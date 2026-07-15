#!/usr/bin/env python3
"""Calibrate frozen FastQC quality severities on the calibration split."""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import fastqc_runtime as fastqc_parser
    import paper_runtime as runtime
    import run_complementarity as execution_runner
except ModuleNotFoundError:
    from scripts import fastqc_runtime as fastqc_parser
    from scripts import paper_runtime as runtime
    from scripts import run_complementarity as execution_runner


SCHEMA_VERSION = "0.1.0"
ANALYZER_VERSION = "0.1.0"
CALL_FIELDS = (
    "analysis_id",
    "execution_id",
    "condition_id",
    "configuration_accession",
    "operator_id",
    "variant",
    "severity_index",
    "case_id",
    "module_names",
    "clean_statuses",
    "observed_statuses",
    "status_worsened",
    "clean_score",
    "observed_score",
    "score_decrease",
    "score_threshold",
    "detected",
)
SEVERITY_FIELDS = (
    "analysis_id",
    "operator_id",
    "variant",
    "severity_index",
    "severity_json",
    "independent_configurations",
    "conditions",
    "detected",
    "detection_rate",
    "eligible",
)
INVARIANCE_FIELDS = (
    "analysis_id",
    "execution_id",
    "condition_id",
    "configuration_accession",
    "operator_id",
    "severity_index",
    "metrics_equal",
    "assessments_equal",
    "passed",
)
CONTEXT_FIELDS = {
    "execution_id",
    "materialization_id",
    "selection_id",
    "condition_id",
    "configuration_accession",
    "family_id",
    "modality",
    "condition_kind",
    "operator_id",
    "variant",
    "event_fraction",
    "mutation_seed",
    "mutated_record_count",
    "total_record_count",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calibrate FastQC quality-control severities."
    )
    parser.add_argument("--execution-manifest", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = analyze_complementarity_calibration(
            execution_path=args.execution_manifest.resolve(),
            protocol_path=args.protocol.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"analyze_complementarity_calibration: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def analyze_complementarity_calibration(
    *, execution_path: Path, protocol_path: Path, output_root: Path
) -> dict[str, Any]:
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    loaded = load_execution(execution_path)
    execution = loaded["execution"]
    if execution.get("cohort_split") != "calibration":
        raise ValueError("quality calibration requires the calibration split")
    protocol = runtime.load_json(protocol_path)
    validate_protocol(protocol)
    if execution["protocol_sha256"] != runtime.file_sha256(protocol_path):
        raise ValueError("analysis protocol differs from execution")
    tools = {
        "analyzer": runtime.script_identity(
            Path(__file__).resolve(), version=ANALYZER_VERSION
        ),
        "fastqc_parser": runtime.script_identity(
            Path(fastqc_parser.__file__).resolve(),
            version=fastqc_parser.RUNTIME_VERSION,
        ),
        "runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
    }
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "execution_id": execution["execution_id"],
        "execution_sha256": runtime.file_sha256(execution_path),
        "protocol_sha256": runtime.file_sha256(protocol_path),
        "analyzer": runtime.functional_script_identity(tools["analyzer"]),
        "fastqc_parser": runtime.functional_script_identity(tools["fastqc_parser"]),
        "runtime": runtime.functional_script_identity(tools["runtime"]),
    }
    analysis_id = runtime.sha256_json(stable_identity)[:16]
    module_map = build_fastqc_module_map(
        loaded["fastqc_runs"], loaded["fastqc_modules"]
    )
    quality_calls = build_quality_calls(
        analysis_id=analysis_id,
        execution_id=execution["execution_id"],
        base_conditions=loaded["base_conditions"],
        controls=loaded["controls"],
        module_map=module_map,
        protocol=protocol,
    )
    severity_rows = summarize_severities(
        analysis_id=analysis_id,
        calls=quality_calls,
        controls=loaded["controls"],
        protocol=protocol,
    )
    invariance = build_seqcheck_invariance(
        analysis_id=analysis_id,
        execution_id=execution["execution_id"],
        base_conditions=loaded["base_conditions"],
        controls=loaded["controls"],
        metrics=loaded["seqcheck_metrics"],
        assessments=loaded["seqcheck_assessments"],
    )
    policy = build_quality_policy(
        analysis_id=analysis_id,
        execution_id=execution["execution_id"],
        protocol_path=protocol_path,
        protocol=protocol,
        severity_rows=severity_rows,
    )
    validation = validate_outputs(
        analysis_id=analysis_id,
        controls=loaded["controls"],
        calls=quality_calls,
        severity_rows=severity_rows,
        invariance=invariance,
        policy=policy,
        protocol=protocol,
    )
    if not validation["valid"]:
        raise ValueError(
            "complementarity calibration analysis failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )

    calls_path = output_root / "tables" / "quality_calls.csv"
    severity_path = output_root / "tables" / "severity_summary.csv"
    invariance_path = output_root / "tables" / "seqcheck_invariance.csv"
    policy_path = output_root / "policy" / "quality_policy.json"
    validation_path = output_root / "validation" / "calibration_analysis.json"
    manifest_path = output_root / "manifests" / "analysis.json"
    runtime.write_csv(calls_path, quality_calls, list(CALL_FIELDS))
    runtime.write_csv(severity_path, severity_rows, list(SEVERITY_FIELDS))
    runtime.write_csv(invariance_path, invariance, list(INVARIANCE_FIELDS))
    runtime.write_json(policy_path, policy)
    runtime.write_json(validation_path, validation)
    manifest = {
        **stable_identity,
        "analysis_id": analysis_id,
        "created_at": utc_now(),
        "valid": True,
        "inputs": {
            "execution": runtime.file_identity(execution_path),
            "protocol": runtime.file_identity(protocol_path),
        },
        "tools": tools,
        "counts": validation["counts"],
        "policy_frozen": policy["frozen"],
        "outputs": {
            "quality_calls": runtime.file_identity(calls_path),
            "severity_summary": runtime.file_identity(severity_path),
            "seqcheck_invariance": runtime.file_identity(invariance_path),
            "quality_policy": runtime.file_identity(policy_path),
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_execution(path: Path) -> dict[str, Any]:
    execution = runtime.load_json(path)
    if execution.get("valid") is not True:
        raise ValueError("complementarity execution is not valid")
    stable = {
        field: execution.get(field)
        for field in execution_runner.EXECUTION_IDENTITY_FIELDS
    }
    if any(value is None for value in stable.values()):
        raise ValueError("complementarity execution identity is incomplete")
    if runtime.sha256_json(stable)[:16] != execution.get("execution_id"):
        raise ValueError("complementarity execution id is not content-addressed")
    materialization_path = verified_path(
        execution.get("inputs", {}).get("materialization", {}), "materialization"
    )
    complementarity, _, base_conditions, controls = (
        execution_runner.load_materialization(materialization_path)
    )
    if complementarity["complementarity_id"] != execution["complementarity_id"]:
        raise ValueError("execution materialization id differs")
    tables = {
        label: runtime.read_csv(
            verified_path(execution.get("outputs", {}).get(label, {}), label)
        )
        for label in (
            "seqspec_runs",
            "seqcheck_runs",
            "seqcheck_metrics",
            "seqcheck_assessments",
            "fastqc_runs",
            "fastqc_modules",
        )
    }
    for label, rows in tables.items():
        if any(
            value.get("execution_id") != execution["execution_id"] for value in rows
        ):
            raise ValueError(f"{label} execution identifiers differ")
    condition_ids = {value["condition_id"] for value in [*base_conditions, *controls]}
    for label, rows in tables.items():
        if any(value.get("condition_id") not in condition_ids for value in rows):
            raise ValueError(f"{label} condition identifiers differ")
    validate_execution_tables(
        execution=execution,
        conditions=[*base_conditions, *controls],
        tables=tables,
    )
    for row in tables["fastqc_runs"]:
        for field in ("data", "summary"):
            if (
                runtime.file_sha256(Path(row[f"{field}_path"]))
                != row[f"{field}_sha256"]
            ):
                raise ValueError(f"FastQC {field} hash changed")
    return {
        "execution": execution,
        "complementarity": complementarity,
        "base_conditions": base_conditions,
        "controls": controls,
        **tables,
    }


def validate_execution_tables(
    *,
    execution: dict[str, Any],
    conditions: list[dict[str, Any]],
    tables: dict[str, list[dict[str, str]]],
) -> None:
    counts = execution.get("counts", {})
    count_fields = {
        "seqspec_runs": "seqspec_runs",
        "seqcheck_runs": "seqcheck_runs",
        "seqcheck_metrics": "seqcheck_metrics",
        "seqcheck_assessments": "seqcheck_assessments",
        "fastqc_runs": "fastqc_input_runs",
        "fastqc_modules": "fastqc_modules",
    }
    for table, count_field in count_fields.items():
        if counts.get(count_field) != len(tables[table]):
            raise ValueError(f"{table} count differs from execution manifest")

    condition_ids = Counter(value["condition_id"] for value in conditions)
    for table in ("seqspec_runs", "seqcheck_runs"):
        if Counter(value["condition_id"] for value in tables[table]) != condition_ids:
            raise ValueError(f"{table} does not reconcile one-to-one with conditions")

    metric_counts = Counter(
        value["condition_id"] for value in tables["seqcheck_metrics"]
    )
    assessment_counts = Counter(
        value["condition_id"] for value in tables["seqcheck_assessments"]
    )
    for row in tables["seqcheck_runs"]:
        condition_id = row["condition_id"]
        if (
            nonnegative_int(row.get("metric_count"), "seqcheck metric count")
            != metric_counts[condition_id]
        ):
            raise ValueError(f"{condition_id}: seqcheck metric count differs")
        if (
            nonnegative_int(row.get("assessment_count"), "seqcheck assessment count")
            != assessment_counts[condition_id]
        ):
            raise ValueError(f"{condition_id}: seqcheck assessment count differs")

    expected_fastqc = Counter(
        condition["condition_id"]
        for condition in conditions
        for _ in condition["inputs"]
    )
    if (
        Counter(value["condition_id"] for value in tables["fastqc_runs"])
        != expected_fastqc
    ):
        raise ValueError("fastqc_runs does not reconcile with condition inputs")
    module_counts = Counter(value["condition_id"] for value in tables["fastqc_modules"])
    expected_modules = Counter()
    for row in tables["fastqc_runs"]:
        expected_modules[row["condition_id"]] += nonnegative_int(
            row.get("module_count"), "FastQC module count"
        )
    if module_counts != expected_modules:
        raise ValueError("fastqc_modules does not reconcile with FastQC runs")


def build_fastqc_module_map(
    fastqc_runs: list[dict[str, str]], fastqc_modules: list[dict[str, str]]
) -> dict[tuple[str, str], dict[str, dict[str, Any]]]:
    rows_by_key = defaultdict(list)
    for row in fastqc_modules:
        rows_by_key[(row["condition_id"], row["case_id"])].append(row)
    run_keys = [(value["condition_id"], value["case_id"]) for value in fastqc_runs]
    if len(run_keys) != len(set(run_keys)):
        raise ValueError("FastQC run keys are not unique")
    if set(rows_by_key) != set(run_keys):
        raise ValueError("FastQC module keys do not reconcile with runs")
    result = {}
    parsed_cache = {}
    for run in fastqc_runs:
        key = (run["condition_id"], run["case_id"])
        data_path = Path(run["data_path"])
        parsed = parsed_cache.get(data_path)
        if parsed is None:
            parsed = fastqc_parser.parse_fastqc_output(data_path.parent)
            parsed_cache[data_path] = parsed
        modules = {value["name"]: value for value in parsed["modules"]}
        rows = rows_by_key[key]
        observed = {(value["module_name"], value["status"]) for value in rows}
        expected = {(name, value["status"]) for name, value in modules.items()}
        if observed != expected or len(rows) != len(modules):
            raise ValueError(f"FastQC module rows differ for {key}")
        result[key] = modules
    return result


def build_quality_calls(
    *,
    analysis_id: str,
    execution_id: str,
    base_conditions: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    module_map: dict[tuple[str, str], dict[str, dict[str, Any]]],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    clean = {
        value["configuration_accession"]: value
        for value in base_conditions
        if value["condition_kind"] == "clean"
    }
    expected_clean = sum(
        value["condition_kind"] == "clean" for value in base_conditions
    )
    if len(clean) != expected_clean:
        raise ValueError("clean configuration accessions are not unique")
    threshold = float(protocol["calibration"]["minimum_quality_score_decrease"])
    result = []
    for condition in controls:
        if condition["condition_kind"] != "quality":
            continue
        configuration = condition["configuration_accession"]
        clean_condition = clean[configuration]
        case_id = only(condition["target"]["case_ids"], "quality target case")
        clean_modules = module_map[(clean_condition["condition_id"], case_id)]
        observed_modules = module_map[(condition["condition_id"], case_id)]
        module_names = condition["expected"]["fastqc_modules"]
        missing = [
            name
            for name in module_names
            if name not in clean_modules or name not in observed_modules
        ]
        if missing:
            raise ValueError(
                f"{condition['condition_id']}: missing FastQC modules {missing}"
            )
        status_worsened = any(
            fastqc_parser.status_worsened(
                clean_modules[name]["status"], observed_modules[name]["status"]
            )
            for name in module_names
        )
        score_pairs = [
            (module_score(clean_modules[name]), module_score(observed_modules[name]))
            for name in module_names
        ]
        decreases = [
            clean_score - observed_score
            for clean_score, observed_score in score_pairs
            if clean_score is not None and observed_score is not None
        ]
        score_decrease = max(decreases, default=None)
        best_pair = (None, None)
        if score_decrease is not None:
            best_pair = next(
                (
                    pair
                    for pair in score_pairs
                    if pair[0] is not None
                    and pair[1] is not None
                    and math.isclose(pair[0] - pair[1], score_decrease)
                ),
                best_pair,
            )
        detected = status_worsened or (
            score_decrease is not None and score_decrease >= threshold
        )
        result.append(
            {
                "analysis_id": analysis_id,
                "execution_id": execution_id,
                "condition_id": condition["condition_id"],
                "configuration_accession": configuration,
                "operator_id": condition["operator_id"],
                "variant": condition["variant"],
                "severity_index": condition["severity_index"],
                "case_id": case_id,
                "module_names": ";".join(module_names),
                "clean_statuses": ";".join(
                    clean_modules[name]["status"] for name in module_names
                ),
                "observed_statuses": ";".join(
                    observed_modules[name]["status"] for name in module_names
                ),
                "status_worsened": status_worsened,
                "clean_score": best_pair[0],
                "observed_score": best_pair[1],
                "score_decrease": score_decrease,
                "score_threshold": threshold,
                "detected": detected,
            }
        )
    return sorted(
        result,
        key=lambda value: (
            value["operator_id"],
            value["severity_index"],
            value["configuration_accession"],
        ),
    )


def module_score(module: dict[str, Any]) -> float | None:
    name = module["name"]
    if name == "Per base sequence quality":
        header = table_header(module, "Base")
        if "Mean" not in header:
            raise ValueError("FastQC per-base quality has no Mean column")
        index = header.index("Mean")
        values = [
            finite_float(row[index]) for row in module["rows"] if len(row) > index
        ]
        values = [value for value in values if value is not None]
        return min(values) if values else None
    if name == "Per sequence quality scores":
        header = table_header(module, "Quality")
        if "Count" not in header:
            raise ValueError("FastQC sequence quality has no Count column")
        quality_index = header.index("Quality")
        count_index = header.index("Count")
        values = []
        for row in module["rows"]:
            if len(row) <= max(quality_index, count_index):
                continue
            quality = finite_float(row[quality_index])
            count = finite_float(row[count_index])
            if quality is not None and count is not None and count >= 0:
                values.append((quality, count))
        total = sum(value[1] for value in values)
        return (
            sum(quality * count for quality, count in values) / total if total else None
        )
    return None


def table_header(module: dict[str, Any], first_field: str) -> list[str]:
    matches = [
        value for value in module["headers"] if value and value[0] == first_field
    ]
    if len(matches) != 1:
        raise ValueError(f"FastQC {module['name']} has no unique table header")
    return matches[0]


def summarize_severities(
    *,
    analysis_id: str,
    calls: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    severities = {
        (value["operator_id"], value["severity_index"]): value["severity"]
        for value in controls
        if value["condition_kind"] == "quality"
    }
    grouped = defaultdict(list)
    for call in calls:
        grouped[(call["operator_id"], call["severity_index"])].append(call)
    minimum = protocol["calibration"]["minimum_independent_configurations"]
    rate_min = protocol["calibration"]["minimum_fastqc_detection_rate"]
    variants = {
        value["operator_id"]: value["variant"] for value in protocol["quality_controls"]
    }
    result = []
    for key, values in sorted(grouped.items()):
        configurations = {value["configuration_accession"] for value in values}
        if len(configurations) != len(values):
            raise ValueError(f"{key[0]} severity {key[1]} has duplicate configurations")
        detected = sum(value["detected"] for value in values)
        rate = detected / len(configurations)
        result.append(
            {
                "analysis_id": analysis_id,
                "operator_id": key[0],
                "variant": variants[key[0]],
                "severity_index": key[1],
                "severity_json": runtime.canonical_json(severities[key]),
                "independent_configurations": len(configurations),
                "conditions": len(values),
                "detected": detected,
                "detection_rate": rate,
                "eligible": len(configurations) >= minimum and rate >= rate_min,
            }
        )
    return result


def build_quality_policy(
    *,
    analysis_id: str,
    execution_id: str,
    protocol_path: Path,
    protocol: dict[str, Any],
    severity_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    by_operator = defaultdict(list)
    for row in severity_rows:
        by_operator[row["operator_id"]].append(row)
    entries = []
    for operator in protocol["quality_controls"]:
        eligible = sorted(
            (
                value
                for value in by_operator[operator["operator_id"]]
                if value["eligible"]
            ),
            key=lambda value: value["severity_index"],
        )
        selected = eligible[0] if eligible else None
        entries.append(
            {
                "operator_id": operator["operator_id"],
                "variant": operator["variant"],
                "severity_index": selected["severity_index"] if selected else None,
                "severity": (
                    operator["severity_candidates"][selected["severity_index"]]
                    if selected
                    else None
                ),
                "independent_configurations": (
                    selected["independent_configurations"] if selected else 0
                ),
                "detection_rate": selected["detection_rate"] if selected else None,
                "fastqc_modules": operator["fastqc_modules"],
                "detection_rule": {
                    "status_worsening": True,
                    "minimum_quality_score_decrease": protocol["calibration"][
                        "minimum_quality_score_decrease"
                    ],
                },
                "ready": selected is not None,
            }
        )
    frozen = all(value["ready"] for value in entries)
    stable = {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": analysis_id,
        "execution_id": execution_id,
        "protocol_sha256": runtime.file_sha256(protocol_path),
        "frozen": frozen,
        "entries": entries,
    }
    return {
        **stable,
        "policy_id": runtime.sha256_json(stable)[:16],
        "created_at": utc_now(),
    }


def build_seqcheck_invariance(
    *,
    analysis_id: str,
    execution_id: str,
    base_conditions: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    metrics: list[dict[str, str]],
    assessments: list[dict[str, str]],
) -> list[dict[str, Any]]:
    clean = {
        value["configuration_accession"]: value["condition_id"]
        for value in base_conditions
        if value["condition_kind"] == "clean"
    }
    metric_map = canonical_rows_by_condition(metrics)
    assessment_map = canonical_rows_by_condition(assessments)
    result = []
    for condition in controls:
        if condition["condition_kind"] != "quality":
            continue
        clean_id = clean[condition["configuration_accession"]]
        metrics_equal = metric_map.get(condition["condition_id"], []) == metric_map.get(
            clean_id, []
        )
        assessments_equal = assessment_map.get(
            condition["condition_id"], []
        ) == assessment_map.get(clean_id, [])
        result.append(
            {
                "analysis_id": analysis_id,
                "execution_id": execution_id,
                "condition_id": condition["condition_id"],
                "configuration_accession": condition["configuration_accession"],
                "operator_id": condition["operator_id"],
                "severity_index": condition["severity_index"],
                "metrics_equal": metrics_equal,
                "assessments_equal": assessments_equal,
                "passed": metrics_equal and assessments_equal,
            }
        )
    return result


def canonical_rows_by_condition(
    rows: list[dict[str, str]],
) -> dict[str, list[str]]:
    result = defaultdict(list)
    for row in rows:
        stable = {key: value for key, value in row.items() if key not in CONTEXT_FIELDS}
        result[row["condition_id"]].append(runtime.canonical_json(stable))
    return {key: sorted(values) for key, values in result.items()}


def validate_protocol(protocol: dict[str, Any]) -> None:
    execution_runner.validate_protocol(protocol)
    calibration = protocol.get("calibration", {})
    minimum = calibration.get("minimum_independent_configurations")
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum <= 0:
        raise ValueError("minimum quality calibration configurations is invalid")
    rate = calibration.get("minimum_fastqc_detection_rate")
    if (
        isinstance(rate, bool)
        or not isinstance(rate, (int, float))
        or not 0 < rate <= 1
    ):
        raise ValueError("minimum FastQC detection rate is invalid")
    decrease = calibration.get("minimum_quality_score_decrease")
    if (
        isinstance(decrease, bool)
        or not isinstance(decrease, (int, float))
        or decrease <= 0
    ):
        raise ValueError("minimum quality score decrease is invalid")
    if calibration.get("severity_selection") != "least_severe_candidate_meeting_target":
        raise ValueError("quality severity selection rule is invalid")


def validate_outputs(
    *,
    analysis_id: str,
    controls: list[dict[str, Any]],
    calls: list[dict[str, Any]],
    severity_rows: list[dict[str, Any]],
    invariance: list[dict[str, Any]],
    policy: dict[str, Any],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    errors = []
    quality = [value for value in controls if value["condition_kind"] == "quality"]
    expected_ids = Counter(value["condition_id"] for value in quality)
    if Counter(value["condition_id"] for value in calls) != expected_ids:
        errors.append("quality calls do not reconcile with controls")
    if Counter(value["condition_id"] for value in invariance) != expected_ids:
        errors.append("seqcheck invariance rows do not reconcile with controls")
    if any(value["analysis_id"] != analysis_id for value in calls + invariance):
        errors.append("analysis identifiers differ")
    if any(value["analysis_id"] != analysis_id for value in severity_rows):
        errors.append("severity analysis identifiers differ")
    expected_severities = Counter(
        (operator["operator_id"], index)
        for operator in protocol["quality_controls"]
        for index in range(len(operator["severity_candidates"]))
    )
    observed_severities = Counter(
        (value["operator_id"], value["severity_index"]) for value in severity_rows
    )
    if observed_severities != expected_severities:
        errors.append("severity summaries do not reconcile with protocol")
    expected_operators = {
        value["operator_id"] for value in protocol["quality_controls"]
    }
    if {value["operator_id"] for value in policy["entries"]} != expected_operators:
        errors.append("quality policy entries do not reconcile with protocol")
    stable_policy = {
        key: value
        for key, value in policy.items()
        if key not in {"policy_id", "created_at"}
    }
    if runtime.sha256_json(stable_policy)[:16] != policy["policy_id"]:
        errors.append("quality policy id is not content-addressed")
    return {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": analysis_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "policy_frozen": policy["frozen"],
        "quality_seqcheck_invariance_passed": all(
            value["passed"] for value in invariance
        ),
        "counts": {
            "quality_conditions": len(quality),
            "quality_calls": len(calls),
            "severity_summaries": len(severity_rows),
            "eligible_severities": sum(value["eligible"] for value in severity_rows),
            "policy_entries": len(policy["entries"]),
            "ready_policy_entries": sum(value["ready"] for value in policy["entries"]),
            "seqcheck_invariance_rows": len(invariance),
            "seqcheck_invariance_failures": sum(
                not value["passed"] for value in invariance
            ),
        },
    }


def finite_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def only(values: list[Any], label: str) -> Any:
    if len(values) != 1:
        raise ValueError(f"{label} must contain exactly one value")
    return values[0]


def nonnegative_int(value: Any, label: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} is invalid") from error
    if result < 0 or str(result) != str(value):
        raise ValueError(f"{label} is invalid")
    return result


def verified_path(identity: Any, label: str) -> Path:
    if not isinstance(identity, dict):
        raise ValueError(f"{label} identity is missing")
    path = Path(str(identity.get("path", ""))).resolve()
    if runtime.file_sha256(path) != identity.get("sha256"):
        raise ValueError(f"{label} hash changed")
    if path.stat().st_size != identity.get("size_bytes"):
        raise ValueError(f"{label} size changed")
    return path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
