#!/usr/bin/env python3
"""Run and export every materialized perturbation condition."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import materialize_perturbations as materializer
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import materialize_perturbations as materializer
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
RUNNER_VERSION = "0.1.0"
MATERIALIZATION_IDENTITY_FIELDS = (
    "schema_version",
    "selection_id",
    "inventory_id",
    "study_run_id",
    "sampling_policy_id",
    "sample_source",
    "base_samples",
    "perturbation_protocol_sha256",
    "seqspec",
    "yq",
    "materializer",
    "mutation_runtime",
    "paper_runtime",
    "fastq_runtime",
)
RUN_FIELDS = (
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
    "expected_process_outcome",
    "observed_process_outcome",
    "process_outcome_matches",
    "exit_code",
    "timed_out",
    "report_schema_version",
    "result_count",
    "metric_count",
    "assessment_count",
    "report_path",
    "report_sha256",
    "stdout_path",
    "stdout_sha256",
    "stderr_path",
    "stderr_sha256",
    "command_json",
)
METRIC_FIELDS = (
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
    "result_index",
    "check",
    "files",
    "reads",
    "regions",
    "ontology_terms",
    "metric_side",
    "metric_id",
    "metric_name",
    "metric_description",
    "data_kind",
    "unit",
    "value_json",
)
ASSESSMENT_FIELDS = (
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
    "result_index",
    "assessment_index",
    "check",
    "files",
    "reads",
    "regions",
    "ontology_terms",
    "assessment_type",
    "assessment_code",
    "assessment_description",
    "expected_metric_ids",
    "observed_metric_ids",
)
PERFORMANCE_FIELDS = (
    "execution_id",
    "materialization_id",
    "condition_id",
    "configuration_accession",
    "operator_id",
    "variant",
    "runtime_schema_version",
    "command_json",
    "exit_code",
    "timed_out",
    "wall_time_seconds",
    "user_cpu_seconds",
    "system_cpu_seconds",
    "peak_resident_memory_bytes",
    "records_processed",
    "compressed_bytes_read",
    "report_size_bytes",
    "stdout_path",
    "stdout_sha256",
    "stderr_path",
    "stderr_sha256",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run every validated controlled perturbation condition."
    )
    parser.add_argument("--materialization-manifest", required=True, type=Path)
    parser.add_argument("--execution-protocol", required=True, type=Path)
    parser.add_argument("--seqcheck-bin", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=300)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = run_perturbation_calibration(
            materialization_manifest_path=args.materialization_manifest.resolve(),
            execution_protocol_path=args.execution_protocol.resolve(),
            seqcheck_bin=args.seqcheck_bin.resolve(),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
        )
    except (OSError, ValueError) as error:
        print(f"run_perturbation_calibration: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def run_perturbation_calibration(
    *,
    materialization_manifest_path: Path,
    execution_protocol_path: Path,
    seqcheck_bin: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    materialization_manifest_path = materialization_manifest_path.resolve()
    execution_protocol_path = execution_protocol_path.resolve()
    seqcheck_bin = seqcheck_bin.resolve()
    output_root = output_root.resolve()
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    materialization, conditions = load_materialization(materialization_manifest_path)
    protocol = runtime.load_json(execution_protocol_path)
    validate_protocol(protocol)
    seqcheck = runtime.executable_identity(
        seqcheck_bin, timeout_seconds=min(timeout_seconds, 30)
    )
    tools = {
        "seqcheck": seqcheck,
        "runner": runtime.script_identity(
            Path(__file__).resolve(), version=RUNNER_VERSION
        ),
        "runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
    }
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "materialization_id": materialization["materialization_id"],
        "materialization_sha256": runtime.file_sha256(materialization_manifest_path),
        "execution_protocol_sha256": runtime.file_sha256(execution_protocol_path),
        "seqcheck": runtime.functional_executable_identity(seqcheck),
        "runner": runtime.functional_script_identity(tools["runner"]),
        "runtime": runtime.functional_script_identity(tools["runtime"]),
    }
    execution_id = runtime.sha256_json(stable_identity)[:16]
    output_root.mkdir(parents=True)
    run_rows = []
    metric_rows = []
    assessment_rows = []
    performance_rows = []
    for condition in conditions:
        result = run_condition(
            condition=condition,
            protocol=protocol,
            execution_id=execution_id,
            materialization_id=materialization["materialization_id"],
            seqcheck_bin=seqcheck_bin,
            output_root=output_root,
            timeout_seconds=timeout_seconds,
        )
        run_rows.append(result["run"])
        metric_rows.extend(result["metrics"])
        assessment_rows.extend(result["assessments"])
        performance_rows.append(result["performance"])

    tables = output_root / "tables"
    runs_path = tables / "runs.csv"
    metrics_path = tables / "metrics.csv"
    assessments_path = tables / "assessments.csv"
    performance_path = tables / "performance.csv"
    runtime.write_csv(runs_path, run_rows, list(RUN_FIELDS))
    runtime.write_csv(metrics_path, metric_rows, list(METRIC_FIELDS))
    runtime.write_csv(assessments_path, assessment_rows, list(ASSESSMENT_FIELDS))
    runtime.write_csv(performance_path, performance_rows, list(PERFORMANCE_FIELDS))
    validation = validate_outputs(
        execution_id=execution_id,
        conditions=conditions,
        runs=run_rows,
        metrics=metric_rows,
        assessments=assessment_rows,
        performance=performance_rows,
        output_root=output_root,
    )
    validation_path = output_root / "validation" / "perturbation_execution.json"
    runtime.write_json(validation_path, validation)
    if not validation["valid"]:
        raise ValueError(
            "perturbation execution failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )
    manifest_path = output_root / "manifests" / "execution.json"
    manifest = {
        **stable_identity,
        "execution_id": execution_id,
        "selection_id": materialization["selection_id"],
        "created_at": utc_now(),
        "valid": True,
        "inputs": {
            "materialization": runtime.file_identity(materialization_manifest_path),
            "execution_protocol": runtime.file_identity(execution_protocol_path),
        },
        "tools": tools,
        "counts": validation["counts"],
        "outputs": {
            "runs": runtime.file_identity(runs_path),
            "metrics": runtime.file_identity(metrics_path),
            "assessments": runtime.file_identity(assessments_path),
            "performance": runtime.file_identity(performance_path),
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_materialization(
    path: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    value = runtime.load_json(path)
    if value.get("valid") is not True:
        raise ValueError("perturbation materialization is not valid")
    stable = {field: value.get(field) for field in MATERIALIZATION_IDENTITY_FIELDS}
    if any(stable[field] is None for field in MATERIALIZATION_IDENTITY_FIELDS):
        raise ValueError("materialization identity is incomplete")
    if runtime.sha256_json(stable)[:16] != value.get("materialization_id"):
        raise ValueError("materialization id is not content-addressed")
    conditions_path = verified_path(
        value.get("outputs", {}).get("conditions", {}), "condition manifest"
    )
    payload = runtime.load_json(conditions_path)
    if payload.get("materialization_id") != value["materialization_id"]:
        raise ValueError("condition manifest materialization id differs")
    conditions = payload.get("conditions")
    if not isinstance(conditions, list) or not conditions:
        raise ValueError("materialization has no conditions")
    for condition in conditions:
        validate_condition(condition, value["materialization_id"])
    if value.get("counts", {}).get("conditions") != len(conditions):
        raise ValueError("materialized condition count does not reconcile")
    return value, conditions


def validate_condition(condition: dict[str, Any], materialization_id: str) -> None:
    if condition.get("materialization_id") != materialization_id:
        raise ValueError("condition materialization id differs")
    stable = {
        key: value
        for key, value in condition.items()
        if key not in {"condition_id", "spec", "inputs", "mutation"}
    }
    stable["spec"] = materializer.functional_file_identity(condition["spec"])
    stable["inputs"] = [
        materializer.functional_input_identity(value) for value in condition["inputs"]
    ]
    stable["mutation"] = materializer.functional_mutation_details(condition["mutation"])
    if runtime.sha256_json(stable)[:16] != condition.get("condition_id"):
        raise ValueError("condition id is not content-addressed")
    identities = [
        condition["spec"],
        *condition["inputs"],
        *materializer.nested_file_identities(condition["mutation"]),
    ]
    for identity in identities:
        verified_path(identity, f"condition {condition['condition_id']} artifact")


def validate_protocol(protocol: dict[str, Any]) -> None:
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("execution protocol schema is unsupported")
    command = protocol.get("command", {})
    if command.get("n_reads") != 0 or command.get("format") != "json":
        raise ValueError(
            "execution protocol must inspect all materialized reads as JSON"
        )
    performance = protocol.get("performance", {})
    if performance.get("warmup_runs") != 0:
        raise ValueError("perturbation execution does not use warm-up runs")
    if performance.get("measured_runs") != 1:
        raise ValueError("perturbation execution requires one measured run")
    validation = protocol.get("validation", {})
    for field in (
        "require_expected_process_outcome",
        "require_json_report_on_success",
        "require_no_timeout",
        "require_one_run_per_condition",
    ):
        if validation.get(field) is not True:
            raise ValueError(f"execution validation is not enabled: {field}")


def run_condition(
    *,
    condition: dict[str, Any],
    protocol: dict[str, Any],
    execution_id: str,
    materialization_id: str,
    seqcheck_bin: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    condition_id = condition["condition_id"]
    report_path = output_root / "reports" / f"{condition_id}.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    log_root = output_root / "logs" / condition_id
    argv = [
        str(seqcheck_bin),
        "check",
        "--modality",
        condition["modality"],
        "--spec",
        condition["spec"]["path"],
        "--n-reads",
        str(protocol["command"]["n_reads"]),
        "--format",
        protocol["command"]["format"],
        "--output",
        str(report_path),
        *[value["path"] for value in condition["inputs"]],
    ]
    measurement = runtime.run_measured_command(
        argv,
        stdout_path=log_root / "stdout.txt",
        stderr_path=log_root / "stderr.txt",
        timeout_seconds=timeout_seconds,
    )
    expected_outcome = condition["expected"]["process_outcome"]
    observed_outcome = (
        "success"
        if measurement["exit_code"] == 0 and not measurement["timed_out"]
        else "failure"
    )
    if measurement["timed_out"]:
        raise ValueError(f"{condition_id}: seqcheck timed out")
    if observed_outcome != expected_outcome:
        stderr = Path(measurement["stderr"]["path"]).read_text(
            encoding="utf-8", errors="replace"
        )
        raise ValueError(
            f"{condition_id}: expected process {expected_outcome}, observed "
            f"{observed_outcome}: {stderr.strip()}"
        )
    context = condition_context(
        condition,
        execution_id=execution_id,
        materialization_id=materialization_id,
    )
    payload = None
    report_identity = {"path": "", "sha256": "", "size_bytes": 0}
    metrics = []
    assessments = []
    if observed_outcome == "success":
        if not report_path.is_file():
            raise ValueError(f"{condition_id}: successful run wrote no JSON report")
        payload = runtime.load_json(report_path)
        if not str(payload.get("report_schema_version", "")).strip():
            raise ValueError(f"{condition_id}: report schema version is missing")
        report_identity = runtime.file_identity(report_path)
        metrics = runtime.flatten_report_metrics(payload, context)
        assessments = runtime.flatten_report_assessments(payload, context)
    performance = performance_row(
        condition=condition,
        measurement=measurement,
        execution_id=execution_id,
        materialization_id=materialization_id,
        report_identity=report_identity,
    )
    run = {
        **context,
        "expected_process_outcome": expected_outcome,
        "observed_process_outcome": observed_outcome,
        "process_outcome_matches": observed_outcome == expected_outcome,
        "exit_code": measurement["exit_code"],
        "timed_out": measurement["timed_out"],
        "report_schema_version": (
            payload.get("report_schema_version", "") if payload else ""
        ),
        "result_count": len(payload.get("results", [])) if payload else 0,
        "metric_count": len(metrics),
        "assessment_count": len(assessments),
        "report_path": report_identity["path"],
        "report_sha256": report_identity["sha256"],
        "stdout_path": measurement["stdout"]["path"],
        "stdout_sha256": measurement["stdout"]["sha256"],
        "stderr_path": measurement["stderr"]["path"],
        "stderr_sha256": measurement["stderr"]["sha256"],
        "command_json": runtime.canonical_json(argv),
    }
    return {
        "run": run,
        "metrics": metrics,
        "assessments": assessments,
        "performance": performance,
    }


def condition_context(
    condition: dict[str, Any], *, execution_id: str, materialization_id: str
) -> dict[str, Any]:
    return {
        "execution_id": execution_id,
        "materialization_id": materialization_id,
        "selection_id": condition["selection_id"],
        "condition_id": condition["condition_id"],
        "configuration_accession": condition["configuration_accession"],
        "family_id": condition["family_id"],
        "modality": condition["modality"],
        "condition_kind": condition["condition_kind"],
        "operator_id": condition["operator_id"],
        "variant": condition["variant"],
        "event_fraction": condition["event_fraction"],
        "mutation_seed": condition["mutation_seed"],
        "mutated_record_count": condition["mutated_record_count"],
        "total_record_count": condition["total_record_count"],
    }


def performance_row(
    *,
    condition: dict[str, Any],
    measurement: dict[str, Any],
    execution_id: str,
    materialization_id: str,
    report_identity: dict[str, Any],
) -> dict[str, Any]:
    return {
        "execution_id": execution_id,
        "materialization_id": materialization_id,
        "condition_id": condition["condition_id"],
        "configuration_accession": condition["configuration_accession"],
        "operator_id": condition["operator_id"],
        "variant": condition["variant"],
        "runtime_schema_version": measurement["runtime_schema_version"],
        "command_json": runtime.canonical_json(measurement["argv"]),
        "exit_code": measurement["exit_code"],
        "timed_out": measurement["timed_out"],
        "wall_time_seconds": measurement["wall_time_seconds"],
        "user_cpu_seconds": measurement["user_cpu_seconds"],
        "system_cpu_seconds": measurement["system_cpu_seconds"],
        "peak_resident_memory_bytes": measurement["peak_resident_memory_bytes"],
        "records_processed": condition["total_record_count"],
        "compressed_bytes_read": sum(
            value["size_bytes"] for value in condition["inputs"]
        ),
        "report_size_bytes": report_identity["size_bytes"],
        "stdout_path": measurement["stdout"]["path"],
        "stdout_sha256": measurement["stdout"]["sha256"],
        "stderr_path": measurement["stderr"]["path"],
        "stderr_sha256": measurement["stderr"]["sha256"],
    }


def validate_outputs(
    *,
    execution_id: str,
    conditions: list[dict[str, Any]],
    runs: list[dict[str, Any]],
    metrics: list[dict[str, Any]],
    assessments: list[dict[str, Any]],
    performance: list[dict[str, Any]],
    output_root: Path,
) -> dict[str, Any]:
    errors = []
    condition_ids = [value["condition_id"] for value in conditions]
    run_ids = [value["condition_id"] for value in runs]
    if Counter(run_ids) != Counter(condition_ids):
        errors.append("run rows do not reconcile one-to-one with conditions")
    if len(performance) != len(conditions):
        errors.append("performance rows do not reconcile with conditions")
    if any(value["execution_id"] != execution_id for value in runs):
        errors.append("run execution identifiers differ")
    if any(not value["process_outcome_matches"] for value in runs):
        errors.append("observed process outcomes differ from expected outcomes")
    metric_counts = Counter(value["condition_id"] for value in metrics)
    assessment_counts = Counter(value["condition_id"] for value in assessments)
    for run in runs:
        condition_id = run["condition_id"]
        if run["metric_count"] != metric_counts[condition_id]:
            errors.append(f"{condition_id}: metric count does not reconcile")
        if run["assessment_count"] != assessment_counts[condition_id]:
            errors.append(f"{condition_id}: assessment count does not reconcile")
        for field in ("stdout", "stderr"):
            path = Path(run[f"{field}_path"])
            if runtime.file_sha256(path) != run[f"{field}_sha256"]:
                errors.append(f"{condition_id}: {field} hash changed")
        if run["observed_process_outcome"] == "success":
            report_path = Path(run["report_path"])
            if runtime.file_sha256(report_path) != run["report_sha256"]:
                errors.append(f"{condition_id}: report hash changed")
            if not report_path.is_relative_to(output_root):
                errors.append(f"{condition_id}: report is outside execution root")
    return {
        "schema_version": SCHEMA_VERSION,
        "execution_id": execution_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "counts": {
            "conditions": len(conditions),
            "runs": len(runs),
            "successful_runs": sum(
                value["observed_process_outcome"] == "success" for value in runs
            ),
            "expected_failure_runs": sum(
                value["observed_process_outcome"] == "failure" for value in runs
            ),
            "metrics": len(metrics),
            "assessments": len(assessments),
            "performance_rows": len(performance),
            "runs_by_operator": dict(
                sorted(Counter(value["operator_id"] for value in runs).items())
            ),
        },
    }


def verified_path(identity: Any, label: str) -> Path:
    if not isinstance(identity, dict):
        raise ValueError(f"{label} identity is missing")
    path = Path(str(identity.get("path", ""))).resolve()
    if runtime.file_sha256(path) != identity.get("sha256"):
        raise ValueError(f"{label} hash changed")
    return path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
