#!/usr/bin/env python3
"""Run seqspec, FastQC, and seqcheck on complementarity conditions."""

from __future__ import annotations

import argparse
import copy
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import fastqc_runtime as fastqc_parser
    import materialize_complementarity as materializer
    import materialize_perturbations as base_materializer
    import paper_runtime as runtime
    import run_perturbation_calibration as seqcheck_runner
except ModuleNotFoundError:
    from scripts import fastqc_runtime as fastqc_parser
    from scripts import materialize_complementarity as materializer
    from scripts import materialize_perturbations as base_materializer
    from scripts import paper_runtime as runtime
    from scripts import run_perturbation_calibration as seqcheck_runner


SCHEMA_VERSION = "0.1.0"
RUNNER_VERSION = "0.1.0"
SEQSPEC_RUN_FIELDS = (
    "execution_id",
    "condition_source",
    "condition_id",
    "configuration_accession",
    "condition_kind",
    "operator_id",
    "variant",
    "expected_status",
    "observed_status",
    "runtime_schema_version",
    "exit_code",
    "timed_out",
    "wall_time_seconds",
    "user_cpu_seconds",
    "system_cpu_seconds",
    "peak_resident_memory_bytes",
    "command_json",
    "stdout_path",
    "stdout_sha256",
    "stderr_path",
    "stderr_sha256",
)
FASTQC_RUN_FIELDS = (
    "execution_id",
    "condition_source",
    "condition_id",
    "configuration_accession",
    "condition_kind",
    "operator_id",
    "variant",
    "case_id",
    "fastq_accession",
    "read_id",
    "input_path",
    "input_sha256",
    "cache_key",
    "cache_hit",
    "cache_source_condition_id",
    "fastqc_version",
    "module_count",
    "command_json",
    "data_path",
    "data_sha256",
    "summary_path",
    "summary_sha256",
)
FASTQC_MODULE_FIELDS = (
    "execution_id",
    "condition_source",
    "condition_id",
    "configuration_accession",
    "condition_kind",
    "operator_id",
    "variant",
    "case_id",
    "fastq_accession",
    "read_id",
    "module_name",
    "status",
    "status_rank",
    "data_path",
    "data_sha256",
    "summary_path",
    "summary_sha256",
)
FASTQC_PERFORMANCE_FIELDS = (
    "execution_id",
    "cache_key",
    "source_condition_id",
    "case_id",
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
    "stdout_path",
    "stdout_sha256",
    "stderr_path",
    "stderr_sha256",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the three-tool QC complementarity experiment."
    )
    parser.add_argument("--materialization-manifest", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--seqcheck-bin", required=True, type=Path)
    parser.add_argument("--seqspec-bin", required=True, type=Path)
    parser.add_argument("--fastqc-bin", required=True, type=Path)
    parser.add_argument("--fastqc-limits", type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=300)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = run_complementarity(
            materialization_path=args.materialization_manifest.resolve(),
            protocol_path=args.protocol.resolve(),
            seqcheck_bin=args.seqcheck_bin.resolve(),
            seqspec_bin=args.seqspec_bin.resolve(),
            fastqc_bin=args.fastqc_bin.resolve(),
            fastqc_limits=(
                args.fastqc_limits.resolve() if args.fastqc_limits else None
            ),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
        )
    except (OSError, ValueError) as error:
        print(f"run_complementarity: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def run_complementarity(
    *,
    materialization_path: Path,
    protocol_path: Path,
    seqcheck_bin: Path,
    seqspec_bin: Path,
    fastqc_bin: Path,
    output_root: Path,
    timeout_seconds: int,
    fastqc_limits: Path | None = None,
) -> dict[str, Any]:
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    complementarity, base, base_conditions, controls = load_materialization(
        materialization_path
    )
    protocol = runtime.load_json(protocol_path)
    validate_protocol(protocol)
    if runtime.file_sha256(protocol_path) != complementarity["protocol_sha256"]:
        raise ValueError("execution protocol differs from materialization")
    tools = tool_identities(
        seqcheck_bin=seqcheck_bin,
        seqspec_bin=seqspec_bin,
        fastqc_bin=fastqc_bin,
        fastqc_limits=fastqc_limits,
        timeout_seconds=timeout_seconds,
    )
    if (
        runtime.functional_executable_identity(tools["seqspec"])
        != complementarity["seqspec"]
    ):
        raise ValueError("seqspec executable differs from materialization")
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "complementarity_id": complementarity["complementarity_id"],
        "materialization_sha256": runtime.file_sha256(materialization_path),
        "protocol_sha256": runtime.file_sha256(protocol_path),
        "seqcheck": runtime.functional_executable_identity(tools["seqcheck"]),
        "seqspec": runtime.functional_executable_identity(tools["seqspec"]),
        "fastqc": functional_fastqc_identity(tools),
        "runner": runtime.functional_script_identity(tools["runner"]),
        "fastqc_parser": runtime.functional_script_identity(tools["fastqc_parser"]),
        "seqcheck_runner": runtime.functional_script_identity(tools["seqcheck_runner"]),
        "runtime": runtime.functional_script_identity(tools["runtime"]),
    }
    execution_id = runtime.sha256_json(stable_identity)[:16]
    combined = [("structural", value) for value in base_conditions] + [
        ("control", value) for value in controls
    ]
    if len({value["condition_id"] for _, value in combined}) != len(combined):
        raise ValueError("combined condition identifiers are not unique")

    output_root.mkdir(parents=True)
    seqspec_rows = []
    seqcheck_runs = []
    seqcheck_metrics = []
    seqcheck_assessments = []
    seqcheck_performance = []
    fastqc_runs = []
    fastqc_modules = []
    fastqc_performance = []
    fastqc_cache = {}
    for source, condition in combined:
        seqspec_rows.append(
            run_seqspec_condition(
                execution_id=execution_id,
                source=source,
                condition=condition,
                seqspec_bin=seqspec_bin,
                output_root=output_root,
                timeout_seconds=timeout_seconds,
            )
        )
        adapted = adapt_seqcheck_condition(condition, source=source)
        seqcheck_result = seqcheck_runner.run_condition(
            condition=adapted,
            protocol={"command": protocol["seqcheck"]},
            execution_id=execution_id,
            materialization_id=complementarity["complementarity_id"],
            seqcheck_bin=seqcheck_bin,
            output_root=output_root / "seqcheck",
            timeout_seconds=timeout_seconds,
        )
        seqcheck_result["run"]["condition_source"] = source
        seqcheck_runs.append(seqcheck_result["run"])
        seqcheck_metrics.extend(seqcheck_result["metrics"])
        seqcheck_assessments.extend(seqcheck_result["assessments"])
        seqcheck_performance.append(seqcheck_result["performance"])
        for input_value in condition["inputs"]:
            result = run_fastqc_input(
                execution_id=execution_id,
                source=source,
                condition=condition,
                input_value=input_value,
                protocol=protocol,
                tools=tools,
                fastqc_bin=fastqc_bin,
                fastqc_limits=fastqc_limits,
                output_root=output_root,
                timeout_seconds=timeout_seconds,
                cache=fastqc_cache,
            )
            fastqc_runs.append(result["run"])
            fastqc_modules.extend(result["modules"])
            if result["performance"] is not None:
                fastqc_performance.append(result["performance"])

    tables = output_root / "tables"
    seqspec_path = tables / "seqspec_runs.csv"
    seqcheck_runs_path = tables / "seqcheck_runs.csv"
    seqcheck_metrics_path = tables / "seqcheck_metrics.csv"
    seqcheck_assessments_path = tables / "seqcheck_assessments.csv"
    seqcheck_performance_path = tables / "seqcheck_performance.csv"
    fastqc_runs_path = tables / "fastqc_runs.csv"
    fastqc_modules_path = tables / "fastqc_modules.csv"
    fastqc_performance_path = tables / "fastqc_performance.csv"
    runtime.write_csv(seqspec_path, seqspec_rows, list(SEQSPEC_RUN_FIELDS))
    runtime.write_csv(
        seqcheck_runs_path,
        seqcheck_runs,
        ["condition_source", *seqcheck_runner.RUN_FIELDS],
    )
    runtime.write_csv(
        seqcheck_metrics_path, seqcheck_metrics, list(seqcheck_runner.METRIC_FIELDS)
    )
    runtime.write_csv(
        seqcheck_assessments_path,
        seqcheck_assessments,
        list(seqcheck_runner.ASSESSMENT_FIELDS),
    )
    runtime.write_csv(
        seqcheck_performance_path,
        seqcheck_performance,
        list(seqcheck_runner.PERFORMANCE_FIELDS),
    )
    runtime.write_csv(fastqc_runs_path, fastqc_runs, list(FASTQC_RUN_FIELDS))
    runtime.write_csv(fastqc_modules_path, fastqc_modules, list(FASTQC_MODULE_FIELDS))
    runtime.write_csv(
        fastqc_performance_path,
        fastqc_performance,
        list(FASTQC_PERFORMANCE_FIELDS),
    )
    validation = validate_outputs(
        execution_id=execution_id,
        combined=combined,
        seqspec_runs=seqspec_rows,
        seqcheck_runs=seqcheck_runs,
        seqcheck_metrics=seqcheck_metrics,
        seqcheck_assessments=seqcheck_assessments,
        fastqc_runs=fastqc_runs,
        fastqc_modules=fastqc_modules,
        fastqc_performance=fastqc_performance,
    )
    validation_path = output_root / "validation" / "execution.json"
    runtime.write_json(validation_path, validation)
    if not validation["valid"]:
        raise ValueError(
            "complementarity execution failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )
    manifest_path = output_root / "manifests" / "execution.json"
    manifest = {
        **stable_identity,
        "execution_id": execution_id,
        "selection_id": complementarity["selection_id"],
        "cohort_split": complementarity["cohort_split"],
        "created_at": utc_now(),
        "valid": True,
        "inputs": {
            "materialization": runtime.file_identity(materialization_path),
            "protocol": runtime.file_identity(protocol_path),
        },
        "tools": tools,
        "counts": validation["counts"],
        "outputs": {
            "seqspec_runs": runtime.file_identity(seqspec_path),
            "seqcheck_runs": runtime.file_identity(seqcheck_runs_path),
            "seqcheck_metrics": runtime.file_identity(seqcheck_metrics_path),
            "seqcheck_assessments": runtime.file_identity(seqcheck_assessments_path),
            "seqcheck_performance": runtime.file_identity(seqcheck_performance_path),
            "fastqc_runs": runtime.file_identity(fastqc_runs_path),
            "fastqc_modules": runtime.file_identity(fastqc_modules_path),
            "fastqc_performance": runtime.file_identity(fastqc_performance_path),
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_materialization(
    path: Path,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    value = runtime.load_json(path)
    if value.get("valid") is not True:
        raise ValueError("complementarity materialization is not valid")
    stable = {
        field: value.get(field)
        for field in materializer.MATERIALIZATION_IDENTITY_FIELDS
    }
    if any(stable[field] is None and field != "quality_policy" for field in stable):
        raise ValueError("complementarity materialization identity is incomplete")
    if runtime.sha256_json(stable)[:16] != value.get("complementarity_id"):
        raise ValueError("complementarity materialization id is not content-addressed")
    base_path = materializer.verified_path(
        value.get("inputs", {}).get("base_materialization", {}),
        "base materialization",
    )
    if runtime.file_sha256(base_path) != value["base_materialization_sha256"]:
        raise ValueError("base materialization identities differ")
    base, base_conditions = seqcheck_runner.load_materialization(base_path)
    if base["materialization_id"] != value["base_materialization_id"]:
        raise ValueError("base materialization id differs")
    controls_path = materializer.verified_path(
        value.get("outputs", {}).get("controls", {}), "control manifest"
    )
    payload = runtime.load_json(controls_path)
    controls = payload.get("controls")
    if (
        payload.get("complementarity_id") != value["complementarity_id"]
        or not isinstance(controls, list)
        or not controls
    ):
        raise ValueError("control manifest does not reconcile")
    if [materializer.control_identity(item) for item in controls] != value["controls"]:
        raise ValueError("control records differ from materialization identity")
    for control in controls:
        if control.get("complementarity_id") != value["complementarity_id"]:
            raise ValueError("control complementarity id differs")
        if runtime.sha256_json(materializer.control_identity(control))[
            :16
        ] != control.get("condition_id"):
            raise ValueError("control condition id is not content-addressed")
        identities = [
            control["spec"],
            *control["inputs"],
            *base_materializer.nested_file_identities(control["mutation"]),
        ]
        for identity in identities:
            materializer.verified_path(identity, "control artifact")
    if value.get("counts", {}).get("controls") != len(controls):
        raise ValueError("control count does not reconcile")
    return value, base, base_conditions, controls


def run_seqspec_condition(
    *,
    execution_id: str,
    source: str,
    condition: dict[str, Any],
    seqspec_bin: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    condition_id = condition["condition_id"]
    log_root = output_root / "logs" / "seqspec" / condition_id
    argv = [str(seqspec_bin), "check", condition["spec"]["path"]]
    measurement = runtime.run_measured_command(
        argv,
        stdout_path=log_root / "stdout.txt",
        stderr_path=log_root / "stderr.txt",
        timeout_seconds=timeout_seconds,
    )
    expected = "failure" if condition["expected_seqspec_check"] == "failure" else "pass"
    observed = (
        "pass"
        if measurement["exit_code"] == 0 and not measurement["timed_out"]
        else "failure"
    )
    if measurement["timed_out"] or observed != expected:
        raise ValueError(
            f"{condition_id}: expected seqspec {expected}, observed {observed}"
        )
    return {
        "execution_id": execution_id,
        "condition_source": source,
        "condition_id": condition_id,
        "configuration_accession": condition["configuration_accession"],
        "condition_kind": condition["condition_kind"],
        "operator_id": condition["operator_id"],
        "variant": condition["variant"],
        "expected_status": expected,
        "observed_status": observed,
        "runtime_schema_version": measurement["runtime_schema_version"],
        "exit_code": measurement["exit_code"],
        "timed_out": measurement["timed_out"],
        "wall_time_seconds": measurement["wall_time_seconds"],
        "user_cpu_seconds": measurement["user_cpu_seconds"],
        "system_cpu_seconds": measurement["system_cpu_seconds"],
        "peak_resident_memory_bytes": measurement["peak_resident_memory_bytes"],
        "command_json": runtime.canonical_json(argv),
        "stdout_path": measurement["stdout"]["path"],
        "stdout_sha256": measurement["stdout"]["sha256"],
        "stderr_path": measurement["stderr"]["path"],
        "stderr_sha256": measurement["stderr"]["sha256"],
    }


def adapt_seqcheck_condition(
    condition: dict[str, Any], *, source: str
) -> dict[str, Any]:
    if source == "structural":
        return condition
    result = copy.deepcopy(condition)
    result["expected"] = {
        **result["expected"],
        "process_outcome": result["expected"]["seqcheck_process"],
    }
    return result


def run_fastqc_input(
    *,
    execution_id: str,
    source: str,
    condition: dict[str, Any],
    input_value: dict[str, Any],
    protocol: dict[str, Any],
    tools: dict[str, Any],
    fastqc_bin: Path,
    fastqc_limits: Path | None,
    output_root: Path,
    timeout_seconds: int,
    cache: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    cache_key = runtime.sha256_json(
        {
            "sha256": input_value["sha256"],
            "size_bytes": input_value["size_bytes"],
            "basename": Path(input_value["path"]).name,
            "fastqc": functional_fastqc_identity(tools),
        }
    )[:16]
    cached = cache.get(cache_key)
    performance = None
    command_json = ""
    if cached is None:
        report_root = output_root / "reports" / "fastqc" / cache_key
        report_root.mkdir(parents=True)
        argv = [
            str(fastqc_bin),
            "--extract",
            "--threads",
            str(protocol["fastqc"]["threads"]),
            "--outdir",
            str(report_root),
        ]
        if fastqc_limits is not None:
            argv.extend(["--limits", str(fastqc_limits)])
        argv.append(input_value["path"])
        measurement = runtime.run_measured_command(
            argv,
            stdout_path=output_root / "logs" / "fastqc" / cache_key / "stdout.txt",
            stderr_path=output_root / "logs" / "fastqc" / cache_key / "stderr.txt",
            timeout_seconds=timeout_seconds,
        )
        runtime.require_success(measurement)
        parsed = fastqc_parser.parse_fastqc_output(
            report_root / fastqc_parser.output_stem(Path(input_value["path"]))
        )
        data_identity = runtime.file_identity(Path(parsed["data_path"]))
        summary_identity = runtime.file_identity(Path(parsed["summary_path"]))
        cached = {
            "source_condition_id": condition["condition_id"],
            "parsed": parsed,
            "data": data_identity,
            "summary": summary_identity,
        }
        cache[cache_key] = cached
        command_json = runtime.canonical_json(argv)
        performance = fastqc_performance_row(
            execution_id=execution_id,
            cache_key=cache_key,
            condition=condition,
            input_value=input_value,
            measurement=measurement,
        )
    parsed = cached["parsed"]
    context = fastqc_context(
        execution_id=execution_id,
        source=source,
        condition=condition,
        input_value=input_value,
    )
    cache_hit = (
        cached["source_condition_id"] != condition["condition_id"] or not command_json
    )
    run = {
        **context,
        "input_path": input_value["path"],
        "input_sha256": input_value["sha256"],
        "cache_key": cache_key,
        "cache_hit": cache_hit,
        "cache_source_condition_id": cached["source_condition_id"],
        "fastqc_version": parsed["fastqc_version"],
        "module_count": len(parsed["modules"]),
        "command_json": command_json,
        "data_path": cached["data"]["path"],
        "data_sha256": cached["data"]["sha256"],
        "summary_path": cached["summary"]["path"],
        "summary_sha256": cached["summary"]["sha256"],
    }
    modules = [
        {
            **context,
            "module_name": module["name"],
            "status": module["status"],
            "status_rank": fastqc_parser.STATUS_RANK[module["status"]],
            "data_path": cached["data"]["path"],
            "data_sha256": cached["data"]["sha256"],
            "summary_path": cached["summary"]["path"],
            "summary_sha256": cached["summary"]["sha256"],
        }
        for module in parsed["modules"]
    ]
    return {"run": run, "modules": modules, "performance": performance}


def fastqc_context(
    *,
    execution_id: str,
    source: str,
    condition: dict[str, Any],
    input_value: dict[str, Any],
) -> dict[str, Any]:
    return {
        "execution_id": execution_id,
        "condition_source": source,
        "condition_id": condition["condition_id"],
        "configuration_accession": condition["configuration_accession"],
        "condition_kind": condition["condition_kind"],
        "operator_id": condition["operator_id"],
        "variant": condition["variant"],
        "case_id": input_value["case_id"],
        "fastq_accession": input_value["fastq_accession"],
        "read_id": input_value["read_id"],
    }


def fastqc_performance_row(
    *,
    execution_id: str,
    cache_key: str,
    condition: dict[str, Any],
    input_value: dict[str, Any],
    measurement: dict[str, Any],
) -> dict[str, Any]:
    return {
        "execution_id": execution_id,
        "cache_key": cache_key,
        "source_condition_id": condition["condition_id"],
        "case_id": input_value["case_id"],
        "runtime_schema_version": measurement["runtime_schema_version"],
        "command_json": runtime.canonical_json(measurement["argv"]),
        "exit_code": measurement["exit_code"],
        "timed_out": measurement["timed_out"],
        "wall_time_seconds": measurement["wall_time_seconds"],
        "user_cpu_seconds": measurement["user_cpu_seconds"],
        "system_cpu_seconds": measurement["system_cpu_seconds"],
        "peak_resident_memory_bytes": measurement["peak_resident_memory_bytes"],
        "records_processed": input_value["records"],
        "compressed_bytes_read": input_value["size_bytes"],
        "stdout_path": measurement["stdout"]["path"],
        "stdout_sha256": measurement["stdout"]["sha256"],
        "stderr_path": measurement["stderr"]["path"],
        "stderr_sha256": measurement["stderr"]["sha256"],
    }


def validate_protocol(protocol: dict[str, Any]) -> None:
    materializer.validate_protocol(protocol)
    fastqc = protocol.get("fastqc", {})
    if fastqc.get("extract") is not True or fastqc.get("threads") != 1:
        raise ValueError("FastQC execution must extract output with one thread")
    if fastqc.get("status_order") != list(fastqc_parser.VALID_STATUSES):
        raise ValueError("FastQC status order differs from parser semantics")
    seqcheck = protocol.get("seqcheck", {})
    if seqcheck.get("n_reads") != 0 or seqcheck.get("format") != "json":
        raise ValueError("seqcheck must inspect all materialized reads as JSON")


def validate_outputs(
    *,
    execution_id: str,
    combined: list[tuple[str, dict[str, Any]]],
    seqspec_runs: list[dict[str, Any]],
    seqcheck_runs: list[dict[str, Any]],
    seqcheck_metrics: list[dict[str, Any]],
    seqcheck_assessments: list[dict[str, Any]],
    fastqc_runs: list[dict[str, Any]],
    fastqc_modules: list[dict[str, Any]],
    fastqc_performance: list[dict[str, Any]],
) -> dict[str, Any]:
    errors = []
    condition_ids = Counter(value["condition_id"] for _, value in combined)
    if Counter(value["condition_id"] for value in seqspec_runs) != condition_ids:
        errors.append("seqspec runs do not reconcile with conditions")
    if Counter(value["condition_id"] for value in seqcheck_runs) != condition_ids:
        errors.append("seqcheck runs do not reconcile with conditions")
    expected_fastqc = Counter(
        value["condition_id"] for _, value in combined for _ in value["inputs"]
    )
    if Counter(value["condition_id"] for value in fastqc_runs) != expected_fastqc:
        errors.append("FastQC runs do not reconcile with condition inputs")
    module_counts = Counter(value["condition_id"] for value in fastqc_modules)
    expected_modules = Counter()
    for run in fastqc_runs:
        expected_modules[run["condition_id"]] += run["module_count"]
    if module_counts != expected_modules:
        errors.append("FastQC modules do not reconcile with runs")
    if any(
        value["execution_id"] != execution_id for value in seqspec_runs + fastqc_runs
    ):
        errors.append("tool execution identifiers differ")
    for label, rows in (
        ("seqcheck run", seqcheck_runs),
        ("seqcheck metric", seqcheck_metrics),
        ("seqcheck assessment", seqcheck_assessments),
        ("FastQC module", fastqc_modules),
    ):
        if any(value["execution_id"] != execution_id for value in rows):
            errors.append(f"{label} execution identifiers differ")
    if any(
        value["expected_status"] != value["observed_status"] for value in seqspec_runs
    ):
        errors.append("seqspec outcomes differ from expected outcomes")
    if any(not value["process_outcome_matches"] for value in seqcheck_runs):
        errors.append("seqcheck outcomes differ from expected outcomes")
    metric_counts = Counter(value["condition_id"] for value in seqcheck_metrics)
    assessment_counts = Counter(value["condition_id"] for value in seqcheck_assessments)
    for run in seqcheck_runs:
        if run["metric_count"] != metric_counts[run["condition_id"]]:
            errors.append(f"{run['condition_id']}: seqcheck metric count differs")
        if run["assessment_count"] != assessment_counts[run["condition_id"]]:
            errors.append(f"{run['condition_id']}: seqcheck assessment count differs")
        for field in ("stdout", "stderr"):
            if (
                runtime.file_sha256(Path(run[f"{field}_path"]))
                != run[f"{field}_sha256"]
            ):
                errors.append(f"{run['condition_id']}: seqcheck {field} hash changed")
        if (
            run["observed_process_outcome"] == "success"
            and runtime.file_sha256(Path(run["report_path"])) != run["report_sha256"]
        ):
            errors.append(f"{run['condition_id']}: seqcheck report hash changed")
    cache_keys = [value["cache_key"] for value in fastqc_performance]
    if len(cache_keys) != len(set(cache_keys)):
        errors.append("FastQC cache keys were executed more than once")
    if set(cache_keys) != {value["cache_key"] for value in fastqc_runs}:
        errors.append("FastQC cache executions do not cover every run")
    for row in fastqc_runs:
        for field in ("data", "summary"):
            if (
                runtime.file_sha256(Path(row[f"{field}_path"]))
                != row[f"{field}_sha256"]
            ):
                errors.append(f"{row['condition_id']}: FastQC {field} hash changed")
    for row in seqspec_runs:
        for field in ("stdout", "stderr"):
            if (
                runtime.file_sha256(Path(row[f"{field}_path"]))
                != row[f"{field}_sha256"]
            ):
                errors.append(f"{row['condition_id']}: seqspec {field} hash changed")
    for row in fastqc_performance:
        for field in ("stdout", "stderr"):
            if (
                runtime.file_sha256(Path(row[f"{field}_path"]))
                != row[f"{field}_sha256"]
            ):
                errors.append(
                    f"{row['source_condition_id']}: FastQC {field} hash changed"
                )
    return {
        "schema_version": SCHEMA_VERSION,
        "execution_id": execution_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "counts": {
            "conditions": len(combined),
            "structural_conditions": sum(
                source == "structural" for source, _ in combined
            ),
            "control_conditions": sum(source == "control" for source, _ in combined),
            "seqspec_runs": len(seqspec_runs),
            "seqcheck_runs": len(seqcheck_runs),
            "seqcheck_metrics": len(seqcheck_metrics),
            "seqcheck_assessments": len(seqcheck_assessments),
            "fastqc_input_runs": len(fastqc_runs),
            "fastqc_cache_executions": len(fastqc_performance),
            "fastqc_cache_hits": sum(value["cache_hit"] for value in fastqc_runs),
            "fastqc_modules": len(fastqc_modules),
        },
    }


def tool_identities(
    *,
    seqcheck_bin: Path,
    seqspec_bin: Path,
    fastqc_bin: Path,
    fastqc_limits: Path | None,
    timeout_seconds: int,
) -> dict[str, Any]:
    result = {
        "seqcheck": runtime.executable_identity(
            seqcheck_bin, timeout_seconds=min(timeout_seconds, 30)
        ),
        "seqspec": runtime.executable_identity(
            seqspec_bin, timeout_seconds=min(timeout_seconds, 30)
        ),
        "fastqc": runtime.executable_identity(
            fastqc_bin, timeout_seconds=min(timeout_seconds, 30)
        ),
        "fastqc_limits": (
            runtime.file_identity(fastqc_limits) if fastqc_limits else None
        ),
        "runner": runtime.script_identity(
            Path(__file__).resolve(), version=RUNNER_VERSION
        ),
        "fastqc_parser": runtime.script_identity(
            Path(fastqc_parser.__file__).resolve(),
            version=fastqc_parser.RUNTIME_VERSION,
        ),
        "seqcheck_runner": runtime.script_identity(
            Path(seqcheck_runner.__file__).resolve(),
            version=seqcheck_runner.RUNNER_VERSION,
        ),
        "runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
    }
    return result


def functional_fastqc_identity(tools: dict[str, Any]) -> dict[str, Any]:
    return {
        "executable": runtime.functional_executable_identity(tools["fastqc"]),
        "limits": (
            {
                "sha256": tools["fastqc_limits"]["sha256"],
                "size_bytes": tools["fastqc_limits"]["size_bytes"],
            }
            if tools["fastqc_limits"]
            else None
        ),
    }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
