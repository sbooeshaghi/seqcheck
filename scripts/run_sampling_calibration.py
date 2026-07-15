#!/usr/bin/env python3
"""Run the seqcheck sampling calibration from a frozen case selection."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.parse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
RUNNER_VERSION = "0.1.0"
SAMPLER_VERSION = "0.1.0"
SAMPLER_CORE_VERSION = "0.1.0"
SELECTION_IDENTITY_FIELDS = (
    "schema_version",
    "freeze_id",
    "cohort_sha256",
    "sampling_protocol_sha256",
    "selector",
    "seqspec",
    "cases",
)
REPORT_FIELDS = [
    "study_run_id",
    "selection_id",
    "case_id",
    "family_id",
    "configuration_accession",
    "fastq_accession",
    "selection_role",
    "condition_id",
    "sampling_method",
    "requested_records_per_fastq",
    "sampling_seed",
    "records_processed",
    "report_schema_version",
    "report_path",
    "report_sha256",
    "result_sha256",
    "metric_count",
]
PERFORMANCE_FIELDS = [
    "study_run_id",
    "selection_id",
    "case_id",
    "condition_id",
    "sampling_method",
    "requested_records_per_fastq",
    "sampling_seed",
    "execution_scope",
    "replicate",
    "is_warmup",
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
    "report_path",
    "report_sha256",
    "stdout_path",
    "stdout_sha256",
    "stderr_path",
    "stderr_sha256",
]
METRIC_FIELDS = [
    "study_run_id",
    "selection_id",
    "case_id",
    "family_id",
    "configuration_accession",
    "fastq_accession",
    "selection_role",
    "condition_id",
    "sampling_method",
    "requested_records_per_fastq",
    "sampling_seed",
    "report_path",
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
]
ERROR_FIELDS = [
    "study_run_id",
    "selection_id",
    "case_id",
    "family_id",
    "configuration_accession",
    "fastq_accession",
    "selection_role",
    "condition_id",
    "sampling_method",
    "requested_records_per_fastq",
    "sampling_seed",
    "check",
    "files",
    "reads",
    "regions",
    "ontology_terms",
    "metric_id",
    "metric_name",
    "unit",
    "complete_value",
    "sample_value",
    "signed_error",
    "absolute_error",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run and reconcile all seqcheck sampling calibration conditions."
    )
    parser.add_argument("--case-selection-manifest", required=True, type=Path)
    parser.add_argument("--sampling-protocol", required=True, type=Path)
    parser.add_argument("--seqcheck-bin", required=True, type=Path)
    parser.add_argument("--sampler-script", type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=3600)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sampler_script = args.sampler_script or (
        Path(__file__).resolve().parent / "sample_fastq_matrix.py"
    )
    try:
        manifest = run_sampling_calibration(
            case_selection_manifest_path=args.case_selection_manifest.resolve(),
            protocol_path=args.sampling_protocol.resolve(),
            seqcheck_bin=args.seqcheck_bin.resolve(),
            sampler_script=sampler_script.resolve(),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
        )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"run_sampling_calibration: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def run_sampling_calibration(
    *,
    case_selection_manifest_path: Path,
    protocol_path: Path,
    seqcheck_bin: Path,
    sampler_script: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    selection_manifest, protocol, cases = load_and_validate_inputs(
        case_selection_manifest_path,
        protocol_path,
        seqcheck_bin,
        sampler_script,
    )
    output_root.mkdir(parents=True)
    sampler_core = sampler_script.with_name("sample_fastq.py")
    tools = {
        "seqcheck": executable_identity(seqcheck_bin),
        "sampler": script_identity(sampler_script, version=SAMPLER_VERSION),
        "sampler_core": script_identity(sampler_core, version=SAMPLER_CORE_VERSION),
        "runtime": script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
        "runner": script_identity(Path(__file__).resolve(), version=RUNNER_VERSION),
    }
    study_identity = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": selection_manifest["selection_id"],
        "sampling_protocol_sha256": runtime.file_sha256(protocol_path),
        "seqcheck": functional_executable_identity(tools["seqcheck"]),
        "sampler": functional_script_identity(tools["sampler"]),
        "sampler_core": functional_script_identity(tools["sampler_core"]),
        "runtime": functional_script_identity(tools["runtime"]),
        "runner": functional_script_identity(tools["runner"]),
    }
    study_run_id = runtime.sha256_json(study_identity)[:16]
    report_rows = []
    metric_rows = []
    performance_rows = []
    matrix_records = []
    repeat_mismatches = []
    for case in cases:
        result = run_case(
            case=case,
            protocol=protocol,
            study_run_id=study_run_id,
            selection_id=selection_manifest["selection_id"],
            seqcheck_bin=seqcheck_bin,
            sampler_script=sampler_script,
            output_root=output_root,
            timeout_seconds=timeout_seconds,
        )
        report_rows.extend(result["report_rows"])
        metric_rows.extend(result["metric_rows"])
        performance_rows.extend(result["performance_rows"])
        matrix_records.append(result["matrix_record"])
        repeat_mismatches.extend(result["repeat_mismatches"])

    tables_dir = output_root / "tables"
    reports_path = tables_dir / "reports.csv"
    metrics_path = tables_dir / "metrics.csv"
    errors_path = tables_dir / "sampling_errors.csv"
    performance_path = tables_dir / "performance.csv"
    error_rows, error_issues = build_sampling_errors(metric_rows)
    runtime.write_csv(reports_path, report_rows, REPORT_FIELDS)
    runtime.write_csv(metrics_path, metric_rows, METRIC_FIELDS)
    runtime.write_csv(errors_path, error_rows, ERROR_FIELDS)
    runtime.write_csv(performance_path, performance_rows, PERFORMANCE_FIELDS)
    validation = validate_outputs(
        study_run_id=study_run_id,
        selection_id=selection_manifest["selection_id"],
        protocol=protocol,
        cases=cases,
        output_root=output_root,
        report_rows=report_rows,
        metric_rows=metric_rows,
        error_rows=error_rows,
        error_issues=error_issues,
        performance_rows=performance_rows,
        matrix_records=matrix_records,
        repeat_mismatches=repeat_mismatches,
    )
    validation_path = output_root / "validation" / "sampling_calibration.json"
    runtime.write_json(validation_path, validation)
    if not validation["valid"]:
        raise ValueError(f"sampling calibration validation failed: {validation_path}")

    manifest_path = output_root / "manifests" / "study.json"
    manifest = {
        **study_identity,
        "study_run_id": study_run_id,
        "created_at": utc_now(),
        "valid": True,
        "inputs": {
            "case_selection_manifest": runtime.file_identity(
                case_selection_manifest_path
            ),
            "sampling_protocol": runtime.file_identity(protocol_path),
        },
        "tools": tools,
        "counts": validation["counts"],
        "matrices": matrix_records,
        "outputs": {
            "reports": runtime.file_identity(reports_path),
            "metrics": runtime.file_identity(metrics_path),
            "sampling_errors": runtime.file_identity(errors_path),
            "performance": runtime.file_identity(performance_path),
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_and_validate_inputs(
    selection_manifest_path: Path,
    protocol_path: Path,
    seqcheck_bin: Path,
    sampler_script: Path,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    selection = runtime.load_json(selection_manifest_path)
    protocol = runtime.load_json(protocol_path)
    if selection.get("valid") is not True:
        raise ValueError("sampling case selection is not valid")
    if selection.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling case selection does not use schema 0.1.0")
    selection_id = str(selection.get("selection_id", "")).strip()
    if not selection_id:
        raise ValueError("sampling case selection has no selection id")
    if selection.get("sampling_protocol_sha256") != runtime.file_sha256(protocol_path):
        raise ValueError("sampling protocol hash differs from case selection")
    cases_identity = selection.get("outputs", {}).get("cases_json", {})
    if not isinstance(cases_identity, dict):
        raise ValueError("case selection has no cases JSON identity")
    cases_path = Path(str(cases_identity.get("path", ""))).resolve()
    if runtime.file_sha256(cases_path) != cases_identity.get("sha256"):
        raise ValueError("sampling cases JSON hash changed")
    cases_payload = runtime.load_json(cases_path)
    if cases_payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling cases do not use schema 0.1.0")
    if cases_payload.get("selection_id") != selection_id:
        raise ValueError("sampling cases selection id does not match manifest")
    if cases_payload.get("freeze_id") != selection.get("freeze_id"):
        raise ValueError("sampling cases freeze id does not match manifest")
    cases = cases_payload.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("sampling cases must be a non-empty list")
    if selection.get("counts", {}).get("fastq_cases") != len(cases):
        raise ValueError("sampling case count does not reconcile")
    validate_selection_identity(selection, cases)

    validate_protocol(protocol)
    if not seqcheck_bin.is_file():
        raise ValueError(f"seqcheck executable does not exist: {seqcheck_bin}")
    if not sampler_script.is_file():
        raise ValueError(f"sampler script does not exist: {sampler_script}")
    sampler_core = sampler_script.with_name("sample_fastq.py")
    if not sampler_core.is_file():
        raise ValueError(f"sampler core does not exist: {sampler_core}")
    case_ids = []
    for case in cases:
        if not isinstance(case, dict):
            raise ValueError("sampling case entries must be objects")
        case_id_value = required_string(case, "case_id")
        case_ids.append(case_id_value)
        if case.get("selection_id") != selection_id:
            raise ValueError(f"{case_id_value}: selection id differs")
        spec_path = Path(required_string(case, "spec_path")).resolve()
        if runtime.file_sha256(spec_path) != required_string(case, "spec_sha256"):
            raise ValueError(f"{case_id_value}: spec hash changed")
        source = required_string(case, "fastq_url")
        if not is_remote_source(source) and not Path(source).is_file():
            raise ValueError(f"{case_id_value}: local FASTQ does not exist: {source}")
        for field in (
            "family_id",
            "configuration_accession",
            "fastq_accession",
            "read_id",
            "modality",
            "selection_role",
        ):
            required_string(case, field)
        if int(case.get("declared_compressed_bytes", 0)) <= 0:
            raise ValueError(f"{case_id_value}: declared FASTQ size must be positive")
    if len(case_ids) != len(set(case_ids)):
        raise ValueError("sampling case identifiers are not unique")
    return selection, protocol, cases


def validate_selection_identity(
    selection: dict[str, Any], cases: list[dict[str, Any]]
) -> None:
    stable_selection = {
        field: selection.get(field) for field in SELECTION_IDENTITY_FIELDS
    }
    if any(stable_selection[field] is None for field in SELECTION_IDENTITY_FIELDS):
        raise ValueError("sampling case selection identity is incomplete")
    expected_id = runtime.sha256_json(stable_selection)[:16]
    if selection["selection_id"] != expected_id:
        raise ValueError("sampling case selection id is not content-addressed")

    stable_cases = [
        {
            key: value
            for key, value in case.items()
            if key not in {"selection_id", "case_id", "spec_path"}
        }
        for case in cases
    ]
    if stable_cases != stable_selection["cases"]:
        raise ValueError("sampling case rows differ from content-addressed selection")


def validate_protocol(protocol: dict[str, Any]) -> None:
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling protocol does not use schema 0.1.0")
    sampling = protocol.get("sampling", {})
    performance = protocol.get("performance", {})
    sizes = sampling.get("records_per_fastq")
    seeds = sampling.get("reservoir_seeds")
    if (
        not isinstance(sizes, list)
        or not sizes
        or any(type(value) is not int or value <= 0 for value in sizes)
        or len(sizes) != len(set(sizes))
    ):
        raise ValueError("sampling protocol needs unique positive record counts")
    if (
        not isinstance(seeds, list)
        or not seeds
        or any(type(value) is not int or value < 0 for value in seeds)
        or len(seeds) != len(set(seeds))
    ):
        raise ValueError("sampling protocol needs unique nonnegative reservoir seeds")
    if sampling.get("include_prefix") is not True:
        raise ValueError("sampling protocol must include prefix conditions")
    if sampling.get("source_traversals_for_bounded_samples_per_fastq") != 1:
        raise ValueError("bounded samples must use one source traversal")
    if sampling.get("complete_stream_reference") is not True:
        raise ValueError("sampling protocol must include a complete reference")
    warmups = performance.get("warmup_runs")
    measured = performance.get("measured_runs")
    if type(warmups) is not int or warmups < 0:
        raise ValueError("performance warm-up count must be nonnegative")
    if type(measured) is not int or measured <= 0:
        raise ValueError("performance measured count must be positive")
    expected = protocol.get("expected_conditions_per_fastq", {})
    if any(
        type(expected.get(field)) is not int
        for field in ("prefix", "reservoir", "complete_stream", "total")
    ):
        raise ValueError("expected condition counts must be integers")
    if expected.get("prefix") != len(sizes):
        raise ValueError("expected prefix condition count does not reconcile")
    if expected.get("reservoir") != len(sizes) * len(seeds):
        raise ValueError("expected reservoir condition count does not reconcile")
    if expected.get("complete_stream") != 1:
        raise ValueError("expected complete-stream condition count must be one")
    if expected.get("total") != len(sizes) * (len(seeds) + 1) + 1:
        raise ValueError("expected total condition count does not reconcile")


def run_case(
    *,
    case: dict[str, Any],
    protocol: dict[str, Any],
    study_run_id: str,
    selection_id: str,
    seqcheck_bin: Path,
    sampler_script: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    case_id_value = case["case_id"]
    sample_root = output_root / "samples" / case_id_value
    sampling = protocol["sampling"]
    sampler_argv = [
        sys.executable,
        str(sampler_script),
        "--input",
        case["fastq_url"],
        "--output-root",
        str(sample_root),
        "--include-prefix",
        "--configuration-accession",
        case["configuration_accession"],
        "--modality",
        case["modality"],
        "--fastq-accession",
        case["fastq_accession"],
        "--read-id",
        case["read_id"],
    ]
    for size in sampling["records_per_fastq"]:
        sampler_argv.extend(["--n-reads", str(size)])
    for seed in sampling["reservoir_seeds"]:
        sampler_argv.extend(["--seed", str(seed)])
    sampler_measurement = runtime.run_measured_command(
        sampler_argv,
        stdout_path=output_root / "logs" / "sampler" / f"{case_id_value}.stdout",
        stderr_path=output_root / "logs" / "sampler" / f"{case_id_value}.stderr",
        timeout_seconds=timeout_seconds,
    )
    runtime.require_success(sampler_measurement)
    manifest_value = Path(sampler_measurement["stdout"]["path"]).read_text().strip()
    matrix_manifest_path = Path(manifest_value).resolve()
    matrix = runtime.load_json(matrix_manifest_path)
    validate_matrix(matrix, protocol, case)
    source_records = int(matrix["source_records"][0]["records_streamed"])
    matrix_record = {
        "case_id": case_id_value,
        "matrix_id": matrix["matrix_id"],
        "manifest": runtime.file_identity(matrix_manifest_path),
        "source_records": source_records,
        "source_uncompressed_sha256": matrix["source_records"][0][
            "source_uncompressed_sha256"
        ],
        "bounded_conditions": len(matrix["conditions"]),
    }
    performance_rows = [
        performance_row(
            measurement=sampler_measurement,
            study_run_id=study_run_id,
            selection_id=selection_id,
            case_id=case_id_value,
            condition_id="sample-matrix",
            method="matrix",
            n_reads=0,
            seed=None,
            scope=(
                "remote_transfer" if is_remote_source(case["fastq_url"]) else "sampling"
            ),
            replicate=1,
            is_warmup=False,
            records_processed=source_records,
            compressed_bytes=source_compressed_bytes(case, matrix),
            report_path=matrix_manifest_path,
        )
    ]
    report_rows = []
    metric_rows = []
    repeat_mismatches = []

    complete = run_report_condition(
        case=case,
        condition_id="complete",
        method="complete",
        n_reads=source_records,
        seed=None,
        input_path=case["fastq_url"],
        execution_scope=(
            "remote_transfer"
            if is_remote_source(case["fastq_url"])
            else "source_stream"
        ),
        warmup_runs=0,
        measured_runs=1,
        records_processed=source_records,
        compressed_bytes=source_compressed_bytes(case, matrix),
        study_run_id=study_run_id,
        selection_id=selection_id,
        seqcheck_bin=seqcheck_bin,
        output_root=output_root,
        timeout_seconds=timeout_seconds,
    )
    report_rows.append(complete["report_row"])
    metric_rows.extend(complete["metric_rows"])
    performance_rows.extend(complete["performance_rows"])
    repeat_mismatches.extend(complete["repeat_mismatches"])

    for condition in matrix["conditions"]:
        output = condition["outputs"][0]
        bounded = run_report_condition(
            case=case,
            condition_id=condition["condition_id"],
            method=condition["method"],
            n_reads=int(condition["requested_records_per_fastq"]),
            seed=condition["seed"],
            input_path=output["output_path"],
            execution_scope="local_compute",
            warmup_runs=int(protocol["performance"]["warmup_runs"]),
            measured_runs=int(protocol["performance"]["measured_runs"]),
            records_processed=int(output["records_selected"]),
            compressed_bytes=Path(output["output_path"]).stat().st_size,
            study_run_id=study_run_id,
            selection_id=selection_id,
            seqcheck_bin=seqcheck_bin,
            output_root=output_root,
            timeout_seconds=timeout_seconds,
        )
        report_rows.append(bounded["report_row"])
        metric_rows.extend(bounded["metric_rows"])
        performance_rows.extend(bounded["performance_rows"])
        repeat_mismatches.extend(bounded["repeat_mismatches"])
    return {
        "matrix_record": matrix_record,
        "report_rows": report_rows,
        "metric_rows": metric_rows,
        "performance_rows": performance_rows,
        "repeat_mismatches": repeat_mismatches,
    }


def run_report_condition(
    *,
    case: dict[str, Any],
    condition_id: str,
    method: str,
    n_reads: int,
    seed: int | None,
    input_path: str,
    execution_scope: str,
    warmup_runs: int,
    measured_runs: int,
    records_processed: int,
    compressed_bytes: int,
    study_run_id: str,
    selection_id: str,
    seqcheck_bin: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    canonical_path = output_root / "reports" / case["case_id"] / f"{condition_id}.json"
    runs = [("warmup", replicate, True) for replicate in range(1, warmup_runs + 1)] + [
        ("measured", replicate, False) for replicate in range(1, measured_runs + 1)
    ]
    performance_rows = []
    repeat_mismatches = []
    canonical_payload = None
    canonical_result_sha256 = ""
    for run_label, replicate, is_warmup in runs:
        if run_label == "measured" and replicate == 1:
            report_path = canonical_path
        else:
            report_path = (
                output_root
                / "performance_reports"
                / case["case_id"]
                / condition_id
                / f"{run_label}-{replicate:02d}.json"
            )
        report_path.parent.mkdir(parents=True, exist_ok=True)
        log_root = output_root / "logs" / "seqcheck" / case["case_id"] / condition_id
        measurement = runtime.run_measured_command(
            [
                str(seqcheck_bin),
                "check",
                "--spec",
                case["spec_path"],
                "--modality",
                case["modality"],
                "--n-reads",
                "0",
                "--format",
                "json",
                "--output",
                str(report_path),
                input_path,
            ],
            stdout_path=log_root / f"{run_label}-{replicate:02d}.stdout",
            stderr_path=log_root / f"{run_label}-{replicate:02d}.stderr",
            timeout_seconds=timeout_seconds,
        )
        runtime.require_success(measurement)
        payload = runtime.load_json(report_path)
        report_schema = str(payload.get("report_schema_version", "")).strip()
        if not report_schema:
            raise ValueError(f"seqcheck report has no schema version: {report_path}")
        result_sha256 = runtime.sha256_json(payload.get("results", []))
        if run_label == "measured":
            if canonical_payload is None:
                canonical_payload = payload
                canonical_result_sha256 = result_sha256
            elif result_sha256 != canonical_result_sha256:
                repeat_mismatches.append(
                    {
                        "case_id": case["case_id"],
                        "condition_id": condition_id,
                        "replicate": replicate,
                        "expected_result_sha256": canonical_result_sha256,
                        "observed_result_sha256": result_sha256,
                    }
                )
        performance_rows.append(
            performance_row(
                measurement=measurement,
                study_run_id=study_run_id,
                selection_id=selection_id,
                case_id=case["case_id"],
                condition_id=condition_id,
                method=method,
                n_reads=n_reads,
                seed=seed,
                scope=execution_scope,
                replicate=replicate,
                is_warmup=is_warmup,
                records_processed=records_processed,
                compressed_bytes=compressed_bytes,
                report_path=report_path,
            )
        )
    if canonical_payload is None:
        raise ValueError(f"condition {condition_id} produced no measured report")
    report_identity = runtime.file_identity(canonical_path)
    report_row = {
        "study_run_id": study_run_id,
        "selection_id": selection_id,
        "case_id": case["case_id"],
        "family_id": case["family_id"],
        "configuration_accession": case["configuration_accession"],
        "fastq_accession": case["fastq_accession"],
        "selection_role": case["selection_role"],
        "condition_id": condition_id,
        "sampling_method": method,
        "requested_records_per_fastq": n_reads,
        "sampling_seed": seed,
        "records_processed": records_processed,
        "report_schema_version": canonical_payload["report_schema_version"],
        "report_path": report_identity["path"],
        "report_sha256": report_identity["sha256"],
        "result_sha256": canonical_result_sha256,
        "metric_count": runtime.count_report_metrics(canonical_payload),
    }
    metric_context = {
        key: report_row[key]
        for key in (
            "study_run_id",
            "selection_id",
            "case_id",
            "family_id",
            "configuration_accession",
            "fastq_accession",
            "selection_role",
            "condition_id",
            "sampling_method",
            "requested_records_per_fastq",
            "sampling_seed",
            "report_path",
        )
    }
    return {
        "report_row": report_row,
        "metric_rows": runtime.flatten_report_metrics(
            canonical_payload, metric_context
        ),
        "performance_rows": performance_rows,
        "repeat_mismatches": repeat_mismatches,
    }


def performance_row(
    *,
    measurement: dict[str, Any],
    study_run_id: str,
    selection_id: str,
    case_id: str,
    condition_id: str,
    method: str,
    n_reads: int,
    seed: int | None,
    scope: str,
    replicate: int,
    is_warmup: bool,
    records_processed: int,
    compressed_bytes: int,
    report_path: Path,
) -> dict[str, Any]:
    report = runtime.file_identity(report_path)
    return {
        "study_run_id": study_run_id,
        "selection_id": selection_id,
        "case_id": case_id,
        "condition_id": condition_id,
        "sampling_method": method,
        "requested_records_per_fastq": n_reads,
        "sampling_seed": seed,
        "execution_scope": scope,
        "replicate": replicate,
        "is_warmup": is_warmup,
        "runtime_schema_version": measurement["runtime_schema_version"],
        "command_json": runtime.canonical_json(measurement["argv"]),
        "exit_code": measurement["exit_code"],
        "timed_out": measurement["timed_out"],
        "wall_time_seconds": measurement["wall_time_seconds"],
        "user_cpu_seconds": measurement["user_cpu_seconds"],
        "system_cpu_seconds": measurement["system_cpu_seconds"],
        "peak_resident_memory_bytes": measurement["peak_resident_memory_bytes"],
        "records_processed": records_processed,
        "compressed_bytes_read": compressed_bytes,
        "report_size_bytes": report["size_bytes"],
        "report_path": report["path"],
        "report_sha256": report["sha256"],
        "stdout_path": measurement["stdout"]["path"],
        "stdout_sha256": measurement["stdout"]["sha256"],
        "stderr_path": measurement["stderr"]["path"],
        "stderr_sha256": measurement["stderr"]["sha256"],
    }


def validate_matrix(
    matrix: dict[str, Any], protocol: dict[str, Any], case: dict[str, Any]
) -> None:
    sampling = protocol["sampling"]
    if matrix.get("matrix_schema_version") != SCHEMA_VERSION:
        raise ValueError(f"{case['case_id']}: sample matrix schema is unsupported")
    if not str(matrix.get("matrix_id", "")).strip():
        raise ValueError(f"{case['case_id']}: sample matrix has no identity")
    if matrix.get("configuration_accession") != case["configuration_accession"]:
        raise ValueError(f"{case['case_id']}: matrix configuration differs from case")
    if matrix.get("modality") != case["modality"]:
        raise ValueError(f"{case['case_id']}: matrix modality differs from case")
    if matrix.get("sample_sizes") != sorted(sampling["records_per_fastq"]):
        raise ValueError(f"{case['case_id']}: sample sizes differ from protocol")
    if matrix.get("reservoir_seeds") != sorted(sampling["reservoir_seeds"]):
        raise ValueError(f"{case['case_id']}: reservoir seeds differ from protocol")
    if matrix.get("include_prefix") is not True:
        raise ValueError(f"{case['case_id']}: sample matrix omitted prefixes")
    if matrix.get("synchronized_mates") is not False:
        raise ValueError(
            f"{case['case_id']}: single-FASTQ matrix is marked synchronized"
        )
    sources = matrix.get("source_records")
    if not isinstance(sources, list) or len(sources) != 1:
        raise ValueError(f"{case['case_id']}: matrix must have one source")
    source = sources[0]
    if source.get("source") != case["fastq_url"]:
        raise ValueError(f"{case['case_id']}: matrix source differs from case")
    for field in ("configuration_accession", "modality", "fastq_accession", "read_id"):
        if source.get(field) != case[field]:
            raise ValueError(f"{case['case_id']}: matrix source {field} differs")
    if source.get("complete_stream_consumed") is not True:
        raise ValueError(f"{case['case_id']}: matrix did not consume complete source")
    source_sha256 = str(source.get("source_uncompressed_sha256", ""))
    if len(source_sha256) != 64 or any(
        character not in "0123456789abcdef" for character in source_sha256
    ):
        raise ValueError(f"{case['case_id']}: matrix source hash is invalid")
    if int(source.get("records_streamed", 0)) < max(sampling["records_per_fastq"]):
        raise ValueError(f"{case['case_id']}: source is shorter than largest sample")
    expected = {("prefix", size, None) for size in sampling["records_per_fastq"]} | {
        ("reservoir", size, seed)
        for size in sampling["records_per_fastq"]
        for seed in sampling["reservoir_seeds"]
    }
    conditions = matrix.get("conditions")
    if not isinstance(conditions, list):
        raise ValueError(f"{case['case_id']}: matrix conditions are not a list")
    observed = {
        (
            value.get("method"),
            value.get("requested_records_per_fastq"),
            value.get("seed"),
        )
        for value in conditions
        if isinstance(value, dict)
    }
    if observed != expected or len(conditions) != len(expected):
        raise ValueError(f"{case['case_id']}: matrix conditions differ from protocol")
    for condition in conditions:
        outputs = condition.get("outputs")
        if not isinstance(outputs, list) or len(outputs) != 1:
            raise ValueError(f"{case['case_id']}: condition must have one output")
        output = outputs[0]
        for field in ("fastq_accession", "read_id"):
            if output.get(field) != case[field]:
                raise ValueError(f"{case['case_id']}: matrix output {field} differs")
        if output.get("records_selected") != condition.get(
            "requested_records_per_fastq"
        ):
            raise ValueError(f"{case['case_id']}: bounded sample is incomplete")
        path = Path(str(output.get("output_path", "")))
        if not path.is_file() or runtime.file_sha256(path) != output.get(
            "output_sha256"
        ):
            raise ValueError(f"{case['case_id']}: sampled FASTQ identity changed")


def build_sampling_errors(
    metric_rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    numeric_rows = []
    for row in metric_rows:
        if row["metric_side"] != "observed" or row["data_kind"] != "scalar":
            continue
        try:
            value = json.loads(row["value_json"])
        except json.JSONDecodeError:
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        numeric_rows.append((row, value))

    baselines: dict[tuple[str, ...], tuple[dict[str, Any], int | float]] = {}
    issues = []
    for row, value in numeric_rows:
        if row["sampling_method"] != "complete":
            continue
        key = scalar_metric_key(row)
        if key in baselines:
            issues.append(
                {
                    "type": "duplicate_complete_metric",
                    "case_id": row["case_id"],
                    "metric_key": list(key[1:]),
                }
            )
        else:
            baselines[key] = (row, value)

    errors = []
    for row, sample_value in numeric_rows:
        if row["sampling_method"] == "complete":
            continue
        key = scalar_metric_key(row)
        baseline = baselines.get(key)
        if baseline is None:
            issues.append(
                {
                    "type": "missing_complete_metric",
                    "case_id": row["case_id"],
                    "condition_id": row["condition_id"],
                    "metric_key": list(key[1:]),
                }
            )
            continue
        _, complete_value = baseline
        signed_error = sample_value - complete_value
        errors.append(
            {
                **{
                    field: row[field]
                    for field in ERROR_FIELDS
                    if field
                    not in {
                        "complete_value",
                        "sample_value",
                        "signed_error",
                        "absolute_error",
                    }
                },
                "complete_value": complete_value,
                "sample_value": sample_value,
                "signed_error": signed_error,
                "absolute_error": abs(signed_error),
            }
        )
    return errors, issues


def scalar_metric_key(row: dict[str, Any]) -> tuple[str, ...]:
    return tuple(
        str(row[field])
        for field in (
            "case_id",
            "check",
            "files",
            "reads",
            "regions",
            "ontology_terms",
            "metric_id",
            "metric_name",
            "unit",
        )
    )


def validate_outputs(
    *,
    study_run_id: str,
    selection_id: str,
    protocol: dict[str, Any],
    cases: list[dict[str, Any]],
    output_root: Path,
    report_rows: list[dict[str, Any]],
    metric_rows: list[dict[str, Any]],
    error_rows: list[dict[str, Any]],
    error_issues: list[dict[str, Any]],
    performance_rows: list[dict[str, Any]],
    matrix_records: list[dict[str, Any]],
    repeat_mismatches: list[dict[str, Any]],
) -> dict[str, Any]:
    expected = protocol["expected_conditions_per_fastq"]
    expected_reports = len(cases) * int(expected["total"])
    bounded_per_case = int(expected["total"]) - 1
    warmups = int(protocol["performance"]["warmup_runs"])
    measured = int(protocol["performance"]["measured_runs"])
    expected_performance = len(cases) * (
        1 + 1 + bounded_per_case * (warmups + measured)
    )
    expected_repeat_reports = len(cases) * bounded_per_case * (warmups + measured - 1)
    expected_metrics = Counter(
        {row["report_path"]: int(row["metric_count"]) for row in report_rows}
    )
    observed_metrics = Counter(row["report_path"] for row in metric_rows)
    report_catalog = {Path(row["report_path"]).resolve() for row in report_rows}
    actual_reports = {
        path.resolve() for path in (output_root / "reports").glob("**/*.json")
    }
    actual_repeat_reports = list(
        (output_root / "performance_reports").glob("**/*.json")
    )
    conditions_by_case = Counter(row["case_id"] for row in report_rows)
    methods_by_case = {
        case["case_id"]: Counter(
            row["sampling_method"]
            for row in report_rows
            if row["case_id"] == case["case_id"]
        )
        for case in cases
    }
    checks = {
        "one_matrix_per_case": len(matrix_records) == len(cases),
        "all_conditions_reported": (
            len(report_rows) == expected_reports
            and all(
                conditions_by_case[case["case_id"]] == expected["total"]
                for case in cases
            )
        ),
        "condition_methods_reconcile": all(
            methods_by_case[case["case_id"]]
            == Counter(
                {
                    "prefix": int(expected["prefix"]),
                    "reservoir": int(expected["reservoir"]),
                    "complete": 1,
                }
            )
            for case in cases
        ),
        "metric_rows_reconcile": expected_metrics == observed_metrics,
        "scalar_sampling_errors_reconcile": bool(error_rows) and not error_issues,
        "report_catalog_reconciles": report_catalog == actual_reports,
        "report_hashes_reconcile": all(
            runtime.file_sha256(Path(row["report_path"])) == row["report_sha256"]
            for row in report_rows
        ),
        "repeat_results_match": not repeat_mismatches,
        "performance_rows_reconcile": len(performance_rows) == expected_performance,
        "repeat_report_catalog_reconciles": (
            len(actual_repeat_reports) == expected_repeat_reports
        ),
        "performance_commands_succeeded": all(
            row["exit_code"] == 0 and not row["timed_out"] for row in performance_rows
        ),
        "performance_measurements_present": all(
            row["wall_time_seconds"] > 0
            and row["peak_resident_memory_bytes"] > 0
            and row["report_size_bytes"] > 0
            for row in performance_rows
        ),
        "performance_provenance_reconciles": all(
            row["runtime_schema_version"] == runtime.RUNTIME_SCHEMA_VERSION
            and bool(row["command_json"])
            and runtime.file_sha256(Path(row["stdout_path"])) == row["stdout_sha256"]
            and runtime.file_sha256(Path(row["stderr_path"])) == row["stderr_sha256"]
            for row in performance_rows
        ),
        "bounded_samples_complete": all(
            row["sampling_method"] == "complete"
            or row["records_processed"] == row["requested_records_per_fastq"]
            for row in report_rows
        ),
        "report_schemas_present": all(
            str(row["report_schema_version"]).strip() for row in report_rows
        ),
        "study_ids_reconcile": all(
            row["study_run_id"] == study_run_id and row["selection_id"] == selection_id
            for row in report_rows + metric_rows + error_rows + performance_rows
        ),
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "study_run_id": study_run_id,
        "selection_id": selection_id,
        "generated_at": utc_now(),
        "valid": all(checks.values()),
        "checks": checks,
        "counts": {
            "cases": len(cases),
            "matrices": len(matrix_records),
            "reports": len(report_rows),
            "expected_reports": expected_reports,
            "metrics": len(metric_rows),
            "expected_metrics": sum(expected_metrics.values()),
            "sampling_errors": len(error_rows),
            "performance_rows": len(performance_rows),
            "expected_performance_rows": expected_performance,
            "repeat_reports": len(actual_repeat_reports),
            "expected_repeat_reports": expected_repeat_reports,
        },
        "details": {
            "repeat_mismatches": repeat_mismatches,
            "metric_count_mismatches": counter_differences(
                expected_metrics, observed_metrics
            ),
            "sampling_error_issues": error_issues,
            "missing_reports": sorted(
                str(path) for path in report_catalog - actual_reports
            ),
            "orphan_reports": sorted(
                str(path) for path in actual_reports - report_catalog
            ),
        },
    }


def source_compressed_bytes(case: dict[str, Any], matrix: dict[str, Any]) -> int:
    access = matrix["source_records"][0].get("access", {})
    content_length = access.get("content_length")
    try:
        observed = int(content_length)
    except (TypeError, ValueError):
        observed = 0
    if observed > 0:
        return observed
    source = case["fastq_url"]
    if not is_remote_source(source):
        return Path(source).stat().st_size
    return int(case["declared_compressed_bytes"])


def executable_identity(path: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [str(path), "--version"],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    version = completed.stdout.strip() or completed.stderr.strip()
    if completed.returncode != 0 or not version:
        raise ValueError(f"could not identify executable: {path}")
    return {
        "path": str(path),
        "sha256": runtime.file_sha256(path),
        "size_bytes": path.stat().st_size,
        "version": version,
    }


def script_identity(path: Path, *, version: str) -> dict[str, Any]:
    root = path.parents[1]
    status = git_output(
        root,
        "status",
        "--porcelain",
        "--untracked-files=normal",
        "--",
        str(path.relative_to(root)),
    )
    return {
        "version": version,
        "path": str(path),
        "sha256": runtime.file_sha256(path),
        "size_bytes": path.stat().st_size,
        "git_commit": git_output(root, "rev-parse", "HEAD"),
        "git_dirty": bool(status),
        "python": sys.version.split()[0],
    }


def functional_executable_identity(value: dict[str, Any]) -> dict[str, Any]:
    return {"sha256": value["sha256"], "version": value["version"]}


def functional_script_identity(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "version": value["version"],
        "sha256": value["sha256"],
        "python": value["python"],
    }


def git_output(root: Path, *args: str) -> str:
    if not (root / ".git").exists():
        return ""
    completed = subprocess.run(
        ["git", "-C", str(root), *args],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def required_string(value: dict[str, Any], field: str) -> str:
    result = str(value.get(field, "")).strip()
    if not result:
        raise ValueError(f"sampling case is missing {field}")
    return result


def is_remote_source(value: str) -> bool:
    return urllib.parse.urlparse(value).scheme in {"http", "https"}


def counter_differences(
    expected: Counter[str], observed: Counter[str]
) -> list[dict[str, Any]]:
    return [
        {
            "key": key,
            "expected": expected[key],
            "observed": observed[key],
        }
        for key in sorted(set(expected) | set(observed))
        if expected[key] != observed[key]
    ]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
