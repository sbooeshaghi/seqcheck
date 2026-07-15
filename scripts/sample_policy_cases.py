#!/usr/bin/env python3
"""Create frozen-policy FASTQ samples for the locked evaluation split."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import build_perturbation_cases as inventory_builder
    import materialize_perturbations as materializer
    import paper_runtime as runtime
    import sample_fastq as sampler
except ModuleNotFoundError:
    from scripts import build_perturbation_cases as inventory_builder
    from scripts import materialize_perturbations as materializer
    from scripts import paper_runtime as runtime
    from scripts import sample_fastq as sampler


SCHEMA_VERSION = "0.1.0"
RUNNER_VERSION = "0.1.0"
SAMPLE_FIELDS = (
    "bundle_id",
    "selection_id",
    "sampling_policy_id",
    "case_id",
    "family_id",
    "configuration_accession",
    "modality",
    "selection_role",
    "fastq_accession",
    "read_id",
    "sample_id",
    "sampling_method",
    "requested_records_per_fastq",
    "sampling_seed",
    "records_streamed",
    "records_selected",
    "complete_stream_consumed",
    "source_uncompressed_sha256",
    "output_path",
    "output_sha256",
    "output_size_bytes",
    "sample_manifest_path",
    "sample_manifest_sha256",
)
PERFORMANCE_FIELDS = (
    "bundle_id",
    "case_id",
    "configuration_accession",
    "fastq_accession",
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
        description="Sample every locked evaluation FASTQ with a frozen policy."
    )
    parser.add_argument("--case-selection-manifest", required=True, type=Path)
    parser.add_argument("--sampling-study", required=True, type=Path)
    parser.add_argument("--sampling-policy", required=True, type=Path)
    parser.add_argument("--sampler-script", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=3600)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = sample_policy_cases(
            case_selection_manifest_path=args.case_selection_manifest.resolve(),
            sampling_study_path=args.sampling_study.resolve(),
            sampling_policy_path=args.sampling_policy.resolve(),
            sampler_script=args.sampler_script.resolve(),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
        )
    except (OSError, ValueError) as error:
        print(f"sample_policy_cases: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def sample_policy_cases(
    *,
    case_selection_manifest_path: Path,
    sampling_study_path: Path,
    sampling_policy_path: Path,
    sampler_script: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    if sampler_script != Path(sampler.__file__).resolve():
        raise ValueError("sampler script must be the repository sample_fastq.py")
    selection, cases = inventory_builder.load_case_selection(
        case_selection_manifest_path
    )
    if selection.get("cohort_split") != "evaluation":
        raise ValueError("policy samples require an evaluation case selection")
    study = load_sampling_study(sampling_study_path)
    if selection["sampling_protocol_sha256"] != study["sampling_protocol_sha256"]:
        raise ValueError("evaluation selection uses a different sampling protocol")
    policy = materializer.load_sampling_policy(
        sampling_policy_path,
        study_run_id=study["study_run_id"],
        require_frozen=True,
    )
    tools = {
        "sampler": runtime.script_identity(
            sampler_script, version=sampler.SAMPLER_VERSION
        ),
        "runner": runtime.script_identity(
            Path(__file__).resolve(), version=RUNNER_VERSION
        ),
        "runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
    }
    if runtime.functional_script_identity(tools["sampler"]) != study["sampler_core"]:
        raise ValueError("evaluation sampler differs from calibration sampler core")

    output_root.mkdir(parents=True)
    sample_rows = []
    measurements = []
    for case in cases:
        result = sample_case(
            case=case,
            policy=policy,
            sampler_script=sampler_script,
            output_root=output_root,
            timeout_seconds=timeout_seconds,
        )
        sample_rows.append(result["sample"])
        measurements.append(result["measurement"])

    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "cohort_split": "evaluation",
        "selection_id": selection["selection_id"],
        "study_run_id": study["study_run_id"],
        "sampling_policy_id": policy["policy_id"],
        "case_selection_sha256": runtime.file_sha256(case_selection_manifest_path),
        "sampling_study_sha256": runtime.file_sha256(sampling_study_path),
        "sampling_policy_sha256": runtime.file_sha256(sampling_policy_path),
        "sampling_rule": policy["default"],
        "sampler": runtime.functional_script_identity(tools["sampler"]),
        "runner": runtime.functional_script_identity(tools["runner"]),
        "runtime": runtime.functional_script_identity(tools["runtime"]),
        "samples": [stable_sample_row(value) for value in sample_rows],
    }
    bundle_id = runtime.sha256_json(stable_identity)[:16]
    samples = [{"bundle_id": bundle_id, **value} for value in sample_rows]
    performance = [
        performance_row(bundle_id, value["case"], value) for value in measurements
    ]
    validation = validate_outputs(
        bundle_id=bundle_id,
        selection=selection,
        cases=cases,
        samples=samples,
        performance=performance,
        policy=policy,
    )
    if not validation["valid"]:
        raise ValueError(
            "policy sample bundle failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )

    samples_path = output_root / "tables" / "policy_samples.csv"
    performance_path = output_root / "tables" / "performance.csv"
    validation_path = output_root / "validation" / "policy_samples.json"
    manifest_path = output_root / "manifests" / "sample_bundle.json"
    runtime.write_csv(samples_path, samples, list(SAMPLE_FIELDS))
    runtime.write_csv(performance_path, performance, list(PERFORMANCE_FIELDS))
    runtime.write_json(validation_path, validation)
    manifest = {
        **stable_identity,
        "bundle_id": bundle_id,
        "created_at": utc_now(),
        "valid": True,
        "inputs": {
            "case_selection_manifest": runtime.file_identity(
                case_selection_manifest_path
            ),
            "sampling_study": runtime.file_identity(sampling_study_path),
            "sampling_policy": runtime.file_identity(sampling_policy_path),
        },
        "tools": tools,
        "sample_records": samples,
        "counts": validation["counts"],
        "outputs": {
            "samples": runtime.file_identity(samples_path),
            "performance": runtime.file_identity(performance_path),
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_sampling_study(path: Path) -> dict[str, Any]:
    study = runtime.load_json(path)
    if study.get("valid") is not True:
        raise ValueError("sampling calibration study is not valid")
    stable = {field: study.get(field) for field in materializer.STUDY_IDENTITY_FIELDS}
    if any(stable[field] is None for field in materializer.STUDY_IDENTITY_FIELDS):
        raise ValueError("sampling calibration study identity is incomplete")
    if runtime.sha256_json(stable)[:16] != study.get("study_run_id"):
        raise ValueError("sampling calibration study id is not content-addressed")
    protocol_path = verified_path(
        study.get("inputs", {}).get("sampling_protocol", {}),
        "sampling protocol",
    )
    if runtime.file_sha256(protocol_path) != study["sampling_protocol_sha256"]:
        raise ValueError("sampling calibration protocol identities differ")
    return study


def sample_case(
    *,
    case: dict[str, Any],
    policy: dict[str, Any],
    sampler_script: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    rule = policy["default"]
    case_id = case["case_id"]
    argv = [
        sys.executable,
        str(sampler_script),
        "--input",
        case["fastq_url"],
        "--output-root",
        str(output_root / "samples" / case_id),
        "--method",
        rule["sampling_method"],
        "--n-reads",
        str(rule["records_per_fastq"]),
        "--configuration-accession",
        case["configuration_accession"],
        "--modality",
        case["modality"],
        "--fastq-accession",
        case["fastq_accession"],
        "--read-id",
        case["read_id"],
    ]
    if rule["sampling_method"] == "reservoir":
        argv.extend(["--seed", str(rule["sampling_seed"])])
    measurement = runtime.run_measured_command(
        argv,
        stdout_path=output_root / "logs" / case_id / "stdout.txt",
        stderr_path=output_root / "logs" / case_id / "stderr.txt",
        timeout_seconds=timeout_seconds,
    )
    runtime.require_success(measurement)
    manifest_path = Path(
        Path(measurement["stdout"]["path"]).read_text(encoding="utf-8").strip()
    ).resolve()
    sample_manifest = runtime.load_json(manifest_path)
    sample = validate_sample_manifest(sample_manifest, manifest_path, case, rule)
    return {
        "sample": {
            "selection_id": case["selection_id"],
            "sampling_policy_id": policy["policy_id"],
            "case_id": case_id,
            "family_id": case["family_id"],
            "configuration_accession": case["configuration_accession"],
            "modality": case["modality"],
            "selection_role": case["selection_role"],
            "fastq_accession": case["fastq_accession"],
            "read_id": case["read_id"],
            **sample,
        },
        "measurement": {
            "case": case,
            "measurement": measurement,
            "records_processed": sample["records_streamed"],
            "complete_stream_consumed": sample["complete_stream_consumed"],
        },
    }


def validate_sample_manifest(
    value: dict[str, Any],
    path: Path,
    case: dict[str, Any],
    rule: dict[str, Any],
) -> dict[str, Any]:
    expected = {
        "configuration_accession": case["configuration_accession"],
        "modality": case["modality"],
        "sampling_method": rule["sampling_method"],
        "requested_records_per_fastq": rule["records_per_fastq"],
        "sampling_seed": rule["sampling_seed"],
    }
    if value.get("sample_schema_version") != sampler.SAMPLE_SCHEMA_VERSION:
        raise ValueError(f"{case['case_id']}: sample schema is unsupported")
    if any(
        value.get(key) != expected_value for key, expected_value in expected.items()
    ):
        raise ValueError(f"{case['case_id']}: sample rule or context differs")
    inputs = value.get("inputs")
    if not isinstance(inputs, list) or len(inputs) != 1:
        raise ValueError(f"{case['case_id']}: sample must contain exactly one FASTQ")
    row = inputs[0]
    if row.get("fastq_accession") != case["fastq_accession"]:
        raise ValueError(f"{case['case_id']}: sampled FASTQ accession differs")
    if row.get("read_id") != case["read_id"]:
        raise ValueError(f"{case['case_id']}: sampled read id differs")
    if row.get("records_selected") != rule["records_per_fastq"]:
        raise ValueError(f"{case['case_id']}: source has fewer records than policy")
    output_path = Path(str(row.get("output_path", ""))).resolve()
    if runtime.file_sha256(output_path) != row.get("output_sha256"):
        raise ValueError(f"{case['case_id']}: sampled FASTQ hash changed")
    identity = {
        key: value.get(key)
        for key in (
            "sample_schema_version",
            "configuration_accession",
            "modality",
            "sampling_method",
            "requested_records_per_fastq",
            "sampling_seed",
            "synchronized_mates",
        )
    }
    identity["inputs"] = [sampler.identity_input_row(row)]
    identity["sampler"] = {
        "version": value.get("sampler", {}).get("version"),
        "script_sha256": value.get("sampler", {}).get("script_sha256"),
    }
    if sampler.sha256_json(identity)[:16] != value.get("sample_id"):
        raise ValueError(f"{case['case_id']}: sample id is not content-addressed")
    return {
        "sample_id": value["sample_id"],
        "sampling_method": value["sampling_method"],
        "requested_records_per_fastq": value["requested_records_per_fastq"],
        "sampling_seed": value["sampling_seed"],
        "records_streamed": row["records_streamed"],
        "records_selected": row["records_selected"],
        "complete_stream_consumed": row["complete_stream_consumed"],
        "source_uncompressed_sha256": row["source_uncompressed_sha256"],
        "output_path": str(output_path),
        "output_sha256": row["output_sha256"],
        "output_size_bytes": output_path.stat().st_size,
        "sample_manifest_path": str(path),
        "sample_manifest_sha256": runtime.file_sha256(path),
    }


def stable_sample_row(value: dict[str, Any]) -> dict[str, Any]:
    return {
        key: item
        for key, item in value.items()
        if key
        not in {
            "bundle_id",
            "output_path",
            "sample_manifest_path",
            "sample_manifest_sha256",
        }
    }


def performance_row(
    bundle_id: str,
    case: dict[str, Any],
    context: dict[str, Any],
) -> dict[str, Any]:
    measurement = context["measurement"]
    return {
        "bundle_id": bundle_id,
        "case_id": case["case_id"],
        "configuration_accession": case["configuration_accession"],
        "fastq_accession": case["fastq_accession"],
        "runtime_schema_version": measurement["runtime_schema_version"],
        "command_json": runtime.canonical_json(measurement["argv"]),
        "exit_code": measurement["exit_code"],
        "timed_out": measurement["timed_out"],
        "wall_time_seconds": measurement["wall_time_seconds"],
        "user_cpu_seconds": measurement["user_cpu_seconds"],
        "system_cpu_seconds": measurement["system_cpu_seconds"],
        "peak_resident_memory_bytes": measurement["peak_resident_memory_bytes"],
        "records_processed": context["records_processed"],
        "compressed_bytes_read": (
            case["declared_compressed_bytes"]
            if context["complete_stream_consumed"]
            else ""
        ),
        "stdout_path": measurement["stdout"]["path"],
        "stdout_sha256": measurement["stdout"]["sha256"],
        "stderr_path": measurement["stderr"]["path"],
        "stderr_sha256": measurement["stderr"]["sha256"],
    }


def validate_outputs(
    *,
    bundle_id: str,
    selection: dict[str, Any],
    cases: list[dict[str, Any]],
    samples: list[dict[str, Any]],
    performance: list[dict[str, Any]],
    policy: dict[str, Any],
) -> dict[str, Any]:
    errors = []
    case_ids = Counter(value["case_id"] for value in cases)
    if Counter(value["case_id"] for value in samples) != case_ids:
        errors.append("sample rows do not reconcile with selected cases")
    if Counter(value["case_id"] for value in performance) != case_ids:
        errors.append("performance rows do not reconcile with selected cases")
    if any(value["bundle_id"] != bundle_id for value in samples + performance):
        errors.append("sample bundle identifiers differ")
    if any(value["selection_id"] != selection["selection_id"] for value in samples):
        errors.append("sample selection identifiers differ")
    if any(value["sampling_policy_id"] != policy["policy_id"] for value in samples):
        errors.append("sample policy identifiers differ")
    if any(
        runtime.file_sha256(Path(value["output_path"])) != value["output_sha256"]
        for value in samples
    ):
        errors.append("sampled FASTQ hashes changed")
    return {
        "schema_version": SCHEMA_VERSION,
        "bundle_id": bundle_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "counts": {
            "configurations": len(
                {value["configuration_accession"] for value in samples}
            ),
            "fastq_cases": len(cases),
            "samples": len(samples),
            "records_selected": sum(value["records_selected"] for value in samples),
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
