#!/usr/bin/env python3
"""Run the integrated reproducibility acceptance gate for the paper study."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import shutil
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any


SCHEMA_VERSION = "0.1.0"
CURRENT_SEQSPEC_VERSION = "0.5.0"
USER_AGENT = "seqcheck-phase0/0.1.0"
CASE_ID_PATTERN = "abcdefghijklmnopqrstuvwxyz0123456789._-"
RUN_FIELDS = (
    "study_run_id",
    "case_id",
    "spec_variant",
    "fastq_transport",
    "raw_seqspec_version",
    "seqspec_version",
    "requested_reads_per_fastq",
    "fastq_count",
    "report_schema_version",
    "report_path",
    "report_sha256",
    "cache_key",
    "metric_count",
)
METRIC_FIELDS = (
    "study_run_id",
    "cache_key",
    "case_id",
    "spec_variant",
    "fastq_transport",
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
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run raw-versus-normalized seqcheck cases, deterministic sampling, "
            "cache probes, and output reconciliation for the paper Phase 0 gate."
        )
    )
    parser.add_argument("--case-manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--seqcheck-bin", required=True, type=Path)
    parser.add_argument("--seqspec-bin", required=True, type=Path)
    parser.add_argument(
        "--sampler-script",
        type=Path,
        default=Path(__file__).resolve().with_name("sample_fastq.py"),
    )
    parser.add_argument("--fastqc-command", default="fastqc")
    parser.add_argument("--fastqc-limits", type=Path)
    parser.add_argument("--command-timeout-seconds", type=int, default=180)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        return run(args)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"run_phase0: {error}", file=sys.stderr)
        return 1


def run(args: argparse.Namespace) -> int:
    if args.command_timeout_seconds <= 0:
        raise ValueError("command timeout must be positive")

    protocol_path = args.case_manifest.resolve()
    protocol = load_protocol(protocol_path)
    output_root = args.output_root.resolve()
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    output_root.mkdir(parents=True)

    seqcheck = executable_identity(args.seqcheck_bin, args.command_timeout_seconds)
    seqspec = executable_identity(args.seqspec_bin, args.command_timeout_seconds)
    sampler = script_identity(args.sampler_script)
    fastqc = external_tool_identity(
        args.fastqc_command,
        args.fastqc_limits,
        args.command_timeout_seconds,
    )
    runner = script_identity(Path(__file__).resolve())
    tools = {
        "seqcheck": seqcheck,
        "seqspec": seqspec,
        "sampler": sampler,
        "fastqc": fastqc,
        "phase0_runner": runner,
    }

    copied_protocol = output_root / "inputs" / "phase0_cases.json"
    copied_protocol.parent.mkdir(parents=True)
    shutil.copyfile(protocol_path, copied_protocol)

    prepared_cases = [
        prepare_case(
            case,
            protocol_path.parent,
            output_root,
            Path(seqspec["resolved_path"]),
            args.command_timeout_seconds,
        )
        for case in protocol["cases"]
    ]
    sampling_sources = sampling_source_identities(
        protocol["sampling_probe"], protocol_path.parent
    )
    run_identity = {
        "schema_version": SCHEMA_VERSION,
        "protocol_sha256": file_sha256(protocol_path),
        "tools": {
            name: functional_tool_identity(value) for name, value in tools.items()
        },
        "cases": [case["identity"] for case in prepared_cases],
        "sampling_probe": {
            **stable_sampling_probe(protocol["sampling_probe"]),
            "inputs": sampling_sources,
        },
    }
    study_run_id = sha256_json(run_identity)[:16]

    run_rows: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    parity_errors = []
    reports_by_case: dict[str, list[dict[str, Any]]] = {}
    for prepared in prepared_cases:
        case_reports = []
        for variant, spec_path, version in (
            ("raw", prepared["raw_spec_path"], prepared["raw_seqspec_version"]),
            (
                "normalized",
                prepared["normalized_spec_path"],
                prepared["normalized_seqspec_version"],
            ),
        ):
            report = run_case_report(
                prepared=prepared,
                variant=variant,
                spec_path=Path(spec_path),
                seqspec_version=version,
                study_run_id=study_run_id,
                seqcheck=seqcheck,
                tools=tools,
                output_root=output_root,
                timeout_seconds=args.command_timeout_seconds,
            )
            case_reports.append(report)
            run_rows.append(report["run_row"])
            metric_rows.extend(report["metric_rows"])
        reports_by_case[prepared["case_id"]] = case_reports
        raw_results = canonical_report_results(case_reports[0]["payload"], prepared)
        normalized_results = canonical_report_results(
            case_reports[1]["payload"], prepared
        )
        if raw_results != normalized_results:
            parity_errors.append(prepared["case_id"])

    sample_probe = run_sampling_probe(
        protocol["sampling_probe"],
        protocol_path.parent,
        output_root,
        Path(sampler["path"]),
        args.command_timeout_seconds,
    )
    cache_probe = build_cache_probe(run_rows, tools)

    tables_dir = output_root / "tables"
    tables_dir.mkdir()
    runs_path = tables_dir / "runs.csv"
    metrics_path = tables_dir / "metrics.csv"
    write_csv(runs_path, run_rows, RUN_FIELDS)
    write_csv(metrics_path, metric_rows, METRIC_FIELDS)

    validation = validate_outputs(
        study_run_id=study_run_id,
        output_root=output_root,
        prepared_cases=prepared_cases,
        reports_by_case=reports_by_case,
        run_rows=run_rows,
        metric_rows=metric_rows,
        parity_errors=parity_errors,
        sample_probe=sample_probe,
        cache_probe=cache_probe,
        tools=tools,
    )
    validation_path = output_root / "validation" / "phase0.json"
    write_json(validation_path, validation)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "study_run_id": study_run_id,
        "generated_at": utc_now(),
        "valid": validation["valid"],
        "invocation": sys.argv,
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "tooling": tools,
        "inputs": {
            "case_manifest": file_identity(protocol_path),
            "copied_case_manifest": file_identity(copied_protocol),
            "cases": [case["manifest_record"] for case in prepared_cases],
        },
        "sampling_probe": sample_probe,
        "cache_probe": cache_probe,
        "outputs": {
            "runs": file_identity(runs_path),
            "metrics": file_identity(metrics_path),
            "validation": file_identity(validation_path),
            "reports": [
                file_identity(Path(report["run_row"]["report_path"]))
                for reports in reports_by_case.values()
                for report in reports
            ],
        },
    }
    manifest_path = output_root / "manifests" / "study.json"
    write_json(manifest_path, manifest)
    print(manifest_path)
    return 0 if validation["valid"] else 1


def load_protocol(path: Path) -> dict[str, Any]:
    protocol = load_json(path)
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Phase 0 case manifest must use schema 0.1.0")
    cases = protocol.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("Phase 0 case manifest needs at least one case")
    seen = set()
    for case in cases:
        if not isinstance(case, dict):
            raise ValueError("Phase 0 cases must be objects")
        case_id = str(case.get("case_id", "")).strip()
        if (
            not case_id
            or any(character not in CASE_ID_PATTERN for character in case_id)
            or case_id in seen
        ):
            raise ValueError(f"invalid or duplicate Phase 0 case id: {case_id!r}")
        seen.add(case_id)
        if not str(case.get("spec", "")).strip():
            raise ValueError(f"case {case_id} is missing its spec")
        if not str(case.get("modality", "")).strip():
            raise ValueError(f"case {case_id} is missing its modality")
        fastqs = case.get("fastqs")
        if not isinstance(fastqs, list) or not fastqs:
            raise ValueError(f"case {case_id} needs at least one FASTQ")
        n_reads = case.get("n_reads", 0)
        if not isinstance(n_reads, int) or n_reads < 0:
            raise ValueError(f"case {case_id} n_reads must be a nonnegative integer")
        validate_resources(case_id, case.get("resources", []))

    probe = protocol.get("sampling_probe")
    if not isinstance(probe, dict):
        raise ValueError("Phase 0 case manifest needs a sampling_probe object")
    inputs = probe.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError("sampling_probe needs at least one input")
    if not isinstance(probe.get("n_reads"), int) or probe["n_reads"] <= 0:
        raise ValueError("sampling_probe n_reads must be a positive integer")
    if not isinstance(probe.get("seed"), int):
        raise ValueError("sampling_probe seed must be an integer")
    return protocol


def validate_resources(case_id: str, resources: Any) -> None:
    if not isinstance(resources, list):
        raise ValueError(f"case {case_id} resources must be a list")
    destinations = set()
    for resource in resources:
        if (
            not isinstance(resource, dict)
            or not str(resource.get("source", "")).strip()
        ):
            raise ValueError(f"case {case_id} has an invalid resource source")
        value = str(resource.get("path", "")).strip()
        destination = PurePosixPath(value)
        if not value or destination.is_absolute() or ".." in destination.parts:
            raise ValueError(f"case {case_id} has an unsafe resource path: {value!r}")
        if value in destinations:
            raise ValueError(f"case {case_id} has duplicate resource path: {value}")
        destinations.add(value)


def prepare_case(
    case: dict[str, Any],
    protocol_root: Path,
    output_root: Path,
    seqspec_bin: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    case_id = str(case["case_id"])
    raw_root = output_root / "specs" / "raw" / case_id
    normalized_root = output_root / "specs" / "normalized" / case_id
    raw_root.mkdir(parents=True)
    normalized_root.mkdir(parents=True)
    raw_spec = raw_root / "spec.yaml"
    normalized_spec = normalized_root / "spec.yaml"
    source_spec = resolve_source(str(case["spec"]), protocol_root)
    raw_source = materialize_source(source_spec, raw_spec)

    resource_records = []
    for resource in case.get("resources", []):
        source = resolve_source(str(resource["source"]), protocol_root)
        relative = PurePosixPath(str(resource["path"]))
        raw_destination = raw_root.joinpath(*relative.parts)
        normalized_destination = normalized_root.joinpath(*relative.parts)
        record = materialize_source(source, raw_destination)
        normalized_destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(raw_destination, normalized_destination)
        resource_records.append(
            {
                **record,
                "relative_path": str(relative),
                "normalized_copy": file_identity(normalized_destination),
            }
        )

    raw_version = seqspec_file_version(seqspec_bin, raw_spec, timeout_seconds)
    run_command(
        [str(seqspec_bin), "upgrade", str(raw_spec), "-o", str(normalized_spec)],
        timeout_seconds,
    )
    normalized_version = seqspec_file_version(
        seqspec_bin, normalized_spec, timeout_seconds
    )
    run_command(
        [str(seqspec_bin), "check", "--skip", "external", str(normalized_spec)],
        timeout_seconds,
    )

    fastq_sources = [
        resolve_source(str(source), protocol_root) for source in case["fastqs"]
    ]
    fastq_identities = [source_identity(source) for source in fastq_sources]
    transport = classify_transport(fastq_sources)
    identity = {
        "case_id": case_id,
        "modality": str(case["modality"]),
        "n_reads": int(case.get("n_reads", 0)),
        "transport": transport,
        "raw_spec_sha256": file_sha256(raw_spec),
        "normalized_spec_sha256": file_sha256(normalized_spec),
        "raw_seqspec_version": raw_version,
        "normalized_seqspec_version": normalized_version,
        "fastqs": [stable_source_identity(value) for value in fastq_identities],
        "resources": [
            {
                "relative_path": value["relative_path"],
                "sha256": value["materialized_sha256"],
            }
            for value in resource_records
        ],
    }
    return {
        "case_id": case_id,
        "modality": str(case["modality"]),
        "n_reads": int(case.get("n_reads", 0)),
        "fastq_sources": fastq_sources,
        "fastq_identities": fastq_identities,
        "transport": transport,
        "raw_spec_path": str(raw_spec),
        "normalized_spec_path": str(normalized_spec),
        "raw_seqspec_version": raw_version,
        "normalized_seqspec_version": normalized_version,
        "identity": identity,
        "manifest_record": {
            **identity,
            "source_spec": raw_source,
            "raw_spec": file_identity(raw_spec),
            "normalized_spec": file_identity(normalized_spec),
            "fastq_sources": fastq_identities,
            "resources": resource_records,
        },
    }


def run_case_report(
    *,
    prepared: dict[str, Any],
    variant: str,
    spec_path: Path,
    seqspec_version: str,
    study_run_id: str,
    seqcheck: dict[str, Any],
    tools: dict[str, dict[str, Any]],
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    report_path = output_root / "reports" / prepared["case_id"] / f"{variant}.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    argv = [
        seqcheck["resolved_path"],
        "check",
        "--spec",
        str(spec_path),
        "--modality",
        prepared["modality"],
        "--n-reads",
        str(prepared["n_reads"]),
        "--format",
        "json",
        "--output",
        str(report_path),
        *prepared["fastq_sources"],
    ]
    run_command(argv, timeout_seconds)
    payload = load_json(report_path)
    report_schema = str(payload.get("report_schema_version", "")).strip()
    if not report_schema:
        raise ValueError(f"seqcheck report has no schema version: {report_path}")
    cache_payload = {
        "phase0_schema_version": SCHEMA_VERSION,
        "report_schema_version": report_schema,
        "seqcheck": functional_tool_identity(tools["seqcheck"]),
        "seqspec": functional_tool_identity(tools["seqspec"]),
        "spec_sha256": file_sha256(spec_path),
        "fastqs": [
            stable_source_identity(value) for value in prepared["fastq_identities"]
        ],
        "sampling": {
            "method": "prefix",
            "n_reads_per_fastq": prepared["n_reads"],
            "seed": None,
        },
    }
    cache_key = sha256_json(cache_payload)
    metric_rows = flatten_metrics(
        payload=payload,
        study_run_id=study_run_id,
        cache_key=cache_key,
        case_id=prepared["case_id"],
        variant=variant,
        transport=prepared["transport"],
        report_path=report_path,
    )
    return {
        "payload": payload,
        "cache_payload": cache_payload,
        "metric_rows": metric_rows,
        "run_row": {
            "study_run_id": study_run_id,
            "case_id": prepared["case_id"],
            "spec_variant": variant,
            "fastq_transport": prepared["transport"],
            "raw_seqspec_version": prepared["raw_seqspec_version"],
            "seqspec_version": seqspec_version,
            "requested_reads_per_fastq": prepared["n_reads"],
            "fastq_count": len(prepared["fastq_sources"]),
            "report_schema_version": report_schema,
            "report_path": str(report_path),
            "report_sha256": file_sha256(report_path),
            "cache_key": cache_key,
            "metric_count": len(metric_rows),
        },
    }


def flatten_metrics(
    *,
    payload: dict[str, Any],
    study_run_id: str,
    cache_key: str,
    case_id: str,
    variant: str,
    transport: str,
    report_path: Path,
) -> list[dict[str, Any]]:
    rows = []
    for result_index, result in enumerate(payload.get("results", [])):
        if not isinstance(result, dict):
            continue
        ontology_terms = ";".join(sorted(extract_ontology_terms(result)))
        for side in ("expected", "observed"):
            values = result.get(side, [])
            if not isinstance(values, list):
                continue
            for metric in values:
                if not isinstance(metric, dict):
                    continue
                data = metric.get("data", {})
                if not isinstance(data, dict):
                    data = {}
                rows.append(
                    {
                        "study_run_id": study_run_id,
                        "cache_key": cache_key,
                        "case_id": case_id,
                        "spec_variant": variant,
                        "fastq_transport": transport,
                        "report_path": str(report_path),
                        "result_index": result_index,
                        "check": str(result.get("check", "")),
                        "files": join_values(result.get("files", [])),
                        "reads": join_values(result.get("reads", [])),
                        "regions": join_values(result.get("regions", [])),
                        "ontology_terms": ontology_terms,
                        "metric_side": side,
                        "metric_id": str(metric.get("id", "")),
                        "metric_name": str(metric.get("name", "")),
                        "metric_description": str(metric.get("description", "")),
                        "data_kind": str(data.get("kind", "")),
                        "unit": str(data.get("unit", "")),
                        "value_json": canonical_json(data.get("value")),
                    }
                )
    return rows


def run_sampling_probe(
    probe: dict[str, Any],
    protocol_root: Path,
    output_root: Path,
    sampler_script: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    inputs = [resolve_source(str(value), protocol_root) for value in probe["inputs"]]
    if any(is_remote_source(value) for value in inputs):
        raise ValueError("Phase 0 sampling_probe inputs must be local files")
    manifests = {}
    for label, seed in (
        ("repeat_a", int(probe["seed"])),
        ("repeat_b", int(probe["seed"])),
        ("changed_seed", int(probe["seed"]) + 1),
    ):
        sample_root = output_root / "samples" / label
        argv = [
            sys.executable,
            str(sampler_script),
            "--output-root",
            str(sample_root),
            "--method",
            "reservoir",
            "--n-reads",
            str(probe["n_reads"]),
            "--seed",
            str(seed),
        ]
        for value in inputs:
            argv.extend(["--input", value])
        for field, option in (
            ("fastq_accessions", "--fastq-accession"),
            ("read_ids", "--read-id"),
        ):
            for value in probe.get(field, []):
                argv.extend([option, str(value)])
        for field, option in (
            ("configuration_accession", "--configuration-accession"),
            ("modality", "--modality"),
        ):
            value = str(probe.get(field, "")).strip()
            if value:
                argv.extend([option, value])
        if bool(probe.get("synchronize_mates", False)):
            argv.append("--synchronize-mates")
        completed = run_command(argv, timeout_seconds)
        manifest_path = Path(completed.stdout.strip())
        manifest = load_json(manifest_path)
        manifests[label] = {
            "sample_id": str(manifest.get("sample_id", "")),
            "sampling_seed": manifest.get("sampling_seed"),
            "manifest": file_identity(manifest_path),
            "output_sha256": [
                str(value.get("output_sha256", ""))
                for value in manifest.get("inputs", [])
                if isinstance(value, dict)
            ],
            "records_selected": [
                int(value.get("records_selected", 0))
                for value in manifest.get("inputs", [])
                if isinstance(value, dict)
            ],
        }
    return {
        "method": "reservoir",
        "requested_records_per_fastq": int(probe["n_reads"]),
        "seed": int(probe["seed"]),
        "inputs": [source_identity(value) for value in inputs],
        "manifests": manifests,
    }


def build_cache_probe(
    run_rows: list[dict[str, Any]], tools: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    if not run_rows:
        raise ValueError("cannot build a cache probe without a report run")
    baseline = {
        "phase0_schema_version": SCHEMA_VERSION,
        "report_schema_version": run_rows[0]["report_schema_version"],
        "tool": functional_tool_identity(tools["seqcheck"]),
        "input_key": run_rows[0]["cache_key"],
        "sampling": {
            "method": "prefix",
            "n_reads_per_fastq": run_rows[0]["requested_reads_per_fastq"],
            "seed": None,
        },
    }
    changed_tool = deepcopy(baseline)
    changed_tool["tool"]["executable_sha256"] += ":changed"
    changed_sampling = deepcopy(baseline)
    changed_sampling["sampling"] = {
        "method": "reservoir",
        "n_reads_per_fastq": int(run_rows[0]["requested_reads_per_fastq"]) + 1,
        "seed": 1,
    }
    return {
        "baseline": sha256_json(baseline),
        "changed_tool": sha256_json(changed_tool),
        "changed_sampling": sha256_json(changed_sampling),
    }


def validate_outputs(
    *,
    study_run_id: str,
    output_root: Path,
    prepared_cases: list[dict[str, Any]],
    reports_by_case: dict[str, list[dict[str, Any]]],
    run_rows: list[dict[str, Any]],
    metric_rows: list[dict[str, Any]],
    parity_errors: list[str],
    sample_probe: dict[str, Any],
    cache_probe: dict[str, Any],
    tools: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    raw_versions = {value["raw_seqspec_version"] for value in prepared_cases}
    transports = {value["transport"] for value in prepared_cases}
    expected_metrics_by_report = {
        report["run_row"]["report_path"]: count_report_metrics(report["payload"])
        for reports in reports_by_case.values()
        for report in reports
    }
    actual_metrics_by_report = Counter(row["report_path"] for row in metric_rows)
    run_metrics_by_report = {
        row["report_path"]: int(row["metric_count"]) for row in run_rows
    }
    expected_metrics = sum(expected_metrics_by_report.values())
    catalog_reports = {Path(row["report_path"]).resolve() for row in run_rows}
    actual_reports = {
        path.resolve() for path in (output_root / "reports").glob("**/*.json")
    }
    repeat_a = sample_probe["manifests"]["repeat_a"]
    repeat_b = sample_probe["manifests"]["repeat_b"]
    changed_seed = sample_probe["manifests"]["changed_seed"]
    remote_sources = [
        source
        for case in prepared_cases
        for source in case["fastq_identities"]
        if source["kind"] == "remote"
    ]
    checks = {
        "required_seqspec_versions_present": {"0.3.0", "0.5.0"}.issubset(raw_versions),
        "required_fastq_transports_present": {"local", "remote"}.issubset(transports),
        "all_specs_normalized_to_0_5_0": all(
            value["normalized_seqspec_version"] == CURRENT_SEQSPEC_VERSION
            for value in prepared_cases
        ),
        "raw_and_normalized_results_match": not parity_errors,
        "two_reports_per_case": all(
            len(reports_by_case.get(value["case_id"], [])) == 2
            for value in prepared_cases
        ),
        "report_schemas_present": all(row["report_schema_version"] for row in run_rows),
        "metric_rows_reconcile": (
            expected_metrics_by_report == actual_metrics_by_report
        ),
        "run_metric_counts_reconcile": (
            expected_metrics_by_report == run_metrics_by_report
        ),
        "report_catalog_reconciles": catalog_reports == actual_reports,
        "study_run_ids_match": all(
            row["study_run_id"] == study_run_id for row in run_rows + metric_rows
        ),
        "remote_content_identifiers_present": all(
            source.get("etag") or source.get("last_modified")
            for source in remote_sources
        ),
        "identical_sampler_seeds_match": (
            repeat_a["sample_id"] == repeat_b["sample_id"]
            and repeat_a["output_sha256"] == repeat_b["output_sha256"]
        ),
        "changed_sampler_seed_invalidates_identity": (
            repeat_a["sample_id"] != changed_seed["sample_id"]
        ),
        "changed_tool_invalidates_cache": (
            cache_probe["baseline"] != cache_probe["changed_tool"]
        ),
        "changed_sampling_invalidates_cache": (
            cache_probe["baseline"] != cache_probe["changed_sampling"]
        ),
        "sampler_record_counts_match": all(
            count == sample_probe["requested_records_per_fastq"]
            for manifest in sample_probe["manifests"].values()
            for count in manifest["records_selected"]
        ),
        "seqcheck_identity_complete": tool_identity_complete(tools["seqcheck"]),
        "seqspec_identity_complete": tool_identity_complete(tools["seqspec"]),
        "fastqc_identity_complete": tool_identity_complete(tools["fastqc"]),
        "runtime_source_identities_complete": all(
            tools[name].get("runtime_source_sha256") for name in ("seqcheck", "seqspec")
        ),
        "script_identities_complete": all(
            script_identity_complete(tools[name])
            for name in ("sampler", "phase0_runner")
        ),
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "study_run_id": study_run_id,
        "generated_at": utc_now(),
        "valid": all(checks.values()),
        "checks": checks,
        "counts": {
            "cases": len(prepared_cases),
            "reports": len(run_rows),
            "expected_reports": len(prepared_cases) * 2,
            "metrics": len(metric_rows),
            "expected_metrics": expected_metrics,
            "local_cases": sum(
                value["transport"] == "local" for value in prepared_cases
            ),
            "remote_cases": sum(
                value["transport"] == "remote" for value in prepared_cases
            ),
        },
        "details": {
            "raw_seqspec_versions": sorted(raw_versions),
            "fastq_transports": sorted(transports),
            "parity_errors": parity_errors,
            "metric_count_mismatches": counter_differences(
                expected_metrics_by_report, actual_metrics_by_report
            ),
            "run_metric_count_mismatches": counter_differences(
                expected_metrics_by_report, run_metrics_by_report
            ),
            "missing_reports": sorted(
                str(value) for value in catalog_reports - actual_reports
            ),
            "orphan_reports": sorted(
                str(value) for value in actual_reports - catalog_reports
            ),
        },
    }


def count_report_metrics(payload: dict[str, Any]) -> int:
    return sum(
        len(result.get(side, []))
        for result in payload.get("results", [])
        if isinstance(result, dict)
        for side in ("expected", "observed")
        if isinstance(result.get(side, []), list)
    )


def canonical_report_results(payload: dict[str, Any], prepared: dict[str, Any]) -> Any:
    replacements = {
        prepared["raw_spec_path"]: f"spec://{prepared['case_id']}",
        prepared["normalized_spec_path"]: f"spec://{prepared['case_id']}",
    }
    for resource in prepared["manifest_record"]["resources"]:
        label = f"resource://{resource['relative_path']}"
        replacements[resource["materialized_path"]] = label
        replacements[resource["normalized_copy"]["path"]] = label
    return replace_exact_strings(payload.get("results", []), replacements)


def replace_exact_strings(value: Any, replacements: dict[str, str]) -> Any:
    if isinstance(value, dict):
        return {
            key: replace_exact_strings(item, replacements)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [replace_exact_strings(item, replacements) for item in value]
    if isinstance(value, str):
        return replacements.get(value, value)
    return value


def seqspec_file_version(binary: Path, spec: Path, timeout_seconds: int) -> str:
    output = run_command([str(binary), "version", str(spec)], timeout_seconds).stdout
    for line in output.splitlines():
        if line.startswith("seqspec file version:"):
            return line.split(":", 1)[1].strip()
    raise ValueError(f"seqspec version output omitted the file version for {spec}")


def classify_transport(sources: list[str]) -> str:
    kinds = {"remote" if is_remote_source(value) else "local" for value in sources}
    if len(kinds) != 1:
        raise ValueError("Phase 0 cases must not mix local and remote FASTQs")
    return next(iter(kinds))


def resolve_source(value: str, root: Path) -> str:
    if is_remote_source(value):
        return value
    path = Path(value)
    if not path.is_absolute():
        path = root / path
    return str(path.resolve())


def is_remote_source(value: str) -> bool:
    return urllib.parse.urlparse(value).scheme in {"http", "https"}


def materialize_source(source: str, destination: Path) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if is_remote_source(source):
        request = urllib.request.Request(source, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request) as response:
            payload = response.read()
            metadata = remote_response_identity(source, response)
        destination.write_bytes(payload)
    else:
        source_path = Path(source)
        shutil.copyfile(source_path, destination)
        metadata = {
            "source": str(source_path),
            "kind": "local",
            "status": "available",
            "size_bytes": source_path.stat().st_size,
            "sha256": file_sha256(source_path),
        }
    return {
        **metadata,
        "materialized_path": str(destination.resolve()),
        "materialized_sha256": file_sha256(destination),
    }


def source_identity(source: str) -> dict[str, Any]:
    if not is_remote_source(source):
        path = Path(source)
        return {
            "source": str(path),
            "kind": "local",
            "status": "available",
            "size_bytes": path.stat().st_size,
            "sha256": file_sha256(path),
        }
    request = urllib.request.Request(
        source, headers={"User-Agent": USER_AGENT}, method="HEAD"
    )
    try:
        response = urllib.request.urlopen(request)
    except urllib.error.HTTPError as error:
        if error.code not in {405, 501}:
            raise
        response = urllib.request.urlopen(
            urllib.request.Request(
                source,
                headers={"User-Agent": USER_AGENT, "Range": "bytes=0-0"},
            )
        )
    with response:
        return remote_response_identity(source, response)


def remote_response_identity(source: str, response: Any) -> dict[str, Any]:
    return {
        "source": source,
        "kind": "remote",
        "status": int(getattr(response, "status", 200)),
        "etag": response.headers.get("ETag", ""),
        "last_modified": response.headers.get("Last-Modified", ""),
        "content_length": response.headers.get("Content-Length", ""),
        "retrieved_at": utc_now(),
    }


def stable_source_identity(value: dict[str, Any]) -> dict[str, Any]:
    excluded = {"retrieved_at", "status"}
    if value.get("kind") == "local":
        excluded.add("source")
    return {key: item for key, item in value.items() if key not in excluded}


def sampling_source_identities(
    probe: dict[str, Any], protocol_root: Path
) -> list[dict[str, Any]]:
    return [
        stable_source_identity(
            source_identity(resolve_source(str(value), protocol_root))
        )
        for value in probe["inputs"]
    ]


def stable_sampling_probe(probe: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in probe.items() if key not in {"inputs"}}


def executable_identity(path: Path, timeout_seconds: int) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise ValueError(f"executable does not exist: {resolved}")
    runtime_scope = ("Cargo.toml", "Cargo.lock", "src")
    source_identity = git_identity(resolved.parent, runtime_scope)
    source_root = str(source_identity.get("git_root", ""))
    identity = {
        "command": str(path),
        "resolved_path": str(resolved),
        "available": True,
        "version": run_command(
            [str(resolved), "--version"], timeout_seconds
        ).stdout.strip(),
        "executable_sha256": file_sha256(resolved),
        "runtime_source_sha256": (
            path_tree_sha256(Path(source_root), runtime_scope)
            if source_root
            else file_sha256(resolved)
        ),
    }
    return {**identity, **source_identity}


def external_tool_identity(
    command: str, limits: Path | None, timeout_seconds: int
) -> dict[str, Any]:
    resolved_value = shutil.which(command)
    if resolved_value is None:
        raise ValueError(f"external tool is unavailable: {command}")
    resolved = Path(resolved_value).resolve()
    configuration = {
        "mode": "custom" if limits else "packaged_defaults",
        "path": str(limits.resolve()) if limits else "",
        "sha256": file_sha256(limits.resolve()) if limits else "",
    }
    return {
        "command": command,
        "resolved_path": str(resolved),
        "available": True,
        "version": run_command(
            [str(resolved), "--version"], timeout_seconds
        ).stdout.strip(),
        "executable_sha256": file_sha256(resolved),
        "configuration": configuration,
    }


def script_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise ValueError(f"script does not exist: {resolved}")
    return {
        "path": str(resolved),
        "sha256": file_sha256(resolved),
        "python": sys.version.split()[0],
        **git_identity(resolved.parent, (str(resolved),)),
    }


def git_identity(path: Path, scope: tuple[str, ...] | None = None) -> dict[str, Any]:
    root = git_output(path, "rev-parse", "--show-toplevel")
    if not root:
        return {
            "git_root": "",
            "git_commit": "",
            "git_dirty": False,
            "git_scope": [],
        }
    status_args = ["status", "--porcelain", "--untracked-files=normal"]
    if scope:
        status_args.extend(["--", *scope])
    return {
        "git_root": root,
        "git_commit": git_output(Path(root), "rev-parse", "HEAD"),
        "git_dirty": bool(git_output(Path(root), *status_args)),
        "git_scope": list(scope or (".",)),
    }


def git_output(path: Path, *args: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(path), *args],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else ""


def functional_tool_identity(value: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "version",
        "executable_sha256",
        "runtime_source_sha256",
        "sha256",
        "python",
    )
    identity = {key: value[key] for key in keys if key in value}
    configuration = value.get("configuration")
    if isinstance(configuration, dict):
        identity["configuration"] = {
            "mode": configuration.get("mode", ""),
            "sha256": configuration.get("sha256", ""),
        }
    return identity


def tool_identity_complete(value: dict[str, Any]) -> bool:
    return bool(
        value.get("available", True)
        and value.get("version")
        and value.get("executable_sha256")
    )


def script_identity_complete(value: dict[str, Any]) -> bool:
    return bool(value.get("sha256") and value.get("python") and value.get("git_commit"))


def counter_differences(
    expected: dict[str, int], observed: dict[str, int]
) -> list[dict[str, Any]]:
    return [
        {
            "report_path": path,
            "expected": expected.get(path, 0),
            "observed": observed.get(path, 0),
        }
        for path in sorted(set(expected).union(observed))
        if expected.get(path, 0) != observed.get(path, 0)
    ]


def run_command(
    argv: list[str], timeout_seconds: int
) -> subprocess.CompletedProcess[str]:
    try:
        completed = subprocess.run(
            argv,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        raise ValueError(
            f"command exceeded {timeout_seconds} seconds: {' '.join(argv)}"
        ) from error
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip()
        raise ValueError(f"{' '.join(argv)} failed: {message}")
    return completed


def extract_ontology_terms(value: Any) -> set[str]:
    terms = set()
    if isinstance(value, dict):
        for child in value.values():
            terms.update(extract_ontology_terms(child))
    elif isinstance(value, list):
        for child in value:
            terms.update(extract_ontology_terms(child))
    elif isinstance(value, str) and value.startswith("RGN:"):
        terms.add(value)
    return terms


def join_values(value: Any) -> str:
    if not isinstance(value, list):
        return ""
    return ";".join(str(item) for item in value)


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def write_csv(path: Path, rows: list[dict[str, Any]], fields: tuple[str, ...]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "sha256": file_sha256(resolved),
        "size_bytes": resolved.stat().st_size,
    }


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def path_tree_sha256(root: Path, paths: tuple[str, ...]) -> str:
    files = []
    for value in paths:
        path = root / value
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            files.extend(
                candidate for candidate in path.rglob("*") if candidate.is_file()
            )
    digest = hashlib.sha256()
    for path in sorted(files):
        digest.update(str(path.relative_to(root)).encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_sha256(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
