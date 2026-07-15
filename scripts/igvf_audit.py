#!/usr/bin/env python3
"""Audit IGVF seqspec-backed FASTQ datasets with seqcheck."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


USER_AGENT = "seqcheck-igvf-audit/0.1"
DEFAULT_API_ROOT = "https://api.data.igvf.org/"
DEFAULT_PORTAL_ROOT = DEFAULT_API_ROOT
AUDIT_SCHEMA_VERSION = "0.4.1"
PORTAL_MANIFEST_SCHEMA_VERSION = "0.1.0"
SAMPLING_POLICY_SCHEMA_VERSION = "0.1.0"
AUDIT_PROTOCOL_SCHEMA_VERSION = "0.1.0"
CURRENT_SEQSPEC_VERSION = "0.5.0"
SAMPLING_METHOD = "prefix"
ONTOLOGY_TERM_PATTERN = re.compile(r"RGN:[A-Za-z0-9_]+:[A-Za-z0-9_]+")
ACCESS_CLASSES = {"public", "controlled", "mixed", "unknown"}
FAILURE_CATEGORIES = {
    "portal",
    "transport",
    "authentication",
    "specification",
    "tool",
    "eligibility",
    "unclassified",
}


@dataclass(frozen=True)
class SequenceFileRecord:
    accession: str
    href: str
    controlled_access: bool
    read_names: list[str]


@dataclass(frozen=True)
class ConfigurationRecord:
    accession: str
    href: str
    lab: str
    submitted_by: str
    award_component: str
    file_set_accession: str
    assay_term: str
    preferred_assay_titles: list[str]
    aliases: list[str]
    seqspec_of: list[str]
    status: str
    upload_status: str


@dataclass(frozen=True)
class RunRecord:
    study_run_id: str
    cache_key: str
    sampling_method: str
    sampling_seed: int | None
    configuration_accession: str
    modality: str
    modality_count: int
    lab: str
    submitted_by: str
    award_component: str
    file_set_accession: str
    assay_term: str
    preferred_assay_titles: str
    aliases: str
    raw_seqspec_version: str
    normalized_seqspec_version: str
    requested_reads: int
    requested_reads_total: int
    sampled_record_count: int
    expected_fastq_count: int
    supplied_fastq_count: int
    controlled_access: bool
    access_class: str
    report_path: str
    run_status: str
    pass_count: int
    warning_count: int
    error_count: int
    interpretation_count: int
    max_attempt_count: int = 1


@dataclass(frozen=True)
class DiagnosticRecord:
    study_run_id: str
    cache_key: str
    configuration_accession: str
    modality: str
    lab: str
    submitted_by: str
    award_component: str
    file_set_accession: str
    assay_term: str
    preferred_assay_titles: str
    aliases: str
    raw_seqspec_version: str
    normalized_seqspec_version: str
    report_path: str
    check: str
    assessment_type: str
    assessment_code: str
    assessment_description: str
    files: str
    reads: str
    regions: str
    ontology_terms: str
    sequence_types: str
    region_annotations_json: str
    access_class: str


@dataclass(frozen=True)
class MetricRecord:
    study_run_id: str
    cache_key: str
    configuration_accession: str
    modality: str
    lab: str
    submitted_by: str
    award_component: str
    file_set_accession: str
    assay_term: str
    preferred_assay_titles: str
    aliases: str
    raw_seqspec_version: str
    normalized_seqspec_version: str
    report_path: str
    result_index: int
    check: str
    files: str
    reads: str
    regions: str
    ontology_terms: str
    sequence_types: str
    region_annotations_json: str
    access_class: str
    metric_side: str
    metric_id: str
    metric_name: str
    metric_description: str
    data_kind: str
    unit: str
    value_json: str


@dataclass(frozen=True)
class FailureRecord:
    study_run_id: str
    configuration_accession: str
    modality: str
    lab: str
    submitted_by: str
    award_component: str
    file_set_accession: str
    assay_term: str
    preferred_assay_titles: str
    aliases: str
    raw_seqspec_version: str
    normalized_seqspec_version: str
    stage: str
    reason: str
    message: str
    access_class: str = "unknown"
    failure_category: str = "tool"
    attempt_count: int = 1


@dataclass(frozen=True)
class ToolCommand:
    argv: list[str]
    env: dict[str, str]


@dataclass(frozen=True)
class ToolIdentity:
    command: list[str]
    version: str
    git_commit: str
    git_dirty: bool
    runtime_source_sha256: str
    executable_sha256: str


@dataclass(frozen=True)
class AuditContext:
    run_id: str
    sampling_method: str
    sampling_seed: int | None
    seqcheck: ToolIdentity
    seqspec: ToolIdentity
    command_timeout_seconds: int | None = None
    max_attempts: int = 1
    retry_delays_seconds: tuple[float, ...] = ()
    retryable_failure_categories: tuple[str, ...] = ("transport",)
    audit_protocol_id: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit IGVF seqspec-backed FASTQ datasets with seqcheck."
    )
    parser.add_argument("--output-root", required=True, help="Directory for reports.")
    parser.add_argument(
        "--n-reads",
        type=int,
        default=10000,
        help="Reads to sample from each FASTQ (0 means all reads).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Concurrent audit workers.",
    )
    parser.add_argument(
        "--status",
        default="released",
        help="IGVF object status filter for configuration files.",
    )
    parser.add_argument(
        "--upload-status",
        default="validated",
        help="IGVF upload_status filter for configuration files.",
    )
    parser.add_argument("--limit", type=int, help="Limit to the first N seqspecs.")
    parser.add_argument(
        "--configuration-accession",
        action="append",
        default=[],
        help="Specific configuration-file accession to audit. Repeatable.",
    )
    parser.add_argument(
        "--portal-manifest",
        type=Path,
        help="Frozen content-addressed portal manifest to consume instead of live search.",
    )
    parser.add_argument(
        "--sampling-policy",
        type=Path,
        help="Frozen Experiment 1 sampling policy; currently must select prefix sampling.",
    )
    parser.add_argument(
        "--audit-protocol",
        type=Path,
        help="Versioned retry, repeatability, and completion rules for the paper audit.",
    )
    parser.add_argument(
        "--public-only",
        action="store_true",
        help="Skip controlled-access FASTQs.",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Retained for compatibility. Remote-mode audits do not stage local inputs.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun even if a report already exists and parses.",
    )
    parser.add_argument(
        "--auth-profile",
        default="igvf",
        help="Auth profile name to use with seqspec and seqcheck for remote resources.",
    )
    parser.add_argument(
        "--api-root",
        default=DEFAULT_API_ROOT,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--portal-root",
        default=DEFAULT_PORTAL_ROOT,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--fastqc-command",
        default="fastqc",
        help="FastQC executable to record for the comparison study.",
    )
    parser.add_argument(
        "--fastqc-limits",
        help="Optional FastQC limits file to identify in the study manifest.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started_at = utc_now()
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "reports").mkdir(exist_ok=True)
    mpl_dir = output_root / ".mplconfig"
    mpl_dir.mkdir(exist_ok=True)

    workspace_root = Path(__file__).resolve().parents[2]
    seqcheck_root = Path(__file__).resolve().parents[1]
    seqspec_cmd = discover_seqspec_command(workspace_root)
    seqcheck_cmd = discover_seqcheck_command(seqcheck_root)
    seqspec_identity = build_tool_identity(seqspec_cmd, workspace_root / "seqspec")
    seqcheck_identity = build_tool_identity(seqcheck_cmd, seqcheck_root)

    seqspec_auth_profile = resolve_ready_auth_profile(
        args.auth_profile,
        config_env="SEQSPEC_AUTH_CONFIG",
        config_subdir="seqspec",
    )
    seqcheck_auth_profile = resolve_ready_auth_profile(
        args.auth_profile,
        config_env="SEQCHECK_AUTH_CONFIG",
        config_subdir="seqcheck",
    )

    try:
        configurations, sequence_files, portal_manifest, sampling_policy, audit_protocol = (
            load_audit_inputs(args)
        )
    except (OSError, ValueError) as err:
        print(f"igvf_audit: {err}", file=sys.stderr)
        return 2

    if args.configuration_accession:
        try:
            configurations = select_requested_configurations(
                configurations, args.configuration_accession
            )
        except ValueError as err:
            print(f"igvf_audit: {err}", file=sys.stderr)
            return 2

    configurations.sort(key=lambda record: record.accession)
    if args.limit is not None:
        configurations = configurations[: args.limit]
    if not configurations:
        print("igvf_audit: no configurations matched the active filters", file=sys.stderr)
        return 2

    audit_context = write_study_manifest(
        output_root,
        args,
        configurations,
        sequence_files,
        seqcheck_identity,
        seqspec_identity,
        seqcheck_auth_profile,
        seqspec_auth_profile,
        started_at,
        portal_manifest=portal_manifest,
        sampling_policy=sampling_policy,
        audit_protocol=audit_protocol,
    )

    print(
        f"auditing {len(configurations)} configuration files against {len(sequence_files)} FASTQ records",
        file=sys.stderr,
    )
    total_configurations = len(configurations)

    runs: list[RunRecord] = []
    diagnostics: list[DiagnosticRecord] = []
    metrics: list[MetricRecord] = []
    failures: list[FailureRecord] = []

    with ThreadPoolExecutor(max_workers=max(args.workers, 1)) as pool:
        future_map = {
            pool.submit(
                process_configuration,
                record,
                sequence_files,
                args,
                output_root,
                seqspec_cmd,
                seqcheck_cmd,
                seqspec_auth_profile,
                seqcheck_auth_profile,
                args.portal_root,
                mpl_dir,
                audit_context,
            ): record
            for record in configurations
        }

        completed_configurations = 0
        for future in as_completed(future_map):
            record = future_map[future]
            accession = record.accession
            completed_configurations += 1
            remaining_configurations = total_configurations - completed_configurations
            try:
                run_rows, diagnostic_rows, metric_rows, failure_rows = future.result()
            except Exception as err:  # pragma: no cover - hard to force deterministically
                print(
                    (
                        f"[{completed_configurations}/{total_configurations}] "
                        f"{accession} unhandled error: {err} "
                        f"({remaining_configurations} remaining)"
                    ),
                    file=sys.stderr,
                )
                failures.append(
                    build_failure(
                        audit_context.run_id,
                        record,
                        modality="",
                        raw_seqspec_version="",
                        normalized_seqspec_version="",
                        stage="worker",
                        reason="unhandled_exception",
                        message=str(err),
                    )
                )
                continue

            runs.extend(run_rows)
            diagnostics.extend(diagnostic_rows)
            metrics.extend(metric_rows)
            failures.extend(failure_rows)
            print(
                (
                    f"[{completed_configurations}/{total_configurations}] "
                    f"{accession} completed: {len(run_rows)} runs, "
                    f"{len(failure_rows)} failures "
                    f"({remaining_configurations} remaining)"
                ),
                file=sys.stderr,
            )

    runs.sort(key=lambda row: (row.configuration_accession, row.modality))
    diagnostics.sort(
        key=lambda row: (
            row.configuration_accession,
            row.modality,
            row.check,
            row.assessment_type,
            row.assessment_code,
        )
    )
    metrics.sort(
        key=lambda row: (
            row.configuration_accession,
            row.modality,
            row.result_index,
            row.metric_side,
            row.metric_id,
        )
    )
    failures.sort(
        key=lambda row: (row.configuration_accession, row.modality, row.stage, row.reason)
    )

    write_rows(output_root / "runs.csv", output_root / "runs.jsonl", runs)
    write_rows(
        output_root / "diagnostics.csv",
        output_root / "diagnostics.jsonl",
        diagnostics,
    )
    write_rows(
        output_root / "failures.csv",
        output_root / "failures.jsonl",
        failures,
    )
    write_rows(
        output_root / "metrics.csv",
        output_root / "metrics.jsonl",
        metrics,
    )
    write_lab_summary(output_root / "lab_summary.csv", runs, failures)
    reconciliation = reconcile_outputs(
        output_root,
        audit_context,
        configurations,
        runs,
        diagnostics,
        metrics,
        failures,
    )
    write_json(output_root / "validation" / "reconciliation.json", reconciliation)
    return 0 if reconciliation["valid"] else 1


def load_audit_inputs(
    args: argparse.Namespace,
) -> tuple[
    list[ConfigurationRecord],
    dict[str, SequenceFileRecord],
    dict[str, Any] | None,
    dict[str, Any] | None,
    dict[str, Any] | None,
]:
    portal_manifest = None
    if args.portal_manifest:
        configurations, sequence_files, portal_manifest = load_portal_manifest(
            args.portal_manifest.resolve()
        )
        query = portal_manifest["query"]
        args.api_root = query["api_root"]
        args.portal_root = query["portal_root"]
        args.status = query["status"]
        args.upload_status = query["upload_status"]
    else:
        configurations = fetch_configuration_records(
            args.api_root,
            status=args.status,
            upload_status=args.upload_status,
        )
        sequence_files = fetch_sequence_file_records(
            args.api_root,
            status=args.status,
            upload_status=args.upload_status,
        )

    sampling_policy = None
    if args.sampling_policy:
        sampling_policy = load_frozen_sampling_policy(args.sampling_policy.resolve())
        args.n_reads = sampling_policy["default"]["records_per_fastq"]
    audit_protocol = None
    audit_protocol_path = getattr(args, "audit_protocol", None)
    if audit_protocol_path:
        audit_protocol = load_audit_protocol(audit_protocol_path.resolve())
    return (
        configurations,
        sequence_files,
        portal_manifest,
        sampling_policy,
        audit_protocol,
    )


def process_configuration(
    record: ConfigurationRecord,
    sequence_files: dict[str, SequenceFileRecord],
    args: argparse.Namespace,
    output_root: Path,
    seqspec_cmd: ToolCommand,
    seqcheck_cmd: ToolCommand,
    seqspec_auth_profile: str | None,
    seqcheck_auth_profile: str | None,
    portal_root: str,
    mpl_dir: Path,
    audit_context: AuditContext,
) -> tuple[
    list[RunRecord],
    list[DiagnosticRecord],
    list[MetricRecord],
    list[FailureRecord],
]:
    runs: list[RunRecord] = []
    diagnostics: list[DiagnosticRecord] = []
    metrics: list[MetricRecord] = []
    failures: list[FailureRecord] = []

    linked_fastqs = [
        sequence_files[accession]
        for accession in extract_accessions(record.seqspec_of)
        if accession in sequence_files
    ]
    configuration_access_class = access_class_for_records(linked_fastqs)

    if not linked_fastqs:
        failures.append(
            build_failure(
                audit_context.run_id,
                record,
                modality="",
                raw_seqspec_version="",
                normalized_seqspec_version="",
                stage="sequence_files",
                reason="no_fastq_sequence_files",
                message="Configuration file has no linked FASTQ sequence files.",
                access_class=configuration_access_class,
            )
        )
        return runs, diagnostics, metrics, failures

    spec_url = absolute_url(portal_root, record.href)

    try:
        raw_seqspec_version, configuration_max_attempt_count = get_seqspec_file_version(
            seqspec_cmd,
            spec_url,
            seqspec_auth_profile,
            mpl_dir,
            audit_context,
        )
    except subprocess.CalledProcessError as err:
        message = stderr_message(err)
        failures.append(
            build_failure(
                audit_context.run_id,
                record,
                modality="",
                raw_seqspec_version="",
                normalized_seqspec_version="",
                stage="seqspec_version",
                reason=classify_seqspec_failure_reason(
                    default_reason="seqspec_version_error",
                    message=message,
                ),
                message=message,
                access_class=configuration_access_class,
                attempt_count=command_attempt_count(err),
            )
        )
        return runs, diagnostics, metrics, failures

    normalized_seqspec_version = normalize_seqspec_version(raw_seqspec_version)

    try:
        modalities, modality_list_attempt_count = list_modalities(
            seqspec_cmd,
            spec_url,
            seqspec_auth_profile,
            mpl_dir,
            audit_context,
        )
    except subprocess.CalledProcessError as err:
        message = stderr_message(err)
        failures.append(
            build_failure(
                audit_context.run_id,
                record,
                modality="",
                raw_seqspec_version=raw_seqspec_version,
                normalized_seqspec_version=normalized_seqspec_version,
                stage="enumerate_modalities",
                reason=classify_seqspec_failure_reason(
                    default_reason="seqspec_info_error",
                    message=message,
                ),
                message=message,
                access_class=configuration_access_class,
                attempt_count=command_attempt_count(err),
            )
        )
        return runs, diagnostics, metrics, failures

    configuration_max_attempt_count = max(
        configuration_max_attempt_count, modality_list_attempt_count
    )

    if not modalities:
        failures.append(
            build_failure(
                audit_context.run_id,
                record,
                modality="",
                raw_seqspec_version=raw_seqspec_version,
                normalized_seqspec_version=normalized_seqspec_version,
                stage="enumerate_modalities",
                reason="no_modalities",
                message="No modalities found in seqspec.",
                access_class=configuration_access_class,
            )
        )
        return runs, diagnostics, metrics, failures

    modality_count = len(modalities)

    for modality in modalities:
        modality_max_attempt_count = configuration_max_attempt_count
        report_path = output_root / "reports" / record.accession / f"{slugify(modality)}.json"
        report_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            region_annotations, region_attempt_count = list_modality_region_annotations(
                seqspec_cmd,
                spec_url,
                modality,
                seqspec_auth_profile,
                mpl_dir,
                audit_context,
            )
        except (subprocess.CalledProcessError, ValueError) as err:
            message = stderr_message(err) if isinstance(
                err, subprocess.CalledProcessError
            ) else str(err)
            failures.append(
                build_failure(
                    audit_context.run_id,
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="enumerate_regions",
                    reason=classify_seqspec_failure_reason(
                        default_reason="seqspec_region_error",
                        message=message,
                    ),
                    message=message,
                    access_class=configuration_access_class,
                    attempt_count=(
                        command_attempt_count(err)
                        if isinstance(err, subprocess.CalledProcessError)
                        else 1
                    ),
                )
            )
            continue
        modality_max_attempt_count = max(
            modality_max_attempt_count, region_attempt_count
        )

        try:
            expected_files, file_attempt_count = list_modality_files(
                seqspec_cmd,
                spec_url,
                modality,
                seqspec_auth_profile,
                mpl_dir,
                audit_context,
            )
        except subprocess.CalledProcessError as err:
            message = stderr_message(err)
            failures.append(
                build_failure(
                    audit_context.run_id,
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="enumerate_files",
                    reason=classify_seqspec_failure_reason(
                        default_reason="seqspec_file_error",
                        message=message,
                    ),
                    message=message,
                    access_class=configuration_access_class,
                    attempt_count=command_attempt_count(err),
                )
            )
            continue
        modality_max_attempt_count = max(
            modality_max_attempt_count, file_attempt_count
        )

        fastq_expectations = [item for item in expected_files if is_fastq_expectation(item)]
        if not fastq_expectations:
            failures.append(
                build_failure(
                    audit_context.run_id,
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="enumerate_files",
                    reason="no_fastq_files_for_modality",
                    message=f"No FASTQ files found for modality '{modality}'.",
                    access_class=configuration_access_class,
                )
            )
            continue

        modality_sequence_records: list[SequenceFileRecord] = []
        fastq_urls: list[str] = []
        controlled_needed = False
        resolved_for_access: list[SequenceFileRecord] = []
        modality_access_class = configuration_access_class

        skip_reason = None
        skip_message = None

        for item in fastq_expectations:
            sequence_record = resolve_expected_sequence_file(item, sequence_files)
            if sequence_record is None:
                skip_reason = "missing_fastq_metadata"
                skip_message = (
                    f"Could not resolve FASTQ metadata for expected file "
                    f"'{item.get('file_id') or item.get('filename')}'."
                )
                break

            resolved_for_access.append(sequence_record)
            modality_access_class = access_class_for_records(resolved_for_access)

            if sequence_record.controlled_access:
                controlled_needed = True
                if args.public_only:
                    skip_reason = "controlled_fastq_skipped"
                    skip_message = (
                        "Encountered controlled-access FASTQ while --public-only was enabled."
                    )
                    break
                if seqcheck_auth_profile is None:
                    skip_reason = "missing_credentials"
                    skip_message = (
                        "Controlled-access FASTQ requires a ready seqcheck auth profile."
                    )
                    break

            modality_sequence_records.append(sequence_record)
            fastq_urls.append(absolute_url(portal_root, sequence_record.href))

        if skip_reason is not None:
            failures.append(
                build_failure(
                    audit_context.run_id,
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="resolve_fastqs",
                    reason=skip_reason,
                    message=skip_message or "",
                    access_class=modality_access_class,
                )
            )
            continue

        cache_key = build_report_cache_key(
            audit_context,
            record,
            modality,
            raw_seqspec_version,
            normalized_seqspec_version,
            spec_url,
            expected_files,
            modality_sequence_records,
            fastq_urls,
            region_annotations,
            args.n_reads,
            controlled_needed,
        )

        if report_path.exists() and not args.force:
            try:
                report = json.loads(report_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                report = None
            if report is not None and report_cache_matches(report, cache_key):
                report = annotate_report_summary(
                    report,
                    configuration_accession=record.accession,
                    modality=modality,
                    modality_count=modality_count,
                    expected_fastq_count=len(modality_sequence_records),
                    supplied_fastq_count=len(modality_sequence_records),
                    controlled_access=controlled_needed,
                    audit_context=audit_context,
                    cache_key=cache_key,
                    access_class=modality_access_class,
                )
                report_path.write_text(
                    json.dumps(report, indent=2, sort_keys=False) + "\n",
                    encoding="utf-8",
                )
                run_row, diagnostic_rows, metric_rows = flatten_report(
                    report,
                    record,
                    modality,
                    raw_seqspec_version,
                    normalized_seqspec_version,
                    report_path,
                    modality_count=modality_count,
                    expected_fastq_count=len(modality_sequence_records),
                    supplied_fastq_count=len(modality_sequence_records),
                    controlled_access=controlled_needed,
                    run_status="cached",
                    region_annotations=region_annotations,
                    access_class=modality_access_class,
                )
                runs.append(run_row)
                diagnostics.extend(diagnostic_rows)
                metrics.extend(metric_rows)
                continue

        try:
            seqcheck_attempt_count = run_seqcheck(
                seqcheck_cmd,
                spec_url,
                modality,
                args.n_reads,
                fastq_urls,
                report_path,
                seqcheck_auth_profile,
                audit_context,
            )
            report = json.loads(report_path.read_text(encoding="utf-8"))
            report = annotate_report_summary(
                report,
                configuration_accession=record.accession,
                modality=modality,
                modality_count=modality_count,
                expected_fastq_count=len(modality_sequence_records),
                supplied_fastq_count=len(modality_sequence_records),
                controlled_access=controlled_needed,
                audit_context=audit_context,
                cache_key=cache_key,
                access_class=modality_access_class,
                max_attempt_count=max(
                    modality_max_attempt_count, seqcheck_attempt_count
                ),
            )
            report_path.write_text(
                json.dumps(report, indent=2, sort_keys=False) + "\n",
                encoding="utf-8",
            )
        except subprocess.CalledProcessError as err:
            failures.append(
                build_failure(
                    audit_context.run_id,
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="seqcheck",
                    reason="seqcheck_error",
                    message=stderr_message(err),
                    access_class=modality_access_class,
                    attempt_count=command_attempt_count(err),
                )
            )
            continue
        except Exception as err:
            failures.append(
                build_failure(
                    audit_context.run_id,
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="seqcheck",
                    reason="report_parse_error",
                    message=str(err),
                    access_class=modality_access_class,
                )
            )
            continue

        run_row, diagnostic_rows, metric_rows = flatten_report(
            report,
            record,
            modality,
            raw_seqspec_version,
            normalized_seqspec_version,
            report_path,
            modality_count=modality_count,
            expected_fastq_count=len(modality_sequence_records),
            supplied_fastq_count=len(modality_sequence_records),
            controlled_access=controlled_needed,
            run_status="completed",
            region_annotations=region_annotations,
            access_class=modality_access_class,
        )
        runs.append(run_row)
        diagnostics.extend(diagnostic_rows)
        metrics.extend(metric_rows)

    return runs, diagnostics, metrics, failures


def absolute_url(api_root: str, href_or_url: str) -> str:
    return urllib.parse.urljoin(api_root, href_or_url)


def select_requested_configurations(
    configurations: list[ConfigurationRecord], requested: list[str]
) -> list[ConfigurationRecord]:
    requested_accessions = {value.strip() for value in requested if value.strip()}
    available = {record.accession for record in configurations}
    missing = sorted(requested_accessions - available)
    if missing:
        raise ValueError(
            "requested configuration accessions were not found under the active "
            f"portal filters: {', '.join(missing)}"
        )
    return [
        record for record in configurations if record.accession in requested_accessions
    ]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def functional_script_provenance(value: dict[str, Any]) -> dict[str, Any]:
    required = {"version", "sha256", "python"}
    if not isinstance(value, dict) or not required.issubset(value):
        raise ValueError("portal manifest tool identity is incomplete")
    return {key: value[key] for key in sorted(required)}


def portal_manifest_stable(value: dict[str, Any]) -> dict[str, Any]:
    tools = value.get("tools", {})
    if not isinstance(tools, dict) or set(tools) != {"audit_runtime", "freezer"}:
        raise ValueError("portal manifest must identify its freezer and audit runtime")
    return {
        "schema_version": value.get("schema_version"),
        "query": value.get("query"),
        "configurations": value.get("configurations"),
        "sequence_files": value.get("sequence_files"),
        "missing_linked_sequence_file_accessions": value.get(
            "missing_linked_sequence_file_accessions"
        ),
        "tools": {
            name: functional_script_provenance(identity)
            for name, identity in sorted(tools.items())
        },
    }


def sampling_policy_stable(value: dict[str, Any]) -> dict[str, Any]:
    return {
        key: item
        for key, item in value.items()
        if key not in {"policy_id", "created_at"}
    }


def load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as err:
        raise ValueError(f"{label} is not valid JSON: {path}") from err
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return value


def load_portal_manifest(
    path: Path,
) -> tuple[list[ConfigurationRecord], dict[str, SequenceFileRecord], dict[str, Any]]:
    value = load_json_object(path, "portal manifest")
    if value.get("schema_version") != PORTAL_MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            "unsupported portal manifest schema: "
            f"{value.get('schema_version', '')}"
        )
    expected_id = sha256_json(portal_manifest_stable(value))[:16]
    if value.get("portal_manifest_id") != expected_id:
        raise ValueError("portal manifest ID does not reconcile with its contents")

    query = value.get("query")
    query_fields = {"api_root", "portal_root", "status", "upload_status"}
    if not isinstance(query, dict) or set(query) != query_fields:
        raise ValueError("portal manifest query fields are incomplete or unexpected")
    if any(not isinstance(query[field], str) or not query[field] for field in query_fields):
        raise ValueError("portal manifest query values must be non-empty strings")

    configurations = decode_configuration_records(value.get("configurations"))
    sequence_files = decode_sequence_file_records(value.get("sequence_files"))
    linked_accessions = {
        accession
        for record in configurations
        for accession in extract_accessions(record.seqspec_of)
    }
    unexpected_sequence_files = sorted(set(sequence_files) - linked_accessions)
    if unexpected_sequence_files:
        raise ValueError(
            "portal manifest contains unlinked sequence files: "
            + ", ".join(unexpected_sequence_files)
        )
    expected_missing = sorted(linked_accessions - set(sequence_files))
    missing = value.get("missing_linked_sequence_file_accessions")
    if missing != expected_missing:
        raise ValueError("portal manifest missing-sequence-file inventory does not reconcile")

    counts = value.get("counts")
    expected_counts = {
        "configurations": len(configurations),
        "linked_sequence_files": len(sequence_files),
        "missing_linked_sequence_files": len(expected_missing),
    }
    if counts != expected_counts:
        raise ValueError("portal manifest counts do not reconcile")
    return configurations, sequence_files, value


def decode_configuration_records(value: Any) -> list[ConfigurationRecord]:
    if not isinstance(value, list) or not value:
        raise ValueError("portal manifest must contain at least one configuration")
    expected_fields = set(ConfigurationRecord.__dataclass_fields__)
    records = []
    for index, item in enumerate(value):
        if not isinstance(item, dict) or set(item) != expected_fields:
            raise ValueError(f"portal configuration {index} has invalid fields")
        list_fields = {"preferred_assay_titles", "aliases", "seqspec_of"}
        if any(
            not isinstance(item[field], list)
            or any(not isinstance(entry, str) for entry in item[field])
            for field in list_fields
        ):
            raise ValueError(f"portal configuration {index} has invalid list fields")
        if any(
            not isinstance(item[field], str)
            for field in expected_fields - list_fields
        ):
            raise ValueError(f"portal configuration {index} has invalid string fields")
        records.append(ConfigurationRecord(**item))
    accessions = [record.accession for record in records]
    if any(not accession for accession in accessions) or len(accessions) != len(
        set(accessions)
    ):
        raise ValueError("portal configuration accessions must be non-empty and unique")
    if accessions != sorted(accessions):
        raise ValueError("portal configurations must be sorted by accession")
    return records


def decode_sequence_file_records(value: Any) -> dict[str, SequenceFileRecord]:
    if not isinstance(value, list):
        raise ValueError("portal manifest sequence_files must be a list")
    expected_fields = set(SequenceFileRecord.__dataclass_fields__)
    records = []
    for index, item in enumerate(value):
        if not isinstance(item, dict) or set(item) != expected_fields:
            raise ValueError(f"portal sequence file {index} has invalid fields")
        if (
            not isinstance(item["accession"], str)
            or not isinstance(item["href"], str)
            or not isinstance(item["controlled_access"], bool)
            or not isinstance(item["read_names"], list)
            or any(not isinstance(entry, str) for entry in item["read_names"])
        ):
            raise ValueError(f"portal sequence file {index} has invalid values")
        records.append(SequenceFileRecord(**item))
    accessions = [record.accession for record in records]
    if any(not accession for accession in accessions) or len(accessions) != len(
        set(accessions)
    ):
        raise ValueError("portal sequence-file accessions must be non-empty and unique")
    if accessions != sorted(accessions):
        raise ValueError("portal sequence files must be sorted by accession")
    return {record.accession: record for record in records}


def load_frozen_sampling_policy(path: Path) -> dict[str, Any]:
    value = load_json_object(path, "sampling policy")
    if value.get("schema_version") != SAMPLING_POLICY_SCHEMA_VERSION:
        raise ValueError(
            "unsupported sampling policy schema: "
            f"{value.get('schema_version', '')}"
        )
    stable = sampling_policy_stable(value)
    if value.get("policy_id") != sha256_json(stable)[:16]:
        raise ValueError("sampling policy ID does not reconcile with its contents")
    if value.get("frozen") is not True:
        raise ValueError("sampling policy is not frozen")
    default = value.get("default")
    if not isinstance(default, dict):
        raise ValueError("sampling policy has no default rule")
    n_reads = default.get("records_per_fastq")
    if not isinstance(n_reads, int) or isinstance(n_reads, bool) or n_reads <= 0:
        raise ValueError("sampling policy records_per_fastq must be a positive integer")
    method = default.get("sampling_method")
    if method != SAMPLING_METHOD:
        raise ValueError(
            "IGVF remote audit currently supports only prefix sampling; "
            f"the frozen policy selected {method!r}"
        )
    if default.get("sampling_seed") is not None:
        raise ValueError("prefix sampling policy must use a null sampling seed")
    return value


def load_audit_protocol(path: Path) -> dict[str, Any]:
    value = load_json_object(path, "audit protocol")
    if value.get("schema_version") != AUDIT_PROTOCOL_SCHEMA_VERSION:
        raise ValueError(
            f"unsupported audit protocol schema: {value.get('schema_version', '')}"
        )
    expected_fields = {
        "schema_version",
        "experiment_id",
        "execution",
        "repeatability",
        "targets",
    }
    if set(value) != expected_fields or value.get("experiment_id") != (
        "current-igvf-audit"
    ):
        raise ValueError("audit protocol fields or experiment ID are invalid")

    execution = value.get("execution")
    execution_fields = {
        "command_timeout_seconds",
        "max_attempts",
        "retry_delays_seconds",
        "retryable_failure_categories",
    }
    if not isinstance(execution, dict) or set(execution) != execution_fields:
        raise ValueError("audit protocol execution rules are incomplete")
    timeout = execution["command_timeout_seconds"]
    attempts = execution["max_attempts"]
    delays = execution["retry_delays_seconds"]
    categories = execution["retryable_failure_categories"]
    if not isinstance(timeout, int) or isinstance(timeout, bool) or timeout <= 0:
        raise ValueError("audit protocol timeout must be a positive integer")
    if not isinstance(attempts, int) or isinstance(attempts, bool) or attempts <= 0:
        raise ValueError("audit protocol max_attempts must be a positive integer")
    if (
        not isinstance(delays, list)
        or len(delays) != attempts - 1
        or any(
            isinstance(delay, bool)
            or not isinstance(delay, (int, float))
            or delay < 0
            for delay in delays
        )
    ):
        raise ValueError(
            "audit protocol needs one nonnegative retry delay per retry"
        )
    if categories != ["transport"]:
        raise ValueError("audit protocol may retry only transport failures")

    repeatability = value.get("repeatability")
    repeatability_fields = {
        "subset_run_count",
        "selection_seed",
        "canonical_excluded_top_level_fields",
    }
    if not isinstance(repeatability, dict) or set(repeatability) != repeatability_fields:
        raise ValueError("audit protocol repeatability rules are incomplete")
    subset_count = repeatability["subset_run_count"]
    seed = repeatability["selection_seed"]
    excluded = repeatability["canonical_excluded_top_level_fields"]
    if (
        not isinstance(subset_count, int)
        or isinstance(subset_count, bool)
        or subset_count <= 0
        or not isinstance(seed, int)
        or isinstance(seed, bool)
        or not isinstance(excluded, list)
        or not excluded
        or any(not isinstance(field, str) or not field for field in excluded)
    ):
        raise ValueError("audit protocol repeatability values are invalid")

    targets = value.get("targets")
    if not isinstance(targets, dict) or set(targets) != {
        "minimum_completion_fraction",
        "maximum_unclassified_failure_count",
        "required_repeat_match_fraction",
    }:
        raise ValueError("audit protocol targets are incomplete")
    fractions = (
        targets["minimum_completion_fraction"],
        targets["required_repeat_match_fraction"],
    )
    if any(
        isinstance(fraction, bool)
        or not isinstance(fraction, (int, float))
        or not 0 <= fraction <= 1
        for fraction in fractions
    ) or (
        isinstance(targets["maximum_unclassified_failure_count"], bool)
        or targets["maximum_unclassified_failure_count"] != 0
    ):
        raise ValueError("audit protocol target values are invalid")
    return value


def audit_protocol_id(value: dict[str, Any]) -> str:
    return sha256_json(value)[:16]


def command_output(command: list[str], env: dict[str, str] | None = None) -> str:
    result = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    output = result.stdout.strip() or result.stderr.strip()
    return output.splitlines()[0] if output else ""


def git_output(repo_root: Path, *args: str) -> str:
    if not (repo_root / ".git").exists():
        return ""
    return command_output(["git", "-C", str(repo_root), *args])


def build_tool_identity(command: ToolCommand, repo_root: Path) -> ToolIdentity:
    env = os.environ.copy()
    env.update(command.env)
    version = command_output(command.argv + ["--version"], env=env)
    git_commit = git_output(repo_root, "rev-parse", "HEAD")
    runtime_paths = [
        value
        for value in ("Cargo.toml", "Cargo.lock", "src")
        if (repo_root / value).exists()
    ]
    status = git_output(
        repo_root,
        "status",
        "--porcelain",
        "--untracked-files=normal",
        "--",
        *runtime_paths,
    )

    executable = Path(command.argv[0]).resolve()
    executable_sha256 = file_sha256(executable) if executable.is_file() else ""

    return ToolIdentity(
        command=command.argv,
        version=version,
        git_commit=git_commit,
        git_dirty=bool(status),
        runtime_source_sha256=path_tree_sha256(repo_root, runtime_paths),
        executable_sha256=executable_sha256,
    )


def write_study_manifest(
    output_root: Path,
    args: argparse.Namespace,
    configurations: list[ConfigurationRecord],
    sequence_files: dict[str, SequenceFileRecord],
    seqcheck_identity: ToolIdentity,
    seqspec_identity: ToolIdentity,
    seqcheck_auth_profile: str | None,
    seqspec_auth_profile: str | None,
    started_at: str,
    *,
    portal_manifest: dict[str, Any] | None = None,
    sampling_policy: dict[str, Any] | None = None,
    audit_protocol: dict[str, Any] | None = None,
) -> AuditContext:
    linked_accessions = {
        accession
        for record in configurations
        for accession in extract_accessions(record.seqspec_of)
    }
    linked_sequence_files = [
        asdict(sequence_files[accession])
        for accession in sorted(linked_accessions)
        if accession in sequence_files
    ]
    sampling_method = (
        sampling_policy["default"]["sampling_method"]
        if sampling_policy is not None
        else SAMPLING_METHOD
    )
    sampling_seed = (
        sampling_policy["default"]["sampling_seed"]
        if sampling_policy is not None
        else None
    )
    execution = (
        audit_protocol["execution"]
        if audit_protocol is not None
        else {
            "command_timeout_seconds": None,
            "max_attempts": 1,
            "retry_delays_seconds": [],
            "retryable_failure_categories": ["transport"],
        }
    )
    protocol_id = audit_protocol_id(audit_protocol) if audit_protocol else ""
    frozen_inputs, input_provenance = build_frozen_input_identities(
        args,
        portal_manifest,
        sampling_policy,
        audit_protocol,
    )
    audit_runner_identity = build_script_identity(
        Path(__file__).resolve(), version=AUDIT_SCHEMA_VERSION
    )
    stable_manifest = {
        "audit_schema_version": AUDIT_SCHEMA_VERSION,
        "sampling_method": sampling_method,
        "sampling_seed": sampling_seed,
        "execution": execution,
        "query": {
            "api_root": args.api_root,
            "portal_root": args.portal_root,
            "status": args.status,
            "upload_status": args.upload_status,
            "public_only": args.public_only,
            "n_reads_per_fastq": args.n_reads,
        },
        "configurations": [asdict(record) for record in configurations],
        "sequence_files": linked_sequence_files,
        "auth": {
            "requested_profile": args.auth_profile,
            "seqcheck_profile_ready": seqcheck_auth_profile is not None,
            "seqspec_profile_ready": seqspec_auth_profile is not None,
        },
        "frozen_inputs": frozen_inputs,
        "tools": {
            "audit_runner": audit_runner_identity,
            "seqcheck": asdict(seqcheck_identity),
            "seqspec": asdict(seqspec_identity),
            "fastqc": build_external_tool_identity(
                args.fastqc_command,
                Path(args.fastqc_limits).resolve() if args.fastqc_limits else None,
            ),
        },
    }
    run_identity = {
        **stable_manifest,
        "tools": {
            "audit_runner": functional_script_identity(audit_runner_identity),
            "seqcheck": functional_tool_identity(seqcheck_identity),
            "seqspec": functional_tool_identity(seqspec_identity),
            "fastqc": functional_external_tool_identity(
                stable_manifest["tools"]["fastqc"]
            ),
        },
    }
    run_id = sha256_json(run_identity)[:16]
    context = AuditContext(
        run_id=run_id,
        sampling_method=sampling_method,
        sampling_seed=sampling_seed,
        seqcheck=seqcheck_identity,
        seqspec=seqspec_identity,
        command_timeout_seconds=execution["command_timeout_seconds"],
        max_attempts=execution["max_attempts"],
        retry_delays_seconds=tuple(execution["retry_delays_seconds"]),
        retryable_failure_categories=tuple(
            execution["retryable_failure_categories"]
        ),
        audit_protocol_id=protocol_id,
    )
    manifest = {
        **stable_manifest,
        "run_id": run_id,
        "started_at": started_at,
        "invocation": sys.argv,
        "workers": args.workers,
        "force": args.force,
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "selected_configuration_count": len(configurations),
        "linked_sequence_file_count": len(linked_sequence_files),
        "inputs": input_provenance,
    }
    write_json(output_root / "manifests" / "study.json", manifest)
    return context


def build_frozen_input_identities(
    args: argparse.Namespace,
    portal_manifest: dict[str, Any] | None,
    sampling_policy: dict[str, Any] | None,
    audit_protocol: dict[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    stable: dict[str, Any] = {}
    provenance: dict[str, Any] = {}
    portal_path = getattr(args, "portal_manifest", None)
    if portal_manifest is not None:
        if portal_path is None:
            raise ValueError("portal manifest path is required for provenance")
        portal_path = Path(portal_path).resolve()
        portal_stable_sha256 = sha256_json(portal_manifest_stable(portal_manifest))
        stable["portal_manifest"] = {
            "portal_manifest_id": portal_manifest["portal_manifest_id"],
            "stable_sha256": portal_stable_sha256,
        }
        provenance["portal_manifest"] = {
            **file_identity(portal_path),
            **stable["portal_manifest"],
        }
    else:
        stable["portal_manifest"] = {"mode": "live"}
        provenance["portal_manifest"] = {"mode": "live"}

    policy_path = getattr(args, "sampling_policy", None)
    if sampling_policy is not None:
        if policy_path is None:
            raise ValueError("sampling policy path is required for provenance")
        policy_path = Path(policy_path).resolve()
        policy_stable_sha256 = sha256_json(sampling_policy_stable(sampling_policy))
        stable["sampling_policy"] = {
            "policy_id": sampling_policy["policy_id"],
            "stable_sha256": policy_stable_sha256,
        }
        provenance["sampling_policy"] = {
            **file_identity(policy_path),
            **stable["sampling_policy"],
        }
    else:
        stable["sampling_policy"] = {"mode": "command_line"}
        provenance["sampling_policy"] = {"mode": "command_line"}

    protocol_path = getattr(args, "audit_protocol", None)
    if audit_protocol is not None:
        if protocol_path is None:
            raise ValueError("audit protocol path is required for provenance")
        protocol_path = Path(protocol_path).resolve()
        stable["audit_protocol"] = {
            "audit_protocol_id": audit_protocol_id(audit_protocol),
            "stable_sha256": sha256_json(audit_protocol),
        }
        provenance["audit_protocol"] = {
            **file_identity(protocol_path),
            **stable["audit_protocol"],
        }
    else:
        stable["audit_protocol"] = {"mode": "command_line_defaults"}
        provenance["audit_protocol"] = {"mode": "command_line_defaults"}
    return stable, provenance


def build_external_tool_identity(
    command: str, configuration_path: Path | None
) -> dict[str, Any]:
    resolved = shutil.which(command)
    configuration = {
        "mode": "custom" if configuration_path else "packaged_defaults",
        "path": str(configuration_path) if configuration_path else "",
        "sha256": file_sha256(configuration_path) if configuration_path else "",
    }
    return {
        "command": command,
        "resolved_path": resolved or "",
        "available": resolved is not None,
        "version": command_output([resolved, "--version"]) if resolved else "",
        "executable_sha256": file_sha256(Path(resolved)) if resolved else "",
        "configuration": configuration,
    }


def build_script_identity(path: Path, *, version: str) -> dict[str, Any]:
    path = path.resolve()
    repo_root = path.parents[1]
    status = git_output(
        repo_root,
        "status",
        "--porcelain",
        "--untracked-files=normal",
        "--",
        str(path.relative_to(repo_root)),
    )
    return {
        **file_identity(path),
        "version": version,
        "git_commit": git_output(repo_root, "rev-parse", "HEAD"),
        "git_dirty": bool(status),
        "python": sys.version.split()[0],
    }


def functional_script_identity(identity: dict[str, Any]) -> dict[str, Any]:
    return {
        "version": identity["version"],
        "sha256": identity["sha256"],
        "python": identity["python"],
    }


def functional_tool_identity(identity: ToolIdentity) -> dict[str, Any]:
    return {
        "command": identity.command,
        "version": identity.version,
        "runtime_source_sha256": identity.runtime_source_sha256,
        "executable_sha256": identity.executable_sha256,
    }


def functional_external_tool_identity(identity: dict[str, Any]) -> dict[str, Any]:
    return {
        "command": identity.get("command", ""),
        "version": identity.get("version", ""),
        "executable_sha256": identity.get("executable_sha256", ""),
        "configuration": identity.get("configuration", {}),
    }


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def path_tree_sha256(root: Path, relative_paths: list[str]) -> str:
    files = []
    for relative_path in relative_paths:
        path = root / relative_path
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            files.extend(item for item in path.rglob("*") if item.is_file())

    digest = hashlib.sha256()
    for path in sorted(files):
        digest.update(str(path.relative_to(root)).encode("utf-8"))
        digest.update(b"\0")
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest.update(b"\0")
    return digest.hexdigest() if files else ""


def build_report_cache_key(
    audit_context: AuditContext,
    record: ConfigurationRecord,
    modality: str,
    raw_seqspec_version: str,
    normalized_seqspec_version: str,
    spec_url: str,
    expected_files: list[dict[str, Any]],
    sequence_records: list[SequenceFileRecord],
    fastq_urls: list[str],
    region_annotations: dict[str, dict[str, Any]],
    n_reads: int,
    controlled_access: bool,
) -> str:
    payload = {
        "audit_schema_version": AUDIT_SCHEMA_VERSION,
        "sampling_method": audit_context.sampling_method,
        "sampling_seed": audit_context.sampling_seed,
        "configuration": asdict(record),
        "modality": modality,
        "raw_seqspec_version": raw_seqspec_version,
        "normalized_seqspec_version": normalized_seqspec_version,
        "spec_url": spec_url,
        "expected_files": expected_files,
        "sequence_records": [asdict(item) for item in sequence_records],
        "fastq_urls": fastq_urls,
        "region_annotations": region_annotations,
        "n_reads": n_reads,
        "controlled_access": controlled_access,
        "tools": {
            "seqcheck": functional_tool_identity(audit_context.seqcheck),
            "seqspec": functional_tool_identity(audit_context.seqspec),
        },
    }
    return sha256_json(payload)


def report_cache_matches(report: dict[str, Any], expected_cache_key: str) -> bool:
    summary = report.get("audit_summary", {})
    return (
        summary.get("audit_schema_version") == AUDIT_SCHEMA_VERSION
        and summary.get("cache_key") == expected_cache_key
    )


def fetch_configuration_records(
    api_root: str, status: str, upload_status: str
) -> list[ConfigurationRecord]:
    params = [
        ("type", "ConfigurationFile"),
        ("status", status),
        ("upload_status", upload_status),
        ("limit", "all"),
        ("format", "json"),
        ("field", "accession"),
        ("field", "href"),
        ("field", "lab.title"),
        ("field", "submitted_by.title"),
        ("field", "award.component"),
        ("field", "file_set.accession"),
        ("field", "file_set.assay_term.term_name"),
        ("field", "preferred_assay_titles"),
        ("field", "aliases"),
        ("field", "seqspec_of"),
        ("field", "status"),
        ("field", "upload_status"),
    ]
    payload = fetch_json(
        absolute_url(api_root, f"search/?{urllib.parse.urlencode(params, doseq=True)}")
    )
    return [normalize_configuration_record(item) for item in payload.get("@graph", [])]


def fetch_sequence_file_records(
    api_root: str, status: str, upload_status: str
) -> dict[str, SequenceFileRecord]:
    params = [
        ("type", "SequenceFile"),
        ("status", status),
        ("upload_status", upload_status),
        ("file_format", "fastq"),
        ("limit", "all"),
        ("format", "json"),
        ("field", "accession"),
        ("field", "href"),
        ("field", "controlled_access"),
        ("field", "read_names"),
    ]
    payload = fetch_json(
        absolute_url(api_root, f"search/?{urllib.parse.urlencode(params, doseq=True)}")
    )
    records = {}
    for item in payload.get("@graph", []):
        record = normalize_sequence_file_record(item)
        if record is None:
            continue
        records[record.accession] = record
    return records


def normalize_sequence_file_record(item: dict[str, Any]) -> SequenceFileRecord | None:
    accession = str(item.get("accession", "")).strip()
    if not accession:
        return None
    return SequenceFileRecord(
        accession=accession,
        href=str(item.get("href", "")).strip(),
        controlled_access=bool(item.get("controlled_access", False)),
        read_names=[str(value) for value in item.get("read_names", [])],
    )


def normalize_configuration_record(item: dict[str, Any]) -> ConfigurationRecord:
    return ConfigurationRecord(
        accession=str(item.get("accession", "")).strip(),
        href=str(item.get("href", "")).strip(),
        lab=nested_string(item, "lab", "title"),
        submitted_by=nested_string(item, "submitted_by", "title"),
        award_component=nested_string(item, "award", "component"),
        file_set_accession=nested_string(item, "file_set", "accession"),
        assay_term=nested_string(item, "file_set", "assay_term", "term_name"),
        preferred_assay_titles=[str(value) for value in item.get("preferred_assay_titles", [])],
        aliases=[str(value) for value in item.get("aliases", [])],
        seqspec_of=[str(value) for value in item.get("seqspec_of", [])],
        status=str(item.get("status", "")).strip(),
        upload_status=str(item.get("upload_status", "")).strip(),
    )


def nested_string(item: dict[str, Any], *path: str) -> str:
    current: Any = item
    for key in path:
        if not isinstance(current, dict):
            return ""
        current = current.get(key)
    if current is None:
        return ""
    return str(current)


def extract_accessions(paths: list[str]) -> list[str]:
    accessions = []
    for value in paths:
        parts = [piece for piece in str(value).strip("/").split("/") if piece]
        if parts:
            accessions.append(parts[-1])
    return accessions


def fetch_json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )
    with urllib.request.urlopen(request) as response:
        return json.load(io.TextIOWrapper(response, encoding="utf-8"))


def request_headers(auth_header: str | None, accept_json: bool) -> dict[str, str]:
    headers = {"User-Agent": USER_AGENT}
    if accept_json:
        headers["Accept"] = "application/json"
    if auth_header:
        headers["Authorization"] = auth_header
    return headers


def discover_seqcheck_command(repo_root: Path) -> ToolCommand:
    env_bin = os.environ.get("SEQCHECK_BIN")
    if env_bin:
        return ToolCommand([env_bin], {})
    local_bin = repo_root / "target" / "debug" / "seqcheck"
    if local_bin.exists():
        return ToolCommand([str(local_bin)], {})
    path_bin = shutil.which("seqcheck")
    if path_bin:
        return ToolCommand([path_bin], {})
    raise RuntimeError(
        "Could not find seqcheck binary. Build it first or set SEQCHECK_BIN."
    )


def discover_seqspec_command(workspace_root: Path) -> ToolCommand:
    env_bin = os.environ.get("SEQSPEC_BIN")
    if env_bin:
        return ToolCommand([env_bin], {})
    sibling_bin = workspace_root / "seqspec" / "target" / "debug" / "seqspec"
    if sibling_bin.exists():
        return ToolCommand([str(sibling_bin)], {})
    path_bin = shutil.which("seqspec")
    if path_bin:
        return ToolCommand([path_bin], {})
    sibling_repo = workspace_root / "seqspec"
    if sibling_repo.exists():
        return ToolCommand(
            [sys.executable, "-m", "seqspec.main"],
            {"PYTHONPATH": str(sibling_repo)},
        )
    raise RuntimeError(
        "Could not find seqspec command. Build/install seqspec first or set SEQSPEC_BIN."
    )


def get_seqspec_file_version(
    command: ToolCommand,
    spec_source: str,
    auth_profile: str | None,
    mpl_dir: Path,
    audit_context: AuditContext,
) -> tuple[str, int]:
    env = os.environ.copy()
    env.update(command.env)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    argv = command.argv + ["version", spec_source]
    if auth_profile:
        argv.extend(["--auth-profile", auth_profile])
    result, attempts = run_audit_command(
        argv,
        env=env,
        audit_context=audit_context,
        stage="seqspec_version",
        reason="seqspec_version_error",
    )
    return parse_seqspec_version_output(result.stdout), attempts


def list_modality_files(
    command: ToolCommand,
    spec_source: str,
    modality: str,
    auth_profile: str | None,
    mpl_dir: Path,
    audit_context: AuditContext,
) -> tuple[list[dict[str, Any]], int]:
    env = os.environ.copy()
    env.update(command.env)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    argv = command.argv + [
        "file",
        "-m",
        modality,
        "-f",
        "json",
        "-s",
        "read",
        "-k",
        "all",
        spec_source,
    ]
    if auth_profile:
        argv.extend(["--auth-profile", auth_profile])
    result, attempts = run_audit_command(
        argv,
        env=env,
        audit_context=audit_context,
        stage="enumerate_files",
        reason="seqspec_file_error",
    )
    return json.loads(result.stdout), attempts


def list_modalities(
    command: ToolCommand,
    spec_source: str,
    auth_profile: str | None,
    mpl_dir: Path,
    audit_context: AuditContext,
) -> tuple[list[str], int]:
    env = os.environ.copy()
    env.update(command.env)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    argv = command.argv + ["info", "-k", "modalities", "-f", "json", spec_source]
    if auth_profile:
        argv.extend(["--auth-profile", auth_profile])
    result, attempts = run_audit_command(
        argv,
        env=env,
        audit_context=audit_context,
        stage="enumerate_modalities",
        reason="seqspec_info_error",
    )
    payload = json.loads(result.stdout)
    return [str(item) for item in payload], attempts


def list_modality_region_annotations(
    command: ToolCommand,
    spec_source: str,
    modality: str,
    auth_profile: str | None,
    mpl_dir: Path,
    audit_context: AuditContext,
) -> tuple[dict[str, dict[str, Any]], int]:
    env = os.environ.copy()
    env.update(command.env)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    argv = command.argv + [
        "info",
        "-k",
        "library_spec",
        "-f",
        "json",
        spec_source,
    ]
    if auth_profile:
        argv.extend(["--auth-profile", auth_profile])
    result, attempts = run_audit_command(
        argv,
        env=env,
        audit_context=audit_context,
        stage="enumerate_regions",
        reason="seqspec_region_error",
    )
    payload = json.loads(result.stdout)
    if not isinstance(payload, dict) or not isinstance(payload.get(modality), list):
        raise ValueError(f"seqspec has no region list for modality: {modality}")
    return index_region_annotations(payload[modality]), attempts


def index_region_annotations(
    regions: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    indexed = {}
    for region in regions:
        if not isinstance(region, dict):
            raise ValueError("seqspec modality contains a non-object region")
        region_id = str(region.get("region_id", "")).strip()
        sequence_type = str(region.get("sequence_type", "")).strip()
        if not region_id or not sequence_type:
            raise ValueError("seqspec region annotation is incomplete")
        if region_id in indexed:
            raise ValueError(f"seqspec region_id is not unique in modality: {region_id}")
        region_type = region.get("region_type")
        terms = [region_type] if isinstance(region_type, str) else region_type
        if (
            not isinstance(terms, list)
            or not terms
            or any(
                not isinstance(value, str) or not value.startswith("RGN:")
                for value in terms
            )
        ):
            raise ValueError(f"seqspec region_type is invalid for region: {region_id}")
        if len(terms) != len(set(terms)):
            raise ValueError(f"seqspec region_type contains duplicates: {region_id}")
        indexed[region_id] = {
            "region_id": region_id,
            "sequence_type": sequence_type,
            "ontology_terms": sorted(terms),
        }
    return indexed


def resolve_expected_sequence_file(
    file_record: dict[str, Any], sequence_files: dict[str, SequenceFileRecord]
) -> SequenceFileRecord | None:
    file_id = str(file_record.get("file_id", "")).strip()
    if file_id and file_id in sequence_files:
        return sequence_files[file_id]

    filename = str(file_record.get("filename", "")).strip()
    if filename.endswith(".fastq.gz"):
        accession = filename[: -len(".fastq.gz")]
        if accession in sequence_files:
            return sequence_files[accession]
    return None


def is_fastq_expectation(file_record: dict[str, Any]) -> bool:
    filetype = str(file_record.get("filetype", "")).strip().lower().lstrip(".")
    if filetype in {"fastq", "fq", "fastq.gz", "fq.gz"}:
        return True

    filename = str(file_record.get("filename", "")).strip().lower()
    return filename.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))


def run_seqcheck(
    command: ToolCommand,
    spec: str,
    modality: str,
    n_reads: int,
    fastqs: list[str],
    report_path: Path,
    auth_profile: str | None,
    audit_context: AuditContext,
) -> int:
    env = os.environ.copy()
    env.update(command.env)
    argv = command.argv + [
        "check",
        "-s",
        spec,
        "-m",
        modality,
        "-n",
        str(n_reads),
        "--format",
        "json",
        "-o",
        str(report_path),
    ]
    if auth_profile:
        argv.extend(["--auth-profile", auth_profile])
    argv.extend(fastqs)
    _, attempts = run_audit_command(
        argv,
        env=env,
        audit_context=audit_context,
        stage="seqcheck",
        reason="seqcheck_error",
    )
    return attempts


def run_audit_command(
    argv: list[str],
    *,
    env: dict[str, str],
    audit_context: AuditContext,
    stage: str,
    reason: str,
) -> tuple[subprocess.CompletedProcess[str], int]:
    for attempt in range(1, audit_context.max_attempts + 1):
        try:
            completed = subprocess.run(
                argv,
                check=True,
                capture_output=True,
                text=True,
                env=env,
                timeout=audit_context.command_timeout_seconds,
            )
        except subprocess.TimeoutExpired as timeout:
            error = subprocess.CalledProcessError(
                returncode=124,
                cmd=argv,
                output=timeout.stdout or "",
                stderr=(
                    f"command timed out after "
                    f"{audit_context.command_timeout_seconds} seconds"
                ),
            )
        except subprocess.CalledProcessError as called:
            error = called
        else:
            return completed, attempt

        message = stderr_message(error)
        category = classify_failure_category(stage, reason, message)
        should_retry = (
            attempt < audit_context.max_attempts
            and category in audit_context.retryable_failure_categories
        )
        if should_retry:
            time.sleep(audit_context.retry_delays_seconds[attempt - 1])
            continue
        error.attempt_count = attempt
        raise error
    raise AssertionError("unreachable retry loop")


def command_attempt_count(error: subprocess.CalledProcessError) -> int:
    value = getattr(error, "attempt_count", 1)
    return value if isinstance(value, int) and value > 0 else 1


def resolve_ready_auth_profile(
    profile_name: str | None,
    *,
    config_env: str,
    config_subdir: str,
) -> str | None:
    if not profile_name:
        return None
    config_path = resolve_auth_config_path(config_env, config_subdir)
    if config_path is None or not config_path.exists():
        return None
    try:
        config = tomllib.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError):
        return None
    profile = config.get("profiles", {}).get(profile_name)
    if not isinstance(profile, dict):
        return None
    username_env = str(profile.get("username_env", "")).strip()
    password_env = str(profile.get("password_env", "")).strip()
    if username_env and password_env:
        if os.environ.get(username_env) and os.environ.get(password_env):
            return profile_name
    return None


def resolve_auth_config_path(config_env: str, config_subdir: str) -> Path | None:
    config_override = os.environ.get(config_env)
    candidate_paths = []
    if config_override:
        candidate_paths.append(Path(config_override))

    xdg_home = os.environ.get("XDG_CONFIG_HOME")
    if xdg_home:
        candidate_paths.append(Path(xdg_home) / config_subdir / "auth.toml")

    candidate_paths.append(Path.home() / ".config" / config_subdir / "auth.toml")

    for path in candidate_paths:
        if path.exists():
            return path
    return None


def parse_seqspec_version_output(output: str) -> str:
    for line in output.splitlines():
        if line.startswith("seqspec file version:"):
            return line.split(":", 1)[1].strip()
    return ""


def normalize_seqspec_version(raw_version: str) -> str:
    if raw_version in {"0.0.0", "0.1.0", "0.1.1", "0.2.0", "0.3.0", "0.4.0"}:
        return CURRENT_SEQSPEC_VERSION
    return raw_version


def flatten_report(
    report: dict[str, Any],
    record: ConfigurationRecord,
    modality: str,
    raw_seqspec_version: str,
    normalized_seqspec_version: str,
    report_path: Path,
    modality_count: int,
    expected_fastq_count: int,
    supplied_fastq_count: int,
    controlled_access: bool,
    run_status: str,
    region_annotations: dict[str, dict[str, Any]] | None = None,
    access_class: str | None = None,
) -> tuple[RunRecord, list[DiagnosticRecord], list[MetricRecord]]:
    if expected_fastq_count <= 0:
        expected_fastq_count = infer_expected_fastq_count(report)
    if supplied_fastq_count <= 0:
        supplied_fastq_count = infer_supplied_fastq_count(report)
    if modality_count <= 0:
        modality_count = int(report.get("audit_summary", {}).get("modality_count", 0))

    assessments = [
        assessment
        for result in report.get("results", [])
        for assessment in result.get("assessment", [])
    ]
    counts = count_assessment_types(assessments)
    requested_reads = int(report.get("meta", {}).get("requested_reads", 0))
    audit_summary = report.get("audit_summary", {})
    study_run_id = str(audit_summary.get("study_run_id", ""))
    cache_key = str(audit_summary.get("cache_key", ""))
    sampling_method = str(audit_summary.get("sampling_method", SAMPLING_METHOD))
    sampling_seed = audit_summary.get("sampling_seed")
    sampled_record_count = infer_sampled_record_count(report)
    resolved_access_class = access_class or (
        "controlled" if controlled_access else "public"
    )

    run = RunRecord(
        study_run_id=study_run_id,
        cache_key=cache_key,
        sampling_method=sampling_method,
        sampling_seed=sampling_seed,
        configuration_accession=record.accession,
        modality=modality,
        modality_count=modality_count,
        lab=record.lab,
        submitted_by=record.submitted_by,
        award_component=record.award_component,
        file_set_accession=record.file_set_accession,
        assay_term=record.assay_term,
        preferred_assay_titles=";".join(record.preferred_assay_titles),
        aliases=";".join(record.aliases),
        raw_seqspec_version=raw_seqspec_version,
        normalized_seqspec_version=normalized_seqspec_version,
        requested_reads=requested_reads,
        requested_reads_total=(
            requested_reads * supplied_fastq_count if requested_reads > 0 else 0
        ),
        sampled_record_count=sampled_record_count,
        expected_fastq_count=expected_fastq_count,
        supplied_fastq_count=supplied_fastq_count,
        controlled_access=controlled_access,
        access_class=resolved_access_class,
        report_path=str(report_path),
        run_status=run_status,
        pass_count=counts["pass"],
        warning_count=counts["warning"],
        error_count=counts["error"],
        interpretation_count=counts["interpretation"],
        max_attempt_count=int(audit_summary.get("max_attempt_count", 0)),
    )

    diagnostics = []
    metric_rows = []
    for result_index, result in enumerate(report.get("results", [])):
        ontology_terms = ";".join(extract_ontology_terms(result))
        annotations = resolve_result_region_annotations(result, region_annotations)
        sequence_types = ";".join(
            sorted({value["sequence_type"] for value in annotations})
        )
        region_annotations_json = canonical_json(annotations)
        for assessment in result.get("assessment", []):
            diagnostics.append(
                DiagnosticRecord(
                    study_run_id=study_run_id,
                    cache_key=cache_key,
                    configuration_accession=record.accession,
                    modality=modality,
                    lab=record.lab,
                    submitted_by=record.submitted_by,
                    award_component=record.award_component,
                    file_set_accession=record.file_set_accession,
                    assay_term=record.assay_term,
                    preferred_assay_titles=";".join(record.preferred_assay_titles),
                    aliases=";".join(record.aliases),
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    report_path=str(report_path),
                    check=str(result.get("check", "")),
                    assessment_type=str(assessment.get("type", "")),
                    assessment_code=str(assessment.get("code", "")),
                    assessment_description=str(assessment.get("description", "")),
                    files=";".join(str(value) for value in result.get("files", [])),
                    reads=";".join(str(value) for value in result.get("reads", [])),
                    regions=";".join(str(value) for value in result.get("regions", [])),
                    ontology_terms=ontology_terms,
                    sequence_types=sequence_types,
                    region_annotations_json=region_annotations_json,
                    access_class=resolved_access_class,
                )
            )

        for metric_side in ("expected", "observed"):
            for metric in result.get(metric_side, []):
                data = metric.get("data", {})
                metric_rows.append(
                    MetricRecord(
                        study_run_id=study_run_id,
                        cache_key=cache_key,
                        configuration_accession=record.accession,
                        modality=modality,
                        lab=record.lab,
                        submitted_by=record.submitted_by,
                        award_component=record.award_component,
                        file_set_accession=record.file_set_accession,
                        assay_term=record.assay_term,
                        preferred_assay_titles=";".join(record.preferred_assay_titles),
                        aliases=";".join(record.aliases),
                        raw_seqspec_version=raw_seqspec_version,
                        normalized_seqspec_version=normalized_seqspec_version,
                        report_path=str(report_path),
                        result_index=result_index,
                        check=str(result.get("check", "")),
                        files=";".join(str(value) for value in result.get("files", [])),
                        reads=";".join(str(value) for value in result.get("reads", [])),
                        regions=";".join(str(value) for value in result.get("regions", [])),
                        ontology_terms=ontology_terms,
                        sequence_types=sequence_types,
                        region_annotations_json=region_annotations_json,
                        access_class=resolved_access_class,
                        metric_side=metric_side,
                        metric_id=str(metric.get("id", "")),
                        metric_name=str(metric.get("name", "")),
                        metric_description=str(metric.get("description", "")),
                        data_kind=str(data.get("kind", "")),
                        unit=str(data.get("unit") or ""),
                        value_json=canonical_json(data.get("value")),
                    )
                )
    return run, diagnostics, metric_rows


def extract_ontology_terms(value: Any) -> list[str]:
    terms: set[str] = set()

    def visit(item: Any) -> None:
        if isinstance(item, dict):
            for child in item.values():
                visit(child)
        elif isinstance(item, list):
            for child in item:
                visit(child)
        elif isinstance(item, str):
            terms.update(ONTOLOGY_TERM_PATTERN.findall(item))

    visit(value.get("expected", []) if isinstance(value, dict) else value)
    return sorted(terms)


def resolve_result_region_annotations(
    result: dict[str, Any],
    region_annotations: dict[str, dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    if region_annotations is None:
        return []
    region_ids = result.get("regions", [])
    if not isinstance(region_ids, list):
        raise ValueError("seqcheck result regions is not a list")
    resolved = []
    for value in region_ids:
        region_id = str(value).strip()
        if not region_id or region_id not in region_annotations:
            raise ValueError(f"seqcheck result references unknown region: {region_id}")
        resolved.append(region_annotations[region_id])
    return resolved


def infer_sampled_record_count(report: dict[str, Any]) -> int:
    sampled_by_file: dict[str, set[int]] = {}
    for result in report.get("results", []):
        files = result.get("files", [])
        if not isinstance(files, list) or len(files) != 1:
            continue
        for metric in result.get("observed", []):
            if metric.get("name") != "sampled_count":
                continue
            data = metric.get("data", {})
            value = data.get("value") if isinstance(data, dict) else None
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError("seqcheck sampled_count is not a nonnegative integer")
            sampled_by_file.setdefault(str(files[0]), set()).add(value)
    inconsistent = {
        file_id: sorted(values)
        for file_id, values in sampled_by_file.items()
        if len(values) != 1
    }
    if inconsistent:
        raise ValueError(f"seqcheck sampled_count differs within files: {inconsistent}")
    return sum(next(iter(values)) for values in sampled_by_file.values())


def infer_expected_fastq_count(report: dict[str, Any]) -> int:
    for result in report.get("results", []):
        if result.get("check") == "input_check":
            return len(result.get("files", []))
    return 0


def infer_supplied_fastq_count(report: dict[str, Any]) -> int:
    for result in report.get("results", []):
        if result.get("check") != "input_check":
            continue
        for observed in result.get("observed", []):
            if observed.get("name") == "matched_inputs":
                value = observed.get("data", {}).get("value", [])
                if isinstance(value, list):
                    return len(value)
    return infer_expected_fastq_count(report)


def annotate_report_summary(
    report: dict[str, Any],
    *,
    configuration_accession: str,
    modality: str,
    modality_count: int,
    expected_fastq_count: int,
    supplied_fastq_count: int,
    controlled_access: bool,
    audit_context: AuditContext,
    cache_key: str,
    access_class: str | None = None,
    max_attempt_count: int | None = None,
) -> dict[str, Any]:
    requested_reads = int(report.get("meta", {}).get("requested_reads", 0))
    sampled_record_count = infer_sampled_record_count(report)
    resolved_access_class = access_class or (
        "controlled" if controlled_access else "public"
    )
    previous_summary = report.get("audit_summary", {})
    resolved_max_attempt_count = (
        int(previous_summary.get("max_attempt_count", 1))
        if max_attempt_count is None
        else max_attempt_count
    )
    report["audit_summary"] = {
        "audit_schema_version": AUDIT_SCHEMA_VERSION,
        "study_run_id": audit_context.run_id,
        "cache_key": cache_key,
        "sampling_method": audit_context.sampling_method,
        "sampling_seed": audit_context.sampling_seed,
        "audit_protocol_id": audit_context.audit_protocol_id,
        "max_attempt_count": resolved_max_attempt_count,
        "audit_invocation": sys.argv,
        "tools": {
            "seqcheck": asdict(audit_context.seqcheck),
            "seqspec": asdict(audit_context.seqspec),
        },
        "configuration_accession": configuration_accession,
        "modality": modality,
        "modality_count": modality_count,
        "expected_fastq_count": expected_fastq_count,
        "supplied_fastq_count": supplied_fastq_count,
        "requested_reads_per_fastq": requested_reads,
        "requested_reads_total": (
            requested_reads * supplied_fastq_count if requested_reads > 0 else 0
        ),
        "sampled_record_count": sampled_record_count,
        "controlled_access": controlled_access,
        "access_class": resolved_access_class,
    }
    return report


def count_assessment_types(assessments: list[dict[str, Any]]) -> dict[str, int]:
    counts = {"pass": 0, "warning": 0, "error": 0, "interpretation": 0}
    for assessment in assessments:
        key = str(assessment.get("type", ""))
        if key in counts:
            counts[key] += 1
    return counts


def build_failure(
    study_run_id: str,
    record: ConfigurationRecord,
    modality: str,
    raw_seqspec_version: str,
    normalized_seqspec_version: str,
    stage: str,
    reason: str,
    message: str,
    access_class: str = "unknown",
    attempt_count: int = 1,
) -> FailureRecord:
    failure_category = classify_failure_category(stage, reason, message)
    return FailureRecord(
        study_run_id=study_run_id,
        configuration_accession=record.accession,
        modality=modality,
        lab=record.lab,
        submitted_by=record.submitted_by,
        award_component=record.award_component,
        file_set_accession=record.file_set_accession,
        assay_term=record.assay_term,
        preferred_assay_titles=";".join(record.preferred_assay_titles),
        aliases=";".join(record.aliases),
        raw_seqspec_version=raw_seqspec_version,
        normalized_seqspec_version=normalized_seqspec_version,
        stage=stage,
        reason=reason,
        message=message,
        access_class=access_class,
        failure_category=failure_category,
        attempt_count=attempt_count,
    )


def classify_failure_category(stage: str, reason: str, message: str) -> str:
    normalized = f"{reason} {message}".lower()
    if reason == "unhandled_exception":
        return "unclassified"
    if reason == "controlled_fastq_skipped":
        return "eligibility"
    if reason == "missing_credentials" or any(
        token in normalized
        for token in ("unauthorized", "forbidden", "401", "403", "credential")
    ):
        return "authentication"
    if any(
        token in normalized
        for token in (
            "timed out",
            "timeout",
            "connection reset",
            "connection refused",
            "temporary failure",
            "name or service not known",
            "could not resolve host",
            "failed to send http request",
            "error sending request",
            "dns error",
            "http 500",
            "http 502",
            "http 503",
            "http 504",
            "500 internal server error",
            "502 bad gateway",
            "503 service unavailable",
            "504 gateway timeout",
        )
    ):
        return "transport"
    if reason in {"no_fastq_sequence_files", "missing_fastq_metadata"}:
        return "portal"
    if (
        stage.startswith("seqspec")
        or stage.startswith("enumerate_")
        or reason
        in {
            "malformed_seqspec_yaml",
            "no_modalities",
            "no_fastq_files_for_modality",
        }
    ):
        return "specification"
    return "tool"


def access_class_for_records(records: list[SequenceFileRecord]) -> str:
    values = {value.controlled_access for value in records}
    if not values:
        return "unknown"
    if values == {False}:
        return "public"
    if values == {True}:
        return "controlled"
    return "mixed"


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_rows(csv_path: Path, jsonl_path: Path, rows: list[Any]) -> None:
    if not rows:
        csv_path.write_text("", encoding="utf-8")
        jsonl_path.write_text("", encoding="utf-8")
        return

    dictionaries = [asdict(row) for row in rows]
    fieldnames = list(dictionaries[0].keys())

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(dictionaries)

    with jsonl_path.open("w", encoding="utf-8") as handle:
        for row in dictionaries:
            handle.write(json.dumps(row, sort_keys=True))
            handle.write("\n")


def reconcile_outputs(
    output_root: Path,
    audit_context: AuditContext,
    configurations: list[ConfigurationRecord],
    runs: list[RunRecord],
    diagnostics: list[DiagnosticRecord],
    metrics: list[MetricRecord],
    failures: list[FailureRecord],
) -> dict[str, Any]:
    expected_configurations = {record.accession for record in configurations}
    observed_configurations = {
        row.configuration_accession for row in runs
    } | {row.configuration_accession for row in failures}

    run_keys = [(row.configuration_accession, row.modality) for row in runs]
    duplicate_run_keys = sorted(
        key for key, count in Counter(run_keys).items() if count > 1
    )
    outcome_keys = [
        (row.configuration_accession, row.modality) for row in [*runs, *failures]
    ]
    duplicate_outcome_keys = sorted(
        key for key, count in Counter(outcome_keys).items() if count > 1
    )
    catalog_reports = {Path(row.report_path).resolve() for row in runs}
    actual_reports = {
        path.resolve() for path in (output_root / "reports").glob("**/*.json")
    }

    report_parse_errors = []
    report_identity_errors = []
    report_sample_count_errors = []
    expected_assessments_by_report: Counter[str] = Counter()
    expected_metrics_by_report: Counter[str] = Counter()
    for run in runs:
        path = Path(run.report_path)
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as err:
            report_parse_errors.append({"report_path": str(path), "error": str(err)})
            continue

        summary = report.get("audit_summary", {})
        if (
            summary.get("study_run_id") != audit_context.run_id
            or summary.get("cache_key") != run.cache_key
            or summary.get("access_class") != run.access_class
            or summary.get("max_attempt_count") != run.max_attempt_count
            or summary.get("audit_protocol_id") != audit_context.audit_protocol_id
            or not report_cache_matches(report, run.cache_key)
        ):
            report_identity_errors.append(str(path))

        try:
            observed_sample_count = infer_sampled_record_count(report)
        except ValueError as err:
            report_sample_count_errors.append(
                {"report_path": str(path), "error": str(err)}
            )
        else:
            summary_sample_count = summary.get("sampled_record_count")
            if (
                observed_sample_count != run.sampled_record_count
                or summary_sample_count != run.sampled_record_count
            ):
                report_sample_count_errors.append(
                    {
                        "report_path": str(path),
                        "report_count": observed_sample_count,
                        "summary_count": summary_sample_count,
                        "run_count": run.sampled_record_count,
                    }
                )

        for result in report.get("results", []):
            report_key = str(path)
            expected_assessments_by_report[report_key] += len(
                result.get("assessment", [])
            )
            expected_metrics_by_report[report_key] += len(result.get("expected", []))
            expected_metrics_by_report[report_key] += len(result.get("observed", []))

    actual_assessments_by_report = Counter(row.report_path for row in diagnostics)
    actual_metrics_by_report = Counter(row.report_path for row in metrics)
    expected_assessment_count = sum(expected_assessments_by_report.values())
    expected_metric_count = sum(expected_metrics_by_report.values())
    known_cache_keys = {(row.report_path, row.cache_key) for row in runs}
    access_class_by_report = {row.report_path: row.access_class for row in runs}
    flattened_identity_errors = [
        {
            "table": table,
            "report_path": row.report_path,
            "cache_key": row.cache_key,
            "study_run_id": row.study_run_id,
        }
        for table, rows in (("diagnostics", diagnostics), ("metrics", metrics))
        for row in rows
        if row.study_run_id != audit_context.run_id
        or (row.report_path, row.cache_key) not in known_cache_keys
        or row.access_class != access_class_by_report.get(row.report_path)
    ]

    checks = {
        "configuration_outcomes_complete": (
            expected_configurations == observed_configurations
        ),
        "run_keys_unique": not duplicate_run_keys,
        "outcome_keys_unique": not duplicate_outcome_keys,
        "access_classes_valid": all(
            row.access_class in ACCESS_CLASSES for row in [*runs, *failures]
        ),
        "failure_categories_valid": all(
            row.failure_category in FAILURE_CATEGORIES for row in failures
        ),
        "unclassified_failures_absent": all(
            row.failure_category != "unclassified" for row in failures
        ),
        "attempt_counts_valid": (
            all(
                (row.run_status == "cached" and row.max_attempt_count >= 0)
                or 1 <= row.max_attempt_count <= audit_context.max_attempts
                for row in runs
            )
            and all(
                1 <= row.attempt_count <= audit_context.max_attempts
                for row in failures
            )
        ),
        "catalog_reports_exist": not (catalog_reports - actual_reports),
        "no_orphan_reports": not (actual_reports - catalog_reports),
        "reports_parse": not report_parse_errors,
        "report_identities_match": not report_identity_errors,
        "sampled_record_counts_match": not report_sample_count_errors,
        "completed_versions_recorded": all(
            row.raw_seqspec_version
            and row.normalized_seqspec_version == CURRENT_SEQSPEC_VERSION
            for row in runs
        ),
        "assessment_rows_reconcile": (
            expected_assessments_by_report == actual_assessments_by_report
        ),
        "metric_rows_reconcile": expected_metrics_by_report == actual_metrics_by_report,
        "flattened_identities_match": not flattened_identity_errors,
        "run_ids_match": (
            all(row.study_run_id == audit_context.run_id for row in runs)
            and all(row.study_run_id == audit_context.run_id for row in failures)
        ),
    }
    return {
        "audit_schema_version": AUDIT_SCHEMA_VERSION,
        "study_run_id": audit_context.run_id,
        "generated_at": utc_now(),
        "valid": all(checks.values()),
        "checks": checks,
        "counts": {
            "selected_configurations": len(expected_configurations),
            "observed_configurations": len(observed_configurations),
            "runs": len(runs),
            "failures": len(failures),
            "runs_by_access_class": dict(
                sorted(Counter(row.access_class for row in runs).items())
            ),
            "failures_by_access_class": dict(
                sorted(Counter(row.access_class for row in failures).items())
            ),
            "failures_by_category": dict(
                sorted(Counter(row.failure_category for row in failures).items())
            ),
            "runs_by_max_attempt_count": dict(
                sorted(Counter(row.max_attempt_count for row in runs).items())
            ),
            "failures_by_attempt_count": dict(
                sorted(Counter(row.attempt_count for row in failures).items())
            ),
            "diagnostics": len(diagnostics),
            "expected_assessments": expected_assessment_count,
            "metrics": len(metrics),
            "expected_metrics": expected_metric_count,
            "sampled_records": sum(row.sampled_record_count for row in runs),
            "catalog_reports": len(catalog_reports),
            "actual_reports": len(actual_reports),
        },
        "details": {
            "missing_configuration_outcomes": sorted(
                expected_configurations - observed_configurations
            ),
            "duplicate_run_keys": [list(key) for key in duplicate_run_keys],
            "duplicate_outcome_keys": [list(key) for key in duplicate_outcome_keys],
            "missing_reports": sorted(str(path) for path in catalog_reports - actual_reports),
            "orphan_reports": sorted(str(path) for path in actual_reports - catalog_reports),
            "report_parse_errors": report_parse_errors,
            "report_identity_errors": report_identity_errors,
            "report_sample_count_errors": report_sample_count_errors,
            "flattened_identity_errors": flattened_identity_errors,
            "assessment_count_mismatches": counter_differences(
                expected_assessments_by_report, actual_assessments_by_report
            ),
            "metric_count_mismatches": counter_differences(
                expected_metrics_by_report, actual_metrics_by_report
            ),
        },
    }


def counter_differences(expected: Counter[str], observed: Counter[str]) -> list[dict[str, Any]]:
    return [
        {
            "report_path": key,
            "expected": expected[key],
            "observed": observed[key],
        }
        for key in sorted(set(expected) | set(observed))
        if expected[key] != observed[key]
    ]


def write_lab_summary(
    path: Path, runs: list[RunRecord], failures: list[FailureRecord]
) -> None:
    summary: dict[str, dict[str, Any]] = {}
    for run in runs:
        row = summary.setdefault(
            run.lab or "(unknown)",
            {
                "lab": run.lab or "(unknown)",
                "run_count": 0,
                "failure_count": 0,
                "pass_count": 0,
                "warning_count": 0,
                "error_count": 0,
                "interpretation_count": 0,
            },
        )
        row["run_count"] += 1
        row["pass_count"] += run.pass_count
        row["warning_count"] += run.warning_count
        row["error_count"] += run.error_count
        row["interpretation_count"] += run.interpretation_count

    for failure in failures:
        row = summary.setdefault(
            failure.lab or "(unknown)",
            {
                "lab": failure.lab or "(unknown)",
                "run_count": 0,
                "failure_count": 0,
                "pass_count": 0,
                "warning_count": 0,
                "error_count": 0,
                "interpretation_count": 0,
            },
        )
        row["failure_count"] += 1

    rows = [summary[key] for key in sorted(summary)]
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def stderr_message(err: subprocess.CalledProcessError) -> str:
    stderr = (err.stderr or "").strip()
    stdout = (err.stdout or "").strip()
    if stderr:
        return stderr
    if stdout:
        return stdout
    return str(err)


def classify_seqspec_failure_reason(default_reason: str, message: str) -> str:
    normalized = message.strip().lower()
    if normalized == "could not read values.":
        return "malformed_seqspec_yaml"
    return default_reason


def slugify(value: str) -> str:
    cleaned = []
    for char in value:
        if char.isalnum() or char in {"-", "_", "."}:
            cleaned.append(char)
        else:
            cleaned.append("_")
    return "".join(cleaned).strip("_") or "unknown"


if __name__ == "__main__":
    raise SystemExit(main())
