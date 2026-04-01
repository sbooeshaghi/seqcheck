#!/usr/bin/env python3
"""Audit IGVF seqspec-backed FASTQ datasets with seqcheck."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import shutil
import subprocess
import sys
import tomllib
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


USER_AGENT = "seqcheck-igvf-audit/0.1"
DEFAULT_API_ROOT = "https://api.data.igvf.org/"
DEFAULT_PORTAL_ROOT = DEFAULT_API_ROOT
CURRENT_SEQSPEC_VERSION = "0.4.0"


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
    expected_fastq_count: int
    supplied_fastq_count: int
    controlled_access: bool
    report_path: str
    run_status: str
    pass_count: int
    warning_count: int
    error_count: int
    interpretation_count: int


@dataclass(frozen=True)
class DiagnosticRecord:
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


@dataclass(frozen=True)
class FailureRecord:
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


@dataclass(frozen=True)
class ToolCommand:
    argv: list[str]
    env: dict[str, str]


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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "reports").mkdir(exist_ok=True)
    mpl_dir = output_root / ".mplconfig"
    mpl_dir.mkdir(exist_ok=True)

    seqspec_cmd = discover_seqspec_command(Path(__file__).resolve().parents[2])
    seqcheck_cmd = discover_seqcheck_command(Path(__file__).resolve().parents[1])

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

    if args.configuration_accession:
        allowed = {item.strip() for item in args.configuration_accession if item.strip()}
        configurations = [
            record for record in configurations if record.accession in allowed
        ]

    configurations.sort(key=lambda record: record.accession)
    if args.limit is not None:
        configurations = configurations[: args.limit]

    print(
        f"auditing {len(configurations)} configuration files against {len(sequence_files)} FASTQ records",
        file=sys.stderr,
    )
    total_configurations = len(configurations)

    runs: list[RunRecord] = []
    diagnostics: list[DiagnosticRecord] = []
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
            ): record.accession
            for record in configurations
        }

        completed_configurations = 0
        for future in as_completed(future_map):
            accession = future_map[future]
            completed_configurations += 1
            remaining_configurations = total_configurations - completed_configurations
            try:
                run_rows, diagnostic_rows, failure_rows = future.result()
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
                    FailureRecord(
                        configuration_accession=accession,
                        modality="",
                        lab="",
                        submitted_by="",
                        award_component="",
                        file_set_accession="",
                        assay_term="",
                        preferred_assay_titles="",
                        aliases="",
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
    write_lab_summary(output_root / "lab_summary.csv", runs, failures)
    return 0


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
) -> tuple[list[RunRecord], list[DiagnosticRecord], list[FailureRecord]]:
    runs: list[RunRecord] = []
    diagnostics: list[DiagnosticRecord] = []
    failures: list[FailureRecord] = []

    linked_fastqs = [
        sequence_files[accession]
        for accession in extract_accessions(record.seqspec_of)
        if accession in sequence_files
    ]

    if not linked_fastqs:
        failures.append(
            build_failure(
                record,
                modality="",
                raw_seqspec_version="",
                normalized_seqspec_version="",
                stage="sequence_files",
                reason="no_fastq_sequence_files",
                message="Configuration file has no linked FASTQ sequence files.",
            )
        )
        return runs, diagnostics, failures

    spec_url = absolute_url(portal_root, record.href)

    try:
        raw_seqspec_version = get_seqspec_file_version(
            seqspec_cmd,
            spec_url,
            seqspec_auth_profile,
            mpl_dir,
        )
    except subprocess.CalledProcessError as err:
        message = stderr_message(err)
        failures.append(
            build_failure(
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
            )
        )
        return runs, diagnostics, failures

    normalized_seqspec_version = normalize_seqspec_version(raw_seqspec_version)

    try:
        modalities = list_modalities(
            seqspec_cmd,
            spec_url,
            seqspec_auth_profile,
            mpl_dir,
        )
    except subprocess.CalledProcessError as err:
        message = stderr_message(err)
        failures.append(
            build_failure(
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
            )
        )
        return runs, diagnostics, failures

    if not modalities:
        failures.append(
            build_failure(
                record,
                modality="",
                raw_seqspec_version=raw_seqspec_version,
                normalized_seqspec_version=normalized_seqspec_version,
                stage="enumerate_modalities",
                reason="no_modalities",
                message="No modalities found in seqspec.",
            )
        )
        return runs, diagnostics, failures

    modality_count = len(modalities)

    for modality in modalities:
        report_path = output_root / "reports" / record.accession / f"{slugify(modality)}.json"
        report_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            expected_files = list_modality_files(
                seqspec_cmd,
                spec_url,
                modality,
                seqspec_auth_profile,
                mpl_dir,
            )
        except subprocess.CalledProcessError as err:
            message = stderr_message(err)
            failures.append(
                build_failure(
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
                )
            )
            continue

        fastq_expectations = [item for item in expected_files if is_fastq_expectation(item)]
        if not fastq_expectations:
            failures.append(
                build_failure(
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="enumerate_files",
                    reason="no_fastq_files_for_modality",
                    message=f"No FASTQ files found for modality '{modality}'.",
                )
            )
            continue

        modality_sequence_records: list[SequenceFileRecord] = []
        fastq_urls: list[str] = []
        controlled_needed = False

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
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="resolve_fastqs",
                    reason=skip_reason,
                    message=skip_message or "",
                )
            )
            continue

        if report_path.exists() and not args.force:
            try:
                report = json.loads(report_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                report = None
            else:
                report = annotate_report_summary(
                    report,
                    configuration_accession=record.accession,
                    modality=modality,
                    modality_count=modality_count,
                    expected_fastq_count=len(modality_sequence_records),
                    supplied_fastq_count=len(modality_sequence_records),
                    controlled_access=controlled_needed,
                )
                report_path.write_text(
                    json.dumps(report, indent=2, sort_keys=False) + "\n",
                    encoding="utf-8",
                )
                run_row, diagnostic_rows = flatten_report(
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
                )
                runs.append(run_row)
                diagnostics.extend(diagnostic_rows)
                continue

        try:
            run_seqcheck(
                seqcheck_cmd,
                spec_url,
                modality,
                args.n_reads,
                fastq_urls,
                report_path,
                seqcheck_auth_profile,
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
            )
            report_path.write_text(
                json.dumps(report, indent=2, sort_keys=False) + "\n",
                encoding="utf-8",
            )
        except subprocess.CalledProcessError as err:
            failures.append(
                build_failure(
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="seqcheck",
                    reason="seqcheck_error",
                    message=stderr_message(err),
                )
            )
            continue
        except Exception as err:
            failures.append(
                build_failure(
                    record,
                    modality=modality,
                    raw_seqspec_version=raw_seqspec_version,
                    normalized_seqspec_version=normalized_seqspec_version,
                    stage="seqcheck",
                    reason="report_parse_error",
                    message=str(err),
                )
            )
            continue

        run_row, diagnostic_rows = flatten_report(
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
        )
        runs.append(run_row)
        diagnostics.extend(diagnostic_rows)

    return runs, diagnostics, failures


def absolute_url(api_root: str, href_or_url: str) -> str:
    return urllib.parse.urljoin(api_root, href_or_url)


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
        accession = str(item.get("accession", "")).strip()
        if not accession:
            continue
        records[accession] = SequenceFileRecord(
            accession=accession,
            href=str(item.get("href", "")).strip(),
            controlled_access=bool(item.get("controlled_access", False)),
            read_names=[str(value) for value in item.get("read_names", [])],
        )
    return records


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
) -> str:
    env = os.environ.copy()
    env.update(command.env)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    argv = command.argv + ["version", spec_source]
    if auth_profile:
        argv.extend(["--auth-profile", auth_profile])
    result = subprocess.run(
        argv,
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    return parse_seqspec_version_output(result.stdout)


def list_modality_files(
    command: ToolCommand,
    spec_source: str,
    modality: str,
    auth_profile: str | None,
    mpl_dir: Path,
) -> list[dict[str, Any]]:
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
    result = subprocess.run(
        argv,
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    return json.loads(result.stdout)


def list_modalities(
    command: ToolCommand,
    spec_source: str,
    auth_profile: str | None,
    mpl_dir: Path,
) -> list[str]:
    env = os.environ.copy()
    env.update(command.env)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    argv = command.argv + ["info", "-k", "modalities", "-f", "json", spec_source]
    if auth_profile:
        argv.extend(["--auth-profile", auth_profile])
    result = subprocess.run(
        argv,
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    payload = json.loads(result.stdout)
    return [str(item) for item in payload]


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
) -> None:
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
    subprocess.run(
        argv,
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )


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
    if raw_version in {"0.0.0", "0.1.0", "0.1.1", "0.2.0", "0.3.0"}:
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
) -> tuple[RunRecord, list[DiagnosticRecord]]:
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

    run = RunRecord(
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
        expected_fastq_count=expected_fastq_count,
        supplied_fastq_count=supplied_fastq_count,
        controlled_access=controlled_access,
        report_path=str(report_path),
        run_status=run_status,
        pass_count=counts["pass"],
        warning_count=counts["warning"],
        error_count=counts["error"],
        interpretation_count=counts["interpretation"],
    )

    diagnostics = []
    for result in report.get("results", []):
        for assessment in result.get("assessment", []):
            diagnostics.append(
                DiagnosticRecord(
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
                )
            )
    return run, diagnostics


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
) -> dict[str, Any]:
    requested_reads = int(report.get("meta", {}).get("requested_reads", 0))
    report["audit_summary"] = {
        "configuration_accession": configuration_accession,
        "modality": modality,
        "modality_count": modality_count,
        "expected_fastq_count": expected_fastq_count,
        "supplied_fastq_count": supplied_fastq_count,
        "requested_reads_per_fastq": requested_reads,
        "requested_reads_total": (
            requested_reads * supplied_fastq_count if requested_reads > 0 else 0
        ),
        "controlled_access": controlled_access,
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
    record: ConfigurationRecord,
    modality: str,
    raw_seqspec_version: str,
    normalized_seqspec_version: str,
    stage: str,
    reason: str,
    message: str,
) -> FailureRecord:
    return FailureRecord(
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
    )


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
