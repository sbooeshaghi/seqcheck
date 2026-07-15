#!/usr/bin/env python3
"""Freeze the IGVF records consumed by the seqcheck paper audit."""

from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

try:
    import igvf_audit as audit
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import igvf_audit as audit
    from scripts import paper_runtime as runtime


FREEZER_VERSION = "0.1.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Freeze content-addressed IGVF records for the seqcheck audit."
    )
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--api-root", default=audit.DEFAULT_API_ROOT)
    parser.add_argument("--portal-root", default=audit.DEFAULT_PORTAL_ROOT)
    parser.add_argument("--status", default="released")
    parser.add_argument("--upload-status", default="validated")
    parser.add_argument(
        "--configuration-snapshot",
        type=Path,
        help="Existing IGVF configuration search response; requires --sequence-file-snapshot.",
    )
    parser.add_argument(
        "--sequence-file-snapshot",
        type=Path,
        help="Existing IGVF FASTQ search response; requires --configuration-snapshot.",
    )
    parser.add_argument(
        "--configuration-accession",
        action="append",
        default=[],
        help="Retain one configuration accession. Repeatable.",
    )
    parser.add_argument("--limit", type=int, help="Retain the first N configurations.")
    return parser.parse_args()


def main() -> int:
    try:
        path = run(parse_args())
    except (OSError, ValueError, RuntimeError) as err:
        print(f"freeze_igvf_audit_portal: {err}", file=sys.stderr)
        return 1
    print(path)
    return 0


def run(args: argparse.Namespace) -> Path:
    if (args.configuration_snapshot is None) != (args.sequence_file_snapshot is None):
        raise ValueError(
            "configuration and sequence-file snapshots must be supplied together"
        )
    if args.limit is not None and args.limit <= 0:
        raise ValueError("limit must be positive")

    output_root = args.output_root.resolve()
    if output_root.exists() and any(output_root.iterdir()):
        raise ValueError(f"output root is not empty: {output_root}")

    if args.configuration_snapshot is not None:
        configurations = load_configuration_snapshot(args.configuration_snapshot)
        sequence_files = load_sequence_file_snapshot(args.sequence_file_snapshot)
        source = {
            "mode": "snapshots",
            "configurations": runtime.file_identity(args.configuration_snapshot),
            "sequence_files": runtime.file_identity(args.sequence_file_snapshot),
        }
    else:
        configurations = audit.fetch_configuration_records(
            args.api_root,
            status=args.status,
            upload_status=args.upload_status,
        )
        sequence_files = audit.fetch_sequence_file_records(
            args.api_root,
            status=args.status,
            upload_status=args.upload_status,
        )
        source = {"mode": "live"}

    if args.configuration_accession:
        configurations = audit.select_requested_configurations(
            configurations, args.configuration_accession
        )
    configurations.sort(key=lambda record: record.accession)
    if args.limit is not None:
        configurations = configurations[: args.limit]
    if not configurations:
        raise ValueError("no configurations matched the requested frozen input")

    tools = {
        "freezer": runtime.script_identity(
            Path(__file__).resolve(), version=FREEZER_VERSION
        ),
        "audit_runtime": runtime.script_identity(
            Path(audit.__file__).resolve(), version=audit.AUDIT_SCHEMA_VERSION
        ),
    }
    manifest = build_portal_manifest(
        configurations=configurations,
        sequence_files=sequence_files,
        query={
            "api_root": args.api_root,
            "portal_root": args.portal_root,
            "status": args.status,
            "upload_status": args.upload_status,
        },
        tools=tools,
        source=source,
        retrieved_at=audit.utc_now(),
    )
    path = output_root / "manifests" / "portal.json"
    audit.write_json(path, manifest)
    audit.load_portal_manifest(path)
    return path


def build_portal_manifest(
    *,
    configurations: list[audit.ConfigurationRecord],
    sequence_files: dict[str, audit.SequenceFileRecord],
    query: dict[str, str],
    tools: dict[str, dict[str, Any]],
    source: dict[str, Any],
    retrieved_at: str,
) -> dict[str, Any]:
    configurations = sorted(configurations, key=lambda record: record.accession)
    linked_accessions = {
        accession
        for record in configurations
        for accession in audit.extract_accessions(record.seqspec_of)
    }
    linked_sequence_files = [
        sequence_files[accession]
        for accession in sorted(linked_accessions)
        if accession in sequence_files
    ]
    missing = sorted(linked_accessions - set(sequence_files))
    value = {
        "schema_version": audit.PORTAL_MANIFEST_SCHEMA_VERSION,
        "query": query,
        "configurations": [asdict(record) for record in configurations],
        "sequence_files": [asdict(record) for record in linked_sequence_files],
        "missing_linked_sequence_file_accessions": missing,
        "counts": {
            "configurations": len(configurations),
            "linked_sequence_files": len(linked_sequence_files),
            "missing_linked_sequence_files": len(missing),
        },
        "tools": tools,
        "source": source,
        "retrieved_at": retrieved_at,
    }
    return {
        **value,
        "portal_manifest_id": audit.sha256_json(audit.portal_manifest_stable(value))[
            :16
        ],
    }


def load_configuration_snapshot(path: Path) -> list[audit.ConfigurationRecord]:
    records = [
        audit.normalize_configuration_record(item)
        for item in snapshot_graph(path, "configuration")
    ]
    records = [record for record in records if record.accession]
    ensure_unique_accessions(records, "configuration")
    return records


def load_sequence_file_snapshot(
    path: Path,
) -> dict[str, audit.SequenceFileRecord]:
    records = [
        record
        for item in snapshot_graph(path, "sequence-file")
        if (record := audit.normalize_sequence_file_record(item)) is not None
    ]
    ensure_unique_accessions(records, "sequence-file")
    return {record.accession: record for record in records}


def snapshot_graph(path: Path, label: str) -> list[dict[str, Any]]:
    value = audit.load_json_object(path.resolve(), f"{label} snapshot")
    graph = value.get("@graph")
    if not isinstance(graph, list):
        raise ValueError(f"{label} snapshot @graph must be a list")
    records = [item for item in graph if isinstance(item, dict)]
    if not records:
        raise ValueError(f"{label} snapshot contains no records")
    return records


def ensure_unique_accessions(records: list[Any], label: str) -> None:
    accessions = [record.accession for record in records]
    if len(accessions) != len(set(accessions)):
        raise ValueError(f"{label} snapshot contains duplicate accessions")


if __name__ == "__main__":
    raise SystemExit(main())
