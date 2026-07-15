#!/usr/bin/env python3
"""Build a reproducible, review-gated seqcheck paper cohort candidate pool."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "0.2.0"
CURRENT_SEQSPEC_VERSION = "0.5.0"
DEFAULT_PORTAL_ROOT = "https://api.data.igvf.org/"
USER_AGENT = "seqcheck-cohort-builder/0.2"
ONTOLOGY_TERM_PATTERN = re.compile(r"RGN:[A-Za-z0-9_]+:[A-Za-z0-9_]+")


@dataclass(frozen=True)
class FamilyRule:
    family_id: str
    label: str
    preferred_assay_titles: tuple[str, ...]
    assay_terms: tuple[str, ...]
    expected_modalities: tuple[str, ...]
    note: str


@dataclass(frozen=True)
class Configuration:
    accession: str
    href: str
    lab: str
    submitted_by: str
    file_set_accession: str
    assay_term: str
    preferred_assay_titles: tuple[str, ...]
    aliases: tuple[str, ...]
    seqspec_of: tuple[str, ...]
    status: str
    upload_status: str


@dataclass(frozen=True)
class SequenceFile:
    accession: str
    href: str
    controlled_access: bool | None
    read_names: tuple[str, ...]


@dataclass(frozen=True)
class Hydration:
    status: str = "not_run"
    message: str = ""
    raw_spec_sha256: str = ""
    download_attempts: int = 0
    raw_seqspec_version: str = ""
    normalized_spec_sha256: str = ""
    normalized_seqspec_version: str = ""
    modalities: tuple[str, ...] = ()
    ontology_terms: tuple[str, ...] = ()
    structure_sha256: str = ""
    expected_fastq_accessions: tuple[str, ...] = ()
    fastq_mapping_status: str = "not_run"
    spec_only_fastq_accessions: tuple[str, ...] = ()
    portal_only_fastq_accessions: tuple[str, ...] = ()
    structural_check_status: str = "not_run"
    structural_check_message: str = ""
    resource_check_status: str = "not_run"
    resource_check_message: str = ""
    resource_check_attempts: int = 0
    normalized_spec_path: str = ""


PROPOSAL_FIELDS = [
    "selection_id",
    "family_id",
    "family_label",
    "family_rule_titles",
    "family_rule_assay_terms",
    "family_expected_modalities",
    "family_rule_note",
    "proposal_rank",
    "family_rank",
    "review_selection_status",
    "review_exclusion_reason",
    "configuration_accession",
    "configuration_url",
    "lab",
    "submitted_by",
    "file_set_accession",
    "assay_term",
    "preferred_assay_titles",
    "aliases",
    "candidate_families",
    "ambiguous_family",
    "fastq_accessions",
    "fastq_urls",
    "fastq_count",
    "fastq_set_sha256",
    "metadata_eligibility",
    "metadata_reasons",
    "hydration_status",
    "hydration_message",
    "raw_spec_sha256",
    "download_attempts",
    "raw_seqspec_version",
    "normalized_spec_sha256",
    "normalized_seqspec_version",
    "modalities",
    "modality_match_status",
    "ontology_terms",
    "structure_sha256",
    "expected_fastq_accessions",
    "fastq_mapping_status",
    "spec_only_fastq_accessions",
    "portal_only_fastq_accessions",
    "structural_check_status",
    "structural_check_message",
    "resource_check_status",
    "resource_check_message",
    "resource_check_attempts",
    "deduplication_key",
    "normalized_spec_path",
]

CANDIDATE_FIELDS = PROPOSAL_FIELDS

REVIEW_FIELDS = PROPOSAL_FIELDS + [
    "reviewer_1",
    "reviewer_1_decision",
    "reviewer_1_rationale",
    "reviewer_1_protocol_url",
    "reviewer_1_date",
    "reviewer_2",
    "reviewer_2_decision",
    "reviewer_2_rationale",
    "reviewer_2_protocol_url",
    "reviewer_2_date",
    "adjudication_decision",
    "adjudication_rationale",
    "final_family",
    "split",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build deterministic seqcheck paper cohort candidates. Automated family "
            "matches are proposals; this command never approves or freezes a cohort."
        )
    )
    parser.add_argument("--configurations", required=True, type=Path)
    parser.add_argument("--sequence-files", required=True, type=Path)
    parser.add_argument(
        "--rules",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "docs"
        / "cohort_family_rules.json",
    )
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--portal-root", default=DEFAULT_PORTAL_ROOT)
    parser.add_argument("--candidate-count", type=int)
    parser.add_argument(
        "--proposal-count",
        type=int,
        help=(
            "Metadata proposals to hydrate per family before modality screening "
            "(default: rules registry)."
        ),
    )
    parser.add_argument("--seed", type=int)
    parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help="Concurrent seqspec hydration workers (default: 2).",
    )
    parser.add_argument(
        "--structural-check-timeout-seconds",
        type=int,
        default=120,
        help="Maximum time for each local structural seqspec check (default: 120).",
    )
    parser.add_argument(
        "--check-timeout-seconds",
        type=int,
        default=120,
        help="Maximum time for each resource-aware seqspec check (default: 120).",
    )
    parser.add_argument(
        "--network-attempts",
        type=int,
        default=3,
        help="Attempts for transient downloads and resource checks (default: 3).",
    )
    parser.add_argument(
        "--retry-backoff-seconds",
        type=float,
        default=2.0,
        help="Fixed delay between transient network attempts (default: 2).",
    )
    parser.add_argument(
        "--hydrate",
        action="store_true",
        help="Download and normalize proposed seqspec files before modality screening.",
    )
    parser.add_argument(
        "--seqspec-bin",
        help="seqspec executable used with --hydrate (default: discover local binary).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        return run(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"build_cohort_candidates: {error}", file=sys.stderr)
        return 1


def run(args: argparse.Namespace) -> int:
    rules_payload = load_json(args.rules)
    rules = parse_family_rules(rules_payload)
    proposal_count = (
        args.proposal_count
        if args.proposal_count is not None
        else int(rules_payload["proposal_count_per_family"])
    )
    candidate_count = (
        args.candidate_count
        if args.candidate_count is not None
        else int(rules_payload["candidate_count_per_family"])
    )
    seed = args.seed if args.seed is not None else int(rules_payload["selection_seed"])
    minimum_final_count = int(rules_payload.get("minimum_final_count_per_family", 5))
    if proposal_count <= 0 or candidate_count <= 0:
        raise ValueError("proposal and candidate counts must be positive")
    if proposal_count < candidate_count:
        raise ValueError("proposal count must be no less than candidate count")
    if minimum_final_count <= 0 or minimum_final_count > candidate_count:
        raise ValueError(
            "minimum final count must be positive and no greater than candidate count"
        )
    if args.workers <= 0:
        raise ValueError("workers must be positive")
    if args.structural_check_timeout_seconds <= 0 or args.check_timeout_seconds <= 0:
        raise ValueError("structural and resource check timeouts must be positive")
    if args.network_attempts <= 0:
        raise ValueError("network attempts must be positive")
    if args.retry_backoff_seconds < 0:
        raise ValueError("retry backoff must be non-negative")

    configuration_payload = load_json(args.configurations)
    sequence_payload = load_json(args.sequence_files)
    configurations = parse_configurations(configuration_payload)
    sequence_files = parse_sequence_files(sequence_payload)
    if not configurations:
        raise ValueError("configuration snapshot contains no records")
    if not sequence_files:
        raise ValueError("sequence-file snapshot contains no records")

    inventory = build_inventory(configurations, sequence_files, rules, args.portal_root)
    proposals = select_proposals(inventory, rules, proposal_count, seed)
    if not proposals:
        raise ValueError("no eligible cohort proposals matched the family rules")

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    for directory in (
        "manifests",
        "tables",
        "validation",
        "specs/raw",
        "specs/normalized",
    ):
        (output_root / directory).mkdir(parents=True, exist_ok=True)

    hydration_by_accession: dict[str, Hydration] = {}
    seqspec_identity: dict[str, Any] = {}
    if args.hydrate:
        seqspec_command = discover_seqspec_command(args.seqspec_bin)
        seqspec_identity = build_seqspec_identity(seqspec_command)
        rows_by_accession = {row["configuration_accession"]: row for row in proposals}
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    hydrate_configuration,
                    row,
                    output_root,
                    seqspec_command,
                    args.structural_check_timeout_seconds,
                    args.check_timeout_seconds,
                    args.network_attempts,
                    args.retry_backoff_seconds,
                ): accession
                for accession, row in sorted(rows_by_accession.items())
            }
            for future in as_completed(futures):
                accession = futures[future]
                hydration = future.result()
                hydration_by_accession[accession] = hydration
                print(
                    f"hydrated {accession}: {hydration.status}, "
                    f"structural={hydration.structural_check_status}, "
                    f"resources={hydration.resource_check_status}, "
                    f"attempts={hydration.resource_check_attempts}, "
                    f"fastqs={hydration.fastq_mapping_status}",
                    file=sys.stderr,
                )

    proposals = [
        enrich_candidate(
            row, hydration_by_accession.get(row["configuration_accession"])
        )
        for row in proposals
    ]
    proposals, selected = select_review_candidates(
        proposals, rules, candidate_count, seed
    )
    stable_selection = selection_identity_payload(
        proposals,
        selected,
        rules_payload,
        proposal_count,
        candidate_count,
        seed,
        file_sha256(args.configurations),
        file_sha256(args.sequence_files),
        functional_seqspec_identity(seqspec_identity),
    )
    selection_id = sha256_json(stable_selection)[:16]
    for row in proposals:
        row["selection_id"] = selection_id

    inventory_path = output_root / "tables" / "candidate_inventory.csv"
    proposals_path = output_root / "tables" / "cohort_proposals.csv"
    candidates_path = output_root / "tables" / "cohort_candidates.csv"
    review_template_path = output_root / "tables" / "cohort_reviews.template.csv"
    review_path = output_root / "tables" / "cohort_reviews.csv"
    write_csv(inventory_path, inventory, inventory_fields(inventory))
    write_csv(proposals_path, proposals, PROPOSAL_FIELDS)
    write_csv(candidates_path, selected, CANDIDATE_FIELDS)
    review_rows = [
        {**row, **{field: "" for field in REVIEW_FIELDS if field not in row}}
        for row in selected
    ]
    write_csv(review_template_path, review_rows, REVIEW_FIELDS)
    review_file_status = "preserved"
    review_selection_ids: list[str] = []
    if not review_path.exists():
        write_csv(review_path, review_rows, REVIEW_FIELDS)
        review_file_status = "created"
        review_selection_ids = [selection_id] if review_rows else []
    else:
        review_selection_ids = read_review_selection_ids(review_path)
        expected_review_ids = [selection_id] if review_rows else []
        if review_selection_ids != expected_review_ids:
            review_file_status = "stale_preserved"

    validation = build_validation(
        inventory,
        proposals,
        selected,
        rules,
        proposal_count,
        candidate_count,
        minimum_final_count,
        selection_id,
        review_file_status,
    )
    manifest = {
        "cohort_candidate_schema_version": SCHEMA_VERSION,
        "selection_id": selection_id,
        "generated_at": utc_now(),
        "frozen": False,
        "split_assigned": False,
        "selection_uses_seqcheck_results": False,
        "selection_gates": [
            "public metadata eligibility",
            "successful seqspec normalization",
            "expected seqspec modality",
        ],
        "selection_stratification": [
            "submitting laboratory",
            "portal-linked FASTQ count within laboratory",
        ],
        "selection_excludes": [
            "structural check outcome",
            "resource check outcome",
            "FASTQ reconciliation outcome",
            "seqcheck result",
        ],
        "rules": {
            "path": str(args.rules.resolve()),
            "sha256": file_sha256(args.rules),
            "schema_version": rules_payload.get("schema_version", ""),
            "proposal_count_per_family": proposal_count,
            "candidate_count_per_family": candidate_count,
            "minimum_final_count_per_family": minimum_final_count,
            "selection_seed": seed,
        },
        "inputs": {
            "configurations": snapshot_identity(
                args.configurations, configuration_payload, len(configurations)
            ),
            "sequence_files": snapshot_identity(
                args.sequence_files, sequence_payload, len(sequence_files)
            ),
            "portal_root": args.portal_root,
        },
        "hydration": {
            "requested": bool(args.hydrate),
            "seqspec": seqspec_identity,
            "workers": args.workers,
            "structural_check_timeout_seconds": args.structural_check_timeout_seconds,
            "resource_check_timeout_seconds": args.check_timeout_seconds,
            "network_attempts": args.network_attempts,
            "retry_backoff_seconds": args.retry_backoff_seconds,
            "counts": dict(
                sorted(
                    Counter(
                        value.status for value in hydration_by_accession.values()
                    ).items()
                )
            ),
        },
        "outputs": {
            "candidate_inventory": str(inventory_path),
            "cohort_proposals": str(proposals_path),
            "cohort_candidates": str(candidates_path),
            "cohort_review_template": str(review_template_path),
            "cohort_reviews": str(review_path),
            "review_file_status": review_file_status,
            "review_selection_ids": review_selection_ids,
        },
    }
    write_json(output_root / "manifests" / "cohort_candidates.json", manifest)
    write_json(output_root / "validation" / "cohort_selection.json", validation)

    action = "hydrated" if args.hydrate else "prepared"
    print(
        f"{action} {len(proposals)} proposals and selected {len(selected)} "
        f"review candidates across {len(rules)} families "
        f"(selection_id={selection_id}, frozen=false)"
    )
    return 0


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def graph_records(payload: dict[str, Any]) -> list[dict[str, Any]]:
    graph = payload.get("@graph", [])
    if not isinstance(graph, list):
        raise ValueError("portal snapshot @graph must be a list")
    return [item for item in graph if isinstance(item, dict)]


def parse_family_rules(payload: dict[str, Any]) -> list[FamilyRule]:
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"unsupported family-rule schema: {payload.get('schema_version', '')}"
        )
    rules = []
    seen: set[str] = set()
    for item in payload.get("families", []):
        family_id = str(item.get("id", "")).strip()
        titles = tuple(
            sorted(
                {
                    str(value).strip()
                    for value in item.get("preferred_assay_titles", [])
                    if str(value).strip()
                }
            )
        )
        assay_terms = tuple(
            sorted({str(value).strip() for value in item.get("assay_terms", [])})
        )
        expected_modalities = tuple(
            sorted(
                {
                    str(value).strip()
                    for value in item.get("expected_modalities", [])
                    if str(value).strip()
                }
            )
        )
        if not family_id or family_id in seen or not titles:
            raise ValueError(
                "each family needs a unique id and at least one exact title"
            )
        seen.add(family_id)
        rules.append(
            FamilyRule(
                family_id=family_id,
                label=str(item.get("label", "")).strip(),
                preferred_assay_titles=titles,
                assay_terms=assay_terms,
                expected_modalities=expected_modalities,
                note=str(item.get("note", "")).strip(),
            )
        )
    if not rules:
        raise ValueError("family-rule registry contains no families")
    return rules


def parse_configurations(payload: dict[str, Any]) -> list[Configuration]:
    records = []
    for item in graph_records(payload):
        accession = str(item.get("accession", "")).strip()
        if not accession:
            continue
        records.append(
            Configuration(
                accession=accession,
                href=str(item.get("href", "")).strip(),
                lab=nested_string(item, "lab", "title"),
                submitted_by=nested_string(item, "submitted_by", "title"),
                file_set_accession=nested_string(item, "file_set", "accession"),
                assay_term=nested_string(item, "file_set", "assay_term", "term_name"),
                preferred_assay_titles=string_tuple(
                    item.get("preferred_assay_titles", [])
                ),
                aliases=string_tuple(item.get("aliases", [])),
                seqspec_of=string_tuple(item.get("seqspec_of", [])),
                status=str(item.get("status", "")).strip(),
                upload_status=str(item.get("upload_status", "")).strip(),
            )
        )
    return sorted(records, key=lambda record: record.accession)


def parse_sequence_files(payload: dict[str, Any]) -> dict[str, SequenceFile]:
    records = {}
    for item in graph_records(payload):
        accession = str(item.get("accession", "")).strip()
        if not accession:
            continue
        controlled_access = item.get("controlled_access")
        if not isinstance(controlled_access, bool):
            controlled_access = None
        records[accession] = SequenceFile(
            accession=accession,
            href=str(item.get("href", "")).strip(),
            controlled_access=controlled_access,
            read_names=string_tuple(item.get("read_names", [])),
        )
    return records


def nested_string(item: dict[str, Any], *path: str) -> str:
    value: Any = item
    for key in path:
        if not isinstance(value, dict):
            return ""
        value = value.get(key)
    return "" if value is None else str(value).strip()


def string_tuple(value: Any) -> tuple[str, ...]:
    if not isinstance(value, list):
        return ()
    return tuple(str(item).strip() for item in value if str(item).strip())


def extract_accession(path: str) -> str:
    pieces = [piece for piece in path.strip("/").split("/") if piece]
    return pieces[-1] if pieces else ""


def match_families(
    configuration: Configuration, rules: list[FamilyRule]
) -> list[tuple[FamilyRule, tuple[str, ...]]]:
    observed = set(configuration.preferred_assay_titles)
    matches = []
    for rule in rules:
        titles = tuple(sorted(observed.intersection(rule.preferred_assay_titles)))
        assay_matches = not rule.assay_terms or configuration.assay_term in set(
            rule.assay_terms
        )
        if titles and assay_matches:
            matches.append((rule, titles))
    return matches


def metadata_eligibility(
    configuration: Configuration, sequence_files: dict[str, SequenceFile]
) -> tuple[str, list[str], list[SequenceFile]]:
    reasons = []
    if configuration.status != "released":
        reasons.append("configuration_not_released")
    if configuration.upload_status != "validated":
        reasons.append("configuration_not_validated")
    if not configuration.href:
        reasons.append("missing_configuration_href")
    if not configuration.lab:
        reasons.append("missing_lab")
    if not configuration.file_set_accession:
        reasons.append("missing_file_set_accession")
    accessions = tuple(extract_accession(value) for value in configuration.seqspec_of)
    accessions = tuple(value for value in accessions if value)
    if not accessions:
        reasons.append("missing_linked_fastqs")
    missing = sorted(value for value in accessions if value not in sequence_files)
    if missing:
        reasons.append("missing_fastq_metadata:" + ";".join(missing))
    linked = [sequence_files[value] for value in accessions if value in sequence_files]
    if any(item.controlled_access is True for item in linked):
        reasons.append("controlled_fastq")
    if any(item.controlled_access is None for item in linked):
        reasons.append("unknown_fastq_access")
    if any(not item.href for item in linked):
        reasons.append("missing_fastq_href")
    return ("eligible" if not reasons else "ineligible", reasons, linked)


def build_inventory(
    configurations: list[Configuration],
    sequence_files: dict[str, SequenceFile],
    rules: list[FamilyRule],
    portal_root: str,
) -> list[dict[str, Any]]:
    inventory = []
    for configuration in configurations:
        matches = match_families(configuration, rules)
        if not matches:
            continue
        eligibility, reasons, linked = metadata_eligibility(
            configuration, sequence_files
        )
        candidate_families = tuple(rule.family_id for rule, _ in matches)
        fastq_accessions = tuple(item.accession for item in linked)
        fastq_urls = tuple(
            absolute_url(portal_root, item.href) if item.href else "" for item in linked
        )
        for rule, matched_titles in matches:
            inventory.append(
                {
                    "family_id": rule.family_id,
                    "family_label": rule.label,
                    "family_rule_titles": ";".join(matched_titles),
                    "family_rule_assay_terms": ";".join(
                        value or "<missing>" for value in rule.assay_terms
                    ),
                    "family_expected_modalities": ";".join(rule.expected_modalities),
                    "family_rule_note": rule.note,
                    "configuration_accession": configuration.accession,
                    "configuration_url": absolute_url(portal_root, configuration.href),
                    "lab": configuration.lab,
                    "submitted_by": configuration.submitted_by,
                    "file_set_accession": configuration.file_set_accession,
                    "assay_term": configuration.assay_term,
                    "preferred_assay_titles": ";".join(
                        configuration.preferred_assay_titles
                    ),
                    "aliases": ";".join(configuration.aliases),
                    "candidate_families": ";".join(candidate_families),
                    "ambiguous_family": len(candidate_families) > 1,
                    "fastq_accessions": ";".join(fastq_accessions),
                    "fastq_urls": ";".join(fastq_urls),
                    "fastq_count": len(fastq_accessions),
                    "fastq_set_sha256": sha256_json(sorted(fastq_accessions)),
                    "metadata_eligibility": eligibility,
                    "metadata_reasons": ";".join(reasons),
                }
            )
    return sorted(
        inventory,
        key=lambda row: (row["family_id"], row["configuration_accession"]),
    )


def select_proposals(
    inventory: list[dict[str, Any]],
    rules: list[FamilyRule],
    proposal_count: int,
    seed: int,
) -> list[dict[str, Any]]:
    proposals = []
    for rule in rules:
        eligible = [
            row
            for row in inventory
            if row["family_id"] == rule.family_id
            and row["metadata_eligibility"] == "eligible"
        ]
        family_rows = select_lab_balanced(
            eligible, rule.family_id, proposal_count, seed
        )
        for rank, row in enumerate(family_rows, start=1):
            proposals.append(
                {
                    **row,
                    "proposal_rank": rank,
                    "family_rank": "",
                    "review_selection_status": "pending_hydration",
                    "review_exclusion_reason": "",
                }
            )
    return proposals


def select_review_candidates(
    proposals: list[dict[str, Any]],
    rules: list[FamilyRule],
    candidate_count: int,
    seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    selected_ranks: dict[tuple[str, str], int] = {}
    for rule in rules:
        qualified = [
            row
            for row in proposals
            if row["family_id"] == rule.family_id and not review_exclusion_reason(row)
        ]
        family_rows = select_lab_balanced(
            qualified, rule.family_id, candidate_count, seed
        )
        for rank, row in enumerate(family_rows, start=1):
            selected_ranks[(row["family_id"], row["configuration_accession"])] = rank

    annotated = []
    for row in proposals:
        key = (row["family_id"], row["configuration_accession"])
        exclusion_reason = review_exclusion_reason(row)
        if key in selected_ranks:
            status = "selected"
            family_rank: int | str = selected_ranks[key]
        elif exclusion_reason:
            status = "excluded"
            family_rank = ""
        else:
            status = "qualified_not_selected"
            exclusion_reason = "candidate_target_reached"
            family_rank = ""
        annotated.append(
            {
                **row,
                "family_rank": family_rank,
                "review_selection_status": status,
                "review_exclusion_reason": exclusion_reason,
            }
        )

    family_order = {rule.family_id: index for index, rule in enumerate(rules)}
    selected = sorted(
        (row for row in annotated if row["review_selection_status"] == "selected"),
        key=lambda row: (family_order[row["family_id"]], row["family_rank"]),
    )
    return annotated, selected


def review_exclusion_reason(row: dict[str, Any]) -> str:
    match_status = row.get("modality_match_status", "not_run")
    if match_status in {"matched", "not_required"}:
        return ""
    if row.get("hydration_status") != "normalized":
        return f"normalization_{row.get('hydration_status', 'not_run')}"
    if match_status == "mismatch":
        observed = row.get("modalities", "") or "<missing>"
        expected = row.get("family_expected_modalities", "") or "<unspecified>"
        return f"modality_mismatch:observed={observed};expected={expected}"
    return "modality_not_evaluated"


def select_lab_balanced(
    rows: list[dict[str, Any]],
    family_id: str,
    limit: int,
    seed: int,
) -> list[dict[str, Any]]:
    rows_by_lab_and_layout: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in rows:
        layout = str(row.get("fastq_count", ""))
        rows_by_lab_and_layout[row["lab"]][layout].append(row)

    by_lab: dict[str, list[dict[str, Any]]] = {}
    for lab, rows_by_layout in rows_by_lab_and_layout.items():
        for layout, layout_rows in rows_by_layout.items():
            layout_rows.sort(
                key=lambda row: stable_order_key(
                    seed,
                    family_id,
                    lab,
                    layout,
                    row["configuration_accession"],
                )
            )
        layouts = sorted(
            rows_by_layout,
            key=lambda layout: stable_order_key(seed, family_id, lab, layout),
        )
        by_lab[lab] = round_robin(
            rows_by_layout,
            layouts,
            sum(len(layout_rows) for layout_rows in rows_by_layout.values()),
        )
    labs = sorted(
        by_lab,
        key=lambda lab: stable_order_key(seed, family_id, lab),
    )
    return round_robin(by_lab, labs, limit)


def round_robin(
    grouped: dict[str, list[dict[str, Any]]],
    groups: list[str],
    limit: int,
) -> list[dict[str, Any]]:
    selected = []
    offset = 0
    while len(selected) < limit:
        added = False
        for group in groups:
            rows = grouped[group]
            if offset < len(rows):
                selected.append(rows[offset])
                added = True
                if len(selected) == limit:
                    return selected
        if not added:
            break
        offset += 1
    return selected


def stable_order_key(seed: int, *values: str) -> str:
    return sha256_json([seed, *values])


def discover_seqspec_command(requested: str | None) -> list[str]:
    if requested:
        executable = Path(requested).expanduser().resolve()
        if not executable.is_file():
            raise RuntimeError(f"seqspec executable does not exist: {executable}")
        return [str(executable)]
    env_bin = os.environ.get("SEQSPEC_BIN")
    if env_bin:
        return discover_seqspec_command(env_bin)
    sibling = (
        Path(__file__).resolve().parents[2] / "seqspec" / "target" / "debug" / "seqspec"
    )
    if sibling.is_file():
        return [str(sibling)]
    path_bin = shutil.which("seqspec")
    if path_bin:
        return [path_bin]
    raise RuntimeError("could not discover seqspec; pass --seqspec-bin")


def build_seqspec_identity(command: list[str]) -> dict[str, Any]:
    executable = Path(command[0]).resolve()
    return {
        "command": command,
        "version": run_command(command + ["--version"]).strip(),
        "executable_sha256": file_sha256(executable) if executable.is_file() else "",
    }


def functional_seqspec_identity(identity: dict[str, Any]) -> dict[str, Any]:
    return {
        "version": identity.get("version", ""),
        "executable_sha256": identity.get("executable_sha256", ""),
    }


def hydrate_configuration(
    row: dict[str, Any],
    output_root: Path,
    seqspec_command: list[str],
    structural_check_timeout_seconds: int,
    resource_check_timeout_seconds: int,
    network_attempts: int,
    retry_backoff_seconds: float,
) -> Hydration:
    accession = row["configuration_accession"]
    try:
        raw_path = output_root / "specs" / "raw" / f"{accession}.yaml"
        normalized_path = output_root / "specs" / "normalized" / f"{accession}.yaml"
        if raw_path.exists():
            source_bytes = raw_path.read_bytes()
            download_attempts = 0
        else:
            try:
                raw_bytes, download_attempts = download_bytes(
                    row["configuration_url"],
                    network_attempts,
                    retry_backoff_seconds,
                )
            except OSError as error:
                return Hydration(
                    status="failed",
                    message=str(error),
                    download_attempts=network_attempts,
                )
            source_bytes = (
                gzip.decompress(raw_bytes)
                if raw_bytes.startswith(b"\x1f\x8b")
                else raw_bytes
            )
            raw_path.write_bytes(source_bytes)
        raw_version = parse_seqspec_version(
            run_command(seqspec_command + ["version", str(raw_path)])
        )
        run_command(
            seqspec_command + ["upgrade", str(raw_path), "-o", str(normalized_path)]
        )
        version_output = run_command(
            seqspec_command + ["version", str(normalized_path)]
        )
        version = parse_seqspec_version(version_output)
        if version != CURRENT_SEQSPEC_VERSION:
            raise RuntimeError(
                f"normalization produced seqspec {version or 'unknown'}, expected {CURRENT_SEQSPEC_VERSION}"
            )
        modalities = tuple(
            str(value)
            for value in json.loads(
                run_command(
                    seqspec_command
                    + ["info", "-k", "modalities", "-f", "json", str(normalized_path)]
                )
            )
        )
        library_spec = json.loads(
            run_command(
                seqspec_command
                + ["info", "-k", "library_spec", "-f", "json", str(normalized_path)]
            )
        )
        sequence_spec = json.loads(
            run_command(
                seqspec_command
                + ["info", "-k", "sequence_spec", "-f", "json", str(normalized_path)]
            )
        )
        (
            structural_status,
            structural_message,
            resource_status,
            resource_message,
            resource_attempts,
        ) = run_seqspec_checks(
            seqspec_command,
            normalized_path,
            structural_check_timeout_seconds,
            resource_check_timeout_seconds,
            network_attempts,
            retry_backoff_seconds,
        )
        structure = normalized_structure(library_spec, sequence_spec)
        ontology_terms = tuple(sorted(extract_ontology_terms(library_spec)))
        expected_fastqs = tuple(
            sorted(extract_expected_fastq_accessions(sequence_spec))
        )
        linked_fastqs = tuple(
            sorted(value for value in row["fastq_accessions"].split(";") if value)
        )
        expected_set = set(expected_fastqs)
        linked_set = set(linked_fastqs)
        spec_only = tuple(sorted(expected_set - linked_set))
        portal_only = tuple(sorted(linked_set - expected_set))
        if not expected_fastqs:
            fastq_mapping_status = "missing_expected_fastqs"
        elif not spec_only and not portal_only:
            fastq_mapping_status = "matched"
        elif portal_only and not spec_only:
            fastq_mapping_status = "portal_extra_fastqs"
        elif spec_only and not portal_only:
            fastq_mapping_status = "spec_extra_fastqs"
        else:
            fastq_mapping_status = "mismatch"
        return Hydration(
            status="normalized",
            raw_spec_sha256=hashlib.sha256(source_bytes).hexdigest(),
            download_attempts=download_attempts,
            raw_seqspec_version=raw_version,
            normalized_spec_sha256=file_sha256(normalized_path),
            normalized_seqspec_version=version,
            modalities=modalities,
            ontology_terms=ontology_terms,
            structure_sha256=sha256_json(structure),
            expected_fastq_accessions=expected_fastqs,
            fastq_mapping_status=fastq_mapping_status,
            spec_only_fastq_accessions=spec_only,
            portal_only_fastq_accessions=portal_only,
            structural_check_status=structural_status,
            structural_check_message=structural_message,
            resource_check_status=resource_status,
            resource_check_message=resource_message,
            resource_check_attempts=resource_attempts,
            normalized_spec_path=str(normalized_path),
        )
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        return Hydration(status="failed", message=str(error))


def download_bytes(
    url: str, attempts: int, retry_backoff_seconds: float
) -> tuple[bytes, int]:
    last_error: OSError | None = None
    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read(), attempt
        except OSError as error:
            last_error = error
            if attempt < attempts:
                time.sleep(retry_backoff_seconds)
    assert last_error is not None
    raise last_error


def run_command(argv: list[str]) -> str:
    result = subprocess.run(argv, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(f"{' '.join(argv)} failed: {message}")
    return result.stdout


def run_check_command(
    argv: list[str],
    timeout_seconds: int,
    attempts: int,
    retry_backoff_seconds: float,
) -> tuple[str, str, int]:
    for attempt in range(1, attempts + 1):
        try:
            result = subprocess.run(
                argv,
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
            )
            message = (result.stderr.strip() or result.stdout.strip()).replace(
                "\n", " | "
            )
            status = classify_check_status(result.returncode, message)
        except subprocess.TimeoutExpired:
            status = "timeout"
            message = f"seqspec check exceeded {timeout_seconds} seconds"
        if status not in {"unavailable", "timeout"} or attempt == attempts:
            return status, message, attempt
        time.sleep(retry_backoff_seconds)
    raise AssertionError("unreachable retry loop")


def classify_check_status(returncode: int, message: str) -> str:
    if returncode == 0:
        return "passed"
    lowered = message.lower()
    if any(
        marker in lowered
        for marker in (
            "failed to send http request",
            "could not resolve host",
            "connection refused",
            "timed out",
        )
    ):
        return "unavailable"
    return "failed"


def run_seqspec_checks(
    seqspec_command: list[str],
    normalized_path: Path,
    structural_timeout_seconds: int,
    resource_timeout_seconds: int,
    network_attempts: int,
    retry_backoff_seconds: float,
) -> tuple[str, str, str, str, int]:
    structural_status, structural_message, _ = run_check_command(
        seqspec_command + ["check", "--skip", "external", str(normalized_path)],
        structural_timeout_seconds,
        1,
        0,
    )
    if structural_status != "passed":
        return structural_status, structural_message, "not_run", "", 0

    resource_status, resource_message, resource_attempts = run_check_command(
        seqspec_command + ["check", str(normalized_path)],
        resource_timeout_seconds,
        network_attempts,
        retry_backoff_seconds,
    )
    return (
        structural_status,
        structural_message,
        resource_status,
        resource_message,
        resource_attempts,
    )


def parse_seqspec_version(output: str) -> str:
    for line in output.splitlines():
        if line.startswith("seqspec file version:"):
            return line.split(":", 1)[1].strip()
    return ""


def normalized_structure(library_spec: Any, sequence_spec: Any) -> dict[str, Any]:
    reads = []
    if isinstance(sequence_spec, list):
        for value in sequence_spec:
            if isinstance(value, dict):
                reads.append(
                    {key: item for key, item in value.items() if key != "files"}
                )
            else:
                reads.append(value)
    return {"library_spec": library_spec, "sequence_spec": reads}


def extract_expected_fastq_accessions(sequence_spec: Any) -> set[str]:
    accessions = set()
    if not isinstance(sequence_spec, list):
        return accessions
    for read in sequence_spec:
        if not isinstance(read, dict):
            continue
        files = read.get("files", [])
        if not isinstance(files, list):
            continue
        for item in files:
            if not isinstance(item, dict) or not is_fastq_file(item):
                continue
            file_id = str(item.get("file_id", "")).strip()
            url_accession = sequence_accession_from_url(str(item.get("url", "")))
            if url_accession:
                accessions.add(url_accession)
                continue
            if file_id:
                accessions.add(normalize_fastq_accession(file_id))
                continue
            filename = Path(str(item.get("filename", ""))).name
            accessions.add(normalize_fastq_accession(filename))
    return {value for value in accessions if value}


def sequence_accession_from_url(value: str) -> str:
    pieces = [piece for piece in urllib.parse.urlparse(value).path.split("/") if piece]
    for index, piece in enumerate(pieces[:-1]):
        if piece == "sequence-files":
            return normalize_fastq_accession(pieces[index + 1])
    return ""


def normalize_fastq_accession(value: str) -> str:
    normalized = extract_accession(value)
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if normalized.lower().endswith(suffix):
            return normalized[: -len(suffix)]
    return normalized


def is_fastq_file(item: dict[str, Any]) -> bool:
    filetype = str(item.get("filetype", "")).strip().lower().lstrip(".")
    if filetype in {"fastq", "fq", "fastq.gz", "fq.gz"}:
        return True
    filename = str(item.get("filename", "")).strip().lower()
    return filename.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))


def extract_ontology_terms(value: Any) -> set[str]:
    terms: set[str] = set()
    if isinstance(value, dict):
        for child in value.values():
            terms.update(extract_ontology_terms(child))
    elif isinstance(value, list):
        for child in value:
            terms.update(extract_ontology_terms(child))
    elif isinstance(value, str):
        terms.update(ONTOLOGY_TERM_PATTERN.findall(value))
    return terms


def enrich_candidate(
    row: dict[str, Any], hydration: Hydration | None
) -> dict[str, Any]:
    hydration = hydration or Hydration()
    deduplication_key = ""
    if hydration.structure_sha256:
        deduplication_key = sha256_json(
            [hydration.structure_sha256, row["fastq_set_sha256"]]
        )
    expected_modalities = {
        value for value in row["family_expected_modalities"].split(";") if value
    }
    observed_modalities = set(hydration.modalities)
    if not expected_modalities:
        modality_match_status = "not_required"
    elif hydration.status != "normalized":
        modality_match_status = "not_run"
    elif expected_modalities.intersection(observed_modalities):
        modality_match_status = "matched"
    else:
        modality_match_status = "mismatch"
    return {
        "selection_id": "",
        **row,
        "hydration_status": hydration.status,
        "hydration_message": hydration.message,
        "raw_spec_sha256": hydration.raw_spec_sha256,
        "download_attempts": hydration.download_attempts,
        "raw_seqspec_version": hydration.raw_seqspec_version,
        "normalized_spec_sha256": hydration.normalized_spec_sha256,
        "normalized_seqspec_version": hydration.normalized_seqspec_version,
        "modalities": ";".join(hydration.modalities),
        "modality_match_status": modality_match_status,
        "ontology_terms": ";".join(hydration.ontology_terms),
        "structure_sha256": hydration.structure_sha256,
        "expected_fastq_accessions": ";".join(hydration.expected_fastq_accessions),
        "fastq_mapping_status": hydration.fastq_mapping_status,
        "spec_only_fastq_accessions": ";".join(hydration.spec_only_fastq_accessions),
        "portal_only_fastq_accessions": ";".join(
            hydration.portal_only_fastq_accessions
        ),
        "structural_check_status": hydration.structural_check_status,
        "structural_check_message": hydration.structural_check_message,
        "resource_check_status": hydration.resource_check_status,
        "resource_check_message": hydration.resource_check_message,
        "resource_check_attempts": hydration.resource_check_attempts,
        "deduplication_key": deduplication_key,
        "normalized_spec_path": hydration.normalized_spec_path,
    }


def selection_identity_payload(
    proposals: list[dict[str, Any]],
    selected: list[dict[str, Any]],
    rules_payload: dict[str, Any],
    proposal_count: int,
    candidate_count: int,
    seed: int,
    configuration_sha256: str,
    sequence_sha256: str,
    seqspec_identity: dict[str, Any],
) -> dict[str, Any]:
    proposal_identity_fields = (
        "family_id",
        "proposal_rank",
        "configuration_accession",
        "raw_spec_sha256",
        "normalized_spec_sha256",
        "normalized_seqspec_version",
        "modalities",
        "modality_match_status",
        "review_selection_status",
        "review_exclusion_reason",
    )
    candidate_identity_fields = (
        "family_id",
        "family_rank",
        "configuration_accession",
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "rules": rules_payload,
        "proposal_count_per_family": proposal_count,
        "candidate_count_per_family": candidate_count,
        "selection_seed": seed,
        "configuration_snapshot_sha256": configuration_sha256,
        "sequence_snapshot_sha256": sequence_sha256,
        "seqspec": seqspec_identity,
        "proposals": [
            {key: row.get(key, "") for key in proposal_identity_fields}
            for row in proposals
        ],
        "candidates": [
            {key: row.get(key, "") for key in candidate_identity_fields}
            for row in selected
        ],
    }


def build_validation(
    inventory: list[dict[str, Any]],
    proposals: list[dict[str, Any]],
    selected: list[dict[str, Any]],
    rules: list[FamilyRule],
    proposal_target: int,
    candidate_target: int,
    minimum_final_count: int,
    selection_id: str,
    review_file_status: str,
) -> dict[str, Any]:
    families = []
    for rule in rules:
        eligible = [
            row
            for row in inventory
            if row["family_id"] == rule.family_id
            and row["metadata_eligibility"] == "eligible"
        ]
        proposed = [row for row in proposals if row["family_id"] == rule.family_id]
        modality_qualified = [
            row for row in proposed if not review_exclusion_reason(row)
        ]
        chosen = [row for row in selected if row["family_id"] == rule.family_id]
        eligible_labs = sorted({row["lab"] for row in eligible})
        proposed_labs = sorted({row["lab"] for row in proposed})
        qualified_labs = sorted({row["lab"] for row in modality_qualified})
        chosen_labs = sorted({row["lab"] for row in chosen})
        families.append(
            {
                "family_id": rule.family_id,
                "family_label": rule.label,
                "eligible_configuration_count": len(eligible),
                "eligible_lab_count": len(eligible_labs),
                "eligible_labs": eligible_labs,
                "selected_proposal_count": len(proposed),
                "selected_proposal_lab_count": len(proposed_labs),
                "selected_proposal_labs": proposed_labs,
                "proposal_target": proposal_target,
                "proposal_target_met": len(proposed) == proposal_target,
                "modality_qualified_proposal_count": len(modality_qualified),
                "modality_qualified_lab_count": len(qualified_labs),
                "modality_qualified_labs": qualified_labs,
                "selected_candidate_count": len(chosen),
                "selected_lab_count": len(chosen_labs),
                "selected_labs": chosen_labs,
                "candidate_target": candidate_target,
                "candidate_target_met": len(chosen) == candidate_target,
                "minimum_final_count": minimum_final_count,
                "metadata_minimum_pool_met": len(eligible) >= minimum_final_count,
                "minimum_final_pool_met": len(modality_qualified)
                >= minimum_final_count,
                "two_lab_diversity_possible": len(qualified_labs) >= 2,
                "two_lab_diversity_selected": len(chosen_labs) >= 2,
            }
        )
    dedup_counts = Counter(
        row["deduplication_key"] for row in selected if row.get("deduplication_key")
    )
    minimum_pool_met = all(family["minimum_final_pool_met"] for family in families)
    candidate_minimum_met = all(
        family["selected_candidate_count"] >= minimum_final_count for family in families
    )
    normalization_complete = bool(selected) and all(
        row["hydration_status"] == "normalized"
        or row["modality_match_status"] == "not_required"
        for row in selected
    )
    fastq_mapping_complete = all(
        row["fastq_mapping_status"] == "matched" for row in selected
    )
    modality_match_complete = all(
        row["modality_match_status"] in {"matched", "not_required"} for row in selected
    )
    structural_checks_complete = all(
        row["structural_check_status"] == "passed" for row in selected
    )
    resource_checks_complete = all(
        row["resource_check_status"] == "passed" for row in selected
    )
    freeze_blockers = []
    if not minimum_pool_met:
        freeze_blockers.append(
            "at least one family has fewer modality-qualified proposals than the final cohort minimum"
        )
    if not candidate_minimum_met:
        freeze_blockers.append(
            "at least one family has fewer selected review candidates than the final cohort minimum"
        )
    if not normalization_complete:
        freeze_blockers.append(
            "selected seqspec files are not all normalized and structurally fingerprinted"
        )
    if not fastq_mapping_complete:
        freeze_blockers.append(
            "selected seqspec expected FASTQs do not all match portal-linked FASTQs"
        )
    if not modality_match_complete:
        freeze_blockers.append(
            "selected seqspec modalities do not all match their proposed assay families"
        )
    if not structural_checks_complete:
        freeze_blockers.append(
            "selected normalized seqspec files do not all pass schema and structural checks"
        )
    if not resource_checks_complete:
        freeze_blockers.append(
            "selected normalized seqspec files do not all have accessible declared resources"
        )
    freeze_blockers.extend(
        [
            "two independent reviews and any adjudication are incomplete",
            "calibration/evaluation split is intentionally unassigned",
        ]
    )
    if review_file_status == "stale_preserved":
        freeze_blockers.append(
            "preserved human review sheet belongs to a different selection id"
        )
    qualified_by_family = {
        rule.family_id: sum(
            row["family_id"] == rule.family_id
            and row["hydration_status"] == "normalized"
            and row["structural_check_status"] == "passed"
            and row["resource_check_status"] == "passed"
            and row["fastq_mapping_status"] == "matched"
            and row["modality_match_status"] in {"matched", "not_required"}
            for row in selected
        )
        for rule in rules
    }
    for family in families:
        qualified = qualified_by_family[family["family_id"]]
        family["qualified_candidate_count"] = qualified
        family["minimum_qualified_pool_met"] = qualified >= minimum_final_count
    qualified_pools_met = all(
        family["minimum_qualified_pool_met"] for family in families
    )
    if not qualified_pools_met:
        freeze_blockers.insert(
            1,
            "at least one family has fewer fully qualified candidates than the final cohort minimum",
        )
    return {
        "cohort_candidate_schema_version": SCHEMA_VERSION,
        "selection_id": selection_id,
        "frozen": False,
        "split_assigned": False,
        "selection_uses_seqcheck_results": False,
        "human_review_required": True,
        "required_independent_reviewers": 2,
        "review_file_status": review_file_status,
        "inventory_row_count": len(inventory),
        "proposal_row_count": len(proposals),
        "proposal_unique_configuration_count": len(
            {row["configuration_accession"] for row in proposals}
        ),
        "proposal_review_selection_status_counts": dict(
            sorted(Counter(row["review_selection_status"] for row in proposals).items())
        ),
        "proposal_review_exclusion_reason_counts": dict(
            sorted(
                Counter(
                    row["review_exclusion_reason"]
                    for row in proposals
                    if row["review_exclusion_reason"]
                ).items()
            )
        ),
        "proposal_normalization_failure_count": sum(
            row["hydration_status"] == "failed" for row in proposals
        ),
        "proposal_normalization_complete_count": sum(
            row["hydration_status"] == "normalized" for row in proposals
        ),
        "proposal_modality_match_status_counts": dict(
            sorted(Counter(row["modality_match_status"] for row in proposals).items())
        ),
        "selected_candidate_row_count": len(selected),
        "selected_unique_configuration_count": len(
            {row["configuration_accession"] for row in selected}
        ),
        "ambiguous_selected_row_count": sum(
            bool(row["ambiguous_family"]) for row in selected
        ),
        "normalization_failure_count": sum(
            row["hydration_status"] == "failed" for row in selected
        ),
        "normalization_complete_count": sum(
            row["hydration_status"] == "normalized" for row in selected
        ),
        "fastq_mapping_matched_count": sum(
            row["fastq_mapping_status"] == "matched" for row in selected
        ),
        "fastq_mapping_issue_count": sum(
            row["fastq_mapping_status"] not in {"matched", "not_run"}
            for row in selected
        ),
        "modality_match_status_counts": dict(
            sorted(Counter(row["modality_match_status"] for row in selected).items())
        ),
        "structural_check_passed_count": sum(
            row["structural_check_status"] == "passed" for row in selected
        ),
        "structural_check_failed_count": sum(
            row["structural_check_status"] not in {"passed", "not_run"}
            for row in selected
        ),
        "structural_check_status_counts": dict(
            sorted(Counter(row["structural_check_status"] for row in selected).items())
        ),
        "resource_check_passed_count": sum(
            row["resource_check_status"] == "passed" for row in selected
        ),
        "resource_check_failed_count": sum(
            row["resource_check_status"] not in {"passed", "not_run"}
            for row in selected
        ),
        "resource_check_status_counts": dict(
            sorted(Counter(row["resource_check_status"] for row in selected).items())
        ),
        "all_minimum_final_pools_met": minimum_pool_met,
        "all_candidate_minimums_met": candidate_minimum_met,
        "all_minimum_qualified_pools_met": qualified_pools_met,
        "all_proposal_targets_met": all(
            family["proposal_target_met"] for family in families
        ),
        "all_candidate_targets_met": all(
            family["candidate_target_met"] for family in families
        ),
        "duplicate_group_count": sum(count > 1 for count in dedup_counts.values()),
        "families": families,
        "ready_to_freeze": False,
        "freeze_blockers": freeze_blockers,
    }


def inventory_fields(rows: list[dict[str, Any]]) -> list[str]:
    return list(rows[0].keys()) if rows else []


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_review_selection_ids(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return sorted(
            {
                row.get("selection_id", "").strip()
                for row in csv.DictReader(handle)
                if row.get("selection_id", "").strip()
            }
        )


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def snapshot_identity(
    path: Path, payload: dict[str, Any], count: int
) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "sha256": file_sha256(path),
        "record_count": count,
        "portal_total": payload.get("total"),
        "query_id": payload.get("@id", ""),
    }


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def absolute_url(root: str, href: str) -> str:
    return urllib.parse.urljoin(root, href)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
