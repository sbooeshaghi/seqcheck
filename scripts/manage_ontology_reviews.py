#!/usr/bin/env python3
"""Prepare blinded ontology reviews and merge two completed review sheets."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

try:
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
TOOL_VERSION = "0.1.0"
SURVEY_SCHEMA_VERSION = "0.1.0"
UNKNOWN_TERM = "RGN:unknown:unclassified"
PACKAGE_FIELDS = ("review_package_id", "review_slot", "review_item_id")
CONTEXT_FIELDS = (
    "configuration_accession",
    "family_labels",
    "assay_term",
    "preferred_assay_titles",
    "assay_id",
    "assay_name",
    "assay_description",
    "raw_seqspec_version",
    "modality",
    "modality_read_ids",
    "modality_primer_ids",
    "region_id",
    "region_name",
    "sequence_type",
    "min_len",
    "max_len",
    "depth",
    "is_leaf",
    "parent_region_id",
    "parent_region_name",
    "path_region_ids",
    "path_region_names",
    "original_region_type_json",
    "original_region_type_labels",
    "has_onlist",
)
EDITABLE_FIELDS = (
    "reviewer",
    "review_date",
    "reviewed_ontology_terms",
    "context_sufficient",
    "confidence",
    "notes",
)
BLINDED_FIELDS = (
    "survey_id",
    "region_key",
    "ontology_terms",
    "mapping_source",
    "unknown_mapping",
    "multi_term_mapping",
    "review_stratum",
    "stratum_population",
    "stratum_sample_size",
    "inclusion_probability",
    "sampling_weight",
    "selection_rank",
)
PRIVATE_FIELDS = (
    "survey_id",
    "region_key",
    "reviewer_1_item_id",
    "reviewer_2_item_id",
    "ontology_terms",
    "review_stratum",
    "stratum_population",
    "stratum_sample_size",
    "inclusion_probability",
    "sampling_weight",
)
REVIEW_RESULT_FIELDS = (
    "reviewer_1",
    "reviewer_1_date",
    "reviewer_1_ontology_terms",
    "reviewer_1_context_sufficient",
    "reviewer_1_confidence",
    "reviewer_1_notes",
    "reviewer_1_registry_exact_match",
    "reviewer_2",
    "reviewer_2_date",
    "reviewer_2_ontology_terms",
    "reviewer_2_context_sufficient",
    "reviewer_2_confidence",
    "reviewer_2_notes",
    "reviewer_2_registry_exact_match",
    "exact_term_set_agreement",
    "needs_adjudication",
    "consensus_ontology_terms",
    "adjudicator",
    "adjudication_date",
    "adjudicated_ontology_terms",
    "adjudication_rationale",
)
TERM_REFERENCE_FIELDS = (
    "ontology_term",
    "label",
    "role",
    "role_scope",
    "target",
    "status",
    "definition",
    "typical_sequence_types_json",
    "examples_json",
    "does_not_mean_json",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare or merge blinded region ontology reviews."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    add_source_arguments(prepare)
    prepare.add_argument("--output-root", required=True, type=Path)
    merge = subparsers.add_parser("merge")
    add_source_arguments(merge)
    merge.add_argument("--reviewer-1-package", required=True, type=Path)
    merge.add_argument("--reviewer-1-sheet", required=True, type=Path)
    merge.add_argument("--reviewer-2-package", required=True, type=Path)
    merge.add_argument("--reviewer-2-sheet", required=True, type=Path)
    merge.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def add_source_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--survey-manifest", required=True, type=Path)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--yq-bin", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=30)


def main() -> int:
    args = parse_args()
    common = {
        "survey_manifest_path": args.survey_manifest.resolve(),
        "registry_path": args.registry.resolve(),
        "protocol_path": args.protocol.resolve(),
        "yq_bin": args.yq_bin.resolve(),
        "timeout_seconds": args.timeout_seconds,
    }
    try:
        if args.command == "prepare":
            prepare_review_packages(output_root=args.output_root.resolve(), **common)
        else:
            merge_review_packages(
                reviewer_1_package=args.reviewer_1_package.resolve(),
                reviewer_1_sheet=args.reviewer_1_sheet.resolve(),
                reviewer_2_package=args.reviewer_2_package.resolve(),
                reviewer_2_sheet=args.reviewer_2_sheet.resolve(),
                output_root=args.output_root.resolve(),
                **common,
            )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"manage_ontology_reviews: {error}", file=sys.stderr)
        return 1
    return 0


def prepare_review_packages(
    *,
    survey_manifest_path: Path,
    registry_path: Path,
    protocol_path: Path,
    yq_bin: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    refuse_output_root(output_root)
    source = load_source(
        survey_manifest_path=survey_manifest_path,
        registry_path=registry_path,
        protocol_path=protocol_path,
        yq_bin=yq_bin,
        timeout_seconds=timeout_seconds,
    )
    script = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    yq_identity = runtime.executable_identity(yq_bin, timeout_seconds=timeout_seconds)
    package_manifests = []
    item_ids_by_slot = {}
    for slot in (1, 2):
        package, item_ids = prepare_one_package(
            slot=slot,
            source=source,
            script=script,
            yq_identity=yq_identity,
            output_root=output_root,
        )
        package_manifests.append(package)
        item_ids_by_slot[slot] = item_ids

    private_rows = []
    for row in source["sample_rows"]:
        private_rows.append(
            {
                "survey_id": source["survey_id"],
                "region_key": row["region_key"],
                "reviewer_1_item_id": item_ids_by_slot[1][row["region_key"]],
                "reviewer_2_item_id": item_ids_by_slot[2][row["region_key"]],
                **{field: row[field] for field in PRIVATE_FIELDS[4:]},
            }
        )
    private_path = output_root / "internal" / "review_key.csv"
    runtime.write_csv(private_path, private_rows, list(PRIVATE_FIELDS))
    manifest_path = output_root / "manifests" / "ontology_review_packages.json"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "survey_id": source["survey_id"],
        "generated_at": utc_now(),
        "review_rows": len(source["sample_rows"]),
        "review_slots": 2,
        "blinded_fields": list(BLINDED_FIELDS),
        "tool": script,
        "inputs": source["input_identities"],
        "packages": [
            {
                "review_slot": value["review_slot"],
                "package_id": value["package_id"],
                "manifest": runtime.file_identity(Path(value["manifest_path"])),
            }
            for value in package_manifests
        ],
        "private_review_key": runtime.file_identity(private_path),
        "private_review_key_sharing": "do_not_share_with_reviewers",
        "manifest_path": str(manifest_path),
        "frozen": False,
    }
    runtime.write_json(manifest_path, manifest)
    print(
        f"prepared two blinded ontology review packages for "
        f"{len(source['sample_rows'])} regions (survey_id={source['survey_id']})"
    )
    return manifest


def prepare_one_package(
    *,
    slot: int,
    source: dict[str, Any],
    script: dict[str, Any],
    yq_identity: dict[str, Any],
    output_root: Path,
) -> tuple[dict[str, Any], dict[str, str]]:
    seed = source["package_seeds"][slot]
    identity = package_identity(slot, seed, source, script["sha256"])
    package_id = runtime.sha256_json(identity)[:16]
    rows, item_ids = expected_package_rows(
        slot=slot,
        source=source,
        package_id=package_id,
    )
    package_dir = output_root / f"reviewer_{slot}"
    sheet_path = package_dir / "ontology_review.csv"
    reference_path = package_dir / "ontology_terms.csv"
    instructions_path = package_dir / "INSTRUCTIONS.md"
    manifest_path = package_dir / "review_package.json"
    runtime.write_csv(
        sheet_path,
        rows,
        [*PACKAGE_FIELDS, *CONTEXT_FIELDS, *EDITABLE_FIELDS],
    )
    runtime.write_csv(
        reference_path, source["term_reference"], list(TERM_REFERENCE_FIELDS)
    )
    instructions_path.parent.mkdir(parents=True, exist_ok=True)
    instructions_path.write_text(review_instructions(slot), encoding="utf-8")
    evidence_sha256 = evidence_hash(rows, [*PACKAGE_FIELDS, *CONTEXT_FIELDS])
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "package_id": package_id,
        "survey_id": source["survey_id"],
        "generated_at": utc_now(),
        "review_slot": slot,
        "selection_seed": seed,
        "review_rows": len(rows),
        "blinded_fields": list(BLINDED_FIELDS),
        "fields": {
            "package": list(PACKAGE_FIELDS),
            "context": list(CONTEXT_FIELDS),
            "editable": list(EDITABLE_FIELDS),
        },
        "review_evidence_sha256": evidence_sha256,
        "tool": script,
        "tools": {"yq": yq_identity},
        "inputs": source["input_identities"],
        "prepared_sheet": runtime.file_identity(sheet_path),
        "term_reference": runtime.file_identity(reference_path),
        "instructions": runtime.file_identity(instructions_path),
        "manifest_path": str(manifest_path),
        "frozen": False,
    }
    runtime.write_json(manifest_path, manifest)
    return manifest, item_ids


def merge_review_packages(
    *,
    survey_manifest_path: Path,
    registry_path: Path,
    protocol_path: Path,
    yq_bin: Path,
    reviewer_1_package: Path,
    reviewer_1_sheet: Path,
    reviewer_2_package: Path,
    reviewer_2_sheet: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    refuse_output_root(output_root)
    source = load_source(
        survey_manifest_path=survey_manifest_path,
        registry_path=registry_path,
        protocol_path=protocol_path,
        yq_bin=yq_bin,
        timeout_seconds=timeout_seconds,
    )
    first = load_completed_review(
        package_path=reviewer_1_package,
        sheet_path=reviewer_1_sheet,
        expected_slot=1,
        source=source,
    )
    second = load_completed_review(
        package_path=reviewer_2_package,
        sheet_path=reviewer_2_sheet,
        expected_slot=2,
        source=source,
    )
    if first["reviewer"].casefold() == second["reviewer"].casefold():
        raise ValueError("ontology reviewers must be distinct people")

    merged_rows = []
    for sample in source["sample_rows"]:
        left = first["rows_by_region_key"][sample["region_key"]]
        right = second["rows_by_region_key"][sample["region_key"]]
        truth = canonical_terms(sample["ontology_terms"], source)
        left_terms = left["_canonical_terms"]
        right_terms = right["_canonical_terms"]
        agreement = left_terms == right_terms
        merged_rows.append(
            {
                **sample,
                **prefixed_review(left, 1, truth),
                **prefixed_review(right, 2, truth),
                "exact_term_set_agreement": agreement,
                "needs_adjudication": not agreement,
                "consensus_ontology_terms": left_terms if agreement else "",
                "adjudicator": "",
                "adjudication_date": "",
                "adjudicated_ontology_terms": "",
                "adjudication_rationale": "",
            }
        )
    merged_path = output_root / "tables" / "ontology_reviews.csv"
    runtime.write_csv(
        merged_path,
        merged_rows,
        [*source["sample_fields"], *REVIEW_RESULT_FIELDS],
    )
    script = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    disagreements = sum(row["needs_adjudication"] for row in merged_rows)
    manifest_path = output_root / "manifests" / "ontology_reviews.json"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "merge_id": runtime.sha256_json(
            {
                "survey_id": source["survey_id"],
                "reviewer_1_package": runtime.file_sha256(reviewer_1_package),
                "reviewer_1_sheet": runtime.file_sha256(reviewer_1_sheet),
                "reviewer_2_package": runtime.file_sha256(reviewer_2_package),
                "reviewer_2_sheet": runtime.file_sha256(reviewer_2_sheet),
                "tool": runtime.functional_script_identity(script),
            }
        )[:16],
        "survey_id": source["survey_id"],
        "generated_at": utc_now(),
        "review_rows": len(merged_rows),
        "reviewers": {
            "reviewer_1": first["reviewer"],
            "reviewer_2": second["reviewer"],
        },
        "exact_term_set_agreements": len(merged_rows) - disagreements,
        "adjudication_required": disagreements,
        "adjudication_complete": disagreements == 0,
        "frozen": False,
        "tool": script,
        "inputs": {
            **source["input_identities"],
            "reviewer_1_package": runtime.file_identity(reviewer_1_package),
            "reviewer_1_sheet": runtime.file_identity(reviewer_1_sheet),
            "reviewer_2_package": runtime.file_identity(reviewer_2_package),
            "reviewer_2_sheet": runtime.file_identity(reviewer_2_sheet),
        },
        "outputs": {"ontology_reviews": runtime.file_identity(merged_path)},
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    print(
        f"merged {len(merged_rows)} blinded ontology reviews "
        f"(merge_id={manifest['merge_id']}, adjudication_required={disagreements})"
    )
    return manifest


def load_source(
    *,
    survey_manifest_path: Path,
    registry_path: Path,
    protocol_path: Path,
    yq_bin: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    manifest = runtime.load_json(survey_manifest_path)
    if manifest.get("schema_version") != SURVEY_SCHEMA_VERSION:
        raise ValueError("ontology survey manifest schema is unsupported")
    if manifest.get("valid") is not True:
        raise ValueError("ontology survey is not valid")
    survey_id = str(manifest.get("survey_id", "")).strip()
    if not survey_id:
        raise ValueError("ontology survey id is missing")
    verify_manifest_input(manifest, "registry", registry_path)
    verify_manifest_input(manifest, "protocol", protocol_path)
    sample_identity = manifest.get("outputs", {}).get("review_sample")
    if not isinstance(sample_identity, dict):
        raise ValueError("ontology survey review sample identity is missing")
    sample_path = Path(str(sample_identity.get("path", ""))).resolve()
    verify_file_identity(sample_path, sample_identity, "ontology review sample")
    sample_rows, sample_fields = read_csv(sample_path)
    expected_rows = manifest.get("counts", {}).get("review_sample")
    if expected_rows != len(sample_rows) or not sample_rows:
        raise ValueError("ontology review sample row count differs")
    required = {
        "survey_id",
        "region_key",
        "ontology_terms",
        "mapping_source",
        "unknown_mapping",
        "multi_term_mapping",
        "review_stratum",
        "stratum_population",
        "stratum_sample_size",
        "inclusion_probability",
        "sampling_weight",
        "selection_rank",
        *CONTEXT_FIELDS,
    }
    missing = sorted(required - set(sample_fields))
    if missing:
        raise ValueError(
            "ontology review sample fields are missing: " + ", ".join(missing)
        )
    if len({row["region_key"] for row in sample_rows}) != len(sample_rows):
        raise ValueError("ontology review sample region keys are not unique")
    if any(row["survey_id"] != survey_id for row in sample_rows):
        raise ValueError("ontology review sample survey identifiers differ")

    protocol = runtime.load_json(protocol_path)
    review_contract = validate_review_protocol(protocol, len(sample_rows))
    registry = load_yaml_json(registry_path, yq_bin, timeout_seconds)
    term_reference, term_order = build_term_reference(registry)
    declared_terms = set(term_order)
    for row in sample_rows:
        terms = split_terms(row["ontology_terms"])
        if (
            not terms
            or len(terms) != len(set(terms))
            or not set(terms) <= declared_terms
        ):
            raise ValueError(f"{row['region_key']}: survey ontology terms are invalid")
    parity = manifest.get("registry_runtime_parity", {})
    if parity.get("matched_labels") != parity.get("legacy_labels"):
        raise ValueError("ontology survey registry runtime parity is incomplete")
    return {
        "survey_id": survey_id,
        "sample_rows": sample_rows,
        "sample_fields": sample_fields,
        "sample_path": sample_path,
        "term_reference": term_reference,
        "term_order": term_order,
        "declared_terms": declared_terms,
        "package_seeds": review_contract["package_seeds"],
        "context_values": review_contract["context_values"],
        "confidence_values": review_contract["confidence_values"],
        "input_identities": {
            "survey_manifest": runtime.file_identity(survey_manifest_path),
            "review_sample": runtime.file_identity(sample_path),
            "registry": runtime.file_identity(registry_path),
            "protocol": runtime.file_identity(protocol_path),
        },
    }


def validate_review_protocol(
    protocol: dict[str, Any], sample_rows: int
) -> dict[str, Any]:
    if protocol.get("schema_version") != SURVEY_SCHEMA_VERSION:
        raise ValueError("ontology review protocol schema is unsupported")
    review = protocol.get("manual_review", {})
    if review.get("sample_size") != sample_rows:
        raise ValueError("ontology review protocol sample size differs")
    if review.get("reviewers") != 2 or review.get("blind_registry_mapping") is not True:
        raise ValueError("ontology review must use two blinded reviewers")
    if review.get("agreement_endpoint") != "order_insensitive_exact_term_set":
        raise ValueError("ontology review agreement endpoint is invalid")
    if review.get("adjudication_required_for") != "term_set_disagreement":
        raise ValueError("ontology review adjudication policy is invalid")
    seeds = review.get("review_package_seeds", {})
    if set(seeds) != {"reviewer_1", "reviewer_2"}:
        raise ValueError("ontology review package seeds are invalid")
    if any(
        isinstance(value, bool) or not isinstance(value, int)
        for value in seeds.values()
    ):
        raise ValueError("ontology review package seed is invalid")
    if seeds["reviewer_1"] == seeds["reviewer_2"]:
        raise ValueError("ontology review package seeds must differ")
    response = review.get("response_contract", {})
    if response.get("term_delimiter") != ";":
        raise ValueError("ontology review term delimiter is invalid")
    if response.get("unknown_term_must_be_exclusive") is not True:
        raise ValueError("ontology unknown term exclusivity is not enabled")
    if response.get("insufficient_context_requires_unknown") is not True:
        raise ValueError("ontology insufficient-context policy is not enabled")
    context_values = response.get("context_sufficient_values")
    confidence_values = response.get("confidence_values")
    if context_values != ["yes", "no"]:
        raise ValueError("ontology context response values are invalid")
    if confidence_values != ["high", "medium", "low"]:
        raise ValueError("ontology confidence response values are invalid")
    return {
        "package_seeds": {1: seeds["reviewer_1"], 2: seeds["reviewer_2"]},
        "context_values": set(context_values),
        "confidence_values": set(confidence_values),
    }


def build_term_reference(
    registry: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[str]]:
    if registry.get("ontology_id") != "seqspec-region-ontology":
        raise ValueError("ontology registry id is invalid")
    roles = registry.get("roles")
    terms = registry.get("terms")
    if not isinstance(roles, dict) or not isinstance(terms, dict) or not terms:
        raise ValueError("ontology registry term catalog is invalid")
    rows = []
    order = []
    for term, value in terms.items():
        if not isinstance(value, dict) or value.get("status") != "active":
            continue
        role = str(value.get("role", ""))
        if role not in roles:
            raise ValueError(f"ontology term has an undeclared role: {term}")
        rows.append(
            {
                "ontology_term": term,
                "label": value.get("label", ""),
                "role": role,
                "role_scope": roles[role].get("scope", ""),
                "target": value.get("target", ""),
                "status": value.get("status", ""),
                "definition": value.get("definition", ""),
                "typical_sequence_types_json": runtime.canonical_json(
                    value.get("typical_sequence_types", [])
                ),
                "examples_json": runtime.canonical_json(value.get("examples", [])),
                "does_not_mean_json": runtime.canonical_json(
                    value.get("does_not_mean", [])
                ),
            }
        )
        order.append(term)
    if UNKNOWN_TERM not in order:
        raise ValueError("ontology term catalog omits the unknown term")
    return rows, order


def load_completed_review(
    *,
    package_path: Path,
    sheet_path: Path,
    expected_slot: int,
    source: dict[str, Any],
) -> dict[str, Any]:
    package = runtime.load_json(package_path)
    validate_package(package_path, package, expected_slot, source)
    rows, fields = read_csv(sheet_path)
    expected_fields = [*PACKAGE_FIELDS, *CONTEXT_FIELDS, *EDITABLE_FIELDS]
    if fields != expected_fields:
        raise ValueError(f"reviewer {expected_slot} sheet fields changed")
    if len(rows) != len(source["sample_rows"]):
        raise ValueError(f"reviewer {expected_slot} sheet row count changed")
    package_id = package["package_id"]
    seed = source["package_seeds"][expected_slot]
    expected_by_item = {}
    for sample in source["sample_rows"]:
        item_id = review_item_id(
            source["survey_id"], sample["region_key"], expected_slot, seed
        )
        expected_by_item[item_id] = sample
    rows_by_item = index_unique(rows, "review_item_id", f"reviewer {expected_slot}")
    if set(rows_by_item) != set(expected_by_item):
        raise ValueError(f"reviewer {expected_slot} review item ids changed")

    reviewers = set()
    dates = set()
    rows_by_region_key = {}
    for item_id, expected in expected_by_item.items():
        row = rows_by_item[item_id]
        if row["review_package_id"] != package_id:
            raise ValueError(f"reviewer {expected_slot} package id changed in sheet")
        if row["review_slot"] != str(expected_slot):
            raise ValueError(f"reviewer {expected_slot} review slot changed in sheet")
        changed = [field for field in CONTEXT_FIELDS if row[field] != expected[field]]
        if changed:
            raise ValueError(
                f"reviewer {expected_slot} changed blinded context: {', '.join(changed)}"
            )
        reviewer = row["reviewer"].strip()
        if not reviewer:
            raise ValueError(f"reviewer {expected_slot} identity is missing")
        reviewers.add(reviewer)
        review_date = row["review_date"].strip()
        validate_date(review_date, f"reviewer {expected_slot}")
        dates.add(review_date)
        context = row["context_sufficient"].strip().lower()
        confidence = row["confidence"].strip().lower()
        if context not in source["context_values"]:
            raise ValueError(f"reviewer {expected_slot} context response is invalid")
        if confidence not in source["confidence_values"]:
            raise ValueError(f"reviewer {expected_slot} confidence response is invalid")
        terms = canonical_terms(row["reviewed_ontology_terms"], source)
        values = split_terms(terms)
        if UNKNOWN_TERM in values and len(values) != 1:
            raise ValueError(
                f"reviewer {expected_slot} combined unknown with precise terms"
            )
        if context == "no" and values != [UNKNOWN_TERM]:
            raise ValueError(
                f"reviewer {expected_slot} must use unknown when context is insufficient"
            )
        rows_by_region_key[expected["region_key"]] = {
            **row,
            "context_sufficient": context,
            "confidence": confidence,
            "_canonical_terms": terms,
        }
    if len(reviewers) != 1:
        raise ValueError(f"reviewer {expected_slot} sheet names multiple reviewers")
    return {
        "reviewer": next(iter(reviewers)),
        "dates": sorted(dates),
        "rows_by_region_key": rows_by_region_key,
        "confidence_counts": dict(
            sorted(
                Counter(
                    row["confidence"] for row in rows_by_region_key.values()
                ).items()
            )
        ),
    }


def validate_package(
    package_path: Path,
    package: dict[str, Any],
    expected_slot: int,
    source: dict[str, Any],
) -> None:
    if package.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"reviewer {expected_slot} package schema is unsupported")
    if package.get("review_slot") != expected_slot:
        raise ValueError(f"reviewer {expected_slot} package slot changed")
    if package.get("survey_id") != source["survey_id"]:
        raise ValueError(f"reviewer {expected_slot} package survey id changed")
    if package.get("review_rows") != len(source["sample_rows"]):
        raise ValueError(f"reviewer {expected_slot} package row count changed")
    if package.get("blinded_fields") != list(BLINDED_FIELDS):
        raise ValueError(f"reviewer {expected_slot} package blinding contract changed")
    if package.get("fields") != {
        "package": list(PACKAGE_FIELDS),
        "context": list(CONTEXT_FIELDS),
        "editable": list(EDITABLE_FIELDS),
    }:
        raise ValueError(f"reviewer {expected_slot} package field contract changed")
    inputs = package.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError(f"reviewer {expected_slot} package inputs are invalid")
    for name, expected in source["input_identities"].items():
        observed = inputs.get(name)
        if (
            not isinstance(observed, dict)
            or observed.get("sha256") != expected["sha256"]
        ):
            raise ValueError(f"reviewer {expected_slot} package input changed: {name}")
    seed = source["package_seeds"][expected_slot]
    if package.get("selection_seed") != seed:
        raise ValueError(f"reviewer {expected_slot} package seed changed")
    tool = package.get("tool")
    if not isinstance(tool, dict) or not tool.get("sha256"):
        raise ValueError(f"reviewer {expected_slot} package tool identity is invalid")
    expected_id = runtime.sha256_json(
        package_identity(expected_slot, seed, source, tool["sha256"])
    )[:16]
    if package.get("package_id") != expected_id:
        raise ValueError(f"reviewer {expected_slot} package id is invalid")
    expected_rows, _ = expected_package_rows(
        slot=expected_slot,
        source=source,
        package_id=expected_id,
    )
    expected_evidence_hash = evidence_hash(
        expected_rows, [*PACKAGE_FIELDS, *CONTEXT_FIELDS]
    )
    if package.get("review_evidence_sha256") != expected_evidence_hash:
        raise ValueError(f"reviewer {expected_slot} package evidence hash changed")
    reference_identity = package.get("term_reference")
    instructions_identity = package.get("instructions")
    if not isinstance(reference_identity, dict) or not isinstance(
        instructions_identity, dict
    ):
        raise ValueError(f"reviewer {expected_slot} package references are invalid")
    reference_path = package_path.parent / "ontology_terms.csv"
    instructions_path = package_path.parent / "INSTRUCTIONS.md"
    verify_file_identity(
        reference_path, reference_identity, f"reviewer {expected_slot} term reference"
    )
    verify_file_identity(
        instructions_path,
        instructions_identity,
        f"reviewer {expected_slot} instructions",
    )
    reference_rows, reference_fields = read_csv(reference_path)
    if reference_fields != list(TERM_REFERENCE_FIELDS) or reference_rows != [
        {field: str(row[field]) for field in TERM_REFERENCE_FIELDS}
        for row in source["term_reference"]
    ]:
        raise ValueError(f"reviewer {expected_slot} term reference changed")


def expected_package_rows(
    *, slot: int, source: dict[str, Any], package_id: str
) -> tuple[list[dict[str, Any]], dict[str, str]]:
    seed = source["package_seeds"][slot]
    item_ids = {
        row["region_key"]: review_item_id(
            source["survey_id"], row["region_key"], slot, seed
        )
        for row in source["sample_rows"]
    }
    rows = [
        {
            "review_package_id": package_id,
            "review_slot": slot,
            "review_item_id": item_ids[row["region_key"]],
            **{field: row[field] for field in CONTEXT_FIELDS},
            **{field: "" for field in EDITABLE_FIELDS},
        }
        for row in source["sample_rows"]
    ]
    rows.sort(
        key=lambda row: hashlib.sha256(
            f"{seed}:order:{row['review_item_id']}".encode()
        ).hexdigest()
    )
    return rows, item_ids


def package_identity(
    slot: int, seed: int, source: dict[str, Any], script_sha256: str
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "survey_id": source["survey_id"],
        "review_slot": slot,
        "selection_seed": seed,
        "review_sample_sha256": source["input_identities"]["review_sample"]["sha256"],
        "registry_sha256": source["input_identities"]["registry"]["sha256"],
        "protocol_sha256": source["input_identities"]["protocol"]["sha256"],
        "script_sha256": script_sha256,
    }


def review_item_id(survey_id: str, region_key: str, slot: int, seed: int) -> str:
    return runtime.sha256_json(
        {
            "survey_id": survey_id,
            "region_key": region_key,
            "review_slot": slot,
            "seed": seed,
        }
    )[:16]


def prefixed_review(
    row: dict[str, Any], slot: int, registry_terms: str
) -> dict[str, Any]:
    return {
        f"reviewer_{slot}": row["reviewer"].strip(),
        f"reviewer_{slot}_date": row["review_date"].strip(),
        f"reviewer_{slot}_ontology_terms": row["_canonical_terms"],
        f"reviewer_{slot}_context_sufficient": row["context_sufficient"],
        f"reviewer_{slot}_confidence": row["confidence"],
        f"reviewer_{slot}_notes": row["notes"].strip(),
        f"reviewer_{slot}_registry_exact_match": row["_canonical_terms"]
        == registry_terms,
    }


def canonical_terms(value: str, source: dict[str, Any]) -> str:
    terms = split_terms(value)
    if not terms:
        raise ValueError("reviewed ontology terms are missing")
    if len(terms) != len(set(terms)):
        raise ValueError("reviewed ontology terms contain duplicates")
    undeclared = sorted(set(terms) - source["declared_terms"])
    if undeclared:
        raise ValueError(
            "reviewed ontology terms are undeclared: " + ", ".join(undeclared)
        )
    present = set(terms)
    return ";".join(term for term in source["term_order"] if term in present)


def split_terms(value: str) -> list[str]:
    return [term.strip() for term in str(value).split(";") if term.strip()]


def verify_manifest_input(manifest: dict[str, Any], name: str, path: Path) -> None:
    identity = manifest.get("inputs", {}).get(name)
    if not isinstance(identity, dict):
        raise ValueError(f"ontology survey input identity is missing: {name}")
    verify_file_identity(path, identity, f"ontology survey {name}")


def verify_file_identity(path: Path, identity: dict[str, Any], label: str) -> None:
    if not path.is_file():
        raise ValueError(f"{label} file is missing: {path}")
    if runtime.file_sha256(path) != identity.get("sha256"):
        raise ValueError(f"{label} hash changed")
    if path.stat().st_size != identity.get("size_bytes"):
        raise ValueError(f"{label} size changed")


def load_yaml_json(path: Path, yq_bin: Path, timeout_seconds: int) -> dict[str, Any]:
    completed = subprocess.run(
        [str(yq_bin), "-o=json", ".", str(path)],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip()
        raise ValueError(f"could not parse {path}: {message}")
    try:
        value = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise ValueError(f"could not parse yq JSON for {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"expected a YAML object in {path}")
    return value


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    import csv

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if not fields or len(fields) != len(set(fields)):
        raise ValueError(f"CSV fields are missing or duplicated: {path}")
    return rows, fields


def index_unique(
    rows: list[dict[str, str]], field: str, label: str
) -> dict[str, dict[str, str]]:
    result = {}
    for row in rows:
        key = row.get(field, "").strip()
        if not key:
            raise ValueError(f"{label} row is missing {field}")
        if key in result:
            raise ValueError(f"{label} rows duplicate {field}: {key}")
        result[key] = row
    return result


def evidence_hash(rows: list[dict[str, Any]], fields: list[str]) -> str:
    return runtime.sha256_json(
        [[str(row.get(field, "")) for field in fields] for row in rows]
    )


def validate_date(value: str, label: str) -> None:
    try:
        parsed = date.fromisoformat(value)
    except ValueError as error:
        raise ValueError(f"{label} date must use YYYY-MM-DD") from error
    if parsed.isoformat() != value:
        raise ValueError(f"{label} date must use YYYY-MM-DD")


def refuse_output_root(path: Path) -> None:
    if path.exists():
        raise ValueError(f"refusing to overwrite existing output root: {path}")


def review_instructions(slot: int) -> str:
    return f"""# Region ontology review {slot}

Review each row independently. The proposed registry mapping, sampling stratum,
and sampling weight are hidden from this sheet.

1. Use `ontology_terms.csv` as the allowed term catalog.
2. Enter one or more term IDs in `reviewed_ontology_terms`, separated by `;`.
3. Assign terms for the direct nominal role of this read interval. A `classify`
   term describes what an aggregation of these intervals classifies.
4. If the supplied assay and region context cannot support a precise assignment,
   enter only `RGN:unknown:unclassified` and set `context_sufficient` to `no`.
5. Set `confidence` to `high`, `medium`, or `low`. Notes are optional.
6. Use one reviewer name throughout the sheet and an ISO date (`YYYY-MM-DD`).

Do not consult the seqspec region ontology legacy-label mapping while reviewing.
"""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
