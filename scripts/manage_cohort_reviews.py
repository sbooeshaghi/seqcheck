#!/usr/bin/env python3
"""Prepare independent cohort review sheets and merge locked reviews."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

try:
    import cohort_corrections as corrections_module
except ModuleNotFoundError:
    from scripts import cohort_corrections as corrections_module

SUPPLEMENTAL_EVIDENCE_FIELDS = corrections_module.EVIDENCE_FIELDS
load_registry = corrections_module.load_registry


SCHEMA_VERSION = "0.1.0"
CANDIDATE_SCHEMA_VERSION = "0.2.0"
REVIEW_DECISIONS = {"include", "exclude", "inconclusive"}
KEY_FIELDS = ("family_id", "configuration_accession")
PACKAGE_FIELDS = ("review_package_id", "review_slot")
REVIEW_INPUT_FIELDS = (
    "reviewer",
    "decision",
    "rationale",
    "protocol_url",
    "date",
)
COMBINED_REVIEW_FIELDS = (
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
)
CORRECTION_ARTIFACT_FIELDS = ("original", "corrected", "diff")
REVIEW_PROTOCOL_PATH = (
    Path(__file__).resolve().parents[1] / "docs" / "COHORT_REVIEW_PROTOCOL.md"
)
SAFE_PATH_COMPONENT = re.compile(r"^[A-Za-z0-9._-]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare blinded cohort review sheets or merge two completed sheets "
            "without trusting copied candidate fields."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser(
        "prepare", description="Create two independent review packages."
    )
    add_candidate_arguments(prepare)
    prepare.add_argument("--output-root", required=True, type=Path)

    merge = subparsers.add_parser(
        "merge", description="Validate and merge two completed review packages."
    )
    add_candidate_arguments(merge)
    merge.add_argument("--reviewer-1-package", required=True, type=Path)
    merge.add_argument("--reviewer-1-sheet", required=True, type=Path)
    merge.add_argument("--reviewer-2-package", required=True, type=Path)
    merge.add_argument("--reviewer-2-sheet", required=True, type=Path)
    merge.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def add_candidate_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--candidate-manifest", required=True, type=Path)
    parser.add_argument("--candidates", type=Path)
    parser.add_argument("--correction-registry", type=Path)


def main() -> int:
    args = parse_args()
    try:
        if args.command == "prepare":
            prepare_review_packages(
                candidate_manifest_path=args.candidate_manifest,
                candidate_path=args.candidates,
                output_root=args.output_root,
                correction_registry_path=args.correction_registry,
            )
        else:
            merge_review_packages(
                candidate_manifest_path=args.candidate_manifest,
                candidate_path=args.candidates,
                reviewer_1_package=args.reviewer_1_package,
                reviewer_1_sheet=args.reviewer_1_sheet,
                reviewer_2_package=args.reviewer_2_package,
                reviewer_2_sheet=args.reviewer_2_sheet,
                output_root=args.output_root,
                correction_registry_path=args.correction_registry,
            )
    except (OSError, ValueError) as error:
        print(f"manage_cohort_reviews: {error}", file=sys.stderr)
        return 1
    return 0


def prepare_review_packages(
    *,
    candidate_manifest_path: Path,
    candidate_path: Path | None,
    output_root: Path,
    correction_registry_path: Path | None = None,
) -> list[dict[str, Any]]:
    source = load_candidate_source(
        candidate_manifest_path, candidate_path, correction_registry_path
    )
    script_path = Path(__file__).resolve()
    output_paths = [output_root / f"reviewer_{slot}" for slot in (1, 2)]
    refuse_existing(output_paths, "review package")

    manifests = []
    for slot in (1, 2):
        package_dir = output_root / f"reviewer_{slot}"
        package_dir.mkdir(parents=True, exist_ok=True)
        sheet_path = package_dir / "cohort_review.csv"
        manifest_path = package_dir / "review_package.json"
        package_identity = package_identity_payload(
            slot=slot,
            source=source,
            prepare_script_sha256=file_sha256(script_path),
        )
        package_id = sha256_json(package_identity)[:16]
        sheet_rows = [
            {
                "review_package_id": package_id,
                "review_slot": str(slot),
                **row,
                **{field: "" for field in REVIEW_INPUT_FIELDS},
            }
            for row in source["review_evidence_rows"]
        ]
        sheet_fields = [
            *PACKAGE_FIELDS,
            *source["review_evidence_fields"],
            *REVIEW_INPUT_FIELDS,
        ]
        write_csv(sheet_path, sheet_rows, sheet_fields)
        review_materials = materialize_review_materials(package_dir, source, slot)
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "package_id": package_id,
            "generated_at": utc_now(),
            "selection_id": source["selection_id"],
            "review_slot": slot,
            "candidate_row_count": len(source["candidate_rows"]),
            "candidate_evidence_sha256": source["candidate_evidence_sha256"],
            "review_evidence_sha256": source["review_evidence_sha256"],
            "fields": {
                "keys": list(KEY_FIELDS),
                "evidence": source["review_evidence_fields"],
                "editable": list(REVIEW_INPUT_FIELDS),
            },
            "tool": file_identity(script_path),
            "dependencies": {
                "cohort_corrections": file_identity(correction_module_path()),
            },
            "inputs": {
                "candidate_manifest": file_identity(candidate_manifest_path),
                "candidate_table": file_identity(source["candidate_path"]),
                "review_protocol": file_identity(source["review_protocol_path"]),
                **optional_file_identity(
                    "correction_registry", correction_registry_path
                ),
            },
            "review_materials": review_materials,
            "prepared_sheet": file_identity(sheet_path),
        }
        write_json(manifest_path, manifest)
        manifests.append(manifest)

    print(
        f"prepared two independent review packages for "
        f"{len(source['candidate_rows'])} candidates "
        f"(selection_id={source['selection_id']})"
    )
    return manifests


def merge_review_packages(
    *,
    candidate_manifest_path: Path,
    candidate_path: Path | None,
    reviewer_1_package: Path,
    reviewer_1_sheet: Path,
    reviewer_2_package: Path,
    reviewer_2_sheet: Path,
    output_root: Path,
    correction_registry_path: Path | None = None,
) -> dict[str, Any]:
    source = load_candidate_source(
        candidate_manifest_path, candidate_path, correction_registry_path
    )
    combined_path = output_root / "tables" / "cohort_reviews.csv"
    manifest_path = output_root / "manifests" / "cohort_reviews.json"
    refuse_existing([combined_path, manifest_path], "merged review")

    review_1 = load_completed_review(
        package_path=reviewer_1_package,
        sheet_path=reviewer_1_sheet,
        expected_slot=1,
        source=source,
    )
    review_2 = load_completed_review(
        package_path=reviewer_2_package,
        sheet_path=reviewer_2_sheet,
        expected_slot=2,
        source=source,
    )
    if review_1["reviewer"].casefold() == review_2["reviewer"].casefold():
        raise ValueError("reviewer 1 and reviewer 2 must be distinct people")

    combined_rows = []
    for evidence in source["review_evidence_rows"]:
        key = row_key(evidence)
        first = review_1["rows_by_key"][key]
        second = review_2["rows_by_key"][key]
        combined_rows.append(
            {
                **evidence,
                **prefixed_review(first, 1),
                **prefixed_review(second, 2),
                "adjudication_decision": "",
                "adjudication_rationale": "",
                "final_family": "",
                "split": "",
            }
        )

    combined_fields = [*source["review_evidence_fields"], *COMBINED_REVIEW_FIELDS]
    combined_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    write_csv(combined_path, combined_rows, combined_fields)
    script_path = Path(__file__).resolve()
    merge_identity = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": source["selection_id"],
        "candidate_manifest_sha256": source["candidate_manifest_sha256"],
        "candidate_table_sha256": source["candidate_table_sha256"],
        "correction_registry_sha256": source["correction_registry_sha256"],
        "correction_module_sha256": source["correction_module_sha256"],
        "reviewer_1_package_sha256": file_sha256(reviewer_1_package),
        "reviewer_1_sheet_sha256": file_sha256(reviewer_1_sheet),
        "reviewer_2_package_sha256": file_sha256(reviewer_2_package),
        "reviewer_2_sheet_sha256": file_sha256(reviewer_2_sheet),
        "merge_script_sha256": file_sha256(script_path),
    }
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "merge_id": sha256_json(merge_identity)[:16],
        "generated_at": utc_now(),
        "selection_id": source["selection_id"],
        "candidate_row_count": len(source["candidate_rows"]),
        "candidate_evidence_sha256": source["candidate_evidence_sha256"],
        "review_evidence_sha256": source["review_evidence_sha256"],
        "tool": file_identity(script_path),
        "dependencies": {
            "cohort_corrections": file_identity(correction_module_path()),
        },
        "inputs": {
            "candidate_manifest": file_identity(candidate_manifest_path),
            "candidate_table": file_identity(source["candidate_path"]),
            **optional_file_identity("correction_registry", correction_registry_path),
            "reviewer_1_package": file_identity(reviewer_1_package),
            "reviewer_1_sheet": file_identity(reviewer_1_sheet),
            "reviewer_2_package": file_identity(reviewer_2_package),
            "reviewer_2_sheet": file_identity(reviewer_2_sheet),
        },
        "reviewers": {
            "reviewer_1": review_1["reviewer"],
            "reviewer_2": review_2["reviewer"],
        },
        "decision_counts": {
            "reviewer_1": review_1["decision_counts"],
            "reviewer_2": review_2["decision_counts"],
        },
        "outputs": {"cohort_reviews": file_identity(combined_path)},
        "adjudication_complete": False,
        "frozen": False,
    }
    write_json(manifest_path, manifest)
    print(
        f"merged {len(combined_rows)} independently reviewed candidates "
        f"(merge_id={manifest['merge_id']}, adjudication_complete=false)"
    )
    return manifest


def load_candidate_source(
    candidate_manifest_path: Path,
    candidate_path: Path | None,
    correction_registry_path: Path | None = None,
) -> dict[str, Any]:
    candidate_manifest = load_json(candidate_manifest_path)
    if candidate_manifest.get("cohort_candidate_schema_version") != (
        CANDIDATE_SCHEMA_VERSION
    ):
        raise ValueError("candidate manifest does not use schema 0.2.0")
    if candidate_manifest.get("frozen") is not False:
        raise ValueError("candidate manifest must describe an unfrozen review pool")
    selection_id = str(candidate_manifest.get("selection_id", "")).strip()
    if not selection_id:
        raise ValueError("candidate selection id is missing")
    candidate_value = str(
        candidate_manifest.get("outputs", {}).get("cohort_candidates", "")
    ).strip()
    if candidate_path is None and not candidate_value:
        raise ValueError("candidate table path is missing")
    resolved_candidate_path = candidate_path or Path(candidate_value)
    candidate_rows, candidate_fields = read_csv(resolved_candidate_path)
    if not candidate_rows:
        raise ValueError("candidate table has no rows")
    if len(candidate_fields) != len(set(candidate_fields)):
        raise ValueError("candidate table has duplicate fields")
    collisions = sorted(
        set(candidate_fields).intersection(
            {
                *PACKAGE_FIELDS,
                *SUPPLEMENTAL_EVIDENCE_FIELDS,
                *REVIEW_INPUT_FIELDS,
                *COMBINED_REVIEW_FIELDS,
            }
        )
    )
    if collisions:
        raise ValueError(
            "candidate table uses reserved review fields: " + ", ".join(collisions)
        )
    index_rows(candidate_rows, "candidate")
    mismatched = [
        row_key(row)
        for row in candidate_rows
        if row.get("selection_id", "").strip() != selection_id
    ]
    if mismatched:
        raise ValueError(
            f"{len(mismatched)} candidate rows do not match selection id {selection_id}"
        )
    correction_by_key = load_registry(
        correction_registry_path,
        selection_id,
        set(index_rows(candidate_rows, "candidate")),
    )
    review_evidence_fields = [*candidate_fields, *SUPPLEMENTAL_EVIDENCE_FIELDS]
    review_evidence_rows = []
    for row in candidate_rows:
        correction = correction_by_key.get(row_key(row), {})
        review_evidence_rows.append(
            {
                **row,
                "proposed_correction_manifest": correction.get("path", ""),
                "proposed_correction_sha256": correction.get("sha256", ""),
            }
        )
    normalized_specs = collect_normalized_specs(
        candidate_rows, resolved_candidate_path.parent
    )
    correction_materials = collect_correction_materials(correction_by_key)
    if not REVIEW_PROTOCOL_PATH.is_file():
        raise ValueError(f"cohort review protocol is missing: {REVIEW_PROTOCOL_PATH}")
    evidence_sha256 = evidence_sha(candidate_rows, candidate_fields)
    review_evidence_sha256 = evidence_sha(review_evidence_rows, review_evidence_fields)
    return {
        "selection_id": selection_id,
        "candidate_path": resolved_candidate_path,
        "candidate_rows": candidate_rows,
        "candidate_fields": candidate_fields,
        "review_evidence_rows": review_evidence_rows,
        "review_evidence_fields": review_evidence_fields,
        "normalized_specs": normalized_specs,
        "correction_materials": correction_materials,
        "review_protocol_path": REVIEW_PROTOCOL_PATH,
        "review_protocol_sha256": file_sha256(REVIEW_PROTOCOL_PATH),
        "candidate_manifest_sha256": file_sha256(candidate_manifest_path),
        "candidate_table_sha256": file_sha256(resolved_candidate_path),
        "candidate_evidence_sha256": evidence_sha256,
        "review_evidence_sha256": review_evidence_sha256,
        "correction_registry_sha256": (
            file_sha256(correction_registry_path) if correction_registry_path else ""
        ),
        "correction_module_sha256": file_sha256(correction_module_path()),
    }


def load_completed_review(
    *,
    package_path: Path,
    sheet_path: Path,
    expected_slot: int,
    source: dict[str, Any],
) -> dict[str, Any]:
    package = load_json(package_path)
    if package.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"reviewer {expected_slot} package uses an unknown schema")
    if package.get("review_slot") != expected_slot:
        raise ValueError(f"reviewer {expected_slot} package has the wrong review slot")
    if package.get("selection_id") != source["selection_id"]:
        raise ValueError(f"reviewer {expected_slot} package selection id changed")
    expected_field_contract = {
        "keys": list(KEY_FIELDS),
        "evidence": source["review_evidence_fields"],
        "editable": list(REVIEW_INPUT_FIELDS),
    }
    if package.get("fields") != expected_field_contract:
        raise ValueError(f"reviewer {expected_slot} package field contract changed")
    package_inputs = package.get("inputs")
    package_tool = package.get("tool")
    if not isinstance(package_inputs, dict) or not isinstance(package_tool, dict):
        raise ValueError(f"reviewer {expected_slot} package identities are malformed")
    expected_hashes = {
        "candidate_manifest_sha256": source["candidate_manifest_sha256"],
        "candidate_table_sha256": source["candidate_table_sha256"],
        "candidate_evidence_sha256": source["candidate_evidence_sha256"],
        "review_evidence_sha256": source["review_evidence_sha256"],
        "review_protocol_sha256": source["review_protocol_sha256"],
        "correction_registry_sha256": source["correction_registry_sha256"],
        "correction_module_sha256": source["correction_module_sha256"],
    }
    for field, expected in expected_hashes.items():
        if package_hash(package, field, expected_slot) != expected:
            raise ValueError(f"reviewer {expected_slot} package {field} changed")

    rows, fields = read_csv(sheet_path)
    expected_fields = [
        *PACKAGE_FIELDS,
        *source["review_evidence_fields"],
        *REVIEW_INPUT_FIELDS,
    ]
    if fields != expected_fields:
        raise ValueError(
            f"reviewer {expected_slot} sheet fields do not match its package"
        )
    if len(rows) != package.get("candidate_row_count"):
        raise ValueError(f"reviewer {expected_slot} sheet row count changed")
    package_id = str(package.get("package_id", "")).strip()
    if not package_id:
        raise ValueError(f"reviewer {expected_slot} package id is missing")
    prepare_script_sha256 = str(package_tool.get("sha256", "")).strip()
    expected_package_id = sha256_json(
        package_identity_payload(
            slot=expected_slot,
            source=source,
            prepare_script_sha256=prepare_script_sha256,
        )
    )[:16]
    if package_id != expected_package_id:
        raise ValueError(f"reviewer {expected_slot} package id is invalid")
    verify_review_materials(
        package=package,
        package_dir=package_path.parent,
        source=source,
        slot=expected_slot,
    )

    rows_by_key = index_rows(rows, f"reviewer {expected_slot}")
    evidence_by_key = index_rows(source["review_evidence_rows"], "candidate")
    if set(rows_by_key) != set(evidence_by_key):
        missing = len(set(evidence_by_key) - set(rows_by_key))
        extra = len(set(rows_by_key) - set(evidence_by_key))
        raise ValueError(
            f"reviewer {expected_slot} sheet candidate keys changed "
            f"(missing={missing}, extra={extra})"
        )

    evidence_rows = []
    reviewers = set()
    decision_counts: Counter[str] = Counter()
    for evidence in source["review_evidence_rows"]:
        key = row_key(evidence)
        row = rows_by_key[key]
        if row.get("review_package_id", "").strip() != package_id:
            raise ValueError(f"reviewer {expected_slot} sheet package id changed")
        if row.get("review_slot", "").strip() != str(expected_slot):
            raise ValueError(f"reviewer {expected_slot} sheet review slot changed")
        changed = [
            field
            for field in source["review_evidence_fields"]
            if row.get(field, "") != evidence.get(field, "")
        ]
        if changed:
            raise ValueError(
                f"reviewer {expected_slot} changed candidate evidence for "
                f"{'/'.join(key)}: {', '.join(changed)}"
            )
        validate_review_input(row, expected_slot, key)
        reviewers.add(row["reviewer"].strip())
        decision_counts[row["decision"].strip()] += 1
        evidence_rows.append(row)

    if len(reviewers) != 1:
        raise ValueError(
            f"reviewer {expected_slot} sheet must use exactly one reviewer identity"
        )
    observed_evidence_sha = evidence_sha(
        evidence_rows, source["review_evidence_fields"]
    )
    if observed_evidence_sha != package.get("review_evidence_sha256"):
        raise ValueError(f"reviewer {expected_slot} sheet evidence hash changed")
    return {
        "reviewer": next(iter(reviewers)),
        "rows_by_key": rows_by_key,
        "decision_counts": dict(sorted(decision_counts.items())),
    }


def collect_normalized_specs(
    rows: list[dict[str, str]], candidate_root: Path
) -> dict[str, dict[str, Any]]:
    specs: dict[str, dict[str, Any]] = {}
    for row in rows:
        accession = safe_path_component(
            row.get("configuration_accession", ""), "configuration accession"
        )
        value = row.get("normalized_spec_path", "").strip()
        declared_sha256 = row.get("normalized_spec_sha256", "").strip()
        if not value or not declared_sha256:
            raise ValueError(f"{accession}: normalized spec identity is incomplete")
        path = Path(value)
        if not path.is_absolute():
            path = candidate_root / path
        path = path.resolve()
        if not path.is_file():
            raise ValueError(f"{accession}: normalized spec is missing: {path}")
        observed_sha256 = file_sha256(path)
        if observed_sha256 != declared_sha256:
            raise ValueError(f"{accession}: normalized spec hash changed")
        observed = {
            "path": path,
            "sha256": observed_sha256,
            "size_bytes": path.stat().st_size,
        }
        prior = specs.get(accession)
        if prior is not None and prior["sha256"] != observed_sha256:
            raise ValueError(
                f"{accession}: candidate rows reference different normalized specs"
            )
        specs.setdefault(accession, observed)
    return dict(sorted(specs.items()))


def collect_correction_materials(
    corrections: dict[tuple[str, str], dict[str, str]],
) -> list[dict[str, Any]]:
    materials = []
    for (family_id, accession), identity in sorted(corrections.items()):
        family_id = safe_path_component(family_id, "correction family id")
        accession = safe_path_component(accession, "correction accession")
        manifest_path = Path(identity["path"]).resolve()
        manifest_name = safe_path_component(
            manifest_path.name, "correction manifest filename"
        )
        correction = load_json(manifest_path)
        artifacts = {}
        occupied = {Path(manifest_name)}
        for field in CORRECTION_ARTIFACT_FIELDS:
            value = correction.get(field)
            if not isinstance(value, dict):
                raise ValueError(
                    f"{family_id}/{accession}: correction {field} is missing"
                )
            raw_path = str(value.get("path", "")).strip()
            relative_path = Path(raw_path)
            if (
                not raw_path
                or relative_path.is_absolute()
                or ".." in relative_path.parts
            ):
                raise ValueError(
                    f"{family_id}/{accession}: correction {field} path is not portable"
                )
            artifact_path, artifact_sha256 = corrections_module.resolve_correction_file(
                correction, field, manifest_path
            )
            if relative_path in occupied:
                raise ValueError(
                    f"{family_id}/{accession}: correction files share a package path"
                )
            occupied.add(relative_path)
            artifacts[field] = {
                "path": artifact_path,
                "relative_path": relative_path,
                "sha256": artifact_sha256,
                "size_bytes": artifact_path.stat().st_size,
            }
        materials.append(
            {
                "family_id": family_id,
                "configuration_accession": accession,
                "manifest": {
                    "path": manifest_path,
                    "relative_path": Path(manifest_name),
                    "sha256": identity["sha256"],
                    "size_bytes": manifest_path.stat().st_size,
                },
                "artifacts": artifacts,
            }
        )
    return materials


def materialize_review_materials(
    package_dir: Path, source: dict[str, Any], slot: int
) -> dict[str, Any]:
    contract = expected_review_materials(source, slot)
    instructions_path = package_dir / contract["instructions"]["path"]
    instructions_path.write_text(review_instructions(slot), encoding="utf-8")
    copy_review_file(
        source["review_protocol_path"],
        package_dir / contract["review_protocol"]["path"],
    )
    for accession, spec in source["normalized_specs"].items():
        destination = package_dir / "evidence" / "specs" / f"{accession}.yaml"
        copy_review_file(spec["path"], destination)
    for correction in source["correction_materials"]:
        root = (
            package_dir
            / "evidence"
            / "corrections"
            / correction["family_id"]
            / correction["configuration_accession"]
        )
        copy_review_file(
            correction["manifest"]["path"],
            root / correction["manifest"]["relative_path"],
        )
        for artifact in correction["artifacts"].values():
            copy_review_file(artifact["path"], root / artifact["relative_path"])
    verify_material_files(package_dir, contract, slot)
    return contract


def expected_review_materials(source: dict[str, Any], slot: int) -> dict[str, Any]:
    instructions = review_instructions(slot).encode("utf-8")
    specs = [
        {
            "configuration_accession": accession,
            "file": portable_identity(
                Path("evidence") / "specs" / f"{accession}.yaml",
                spec["sha256"],
                spec["size_bytes"],
            ),
        }
        for accession, spec in source["normalized_specs"].items()
    ]
    corrections = []
    for correction in source["correction_materials"]:
        root = (
            Path("evidence")
            / "corrections"
            / correction["family_id"]
            / correction["configuration_accession"]
        )
        corrections.append(
            {
                "family_id": correction["family_id"],
                "configuration_accession": correction["configuration_accession"],
                "manifest": portable_identity(
                    root / correction["manifest"]["relative_path"],
                    correction["manifest"]["sha256"],
                    correction["manifest"]["size_bytes"],
                ),
                "artifacts": {
                    field: portable_identity(
                        root / correction["artifacts"][field]["relative_path"],
                        correction["artifacts"][field]["sha256"],
                        correction["artifacts"][field]["size_bytes"],
                    )
                    for field in CORRECTION_ARTIFACT_FIELDS
                },
            }
        )
    return {
        "instructions": portable_identity(
            Path("INSTRUCTIONS.md"),
            hashlib.sha256(instructions).hexdigest(),
            len(instructions),
        ),
        "review_protocol": portable_identity(
            Path("REVIEW_PROTOCOL.md"),
            source["review_protocol_sha256"],
            source["review_protocol_path"].stat().st_size,
        ),
        "normalized_specs": specs,
        "corrections": corrections,
    }


def verify_review_materials(
    *,
    package: dict[str, Any],
    package_dir: Path,
    source: dict[str, Any],
    slot: int,
) -> None:
    expected = expected_review_materials(source, slot)
    if package.get("review_materials") != expected:
        raise ValueError(f"reviewer {slot} package review material contract changed")
    verify_material_files(package_dir, expected, slot)


def verify_material_files(
    package_dir: Path, contract: dict[str, Any], slot: int
) -> None:
    package_root = package_dir.resolve()
    for identity in iter_file_identities(contract):
        relative_path = Path(str(identity["path"]))
        if relative_path.is_absolute() or ".." in relative_path.parts:
            raise ValueError(f"reviewer {slot} package material path is unsafe")
        unresolved_path = package_root / relative_path
        path = unresolved_path.resolve()
        if (
            not path.is_relative_to(package_root)
            or any(
                (package_root / Path(*relative_path.parts[:index])).is_symlink()
                for index in range(1, len(relative_path.parts) + 1)
            )
            or not path.is_file()
        ):
            raise ValueError(
                f"reviewer {slot} package review material is missing: {relative_path}"
            )
        if (
            file_sha256(path) != identity["sha256"]
            or path.stat().st_size != identity["size_bytes"]
        ):
            raise ValueError(
                f"reviewer {slot} package review material hash changed: {relative_path}"
            )


def iter_file_identities(value: Any):
    if isinstance(value, dict):
        if set(value) == {"path", "sha256", "size_bytes"}:
            yield value
            return
        for child in value.values():
            yield from iter_file_identities(child)
    elif isinstance(value, list):
        for child in value:
            yield from iter_file_identities(child)


def portable_identity(path: Path, sha256: str, size_bytes: int) -> dict[str, Any]:
    return {
        "path": path.as_posix(),
        "sha256": sha256,
        "size_bytes": size_bytes,
    }


def copy_review_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def review_instructions(slot: int) -> str:
    return f"""# Cohort review {slot}

Review every row independently. Do not consult the other reviewer's package.

1. Read `REVIEW_PROTOCOL.md` before starting.
2. Open each normalized spec at `evidence/specs/<configuration_accession>.yaml`.
   The absolute `normalized_spec_path` column records provenance only.
3. Compare the assay protocol, portal FASTQs, declared FASTQs, modality, and read
   structure as described in the protocol.
4. When `proposed_correction_manifest` is populated, inspect
   `evidence/corrections/<family_id>/<configuration_accession>/` and compare the
   original spec, corrected spec, and diff before deciding.
5. Edit only `reviewer`, `decision`, `rationale`, `protocol_url`, and `date` in
   `cohort_review.csv`. Use one reviewer name and ISO date (`YYYY-MM-DD`)
   throughout the sheet.
6. Complete every row. Allowed decisions are `include`, `exclude`, and
   `inconclusive`.

Do not edit copied evidence, package identifiers, or files under `evidence/`.
"""


def safe_path_component(value: str, label: str) -> str:
    value = str(value).strip()
    if (
        not value
        or value in {".", ".."}
        or SAFE_PATH_COMPONENT.fullmatch(value) is None
    ):
        raise ValueError(f"{label} is not safe for a review package path: {value}")
    return value


def package_hash(package: dict[str, Any], field: str, slot: int) -> str:
    if field in {
        "candidate_manifest_sha256",
        "candidate_table_sha256",
        "review_protocol_sha256",
        "correction_registry_sha256",
        "correction_module_sha256",
    }:
        input_name = field.removesuffix("_sha256")
        if field == "correction_module_sha256":
            dependencies = package.get("dependencies")
            identity = (
                dependencies.get("cohort_corrections")
                if isinstance(dependencies, dict)
                else None
            )
            if not isinstance(identity, dict):
                raise ValueError(
                    f"reviewer {slot} package correction module identity is malformed"
                )
            return str(identity.get("sha256", ""))
        identity = package["inputs"].get(input_name)
        if not identity and field == "correction_registry_sha256":
            return ""
        if not isinstance(identity, dict):
            raise ValueError(
                f"reviewer {slot} package {input_name} identity is malformed"
            )
        return str(identity.get("sha256", ""))
    return str(package.get(field, ""))


def package_identity_payload(
    *, slot: int, source: dict[str, Any], prepare_script_sha256: str
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "selection_id": source["selection_id"],
        "review_slot": slot,
        "candidate_manifest_sha256": source["candidate_manifest_sha256"],
        "candidate_table_sha256": source["candidate_table_sha256"],
        "candidate_evidence_sha256": source["candidate_evidence_sha256"],
        "review_evidence_sha256": source["review_evidence_sha256"],
        "review_protocol_sha256": source["review_protocol_sha256"],
        "correction_registry_sha256": source["correction_registry_sha256"],
        "correction_module_sha256": source["correction_module_sha256"],
        "prepare_script_sha256": prepare_script_sha256,
    }


def optional_file_identity(name: str, path: Path | None) -> dict[str, Any]:
    return {name: file_identity(path)} if path is not None else {}


def correction_module_path() -> Path:
    value = corrections_module.__file__
    if value is None:
        raise ValueError("could not resolve cohort correction module")
    return Path(value).resolve()


def validate_review_input(row: dict[str, str], slot: int, key: tuple[str, str]) -> None:
    label = f"reviewer {slot} row {'/'.join(key)}"
    if not row.get("reviewer", "").strip():
        raise ValueError(f"{label} is missing reviewer")
    decision = row.get("decision", "").strip()
    if decision not in REVIEW_DECISIONS:
        raise ValueError(f"{label} decision must be include, exclude, or inconclusive")
    for field in ("rationale", "protocol_url"):
        if not row.get(field, "").strip():
            raise ValueError(f"{label} is missing {field}")
    review_date = row.get("date", "").strip()
    if not review_date:
        raise ValueError(f"{label} is missing date")
    if not is_iso_date(review_date):
        raise ValueError(f"{label} date must use YYYY-MM-DD")


def prefixed_review(row: dict[str, str], slot: int) -> dict[str, str]:
    return {
        f"reviewer_{slot}": row["reviewer"].strip(),
        f"reviewer_{slot}_decision": row["decision"].strip(),
        f"reviewer_{slot}_rationale": row["rationale"].strip(),
        f"reviewer_{slot}_protocol_url": row["protocol_url"].strip(),
        f"reviewer_{slot}_date": row["date"].strip(),
    }


def evidence_sha(rows: list[dict[str, str]], fields: list[str]) -> str:
    return sha256_json(
        [{field: row.get(field, "") for field in fields} for row in rows]
    )


def row_key(row: dict[str, str]) -> tuple[str, str]:
    return (
        row.get(KEY_FIELDS[0], "").strip(),
        row.get(KEY_FIELDS[1], "").strip(),
    )


def index_rows(
    rows: list[dict[str, str]], label: str
) -> dict[tuple[str, str], dict[str, str]]:
    indexed = {}
    for row in rows:
        key = row_key(row)
        if not all(key):
            raise ValueError(f"{label} row has a missing candidate key")
        if key in indexed:
            raise ValueError(f"{label} has duplicate candidate key: {'/'.join(key)}")
        indexed[key] = row
    return indexed


def refuse_existing(paths: list[Path], label: str) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise ValueError(
            f"refusing to overwrite existing {label}: {', '.join(existing)}"
        )


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV header is missing: {path}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"CSV header has duplicate fields: {path}")
        return [dict(row) for row in reader], list(reader.fieldnames)


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def is_iso_date(value: str) -> bool:
    try:
        parsed = date.fromisoformat(value)
    except ValueError:
        return False
    return parsed.isoformat() == value


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


def file_identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
