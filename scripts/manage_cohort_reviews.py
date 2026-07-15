#!/usr/bin/env python3
"""Prepare independent cohort review sheets and merge locked reviews."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any


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


def main() -> int:
    args = parse_args()
    try:
        if args.command == "prepare":
            prepare_review_packages(
                candidate_manifest_path=args.candidate_manifest,
                candidate_path=args.candidates,
                output_root=args.output_root,
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
) -> list[dict[str, Any]]:
    source = load_candidate_source(candidate_manifest_path, candidate_path)
    script_path = Path(__file__).resolve()
    output_paths = [
        output_root / f"reviewer_{slot}" / filename
        for slot in (1, 2)
        for filename in ("cohort_review.csv", "review_package.json")
    ]
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
            for row in source["candidate_rows"]
        ]
        sheet_fields = [
            *PACKAGE_FIELDS,
            *source["candidate_fields"],
            *REVIEW_INPUT_FIELDS,
        ]
        write_csv(sheet_path, sheet_rows, sheet_fields)
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "package_id": package_id,
            "generated_at": utc_now(),
            "selection_id": source["selection_id"],
            "review_slot": slot,
            "candidate_row_count": len(source["candidate_rows"]),
            "candidate_evidence_sha256": source["candidate_evidence_sha256"],
            "fields": {
                "keys": list(KEY_FIELDS),
                "evidence": source["candidate_fields"],
                "editable": list(REVIEW_INPUT_FIELDS),
            },
            "tool": file_identity(script_path),
            "inputs": {
                "candidate_manifest": file_identity(candidate_manifest_path),
                "candidate_table": file_identity(source["candidate_path"]),
            },
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
) -> dict[str, Any]:
    source = load_candidate_source(candidate_manifest_path, candidate_path)
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
    for candidate in source["candidate_rows"]:
        key = row_key(candidate)
        first = review_1["rows_by_key"][key]
        second = review_2["rows_by_key"][key]
        combined_rows.append(
            {
                **candidate,
                **prefixed_review(first, 1),
                **prefixed_review(second, 2),
                "adjudication_decision": "",
                "adjudication_rationale": "",
                "final_family": "",
                "split": "",
            }
        )

    combined_fields = [*source["candidate_fields"], *COMBINED_REVIEW_FIELDS]
    combined_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    write_csv(combined_path, combined_rows, combined_fields)
    script_path = Path(__file__).resolve()
    merge_identity = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": source["selection_id"],
        "candidate_manifest_sha256": source["candidate_manifest_sha256"],
        "candidate_table_sha256": source["candidate_table_sha256"],
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
        "tool": file_identity(script_path),
        "inputs": {
            "candidate_manifest": file_identity(candidate_manifest_path),
            "candidate_table": file_identity(source["candidate_path"]),
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
    candidate_manifest_path: Path, candidate_path: Path | None
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
            {*PACKAGE_FIELDS, *REVIEW_INPUT_FIELDS, *COMBINED_REVIEW_FIELDS}
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
    evidence_sha256 = evidence_sha(candidate_rows, candidate_fields)
    return {
        "selection_id": selection_id,
        "candidate_path": resolved_candidate_path,
        "candidate_rows": candidate_rows,
        "candidate_fields": candidate_fields,
        "candidate_manifest_sha256": file_sha256(candidate_manifest_path),
        "candidate_table_sha256": file_sha256(resolved_candidate_path),
        "candidate_evidence_sha256": evidence_sha256,
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
        "evidence": source["candidate_fields"],
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
    }
    for field, expected in expected_hashes.items():
        if package_hash(package, field, expected_slot) != expected:
            raise ValueError(f"reviewer {expected_slot} package {field} changed")

    rows, fields = read_csv(sheet_path)
    expected_fields = [
        *PACKAGE_FIELDS,
        *source["candidate_fields"],
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

    rows_by_key = index_rows(rows, f"reviewer {expected_slot}")
    candidate_by_key = index_rows(source["candidate_rows"], "candidate")
    if set(rows_by_key) != set(candidate_by_key):
        missing = len(set(candidate_by_key) - set(rows_by_key))
        extra = len(set(rows_by_key) - set(candidate_by_key))
        raise ValueError(
            f"reviewer {expected_slot} sheet candidate keys changed "
            f"(missing={missing}, extra={extra})"
        )

    evidence_rows = []
    reviewers = set()
    decision_counts: Counter[str] = Counter()
    for candidate in source["candidate_rows"]:
        key = row_key(candidate)
        row = rows_by_key[key]
        if row.get("review_package_id", "").strip() != package_id:
            raise ValueError(f"reviewer {expected_slot} sheet package id changed")
        if row.get("review_slot", "").strip() != str(expected_slot):
            raise ValueError(f"reviewer {expected_slot} sheet review slot changed")
        changed = [
            field
            for field in source["candidate_fields"]
            if row.get(field, "") != candidate.get(field, "")
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
        evidence_rows, source["candidate_fields"]
    )
    if observed_evidence_sha != package.get("candidate_evidence_sha256"):
        raise ValueError(f"reviewer {expected_slot} sheet evidence hash changed")
    return {
        "reviewer": next(iter(reviewers)),
        "rows_by_key": rows_by_key,
        "decision_counts": dict(sorted(decision_counts.items())),
    }


def package_hash(package: dict[str, Any], field: str, slot: int) -> str:
    if field in {"candidate_manifest_sha256", "candidate_table_sha256"}:
        input_name = field.removesuffix("_sha256")
        identity = package["inputs"].get(input_name)
        if not isinstance(identity, dict):
            raise ValueError(f"reviewer {slot} package {input_name} identity is malformed")
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
        "prepare_script_sha256": prepare_script_sha256,
    }


def validate_review_input(
    row: dict[str, str], slot: int, key: tuple[str, str]
) -> None:
    label = f"reviewer {slot} row {'/'.join(key)}"
    if not row.get("reviewer", "").strip():
        raise ValueError(f"{label} is missing reviewer")
    decision = row.get("decision", "").strip()
    if decision not in REVIEW_DECISIONS:
        raise ValueError(
            f"{label} decision must be include, exclude, or inconclusive"
        )
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
    return sha256_json([{field: row.get(field, "") for field in fields} for row in rows])


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
        raise ValueError(f"refusing to overwrite existing {label}: {', '.join(existing)}")


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
