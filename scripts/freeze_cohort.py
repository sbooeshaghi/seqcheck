#!/usr/bin/env python3
"""Validate two-reviewer cohort decisions and freeze a deterministic split."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

try:
    import cohort_corrections as corrections_module
except ModuleNotFoundError:
    from scripts import cohort_corrections as corrections_module

EFFECTIVE_FIELDS = corrections_module.EFFECTIVE_FIELDS
command_identity = corrections_module.command_identity
effective_candidate_fields = corrections_module.effective_candidate_fields
included_keys = corrections_module.included_keys
load_registry = corrections_module.load_registry
reconcile_evidence = corrections_module.reconcile_evidence
resolve_seqspec_command = corrections_module.resolve_seqspec_command
validate_correction = corrections_module.validate_correction


SCHEMA_VERSION = "0.1.0"
CANDIDATE_SCHEMA_VERSION = "0.2.0"
REVIEW_DECISIONS = {"include", "exclude", "inconclusive"}
FINAL_DECISIONS = {"include", "exclude"}
REVIEW_FIELDS = (
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
BASELINE_REQUIREMENTS = {
    "hydration_status": "normalized",
    "normalized_seqspec_version": "0.5.0",
    "modality_match_status": "matched",
    "structural_check_status": "passed",
    "resource_check_status": "passed",
    "fastq_mapping_status": "matched",
}
VOLATILE_CANDIDATE_FIELDS = {
    "download_attempts",
    "hydration_message",
    "normalized_spec_path",
    "resource_check_attempts",
    "resource_check_message",
    "structural_check_message",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate independent cohort reviews and write a frozen cohort only "
            "when every policy check passes."
        )
    )
    parser.add_argument("--candidate-manifest", required=True, type=Path)
    parser.add_argument("--reviews", required=True, type=Path)
    parser.add_argument("--family-rules", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--candidates", type=Path)
    parser.add_argument("--target-per-family", type=int, default=5)
    parser.add_argument("--calibration-per-family", type=int, default=2)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--correction-registry", type=Path)
    parser.add_argument("--seqspec-bin", type=Path)
    parser.add_argument("--structural-check-timeout-seconds", type=int, default=120)
    parser.add_argument("--resource-check-timeout-seconds", type=int, default=20)
    parser.add_argument("--network-attempts", type=int, default=2)
    parser.add_argument("--retry-backoff-seconds", type=float, default=1.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        return run(args)
    except (OSError, ValueError) as error:
        print(f"freeze_cohort: {error}", file=sys.stderr)
        return 1


def run(args: argparse.Namespace) -> int:
    if args.target_per_family <= 0:
        raise ValueError("target per family must be positive")
    if not 0 < args.calibration_per_family < args.target_per_family:
        raise ValueError(
            "calibration count must be positive and less than the family target"
        )
    if args.structural_check_timeout_seconds <= 0:
        raise ValueError("structural check timeout must be positive")
    if args.resource_check_timeout_seconds <= 0:
        raise ValueError("resource check timeout must be positive")
    if args.network_attempts <= 0:
        raise ValueError("network attempts must be positive")
    if args.retry_backoff_seconds < 0:
        raise ValueError("retry backoff must be nonnegative")

    candidate_manifest = load_json(args.candidate_manifest)
    rules_payload = load_json(args.family_rules)
    if candidate_manifest.get("cohort_candidate_schema_version") != (
        CANDIDATE_SCHEMA_VERSION
    ):
        raise ValueError("candidate manifest does not use schema 0.2.0")
    if candidate_manifest.get("frozen") is not False:
        raise ValueError("candidate manifest must describe an unfrozen review pool")

    candidate_value = str(
        candidate_manifest.get("outputs", {}).get("cohort_candidates", "")
    ).strip()
    if args.candidates is None and not candidate_value:
        raise ValueError("candidate table path is missing")
    candidate_path = args.candidates or Path(candidate_value)
    candidate_rows, candidate_fields = read_csv(candidate_path)
    review_rows, review_fields = read_csv(args.reviews)
    family_ids = parse_family_ids(rules_payload)
    selection_id = str(candidate_manifest.get("selection_id", "")).strip()
    seed = (
        args.seed
        if args.seed is not None
        else int(candidate_manifest.get("rules", {}).get("selection_seed", 0))
    )

    correction_registry = load_registry(
        args.correction_registry,
        selection_id,
        set(strict_index(candidate_rows, "candidate")),
    )
    correction_evidence_errors = reconcile_evidence(
        review_rows, correction_registry
    )
    correction_overlays: dict[tuple[str, str], dict[str, str]] = {}
    correction_validation_errors: dict[tuple[str, str], list[str]] = {}
    seqspec_identity: dict[str, Any] | None = None
    correction_keys_to_validate = included_keys(
        review_rows, correction_registry, effective_review_decision
    )
    if correction_keys_to_validate:
        seqspec_command = resolve_seqspec_command(args.seqspec_bin, candidate_manifest)
        seqspec_identity = command_identity(seqspec_command)
        candidate_by_key = strict_index(candidate_rows, "candidate")
        for key in sorted(correction_keys_to_validate):
            try:
                correction_overlays[key] = validate_correction(
                    candidate=candidate_by_key[key],
                    registry_entry=correction_registry[key],
                    seqspec_command=seqspec_command,
                    structural_timeout_seconds=args.structural_check_timeout_seconds,
                    resource_timeout_seconds=args.resource_check_timeout_seconds,
                    network_attempts=args.network_attempts,
                    retry_backoff_seconds=args.retry_backoff_seconds,
                )
            except (OSError, ValueError, subprocess.SubprocessError) as error:
                correction_validation_errors[key] = [str(error)]

    validation, included = validate_reviews(
        selection_id=selection_id,
        candidate_rows=candidate_rows,
        candidate_fields=candidate_fields,
        review_rows=review_rows,
        review_fields=review_fields,
        family_ids=family_ids,
        target_per_family=args.target_per_family,
        correction_evidence_errors=correction_evidence_errors,
        correction_validation_errors=correction_validation_errors,
        correction_overlays=correction_overlays,
    )
    output_root = args.output_root.resolve()
    validation_dir = output_root / "validation"
    manifest_dir = output_root / "manifests"
    table_dir = output_root / "tables"
    frozen_manifest_path = manifest_dir / "cohort_frozen.json"
    cohort_path = table_dir / "cohort.csv"
    if frozen_manifest_path.exists() or cohort_path.exists():
        raise ValueError("refusing to modify an existing frozen cohort")
    validation_dir.mkdir(parents=True, exist_ok=True)
    validation_path = validation_dir / "cohort_freeze.json"
    validation.update(
        {
            "schema_version": SCHEMA_VERSION,
            "selection_id": selection_id,
            "freeze_script_sha256": file_sha256(Path(__file__).resolve()),
            "correction_module_sha256": correction_module_sha256(),
            "candidate_manifest_sha256": file_sha256(args.candidate_manifest),
            "candidate_table_sha256": file_sha256(candidate_path),
            "reviews_sha256": file_sha256(args.reviews),
            "family_rules_sha256": file_sha256(args.family_rules),
            "correction_registry_sha256": optional_file_sha256(
                args.correction_registry
            ),
            "seqspec": seqspec_identity,
            "correction_count": len(correction_registry),
            "validated_correction_count": len(correction_overlays),
            "target_per_family": args.target_per_family,
            "calibration_per_family": args.calibration_per_family,
            "split_seed": seed,
        }
    )
    write_json(validation_path, validation)
    if not validation["ready_to_freeze"]:
        print(
            f"cohort remains unfrozen: {len(validation['errors'])} validation errors "
            f"({validation_path})",
            file=sys.stderr,
        )
        return 1

    manifest_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    frozen_rows = assign_splits(
        included,
        family_ids,
        args.calibration_per_family,
        seed,
        selection_id,
    )
    output_fields = list(review_fields)
    for field in EFFECTIVE_FIELDS:
        if field not in output_fields:
            output_fields.append(field)
    if "effective_review_decision" not in output_fields:
        output_fields.append("effective_review_decision")
    write_csv(cohort_path, frozen_rows, output_fields)
    freeze_identity = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": selection_id,
        "freeze_script_sha256": file_sha256(Path(__file__).resolve()),
        "correction_module_sha256": correction_module_sha256(),
        "candidate_manifest_sha256": file_sha256(args.candidate_manifest),
        "candidate_table_sha256": file_sha256(candidate_path),
        "reviews_sha256": file_sha256(args.reviews),
        "family_rules_sha256": file_sha256(args.family_rules),
        "correction_registry_sha256": optional_file_sha256(
            args.correction_registry
        ),
        "seqspec": seqspec_identity,
        "target_per_family": args.target_per_family,
        "calibration_per_family": args.calibration_per_family,
        "split_seed": seed,
        "assignments": [
            {
                "family_id": row["final_family"],
                "configuration_accession": row["configuration_accession"],
                "split": row["split"],
                "correction_manifest_sha256": row["correction_manifest_sha256"],
                "effective_spec_sha256": row["effective_spec_sha256"],
                "effective_deduplication_key": row[
                    "effective_deduplication_key"
                ],
            }
            for row in frozen_rows
        ],
    }
    freeze_id = sha256_json(freeze_identity)[:16]
    frozen_manifest = {
        "schema_version": SCHEMA_VERSION,
        "freeze_id": freeze_id,
        "selection_id": selection_id,
        "generated_at": utc_now(),
        "frozen": True,
        "split_assigned": True,
        "tool": {
            "script": file_identity(Path(__file__).resolve()),
            "dependencies": {
                "cohort_corrections": file_identity(correction_module_path()),
            },
        },
        "policy": {
            "target_per_family": args.target_per_family,
            "calibration_per_family": args.calibration_per_family,
            "evaluation_per_family": (
                args.target_per_family - args.calibration_per_family
            ),
            "split_seed": seed,
            "required_independent_reviewers": 2,
            "corrected_baseline_requires_unanimous_review": True,
            "baseline_requirements": BASELINE_REQUIREMENTS,
        },
        "inputs": {
            "candidate_manifest": file_identity(args.candidate_manifest),
            "candidate_table": file_identity(candidate_path),
            "reviews": file_identity(args.reviews),
            "family_rules": file_identity(args.family_rules),
            **optional_file_identity(
                "correction_registry", args.correction_registry
            ),
        },
        "outputs": {
            "cohort": file_identity(cohort_path),
            "validation": file_identity(validation_path),
        },
        "corrections": [
            {
                "family_id": row["final_family"],
                "configuration_accession": row["configuration_accession"],
                "manifest": file_identity(Path(row["correction_manifest"])),
                "effective_spec": file_identity(Path(row["effective_spec_path"])),
                "effective_structure_sha256": row["effective_structure_sha256"],
                "effective_deduplication_key": row[
                    "effective_deduplication_key"
                ],
            }
            for row in frozen_rows
            if row.get("correction_applied") == "true"
        ],
        "counts": {
            "total": len(frozen_rows),
            "calibration": sum(row["split"] == "calibration" for row in frozen_rows),
            "evaluation": sum(row["split"] == "evaluation" for row in frozen_rows),
            "corrected": sum(
                row.get("correction_applied") == "true" for row in frozen_rows
            ),
            "families": dict(
                sorted(Counter(row["final_family"] for row in frozen_rows).items())
            ),
        },
    }
    write_json(frozen_manifest_path, frozen_manifest)
    print(
        f"froze {len(frozen_rows)} configurations "
        f"(freeze_id={freeze_id}, calibration={frozen_manifest['counts']['calibration']}, "
        f"evaluation={frozen_manifest['counts']['evaluation']})"
    )
    return 0


def validate_reviews(
    *,
    selection_id: str,
    candidate_rows: list[dict[str, str]],
    candidate_fields: list[str],
    review_rows: list[dict[str, str]],
    review_fields: list[str],
    family_ids: list[str],
    target_per_family: int,
    correction_evidence_errors: dict[tuple[str, str], list[str]] | None = None,
    correction_validation_errors: dict[tuple[str, str], list[str]] | None = None,
    correction_overlays: dict[tuple[str, str], dict[str, str]] | None = None,
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    correction_evidence_errors = correction_evidence_errors or {}
    correction_validation_errors = correction_validation_errors or {}
    correction_overlays = correction_overlays or {}
    errors: list[str] = []
    row_errors: list[dict[str, Any]] = []
    warnings: list[str] = []
    key_fields = ("family_id", "configuration_accession")
    missing_review_fields = sorted(set(REVIEW_FIELDS) - set(review_fields))
    if missing_review_fields:
        errors.append(
            "review table is missing fields: " + ", ".join(missing_review_fields)
        )
    if not selection_id:
        errors.append("candidate selection id is missing")

    candidates_by_key = index_rows(candidate_rows, key_fields, "candidate", errors)
    reviews_by_key = index_rows(review_rows, key_fields, "review", errors)
    missing_reviews = sorted(set(candidates_by_key) - set(reviews_by_key))
    extra_reviews = sorted(set(reviews_by_key) - set(candidates_by_key))
    if missing_reviews:
        errors.append(f"review table is missing {len(missing_reviews)} candidate rows")
    if extra_reviews:
        errors.append(f"review table has {len(extra_reviews)} rows not in candidates")

    included: list[dict[str, str]] = []
    effective_counts: Counter[str] = Counter()
    for key in sorted(set(candidates_by_key).intersection(reviews_by_key)):
        candidate = candidates_by_key[key]
        review = reviews_by_key[key]
        problems = compare_candidate_fields(candidate, review, candidate_fields)
        problems.extend(correction_evidence_errors.get(key, []))
        if candidate.get("selection_id", "") != selection_id:
            problems.append("candidate selection_id does not match manifest")
        if review.get("selection_id", "") != selection_id:
            problems.append("review selection_id does not match manifest")
        decision, decision_errors = effective_review_decision(review)
        problems.extend(decision_errors)
        if review.get("split", "").strip():
            problems.append(
                "split must remain blank until the freeze command assigns it"
            )
        if decision:
            effective_counts[decision] += 1
        if decision == "include":
            if review.get("final_family", "").strip() != candidate.get("family_id", ""):
                problems.append(
                    "included row final_family must equal its proposed family"
                )
            correction_expected = (
                key in correction_overlays or key in correction_validation_errors
            )
            if correction_expected and any(
                review.get(f"reviewer_{number}_decision", "").strip() != "include"
                for number in (1, 2)
            ):
                problems.append(
                    "corrected baseline requires unanimous include decisions from both reviewers"
                )
            problems.extend(correction_validation_errors.get(key, []))
            effective = correction_overlays.get(key)
            if correction_expected and effective is None:
                problems.append("proposed correction did not pass freeze-time validation")
            if effective is None:
                effective = effective_candidate_fields(candidate)
            for field, expected in BASELINE_REQUIREMENTS.items():
                effective_field = f"effective_{field}"
                if effective.get(effective_field, "") != expected:
                    problems.append(
                        f"included row requires {field}={expected}, observed "
                        f"{effective.get(effective_field, '') or '<missing>'}"
                    )
            if not effective.get("effective_deduplication_key", "").strip():
                problems.append("included row is missing a deduplication key")
            if not problems:
                included.append(
                    {
                        **review,
                        **effective,
                        "effective_review_decision": decision,
                    }
                )
        elif decision == "exclude" and review.get("final_family", "").strip():
            problems.append("excluded row must not have final_family")
        if problems:
            row_errors.append(
                {
                    "family_id": key[0],
                    "configuration_accession": key[1],
                    "errors": problems,
                }
            )

    if row_errors:
        errors.append(f"{len(row_errors)} review rows failed validation")
    included_by_family = Counter(row["final_family"] for row in included)
    for family_id in family_ids:
        observed = included_by_family[family_id]
        if observed != target_per_family:
            errors.append(
                f"family {family_id} requires {target_per_family} included rows, "
                f"observed {observed}"
            )
    unknown_families = sorted(set(included_by_family) - set(family_ids))
    if unknown_families:
        errors.append(
            "included rows use unknown families: " + ", ".join(unknown_families)
        )
    expected_total = target_per_family * len(family_ids)
    if len(included) != expected_total:
        errors.append(
            f"cohort requires {expected_total} included rows, observed {len(included)}"
        )

    duplicate_accessions = duplicate_values(
        row["configuration_accession"] for row in included
    )
    if duplicate_accessions:
        errors.append(
            "included configurations are not unique: " + ", ".join(duplicate_accessions)
        )
    duplicate_structures = duplicate_values(
        row["effective_deduplication_key"] for row in included
    )
    if duplicate_structures:
        errors.append(
            f"included cohort has {len(duplicate_structures)} duplicate structure/FASTQ keys"
        )

    for family_id in family_ids:
        labs = {
            row.get("lab", "")
            for row in included
            if row["final_family"] == family_id and row.get("lab", "")
        }
        if included_by_family[family_id] and len(labs) < 2:
            warnings.append(
                f"family {family_id} has {len(labs)} submitting laboratory in the final cohort"
            )

    return (
        {
            "ready_to_freeze": not errors,
            "candidate_row_count": len(candidate_rows),
            "review_row_count": len(review_rows),
            "included_row_count": len(included),
            "included_correction_count": sum(
                row.get("correction_applied") == "true" for row in included
            ),
            "included_by_family": dict(sorted(included_by_family.items())),
            "effective_decision_counts": dict(sorted(effective_counts.items())),
            "errors": errors,
            "row_errors": row_errors,
            "warnings": warnings,
        },
        included,
    )


def strict_index(
    rows: list[dict[str, str]], label: str
) -> dict[tuple[str, str], dict[str, str]]:
    errors: list[str] = []
    indexed = index_rows(rows, ("family_id", "configuration_accession"), label, errors)
    if errors:
        raise ValueError("; ".join(errors))
    return indexed


def row_key(row: dict[str, str]) -> tuple[str, str]:
    return (
        row.get("family_id", "").strip(),
        row.get("configuration_accession", "").strip(),
    )


def effective_review_decision(row: dict[str, str]) -> tuple[str, list[str]]:
    errors = []
    reviewer_1 = row.get("reviewer_1", "").strip()
    reviewer_2 = row.get("reviewer_2", "").strip()
    if not reviewer_1:
        errors.append("reviewer_1 is required")
    if not reviewer_2:
        errors.append("reviewer_2 is required")
    if reviewer_1 and reviewer_1 == reviewer_2:
        errors.append("reviewers must be distinct")
    decisions = []
    for number in (1, 2):
        decision = row.get(f"reviewer_{number}_decision", "").strip()
        decisions.append(decision)
        if decision not in REVIEW_DECISIONS:
            errors.append(
                f"reviewer_{number}_decision must be include, exclude, or inconclusive"
            )
        for suffix in ("rationale", "protocol_url"):
            if not row.get(f"reviewer_{number}_{suffix}", "").strip():
                errors.append(f"reviewer_{number}_{suffix} is required")
        review_date = row.get(f"reviewer_{number}_date", "").strip()
        if not review_date:
            errors.append(f"reviewer_{number}_date is required")
        elif not is_iso_date(review_date):
            errors.append(f"reviewer_{number}_date must use YYYY-MM-DD")

    adjudication = row.get("adjudication_decision", "").strip()
    needs_adjudication = len(decisions) == 2 and (
        decisions[0] != decisions[1] or "inconclusive" in decisions
    )
    if needs_adjudication and adjudication not in FINAL_DECISIONS:
        errors.append(
            "review disagreement or inconclusive decision requires adjudication"
        )
    if adjudication:
        if not needs_adjudication:
            errors.append(
                "adjudication must remain blank when reviewer decisions agree"
            )
        if adjudication not in FINAL_DECISIONS:
            errors.append("adjudication_decision must be include or exclude")
        if not row.get("adjudication_rationale", "").strip():
            errors.append("adjudication_rationale is required when adjudication is set")

    if errors:
        return "", errors
    if adjudication:
        return adjudication, []
    if decisions[0] == decisions[1] and decisions[0] in FINAL_DECISIONS:
        return decisions[0], []
    return "", ["reviews do not resolve to a final include or exclude decision"]


def assign_splits(
    included: list[dict[str, str]],
    family_ids: list[str],
    calibration_per_family: int,
    seed: int,
    selection_id: str,
) -> list[dict[str, str]]:
    assigned = []
    for family_id in family_ids:
        family_rows = sorted(
            (row for row in included if row["final_family"] == family_id),
            key=lambda row: stable_order_key(
                seed,
                selection_id,
                family_id,
                row["configuration_accession"],
            ),
        )
        for index, row in enumerate(family_rows):
            assigned.append(
                {
                    **row,
                    "split": (
                        "calibration"
                        if index < calibration_per_family
                        else "evaluation"
                    ),
                }
            )
    return assigned


def compare_candidate_fields(
    candidate: dict[str, str], review: dict[str, str], candidate_fields: list[str]
) -> list[str]:
    changed = [
        field
        for field in candidate_fields
        if field not in VOLATILE_CANDIDATE_FIELDS
        if review.get(field, "") != candidate.get(field, "")
    ]
    if not changed:
        return []
    return ["candidate fields changed in review table: " + ", ".join(changed)]


def index_rows(
    rows: list[dict[str, str]],
    key_fields: tuple[str, ...],
    label: str,
    errors: list[str],
) -> dict[tuple[str, ...], dict[str, str]]:
    indexed = {}
    for row in rows:
        key = tuple(row.get(field, "").strip() for field in key_fields)
        if not all(key):
            errors.append(f"{label} row has a missing key field")
            continue
        if key in indexed:
            errors.append(f"{label} table has duplicate key: {'/'.join(key)}")
            continue
        indexed[key] = row
    return indexed


def duplicate_values(values: Any) -> list[str]:
    counts = Counter(value for value in values if value)
    return sorted(value for value, count in counts.items() if count > 1)


def parse_family_ids(payload: dict[str, Any]) -> list[str]:
    family_ids = [
        str(item.get("id", "")).strip()
        for item in payload.get("families", [])
        if isinstance(item, dict) and str(item.get("id", "")).strip()
    ]
    if not family_ids or len(family_ids) != len(set(family_ids)):
        raise ValueError("family rules need unique family ids")
    return family_ids


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV header is missing: {path}")
        rows = [dict(row) for row in reader]
        return rows, list(reader.fieldnames)


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


def stable_order_key(seed: int, *values: str) -> str:
    return sha256_json([seed, *values])


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


def optional_file_sha256(path: Path | None) -> str:
    return file_sha256(path) if path is not None else ""


def optional_file_identity(name: str, path: Path | None) -> dict[str, Any]:
    return {name: file_identity(path)} if path is not None else {}


def correction_module_path() -> Path:
    value = corrections_module.__file__
    if value is None:
        raise ValueError("could not resolve cohort correction module")
    return Path(value).resolve()


def correction_module_sha256() -> str:
    return file_sha256(correction_module_path())


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
