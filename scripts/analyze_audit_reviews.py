#!/usr/bin/env python3
"""Analyze completed, adjudicated reviews of IGVF audit findings."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from statistics import NormalDist
from typing import Any

try:
    import manage_audit_reviews as reviews
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import manage_audit_reviews as reviews
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
TOOL_VERSION = "0.1.0"
SUBMITTER_CONTACT_VALUES = {"yes", "no", "not_attempted"}
ADJUDICATION_EDITABLE_FIELDS = (
    "adjudicator",
    "adjudication_date",
    "final_classification",
    "adjudication_rationale",
    "submitter_contacted",
    "submitter_response",
    "proposed_correction",
    "final_status_note",
)
FINAL_FIELDS = (
    "final_classification_source",
    "confirmed_problem",
    "inconclusive",
    "pass_control_consistent",
)
PRECISION_FIELDS = (
    "analysis_id",
    "dimension",
    "group_id",
    "candidate_cases",
    "confirmed_cases",
    "inconclusive_cases",
    "conservative_denominator",
    "conservative_precision",
    "conservative_ci_lower",
    "conservative_ci_upper",
    "evaluable_denominator",
    "evaluable_precision",
    "evaluable_ci_lower",
    "evaluable_ci_upper",
)
AGREEMENT_FIELDS = (
    "analysis_id",
    "source_kind",
    "reviewed_cases",
    "exact_agreements",
    "exact_agreement_fraction",
    "cohen_kappa",
)
ENDPOINT_FIELDS = (
    "analysis_id",
    "endpoint",
    "value",
    "target",
    "comparison",
    "target_met",
    "endpoint_class",
    "description",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze completed adjudication of blinded audit reviews."
    )
    parser.add_argument("--review-manifest", required=True, type=Path)
    parser.add_argument("--adjudicated-reviews", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        path = analyze_reviews(
            review_manifest_path=args.review_manifest.resolve(),
            adjudicated_reviews_path=args.adjudicated_reviews.resolve(),
            protocol_path=args.protocol.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"analyze_audit_reviews: {error}", file=sys.stderr)
        return 1
    print(path)
    return 0


def analyze_reviews(
    *,
    review_manifest_path: Path,
    adjudicated_reviews_path: Path,
    protocol_path: Path,
    output_root: Path,
) -> Path:
    refuse_output_root(output_root)
    manifest, source_rows, source_fields = load_review_source(
        review_manifest_path, protocol_path
    )
    protocol = reviews.load_protocol(protocol_path)
    completed_rows, completed_fields = read_csv(adjudicated_reviews_path)
    if completed_fields != source_fields:
        raise ValueError("adjudicated audit review fields changed")
    outcomes = validate_adjudication(
        source_rows=source_rows,
        completed_rows=completed_rows,
        manifest=manifest,
        protocol=protocol,
    )
    tool = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    stable = {
        "schema_version": SCHEMA_VERSION,
        "merge_id": manifest["merge_id"],
        "selection_id": manifest["selection_id"],
        "tool": runtime.functional_script_identity(tool),
        "inputs": {
            "review_manifest": runtime.file_sha256(review_manifest_path),
            "adjudicated_reviews": runtime.file_sha256(adjudicated_reviews_path),
            "protocol": runtime.file_sha256(protocol_path),
        },
    }
    analysis_id = runtime.sha256_json(stable)[:16]
    precision = build_precision_summary(analysis_id, outcomes, protocol)
    agreement = build_agreement_summary(analysis_id, outcomes)
    endpoints = build_endpoints(
        analysis_id=analysis_id,
        outcomes=outcomes,
        agreement=agreement,
        protocol=protocol,
    )
    output_fields = [*source_fields, *FINAL_FIELDS]
    outcome_path = output_root / "tables" / "adjudicated_outcomes.csv"
    precision_path = output_root / "tables" / "confirmation_precision.csv"
    agreement_path = output_root / "tables" / "reviewer_agreement.csv"
    endpoint_path = output_root / "tables" / "review_endpoints.csv"
    runtime.write_csv(outcome_path, outcomes, output_fields)
    runtime.write_csv(precision_path, precision, list(PRECISION_FIELDS))
    runtime.write_csv(agreement_path, agreement, list(AGREEMENT_FIELDS))
    runtime.write_csv(endpoint_path, endpoints, list(ENDPOINT_FIELDS))

    candidate_count = sum(row["source_kind"] == "candidate" for row in outcomes)
    control_count = sum(row["source_kind"] == "pass_control" for row in outcomes)
    selection = protocol["selection"]
    checks = {
        "review_row_count_reconciles": len(outcomes) == manifest["review_rows"],
        "all_rows_have_final_classification": all(
            row["final_classification"] in reviews.CLASSIFICATIONS for row in outcomes
        ),
        "adjudications_reconcile": sum(
            row["final_classification_source"] == "human_adjudication"
            for row in outcomes
        )
        == manifest["adjudication_required"],
        "candidate_count_complete": candidate_count
        >= selection["minimum_candidate_findings"],
        "pass_control_count_complete": control_count
        == selection["pass_control_findings"],
        "source_classes_reconcile": candidate_count + control_count == len(outcomes),
    }
    review_target_rows = [
        row for row in endpoints if row["endpoint_class"] == "review_scientific_target"
    ]
    validation = {
        **stable,
        "analysis_id": analysis_id,
        "generated_at": utc_now(),
        "valid": all(checks.values()),
        "review_scientific_targets_met": all(
            row["target_met"] for row in review_target_rows
        ),
        "downstream_targets_status": "not_evaluated_by_review_analysis",
        "checks": checks,
        "counts": {
            "reviewed_cases": len(outcomes),
            "candidate_cases": candidate_count,
            "pass_control_cases": control_count,
            "high_confidence_candidate_cases": sum(
                row["source_kind"] == "candidate"
                and parse_bool(row["high_confidence"], "high confidence")
                for row in outcomes
            ),
            "confirmed_candidate_cases": sum(
                row["source_kind"] == "candidate" and row["confirmed_problem"]
                for row in outcomes
            ),
            "inconclusive_candidate_cases": sum(
                row["source_kind"] == "candidate" and row["inconclusive"]
                for row in outcomes
            ),
            "human_adjudications": sum(
                row["final_classification_source"] == "human_adjudication"
                for row in outcomes
            ),
        },
        "outputs": {
            "outcomes": runtime.file_identity(outcome_path),
            "precision": runtime.file_identity(precision_path),
            "agreement": runtime.file_identity(agreement_path),
            "endpoints": runtime.file_identity(endpoint_path),
        },
        "tools": {"analyzer": tool},
    }
    validation_path = output_root / "validation" / "audit_reviews.json"
    runtime.write_json(validation_path, validation)
    if not validation["valid"]:
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise ValueError(f"audit review analysis validation failed: {failed}")
    analysis_manifest_path = output_root / "manifests" / "audit_review_analysis.json"
    analysis_manifest = {
        **stable,
        "analysis_id": analysis_id,
        "generated_at": utc_now(),
        "valid": True,
        "review_scientific_targets_met": validation["review_scientific_targets_met"],
        "inputs": {
            "review_manifest": runtime.file_identity(review_manifest_path),
            "adjudicated_reviews": runtime.file_identity(adjudicated_reviews_path),
            "protocol": runtime.file_identity(protocol_path),
        },
        "outputs": {
            **validation["outputs"],
            "validation": runtime.file_identity(validation_path),
        },
        "tools": {"analyzer": tool},
        "manifest_path": str(analysis_manifest_path),
    }
    runtime.write_json(analysis_manifest_path, analysis_manifest)
    return validation_path


def load_review_source(
    manifest_path: Path, protocol_path: Path
) -> tuple[dict[str, Any], list[dict[str, str]], list[str]]:
    manifest = runtime.load_json(manifest_path)
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("audit review manifest schema is unsupported")
    stable = {
        key: manifest.get(key)
        for key in (
            "schema_version",
            "selection_id",
            "tool",
            "reviewer_1_package_sha256",
            "reviewer_1_sheet_sha256",
            "reviewer_2_package_sha256",
            "reviewer_2_sheet_sha256",
        )
    }
    if runtime.sha256_json(stable)[:16] != manifest.get("merge_id"):
        raise ValueError("audit review merge id is not content-addressed")
    protocol_identity = manifest.get("inputs", {}).get("protocol", {})
    verify_identity(protocol_path, protocol_identity, "audit adjudication protocol")
    output_identity = manifest.get("outputs", {}).get("audit_reviews", {})
    source_path = identity_path(output_identity, "merged audit reviews")
    verify_identity(source_path, output_identity, "merged audit reviews")
    rows, fields = read_csv(source_path)
    if not rows or len(rows) != manifest.get("review_rows"):
        raise ValueError("merged audit review row count changed")
    if manifest.get("adjudication_complete") not in {True, False}:
        raise ValueError("audit review adjudication status is invalid")
    return manifest, rows, fields


def validate_adjudication(
    *,
    source_rows: list[dict[str, str]],
    completed_rows: list[dict[str, str]],
    manifest: dict[str, Any],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    source = index_unique(source_rows, "case_id", "merged audit reviews")
    completed = index_unique(completed_rows, "case_id", "adjudicated audit reviews")
    if set(source) != set(completed):
        raise ValueError("adjudicated audit review case ids changed")
    immutable_fields = [
        field for field in source_rows[0] if field not in ADJUDICATION_EDITABLE_FIELDS
    ]
    confirmed = set(protocol["analysis"]["confirmed_classifications"])
    inconclusive = protocol["analysis"]["inconclusive_classification"]
    pass_consistent = set(protocol["analysis"]["pass_consistent_classifications"])
    outcomes = []
    adjudications = 0
    for case_id in sorted(source):
        original = source[case_id]
        row = completed[case_id]
        changed = [field for field in immutable_fields if row[field] != original[field]]
        if changed:
            raise ValueError(
                f"{case_id}: adjudication changed review evidence: {', '.join(changed)}"
            )
        left = row["reviewer_1_classification"]
        right = row["reviewer_2_classification"]
        if left not in reviews.CLASSIFICATIONS or right not in reviews.CLASSIFICATIONS:
            raise ValueError(f"{case_id}: reviewer classification is invalid")
        agreement = parse_bool(
            row["exact_classification_agreement"], "classification agreement"
        )
        needs_adjudication = parse_bool(
            row["needs_adjudication"], "adjudication requirement"
        )
        if agreement == needs_adjudication or agreement != (left == right):
            raise ValueError(f"{case_id}: agreement fields do not reconcile")
        final = row["final_classification"].strip()
        if final not in reviews.CLASSIFICATIONS:
            raise ValueError(f"{case_id}: final classification is invalid")
        if agreement:
            if (
                final != left
                or row["consensus_classification"] != left
                or any(
                    row[field].strip()
                    for field in (
                        "adjudicator",
                        "adjudication_date",
                        "adjudication_rationale",
                    )
                )
            ):
                raise ValueError(f"{case_id}: consensus classification was overridden")
            source_label = "reviewer_consensus"
        else:
            adjudicator = row["adjudicator"].strip()
            rationale = row["adjudication_rationale"].strip()
            if not adjudicator or not rationale:
                raise ValueError(f"{case_id}: adjudication is incomplete")
            if adjudicator.casefold() in {
                row["reviewer_1"].casefold(),
                row["reviewer_2"].casefold(),
            }:
                raise ValueError(f"{case_id}: adjudicator must be a third person")
            validate_date(row["adjudication_date"], "adjudication")
            source_label = "human_adjudication"
            adjudications += 1
        contacted = row["submitter_contacted"].strip()
        if contacted not in SUBMITTER_CONTACT_VALUES:
            raise ValueError(f"{case_id}: submitter contact status is invalid")
        if contacted == "yes" and not row["submitter_response"].strip():
            raise ValueError(f"{case_id}: submitter response is required")
        source_kind = row["source_kind"]
        if source_kind not in {"candidate", "pass_control"}:
            raise ValueError(f"{case_id}: hidden source class is invalid")
        outcomes.append(
            {
                **row,
                "final_classification": final,
                "final_classification_source": source_label,
                "confirmed_problem": final in confirmed,
                "inconclusive": final == inconclusive,
                "pass_control_consistent": source_kind == "pass_control"
                and final in pass_consistent,
            }
        )
    if adjudications != manifest.get("adjudication_required"):
        raise ValueError("completed adjudication count differs from review manifest")
    return outcomes


def build_precision_summary(
    analysis_id: str,
    outcomes: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    candidates = [row for row in outcomes if row["source_kind"] == "candidate"]
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in candidates:
        memberships = [
            ("overall", "all"),
            ("high_confidence", row["high_confidence"]),
            ("assessment_type", row["assessment_type"]),
            ("assessment_code", row["assessment_code"]),
            ("assay_family", row["matched_assay_family"]),
            ("access_class", row["access_class"]),
        ]
        for membership in memberships:
            groups[membership].append(row)
    confidence = protocol["analysis"]["confidence_level"]
    rows = []
    for (dimension, group_id), values in sorted(groups.items()):
        confirmed = sum(row["confirmed_problem"] for row in values)
        inconclusive = sum(row["inconclusive"] for row in values)
        conservative_n = len(values)
        evaluable_n = conservative_n - inconclusive
        conservative_interval = wilson_interval(confirmed, conservative_n, confidence)
        evaluable_interval = wilson_interval(confirmed, evaluable_n, confidence)
        rows.append(
            {
                "analysis_id": analysis_id,
                "dimension": dimension,
                "group_id": group_id,
                "candidate_cases": len(values),
                "confirmed_cases": confirmed,
                "inconclusive_cases": inconclusive,
                "conservative_denominator": conservative_n,
                "conservative_precision": confirmed / conservative_n,
                "conservative_ci_lower": conservative_interval[0],
                "conservative_ci_upper": conservative_interval[1],
                "evaluable_denominator": evaluable_n,
                "evaluable_precision": (
                    confirmed / evaluable_n if evaluable_n else None
                ),
                "evaluable_ci_lower": evaluable_interval[0],
                "evaluable_ci_upper": evaluable_interval[1],
            }
        )
    return rows


def build_agreement_summary(
    analysis_id: str, outcomes: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    groups = {"all": outcomes}
    groups.update(
        {
            kind: [row for row in outcomes if row["source_kind"] == kind]
            for kind in ("candidate", "pass_control")
        }
    )
    result = []
    for source_kind, rows in groups.items():
        exact = sum(
            row["reviewer_1_classification"] == row["reviewer_2_classification"]
            for row in rows
        )
        result.append(
            {
                "analysis_id": analysis_id,
                "source_kind": source_kind,
                "reviewed_cases": len(rows),
                "exact_agreements": exact,
                "exact_agreement_fraction": exact / len(rows) if rows else None,
                "cohen_kappa": cohen_kappa(rows),
            }
        )
    return result


def build_endpoints(
    *,
    analysis_id: str,
    outcomes: list[dict[str, Any]],
    agreement: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    targets = protocol["scientific_targets"]
    high = [
        row
        for row in outcomes
        if row["source_kind"] == "candidate"
        and parse_bool(row["high_confidence"], "high confidence")
    ]
    high_precision = (
        sum(row["confirmed_problem"] for row in high) / len(high) if high else None
    )
    overall_kappa = next(
        row["cohen_kappa"] for row in agreement if row["source_kind"] == "all"
    )
    controls = [row for row in outcomes if row["source_kind"] == "pass_control"]
    control_consistency = (
        sum(row["pass_control_consistent"] for row in controls) / len(controls)
        if controls
        else None
    )
    values = [
        (
            "high_confidence_conservative_confirmation_precision",
            high_precision,
            targets["minimum_high_confidence_confirmation_precision"],
            "At least",
            "Conservative precision among policy-supported candidates.",
        ),
        (
            "unweighted_cohen_kappa",
            overall_kappa,
            targets["minimum_cohen_kappa"],
            "At least",
            "Agreement between the two blinded reviewers.",
        ),
        (
            "pass_control_consistency",
            control_consistency,
            targets["minimum_pass_control_consistency"],
            "At least",
            "Pass controls adjudicated as intended assay design.",
        ),
    ]
    return [
        {
            "analysis_id": analysis_id,
            "endpoint": endpoint,
            "value": value,
            "target": target,
            "comparison": comparison,
            "target_met": value is not None and value >= target,
            "endpoint_class": "review_scientific_target",
            "description": description,
        }
        for endpoint, value, target, comparison, description in values
    ]


def cohen_kappa(rows: list[dict[str, Any]]) -> float | None:
    if not rows:
        return None
    left = Counter(row["reviewer_1_classification"] for row in rows)
    right = Counter(row["reviewer_2_classification"] for row in rows)
    observed = sum(
        row["reviewer_1_classification"] == row["reviewer_2_classification"]
        for row in rows
    ) / len(rows)
    expected = sum(left[label] * right[label] for label in reviews.CLASSIFICATIONS) / (
        len(rows) ** 2
    )
    if math.isclose(expected, 1.0):
        return None
    return (observed - expected) / (1 - expected)


def wilson_interval(
    successes: int, total: int, confidence: float
) -> tuple[float | None, float | None]:
    if total <= 0:
        return None, None
    z = NormalDist().inv_cdf(0.5 + confidence / 2)
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    margin = (
        z
        * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total))
        / denominator
    )
    return max(0.0, center - margin), min(1.0, center + margin)


def verify_identity(path: Path, identity: Any, label: str) -> None:
    if not isinstance(identity, dict):
        raise ValueError(f"{label} identity is missing")
    if runtime.file_sha256(path) != identity.get("sha256"):
        raise ValueError(f"{label} hash changed")


def identity_path(identity: Any, label: str) -> Path:
    if not isinstance(identity, dict) or not str(identity.get("path", "")).strip():
        raise ValueError(f"{label} path is missing")
    return Path(str(identity["path"])).resolve()


def index_unique(
    rows: list[dict[str, str]], field: str, label: str
) -> dict[str, dict[str, str]]:
    result = {}
    for row in rows:
        key = row.get(field, "")
        if not key or key in result:
            raise ValueError(f"{label} {field} values are incomplete or duplicated")
        result[key] = row
    return result


def parse_bool(value: Any, label: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise ValueError(f"{label} is not boolean")


def validate_date(value: str, label: str) -> None:
    try:
        parsed = date.fromisoformat(value)
    except ValueError as error:
        raise ValueError(f"{label} date is not ISO YYYY-MM-DD") from error
    if parsed.isoformat() != value:
        raise ValueError(f"{label} date is not canonical ISO YYYY-MM-DD")


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def refuse_output_root(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"output root is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
