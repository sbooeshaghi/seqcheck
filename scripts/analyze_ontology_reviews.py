#!/usr/bin/env python3
"""Analyze completed blinded ontology reviews with locked weighted endpoints."""

from __future__ import annotations

import argparse
import math
import random
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import manage_ontology_reviews as reviews
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import manage_ontology_reviews as reviews
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
TOOL_VERSION = "0.1.0"
ADJUDICATION_FIELDS = (
    "adjudicator",
    "adjudication_date",
    "adjudicated_ontology_terms",
    "adjudication_category",
    "adjudication_rationale",
)
OUTCOME_FIELDS = (
    "final_ontology_terms",
    "final_assignment_source",
    "registry_exact_match",
)
SUMMARY_FIELDS = (
    "analysis_id",
    "metric",
    "estimate",
    "ci_lower",
    "ci_upper",
    "target",
    "target_met",
    "weighted",
    "description",
)
TERM_FIELDS = (
    "analysis_id",
    "ontology_term",
    "weighted_true_positive",
    "weighted_false_positive",
    "weighted_false_negative",
    "weighted_precision",
    "weighted_recall",
    "reviewer_weighted_kappa",
    "registry_sample_regions",
    "final_sample_regions",
    "sample_configurations",
)
DISAGREEMENT_FIELDS = (
    "analysis_id",
    "adjudication_category",
    "sample_regions",
    "weighted_regions",
    "weighted_fraction_of_disagreements",
)
BOOTSTRAP_FIELDS = (
    "analysis_id",
    "replicate",
    "sampled_configuration_clusters",
    "weighted_mapping_precision",
    "weighted_exact_set_agreement",
    "weighted_exact_set_cohen_kappa",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze a completed ontology review adjudication table."
    )
    parser.add_argument("--reviews-manifest", required=True, type=Path)
    parser.add_argument("--adjudicated-reviews", required=True, type=Path)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--yq-bin", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=30)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = analyze_reviews(
            reviews_manifest_path=args.reviews_manifest.resolve(),
            adjudicated_reviews_path=args.adjudicated_reviews.resolve(),
            registry_path=args.registry.resolve(),
            protocol_path=args.protocol.resolve(),
            yq_bin=args.yq_bin.resolve(),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
        )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"analyze_ontology_reviews: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def analyze_reviews(
    *,
    reviews_manifest_path: Path,
    adjudicated_reviews_path: Path,
    registry_path: Path,
    protocol_path: Path,
    yq_bin: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    review_manifest = runtime.load_json(reviews_manifest_path)
    if review_manifest.get("schema_version") != reviews.SCHEMA_VERSION:
        raise ValueError("ontology review manifest schema is unsupported")
    survey_id = str(review_manifest.get("survey_id", "")).strip()
    if not survey_id:
        raise ValueError("ontology review survey id is missing")
    verify_manifest_input(review_manifest, "registry", registry_path)
    verify_manifest_input(review_manifest, "protocol", protocol_path)
    template_identity = review_manifest.get("outputs", {}).get("ontology_reviews")
    if not isinstance(template_identity, dict):
        raise ValueError("ontology review template identity is missing")
    template_path = Path(str(template_identity.get("path", ""))).resolve()
    reviews.verify_file_identity(
        template_path, template_identity, "ontology review template"
    )
    template_rows, template_fields = reviews.read_csv(template_path)
    completed_rows, completed_fields = reviews.read_csv(adjudicated_reviews_path)
    if completed_fields != template_fields:
        raise ValueError("adjudicated ontology review fields changed")
    expected_rows = review_manifest.get("review_rows")
    if len(template_rows) != expected_rows or len(completed_rows) != expected_rows:
        raise ValueError("adjudicated ontology review row count changed")

    protocol = runtime.load_json(protocol_path)
    review_contract = reviews.validate_review_protocol(protocol, len(template_rows))
    registry = reviews.load_yaml_json(registry_path, yq_bin, timeout_seconds)
    _, term_order = reviews.build_term_reference(
        registry, review_contract["reviewable_term_statuses"]
    )
    source = {"term_order": term_order, "declared_terms": set(term_order)}
    outcomes = validate_and_finalize_rows(
        survey_id=survey_id,
        template_rows=template_rows,
        completed_rows=completed_rows,
        source=source,
        adjudication_categories=review_contract["adjudication_categories"],
    )
    script = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    dependency = runtime.script_identity(
        Path(reviews.__file__).resolve(), version=reviews.TOOL_VERSION
    )
    yq_identity = runtime.executable_identity(
        yq_bin, timeout_seconds=min(timeout_seconds, 30)
    )
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "survey_id": survey_id,
        "reviews_manifest_sha256": runtime.file_sha256(reviews_manifest_path),
        "review_template_sha256": runtime.file_sha256(template_path),
        "adjudicated_reviews_sha256": runtime.file_sha256(adjudicated_reviews_path),
        "registry_sha256": runtime.file_sha256(registry_path),
        "protocol_sha256": runtime.file_sha256(protocol_path),
        "tool": runtime.functional_script_identity(script),
        "review_runtime": runtime.functional_script_identity(dependency),
        "yq": runtime.functional_executable_identity(yq_identity),
    }
    analysis_id = runtime.sha256_json(stable_identity)[:16]
    analysis = review_contract["analysis"]
    bootstrap = bootstrap_metrics(analysis_id, outcomes, analysis)
    summary = summarize_primary_metrics(
        analysis_id=analysis_id,
        outcomes=outcomes,
        bootstrap=bootstrap,
        confidence_level=analysis["confidence_level"],
        targets=protocol["scientific_targets"],
    )
    term_summary = summarize_terms(analysis_id, outcomes, term_order)
    disagreement_summary = summarize_disagreements(analysis_id, outcomes)
    validation = validate_outputs(
        analysis_id=analysis_id,
        outcomes=outcomes,
        bootstrap=bootstrap,
        summary=summary,
        review_manifest=review_manifest,
        analysis=analysis,
    )
    if not validation["valid"]:
        raise ValueError(
            "ontology review analysis failed: " + "; ".join(validation["errors"])
        )

    tables = output_root / "tables"
    paths = {
        "outcomes": tables / "ontology_review_outcomes.csv",
        "summary": tables / "ontology_review_summary.csv",
        "terms": tables / "ontology_term_review_metrics.csv",
        "disagreements": tables / "ontology_disagreement_categories.csv",
        "bootstrap": output_root / "bootstrap" / "ontology_review_bootstrap.csv",
    }
    runtime.write_csv(paths["outcomes"], outcomes, [*template_fields, *OUTCOME_FIELDS])
    runtime.write_csv(paths["summary"], summary, list(SUMMARY_FIELDS))
    runtime.write_csv(paths["terms"], term_summary, list(TERM_FIELDS))
    runtime.write_csv(
        paths["disagreements"], disagreement_summary, list(DISAGREEMENT_FIELDS)
    )
    runtime.write_csv(paths["bootstrap"], bootstrap, list(BOOTSTRAP_FIELDS))
    validation_path = output_root / "validation" / "ontology_reviews.json"
    runtime.write_json(validation_path, validation)
    manifest_path = output_root / "manifests" / "ontology_review_analysis.json"
    manifest = {
        **stable_identity,
        "analysis_id": analysis_id,
        "generated_at": utc_now(),
        "valid": True,
        "frozen": True,
        "analysis_contract": analysis,
        "counts": validation["counts"],
        "targets": validation["targets"],
        "tools": {
            "analyzer": script,
            "review_runtime": dependency,
            "yq": yq_identity,
        },
        "inputs": {
            "reviews_manifest": runtime.file_identity(reviews_manifest_path),
            "review_template": runtime.file_identity(template_path),
            "adjudicated_reviews": runtime.file_identity(adjudicated_reviews_path),
            "registry": runtime.file_identity(registry_path),
            "protocol": runtime.file_identity(protocol_path),
        },
        "outputs": {
            **{name: runtime.file_identity(path) for name, path in paths.items()},
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def validate_and_finalize_rows(
    *,
    survey_id: str,
    template_rows: list[dict[str, str]],
    completed_rows: list[dict[str, str]],
    source: dict[str, Any],
    adjudication_categories: set[str],
) -> list[dict[str, Any]]:
    template_by_key = reviews.index_unique(
        template_rows, "region_key", "ontology review template"
    )
    completed_by_key = reviews.index_unique(
        completed_rows, "region_key", "adjudicated ontology review"
    )
    if set(template_by_key) != set(completed_by_key):
        raise ValueError("adjudicated ontology review region keys changed")
    immutable = [
        field for field in template_rows[0] if field not in ADJUDICATION_FIELDS
    ]
    result = []
    for key, template in template_by_key.items():
        row = completed_by_key[key]
        changed = [field for field in immutable if row[field] != template[field]]
        if changed:
            raise ValueError(
                f"{key}: adjudication changed reviewer evidence: {', '.join(changed)}"
            )
        if row.get("survey_id") != survey_id:
            raise ValueError(f"{key}: adjudication survey id changed")
        agreement = parse_bool(row["exact_term_set_agreement"], "agreement")
        needs_adjudication = parse_bool(row["needs_adjudication"], "adjudication")
        if agreement == needs_adjudication:
            raise ValueError(f"{key}: agreement and adjudication flags conflict")
        left = canonical_terms(row["reviewer_1_ontology_terms"], source)
        right = canonical_terms(row["reviewer_2_ontology_terms"], source)
        truth = canonical_terms(row["ontology_terms"], source)
        if agreement:
            consensus = canonical_terms(row["consensus_ontology_terms"], source)
            if left != right or consensus != left:
                raise ValueError(f"{key}: reviewer consensus is inconsistent")
            if any(row[field].strip() for field in ADJUDICATION_FIELDS):
                raise ValueError(f"{key}: agreed row contains adjudication values")
            final_terms = consensus
            final_source = "reviewer_consensus"
        else:
            if left == right or row["consensus_ontology_terms"].strip():
                raise ValueError(f"{key}: disagreement row contains consensus")
            missing = [field for field in ADJUDICATION_FIELDS if not row[field].strip()]
            if missing:
                raise ValueError(
                    f"{key}: adjudication is incomplete: {', '.join(missing)}"
                )
            reviews.validate_date(row["adjudication_date"], "adjudication")
            category = row["adjudication_category"].strip()
            if category not in adjudication_categories:
                raise ValueError(f"{key}: adjudication category is invalid")
            final_terms = canonical_terms(row["adjudicated_ontology_terms"], source)
            final_source = "human_adjudication"
        final_values = reviews.split_terms(final_terms)
        if reviews.UNKNOWN_TERM in final_values and len(final_values) != 1:
            raise ValueError(
                f"{key}: final assignment combines unknown and precise terms"
            )
        try:
            weight = float(row["sampling_weight"])
        except ValueError as error:
            raise ValueError(f"{key}: sampling weight is invalid") from error
        if not math.isfinite(weight) or weight <= 0:
            raise ValueError(f"{key}: sampling weight must be finite and positive")
        if not row["configuration_accession"].strip():
            raise ValueError(f"{key}: configuration accession is missing")
        result.append(
            {
                **row,
                "final_ontology_terms": final_terms,
                "final_assignment_source": final_source,
                "registry_exact_match": final_terms == truth,
                "_weight": weight,
                "_reviewer_1_terms": left,
                "_reviewer_2_terms": right,
                "_registry_terms": truth,
            }
        )
    return sorted(result, key=lambda row: row["region_key"])


def bootstrap_metrics(
    analysis_id: str, rows: list[dict[str, Any]], analysis: dict[str, Any]
) -> list[dict[str, Any]]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["configuration_accession"]].append(row)
    configurations = sorted(grouped)
    if len(configurations) < 2:
        raise ValueError("ontology review bootstrap needs at least two configurations")
    rng = random.Random(analysis["bootstrap_seed"])
    result = []
    for replicate in range(1, analysis["bootstrap_replicates"] + 1):
        sampled = [rng.choice(configurations) for _ in configurations]
        replicate_rows = [row for accession in sampled for row in grouped[accession]]
        result.append(
            {
                "analysis_id": analysis_id,
                "replicate": replicate,
                "sampled_configuration_clusters": len(sampled),
                "weighted_mapping_precision": weighted_mean(
                    replicate_rows, "registry_exact_match"
                ),
                "weighted_exact_set_agreement": weighted_agreement(replicate_rows),
                "weighted_exact_set_cohen_kappa": weighted_kappa(
                    replicate_rows,
                    "_reviewer_1_terms",
                    "_reviewer_2_terms",
                ),
            }
        )
    return result


def summarize_primary_metrics(
    *,
    analysis_id: str,
    outcomes: list[dict[str, Any]],
    bootstrap: list[dict[str, Any]],
    confidence_level: float,
    targets: dict[str, Any],
) -> list[dict[str, Any]]:
    specifications = (
        (
            "weighted_mapping_precision",
            weighted_mean(outcomes, "registry_exact_match"),
            "weighted_mapping_precision",
            targets["weighted_mapping_precision_min"],
            "Inverse-inclusion weighted exact registry-to-final term-set precision.",
        ),
        (
            "weighted_exact_set_agreement",
            weighted_agreement(outcomes),
            "weighted_exact_set_agreement",
            "",
            "Inverse-inclusion weighted reviewer exact term-set agreement.",
        ),
        (
            "weighted_exact_set_cohen_kappa",
            weighted_kappa(outcomes, "_reviewer_1_terms", "_reviewer_2_terms"),
            "weighted_exact_set_cohen_kappa",
            targets["cohen_kappa_min"],
            "Population-weighted Cohen's kappa over exact ontology term sets.",
        ),
        (
            "unweighted_mapping_precision",
            unweighted_mean(outcomes, "registry_exact_match"),
            None,
            "",
            "Sample exact registry-to-final term-set precision without weights.",
        ),
    )
    rows = []
    for metric, estimate, bootstrap_field, target, description in specifications:
        values = (
            [
                row[bootstrap_field]
                for row in bootstrap
                if row[bootstrap_field] not in ("", None)
            ]
            if bootstrap_field
            else []
        )
        lower, upper = confidence_interval(values, confidence_level)
        rows.append(
            {
                "analysis_id": analysis_id,
                "metric": metric,
                "estimate": optional_number(estimate),
                "ci_lower": optional_number(lower),
                "ci_upper": optional_number(upper),
                "target": target,
                "target_met": (
                    "" if target == "" or estimate is None else estimate >= target
                ),
                "weighted": metric.startswith("weighted_"),
                "description": description,
            }
        )
    return rows


def summarize_terms(
    analysis_id: str, rows: list[dict[str, Any]], term_order: list[str]
) -> list[dict[str, Any]]:
    result = []
    for term in term_order:
        predicted = [
            term in reviews.split_terms(row["_registry_terms"]) for row in rows
        ]
        final = [
            term in reviews.split_terms(row["final_ontology_terms"]) for row in rows
        ]
        if not any(predicted) and not any(final):
            continue
        tp = weighted_condition(rows, lambda index: predicted[index] and final[index])
        fp = weighted_condition(
            rows, lambda index: predicted[index] and not final[index]
        )
        fn = weighted_condition(
            rows, lambda index: not predicted[index] and final[index]
        )
        reviewer_rows = [
            {
                **row,
                "_left_binary": "present"
                if term in reviews.split_terms(row["_reviewer_1_terms"])
                else "absent",
                "_right_binary": "present"
                if term in reviews.split_terms(row["_reviewer_2_terms"])
                else "absent",
            }
            for row in rows
        ]
        configurations = {
            row["configuration_accession"]
            for index, row in enumerate(rows)
            if predicted[index] or final[index]
        }
        result.append(
            {
                "analysis_id": analysis_id,
                "ontology_term": term,
                "weighted_true_positive": tp,
                "weighted_false_positive": fp,
                "weighted_false_negative": fn,
                "weighted_precision": optional_ratio(tp, tp + fp),
                "weighted_recall": optional_ratio(tp, tp + fn),
                "reviewer_weighted_kappa": optional_number(
                    weighted_kappa(reviewer_rows, "_left_binary", "_right_binary")
                ),
                "registry_sample_regions": sum(predicted),
                "final_sample_regions": sum(final),
                "sample_configurations": len(configurations),
            }
        )
    return result


def summarize_disagreements(
    analysis_id: str, rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    disagreements = [
        row for row in rows if row["final_assignment_source"] == "human_adjudication"
    ]
    total_weight = sum(row["_weight"] for row in disagreements)
    grouped = defaultdict(list)
    for row in disagreements:
        grouped[row["adjudication_category"]].append(row)
    return [
        {
            "analysis_id": analysis_id,
            "adjudication_category": category,
            "sample_regions": len(values),
            "weighted_regions": sum(row["_weight"] for row in values),
            "weighted_fraction_of_disagreements": optional_ratio(
                sum(row["_weight"] for row in values), total_weight
            ),
        }
        for category, values in sorted(grouped.items())
    ]


def validate_outputs(
    *,
    analysis_id: str,
    outcomes: list[dict[str, Any]],
    bootstrap: list[dict[str, Any]],
    summary: list[dict[str, Any]],
    review_manifest: dict[str, Any],
    analysis: dict[str, Any],
) -> dict[str, Any]:
    errors = []
    if len(outcomes) != review_manifest.get("review_rows"):
        errors.append("outcome rows do not reconcile with review manifest")
    if len({row["region_key"] for row in outcomes}) != len(outcomes):
        errors.append("outcome region keys are not unique")
    if any(row["survey_id"] != review_manifest.get("survey_id") for row in outcomes):
        errors.append("outcome survey identifiers differ")
    completed_adjudications = sum(
        row["final_assignment_source"] == "human_adjudication" for row in outcomes
    )
    if completed_adjudications != review_manifest.get("adjudication_required"):
        errors.append("completed adjudications do not reconcile with review manifest")
    if review_manifest.get("exact_term_set_agreements") not in (
        None,
        len(outcomes) - completed_adjudications,
    ):
        errors.append("reviewer agreements do not reconcile with review manifest")
    if len(bootstrap) != analysis["bootstrap_replicates"]:
        errors.append("bootstrap replicate count differs from protocol")
    if any(row["analysis_id"] != analysis_id for row in bootstrap + summary):
        errors.append("analysis identifiers differ")
    metrics = {row["metric"]: row for row in summary}
    required = {
        "weighted_mapping_precision",
        "weighted_exact_set_agreement",
        "weighted_exact_set_cohen_kappa",
        "unweighted_mapping_precision",
    }
    if set(metrics) != required:
        errors.append("primary metric rows are incomplete")
    return {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": analysis_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "counts": {
            "reviewed_regions": len(outcomes),
            "configurations": len({row["configuration_accession"] for row in outcomes}),
            "reviewer_agreements": sum(
                parse_bool(row["exact_term_set_agreement"], "agreement")
                for row in outcomes
            ),
            "human_adjudications": sum(
                row["final_assignment_source"] == "human_adjudication"
                for row in outcomes
            ),
            "registry_exact_matches": sum(
                row["registry_exact_match"] for row in outcomes
            ),
            "bootstrap_replicates": len(bootstrap),
        },
        "targets": {
            row["metric"]: {
                "estimate": row["estimate"],
                "target": row["target"],
                "met": row["target_met"],
            }
            for row in summary
            if row["target"] != ""
        },
    }


def weighted_mean(rows: list[dict[str, Any]], field: str) -> float:
    denominator = sum(row["_weight"] for row in rows)
    return sum(row["_weight"] * bool(row[field]) for row in rows) / denominator


def unweighted_mean(rows: list[dict[str, Any]], field: str) -> float:
    return sum(bool(row[field]) for row in rows) / len(rows)


def weighted_agreement(rows: list[dict[str, Any]]) -> float:
    denominator = sum(row["_weight"] for row in rows)
    return (
        sum(
            row["_weight"] * (row["_reviewer_1_terms"] == row["_reviewer_2_terms"])
            for row in rows
        )
        / denominator
    )


def weighted_kappa(
    rows: list[dict[str, Any]], left_field: str, right_field: str
) -> float | None:
    total = sum(row["_weight"] for row in rows)
    left = Counter()
    right = Counter()
    observed = 0.0
    for row in rows:
        weight = row["_weight"]
        left[row[left_field]] += weight
        right[row[right_field]] += weight
        observed += weight * (row[left_field] == row[right_field])
    observed /= total
    expected = sum(
        (left[category] / total) * (right[category] / total)
        for category in set(left) | set(right)
    )
    if math.isclose(expected, 1.0):
        return None
    return (observed - expected) / (1 - expected)


def weighted_condition(rows: list[dict[str, Any]], predicate: Any) -> float:
    return sum(row["_weight"] for index, row in enumerate(rows) if predicate(index))


def confidence_interval(
    values: list[float], confidence_level: float
) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    alpha = (1 - confidence_level) / 2
    ordered = sorted(float(value) for value in values)
    return percentile(ordered, alpha), percentile(ordered, 1 - alpha)


def percentile(values: list[float], probability: float) -> float:
    position = (len(values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    fraction = position - lower
    return values[lower] * (1 - fraction) + values[upper] * fraction


def optional_ratio(numerator: float, denominator: float) -> float | str:
    return numerator / denominator if denominator else ""


def optional_number(value: float | None) -> float | str:
    return "" if value is None else value


def canonical_terms(value: str, source: dict[str, Any]) -> str:
    result = reviews.canonical_terms(value, source)
    terms = reviews.split_terms(result)
    if reviews.UNKNOWN_TERM in terms and len(terms) != 1:
        raise ValueError("ontology assignment combines unknown and precise terms")
    return result


def parse_bool(value: str, label: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise ValueError(f"{label} boolean is invalid")


def verify_manifest_input(manifest: dict[str, Any], name: str, path: Path) -> None:
    identity = manifest.get("inputs", {}).get(name)
    if not isinstance(identity, dict):
        raise ValueError(f"ontology review input identity is missing: {name}")
    reviews.verify_file_identity(path, identity, f"ontology review {name}")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
