#!/usr/bin/env python3
"""Compare region-reference models with and without ontology semantics."""

from __future__ import annotations

import argparse
import json
import math
import random
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import analyze_perturbation_calibration as perturbation_analysis
    import paper_runtime as runtime
    import run_perturbation_calibration as execution_runner
except ModuleNotFoundError:
    from scripts import analyze_perturbation_calibration as perturbation_analysis
    from scripts import paper_runtime as runtime
    from scripts import run_perturbation_calibration as execution_runner


SCHEMA_VERSION = "0.1.0"
ANALYZER_VERSION = "0.1.0"
OBSERVATION_FIELDS = (
    "analysis_id",
    "split",
    "observation_id",
    "label",
    "condition_class",
    "configuration_accession",
    "family_id",
    "modality",
    "operator_id",
    "condition_count",
    "condition_ids",
    "check",
    "metric_name",
    "unit",
    "files",
    "reads",
    "region_id",
    "sequence_type",
    "ontology_signature",
    "ontology_terms",
    "value",
    "eligible",
)
TERM_COVERAGE_FIELDS = (
    "analysis_id",
    "metric_name",
    "ontology_term",
    "calibration_configurations",
    "evaluation_configurations",
    "minimum_configurations_per_split",
    "eligible",
)
REFERENCE_FIELDS = (
    "analysis_id",
    "model_id",
    "model_name",
    "group_key_json",
    "group_fields",
    "independent_configurations",
    "clean_observations",
    "reference_value",
)
SCORE_FIELDS = (
    "analysis_id",
    "split",
    "model_id",
    "model_name",
    "observation_id",
    "configuration_accession",
    "label",
    "metric_name",
    "region_id",
    "ontology_signature",
    "group_key_json",
    "reference_value",
    "observed_value",
    "anomaly_score",
    "threshold",
    "detected",
)
MODEL_SUMMARY_FIELDS = (
    "analysis_id",
    "split",
    "model_id",
    "model_name",
    "eligible_observations",
    "scored_observations",
    "coverage",
    "independent_configurations",
    "positive_observations",
    "negative_observations",
    "threshold",
    "sensitivity",
    "false_positive_rate",
    "average_precision",
    "sensitivity_calibration_error",
)
COMPARISON_FIELDS = (
    "analysis_id",
    "comparison_id",
    "primary_model_id",
    "baseline_model_id",
    "common_observations",
    "independent_configurations",
    "positive_observations",
    "negative_observations",
    "primary_sensitivity",
    "baseline_sensitivity",
    "primary_false_positive_rate",
    "baseline_false_positive_rate",
    "false_positive_rate_difference",
    "relative_false_positive_rate_reduction",
    "primary_average_precision",
    "baseline_average_precision",
    "average_precision_difference",
    "bootstrap_ci_lower",
    "bootstrap_ci_upper",
    "confidence_level",
    "bootstrap_replicates",
    "relative_reduction_target",
    "upper_ci_target",
    "scientific_target_met",
)
TERM_COMPARISON_FIELDS = (
    "analysis_id",
    "metric_name",
    "ontology_term",
    "independent_configurations",
    "positive_observations",
    "negative_observations",
    "primary_false_positive_rate",
    "baseline_false_positive_rate",
    "false_positive_rate_difference",
    "primary_sensitivity",
    "baseline_sensitivity",
    "primary_average_precision",
    "baseline_average_precision",
)
BOOTSTRAP_FIELDS = (
    "analysis_id",
    "replicate",
    "sampled_configurations",
    "primary_false_positive_rate",
    "baseline_false_positive_rate",
    "false_positive_rate_difference",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Measure the held-out QC value added by region ontology terms."
    )
    parser.add_argument("--calibration-execution", required=True, type=Path)
    parser.add_argument("--evaluation-execution", required=True, type=Path)
    parser.add_argument("--analysis-protocol", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = analyze_ontology_utility(
            calibration_execution_path=args.calibration_execution.resolve(),
            evaluation_execution_path=args.evaluation_execution.resolve(),
            analysis_protocol_path=args.analysis_protocol.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"analyze_ontology_utility: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def analyze_ontology_utility(
    *,
    calibration_execution_path: Path,
    evaluation_execution_path: Path,
    analysis_protocol_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    protocol = runtime.load_json(analysis_protocol_path)
    validate_protocol(protocol)
    calibration = load_execution(calibration_execution_path, "sampling_study")
    evaluation = load_execution(evaluation_execution_path, "policy_sample_bundle")
    calibration_configurations = {
        value["configuration_accession"] for value in calibration["conditions"]
    }
    evaluation_configurations = {
        value["configuration_accession"] for value in evaluation["conditions"]
    }
    if protocol["evaluation"]["require_disjoint_configurations"] and (
        calibration_configurations & evaluation_configurations
    ):
        raise ValueError("calibration and evaluation configurations overlap")

    tools = {
        "analyzer": runtime.script_identity(
            Path(__file__).resolve(), version=ANALYZER_VERSION
        ),
        "perturbation_analyzer": runtime.script_identity(
            Path(perturbation_analysis.__file__).resolve(),
            version=perturbation_analysis.ANALYZER_VERSION,
        ),
        "execution_runner": runtime.script_identity(
            Path(execution_runner.__file__).resolve(),
            version=execution_runner.RUNNER_VERSION,
        ),
        "runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
    }
    stable = {
        "schema_version": SCHEMA_VERSION,
        "calibration_execution_id": calibration["execution"]["execution_id"],
        "calibration_execution_sha256": runtime.file_sha256(calibration_execution_path),
        "evaluation_execution_id": evaluation["execution"]["execution_id"],
        "evaluation_execution_sha256": runtime.file_sha256(evaluation_execution_path),
        "analysis_protocol_sha256": runtime.file_sha256(analysis_protocol_path),
        "analyzer": runtime.functional_script_identity(tools["analyzer"]),
        "perturbation_analyzer": runtime.functional_script_identity(
            tools["perturbation_analyzer"]
        ),
        "execution_runner": runtime.functional_script_identity(
            tools["execution_runner"]
        ),
        "runtime": runtime.functional_script_identity(tools["runtime"]),
    }
    analysis_id = runtime.sha256_json(stable)[:16]

    observations = [
        *build_observations(
            split="calibration",
            conditions=calibration["conditions"],
            metrics=calibration["metrics"],
            protocol=protocol,
        ),
        *build_observations(
            split="evaluation",
            conditions=evaluation["conditions"],
            metrics=evaluation["metrics"],
            protocol=protocol,
        ),
    ]
    coverage, eligible_terms = build_term_coverage(observations, protocol, analysis_id)
    for value in observations:
        terms = split_terms(value["ontology_terms"])
        value["eligible"] = bool(terms) and all(
            (value["metric_name"], term) in eligible_terms for term in terms
        )
        value["analysis_id"] = analysis_id
    eligible = [value for value in observations if value["eligible"]]
    references, scores, summaries = fit_and_score_models(
        eligible, protocol, analysis_id
    )
    comparison, bootstrap = compare_primary_models(
        scores=scores,
        protocol=protocol,
        analysis_id=analysis_id,
    )
    term_comparisons = compare_terms(
        scores=scores,
        observations=eligible,
        comparison=comparison,
        protocol=protocol,
        analysis_id=analysis_id,
    )
    validation = validate_outputs(
        analysis_id=analysis_id,
        observations=observations,
        eligible=eligible,
        coverage=coverage,
        references=references,
        scores=scores,
        summaries=summaries,
        comparison=comparison,
        bootstrap=bootstrap,
        protocol=protocol,
    )
    if not validation["valid"]:
        raise ValueError(
            "ontology utility analysis failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )

    output_root.mkdir(parents=True)
    paths = {
        "observations": output_root / "tables" / "observations.csv",
        "term_coverage": output_root / "tables" / "term_coverage.csv",
        "references": output_root / "tables" / "model_references.csv",
        "scores": output_root / "tables" / "model_scores.csv",
        "model_summary": output_root / "tables" / "model_summary.csv",
        "primary_comparison": output_root / "tables" / "primary_comparison.csv",
        "term_comparisons": output_root / "tables" / "term_comparisons.csv",
        "bootstrap": output_root / "tables" / "bootstrap_fpr_difference.csv",
        "validation": output_root / "validation" / "ontology_utility.json",
    }
    runtime.write_csv(paths["observations"], observations, list(OBSERVATION_FIELDS))
    runtime.write_csv(paths["term_coverage"], coverage, list(TERM_COVERAGE_FIELDS))
    runtime.write_csv(paths["references"], references, list(REFERENCE_FIELDS))
    runtime.write_csv(paths["scores"], scores, list(SCORE_FIELDS))
    runtime.write_csv(paths["model_summary"], summaries, list(MODEL_SUMMARY_FIELDS))
    runtime.write_csv(
        paths["primary_comparison"], [comparison], list(COMPARISON_FIELDS)
    )
    runtime.write_csv(
        paths["term_comparisons"], term_comparisons, list(TERM_COMPARISON_FIELDS)
    )
    runtime.write_csv(paths["bootstrap"], bootstrap, list(BOOTSTRAP_FIELDS))
    runtime.write_json(paths["validation"], validation)
    manifest_path = output_root / "manifests" / "analysis.json"
    manifest = {
        **stable,
        "analysis_id": analysis_id,
        "created_at": utc_now(),
        "valid": True,
        "scientific_target_met": comparison["scientific_target_met"],
        "inputs": {
            "calibration_execution": runtime.file_identity(calibration_execution_path),
            "evaluation_execution": runtime.file_identity(evaluation_execution_path),
            "analysis_protocol": runtime.file_identity(analysis_protocol_path),
        },
        "tools": tools,
        "counts": validation["counts"],
        "outputs": {key: runtime.file_identity(path) for key, path in paths.items()},
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_execution(path: Path, source_kind: str) -> dict[str, Any]:
    execution, conditions, runs, metrics, assessments = (
        perturbation_analysis.load_execution(path)
    )
    materialization_path = perturbation_analysis.verified_path(
        execution.get("inputs", {}).get("materialization", {}),
        "materialization manifest",
    )
    materialization, _ = execution_runner.load_materialization(materialization_path)
    observed_kind = materialization.get("sample_source", {}).get("kind")
    if observed_kind != source_kind:
        raise ValueError(
            f"expected {source_kind} execution, observed {observed_kind or 'missing'}"
        )
    if execution["execution_id"] == "":
        raise ValueError("execution identifier is empty")
    return {
        "execution": execution,
        "materialization": materialization,
        "conditions": conditions,
        "runs": runs,
        "metrics": metrics,
        "assessments": assessments,
    }


def validate_protocol(protocol: dict[str, Any]) -> None:
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("ontology utility protocol schema is unsupported")
    observations = protocol.get("observations", {})
    if observations.get("negative_condition_kind") != "clean":
        raise ValueError("ontology utility negatives must be clean conditions")
    if set(observations.get("positive_condition_kinds", [])) != {
        "deterministic",
        "stochastic",
    }:
        raise ValueError("ontology utility positive condition kinds are invalid")
    anchor = finite_number(
        observations.get("stochastic_anchor_fraction"), "stochastic anchor"
    )
    if not 0 < anchor < 1:
        raise ValueError("stochastic anchor must be between zero and one")
    if observations.get("positive_aggregation") != (
        "configuration_region_metric_operator_median"
    ):
        raise ValueError("ontology utility positive aggregation is unsupported")
    for field in (
        "positive_requires_target_region",
        "positive_requires_expected_metric",
    ):
        if observations.get(field) is not True:
            raise ValueError(f"ontology utility observation rule is disabled: {field}")
    metrics = observations.get("metrics")
    if not isinstance(metrics, list) or not metrics:
        raise ValueError("ontology utility metrics are missing")
    keys = set()
    for metric in metrics:
        key = (required_string(metric, "check"), required_string(metric, "metric_name"))
        if key in keys:
            raise ValueError(f"duplicate ontology utility metric: {key}")
        keys.add(key)
        required_string(metric, "unit")
        required_string(metric, "sequence_type")
        if metric.get("anomaly_direction") != "decrease":
            raise ValueError("only decreasing ontology utility metrics are supported")
        terms = metric.get("ontology_terms")
        if (
            not isinstance(terms, list)
            or not terms
            or any(
                not isinstance(value, str) or not value.startswith("RGN:")
                for value in terms
            )
        ):
            raise ValueError(f"ontology utility metric terms are invalid: {key}")
    minimum = positive_int(
        protocol.get("eligibility", {}).get(
            "minimum_independent_configurations_per_term_per_split"
        ),
        "minimum term configurations",
    )
    if minimum < 1:
        raise ValueError("minimum term configurations must be positive")
    if protocol["eligibility"].get("require_all_signature_terms_eligible") is not True:
        raise ValueError("all ontology signature terms must meet eligibility")
    models = protocol.get("models")
    if not isinstance(models, list) or {value.get("model_id") for value in models} != {
        "M1",
        "M2",
        "M3",
        "M4",
    }:
        raise ValueError("ontology utility models must be M1 through M4")
    expected_groups = {
        "M1": ["metric_name"],
        "M2": ["metric_name", "sequence_type"],
        "M3": ["metric_name", "sequence_type", "ontology_signature"],
        "M4": ["metric_name", "family_id", "ontology_signature"],
    }
    for model in models:
        fields = model.get("group_fields")
        if fields != expected_groups[model["model_id"]]:
            raise ValueError(f"model group fields are invalid: {model.get('model_id')}")
        positive_int(
            model.get("minimum_reference_configurations"),
            "minimum reference configurations",
        )
    fit = protocol.get("fit", {})
    sensitivity = finite_number(fit.get("sensitivity_target"), "sensitivity target")
    if not 0 < sensitivity <= 1:
        raise ValueError("sensitivity target must be in (0, 1]")
    if fit.get("reference_statistic") != "median_of_configuration_medians":
        raise ValueError("reference statistic is unsupported")
    expected_fit = {
        "anomaly_score": "reference_minus_observed",
        "threshold_selection": (
            "largest_score_cutoff_reaching_target_with_ties_included"
        ),
        "unseen_group_policy": "exclude_and_report_coverage",
    }
    for field, expected in expected_fit.items():
        if fit.get(field) != expected:
            raise ValueError(f"ontology utility fit rule is unsupported: {field}")
    evaluation = protocol.get("evaluation", {})
    if (
        evaluation.get("primary_model") != "M3"
        or evaluation.get("baseline_model") != "M2"
    ):
        raise ValueError("primary ontology comparison must be M3 versus M2")
    expected_evaluation = {
        "operating_point": "calibration_threshold_fixed_at_target_sensitivity",
        "calibration_metric": "absolute_held_out_sensitivity_minus_target",
        "ranking_metric": "tie_aware_average_precision",
    }
    for field, expected in expected_evaluation.items():
        if evaluation.get(field) != expected:
            raise ValueError(
                f"ontology utility evaluation rule is unsupported: {field}"
            )
    if evaluation.get("require_disjoint_configurations") is not True:
        raise ValueError("calibration and evaluation configurations must be disjoint")
    statistics_value = protocol.get("statistics", {})
    if statistics_value.get("bootstrap_unit") != "configuration_accession" or (
        statistics_value.get("bootstrap_method") != "percentile_cluster_bootstrap"
    ):
        raise ValueError("ontology utility bootstrap rule is unsupported")
    positive_int(statistics_value.get("bootstrap_replicates"), "bootstrap replicates")
    seed = statistics_value.get("bootstrap_seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("bootstrap seed must be a nonnegative integer")
    confidence = finite_number(statistics_value.get("confidence_level"), "confidence")
    if not 0 < confidence < 1:
        raise ValueError("confidence level must be between zero and one")
    targets = protocol.get("scientific_targets", {})
    reduction = finite_number(
        targets.get("relative_false_positive_rate_reduction_min"),
        "relative false-positive-rate reduction target",
    )
    if not 0 <= reduction <= 1:
        raise ValueError("relative false-positive-rate reduction target is invalid")
    finite_number(
        targets.get("false_positive_rate_difference_upper_ci_max"),
        "false-positive-rate confidence bound target",
    )


def build_observations(
    *,
    split: str,
    conditions: list[dict[str, Any]],
    metrics: list[dict[str, str]],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    metric_specs = {
        (value["check"], value["metric_name"]): value
        for value in protocol["observations"]["metrics"]
    }
    conditions_by_id = {value["condition_id"]: value for value in conditions}
    if len(conditions_by_id) != len(conditions):
        raise ValueError(f"{split}: condition identifiers are not unique")
    clean = []
    positives = defaultdict(list)
    anchor = protocol["observations"]["stochastic_anchor_fraction"]
    for row in metrics:
        spec = metric_specs.get((row.get("check"), row.get("metric_name")))
        if spec is None or row.get("metric_side") != "observed":
            continue
        if row.get("data_kind") != "scalar" or row.get("unit") != spec["unit"]:
            continue
        value = numeric_json(row.get("value_json", ""))
        if value is None:
            continue
        annotation = single_region_annotation(row, split)
        if annotation["sequence_type"] != spec["sequence_type"]:
            continue
        report_terms = split_terms(row.get("ontology_terms", ""))
        annotation_terms = set(annotation["ontology_terms"])
        if report_terms != annotation_terms:
            raise ValueError(
                f"{split}: report and seqspec ontology terms differ for "
                f"{annotation['region_id']}"
            )
        relevant_terms = sorted(annotation_terms & set(spec["ontology_terms"]))
        if not relevant_terms:
            continue
        condition_id = row.get("condition_id", "")
        condition = conditions_by_id.get(condition_id)
        if condition is None:
            raise ValueError(f"{split}: metric has unknown condition {condition_id}")
        base = {
            "split": split,
            "configuration_accession": condition["configuration_accession"],
            "family_id": condition["family_id"],
            "modality": condition["modality"],
            "check": spec["check"],
            "metric_name": spec["metric_name"],
            "unit": spec["unit"],
            "files": row.get("files", ""),
            "reads": row.get("reads", ""),
            "region_id": annotation["region_id"],
            "sequence_type": annotation["sequence_type"],
            "ontology_signature": "+".join(relevant_terms),
            "ontology_terms": ";".join(relevant_terms),
            "value": value,
        }
        kind = condition["condition_kind"]
        if kind == protocol["observations"]["negative_condition_kind"]:
            clean.append(
                finalize_observation(
                    {
                        **base,
                        "label": 0,
                        "condition_class": "clean",
                        "operator_id": "CLEAN",
                        "condition_count": 1,
                        "condition_ids": condition_id,
                    }
                )
            )
            continue
        if kind not in protocol["observations"]["positive_condition_kinds"]:
            continue
        if kind == "stochastic" and not math.isclose(
            float(condition["event_fraction"]), anchor
        ):
            continue
        if annotation["region_id"] not in set(
            condition["target"].get("region_ids", [])
        ):
            continue
        if spec["metric_name"] not in set(
            condition["expected"].get("metric_names", [])
        ):
            continue
        key = (
            base["configuration_accession"],
            base["family_id"],
            base["modality"],
            base["check"],
            base["metric_name"],
            base["unit"],
            base["files"],
            base["reads"],
            base["region_id"],
            base["sequence_type"],
            base["ontology_signature"],
            base["ontology_terms"],
            condition["operator_id"],
        )
        positives[key].append((value, condition_id, base))

    result = clean
    for key, values in positives.items():
        base = values[0][2]
        condition_ids = sorted({value[1] for value in values})
        result.append(
            finalize_observation(
                {
                    **base,
                    "label": 1,
                    "condition_class": "targeted_perturbation",
                    "operator_id": key[-1],
                    "condition_count": len(condition_ids),
                    "condition_ids": ";".join(condition_ids),
                    "value": statistics.median(value[0] for value in values),
                }
            )
        )
    identifiers = [value["observation_id"] for value in result]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"{split}: observation identifiers are not unique")
    return sorted(result, key=lambda value: value["observation_id"])


def single_region_annotation(row: dict[str, str], split: str) -> dict[str, Any]:
    try:
        annotations = json.loads(row.get("region_annotations_json", ""))
    except json.JSONDecodeError as error:
        raise ValueError(
            f"{split}: metric region annotations are invalid JSON"
        ) from error
    if not isinstance(annotations, list) or len(annotations) != 1:
        raise ValueError(f"{split}: role metric does not resolve to exactly one region")
    annotation = annotations[0]
    if not isinstance(annotation, dict):
        raise ValueError(f"{split}: role metric region annotation is not an object")
    region_id = required_string(annotation, "region_id")
    sequence_type = required_string(annotation, "sequence_type")
    terms = annotation.get("ontology_terms")
    if (
        not isinstance(terms, list)
        or not terms
        or any(
            not isinstance(value, str) or not value.startswith("RGN:")
            for value in terms
        )
    ):
        raise ValueError(f"{split}: role metric ontology annotation is invalid")
    if len(terms) != len(set(terms)):
        raise ValueError(f"{split}: role metric ontology terms are duplicated")
    if row.get("regions") != region_id or row.get("sequence_types") != sequence_type:
        raise ValueError(f"{split}: flattened region annotation does not reconcile")
    return {
        "region_id": region_id,
        "sequence_type": sequence_type,
        "ontology_terms": sorted(terms),
    }


def finalize_observation(value: dict[str, Any]) -> dict[str, Any]:
    stable = {key: item for key, item in value.items() if key != "eligible"}
    return {
        **stable,
        "observation_id": runtime.sha256_json(stable)[:20],
        "eligible": False,
    }


def build_term_coverage(
    observations: list[dict[str, Any]],
    protocol: dict[str, Any],
    analysis_id: str,
) -> tuple[list[dict[str, Any]], set[tuple[str, str]]]:
    configurations = defaultdict(set)
    declared = set()
    for metric in protocol["observations"]["metrics"]:
        for term in metric["ontology_terms"]:
            declared.add((metric["metric_name"], term))
    for value in observations:
        if value["label"] != 0:
            continue
        for term in split_terms(value["ontology_terms"]):
            configurations[(value["split"], value["metric_name"], term)].add(
                value["configuration_accession"]
            )
    minimum = protocol["eligibility"][
        "minimum_independent_configurations_per_term_per_split"
    ]
    rows = []
    eligible = set()
    for metric_name, term in sorted(declared):
        calibration = len(configurations[("calibration", metric_name, term)])
        evaluation = len(configurations[("evaluation", metric_name, term)])
        is_eligible = calibration >= minimum and evaluation >= minimum
        if is_eligible:
            eligible.add((metric_name, term))
        rows.append(
            {
                "analysis_id": analysis_id,
                "metric_name": metric_name,
                "ontology_term": term,
                "calibration_configurations": calibration,
                "evaluation_configurations": evaluation,
                "minimum_configurations_per_split": minimum,
                "eligible": is_eligible,
            }
        )
    return rows, eligible


def fit_and_score_models(
    observations: list[dict[str, Any]],
    protocol: dict[str, Any],
    analysis_id: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    calibration = [value for value in observations if value["split"] == "calibration"]
    target = protocol["fit"]["sensitivity_target"]
    references_out = []
    scores_out = []
    summaries = []
    for model in protocol["models"]:
        fields = model["group_fields"]
        grouped = defaultdict(lambda: defaultdict(list))
        for value in calibration:
            if value["label"] == 0:
                grouped[group_key(value, fields)][
                    value["configuration_accession"]
                ].append(value["value"])
        references = {}
        for key, by_configuration in sorted(grouped.items()):
            if len(by_configuration) < model["minimum_reference_configurations"]:
                continue
            configuration_medians = [
                statistics.median(values) for values in by_configuration.values()
            ]
            references[key] = statistics.median(configuration_medians)
            references_out.append(
                {
                    "analysis_id": analysis_id,
                    "model_id": model["model_id"],
                    "model_name": model["name"],
                    "group_key_json": runtime.canonical_json(list(key)),
                    "group_fields": ";".join(fields),
                    "independent_configurations": len(by_configuration),
                    "clean_observations": sum(
                        len(values) for values in by_configuration.values()
                    ),
                    "reference_value": references[key],
                }
            )
        calibration_scores = score_observations(calibration, fields, references)
        positive_scores = [
            value["anomaly_score"]
            for value in calibration_scores
            if value["label"] == 1
        ]
        if not positive_scores:
            raise ValueError(
                f"{model['model_id']}: no calibration positives are scorable"
            )
        threshold = sensitivity_threshold(positive_scores, target)
        for split in ("calibration", "evaluation"):
            eligible_split = [
                value for value in observations if value["split"] == split
            ]
            scored = score_observations(eligible_split, fields, references)
            serialized = serialize_scores(scored, model, threshold, analysis_id, split)
            scores_out.extend(serialized)
            summaries.append(
                summarize_model(
                    eligible=eligible_split,
                    scores=serialized,
                    model=model,
                    threshold=threshold,
                    target=target,
                    analysis_id=analysis_id,
                    split=split,
                )
            )
    return references_out, scores_out, summaries


def score_observations(
    observations: list[dict[str, Any]],
    fields: list[str],
    references: dict[tuple[str, ...], float],
) -> list[dict[str, Any]]:
    result = []
    for value in observations:
        key = group_key(value, fields)
        if key not in references:
            continue
        result.append(
            {
                **value,
                "group_key": key,
                "reference_value": references[key],
                "anomaly_score": references[key] - value["value"],
            }
        )
    return result


def serialize_scores(
    scores: list[dict[str, Any]],
    model: dict[str, Any],
    threshold: float,
    analysis_id: str,
    split: str,
) -> list[dict[str, Any]]:
    return [
        {
            "analysis_id": analysis_id,
            "split": split,
            "model_id": model["model_id"],
            "model_name": model["name"],
            "observation_id": value["observation_id"],
            "configuration_accession": value["configuration_accession"],
            "label": value["label"],
            "metric_name": value["metric_name"],
            "region_id": value["region_id"],
            "ontology_signature": value["ontology_signature"],
            "group_key_json": runtime.canonical_json(list(value["group_key"])),
            "reference_value": value["reference_value"],
            "observed_value": value["value"],
            "anomaly_score": value["anomaly_score"],
            "threshold": threshold,
            "detected": value["anomaly_score"] >= threshold,
        }
        for value in scores
    ]


def summarize_model(
    *,
    eligible: list[dict[str, Any]],
    scores: list[dict[str, Any]],
    model: dict[str, Any],
    threshold: float,
    target: float,
    analysis_id: str,
    split: str,
) -> dict[str, Any]:
    sensitivity = detection_rate(scores, label=1)
    return {
        "analysis_id": analysis_id,
        "split": split,
        "model_id": model["model_id"],
        "model_name": model["name"],
        "eligible_observations": len(eligible),
        "scored_observations": len(scores),
        "coverage": len(scores) / len(eligible) if eligible else None,
        "independent_configurations": len(
            {value["configuration_accession"] for value in scores}
        ),
        "positive_observations": sum(value["label"] == 1 for value in scores),
        "negative_observations": sum(value["label"] == 0 for value in scores),
        "threshold": threshold,
        "sensitivity": sensitivity,
        "false_positive_rate": detection_rate(scores, label=0),
        "average_precision": average_precision(scores),
        "sensitivity_calibration_error": (
            abs(sensitivity - target) if sensitivity is not None else None
        ),
    }


def compare_primary_models(
    *,
    scores: list[dict[str, Any]],
    protocol: dict[str, Any],
    analysis_id: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    primary_id = protocol["evaluation"]["primary_model"]
    baseline_id = protocol["evaluation"]["baseline_model"]
    by_model = {}
    for model_id in (primary_id, baseline_id):
        by_model[model_id] = {
            value["observation_id"]: value
            for value in scores
            if value["split"] == "evaluation" and value["model_id"] == model_id
        }
    common_ids = sorted(set(by_model[primary_id]) & set(by_model[baseline_id]))
    if not common_ids:
        raise ValueError("primary models have no common evaluation observations")
    primary = [by_model[primary_id][value] for value in common_ids]
    baseline = [by_model[baseline_id][value] for value in common_ids]
    primary_fpr = detection_rate(primary, label=0)
    baseline_fpr = detection_rate(baseline, label=0)
    if primary_fpr is None or baseline_fpr is None:
        raise ValueError("primary comparison has no clean observations")
    difference = primary_fpr - baseline_fpr
    relative = (baseline_fpr - primary_fpr) / baseline_fpr if baseline_fpr > 0 else None
    bootstrap = bootstrap_fpr_difference(
        analysis_id=analysis_id,
        primary=primary,
        baseline=baseline,
        replicates=protocol["statistics"]["bootstrap_replicates"],
        seed=protocol["statistics"]["bootstrap_seed"],
    )
    confidence = protocol["statistics"]["confidence_level"]
    differences = [value["false_positive_rate_difference"] for value in bootstrap]
    alpha = 1 - confidence
    lower = percentile(differences, alpha / 2)
    upper = percentile(differences, 1 - alpha / 2)
    targets = protocol["scientific_targets"]
    target_met = (
        relative is not None
        and relative >= targets["relative_false_positive_rate_reduction_min"]
        and upper < targets["false_positive_rate_difference_upper_ci_max"]
    )
    primary_ap = average_precision(primary)
    baseline_ap = average_precision(baseline)
    comparison = {
        "analysis_id": analysis_id,
        "comparison_id": "M3_vs_M2",
        "primary_model_id": primary_id,
        "baseline_model_id": baseline_id,
        "common_observations": len(common_ids),
        "independent_configurations": len(
            {value["configuration_accession"] for value in primary}
        ),
        "positive_observations": sum(value["label"] == 1 for value in primary),
        "negative_observations": sum(value["label"] == 0 for value in primary),
        "primary_sensitivity": detection_rate(primary, label=1),
        "baseline_sensitivity": detection_rate(baseline, label=1),
        "primary_false_positive_rate": primary_fpr,
        "baseline_false_positive_rate": baseline_fpr,
        "false_positive_rate_difference": difference,
        "relative_false_positive_rate_reduction": relative,
        "primary_average_precision": primary_ap,
        "baseline_average_precision": baseline_ap,
        "average_precision_difference": primary_ap - baseline_ap,
        "bootstrap_ci_lower": lower,
        "bootstrap_ci_upper": upper,
        "confidence_level": confidence,
        "bootstrap_replicates": len(bootstrap),
        "relative_reduction_target": targets[
            "relative_false_positive_rate_reduction_min"
        ],
        "upper_ci_target": targets["false_positive_rate_difference_upper_ci_max"],
        "scientific_target_met": target_met,
    }
    return comparison, bootstrap


def bootstrap_fpr_difference(
    *,
    analysis_id: str,
    primary: list[dict[str, Any]],
    baseline: list[dict[str, Any]],
    replicates: int,
    seed: int,
) -> list[dict[str, Any]]:
    baseline_by_id = {value["observation_id"]: value for value in baseline}
    paired = [(value, baseline_by_id[value["observation_id"]]) for value in primary]
    by_configuration = defaultdict(list)
    for pair in paired:
        by_configuration[pair[0]["configuration_accession"]].append(pair)
    configurations = sorted(by_configuration)
    if not configurations:
        raise ValueError("bootstrap has no configurations")
    rng = random.Random(seed)
    rows = []
    for replicate in range(replicates):
        sampled = [rng.choice(configurations) for _ in configurations]
        primary_sample = []
        baseline_sample = []
        for configuration in sampled:
            for left, right in by_configuration[configuration]:
                primary_sample.append(left)
                baseline_sample.append(right)
        left_fpr = detection_rate(primary_sample, label=0)
        right_fpr = detection_rate(baseline_sample, label=0)
        if left_fpr is None or right_fpr is None:
            raise ValueError("bootstrap replicate has no clean observations")
        rows.append(
            {
                "analysis_id": analysis_id,
                "replicate": replicate,
                "sampled_configurations": ";".join(sampled),
                "primary_false_positive_rate": left_fpr,
                "baseline_false_positive_rate": right_fpr,
                "false_positive_rate_difference": left_fpr - right_fpr,
            }
        )
    return rows


def compare_terms(
    *,
    scores: list[dict[str, Any]],
    observations: list[dict[str, Any]],
    comparison: dict[str, Any],
    protocol: dict[str, Any],
    analysis_id: str,
) -> list[dict[str, Any]]:
    observation_by_id = {value["observation_id"]: value for value in observations}
    primary_id = comparison["primary_model_id"]
    baseline_id = comparison["baseline_model_id"]
    primary = {
        value["observation_id"]: value
        for value in scores
        if value["split"] == "evaluation" and value["model_id"] == primary_id
    }
    baseline = {
        value["observation_id"]: value
        for value in scores
        if value["split"] == "evaluation" and value["model_id"] == baseline_id
    }
    common = set(primary) & set(baseline)
    minimum = protocol["eligibility"][
        "minimum_independent_configurations_per_term_per_split"
    ]
    rows = []
    declared = sorted(
        {
            (value["metric_name"], term)
            for value in observations
            for term in split_terms(value["ontology_terms"])
        }
    )
    for metric_name, term in declared:
        identifiers = [
            identifier
            for identifier in common
            if observation_by_id[identifier]["metric_name"] == metric_name
            and term in split_terms(observation_by_id[identifier]["ontology_terms"])
        ]
        configurations = {
            observation_by_id[value]["configuration_accession"] for value in identifiers
        }
        if len(configurations) < minimum:
            continue
        left = [primary[value] for value in identifiers]
        right = [baseline[value] for value in identifiers]
        rows.append(
            {
                "analysis_id": analysis_id,
                "metric_name": metric_name,
                "ontology_term": term,
                "independent_configurations": len(configurations),
                "positive_observations": sum(value["label"] == 1 for value in left),
                "negative_observations": sum(value["label"] == 0 for value in left),
                "primary_false_positive_rate": detection_rate(left, label=0),
                "baseline_false_positive_rate": detection_rate(right, label=0),
                "false_positive_rate_difference": rate_difference(left, right, label=0),
                "primary_sensitivity": detection_rate(left, label=1),
                "baseline_sensitivity": detection_rate(right, label=1),
                "primary_average_precision": average_precision(left),
                "baseline_average_precision": average_precision(right),
            }
        )
    return rows


def validate_outputs(
    *,
    analysis_id: str,
    observations: list[dict[str, Any]],
    eligible: list[dict[str, Any]],
    coverage: list[dict[str, Any]],
    references: list[dict[str, Any]],
    scores: list[dict[str, Any]],
    summaries: list[dict[str, Any]],
    comparison: dict[str, Any],
    bootstrap: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    errors = []
    if not observations:
        errors.append("no role metric observations were selected")
    if not eligible:
        errors.append("no ontology terms meet the split eligibility rule")
    for split in ("calibration", "evaluation"):
        split_rows = [value for value in eligible if value["split"] == split]
        if not any(value["label"] == 0 for value in split_rows):
            errors.append(f"{split} has no eligible clean observations")
        if not any(value["label"] == 1 for value in split_rows):
            errors.append(f"{split} has no eligible targeted perturbations")
    if len(summaries) != 2 * len(protocol["models"]):
        errors.append("model summaries do not reconcile")
    if not references or not scores:
        errors.append("model references or scores are empty")
    if comparison["common_observations"] <= 0:
        errors.append("primary model comparison is empty")
    if len(bootstrap) != protocol["statistics"]["bootstrap_replicates"]:
        errors.append("bootstrap replicate count does not reconcile")
    if any(value.get("analysis_id") != analysis_id for value in coverage):
        errors.append("term coverage analysis identifiers differ")
    return {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": analysis_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "scientific_target_met": comparison["scientific_target_met"],
        "counts": {
            "observations": len(observations),
            "eligible_observations": len(eligible),
            "calibration_configurations": len(
                {
                    value["configuration_accession"]
                    for value in eligible
                    if value["split"] == "calibration"
                }
            ),
            "evaluation_configurations": len(
                {
                    value["configuration_accession"]
                    for value in eligible
                    if value["split"] == "evaluation"
                }
            ),
            "eligible_terms": sum(value["eligible"] for value in coverage),
            "references": len(references),
            "scores": len(scores),
            "model_summaries": len(summaries),
            "bootstrap_replicates": len(bootstrap),
        },
    }


def group_key(value: dict[str, Any], fields: list[str]) -> tuple[str, ...]:
    return tuple(str(value[field]) for field in fields)


def sensitivity_threshold(scores: list[float], target: float) -> float:
    if not scores:
        raise ValueError("cannot fit a threshold without positive scores")
    ranked = sorted(scores, reverse=True)
    required = math.ceil(target * len(ranked))
    return ranked[required - 1]


def detection_rate(rows: list[dict[str, Any]], *, label: int) -> float | None:
    selected = [value for value in rows if int(value["label"]) == label]
    if not selected:
        return None
    return sum(as_bool(value["detected"]) for value in selected) / len(selected)


def rate_difference(
    primary: list[dict[str, Any]], baseline: list[dict[str, Any]], *, label: int
) -> float | None:
    left = detection_rate(primary, label=label)
    right = detection_rate(baseline, label=label)
    return left - right if left is not None and right is not None else None


def average_precision(rows: list[dict[str, Any]]) -> float:
    positives = sum(int(value["label"]) == 1 for value in rows)
    if positives == 0:
        return 0.0
    grouped = defaultdict(list)
    for value in rows:
        grouped[float(value["anomaly_score"])].append(int(value["label"]))
    true_positives = 0
    false_positives = 0
    previous_recall = 0.0
    result = 0.0
    for score in sorted(grouped, reverse=True):
        labels = grouped[score]
        true_positives += sum(labels)
        false_positives += len(labels) - sum(labels)
        recall = true_positives / positives
        precision = true_positives / (true_positives + false_positives)
        result += (recall - previous_recall) * precision
        previous_recall = recall
    return result


def percentile(values: list[float], probability: float) -> float:
    if not values:
        raise ValueError("cannot compute a percentile of no values")
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def numeric_json(value: str) -> float | None:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError:
        return None
    if isinstance(parsed, bool) or not isinstance(parsed, (int, float)):
        return None
    result = float(parsed)
    return result if math.isfinite(result) else None


def split_terms(value: str) -> set[str]:
    return {item for item in value.split(";") if item}


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value in {"True", "true", "1", 1}:
        return True
    if value in {"False", "false", "0", 0}:
        return False
    raise ValueError(f"invalid boolean value: {value}")


def required_string(value: dict[str, Any], field: str) -> str:
    result = value.get(field)
    if not isinstance(result, str) or not result.strip():
        raise ValueError(f"required string is missing: {field}")
    return result.strip()


def positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
