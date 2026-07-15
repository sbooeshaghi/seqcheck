#!/usr/bin/env python3
"""Evaluate locked seqspec, FastQC, and seqcheck complementarity outputs."""

from __future__ import annotations

import argparse
import math
import statistics
import sys
from collections import Counter, defaultdict
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import analyze_complementarity_calibration as quality_calibration
    import analyze_perturbation_calibration as perturbation_calibration
    import evaluate_perturbation_execution as locked_perturbation
    import materialize_complementarity as materializer
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import analyze_complementarity_calibration as quality_calibration
    from scripts import analyze_perturbation_calibration as perturbation_calibration
    from scripts import evaluate_perturbation_execution as locked_perturbation
    from scripts import materialize_complementarity as materializer
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
EVALUATOR_VERSION = "0.1.0"
TOOL_CALL_FIELDS = (
    "evaluation_id",
    "execution_id",
    "condition_source",
    "condition_id",
    "configuration_accession",
    "family_id",
    "modality",
    "condition_kind",
    "operator_id",
    "variant",
    "tool",
    "baseline_condition_id",
    "detected",
    "detection_source",
    "localized",
    "evidence_json",
)
FASTQC_DELTA_FIELDS = (
    "evaluation_id",
    "execution_id",
    "condition_source",
    "condition_id",
    "configuration_accession",
    "condition_kind",
    "operator_id",
    "variant",
    "case_id",
    "module_name",
    "clean_status",
    "observed_status",
    "status_worsened",
    "clean_score",
    "observed_score",
    "score_decrease",
    "in_predeclared_scope",
)
SUMMARY_FIELDS = (
    "evaluation_id",
    "operator_id",
    "variant",
    "condition_kind",
    "tool",
    "independent_configurations",
    "conditions",
    "detected",
    "detection_rate",
    "localized_detections",
    "localization_accuracy",
)
ENDPOINT_FIELDS = (
    "evaluation_id",
    "endpoint",
    "independent_configurations",
    "observations",
    "estimate",
    "confidence_level",
    "bootstrap_ci_lower",
    "bootstrap_ci_upper",
    "target_direction",
    "target",
    "target_met",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate the locked three-tool complementarity experiment."
    )
    parser.add_argument("--execution-manifest", required=True, type=Path)
    parser.add_argument("--quality-policy", required=True, type=Path)
    parser.add_argument("--detection-policy", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = evaluate_complementarity(
            execution_path=args.execution_manifest.resolve(),
            quality_policy_path=args.quality_policy.resolve(),
            detection_policy_path=args.detection_policy.resolve(),
            protocol_path=args.protocol.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"evaluate_complementarity: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def evaluate_complementarity(
    *,
    execution_path: Path,
    quality_policy_path: Path,
    detection_policy_path: Path,
    protocol_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    loaded = quality_calibration.load_execution(execution_path)
    execution = loaded["execution"]
    complementarity = loaded["complementarity"]
    base = loaded["base"]
    base_conditions = loaded["base_conditions"]
    controls = loaded["controls"]
    conditions = [*base_conditions, *controls]
    if execution.get("cohort_split") != "evaluation":
        raise ValueError("complementarity evaluation requires the evaluation split")
    if base.get("sample_source", {}).get("kind") != "policy_sample_bundle":
        raise ValueError("complementarity evaluation requires locked sampled reads")

    protocol = runtime.load_json(protocol_path)
    validate_protocol(protocol)
    protocol_sha256 = runtime.file_sha256(protocol_path)
    if execution["protocol_sha256"] != protocol_sha256:
        raise ValueError("evaluation protocol differs from execution")
    quality_policy = materializer.load_quality_policy(
        quality_policy_path, protocol_path, protocol
    )
    expected_quality_policy = complementarity.get("quality_policy")
    if expected_quality_policy != {
        "policy_id": quality_policy["policy_id"],
        "sha256": runtime.file_sha256(quality_policy_path),
    }:
        raise ValueError("evaluation quality policy differs from materialization")
    detection_policy = locked_perturbation.load_detection_policy(detection_policy_path)
    if detection_policy["execution_id"] == execution["execution_id"]:
        raise ValueError("evaluation cannot reuse its own execution as calibration")
    locked_perturbation.validate_policy_coverage(detection_policy, base_conditions)

    tools = tool_identities()
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "execution_id": execution["execution_id"],
        "execution_sha256": runtime.file_sha256(execution_path),
        "quality_policy_id": quality_policy["policy_id"],
        "quality_policy_sha256": runtime.file_sha256(quality_policy_path),
        "detection_policy_id": detection_policy["policy_id"],
        "detection_policy_sha256": runtime.file_sha256(detection_policy_path),
        "protocol_sha256": protocol_sha256,
        **{
            name: runtime.functional_script_identity(identity)
            for name, identity in tools.items()
        },
    }
    evaluation_id = runtime.sha256_json(stable_identity)[:16]

    seqcheck_run_map = unique_rows(
        loaded["seqcheck_runs"], "condition_id", "seqcheck run"
    )
    seqspec_run_map = unique_rows(loaded["seqspec_runs"], "condition_id", "seqspec run")
    structural_effects = perturbation_calibration.build_effect_rows(
        analysis_id=evaluation_id,
        execution_id=execution["execution_id"],
        conditions=base_conditions,
        metrics=loaded["seqcheck_metrics"],
    )
    structural_assessments = perturbation_calibration.build_assessment_evidence(
        conditions=base_conditions,
        assessments=loaded["seqcheck_assessments"],
    )
    structural_seqcheck_calls = perturbation_calibration.call_conditions(
        analysis_id=evaluation_id,
        execution_id=execution["execution_id"],
        conditions=base_conditions,
        runs={
            value["condition_id"]: seqcheck_run_map[value["condition_id"]]
            for value in base_conditions
        },
        effects=structural_effects,
        assessment_evidence=structural_assessments,
        policy=detection_policy,
    )

    module_map = quality_calibration.build_fastqc_module_map(
        loaded["fastqc_runs"], loaded["fastqc_modules"]
    )
    quality_fastqc_calls = quality_calibration.build_quality_calls(
        analysis_id=evaluation_id,
        execution_id=execution["execution_id"],
        base_conditions=base_conditions,
        controls=controls,
        module_map=module_map,
        protocol=protocol,
    )
    seqcheck_invariance = quality_calibration.build_seqcheck_invariance(
        analysis_id=evaluation_id,
        execution_id=execution["execution_id"],
        base_conditions=base_conditions,
        controls=controls,
        metrics=loaded["seqcheck_metrics"],
        assessments=loaded["seqcheck_assessments"],
    )
    fastqc_deltas = build_fastqc_deltas(
        evaluation_id=evaluation_id,
        execution_id=execution["execution_id"],
        base_conditions=base_conditions,
        controls=controls,
        module_map=module_map,
    )
    tool_calls = build_tool_calls(
        evaluation_id=evaluation_id,
        execution_id=execution["execution_id"],
        base_conditions=base_conditions,
        controls=controls,
        seqspec_runs=seqspec_run_map,
        seqcheck_runs=seqcheck_run_map,
        seqcheck_metrics=loaded["seqcheck_metrics"],
        seqcheck_assessments=loaded["seqcheck_assessments"],
        structural_seqcheck_calls=structural_seqcheck_calls,
        quality_fastqc_calls=quality_fastqc_calls,
        fastqc_deltas=fastqc_deltas,
    )
    summary = summarize_tool_calls(evaluation_id, tool_calls)
    endpoints = build_endpoints(
        evaluation_id=evaluation_id,
        base_conditions=base_conditions,
        controls=controls,
        calls=tool_calls,
        protocol=protocol,
    )
    validation = validate_outputs(
        evaluation_id=evaluation_id,
        conditions=conditions,
        structural_conditions=base_conditions,
        controls=controls,
        tool_calls=tool_calls,
        structural_seqcheck_calls=structural_seqcheck_calls,
        quality_fastqc_calls=quality_fastqc_calls,
        seqcheck_invariance=seqcheck_invariance,
        fastqc_deltas=fastqc_deltas,
        endpoints=endpoints,
    )
    if not validation["valid"]:
        raise ValueError(
            "complementarity evaluation failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )

    tables = output_root / "tables"
    paths = {
        "tool_calls": tables / "tool_calls.csv",
        "tool_summary": tables / "tool_summary.csv",
        "fastqc_deltas": tables / "fastqc_module_deltas.csv",
        "quality_fastqc_calls": tables / "quality_fastqc_calls.csv",
        "seqcheck_invariance": tables / "seqcheck_invariance.csv",
        "structural_effects": tables / "structural_metric_effects.csv",
        "structural_seqcheck_calls": tables / "structural_seqcheck_calls.csv",
        "endpoints": tables / "scientific_endpoints.csv",
    }
    runtime.write_csv(paths["tool_calls"], tool_calls, list(TOOL_CALL_FIELDS))
    runtime.write_csv(paths["tool_summary"], summary, list(SUMMARY_FIELDS))
    runtime.write_csv(paths["fastqc_deltas"], fastqc_deltas, list(FASTQC_DELTA_FIELDS))
    runtime.write_csv(
        paths["quality_fastqc_calls"],
        quality_fastqc_calls,
        list(quality_calibration.CALL_FIELDS),
    )
    runtime.write_csv(
        paths["seqcheck_invariance"],
        seqcheck_invariance,
        list(quality_calibration.INVARIANCE_FIELDS),
    )
    runtime.write_csv(
        paths["structural_effects"],
        structural_effects,
        list(perturbation_calibration.EFFECT_FIELDS),
    )
    runtime.write_csv(
        paths["structural_seqcheck_calls"],
        structural_seqcheck_calls,
        list(perturbation_calibration.CALL_FIELDS),
    )
    runtime.write_csv(paths["endpoints"], endpoints, list(ENDPOINT_FIELDS))
    validation_path = output_root / "validation" / "complementarity_evaluation.json"
    runtime.write_json(validation_path, validation)
    manifest_path = output_root / "manifests" / "evaluation.json"
    manifest = {
        **stable_identity,
        "evaluation_id": evaluation_id,
        "selection_id": execution["selection_id"],
        "created_at": utc_now(),
        "valid": True,
        "inputs": {
            "execution": runtime.file_identity(execution_path),
            "quality_policy": runtime.file_identity(quality_policy_path),
            "detection_policy": runtime.file_identity(detection_policy_path),
            "protocol": runtime.file_identity(protocol_path),
        },
        "tools": tools,
        "counts": validation["counts"],
        "scientific_targets_met": validation["scientific_targets_met"],
        "outputs": {
            **{name: runtime.file_identity(path) for name, path in paths.items()},
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def tool_identities() -> dict[str, dict[str, Any]]:
    return {
        "evaluator": runtime.script_identity(
            Path(__file__).resolve(), version=EVALUATOR_VERSION
        ),
        "quality_calibration": runtime.script_identity(
            Path(quality_calibration.__file__).resolve(),
            version=quality_calibration.ANALYZER_VERSION,
        ),
        "perturbation_calibration": runtime.script_identity(
            Path(perturbation_calibration.__file__).resolve(),
            version=perturbation_calibration.ANALYZER_VERSION,
        ),
        "locked_perturbation": runtime.script_identity(
            Path(locked_perturbation.__file__).resolve(),
            version=locked_perturbation.EVALUATOR_VERSION,
        ),
        "materializer": runtime.script_identity(
            Path(materializer.__file__).resolve(),
            version=materializer.MATERIALIZER_VERSION,
        ),
        "fastqc_parser": runtime.script_identity(
            Path(quality_calibration.fastqc_parser.__file__).resolve(),
            version=quality_calibration.fastqc_parser.RUNTIME_VERSION,
        ),
        "runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
    }


def build_fastqc_deltas(
    *,
    evaluation_id: str,
    execution_id: str,
    base_conditions: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    module_map: dict[tuple[str, str], dict[str, dict[str, Any]]],
) -> list[dict[str, Any]]:
    clean = clean_conditions(base_conditions)
    result = []
    combined = [("structural", value) for value in base_conditions] + [
        ("control", value) for value in controls
    ]
    for source, condition in combined:
        if condition["condition_kind"] == "clean":
            continue
        baseline = clean[condition["configuration_accession"]]
        scoped = set(condition.get("expected", {}).get("fastqc_modules", []))
        for input_value in condition["inputs"]:
            key = (condition["condition_id"], input_value["case_id"])
            clean_key = (baseline["condition_id"], input_value["case_id"])
            if key not in module_map or clean_key not in module_map:
                raise ValueError(
                    f"{condition['condition_id']}: FastQC baseline is missing"
                )
            observed_modules = module_map[key]
            clean_modules = module_map[clean_key]
            if set(observed_modules) != set(clean_modules):
                raise ValueError(
                    f"{condition['condition_id']}: FastQC module names differ from clean"
                )
            for name in sorted(observed_modules):
                clean_module = clean_modules[name]
                observed_module = observed_modules[name]
                clean_score = quality_calibration.module_score(clean_module)
                observed_score = quality_calibration.module_score(observed_module)
                result.append(
                    {
                        "evaluation_id": evaluation_id,
                        "execution_id": execution_id,
                        "condition_source": source,
                        "condition_id": condition["condition_id"],
                        "configuration_accession": condition["configuration_accession"],
                        "condition_kind": condition["condition_kind"],
                        "operator_id": condition["operator_id"],
                        "variant": condition["variant"],
                        "case_id": input_value["case_id"],
                        "module_name": name,
                        "clean_status": clean_module["status"],
                        "observed_status": observed_module["status"],
                        "status_worsened": quality_calibration.fastqc_parser.status_worsened(
                            clean_module["status"], observed_module["status"]
                        ),
                        "clean_score": clean_score,
                        "observed_score": observed_score,
                        "score_decrease": (
                            clean_score - observed_score
                            if clean_score is not None and observed_score is not None
                            else None
                        ),
                        "in_predeclared_scope": not scoped or name in scoped,
                    }
                )
    return result


def build_tool_calls(
    *,
    evaluation_id: str,
    execution_id: str,
    base_conditions: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    seqspec_runs: dict[str, dict[str, str]],
    seqcheck_runs: dict[str, dict[str, str]],
    seqcheck_metrics: list[dict[str, str]],
    seqcheck_assessments: list[dict[str, str]],
    structural_seqcheck_calls: list[dict[str, Any]],
    quality_fastqc_calls: list[dict[str, Any]],
    fastqc_deltas: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    clean = clean_conditions(base_conditions)
    structural_map = {
        value["condition_id"]: value for value in structural_seqcheck_calls
    }
    quality_map = {value["condition_id"]: value for value in quality_fastqc_calls}
    deltas = defaultdict(list)
    for row in fastqc_deltas:
        deltas[row["condition_id"]].append(row)
    metric_map = quality_calibration.canonical_rows_by_condition(seqcheck_metrics)
    assessment_map = quality_calibration.canonical_rows_by_condition(
        seqcheck_assessments
    )
    result = []
    combined = [("structural", value) for value in base_conditions] + [
        ("control", value) for value in controls
    ]
    for source, condition in combined:
        baseline = clean[condition["configuration_accession"]]
        baseline_id = baseline["condition_id"]
        clean_condition = condition["condition_kind"] == "clean"
        seqspec_detected = (
            not clean_condition
            and seqspec_runs[baseline_id]["observed_status"] == "pass"
            and seqspec_runs[condition["condition_id"]]["observed_status"] == "failure"
        )
        result.append(
            tool_call(
                evaluation_id,
                execution_id,
                source,
                condition,
                "seqspec",
                baseline_id,
                seqspec_detected,
                "check_failure" if seqspec_detected else "",
                False,
                {
                    "clean_status": seqspec_runs[baseline_id]["observed_status"],
                    "observed_status": seqspec_runs[condition["condition_id"]][
                        "observed_status"
                    ],
                },
            )
        )

        if source == "structural":
            structural = structural_map[condition["condition_id"]]
            seqcheck_detected = structural["detected"]
            seqcheck_source = structural["detection_source"]
            localized = structural["localized"]
            seqcheck_evidence = {
                "selected_check": structural["selected_check"],
                "selected_metric_name": structural["selected_metric_name"],
                "selected_threshold": structural["selected_threshold"],
                "selected_effect": structural["selected_effect"],
                "new_assessment_codes": structural["new_assessment_codes"],
            }
        else:
            condition_id = condition["condition_id"]
            process_detected = (
                seqcheck_runs[baseline_id]["observed_process_outcome"] == "success"
                and seqcheck_runs[condition_id]["observed_process_outcome"] == "failure"
            )
            metrics_differ = metric_map.get(condition_id, []) != metric_map.get(
                baseline_id, []
            )
            assessments_differ = assessment_map.get(
                condition_id, []
            ) != assessment_map.get(baseline_id, [])
            seqcheck_detected = process_detected or metrics_differ or assessments_differ
            sources = []
            if process_detected:
                sources.append("process")
            if metrics_differ:
                sources.append("metric")
            if assessments_differ:
                sources.append("assessment")
            seqcheck_source = ";".join(sources)
            localized = False
            seqcheck_evidence = {
                "process_detected": process_detected,
                "metrics_differ": metrics_differ,
                "assessments_differ": assessments_differ,
            }
        result.append(
            tool_call(
                evaluation_id,
                execution_id,
                source,
                condition,
                "seqcheck",
                baseline_id,
                seqcheck_detected,
                seqcheck_source,
                localized,
                seqcheck_evidence,
            )
        )

        if clean_condition:
            fastqc_detected = False
            fastqc_source = ""
            fastqc_evidence = {"module_names": []}
        elif condition["condition_kind"] == "quality":
            quality = quality_map[condition["condition_id"]]
            fastqc_detected = quality["detected"]
            sources = []
            if quality["status_worsened"]:
                sources.append("module_status")
            if (
                quality["score_decrease"] is not None
                and quality["score_decrease"] >= quality["score_threshold"]
            ):
                sources.append("raw_quality_score")
            fastqc_source = ";".join(sources)
            fastqc_evidence = {
                "module_names": quality["module_names"],
                "score_decrease": quality["score_decrease"],
                "score_threshold": quality["score_threshold"],
            }
        else:
            worsened = [
                value
                for value in deltas[condition["condition_id"]]
                if value["in_predeclared_scope"] and value["status_worsened"]
            ]
            fastqc_detected = bool(worsened)
            fastqc_source = "module_status" if worsened else ""
            fastqc_evidence = {
                "module_names": sorted({value["module_name"] for value in worsened})
            }
        result.append(
            tool_call(
                evaluation_id,
                execution_id,
                source,
                condition,
                "fastqc",
                baseline_id,
                fastqc_detected,
                fastqc_source,
                False,
                fastqc_evidence,
            )
        )
    return sorted(
        result,
        key=lambda value: (
            value["operator_id"],
            value["variant"],
            value["configuration_accession"],
            value["condition_id"],
            value["tool"],
        ),
    )


def tool_call(
    evaluation_id: str,
    execution_id: str,
    source: str,
    condition: dict[str, Any],
    tool: str,
    baseline_condition_id: str,
    detected: bool,
    detection_source: str,
    localized: bool,
    evidence: dict[str, Any],
) -> dict[str, Any]:
    return {
        "evaluation_id": evaluation_id,
        "execution_id": execution_id,
        "condition_source": source,
        "condition_id": condition["condition_id"],
        "configuration_accession": condition["configuration_accession"],
        "family_id": condition["family_id"],
        "modality": condition["modality"],
        "condition_kind": condition["condition_kind"],
        "operator_id": condition["operator_id"],
        "variant": condition["variant"],
        "tool": tool,
        "baseline_condition_id": baseline_condition_id,
        "detected": bool(detected),
        "detection_source": detection_source,
        "localized": bool(localized),
        "evidence_json": runtime.canonical_json(evidence),
    }


def summarize_tool_calls(
    evaluation_id: str, calls: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    grouped = defaultdict(list)
    for row in calls:
        grouped[
            (
                row["operator_id"],
                row["variant"],
                row["condition_kind"],
                row["tool"],
            )
        ].append(row)
    result = []
    for key, rows in sorted(grouped.items()):
        detected = [value for value in rows if value["detected"]]
        localized = sum(value["localized"] for value in detected)
        configuration_rates = grouped_values(rows, lambda value: value["detected"])
        result.append(
            {
                "evaluation_id": evaluation_id,
                "operator_id": key[0],
                "variant": key[1],
                "condition_kind": key[2],
                "tool": key[3],
                "independent_configurations": len(configuration_rates),
                "conditions": len(rows),
                "detected": len(detected),
                "detection_rate": (
                    statistics.mean(configuration_rates)
                    if configuration_rates
                    else None
                ),
                "localized_detections": localized,
                "localization_accuracy": (
                    localized / len(detected) if detected else None
                ),
            }
        )
    return result


def build_endpoints(
    *,
    evaluation_id: str,
    base_conditions: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    calls: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    call_map = {(value["condition_id"], value["tool"]): value for value in calls}
    primary = primary_structural_conditions(base_conditions, protocol)
    invalid = [
        value for value in controls if value["condition_kind"] == "schema_invalid"
    ]
    quality = [value for value in controls if value["condition_kind"] == "quality"]
    structural_detected = [
        value
        for value in primary
        if call_map[(value["condition_id"], "seqcheck")]["detected"]
    ]
    targets = protocol["scientific_targets"]
    specifications = (
        (
            "seqcheck_structural_advantage",
            primary,
            lambda value: float(
                call_map[(value["condition_id"], "seqcheck")]["detected"]
            )
            - float(call_map[(value["condition_id"], "fastqc")]["detected"]),
            "minimum",
            targets["seqcheck_structural_advantage_min"],
        ),
        (
            "seqcheck_structural_localization",
            structural_detected,
            lambda value: float(
                call_map[(value["condition_id"], "seqcheck")]["localized"]
            ),
            "minimum",
            targets["seqcheck_localization_min"],
        ),
        (
            "seqspec_invalid_rejection",
            invalid,
            lambda value: float(
                call_map[(value["condition_id"], "seqspec")]["detected"]
            ),
            "minimum",
            targets["seqspec_invalid_rejection_min"],
        ),
        (
            "fastqc_quality_sensitivity",
            quality,
            lambda value: float(
                call_map[(value["condition_id"], "fastqc")]["detected"]
            ),
            "minimum",
            targets["fastqc_quality_sensitivity_min"],
        ),
        (
            "seqcheck_quality_difference_rate",
            quality,
            lambda value: float(
                call_map[(value["condition_id"], "seqcheck")]["detected"]
            ),
            "maximum",
            targets["quality_seqcheck_difference_max"],
        ),
    )
    statistics_value = protocol["statistics"]
    result = []
    for index, (name, rows, outcome, direction, target) in enumerate(specifications):
        values = grouped_values(rows, outcome)
        estimate = statistics.mean(values) if values else None
        lower, upper = locked_perturbation.bootstrap_mean_interval(
            values,
            confidence=statistics_value["confidence_level"],
            resamples=statistics_value["bootstrap_resamples"],
            seed=statistics_value["bootstrap_seed"] + index,
        )
        target_met = None
        if estimate is not None:
            target_met = (
                estimate >= target if direction == "minimum" else estimate <= target
            )
        result.append(
            {
                "evaluation_id": evaluation_id,
                "endpoint": name,
                "independent_configurations": len(values),
                "observations": len(rows),
                "estimate": estimate,
                "confidence_level": statistics_value["confidence_level"],
                "bootstrap_ci_lower": lower,
                "bootstrap_ci_upper": upper,
                "target_direction": direction,
                "target": target,
                "target_met": target_met,
            }
        )
    return result


def primary_structural_conditions(
    conditions: list[dict[str, Any]], protocol: dict[str, Any]
) -> list[dict[str, Any]]:
    scope = protocol["primary_structural_scope"]
    anchor = scope["stochastic_event_fraction"]
    return [
        value
        for value in conditions
        if value["condition_kind"] != "clean"
        and value["expected_seqspec_check"] != "failure"
        and (
            (
                scope["include_deterministic"]
                and value["condition_kind"] == "deterministic"
            )
            or (
                value["condition_kind"] == "stochastic"
                and math.isclose(float(value["event_fraction"]), anchor)
            )
        )
    ]


def grouped_values(
    rows: list[dict[str, Any]], outcome: Callable[[dict[str, Any]], float | bool]
) -> list[float]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["configuration_accession"]].append(float(outcome(row)))
    return [statistics.mean(values) for _, values in sorted(grouped.items())]


def validate_protocol(protocol: dict[str, Any]) -> None:
    quality_calibration.validate_protocol(protocol)
    statistics_value = protocol.get("statistics", {})
    confidence = finite_number(statistics_value.get("confidence_level"), "confidence")
    if not 0 < confidence < 1:
        raise ValueError("bootstrap confidence level must be between zero and one")
    positive_int(statistics_value.get("bootstrap_resamples"), "bootstrap resamples")
    seed = statistics_value.get("bootstrap_seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("bootstrap seed must be a nonnegative integer")
    targets = protocol.get("scientific_targets", {})
    for field in (
        "seqcheck_structural_advantage_min",
        "seqcheck_localization_min",
        "seqspec_invalid_rejection_min",
        "fastqc_quality_sensitivity_min",
        "quality_seqcheck_difference_max",
    ):
        value = finite_number(targets.get(field), field)
        if not 0 <= value <= 1:
            raise ValueError(f"scientific target is outside [0, 1]: {field}")


def validate_outputs(
    *,
    evaluation_id: str,
    conditions: list[dict[str, Any]],
    structural_conditions: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    tool_calls: list[dict[str, Any]],
    structural_seqcheck_calls: list[dict[str, Any]],
    quality_fastqc_calls: list[dict[str, Any]],
    seqcheck_invariance: list[dict[str, Any]],
    fastqc_deltas: list[dict[str, Any]],
    endpoints: list[dict[str, Any]],
) -> dict[str, Any]:
    errors = []
    expected_keys = {
        (value["condition_id"], tool)
        for value in conditions
        for tool in ("seqspec", "seqcheck", "fastqc")
    }
    observed_keys = [(value["condition_id"], value["tool"]) for value in tool_calls]
    if set(observed_keys) != expected_keys or len(observed_keys) != len(expected_keys):
        errors.append("tool calls do not reconcile with conditions and tools")
    if any(value["evaluation_id"] != evaluation_id for value in tool_calls):
        errors.append("tool call evaluation identifiers differ")
    if Counter(value["condition_id"] for value in structural_seqcheck_calls) != Counter(
        value["condition_id"] for value in structural_conditions
    ):
        errors.append("structural seqcheck calls do not reconcile")
    quality_ids = Counter(
        value["condition_id"]
        for value in controls
        if value["condition_kind"] == "quality"
    )
    if Counter(value["condition_id"] for value in quality_fastqc_calls) != quality_ids:
        errors.append("quality FastQC calls do not reconcile")
    if Counter(value["condition_id"] for value in seqcheck_invariance) != quality_ids:
        errors.append("seqcheck quality invariance rows do not reconcile")
    call_map = {(value["condition_id"], value["tool"]): value for value in tool_calls}
    if any(
        bool(call_map[(value["condition_id"], "fastqc")]["detected"])
        != bool(value["detected"])
        for value in quality_fastqc_calls
    ):
        errors.append("quality FastQC calls differ from tool calls")
    if any(
        bool(call_map[(value["condition_id"], "seqcheck")]["detected"])
        == bool(value["passed"])
        for value in seqcheck_invariance
    ):
        errors.append("seqcheck invariance differs from quality tool calls")
    if any(
        bool(call_map[(value["condition_id"], "seqcheck")]["detected"])
        != bool(value["detected"])
        for value in structural_seqcheck_calls
    ):
        errors.append("structural seqcheck calls differ from tool calls")
    if not fastqc_deltas:
        errors.append("FastQC module delta evidence is empty")
    if any(value["evaluation_id"] != evaluation_id for value in fastqc_deltas):
        errors.append("FastQC delta evaluation identifiers differ")
    expected_endpoints = {
        "seqcheck_structural_advantage",
        "seqcheck_structural_localization",
        "seqspec_invalid_rejection",
        "fastqc_quality_sensitivity",
        "seqcheck_quality_difference_rate",
    }
    if {value["endpoint"] for value in endpoints} != expected_endpoints:
        errors.append("scientific endpoints do not reconcile")
    if any(value["evaluation_id"] != evaluation_id for value in endpoints):
        errors.append("scientific endpoint evaluation identifiers differ")
    complete = all(value["estimate"] is not None for value in endpoints)
    targets_met = complete and all(value["target_met"] for value in endpoints)
    return {
        "schema_version": SCHEMA_VERSION,
        "evaluation_id": evaluation_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "scientific_endpoints_complete": complete,
        "scientific_targets_met": targets_met,
        "counts": {
            "configurations": len(
                {value["configuration_accession"] for value in conditions}
            ),
            "conditions": len(conditions),
            "tool_calls": len(tool_calls),
            "fastqc_module_deltas": len(fastqc_deltas),
            "quality_fastqc_calls": len(quality_fastqc_calls),
            "seqcheck_invariance_rows": len(seqcheck_invariance),
            "structural_seqcheck_calls": len(structural_seqcheck_calls),
            "scientific_endpoints": len(endpoints),
        },
    }


def clean_conditions(
    conditions: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    grouped = defaultdict(list)
    for value in conditions:
        if value["condition_kind"] == "clean":
            grouped[value["configuration_accession"]].append(value)
    if not grouped or any(len(values) != 1 for values in grouped.values()):
        raise ValueError("structural conditions need one clean row per configuration")
    return {key: values[0] for key, values in grouped.items()}


def unique_rows(
    rows: list[dict[str, str]], key: str, label: str
) -> dict[str, dict[str, str]]:
    result = {value[key]: value for value in rows}
    if len(result) != len(rows):
        raise ValueError(f"{label} identifiers are not unique")
    return result


def finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
