#!/usr/bin/env python3
"""Apply a frozen detection policy to a locked perturbation execution."""

from __future__ import annotations

import argparse
import math
import random
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

try:
    import analyze_perturbation_calibration as calibration
    import paper_runtime as runtime
    import run_perturbation_calibration as execution_runner
except ModuleNotFoundError:
    from scripts import analyze_perturbation_calibration as calibration
    from scripts import paper_runtime as runtime
    from scripts import run_perturbation_calibration as execution_runner


SCHEMA_VERSION = "0.1.0"
EVALUATOR_VERSION = "0.1.0"
ENDPOINT_FIELDS = (
    "evaluation_id",
    "endpoint",
    "independent_configurations",
    "observations",
    "estimate",
    "confidence_level",
    "bootstrap_ci_lower",
    "bootstrap_ci_upper",
    "target_minimum",
    "target_met",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate a locked execution with a frozen detection policy."
    )
    parser.add_argument("--execution-manifest", required=True, type=Path)
    parser.add_argument("--detection-policy", required=True, type=Path)
    parser.add_argument("--analysis-protocol", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = evaluate_perturbation_execution(
            execution_manifest_path=args.execution_manifest.resolve(),
            detection_policy_path=args.detection_policy.resolve(),
            analysis_protocol_path=args.analysis_protocol.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"evaluate_perturbation_execution: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def evaluate_perturbation_execution(
    *,
    execution_manifest_path: Path,
    detection_policy_path: Path,
    analysis_protocol_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    execution_manifest_path = execution_manifest_path.resolve()
    detection_policy_path = detection_policy_path.resolve()
    analysis_protocol_path = analysis_protocol_path.resolve()
    output_root = output_root.resolve()
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    execution, conditions, runs, metrics, assessments = calibration.load_execution(
        execution_manifest_path
    )
    materialization_path = calibration.verified_path(
        execution.get("inputs", {}).get("materialization", {}),
        "materialization manifest",
    )
    materialization, _ = execution_runner.load_materialization(materialization_path)
    if materialization.get("sample_source", {}).get("kind") != "policy_sample_bundle":
        raise ValueError("locked evaluation requires a policy sample bundle")
    policy = load_detection_policy(detection_policy_path)
    if policy["execution_id"] == execution["execution_id"]:
        raise ValueError("evaluation execution cannot be its calibration execution")
    protocol = runtime.load_json(analysis_protocol_path)
    validate_protocol(protocol)
    if policy["calibration"] != protocol["calibration"]:
        raise ValueError("detection policy calibration differs from analysis protocol")
    validate_policy_coverage(policy, conditions)

    tools = {
        "evaluator": runtime.script_identity(
            Path(__file__).resolve(), version=EVALUATOR_VERSION
        ),
        "calibration_analyzer": runtime.script_identity(
            Path(calibration.__file__).resolve(),
            version=calibration.ANALYZER_VERSION,
        ),
        "runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
    }
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "execution_id": execution["execution_id"],
        "execution_sha256": runtime.file_sha256(execution_manifest_path),
        "policy_id": policy["policy_id"],
        "policy_sha256": runtime.file_sha256(detection_policy_path),
        "analysis_protocol_sha256": runtime.file_sha256(analysis_protocol_path),
        "evaluator": runtime.functional_script_identity(tools["evaluator"]),
        "calibration_analyzer": runtime.functional_script_identity(
            tools["calibration_analyzer"]
        ),
        "runtime": runtime.functional_script_identity(tools["runtime"]),
    }
    evaluation_id = runtime.sha256_json(stable_identity)[:16]
    run_map = {value["condition_id"]: value for value in runs}
    effects = calibration.build_effect_rows(
        analysis_id=evaluation_id,
        execution_id=execution["execution_id"],
        conditions=conditions,
        metrics=metrics,
    )
    assessment_evidence = calibration.build_assessment_evidence(
        conditions=conditions,
        assessments=assessments,
    )
    calls = calibration.call_conditions(
        analysis_id=evaluation_id,
        execution_id=execution["execution_id"],
        conditions=conditions,
        runs=run_map,
        effects=effects,
        assessment_evidence=assessment_evidence,
        policy=policy,
    )
    summary = calibration.summarize_calls(evaluation_id, calls)
    monotonic = calibration.monotonic_rows(
        analysis_id=evaluation_id,
        effects=effects,
        policy=policy,
        epsilon=protocol["calibration"]["numeric_epsilon"],
    )
    endpoint_rows = evaluation_endpoints(
        evaluation_id=evaluation_id,
        conditions=conditions,
        calls=calls,
        monotonic=monotonic,
        protocol=protocol,
    )
    validation = validate_outputs(
        evaluation_id=evaluation_id,
        conditions=conditions,
        calls=calls,
        effects=effects,
        monotonic=monotonic,
        endpoints=endpoint_rows,
        policy=policy,
    )
    if not validation["valid"]:
        raise ValueError(
            "locked perturbation evaluation failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )

    effects_path = output_root / "tables" / "metric_effects.csv"
    calls_path = output_root / "tables" / "detection_calls.csv"
    summary_path = output_root / "tables" / "operator_summary.csv"
    monotonic_path = output_root / "tables" / "monotonicity.csv"
    endpoints_path = output_root / "tables" / "evaluation_endpoints.csv"
    validation_path = output_root / "validation" / "perturbation_evaluation.json"
    manifest_path = output_root / "manifests" / "evaluation.json"
    runtime.write_csv(effects_path, effects, list(calibration.EFFECT_FIELDS))
    runtime.write_csv(calls_path, calls, list(calibration.CALL_FIELDS))
    runtime.write_csv(summary_path, summary, list(calibration.SUMMARY_FIELDS))
    runtime.write_csv(monotonic_path, monotonic, list(calibration.MONOTONIC_FIELDS))
    runtime.write_csv(endpoints_path, endpoint_rows, list(ENDPOINT_FIELDS))
    runtime.write_json(validation_path, validation)
    manifest = {
        **stable_identity,
        "evaluation_id": evaluation_id,
        "selection_id": execution["selection_id"],
        "created_at": utc_now(),
        "valid": True,
        "inputs": {
            "execution": runtime.file_identity(execution_manifest_path),
            "detection_policy": runtime.file_identity(detection_policy_path),
            "analysis_protocol": runtime.file_identity(analysis_protocol_path),
        },
        "tools": tools,
        "counts": validation["counts"],
        "scientific_targets_met": validation["scientific_targets_met"],
        "outputs": {
            "effects": runtime.file_identity(effects_path),
            "calls": runtime.file_identity(calls_path),
            "summary": runtime.file_identity(summary_path),
            "monotonicity": runtime.file_identity(monotonic_path),
            "endpoints": runtime.file_identity(endpoints_path),
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_detection_policy(path: Path) -> dict[str, Any]:
    policy = runtime.load_json(path)
    if policy.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("detection policy schema is unsupported")
    stable = {
        key: value
        for key, value in policy.items()
        if key not in {"policy_id", "created_at"}
    }
    if runtime.sha256_json(stable)[:16] != policy.get("policy_id"):
        raise ValueError("detection policy id is not content-addressed")
    entries = policy.get("entries")
    if policy.get("frozen") is not True or not isinstance(entries, list) or not entries:
        raise ValueError("detection policy is not frozen")
    keys = [(value.get("operator_id"), value.get("variant")) for value in entries]
    identifiers_invalid = any(
        not isinstance(identifier, str) or not identifier.strip()
        for key in keys
        for identifier in key
    )
    if len(keys) != len(set(keys)) or identifiers_invalid:
        raise ValueError("detection policy entry identifiers are invalid")
    if any(value.get("ready") is not True for value in entries):
        raise ValueError("detection policy contains an unready entry")
    return policy


def validate_protocol(protocol: dict[str, Any]) -> None:
    calibration.validate_protocol(protocol)
    statistics_value = protocol.get("statistics", {})
    confidence = finite_number(statistics_value.get("confidence_level"), "confidence")
    if not 0 < confidence < 1:
        raise ValueError("bootstrap confidence level must be between zero and one")
    positive_int(statistics_value.get("bootstrap_resamples"), "bootstrap resamples")
    seed = statistics_value.get("bootstrap_seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("bootstrap seed must be a nonnegative integer")


def validate_policy_coverage(
    policy: dict[str, Any], conditions: list[dict[str, Any]]
) -> None:
    available = {
        (value["operator_id"], value["variant"]) for value in policy["entries"]
    }
    required = {
        (value["operator_id"], value["variant"])
        for value in conditions
        if value["operator_id"] != "CLEAN"
    }
    missing = sorted(required - available)
    if missing:
        raise ValueError(
            f"detection policy has no rule for evaluation variants: {missing}"
        )


def evaluation_endpoints(
    *,
    evaluation_id: str,
    conditions: list[dict[str, Any]],
    calls: list[dict[str, Any]],
    monotonic: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    targets = protocol["scientific_targets"]
    anchor = protocol["calibration"]["stochastic_anchor_fraction"]
    primary_condition_ids = {
        value["condition_id"]
        for value in conditions
        if value["expected_seqspec_check"] != "failure"
    }
    primary_calls = [
        value for value in calls if value["condition_id"] in primary_condition_ids
    ]
    specifications = (
        (
            "deterministic_sensitivity",
            [
                value
                for value in primary_calls
                if value["operator_id"] != "CLEAN"
                and value["condition_kind"] == "deterministic"
            ],
            lambda value: value["detected"],
            targets["deterministic_sensitivity_min"],
        ),
        (
            "clean_specificity",
            [value for value in primary_calls if value["operator_id"] == "CLEAN"],
            lambda value: not value["detected"],
            targets["clean_specificity_min"],
        ),
        (
            "localization_accuracy",
            [
                value
                for value in primary_calls
                if value["operator_id"] != "CLEAN" and value["detected"]
            ],
            lambda value: value["localized"],
            targets["localization_accuracy_min"],
        ),
        (
            "stochastic_one_percent_sensitivity",
            [
                value
                for value in primary_calls
                if value["condition_kind"] == "stochastic"
                and math.isclose(float(value["event_fraction"]), anchor)
            ],
            lambda value: value["detected"],
            targets["stochastic_one_percent_sensitivity_min"],
        ),
        (
            "monotonic_series_fraction",
            monotonic,
            lambda value: value["monotonic"],
            targets["monotonic_series_fraction_min"],
        ),
    )
    statistics_value = protocol["statistics"]
    result = []
    for index, (name, rows, outcome, target) in enumerate(specifications):
        rates = configuration_rates(rows, outcome)
        estimate = statistics.mean(rates) if rates else None
        lower, upper = bootstrap_mean_interval(
            rates,
            confidence=statistics_value["confidence_level"],
            resamples=statistics_value["bootstrap_resamples"],
            seed=statistics_value["bootstrap_seed"] + index,
        )
        result.append(
            {
                "evaluation_id": evaluation_id,
                "endpoint": name,
                "independent_configurations": len(rates),
                "observations": len(rows),
                "estimate": estimate,
                "confidence_level": statistics_value["confidence_level"],
                "bootstrap_ci_lower": lower,
                "bootstrap_ci_upper": upper,
                "target_minimum": target,
                "target_met": estimate >= target if estimate is not None else None,
            }
        )
    return result


def configuration_rates(
    rows: list[dict[str, Any]],
    outcome: Callable[[dict[str, Any]], bool],
) -> list[float]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["configuration_accession"]].append(bool(outcome(row)))
    return [sum(values) / len(values) for _, values in sorted(grouped.items())]


def bootstrap_mean_interval(
    values: list[float], *, confidence: float, resamples: int, seed: int
) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    randomizer = random.Random(seed)
    size = len(values)
    estimates = sorted(
        statistics.mean(randomizer.choice(values) for _ in range(size))
        for _ in range(resamples)
    )
    alpha = (1 - confidence) / 2
    return percentile(estimates, alpha), percentile(estimates, 1 - alpha)


def percentile(values: list[float], probability: float) -> float:
    position = probability * (len(values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    weight = position - lower
    return values[lower] * (1 - weight) + values[upper] * weight


def validate_outputs(
    *,
    evaluation_id: str,
    conditions: list[dict[str, Any]],
    calls: list[dict[str, Any]],
    effects: list[dict[str, Any]],
    monotonic: list[dict[str, Any]],
    endpoints: list[dict[str, Any]],
    policy: dict[str, Any],
) -> dict[str, Any]:
    errors = []
    if Counter(value["condition_id"] for value in calls) != Counter(
        value["condition_id"] for value in conditions
    ):
        errors.append("evaluation calls do not reconcile with conditions")
    if any(value["analysis_id"] != evaluation_id for value in calls):
        errors.append("evaluation call identifiers differ")
    if len(endpoints) != 5 or len({value["endpoint"] for value in endpoints}) != 5:
        errors.append("scientific endpoints do not reconcile")
    complete = all(value["estimate"] is not None for value in endpoints)
    targets_met = complete and all(value["target_met"] for value in endpoints)
    return {
        "schema_version": SCHEMA_VERSION,
        "evaluation_id": evaluation_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "policy_id": policy["policy_id"],
        "scientific_endpoints_complete": complete,
        "scientific_targets_met": targets_met,
        "counts": {
            "configurations": len(
                {value["configuration_accession"] for value in conditions}
            ),
            "conditions": len(conditions),
            "calls": len(calls),
            "effects": len(effects),
            "monotonic_series": len(monotonic),
            "scientific_endpoints": len(endpoints),
        },
    }


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
