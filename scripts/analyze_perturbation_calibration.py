#!/usr/bin/env python3
"""Calibrate detection endpoints from controlled perturbation executions."""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import paper_runtime as runtime
    import run_perturbation_calibration as execution_runner
except ModuleNotFoundError:
    from scripts import paper_runtime as runtime
    from scripts import run_perturbation_calibration as execution_runner


SCHEMA_VERSION = "0.1.0"
ANALYZER_VERSION = "0.1.0"
EXECUTION_IDENTITY_FIELDS = (
    "schema_version",
    "materialization_id",
    "materialization_sha256",
    "execution_protocol_sha256",
    "seqcheck",
    "runner",
    "runtime",
)
EFFECT_FIELDS = (
    "analysis_id",
    "execution_id",
    "condition_id",
    "configuration_accession",
    "family_id",
    "modality",
    "operator_id",
    "variant",
    "condition_kind",
    "event_fraction",
    "mutation_seed",
    "check",
    "metric_name",
    "unit",
    "reads",
    "regions",
    "clean_value",
    "perturbed_value",
    "signed_effect",
    "absolute_effect",
)
CALL_FIELDS = (
    "analysis_id",
    "execution_id",
    "condition_id",
    "configuration_accession",
    "family_id",
    "modality",
    "operator_id",
    "variant",
    "condition_kind",
    "event_fraction",
    "mutation_seed",
    "detected",
    "detection_source",
    "localized",
    "process_detected",
    "assessment_detected",
    "metric_detected",
    "selected_check",
    "selected_metric_name",
    "selected_unit",
    "selected_direction",
    "selected_threshold",
    "selected_effect",
    "new_assessment_codes",
)
SUMMARY_FIELDS = (
    "analysis_id",
    "operator_id",
    "variant",
    "condition_kind",
    "event_fraction",
    "configurations",
    "conditions",
    "detected",
    "sensitivity",
    "localized_detections",
    "localization_accuracy",
    "end_to_end_localization",
)
ENDPOINT_FIELDS = (
    "analysis_id",
    "operator_id",
    "variant",
    "check",
    "metric_name",
    "unit",
    "direction",
    "configurations",
    "direction_consistency",
    "median_absolute_effect",
    "p10_absolute_effect",
    "threshold",
    "eligible",
    "rank",
)
MONOTONIC_FIELDS = (
    "analysis_id",
    "operator_id",
    "variant",
    "configuration_accession",
    "fractions",
    "median_directional_effects",
    "monotonic",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calibrate controlled perturbation detection endpoints."
    )
    parser.add_argument("--execution-manifest", required=True, type=Path)
    parser.add_argument("--analysis-protocol", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = analyze_perturbation_calibration(
            execution_manifest_path=args.execution_manifest.resolve(),
            analysis_protocol_path=args.analysis_protocol.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"analyze_perturbation_calibration: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def analyze_perturbation_calibration(
    *,
    execution_manifest_path: Path,
    analysis_protocol_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    execution_manifest_path = execution_manifest_path.resolve()
    analysis_protocol_path = analysis_protocol_path.resolve()
    output_root = output_root.resolve()
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    execution, conditions, runs, metrics, assessments = load_execution(
        execution_manifest_path
    )
    protocol = runtime.load_json(analysis_protocol_path)
    validate_protocol(protocol)
    tools = {
        "analyzer": runtime.script_identity(
            Path(__file__).resolve(), version=ANALYZER_VERSION
        ),
        "runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
    }
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "execution_id": execution["execution_id"],
        "execution_sha256": runtime.file_sha256(execution_manifest_path),
        "analysis_protocol_sha256": runtime.file_sha256(analysis_protocol_path),
        "analyzer": runtime.functional_script_identity(tools["analyzer"]),
        "runtime": runtime.functional_script_identity(tools["runtime"]),
    }
    analysis_id = runtime.sha256_json(stable_identity)[:16]
    run_map = {value["condition_id"]: value for value in runs}
    effects = build_effect_rows(
        analysis_id=analysis_id,
        execution_id=execution["execution_id"],
        conditions=conditions,
        metrics=metrics,
    )
    assessment_evidence = build_assessment_evidence(
        conditions=conditions,
        assessments=assessments,
    )
    endpoints = calibrate_metric_endpoints(
        analysis_id=analysis_id,
        conditions=conditions,
        effects=effects,
        protocol=protocol,
    )
    policy = build_policy(
        analysis_id=analysis_id,
        execution=execution,
        conditions=conditions,
        runs=run_map,
        endpoints=endpoints,
        assessment_evidence=assessment_evidence,
        protocol=protocol,
    )
    calls = call_conditions(
        analysis_id=analysis_id,
        execution_id=execution["execution_id"],
        conditions=conditions,
        runs=run_map,
        effects=effects,
        assessment_evidence=assessment_evidence,
        policy=policy,
    )
    summary = summarize_calls(analysis_id, calls)
    monotonic = monotonic_rows(
        analysis_id=analysis_id,
        effects=effects,
        policy=policy,
        epsilon=protocol["calibration"]["numeric_epsilon"],
    )
    validation = validate_outputs(
        analysis_id=analysis_id,
        conditions=conditions,
        calls=calls,
        effects=effects,
        endpoints=endpoints,
        monotonic=monotonic,
        policy=policy,
        protocol=protocol,
    )
    output_root.mkdir(parents=True)
    effects_path = output_root / "tables" / "metric_effects.csv"
    calls_path = output_root / "tables" / "detection_calls.csv"
    summary_path = output_root / "tables" / "operator_summary.csv"
    endpoints_path = output_root / "tables" / "endpoint_candidates.csv"
    monotonic_path = output_root / "tables" / "monotonicity.csv"
    policy_path = output_root / "policy" / "detection_policy.json"
    validation_path = output_root / "validation" / "perturbation_analysis.json"
    runtime.write_csv(effects_path, effects, list(EFFECT_FIELDS))
    runtime.write_csv(calls_path, calls, list(CALL_FIELDS))
    runtime.write_csv(summary_path, summary, list(SUMMARY_FIELDS))
    runtime.write_csv(endpoints_path, endpoints, list(ENDPOINT_FIELDS))
    runtime.write_csv(monotonic_path, monotonic, list(MONOTONIC_FIELDS))
    runtime.write_json(policy_path, policy)
    runtime.write_json(validation_path, validation)
    if not validation["valid"]:
        raise ValueError(
            "perturbation analysis failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )
    manifest_path = output_root / "manifests" / "analysis.json"
    manifest = {
        **stable_identity,
        "analysis_id": analysis_id,
        "selection_id": execution["selection_id"],
        "created_at": utc_now(),
        "valid": True,
        "inputs": {
            "execution": runtime.file_identity(execution_manifest_path),
            "analysis_protocol": runtime.file_identity(analysis_protocol_path),
        },
        "tools": tools,
        "counts": validation["counts"],
        "outputs": {
            "effects": runtime.file_identity(effects_path),
            "calls": runtime.file_identity(calls_path),
            "summary": runtime.file_identity(summary_path),
            "endpoints": runtime.file_identity(endpoints_path),
            "monotonicity": runtime.file_identity(monotonic_path),
            "policy": runtime.file_identity(policy_path),
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_execution(
    path: Path,
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    execution = runtime.load_json(path)
    if execution.get("valid") is not True:
        raise ValueError("perturbation execution is not valid")
    stable = {field: execution.get(field) for field in EXECUTION_IDENTITY_FIELDS}
    if any(stable[field] is None for field in EXECUTION_IDENTITY_FIELDS):
        raise ValueError("execution identity is incomplete")
    if runtime.sha256_json(stable)[:16] != execution.get("execution_id"):
        raise ValueError("execution id is not content-addressed")
    inputs = execution.get("inputs", {})
    materialization_path = verified_path(
        inputs.get("materialization", {}),
        "materialization manifest",
    )
    if runtime.file_sha256(materialization_path) != execution["materialization_sha256"]:
        raise ValueError("execution materialization hash differs")
    execution_protocol_path = verified_path(
        inputs.get("execution_protocol", {}),
        "execution protocol",
    )
    if (
        runtime.file_sha256(execution_protocol_path)
        != execution["execution_protocol_sha256"]
    ):
        raise ValueError("execution protocol hash differs")
    materialization, conditions = execution_runner.load_materialization(
        materialization_path
    )
    if materialization["materialization_id"] != execution["materialization_id"]:
        raise ValueError("execution materialization id differs")
    runs = runtime.read_csv(
        verified_path(execution.get("outputs", {}).get("runs", {}), "run table")
    )
    metrics = runtime.read_csv(
        verified_path(execution.get("outputs", {}).get("metrics", {}), "metric table")
    )
    assessments = runtime.read_csv(
        verified_path(
            execution.get("outputs", {}).get("assessments", {}),
            "assessment table",
        )
    )
    condition_ids = Counter(value["condition_id"] for value in conditions)
    run_ids = Counter(value.get("condition_id", "") for value in runs)
    if run_ids != condition_ids:
        raise ValueError("execution run condition identifiers do not reconcile")
    valid_condition_ids = set(condition_ids)
    for label, rows in (
        ("run", runs),
        ("metric", metrics),
        ("assessment", assessments),
    ):
        if any(
            value.get("execution_id") != execution["execution_id"] for value in rows
        ):
            raise ValueError(f"execution {label} identifiers differ")
        if any(value.get("condition_id") not in valid_condition_ids for value in rows):
            raise ValueError(f"execution {label} condition identifiers differ")
    return execution, conditions, runs, metrics, assessments


def validate_protocol(protocol: dict[str, Any]) -> None:
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("perturbation analysis schema is unsupported")
    calibration = protocol.get("calibration", {})
    positive_int(
        calibration.get("minimum_independent_configurations_per_variant"),
        "minimum configurations",
    )
    anchor = calibration.get("stochastic_anchor_fraction")
    if (
        isinstance(anchor, bool)
        or not isinstance(anchor, (int, float))
        or not 0 < anchor < 1
    ):
        raise ValueError("stochastic anchor fraction is invalid")
    for field in (
        "minimum_direction_consistency",
        "minimum_assessment_support",
        "threshold_fraction_of_p10_effect",
    ):
        value = calibration.get(field)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not 0 < value <= 1
        ):
            raise ValueError(f"analysis calibration value is invalid: {field}")
    epsilon = calibration.get("numeric_epsilon")
    if (
        isinstance(epsilon, bool)
        or not isinstance(epsilon, (int, float))
        or epsilon <= 0
    ):
        raise ValueError("numeric epsilon is invalid")
    floors = calibration.get("unit_threshold_floors")
    if not isinstance(floors, dict) or any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0
        for value in floors.values()
    ):
        raise ValueError("unit threshold floors are invalid")
    targets = protocol.get("scientific_targets", {})
    for field in (
        "deterministic_sensitivity_min",
        "clean_specificity_min",
        "localization_accuracy_min",
        "stochastic_one_percent_sensitivity_min",
        "monotonic_series_fraction_min",
    ):
        value = targets.get(field)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not 0 < value <= 1
        ):
            raise ValueError(f"scientific target is invalid: {field}")


def build_effect_rows(
    *,
    analysis_id: str,
    execution_id: str,
    conditions: list[dict[str, Any]],
    metrics: list[dict[str, str]],
) -> list[dict[str, Any]]:
    clean_by_configuration = {
        value["configuration_accession"]: value["condition_id"]
        for value in conditions
        if value["condition_kind"] == "clean"
    }
    metric_maps = defaultdict(dict)
    for row in metrics:
        if row.get("metric_side") != "observed":
            continue
        value = numeric_json(row.get("value_json", ""))
        if value is None:
            continue
        key = metric_pair_key(row)
        condition_id = row["condition_id"]
        if key in metric_maps[condition_id]:
            raise ValueError(f"{condition_id}: duplicate numeric metric key {key}")
        metric_maps[condition_id][key] = (value, row)
    result = []
    for condition in conditions:
        if condition["condition_kind"] == "clean":
            continue
        clean_id = clean_by_configuration.get(condition["configuration_accession"])
        if clean_id is None:
            raise ValueError(
                f"{condition['configuration_accession']}: clean condition is missing"
            )
        expected_metrics = set(condition["expected"].get("metric_names", []))
        for key, (perturbed_value, row) in metric_maps[
            condition["condition_id"]
        ].items():
            if row["metric_name"] not in expected_metrics:
                continue
            clean = metric_maps[clean_id].get(key)
            if clean is None or not target_scope_matches(row, condition):
                continue
            clean_value = clean[0]
            signed = perturbed_value - clean_value
            result.append(
                {
                    "analysis_id": analysis_id,
                    "execution_id": execution_id,
                    "condition_id": condition["condition_id"],
                    "configuration_accession": condition["configuration_accession"],
                    "family_id": condition["family_id"],
                    "modality": condition["modality"],
                    "operator_id": condition["operator_id"],
                    "variant": condition["variant"],
                    "condition_kind": condition["condition_kind"],
                    "event_fraction": condition["event_fraction"],
                    "mutation_seed": condition["mutation_seed"],
                    "check": row["check"],
                    "metric_name": row["metric_name"],
                    "unit": row["unit"],
                    "reads": row["reads"],
                    "regions": row["regions"],
                    "clean_value": clean_value,
                    "perturbed_value": perturbed_value,
                    "signed_effect": signed,
                    "absolute_effect": abs(signed),
                }
            )
    return sorted(
        result,
        key=lambda value: (
            value["operator_id"],
            value["variant"],
            value["event_fraction"] if value["event_fraction"] is not None else -1,
            value["mutation_seed"] if value["mutation_seed"] is not None else -1,
            value["configuration_accession"],
            value["check"],
            value["metric_name"],
            value["regions"],
        ),
    )


def build_assessment_evidence(
    *,
    conditions: list[dict[str, Any]],
    assessments: list[dict[str, str]],
) -> dict[str, dict[str, Any]]:
    clean_by_configuration = {
        value["configuration_accession"]: value["condition_id"]
        for value in conditions
        if value["condition_kind"] == "clean"
    }
    rows_by_condition = defaultdict(list)
    for row in assessments:
        rows_by_condition[row["condition_id"]].append(row)
    result = {}
    for condition in conditions:
        expected = set(condition["expected"].get("assessment_codes", []))
        if not expected or condition["condition_kind"] == "clean":
            result[condition["condition_id"]] = {
                "new_codes": [],
                "localized_codes": [],
            }
            continue
        current = {
            row["assessment_code"]
            for row in rows_by_condition[condition["condition_id"]]
            if row["assessment_code"] in expected
            and target_scope_matches(row, condition)
        }
        clean_id = clean_by_configuration[condition["configuration_accession"]]
        baseline = {
            row["assessment_code"]
            for row in rows_by_condition[clean_id]
            if row["assessment_code"] in expected
            and target_scope_matches(row, condition)
        }
        new = sorted(current - baseline)
        result[condition["condition_id"]] = {
            "new_codes": new,
            "localized_codes": new,
        }
    return result


def calibrate_metric_endpoints(
    *,
    analysis_id: str,
    conditions: list[dict[str, Any]],
    effects: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    calibration = protocol["calibration"]
    anchor = calibration["stochastic_anchor_fraction"]
    minimum = calibration["minimum_independent_configurations_per_variant"]
    consistency_min = calibration["minimum_direction_consistency"]
    epsilon = calibration["numeric_epsilon"]
    kind_by_variant = {
        (value["operator_id"], value["variant"]): value["condition_kind"]
        for value in conditions
        if value["operator_id"] != "CLEAN"
    }
    grouped = defaultdict(lambda: defaultdict(list))
    for row in effects:
        key = (row["operator_id"], row["variant"])
        if kind_by_variant[key] == "stochastic" and not math.isclose(
            float(row["event_fraction"]), anchor
        ):
            continue
        selector = (row["check"], row["metric_name"], row["unit"])
        grouped[(key, selector)][row["configuration_accession"]].append(
            row["signed_effect"]
        )
    rows = []
    by_variant = defaultdict(list)
    for (variant_key, selector), by_configuration in grouped.items():
        values = [statistics.median(items) for items in by_configuration.values()]
        nonzero = [value for value in values if abs(value) > epsilon]
        positive = sum(value > 0 for value in nonzero)
        negative = sum(value < 0 for value in nonzero)
        direction = "increase" if positive >= negative else "decrease"
        consistency = max(positive, negative) / len(values) if values else 0.0
        absolute = sorted(abs(value) for value in values)
        p10 = lower_quantile(absolute, 0.1) if absolute else 0.0
        median = statistics.median(absolute) if absolute else 0.0
        floor = unit_floor(calibration["unit_threshold_floors"], selector[2])
        threshold = max(
            floor,
            p10 * calibration["threshold_fraction_of_p10_effect"],
        )
        eligible = (
            len(values) >= minimum and consistency >= consistency_min and p10 > epsilon
        )
        row = {
            "analysis_id": analysis_id,
            "operator_id": variant_key[0],
            "variant": variant_key[1],
            "check": selector[0],
            "metric_name": selector[1],
            "unit": selector[2],
            "direction": direction,
            "configurations": len(values),
            "direction_consistency": consistency,
            "median_absolute_effect": median,
            "p10_absolute_effect": p10,
            "threshold": threshold,
            "eligible": eligible,
            "rank": 0,
        }
        by_variant[variant_key].append(row)
    for values in by_variant.values():
        ordered = sorted(
            values,
            key=lambda value: (
                not value["eligible"],
                -value["configurations"],
                -value["direction_consistency"],
                -value["p10_absolute_effect"],
                -value["median_absolute_effect"],
                value["check"],
                value["metric_name"],
                value["unit"],
            ),
        )
        for rank, row in enumerate(ordered, 1):
            row["rank"] = rank
            rows.append(row)
    return sorted(
        rows, key=lambda value: (value["operator_id"], value["variant"], value["rank"])
    )


def build_policy(
    *,
    analysis_id: str,
    execution: dict[str, Any],
    conditions: list[dict[str, Any]],
    runs: dict[str, dict[str, str]],
    endpoints: list[dict[str, Any]],
    assessment_evidence: dict[str, dict[str, Any]],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    calibration = protocol["calibration"]
    minimum = calibration["minimum_independent_configurations_per_variant"]
    support_min = calibration["minimum_assessment_support"]
    best_endpoints = {
        (value["operator_id"], value["variant"]): value
        for value in endpoints
        if value["rank"] == 1 and value["eligible"]
    }
    grouped = defaultdict(list)
    for condition in conditions:
        if condition["operator_id"] != "CLEAN":
            grouped[(condition["operator_id"], condition["variant"])].append(condition)
    entries = []
    for key, values in sorted(grouped.items()):
        kind = values[0]["condition_kind"]
        selected = policy_anchor_conditions(values, calibration)
        configurations = {value["configuration_accession"] for value in selected}
        declared_process_failure = all(
            value["expected"]["process_outcome"] == "failure" for value in values
        )
        stderr_patterns = sorted(
            {
                pattern
                for value in values
                for pattern in value["expected"].get("stderr_patterns", [])
            }
        )
        process_rule = (
            declared_process_failure
            and bool(stderr_patterns)
            and all(
                process_failure_evidence(
                    runs[value["condition_id"]], value, stderr_patterns
                )[0]
                for value in selected
            )
        )
        selected_by_configuration = defaultdict(list)
        for condition in selected:
            selected_by_configuration[condition["configuration_accession"]].append(
                condition
            )
        all_codes = {
            code
            for condition in selected
            for code in assessment_evidence[condition["condition_id"]]["new_codes"]
        }
        code_configurations = defaultdict(set)
        for code in all_codes:
            for (
                configuration,
                configuration_conditions,
            ) in selected_by_configuration.items():
                supported = sum(
                    code in assessment_evidence[value["condition_id"]]["new_codes"]
                    for value in configuration_conditions
                )
                if supported / len(configuration_conditions) >= support_min:
                    code_configurations[code].add(configuration)
        assessment_codes = sorted(
            code
            for code, supported in code_configurations.items()
            if len(configurations) >= minimum
            and len(supported) / len(configurations) >= support_min
        )
        endpoint = best_endpoints.get(key)
        ready = len(configurations) >= minimum and bool(
            process_rule or assessment_codes or endpoint
        )
        entries.append(
            {
                "operator_id": key[0],
                "variant": key[1],
                "condition_kind": kind,
                "calibration_configurations": len(configurations),
                "process_failure": process_rule,
                "stderr_patterns": stderr_patterns,
                "assessment_codes": assessment_codes,
                "metric_endpoint": (
                    {
                        field: endpoint[field]
                        for field in (
                            "check",
                            "metric_name",
                            "unit",
                            "direction",
                            "threshold",
                            "direction_consistency",
                            "p10_absolute_effect",
                        )
                    }
                    if endpoint
                    else None
                ),
                "ready": ready,
            }
        )
    frozen = bool(entries) and all(value["ready"] for value in entries)
    stable = {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": analysis_id,
        "execution_id": execution["execution_id"],
        "frozen": frozen,
        "calibration": protocol["calibration"],
        "entries": entries,
    }
    return {
        **stable,
        "policy_id": runtime.sha256_json(stable)[:16],
        "created_at": utc_now(),
    }


def call_conditions(
    *,
    analysis_id: str,
    execution_id: str,
    conditions: list[dict[str, Any]],
    runs: dict[str, dict[str, str]],
    effects: list[dict[str, Any]],
    assessment_evidence: dict[str, dict[str, Any]],
    policy: dict[str, Any],
) -> list[dict[str, Any]]:
    effects_by_condition = defaultdict(list)
    for row in effects:
        effects_by_condition[row["condition_id"]].append(row)
    policy_map = {
        (value["operator_id"], value["variant"]): value for value in policy["entries"]
    }
    result = []
    for condition in conditions:
        if condition["operator_id"] == "CLEAN":
            result.append(
                call_row(
                    analysis_id,
                    execution_id,
                    condition,
                    detected=False,
                    source="",
                    localized=False,
                )
            )
            continue
        entry = policy_map[(condition["operator_id"], condition["variant"])]
        process_detected, process_localized = process_failure_evidence(
            runs[condition["condition_id"]],
            condition,
            entry["stderr_patterns"] if entry["process_failure"] else [],
        )
        assessment_codes = set(
            assessment_evidence[condition["condition_id"]]["new_codes"]
        )
        assessment_detected = bool(assessment_codes & set(entry["assessment_codes"]))
        endpoint = entry["metric_endpoint"]
        selected_effect = None
        metric_detected = False
        if endpoint:
            candidates = [
                row
                for row in effects_by_condition[condition["condition_id"]]
                if row["check"] == endpoint["check"]
                and row["metric_name"] == endpoint["metric_name"]
                and row["unit"] == endpoint["unit"]
            ]
            if candidates:
                selected_effect = max(
                    candidates,
                    key=lambda value: directional_effect(
                        value["signed_effect"], endpoint["direction"]
                    ),
                )["signed_effect"]
                metric_detected = (
                    directional_effect(selected_effect, endpoint["direction"])
                    >= endpoint["threshold"]
                )
        detected = process_detected or assessment_detected or metric_detected
        sources = []
        if process_detected:
            sources.append("process")
        if assessment_detected:
            sources.append("assessment")
        if metric_detected:
            sources.append("metric")
        localized_assessment_codes = set(
            assessment_evidence[condition["condition_id"]]["localized_codes"]
        )
        assessment_localized = bool(
            localized_assessment_codes & set(entry["assessment_codes"])
        )
        localized = bool(
            detected and (process_localized or assessment_localized or metric_detected)
        )
        result.append(
            call_row(
                analysis_id,
                execution_id,
                condition,
                detected=detected,
                source=";".join(sources),
                localized=localized,
                process_detected=process_detected,
                assessment_detected=assessment_detected,
                metric_detected=metric_detected,
                endpoint=endpoint,
                selected_effect=selected_effect,
                new_codes=sorted(assessment_codes),
            )
        )
    return result


def call_row(
    analysis_id: str,
    execution_id: str,
    condition: dict[str, Any],
    *,
    detected: bool,
    source: str,
    localized: bool,
    process_detected: bool = False,
    assessment_detected: bool = False,
    metric_detected: bool = False,
    endpoint: dict[str, Any] | None = None,
    selected_effect: float | None = None,
    new_codes: list[str] | None = None,
) -> dict[str, Any]:
    endpoint = endpoint or {}
    return {
        "analysis_id": analysis_id,
        "execution_id": execution_id,
        "condition_id": condition["condition_id"],
        "configuration_accession": condition["configuration_accession"],
        "family_id": condition["family_id"],
        "modality": condition["modality"],
        "operator_id": condition["operator_id"],
        "variant": condition["variant"],
        "condition_kind": condition["condition_kind"],
        "event_fraction": condition["event_fraction"],
        "mutation_seed": condition["mutation_seed"],
        "detected": detected,
        "detection_source": source,
        "localized": localized,
        "process_detected": process_detected,
        "assessment_detected": assessment_detected,
        "metric_detected": metric_detected,
        "selected_check": endpoint.get("check", ""),
        "selected_metric_name": endpoint.get("metric_name", ""),
        "selected_unit": endpoint.get("unit", ""),
        "selected_direction": endpoint.get("direction", ""),
        "selected_threshold": endpoint.get("threshold"),
        "selected_effect": selected_effect,
        "new_assessment_codes": ";".join(new_codes or []),
    }


def summarize_calls(
    analysis_id: str, calls: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    grouped = defaultdict(list)
    for row in calls:
        key = (
            row["operator_id"],
            row["variant"],
            row["condition_kind"],
            row["event_fraction"],
        )
        grouped[key].append(row)
    result = []
    for key, values in sorted(
        grouped.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            item[0][3] if item[0][3] is not None else -1,
        ),
    ):
        detected = sum(value["detected"] for value in values)
        localized = sum(value["detected"] and value["localized"] for value in values)
        result.append(
            {
                "analysis_id": analysis_id,
                "operator_id": key[0],
                "variant": key[1],
                "condition_kind": key[2],
                "event_fraction": key[3],
                "configurations": len(
                    {value["configuration_accession"] for value in values}
                ),
                "conditions": len(values),
                "detected": detected,
                "sensitivity": detected / len(values),
                "localized_detections": localized,
                "localization_accuracy": localized / detected if detected else None,
                "end_to_end_localization": localized / len(values),
            }
        )
    return result


def monotonic_rows(
    *,
    analysis_id: str,
    effects: list[dict[str, Any]],
    policy: dict[str, Any],
    epsilon: float,
) -> list[dict[str, Any]]:
    policy_map = {
        (value["operator_id"], value["variant"]): value
        for value in policy["entries"]
        if value["condition_kind"] == "stochastic" and value["metric_endpoint"]
    }
    grouped = defaultdict(lambda: defaultdict(list))
    for row in effects:
        key = (row["operator_id"], row["variant"])
        entry = policy_map.get(key)
        if not entry:
            continue
        endpoint = entry["metric_endpoint"]
        if (
            row["check"] == endpoint["check"]
            and row["metric_name"] == endpoint["metric_name"]
            and row["unit"] == endpoint["unit"]
        ):
            grouped[(key, row["configuration_accession"])][
                float(row["event_fraction"])
            ].append(directional_effect(row["signed_effect"], endpoint["direction"]))
    result = []
    for (key, configuration), by_fraction in sorted(grouped.items()):
        fractions = sorted(by_fraction)
        values = [statistics.median(by_fraction[value]) for value in fractions]
        monotonic = all(
            right + epsilon >= left for left, right in zip(values, values[1:])
        )
        result.append(
            {
                "analysis_id": analysis_id,
                "operator_id": key[0],
                "variant": key[1],
                "configuration_accession": configuration,
                "fractions": ";".join(str(value) for value in fractions),
                "median_directional_effects": ";".join(str(value) for value in values),
                "monotonic": monotonic,
            }
        )
    return result


def validate_outputs(
    *,
    analysis_id: str,
    conditions: list[dict[str, Any]],
    calls: list[dict[str, Any]],
    effects: list[dict[str, Any]],
    endpoints: list[dict[str, Any]],
    monotonic: list[dict[str, Any]],
    policy: dict[str, Any],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    errors = []
    if Counter(value["condition_id"] for value in calls) != Counter(
        value["condition_id"] for value in conditions
    ):
        errors.append("detection calls do not reconcile with conditions")
    if any(value["analysis_id"] != analysis_id for value in calls):
        errors.append("detection call analysis ids differ")
    stable_policy = {
        key: value
        for key, value in policy.items()
        if key not in {"policy_id", "created_at"}
    }
    if runtime.sha256_json(stable_policy)[:16] != policy["policy_id"]:
        errors.append("detection policy id is not content-addressed")
    eligible = [value for value in endpoints if value["eligible"]]
    scientific = scientific_diagnostics(
        calls=calls,
        monotonic=monotonic,
        protocol=protocol,
        policy_frozen=policy["frozen"],
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": analysis_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "policy_frozen": policy["frozen"],
        "scientific_diagnostics": scientific,
        "counts": {
            "conditions": len(conditions),
            "calls": len(calls),
            "effects": len(effects),
            "endpoint_candidates": len(endpoints),
            "eligible_endpoints": len(eligible),
            "policy_entries": len(policy["entries"]),
            "ready_policy_entries": sum(value["ready"] for value in policy["entries"]),
            "minimum_configurations_required": protocol["calibration"][
                "minimum_independent_configurations_per_variant"
            ],
        },
    }


def scientific_diagnostics(
    *,
    calls: list[dict[str, Any]],
    monotonic: list[dict[str, Any]],
    protocol: dict[str, Any],
    policy_frozen: bool,
) -> dict[str, Any]:
    targets = protocol["scientific_targets"]
    anchor = protocol["calibration"]["stochastic_anchor_fraction"]
    deterministic = [
        value
        for value in calls
        if value["operator_id"] != "CLEAN"
        and value["condition_kind"] == "deterministic"
    ]
    clean = [value for value in calls if value["operator_id"] == "CLEAN"]
    detected = [
        value
        for value in calls
        if value["operator_id"] != "CLEAN" and value["detected"]
    ]
    stochastic_anchor = [
        value
        for value in calls
        if value["condition_kind"] == "stochastic"
        and math.isclose(float(value["event_fraction"]), anchor)
    ]
    measures = {
        "deterministic_sensitivity": diagnostic_measure(
            deterministic,
            outcome=lambda value: value["detected"],
            target=targets["deterministic_sensitivity_min"],
        ),
        "clean_specificity": diagnostic_measure(
            clean,
            outcome=lambda value: not value["detected"],
            target=targets["clean_specificity_min"],
        ),
        "localization_accuracy": diagnostic_measure(
            detected,
            outcome=lambda value: value["localized"],
            target=targets["localization_accuracy_min"],
        ),
        "stochastic_one_percent_sensitivity": diagnostic_measure(
            stochastic_anchor,
            outcome=lambda value: value["detected"],
            target=targets["stochastic_one_percent_sensitivity_min"],
        ),
        "monotonic_series_fraction": diagnostic_measure(
            monotonic,
            outcome=lambda value: value["monotonic"],
            target=targets["monotonic_series_fraction_min"],
        ),
    }
    complete = all(value["value"] is not None for value in measures.values())
    targets_met = complete and all(value["met"] for value in measures.values())
    if not policy_frozen:
        status = "policy_not_frozen"
    elif not complete:
        status = "insufficient_applicable_conditions"
    elif targets_met:
        status = "targets_met"
    else:
        status = "targets_not_met"
    return {
        "scope": "calibration",
        "aggregation": "mean_of_configuration_level_rates",
        "status": status,
        "complete": complete,
        "targets_met": targets_met,
        "measures": measures,
    }


def diagnostic_measure(
    rows: list[dict[str, Any]],
    *,
    outcome: Callable[[dict[str, Any]], bool],
    target: float,
) -> dict[str, Any]:
    by_configuration = defaultdict(list)
    for row in rows:
        by_configuration[row["configuration_accession"]].append(bool(outcome(row)))
    rates = [sum(values) / len(values) for values in by_configuration.values()]
    value = statistics.mean(rates) if rates else None
    return {
        "independent_configurations": len(by_configuration),
        "observations": len(rows),
        "value": value,
        "target_minimum": target,
        "met": value >= target if value is not None else None,
    }


def policy_anchor_conditions(
    values: list[dict[str, Any]], calibration: dict[str, Any]
) -> list[dict[str, Any]]:
    if values[0]["condition_kind"] != "stochastic":
        return values
    anchor = calibration["stochastic_anchor_fraction"]
    return [
        value
        for value in values
        if math.isclose(float(value["event_fraction"]), anchor)
    ]


def target_scope_matches(row: dict[str, Any], condition: dict[str, Any]) -> bool:
    localization = set(condition["expected"].get("localization", []))
    target = condition["target"]
    if "read" in localization and target.get("read_ids"):
        if not set(split_values(row.get("reads", ""))) & set(target["read_ids"]):
            return False
    if "region" in localization and target.get("region_ids"):
        if not set(split_values(row.get("regions", ""))) & set(target["region_ids"]):
            return False
    if "file" in localization and target.get("fastq_accessions"):
        files = normalized_file_tokens(split_values(row.get("files", "")))
        targets = normalized_file_tokens(target["fastq_accessions"])
        if not files & targets:
            return False
    return True


def process_failure_evidence(
    run: dict[str, str],
    condition: dict[str, Any],
    stderr_patterns: list[str],
) -> tuple[bool, bool]:
    if run["observed_process_outcome"] != "failure" or not stderr_patterns:
        return False, False
    path = Path(run["stderr_path"]).resolve()
    if runtime.file_sha256(path) != run["stderr_sha256"]:
        raise ValueError(f"{condition['condition_id']}: stderr hash changed")
    stderr = path.read_text(encoding="utf-8", errors="replace").lower()
    detected = any(pattern.lower() in stderr for pattern in stderr_patterns)
    input_names = [Path(value["path"]).name.lower() for value in condition["inputs"]]
    localized = detected and any(value in stderr for value in input_names)
    return detected, localized


def normalized_file_tokens(values: list[str]) -> set[str]:
    result = set()
    for value in values:
        basename = Path(value).name
        result.update((value, basename))
        if basename.endswith(".gz"):
            result.add(basename[:-3])
    return result


def metric_pair_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row.get("check", ""),
        row.get("files", ""),
        row.get("reads", ""),
        row.get("regions", ""),
        row.get("metric_name", ""),
        row.get("unit", ""),
    )


def numeric_json(value: str) -> float | None:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError:
        return None
    if isinstance(parsed, bool) or not isinstance(parsed, (int, float)):
        return None
    result = float(parsed)
    return result if math.isfinite(result) else None


def directional_effect(value: float, direction: str) -> float:
    return value if direction == "increase" else -value


def lower_quantile(values: list[float], probability: float) -> float:
    if not values:
        raise ValueError("cannot compute a quantile of an empty list")
    index = max(0, math.ceil(probability * len(values)) - 1)
    return sorted(values)[index]


def unit_floor(floors: dict[str, Any], unit: str) -> float:
    value = floors.get(unit, floors.get(""))
    if value is None:
        raise ValueError(f"no threshold floor is defined for unit '{unit}'")
    return float(value)


def split_values(value: str) -> list[str]:
    return [item for item in str(value).split(";") if item]


def positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def verified_path(identity: Any, label: str) -> Path:
    if not isinstance(identity, dict):
        raise ValueError(f"{label} identity is missing")
    path = Path(str(identity.get("path", ""))).resolve()
    if runtime.file_sha256(path) != identity.get("sha256"):
        raise ValueError(f"{label} hash changed")
    return path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
