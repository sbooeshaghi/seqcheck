#!/usr/bin/env python3
"""Analyze a reconciled sampling calibration and freeze its sampling policy."""

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
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
ANALYZER_VERSION = "0.1.0"
INSTANCE_FIELDS = (
    "check",
    "files",
    "reads",
    "regions",
    "ontology_terms",
    "metric_id",
    "metric_name",
    "unit",
)
CASE_ACCURACY_FIELDS = (
    "analysis_run_id",
    "study_run_id",
    "endpoint_id",
    "check",
    "metric_name",
    "unit",
    "case_id",
    "family_id",
    "configuration_accession",
    "fastq_accession",
    "selection_role",
    "requested_records_per_fastq",
    "sampling_method",
    "metric_instances",
    "observations",
    "case_mean_signed_error",
    "case_median_absolute_error",
    "case_max_absolute_error",
)
ACCURACY_SUMMARY_FIELDS = (
    "analysis_run_id",
    "study_run_id",
    "endpoint_id",
    "check",
    "metric_name",
    "unit",
    "requested_records_per_fastq",
    "sampling_method",
    "independent_configurations",
    "fastqs",
    "families",
    "metric_instances",
    "median_configuration_max_absolute_error",
    "p95_configuration_max_absolute_error",
    "mean_configuration_mean_signed_error",
    "minimum_configurations_required",
    "sufficient_evidence",
    "accuracy_pass",
)
CASE_BIAS_FIELDS = (
    "analysis_run_id",
    "study_run_id",
    "endpoint_id",
    "check",
    "metric_name",
    "unit",
    "case_id",
    "family_id",
    "configuration_accession",
    "fastq_accession",
    "selection_role",
    "requested_records_per_fastq",
    "metric_instances",
    "reservoir_seeds",
    "case_mean_prefix_minus_reservoir",
    "case_max_absolute_prefix_minus_reservoir",
)
BIAS_SUMMARY_FIELDS = (
    "analysis_run_id",
    "study_run_id",
    "endpoint_id",
    "check",
    "metric_name",
    "unit",
    "requested_records_per_fastq",
    "independent_configurations",
    "fastqs",
    "families",
    "metric_instances",
    "mean_prefix_minus_reservoir",
    "median_prefix_minus_reservoir",
    "bootstrap_ci_lower",
    "bootstrap_ci_upper",
    "confidence_level",
    "bootstrap_resamples",
    "minimum_configurations_required",
    "sufficient_evidence",
    "bias_pass",
)
PERFORMANCE_FIELDS = (
    "analysis_run_id",
    "study_run_id",
    "requested_records_per_fastq",
    "sampling_method",
    "measured_runs",
    "independent_fastqs",
    "median_wall_time_seconds",
    "p95_wall_time_seconds",
    "median_user_cpu_seconds",
    "median_system_cpu_seconds",
    "max_peak_resident_memory_bytes",
    "median_records_per_second",
    "median_report_size_bytes",
)
DESCRIPTIVE_FIELDS = (
    "analysis_run_id",
    "study_run_id",
    "descriptive_metric_id",
    "case_id",
    "family_id",
    "configuration_accession",
    "fastq_accession",
    "selection_role",
    "condition_id",
    "sampling_method",
    "requested_records_per_fastq",
    "sampling_seed",
    "check",
    "files",
    "reads",
    "regions",
    "ontology_terms",
    "metric_id",
    "metric_name",
    "unit",
    "value",
)
DESCRIPTIVE_PAIR_FIELDS = (
    "analysis_run_id",
    "study_run_id",
    "descriptive_metric_id",
    "case_id",
    "family_id",
    "configuration_accession",
    "fastq_accession",
    "selection_role",
    "requested_records_per_fastq",
    "check",
    "files",
    "reads",
    "regions",
    "ontology_terms",
    "metric_id",
    "metric_name",
    "unit",
    "prefix_value",
    "mean_reservoir_value",
    "prefix_minus_reservoir",
    "reservoir_seeds",
)
EVENT_FIELDS = (
    "analysis_run_id",
    "records_per_fastq",
    "event_prevalence",
    "probability_at_least_one_event",
)
DECISION_FIELDS = (
    "analysis_run_id",
    "study_run_id",
    "endpoint_id",
    "check",
    "metric_name",
    "unit",
    "evidence_complete",
    "preferred_prefix_pass",
    "preferred_reservoir_pass",
    "escalation_prefix_pass",
    "escalation_reservoir_pass",
    "decision",
    "recommended_records_per_fastq",
    "recommended_sampling_method",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze sampling accuracy, bias, and measured performance."
    )
    parser.add_argument("--study-manifest", required=True, type=Path)
    parser.add_argument("--analysis-protocol", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = analyze_sampling_calibration(
            study_manifest_path=args.study_manifest.resolve(),
            analysis_protocol_path=args.analysis_protocol.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"analyze_sampling_calibration: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def analyze_sampling_calibration(
    *,
    study_manifest_path: Path,
    analysis_protocol_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    inputs = load_and_validate_inputs(study_manifest_path, analysis_protocol_path)
    study = inputs["study"]
    sampling_protocol = inputs["sampling_protocol"]
    analysis_protocol = inputs["analysis_protocol"]
    analyzer = runtime.script_identity(
        Path(__file__).resolve(), version=ANALYZER_VERSION
    )
    runtime_tool = runtime.script_identity(
        Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
    )
    analysis_identity = {
        "schema_version": SCHEMA_VERSION,
        "study_manifest_sha256": runtime.file_sha256(study_manifest_path),
        "study_run_id": study["study_run_id"],
        "analysis_protocol_sha256": runtime.file_sha256(analysis_protocol_path),
        "analyzer": runtime.functional_script_identity(analyzer),
        "runtime": runtime.functional_script_identity(runtime_tool),
    }
    analysis_run_id = runtime.sha256_json(analysis_identity)[:16]

    selected_errors = select_endpoint_errors(
        inputs["error_rows"], analysis_protocol, study["study_run_id"]
    )
    case_accuracy = build_case_accuracy(selected_errors, analysis_run_id)
    accuracy_summary = summarize_accuracy(
        case_accuracy, analysis_protocol, analysis_run_id
    )
    case_bias, bias_issues = build_case_bias(
        selected_errors,
        sampling_protocol,
        analysis_run_id,
    )
    bias_summary = summarize_bias(case_bias, analysis_protocol, analysis_run_id)
    performance_summary, memory = summarize_performance(
        inputs["performance_rows"],
        analysis_protocol,
        study["study_run_id"],
        analysis_run_id,
    )
    descriptive_rows = select_descriptive_metrics(
        inputs["metric_rows"],
        analysis_protocol,
        study["study_run_id"],
        analysis_run_id,
    )
    descriptive_pairs, descriptive_issues = build_descriptive_pairs(
        descriptive_rows,
        sampling_protocol,
    )
    event_rows = build_event_probabilities(
        sampling_protocol, analysis_protocol, analysis_run_id
    )
    decisions = build_endpoint_decisions(
        accuracy_summary,
        bias_summary,
        analysis_protocol,
        study["study_run_id"],
        analysis_run_id,
    )
    policy = build_policy(
        study_run_id=study["study_run_id"],
        analysis_run_id=analysis_run_id,
        analysis_protocol_path=analysis_protocol_path,
        analysis_protocol=analysis_protocol,
        decisions=decisions,
        memory=memory,
    )

    output_root.mkdir(parents=True)
    tables = {
        "case_accuracy": (
            output_root / "tables" / "case_accuracy.csv",
            case_accuracy,
            CASE_ACCURACY_FIELDS,
        ),
        "accuracy_summary": (
            output_root / "tables" / "accuracy_summary.csv",
            accuracy_summary,
            ACCURACY_SUMMARY_FIELDS,
        ),
        "case_prefix_bias": (
            output_root / "tables" / "case_prefix_bias.csv",
            case_bias,
            CASE_BIAS_FIELDS,
        ),
        "prefix_bias_summary": (
            output_root / "tables" / "prefix_bias_summary.csv",
            bias_summary,
            BIAS_SUMMARY_FIELDS,
        ),
        "performance_summary": (
            output_root / "tables" / "performance_summary.csv",
            performance_summary,
            PERFORMANCE_FIELDS,
        ),
        "descriptive_metrics": (
            output_root / "tables" / "descriptive_metrics.csv",
            descriptive_rows,
            DESCRIPTIVE_FIELDS,
        ),
        "descriptive_prefix_reservoir": (
            output_root / "tables" / "descriptive_prefix_reservoir.csv",
            descriptive_pairs,
            DESCRIPTIVE_PAIR_FIELDS,
        ),
        "event_detection_probability": (
            output_root / "tables" / "event_detection_probability.csv",
            event_rows,
            EVENT_FIELDS,
        ),
        "endpoint_decisions": (
            output_root / "tables" / "endpoint_decisions.csv",
            decisions,
            DECISION_FIELDS,
        ),
    }
    for path, rows, fields in tables.values():
        runtime.write_csv(path, rows, list(fields))
    policy_path = output_root / "policy" / "sampling_policy.json"
    runtime.write_json(policy_path, policy)

    validation = validate_outputs(
        analysis_run_id=analysis_run_id,
        study_run_id=study["study_run_id"],
        analysis_protocol=analysis_protocol,
        sampling_protocol=sampling_protocol,
        selected_errors=selected_errors,
        case_accuracy=case_accuracy,
        accuracy_summary=accuracy_summary,
        case_bias=case_bias,
        bias_summary=bias_summary,
        bias_issues=bias_issues,
        performance_summary=performance_summary,
        descriptive_rows=descriptive_rows,
        descriptive_pairs=descriptive_pairs,
        descriptive_issues=descriptive_issues,
        event_rows=event_rows,
        decisions=decisions,
        policy=policy,
    )
    validation_path = output_root / "validation" / "sampling_analysis.json"
    runtime.write_json(validation_path, validation)
    if not validation["valid"]:
        raise ValueError(f"sampling analysis validation failed: {validation_path}")

    manifest_path = output_root / "manifests" / "sampling_analysis.json"
    manifest = {
        **analysis_identity,
        "analysis_run_id": analysis_run_id,
        "created_at": utc_now(),
        "valid": True,
        "policy_frozen": policy["frozen"],
        "inputs": {
            "study_manifest": runtime.file_identity(study_manifest_path),
            "analysis_protocol": runtime.file_identity(analysis_protocol_path),
            "sampling_protocol": runtime.file_identity(
                inputs["sampling_protocol_path"]
            ),
        },
        "tools": {"analyzer": analyzer, "runtime": runtime_tool},
        "counts": validation["counts"],
        "outputs": {
            key: runtime.file_identity(path) for key, (path, _, _) in tables.items()
        }
        | {
            "policy": runtime.file_identity(policy_path),
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_and_validate_inputs(
    study_manifest_path: Path, analysis_protocol_path: Path
) -> dict[str, Any]:
    study = runtime.load_json(study_manifest_path)
    analysis_protocol = runtime.load_json(analysis_protocol_path)
    if study.get("valid") is not True:
        raise ValueError("sampling calibration study is not valid")
    if study.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling calibration study does not use schema 0.1.0")
    study_run_id = required_string(study, "study_run_id")
    validation_path = verified_output(study, "validation")
    source_validation = runtime.load_json(validation_path)
    if source_validation.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling calibration validation schema is unsupported")
    if source_validation.get("valid") is not True:
        raise ValueError("sampling calibration validation is not valid")
    if source_validation.get("study_run_id") != study_run_id:
        raise ValueError("sampling calibration validation study id differs")

    sampling_identity = study.get("inputs", {}).get("sampling_protocol", {})
    if not isinstance(sampling_identity, dict):
        raise ValueError("sampling calibration has no protocol identity")
    sampling_protocol_path = Path(required_string(sampling_identity, "path")).resolve()
    if runtime.file_sha256(sampling_protocol_path) != required_string(
        sampling_identity, "sha256"
    ):
        raise ValueError("sampling calibration protocol hash changed")
    if study.get("sampling_protocol_sha256") != sampling_identity["sha256"]:
        raise ValueError("sampling calibration protocol identities differ")
    sampling_protocol = runtime.load_json(sampling_protocol_path)
    validate_analysis_protocol(analysis_protocol, sampling_protocol)

    return {
        "study": study,
        "analysis_protocol": analysis_protocol,
        "sampling_protocol": sampling_protocol,
        "sampling_protocol_path": sampling_protocol_path,
        "error_rows": runtime.read_csv(verified_output(study, "sampling_errors")),
        "performance_rows": runtime.read_csv(verified_output(study, "performance")),
        "metric_rows": runtime.read_csv(verified_output(study, "metrics")),
    }


def verified_output(study: dict[str, Any], key: str) -> Path:
    identity = study.get("outputs", {}).get(key, {})
    if not isinstance(identity, dict):
        raise ValueError(f"sampling calibration has no {key} output identity")
    path = Path(required_string(identity, "path")).resolve()
    if runtime.file_sha256(path) != required_string(identity, "sha256"):
        raise ValueError(f"sampling calibration {key} hash changed")
    return path


def validate_analysis_protocol(
    protocol: dict[str, Any], sampling_protocol: dict[str, Any]
) -> None:
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling analysis protocol does not use schema 0.1.0")
    required_string(protocol, "analysis_id")
    if sampling_protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling protocol does not use schema 0.1.0")
    sampling = sampling_protocol.get("sampling", {})
    sizes = sampling.get("records_per_fastq")
    seeds = sampling.get("reservoir_seeds")
    if (
        not isinstance(sizes, list)
        or not sizes
        or any(type(value) is not int or value <= 0 for value in sizes)
        or len(sizes) != len(set(sizes))
    ):
        raise ValueError("sampling protocol record counts are invalid")
    if (
        not isinstance(seeds, list)
        or not seeds
        or any(type(value) is not int or value < 0 for value in seeds)
        or len(seeds) != len(set(seeds))
    ):
        raise ValueError("sampling protocol reservoir seeds are invalid")
    if protocol.get("accuracy_aggregation") != (
        "configuration_max_of_fastq_case_max_absolute_error"
    ):
        raise ValueError("unsupported sampling accuracy aggregation")
    if protocol.get("prefix_bias_aggregation") != (
        "configuration_mean_of_fastq_case_means"
    ):
        raise ValueError("unsupported prefix bias aggregation")
    policy = protocol.get("policy", {})
    statistics_value = protocol.get("statistics", {})
    preferred = positive_int(policy.get("preferred_records_per_fastq"), "preferred")
    escalation = positive_int(policy.get("escalation_records_per_fastq"), "escalation")
    if preferred not in sizes or escalation not in sizes or preferred >= escalation:
        raise ValueError("analysis policy sizes do not match the sampling protocol")
    reservoir_seed = nonnegative_int(policy.get("reservoir_seed"), "reservoir seed")
    if reservoir_seed not in seeds:
        raise ValueError(
            "analysis policy reservoir seed is not in the sampling protocol"
        )
    thresholds = {}
    for field in (
        "median_absolute_error_max",
        "p95_absolute_error_max",
        "prefix_bias_ci_lower_min",
        "prefix_bias_ci_upper_max",
    ):
        thresholds[field] = finite_float(policy.get(field), field)
    if (
        thresholds["median_absolute_error_max"] < 0
        or thresholds["p95_absolute_error_max"] < 0
        or thresholds["median_absolute_error_max"]
        > thresholds["p95_absolute_error_max"]
    ):
        raise ValueError("sampling accuracy thresholds are inconsistent")
    if not (
        thresholds["prefix_bias_ci_lower_min"]
        < 0
        < thresholds["prefix_bias_ci_upper_max"]
    ):
        raise ValueError("prefix bias interval must span zero")
    positive_int(policy.get("peak_resident_memory_bytes_max"), "peak memory threshold")
    minimum = positive_int(
        statistics_value.get("minimum_independent_configurations_per_endpoint"),
        "minimum independent configurations",
    )
    if minimum < 2:
        raise ValueError("minimum independent configurations must be at least two")
    confidence = finite_float(
        statistics_value.get("confidence_level"), "confidence level"
    )
    if not 0 < confidence < 1:
        raise ValueError("confidence level must be between zero and one")
    positive_int(statistics_value.get("bootstrap_resamples"), "bootstrap resamples")
    nonnegative_int(statistics_value.get("bootstrap_seed"), "bootstrap seed")

    endpoints = protocol.get("endpoint_metrics")
    if not isinstance(endpoints, list) or not endpoints:
        raise ValueError("sampling analysis needs endpoint metrics")
    endpoint_ids = []
    selectors = []
    for endpoint in endpoints:
        if not isinstance(endpoint, dict):
            raise ValueError("endpoint metrics must be objects")
        endpoint_ids.append(required_string(endpoint, "endpoint_id"))
        selectors.append(
            (
                required_string(endpoint, "check"),
                required_string(endpoint, "metric_name"),
                required_string(endpoint, "unit"),
            )
        )
        if selectors[-1][2] != "fraction":
            raise ValueError("sampling policy endpoints must use fraction units")
    if len(endpoint_ids) != len(set(endpoint_ids)):
        raise ValueError("endpoint identifiers are not unique")
    if len(selectors) != len(set(selectors)):
        raise ValueError("endpoint metric selectors are not unique")

    descriptors = protocol.get("descriptive_metrics", [])
    if not isinstance(descriptors, list):
        raise ValueError("descriptive metrics must be a list")
    descriptor_ids = []
    for descriptor in descriptors:
        if not isinstance(descriptor, dict):
            raise ValueError("descriptive metrics must be objects")
        descriptor_ids.append(required_string(descriptor, "metric_id"))
        for field in ("check", "metric_name", "unit", "comparison"):
            required_string(descriptor, field)
        if descriptor["comparison"] != "prefix_vs_reservoir_at_equal_sample_size":
            raise ValueError("unsupported descriptive metric comparison")
    if len(descriptor_ids) != len(set(descriptor_ids)):
        raise ValueError("descriptive metric identifiers are not unique")

    prevalences = protocol.get("event_prevalences")
    if not isinstance(prevalences, list) or not prevalences:
        raise ValueError("event prevalence grid is empty")
    if any(
        not 0 < finite_float(value, "event prevalence") < 1 for value in prevalences
    ):
        raise ValueError("event prevalences must be between zero and one")
    if len(prevalences) != len(set(prevalences)):
        raise ValueError("event prevalences are not unique")


def select_endpoint_errors(
    rows: list[dict[str, str]], protocol: dict[str, Any], study_run_id: str
) -> list[dict[str, Any]]:
    selectors = {
        (value["check"], value["metric_name"], value["unit"]): value
        for value in protocol["endpoint_metrics"]
    }
    selected = []
    case_metadata: dict[str, tuple[str, ...]] = {}
    for row in rows:
        if row.get("study_run_id") != study_run_id:
            raise ValueError("sampling error row study id differs")
        case_id = required_string(row, "case_id")
        metadata = tuple(
            required_string(row, field)
            for field in (
                "family_id",
                "configuration_accession",
                "fastq_accession",
                "selection_role",
            )
        )
        if case_id in case_metadata and case_metadata[case_id] != metadata:
            raise ValueError(f"{case_id}: sampling error metadata differs across rows")
        case_metadata[case_id] = metadata
        endpoint = selectors.get(
            (row.get("check", ""), row.get("metric_name", ""), row.get("unit", ""))
        )
        if endpoint is None:
            continue
        method = row.get("sampling_method", "")
        if method not in {"prefix", "reservoir"}:
            raise ValueError("sampling error row has unsupported method")
        signed = finite_float(row.get("signed_error"), "signed error")
        absolute = finite_float(row.get("absolute_error"), "absolute error")
        complete = finite_float(row.get("complete_value"), "complete value")
        sample = finite_float(row.get("sample_value"), "sample value")
        if not math.isclose(absolute, abs(signed), rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("sampling absolute error does not match signed error")
        if not math.isclose(sample - complete, signed, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("sampling signed error does not match metric values")
        seed = optional_int(row.get("sampling_seed", ""), "sampling seed")
        if method == "prefix" and seed is not None:
            raise ValueError("prefix sampling error unexpectedly has a seed")
        selected.append(
            {
                **row,
                "endpoint_id": endpoint["endpoint_id"],
                "requested_records_per_fastq": positive_int(
                    row.get("requested_records_per_fastq"), "sample size"
                ),
                "sampling_seed": seed,
                "signed_error": signed,
                "absolute_error": absolute,
                "complete_value": complete,
                "sample_value": sample,
                "instance_key": tuple(row.get(field, "") for field in INSTANCE_FIELDS),
            }
        )
    return selected


def build_case_accuracy(
    rows: list[dict[str, Any]], analysis_run_id: str
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = (
            row["case_id"],
            row["endpoint_id"],
            str(row["requested_records_per_fastq"]),
            row["sampling_method"],
        )
        grouped[key].append(row)
    result = []
    for values in grouped.values():
        first = values[0]
        signed = [value["signed_error"] for value in values]
        absolute = [value["absolute_error"] for value in values]
        result.append(
            {
                **analysis_context(first, analysis_run_id),
                "requested_records_per_fastq": first["requested_records_per_fastq"],
                "sampling_method": first["sampling_method"],
                "metric_instances": len({value["instance_key"] for value in values}),
                "observations": len(values),
                "case_mean_signed_error": statistics.fmean(signed),
                "case_median_absolute_error": statistics.median(absolute),
                "case_max_absolute_error": max(absolute),
            }
        )
    return sorted(
        result,
        key=lambda row: (
            row["endpoint_id"],
            row["requested_records_per_fastq"],
            row["sampling_method"],
            row["case_id"],
        ),
    )


def summarize_accuracy(
    rows: list[dict[str, Any]], protocol: dict[str, Any], analysis_run_id: str
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, int, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[
            (
                row["endpoint_id"],
                row["requested_records_per_fastq"],
                row["sampling_method"],
            )
        ].append(row)
    minimum = protocol["statistics"]["minimum_independent_configurations_per_endpoint"]
    policy = protocol["policy"]
    result = []
    for values in grouped.values():
        first = values[0]
        by_configuration: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for value in values:
            by_configuration[value["configuration_accession"]].append(value)
        configuration_maxima = [
            max(value["case_max_absolute_error"] for value in configuration_values)
            for configuration_values in by_configuration.values()
        ]
        configuration_means = [
            statistics.fmean(
                value["case_mean_signed_error"] for value in configuration_values
            )
            for configuration_values in by_configuration.values()
        ]
        sufficient = len(by_configuration) >= minimum
        median_error = statistics.median(configuration_maxima)
        p95_error = percentile(configuration_maxima, 0.95)
        result.append(
            {
                "analysis_run_id": analysis_run_id,
                "study_run_id": first["study_run_id"],
                "endpoint_id": first["endpoint_id"],
                "check": first["check"],
                "metric_name": first["metric_name"],
                "unit": first["unit"],
                "requested_records_per_fastq": first["requested_records_per_fastq"],
                "sampling_method": first["sampling_method"],
                "independent_configurations": len(by_configuration),
                "fastqs": len(values),
                "families": len({value["family_id"] for value in values}),
                "metric_instances": sum(value["metric_instances"] for value in values),
                "median_configuration_max_absolute_error": median_error,
                "p95_configuration_max_absolute_error": p95_error,
                "mean_configuration_mean_signed_error": statistics.fmean(
                    configuration_means
                ),
                "minimum_configurations_required": minimum,
                "sufficient_evidence": sufficient,
                "accuracy_pass": (
                    sufficient
                    and median_error <= policy["median_absolute_error_max"]
                    and p95_error <= policy["p95_absolute_error_max"]
                ),
            }
        )
    return sorted(
        result,
        key=lambda row: (
            row["endpoint_id"],
            row["requested_records_per_fastq"],
            row["sampling_method"],
        ),
    )


def build_case_bias(
    rows: list[dict[str, Any]],
    sampling_protocol: dict[str, Any],
    analysis_run_id: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    expected_seeds = set(sampling_protocol["sampling"]["reservoir_seeds"])
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = (
            row["case_id"],
            row["endpoint_id"],
            row["requested_records_per_fastq"],
            row["instance_key"],
        )
        grouped[key].append(row)
    paired_by_case: dict[tuple[str, str, int], list[tuple[dict[str, Any], float]]] = (
        defaultdict(list)
    )
    issues = []
    for key, values in grouped.items():
        prefix = [value for value in values if value["sampling_method"] == "prefix"]
        reservoir = [
            value for value in values if value["sampling_method"] == "reservoir"
        ]
        observed_seeds = [value["sampling_seed"] for value in reservoir]
        if (
            len(prefix) != 1
            or set(observed_seeds) != expected_seeds
            or len(observed_seeds) != len(expected_seeds)
        ):
            issues.append(
                {
                    "type": "incomplete_prefix_reservoir_pair",
                    "case_id": key[0],
                    "endpoint_id": key[1],
                    "requested_records_per_fastq": key[2],
                    "prefix_rows": len(prefix),
                    "reservoir_seeds": sorted(
                        seed for seed in observed_seeds if seed is not None
                    ),
                }
            )
            continue
        complete_values = {value["complete_value"] for value in values}
        if len(complete_values) != 1:
            issues.append(
                {
                    "type": "complete_value_mismatch",
                    "case_id": key[0],
                    "endpoint_id": key[1],
                    "requested_records_per_fastq": key[2],
                }
            )
            continue
        difference = prefix[0]["sample_value"] - statistics.fmean(
            value["sample_value"] for value in reservoir
        )
        paired_by_case[(key[0], key[1], key[2])].append((prefix[0], difference))

    result = []
    for values in paired_by_case.values():
        first = values[0][0]
        differences = [difference for _, difference in values]
        result.append(
            {
                **analysis_context(first, analysis_run_id),
                "requested_records_per_fastq": first["requested_records_per_fastq"],
                "metric_instances": len(values),
                "reservoir_seeds": ";".join(
                    str(seed) for seed in sorted(expected_seeds)
                ),
                "case_mean_prefix_minus_reservoir": statistics.fmean(differences),
                "case_max_absolute_prefix_minus_reservoir": max(
                    abs(value) for value in differences
                ),
            }
        )
    return (
        sorted(
            result,
            key=lambda row: (
                row["endpoint_id"],
                row["requested_records_per_fastq"],
                row["case_id"],
            ),
        ),
        issues,
    )


def summarize_bias(
    rows: list[dict[str, Any]], protocol: dict[str, Any], analysis_run_id: str
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(row["endpoint_id"], row["requested_records_per_fastq"])].append(row)
    statistics_protocol = protocol["statistics"]
    policy = protocol["policy"]
    minimum = statistics_protocol["minimum_independent_configurations_per_endpoint"]
    result = []
    for values in grouped.values():
        first = values[0]
        by_configuration: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for value in values:
            by_configuration[value["configuration_accession"]].append(value)
        differences = sorted(
            statistics.fmean(
                value["case_mean_prefix_minus_reservoir"]
                for value in configuration_values
            )
            for configuration_values in by_configuration.values()
        )
        lower, upper = bootstrap_mean_interval(
            differences,
            confidence=statistics_protocol["confidence_level"],
            resamples=statistics_protocol["bootstrap_resamples"],
            seed=bootstrap_seed(
                statistics_protocol["bootstrap_seed"],
                first["endpoint_id"],
                first["requested_records_per_fastq"],
            ),
        )
        sufficient = len(by_configuration) >= minimum
        result.append(
            {
                "analysis_run_id": analysis_run_id,
                "study_run_id": first["study_run_id"],
                "endpoint_id": first["endpoint_id"],
                "check": first["check"],
                "metric_name": first["metric_name"],
                "unit": first["unit"],
                "requested_records_per_fastq": first["requested_records_per_fastq"],
                "independent_configurations": len(by_configuration),
                "fastqs": len(values),
                "families": len({value["family_id"] for value in values}),
                "metric_instances": sum(value["metric_instances"] for value in values),
                "mean_prefix_minus_reservoir": statistics.fmean(differences),
                "median_prefix_minus_reservoir": statistics.median(differences),
                "bootstrap_ci_lower": lower,
                "bootstrap_ci_upper": upper,
                "confidence_level": statistics_protocol["confidence_level"],
                "bootstrap_resamples": statistics_protocol["bootstrap_resamples"],
                "minimum_configurations_required": minimum,
                "sufficient_evidence": sufficient,
                "bias_pass": (
                    sufficient
                    and lower >= policy["prefix_bias_ci_lower_min"]
                    and upper <= policy["prefix_bias_ci_upper_max"]
                ),
            }
        )
    return sorted(
        result,
        key=lambda row: (
            row["endpoint_id"],
            row["requested_records_per_fastq"],
        ),
    )


def summarize_performance(
    rows: list[dict[str, str]],
    protocol: dict[str, Any],
    study_run_id: str,
    analysis_run_id: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    grouped: dict[tuple[int, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row.get("study_run_id") != study_run_id:
            raise ValueError("performance row study id differs")
        if parse_bool(row.get("is_warmup", "")):
            continue
        if row.get("execution_scope") != "local_compute":
            continue
        method = row.get("sampling_method", "")
        if method not in {"prefix", "reservoir"}:
            continue
        n_reads = positive_int(row.get("requested_records_per_fastq"), "sample size")
        wall = positive_float(row.get("wall_time_seconds"), "wall time")
        records = positive_int(row.get("records_processed"), "records processed")
        grouped[(n_reads, method)].append(
            {
                "case_id": required_string(row, "case_id"),
                "wall": wall,
                "user": nonnegative_float(row.get("user_cpu_seconds"), "user CPU"),
                "system": nonnegative_float(
                    row.get("system_cpu_seconds"), "system CPU"
                ),
                "rss": positive_int(
                    row.get("peak_resident_memory_bytes"), "peak memory"
                ),
                "throughput": records / wall,
                "report_size": positive_int(
                    row.get("report_size_bytes"), "report size"
                ),
            }
        )
    result = []
    for (n_reads, method), values in grouped.items():
        walls = [value["wall"] for value in values]
        result.append(
            {
                "analysis_run_id": analysis_run_id,
                "study_run_id": study_run_id,
                "requested_records_per_fastq": n_reads,
                "sampling_method": method,
                "measured_runs": len(values),
                "independent_fastqs": len({value["case_id"] for value in values}),
                "median_wall_time_seconds": statistics.median(walls),
                "p95_wall_time_seconds": percentile(walls, 0.95),
                "median_user_cpu_seconds": statistics.median(
                    value["user"] for value in values
                ),
                "median_system_cpu_seconds": statistics.median(
                    value["system"] for value in values
                ),
                "max_peak_resident_memory_bytes": max(value["rss"] for value in values),
                "median_records_per_second": statistics.median(
                    value["throughput"] for value in values
                ),
                "median_report_size_bytes": statistics.median(
                    value["report_size"] for value in values
                ),
            }
        )
    escalation = protocol["policy"]["escalation_records_per_fastq"]
    escalation_rows = [
        row for row in result if row["requested_records_per_fastq"] == escalation
    ]
    peak = (
        max(row["max_peak_resident_memory_bytes"] for row in escalation_rows)
        if escalation_rows
        else None
    )
    memory = {
        "records_per_fastq": escalation,
        "evaluable": peak is not None,
        "observed_peak_resident_memory_bytes": peak,
        "threshold_bytes": protocol["policy"]["peak_resident_memory_bytes_max"],
        "pass": (
            peak is not None
            and peak < protocol["policy"]["peak_resident_memory_bytes_max"]
        ),
    }
    return (
        sorted(
            result,
            key=lambda row: (
                row["requested_records_per_fastq"],
                row["sampling_method"],
            ),
        ),
        memory,
    )


def select_descriptive_metrics(
    rows: list[dict[str, str]],
    protocol: dict[str, Any],
    study_run_id: str,
    analysis_run_id: str,
) -> list[dict[str, Any]]:
    selectors = {
        (value["check"], value["metric_name"], value["unit"]): value
        for value in protocol["descriptive_metrics"]
    }
    result = []
    for row in rows:
        if row.get("study_run_id") != study_run_id:
            raise ValueError("metric row study id differs")
        if row.get("metric_side") != "observed" or row.get("data_kind") != "scalar":
            continue
        descriptor = selectors.get(
            (row.get("check", ""), row.get("metric_name", ""), row.get("unit", ""))
        )
        if descriptor is None:
            continue
        value = json.loads(row.get("value_json", "null"))
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError("descriptive metric value is not numeric")
        result.append(
            {
                "analysis_run_id": analysis_run_id,
                "study_run_id": study_run_id,
                "descriptive_metric_id": descriptor["metric_id"],
                **{
                    field: row.get(field, "")
                    for field in (
                        "case_id",
                        "family_id",
                        "configuration_accession",
                        "fastq_accession",
                        "selection_role",
                        "condition_id",
                        "sampling_method",
                        "sampling_seed",
                        "check",
                        "files",
                        "reads",
                        "regions",
                        "ontology_terms",
                        "metric_id",
                        "metric_name",
                        "unit",
                    )
                },
                "requested_records_per_fastq": positive_int(
                    row.get("requested_records_per_fastq"), "sample size"
                ),
                "value": finite_float(value, "descriptive metric value"),
            }
        )
    return sorted(
        result,
        key=lambda row: (
            row["descriptive_metric_id"],
            row["case_id"],
            row["requested_records_per_fastq"],
            row["sampling_method"],
            row["sampling_seed"],
            row["regions"],
        ),
    )


def build_descriptive_pairs(
    rows: list[dict[str, Any]], sampling_protocol: dict[str, Any]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    expected_seeds = set(sampling_protocol["sampling"]["reservoir_seeds"])
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["sampling_method"] not in {"prefix", "reservoir"}:
            continue
        key = (
            row["case_id"],
            row["descriptive_metric_id"],
            row["requested_records_per_fastq"],
            *(row[field] for field in INSTANCE_FIELDS),
        )
        grouped[key].append(row)
    result = []
    issues = []
    for values in grouped.values():
        prefix = [value for value in values if value["sampling_method"] == "prefix"]
        reservoir = [
            value for value in values if value["sampling_method"] == "reservoir"
        ]
        observed_seeds = [
            optional_int(value["sampling_seed"], "sampling seed") for value in reservoir
        ]
        if (
            len(prefix) != 1
            or set(observed_seeds) != expected_seeds
            or len(observed_seeds) != len(expected_seeds)
        ):
            issues.append(
                {
                    "type": "incomplete_descriptive_pair",
                    "case_id": values[0]["case_id"],
                    "descriptive_metric_id": values[0]["descriptive_metric_id"],
                    "requested_records_per_fastq": values[0][
                        "requested_records_per_fastq"
                    ],
                }
            )
            continue
        first = prefix[0]
        reservoir_mean = statistics.fmean(value["value"] for value in reservoir)
        result.append(
            {
                **{
                    field: first[field]
                    for field in DESCRIPTIVE_PAIR_FIELDS
                    if field
                    not in {
                        "prefix_value",
                        "mean_reservoir_value",
                        "prefix_minus_reservoir",
                        "reservoir_seeds",
                    }
                },
                "prefix_value": first["value"],
                "mean_reservoir_value": reservoir_mean,
                "prefix_minus_reservoir": first["value"] - reservoir_mean,
                "reservoir_seeds": ";".join(
                    str(seed) for seed in sorted(expected_seeds)
                ),
            }
        )
    return result, issues


def build_event_probabilities(
    sampling_protocol: dict[str, Any],
    analysis_protocol: dict[str, Any],
    analysis_run_id: str,
) -> list[dict[str, Any]]:
    return [
        {
            "analysis_run_id": analysis_run_id,
            "records_per_fastq": size,
            "event_prevalence": prevalence,
            "probability_at_least_one_event": -math.expm1(
                size * math.log1p(-prevalence)
            ),
        }
        for size in sorted(sampling_protocol["sampling"]["records_per_fastq"])
        for prevalence in sorted(analysis_protocol["event_prevalences"])
    ]


def build_endpoint_decisions(
    accuracy_rows: list[dict[str, Any]],
    bias_rows: list[dict[str, Any]],
    protocol: dict[str, Any],
    study_run_id: str,
    analysis_run_id: str,
) -> list[dict[str, Any]]:
    accuracy = {
        (
            row["endpoint_id"],
            row["requested_records_per_fastq"],
            row["sampling_method"],
        ): row
        for row in accuracy_rows
    }
    bias = {
        (row["endpoint_id"], row["requested_records_per_fastq"]): row
        for row in bias_rows
    }
    preferred = protocol["policy"]["preferred_records_per_fastq"]
    escalation = protocol["policy"]["escalation_records_per_fastq"]
    result = []
    for endpoint in protocol["endpoint_metrics"]:
        endpoint_id = endpoint["endpoint_id"]
        size_state = {}
        evidence = []
        for size in (preferred, escalation):
            prefix = accuracy.get((endpoint_id, size, "prefix"))
            reservoir = accuracy.get((endpoint_id, size, "reservoir"))
            bias_value = bias.get((endpoint_id, size))
            sufficient = all(
                value is not None and value["sufficient_evidence"]
                for value in (prefix, reservoir, bias_value)
            )
            evidence.append(sufficient)
            size_state[size] = {
                "prefix": bool(
                    sufficient and prefix["accuracy_pass"] and bias_value["bias_pass"]
                ),
                "reservoir": bool(sufficient and reservoir["accuracy_pass"]),
            }
        evidence_complete = all(evidence)
        recommendation = ("", None, "")
        if evidence_complete:
            for size, label in ((preferred, "default"), (escalation, "escalation")):
                for method in ("prefix", "reservoir"):
                    if size_state[size][method]:
                        recommendation = (f"{label}_{method}", size, method)
                        break
                if recommendation[0]:
                    break
            if not recommendation[0]:
                recommendation = ("quantitative_only", None, "")
        else:
            recommendation = ("insufficient_evidence", None, "")
        result.append(
            {
                "analysis_run_id": analysis_run_id,
                "study_run_id": study_run_id,
                **endpoint,
                "evidence_complete": evidence_complete,
                "preferred_prefix_pass": size_state[preferred]["prefix"],
                "preferred_reservoir_pass": size_state[preferred]["reservoir"],
                "escalation_prefix_pass": size_state[escalation]["prefix"],
                "escalation_reservoir_pass": size_state[escalation]["reservoir"],
                "decision": recommendation[0],
                "recommended_records_per_fastq": recommendation[1],
                "recommended_sampling_method": recommendation[2],
            }
        )
    return result


def build_policy(
    *,
    study_run_id: str,
    analysis_run_id: str,
    analysis_protocol_path: Path,
    analysis_protocol: dict[str, Any],
    decisions: list[dict[str, Any]],
    memory: dict[str, Any],
) -> dict[str, Any]:
    frozen = all(row["evidence_complete"] for row in decisions) and memory["evaluable"]
    stable_decisions = [
        {
            key: value
            for key, value in row.items()
            if key not in {"analysis_run_id", "study_run_id"}
        }
        for row in decisions
    ]
    default_reservoir = any(row["decision"] == "default_reservoir" for row in decisions)
    stable_policy = {
        "schema_version": SCHEMA_VERSION,
        "study_run_id": study_run_id,
        "analysis_run_id": analysis_run_id,
        "analysis_protocol_sha256": runtime.file_sha256(analysis_protocol_path),
        "frozen": frozen,
        "scientific_targets_met": (
            frozen
            and memory["pass"]
            and all(
                row["decision"] not in {"quantitative_only", "insufficient_evidence"}
                for row in decisions
            )
        ),
        "default": {
            "records_per_fastq": analysis_protocol["policy"][
                "preferred_records_per_fastq"
            ],
            "sampling_method": "reservoir" if default_reservoir else "prefix",
            "sampling_seed": (
                analysis_protocol["policy"]["reservoir_seed"]
                if default_reservoir
                else None
            ),
        },
        "escalation_records_per_fastq": analysis_protocol["policy"][
            "escalation_records_per_fastq"
        ],
        "endpoint_decisions": stable_decisions,
        "memory": memory,
    }
    return {
        **stable_policy,
        "policy_id": runtime.sha256_json(stable_policy)[:16],
        "created_at": utc_now(),
    }


def validate_outputs(
    *,
    analysis_run_id: str,
    study_run_id: str,
    analysis_protocol: dict[str, Any],
    sampling_protocol: dict[str, Any],
    selected_errors: list[dict[str, Any]],
    case_accuracy: list[dict[str, Any]],
    accuracy_summary: list[dict[str, Any]],
    case_bias: list[dict[str, Any]],
    bias_summary: list[dict[str, Any]],
    bias_issues: list[dict[str, Any]],
    performance_summary: list[dict[str, Any]],
    descriptive_rows: list[dict[str, Any]],
    descriptive_pairs: list[dict[str, Any]],
    descriptive_issues: list[dict[str, Any]],
    event_rows: list[dict[str, Any]],
    decisions: list[dict[str, Any]],
    policy: dict[str, Any],
) -> dict[str, Any]:
    stable_policy = {
        key: value
        for key, value in policy.items()
        if key not in {"policy_id", "created_at"}
    }
    all_rows = (
        case_accuracy
        + accuracy_summary
        + case_bias
        + bias_summary
        + performance_summary
        + descriptive_rows
        + descriptive_pairs
        + event_rows
        + decisions
    )
    expected_events = len(sampling_protocol["sampling"]["records_per_fastq"]) * len(
        analysis_protocol["event_prevalences"]
    )
    checks = {
        "endpoint_errors_present": bool(selected_errors),
        "case_accuracy_present": bool(case_accuracy),
        "prefix_reservoir_pairs_reconcile": not bias_issues,
        "performance_summary_present": bool(performance_summary),
        "descriptive_pairs_reconcile": not descriptive_issues,
        "event_grid_reconciles": len(event_rows) == expected_events,
        "endpoint_decisions_reconcile": len(decisions)
        == len(analysis_protocol["endpoint_metrics"]),
        "policy_id_reconciles": policy["policy_id"]
        == runtime.sha256_json(stable_policy)[:16],
        "analysis_ids_reconcile": all(
            row["analysis_run_id"] == analysis_run_id for row in all_rows
        ),
        "study_ids_reconcile": all(
            row.get("study_run_id", study_run_id) == study_run_id for row in all_rows
        ),
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "analysis_run_id": analysis_run_id,
        "study_run_id": study_run_id,
        "generated_at": utc_now(),
        "valid": all(checks.values()),
        "policy_frozen": policy["frozen"],
        "checks": checks,
        "counts": {
            "selected_error_rows": len(selected_errors),
            "case_accuracy_rows": len(case_accuracy),
            "accuracy_summary_rows": len(accuracy_summary),
            "case_bias_rows": len(case_bias),
            "bias_summary_rows": len(bias_summary),
            "performance_summary_rows": len(performance_summary),
            "descriptive_metric_rows": len(descriptive_rows),
            "descriptive_pair_rows": len(descriptive_pairs),
            "event_probability_rows": len(event_rows),
            "endpoint_decisions": len(decisions),
        },
        "details": {
            "prefix_reservoir_issues": bias_issues,
            "descriptive_pair_issues": descriptive_issues,
        },
    }


def analysis_context(row: dict[str, Any], analysis_run_id: str) -> dict[str, Any]:
    return {
        "analysis_run_id": analysis_run_id,
        **{
            field: row[field]
            for field in (
                "study_run_id",
                "endpoint_id",
                "check",
                "metric_name",
                "unit",
                "case_id",
                "family_id",
                "configuration_accession",
                "fastq_accession",
                "selection_role",
            )
        },
    }


def bootstrap_mean_interval(
    values: list[float], *, confidence: float, resamples: int, seed: int
) -> tuple[float, float]:
    if not values:
        raise ValueError("cannot bootstrap an empty value list")
    generator = random.Random(seed)
    count = len(values)
    estimates = sorted(
        statistics.fmean(values[generator.randrange(count)] for _ in range(count))
        for _ in range(resamples)
    )
    alpha = (1 - confidence) / 2
    return percentile(estimates, alpha), percentile(estimates, 1 - alpha)


def bootstrap_seed(base: int, endpoint_id: str, n_reads: int) -> int:
    offset = int(runtime.sha256_json([endpoint_id, n_reads])[:16], 16)
    return base + offset


def percentile(values: list[float], probability: float) -> float:
    if not values:
        raise ValueError("cannot compute a percentile of an empty value list")
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def required_string(value: dict[str, Any], field: str) -> str:
    result = str(value.get(field, "")).strip()
    if not result:
        raise ValueError(f"missing required field: {field}")
    return result


def finite_float(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be numeric")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric") from error
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def positive_float(value: Any, label: str) -> float:
    result = finite_float(value, label)
    if result <= 0:
        raise ValueError(f"{label} must be positive")
    return result


def nonnegative_float(value: Any, label: str) -> float:
    result = finite_float(value, label)
    if result < 0:
        raise ValueError(f"{label} must be nonnegative")
    return result


def positive_int(value: Any, label: str) -> int:
    result = strict_int(value, label)
    if result <= 0:
        raise ValueError(f"{label} must be positive")
    return result


def nonnegative_int(value: Any, label: str) -> int:
    result = strict_int(value, label)
    if result < 0:
        raise ValueError(f"{label} must be nonnegative")
    return result


def optional_int(value: Any, label: str) -> int | None:
    if value is None or str(value).strip() == "":
        return None
    return strict_int(value, label)


def strict_int(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer")
    text = str(value).strip()
    try:
        result = int(text)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer") from error
    if text not in {str(result), f"+{result}"}:
        raise ValueError(f"{label} must be an integer")
    return result


def parse_bool(value: Any) -> bool:
    if value is True or value == "True" or value == "true" or value == "1":
        return True
    if value is False or value == "False" or value == "false" or value == "0":
        return False
    raise ValueError(f"invalid boolean value: {value}")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
