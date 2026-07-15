#!/usr/bin/env python3
"""Analyze the frozen current-IGVF seqcheck audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import igvf_audit as audit
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import igvf_audit as audit
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
TOOL_VERSION = "0.1.0"
COMPLETION_FIELDS = (
    "analysis_id",
    "dimension",
    "group_id",
    "group_label",
    "configuration_count",
    "completed_runs",
    "classified_failures",
    "eligibility_exclusions",
    "eligible_outcomes",
    "completion_fraction",
    "minimum_completion_fraction",
    "completion_target_met",
    "sampled_records",
)
FAILURE_FIELDS = (
    "analysis_id",
    "failure_category",
    "access_class",
    "stage",
    "reason",
    "assay_family_ids",
    "failure_count",
    "configuration_count",
    "maximum_attempt_count",
)
METRIC_FIELDS = (
    "analysis_id",
    "dimension",
    "group_id",
    "group_label",
    "access_class",
    "check",
    "metric_name",
    "unit",
    "raw_observations",
    "independent_configurations",
    "minimum",
    "q1",
    "median",
    "mean",
    "q3",
    "maximum",
)
ASSESSMENT_FIELDS = (
    "analysis_id",
    "assay_family_id",
    "assay_family_label",
    "access_class",
    "check",
    "assessment_type",
    "assessment_code",
    "candidate_inconsistency",
    "assessment_count",
    "configuration_count",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze the frozen IGVF audit.")
    parser.add_argument("--audit-root", required=True, type=Path)
    parser.add_argument("--repeatability-validation", required=True, type=Path)
    parser.add_argument("--audit-protocol", required=True, type=Path)
    parser.add_argument("--family-rules", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        path = analyze_igvf_audit(
            audit_root=args.audit_root,
            repeatability_path=args.repeatability_validation,
            protocol_path=args.audit_protocol,
            family_rules_path=args.family_rules,
            output_root=args.output_root,
        )
    except (OSError, ValueError) as err:
        print(f"analyze_igvf_audit: {err}", file=sys.stderr)
        return 1
    print(path)
    return 0


def analyze_igvf_audit(
    *,
    audit_root: Path,
    repeatability_path: Path,
    protocol_path: Path,
    family_rules_path: Path,
    output_root: Path,
) -> Path:
    audit_root = audit_root.resolve()
    repeatability_path = repeatability_path.resolve()
    protocol_path = protocol_path.resolve()
    family_rules_path = family_rules_path.resolve()
    output_root = prepare_output_root(output_root)

    study = load_json(audit_root / "manifests" / "study.json", "study manifest")
    reconciliation = load_json(
        audit_root / "validation" / "reconciliation.json", "audit reconciliation"
    )
    repeatability = load_json(repeatability_path, "repeatability validation")
    protocol = audit.load_audit_protocol(protocol_path)
    family_rules = load_family_rules(family_rules_path)
    if reconciliation.get("valid") is not True:
        raise ValueError("audit reconciliation is not valid")
    if repeatability.get("valid") is not True:
        raise ValueError("repeatability validation is not valid")
    study_run_id = str(study.get("run_id", ""))
    if not study_run_id or reconciliation.get("study_run_id") != study_run_id:
        raise ValueError("study and reconciliation IDs differ")
    if repeatability.get("primary_study_run_id") != study_run_id:
        raise ValueError("repeatability validation belongs to another audit")
    protocol_id = audit.audit_protocol_id(protocol)
    study_protocol = study.get("frozen_inputs", {}).get("audit_protocol", {})
    if study_protocol.get("audit_protocol_id") != protocol_id:
        raise ValueError("audit study belongs to another audit protocol")
    repeat_protocol = repeatability.get("inputs", {}).get("protocol", {})
    if repeat_protocol.get("sha256") != runtime.file_sha256(protocol_path):
        raise ValueError("repeatability validation belongs to another audit protocol")
    reconciliation_counts = reconciliation.get("counts")
    if not isinstance(reconciliation_counts, dict):
        raise ValueError("audit reconciliation does not contain table counts")

    paths = {
        "runs": audit_root / "runs.csv",
        "failures": audit_root / "failures.csv",
        "diagnostics": audit_root / "diagnostics.csv",
        "metrics": audit_root / "metrics.csv",
    }
    runs = read_csv(paths["runs"])
    failures = read_csv(paths["failures"])
    diagnostics = read_csv(paths["diagnostics"])
    metrics = read_csv(paths["metrics"])
    validate_input_rows(study_run_id, runs, failures, diagnostics, metrics)
    if not runs and not failures:
        raise ValueError("audit contains no configuration outcomes")

    tool = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "study_run_id": study_run_id,
        "audit_protocol_id": protocol_id,
        "tool": runtime.functional_script_identity(tool),
        "inputs": {
            "study": runtime.file_sha256(audit_root / "manifests" / "study.json"),
            "reconciliation": runtime.file_sha256(
                audit_root / "validation" / "reconciliation.json"
            ),
            "repeatability": runtime.file_sha256(repeatability_path),
            "protocol": runtime.file_sha256(protocol_path),
            "family_rules": runtime.file_sha256(family_rules_path),
            **{name: runtime.file_sha256(path) for name, path in paths.items()},
        },
    }
    analysis_id = runtime.sha256_json(stable_identity)[:16]
    completion = build_completion_summary(
        runs=runs,
        failures=failures,
        families=family_rules,
        minimum_completion=protocol["targets"]["minimum_completion_fraction"],
        analysis_id=analysis_id,
    )
    failure_summary = build_failure_summary(failures, family_rules, analysis_id)
    metric_summary, numeric_metric_count = build_metric_distributions(
        metrics, family_rules, analysis_id
    )
    assessment_summary = build_assessment_summary(
        diagnostics, family_rules, analysis_id
    )

    tables = {
        "completion": (
            output_root / "tables" / "completion_summary.csv",
            completion,
            COMPLETION_FIELDS,
        ),
        "failures": (
            output_root / "tables" / "failure_summary.csv",
            failure_summary,
            FAILURE_FIELDS,
        ),
        "metrics": (
            output_root / "tables" / "metric_distributions.csv",
            metric_summary,
            METRIC_FIELDS,
        ),
        "assessments": (
            output_root / "tables" / "assessment_summary.csv",
            assessment_summary,
            ASSESSMENT_FIELDS,
        ),
    }
    for _, (path, rows, fields) in tables.items():
        write_csv(path, rows, fields)

    overall = next(
        row
        for row in completion
        if row["dimension"] == "overall" and row["group_id"] == "all"
    )
    unclassified = sum(
        1 for row in failures if row.get("failure_category") == "unclassified"
    )
    repeat_target_met = (
        parse_float(repeatability.get("match_fraction"), "repeat match fraction")
        >= protocol["targets"]["required_repeat_match_fraction"]
    )
    scientific_targets = {
        "completion_target_met": overall["completion_target_met"],
        "unclassified_failure_target_met": unclassified
        <= protocol["targets"]["maximum_unclassified_failure_count"],
        "repeatability_target_met": repeat_target_met,
    }
    run_keys = [(row["configuration_accession"], row["modality"]) for row in runs]
    run_key_set = set(run_keys)
    outcome_keys = [
        (row["configuration_accession"], row["modality"]) for row in [*runs, *failures]
    ]
    sampled_records = sum(parse_int(row["sampled_record_count"]) for row in runs)
    table_counts = {
        "selected_configurations": overall["configuration_count"],
        "runs": len(runs),
        "failures": len(failures),
        "diagnostics": len(diagnostics),
        "metrics": len(metrics),
        "sampled_records": sampled_records,
    }
    checks = {
        "audit_reconciliation_valid": reconciliation["valid"] is True,
        "repeatability_valid": repeatability["valid"] is True,
        "study_ids_reconcile": all(
            row.get("study_run_id") == study_run_id
            for row in [*runs, *failures, *diagnostics, *metrics]
        ),
        "configuration_outcomes_reconcile": overall["configuration_count"]
        == int(study.get("selected_configuration_count", 0)),
        "audit_table_counts_reconcile": all(
            parse_int(reconciliation_counts.get(name)) == count
            for name, count in table_counts.items()
        ),
        "run_keys_unique": len(run_keys) == len(run_key_set),
        "outcome_keys_unique": len(outcome_keys) == len(set(outcome_keys)),
        "flattened_rows_match_completed_runs": all(
            (row["configuration_accession"], row["modality"]) in run_key_set
            for row in [*diagnostics, *metrics]
        ),
        "access_classes_valid": all(
            row["access_class"] in audit.ACCESS_CLASSES
            for row in [*runs, *failures, *diagnostics, *metrics]
        ),
        "failure_categories_valid": all(
            row["failure_category"] in audit.FAILURE_CATEGORIES for row in failures
        ),
        "attempt_counts_valid": all(
            (row["run_status"] == "cached" and parse_int(row["max_attempt_count"]) >= 0)
            or 1
            <= parse_int(row["max_attempt_count"])
            <= protocol["execution"]["max_attempts"]
            for row in runs
        )
        and all(
            1
            <= parse_int(row["attempt_count"])
            <= protocol["execution"]["max_attempts"]
            for row in failures
        ),
        "numeric_metrics_reconcile": numeric_metric_count
        == sum(
            row["raw_observations"]
            for row in metric_summary
            if row["dimension"] == "overall"
        ),
        "metric_groups_have_configurations": all(
            row["independent_configurations"] > 0 for row in metric_summary
        ),
        "unclassified_failures_absent": scientific_targets[
            "unclassified_failure_target_met"
        ],
    }
    validation = {
        **stable_identity,
        "analysis_id": analysis_id,
        "generated_at": utc_now(),
        "valid": all(checks.values()),
        "scientific_target_met": all(scientific_targets.values()),
        "checks": checks,
        "scientific_targets": scientific_targets,
        "counts": {
            "runs": len(runs),
            "failures": len(failures),
            "diagnostics": len(diagnostics),
            "metrics": len(metrics),
            "numeric_observed_scalar_metrics": numeric_metric_count,
            "metric_distribution_rows": len(metric_summary),
            "assessment_summary_rows": len(assessment_summary),
            "sampled_records": sampled_records,
        },
        "overall_completion_fraction": overall["completion_fraction"],
        "outputs": {
            name: runtime.file_identity(path) for name, (path, _, _) in tables.items()
        },
        "tools": {"analyzer": tool},
    }
    validation_path = output_root / "validation" / "igvf_audit_analysis.json"
    write_json(validation_path, validation)
    if not validation["valid"]:
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise ValueError(f"IGVF audit analysis validation failed: {failed}")
    manifest = {
        **stable_identity,
        "analysis_id": analysis_id,
        "generated_at": utc_now(),
        "aggregation": {
            "metric_unit": "configuration median",
            "summary_unit": "independent configuration",
            "quantile_method": "linear interpolation, type 7",
            "assay_family_membership": "exact multi-membership with unclassified fallback",
            "warning_interpretation": "candidate inconsistency, not confirmed error",
        },
        "inputs": {
            "study": runtime.file_identity(audit_root / "manifests" / "study.json"),
            "reconciliation": runtime.file_identity(
                audit_root / "validation" / "reconciliation.json"
            ),
            "repeatability": runtime.file_identity(repeatability_path),
            "protocol": runtime.file_identity(protocol_path),
            "family_rules": runtime.file_identity(family_rules_path),
            **{name: runtime.file_identity(path) for name, path in paths.items()},
        },
        "outputs": {
            **validation["outputs"],
            "validation": runtime.file_identity(validation_path),
        },
        "tools": {"analyzer": tool},
    }
    manifest_path = output_root / "manifests" / "igvf_audit_analysis.json"
    write_json(manifest_path, manifest)
    return validation_path


def build_completion_summary(
    *,
    runs: list[dict[str, str]],
    failures: list[dict[str, str]],
    families: list[dict[str, Any]],
    minimum_completion: float,
    analysis_id: str,
) -> list[dict[str, Any]]:
    outcomes = []
    for row in runs:
        outcomes.append(
            {
                **row,
                "completed": True,
                "eligibility_exclusion": False,
                "sampled_records": parse_int(row["sampled_record_count"]),
                "family_ids": family_ids(row, families),
            }
        )
    for row in failures:
        outcomes.append(
            {
                **row,
                "completed": False,
                "eligibility_exclusion": row["failure_category"] == "eligibility",
                "sampled_records": 0,
                "family_ids": family_ids(row, families),
            }
        )
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    family_labels = {rule["id"]: rule["label"] for rule in families}
    family_labels["unclassified"] = "Unclassified assay family"
    for row in outcomes:
        memberships = [
            ("overall", "all", "All outcomes"),
            ("access_class", row["access_class"], row["access_class"]),
        ]
        for family_id in row["family_ids"]:
            memberships.append(("assay_family", family_id, family_labels[family_id]))
            memberships.append(
                (
                    "assay_family_access_class",
                    f"{family_id}|{row['access_class']}",
                    f"{family_labels[family_id]} | {row['access_class']}",
                )
            )
        for membership in memberships:
            grouped[membership].append(row)

    result = []
    for (dimension, group_id, group_label), rows in sorted(grouped.items()):
        completed = sum(row["completed"] for row in rows)
        exclusions = sum(row["eligibility_exclusion"] for row in rows)
        failures = len(rows) - completed - exclusions
        denominator = completed + failures
        completion = completed / denominator if denominator else None
        result.append(
            {
                "analysis_id": analysis_id,
                "dimension": dimension,
                "group_id": group_id,
                "group_label": group_label,
                "configuration_count": len(
                    {row["configuration_accession"] for row in rows}
                ),
                "completed_runs": completed,
                "classified_failures": failures,
                "eligibility_exclusions": exclusions,
                "eligible_outcomes": denominator,
                "completion_fraction": completion,
                "minimum_completion_fraction": minimum_completion,
                "completion_target_met": completion is not None
                and completion >= minimum_completion,
                "sampled_records": sum(row["sampled_records"] for row in rows),
            }
        )
    return result


def build_failure_summary(
    failures: list[dict[str, str]],
    families: list[dict[str, Any]],
    analysis_id: str,
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in failures:
        ids = family_ids(row, families)
        key = (
            row["failure_category"],
            row["access_class"],
            row["stage"],
            row["reason"],
            ";".join(ids),
        )
        grouped[key].append(row)
    return [
        {
            "analysis_id": analysis_id,
            "failure_category": key[0],
            "access_class": key[1],
            "stage": key[2],
            "reason": key[3],
            "assay_family_ids": key[4],
            "failure_count": len(rows),
            "configuration_count": len(
                {row["configuration_accession"] for row in rows}
            ),
            "maximum_attempt_count": max(
                parse_int(row["attempt_count"]) for row in rows
            ),
        }
        for key, rows in sorted(grouped.items())
    ]


def build_metric_distributions(
    metrics: list[dict[str, str]],
    families: list[dict[str, Any]],
    analysis_id: str,
) -> tuple[list[dict[str, Any]], int]:
    observations: dict[tuple[str, ...], list[tuple[str, float]]] = defaultdict(list)
    numeric_count = 0
    family_labels = {rule["id"]: rule["label"] for rule in families}
    family_labels["unclassified"] = "Unclassified assay family"
    for row in metrics:
        if row["metric_side"] != "observed" or row["data_kind"] != "scalar":
            continue
        try:
            value = json.loads(row["value_json"])
        except json.JSONDecodeError:
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        value = float(value)
        if not math.isfinite(value):
            continue
        numeric_count += 1
        groups = [("overall", "all", "All annotated metrics")]
        groups.extend(
            ("assay_family", family_id, family_labels[family_id])
            for family_id in family_ids(row, families)
        )
        groups.extend(
            ("sequence_type", sequence_type, sequence_type)
            for sequence_type in split_terms(row.get("sequence_types", ""))
            or ["unannotated"]
        )
        groups.extend(
            ("ontology_term", ontology_term, ontology_term)
            for ontology_term in split_terms(row.get("ontology_terms", ""))
            or ["unannotated"]
        )
        for dimension, group_id, group_label in groups:
            key = (
                dimension,
                group_id,
                group_label,
                row["access_class"],
                row["check"],
                row["metric_name"],
                row["unit"],
            )
            observations[key].append((row["configuration_accession"], value))

    result = []
    for key, raw_values in sorted(observations.items()):
        by_configuration: dict[str, list[float]] = defaultdict(list)
        for configuration, value in raw_values:
            by_configuration[configuration].append(value)
        values = sorted(
            statistics.median(configuration_values)
            for configuration_values in by_configuration.values()
        )
        result.append(
            {
                "analysis_id": analysis_id,
                "dimension": key[0],
                "group_id": key[1],
                "group_label": key[2],
                "access_class": key[3],
                "check": key[4],
                "metric_name": key[5],
                "unit": key[6],
                "raw_observations": len(raw_values),
                "independent_configurations": len(values),
                "minimum": values[0],
                "q1": linear_quantile(values, 0.25),
                "median": statistics.median(values),
                "mean": statistics.fmean(values),
                "q3": linear_quantile(values, 0.75),
                "maximum": values[-1],
            }
        )
    return result, numeric_count


def build_assessment_summary(
    diagnostics: list[dict[str, str]],
    families: list[dict[str, Any]],
    analysis_id: str,
) -> list[dict[str, Any]]:
    labels = {rule["id"]: rule["label"] for rule in families}
    labels["unclassified"] = "Unclassified assay family"
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in diagnostics:
        for family_id in family_ids(row, families):
            key = (
                family_id,
                labels[family_id],
                row["access_class"],
                row["check"],
                row["assessment_type"],
                row["assessment_code"],
            )
            grouped[key].append(row)
    return [
        {
            "analysis_id": analysis_id,
            "assay_family_id": key[0],
            "assay_family_label": key[1],
            "access_class": key[2],
            "check": key[3],
            "assessment_type": key[4],
            "assessment_code": key[5],
            "candidate_inconsistency": key[4] in {"warning", "error"},
            "assessment_count": len(rows),
            "configuration_count": len(
                {row["configuration_accession"] for row in rows}
            ),
        }
        for key, rows in sorted(grouped.items())
    ]


def load_family_rules(path: Path) -> list[dict[str, Any]]:
    value = load_json(path, "assay family rules")
    if value.get("schema_version") != "0.2.0":
        raise ValueError("unsupported assay-family rule schema")
    rules = value.get("families")
    if not isinstance(rules, list) or not rules:
        raise ValueError("assay-family registry contains no rules")
    result = []
    seen = set()
    for rule in rules:
        if not isinstance(rule, dict) or not rule.get("id") or not rule.get("label"):
            raise ValueError("assay-family rule is incomplete")
        rule_id = str(rule["id"])
        if rule_id in seen:
            raise ValueError(f"duplicate assay-family rule: {rule_id}")
        seen.add(rule_id)
        for field in (
            "preferred_assay_titles",
            "assay_terms",
            "expected_modalities",
        ):
            values = rule.get(field, [])
            if not isinstance(values, list) or any(
                not isinstance(value, str) for value in values
            ):
                raise ValueError(f"assay-family {rule_id} has invalid {field}")
        if not rule["preferred_assay_titles"]:
            raise ValueError(f"assay-family {rule_id} has no exact title")
        result.append(rule)
    return result


def family_ids(row: dict[str, Any], families: list[dict[str, Any]]) -> list[str]:
    titles = set(split_terms(str(row.get("preferred_assay_titles", ""))))
    assay_term = str(row.get("assay_term", ""))
    modality = str(row.get("modality", ""))
    result = []
    for rule in families:
        rule_titles = {str(value) for value in rule.get("preferred_assay_titles", [])}
        assay_terms = {str(value) for value in rule.get("assay_terms", [])}
        modalities = {str(value) for value in rule.get("expected_modalities", [])}
        if not titles.intersection(rule_titles):
            continue
        if assay_terms and assay_term not in assay_terms:
            continue
        if modality and modalities and modality not in modalities:
            continue
        result.append(str(rule["id"]))
    return sorted(result) or ["unclassified"]


def validate_input_rows(
    study_run_id: str,
    runs: list[dict[str, str]],
    failures: list[dict[str, str]],
    diagnostics: list[dict[str, str]],
    metrics: list[dict[str, str]],
) -> None:
    required = {
        "runs": {
            "study_run_id",
            "configuration_accession",
            "modality",
            "access_class",
            "sampled_record_count",
            "run_status",
            "max_attempt_count",
            "preferred_assay_titles",
            "assay_term",
        },
        "failures": {
            "study_run_id",
            "configuration_accession",
            "modality",
            "access_class",
            "failure_category",
            "stage",
            "reason",
            "attempt_count",
            "preferred_assay_titles",
            "assay_term",
        },
        "diagnostics": {
            "study_run_id",
            "configuration_accession",
            "modality",
            "access_class",
            "check",
            "assessment_type",
            "assessment_code",
            "preferred_assay_titles",
            "assay_term",
        },
        "metrics": {
            "study_run_id",
            "configuration_accession",
            "modality",
            "access_class",
            "check",
            "metric_side",
            "metric_name",
            "data_kind",
            "unit",
            "value_json",
            "sequence_types",
            "ontology_terms",
            "preferred_assay_titles",
            "assay_term",
        },
    }
    for name, rows in (
        ("runs", runs),
        ("failures", failures),
        ("diagnostics", diagnostics),
        ("metrics", metrics),
    ):
        if rows and any(not required[name].issubset(row) for row in rows):
            raise ValueError(f"{name} table fields are incomplete")
        if any(row.get("study_run_id") != study_run_id for row in rows):
            raise ValueError(f"{name} rows belong to another study")


def split_terms(value: str) -> list[str]:
    return sorted({item.strip() for item in value.split(";") if item.strip()})


def linear_quantile(values: list[float], probability: float) -> float:
    if not values:
        raise ValueError("cannot compute a quantile of an empty sequence")
    position = (len(values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    fraction = position - lower
    return values[lower] + fraction * (values[upper] - values[lower])


def parse_int(value: Any) -> int:
    if isinstance(value, bool):
        raise ValueError(f"expected an integer, found {value!r}")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as err:
        raise ValueError(f"expected an integer, found {value!r}") from err
    if parsed < 0:
        raise ValueError(f"expected a nonnegative integer, found {parsed}")
    return parsed


def parse_float(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} is not numeric")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as err:
        raise ValueError(f"{label} is not numeric") from err
    if not math.isfinite(parsed):
        raise ValueError(f"{label} is not finite")
    return parsed


def prepare_output_root(path: Path) -> Path:
    path = path.resolve()
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"output root is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as err:
        raise ValueError(f"{label} is not valid JSON: {path}") from err
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
