#!/usr/bin/env python3
"""Select, blind, and merge independent reviews of IGVF audit findings."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

try:
    import analyze_igvf_audit as audit_analysis
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import analyze_igvf_audit as audit_analysis
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
TOOL_VERSION = "0.1.0"
CLASSIFICATIONS = (
    "confirmed_specification_problem",
    "confirmed_read_or_resource_problem",
    "intended_assay_design",
    "seqcheck_problem",
    "inconclusive",
)
CONFIDENCE_VALUES = ("high", "medium", "low")
PACKAGE_FIELDS = ("review_package_id", "review_slot", "review_item_id")
CONTEXT_FIELDS = (
    "configuration_accession",
    "configuration_href",
    "modality",
    "access_class",
    "assay_family_labels",
    "assay_term",
    "preferred_assay_titles",
    "lab",
    "normalized_seqspec_version",
    "check",
    "files",
    "reads",
    "regions",
    "ontology_terms",
    "sequence_types",
    "region_annotations_json",
    "assessment_description",
    "observed_metrics_json",
)
EDITABLE_FIELDS = (
    "reviewer",
    "review_date",
    "classification",
    "confidence",
    "rationale",
)
HIDDEN_FIELDS = (
    "source_kind",
    "assessment_type",
    "assessment_code",
    "high_confidence",
    "high_confidence_basis",
    "matched_assay_family",
    "metric_value_stratum",
    "finding_key",
    "source_row_count",
)
MERGED_REVIEW_FIELDS = (
    "reviewer_1",
    "reviewer_1_date",
    "reviewer_1_classification",
    "reviewer_1_confidence",
    "reviewer_1_rationale",
    "reviewer_2",
    "reviewer_2_date",
    "reviewer_2_classification",
    "reviewer_2_confidence",
    "reviewer_2_rationale",
    "exact_classification_agreement",
    "needs_adjudication",
    "consensus_classification",
    "adjudicator",
    "adjudication_date",
    "final_classification",
    "adjudication_rationale",
    "submitter_contacted",
    "submitter_response",
    "proposed_correction",
    "final_status_note",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare or merge blinded reviews of IGVF audit findings."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    add_source_arguments(prepare)
    prepare.add_argument("--output-root", required=True, type=Path)
    merge = subparsers.add_parser("merge")
    add_source_arguments(merge)
    merge.add_argument("--reviewer-1-package", required=True, type=Path)
    merge.add_argument("--reviewer-1-sheet", required=True, type=Path)
    merge.add_argument("--reviewer-2-package", required=True, type=Path)
    merge.add_argument("--reviewer-2-sheet", required=True, type=Path)
    merge.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def add_source_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--audit-root", required=True, type=Path)
    parser.add_argument("--audit-analysis", required=True, type=Path)
    parser.add_argument("--detection-policy", required=True, type=Path)
    parser.add_argument("--family-rules", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)


def main() -> int:
    args = parse_args()
    common = {
        "audit_root": args.audit_root.resolve(),
        "audit_analysis_path": args.audit_analysis.resolve(),
        "detection_policy_path": args.detection_policy.resolve(),
        "family_rules_path": args.family_rules.resolve(),
        "protocol_path": args.protocol.resolve(),
    }
    try:
        if args.command == "prepare":
            prepare_review_packages(output_root=args.output_root.resolve(), **common)
        else:
            merge_review_packages(
                reviewer_1_package=args.reviewer_1_package.resolve(),
                reviewer_1_sheet=args.reviewer_1_sheet.resolve(),
                reviewer_2_package=args.reviewer_2_package.resolve(),
                reviewer_2_sheet=args.reviewer_2_sheet.resolve(),
                output_root=args.output_root.resolve(),
                **common,
            )
    except (OSError, ValueError) as error:
        print(f"manage_audit_reviews: {error}", file=sys.stderr)
        return 1
    return 0


def prepare_review_packages(
    *,
    audit_root: Path,
    audit_analysis_path: Path,
    detection_policy_path: Path,
    family_rules_path: Path,
    protocol_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    refuse_output_root(output_root)
    source = load_source(
        audit_root=audit_root,
        audit_analysis_path=audit_analysis_path,
        detection_policy_path=detection_policy_path,
        family_rules_path=family_rules_path,
        protocol_path=protocol_path,
    )
    tool = source["tool"]
    package_manifests = []
    item_ids_by_slot = {}
    for slot in (1, 2):
        package, item_ids = prepare_one_package(
            slot=slot,
            source=source,
            output_root=output_root,
        )
        package_manifests.append(package)
        item_ids_by_slot[slot] = item_ids

    private_rows = []
    for row in source["selected_cases"]:
        private_rows.append(
            {
                "selection_id": source["selection_id"],
                "case_id": row["case_id"],
                "reviewer_1_item_id": item_ids_by_slot[1][row["case_id"]],
                "reviewer_2_item_id": item_ids_by_slot[2][row["case_id"]],
                **{field: row[field] for field in HIDDEN_FIELDS},
                **{field: row[field] for field in CONTEXT_FIELDS},
            }
        )
    private_fields = (
        "selection_id",
        "case_id",
        "reviewer_1_item_id",
        "reviewer_2_item_id",
        *HIDDEN_FIELDS,
        *CONTEXT_FIELDS,
    )
    private_path = output_root / "internal" / "review_key.csv"
    runtime.write_csv(private_path, private_rows, list(private_fields))
    selected_path = output_root / "internal" / "selected_cases.json"
    runtime.write_json(
        selected_path,
        {
            "schema_version": SCHEMA_VERSION,
            "selection_id": source["selection_id"],
            "cases": source["selected_cases"],
        },
    )
    stable = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": source["selection_id"],
        "tool": runtime.functional_script_identity(tool),
        "packages": [
            {"review_slot": value["review_slot"], "package_id": value["package_id"]}
            for value in package_manifests
        ],
    }
    preparation_id = runtime.sha256_json(stable)[:16]
    manifest_path = output_root / "manifests" / "audit_review_packages.json"
    manifest = {
        **stable,
        "preparation_id": preparation_id,
        "generated_at": utc_now(),
        "review_rows": len(source["selected_cases"]),
        "candidate_rows": sum(
            row["source_kind"] == "candidate" for row in source["selected_cases"]
        ),
        "pass_control_rows": sum(
            row["source_kind"] == "pass_control" for row in source["selected_cases"]
        ),
        "high_confidence_candidate_rows": sum(
            row["source_kind"] == "candidate" and row["high_confidence"] == "true"
            for row in source["selected_cases"]
        ),
        "review_slots": 2,
        "blinded_fields": [*HIDDEN_FIELDS],
        "inputs": source["input_identities"],
        "packages": [
            {
                "review_slot": value["review_slot"],
                "package_id": value["package_id"],
                "manifest": runtime.file_identity(Path(value["manifest_path"])),
            }
            for value in package_manifests
        ],
        "private_review_key": runtime.file_identity(private_path),
        "private_selected_cases": runtime.file_identity(selected_path),
        "private_review_key_sharing": "do_not_share_with_reviewers",
        "manifest_path": str(manifest_path),
        "frozen": False,
    }
    runtime.write_json(manifest_path, manifest)
    print(
        f"prepared two blinded audit review packages for "
        f"{len(source['selected_cases'])} cases "
        f"(selection_id={source['selection_id']})"
    )
    return manifest


def prepare_one_package(
    *, slot: int, source: dict[str, Any], output_root: Path
) -> tuple[dict[str, Any], dict[str, str]]:
    seed = source["package_seeds"][slot]
    identity = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": source["selection_id"],
        "review_slot": slot,
        "package_seed": seed,
        "tool": runtime.functional_script_identity(source["tool"]),
    }
    package_id = runtime.sha256_json(identity)[:16]
    rows, item_ids = expected_package_rows(
        slot=slot,
        seed=seed,
        package_id=package_id,
        source=source,
    )
    package_dir = output_root / f"reviewer_{slot}"
    sheet_path = package_dir / "audit_review.csv"
    instructions_path = package_dir / "INSTRUCTIONS.md"
    manifest_path = package_dir / "review_package.json"
    runtime.write_csv(
        sheet_path,
        rows,
        [*PACKAGE_FIELDS, *CONTEXT_FIELDS, *EDITABLE_FIELDS],
    )
    instructions_path.parent.mkdir(parents=True, exist_ok=True)
    instructions_path.write_text(review_instructions(slot), encoding="utf-8")
    evidence_sha256 = evidence_hash(rows)
    manifest = {
        **identity,
        "package_id": package_id,
        "generated_at": utc_now(),
        "review_rows": len(rows),
        "blinded_fields": [*HIDDEN_FIELDS],
        "fields": {
            "package": [*PACKAGE_FIELDS],
            "context": [*CONTEXT_FIELDS],
            "editable": [*EDITABLE_FIELDS],
        },
        "review_evidence_sha256": evidence_sha256,
        "inputs": source["input_identities"],
        "prepared_sheet": runtime.file_identity(sheet_path),
        "instructions": runtime.file_identity(instructions_path),
        "manifest_path": str(manifest_path),
        "frozen": False,
    }
    runtime.write_json(manifest_path, manifest)
    return manifest, item_ids


def merge_review_packages(
    *,
    audit_root: Path,
    audit_analysis_path: Path,
    detection_policy_path: Path,
    family_rules_path: Path,
    protocol_path: Path,
    reviewer_1_package: Path,
    reviewer_1_sheet: Path,
    reviewer_2_package: Path,
    reviewer_2_sheet: Path,
    output_root: Path,
) -> dict[str, Any]:
    refuse_output_root(output_root)
    source = load_source(
        audit_root=audit_root,
        audit_analysis_path=audit_analysis_path,
        detection_policy_path=detection_policy_path,
        family_rules_path=family_rules_path,
        protocol_path=protocol_path,
    )
    first = load_completed_review(
        package_path=reviewer_1_package,
        sheet_path=reviewer_1_sheet,
        expected_slot=1,
        source=source,
    )
    second = load_completed_review(
        package_path=reviewer_2_package,
        sheet_path=reviewer_2_sheet,
        expected_slot=2,
        source=source,
    )
    if first["reviewer"].casefold() == second["reviewer"].casefold():
        raise ValueError("audit reviewers must be distinct people")

    merged_rows = []
    for case in source["selected_cases"]:
        left = first["rows_by_case_id"][case["case_id"]]
        right = second["rows_by_case_id"][case["case_id"]]
        agreement = left["classification"] == right["classification"]
        consensus = left["classification"] if agreement else ""
        merged_rows.append(
            {
                "selection_id": source["selection_id"],
                **case,
                **prefixed_review(left, 1),
                **prefixed_review(right, 2),
                "exact_classification_agreement": agreement,
                "needs_adjudication": not agreement,
                "consensus_classification": consensus,
                "adjudicator": "",
                "adjudication_date": "",
                "final_classification": consensus,
                "adjudication_rationale": "",
                "submitter_contacted": "",
                "submitter_response": "",
                "proposed_correction": "",
                "final_status_note": "",
            }
        )
    output_fields = (
        "selection_id",
        "case_id",
        *HIDDEN_FIELDS,
        *CONTEXT_FIELDS,
        *MERGED_REVIEW_FIELDS,
    )
    merged_path = output_root / "tables" / "audit_reviews.csv"
    runtime.write_csv(merged_path, merged_rows, list(output_fields))
    disagreements = sum(row["needs_adjudication"] for row in merged_rows)
    stable = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": source["selection_id"],
        "tool": runtime.functional_script_identity(source["tool"]),
        "reviewer_1_package_sha256": runtime.file_sha256(reviewer_1_package),
        "reviewer_1_sheet_sha256": runtime.file_sha256(reviewer_1_sheet),
        "reviewer_2_package_sha256": runtime.file_sha256(reviewer_2_package),
        "reviewer_2_sheet_sha256": runtime.file_sha256(reviewer_2_sheet),
    }
    merge_id = runtime.sha256_json(stable)[:16]
    manifest_path = output_root / "manifests" / "audit_reviews.json"
    manifest = {
        **stable,
        "merge_id": merge_id,
        "generated_at": utc_now(),
        "review_rows": len(merged_rows),
        "reviewers": {
            "reviewer_1": first["reviewer"],
            "reviewer_2": second["reviewer"],
        },
        "exact_agreements": len(merged_rows) - disagreements,
        "adjudication_required": disagreements,
        "adjudication_complete": disagreements == 0,
        "frozen": False,
        "inputs": {
            **source["input_identities"],
            "reviewer_1_package": runtime.file_identity(reviewer_1_package),
            "reviewer_1_sheet": runtime.file_identity(reviewer_1_sheet),
            "reviewer_2_package": runtime.file_identity(reviewer_2_package),
            "reviewer_2_sheet": runtime.file_identity(reviewer_2_sheet),
        },
        "outputs": {"audit_reviews": runtime.file_identity(merged_path)},
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    print(
        f"merged {len(merged_rows)} blinded audit reviews "
        f"(merge_id={merge_id}, adjudication_required={disagreements})"
    )
    return manifest


def load_source(
    *,
    audit_root: Path,
    audit_analysis_path: Path,
    detection_policy_path: Path,
    family_rules_path: Path,
    protocol_path: Path,
) -> dict[str, Any]:
    audit_root = audit_root.resolve()
    paths = {
        "study": audit_root / "manifests" / "study.json",
        "reconciliation": audit_root / "validation" / "reconciliation.json",
        "runs": audit_root / "runs.csv",
        "failures": audit_root / "failures.csv",
        "diagnostics": audit_root / "diagnostics.csv",
        "metrics": audit_root / "metrics.csv",
    }
    study = runtime.load_json(paths["study"])
    reconciliation = runtime.load_json(paths["reconciliation"])
    analysis = runtime.load_json(audit_analysis_path)
    if reconciliation.get("valid") is not True or analysis.get("valid") is not True:
        raise ValueError("audit and audit analysis must both be valid")
    if analysis.get("schema_version") != audit_analysis.SCHEMA_VERSION:
        raise ValueError("audit analysis schema is unsupported")
    study_run_id = str(study.get("run_id", ""))
    if (
        not study_run_id
        or reconciliation.get("study_run_id") != study_run_id
        or analysis.get("study_run_id") != study_run_id
    ):
        raise ValueError("audit study identifiers differ")
    analysis_stable = {
        key: analysis.get(key)
        for key in (
            "schema_version",
            "study_run_id",
            "audit_protocol_id",
            "tool",
            "inputs",
        )
    }
    if runtime.sha256_json(analysis_stable)[:16] != analysis.get("analysis_id"):
        raise ValueError("audit analysis id is not content-addressed")
    for name, path in paths.items():
        expected = analysis.get("inputs", {}).get(name)
        if expected != runtime.file_sha256(path):
            raise ValueError(f"audit analysis {name} input changed")
    if analysis.get("inputs", {}).get("family_rules") != runtime.file_sha256(
        family_rules_path
    ):
        raise ValueError("audit analysis family rules changed")

    runs = audit_analysis.read_csv(paths["runs"])
    failures = audit_analysis.read_csv(paths["failures"])
    diagnostics = audit_analysis.read_csv(paths["diagnostics"])
    metrics = audit_analysis.read_csv(paths["metrics"])
    audit_analysis.validate_input_rows(
        study_run_id, runs, failures, diagnostics, metrics
    )
    if not diagnostics:
        raise ValueError("audit has no diagnostic rows to review")
    run_index = index_unique(
        runs, ("configuration_accession", "modality"), "audit runs"
    )
    if any(
        (row["configuration_accession"], row["modality"]) not in run_index
        for row in diagnostics
    ):
        raise ValueError("audit diagnostics do not match completed runs")

    policy = load_detection_policy(detection_policy_path)
    family_rules = audit_analysis.load_family_rules(family_rules_path)
    protocol = load_protocol(protocol_path)
    tool = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    high_confidence_codes = sorted(
        {
            str(code)
            for entry in policy["entries"]
            for code in entry.get("assessment_codes", [])
            if str(code)
        }
    )
    findings = build_findings(
        study=study,
        diagnostics=diagnostics,
        metrics=metrics,
        families=family_rules,
        high_confidence_codes=set(high_confidence_codes),
        protocol=protocol,
    )
    selected = select_cases(findings, protocol)
    stable_inputs = {
        "audit_analysis": runtime.file_sha256(audit_analysis_path),
        "detection_policy": runtime.file_sha256(detection_policy_path),
        "family_rules": runtime.file_sha256(family_rules_path),
        "protocol": runtime.file_sha256(protocol_path),
        **{name: runtime.file_sha256(path) for name, path in paths.items()},
    }
    selection_stable = {
        "schema_version": SCHEMA_VERSION,
        "study_run_id": study_run_id,
        "audit_analysis_id": analysis.get("analysis_id"),
        "detection_policy_id": policy["policy_id"],
        "tool": runtime.functional_script_identity(tool),
        "inputs": stable_inputs,
        "selected": [
            {
                "finding_key": row["finding_key"],
                "source_kind": row["source_kind"],
                "high_confidence": row["high_confidence"],
                "matched_assay_family": row["matched_assay_family"],
            }
            for row in selected
        ],
    }
    selection_id = runtime.sha256_json(selection_stable)[:16]
    for row in selected:
        row["case_id"] = runtime.sha256_json(
            [selection_id, row["finding_key"], row["source_kind"]]
        )[:16]
    input_identities = {
        "audit_analysis": runtime.file_identity(audit_analysis_path),
        "detection_policy": runtime.file_identity(detection_policy_path),
        "family_rules": runtime.file_identity(family_rules_path),
        "protocol": runtime.file_identity(protocol_path),
        **{name: runtime.file_identity(path) for name, path in paths.items()},
    }
    return {
        "selection_id": selection_id,
        "study_run_id": study_run_id,
        "selected_cases": selected,
        "package_seeds": protocol["review"]["package_seeds"],
        "classifications": set(protocol["review"]["classifications"]),
        "confidence_values": set(protocol["review"]["confidence_values"]),
        "tool": tool,
        "input_identities": input_identities,
    }


def build_findings(
    *,
    study: dict[str, Any],
    diagnostics: list[dict[str, str]],
    metrics: list[dict[str, str]],
    families: list[dict[str, Any]],
    high_confidence_codes: set[str],
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    configuration_hrefs = {
        str(row.get("accession", "")): str(row.get("href", ""))
        for row in study.get("configurations", [])
        if isinstance(row, dict)
    }
    metric_index: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in metrics:
        if row.get("metric_side") != "observed":
            continue
        metric_index[scope_key(row)].append(row)
    eligible_high_types = set(protocol["high_confidence"]["eligible_assessment_types"])
    family_labels = {row["id"]: row["label"] for row in families}
    family_labels["unclassified"] = "Unclassified assay family"
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in diagnostics:
        key_payload = {field: row.get(field, "") for field in finding_key_fields()}
        finding_key = runtime.sha256_json(key_payload)
        grouped[finding_key].append(row)

    findings = []
    for finding_key, rows in sorted(grouped.items()):
        row = rows[0]
        family_ids = audit_analysis.family_ids(row, families)
        primary_family = family_ids[0]
        observed = [
            {
                "metric_name": value["metric_name"],
                "unit": value["unit"],
                "data_kind": value["data_kind"],
                "value": json.loads(value["value_json"]),
            }
            for value in sorted(
                metric_index.get(scope_key(row), []),
                key=lambda value: (
                    value["metric_name"],
                    value["unit"],
                    value["value_json"],
                ),
            )
        ]
        high_confidence = (
            row["assessment_type"] in eligible_high_types
            and row["assessment_code"] in high_confidence_codes
        )
        findings.append(
            {
                "finding_key": finding_key,
                "source_row_count": str(len(rows)),
                "assessment_type": row["assessment_type"],
                "assessment_code": row["assessment_code"],
                "high_confidence": bool_text(high_confidence),
                "high_confidence_basis": (
                    "controlled_policy_assessment_code" if high_confidence else ""
                ),
                "metric_value_stratum": metric_value_stratum(observed),
                "assay_family_ids": family_ids,
                "primary_assay_family": primary_family,
                "configuration_accession": row["configuration_accession"],
                "configuration_href": configuration_hrefs.get(
                    row["configuration_accession"], ""
                ),
                "modality": row["modality"],
                "access_class": row["access_class"],
                "assay_family_labels": ";".join(
                    family_labels[value] for value in family_ids
                ),
                "assay_term": row["assay_term"],
                "preferred_assay_titles": row["preferred_assay_titles"],
                "lab": row.get("lab", ""),
                "normalized_seqspec_version": row.get("normalized_seqspec_version", ""),
                "check": row["check"],
                "files": row.get("files", ""),
                "reads": row.get("reads", ""),
                "regions": row.get("regions", ""),
                "ontology_terms": row.get("ontology_terms", ""),
                "sequence_types": row.get("sequence_types", ""),
                "region_annotations_json": row.get("region_annotations_json", "[]"),
                "assessment_description": row.get("assessment_description", ""),
                "observed_metrics_json": runtime.canonical_json(observed),
            }
        )
    return findings


def select_cases(
    findings: list[dict[str, Any]], protocol: dict[str, Any]
) -> list[dict[str, Any]]:
    selection = protocol["selection"]
    candidate_types = set(selection["candidate_assessment_types"])
    pass_type = selection["pass_control_assessment_type"]
    candidates = [row for row in findings if row["assessment_type"] in candidate_types]
    controls = [row for row in findings if row["assessment_type"] == pass_type]
    target = selection["minimum_candidate_findings"]
    control_target = selection["pass_control_findings"]
    if len(candidates) < target:
        raise ValueError(
            f"audit has {len(candidates)} candidate findings; {target} are required"
        )
    mandatory_types = set(selection["include_all_assessment_types"])
    mandatory = [
        row
        for row in candidates
        if row["assessment_type"] in mandatory_types
        or (
            protocol["high_confidence"]["include_all_high_confidence"]
            and row["high_confidence"] == "true"
        )
    ]
    selected_candidates = select_stratified(
        candidates=candidates,
        mandatory=mandatory,
        target=max(target, len(mandatory)),
        fields=selection["stratification_fields"],
        seed=selection["selection_seed"],
    )
    for row in selected_candidates:
        row["source_kind"] = "candidate"
        row["matched_assay_family"] = row["primary_assay_family"]
    selected_controls = select_matched_controls(
        controls=controls,
        candidates=selected_candidates,
        target=control_target,
        seed=selection["selection_seed"],
    )
    return sorted(
        [*selected_candidates, *selected_controls],
        key=lambda row: (row["source_kind"], row["finding_key"]),
    )


def select_stratified(
    *,
    candidates: list[dict[str, Any]],
    mandatory: list[dict[str, Any]],
    target: int,
    fields: list[str],
    seed: int,
) -> list[dict[str, Any]]:
    selected = {row["finding_key"]: dict(row) for row in mandatory}
    grouped: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in candidates:
        if row["finding_key"] in selected:
            continue
        grouped[tuple(str(row[field]) for field in fields)].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: stable_rank(seed, "candidate", row["finding_key"]))
    strata = sorted(grouped, key=lambda value: stable_rank(seed, "stratum", *value))
    while len(selected) < target:
        added = False
        for stratum in strata:
            if grouped[stratum] and len(selected) < target:
                row = grouped[stratum].pop(0)
                selected[row["finding_key"]] = dict(row)
                added = True
        if not added:
            break
    if len(selected) < target:
        raise ValueError("stratified candidate selection did not reach its target")
    return sorted(
        selected.values(),
        key=lambda row: stable_rank(seed, "selected-candidate", row["finding_key"]),
    )


def select_matched_controls(
    *,
    controls: list[dict[str, Any]],
    candidates: list[dict[str, Any]],
    target: int,
    seed: int,
) -> list[dict[str, Any]]:
    candidate_counts = Counter(row["primary_assay_family"] for row in candidates)
    allocations = proportional_allocations(candidate_counts, target)
    available = {
        family: sorted(
            [row for row in controls if family in row["assay_family_ids"]],
            key=lambda row: stable_rank(seed, "control", family, row["finding_key"]),
        )
        for family in allocations
    }
    selected = []
    used = set()
    order = sorted(
        allocations,
        key=lambda family: (
            len(available[family]) / allocations[family]
            if allocations[family]
            else float("inf"),
            family,
        ),
    )
    for family in order:
        chosen = [row for row in available[family] if row["finding_key"] not in used][
            : allocations[family]
        ]
        if len(chosen) != allocations[family]:
            raise ValueError(f"not enough pass controls matched assay family {family}")
        for row in chosen:
            value = dict(row)
            value["source_kind"] = "pass_control"
            value["matched_assay_family"] = family
            value["high_confidence"] = "false"
            value["high_confidence_basis"] = ""
            selected.append(value)
            used.add(row["finding_key"])
    if len(selected) != target:
        raise ValueError("matched pass-control selection did not reach its target")
    return selected


def proportional_allocations(counts: Counter[str], target: int) -> dict[str, int]:
    if target <= 0 or not counts:
        raise ValueError(
            "pass-control allocation target and candidate families are required"
        )
    total = sum(counts.values())
    raw = {key: target * value / total for key, value in counts.items()}
    result = {key: int(value) for key, value in raw.items()}
    remaining = target - sum(result.values())
    order = sorted(counts, key=lambda key: (-(raw[key] - result[key]), key))
    for key in order[:remaining]:
        result[key] += 1
    return {key: value for key, value in result.items() if value}


def load_completed_review(
    *,
    package_path: Path,
    sheet_path: Path,
    expected_slot: int,
    source: dict[str, Any],
) -> dict[str, Any]:
    package = runtime.load_json(package_path)
    seed = source["package_seeds"][expected_slot]
    expected_identity = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": source["selection_id"],
        "review_slot": expected_slot,
        "package_seed": seed,
        "tool": runtime.functional_script_identity(source["tool"]),
    }
    expected_package_id = runtime.sha256_json(expected_identity)[:16]
    if any(package.get(key) != value for key, value in expected_identity.items()):
        raise ValueError(f"reviewer {expected_slot} package identity changed")
    if package.get("package_id") != expected_package_id:
        raise ValueError(f"reviewer {expected_slot} package id changed")
    rows, fields = read_csv(sheet_path)
    expected_fields = [*PACKAGE_FIELDS, *CONTEXT_FIELDS, *EDITABLE_FIELDS]
    if fields != expected_fields:
        raise ValueError(f"reviewer {expected_slot} sheet fields changed")
    if len(rows) != len(source["selected_cases"]):
        raise ValueError(f"reviewer {expected_slot} sheet row count changed")
    expected_rows, item_ids = expected_package_rows(
        slot=expected_slot,
        seed=seed,
        package_id=expected_package_id,
        source=source,
    )
    expected_by_item = {row["review_item_id"]: row for row in expected_rows}
    case_id_by_item = {item_id: case_id for case_id, item_id in item_ids.items()}
    observed_by_item = index_unique(
        rows, ("review_item_id",), f"reviewer {expected_slot} rows"
    )
    observed_by_item = {key[0]: value for key, value in observed_by_item.items()}
    if set(observed_by_item) != set(expected_by_item):
        raise ValueError(f"reviewer {expected_slot} item ids changed")
    reviewers = set()
    dates = set()
    rows_by_case_id = {}
    for item_id, expected in expected_by_item.items():
        row = observed_by_item[item_id]
        for field in (*PACKAGE_FIELDS, *CONTEXT_FIELDS):
            if row[field] != str(expected[field]):
                raise ValueError(
                    f"reviewer {expected_slot} changed blinded evidence: {field}"
                )
        reviewer = row["reviewer"].strip()
        review_date = row["review_date"].strip()
        classification = row["classification"].strip()
        confidence = row["confidence"].strip()
        rationale = row["rationale"].strip()
        if not reviewer or not rationale:
            raise ValueError(f"reviewer {expected_slot} response is incomplete")
        validate_date(review_date, f"reviewer {expected_slot}")
        if classification not in source["classifications"]:
            raise ValueError(f"reviewer {expected_slot} classification is invalid")
        if confidence not in source["confidence_values"]:
            raise ValueError(f"reviewer {expected_slot} confidence is invalid")
        reviewers.add(reviewer)
        dates.add(review_date)
        rows_by_case_id[case_id_by_item[item_id]] = {
            **row,
            "reviewer": reviewer,
            "review_date": review_date,
            "classification": classification,
            "confidence": confidence,
            "rationale": rationale,
        }
    if len(reviewers) != 1 or len(dates) != 1:
        raise ValueError(f"reviewer {expected_slot} name and date must be consistent")
    if evidence_hash(rows) != package.get("review_evidence_sha256"):
        raise ValueError(f"reviewer {expected_slot} blinded evidence hash changed")
    return {"reviewer": next(iter(reviewers)), "rows_by_case_id": rows_by_case_id}


def expected_package_rows(
    *, slot: int, seed: int, package_id: str, source: dict[str, Any]
) -> tuple[list[dict[str, str]], dict[str, str]]:
    rows = []
    item_ids = {}
    for case in source["selected_cases"]:
        item_id = runtime.sha256_json(
            [source["selection_id"], case["case_id"], slot, seed]
        )[:16]
        item_ids[case["case_id"]] = item_id
        rows.append(
            {
                "review_package_id": package_id,
                "review_slot": str(slot),
                "review_item_id": item_id,
                **{field: str(case[field]) for field in CONTEXT_FIELDS},
                **{field: "" for field in EDITABLE_FIELDS},
            }
        )
    rows.sort(key=lambda row: stable_rank(seed, "package-order", row["review_item_id"]))
    return rows, item_ids


def evidence_hash(rows: list[dict[str, str]]) -> str:
    fields = [*PACKAGE_FIELDS, *CONTEXT_FIELDS]
    return runtime.sha256_json(
        [{field: row[field] for field in fields} for row in rows]
    )


def prefixed_review(row: dict[str, str], slot: int) -> dict[str, str]:
    return {
        f"reviewer_{slot}": row["reviewer"],
        f"reviewer_{slot}_date": row["review_date"],
        f"reviewer_{slot}_classification": row["classification"],
        f"reviewer_{slot}_confidence": row["confidence"],
        f"reviewer_{slot}_rationale": row["rationale"],
    }


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
    if any(
        not isinstance(entry, dict)
        or entry.get("ready") is not True
        or not isinstance(entry.get("assessment_codes", []), list)
        for entry in entries
    ):
        raise ValueError("detection policy entries are invalid")
    return policy


def load_protocol(path: Path) -> dict[str, Any]:
    protocol = runtime.load_json(path)
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("audit adjudication protocol schema is unsupported")
    if protocol.get("experiment_id") != "audit-adjudication":
        raise ValueError("audit adjudication experiment id is invalid")
    selection = protocol.get("selection", {})
    for field in ("minimum_candidate_findings", "pass_control_findings"):
        positive_int(selection.get(field), field)
    nonnegative_int(selection.get("selection_seed"), "selection seed")
    candidate_types = selection.get("candidate_assessment_types")
    if candidate_types != ["error", "warning", "interpretation"]:
        raise ValueError("candidate assessment types are invalid")
    if selection.get("pass_control_assessment_type") != "pass":
        raise ValueError("pass-control assessment type is invalid")
    include_all = selection.get("include_all_assessment_types")
    if not isinstance(include_all, list) or not set(include_all) <= set(
        candidate_types
    ):
        raise ValueError("include-all assessment types are invalid")
    fields = selection.get("stratification_fields")
    allowed_fields = {
        "assessment_type",
        "check",
        "primary_assay_family",
        "ontology_terms",
        "metric_value_stratum",
        "lab",
    }
    if not isinstance(fields, list) or not fields or not set(fields) <= allowed_fields:
        raise ValueError("candidate stratification fields are invalid")
    if selection.get("control_matching_field") != "assay_family":
        raise ValueError("pass-control matching field is invalid")

    high = protocol.get("high_confidence", {})
    if high.get("source") != ("frozen_controlled_perturbation_policy_assessment_codes"):
        raise ValueError("high-confidence source is invalid")
    if high.get("eligible_assessment_types") != ["error", "warning"]:
        raise ValueError("high-confidence assessment types are invalid")
    if high.get("include_all_high_confidence") is not True:
        raise ValueError("high-confidence findings must all be selected")

    review = protocol.get("review", {})
    if (
        review.get("reviewers") != 2
        or review.get("blind_source_kind") is not True
        or review.get("blind_original_assessment") is not True
        or review.get("blind_high_confidence") is not True
    ):
        raise ValueError("audit reviews must use two fully blinded reviewers")
    if review.get("classifications") != list(CLASSIFICATIONS):
        raise ValueError("audit review classifications are invalid")
    if review.get("confidence_values") != list(CONFIDENCE_VALUES):
        raise ValueError("audit review confidence values are invalid")
    if review.get("adjudication_required_for") != "classification_disagreement":
        raise ValueError("audit review adjudication rule is invalid")
    seeds = review.get("package_seeds", {})
    if set(seeds) != {"reviewer_1", "reviewer_2"}:
        raise ValueError("audit review package seeds are invalid")
    for value in seeds.values():
        nonnegative_int(value, "package seed")
    if seeds["reviewer_1"] == seeds["reviewer_2"]:
        raise ValueError("audit review package seeds must differ")
    review["package_seeds"] = {1: seeds["reviewer_1"], 2: seeds["reviewer_2"]}

    analysis = protocol.get("analysis", {})
    if analysis.get("confirmed_classifications") != list(CLASSIFICATIONS[:2]):
        raise ValueError("confirmed audit classifications are invalid")
    if analysis.get("inconclusive_classification") != "inconclusive":
        raise ValueError("audit inconclusive classification is invalid")
    if analysis.get("pass_consistent_classifications") != ["intended_assay_design"]:
        raise ValueError("pass-control consistency classifications are invalid")
    if (
        analysis.get("agreement_endpoint") != "unweighted_cohen_kappa"
        or analysis.get("confidence_interval") != "wilson"
    ):
        raise ValueError("audit review analysis endpoint is invalid")
    confidence = probability(analysis.get("confidence_level"), "confidence level")
    if not 0 < confidence < 1:
        raise ValueError("confidence level must be inside (0, 1)")
    targets = protocol.get("scientific_targets", {})
    for field in (
        "minimum_high_confidence_confirmation_precision",
        "minimum_cohen_kappa",
        "minimum_pass_control_consistency",
    ):
        probability(targets.get(field), field)
    for field in (
        "minimum_presentable_confirmed_cases",
        "minimum_downstream_cases_with_expected_change",
    ):
        positive_int(targets.get(field), field)
    return protocol


def metric_value_stratum(observed: list[dict[str, Any]]) -> str:
    strata = set()
    for metric in observed:
        value = metric.get("value")
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            strata.add("nonnumeric")
            continue
        magnitude = abs(float(value))
        if magnitude == 0:
            strata.add("zero")
        elif magnitude < 0.001:
            strata.add("lt_1e-3")
        elif magnitude < 0.01:
            strata.add("1e-3_to_1e-2")
        elif magnitude < 0.1:
            strata.add("1e-2_to_1e-1")
        elif magnitude < 1:
            strata.add("1e-1_to_1")
        elif magnitude < 10:
            strata.add("1_to_10")
        elif magnitude < 100:
            strata.add("10_to_100")
        else:
            strata.add("ge_100")
    return ";".join(sorted(strata)) or "no_observed_metric"


def finding_key_fields() -> tuple[str, ...]:
    return (
        "configuration_accession",
        "modality",
        "check",
        "files",
        "reads",
        "regions",
        "ontology_terms",
        "sequence_types",
        "region_annotations_json",
        "assessment_type",
        "assessment_code",
        "assessment_description",
    )


def scope_key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(
        row.get(field, "")
        for field in (
            "configuration_accession",
            "modality",
            "check",
            "files",
            "reads",
            "regions",
        )
    )


def index_unique(
    rows: list[dict[str, str]], fields: tuple[str, ...], label: str
) -> dict[tuple[str, ...], dict[str, str]]:
    indexed = {}
    for row in rows:
        key = tuple(row.get(field, "") for field in fields)
        if not all(key) or key in indexed:
            raise ValueError(f"{label} keys are incomplete or duplicated")
        indexed[key] = row
    return indexed


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        return list(reader), fields


def review_instructions(slot: int) -> str:
    labels = "\n".join(f"- `{value}`" for value in CLASSIFICATIONS)
    return f"""# Audit finding review {slot}

Review every row independently. Candidate/pass-control status, the original
seqcheck assessment, and the policy-derived high-confidence flag are hidden.

1. Use the assay, seqspec region context, observed metrics, and linked portal
   record to classify the case.
2. Enter one classification from this list:
{labels}
3. Set `confidence` to `high`, `medium`, or `low`.
4. Give a concrete rationale. Do not infer the hidden source class.
5. Use one reviewer name and one ISO date (`YYYY-MM-DD`) throughout the sheet.

Do not compare sheets with the other reviewer before both reviews are complete.
"""


def stable_rank(seed: int, *values: str) -> str:
    return runtime.sha256_json([seed, *values])


def validate_date(value: str, label: str) -> None:
    try:
        parsed = date.fromisoformat(value)
    except ValueError as error:
        raise ValueError(f"{label} date is not ISO YYYY-MM-DD") from error
    if parsed.isoformat() != value:
        raise ValueError(f"{label} date is not canonical ISO YYYY-MM-DD")


def positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def nonnegative_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def probability(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    parsed = float(value)
    if not 0 <= parsed <= 1:
        raise ValueError(f"{label} must be between zero and one")
    return parsed


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def refuse_output_root(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"output root is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
