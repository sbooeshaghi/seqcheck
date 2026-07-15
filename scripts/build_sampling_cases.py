#!/usr/bin/env python3
"""Select ontology-diverse FASTQs from one frozen cohort split."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import urllib.parse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "0.1.0"
SELECTOR_VERSION = "0.1.0"
COHORT_SPLITS = ("calibration", "evaluation")
ROLE_PREFIXES = {
    "measure": "RGN:measure:",
    "partition": "RGN:partition:",
    "technical": "RGN:technical:",
}
CASE_FIELDS = (
    "selection_id",
    "case_id",
    "family_id",
    "configuration_accession",
    "modality",
    "selection_role",
    "read_id",
    "fastq_accession",
    "fastq_url",
    "declared_compressed_bytes",
    "indexed_bases",
    "measure_bases",
    "partition_bases",
    "technical_bases",
    "role_classes",
    "ontology_terms",
    "spec_sha256",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select one or two ontology-diverse FASTQs from every configuration "
            "in one frozen cohort split."
        )
    )
    parser.add_argument("--cohort-manifest", required=True, type=Path)
    parser.add_argument("--sampling-protocol", required=True, type=Path)
    parser.add_argument("--seqspec-bin", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument(
        "--cohort-split",
        choices=COHORT_SPLITS,
        default="calibration",
    )
    parser.add_argument("--timeout-seconds", type=int, default=120)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = build_sampling_cases(
            cohort_manifest_path=args.cohort_manifest.resolve(),
            protocol_path=args.sampling_protocol.resolve(),
            seqspec_bin=args.seqspec_bin.resolve(),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
            cohort_split=args.cohort_split,
        )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"build_sampling_cases: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def build_sampling_cases(
    *,
    cohort_manifest_path: Path,
    protocol_path: Path,
    seqspec_bin: Path,
    output_root: Path,
    timeout_seconds: int,
    cohort_split: str = "calibration",
) -> dict[str, Any]:
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")

    cohort_manifest = load_json(cohort_manifest_path)
    protocol = load_json(protocol_path)
    cohort_path, cohort_rows = validate_inputs(
        cohort_manifest,
        protocol,
        seqspec_bin,
        cohort_split,
    )
    selection = protocol["case_selection"]
    max_fastqs = int(selection["max_fastqs_per_configuration"])
    selected_rows = []
    configuration_records = []
    for row in sorted(
        (value for value in cohort_rows if value.get("split") == cohort_split),
        key=lambda value: (
            value.get("final_family", ""),
            value.get("configuration_accession", ""),
        ),
    ):
        record = prepare_configuration(
            row=row,
            seqspec_bin=seqspec_bin,
            timeout_seconds=timeout_seconds,
            max_fastqs=max_fastqs,
            require_positive_size=bool(
                selection["require_positive_declared_compressed_size"]
            ),
        )
        configuration_records.append(record)
        selected_rows.extend(record["selected"])

    selector = selector_identity()
    seqspec = executable_identity(seqspec_bin)
    stable_cases = [stable_case_row(value) for value in selected_rows]
    stable_manifest = {
        "schema_version": SCHEMA_VERSION,
        "freeze_id": str(cohort_manifest["freeze_id"]),
        "cohort_split": cohort_split,
        "cohort_sha256": file_sha256(cohort_path),
        "sampling_protocol_sha256": file_sha256(protocol_path),
        "selector": functional_selector_identity(selector),
        "seqspec": functional_executable_identity(seqspec),
        "cases": stable_cases,
    }
    selection_id = sha256_json(stable_manifest)[:16]
    cases = [
        {
            **public_case_row(value),
            "selection_id": selection_id,
            "case_id": case_id(value),
        }
        for value in selected_rows
    ]
    validation = validate_selection(
        cohort_manifest=cohort_manifest,
        protocol=protocol,
        cohort_rows=cohort_rows,
        configurations=configuration_records,
        cases=cases,
        cohort_split=cohort_split,
    )
    if not validation["valid"]:
        raise ValueError(
            "sampling case selection failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )

    output_root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_root.name}-", dir=output_root.parent
    ) as tmpdir:
        temporary_root = Path(tmpdir)
        tables_dir = temporary_root / "tables"
        manifests_dir = temporary_root / "manifests"
        validation_dir = temporary_root / "validation"
        tables_dir.mkdir()
        manifests_dir.mkdir()
        validation_dir.mkdir()
        cases_csv_path = tables_dir / "sampling_cases.csv"
        cases_json_path = manifests_dir / "sampling_cases.json"
        validation_path = validation_dir / "sampling_cases.json"
        manifest_path = manifests_dir / "sampling_case_selection.json"

        write_csv(cases_csv_path, cases, CASE_FIELDS)
        write_json(
            cases_json_path,
            {
                "schema_version": SCHEMA_VERSION,
                "selection_id": selection_id,
                "freeze_id": cohort_manifest["freeze_id"],
                "cohort_split": cohort_split,
                "cases": cases,
            },
        )
        write_json(validation_path, validation)
        final_manifest_path = output_root / manifest_path.relative_to(temporary_root)
        manifest = {
            **stable_manifest,
            "selection_id": selection_id,
            "created_at": utc_now(),
            "valid": True,
            "inputs": {
                "cohort_manifest": file_identity(cohort_manifest_path),
                "cohort": file_identity(cohort_path),
                "sampling_protocol": file_identity(protocol_path),
            },
            "tooling": {"selector": selector, "seqspec": seqspec},
            "counts": validation["counts"],
            "declared_transfer": validation["declared_transfer"],
            "outputs": {
                "cases_csv": final_file_identity(
                    output_root, temporary_root, cases_csv_path
                ),
                "cases_json": final_file_identity(
                    output_root, temporary_root, cases_json_path
                ),
                "validation": final_file_identity(
                    output_root, temporary_root, validation_path
                ),
            },
            "manifest_path": str(final_manifest_path),
        }
        write_json(manifest_path, manifest)
        os.replace(temporary_root, output_root)
    return manifest


def validate_inputs(
    cohort_manifest: dict[str, Any],
    protocol: dict[str, Any],
    seqspec_bin: Path,
    cohort_split: str,
) -> tuple[Path, list[dict[str, str]]]:
    if cohort_split not in COHORT_SPLITS:
        raise ValueError(f"unsupported cohort split: {cohort_split}")
    if cohort_manifest.get("frozen") is not True:
        raise ValueError("cohort manifest is not frozen")
    if not str(cohort_manifest.get("freeze_id", "")).strip():
        raise ValueError("cohort manifest has no freeze id")
    cohort_identity = cohort_manifest.get("outputs", {}).get("cohort", {})
    if not isinstance(cohort_identity, dict):
        raise ValueError("cohort manifest has no cohort output identity")
    cohort_path = Path(str(cohort_identity.get("path", ""))).resolve()
    if file_sha256(cohort_path) != str(cohort_identity.get("sha256", "")):
        raise ValueError("frozen cohort table hash changed")
    cohort_rows = read_csv(cohort_path)
    expected_split = int(cohort_manifest.get("counts", {}).get(cohort_split, -1))
    observed_split = sum(row.get("split") == cohort_split for row in cohort_rows)
    if expected_split <= 0 or observed_split != expected_split:
        raise ValueError(f"frozen cohort {cohort_split} count does not reconcile")

    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling protocol does not use schema 0.1.0")
    limits = cohort_limits(protocol, cohort_split)
    selection = protocol.get("case_selection", {})
    if int(limits.get("configurations", -1)) != observed_split:
        raise ValueError("sampling protocol configuration count does not match cohort")
    max_fastqs = int(selection.get("max_fastqs_per_configuration", 0))
    if max_fastqs not in {1, 2}:
        raise ValueError(
            "sampling protocol must select one or two FASTQs per configuration"
        )
    if int(limits.get("fastqs", -1)) != observed_split * max_fastqs:
        raise ValueError(
            "sampling protocol FASTQ cap does not match its selection rule"
        )
    if selection.get("primary_read_rule") != "largest_measurement_span":
        raise ValueError("unsupported primary read selection rule")
    if selection.get("secondary_read_rule") != (
        "maximum_additional_ontology_role_coverage"
    ):
        raise ValueError("unsupported secondary read selection rule")
    if not isinstance(selection.get("require_positive_declared_compressed_size"), bool):
        raise ValueError("sampling protocol size requirement must be boolean")
    if not seqspec_bin.is_file():
        raise ValueError(f"seqspec executable does not exist: {seqspec_bin}")
    return cohort_path, cohort_rows


def prepare_configuration(
    *,
    row: dict[str, str],
    seqspec_bin: Path,
    timeout_seconds: int,
    max_fastqs: int,
    require_positive_size: bool,
) -> dict[str, Any]:
    family_id = required_field(row, "final_family")
    configuration_accession = required_field(row, "configuration_accession")
    modalities = split_values(required_field(row, "modalities"))
    if len(modalities) != 1:
        raise ValueError(
            f"{configuration_accession}: selected cohort row must have one modality"
        )
    modality = modalities[0]
    spec_path = Path(required_field(row, "effective_spec_path")).resolve()
    declared_spec_sha256 = required_field(row, "effective_spec_sha256")
    if file_sha256(spec_path) != declared_spec_sha256:
        raise ValueError(f"{configuration_accession}: effective spec hash changed")

    accessions = split_values(required_field(row, "fastq_accessions"))
    urls = split_values(required_field(row, "fastq_urls"))
    if len(accessions) != len(urls) or len(accessions) != len(set(accessions)):
        raise ValueError(
            f"{configuration_accession}: FASTQ accessions and URLs do not align"
        )
    expected_accessions = set(
        split_values(required_field(row, "effective_expected_fastq_accessions"))
    )
    if set(accessions) != expected_accessions:
        raise ValueError(
            f"{configuration_accession}: frozen FASTQ mapping is not exact"
        )

    sequence_spec = run_seqspec_json(
        seqspec_bin,
        ["info", "-k", "sequence_spec", "-f", "json", str(spec_path)],
        timeout_seconds,
    )
    library_spec = run_seqspec_json(
        seqspec_bin,
        ["info", "-k", "library_spec", "-f", "json", str(spec_path)],
        timeout_seconds,
    )
    index_output = run_command(
        [
            str(seqspec_bin),
            "index",
            "--modality",
            modality,
            str(spec_path),
        ],
        timeout_seconds,
    ).stdout
    region_index_output = run_command(
        [
            str(seqspec_bin),
            "index",
            "--selector",
            "region",
            "--modality",
            modality,
            str(spec_path),
        ],
        timeout_seconds,
    ).stdout
    leaf_regions = extract_leaf_regions(library_spec, modality)
    label_terms = parse_region_labels(
        region_index_output, leaf_regions, configuration_accession
    )
    read_spans = parse_index(index_output, label_terms, configuration_accession)
    candidates = build_candidates(
        sequence_spec=sequence_spec,
        read_spans=read_spans,
        accessions=accessions,
        urls=urls,
        family_id=family_id,
        configuration_accession=configuration_accession,
        modality=modality,
        spec_path=spec_path,
        spec_sha256=declared_spec_sha256,
        require_positive_size=require_positive_size,
    )
    selected = select_candidates(candidates, max_fastqs)
    return {
        "family_id": family_id,
        "configuration_accession": configuration_accession,
        "modality": modality,
        "candidate_count": len(candidates),
        "candidate_role_classes": sorted(
            {role for value in candidates for role in value["role_class_values"]}
        ),
        "candidate_ontology_terms": sorted(
            {term for value in candidates for term in value["ontology_term_values"]}
        ),
        "selected": selected,
    }


def build_candidates(
    *,
    sequence_spec: Any,
    read_spans: dict[str, dict[str, Any]],
    accessions: list[str],
    urls: list[str],
    family_id: str,
    configuration_accession: str,
    modality: str,
    spec_path: Path,
    spec_sha256: str,
    require_positive_size: bool,
) -> list[dict[str, Any]]:
    if not isinstance(sequence_spec, list):
        raise ValueError(f"{configuration_accession}: sequence_spec is not a list")
    url_by_accession = dict(zip(accessions, urls, strict=True))
    candidates = []
    observed_accessions = set()
    for read in sequence_spec:
        if not isinstance(read, dict) or str(read.get("modality", "")) != modality:
            continue
        read_id = str(read.get("read_id", "")).strip()
        if not read_id:
            raise ValueError(f"{configuration_accession}: read has no read_id")
        profile = read_spans.get(read_id, empty_profile())
        files = read.get("files", [])
        if not isinstance(files, list):
            raise ValueError(f"{configuration_accession}: read files are not a list")
        for file_value in files:
            if not isinstance(file_value, dict):
                continue
            accession = sequence_file_accession(file_value)
            if accession not in url_by_accession:
                continue
            if accession in observed_accessions:
                raise ValueError(
                    f"{configuration_accession}: FASTQ {accession} maps to multiple reads"
                )
            size = int(file_value.get("filesize") or 0)
            if require_positive_size and size <= 0:
                raise ValueError(
                    f"{configuration_accession}: FASTQ {accession} has no positive "
                    "declared compressed size"
                )
            observed_accessions.add(accession)
            role_classes = sorted(
                role for role in ROLE_PREFIXES if profile[f"{role}_bases"] > 0
            )
            ontology_terms = sorted(profile["ontology_terms"])
            candidates.append(
                {
                    "family_id": family_id,
                    "configuration_accession": configuration_accession,
                    "modality": modality,
                    "selection_role": "",
                    "read_id": read_id,
                    "fastq_accession": accession,
                    "fastq_url": url_by_accession[accession],
                    "declared_compressed_bytes": size,
                    "indexed_bases": profile["indexed_bases"],
                    "measure_bases": profile["measure_bases"],
                    "partition_bases": profile["partition_bases"],
                    "technical_bases": profile["technical_bases"],
                    "role_classes": ";".join(role_classes),
                    "role_class_values": role_classes,
                    "ontology_terms": ";".join(ontology_terms),
                    "ontology_term_values": ontology_terms,
                    "spec_path": str(spec_path),
                    "spec_sha256": spec_sha256,
                }
            )
    missing = sorted(set(accessions) - observed_accessions)
    if missing:
        raise ValueError(
            f"{configuration_accession}: FASTQs do not map to modality {modality}: "
            + ", ".join(missing)
        )
    return candidates


def select_candidates(
    candidates: list[dict[str, Any]], max_fastqs: int
) -> list[dict[str, Any]]:
    if not candidates:
        raise ValueError("configuration has no FASTQ candidates")
    has_measurement = any(value["measure_bases"] > 0 for value in candidates)
    primary = min(
        candidates,
        key=(primary_measure_key if has_measurement else primary_fallback_key),
    )
    selected = [{**primary, "selection_role": "primary"}]
    if max_fastqs == 1:
        return selected

    remaining = [value for value in candidates if value is not primary]
    if not remaining:
        return selected
    primary_roles = set(primary["role_class_values"])
    primary_terms = selection_ontology_terms(primary)
    secondary = min(
        remaining,
        key=lambda value: secondary_key(value, primary_roles, primary_terms),
    )
    adds_role = bool(set(secondary["role_class_values"]) - primary_roles)
    adds_term = bool(selection_ontology_terms(secondary) - primary_terms)
    if adds_role or adds_term:
        selected.append({**secondary, "selection_role": "secondary"})
    return selected


def primary_measure_key(value: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -int(value["measure_bases"]),
        -len(
            [
                term
                for term in value["ontology_term_values"]
                if term.startswith(ROLE_PREFIXES["measure"])
            ]
        ),
        -int(value["indexed_bases"]),
        int(value["declared_compressed_bytes"]),
        str(value["fastq_accession"]),
    )


def primary_fallback_key(value: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -int(value["partition_bases"]),
        -int(value["technical_bases"]),
        -int(value["indexed_bases"]),
        int(value["declared_compressed_bytes"]),
        str(value["fastq_accession"]),
    )


def secondary_key(
    value: dict[str, Any], primary_roles: set[str], primary_terms: set[str]
) -> tuple[Any, ...]:
    new_roles = set(value["role_class_values"]) - primary_roles
    new_terms = selection_ontology_terms(value) - primary_terms
    return (
        -len(new_roles),
        -len(new_terms),
        -int(value["partition_bases"]),
        -int(value["technical_bases"]),
        -int(value["indexed_bases"]),
        int(value["declared_compressed_bytes"]),
        str(value["fastq_accession"]),
    )


def selection_ontology_terms(value: dict[str, Any]) -> set[str]:
    return {
        term
        for term in value["ontology_term_values"]
        if any(term.startswith(prefix) for prefix in ROLE_PREFIXES.values())
    }


def extract_leaf_regions(library_spec: Any, modality: str) -> list[dict[str, Any]]:
    if isinstance(library_spec, dict) and modality in library_spec:
        roots = library_spec[modality]
    else:
        roots = library_spec
    if isinstance(roots, dict):
        roots = [roots]
    if not isinstance(roots, list):
        raise ValueError(f"library_spec modality {modality} is not a region list")
    leaves = []
    stack = list(reversed(roots))
    while stack:
        region = stack.pop()
        if not isinstance(region, dict):
            continue
        children = region.get("regions", [])
        if isinstance(children, list) and children:
            stack.extend(reversed(children))
            continue
        region_id = str(region.get("region_id", "")).strip()
        name = str(region.get("name", "")).strip()
        if not region_id or not name:
            raise ValueError("library_spec leaf is missing region_id or name")
        values = region.get("region_type", [])
        if isinstance(values, str):
            values = [values]
        if not isinstance(values, list):
            raise ValueError(f"region {region_id} has invalid region_type")
        leaves.append(
            {
                "region_id": region_id,
                "name": name,
                "ontology_terms": {
                    str(value).strip() for value in values if str(value).strip()
                },
            }
        )
    return leaves


def parse_region_labels(
    output: str,
    leaf_regions: list[dict[str, Any]],
    configuration_accession: str,
) -> dict[tuple[str, str], set[str]]:
    rows = [line.split("\t") for line in output.splitlines() if line.strip()]
    if any(len(fields) != 5 for fields in rows):
        raise ValueError(f"{configuration_accession}: invalid region-level index")
    if len(rows) != len(leaf_regions):
        raise ValueError(
            f"{configuration_accession}: region index has {len(rows)} leaves, "
            f"library_spec has {len(leaf_regions)}"
        )
    label_terms = {}
    for fields, region in zip(rows, leaf_regions, strict=True):
        _, name, tool_label, _, _ = fields
        if name != region["name"]:
            raise ValueError(
                f"{configuration_accession}: region index order differs at "
                f"{region['region_id']}"
            )
        key = (name, tool_label)
        terms = region["ontology_terms"]
        if key in label_terms and label_terms[key] != terms:
            raise ValueError(
                f"{configuration_accession}: tab index label {name}/{tool_label} "
                "maps to multiple ontology term sets"
            )
        label_terms[key] = terms
    return label_terms


def parse_index(
    output: str,
    label_terms: dict[tuple[str, str], set[str]],
    configuration_accession: str,
) -> dict[str, dict[str, Any]]:
    profiles: dict[str, dict[str, Any]] = {}
    for line_number, line in enumerate(output.splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 5:
            raise ValueError(
                f"{configuration_accession}: invalid seqspec index line {line_number}"
            )
        read_id, name, tool_label, start_value, end_value = fields
        key = (name, tool_label)
        if key not in label_terms:
            raise ValueError(
                f"{configuration_accession}: indexed label {name}/{tool_label} is "
                "not in the region-level index"
            )
        try:
            start = int(start_value)
            end = int(end_value)
        except ValueError as error:
            raise ValueError(
                f"{configuration_accession}: invalid index coordinates on line "
                f"{line_number}"
            ) from error
        if start < 0 or end <= start:
            raise ValueError(
                f"{configuration_accession}: non-positive index span on line "
                f"{line_number}"
            )
        span = end - start
        profile = profiles.setdefault(read_id, empty_profile())
        profile["indexed_bases"] += span
        terms = label_terms[key]
        profile["ontology_terms"].update(terms)
        for role, prefix in ROLE_PREFIXES.items():
            if any(term.startswith(prefix) for term in terms):
                profile[f"{role}_bases"] += span
    return profiles


def empty_profile() -> dict[str, Any]:
    return {
        "indexed_bases": 0,
        "measure_bases": 0,
        "partition_bases": 0,
        "technical_bases": 0,
        "ontology_terms": set(),
    }


def validate_selection(
    *,
    cohort_manifest: dict[str, Any],
    protocol: dict[str, Any],
    cohort_rows: list[dict[str, str]],
    configurations: list[dict[str, Any]],
    cases: list[dict[str, Any]],
    cohort_split: str,
) -> dict[str, Any]:
    errors = []
    limits = cohort_limits(protocol, cohort_split)
    expected_configurations = int(limits["configurations"])
    max_fastqs = int(protocol["case_selection"]["max_fastqs_per_configuration"])
    expected_split = int(cohort_manifest["counts"][cohort_split])
    if len(configurations) != expected_configurations:
        errors.append(
            f"expected {expected_configurations} configurations, "
            f"observed {len(configurations)}"
        )
    if len(configurations) != expected_split:
        errors.append(
            f"selected configurations do not match frozen {cohort_split} count"
        )
    if sum(row.get("split") == cohort_split for row in cohort_rows) != len(
        configurations
    ):
        errors.append(f"not every {cohort_split} row produced a configuration")

    for configuration in configurations:
        selected = configuration["selected"]
        accession = configuration["configuration_accession"]
        if not 1 <= len(selected) <= max_fastqs:
            errors.append(f"{accession}: selected FASTQ count is outside the protocol")
            continue
        if selected[0]["selection_role"] != "primary":
            errors.append(f"{accession}: first selected FASTQ is not primary")
        candidate_has_measurement = any(
            role == "measure" for role in configuration["candidate_role_classes"]
        )
        if candidate_has_measurement and selected[0]["measure_bases"] <= 0:
            errors.append(
                f"{accession}: primary FASTQ omits available measurement span"
            )
        if len(selected) == 2:
            first_roles = set(selected[0]["role_class_values"])
            first_terms = selection_ontology_terms(selected[0])
            second_roles = set(selected[1]["role_class_values"])
            second_terms = selection_ontology_terms(selected[1])
            if not (second_roles - first_roles or second_terms - first_terms):
                errors.append(f"{accession}: secondary FASTQ adds no ontology coverage")

    case_ids = [value["case_id"] for value in cases]
    fastq_keys = [
        (value["configuration_accession"], value["fastq_accession"]) for value in cases
    ]
    if len(case_ids) != len(set(case_ids)):
        errors.append("case identifiers are not unique")
    if len(fastq_keys) != len(set(fastq_keys)):
        errors.append("selected configuration/FASTQ pairs are not unique")
    if len(cases) > int(limits["fastqs"]):
        errors.append("selected FASTQ count exceeds protocol cap")

    source_bytes = sum(int(value["declared_compressed_bytes"]) for value in cases)
    calibration = cohort_split == "calibration"
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": utc_now(),
        "cohort_split": cohort_split,
        "valid": not errors,
        "errors": errors,
        "counts": {
            "configurations": len(configurations),
            "fastq_cases": len(cases),
            "configuration_families": dict(
                sorted(Counter(value["family_id"] for value in configurations).items())
            ),
            "fastq_case_families": dict(
                sorted(Counter(value["family_id"] for value in cases).items())
            ),
        },
        "declared_transfer": {
            "selected_source_bytes": source_bytes,
            "minimum_source_traversals_per_fastq": 2 if calibration else 1,
            "minimum_remote_bytes": source_bytes * (2 if calibration else 1),
            "explanation": (
                "One traversal creates bounded samples and one traversal creates "
                "the complete-stream seqcheck reference."
                if calibration
                else "One traversal creates the frozen-policy bounded sample."
            ),
        },
    }


def cohort_limits(protocol: dict[str, Any], cohort_split: str) -> dict[str, Any]:
    limits = protocol.get("cohort_limits", {})
    value = limits.get(cohort_split) if isinstance(limits, dict) else None
    if not isinstance(value, dict):
        raise ValueError(f"sampling protocol has no {cohort_split} cohort limits")
    return value


def public_case_row(value: dict[str, Any]) -> dict[str, Any]:
    return {
        key: field_value
        for key, field_value in value.items()
        if key not in {"role_class_values", "ontology_term_values"}
    }


def stable_case_row(value: dict[str, Any]) -> dict[str, Any]:
    return {
        key: field_value
        for key, field_value in public_case_row(value).items()
        if key != "spec_path"
    }


def case_id(value: dict[str, Any]) -> str:
    return (f"{value['configuration_accession']}--{value['fastq_accession']}").lower()


def run_seqspec_json(seqspec_bin: Path, args: list[str], timeout_seconds: int) -> Any:
    completed = run_command([str(seqspec_bin), *args], timeout_seconds)
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise ValueError(
            f"seqspec returned invalid JSON for {' '.join(args)}"
        ) from error


def run_command(
    argv: list[str], timeout_seconds: int
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        argv,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip()
        raise ValueError(f"{' '.join(argv)} failed: {message}")
    return completed


def selector_identity() -> dict[str, Any]:
    path = Path(__file__).resolve()
    root = path.parents[1]
    status = git_output(
        root,
        "status",
        "--porcelain",
        "--untracked-files=normal",
        "--",
        str(path.relative_to(root)),
    )
    return {
        "version": SELECTOR_VERSION,
        "script_path": str(path),
        "script_sha256": file_sha256(path),
        "git_commit": git_output(root, "rev-parse", "HEAD"),
        "git_dirty": bool(status),
        "python": sys.version.split()[0],
    }


def executable_identity(path: Path) -> dict[str, Any]:
    completed = run_command([str(path), "--version"], 10)
    version = completed.stdout.strip() or completed.stderr.strip()
    if not version:
        raise ValueError(f"seqspec executable returned no version: {path}")
    return {
        "path": str(path),
        "sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
        "version": version,
    }


def functional_selector_identity(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "version": value["version"],
        "script_sha256": value["script_sha256"],
        "python": value["python"],
    }


def functional_executable_identity(value: dict[str, Any]) -> dict[str, Any]:
    return {"sha256": value["sha256"], "version": value["version"]}


def git_output(root: Path, *args: str) -> str:
    if not (root / ".git").exists():
        return ""
    completed = subprocess.run(
        ["git", "-C", str(root), *args],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def required_field(row: dict[str, str], field: str) -> str:
    value = row.get(field, "").strip()
    if not value:
        raise ValueError(f"cohort row is missing {field}")
    return value


def split_values(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def sequence_file_accession(value: dict[str, Any]) -> str:
    url = str(value.get("url", "")).strip()
    pieces = [piece for piece in urllib.parse.urlparse(url).path.split("/") if piece]
    for index, piece in enumerate(pieces[:-1]):
        if piece == "sequence-files":
            return normalize_fastq_accession(pieces[index + 1])
    for field in ("file_id", "filename"):
        candidate = str(value.get(field, "")).strip()
        if candidate:
            return normalize_fastq_accession(Path(candidate).name)
    return ""


def normalize_fastq_accession(value: str) -> str:
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if value.lower().endswith(suffix):
            return value[: -len(suffix)]
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV header is missing: {path}")
        return [dict(row) for row in reader]


def write_csv(path: Path, rows: list[dict[str, Any]], fields: tuple[str, ...]) -> None:
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


def final_file_identity(
    output_root: Path, temporary_root: Path, temporary_path: Path
) -> dict[str, Any]:
    relative = temporary_path.relative_to(temporary_root)
    return {
        "path": str(output_root / relative),
        "sha256": file_sha256(temporary_path),
        "size_bytes": temporary_path.stat().st_size,
    }


def file_identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
