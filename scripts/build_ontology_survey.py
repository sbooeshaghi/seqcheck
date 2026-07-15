#!/usr/bin/env python3
"""Build the paired seqspec region ontology survey and weighted review sample."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
BUILDER_VERSION = "0.1.0"
REGION_FIELDS = (
    "survey_id",
    "region_key",
    "configuration_accession",
    "family_ids",
    "family_labels",
    "assay_term",
    "preferred_assay_titles",
    "assay_id",
    "assay_name",
    "assay_description",
    "raw_seqspec_version",
    "normalized_seqspec_version",
    "raw_spec_path",
    "raw_spec_sha256",
    "normalized_spec_path",
    "normalized_spec_sha256",
    "modality",
    "modality_read_ids",
    "modality_primer_ids",
    "region_id",
    "region_name",
    "sequence_type",
    "min_len",
    "max_len",
    "depth",
    "is_leaf",
    "parent_region_id",
    "parent_region_name",
    "path_region_ids",
    "path_region_names",
    "original_region_type_json",
    "original_region_type_labels",
    "ontology_terms",
    "mapping_source",
    "unknown_mapping",
    "multi_term_mapping",
    "has_onlist",
)
SAMPLE_FIELDS = (
    *REGION_FIELDS,
    "review_stratum",
    "stratum_population",
    "stratum_sample_size",
    "inclusion_probability",
    "sampling_weight",
    "selection_rank",
)
LABEL_FIELDS = (
    "survey_id",
    "original_region_type",
    "regions",
    "configurations",
    "ontology_terms",
    "unknown_regions",
    "multi_term_regions",
)
TERM_FIELDS = (
    "survey_id",
    "ontology_term",
    "regions",
    "configurations",
    "original_region_types",
    "sequence_types",
)
SPEC_FIELDS = (
    "survey_id",
    "configuration_accession",
    "regions",
    "unknown_regions",
    "multi_term_regions",
    "original_region_types",
    "ontology_terms",
)
QUERY_FIELDS = (
    "survey_id",
    "configuration_accession",
    "modality",
    "selector",
    "selector_kind",
    "query_terms",
    "raw_region_keys",
    "normalized_region_keys",
    "raw_region_count",
    "normalized_region_count",
    "added_region_keys",
    "removed_region_keys",
    "change_class",
    "unexpected_change",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the paired region ontology mapping survey."
    )
    parser.add_argument("--cohort-manifest", required=True, type=Path)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--seqspec-bin", required=True, type=Path)
    parser.add_argument("--yq-bin", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=30)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = build_ontology_survey(
            cohort_manifest_path=args.cohort_manifest.resolve(),
            registry_path=args.registry.resolve(),
            protocol_path=args.protocol.resolve(),
            seqspec_bin=args.seqspec_bin.resolve(),
            yq_bin=args.yq_bin.resolve(),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
        )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"build_ontology_survey: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def build_ontology_survey(
    *,
    cohort_manifest_path: Path,
    registry_path: Path,
    protocol_path: Path,
    seqspec_bin: Path,
    yq_bin: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    protocol = runtime.load_json(protocol_path)
    validate_protocol(protocol)
    cohort, proposal_path, specifications = load_population(
        cohort_manifest_path, protocol
    )
    registry = load_yaml_json(registry_path, yq_bin, timeout_seconds)
    legacy_map, registry_terms = validate_registry(registry, protocol)
    seqspec_identity = runtime.executable_identity(
        seqspec_bin, timeout_seconds=min(timeout_seconds, 30)
    )
    runtime_parity = validate_runtime_registry_parity(
        seqspec_bin=seqspec_bin,
        yq_bin=yq_bin,
        legacy_map=legacy_map,
        timeout_seconds=timeout_seconds,
    )
    yq_identity = runtime.executable_identity(
        yq_bin, timeout_seconds=min(timeout_seconds, 30)
    )
    script = runtime.script_identity(Path(__file__).resolve(), version=BUILDER_VERSION)
    runtime_identity = runtime.script_identity(
        Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
    )
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "cohort_selection_id": cohort["selection_id"],
        "cohort_manifest_sha256": runtime.file_sha256(cohort_manifest_path),
        "proposal_table_sha256": runtime.file_sha256(proposal_path),
        "registry_sha256": runtime.file_sha256(registry_path),
        "protocol_sha256": runtime.file_sha256(protocol_path),
        "specifications": [
            {
                "configuration_accession": value["configuration_accession"],
                "raw_spec_sha256": value["raw_spec_sha256"],
                "normalized_spec_sha256": value["normalized_spec_sha256"],
            }
            for value in specifications
        ],
        "yq": runtime.functional_executable_identity(yq_identity),
        "seqspec": runtime.functional_executable_identity(seqspec_identity),
        "registry_runtime_parity": runtime_parity,
        "builder": runtime.functional_script_identity(script),
        "runtime": runtime.functional_script_identity(runtime_identity),
    }
    survey_id = runtime.sha256_json(stable_identity)[:16]
    rows = []
    query_rows = []
    for specification in specifications:
        spec_rows, spec_queries = survey_specification(
            survey_id=survey_id,
            specification=specification,
            legacy_map=legacy_map,
            registry_terms=registry_terms,
            protocol=protocol,
            yq_bin=yq_bin,
            timeout_seconds=timeout_seconds,
        )
        rows.extend(spec_rows)
        query_rows.extend(spec_queries)
    rows.sort(
        key=lambda value: (
            value["configuration_accession"],
            value["modality"],
            value["path_region_ids"],
        )
    )
    query_rows.sort(
        key=lambda value: (
            value["configuration_accession"],
            value["modality"],
            value["selector"],
        )
    )
    if len({value["region_key"] for value in rows}) != len(rows):
        raise ValueError("region survey keys are not unique")
    sample = select_review_sample(rows, protocol)
    label_summary = summarize_labels(survey_id, rows)
    term_summary = summarize_terms(survey_id, rows)
    spec_summary = summarize_specs(survey_id, rows)
    validation = validate_outputs(
        survey_id=survey_id,
        specifications=specifications,
        rows=rows,
        query_rows=query_rows,
        sample=sample,
        protocol=protocol,
        runtime_parity=runtime_parity,
    )
    if not validation["valid"]:
        raise ValueError(
            "ontology survey failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )

    tables = output_root / "tables"
    paths = {
        "regions": tables / "ontology_regions.csv",
        "review_sample": tables / "review_sample.csv",
        "original_labels": tables / "original_label_summary.csv",
        "ontology_terms": tables / "ontology_term_summary.csv",
        "specifications": tables / "specification_summary.csv",
        "query_parity": tables / "query_parity.csv",
    }
    runtime.write_csv(paths["regions"], rows, list(REGION_FIELDS))
    runtime.write_csv(paths["review_sample"], sample, list(SAMPLE_FIELDS))
    runtime.write_csv(paths["original_labels"], label_summary, list(LABEL_FIELDS))
    runtime.write_csv(paths["ontology_terms"], term_summary, list(TERM_FIELDS))
    runtime.write_csv(paths["specifications"], spec_summary, list(SPEC_FIELDS))
    runtime.write_csv(paths["query_parity"], query_rows, list(QUERY_FIELDS))
    validation_path = output_root / "validation" / "ontology_survey.json"
    runtime.write_json(validation_path, validation)
    manifest_path = output_root / "manifests" / "ontology_survey.json"
    manifest = {
        **stable_identity,
        "survey_id": survey_id,
        "created_at": utc_now(),
        "valid": True,
        "population_scope": "successfully normalized specification pairs in cohort proposal snapshot",
        "inputs": {
            "cohort_manifest": runtime.file_identity(cohort_manifest_path),
            "proposal_table": runtime.file_identity(proposal_path),
            "registry": runtime.file_identity(registry_path),
            "protocol": runtime.file_identity(protocol_path),
        },
        "tools": {
            "builder": script,
            "runtime": runtime_identity,
            "seqspec": seqspec_identity,
            "yq": yq_identity,
        },
        "counts": validation["counts"],
        "outputs": {
            **{name: runtime.file_identity(path) for name, path in paths.items()},
            "validation": runtime.file_identity(validation_path),
        },
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest


def load_population(
    cohort_manifest_path: Path, protocol: dict[str, Any]
) -> tuple[dict[str, Any], Path, list[dict[str, Any]]]:
    cohort = runtime.load_json(cohort_manifest_path)
    if cohort.get("cohort_candidate_schema_version") != "0.2.0":
        raise ValueError("cohort candidate schema is unsupported")
    selection_id = str(cohort.get("selection_id", "")).strip()
    if not selection_id:
        raise ValueError("cohort selection id is missing")
    proposal_value = str(cohort.get("outputs", {}).get("cohort_proposals", ""))
    proposal_path = Path(proposal_value).resolve()
    proposals = runtime.read_csv(proposal_path)
    required_status = protocol["population"]["required_hydration_status"]
    grouped = defaultdict(list)
    for row in proposals:
        if row.get("hydration_status") == required_status:
            grouped[row["configuration_accession"]].append(row)
    if not grouped:
        raise ValueError("cohort proposal table has no normalized specifications")
    run_root = cohort_manifest_path.resolve().parents[1]
    result = []
    for accession, values in sorted(grouped.items()):
        if any(value.get("selection_id") != selection_id for value in values):
            raise ValueError(f"{accession}: proposal selection identifiers differ")
        invariant_fields = (
            "raw_spec_sha256",
            "raw_seqspec_version",
            "normalized_spec_sha256",
            "normalized_seqspec_version",
            "normalized_spec_path",
            "assay_term",
            "preferred_assay_titles",
        )
        for field in invariant_fields:
            if len({value.get(field, "") for value in values}) != 1:
                raise ValueError(f"{accession}: proposal field differs: {field}")
        row = values[0]
        normalized_path = Path(row["normalized_spec_path"]).resolve()
        raw_path = run_root / "specs" / "raw" / f"{accession}.yaml"
        verify_hash(raw_path, row["raw_spec_sha256"], f"{accession} raw spec")
        verify_hash(
            normalized_path,
            row["normalized_spec_sha256"],
            f"{accession} normalized spec",
        )
        required_version = protocol["population"]["required_normalized_seqspec_version"]
        if row["normalized_seqspec_version"] != required_version:
            raise ValueError(f"{accession}: normalized seqspec version differs")
        result.append(
            {
                **row,
                "family_ids": ";".join(
                    sorted({value["family_id"] for value in values})
                ),
                "family_labels": ";".join(
                    sorted({value["family_label"] for value in values})
                ),
                "raw_spec_path": str(raw_path.resolve()),
                "normalized_spec_path": str(normalized_path),
            }
        )
    return cohort, proposal_path, result


def survey_specification(
    *,
    survey_id: str,
    specification: dict[str, Any],
    legacy_map: dict[str, list[str]],
    registry_terms: set[str],
    protocol: dict[str, Any],
    yq_bin: Path,
    timeout_seconds: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    raw = load_yaml_json(Path(specification["raw_spec_path"]), yq_bin, timeout_seconds)
    normalized = load_yaml_json(
        Path(specification["normalized_spec_path"]), yq_bin, timeout_seconds
    )
    raw_regions = flatten_regions(raw, protocol)
    normalized_regions = flatten_regions(normalized, protocol)
    if set(raw_regions) != set(normalized_regions):
        raise ValueError(
            f"{specification['configuration_accession']}: raw and normalized region paths differ"
        )
    compare_read_geometry(specification["configuration_accession"], raw, normalized)
    read_context = modality_read_context(raw)
    result = []
    for path, raw_region in sorted(raw_regions.items()):
        normalized_region = normalized_regions[path]
        for field in (
            "region_id",
            "name",
            "sequence_type",
            "sequence",
            "min_len",
            "max_len",
        ):
            if raw_region.get(field) != normalized_region.get(field):
                raise ValueError(
                    f"{specification['configuration_accession']} {path}: {field} changed"
                )
        if normalize_resource(raw_region.get("onlist")) != normalize_resource(
            normalized_region.get("onlist")
        ):
            raise ValueError(
                f"{specification['configuration_accession']} {path}: onlist changed"
            )
        labels = region_type_values(raw_region.get("region_type"))
        observed = region_type_values(normalized_region.get("region_type"))
        expected, source = map_region_types(labels, legacy_map)
        if (
            protocol["mapping"]["require_exact_upgrade_mapping"]
            and observed != expected
        ):
            raise ValueError(
                f"{specification['configuration_accession']} {path}: ontology mapping differs; "
                f"expected {expected}, observed {observed}"
            )
        undeclared = sorted(set(observed) - registry_terms)
        if undeclared:
            raise ValueError(
                f"{specification['configuration_accession']} {path}: undeclared terms {undeclared}"
            )
        modality = path[0]
        context = read_context.get(modality, {"read_ids": [], "primer_ids": []})
        parent = raw_region["_parent"]
        region_key = make_region_key(
            specification["configuration_accession"], modality, path[1:]
        )
        result.append(
            {
                "survey_id": survey_id,
                "region_key": region_key,
                "configuration_accession": specification["configuration_accession"],
                "family_ids": specification["family_ids"],
                "family_labels": specification["family_labels"],
                "assay_term": specification["assay_term"],
                "preferred_assay_titles": specification["preferred_assay_titles"],
                "assay_id": str(raw.get("assay_id", "")),
                "assay_name": str(raw.get("name", "")),
                "assay_description": str(raw.get("description", "")),
                "raw_seqspec_version": specification["raw_seqspec_version"],
                "normalized_seqspec_version": specification[
                    "normalized_seqspec_version"
                ],
                "raw_spec_path": specification["raw_spec_path"],
                "raw_spec_sha256": specification["raw_spec_sha256"],
                "normalized_spec_path": specification["normalized_spec_path"],
                "normalized_spec_sha256": specification["normalized_spec_sha256"],
                "modality": modality,
                "modality_read_ids": ";".join(context["read_ids"]),
                "modality_primer_ids": ";".join(context["primer_ids"]),
                "region_id": raw_region["region_id"],
                "region_name": str(raw_region.get("name", "")),
                "sequence_type": str(raw_region.get("sequence_type", "")),
                "min_len": raw_region.get("min_len"),
                "max_len": raw_region.get("max_len"),
                "depth": raw_region["_depth"],
                "is_leaf": not raw_region["_children"],
                "parent_region_id": parent.get("region_id", "") if parent else "",
                "parent_region_name": parent.get("name", "") if parent else "",
                "path_region_ids": ";".join(path[1:]),
                "path_region_names": ";".join(raw_region["_path_names"]),
                "original_region_type_json": runtime.canonical_json(
                    raw_region.get("region_type")
                ),
                "original_region_type_labels": ";".join(labels),
                "ontology_terms": ";".join(observed),
                "mapping_source": source,
                "unknown_mapping": protocol["mapping"]["unknown_term"] in observed,
                "multi_term_mapping": len(observed) > 1,
                "has_onlist": raw_region.get("onlist") is not None,
            }
        )
    queries = build_query_parity(
        survey_id=survey_id,
        configuration_accession=specification["configuration_accession"],
        raw_regions=raw_regions,
        normalized_regions=normalized_regions,
        region_rows=result,
        legacy_map=legacy_map,
        unknown_term=protocol["mapping"]["unknown_term"],
    )
    return result, queries


def flatten_regions(
    spec: dict[str, Any], protocol: dict[str, Any]
) -> dict[tuple[str, ...], dict[str, Any]]:
    roots = spec.get("library_spec")
    if not isinstance(roots, list) or not roots:
        raise ValueError("specification library_spec is empty")
    result = {}

    def visit(
        region: dict[str, Any],
        modality: str,
        path: tuple[str, ...],
        names: tuple[str, ...],
        parent: dict[str, Any] | None,
        depth: int,
    ) -> None:
        region_id = str(region.get("region_id", "")).strip()
        if not region_id:
            raise ValueError("region id is missing")
        current = (*path, region_id)
        children = region.get("regions") or []
        if not isinstance(children, list):
            raise ValueError(f"{region_id}: child regions are invalid")
        key = (modality, *current)
        if key in result:
            raise ValueError(f"duplicate region path: {key}")
        result[key] = {
            **region,
            "_parent": parent,
            "_depth": depth,
            "_children": children,
            "_path_names": (*names, str(region.get("name", ""))),
        }
        for child in children:
            if not isinstance(child, dict):
                raise ValueError(f"{region_id}: child region is invalid")
            visit(
                child,
                modality,
                current,
                result[key]["_path_names"],
                region,
                depth + 1,
            )

    include_containers = protocol["population"]["include_container_regions"]
    for root in roots:
        if not isinstance(root, dict):
            raise ValueError("library root is invalid")
        modality = str(root.get("region_id", "")).strip()
        visit(root, modality, (), (), None, 0)
    if not include_containers:
        result = {key: value for key, value in result.items() if not value["_children"]}
    return result


def modality_read_context(spec: dict[str, Any]) -> dict[str, dict[str, list[str]]]:
    grouped = defaultdict(lambda: {"read_ids": [], "primer_ids": []})
    for read in spec.get("sequence_spec") or []:
        modality = str(read.get("modality", ""))
        grouped[modality]["read_ids"].append(str(read.get("read_id", "")))
        grouped[modality]["primer_ids"].append(str(read.get("primer_id", "")))
    return {
        key: {
            "read_ids": sorted(set(value["read_ids"])),
            "primer_ids": sorted(set(value["primer_ids"])),
        }
        for key, value in grouped.items()
    }


def compare_read_geometry(
    configuration_accession: str,
    raw: dict[str, Any],
    normalized: dict[str, Any],
) -> None:
    raw_reads = functional_read_geometry(raw)
    normalized_reads = functional_read_geometry(normalized)
    if raw_reads != normalized_reads:
        raise ValueError(
            f"{configuration_accession}: read or file query geometry changed"
        )


def functional_read_geometry(spec: dict[str, Any]) -> list[dict[str, Any]]:
    reads = spec.get("sequence_spec") or []
    if not isinstance(reads, list):
        raise ValueError("specification sequence_spec is invalid")
    result = []
    for read in reads:
        if not isinstance(read, dict):
            raise ValueError("specification read is invalid")
        read_id = str(read.get("read_id", ""))
        files = read.get("files") or []
        if not files:
            files = [
                {
                    "file_id": read_id,
                    "filename": read_id,
                    "filetype": "",
                    "filesize": 0,
                    "url": "",
                    "urltype": "",
                    "md5": "",
                }
            ]
        result.append(
            {
                field: read.get(field)
                for field in (
                    "read_id",
                    "name",
                    "modality",
                    "primer_id",
                    "min_len",
                    "max_len",
                    "strand",
                )
            }
            | {"files": [normalize_resource(value) for value in files]}
        )
    return result


def normalize_resource(value: Any) -> dict[str, Any] | None:
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError("file or onlist resource is invalid")
    defaults = {
        "file_id": "",
        "filename": "",
        "filetype": "",
        "filesize": 0,
        "url": "",
        "urltype": "",
        "md5": "",
    }
    if "file_id" not in value and "location" in value:
        filename = value.get("filename", "")
        return {
            **defaults,
            "file_id": filename,
            "filename": filename,
            "md5": value.get("md5", ""),
        }
    return {field: value.get(field, default) for field, default in defaults.items()}


def build_query_parity(
    *,
    survey_id: str,
    configuration_accession: str,
    raw_regions: dict[tuple[str, ...], dict[str, Any]],
    normalized_regions: dict[tuple[str, ...], dict[str, Any]],
    region_rows: list[dict[str, Any]],
    legacy_map: dict[str, list[str]],
    unknown_term: str,
) -> list[dict[str, Any]]:
    row_keys = {value["region_key"] for value in region_rows}
    paths_by_modality = defaultdict(list)
    for path in raw_regions:
        paths_by_modality[path[0]].append(path)
    result = []
    for modality, paths in sorted(paths_by_modality.items()):
        selectors = {
            label
            for path in paths
            for label in region_type_values(raw_regions[path].get("region_type"))
        } | {
            term
            for path in paths
            for term in region_type_values(normalized_regions[path].get("region_type"))
        }
        for selector in sorted(selectors):
            query_terms = runtime_region_terms([selector], legacy_map)
            raw_keys = query_region_keys(
                configuration_accession,
                paths,
                raw_regions,
                selector,
                query_terms,
                legacy_map,
            )
            normalized_keys = query_region_keys(
                configuration_accession,
                paths,
                normalized_regions,
                selector,
                query_terms,
                legacy_map,
            )
            if not raw_keys | normalized_keys <= row_keys:
                raise ValueError("query parity refers to an unknown region key")
            added = sorted(normalized_keys - raw_keys)
            removed = sorted(raw_keys - normalized_keys)
            if not added and not removed:
                change_class = "stable"
            elif (
                unknown_term in query_terms
                or selector not in legacy_map
                and not selector.startswith("RGN:")
            ):
                change_class = "expected_unknown_collapse"
            else:
                change_class = "unexpected"
            result.append(
                {
                    "survey_id": survey_id,
                    "configuration_accession": configuration_accession,
                    "modality": modality,
                    "selector": selector,
                    "selector_kind": (
                        "ontology_term"
                        if selector.startswith("RGN:")
                        else "legacy_label"
                    ),
                    "query_terms": ";".join(query_terms),
                    "raw_region_keys": ";".join(sorted(raw_keys)),
                    "normalized_region_keys": ";".join(sorted(normalized_keys)),
                    "raw_region_count": len(raw_keys),
                    "normalized_region_count": len(normalized_keys),
                    "added_region_keys": ";".join(added),
                    "removed_region_keys": ";".join(removed),
                    "change_class": change_class,
                    "unexpected_change": change_class == "unexpected",
                }
            )
    return result


def query_region_keys(
    configuration_accession: str,
    paths: list[tuple[str, ...]],
    regions: dict[tuple[str, ...], dict[str, Any]],
    selector: str,
    query_terms: list[str],
    legacy_map: dict[str, list[str]],
) -> set[str]:
    result = set()
    query_set = set(query_terms)
    for path in paths:
        values = region_type_values(regions[path].get("region_type"))
        terms = set(runtime_region_terms(values, legacy_map))
        if selector in values or terms & query_set:
            result.add(make_region_key(configuration_accession, path[0], path[1:]))
    return result


def runtime_region_terms(
    values: list[str], legacy_map: dict[str, list[str]]
) -> list[str]:
    terms = []
    for value in values:
        mapped = legacy_map.get(value, [value])
        for term in mapped:
            if term not in terms:
                terms.append(term)
    return terms


def make_region_key(
    configuration_accession: str, modality: str, path: tuple[str, ...]
) -> str:
    return runtime.sha256_json(
        {
            "configuration_accession": configuration_accession,
            "modality": modality,
            "path": list(path),
        }
    )[:16]


def select_review_sample(
    rows: list[dict[str, Any]], protocol: dict[str, Any]
) -> list[dict[str, Any]]:
    review = protocol["manual_review"]
    label_counts = Counter(
        label
        for row in rows
        for label in split_values(row["original_region_type_labels"])
    )
    grouped = defaultdict(list)
    for row in rows:
        if row["unknown_mapping"]:
            stratum = "unknown"
        elif row["multi_term_mapping"]:
            stratum = "multi_term"
        elif any(
            label_counts[label] <= review["rare_label_max_population_count"]
            for label in split_values(row["original_region_type_labels"])
        ):
            stratum = "rare_label"
        else:
            stratum = "common"
        grouped[stratum].append(row)
    sample_size = min(review["sample_size"], len(rows))
    allocations = {
        stratum: min(review["target_allocations"][stratum], len(grouped[stratum]))
        for stratum in review["stratum_precedence"]
    }
    remaining = sample_size - sum(allocations.values())
    while remaining:
        progressed = False
        for stratum in review["stratum_precedence"]:
            if allocations[stratum] < len(grouped[stratum]):
                allocations[stratum] += 1
                remaining -= 1
                progressed = True
                if not remaining:
                    break
        if not progressed:
            raise ValueError("review allocation cannot fill the requested sample")
    seed = review["selection_seed"]
    result = []
    for stratum in review["stratum_precedence"]:
        population = sorted(
            grouped[stratum],
            key=lambda value: hashlib.sha256(
                f"{seed}:{value['region_key']}".encode("utf-8")
            ).hexdigest(),
        )
        selected = population[: allocations[stratum]]
        probability = len(selected) / len(population) if population else 0
        for rank, row in enumerate(selected, 1):
            result.append(
                {
                    **row,
                    "review_stratum": stratum,
                    "stratum_population": len(population),
                    "stratum_sample_size": len(selected),
                    "inclusion_probability": probability,
                    "sampling_weight": 1 / probability,
                    "selection_rank": rank,
                }
            )
    return sorted(result, key=lambda value: value["region_key"])


def summarize_labels(
    survey_id: str, rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    grouped = defaultdict(list)
    for row in rows:
        for label in split_values(row["original_region_type_labels"]):
            grouped[label].append(row)
    return [
        {
            "survey_id": survey_id,
            "original_region_type": label,
            "regions": len(values),
            "configurations": len(
                {value["configuration_accession"] for value in values}
            ),
            "ontology_terms": ";".join(
                sorted(
                    {
                        term
                        for value in values
                        for term in split_values(value["ontology_terms"])
                    }
                )
            ),
            "unknown_regions": sum(value["unknown_mapping"] for value in values),
            "multi_term_regions": sum(value["multi_term_mapping"] for value in values),
        }
        for label, values in sorted(grouped.items())
    ]


def summarize_terms(survey_id: str, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped = defaultdict(list)
    for row in rows:
        for term in split_values(row["ontology_terms"]):
            grouped[term].append(row)
    return [
        {
            "survey_id": survey_id,
            "ontology_term": term,
            "regions": len(values),
            "configurations": len(
                {value["configuration_accession"] for value in values}
            ),
            "original_region_types": ";".join(
                sorted(
                    {
                        label
                        for value in values
                        for label in split_values(value["original_region_type_labels"])
                    }
                )
            ),
            "sequence_types": ";".join(
                sorted({value["sequence_type"] for value in values})
            ),
        }
        for term, values in sorted(grouped.items())
    ]


def summarize_specs(survey_id: str, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["configuration_accession"]].append(row)
    return [
        {
            "survey_id": survey_id,
            "configuration_accession": accession,
            "regions": len(values),
            "unknown_regions": sum(value["unknown_mapping"] for value in values),
            "multi_term_regions": sum(value["multi_term_mapping"] for value in values),
            "original_region_types": len(
                {
                    label
                    for value in values
                    for label in split_values(value["original_region_type_labels"])
                }
            ),
            "ontology_terms": len(
                {
                    term
                    for value in values
                    for term in split_values(value["ontology_terms"])
                }
            ),
        }
        for accession, values in sorted(grouped.items())
    ]


def validate_protocol(protocol: dict[str, Any]) -> None:
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("ontology mapping protocol schema is unsupported")
    if protocol.get("experiment_id") != "region-ontology-mapping-survey":
        raise ValueError("ontology mapping experiment id is invalid")
    population = protocol.get("population", {})
    if population.get("required_hydration_status") != "normalized":
        raise ValueError("ontology population hydration status is invalid")
    if population.get("required_normalized_seqspec_version") != "0.5.0":
        raise ValueError("ontology population seqspec version is invalid")
    if population.get("pair_regions_by") != "modality_and_region_id_path":
        raise ValueError("region pairing rule is invalid")
    if population.get("include_container_regions") is not True:
        raise ValueError("mapping survey must include container regions")
    mapping = protocol.get("mapping", {})
    if mapping.get("unknown_term") != "RGN:unknown:unclassified":
        raise ValueError("unknown ontology term is invalid")
    if (
        mapping.get("require_registry_runtime_parity") is not True
        or mapping.get("require_exact_upgrade_mapping") is not True
    ):
        raise ValueError("mapping parity gates must be enabled")
    if mapping.get("query_semantics") != "seqspec_region_type_matches":
        raise ValueError("mapping query semantics are invalid")
    if mapping.get("classify_unknown_collapse_as_expected") is not True:
        raise ValueError("unknown mapping query policy must be explicit")
    review = protocol.get("manual_review", {})
    positive_int(review.get("sample_size"), "review sample size")
    positive_int(review.get("rare_label_max_population_count"), "rare label count")
    seed = review.get("selection_seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("review selection seed is invalid")
    precedence = review.get("stratum_precedence")
    expected = ["unknown", "multi_term", "rare_label", "common"]
    if precedence != expected:
        raise ValueError("review stratum precedence is invalid")
    allocations = review.get("target_allocations", {})
    if set(allocations) != set(expected) or any(
        isinstance(value, bool) or not isinstance(value, int) or value < 0
        for value in allocations.values()
    ):
        raise ValueError("review target allocations are invalid")
    if sum(allocations.values()) != review["sample_size"]:
        raise ValueError("review target allocations do not sum to sample size")
    if review.get("allocation_shortfall") != "round_robin_stratum_precedence":
        raise ValueError("review allocation shortfall policy is invalid")
    if review.get("reviewers") != 2 or review.get("blind_registry_mapping") is not True:
        raise ValueError("manual review must use two blinded reviewers")
    targets = protocol.get("scientific_targets", {})
    if targets.get("unexpected_query_changes_max") != 0:
        raise ValueError("unexpected query change target must be zero")
    for key in ("weighted_mapping_precision_min", "cohen_kappa_min"):
        value = targets.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"scientific target is invalid: {key}")
        if not 0 <= value <= 1:
            raise ValueError(f"scientific target is outside [0, 1]: {key}")


def validate_registry(
    registry: dict[str, Any], protocol: dict[str, Any]
) -> tuple[dict[str, list[str]], set[str]]:
    if registry.get("ontology_id") != "seqspec-region-ontology":
        raise ValueError("region ontology registry id is invalid")
    legacy = registry.get("legacy_region_types")
    terms = registry.get("terms")
    if not isinstance(legacy, dict) or not legacy or not isinstance(terms, dict):
        raise ValueError("region ontology registry mappings are invalid")
    term_ids = set(terms)
    unknown = protocol["mapping"]["unknown_term"]
    if unknown not in term_ids:
        raise ValueError("registry does not declare the unknown term")
    for label, mapped in legacy.items():
        if (
            not isinstance(label, str)
            or not label
            or not isinstance(mapped, list)
            or not mapped
            or len(mapped) != len(set(mapped))
            or any(value not in term_ids for value in mapped)
        ):
            raise ValueError(f"registry legacy mapping is invalid: {label}")
    return legacy, term_ids


def validate_runtime_registry_parity(
    *,
    seqspec_bin: Path,
    yq_bin: Path,
    legacy_map: dict[str, list[str]],
    timeout_seconds: int,
) -> dict[str, Any]:
    labels = sorted(legacy_map)
    children = [
        {
            "region_id": f"region_{index:03d}",
            "region_type": label,
            "name": label,
            "sequence_type": "fixed",
            "sequence": "A",
            "min_len": 1,
            "max_len": 1,
            "onlist": None,
            "regions": [],
        }
        for index, label in enumerate(labels)
    ]
    probe = {
        "seqspec_version": "0.4.0",
        "assay_id": "ontology-runtime-parity",
        "name": "Ontology runtime parity probe",
        "doi": "",
        "date": "",
        "description": "Generated registry parity probe",
        "modalities": ["probe"],
        "lib_struct": "",
        "sequence_protocol": None,
        "sequence_kit": None,
        "library_protocol": None,
        "library_kit": None,
        "sequence_spec": [],
        "library_spec": [
            {
                "region_id": "probe",
                "region_type": "RGN:unknown:unclassified",
                "name": "probe",
                "sequence_type": "joined",
                "sequence": "A" * len(children),
                "min_len": len(children),
                "max_len": len(children),
                "onlist": None,
                "regions": children,
            }
        ],
    }
    with tempfile.TemporaryDirectory(prefix="seqcheck-ontology-parity-") as tmpdir:
        source = Path(tmpdir) / "probe.yaml"
        upgraded_path = Path(tmpdir) / "upgraded.yaml"
        runtime.write_json(source, probe)
        completed = subprocess.run(
            [str(seqspec_bin), "upgrade", "-o", str(upgraded_path), str(source)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
        if completed.returncode != 0 or not upgraded_path.is_file():
            message = completed.stderr.strip() or completed.stdout.strip()
            raise ValueError(f"seqspec runtime parity probe failed: {message}")
        upgraded = load_yaml_json(upgraded_path, yq_bin, timeout_seconds)
        roots = upgraded.get("library_spec") or []
        observed_children = roots[0].get("regions", []) if len(roots) == 1 else []
        if upgraded.get("seqspec_version") != "0.5.0":
            raise ValueError("seqspec runtime parity probe did not upgrade to 0.5.0")
        if len(observed_children) != len(children):
            raise ValueError(
                "seqspec runtime parity probe changed the registry population"
            )
        for index, observed in enumerate(observed_children):
            label = labels[index]
            terms = region_type_values(observed.get("region_type"))
            if observed.get("region_id") != f"region_{index:03d}":
                raise ValueError("seqspec runtime parity probe changed region order")
            if terms != legacy_map[label]:
                raise ValueError(
                    f"seqspec runtime mapping differs from registry: {label}; "
                    f"expected {legacy_map[label]}, observed {terms}"
                )
        return {
            "legacy_labels": len(labels),
            "matched_labels": len(labels),
            "probe_sha256": runtime.file_sha256(source),
            "upgraded_probe_sha256": runtime.file_sha256(upgraded_path),
        }


def validate_outputs(
    *,
    survey_id: str,
    specifications: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    query_rows: list[dict[str, Any]],
    sample: list[dict[str, Any]],
    protocol: dict[str, Any],
    runtime_parity: dict[str, Any],
) -> dict[str, Any]:
    errors = []
    if any(value["survey_id"] != survey_id for value in rows + query_rows + sample):
        errors.append("survey identifiers differ")
    if runtime_parity["matched_labels"] != runtime_parity["legacy_labels"]:
        errors.append("seqspec runtime and registry mappings differ")
    expected_configs = {value["configuration_accession"] for value in specifications}
    observed_configs = {value["configuration_accession"] for value in rows}
    if observed_configs != expected_configs:
        errors.append("region rows do not cover every paired specification")
    query_configs = {value["configuration_accession"] for value in query_rows}
    if query_configs != expected_configs:
        errors.append("query rows do not cover every paired specification")
    unexpected_queries = sum(value["unexpected_change"] for value in query_rows)
    if (
        unexpected_queries
        > protocol["scientific_targets"]["unexpected_query_changes_max"]
    ):
        errors.append("unexpected query changes exceed the scientific target")
    if len({value["region_key"] for value in sample}) != len(sample):
        errors.append("review sample region keys are not unique")
    if not {value["region_key"] for value in sample} <= {
        value["region_key"] for value in rows
    }:
        errors.append("review sample is not a subset of region rows")
    expected_sample = min(protocol["manual_review"]["sample_size"], len(rows))
    if len(sample) != expected_sample:
        errors.append("review sample size differs from protocol")
    if any(
        not 0 < value["inclusion_probability"] <= 1
        or not math.isclose(
            value["sampling_weight"], 1 / value["inclusion_probability"]
        )
        for value in sample
    ):
        errors.append("review sampling weights are invalid")
    return {
        "schema_version": SCHEMA_VERSION,
        "survey_id": survey_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "counts": {
            "specifications": len(specifications),
            "registry_legacy_labels": runtime_parity["legacy_labels"],
            "runtime_registry_mappings": runtime_parity["matched_labels"],
            "regions": len(rows),
            "original_region_types": len(
                {
                    label
                    for value in rows
                    for label in split_values(value["original_region_type_labels"])
                }
            ),
            "ontology_terms": len(
                {
                    term
                    for value in rows
                    for term in split_values(value["ontology_terms"])
                }
            ),
            "unknown_regions": sum(value["unknown_mapping"] for value in rows),
            "multi_term_regions": sum(value["multi_term_mapping"] for value in rows),
            "query_comparisons": len(query_rows),
            "query_changes": sum(
                value["change_class"] != "stable" for value in query_rows
            ),
            "expected_unknown_query_changes": sum(
                value["change_class"] == "expected_unknown_collapse"
                for value in query_rows
            ),
            "unexpected_query_changes": unexpected_queries,
            "review_sample": len(sample),
            "review_strata": dict(
                sorted(Counter(value["review_stratum"] for value in sample).items())
            ),
        },
    }


def load_yaml_json(path: Path, yq_bin: Path, timeout_seconds: int) -> dict[str, Any]:
    completed = subprocess.run(
        [str(yq_bin), "-o=json", ".", str(path)],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip()
        raise ValueError(f"could not parse {path}: {message}")
    try:
        value = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise ValueError(f"could not parse yq JSON for {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"expected a YAML object in {path}")
    return value


def map_region_types(
    labels: list[str], legacy_map: dict[str, list[str]]
) -> tuple[list[str], str]:
    terms = []
    sources = set()
    for label in labels:
        if label.startswith("RGN:"):
            mapped = [label]
            sources.add("ontology_passthrough")
        elif label in legacy_map:
            mapped = legacy_map[label]
            sources.add("deterministic_registry")
        else:
            mapped = ["RGN:unknown:unclassified"]
            sources.add("unknown_fallback")
        for term in mapped:
            if term not in terms:
                terms.append(term)
    source = next(iter(sources)) if len(sources) == 1 else "mixed"
    return terms, source


def region_type_values(value: Any) -> list[str]:
    values = value if isinstance(value, list) else [value]
    if not values or any(not isinstance(item, str) or not item for item in values):
        raise ValueError("region_type is empty or invalid")
    return values


def split_values(value: str) -> list[str]:
    return [item for item in str(value).split(";") if item]


def verify_hash(path: Path, expected: str, label: str) -> None:
    if runtime.file_sha256(path) != expected:
        raise ValueError(f"{label} hash changed")


def positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
