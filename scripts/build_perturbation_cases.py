#!/usr/bin/env python3
"""Build deterministic perturbation applicability records from selected cases."""

from __future__ import annotations

import argparse
import json
import os
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
OPERATOR_IDS = tuple(
    [f"S{index:02d}" for index in range(1, 11)]
    + [f"D{index:02d}" for index in range(1, 5)]
)
SELECTION_IDENTITY_FIELDS = (
    "schema_version",
    "freeze_id",
    "cohort_split",
    "cohort_sha256",
    "sampling_protocol_sha256",
    "selector",
    "seqspec",
    "cases",
)
ROW_FIELDS = (
    "inventory_id",
    "record_id",
    "selection_id",
    "family_id",
    "configuration_accession",
    "modality",
    "operator_id",
    "operator_name",
    "layer",
    "stochastic",
    "applicable",
    "reason",
    "variants",
    "case_ids",
    "fastq_accessions",
    "read_ids",
    "region_ids",
    "target_json",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build one perturbation applicability row per configuration/operator."
    )
    parser.add_argument("--case-selection-manifest", required=True, type=Path)
    parser.add_argument("--perturbation-protocol", required=True, type=Path)
    parser.add_argument("--seqspec-bin", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=120)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = build_perturbation_cases(
            case_selection_manifest_path=args.case_selection_manifest.resolve(),
            perturbation_protocol_path=args.perturbation_protocol.resolve(),
            seqspec_bin=args.seqspec_bin.resolve(),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
        )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"build_perturbation_cases: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def build_perturbation_cases(
    *,
    case_selection_manifest_path: Path,
    perturbation_protocol_path: Path,
    seqspec_bin: Path,
    output_root: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    selection, cases = load_case_selection(case_selection_manifest_path)
    protocol = runtime.load_json(perturbation_protocol_path)
    validate_protocol(protocol)
    if not seqspec_bin.is_file():
        raise ValueError(f"seqspec executable does not exist: {seqspec_bin}")
    seqspec = runtime.executable_identity(
        seqspec_bin, timeout_seconds=min(timeout_seconds, 30)
    )
    builder = runtime.script_identity(Path(__file__).resolve(), version=BUILDER_VERSION)
    runtime_tool = runtime.script_identity(
        Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
    )

    grouped = group_configuration_cases(cases)
    configuration_records = [
        inspect_configuration(
            configuration_cases,
            protocol=protocol,
            seqspec_bin=seqspec_bin,
            timeout_seconds=timeout_seconds,
            selection_id=selection["selection_id"],
        )
        for configuration_cases in grouped
    ]
    stable_records = [
        record
        for configuration in configuration_records
        for record in configuration["operators"]
    ]
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": selection["selection_id"],
        "perturbation_protocol_sha256": runtime.file_sha256(perturbation_protocol_path),
        "seqspec": runtime.functional_executable_identity(seqspec),
        "builder": runtime.functional_script_identity(builder),
        "runtime": runtime.functional_script_identity(runtime_tool),
        "records": stable_records,
    }
    inventory_id = runtime.sha256_json(stable_identity)[:16]
    records = [{**record, "inventory_id": inventory_id} for record in stable_records]
    rows = [public_row(record, protocol) for record in records]
    validation = validate_inventory(
        inventory_id=inventory_id,
        selection_id=selection["selection_id"],
        configurations=configuration_records,
        records=records,
        protocol=protocol,
    )
    if not validation["valid"]:
        raise ValueError(
            "perturbation applicability inventory failed: "
            + "; ".join(validation["errors"] or ["unknown validation error"])
        )

    output_root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_root.name}-", dir=output_root.parent
    ) as tmpdir:
        temporary_root = Path(tmpdir)
        table_path = temporary_root / "tables" / "perturbation_applicability.csv"
        cases_path = temporary_root / "manifests" / "perturbation_cases.json"
        validation_path = temporary_root / "validation" / "perturbation_cases.json"
        manifest_path = temporary_root / "manifests" / "perturbation_inventory.json"
        runtime.write_csv(table_path, rows, list(ROW_FIELDS))
        runtime.write_json(
            cases_path,
            {
                "schema_version": SCHEMA_VERSION,
                "inventory_id": inventory_id,
                "selection_id": selection["selection_id"],
                "records": records,
            },
        )
        runtime.write_json(validation_path, validation)
        final_manifest_path = output_root / manifest_path.relative_to(temporary_root)
        manifest = {
            **stable_identity,
            "inventory_id": inventory_id,
            "created_at": utc_now(),
            "valid": True,
            "inputs": {
                "case_selection_manifest": runtime.file_identity(
                    case_selection_manifest_path
                ),
                "perturbation_protocol": runtime.file_identity(
                    perturbation_protocol_path
                ),
            },
            "tools": {
                "seqspec": seqspec,
                "builder": builder,
                "runtime": runtime_tool,
            },
            "counts": validation["counts"],
            "outputs": {
                "applicability": final_file_identity(
                    output_root, temporary_root, table_path
                ),
                "cases": final_file_identity(output_root, temporary_root, cases_path),
                "validation": final_file_identity(
                    output_root, temporary_root, validation_path
                ),
            },
            "manifest_path": str(final_manifest_path),
        }
        runtime.write_json(manifest_path, manifest)
        os.replace(temporary_root, output_root)
    return manifest


def load_case_selection(
    manifest_path: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    selection = runtime.load_json(manifest_path)
    if selection.get("valid") is not True:
        raise ValueError("sampling case selection is not valid")
    if selection.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling case selection schema is unsupported")
    selection_id = required_string(selection, "selection_id")
    stable_selection = {
        field: selection.get(field) for field in SELECTION_IDENTITY_FIELDS
    }
    if any(stable_selection[field] is None for field in SELECTION_IDENTITY_FIELDS):
        raise ValueError("sampling case selection identity is incomplete")
    if runtime.sha256_json(stable_selection)[:16] != selection_id:
        raise ValueError("sampling case selection id is not content-addressed")

    identity = selection.get("outputs", {}).get("cases_json", {})
    if not isinstance(identity, dict):
        raise ValueError("sampling case selection has no cases identity")
    cases_path = Path(required_string(identity, "path")).resolve()
    if runtime.file_sha256(cases_path) != required_string(identity, "sha256"):
        raise ValueError("sampling cases JSON hash changed")
    payload = runtime.load_json(cases_path)
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling cases schema is unsupported")
    if payload.get("selection_id") != selection_id:
        raise ValueError("sampling cases selection id differs")
    if payload.get("freeze_id") != selection.get("freeze_id"):
        raise ValueError("sampling cases freeze id differs")
    cases = payload.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("sampling cases must be a non-empty list")
    stable_cases = [
        {
            key: value
            for key, value in case.items()
            if key not in {"selection_id", "case_id", "spec_path"}
        }
        for case in cases
        if isinstance(case, dict)
    ]
    if len(stable_cases) != len(cases) or stable_cases != selection["cases"]:
        raise ValueError("sampling case rows differ from content-addressed selection")
    if selection.get("counts", {}).get("fastq_cases") != len(cases):
        raise ValueError("sampling case count does not reconcile")

    case_ids = []
    for case in cases:
        if case.get("selection_id") != selection_id:
            raise ValueError("sampling case selection id differs")
        case_id = required_string(case, "case_id")
        case_ids.append(case_id)
        for field in (
            "family_id",
            "configuration_accession",
            "modality",
            "selection_role",
            "read_id",
            "fastq_accession",
            "fastq_url",
            "spec_path",
            "spec_sha256",
        ):
            required_string(case, field)
        spec_path = Path(case["spec_path"]).resolve()
        if runtime.file_sha256(spec_path) != case["spec_sha256"]:
            raise ValueError(f"{case_id}: spec hash changed")
    if len(case_ids) != len(set(case_ids)):
        raise ValueError("sampling case identifiers are not unique")
    return selection, cases


def validate_protocol(protocol: dict[str, Any]) -> None:
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("perturbation protocol does not use schema 0.1.0")
    required_string(protocol, "experiment_id")
    principles = protocol.get("principles", {})
    for field in (
        "one_conceptual_property_per_condition",
        "preserve_fastq_record_names",
        "require_declared_seqspec_check_outcome",
        "require_target_file_read_region",
    ):
        if principles.get(field) is not True:
            raise ValueError(f"perturbation principle is not enabled: {field}")
    if principles.get("preserve_read_length_except_for_operator") != "D01":
        raise ValueError("only D01 may alter read length")
    stochastic = protocol.get("stochastic_conditions", {})
    fractions = stochastic.get("event_fractions")
    seeds = stochastic.get("seeds")
    if (
        not isinstance(fractions, list)
        or not fractions
        or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not 0 < value < 1
            for value in fractions
        )
        or len(fractions) != len(set(fractions))
    ):
        raise ValueError("perturbation event fractions are invalid")
    if (
        not isinstance(seeds, list)
        or not seeds
        or any(type(value) is not int or value < 0 for value in seeds)
        or len(seeds) != len(set(seeds))
    ):
        raise ValueError("perturbation seeds are invalid")
    if stochastic.get("selection_method") != "hash_ranked_record_ids":
        raise ValueError("unsupported stochastic record selection method")
    base_sample = protocol.get("base_sample", {})
    if base_sample.get("require_frozen_sampling_policy") is not True:
        raise ValueError("perturbations must require a frozen sampling policy")
    if (
        type(base_sample.get("reservoir_seed")) is not int
        or base_sample["reservoir_seed"] < 0
    ):
        raise ValueError("perturbation base reservoir seed is invalid")

    operators = protocol.get("operators")
    if not isinstance(operators, list):
        raise ValueError("perturbation operators must be a list")
    ids = [required_string(operator, "operator_id") for operator in operators]
    if tuple(ids) != OPERATOR_IDS:
        raise ValueError("perturbation operator ids or ordering differ")
    for operator in operators:
        operator_id = operator["operator_id"]
        for field in ("name", "layer", "mutation", "expected_seqspec_check"):
            required_string(operator, field)
        if operator["expected_seqspec_check"] not in {
            "pass",
            "failure",
            "not_applicable",
        }:
            raise ValueError(f"{operator_id}: seqspec check outcome is invalid")
        if operator.get("stochastic") is not operator_id.startswith("D"):
            raise ValueError(f"{operator_id}: stochastic classification differs")
        variants = operator.get("variants")
        if (
            not isinstance(variants, list)
            or not variants
            or len(variants) != len(set(variants))
        ):
            raise ValueError(f"{operator_id}: variants are invalid")
        applicability = operator.get("applicability", {})
        if not isinstance(applicability.get("required"), list):
            raise ValueError(f"{operator_id}: applicability requirements are invalid")
        target = operator.get("target", {})
        if target.get("process_outcome") not in {"success", "failure"}:
            raise ValueError(f"{operator_id}: process outcome is invalid")
        for field in (
            "checks",
            "assessment_codes",
            "metric_names",
            "localization",
        ):
            if not isinstance(target.get(field), list):
                raise ValueError(f"{operator_id}: target {field} is invalid")


def group_configuration_cases(
    cases: list[dict[str, Any]],
) -> list[list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for case in cases:
        grouped[case["configuration_accession"]].append(case)
    result = []
    for accession in sorted(grouped):
        values = sorted(
            grouped[accession],
            key=lambda value: (
                0 if value["selection_role"] == "primary" else 1,
                value["fastq_accession"],
            ),
        )
        for field in ("family_id", "modality", "spec_sha256", "spec_path"):
            if len({value[field] for value in values}) != 1:
                raise ValueError(f"{accession}: selected cases disagree on {field}")
        result.append(values)
    return result


def inspect_configuration(
    cases: list[dict[str, Any]],
    *,
    protocol: dict[str, Any],
    seqspec_bin: Path,
    timeout_seconds: int,
    selection_id: str,
) -> dict[str, Any]:
    first = cases[0]
    spec_path = Path(first["spec_path"]).resolve()
    run_command([str(seqspec_bin), "check", str(spec_path)], timeout_seconds)
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
    modality = first["modality"]
    read_index = run_command(
        [str(seqspec_bin), "index", "-m", modality, str(spec_path)],
        timeout_seconds,
    ).stdout
    region_index = run_command(
        [
            str(seqspec_bin),
            "index",
            "-s",
            "region",
            "-m",
            modality,
            str(spec_path),
        ],
        timeout_seconds,
    ).stdout
    reads = parse_reads(sequence_spec, modality, first["configuration_accession"])
    leaves = extract_leaf_regions(
        library_spec, modality, spec_path, first["configuration_accession"]
    )
    label_regions = parse_region_labels(
        region_index, leaves, first["configuration_accession"]
    )
    projected = parse_read_index(
        read_index, label_regions, first["configuration_accession"]
    )
    profiles = []
    for case in cases:
        read_id = case["read_id"]
        if read_id not in reads:
            raise ValueError(
                f"{first['configuration_accession']}: selected read {read_id} is absent"
            )
        read = reads[read_id]
        profiles.append(
            {
                "case": case,
                "read": read,
                "projected": projected.get(read_id, []),
                "primer": next(
                    (leaf for leaf in leaves if leaf["region_id"] == read["primer_id"]),
                    None,
                ),
            }
        )
    context = {
        "selection_id": selection_id,
        "family_id": first["family_id"],
        "configuration_accession": first["configuration_accession"],
        "modality": modality,
        "spec_sha256": first["spec_sha256"],
        "reads": reads,
        "profiles": profiles,
        "leaves": leaves,
    }
    operators = [
        build_operator_record(operator, context) for operator in protocol["operators"]
    ]
    return {
        "family_id": first["family_id"],
        "configuration_accession": first["configuration_accession"],
        "modality": modality,
        "operators": operators,
    }


def parse_reads(
    value: Any, modality: str, configuration_accession: str
) -> dict[str, dict[str, Any]]:
    if not isinstance(value, list):
        raise ValueError(f"{configuration_accession}: sequence_spec is not a list")
    reads = {}
    for read in value:
        if not isinstance(read, dict) or read.get("modality") != modality:
            continue
        read_id = required_string(read, "read_id")
        if read_id in reads:
            raise ValueError(f"{configuration_accession}: duplicate read id {read_id}")
        if read.get("strand") not in {"pos", "neg"}:
            raise ValueError(f"{configuration_accession}: unsupported read strand")
        reads[read_id] = read
    if not reads:
        raise ValueError(f"{configuration_accession}: modality has no reads")
    return reads


def extract_leaf_regions(
    value: Any,
    modality: str,
    spec_path: Path,
    configuration_accession: str,
) -> list[dict[str, Any]]:
    roots = value.get(modality) if isinstance(value, dict) else value
    if isinstance(roots, dict):
        roots = [roots]
    if not isinstance(roots, list):
        raise ValueError(f"{configuration_accession}: library_spec is not a list")
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
        result = dict(region)
        result["region_id"] = required_string(region, "region_id")
        result["name"] = required_string(region, "name")
        result["sequence_type"] = required_string(region, "sequence_type")
        result["min_len"] = strict_int(region.get("min_len"), "region min_len")
        result["max_len"] = strict_int(region.get("max_len"), "region max_len")
        result["resource"] = onlist_resource(result, spec_path)
        leaves.append(result)
    return leaves


def parse_region_labels(
    output: str,
    leaves: list[dict[str, Any]],
    configuration_accession: str,
) -> dict[tuple[str, str], dict[str, Any]]:
    rows = [line.split("\t") for line in output.splitlines() if line.strip()]
    if any(len(fields) != 5 for fields in rows) or len(rows) != len(leaves):
        raise ValueError(f"{configuration_accession}: region index does not reconcile")
    result = {}
    for fields, leaf in zip(rows, leaves, strict=True):
        _, name, tool_label, _, _ = fields
        if name != leaf["name"]:
            raise ValueError(f"{configuration_accession}: region index order differs")
        key = (name, tool_label)
        if key in result:
            raise ValueError(f"{configuration_accession}: ambiguous region index label")
        result[key] = leaf
    return result


def parse_read_index(
    output: str,
    label_regions: dict[tuple[str, str], dict[str, Any]],
    configuration_accession: str,
) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for line_number, line in enumerate(output.splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 5:
            raise ValueError(
                f"{configuration_accession}: invalid read index line {line_number}"
            )
        read_id, name, tool_label, start_value, stop_value = fields
        region = label_regions.get((name, tool_label))
        if region is None:
            raise ValueError(f"{configuration_accession}: read index label is unknown")
        start = strict_int(start_value, "region start")
        stop = strict_int(stop_value, "region stop")
        if start < 0 or stop <= start:
            raise ValueError(f"{configuration_accession}: projected span is invalid")
        result[read_id].append({**region, "start": start, "stop": stop})
    return result


def build_operator_record(
    operator: dict[str, Any], context: dict[str, Any]
) -> dict[str, Any]:
    operator_id = operator["operator_id"]
    applicable, reason, variants, target = operator_target(operator_id, context)
    stable = {
        "selection_id": context["selection_id"],
        "family_id": context["family_id"],
        "configuration_accession": context["configuration_accession"],
        "modality": context["modality"],
        "spec_sha256": context["spec_sha256"],
        "operator_id": operator_id,
        "operator_name": operator["name"],
        "layer": operator["layer"],
        "stochastic": operator["stochastic"],
        "applicable": applicable,
        "reason": reason,
        "variants": variants,
        "target": target,
    }
    return {**stable, "record_id": runtime.sha256_json(stable)[:16]}


def operator_target(
    operator_id: str, context: dict[str, Any]
) -> tuple[bool, str, list[str], dict[str, Any]]:
    profiles = context["profiles"]
    leaves = context["leaves"]
    first = profiles[0]
    if operator_id == "S01":
        variants = ["unexpected_name"]
        ambiguous_basename = find_ambiguous_input_name(context["reads"])
        if ambiguous_basename is not None:
            variants.append("ambiguous_name")
        return applicable_target(
            variants,
            profile_target(
                [first],
                [],
                {
                    "ambiguous_basename": ambiguous_basename or "",
                    "unexpected_basename": "seqcheck_unexpected_input.fastq.gz",
                },
            ),
            "selected FASTQ available",
        )
    if operator_id in {"S02", "S03"}:
        if len(profiles) < 2:
            return inapplicable("fewer than two selected FASTQs")
        if (
            operator_id == "S02"
            and len({profile["read"]["read_id"] for profile in profiles}) < 2
        ):
            return inapplicable("selected FASTQs do not map to distinct reads")
        variants = ["pairwise_swap"] if operator_id == "S02" else ["omit_one"]
        target_profiles = profiles[:2] if operator_id == "S02" else [profiles[-1]]
        return applicable_target(
            variants,
            profile_target(target_profiles, []),
            "multiple selected FASTQs available",
        )
    if operator_id == "S04":
        return applicable_target(
            ["exclude_observed_length"],
            profile_target(
                [first],
                [],
                {
                    "original_min_len": strict_int(
                        first["read"].get("min_len"), "read min_len"
                    ),
                    "original_max_len": strict_int(
                        first["read"].get("max_len"), "read max_len"
                    ),
                },
            ),
            "selected read has mutable length bounds",
        )
    if operator_id == "S05":
        boundary = select_boundary(profiles)
        if boundary is None:
            return inapplicable("no adjacent projected leaves support a positive shift")
        profile, left, right, direction, deltas = boundary
        return applicable_target(
            [f"shift_{delta}" for delta in deltas],
            profile_target(
                [profile],
                [left, right],
                {
                    "direction": direction,
                    "deltas": deltas,
                    "combined_length": (left["stop"] - left["start"])
                    + (right["stop"] - right["start"]),
                },
            ),
            "adjacent projected leaves support a length-preserving boundary shift",
        )
    if operator_id == "S06":
        selected = select_profile_region(
            profiles,
            concrete_fixed,
            require_non_palindromic=True,
            profile_predicate=lambda profile: strand_flip_supported(profile, leaves),
        )
        if selected is None:
            return inapplicable(
                "no read has opposite-strand range and a projected non-palindromic fixed region"
            )
        profile, region = selected
        variant = "pos_to_neg" if profile["read"]["strand"] == "pos" else "neg_to_pos"
        return applicable_target(
            [variant],
            profile_target([profile], [region]),
            "read strand and orientation-sensitive fixed region available",
        )
    if operator_id == "S07":
        for profile in profiles:
            alternatives = sorted(
                (
                    leaf
                    for leaf in leaves
                    if concrete_fixed(leaf)
                    and leaf["region_id"] != profile["read"]["primer_id"]
                    and primer_anchor_supported(profile, leaves, leaf)
                ),
                key=lambda value: value["region_id"],
            )
            if alternatives:
                return applicable_target(
                    ["alternate_fixed_primer"],
                    profile_target(
                        [profile],
                        [alternatives[0]],
                        {
                            "original_primer_id": profile["read"]["primer_id"],
                            "replacement_primer_id": alternatives[0]["region_id"],
                        },
                    ),
                    "alternate fixed scannable primer is available",
                )
        return inapplicable(
            "no alternate fixed primer preserves the read sequence-able range"
        )
    if operator_id == "S08":
        selected = select_profile_region(
            profiles, concrete_fixed, require_non_palindromic=True
        )
        if selected is None:
            return inapplicable("no projected non-palindromic fixed region")
        profile, region = selected
        return applicable_target(
            ["deterministic_replacement", "reverse_complement"],
            profile_target([profile], [region]),
            "projected concrete fixed region available",
        )
    if operator_id in {"S09", "S10", "D03"}:
        local_only = operator_id == "S10"
        selected = select_profile_region(
            profiles,
            lambda region: usable_onlist(region, local_only=local_only),
        )
        if selected is None:
            requirement = "loadable local onlist" if local_only else "loadable onlist"
            return inapplicable(f"no projected {requirement} region")
        profile, region = selected
        variants = {
            "S09": ["disjoint_same_length_onlist"],
            "S10": ["missing_local_resource"],
            "D03": ["deterministic_offlist_substitution"],
        }[operator_id]
        return applicable_target(
            variants,
            profile_target([profile], [region]),
            "projected onlist region has a resolvable resource",
        )
    if operator_id == "D01":
        selected = select_truncation(profiles)
        if selected is None:
            return inapplicable("no projected region can be crossed by truncation")
        profile, region, truncate_to = selected
        return applicable_target(
            ["truncate_before_expected_end"],
            profile_target([profile], [region], {"truncate_to_bases": truncate_to}),
            "selected read has a positive within-region truncation point",
        )
    if operator_id == "D02":
        selected = select_profile_region(profiles, concrete_fixed)
        if selected is None:
            return inapplicable("no projected concrete fixed region")
        profile, region = selected
        return applicable_target(
            ["single_base_substitution"],
            profile_target([profile], [region], {"target_position": region["start"]}),
            "projected fixed region has a covered substitution position",
        )
    if operator_id == "D04":
        for profile in profiles:
            primer = profile["primer"]
            if (
                primer is not None
                and concrete_fixed(primer)
                and strict_int(profile["read"].get("max_len"), "read max_len")
                >= len(primer["sequence"])
            ):
                variants = ["introduce_by_substitution"]
                if non_palindromic(primer):
                    variants.append("reverse_complement_by_substitution")
                return applicable_target(
                    variants,
                    profile_target(
                        [profile],
                        [primer],
                        {
                            "conditional_variant": "remove_by_substitution_if_baseline_hit"
                        },
                    ),
                    "selected read has a fixed scannable primer and substitution window",
                )
        return inapplicable("selected reads have no fixed scannable primer")
    raise ValueError(f"unsupported perturbation operator: {operator_id}")


def select_boundary(
    profiles: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], str, list[int]] | None:
    candidates = []
    for profile in profiles:
        for left, right in zip(profile["projected"], profile["projected"][1:]):
            if not (
                left["sequence_type"] in {"fixed", "onlist", "random"}
                or right["sequence_type"] in {"fixed", "onlist", "random"}
            ):
                continue
            for donor, direction in ((left, "left_to_right"), (right, "right_to_left")):
                donor_len = donor["stop"] - donor["start"]
                donor_min_len = strict_int(donor.get("min_len"), "region min_len")
                deltas = [
                    delta
                    for delta in (1, 4, 8)
                    if donor_len > delta and donor_min_len > delta
                ]
                if deltas:
                    candidates.append(
                        (profile, left, right, direction, deltas, donor_len)
                    )
    if not candidates:
        return None
    profile, left, right, direction, deltas, _ = min(
        candidates,
        key=lambda value: (
            -len(value[4]),
            -value[5],
            0 if value[0]["case"]["selection_role"] == "primary" else 1,
            value[0]["case"]["fastq_accession"],
            value[1]["start"],
            value[1]["region_id"],
            value[2]["region_id"],
            value[3],
        ),
    )
    return profile, left, right, direction, deltas


def select_profile_region(
    profiles: list[dict[str, Any]],
    predicate: Any,
    *,
    require_non_palindromic: bool = False,
    profile_predicate: Any = None,
) -> tuple[dict[str, Any], dict[str, Any]] | None:
    candidates = []
    for profile in profiles:
        if profile_predicate is not None and not profile_predicate(profile):
            continue
        for region in profile["projected"]:
            if predicate(region) and (
                not require_non_palindromic or non_palindromic(region)
            ):
                candidates.append((profile, region))
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda value: (
            0 if value[0]["case"]["selection_role"] == "primary" else 1,
            value[0]["case"]["fastq_accession"],
            value[1]["start"],
            value[1]["region_id"],
        ),
    )


def select_truncation(
    profiles: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any], int] | None:
    candidates = []
    for profile in profiles:
        for region in profile["projected"]:
            truncate_to = region["stop"] - 1
            if truncate_to > 0 and truncate_to >= region["start"]:
                candidates.append((profile, region, truncate_to))
    if not candidates:
        return None
    return max(
        candidates,
        key=lambda value: (
            value[2],
            1 if value[0]["case"]["selection_role"] == "primary" else 0,
            value[1]["region_id"],
        ),
    )


def concrete_fixed(region: dict[str, Any]) -> bool:
    sequence = str(region.get("sequence", "")).upper()
    return (
        region.get("sequence_type") == "fixed"
        and bool(sequence)
        and all(base in "ACGT" for base in sequence)
    )


def strand_flip_supported(
    profile: dict[str, Any], leaves: list[dict[str, Any]]
) -> bool:
    primer_id = profile["read"]["primer_id"]
    positions = [
        index for index, leaf in enumerate(leaves) if leaf["region_id"] == primer_id
    ]
    if len(positions) != 1:
        return False
    primer_index = positions[0]
    strand = profile["read"]["strand"]
    opposite_leaves = (
        leaves[:primer_index] if strand == "pos" else leaves[primer_index + 1 :]
    )
    available = sum(
        strict_int(leaf.get("max_len"), "region max_len") for leaf in opposite_leaves
    )
    return available >= strict_int(profile["read"].get("max_len"), "read max_len")


def primer_anchor_supported(
    profile: dict[str, Any],
    leaves: list[dict[str, Any]],
    primer: dict[str, Any],
) -> bool:
    positions = [
        index
        for index, leaf in enumerate(leaves)
        if leaf["region_id"] == primer["region_id"]
    ]
    if len(positions) != 1:
        return False
    primer_index = positions[0]
    strand = profile["read"]["strand"]
    sequenceable = (
        leaves[primer_index + 1 :] if strand == "pos" else leaves[:primer_index]
    )
    available = sum(
        strict_int(leaf.get("max_len"), "region max_len") for leaf in sequenceable
    )
    return available >= strict_int(profile["read"].get("max_len"), "read max_len")


def non_palindromic(region: dict[str, Any]) -> bool:
    sequence = str(region.get("sequence", "")).upper()
    return sequence != reverse_complement(sequence)


def find_ambiguous_input_name(reads: dict[str, dict[str, Any]]) -> str | None:
    possible_names = set()
    for read_id, read in reads.items():
        possible_names.add(read_id)
        possible_names.add(normalize_fastq_name(read_id))
        for file in read.get("files", []):
            if not isinstance(file, dict):
                continue
            url_basename = Path(str(file.get("url", ""))).name
            for name in (
                str(file.get("file_id", "")).strip(),
                str(file.get("filename", "")).strip(),
                url_basename,
            ):
                if name:
                    possible_names.add(name)
                    possible_names.add(normalize_fastq_name(name))
    candidates = []
    for basename in sorted(name for name in possible_names if name):
        ranked = input_match_ranks(reads, basename)
        if not ranked:
            continue
        best_rank = min(ranked.values())
        if sum(rank == best_rank for rank in ranked.values()) > 1:
            candidates.append((best_rank, basename))
    return min(candidates)[1] if candidates else None


def input_match_ranks(
    reads: dict[str, dict[str, Any]], basename: str
) -> dict[str, int]:
    result = {}
    normalized = normalize_fastq_name(basename)
    for read_id, read in reads.items():
        ranks = []
        for file in read.get("files", []):
            if not isinstance(file, dict):
                continue
            names = (
                str(file.get("file_id", "")).strip(),
                str(file.get("filename", "")).strip(),
                Path(str(file.get("url", ""))).name,
            )
            ranks.extend(rank for rank, name in enumerate(names) if name == basename)
            if any(name and normalize_fastq_name(name) == normalized for name in names):
                ranks.append(4)
        if read_id == basename:
            ranks.append(3)
        if normalize_fastq_name(read_id) == normalized:
            ranks.append(5)
        if ranks:
            result[read_id] = min(ranks)
    return result


def normalize_fastq_name(value: str) -> str:
    result = value.lower()
    if result.endswith(".gz"):
        result = result[: -len(".gz")]
    for suffix in (".fastq", ".fq"):
        if result.endswith(suffix):
            return result[: -len(suffix)]
    return result


def usable_onlist(region: dict[str, Any], *, local_only: bool) -> bool:
    resource = region.get("resource", {})
    if region.get("sequence_type") != "onlist" or not resource.get("locator"):
        return False
    if local_only and resource.get("urltype") != "local":
        return False
    return resource.get("available") is True


def onlist_resource(region: dict[str, Any], spec_path: Path) -> dict[str, Any]:
    value = region.get("onlist")
    if not isinstance(value, dict):
        return {"urltype": "", "locator": "", "available": False}
    urltype = str(value.get("urltype", "")).strip()
    locator = str(value.get("url") or value.get("filename") or "").strip()
    if urltype == "local":
        path = (spec_path.parent / locator).resolve()
        result = {
            "urltype": urltype,
            "locator": locator,
            "available": path.is_file(),
        }
        if path.is_file():
            result["sha256"] = runtime.file_sha256(path)
            result["size_bytes"] = path.stat().st_size
        return result
    return {
        "urltype": urltype,
        "locator": locator,
        "available": bool(locator and urltype in {"http", "https", "ftp"}),
    }


def profile_target(
    profiles: list[dict[str, Any]],
    regions: list[dict[str, Any]],
    details: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "case_ids": [profile["case"]["case_id"] for profile in profiles],
        "fastq_accessions": [
            profile["case"]["fastq_accession"] for profile in profiles
        ],
        "read_ids": [profile["read"]["read_id"] for profile in profiles],
        "region_ids": [region["region_id"] for region in regions],
        "reads": [read_target(profile) for profile in profiles],
        "regions": [region_target(region) for region in regions],
        "details": details or {},
    }


def read_target(profile: dict[str, Any]) -> dict[str, Any]:
    read = profile["read"]
    return {
        "case_id": profile["case"]["case_id"],
        "fastq_accession": profile["case"]["fastq_accession"],
        "read_id": read["read_id"],
        "primer_id": read["primer_id"],
        "strand": read["strand"],
        "min_len": strict_int(read.get("min_len"), "read min_len"),
        "max_len": strict_int(read.get("max_len"), "read max_len"),
    }


def region_target(region: dict[str, Any]) -> dict[str, Any]:
    result = {
        "region_id": region["region_id"],
        "name": region["name"],
        "sequence_type": region["sequence_type"],
        "sequence": str(region.get("sequence", "")),
        "min_len": strict_int(region.get("min_len"), "region min_len"),
        "max_len": strict_int(region.get("max_len"), "region max_len"),
    }
    if "start" in region or "stop" in region:
        result["start"] = strict_int(region.get("start"), "region start")
        result["stop"] = strict_int(region.get("stop"), "region stop")
    resource = region.get("resource")
    if isinstance(resource, dict) and resource.get("locator"):
        result["resource"] = resource
    return result


def applicable_target(
    variants: list[str], target: dict[str, Any], reason: str
) -> tuple[bool, str, list[str], dict[str, Any]]:
    return True, reason, variants, target


def inapplicable(reason: str) -> tuple[bool, str, list[str], dict[str, Any]]:
    return False, reason, [], profile_target([], [])


def public_row(record: dict[str, Any], protocol: dict[str, Any]) -> dict[str, Any]:
    operator = next(
        value
        for value in protocol["operators"]
        if value["operator_id"] == record["operator_id"]
    )
    target = record["target"]
    return {
        **record,
        "variants": ";".join(record["variants"]),
        "case_ids": ";".join(target["case_ids"]),
        "fastq_accessions": ";".join(target["fastq_accessions"]),
        "read_ids": ";".join(target["read_ids"]),
        "region_ids": ";".join(target["region_ids"]),
        "target_json": runtime.canonical_json(
            {"selected": target, "expected": operator["target"]}
        ),
    }


def validate_inventory(
    *,
    inventory_id: str,
    selection_id: str,
    configurations: list[dict[str, Any]],
    records: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    errors = []
    expected_count = len(configurations) * len(OPERATOR_IDS)
    if len(records) != expected_count:
        errors.append(f"expected {expected_count} records, observed {len(records)}")
    expected_pairs = {
        (configuration["configuration_accession"], operator_id)
        for configuration in configurations
        for operator_id in OPERATOR_IDS
    }
    observed_pairs = {
        (record["configuration_accession"], record["operator_id"]) for record in records
    }
    if observed_pairs != expected_pairs or len(records) != len(observed_pairs):
        errors.append("configuration/operator matrix does not reconcile")
    if len({record["record_id"] for record in records}) != len(records):
        errors.append("perturbation record identifiers are not unique")
    operators = {
        operator["operator_id"]: operator for operator in protocol["operators"]
    }
    for record in records:
        if record.get("inventory_id") != inventory_id:
            errors.append(f"{record['record_id']}: inventory id differs")
        if record.get("selection_id") != selection_id:
            errors.append(f"{record['record_id']}: selection id differs")
        if record["applicable"] and not record["target"]["case_ids"]:
            errors.append(
                f"{record['record_id']}: applicable record has no target case"
            )
        if not record["applicable"] and not record["reason"]:
            errors.append(f"{record['record_id']}: inapplicable record has no reason")
        if not set(record["variants"]).issubset(
            operators[record["operator_id"]]["variants"]
        ):
            errors.append(f"{record['record_id']}: selected variant is undeclared")
        if record["applicable"]:
            localization = operators[record["operator_id"]]["target"]["localization"]
            for scope, target_field in (
                ("file", "fastq_accessions"),
                ("read", "read_ids"),
                ("region", "region_ids"),
            ):
                if scope in localization and not record["target"][target_field]:
                    errors.append(
                        f"{record['record_id']}: applicable target has no {scope}"
                    )
    applicable = [record for record in records if record["applicable"]]
    return {
        "schema_version": SCHEMA_VERSION,
        "inventory_id": inventory_id,
        "selection_id": selection_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "counts": {
            "configurations": len(configurations),
            "operators": len(protocol["operators"]),
            "records": len(records),
            "applicable": len(applicable),
            "inapplicable": len(records) - len(applicable),
            "applicable_by_operator": dict(
                sorted(Counter(record["operator_id"] for record in applicable).items())
            ),
            "applicable_by_family": dict(
                sorted(Counter(record["family_id"] for record in applicable).items())
            ),
        },
    }


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
        raise ValueError(
            f"command exited {completed.returncode}: {' '.join(argv)}"
            + (f": {message}" if message else "")
        )
    return completed


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


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


def required_string(value: dict[str, Any], field: str) -> str:
    result = str(value.get(field, "")).strip()
    if not result:
        raise ValueError(f"missing required field: {field}")
    return result


def final_file_identity(
    output_root: Path, temporary_root: Path, path: Path
) -> dict[str, Any]:
    return {
        "path": str(output_root / path.relative_to(temporary_root)),
        "sha256": runtime.file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
