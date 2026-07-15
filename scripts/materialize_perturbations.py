#!/usr/bin/env python3
"""Materialize controlled seqspec and FASTQ perturbation conditions."""

from __future__ import annotations

import argparse
import copy
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import build_perturbation_cases as inventory_builder
    import paper_runtime as runtime
    import perturbation_runtime as mutation
except ModuleNotFoundError:
    from scripts import build_perturbation_cases as inventory_builder
    from scripts import paper_runtime as runtime
    from scripts import perturbation_runtime as mutation


SCHEMA_VERSION = "0.1.0"
MATERIALIZER_VERSION = "0.1.0"
INVENTORY_IDENTITY_FIELDS = (
    "schema_version",
    "selection_id",
    "perturbation_protocol_sha256",
    "seqspec",
    "builder",
    "runtime",
    "records",
)
STUDY_IDENTITY_FIELDS = (
    "schema_version",
    "selection_id",
    "sampling_protocol_sha256",
    "seqcheck",
    "sampler",
    "sampler_core",
    "runtime",
    "runner",
)
BUNDLE_IDENTITY_FIELDS = (
    "schema_version",
    "cohort_split",
    "selection_id",
    "study_run_id",
    "sampling_policy_id",
    "case_selection_sha256",
    "sampling_study_sha256",
    "sampling_policy_sha256",
    "sampling_rule",
    "sampler",
    "runner",
    "runtime",
    "samples",
)
CONDITION_FIELDS = (
    "materialization_id",
    "condition_id",
    "selection_id",
    "configuration_accession",
    "family_id",
    "modality",
    "condition_kind",
    "operator_id",
    "variant",
    "event_fraction",
    "mutation_seed",
    "mutated_record_count",
    "total_record_count",
    "spec_path",
    "input_paths",
    "target_case_ids",
    "target_read_ids",
    "target_region_ids",
    "expected_seqspec_check",
    "observed_seqspec_check",
    "expected_process_outcome",
    "expected_checks",
    "expected_assessment_codes",
    "expected_metric_names",
    "condition_json",
)
SKIP_FIELDS = (
    "materialization_id",
    "selection_id",
    "configuration_accession",
    "family_id",
    "modality",
    "operator_id",
    "variant",
    "event_fraction",
    "mutation_seed",
    "reason_code",
    "reason",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Materialize clean and controlled perturbation conditions."
    )
    parser.add_argument("--perturbation-inventory", required=True, type=Path)
    samples = parser.add_mutually_exclusive_group(required=True)
    samples.add_argument("--sampling-study", type=Path)
    samples.add_argument("--sample-bundle", type=Path)
    parser.add_argument("--sampling-policy", required=True, type=Path)
    parser.add_argument("--perturbation-protocol", required=True, type=Path)
    parser.add_argument("--seqspec-bin", required=True, type=Path)
    parser.add_argument("--yq-bin", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=120)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = materialize_perturbations(
            inventory_manifest_path=args.perturbation_inventory.resolve(),
            study_manifest_path=(
                args.sampling_study.resolve() if args.sampling_study else None
            ),
            sample_bundle_path=(
                args.sample_bundle.resolve() if args.sample_bundle else None
            ),
            sampling_policy_path=args.sampling_policy.resolve(),
            perturbation_protocol_path=args.perturbation_protocol.resolve(),
            seqspec_bin=args.seqspec_bin.resolve(),
            yq_bin=args.yq_bin.resolve(),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
        )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"materialize_perturbations: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def materialize_perturbations(
    *,
    inventory_manifest_path: Path,
    sampling_policy_path: Path,
    perturbation_protocol_path: Path,
    seqspec_bin: Path,
    yq_bin: Path,
    output_root: Path,
    timeout_seconds: int,
    study_manifest_path: Path | None = None,
    sample_bundle_path: Path | None = None,
) -> dict[str, Any]:
    inventory_manifest_path = inventory_manifest_path.resolve()
    study_manifest_path = study_manifest_path.resolve() if study_manifest_path else None
    sample_bundle_path = sample_bundle_path.resolve() if sample_bundle_path else None
    sampling_policy_path = sampling_policy_path.resolve()
    perturbation_protocol_path = perturbation_protocol_path.resolve()
    seqspec_bin = seqspec_bin.resolve()
    yq_bin = yq_bin.resolve()
    output_root = output_root.resolve()
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    if (study_manifest_path is None) == (sample_bundle_path is None):
        raise ValueError("provide exactly one sampling study or sample bundle")
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    inventory, inventory_records = load_inventory(inventory_manifest_path)
    protocol = runtime.load_json(perturbation_protocol_path)
    inventory_builder.validate_protocol(protocol)
    if inventory["perturbation_protocol_sha256"] != runtime.file_sha256(
        perturbation_protocol_path
    ):
        raise ValueError("perturbation protocol differs from applicability inventory")
    if study_manifest_path is not None:
        study, cases, clean_inputs = load_sampling_study(
            study_manifest_path,
            selection_id=inventory["selection_id"],
        )
        study_run_id = study["study_run_id"]
        sample_source = {
            "kind": "sampling_study",
            "id": study_run_id,
            "sha256": runtime.file_sha256(study_manifest_path),
        }
        selected_inputs = None
    else:
        bundle, cases, selected_inputs = load_sample_bundle(
            sample_bundle_path,
            selection_id=inventory["selection_id"],
        )
        study_run_id = bundle["study_run_id"]
        sample_source = {
            "kind": "policy_sample_bundle",
            "id": bundle["bundle_id"],
            "sha256": runtime.file_sha256(sample_bundle_path),
        }
    policy = load_sampling_policy(
        sampling_policy_path,
        study_run_id=study_run_id,
        require_frozen=protocol["base_sample"]["require_frozen_sampling_policy"],
    )
    if sample_bundle_path is not None:
        if bundle["sampling_policy_id"] != policy["policy_id"]:
            raise ValueError("sample bundle uses a different sampling policy")
        if bundle["sampling_policy_sha256"] != runtime.file_sha256(
            sampling_policy_path
        ):
            raise ValueError("sample bundle sampling policy hash differs")
        if bundle["sampling_rule"] != policy["default"]:
            raise ValueError("sample bundle sampling rule differs from policy")
    policy_seed = policy["default"]["sampling_seed"]
    expected_seed = (
        protocol["base_sample"]["reservoir_seed"]
        if policy["default"]["sampling_method"] == "reservoir"
        else None
    )
    if policy_seed != expected_seed:
        raise ValueError("perturbation base-sample seed differs from sampling policy")
    if selected_inputs is None:
        selected_inputs = select_clean_inputs(clean_inputs, policy=policy)
    case_map = {case["case_id"]: case for case in cases}
    if set(selected_inputs) != set(case_map):
        raise ValueError("sample source does not contain every selected FASTQ case")

    seqspec_identity = require_executable(seqspec_bin, timeout_seconds)
    if runtime.functional_executable_identity(seqspec_identity) != inventory["seqspec"]:
        raise ValueError("seqspec executable differs from applicability inventory")
    yq_identity = require_executable(yq_bin, timeout_seconds)
    tools = tool_identities(seqspec_identity, yq_identity)
    stable_identity = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": inventory["selection_id"],
        "inventory_id": inventory["inventory_id"],
        "study_run_id": study_run_id,
        "sampling_policy_id": policy["policy_id"],
        "sample_source": sample_source,
        "base_samples": [
            {
                "case_id": case_id,
                "source_id": value["source_id"],
                "condition_id": value["condition_id"],
                "sha256": value["sha256"],
                "size_bytes": value["size_bytes"],
                "records_selected": value["records_selected"],
            }
            for case_id, value in sorted(selected_inputs.items())
        ],
        "perturbation_protocol_sha256": runtime.file_sha256(perturbation_protocol_path),
        "seqspec": runtime.functional_executable_identity(seqspec_identity),
        "yq": runtime.functional_executable_identity(yq_identity),
        "materializer": runtime.functional_script_identity(tools["materializer"]),
        "mutation_runtime": runtime.functional_script_identity(
            tools["mutation_runtime"]
        ),
        "paper_runtime": runtime.functional_script_identity(tools["paper_runtime"]),
        "fastq_runtime": runtime.functional_script_identity(tools["fastq_runtime"]),
    }
    materialization_id = runtime.sha256_json(stable_identity)[:16]
    operators = {value["operator_id"]: value for value in protocol["operators"]}

    output_root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_root.name}-", dir=output_root.parent
    ) as tmpdir:
        temporary_root = Path(tmpdir).resolve()
        condition_rows = []
        skip_rows = []
        grouped_records = group_inventory_records(inventory_records)
        for configuration_accession, records in grouped_records.items():
            configuration_cases = sorted(
                (
                    case
                    for case in cases
                    if case["configuration_accession"] == configuration_accession
                ),
                key=lambda value: (
                    0 if value["selection_role"] == "primary" else 1,
                    value["fastq_accession"],
                ),
            )
            if not configuration_cases:
                raise ValueError(f"{configuration_accession}: no sampling cases")
            context = build_configuration_context(
                configuration_cases,
                selected_inputs=selected_inputs,
                yq_bin=yq_bin,
                timeout_seconds=timeout_seconds,
            )
            condition_rows.append(
                clean_condition(
                    context,
                    materialization_id=materialization_id,
                )
            )
            for record in records:
                operator = operators[record["operator_id"]]
                generated = generate_operator_conditions(
                    record=record,
                    operator=operator,
                    protocol=protocol,
                    context=context,
                    materialization_id=materialization_id,
                    temporary_root=temporary_root,
                    output_root=output_root,
                    seqspec_bin=seqspec_bin,
                    timeout_seconds=timeout_seconds,
                )
                condition_rows.extend(generated["conditions"])
                skip_rows.extend(generated["skips"])

        validation = validate_materialization(
            materialization_id=materialization_id,
            inventory_records=inventory_records,
            protocol=protocol,
            conditions=condition_rows,
            skips=skip_rows,
            temporary_root=temporary_root,
            output_root=output_root,
        )
        if not validation["valid"]:
            raise ValueError(
                "perturbation materialization failed: "
                + "; ".join(validation["errors"] or ["unknown validation error"])
            )
        conditions_path = temporary_root / "manifests" / "conditions.json"
        conditions_table_path = temporary_root / "tables" / "conditions.csv"
        skips_path = temporary_root / "tables" / "skipped_conditions.csv"
        validation_path = temporary_root / "validation" / "materialization.json"
        manifest_path = temporary_root / "manifests" / "materialization.json"
        runtime.write_json(
            conditions_path,
            {
                "schema_version": SCHEMA_VERSION,
                "materialization_id": materialization_id,
                "conditions": condition_rows,
                "skips": skip_rows,
            },
        )
        runtime.write_csv(
            conditions_table_path,
            [public_condition_row(row) for row in condition_rows],
            list(CONDITION_FIELDS),
        )
        runtime.write_csv(skips_path, skip_rows, list(SKIP_FIELDS))
        runtime.write_json(validation_path, validation)
        final_manifest_path = output_root / manifest_path.relative_to(temporary_root)
        manifest = {
            **stable_identity,
            "materialization_id": materialization_id,
            "created_at": utc_now(),
            "valid": True,
            "inputs": {
                "perturbation_inventory": runtime.file_identity(
                    inventory_manifest_path
                ),
                "sample_source": runtime.file_identity(
                    study_manifest_path or sample_bundle_path
                ),
                "sampling_policy": runtime.file_identity(sampling_policy_path),
                "perturbation_protocol": runtime.file_identity(
                    perturbation_protocol_path
                ),
            },
            "tools": tools,
            "base_sample": {
                **policy["default"],
            },
            "counts": validation["counts"],
            "outputs": {
                "conditions": final_file_identity(
                    output_root, temporary_root, conditions_path
                ),
                "conditions_table": final_file_identity(
                    output_root, temporary_root, conditions_table_path
                ),
                "skips": final_file_identity(output_root, temporary_root, skips_path),
                "validation": final_file_identity(
                    output_root, temporary_root, validation_path
                ),
            },
            "manifest_path": str(final_manifest_path),
        }
        runtime.write_json(manifest_path, manifest)
        os.replace(temporary_root, output_root)
    return manifest


def load_inventory(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    value = runtime.load_json(path)
    if value.get("valid") is not True:
        raise ValueError("perturbation inventory is not valid")
    stable = {field: value.get(field) for field in INVENTORY_IDENTITY_FIELDS}
    if any(stable[field] is None for field in INVENTORY_IDENTITY_FIELDS):
        raise ValueError("perturbation inventory identity is incomplete")
    if runtime.sha256_json(stable)[:16] != value.get("inventory_id"):
        raise ValueError("perturbation inventory id is not content-addressed")
    output = value.get("outputs", {}).get("cases", {})
    cases_path = verified_path(output, "perturbation cases")
    payload = runtime.load_json(cases_path)
    records = payload.get("records")
    if not isinstance(records, list) or not records:
        raise ValueError("perturbation cases are empty")
    if payload.get("inventory_id") != value["inventory_id"]:
        raise ValueError("perturbation cases inventory id differs")
    stable_records = [
        {key: item for key, item in record.items() if key != "inventory_id"}
        for record in records
        if isinstance(record, dict)
    ]
    if len(stable_records) != len(records) or stable_records != value["records"]:
        raise ValueError("perturbation case records differ from inventory")
    return value, records


def load_sampling_study(
    path: Path, *, selection_id: str
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, dict[str, Any]]]:
    study = runtime.load_json(path)
    if study.get("valid") is not True:
        raise ValueError("sampling study is not valid")
    stable = {field: study.get(field) for field in STUDY_IDENTITY_FIELDS}
    if any(stable[field] is None for field in STUDY_IDENTITY_FIELDS):
        raise ValueError("sampling study identity is incomplete")
    if runtime.sha256_json(stable)[:16] != study.get("study_run_id"):
        raise ValueError("sampling study id is not content-addressed")
    if study.get("selection_id") != selection_id:
        raise ValueError("sampling study selection differs from perturbation inventory")
    selection_identity = study.get("inputs", {}).get("case_selection_manifest", {})
    selection_path = verified_path(selection_identity, "sampling case selection")
    selection, cases = inventory_builder.load_case_selection(selection_path)
    if selection["selection_id"] != selection_id:
        raise ValueError("sampling case selection id differs")
    matrices = study.get("matrices")
    if not isinstance(matrices, list) or len(matrices) != len(cases):
        raise ValueError("sampling matrices do not reconcile with selected cases")
    result = {}
    for matrix_row in matrices:
        case_id = str(matrix_row.get("case_id", ""))
        if not case_id or case_id in result:
            raise ValueError("sampling matrix case identifiers are invalid")
        matrix_path = verified_path(matrix_row.get("manifest", {}), "sample matrix")
        matrix = runtime.load_json(matrix_path)
        if matrix.get("matrix_id") != matrix_row.get("matrix_id"):
            raise ValueError(f"{case_id}: sampling matrix id differs")
        for condition in matrix.get("conditions", []):
            for output in condition.get("outputs", []):
                output_path = Path(str(output.get("output_path", ""))).resolve()
                if runtime.file_sha256(output_path) != output.get("output_sha256"):
                    raise ValueError(f"{case_id}: sampled FASTQ hash changed")
        result[case_id] = matrix
    return study, cases, result


def load_sampling_policy(
    path: Path, *, study_run_id: str, require_frozen: bool
) -> dict[str, Any]:
    policy = runtime.load_json(path)
    if policy.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("sampling policy schema is unsupported")
    stable = {
        key: value
        for key, value in policy.items()
        if key not in {"policy_id", "created_at"}
    }
    if runtime.sha256_json(stable)[:16] != policy.get("policy_id"):
        raise ValueError("sampling policy id is not content-addressed")
    if policy.get("study_run_id") != study_run_id:
        raise ValueError("sampling policy refers to a different study")
    if require_frozen and policy.get("frozen") is not True:
        raise ValueError("sampling policy is not frozen")
    default = policy.get("default", {})
    if default.get("sampling_method") not in {"prefix", "reservoir"}:
        raise ValueError("sampling policy method is invalid")
    if positive_int(default.get("records_per_fastq"), "policy sample size") <= 0:
        raise ValueError("sampling policy size is invalid")
    seed = default.get("sampling_seed")
    if default["sampling_method"] == "reservoir":
        if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
            raise ValueError("sampling policy reservoir seed is invalid")
    elif seed is not None:
        raise ValueError("prefix sampling policy must not declare a seed")
    return policy


def load_sample_bundle(
    path: Path | None, *, selection_id: str
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, dict[str, Any]]]:
    if path is None:
        raise ValueError("sample bundle path is missing")
    bundle = runtime.load_json(path)
    if bundle.get("valid") is not True:
        raise ValueError("policy sample bundle is not valid")
    stable = {field: bundle.get(field) for field in BUNDLE_IDENTITY_FIELDS}
    if any(stable[field] is None for field in BUNDLE_IDENTITY_FIELDS):
        raise ValueError("policy sample bundle identity is incomplete")
    if runtime.sha256_json(stable)[:16] != bundle.get("bundle_id"):
        raise ValueError("policy sample bundle id is not content-addressed")
    if bundle.get("cohort_split") != "evaluation":
        raise ValueError("policy sample bundle is not from the evaluation split")
    if bundle.get("selection_id") != selection_id:
        raise ValueError("sample bundle selection differs from perturbation inventory")
    selection_path = verified_path(
        bundle.get("inputs", {}).get("case_selection_manifest", {}),
        "sample bundle case selection",
    )
    selection, cases = inventory_builder.load_case_selection(selection_path)
    if selection["selection_id"] != selection_id:
        raise ValueError("sample bundle case selection id differs")
    records = bundle.get("sample_records")
    if not isinstance(records, list) or len(records) != len(cases):
        raise ValueError("sample bundle records do not reconcile with selected cases")
    if [stable_bundle_sample(value) for value in records] != bundle["samples"]:
        raise ValueError("sample bundle records differ from content-addressed samples")
    verified_path(bundle.get("outputs", {}).get("samples", {}), "sample bundle table")
    validation_path = verified_path(
        bundle.get("outputs", {}).get("validation", {}),
        "sample bundle validation",
    )
    validation = runtime.load_json(validation_path)
    if (
        validation.get("valid") is not True
        or validation.get("bundle_id") != bundle["bundle_id"]
    ):
        raise ValueError("sample bundle validation is not valid")
    result = {}
    for record in records:
        case_id = str(record.get("case_id", ""))
        if not case_id or case_id in result:
            raise ValueError("sample bundle case identifiers are invalid")
        output_path = Path(str(record.get("output_path", ""))).resolve()
        if runtime.file_sha256(output_path) != record.get("output_sha256"):
            raise ValueError(f"{case_id}: policy sample FASTQ hash changed")
        manifest_path = Path(str(record.get("sample_manifest_path", ""))).resolve()
        if runtime.file_sha256(manifest_path) != record.get("sample_manifest_sha256"):
            raise ValueError(f"{case_id}: policy sample manifest hash changed")
        result[case_id] = {
            "source_id": bundle["bundle_id"],
            "condition_id": record["sample_id"],
            "method": record["sampling_method"],
            "requested_records_per_fastq": record["requested_records_per_fastq"],
            "seed": record["sampling_seed"],
            "path": str(output_path),
            "sha256": record["output_sha256"],
            "size_bytes": record["output_size_bytes"],
            "records_selected": record["records_selected"],
        }
    if set(result) != {case["case_id"] for case in cases}:
        raise ValueError("sample bundle does not contain every selected case")
    return bundle, cases, result


def stable_bundle_sample(value: dict[str, Any]) -> dict[str, Any]:
    return {
        key: item
        for key, item in value.items()
        if key
        not in {
            "bundle_id",
            "output_path",
            "sample_manifest_path",
            "sample_manifest_sha256",
        }
    }


def select_clean_inputs(
    matrices: dict[str, dict[str, Any]],
    *,
    policy: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    method = policy["default"]["sampling_method"]
    size = positive_int(policy["default"]["records_per_fastq"], "sample size")
    seed = policy["default"]["sampling_seed"]
    result = {}
    for case_id, matrix in matrices.items():
        matches = [
            condition
            for condition in matrix.get("conditions", [])
            if condition.get("method") == method
            and condition.get("requested_records_per_fastq") == size
            and condition.get("seed") == seed
        ]
        if len(matches) != 1 or len(matches[0].get("outputs", [])) != 1:
            raise ValueError(
                f"{case_id}: sampling matrix has no unique frozen base condition"
            )
        output = matches[0]["outputs"][0]
        path = Path(output["output_path"]).resolve()
        if int(output.get("records_selected", 0)) != size:
            raise ValueError(f"{case_id}: base sample has fewer records than required")
        result[case_id] = {
            "source_id": matrix["matrix_id"],
            "condition_id": matches[0]["condition_id"],
            "method": method,
            "requested_records_per_fastq": size,
            "seed": seed,
            "path": str(path),
            "sha256": runtime.file_sha256(path),
            "size_bytes": path.stat().st_size,
            "records_selected": size,
        }
    return result


def group_inventory_records(
    records: list[dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    grouped = defaultdict(list)
    for record in records:
        grouped[record["configuration_accession"]].append(record)
    return {
        key: sorted(values, key=lambda value: value["operator_id"])
        for key, values in sorted(grouped.items())
    }


def build_configuration_context(
    cases: list[dict[str, Any]],
    *,
    selected_inputs: dict[str, dict[str, Any]],
    yq_bin: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    first = cases[0]
    spec_path = Path(first["spec_path"]).resolve()
    spec = mutation.load_seqspec_json(
        spec_path, yq_bin=yq_bin, timeout_seconds=timeout_seconds
    )
    clean = []
    records = {}
    for case in cases:
        sample = selected_inputs[case["case_id"]]
        path = Path(sample["path"])
        values = mutation.read_fastq(path)
        if len(values) != sample["records_selected"]:
            raise ValueError(f"{case['case_id']}: sampled FASTQ record count differs")
        records[case["case_id"]] = values
        clean.append(clean_input_row(case, sample))
    if len({Path(value["path"]).name for value in clean}) != len(clean):
        raise ValueError(
            f"{first['configuration_accession']}: clean input basenames collide"
        )
    return {
        "selection_id": first["selection_id"],
        "configuration_accession": first["configuration_accession"],
        "family_id": first["family_id"],
        "modality": first["modality"],
        "spec_path": spec_path,
        "spec": spec,
        "cases": {case["case_id"]: case for case in cases},
        "clean_inputs": clean,
        "records": records,
        "onlists": {},
    }


def clean_input_row(case: dict[str, Any], sample: dict[str, Any]) -> dict[str, Any]:
    return {
        "case_id": case["case_id"],
        "fastq_accession": case["fastq_accession"],
        "read_id": case["read_id"],
        "path": sample["path"],
        "sha256": sample["sha256"],
        "size_bytes": sample["size_bytes"],
        "records": sample["records_selected"],
        "sample_source_id": sample["source_id"],
        "sample_condition_id": sample["condition_id"],
    }


def clean_condition(
    context: dict[str, Any], *, materialization_id: str
) -> dict[str, Any]:
    return build_condition(
        context=context,
        materialization_id=materialization_id,
        condition_kind="clean",
        operator_id="CLEAN",
        variant="clean",
        event_fraction=None,
        mutation_seed=None,
        mutated_record_count=0,
        total_record_count=sum(value["records"] for value in context["clean_inputs"]),
        spec=runtime.file_identity(context["spec_path"]),
        inputs=context["clean_inputs"],
        target={
            "case_ids": [],
            "read_ids": [],
            "region_ids": [],
            "details": {},
        },
        mutation_details={},
        expected_seqspec_check="pass",
        observed_seqspec_check="pass",
        expected={
            "process_outcome": "success",
            "checks": [],
            "assessment_codes": [],
            "metric_names": [],
            "localization": [],
        },
    )


def generate_operator_conditions(
    *,
    record: dict[str, Any],
    operator: dict[str, Any],
    protocol: dict[str, Any],
    context: dict[str, Any],
    materialization_id: str,
    temporary_root: Path,
    output_root: Path,
    seqspec_bin: Path,
    timeout_seconds: int,
) -> dict[str, list[dict[str, Any]]]:
    if not record["applicable"]:
        return {
            "conditions": [],
            "skips": [
                skip_row(
                    context,
                    materialization_id=materialization_id,
                    operator_id=record["operator_id"],
                    variant="",
                    event_fraction=None,
                    mutation_seed=None,
                    reason_code="operator_inapplicable",
                    reason=record["reason"],
                )
            ],
        }
    if operator["stochastic"]:
        return generate_read_conditions(
            record=record,
            operator=operator,
            protocol=protocol,
            context=context,
            materialization_id=materialization_id,
            temporary_root=temporary_root,
            output_root=output_root,
        )
    conditions = []
    for variant in record["variants"]:
        conditions.append(
            generate_deterministic_condition(
                record=record,
                operator=operator,
                variant=variant,
                context=context,
                materialization_id=materialization_id,
                temporary_root=temporary_root,
                output_root=output_root,
                seqspec_bin=seqspec_bin,
                timeout_seconds=timeout_seconds,
            )
        )
    return {"conditions": conditions, "skips": []}


def generate_deterministic_condition(
    *,
    record: dict[str, Any],
    operator: dict[str, Any],
    variant: str,
    context: dict[str, Any],
    materialization_id: str,
    temporary_root: Path,
    output_root: Path,
    seqspec_bin: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    operator_id = record["operator_id"]
    relative_root = condition_directory(
        context["configuration_accession"], operator_id, variant, None, None
    )
    temporary_condition_root = temporary_root / relative_root
    final_condition_root = output_root / relative_root
    inputs = copy.deepcopy(context["clean_inputs"])
    spec_identity = runtime.file_identity(context["spec_path"])
    mutation_details = {}
    observed_check = "not_applicable"
    if operator_id in {"S01", "S02", "S03"}:
        inputs, mutation_details = mutate_input_contract(
            operator_id,
            variant,
            record["target"],
            inputs,
            temporary_condition_root=temporary_condition_root,
            final_condition_root=final_condition_root,
        )
    else:
        target_case_id = record["target"]["case_ids"][0]
        lengths = mutation.observed_lengths(context["records"][target_case_id])
        spec, mutation_details = mutation.mutate_spec(
            context["spec"],
            operator_id=operator_id,
            variant=variant,
            target=record["target"],
            observed_lengths=lengths,
        )
        replacements = {}
        omitted = set()
        if operator_id == "S09":
            entries = onlist_entries(context, record["target"])
            region_id = record["target"]["region_ids"][0]
            length = record["target"]["regions"][0]["max_len"]
            replacement = mutation.first_offlist_sequence(length, entries)
            replacements[region_id] = replacement + "\n"
            mutation_details["replacement_sequence"] = replacement
            mutation_details["source_onlist_entries"] = len(entries)
        elif operator_id == "S10":
            omitted.add(record["target"]["region_ids"][0])
        resources = mutation.bundle_local_resources(
            spec,
            source_spec_path=context["spec_path"],
            output_root=temporary_condition_root,
            omit_onlist_region_ids=omitted,
            replacement_onlists=replacements,
        )
        temporary_spec_path = temporary_condition_root / "spec.json"
        final_spec_path = final_condition_root / "spec.json"
        spec_identity = final_file_identity(
            output_root,
            temporary_root,
            temporary_spec_path,
            writer=lambda: mutation.write_spec(temporary_spec_path, spec),
        )
        observed_check, check_details = seqspec_check(
            seqspec_bin,
            temporary_spec_path,
            timeout_seconds=timeout_seconds,
        )
        expected_check = operator["expected_seqspec_check"]
        if observed_check != expected_check:
            raise ValueError(
                f"{operator_id}/{variant}: expected seqspec check {expected_check}, "
                f"observed {observed_check}: "
                f"{check_details['stderr'] or check_details['stdout']}"
            )
        spec_identity["path"] = str(final_spec_path)
        mutation_details["resources"] = [
            relocate_identity(value, temporary_root, output_root) for value in resources
        ]
        mutation_details["seqspec_check"] = check_details
    return build_condition(
        context=context,
        materialization_id=materialization_id,
        condition_kind="deterministic",
        operator_id=operator_id,
        variant=variant,
        event_fraction=None,
        mutation_seed=None,
        mutated_record_count=0,
        total_record_count=sum(value["records"] for value in inputs),
        spec=spec_identity,
        inputs=inputs,
        target=record["target"],
        mutation_details=mutation_details,
        expected_seqspec_check=operator["expected_seqspec_check"],
        observed_seqspec_check=observed_check,
        expected=operator["target"],
    )


def mutate_input_contract(
    operator_id: str,
    variant: str,
    target: dict[str, Any],
    inputs: list[dict[str, Any]],
    *,
    temporary_condition_root: Path,
    final_condition_root: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    by_case = {value["case_id"]: value for value in inputs}
    case_ids = target["case_ids"]
    if operator_id == "S01":
        case_id = case_ids[0]
        source = Path(by_case[case_id]["path"])
        basename = target["details"][
            "unexpected_basename"
            if variant == "unexpected_name"
            else "ambiguous_basename"
        ]
        if not basename:
            raise ValueError(f"S01 variant {variant} has no target basename")
        destination = temporary_condition_root / "inputs" / basename
        replacement = copied_input(
            by_case[case_id],
            source=source,
            destination=destination,
            final_path=final_condition_root / "inputs" / basename,
        )
        return replace_inputs(inputs, {case_id: replacement}), {"basename": basename}
    if operator_id == "S02":
        if len(case_ids) != 2:
            raise ValueError("S02 requires two target inputs")
        replacements = {}
        for target_case, donor_case in zip(case_ids, reversed(case_ids), strict=True):
            target_input = by_case[target_case]
            donor_input = by_case[donor_case]
            basename = Path(target_input["path"]).name
            replacements[target_case] = copied_input(
                target_input,
                source=Path(donor_input["path"]),
                destination=temporary_condition_root / "inputs" / basename,
                final_path=final_condition_root / "inputs" / basename,
            )
            replacements[target_case]["content_source_case_id"] = donor_case
        return replace_inputs(inputs, replacements), {"swapped_case_ids": case_ids}
    if operator_id == "S03":
        omitted = case_ids[0]
        result = [value for value in inputs if value["case_id"] != omitted]
        if len(result) == len(inputs):
            raise ValueError("S03 target input was not present")
        return result, {"omitted_case_id": omitted}
    raise ValueError(f"unsupported input operator: {operator_id}")


def generate_read_conditions(
    *,
    record: dict[str, Any],
    operator: dict[str, Any],
    protocol: dict[str, Any],
    context: dict[str, Any],
    materialization_id: str,
    temporary_root: Path,
    output_root: Path,
) -> dict[str, list[dict[str, Any]]]:
    operator_id = record["operator_id"]
    case_id = record["target"]["case_ids"][0]
    source_records = context["records"][case_id]
    variants = list(record["variants"])
    if operator_id == "D04" and record["target"]["details"].get("conditional_variant"):
        variants.append("remove_by_substitution")
    entries = (
        onlist_entries(context, record["target"]) if operator_id == "D03" else None
    )
    conditions = []
    skips = []
    for variant in variants:
        mutate_record, eligible_record = mutation.mutation_for_read_operator(
            operator_id,
            variant,
            record["target"],
            onlist_entries=entries,
        )
        eligible = {
            index
            for index, value in enumerate(source_records)
            if eligible_record(value[1].rstrip(b"\r\n"))
        }
        for event_fraction in protocol["stochastic_conditions"]["event_fractions"]:
            target_count = mutation.target_record_count(
                len(source_records), event_fraction
            )
            for seed in protocol["stochastic_conditions"]["seeds"]:
                if target_count == 0:
                    skips.append(
                        skip_row(
                            context,
                            materialization_id=materialization_id,
                            operator_id=operator_id,
                            variant=variant,
                            event_fraction=event_fraction,
                            mutation_seed=seed,
                            reason_code="event_fraction_unrepresentable",
                            reason=(
                                f"fraction {event_fraction} does not produce an integer "
                                f"event count in {len(source_records)} records"
                            ),
                        )
                    )
                    continue
                if len(eligible) < target_count:
                    skips.append(
                        skip_row(
                            context,
                            materialization_id=materialization_id,
                            operator_id=operator_id,
                            variant=variant,
                            event_fraction=event_fraction,
                            mutation_seed=seed,
                            reason_code="insufficient_eligible_records",
                            reason=(
                                f"{len(eligible)} eligible records for {target_count} "
                                "requested mutations"
                            ),
                        )
                    )
                    continue
                selected = mutation.select_record_indices(
                    source_records,
                    event_fraction=event_fraction,
                    seed=seed,
                    eligible=eligible,
                )
                changed = mutation.mutate_fastq_records(
                    source_records, selected, mutate_record
                )
                conditions.append(
                    write_read_condition(
                        record=record,
                        operator=operator,
                        variant=variant,
                        event_fraction=event_fraction,
                        seed=seed,
                        selected=selected,
                        eligible_count=len(eligible),
                        changed=changed,
                        context=context,
                        materialization_id=materialization_id,
                        temporary_root=temporary_root,
                        output_root=output_root,
                    )
                )
    return {"conditions": conditions, "skips": skips}


def write_read_condition(
    *,
    record: dict[str, Any],
    operator: dict[str, Any],
    variant: str,
    event_fraction: float,
    seed: int,
    selected: list[int],
    eligible_count: int,
    changed: list[mutation.FastqRecord],
    context: dict[str, Any],
    materialization_id: str,
    temporary_root: Path,
    output_root: Path,
) -> dict[str, Any]:
    operator_id = record["operator_id"]
    case_id = record["target"]["case_ids"][0]
    relative_root = condition_directory(
        context["configuration_accession"],
        operator_id,
        variant,
        event_fraction,
        seed,
    )
    temporary_condition_root = temporary_root / relative_root
    final_condition_root = output_root / relative_root
    original_input = next(
        value for value in context["clean_inputs"] if value["case_id"] == case_id
    )
    basename = Path(original_input["path"]).name
    temporary_fastq = temporary_condition_root / "inputs" / basename
    final_fastq = final_condition_root / "inputs" / basename
    mutation.write_fastq(temporary_fastq, changed)
    replacement = {
        **original_input,
        **relocated_file_identity(temporary_fastq, final_fastq),
        "content_source_case_id": case_id,
    }
    inputs = replace_inputs(context["clean_inputs"], {case_id: replacement})
    selected_ids = [
        mutation.fastq_record_id(context["records"][case_id][index])
        for index in selected
    ]
    selection_path = temporary_condition_root / "selected_record_ids.txt.gz"
    write_gzip_lines(selection_path, selected_ids)
    selection_identity = relocated_file_identity(
        selection_path, final_condition_root / selection_path.name
    )
    lengths_before = mutation.observed_lengths(context["records"][case_id])
    lengths_after = mutation.observed_lengths(changed)
    if operator_id != "D01" and lengths_before != lengths_after:
        raise ValueError(f"{operator_id} changed read lengths")
    return build_condition(
        context=context,
        materialization_id=materialization_id,
        condition_kind="stochastic",
        operator_id=operator_id,
        variant=variant,
        event_fraction=event_fraction,
        mutation_seed=seed,
        mutated_record_count=len(selected),
        total_record_count=len(changed),
        spec=runtime.file_identity(context["spec_path"]),
        inputs=inputs,
        target=record["target"],
        mutation_details={
            "selection_method": "hash_ranked_record_ids",
            "selected_record_ids": selection_identity,
            "selected_record_ids_sha256": runtime.sha256_json(selected_ids),
            "eligible_records": eligible_count,
            "observed_event_fraction": len(selected) / len(changed),
        },
        expected_seqspec_check=operator["expected_seqspec_check"],
        observed_seqspec_check="not_applicable",
        expected=operator["target"],
    )


def onlist_entries(context: dict[str, Any], target: dict[str, Any]) -> set[str]:
    region = target["regions"][0]
    resource = region.get("resource", {})
    key = (
        str(resource.get("urltype", "")),
        str(resource.get("locator", "")),
        int(region["max_len"]),
    )
    if key not in context["onlists"]:
        locator = key[1]
        source = (
            str((context["spec_path"].parent / locator).resolve())
            if key[0] == "local"
            else locator
        )
        entries = mutation.read_onlist_entries(source, key[2])
        if not entries:
            raise ValueError(
                f"onlist parser found no length-{key[2]} DNA entries in {locator}"
            )
        context["onlists"][key] = entries
    return context["onlists"][key]


def build_condition(
    *,
    context: dict[str, Any],
    materialization_id: str,
    condition_kind: str,
    operator_id: str,
    variant: str,
    event_fraction: float | None,
    mutation_seed: int | None,
    mutated_record_count: int,
    total_record_count: int,
    spec: dict[str, Any],
    inputs: list[dict[str, Any]],
    target: dict[str, Any],
    mutation_details: dict[str, Any],
    expected_seqspec_check: str,
    observed_seqspec_check: str,
    expected: dict[str, Any],
) -> dict[str, Any]:
    stable = {
        "schema_version": SCHEMA_VERSION,
        "materialization_id": materialization_id,
        "selection_id": context["selection_id"],
        "configuration_accession": context["configuration_accession"],
        "family_id": context["family_id"],
        "modality": context["modality"],
        "condition_kind": condition_kind,
        "operator_id": operator_id,
        "variant": variant,
        "event_fraction": event_fraction,
        "mutation_seed": mutation_seed,
        "mutated_record_count": mutated_record_count,
        "total_record_count": total_record_count,
        "spec": functional_file_identity(spec),
        "inputs": [functional_input_identity(value) for value in inputs],
        "target": target,
        "mutation": functional_mutation_details(mutation_details),
        "expected_seqspec_check": expected_seqspec_check,
        "observed_seqspec_check": observed_seqspec_check,
        "expected": expected,
    }
    return {
        **stable,
        "condition_id": runtime.sha256_json(stable)[:16],
        "spec": spec,
        "inputs": inputs,
        "mutation": mutation_details,
    }


def skip_row(
    context: dict[str, Any],
    *,
    materialization_id: str,
    operator_id: str,
    variant: str,
    event_fraction: float | None,
    mutation_seed: int | None,
    reason_code: str,
    reason: str,
) -> dict[str, Any]:
    return {
        "materialization_id": materialization_id,
        "selection_id": context["selection_id"],
        "configuration_accession": context["configuration_accession"],
        "family_id": context["family_id"],
        "modality": context["modality"],
        "operator_id": operator_id,
        "variant": variant,
        "event_fraction": event_fraction,
        "mutation_seed": mutation_seed,
        "reason_code": reason_code,
        "reason": reason,
    }


def validate_materialization(
    *,
    materialization_id: str,
    inventory_records: list[dict[str, Any]],
    protocol: dict[str, Any],
    conditions: list[dict[str, Any]],
    skips: list[dict[str, Any]],
    temporary_root: Path,
    output_root: Path,
) -> dict[str, Any]:
    errors = []
    condition_ids = [value["condition_id"] for value in conditions]
    if len(condition_ids) != len(set(condition_ids)):
        errors.append("condition identifiers are not unique")
    configurations = {value["configuration_accession"] for value in inventory_records}
    clean_counts = Counter(
        value["configuration_accession"]
        for value in conditions
        if value["condition_kind"] == "clean"
    )
    if clean_counts != Counter({value: 1 for value in configurations}):
        errors.append("clean condition count does not reconcile")
    declared_variants = {
        operator["operator_id"]: set(operator["variants"])
        for operator in protocol["operators"]
    }
    for condition in conditions:
        if condition["materialization_id"] != materialization_id:
            errors.append(f"{condition['condition_id']}: materialization id differs")
        if (
            condition["operator_id"] != "CLEAN"
            and condition["variant"] not in declared_variants[condition["operator_id"]]
        ):
            errors.append(f"{condition['condition_id']}: variant is undeclared")
        file_identities = [
            condition["spec"],
            *condition["inputs"],
            *nested_file_identities(condition["mutation"]),
        ]
        for value in file_identities:
            path = Path(value["path"])
            if path.is_relative_to(output_root):
                path = temporary_root / path.relative_to(output_root)
            if runtime.file_sha256(path) != value["sha256"]:
                errors.append(f"{condition['condition_id']}: output hash changed")
        if condition["expected_seqspec_check"] != condition["observed_seqspec_check"]:
            errors.append(f"{condition['condition_id']}: seqspec check outcome differs")
    covered_pairs = {
        (value["configuration_accession"], value["operator_id"])
        for value in conditions
        if value["operator_id"] != "CLEAN"
    } | {(value["configuration_accession"], value["operator_id"]) for value in skips}
    expected_pairs = {
        (value["configuration_accession"], value["operator_id"])
        for value in inventory_records
    }
    if covered_pairs != expected_pairs:
        errors.append("inventory configuration/operator coverage does not reconcile")
    return {
        "schema_version": SCHEMA_VERSION,
        "materialization_id": materialization_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "counts": {
            "configurations": len(configurations),
            "inventory_records": len(inventory_records),
            "conditions": len(conditions),
            "clean_conditions": sum(
                value["condition_kind"] == "clean" for value in conditions
            ),
            "deterministic_conditions": sum(
                value["condition_kind"] == "deterministic" for value in conditions
            ),
            "stochastic_conditions": sum(
                value["condition_kind"] == "stochastic" for value in conditions
            ),
            "skips": len(skips),
            "skips_by_reason": dict(
                sorted(Counter(value["reason_code"] for value in skips).items())
            ),
            "conditions_by_operator": dict(
                sorted(Counter(value["operator_id"] for value in conditions).items())
            ),
        },
    }


def public_condition_row(value: dict[str, Any]) -> dict[str, Any]:
    expected = value["expected"]
    target = value["target"]
    return {
        **value,
        "spec_path": value["spec"]["path"],
        "input_paths": ";".join(item["path"] for item in value["inputs"]),
        "target_case_ids": ";".join(target.get("case_ids", [])),
        "target_read_ids": ";".join(target.get("read_ids", [])),
        "target_region_ids": ";".join(target.get("region_ids", [])),
        "expected_process_outcome": expected.get("process_outcome", ""),
        "expected_checks": ";".join(expected.get("checks", [])),
        "expected_assessment_codes": ";".join(expected.get("assessment_codes", [])),
        "expected_metric_names": ";".join(expected.get("metric_names", [])),
        "condition_json": runtime.canonical_json(value),
    }


def copied_input(
    original: dict[str, Any],
    *,
    source: Path,
    destination: Path,
    final_path: Path,
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    return {**original, **relocated_file_identity(destination, final_path)}


def replace_inputs(
    inputs: list[dict[str, Any]], replacements: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    return [replacements.get(value["case_id"], value) for value in inputs]


def condition_directory(
    configuration_accession: str,
    operator_id: str,
    variant: str,
    event_fraction: float | None,
    seed: int | None,
) -> Path:
    parts = [
        "conditions",
        mutation.safe_name(configuration_accession),
        operator_id,
        mutation.safe_name(variant),
    ]
    if event_fraction is not None:
        parts.append(f"fraction-{str(event_fraction).replace('.', '_')}")
    if seed is not None:
        parts.append(f"seed-{seed:09d}")
    return Path(*parts)


def seqspec_check(
    seqspec_bin: Path, spec_path: Path, *, timeout_seconds: int
) -> tuple[str, dict[str, Any]]:
    completed = subprocess.run(
        [str(seqspec_bin), "check", str(spec_path)],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    return (
        "pass" if completed.returncode == 0 else "failure",
        {
            "exit_code": completed.returncode,
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
        },
    )


def tool_identities(
    seqspec_identity: dict[str, Any], yq_identity: dict[str, Any]
) -> dict[str, Any]:
    return {
        "seqspec": seqspec_identity,
        "yq": yq_identity,
        "materializer": runtime.script_identity(
            Path(__file__).resolve(), version=MATERIALIZER_VERSION
        ),
        "mutation_runtime": runtime.script_identity(
            Path(mutation.__file__).resolve(), version=mutation.RUNTIME_VERSION
        ),
        "paper_runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
        "fastq_runtime": runtime.script_identity(
            Path(mutation.sampler.__file__).resolve(),
            version=mutation.sampler.SAMPLER_VERSION,
        ),
    }


def require_executable(path: Path, timeout_seconds: int) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"executable does not exist: {path}")
    return runtime.executable_identity(path, timeout_seconds=min(timeout_seconds, 30))


def verified_path(identity: Any, label: str) -> Path:
    if not isinstance(identity, dict):
        raise ValueError(f"{label} identity is missing")
    path = Path(str(identity.get("path", ""))).resolve()
    if runtime.file_sha256(path) != identity.get("sha256"):
        raise ValueError(f"{label} hash changed")
    return path


def functional_file_identity(value: dict[str, Any]) -> dict[str, Any]:
    return {"sha256": value["sha256"], "size_bytes": value["size_bytes"]}


def functional_input_identity(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "case_id": value["case_id"],
        "fastq_accession": value["fastq_accession"],
        "read_id": value["read_id"],
        "sha256": value["sha256"],
        "size_bytes": value["size_bytes"],
        "records": value["records"],
        "content_source_case_id": value.get("content_source_case_id", value["case_id"]),
    }


def functional_mutation_details(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: functional_mutation_details(item)
            for key, item in value.items()
            if key not in {"path", "seqspec_check"}
        }
    if isinstance(value, list):
        return [functional_mutation_details(item) for item in value]
    return value


def nested_file_identities(value: Any) -> list[dict[str, Any]]:
    result = []
    if isinstance(value, dict):
        if {"path", "sha256", "size_bytes"}.issubset(value):
            result.append(value)
        for item in value.values():
            result.extend(nested_file_identities(item))
    elif isinstance(value, list):
        for item in value:
            result.extend(nested_file_identities(item))
    return result


def relocated_file_identity(path: Path, final_path: Path) -> dict[str, Any]:
    return {
        "path": str(final_path),
        "sha256": runtime.file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def relocate_identity(
    value: dict[str, Any], temporary_root: Path, output_root: Path
) -> dict[str, Any]:
    result = dict(value)
    if "path" in result:
        path = Path(result["path"]).resolve()
        temporary_root = temporary_root.resolve()
        output_root = output_root.resolve()
        if path.is_relative_to(temporary_root):
            result["path"] = str(output_root / path.relative_to(temporary_root))
    return result


def final_file_identity(
    output_root: Path,
    temporary_root: Path,
    path: Path,
    writer: Any | None = None,
) -> dict[str, Any]:
    if writer is not None:
        writer()
    return {
        "path": str(output_root / path.relative_to(temporary_root)),
        "sha256": runtime.file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def write_gzip_lines(path: Path, values: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as handle:
            for value in values:
                handle.write(value.encode() + b"\n")


def positive_int(value: Any, label: str) -> int:
    result = mutation.strict_int(value, label)
    if result <= 0:
        raise ValueError(f"{label} must be positive")
    return result


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
