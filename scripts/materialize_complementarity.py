#!/usr/bin/env python3
"""Materialize quality, composition, and schema controls for Experiment 3."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import build_perturbation_cases as inventory_builder
    import complementarity_runtime as controls
    import materialize_perturbations as base_materializer
    import paper_runtime as runtime
    import perturbation_runtime as perturb
    import run_perturbation_calibration as base_runner
except ModuleNotFoundError:
    from scripts import build_perturbation_cases as inventory_builder
    from scripts import complementarity_runtime as controls
    from scripts import materialize_perturbations as base_materializer
    from scripts import paper_runtime as runtime
    from scripts import perturbation_runtime as perturb
    from scripts import run_perturbation_calibration as base_runner


SCHEMA_VERSION = "0.1.0"
MATERIALIZER_VERSION = "0.1.0"
COHORT_SPLITS = ("calibration", "evaluation")
MATERIALIZATION_IDENTITY_FIELDS = (
    "schema_version",
    "cohort_split",
    "selection_id",
    "base_materialization_id",
    "base_materialization_sha256",
    "protocol_sha256",
    "quality_policy",
    "seqspec",
    "yq",
    "materializer",
    "control_runtime",
    "perturbation_runtime",
    "paper_runtime",
    "controls",
    "skips",
)
CONTROL_FIELDS = (
    "complementarity_id",
    "condition_id",
    "base_materialization_id",
    "selection_id",
    "configuration_accession",
    "family_id",
    "modality",
    "condition_kind",
    "operator_id",
    "variant",
    "severity_index",
    "event_fraction",
    "mutation_seed",
    "mutated_record_count",
    "total_record_count",
    "target_case_ids",
    "target_read_ids",
    "target_region_ids",
    "expected_seqspec_check",
    "observed_seqspec_check",
)
SKIP_FIELDS = (
    "configuration_accession",
    "operator_id",
    "variant",
    "reason_code",
    "reason",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Materialize frozen QC complementarity controls."
    )
    parser.add_argument("--base-materialization", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--seqspec-bin", required=True, type=Path)
    parser.add_argument("--yq-bin", required=True, type=Path)
    parser.add_argument("--cohort-split", required=True, choices=COHORT_SPLITS)
    parser.add_argument("--quality-policy", type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=120)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = materialize_complementarity(
            base_materialization_path=args.base_materialization.resolve(),
            protocol_path=args.protocol.resolve(),
            seqspec_bin=args.seqspec_bin.resolve(),
            yq_bin=args.yq_bin.resolve(),
            cohort_split=args.cohort_split,
            quality_policy_path=(
                args.quality_policy.resolve() if args.quality_policy else None
            ),
            output_root=args.output_root.resolve(),
            timeout_seconds=args.timeout_seconds,
        )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"materialize_complementarity: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def materialize_complementarity(
    *,
    base_materialization_path: Path,
    protocol_path: Path,
    seqspec_bin: Path,
    yq_bin: Path,
    cohort_split: str,
    output_root: Path,
    timeout_seconds: int,
    quality_policy_path: Path | None = None,
) -> dict[str, Any]:
    if timeout_seconds <= 0:
        raise ValueError("timeout must be positive")
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    if cohort_split not in COHORT_SPLITS:
        raise ValueError(f"unsupported cohort split: {cohort_split}")
    if (cohort_split == "evaluation") != (quality_policy_path is not None):
        raise ValueError(
            "evaluation requires one quality policy; calibration forbids it"
        )
    base, base_conditions = base_runner.load_materialization(base_materialization_path)
    protocol = runtime.load_json(protocol_path)
    validate_protocol(protocol)
    cases = load_cases(base, cohort_split=cohort_split)
    validate_base_source(base, cohort_split)
    policy = (
        load_quality_policy(quality_policy_path, protocol_path, protocol)
        if quality_policy_path
        else None
    )
    seqspec = base_materializer.require_executable(seqspec_bin, timeout_seconds)
    if runtime.functional_executable_identity(seqspec) != base["seqspec"]:
        raise ValueError("seqspec executable differs from base materialization")
    yq = base_materializer.require_executable(yq_bin, timeout_seconds)
    if runtime.functional_executable_identity(yq) != base["yq"]:
        raise ValueError("yq executable differs from base materialization")
    tools = tool_identities(seqspec, yq)
    clean_by_configuration = clean_conditions(base_conditions)
    configuration_conditions = conditions_by_configuration(base_conditions)

    output_root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_root.name}-", dir=output_root.parent
    ) as tmpdir:
        temporary_root = Path(tmpdir).resolve()
        control_rows = []
        skips = []
        for configuration, clean in sorted(clean_by_configuration.items()):
            case_map = {
                case["case_id"]: case
                for case in cases
                if case["configuration_accession"] == configuration
            }
            primary = primary_input(clean, case_map)
            records = perturb.read_fastq(Path(primary["path"]))
            quality = materialize_quality_controls(
                clean=clean,
                primary=primary,
                records=records,
                protocol=protocol,
                policy=policy,
                base_materialization_id=base["materialization_id"],
                temporary_root=temporary_root,
                output_root=output_root,
            )
            control_rows.extend(quality)
            composition = materialize_composition_controls(
                clean=clean,
                primary=primary,
                records=records,
                base_conditions=configuration_conditions[configuration],
                protocol=protocol,
                yq_bin=yq_bin,
                timeout_seconds=timeout_seconds,
                base_materialization_id=base["materialization_id"],
                temporary_root=temporary_root,
                output_root=output_root,
            )
            control_rows.extend(composition["conditions"])
            skips.extend(composition["skips"])
            control_rows.extend(
                materialize_schema_controls(
                    clean=clean,
                    protocol=protocol,
                    seqspec_bin=seqspec_bin,
                    yq_bin=yq_bin,
                    timeout_seconds=timeout_seconds,
                    base_materialization_id=base["materialization_id"],
                    temporary_root=temporary_root,
                    output_root=output_root,
                )
            )

        stable_identity = {
            "schema_version": SCHEMA_VERSION,
            "cohort_split": cohort_split,
            "selection_id": base["selection_id"],
            "base_materialization_id": base["materialization_id"],
            "base_materialization_sha256": runtime.file_sha256(
                base_materialization_path
            ),
            "protocol_sha256": runtime.file_sha256(protocol_path),
            "quality_policy": (
                {
                    "policy_id": policy["policy_id"],
                    "sha256": runtime.file_sha256(quality_policy_path),
                }
                if policy
                else None
            ),
            "seqspec": runtime.functional_executable_identity(seqspec),
            "yq": runtime.functional_executable_identity(yq),
            "materializer": runtime.functional_script_identity(tools["materializer"]),
            "control_runtime": runtime.functional_script_identity(
                tools["control_runtime"]
            ),
            "perturbation_runtime": runtime.functional_script_identity(
                tools["perturbation_runtime"]
            ),
            "paper_runtime": runtime.functional_script_identity(tools["paper_runtime"]),
            "controls": [control_identity(value) for value in control_rows],
            "skips": skips,
        }
        complementarity_id = runtime.sha256_json(stable_identity)[:16]
        controls_with_id = [
            {"complementarity_id": complementarity_id, **value}
            for value in control_rows
        ]
        validation = validate_outputs(
            complementarity_id=complementarity_id,
            cohort_split=cohort_split,
            configurations=set(clean_by_configuration),
            controls=controls_with_id,
            skips=skips,
            protocol=protocol,
            policy=policy,
        )
        if not validation["valid"]:
            raise ValueError(
                "complementarity materialization failed: "
                + "; ".join(validation["errors"] or ["unknown validation error"])
            )

        controls_path = temporary_root / "manifests" / "controls.json"
        table_path = temporary_root / "tables" / "controls.csv"
        skips_path = temporary_root / "tables" / "skips.csv"
        validation_path = temporary_root / "validation" / "materialization.json"
        manifest_path = temporary_root / "manifests" / "materialization.json"
        runtime.write_json(
            controls_path,
            {
                "schema_version": SCHEMA_VERSION,
                "complementarity_id": complementarity_id,
                "controls": controls_with_id,
            },
        )
        runtime.write_csv(
            table_path,
            [public_control_row(value) for value in controls_with_id],
            list(CONTROL_FIELDS),
        )
        runtime.write_csv(skips_path, skips, list(SKIP_FIELDS))
        runtime.write_json(validation_path, validation)
        final_manifest_path = output_root / manifest_path.relative_to(temporary_root)
        manifest = {
            **stable_identity,
            "complementarity_id": complementarity_id,
            "created_at": utc_now(),
            "valid": True,
            "inputs": {
                "base_materialization": runtime.file_identity(
                    base_materialization_path
                ),
                "protocol": runtime.file_identity(protocol_path),
                "quality_policy": (
                    runtime.file_identity(quality_policy_path)
                    if quality_policy_path
                    else None
                ),
            },
            "tools": tools,
            "counts": validation["counts"],
            "outputs": {
                "controls": base_materializer.final_file_identity(
                    output_root, temporary_root, controls_path
                ),
                "table": base_materializer.final_file_identity(
                    output_root, temporary_root, table_path
                ),
                "skips": base_materializer.final_file_identity(
                    output_root, temporary_root, skips_path
                ),
                "validation": base_materializer.final_file_identity(
                    output_root, temporary_root, validation_path
                ),
            },
            "manifest_path": str(final_manifest_path),
        }
        runtime.write_json(manifest_path, manifest)
        os.replace(temporary_root, output_root)
    return manifest


def materialize_quality_controls(
    *,
    clean: dict[str, Any],
    primary: dict[str, Any],
    records: list[perturb.FastqRecord],
    protocol: dict[str, Any],
    policy: dict[str, Any] | None,
    base_materialization_id: str,
    temporary_root: Path,
    output_root: Path,
) -> list[dict[str, Any]]:
    policy_map = (
        {value["operator_id"]: value for value in policy["entries"]} if policy else {}
    )
    result = []
    for operator in protocol["quality_controls"]:
        candidates = operator["severity_candidates"]
        indexed = (
            [
                (
                    policy_map[operator["operator_id"]]["severity_index"],
                    policy_map[operator["operator_id"]]["severity"],
                )
            ]
            if policy
            else list(enumerate(candidates))
        )
        for severity_index, severity in indexed:
            mutated, details = controls.mutate_quality_records(
                records,
                operator_id=operator["operator_id"],
                severity=severity,
                seed=protocol["quality_seed"],
            )
            relative = control_directory(
                clean["configuration_accession"],
                operator["operator_id"],
                operator["variant"],
                severity_index,
            )
            temporary_path = (
                temporary_root / relative / "inputs" / Path(primary["path"]).name
            )
            identity = base_materializer.final_file_identity(
                output_root,
                temporary_root,
                temporary_path,
                writer=lambda path=temporary_path, values=mutated: controls.write_fastq(
                    path, values
                ),
            )
            inputs = replace_input(clean["inputs"], primary["case_id"], identity)
            event_fraction = severity.get(
                "event_fraction", severity.get("low_fraction")
            )
            result.append(
                build_control(
                    clean=clean,
                    base_materialization_id=base_materialization_id,
                    condition_kind="quality",
                    operator_id=operator["operator_id"],
                    variant=operator["variant"],
                    severity_index=severity_index,
                    severity=severity,
                    event_fraction=event_fraction,
                    mutation_seed=protocol["quality_seed"],
                    mutated_record_count=details["records_changed"],
                    spec=clean["spec"],
                    inputs=inputs,
                    target={
                        "case_ids": [primary["case_id"]],
                        "read_ids": [primary["read_id"]],
                        "region_ids": [],
                    },
                    mutation=details,
                    expected_seqspec_check="pass",
                    observed_seqspec_check="pass",
                    expected={
                        "seqcheck_process": "success",
                        "fastqc_modules": operator["fastqc_modules"],
                        "sequence_invariant": True,
                    },
                )
            )
    return result


def materialize_composition_controls(
    *,
    clean: dict[str, Any],
    primary: dict[str, Any],
    records: list[perturb.FastqRecord],
    base_conditions: list[dict[str, Any]],
    protocol: dict[str, Any],
    yq_bin: Path,
    timeout_seconds: int,
    base_materialization_id: str,
    temporary_root: Path,
    output_root: Path,
) -> dict[str, list[dict[str, Any]]]:
    spec = perturb.load_seqspec_json(
        Path(clean["spec"]["path"]), yq_bin=yq_bin, timeout_seconds=timeout_seconds
    )
    measurement = measurement_target(
        spec, base_conditions=base_conditions, primary=primary
    )
    conditions = []
    skips = []
    for operator in protocol["composition_controls"]:
        operator_id = operator["operator_id"]
        span = None
        target = {
            "case_ids": [primary["case_id"]],
            "read_ids": [primary["read_id"]],
            "region_ids": [],
        }
        if operator["target"] == "projected_measurement_payload":
            minimum = len(operator.get("motif", "")) if operator_id == "C01" else 4
            if (
                measurement is None
                or measurement["stop"] - measurement["start"] < minimum
            ):
                skips.append(
                    skip_row(
                        clean,
                        operator,
                        "no_eligible_measurement_payload",
                        f"requires a projected measurement payload of at least {minimum} bases",
                    )
                )
                continue
            span = (measurement["start"], measurement["stop"])
            target["region_ids"] = [measurement["region_id"]]
        try:
            mutated, details = controls.mutate_composition_records(
                records,
                operator_id=operator_id,
                event_fraction=operator["event_fraction"],
                seed=operator["seed"],
                target_span=span,
                motif=operator.get("motif"),
            )
        except ValueError as error:
            skips.append(
                skip_row(
                    clean,
                    operator,
                    "insufficient_eligible_records",
                    str(error),
                )
            )
            continue
        relative = control_directory(
            clean["configuration_accession"], operator_id, operator["variant"], None
        )
        temporary_path = (
            temporary_root / relative / "inputs" / Path(primary["path"]).name
        )
        identity = base_materializer.final_file_identity(
            output_root,
            temporary_root,
            temporary_path,
            writer=lambda path=temporary_path, values=mutated: controls.write_fastq(
                path, values
            ),
        )
        inputs = replace_input(clean["inputs"], primary["case_id"], identity)
        conditions.append(
            build_control(
                clean=clean,
                base_materialization_id=base_materialization_id,
                condition_kind="composition",
                operator_id=operator_id,
                variant=operator["variant"],
                severity_index=None,
                severity=None,
                event_fraction=operator["event_fraction"],
                mutation_seed=operator["seed"],
                mutated_record_count=details["records_selected"],
                spec=clean["spec"],
                inputs=inputs,
                target=target,
                mutation=details,
                expected_seqspec_check="pass",
                observed_seqspec_check="pass",
                expected={
                    "seqcheck_process": "success",
                    "fastqc_modules": operator["fastqc_modules"],
                    "sequence_invariant": False,
                },
            )
        )
    return {"conditions": conditions, "skips": skips}


def materialize_schema_controls(
    *,
    clean: dict[str, Any],
    protocol: dict[str, Any],
    seqspec_bin: Path,
    yq_bin: Path,
    timeout_seconds: int,
    base_materialization_id: str,
    temporary_root: Path,
    output_root: Path,
) -> list[dict[str, Any]]:
    original = perturb.load_seqspec_json(
        Path(clean["spec"]["path"]), yq_bin=yq_bin, timeout_seconds=timeout_seconds
    )
    result = []
    for operator in protocol["schema_invalid_controls"]:
        spec, details = controls.mutate_invalid_spec(
            original, variant=operator["variant"]
        )
        relative = control_directory(
            clean["configuration_accession"],
            operator["operator_id"],
            operator["variant"],
            None,
        )
        temporary_condition_root = temporary_root / relative
        resources = perturb.bundle_local_resources(
            spec,
            source_spec_path=Path(clean["spec"]["path"]),
            output_root=temporary_condition_root,
        )
        temporary_spec_path = temporary_condition_root / "spec.json"
        spec_identity = base_materializer.final_file_identity(
            output_root,
            temporary_root,
            temporary_spec_path,
            writer=lambda path=temporary_spec_path, value=spec: perturb.write_spec(
                path, value
            ),
        )
        observed, check = base_materializer.seqspec_check(
            seqspec_bin, temporary_spec_path, timeout_seconds=timeout_seconds
        )
        if observed != operator["expected_seqspec_check"]:
            raise ValueError(
                f"{operator['operator_id']}: expected seqspec check "
                f"{operator['expected_seqspec_check']}, observed {observed}"
            )
        details["resources"] = [
            base_materializer.relocate_identity(value, temporary_root, output_root)
            for value in resources
        ]
        details["seqspec_check"] = check
        result.append(
            build_control(
                clean=clean,
                base_materialization_id=base_materialization_id,
                condition_kind="schema_invalid",
                operator_id=operator["operator_id"],
                variant=operator["variant"],
                severity_index=None,
                severity=None,
                event_fraction=None,
                mutation_seed=None,
                mutated_record_count=0,
                spec=spec_identity,
                inputs=clean["inputs"],
                target={"case_ids": [], "read_ids": [], "region_ids": []},
                mutation=details,
                expected_seqspec_check=operator["expected_seqspec_check"],
                observed_seqspec_check=observed,
                expected={
                    "seqcheck_process": "failure",
                    "fastqc_modules": [],
                    "sequence_invariant": True,
                },
            )
        )
    return result


def build_control(
    *,
    clean: dict[str, Any],
    base_materialization_id: str,
    condition_kind: str,
    operator_id: str,
    variant: str,
    severity_index: int | None,
    severity: dict[str, Any] | None,
    event_fraction: float | None,
    mutation_seed: int | None,
    mutated_record_count: int,
    spec: dict[str, Any],
    inputs: list[dict[str, Any]],
    target: dict[str, Any],
    mutation: dict[str, Any],
    expected_seqspec_check: str,
    observed_seqspec_check: str,
    expected: dict[str, Any],
) -> dict[str, Any]:
    stable = {
        "schema_version": SCHEMA_VERSION,
        "base_materialization_id": base_materialization_id,
        "selection_id": clean["selection_id"],
        "configuration_accession": clean["configuration_accession"],
        "family_id": clean["family_id"],
        "modality": clean["modality"],
        "condition_kind": condition_kind,
        "operator_id": operator_id,
        "variant": variant,
        "severity_index": severity_index,
        "severity": severity,
        "event_fraction": event_fraction,
        "mutation_seed": mutation_seed,
        "mutated_record_count": mutated_record_count,
        "total_record_count": sum(value["records"] for value in inputs),
        "spec": base_materializer.functional_file_identity(spec),
        "inputs": [
            base_materializer.functional_input_identity(value) for value in inputs
        ],
        "target": target,
        "mutation": base_materializer.functional_mutation_details(mutation),
        "expected_seqspec_check": expected_seqspec_check,
        "observed_seqspec_check": observed_seqspec_check,
        "expected": expected,
    }
    return {
        **stable,
        "condition_id": runtime.sha256_json(stable)[:16],
        "spec": spec,
        "inputs": inputs,
        "mutation": mutation,
    }


def control_identity(value: dict[str, Any]) -> dict[str, Any]:
    stable = {
        key: item
        for key, item in value.items()
        if key
        not in {"complementarity_id", "condition_id", "spec", "inputs", "mutation"}
    }
    stable["spec"] = base_materializer.functional_file_identity(value["spec"])
    stable["inputs"] = [
        base_materializer.functional_input_identity(item) for item in value["inputs"]
    ]
    stable["mutation"] = base_materializer.functional_mutation_details(
        value["mutation"]
    )
    return stable


def load_cases(base: dict[str, Any], *, cohort_split: str) -> list[dict[str, Any]]:
    inventory_path = verified_path(
        base.get("inputs", {}).get("perturbation_inventory", {}),
        "base perturbation inventory",
    )
    inventory = runtime.load_json(inventory_path)
    selection_path = verified_path(
        inventory.get("inputs", {}).get("case_selection_manifest", {}),
        "base case selection",
    )
    selection, cases = inventory_builder.load_case_selection(selection_path)
    if selection["selection_id"] != base["selection_id"]:
        raise ValueError("base case selection differs from materialization")
    if selection.get("cohort_split") != cohort_split:
        raise ValueError("base case selection uses a different cohort split")
    return cases


def validate_base_source(base: dict[str, Any], cohort_split: str) -> None:
    expected = (
        "sampling_study" if cohort_split == "calibration" else "policy_sample_bundle"
    )
    if base.get("sample_source", {}).get("kind") != expected:
        raise ValueError(f"{cohort_split} requires base sample source {expected}")


def clean_conditions(
    conditions: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    grouped = defaultdict(list)
    for condition in conditions:
        if condition["condition_kind"] == "clean":
            grouped[condition["configuration_accession"]].append(condition)
    if not grouped or any(len(values) != 1 for values in grouped.values()):
        raise ValueError(
            "base materialization must have one clean condition per configuration"
        )
    return {key: values[0] for key, values in grouped.items()}


def conditions_by_configuration(
    conditions: list[dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    grouped = defaultdict(list)
    for condition in conditions:
        grouped[condition["configuration_accession"]].append(condition)
    return grouped


def primary_input(
    clean: dict[str, Any], case_map: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    matches = [
        value
        for value in clean["inputs"]
        if case_map.get(value["case_id"], {}).get("selection_role") == "primary"
    ]
    if len(matches) != 1:
        raise ValueError(
            f"{clean['configuration_accession']}: no unique primary FASTQ input"
        )
    return matches[0]


def measurement_target(
    spec: dict[str, Any],
    *,
    base_conditions: list[dict[str, Any]],
    primary: dict[str, Any],
) -> dict[str, Any] | None:
    candidates = {}
    for condition in base_conditions:
        target = condition.get("target", {})
        if primary["case_id"] not in target.get("case_ids", []):
            continue
        for region in target.get("regions", []):
            region_id = str(region.get("region_id", ""))
            if not region_id:
                continue
            declared = perturb.find_region(spec, region_id)
            terms = declared.get("region_type", [])
            if isinstance(terms, str):
                terms = [terms]
            if (
                declared.get("sequence_type") == "random"
                and any(str(value).startswith("RGN:measure:") for value in terms)
                and isinstance(region.get("start"), int)
                and isinstance(region.get("stop"), int)
            ):
                candidates[region_id] = {
                    "region_id": region_id,
                    "start": region["start"],
                    "stop": region["stop"],
                }
    return max(
        candidates.values(),
        key=lambda value: (value["stop"] - value["start"], value["region_id"]),
        default=None,
    )


def replace_input(
    inputs: list[dict[str, Any]], case_id: str, identity: dict[str, Any]
) -> list[dict[str, Any]]:
    result = []
    for value in inputs:
        result.append({**value, **identity} if value["case_id"] == case_id else value)
    return result


def load_quality_policy(
    path: Path, protocol_path: Path, protocol: dict[str, Any]
) -> dict[str, Any]:
    policy = runtime.load_json(path)
    stable = {
        key: value
        for key, value in policy.items()
        if key not in {"policy_id", "created_at"}
    }
    if policy.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("quality policy schema is unsupported")
    if runtime.sha256_json(stable)[:16] != policy.get("policy_id"):
        raise ValueError("quality policy id is not content-addressed")
    if policy.get("frozen") is not True:
        raise ValueError("quality policy is not frozen")
    if policy.get("protocol_sha256") != runtime.file_sha256(protocol_path):
        raise ValueError("quality policy uses a different protocol")
    entries = policy.get("entries")
    expected = {value["operator_id"]: value for value in protocol["quality_controls"]}
    entry_ids = (
        [value.get("operator_id") for value in entries]
        if isinstance(entries, list)
        else []
    )
    if (
        not isinstance(entries, list)
        or len(entry_ids) != len(set(entry_ids))
        or set(entry_ids) != set(expected)
    ):
        raise ValueError("quality policy entries do not reconcile")
    for entry in entries:
        index = entry.get("severity_index")
        if isinstance(index, bool) or not isinstance(index, int):
            raise ValueError("quality policy severity index is invalid")
        candidates = expected[entry["operator_id"]]["severity_candidates"]
        if (
            not 0 <= index < len(candidates)
            or entry.get("severity") != candidates[index]
        ):
            raise ValueError("quality policy severity differs from protocol")
    return policy


def validate_protocol(protocol: dict[str, Any]) -> None:
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("complementarity protocol schema is unsupported")
    quality = protocol.get("quality_controls")
    composition = protocol.get("composition_controls")
    invalid = protocol.get("schema_invalid_controls")
    if not all(
        isinstance(value, list) and value for value in (quality, composition, invalid)
    ):
        raise ValueError("complementarity control lists must be nonempty")
    expected_ids = {
        "quality": {"Q01", "Q02", "Q03", "Q04"},
        "composition": {"C01", "C02", "C03"},
        "invalid": {"X01", "X02"},
    }
    for label, values in (
        ("quality", quality),
        ("composition", composition),
        ("invalid", invalid),
    ):
        ids = [required_string(value, "operator_id") for value in values]
        if set(ids) != expected_ids[label] or len(ids) != len(set(ids)):
            raise ValueError(f"{label} control identifiers are invalid")
        for value in values:
            required_string(value, "variant")
    for value in quality:
        candidates = value.get("severity_candidates")
        if not isinstance(candidates, list) or not candidates:
            raise ValueError(f"{value['operator_id']}: severity candidates are empty")
        validate_modules(value)
    seed = protocol.get("quality_seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("quality seed is invalid")
    for value in composition:
        controls.unit_fraction(value.get("event_fraction"), "event fraction")
        seed = value.get("seed")
        if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
            raise ValueError(f"{value['operator_id']}: seed is invalid")
        if value.get("target") not in {
            "selected_fastq",
            "projected_measurement_payload",
        }:
            raise ValueError(f"{value['operator_id']}: target is invalid")
        validate_modules(value)
    if any(value.get("expected_seqspec_check") != "failure" for value in invalid):
        raise ValueError("schema-invalid controls must expect seqspec failure")
    scope = protocol.get("primary_structural_scope", {})
    if scope.get("include_deterministic") is not True:
        raise ValueError("primary structural scope must include deterministic defects")
    controls.unit_fraction(
        scope.get("stochastic_event_fraction"), "stochastic anchor fraction"
    )
    if scope.get("exclude_expected_seqspec_failure") is not True:
        raise ValueError("primary structural scope must exclude invalid specs")
    if scope.get("fastqc_module_scope") != "all_modules":
        raise ValueError("primary structural FastQC scope must include all modules")


def validate_modules(value: dict[str, Any]) -> None:
    modules = value.get("fastqc_modules")
    if (
        not isinstance(modules, list)
        or not modules
        or any(not isinstance(item, str) or not item.strip() for item in modules)
    ):
        raise ValueError(f"{value['operator_id']}: FastQC modules are invalid")


def validate_outputs(
    *,
    complementarity_id: str,
    cohort_split: str,
    configurations: set[str],
    controls: list[dict[str, Any]],
    skips: list[dict[str, Any]],
    protocol: dict[str, Any],
    policy: dict[str, Any] | None,
) -> dict[str, Any]:
    errors = []
    if any(value["complementarity_id"] != complementarity_id for value in controls):
        errors.append("control complementarity identifiers differ")
    condition_ids = [value["condition_id"] for value in controls]
    if len(condition_ids) != len(set(condition_ids)):
        errors.append("control condition identifiers are not unique")
    if any(
        runtime.sha256_json(control_identity(value))[:16] != value["condition_id"]
        for value in controls
    ):
        errors.append("control condition identifiers are not content-addressed")
    represented = {
        (value["configuration_accession"], value["operator_id"]) for value in controls
    } | {(value["configuration_accession"], value["operator_id"]) for value in skips}
    operator_ids = {
        value["operator_id"]
        for field in (
            "quality_controls",
            "composition_controls",
            "schema_invalid_controls",
        )
        for value in protocol[field]
    }
    expected = {
        (configuration, operator)
        for configuration in configurations
        for operator in operator_ids
    }
    if represented != expected:
        errors.append("controls and skips do not cover every configuration/operator")
    quality_counts = Counter(
        (value["configuration_accession"], value["operator_id"])
        for value in controls
        if value["condition_kind"] == "quality"
    )
    expected_quality = {
        value["operator_id"]: (
            1 if cohort_split == "evaluation" else len(value["severity_candidates"])
        )
        for value in protocol["quality_controls"]
    }
    if any(
        quality_counts[(configuration, operator)] != count
        for configuration in configurations
        for operator, count in expected_quality.items()
    ):
        errors.append("quality control severity counts do not reconcile")
    if (cohort_split == "evaluation") != (policy is not None):
        errors.append("quality policy use differs from cohort split")
    return {
        "schema_version": SCHEMA_VERSION,
        "complementarity_id": complementarity_id,
        "generated_at": utc_now(),
        "valid": not errors,
        "errors": errors,
        "counts": {
            "configurations": len(configurations),
            "controls": len(controls),
            "skips": len(skips),
            "controls_by_kind": dict(
                sorted(Counter(value["condition_kind"] for value in controls).items())
            ),
            "controls_by_operator": dict(
                sorted(Counter(value["operator_id"] for value in controls).items())
            ),
        },
    }


def public_control_row(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "complementarity_id": value["complementarity_id"],
        "condition_id": value["condition_id"],
        "base_materialization_id": value["base_materialization_id"],
        "selection_id": value["selection_id"],
        "configuration_accession": value["configuration_accession"],
        "family_id": value["family_id"],
        "modality": value["modality"],
        "condition_kind": value["condition_kind"],
        "operator_id": value["operator_id"],
        "variant": value["variant"],
        "severity_index": value["severity_index"],
        "event_fraction": value["event_fraction"],
        "mutation_seed": value["mutation_seed"],
        "mutated_record_count": value["mutated_record_count"],
        "total_record_count": value["total_record_count"],
        "target_case_ids": ";".join(value["target"]["case_ids"]),
        "target_read_ids": ";".join(value["target"]["read_ids"]),
        "target_region_ids": ";".join(value["target"]["region_ids"]),
        "expected_seqspec_check": value["expected_seqspec_check"],
        "observed_seqspec_check": value["observed_seqspec_check"],
    }


def skip_row(
    clean: dict[str, Any],
    operator: dict[str, Any],
    reason_code: str,
    reason: str,
) -> dict[str, str]:
    return {
        "configuration_accession": clean["configuration_accession"],
        "operator_id": operator["operator_id"],
        "variant": operator["variant"],
        "reason_code": reason_code,
        "reason": reason,
    }


def control_directory(
    configuration: str, operator_id: str, variant: str, severity_index: int | None
) -> Path:
    parts = [
        "controls",
        perturb.safe_name(configuration),
        operator_id,
        perturb.safe_name(variant),
    ]
    if severity_index is not None:
        parts.append(f"severity-{severity_index:02d}")
    return Path(*parts)


def verified_path(identity: Any, label: str) -> Path:
    if not isinstance(identity, dict):
        raise ValueError(f"{label} identity is missing")
    path = Path(str(identity.get("path", ""))).resolve()
    if runtime.file_sha256(path) != identity.get("sha256"):
        raise ValueError(f"{label} hash changed")
    return path


def required_string(value: dict[str, Any], field: str) -> str:
    result = value.get(field)
    if not isinstance(result, str) or not result.strip():
        raise ValueError(f"{field} must be a nonempty string")
    return result


def tool_identities(
    seqspec: dict[str, Any], yq: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    return {
        "seqspec": seqspec,
        "yq": yq,
        "materializer": runtime.script_identity(
            Path(__file__).resolve(), version=MATERIALIZER_VERSION
        ),
        "control_runtime": runtime.script_identity(
            Path(controls.__file__).resolve(), version=controls.RUNTIME_VERSION
        ),
        "perturbation_runtime": runtime.script_identity(
            Path(perturb.__file__).resolve(), version=perturb.RUNTIME_VERSION
        ),
        "paper_runtime": runtime.script_identity(
            Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
        ),
    }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
