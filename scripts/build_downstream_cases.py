#!/usr/bin/env python3
"""Freeze confirmed downstream before/after case studies before execution."""

from __future__ import annotations

import argparse
import math
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import manage_audit_reviews as reviews
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import manage_audit_reviews as reviews
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
TOOL_VERSION = "0.1.0"
CASE_FIELDS = (
    "selection_id",
    "selection_rank",
    "case_id",
    "configuration_accession",
    "modality",
    "final_classification",
    "matched_assay_family",
    "access_class",
    "pipeline_id",
    "pipeline_version",
    "endpoint_id",
    "endpoint_label",
    "endpoint_unit",
    "endpoint_direction",
    "minimum_meaningful_change",
    "shared_input_roles",
    "before_input_count",
    "after_input_count",
    "changed_input_roles",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Freeze confirmed downstream before/after case studies."
    )
    parser.add_argument("--review-analysis-manifest", required=True, type=Path)
    parser.add_argument("--case-registry", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        path = build_cases(
            review_analysis_path=args.review_analysis_manifest.resolve(),
            registry_path=args.case_registry.resolve(),
            protocol_path=args.protocol.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"build_downstream_cases: {error}", file=sys.stderr)
        return 1
    print(path)
    return 0


def build_cases(
    *,
    review_analysis_path: Path,
    registry_path: Path,
    protocol_path: Path,
    output_root: Path,
) -> Path:
    refuse_output_root(output_root)
    analysis, outcomes = load_review_analysis(review_analysis_path)
    protocol = load_protocol(protocol_path)
    registry = load_registry(registry_path, analysis["analysis_id"])
    outcomes_by_case = index_unique(outcomes, "case_id", "review outcomes")
    materialized = []
    for entry in registry["cases"]:
        case_id = entry["case_id"]
        outcome = outcomes_by_case.get(case_id)
        if outcome is None:
            raise ValueError(f"registry case is not a reviewed case: {case_id}")
        materialized.append(
            materialize_case(
                entry=entry,
                outcome=outcome,
                protocol=protocol,
                registry_root=registry_path.parent,
            )
        )
    eligible = [row for row in materialized if row["eligible"]]
    selection = protocol["selection"]
    minimum = selection["minimum_case_count"]
    if len(eligible) < minimum:
        raise ValueError(
            f"registry has {len(eligible)} eligible downstream cases; "
            f"{minimum} are required"
        )
    target = min(selection["target_case_count"], len(eligible))
    selected = select_stratified(
        eligible,
        target=target,
        fields=selection["stratification_fields"],
        seed=selection["selection_seed"],
    )
    for rank, row in enumerate(selected, 1):
        row["selection_rank"] = rank
    tool = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    stable = {
        "schema_version": SCHEMA_VERSION,
        "review_analysis_id": analysis["analysis_id"],
        "tool": runtime.functional_script_identity(tool),
        "inputs": {
            "review_analysis": runtime.file_sha256(review_analysis_path),
            "case_registry": runtime.file_sha256(registry_path),
            "protocol": runtime.file_sha256(protocol_path),
        },
        "selected_cases": selected,
    }
    selection_id = runtime.sha256_json(stable)[:16]
    rows = [case_table_row(selection_id, row) for row in selected]
    table_path = output_root / "tables" / "downstream_cases.csv"
    case_path = output_root / "cases" / "selected_cases.json"
    runtime.write_csv(table_path, rows, list(CASE_FIELDS))
    runtime.write_json(
        case_path,
        {
            "schema_version": SCHEMA_VERSION,
            "selection_id": selection_id,
            "cases": selected,
        },
    )
    checks = {
        "minimum_case_count_met": len(selected) >= minimum,
        "maximum_case_count_met": len(selected) <= selection["maximum_case_count"],
        "case_ids_unique": len({row["case_id"] for row in selected}) == len(selected),
        "all_cases_confirmed": all(row["confirmed_problem"] for row in selected),
        "all_cases_ready": all(row["eligible"] for row in selected),
        "all_cases_pair_before_after": all(
            set(row["conditions"]) == {"before", "after"} for row in selected
        ),
        "all_endpoints_predeclared": all(
            row["endpoint"]["direction"]
            in protocol["case_contract"]["endpoint_directions"]
            for row in selected
        ),
    }
    validation = {
        **stable,
        "selection_id": selection_id,
        "generated_at": utc_now(),
        "valid": all(checks.values()),
        "checks": checks,
        "counts": {
            "registered_cases": len(materialized),
            "eligible_cases": len(eligible),
            "selected_cases": len(selected),
            "selected_by_defect_class": dict(
                sorted(row_counts(selected, "final_classification").items())
            ),
            "selected_by_assay_family": dict(
                sorted(row_counts(selected, "matched_assay_family").items())
            ),
            "selected_by_access_class": dict(
                sorted(row_counts(selected, "access_class").items())
            ),
            "selected_by_pipeline": dict(
                sorted(row_counts(selected, "pipeline_id").items())
            ),
        },
        "outputs": {
            "case_table": runtime.file_identity(table_path),
            "selected_cases": runtime.file_identity(case_path),
        },
        "tools": {"builder": tool},
    }
    validation_path = output_root / "validation" / "downstream_cases.json"
    runtime.write_json(validation_path, validation)
    if not validation["valid"]:
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise ValueError(f"downstream case selection validation failed: {failed}")
    manifest_path = output_root / "manifests" / "downstream_cases.json"
    manifest = {
        **stable,
        "selection_id": selection_id,
        "generated_at": utc_now(),
        "valid": True,
        "frozen": True,
        "counts": validation["counts"],
        "inputs": {
            "review_analysis": runtime.file_identity(review_analysis_path),
            "case_registry": runtime.file_identity(registry_path),
            "protocol": runtime.file_identity(protocol_path),
        },
        "outputs": {
            **validation["outputs"],
            "validation": runtime.file_identity(validation_path),
        },
        "tools": {"builder": tool},
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest_path


def load_review_analysis(path: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    analysis = runtime.load_json(path)
    if (
        analysis.get("schema_version") != SCHEMA_VERSION
        or analysis.get("valid") is not True
    ):
        raise ValueError("audit review analysis is not valid")
    analysis_id = str(analysis.get("analysis_id", ""))
    if not analysis_id:
        raise ValueError("audit review analysis id is missing")
    validation_identity = analysis.get("outputs", {}).get("validation", {})
    validation_path = identity_path(validation_identity, "review validation")
    verify_identity(validation_path, validation_identity, "review validation")
    validation = runtime.load_json(validation_path)
    if (
        validation.get("valid") is not True
        or validation.get("analysis_id") != analysis_id
    ):
        raise ValueError("audit review validation does not match its analysis")
    validation_stable = {
        key: validation.get(key)
        for key in (
            "schema_version",
            "merge_id",
            "selection_id",
            "tool",
            "inputs",
        )
    }
    if runtime.sha256_json(validation_stable)[:16] != analysis_id:
        raise ValueError("audit review validation id is not content-addressed")
    outcome_identity = analysis.get("outputs", {}).get("outcomes", {})
    outcome_path = identity_path(outcome_identity, "review outcomes")
    verify_identity(outcome_path, outcome_identity, "review outcomes")
    outcomes = runtime.read_csv(outcome_path)
    if not outcomes:
        raise ValueError("audit review analysis has no outcomes")
    if any(row.get("analysis_id") not in {None, "", analysis_id} for row in outcomes):
        raise ValueError("audit review outcome analysis ids differ")
    return analysis, outcomes


def load_registry(path: Path, analysis_id: str) -> dict[str, Any]:
    registry = runtime.load_json(path)
    if registry.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("downstream case registry schema is unsupported")
    if registry.get("review_analysis_id") != analysis_id:
        raise ValueError("downstream registry belongs to another review analysis")
    if set(registry) != {"schema_version", "review_analysis_id", "cases"}:
        raise ValueError("downstream case registry fields are invalid")
    cases = registry.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("downstream case registry has no cases")
    case_ids = [entry.get("case_id") for entry in cases if isinstance(entry, dict)]
    if len(case_ids) != len(cases) or len(case_ids) != len(set(case_ids)):
        raise ValueError("downstream registry case ids are missing or duplicated")
    return registry


def materialize_case(
    *,
    entry: dict[str, Any],
    outcome: dict[str, str],
    protocol: dict[str, Any],
    registry_root: Path,
) -> dict[str, Any]:
    required_entry_fields = {
        "case_id",
        "access_available",
        "pipeline_available",
        "presentation_ready",
        "pipeline",
        "endpoint",
        "shared_input_roles",
        "before",
        "after",
    }
    if set(entry) != required_entry_fields:
        raise ValueError(
            f"{entry.get('case_id', '<unknown>')}: registry fields are invalid"
        )
    case_id = str(entry["case_id"])
    if outcome.get("source_kind") != "candidate" or not parse_bool(
        outcome.get("confirmed_problem"), "confirmed problem"
    ):
        raise ValueError(f"{case_id}: downstream case is not a confirmed candidate")
    required_classes = set(protocol["selection"]["required_final_classifications"])
    final_classification = outcome.get("final_classification", "")
    if final_classification not in required_classes:
        raise ValueError(f"{case_id}: final classification is not eligible")
    readiness = {
        field: entry.get(field)
        for field in protocol["selection"]["required_readiness_flags"]
    }
    if any(not isinstance(value, bool) for value in readiness.values()):
        raise ValueError(f"{case_id}: readiness flags must be boolean")
    pipeline = entry.get("pipeline")
    if not isinstance(pipeline, dict) or set(pipeline) != {"id", "version"}:
        raise ValueError(f"{case_id}: pipeline identity is invalid")
    pipeline_id = nonempty_string(pipeline.get("id"), f"{case_id} pipeline id")
    pipeline_version = nonempty_string(
        pipeline.get("version"), f"{case_id} pipeline version"
    )
    endpoint = validate_endpoint(entry.get("endpoint"), case_id, protocol)
    shared_roles = entry.get("shared_input_roles")
    if (
        not isinstance(shared_roles, list)
        or len(shared_roles) < protocol["case_contract"]["minimum_shared_input_roles"]
        or any(not isinstance(role, str) or not role for role in shared_roles)
        or len(shared_roles) != len(set(shared_roles))
    ):
        raise ValueError(f"{case_id}: shared input roles are invalid")
    eligible = all(readiness.values())
    if not eligible:
        for condition in protocol["case_contract"]["conditions"]:
            validate_unready_condition(entry.get(condition), condition, case_id)
        return {
            "case_id": case_id,
            "configuration_accession": outcome.get("configuration_accession", ""),
            "modality": outcome.get("modality", ""),
            "final_classification": final_classification,
            "matched_assay_family": outcome.get("matched_assay_family", ""),
            "access_class": outcome.get("access_class", ""),
            "confirmed_problem": True,
            "pipeline_id": pipeline_id,
            "pipeline_version": pipeline_version,
            "endpoint": endpoint,
            "shared_input_roles": sorted(shared_roles),
            "changed_input_roles": [],
            "readiness": readiness,
            "eligible": False,
            "conditions": {},
        }
    conditions = {
        condition: materialize_condition(
            value=entry.get(condition),
            condition=condition,
            case_id=case_id,
            registry_root=registry_root,
            protocol=protocol,
        )
        for condition in protocol["case_contract"]["conditions"]
    }
    before_inputs = {value["role"]: value for value in conditions["before"]["inputs"]}
    after_inputs = {value["role"]: value for value in conditions["after"]["inputs"]}
    if set(before_inputs) != set(after_inputs):
        raise ValueError(f"{case_id}: before and after input roles differ")
    for role in shared_roles:
        if role not in before_inputs or role not in after_inputs:
            raise ValueError(f"{case_id}: shared input role is missing: {role}")
        if before_inputs[role]["sha256"] != after_inputs[role]["sha256"]:
            raise ValueError(f"{case_id}: shared input changed: {role}")
    all_roles = set(before_inputs) | set(after_inputs)
    changed_roles = sorted(
        role
        for role in all_roles - set(shared_roles)
        if role not in before_inputs
        or role not in after_inputs
        or before_inputs[role]["sha256"] != after_inputs[role]["sha256"]
    )
    if (
        protocol["case_contract"]["require_changed_nonshared_input"]
        and not changed_roles
    ):
        raise ValueError(f"{case_id}: before and after have no changed input")
    if (
        conditions["before"]["executable"]["sha256"]
        != conditions["after"]["executable"]["sha256"]
    ):
        raise ValueError(f"{case_id}: before and after use different executables")
    if protocol["case_contract"]["require_condition_command_parity"]:
        if normalized_condition_command(
            conditions["before"]
        ) != normalized_condition_command(conditions["after"]):
            raise ValueError(
                f"{case_id}: before and after commands differ beyond declared inputs"
            )
    return {
        "case_id": case_id,
        "configuration_accession": outcome.get("configuration_accession", ""),
        "modality": outcome.get("modality", ""),
        "final_classification": final_classification,
        "matched_assay_family": outcome.get("matched_assay_family", ""),
        "access_class": outcome.get("access_class", ""),
        "confirmed_problem": True,
        "pipeline_id": pipeline_id,
        "pipeline_version": pipeline_version,
        "endpoint": endpoint,
        "shared_input_roles": sorted(shared_roles),
        "changed_input_roles": changed_roles,
        "readiness": readiness,
        "eligible": True,
        "conditions": conditions,
    }


def validate_unready_condition(value: Any, condition: str, case_id: str) -> None:
    if not isinstance(value, dict) or set(value) != {"inputs", "command"}:
        raise ValueError(f"{case_id}/{condition}: condition fields are invalid")
    inputs = value.get("inputs")
    command = value.get("command")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError(f"{case_id}/{condition}: inputs are missing")
    roles = []
    for item in inputs:
        if not isinstance(item, dict) or set(item) != {"role", "path"}:
            raise ValueError(f"{case_id}/{condition}: input fields are invalid")
        roles.append(nonempty_string(item.get("role"), "input role"))
        nonempty_string(item.get("path"), "input path")
    if len(roles) != len(set(roles)):
        raise ValueError(f"{case_id}/{condition}: input roles are duplicated")
    if (
        not isinstance(command, list)
        or not command
        or any(not isinstance(argument, str) or not argument for argument in command)
    ):
        raise ValueError(f"{case_id}/{condition}: command is invalid")


def materialize_condition(
    *,
    value: Any,
    condition: str,
    case_id: str,
    registry_root: Path,
    protocol: dict[str, Any],
) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != {"inputs", "command"}:
        raise ValueError(f"{case_id}/{condition}: condition fields are invalid")
    inputs = value.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError(f"{case_id}/{condition}: inputs are missing")
    materialized_inputs = []
    input_path_aliases = []
    roles = set()
    for item in inputs:
        if not isinstance(item, dict) or set(item) != {"role", "path"}:
            raise ValueError(f"{case_id}/{condition}: input fields are invalid")
        role = nonempty_string(item.get("role"), f"{case_id}/{condition} input role")
        if role in roles:
            raise ValueError(f"{case_id}/{condition}: duplicate input role: {role}")
        roles.add(role)
        raw_path = Path(nonempty_string(item.get("path"), "input path"))
        if not raw_path.is_absolute():
            raw_path = registry_root / raw_path
        declared_path = str(raw_path.absolute())
        path = raw_path.resolve()
        if not path.is_file():
            raise ValueError(f"{case_id}/{condition}: input is not a file: {path}")
        materialized_inputs.append({"role": role, **runtime.file_identity(path)})
        input_path_aliases.append((declared_path, str(path)))
    command = value.get("command")
    if (
        not isinstance(command, list)
        or not command
        or any(not isinstance(arg, str) or not arg for arg in command)
    ):
        raise ValueError(f"{case_id}/{condition}: command is invalid")
    executable_path = resolve_executable(command[0])
    resolved_command = [str(executable_path), *command[1:]]
    for declared_path, resolved_path in sorted(
        input_path_aliases, key=lambda value: -len(value[0])
    ):
        resolved_command = [
            argument.replace(declared_path, resolved_path)
            for argument in resolved_command
        ]
    if protocol["case_contract"]["require_output_dir_placeholder"] and not any(
        "{output_dir}" in arg for arg in resolved_command
    ):
        raise ValueError(f"{case_id}/{condition}: command omits {{output_dir}}")
    if protocol["case_contract"]["require_input_paths_in_argv"]:
        for item in materialized_inputs:
            if not any(item["path"] in arg for arg in resolved_command):
                raise ValueError(
                    f"{case_id}/{condition}: command omits input role {item['role']}"
                )
    if protocol["case_contract"]["forbid_shell_command_strings"]:
        shell = executable_path.name in {"sh", "bash", "zsh", "dash"}
        if shell and "-c" in resolved_command[1:]:
            raise ValueError(
                f"{case_id}/{condition}: shell command strings are forbidden"
            )
    return {
        "inputs": materialized_inputs,
        "command_template": resolved_command,
        "executable": runtime.file_identity(executable_path),
    }


def validate_endpoint(
    value: Any, case_id: str, protocol: dict[str, Any]
) -> dict[str, Any]:
    fields = {"id", "label", "unit", "direction", "minimum_meaningful_change"}
    if not isinstance(value, dict) or set(value) != fields:
        raise ValueError(f"{case_id}: endpoint fields are invalid")
    endpoint_id = nonempty_string(value.get("id"), f"{case_id} endpoint id")
    label = nonempty_string(value.get("label"), f"{case_id} endpoint label")
    unit = nonempty_string(value.get("unit"), f"{case_id} endpoint unit")
    direction = value.get("direction")
    if direction not in protocol["case_contract"]["endpoint_directions"]:
        raise ValueError(f"{case_id}: endpoint direction is invalid")
    threshold = finite_positive(
        value.get("minimum_meaningful_change"), f"{case_id} endpoint threshold"
    )
    return {
        "id": endpoint_id,
        "label": label,
        "unit": unit,
        "direction": direction,
        "minimum_meaningful_change": threshold,
        "output_filename": protocol["case_contract"]["endpoint_output_filename"],
    }


def normalized_condition_command(condition: dict[str, Any]) -> list[str]:
    replacements = sorted(
        (
            (value["path"], f"{{input:{value['role']}}}")
            for value in condition["inputs"]
        ),
        key=lambda value: -len(value[0]),
    )
    normalized = []
    for argument in condition["command_template"]:
        value = argument
        for path, token in replacements:
            value = value.replace(path, token)
        normalized.append(value)
    return normalized


def select_stratified(
    cases: list[dict[str, Any]], *, target: int, fields: list[str], seed: int
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in cases:
        grouped[tuple(str(row[field]) for field in fields)].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: stable_rank(seed, "case", row["case_id"]))
    strata = sorted(grouped, key=lambda value: stable_rank(seed, "stratum", *value))
    selected = []
    while len(selected) < target:
        added = False
        for stratum in strata:
            if grouped[stratum] and len(selected) < target:
                selected.append(grouped[stratum].pop(0))
                added = True
        if not added:
            break
    if len(selected) != target:
        raise ValueError("downstream case selection did not reach its target")
    return selected


def case_table_row(selection_id: str, row: dict[str, Any]) -> dict[str, Any]:
    before = {value["role"]: value for value in row["conditions"]["before"]["inputs"]}
    after = {value["role"]: value for value in row["conditions"]["after"]["inputs"]}
    return {
        "selection_id": selection_id,
        "selection_rank": row["selection_rank"],
        "case_id": row["case_id"],
        "configuration_accession": row["configuration_accession"],
        "modality": row["modality"],
        "final_classification": row["final_classification"],
        "matched_assay_family": row["matched_assay_family"],
        "access_class": row["access_class"],
        "pipeline_id": row["pipeline_id"],
        "pipeline_version": row["pipeline_version"],
        "endpoint_id": row["endpoint"]["id"],
        "endpoint_label": row["endpoint"]["label"],
        "endpoint_unit": row["endpoint"]["unit"],
        "endpoint_direction": row["endpoint"]["direction"],
        "minimum_meaningful_change": row["endpoint"]["minimum_meaningful_change"],
        "shared_input_roles": ";".join(row["shared_input_roles"]),
        "before_input_count": len(before),
        "after_input_count": len(after),
        "changed_input_roles": ";".join(row["changed_input_roles"]),
    }


def load_protocol(path: Path) -> dict[str, Any]:
    protocol = runtime.load_json(path)
    if protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("downstream case protocol schema is unsupported")
    if protocol.get("experiment_id") != "downstream-case-studies":
        raise ValueError("downstream case experiment id is invalid")
    selection = protocol.get("selection", {})
    minimum = positive_int(selection.get("minimum_case_count"), "minimum cases")
    target = positive_int(selection.get("target_case_count"), "target cases")
    maximum = positive_int(selection.get("maximum_case_count"), "maximum cases")
    if not minimum <= target <= maximum:
        raise ValueError("downstream case count bounds are inconsistent")
    nonnegative_int(selection.get("selection_seed"), "selection seed")
    if selection.get("required_final_classifications") != list(
        reviews.CLASSIFICATIONS[:2]
    ):
        raise ValueError("downstream eligible classifications are invalid")
    expected_strata = [
        "final_classification",
        "matched_assay_family",
        "access_class",
        "pipeline_id",
    ]
    if selection.get("stratification_fields") != expected_strata:
        raise ValueError("downstream stratification fields are invalid")
    if selection.get("required_readiness_flags") != [
        "access_available",
        "pipeline_available",
        "presentation_ready",
    ]:
        raise ValueError("downstream readiness flags are invalid")
    contract = protocol.get("case_contract", {})
    if contract.get("conditions") != ["before", "after"]:
        raise ValueError("downstream conditions are invalid")
    if contract.get("endpoint_output_filename") != "endpoint.json":
        raise ValueError("downstream endpoint filename is invalid")
    if contract.get("endpoint_directions") != ["increase", "decrease"]:
        raise ValueError("downstream endpoint directions are invalid")
    positive_int(contract.get("minimum_shared_input_roles"), "shared input roles")
    for field in (
        "require_changed_nonshared_input",
        "require_input_paths_in_argv",
        "require_output_dir_placeholder",
        "require_condition_command_parity",
        "forbid_shell_command_strings",
    ):
        if contract.get(field) is not True:
            raise ValueError(f"downstream contract must enable {field}")
    positive_int(
        protocol.get("execution", {}).get("condition_timeout_seconds"),
        "condition timeout",
    )
    targets = protocol.get("scientific_targets", {})
    positive_int(
        targets.get("minimum_presentable_confirmed_cases"), "presentable cases"
    )
    positive_int(
        targets.get("minimum_cases_with_expected_endpoint_change"),
        "changed endpoint cases",
    )
    return protocol


def resolve_executable(value: str) -> Path:
    path = Path(value)
    resolved = path.resolve() if path.is_absolute() else None
    if resolved is None or not resolved.is_file():
        found = shutil.which(value)
        if not found:
            raise ValueError(f"pipeline executable is unavailable: {value}")
        resolved = Path(found).resolve()
    return resolved


def verify_identity(path: Path, identity: Any, label: str) -> None:
    if not isinstance(identity, dict) or runtime.file_sha256(path) != identity.get(
        "sha256"
    ):
        raise ValueError(f"{label} hash changed")


def identity_path(identity: Any, label: str) -> Path:
    if not isinstance(identity, dict) or not str(identity.get("path", "")).strip():
        raise ValueError(f"{label} path is missing")
    return Path(str(identity["path"])).resolve()


def index_unique(
    rows: list[dict[str, str]], field: str, label: str
) -> dict[str, dict[str, str]]:
    result = {}
    for row in rows:
        key = row.get(field, "")
        if not key or key in result:
            raise ValueError(f"{label} {field} values are incomplete or duplicated")
        result[key] = row
    return result


def row_counts(rows: list[dict[str, Any]], field: str) -> Counter[str]:
    return Counter(str(row[field]) for row in rows)


def stable_rank(seed: int, *values: str) -> str:
    return runtime.sha256_json([seed, *values])


def parse_bool(value: Any, label: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise ValueError(f"{label} is not boolean")


def nonempty_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty string")
    return value.strip()


def finite_positive(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise ValueError(f"{label} must be finite and positive")
    return parsed


def positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def nonnegative_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def refuse_output_root(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"output root is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
