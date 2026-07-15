#!/usr/bin/env python3
"""Execute frozen downstream before/after case studies."""

from __future__ import annotations

import argparse
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import build_downstream_cases as builder
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import build_downstream_cases as builder
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
TOOL_VERSION = "0.1.0"
RUN_FIELDS = (
    "execution_id",
    "selection_id",
    "case_id",
    "condition",
    "pipeline_id",
    "pipeline_version",
    "command_json",
    "status",
    "message",
    "exit_code",
    "timed_out",
    "wall_time_seconds",
    "user_cpu_seconds",
    "system_cpu_seconds",
    "peak_resident_memory_bytes",
    "endpoint_id",
    "endpoint_unit",
    "endpoint_value",
    "endpoint_path",
    "endpoint_sha256",
    "output_inventory_json",
    "stdout_path",
    "stdout_sha256",
    "stderr_path",
    "stderr_sha256",
)
RESULT_FIELDS = (
    "execution_id",
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
    "before_value",
    "after_value",
    "after_minus_before",
    "directional_change",
    "expected_change_met",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Execute frozen downstream before/after case studies."
    )
    parser.add_argument("--case-manifest", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        path = run_cases(
            case_manifest_path=args.case_manifest.resolve(),
            protocol_path=args.protocol.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"run_downstream_cases: {error}", file=sys.stderr)
        return 1
    print(path)
    return 0


def run_cases(
    *, case_manifest_path: Path, protocol_path: Path, output_root: Path
) -> Path:
    refuse_output_root(output_root)
    case_manifest, cases = load_case_manifest(case_manifest_path, protocol_path)
    protocol = builder.load_protocol(protocol_path)
    verify_case_inputs(cases, protocol)
    tool = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    stable = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": case_manifest["selection_id"],
        "tool": runtime.functional_script_identity(tool),
        "inputs": {
            "case_manifest": runtime.file_sha256(case_manifest_path),
            "protocol": runtime.file_sha256(protocol_path),
        },
    }
    execution_id = runtime.sha256_json(stable)[:16]
    timeout = protocol["execution"]["condition_timeout_seconds"]
    run_rows = []
    values: dict[tuple[str, str], float] = {}
    for case in cases:
        for condition in protocol["case_contract"]["conditions"]:
            row, value = execute_condition(
                execution_id=execution_id,
                selection_id=case_manifest["selection_id"],
                case=case,
                condition=condition,
                timeout_seconds=timeout,
                output_root=output_root,
            )
            run_rows.append(row)
            if value is not None:
                values[(case["case_id"], condition)] = value
    result_rows = build_results(
        execution_id=execution_id,
        selection_id=case_manifest["selection_id"],
        cases=cases,
        values=values,
    )
    run_path = output_root / "tables" / "condition_runs.csv"
    result_path = output_root / "tables" / "paired_endpoints.csv"
    runtime.write_csv(run_path, run_rows, list(RUN_FIELDS))
    runtime.write_csv(result_path, result_rows, list(RESULT_FIELDS))

    expected_run_count = len(cases) * len(protocol["case_contract"]["conditions"])
    checks = {
        "condition_run_count_reconciles": len(run_rows) == expected_run_count,
        "all_conditions_completed": all(
            row["status"] == "completed" for row in run_rows
        ),
        "paired_result_count_reconciles": len(result_rows) == len(cases),
        "endpoint_values_finite": all(
            math.isfinite(float(row["before_value"]))
            and math.isfinite(float(row["after_value"]))
            for row in result_rows
        ),
        "case_ids_unique": len({row["case_id"] for row in result_rows})
        == len(result_rows),
    }
    targets = protocol["scientific_targets"]
    changed_cases = sum(row["expected_change_met"] for row in result_rows)
    scientific_targets = {
        "presentable_confirmed_case_target_met": len(result_rows)
        >= targets["minimum_presentable_confirmed_cases"],
        "expected_endpoint_change_target_met": changed_cases
        >= targets["minimum_cases_with_expected_endpoint_change"],
    }
    validation = {
        **stable,
        "execution_id": execution_id,
        "generated_at": utc_now(),
        "valid": all(checks.values()),
        "scientific_targets_met": all(scientific_targets.values()),
        "checks": checks,
        "scientific_targets": scientific_targets,
        "counts": {
            "selected_cases": len(cases),
            "condition_runs": len(run_rows),
            "completed_conditions": sum(
                row["status"] == "completed" for row in run_rows
            ),
            "paired_results": len(result_rows),
            "cases_with_expected_endpoint_change": changed_cases,
        },
        "outputs": {
            "condition_runs": runtime.file_identity(run_path),
            "paired_endpoints": runtime.file_identity(result_path),
        },
        "tools": {"runner": tool},
    }
    validation_path = output_root / "validation" / "downstream_execution.json"
    runtime.write_json(validation_path, validation)
    if not validation["valid"]:
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise ValueError(f"downstream execution validation failed: {failed}")
    manifest_path = output_root / "manifests" / "downstream_execution.json"
    manifest = {
        **stable,
        "execution_id": execution_id,
        "generated_at": utc_now(),
        "valid": True,
        "scientific_targets_met": validation["scientific_targets_met"],
        "inputs": {
            "case_manifest": runtime.file_identity(case_manifest_path),
            "protocol": runtime.file_identity(protocol_path),
        },
        "counts": validation["counts"],
        "outputs": {
            **validation["outputs"],
            "validation": runtime.file_identity(validation_path),
        },
        "tools": {"runner": tool},
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest_path


def load_case_manifest(
    manifest_path: Path, protocol_path: Path
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    manifest = runtime.load_json(manifest_path)
    if (
        manifest.get("schema_version") != SCHEMA_VERSION
        or manifest.get("valid") is not True
        or manifest.get("frozen") is not True
    ):
        raise ValueError("downstream case manifest is not valid and frozen")
    stable = {
        "schema_version": manifest.get("schema_version"),
        "review_analysis_id": manifest.get("review_analysis_id"),
        "tool": manifest.get("tool"),
        "inputs": {
            name: manifest.get("inputs", {}).get(name, {}).get("sha256")
            for name in ("review_analysis", "case_registry", "protocol")
        },
        "selected_cases": manifest.get("selected_cases"),
    }
    if runtime.sha256_json(stable)[:16] != manifest.get("selection_id"):
        raise ValueError("downstream case selection id is not content-addressed")
    protocol_identity = manifest.get("inputs", {}).get("protocol", {})
    if runtime.file_sha256(protocol_path) != protocol_identity.get("sha256"):
        raise ValueError("downstream execution protocol differs from selection")
    validation_identity = manifest.get("outputs", {}).get("validation", {})
    validation_path = identity_path(validation_identity, "case selection validation")
    verify_identity(validation_path, validation_identity, "case selection validation")
    validation = runtime.load_json(validation_path)
    if validation.get("valid") is not True or validation.get(
        "selection_id"
    ) != manifest.get("selection_id"):
        raise ValueError("downstream case selection validation differs")
    case_identity = manifest.get("outputs", {}).get("selected_cases", {})
    case_path = identity_path(case_identity, "selected case file")
    verify_identity(case_path, case_identity, "selected case file")
    case_file = runtime.load_json(case_path)
    cases = case_file.get("cases")
    if (
        case_file.get("selection_id") != manifest.get("selection_id")
        or cases != manifest.get("selected_cases")
        or not isinstance(cases, list)
        or not cases
    ):
        raise ValueError("selected downstream cases differ from manifest")
    return manifest, cases


def verify_case_inputs(cases: list[dict[str, Any]], protocol: dict[str, Any]) -> None:
    ranks = [case.get("selection_rank") for case in cases]
    if ranks != list(range(1, len(cases) + 1)):
        raise ValueError("downstream case selection ranks do not reconcile")
    for case in cases:
        if set(case.get("conditions", {})) != {"before", "after"}:
            raise ValueError(f"{case.get('case_id', '')}: conditions are incomplete")
        before = {item["role"]: item for item in case["conditions"]["before"]["inputs"]}
        after = {item["role"]: item for item in case["conditions"]["after"]["inputs"]}
        if set(before) != set(after):
            raise ValueError(f"{case['case_id']}: before and after input roles differ")
        shared = set(case.get("shared_input_roles", []))
        if not shared or not shared <= set(before):
            raise ValueError(f"{case['case_id']}: shared input roles are invalid")
        if any(before[role]["sha256"] != after[role]["sha256"] for role in shared):
            raise ValueError(f"{case['case_id']}: shared input identity changed")
        changed = sorted(
            role
            for role in set(before) - shared
            if before[role]["sha256"] != after[role]["sha256"]
        )
        if changed != case.get("changed_input_roles") or not changed:
            raise ValueError(f"{case['case_id']}: changed input roles do not reconcile")
        if builder.normalized_condition_command(
            case["conditions"]["before"]
        ) != builder.normalized_condition_command(case["conditions"]["after"]):
            raise ValueError(f"{case['case_id']}: condition commands do not reconcile")
        if (
            case.get("endpoint", {}).get("direction")
            not in protocol["case_contract"]["endpoint_directions"]
        ):
            raise ValueError(f"{case['case_id']}: endpoint direction is invalid")
        for condition in ("before", "after"):
            value = case["conditions"][condition]
            if not any("{output_dir}" in arg for arg in value["command_template"]):
                raise ValueError(
                    f"{case['case_id']}/{condition}: output placeholder is missing"
                )
            executable = value["executable"]
            path = identity_path(executable, "pipeline executable")
            verify_identity(path, executable, "pipeline executable")
            if str(path) != value["command_template"][0]:
                raise ValueError(
                    f"{case['case_id']}/{condition}: executable command changed"
                )
            for item in value["inputs"]:
                path = identity_path(item, f"{case['case_id']}/{condition} input")
                verify_identity(path, item, f"{case['case_id']}/{condition} input")


def execute_condition(
    *,
    execution_id: str,
    selection_id: str,
    case: dict[str, Any],
    condition: str,
    timeout_seconds: int,
    output_root: Path,
) -> tuple[dict[str, Any], float | None]:
    condition_root = output_root / "artifacts" / case["case_id"] / condition
    condition_root.mkdir(parents=True, exist_ok=False)
    log_root = output_root / "logs" / case["case_id"]
    stdout_path = log_root / f"{condition}.stdout.txt"
    stderr_path = log_root / f"{condition}.stderr.txt"
    command = [
        argument.replace("{output_dir}", str(condition_root))
        .replace("{case_id}", case["case_id"])
        .replace("{condition}", condition)
        for argument in case["conditions"][condition]["command_template"]
    ]
    measurement = runtime.run_measured_command(
        command,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        timeout_seconds=timeout_seconds,
    )
    status = "completed"
    message = ""
    endpoint_value = None
    endpoint_identity: dict[str, Any] = {}
    if measurement["timed_out"]:
        status = "timed_out"
        message = "condition exceeded the frozen timeout"
    elif measurement["exit_code"] != 0:
        status = "process_failure"
        message = stderr_path.read_text(encoding="utf-8", errors="replace").strip()
    else:
        try:
            endpoint_path = condition_root / case["endpoint"]["output_filename"]
            endpoint_value = load_endpoint(endpoint_path, case["endpoint"])
            endpoint_identity = runtime.file_identity(endpoint_path)
        except (OSError, ValueError) as error:
            status = "invalid_endpoint"
            message = str(error)
    try:
        inventory = output_inventory(condition_root)
    except ValueError as error:
        status = "invalid_artifact"
        message = str(error)
        inventory = []
        endpoint_value = None
        endpoint_identity = {}
    row = {
        "execution_id": execution_id,
        "selection_id": selection_id,
        "case_id": case["case_id"],
        "condition": condition,
        "pipeline_id": case["pipeline_id"],
        "pipeline_version": case["pipeline_version"],
        "command_json": runtime.canonical_json(command),
        "status": status,
        "message": message,
        "exit_code": measurement["exit_code"],
        "timed_out": measurement["timed_out"],
        "wall_time_seconds": measurement["wall_time_seconds"],
        "user_cpu_seconds": measurement["user_cpu_seconds"],
        "system_cpu_seconds": measurement["system_cpu_seconds"],
        "peak_resident_memory_bytes": measurement["peak_resident_memory_bytes"],
        "endpoint_id": case["endpoint"]["id"] if endpoint_value is not None else "",
        "endpoint_unit": (
            case["endpoint"]["unit"] if endpoint_value is not None else ""
        ),
        "endpoint_value": endpoint_value,
        "endpoint_path": endpoint_identity.get("path", ""),
        "endpoint_sha256": endpoint_identity.get("sha256", ""),
        "output_inventory_json": runtime.canonical_json(inventory),
        "stdout_path": measurement["stdout"]["path"],
        "stdout_sha256": measurement["stdout"]["sha256"],
        "stderr_path": measurement["stderr"]["path"],
        "stderr_sha256": measurement["stderr"]["sha256"],
    }
    return row, endpoint_value


def load_endpoint(path: Path, expected: dict[str, Any]) -> float:
    value = runtime.load_json(path)
    if set(value) != {"endpoint_id", "value", "unit"}:
        raise ValueError("endpoint output fields are invalid")
    if value["endpoint_id"] != expected["id"] or value["unit"] != expected["unit"]:
        raise ValueError("endpoint output identity differs from the frozen endpoint")
    observed = value["value"]
    if isinstance(observed, bool) or not isinstance(observed, (int, float)):
        raise ValueError("endpoint value is not numeric")
    observed = float(observed)
    if not math.isfinite(observed):
        raise ValueError("endpoint value is not finite")
    return observed


def output_inventory(root: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"condition output contains a symbolic link: {path}")
        if path.is_file():
            rows.append(
                {
                    "relative_path": str(path.relative_to(root)),
                    **runtime.file_identity(path),
                }
            )
    return rows


def build_results(
    *,
    execution_id: str,
    selection_id: str,
    cases: list[dict[str, Any]],
    values: dict[tuple[str, str], float],
) -> list[dict[str, Any]]:
    rows = []
    for case in cases:
        key_before = (case["case_id"], "before")
        key_after = (case["case_id"], "after")
        if key_before not in values or key_after not in values:
            continue
        before = values[key_before]
        after = values[key_after]
        delta = after - before
        directional = delta if case["endpoint"]["direction"] == "increase" else -delta
        rows.append(
            {
                "execution_id": execution_id,
                "selection_id": selection_id,
                "selection_rank": case["selection_rank"],
                "case_id": case["case_id"],
                "configuration_accession": case["configuration_accession"],
                "modality": case["modality"],
                "final_classification": case["final_classification"],
                "matched_assay_family": case["matched_assay_family"],
                "access_class": case["access_class"],
                "pipeline_id": case["pipeline_id"],
                "pipeline_version": case["pipeline_version"],
                "endpoint_id": case["endpoint"]["id"],
                "endpoint_label": case["endpoint"]["label"],
                "endpoint_unit": case["endpoint"]["unit"],
                "endpoint_direction": case["endpoint"]["direction"],
                "minimum_meaningful_change": case["endpoint"][
                    "minimum_meaningful_change"
                ],
                "before_value": before,
                "after_value": after,
                "after_minus_before": delta,
                "directional_change": directional,
                "expected_change_met": directional
                >= case["endpoint"]["minimum_meaningful_change"],
            }
        )
    return rows


def verify_identity(path: Path, identity: Any, label: str) -> None:
    if not isinstance(identity, dict) or runtime.file_sha256(path) != identity.get(
        "sha256"
    ):
        raise ValueError(f"{label} hash changed")


def identity_path(identity: Any, label: str) -> Path:
    if not isinstance(identity, dict) or not str(identity.get("path", "")).strip():
        raise ValueError(f"{label} path is missing")
    return Path(str(identity["path"])).resolve()


def refuse_output_root(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"output root is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
