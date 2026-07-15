#!/usr/bin/env python3
"""Regenerate paper analyses twice and freeze matching canonical results."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
TOOL_VERSION = "0.1.0"
CANONICAL_KINDS = ("table", "figure_data")
OUTPUT_KINDS = {*CANONICAL_KINDS, "figure"}
RESULT_DIRECTORIES = ("tables", "analysis", "figure_data", "figures")
SAFE_ID = re.compile(r"^[A-Za-z0-9._-]+$")
ENVIRONMENT_KEY = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
RUN_FIELDS = (
    "freeze_plan_id",
    "reproduction",
    "step_id",
    "status",
    "message",
    "command_json",
    "working_directory",
    "exit_code",
    "timed_out",
    "wall_time_seconds",
    "user_cpu_seconds",
    "system_cpu_seconds",
    "peak_resident_memory_bytes",
    "stdout_path",
    "stdout_sha256",
    "stderr_path",
    "stderr_sha256",
    "output_file_count",
    "output_tree_sha256",
)
COMPARISON_FIELDS = (
    "freeze_plan_id",
    "output_id",
    "step_id",
    "kind",
    "relative_path",
    "canonical",
    "reproduction_1_path",
    "reproduction_1_sha256",
    "reproduction_1_size_bytes",
    "reproduction_2_path",
    "reproduction_2_sha256",
    "reproduction_2_size_bytes",
    "byte_identical",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run two clean paper-analysis regenerations and freeze matching results."
        )
    )
    parser.add_argument("--study-manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        path = run_analysis_freeze(
            study_manifest_path=args.study_manifest.resolve(),
            output_root=args.output_root.resolve(),
        )
    except (OSError, ValueError) as error:
        print(f"run_analysis_freeze: {error}", file=sys.stderr)
        return 1
    print(path)
    return 0


def run_analysis_freeze(*, study_manifest_path: Path, output_root: Path) -> Path:
    study_manifest_path = study_manifest_path.resolve()
    output_root = output_root.resolve()
    study = load_study_manifest(study_manifest_path)
    ensure_output_separate(output_root, study)
    tool = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    runtime_tool = runtime.script_identity(
        Path(runtime.__file__).resolve(), version=runtime.RUNTIME_SCHEMA_VERSION
    )
    artifact_snapshot = snapshot_artifacts(study["inputs"])
    repository_snapshot = snapshot_repositories(study["repositories"])
    software_entries = [
        *study["software"],
        *internal_software_entries(study, tool, runtime_tool),
    ]
    software_snapshot = snapshot_software(software_entries, repository_snapshot)
    refuse_output_root(output_root)
    study_manifest_sha256 = runtime.file_sha256(study_manifest_path)
    stable = {
        "schema_version": SCHEMA_VERSION,
        "study_id": study["study_id"],
        "tools": {
            "analysis_freezer": runtime.functional_script_identity(tool),
            "paper_runtime": runtime.functional_script_identity(runtime_tool),
        },
        "study_manifest_sha256": study_manifest_sha256,
        "artifacts": stable_artifacts(artifact_snapshot),
        "repositories": stable_repositories(repository_snapshot),
        "software": stable_software(software_snapshot),
        "environment": study["environment"],
        "steps": stable_steps(study["steps"]),
    }
    freeze_plan_id = runtime.sha256_json(stable)[:16]

    run_rows: list[dict[str, Any]] = []
    outputs_by_reproduction: dict[int, dict[str, dict[str, Any]]] = {}
    for reproduction in range(1, study["repetitions"] + 1):
        outputs_by_reproduction[reproduction] = execute_reproduction(
            reproduction=reproduction,
            freeze_plan_id=freeze_plan_id,
            study=study,
            output_root=output_root,
            run_rows=run_rows,
        )

    run_path = output_root / "tables" / "analysis_runs.csv"
    runtime.write_csv(run_path, run_rows, list(RUN_FIELDS))
    comparison_rows = compare_outputs(
        freeze_plan_id=freeze_plan_id,
        study=study,
        outputs_by_reproduction=outputs_by_reproduction,
    )
    comparison_path = output_root / "tables" / "result_hashes.csv"
    runtime.write_csv(comparison_path, comparison_rows, list(COMPARISON_FIELDS))

    postflight_errors = {}
    artifacts_after = None
    repositories_after = None
    software_after = None
    study_manifest_after_sha256 = None
    try:
        study_manifest_after_sha256 = runtime.file_sha256(study_manifest_path)
    except OSError as error:
        postflight_errors["study_manifest"] = str(error)
    try:
        artifacts_after = snapshot_artifacts(study["inputs"])
    except ValueError as error:
        postflight_errors["artifacts"] = str(error)
    try:
        repositories_after = snapshot_repositories(
            study["repositories"], require_clean=False
        )
    except ValueError as error:
        postflight_errors["repositories"] = str(error)
    try:
        software_after = snapshot_software(
            software_entries, repositories_after or repository_snapshot
        )
    except ValueError as error:
        postflight_errors["software"] = str(error)
    expected_run_count = study["repetitions"] * len(study["steps"])
    expected_output_count = sum(len(step["outputs"]) for step in study["steps"])
    checks = {
        "two_reproductions_requested": study["repetitions"] == 2,
        "run_count_reconciles": len(run_rows) == expected_run_count,
        "all_analysis_steps_completed": all(
            row["status"] == "completed" for row in run_rows
        ),
        "reproduction_1_output_count_reconciles": len(
            outputs_by_reproduction.get(1, {})
        )
        == expected_output_count,
        "reproduction_2_output_count_reconciles": len(
            outputs_by_reproduction.get(2, {})
        )
        == expected_output_count,
        "comparison_count_reconciles": len(comparison_rows) == expected_output_count,
        "all_canonical_outputs_byte_identical": all(
            row["byte_identical"] for row in comparison_rows if row["canonical"]
        ),
        "study_manifest_unchanged": study_manifest_after_sha256
        == study_manifest_sha256,
        "all_source_artifacts_unchanged": artifacts_after is not None
        and stable_artifacts(artifact_snapshot) == stable_artifacts(artifacts_after),
        "all_software_unchanged": software_after is not None
        and stable_software(software_snapshot) == stable_software(software_after),
        "all_repositories_unchanged": repositories_after is not None
        and stable_repositories(repository_snapshot)
        == stable_repositories(repositories_after),
    }
    prebundle_valid = all(checks.values())
    bundle: dict[str, Any] = {}
    analysis_freeze_id = ""
    if prebundle_valid:
        canonical_outputs = [
            {
                "output_id": row["output_id"],
                "sha256": row["reproduction_1_sha256"],
                "size_bytes": row["reproduction_1_size_bytes"],
            }
            for row in comparison_rows
            if row["canonical"]
        ]
        analysis_freeze_id = runtime.sha256_json(
            {"freeze_plan_id": freeze_plan_id, "canonical_outputs": canonical_outputs}
        )[:16]
        try:
            bundle = build_bundle(
                analysis_freeze_id=analysis_freeze_id,
                freeze_plan_id=freeze_plan_id,
                study_manifest_path=study_manifest_path,
                study=study,
                artifact_snapshot=artifact_snapshot,
                repository_snapshot=repository_snapshot,
                software_snapshot=software_snapshot,
                comparison_rows=comparison_rows,
                run_path=run_path,
                comparison_path=comparison_path,
                output_root=output_root,
            )
        except (OSError, ValueError) as error:
            postflight_errors["bundle"] = str(error)
            analysis_freeze_id = ""
    checks["bundle_archived_when_eligible"] = not prebundle_valid or bool(bundle)
    valid = all(checks.values())

    validation = {
        **stable,
        "freeze_plan_id": freeze_plan_id,
        "analysis_freeze_id": analysis_freeze_id,
        "generated_at": utc_now(),
        "valid": valid,
        "frozen": valid,
        "checks": checks,
        "postflight_errors": postflight_errors,
        "counts": {
            "reproductions": study["repetitions"],
            "analysis_steps": len(study["steps"]),
            "run_rows": len(run_rows),
            "declared_outputs": expected_output_count,
            "canonical_outputs": sum(
                output["kind"] in CANONICAL_KINDS
                for step in study["steps"]
                for output in step["outputs"]
            ),
            "matching_canonical_outputs": sum(
                row["canonical"] and row["byte_identical"] for row in comparison_rows
            ),
            "source_artifacts": len(artifact_snapshot),
            "source_files": sum(
                artifact["file_count"] for artifact in artifact_snapshot
            ),
            "software_files": len(software_snapshot),
            "repositories": len(repository_snapshot),
        },
        "outputs": {
            "analysis_runs": runtime.file_identity(run_path),
            "result_hashes": runtime.file_identity(comparison_path),
            **bundle,
        },
        "tools": {"analysis_freezer": tool, "paper_runtime": runtime_tool},
    }
    validation_path = output_root / "validation" / "analysis_freeze.json"
    runtime.write_json(validation_path, validation)
    if not valid:
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise ValueError(f"analysis freeze validation failed: {failed}")

    manifest_path = output_root / "manifests" / "paper_results.json"
    manifest = {
        **stable,
        "freeze_plan_id": freeze_plan_id,
        "analysis_freeze_id": analysis_freeze_id,
        "generated_at": utc_now(),
        "valid": True,
        "frozen": True,
        "counts": validation["counts"],
        "inputs": {
            "study_manifest": runtime.file_identity(study_manifest_path),
            "artifact_manifest": bundle["artifact_manifest"],
            "repository_manifest": bundle["repository_manifest"],
            "software_manifest": bundle["software_manifest"],
        },
        "outputs": {
            "analysis_runs": runtime.file_identity(run_path),
            "result_hashes": runtime.file_identity(comparison_path),
            "result_bundle": bundle["result_bundle"],
            "archived_analysis_runs": bundle["archived_analysis_runs"],
            "archived_result_hashes": bundle["archived_result_hashes"],
            "log_manifest": bundle["log_manifest"],
            "bundle_manifest": bundle["bundle_manifest"],
            "validation": runtime.file_identity(validation_path),
        },
        "tools": {"analysis_freezer": tool, "paper_runtime": runtime_tool},
        "manifest_path": str(manifest_path),
    }
    runtime.write_json(manifest_path, manifest)
    return manifest_path


def execute_reproduction(
    *,
    reproduction: int,
    freeze_plan_id: str,
    study: dict[str, Any],
    output_root: Path,
    run_rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    analysis_root = output_root / "reproductions" / f"reproduction_{reproduction}"
    analysis_root.mkdir(parents=True, exist_ok=False)
    status_by_step: dict[str, str] = {}
    outputs: dict[str, dict[str, Any]] = {}
    for step in study["steps"]:
        step_id = step["id"]
        failed_dependencies = [
            dependency
            for dependency in step["depends_on"]
            if status_by_step.get(dependency) != "completed"
        ]
        if failed_dependencies:
            status_by_step[step_id] = "skipped_dependency"
            run_rows.append(
                empty_run_row(
                    freeze_plan_id=freeze_plan_id,
                    reproduction=reproduction,
                    step=step,
                    status="skipped_dependency",
                    message="failed dependencies: " + ", ".join(failed_dependencies),
                )
            )
            continue

        step_root = analysis_root / step_id
        command = expand_command(
            step["command"],
            analysis_root=analysis_root,
            step_root=step_root,
            reproduction=reproduction,
        )
        stdout_path = (
            output_root
            / "logs"
            / f"reproduction_{reproduction}"
            / f"{step_id}.stdout.txt"
        )
        stderr_path = (
            output_root
            / "logs"
            / f"reproduction_{reproduction}"
            / f"{step_id}.stderr.txt"
        )
        try:
            measurement = runtime.run_measured_command(
                command,
                stdout_path=stdout_path,
                stderr_path=stderr_path,
                timeout_seconds=step["timeout_seconds"],
                cwd=step["working_directory"],
                env=study["environment"],
            )
        except OSError as error:
            status_by_step[step_id] = "spawn_failure"
            run_rows.append(
                empty_run_row(
                    freeze_plan_id=freeze_plan_id,
                    reproduction=reproduction,
                    step=step,
                    status="spawn_failure",
                    message=str(error),
                    command=command,
                )
            )
            continue

        status = "completed"
        message = ""
        if measurement["timed_out"]:
            status = "timed_out"
            message = "analysis step exceeded its frozen timeout"
        elif measurement["exit_code"] != 0:
            status = "process_failure"
            message = stderr_path.read_text(encoding="utf-8", errors="replace").strip()
        output_file_count = 0
        output_tree_sha256 = ""
        step_outputs: dict[str, dict[str, Any]] = {}
        if status == "completed":
            try:
                step_outputs, inventory = inspect_step_outputs(
                    step=step,
                    step_root=step_root,
                    result_directories=study["result_directory_names"],
                )
                output_file_count = len(inventory)
                output_tree_sha256 = runtime.sha256_json(inventory)
            except ValueError as error:
                status = "invalid_outputs"
                message = str(error)
        status_by_step[step_id] = status
        if status == "completed":
            outputs.update(step_outputs)
        run_rows.append(
            measured_run_row(
                freeze_plan_id=freeze_plan_id,
                reproduction=reproduction,
                step=step,
                status=status,
                message=message,
                command=command,
                measurement=measurement,
                output_file_count=output_file_count,
                output_tree_sha256=output_tree_sha256,
            )
        )
    return outputs


def inspect_step_outputs(
    *, step: dict[str, Any], step_root: Path, result_directories: list[str]
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    if not step_root.is_dir():
        raise ValueError(f"{step['id']}: analysis step did not create its output root")
    inventory = directory_inventory(step_root)
    declared_paths = {output["path"] for output in step["outputs"]}
    result_paths = {
        row["relative_path"]
        for row in inventory
        if Path(row["relative_path"]).parts[0] in result_directories
    }
    undeclared = sorted(result_paths - declared_paths)
    if undeclared:
        raise ValueError(
            f"{step['id']}: undeclared result files: {', '.join(undeclared)}"
        )
    outputs = {}
    inventory_by_path = {row["relative_path"]: row for row in inventory}
    for output in step["outputs"]:
        relative_path = output["path"]
        identity = inventory_by_path.get(relative_path)
        if identity is None:
            raise ValueError(
                f"{step['id']}: declared output is missing: {relative_path}"
            )
        outputs[output["id"]] = {
            "output_id": output["id"],
            "step_id": step["id"],
            "kind": output["kind"],
            "relative_path": relative_path,
            "path": str((step_root / relative_path).resolve()),
            "sha256": identity["sha256"],
            "size_bytes": identity["size_bytes"],
        }
    return outputs, inventory


def compare_outputs(
    *,
    freeze_plan_id: str,
    study: dict[str, Any],
    outputs_by_reproduction: dict[int, dict[str, dict[str, Any]]],
) -> list[dict[str, Any]]:
    first = outputs_by_reproduction.get(1, {})
    second = outputs_by_reproduction.get(2, {})
    rows = []
    for step in study["steps"]:
        for output in step["outputs"]:
            left = first.get(output["id"], {})
            right = second.get(output["id"], {})
            identical = bool(left and right) and (
                left.get("sha256") == right.get("sha256")
                and left.get("size_bytes") == right.get("size_bytes")
            )
            rows.append(
                {
                    "freeze_plan_id": freeze_plan_id,
                    "output_id": output["id"],
                    "step_id": step["id"],
                    "kind": output["kind"],
                    "relative_path": output["path"],
                    "canonical": output["kind"] in CANONICAL_KINDS,
                    "reproduction_1_path": left.get("path", ""),
                    "reproduction_1_sha256": left.get("sha256", ""),
                    "reproduction_1_size_bytes": left.get("size_bytes", ""),
                    "reproduction_2_path": right.get("path", ""),
                    "reproduction_2_sha256": right.get("sha256", ""),
                    "reproduction_2_size_bytes": right.get("size_bytes", ""),
                    "byte_identical": identical,
                }
            )
    return rows


def build_bundle(
    *,
    analysis_freeze_id: str,
    freeze_plan_id: str,
    study_manifest_path: Path,
    study: dict[str, Any],
    artifact_snapshot: list[dict[str, Any]],
    repository_snapshot: list[dict[str, Any]],
    software_snapshot: list[dict[str, Any]],
    comparison_rows: list[dict[str, Any]],
    run_path: Path,
    comparison_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    bundle_root = output_root / "bundle"
    protocol_path = bundle_root / "protocol" / "study_manifest.json"
    protocol_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(study_manifest_path, protocol_path)

    evidence_root = bundle_root / "evidence"
    archived_run_path = evidence_root / "analysis_runs.csv"
    archived_comparison_path = evidence_root / "result_hashes.csv"
    evidence_root.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(run_path, archived_run_path)
    shutil.copyfile(comparison_path, archived_comparison_path)
    if not same_file_content(archived_run_path, run_path):
        raise ValueError("bundle analysis-run copy differs from source evidence")
    if not same_file_content(archived_comparison_path, comparison_path):
        raise ValueError("bundle result-hash copy differs from source evidence")
    source_log_root = output_root / "logs"
    log_inventory = directory_inventory(source_log_root)
    archived_log_root = bundle_root / "logs"
    shutil.copytree(source_log_root, archived_log_root)
    archived_log_inventory = directory_inventory(archived_log_root)
    if archived_log_inventory != log_inventory:
        raise ValueError("bundle log copy differs from execution logs")

    artifact_path = bundle_root / "manifests" / "artifacts.json"
    repository_path = bundle_root / "manifests" / "repositories.json"
    software_path = bundle_root / "manifests" / "software.json"
    log_path = bundle_root / "manifests" / "logs.json"
    runtime.write_json(
        artifact_path,
        {
            "schema_version": SCHEMA_VERSION,
            "analysis_freeze_id": analysis_freeze_id,
            "artifacts": artifact_snapshot,
        },
    )
    runtime.write_json(
        repository_path,
        {
            "schema_version": SCHEMA_VERSION,
            "analysis_freeze_id": analysis_freeze_id,
            "repositories": repository_snapshot,
        },
    )
    runtime.write_json(
        software_path,
        {
            "schema_version": SCHEMA_VERSION,
            "analysis_freeze_id": analysis_freeze_id,
            "software": software_snapshot,
        },
    )
    runtime.write_json(
        log_path,
        {
            "schema_version": SCHEMA_VERSION,
            "analysis_freeze_id": analysis_freeze_id,
            "file_count": len(archived_log_inventory),
            "tree_sha256": runtime.sha256_json(archived_log_inventory),
            "files": archived_log_inventory,
        },
    )

    frozen_results = []
    for row in comparison_rows:
        source = Path(str(row["reproduction_1_path"]))
        destination = bundle_root / "results" / row["step_id"] / row["relative_path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        identity = runtime.file_identity(destination)
        if (
            identity["sha256"] != row["reproduction_1_sha256"]
            or identity["size_bytes"] != row["reproduction_1_size_bytes"]
        ):
            raise ValueError(f"bundle copy changed result: {row['output_id']}")
        frozen_results.append(
            {
                "output_id": row["output_id"],
                "step_id": row["step_id"],
                "kind": row["kind"],
                "path": str(destination.relative_to(bundle_root)),
                "sha256": identity["sha256"],
                "size_bytes": identity["size_bytes"],
            }
        )
    result_path = bundle_root / "manifests" / "results.json"
    runtime.write_json(
        result_path,
        {
            "schema_version": SCHEMA_VERSION,
            "analysis_freeze_id": analysis_freeze_id,
            "results": frozen_results,
        },
    )
    bundle_stable = {
        "schema_version": SCHEMA_VERSION,
        "analysis_freeze_id": analysis_freeze_id,
        "freeze_plan_id": freeze_plan_id,
        "study_id": study["study_id"],
        "protocol_sha256": runtime.file_sha256(protocol_path),
        "artifact_manifest_sha256": runtime.file_sha256(artifact_path),
        "repository_manifest_sha256": runtime.file_sha256(repository_path),
        "software_manifest_sha256": runtime.file_sha256(software_path),
        "analysis_runs_sha256": runtime.file_sha256(archived_run_path),
        "result_hashes_sha256": runtime.file_sha256(archived_comparison_path),
        "log_manifest_sha256": runtime.file_sha256(log_path),
        "result_manifest_sha256": runtime.file_sha256(result_path),
    }
    bundle_id = runtime.sha256_json(bundle_stable)[:16]
    bundle_manifest_path = bundle_root / "manifests" / "bundle.json"
    runtime.write_json(
        bundle_manifest_path,
        {
            **bundle_stable,
            "bundle_id": bundle_id,
            "generated_at": utc_now(),
            "files": {
                "protocol": bundle_identity(protocol_path, bundle_root),
                "artifacts": bundle_identity(artifact_path, bundle_root),
                "repositories": bundle_identity(repository_path, bundle_root),
                "software": bundle_identity(software_path, bundle_root),
                "analysis_runs": bundle_identity(archived_run_path, bundle_root),
                "result_hashes": bundle_identity(archived_comparison_path, bundle_root),
                "logs": bundle_identity(log_path, bundle_root),
                "results": bundle_identity(result_path, bundle_root),
            },
        },
    )
    return {
        "artifact_manifest": runtime.file_identity(artifact_path),
        "repository_manifest": runtime.file_identity(repository_path),
        "software_manifest": runtime.file_identity(software_path),
        "result_bundle": runtime.file_identity(result_path),
        "archived_analysis_runs": runtime.file_identity(archived_run_path),
        "archived_result_hashes": runtime.file_identity(archived_comparison_path),
        "log_manifest": runtime.file_identity(log_path),
        "bundle_manifest": runtime.file_identity(bundle_manifest_path),
    }


def load_study_manifest(path: Path) -> dict[str, Any]:
    value = runtime.load_json(path)
    required = {
        "schema_version",
        "study_id",
        "repetitions",
        "canonical_output_kinds",
        "result_directory_names",
        "environment",
        "inputs",
        "repositories",
        "software",
        "steps",
    }
    if (
        not isinstance(value, dict)
        or set(value) != required
        or value.get("schema_version") != SCHEMA_VERSION
    ):
        raise ValueError("analysis freeze study manifest fields are invalid")
    study_id = safe_id(value.get("study_id"), "study id")
    if value.get("repetitions") != 2:
        raise ValueError("analysis freeze requires exactly two reproductions")
    if value.get("canonical_output_kinds") != list(CANONICAL_KINDS):
        raise ValueError("analysis freeze canonical output kinds are invalid")
    if value.get("result_directory_names") != list(RESULT_DIRECTORIES):
        raise ValueError("analysis freeze result directory names are invalid")
    environment = normalize_environment(value.get("environment"))
    root = path.parent
    inputs = normalize_path_entries(value.get("inputs"), root, "input")
    repositories = normalize_path_entries(
        value.get("repositories"), root, "repository", allow_empty=True
    )
    repository_ids = {entry["id"] for entry in repositories}
    software = normalize_software(value.get("software"), root, repository_ids)
    steps = normalize_steps(
        value.get("steps"),
        root=root,
        inputs=inputs,
        repositories=repositories,
        software=software,
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "study_id": study_id,
        "repetitions": 2,
        "canonical_output_kinds": list(CANONICAL_KINDS),
        "result_directory_names": list(RESULT_DIRECTORIES),
        "environment": environment,
        "inputs": inputs,
        "repositories": repositories,
        "software": software,
        "steps": steps,
    }


def normalize_environment(value: Any) -> dict[str, str]:
    if not isinstance(value, dict) or not value:
        raise ValueError("analysis freeze environment is missing")
    environment = {}
    for key, item in value.items():
        if (
            not isinstance(key, str)
            or ENVIRONMENT_KEY.fullmatch(key) is None
            or not isinstance(item, str)
            or "\x00" in item
        ):
            raise ValueError("analysis freeze environment fields are invalid")
        environment[key] = item
    required = {
        "LC_ALL": "C",
        "TZ": "UTC",
        "PYTHONHASHSEED": "0",
    }
    if any(environment.get(key) != expected for key, expected in required.items()):
        raise ValueError(
            "analysis freeze environment must set LC_ALL=C, TZ=UTC, "
            "and PYTHONHASHSEED=0"
        )
    return dict(sorted(environment.items()))


def normalize_path_entries(
    value: Any, root: Path, label: str, *, allow_empty: bool = False
) -> list[dict[str, Any]]:
    if not isinstance(value, list) or (not value and not allow_empty):
        raise ValueError(f"analysis freeze {label}s are missing")
    rows = []
    ids = set()
    for entry in value:
        if not isinstance(entry, dict) or set(entry) != {"id", "path"}:
            raise ValueError(f"analysis freeze {label} fields are invalid")
        entry_id = safe_id(entry.get("id"), f"{label} id")
        if entry_id in ids:
            raise ValueError(f"analysis freeze {label} ids are duplicated")
        ids.add(entry_id)
        raw_path = Path(nonempty_string(entry["path"], f"{label} path"))
        if not raw_path.is_absolute():
            raw_path = root / raw_path
        rows.append(
            {
                "id": entry_id,
                "declared_path": str(raw_path.absolute()),
                "path": raw_path.resolve(),
            }
        )
    return rows


def normalize_software(
    value: Any, root: Path, repository_ids: set[str]
) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not value:
        raise ValueError("analysis freeze software entries are missing")
    rows = []
    ids = set()
    for entry in value:
        fields = {"id", "path", "version", "repository_id"}
        if not isinstance(entry, dict) or set(entry) != fields:
            raise ValueError("analysis freeze software fields are invalid")
        software_id = safe_id(entry.get("id"), "software id")
        if software_id in ids:
            raise ValueError("analysis freeze software ids are duplicated")
        ids.add(software_id)
        version = nonempty_string(entry.get("version"), f"{software_id} version")
        repository_id = entry.get("repository_id")
        if repository_id is not None and repository_id not in repository_ids:
            raise ValueError(f"{software_id}: software repository id is unknown")
        rows.append(
            {
                "id": software_id,
                "path": resolve_path(entry["path"], root),
                "version": version,
                "repository_id": repository_id,
            }
        )
    return rows


def normalize_steps(
    value: Any,
    *,
    root: Path,
    inputs: list[dict[str, Any]],
    repositories: list[dict[str, Any]],
    software: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not value:
        raise ValueError("analysis freeze steps are missing")
    input_by_id = {entry["id"]: entry for entry in inputs}
    software_by_id = {entry["id"]: entry for entry in software}
    seen_steps = set()
    seen_outputs = set()
    used_inputs = set()
    used_software = set()
    observed_kinds = set()
    steps = []
    fields = {
        "id",
        "working_directory",
        "timeout_seconds",
        "depends_on",
        "input_ids",
        "software_ids",
        "command",
        "outputs",
    }
    for entry in value:
        if not isinstance(entry, dict) or set(entry) != fields:
            raise ValueError("analysis freeze step fields are invalid")
        step_id = safe_id(entry.get("id"), "step id")
        if step_id in seen_steps:
            raise ValueError("analysis freeze step ids are duplicated")
        dependencies = string_list(entry.get("depends_on"), f"{step_id} dependencies")
        if any(dependency not in seen_steps for dependency in dependencies):
            raise ValueError(f"{step_id}: dependencies must name earlier steps")
        input_ids = string_list(entry.get("input_ids"), f"{step_id} input ids")
        software_ids = string_list(entry.get("software_ids"), f"{step_id} software ids")
        if any(input_id not in input_by_id for input_id in input_ids):
            raise ValueError(f"{step_id}: input id is unknown")
        if any(software_id not in software_by_id for software_id in software_ids):
            raise ValueError(f"{step_id}: software id is unknown")
        if not input_ids and not dependencies:
            raise ValueError(f"{step_id}: step has no source input or dependency")
        workdir = resolve_path(entry.get("working_directory"), root)
        if not workdir.is_dir():
            raise ValueError(f"{step_id}: working directory does not exist: {workdir}")
        command = normalize_command(entry.get("command"), root, step_id)
        if not any("{step_root}" in argument for argument in command):
            raise ValueError(f"{step_id}: command omits {{step_root}}")
        executable = Path(command[0]).resolve()
        executable_ids = {
            software_id
            for software_id in software_ids
            if software_by_id[software_id]["path"] == executable
        }
        if not executable_ids:
            raise ValueError(f"{step_id}: command executable is not declared software")
        command_paths = command_path_candidates(command, workdir)
        for input_id in input_ids:
            input_path = input_by_id[input_id]["path"]
            if not any(
                candidate == input_path or candidate.is_relative_to(input_path)
                for candidate in command_paths
            ):
                raise ValueError(f"{step_id}: command omits input {input_id}")
        validate_command_paths(
            command=command,
            step_id=step_id,
            inputs=inputs,
            repositories=repositories,
            software=software,
            working_directory=workdir,
        )
        outputs = normalize_outputs(entry.get("outputs"), step_id, seen_outputs)
        observed_kinds.update(output["kind"] for output in outputs)
        timeout = positive_int(entry.get("timeout_seconds"), f"{step_id} timeout")
        steps.append(
            {
                "id": step_id,
                "working_directory": workdir,
                "timeout_seconds": timeout,
                "depends_on": dependencies,
                "input_ids": input_ids,
                "software_ids": software_ids,
                "command": command,
                "outputs": outputs,
            }
        )
        seen_steps.add(step_id)
        used_inputs.update(input_ids)
        used_software.update(software_ids)
    if used_inputs != set(input_by_id):
        raise ValueError("analysis freeze contains unused source inputs")
    if used_software != set(software_by_id):
        raise ValueError("analysis freeze contains unused software entries")
    if not set(CANONICAL_KINDS).issubset(observed_kinds):
        raise ValueError("analysis freeze requires table and figure_data outputs")
    return steps


def normalize_command(value: Any, root: Path, step_id: str) -> list[str]:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(argument, str) or not argument for argument in value)
    ):
        raise ValueError(f"{step_id}: command is invalid")
    executable = Path(value[0])
    if executable.is_absolute():
        executable = executable.resolve()
    else:
        found = shutil.which(value[0])
        if found is None:
            candidate = (root / executable).resolve()
            if not candidate.is_file():
                raise ValueError(f"{step_id}: command executable is unavailable")
            executable = candidate
        else:
            executable = Path(found).resolve()
    if not executable.is_file():
        raise ValueError(f"{step_id}: command executable is unavailable")
    command = [str(executable), *value[1:]]
    if executable.name in {"sh", "bash", "zsh", "dash"} and "-c" in command[1:]:
        raise ValueError(f"{step_id}: shell command strings are forbidden")
    return command


def validate_command_paths(
    *,
    command: list[str],
    step_id: str,
    inputs: list[dict[str, Any]],
    repositories: list[dict[str, Any]],
    software: list[dict[str, Any]],
    working_directory: Path,
) -> None:
    input_paths = [entry["path"] for entry in inputs]
    software_paths = {entry["path"] for entry in software}
    for candidate in command_path_candidates(command, working_directory):
        if candidate in software_paths or any(
            candidate == path or candidate.is_relative_to(path) for path in input_paths
        ):
            continue
        if any(
            tracked_repository_path(candidate, repository)
            for repository in repositories
        ):
            continue
        raise ValueError(f"{step_id}: command uses undeclared path: {candidate}")


def command_path_candidates(command: list[str], working_directory: Path) -> list[Path]:
    candidates = []
    for argument in command:
        if "{analysis_root}" in argument or "{step_root}" in argument:
            continue
        candidate = Path(argument.split("=", 1)[-1])
        if not candidate.is_absolute():
            candidate = working_directory / candidate
            if not candidate.exists():
                continue
        candidates.append(candidate.resolve())
    return candidates


def tracked_repository_path(candidate: Path, repository: dict[str, Any]) -> bool:
    root = repository["path"]
    if not candidate.is_file() or not candidate.is_relative_to(root):
        return False
    relative_path = candidate.relative_to(root).as_posix()
    completed = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--error-unmatch", relative_path],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.returncode == 0


def normalize_outputs(
    value: Any, step_id: str, seen_outputs: set[str]
) -> list[dict[str, str]]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{step_id}: outputs are missing")
    outputs = []
    paths = set()
    for entry in value:
        if not isinstance(entry, dict) or set(entry) != {"id", "path", "kind"}:
            raise ValueError(f"{step_id}: output fields are invalid")
        output_id = safe_id(entry.get("id"), "output id")
        if output_id in seen_outputs:
            raise ValueError("analysis freeze output ids are duplicated")
        relative_path = safe_relative_path(entry.get("path"), f"{output_id} path")
        if relative_path in paths:
            raise ValueError(f"{step_id}: output paths are duplicated")
        if Path(relative_path).parts[0] not in RESULT_DIRECTORIES:
            raise ValueError(f"{output_id}: output is outside a result directory")
        kind = entry.get("kind")
        if kind not in OUTPUT_KINDS:
            raise ValueError(f"{output_id}: output kind is invalid")
        outputs.append({"id": output_id, "path": relative_path, "kind": kind})
        seen_outputs.add(output_id)
        paths.add(relative_path)
    return outputs


def snapshot_artifacts(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [snapshot_artifact(entry) for entry in entries]


def snapshot_artifact(entry: dict[str, Any]) -> dict[str, Any]:
    path = entry["path"]
    declared_path = Path(entry["declared_path"])
    if declared_path.is_symlink() or not path.exists():
        raise ValueError(f"source artifact is missing or symbolic: {path}")
    if path.is_file():
        files = [
            {
                "relative_path": ".",
                "sha256": runtime.file_sha256(path),
                "size_bytes": path.stat().st_size,
            }
        ]
        kind = "file"
    elif path.is_dir():
        files = directory_inventory(path)
        kind = "directory"
    else:
        raise ValueError(f"source artifact is not a file or directory: {path}")
    if not files:
        raise ValueError(f"source artifact has no files: {path}")
    return {
        "id": entry["id"],
        "path": str(path),
        "kind": kind,
        "file_count": len(files),
        "size_bytes": sum(row["size_bytes"] for row in files),
        "tree_sha256": runtime.sha256_json(files),
        "files": files,
    }


def directory_inventory(root: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"artifact tree contains a symbolic link: {path}")
        if path.is_file():
            rows.append(
                {
                    "relative_path": path.relative_to(root).as_posix(),
                    "sha256": runtime.file_sha256(path),
                    "size_bytes": path.stat().st_size,
                }
            )
    return rows


def internal_software_entries(
    study: dict[str, Any], tool: dict[str, Any], runtime_tool: dict[str, Any]
) -> list[dict[str, Any]]:
    declared_ids = {entry["id"] for entry in study["software"]}
    internal = []
    for software_id, identity in (
        ("internal-analysis-freezer", tool),
        ("internal-paper-runtime", runtime_tool),
    ):
        if software_id in declared_ids:
            raise ValueError(
                f"reserved internal software id is declared: {software_id}"
            )
        path = Path(identity["path"])
        repository_id = next(
            (
                repository["id"]
                for repository in study["repositories"]
                if path.is_relative_to(repository["path"])
            ),
            None,
        )
        internal.append(
            {
                "id": software_id,
                "path": path,
                "version": identity["version"],
                "repository_id": repository_id,
            }
        )
    return internal


def snapshot_software(
    entries: list[dict[str, Any]], repositories: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    repository_by_id = {entry["id"]: entry for entry in repositories}
    rows = []
    for entry in entries:
        path = entry["path"]
        if not path.is_file():
            raise ValueError(f"software file is unavailable: {path}")
        repository_id = entry["repository_id"]
        if repository_id is not None:
            repository_path = Path(repository_by_id[repository_id]["path"])
            if not path.is_relative_to(repository_path):
                raise ValueError(
                    f"{entry['id']}: software is outside its declared repository"
                )
            if not tracked_repository_path(path, {"path": repository_path}):
                raise ValueError(
                    f"{entry['id']}: repository software is not tracked by Git"
                )
        rows.append(
            {
                "id": entry["id"],
                "path": str(path),
                "version": entry["version"],
                "repository_id": repository_id,
                "sha256": runtime.file_sha256(path),
                "size_bytes": path.stat().st_size,
            }
        )
    return rows


def snapshot_repositories(
    entries: list[dict[str, Any]], *, require_clean: bool = True
) -> list[dict[str, Any]]:
    rows = []
    for entry in entries:
        path = entry["path"]
        if not path.is_dir():
            raise ValueError(f"repository is unavailable: {path}")
        top = git_command(path, "rev-parse", "--show-toplevel")
        if Path(top).resolve() != path:
            raise ValueError(f"repository path is not a git root: {path}")
        tracked_status = git_command(
            path, "status", "--porcelain", "--untracked-files=no"
        )
        if tracked_status and require_clean:
            raise ValueError(f"repository has tracked changes: {path}")
        untracked = git_command(path, "status", "--porcelain", "--untracked-files=all")
        remote = git_command(path, "remote", "get-url", "origin", required=False)
        rows.append(
            {
                "id": entry["id"],
                "path": str(path),
                "commit": git_command(path, "rev-parse", "HEAD"),
                "remote": remote,
                "tracked_clean": not tracked_status,
                "untracked_paths": sorted(
                    line[3:]
                    for line in untracked.splitlines()
                    if line.startswith("?? ")
                ),
            }
        )
    return rows


def stable_artifacts(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            key: row[key]
            for key in ("id", "path", "kind", "file_count", "size_bytes", "tree_sha256")
        }
        for row in rows
    ]


def stable_software(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            key: row[key]
            for key in (
                "id",
                "path",
                "version",
                "repository_id",
                "sha256",
                "size_bytes",
            )
        }
        for row in rows
    ]


def stable_repositories(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {key: row[key] for key in ("id", "path", "commit", "remote", "tracked_clean")}
        for row in rows
    ]


def stable_steps(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            **row,
            "working_directory": str(row["working_directory"]),
        }
        for row in rows
    ]


def expand_command(
    command: list[str], *, analysis_root: Path, step_root: Path, reproduction: int
) -> list[str]:
    replacements = {
        "{analysis_root}": str(analysis_root),
        "{step_root}": str(step_root),
        "{reproduction}": str(reproduction),
    }
    expanded = []
    for argument in command:
        value = argument
        for token, replacement in replacements.items():
            value = value.replace(token, replacement)
        expanded.append(value)
    return expanded


def measured_run_row(
    *,
    freeze_plan_id: str,
    reproduction: int,
    step: dict[str, Any],
    status: str,
    message: str,
    command: list[str],
    measurement: dict[str, Any],
    output_file_count: int,
    output_tree_sha256: str,
) -> dict[str, Any]:
    return {
        "freeze_plan_id": freeze_plan_id,
        "reproduction": reproduction,
        "step_id": step["id"],
        "status": status,
        "message": message,
        "command_json": runtime.canonical_json(command),
        "working_directory": str(step["working_directory"]),
        "exit_code": measurement["exit_code"],
        "timed_out": measurement["timed_out"],
        "wall_time_seconds": measurement["wall_time_seconds"],
        "user_cpu_seconds": measurement["user_cpu_seconds"],
        "system_cpu_seconds": measurement["system_cpu_seconds"],
        "peak_resident_memory_bytes": measurement["peak_resident_memory_bytes"],
        "stdout_path": measurement["stdout"]["path"],
        "stdout_sha256": measurement["stdout"]["sha256"],
        "stderr_path": measurement["stderr"]["path"],
        "stderr_sha256": measurement["stderr"]["sha256"],
        "output_file_count": output_file_count,
        "output_tree_sha256": output_tree_sha256,
    }


def empty_run_row(
    *,
    freeze_plan_id: str,
    reproduction: int,
    step: dict[str, Any],
    status: str,
    message: str,
    command: list[str] | None = None,
) -> dict[str, Any]:
    return {
        field: value
        for field, value in {
            "freeze_plan_id": freeze_plan_id,
            "reproduction": reproduction,
            "step_id": step["id"],
            "status": status,
            "message": message,
            "command_json": runtime.canonical_json(command or step["command"]),
            "working_directory": str(step["working_directory"]),
        }.items()
    }


def bundle_identity(path: Path, bundle_root: Path) -> dict[str, Any]:
    return {
        "path": path.relative_to(bundle_root).as_posix(),
        "sha256": runtime.file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def same_file_content(left: Path, right: Path) -> bool:
    return left.stat().st_size == right.stat().st_size and runtime.file_sha256(
        left
    ) == runtime.file_sha256(right)


def ensure_output_separate(output_root: Path, study: dict[str, Any]) -> None:
    for entry in study["inputs"]:
        path = entry["path"]
        if output_root == path or output_root.is_relative_to(path):
            raise ValueError(
                f"analysis freeze output root is inside immutable input: {path}"
            )


def refuse_output_root(path: Path) -> None:
    if path.exists():
        raise ValueError(f"refusing to overwrite existing output root: {path}")
    path.mkdir(parents=True)


def resolve_path(value: Any, root: Path) -> Path:
    text = nonempty_string(value, "path")
    path = Path(text)
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def safe_relative_path(value: Any, label: str) -> str:
    text = nonempty_string(value, label)
    path = Path(text)
    if path.is_absolute() or ".." in path.parts or text in {".", ".."}:
        raise ValueError(f"{label} must be a safe relative path")
    return path.as_posix()


def safe_id(value: Any, label: str) -> str:
    text = nonempty_string(value, label)
    if text in {".", ".."} or SAFE_ID.fullmatch(text) is None:
        raise ValueError(f"{label} contains unsafe characters")
    return text


def string_list(value: Any, label: str) -> list[str]:
    if (
        not isinstance(value, list)
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise ValueError(f"{label} must be a unique string list")
    return list(value)


def nonempty_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty string")
    return value.strip()


def positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def git_command(root: Path, *args: str, required: bool = True) -> str:
    completed = subprocess.run(
        ["git", "-C", str(root), *args],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        if not required:
            return ""
        raise ValueError(
            f"git {' '.join(args)} failed for {root}: {completed.stderr.strip()}"
        )
    return completed.stdout.strip()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
