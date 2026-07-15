"""Shared runtime and report utilities for seqcheck paper experiments."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


RUNTIME_SCHEMA_VERSION = "0.1.0"


def run_measured_command(
    argv: list[str],
    *,
    stdout_path: Path,
    stderr_path: Path,
    timeout_seconds: int,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
) -> dict[str, Any]:
    if not argv:
        raise ValueError("measured command cannot be empty")
    if timeout_seconds <= 0:
        raise ValueError("measured command timeout must be positive")
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    timed_out = False
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(argv, stdout=stdout, stderr=stderr, cwd=cwd, env=env)
        while True:
            waited_pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if waited_pid == process.pid:
                break
            if time.perf_counter() - started >= timeout_seconds:
                timed_out = True
                process.kill()
                waited_pid, status, usage = os.wait4(process.pid, 0)
                if waited_pid != process.pid:
                    raise RuntimeError("failed to reap timed-out measured command")
                break
            time.sleep(0.01)
        process.returncode = os.waitstatus_to_exitcode(status)
    wall_time = time.perf_counter() - started
    return {
        "runtime_schema_version": RUNTIME_SCHEMA_VERSION,
        "argv": argv,
        "exit_code": process.returncode,
        "timed_out": timed_out,
        "wall_time_seconds": wall_time,
        "user_cpu_seconds": usage.ru_utime,
        "system_cpu_seconds": usage.ru_stime,
        "peak_resident_memory_bytes": max_rss_bytes(usage.ru_maxrss),
        "stdout": file_identity(stdout_path),
        "stderr": file_identity(stderr_path),
    }


def require_success(measurement: dict[str, Any]) -> None:
    if measurement["exit_code"] == 0 and not measurement["timed_out"]:
        return
    stderr_path = Path(measurement["stderr"]["path"])
    message = stderr_path.read_text(encoding="utf-8", errors="replace").strip()
    if measurement["timed_out"]:
        raise ValueError(
            f"command exceeded its timeout: {' '.join(measurement['argv'])}"
        )
    raise ValueError(
        f"command exited {measurement['exit_code']}: {' '.join(measurement['argv'])}"
        + (f": {message}" if message else "")
    )


def flatten_report_metrics(
    payload: dict[str, Any],
    context: dict[str, Any],
    *,
    region_annotations: dict[str, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    rows = []
    for result_index, result in enumerate(payload.get("results", [])):
        if not isinstance(result, dict):
            continue
        ontology_terms = ";".join(sorted(extract_ontology_terms(result)))
        annotations = resolve_result_region_annotations(result, region_annotations)
        for side in ("expected", "observed"):
            values = result.get(side, [])
            if not isinstance(values, list):
                continue
            for metric in values:
                if not isinstance(metric, dict):
                    continue
                data = metric.get("data", {})
                if not isinstance(data, dict):
                    data = {}
                rows.append(
                    {
                        **context,
                        "result_index": result_index,
                        "check": str(result.get("check", "")),
                        "files": join_values(result.get("files", [])),
                        "reads": join_values(result.get("reads", [])),
                        "regions": join_values(result.get("regions", [])),
                        "ontology_terms": ontology_terms,
                        "sequence_types": ";".join(
                            sorted({value["sequence_type"] for value in annotations})
                        ),
                        "region_annotations_json": canonical_json(annotations),
                        "metric_side": side,
                        "metric_id": str(metric.get("id", "")),
                        "metric_name": str(metric.get("name", "")),
                        "metric_description": str(metric.get("description", "")),
                        "data_kind": str(data.get("kind", "")),
                        "unit": str(data.get("unit", "")),
                        "value_json": canonical_json(data.get("value")),
                    }
                )
    return rows


def count_report_metrics(payload: dict[str, Any]) -> int:
    return sum(
        len(result.get(side, []))
        for result in payload.get("results", [])
        if isinstance(result, dict)
        for side in ("expected", "observed")
        if isinstance(result.get(side, []), list)
    )


def flatten_report_assessments(
    payload: dict[str, Any],
    context: dict[str, Any],
    *,
    region_annotations: dict[str, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    rows = []
    for result_index, result in enumerate(payload.get("results", [])):
        if not isinstance(result, dict):
            continue
        ontology_terms = ";".join(sorted(extract_ontology_terms(result)))
        annotations = resolve_result_region_annotations(result, region_annotations)
        assessments = result.get("assessment", [])
        if not isinstance(assessments, list):
            continue
        for assessment_index, assessment in enumerate(assessments):
            if not isinstance(assessment, dict):
                continue
            rows.append(
                {
                    **context,
                    "result_index": result_index,
                    "assessment_index": assessment_index,
                    "check": str(result.get("check", "")),
                    "files": join_values(result.get("files", [])),
                    "reads": join_values(result.get("reads", [])),
                    "regions": join_values(result.get("regions", [])),
                    "ontology_terms": ontology_terms,
                    "sequence_types": ";".join(
                        sorted({value["sequence_type"] for value in annotations})
                    ),
                    "region_annotations_json": canonical_json(annotations),
                    "assessment_type": str(assessment.get("type", "")),
                    "assessment_code": str(assessment.get("code", "")),
                    "assessment_description": str(assessment.get("description", "")),
                    "expected_metric_ids": join_values(
                        assessment.get("expected_ids", [])
                    ),
                    "observed_metric_ids": join_values(
                        assessment.get("observed_ids", [])
                    ),
                }
            )
    return rows


def extract_ontology_terms(value: Any) -> set[str]:
    terms = set()
    if isinstance(value, dict):
        for child in value.values():
            terms.update(extract_ontology_terms(child))
    elif isinstance(value, list):
        for child in value:
            terms.update(extract_ontology_terms(child))
    elif isinstance(value, str) and value.startswith("RGN:"):
        terms.add(value)
    return terms


def index_seqspec_regions(
    payload: dict[str, Any], *, modality: str | None = None
) -> dict[str, dict[str, Any]]:
    """Index nominal region annotations without depending on free-text labels."""
    roots = payload.get("library_spec", [])
    if not isinstance(roots, list):
        raise ValueError("seqspec library_spec is not a list")
    if modality is not None:
        roots = [
            value
            for value in roots
            if isinstance(value, dict) and value.get("region_id") == modality
        ]
        if len(roots) != 1:
            raise ValueError(
                f"seqspec has no unique library root for modality: {modality}"
            )
    indexed = {}
    stack = list(reversed(roots))
    while stack:
        region = stack.pop()
        if not isinstance(region, dict):
            raise ValueError("seqspec library_spec contains a non-object region")
        region_id = str(region.get("region_id", "")).strip()
        sequence_type = str(region.get("sequence_type", "")).strip()
        if not region_id or not sequence_type:
            raise ValueError("seqspec region annotation is incomplete")
        if region_id in indexed:
            raise ValueError(f"seqspec region_id is not unique: {region_id}")
        indexed[region_id] = {
            "region_id": region_id,
            "sequence_type": sequence_type,
            "ontology_terms": region_ontology_terms(region.get("region_type")),
        }
        children = region.get("regions", [])
        if not isinstance(children, list):
            raise ValueError(f"seqspec region {region_id} has invalid children")
        stack.extend(reversed(children))
    return indexed


def region_ontology_terms(value: Any) -> list[str]:
    values = [value] if isinstance(value, str) else value
    if not isinstance(values, list) or not values:
        raise ValueError("seqspec region_type is not a string or nonempty list")
    terms = []
    for term in values:
        if not isinstance(term, str) or not term.startswith("RGN:"):
            raise ValueError(f"seqspec region_type is not an ontology term: {term}")
        terms.append(term)
    if len(terms) != len(set(terms)):
        raise ValueError("seqspec region_type contains duplicate ontology terms")
    return sorted(terms)


def resolve_result_region_annotations(
    result: dict[str, Any],
    region_annotations: dict[str, dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    if region_annotations is None:
        return []
    region_ids = result.get("regions", [])
    if not isinstance(region_ids, list):
        raise ValueError("seqcheck result regions is not a list")
    resolved = []
    for value in region_ids:
        region_id = str(value).strip()
        if not region_id or region_id not in region_annotations:
            raise ValueError(f"seqcheck result references unknown region: {region_id}")
        resolved.append(region_annotations[region_id])
    return resolved


def join_values(value: Any) -> str:
    if not isinstance(value, list):
        return ""
    return ";".join(str(item) for item in value)


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def max_rss_bytes(value: int) -> int:
    return value if sys.platform == "darwin" else value * 1024


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


def script_identity(path: Path, *, version: str) -> dict[str, Any]:
    path = path.resolve()
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
        "version": version,
        "path": str(path),
        "sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
        "git_commit": git_output(root, "rev-parse", "HEAD"),
        "git_dirty": bool(status),
        "python": sys.version.split()[0],
    }


def executable_identity(path: Path, *, timeout_seconds: int = 30) -> dict[str, Any]:
    path = path.resolve()
    completed = subprocess.run(
        [str(path), "--version"],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    version = completed.stdout.strip() or completed.stderr.strip()
    if completed.returncode != 0 or not version:
        raise ValueError(f"could not identify executable: {path}")
    return {
        "path": str(path),
        "sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
        "version": version,
    }


def functional_script_identity(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "version": value["version"],
        "sha256": value["sha256"],
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
