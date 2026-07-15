#!/usr/bin/env python3
"""Prepare and verify the deterministic IGVF audit repeatability subset."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import igvf_audit as audit
    import paper_runtime as runtime
except ModuleNotFoundError:
    from scripts import igvf_audit as audit
    from scripts import paper_runtime as runtime


SCHEMA_VERSION = "0.1.0"
TOOL_VERSION = "0.1.0"
SELECTION_FIELDS = (
    "selection_id",
    "selection_rank",
    "configuration_accession",
    "modality",
    "access_class",
    "primary_core_sha256",
)
COMPARISON_FIELDS = (
    "verification_id",
    "selection_id",
    "selection_rank",
    "configuration_accession",
    "modality",
    "primary_core_sha256",
    "repeat_core_sha256",
    "primary_unchanged",
    "repeat_matches_primary",
    "status",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare or verify the frozen IGVF audit repeatability subset."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--primary-root", required=True, type=Path)
    prepare.add_argument("--audit-protocol", required=True, type=Path)
    prepare.add_argument("--output-root", required=True, type=Path)

    verify = subparsers.add_parser("verify")
    verify.add_argument("--primary-root", required=True, type=Path)
    verify.add_argument("--repeat-root", required=True, type=Path)
    verify.add_argument("--selection-manifest", required=True, type=Path)
    verify.add_argument("--audit-protocol", required=True, type=Path)
    verify.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.command == "prepare":
            path = prepare_repeatability_subset(
                primary_root=args.primary_root,
                protocol_path=args.audit_protocol,
                output_root=args.output_root,
            )
        else:
            path = verify_repeatability(
                primary_root=args.primary_root,
                repeat_root=args.repeat_root,
                selection_path=args.selection_manifest,
                protocol_path=args.audit_protocol,
                output_root=args.output_root,
            )
    except (OSError, ValueError) as err:
        print(f"verify_igvf_audit_repeatability: {err}", file=sys.stderr)
        return 1
    print(path)
    return 0


def prepare_repeatability_subset(
    *,
    primary_root: Path,
    protocol_path: Path,
    output_root: Path,
) -> Path:
    primary_root = primary_root.resolve()
    protocol_path = protocol_path.resolve()
    output_root = prepare_output_root(output_root)
    protocol = audit.load_audit_protocol(protocol_path)
    tool = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    study, reconciliation, runs = load_audit_root(primary_root)
    repeatability = protocol["repeatability"]
    requested = repeatability["subset_run_count"]
    if len(runs) < requested:
        raise ValueError(
            f"primary audit has {len(runs)} completed runs; {requested} are required"
        )

    ranked = []
    for row in runs:
        key = (row["configuration_accession"], row["modality"])
        report = load_report(primary_root, *key)
        core_sha256 = canonical_report_sha256(
            report,
            repeatability["canonical_excluded_top_level_fields"],
        )
        ranked.append(
            {
                "selection_rank_sha256": runtime.sha256_json(
                    [repeatability["selection_seed"], *key]
                ),
                "configuration_accession": key[0],
                "modality": key[1],
                "access_class": row["access_class"],
                "primary_core_sha256": core_sha256,
            }
        )
    ranked.sort(
        key=lambda row: (
            row["selection_rank_sha256"],
            row["configuration_accession"],
            row["modality"],
        )
    )
    selected = ranked[:requested]
    for rank, row in enumerate(selected, start=1):
        row["selection_rank"] = rank

    stable = {
        "schema_version": SCHEMA_VERSION,
        "tool": runtime.functional_script_identity(tool),
        "protocol": {
            "audit_protocol_id": audit.audit_protocol_id(protocol),
            "stable_sha256": audit.sha256_json(protocol),
        },
        "primary": {
            "study_run_id": study["run_id"],
            "functional_study_sha256": runtime.sha256_json(
                functional_study_identity(study)
            ),
            "runs_sha256": runtime.file_sha256(primary_root / "runs.csv"),
            "reconciliation_sha256": runtime.file_sha256(
                primary_root / "validation" / "reconciliation.json"
            ),
        },
        "selection_seed": repeatability["selection_seed"],
        "requested_run_count": requested,
        "canonical_excluded_top_level_fields": repeatability[
            "canonical_excluded_top_level_fields"
        ],
        "selected_runs": selected,
    }
    selection_id = runtime.sha256_json(stable)[:16]
    table_rows = [
        {
            "selection_id": selection_id,
            "selection_rank": row["selection_rank"],
            "configuration_accession": row["configuration_accession"],
            "modality": row["modality"],
            "access_class": row["access_class"],
            "primary_core_sha256": row["primary_core_sha256"],
        }
        for row in selected
    ]
    table_path = output_root / "tables" / "repeat_selection.csv"
    write_csv(table_path, table_rows, SELECTION_FIELDS)
    configuration_path = output_root / "inputs" / "configuration_accessions.txt"
    configuration_path.parent.mkdir(parents=True, exist_ok=True)
    configuration_path.write_text(
        "".join(
            f"{accession}\n"
            for accession in sorted(
                {row["configuration_accession"] for row in selected}
            )
        ),
        encoding="utf-8",
    )
    manifest = {
        **stable,
        "selection_id": selection_id,
        "generated_at": utc_now(),
        "inputs": {
            "primary_study": runtime.file_identity(
                primary_root / "manifests" / "study.json"
            ),
            "primary_reconciliation": runtime.file_identity(
                primary_root / "validation" / "reconciliation.json"
            ),
            "primary_runs": runtime.file_identity(primary_root / "runs.csv"),
            "audit_protocol": runtime.file_identity(protocol_path),
        },
        "tools": {"repeatability": tool},
        "outputs": {
            "selection_table": runtime.file_identity(table_path),
            "configuration_accessions": runtime.file_identity(configuration_path),
        },
        "primary_reconciliation_valid": reconciliation["valid"],
    }
    manifest_path = output_root / "manifests" / "repeat_selection.json"
    write_json(manifest_path, manifest)
    load_selection_manifest(manifest_path)
    return manifest_path


def verify_repeatability(
    *,
    primary_root: Path,
    repeat_root: Path,
    selection_path: Path,
    protocol_path: Path,
    output_root: Path,
) -> Path:
    primary_root = primary_root.resolve()
    repeat_root = repeat_root.resolve()
    selection_path = selection_path.resolve()
    protocol_path = protocol_path.resolve()
    output_root = prepare_output_root(output_root)
    protocol = audit.load_audit_protocol(protocol_path)
    tool = runtime.script_identity(Path(__file__).resolve(), version=TOOL_VERSION)
    selection = load_selection_manifest(selection_path)
    if selection["protocol"]["audit_protocol_id"] != audit.audit_protocol_id(protocol):
        raise ValueError("repeat selection and audit protocol IDs differ")
    if selection["protocol"]["stable_sha256"] != audit.sha256_json(protocol):
        raise ValueError("repeat selection and audit protocol hashes differ")

    primary_study, _, primary_runs = load_audit_root(primary_root)
    repeat_study, _, repeat_runs = load_audit_root(repeat_root)
    if primary_study["run_id"] != selection["primary"]["study_run_id"]:
        raise ValueError("primary audit study ID differs from repeat selection")
    primary_by_key = index_runs(primary_runs, "primary")
    repeat_by_key = index_runs(repeat_runs, "repeat")

    stable_rows = []
    comparison_rows = []
    excluded = selection["canonical_excluded_top_level_fields"]
    for selected in selection["selected_runs"]:
        key = (selected["configuration_accession"], selected["modality"])
        primary_row = primary_by_key.get(key)
        repeat_row = repeat_by_key.get(key)
        primary_sha = ""
        repeat_sha = ""
        status = "matched"
        if primary_row is None:
            status = "missing_primary_run"
        else:
            primary_sha = canonical_report_sha256(
                load_report(primary_root, *key), excluded
            )
        if repeat_row is None:
            status = "missing_repeat_run"
        else:
            repeat_sha = canonical_report_sha256(
                load_report(repeat_root, *key), excluded
            )
        primary_unchanged = primary_sha == selected["primary_core_sha256"]
        repeat_matches = bool(primary_sha) and primary_sha == repeat_sha
        if status == "matched" and not primary_unchanged:
            status = "primary_changed"
        elif status == "matched" and not repeat_matches:
            status = "report_mismatch"
        stable_rows.append(
            {
                "selection_rank": selected["selection_rank"],
                "configuration_accession": key[0],
                "modality": key[1],
                "primary_core_sha256": primary_sha,
                "repeat_core_sha256": repeat_sha,
                "primary_unchanged": primary_unchanged,
                "repeat_matches_primary": repeat_matches,
                "status": status,
            }
        )

    verification_stable = {
        "schema_version": SCHEMA_VERSION,
        "tool": runtime.functional_script_identity(tool),
        "selection_id": selection["selection_id"],
        "primary_study_run_id": primary_study["run_id"],
        "repeat_study_run_id": repeat_study["run_id"],
        "primary_functional_study_sha256": runtime.sha256_json(
            functional_study_identity(primary_study)
        ),
        "repeat_functional_study_sha256": runtime.sha256_json(
            functional_study_identity(repeat_study)
        ),
        "rows": stable_rows,
    }
    verification_id = runtime.sha256_json(verification_stable)[:16]
    for row in stable_rows:
        comparison_rows.append(
            {
                "verification_id": verification_id,
                "selection_id": selection["selection_id"],
                **row,
            }
        )
    comparison_path = output_root / "tables" / "repeat_comparison.csv"
    write_csv(comparison_path, comparison_rows, COMPARISON_FIELDS)

    requested = selection["requested_run_count"]
    matched = sum(row["status"] == "matched" for row in stable_rows)
    match_fraction = matched / requested if requested else 0.0
    target = protocol["targets"]["required_repeat_match_fraction"]
    checks = {
        "selection_count_reconciles": len(stable_rows) == requested,
        "primary_study_identity_reconciles": verification_stable[
            "primary_functional_study_sha256"
        ]
        == selection["primary"]["functional_study_sha256"],
        "repeat_study_identity_matches": verification_stable[
            "repeat_functional_study_sha256"
        ]
        == verification_stable["primary_functional_study_sha256"],
        "primary_runs_hash_reconciles": runtime.file_sha256(primary_root / "runs.csv")
        == selection["primary"]["runs_sha256"],
        "primary_reconciliation_hash_reconciles": runtime.file_sha256(
            primary_root / "validation" / "reconciliation.json"
        )
        == selection["primary"]["reconciliation_sha256"],
        "all_selected_runs_present": all(
            row["status"] not in {"missing_primary_run", "missing_repeat_run"}
            for row in stable_rows
        ),
        "primary_reports_unchanged": all(
            row["primary_unchanged"] for row in stable_rows
        ),
        "repeat_match_target_met": match_fraction >= target,
    }
    validation = {
        **verification_stable,
        "verification_id": verification_id,
        "generated_at": utc_now(),
        "valid": all(checks.values()),
        "checks": checks,
        "counts": {
            "requested_runs": requested,
            "compared_runs": len(stable_rows),
            "matched_runs": matched,
            "mismatched_runs": len(stable_rows) - matched,
        },
        "match_fraction": match_fraction,
        "required_match_fraction": target,
        "inputs": {
            "selection": runtime.file_identity(selection_path),
            "protocol": runtime.file_identity(protocol_path),
            "primary_study": runtime.file_identity(
                primary_root / "manifests" / "study.json"
            ),
            "repeat_study": runtime.file_identity(
                repeat_root / "manifests" / "study.json"
            ),
        },
        "outputs": {"comparison": runtime.file_identity(comparison_path)},
        "tools": {"repeatability": tool},
    }
    validation_path = output_root / "validation" / "repeatability.json"
    write_json(validation_path, validation)
    return validation_path


def load_audit_root(
    root: Path,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, str]]]:
    study = load_json(root / "manifests" / "study.json", "study manifest")
    reconciliation = load_json(
        root / "validation" / "reconciliation.json", "audit reconciliation"
    )
    if reconciliation.get("valid") is not True:
        raise ValueError(f"audit reconciliation is not valid: {root}")
    if reconciliation.get("study_run_id") != study.get("run_id"):
        raise ValueError(f"study and reconciliation IDs differ: {root}")
    runs = read_csv(root / "runs.csv")
    if not runs:
        raise ValueError(f"audit root contains no completed runs: {root}")
    required_run_fields = {
        "study_run_id",
        "configuration_accession",
        "modality",
        "access_class",
    }
    if any(not required_run_fields.issubset(row) for row in runs):
        raise ValueError(f"audit run table fields are incomplete: {root}")
    if any(row.get("study_run_id") != study.get("run_id") for row in runs):
        raise ValueError(f"run rows do not match the study ID: {root}")
    index_runs(runs, str(root))
    return study, reconciliation, runs


def load_selection_manifest(path: Path) -> dict[str, Any]:
    value = load_json(path, "repeat selection")
    if value.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("unsupported repeat selection schema")
    stable = {
        key: value.get(key)
        for key in (
            "schema_version",
            "tool",
            "protocol",
            "primary",
            "selection_seed",
            "requested_run_count",
            "canonical_excluded_top_level_fields",
            "selected_runs",
        )
    }
    if value.get("selection_id") != runtime.sha256_json(stable)[:16]:
        raise ValueError("repeat selection ID does not reconcile")
    selected = value.get("selected_runs")
    requested = value.get("requested_run_count")
    if (
        not isinstance(requested, int)
        or isinstance(requested, bool)
        or requested <= 0
        or not isinstance(selected, list)
        or len(selected) != requested
    ):
        raise ValueError("repeat selection row count does not reconcile")
    required_selected_fields = {
        "selection_rank_sha256",
        "selection_rank",
        "configuration_accession",
        "modality",
        "access_class",
        "primary_core_sha256",
    }
    if any(
        not isinstance(row, dict) or set(row) != required_selected_fields
        for row in selected
    ):
        raise ValueError("repeat selection fields are invalid")
    keys = [
        (row.get("configuration_accession"), row.get("modality")) for row in selected
    ]
    if (
        any(
            not isinstance(configuration, str)
            or not configuration
            or not isinstance(modality, str)
            or not modality
            for configuration, modality in keys
        )
        or len(keys) != len(set(keys))
        or any(
            not isinstance(row["access_class"], str)
            or row["access_class"] not in audit.ACCESS_CLASSES
            or not isinstance(row["selection_rank_sha256"], str)
            or len(row["selection_rank_sha256"]) != 64
            or not isinstance(row["primary_core_sha256"], str)
            or len(row["primary_core_sha256"]) != 64
            for row in selected
        )
    ):
        raise ValueError("repeat selection keys must be unique objects")
    if [row.get("selection_rank") for row in selected] != list(
        range(1, len(selected) + 1)
    ):
        raise ValueError("repeat selection ranks do not reconcile")
    return value


def index_runs(
    runs: list[dict[str, str]], label: str
) -> dict[tuple[str, str], dict[str, str]]:
    indexed = {}
    for row in runs:
        key = (row.get("configuration_accession", ""), row.get("modality", ""))
        if not all(key) or key in indexed:
            raise ValueError(f"{label} audit run keys are incomplete or duplicated")
        indexed[key] = row
    return indexed


def load_report(
    root: Path, configuration_accession: str, modality: str
) -> dict[str, Any]:
    path = (
        root / "reports" / configuration_accession / f"{audit.slugify(modality)}.json"
    )
    return load_json(path, "audit report")


def canonical_report_sha256(
    report: dict[str, Any], excluded_top_level_fields: list[str]
) -> str:
    core = {
        key: value
        for key, value in report.items()
        if key not in excluded_top_level_fields
    }
    return runtime.sha256_json(core)


def functional_study_identity(study: dict[str, Any]) -> dict[str, Any]:
    tools = study.get("tools", {})
    tools = tools if isinstance(tools, dict) else {}
    runner = tools.get("audit_runner", {})
    seqcheck = tools.get("seqcheck", {})
    seqspec = tools.get("seqspec", {})
    fastqc = tools.get("fastqc", {})
    return {
        "audit_schema_version": study.get("audit_schema_version"),
        "sampling_method": study.get("sampling_method"),
        "sampling_seed": study.get("sampling_seed"),
        "execution": study.get("execution"),
        "query": study.get("query"),
        "auth": study.get("auth"),
        "frozen_inputs": study.get("frozen_inputs"),
        "tools": {
            "audit_runner": select_fields(runner, ("version", "sha256", "python")),
            "seqcheck": select_fields(
                seqcheck,
                (
                    "command",
                    "version",
                    "runtime_source_sha256",
                    "executable_sha256",
                ),
            ),
            "seqspec": select_fields(
                seqspec,
                (
                    "command",
                    "version",
                    "runtime_source_sha256",
                    "executable_sha256",
                ),
            ),
            "fastqc": select_fields(
                fastqc,
                ("command", "version", "executable_sha256", "configuration"),
            ),
        },
    }


def select_fields(value: Any, fields: tuple[str, ...]) -> dict[str, Any]:
    if not isinstance(value, dict):
        return {}
    return {field: value.get(field) for field in fields}


def prepare_output_root(path: Path) -> Path:
    path = path.resolve()
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"output root is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as err:
        raise ValueError(f"{label} is not valid JSON: {path}") from err
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
