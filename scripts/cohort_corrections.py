"""Validate immutable seqspec correction proposals for paper cohorts."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import time
import urllib.parse
from pathlib import Path
from typing import Any, Callable


SCHEMA_VERSION = "0.1.0"
HUNK_HEADER = re.compile(
    r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@"
)
EVIDENCE_FIELDS = (
    "proposed_correction_manifest",
    "proposed_correction_sha256",
)
EFFECTIVE_FIELDS = (
    "correction_applied",
    "correction_manifest",
    "correction_manifest_sha256",
    "effective_spec_path",
    "effective_spec_sha256",
    "effective_hydration_status",
    "effective_normalized_seqspec_version",
    "effective_modality_match_status",
    "effective_structural_check_status",
    "effective_resource_check_status",
    "effective_fastq_mapping_status",
    "effective_expected_fastq_accessions",
    "effective_structure_sha256",
    "effective_deduplication_key",
)


def load_registry(
    path: Path | None,
    selection_id: str,
    candidate_keys: set[tuple[str, str]],
) -> dict[tuple[str, str], dict[str, str]]:
    if path is None:
        return {}
    registry = load_json(path)
    if registry.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("correction registry does not use schema 0.1.0")
    if registry.get("selection_id") != selection_id:
        raise ValueError("correction registry selection id does not match candidates")
    entries = registry.get("corrections")
    if not isinstance(entries, list):
        raise ValueError("correction registry corrections must be a list")

    corrections = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("correction registry entries must be objects")
        key = (
            str(entry.get("family_id", "")).strip(),
            str(entry.get("configuration_accession", "")).strip(),
        )
        if not all(key):
            raise ValueError("correction registry entry has a missing candidate key")
        if key in corrections:
            raise ValueError(f"correction registry has duplicate key: {'/'.join(key)}")
        if key not in candidate_keys:
            raise ValueError(
                f"correction registry key is not a review candidate: {'/'.join(key)}"
            )
        manifest_value = str(entry.get("manifest", "")).strip()
        declared_sha256 = str(entry.get("sha256", "")).strip()
        if not manifest_value or not declared_sha256:
            raise ValueError(
                f"correction registry entry {'/'.join(key)} needs manifest and sha256"
            )
        manifest_path = Path(manifest_value)
        if not manifest_path.is_absolute():
            manifest_path = path.parent / manifest_path
        manifest_path = manifest_path.resolve()
        observed_sha256 = file_sha256(manifest_path)
        if observed_sha256 != declared_sha256:
            raise ValueError(
                f"correction manifest hash changed for {'/'.join(key)}"
            )
        correction = load_json(manifest_path)
        validate_proposal(correction, key, selection_id)
        corrections[key] = {
            "path": str(manifest_path),
            "sha256": observed_sha256,
        }
    return corrections


def validate_proposal(
    correction: dict[str, Any], key: tuple[str, str], selection_id: str
) -> None:
    if correction.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"correction {'/'.join(key)} does not use schema 0.1.0")
    if correction.get("status") != "proposed_unapproved":
        raise ValueError(f"correction {'/'.join(key)} must remain proposed_unapproved")
    observed_key = (
        str(correction.get("family_id", "")).strip(),
        str(correction.get("configuration_accession", "")).strip(),
    )
    if observed_key != key or correction.get("selection_id") != selection_id:
        raise ValueError(f"correction {'/'.join(key)} candidate identity changed")
    if not str(correction.get("rationale", "")).strip():
        raise ValueError(f"correction {'/'.join(key)} is missing a rationale")
    approval = correction.get("approval", {})
    if not isinstance(approval, dict) or any(
        str(value).strip() for value in approval.values()
    ):
        raise ValueError(
            f"correction {'/'.join(key)} must not contain prefilled approval fields"
        )


def effective_candidate_fields(candidate: dict[str, str]) -> dict[str, str]:
    return {
        "correction_applied": "false",
        "correction_manifest": "",
        "correction_manifest_sha256": "",
        "effective_spec_path": candidate.get("normalized_spec_path", ""),
        "effective_spec_sha256": candidate.get("normalized_spec_sha256", ""),
        "effective_hydration_status": candidate.get("hydration_status", ""),
        "effective_normalized_seqspec_version": candidate.get(
            "normalized_seqspec_version", ""
        ),
        "effective_modality_match_status": candidate.get(
            "modality_match_status", ""
        ),
        "effective_structural_check_status": candidate.get(
            "structural_check_status", ""
        ),
        "effective_resource_check_status": candidate.get(
            "resource_check_status", ""
        ),
        "effective_fastq_mapping_status": candidate.get(
            "fastq_mapping_status", ""
        ),
        "effective_expected_fastq_accessions": candidate.get(
            "expected_fastq_accessions", ""
        ),
        "effective_structure_sha256": candidate.get("structure_sha256", ""),
        "effective_deduplication_key": candidate.get("deduplication_key", ""),
    }


def reconcile_evidence(
    review_rows: list[dict[str, str]],
    correction_registry: dict[tuple[str, str], dict[str, str]],
) -> dict[tuple[str, str], list[str]]:
    errors: dict[tuple[str, str], list[str]] = {}
    for row in review_rows:
        key = row_key(row)
        if not all(key):
            continue
        expected = correction_registry.get(key, {})
        observed_path = row.get("proposed_correction_manifest", "").strip()
        observed_sha256 = row.get("proposed_correction_sha256", "").strip()
        problems = []
        if observed_path != expected.get("path", ""):
            problems.append("proposed correction manifest changed after review packaging")
        if observed_sha256 != expected.get("sha256", ""):
            problems.append("proposed correction hash changed after review packaging")
        if problems:
            errors[key] = problems
    return errors


def included_keys(
    review_rows: list[dict[str, str]],
    correction_registry: dict[tuple[str, str], dict[str, str]],
    decision_fn: Callable[[dict[str, str]], tuple[str, list[str]]],
) -> set[tuple[str, str]]:
    keys = set()
    for row in review_rows:
        key = row_key(row)
        if key not in correction_registry:
            continue
        decision, _ = decision_fn(row)
        if decision == "include":
            keys.add(key)
    return keys


def resolve_seqspec_command(
    seqspec_bin: Path | None, candidate_manifest: dict[str, Any]
) -> list[str]:
    if seqspec_bin is not None:
        command = [str(seqspec_bin.resolve())]
    else:
        value = candidate_manifest.get("hydration", {}).get("seqspec", {}).get(
            "command", []
        )
        command = [str(item) for item in value] if isinstance(value, list) else []
    if not command:
        raise ValueError(
            "a seqspec executable is required to validate an included correction"
        )
    executable = Path(command[0])
    if not executable.is_file():
        raise ValueError(f"seqspec executable does not exist: {executable}")
    return command


def command_identity(command: list[str]) -> dict[str, Any]:
    return {
        "command": command,
        "version": run_seqspec_command(command + ["--version"], 30).strip(),
        "executable_sha256": file_sha256(Path(command[0])),
    }


def validate_correction(
    *,
    candidate: dict[str, str],
    registry_entry: dict[str, str],
    seqspec_command: list[str],
    structural_timeout_seconds: int,
    resource_timeout_seconds: int,
    network_attempts: int,
    retry_backoff_seconds: float,
) -> dict[str, str]:
    manifest_path = Path(registry_entry["path"])
    if file_sha256(manifest_path) != registry_entry["sha256"]:
        raise ValueError("correction manifest changed after registry validation")
    correction = load_json(manifest_path)
    key = row_key(candidate)
    validate_proposal(correction, key, candidate.get("selection_id", ""))
    original_path, original_sha256 = resolve_correction_file(
        correction, "original", manifest_path
    )
    corrected_path, corrected_sha256 = resolve_correction_file(
        correction, "corrected", manifest_path
    )
    diff_path, _ = resolve_correction_file(correction, "diff", manifest_path)
    if original_sha256 != candidate.get("normalized_spec_sha256", ""):
        raise ValueError("correction original hash does not match the candidate spec")
    if corrected_sha256 == original_sha256:
        raise ValueError("correction does not change the candidate spec")
    verify_unified_diff(original_path, corrected_path, diff_path)

    version = parse_seqspec_version(
        run_seqspec_command(
            seqspec_command + ["version", str(corrected_path)],
            structural_timeout_seconds,
        )
    )
    if version != "0.5.0":
        raise ValueError(
            f"corrected seqspec must use version 0.5.0, observed {version or '<missing>'}"
        )
    modalities = json.loads(
        run_seqspec_command(
            seqspec_command
            + ["info", "-k", "modalities", "-f", "json", str(corrected_path)],
            structural_timeout_seconds,
        )
    )
    library_spec = json.loads(
        run_seqspec_command(
            seqspec_command
            + ["info", "-k", "library_spec", "-f", "json", str(corrected_path)],
            structural_timeout_seconds,
        )
    )
    sequence_spec = json.loads(
        run_seqspec_command(
            seqspec_command
            + ["info", "-k", "sequence_spec", "-f", "json", str(corrected_path)],
            structural_timeout_seconds,
        )
    )
    if not isinstance(modalities, list):
        raise ValueError("corrected seqspec modalities output is not a list")
    expected_modalities = split_values(candidate.get("family_expected_modalities", ""))
    observed_modalities = {str(value) for value in modalities}
    if expected_modalities and not expected_modalities.intersection(
        observed_modalities
    ):
        raise ValueError("corrected seqspec does not match the candidate family modality")

    run_seqspec_check(
        seqspec_command + ["check", "--skip", "external", str(corrected_path)],
        structural_timeout_seconds,
        1,
        0,
        "structural",
    )
    run_seqspec_check(
        seqspec_command + ["check", str(corrected_path)],
        resource_timeout_seconds,
        network_attempts,
        retry_backoff_seconds,
        "resource",
    )
    expected_fastqs = extract_expected_fastq_accessions(sequence_spec)
    portal_fastqs = split_values(candidate.get("fastq_accessions", ""))
    if expected_fastqs != portal_fastqs:
        missing = ";".join(sorted(portal_fastqs - expected_fastqs)) or "<none>"
        extra = ";".join(sorted(expected_fastqs - portal_fastqs)) or "<none>"
        raise ValueError(
            f"corrected FASTQs do not match portal FASTQs (missing={missing}, extra={extra})"
        )
    structure_sha256 = sha256_json(
        normalized_structure(library_spec, sequence_spec)
    )
    fastq_set_sha256 = candidate.get("fastq_set_sha256", "").strip()
    if not fastq_set_sha256:
        raise ValueError("candidate is missing its portal FASTQ-set hash")
    deduplication_key = sha256_json([structure_sha256, fastq_set_sha256])
    return {
        "correction_applied": "true",
        "correction_manifest": str(manifest_path.resolve()),
        "correction_manifest_sha256": registry_entry["sha256"],
        "effective_spec_path": str(corrected_path.resolve()),
        "effective_spec_sha256": corrected_sha256,
        "effective_hydration_status": "normalized",
        "effective_normalized_seqspec_version": version,
        "effective_modality_match_status": "matched",
        "effective_structural_check_status": "passed",
        "effective_resource_check_status": "passed",
        "effective_fastq_mapping_status": "matched",
        "effective_expected_fastq_accessions": ";".join(sorted(expected_fastqs)),
        "effective_structure_sha256": structure_sha256,
        "effective_deduplication_key": deduplication_key,
    }


def resolve_correction_file(
    correction: dict[str, Any], field: str, manifest_path: Path
) -> tuple[Path, str]:
    value = correction.get(field)
    if not isinstance(value, dict):
        raise ValueError(f"correction {field} identity is missing")
    path_value = str(value.get("path", "")).strip()
    declared_sha256 = str(value.get("sha256", "")).strip()
    if not path_value or not declared_sha256:
        raise ValueError(f"correction {field} identity needs path and sha256")
    path = Path(path_value)
    if not path.is_absolute():
        path = manifest_path.parent / path
    path = path.resolve()
    if file_sha256(path) != declared_sha256:
        raise ValueError(f"correction {field} file hash changed")
    return path, declared_sha256


def verify_unified_diff(original: Path, corrected: Path, diff_path: Path) -> None:
    original_lines = original.read_text(encoding="utf-8").splitlines(keepends=True)
    corrected_lines = corrected.read_text(encoding="utf-8").splitlines(keepends=True)
    diff_lines = diff_path.read_text(encoding="utf-8").splitlines(keepends=True)
    applied = apply_unified_diff(original_lines, diff_lines)
    if applied != corrected_lines:
        raise ValueError("correction diff does not reproduce original-to-corrected changes")


def apply_unified_diff(original: list[str], diff_lines: list[str]) -> list[str]:
    result = []
    original_index = 0
    index = 0
    hunk_count = 0
    while index < len(diff_lines):
        match = HUNK_HEADER.match(diff_lines[index])
        if match is None:
            index += 1
            continue
        hunk_count += 1
        old_start = int(match.group(1))
        old_count = int(match.group(2) or "1")
        new_start = int(match.group(3))
        new_count = int(match.group(4) or "1")
        hunk_start = 0 if old_start == 0 else old_start - 1
        new_hunk_start = 0 if new_start == 0 else new_start - 1
        if hunk_start < original_index or hunk_start > len(original):
            raise ValueError("correction diff has an invalid hunk location")
        result.extend(original[original_index:hunk_start])
        if new_hunk_start != len(result):
            raise ValueError("correction diff has an invalid corrected-file location")
        original_index = hunk_start
        old_used = 0
        new_used = 0
        index += 1
        while index < len(diff_lines) and not diff_lines[index].startswith("@@"):
            line = diff_lines[index]
            if line.startswith("\\ No newline at end of file"):
                index += 1
                continue
            if not line or line[0] not in {" ", "+", "-"}:
                raise ValueError("correction diff contains an invalid hunk line")
            marker, content = line[0], line[1:]
            if marker in {" ", "-"}:
                if original_index >= len(original) or original[original_index] != content:
                    raise ValueError("correction diff does not match the original file")
                original_index += 1
                old_used += 1
            if marker in {" ", "+"}:
                result.append(content)
                new_used += 1
            index += 1
        if old_used != old_count or new_used != new_count:
            raise ValueError("correction diff hunk counts are inconsistent")
    if hunk_count == 0:
        raise ValueError("correction diff has no change hunks")
    result.extend(original[original_index:])
    return result


def run_seqspec_command(argv: list[str], timeout_seconds: int) -> str:
    try:
        result = subprocess.run(
            argv,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        raise ValueError(
            f"seqspec command exceeded {timeout_seconds} seconds: {' '.join(argv)}"
        ) from error
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip()
        raise ValueError(f"{' '.join(argv)} failed: {message}")
    return result.stdout


def run_seqspec_check(
    argv: list[str],
    timeout_seconds: int,
    attempts: int,
    retry_backoff_seconds: float,
    label: str,
) -> None:
    last_message = ""
    for attempt in range(1, attempts + 1):
        try:
            result = subprocess.run(
                argv,
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
            )
            if result.returncode == 0:
                return
            last_message = result.stderr.strip() or result.stdout.strip()
        except subprocess.TimeoutExpired:
            last_message = f"exceeded {timeout_seconds} seconds"
        if attempt < attempts:
            time.sleep(retry_backoff_seconds)
    raise ValueError(f"corrected seqspec {label} check failed: {last_message}")


def parse_seqspec_version(output: str) -> str:
    for line in output.splitlines():
        if line.startswith("seqspec file version:"):
            return line.split(":", 1)[1].strip()
    return ""


def normalized_structure(library_spec: Any, sequence_spec: Any) -> dict[str, Any]:
    reads = []
    if isinstance(sequence_spec, list):
        for value in sequence_spec:
            if isinstance(value, dict):
                reads.append(
                    {key: item for key, item in value.items() if key != "files"}
                )
            else:
                reads.append(value)
    return {"library_spec": library_spec, "sequence_spec": reads}


def extract_expected_fastq_accessions(sequence_spec: Any) -> set[str]:
    accessions = set()
    if not isinstance(sequence_spec, list):
        return accessions
    for read in sequence_spec:
        if not isinstance(read, dict):
            continue
        files = read.get("files", [])
        if not isinstance(files, list):
            continue
        for item in files:
            if not isinstance(item, dict) or not is_fastq_file(item):
                continue
            file_id = str(item.get("file_id", "")).strip()
            url_accession = sequence_accession_from_url(str(item.get("url", "")))
            if url_accession:
                accessions.add(url_accession)
            elif file_id:
                accessions.add(normalize_fastq_accession(file_id))
            else:
                filename = Path(str(item.get("filename", ""))).name
                accessions.add(normalize_fastq_accession(filename))
    return {value for value in accessions if value}


def sequence_accession_from_url(value: str) -> str:
    pieces = [piece for piece in urllib.parse.urlparse(value).path.split("/") if piece]
    for index, piece in enumerate(pieces[:-1]):
        if piece == "sequence-files":
            return normalize_fastq_accession(pieces[index + 1])
    return ""


def normalize_fastq_accession(value: str) -> str:
    normalized = value.strip().rstrip("/").split("/")[-1]
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if normalized.lower().endswith(suffix):
            return normalized[: -len(suffix)]
    return normalized


def is_fastq_file(item: dict[str, Any]) -> bool:
    filetype = str(item.get("filetype", "")).strip().lower().lstrip(".")
    if filetype in {"fastq", "fq", "fastq.gz", "fq.gz"}:
        return True
    filename = str(item.get("filename", "")).strip().lower()
    return filename.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))


def split_values(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def row_key(row: dict[str, str]) -> tuple[str, str]:
    return (
        row.get("family_id", "").strip(),
        row.get("configuration_accession", "").strip(),
    )


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
