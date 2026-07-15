"""Deterministic seqspec and FASTQ mutation helpers for paper experiments."""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import urllib.request
from decimal import Decimal
from pathlib import Path
from typing import Any, Callable, Iterator

try:
    import paper_runtime as runtime
    import sample_fastq as sampler
except ModuleNotFoundError:
    from scripts import paper_runtime as runtime
    from scripts import sample_fastq as sampler


RUNTIME_VERSION = "0.1.0"
FastqRecord = sampler.FastqRecord


def load_seqspec_json(
    path: Path, *, yq_bin: Path, timeout_seconds: int
) -> dict[str, Any]:
    completed = subprocess.run(
        [str(yq_bin), "-o=json", ".", str(path)],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip()
        raise ValueError(f"could not parse seqspec {path}: {message}")
    try:
        value = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise ValueError(f"yq returned invalid JSON for seqspec {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"seqspec root is not an object: {path}")
    return value


def iter_regions(spec: dict[str, Any]) -> Iterator[dict[str, Any]]:
    roots = spec.get("library_spec", [])
    if not isinstance(roots, list):
        return
    stack = list(reversed(roots))
    while stack:
        region = stack.pop()
        if not isinstance(region, dict):
            continue
        yield region
        children = region.get("regions", [])
        if isinstance(children, list):
            stack.extend(reversed(children))


def find_region(spec: dict[str, Any], region_id: str) -> dict[str, Any]:
    matches = [
        region for region in iter_regions(spec) if region.get("region_id") == region_id
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one region '{region_id}', observed {len(matches)}")
    return matches[0]


def find_read(spec: dict[str, Any], read_id: str) -> dict[str, Any]:
    reads = spec.get("sequence_spec", [])
    matches = [
        read
        for read in reads
        if isinstance(read, dict) and read.get("read_id") == read_id
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one read '{read_id}', observed {len(matches)}")
    return matches[0]


def mutate_spec(
    original: dict[str, Any],
    *,
    operator_id: str,
    variant: str,
    target: dict[str, Any],
    observed_lengths: list[int] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    spec = copy.deepcopy(original)
    read_ids = target.get("read_ids", [])
    region_ids = target.get("region_ids", [])
    details = target.get("details", {})
    if operator_id == "S04":
        if not observed_lengths:
            raise ValueError("S04 requires observed FASTQ lengths")
        read = find_read(spec, only(read_ids, "S04 read"))
        original_max = strict_int(read.get("max_len"), "S04 read max_len")
        observed = set(observed_lengths)
        replacement = next(
            (value for value in range(original_max, 0, -1) if value not in observed),
            None,
        )
        if replacement is None:
            raise ValueError("S04 has no schema-valid unobserved read length")
        before = {"min_len": read.get("min_len"), "max_len": read.get("max_len")}
        read["min_len"] = replacement
        read["max_len"] = replacement
        return spec, {
            "before": before,
            "after": {"min_len": replacement, "max_len": replacement},
        }
    if operator_id == "S05":
        if len(region_ids) != 2:
            raise ValueError("S05 requires two target regions")
        delta = parse_shift_variant(variant)
        declared = details.get("deltas", [])
        if delta not in declared:
            raise ValueError(f"S05 shift {delta} was not declared applicable")
        left = find_region(spec, region_ids[0])
        right = find_region(spec, region_ids[1])
        direction = str(details.get("direction", ""))
        before = {
            region_ids[0]: region_geometry(left),
            region_ids[1]: region_geometry(right),
        }
        shift_region_boundary(left, right, direction=direction, delta=delta)
        recompute_joined_regions(spec)
        return spec, {
            "direction": direction,
            "delta": delta,
            "before": before,
            "after": {
                region_ids[0]: region_geometry(left),
                region_ids[1]: region_geometry(right),
            },
        }
    if operator_id == "S06":
        read = find_read(spec, only(read_ids, "S06 read"))
        before = str(read.get("strand", ""))
        after = {"pos": "neg", "neg": "pos"}.get(before)
        if after is None or variant != f"{before}_to_{after}":
            raise ValueError("S06 strand variant does not match the target read")
        read["strand"] = after
        return spec, {"before": before, "after": after}
    if operator_id == "S07":
        read = find_read(spec, only(read_ids, "S07 read"))
        replacement = str(details.get("replacement_primer_id", ""))
        if not replacement or replacement != only(region_ids, "S07 replacement region"):
            raise ValueError("S07 replacement primer is incomplete")
        before = str(read.get("primer_id", ""))
        read["primer_id"] = replacement
        return spec, {"before": before, "after": replacement}
    if operator_id == "S08":
        region = find_region(spec, only(region_ids, "S08 region"))
        before = concrete_dna(region.get("sequence"), "S08 fixed sequence")
        if variant == "deterministic_replacement":
            after = substitute_sequence(before)
        elif variant == "reverse_complement":
            after = reverse_complement(before)
        else:
            raise ValueError(f"unsupported S08 variant: {variant}")
        if after == before:
            raise ValueError("S08 mutation did not change the fixed sequence")
        region["sequence"] = after
        recompute_joined_regions(spec)
        return spec, {"before": before, "after": after}
    if operator_id in {"S09", "S10"}:
        only(region_ids, f"{operator_id} region")
        return spec, {}
    raise ValueError(f"operator does not mutate a seqspec: {operator_id}")


def shift_region_boundary(
    left: dict[str, Any],
    right: dict[str, Any],
    *,
    direction: str,
    delta: int,
) -> None:
    if delta <= 0:
        raise ValueError("boundary shift must be positive")
    left_sequence = str(left.get("sequence", ""))
    right_sequence = str(right.get("sequence", ""))
    if direction == "left_to_right":
        donor, recipient = left, right
        if len(left_sequence) <= delta:
            raise ValueError("left donor sequence is too short for boundary shift")
        moved = left_sequence[-delta:]
        left["sequence"] = left_sequence[:-delta]
        right["sequence"] = moved + right_sequence
    elif direction == "right_to_left":
        donor, recipient = right, left
        if len(right_sequence) <= delta:
            raise ValueError("right donor sequence is too short for boundary shift")
        moved = right_sequence[:delta]
        right["sequence"] = right_sequence[delta:]
        left["sequence"] = left_sequence + moved
    else:
        raise ValueError(f"unsupported boundary direction: {direction}")
    for field in ("min_len", "max_len"):
        donor_value = strict_int(donor.get(field), f"donor {field}")
        recipient_value = strict_int(recipient.get(field), f"recipient {field}")
        if donor_value <= delta:
            raise ValueError(f"boundary shift makes donor {field} nonpositive")
        donor[field] = donor_value - delta
        recipient[field] = recipient_value + delta


def region_geometry(region: dict[str, Any]) -> dict[str, Any]:
    return {
        "sequence": str(region.get("sequence", "")),
        "min_len": strict_int(region.get("min_len"), "region min_len"),
        "max_len": strict_int(region.get("max_len"), "region max_len"),
    }


def recompute_joined_regions(spec: dict[str, Any]) -> None:
    def visit(region: dict[str, Any]) -> None:
        children = [
            value for value in region.get("regions", []) if isinstance(value, dict)
        ]
        for child in children:
            visit(child)
        if children and region.get("sequence_type") == "joined":
            region["sequence"] = "".join(
                str(child.get("sequence", "")) for child in children
            )
            region["min_len"] = sum(
                strict_int(child.get("min_len"), "child min_len") for child in children
            )
            region["max_len"] = sum(
                strict_int(child.get("max_len"), "child max_len") for child in children
            )

    roots = spec.get("library_spec", [])
    if isinstance(roots, list):
        for root in roots:
            if isinstance(root, dict):
                visit(root)


def bundle_local_resources(
    spec: dict[str, Any],
    *,
    source_spec_path: Path,
    output_root: Path,
    omit_onlist_region_ids: set[str] | None = None,
    replacement_onlists: dict[str, str] | None = None,
) -> list[dict[str, Any]]:
    omitted = omit_onlist_region_ids or set()
    replacements = replacement_onlists or {}
    resources = []
    for region in iter_regions(spec):
        region_id = str(region.get("region_id", ""))
        onlist = region.get("onlist")
        if not isinstance(onlist, dict):
            continue
        if region_id in replacements:
            relative = Path("resources") / f"{safe_name(region_id)}-replacement.txt"
            destination = output_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(replacements[region_id], encoding="ascii")
            rewrite_local_resource(onlist, relative, destination)
            resources.append(
                {"kind": "onlist_replacement", **runtime.file_identity(destination)}
            )
            continue
        if region_id in omitted:
            relative = Path("resources") / f"{safe_name(region_id)}-missing.txt"
            rewrite_local_resource(onlist, relative, None)
            resources.append(
                {
                    "kind": "onlist_omission",
                    "region_id": region_id,
                    "path": str(output_root / relative),
                }
            )
            continue
        identity = bundle_declared_resource(
            onlist,
            source_spec_path=source_spec_path,
            output_root=output_root,
            kind="onlist",
        )
        if identity is not None:
            resources.append(identity)
    for read in spec.get("sequence_spec", []):
        if not isinstance(read, dict):
            continue
        for file in read.get("files", []):
            if not isinstance(file, dict):
                continue
            identity = bundle_declared_resource(
                file,
                source_spec_path=source_spec_path,
                output_root=output_root,
                kind="sequence_file",
            )
            if identity is not None:
                resources.append(identity)
    return resources


def bundle_declared_resource(
    value: dict[str, Any],
    *,
    source_spec_path: Path,
    output_root: Path,
    kind: str,
) -> dict[str, Any] | None:
    if str(value.get("urltype", "")) != "local":
        return None
    locator = str(value.get("url", "")).strip()
    if not locator:
        raise ValueError(f"local {kind} has no URL")
    source = Path(locator)
    if not source.is_absolute():
        source = (source_spec_path.parent / source).resolve()
    if not source.is_file():
        raise ValueError(f"local {kind} does not exist: {source}")
    sha256 = runtime.file_sha256(source)
    relative = Path("resources") / f"{sha256[:16]}-{safe_name(source.name)}"
    destination = output_root / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        shutil.copyfile(source, destination)
    rewrite_local_resource(value, relative, destination)
    return {"kind": kind, **runtime.file_identity(destination)}


def rewrite_local_resource(
    value: dict[str, Any], relative: Path, path: Path | None
) -> None:
    locator = relative.as_posix()
    value["urltype"] = "local"
    value["url"] = locator
    value["filename"] = relative.name
    value["file_id"] = str(value.get("file_id") or relative.name)
    value["filetype"] = str(value.get("filetype") or relative.suffix.lstrip("."))
    value["filesize"] = path.stat().st_size if path is not None else 0
    value["md5"] = md5_file(path) if path is not None else ""


def write_spec(path: Path, spec: dict[str, Any]) -> dict[str, Any]:
    runtime.write_json(path, spec)
    return runtime.file_identity(path)


def read_fastq(path: Path) -> list[FastqRecord]:
    digest = hashlib.sha256()
    with sampler.open_source(str(path)) as (reader, _):
        return list(sampler.iter_fastq(reader, str(path), digest))


def observed_lengths(records: list[FastqRecord]) -> list[int]:
    return [len(record[1].rstrip(b"\r\n")) for record in records]


def fastq_record_id(record: FastqRecord) -> str:
    header = record[0].rstrip(b"\r\n")
    if not header.startswith(b"@"):
        raise ValueError("FASTQ header does not start with '@'")
    value = header[1:].split(maxsplit=1)[0].decode("utf-8", errors="strict")
    if not value:
        raise ValueError("FASTQ record has an empty identifier")
    return value


def target_record_count(total_records: int, event_fraction: float) -> int:
    if total_records <= 0:
        raise ValueError("FASTQ must contain at least one record")
    fraction = Decimal(str(event_fraction))
    if not Decimal(0) < fraction < Decimal(1):
        raise ValueError("event fraction must be between zero and one")
    count = fraction * total_records
    integral = count.to_integral_value()
    return int(integral) if count == integral else 0


def select_record_indices(
    records: list[FastqRecord],
    *,
    event_fraction: float,
    seed: int,
    eligible: set[int] | None = None,
) -> list[int]:
    target_count = target_record_count(len(records), event_fraction)
    if target_count == 0:
        return []
    identifiers = [fastq_record_id(record) for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("hash-ranked FASTQ record identifiers are not unique")
    allowed = set(range(len(records))) if eligible is None else set(eligible)
    if len(allowed) < target_count:
        raise ValueError(
            f"only {len(allowed)} eligible records for {target_count} requested mutations"
        )
    ranked = sorted(
        allowed,
        key=lambda index: (
            hashlib.sha256(f"{seed}\0{identifiers[index]}".encode()).digest(),
            identifiers[index],
        ),
    )
    return sorted(ranked[:target_count])


def mutate_fastq_records(
    records: list[FastqRecord],
    indices: list[int],
    mutation: Callable[[bytes], bytes],
) -> list[FastqRecord]:
    selected = set(indices)
    result = []
    for index, record in enumerate(records):
        if index not in selected:
            result.append(record)
            continue
        sequence = record[1].rstrip(b"\r\n")
        quality = record[3].rstrip(b"\r\n")
        changed = mutation(sequence)
        if changed == sequence:
            raise ValueError(
                f"FASTQ mutation did not change record {fastq_record_id(record)}"
            )
        if len(changed) > len(quality):
            raise ValueError(
                "FASTQ mutation made sequence longer than its quality string"
            )
        newline = b"\r\n" if record[1].endswith(b"\r\n") else b"\n"
        result.append(
            (record[0], changed + newline, record[2], quality[: len(changed)] + newline)
        )
    return result


def write_fastq(path: Path, records: list[FastqRecord]) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    sampler.write_fastq_gzip(path, records)
    return runtime.file_identity(path)


def mutation_for_read_operator(
    operator_id: str,
    variant: str,
    target: dict[str, Any],
    *,
    onlist_entries: set[str] | None = None,
) -> tuple[Callable[[bytes], bytes], Callable[[bytes], bool]]:
    regions = target.get("regions", [])
    if operator_id == "D01":
        truncate_to = strict_int(
            target.get("details", {}).get("truncate_to_bases"), "D01 truncation length"
        )
        return lambda sequence: sequence[:truncate_to], lambda sequence: len(
            sequence
        ) > truncate_to
    if operator_id == "D02":
        region = only(regions, "D02 region")
        position = strict_int(
            target.get("details", {}).get("target_position"), "D02 target position"
        )

        def mutate(sequence: bytes) -> bytes:
            if position >= len(sequence):
                return sequence
            replacement = substitute_base(chr(sequence[position])).encode()
            return sequence[:position] + replacement + sequence[position + 1 :]

        return mutate, lambda sequence: position < len(sequence)
    if operator_id == "D03":
        region = only(regions, "D03 region")
        start, stop = region_span(region, "D03")
        entries = onlist_entries or set()
        replacement = first_offlist_sequence(stop - start, entries).encode()
        return (
            lambda sequence: sequence[:start] + replacement + sequence[stop:]
            if len(sequence) >= stop
            else sequence,
            lambda sequence: len(sequence) >= stop
            and sequence[start:stop] != replacement,
        )
    if operator_id == "D04":
        primer = concrete_dna(only(regions, "D04 primer").get("sequence"), "D04 primer")
        motif = (
            reverse_complement(primer)
            if variant == "reverse_complement_by_substitution"
            else primer
        )
        motif_bytes = motif.encode()
        if variant in {
            "introduce_by_substitution",
            "reverse_complement_by_substitution",
        }:

            def introduce(sequence: bytes) -> bytes:
                position = insertion_window(sequence, motif_bytes)
                if position is None:
                    return sequence
                return (
                    sequence[:position]
                    + motif_bytes
                    + sequence[position + len(motif_bytes) :]
                )

            return introduce, lambda sequence: insertion_window(
                sequence, motif_bytes
            ) is not None
        if variant == "remove_by_substitution":

            def remove(sequence: bytes) -> bytes:
                position = sequence.find(motif_bytes)
                if position < 0:
                    return sequence
                base = substitute_base(chr(sequence[position])).encode()
                return sequence[:position] + base + sequence[position + 1 :]

            return remove, lambda sequence: sequence.find(motif_bytes) >= 0
        raise ValueError(f"unsupported D04 variant: {variant}")
    raise ValueError(f"unsupported read mutation operator: {operator_id}")


def insertion_window(sequence: bytes, motif: bytes) -> int | None:
    if len(sequence) < len(motif):
        return None
    for position in range(len(sequence) - len(motif) + 1):
        if sequence[position : position + len(motif)] != motif:
            return position
    return None


def region_span(region: dict[str, Any], label: str) -> tuple[int, int]:
    start = strict_int(region.get("start"), f"{label} region start")
    stop = strict_int(region.get("stop"), f"{label} region stop")
    if start < 0 or stop <= start:
        raise ValueError(f"{label} region span is invalid")
    return start, stop


def read_onlist_entries(source: str | Path, sequence_length: int) -> set[str]:
    if sequence_length <= 0:
        raise ValueError("onlist sequence length must be positive")
    value = str(source)
    if value.startswith(("http://", "https://", "ftp://")):
        with urllib.request.urlopen(value) as response:
            data = response.read()
    else:
        data = Path(value).read_bytes()
    if data.startswith(b"\x1f\x8b") or value.lower().endswith(".gz"):
        data = gzip.decompress(data)
    text = data.decode("utf-8")
    pattern = re.compile(rf"(?<![ACGT])[ACGT]{{{sequence_length}}}(?![ACGT])")
    return set(pattern.findall(text.upper()))


def first_offlist_sequence(length: int, forbidden: set[str]) -> str:
    if length <= 0:
        raise ValueError("offlist sequence length must be positive")
    normalized = {value.upper() for value in forbidden}
    value = 0
    limit = 4**length
    alphabet = "ACGT"
    while value < limit:
        current = value
        bases = ["A"] * length
        for index in range(length - 1, -1, -1):
            bases[index] = alphabet[current % 4]
            current //= 4
        candidate = "".join(bases)
        if candidate not in normalized:
            return candidate
        value += 1
    raise ValueError(f"onlist contains every possible DNA sequence of length {length}")


def substitute_sequence(sequence: str) -> str:
    return "".join(substitute_base(base) for base in sequence)


def substitute_base(base: str) -> str:
    value = base.upper()
    replacement = {"A": "C", "C": "G", "G": "T", "T": "A"}.get(value)
    if replacement is None:
        raise ValueError(f"cannot substitute non-ACGT base: {base}")
    return replacement


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def concrete_dna(value: Any, label: str) -> str:
    sequence = str(value or "").upper()
    if not sequence or any(base not in "ACGT" for base in sequence):
        raise ValueError(f"{label} is not concrete DNA")
    return sequence


def parse_shift_variant(value: str) -> int:
    if not value.startswith("shift_"):
        raise ValueError(f"invalid boundary shift variant: {value}")
    return strict_int(value.removeprefix("shift_"), "boundary shift")


def only(values: list[Any], label: str) -> Any:
    if len(values) != 1:
        raise ValueError(f"{label} must contain exactly one value")
    return values[0]


def strict_int(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer")
    try:
        result = int(str(value).strip())
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer") from error
    if str(value).strip() not in {str(result), f"+{result}"}:
        raise ValueError(f"{label} must be an integer")
    return result


def safe_name(value: str) -> str:
    result = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return result or "resource"


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
