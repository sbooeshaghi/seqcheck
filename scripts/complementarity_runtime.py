"""Deterministic mutation helpers for the QC complementarity experiment."""

from __future__ import annotations

import copy
import hashlib
import math
from pathlib import Path
from typing import Any

try:
    import perturbation_runtime as perturb
except ModuleNotFoundError:
    from scripts import perturbation_runtime as perturb


RUNTIME_VERSION = "0.1.0"
FastqRecord = perturb.FastqRecord


def mutate_quality_records(
    records: list[FastqRecord],
    *,
    operator_id: str,
    severity: dict[str, Any],
    seed: int,
) -> tuple[list[FastqRecord], dict[str, Any]]:
    require_records(records)
    if operator_id == "Q01":
        selected = list(range(len(records)))
        result = replace_quality(records, selected, phred(severity, "phred"))
    elif operator_id == "Q02":
        score = phred(severity, "phred")
        fraction = unit_fraction(severity.get("cycle_fraction"), "cycle fraction")
        selected = list(range(len(records)))
        result = replace_quality_cycles(records, score=score, fraction=fraction)
    elif operator_id == "Q03":
        score = phred(severity, "phred")
        fraction = unit_fraction(severity.get("event_fraction"), "event fraction")
        selected = perturb.select_record_indices(
            records, event_fraction=fraction, seed=seed
        )
        require_selected(selected, operator_id)
        result = replace_quality(records, selected, score)
    elif operator_id == "Q04":
        low = phred(severity, "low_phred")
        high = phred(severity, "high_phred")
        if low >= high:
            raise ValueError("Q04 low quality must be below high quality")
        fraction = unit_fraction(severity.get("low_fraction"), "low fraction")
        selected = perturb.select_record_indices(
            records, event_fraction=fraction, seed=seed
        )
        require_selected(selected, operator_id)
        result = replace_quality(records, list(range(len(records))), high)
        result = replace_quality(result, selected, low)
    else:
        raise ValueError(f"unsupported quality operator: {operator_id}")
    validate_quality_only(records, result)
    changed_records = sum(
        before[3] != after[3] for before, after in zip(records, result)
    )
    if changed_records == 0:
        raise ValueError(f"{operator_id} did not change any quality strings")
    return result, {
        "operator_id": operator_id,
        "severity": severity,
        "seed": seed,
        "records_selected": len(selected),
        "records_changed": changed_records,
        "total_records": len(records),
    }


def mutate_composition_records(
    records: list[FastqRecord],
    *,
    operator_id: str,
    event_fraction: float,
    seed: int,
    target_span: tuple[int, int] | None = None,
    motif: str | None = None,
) -> tuple[list[FastqRecord], dict[str, Any]]:
    require_records(records)
    fraction = unit_fraction(event_fraction, "event fraction")
    if operator_id == "C02":
        result, selected, details = duplicate_sequences(
            records, event_fraction=fraction, seed=seed
        )
    else:
        start, stop = require_span(target_span)
        if operator_id == "C01":
            replacement = repeated_motif(motif, stop - start)
        elif operator_id == "C03":
            replacement = alternating_gc(stop - start)
        else:
            raise ValueError(f"unsupported composition operator: {operator_id}")
        eligible = {
            index
            for index, record in enumerate(records)
            if len(sequence(record)) >= stop
            and sequence(record)[start:stop] != replacement
        }
        selected = perturb.select_record_indices(
            records,
            event_fraction=fraction,
            seed=seed,
            eligible=eligible,
        )
        require_selected(selected, operator_id)
        result = perturb.mutate_fastq_records(
            records,
            selected,
            lambda value: value[:start] + replacement + value[stop:],
        )
        details = {"target_start": start, "target_stop": stop}
    validate_sequence_mutation(records, result)
    return result, {
        "operator_id": operator_id,
        "event_fraction": fraction,
        "seed": seed,
        "records_selected": len(selected),
        "total_records": len(records),
        **details,
    }


def mutate_invalid_spec(
    original: dict[str, Any], *, variant: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    spec = copy.deepcopy(original)
    if variant == "missing_required_root_region_id":
        roots = spec.get("library_spec")
        if not isinstance(roots, list) or not roots or not isinstance(roots[0], dict):
            raise ValueError("seqspec has no root library region")
        if "region_id" not in roots[0]:
            raise ValueError("root library region already lacks region_id")
        removed = roots[0].pop("region_id")
        return spec, {"removed_field": "region_id", "removed_value": removed}
    if variant == "region_min_exceeds_max":
        region = next(
            (
                value
                for value in perturb.iter_regions(spec)
                if isinstance(value.get("max_len"), int)
            ),
            None,
        )
        if region is None:
            raise ValueError("seqspec has no region with an integer max_len")
        before = {"min_len": region.get("min_len"), "max_len": region["max_len"]}
        region["min_len"] = region["max_len"] + 1
        return spec, {
            "region_id": region.get("region_id"),
            "before": before,
            "after": {"min_len": region["min_len"], "max_len": region["max_len"]},
        }
    raise ValueError(f"unsupported schema-invalid variant: {variant}")


def replace_quality(
    records: list[FastqRecord], indices: list[int], score: int
) -> list[FastqRecord]:
    selected = set(indices)
    result = []
    for index, record in enumerate(records):
        if index not in selected:
            result.append(record)
            continue
        value = bytes([score + 33]) * len(sequence(record))
        result.append((record[0], record[1], record[2], value + newline(record[3])))
    return result


def replace_quality_cycles(
    records: list[FastqRecord], *, score: int, fraction: float
) -> list[FastqRecord]:
    result = []
    for record in records:
        quality = record[3].rstrip(b"\r\n")
        count = max(1, math.ceil(len(quality) * fraction))
        changed = quality[:-count] + bytes([score + 33]) * count
        result.append((record[0], record[1], record[2], changed + newline(record[3])))
    return result


def duplicate_sequences(
    records: list[FastqRecord], *, event_fraction: float, seed: int
) -> tuple[list[FastqRecord], list[int], dict[str, Any]]:
    ranked = sorted(
        range(len(records)),
        key=lambda index: (
            hashlib.sha256(
                f"{seed}\0{perturb.fastq_record_id(records[index])}".encode()
            ).digest(),
            perturb.fastq_record_id(records[index]),
        ),
    )
    donor_index = ranked[0]
    donor_sequence = sequence(records[donor_index])
    donor_quality = records[donor_index][3].rstrip(b"\r\n")
    eligible = {
        index
        for index, record in enumerate(records)
        if index != donor_index
        and len(sequence(record)) == len(donor_sequence)
        and sequence(record) != donor_sequence
    }
    selected = perturb.select_record_indices(
        records,
        event_fraction=event_fraction,
        seed=seed,
        eligible=eligible,
    )
    require_selected(selected, "C02")
    result = list(records)
    for index in selected:
        record = records[index]
        result[index] = (
            record[0],
            donor_sequence + newline(record[1]),
            record[2],
            donor_quality + newline(record[3]),
        )
    return (
        result,
        selected,
        {"donor_record_id": perturb.fastq_record_id(records[donor_index])},
    )


def validate_quality_only(before: list[FastqRecord], after: list[FastqRecord]) -> None:
    if len(before) != len(after):
        raise ValueError("quality mutation changed the FASTQ record count")
    for original, changed in zip(before, after):
        if original[:3] != changed[:3]:
            raise ValueError("quality mutation changed a name, sequence, or separator")
        if len(sequence(changed)) != len(changed[3].rstrip(b"\r\n")):
            raise ValueError("quality mutation broke sequence/quality length parity")


def validate_sequence_mutation(
    before: list[FastqRecord], after: list[FastqRecord]
) -> None:
    if len(before) != len(after):
        raise ValueError("composition mutation changed the FASTQ record count")
    changed = 0
    for original, mutated in zip(before, after):
        if original[0] != mutated[0] or original[2] != mutated[2]:
            raise ValueError("composition mutation changed a name or separator")
        if len(sequence(original)) != len(sequence(mutated)):
            raise ValueError("composition mutation changed read length")
        if len(sequence(mutated)) != len(mutated[3].rstrip(b"\r\n")):
            raise ValueError(
                "composition mutation broke sequence/quality length parity"
            )
        changed += sequence(original) != sequence(mutated)
    if changed == 0:
        raise ValueError("composition mutation changed no sequences")


def require_records(records: list[FastqRecord]) -> None:
    if not records:
        raise ValueError("FASTQ must contain at least one record")
    if any(
        len(sequence(record)) != len(record[3].rstrip(b"\r\n")) for record in records
    ):
        raise ValueError("FASTQ sequence and quality lengths differ")


def phred(value: dict[str, Any], field: str) -> int:
    score = value.get(field)
    if isinstance(score, bool) or not isinstance(score, int) or not 0 <= score <= 41:
        raise ValueError(f"{field} must be an integer from 0 to 41")
    return score


def unit_fraction(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    result = float(value)
    if not 0 < result < 1:
        raise ValueError(f"{label} must be between zero and one")
    return result


def require_span(value: tuple[int, int] | None) -> tuple[int, int]:
    if value is None:
        raise ValueError("composition mutation requires a target span")
    start, stop = value
    if start < 0 or stop <= start:
        raise ValueError("composition target span is invalid")
    return start, stop


def repeated_motif(value: str | None, length: int) -> bytes:
    motif = str(value or "").upper()
    if not motif or any(base not in "ACGT" for base in motif):
        raise ValueError("adapter motif must contain only A, C, G, and T")
    return (motif * math.ceil(length / len(motif)))[:length].encode()


def alternating_gc(length: int) -> bytes:
    return ("GC" * math.ceil(length / 2))[:length].encode()


def sequence(record: FastqRecord) -> bytes:
    return record[1].rstrip(b"\r\n")


def newline(value: bytes) -> bytes:
    return b"\r\n" if value.endswith(b"\r\n") else b"\n"


def require_selected(indices: list[int], operator_id: str) -> None:
    if not indices:
        raise ValueError(f"{operator_id} selected no records at this sample size")


def write_fastq(path: Path, records: list[FastqRecord]) -> dict[str, Any]:
    return perturb.write_fastq(path, records)
