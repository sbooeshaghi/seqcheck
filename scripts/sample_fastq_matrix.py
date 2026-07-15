#!/usr/bin/env python3
"""Create a deterministic FASTQ sampling matrix in one complete stream."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
import tempfile
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

try:
    import sample_fastq as sampler
except ModuleNotFoundError:
    from scripts import sample_fastq as sampler


SCHEMA_VERSION = "0.1.0"
MATRIX_VERSION = "0.1.0"
Condition = tuple[str, int, int | None]
IndexedGroup = tuple[int, tuple[sampler.FastqRecord, ...]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create all prefix and seeded reservoir FASTQ samples during one "
            "complete stream through each source."
        )
    )
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="Local FASTQ path or HTTP(S) URL. Repeat for multiple reads.",
    )
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument(
        "--n-reads",
        action="append",
        required=True,
        type=int,
        help="Records per sample. Repeat for each planned sample size.",
    )
    parser.add_argument(
        "--seed",
        action="append",
        required=True,
        type=int,
        help="Reservoir seed. Repeat for independent samples.",
    )
    parser.add_argument(
        "--include-prefix",
        action="store_true",
        help="Also write one prefix sample at each requested size.",
    )
    parser.add_argument(
        "--synchronize-mates",
        action="store_true",
        help="Validate read names and select the same indices from every input.",
    )
    parser.add_argument("--configuration-accession", default="")
    parser.add_argument("--modality", default="")
    parser.add_argument("--fastq-accession", action="append", default=[])
    parser.add_argument("--read-id", action="append", default=[])
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = sample_fastq_matrix(
            inputs=args.input,
            output_root=args.output_root.resolve(),
            sample_sizes=args.n_reads,
            seeds=args.seed,
            include_prefix=args.include_prefix,
            synchronize_mates=args.synchronize_mates,
            configuration_accession=args.configuration_accession,
            modality=args.modality,
            fastq_accessions=args.fastq_accession,
            read_ids=args.read_id,
        )
    except (OSError, ValueError, sampler.SampleError) as error:
        print(f"sample_fastq_matrix: {error}", file=sys.stderr)
        return 1
    print(manifest["manifest_path"])
    return 0


def sample_fastq_matrix(
    *,
    inputs: list[str],
    output_root: Path,
    sample_sizes: list[int],
    seeds: list[int],
    include_prefix: bool,
    synchronize_mates: bool,
    configuration_accession: str = "",
    modality: str = "",
    fastq_accessions: list[str] | None = None,
    read_ids: list[str] | None = None,
) -> dict[str, Any]:
    validate_inputs(
        inputs,
        sample_sizes,
        seeds,
        synchronize_mates,
    )
    sample_sizes = sorted(sample_sizes)
    seeds = sorted(seeds)
    if output_root.exists():
        raise ValueError(f"refusing to overwrite existing output root: {output_root}")
    accessions = sampler.aligned_metadata(
        fastq_accessions or [], len(inputs), "accession"
    )
    seqspec_read_ids = sampler.aligned_metadata(read_ids or [], len(inputs), "read id")
    filenames = sampler.output_filenames(inputs, accessions)

    with ExitStack() as stack:
        opened = [stack.enter_context(sampler.open_source(source)) for source in inputs]
        readers = [item[0] for item in opened]
        source_metadata = [item[1] for item in opened]
        digests = [hashlib.sha256() for _ in inputs]
        iterators = [
            sampler.iter_fastq(reader, source, digest)
            for reader, source, digest in zip(readers, inputs, digests, strict=True)
        ]
        if synchronize_mates:
            selected, streamed = sample_synchronized_matrix(
                iterators,
                inputs,
                sample_sizes,
                seeds,
                include_prefix,
            )
        else:
            selected = []
            streamed = []
            for input_index, iterator in enumerate(iterators):
                input_selected, input_streamed = sample_single_matrix(
                    iterator,
                    sample_sizes,
                    seeds,
                    include_prefix,
                    input_index,
                )
                selected.append(input_selected)
                streamed.append(input_streamed)

    source_rows = [
        {
            "source": inputs[index],
            "configuration_accession": configuration_accession,
            "modality": modality,
            "fastq_accession": accessions[index],
            "read_id": seqspec_read_ids[index],
            "access": source_metadata[index],
            "records_streamed": streamed[index],
            "complete_stream_consumed": True,
            "source_uncompressed_sha256": digests[index].hexdigest(),
        }
        for index in range(len(inputs))
    ]
    matrix_tool = matrix_identity()
    sampler_tool = sampler.sampler_identity()
    conditions = condition_order(sample_sizes, seeds, include_prefix)
    output_root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_root.name}-", dir=output_root.parent
    ) as tmpdir:
        temporary_root = Path(tmpdir)
        condition_rows = []
        for condition in conditions:
            method, n_reads, seed = condition
            relative_root = condition_directory(condition)
            output_rows = []
            for input_index in range(len(inputs)):
                relative_path = relative_root / filenames[input_index]
                temporary_path = temporary_root / relative_path
                temporary_path.parent.mkdir(parents=True, exist_ok=True)
                records = selected[input_index][condition]
                sampler.write_fastq_gzip(temporary_path, records)
                output_rows.append(
                    {
                        "input_index": input_index,
                        "fastq_accession": accessions[input_index],
                        "read_id": seqspec_read_ids[input_index],
                        "records_selected": len(records),
                        "effective_seed": effective_seed(
                            method, seed, input_index, synchronize_mates
                        ),
                        "output_path": str(output_root / relative_path),
                        "output_sha256": sampler.file_sha256(temporary_path),
                    }
                )
            condition_identity = {
                "method": method,
                "requested_records_per_fastq": n_reads,
                "seed": seed,
                "inputs": [stable_source_row(row) for row in source_rows],
                "outputs": [stable_output_row(row) for row in output_rows],
                "tools": {
                    "matrix_script_sha256": matrix_tool["script_sha256"],
                    "sampler_script_sha256": sampler_tool["script_sha256"],
                },
            }
            condition_rows.append(
                {
                    **condition_identity,
                    "condition_id": sha256_json(condition_identity)[:16],
                    "outputs": output_rows,
                }
            )

        stable_manifest = {
            "matrix_schema_version": SCHEMA_VERSION,
            "configuration_accession": configuration_accession,
            "modality": modality,
            "sample_sizes": sample_sizes,
            "reservoir_seeds": seeds,
            "include_prefix": include_prefix,
            "synchronized_mates": synchronize_mates,
            "sources": [stable_source_row(row) for row in source_rows],
            "conditions": [
                {
                    **{key: value for key, value in row.items() if key != "outputs"},
                    "outputs": [stable_output_row(value) for value in row["outputs"]],
                }
                for row in condition_rows
            ],
            "tools": {
                "matrix": functional_matrix_identity(matrix_tool),
                "sampler": functional_sampler_identity(sampler_tool),
            },
        }
        matrix_id = sha256_json(stable_manifest)[:16]
        manifest_path = output_root / "sample_matrix_manifest.json"
        manifest = {
            **stable_manifest,
            "matrix_id": matrix_id,
            "created_at": utc_now(),
            "invocation": sys.argv,
            "source_records": source_rows,
            "conditions": condition_rows,
            "tooling": {
                "matrix": matrix_tool,
                "sampler": sampler_tool,
            },
            "manifest_path": str(manifest_path),
        }
        (temporary_root / manifest_path.name).write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary_root, output_root)
    return manifest


def validate_inputs(
    inputs: list[str],
    sample_sizes: list[int],
    seeds: list[int],
    synchronize_mates: bool,
) -> None:
    if not inputs:
        raise ValueError("at least one input is required")
    if not sample_sizes or any(value <= 0 for value in sample_sizes):
        raise ValueError("sample sizes must be positive")
    if len(sample_sizes) != len(set(sample_sizes)):
        raise ValueError("sample sizes must be unique")
    if not seeds or any(value < 0 for value in seeds):
        raise ValueError("reservoir seeds must be nonnegative")
    if len(seeds) != len(set(seeds)):
        raise ValueError("reservoir seeds must be unique")
    if synchronize_mates and len(inputs) < 2:
        raise ValueError("--synchronize-mates requires at least two inputs")


def sample_single_matrix(
    records: Iterator[sampler.FastqRecord],
    sample_sizes: list[int],
    seeds: list[int],
    include_prefix: bool,
    input_index: int,
) -> tuple[dict[Condition, list[sampler.FastqRecord]], int]:
    selected, randomizers = initialize_matrix(
        sample_sizes,
        seeds,
        include_prefix,
        seed_offset=input_index,
    )
    streamed = 0
    for record in records:
        update_matrix(selected, randomizers, streamed, record)
        streamed += 1
    return finalize_matrix(selected), streamed


def sample_synchronized_matrix(
    iterators: list[Iterator[sampler.FastqRecord]],
    sources: list[str],
    sample_sizes: list[int],
    seeds: list[int],
    include_prefix: bool,
) -> tuple[list[dict[Condition, list[sampler.FastqRecord]]], list[int]]:
    conditions = condition_order(sample_sizes, seeds, include_prefix)
    grouped: dict[Condition, list[IndexedGroup]] = {
        condition: [] for condition in conditions
    }
    group_randomizers = {
        condition: random.Random(condition[2])
        for condition in conditions
        if condition[0] == "reservoir"
    }
    streamed = 0
    sentinel = object()
    while True:
        values = [next(iterator, sentinel) for iterator in iterators]
        if all(value is sentinel for value in values):
            break
        if any(value is sentinel for value in values):
            raise sampler.SampleError(
                "synchronized FASTQ inputs contain different record counts"
            )
        records = tuple(value for value in values if value is not sentinel)
        names = [sampler.read_name(record) for record in records]
        if len(set(names)) != 1:
            rendered = ", ".join(
                f"{source}={name.decode('utf-8', errors='replace')}"
                for source, name in zip(sources, names, strict=True)
            )
            raise sampler.SampleError(
                f"synchronized FASTQ read names differ at record {streamed + 1}: "
                f"{rendered}"
            )
        update_matrix(grouped, group_randomizers, streamed, records)
        streamed += 1

    finalized_groups = finalize_matrix(grouped)
    per_input = []
    for input_index in range(len(iterators)):
        per_input.append(
            {
                condition: [group[input_index] for group in groups]
                for condition, groups in finalized_groups.items()
            }
        )
    return per_input, [streamed] * len(iterators)


def initialize_matrix(
    sample_sizes: list[int],
    seeds: list[int],
    include_prefix: bool,
    seed_offset: int,
) -> tuple[dict[Condition, list[Any]], dict[Condition, random.Random]]:
    selected: dict[Condition, list[Any]] = {}
    if include_prefix:
        for n_reads in sample_sizes:
            selected[("prefix", n_reads, None)] = []
    randomizers = {}
    for n_reads in sample_sizes:
        for seed in seeds:
            condition = ("reservoir", n_reads, seed)
            selected[condition] = []
            randomizers[condition] = random.Random(seed + seed_offset)
    return selected, randomizers


def update_matrix(
    selected: dict[Condition, list[Any]],
    randomizers: dict[Condition, random.Random],
    index: int,
    record: Any,
) -> None:
    for condition, values in selected.items():
        method, n_reads, _ = condition
        if method == "prefix":
            if len(values) < n_reads:
                values.append((index, record))
            continue
        sampler.update_reservoir(
            values,
            index,
            record,
            method,
            n_reads,
            randomizers[condition],
        )


def finalize_matrix(
    selected: dict[Condition, list[Any]],
) -> dict[Condition, list[Any]]:
    result = {}
    for condition, values in selected.items():
        values.sort(key=lambda item: item[0])
        result[condition] = [record for _, record in values]
    return result


def condition_order(
    sample_sizes: list[int], seeds: list[int], include_prefix: bool
) -> list[Condition]:
    conditions = []
    for n_reads in sample_sizes:
        if include_prefix:
            conditions.append(("prefix", n_reads, None))
        conditions.extend(("reservoir", n_reads, seed) for seed in seeds)
    return conditions


def condition_directory(condition: Condition) -> Path:
    method, n_reads, seed = condition
    root = Path(method) / f"n-{n_reads:09d}"
    return root if seed is None else root / f"seed-{seed:09d}"


def effective_seed(
    method: str,
    seed: int | None,
    input_index: int,
    synchronize_mates: bool,
) -> int | None:
    if method == "prefix":
        return None
    if seed is None:
        raise ValueError("reservoir conditions require a seed")
    return seed if synchronize_mates else seed + input_index


def stable_source_row(row: dict[str, Any]) -> dict[str, Any]:
    access = {
        key: value
        for key, value in row["access"].items()
        if key not in {"retrieved_at", "path", "modified_time_ns", "status"}
    }
    source = row["source"] if row["access"].get("kind") == "remote" else ""
    return {
        "source": source,
        "kind": row["access"].get("kind", ""),
        "configuration_accession": row["configuration_accession"],
        "modality": row["modality"],
        "fastq_accession": row["fastq_accession"],
        "read_id": row["read_id"],
        "access": access,
        "records_streamed": row["records_streamed"],
        "source_uncompressed_sha256": row["source_uncompressed_sha256"],
    }


def stable_output_row(row: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in row.items() if key != "output_path"}


def matrix_identity() -> dict[str, Any]:
    path = Path(__file__).resolve()
    root = path.parents[1]
    status = sampler.git_output(
        root,
        "status",
        "--porcelain",
        "--untracked-files=normal",
        "--",
        str(path.relative_to(root)),
    )
    return {
        "version": MATRIX_VERSION,
        "script_path": str(path),
        "script_sha256": sampler.file_sha256(path),
        "git_commit": sampler.git_output(root, "rev-parse", "HEAD"),
        "git_dirty": bool(status),
        "python": sys.version.split()[0],
    }


def functional_matrix_identity(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "version": value["version"],
        "script_sha256": value["script_sha256"],
        "python": value["python"],
    }


def functional_sampler_identity(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "version": value["version"],
        "script_sha256": value["script_sha256"],
    }


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
