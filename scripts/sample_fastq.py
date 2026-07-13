#!/usr/bin/env python3
"""Create deterministic FASTQ samples for seqcheck experiments."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
from contextlib import ExitStack, contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, BinaryIO, Iterator


SAMPLE_SCHEMA_VERSION = "0.1.0"
SAMPLER_VERSION = "0.1.0"
USER_AGENT = f"seqcheck-fastq-sampler/{SAMPLER_VERSION}"
FastqRecord = tuple[bytes, bytes, bytes, bytes]


class SampleError(RuntimeError):
    """A FASTQ input cannot produce a valid study sample."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create deterministic prefix or reservoir FASTQ samples."
    )
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="Local FASTQ path or HTTP(S) URL. Repeat for multiple reads.",
    )
    parser.add_argument("--output-root", required=True, help="Sample output directory.")
    parser.add_argument(
        "--method",
        choices=("prefix", "reservoir"),
        default="prefix",
        help="Sampling method.",
    )
    parser.add_argument(
        "--n-reads",
        type=int,
        required=True,
        help="Records per FASTQ (0 retains the complete stream).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Reservoir seed. Recorded as null for prefix sampling.",
    )
    parser.add_argument(
        "--synchronize-mates",
        action="store_true",
        help="Validate read names and select the same record indices from every input.",
    )
    parser.add_argument("--configuration-accession", default="")
    parser.add_argument("--modality", default="")
    parser.add_argument(
        "--fastq-accession",
        action="append",
        default=[],
        help="FASTQ accession aligned with each --input. Repeatable.",
    )
    parser.add_argument(
        "--read-id",
        action="append",
        default=[],
        help="Seqspec read id aligned with each --input. Repeatable.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest = sample_fastqs(
            inputs=args.input,
            output_root=Path(args.output_root).resolve(),
            method=args.method,
            n_reads=args.n_reads,
            seed=args.seed,
            synchronize_mates=args.synchronize_mates,
            configuration_accession=args.configuration_accession,
            modality=args.modality,
            fastq_accessions=args.fastq_accession,
            read_ids=args.read_id,
        )
    except (OSError, SampleError) as err:
        print(f"sample_fastq: {err}", file=sys.stderr)
        return 1

    print(manifest["manifest_path"])
    return 0


def sample_fastqs(
    *,
    inputs: list[str],
    output_root: Path,
    method: str,
    n_reads: int,
    seed: int,
    synchronize_mates: bool,
    configuration_accession: str = "",
    modality: str = "",
    fastq_accessions: list[str] | None = None,
    read_ids: list[str] | None = None,
) -> dict[str, Any]:
    if not inputs:
        raise SampleError("at least one input is required")
    if method not in {"prefix", "reservoir"}:
        raise SampleError(f"unsupported sampling method: {method}")
    if n_reads < 0:
        raise SampleError("--n-reads must be non-negative")
    if synchronize_mates and len(inputs) < 2:
        raise SampleError("--synchronize-mates requires at least two inputs")

    accessions = aligned_metadata(fastq_accessions or [], len(inputs), "accession")
    seqspec_read_ids = aligned_metadata(read_ids or [], len(inputs), "read id")
    output_root.mkdir(parents=True, exist_ok=True)

    with ExitStack() as stack:
        opened = [stack.enter_context(open_source(source)) for source in inputs]
        readers = [item[0] for item in opened]
        source_metadata = [item[1] for item in opened]
        digests = [hashlib.sha256() for _ in inputs]
        iterators = [
            iter_fastq(reader, source, digest)
            for reader, source, digest in zip(readers, inputs, digests, strict=True)
        ]

        if synchronize_mates:
            selected, records_streamed, complete = sample_synchronized(
                iterators, inputs, method, n_reads, seed
            )
        else:
            selected = []
            records_streamed = []
            complete = []
            for input_index, iterator in enumerate(iterators):
                records, streamed, exhausted = sample_records(
                    iterator, method, n_reads, seed + input_index
                )
                selected.append(records)
                records_streamed.append(streamed)
                complete.append(exhausted)

    effective_seed = seed if method == "reservoir" else None
    output_rows = []
    with tempfile.TemporaryDirectory(prefix=".sample-", dir=output_root) as tmpdir:
        temporary_root = Path(tmpdir)
        for index, records in enumerate(selected):
            filename = output_filename(index, inputs[index], accessions[index])
            temporary_path = temporary_root / filename
            write_fastq_gzip(temporary_path, records)
            final_path = output_root / filename
            os.replace(temporary_path, final_path)
            output_rows.append(
                {
                    "source": inputs[index],
                    "configuration_accession": configuration_accession,
                    "modality": modality,
                    "fastq_accession": accessions[index],
                    "read_id": seqspec_read_ids[index],
                    "access": source_metadata[index],
                    "records_streamed": records_streamed[index],
                    "records_selected": len(records),
                    "complete_stream_consumed": complete[index],
                    "streamed_uncompressed_sha256": digests[index].hexdigest(),
                    "source_uncompressed_sha256": (
                        digests[index].hexdigest() if complete[index] else ""
                    ),
                    "output_path": str(final_path),
                    "output_sha256": file_sha256(final_path),
                }
            )

    stable_manifest = {
        "sample_schema_version": SAMPLE_SCHEMA_VERSION,
        "configuration_accession": configuration_accession,
        "modality": modality,
        "sampling_method": method,
        "requested_records_per_fastq": n_reads,
        "sampling_seed": effective_seed,
        "synchronized_mates": synchronize_mates,
        "inputs": output_rows,
        "sampler": sampler_identity(),
    }
    identity_manifest = {
        **stable_manifest,
        "inputs": [
            identity_input_row(row)
            for row in output_rows
        ],
        "sampler": {
            "version": stable_manifest["sampler"]["version"],
            "script_sha256": stable_manifest["sampler"]["script_sha256"],
        },
    }
    sample_id = sha256_json(identity_manifest)[:16]
    manifest_path = output_root / "sample_manifest.json"
    manifest = {
        **stable_manifest,
        "sample_id": sample_id,
        "created_at": utc_now(),
        "invocation": sys.argv,
        "manifest_path": str(manifest_path),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def aligned_metadata(values: list[str], size: int, label: str) -> list[str]:
    if not values:
        return [""] * size
    if len(values) != size:
        raise SampleError(
            f"received {len(values)} {label} values for {size} FASTQ inputs"
        )
    return values


def identity_input_row(row: dict[str, Any]) -> dict[str, Any]:
    identity = {key: value for key, value in row.items() if key != "output_path"}
    identity["access"] = {
        key: value
        for key, value in row["access"].items()
        if key != "retrieved_at"
    }
    return identity


@contextmanager
def open_source(source: str) -> Iterator[tuple[BinaryIO, dict[str, Any]]]:
    parsed = urllib.parse.urlparse(source)
    is_remote = parsed.scheme in {"http", "https"}
    if is_remote:
        request = urllib.request.Request(source, headers={"User-Agent": USER_AGENT})
        response = urllib.request.urlopen(request)
        raw: BinaryIO = response
        metadata = {
            "kind": "remote",
            "retrieved_at": utc_now(),
            "status": int(getattr(response, "status", 200)),
            "etag": response.headers.get("ETag", ""),
            "last_modified": response.headers.get("Last-Modified", ""),
            "content_length": response.headers.get("Content-Length", ""),
        }
        path = parsed.path
    else:
        path_obj = Path(source).resolve()
        raw = path_obj.open("rb")
        stat = path_obj.stat()
        metadata = {
            "kind": "local",
            "retrieved_at": utc_now(),
            "status": "available",
            "path": str(path_obj),
            "size_bytes": stat.st_size,
            "modified_time_ns": stat.st_mtime_ns,
        }
        path = str(path_obj)

    reader: BinaryIO
    if path.lower().endswith(".gz"):
        reader = gzip.GzipFile(fileobj=raw, mode="rb")
    else:
        reader = raw

    try:
        yield reader, metadata
    finally:
        reader.close()
        if reader is not raw:
            raw.close()


def iter_fastq(
    reader: BinaryIO, source: str, digest: Any
) -> Iterator[FastqRecord]:
    record_index = 0
    while True:
        header = reader.readline()
        if not header:
            return
        sequence = reader.readline()
        separator = reader.readline()
        quality = reader.readline()
        record_index += 1

        if not sequence or not separator or not quality:
            raise SampleError(f"{source}: incomplete FASTQ record {record_index}")
        if not header.startswith(b"@"):
            raise SampleError(f"{source}: record {record_index} header does not start with @")
        if not separator.startswith(b"+"):
            raise SampleError(f"{source}: record {record_index} separator does not start with +")
        if len(strip_newline(sequence)) != len(strip_newline(quality)):
            raise SampleError(
                f"{source}: record {record_index} sequence and quality lengths differ"
            )

        record = (header, sequence, separator, quality)
        for line in record:
            digest.update(line)
        yield record


def strip_newline(value: bytes) -> bytes:
    return value.rstrip(b"\r\n")


def read_name(record: FastqRecord) -> bytes:
    name = strip_newline(record[0])[1:].split(maxsplit=1)[0]
    if name.endswith((b"/1", b"/2")):
        return name[:-2]
    return name


def sample_records(
    records: Iterator[FastqRecord], method: str, n_reads: int, seed: int
) -> tuple[list[FastqRecord], int, bool]:
    selected: list[tuple[int, FastqRecord]] = []
    rng = random.Random(seed)
    streamed = 0
    exhausted = False

    while True:
        try:
            record = next(records)
        except StopIteration:
            exhausted = True
            break

        index = streamed
        streamed += 1
        update_reservoir(selected, index, record, method, n_reads, rng)
        if method == "prefix" and n_reads > 0 and len(selected) == n_reads:
            break

    selected.sort(key=lambda item: item[0])
    return [record for _, record in selected], streamed, exhausted


def sample_synchronized(
    iterators: list[Iterator[FastqRecord]],
    sources: list[str],
    method: str,
    n_reads: int,
    seed: int,
) -> tuple[list[list[FastqRecord]], list[int], list[bool]]:
    selected: list[tuple[int, tuple[FastqRecord, ...]]] = []
    rng = random.Random(seed)
    streamed = 0
    exhausted = False
    sentinel = object()

    while True:
        group = [next(iterator, sentinel) for iterator in iterators]
        if all(record is sentinel for record in group):
            exhausted = True
            break
        if any(record is sentinel for record in group):
            raise SampleError("synchronized FASTQ inputs contain different record counts")

        records = tuple(record for record in group if record is not sentinel)
        names = [read_name(record) for record in records]
        if len(set(names)) != 1:
            rendered = ", ".join(
                f"{source}={name.decode('utf-8', errors='replace')}"
                for source, name in zip(sources, names, strict=True)
            )
            raise SampleError(
                f"synchronized FASTQ read names differ at record {streamed + 1}: {rendered}"
            )

        index = streamed
        streamed += 1
        update_reservoir(selected, index, records, method, n_reads, rng)
        if method == "prefix" and n_reads > 0 and len(selected) == n_reads:
            break

    selected.sort(key=lambda item: item[0])
    per_input = [
        [group[input_index] for _, group in selected]
        for input_index in range(len(iterators))
    ]
    return per_input, [streamed] * len(iterators), [exhausted] * len(iterators)


def update_reservoir(
    selected: list[tuple[int, Any]],
    index: int,
    record: Any,
    method: str,
    n_reads: int,
    rng: random.Random,
) -> None:
    if method == "prefix" or n_reads == 0:
        selected.append((index, record))
        return
    if len(selected) < n_reads:
        selected.append((index, record))
        return

    replacement = rng.randrange(index + 1)
    if replacement < n_reads:
        selected[replacement] = (index, record)


def output_filename(index: int, source: str, accession: str) -> str:
    source_name = Path(urllib.parse.urlparse(source).path).name
    label = accession or source_name or f"input-{index + 1}"
    label = re.sub(r"(?:\.fastq|\.fq)?\.gz$", "", label, flags=re.IGNORECASE)
    label = re.sub(r"\.(?:fastq|fq)$", "", label, flags=re.IGNORECASE)
    safe_label = re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("._")
    return f"{index + 1:02d}_{safe_label or f'input-{index + 1}'}.sample.fastq.gz"


def write_fastq_gzip(path: Path, records: list[FastqRecord]) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as handle:
            for record in records:
                for line in record:
                    handle.write(line)


def sampler_identity() -> dict[str, Any]:
    root = Path(__file__).resolve().parents[1]
    commit = git_output(root, "rev-parse", "HEAD")
    status = git_output(
        root,
        "status",
        "--porcelain",
        "--untracked-files=normal",
        "--",
        "scripts/sample_fastq.py",
    )
    diff = subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "diff",
            "--binary",
            "HEAD",
            "--",
            "scripts/sample_fastq.py",
        ],
        check=False,
        capture_output=True,
    ).stdout
    return {
        "version": SAMPLER_VERSION,
        "git_commit": commit,
        "git_dirty": bool(status),
        "working_tree_hash": hashlib.sha256(diff).hexdigest() if diff else "",
        "script_sha256": file_sha256(Path(__file__).resolve()),
    }


def git_output(root: Path, *args: str) -> str:
    if not (root / ".git").exists():
        return ""
    result = subprocess.run(
        ["git", "-C", str(root), *args],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_json(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
