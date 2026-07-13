import gzip
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "sample_fastq.py"
SPEC = importlib.util.spec_from_file_location("sample_fastq", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def fastq_text(count: int, mate: int | None = None) -> str:
    records = []
    for index in range(count):
        suffix = f"/{mate}" if mate is not None else ""
        records.extend(
            [
                f"@read-{index}{suffix} metadata\n",
                f"ACG{index % 10}\n",
                "+\n",
                "!#I$\n",
            ]
        )
    return "".join(records)


def write_fastq(path: Path, count: int, mate: int | None = None) -> str:
    text = fastq_text(count, mate)
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="ascii") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="ascii")
    return text


def read_sample(path: str) -> str:
    with gzip.open(path, "rt", encoding="ascii") as handle:
        return handle.read()


def sampled_names(path: str) -> list[str]:
    lines = read_sample(path).splitlines()
    return [lines[index][1:].split()[0] for index in range(0, len(lines), 4)]


class SampleFastqTests(unittest.TestCase):
    @staticmethod
    def sample(**kwargs: object) -> dict[str, object]:
        with patch.object(
            MODULE,
            "sampler_identity",
            return_value={
                "version": "0.1.0",
                "git_commit": "abc123",
                "git_dirty": False,
                "working_tree_hash": "",
                "script_sha256": "script-hash",
            },
        ):
            return MODULE.sample_fastqs(**kwargs)

    def test_prefix_preserves_complete_records_and_quality(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "reads.fastq"
            original = write_fastq(source, 5)
            manifest = self.sample(
                inputs=[str(source)],
                output_root=root / "sample",
                method="prefix",
                n_reads=2,
                seed=99,
                synchronize_mates=False,
            )

            row = manifest["inputs"][0]
            observed = read_sample(row["output_path"])

        self.assertEqual(observed, "".join(original.splitlines(keepends=True)[:8]))
        self.assertEqual(row["records_streamed"], 2)
        self.assertEqual(row["records_selected"], 2)
        self.assertFalse(row["complete_stream_consumed"])
        self.assertEqual(manifest["sampling_seed"], None)
        self.assertEqual(row["source_uncompressed_sha256"], "")

    def test_reservoir_is_deterministic_and_consumes_complete_stream(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "reads.fastq.gz"
            write_fastq(source, 30)

            common = {
                "inputs": [str(source)],
                "method": "reservoir",
                "n_reads": 6,
                "seed": 17,
                "synchronize_mates": False,
            }
            first = self.sample(output_root=root / "first", **common)
            second = self.sample(output_root=root / "second", **common)
            changed = self.sample(
                output_root=root / "changed",
                **{**common, "seed": 18},
            )

            first_row = first["inputs"][0]
            second_row = second["inputs"][0]
            changed_row = changed["inputs"][0]

        self.assertEqual(first_row["output_sha256"], second_row["output_sha256"])
        self.assertNotEqual(first_row["output_sha256"], changed_row["output_sha256"])
        self.assertEqual(first["sample_id"], second["sample_id"])
        self.assertNotEqual(first["sample_id"], changed["sample_id"])
        self.assertEqual(first_row["records_streamed"], 30)
        self.assertEqual(first_row["records_selected"], 6)
        self.assertTrue(first_row["complete_stream_consumed"])
        self.assertEqual(
            first_row["streamed_uncompressed_sha256"],
            first_row["source_uncompressed_sha256"],
        )
        self.assertEqual(first["sampling_seed"], 17)

    def test_synchronized_reservoir_selects_matching_mate_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            read1 = root / "R1.fastq"
            read2 = root / "R2.fastq"
            write_fastq(read1, 12, mate=1)
            write_fastq(read2, 12, mate=2)

            manifest = self.sample(
                inputs=[str(read1), str(read2)],
                output_root=root / "sample",
                method="reservoir",
                n_reads=4,
                seed=7,
                synchronize_mates=True,
                fastq_accessions=["FASTQ1", "FASTQ2"],
                read_ids=["Read1", "Read2"],
            )
            first_names = sampled_names(manifest["inputs"][0]["output_path"])
            second_names = sampled_names(manifest["inputs"][1]["output_path"])

        self.assertEqual(
            [name.removesuffix("/1") for name in first_names],
            [name.removesuffix("/2") for name in second_names],
        )
        self.assertTrue(manifest["synchronized_mates"])
        self.assertEqual(manifest["inputs"][0]["read_id"], "Read1")
        self.assertEqual(manifest["inputs"][1]["fastq_accession"], "FASTQ2")

    def test_synchronized_sampling_rejects_mismatched_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            read1 = root / "R1.fastq"
            read2 = root / "R2.fastq"
            write_fastq(read1, 4, mate=1)
            read2.write_text(
                fastq_text(2, mate=2) + fastq_text(1, mate=2).replace("read-0", "wrong"),
                encoding="ascii",
            )

            with self.assertRaisesRegex(MODULE.SampleError, "read names differ"):
                self.sample(
                    inputs=[str(read1), str(read2)],
                    output_root=root / "sample",
                    method="reservoir",
                    n_reads=2,
                    seed=1,
                    synchronize_mates=True,
                )

    def test_sampling_rejects_truncated_quality(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "bad.fastq"
            source.write_text("@read\nACGT\n+\n!!!\n", encoding="ascii")

            with self.assertRaisesRegex(MODULE.SampleError, "lengths differ"):
                self.sample(
                    inputs=[str(source)],
                    output_root=root / "sample",
                    method="prefix",
                    n_reads=1,
                    seed=0,
                    synchronize_mates=False,
                )

    def test_cli_writes_manifest_with_study_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "reads.fastq"
            write_fastq(source, 8)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--input",
                    str(source),
                    "--output-root",
                    str(root / "sample"),
                    "--method",
                    "reservoir",
                    "--n-reads",
                    "3",
                    "--seed",
                    "5",
                    "--configuration-accession",
                    "IGVFFI1234TEST",
                    "--modality",
                    "rna",
                    "--fastq-accession",
                    "IGVFFI0001TEST",
                    "--read-id",
                    "Read1",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest_path = Path(completed.stdout.strip())
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

        self.assertEqual(manifest["configuration_accession"], "IGVFFI1234TEST")
        self.assertEqual(manifest["inputs"][0]["fastq_accession"], "IGVFFI0001TEST")
        self.assertEqual(manifest["inputs"][0]["records_selected"], 3)


if __name__ == "__main__":
    unittest.main()
