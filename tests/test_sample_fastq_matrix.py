import gzip
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "sample_fastq_matrix.py"
SPEC = importlib.util.spec_from_file_location("sample_fastq_matrix", SCRIPT_PATH)
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


def sampled_names(path: str) -> list[str]:
    with gzip.open(path, "rt", encoding="ascii") as handle:
        lines = handle.read().splitlines()
    return [lines[index][1:].split()[0] for index in range(0, len(lines), 4)]


def condition_hashes(manifest: dict[str, object]) -> dict[str, list[str]]:
    return {
        row["condition_id"]: [output["output_sha256"] for output in row["outputs"]]
        for row in manifest["conditions"]
    }


class SampleFastqMatrixTests(unittest.TestCase):
    @staticmethod
    def sample(**kwargs: object) -> dict[str, object]:
        with (
            patch.object(
                MODULE,
                "matrix_identity",
                return_value={
                    "version": "0.1.0",
                    "script_path": "/source/sample_fastq_matrix.py",
                    "script_sha256": "matrix-script-hash",
                    "git_commit": "abc123",
                    "git_dirty": False,
                    "python": "3.12.0",
                },
            ),
            patch.object(
                MODULE.sampler,
                "sampler_identity",
                return_value={
                    "version": "0.1.0",
                    "git_commit": "abc123",
                    "git_dirty": False,
                    "working_tree_hash": "",
                    "script_sha256": "sampler-script-hash",
                },
            ),
        ):
            return MODULE.sample_fastq_matrix(**kwargs)

    def test_matrix_produces_every_condition_during_one_complete_stream(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "reads.fastq.gz"
            source_text = write_fastq(source, 30)
            manifest = self.sample(
                inputs=[str(source)],
                output_root=root / "matrix",
                sample_sizes=[6, 3],
                seeds=[18, 17],
                include_prefix=True,
                synchronize_mates=False,
                configuration_accession="IGVFFI1234TEST",
                modality="rna",
                fastq_accessions=["IGVFFI0001TEST"],
                read_ids=["Read1"],
            )
            condition_rows = manifest["conditions"]
            selected_counts = sorted(
                row["outputs"][0]["records_selected"] for row in condition_rows
            )
            output_paths = [
                Path(row["outputs"][0]["output_path"]) for row in condition_rows
            ]

            self.assertTrue(all(path.is_file() for path in output_paths))
            persisted = json.loads(
                Path(manifest["manifest_path"]).read_text(encoding="utf-8")
            )

        self.assertEqual(manifest["sample_sizes"], [3, 6])
        self.assertEqual(manifest["reservoir_seeds"], [17, 18])
        self.assertEqual(len(condition_rows), 6)
        self.assertEqual(selected_counts, [3, 3, 3, 6, 6, 6])
        self.assertEqual(manifest["source_records"][0]["records_streamed"], 30)
        self.assertTrue(manifest["source_records"][0]["complete_stream_consumed"])
        self.assertEqual(
            manifest["source_records"][0]["source_uncompressed_sha256"],
            hashlib.sha256(source_text.encode("ascii")).hexdigest(),
        )
        self.assertEqual(persisted["matrix_id"], manifest["matrix_id"])

    def test_matrix_identity_is_stable_across_output_roots_and_seed_sensitive(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "reads.fastq"
            write_fastq(source, 30)
            common = {
                "inputs": [str(source)],
                "sample_sizes": [6],
                "seeds": [17],
                "include_prefix": True,
                "synchronize_mates": False,
            }
            first = self.sample(output_root=root / "first", **common)
            second = self.sample(output_root=root / "second", **common)
            changed = self.sample(
                output_root=root / "changed",
                **{**common, "seeds": [18]},
            )

        self.assertEqual(first["matrix_id"], second["matrix_id"])
        self.assertEqual(condition_hashes(first), condition_hashes(second))
        self.assertNotEqual(first["matrix_id"], changed["matrix_id"])
        first_reservoir = next(
            row for row in first["conditions"] if row["method"] == "reservoir"
        )
        changed_reservoir = next(
            row for row in changed["conditions"] if row["method"] == "reservoir"
        )
        self.assertNotEqual(
            first_reservoir["outputs"][0]["output_sha256"],
            changed_reservoir["outputs"][0]["output_sha256"],
        )

    def test_synchronized_matrix_preserves_mate_names_for_every_condition(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            read1 = root / "R1.fastq"
            read2 = root / "R2.fastq"
            write_fastq(read1, 20, mate=1)
            write_fastq(read2, 20, mate=2)
            manifest = self.sample(
                inputs=[str(read1), str(read2)],
                output_root=root / "matrix",
                sample_sizes=[4, 8],
                seeds=[7, 11],
                include_prefix=True,
                synchronize_mates=True,
                fastq_accessions=["R1", "R2"],
                read_ids=["Read1", "Read2"],
            )
            condition_names = []
            for condition in manifest["conditions"]:
                first = sampled_names(condition["outputs"][0]["output_path"])
                second = sampled_names(condition["outputs"][1]["output_path"])
                condition_names.append(
                    (
                        [name.removesuffix("/1") for name in first],
                        [name.removesuffix("/2") for name in second],
                    )
                )

        self.assertEqual(len(condition_names), 6)
        for first, second in condition_names:
            self.assertEqual(first, second)
        for condition in manifest["conditions"]:
            expected_seed = condition["seed"]
            self.assertEqual(
                [output["effective_seed"] for output in condition["outputs"]],
                [expected_seed, expected_seed],
            )
        self.assertEqual(
            [row["records_streamed"] for row in manifest["source_records"]],
            [20, 20],
        )

    def test_unsynchronized_matrix_records_input_seed_offsets(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            read1 = root / "R1.fastq"
            read2 = root / "I1.fastq"
            write_fastq(read1, 10)
            write_fastq(read2, 10)
            manifest = self.sample(
                inputs=[str(read1), str(read2)],
                output_root=root / "matrix",
                sample_sizes=[4],
                seeds=[7],
                include_prefix=True,
                synchronize_mates=False,
            )

        prefix, reservoir = manifest["conditions"]
        self.assertEqual(prefix["method"], "prefix")
        self.assertEqual(
            [output["effective_seed"] for output in prefix["outputs"]],
            [None, None],
        )
        self.assertEqual(reservoir["method"], "reservoir")
        self.assertEqual(
            [output["effective_seed"] for output in reservoir["outputs"]],
            [7, 8],
        )

    def test_synchronized_matrix_rejects_mismatched_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            read1 = root / "R1.fastq"
            read2 = root / "R2.fastq"
            write_fastq(read1, 4, mate=1)
            read2.write_text(
                fastq_text(4, mate=2).replace("read-2", "wrong", 1),
                encoding="ascii",
            )
            output_root = root / "matrix"

            with self.assertRaisesRegex(MODULE.sampler.SampleError, "names differ"):
                self.sample(
                    inputs=[str(read1), str(read2)],
                    output_root=output_root,
                    sample_sizes=[2],
                    seeds=[1],
                    include_prefix=True,
                    synchronize_mates=True,
                )

            self.assertFalse(output_root.exists())

    def test_synchronized_matrix_rejects_mismatched_record_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            read1 = root / "R1.fastq"
            read2 = root / "R2.fastq"
            write_fastq(read1, 4, mate=1)
            write_fastq(read2, 3, mate=2)
            output_root = root / "matrix"

            with self.assertRaisesRegex(
                MODULE.sampler.SampleError, "different record counts"
            ):
                self.sample(
                    inputs=[str(read1), str(read2)],
                    output_root=output_root,
                    sample_sizes=[2],
                    seeds=[1],
                    include_prefix=True,
                    synchronize_mates=True,
                )

            self.assertFalse(output_root.exists())

    def test_matrix_rejects_duplicate_protocol_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "reads.fastq"
            write_fastq(source, 4)
            common = {
                "inputs": [str(source)],
                "output_root": root / "matrix",
                "include_prefix": True,
                "synchronize_mates": False,
            }
            with self.assertRaisesRegex(ValueError, "sample sizes must be unique"):
                self.sample(sample_sizes=[2, 2], seeds=[1], **common)
            with self.assertRaisesRegex(ValueError, "seeds must be unique"):
                self.sample(sample_sizes=[2], seeds=[1, 1], **common)

    def test_cli_writes_manifest_and_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "reads.fastq"
            write_fastq(source, 8)
            command = [
                sys.executable,
                str(SCRIPT_PATH),
                "--input",
                str(source),
                "--output-root",
                str(root / "matrix"),
                "--n-reads",
                "3",
                "--seed",
                "5",
                "--include-prefix",
                "--configuration-accession",
                "IGVFFI1234TEST",
            ]
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
            )
            repeated = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest_path = Path(completed.stdout.strip())
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

        self.assertEqual(len(manifest["conditions"]), 2)
        self.assertEqual(manifest["configuration_accession"], "IGVFFI1234TEST")
        self.assertNotEqual(repeated.returncode, 0)
        self.assertIn("refusing to overwrite", repeated.stderr)


if __name__ == "__main__":
    unittest.main()
