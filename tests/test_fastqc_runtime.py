import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "fastqc_runtime.py"
FIXTURE = ROOT / "tests" / "fixtures" / "fastqc"
SPEC = importlib.util.spec_from_file_location("fastqc_runtime", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class FastqcRuntimeTests(unittest.TestCase):
    def test_parser_preserves_modules_and_reconciles_summary(self) -> None:
        parsed = MODULE.parse_fastqc_output(FIXTURE)

        self.assertEqual(parsed["fastqc_version"], "0.12.1")
        self.assertEqual(parsed["filename"], "reads.fastq")
        self.assertEqual(
            [(value["name"], value["status"]) for value in parsed["modules"]],
            [
                ("Basic Statistics", "pass"),
                ("Per base sequence quality", "warn"),
                ("Adapter Content", "fail"),
            ],
        )
        self.assertEqual(parsed["modules"][1]["rows"][0], ["1", "20.0", "20.0"])

    def test_parser_rejects_summary_drift_and_truncated_modules(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            root.joinpath("fastqc_data.txt").write_text(
                FIXTURE.joinpath("fastqc_data.txt").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            root.joinpath("summary.txt").write_text(
                "PASS\tBasic Statistics\treads.fastq\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "does not reconcile"):
                MODULE.parse_fastqc_output(root)

            root.joinpath("fastqc_data.txt").write_text(
                ">>Basic Statistics\tpass\nFilename\treads.fastq\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "header is missing"):
                MODULE.parse_fastqc_data(root / "fastqc_data.txt")

    def test_status_worsening_and_output_names_are_explicit(self) -> None:
        self.assertTrue(MODULE.status_worsened("pass", "warn"))
        self.assertTrue(MODULE.status_worsened("warn", "fail"))
        self.assertFalse(MODULE.status_worsened("fail", "fail"))
        self.assertEqual(MODULE.output_stem(Path("reads.fastq.gz")), "reads_fastqc")
        self.assertEqual(MODULE.output_stem(Path("reads.bam")), "reads_fastqc")


if __name__ == "__main__":
    unittest.main()
