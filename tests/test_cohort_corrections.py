import importlib.util
import json
import sys
import tempfile
import unittest
from difflib import unified_diff
from pathlib import Path
from unittest.mock import patch


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "cohort_corrections.py"
SPEC = importlib.util.spec_from_file_location("cohort_corrections_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def candidate(original: Path) -> dict[str, str]:
    return {
        "selection_id": "selection",
        "family_id": "rna",
        "configuration_accession": "C1",
        "family_expected_modalities": "rna",
        "normalized_spec_sha256": MODULE.file_sha256(original),
        "fastq_accessions": "F1",
        "fastq_set_sha256": "fastq-set",
    }


def write_correction(root: Path) -> tuple[Path, Path, Path, Path]:
    original = root / "original.yaml"
    corrected = root / "corrected.yaml"
    diff = root / "correction.diff"
    original.write_text("original\n", encoding="utf-8")
    corrected.write_text("corrected\n", encoding="utf-8")
    diff.write_text(
        "".join(
            unified_diff(
                original.read_text().splitlines(keepends=True),
                corrected.read_text().splitlines(keepends=True),
            )
        ),
        encoding="utf-8",
    )
    manifest_path = root / "correction.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": "0.1.0",
                "status": "proposed_unapproved",
                "selection_id": "selection",
                "family_id": "rna",
                "configuration_accession": "C1",
                "rationale": "Add omitted read metadata.",
                "original": {
                    "path": original.name,
                    "sha256": MODULE.file_sha256(original),
                },
                "corrected": {
                    "path": corrected.name,
                    "sha256": MODULE.file_sha256(corrected),
                },
                "diff": {
                    "path": diff.name,
                    "sha256": MODULE.file_sha256(diff),
                },
                "approval": {"reviewer_1": "", "reviewer_2": ""},
            }
        ),
        encoding="utf-8",
    )
    return original, corrected, diff, manifest_path


class CohortCorrectionTests(unittest.TestCase):
    def test_validate_correction_recomputes_effective_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            original, corrected, _, manifest_path = write_correction(root)

            def seqspec_output(argv: list[str], _: int) -> str:
                if "version" in argv:
                    return "seqspec file version: 0.5.0\n"
                key = argv[argv.index("-k") + 1]
                values = {
                    "modalities": ["rna"],
                    "library_spec": {"rna": []},
                    "sequence_spec": [
                        {
                            "read_id": "R1",
                            "files": [
                                {
                                    "file_id": "F1",
                                    "filename": "F1.fastq.gz",
                                    "filetype": "fastq",
                                }
                            ],
                        }
                    ],
                }
                return json.dumps(values[key])

            with patch.object(
                MODULE, "run_seqspec_command", side_effect=seqspec_output
            ), patch.object(MODULE, "run_seqspec_check") as run_check:
                effective = MODULE.validate_correction(
                    candidate=candidate(original),
                    registry_entry={
                        "path": str(manifest_path),
                        "sha256": MODULE.file_sha256(manifest_path),
                    },
                    seqspec_command=["seqspec"],
                    structural_timeout_seconds=5,
                    resource_timeout_seconds=5,
                    network_attempts=2,
                    retry_backoff_seconds=0,
                )

            self.assertEqual(effective["effective_fastq_mapping_status"], "matched")
            self.assertEqual(
                effective["effective_spec_sha256"], MODULE.file_sha256(corrected)
            )
            self.assertTrue(effective["effective_deduplication_key"])
            self.assertEqual(run_check.call_count, 2)

    def test_unified_diff_must_apply_at_recorded_location(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            original, corrected, diff, _ = write_correction(root)
            diff.write_text(
                "--- original\n+++ corrected\n@@ -2 +2 @@\n-original\n+corrected\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "original file|hunk"):
                MODULE.verify_unified_diff(original, corrected, diff)


if __name__ == "__main__":
    unittest.main()
