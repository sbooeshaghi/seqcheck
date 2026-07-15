import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from difflib import unified_diff
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "freeze_cohort.py"
SPEC = importlib.util.spec_from_file_location("freeze_cohort", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def candidate(accession: str, *, lab: str = "Lab A") -> dict[str, str]:
    return {
        "selection_id": "selection",
        "family_id": "rna",
        "configuration_accession": accession,
        "lab": lab,
        "hydration_status": "normalized",
        "normalized_seqspec_version": "0.5.0",
        "modality_match_status": "matched",
        "structural_check_status": "passed",
        "resource_check_status": "passed",
        "fastq_mapping_status": "matched",
        "deduplication_key": f"structure-{accession}",
    }


def review(row: dict[str, str], decision: str = "include") -> dict[str, str]:
    return {
        **row,
        "reviewer_1": "Reviewer A",
        "reviewer_1_decision": decision,
        "reviewer_1_rationale": "Protocol and read layout agree.",
        "reviewer_1_protocol_url": "https://example.org/protocol",
        "reviewer_1_date": "2026-07-14",
        "reviewer_2": "Reviewer B",
        "reviewer_2_decision": decision,
        "reviewer_2_rationale": "FASTQ roles and modality agree.",
        "reviewer_2_protocol_url": "https://example.org/protocol",
        "reviewer_2_date": "2026-07-14",
        "adjudication_decision": "",
        "adjudication_rationale": "",
        "final_family": row["family_id"] if decision == "include" else "",
        "split": "",
    }


class FreezeCohortTests(unittest.TestCase):
    def test_valid_reviews_freeze_with_deterministic_split(self) -> None:
        candidates = [candidate("C1"), candidate("C2", lab="Lab B")]
        reviews = [review(row) for row in candidates]
        review_fields = list(reviews[0])

        validation, included = MODULE.validate_reviews(
            selection_id="selection",
            candidate_rows=candidates,
            candidate_fields=list(candidates[0]),
            review_rows=reviews,
            review_fields=review_fields,
            family_ids=["rna"],
            target_per_family=2,
        )
        forward = MODULE.assign_splits(included, ["rna"], 1, 7, "selection")
        reverse = MODULE.assign_splits(
            list(reversed(included)), ["rna"], 1, 7, "selection"
        )

        self.assertTrue(validation["ready_to_freeze"])
        self.assertEqual(forward, reverse)
        self.assertEqual(
            {row["split"] for row in forward}, {"calibration", "evaluation"}
        )

    def test_invalid_reviews_and_baselines_remain_unfrozen(self) -> None:
        first = candidate("C1")
        second = candidate("C2")
        second["resource_check_status"] = "timeout"
        reviews = [review(first), review(second)]
        reviews[0]["reviewer_2"] = reviews[0]["reviewer_1"]
        reviews[1]["lab"] = "Changed Lab"

        validation, included = MODULE.validate_reviews(
            selection_id="selection",
            candidate_rows=[first, second],
            candidate_fields=list(first),
            review_rows=reviews,
            review_fields=list(reviews[0]),
            family_ids=["rna"],
            target_per_family=2,
        )
        messages = " ".join(
            problem for row in validation["row_errors"] for problem in row["errors"]
        )

        self.assertFalse(validation["ready_to_freeze"])
        self.assertEqual(included, [])
        self.assertIn("reviewers must be distinct", messages)
        self.assertIn("candidate fields changed", messages)
        self.assertIn("resource_check_status=passed", messages)

    def test_disagreement_requires_adjudication(self) -> None:
        row = candidate("C1")
        reviewed = review(row)
        reviewed["reviewer_2_decision"] = "exclude"

        decision, errors = MODULE.effective_review_decision(reviewed)
        self.assertEqual(decision, "")
        self.assertTrue(any("requires adjudication" in error for error in errors))

        reviewed["adjudication_decision"] = "include"
        reviewed["adjudication_rationale"] = "Protocol evidence resolves the conflict."
        self.assertEqual(
            MODULE.effective_review_decision(reviewed),
            ("include", []),
        )

    def test_adjudication_cannot_override_agreeing_reviewers(self) -> None:
        reviewed = review(candidate("C1"), decision="exclude")
        reviewed["adjudication_decision"] = "include"
        reviewed["adjudication_rationale"] = "Override both reviewers."

        decision, errors = MODULE.effective_review_decision(reviewed)

        self.assertEqual(decision, "")
        self.assertIn(
            "adjudication must remain blank when reviewer decisions agree", errors
        )

    def test_review_dates_require_extended_iso_format(self) -> None:
        reviewed = review(candidate("C1"))
        reviewed["reviewer_1_date"] = "20260714"

        _, errors = MODULE.effective_review_decision(reviewed)

        self.assertIn("reviewer_1_date must use YYYY-MM-DD", errors)

    def test_corrected_baseline_uses_effective_fields_and_requires_unanimity(
        self,
    ) -> None:
        row = candidate("C1") | {
            "fastq_mapping_status": "portal_extra_fastqs",
            "proposed_correction_manifest": "/correction.json",
            "proposed_correction_sha256": "correction-hash",
        }
        reviewed = review(row)
        overlay = MODULE.effective_candidate_fields(row) | {
            "correction_applied": "true",
            "correction_manifest": "/correction.json",
            "correction_manifest_sha256": "correction-hash",
            "effective_fastq_mapping_status": "matched",
            "effective_deduplication_key": "corrected-structure",
        }
        key = ("rna", "C1")

        validation, included = MODULE.validate_reviews(
            selection_id="selection",
            candidate_rows=[row],
            candidate_fields=list(candidate("C1")),
            review_rows=[reviewed],
            review_fields=list(reviewed),
            family_ids=["rna"],
            target_per_family=1,
            correction_overlays={key: overlay},
        )

        self.assertTrue(validation["ready_to_freeze"])
        self.assertEqual(included[0]["correction_applied"], "true")
        self.assertEqual(validation["included_correction_count"], 1)

        reviewed["reviewer_2_decision"] = "exclude"
        reviewed["adjudication_decision"] = "include"
        reviewed["adjudication_rationale"] = "Resolve disagreement."
        validation, _ = MODULE.validate_reviews(
            selection_id="selection",
            candidate_rows=[row],
            candidate_fields=list(candidate("C1")),
            review_rows=[reviewed],
            review_fields=list(reviewed),
            family_ids=["rna"],
            target_per_family=1,
            correction_overlays={key: overlay},
        )
        messages = " ".join(
            error
            for item in validation["row_errors"]
            for error in item["errors"]
        )
        self.assertFalse(validation["ready_to_freeze"])
        self.assertIn("requires unanimous include decisions", messages)

    def test_volatile_retry_fields_do_not_stale_a_review(self) -> None:
        original = candidate("C1") | {
            "download_attempts": "1",
            "resource_check_attempts": "1",
            "resource_check_message": "timeout after one second",
            "normalized_spec_path": "/first/spec.yaml",
        }
        rerun = original | {
            "download_attempts": "0",
            "resource_check_attempts": "2",
            "resource_check_message": "timeout after two seconds",
            "normalized_spec_path": "/second/spec.yaml",
        }

        self.assertEqual(
            MODULE.compare_candidate_fields(rerun, original, list(rerun)),
            [],
        )

    def test_cli_writes_immutable_frozen_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            candidates_path = root / "candidates.csv"
            reviews_path = root / "reviews.csv"
            rules_path = root / "rules.json"
            candidate_manifest_path = root / "candidate-manifest.json"
            output_root = root / "output"
            candidates = [candidate("C1"), candidate("C2", lab="Lab B")]
            reviews = [review(row) for row in candidates]
            write_csv(candidates_path, candidates, list(candidates[0]))
            write_csv(reviews_path, reviews, list(reviews[0]))
            rules_path.write_text(
                json.dumps({"families": [{"id": "rna"}]}), encoding="utf-8"
            )
            candidate_manifest_path.write_text(
                json.dumps(
                    {
                        "cohort_candidate_schema_version": "0.2.0",
                        "selection_id": "selection",
                        "frozen": False,
                        "rules": {"selection_seed": 7},
                        "outputs": {"cohort_candidates": str(candidates_path)},
                    }
                ),
                encoding="utf-8",
            )
            argv = [
                sys.executable,
                str(SCRIPT_PATH),
                "--candidate-manifest",
                str(candidate_manifest_path),
                "--reviews",
                str(reviews_path),
                "--family-rules",
                str(rules_path),
                "--output-root",
                str(output_root),
                "--target-per-family",
                "2",
                "--calibration-per-family",
                "1",
            ]

            result = subprocess.run(argv, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads(
                (output_root / "manifests" / "cohort_frozen.json").read_text()
            )
            self.assertTrue(manifest["frozen"])
            self.assertEqual(manifest["counts"]["calibration"], 1)
            self.assertEqual(manifest["counts"]["evaluation"], 1)

            rerun = subprocess.run(argv, capture_output=True, text=True, check=False)
            self.assertEqual(rerun.returncode, 1)
            self.assertIn("existing frozen cohort", rerun.stderr)

    def test_cli_freezes_a_unanimously_approved_correction(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
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
            correction_path = root / "correction.json"
            correction_path.write_text(
                json.dumps(
                    {
                        "schema_version": "0.1.0",
                        "status": "proposed_unapproved",
                        "selection_id": "selection",
                        "family_id": "rna",
                        "configuration_accession": "C1",
                        "rationale": "Add the omitted read.",
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
            registry_path = root / "correction_registry.json"
            registry_path.write_text(
                json.dumps(
                    {
                        "schema_version": "0.1.0",
                        "selection_id": "selection",
                        "corrections": [
                            {
                                "family_id": "rna",
                                "configuration_accession": "C1",
                                "manifest": correction_path.name,
                                "sha256": MODULE.file_sha256(correction_path),
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            seqspec = root / "seqspec"
            seqspec.write_text(
                """#!/usr/bin/env python3
import json
import sys

args = sys.argv[1:]
if args == ["--version"]:
    print("seqspec test")
elif args[0] == "version":
    print("seqspec file version: 0.5.0")
elif args[0] == "info":
    key = args[args.index("-k") + 1]
    values = {
        "modalities": ["rna"],
        "library_spec": {"rna": []},
        "sequence_spec": [{"read_id": "R1", "files": [{"file_id": "F1", "filetype": "fastq"}]}],
    }
    print(json.dumps(values[key]))
elif args[0] == "check":
    pass
else:
    raise SystemExit(2)
""",
                encoding="utf-8",
            )
            seqspec.chmod(0o755)

            first = candidate("C1") | {
                "family_expected_modalities": "rna",
                "normalized_spec_path": str(original),
                "normalized_spec_sha256": MODULE.file_sha256(original),
                "structure_sha256": "original-structure",
                "fastq_accessions": "F1",
                "fastq_set_sha256": "fastq-set-1",
                "expected_fastq_accessions": "",
                "fastq_mapping_status": "portal_extra_fastqs",
            }
            second = candidate("C2", lab="Lab B") | {
                "family_expected_modalities": "rna",
                "normalized_spec_path": "/specs/C2.yaml",
                "normalized_spec_sha256": "normalized-C2",
                "structure_sha256": "structure-C2",
                "fastq_accessions": "F2",
                "fastq_set_sha256": "fastq-set-2",
                "expected_fastq_accessions": "F2",
            }
            candidates = [first, second]
            correction_sha = MODULE.file_sha256(correction_path)
            reviews = [
                review(first)
                | {
                    "proposed_correction_manifest": str(correction_path.resolve()),
                    "proposed_correction_sha256": correction_sha,
                },
                review(second)
                | {
                    "proposed_correction_manifest": "",
                    "proposed_correction_sha256": "",
                },
            ]
            candidates_path = root / "candidates.csv"
            reviews_path = root / "reviews.csv"
            rules_path = root / "rules.json"
            candidate_manifest_path = root / "candidate-manifest.json"
            output_root = root / "output"
            write_csv(candidates_path, candidates, list(first))
            write_csv(reviews_path, reviews, list(reviews[0]))
            rules_path.write_text(
                json.dumps({"families": [{"id": "rna"}]}), encoding="utf-8"
            )
            candidate_manifest_path.write_text(
                json.dumps(
                    {
                        "cohort_candidate_schema_version": "0.2.0",
                        "selection_id": "selection",
                        "frozen": False,
                        "rules": {"selection_seed": 7},
                        "outputs": {"cohort_candidates": str(candidates_path)},
                    }
                ),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--candidate-manifest",
                    str(candidate_manifest_path),
                    "--reviews",
                    str(reviews_path),
                    "--family-rules",
                    str(rules_path),
                    "--correction-registry",
                    str(registry_path),
                    "--seqspec-bin",
                    str(seqspec),
                    "--output-root",
                    str(output_root),
                    "--target-per-family",
                    "2",
                    "--calibration-per-family",
                    "1",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            validation_path = output_root / "validation" / "cohort_freeze.json"
            failure_detail = (
                validation_path.read_text() if validation_path.exists() else result.stderr
            )
            self.assertEqual(result.returncode, 0, failure_detail)
            rows, _ = MODULE.read_csv(output_root / "tables" / "cohort.csv")
            by_accession = {row["configuration_accession"]: row for row in rows}
            self.assertEqual(by_accession["C1"]["correction_applied"], "true")
            self.assertEqual(
                by_accession["C1"]["effective_fastq_mapping_status"], "matched"
            )
            self.assertEqual(by_accession["C2"]["correction_applied"], "false")
            validation = json.loads(validation_path.read_text())
            self.assertEqual(validation["included_correction_count"], 1)
            frozen = json.loads(
                (output_root / "manifests" / "cohort_frozen.json").read_text()
            )
            self.assertEqual(frozen["counts"]["corrected"], 1)
            self.assertEqual(
                frozen["corrections"][0]["configuration_accession"], "C1"
            )


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
