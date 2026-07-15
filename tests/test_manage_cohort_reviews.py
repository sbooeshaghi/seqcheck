import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "manage_cohort_reviews.py"
)
SPEC = importlib.util.spec_from_file_location("manage_cohort_reviews", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

FREEZE_SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "freeze_cohort.py"
FREEZE_SPEC = importlib.util.spec_from_file_location(
    "freeze_cohort_for_review_test", FREEZE_SCRIPT_PATH
)
FREEZE_MODULE = importlib.util.module_from_spec(FREEZE_SPEC)
assert FREEZE_SPEC.loader is not None
sys.modules[FREEZE_SPEC.name] = FREEZE_MODULE
FREEZE_SPEC.loader.exec_module(FREEZE_MODULE)


def candidate(accession: str, family: str = "rna") -> dict[str, str]:
    return {
        "selection_id": "selection",
        "family_id": family,
        "configuration_accession": accession,
        "configuration_url": f"https://example.org/{accession}",
        "lab": "Lab A" if accession == "C1" else "Lab B",
        "normalized_spec_path": f"/specs/{accession}.yaml",
        "normalized_spec_sha256": f"normalized-{accession}",
        "structural_check_status": "passed",
        "resource_check_status": "passed",
        "fastq_mapping_status": "matched",
    }


def write_csv(
    path: Path, rows: list[dict[str, str]], fields: list[str] | None = None
) -> None:
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        return [dict(row) for row in reader], list(reader.fieldnames)


def make_candidate_source(root: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    rows = [candidate("C1"), candidate("C2", "atac")]
    candidates_path = root / "candidates.csv"
    manifest_path = root / "candidate_manifest.json"
    write_csv(candidates_path, rows)
    manifest_path.write_text(
        json.dumps(
            {
                "cohort_candidate_schema_version": "0.2.0",
                "selection_id": "selection",
                "frozen": False,
                "outputs": {"cohort_candidates": str(candidates_path)},
            }
        ),
        encoding="utf-8",
    )
    return manifest_path, candidates_path, rows


def make_correction_registry(root: Path) -> Path:
    correction_dir = root / "corrections" / "C1"
    correction_dir.mkdir(parents=True)
    correction_path = correction_dir / "correction.json"
    correction_path.write_text(
        json.dumps(
            {
                "schema_version": "0.1.0",
                "status": "proposed_unapproved",
                "selection_id": "selection",
                "family_id": "rna",
                "configuration_accession": "C1",
                "rationale": "Add the omitted index read.",
                "approval": {
                    "reviewer_1": "",
                    "reviewer_1_decision": "",
                    "reviewer_2": "",
                    "reviewer_2_decision": "",
                },
            }
        ),
        encoding="utf-8",
    )
    registry_path = root / "corrections" / "correction_registry.json"
    registry_path.write_text(
        json.dumps(
            {
                "schema_version": "0.1.0",
                "selection_id": "selection",
                "corrections": [
                    {
                        "family_id": "rna",
                        "configuration_accession": "C1",
                        "manifest": "C1/correction.json",
                        "sha256": MODULE.file_sha256(correction_path),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return registry_path


def complete_sheet(
    path: Path,
    reviewer: str,
    *,
    decisions: tuple[str, ...] = ("include", "exclude"),
) -> None:
    rows, fields = read_csv(path)
    for row, decision in zip(rows, decisions, strict=True):
        row.update(
            {
                "reviewer": reviewer,
                "decision": decision,
                "rationale": "Protocol and read layout agree.",
                "protocol_url": "https://example.org/protocol",
                "date": "2026-07-14",
            }
        )
    write_csv(path, rows, fields)


class ManageCohortReviewsTests(unittest.TestCase):
    def test_combined_fields_match_freeze_contract(self) -> None:
        self.assertEqual(MODULE.COMBINED_REVIEW_FIELDS, FREEZE_MODULE.REVIEW_FIELDS)

    def test_prepare_creates_two_blinded_exact_packages(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _, candidates = make_candidate_source(root)
            output_root = root / "packages"

            manifests = MODULE.prepare_review_packages(
                candidate_manifest_path=manifest_path,
                candidate_path=None,
                output_root=output_root,
            )

            self.assertEqual(len(manifests), 2)
            self.assertNotEqual(manifests[0]["package_id"], manifests[1]["package_id"])
            for slot in (1, 2):
                sheet = output_root / f"reviewer_{slot}" / "cohort_review.csv"
                package = json.loads(
                    (
                        output_root
                        / f"reviewer_{slot}"
                        / "review_package.json"
                    ).read_text()
                )
                rows, fields = read_csv(sheet)
                self.assertEqual(package["review_slot"], slot)
                self.assertEqual(package["candidate_row_count"], len(candidates))
                self.assertEqual(
                    [
                        (row["family_id"], row["configuration_accession"])
                        for row in rows
                    ],
                    [
                        (row["family_id"], row["configuration_accession"])
                        for row in candidates
                    ],
                )
                self.assertEqual(
                    fields[-len(MODULE.REVIEW_INPUT_FIELDS) :],
                    list(MODULE.REVIEW_INPUT_FIELDS),
                )
                self.assertFalse(any("reviewer_1_decision" == field for field in fields))
                self.assertFalse(any("reviewer_2_decision" == field for field in fields))

            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                MODULE.prepare_review_packages(
                    candidate_manifest_path=manifest_path,
                    candidate_path=None,
                    output_root=output_root,
                )

    def test_merge_is_deterministic_and_uses_authoritative_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _, candidates = make_candidate_source(root)
            packages = root / "packages"
            MODULE.prepare_review_packages(
                candidate_manifest_path=manifest_path,
                candidate_path=None,
                output_root=packages,
            )
            sheet_1 = packages / "reviewer_1" / "cohort_review.csv"
            sheet_2 = packages / "reviewer_2" / "cohort_review.csv"
            package_1 = packages / "reviewer_1" / "review_package.json"
            package_2 = packages / "reviewer_2" / "review_package.json"
            complete_sheet(sheet_1, "Reviewer A")
            complete_sheet(sheet_2, "Reviewer B")

            outputs = []
            for name in ("merge-1", "merge-2"):
                output = root / name
                MODULE.merge_review_packages(
                    candidate_manifest_path=manifest_path,
                    candidate_path=None,
                    reviewer_1_package=package_1,
                    reviewer_1_sheet=sheet_1,
                    reviewer_2_package=package_2,
                    reviewer_2_sheet=sheet_2,
                    output_root=output,
                )
                outputs.append(output / "tables" / "cohort_reviews.csv")

            self.assertEqual(outputs[0].read_bytes(), outputs[1].read_bytes())
            merged, fields = read_csv(outputs[0])
            self.assertEqual(
                [{field: row[field] for field in candidates[0]} for row in merged],
                candidates,
            )
            self.assertEqual(merged[0]["reviewer_1"], "Reviewer A")
            self.assertEqual(merged[0]["reviewer_2"], "Reviewer B")
            self.assertEqual(merged[0]["reviewer_1_decision"], "include")
            self.assertTrue(all(row["final_family"] == "" for row in merged))
            self.assertEqual(
                fields[-len(MODULE.COMBINED_REVIEW_FIELDS) :],
                list(MODULE.COMBINED_REVIEW_FIELDS),
            )
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                MODULE.merge_review_packages(
                    candidate_manifest_path=manifest_path,
                    candidate_path=None,
                    reviewer_1_package=package_1,
                    reviewer_1_sheet=sheet_1,
                    reviewer_2_package=package_2,
                    reviewer_2_sheet=sheet_2,
                    output_root=root / "merge-1",
                )

    def test_merge_rejects_changed_candidate_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _, _ = make_candidate_source(root)
            packages = root / "packages"
            MODULE.prepare_review_packages(
                candidate_manifest_path=manifest_path,
                candidate_path=None,
                output_root=packages,
            )
            sheet_1 = packages / "reviewer_1" / "cohort_review.csv"
            sheet_2 = packages / "reviewer_2" / "cohort_review.csv"
            complete_sheet(sheet_1, "Reviewer A")
            complete_sheet(sheet_2, "Reviewer B")
            rows, fields = read_csv(sheet_1)
            rows[0]["lab"] = "Changed Lab"
            write_csv(sheet_1, rows, fields)

            with self.assertRaisesRegex(ValueError, "changed candidate evidence"):
                self.merge(root, manifest_path, packages)

    def test_correction_registry_is_read_only_review_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _, _ = make_candidate_source(root)
            registry_path = make_correction_registry(root)
            packages = root / "packages"
            MODULE.prepare_review_packages(
                candidate_manifest_path=manifest_path,
                candidate_path=None,
                output_root=packages,
                correction_registry_path=registry_path,
            )
            sheet_1 = packages / "reviewer_1" / "cohort_review.csv"
            sheet_2 = packages / "reviewer_2" / "cohort_review.csv"
            rows, fields = read_csv(sheet_1)
            self.assertTrue(rows[0]["proposed_correction_manifest"].endswith(
                "C1/correction.json"
            ))
            self.assertTrue(rows[0]["proposed_correction_sha256"])
            self.assertEqual(rows[1]["proposed_correction_manifest"], "")
            complete_sheet(sheet_1, "Reviewer A")
            complete_sheet(sheet_2, "Reviewer B")

            merged_root = root / "merged"
            MODULE.merge_review_packages(
                candidate_manifest_path=manifest_path,
                candidate_path=None,
                reviewer_1_package=packages
                / "reviewer_1"
                / "review_package.json",
                reviewer_1_sheet=sheet_1,
                reviewer_2_package=packages
                / "reviewer_2"
                / "review_package.json",
                reviewer_2_sheet=sheet_2,
                output_root=merged_root,
                correction_registry_path=registry_path,
            )
            merged, _ = read_csv(merged_root / "tables" / "cohort_reviews.csv")
            self.assertEqual(
                merged[0]["proposed_correction_manifest"],
                rows[0]["proposed_correction_manifest"],
            )

            tampered_root = root / "tampered"
            rows, fields = read_csv(sheet_1)
            rows[0]["proposed_correction_sha256"] = "changed"
            write_csv(sheet_1, rows, fields)
            with self.assertRaisesRegex(ValueError, "changed candidate evidence"):
                MODULE.merge_review_packages(
                    candidate_manifest_path=manifest_path,
                    candidate_path=None,
                    reviewer_1_package=packages
                    / "reviewer_1"
                    / "review_package.json",
                    reviewer_1_sheet=sheet_1,
                    reviewer_2_package=packages
                    / "reviewer_2"
                    / "review_package.json",
                    reviewer_2_sheet=sheet_2,
                    output_root=tampered_root,
                    correction_registry_path=registry_path,
                )

    def test_merge_rejects_a_changed_package_identity(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _, _ = make_candidate_source(root)
            packages = root / "packages"
            MODULE.prepare_review_packages(
                candidate_manifest_path=manifest_path,
                candidate_path=None,
                output_root=packages,
            )
            sheet_1 = packages / "reviewer_1" / "cohort_review.csv"
            sheet_2 = packages / "reviewer_2" / "cohort_review.csv"
            complete_sheet(sheet_1, "Reviewer A")
            complete_sheet(sheet_2, "Reviewer B")
            package_path = packages / "reviewer_1" / "review_package.json"
            package = json.loads(package_path.read_text())
            package["package_id"] = "changed-package"
            package_path.write_text(json.dumps(package), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "package id is invalid"):
                self.merge(root, manifest_path, packages)

    def test_merge_rejects_missing_extra_and_duplicate_rows(self) -> None:
        mutations = {
            "missing": lambda rows: rows[:1],
            "extra": lambda rows: [*rows, {**rows[0], "configuration_accession": "C3"}],
            "duplicate": lambda rows: [rows[0], rows[0]],
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                root = Path(tmpdir)
                manifest_path, _, _ = make_candidate_source(root)
                packages = root / "packages"
                MODULE.prepare_review_packages(
                    candidate_manifest_path=manifest_path,
                    candidate_path=None,
                    output_root=packages,
                )
                sheet_1 = packages / "reviewer_1" / "cohort_review.csv"
                sheet_2 = packages / "reviewer_2" / "cohort_review.csv"
                complete_sheet(sheet_1, "Reviewer A")
                complete_sheet(sheet_2, "Reviewer B")
                rows, fields = read_csv(sheet_1)
                write_csv(sheet_1, mutate(rows), fields)

                with self.assertRaisesRegex(ValueError, "row count|duplicate"):
                    self.merge(root, manifest_path, packages)

    def test_merge_requires_complete_distinct_reviewer_identities(self) -> None:
        for scenario in ("incomplete", "same-reviewer"):
            with self.subTest(scenario=scenario), tempfile.TemporaryDirectory() as tmpdir:
                root = Path(tmpdir)
                manifest_path, _, _ = make_candidate_source(root)
                packages = root / "packages"
                MODULE.prepare_review_packages(
                    candidate_manifest_path=manifest_path,
                    candidate_path=None,
                    output_root=packages,
                )
                sheet_1 = packages / "reviewer_1" / "cohort_review.csv"
                sheet_2 = packages / "reviewer_2" / "cohort_review.csv"
                complete_sheet(sheet_1, "Reviewer A")
                reviewer_2 = (
                    "Reviewer A" if scenario == "same-reviewer" else "Reviewer B"
                )
                complete_sheet(sheet_2, reviewer_2)
                if scenario == "incomplete":
                    rows, fields = read_csv(sheet_1)
                    rows[0]["rationale"] = ""
                    write_csv(sheet_1, rows, fields)

                pattern = "missing rationale|distinct people"
                with self.assertRaisesRegex(ValueError, pattern):
                    self.merge(root, manifest_path, packages)

    def test_merge_requires_extended_iso_review_dates(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _, _ = make_candidate_source(root)
            packages = root / "packages"
            MODULE.prepare_review_packages(
                candidate_manifest_path=manifest_path,
                candidate_path=None,
                output_root=packages,
            )
            sheet_1 = packages / "reviewer_1" / "cohort_review.csv"
            sheet_2 = packages / "reviewer_2" / "cohort_review.csv"
            complete_sheet(sheet_1, "Reviewer A")
            complete_sheet(sheet_2, "Reviewer B")
            rows, fields = read_csv(sheet_1)
            rows[0]["date"] = "20260714"
            write_csv(sheet_1, rows, fields)

            with self.assertRaisesRegex(ValueError, "date must use YYYY-MM-DD"):
                self.merge(root, manifest_path, packages)

    def test_cli_prepares_and_merges_review_packages(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _, _ = make_candidate_source(root)
            packages = root / "packages"
            prepare = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "prepare",
                    "--candidate-manifest",
                    str(manifest_path),
                    "--output-root",
                    str(packages),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(prepare.returncode, 0, prepare.stderr)
            complete_sheet(
                packages / "reviewer_1" / "cohort_review.csv", "Reviewer A"
            )
            complete_sheet(
                packages / "reviewer_2" / "cohort_review.csv", "Reviewer B"
            )
            merge = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "merge",
                    "--candidate-manifest",
                    str(manifest_path),
                    "--reviewer-1-package",
                    str(packages / "reviewer_1" / "review_package.json"),
                    "--reviewer-1-sheet",
                    str(packages / "reviewer_1" / "cohort_review.csv"),
                    "--reviewer-2-package",
                    str(packages / "reviewer_2" / "review_package.json"),
                    "--reviewer-2-sheet",
                    str(packages / "reviewer_2" / "cohort_review.csv"),
                    "--output-root",
                    str(root / "merged"),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(merge.returncode, 0, merge.stderr)
            self.assertTrue((root / "merged/tables/cohort_reviews.csv").is_file())
            self.assertTrue((root / "merged/manifests/cohort_reviews.json").is_file())

    def merge(self, root: Path, manifest_path: Path, packages: Path) -> None:
        MODULE.merge_review_packages(
            candidate_manifest_path=manifest_path,
            candidate_path=None,
            reviewer_1_package=packages / "reviewer_1" / "review_package.json",
            reviewer_1_sheet=packages / "reviewer_1" / "cohort_review.csv",
            reviewer_2_package=packages / "reviewer_2" / "review_package.json",
            reviewer_2_sheet=packages / "reviewer_2" / "cohort_review.csv",
            output_root=root / "merged",
        )


if __name__ == "__main__":
    unittest.main()
