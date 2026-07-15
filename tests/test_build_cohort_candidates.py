import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import call, patch


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "build_cohort_candidates.py"
)
SPEC = importlib.util.spec_from_file_location("build_cohort_candidates", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def rules_payload() -> dict:
    return {
        "schema_version": "0.2.0",
        "proposal_count_per_family": 4,
        "candidate_count_per_family": 2,
        "minimum_final_count_per_family": 1,
        "selection_seed": 7,
        "families": [
            {
                "id": "rna",
                "label": "RNA",
                "preferred_assay_titles": ["RNA-seq", "SHARE-seq"],
                "assay_terms": ["RNA assay"],
                "expected_modalities": ["rna"],
                "note": "review",
            },
            {
                "id": "multiome",
                "label": "Multiome",
                "preferred_assay_titles": ["SHARE-seq"],
                "assay_terms": ["RNA assay"],
                "expected_modalities": ["rna", "atac"],
                "note": "review both modalities",
            },
        ],
    }


def configuration(
    accession: str,
    title: str,
    lab: str,
    fastq: str,
    *,
    status: str = "released",
) -> dict:
    return {
        "accession": accession,
        "href": f"/configuration-files/{accession}/@@download/{accession}.yaml.gz",
        "lab": {"title": lab},
        "submitted_by": {"title": "Submitter"},
        "file_set": {
            "accession": f"DS-{accession}",
            "assay_term": {"term_name": "RNA assay"},
        },
        "preferred_assay_titles": [title],
        "aliases": [f"test:{accession}"],
        "seqspec_of": [f"/sequence-files/{fastq}/"],
        "status": status,
        "upload_status": "validated",
    }


def sequence(accession: str, controlled: bool = False) -> dict:
    return {
        "accession": accession,
        "href": f"/sequence-files/{accession}/@@download/{accession}.fastq.gz",
        "controlled_access": controlled,
        "read_names": ["R1"],
    }


class CohortCandidateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.rules = MODULE.parse_family_rules(rules_payload())

    def test_ambiguous_exact_title_matches_are_preserved(self) -> None:
        config = MODULE.parse_configurations(
            {"@graph": [configuration("C1", "SHARE-seq", "Lab A", "F1")]}
        )[0]
        matches = MODULE.match_families(config, self.rules)
        self.assertEqual([rule.family_id for rule, _ in matches], ["rna", "multiome"])

        inventory = MODULE.build_inventory(
            [config],
            MODULE.parse_sequence_files({"@graph": [sequence("F1")]}),
            self.rules,
            "https://example.org/",
        )
        self.assertEqual(len(inventory), 2)
        self.assertTrue(all(row["ambiguous_family"] for row in inventory))
        self.assertEqual(
            {row["candidate_families"] for row in inventory}, {"rna;multiome"}
        )

    def test_metadata_eligibility_rejects_controlled_missing_and_unknown_access(
        self,
    ) -> None:
        configs = MODULE.parse_configurations(
            {
                "@graph": [
                    configuration("PUBLIC", "RNA-seq", "Lab A", "F1"),
                    configuration("CONTROLLED", "RNA-seq", "Lab A", "F2"),
                    configuration("MISSING", "RNA-seq", "Lab A", "F3"),
                    configuration("UNKNOWN", "RNA-seq", "Lab A", "F4"),
                ]
            }
        )
        files = MODULE.parse_sequence_files(
            {
                "@graph": [
                    sequence("F1"),
                    sequence("F2", controlled=True),
                    {"accession": "F4"},
                ]
            }
        )
        observed = {}
        for config in configs:
            status, reasons, _ = MODULE.metadata_eligibility(config, files)
            observed[config.accession] = (status, reasons)
        self.assertEqual(observed["PUBLIC"][0], "eligible")
        self.assertIn("controlled_fastq", observed["CONTROLLED"][1])
        self.assertTrue(observed["MISSING"][1][0].startswith("missing_fastq_metadata"))
        self.assertIn("unknown_fastq_access", observed["UNKNOWN"][1])

    def test_proposals_are_deterministic_and_round_robin_labs(self) -> None:
        configs = MODULE.parse_configurations(
            {
                "@graph": [
                    configuration("A1", "RNA-seq", "Lab A", "F1"),
                    configuration("A2", "RNA-seq", "Lab A", "F2"),
                    configuration("A3", "RNA-seq", "Lab A", "F3"),
                    configuration("B1", "RNA-seq", "Lab B", "F4"),
                ]
            }
        )
        files = MODULE.parse_sequence_files(
            {"@graph": [sequence(f"F{index}") for index in range(1, 5)]}
        )
        inventory = MODULE.build_inventory(
            configs, files, self.rules, "https://example.org/"
        )
        forward = MODULE.select_proposals(inventory, self.rules, 2, 19)
        reverse = MODULE.select_proposals(list(reversed(inventory)), self.rules, 2, 19)
        self.assertEqual(forward, reverse)
        rna = [row for row in forward if row["family_id"] == "rna"]
        self.assertEqual({row["lab"] for row in rna}, {"Lab A", "Lab B"})
        self.assertEqual([row["proposal_rank"] for row in rna], [1, 2])

    def test_selection_represents_fastq_counts_within_a_lab(self) -> None:
        rows = [
            {
                "family_id": "rna",
                "configuration_accession": f"C{index}",
                "lab": "Lab A",
                "fastq_count": 4,
            }
            for index in range(5)
        ]
        rows.append(
            {
                "family_id": "rna",
                "configuration_accession": "RARE",
                "lab": "Lab A",
                "fastq_count": 3,
            }
        )

        selected = MODULE.select_lab_balanced(rows, "rna", 2, 19)

        self.assertEqual({row["fastq_count"] for row in selected}, {3, 4})

    def test_modality_screen_fills_review_pool_after_hydration(self) -> None:
        configs = MODULE.parse_configurations(
            {
                "@graph": [
                    configuration("A1", "RNA-seq", "Lab A", "F1"),
                    configuration("A2", "RNA-seq", "Lab A", "F2"),
                    configuration("A3", "RNA-seq", "Lab A", "F3"),
                    configuration("B1", "RNA-seq", "Lab B", "F4"),
                ]
            }
        )
        files = MODULE.parse_sequence_files(
            {"@graph": [sequence(f"F{index}") for index in range(1, 5)]}
        )
        inventory = MODULE.build_inventory(
            configs, files, self.rules, "https://example.org/"
        )
        proposals = MODULE.select_proposals(inventory, self.rules, 4, 19)
        modalities = {"A1": ("atac",), "A2": ("rna",), "A3": ("rna",), "B1": ("rna",)}
        proposals = [
            MODULE.enrich_candidate(
                row,
                MODULE.Hydration(
                    status="normalized",
                    modalities=modalities[row["configuration_accession"]],
                ),
            )
            for row in proposals
        ]

        annotated, selected = MODULE.select_review_candidates(
            proposals, self.rules, 2, 19
        )
        rna_selected = [row for row in selected if row["family_id"] == "rna"]
        self.assertEqual(len(rna_selected), 2)
        self.assertEqual({row["lab"] for row in rna_selected}, {"Lab A", "Lab B"})
        self.assertTrue(all(row["modalities"] == "rna" for row in rna_selected))
        mismatch = next(
            row
            for row in annotated
            if row["family_id"] == "rna" and row["configuration_accession"] == "A1"
        )
        self.assertEqual(mismatch["review_selection_status"], "excluded")
        self.assertIn("modality_mismatch", mismatch["review_exclusion_reason"])

    def test_diagnostics_do_not_change_review_candidate_membership(self) -> None:
        base = [
            {
                "family_id": "rna",
                "configuration_accession": accession,
                "lab": lab,
                "hydration_status": "normalized",
                "modality_match_status": "matched",
                "modalities": "rna",
                "family_expected_modalities": "rna",
                "proposal_rank": rank,
            }
            for rank, (accession, lab) in enumerate(
                [("A1", "Lab A"), ("A2", "Lab A"), ("B1", "Lab B")],
                start=1,
            )
        ]
        favorable = [
            row
            | {
                "structural_check_status": "passed",
                "resource_check_status": "passed",
                "fastq_mapping_status": "matched",
            }
            for row in base
        ]
        unfavorable = [
            row
            | {
                "structural_check_status": "failed",
                "resource_check_status": "unavailable",
                "fastq_mapping_status": "portal_extra",
            }
            for row in base
        ]
        _, first = MODULE.select_review_candidates(favorable, self.rules, 2, 19)
        _, second = MODULE.select_review_candidates(unfavorable, self.rules, 2, 19)

        def identity(rows: list[dict]) -> list[tuple[str, str]]:
            return [(row["family_id"], row["configuration_accession"]) for row in rows]

        self.assertEqual(identity(first), identity(second))

    def test_validation_reports_an_underfilled_modality_pool(self) -> None:
        payload = rules_payload()
        rules = MODULE.parse_family_rules(
            {**payload, "families": [payload["families"][0]]}
        )
        configs = MODULE.parse_configurations(
            {
                "@graph": [
                    configuration("C1", "RNA-seq", "Lab A", "F1"),
                    configuration("C2", "RNA-seq", "Lab B", "F2"),
                ]
            }
        )
        inventory = MODULE.build_inventory(
            configs,
            MODULE.parse_sequence_files({"@graph": [sequence("F1"), sequence("F2")]}),
            rules,
            "https://example.org/",
        )
        proposals = MODULE.select_proposals(inventory, rules, 2, 7)
        proposals = [
            MODULE.enrich_candidate(
                row,
                MODULE.Hydration(
                    status="normalized",
                    modalities=(
                        "rna" if row["configuration_accession"] == "C1" else "atac",
                    ),
                ),
            )
            for row in proposals
        ]
        proposals, selected = MODULE.select_review_candidates(proposals, rules, 2, 7)

        validation = MODULE.build_validation(
            inventory,
            proposals,
            selected,
            rules,
            2,
            2,
            2,
            "selection",
            "created",
        )
        family = validation["families"][0]
        self.assertEqual(family["modality_qualified_proposal_count"], 1)
        self.assertEqual(family["selected_candidate_count"], 1)
        self.assertFalse(family["minimum_final_pool_met"])
        self.assertFalse(validation["all_candidate_minimums_met"])
        self.assertTrue(
            any(
                "modality-qualified proposals" in blocker
                for blocker in validation["freeze_blockers"]
            )
        )

    def test_normalized_structure_excludes_fastq_file_bindings(self) -> None:
        library = {"rna": [{"region_type": ["RGN:partition:cell"]}]}
        reads_a = [
            {
                "read_id": "R1",
                "files": [{"file_id": "A", "filetype": "fastq"}],
                "min_len": 10,
            }
        ]
        reads_b = [
            {
                "read_id": "R1",
                "files": [{"file_id": "B", "filetype": "fastq"}],
                "min_len": 10,
            }
        ]
        structure_a = MODULE.normalized_structure(library, reads_a)
        structure_b = MODULE.normalized_structure(library, reads_b)
        self.assertEqual(structure_a, structure_b)
        self.assertEqual(
            MODULE.extract_ontology_terms(structure_a), {"RGN:partition:cell"}
        )
        self.assertEqual(MODULE.extract_expected_fastq_accessions(reads_a), {"A"})
        self.assertEqual(
            MODULE.normalize_fastq_accession("IGVFFI123.fastq.gz"), "IGVFFI123"
        )
        url_reads = [
            {
                "files": [
                    {
                        "file_id": "descriptive_R1",
                        "filetype": "fastq",
                        "url": "https://example.org/sequence-files/IGVFFI999/@@download/read.fastq.gz",
                    }
                ]
            }
        ]
        self.assertEqual(
            MODULE.extract_expected_fastq_accessions(url_reads), {"IGVFFI999"}
        )

    def test_validation_never_marks_unreviewed_candidates_frozen(self) -> None:
        config = MODULE.parse_configurations(
            {"@graph": [configuration("C1", "RNA-seq", "Lab A", "F1")]}
        )[0]
        inventory = MODULE.build_inventory(
            [config],
            MODULE.parse_sequence_files({"@graph": [sequence("F1")]}),
            self.rules,
            "https://example.org/",
        )
        proposals = MODULE.select_proposals(inventory, self.rules, 1, 7)
        proposals = [
            MODULE.enrich_candidate(
                row,
                MODULE.Hydration(
                    status="normalized",
                    modalities=("rna",),
                    structural_check_status="passed",
                    resource_check_status="passed",
                    fastq_mapping_status="matched",
                ),
            )
            for row in proposals
        ]
        proposals, selected = MODULE.select_review_candidates(
            proposals, self.rules, 1, 7
        )
        validation = MODULE.build_validation(
            inventory,
            proposals,
            selected,
            self.rules,
            1,
            1,
            1,
            "selection",
            "created",
        )
        self.assertFalse(validation["frozen"])
        self.assertFalse(validation["ready_to_freeze"])
        self.assertEqual(validation["required_independent_reviewers"], 2)
        self.assertFalse(validation["selection_uses_seqcheck_results"])

    def test_check_status_separates_check_and_transport_failures(self) -> None:
        self.assertEqual(MODULE.classify_check_status(0, ""), "passed")
        self.assertEqual(
            MODULE.classify_check_status(1, "[error 1] invalid region"),
            "failed",
        )
        self.assertEqual(
            MODULE.classify_check_status(1, "failed to send HTTP request"),
            "unavailable",
        )
        with (
            patch.object(
                MODULE.subprocess,
                "run",
                side_effect=[
                    subprocess.TimeoutExpired(["seqspec", "check"], 5),
                    subprocess.CompletedProcess(
                        ["seqspec", "check"], 0, stdout="valid", stderr=""
                    ),
                ],
            ),
            patch.object(MODULE.time, "sleep"),
        ):
            self.assertEqual(
                MODULE.run_check_command(["seqspec", "check"], 5, 2, 0),
                ("passed", "valid", 2),
            )

    def test_resource_check_is_not_run_after_structural_failure(self) -> None:
        with patch.object(
            MODULE,
            "run_check_command",
            return_value=("failed", "invalid structure", 1),
        ) as run_check:
            result = MODULE.run_seqspec_checks(
                ["seqspec"], Path("spec.yaml"), 5, 120, 3, 0
            )

        self.assertEqual(
            result,
            ("failed", "invalid structure", "not_run", "", 0),
        )
        run_check.assert_called_once_with(
            ["seqspec", "check", "--skip", "external", "spec.yaml"],
            5,
            1,
            0,
        )

    def test_resource_check_runs_after_structural_pass(self) -> None:
        with patch.object(
            MODULE,
            "run_check_command",
            side_effect=[
                ("passed", "", 1),
                ("unavailable", "endpoint unavailable", 3),
            ],
        ) as run_check:
            result = MODULE.run_seqspec_checks(
                ["seqspec"], Path("spec.yaml"), 5, 120, 3, 2
            )

        self.assertEqual(
            result,
            ("passed", "", "unavailable", "endpoint unavailable", 3),
        )
        self.assertEqual(
            run_check.call_args_list,
            [
                call(
                    ["seqspec", "check", "--skip", "external", "spec.yaml"],
                    5,
                    1,
                    0,
                ),
                call(
                    ["seqspec", "check", "spec.yaml"],
                    120,
                    3,
                    2,
                ),
            ],
        )

    def test_attempt_counts_do_not_change_selection_identity(self) -> None:
        common = ({}, 1, 1, 7, "configuration", "sequence", {})
        first_row = {"download_attempts": 1, "resource_check_attempts": 1}
        retried_row = {"download_attempts": 3, "resource_check_attempts": 2}
        first = MODULE.selection_identity_payload([first_row], [first_row], *common)
        retried = MODULE.selection_identity_payload(
            [retried_row], [retried_row], *common
        )
        self.assertEqual(first, retried)

    def test_resource_outcomes_do_not_change_selection_identity(self) -> None:
        common = ({}, 1, 1, 7, "configuration", "sequence", {})
        passed_row = {"resource_check_status": "passed"}
        unavailable_row = {"resource_check_status": "unavailable"}
        passed = MODULE.selection_identity_payload([passed_row], [passed_row], *common)
        unavailable = MODULE.selection_identity_payload(
            [unavailable_row], [unavailable_row], *common
        )
        self.assertEqual(passed, unavailable)

    def test_structural_and_fastq_outcomes_do_not_change_selection_identity(
        self,
    ) -> None:
        common = ({}, 1, 1, 7, "configuration", "sequence", {})
        passed_row = {
            "structural_check_status": "passed",
            "fastq_mapping_status": "matched",
        }
        failed_row = {
            "structural_check_status": "failed",
            "fastq_mapping_status": "portal_extra",
        }
        passed = MODULE.selection_identity_payload([passed_row], [passed_row], *common)
        failed = MODULE.selection_identity_payload([failed_row], [failed_row], *common)
        self.assertEqual(passed, failed)

    def test_cli_writes_blank_review_sheet_and_preserves_it_on_rerun(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            configurations = root / "configurations.json"
            sequences = root / "sequences.json"
            rules = root / "rules.json"
            output = root / "output"
            configurations.write_text(
                json.dumps(
                    {
                        "@graph": [
                            configuration("C1", "RNA-seq", "Lab A", "F1"),
                            configuration("C2", "SHARE-seq", "Lab B", "F2"),
                        ]
                    }
                ),
                encoding="utf-8",
            )
            sequences.write_text(
                json.dumps({"@graph": [sequence("F1"), sequence("F2")]}),
                encoding="utf-8",
            )
            payload = rules_payload()
            for family in payload["families"]:
                family["expected_modalities"] = []
            rules.write_text(json.dumps(payload), encoding="utf-8")
            argv = [
                sys.executable,
                str(SCRIPT_PATH),
                "--configurations",
                str(configurations),
                "--sequence-files",
                str(sequences),
                "--rules",
                str(rules),
                "--output-root",
                str(output),
            ]
            first = subprocess.run(argv, capture_output=True, text=True, check=False)
            self.assertEqual(first.returncode, 0, first.stderr)
            review_path = output / "tables" / "cohort_reviews.csv"
            with review_path.open(newline="", encoding="utf-8") as handle:
                review_rows = list(csv.DictReader(handle))
            self.assertTrue(review_rows)
            self.assertTrue(all(not row["reviewer_1_decision"] for row in review_rows))
            self.assertTrue(all(not row["split"] for row in review_rows))
            with (output / "tables" / "cohort_proposals.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                proposal_rows = list(csv.DictReader(handle))
            self.assertEqual(
                len(review_rows),
                sum(
                    row["review_selection_status"] == "selected"
                    for row in proposal_rows
                ),
            )

            review_rows[0]["reviewer_1_decision"] = "include"
            MODULE.write_csv(review_path, review_rows, MODULE.REVIEW_FIELDS)
            second = subprocess.run(argv, capture_output=True, text=True, check=False)
            self.assertEqual(second.returncode, 0, second.stderr)
            with review_path.open(newline="", encoding="utf-8") as handle:
                preserved = list(csv.DictReader(handle))
            self.assertEqual(preserved[0]["reviewer_1_decision"], "include")
            manifest = json.loads(
                (output / "manifests" / "cohort_candidates.json").read_text()
            )
            self.assertEqual(manifest["outputs"]["review_file_status"], "preserved")


if __name__ == "__main__":
    unittest.main()
