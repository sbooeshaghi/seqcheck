import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "manage_audit_reviews.py"
ANALYSIS_SCRIPT_PATH = ROOT / "scripts" / "analyze_audit_reviews.py"
PROTOCOL_PATH = ROOT / "experiments" / "paper" / "protocol" / "audit_adjudication.json"
SPEC = importlib.util.spec_from_file_location("manage_audit_reviews_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
ANALYSIS_SPEC = importlib.util.spec_from_file_location(
    "analyze_audit_reviews_test", ANALYSIS_SCRIPT_PATH
)
ANALYSIS_MODULE = importlib.util.module_from_spec(ANALYSIS_SPEC)
assert ANALYSIS_SPEC.loader is not None
sys.modules[ANALYSIS_SPEC.name] = ANALYSIS_MODULE
ANALYSIS_SPEC.loader.exec_module(ANALYSIS_MODULE)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def family_rules() -> dict[str, object]:
    return {
        "schema_version": "0.2.0",
        "families": [
            {
                "id": "rna",
                "label": "RNA",
                "preferred_assay_titles": ["RNA assay"],
                "assay_terms": ["RNA sequencing"],
                "expected_modalities": ["rna"],
            },
            {
                "id": "atac",
                "label": "ATAC",
                "preferred_assay_titles": ["ATAC assay"],
                "assay_terms": ["ATAC sequencing"],
                "expected_modalities": ["atac"],
            },
        ],
    }


def adjudication_protocol() -> dict[str, object]:
    value = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    value["selection"]["minimum_candidate_findings"] = 3
    value["selection"]["pass_control_findings"] = 2
    return value


def run_row(
    configuration: str, modality: str, title: str, assay_term: str
) -> dict[str, object]:
    return {
        "study_run_id": "audit-study",
        "configuration_accession": configuration,
        "modality": modality,
        "access_class": "public",
        "sampled_record_count": 100,
        "run_status": "completed",
        "max_attempt_count": 1,
        "preferred_assay_titles": title,
        "assay_term": assay_term,
    }


def diagnostic_row(
    configuration: str,
    modality: str,
    title: str,
    assay_term: str,
    assessment_type: str,
    assessment_code: str,
) -> dict[str, object]:
    region = f"{configuration}-region"
    return {
        "study_run_id": "audit-study",
        "configuration_accession": configuration,
        "modality": modality,
        "access_class": "public",
        "lab": f"lab-{configuration}",
        "normalized_seqspec_version": "0.5.0",
        "check": "length",
        "assessment_type": assessment_type,
        "assessment_code": assessment_code,
        "assessment_description": f"Assessment for {configuration}",
        "files": f"{configuration}.fastq.gz",
        "reads": "Read1",
        "regions": region,
        "ontology_terms": "RGN:measure:transcript"
        if modality == "rna"
        else "RGN:measure:genomic",
        "sequence_types": "joined",
        "region_annotations_json": json.dumps([{"region_id": region}]),
        "preferred_assay_titles": title,
        "assay_term": assay_term,
    }


def metric_row(diagnostic: dict[str, object], value: float) -> dict[str, object]:
    return {
        "study_run_id": diagnostic["study_run_id"],
        "configuration_accession": diagnostic["configuration_accession"],
        "modality": diagnostic["modality"],
        "access_class": diagnostic["access_class"],
        "check": diagnostic["check"],
        "metric_side": "observed",
        "metric_name": "in_range_fraction",
        "data_kind": "scalar",
        "unit": "fraction",
        "value_json": json.dumps(value),
        "files": diagnostic["files"],
        "reads": diagnostic["reads"],
        "regions": diagnostic["regions"],
        "sequence_types": diagnostic["sequence_types"],
        "ontology_terms": diagnostic["ontology_terms"],
        "preferred_assay_titles": diagnostic["preferred_assay_titles"],
        "assay_term": diagnostic["assay_term"],
    }


def write_fixture(root: Path) -> dict[str, Path]:
    audit_root = root / "audit"
    family_path = root / "families.json"
    protocol_path = root / "protocol.json"
    policy_path = root / "detection_policy.json"
    analysis_path = root / "audit_analysis.json"
    write_json(family_path, family_rules())
    write_json(protocol_path, adjudication_protocol())

    definitions = [
        ("C1", "rna", "RNA assay", "RNA sequencing", "warning", "known_warning"),
        ("C2", "rna", "RNA assay", "RNA sequencing", "error", "hard_error"),
        ("C3", "atac", "ATAC assay", "ATAC sequencing", "interpretation", "info"),
        ("C4", "atac", "ATAC assay", "ATAC sequencing", "warning", "other_warning"),
        ("C5", "rna", "RNA assay", "RNA sequencing", "pass", "length_pass"),
        ("C6", "atac", "ATAC assay", "ATAC sequencing", "pass", "length_pass"),
    ]
    runs = [run_row(*value[:4]) for value in definitions]
    diagnostics = [diagnostic_row(*value) for value in definitions]
    metrics = [
        metric_row(value, 0.5 + index / 100) for index, value in enumerate(diagnostics)
    ]
    study = {
        "audit_schema_version": "0.4.1",
        "run_id": "audit-study",
        "selected_configuration_count": len(runs),
        "configurations": [
            {"accession": value[0], "href": f"/configuration-files/{value[0]}/"}
            for value in definitions
        ],
    }
    reconciliation = {"study_run_id": "audit-study", "valid": True}
    paths = {
        "study": audit_root / "manifests" / "study.json",
        "reconciliation": audit_root / "validation" / "reconciliation.json",
        "runs": audit_root / "runs.csv",
        "failures": audit_root / "failures.csv",
        "diagnostics": audit_root / "diagnostics.csv",
        "metrics": audit_root / "metrics.csv",
    }
    write_json(paths["study"], study)
    write_json(paths["reconciliation"], reconciliation)
    write_csv(paths["runs"], runs)
    write_csv(paths["failures"], [])
    write_csv(paths["diagnostics"], diagnostics)
    write_csv(paths["metrics"], metrics)

    policy_stable = {
        "schema_version": "0.1.0",
        "analysis_id": "calibration",
        "execution_id": "execution",
        "frozen": True,
        "calibration": {},
        "entries": [
            {
                "operator_id": "L01",
                "variant": "length",
                "ready": True,
                "assessment_codes": ["known_warning"],
            }
        ],
    }
    write_json(
        policy_path,
        {
            **policy_stable,
            "policy_id": MODULE.runtime.sha256_json(policy_stable)[:16],
            "created_at": "2026-07-14T00:00:00+00:00",
        },
    )
    analysis_inputs = {
        name: MODULE.runtime.file_sha256(path) for name, path in paths.items()
    }
    analysis_inputs.update(
        {
            "repeatability": "1" * 64,
            "protocol": "2" * 64,
            "family_rules": MODULE.runtime.file_sha256(family_path),
        }
    )
    analysis_stable = {
        "schema_version": "0.1.0",
        "study_run_id": "audit-study",
        "audit_protocol_id": "audit-protocol",
        "tool": {"version": "0.1.0", "sha256": "3" * 64, "python": "3.12"},
        "inputs": analysis_inputs,
    }
    write_json(
        analysis_path,
        {
            **analysis_stable,
            "analysis_id": MODULE.runtime.sha256_json(analysis_stable)[:16],
            "valid": True,
        },
    )
    return {
        "audit_root": audit_root,
        "audit_analysis_path": analysis_path,
        "detection_policy_path": policy_path,
        "family_rules_path": family_path,
        "protocol_path": protocol_path,
    }


def fill_review(path: Path, reviewer: str, *, disagree_first: bool = False) -> None:
    rows = read_csv(path)
    for index, row in enumerate(rows):
        row["reviewer"] = reviewer
        row["review_date"] = "2026-07-14"
        row["classification"] = (
            "confirmed_specification_problem"
            if disagree_first and index == 0
            else "intended_assay_design"
        )
        row["confidence"] = "high"
        row["rationale"] = f"Independent assessment by {reviewer}"
    write_csv(path, rows)


def fill_scientific_review(path: Path, reviewer: str, slot: int) -> None:
    rows = read_csv(path)
    for row in rows:
        configuration = row["configuration_accession"]
        if configuration == "C1":
            classification = "confirmed_specification_problem"
        elif configuration == "C2":
            classification = (
                "confirmed_read_or_resource_problem" if slot == 1 else "inconclusive"
            )
        elif configuration in {"C3", "C4"}:
            classification = "inconclusive"
        else:
            classification = "intended_assay_design"
        row["reviewer"] = reviewer
        row["review_date"] = "2026-07-14"
        row["classification"] = classification
        row["confidence"] = "high"
        row["rationale"] = f"Independent assessment by {reviewer}"
    write_csv(path, rows)


def prepare_merged_reviews(root: Path) -> tuple[dict[str, Path], Path, Path]:
    inputs = write_fixture(root)
    MODULE.prepare_review_packages(output_root=root / "packages", **inputs)
    first_sheet = root / "packages" / "reviewer_1" / "audit_review.csv"
    second_sheet = root / "packages" / "reviewer_2" / "audit_review.csv"
    fill_scientific_review(first_sheet, "Reviewer One", 1)
    fill_scientific_review(second_sheet, "Reviewer Two", 2)
    manifest = MODULE.merge_review_packages(
        reviewer_1_package=root / "packages" / "reviewer_1" / "review_package.json",
        reviewer_1_sheet=first_sheet,
        reviewer_2_package=root / "packages" / "reviewer_2" / "review_package.json",
        reviewer_2_sheet=second_sheet,
        output_root=root / "merged",
        **inputs,
    )
    merged_path = root / "merged" / "tables" / "audit_reviews.csv"
    completed_path = root / "adjudicated.csv"
    rows = read_csv(merged_path)
    for row in rows:
        row["submitter_contacted"] = "not_attempted"
        if row["needs_adjudication"] == "True":
            row["adjudicator"] = "Third Reviewer"
            row["adjudication_date"] = "2026-07-15"
            row["final_classification"] = "confirmed_read_or_resource_problem"
            row["adjudication_rationale"] = "Independent adjudication"
    write_csv(completed_path, rows)
    return inputs, Path(manifest["manifest_path"]), completed_path


class ManageAuditReviewsTests(unittest.TestCase):
    def test_prepare_blinds_source_and_selects_matched_controls(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = write_fixture(root)
            first = MODULE.prepare_review_packages(output_root=root / "first", **inputs)
            second = MODULE.prepare_review_packages(
                output_root=root / "second", **inputs
            )
            first_rows = read_csv(root / "first" / "reviewer_1" / "audit_review.csv")
            second_slot_rows = read_csv(
                root / "first" / "reviewer_2" / "audit_review.csv"
            )
            private = read_csv(root / "first" / "internal" / "review_key.csv")

        self.assertEqual(first["selection_id"], second["selection_id"])
        self.assertEqual(first["candidate_rows"], 3)
        self.assertEqual(first["pass_control_rows"], 2)
        self.assertEqual(first["high_confidence_candidate_rows"], 1)
        self.assertNotIn("source_kind", first_rows[0])
        self.assertNotIn("assessment_type", first_rows[0])
        self.assertNotIn("assessment_code", first_rows[0])
        self.assertNotIn("high_confidence", first_rows[0])
        self.assertNotIn("case_id", first_rows[0])
        self.assertTrue(
            set(row["review_item_id"] for row in first_rows).isdisjoint(
                row["review_item_id"] for row in second_slot_rows
            )
        )
        controls = [row for row in private if row["source_kind"] == "pass_control"]
        self.assertEqual(
            {row["matched_assay_family"] for row in controls}, {"rna", "atac"}
        )

    def test_merge_preserves_hidden_truth_and_marks_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = write_fixture(root)
            MODULE.prepare_review_packages(output_root=root / "packages", **inputs)
            first_sheet = root / "packages" / "reviewer_1" / "audit_review.csv"
            second_sheet = root / "packages" / "reviewer_2" / "audit_review.csv"
            fill_review(first_sheet, "Reviewer One")
            fill_review(second_sheet, "Reviewer Two", disagree_first=True)
            manifest = MODULE.merge_review_packages(
                reviewer_1_package=root
                / "packages"
                / "reviewer_1"
                / "review_package.json",
                reviewer_1_sheet=first_sheet,
                reviewer_2_package=root
                / "packages"
                / "reviewer_2"
                / "review_package.json",
                reviewer_2_sheet=second_sheet,
                output_root=root / "merged",
                **inputs,
            )
            rows = read_csv(root / "merged" / "tables" / "audit_reviews.csv")

        self.assertEqual(manifest["adjudication_required"], 1)
        self.assertFalse(manifest["adjudication_complete"])
        self.assertEqual(sum(row["needs_adjudication"] == "True" for row in rows), 1)
        self.assertEqual(sum(bool(row["final_classification"]) for row in rows), 4)
        self.assertEqual(
            {row["source_kind"] for row in rows}, {"candidate", "pass_control"}
        )

    def test_merge_rejects_changed_context_and_same_reviewer(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = write_fixture(root)
            MODULE.prepare_review_packages(output_root=root / "packages", **inputs)
            first_sheet = root / "packages" / "reviewer_1" / "audit_review.csv"
            second_sheet = root / "packages" / "reviewer_2" / "audit_review.csv"
            fill_review(first_sheet, "Same Reviewer")
            fill_review(second_sheet, "Same Reviewer")
            with self.assertRaisesRegex(ValueError, "distinct people"):
                MODULE.merge_review_packages(
                    reviewer_1_package=root
                    / "packages"
                    / "reviewer_1"
                    / "review_package.json",
                    reviewer_1_sheet=first_sheet,
                    reviewer_2_package=root
                    / "packages"
                    / "reviewer_2"
                    / "review_package.json",
                    reviewer_2_sheet=second_sheet,
                    output_root=root / "same-reviewer",
                    **inputs,
                )

            fill_review(second_sheet, "Reviewer Two")
            rows = read_csv(second_sheet)
            rows[0]["observed_metrics_json"] = "[]"
            write_csv(second_sheet, rows)
            with self.assertRaisesRegex(ValueError, "changed blinded evidence"):
                MODULE.merge_review_packages(
                    reviewer_1_package=root
                    / "packages"
                    / "reviewer_1"
                    / "review_package.json",
                    reviewer_1_sheet=first_sheet,
                    reviewer_2_package=root
                    / "packages"
                    / "reviewer_2"
                    / "review_package.json",
                    reviewer_2_sheet=second_sheet,
                    output_root=root / "tampered",
                    **inputs,
                )

    def test_prepare_rejects_unfrozen_policy(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = write_fixture(root)
            policy = json.loads(inputs["detection_policy_path"].read_text())
            policy["frozen"] = False
            write_json(inputs["detection_policy_path"], policy)
            with self.assertRaisesRegex(ValueError, "not content-addressed|not frozen"):
                MODULE.prepare_review_packages(output_root=root / "unfrozen", **inputs)

    def test_matched_control_selection_rejects_family_shortfall(self) -> None:
        candidates = [
            {"primary_assay_family": "rna"},
            {"primary_assay_family": "atac"},
        ]
        controls = [{"finding_key": "rna-control", "assay_family_ids": ["rna"]}]
        with self.assertRaisesRegex(ValueError, "assay family atac"):
            MODULE.select_matched_controls(
                controls=controls,
                candidates=candidates,
                target=2,
                seed=17,
            )

    def test_cli_prepare_writes_review_packages(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = write_fixture(root)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "prepare",
                    "--audit-root",
                    str(inputs["audit_root"]),
                    "--audit-analysis",
                    str(inputs["audit_analysis_path"]),
                    "--detection-policy",
                    str(inputs["detection_policy_path"]),
                    "--family-rules",
                    str(inputs["family_rules_path"]),
                    "--protocol",
                    str(inputs["protocol_path"]),
                    "--output-root",
                    str(root / "packages"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            manifest = json.loads(
                (
                    root / "packages" / "manifests" / "audit_review_packages.json"
                ).read_text()
            )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(manifest["review_rows"], 5)


class AnalyzeAuditReviewsTests(unittest.TestCase):
    def test_analysis_reports_conservative_and_evaluable_precision(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs, manifest_path, completed_path = prepare_merged_reviews(root)
            validation_path = ANALYSIS_MODULE.analyze_reviews(
                review_manifest_path=manifest_path,
                adjudicated_reviews_path=completed_path,
                protocol_path=inputs["protocol_path"],
                output_root=root / "analysis",
            )
            validation = json.loads(validation_path.read_text())
            precision = read_csv(
                root / "analysis" / "tables" / "confirmation_precision.csv"
            )
            endpoints = read_csv(root / "analysis" / "tables" / "review_endpoints.csv")

        overall = next(
            row
            for row in precision
            if row["dimension"] == "overall" and row["group_id"] == "all"
        )
        endpoint_values = {row["endpoint"]: row for row in endpoints}
        self.assertTrue(validation["valid"])
        self.assertEqual(validation["counts"]["candidate_cases"], 3)
        self.assertEqual(validation["counts"]["pass_control_cases"], 2)
        self.assertEqual(overall["confirmed_cases"], "2")
        self.assertEqual(overall["inconclusive_cases"], "1")
        self.assertAlmostEqual(float(overall["conservative_precision"]), 2 / 3)
        self.assertEqual(float(overall["evaluable_precision"]), 1.0)
        self.assertEqual(
            float(
                endpoint_values["high_confidence_conservative_confirmation_precision"][
                    "value"
                ]
            ),
            1.0,
        )
        self.assertEqual(
            float(endpoint_values["pass_control_consistency"]["value"]), 1.0
        )

    def test_analysis_rejects_changed_evidence_and_incomplete_adjudication(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs, manifest_path, completed_path = prepare_merged_reviews(root)
            rows = read_csv(completed_path)
            rows[0]["assessment_code"] = "tampered"
            write_csv(completed_path, rows)
            with self.assertRaisesRegex(ValueError, "changed review evidence"):
                ANALYSIS_MODULE.analyze_reviews(
                    review_manifest_path=manifest_path,
                    adjudicated_reviews_path=completed_path,
                    protocol_path=inputs["protocol_path"],
                    output_root=root / "tampered",
                )

            rows[0]["assessment_code"] = read_csv(
                root / "merged" / "tables" / "audit_reviews.csv"
            )[0]["assessment_code"]
            disagreement = next(
                row for row in rows if row["needs_adjudication"] == "True"
            )
            disagreement["adjudicator"] = ""
            write_csv(completed_path, rows)
            with self.assertRaisesRegex(ValueError, "adjudication is incomplete"):
                ANALYSIS_MODULE.analyze_reviews(
                    review_manifest_path=manifest_path,
                    adjudicated_reviews_path=completed_path,
                    protocol_path=inputs["protocol_path"],
                    output_root=root / "incomplete",
                )

    def test_analysis_cli_writes_valid_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs, manifest_path, completed_path = prepare_merged_reviews(root)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(ANALYSIS_SCRIPT_PATH),
                    "--review-manifest",
                    str(manifest_path),
                    "--adjudicated-reviews",
                    str(completed_path),
                    "--protocol",
                    str(inputs["protocol_path"]),
                    "--output-root",
                    str(root / "analysis"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            validation = json.loads(Path(completed.stdout.strip()).read_text())

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(validation["valid"])


if __name__ == "__main__":
    unittest.main()
