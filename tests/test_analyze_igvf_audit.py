import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "analyze_igvf_audit.py"
PROTOCOL_PATH = ROOT / "experiments" / "paper" / "protocol" / "igvf_audit.json"
FAMILY_RULES_PATH = ROOT / "docs" / "cohort_family_rules.json"
SPEC = importlib.util.spec_from_file_location("analyze_igvf_audit_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


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


def completed_run(
    configuration: str,
    modality: str,
    access_class: str,
    title: str,
    assay_term: str,
    sampled_records: int,
) -> dict[str, object]:
    return {
        "study_run_id": "audit-study",
        "configuration_accession": configuration,
        "modality": modality,
        "access_class": access_class,
        "sampled_record_count": sampled_records,
        "run_status": "completed",
        "max_attempt_count": 1,
        "preferred_assay_titles": title,
        "assay_term": assay_term,
    }


def failure(
    configuration: str,
    category: str = "transport",
) -> dict[str, object]:
    return {
        "study_run_id": "audit-study",
        "configuration_accession": configuration,
        "modality": "guide",
        "access_class": "public",
        "failure_category": category,
        "stage": "seqcheck",
        "reason": "command_timeout",
        "attempt_count": 3,
        "preferred_assay_titles": "Perturb-seq",
        "assay_term": "single-cell RNA sequencing assay",
    }


def diagnostic() -> dict[str, object]:
    return {
        "study_run_id": "audit-study",
        "configuration_accession": "C1",
        "modality": "rna",
        "access_class": "public",
        "check": "length",
        "assessment_type": "warning",
        "assessment_code": "length_out_of_range",
        "preferred_assay_titles": "10x multiome",
        "assay_term": "single-cell RNA sequencing assay",
    }


def metric(
    configuration: str,
    modality: str,
    access_class: str,
    title: str,
    assay_term: str,
    value: float,
    sequence_type: str,
    ontology_term: str,
) -> dict[str, object]:
    return {
        "study_run_id": "audit-study",
        "configuration_accession": configuration,
        "modality": modality,
        "access_class": access_class,
        "check": "length",
        "metric_side": "observed",
        "metric_name": "in_range_fraction",
        "data_kind": "scalar",
        "unit": "fraction",
        "value_json": json.dumps(value),
        "sequence_types": sequence_type,
        "ontology_terms": ontology_term,
        "preferred_assay_titles": title,
        "assay_term": assay_term,
    }


def write_audit_fixture(root: Path, *, include_runs: bool = True) -> Path:
    runs = []
    metrics = []
    diagnostics = []
    if include_runs:
        runs = [
            completed_run(
                "C1",
                "rna",
                "public",
                "10x multiome",
                "single-cell RNA sequencing assay",
                100,
            ),
            completed_run(
                "C2",
                "atac",
                "controlled",
                "snATAC-seq",
                "single-nucleus ATAC-seq",
                80,
            ),
        ]
        diagnostics = [diagnostic()]
        metrics = [
            metric(
                "C1",
                "rna",
                "public",
                "10x multiome",
                "single-cell RNA sequencing assay",
                0.8,
                "random",
                "RGN:partition:molecule",
            ),
            metric(
                "C1",
                "rna",
                "public",
                "10x multiome",
                "single-cell RNA sequencing assay",
                0.6,
                "random",
                "RGN:partition:molecule",
            ),
            metric(
                "C2",
                "atac",
                "controlled",
                "snATAC-seq",
                "single-nucleus ATAC-seq",
                0.4,
                "genomic",
                "RGN:measure:genomic",
            ),
        ]

    failures = [failure("C3")]
    selected_count = 3 if include_runs else 1
    protocol_id = MODULE.audit.audit_protocol_id(
        MODULE.audit.load_audit_protocol(PROTOCOL_PATH)
    )
    write_json(
        root / "manifests" / "study.json",
        {
            "audit_schema_version": "0.4.1",
            "run_id": "audit-study",
            "selected_configuration_count": selected_count,
            "frozen_inputs": {"audit_protocol": {"audit_protocol_id": protocol_id}},
        },
    )
    write_json(
        root / "validation" / "reconciliation.json",
        {
            "audit_schema_version": "0.4.1",
            "study_run_id": "audit-study",
            "valid": True,
            "counts": {
                "selected_configurations": selected_count,
                "runs": len(runs),
                "failures": len(failures),
                "diagnostics": len(diagnostics),
                "metrics": len(metrics),
                "sampled_records": sum(
                    int(row["sampled_record_count"]) for row in runs
                ),
            },
        },
    )
    write_csv(root / "runs.csv", runs)
    write_csv(root / "failures.csv", failures)
    write_csv(root / "diagnostics.csv", diagnostics)
    write_csv(root / "metrics.csv", metrics)
    repeatability_path = root / "repeatability.json"
    write_json(
        repeatability_path,
        {
            "valid": True,
            "primary_study_run_id": "audit-study",
            "match_fraction": 1.0,
            "inputs": {
                "protocol": {"sha256": MODULE.runtime.file_sha256(PROTOCOL_PATH)}
            },
        },
    )
    return repeatability_path


class AnalyzeIgvfAuditTests(unittest.TestCase):
    def analyze(self, root: Path, output_name: str = "analysis") -> Path:
        return MODULE.analyze_igvf_audit(
            audit_root=root / "audit",
            repeatability_path=root / "audit" / "repeatability.json",
            protocol_path=PROTOCOL_PATH,
            family_rules_path=FAMILY_RULES_PATH,
            output_root=root / output_name,
        )

    def test_analysis_summarizes_outcomes_metrics_and_assessments(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_audit_fixture(root / "audit")
            validation_path = self.analyze(root)
            validation = json.loads(validation_path.read_text())
            completion = MODULE.read_csv(
                root / "analysis" / "tables" / "completion_summary.csv"
            )
            metric_rows = MODULE.read_csv(
                root / "analysis" / "tables" / "metric_distributions.csv"
            )
            assessments = MODULE.read_csv(
                root / "analysis" / "tables" / "assessment_summary.csv"
            )
            failures = MODULE.read_csv(
                root / "analysis" / "tables" / "failure_summary.csv"
            )

        overall = next(
            row
            for row in completion
            if row["dimension"] == "overall" and row["group_id"] == "all"
        )
        molecule = next(
            row
            for row in metric_rows
            if row["dimension"] == "ontology_term"
            and row["group_id"] == "RGN:partition:molecule"
        )
        c1_families = {
            row["group_id"]
            for row in completion
            if row["dimension"] == "assay_family"
            and row["completed_runs"] == "1"
            and row["sampled_records"] == "100"
        }

        self.assertTrue(validation["valid"])
        self.assertFalse(validation["scientific_target_met"])
        self.assertEqual(overall["configuration_count"], "3")
        self.assertEqual(overall["completed_runs"], "2")
        self.assertEqual(overall["classified_failures"], "1")
        self.assertAlmostEqual(float(overall["completion_fraction"]), 2 / 3)
        self.assertEqual(c1_families, {"droplet_rna", "multiome"})
        self.assertEqual(molecule["raw_observations"], "2")
        self.assertEqual(molecule["independent_configurations"], "1")
        self.assertEqual(float(molecule["median"]), 0.7)
        self.assertEqual(assessments[0]["candidate_inconsistency"], "True")
        self.assertEqual(failures[0]["failure_category"], "transport")
        self.assertEqual(failures[0]["maximum_attempt_count"], "3")

    def test_analysis_id_is_independent_of_output_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_audit_fixture(root / "audit")
            first = json.loads(self.analyze(root, "first").read_text())
            second = json.loads(self.analyze(root, "second").read_text())

        self.assertEqual(first["analysis_id"], second["analysis_id"])
        self.assertEqual(first["inputs"], second["inputs"])

    def test_zero_completed_runs_is_analyzable_but_misses_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_audit_fixture(root / "audit", include_runs=False)
            validation = json.loads(self.analyze(root).read_text())
            completion = MODULE.read_csv(
                root / "analysis" / "tables" / "completion_summary.csv"
            )

        overall = next(row for row in completion if row["dimension"] == "overall")
        self.assertTrue(validation["valid"])
        self.assertFalse(validation["scientific_target_met"])
        self.assertEqual(overall["completed_runs"], "0")
        self.assertEqual(overall["completion_fraction"], "0.0")
        self.assertEqual(validation["counts"]["numeric_observed_scalar_metrics"], 0)

    def test_eligibility_exclusions_do_not_enter_completion_denominator(self) -> None:
        rules = MODULE.load_family_rules(FAMILY_RULES_PATH)
        rows = MODULE.build_completion_summary(
            runs=[
                completed_run(
                    "C1",
                    "rna",
                    "public",
                    "10x multiome",
                    "single-cell RNA sequencing assay",
                    100,
                )
            ],
            failures=[failure("C2", "eligibility"), failure("C3")],
            families=rules,
            minimum_completion=0.95,
            analysis_id="analysis",
        )

        overall = next(row for row in rows if row["dimension"] == "overall")
        self.assertEqual(overall["configuration_count"], 3)
        self.assertEqual(overall["eligibility_exclusions"], 1)
        self.assertEqual(overall["eligible_outcomes"], 2)
        self.assertEqual(overall["completion_fraction"], 0.5)

    def test_protocol_identity_must_reconcile_across_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            repeatability_path = write_audit_fixture(root / "audit")
            study_path = root / "audit" / "manifests" / "study.json"
            study = json.loads(study_path.read_text())
            study["frozen_inputs"]["audit_protocol"]["audit_protocol_id"] = "other"
            write_json(study_path, study)
            with self.assertRaisesRegex(ValueError, "another audit protocol"):
                self.analyze(root, "invalid-study-protocol")

            study["frozen_inputs"]["audit_protocol"]["audit_protocol_id"] = (
                MODULE.audit.audit_protocol_id(
                    MODULE.audit.load_audit_protocol(PROTOCOL_PATH)
                )
            )
            write_json(study_path, study)
            repeatability = json.loads(repeatability_path.read_text())
            repeatability["inputs"]["protocol"]["sha256"] = "0" * 64
            write_json(repeatability_path, repeatability)
            with self.assertRaisesRegex(ValueError, "another audit protocol"):
                self.analyze(root, "invalid-repeat-protocol")

    def test_stale_reconciliation_counts_make_analysis_invalid(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            write_audit_fixture(root / "audit")
            reconciliation_path = root / "audit" / "validation" / "reconciliation.json"
            reconciliation = json.loads(reconciliation_path.read_text())
            reconciliation["counts"]["metrics"] += 1
            write_json(reconciliation_path, reconciliation)
            with self.assertRaisesRegex(ValueError, "audit_table_counts_reconcile"):
                self.analyze(root)
            validation = json.loads(
                (
                    root / "analysis" / "validation" / "igvf_audit_analysis.json"
                ).read_text()
            )

        self.assertFalse(validation["valid"])
        self.assertFalse(validation["checks"]["audit_table_counts_reconcile"])

    def test_invalid_repeatability_and_mismatched_rows_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            repeatability_path = write_audit_fixture(root / "audit")
            repeatability = json.loads(repeatability_path.read_text())
            repeatability["valid"] = False
            write_json(repeatability_path, repeatability)
            with self.assertRaisesRegex(ValueError, "repeatability validation"):
                self.analyze(root, "invalid-repeat")

            repeatability["valid"] = True
            write_json(repeatability_path, repeatability)
            rows = MODULE.read_csv(root / "audit" / "runs.csv")
            rows[0]["study_run_id"] = "another-study"
            write_csv(root / "audit" / "runs.csv", rows)
            with self.assertRaisesRegex(ValueError, "another study"):
                self.analyze(root, "invalid-row")

    def test_cli_writes_valid_analysis(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            repeatability_path = write_audit_fixture(root / "audit")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--audit-root",
                    str(root / "audit"),
                    "--repeatability-validation",
                    str(repeatability_path),
                    "--audit-protocol",
                    str(PROTOCOL_PATH),
                    "--family-rules",
                    str(FAMILY_RULES_PATH),
                    "--output-root",
                    str(root / "analysis"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            validation_path = Path(completed.stdout.strip())
            validation = json.loads(validation_path.read_text())

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(validation["valid"])

    def test_cli_returns_nonzero_for_invalid_analysis(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            repeatability_path = write_audit_fixture(root / "audit")
            reconciliation_path = root / "audit" / "validation" / "reconciliation.json"
            reconciliation = json.loads(reconciliation_path.read_text())
            reconciliation["counts"]["metrics"] += 1
            write_json(reconciliation_path, reconciliation)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--audit-root",
                    str(root / "audit"),
                    "--repeatability-validation",
                    str(repeatability_path),
                    "--audit-protocol",
                    str(PROTOCOL_PATH),
                    "--family-rules",
                    str(FAMILY_RULES_PATH),
                    "--output-root",
                    str(root / "analysis"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            validation = json.loads(
                (
                    root / "analysis" / "validation" / "igvf_audit_analysis.json"
                ).read_text()
            )

        self.assertEqual(completed.returncode, 1)
        self.assertIn("audit_table_counts_reconcile", completed.stderr)
        self.assertFalse(validation["valid"])


if __name__ == "__main__":
    unittest.main()
