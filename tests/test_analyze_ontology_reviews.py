import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "analyze_ontology_reviews.py"
SPEC = importlib.util.spec_from_file_location(
    "analyze_ontology_reviews_test", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fake_yq(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if sys.argv[1:] == ["--version"]:
    print("yq test 4.0")
else:
    print(Path(sys.argv[-1]).read_text())
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def registry() -> dict[str, object]:
    roles = {
        "partition": {"scope": "direct_read_interval"},
        "technical": {"scope": "direct_read_interval"},
        "unknown": {"scope": "direct_read_interval"},
    }
    definitions = {
        "RGN:partition:cell": ("partition", "cell"),
        "RGN:partition:molecule": ("partition", "molecule"),
        "RGN:partition:sample": ("partition", "sample"),
        "RGN:technical:index7": ("technical", "index7"),
        "RGN:unknown:unclassified": ("unknown", "unclassified"),
    }
    terms = {
        term: {
            "label": term,
            "role": role,
            "target": target,
            "status": "active",
            "definition": term,
            "typical_sequence_types": [],
            "examples": [],
            "does_not_mean": [],
        }
        for term, (role, target) in definitions.items()
    }
    return {
        "ontology_id": "seqspec-region-ontology",
        "roles": roles,
        "terms": terms,
    }


def review_row(
    *,
    index: int,
    configuration: str,
    truth: str,
    reviewer_1: str,
    reviewer_2: str,
    weight: float,
) -> dict[str, object]:
    agreement = set(reviewer_1.split(";")) == set(reviewer_2.split(";"))
    row = {
        "survey_id": "survey-1",
        "region_key": f"region-{index}",
        "configuration_accession": configuration,
        "ontology_terms": truth,
        "sampling_weight": weight,
        "review_stratum": "common",
        "reviewer_1": "Reviewer A",
        "reviewer_1_date": "2026-07-14",
        "reviewer_1_ontology_terms": reviewer_1,
        "reviewer_1_context_sufficient": "yes",
        "reviewer_1_confidence": "high",
        "reviewer_1_notes": "",
        "reviewer_1_registry_exact_match": set(reviewer_1.split(";"))
        == set(truth.split(";")),
        "reviewer_2": "Reviewer B",
        "reviewer_2_date": "2026-07-14",
        "reviewer_2_ontology_terms": reviewer_2,
        "reviewer_2_context_sufficient": "yes",
        "reviewer_2_confidence": "high",
        "reviewer_2_notes": "",
        "reviewer_2_registry_exact_match": set(reviewer_2.split(";"))
        == set(truth.split(";")),
        "exact_term_set_agreement": agreement,
        "needs_adjudication": not agreement,
        "consensus_ontology_terms": reviewer_1 if agreement else "",
        "adjudicator": "",
        "adjudication_date": "",
        "adjudicated_ontology_terms": "",
        "adjudication_category": "",
        "adjudication_rationale": "",
    }
    return row


def create_source(root: Path) -> dict[str, Path]:
    yq = root / "bin" / "yq"
    yq.parent.mkdir(parents=True)
    fake_yq(yq)
    registry_path = root / "registry.yaml"
    write_json(registry_path, registry())
    protocol = MODULE.runtime.load_json(
        ROOT / "experiments" / "paper" / "protocol" / "ontology_mapping.json"
    )
    protocol["manual_review"]["sample_size"] = 4
    protocol["manual_review"]["analysis"]["bootstrap_replicates"] = 50
    protocol_path = root / "protocol.json"
    write_json(protocol_path, protocol)

    multi = "RGN:partition:sample;RGN:technical:index7"
    rows = [
        review_row(
            index=1,
            configuration="CONFIG1",
            truth="RGN:partition:cell",
            reviewer_1="RGN:partition:cell",
            reviewer_2="RGN:partition:cell",
            weight=1,
        ),
        review_row(
            index=2,
            configuration="CONFIG1",
            truth="RGN:partition:molecule",
            reviewer_1="RGN:partition:molecule",
            reviewer_2="RGN:partition:molecule",
            weight=2,
        ),
        review_row(
            index=3,
            configuration="CONFIG2",
            truth=multi,
            reviewer_1=multi,
            reviewer_2="RGN:partition:sample",
            weight=1,
        ),
        review_row(
            index=4,
            configuration="CONFIG3",
            truth="RGN:unknown:unclassified",
            reviewer_1="RGN:unknown:unclassified",
            reviewer_2="RGN:unknown:unclassified",
            weight=3,
        ),
    ]
    sample_fields = [
        "survey_id",
        "region_key",
        "configuration_accession",
        "ontology_terms",
        "sampling_weight",
        "review_stratum",
    ]
    fields = [*sample_fields, *MODULE.reviews.REVIEW_RESULT_FIELDS]
    template_path = root / "tables" / "ontology_reviews.csv"
    write_csv(template_path, rows, fields)
    completed = [{**row} for row in rows]
    disagreement = next(row for row in completed if row["needs_adjudication"])
    disagreement.update(
        {
            "adjudicator": "Reviewer C",
            "adjudication_date": "2026-07-15",
            "adjudicated_ontology_terms": multi,
            "adjudication_category": "missing_role",
            "adjudication_rationale": "The index has both direct roles.",
        }
    )
    completed_path = root / "tables" / "ontology_reviews.adjudicated.csv"
    write_csv(completed_path, completed, fields)
    manifest_path = root / "manifests" / "ontology_reviews.json"
    write_json(
        manifest_path,
        {
            "schema_version": "0.1.0",
            "survey_id": "survey-1",
            "review_rows": 4,
            "reviewers": {
                "reviewer_1": "Reviewer A",
                "reviewer_2": "Reviewer B",
            },
            "adjudication_required": 1,
            "adjudication_complete": False,
            "inputs": {
                "registry": MODULE.runtime.file_identity(registry_path),
                "protocol": MODULE.runtime.file_identity(protocol_path),
            },
            "outputs": {
                "ontology_reviews": MODULE.runtime.file_identity(template_path)
            },
        },
    )
    return {
        "manifest": manifest_path,
        "completed": completed_path,
        "template": template_path,
        "registry": registry_path,
        "protocol": protocol_path,
        "yq": yq,
    }


def analyze(
    root: Path, source: dict[str, Path], name: str = "analysis"
) -> dict[str, object]:
    return MODULE.analyze_reviews(
        reviews_manifest_path=source["manifest"],
        adjudicated_reviews_path=source["completed"],
        registry_path=source["registry"],
        protocol_path=source["protocol"],
        yq_bin=source["yq"],
        output_root=root / name,
        timeout_seconds=10,
    )


class AnalyzeOntologyReviewTests(unittest.TestCase):
    def test_analyzes_weighted_reviews_and_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = create_source(root / "source")
            first = analyze(root, source, "first")
            second = analyze(root, source, "second")

            self.assertEqual(first["analysis_id"], second["analysis_id"])
            self.assertTrue(first["valid"])
            self.assertEqual(first["counts"]["reviewed_regions"], 4)
            self.assertEqual(first["counts"]["human_adjudications"], 1)
            self.assertEqual(first["counts"]["registry_exact_matches"], 4)
            summary = {
                row["metric"]: row
                for row in read_csv(root / "first/tables/ontology_review_summary.csv")
            }
            self.assertEqual(
                float(summary["weighted_mapping_precision"]["estimate"]), 1
            )
            self.assertEqual(
                summary["weighted_mapping_precision"]["target_met"], "True"
            )
            self.assertLess(
                float(summary["weighted_exact_set_cohen_kappa"]["estimate"]), 1
            )
            outcomes = read_csv(root / "first/tables/ontology_review_outcomes.csv")
            adjudicated = next(
                row
                for row in outcomes
                if row["final_assignment_source"] == "human_adjudication"
            )
            self.assertEqual(adjudicated["registry_exact_match"], "True")
            self.assertEqual(
                (root / "first/bootstrap/ontology_review_bootstrap.csv").read_bytes(),
                (root / "second/bootstrap/ontology_review_bootstrap.csv").read_bytes(),
            )
            self.assertEqual(
                (root / "first/tables/ontology_review_summary.csv").read_bytes(),
                (root / "second/tables/ontology_review_summary.csv").read_bytes(),
            )

    def test_rejects_incomplete_tampered_and_invalid_adjudication(self) -> None:
        scenarios = {
            "incomplete": (lambda rows: rows[2].update(adjudicator=""), "incomplete"),
            "tampered": (
                lambda rows: rows[0].update(
                    reviewer_1_ontology_terms=("RGN:partition:molecule")
                ),
                "changed reviewer evidence",
            ),
            "invalid_category": (
                lambda rows: rows[2].update(adjudication_category="not_allowed"),
                "category is invalid",
            ),
            "unknown_mixed": (
                lambda rows: rows[2].update(
                    adjudicated_ontology_terms=(
                        "RGN:unknown:unclassified;RGN:partition:sample"
                    )
                ),
                "combines unknown",
            ),
        }
        for name, (mutate, pattern) in scenarios.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                root = Path(tmpdir)
                source = create_source(root / "source")
                rows = read_csv(source["completed"])
                fields = list(rows[0])
                mutate(rows)
                write_csv(source["completed"], rows, fields)
                with self.assertRaisesRegex(ValueError, pattern):
                    analyze(root, source)

    def test_cli_analyzes_completed_reviews(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = create_source(root / "source")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--reviews-manifest",
                    str(source["manifest"]),
                    "--adjudicated-reviews",
                    str(source["completed"]),
                    "--registry",
                    str(source["registry"]),
                    "--protocol",
                    str(source["protocol"]),
                    "--yq-bin",
                    str(source["yq"]),
                    "--output-root",
                    str(root / "analysis"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(
                (root / "analysis/manifests/ontology_review_analysis.json").is_file()
            )
            self.assertTrue(
                (root / "analysis/validation/ontology_reviews.json").is_file()
            )


if __name__ == "__main__":
    unittest.main()
