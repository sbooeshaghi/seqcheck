import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import test_materialize_perturbations as materialize_fixture
from tests import test_run_perturbation_calibration as execution_fixture


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "analyze_ontology_utility.py"
SPEC = importlib.util.spec_from_file_location("analyze_ontology_utility", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
PROTOCOL_PATH = ROOT / "experiments" / "paper" / "protocol" / "ontology_utility.json"


def protocol() -> dict[str, object]:
    return MODULE.runtime.load_json(PROTOCOL_PATH)


def condition(
    condition_id: str,
    *,
    configuration: str = "CONFIG1",
    kind: str = "clean",
    operator_id: str = "CLEAN",
    event_fraction: float | None = None,
    region_id: str = "umi",
    metric_name: str = "sequence_entropy_fraction",
) -> dict[str, object]:
    return {
        "condition_id": condition_id,
        "configuration_accession": configuration,
        "family_id": "rna_family",
        "modality": "rna",
        "condition_kind": kind,
        "operator_id": operator_id,
        "event_fraction": event_fraction,
        "target": {"region_ids": [] if kind == "clean" else [region_id]},
        "expected": {"metric_names": [] if kind == "clean" else [metric_name]},
    }


def metric(
    condition_id: str,
    value: float,
    *,
    ontology_terms: list[str] | None = None,
) -> dict[str, str]:
    terms = ontology_terms or ["RGN:partition:molecule"]
    annotation = [
        {
            "region_id": "umi",
            "sequence_type": "random",
            "ontology_terms": terms,
        }
    ]
    return {
        "condition_id": condition_id,
        "metric_side": "observed",
        "data_kind": "scalar",
        "check": "random",
        "metric_name": "sequence_entropy_fraction",
        "unit": "fraction",
        "files": "R1.fastq.gz",
        "reads": "rna_R1",
        "regions": "umi",
        "ontology_terms": ";".join(sorted(terms)),
        "sequence_types": "random",
        "region_annotations_json": MODULE.runtime.canonical_json(annotation),
        "value_json": json.dumps(value),
    }


def synthetic_observation(
    *,
    split: str,
    configuration: str,
    term: str,
    label: int,
    value: float,
) -> dict[str, object]:
    result = MODULE.finalize_observation(
        {
            "split": split,
            "configuration_accession": configuration,
            "family_id": "rna_family",
            "modality": "rna",
            "operator_id": "D00" if label else "CLEAN",
            "condition_count": 1,
            "condition_ids": f"{split}-{configuration}-{term}-{label}",
            "condition_class": "targeted_perturbation" if label else "clean",
            "label": label,
            "check": "random",
            "metric_name": "sequence_entropy_fraction",
            "unit": "fraction",
            "files": "R1.fastq.gz",
            "reads": "rna_R1",
            "region_id": f"region-{term}",
            "sequence_type": "random",
            "ontology_signature": term,
            "ontology_terms": term,
            "value": value,
        }
    )
    result["eligible"] = True
    result["analysis_id"] = "analysis"
    return result


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def create_execution(
    root: Path, *, evaluation: bool, configuration_accession: str
) -> Path:
    inputs = materialize_fixture.create_inputs(
        root / "inputs",
        cohort_split="evaluation" if evaluation else "calibration",
        configuration_accession=configuration_accession,
    )
    kwargs = {}
    if evaluation:
        kwargs["sample_bundle_path"] = materialize_fixture.create_sample_bundle(
            inputs, root / "bundle" / "bundle.json"
        )
    else:
        kwargs["study_manifest_path"] = inputs["study"]
    materialization = materialize_fixture.MODULE.materialize_perturbations(
        inventory_manifest_path=inputs["inventory"],
        sampling_policy_path=inputs["policy"],
        perturbation_protocol_path=inputs["protocol"],
        seqspec_bin=inputs["seqspec"],
        yq_bin=inputs["yq"],
        output_root=root / "materialized",
        timeout_seconds=10,
        **kwargs,
    )
    seqcheck = root / "seqcheck"
    execution_fixture.fake_seqcheck(seqcheck)
    execution = execution_fixture.MODULE.run_perturbation_calibration(
        materialization_manifest_path=Path(materialization["manifest_path"]),
        execution_protocol_path=(
            ROOT / "experiments" / "paper" / "protocol" / "perturbation_execution.json"
        ),
        seqcheck_bin=seqcheck,
        output_root=root / "execution",
        timeout_seconds=10,
    )
    return Path(execution["manifest_path"])


class AnalyzeOntologyUtilityTests(unittest.TestCase):
    def test_protocol_rejects_a_model_that_does_not_match_its_label(self) -> None:
        value = protocol()
        value["models"][2]["group_fields"] = ["metric_name", "sequence_type"]

        with self.assertRaisesRegex(ValueError, "model group fields are invalid"):
            MODULE.validate_protocol(value)

    def test_observations_aggregate_anchor_seeds_by_configuration(self) -> None:
        conditions = [condition("clean")]
        metrics = [metric("clean", 0.9)]
        for seed, value in ((101, 0.5), (211, 0.7)):
            condition_id = f"positive-{seed}"
            conditions.append(
                condition(
                    condition_id,
                    kind="stochastic",
                    operator_id="D00",
                    event_fraction=0.01,
                )
            )
            metrics.append(metric(condition_id, value))
        conditions.append(
            condition(
                "other-fraction",
                kind="stochastic",
                operator_id="D00",
                event_fraction=0.1,
            )
        )
        metrics.append(metric("other-fraction", 0.1))

        observations = MODULE.build_observations(
            split="calibration",
            conditions=conditions,
            metrics=metrics,
            protocol=protocol(),
        )

        self.assertEqual(len(observations), 2)
        positive = next(value for value in observations if value["label"] == 1)
        self.assertEqual(positive["condition_count"], 2)
        self.assertAlmostEqual(positive["value"], 0.6)
        self.assertEqual(positive["ontology_signature"], "RGN:partition:molecule")

    def test_observations_reject_report_and_seqspec_ontology_disagreement(self) -> None:
        row = metric("clean", 0.9)
        row["ontology_terms"] = "RGN:measure:transcript"

        with self.assertRaisesRegex(ValueError, "ontology terms differ"):
            MODULE.build_observations(
                split="calibration",
                conditions=[condition("clean")],
                metrics=[row],
                protocol=protocol(),
            )

    def test_ontology_model_removes_role_pooling_false_positives(self) -> None:
        observations = []
        terms = (
            ("RGN:partition:molecule", 0.95, 0.75),
            ("RGN:measure:transcript", 0.45, 0.25),
        )
        for split, prefix in (("calibration", "C"), ("evaluation", "E")):
            for term_index, (term, clean_value, positive_value) in enumerate(terms):
                for index in range(3):
                    configuration = f"{prefix}{term_index}{index}"
                    observations.append(
                        synthetic_observation(
                            split=split,
                            configuration=configuration,
                            term=term,
                            label=0,
                            value=clean_value,
                        )
                    )
                    observations.append(
                        synthetic_observation(
                            split=split,
                            configuration=configuration,
                            term=term,
                            label=1,
                            value=positive_value,
                        )
                    )
        value = protocol()
        value["statistics"]["bootstrap_replicates"] = 100
        _, scores, _ = MODULE.fit_and_score_models(observations, value, "analysis")
        comparison, bootstrap = MODULE.compare_primary_models(
            scores=scores,
            protocol=value,
            analysis_id="analysis",
        )

        self.assertEqual(len(bootstrap), 100)
        self.assertEqual(comparison["primary_false_positive_rate"], 0.0)
        self.assertEqual(comparison["baseline_false_positive_rate"], 0.5)
        self.assertEqual(comparison["relative_false_positive_rate_reduction"], 1.0)
        self.assertTrue(comparison["scientific_target_met"])

    def test_cli_runs_locked_disjoint_execution_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            calibration = create_execution(
                root / "calibration",
                evaluation=False,
                configuration_accession="CAL1",
            )
            evaluation = create_execution(
                root / "evaluation",
                evaluation=True,
                configuration_accession="EVAL1",
            )
            value = protocol()
            value["eligibility"][
                "minimum_independent_configurations_per_term_per_split"
            ] = 1
            for model in value["models"]:
                model["minimum_reference_configurations"] = 1
            value["statistics"]["bootstrap_replicates"] = 20
            analysis_protocol = root / "ontology_utility.json"
            MODULE.runtime.write_json(analysis_protocol, value)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--calibration-execution",
                    str(calibration),
                    "--evaluation-execution",
                    str(evaluation),
                    "--analysis-protocol",
                    str(analysis_protocol),
                    "--output-root",
                    str(root / "analysis"),
                ],
                capture_output=True,
                text=True,
            )
            manifest = MODULE.runtime.load_json(Path(completed.stdout.strip()))
            comparison = read_csv(manifest["outputs"]["primary_comparison"]["path"])

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(manifest["valid"])
        self.assertEqual(len(comparison), 1)
        self.assertEqual(manifest["counts"]["calibration_configurations"], 1)
        self.assertEqual(manifest["counts"]["evaluation_configurations"], 1)


if __name__ == "__main__":
    unittest.main()
