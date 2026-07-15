import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import test_materialize_perturbations as materialize_fixture


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "run_perturbation_calibration.py"
SPEC = importlib.util.spec_from_file_location(
    "run_perturbation_calibration", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fake_seqcheck(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import json
import sys
from pathlib import Path

if sys.argv[1] == "--version":
    print("seqcheck 0.2.0-test")
    raise SystemExit(0)
inputs = [Path(value) for value in sys.argv[sys.argv.index("--output") + 2:]]
if any(value.name == "seqcheck_unexpected_input.fastq.gz" for value in inputs):
    print("could not match unexpected FASTQ", file=sys.stderr)
    raise SystemExit(2)
output = Path(sys.argv[sys.argv.index("--output") + 1])
payload = {
    "report_schema_version": "0.1.0",
    "meta": {"command": "check", "requested_reads": 0},
    "results": [
        {
            "check": "fixed",
            "files": [inputs[0].name],
            "reads": ["synthetic_R1"],
            "regions": ["linker"],
            "ontology": ["RGN:technical:linker"],
            "expected": [
                {
                    "id": "e1",
                    "name": "expected_sequence",
                    "description": "Expected fixed sequence.",
                    "data": {"kind": "text", "value": "TT", "unit": "sequence"},
                }
            ],
            "observed": [
                {
                    "id": "o1",
                    "name": "exact_match_fraction",
                    "description": "Exact match fraction.",
                    "data": {"kind": "scalar", "value": 0.5, "unit": "fraction"},
                }
            ],
            "assessment": [
                {
                    "type": "interpretation",
                    "code": "fixed_exact_match_partial",
                    "description": "Some reads match.",
                    "expected_ids": ["e1"],
                    "observed_ids": ["o1"],
                }
            ],
        }
    ],
}
output.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def create_materialization(root: Path) -> dict[str, Path]:
    inputs = materialize_fixture.create_inputs(root / "inputs")
    materialization = materialize_fixture.MODULE.materialize_perturbations(
        inventory_manifest_path=inputs["inventory"],
        study_manifest_path=inputs["study"],
        sampling_policy_path=inputs["policy"],
        perturbation_protocol_path=inputs["protocol"],
        seqspec_bin=inputs["seqspec"],
        yq_bin=inputs["yq"],
        output_root=root / "materialized",
        timeout_seconds=10,
    )
    seqcheck = root / "seqcheck"
    fake_seqcheck(seqcheck)
    return {
        "materialization": Path(materialization["manifest_path"]),
        "conditions": Path(materialization["outputs"]["conditions"]["path"]),
        "seqcheck": seqcheck,
        "protocol": ROOT
        / "experiments"
        / "paper"
        / "protocol"
        / "perturbation_execution.json",
    }


class RunPerturbationCalibrationTests(unittest.TestCase):
    def test_runner_exports_reports_metrics_assessments_and_failures(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_materialization(root)
            manifest = MODULE.run_perturbation_calibration(
                materialization_manifest_path=inputs["materialization"],
                execution_protocol_path=inputs["protocol"],
                seqcheck_bin=inputs["seqcheck"],
                output_root=root / "execution",
                timeout_seconds=10,
            )
            validation = MODULE.runtime.load_json(
                Path(manifest["outputs"]["validation"]["path"])
            )
            runs = read_csv(manifest["outputs"]["runs"]["path"])
            metrics = read_csv(manifest["outputs"]["metrics"]["path"])
            assessments = read_csv(manifest["outputs"]["assessments"]["path"])
            performance = read_csv(manifest["outputs"]["performance"]["path"])

        self.assertTrue(validation["valid"])
        self.assertEqual(validation["counts"]["runs"], 13)
        self.assertEqual(validation["counts"]["successful_runs"], 12)
        self.assertEqual(validation["counts"]["expected_failure_runs"], 1)
        self.assertEqual(len(runs), 13)
        self.assertEqual(len(metrics), 24)
        self.assertEqual(len(assessments), 12)
        self.assertEqual(len(performance), 13)
        failure = next(value for value in runs if value["operator_id"] == "S01")
        self.assertEqual(failure["observed_process_outcome"], "failure")
        self.assertEqual(failure["report_path"], "")
        self.assertTrue(
            all(value["ontology_terms"] == "RGN:technical:linker" for value in metrics)
        )
        self.assertTrue(
            all(
                value["assessment_code"] == "fixed_exact_match_partial"
                for value in assessments
            )
        )

    def test_runner_rejects_changed_condition_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_materialization(root)
            payload = MODULE.runtime.load_json(inputs["conditions"])
            artifact = next(
                Path(condition["inputs"][0]["path"])
                for condition in payload["conditions"]
                if condition["condition_kind"] == "stochastic"
            )
            artifact.write_bytes(artifact.read_bytes() + b"changed")

            with self.assertRaisesRegex(ValueError, "artifact hash changed"):
                MODULE.run_perturbation_calibration(
                    materialization_manifest_path=inputs["materialization"],
                    execution_protocol_path=inputs["protocol"],
                    seqcheck_bin=inputs["seqcheck"],
                    output_root=root / "execution",
                    timeout_seconds=10,
                )

    def test_cli_writes_valid_execution_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_materialization(root)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--materialization-manifest",
                    str(inputs["materialization"]),
                    "--execution-protocol",
                    str(inputs["protocol"]),
                    "--seqcheck-bin",
                    str(inputs["seqcheck"]),
                    "--output-root",
                    str(root / "execution"),
                    "--timeout-seconds",
                    "10",
                ],
                capture_output=True,
                text=True,
            )
            manifest = MODULE.runtime.load_json(Path(completed.stdout.strip()))

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(manifest["valid"])


if __name__ == "__main__":
    unittest.main()
