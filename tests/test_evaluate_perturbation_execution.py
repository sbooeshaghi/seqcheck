import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import test_run_perturbation_calibration as execution_fixture


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "evaluate_perturbation_execution.py"
PROTOCOL_PATH = (
    ROOT / "experiments" / "paper" / "protocol" / "perturbation_analysis.json"
)
SPEC = importlib.util.spec_from_file_location(
    "evaluate_perturbation_execution", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def detection_policy(path: Path, conditions: list[dict[str, object]]) -> Path:
    entries = []
    keys = sorted(
        {
            (str(value["operator_id"]), str(value["variant"]))
            for value in conditions
            if value["operator_id"] != "CLEAN"
        }
    )
    for operator_id, variant in keys:
        process_failure = operator_id == "S01"
        entries.append(
            {
                "operator_id": operator_id,
                "variant": variant,
                "condition_kind": next(
                    value["condition_kind"]
                    for value in conditions
                    if value["operator_id"] == operator_id
                    and value["variant"] == variant
                ),
                "calibration_configurations": 6,
                "process_failure": process_failure,
                "stderr_patterns": ["could not match"] if process_failure else [],
                "assessment_codes": [],
                "metric_endpoint": (
                    None
                    if process_failure
                    else {
                        "check": "fixed",
                        "metric_name": "exact_match_fraction",
                        "unit": "fraction",
                        "direction": "decrease",
                        "threshold": 0.1,
                        "direction_consistency": 1.0,
                        "p10_absolute_effect": 0.2,
                    }
                ),
                "ready": True,
            }
        )
    protocol = MODULE.runtime.load_json(PROTOCOL_PATH)
    stable = {
        "schema_version": "0.1.0",
        "analysis_id": "calibration-analysis",
        "execution_id": "different-calibration-execution",
        "frozen": True,
        "calibration": protocol["calibration"],
        "entries": entries,
    }
    MODULE.runtime.write_json(
        path,
        {
            **stable,
            "policy_id": MODULE.runtime.sha256_json(stable)[:16],
            "created_at": "2026-07-14T00:00:00+00:00",
        },
    )
    return path


class EvaluatePerturbationExecutionTests(unittest.TestCase):
    def test_cli_applies_frozen_policy_without_refitting(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = execution_fixture.create_materialization(root, evaluation=True)
            execution = execution_fixture.MODULE.run_perturbation_calibration(
                materialization_manifest_path=inputs["materialization"],
                execution_protocol_path=inputs["protocol"],
                seqcheck_bin=inputs["seqcheck"],
                output_root=root / "execution",
                timeout_seconds=10,
            )
            materialization = MODULE.runtime.load_json(inputs["materialization"])
            conditions = MODULE.runtime.load_json(
                Path(materialization["outputs"]["conditions"]["path"])
            )["conditions"]
            policy_path = detection_policy(root / "detection_policy.json", conditions)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--execution-manifest",
                    execution["manifest_path"],
                    "--detection-policy",
                    str(policy_path),
                    "--analysis-protocol",
                    str(PROTOCOL_PATH),
                    "--output-root",
                    str(root / "evaluation"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            manifest = json.loads(
                Path(completed.stdout.strip()).read_text(encoding="utf-8")
            )
            calls = read_csv(manifest["outputs"]["calls"]["path"])
            endpoints = read_csv(manifest["outputs"]["endpoints"]["path"])
            validation = MODULE.runtime.load_json(
                Path(manifest["outputs"]["validation"]["path"])
            )
            repeated = MODULE.evaluate_perturbation_execution(
                execution_manifest_path=Path(execution["manifest_path"]),
                detection_policy_path=policy_path,
                analysis_protocol_path=PROTOCOL_PATH,
                output_root=root / "evaluation-repeated",
            )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(manifest["valid"])
        self.assertEqual(manifest["evaluation_id"], repeated["evaluation_id"])
        self.assertEqual(len(calls), manifest["counts"]["conditions"])
        self.assertEqual(len(endpoints), 5)
        deterministic = next(
            value
            for value in endpoints
            if value["endpoint"] == "deterministic_sensitivity"
        )
        expected_deterministic = sum(
            value["condition_kind"] == "deterministic"
            and value["expected_seqspec_check"] != "failure"
            for value in conditions
        )
        self.assertEqual(int(deterministic["observations"]), expected_deterministic)
        self.assertFalse(validation["scientific_endpoints_complete"])
        s01 = next(value for value in calls if value["operator_id"] == "S01")
        self.assertEqual(s01["detected"], "True")
        self.assertEqual(s01["detection_source"], "process")

    def test_policy_coverage_rejects_an_unseen_evaluation_variant(self) -> None:
        policy = {
            "entries": [{"operator_id": "S04", "variant": "exclude_observed_length"}]
        }
        conditions = [
            {
                "operator_id": "S05",
                "variant": "shift_1",
            }
        ]
        with self.assertRaisesRegex(ValueError, "no rule"):
            MODULE.validate_policy_coverage(policy, conditions)

    def test_policy_loader_rejects_a_partially_missing_identifier(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "policy.json"
            stable = {
                "schema_version": "0.1.0",
                "analysis_id": "calibration-analysis",
                "execution_id": "calibration-execution",
                "frozen": True,
                "calibration": {},
                "entries": [
                    {
                        "operator_id": "S04",
                        "variant": None,
                        "ready": True,
                    }
                ],
            }
            MODULE.runtime.write_json(
                path,
                {
                    **stable,
                    "policy_id": MODULE.runtime.sha256_json(stable)[:16],
                    "created_at": "2026-07-14T00:00:00+00:00",
                },
            )

            with self.assertRaisesRegex(ValueError, "identifiers are invalid"):
                MODULE.load_detection_policy(path)

    def test_configuration_cluster_interval_is_deterministic(self) -> None:
        first = MODULE.bootstrap_mean_interval(
            [0.0, 0.5, 1.0], confidence=0.95, resamples=1000, seed=17
        )
        second = MODULE.bootstrap_mean_interval(
            [0.0, 0.5, 1.0], confidence=0.95, resamples=1000, seed=17
        )
        self.assertEqual(first, second)
        self.assertLessEqual(first[0], 0.5)
        self.assertGreaterEqual(first[1], 0.5)


if __name__ == "__main__":
    unittest.main()
