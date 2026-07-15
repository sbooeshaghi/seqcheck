import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import test_evaluate_perturbation_execution as policy_fixture
from tests import test_materialize_complementarity as materialization_fixture
from tests import test_run_complementarity as execution_fixture
from tests import test_run_perturbation_calibration as seqcheck_fixture


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "evaluate_complementarity.py"
SPEC = importlib.util.spec_from_file_location("evaluate_complementarity", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def create_evaluation_execution(root: Path) -> dict[str, object]:
    base = materialization_fixture.create_base(
        root / "study", cohort_split="evaluation"
    )
    protocol = materialization_fixture.create_protocol(root / "protocol.json")
    quality_policy = materialization_fixture.quality_policy(
        root / "quality_policy.json", protocol
    )
    materialization = materialization_fixture.MODULE.materialize_complementarity(
        base_materialization_path=base["base"],
        protocol_path=protocol,
        seqspec_bin=base["seqspec"],
        yq_bin=base["yq"],
        cohort_split="evaluation",
        quality_policy_path=quality_policy,
        output_root=root / "controls",
        timeout_seconds=10,
    )
    seqcheck = root / "seqcheck"
    fastqc = root / "fastqc"
    seqcheck_fixture.fake_seqcheck(seqcheck)
    execution_fixture.fake_fastqc(fastqc)
    execution = execution_fixture.MODULE.run_complementarity(
        materialization_path=Path(materialization["manifest_path"]),
        protocol_path=protocol,
        seqcheck_bin=seqcheck,
        seqspec_bin=base["seqspec"],
        fastqc_bin=fastqc,
        output_root=root / "execution",
        timeout_seconds=10,
    )
    _, _, base_conditions, controls = execution_fixture.MODULE.load_materialization(
        Path(materialization["manifest_path"])
    )
    detection_policy = policy_fixture.detection_policy(
        root / "detection_policy.json", base_conditions
    )
    return {
        "execution": execution,
        "protocol": protocol,
        "quality_policy": quality_policy,
        "detection_policy": detection_policy,
        "base_conditions": base_conditions,
        "controls": controls,
    }


class EvaluateComplementarityTests(unittest.TestCase):
    def test_cli_applies_frozen_policies_and_reconciles_endpoints(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fixture = create_evaluation_execution(root)
            execution = fixture["execution"]
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--execution-manifest",
                    str(execution["manifest_path"]),
                    "--quality-policy",
                    str(fixture["quality_policy"]),
                    "--detection-policy",
                    str(fixture["detection_policy"]),
                    "--protocol",
                    str(fixture["protocol"]),
                    "--output-root",
                    str(root / "evaluation"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest = json.loads(
                Path(completed.stdout.strip()).read_text(encoding="utf-8")
            )
            calls = MODULE.runtime.read_csv(
                Path(manifest["outputs"]["tool_calls"]["path"])
            )
            endpoints = MODULE.runtime.read_csv(
                Path(manifest["outputs"]["endpoints"]["path"])
            )
            repeated = MODULE.evaluate_complementarity(
                execution_path=Path(execution["manifest_path"]),
                quality_policy_path=Path(fixture["quality_policy"]),
                detection_policy_path=Path(fixture["detection_policy"]),
                protocol_path=Path(fixture["protocol"]),
                output_root=root / "evaluation-repeated",
            )

            conditions = [*fixture["base_conditions"], *fixture["controls"]]
            self.assertTrue(manifest["valid"])
            self.assertEqual(manifest["evaluation_id"], repeated["evaluation_id"])
            self.assertEqual(len(calls), len(conditions) * 3)
            self.assertEqual(len(endpoints), 5)
            quality_ids = {
                value["condition_id"]
                for value in fixture["controls"]
                if value["condition_kind"] == "quality"
            }
            quality_fastqc = [
                value
                for value in calls
                if value["condition_id"] in quality_ids and value["tool"] == "fastqc"
            ]
            quality_seqcheck = [
                value
                for value in calls
                if value["condition_id"] in quality_ids and value["tool"] == "seqcheck"
            ]
            invalid_ids = {
                value["condition_id"]
                for value in fixture["controls"]
                if value["condition_kind"] == "schema_invalid"
            }
            invalid_seqspec = [
                value
                for value in calls
                if value["condition_id"] in invalid_ids and value["tool"] == "seqspec"
            ]
            self.assertTrue(
                all(value["detected"] == "True" for value in quality_fastqc)
            )
            self.assertTrue(
                all(value["detected"] == "False" for value in quality_seqcheck)
            )
            self.assertTrue(
                all(value["detected"] == "True" for value in invalid_seqspec)
            )
            primary = MODULE.primary_structural_conditions(
                fixture["base_conditions"],
                MODULE.runtime.load_json(Path(fixture["protocol"])),
            )
            advantage = next(
                value
                for value in endpoints
                if value["endpoint"] == "seqcheck_structural_advantage"
            )
            self.assertEqual(int(advantage["observations"]), len(primary))

            alternate_policy = root / "alternate_quality_policy.json"
            policy = MODULE.runtime.load_json(Path(fixture["quality_policy"]))
            policy["created_at"] = "2026-07-15T00:00:00+00:00"
            MODULE.runtime.write_json(alternate_policy, policy)
            with self.assertRaisesRegex(ValueError, "differs from materialization"):
                MODULE.evaluate_complementarity(
                    execution_path=Path(execution["manifest_path"]),
                    quality_policy_path=alternate_policy,
                    detection_policy_path=Path(fixture["detection_policy"]),
                    protocol_path=Path(fixture["protocol"]),
                    output_root=root / "rejected",
                )

    def test_primary_scope_uses_deterministic_and_stochastic_anchor(self) -> None:
        protocol = MODULE.runtime.load_json(
            ROOT / "experiments" / "paper" / "protocol" / "complementarity.json"
        )
        common = {"expected_seqspec_check": "pass"}
        conditions = [
            {**common, "condition_kind": "clean", "event_fraction": None},
            {**common, "condition_kind": "deterministic", "event_fraction": None},
            {**common, "condition_kind": "stochastic", "event_fraction": 0.01},
            {**common, "condition_kind": "stochastic", "event_fraction": 0.1},
            {
                "condition_kind": "deterministic",
                "event_fraction": None,
                "expected_seqspec_check": "failure",
            },
        ]

        selected = MODULE.primary_structural_conditions(conditions, protocol)

        self.assertEqual(selected, conditions[1:3])


if __name__ == "__main__":
    unittest.main()
