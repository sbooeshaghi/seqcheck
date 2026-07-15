import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import test_run_complementarity as execution_fixture


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "analyze_complementarity_calibration.py"
SOURCE_PROTOCOL = ROOT / "experiments" / "paper" / "protocol" / "complementarity.json"
SPEC = importlib.util.spec_from_file_location(
    "analyze_complementarity_calibration", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def create_execution(root: Path) -> tuple[dict[str, Path], dict[str, object]]:
    inputs = execution_fixture.create_materialization(root)
    manifest = execution_fixture.MODULE.run_complementarity(
        materialization_path=inputs["materialization"],
        protocol_path=inputs["protocol"],
        seqcheck_bin=inputs["seqcheck"],
        seqspec_bin=inputs["seqspec"],
        fastqc_bin=inputs["fastqc"],
        output_root=root / "execution",
        timeout_seconds=10,
    )
    return inputs, manifest


class AnalyzeComplementarityCalibrationTests(unittest.TestCase):
    def test_cli_calibrates_quality_and_rejects_changed_raw_fastqc(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs, execution = create_execution(root)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--execution-manifest",
                    str(execution["manifest_path"]),
                    "--protocol",
                    str(inputs["protocol"]),
                    "--output-root",
                    str(root / "analysis"),
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
                Path(manifest["outputs"]["quality_calls"]["path"])
            )
            severities = MODULE.runtime.read_csv(
                Path(manifest["outputs"]["severity_summary"]["path"])
            )
            invariance = MODULE.runtime.read_csv(
                Path(manifest["outputs"]["seqcheck_invariance"]["path"])
            )
            policy = MODULE.runtime.load_json(
                Path(manifest["outputs"]["quality_policy"]["path"])
            )

            self.assertTrue(manifest["valid"])
            self.assertEqual(len(calls), 4)
            self.assertTrue(all(value["detected"] == "True" for value in calls))
            self.assertEqual(len(severities), 4)
            self.assertTrue(all(value["eligible"] == "False" for value in severities))
            self.assertEqual(len(invariance), 4)
            self.assertTrue(all(value["passed"] == "True" for value in invariance))
            self.assertFalse(policy["frozen"])
            self.assertTrue(all(not value["ready"] for value in policy["entries"]))

            fastqc_runs = MODULE.runtime.read_csv(
                Path(execution["outputs"]["fastqc_runs"]["path"])
            )
            raw_path = Path(fastqc_runs[0]["data_path"])
            raw_path.write_text(
                raw_path.read_text(encoding="utf-8") + "\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "FastQC data hash changed"):
                MODULE.load_execution(Path(execution["manifest_path"]))

    def test_status_worsening_detects_quality_without_numeric_score(self) -> None:
        module_name = "Per base sequence quality"
        module_map = {
            ("clean", "case"): {
                module_name: {
                    "name": module_name,
                    "status": "pass",
                    "headers": [["Base", "Mean"]],
                    "rows": [],
                }
            },
            ("quality", "case"): {
                module_name: {
                    "name": module_name,
                    "status": "fail",
                    "headers": [["Base", "Mean"]],
                    "rows": [],
                }
            },
        }
        calls = MODULE.build_quality_calls(
            analysis_id="analysis",
            execution_id="execution",
            base_conditions=[
                {
                    "condition_id": "clean",
                    "condition_kind": "clean",
                    "configuration_accession": "config",
                }
            ],
            controls=[
                {
                    "condition_id": "quality",
                    "condition_kind": "quality",
                    "configuration_accession": "config",
                    "operator_id": "Q01",
                    "variant": "all_cycles_low",
                    "severity_index": 0,
                    "target": {"case_ids": ["case"]},
                    "expected": {"fastqc_modules": [module_name]},
                }
            ],
            module_map=module_map,
            protocol={"calibration": {"minimum_quality_score_decrease": 1.0}},
        )

        self.assertTrue(calls[0]["status_worsened"])
        self.assertTrue(calls[0]["detected"])
        self.assertIsNone(calls[0]["score_decrease"])
        self.assertIsNone(calls[0]["clean_score"])

    def test_policy_selects_least_severe_eligible_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            protocol_path = Path(tmpdir) / "protocol.json"
            protocol = MODULE.runtime.load_json(SOURCE_PROTOCOL)
            MODULE.runtime.write_json(protocol_path, protocol)
            severity_rows = []
            for operator in protocol["quality_controls"]:
                severity_rows.extend(
                    [
                        {
                            "operator_id": operator["operator_id"],
                            "severity_index": 0,
                            "eligible": False,
                        },
                        {
                            "operator_id": operator["operator_id"],
                            "severity_index": 1,
                            "eligible": True,
                            "independent_configurations": 6,
                            "detection_rate": 1.0,
                        },
                    ]
                )
            policy = MODULE.build_quality_policy(
                analysis_id="analysis",
                execution_id="execution",
                protocol_path=protocol_path,
                protocol=protocol,
                severity_rows=severity_rows,
            )

        self.assertTrue(policy["frozen"])
        self.assertTrue(
            all(value["severity_index"] == 1 for value in policy["entries"])
        )
        stable = {
            key: value
            for key, value in policy.items()
            if key not in {"policy_id", "created_at"}
        }
        self.assertEqual(policy["policy_id"], MODULE.runtime.sha256_json(stable)[:16])


if __name__ == "__main__":
    unittest.main()
