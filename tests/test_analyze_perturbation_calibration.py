import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import test_run_perturbation_calibration as execution_fixture


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "analyze_perturbation_calibration.py"
SPEC = importlib.util.spec_from_file_location(
    "analyze_perturbation_calibration", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
PROTOCOL_PATH = (
    ROOT / "experiments" / "paper" / "protocol" / "perturbation_analysis.json"
)


def protocol() -> dict[str, object]:
    return MODULE.runtime.load_json(PROTOCOL_PATH)


def condition(
    configuration: int,
    *,
    operator_id: str,
    variant: str,
    kind: str = "deterministic",
    event_fraction: float | None = None,
    seed: int | None = None,
    expected: dict[str, object] | None = None,
    region_ids: list[str] | None = None,
) -> dict[str, object]:
    fastq = f"FASTQ{configuration}.fastq"
    condition_id = f"{operator_id}-{variant}-{configuration}-{event_fraction}-{seed}"
    return {
        "condition_id": condition_id,
        "configuration_accession": f"CONFIG{configuration}",
        "family_id": "rna_family",
        "modality": "rna",
        "operator_id": operator_id,
        "variant": variant,
        "condition_kind": kind,
        "event_fraction": event_fraction,
        "mutation_seed": seed,
        "target": {
            "fastq_accessions": [fastq],
            "read_ids": ["read1"],
            "region_ids": region_ids or [],
        },
        "expected": expected
        or {
            "process_outcome": "success",
            "metric_names": [],
            "assessment_codes": [],
            "localization": ["file", "read"],
        },
        "inputs": [{"path": f"/tmp/{fastq}.gz"}],
    }


def clean_condition(configuration: int) -> dict[str, object]:
    return condition(
        configuration,
        operator_id="CLEAN",
        variant="clean",
        kind="clean",
        expected={
            "process_outcome": "success",
            "metric_names": [],
            "assessment_codes": [],
            "localization": [],
        },
    )


def metric_row(
    value: dict[str, object],
    *,
    metric_value: float,
    metric_name: str = "out_of_range_fraction",
    check: str = "length",
    unit: str = "fraction",
    regions: str = "",
) -> dict[str, str]:
    configuration = int(str(value["configuration_accession"]).replace("CONFIG", ""))
    return {
        "condition_id": str(value["condition_id"]),
        "metric_side": "observed",
        "check": check,
        "files": f"FASTQ{configuration}.fastq.gz",
        "reads": "read1",
        "regions": regions,
        "metric_name": metric_name,
        "unit": unit,
        "value_json": str(metric_value),
    }


def success_runs(conditions: list[dict[str, object]]) -> dict[str, dict[str, str]]:
    return {
        str(value["condition_id"]): {
            "observed_process_outcome": "success",
            "stderr_path": "",
            "stderr_sha256": "",
        }
        for value in conditions
    }


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class AnalyzePerturbationCalibrationTests(unittest.TestCase):
    def test_metric_policy_requires_and_detects_six_configurations(self) -> None:
        conditions = []
        metrics = []
        for configuration in range(6):
            clean = clean_condition(configuration)
            perturbed = condition(
                configuration,
                operator_id="S04",
                variant="exclude_observed_length",
                expected={
                    "process_outcome": "success",
                    "metric_names": ["out_of_range_fraction"],
                    "assessment_codes": [],
                    "localization": ["file", "read"],
                },
            )
            conditions.extend((clean, perturbed))
            metrics.extend(
                (
                    metric_row(clean, metric_value=0.0),
                    metric_row(perturbed, metric_value=0.2 + configuration / 100),
                )
            )

        effects = MODULE.build_effect_rows(
            analysis_id="analysis",
            execution_id="execution",
            conditions=conditions,
            metrics=metrics,
        )
        endpoints = MODULE.calibrate_metric_endpoints(
            analysis_id="analysis",
            conditions=conditions,
            effects=effects,
            protocol=protocol(),
        )
        evidence = MODULE.build_assessment_evidence(
            conditions=conditions, assessments=[]
        )
        runs = success_runs(conditions)
        policy = MODULE.build_policy(
            analysis_id="analysis",
            execution={"execution_id": "execution"},
            conditions=conditions,
            runs=runs,
            endpoints=endpoints,
            assessment_evidence=evidence,
            protocol=protocol(),
        )
        calls = MODULE.call_conditions(
            analysis_id="analysis",
            execution_id="execution",
            conditions=conditions,
            runs=runs,
            effects=effects,
            assessment_evidence=evidence,
            policy=policy,
        )

        endpoint = next(value for value in endpoints if value["eligible"])
        self.assertEqual(endpoint["direction"], "increase")
        self.assertEqual(endpoint["configurations"], 6)
        self.assertAlmostEqual(endpoint["threshold"], 0.1)
        self.assertTrue(policy["frozen"])
        perturbed_calls = [value for value in calls if value["operator_id"] == "S04"]
        self.assertTrue(all(value["detected"] for value in perturbed_calls))
        self.assertTrue(all(value["localized"] for value in perturbed_calls))
        self.assertTrue(
            all(
                not value["localized"]
                for value in calls
                if value["operator_id"] == "CLEAN"
            )
        )
        diagnostics = MODULE.scientific_diagnostics(
            calls=calls,
            monotonic=[],
            protocol=protocol(),
            policy_frozen=policy["frozen"],
        )
        self.assertEqual(
            diagnostics["measures"]["deterministic_sensitivity"]["value"], 1.0
        )
        self.assertEqual(diagnostics["measures"]["clean_specificity"]["value"], 1.0)
        self.assertEqual(diagnostics["status"], "insufficient_applicable_conditions")

    def test_assessment_policy_requires_support_within_each_configuration(self) -> None:
        conditions = []
        assessments = []
        expected = {
            "process_outcome": "success",
            "metric_names": [],
            "assessment_codes": ["missing_onlist_resource"],
            "localization": ["file", "read", "region"],
        }
        for configuration in range(6):
            conditions.append(clean_condition(configuration))
            for seed in (101, 211, 307):
                perturbed = condition(
                    configuration,
                    operator_id="D03",
                    variant="offlist_substitution",
                    kind="stochastic",
                    event_fraction=0.01,
                    seed=seed,
                    expected=expected,
                    region_ids=["barcode"],
                )
                conditions.append(perturbed)
                if seed == 101:
                    assessments.append(
                        {
                            "condition_id": perturbed["condition_id"],
                            "assessment_code": "missing_onlist_resource",
                            "files": f"FASTQ{configuration}.fastq",
                            "reads": "read1",
                            "regions": "barcode",
                        }
                    )

        evidence = MODULE.build_assessment_evidence(
            conditions=conditions, assessments=assessments
        )
        policy = MODULE.build_policy(
            analysis_id="analysis",
            execution={"execution_id": "execution"},
            conditions=conditions,
            runs=success_runs(conditions),
            endpoints=[],
            assessment_evidence=evidence,
            protocol=protocol(),
        )

        self.assertEqual(policy["entries"][0]["assessment_codes"], [])
        self.assertFalse(policy["entries"][0]["ready"])
        for configuration in range(6):
            for seed in (211, 307):
                value = next(
                    item
                    for item in conditions
                    if item["configuration_accession"] == f"CONFIG{configuration}"
                    and item["mutation_seed"] == seed
                )
                evidence[str(value["condition_id"])] = {
                    "new_codes": ["missing_onlist_resource"],
                    "localized_codes": ["missing_onlist_resource"],
                }
        supported = MODULE.build_policy(
            analysis_id="analysis",
            execution={"execution_id": "execution"},
            conditions=conditions,
            runs=success_runs(conditions),
            endpoints=[],
            assessment_evidence=evidence,
            protocol=protocol(),
        )
        self.assertEqual(
            supported["entries"][0]["assessment_codes"],
            ["missing_onlist_resource"],
        )
        self.assertTrue(supported["frozen"])

    def test_process_detection_requires_matching_hashed_stderr_and_filename(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            conditions = []
            runs = {}
            expected = {
                "process_outcome": "failure",
                "stderr_patterns": ["could not match"],
                "metric_names": [],
                "assessment_codes": [],
                "localization": ["file"],
            }
            for configuration in range(6):
                clean = clean_condition(configuration)
                perturbed = condition(
                    configuration,
                    operator_id="S01",
                    variant="unexpected_name",
                    expected=expected,
                )
                input_path = root / f"renamed-{configuration}.fastq.gz"
                perturbed["inputs"] = [{"path": str(input_path)}]
                stderr = root / f"stderr-{configuration}.txt"
                stderr.write_text(
                    f"error: could not match '{input_path.name}' to any read\n",
                    encoding="utf-8",
                )
                conditions.extend((clean, perturbed))
                runs[str(clean["condition_id"])] = {
                    "observed_process_outcome": "success",
                    "stderr_path": "",
                    "stderr_sha256": "",
                }
                runs[str(perturbed["condition_id"])] = {
                    "observed_process_outcome": "failure",
                    "stderr_path": str(stderr),
                    "stderr_sha256": MODULE.runtime.file_sha256(stderr),
                }
            evidence = MODULE.build_assessment_evidence(
                conditions=conditions, assessments=[]
            )
            policy = MODULE.build_policy(
                analysis_id="analysis",
                execution={"execution_id": "execution"},
                conditions=conditions,
                runs=runs,
                endpoints=[],
                assessment_evidence=evidence,
                protocol=protocol(),
            )
            calls = MODULE.call_conditions(
                analysis_id="analysis",
                execution_id="execution",
                conditions=conditions,
                runs=runs,
                effects=[],
                assessment_evidence=evidence,
                policy=policy,
            )
            perturbed_calls = [
                value for value in calls if value["operator_id"] == "S01"
            ]
            self.assertTrue(policy["frozen"])
            self.assertTrue(all(value["process_detected"] for value in perturbed_calls))
            self.assertTrue(all(value["localized"] for value in perturbed_calls))

            stderr_path = Path(runs[str(conditions[1]["condition_id"])]["stderr_path"])
            stderr_path.write_text("changed", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "stderr hash changed"):
                MODULE.call_conditions(
                    analysis_id="analysis",
                    execution_id="execution",
                    conditions=conditions,
                    runs=runs,
                    effects=[],
                    assessment_evidence=evidence,
                    policy=policy,
                )

    def test_stochastic_metric_response_is_monotonic_after_seed_collapse(self) -> None:
        conditions = []
        metrics = []
        expected = {
            "process_outcome": "success",
            "metric_names": ["exact_match_fraction"],
            "assessment_codes": [],
            "localization": ["file", "read", "region"],
        }
        for configuration in range(6):
            clean = clean_condition(configuration)
            conditions.append(clean)
            metrics.append(
                metric_row(
                    clean,
                    metric_value=1.0,
                    metric_name="exact_match_fraction",
                    check="fixed",
                    regions="linker",
                )
            )
            for fraction in (0.01, 0.05, 0.1):
                for seed in (101, 211, 307):
                    perturbed = condition(
                        configuration,
                        operator_id="D02",
                        variant="fixed_substitution",
                        kind="stochastic",
                        event_fraction=fraction,
                        seed=seed,
                        expected=expected,
                        region_ids=["linker"],
                    )
                    conditions.append(perturbed)
                    metrics.append(
                        metric_row(
                            perturbed,
                            metric_value=1.0 - fraction,
                            metric_name="exact_match_fraction",
                            check="fixed",
                            regions="linker",
                        )
                    )
        effects = MODULE.build_effect_rows(
            analysis_id="analysis",
            execution_id="execution",
            conditions=conditions,
            metrics=metrics,
        )
        endpoints = MODULE.calibrate_metric_endpoints(
            analysis_id="analysis",
            conditions=conditions,
            effects=effects,
            protocol=protocol(),
        )
        evidence = MODULE.build_assessment_evidence(
            conditions=conditions, assessments=[]
        )
        policy = MODULE.build_policy(
            analysis_id="analysis",
            execution={"execution_id": "execution"},
            conditions=conditions,
            runs=success_runs(conditions),
            endpoints=endpoints,
            assessment_evidence=evidence,
            protocol=protocol(),
        )
        monotonic = MODULE.monotonic_rows(
            analysis_id="analysis",
            effects=effects,
            policy=policy,
            epsilon=protocol()["calibration"]["numeric_epsilon"],
        )

        endpoint = next(value for value in endpoints if value["eligible"])
        self.assertEqual(endpoint["direction"], "decrease")
        self.assertTrue(policy["frozen"])
        self.assertEqual(len(monotonic), 6)
        self.assertTrue(all(value["monotonic"] for value in monotonic))

        reversed_effects = [dict(value) for value in effects]
        for value in reversed_effects:
            if (
                value["configuration_accession"] == "CONFIG0"
                and value["event_fraction"] == 0.1
            ):
                value["signed_effect"] = 0.2
                value["absolute_effect"] = 0.2
        reversed_monotonic = MODULE.monotonic_rows(
            analysis_id="analysis",
            effects=reversed_effects,
            policy=policy,
            epsilon=protocol()["calibration"]["numeric_epsilon"],
        )
        configuration_zero = next(
            value
            for value in reversed_monotonic
            if value["configuration_accession"] == "CONFIG0"
        )
        self.assertFalse(configuration_zero["monotonic"])

    def test_cli_analyzes_reconciled_execution_without_premature_freeze(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = execution_fixture.create_materialization(root)
            execution = execution_fixture.MODULE.run_perturbation_calibration(
                materialization_manifest_path=inputs["materialization"],
                execution_protocol_path=inputs["protocol"],
                seqcheck_bin=inputs["seqcheck"],
                output_root=root / "execution",
                timeout_seconds=10,
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--execution-manifest",
                    execution["manifest_path"],
                    "--analysis-protocol",
                    str(PROTOCOL_PATH),
                    "--output-root",
                    str(root / "analysis"),
                ],
                capture_output=True,
                text=True,
            )
            manifest = MODULE.runtime.load_json(Path(completed.stdout.strip()))
            policy = MODULE.runtime.load_json(
                Path(manifest["outputs"]["policy"]["path"])
            )
            calls = read_csv(manifest["outputs"]["calls"]["path"])
            validation = MODULE.runtime.load_json(
                Path(manifest["outputs"]["validation"]["path"])
            )

            execution_path = Path(execution["manifest_path"])
            execution_payload = MODULE.runtime.load_json(execution_path)
            run_path = Path(execution_payload["outputs"]["runs"]["path"])
            run_rows = read_csv(run_path)
            run_rows[0]["condition_id"] = "not-a-condition"
            MODULE.runtime.write_csv(
                run_path,
                run_rows,
                list(execution_fixture.MODULE.RUN_FIELDS),
            )
            execution_payload["outputs"]["runs"] = MODULE.runtime.file_identity(
                run_path
            )
            MODULE.runtime.write_json(execution_path, execution_payload)
            with self.assertRaisesRegex(ValueError, "condition identifiers"):
                MODULE.load_execution(execution_path)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(manifest["valid"])
        self.assertFalse(policy["frozen"])
        self.assertEqual(len(calls), manifest["counts"]["conditions"])
        self.assertEqual(
            validation["scientific_diagnostics"]["status"], "policy_not_frozen"
        )


if __name__ == "__main__":
    unittest.main()
