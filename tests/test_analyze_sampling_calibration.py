import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "analyze_sampling_calibration.py"
)
SPEC = importlib.util.spec_from_file_location(
    "analyze_sampling_calibration", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

ERROR_FIELDS = (
    "study_run_id",
    "selection_id",
    "case_id",
    "family_id",
    "configuration_accession",
    "fastq_accession",
    "selection_role",
    "condition_id",
    "sampling_method",
    "requested_records_per_fastq",
    "sampling_seed",
    "check",
    "files",
    "reads",
    "regions",
    "ontology_terms",
    "metric_id",
    "metric_name",
    "unit",
    "complete_value",
    "sample_value",
    "signed_error",
    "absolute_error",
)
PERFORMANCE_FIELDS = (
    "study_run_id",
    "case_id",
    "sampling_method",
    "requested_records_per_fastq",
    "execution_scope",
    "is_warmup",
    "wall_time_seconds",
    "user_cpu_seconds",
    "system_cpu_seconds",
    "peak_resident_memory_bytes",
    "records_processed",
    "report_size_bytes",
)
METRIC_FIELDS = (
    "study_run_id",
    "case_id",
    "family_id",
    "configuration_accession",
    "fastq_accession",
    "selection_role",
    "condition_id",
    "sampling_method",
    "requested_records_per_fastq",
    "sampling_seed",
    "check",
    "files",
    "reads",
    "regions",
    "ontology_terms",
    "metric_side",
    "metric_id",
    "metric_name",
    "data_kind",
    "unit",
    "value_json",
)


def write_json(path: Path, value: object) -> None:
    MODULE.runtime.write_json(path, value)


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def analysis_protocol() -> dict[str, object]:
    return {
        "schema_version": "0.1.0",
        "analysis_id": "test-analysis",
        "accuracy_aggregation": "configuration_max_of_fastq_case_max_absolute_error",
        "prefix_bias_aggregation": "configuration_mean_of_fastq_case_means",
        "policy": {
            "preferred_records_per_fastq": 10000,
            "escalation_records_per_fastq": 100000,
            "reservoir_seed": 17,
            "median_absolute_error_max": 0.01,
            "p95_absolute_error_max": 0.05,
            "prefix_bias_ci_lower_min": -0.01,
            "prefix_bias_ci_upper_max": 0.01,
            "peak_resident_memory_bytes_max": 1073741824,
        },
        "statistics": {
            "minimum_independent_configurations_per_endpoint": 3,
            "confidence_level": 0.95,
            "bootstrap_resamples": 1000,
            "bootstrap_seed": 1234,
        },
        "endpoint_metrics": [
            {
                "endpoint_id": "stable_fraction",
                "check": "coverage",
                "metric_name": "covered_fraction",
                "unit": "fraction",
            },
            {
                "endpoint_id": "escalating_fraction",
                "check": "fixed",
                "metric_name": "exact_match_fraction",
                "unit": "fraction",
            },
        ],
        "descriptive_metrics": [
            {
                "metric_id": "unique_count",
                "check": "random",
                "metric_name": "unique_sequence_count",
                "unit": "count",
                "comparison": "prefix_vs_reservoir_at_equal_sample_size",
            }
        ],
        "event_prevalences": [0.001, 0.01],
    }


def row_context(case_index: int) -> dict[str, object]:
    return {
        "case_id": f"case-{case_index}",
        "family_id": "rna_family",
        "configuration_accession": f"CONFIG{case_index}",
        "fastq_accession": f"FASTQ{case_index}",
        "selection_role": "primary",
        "files": f"FASTQ{case_index}",
        "reads": "read1",
        "regions": "insert",
        "ontology_terms": "RGN:measure:transcript",
    }


def error_row(
    *,
    study_run_id: str,
    case_index: int,
    endpoint: str,
    n_reads: int,
    method: str,
    seed: int | None,
    error: float,
) -> dict[str, object]:
    endpoint_values = {
        "stable": ("coverage", "covered_fraction", "o1", 0.8),
        "escalating": ("fixed", "exact_match_fraction", "o2", 0.9),
    }
    check, metric_name, metric_id, complete = endpoint_values[endpoint]
    return {
        "study_run_id": study_run_id,
        "selection_id": "selection-test",
        **row_context(case_index),
        "condition_id": f"{method}-{n_reads}-{seed}",
        "sampling_method": method,
        "requested_records_per_fastq": n_reads,
        "sampling_seed": seed,
        "check": check,
        "metric_id": metric_id,
        "metric_name": metric_name,
        "unit": "fraction",
        "complete_value": complete,
        "sample_value": complete + error,
        "signed_error": error,
        "absolute_error": abs(error),
    }


def create_study(root: Path) -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    study_run_id = "study-test"
    sampling_protocol = root / "sampling_protocol.json"
    write_json(
        sampling_protocol,
        {
            "schema_version": "0.1.0",
            "sampling": {
                "records_per_fastq": [10000, 100000],
                "reservoir_seeds": [17, 29, 43],
            },
        },
    )
    analysis_protocol_path = root / "analysis_protocol.json"
    write_json(analysis_protocol_path, analysis_protocol())

    errors = []
    performance = []
    metrics = []
    reservoir_offsets = {17: -0.001, 29: 0.0, 43: 0.001}
    for case_index in range(3):
        for n_reads in (10000, 100000):
            for endpoint in ("stable", "escalating"):
                base_error = (
                    0.005
                    if endpoint == "stable"
                    else (0.08 if n_reads == 10000 else 0.005)
                )
                errors.append(
                    error_row(
                        study_run_id=study_run_id,
                        case_index=case_index,
                        endpoint=endpoint,
                        n_reads=n_reads,
                        method="prefix",
                        seed=None,
                        error=base_error,
                    )
                )
                for seed, offset in reservoir_offsets.items():
                    errors.append(
                        error_row(
                            study_run_id=study_run_id,
                            case_index=case_index,
                            endpoint=endpoint,
                            n_reads=n_reads,
                            method="reservoir",
                            seed=seed,
                            error=base_error + offset,
                        )
                    )

            for method, seeds in (("prefix", [None]), ("reservoir", [17, 29, 43])):
                for seed in seeds:
                    performance.append(
                        {
                            "study_run_id": study_run_id,
                            "case_id": f"case-{case_index}",
                            "sampling_method": method,
                            "requested_records_per_fastq": n_reads,
                            "execution_scope": "local_compute",
                            "is_warmup": False,
                            "wall_time_seconds": n_reads / 100000,
                            "user_cpu_seconds": n_reads / 200000,
                            "system_cpu_seconds": 0.01,
                            "peak_resident_memory_bytes": 100000000 + n_reads,
                            "records_processed": n_reads,
                            "report_size_bytes": 1000,
                        }
                    )
                    metrics.append(
                        {
                            "study_run_id": study_run_id,
                            **row_context(case_index),
                            "condition_id": f"{method}-{n_reads}-{seed}",
                            "sampling_method": method,
                            "requested_records_per_fastq": n_reads,
                            "sampling_seed": seed,
                            "check": "random",
                            "metric_side": "observed",
                            "metric_id": "o5",
                            "metric_name": "unique_sequence_count",
                            "data_kind": "scalar",
                            "unit": "count",
                            "value_json": json.dumps(
                                n_reads // 2 + (0 if seed is None else seed)
                            ),
                        }
                    )

    tables = root / "source_tables"
    errors_path = tables / "sampling_errors.csv"
    performance_path = tables / "performance.csv"
    metrics_path = tables / "metrics.csv"
    MODULE.runtime.write_csv(errors_path, errors, list(ERROR_FIELDS))
    MODULE.runtime.write_csv(performance_path, performance, list(PERFORMANCE_FIELDS))
    MODULE.runtime.write_csv(metrics_path, metrics, list(METRIC_FIELDS))
    validation_path = root / "source_validation.json"
    write_json(
        validation_path,
        {"schema_version": "0.1.0", "study_run_id": study_run_id, "valid": True},
    )
    study_manifest = root / "study.json"
    write_json(
        study_manifest,
        {
            "schema_version": "0.1.0",
            "study_run_id": study_run_id,
            "sampling_protocol_sha256": MODULE.runtime.file_sha256(sampling_protocol),
            "valid": True,
            "inputs": {
                "sampling_protocol": MODULE.runtime.file_identity(sampling_protocol)
            },
            "outputs": {
                "sampling_errors": MODULE.runtime.file_identity(errors_path),
                "performance": MODULE.runtime.file_identity(performance_path),
                "metrics": MODULE.runtime.file_identity(metrics_path),
                "validation": MODULE.runtime.file_identity(validation_path),
            },
        },
    )
    return {
        "study": study_manifest,
        "analysis_protocol": analysis_protocol_path,
        "errors": errors_path,
    }


class AnalyzeSamplingCalibrationTests(unittest.TestCase):
    def test_summaries_cluster_multiple_fastqs_by_configuration(self) -> None:
        protocol = analysis_protocol()
        accuracy_rows = []
        bias_rows = []
        for case_id, configuration in (
            ("case-a1", "CONFIG-A"),
            ("case-a2", "CONFIG-A"),
            ("case-b", "CONFIG-B"),
            ("case-c", "CONFIG-C"),
        ):
            context = {
                "analysis_run_id": "analysis",
                "study_run_id": "study",
                "endpoint_id": "stable_fraction",
                "check": "coverage",
                "metric_name": "covered_fraction",
                "unit": "fraction",
                "case_id": case_id,
                "family_id": "rna_family",
                "configuration_accession": configuration,
                "fastq_accession": case_id,
                "selection_role": "primary",
                "requested_records_per_fastq": 10000,
                "metric_instances": 1,
            }
            accuracy_rows.append(
                {
                    **context,
                    "sampling_method": "prefix",
                    "case_mean_signed_error": 0.001,
                    "case_max_absolute_error": 0.001,
                }
            )
            bias_rows.append(
                {
                    **context,
                    "case_mean_prefix_minus_reservoir": 0.0,
                }
            )

        accuracy = MODULE.summarize_accuracy(accuracy_rows, protocol, "analysis")
        bias = MODULE.summarize_bias(bias_rows, protocol, "analysis")

        self.assertEqual(accuracy[0]["independent_configurations"], 3)
        self.assertEqual(accuracy[0]["fastqs"], 4)
        self.assertEqual(bias[0]["independent_configurations"], 3)
        self.assertEqual(bias[0]["fastqs"], 4)

    def test_analysis_freezes_default_and_escalation_decisions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_study(root)
            manifest = MODULE.analyze_sampling_calibration(
                study_manifest_path=inputs["study"],
                analysis_protocol_path=inputs["analysis_protocol"],
                output_root=root / "analysis",
            )
            validation = json.loads(
                Path(manifest["outputs"]["validation"]["path"]).read_text()
            )
            policy = json.loads(Path(manifest["outputs"]["policy"]["path"]).read_text())
            decisions = read_csv(manifest["outputs"]["endpoint_decisions"]["path"])

        self.assertTrue(manifest["valid"])
        self.assertTrue(manifest["policy_frozen"])
        self.assertTrue(validation["valid"])
        self.assertTrue(policy["frozen"])
        self.assertTrue(policy["scientific_targets_met"])
        self.assertIsNone(policy["default"]["sampling_seed"])
        self.assertEqual(validation["counts"]["selected_error_rows"], 48)
        self.assertEqual(validation["counts"]["case_accuracy_rows"], 24)
        self.assertEqual(validation["counts"]["case_bias_rows"], 12)
        self.assertEqual(validation["counts"]["descriptive_pair_rows"], 6)
        self.assertEqual(validation["counts"]["event_probability_rows"], 4)
        self.assertEqual(
            {row["endpoint_id"]: row["decision"] for row in decisions},
            {
                "stable_fraction": "default_prefix",
                "escalating_fraction": "escalation_prefix",
            },
        )

    def test_reservoir_policy_records_the_predeclared_seed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            protocol_path = Path(tmpdir) / "analysis.json"
            write_json(protocol_path, analysis_protocol())
            policy = MODULE.build_policy(
                study_run_id="study",
                analysis_run_id="analysis",
                analysis_protocol_path=protocol_path,
                analysis_protocol=analysis_protocol(),
                decisions=[
                    {
                        "evidence_complete": True,
                        "decision": "default_reservoir",
                    }
                ],
                memory={"evaluable": True, "pass": True},
            )

        self.assertEqual(
            policy["default"],
            {
                "records_per_fastq": 10000,
                "sampling_method": "reservoir",
                "sampling_seed": 17,
            },
        )

    def test_analysis_is_deterministic_and_rejects_tampering(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_study(root)
            first = MODULE.analyze_sampling_calibration(
                study_manifest_path=inputs["study"],
                analysis_protocol_path=inputs["analysis_protocol"],
                output_root=root / "first",
            )
            second = MODULE.analyze_sampling_calibration(
                study_manifest_path=inputs["study"],
                analysis_protocol_path=inputs["analysis_protocol"],
                output_root=root / "second",
            )
            first_policy = json.loads(
                Path(first["outputs"]["policy"]["path"]).read_text()
            )
            second_policy = json.loads(
                Path(second["outputs"]["policy"]["path"]).read_text()
            )
            first_summary = first["outputs"]["prefix_bias_summary"]["sha256"]
            second_summary = second["outputs"]["prefix_bias_summary"]["sha256"]

            inputs["errors"].write_text(
                inputs["errors"].read_text() + "tampered\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "hash changed"):
                MODULE.analyze_sampling_calibration(
                    study_manifest_path=inputs["study"],
                    analysis_protocol_path=inputs["analysis_protocol"],
                    output_root=root / "tampered",
                )

        self.assertEqual(first["analysis_run_id"], second["analysis_run_id"])
        self.assertEqual(first_policy["policy_id"], second_policy["policy_id"])
        self.assertEqual(first_summary, second_summary)

    def test_cli_writes_valid_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_study(root)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--study-manifest",
                    str(inputs["study"]),
                    "--analysis-protocol",
                    str(inputs["analysis_protocol"]),
                    "--output-root",
                    str(root / "analysis"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest = json.loads(Path(completed.stdout.strip()).read_text())

        self.assertTrue(manifest["valid"])
        self.assertTrue(manifest["policy_frozen"])


if __name__ == "__main__":
    unittest.main()
