import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "paper_runtime.py"
SPEC = importlib.util.spec_from_file_location("paper_runtime", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class PaperRuntimeTests(unittest.TestCase):
    def test_flatten_report_assessments_preserves_scope_and_ontology(self) -> None:
        payload = {
            "results": [
                {
                    "check": "fixed",
                    "files": ["R1"],
                    "reads": ["read1"],
                    "regions": ["linker"],
                    "ontology": ["RGN:technical:linker"],
                    "assessment": [
                        {
                            "type": "warning",
                            "code": "fixed_exact_match_absent",
                            "description": "No exact matches.",
                            "expected_ids": ["e1"],
                            "observed_ids": ["o1", "o2"],
                        }
                    ],
                }
            ]
        }

        rows = MODULE.flatten_report_assessments(payload, {"condition_id": "c1"})

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["condition_id"], "c1")
        self.assertEqual(rows[0]["regions"], "linker")
        self.assertEqual(rows[0]["ontology_terms"], "RGN:technical:linker")
        self.assertEqual(rows[0]["assessment_code"], "fixed_exact_match_absent")
        self.assertEqual(rows[0]["observed_metric_ids"], "o1;o2")

    def test_measured_command_captures_outputs_and_resource_usage(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            measurement = MODULE.run_measured_command(
                [
                    sys.executable,
                    "-c",
                    (
                        "import sys,time; data=bytearray(1000000); "
                        "print(len(data)); print('note', file=sys.stderr); "
                        "time.sleep(0.02)"
                    ),
                ],
                stdout_path=root / "stdout.txt",
                stderr_path=root / "stderr.txt",
                timeout_seconds=5,
            )
            MODULE.require_success(measurement)
            stdout = Path(measurement["stdout"]["path"]).read_text().strip()
            stderr = Path(measurement["stderr"]["path"]).read_text().strip()

        self.assertEqual(measurement["exit_code"], 0)
        self.assertFalse(measurement["timed_out"])
        self.assertEqual(stdout, "1000000")
        self.assertEqual(stderr, "note")
        self.assertGreater(measurement["wall_time_seconds"], 0)
        self.assertGreater(measurement["peak_resident_memory_bytes"], 0)

    def test_measured_command_reports_failure_and_timeout(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            failed = MODULE.run_measured_command(
                [sys.executable, "-c", "import sys; sys.exit(7)"],
                stdout_path=root / "failed.out",
                stderr_path=root / "failed.err",
                timeout_seconds=5,
            )
            with self.assertRaisesRegex(ValueError, "exited 7"):
                MODULE.require_success(failed)

            timed_out = MODULE.run_measured_command(
                [sys.executable, "-c", "import time; time.sleep(5)"],
                stdout_path=root / "timeout.out",
                stderr_path=root / "timeout.err",
                timeout_seconds=1,
            )
            with self.assertRaisesRegex(ValueError, "exceeded"):
                MODULE.require_success(timed_out)

        self.assertEqual(failed["exit_code"], 7)
        self.assertTrue(timed_out["timed_out"])

    def test_report_metrics_flatten_and_reconcile(self) -> None:
        payload = {
            "report_schema_version": "0.1.0",
            "results": [
                {
                    "check": "coverage",
                    "files": ["R1"],
                    "reads": ["read1"],
                    "regions": ["barcode"],
                    "expected": [
                        {
                            "id": "e1",
                            "name": "expected_region",
                            "description": "Expected region.",
                            "data": {
                                "kind": "records",
                                "value": [{"region_type": ["RGN:partition:cell"]}],
                            },
                        }
                    ],
                    "observed": [
                        {
                            "id": "o1",
                            "name": "covered_fraction",
                            "description": "Covered fraction.",
                            "data": {
                                "kind": "scalar",
                                "value": 0.75,
                                "unit": "fraction",
                            },
                        }
                    ],
                }
            ],
        }
        rows = MODULE.flatten_report_metrics(payload, {"case_id": "case"})

        self.assertEqual(len(rows), MODULE.count_report_metrics(payload))
        self.assertEqual(rows[0]["case_id"], "case")
        self.assertEqual(rows[0]["ontology_terms"], "RGN:partition:cell")
        self.assertEqual(rows[1]["value_json"], "0.75")
        self.assertEqual(json.loads(rows[1]["value_json"]), 0.75)

    def test_report_flatten_preserves_nominal_region_annotations(self) -> None:
        payload = {
            "results": [
                {
                    "check": "random",
                    "files": ["R1"],
                    "reads": ["read1"],
                    "regions": ["umi"],
                    "expected": [],
                    "observed": [
                        {
                            "id": "o1",
                            "name": "sequence_entropy_fraction",
                            "description": "Entropy fraction.",
                            "data": {"kind": "scalar", "value": 0.9},
                        }
                    ],
                    "ontology": ["RGN:partition:molecule"],
                }
            ]
        }
        spec = {
            "library_spec": [
                {
                    "region_id": "rna",
                    "region_type": ["RGN:unknown:unclassified"],
                    "sequence_type": "joined",
                    "regions": [
                        {
                            "region_id": "umi",
                            "region_type": ["RGN:partition:molecule"],
                            "sequence_type": "random",
                            "regions": [],
                        }
                    ],
                }
            ]
        }

        rows = MODULE.flatten_report_metrics(
            payload,
            {"case_id": "case"},
            region_annotations=MODULE.index_seqspec_regions(spec),
        )

        self.assertEqual(rows[0]["sequence_types"], "random")
        self.assertEqual(
            json.loads(rows[0]["region_annotations_json"]),
            [
                {
                    "ontology_terms": ["RGN:partition:molecule"],
                    "region_id": "umi",
                    "sequence_type": "random",
                }
            ],
        )

    def test_report_flatten_rejects_unknown_region_annotation(self) -> None:
        with self.assertRaisesRegex(ValueError, "unknown region"):
            MODULE.flatten_report_metrics(
                {"results": [{"regions": ["missing"]}]},
                {},
                region_annotations={},
            )

    def test_seqspec_region_index_is_scoped_to_modality(self) -> None:
        spec = {
            "library_spec": [
                {
                    "region_id": "rna",
                    "region_type": ["RGN:unknown:unclassified"],
                    "sequence_type": "joined",
                    "regions": [
                        {
                            "region_id": "barcode",
                            "region_type": ["RGN:partition:cell"],
                            "sequence_type": "onlist",
                            "regions": [],
                        }
                    ],
                },
                {
                    "region_id": "atac",
                    "region_type": ["RGN:unknown:unclassified"],
                    "sequence_type": "joined",
                    "regions": [
                        {
                            "region_id": "barcode",
                            "region_type": ["RGN:partition:nucleus"],
                            "sequence_type": "onlist",
                            "regions": [],
                        }
                    ],
                },
            ]
        }

        indexed = MODULE.index_seqspec_regions(spec, modality="atac")

        self.assertEqual(
            indexed["barcode"]["ontology_terms"], ["RGN:partition:nucleus"]
        )
        with self.assertRaisesRegex(ValueError, "not unique"):
            MODULE.index_seqspec_regions(spec)

    def test_csv_and_script_identity_are_content_addressed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            table = root / "table.csv"
            script = root / "tools" / "script.py"
            executable = root / "tools" / "executable"
            script.parent.mkdir()
            script.write_text("print('first')\n", encoding="utf-8")
            executable.write_text(
                "#!/bin/sh\necho 'paper-tool test'\n", encoding="utf-8"
            )
            executable.chmod(0o755)
            MODULE.write_csv(table, [{"value": 1}], ["value"])

            first = MODULE.script_identity(script, version="test")
            executable_value = MODULE.executable_identity(executable)
            rows = MODULE.read_csv(table)
            script.write_text("print('second')\n", encoding="utf-8")
            second = MODULE.script_identity(script, version="test")

        self.assertEqual(rows, [{"value": "1"}])
        self.assertNotEqual(first["sha256"], second["sha256"])
        self.assertEqual(executable_value["version"], "paper-tool test")
        self.assertEqual(
            MODULE.functional_executable_identity(executable_value),
            {
                "version": "paper-tool test",
                "sha256": executable_value["sha256"],
            },
        )
        self.assertEqual(
            MODULE.functional_script_identity(first),
            {
                "version": "test",
                "sha256": first["sha256"],
                "python": first["python"],
            },
        )


if __name__ == "__main__":
    unittest.main()
