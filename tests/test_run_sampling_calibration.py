import csv
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "run_sampling_calibration.py"
)
SAMPLER_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "sample_fastq_matrix.py"
)
SPEC = importlib.util.spec_from_file_location("run_sampling_calibration", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def write_fastq(path: Path, count: int) -> None:
    records = []
    for index in range(count):
        records.extend(
            [
                f"@read-{index}\n",
                "ACGT\n",
                "+\n",
                "!!!!\n",
            ]
        )
    path.write_text("".join(records), encoding="ascii")


def fake_seqcheck(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import gzip
import json
import sys
from pathlib import Path

if sys.argv[1] == "--version":
    print("seqcheck 0.2.0-test")
    raise SystemExit(0)
output = Path(sys.argv[sys.argv.index("--output") + 1])
source = Path(sys.argv[-1])
opener = gzip.open if source.suffix == ".gz" else open
with opener(source, "rt", encoding="ascii") as handle:
    count = sum(1 for _ in handle) // 4
payload = {
    "report_schema_version": "0.1.0",
    "meta": {"command": "check", "requested_reads": 0},
    "results": [
        {
            "check": "length",
            "files": ["FASTQ1"],
            "reads": ["read1"],
            "regions": ["insert"],
            "expected": [
                {
                    "id": "e1",
                    "name": "expected_fraction",
                    "description": "Expected fraction.",
                    "data": {
                        "kind": "scalar",
                        "value": 1.0,
                        "unit": "fraction",
                    },
                }
            ],
            "observed": [
                {
                    "id": "o1",
                    "name": "sampled_count",
                    "description": "Sampled count.",
                    "data": {
                        "kind": "scalar",
                        "value": count,
                        "unit": "count",
                    },
                },
                {
                    "id": "o2",
                    "name": "covered_fraction",
                    "description": "Covered fraction.",
                    "data": {
                        "kind": "scalar",
                        "value": count / 12,
                        "unit": "fraction",
                    },
                },
            ],
            "assessment": [],
            "ontology": ["RGN:measure:transcript"],
        }
    ],
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def create_inputs(root: Path) -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    seqcheck = root / "seqcheck"
    fake_seqcheck(seqcheck)
    fastq = root / "FASTQ1.fastq"
    write_fastq(fastq, 12)
    spec = root / "spec.yaml"
    spec.write_text("seqspec_version: 0.5.0\n", encoding="utf-8")
    protocol = root / "protocol.json"
    write_json(
        protocol,
        {
            "schema_version": "0.1.0",
            "sampling": {
                "records_per_fastq": [3],
                "include_prefix": True,
                "reservoir_seeds": [7],
                "source_traversals_for_bounded_samples_per_fastq": 1,
                "synchronize_mates_when_paired": True,
                "complete_stream_reference": True,
            },
            "performance": {"warmup_runs": 1, "measured_runs": 2},
            "expected_conditions_per_fastq": {
                "prefix": 1,
                "reservoir": 1,
                "complete_stream": 1,
                "total": 3,
            },
        },
    )
    selected_case = {
        "family_id": "rna_family",
        "configuration_accession": "CONFIG1",
        "modality": "rna",
        "selection_role": "primary",
        "read_id": "read1",
        "fastq_accession": "FASTQ1",
        "fastq_url": str(fastq),
        "declared_compressed_bytes": fastq.stat().st_size,
        "indexed_bases": 4,
        "measure_bases": 4,
        "partition_bases": 0,
        "technical_bases": 0,
        "role_classes": "measure",
        "ontology_terms": "RGN:measure:transcript",
        "spec_path": str(spec),
        "spec_sha256": MODULE.runtime.file_sha256(spec),
    }
    stable_case = {
        key: value for key, value in selected_case.items() if key != "spec_path"
    }
    stable_selection = {
        "schema_version": "0.1.0",
        "freeze_id": "freeze-test",
        "cohort_sha256": "cohort-sha256-test",
        "sampling_protocol_sha256": MODULE.runtime.file_sha256(protocol),
        "selector": {"version": "0.1.0", "sha256": "selector-sha256-test"},
        "seqspec": {"version": "seqspec 0.5.0-test", "sha256": "seqspec-test"},
        "cases": [stable_case],
    }
    selection_id = MODULE.runtime.sha256_json(stable_selection)[:16]
    cases_path = root / "cases.json"
    write_json(
        cases_path,
        {
            "schema_version": "0.1.0",
            "selection_id": selection_id,
            "freeze_id": "freeze-test",
            "cases": [
                {
                    **selected_case,
                    "selection_id": selection_id,
                    "case_id": "config1--fastq1",
                }
            ],
        },
    )
    selection_manifest = root / "selection.json"
    write_json(
        selection_manifest,
        {
            **stable_selection,
            "selection_id": selection_id,
            "valid": True,
            "counts": {"fastq_cases": 1},
            "outputs": {"cases_json": MODULE.runtime.file_identity(cases_path)},
        },
    )
    return {
        "seqcheck": seqcheck,
        "fastq": fastq,
        "spec": spec,
        "protocol": protocol,
        "selection": selection_manifest,
    }


def read_csv(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class RunSamplingCalibrationTests(unittest.TestCase):
    def test_runner_executes_and_reconciles_all_conditions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root)
            manifest = MODULE.run_sampling_calibration(
                case_selection_manifest_path=inputs["selection"],
                protocol_path=inputs["protocol"],
                seqcheck_bin=inputs["seqcheck"],
                sampler_script=SAMPLER_PATH,
                output_root=root / "output",
                timeout_seconds=10,
            )
            validation = json.loads(
                Path(manifest["outputs"]["validation"]["path"]).read_text()
            )
            reports = read_csv(manifest["outputs"]["reports"]["path"])
            metrics = read_csv(manifest["outputs"]["metrics"]["path"])
            errors = read_csv(manifest["outputs"]["sampling_errors"]["path"])
            performance = read_csv(manifest["outputs"]["performance"]["path"])

        self.assertTrue(manifest["valid"])
        self.assertTrue(validation["valid"])
        self.assertEqual(validation["counts"]["reports"], 3)
        self.assertEqual(validation["counts"]["metrics"], 9)
        self.assertEqual(validation["counts"]["sampling_errors"], 4)
        self.assertEqual(validation["counts"]["performance_rows"], 8)
        self.assertEqual(validation["counts"]["repeat_reports"], 4)
        self.assertEqual(
            {row["sampling_method"] for row in reports},
            {"prefix", "reservoir", "complete"},
        )
        self.assertEqual(len(metrics), 9)
        self.assertEqual(len(errors), 4)
        self.assertEqual(len(performance), 8)
        self.assertTrue(all(json.loads(row["command_json"])[0] for row in performance))
        self.assertTrue(all(len(row["stdout_sha256"]) == 64 for row in performance))
        self.assertTrue(all(len(row["stderr_sha256"]) == 64 for row in performance))
        self.assertEqual(manifest["tools"]["runner"]["version"], MODULE.RUNNER_VERSION)
        self.assertEqual(
            manifest["tools"]["sampler"]["version"], MODULE.SAMPLER_VERSION
        )
        self.assertEqual(
            manifest["tools"]["sampler"]["sha256"],
            MODULE.runtime.file_sha256(SAMPLER_PATH),
        )
        self.assertEqual(
            manifest["tools"]["sampler_core"]["sha256"],
            MODULE.runtime.file_sha256(SAMPLER_PATH.with_name("sample_fastq.py")),
        )
        self.assertEqual(
            manifest["tools"]["runtime"]["sha256"],
            MODULE.runtime.file_sha256(Path(MODULE.runtime.__file__)),
        )
        complete = next(row for row in reports if row["sampling_method"] == "complete")
        self.assertEqual(complete["records_processed"], "12")
        covered = next(
            row for row in errors if row["metric_name"] == "covered_fraction"
        )
        self.assertEqual(float(covered["complete_value"]), 1.0)
        self.assertEqual(float(covered["sample_value"]), 0.25)
        self.assertEqual(float(covered["signed_error"]), -0.75)
        self.assertEqual(float(covered["absolute_error"]), 0.75)
        self.assertEqual(covered["metric_id"], "o2")

    def test_runner_rejects_changed_protocol_and_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root)
            output_root = root / "output"
            output_root.mkdir()
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                MODULE.run_sampling_calibration(
                    case_selection_manifest_path=inputs["selection"],
                    protocol_path=inputs["protocol"],
                    seqcheck_bin=inputs["seqcheck"],
                    sampler_script=SAMPLER_PATH,
                    output_root=output_root,
                    timeout_seconds=10,
                )

            protocol = json.loads(inputs["protocol"].read_text())
            protocol["sampling"]["reservoir_seeds"] = [8]
            write_json(inputs["protocol"], protocol)
            with self.assertRaisesRegex(ValueError, "hash differs"):
                MODULE.run_sampling_calibration(
                    case_selection_manifest_path=inputs["selection"],
                    protocol_path=inputs["protocol"],
                    seqcheck_bin=inputs["seqcheck"],
                    sampler_script=SAMPLER_PATH,
                    output_root=root / "changed",
                    timeout_seconds=10,
                )

    def test_runner_rejects_tampered_content_addressed_selection(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root)
            selection = json.loads(inputs["selection"].read_text())
            selection["cases"][0]["fastq_accession"] = "CHANGED"
            write_json(inputs["selection"], selection)

            with self.assertRaisesRegex(ValueError, "not content-addressed"):
                MODULE.run_sampling_calibration(
                    case_selection_manifest_path=inputs["selection"],
                    protocol_path=inputs["protocol"],
                    seqcheck_bin=inputs["seqcheck"],
                    sampler_script=SAMPLER_PATH,
                    output_root=root / "changed-selection",
                    timeout_seconds=10,
                )

            inputs = create_inputs(root / "case-tamper")
            selection = json.loads(inputs["selection"].read_text())
            cases_path = Path(selection["outputs"]["cases_json"]["path"])
            cases = json.loads(cases_path.read_text())
            cases["cases"][0]["fastq_accession"] = "CHANGED"
            write_json(cases_path, cases)
            selection["outputs"]["cases_json"] = MODULE.runtime.file_identity(
                cases_path
            )
            write_json(inputs["selection"], selection)

            with self.assertRaisesRegex(ValueError, "rows differ"):
                MODULE.run_sampling_calibration(
                    case_selection_manifest_path=inputs["selection"],
                    protocol_path=inputs["protocol"],
                    seqcheck_bin=inputs["seqcheck"],
                    sampler_script=SAMPLER_PATH,
                    output_root=root / "changed-cases",
                    timeout_seconds=10,
                )

    def test_sampler_dependency_changes_study_identity(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root)
            first_tools = root / "first-tools"
            second_tools = root / "second-tools"
            first_tools.mkdir()
            second_tools.mkdir()
            for destination in (first_tools, second_tools):
                shutil.copy2(SAMPLER_PATH, destination / SAMPLER_PATH.name)
                shutil.copy2(
                    SAMPLER_PATH.with_name("sample_fastq.py"),
                    destination / "sample_fastq.py",
                )
            with (second_tools / "sample_fastq.py").open(
                "a", encoding="utf-8"
            ) as handle:
                handle.write("\n# Distinct sampling implementation identity.\n")

            first = MODULE.run_sampling_calibration(
                case_selection_manifest_path=inputs["selection"],
                protocol_path=inputs["protocol"],
                seqcheck_bin=inputs["seqcheck"],
                sampler_script=first_tools / SAMPLER_PATH.name,
                output_root=root / "first-output",
                timeout_seconds=10,
            )
            second = MODULE.run_sampling_calibration(
                case_selection_manifest_path=inputs["selection"],
                protocol_path=inputs["protocol"],
                seqcheck_bin=inputs["seqcheck"],
                sampler_script=second_tools / SAMPLER_PATH.name,
                output_root=root / "second-output",
                timeout_seconds=10,
            )

        self.assertNotEqual(first["study_run_id"], second["study_run_id"])
        self.assertNotEqual(
            first["sampler_core"]["sha256"], second["sampler_core"]["sha256"]
        )

    def test_cli_writes_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--case-selection-manifest",
                    str(inputs["selection"]),
                    "--sampling-protocol",
                    str(inputs["protocol"]),
                    "--seqcheck-bin",
                    str(inputs["seqcheck"]),
                    "--sampler-script",
                    str(SAMPLER_PATH),
                    "--output-root",
                    str(root / "output"),
                    "--timeout-seconds",
                    "10",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest = json.loads(Path(completed.stdout.strip()).read_text())

        self.assertTrue(manifest["valid"])
        self.assertEqual(manifest["counts"]["reports"], 3)


if __name__ == "__main__":
    unittest.main()
