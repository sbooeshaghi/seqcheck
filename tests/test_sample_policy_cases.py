import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import test_sample_fastq as fastq_fixture


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "sample_policy_cases.py"
SAMPLER_PATH = ROOT / "scripts" / "sample_fastq.py"
SPEC = importlib.util.spec_from_file_location("sample_policy_cases", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, value: object) -> None:
    MODULE.runtime.write_json(path, value)


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def create_inputs(root: Path, *, method: str = "prefix") -> dict[str, Path]:
    root.mkdir(parents=True)
    source = root / "reads.fastq.gz"
    fastq_fixture.write_fastq(source, 20)
    spec = root / "spec.json"
    write_json(spec, {"seqspec_version": "0.5.0"})
    sampling_protocol = root / "sampling_protocol.json"
    write_json(sampling_protocol, {"schema_version": "0.1.0"})
    protocol_sha256 = MODULE.runtime.file_sha256(sampling_protocol)

    selected_case = {
        "family_id": "rna_family",
        "configuration_accession": "EVAL1",
        "modality": "rna",
        "selection_role": "primary",
        "read_id": "read1",
        "fastq_accession": "FASTQ1",
        "fastq_url": str(source),
        "declared_compressed_bytes": source.stat().st_size,
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
        "cohort_split": "evaluation",
        "cohort_sha256": "cohort-test",
        "sampling_protocol_sha256": protocol_sha256,
        "selector": {"version": "test", "sha256": "selector-test"},
        "seqspec": {"version": "test", "sha256": "seqspec-test"},
        "cases": [stable_case],
    }
    selection_id = MODULE.runtime.sha256_json(stable_selection)[:16]
    case = {
        **selected_case,
        "selection_id": selection_id,
        "case_id": "eval1--fastq1",
    }
    cases_path = root / "cases.json"
    write_json(
        cases_path,
        {
            "schema_version": "0.1.0",
            "selection_id": selection_id,
            "freeze_id": "freeze-test",
            "cohort_split": "evaluation",
            "cases": [case],
        },
    )
    selection_path = root / "selection.json"
    write_json(
        selection_path,
        {
            **stable_selection,
            "selection_id": selection_id,
            "valid": True,
            "counts": {"fastq_cases": 1},
            "outputs": {"cases_json": MODULE.runtime.file_identity(cases_path)},
        },
    )

    sampler_identity = MODULE.runtime.script_identity(
        SAMPLER_PATH, version=MODULE.sampler.SAMPLER_VERSION
    )
    stable_study = {
        "schema_version": "0.1.0",
        "selection_id": "calibration-selection",
        "sampling_protocol_sha256": protocol_sha256,
        "seqcheck": {"version": "test", "sha256": "seqcheck-test"},
        "sampler": {"version": "test", "sha256": "matrix-test", "python": "test"},
        "sampler_core": MODULE.runtime.functional_script_identity(sampler_identity),
        "runtime": {"version": "test", "sha256": "runtime-test", "python": "test"},
        "runner": {"version": "test", "sha256": "runner-test", "python": "test"},
    }
    study_run_id = MODULE.runtime.sha256_json(stable_study)[:16]
    study_path = root / "study.json"
    write_json(
        study_path,
        {
            **stable_study,
            "study_run_id": study_run_id,
            "valid": True,
            "inputs": {
                "sampling_protocol": MODULE.runtime.file_identity(sampling_protocol)
            },
        },
    )

    stable_policy = {
        "schema_version": "0.1.0",
        "study_run_id": study_run_id,
        "analysis_run_id": "analysis-test",
        "analysis_protocol_sha256": "analysis-test",
        "frozen": True,
        "scientific_targets_met": True,
        "default": {
            "records_per_fastq": 5,
            "sampling_method": method,
            "sampling_seed": 17 if method == "reservoir" else None,
        },
        "escalation_records_per_fastq": 10,
        "endpoint_decisions": [],
        "memory": {"evaluable": True, "pass": True},
    }
    policy_path = root / "policy.json"
    write_json(
        policy_path,
        {
            **stable_policy,
            "policy_id": MODULE.runtime.sha256_json(stable_policy)[:16],
        },
    )
    return {
        "selection": selection_path,
        "study": study_path,
        "policy": policy_path,
        "source": source,
    }


class SamplePolicyCasesTests(unittest.TestCase):
    def test_bundle_is_content_addressed_and_samples_every_case_once(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "inputs")
            first = MODULE.sample_policy_cases(
                case_selection_manifest_path=inputs["selection"],
                sampling_study_path=inputs["study"],
                sampling_policy_path=inputs["policy"],
                sampler_script=SAMPLER_PATH,
                output_root=root / "first",
                timeout_seconds=10,
            )
            second = MODULE.sample_policy_cases(
                case_selection_manifest_path=inputs["selection"],
                sampling_study_path=inputs["study"],
                sampling_policy_path=inputs["policy"],
                sampler_script=SAMPLER_PATH,
                output_root=root / "second",
                timeout_seconds=10,
            )
            rows = read_csv(first["outputs"]["samples"]["path"])
            validation = MODULE.runtime.load_json(
                Path(first["outputs"]["validation"]["path"])
            )

        self.assertEqual(first["bundle_id"], second["bundle_id"])
        self.assertTrue(validation["valid"])
        self.assertEqual(first["counts"]["samples"], 1)
        self.assertEqual(first["counts"]["records_selected"], 5)
        self.assertEqual(rows[0]["sampling_method"], "prefix")
        self.assertEqual(rows[0]["records_streamed"], "5")
        self.assertEqual(rows[0]["records_selected"], "5")

    def test_cli_applies_the_frozen_reservoir_seed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "inputs", method="reservoir")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--case-selection-manifest",
                    str(inputs["selection"]),
                    "--sampling-study",
                    str(inputs["study"]),
                    "--sampling-policy",
                    str(inputs["policy"]),
                    "--sampler-script",
                    str(SAMPLER_PATH),
                    "--output-root",
                    str(root / "bundle"),
                    "--timeout-seconds",
                    "10",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            manifest = json.loads(
                Path(completed.stdout.strip()).read_text(encoding="utf-8")
            )
            rows = read_csv(manifest["outputs"]["samples"]["path"])

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(rows[0]["sampling_method"], "reservoir")
        self.assertEqual(rows[0]["sampling_seed"], "17")
        self.assertEqual(rows[0]["records_streamed"], "20")

    def test_runner_rejects_a_non_evaluation_selection(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "inputs")
            selection = MODULE.runtime.load_json(inputs["selection"])
            selection["cohort_split"] = "calibration"
            stable = {
                field: selection[field]
                for field in MODULE.inventory_builder.SELECTION_IDENTITY_FIELDS
            }
            selection["selection_id"] = MODULE.runtime.sha256_json(stable)[:16]
            cases_path = Path(selection["outputs"]["cases_json"]["path"])
            cases = MODULE.runtime.load_json(cases_path)
            cases["selection_id"] = selection["selection_id"]
            cases["cohort_split"] = "calibration"
            cases["cases"][0]["selection_id"] = selection["selection_id"]
            MODULE.runtime.write_json(cases_path, cases)
            selection["outputs"]["cases_json"] = MODULE.runtime.file_identity(
                cases_path
            )
            MODULE.runtime.write_json(inputs["selection"], selection)

            with self.assertRaisesRegex(ValueError, "require an evaluation"):
                MODULE.sample_policy_cases(
                    case_selection_manifest_path=inputs["selection"],
                    sampling_study_path=inputs["study"],
                    sampling_policy_path=inputs["policy"],
                    sampler_script=SAMPLER_PATH,
                    output_root=root / "bundle",
                    timeout_seconds=10,
                )


if __name__ == "__main__":
    unittest.main()
