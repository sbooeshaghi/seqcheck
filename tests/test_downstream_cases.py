import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER_PATH = ROOT / "scripts" / "build_downstream_cases.py"
RUNNER_PATH = ROOT / "scripts" / "run_downstream_cases.py"
PROTOCOL_PATH = ROOT / "experiments" / "paper" / "protocol" / "downstream_cases.json"
SPEC = importlib.util.spec_from_file_location(
    "build_downstream_cases_test", BUILDER_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
RUNNER_SPEC = importlib.util.spec_from_file_location(
    "run_downstream_cases_test", RUNNER_PATH
)
RUNNER = importlib.util.module_from_spec(RUNNER_SPEC)
assert RUNNER_SPEC.loader is not None
sys.modules[RUNNER_SPEC.name] = RUNNER
RUNNER_SPEC.loader.exec_module(RUNNER)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def write_review_analysis(
    root: Path, count: int = 4
) -> tuple[Path, list[dict[str, object]], str]:
    outcomes = []
    for index in range(count):
        outcomes.append(
            {
                "case_id": f"CASE{index + 1}",
                "source_kind": "candidate",
                "confirmed_problem": True,
                "final_classification": (
                    "confirmed_specification_problem"
                    if index % 2 == 0
                    else "confirmed_read_or_resource_problem"
                ),
                "configuration_accession": f"C{index + 1}",
                "modality": "rna" if index % 2 == 0 else "atac",
                "matched_assay_family": "rna" if index % 2 == 0 else "atac",
                "access_class": "public" if index < 3 else "controlled",
            }
        )
    outcomes_path = root / "review" / "tables" / "outcomes.csv"
    MODULE.runtime.write_csv(outcomes_path, outcomes, list(outcomes[0]))
    validation_stable = {
        "schema_version": "0.1.0",
        "merge_id": "merge",
        "selection_id": "review-selection",
        "tool": {"version": "0.1.0", "sha256": "1" * 64, "python": "3.13"},
        "inputs": {
            "review_manifest": "2" * 64,
            "adjudicated_reviews": "3" * 64,
            "protocol": "4" * 64,
        },
    }
    analysis_id = MODULE.runtime.sha256_json(validation_stable)[:16]
    validation_path = root / "review" / "validation" / "audit_reviews.json"
    write_json(
        validation_path,
        {**validation_stable, "analysis_id": analysis_id, "valid": True},
    )
    manifest_path = root / "review" / "manifests" / "analysis.json"
    write_json(
        manifest_path,
        {
            "schema_version": "0.1.0",
            "analysis_id": analysis_id,
            "valid": True,
            "outputs": {
                "outcomes": MODULE.runtime.file_identity(outcomes_path),
                "validation": MODULE.runtime.file_identity(validation_path),
            },
        },
    )
    return manifest_path, outcomes, analysis_id


def pipeline_source() -> str:
    return """#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--reads", required=True)
parser.add_argument("--spec", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()
value = float(Path(args.spec).read_text().strip())
output = Path(args.output)
output.mkdir(parents=True, exist_ok=True)
(output / "endpoint.json").write_text(json.dumps({
    "endpoint_id": "usable_fraction",
    "value": value,
    "unit": "fraction"
}) + "\\n")
"""


def write_case_registry(
    root: Path, outcomes: list[dict[str, object]], analysis_id: str
) -> tuple[Path, dict[str, object]]:
    script = root / "inputs" / "pipeline.py"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(pipeline_source(), encoding="utf-8")
    cases = []
    for index, outcome in enumerate(outcomes):
        case_id = str(outcome["case_id"])
        case_root = root / "inputs" / case_id
        case_root.mkdir(parents=True)
        reads = case_root / "reads.fastq"
        before_spec = case_root / "before.spec"
        after_spec = case_root / "after.spec"
        reads.write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
        before_spec.write_text(f"{0.4 + index / 100}\n", encoding="utf-8")
        after_spec.write_text(f"{0.6 + index / 100}\n", encoding="utf-8")

        def condition(spec_path: Path) -> dict[str, object]:
            return {
                "inputs": [
                    {"role": "reads", "path": str(reads)},
                    {"role": "seqspec", "path": str(spec_path)},
                    {"role": "pipeline_script", "path": str(script)},
                ],
                "command": [
                    sys.executable,
                    str(script),
                    "--reads",
                    str(reads),
                    "--spec",
                    str(spec_path),
                    "--output",
                    "{output_dir}",
                ],
            }

        cases.append(
            {
                "case_id": case_id,
                "access_available": True,
                "pipeline_available": True,
                "presentation_ready": True,
                "pipeline": {"id": f"pipeline-{index % 2}", "version": "1.0"},
                "endpoint": {
                    "id": "usable_fraction",
                    "label": "Usable fraction",
                    "unit": "fraction",
                    "direction": "increase",
                    "minimum_meaningful_change": 0.1,
                },
                "shared_input_roles": ["reads", "pipeline_script"],
                "before": condition(before_spec),
                "after": condition(after_spec),
            }
        )
    registry = {
        "schema_version": "0.1.0",
        "review_analysis_id": analysis_id,
        "cases": cases,
    }
    path = root / "case_registry.json"
    write_json(path, registry)
    return path, registry


def write_protocol(root: Path, minimum: int = 3, target: int = 3) -> Path:
    value = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    value["selection"]["minimum_case_count"] = minimum
    value["selection"]["target_case_count"] = target
    path = root / "protocol.json"
    write_json(path, value)
    return path


def fixture(root: Path) -> dict[str, Path]:
    analysis_path, outcomes, analysis_id = write_review_analysis(root)
    registry_path, _ = write_case_registry(root, outcomes, analysis_id)
    protocol_path = write_protocol(root)
    return {
        "review_analysis_path": analysis_path,
        "registry_path": registry_path,
        "protocol_path": protocol_path,
    }


def build_selection(root: Path) -> tuple[dict[str, Path], Path]:
    inputs = fixture(root)
    path = MODULE.build_cases(output_root=root / "selection", **inputs)
    return inputs, path


class BuildDownstreamCasesTests(unittest.TestCase):
    def test_builder_freezes_deterministic_diverse_cases(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = fixture(root)
            first_path = MODULE.build_cases(output_root=root / "first", **inputs)
            second_path = MODULE.build_cases(output_root=root / "second", **inputs)
            first = json.loads(first_path.read_text())
            second = json.loads(second_path.read_text())
            selected = json.loads(
                (root / "first" / "cases" / "selected_cases.json").read_text()
            )["cases"]

        self.assertEqual(first["selection_id"], second["selection_id"])
        self.assertEqual(first["counts"]["selected_cases"], 3)
        self.assertEqual(
            {row["final_classification"] for row in selected},
            {
                "confirmed_specification_problem",
                "confirmed_read_or_resource_problem",
            },
        )
        self.assertTrue(
            all(row["changed_input_roles"] == ["seqspec"] for row in selected)
        )
        self.assertTrue(
            all(
                MODULE.normalized_condition_command(row["conditions"]["before"])
                == MODULE.normalized_condition_command(row["conditions"]["after"])
                for row in selected
            )
        )

    def test_builder_rejects_command_confound_and_changed_shared_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = fixture(root)
            registry = json.loads(inputs["registry_path"].read_text())
            registry["cases"][0]["after"]["command"].extend(["--mode", "different"])
            write_json(inputs["registry_path"], registry)
            with self.assertRaisesRegex(ValueError, "commands differ"):
                MODULE.build_cases(output_root=root / "command-confound", **inputs)

            registry["cases"][0]["after"]["command"] = registry["cases"][0]["after"][
                "command"
            ][:-2]
            after_reads = root / "inputs" / "CASE1" / "after-reads.fastq"
            after_reads.write_text("@r\nTGCA\n+\nIIII\n", encoding="utf-8")
            for item in registry["cases"][0]["after"]["inputs"]:
                if item["role"] == "reads":
                    item["path"] = str(after_reads)
            registry["cases"][0]["after"]["command"] = [
                str(after_reads) if value.endswith("reads.fastq") else value
                for value in registry["cases"][0]["after"]["command"]
            ]
            write_json(inputs["registry_path"], registry)
            with self.assertRaisesRegex(ValueError, "shared input changed"):
                MODULE.build_cases(output_root=root / "shared-change", **inputs)

    def test_builder_rejects_unconfirmed_case_and_insufficient_pool(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = fixture(root)
            outcomes_identity = json.loads(inputs["review_analysis_path"].read_text())[
                "outputs"
            ]["outcomes"]
            outcomes_path = Path(outcomes_identity["path"])
            outcomes = MODULE.runtime.read_csv(outcomes_path)
            outcomes[0]["confirmed_problem"] = "False"
            MODULE.runtime.write_csv(outcomes_path, outcomes, list(outcomes[0]))
            analysis = json.loads(inputs["review_analysis_path"].read_text())
            analysis["outputs"]["outcomes"] = MODULE.runtime.file_identity(
                outcomes_path
            )
            write_json(inputs["review_analysis_path"], analysis)
            with self.assertRaisesRegex(ValueError, "not a confirmed candidate"):
                MODULE.build_cases(output_root=root / "unconfirmed", **inputs)

    def test_builder_records_but_excludes_unready_case(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = fixture(root)
            registry = json.loads(inputs["registry_path"].read_text())
            excluded_id = registry["cases"][0]["case_id"]
            registry["cases"][0]["access_available"] = False
            registry["cases"][0]["before"]["inputs"][0]["path"] = "/missing/reads"
            registry["cases"][0]["after"]["inputs"][0]["path"] = "/missing/reads"
            registry["cases"][0]["before"]["command"][0] = "/missing/pipeline"
            registry["cases"][0]["after"]["command"][0] = "/missing/pipeline"
            write_json(inputs["registry_path"], registry)
            manifest_path = MODULE.build_cases(output_root=root / "selection", **inputs)
            manifest = json.loads(manifest_path.read_text())

        self.assertEqual(manifest["counts"]["registered_cases"], 4)
        self.assertEqual(manifest["counts"]["eligible_cases"], 3)
        self.assertNotIn(
            excluded_id, {row["case_id"] for row in manifest["selected_cases"]}
        )

    def test_builder_cli_writes_valid_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = fixture(root)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(BUILDER_PATH),
                    "--review-analysis-manifest",
                    str(inputs["review_analysis_path"]),
                    "--case-registry",
                    str(inputs["registry_path"]),
                    "--protocol",
                    str(inputs["protocol_path"]),
                    "--output-root",
                    str(root / "selection"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            manifest = json.loads(Path(completed.stdout.strip()).read_text())

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(manifest["valid"])


class RunDownstreamCasesTests(unittest.TestCase):
    def test_runner_executes_paired_cases_and_meets_expected_direction(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs, selection_path = build_selection(root)
            manifest_path = RUNNER.run_cases(
                case_manifest_path=selection_path,
                protocol_path=inputs["protocol_path"],
                output_root=root / "execution",
            )
            manifest = json.loads(manifest_path.read_text())
            validation = json.loads(
                (
                    root / "execution" / "validation" / "downstream_execution.json"
                ).read_text()
            )
            results = MODULE.runtime.read_csv(
                root / "execution" / "tables" / "paired_endpoints.csv"
            )

        self.assertTrue(manifest["valid"])
        self.assertTrue(manifest["scientific_targets_met"])
        self.assertEqual(validation["counts"]["condition_runs"], 6)
        self.assertEqual(validation["counts"]["paired_results"], 3)
        self.assertEqual(validation["counts"]["cases_with_expected_endpoint_change"], 3)
        self.assertTrue(all(row["expected_change_met"] == "True" for row in results))

    def test_runner_separates_valid_execution_from_missed_scientific_target(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = fixture(root)
            registry = json.loads(inputs["registry_path"].read_text())
            for entry in registry["cases"]:
                before_spec = next(
                    Path(item["path"])
                    for item in entry["before"]["inputs"]
                    if item["role"] == "seqspec"
                )
                after_spec = next(
                    Path(item["path"])
                    for item in entry["after"]["inputs"]
                    if item["role"] == "seqspec"
                )
                before_spec.write_text("0.7\n", encoding="utf-8")
                after_spec.write_text("0.6\n", encoding="utf-8")
            selection_path = MODULE.build_cases(
                output_root=root / "selection", **inputs
            )
            manifest_path = RUNNER.run_cases(
                case_manifest_path=selection_path,
                protocol_path=inputs["protocol_path"],
                output_root=root / "execution",
            )
            manifest = json.loads(manifest_path.read_text())
            validation = json.loads(
                (
                    root / "execution" / "validation" / "downstream_execution.json"
                ).read_text()
            )

        self.assertTrue(validation["valid"])
        self.assertFalse(manifest["scientific_targets_met"])
        self.assertEqual(validation["counts"]["cases_with_expected_endpoint_change"], 0)

    def test_runner_rejects_inputs_changed_after_selection(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs, selection_path = build_selection(root)
            reads = root / "inputs" / "CASE1" / "reads.fastq"
            reads.write_text("@changed\nAAAA\n+\nIIII\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "input hash changed"):
                RUNNER.run_cases(
                    case_manifest_path=selection_path,
                    protocol_path=inputs["protocol_path"],
                    output_root=root / "execution",
                )

    def test_runner_records_invalid_endpoint_and_returns_nonzero(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = fixture(root)
            script = root / "inputs" / "pipeline.py"
            script.write_text(
                pipeline_source().replace('"unit": "fraction"', '"unit": "count"'),
                encoding="utf-8",
            )
            selection_path = MODULE.build_cases(
                output_root=root / "selection", **inputs
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER_PATH),
                    "--case-manifest",
                    str(selection_path),
                    "--protocol",
                    str(inputs["protocol_path"]),
                    "--output-root",
                    str(root / "execution"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            validation = json.loads(
                (
                    root / "execution" / "validation" / "downstream_execution.json"
                ).read_text()
            )
            runs = MODULE.runtime.read_csv(
                root / "execution" / "tables" / "condition_runs.csv"
            )

        self.assertEqual(completed.returncode, 1)
        self.assertFalse(validation["valid"])
        self.assertEqual({row["status"] for row in runs}, {"invalid_endpoint"})


if __name__ == "__main__":
    unittest.main()
