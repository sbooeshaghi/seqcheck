import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "run_analysis_freeze.py"
SPEC = importlib.util.spec_from_file_location("run_analysis_freeze_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def analysis_source() -> str:
    return """#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output-root", required=True)
parser.add_argument("--mode", required=True)
parser.add_argument("--tag", required=True)
args = parser.parse_args()
source = Path(args.input)
value = source.read_text(encoding="utf-8").strip()
root = Path(args.output_root)
(root / "tables").mkdir(parents=True)
(root / "figure_data").mkdir(parents=True)
(root / "figures").mkdir(parents=True)
(root / "tables" / "summary.csv").write_text(
    "metric,value\\nobserved," + value + "\\n", encoding="utf-8"
)
figure_value = args.tag if args.mode == "nondeterministic" else value
(root / "figure_data" / "summary.json").write_text(
    json.dumps({"value": figure_value}, sort_keys=True) + "\\n", encoding="utf-8"
)
figure = args.tag if args.mode == "figure-nondeterministic" else value
(root / "figures" / "summary.txt").write_text(figure + "\\n", encoding="utf-8")
if args.mode == "extra":
    (root / "tables" / "undeclared.csv").write_text("x\\n1\\n", encoding="utf-8")
elif args.mode == "symlink":
    (root / "tables" / "summary.csv").unlink()
    (root / "tables" / "summary.csv").symlink_to(source)
elif args.mode == "mutate-input":
    source.write_text("changed\\n", encoding="utf-8")
elif args.mode == "delete-input":
    source.unlink()
elif args.mode == "mutate-study-manifest":
    Path(__file__).with_name("study.json").write_text("{}\\n", encoding="utf-8")
"""


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def make_study(root: Path, mode: str = "deterministic") -> tuple[Path, Path]:
    source = root / "source.txt"
    source.write_text("7\n", encoding="utf-8")
    script = root / "analysis.py"
    script.write_text(analysis_source(), encoding="utf-8")
    python = Path(sys.executable).resolve()
    manifest = {
        "schema_version": "0.1.0",
        "study_id": "paper-test",
        "repetitions": 2,
        "canonical_output_kinds": ["table", "figure_data"],
        "result_directory_names": [
            "tables",
            "analysis",
            "figure_data",
            "figures",
        ],
        "environment": {
            "LC_ALL": "C",
            "PYTHONHASHSEED": "0",
            "TZ": "UTC",
        },
        "inputs": [{"id": "source", "path": str(source)}],
        "repositories": [],
        "software": [
            {
                "id": "python",
                "path": str(python),
                "version": sys.version.split()[0],
                "repository_id": None,
            },
            {
                "id": "analysis-script",
                "path": str(script),
                "version": "test",
                "repository_id": None,
            },
        ],
        "steps": [
            {
                "id": "analysis",
                "working_directory": str(root),
                "timeout_seconds": 10,
                "depends_on": [],
                "input_ids": ["source"],
                "software_ids": ["python", "analysis-script"],
                "command": [
                    str(python),
                    str(script),
                    "--input",
                    str(source),
                    "--output-root",
                    "{step_root}",
                    "--mode",
                    mode,
                    "--tag",
                    "{reproduction}",
                ],
                "outputs": [
                    {
                        "id": "summary-table",
                        "path": "tables/summary.csv",
                        "kind": "table",
                    },
                    {
                        "id": "summary-figure-data",
                        "path": "figure_data/summary.json",
                        "kind": "figure_data",
                    },
                    {
                        "id": "summary-figure",
                        "path": "figures/summary.txt",
                        "kind": "figure",
                    },
                ],
            }
        ],
    }
    manifest_path = root / "study.json"
    write_json(manifest_path, manifest)
    return manifest_path, source


class RunAnalysisFreezeTests(unittest.TestCase):
    def test_freeze_runs_twice_and_archives_matching_results(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root)
            output_root = root / "freeze"
            result_path = MODULE.run_analysis_freeze(
                study_manifest_path=manifest_path,
                output_root=output_root,
            )
            result = json.loads(result_path.read_text())
            validation = json.loads(
                (output_root / "validation" / "analysis_freeze.json").read_text()
            )
            comparisons = MODULE.runtime.read_csv(
                output_root / "tables" / "result_hashes.csv"
            )
            bundle = json.loads(
                (output_root / "bundle" / "manifests" / "bundle.json").read_text()
            )
            archived_runs = output_root / "bundle" / "evidence" / "analysis_runs.csv"
            archived_hashes = output_root / "bundle" / "evidence" / "result_hashes.csv"
            archived_logs = output_root / "bundle" / "manifests" / "logs.json"
            archived_evidence_exists = all(
                path.is_file()
                for path in (archived_runs, archived_hashes, archived_logs)
            )

        self.assertTrue(result["valid"])
        self.assertTrue(result["frozen"])
        self.assertTrue(validation["checks"]["all_canonical_outputs_byte_identical"])
        self.assertEqual(validation["counts"]["run_rows"], 2)
        self.assertEqual(validation["counts"]["canonical_outputs"], 2)
        self.assertEqual(len(comparisons), 3)
        self.assertTrue(all(row["byte_identical"] == "True" for row in comparisons))
        self.assertEqual(bundle["analysis_freeze_id"], result["analysis_freeze_id"])
        self.assertTrue(archived_evidence_exists)

    def test_noncanonical_figure_may_differ_without_failing_freeze(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root, "figure-nondeterministic")
            output_root = root / "freeze"
            MODULE.run_analysis_freeze(
                study_manifest_path=manifest_path,
                output_root=output_root,
            )
            rows = MODULE.runtime.read_csv(output_root / "tables" / "result_hashes.csv")

        figure = next(row for row in rows if row["kind"] == "figure")
        self.assertEqual(figure["byte_identical"], "False")
        self.assertEqual(figure["canonical"], "False")

    def test_ordered_step_can_consume_prior_reproduction_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, source = make_study(root)
            manifest = json.loads(manifest_path.read_text())
            second = json.loads(json.dumps(manifest["steps"][0]))
            second["id"] = "figures"
            second["depends_on"] = ["analysis"]
            second["input_ids"] = []
            second["command"][second["command"].index(str(source))] = (
                "{analysis_root}/analysis/tables/summary.csv"
            )
            for output in second["outputs"]:
                output["id"] = "second-" + output["id"]
            manifest["steps"].append(second)
            write_json(manifest_path, manifest)
            output_root = root / "freeze"
            result_path = MODULE.run_analysis_freeze(
                study_manifest_path=manifest_path,
                output_root=output_root,
            )
            result = json.loads(result_path.read_text())

        self.assertEqual(result["counts"]["analysis_steps"], 2)
        self.assertEqual(result["counts"]["run_rows"], 4)
        self.assertEqual(result["counts"]["declared_outputs"], 6)

    def test_canonical_difference_fails_with_validation_record(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root, "nondeterministic")
            output_root = root / "freeze"
            with self.assertRaisesRegex(ValueError, "canonical_outputs_byte_identical"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=output_root,
                )
            validation = json.loads(
                (output_root / "validation" / "analysis_freeze.json").read_text()
            )

        self.assertFalse(validation["valid"])
        self.assertFalse(validation["frozen"])
        self.assertFalse(validation["checks"]["all_canonical_outputs_byte_identical"])

    def test_bundle_failure_writes_fail_closed_validation(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root)
            output_root = root / "freeze"
            with (
                mock.patch.object(
                    MODULE,
                    "build_bundle",
                    side_effect=OSError("archive unavailable"),
                ),
                self.assertRaisesRegex(ValueError, "bundle_archived_when_eligible"),
            ):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=output_root,
                )
            validation = json.loads(
                (output_root / "validation" / "analysis_freeze.json").read_text()
            )

        self.assertFalse(validation["valid"])
        self.assertFalse(validation["frozen"])
        self.assertEqual(validation["analysis_freeze_id"], "")
        self.assertFalse(validation["checks"]["bundle_archived_when_eligible"])
        self.assertEqual(
            validation["postflight_errors"]["bundle"], "archive unavailable"
        )

    def test_undeclared_result_file_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root, "extra")
            output_root = root / "freeze"
            with self.assertRaisesRegex(ValueError, "all_analysis_steps_completed"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=output_root,
                )
            rows = MODULE.runtime.read_csv(output_root / "tables" / "analysis_runs.csv")

        self.assertEqual({row["status"] for row in rows}, {"invalid_outputs"})
        self.assertTrue(
            all("undeclared result files" in row["message"] for row in rows)
        )

    def test_output_symlink_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root, "symlink")
            output_root = root / "freeze"
            with self.assertRaisesRegex(ValueError, "all_analysis_steps_completed"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=output_root,
                )
            rows = MODULE.runtime.read_csv(output_root / "tables" / "analysis_runs.csv")

        self.assertEqual({row["status"] for row in rows}, {"invalid_outputs"})
        self.assertTrue(all("symbolic link" in row["message"] for row in rows))

    def test_changed_source_artifact_fails_post_run_identity_check(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, source = make_study(root, "mutate-input")
            output_root = root / "freeze"
            with self.assertRaisesRegex(ValueError, "source_artifacts_unchanged"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=output_root,
                )
            validation = json.loads(
                (output_root / "validation" / "analysis_freeze.json").read_text()
            )
            observed_source = source.read_text()

        self.assertEqual(observed_source, "changed\n")
        self.assertFalse(validation["checks"]["all_source_artifacts_unchanged"])

    def test_changed_study_manifest_fails_post_run_identity_check(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root, "mutate-study-manifest")
            output_root = root / "freeze"
            with self.assertRaisesRegex(ValueError, "study_manifest_unchanged"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=output_root,
                )
            validation = json.loads(
                (output_root / "validation" / "analysis_freeze.json").read_text()
            )

        self.assertFalse(validation["valid"])
        self.assertFalse(validation["checks"]["study_manifest_unchanged"])

    def test_deleted_source_still_writes_fail_closed_validation(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root, "delete-input")
            output_root = root / "freeze"
            with self.assertRaisesRegex(ValueError, "source_artifacts_unchanged"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=output_root,
                )
            validation = json.loads(
                (output_root / "validation" / "analysis_freeze.json").read_text()
            )

        self.assertFalse(validation["valid"])
        self.assertIn(
            "source artifact is missing", validation["postflight_errors"]["artifacts"]
        )

    def test_manifest_rejects_top_level_symlink_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, source = make_study(root)
            symlink = root / "source-link.txt"
            symlink.symlink_to(source)
            manifest = json.loads(manifest_path.read_text())
            manifest["inputs"][0]["path"] = str(symlink)
            command = manifest["steps"][0]["command"]
            command[command.index(str(source))] = str(symlink)
            write_json(manifest_path, manifest)

            with self.assertRaisesRegex(ValueError, "missing or symbolic"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=root / "freeze",
                )

        self.assertFalse((root / "freeze").exists())

    def test_manifest_rejects_non_object_document(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path = root / "study.json"
            write_json(manifest_path, [])

            with self.assertRaisesRegex(ValueError, "expected a JSON object"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=root / "freeze",
                )

        self.assertFalse((root / "freeze").exists())

    def test_manifest_rejects_command_that_omits_declared_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, source = make_study(root)
            manifest = json.loads(manifest_path.read_text())
            command = manifest["steps"][0]["command"]
            command[command.index(str(source))] = str(root / "another-input")
            write_json(manifest_path, manifest)

            with self.assertRaisesRegex(ValueError, "command omits input"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=root / "freeze",
                )
            self.assertFalse((root / "freeze").exists())

    def test_manifest_rejects_unregistered_absolute_command_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root)
            unregistered = root / "unregistered.json"
            unregistered.write_text("{}\n", encoding="utf-8")
            manifest = json.loads(manifest_path.read_text())
            manifest["steps"][0]["command"].extend(
                ["--unregistered", str(unregistered)]
            )
            write_json(manifest_path, manifest)

            with self.assertRaisesRegex(ValueError, "command uses undeclared path"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=root / "freeze",
                )
            self.assertFalse((root / "freeze").exists())

    def test_manifest_requires_deterministic_environment(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root)
            manifest = json.loads(manifest_path.read_text())
            manifest["environment"]["PYTHONHASHSEED"] = "random"
            write_json(manifest_path, manifest)

            with self.assertRaisesRegex(ValueError, "environment must set"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=root / "freeze",
                )
            self.assertFalse((root / "freeze").exists())

    def test_output_root_cannot_be_inside_immutable_input_tree(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root)
            manifest = json.loads(manifest_path.read_text())
            manifest["inputs"][0]["path"] = str(root)
            write_json(manifest_path, manifest)

            with self.assertRaisesRegex(ValueError, "inside immutable input"):
                MODULE.run_analysis_freeze(
                    study_manifest_path=manifest_path,
                    output_root=root / "freeze",
                )
            self.assertFalse((root / "freeze").exists())

    def test_repository_snapshot_rejects_tracked_changes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir).resolve()
            repository = root / "repository"
            repository.mkdir()
            subprocess.run(["git", "init", "-q", str(repository)], check=True)
            subprocess.run(
                ["git", "-C", str(repository), "config", "user.name", "Test"],
                check=True,
            )
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(repository),
                    "config",
                    "user.email",
                    "test@example.org",
                ],
                check=True,
            )
            tracked = repository / "analysis.py"
            tracked.write_text("print('stable')\n", encoding="utf-8")
            subprocess.run(
                ["git", "-C", str(repository), "add", "analysis.py"], check=True
            )
            subprocess.run(
                ["git", "-C", str(repository), "commit", "-q", "-m", "initial"],
                check=True,
            )
            clean = MODULE.snapshot_repositories(
                [{"id": "analysis", "path": repository}]
            )
            untracked = repository / "untracked.py"
            untracked.write_text("print('untracked')\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not tracked by Git"):
                MODULE.snapshot_software(
                    [
                        {
                            "id": "untracked",
                            "path": untracked,
                            "version": "test",
                            "repository_id": "analysis",
                        }
                    ],
                    clean,
                )
            tracked.write_text("print('changed')\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "tracked changes"):
                MODULE.snapshot_repositories([{"id": "analysis", "path": repository}])

        self.assertTrue(clean[0]["tracked_clean"])

    def test_cli_returns_nonzero_and_preserves_invalid_validation(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root, "nondeterministic")
            output_root = root / "freeze"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--study-manifest",
                    str(manifest_path),
                    "--output-root",
                    str(output_root),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            validation = json.loads(
                (output_root / "validation" / "analysis_freeze.json").read_text()
            )

        self.assertEqual(completed.returncode, 1)
        self.assertFalse(validation["valid"])

    def test_cli_freezes_valid_study(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_path, _ = make_study(root)
            output_root = root / "freeze"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--study-manifest",
                    str(manifest_path),
                    "--output-root",
                    str(output_root),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            result_path = Path(completed.stdout.strip())
            result = json.loads(result_path.read_text())

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(result["valid"])
        self.assertTrue(result["frozen"])


if __name__ == "__main__":
    unittest.main()
