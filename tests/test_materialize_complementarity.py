import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import test_materialize_perturbations as base_fixture


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "materialize_complementarity.py"
SOURCE_PROTOCOL = ROOT / "experiments" / "paper" / "protocol" / "complementarity.json"
SPEC = importlib.util.spec_from_file_location(
    "materialize_complementarity", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def create_protocol(path: Path) -> Path:
    protocol = MODULE.runtime.load_json(SOURCE_PROTOCOL)
    severities = {
        "Q01": {"phred": 2},
        "Q02": {"phred": 2, "cycle_fraction": 0.25},
        "Q03": {"phred": 2, "event_fraction": 0.1},
        "Q04": {"low_phred": 2, "high_phred": 40, "low_fraction": 0.1},
    }
    for operator in protocol["quality_controls"]:
        operator["severity_candidates"] = [severities[operator["operator_id"]]]
    for operator in protocol["composition_controls"]:
        operator["event_fraction"] = 0.1
    MODULE.runtime.write_json(path, protocol)
    return path


def create_base(root: Path, *, cohort_split: str) -> dict[str, Path]:
    inputs = base_fixture.create_inputs(root / "inputs", cohort_split=cohort_split)
    sample_args = {}
    if cohort_split == "evaluation":
        sample_args["sample_bundle_path"] = base_fixture.create_sample_bundle(
            inputs, root / "bundle" / "bundle.json"
        )
    else:
        sample_args["study_manifest_path"] = inputs["study"]
    manifest = base_fixture.MODULE.materialize_perturbations(
        inventory_manifest_path=inputs["inventory"],
        sampling_policy_path=inputs["policy"],
        perturbation_protocol_path=inputs["protocol"],
        seqspec_bin=inputs["seqspec"],
        yq_bin=inputs["yq"],
        output_root=root / "base",
        timeout_seconds=10,
        **sample_args,
    )
    return {
        "base": Path(manifest["manifest_path"]),
        "seqspec": inputs["seqspec"],
        "yq": inputs["yq"],
    }


def quality_policy(path: Path, protocol_path: Path) -> Path:
    protocol = MODULE.runtime.load_json(protocol_path)
    stable = {
        "schema_version": "0.1.0",
        "analysis_id": "quality-calibration",
        "protocol_sha256": MODULE.runtime.file_sha256(protocol_path),
        "frozen": True,
        "entries": [
            {
                "operator_id": value["operator_id"],
                "severity_index": 0,
                "severity": value["severity_candidates"][0],
            }
            for value in protocol["quality_controls"]
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
    return path


def controls(manifest: dict[str, object]) -> list[dict[str, object]]:
    payload = MODULE.runtime.load_json(Path(manifest["outputs"]["controls"]["path"]))
    return payload["controls"]


class MaterializeComplementarityTests(unittest.TestCase):
    def test_calibration_materialization_is_deterministic_and_preserves_sequences(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            base = create_base(root / "study", cohort_split="calibration")
            protocol = create_protocol(root / "protocol.json")
            first = MODULE.materialize_complementarity(
                base_materialization_path=base["base"],
                protocol_path=protocol,
                seqspec_bin=base["seqspec"],
                yq_bin=base["yq"],
                cohort_split="calibration",
                output_root=root / "first",
                timeout_seconds=10,
            )
            second = MODULE.materialize_complementarity(
                base_materialization_path=base["base"],
                protocol_path=protocol,
                seqspec_bin=base["seqspec"],
                yq_bin=base["yq"],
                cohort_split="calibration",
                output_root=root / "second",
                timeout_seconds=10,
            )
            rows = controls(first)
            base_manifest, base_conditions = MODULE.base_runner.load_materialization(
                base["base"]
            )
            clean = next(
                value for value in base_conditions if value["condition_kind"] == "clean"
            )
            clean_records = MODULE.perturb.read_fastq(Path(clean["inputs"][0]["path"]))
            quality_sequences_match = []
            for row in (
                value for value in rows if value["condition_kind"] == "quality"
            ):
                observed = MODULE.perturb.read_fastq(Path(row["inputs"][0]["path"]))
                quality_sequences_match.append(
                    [value[1] for value in clean_records]
                    == [value[1] for value in observed]
                )

        self.assertEqual(first["complementarity_id"], second["complementarity_id"])
        self.assertEqual(first["counts"]["controls_by_kind"]["quality"], 4)
        self.assertEqual(first["counts"]["controls_by_kind"]["composition"], 1)
        self.assertEqual(first["counts"]["controls_by_kind"]["schema_invalid"], 2)
        self.assertEqual(first["counts"]["skips"], 2)
        self.assertEqual(
            base_manifest["materialization_id"], first["base_materialization_id"]
        )
        self.assertEqual(quality_sequences_match, [True] * 4)
        self.assertTrue(
            all(
                value["observed_seqspec_check"] == "failure"
                for value in rows
                if value["condition_kind"] == "schema_invalid"
            )
        )

    def test_cli_applies_only_frozen_evaluation_severities(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            base = create_base(root / "study", cohort_split="evaluation")
            protocol = create_protocol(root / "protocol.json")
            policy = quality_policy(root / "quality_policy.json", protocol)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--base-materialization",
                    str(base["base"]),
                    "--protocol",
                    str(protocol),
                    "--seqspec-bin",
                    str(base["seqspec"]),
                    "--yq-bin",
                    str(base["yq"]),
                    "--cohort-split",
                    "evaluation",
                    "--quality-policy",
                    str(policy),
                    "--output-root",
                    str(root / "controls"),
                    "--timeout-seconds",
                    "10",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest = json.loads(
                Path(completed.stdout.strip()).read_text(encoding="utf-8")
            )
            rows = controls(manifest)

        quality = [value for value in rows if value["condition_kind"] == "quality"]
        self.assertEqual(len(quality), 4)
        self.assertEqual({value["severity_index"] for value in quality}, {0})

    def test_evaluation_rejects_a_missing_quality_policy(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            base = create_base(root / "study", cohort_split="evaluation")
            protocol = create_protocol(root / "protocol.json")
            with self.assertRaisesRegex(ValueError, "requires one quality policy"):
                MODULE.materialize_complementarity(
                    base_materialization_path=base["base"],
                    protocol_path=protocol,
                    seqspec_bin=base["seqspec"],
                    yq_bin=base["yq"],
                    cohort_split="evaluation",
                    output_root=root / "controls",
                    timeout_seconds=10,
                )


if __name__ == "__main__":
    unittest.main()
