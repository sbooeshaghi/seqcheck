import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "verify_igvf_audit_repeatability.py"
PROTOCOL_PATH = ROOT / "experiments" / "paper" / "protocol" / "igvf_audit.json"
SPEC = importlib.util.spec_from_file_location(
    "verify_igvf_audit_repeatability_test", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def write_audit_root(
    root: Path, run_id: str, values: dict[tuple[str, str], float]
) -> None:
    write_json(
        root / "manifests" / "study.json",
        {"audit_schema_version": "0.4.1", "run_id": run_id},
    )
    write_json(
        root / "validation" / "reconciliation.json",
        {"audit_schema_version": "0.4.1", "study_run_id": run_id, "valid": True},
    )
    rows = []
    for configuration, modality in sorted(values):
        rows.append(
            {
                "study_run_id": run_id,
                "configuration_accession": configuration,
                "modality": modality,
                "access_class": "public",
            }
        )
        write_json(
            root / "reports" / configuration / f"{MODULE.audit.slugify(modality)}.json",
            {
                "meta": {
                    "command": "check",
                    "modality": modality,
                    "requested_reads": 10000,
                },
                "results": [
                    {
                        "check": "length",
                        "observed": [
                            {
                                "name": "out_of_range_fraction",
                                "data": {
                                    "kind": "scalar",
                                    "unit": "fraction",
                                    "value": values[(configuration, modality)],
                                },
                            }
                        ],
                    }
                ],
                "audit_summary": {
                    "study_run_id": run_id,
                    "audit_invocation": [run_id],
                },
            },
        )
    with (root / "runs.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_protocol(path: Path, subset_count: int = 2) -> None:
    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    protocol["repeatability"]["subset_run_count"] = subset_count
    write_json(path, protocol)


class VerifyIgvfAuditRepeatabilityTests(unittest.TestCase):
    def test_cli_prepare_and_verify_ignores_audit_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            values = {("C1", "rna"): 0.1, ("C2", "atac"): 0.2, ("C3", "guide"): 0.3}
            write_audit_root(root / "primary", "primary-study", values)
            write_audit_root(root / "repeat", "repeat-study", values)
            protocol = root / "protocol.json"
            write_protocol(protocol)

            prepared = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "prepare",
                    "--primary-root",
                    str(root / "primary"),
                    "--audit-protocol",
                    str(protocol),
                    "--output-root",
                    str(root / "selection"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(prepared.returncode, 0, prepared.stderr)
            selection_path = Path(prepared.stdout.strip())

            verified = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "verify",
                    "--primary-root",
                    str(root / "primary"),
                    "--repeat-root",
                    str(root / "repeat"),
                    "--selection-manifest",
                    str(selection_path),
                    "--audit-protocol",
                    str(protocol),
                    "--output-root",
                    str(root / "verification"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(verified.returncode, 0, verified.stderr)
            validation = json.loads(Path(verified.stdout.strip()).read_text())
            configurations = (
                root / "selection" / "inputs" / "configuration_accessions.txt"
            ).read_text()

        self.assertTrue(validation["valid"])
        self.assertEqual(validation["counts"]["requested_runs"], 2)
        self.assertEqual(validation["counts"]["matched_runs"], 2)
        self.assertEqual(validation["match_fraction"], 1.0)
        self.assertEqual(len(configurations.splitlines()), 2)

    def test_selection_is_deterministic_and_metric_changes_fail(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            values = {("C1", "rna"): 0.1, ("C2", "atac"): 0.2, ("C3", "guide"): 0.3}
            write_audit_root(root / "primary", "primary-study", values)
            write_audit_root(root / "repeat", "repeat-study", values)
            protocol = root / "protocol.json"
            write_protocol(protocol)
            first_path = MODULE.prepare_repeatability_subset(
                primary_root=root / "primary",
                protocol_path=protocol,
                output_root=root / "first",
            )
            second_path = MODULE.prepare_repeatability_subset(
                primary_root=root / "primary",
                protocol_path=protocol,
                output_root=root / "second",
            )
            first = json.loads(first_path.read_text())
            second = json.loads(second_path.read_text())
            repeat_study_path = root / "repeat" / "manifests" / "study.json"
            repeat_study = json.loads(repeat_study_path.read_text())
            repeat_study["sampling_method"] = "reservoir"
            write_json(repeat_study_path, repeat_study)
            identity_validation_path = MODULE.verify_repeatability(
                primary_root=root / "primary",
                repeat_root=root / "repeat",
                selection_path=first_path,
                protocol_path=protocol,
                output_root=root / "identity-mismatch",
            )
            identity_validation = json.loads(identity_validation_path.read_text())
            del repeat_study["sampling_method"]
            write_json(repeat_study_path, repeat_study)
            selected = first["selected_runs"][0]
            repeat_report = (
                root
                / "repeat"
                / "reports"
                / selected["configuration_accession"]
                / f"{MODULE.audit.slugify(selected['modality'])}.json"
            )
            changed = json.loads(repeat_report.read_text())
            changed["results"][0]["observed"][0]["data"]["value"] = 0.9
            write_json(repeat_report, changed)
            validation_path = MODULE.verify_repeatability(
                primary_root=root / "primary",
                repeat_root=root / "repeat",
                selection_path=first_path,
                protocol_path=protocol,
                output_root=root / "mismatch",
            )
            validation = json.loads(validation_path.read_text())

        self.assertEqual(first["selection_id"], second["selection_id"])
        self.assertEqual(first["selected_runs"], second["selected_runs"])
        self.assertFalse(identity_validation["valid"])
        self.assertFalse(identity_validation["checks"]["repeat_study_identity_matches"])
        self.assertEqual(identity_validation["match_fraction"], 1.0)
        self.assertFalse(validation["valid"])
        self.assertEqual(validation["counts"]["mismatched_runs"], 1)
        self.assertFalse(validation["checks"]["repeat_match_target_met"])

    def test_tampered_selection_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            values = {("C1", "rna"): 0.1, ("C2", "atac"): 0.2}
            write_audit_root(root / "primary", "primary-study", values)
            protocol = root / "protocol.json"
            write_protocol(protocol)
            path = MODULE.prepare_repeatability_subset(
                primary_root=root / "primary",
                protocol_path=protocol,
                output_root=root / "selection",
            )
            value = json.loads(path.read_text())
            value["selected_runs"][0]["modality"] = "tampered"
            write_json(path, value)

            with self.assertRaisesRegex(ValueError, "ID does not reconcile"):
                MODULE.load_selection_manifest(path)


if __name__ == "__main__":
    unittest.main()
