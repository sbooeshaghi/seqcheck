import argparse
import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_phase0.py"
SAMPLER_PATH = Path(__file__).resolve().parents[1] / "scripts" / "sample_fastq.py"
SPEC = importlib.util.spec_from_file_location("run_phase0_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_executable(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")
    path.chmod(0o755)


def write_fastq(path: Path, count: int = 20) -> None:
    with path.open("w", encoding="ascii") as handle:
        for index in range(count):
            handle.write(f"@read-{index}\nACGT\n+\nIIII\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class RunPhase0Tests(unittest.TestCase):
    def test_stable_local_identity_uses_content_not_location(self) -> None:
        first = {
            "source": "/first/reads.fastq",
            "kind": "local",
            "status": "available",
            "size_bytes": 100,
            "sha256": "content",
        }
        second = {
            **first,
            "source": "/second/reads.fastq",
            "status": "cached",
        }

        self.assertEqual(
            MODULE.stable_source_identity(first),
            MODULE.stable_source_identity(second),
        )

    def test_functional_tool_identity_uses_configuration_content(self) -> None:
        first = {
            "version": "FastQC test",
            "executable_sha256": "executable",
            "configuration": {
                "mode": "custom",
                "path": "/first/limits.txt",
                "sha256": "limits",
            },
        }
        second = {
            **first,
            "configuration": {
                **first["configuration"],
                "path": "/second/limits.txt",
            },
        }

        self.assertEqual(
            MODULE.functional_tool_identity(first),
            MODULE.functional_tool_identity(second),
        )

    def test_integrated_gate_reconciles_cases_samples_and_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            local_spec = root / "legacy.yaml"
            local_fastq = root / "local.fastq"
            remote_spec = root / "current.yaml"
            remote_fastq = root / "remote.fastq"
            remote_resource = root / "resource.txt"
            local_spec.write_text("seqspec_version: 0.3.0\n", encoding="utf-8")
            remote_spec.write_text("seqspec_version: 0.5.0\n", encoding="utf-8")
            remote_resource.write_text("ACGT\n", encoding="utf-8")
            write_fastq(local_fastq)
            write_fastq(remote_fastq)

            seqspec = root / "seqspec"
            write_executable(
                seqspec,
                """#!/usr/bin/env python3
import re
import sys
from pathlib import Path

args = sys.argv[1:]
if args == ["--version"]:
    print("seqspec test 0.1")
elif args[0] == "version":
    text = Path(args[1]).read_text()
    version = re.search(r"seqspec_version: ([0-9.]+)", text).group(1)
    print(f"seqspec file version: {version}")
elif args[0] == "upgrade":
    source = Path(args[1])
    output = Path(args[args.index("-o") + 1])
    output.write_text(source.read_text().replace("0.3.0", "0.5.0"))
elif args[0] != "check":
    raise SystemExit(2)
""",
            )
            seqcheck = root / "seqcheck"
            write_executable(
                seqcheck,
                """#!/usr/bin/env python3
import json
import sys
from pathlib import Path

args = sys.argv[1:]
if args == ["--version"]:
    print("seqcheck test 0.1")
    raise SystemExit()
spec = Path(args[args.index("--spec") + 1])
output = Path(args[args.index("--output") + 1])
resource = spec.parent / "resource.txt"
expected = []
if resource.exists():
    expected.append({
        "id": "e1",
        "name": "onlist_source",
        "description": "Declared resource.",
        "data": {"kind": "scalar", "value": str(resource)},
    })
payload = {
    "report_schema_version": "0.1.0",
    "meta": {"spec": str(spec), "requested_reads": 0},
    "results": [{
        "check": "coverage",
        "files": ["reads.fastq"],
        "reads": ["R1"],
        "regions": ["region"],
        "expected": expected,
        "observed": [{
            "id": "o1",
            "name": "covered_fraction",
            "description": "Covered fraction.",
            "data": {"kind": "scalar", "value": 1.0, "unit": "fraction"},
        }],
        "assessment": [{
            "type": "pass",
            "code": "full_coverage",
            "description": "All reads cover the region.",
        }],
    }],
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload))
""",
            )
            fastqc = root / "fastqc"
            write_executable(
                fastqc,
                "#!/bin/sh\necho 'FastQC test 0.1'\n",
            )

            protocol_path = root / "cases.json"
            protocol_path.write_text(
                json.dumps(
                    {
                        "schema_version": "0.1.0",
                        "cases": [
                            {
                                "case_id": "legacy-local",
                                "spec": local_spec.name,
                                "modality": "rna",
                                "fastqs": [local_fastq.name],
                                "n_reads": 0,
                                "resources": [],
                            },
                            {
                                "case_id": "current-remote",
                                "spec": "https://fixture/spec.yaml",
                                "modality": "rna",
                                "fastqs": ["https://fixture/remote.fastq"],
                                "n_reads": 0,
                                "resources": [
                                    {
                                        "source": "https://fixture/resource.txt",
                                        "path": "resource.txt",
                                    }
                                ],
                            },
                        ],
                        "sampling_probe": {
                            "inputs": [local_fastq.name],
                            "n_reads": 5,
                            "seed": 11,
                            "synchronize_mates": False,
                        },
                    }
                ),
                encoding="utf-8",
            )
            remote_payloads = {
                "https://fixture/spec.yaml": remote_spec,
                "https://fixture/resource.txt": remote_resource,
            }
            original_materialize = MODULE.materialize_source
            original_identity = MODULE.source_identity

            def materialize(source: str, destination: Path) -> dict[str, object]:
                if source not in remote_payloads:
                    return original_materialize(source, destination)
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(remote_payloads[source].read_bytes())
                return {
                    "source": source,
                    "kind": "remote",
                    "status": 200,
                    "etag": '"fixture"',
                    "last_modified": "",
                    "content_length": str(destination.stat().st_size),
                    "retrieved_at": "2026-07-14T00:00:00+00:00",
                    "materialized_path": str(destination.resolve()),
                    "materialized_sha256": MODULE.file_sha256(destination),
                }

            def source_identity(source: str) -> dict[str, object]:
                if source == "https://fixture/remote.fastq":
                    return {
                        "source": source,
                        "kind": "remote",
                        "status": 200,
                        "etag": '"fastq-fixture"',
                        "last_modified": "",
                        "content_length": str(remote_fastq.stat().st_size),
                        "retrieved_at": "2026-07-14T00:00:00+00:00",
                    }
                return original_identity(source)

            output_root = root / "output"
            args = argparse.Namespace(
                case_manifest=protocol_path,
                output_root=output_root,
                seqcheck_bin=seqcheck,
                seqspec_bin=seqspec,
                sampler_script=SAMPLER_PATH,
                fastqc_command=str(fastqc),
                fastqc_limits=None,
                command_timeout_seconds=30,
            )
            with (
                patch.object(MODULE, "materialize_source", side_effect=materialize),
                patch.object(MODULE, "source_identity", side_effect=source_identity),
            ):
                result = MODULE.run(args)
                second_output = root / "second-output"
                second_result = MODULE.run(
                    argparse.Namespace(**{**vars(args), "output_root": second_output})
                )

            self.assertEqual(result, 0)
            self.assertEqual(second_result, 0)
            validation = json.loads(
                (output_root / "validation" / "phase0.json").read_text()
            )
            self.assertTrue(validation["valid"])
            self.assertTrue(all(validation["checks"].values()))
            self.assertEqual(validation["counts"]["reports"], 4)
            runs = read_csv(output_root / "tables" / "runs.csv")
            metrics = read_csv(output_root / "tables" / "metrics.csv")
            self.assertEqual(len(runs), 4)
            self.assertEqual(
                len(metrics), sum(int(row["metric_count"]) for row in runs)
            )
            sample_probe = json.loads(
                (output_root / "manifests" / "study.json").read_text()
            )["sampling_probe"]
            self.assertEqual(
                sample_probe["manifests"]["repeat_a"]["output_sha256"],
                sample_probe["manifests"]["repeat_b"]["output_sha256"],
            )
            first_study = json.loads(
                (output_root / "manifests" / "study.json").read_text()
            )
            second_study = json.loads(
                (second_output / "manifests" / "study.json").read_text()
            )
            self.assertEqual(first_study["study_run_id"], second_study["study_run_id"])

    def test_protocol_rejects_unsafe_resource_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "cases.json"
            path.write_text(
                json.dumps(
                    {
                        "schema_version": "0.1.0",
                        "cases": [
                            {
                                "case_id": "unsafe",
                                "spec": "spec.yaml",
                                "modality": "rna",
                                "fastqs": ["reads.fastq"],
                                "resources": [
                                    {"source": "resource.txt", "path": "../escape"}
                                ],
                            }
                        ],
                        "sampling_probe": {
                            "inputs": ["reads.fastq"],
                            "n_reads": 1,
                            "seed": 1,
                        },
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "unsafe resource path"):
                MODULE.load_protocol(path)

    def test_cli_help_exposes_required_tool_identities(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), "--help"],
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--seqcheck-bin", completed.stdout)
        self.assertIn("--seqspec-bin", completed.stdout)
        self.assertIn("--fastqc-command", completed.stdout)


if __name__ == "__main__":
    unittest.main()
