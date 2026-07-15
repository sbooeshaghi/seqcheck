import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "freeze_igvf_audit_portal.py"
SPEC = importlib.util.spec_from_file_location(
    "freeze_igvf_audit_portal_test", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class FreezeIgvfAuditPortalTests(unittest.TestCase):
    def test_cli_freezes_only_linked_records_from_snapshots(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            configurations = root / "configurations.json"
            sequence_files = root / "sequence-files.json"
            configurations.write_text(
                json.dumps(
                    {
                        "@graph": [
                            {
                                "accession": "IGVFFI0001TEST",
                                "href": "/configuration-files/IGVFFI0001TEST/spec.yaml.gz",
                                "lab": {"title": "Test Lab"},
                                "submitted_by": {"title": "Tester"},
                                "award": {"component": "mapping"},
                                "file_set": {
                                    "accession": "IGVFDS0001TEST",
                                    "assay_term": {"term_name": "RNA sequencing"},
                                },
                                "preferred_assay_titles": ["scRNA-seq"],
                                "aliases": ["test:one"],
                                "seqspec_of": [
                                    "/sequence-files/IGVFFI0002TEST/",
                                    "/sequence-files/IGVFFI0003MISS/",
                                ],
                                "status": "released",
                                "upload_status": "validated",
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            sequence_files.write_text(
                json.dumps(
                    {
                        "@graph": [
                            {
                                "accession": "IGVFFI0002TEST",
                                "href": "/sequence-files/IGVFFI0002TEST/read.fastq.gz",
                                "controlled_access": False,
                                "read_names": ["Read1"],
                            },
                            {
                                "accession": "IGVFFI9999OTHER",
                                "href": "/sequence-files/IGVFFI9999OTHER/read.fastq.gz",
                                "controlled_access": False,
                                "read_names": ["Read2"],
                            },
                        ]
                    }
                ),
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--configuration-snapshot",
                    str(configurations),
                    "--sequence-file-snapshot",
                    str(sequence_files),
                    "--output-root",
                    str(root / "frozen"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            path = Path(completed.stdout.strip())
            manifest = json.loads(path.read_text(encoding="utf-8"))
            records, fastqs, _ = MODULE.audit.load_portal_manifest(path)

        self.assertEqual([record.accession for record in records], ["IGVFFI0001TEST"])
        self.assertEqual(list(fastqs), ["IGVFFI0002TEST"])
        self.assertEqual(
            manifest["missing_linked_sequence_file_accessions"],
            ["IGVFFI0003MISS"],
        )
        self.assertEqual(manifest["counts"]["linked_sequence_files"], 1)

    def test_manifest_id_ignores_source_paths_and_retrieval_time(self) -> None:
        configuration = MODULE.audit.ConfigurationRecord(
            accession="C",
            href="/c",
            lab="lab",
            submitted_by="user",
            award_component="mapping",
            file_set_accession="D",
            assay_term="assay",
            preferred_assay_titles=[],
            aliases=[],
            seqspec_of=["/sequence-files/F/"],
            status="released",
            upload_status="validated",
        )
        sequence = MODULE.audit.SequenceFileRecord("F", "/f", False, ["Read1"])
        tools = {
            "audit_runtime": {
                "version": "0.4.1",
                "sha256": "a" * 64,
                "python": "3.12.0",
            },
            "freezer": {
                "version": "0.1.0",
                "sha256": "b" * 64,
                "python": "3.12.0",
            },
        }
        query = {
            "api_root": "https://api.example.org/",
            "portal_root": "https://example.org/",
            "status": "released",
            "upload_status": "validated",
        }
        first = MODULE.build_portal_manifest(
            configurations=[configuration],
            sequence_files={"F": sequence},
            query=query,
            tools=tools,
            source={"path": "/first"},
            retrieved_at="first",
        )
        second = MODULE.build_portal_manifest(
            configurations=[configuration],
            sequence_files={"F": sequence},
            query=query,
            tools=tools,
            source={"path": "/second"},
            retrieved_at="second",
        )

        self.assertEqual(first["portal_manifest_id"], second["portal_manifest_id"])


if __name__ == "__main__":
    unittest.main()
