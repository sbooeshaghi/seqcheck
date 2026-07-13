import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "igvf_audit.py"
SPEC = importlib.util.spec_from_file_location("igvf_audit", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class IgvfAuditTests(unittest.TestCase):
    @staticmethod
    def configuration_record() -> object:
        return MODULE.ConfigurationRecord(
            accession="IGVFFI1234TEST",
            href="/configuration-files/IGVFFI1234TEST/@@download/IGVFFI1234TEST.yaml.gz",
            lab="Test Lab",
            submitted_by="Tester",
            award_component="mapping",
            file_set_accession="IGVFDS1234TEST",
            assay_term="single-cell RNA sequencing assay",
            preferred_assay_titles=["scRNA-seq"],
            aliases=["test:seqspec"],
            seqspec_of=["/sequence-files/IGVFFI0001TEST/"],
            status="released",
            upload_status="validated",
        )

    @staticmethod
    def tool_identity(version: str = "tool 1.0") -> object:
        return MODULE.ToolIdentity(
            command=["tool"],
            version=version,
            git_commit="abc123",
            git_dirty=False,
            runtime_source_sha256="source-hash",
            executable_sha256="binary-hash",
        )

    @classmethod
    def audit_context(cls, run_id: str = "study-1") -> object:
        return MODULE.AuditContext(
            run_id=run_id,
            sampling_method="prefix",
            sampling_seed=None,
            seqcheck=cls.tool_identity("seqcheck 0.2.0"),
            seqspec=cls.tool_identity("seqspec 0.4.0"),
        )

    def test_parse_seqspec_version_output(self) -> None:
        output = "seqspec version: 0.4.2\nseqspec file version: 0.3.0\n"
        self.assertEqual(MODULE.parse_seqspec_version_output(output), "0.3.0")

    def test_normalize_seqspec_version_maps_legacy_versions(self) -> None:
        self.assertEqual(MODULE.normalize_seqspec_version("0.3.0"), "0.5.0")
        self.assertEqual(MODULE.normalize_seqspec_version("0.4.0"), "0.5.0")
        self.assertEqual(MODULE.normalize_seqspec_version("0.5.0"), "0.5.0")

    def test_select_requested_configurations_rejects_missing_accessions(self) -> None:
        record = self.configuration_record()
        selected = MODULE.select_requested_configurations(
            [record], [record.accession]
        )
        self.assertEqual(selected, [record])

        with self.assertRaisesRegex(ValueError, "IGVFFI_MISSING"):
            MODULE.select_requested_configurations(
                [record], [record.accession, "IGVFFI_MISSING"]
            )

    def test_classify_seqspec_failure_reason_marks_malformed_yaml(self) -> None:
        self.assertEqual(
            MODULE.classify_seqspec_failure_reason(
                "seqspec_version_error", "Could not read values."
            ),
            "malformed_seqspec_yaml",
        )
        self.assertEqual(
            MODULE.classify_seqspec_failure_reason(
                "seqspec_version_error", "Input source does not exist: https://example.org/spec.yaml"
            ),
            "seqspec_version_error",
        )

    def test_resolve_ready_auth_profile_requires_config_and_env(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "auth.toml"
            config_path.write_text(
                "[profiles.igvf]\n"
                'hosts = ["data.igvf.org"]\n'
                'kind = "basic"\n'
                'username_env = "TEST_USER"\n'
                'password_env = "TEST_PASS"\n',
                encoding="utf-8",
            )

            with patch.dict(
                os.environ,
                {
                    "SEQCHECK_AUTH_CONFIG": str(config_path),
                    "TEST_USER": "user",
                    "TEST_PASS": "pass",
                },
                clear=False,
            ):
                self.assertEqual(
                    MODULE.resolve_ready_auth_profile(
                        "igvf",
                        config_env="SEQCHECK_AUTH_CONFIG",
                        config_subdir="seqcheck",
                    ),
                    "igvf",
                )

            with patch.dict(
                os.environ,
                {"SEQCHECK_AUTH_CONFIG": str(config_path)},
                clear=False,
            ):
                self.assertIsNone(
                    MODULE.resolve_ready_auth_profile(
                        "igvf",
                        config_env="SEQCHECK_AUTH_CONFIG",
                        config_subdir="seqcheck",
                    )
                )

    def test_is_fastq_expectation_accepts_compressed_fastq_filetypes(self) -> None:
        self.assertTrue(
            MODULE.is_fastq_expectation(
                {"filetype": "fastq.gz", "filename": "IGVFFI0001.fastq.gz"}
            )
        )
        self.assertTrue(
            MODULE.is_fastq_expectation(
                {"filetype": ".fq.gz", "filename": "IGVFFI0002.fq.gz"}
            )
        )
        self.assertTrue(
            MODULE.is_fastq_expectation(
                {"filetype": "", "filename": "IGVFFI0003.fastq.gz"}
            )
        )
        self.assertFalse(
            MODULE.is_fastq_expectation(
                {"filetype": "bam", "filename": "IGVFFI0004.bam"}
            )
        )

    def test_flatten_report_counts_assessments(self) -> None:
        record = self.configuration_record()
        report = {
            "meta": {"requested_reads": 10000},
            "audit_summary": {
                "study_run_id": "study-1",
                "cache_key": "cache-1",
                "sampling_method": "prefix",
                "sampling_seed": None,
            },
            "results": [
                {
                    "check": "input_check",
                    "files": ["IGVFFI0001TEST"],
                    "reads": ["Read1"],
                    "regions": [],
                    "expected": [
                        {
                            "id": "e1",
                            "name": "expected_regions",
                            "description": "Expected region annotations.",
                            "data": {
                                "kind": "records",
                                "value": [
                                    {
                                        "region_type": [
                                            "RGN:partition:cell",
                                            "RGN:classify:cell_identity",
                                        ]
                                    }
                                ],
                            },
                        }
                    ],
                    "observed": [
                        {
                            "id": "o1",
                            "name": "matched_count",
                            "description": "Matched input count.",
                            "data": {"kind": "scalar", "value": 1, "unit": "count"},
                        }
                    ],
                    "assessment": [
                        {"type": "pass", "code": "all_expected_files_matched", "description": "ok"}
                    ],
                },
                {
                    "check": "length",
                    "files": ["IGVFFI0001TEST"],
                    "reads": ["Read1"],
                    "regions": [],
                    "assessment": [
                        {"type": "warning", "code": "length_out_of_range", "description": "warn"},
                        {"type": "interpretation", "code": "note", "description": "info"},
                    ],
                },
            ],
        }

        run, diagnostics, metrics = MODULE.flatten_report(
            report=report,
            record=record,
            modality="rna",
            raw_seqspec_version="0.3.0",
            normalized_seqspec_version="0.5.0",
            report_path=Path("/tmp/report.json"),
            modality_count=2,
            expected_fastq_count=1,
            supplied_fastq_count=1,
            controlled_access=False,
            run_status="completed",
        )

        self.assertEqual(run.modality_count, 2)
        self.assertEqual(run.pass_count, 1)
        self.assertEqual(run.warning_count, 1)
        self.assertEqual(run.error_count, 0)
        self.assertEqual(run.interpretation_count, 1)
        self.assertEqual(run.requested_reads_total, 10000)
        self.assertEqual(run.supplied_fastq_count, 1)
        self.assertEqual(run.study_run_id, "study-1")
        self.assertEqual(len(diagnostics), 3)
        self.assertEqual(diagnostics[1].assessment_code, "length_out_of_range")
        self.assertEqual(len(metrics), 2)
        self.assertEqual(metrics[0].metric_side, "expected")
        self.assertEqual(metrics[0].data_kind, "records")
        self.assertEqual(
            metrics[0].ontology_terms,
            "RGN:classify:cell_identity;RGN:partition:cell",
        )
        self.assertEqual(metrics[1].value_json, "1")

    def test_annotate_report_summary_adds_counts(self) -> None:
        report = {"meta": {"requested_reads": 10000}, "results": []}
        context = self.audit_context()
        with patch.object(sys, "argv", ["igvf_audit.py", "--limit", "1"]):
            annotated = MODULE.annotate_report_summary(
                report,
                configuration_accession="IGVFFI1234TEST",
                modality="rna",
                modality_count=3,
                expected_fastq_count=2,
                supplied_fastq_count=2,
                controlled_access=True,
                audit_context=context,
                cache_key="cache-1",
            )

        self.assertEqual(
            annotated["audit_summary"],
            {
                "audit_schema_version": "0.2.0",
                "study_run_id": "study-1",
                "cache_key": "cache-1",
                "sampling_method": "prefix",
                "sampling_seed": None,
                "audit_invocation": ["igvf_audit.py", "--limit", "1"],
                "tools": {
                    "seqcheck": {
                        "command": ["tool"],
                        "version": "seqcheck 0.2.0",
                        "git_commit": "abc123",
                        "git_dirty": False,
                        "runtime_source_sha256": "source-hash",
                        "executable_sha256": "binary-hash",
                    },
                    "seqspec": {
                        "command": ["tool"],
                        "version": "seqspec 0.4.0",
                        "git_commit": "abc123",
                        "git_dirty": False,
                        "runtime_source_sha256": "source-hash",
                        "executable_sha256": "binary-hash",
                    },
                },
                "configuration_accession": "IGVFFI1234TEST",
                "modality": "rna",
                "modality_count": 3,
                "expected_fastq_count": 2,
                "supplied_fastq_count": 2,
                "requested_reads_per_fastq": 10000,
                "requested_reads_total": 20000,
                "controlled_access": True,
            },
        )

    def test_cache_key_changes_with_sampling_and_tool_identity(self) -> None:
        record = self.configuration_record()
        sequence = MODULE.SequenceFileRecord(
            accession="IGVFFI0001TEST",
            href="/sequence-files/IGVFFI0001TEST/@@download/test.fastq.gz",
            controlled_access=False,
            read_names=["Read1"],
        )
        expected_files = [{"file_id": sequence.accession, "read_id": "Read1"}]

        def key(context: object, n_reads: int) -> str:
            return MODULE.build_report_cache_key(
                context,
                record,
                "rna",
                "0.4.0",
                "0.5.0",
                "https://example.org/spec.yaml",
                expected_files,
                [sequence],
                ["https://example.org/test.fastq.gz"],
                n_reads,
                False,
            )

        context = self.audit_context()
        self.assertEqual(key(context, 10000), key(context, 10000))
        self.assertNotEqual(key(context, 10000), key(context, 100000))

        changed_tool = MODULE.AuditContext(
            run_id=context.run_id,
            sampling_method=context.sampling_method,
            sampling_seed=context.sampling_seed,
            seqcheck=self.tool_identity("seqcheck 0.2.1"),
            seqspec=context.seqspec,
        )
        self.assertNotEqual(key(context, 10000), key(changed_tool, 10000))

        provenance_only = MODULE.AuditContext(
            run_id=context.run_id,
            sampling_method=context.sampling_method,
            sampling_seed=context.sampling_seed,
            seqcheck=MODULE.ToolIdentity(
                command=context.seqcheck.command,
                version=context.seqcheck.version,
                git_commit="different-commit",
                git_dirty=True,
                runtime_source_sha256=context.seqcheck.runtime_source_sha256,
                executable_sha256=context.seqcheck.executable_sha256,
            ),
            seqspec=context.seqspec,
        )
        self.assertEqual(key(context, 10000), key(provenance_only, 10000))

        report = {
            "audit_summary": {
                "audit_schema_version": MODULE.AUDIT_SCHEMA_VERSION,
                "cache_key": key(context, 10000),
            }
        }
        self.assertTrue(MODULE.report_cache_matches(report, key(context, 10000)))
        self.assertFalse(MODULE.report_cache_matches(report, key(context, 100000)))

    def test_runtime_source_hash_ignores_non_runtime_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "src").mkdir()
            (root / "docs").mkdir()
            (root / "Cargo.toml").write_text("[package]\n", encoding="utf-8")
            source = root / "src" / "lib.rs"
            source.write_text("pub fn value() -> u8 { 1 }\n", encoding="utf-8")
            notes = root / "docs" / "notes.md"
            notes.write_text("first\n", encoding="utf-8")

            first = MODULE.path_tree_sha256(root, ["Cargo.toml", "src"])
            notes.write_text("second\n", encoding="utf-8")
            self.assertEqual(
                first, MODULE.path_tree_sha256(root, ["Cargo.toml", "src"])
            )

            source.write_text("pub fn value() -> u8 { 2 }\n", encoding="utf-8")
            self.assertNotEqual(
                first, MODULE.path_tree_sha256(root, ["Cargo.toml", "src"])
            )

    def test_external_tool_identity_hashes_executable_and_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            limits = Path(tmpdir) / "limits.txt"
            limits.write_text("duplication warn = 70\n", encoding="utf-8")
            with patch.object(MODULE.shutil, "which", return_value=sys.executable):
                identity = MODULE.build_external_tool_identity("fastqc", limits)

        self.assertTrue(identity["available"])
        self.assertTrue(identity["version"].startswith("Python"))
        self.assertEqual(identity["configuration"]["mode"], "custom")
        self.assertEqual(identity["configuration"]["path"], str(limits))
        self.assertEqual(len(identity["configuration"]["sha256"]), 64)
        self.assertEqual(len(identity["executable_sha256"]), 64)

    def test_study_manifest_run_id_ignores_runtime_only_fields(self) -> None:
        record = self.configuration_record()
        sequence = MODULE.SequenceFileRecord(
            accession="IGVFFI0001TEST",
            href="/sequence-files/IGVFFI0001TEST/@@download/test.fastq.gz",
            controlled_access=False,
            read_names=["Read1"],
        )
        args = argparse.Namespace(
            api_root="https://api.example.org/",
            portal_root="https://example.org/",
            status="released",
            upload_status="validated",
            public_only=True,
            n_reads=10000,
            fastqc_command="fastqc",
            fastqc_limits=None,
            workers=2,
            force=False,
            auth_profile="igvf",
        )
        external_tool = {
            "command": "fastqc",
            "resolved_path": "/usr/bin/fastqc",
            "available": True,
            "version": "FastQC v0.12.1",
            "configuration": {
                "mode": "packaged_defaults",
                "path": "",
                "sha256": "",
            },
        }

        with tempfile.TemporaryDirectory() as tmpdir, patch.object(
            MODULE, "build_external_tool_identity", return_value=external_tool
        ):
            first = MODULE.write_study_manifest(
                Path(tmpdir) / "first",
                args,
                [record],
                {sequence.accession: sequence},
                self.tool_identity("seqcheck 0.2.0"),
                self.tool_identity("seqspec 0.4.0"),
                None,
                None,
                "2026-07-13T12:00:00+00:00",
            )
            args.workers = 8
            args.force = True
            second = MODULE.write_study_manifest(
                Path(tmpdir) / "second",
                args,
                [record],
                {sequence.accession: sequence},
                self.tool_identity("seqcheck 0.2.0"),
                self.tool_identity("seqspec 0.4.0"),
                None,
                None,
                "2026-07-14T12:00:00+00:00",
            )

            self.assertEqual(first.run_id, second.run_id)
            manifest = json.loads(
                (Path(tmpdir) / "second" / "manifests" / "study.json").read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(manifest["run_id"], first.run_id)
        self.assertEqual(manifest["tools"]["fastqc"]["version"], "FastQC v0.12.1")
        self.assertEqual(manifest["auth"]["requested_profile"], "igvf")
        self.assertFalse(manifest["auth"]["seqcheck_profile_ready"])

    def test_reconciliation_checks_each_report_and_rejects_orphans(self) -> None:
        record = self.configuration_record()
        context = self.audit_context()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_root = Path(tmpdir)
            report_path = output_root / "reports" / record.accession / "rna.json"
            report_path.parent.mkdir(parents=True)
            report = {
                "meta": {"requested_reads": 10},
                "results": [
                    {
                        "check": "coverage",
                        "files": ["IGVFFI0001TEST"],
                        "reads": ["Read1"],
                        "regions": ["barcode"],
                        "expected": [
                            {
                                "id": "e1",
                                "name": "expected_regions",
                                "description": "Expected regions.",
                                "data": {
                                    "kind": "records",
                                    "value": [
                                        {"region_type": ["RGN:partition:cell"]}
                                    ],
                                },
                            }
                        ],
                        "observed": [
                            {
                                "id": "o1",
                                "name": "sampled_count",
                                "description": "Sampled records.",
                                "data": {
                                    "kind": "scalar",
                                    "value": 10,
                                    "unit": "count",
                                },
                            }
                        ],
                        "assessment": [
                            {
                                "type": "pass",
                                "code": "full_coverage",
                                "description": "All records cover the region.",
                            }
                        ],
                    }
                ],
            }
            report = MODULE.annotate_report_summary(
                report,
                configuration_accession=record.accession,
                modality="rna",
                modality_count=1,
                expected_fastq_count=1,
                supplied_fastq_count=1,
                controlled_access=False,
                audit_context=context,
                cache_key="cache-1",
            )
            report_path.write_text(json.dumps(report), encoding="utf-8")
            run, diagnostics, metrics = MODULE.flatten_report(
                report,
                record,
                "rna",
                "0.4.0",
                "0.5.0",
                report_path,
                1,
                1,
                1,
                False,
                "completed",
            )

            reconciled = MODULE.reconcile_outputs(
                output_root,
                context,
                [record],
                [run],
                diagnostics,
                metrics,
                [],
            )
            self.assertTrue(reconciled["valid"])

            missing_metric = MODULE.reconcile_outputs(
                output_root,
                context,
                [record],
                [run],
                diagnostics,
                metrics[:-1],
                [],
            )
            self.assertFalse(missing_metric["checks"]["metric_rows_reconcile"])

            orphan = output_root / "reports" / "orphan.json"
            orphan.write_text("{}", encoding="utf-8")
            with_orphan = MODULE.reconcile_outputs(
                output_root,
                context,
                [record],
                [run],
                diagnostics,
                metrics,
                [],
            )
            self.assertFalse(with_orphan["checks"]["no_orphan_reports"])

    def test_cli_help_documents_reproducibility_options(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), "--help"],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--fastqc-command", completed.stdout)
        self.assertIn("--fastqc-limits", completed.stdout)

    def test_write_lab_summary_aggregates_runs_and_failures(self) -> None:
        runs = [
            MODULE.RunRecord(
                study_run_id="study-1",
                cache_key="cache-1",
                sampling_method="prefix",
                sampling_seed=None,
                configuration_accession="A",
                modality="rna",
                modality_count=1,
                lab="Lab 1",
                submitted_by="User 1",
                award_component="mapping",
                file_set_accession="DS1",
                assay_term="rna",
                preferred_assay_titles="RNA-seq",
                aliases="lab1:a",
                raw_seqspec_version="0.3.0",
                normalized_seqspec_version="0.5.0",
                requested_reads=10000,
                requested_reads_total=10000,
                expected_fastq_count=1,
                supplied_fastq_count=1,
                controlled_access=False,
                report_path="/tmp/a.json",
                run_status="completed",
                pass_count=2,
                warning_count=1,
                error_count=0,
                interpretation_count=1,
            )
        ]
        failures = [
            MODULE.FailureRecord(
                study_run_id="study-1",
                configuration_accession="B",
                modality="rna",
                lab="Lab 1",
                submitted_by="User 2",
                award_component="mapping",
                file_set_accession="DS2",
                assay_term="rna",
                preferred_assay_titles="RNA-seq",
                aliases="lab1:b",
                raw_seqspec_version="0.3.0",
                normalized_seqspec_version="0.5.0",
                stage="seqcheck",
                reason="seqcheck_error",
                message="boom",
            )
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "lab_summary.csv"
            MODULE.write_lab_summary(output, runs, failures)
            text = output.read_text(encoding="utf-8")

        self.assertIn("Lab 1", text)
        self.assertIn("1,1,2,1,0,1", text)


if __name__ == "__main__":
    unittest.main()
