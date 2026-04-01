import importlib.util
import os
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
    def test_parse_seqspec_version_output(self) -> None:
        output = "seqspec version: 0.4.2\nseqspec file version: 0.3.0\n"
        self.assertEqual(MODULE.parse_seqspec_version_output(output), "0.3.0")

    def test_normalize_seqspec_version_maps_legacy_versions(self) -> None:
        self.assertEqual(MODULE.normalize_seqspec_version("0.3.0"), "0.4.0")
        self.assertEqual(MODULE.normalize_seqspec_version("0.4.0"), "0.4.0")
        self.assertEqual(MODULE.normalize_seqspec_version("0.5.0"), "0.5.0")

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
        record = MODULE.ConfigurationRecord(
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
        report = {
            "meta": {"requested_reads": 10000},
            "results": [
                {
                    "check": "input_check",
                    "files": ["IGVFFI0001TEST"],
                    "reads": ["Read1"],
                    "regions": [],
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

        run, diagnostics = MODULE.flatten_report(
            report=report,
            record=record,
            modality="rna",
            raw_seqspec_version="0.3.0",
            normalized_seqspec_version="0.4.0",
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
        self.assertEqual(len(diagnostics), 3)
        self.assertEqual(diagnostics[1].assessment_code, "length_out_of_range")

    def test_annotate_report_summary_adds_counts(self) -> None:
        report = {"meta": {"requested_reads": 10000}, "results": []}
        annotated = MODULE.annotate_report_summary(
            report,
            configuration_accession="IGVFFI1234TEST",
            modality="rna",
            modality_count=3,
            expected_fastq_count=2,
            supplied_fastq_count=2,
            controlled_access=True,
        )

        self.assertEqual(
            annotated["audit_summary"],
            {
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

    def test_write_lab_summary_aggregates_runs_and_failures(self) -> None:
        runs = [
            MODULE.RunRecord(
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
                normalized_seqspec_version="0.4.0",
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
                normalized_seqspec_version="0.4.0",
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
