import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "build_sampling_cases.py"
)
SPEC = importlib.util.spec_from_file_location("build_sampling_cases", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def fake_seqspec(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import json
import sys

if sys.argv[1] == "--version":
    print("seqspec 0.5.0-test")
    raise SystemExit(0)

with open(sys.argv[-1], encoding="utf-8") as handle:
    payload = json.load(handle)
if sys.argv[1] == "info":
    key = sys.argv[sys.argv.index("-k") + 1]
    print(json.dumps(payload[key]))
elif sys.argv[1] == "index":
    key = "region_index" if "--selector" in sys.argv else "index"
    print(payload[key], end="")
else:
    raise SystemExit(2)
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def fixture_spec() -> dict[str, object]:
    return {
        "sequence_spec": [
            {
                "read_id": "read-measure",
                "modality": "rna",
                "files": [
                    {
                        "file_id": "FASTQ_MEASURE",
                        "filesize": 100,
                    }
                ],
            },
            {
                "read_id": "read-complement",
                "modality": "rna",
                "files": [
                    {
                        "file_id": "FASTQ_COMPLEMENT",
                        "filesize": 80,
                    }
                ],
            },
        ],
        "library_spec": {
            "rna": [
                {
                    "region_id": "rna",
                    "region_type": ["RGN:unknown:unclassified"],
                    "regions": [
                        {
                            "region_id": "transcript",
                            "name": "Transcript",
                            "region_type": ["RGN:measure:transcript"],
                            "regions": [],
                        },
                        {
                            "region_id": "barcode",
                            "name": "Barcode",
                            "region_type": ["RGN:partition:cell"],
                            "regions": [],
                        },
                        {
                            "region_id": "linker",
                            "name": "Linker",
                            "region_type": ["RGN:technical:linker"],
                            "regions": [],
                        },
                    ],
                }
            ]
        },
        "index": (
            "read-measure\tTranscript\tcdna\t0\t90\n"
            "read-complement\tTranscript\tcdna\t0\t30\n"
            "read-complement\tBarcode\tbarcode\t30\t46\n"
            "read-complement\tLinker\tlinker\t46\t50\n"
        ),
        "region_index": (
            "rna\tTranscript\tcdna\t0\t90\n"
            "rna\tBarcode\tbarcode\t90\t106\n"
            "rna\tLinker\tlinker\t106\t110\n"
        ),
    }


def create_inputs(root: Path) -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    seqspec = root / "seqspec"
    fake_seqspec(seqspec)
    spec = root / "spec.json"
    write_json(spec, fixture_spec())
    cohort = root / "cohort.csv"
    write_csv(
        cohort,
        [
            {
                "split": "calibration",
                "final_family": "rna_family",
                "configuration_accession": "CONFIG1",
                "modalities": "rna",
                "effective_spec_path": str(spec),
                "effective_spec_sha256": MODULE.file_sha256(spec),
                "fastq_accessions": "FASTQ_MEASURE;FASTQ_COMPLEMENT",
                "fastq_urls": "https://example.test/measure.fastq.gz;https://example.test/complement.fastq.gz",
                "effective_expected_fastq_accessions": (
                    "FASTQ_COMPLEMENT;FASTQ_MEASURE"
                ),
            }
        ],
    )
    cohort_manifest = root / "cohort_manifest.json"
    write_json(
        cohort_manifest,
        {
            "schema_version": "0.1.0",
            "freeze_id": "freeze-test",
            "frozen": True,
            "counts": {"calibration": 1},
            "outputs": {"cohort": MODULE.file_identity(cohort)},
        },
    )
    protocol = root / "protocol.json"
    write_json(
        protocol,
        {
            "schema_version": "0.1.0",
            "cohort_limits": {"configurations": 1, "fastqs": 2},
            "case_selection": {
                "max_fastqs_per_configuration": 2,
                "primary_read_rule": "largest_measurement_span",
                "secondary_read_rule": ("maximum_additional_ontology_role_coverage"),
                "require_positive_declared_compressed_size": True,
            },
        },
    )
    return {
        "seqspec": seqspec,
        "spec": spec,
        "cohort": cohort,
        "cohort_manifest": cohort_manifest,
        "protocol": protocol,
    }


class BuildSamplingCasesTests(unittest.TestCase):
    def test_selects_measurement_primary_and_ontology_complement(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root)
            first = MODULE.build_sampling_cases(
                cohort_manifest_path=inputs["cohort_manifest"],
                protocol_path=inputs["protocol"],
                seqspec_bin=inputs["seqspec"],
                output_root=root / "first",
                timeout_seconds=10,
            )
            second = MODULE.build_sampling_cases(
                cohort_manifest_path=inputs["cohort_manifest"],
                protocol_path=inputs["protocol"],
                seqspec_bin=inputs["seqspec"],
                output_root=root / "second",
                timeout_seconds=10,
            )
            cases = json.loads(
                Path(first["outputs"]["cases_json"]["path"]).read_text(encoding="utf-8")
            )["cases"]

        self.assertEqual(first["selection_id"], second["selection_id"])
        self.assertEqual(first["counts"]["configurations"], 1)
        self.assertEqual(first["counts"]["fastq_cases"], 2)
        self.assertEqual(first["counts"]["configuration_families"], {"rna_family": 1})
        self.assertEqual(first["counts"]["fastq_case_families"], {"rna_family": 2})
        self.assertEqual(
            first["declared_transfer"],
            {
                "selected_source_bytes": 180,
                "minimum_source_traversals_per_fastq": 2,
                "minimum_remote_bytes": 360,
                "explanation": (
                    "One traversal creates bounded samples and one traversal "
                    "creates the complete-stream seqcheck reference."
                ),
            },
        )
        primary, secondary = cases
        self.assertEqual(primary["selection_role"], "primary")
        self.assertEqual(primary["fastq_accession"], "FASTQ_MEASURE")
        self.assertEqual(primary["measure_bases"], 90)
        self.assertEqual(secondary["selection_role"], "secondary")
        self.assertEqual(secondary["fastq_accession"], "FASTQ_COMPLEMENT")
        self.assertEqual(secondary["partition_bases"], 16)
        self.assertEqual(secondary["technical_bases"], 4)

    def test_redundant_secondary_is_not_selected(self) -> None:
        base = {
            "measure_bases": 10,
            "partition_bases": 0,
            "technical_bases": 0,
            "indexed_bases": 10,
            "declared_compressed_bytes": 20,
            "role_class_values": ["measure"],
            "ontology_term_values": ["RGN:measure:transcript"],
        }
        selected = MODULE.select_candidates(
            [
                {**base, "fastq_accession": "A"},
                {
                    **base,
                    "fastq_accession": "B",
                    "ontology_term_values": [
                        "RGN:measure:transcript",
                        "RGN:unknown:unclassified",
                    ],
                },
            ],
            2,
        )
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["fastq_accession"], "A")

    def test_sequence_file_accession_prefers_url_and_normalizes_suffix(self) -> None:
        self.assertEqual(
            MODULE.sequence_file_accession(
                {
                    "file_id": "different.fastq.gz",
                    "url": (
                        "https://api.data.igvf.org/sequence-files/"
                        "IGVFFI1234TEST/@@download/reads.fastq.gz"
                    ),
                }
            ),
            "IGVFFI1234TEST",
        )
        self.assertEqual(
            MODULE.sequence_file_accession({"file_id": "reads.fastq.gz"}),
            "reads",
        )

    def test_tampered_or_unfrozen_cohort_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root)
            inputs["cohort"].write_text("changed\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash changed"):
                MODULE.build_sampling_cases(
                    cohort_manifest_path=inputs["cohort_manifest"],
                    protocol_path=inputs["protocol"],
                    seqspec_bin=inputs["seqspec"],
                    output_root=root / "tampered",
                    timeout_seconds=10,
                )

            inputs = create_inputs(root / "unfrozen")
            manifest = json.loads(inputs["cohort_manifest"].read_text())
            manifest["frozen"] = False
            write_json(inputs["cohort_manifest"], manifest)
            with self.assertRaisesRegex(ValueError, "not frozen"):
                MODULE.build_sampling_cases(
                    cohort_manifest_path=inputs["cohort_manifest"],
                    protocol_path=inputs["protocol"],
                    seqspec_bin=inputs["seqspec"],
                    output_root=root / "unfrozen-output",
                    timeout_seconds=10,
                )

    def test_cli_writes_manifest_and_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root)
            command = [
                sys.executable,
                str(SCRIPT_PATH),
                "--cohort-manifest",
                str(inputs["cohort_manifest"]),
                "--sampling-protocol",
                str(inputs["protocol"]),
                "--seqspec-bin",
                str(inputs["seqspec"]),
                "--output-root",
                str(root / "output"),
            ]
            completed = subprocess.run(
                command, check=False, capture_output=True, text=True
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest_path = Path(completed.stdout.strip())
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            repeated = subprocess.run(
                command, check=False, capture_output=True, text=True
            )

        self.assertTrue(manifest["valid"])
        self.assertEqual(manifest["counts"]["fastq_cases"], 2)
        self.assertNotEqual(repeated.returncode, 0)
        self.assertIn("refusing to overwrite", repeated.stderr)


if __name__ == "__main__":
    unittest.main()
