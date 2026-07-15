import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "build_perturbation_cases.py"
PROTOCOL_PATH = ROOT / "experiments" / "paper" / "protocol" / "perturbations.json"
SPEC = importlib.util.spec_from_file_location("build_perturbation_cases", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fake_seqspec(path: Path) -> None:
    sequence_spec = [
        {
            "read_id": "synthetic_R1",
            "name": "Synthetic Read 1",
            "modality": "rna",
            "primer_id": "primer",
            "min_len": 12,
            "max_len": 12,
            "strand": "pos",
            "files": [
                {
                    "file_id": "synthetic_R1.fastq",
                    "filename": "synthetic_R1.fastq",
                    "url": "fastqs/synthetic_R1.fastq",
                }
            ],
        }
    ]
    regions = [
        {
            "region_id": "primer",
            "region_type": "truseq_read1",
            "name": "Primer",
            "sequence_type": "fixed",
            "sequence": "AAAA",
            "min_len": 4,
            "max_len": 4,
            "onlist": None,
            "regions": [],
        },
        {
            "region_id": "barcode",
            "region_type": "barcode",
            "name": "Barcode",
            "sequence_type": "onlist",
            "sequence": "NNNN",
            "min_len": 4,
            "max_len": 4,
            "onlist": {
                "filename": "synthetic_barcodes.txt",
                "url": "onlists/synthetic_barcodes.txt",
                "urltype": "local",
            },
            "regions": [],
        },
        {
            "region_id": "linker",
            "region_type": "linker",
            "name": "Linker",
            "sequence_type": "fixed",
            "sequence": "TT",
            "min_len": 2,
            "max_len": 2,
            "onlist": None,
            "regions": [],
        },
        {
            "region_id": "umi",
            "region_type": "umi",
            "name": "UMI",
            "sequence_type": "random",
            "sequence": "XX",
            "min_len": 2,
            "max_len": 2,
            "onlist": None,
            "regions": [],
        },
        {
            "region_id": "cdna",
            "region_type": "cdna",
            "name": "cDNA",
            "sequence_type": "random",
            "sequence": "XXXX",
            "min_len": 4,
            "max_len": 4,
            "onlist": None,
            "regions": [],
        },
    ]
    script = (
        f"""#!/usr/bin/env python3
import json
import sys

SEQUENCE_SPEC = {sequence_spec!r}
LIBRARY_SPEC = {{"rna": {regions!r}}}
READ_INDEX = """
        + repr(
            "synthetic_R1\tBarcode\tbarcode\t0\t4\n"
            "synthetic_R1\tLinker\tlinker\t4\t6\n"
            "synthetic_R1\tUMI\tumi\t6\t8\n"
            "synthetic_R1\tcDNA\tcdna\t8\t12\n"
        )
        + """
REGION_INDEX = """
        + repr(
            "rna\tPrimer\ttruseq_read1\t0\t4\n"
            "rna\tBarcode\tbarcode\t4\t8\n"
            "rna\tLinker\tlinker\t8\t10\n"
            "rna\tUMI\tumi\t10\t12\n"
            "rna\tcDNA\tcdna\t12\t16\n"
        )
        + """

if sys.argv[1] == "--version":
    print("seqspec 0.5.0-test")
elif sys.argv[1] == "check":
    print("valid")
elif sys.argv[1] == "info":
    key = sys.argv[sys.argv.index("-k") + 1]
    print(json.dumps(SEQUENCE_SPEC if key == "sequence_spec" else LIBRARY_SPEC))
elif sys.argv[1] == "index":
    print(REGION_INDEX if "-s" in sys.argv else READ_INDEX, end="")
else:
    raise SystemExit(2)
"""
    )
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def create_inputs(root: Path) -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    seqspec = root / "seqspec"
    fake_seqspec(seqspec)
    spec = root / "spec.yaml"
    spec.write_text("seqspec_version: 0.5.0\n", encoding="utf-8")
    onlist = root / "onlists" / "synthetic_barcodes.txt"
    onlist.parent.mkdir()
    onlist.write_text("ACGT\nTGCA\n", encoding="ascii")
    fastq = root / "synthetic_R1.fastq"
    fastq.write_text("@read-1\nACGTTTXXXXXX\n+\n!!!!!!!!!!!!\n", encoding="ascii")

    selected_case = {
        "family_id": "rna_family",
        "configuration_accession": "CONFIG1",
        "modality": "rna",
        "selection_role": "primary",
        "read_id": "synthetic_R1",
        "fastq_accession": "FASTQ1",
        "fastq_url": str(fastq),
        "declared_compressed_bytes": fastq.stat().st_size,
        "indexed_bases": 12,
        "measure_bases": 4,
        "partition_bases": 4,
        "technical_bases": 2,
        "role_classes": "measure;partition;technical",
        "ontology_terms": "RGN:measure:transcript;RGN:partition:cell",
        "spec_path": str(spec),
        "spec_sha256": MODULE.runtime.file_sha256(spec),
    }
    stable_case = {
        key: value for key, value in selected_case.items() if key != "spec_path"
    }
    stable_selection = {
        "schema_version": "0.1.0",
        "freeze_id": "freeze-test",
        "cohort_sha256": "cohort-sha256-test",
        "sampling_protocol_sha256": "sampling-protocol-sha256-test",
        "selector": {"version": "0.1.0", "sha256": "selector-sha256-test"},
        "seqspec": {"version": "seqspec 0.5.0-test", "sha256": "seqspec-test"},
        "cases": [stable_case],
    }
    selection_id = MODULE.runtime.sha256_json(stable_selection)[:16]
    cases_path = root / "cases.json"
    write_json(
        cases_path,
        {
            "schema_version": "0.1.0",
            "selection_id": selection_id,
            "freeze_id": "freeze-test",
            "cases": [
                {
                    **selected_case,
                    "selection_id": selection_id,
                    "case_id": "config1--fastq1",
                }
            ],
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
    return {
        "seqspec": seqspec,
        "spec": spec,
        "onlist": onlist,
        "fastq": fastq,
        "cases": cases_path,
        "selection": selection_path,
    }


class BuildPerturbationCasesTests(unittest.TestCase):
    def test_ambiguous_input_name_matches_at_same_priority(self) -> None:
        reads = {
            "read_a": {
                "files": [
                    {
                        "file_id": "shared.fastq",
                        "filename": "a.fastq",
                        "url": "a.fastq",
                    }
                ]
            },
            "read_b": {
                "files": [
                    {
                        "file_id": "shared.fastq",
                        "filename": "b.fastq",
                        "url": "b.fastq",
                    }
                ]
            },
        }

        self.assertEqual(MODULE.find_ambiguous_input_name(reads), "shared.fastq")

    def test_boundary_shift_respects_declared_minimum_length(self) -> None:
        profile = {
            "case": {"selection_role": "primary", "fastq_accession": "FASTQ1"},
            "projected": [
                {
                    "region_id": "variable",
                    "sequence_type": "random",
                    "min_len": 1,
                    "start": 0,
                    "stop": 10,
                },
                {
                    "region_id": "fixed",
                    "sequence_type": "fixed",
                    "min_len": 1,
                    "start": 10,
                    "stop": 11,
                },
            ],
        }

        self.assertIsNone(MODULE.select_boundary([profile]))

    def test_builder_inventories_applicable_targets_deterministically(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "inputs")
            first = MODULE.build_perturbation_cases(
                case_selection_manifest_path=inputs["selection"],
                perturbation_protocol_path=PROTOCOL_PATH,
                seqspec_bin=inputs["seqspec"],
                output_root=root / "first",
                timeout_seconds=10,
            )
            second = MODULE.build_perturbation_cases(
                case_selection_manifest_path=inputs["selection"],
                perturbation_protocol_path=PROTOCOL_PATH,
                seqspec_bin=inputs["seqspec"],
                output_root=root / "second",
                timeout_seconds=10,
            )
            onlist_sha256 = MODULE.runtime.file_sha256(inputs["onlist"])
            validation = json.loads(
                Path(first["outputs"]["validation"]["path"]).read_text()
            )
            records = json.loads(Path(first["outputs"]["cases"]["path"]).read_text())[
                "records"
            ]
            rows = read_csv(first["outputs"]["applicability"]["path"])

        by_operator = {record["operator_id"]: record for record in records}
        self.assertTrue(first["valid"])
        self.assertTrue(validation["valid"])
        self.assertEqual(validation["counts"]["records"], 14)
        self.assertEqual(validation["counts"]["applicable"], 10)
        self.assertEqual(len(rows), 14)
        self.assertFalse(by_operator["S02"]["applicable"])
        self.assertFalse(by_operator["S03"]["applicable"])
        self.assertFalse(by_operator["S06"]["applicable"])
        self.assertFalse(by_operator["S07"]["applicable"])
        self.assertEqual(by_operator["S01"]["variants"], ["unexpected_name"])
        self.assertEqual(by_operator["S05"]["variants"], ["shift_1"])
        self.assertEqual(
            by_operator["S05"]["target"]["region_ids"], ["barcode", "linker"]
        )
        self.assertEqual(by_operator["S09"]["target"]["region_ids"], ["barcode"])
        self.assertEqual(by_operator["D01"]["target"]["region_ids"], ["cdna"])
        self.assertEqual(by_operator["D04"]["target"]["region_ids"], ["primer"])
        self.assertEqual(
            by_operator["S09"]["target"]["regions"][0]["resource"]["sha256"],
            onlist_sha256,
        )
        self.assertEqual(first["inventory_id"], second["inventory_id"])
        self.assertEqual(
            first["outputs"]["applicability"]["sha256"],
            second["outputs"]["applicability"]["sha256"],
        )
        self.assertEqual(
            first["outputs"]["cases"]["sha256"],
            second["outputs"]["cases"]["sha256"],
        )

    def test_builder_rejects_changed_content_addressed_cases(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "inputs")
            payload = json.loads(inputs["cases"].read_text())
            payload["cases"][0]["fastq_accession"] = "CHANGED"
            write_json(inputs["cases"], payload)

            with self.assertRaisesRegex(ValueError, "cases JSON hash changed"):
                MODULE.build_perturbation_cases(
                    case_selection_manifest_path=inputs["selection"],
                    perturbation_protocol_path=PROTOCOL_PATH,
                    seqspec_bin=inputs["seqspec"],
                    output_root=root / "output",
                    timeout_seconds=10,
                )

    def test_cli_writes_manifest_and_refuses_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "inputs")
            output = root / "output"
            command = [
                sys.executable,
                str(SCRIPT_PATH),
                "--case-selection-manifest",
                str(inputs["selection"]),
                "--perturbation-protocol",
                str(PROTOCOL_PATH),
                "--seqspec-bin",
                str(inputs["seqspec"]),
                "--output-root",
                str(output),
                "--timeout-seconds",
                "10",
            ]
            completed = subprocess.run(command, capture_output=True, text=True)
            repeated = subprocess.run(command, capture_output=True, text=True)

            manifest_path = Path(completed.stdout.strip())
            manifest = json.loads(manifest_path.read_text())

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(manifest["valid"])
        self.assertEqual(repeated.returncode, 1)
        self.assertIn("refusing to overwrite", repeated.stderr)


if __name__ == "__main__":
    unittest.main()
