import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "materialize_perturbations.py"
SPEC = importlib.util.spec_from_file_location("materialize_perturbations", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, value: object) -> None:
    MODULE.runtime.write_json(path, value)


def seqspec_value() -> dict[str, object]:
    return {
        "seqspec_version": "0.5.0",
        "assay_id": "synthetic",
        "name": "Synthetic",
        "doi": "https://doi.org/10.0000/synthetic",
        "date": "2026-07-14",
        "description": "Synthetic perturbation fixture",
        "modalities": ["rna"],
        "lib_struct": "",
        "sequence_protocol": "Custom",
        "sequence_kit": "Custom",
        "library_protocol": "Custom",
        "library_kit": "Custom",
        "sequence_spec": [
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
                        "filetype": "fastq",
                        "filesize": 0,
                        "url": "fastqs/synthetic_R1.fastq",
                        "urltype": "local",
                        "md5": "",
                    }
                ],
            }
        ],
        "library_spec": [
            {
                "region_id": "rna",
                "region_type": ["RGN:unknown:unclassified"],
                "name": "Synthetic RNA",
                "sequence_type": "joined",
                "sequence": "AAAANNNNTTXXXXXX",
                "min_len": 16,
                "max_len": 16,
                "onlist": None,
                "regions": [
                    {
                        "region_id": "primer",
                        "region_type": ["RGN:technical:primer"],
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
                        "region_type": ["RGN:partition:cell"],
                        "name": "Barcode",
                        "sequence_type": "onlist",
                        "sequence": "NNNN",
                        "min_len": 4,
                        "max_len": 4,
                        "onlist": {
                            "file_id": "synthetic_barcodes.txt",
                            "filename": "synthetic_barcodes.txt",
                            "filetype": "txt",
                            "filesize": 0,
                            "url": "onlists/synthetic_barcodes.txt",
                            "urltype": "local",
                            "md5": "",
                        },
                        "regions": [],
                    },
                    {
                        "region_id": "linker",
                        "region_type": ["RGN:technical:linker"],
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
                        "region_type": ["RGN:partition:molecule"],
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
                        "region_type": ["RGN:measure:transcript"],
                        "name": "cDNA",
                        "sequence_type": "random",
                        "sequence": "XXXX",
                        "min_len": 4,
                        "max_len": 4,
                        "onlist": None,
                        "regions": [],
                    },
                ],
            }
        ],
    }


def fake_seqspec(path: Path) -> None:
    sequence_spec = seqspec_value()["sequence_spec"]
    leaves = seqspec_value()["library_spec"][0]["regions"]
    script = f"""#!/usr/bin/env python3
import json
import sys
from pathlib import Path

SEQUENCE_SPEC = {sequence_spec!r}
LIBRARY_SPEC = {{"rna": {leaves!r}}}
READ_INDEX = {("synthetic_R1\tBarcode\tRGN:partition:cell\t0\t4\nsynthetic_R1\tLinker\tRGN:technical:linker\t4\t6\nsynthetic_R1\tUMI\tRGN:partition:molecule\t6\t8\nsynthetic_R1\tcDNA\tRGN:measure:transcript\t8\t12\n")!r}
REGION_INDEX = {("rna\tPrimer\tRGN:technical:primer\t0\t4\nrna\tBarcode\tRGN:partition:cell\t4\t8\nrna\tLinker\tRGN:technical:linker\t8\t10\nrna\tUMI\tRGN:partition:molecule\t10\t12\nrna\tcDNA\tRGN:measure:transcript\t12\t16\n")!r}

def regions(value):
    stack = list(value.get("library_spec", []))
    while stack:
        region = stack.pop()
        yield region
        stack.extend(region.get("regions", []))

if sys.argv[1] == "--version":
    print("seqspec 0.5.0-test")
elif sys.argv[1] == "check":
    spec_path = Path(sys.argv[-1])
    try:
        value = json.loads(spec_path.read_text())
    except json.JSONDecodeError:
        raise SystemExit(0)
    for region in regions(value):
        onlist = region.get("onlist")
        if isinstance(onlist, dict) and onlist.get("urltype") == "local":
            resource = spec_path.parent / onlist.get("url", "")
            if not resource.is_file():
                print(f"missing onlist: {{resource}}", file=sys.stderr)
                raise SystemExit(1)
    for read in value.get("sequence_spec", []):
        for file in read.get("files", []):
            if file.get("urltype") == "local":
                resource = spec_path.parent / file.get("url", "")
                if not resource.is_file():
                    print(f"missing FASTQ: {{resource}}", file=sys.stderr)
                    raise SystemExit(1)
elif sys.argv[1] == "info":
    key = sys.argv[sys.argv.index("-k") + 1]
    print(json.dumps(SEQUENCE_SPEC if key == "sequence_spec" else LIBRARY_SPEC))
elif sys.argv[1] == "index":
    print(REGION_INDEX if "-s" in sys.argv else READ_INDEX, end="")
else:
    raise SystemExit(2)
"""
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def fake_yq(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import sys
from pathlib import Path
if sys.argv[1] == "--version":
    print("yq test")
else:
    print(Path(sys.argv[-1]).read_text())
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def fastq_records(count: int) -> list[MODULE.mutation.FastqRecord]:
    return [
        (
            f"@record-{index}\n".encode(),
            b"ACGTTTGGCCCC\n",
            b"+\n",
            b"IIIIIIIIIIII\n",
        )
        for index in range(count)
    ]


def create_inputs(root: Path) -> dict[str, Path]:
    root.mkdir(parents=True)
    seqspec = root / "seqspec"
    yq = root / "yq"
    fake_seqspec(seqspec)
    fake_yq(yq)
    spec = root / "spec.yaml"
    write_json(spec, seqspec_value())
    (root / "onlists").mkdir()
    (root / "onlists" / "synthetic_barcodes.txt").write_text(
        "ACGT\nTGCA\n", encoding="ascii"
    )
    (root / "fastqs").mkdir()
    (root / "fastqs" / "synthetic_R1.fastq").write_text(
        "@source\nACGTTTGGCCCC\n+\nIIIIIIIIIIII\n", encoding="ascii"
    )
    sample = root / "sample" / "synthetic_R1.fastq.gz"
    MODULE.mutation.write_fastq(sample, fastq_records(20))

    protocol = root / "perturbations.json"
    value = MODULE.runtime.load_json(
        ROOT / "experiments" / "paper" / "protocol" / "perturbations.json"
    )
    value["stochastic_conditions"]["event_fractions"] = [0.1]
    value["stochastic_conditions"]["seeds"] = [7]
    write_json(protocol, value)
    selected_case = {
        "family_id": "rna_family",
        "configuration_accession": "CONFIG1",
        "modality": "rna",
        "selection_role": "primary",
        "read_id": "synthetic_R1",
        "fastq_accession": "synthetic_R1.fastq",
        "fastq_url": str(root / "fastqs" / "synthetic_R1.fastq"),
        "declared_compressed_bytes": sample.stat().st_size,
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
        key: item for key, item in selected_case.items() if key != "spec_path"
    }
    stable_selection = {
        "schema_version": "0.1.0",
        "freeze_id": "freeze-test",
        "cohort_sha256": "cohort-test",
        "sampling_protocol_sha256": "sampling-test",
        "selector": {"version": "test", "sha256": "selector-test"},
        "seqspec": {"version": "seqspec 0.5.0-test", "sha256": "selection-test"},
        "cases": [stable_case],
    }
    selection_id = MODULE.runtime.sha256_json(stable_selection)[:16]
    cases_path = root / "cases.json"
    case = {
        **selected_case,
        "selection_id": selection_id,
        "case_id": "config1--r1",
    }
    write_json(
        cases_path,
        {
            "schema_version": "0.1.0",
            "selection_id": selection_id,
            "freeze_id": "freeze-test",
            "cases": [case],
        },
    )
    selection = root / "selection.json"
    write_json(
        selection,
        {
            **stable_selection,
            "selection_id": selection_id,
            "valid": True,
            "counts": {"fastq_cases": 1},
            "outputs": {"cases_json": MODULE.runtime.file_identity(cases_path)},
        },
    )
    inventory = MODULE.inventory_builder.build_perturbation_cases(
        case_selection_manifest_path=selection,
        perturbation_protocol_path=protocol,
        seqspec_bin=seqspec,
        output_root=root / "inventory",
        timeout_seconds=10,
    )

    matrix = root / "matrix.json"
    matrix_id = "matrix-test"
    write_json(
        matrix,
        {
            "matrix_schema_version": "0.1.0",
            "matrix_id": matrix_id,
            "conditions": [
                {
                    "condition_id": "base-sample",
                    "method": "reservoir",
                    "requested_records_per_fastq": 20,
                    "seed": 17,
                    "outputs": [
                        {
                            "output_path": str(sample),
                            "output_sha256": MODULE.runtime.file_sha256(sample),
                            "records_selected": 20,
                        }
                    ],
                }
            ],
        },
    )
    study_identity = {
        "schema_version": "0.1.0",
        "selection_id": selection_id,
        "sampling_protocol_sha256": "sampling-test",
        "seqcheck": {"version": "seqcheck test", "sha256": "seqcheck-test"},
        "sampler": {"version": "test", "sha256": "sampler-test", "size_bytes": 1},
        "sampler_core": {"version": "test", "sha256": "core-test", "size_bytes": 1},
        "runtime": {"version": "test", "sha256": "runtime-test", "size_bytes": 1},
        "runner": {"version": "test", "sha256": "runner-test", "size_bytes": 1},
    }
    study_run_id = MODULE.runtime.sha256_json(study_identity)[:16]
    study = root / "study.json"
    write_json(
        study,
        {
            **study_identity,
            "study_run_id": study_run_id,
            "valid": True,
            "inputs": {
                "case_selection_manifest": MODULE.runtime.file_identity(selection)
            },
            "matrices": [
                {
                    "case_id": case["case_id"],
                    "matrix_id": matrix_id,
                    "manifest": MODULE.runtime.file_identity(matrix),
                }
            ],
        },
    )
    stable_policy = {
        "schema_version": "0.1.0",
        "study_run_id": study_run_id,
        "analysis_run_id": "analysis-test",
        "analysis_protocol_sha256": "analysis-test",
        "frozen": True,
        "scientific_targets_met": True,
        "default": {"records_per_fastq": 20, "sampling_method": "reservoir"},
        "escalation_records_per_fastq": 20,
        "endpoint_decisions": [],
        "memory": {},
    }
    policy = root / "policy.json"
    write_json(
        policy,
        {
            **stable_policy,
            "policy_id": MODULE.runtime.sha256_json(stable_policy)[:16],
            "created_at": "2026-07-14T00:00:00+00:00",
        },
    )
    return {
        "inventory": Path(inventory["manifest_path"]),
        "study": study,
        "policy": policy,
        "protocol": protocol,
        "seqspec": seqspec,
        "yq": yq,
    }


class MaterializePerturbationsTests(unittest.TestCase):
    def test_materializer_emits_complete_deterministic_conditions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "inputs")
            first = MODULE.materialize_perturbations(
                inventory_manifest_path=inputs["inventory"],
                study_manifest_path=inputs["study"],
                sampling_policy_path=inputs["policy"],
                perturbation_protocol_path=inputs["protocol"],
                seqspec_bin=inputs["seqspec"],
                yq_bin=inputs["yq"],
                output_root=root / "first",
                timeout_seconds=10,
            )
            second = MODULE.materialize_perturbations(
                inventory_manifest_path=inputs["inventory"],
                study_manifest_path=inputs["study"],
                sampling_policy_path=inputs["policy"],
                perturbation_protocol_path=inputs["protocol"],
                seqspec_bin=inputs["seqspec"],
                yq_bin=inputs["yq"],
                output_root=root / "second",
                timeout_seconds=10,
            )
            payload = MODULE.runtime.load_json(
                Path(first["outputs"]["conditions"]["path"])
            )
            validation = MODULE.runtime.load_json(
                Path(first["outputs"]["validation"]["path"])
            )
            second_condition_ids = [
                value["condition_id"]
                for value in MODULE.runtime.load_json(
                    Path(second["outputs"]["conditions"]["path"])
                )["conditions"]
            ]

        conditions = payload["conditions"]
        by_operator = {}
        for condition in conditions:
            by_operator.setdefault(condition["operator_id"], []).append(condition)
        self.assertTrue(validation["valid"])
        self.assertEqual(validation["counts"]["conditions"], 13)
        self.assertEqual(validation["counts"]["skips"], 5)
        self.assertEqual(len(by_operator["CLEAN"]), 1)
        self.assertEqual(len(by_operator["S08"]), 2)
        self.assertEqual(len(by_operator["D04"]), 2)
        self.assertEqual(by_operator["S10"][0]["observed_seqspec_check"], "failure")
        self.assertTrue(
            all(
                condition["mutated_record_count"] == 2
                for condition in by_operator["D01"]
            )
        )
        self.assertEqual(first["materialization_id"], second["materialization_id"])
        self.assertEqual(
            [value["condition_id"] for value in conditions],
            second_condition_ids,
        )

    def test_materializer_rejects_unfrozen_policy(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "inputs")
            policy = MODULE.runtime.load_json(inputs["policy"])
            stable = {
                key: value
                for key, value in policy.items()
                if key not in {"policy_id", "created_at"}
            }
            stable["frozen"] = False
            write_json(
                inputs["policy"],
                {
                    **stable,
                    "policy_id": MODULE.runtime.sha256_json(stable)[:16],
                    "created_at": "2026-07-14T00:00:00+00:00",
                },
            )

            with self.assertRaisesRegex(ValueError, "policy is not frozen"):
                MODULE.materialize_perturbations(
                    inventory_manifest_path=inputs["inventory"],
                    study_manifest_path=inputs["study"],
                    sampling_policy_path=inputs["policy"],
                    perturbation_protocol_path=inputs["protocol"],
                    seqspec_bin=inputs["seqspec"],
                    yq_bin=inputs["yq"],
                    output_root=root / "output",
                    timeout_seconds=10,
                )

    def test_cli_writes_a_valid_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "inputs")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--perturbation-inventory",
                    str(inputs["inventory"]),
                    "--sampling-study",
                    str(inputs["study"]),
                    "--sampling-policy",
                    str(inputs["policy"]),
                    "--perturbation-protocol",
                    str(inputs["protocol"]),
                    "--seqspec-bin",
                    str(inputs["seqspec"]),
                    "--yq-bin",
                    str(inputs["yq"]),
                    "--output-root",
                    str(root / "output"),
                    "--timeout-seconds",
                    "10",
                ],
                capture_output=True,
                text=True,
            )
            manifest = MODULE.runtime.load_json(Path(completed.stdout.strip()))

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(manifest["valid"])


if __name__ == "__main__":
    unittest.main()
