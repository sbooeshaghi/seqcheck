import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "perturbation_runtime.py"
)
SPEC = importlib.util.spec_from_file_location("perturbation_runtime", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def assay() -> dict[str, object]:
    return {
        "seqspec_version": "0.5.0",
        "sequence_spec": [
            {
                "read_id": "read1",
                "primer_id": "primer",
                "strand": "pos",
                "min_len": 8,
                "max_len": 8,
                "files": [],
            }
        ],
        "library_spec": [
            {
                "region_id": "rna",
                "sequence_type": "joined",
                "sequence": "AAAACCGG",
                "min_len": 8,
                "max_len": 8,
                "regions": [
                    {
                        "region_id": "primer",
                        "name": "Primer",
                        "sequence_type": "fixed",
                        "sequence": "AAAA",
                        "min_len": 4,
                        "max_len": 4,
                        "regions": [],
                    },
                    {
                        "region_id": "left",
                        "name": "Left",
                        "sequence_type": "fixed",
                        "sequence": "CCGG",
                        "min_len": 4,
                        "max_len": 4,
                        "regions": [],
                    },
                    {
                        "region_id": "right",
                        "name": "Right",
                        "sequence_type": "random",
                        "sequence": "XX",
                        "min_len": 2,
                        "max_len": 2,
                        "regions": [],
                    },
                ],
            }
        ],
    }


def target(*, regions: list[dict[str, object]] | None = None) -> dict[str, object]:
    return {
        "read_ids": ["read1"],
        "region_ids": [region["region_id"] for region in regions or []],
        "regions": regions or [],
        "details": {},
    }


def fastq_record(identifier: str, sequence: str) -> MODULE.FastqRecord:
    quality = "I" * len(sequence)
    return (
        f"@{identifier}\n".encode(),
        f"{sequence}\n".encode(),
        b"+\n",
        f"{quality}\n".encode(),
    )


class PerturbationRuntimeTests(unittest.TestCase):
    def test_mutate_spec_preserves_boundary_total_and_changes_one_contract(
        self,
    ) -> None:
        original = assay()
        boundary_target = target(
            regions=[
                {"region_id": "left"},
                {"region_id": "right"},
            ]
        )
        boundary_target["details"] = {
            "direction": "left_to_right",
            "deltas": [1],
        }
        shifted, change = MODULE.mutate_spec(
            original,
            operator_id="S05",
            variant="shift_1",
            target=boundary_target,
        )
        left = MODULE.find_region(shifted, "left")
        right = MODULE.find_region(shifted, "right")

        self.assertEqual(left["sequence"], "CCG")
        self.assertEqual(right["sequence"], "GXX")
        self.assertEqual(left["max_len"] + right["max_len"], 6)
        self.assertEqual(change["delta"], 1)
        self.assertEqual(MODULE.find_region(original, "left")["sequence"], "CCGG")

    def test_mutate_spec_handles_length_strand_primer_and_fixed_sequence(self) -> None:
        original = assay()
        length, _ = MODULE.mutate_spec(
            original,
            operator_id="S04",
            variant="exclude_observed_length",
            target=target(),
            observed_lengths=[8, 9],
        )
        strand_target = target(regions=[{"region_id": "left"}])
        strand, _ = MODULE.mutate_spec(
            original,
            operator_id="S06",
            variant="pos_to_neg",
            target=strand_target,
        )
        primer_target = target(regions=[{"region_id": "left"}])
        primer_target["details"] = {"replacement_primer_id": "left"}
        primer, _ = MODULE.mutate_spec(
            original,
            operator_id="S07",
            variant="alternate_fixed_primer",
            target=primer_target,
        )
        fixed_target = target(regions=[{"region_id": "left"}])
        fixed, _ = MODULE.mutate_spec(
            original,
            operator_id="S08",
            variant="deterministic_replacement",
            target=fixed_target,
        )

        self.assertEqual(MODULE.find_read(length, "read1")["min_len"], 7)
        self.assertEqual(MODULE.find_read(strand, "read1")["strand"], "neg")
        self.assertEqual(MODULE.find_read(primer, "read1")["primer_id"], "left")
        self.assertEqual(MODULE.find_region(fixed, "left")["sequence"], "GGTT")
        self.assertEqual(MODULE.find_region(fixed, "rna")["sequence"], "AAAAGGTTXX")
        with self.assertRaisesRegex(ValueError, "did not change"):
            MODULE.mutate_spec(
                original,
                operator_id="S08",
                variant="reverse_complement",
                target=fixed_target,
            )

    def test_record_selection_and_mutation_are_deterministic(self) -> None:
        records = [fastq_record(f"record-{index}", "ACGT") for index in range(10)]
        first = MODULE.select_record_indices(records, event_fraction=0.2, seed=17)
        second = MODULE.select_record_indices(records, event_fraction=0.2, seed=17)
        changed = MODULE.mutate_fastq_records(
            records,
            first,
            lambda sequence: sequence[:2],
        )

        self.assertEqual(first, second)
        self.assertEqual(len(first), 2)
        self.assertEqual(sum(len(record[1].strip()) == 2 for record in changed), 2)
        self.assertTrue(
            all(len(record[1].strip()) == len(record[3].strip()) for record in changed)
        )
        self.assertEqual(MODULE.target_record_count(10, 0.01), 0)

    def test_read_mutations_preserve_length_except_truncation(self) -> None:
        region = {
            "region_id": "target",
            "sequence": "AAAA",
            "start": 2,
            "stop": 6,
        }
        read_target = {
            "regions": [region],
            "details": {"target_position": 2, "truncate_to_bases": 5},
        }
        sequence = b"GGAAAATT"
        d01, d01_eligible = MODULE.mutation_for_read_operator(
            "D01", "truncate_before_expected_end", read_target
        )
        d02, d02_eligible = MODULE.mutation_for_read_operator(
            "D02", "single_base_substitution", read_target
        )
        d03, d03_eligible = MODULE.mutation_for_read_operator(
            "D03",
            "deterministic_offlist_substitution",
            read_target,
            onlist_entries={"AAAA"},
        )
        d04, d04_eligible = MODULE.mutation_for_read_operator(
            "D04", "introduce_by_substitution", read_target
        )

        self.assertTrue(d01_eligible(sequence))
        self.assertEqual(len(d01(sequence)), 5)
        self.assertTrue(d02_eligible(sequence))
        self.assertEqual(len(d02(sequence)), len(sequence))
        self.assertTrue(d03_eligible(sequence))
        self.assertEqual(d03(sequence)[2:6], b"AAAC")
        self.assertTrue(d04_eligible(sequence))
        self.assertEqual(len(d04(sequence)), len(sequence))

    def test_onlist_reader_finds_tokens_and_offlist_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "onlist.csv.gz"
            path.write_bytes(MODULE.gzip.compress(b"name,index\none,AAAA\ntwo,AAAC\n"))
            entries = MODULE.read_onlist_entries(path, 4)

        self.assertEqual(entries, {"AAAA", "AAAC"})
        self.assertEqual(MODULE.first_offlist_sequence(4, entries), "AAAG")

    def test_bundle_local_resources_rewrites_and_hashes_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "source"
            source.mkdir()
            (source / "onlists").mkdir()
            (source / "fastqs").mkdir()
            (source / "onlists" / "barcodes.txt").write_text("AAAA\n")
            (source / "fastqs" / "read.fastq").write_text("@r\nAAAA\n+\nIIII\n")
            spec_path = source / "spec.yaml"
            spec_path.write_text("seqspec_version: 0.5.0\n")
            value = assay()
            MODULE.find_region(value, "left")["onlist"] = {
                "urltype": "local",
                "url": "onlists/barcodes.txt",
                "filename": "barcodes.txt",
            }
            MODULE.find_read(value, "read1")["files"] = [
                {
                    "urltype": "local",
                    "url": "fastqs/read.fastq",
                    "filename": "read.fastq",
                }
            ]
            output = root / "output"
            resources = MODULE.bundle_local_resources(
                value,
                source_spec_path=spec_path,
                output_root=output,
            )

        self.assertEqual(len(resources), 2)
        self.assertTrue(all(len(resource["sha256"]) == 64 for resource in resources))
        self.assertTrue(
            MODULE.find_region(value, "left")["onlist"]["url"].startswith("resources/")
        )
        self.assertTrue(
            MODULE.find_read(value, "read1")["files"][0]["url"].startswith("resources/")
        )

    def test_yq_conversion_requires_json_object(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            yq = root / "yq"
            yq.write_text('#!/bin/sh\necho \'{"seqspec_version": "0.5.0"}\'\n')
            yq.chmod(0o755)
            path = root / "spec.yaml"
            path.write_text("seqspec_version: 0.5.0\n")
            value = MODULE.load_seqspec_json(path, yq_bin=yq, timeout_seconds=5)

        self.assertEqual(value["seqspec_version"], "0.5.0")


if __name__ == "__main__":
    unittest.main()
