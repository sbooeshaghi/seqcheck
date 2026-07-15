import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "complementarity_runtime.py"
)
SPEC = importlib.util.spec_from_file_location("complementarity_runtime", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def record(index: int, sequence: str = "AACCTTGG") -> MODULE.FastqRecord:
    return (
        f"@record-{index}\n".encode(),
        f"{sequence}\n".encode(),
        b"+\n",
        ("I" * len(sequence) + "\n").encode(),
    )


class ComplementarityRuntimeTests(unittest.TestCase):
    def test_quality_operators_preserve_names_and_sequences(self) -> None:
        records = [record(index) for index in range(100)]
        cases = [
            ("Q01", {"phred": 2}),
            ("Q02", {"phred": 2, "cycle_fraction": 0.25}),
            ("Q03", {"phred": 2, "event_fraction": 0.1}),
            (
                "Q04",
                {"low_phred": 2, "high_phred": 40, "low_fraction": 0.25},
            ),
        ]

        for operator_id, severity in cases:
            with self.subTest(operator_id=operator_id):
                mutated, details = MODULE.mutate_quality_records(
                    records,
                    operator_id=operator_id,
                    severity=severity,
                    seed=401,
                )
                self.assertEqual(
                    [(value[0], value[1], value[2]) for value in records],
                    [(value[0], value[1], value[2]) for value in mutated],
                )
                self.assertGreater(details["records_changed"], 0)
                self.assertTrue(
                    all(
                        len(value[1].strip()) == len(value[3].strip())
                        for value in mutated
                    )
                )

    def test_quality_selection_is_seeded_and_deterministic(self) -> None:
        records = [record(index) for index in range(100)]
        severity = {"phred": 2, "event_fraction": 0.1}
        first, first_details = MODULE.mutate_quality_records(
            records, operator_id="Q03", severity=severity, seed=401
        )
        second, second_details = MODULE.mutate_quality_records(
            records, operator_id="Q03", severity=severity, seed=401
        )

        self.assertEqual(first, second)
        self.assertEqual(first_details, second_details)
        self.assertEqual(first_details["records_selected"], 10)

    def test_composition_mutations_preserve_names_lengths_and_coordinates(self) -> None:
        bases = ("A", "C", "G", "T")
        records = [
            record(index, f"AA{bases[index % len(bases)]}CTTGG") for index in range(100)
        ]
        for operator_id, options in (
            ("C01", {"target_span": (2, 6), "motif": "AGAT"}),
            ("C03", {"target_span": (2, 6)}),
        ):
            with self.subTest(operator_id=operator_id):
                mutated, details = MODULE.mutate_composition_records(
                    records,
                    operator_id=operator_id,
                    event_fraction=0.1,
                    seed=503,
                    **options,
                )
                changed = [
                    index
                    for index, (before, after) in enumerate(zip(records, mutated))
                    if before[1] != after[1]
                ]
                self.assertEqual(len(changed), 10)
                self.assertEqual(details["records_selected"], 10)
                for index in changed:
                    self.assertEqual(
                        records[index][1].strip()[:2], mutated[index][1].strip()[:2]
                    )
                    self.assertEqual(
                        records[index][1].strip()[6:], mutated[index][1].strip()[6:]
                    )

    def test_duplication_preserves_record_names_and_count(self) -> None:
        bases = ("A", "C", "G", "T")
        records = [
            record(index, "".join(bases[(index + offset) % 4] for offset in range(8)))
            for index in range(100)
        ]
        mutated, details = MODULE.mutate_composition_records(
            records,
            operator_id="C02",
            event_fraction=0.1,
            seed=509,
        )

        self.assertEqual(
            [value[0] for value in records], [value[0] for value in mutated]
        )
        self.assertEqual(len(records), len(mutated))
        self.assertEqual(details["records_selected"], 10)

    def test_schema_invalid_mutations_change_one_declared_constraint(self) -> None:
        spec = {
            "assay_id": "example",
            "library_spec": [
                {
                    "region_id": "rna",
                    "min_len": 8,
                    "max_len": 8,
                    "regions": [],
                }
            ],
        }
        missing, missing_details = MODULE.mutate_invalid_spec(
            spec, variant="missing_required_root_region_id"
        )
        geometry, geometry_details = MODULE.mutate_invalid_spec(
            spec, variant="region_min_exceeds_max"
        )

        self.assertNotIn("region_id", missing["library_spec"][0])
        self.assertEqual(missing_details["removed_field"], "region_id")
        self.assertGreater(
            geometry["library_spec"][0]["min_len"],
            geometry["library_spec"][0]["max_len"],
        )
        self.assertEqual(geometry_details["region_id"], "rna")


if __name__ == "__main__":
    unittest.main()
