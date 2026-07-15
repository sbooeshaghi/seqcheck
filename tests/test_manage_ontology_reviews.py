import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "manage_ontology_reviews.py"
SPEC = importlib.util.spec_from_file_location(
    "manage_ontology_reviews_test", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def fake_yq(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if sys.argv[1:] == ["--version"]:
    print("yq test 4.0")
else:
    print(Path(sys.argv[-1]).read_text())
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def registry() -> dict[str, object]:
    roles = {
        "partition": {"scope": "direct_read_interval"},
        "technical": {"scope": "direct_read_interval"},
        "unknown": {"scope": "direct_read_interval"},
    }
    values = {
        "RGN:partition:cell": ("Cell partition barcode", "partition", "cell"),
        "RGN:partition:molecule": (
            "Molecule partition barcode",
            "partition",
            "molecule",
        ),
        "RGN:partition:sample": (
            "Sample partition barcode",
            "partition",
            "sample",
        ),
        "RGN:technical:index7": ("Index 7", "technical", "index7"),
        "RGN:unknown:unclassified": (
            "Unclassified region",
            "unknown",
            "unclassified",
        ),
    }
    terms = {}
    for term, (label, role, target) in values.items():
        terms[term] = {
            "label": label,
            "role": role,
            "target": target,
            "status": "active",
            "definition": f"Definition for {label}.",
            "typical_sequence_types": ["onlist"],
            "examples": [f"Example for {label}."],
            "does_not_mean": [f"Not another role for {label}."],
            "aliases": ["hidden legacy hint"],
        }
    return {
        "ontology_id": "seqspec-region-ontology",
        "roles": roles,
        "terms": terms,
        "legacy_region_types": {"barcode": ["RGN:partition:cell"]},
    }


def create_source(root: Path) -> dict[str, Path]:
    yq = root / "bin" / "yq"
    yq.parent.mkdir(parents=True)
    fake_yq(yq)
    registry_path = root / "registry.yaml"
    write_json(registry_path, registry())
    protocol = MODULE.runtime.load_json(
        ROOT / "experiments" / "paper" / "protocol" / "ontology_mapping.json"
    )
    protocol["manual_review"]["sample_size"] = 4
    protocol["manual_review"]["target_allocations"] = {
        "unknown": 1,
        "multi_term": 1,
        "rare_label": 1,
        "common": 1,
    }
    protocol_path = root / "protocol.json"
    write_json(protocol_path, protocol)

    truths = [
        "RGN:partition:cell",
        "RGN:partition:molecule",
        "RGN:partition:sample;RGN:technical:index7",
        "RGN:unknown:unclassified",
    ]
    labels = ["barcode", "umi", "index7", "named"]
    strata = ["common", "rare_label", "multi_term", "unknown"]
    rows = []
    for index, (truth, label, stratum) in enumerate(
        zip(truths, labels, strata, strict=True), 1
    ):
        row = {
            "survey_id": "survey-1",
            "region_key": f"region-key-{index}",
            "ontology_terms": truth,
            "mapping_source": "deterministic_registry",
            "unknown_mapping": str("unknown" in truth),
            "multi_term_mapping": str(";" in truth),
            "review_stratum": stratum,
            "stratum_population": str(index * 10),
            "stratum_sample_size": "1",
            "inclusion_probability": str(1 / (index * 10)),
            "sampling_weight": str(index * 10),
            "selection_rank": "1",
        }
        row.update(
            {
                field: value
                for field, value in zip(
                    MODULE.CONTEXT_FIELDS,
                    (
                        "CONFIG1",
                        "RNA family",
                        "single-cell RNA sequencing assay",
                        "synthetic",
                        "assay-1",
                        "Synthetic assay",
                        "Synthetic assay context",
                        "0.3.0",
                        "rna",
                        "R1",
                        "primer",
                        f"region-{index}",
                        f"Region {index}",
                        "onlist",
                        "4",
                        "4",
                        "1",
                        "True",
                        "rna",
                        "RNA",
                        f"rna;region-{index}",
                        f"RNA;Region {index}",
                        json.dumps(label),
                        label,
                        "True",
                    ),
                    strict=True,
                )
            }
        )
        rows.append(row)
    sample_fields = [
        "survey_id",
        "region_key",
        *MODULE.CONTEXT_FIELDS,
        "ontology_terms",
        "mapping_source",
        "unknown_mapping",
        "multi_term_mapping",
        "review_stratum",
        "stratum_population",
        "stratum_sample_size",
        "inclusion_probability",
        "sampling_weight",
        "selection_rank",
    ]
    sample_path = root / "tables" / "review_sample.csv"
    write_csv(sample_path, rows, sample_fields)
    manifest_path = root / "manifests" / "survey.json"
    write_json(
        manifest_path,
        {
            "schema_version": "0.1.0",
            "survey_id": "survey-1",
            "valid": True,
            "counts": {"review_sample": 4},
            "registry_runtime_parity": {"legacy_labels": 40, "matched_labels": 40},
            "inputs": {
                "registry": MODULE.runtime.file_identity(registry_path),
                "protocol": MODULE.runtime.file_identity(protocol_path),
            },
            "outputs": {"review_sample": MODULE.runtime.file_identity(sample_path)},
        },
    )
    return {
        "manifest": manifest_path,
        "registry": registry_path,
        "protocol": protocol_path,
        "sample": sample_path,
        "yq": yq,
    }


def prepare(root: Path, source: dict[str, Path]) -> Path:
    output = root / "packages"
    MODULE.prepare_review_packages(
        survey_manifest_path=source["manifest"],
        registry_path=source["registry"],
        protocol_path=source["protocol"],
        yq_bin=source["yq"],
        output_root=output,
        timeout_seconds=10,
    )
    return output


def complete_sheet(
    packages: Path,
    slot: int,
    reviewer: str,
    mutate: Callable[[list[dict[str, str]]], None] | None = None,
) -> Path:
    key_rows, _ = read_csv(packages / "internal" / "review_key.csv")
    key_field = f"reviewer_{slot}_item_id"
    truth_by_item = {row[key_field]: row["ontology_terms"] for row in key_rows}
    sheet = packages / f"reviewer_{slot}" / "ontology_review.csv"
    rows, fields = read_csv(sheet)
    for row in rows:
        terms = truth_by_item[row["review_item_id"]]
        if slot == 2 and ";" in terms:
            terms = ";".join(reversed(terms.split(";")))
        row.update(
            {
                "reviewer": reviewer,
                "review_date": "2026-07-14",
                "reviewed_ontology_terms": terms,
                "context_sufficient": "yes",
                "confidence": "high",
                "notes": "",
            }
        )
    if mutate is not None:
        mutate(rows)
    write_csv(sheet, rows, fields)
    return sheet


def merge(root: Path, source: dict[str, Path], packages: Path) -> dict[str, object]:
    return MODULE.merge_review_packages(
        survey_manifest_path=source["manifest"],
        registry_path=source["registry"],
        protocol_path=source["protocol"],
        yq_bin=source["yq"],
        reviewer_1_package=packages / "reviewer_1" / "review_package.json",
        reviewer_1_sheet=packages / "reviewer_1" / "ontology_review.csv",
        reviewer_2_package=packages / "reviewer_2" / "review_package.json",
        reviewer_2_sheet=packages / "reviewer_2" / "ontology_review.csv",
        output_root=root / "merged",
        timeout_seconds=10,
    )


class ManageOntologyReviewTests(unittest.TestCase):
    def test_prepare_blinds_truth_and_builds_independent_packages(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = create_source(root / "source")
            packages = prepare(root, source)
            first, first_fields = read_csv(
                packages / "reviewer_1" / "ontology_review.csv"
            )
            second, second_fields = read_csv(
                packages / "reviewer_2" / "ontology_review.csv"
            )
            self.assertEqual(first_fields, second_fields)
            self.assertFalse(set(MODULE.BLINDED_FIELDS) & set(first_fields))
            self.assertTrue(set(MODULE.EDITABLE_FIELDS) <= set(first_fields))
            self.assertNotEqual(
                [row["review_item_id"] for row in first],
                [row["review_item_id"] for row in second],
            )
            self.assertTrue(
                set(row["review_item_id"] for row in first).isdisjoint(
                    row["review_item_id"] for row in second
                )
            )
            self.assertNotEqual(
                [row["region_id"] for row in first],
                [row["region_id"] for row in second],
            )
            reference, reference_fields = read_csv(
                packages / "reviewer_1" / "ontology_terms.csv"
            )
            self.assertEqual(reference_fields, list(MODULE.TERM_REFERENCE_FIELDS))
            self.assertNotIn("aliases", reference_fields)
            self.assertNotIn("legacy_region_types", reference_fields)
            self.assertEqual(len(reference), 5)
            private, private_fields = read_csv(packages / "internal" / "review_key.csv")
            self.assertEqual(private_fields, list(MODULE.PRIVATE_FIELDS))
            self.assertTrue(all(row["ontology_terms"] for row in private))
            self.assertTrue(all(float(row["sampling_weight"]) > 0 for row in private))
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                MODULE.prepare_review_packages(
                    survey_manifest_path=source["manifest"],
                    registry_path=source["registry"],
                    protocol_path=source["protocol"],
                    yq_bin=source["yq"],
                    output_root=packages,
                    timeout_seconds=10,
                )

    def test_merge_canonicalizes_term_order_and_recovers_truth(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = create_source(root / "source")
            packages = prepare(root, source)
            complete_sheet(packages, 1, "Reviewer A")
            complete_sheet(packages, 2, "Reviewer B")
            manifest = merge(root, source, packages)

            self.assertEqual(manifest["exact_term_set_agreements"], 4)
            self.assertEqual(manifest["adjudication_required"], 0)
            self.assertTrue(manifest["adjudication_complete"])
            merged, fields = read_csv(
                root / "merged" / "tables" / "ontology_reviews.csv"
            )
            self.assertEqual(
                fields[-len(MODULE.REVIEW_RESULT_FIELDS) :],
                list(MODULE.REVIEW_RESULT_FIELDS),
            )
            self.assertTrue(
                all(row["exact_term_set_agreement"] == "True" for row in merged)
            )
            self.assertTrue(
                all(row["reviewer_1_registry_exact_match"] == "True" for row in merged)
            )
            self.assertTrue(
                all(row["reviewer_2_registry_exact_match"] == "True" for row in merged)
            )
            index7 = next(
                row for row in merged if row["original_region_type_labels"] == "index7"
            )
            self.assertEqual(
                index7["reviewer_2_ontology_terms"],
                "RGN:partition:sample;RGN:technical:index7",
            )

    def test_merge_routes_term_disagreement_to_adjudication(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = create_source(root / "source")
            packages = prepare(root, source)
            complete_sheet(packages, 1, "Reviewer A")

            def disagree(rows: list[dict[str, str]]) -> None:
                rows[0]["reviewed_ontology_terms"] = "RGN:partition:molecule"

            complete_sheet(packages, 2, "Reviewer B", disagree)
            manifest = merge(root, source, packages)
            self.assertEqual(manifest["adjudication_required"], 1)
            self.assertFalse(manifest["adjudication_complete"])
            rows, _ = read_csv(root / "merged" / "tables" / "ontology_reviews.csv")
            disagreement = [row for row in rows if row["needs_adjudication"] == "True"]
            self.assertEqual(len(disagreement), 1)
            self.assertEqual(disagreement[0]["consensus_ontology_terms"], "")

    def test_merge_rejects_changed_context_and_package_identity(self) -> None:
        for scenario in ("context", "package", "reference"):
            with (
                self.subTest(scenario=scenario),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                root = Path(tmpdir)
                source = create_source(root / "source")
                packages = prepare(root, source)
                complete_sheet(packages, 1, "Reviewer A")
                complete_sheet(packages, 2, "Reviewer B")
                if scenario == "context":
                    sheet = packages / "reviewer_1" / "ontology_review.csv"
                    rows, fields = read_csv(sheet)
                    rows[0]["assay_name"] = "Changed assay"
                    write_csv(sheet, rows, fields)
                    pattern = "changed blinded context"
                elif scenario == "package":
                    package_path = packages / "reviewer_1" / "review_package.json"
                    package = json.loads(package_path.read_text())
                    package["package_id"] = "changed"
                    write_json(package_path, package)
                    pattern = "package id is invalid"
                else:
                    reference_path = packages / "reviewer_1" / "ontology_terms.csv"
                    reference, fields = read_csv(reference_path)
                    reference[0]["definition"] = "Changed definition"
                    write_csv(reference_path, reference, fields)
                    pattern = "term reference hash changed"
                with self.assertRaisesRegex(ValueError, pattern):
                    merge(root, source, packages)

    def test_merge_rejects_invalid_review_responses(self) -> None:
        scenarios = {
            "duplicate": (
                lambda rows: rows[0].update(
                    reviewed_ontology_terms=("RGN:partition:cell;RGN:partition:cell")
                ),
                "duplicates",
            ),
            "undeclared": (
                lambda rows: rows[0].update(reviewed_ontology_terms="RGN:not:declared"),
                "undeclared",
            ),
            "unknown_mixed": (
                lambda rows: rows[0].update(
                    reviewed_ontology_terms=(
                        "RGN:unknown:unclassified;RGN:partition:cell"
                    )
                ),
                "combined unknown",
            ),
            "insufficient_precise": (
                lambda rows: rows[0].update(context_sufficient="no"),
                "must use unknown",
            ),
            "missing_reviewer": (
                lambda rows: rows[0].update(reviewer=""),
                "identity is missing",
            ),
        }
        for name, (mutate, pattern) in scenarios.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                root = Path(tmpdir)
                source = create_source(root / "source")
                packages = prepare(root, source)
                complete_sheet(packages, 1, "Reviewer A", mutate)
                complete_sheet(packages, 2, "Reviewer B")
                with self.assertRaisesRegex(ValueError, pattern):
                    merge(root, source, packages)

    def test_merge_requires_distinct_reviewers(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = create_source(root / "source")
            packages = prepare(root, source)
            complete_sheet(packages, 1, "Reviewer A")
            complete_sheet(packages, 2, "Reviewer A")
            with self.assertRaisesRegex(ValueError, "distinct people"):
                merge(root, source, packages)

    def test_cli_prepares_and_merges_reviews(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = create_source(root / "source")
            packages = root / "packages"
            common = [
                "--survey-manifest",
                str(source["manifest"]),
                "--registry",
                str(source["registry"]),
                "--protocol",
                str(source["protocol"]),
                "--yq-bin",
                str(source["yq"]),
            ]
            prepared = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "prepare",
                    *common,
                    "--output-root",
                    str(packages),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(prepared.returncode, 0, prepared.stderr)
            complete_sheet(packages, 1, "Reviewer A")
            complete_sheet(packages, 2, "Reviewer B")
            merged = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "merge",
                    *common,
                    "--reviewer-1-package",
                    str(packages / "reviewer_1" / "review_package.json"),
                    "--reviewer-1-sheet",
                    str(packages / "reviewer_1" / "ontology_review.csv"),
                    "--reviewer-2-package",
                    str(packages / "reviewer_2" / "review_package.json"),
                    "--reviewer-2-sheet",
                    str(packages / "reviewer_2" / "ontology_review.csv"),
                    "--output-root",
                    str(root / "merged"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(merged.returncode, 0, merged.stderr)
            self.assertTrue((root / "merged/tables/ontology_reviews.csv").is_file())
            self.assertTrue((root / "merged/manifests/ontology_reviews.json").is_file())


if __name__ == "__main__":
    unittest.main()
