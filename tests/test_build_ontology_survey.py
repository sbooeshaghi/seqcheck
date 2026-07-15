import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "build_ontology_survey.py"
SPEC = importlib.util.spec_from_file_location("build_ontology_survey_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


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


def fake_seqspec(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import json
import sys
from pathlib import Path

MAPPING = {
    "rna": ["RGN:unknown:unclassified"],
    "barcode": ["RGN:partition:cell"],
    "umi": ["RGN:partition:molecule"],
    "index7": ["RGN:partition:sample", "RGN:technical:index7"],
}

if sys.argv[1:] == ["--version"]:
    print("seqspec test 0.5.0")
elif sys.argv[1] == "upgrade":
    output = Path(sys.argv[sys.argv.index("-o") + 1])
    source = Path(sys.argv[-1])
    value = json.loads(source.read_text())
    value["seqspec_version"] = "0.5.0"
    for child in value["library_spec"][0]["regions"]:
        child["region_type"] = MAPPING[child["region_type"]]
    output.write_text(json.dumps(value))
else:
    raise SystemExit(2)
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def region(
    region_id: str,
    region_type: str | list[str],
    *,
    sequence_type: str = "fixed",
    children: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    length = sum(child["max_len"] for child in children or []) or 4
    return {
        "region_id": region_id,
        "region_type": region_type,
        "name": region_id.replace("_", " ").title(),
        "sequence_type": "joined" if children else sequence_type,
        "sequence": "X" * length,
        "min_len": length,
        "max_len": length,
        "onlist": None,
        "regions": children,
    }


def specs() -> tuple[dict[str, object], dict[str, object]]:
    raw_children = [
        region("barcode_a", "barcode", sequence_type="onlist"),
        region("barcode_b", "barcode", sequence_type="onlist"),
        region("umi", "umi", sequence_type="random"),
        region("index7", "index7", sequence_type="onlist"),
        region("custom", "custom_unmapped"),
    ]
    normalized_children = [
        region("barcode_a", ["RGN:partition:cell"], sequence_type="onlist"),
        region("barcode_b", ["RGN:partition:cell"], sequence_type="onlist"),
        region("umi", ["RGN:partition:molecule"], sequence_type="random"),
        region(
            "index7",
            ["RGN:partition:sample", "RGN:technical:index7"],
            sequence_type="onlist",
        ),
        region("custom", ["RGN:unknown:unclassified"]),
    ]
    raw_children[0]["onlist"] = {
        "location": "local",
        "filename": "barcodes.txt",
        "md5": "abc123",
    }
    normalized_children[0]["onlist"] = {
        "file_id": "barcodes.txt",
        "filename": "barcodes.txt",
        "filetype": "",
        "filesize": 0,
        "url": "",
        "urltype": "",
        "md5": "abc123",
    }
    read = {
        "read_id": "R1",
        "name": "Read 1",
        "modality": "rna",
        "primer_id": "barcode_a",
        "min_len": 20,
        "max_len": 20,
        "strand": "pos",
    }
    common = {
        "assay_id": "synthetic",
        "name": "Synthetic assay",
        "description": "Ontology survey fixture",
        "modalities": ["rna"],
        "sequence_spec": [read],
    }
    raw = {
        **common,
        "seqspec_version": "0.2.0",
        "library_spec": [region("rna", "rna", children=raw_children)],
    }
    normalized = {
        **common,
        "seqspec_version": "0.5.0",
        "library_spec": [
            region("rna", ["RGN:unknown:unclassified"], children=normalized_children)
        ],
    }
    return raw, normalized


def registry() -> dict[str, object]:
    term_ids = [
        "RGN:unknown:unclassified",
        "RGN:partition:cell",
        "RGN:partition:molecule",
        "RGN:partition:sample",
        "RGN:technical:index7",
    ]
    return {
        "ontology_id": "seqspec-region-ontology",
        "legacy_region_types": {
            "rna": ["RGN:unknown:unclassified"],
            "barcode": ["RGN:partition:cell"],
            "umi": ["RGN:partition:molecule"],
            "index7": ["RGN:partition:sample", "RGN:technical:index7"],
        },
        "terms": {term: {"status": "active"} for term in term_ids},
    }


def create_inputs(root: Path) -> dict[str, Path]:
    yq = root / "bin" / "yq"
    seqspec = root / "bin" / "seqspec"
    yq.parent.mkdir(parents=True)
    fake_yq(yq)
    fake_seqspec(seqspec)
    raw, normalized = specs()
    raw_path = root / "specs" / "raw" / "CONFIG1.yaml"
    normalized_path = root / "specs" / "normalized" / "CONFIG1.yaml"
    write_json(raw_path, raw)
    write_json(normalized_path, normalized)

    proposals = root / "tables" / "cohort_proposals.csv"
    proposals.parent.mkdir(parents=True)
    fields = [
        "selection_id",
        "family_id",
        "family_label",
        "configuration_accession",
        "assay_term",
        "preferred_assay_titles",
        "hydration_status",
        "raw_spec_sha256",
        "raw_seqspec_version",
        "normalized_spec_sha256",
        "normalized_seqspec_version",
        "normalized_spec_path",
    ]
    with proposals.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(
            {
                "selection_id": "selection-1",
                "family_id": "rna_family",
                "family_label": "RNA family",
                "configuration_accession": "CONFIG1",
                "assay_term": "single-cell RNA sequencing assay",
                "preferred_assay_titles": "synthetic",
                "hydration_status": "normalized",
                "raw_spec_sha256": MODULE.runtime.file_sha256(raw_path),
                "raw_seqspec_version": "0.2.0",
                "normalized_spec_sha256": MODULE.runtime.file_sha256(normalized_path),
                "normalized_seqspec_version": "0.5.0",
                "normalized_spec_path": str(normalized_path),
            }
        )

    cohort = root / "manifests" / "cohort.json"
    write_json(
        cohort,
        {
            "cohort_candidate_schema_version": "0.2.0",
            "selection_id": "selection-1",
            "outputs": {"cohort_proposals": str(proposals)},
        },
    )
    registry_path = root / "registry.yaml"
    write_json(registry_path, registry())
    protocol = MODULE.runtime.load_json(
        ROOT / "experiments" / "paper" / "protocol" / "ontology_mapping.json"
    )
    protocol["manual_review"] |= {
        "sample_size": 4,
        "rare_label_max_population_count": 1,
        "target_allocations": {
            "unknown": 1,
            "multi_term": 1,
            "rare_label": 1,
            "common": 1,
        },
    }
    protocol_path = root / "protocol.json"
    write_json(protocol_path, protocol)
    return {
        "cohort": cohort,
        "proposals": proposals,
        "raw": raw_path,
        "normalized": normalized_path,
        "registry": registry_path,
        "protocol": protocol_path,
        "seqspec": seqspec,
        "yq": yq,
    }


class BuildOntologySurveyTests(unittest.TestCase):
    def test_builds_complete_deterministic_weighted_survey(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            inputs = create_inputs(Path(tmpdir) / "run")
            first = Path(tmpdir) / "survey-first"
            second = Path(tmpdir) / "survey-second"
            first_manifest = MODULE.build_ontology_survey(
                cohort_manifest_path=inputs["cohort"],
                registry_path=inputs["registry"],
                protocol_path=inputs["protocol"],
                seqspec_bin=inputs["seqspec"],
                yq_bin=inputs["yq"],
                output_root=first,
                timeout_seconds=10,
            )
            second_manifest = MODULE.build_ontology_survey(
                cohort_manifest_path=inputs["cohort"],
                registry_path=inputs["registry"],
                protocol_path=inputs["protocol"],
                seqspec_bin=inputs["seqspec"],
                yq_bin=inputs["yq"],
                output_root=second,
                timeout_seconds=10,
            )

            self.assertEqual(first_manifest["survey_id"], second_manifest["survey_id"])
            self.assertEqual(first_manifest["counts"]["specifications"], 1)
            self.assertEqual(first_manifest["counts"]["regions"], 6)
            self.assertEqual(first_manifest["counts"]["review_sample"], 4)
            self.assertEqual(
                first_manifest["counts"]["review_strata"],
                {"common": 1, "multi_term": 1, "rare_label": 1, "unknown": 1},
            )
            self.assertEqual(first_manifest["counts"]["unexpected_query_changes"], 0)
            self.assertGreater(
                first_manifest["counts"]["expected_unknown_query_changes"], 0
            )

            regions = read_csv(first / "tables" / "ontology_regions.csv")
            root = next(row for row in regions if row["region_id"] == "rna")
            index7 = next(row for row in regions if row["region_id"] == "index7")
            self.assertEqual(root["depth"], "0")
            self.assertEqual(root["path_region_ids"], "rna")
            self.assertEqual(index7["parent_region_id"], "rna")
            self.assertEqual(
                index7["ontology_terms"],
                "RGN:partition:sample;RGN:technical:index7",
            )
            sample = read_csv(first / "tables" / "review_sample.csv")
            self.assertEqual(
                {row["review_stratum"] for row in sample},
                {
                    "unknown",
                    "multi_term",
                    "rare_label",
                    "common",
                },
            )
            for row in sample:
                probability = float(row["inclusion_probability"])
                self.assertAlmostEqual(float(row["sampling_weight"]), 1 / probability)
            self.assertEqual(
                [row["region_key"] for row in sample],
                [
                    row["region_key"]
                    for row in read_csv(second / "tables" / "review_sample.csv")
                ],
            )

            queries = read_csv(first / "tables" / "query_parity.csv")
            custom = next(
                row for row in queries if row["selector"] == "custom_unmapped"
            )
            self.assertEqual(custom["change_class"], "expected_unknown_collapse")
            self.assertEqual(custom["raw_region_count"], "1")
            self.assertEqual(custom["normalized_region_count"], "0")
            barcode = next(row for row in queries if row["selector"] == "barcode")
            self.assertEqual(barcode["change_class"], "stable")
            self.assertEqual(barcode["raw_region_count"], "2")

    def test_cli_rejects_hash_tampering(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "run")
            inputs["raw"].write_text("tampered\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--cohort-manifest",
                    str(inputs["cohort"]),
                    "--registry",
                    str(inputs["registry"]),
                    "--protocol",
                    str(inputs["protocol"]),
                    "--seqspec-bin",
                    str(inputs["seqspec"]),
                    "--yq-bin",
                    str(inputs["yq"]),
                    "--output-root",
                    str(root / "survey"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 1)
            self.assertIn("raw spec hash changed", completed.stderr)

    def test_rejects_normalized_mapping_that_differs_from_registry(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "run")
            normalized = json.loads(inputs["normalized"].read_text())
            normalized["library_spec"][0]["regions"][0]["region_type"] = [
                "RGN:partition:molecule"
            ]
            write_json(inputs["normalized"], normalized)
            proposals = read_csv(inputs["proposals"])
            proposals[0]["normalized_spec_sha256"] = MODULE.runtime.file_sha256(
                inputs["normalized"]
            )
            with inputs["proposals"].open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(proposals[0]))
                writer.writeheader()
                writer.writerows(proposals)

            with self.assertRaisesRegex(ValueError, "ontology mapping differs"):
                MODULE.build_ontology_survey(
                    cohort_manifest_path=inputs["cohort"],
                    registry_path=inputs["registry"],
                    protocol_path=inputs["protocol"],
                    seqspec_bin=inputs["seqspec"],
                    yq_bin=inputs["yq"],
                    output_root=root / "survey",
                    timeout_seconds=10,
                )

    def test_rejects_seqspec_runtime_mapping_that_differs_from_registry(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_inputs(root / "run")
            executable = inputs["seqspec"].read_text()
            inputs["seqspec"].write_text(
                executable.replace(
                    '"barcode": ["RGN:partition:cell"]',
                    '"barcode": ["RGN:partition:molecule"]',
                )
            )

            with self.assertRaisesRegex(ValueError, "runtime mapping differs"):
                MODULE.build_ontology_survey(
                    cohort_manifest_path=inputs["cohort"],
                    registry_path=inputs["registry"],
                    protocol_path=inputs["protocol"],
                    seqspec_bin=inputs["seqspec"],
                    yq_bin=inputs["yq"],
                    output_root=root / "survey",
                    timeout_seconds=10,
                )


if __name__ == "__main__":
    unittest.main()
