import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests import test_materialize_complementarity as control_fixture
from tests import test_run_perturbation_calibration as execution_fixture


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "run_complementarity.py"
SPEC = importlib.util.spec_from_file_location("run_complementarity", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fake_fastqc(path: Path) -> None:
    path.write_text(
        r"""#!/usr/bin/env python3
import gzip
import sys
from pathlib import Path

if "--version" in sys.argv:
    print("FastQC v0.12.1-test")
    raise SystemExit(0)

source = Path(sys.argv[-1])
output_root = Path(sys.argv[sys.argv.index("--outdir") + 1])
name = source.name
lower = name.lower()
for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
    if lower.endswith(suffix):
        stem = name[:-len(suffix)] + "_fastqc"
        break
else:
    stem = name + "_fastqc"
opener = gzip.open if lower.endswith(".gz") else open
with opener(source, "rt", encoding="ascii") as handle:
    lines = handle.read().splitlines()
qualities = lines[3::4]
scores = [ord(value) - 33 for quality in qualities for value in quality]
mean = sum(scores) / len(scores)
quality_status = "fail" if mean < 15 else "pass"
modules = [
    (
        "Basic Statistics",
        "pass",
        ["#Measure\tValue", f"Filename\t{name}", f"Total Sequences\t{len(qualities)}"],
    ),
    (
        "Per base sequence quality",
        quality_status,
        ["#Base\tMean\tMedian", f"1\t{mean:.3f}\t{mean:.3f}"],
    ),
    (
        "Per sequence quality scores",
        quality_status,
        ["#Quality\tCount", f"{round(mean)}\t{len(qualities)}"],
    ),
    ("Per sequence GC content", "pass", ["#GC Content\tCount", "50\t1"]),
    (
        "Sequence Duplication Levels",
        "pass",
        ["#Total Deduplicated Percentage\t100.0", "#Duplication Level\tPercentage of deduplicated\tPercentage of total", "1\t100\t100"],
    ),
    (
        "Overrepresented sequences",
        "pass",
        ["#Sequence\tCount\tPercentage\tPossible Source"],
    ),
    (
        "Adapter Content",
        "pass",
        ["#Position\tIllumina Universal Adapter", "1\t0.0"],
    ),
]
report = output_root / stem
report.mkdir(parents=True)
data = ["##FastQC\t0.12.1-test"]
summary = []
for module, status, rows in modules:
    data.append(f">>{module}\t{status}")
    data.extend(rows)
    data.append(">>END_MODULE")
    summary.append(f"{status.upper()}\t{module}\t{name}")
(report / "fastqc_data.txt").write_text("\n".join(data) + "\n", encoding="utf-8")
(report / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def create_materialization(root: Path) -> dict[str, Path]:
    base = control_fixture.create_base(root / "study", cohort_split="calibration")
    protocol = control_fixture.create_protocol(root / "protocol.json")
    manifest = control_fixture.MODULE.materialize_complementarity(
        base_materialization_path=base["base"],
        protocol_path=protocol,
        seqspec_bin=base["seqspec"],
        yq_bin=base["yq"],
        cohort_split="calibration",
        output_root=root / "controls",
        timeout_seconds=10,
    )
    seqcheck = root / "seqcheck"
    fastqc = root / "fastqc"
    execution_fixture.fake_seqcheck(seqcheck)
    fake_fastqc(fastqc)
    return {
        "materialization": Path(manifest["manifest_path"]),
        "protocol": protocol,
        "seqcheck": seqcheck,
        "seqspec": base["seqspec"],
        "fastqc": fastqc,
    }


class RunComplementarityTests(unittest.TestCase):
    def test_cli_runs_three_tools_and_reuses_fastqc_by_content(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_materialization(root)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--materialization-manifest",
                    str(inputs["materialization"]),
                    "--protocol",
                    str(inputs["protocol"]),
                    "--seqcheck-bin",
                    str(inputs["seqcheck"]),
                    "--seqspec-bin",
                    str(inputs["seqspec"]),
                    "--fastqc-bin",
                    str(inputs["fastqc"]),
                    "--output-root",
                    str(root / "execution"),
                    "--timeout-seconds",
                    "10",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest = json.loads(
                Path(completed.stdout.strip()).read_text(encoding="utf-8")
            )
            seqspec = read_csv(manifest["outputs"]["seqspec_runs"]["path"])
            seqcheck = read_csv(manifest["outputs"]["seqcheck_runs"]["path"])
            fastqc = read_csv(manifest["outputs"]["fastqc_runs"]["path"])
            modules = read_csv(manifest["outputs"]["fastqc_modules"]["path"])
            performance = read_csv(manifest["outputs"]["fastqc_performance"]["path"])

        self.assertTrue(manifest["valid"])
        self.assertEqual(len(seqspec), manifest["counts"]["conditions"])
        self.assertEqual(len(seqcheck), manifest["counts"]["conditions"])
        self.assertEqual(len(fastqc), manifest["counts"]["fastqc_input_runs"])
        self.assertEqual(len(modules), manifest["counts"]["fastqc_modules"])
        self.assertEqual(
            len(performance), manifest["counts"]["fastqc_cache_executions"]
        )
        self.assertGreater(manifest["counts"]["fastqc_cache_hits"], 0)
        invalid = [value for value in seqspec if value["operator_id"].startswith("X")]
        self.assertTrue(all(value["observed_status"] == "failure" for value in invalid))
        quality_modules = [
            value
            for value in modules
            if value["operator_id"].startswith("Q")
            and value["module_name"]
            in {
                "Per base sequence quality",
                "Per sequence quality scores",
            }
        ]
        self.assertTrue(any(value["status"] == "fail" for value in quality_modules))

    def test_loader_rejects_changed_control_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = create_materialization(root)
            materialization = MODULE.runtime.load_json(inputs["materialization"])
            payload = MODULE.runtime.load_json(
                Path(materialization["outputs"]["controls"]["path"])
            )
            artifact = Path(payload["controls"][0]["inputs"][0]["path"])
            artifact.write_bytes(artifact.read_bytes() + b"changed")

            with self.assertRaisesRegex(ValueError, "artifact hash changed"):
                MODULE.load_materialization(inputs["materialization"])


if __name__ == "__main__":
    unittest.main()
