"""Strict FastQC output parsing for paper experiments."""

from __future__ import annotations

from pathlib import Path
from typing import Any


RUNTIME_VERSION = "0.1.0"
VALID_STATUSES = ("pass", "warn", "fail")
STATUS_RANK = {value: index for index, value in enumerate(VALID_STATUSES)}


def parse_fastqc_output(root: Path) -> dict[str, Any]:
    data_path = root / "fastqc_data.txt"
    summary_path = root / "summary.txt"
    data = parse_fastqc_data(data_path)
    summary = parse_summary(summary_path)
    data_statuses = [
        (value["status"], value["name"], data["filename"]) for value in data["modules"]
    ]
    summary_statuses = [
        (value["status"], value["name"], value["filename"]) for value in summary
    ]
    if data_statuses != summary_statuses:
        raise ValueError("FastQC summary does not reconcile with module data")
    return {
        **data,
        "data_path": str(data_path.resolve()),
        "summary_path": str(summary_path.resolve()),
    }


def parse_fastqc_data(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    if not lines or not lines[0].startswith("##FastQC\t"):
        raise ValueError("FastQC data header is missing")
    version = lines[0].split("\t", 1)[1].strip()
    if not version:
        raise ValueError("FastQC version is missing")
    modules = []
    current = None
    for line_number, line in enumerate(lines[1:], 2):
        if line.startswith(">>"):
            if line == ">>END_MODULE":
                if current is None:
                    raise ValueError(
                        f"FastQC data has an unmatched module end at line {line_number}"
                    )
                modules.append(current)
                current = None
                continue
            if current is not None:
                raise ValueError(
                    f"FastQC data has a nested module at line {line_number}"
                )
            fields = line[2:].split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"FastQC module header is invalid at line {line_number}"
                )
            name, status = fields[0].strip(), fields[1].strip().lower()
            if not name or status not in VALID_STATUSES:
                raise ValueError(
                    f"FastQC module name or status is invalid at line {line_number}"
                )
            current = {"name": name, "status": status, "headers": [], "rows": []}
            continue
        if current is None:
            if line.strip():
                raise ValueError(f"FastQC data outside a module at line {line_number}")
            continue
        fields = line.split("\t")
        if line.startswith("#"):
            current["headers"].append([fields[0][1:], *fields[1:]])
        else:
            current["rows"].append(fields)
    if current is not None:
        raise ValueError("FastQC data ends before END_MODULE")
    names = [value["name"] for value in modules]
    if not modules or len(names) != len(set(names)):
        raise ValueError("FastQC modules are empty or duplicated")
    filename = basic_statistic(modules, "Filename")
    return {"fastqc_version": version, "filename": filename, "modules": modules}


def parse_summary(path: Path) -> list[dict[str, str]]:
    rows = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8", errors="strict").splitlines(), 1
    ):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError(f"FastQC summary row is invalid at line {line_number}")
        status = fields[0].strip().lower()
        name = fields[1].strip()
        filename = fields[2].strip()
        if status not in VALID_STATUSES or not name or not filename:
            raise ValueError(f"FastQC summary value is invalid at line {line_number}")
        rows.append({"status": status, "name": name, "filename": filename})
    if not rows:
        raise ValueError("FastQC summary is empty")
    return rows


def basic_statistic(modules: list[dict[str, Any]], measure: str) -> str:
    basic = next(
        (value for value in modules if value["name"] == "Basic Statistics"), None
    )
    if basic is None:
        raise ValueError("FastQC Basic Statistics module is missing")
    matches = [
        row[1]
        for row in basic["rows"]
        if len(row) == 2 and row[0] == measure and row[1].strip()
    ]
    if len(matches) != 1:
        raise ValueError(f"FastQC Basic Statistics has no unique {measure}")
    return matches[0]


def status_worsened(clean: str, observed: str) -> bool:
    if clean not in STATUS_RANK or observed not in STATUS_RANK:
        raise ValueError("FastQC status is invalid")
    return STATUS_RANK[observed] > STATUS_RANK[clean]


def output_stem(path: Path) -> str:
    name = path.name
    lowered = name.lower()
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam", ".sam"):
        if lowered.endswith(suffix):
            return name[: -len(suffix)] + "_fastqc"
    return name + "_fastqc"
