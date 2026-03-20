use serde_json::Value;
use std::path::Path;
use std::process::Command;

fn manifest_dir() -> &'static str {
    env!("CARGO_MANIFEST_DIR")
}

fn run_success(args: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_seqcheck"))
        .current_dir(manifest_dir())
        .args(args)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    String::from_utf8(output.stdout).unwrap()
}

fn run_failure(args: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_seqcheck"))
        .current_dir(manifest_dir())
        .args(args)
        .output()
        .unwrap();

    assert!(
        !output.status.success(),
        "command unexpectedly succeeded\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    String::from_utf8(output.stderr).unwrap()
}

#[test]
fn test_length_json_synthetic_golden() {
    let observed = run_success(&[
        "length",
        "--format",
        "json",
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);

    let expected = r#"{
  "spec": "tests/fixtures/synthetic/spec.yaml",
  "modality": "rna",
  "command": "length",
  "n_reads": 0,
  "input_check": {
    "expected_files": [
      {
        "read_id": "synthetic_R1",
        "file_id": "synthetic_R1.fastq",
        "filename": "synthetic_R1.fastq",
        "url_basename": "synthetic_R1.fastq"
      }
    ],
    "supplied_inputs": [
      "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq"
    ],
    "matched_inputs": [
      {
        "input_path": "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
        "read_id": "synthetic_R1",
        "file_id": "synthetic_R1.fastq",
        "matched_by": "file_id"
      }
    ],
    "missing_expected_files": []
  },
  "files": [
    {
      "input_path": "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
      "read_id": "synthetic_R1",
      "file_id": "synthetic_R1.fastq",
      "matched_by": "file_id",
      "results": {
        "sampled_count": 5,
        "expected_min_len": 12,
        "expected_max_len": 12,
        "observed_min_len": 2,
        "observed_max_len": 12,
        "out_of_range_count": 1,
        "out_of_range_fraction": 0.2
      }
    }
  ]
}
"#;

    assert_eq!(observed, expected);
}

#[test]
fn test_fixed_text_synthetic_golden() {
    let observed = run_success(&[
        "fixed",
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);

    let expected = r#"seqcheck fixed
spec: tests/fixtures/synthetic/spec.yaml
modality: rna
requested_reads: 0
input_check:
  expected_files:
    - synthetic_R1.fastq (read_id synthetic_R1, filename synthetic_R1.fastq, url_basename synthetic_R1.fastq)
  supplied_inputs:
    - tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
  matched_inputs:
    - tests/fixtures/synthetic/fastqs/synthetic_R1.fastq -> synthetic_R1.fastq (read_id synthetic_R1 via file_id)
  missing_expected_files:
    (none)

file: tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
read_id: synthetic_R1
file_id: synthetic_R1.fastq
matched_by: file_id
sampled_reads: 5
region: linker [4:6] expected=TT
  covered: 4 (0.8000)
  short_reads: 1
  exact_matches: 3 (0.7500)
  orientation_matches:
    forward: 3
    reverse: 3
    complement: 0
    reverse_complement: 0
  offset_histogram:
    -1 1
    +0 3
    +2 1
  absent_reads: 1
  multi_hit_reads: 1
  top_nonmatching_sequences:
    AC 1
"#;

    assert_eq!(observed, expected);
}

#[test]
fn test_onlist_text_synthetic_golden() {
    let observed = run_success(&[
        "onlist",
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);

    let expected = r#"seqcheck onlist
spec: tests/fixtures/synthetic/spec.yaml
modality: rna
requested_reads: 0
input_check:
  expected_files:
    - synthetic_R1.fastq (read_id synthetic_R1, filename synthetic_R1.fastq, url_basename synthetic_R1.fastq)
  supplied_inputs:
    - tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
  matched_inputs:
    - tests/fixtures/synthetic/fastqs/synthetic_R1.fastq -> synthetic_R1.fastq (read_id synthetic_R1 via file_id)
  missing_expected_files:
    (none)

file: tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
read_id: synthetic_R1
file_id: synthetic_R1.fastq
matched_by: file_id
sampled_reads: 5
region: barcode [0:4] onlist=tests/fixtures/synthetic/onlists/synthetic_barcodes.txt entries=2
  covered: 4 (0.8000)
  short_reads: 1
  onlist_matches: 3 (0.7500)
  offlist_reads: 1
  top_offlist_sequences:
    GGGG 1
"#;

    assert_eq!(observed, expected);
}

#[test]
fn test_cut_and_hist_text_synthetic_golden() {
    let cut = run_success(&[
        "cut",
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
        "-m",
        "rna",
        "-r",
        "barcode",
        "-n",
        "0",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);
    let expected_cut = r#"seqcheck cut
spec: tests/fixtures/synthetic/spec.yaml
modality: rna
requested_reads: 0
input_check:
  expected_files:
    - synthetic_R1.fastq (read_id synthetic_R1, filename synthetic_R1.fastq, url_basename synthetic_R1.fastq)
  supplied_inputs:
    - tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
  matched_inputs:
    - tests/fixtures/synthetic/fastqs/synthetic_R1.fastq -> synthetic_R1.fastq (read_id synthetic_R1 via file_id)
  missing_expected_files:
    (none)

file: tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
read_id: synthetic_R1
file_id: synthetic_R1.fastq
matched_by: file_id
region: barcode
sampled_reads: 5
covered_reads: 4
short_reads: 1
  1 ACGT
  2 TGCA
  3 ACGT
  4 GGGG
"#;
    assert_eq!(cut, expected_cut);

    let hist = run_success(&[
        "hist",
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
        "-m",
        "rna",
        "-r",
        "barcode",
        "-n",
        "0",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);
    let expected_hist = r#"seqcheck hist
spec: tests/fixtures/synthetic/spec.yaml
modality: rna
requested_reads: 0
input_check:
  expected_files:
    - synthetic_R1.fastq (read_id synthetic_R1, filename synthetic_R1.fastq, url_basename synthetic_R1.fastq)
  supplied_inputs:
    - tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
  matched_inputs:
    - tests/fixtures/synthetic/fastqs/synthetic_R1.fastq -> synthetic_R1.fastq (read_id synthetic_R1 via file_id)
  missing_expected_files:
    (none)

file: tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
read_id: synthetic_R1
file_id: synthetic_R1.fastq
matched_by: file_id
region: barcode
sampled_reads: 5
covered_reads: 4
short_reads: 1
  ACGT 2
  GGGG 1
  TGCA 1
"#;
    assert_eq!(hist, expected_hist);
}

#[test]
fn test_coverage_json_uses_current_seqspec_fixture() {
    let observed = run_success(&[
        "coverage",
        "--format",
        "json",
        "-s",
        "../seqspec/tests/fixtures/spec.yaml",
        "-m",
        "rna",
        "-n",
        "100",
        "../seqspec/tests/fixtures/fastqs/rna_R1_SRR18677638.fastq.gz",
        "../seqspec/tests/fixtures/fastqs/rna_R2_SRR18677638.fastq.gz",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();

    assert_eq!(parsed["modality"], "rna");
    assert_eq!(parsed["command"], "coverage");
    assert_eq!(
        parsed["input_check"]["missing_expected_files"]
            .as_array()
            .unwrap()
            .len(),
        0
    );
    assert_eq!(parsed["files"].as_array().unwrap().len(), 2);
    assert_eq!(parsed["files"][0]["read_id"], "rna_R1");
    assert_eq!(
        parsed["files"][0]["results"]["expected_regions"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert_eq!(
        parsed["files"][0]["results"]["full_read_coverage_count"],
        100
    );
    assert_eq!(parsed["files"][1]["read_id"], "rna_R2");
    assert_eq!(
        parsed["files"][1]["results"]["expected_regions"][0]["region_id"],
        "cdna"
    );
    assert_eq!(
        parsed["files"][1]["results"]["regions"][0]["covered_count"],
        100
    );
}

#[test]
fn test_coverage_surfaces_duplicate_assignment_warnings() {
    let observed = run_success(&[
        "coverage",
        "--format",
        "json",
        "-s",
        "tests/fixtures/bad_geometry/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/bad_geometry/fastqs/bad_R1.fastq",
        "tests/fixtures/bad_geometry/fastqs/bad_R2.fastq",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();

    let warnings = parsed["warnings"].as_array().unwrap();
    assert!(!warnings.is_empty());
    assert!(warnings.iter().any(|warning| {
        warning
            .as_str()
            .unwrap()
            .contains("region_id 'barcode' appears in multiple reads/files")
    }));
    assert!(warnings.iter().any(|warning| {
        warning
            .as_str()
            .unwrap()
            .contains("intended overlapping paired-end reads")
    }));
    assert!(warnings.iter().any(|warning| {
        warning
            .as_str()
            .unwrap()
            .contains("region_type 'barcode' appears in multiple reads/files")
    }));
    assert!(warnings.iter().any(|warning| {
        warning
            .as_str()
            .unwrap()
            .contains("region_id 'umi' appears in multiple reads/files")
    }));
}

#[test]
fn test_subset_input_check_warns_and_succeeds() {
    let observed = run_success(&[
        "length",
        "--format",
        "json",
        "-s",
        "../seqspec/tests/fixtures/spec.yaml",
        "-m",
        "rna",
        "-n",
        "10",
        "../seqspec/tests/fixtures/fastqs/rna_R1_SRR18677638.fastq.gz",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();

    assert_eq!(
        parsed["input_check"]["matched_inputs"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        parsed["input_check"]["missing_expected_files"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        parsed["input_check"]["missing_expected_files"][0]["read_id"],
        "rna_R2"
    );
    assert!(parsed["warnings"]
        .as_array()
        .unwrap()
        .iter()
        .any(|warning| warning
            .as_str()
            .unwrap()
            .contains("missing expected modality files")));
}

#[test]
fn test_duplicate_resolved_input_fails_before_scanning() {
    let stderr = run_failure(&[
        "length",
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
        "-m",
        "rna",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);

    assert!(stderr.contains("both resolved to read 'synthetic_R1'"));
}

#[test]
fn test_unmatched_input_fails_before_scanning() {
    let stderr = run_failure(&[
        "length",
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
        "-m",
        "rna",
        "tests/fixtures/synthetic/fastqs/unmatched.fastq",
    ]);

    assert!(stderr.contains("could not match 'unmatched.fastq'"));
}

#[test]
fn test_random_json_reports_sequence_entropy_against_max() {
    let observed = run_success(&[
        "random",
        "--format",
        "json",
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();

    let regions = parsed["files"][0]["results"]["regions"].as_array().unwrap();
    let umi = regions
        .iter()
        .find(|region| region["region_id"] == "umi")
        .unwrap();
    let cdna = regions
        .iter()
        .find(|region| region["region_id"] == "cdna")
        .unwrap();

    assert_eq!(umi["covered_count"], 4);
    assert_eq!(umi["short_read_count"], 1);
    assert_eq!(umi["unique_sequence_count"], 3);
    assert!((umi["sequence_entropy_bits"].as_f64().unwrap() - 1.5).abs() < 1e-9);
    assert_eq!(umi["max_entropy_bits"].as_f64().unwrap(), 4.0);
    assert!((umi["sequence_entropy_fraction"].as_f64().unwrap() - 0.375).abs() < 1e-9);
    assert_eq!(umi["top_sequences"][0]["sequence"], "AA");
    assert_eq!(umi["top_sequences"][0]["count"], 2);
    assert_eq!(umi["top_sequences"][1]["sequence"], "GG");
    assert_eq!(umi["top_sequences"][2]["sequence"], "TT");

    assert_eq!(cdna["covered_count"], 4);
    assert_eq!(cdna["unique_sequence_count"], 2);
    assert_eq!(cdna["top_sequences"][0]["sequence"], "CCCC");
    assert_eq!(cdna["top_sequences"][0]["count"], 3);
    assert_eq!(cdna["top_sequences"][1]["sequence"], "GGGG");
    assert!((cdna["sequence_entropy_bits"].as_f64().unwrap() - 0.8112781244591328).abs() < 1e-12);
    assert!(
        (cdna["sequence_entropy_fraction"].as_f64().unwrap() - 0.1014097655573916).abs() < 1e-12
    );
}

#[test]
fn test_fixed_json_matches_reverse_complement_on_negative_strand() {
    let observed = run_success(&[
        "fixed",
        "--format",
        "json",
        "-s",
        "tests/fixtures/neg_fixed/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/neg_fixed/fastqs/neg.fastq",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();
    let regions = parsed["files"][0]["results"]["regions"].as_array().unwrap();
    let target = regions
        .iter()
        .find(|region| region["region_id"] == "target")
        .unwrap();

    assert_eq!(target["exact_match_count"], 3);
    assert_eq!(
        target["top_nonmatching_sequences"]
            .as_array()
            .unwrap()
            .len(),
        0
    );
    assert_eq!(target["orientation_counts"]["forward"], 0);
    assert_eq!(target["orientation_counts"]["reverse_complement"], 3);
    let offset_histogram = target["offset_histogram"].as_array().unwrap();
    assert_eq!(offset_histogram.len(), 1);
    assert_eq!(offset_histogram[0]["offset"], 0);
    assert_eq!(offset_histogram[0]["count"], 3);
    assert_eq!(target["absent_count"], 0);
    assert_eq!(target["multi_hit_count"], 0);
}

#[test]
fn test_primer_json_reports_hits_and_classification() {
    let observed = run_success(&[
        "primer",
        "--format",
        "json",
        "-s",
        "tests/fixtures/primer_cases/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/primer_cases/fastqs/fixed.fastq",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();
    let results = &parsed["files"][0]["results"];

    assert_eq!(results["primer_classification"]["kind"], "fixed_scannable");
    assert_eq!(results["primer_classification"]["scannable"], true);
    assert_eq!(results["forward_start_hit_count"], 1);
    assert_eq!(results["forward_internal_hit_count"], 1);
    assert_eq!(results["reverse_complement_start_hit_count"], 1);
    assert_eq!(results["reverse_complement_internal_hit_count"], 1);
    assert_eq!(results["absent_count"], 1);
    assert_eq!(results["forward_hit_positions"][0]["position"], 0);
    assert_eq!(results["forward_hit_positions"][1]["position"], 2);
    assert_eq!(
        parsed["input_check"]["missing_expected_files"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn test_primer_json_reports_ghost_primer() {
    let observed = run_success(&[
        "primer",
        "--format",
        "json",
        "-s",
        "tests/fixtures/primer_cases/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/primer_cases/fastqs/ghost.fastq",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();
    let results = &parsed["files"][0]["results"];

    assert_eq!(results["primer_classification"]["kind"], "ghost_primer");
    assert_eq!(results["primer_classification"]["scannable"], false);
    assert_eq!(results["sampled_count"], 3);
}

#[test]
fn test_primer_json_reports_non_scannable_primer() {
    let observed = run_success(&[
        "primer",
        "--format",
        "json",
        "-s",
        "tests/fixtures/primer_cases/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/primer_cases/fastqs/bad.fastq",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();
    let results = &parsed["files"][0]["results"];

    assert_eq!(
        results["primer_classification"]["kind"],
        "non_scannable_primer"
    );
    assert_eq!(results["primer_classification"]["scannable"], false);
    assert!(parsed["warnings"]
        .as_array()
        .unwrap()
        .iter()
        .any(|warning| warning
            .as_str()
            .unwrap()
            .contains("non-scannable primer region 'bad_primer'")));
}

#[test]
fn test_non_rna_length_json_uses_seqspec_fixture() {
    let observed = run_success(&[
        "length",
        "--format",
        "json",
        "-s",
        "../seqspec/tests/fixtures/spec.yaml",
        "-m",
        "protein",
        "-n",
        "100",
        "../seqspec/tests/fixtures/fastqs/protein_R2_SRR18677644.fastq.gz",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();

    assert_eq!(parsed["files"][0]["read_id"], "protein_R2");
    assert_eq!(parsed["files"][0]["results"]["expected_min_len"], 15);
    assert_eq!(parsed["files"][0]["results"]["observed_max_len"], 15);
    assert_eq!(parsed["files"][0]["results"]["out_of_range_count"], 0);
}

#[test]
fn test_onlist_json_uses_upgraded_10x_index_fixture() {
    assert!(Path::new(manifest_dir())
        .join("examples/10xv3/spec.yaml")
        .exists());

    let observed = run_success(&[
        "onlist",
        "--format",
        "json",
        "-s",
        "examples/10xv3/spec.yaml",
        "-m",
        "rna",
        "-n",
        "100",
        "examples/10xv3/fastqs/I1.fastq.gz",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();

    assert_eq!(parsed["files"][0]["read_id"], "I1.fastq.gz");
    assert_eq!(
        parsed["files"][0]["results"]["regions"][0]["region_id"],
        "index7"
    );
    assert_eq!(
        parsed["files"][0]["results"]["regions"][0]["onlist_source"],
        "examples/10xv3/index7_onlist.txt"
    );
    assert_eq!(
        parsed["files"][0]["results"]["regions"][0]["covered_count"],
        100
    );
}

#[test]
fn test_version_json_reports_seqspec_version() {
    let observed = run_success(&[
        "version",
        "--format",
        "json",
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
    ]);
    let parsed: Value = serde_json::from_str(&observed).unwrap();

    assert_eq!(parsed["seqcheck_version"], env!("CARGO_PKG_VERSION"));
    assert_eq!(parsed["seqspec_file_version"], "0.4.0");
    assert_eq!(parsed["assay_id"], "synthetic");
}
