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

file: tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
read_id: synthetic_R1
file_id: synthetic_R1.fastq
matched_by: file_id
sampled_reads: 5
region: linker [4:6] expected=TT
  covered: 4 (0.8000)
  short_reads: 1
  exact_matches: 3 (0.7500)
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
            .contains("region_id 'barcode' is assigned to multiple reads/files")
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
            .contains("region_id 'umi' is assigned to multiple reads/files")
    }));
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
