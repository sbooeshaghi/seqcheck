use serde_json::Value;
use std::fs;
use std::path::Path;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

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

fn parse_json(args: &[&str]) -> Value {
    serde_json::from_str(&run_success(args)).unwrap()
}

fn find_result<F>(report: &Value, predicate: F) -> &Value
where
    F: Fn(&Value) -> bool,
{
    report["results"]
        .as_array()
        .unwrap()
        .iter()
        .find(|result| predicate(result))
        .unwrap()
}

fn find_metric<'a>(result: &'a Value, section: &str, name: &str) -> &'a Value {
    result[section]
        .as_array()
        .unwrap()
        .iter()
        .find(|item| item["name"] == name)
        .unwrap()
}

fn metric_value<'a>(result: &'a Value, section: &str, name: &str) -> &'a Value {
    &find_metric(result, section, name)["data"]["value"]
}

fn has_assessment(result: &Value, code: &str) -> bool {
    result["assessment"]
        .as_array()
        .unwrap_or(&Vec::new())
        .iter()
        .any(|assessment| assessment["code"] == code)
}

fn temp_path(extension: &str) -> String {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir()
        .join(format!(
            "seqcheck-cli-{}-{}.{}",
            std::process::id(),
            nanos,
            extension
        ))
        .display()
        .to_string()
}

#[test]
fn test_length_json_uses_atomic_schema() {
    let parsed = parse_json(&[
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

    assert_eq!(parsed["report_schema_version"], "0.1.0");
    assert_eq!(parsed["meta"]["command"], "length");
    assert_eq!(parsed["meta"]["modality"], "rna");

    let input_check = find_result(&parsed, |result| result["check"] == "input_check");
    assert_eq!(input_check["files"][0], "synthetic_R1.fastq");
    assert!(has_assessment(input_check, "all_expected_files_matched"));

    let length = find_result(&parsed, |result| {
        result["check"] == "length"
            && result["files"] == serde_json::json!(["synthetic_R1.fastq"])
            && result["regions"].as_array().unwrap().is_empty()
    });
    assert_eq!(metric_value(length, "expected", "expected_min_len"), 12);
    assert_eq!(metric_value(length, "expected", "expected_max_len"), 12);
    assert_eq!(metric_value(length, "observed", "sampled_count"), 5);
    assert_eq!(metric_value(length, "observed", "observed_min_len"), 2);
    assert_eq!(metric_value(length, "observed", "observed_max_len"), 12);
    assert_eq!(metric_value(length, "observed", "out_of_range_count"), 1);
    assert_eq!(
        metric_value(length, "observed", "out_of_range_fraction"),
        0.2
    );
    assert!(has_assessment(length, "length_out_of_range"));
}

#[test]
fn test_check_json_aggregates_core_checks() {
    let parsed = parse_json(&[
        "check",
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

    assert_eq!(parsed["report_schema_version"], "0.1.0");
    assert_eq!(parsed["meta"]["command"], "check");
    assert_eq!(parsed["results"][0]["check"], "input_check");

    let checks = parsed["results"]
        .as_array()
        .unwrap()
        .iter()
        .map(|result| result["check"].as_str().unwrap())
        .collect::<Vec<_>>();

    assert!(checks.contains(&"input_check"));
    assert!(checks.contains(&"length"));
    assert!(checks.contains(&"coverage"));
    assert!(checks.contains(&"primer"));
    assert!(checks.contains(&"fixed"));
    assert!(checks.contains(&"onlist"));
    assert!(checks.contains(&"random"));
    assert!(!checks.contains(&"cut"));
    assert!(!checks.contains(&"hist"));
}

#[test]
fn test_check_json_accepts_ontology_region_type_lists() {
    let parsed = parse_json(&[
        "check",
        "--format",
        "json",
        "-s",
        "tests/fixtures/synthetic/spec_0_5.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);

    let barcode = find_result(&parsed, |result| {
        result["check"] == "coverage" && result["regions"] == serde_json::json!(["barcode"])
    });
    assert_eq!(
        metric_value(barcode, "expected", "expected_region")[0]["region_type"],
        serde_json::json!(["RGN:partition:cell"])
    );

    let primer = find_result(&parsed, |result| {
        result["check"] == "primer" && result["regions"] == serde_json::json!(["primer"])
    });
    assert_eq!(
        metric_value(primer, "expected", "primer_region_type"),
        &serde_json::json!(["RGN:technical:primer"])
    );

    let onlist = find_result(&parsed, |result| {
        result["check"] == "onlist" && result["regions"] == serde_json::json!(["barcode"])
    });
    assert_eq!(metric_value(onlist, "observed", "covered_count"), 4);
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

check: input_check
files: synthetic_R1.fastq
reads: synthetic_R1
regions: (none)
assessment:
  - [pass] all_expected_files_matched: All expected modality files were supplied and matched uniquely.
expected:
  - expected_file_mappings:
    - file_id=synthetic_R1.fastq filename=synthetic_R1.fastq read_id=synthetic_R1 url_basename=synthetic_R1.fastq
observed:
  - supplied_input_paths:
    - input_path=tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
  - matched_inputs:
    - file_id=synthetic_R1.fastq input_path=tests/fixtures/synthetic/fastqs/synthetic_R1.fastq matched_by=file_id read_id=synthetic_R1
  - missing_expected_file_ids:
    (none)
  - primer_classifications:
    - classification=fixed_scannable file_id=synthetic_R1.fastq input_path=tests/fixtures/synthetic/fastqs/synthetic_R1.fastq primer_id=primer primer_region_id=primer read_id=synthetic_R1 reason=null scannable=true

check: fixed
files: synthetic_R1.fastq
reads: synthetic_R1
regions: linker
assessment:
  - [interpretation] fixed_exact_match_partial: Fixed region 'linker' matches exactly in some covered reads, with the dominant orientation 'forward'.
expected:
  - expected_sequence: TT bases
  - expected_coordinates:
    - name=Linker region_type=linker start=4 stop=6
  - primary_orientation: forward
observed:
  - sampled_count: 5 count
  - covered_count: 4 count
  - covered_fraction: 0.8 fraction
  - short_read_count: 1 count
  - exact_match_count: 3 count
  - exact_match_fraction: 0.75 fraction
  - orientation_counts:
    - count=3 orientation=forward
    - count=3 orientation=reverse
    - count=0 orientation=complement
    - count=0 orientation=reverse_complement
  - offset_histogram:
    -1 1 count
    +0 3 count
    +2 1 count
  - absent_count: 1 count
  - multi_hit_count: 1 count
  - top_nonmatching_sequences:
    - count=1 sequence=AC
"#;

    assert_eq!(observed, expected);
}

#[test]
fn test_onlist_json_synthetic_reports_offlist_sequences() {
    let parsed = parse_json(&[
        "onlist",
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

    let onlist = find_result(&parsed, |result| {
        result["check"] == "onlist" && result["regions"] == serde_json::json!(["barcode"])
    });
    assert_eq!(
        metric_value(onlist, "expected", "onlist_source"),
        "tests/fixtures/synthetic/onlists/synthetic_barcodes.txt"
    );
    assert_eq!(metric_value(onlist, "observed", "covered_count"), 4);
    assert_eq!(metric_value(onlist, "observed", "short_read_count"), 1);
    assert_eq!(metric_value(onlist, "observed", "exact_onlist_count"), 3);
    assert_eq!(metric_value(onlist, "observed", "offlist_count"), 1);
    assert!(has_assessment(onlist, "offlist_sequences_detected"));
}

#[test]
fn test_cut_and_hist_json_report_atomic_results() {
    let cut = parse_json(&[
        "cut",
        "--format",
        "json",
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
    let cut_result = find_result(&cut, |result| {
        result["check"] == "cut" && result["regions"] == serde_json::json!(["barcode"])
    });
    assert_eq!(metric_value(cut_result, "observed", "covered_count"), 4);
    assert_eq!(metric_value(cut_result, "observed", "short_read_count"), 1);
    assert_eq!(
        metric_value(cut_result, "observed", "extracted_sequences")[0]["sequence"],
        "ACGT"
    );

    let hist = parse_json(&[
        "hist",
        "--format",
        "json",
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
    let hist_result = find_result(&hist, |result| {
        result["check"] == "hist" && result["regions"] == serde_json::json!(["barcode"])
    });
    assert_eq!(metric_value(hist_result, "observed", "covered_count"), 4);
    assert_eq!(
        metric_value(hist_result, "observed", "sequence_histogram")[0]["sequence"],
        "ACGT"
    );
    assert_eq!(
        metric_value(hist_result, "observed", "sequence_histogram")[0]["count"],
        2
    );
}

#[test]
fn test_coverage_json_uses_current_seqspec_fixture() {
    let parsed = parse_json(&[
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

    assert_eq!(parsed["meta"]["modality"], "rna");
    assert_eq!(parsed["meta"]["command"], "coverage");

    let input_check = find_result(&parsed, |result| result["check"] == "input_check");
    assert!(has_assessment(input_check, "all_expected_files_matched"));

    let r1_file = find_result(&parsed, |result| {
        result["check"] == "coverage"
            && result["files"] == serde_json::json!(["rna_R1_SRR18677638.fastq.gz"])
            && result["regions"].as_array().unwrap().is_empty()
    });
    assert_eq!(
        metric_value(r1_file, "observed", "full_read_coverage_count"),
        100
    );
    assert_eq!(
        metric_value(r1_file, "expected", "expected_regions")
            .as_array()
            .unwrap()
            .len(),
        2
    );

    let cdna_region = find_result(&parsed, |result| {
        result["check"] == "coverage"
            && result["reads"] == serde_json::json!(["rna_R2"])
            && result["regions"] == serde_json::json!(["cdna"])
    });
    assert_eq!(metric_value(cdna_region, "observed", "covered_count"), 100);
}

#[test]
fn test_coverage_surfaces_duplicate_assignment_as_results() {
    let parsed = parse_json(&[
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

    let shared_barcode = find_result(&parsed, |result| {
        result["check"] == "coverage"
            && result["regions"] == serde_json::json!(["barcode"])
            && has_assessment(result, "shared_region_id_visible_in_multiple_reads")
    });
    assert!(shared_barcode["assessment"][0]["description"]
        .as_str()
        .unwrap()
        .contains("intended paired-end overlap"));

    let shared_region_type = find_result(&parsed, |result| {
        result["check"] == "coverage"
            && result["regions"].as_array().unwrap().is_empty()
            && has_assessment(result, "shared_region_type_visible_in_multiple_reads")
            && metric_value(result, "expected", "projected_coordinates_by_read")[0]["region_type"]
                == "barcode"
    });
    assert!(has_assessment(
        shared_region_type,
        "shared_region_type_visible_in_multiple_reads"
    ));
    assert_eq!(
        metric_value(
            shared_region_type,
            "expected",
            "shared_region_type_term"
        ),
        "RGN:partition:cell"
    );
}

#[test]
fn test_coverage_groups_ontology_terms_and_preserves_multiple_roles() {
    let parsed = parse_json(&[
        "coverage",
        "--format",
        "json",
        "-s",
        "tests/fixtures/bad_geometry/spec_0_5.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/bad_geometry/fastqs/bad_R1.fastq",
        "tests/fixtures/bad_geometry/fastqs/bad_R2.fastq",
    ]);

    let shared_cell_partition = find_result(&parsed, |result| {
        result["check"] == "coverage"
            && result["regions"].as_array().unwrap().is_empty()
            && has_assessment(result, "shared_region_type_visible_in_multiple_reads")
            && metric_value(result, "expected", "projected_coordinates_by_read")[0]["region_type"]
                == serde_json::json!(["RGN:partition:cell"])
    });
    assert!(shared_cell_partition["assessment"][0]["description"]
        .as_str()
        .unwrap()
        .contains("RGN:partition:cell"));
    assert_eq!(
        metric_value(
            shared_cell_partition,
            "expected",
            "shared_region_type_term"
        ),
        "RGN:partition:cell"
    );

    let guide = find_result(&parsed, |result| {
        result["check"] == "coverage" && result["regions"] == serde_json::json!(["feature"])
    });
    assert_eq!(
        metric_value(guide, "expected", "expected_region")[0]["region_type"],
        serde_json::json!(["RGN:measure:guide", "RGN:classify:perturbation"])
    );
}

#[test]
fn test_subset_input_check_warns_and_succeeds() {
    let parsed = parse_json(&[
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

    let input_check = find_result(&parsed, |result| result["check"] == "input_check");
    assert_eq!(
        metric_value(input_check, "observed", "matched_inputs")
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        metric_value(input_check, "observed", "missing_expected_file_ids")
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        metric_value(input_check, "observed", "missing_expected_file_ids")[0]["read_id"],
        "rna_R2"
    );
    assert!(has_assessment(input_check, "missing_expected_files"));
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
    let parsed = parse_json(&[
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

    let umi = find_result(&parsed, |result| {
        result["check"] == "random" && result["regions"] == serde_json::json!(["umi"])
    });
    let cdna = find_result(&parsed, |result| {
        result["check"] == "random" && result["regions"] == serde_json::json!(["cdna"])
    });

    assert_eq!(metric_value(umi, "observed", "covered_count"), 4);
    assert_eq!(metric_value(umi, "observed", "short_read_count"), 1);
    assert_eq!(metric_value(umi, "observed", "unique_sequence_count"), 3);
    assert!(
        (metric_value(umi, "observed", "sequence_entropy_bits")
            .as_f64()
            .unwrap()
            - 1.5)
            .abs()
            < 1e-9
    );
    assert_eq!(
        metric_value(umi, "expected", "max_entropy_bits")
            .as_f64()
            .unwrap(),
        4.0
    );
    assert!(
        (metric_value(umi, "observed", "sequence_entropy_fraction")
            .as_f64()
            .unwrap()
            - 0.375)
            .abs()
            < 1e-9
    );
    assert_eq!(
        metric_value(umi, "observed", "top_sequences")[0]["sequence"],
        "AA"
    );

    assert_eq!(metric_value(cdna, "observed", "covered_count"), 4);
    assert_eq!(metric_value(cdna, "observed", "unique_sequence_count"), 2);
    assert_eq!(
        metric_value(cdna, "observed", "top_sequences")[0]["sequence"],
        "CCCC"
    );
    assert!(
        (metric_value(cdna, "observed", "sequence_entropy_bits")
            .as_f64()
            .unwrap()
            - 0.8112781244591328)
            .abs()
            < 1e-12
    );
}

#[test]
fn test_fixed_json_matches_reverse_complement_on_negative_strand() {
    let parsed = parse_json(&[
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
    let target = find_result(&parsed, |result| {
        result["check"] == "fixed" && result["regions"] == serde_json::json!(["target"])
    });

    assert_eq!(metric_value(target, "observed", "exact_match_count"), 3);
    assert_eq!(
        metric_value(target, "observed", "top_nonmatching_sequences")
            .as_array()
            .unwrap()
            .len(),
        0
    );
    assert_eq!(
        metric_value(target, "observed", "orientation_counts")[0]["count"],
        0
    );
    assert_eq!(
        metric_value(target, "observed", "orientation_counts")[3]["count"],
        3
    );
    let offset_histogram = metric_value(target, "observed", "offset_histogram")
        .as_array()
        .unwrap();
    assert_eq!(offset_histogram.len(), 1);
    assert_eq!(offset_histogram[0]["key"], "+0");
    assert_eq!(offset_histogram[0]["value"], 3);
    assert_eq!(metric_value(target, "observed", "absent_count"), 0);
    assert_eq!(metric_value(target, "observed", "multi_hit_count"), 0);
}

#[test]
fn test_primer_json_reports_hits_and_classification() {
    let parsed = parse_json(&[
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
    let result = find_result(&parsed, |result| {
        result["check"] == "primer" && result["files"] == serde_json::json!(["fixed.fastq"])
    });

    assert_eq!(
        metric_value(result, "observed", "primer_classification_kind"),
        "fixed_scannable"
    );
    assert_eq!(metric_value(result, "observed", "primer_scannable"), true);
    assert_eq!(
        metric_value(result, "observed", "forward_start_hit_count"),
        1
    );
    assert_eq!(
        metric_value(result, "observed", "forward_start_hit_fraction"),
        0.2
    );
    assert_eq!(
        metric_value(result, "observed", "forward_internal_hit_count"),
        1
    );
    assert_eq!(
        metric_value(result, "observed", "reverse_complement_start_hit_count"),
        1
    );
    assert_eq!(
        metric_value(result, "observed", "reverse_complement_internal_hit_count"),
        1
    );
    assert_eq!(metric_value(result, "observed", "absent_count"), 1);
    assert_eq!(metric_value(result, "observed", "absent_fraction"), 0.2);
    assert_eq!(
        metric_value(result, "observed", "forward_hit_positions")[0]["position"],
        0
    );
    assert_eq!(
        metric_value(result, "observed", "forward_hit_positions")[1]["position"],
        2
    );

    let input_check = find_result(&parsed, |result| result["check"] == "input_check");
    assert_eq!(
        metric_value(input_check, "observed", "missing_expected_file_ids")
            .as_array()
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn test_primer_json_reports_ghost_primer() {
    let parsed = parse_json(&[
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
    let result = find_result(&parsed, |result| {
        result["check"] == "primer" && result["files"] == serde_json::json!(["ghost.fastq"])
    });

    assert_eq!(
        metric_value(result, "observed", "primer_classification_kind"),
        "ghost_primer"
    );
    assert_eq!(metric_value(result, "observed", "primer_scannable"), false);
    assert_eq!(metric_value(result, "observed", "sampled_count"), 3);
    assert!(has_assessment(result, "ghost_primer_anchor"));
}

#[test]
fn test_primer_json_reports_non_scannable_primer() {
    let parsed = parse_json(&[
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
    let result = find_result(&parsed, |result| {
        result["check"] == "primer" && result["files"] == serde_json::json!(["bad.fastq"])
    });

    assert_eq!(
        metric_value(result, "observed", "primer_classification_kind"),
        "non_scannable_primer"
    );
    assert_eq!(metric_value(result, "observed", "primer_scannable"), false);
    assert!(has_assessment(result, "non_scannable_primer"));

    let input_check = find_result(&parsed, |result| result["check"] == "input_check");
    assert!(has_assessment(input_check, "non_scannable_primer"));
}

#[test]
fn test_non_rna_length_json_uses_seqspec_fixture() {
    let parsed = parse_json(&[
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

    let result = find_result(&parsed, |result| {
        result["check"] == "length"
            && result["files"] == serde_json::json!(["protein_R2_SRR18677644.fastq.gz"])
    });
    assert_eq!(result["reads"][0], "protein_R2");
    assert_eq!(metric_value(result, "expected", "expected_min_len"), 15);
    assert_eq!(metric_value(result, "observed", "observed_max_len"), 15);
    assert_eq!(metric_value(result, "observed", "out_of_range_count"), 0);
}

#[test]
fn test_onlist_json_uses_upgraded_10x_index_fixture() {
    assert!(Path::new(manifest_dir())
        .join("examples/10xv3/spec.yaml")
        .exists());

    let parsed = parse_json(&[
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

    let result = find_result(&parsed, |result| {
        result["check"] == "onlist" && result["regions"] == serde_json::json!(["index7"])
    });
    assert_eq!(result["reads"][0], "I1.fastq.gz");
    assert_eq!(
        metric_value(result, "expected", "onlist_source"),
        "examples/10xv3/index7_onlist.txt"
    );
    assert_eq!(metric_value(result, "observed", "covered_count"), 100);
}

#[test]
fn test_onlist_missing_resource_is_structured_error() {
    let parsed = parse_json(&[
        "onlist",
        "--format",
        "json",
        "-s",
        "tests/fixtures/broken_onlist/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);

    let result = find_result(&parsed, |result| {
        result["check"] == "onlist" && result["regions"] == serde_json::json!(["barcode"])
    });
    assert_eq!(
        metric_value(result, "observed", "fetch_load_status"),
        "error"
    );
    assert!(has_assessment(result, "missing_onlist_resource"));
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

#[test]
fn test_report_command_writes_html_with_embedded_payloads() {
    let report_json = temp_path("json");
    let report_html = temp_path("html");

    run_success(&[
        "check",
        "--format",
        "json",
        "-o",
        &report_json,
        "-s",
        "tests/fixtures/synthetic/spec.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/synthetic/fastqs/synthetic_R1.fastq",
    ]);

    run_success(&["report", "-i", &report_json, "-o", &report_html]);

    let html = fs::read_to_string(&report_html).unwrap();
    assert!(html.contains("seqcheck report"));
    assert!(html.contains("id=\"seqcheck-report-data\""));
    assert!(html.contains("id=\"seqspec-lib-data\""));
    assert!(html.contains("\"assay_name\": \"Synthetic Seqcheck Fixture\""));
    assert!(html.contains("\"modality\": \"rna\""));

    let _ = fs::remove_file(report_json);
    let _ = fs::remove_file(report_html);
}

#[test]
fn test_report_command_displays_multiple_ontology_terms() {
    let report_json = temp_path("json");
    let report_html = temp_path("html");

    run_success(&[
        "coverage",
        "--format",
        "json",
        "-o",
        &report_json,
        "-s",
        "tests/fixtures/bad_geometry/spec_0_5.yaml",
        "-m",
        "rna",
        "-n",
        "0",
        "tests/fixtures/bad_geometry/fastqs/bad_R1.fastq",
        "tests/fixtures/bad_geometry/fastqs/bad_R2.fastq",
    ]);
    run_success(&["report", "-i", &report_json, "-o", &report_html]);

    let html = fs::read_to_string(&report_html).unwrap();
    assert!(html.contains("\"region_type\": \"RGN:measure:guide+RGN:classify:perturbation\""));

    let _ = fs::remove_file(report_json);
    let _ = fs::remove_file(report_html);
}

#[test]
fn test_report_command_generates_html_without_library_payload_when_spec_missing() {
    let report_json = temp_path("json");
    let report_html = temp_path("html");

    let mut parsed = parse_json(&[
        "check",
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
    parsed["meta"]["spec"] = Value::String("tests/fixtures/does_not_exist.yaml".to_string());
    fs::write(&report_json, serde_json::to_string_pretty(&parsed).unwrap()).unwrap();

    run_success(&["report", "-i", &report_json, "-o", &report_html]);

    let html = fs::read_to_string(&report_html).unwrap();
    assert!(html.contains("id=\"seqspec-lib-data\""));
    assert!(html.contains(">null</script>") || html.contains(">\nnull\n</script>"));

    let _ = fs::remove_file(report_json);
    let _ = fs::remove_file(report_html);
}

#[test]
fn test_report_command_rejects_invalid_report_json() {
    let report_json = temp_path("json");
    let report_html = temp_path("html");

    fs::write(&report_json, "{ not valid json").unwrap();

    let stderr = run_failure(&["report", "-i", &report_json, "-o", &report_html]);
    assert!(stderr.contains("failed to parse JSON report"));

    let _ = fs::remove_file(report_json);
    let _ = fs::remove_file(report_html);
}
