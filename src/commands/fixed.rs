use crate::context::{filter_regions_by_sequence_type, load_resolved_inputs};
use crate::report::{
    format_fraction, render_report_prelude, top_sequences, write_report, FileReport,
    ReportEnvelope, SequenceCount,
};
use crate::scan::{extract_region, scan_fastq};
use crate::sequence::{complement_sequence, reverse_complement_sequence, reverse_sequence};
use crate::CommonMetricArgs;
use anyhow::Result;
use serde::Serialize;
use std::collections::HashMap;

const TOP_SEQUENCE_LIMIT: usize = 10;

#[derive(Debug, clap::Args)]
pub struct FixedArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct OrientationCounts {
    pub forward: usize,
    pub reverse: usize,
    pub complement: usize,
    pub reverse_complement: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct FixedRegionResult {
    pub region_id: String,
    pub name: String,
    pub region_type: String,
    pub start: usize,
    pub stop: usize,
    pub expected_sequence: String,
    pub sampled_count: usize,
    pub covered_count: usize,
    pub covered_fraction: f64,
    pub short_read_count: usize,
    pub exact_match_count: usize,
    pub exact_match_fraction: f64,
    pub orientation_counts: OrientationCounts,
    pub top_nonmatching_sequences: Vec<SequenceCount>,
}

#[derive(Debug, Clone, Serialize)]
pub struct FixedResult {
    pub sampled_count: usize,
    pub regions: Vec<FixedRegionResult>,
}

#[derive(Debug, Clone, Default)]
struct FixedRegionState {
    covered_count: usize,
    short_read_count: usize,
    exact_match_count: usize,
    orientation_counts: OrientationCounts,
    mismatches: HashMap<String, usize>,
}

pub fn run(args: &FixedArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
    )?;
    let crate::context::LoadedInputs {
        input_check,
        warnings,
        inputs,
        ..
    } = loaded;
    let mut files = Vec::new();

    for input in inputs {
        let input_path = input.input_path.clone();
        let read_id = input.read.read_id.clone();
        let file_id = input.file_id();
        let matched_by = input.matched_by.clone();
        let primary_orientation = if input.read.strand == "neg" {
            ExpectedOrientation::ReverseComplement
        } else {
            ExpectedOrientation::Forward
        };
        let fixed_regions = filter_regions_by_sequence_type(&input, "fixed");
        let state = vec![FixedRegionState::default(); fixed_regions.len()];
        let regions = fixed_regions.clone();

        let (state, sampled_count) = scan_fastq(
            &input,
            args.common.n_reads,
            state,
            |state, _, sequence, _| {
                for (idx, region) in regions.iter().enumerate() {
                    match extract_region(sequence, region) {
                        Some(observed) => {
                            state[idx].covered_count += 1;
                            let expected = expected_sequences(&region.region.sequence);
                            update_orientation_counts(
                                &mut state[idx].orientation_counts,
                                &observed,
                                &expected,
                            );

                            let primary_expected = primary_orientation.sequence(&expected);
                            if observed == primary_expected {
                                state[idx].exact_match_count += 1;
                            } else {
                                *state[idx].mismatches.entry(observed).or_insert(0) += 1;
                            }
                        }
                        None => {
                            state[idx].short_read_count += 1;
                        }
                    }
                }
                Ok(())
            },
        )?;

        let results = fixed_regions
            .iter()
            .enumerate()
            .map(|(idx, region)| FixedRegionResult {
                region_id: region.region.region_id.clone(),
                name: region.region.name.clone(),
                region_type: region.region.region_type.clone(),
                start: usize::try_from(region.start).unwrap_or_default(),
                stop: usize::try_from(region.stop).unwrap_or_default(),
                expected_sequence: region.region.sequence.clone(),
                sampled_count,
                covered_count: state[idx].covered_count,
                covered_fraction: crate::report::fraction(state[idx].covered_count, sampled_count),
                short_read_count: state[idx].short_read_count,
                exact_match_count: state[idx].exact_match_count,
                exact_match_fraction: crate::report::fraction(
                    state[idx].exact_match_count,
                    state[idx].covered_count,
                ),
                orientation_counts: state[idx].orientation_counts.clone(),
                top_nonmatching_sequences: top_sequences(
                    &state[idx].mismatches,
                    TOP_SEQUENCE_LIMIT,
                ),
            })
            .collect();

        files.push(FileReport {
            input_path,
            read_id,
            file_id,
            matched_by,
            results: FixedResult {
                sampled_count,
                regions: results,
            },
        });
    }

    let report = ReportEnvelope {
        spec: args.common.spec.clone(),
        modality: args.common.modality.clone(),
        command: "fixed".to_string(),
        n_reads: args.common.n_reads,
        input_check,
        warnings,
        files,
    };

    write_report(
        &args.common.output,
        args.common.format,
        &report,
        render_text(&report),
    )
}

fn render_text(report: &ReportEnvelope<FixedResult>) -> String {
    let mut out = render_report_prelude("fixed", report);

    for file in &report.files {
        out.push_str(&format!(
            "\nfile: {}\nread_id: {}\nfile_id: {}\nmatched_by: {}\nsampled_reads: {}\n",
            file.input_path.display(),
            file.read_id,
            file.file_id,
            file.matched_by,
            file.results.sampled_count,
        ));
        for region in &file.results.regions {
            out.push_str(&format!(
                "region: {} [{}:{}] expected={}\n  covered: {} ({})\n  short_reads: {}\n  exact_matches: {} ({})\n",
                region.region_id,
                region.start,
                region.stop,
                region.expected_sequence,
                region.covered_count,
                format_fraction(region.covered_count, region.sampled_count),
                region.short_read_count,
                region.exact_match_count,
                format_fraction(region.exact_match_count, region.covered_count),
            ));
            out.push_str(&format!(
                "  orientation_matches:\n    forward: {}\n    reverse: {}\n    complement: {}\n    reverse_complement: {}\n",
                region.orientation_counts.forward,
                region.orientation_counts.reverse,
                region.orientation_counts.complement,
                region.orientation_counts.reverse_complement,
            ));
            if !region.top_nonmatching_sequences.is_empty() {
                out.push_str("  top_nonmatching_sequences:\n");
                for entry in &region.top_nonmatching_sequences {
                    out.push_str(&format!("    {} {}\n", entry.sequence, entry.count));
                }
            }
        }
    }

    out
}

#[derive(Debug, Clone)]
struct ExpectedSequences {
    forward: String,
    reverse: String,
    complement: String,
    reverse_complement: String,
}

#[derive(Debug, Clone, Copy)]
enum ExpectedOrientation {
    Forward,
    ReverseComplement,
}

impl ExpectedOrientation {
    fn sequence<'a>(&self, expected: &'a ExpectedSequences) -> &'a str {
        match self {
            ExpectedOrientation::Forward => &expected.forward,
            ExpectedOrientation::ReverseComplement => &expected.reverse_complement,
        }
    }
}

fn expected_sequences(sequence: &str) -> ExpectedSequences {
    ExpectedSequences {
        forward: sequence.to_string(),
        reverse: reverse_sequence(sequence),
        complement: complement_sequence(sequence),
        reverse_complement: reverse_complement_sequence(sequence),
    }
}

fn update_orientation_counts(
    counts: &mut OrientationCounts,
    observed: &str,
    expected: &ExpectedSequences,
) {
    if observed == expected.forward {
        counts.forward += 1;
    }
    if observed == expected.reverse {
        counts.reverse += 1;
    }
    if observed == expected.complement {
        counts.complement += 1;
    }
    if observed == expected.reverse_complement {
        counts.reverse_complement += 1;
    }
}
