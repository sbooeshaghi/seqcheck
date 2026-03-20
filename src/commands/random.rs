use crate::context::{filter_regions_by_sequence_type, load_resolved_inputs};
use crate::report::{
    format_fraction, render_report_prelude, top_sequences, write_report, FileReport,
    ReportEnvelope, SequenceCount,
};
use crate::scan::{extract_region, scan_fastq};
use crate::CommonMetricArgs;
use anyhow::Result;
use serde::Serialize;
use std::collections::HashMap;

const DNA_BITS_PER_BASE: f64 = 2.0;
const TOP_SEQUENCE_LIMIT: usize = 10;

#[derive(Debug, clap::Args)]
pub struct RandomArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Serialize)]
pub struct RandomRegionResult {
    pub region_id: String,
    pub name: String,
    pub region_type: String,
    pub start: usize,
    pub stop: usize,
    pub sampled_count: usize,
    pub covered_count: usize,
    pub covered_fraction: f64,
    pub short_read_count: usize,
    pub unique_sequence_count: usize,
    pub sequence_entropy_bits: f64,
    pub max_entropy_bits: f64,
    pub sequence_entropy_fraction: f64,
    pub top_sequences: Vec<SequenceCount>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RandomResult {
    pub sampled_count: usize,
    pub regions: Vec<RandomRegionResult>,
}

#[derive(Debug, Clone, Default)]
struct RandomRegionState {
    covered_count: usize,
    short_read_count: usize,
    sequences: HashMap<String, usize>,
}

pub fn run(args: &RandomArgs) -> Result<()> {
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
        let random_regions = filter_regions_by_sequence_type(&input, "random");
        let state = vec![RandomRegionState::default(); random_regions.len()];
        let regions = random_regions.clone();

        let (state, sampled_count) = scan_fastq(
            &input,
            args.common.n_reads,
            state,
            |state, _, sequence, _| {
                for (idx, region) in regions.iter().enumerate() {
                    match extract_region(sequence, region) {
                        Some(observed) => {
                            state[idx].covered_count += 1;
                            *state[idx].sequences.entry(observed).or_insert(0) += 1;
                        }
                        None => {
                            state[idx].short_read_count += 1;
                        }
                    }
                }
                Ok(())
            },
        )?;

        let regions = random_regions
            .iter()
            .enumerate()
            .map(|(idx, region)| build_region_result(region, sampled_count, &state[idx]))
            .collect();

        files.push(FileReport {
            input_path,
            read_id,
            file_id,
            matched_by,
            results: RandomResult {
                sampled_count,
                regions,
            },
        });
    }

    let report = ReportEnvelope {
        spec: args.common.spec.clone(),
        modality: args.common.modality.clone(),
        command: "random".to_string(),
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

fn build_region_result(
    region: &seqspec::region::RegionCoordinate,
    sampled_count: usize,
    state: &RandomRegionState,
) -> RandomRegionResult {
    let region_len = usize::try_from(region.stop - region.start).unwrap_or_default();
    let max_entropy_bits = region_len as f64 * DNA_BITS_PER_BASE;
    let sequence_entropy_bits = sequence_entropy_bits(&state.sequences);

    RandomRegionResult {
        region_id: region.region.region_id.clone(),
        name: region.region.name.clone(),
        region_type: region.region.region_type.clone(),
        start: usize::try_from(region.start).unwrap_or_default(),
        stop: usize::try_from(region.stop).unwrap_or_default(),
        sampled_count,
        covered_count: state.covered_count,
        covered_fraction: crate::report::fraction(state.covered_count, sampled_count),
        short_read_count: state.short_read_count,
        unique_sequence_count: state.sequences.len(),
        sequence_entropy_bits,
        max_entropy_bits,
        sequence_entropy_fraction: if max_entropy_bits == 0.0 {
            0.0
        } else {
            sequence_entropy_bits / max_entropy_bits
        },
        top_sequences: top_sequences(&state.sequences, TOP_SEQUENCE_LIMIT),
    }
}

fn sequence_entropy_bits(sequences: &HashMap<String, usize>) -> f64 {
    let total = sequences.values().sum::<usize>() as f64;
    if total == 0.0 {
        return 0.0;
    }

    sequences
        .values()
        .filter(|count| **count > 0)
        .map(|count| {
            let probability = *count as f64 / total;
            -(probability * probability.log2())
        })
        .sum()
}

fn render_text(report: &ReportEnvelope<RandomResult>) -> String {
    let mut out = render_report_prelude("random", report);

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
                "region: {} [{}:{}]\n  covered: {} ({})\n  short_reads: {}\n  unique_sequences: {}\n  sequence_entropy_bits: {:.4} / {:.4} ({:.4})\n",
                region.region_id,
                region.start,
                region.stop,
                region.covered_count,
                format_fraction(region.covered_count, region.sampled_count),
                region.short_read_count,
                region.unique_sequence_count,
                region.sequence_entropy_bits,
                region.max_entropy_bits,
                region.sequence_entropy_fraction,
            ));
            if !region.top_sequences.is_empty() {
                out.push_str("  top_sequences:\n");
                for entry in &region.top_sequences {
                    out.push_str(&format!("    {} {}\n", entry.sequence, entry.count));
                }
            }
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_entropy_bits_uses_exact_sequence_distribution() {
        let mut sequences = HashMap::new();
        sequences.insert("AAA".to_string(), 2);
        sequences.insert("AAT".to_string(), 1);

        let entropy = sequence_entropy_bits(&sequences);

        assert!((entropy - 0.9182958340544896).abs() < 1e-12);
    }

    #[test]
    fn test_sequence_entropy_bits_is_zero_for_empty_distribution() {
        let entropy = sequence_entropy_bits(&HashMap::new());
        assert_eq!(entropy, 0.0);
    }
}
