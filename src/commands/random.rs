use crate::context::{filter_regions_by_sequence_type, load_resolved_inputs};
use crate::report::{format_fraction, write_report, FileReport, ReportEnvelope};
use crate::scan::{extract_region, scan_fastq};
use crate::CommonMetricArgs;
use anyhow::Result;
use serde::Serialize;

const MAX_DNA_ENTROPY_BITS: f64 = 2.0;

#[derive(Debug, clap::Args)]
pub struct RandomArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Serialize)]
pub struct BaseCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
    pub other: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct RandomPositionResult {
    pub offset: usize,
    pub counts: BaseCounts,
    pub valid_base_count: usize,
    pub entropy_bits: f64,
    pub max_entropy_bits: f64,
    pub entropy_fraction: f64,
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
    pub mean_entropy_bits: f64,
    pub max_entropy_bits: f64,
    pub mean_entropy_fraction: f64,
    pub positions: Vec<RandomPositionResult>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RandomResult {
    pub sampled_count: usize,
    pub regions: Vec<RandomRegionResult>,
}

#[derive(Debug, Clone, Default)]
struct PositionState {
    a: usize,
    c: usize,
    g: usize,
    t: usize,
    n: usize,
    other: usize,
}

#[derive(Debug, Clone, Default)]
struct RandomRegionState {
    covered_count: usize,
    short_read_count: usize,
    positions: Vec<PositionState>,
}

pub fn run(args: &RandomArgs) -> Result<()> {
    let (_, inputs) = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
    )?;
    let mut files = Vec::new();

    for input in inputs {
        let input_path = input.input_path.clone();
        let read_id = input.read.read_id.clone();
        let file_id = input.file_id();
        let matched_by = input.matched_by.clone();
        let random_regions = filter_regions_by_sequence_type(&input, "random");
        let state = random_regions
            .iter()
            .map(|region| RandomRegionState {
                covered_count: 0,
                short_read_count: 0,
                positions: vec![
                    PositionState::default();
                    usize::try_from(region.stop - region.start).unwrap_or_default()
                ],
            })
            .collect::<Vec<_>>();
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
                            for (offset, base) in observed.chars().enumerate() {
                                increment_base(&mut state[idx].positions[offset], base);
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
        warnings: Vec::new(),
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
    let positions = state
        .positions
        .iter()
        .enumerate()
        .map(|(offset, position)| {
            let valid_base_count = position.a + position.c + position.g + position.t;
            let entropy_bits = entropy_bits(position);
            RandomPositionResult {
                offset,
                counts: BaseCounts {
                    a: position.a,
                    c: position.c,
                    g: position.g,
                    t: position.t,
                    n: position.n,
                    other: position.other,
                },
                valid_base_count,
                entropy_bits,
                max_entropy_bits: MAX_DNA_ENTROPY_BITS,
                entropy_fraction: if MAX_DNA_ENTROPY_BITS == 0.0 {
                    0.0
                } else {
                    entropy_bits / MAX_DNA_ENTROPY_BITS
                },
            }
        })
        .collect::<Vec<_>>();

    let mean_entropy_bits = if positions.is_empty() {
        0.0
    } else {
        positions
            .iter()
            .map(|position| position.entropy_bits)
            .sum::<f64>()
            / positions.len() as f64
    };

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
        mean_entropy_bits,
        max_entropy_bits: MAX_DNA_ENTROPY_BITS,
        mean_entropy_fraction: if MAX_DNA_ENTROPY_BITS == 0.0 {
            0.0
        } else {
            mean_entropy_bits / MAX_DNA_ENTROPY_BITS
        },
        positions,
    }
}

fn increment_base(position: &mut PositionState, base: char) {
    match base.to_ascii_uppercase() {
        'A' => position.a += 1,
        'C' => position.c += 1,
        'G' => position.g += 1,
        'T' => position.t += 1,
        'N' => position.n += 1,
        _ => position.other += 1,
    }
}

fn entropy_bits(position: &PositionState) -> f64 {
    let total = (position.a + position.c + position.g + position.t) as f64;
    if total == 0.0 {
        return 0.0;
    }

    [position.a, position.c, position.g, position.t]
        .into_iter()
        .filter(|count| *count > 0)
        .map(|count| {
            let probability = count as f64 / total;
            -(probability * probability.log2())
        })
        .sum()
}

fn render_text(report: &ReportEnvelope<RandomResult>) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "seqcheck random\nspec: {}\nmodality: {}\nrequested_reads: {}\n",
        report.spec.display(),
        report.modality,
        report.n_reads
    ));

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
                "region: {} [{}:{}]\n  covered: {} ({})\n  short_reads: {}\n  mean_entropy_bits: {:.4} / {:.4} ({:.4})\n",
                region.region_id,
                region.start,
                region.stop,
                region.covered_count,
                format_fraction(region.covered_count, region.sampled_count),
                region.short_read_count,
                region.mean_entropy_bits,
                region.max_entropy_bits,
                region.mean_entropy_fraction,
            ));
            out.push_str("  positions:\n");
            for position in &region.positions {
                out.push_str(&format!(
                    "    {} A={} C={} G={} T={} N={} other={} entropy={:.4} fraction={:.4}\n",
                    position.offset,
                    position.counts.a,
                    position.counts.c,
                    position.counts.g,
                    position.counts.t,
                    position.counts.n,
                    position.counts.other,
                    position.entropy_bits,
                    position.entropy_fraction,
                ));
            }
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_entropy_bits() {
        let position = PositionState {
            a: 1,
            c: 1,
            g: 1,
            t: 1,
            n: 0,
            other: 0,
        };
        assert!((entropy_bits(&position) - 2.0).abs() < 1e-9);

        let position = PositionState {
            a: 2,
            c: 0,
            g: 1,
            t: 1,
            n: 0,
            other: 0,
        };
        assert!((entropy_bits(&position) - 1.5).abs() < 1e-9);
    }
}
