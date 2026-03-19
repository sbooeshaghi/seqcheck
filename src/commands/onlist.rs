use crate::context::{filter_regions_by_sequence_type, load_onlist, load_resolved_inputs};
use crate::report::{
    format_fraction, top_sequences, write_report, FileReport, ReportEnvelope, SequenceCount,
};
use crate::scan::{extract_region, scan_fastq};
use crate::CommonMetricArgs;
use anyhow::Result;
use serde::Serialize;
use std::collections::HashMap;

const TOP_SEQUENCE_LIMIT: usize = 10;

#[derive(Debug, clap::Args)]
pub struct OnlistArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Serialize)]
pub struct OnlistRegionResult {
    pub region_id: String,
    pub name: String,
    pub region_type: String,
    pub start: usize,
    pub stop: usize,
    pub onlist_source: String,
    pub onlist_size: usize,
    pub sampled_count: usize,
    pub covered_count: usize,
    pub covered_fraction: f64,
    pub short_read_count: usize,
    pub exact_onlist_count: usize,
    pub exact_onlist_fraction: f64,
    pub offlist_count: usize,
    pub top_offlist_sequences: Vec<SequenceCount>,
}

#[derive(Debug, Clone, Serialize)]
pub struct OnlistResult {
    pub sampled_count: usize,
    pub regions: Vec<OnlistRegionResult>,
}

#[derive(Debug, Clone, Default)]
struct OnlistRegionState {
    covered_count: usize,
    short_read_count: usize,
    exact_onlist_count: usize,
    offlist_count: usize,
    offlist_sequences: HashMap<String, usize>,
}

pub fn run(args: &OnlistArgs) -> Result<()> {
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
        let onlist_regions = filter_regions_by_sequence_type(&input, "onlist");
        let loaded_onlists = onlist_regions
            .iter()
            .map(|region| load_onlist(&input.spec_base, &region.region))
            .collect::<Result<Vec<_>>>()?;
        let state = vec![OnlistRegionState::default(); onlist_regions.len()];
        let regions = onlist_regions.clone();

        let (state, sampled_count) = scan_fastq(
            &input,
            args.common.n_reads,
            state,
            |state, _, sequence, _| {
                for (idx, region) in regions.iter().enumerate() {
                    match extract_region(sequence, region) {
                        Some(observed) => {
                            state[idx].covered_count += 1;
                            if loaded_onlists[idx].entries.contains(&observed) {
                                state[idx].exact_onlist_count += 1;
                            } else {
                                state[idx].offlist_count += 1;
                                *state[idx].offlist_sequences.entry(observed).or_insert(0) += 1;
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

        let results = onlist_regions
            .iter()
            .enumerate()
            .map(|(idx, region)| OnlistRegionResult {
                region_id: region.region.region_id.clone(),
                name: region.region.name.clone(),
                region_type: region.region.region_type.clone(),
                start: usize::try_from(region.start).unwrap_or_default(),
                stop: usize::try_from(region.stop).unwrap_or_default(),
                onlist_source: loaded_onlists[idx].source.clone(),
                onlist_size: loaded_onlists[idx].entries.len(),
                sampled_count,
                covered_count: state[idx].covered_count,
                covered_fraction: crate::report::fraction(state[idx].covered_count, sampled_count),
                short_read_count: state[idx].short_read_count,
                exact_onlist_count: state[idx].exact_onlist_count,
                exact_onlist_fraction: crate::report::fraction(
                    state[idx].exact_onlist_count,
                    state[idx].covered_count,
                ),
                offlist_count: state[idx].offlist_count,
                top_offlist_sequences: top_sequences(
                    &state[idx].offlist_sequences,
                    TOP_SEQUENCE_LIMIT,
                ),
            })
            .collect();

        files.push(FileReport {
            input_path,
            read_id,
            file_id,
            matched_by,
            results: OnlistResult {
                sampled_count,
                regions: results,
            },
        });
    }

    let report = ReportEnvelope {
        spec: args.common.spec.clone(),
        modality: args.common.modality.clone(),
        command: "onlist".to_string(),
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

fn render_text(report: &ReportEnvelope<OnlistResult>) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "seqcheck onlist\nspec: {}\nmodality: {}\nrequested_reads: {}\n",
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
                "region: {} [{}:{}] onlist={} entries={}\n  covered: {} ({})\n  short_reads: {}\n  onlist_matches: {} ({})\n  offlist_reads: {}\n",
                region.region_id,
                region.start,
                region.stop,
                region.onlist_source,
                region.onlist_size,
                region.covered_count,
                format_fraction(region.covered_count, region.sampled_count),
                region.short_read_count,
                region.exact_onlist_count,
                format_fraction(region.exact_onlist_count, region.covered_count),
                region.offlist_count,
            ));
            if !region.top_offlist_sequences.is_empty() {
                out.push_str("  top_offlist_sequences:\n");
                for entry in &region.top_offlist_sequences {
                    out.push_str(&format!("    {} {}\n", entry.sequence, entry.count));
                }
            }
        }
    }

    out
}
