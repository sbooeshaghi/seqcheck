use crate::context::{filter_regions_by_id, load_resolved_inputs};
use crate::report::{
    render_report_prelude, top_sequences, write_report, FileReport, ReportEnvelope, SequenceCount,
};
use crate::scan::{extract_region, scan_fastq};
use crate::CommonMetricArgs;
use anyhow::Result;
use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, clap::Args)]
pub struct HistArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,

    #[arg(short = 'r', long, help = "Region ID to summarize", required = true)]
    pub region_id: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct HistResult {
    pub region_id: String,
    pub region_type: String,
    pub sampled_count: usize,
    pub covered_count: usize,
    pub short_read_count: usize,
    pub histogram: Vec<SequenceCount>,
}

#[derive(Debug, Default)]
struct HistState {
    covered_count: usize,
    short_read_count: usize,
    counts: HashMap<String, usize>,
}

pub fn run(args: &HistArgs) -> Result<()> {
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
        let regions = filter_regions_by_id(&input, &args.region_id)?;
        let region = regions.first().unwrap().clone();
        let target = region.clone();

        let (state, sampled_count) = scan_fastq(
            &input,
            args.common.n_reads,
            HistState::default(),
            |state, _, sequence, _| {
                match extract_region(sequence, &target) {
                    Some(observed) => {
                        state.covered_count += 1;
                        *state.counts.entry(observed).or_insert(0) += 1;
                    }
                    None => {
                        state.short_read_count += 1;
                    }
                }
                Ok(())
            },
        )?;

        files.push(FileReport {
            input_path,
            read_id,
            file_id,
            matched_by,
            results: HistResult {
                region_id: region.region.region_id.clone(),
                region_type: region.region.region_type.clone(),
                sampled_count,
                covered_count: state.covered_count,
                short_read_count: state.short_read_count,
                histogram: top_sequences(&state.counts, usize::MAX),
            },
        });
    }

    let report = ReportEnvelope {
        spec: args.common.spec.clone(),
        modality: args.common.modality.clone(),
        command: "hist".to_string(),
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

fn render_text(report: &ReportEnvelope<HistResult>) -> String {
    let mut out = render_report_prelude("hist", report);

    for file in &report.files {
        out.push_str(&format!(
            "\nfile: {}\nread_id: {}\nfile_id: {}\nmatched_by: {}\nregion: {}\nsampled_reads: {}\ncovered_reads: {}\nshort_reads: {}\n",
            file.input_path.display(),
            file.read_id,
            file.file_id,
            file.matched_by,
            file.results.region_id,
            file.results.sampled_count,
            file.results.covered_count,
            file.results.short_read_count,
        ));
        for entry in &file.results.histogram {
            out.push_str(&format!("  {} {}\n", entry.sequence, entry.count));
        }
    }

    out
}
