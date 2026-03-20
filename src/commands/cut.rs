use crate::context::{filter_regions_by_id, load_resolved_inputs};
use crate::report::{render_report_prelude, write_report, FileReport, ReportEnvelope};
use crate::scan::{extract_region, scan_fastq};
use crate::CommonMetricArgs;
use anyhow::Result;
use serde::Serialize;

#[derive(Debug, clap::Args)]
pub struct CutArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,

    #[arg(short = 'r', long, help = "Region ID to extract", required = true)]
    pub region_id: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct CutSample {
    pub read_index: usize,
    pub sequence: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct CutResult {
    pub region_id: String,
    pub region_type: String,
    pub sampled_count: usize,
    pub covered_count: usize,
    pub short_read_count: usize,
    pub sequences: Vec<CutSample>,
}

#[derive(Debug, Default)]
struct CutState {
    covered_count: usize,
    short_read_count: usize,
    sequences: Vec<CutSample>,
}

pub fn run(args: &CutArgs) -> Result<()> {
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
            CutState::default(),
            |state, read_index, sequence, _| {
                match extract_region(sequence, &target) {
                    Some(observed) => {
                        state.covered_count += 1;
                        state.sequences.push(CutSample {
                            read_index,
                            sequence: observed,
                        });
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
            results: CutResult {
                region_id: region.region.region_id.clone(),
                region_type: region.region.region_type.clone(),
                sampled_count,
                covered_count: state.covered_count,
                short_read_count: state.short_read_count,
                sequences: state.sequences,
            },
        });
    }

    let report = ReportEnvelope {
        spec: args.common.spec.clone(),
        modality: args.common.modality.clone(),
        command: "cut".to_string(),
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

fn render_text(report: &ReportEnvelope<CutResult>) -> String {
    let mut out = render_report_prelude("cut", report);

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
        for sample in &file.results.sequences {
            out.push_str(&format!("  {} {}\n", sample.read_index, sample.sequence));
        }
    }

    out
}
