use crate::context::load_resolved_inputs;
use crate::report::{
    format_fraction, render_report_prelude, write_report, FileReport, ReportEnvelope,
};
use crate::scan::scan_fastq;
use crate::CommonMetricArgs;
use anyhow::Result;
use serde::Serialize;

#[derive(Debug, clap::Args)]
pub struct LengthArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Serialize)]
pub struct LengthResult {
    pub sampled_count: usize,
    pub expected_min_len: i64,
    pub expected_max_len: i64,
    pub observed_min_len: usize,
    pub observed_max_len: usize,
    pub out_of_range_count: usize,
    pub out_of_range_fraction: f64,
}

#[derive(Debug, Default)]
struct LengthState {
    min_len: Option<usize>,
    max_len: usize,
    out_of_range_count: usize,
}

pub fn run(args: &LengthArgs) -> Result<()> {
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
        let expected_min_len = input.read.min_len;
        let expected_max_len = input.read.max_len;
        let (state, sampled_count) = scan_fastq(
            &input,
            args.common.n_reads,
            LengthState::default(),
            |state, _, _, len| {
                state.min_len = Some(state.min_len.map_or(len, |current| current.min(len)));
                state.max_len = state.max_len.max(len);
                if (len as i64) < expected_min_len || (len as i64) > expected_max_len {
                    state.out_of_range_count += 1;
                }
                Ok(())
            },
        )?;

        files.push(FileReport {
            input_path,
            read_id,
            file_id,
            matched_by,
            results: LengthResult {
                sampled_count,
                expected_min_len,
                expected_max_len,
                observed_min_len: state.min_len.unwrap_or_default(),
                observed_max_len: state.max_len,
                out_of_range_count: state.out_of_range_count,
                out_of_range_fraction: crate::report::fraction(
                    state.out_of_range_count,
                    sampled_count,
                ),
            },
        });
    }

    let report = ReportEnvelope {
        spec: args.common.spec.clone(),
        modality: args.common.modality.clone(),
        command: "length".to_string(),
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

fn render_text(report: &ReportEnvelope<LengthResult>) -> String {
    let mut out = render_report_prelude("length", report);

    for file in &report.files {
        let results = &file.results;
        out.push_str(&format!(
            "\nfile: {}\nread_id: {}\nfile_id: {}\nmatched_by: {}\nsampled_reads: {}\nexpected_length_range: {}-{}\nobserved_length_range: {}-{}\nout_of_range_reads: {} ({})\n",
            file.input_path.display(),
            file.read_id,
            file.file_id,
            file.matched_by,
            results.sampled_count,
            results.expected_min_len,
            results.expected_max_len,
            results.observed_min_len,
            results.observed_max_len,
            results.out_of_range_count,
            format_fraction(results.out_of_range_count, results.sampled_count),
        ));
    }

    out
}
