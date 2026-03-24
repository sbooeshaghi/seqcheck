use crate::html_report::{
    generated_timestamp_utc, load_optional_seqspec_lib_data, load_report_json, render_report_html,
};
use anyhow::{Context, Result};
use clap::Args;
use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone, Args)]
pub struct ReportArgs {
    #[arg(
        short = 'i',
        long = "input",
        help = "Path to a seqcheck JSON report",
        value_name = "REPORT",
        required = true
    )]
    pub input: PathBuf,

    #[arg(
        short = 'o',
        long = "output",
        help = "Path to output HTML report",
        value_name = "HTML",
        required = true
    )]
    pub output: PathBuf,

    #[arg(
        long = "spec",
        visible_alias = "yaml",
        help = "Optional seqspec YAML path to use for the library diagram",
        value_name = "SPEC"
    )]
    pub spec: Option<PathBuf>,
}

pub fn run(args: &ReportArgs) -> Result<()> {
    let report = load_report_json(&args.input)?;
    let lib_data = load_optional_seqspec_lib_data(&args.input, &report, args.spec.as_deref());
    let generated_at = generated_timestamp_utc();
    let html = render_report_html(
        &report,
        lib_data.as_ref(),
        &generated_at,
        &args.input,
        &args.output,
    )?;

    fs::write(&args.output, html)
        .with_context(|| format!("failed to write {}", args.output.display()))?;

    Ok(())
}
