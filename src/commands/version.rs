use crate::context;
use crate::report::OutputFormat;
use anyhow::{bail, Result};
use clap::Args;
use serde::Serialize;
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug, Args)]
pub struct VersionArgs {
    #[arg(short = 'o', long, help = "Path to output file", value_name = "OUT")]
    pub output: Option<PathBuf>,

    #[arg(
        short = 's',
        long = "spec",
        visible_alias = "yaml",
        help = "Path to seqspec YAML file",
        value_name = "SPEC",
        required = true
    )]
    pub spec: PathBuf,

    #[arg(long, help = "Output format", value_enum, default_value_t = OutputFormat::Text)]
    pub format: OutputFormat,
}

#[derive(Debug, Clone, Serialize)]
struct VersionReport {
    seqcheck_version: String,
    seqspec_file_version: Option<String>,
    assay_id: String,
}

pub fn run(args: &VersionArgs) -> Result<()> {
    if !args.spec.exists() {
        bail!("spec file does not exist: {}", args.spec.display());
    }

    let spec = context::load_spec(&args.spec)?;
    let report = VersionReport {
        seqcheck_version: env!("CARGO_PKG_VERSION").to_string(),
        seqspec_file_version: spec.seqspec_version.clone(),
        assay_id: spec.assay_id,
    };

    let payload = match args.format {
        OutputFormat::Text => format!(
            "seqcheck version: {}\nseqspec file version: {}\nassay_id: {}",
            report.seqcheck_version,
            report
                .seqspec_file_version
                .clone()
                .unwrap_or_else(|| "unknown".to_string()),
            report.assay_id,
        ),
        OutputFormat::Json => serde_json::to_string_pretty(&report)?,
    };

    match &args.output {
        Some(path) => {
            let mut file = std::fs::File::create(path)?;
            file.write_all(payload.as_bytes())?;
            if !payload.ends_with('\n') {
                file.write_all(b"\n")?;
            }
        }
        None => {
            let mut stdout = std::io::stdout();
            stdout.write_all(payload.as_bytes())?;
            if !payload.ends_with('\n') {
                stdout.write_all(b"\n")?;
            }
        }
    }

    Ok(())
}
