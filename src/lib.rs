pub mod auth;
pub mod commands;
pub mod context;
pub mod html_report;
pub mod report;
pub mod scan;
pub mod sequence;

use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use report::OutputFormat;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "seqcheck",
    version,
    about = "Validate FASTQ reads against a seqspec file"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    #[command(hide = true)]
    Auth(commands::auth::AuthArgs),
    Check(commands::check::CheckArgs),
    Coverage(commands::coverage::CoverageArgs),
    Cut(commands::cut::CutArgs),
    Fixed(commands::fixed::FixedArgs),
    Hist(commands::hist::HistArgs),
    Length(commands::length::LengthArgs),
    Onlist(commands::onlist::OnlistArgs),
    Primer(commands::primer::PrimerArgs),
    Random(commands::random::RandomArgs),
    Report(commands::report::ReportArgs),
    Version(commands::version::VersionArgs),
}

#[derive(Debug, Clone, Args)]
pub struct CommonMetricArgs {
    #[arg(short = 'o', long, help = "Path to output file", value_name = "OUT")]
    pub output: Option<PathBuf>,

    #[arg(
        short = 'm',
        long,
        help = "Modality to inspect",
        value_name = "MODALITY",
        required = true
    )]
    pub modality: String,

    #[arg(
        short = 's',
        long = "spec",
        visible_alias = "yaml",
        help = "Path or URL to seqspec YAML file",
        value_name = "SPEC",
        required = true
    )]
    pub spec: String,

    #[arg(
        short = 'n',
        long = "n-reads",
        help = "Number of reads to inspect per FASTQ (0 means all reads)",
        default_value_t = 10000,
        value_name = "N"
    )]
    pub n_reads: usize,

    #[arg(
        long,
        help = "Output format",
        value_enum,
        default_value_t = OutputFormat::Text
    )]
    pub format: OutputFormat,

    #[arg(
        long,
        env = "SEQCHECK_AUTH_PROFILE",
        help = "Auth profile name for remote resources declared in the seqspec",
        value_name = "PROFILE"
    )]
    pub auth_profile: Option<String>,

    #[arg(
        help = "FASTQ files or URLs to inspect",
        required = true,
        value_name = "FASTQ"
    )]
    pub fastqs: Vec<String>,
}

pub fn run() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Auth(args) => commands::auth::run(&args),
        Commands::Check(args) => commands::check::run(&args),
        Commands::Coverage(args) => commands::coverage::run(&args),
        Commands::Cut(args) => commands::cut::run(&args),
        Commands::Fixed(args) => commands::fixed::run(&args),
        Commands::Hist(args) => commands::hist::run(&args),
        Commands::Length(args) => commands::length::run(&args),
        Commands::Onlist(args) => commands::onlist::run(&args),
        Commands::Primer(args) => commands::primer::run(&args),
        Commands::Random(args) => commands::random::run(&args),
        Commands::Report(args) => commands::report::run(&args),
        Commands::Version(args) => commands::version::run(&args),
    }
}
