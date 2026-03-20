pub mod commands;
pub mod context;
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
    Coverage(commands::coverage::CoverageArgs),
    Cut(commands::cut::CutArgs),
    Fixed(commands::fixed::FixedArgs),
    Hist(commands::hist::HistArgs),
    Length(commands::length::LengthArgs),
    Onlist(commands::onlist::OnlistArgs),
    Primer(commands::primer::PrimerArgs),
    Random(commands::random::RandomArgs),
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
        help = "Path to seqspec YAML file",
        value_name = "SPEC",
        required = true
    )]
    pub spec: PathBuf,

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

    #[arg(help = "FASTQ files to inspect", required = true, value_name = "FASTQ")]
    pub fastqs: Vec<PathBuf>,
}

pub fn run() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Coverage(args) => commands::coverage::run(&args),
        Commands::Cut(args) => commands::cut::run(&args),
        Commands::Fixed(args) => commands::fixed::run(&args),
        Commands::Hist(args) => commands::hist::run(&args),
        Commands::Length(args) => commands::length::run(&args),
        Commands::Onlist(args) => commands::onlist::run(&args),
        Commands::Primer(args) => commands::primer::run(&args),
        Commands::Random(args) => commands::random::run(&args),
        Commands::Version(args) => commands::version::run(&args),
    }
}
