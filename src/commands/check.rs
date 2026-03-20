use crate::commands::coverage::{self, CoverageCollector};
use crate::commands::fixed::FixedCollector;
use crate::commands::length::LengthCollector;
use crate::commands::onlist::OnlistCollector;
use crate::commands::primer::PrimerCollector;
use crate::commands::random::RandomCollector;
use crate::context::load_resolved_inputs;
use crate::report::{input_check_result, write_report, AtomicResult, Report};
use crate::scan::{run_collector, FastqCollector};
use crate::CommonMetricArgs;
use anyhow::Result;

#[derive(Debug, clap::Args)]
pub struct CheckArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

struct CheckCollectorBundle {
    length: LengthCollector,
    coverage: CoverageCollector,
    primer: PrimerCollector,
    fixed: FixedCollector,
    onlist: OnlistCollector,
    random: RandomCollector,
}

struct CheckCollectorOutput {
    results: Vec<AtomicResult>,
    coverage_summary: coverage::CoverageFileSummary,
}

pub fn run(args: &CheckArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
        args.common.auth_profile.as_deref(),
    )?;

    let mut report = Report::new(
        "check",
        args.common.spec.clone(),
        &args.common.modality,
        args.common.n_reads,
    );
    report
        .results
        .push(input_check_result(&loaded.input_check, &loaded.inputs));

    let mut coverage_summaries = Vec::new();

    for input in &loaded.inputs {
        let output = run_collector(
            input,
            args.common.n_reads,
            CheckCollectorBundle::new(input, &loaded.remote_access)?,
        )?;
        report.results.extend(output.results);
        coverage_summaries.push(output.coverage_summary);
    }

    report
        .results
        .extend(coverage::build_overlap_results(&coverage_summaries));

    write_report(&args.common.output, args.common.format, &report)
}

impl CheckCollectorBundle {
    fn new(
        input: &crate::context::ResolvedInput,
        remote_access: &crate::auth::RemoteAccess,
    ) -> Result<Self> {
        Ok(Self {
            length: LengthCollector::new(input),
            coverage: CoverageCollector::new(input),
            primer: PrimerCollector::new(input),
            fixed: FixedCollector::new(input),
            onlist: OnlistCollector::new(input, remote_access)?,
            random: RandomCollector::new(input),
        })
    }
}

impl FastqCollector for CheckCollectorBundle {
    type Output = CheckCollectorOutput;

    fn observe(&mut self, sequence: &str, len: usize) -> Result<()> {
        self.length.observe(sequence, len)?;
        self.coverage.observe(sequence, len)?;
        self.primer.observe(sequence, len)?;
        self.fixed.observe(sequence, len)?;
        self.onlist.observe(sequence, len)?;
        self.random.observe(sequence, len)?;
        Ok(())
    }

    fn finish(self, sampled_count: usize) -> CheckCollectorOutput {
        let mut results = Vec::new();
        results.push(self.length.finish(sampled_count));

        let coverage_output = self.coverage.finish(sampled_count);
        results.extend(coverage_output.results);

        results.push(self.primer.finish(sampled_count));
        results.extend(self.fixed.finish(sampled_count));
        results.extend(self.onlist.finish(sampled_count));
        results.extend(self.random.finish(sampled_count));

        CheckCollectorOutput {
            results,
            coverage_summary: coverage_output.summary,
        }
    }
}
