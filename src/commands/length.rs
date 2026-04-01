use crate::context::{load_resolved_inputs, LoadedInputs};
use crate::report::{
    input_check_result, write_report, AssessmentType, AtomicResult, Report, ResultBuilder,
};
use crate::scan::{run_collector, FastqCollector};
use crate::CommonMetricArgs;
use anyhow::Result;

#[derive(Debug, clap::Args)]
pub struct LengthArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Default)]
struct LengthState {
    min_len: Option<usize>,
    max_len: usize,
    out_of_range_count: usize,
}

pub(crate) struct LengthCollector {
    file_id: String,
    read_id: String,
    expected_min_len: i64,
    expected_max_len: i64,
    state: LengthState,
}

pub fn run(args: &LengthArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
        args.common.n_reads,
        args.common.auth_profile.as_deref(),
    )?;
    let mut report = Report::new(
        "length",
        args.common.spec.clone(),
        &args.common.modality,
        args.common.n_reads,
    );
    report
        .results
        .push(input_check_result(&loaded.input_check, &loaded.inputs));
    report
        .results
        .extend(collect_results(&loaded, args.common.n_reads)?);

    write_report(&args.common.output, args.common.format, &report)
}

pub fn collect_results(loaded: &LoadedInputs, n_reads: usize) -> Result<Vec<AtomicResult>> {
    let mut results = Vec::new();

    for input in &loaded.inputs {
        results.push(run_collector(input, n_reads, LengthCollector::new(input))?);
    }

    Ok(results)
}

impl LengthCollector {
    pub(crate) fn new(input: &crate::context::ResolvedInput) -> Self {
        Self {
            file_id: input.file_id(),
            read_id: input.read.read_id.clone(),
            expected_min_len: input.read.min_len,
            expected_max_len: input.read.max_len,
            state: LengthState::default(),
        }
    }
}

impl FastqCollector for LengthCollector {
    type Output = AtomicResult;

    fn observe(&mut self, _sequence: &str, len: usize) -> Result<()> {
        self.state.min_len = Some(self.state.min_len.map_or(len, |current| current.min(len)));
        self.state.max_len = self.state.max_len.max(len);
        if (len as i64) < self.expected_min_len || (len as i64) > self.expected_max_len {
            self.state.out_of_range_count += 1;
        }
        Ok(())
    }

    fn finish(self, sampled_count: usize) -> AtomicResult {
        let observed_min_len = self.state.min_len.unwrap_or_default();
        let observed_max_len = self.state.max_len;
        let out_of_range_fraction =
            crate::report::fraction(self.state.out_of_range_count, sampled_count);

        let mut result =
            ResultBuilder::new("length", vec![self.file_id], vec![self.read_id], Vec::new());
        let expected_min_id = result.expected_scalar(
            "expected_min_len",
            "Minimum read length allowed by the seqspec for this read.",
            self.expected_min_len,
            Some("bp"),
        );
        let expected_max_id = result.expected_scalar(
            "expected_max_len",
            "Maximum read length allowed by the seqspec for this read.",
            self.expected_max_len,
            Some("bp"),
        );
        let sampled_id = result.observed_scalar(
            "sampled_count",
            "Number of reads sampled from the FASTQ file.",
            sampled_count,
            Some("count"),
        );
        let observed_min_id = result.observed_scalar(
            "observed_min_len",
            "Shortest observed read length in the sampled reads.",
            observed_min_len,
            Some("bp"),
        );
        let observed_max_id = result.observed_scalar(
            "observed_max_len",
            "Longest observed read length in the sampled reads.",
            observed_max_len,
            Some("bp"),
        );
        let out_of_range_count_id = result.observed_scalar(
            "out_of_range_count",
            "Number of sampled reads whose length falls outside the seqspec range.",
            self.state.out_of_range_count,
            Some("count"),
        );
        let out_of_range_fraction_id = result.observed_scalar(
            "out_of_range_fraction",
            "Fraction of sampled reads whose length falls outside the seqspec range.",
            out_of_range_fraction,
            Some("fraction"),
        );

        if self.state.out_of_range_count == 0 {
            result.assessment(
                AssessmentType::Pass,
                "length_in_range",
                "All sampled reads fall within the seqspec length range.",
                vec![expected_min_id, expected_max_id],
                vec![
                    sampled_id,
                    observed_min_id,
                    observed_max_id,
                    out_of_range_count_id,
                    out_of_range_fraction_id,
                ],
            );
        } else {
            result.assessment(
                AssessmentType::Warning,
                "length_out_of_range",
                "Observed read lengths fall outside the seqspec range.",
                vec![expected_min_id, expected_max_id],
                vec![
                    sampled_id,
                    observed_min_id,
                    observed_max_id,
                    out_of_range_count_id,
                    out_of_range_fraction_id,
                ],
            );
        }

        result.build()
    }
}
