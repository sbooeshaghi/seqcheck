use crate::context::load_resolved_inputs;
use crate::report::{input_check_result, write_report, AssessmentType, Report, ResultBuilder};
use crate::scan::scan_fastq;
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

pub fn run(args: &LengthArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
        args.common.auth_profile.as_deref(),
    )?;
    let crate::context::LoadedInputs {
        input_check,
        inputs,
        ..
    } = loaded;
    let mut report = Report::new(
        "length",
        args.common.spec.clone(),
        &args.common.modality,
        args.common.n_reads,
    );
    report
        .results
        .push(input_check_result(&input_check, &inputs));

    for input in inputs {
        let read_id = input.read.read_id.clone();
        let file_id = input.file_id();
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

        let observed_min_len = state.min_len.unwrap_or_default();
        let observed_max_len = state.max_len;
        let out_of_range_fraction =
            crate::report::fraction(state.out_of_range_count, sampled_count);

        let mut result = ResultBuilder::new("length", vec![file_id], vec![read_id], Vec::new());
        let expected_min_id = result.expected_scalar(
            "expected_min_len",
            "Minimum read length allowed by the seqspec for this read.",
            expected_min_len,
            Some("bp"),
        );
        let expected_max_id = result.expected_scalar(
            "expected_max_len",
            "Maximum read length allowed by the seqspec for this read.",
            expected_max_len,
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
            state.out_of_range_count,
            Some("count"),
        );
        let out_of_range_fraction_id = result.observed_scalar(
            "out_of_range_fraction",
            "Fraction of sampled reads whose length falls outside the seqspec range.",
            out_of_range_fraction,
            Some("fraction"),
        );

        if state.out_of_range_count == 0 {
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

        report.results.push(result.build());
    }

    write_report(&args.common.output, args.common.format, &report)
}
