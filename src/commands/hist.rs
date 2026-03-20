use crate::context::{filter_regions_by_id, load_resolved_inputs};
use crate::report::{
    input_check_result, top_sequences, write_report, AssessmentType, Report, ResultBuilder,
};
use crate::scan::{extract_region, scan_fastq};
use crate::CommonMetricArgs;
use anyhow::Result;
use std::collections::HashMap;

#[derive(Debug, clap::Args)]
pub struct HistArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,

    #[arg(short = 'r', long, help = "Region ID to summarize", required = true)]
    pub region_id: String,
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
        args.common.auth_profile.as_deref(),
    )?;
    let crate::context::LoadedInputs {
        input_check,
        inputs,
        ..
    } = loaded;
    let mut report = Report::new(
        "hist",
        args.common.spec.clone(),
        &args.common.modality,
        args.common.n_reads,
    );
    report
        .results
        .push(input_check_result(&input_check, &inputs));

    for input in inputs {
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

        let mut result = ResultBuilder::new(
            "hist",
            vec![input.file_id()],
            vec![input.read.read_id.clone()],
            vec![region.region.region_id.clone()],
        );
        let expected_region_id = result.expected_records(
            "expected_region",
            "Expected projected coordinates and annotations for the summarized region.",
            vec![serde_json::json!({
                "start": usize::try_from(region.start).unwrap_or_default(),
                "stop": usize::try_from(region.stop).unwrap_or_default(),
                "region_type": region.region.region_type,
                "name": region.region.name
            })],
        );
        let sampled_id = result.observed_scalar(
            "sampled_count",
            "Number of reads sampled from the FASTQ file.",
            sampled_count,
            Some("count"),
        );
        let covered_id = result.observed_scalar(
            "covered_count",
            "Number of sampled reads that fully cover the summarized region.",
            state.covered_count,
            Some("count"),
        );
        let short_id = result.observed_scalar(
            "short_read_count",
            "Number of sampled reads that do not extend far enough to cover the summarized region.",
            state.short_read_count,
            Some("count"),
        );
        let histogram_id = result.observed_records(
            "sequence_histogram",
            "Exact histogram of extracted sequences for this region.",
            top_sequences(&state.counts, usize::MAX),
        );
        result.assessment(
            AssessmentType::Interpretation,
            "region_histogram_computed",
            format!(
                "Computed the exact sequence histogram for region '{}'.",
                region.region.region_id
            ),
            vec![expected_region_id],
            vec![sampled_id, covered_id, short_id, histogram_id],
        );

        report.results.push(result.build());
    }

    write_report(&args.common.output, args.common.format, &report)
}
