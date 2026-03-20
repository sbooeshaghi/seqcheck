use crate::context::{
    filter_regions_by_sequence_type, load_onlist, load_resolved_inputs, LoadedOnlist,
};
use crate::report::{
    input_check_result, top_sequences, write_report, AssessmentType, Report, ResultBuilder,
};
use crate::scan::{extract_region, scan_fastq};
use crate::CommonMetricArgs;
use anyhow::Result;
use std::collections::HashMap;
use std::path::Path;

const TOP_SEQUENCE_LIMIT: usize = 10;

#[derive(Debug, clap::Args)]
pub struct OnlistArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Default)]
struct OnlistRegionState {
    covered_count: usize,
    short_read_count: usize,
    exact_onlist_count: usize,
    offlist_count: usize,
    offlist_sequences: HashMap<String, usize>,
}

#[derive(Debug, Clone)]
struct OnlistRegionLoad {
    source: String,
    loaded_onlist: Option<LoadedOnlist>,
    error: Option<String>,
}

pub fn run(args: &OnlistArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
        args.common.auth_profile.as_deref(),
    )?;
    let crate::context::LoadedInputs {
        input_check,
        inputs,
        remote_access,
        ..
    } = loaded;
    let mut report = Report::new(
        "onlist",
        args.common.spec.clone(),
        &args.common.modality,
        args.common.n_reads,
    );
    report
        .results
        .push(input_check_result(&input_check, &inputs));

    for input in inputs {
        let onlist_regions = filter_regions_by_sequence_type(&input, "onlist");
        let loaded_onlists = onlist_regions
            .iter()
            .map(
                |region| match load_onlist(&input.spec_base, &region.region, &remote_access) {
                    Ok(loaded_onlist) => OnlistRegionLoad {
                        source: loaded_onlist.source.clone(),
                        loaded_onlist: Some(loaded_onlist),
                        error: None,
                    },
                    Err(error) => OnlistRegionLoad {
                        source: region
                            .region
                            .onlist
                            .as_ref()
                            .map(|onlist| {
                                if onlist.urltype == "local" {
                                    input
                                        .spec_base
                                        .join(Path::new(&onlist.filename))
                                        .display()
                                        .to_string()
                                } else {
                                    onlist.filename.clone()
                                }
                            })
                            .unwrap_or_default(),
                        loaded_onlist: None,
                        error: Some(error.to_string()),
                    },
                },
            )
            .collect::<Vec<_>>();
        let state = vec![OnlistRegionState::default(); onlist_regions.len()];
        let regions = onlist_regions.clone();
        let loads = loaded_onlists.clone();

        let (state, sampled_count) = scan_fastq(
            &input,
            args.common.n_reads,
            state,
            |state, _, sequence, _| {
                for (idx, region) in regions.iter().enumerate() {
                    match extract_region(sequence, region) {
                        Some(observed) => {
                            state[idx].covered_count += 1;
                            if let Some(loaded_onlist) = loads[idx].loaded_onlist.as_ref() {
                                if loaded_onlist.entries.contains(&observed) {
                                    state[idx].exact_onlist_count += 1;
                                } else {
                                    state[idx].offlist_count += 1;
                                    *state[idx].offlist_sequences.entry(observed).or_insert(0) += 1;
                                }
                            }
                        }
                        None => {
                            state[idx].short_read_count += 1;
                        }
                    }
                }
                Ok(())
            },
        )?;

        for (idx, region) in onlist_regions.iter().enumerate() {
            report.results.push(build_onlist_result(
                &input.file_id(),
                &input.read.read_id,
                sampled_count,
                region,
                &loaded_onlists[idx],
                &state[idx],
            ));
        }
    }

    write_report(&args.common.output, args.common.format, &report)
}

fn build_onlist_result(
    file_id: &str,
    read_id: &str,
    sampled_count: usize,
    region: &seqspec::region::RegionCoordinate,
    load: &OnlistRegionLoad,
    state: &OnlistRegionState,
) -> crate::report::AtomicResult {
    let mut result = ResultBuilder::new(
        "onlist",
        vec![file_id.to_string()],
        vec![read_id.to_string()],
        vec![region.region.region_id.clone()],
    );

    let source_id = result.expected_scalar(
        "onlist_source",
        "Remote or local source declared for this whitelist.",
        load.source.clone(),
        None,
    );
    let region_id = result.expected_records(
        "expected_region",
        "Expected projected coordinates and annotations for this onlist region.",
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
        "Number of sampled reads that fully cover the expected region coordinates.",
        state.covered_count,
        Some("count"),
    );
    let covered_fraction_id = result.observed_scalar(
        "covered_fraction",
        "Fraction of sampled reads that fully cover the expected region coordinates.",
        crate::report::fraction(state.covered_count, sampled_count),
        Some("fraction"),
    );
    let short_id = result.observed_scalar(
        "short_read_count",
        "Number of sampled reads that do not extend far enough to cover the expected region coordinates.",
        state.short_read_count,
        Some("count"),
    );
    let status_id = result.observed_scalar(
        "fetch_load_status",
        "Whether the whitelist source was loaded successfully.",
        if load.loaded_onlist.is_some() {
            "loaded"
        } else {
            "error"
        },
        None,
    );
    let onlist_size_id = result.observed_scalar(
        "onlist_entry_count",
        "Number of entries available in the loaded whitelist.",
        load.loaded_onlist
            .as_ref()
            .map(|loaded_onlist| loaded_onlist.entries.len()),
        Some("count"),
    );
    let onlist_count_id = result.observed_scalar(
        "exact_onlist_count",
        "Number of covered reads whose extracted sequence is present in the whitelist.",
        load.loaded_onlist
            .as_ref()
            .map(|_| state.exact_onlist_count),
        Some("count"),
    );
    let onlist_fraction_id = result.observed_scalar(
        "exact_onlist_fraction",
        "Fraction of covered reads whose extracted sequence is present in the whitelist.",
        load.loaded_onlist
            .as_ref()
            .map(|_| crate::report::fraction(state.exact_onlist_count, state.covered_count)),
        Some("fraction"),
    );
    let offlist_count_id = result.observed_scalar(
        "offlist_count",
        "Number of covered reads whose extracted sequence is absent from the whitelist.",
        load.loaded_onlist.as_ref().map(|_| state.offlist_count),
        Some("count"),
    );
    let offlist_sequences_id = result.observed_records(
        "top_offlist_sequences",
        "Most frequent offlist sequences observed for this region.",
        if load.loaded_onlist.is_some() {
            top_sequences(&state.offlist_sequences, TOP_SEQUENCE_LIMIT)
        } else {
            Vec::new()
        },
    );

    if let Some(error) = &load.error {
        result.assessment(
            AssessmentType::Error,
            "missing_onlist_resource",
            format!("The whitelist source could not be loaded: {}", error),
            vec![source_id, region_id],
            vec![
                sampled_id,
                covered_id,
                covered_fraction_id,
                short_id,
                status_id,
                onlist_size_id,
                onlist_count_id,
                onlist_fraction_id,
                offlist_count_id,
                offlist_sequences_id,
            ],
        );
    } else if state.offlist_count == 0 {
        result.assessment(
            AssessmentType::Pass,
            "all_sequences_onlist",
            "All covered sequences are present in the whitelist.",
            vec![source_id, region_id],
            vec![
                sampled_id,
                covered_id,
                covered_fraction_id,
                short_id,
                status_id,
                onlist_size_id,
                onlist_count_id,
                onlist_fraction_id,
                offlist_count_id,
                offlist_sequences_id,
            ],
        );
    } else {
        result.assessment(
            AssessmentType::Warning,
            "offlist_sequences_detected",
            "Some covered sequences are absent from the whitelist.",
            vec![source_id, region_id],
            vec![
                sampled_id,
                covered_id,
                covered_fraction_id,
                short_id,
                status_id,
                onlist_size_id,
                onlist_count_id,
                onlist_fraction_id,
                offlist_count_id,
                offlist_sequences_id,
            ],
        );
    }

    result.build()
}
