use crate::context::{
    filter_regions_by_sequence_type, load_onlist, load_resolved_inputs, LoadedInputs, LoadedOnlist,
};
use crate::report::{
    input_check_result, top_sequences, write_report, AssessmentType, AtomicResult, Report,
    ResultBuilder,
};
use crate::scan::{extract_region, run_collector, FastqCollector};
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

fn onlist_error_source(spec_base: &Path, onlist: &seqspec::onlist::Onlist) -> String {
    if onlist.urltype == "local" {
        match seqspec::utils::local_onlist_locator(onlist) {
            Ok(locator) => spec_base.join(Path::new(locator)).display().to_string(),
            Err(err) => err,
        }
    } else {
        onlist.url.clone()
    }
}

pub(crate) struct OnlistCollector {
    file_id: String,
    read_id: String,
    regions: Vec<seqspec::region::RegionCoordinate>,
    loads: Vec<OnlistRegionLoad>,
    state: Vec<OnlistRegionState>,
}

pub fn run(args: &OnlistArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
        args.common.auth_profile.as_deref(),
    )?;
    let mut report = Report::new(
        "onlist",
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
        results.extend(run_collector(
            input,
            n_reads,
            OnlistCollector::new(input, &loaded.remote_access)?,
        )?);
    }

    Ok(results)
}

impl OnlistCollector {
    pub(crate) fn new(
        input: &crate::context::ResolvedInput,
        remote_access: &crate::auth::RemoteAccess,
    ) -> Result<Self> {
        let regions = filter_regions_by_sequence_type(input, "onlist");
        let loads = regions
            .iter()
            .map(
                |region| match load_onlist(&input.spec_base, &region.region, remote_access) {
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
                            .map(|onlist| onlist_error_source(&input.spec_base, onlist))
                            .unwrap_or_default(),
                        loaded_onlist: None,
                        error: Some(error.to_string()),
                    },
                },
            )
            .collect::<Vec<_>>();

        Ok(Self {
            file_id: input.file_id(),
            read_id: input.read.read_id.clone(),
            state: vec![OnlistRegionState::default(); regions.len()],
            regions,
            loads,
        })
    }
}

impl FastqCollector for OnlistCollector {
    type Output = Vec<AtomicResult>;

    fn observe(&mut self, sequence: &str, _len: usize) -> Result<()> {
        for (idx, region) in self.regions.iter().enumerate() {
            match extract_region(sequence, region) {
                Some(observed) => {
                    self.state[idx].covered_count += 1;
                    if let Some(loaded_onlist) = self.loads[idx].loaded_onlist.as_ref() {
                        if loaded_onlist.entries.contains(&observed) {
                            self.state[idx].exact_onlist_count += 1;
                        } else {
                            self.state[idx].offlist_count += 1;
                            *self.state[idx]
                                .offlist_sequences
                                .entry(observed)
                                .or_insert(0) += 1;
                        }
                    }
                }
                None => {
                    self.state[idx].short_read_count += 1;
                }
            }
        }
        Ok(())
    }

    fn finish(self, sampled_count: usize) -> Vec<AtomicResult> {
        let mut results = Vec::new();
        for (idx, region) in self.regions.iter().enumerate() {
            results.push(build_onlist_result(
                &self.file_id,
                &self.read_id,
                sampled_count,
                region,
                &self.loads[idx],
                &self.state[idx],
            ));
        }
        results
    }
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

#[cfg(test)]
mod tests {
    use super::onlist_error_source;
    use seqspec::onlist::Onlist;
    use std::path::Path;

    #[test]
    fn test_onlist_error_source_prefers_local_url() {
        let onlist = Onlist::new(
            "ol".to_string(),
            "display.txt".to_string(),
            "txt".to_string(),
            0,
            "nested/whitelist.txt".to_string(),
            "local".to_string(),
            String::new(),
        );

        let source = onlist_error_source(Path::new("/tmp/spec-root"), &onlist);
        assert_eq!(source, "/tmp/spec-root/nested/whitelist.txt");
    }

    #[test]
    fn test_onlist_error_source_reports_empty_local_url() {
        let onlist = Onlist::new(
            "ol".to_string(),
            "display.txt".to_string(),
            "txt".to_string(),
            0,
            String::new(),
            "local".to_string(),
            String::new(),
        );

        let source = onlist_error_source(Path::new("/tmp/spec-root"), &onlist);
        assert_eq!(source, "local onlist 'display.txt' has empty url");
    }
}
