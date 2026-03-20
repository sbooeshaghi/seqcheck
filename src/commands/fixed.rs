use crate::context::{filter_regions_by_sequence_type, load_resolved_inputs, LoadedInputs};
use crate::report::{
    input_check_result, top_sequences, write_report, AssessmentType, AtomicResult, Report,
    ResultBuilder,
};
use crate::scan::{extract_region, run_collector, FastqCollector};
use crate::sequence::{
    complement_sequence, find_all_exact_hits, reverse_complement_sequence, reverse_sequence,
};
use crate::CommonMetricArgs;
use anyhow::Result;
use std::collections::{BTreeMap, HashMap};

const TOP_SEQUENCE_LIMIT: usize = 10;

#[derive(Debug, clap::Args)]
pub struct FixedArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Default)]
pub struct OrientationCounts {
    pub forward: usize,
    pub reverse: usize,
    pub complement: usize,
    pub reverse_complement: usize,
}

#[derive(Debug, Clone, Default)]
struct FixedRegionState {
    covered_count: usize,
    short_read_count: usize,
    exact_match_count: usize,
    orientation_counts: OrientationCounts,
    offset_histogram: BTreeMap<i64, usize>,
    absent_count: usize,
    multi_hit_count: usize,
    mismatches: HashMap<String, usize>,
}

pub(crate) struct FixedCollector {
    file_id: String,
    read_id: String,
    primary_orientation: ExpectedOrientation,
    fixed_regions: Vec<seqspec::region::RegionCoordinate>,
    expected_sequences: Vec<ExpectedSequences>,
    state: Vec<FixedRegionState>,
}

pub fn run(args: &FixedArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
        args.common.auth_profile.as_deref(),
    )?;
    let mut report = Report::new(
        "fixed",
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
        results.extend(run_collector(input, n_reads, FixedCollector::new(input))?);
    }

    Ok(results)
}

impl FixedCollector {
    pub(crate) fn new(input: &crate::context::ResolvedInput) -> Self {
        let primary_orientation = if input.read.strand == "neg" {
            ExpectedOrientation::ReverseComplement
        } else {
            ExpectedOrientation::Forward
        };
        let fixed_regions = filter_regions_by_sequence_type(input, "fixed");
        let expected_sequences = fixed_regions
            .iter()
            .map(|region| expected_sequences(&region.region.sequence))
            .collect::<Vec<_>>();
        let state_len = expected_sequences.len();

        Self {
            file_id: input.file_id(),
            read_id: input.read.read_id.clone(),
            primary_orientation,
            fixed_regions,
            expected_sequences,
            state: vec![FixedRegionState::default(); state_len],
        }
    }
}

impl FastqCollector for FixedCollector {
    type Output = Vec<AtomicResult>;

    fn observe(&mut self, sequence: &str, _len: usize) -> Result<()> {
        for (idx, region) in self.fixed_regions.iter().enumerate() {
            let expected = &self.expected_sequences[idx];
            let primary_expected = self.primary_orientation.sequence(expected);
            let hits = find_all_exact_hits(sequence, primary_expected);
            if hits.is_empty() {
                self.state[idx].absent_count += 1;
            } else {
                if hits.len() > 1 {
                    self.state[idx].multi_hit_count += 1;
                }
                for hit in hits {
                    let offset = i64::try_from(hit).unwrap_or_default() - region.start;
                    *self.state[idx].offset_histogram.entry(offset).or_insert(0) += 1;
                }
            }

            match extract_region(sequence, region) {
                Some(observed) => {
                    self.state[idx].covered_count += 1;
                    update_orientation_counts(
                        &mut self.state[idx].orientation_counts,
                        &observed,
                        expected,
                    );

                    if observed == primary_expected {
                        self.state[idx].exact_match_count += 1;
                    } else {
                        *self.state[idx].mismatches.entry(observed).or_insert(0) += 1;
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
        for (idx, region) in self.fixed_regions.iter().enumerate() {
            results.push(build_fixed_region_result(
                &self.file_id,
                &self.read_id,
                sampled_count,
                region,
                &self.expected_sequences[idx],
                self.primary_orientation,
                &self.state[idx],
            ));
        }
        results
    }
}

fn build_fixed_region_result(
    file_id: &str,
    read_id: &str,
    sampled_count: usize,
    region: &seqspec::region::RegionCoordinate,
    _expected: &ExpectedSequences,
    primary_orientation: ExpectedOrientation,
    state: &FixedRegionState,
) -> crate::report::AtomicResult {
    let mut result = ResultBuilder::new(
        "fixed",
        vec![file_id.to_string()],
        vec![read_id.to_string()],
        vec![region.region.region_id.clone()],
    );

    let expected_sequence_id = result.expected_scalar(
        "expected_sequence",
        "Expected fixed sequence for this region from the seqspec.",
        region.region.sequence.clone(),
        Some("bases"),
    );
    let expected_coordinates_id = result.expected_records(
        "expected_coordinates",
        "Expected projected start and stop coordinates for this fixed region.",
        vec![serde_json::json!({
            "start": usize::try_from(region.start).unwrap_or_default(),
            "stop": usize::try_from(region.stop).unwrap_or_default(),
            "region_type": region.region.region_type,
            "name": region.region.name
        })],
    );
    let primary_orientation_label = primary_orientation.label();
    let primary_orientation_id = result.expected_scalar(
        "primary_orientation",
        "Primary biological orientation used for exact matching on this read.",
        primary_orientation_label,
        None,
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
    let exact_count_id = result.observed_scalar(
        "exact_match_count",
        "Number of covered reads whose extracted sequence exactly matches the expected motif in the primary orientation.",
        state.exact_match_count,
        Some("count"),
    );
    let exact_fraction_id = result.observed_scalar(
        "exact_match_fraction",
        "Fraction of covered reads whose extracted sequence exactly matches the expected motif in the primary orientation.",
        crate::report::fraction(state.exact_match_count, state.covered_count),
        Some("fraction"),
    );
    let orientation_id = result.observed_records(
        "orientation_counts",
        "Counts of exact motif matches under each orientation transform.",
        vec![
            serde_json::json!({ "orientation": "forward", "count": state.orientation_counts.forward }),
            serde_json::json!({ "orientation": "reverse", "count": state.orientation_counts.reverse }),
            serde_json::json!({ "orientation": "complement", "count": state.orientation_counts.complement }),
            serde_json::json!({ "orientation": "reverse_complement", "count": state.orientation_counts.reverse_complement }),
        ],
    );
    let offset_id = result.observed_series(
        "offset_histogram",
        "Distribution of exact-hit offsets relative to the expected region start.",
        state
            .offset_histogram
            .iter()
            .map(|(offset, count)| (format_offset(*offset), *count)),
        Some("count"),
    );
    let absent_id = result.observed_scalar(
        "absent_count",
        "Number of sampled reads with no exact motif hit anywhere in the read.",
        state.absent_count,
        Some("count"),
    );
    let multi_hit_id = result.observed_scalar(
        "multi_hit_count",
        "Number of sampled reads with more than one exact motif hit in the read.",
        state.multi_hit_count,
        Some("count"),
    );
    let top_mismatches_id = result.observed_records(
        "top_nonmatching_sequences",
        "Most frequent nonmatching extracted sequences observed at the expected coordinates.",
        top_sequences(&state.mismatches, TOP_SEQUENCE_LIMIT),
    );

    let (assessment_type, assessment_code, description) = if state.exact_match_count
        == state.covered_count
    {
        (
            AssessmentType::Pass,
            "fixed_exact_match_complete",
            format!(
                "All covered reads exactly match fixed region '{}' in {} orientation.",
                region.region.region_id, primary_orientation_label
            ),
        )
    } else if state.exact_match_count > 0 {
        (
            AssessmentType::Interpretation,
            "fixed_exact_match_partial",
            format!(
                "Fixed region '{}' matches exactly in some covered reads, with the dominant orientation '{}'.",
                region.region.region_id, primary_orientation_label
            ),
        )
    } else {
        (
            AssessmentType::Warning,
            "fixed_exact_match_absent",
            format!(
                "Fixed region '{}' does not exactly match at the expected coordinates in the primary orientation.",
                region.region.region_id
            ),
        )
    };
    result.assessment(
        assessment_type,
        assessment_code,
        description,
        vec![
            expected_sequence_id,
            expected_coordinates_id,
            primary_orientation_id,
        ],
        vec![
            sampled_id,
            covered_id,
            covered_fraction_id,
            short_id,
            exact_count_id,
            exact_fraction_id,
            orientation_id,
            offset_id,
            absent_id,
            multi_hit_id,
            top_mismatches_id,
        ],
    );

    result.build()
}

#[derive(Debug, Clone)]
struct ExpectedSequences {
    forward: String,
    reverse: String,
    complement: String,
    reverse_complement: String,
}

#[derive(Debug, Clone, Copy)]
enum ExpectedOrientation {
    Forward,
    ReverseComplement,
}

impl ExpectedOrientation {
    fn sequence<'a>(&self, expected: &'a ExpectedSequences) -> &'a str {
        match self {
            ExpectedOrientation::Forward => &expected.forward,
            ExpectedOrientation::ReverseComplement => &expected.reverse_complement,
        }
    }

    fn label(&self) -> &'static str {
        match self {
            ExpectedOrientation::Forward => "forward",
            ExpectedOrientation::ReverseComplement => "reverse_complement",
        }
    }
}

fn expected_sequences(sequence: &str) -> ExpectedSequences {
    ExpectedSequences {
        forward: sequence.to_string(),
        reverse: reverse_sequence(sequence),
        complement: complement_sequence(sequence),
        reverse_complement: reverse_complement_sequence(sequence),
    }
}

fn update_orientation_counts(
    counts: &mut OrientationCounts,
    observed: &str,
    expected: &ExpectedSequences,
) {
    if observed == expected.forward {
        counts.forward += 1;
    }
    if observed == expected.reverse {
        counts.reverse += 1;
    }
    if observed == expected.complement {
        counts.complement += 1;
    }
    if observed == expected.reverse_complement {
        counts.reverse_complement += 1;
    }
}

fn format_offset(offset: i64) -> String {
    if offset >= 0 {
        format!("+{}", offset)
    } else {
        offset.to_string()
    }
}
