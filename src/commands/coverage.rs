use crate::context::{load_resolved_inputs, ExpectedRegion, LoadedInputs};
use crate::report::{
    input_check_result, write_report, AssessmentType, AtomicResult, Report, ResultBuilder,
};
use crate::scan::{is_region_covered, run_collector, FastqCollector};
use crate::CommonMetricArgs;
use anyhow::Result;
use std::collections::BTreeMap;

const SHARED_REGION_TYPES: &[&str] = &[
    "barcode",
    "umi",
    "cdna",
    "gdna",
    "protein",
    "tag",
    "sgrna_target",
];

#[derive(Debug, clap::Args)]
pub struct CoverageArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone)]
pub(crate) struct CoverageRegionSummary {
    region_id: String,
    name: String,
    region_type: String,
    sequence_type: String,
    start: usize,
    stop: usize,
    covered_count: usize,
    covered_fraction: f64,
}

#[derive(Debug, Clone)]
pub(crate) struct CoverageFileSummary {
    file_id: String,
    read_id: String,
    expected_regions: Vec<ExpectedRegion>,
    sampled_count: usize,
    full_read_coverage_count: usize,
    full_read_coverage_fraction: f64,
    regions: Vec<CoverageRegionSummary>,
}

#[derive(Debug)]
struct CoverageState {
    full_read_coverage_count: usize,
    region_counts: Vec<usize>,
}

pub(crate) struct CoverageCollector {
    file_id: String,
    read_id: String,
    expected_regions: Vec<ExpectedRegion>,
    expected_stop: usize,
    coordinates: Vec<seqspec::region::RegionCoordinate>,
    state: CoverageState,
}

pub(crate) struct CoverageCollectorOutput {
    pub(crate) results: Vec<AtomicResult>,
    pub(crate) summary: CoverageFileSummary,
}

pub fn run(args: &CoverageArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
        args.common.auth_profile.as_deref(),
    )?;
    let mut report = Report::new(
        "coverage",
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
    let mut summaries = Vec::new();

    for input in &loaded.inputs {
        let output = run_collector(input, n_reads, CoverageCollector::new(input))?;
        results.extend(output.results);
        summaries.push(output.summary);
    }

    results.extend(build_overlap_results(&summaries));
    Ok(results)
}

impl CoverageCollector {
    pub(crate) fn new(input: &crate::context::ResolvedInput) -> Self {
        Self {
            file_id: input.file_id(),
            read_id: input.read.read_id.clone(),
            expected_regions: input.expected_regions(),
            expected_stop: input.expected_stop(),
            coordinates: input.coordinates.clone(),
            state: CoverageState {
                full_read_coverage_count: 0,
                region_counts: vec![0; input.coordinates.len()],
            },
        }
    }
}

impl FastqCollector for CoverageCollector {
    type Output = CoverageCollectorOutput;

    fn observe(&mut self, _sequence: &str, len: usize) -> Result<()> {
        if len >= self.expected_stop {
            self.state.full_read_coverage_count += 1;
        }

        for (idx, region) in self.coordinates.iter().enumerate() {
            if is_region_covered(len, region) {
                self.state.region_counts[idx] += 1;
            }
        }

        Ok(())
    }

    fn finish(self, sampled_count: usize) -> CoverageCollectorOutput {
        let regions = self
            .coordinates
            .iter()
            .enumerate()
            .map(|(idx, region)| CoverageRegionSummary {
                region_id: region.region.region_id.clone(),
                name: region.region.name.clone(),
                region_type: region.region.region_type.clone(),
                sequence_type: region.region.sequence_type.clone(),
                start: usize::try_from(region.start).unwrap_or_default(),
                stop: usize::try_from(region.stop).unwrap_or_default(),
                covered_count: self.state.region_counts[idx],
                covered_fraction: crate::report::fraction(
                    self.state.region_counts[idx],
                    sampled_count,
                ),
            })
            .collect::<Vec<_>>();

        let summary = CoverageFileSummary {
            file_id: self.file_id,
            read_id: self.read_id,
            expected_regions: self.expected_regions,
            sampled_count,
            full_read_coverage_count: self.state.full_read_coverage_count,
            full_read_coverage_fraction: crate::report::fraction(
                self.state.full_read_coverage_count,
                sampled_count,
            ),
            regions,
        };

        let mut results = vec![build_file_coverage_result(&summary)];
        for region in &summary.regions {
            results.push(build_region_coverage_result(&summary, region));
        }

        CoverageCollectorOutput { results, summary }
    }
}

fn build_file_coverage_result(summary: &CoverageFileSummary) -> crate::report::AtomicResult {
    let mut result = ResultBuilder::new(
        "coverage",
        vec![summary.file_id.clone()],
        vec![summary.read_id.clone()],
        Vec::new(),
    );

    let expected_regions_id = result.expected_records(
        "expected_regions",
        "Expected projected regions for this resolved read.",
        summary
            .expected_regions
            .iter()
            .map(|region| {
                serde_json::json!({
                    "region_id": region.region_id,
                    "name": region.name,
                    "region_type": region.region_type,
                    "sequence_type": region.sequence_type,
                    "start": region.start,
                    "stop": region.stop
                })
            })
            .collect::<Vec<_>>(),
    );
    let sampled_id = result.observed_scalar(
        "sampled_count",
        "Number of reads sampled from the FASTQ file.",
        summary.sampled_count,
        Some("count"),
    );
    let covered_count_id = result.observed_scalar(
        "full_read_coverage_count",
        "Number of sampled reads that are long enough to cover the full projected read geometry.",
        summary.full_read_coverage_count,
        Some("count"),
    );
    let covered_fraction_id = result.observed_scalar(
        "full_read_coverage_fraction",
        "Fraction of sampled reads that are long enough to cover the full projected read geometry.",
        summary.full_read_coverage_fraction,
        Some("fraction"),
    );

    let assessment_type = if summary.full_read_coverage_count == summary.sampled_count {
        AssessmentType::Pass
    } else {
        AssessmentType::Warning
    };
    let assessment_code = if summary.full_read_coverage_count == summary.sampled_count {
        "full_read_coverage_complete"
    } else {
        "full_read_coverage_partial"
    };
    let description = if summary.full_read_coverage_count == summary.sampled_count {
        "All sampled reads fully cover the projected read geometry."
    } else {
        "Some sampled reads are shorter than the projected read geometry."
    };
    result.assessment(
        assessment_type,
        assessment_code,
        description,
        vec![expected_regions_id],
        vec![sampled_id, covered_count_id, covered_fraction_id],
    );

    result.build()
}

fn build_region_coverage_result(
    summary: &CoverageFileSummary,
    region: &CoverageRegionSummary,
) -> crate::report::AtomicResult {
    let mut result = ResultBuilder::new(
        "coverage",
        vec![summary.file_id.clone()],
        vec![summary.read_id.clone()],
        vec![region.region_id.clone()],
    );

    let expected_region_id = result.expected_records(
        "expected_region",
        "Expected coordinates and seqspec annotations for this projected region.",
        vec![serde_json::json!({
            "region_id": region.region_id,
            "name": region.name,
            "region_type": region.region_type,
            "sequence_type": region.sequence_type,
            "start": region.start,
            "stop": region.stop
        })],
    );
    let sampled_id = result.observed_scalar(
        "sampled_count",
        "Number of reads sampled from the FASTQ file.",
        summary.sampled_count,
        Some("count"),
    );
    let covered_count_id = result.observed_scalar(
        "covered_count",
        "Number of sampled reads that fully cover this projected region.",
        region.covered_count,
        Some("count"),
    );
    let covered_fraction_id = result.observed_scalar(
        "covered_fraction",
        "Fraction of sampled reads that fully cover this projected region.",
        region.covered_fraction,
        Some("fraction"),
    );

    let assessment_type = if region.covered_count == summary.sampled_count {
        AssessmentType::Pass
    } else {
        AssessmentType::Warning
    };
    let assessment_code = if region.covered_count == summary.sampled_count {
        "region_coverage_complete"
    } else {
        "region_coverage_partial"
    };
    let description = if region.covered_count == summary.sampled_count {
        format!("All sampled reads cover region '{}'.", region.region_id)
    } else {
        format!(
            "Some sampled reads do not extend far enough to cover region '{}'.",
            region.region_id
        )
    };
    result.assessment(
        assessment_type,
        assessment_code,
        description,
        vec![expected_region_id],
        vec![sampled_id, covered_count_id, covered_fraction_id],
    );

    result.build()
}

pub(crate) fn build_overlap_results(
    summaries: &[CoverageFileSummary],
) -> Vec<crate::report::AtomicResult> {
    let mut by_region_id: BTreeMap<String, Vec<OverlapProjection>> = BTreeMap::new();
    let mut by_region_type: BTreeMap<String, Vec<OverlapProjection>> = BTreeMap::new();

    for summary in summaries {
        for region in &summary.regions {
            let projection = OverlapProjection {
                file_id: summary.file_id.clone(),
                read_id: summary.read_id.clone(),
                region_id: region.region_id.clone(),
                region_type: region.region_type.clone(),
                start: region.start,
                stop: region.stop,
                covered_count: region.covered_count,
                covered_fraction: region.covered_fraction,
            };
            by_region_id
                .entry(region.region_id.clone())
                .or_default()
                .push(projection.clone());
            by_region_type
                .entry(region.region_type.clone())
                .or_default()
                .push(projection);
        }
    }

    let mut results = Vec::new();

    for (region_id, projections) in by_region_id {
        if projections.len() > 1 {
            results.push(build_region_id_overlap_result(&region_id, &projections));
        }
    }

    for (region_type, projections) in by_region_type {
        if projections.len() > 1 && SHARED_REGION_TYPES.contains(&region_type.as_str()) {
            results.push(build_region_type_overlap_result(&region_type, &projections));
        }
    }

    results
}

#[derive(Debug, Clone)]
struct OverlapProjection {
    file_id: String,
    read_id: String,
    region_id: String,
    region_type: String,
    start: usize,
    stop: usize,
    covered_count: usize,
    covered_fraction: f64,
}

fn build_region_id_overlap_result(
    region_id: &str,
    projections: &[OverlapProjection],
) -> crate::report::AtomicResult {
    let mut result = ResultBuilder::new(
        "coverage",
        projections
            .iter()
            .map(|projection| projection.file_id.clone())
            .collect(),
        projections
            .iter()
            .map(|projection| projection.read_id.clone())
            .collect(),
        vec![region_id.to_string()],
    );
    let expected_id = result.expected_records(
        "projected_coordinates_by_read",
        "Projected coordinates for this shared region across matched reads.",
        projections
            .iter()
            .map(|projection| {
                serde_json::json!({
                    "file_id": projection.file_id,
                    "read_id": projection.read_id,
                    "region_id": projection.region_id,
                    "start": projection.start,
                    "stop": projection.stop
                })
            })
            .collect::<Vec<_>>(),
    );
    let observed_id = result.observed_records(
        "covered_counts_by_read",
        "Observed coverage counts for this shared region across matched reads.",
        projections
            .iter()
            .map(|projection| {
                serde_json::json!({
                    "file_id": projection.file_id,
                    "read_id": projection.read_id,
                    "covered_count": projection.covered_count,
                    "covered_fraction": projection.covered_fraction
                })
            })
            .collect::<Vec<_>>(),
    );
    result.assessment(
        AssessmentType::Interpretation,
        "shared_region_id_visible_in_multiple_reads",
        format!(
            "Region '{}' is visible in multiple reads. This can reflect intended paired-end overlap or a read-geometry mismatch when observed reads extend farther than expected.",
            region_id
        ),
        vec![expected_id],
        vec![observed_id],
    );
    result.build()
}

fn build_region_type_overlap_result(
    region_type: &str,
    projections: &[OverlapProjection],
) -> crate::report::AtomicResult {
    let mut result = ResultBuilder::new(
        "coverage",
        projections
            .iter()
            .map(|projection| projection.file_id.clone())
            .collect(),
        projections
            .iter()
            .map(|projection| projection.read_id.clone())
            .collect(),
        Vec::new(),
    );
    let expected_id = result.expected_records(
        "projected_coordinates_by_read",
        "Projected coordinates for this shared biological region type across matched reads.",
        projections
            .iter()
            .map(|projection| {
                serde_json::json!({
                    "file_id": projection.file_id,
                    "read_id": projection.read_id,
                    "region_id": projection.region_id,
                    "region_type": projection.region_type,
                    "start": projection.start,
                    "stop": projection.stop
                })
            })
            .collect::<Vec<_>>(),
    );
    let observed_id = result.observed_records(
        "covered_counts_by_read",
        "Observed coverage counts for this shared biological region type across matched reads.",
        projections
            .iter()
            .map(|projection| {
                serde_json::json!({
                    "file_id": projection.file_id,
                    "read_id": projection.read_id,
                    "region_id": projection.region_id,
                    "covered_count": projection.covered_count,
                    "covered_fraction": projection.covered_fraction
                })
            })
            .collect::<Vec<_>>(),
    );
    result.assessment(
        AssessmentType::Interpretation,
        "shared_region_type_visible_in_multiple_reads",
        format!(
            "Region type '{}' is visible in multiple reads. This can reflect intended paired-end overlap or a read-geometry mismatch when observed reads extend farther than expected.",
            region_type
        ),
        vec![expected_id],
        vec![observed_id],
    );
    result.build()
}
