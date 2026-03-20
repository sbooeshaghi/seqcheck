use crate::context::load_resolved_inputs;
use crate::report::{
    format_fraction, render_report_prelude, write_report, FileReport, ReportEnvelope,
};
use crate::scan::{is_region_covered, scan_fastq};
use crate::CommonMetricArgs;
use anyhow::Result;
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Debug, clap::Args)]
pub struct CoverageArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Serialize)]
pub struct CoverageRegionResult {
    pub region_id: String,
    pub name: String,
    pub region_type: String,
    pub sequence_type: String,
    pub start: usize,
    pub stop: usize,
    pub covered_count: usize,
    pub covered_fraction: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct CoverageResult {
    pub sampled_count: usize,
    pub expected_regions: Vec<crate::context::ExpectedRegion>,
    pub full_read_coverage_count: usize,
    pub full_read_coverage_fraction: f64,
    pub regions: Vec<CoverageRegionResult>,
}

#[derive(Debug)]
struct CoverageState {
    full_read_coverage_count: usize,
    region_counts: Vec<usize>,
}

pub fn run(args: &CoverageArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
    )?;
    let crate::context::LoadedInputs {
        input_check,
        warnings: base_warnings,
        inputs,
        ..
    } = loaded;
    let mut files = Vec::new();

    for input in inputs {
        let input_path = input.input_path.clone();
        let read_id = input.read.read_id.clone();
        let file_id = input.file_id();
        let matched_by = input.matched_by.clone();
        let expected_regions = input.expected_regions();
        let expected_stop = input.expected_stop();
        let state = CoverageState {
            full_read_coverage_count: 0,
            region_counts: vec![0; input.coordinates.len()],
        };
        let coordinates = input.coordinates.clone();

        let (state, sampled_count) =
            scan_fastq(&input, args.common.n_reads, state, |state, _, _, len| {
                if len >= expected_stop {
                    state.full_read_coverage_count += 1;
                }

                for (idx, region) in coordinates.iter().enumerate() {
                    if is_region_covered(len, region) {
                        state.region_counts[idx] += 1;
                    }
                }

                Ok(())
            })?;

        let regions = input
            .coordinates
            .iter()
            .enumerate()
            .map(|(idx, region)| CoverageRegionResult {
                region_id: region.region.region_id.clone(),
                name: region.region.name.clone(),
                region_type: region.region.region_type.clone(),
                sequence_type: region.region.sequence_type.clone(),
                start: usize::try_from(region.start).unwrap_or_default(),
                stop: usize::try_from(region.stop).unwrap_or_default(),
                covered_count: state.region_counts[idx],
                covered_fraction: crate::report::fraction(state.region_counts[idx], sampled_count),
            })
            .collect();

        files.push(FileReport {
            input_path,
            read_id,
            file_id,
            matched_by,
            results: CoverageResult {
                sampled_count,
                expected_regions,
                full_read_coverage_count: state.full_read_coverage_count,
                full_read_coverage_fraction: crate::report::fraction(
                    state.full_read_coverage_count,
                    sampled_count,
                ),
                regions,
            },
        });
    }

    let mut warnings = base_warnings;
    warnings.extend(build_assignment_warnings(&files));

    let report = ReportEnvelope {
        spec: args.common.spec.clone(),
        modality: args.common.modality.clone(),
        command: "coverage".to_string(),
        n_reads: args.common.n_reads,
        input_check,
        warnings,
        files,
    };

    write_report(
        &args.common.output,
        args.common.format,
        &report,
        render_text(&report),
    )
}

fn render_text(report: &ReportEnvelope<CoverageResult>) -> String {
    let mut out = render_report_prelude("coverage", report);

    for file in &report.files {
        let results = &file.results;
        out.push_str(&format!(
            "\nfile: {}\nread_id: {}\nfile_id: {}\nmatched_by: {}\nsampled_reads: {}\nfull_read_coverage: {} ({})\nexpected_regions:\n",
            file.input_path.display(),
            file.read_id,
            file.file_id,
            file.matched_by,
            results.sampled_count,
            results.full_read_coverage_count,
            format_fraction(results.full_read_coverage_count, results.sampled_count),
        ));
        for region in &results.regions {
            out.push_str(&format!(
                "  - {} [{}:{}] {} {} covered {} ({})\n",
                region.region_id,
                region.start,
                region.stop,
                region.region_type,
                region.sequence_type,
                region.covered_count,
                format_fraction(region.covered_count, results.sampled_count),
            ));
        }
    }

    out
}

fn build_assignment_warnings(files: &[FileReport<CoverageResult>]) -> Vec<String> {
    let mut by_region_id: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut by_region_type: BTreeMap<String, Vec<String>> = BTreeMap::new();

    for file in files {
        let location = format!("{}:{}", file.read_id, file.file_id);
        for region in &file.results.regions {
            by_region_id
                .entry(region.region_id.clone())
                .or_default()
                .push(location.clone());
            by_region_type
                .entry(region.region_type.clone())
                .or_default()
                .push(location.clone());
        }
    }

    let mut warnings = Vec::new();

    for (region_id, locations) in by_region_id {
        if locations.len() > 1 {
            warnings.push(format!(
                "region_id '{}' appears in multiple reads/files: {}. This can reflect intended overlapping paired-end reads, or a seqspec/read-geometry mismatch when observed reads extend farther than expected.",
                region_id,
                locations.join(", ")
            ));
        }
    }

    for (region_type, locations) in by_region_type {
        if locations.len() > 1
            && matches!(
                region_type.as_str(),
                "barcode" | "umi" | "cdna" | "gdna" | "protein" | "tag" | "sgrna_target"
            )
        {
            warnings.push(format!(
                "region_type '{}' appears in multiple reads/files: {}. This can reflect intended overlapping paired-end reads, or a seqspec/read-geometry mismatch when observed reads extend farther than expected.",
                region_type,
                locations.join(", ")
            ));
        }
    }

    warnings
}
