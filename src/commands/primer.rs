use crate::context::{
    load_resolved_inputs, LoadedInputs, PrimerClassification, PrimerClassificationKind,
};
use crate::report::{
    input_check_result, write_report, AssessmentType, AtomicResult, Report, ResultBuilder,
};
use crate::scan::{run_collector, FastqCollector};
use crate::sequence::find_all_exact_hits;
use crate::sequence::reverse_complement_sequence;
use crate::CommonMetricArgs;
use anyhow::Result;
use std::collections::HashMap;

const TOP_POSITION_LIMIT: usize = 10;

#[derive(Debug, clap::Args)]
pub struct PrimerArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct PositionCount {
    pub position: usize,
    pub count: usize,
}

#[derive(Debug, Default)]
struct PrimerState {
    forward_start_hit_count: usize,
    forward_internal_hit_count: usize,
    reverse_complement_start_hit_count: usize,
    reverse_complement_internal_hit_count: usize,
    absent_count: usize,
    forward_positions: HashMap<usize, usize>,
    reverse_complement_positions: HashMap<usize, usize>,
}

pub(crate) struct PrimerCollector {
    file_id: String,
    read_id: String,
    primer_id: String,
    primer_region_type: String,
    primer_sequence_type: String,
    primer_sequence: String,
    reverse_complement: String,
    primer_classification: PrimerClassification,
    state: PrimerState,
}

pub fn run(args: &PrimerArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
        args.common.n_reads,
        args.common.auth_profile.as_deref(),
    )?;
    let mut report = Report::new(
        "primer",
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
        results.push(run_collector(input, n_reads, PrimerCollector::new(input))?);
    }

    Ok(results)
}

impl PrimerCollector {
    pub(crate) fn new(input: &crate::context::ResolvedInput) -> Self {
        let primer_region = input.primer_region.clone();
        let primer_sequence = primer_region.sequence.clone();
        Self {
            file_id: input.file_id(),
            read_id: input.read.read_id.clone(),
            primer_id: primer_region.region_id,
            primer_region_type: primer_region.region_type,
            primer_sequence_type: primer_region.sequence_type,
            reverse_complement: reverse_complement_sequence(&primer_sequence),
            primer_sequence,
            primer_classification: input.primer_classification.clone(),
            state: PrimerState::default(),
        }
    }
}

impl FastqCollector for PrimerCollector {
    type Output = AtomicResult;

    fn observe(&mut self, sequence: &str, _len: usize) -> Result<()> {
        if !self.primer_classification.scannable {
            return Ok(());
        }

        let forward_hits = find_all_exact_hits(sequence, &self.primer_sequence);
        let reverse_complement_hits = find_all_exact_hits(sequence, &self.reverse_complement);

        if forward_hits.iter().any(|position| *position == 0) {
            self.state.forward_start_hit_count += 1;
        }
        if forward_hits.iter().any(|position| *position > 0) {
            self.state.forward_internal_hit_count += 1;
        }
        if reverse_complement_hits
            .iter()
            .any(|position| *position == 0)
        {
            self.state.reverse_complement_start_hit_count += 1;
        }
        if reverse_complement_hits.iter().any(|position| *position > 0) {
            self.state.reverse_complement_internal_hit_count += 1;
        }
        if forward_hits.is_empty() && reverse_complement_hits.is_empty() {
            self.state.absent_count += 1;
        }

        for position in forward_hits {
            *self.state.forward_positions.entry(position).or_insert(0) += 1;
        }
        for position in reverse_complement_hits {
            *self
                .state
                .reverse_complement_positions
                .entry(position)
                .or_insert(0) += 1;
        }

        Ok(())
    }

    fn finish(self, sampled_count: usize) -> AtomicResult {
        build_primer_result(
            &self.file_id,
            &self.read_id,
            &self.primer_id,
            &self.primer_region_type,
            &self.primer_sequence_type,
            &self.primer_sequence,
            &self.primer_classification,
            sampled_count,
            &self.state,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn build_primer_result(
    file_id: &str,
    read_id: &str,
    primer_id: &str,
    primer_region_type: &str,
    primer_sequence_type: &str,
    primer_sequence: &str,
    primer_classification: &PrimerClassification,
    sampled_count: usize,
    state: &PrimerState,
) -> crate::report::AtomicResult {
    let mut result = ResultBuilder::new(
        "primer",
        vec![file_id.to_string()],
        vec![read_id.to_string()],
        vec![primer_id.to_string()],
    );

    let region_type_id = result.expected_scalar(
        "primer_region_type",
        "Region type of the seqspec region referenced by primer_id.",
        primer_region_type,
        None,
    );
    let sequence_type_id = result.expected_scalar(
        "primer_sequence_type",
        "Sequence type of the seqspec region referenced by primer_id.",
        primer_sequence_type,
        None,
    );
    let sequence_id = result.expected_scalar(
        "primer_sequence",
        "Sequence stored on the seqspec region referenced by primer_id.",
        primer_sequence,
        Some("bases"),
    );
    let length_id = result.expected_scalar(
        "primer_length",
        "Length of the sequence stored on the seqspec region referenced by primer_id.",
        primer_sequence.len(),
        Some("bp"),
    );
    let sampled_id = result.observed_scalar(
        "sampled_count",
        "Number of reads sampled from the FASTQ file.",
        sampled_count,
        Some("count"),
    );
    let classification_id = result.observed_scalar(
        "primer_classification_kind",
        "Primer classification derived from the referenced seqspec region.",
        primer_kind_label(&primer_classification.kind),
        None,
    );
    let scannable_id = result.observed_scalar(
        "primer_scannable",
        "Whether the primer can be scanned as a concrete fixed motif.",
        primer_classification.scannable,
        None,
    );
    let reason_id = result.observed_scalar(
        "primer_reason",
        "Reason for the computed primer classification when one is available.",
        primer_classification.reason.clone(),
        None,
    );

    let mut observed_ids = vec![sampled_id, classification_id, scannable_id, reason_id];

    if primer_classification.scannable {
        let forward_start_id = result.observed_scalar(
            "forward_start_hit_count",
            "Number of sampled reads with at least one exact forward hit at position 0.",
            state.forward_start_hit_count,
            Some("count"),
        );
        let forward_start_fraction_id = result.observed_scalar(
            "forward_start_hit_fraction",
            "Fraction of sampled reads with at least one exact forward hit at position 0.",
            crate::report::fraction(state.forward_start_hit_count, sampled_count),
            Some("fraction"),
        );
        let forward_internal_id = result.observed_scalar(
            "forward_internal_hit_count",
            "Number of sampled reads with at least one exact forward hit after position 0.",
            state.forward_internal_hit_count,
            Some("count"),
        );
        let forward_internal_fraction_id = result.observed_scalar(
            "forward_internal_hit_fraction",
            "Fraction of sampled reads with at least one exact forward hit after position 0.",
            crate::report::fraction(state.forward_internal_hit_count, sampled_count),
            Some("fraction"),
        );
        let reverse_start_id = result.observed_scalar(
            "reverse_complement_start_hit_count",
            "Number of sampled reads with at least one exact reverse-complement hit at position 0.",
            state.reverse_complement_start_hit_count,
            Some("count"),
        );
        let reverse_start_fraction_id = result.observed_scalar(
            "reverse_complement_start_hit_fraction",
            "Fraction of sampled reads with at least one exact reverse-complement hit at position 0.",
            crate::report::fraction(state.reverse_complement_start_hit_count, sampled_count),
            Some("fraction"),
        );
        let reverse_internal_id = result.observed_scalar(
            "reverse_complement_internal_hit_count",
            "Number of sampled reads with at least one exact reverse-complement hit after position 0.",
            state.reverse_complement_internal_hit_count,
            Some("count"),
        );
        let reverse_internal_fraction_id = result.observed_scalar(
            "reverse_complement_internal_hit_fraction",
            "Fraction of sampled reads with at least one exact reverse-complement hit after position 0.",
            crate::report::fraction(state.reverse_complement_internal_hit_count, sampled_count),
            Some("fraction"),
        );
        let absent_id = result.observed_scalar(
            "absent_count",
            "Number of sampled reads with no exact forward or reverse-complement primer hit.",
            state.absent_count,
            Some("count"),
        );
        let absent_fraction_id = result.observed_scalar(
            "absent_fraction",
            "Fraction of sampled reads with no exact forward or reverse-complement primer hit.",
            crate::report::fraction(state.absent_count, sampled_count),
            Some("fraction"),
        );
        let forward_positions_id = result.observed_records(
            "forward_hit_positions",
            "Most frequent forward primer hit positions.",
            top_positions(&state.forward_positions, TOP_POSITION_LIMIT),
        );
        let reverse_positions_id = result.observed_records(
            "reverse_complement_hit_positions",
            "Most frequent reverse-complement primer hit positions.",
            top_positions(&state.reverse_complement_positions, TOP_POSITION_LIMIT),
        );
        observed_ids.extend([
            forward_start_id,
            forward_start_fraction_id,
            forward_internal_id,
            forward_internal_fraction_id,
            reverse_start_id,
            reverse_start_fraction_id,
            reverse_internal_id,
            reverse_internal_fraction_id,
            absent_id,
            absent_fraction_id,
            forward_positions_id,
            reverse_positions_id,
        ]);
    }

    let expected_ids = vec![region_type_id, sequence_type_id, sequence_id, length_id];

    match primer_classification.kind {
        PrimerClassificationKind::FixedScannable => {
            if state.absent_count == sampled_count {
                result.assessment(
                    AssessmentType::Pass,
                    "primer_absent_from_reads",
                    "The primer motif is absent from the sampled reads, which is consistent with sequencing starting at the primer.",
                    expected_ids,
                    observed_ids,
                );
            } else {
                result.assessment(
                    AssessmentType::Interpretation,
                    "primer_sequence_detected_in_reads",
                    "The primer motif appears in some sampled reads. This can indicate shifted geometry, readthrough, or primer sequence carried into the read.",
                    expected_ids,
                    observed_ids,
                );
            }
        }
        PrimerClassificationKind::GhostPrimer => {
            result.assessment(
                AssessmentType::Interpretation,
                "ghost_primer_anchor",
                "The read is anchored by a zero-length ghost primer. This is allowed for projection but cannot be scanned as a sequence motif.",
                expected_ids,
                observed_ids,
            );
        }
        PrimerClassificationKind::NonScannablePrimer => {
            result.assessment(
                AssessmentType::Warning,
                "non_scannable_primer",
                "The region referenced by primer_id is not a concrete fixed primer and cannot be scanned as an anchor motif.",
                expected_ids,
                observed_ids,
            );
        }
    }

    result.build()
}

fn top_positions(counts: &HashMap<usize, usize>, limit: usize) -> Vec<PositionCount> {
    let mut items = counts
        .iter()
        .map(|(position, count)| PositionCount {
            position: *position,
            count: *count,
        })
        .collect::<Vec<_>>();

    items.sort_by(|left, right| {
        right
            .count
            .cmp(&left.count)
            .then_with(|| left.position.cmp(&right.position))
    });
    items.truncate(limit);
    items
}

fn primer_kind_label(kind: &PrimerClassificationKind) -> &'static str {
    match kind {
        PrimerClassificationKind::FixedScannable => "fixed_scannable",
        PrimerClassificationKind::GhostPrimer => "ghost_primer",
        PrimerClassificationKind::NonScannablePrimer => "non_scannable_primer",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_top_positions_sort_by_count_then_position() {
        let counts = HashMap::from([(2usize, 1usize), (1usize, 2usize), (0usize, 2usize)]);
        let observed = top_positions(&counts, 3);

        assert_eq!(
            observed,
            vec![
                PositionCount {
                    position: 0,
                    count: 2,
                },
                PositionCount {
                    position: 1,
                    count: 2,
                },
                PositionCount {
                    position: 2,
                    count: 1,
                },
            ]
        );
    }
}
