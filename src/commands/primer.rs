use crate::context::{load_resolved_inputs, PrimerClassification, PrimerClassificationKind};
use crate::report::{
    format_fraction, render_report_prelude, write_report, FileReport, ReportEnvelope,
};
use crate::scan::scan_fastq;
use crate::sequence::find_all_exact_hits;
use crate::sequence::reverse_complement_sequence;
use crate::CommonMetricArgs;
use anyhow::Result;
use serde::Serialize;
use std::collections::HashMap;

const TOP_POSITION_LIMIT: usize = 10;

#[derive(Debug, clap::Args)]
pub struct PrimerArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct PositionCount {
    pub position: usize,
    pub count: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct PrimerResult {
    pub primer_id: String,
    pub primer_region_type: String,
    pub primer_sequence_type: String,
    pub primer_sequence: String,
    pub primer_length: usize,
    pub primer_classification: PrimerClassification,
    pub sampled_count: usize,
    pub forward_start_hit_count: usize,
    pub forward_start_hit_fraction: f64,
    pub forward_internal_hit_count: usize,
    pub forward_internal_hit_fraction: f64,
    pub reverse_complement_start_hit_count: usize,
    pub reverse_complement_start_hit_fraction: f64,
    pub reverse_complement_internal_hit_count: usize,
    pub reverse_complement_internal_hit_fraction: f64,
    pub absent_count: usize,
    pub absent_fraction: f64,
    pub forward_hit_positions: Vec<PositionCount>,
    pub reverse_complement_hit_positions: Vec<PositionCount>,
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

pub fn run(args: &PrimerArgs) -> Result<()> {
    let loaded = load_resolved_inputs(
        &args.common.spec,
        &args.common.modality,
        &args.common.fastqs,
    )?;
    let crate::context::LoadedInputs {
        input_check,
        warnings,
        inputs,
        ..
    } = loaded;
    let mut files = Vec::new();

    for input in inputs {
        let input_path = input.input_path.clone();
        let read_id = input.read.read_id.clone();
        let file_id = input.file_id();
        let matched_by = input.matched_by.clone();
        let primer_region = input.primer_region.clone();
        let primer_classification = input.primer_classification.clone();
        let primer_sequence = primer_region.sequence.clone();
        let reverse_complement = reverse_complement_sequence(&primer_sequence);
        let scannable = primer_classification.scannable;

        let (state, sampled_count) = scan_fastq(
            &input,
            args.common.n_reads,
            PrimerState::default(),
            |state, _, sequence, _| {
                if !scannable {
                    return Ok(());
                }

                let forward_hits = find_all_exact_hits(sequence, &primer_sequence);
                let reverse_complement_hits = find_all_exact_hits(sequence, &reverse_complement);

                if forward_hits.iter().any(|position| *position == 0) {
                    state.forward_start_hit_count += 1;
                }
                if forward_hits.iter().any(|position| *position > 0) {
                    state.forward_internal_hit_count += 1;
                }
                if reverse_complement_hits
                    .iter()
                    .any(|position| *position == 0)
                {
                    state.reverse_complement_start_hit_count += 1;
                }
                if reverse_complement_hits.iter().any(|position| *position > 0) {
                    state.reverse_complement_internal_hit_count += 1;
                }
                if forward_hits.is_empty() && reverse_complement_hits.is_empty() {
                    state.absent_count += 1;
                }

                for position in forward_hits {
                    *state.forward_positions.entry(position).or_insert(0) += 1;
                }
                for position in reverse_complement_hits {
                    *state
                        .reverse_complement_positions
                        .entry(position)
                        .or_insert(0) += 1;
                }

                Ok(())
            },
        )?;

        files.push(FileReport {
            input_path,
            read_id,
            file_id,
            matched_by,
            results: PrimerResult {
                primer_id: primer_region.region_id.clone(),
                primer_region_type: primer_region.region_type.clone(),
                primer_sequence_type: primer_region.sequence_type.clone(),
                primer_sequence,
                primer_length: usize::try_from(primer_region.max_len.max(0)).unwrap_or_default(),
                primer_classification,
                sampled_count,
                forward_start_hit_count: state.forward_start_hit_count,
                forward_start_hit_fraction: crate::report::fraction(
                    state.forward_start_hit_count,
                    sampled_count,
                ),
                forward_internal_hit_count: state.forward_internal_hit_count,
                forward_internal_hit_fraction: crate::report::fraction(
                    state.forward_internal_hit_count,
                    sampled_count,
                ),
                reverse_complement_start_hit_count: state.reverse_complement_start_hit_count,
                reverse_complement_start_hit_fraction: crate::report::fraction(
                    state.reverse_complement_start_hit_count,
                    sampled_count,
                ),
                reverse_complement_internal_hit_count: state.reverse_complement_internal_hit_count,
                reverse_complement_internal_hit_fraction: crate::report::fraction(
                    state.reverse_complement_internal_hit_count,
                    sampled_count,
                ),
                absent_count: state.absent_count,
                absent_fraction: crate::report::fraction(state.absent_count, sampled_count),
                forward_hit_positions: top_positions(&state.forward_positions, TOP_POSITION_LIMIT),
                reverse_complement_hit_positions: top_positions(
                    &state.reverse_complement_positions,
                    TOP_POSITION_LIMIT,
                ),
            },
        });
    }

    let report = ReportEnvelope {
        spec: args.common.spec.clone(),
        modality: args.common.modality.clone(),
        command: "primer".to_string(),
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

fn render_text(report: &ReportEnvelope<PrimerResult>) -> String {
    let mut out = render_report_prelude("primer", report);

    for file in &report.files {
        out.push_str(&format!(
            "\nfile: {}\nread_id: {}\nfile_id: {}\nmatched_by: {}\nsampled_reads: {}\nprimer_id: {}\nprimer_region_type: {}\nprimer_sequence_type: {}\nprimer_length: {}\nprimer_classification: {:?}\nscannable: {}\n",
            file.input_path.display(),
            file.read_id,
            file.file_id,
            file.matched_by,
            file.results.sampled_count,
            file.results.primer_id,
            file.results.primer_region_type,
            file.results.primer_sequence_type,
            file.results.primer_length,
            primer_kind_label(&file.results.primer_classification.kind),
            file.results.primer_classification.scannable,
        ));

        if let Some(reason) = &file.results.primer_classification.reason {
            out.push_str(&format!("primer_reason: {}\n", reason));
        }
        if file.results.primer_classification.scannable {
            out.push_str(&format!(
                "primer_sequence: {}\nforward_start_hits: {} ({})\nforward_internal_hits: {} ({})\nreverse_complement_start_hits: {} ({})\nreverse_complement_internal_hits: {} ({})\nprimer_absent_reads: {} ({})\n",
                file.results.primer_sequence,
                file.results.forward_start_hit_count,
                format_fraction(
                    file.results.forward_start_hit_count,
                    file.results.sampled_count,
                ),
                file.results.forward_internal_hit_count,
                format_fraction(
                    file.results.forward_internal_hit_count,
                    file.results.sampled_count,
                ),
                file.results.reverse_complement_start_hit_count,
                format_fraction(
                    file.results.reverse_complement_start_hit_count,
                    file.results.sampled_count,
                ),
                file.results.reverse_complement_internal_hit_count,
                format_fraction(
                    file.results.reverse_complement_internal_hit_count,
                    file.results.sampled_count,
                ),
                file.results.absent_count,
                format_fraction(file.results.absent_count, file.results.sampled_count),
            ));
            render_positions(
                &mut out,
                "forward_hit_positions",
                &file.results.forward_hit_positions,
            );
            render_positions(
                &mut out,
                "reverse_complement_hit_positions",
                &file.results.reverse_complement_hit_positions,
            );
        }
    }

    out
}

fn render_positions(out: &mut String, label: &str, positions: &[PositionCount]) {
    out.push_str(&format!("{}:\n", label));
    if positions.is_empty() {
        out.push_str("  (none)\n");
    } else {
        for position in positions {
            out.push_str(&format!("  {} {}\n", position.position, position.count));
        }
    }
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
