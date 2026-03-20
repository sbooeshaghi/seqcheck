use crate::context::{filter_regions_by_sequence_type, load_resolved_inputs};
use crate::report::{
    input_check_result, top_sequences, write_report, AssessmentType, Report, ResultBuilder,
};
use crate::scan::{extract_region, scan_fastq};
use crate::CommonMetricArgs;
use anyhow::Result;
use std::collections::HashMap;

const DNA_BITS_PER_BASE: f64 = 2.0;
const TOP_SEQUENCE_LIMIT: usize = 10;

#[derive(Debug, clap::Args)]
pub struct RandomArgs {
    #[command(flatten)]
    pub common: CommonMetricArgs,
}

#[derive(Debug, Clone, Default)]
struct RandomRegionState {
    covered_count: usize,
    short_read_count: usize,
    sequences: HashMap<String, usize>,
}

pub fn run(args: &RandomArgs) -> Result<()> {
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
        "random",
        args.common.spec.clone(),
        &args.common.modality,
        args.common.n_reads,
    );
    report
        .results
        .push(input_check_result(&input_check, &inputs));

    for input in inputs {
        let random_regions = filter_regions_by_sequence_type(&input, "random");
        let state = vec![RandomRegionState::default(); random_regions.len()];
        let regions = random_regions.clone();

        let (state, sampled_count) = scan_fastq(
            &input,
            args.common.n_reads,
            state,
            |state, _, sequence, _| {
                for (idx, region) in regions.iter().enumerate() {
                    match extract_region(sequence, region) {
                        Some(observed) => {
                            state[idx].covered_count += 1;
                            *state[idx].sequences.entry(observed).or_insert(0) += 1;
                        }
                        None => {
                            state[idx].short_read_count += 1;
                        }
                    }
                }
                Ok(())
            },
        )?;

        for (idx, region) in random_regions.iter().enumerate() {
            report.results.push(build_random_result(
                &input.file_id(),
                &input.read.read_id,
                sampled_count,
                region,
                &state[idx],
            ));
        }
    }

    write_report(&args.common.output, args.common.format, &report)
}

fn build_random_result(
    file_id: &str,
    read_id: &str,
    sampled_count: usize,
    region: &seqspec::region::RegionCoordinate,
    state: &RandomRegionState,
) -> crate::report::AtomicResult {
    let region_len = usize::try_from(region.stop - region.start).unwrap_or_default();
    let max_entropy_bits = region_len as f64 * DNA_BITS_PER_BASE;
    let sequence_entropy_bits = sequence_entropy_bits(&state.sequences);
    let sequence_entropy_fraction = if max_entropy_bits == 0.0 {
        0.0
    } else {
        sequence_entropy_bits / max_entropy_bits
    };

    let mut result = ResultBuilder::new(
        "random",
        vec![file_id.to_string()],
        vec![read_id.to_string()],
        vec![region.region.region_id.clone()],
    );

    let region_id = result.expected_records(
        "expected_region",
        "Expected projected coordinates and annotations for this random region.",
        vec![serde_json::json!({
            "start": usize::try_from(region.start).unwrap_or_default(),
            "stop": usize::try_from(region.stop).unwrap_or_default(),
            "region_type": region.region.region_type,
            "name": region.region.name
        })],
    );
    let max_entropy_id = result.expected_scalar(
        "max_entropy_bits",
        "Theoretical maximum Shannon entropy for a random DNA sequence of this length.",
        max_entropy_bits,
        Some("bits"),
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
    let unique_id = result.observed_scalar(
        "unique_sequence_count",
        "Number of unique extracted sequences observed for this region.",
        state.sequences.len(),
        Some("count"),
    );
    let entropy_id = result.observed_scalar(
        "sequence_entropy_bits",
        "Observed Shannon entropy of the exact extracted sequence distribution.",
        sequence_entropy_bits,
        Some("bits"),
    );
    let entropy_fraction_id = result.observed_scalar(
        "sequence_entropy_fraction",
        "Observed Shannon entropy divided by the theoretical maximum for this region length.",
        sequence_entropy_fraction,
        Some("fraction"),
    );
    let top_sequences_id = result.observed_records(
        "top_sequences",
        "Most frequent exact sequences observed for this region.",
        top_sequences(&state.sequences, TOP_SEQUENCE_LIMIT),
    );

    result.assessment(
        AssessmentType::Interpretation,
        "random_sequence_entropy",
        format!(
            "Observed sequence entropy for region '{}' is {:.4} of the theoretical maximum.",
            region.region.region_id, sequence_entropy_fraction
        ),
        vec![region_id, max_entropy_id],
        vec![
            sampled_id,
            covered_id,
            covered_fraction_id,
            short_id,
            unique_id,
            entropy_id,
            entropy_fraction_id,
            top_sequences_id,
        ],
    );

    result.build()
}

fn sequence_entropy_bits(sequences: &HashMap<String, usize>) -> f64 {
    let total = sequences.values().sum::<usize>() as f64;
    if total == 0.0 {
        return 0.0;
    }

    sequences
        .values()
        .filter(|count| **count > 0)
        .map(|count| {
            let probability = *count as f64 / total;
            -(probability * probability.log2())
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_entropy_bits_uses_exact_sequence_distribution() {
        let mut sequences = HashMap::new();
        sequences.insert("AAA".to_string(), 2);
        sequences.insert("AAT".to_string(), 1);

        let entropy = sequence_entropy_bits(&sequences);

        assert!((entropy - 0.9182958340544896).abs() < 1e-12);
    }

    #[test]
    fn test_sequence_entropy_bits_is_zero_for_empty_distribution() {
        let entropy = sequence_entropy_bits(&HashMap::new());
        assert_eq!(entropy, 0.0);
    }
}
