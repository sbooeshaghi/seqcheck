use crate::context::ResolvedInput;
use anyhow::{Context, Result};
use seqspec::region::RegionCoordinate;
use seqspec::utils::is_remote_source;

pub trait FastqCollector {
    type Output;

    fn observe(&mut self, sequence: &str, len: usize) -> Result<()>;

    fn finish(self, sampled_count: usize) -> Self::Output;
}

pub fn scan_fastq<S, F>(
    input: &ResolvedInput,
    n_reads: usize,
    mut state: S,
    mut update: F,
) -> Result<(S, usize)>
where
    F: FnMut(&mut S, usize, &str, usize) -> Result<()>,
{
    let sampled = scan_fastq_records(input, n_reads, |record_idx, sequence, len| {
        update(&mut state, record_idx, sequence, len)
    })?;

    Ok((state, sampled))
}

pub fn scan_fastq_records<F>(input: &ResolvedInput, n_reads: usize, mut observe: F) -> Result<usize>
where
    F: FnMut(usize, &str, usize) -> Result<()>,
{
    let limit = if n_reads == 0 { usize::MAX } else { n_reads };
    let mut records = if is_remote_source(&input.input_source) {
        input
            .remote_access
            .with_reader(&input.input_source, |reader| {
                kseq::parse_reader(reader)
                    .with_context(|| format!("failed to open FASTQ {}", input.input_source))
            })?
    } else {
        kseq::parse_path(&input.input_path)
            .with_context(|| format!("failed to open FASTQ {}", input.input_source))?
    };
    let mut sampled = 0usize;

    while sampled < limit {
        let Some(record) = records.iter_record()? else {
            break;
        };

        sampled += 1;
        let sequence = record.seq();
        observe(sampled, sequence, sequence.len())?;
    }

    Ok(sampled)
}

pub fn run_collector<C>(
    input: &ResolvedInput,
    n_reads: usize,
    mut collector: C,
) -> Result<C::Output>
where
    C: FastqCollector,
{
    let sampled_count = scan_fastq_records(input, n_reads, |_, sequence, len| {
        collector.observe(sequence, len)
    })?;
    Ok(collector.finish(sampled_count))
}

pub fn region_stop(region: &RegionCoordinate) -> usize {
    usize::try_from(region.stop).unwrap_or_default()
}

pub fn region_start(region: &RegionCoordinate) -> usize {
    usize::try_from(region.start).unwrap_or_default()
}

pub fn is_region_covered(read_len: usize, region: &RegionCoordinate) -> bool {
    read_len >= region_stop(region)
}

pub fn extract_region(sequence: &str, region: &RegionCoordinate) -> Option<String> {
    let start = region_start(region);
    let stop = region_stop(region);
    if stop > sequence.len() || start > stop {
        return None;
    }
    Some(sequence[start..stop].to_string())
}
