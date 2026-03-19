use crate::context::ResolvedInput;
use anyhow::{Context, Result};
use seqspec::region::RegionCoordinate;

pub fn scan_fastq<S, F>(
    input: &ResolvedInput,
    n_reads: usize,
    mut state: S,
    mut update: F,
) -> Result<(S, usize)>
where
    F: FnMut(&mut S, usize, &str, usize) -> Result<()>,
{
    let limit = if n_reads == 0 { usize::MAX } else { n_reads };
    let mut records = kseq::parse_path(&input.input_path)
        .with_context(|| format!("failed to open FASTQ {}", input.input_path.display()))?;
    let mut sampled = 0usize;

    while sampled < limit {
        let Some(record) = records.iter_record()? else {
            break;
        };

        sampled += 1;
        let sequence = record.seq();
        update(&mut state, sampled, sequence, sequence.len())?;
    }

    Ok((state, sampled))
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
