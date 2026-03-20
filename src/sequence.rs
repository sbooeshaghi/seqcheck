pub fn reverse_sequence(sequence: &str) -> String {
    sequence.chars().rev().collect()
}

pub fn complement_sequence(sequence: &str) -> String {
    seqspec::utils::complement_seq(sequence)
}

pub fn reverse_complement_sequence(sequence: &str) -> String {
    reverse_sequence(&complement_sequence(sequence))
}

pub fn find_all_exact_hits(sequence: &str, pattern: &str) -> Vec<usize> {
    if pattern.is_empty() || pattern.len() > sequence.len() {
        return Vec::new();
    }

    let mut hits = Vec::new();
    let mut start = 0usize;

    while start + pattern.len() <= sequence.len() {
        let window = &sequence[start..];
        let Some(offset) = window.find(pattern) else {
            break;
        };
        let position = start + offset;
        hits.push(position);
        start = position + 1;
    }

    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_sequence() {
        assert_eq!(reverse_complement_sequence("ATGC"), "GCAT");
    }

    #[test]
    fn test_find_all_exact_hits_returns_overlapping_hits() {
        let hits = find_all_exact_hits("AAAA", "AA");
        assert_eq!(hits, vec![0, 1, 2]);
    }
}
