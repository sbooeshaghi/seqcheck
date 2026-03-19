use anyhow::{Context, Result};
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize, ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum OutputFormat {
    #[default]
    Text,
    Json,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SequenceCount {
    pub sequence: String,
    pub count: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct FileReport<T> {
    pub input_path: PathBuf,
    pub read_id: String,
    pub file_id: String,
    pub matched_by: String,
    pub results: T,
}

#[derive(Debug, Clone, Serialize)]
pub struct ReportEnvelope<T> {
    pub spec: PathBuf,
    pub modality: String,
    pub command: String,
    pub n_reads: usize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
    pub files: Vec<FileReport<T>>,
}

pub fn write_report<T: Serialize>(
    output: &Option<PathBuf>,
    format: OutputFormat,
    envelope: &ReportEnvelope<T>,
    text: String,
) -> Result<()> {
    let payload = match format {
        OutputFormat::Text => text,
        OutputFormat::Json => serde_json::to_string_pretty(envelope)?,
    };

    match output {
        Some(path) => {
            let mut file = std::fs::File::create(path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            file.write_all(payload.as_bytes())?;
            if !payload.ends_with('\n') {
                file.write_all(b"\n")?;
            }
        }
        None => {
            let mut stdout = std::io::stdout();
            stdout.write_all(payload.as_bytes())?;
            if !payload.ends_with('\n') {
                stdout.write_all(b"\n")?;
            }
        }
    }

    Ok(())
}

pub fn top_sequences(counts: &HashMap<String, usize>, limit: usize) -> Vec<SequenceCount> {
    let mut items: Vec<SequenceCount> = counts
        .iter()
        .map(|(sequence, count)| SequenceCount {
            sequence: sequence.clone(),
            count: *count,
        })
        .collect();

    items.sort_by(|left, right| {
        right
            .count
            .cmp(&left.count)
            .then_with(|| left.sequence.cmp(&right.sequence))
    });

    items.truncate(limit);
    items
}

pub fn fraction(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

pub fn format_fraction(numerator: usize, denominator: usize) -> String {
    format!("{:.4}", fraction(numerator, denominator))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_top_sequences_are_sorted_and_truncated() {
        let counts = HashMap::from([
            ("CCC".to_string(), 1usize),
            ("AAA".to_string(), 3usize),
            ("BBB".to_string(), 3usize),
        ]);

        let observed = top_sequences(&counts, 2);

        assert_eq!(
            observed,
            vec![
                SequenceCount {
                    sequence: "AAA".to_string(),
                    count: 3,
                },
                SequenceCount {
                    sequence: "BBB".to_string(),
                    count: 3,
                }
            ]
        );
    }

    #[test]
    fn test_report_envelope_serializes() {
        let envelope = ReportEnvelope {
            spec: PathBuf::from("spec.yaml"),
            modality: "rna".to_string(),
            command: "length".to_string(),
            n_reads: 10,
            warnings: Vec::new(),
            files: vec![FileReport {
                input_path: PathBuf::from("R1.fastq.gz"),
                read_id: "rna_R1".to_string(),
                file_id: "R1.fastq.gz".to_string(),
                matched_by: "file_id".to_string(),
                results: serde_json::json!({
                    "sampled_count": 10,
                    "out_of_range_count": 0
                }),
            }],
        };

        let json = serde_json::to_string(&envelope).unwrap();
        assert!(json.contains("\"command\":\"length\""));
        assert!(json.contains("\"matched_by\":\"file_id\""));
    }
}
