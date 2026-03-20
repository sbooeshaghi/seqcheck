use crate::context::{InputCheck, ResolvedInput};
use anyhow::{anyhow, Context, Result};
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeSet, HashMap};
use std::io::Write;
use std::path::PathBuf;

const REPORT_SCHEMA_VERSION: &str = "0.1.0";

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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ReportMeta {
    pub command: String,
    pub spec: PathBuf,
    pub modality: String,
    pub requested_reads: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Report {
    pub report_schema_version: String,
    pub meta: ReportMeta,
    pub results: Vec<AtomicResult>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct AtomicResult {
    pub check: String,
    pub files: Vec<String>,
    pub reads: Vec<String>,
    pub regions: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub expected: Vec<MetricItem>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub observed: Vec<MetricItem>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub assessment: Vec<Assessment>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MetricItem {
    pub id: String,
    pub name: String,
    pub description: String,
    pub data: MetricData,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MetricData {
    pub kind: MetricDataKind,
    pub value: Value,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum MetricDataKind {
    Scalar,
    Series,
    Records,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum AssessmentType {
    Pass,
    Warning,
    Error,
    Interpretation,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Assessment {
    #[serde(rename = "type")]
    pub assessment_type: AssessmentType,
    pub code: String,
    pub description: String,
    pub expected_ids: Vec<String>,
    pub observed_ids: Vec<String>,
}

#[derive(Debug)]
pub struct ResultBuilder {
    result: AtomicResult,
    next_expected: usize,
    next_observed: usize,
}

impl Report {
    pub fn new(command: &str, spec: PathBuf, modality: &str, requested_reads: usize) -> Self {
        Self {
            report_schema_version: REPORT_SCHEMA_VERSION.to_string(),
            meta: ReportMeta {
                command: command.to_string(),
                spec,
                modality: modality.to_string(),
                requested_reads,
            },
            results: Vec::new(),
        }
    }

    pub fn validate(&self) -> Result<()> {
        for result in &self.results {
            validate_result(result)?;
        }
        Ok(())
    }
}

impl ResultBuilder {
    pub fn new(check: &str, files: Vec<String>, reads: Vec<String>, regions: Vec<String>) -> Self {
        Self {
            result: AtomicResult {
                check: check.to_string(),
                files,
                reads,
                regions,
                expected: Vec::new(),
                observed: Vec::new(),
                assessment: Vec::new(),
            },
            next_expected: 1,
            next_observed: 1,
        }
    }

    pub fn expected_scalar<T: Serialize>(
        &mut self,
        name: &str,
        description: &str,
        value: T,
        unit: Option<&str>,
    ) -> String {
        self.push_metric(true, name, description, MetricDataKind::Scalar, value, unit)
    }

    pub fn observed_scalar<T: Serialize>(
        &mut self,
        name: &str,
        description: &str,
        value: T,
        unit: Option<&str>,
    ) -> String {
        self.push_metric(
            false,
            name,
            description,
            MetricDataKind::Scalar,
            value,
            unit,
        )
    }

    pub fn expected_records<T: Serialize>(
        &mut self,
        name: &str,
        description: &str,
        value: T,
    ) -> String {
        self.push_metric(
            true,
            name,
            description,
            MetricDataKind::Records,
            value,
            None,
        )
    }

    pub fn observed_records<T: Serialize>(
        &mut self,
        name: &str,
        description: &str,
        value: T,
    ) -> String {
        self.push_metric(
            false,
            name,
            description,
            MetricDataKind::Records,
            value,
            None,
        )
    }

    pub fn expected_series<K, V, I>(
        &mut self,
        name: &str,
        description: &str,
        values: I,
        unit: Option<&str>,
    ) -> String
    where
        K: ToString,
        V: Serialize,
        I: IntoIterator<Item = (K, V)>,
    {
        self.push_metric(
            true,
            name,
            description,
            MetricDataKind::Series,
            build_series(values),
            unit,
        )
    }

    pub fn observed_series<K, V, I>(
        &mut self,
        name: &str,
        description: &str,
        values: I,
        unit: Option<&str>,
    ) -> String
    where
        K: ToString,
        V: Serialize,
        I: IntoIterator<Item = (K, V)>,
    {
        self.push_metric(
            false,
            name,
            description,
            MetricDataKind::Series,
            build_series(values),
            unit,
        )
    }

    pub fn assessment(
        &mut self,
        assessment_type: AssessmentType,
        code: &str,
        description: impl Into<String>,
        expected_ids: Vec<String>,
        observed_ids: Vec<String>,
    ) {
        self.result.assessment.push(Assessment {
            assessment_type,
            code: code.to_string(),
            description: description.into(),
            expected_ids,
            observed_ids,
        });
    }

    pub fn build(self) -> AtomicResult {
        self.result
    }

    fn push_metric<T: Serialize>(
        &mut self,
        expected: bool,
        name: &str,
        description: &str,
        kind: MetricDataKind,
        value: T,
        unit: Option<&str>,
    ) -> String {
        let id = if expected {
            let id = format!("e{}", self.next_expected);
            self.next_expected += 1;
            id
        } else {
            let id = format!("o{}", self.next_observed);
            self.next_observed += 1;
            id
        };

        let item = MetricItem {
            id: id.clone(),
            name: name.to_string(),
            description: description.to_string(),
            data: MetricData {
                kind,
                value: serde_json::to_value(value).expect("metric item should serialize"),
                unit: unit.map(ToOwned::to_owned),
            },
        };

        if expected {
            self.result.expected.push(item);
        } else {
            self.result.observed.push(item);
        }

        id
    }
}

pub fn write_report(output: &Option<PathBuf>, format: OutputFormat, report: &Report) -> Result<()> {
    report.validate()?;

    let payload = match format {
        OutputFormat::Text => render_report(report),
        OutputFormat::Json => serde_json::to_string_pretty(report)?,
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

pub fn input_check_result(input_check: &InputCheck, inputs: &[ResolvedInput]) -> AtomicResult {
    let matched_files = input_check
        .matched_inputs
        .iter()
        .map(|matched| matched.file_id.clone())
        .collect::<Vec<_>>();
    let matched_reads = input_check
        .matched_inputs
        .iter()
        .map(|matched| matched.read_id.clone())
        .collect::<Vec<_>>();

    let mut builder = ResultBuilder::new("input_check", matched_files, matched_reads, Vec::new());

    let expected_ids = builder.expected_records(
        "expected_file_mappings",
        "Expected seqspec file and read mappings for the requested modality.",
        input_check
            .expected_files
            .iter()
            .map(|expected| {
                json!({
                    "file_id": expected.file_id,
                    "read_id": expected.read_id,
                    "filename": expected.filename,
                    "url_basename": expected.url_basename
                })
            })
            .collect::<Vec<_>>(),
    );
    let supplied_ids = builder.observed_records(
        "supplied_input_paths",
        "FASTQ paths supplied to the command.",
        input_check
            .supplied_inputs
            .iter()
            .map(|path| json!({ "input_path": path }))
            .collect::<Vec<_>>(),
    );
    let matched_ids = builder.observed_records(
        "matched_inputs",
        "FASTQ inputs resolved to seqspec file and read ids.",
        input_check
            .matched_inputs
            .iter()
            .map(|matched| {
                json!({
                    "input_path": matched.input_path,
                    "file_id": matched.file_id,
                    "read_id": matched.read_id,
                    "matched_by": matched.matched_by
                })
            })
            .collect::<Vec<_>>(),
    );
    let missing_ids = builder.observed_records(
        "missing_expected_file_ids",
        "Expected seqspec file ids that were not supplied in this run.",
        input_check
            .missing_expected_files
            .iter()
            .map(|expected| {
                json!({
                    "file_id": expected.file_id,
                    "read_id": expected.read_id
                })
            })
            .collect::<Vec<_>>(),
    );
    let primer_status_ids = builder.observed_records(
        "primer_classifications",
        "Primer classification state for each matched input.",
        inputs
            .iter()
            .map(|input| {
                json!({
                    "input_path": input.input_path,
                    "file_id": input.file_id(),
                    "read_id": input.read.read_id,
                    "primer_id": input.read.primer_id,
                    "primer_region_id": input.primer_region.region_id,
                    "classification": input.primer_classification.kind,
                    "scannable": input.primer_classification.scannable,
                    "reason": input.primer_classification.reason
                })
            })
            .collect::<Vec<_>>(),
    );

    if input_check.missing_expected_files.is_empty() {
        builder.assessment(
            AssessmentType::Pass,
            "all_expected_files_matched",
            "All expected modality files were supplied and matched uniquely.",
            vec![expected_ids.clone()],
            vec![
                supplied_ids.clone(),
                matched_ids.clone(),
                missing_ids.clone(),
            ],
        );
    } else {
        builder.assessment(
            AssessmentType::Warning,
            "missing_expected_files",
            format!(
                "Missing expected modality files were not supplied: {}",
                input_check
                    .missing_expected_files
                    .iter()
                    .map(|expected| format!("{} (read_id {})", expected.file_id, expected.read_id))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            vec![expected_ids.clone()],
            vec![
                supplied_ids.clone(),
                matched_ids.clone(),
                missing_ids.clone(),
            ],
        );
    }

    for input in inputs {
        if !input.primer_classification.scannable {
            builder.assessment(
                AssessmentType::Warning,
                "non_scannable_primer",
                format!(
                    "Read '{}' primer_id '{}' points to primer region '{}' that is not scannable: {}",
                    input.read.read_id,
                    input.read.primer_id,
                    input.primer_region.region_id,
                    input
                        .primer_classification
                        .reason
                        .clone()
                        .unwrap_or_else(|| "unknown reason".to_string())
                ),
                Vec::new(),
                vec![primer_status_ids.clone()],
            );
        }
    }

    builder.build()
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

pub fn render_report(report: &Report) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "seqcheck {}\nspec: {}\nmodality: {}\nrequested_reads: {}\n",
        report.meta.command,
        report.meta.spec.display(),
        report.meta.modality,
        report.meta.requested_reads
    ));

    for result in &report.results {
        out.push_str(&format!("\ncheck: {}\n", result.check));
        append_id_list(&mut out, "files", &result.files);
        append_id_list(&mut out, "reads", &result.reads);
        append_id_list(&mut out, "regions", &result.regions);
        append_assessments(&mut out, &result.assessment);
        append_metrics(&mut out, "expected", &result.expected);
        append_metrics(&mut out, "observed", &result.observed);
    }

    out
}

fn validate_result(result: &AtomicResult) -> Result<()> {
    let mut ids = BTreeSet::new();

    for item in result.expected.iter().chain(result.observed.iter()) {
        if !ids.insert(item.id.clone()) {
            return Err(anyhow!(
                "duplicate metric item id '{}' in result '{}'",
                item.id,
                result.check
            ));
        }
    }

    for assessment in &result.assessment {
        for expected_id in &assessment.expected_ids {
            if !result.expected.iter().any(|item| item.id == *expected_id) {
                return Err(anyhow!(
                    "assessment '{}' in result '{}' references unknown expected id '{}'",
                    assessment.code,
                    result.check,
                    expected_id
                ));
            }
        }

        for observed_id in &assessment.observed_ids {
            if !result.observed.iter().any(|item| item.id == *observed_id) {
                return Err(anyhow!(
                    "assessment '{}' in result '{}' references unknown observed id '{}'",
                    assessment.code,
                    result.check,
                    observed_id
                ));
            }
        }
    }

    Ok(())
}

fn build_series<K, V, I>(values: I) -> Vec<Value>
where
    K: ToString,
    V: Serialize,
    I: IntoIterator<Item = (K, V)>,
{
    values
        .into_iter()
        .map(|(key, value)| {
            json!({
                "key": key.to_string(),
                "value": value
            })
        })
        .collect()
}

fn append_id_list(out: &mut String, label: &str, values: &[String]) {
    if values.is_empty() {
        out.push_str(&format!("{}: (none)\n", label));
    } else {
        out.push_str(&format!("{}: {}\n", label, values.join(", ")));
    }
}

fn append_assessments(out: &mut String, assessments: &[Assessment]) {
    if assessments.is_empty() {
        return;
    }

    out.push_str("assessment:\n");
    for assessment in assessments {
        out.push_str(&format!(
            "  - [{}] {}: {}\n",
            assessment_type_label(&assessment.assessment_type),
            assessment.code,
            assessment.description
        ));
    }
}

fn append_metrics(out: &mut String, label: &str, metrics: &[MetricItem]) {
    if metrics.is_empty() {
        return;
    }

    out.push_str(&format!("{}:\n", label));
    for metric in metrics {
        match metric.data.kind {
            MetricDataKind::Scalar => {
                out.push_str(&format!(
                    "  - {}: {}\n",
                    metric.name,
                    format_scalar_value(&metric.data.value, metric.data.unit.as_deref())
                ));
            }
            MetricDataKind::Series => {
                out.push_str(&format!("  - {}:\n", metric.name));
                append_series_values(out, &metric.data.value, metric.data.unit.as_deref());
            }
            MetricDataKind::Records => {
                out.push_str(&format!("  - {}:\n", metric.name));
                append_record_values(out, &metric.data.value);
            }
        }
    }
}

fn append_series_values(out: &mut String, value: &Value, unit: Option<&str>) {
    let Some(entries) = value.as_array() else {
        out.push_str("    (invalid)\n");
        return;
    };

    if entries.is_empty() {
        out.push_str("    (none)\n");
        return;
    }

    for entry in entries {
        let key = entry
            .get("key")
            .map(format_json_value)
            .unwrap_or_else(|| "?".to_string());
        let value = entry
            .get("value")
            .map(|value| format_scalar_value(value, unit))
            .unwrap_or_else(|| "?".to_string());
        out.push_str(&format!("    {} {}\n", key, value));
    }
}

fn append_record_values(out: &mut String, value: &Value) {
    match value {
        Value::Array(entries) => {
            if entries.is_empty() {
                out.push_str("    (none)\n");
                return;
            }

            for entry in entries {
                out.push_str(&format!("    - {}\n", format_record(entry)));
            }
        }
        Value::Object(_) => {
            out.push_str(&format!("    - {}\n", format_record(value)));
        }
        _ => {
            out.push_str(&format!("    - {}\n", format_json_value(value)));
        }
    }
}

fn format_record(value: &Value) -> String {
    match value {
        Value::Object(map) => {
            let mut keys = map.keys().cloned().collect::<Vec<_>>();
            keys.sort();
            keys.into_iter()
                .map(|key| {
                    let value = map
                        .get(&key)
                        .map(format_json_value)
                        .unwrap_or_else(|| "null".to_string());
                    format!("{}={}", key, value)
                })
                .collect::<Vec<_>>()
                .join(" ")
        }
        _ => format_json_value(value),
    }
}

fn format_scalar_value(value: &Value, unit: Option<&str>) -> String {
    let rendered = format_json_value(value);
    match unit {
        Some(unit) => format!("{} {}", rendered, unit),
        None => rendered,
    }
}

fn format_json_value(value: &Value) -> String {
    match value {
        Value::Null => "null".to_string(),
        Value::Bool(value) => value.to_string(),
        Value::Number(value) => value.to_string(),
        Value::String(value) => value.clone(),
        Value::Array(_) | Value::Object(_) => serde_json::to_string(value).unwrap_or_default(),
    }
}

fn assessment_type_label(assessment_type: &AssessmentType) -> &'static str {
    match assessment_type {
        AssessmentType::Pass => "pass",
        AssessmentType::Warning => "warning",
        AssessmentType::Error => "error",
        AssessmentType::Interpretation => "interpretation",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{ExpectedFile, MatchedInput};

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
    fn test_metric_data_serializes_scalar_series_and_records() {
        let mut builder = ResultBuilder::new(
            "test",
            vec!["file".to_string()],
            vec!["read".to_string()],
            vec![],
        );
        builder.expected_scalar("expected_min_len", "desc", 12usize, Some("bp"));
        builder.observed_series(
            "offset_histogram",
            "desc",
            [("+0", 3usize), ("+1", 1usize)],
            Some("count"),
        );
        builder.observed_records(
            "matched_inputs",
            "desc",
            vec![json!({"file_id": "R1", "read_id": "rna_R1"})],
        );
        let result = builder.build();

        assert_eq!(result.expected[0].data.kind, MetricDataKind::Scalar);
        assert_eq!(result.observed[0].data.kind, MetricDataKind::Series);
        assert_eq!(result.observed[1].data.kind, MetricDataKind::Records);
    }

    #[test]
    fn test_report_validation_rejects_unknown_assessment_ids() {
        let result = AtomicResult {
            check: "length".to_string(),
            files: vec!["R1".to_string()],
            reads: vec!["rna_R1".to_string()],
            regions: Vec::new(),
            expected: vec![MetricItem {
                id: "e1".to_string(),
                name: "expected_min_len".to_string(),
                description: "desc".to_string(),
                data: MetricData {
                    kind: MetricDataKind::Scalar,
                    value: json!(12),
                    unit: Some("bp".to_string()),
                },
            }],
            observed: Vec::new(),
            assessment: vec![Assessment {
                assessment_type: AssessmentType::Warning,
                code: "bad".to_string(),
                description: "desc".to_string(),
                expected_ids: vec!["e1".to_string()],
                observed_ids: vec!["o1".to_string()],
            }],
        };

        let report = Report {
            report_schema_version: REPORT_SCHEMA_VERSION.to_string(),
            meta: ReportMeta {
                command: "length".to_string(),
                spec: PathBuf::from("spec.yaml"),
                modality: "rna".to_string(),
                requested_reads: 10,
            },
            results: vec![result],
        };

        let error = report.validate().unwrap_err().to_string();
        assert!(error.contains("unknown observed id 'o1'"));
    }

    #[test]
    fn test_render_report_orders_results_and_metrics_deterministically() {
        let mut input = ResultBuilder::new(
            "input_check",
            vec!["R1".to_string()],
            vec!["rna_R1".to_string()],
            vec![],
        );
        input.expected_records(
            "expected_file_mappings",
            "desc",
            vec![json!({"file_id": "R1", "read_id": "rna_R1"})],
        );
        input.observed_records(
            "matched_inputs",
            "desc",
            vec![json!({"file_id": "R1", "read_id": "rna_R1"})],
        );
        let mut length = ResultBuilder::new(
            "length",
            vec!["R1".to_string()],
            vec!["rna_R1".to_string()],
            vec![],
        );
        length.expected_scalar("expected_min_len", "desc", 12usize, Some("bp"));
        length.observed_scalar("observed_max_len", "desc", 12usize, Some("bp"));

        let report = Report {
            report_schema_version: REPORT_SCHEMA_VERSION.to_string(),
            meta: ReportMeta {
                command: "length".to_string(),
                spec: PathBuf::from("spec.yaml"),
                modality: "rna".to_string(),
                requested_reads: 10,
            },
            results: vec![input.build(), length.build()],
        };

        let observed = render_report(&report);
        assert!(observed.contains("check: input_check"));
        assert!(observed.contains("check: length"));
        assert!(
            observed.find("check: input_check").unwrap() < observed.find("check: length").unwrap()
        );
        assert!(observed.find("expected:").unwrap() < observed.find("observed:").unwrap());
    }

    #[test]
    fn test_input_check_result_builds_atomic_result() {
        let input_check = InputCheck {
            expected_files: vec![ExpectedFile {
                read_id: "rna_R1".to_string(),
                file_id: "R1.fastq.gz".to_string(),
                filename: "R1.fastq.gz".to_string(),
                url_basename: "R1.fastq.gz".to_string(),
            }],
            supplied_inputs: vec![PathBuf::from("R1.fastq.gz")],
            matched_inputs: vec![MatchedInput {
                input_path: PathBuf::from("R1.fastq.gz"),
                read_id: "rna_R1".to_string(),
                file_id: "R1.fastq.gz".to_string(),
                matched_by: "file_id".to_string(),
            }],
            missing_expected_files: Vec::new(),
        };

        let result = input_check_result(&input_check, &[]);
        assert_eq!(result.check, "input_check");
        assert_eq!(result.files, vec!["R1.fastq.gz".to_string()]);
        assert_eq!(result.reads, vec!["rna_R1".to_string()]);
        assert_eq!(result.expected.len(), 1);
        assert_eq!(result.observed.len(), 4);
        assert_eq!(result.assessment[0].code, "all_expected_files_matched");
    }
}
