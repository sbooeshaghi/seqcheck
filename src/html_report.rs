use crate::report::{AssessmentType, Report};
use anyhow::{Context, Result};
use seqspec::assay::{Assay, LibKit, LibProtocol, SeqKit, SeqProtocol};
use seqspec::file::File as SeqspecFile;
use seqspec::onlist::Onlist as SeqspecOnlist;
use seqspec::region::Region;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

const TEMPLATE_HTML: &str = include_str!("report_assets/template.html");
const STYLE_CSS: &str = include_str!("report_assets/style.css");
const APP_JS: &str = include_str!("report_assets/app.js");
const REPOSITORY_URL: &str = env!("CARGO_PKG_REPOSITORY");
const SEQCHECK_VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SeqspecLibRegion {
    pub region_id: String,
    pub region_type: String,
    pub name: String,
    pub sequence_type: String,
    pub sequence: String,
    pub min_len: i64,
    pub max_len: i64,
    pub len: i64,
    pub bp_start: i64,
    pub bp_end: i64,
    pub depth: usize,
    pub parent_region_id: Option<String>,
    pub path_region_ids: Vec<String>,
    pub path_names: Vec<String>,
    pub is_leaf: bool,
    pub child_region_ids: Vec<String>,
    pub onlist: Option<SeqspecLibOnlist>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SeqspecLibOnlist {
    pub file_id: String,
    pub filename: String,
    pub filetype: String,
    pub filesize: i64,
    pub url: String,
    pub urltype: String,
    pub md5: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SeqspecLibRead {
    pub read_id: String,
    pub name: String,
    pub label: String,
    pub min_len: i64,
    pub max_len: i64,
    pub strand: String,
    pub start: i64,
    pub end: i64,
    pub primer_id: String,
    pub files: Vec<SeqspecLibFile>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SeqspecLibFile {
    pub file_id: String,
    pub filename: String,
    pub filetype: String,
    pub filesize: i64,
    pub url: String,
    pub urltype: String,
    pub md5: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SeqspecMetadataRow {
    pub protocol_id: Option<String>,
    pub kit_id: Option<String>,
    pub name: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SeqspecLibData {
    pub assay_id: String,
    pub assay_name: String,
    pub modality: String,
    pub library_region_id: String,
    pub seqspec_version: Option<String>,
    pub sequence_protocols: Vec<SeqspecMetadataRow>,
    pub sequence_kits: Vec<SeqspecMetadataRow>,
    pub library_protocols: Vec<SeqspecMetadataRow>,
    pub library_kits: Vec<SeqspecMetadataRow>,
    pub region_nodes: Vec<SeqspecLibRegion>,
    pub regions: Vec<SeqspecLibRegion>,
    pub reads: Vec<SeqspecLibRead>,
    pub total_bp: i64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RegionSpan {
    start: i64,
    end: i64,
}

pub fn load_report_json(path: &Path) -> Result<Report> {
    let file = File::open(path).with_context(|| format!("failed to read {}", path.display()))?;
    let report: Report = serde_json::from_reader(file)
        .with_context(|| format!("failed to parse JSON report from {}", path.display()))?;
    report.validate().context("report JSON failed validation")?;
    Ok(report)
}

pub fn render_report_html(
    report: &Report,
    lib_data: Option<&SeqspecLibData>,
    generated_at: &str,
    report_input_path: &Path,
    html_output_path: &Path,
) -> Result<String> {
    report.validate()?;

    let report_json = escape_script_json(&serde_json::to_string_pretty(report)?);
    let lib_json = escape_script_json(&serde_json::to_string_pretty(&lib_data)?);
    let generated_at_json = escape_script_json(&serde_json::to_string(generated_at)?);
    let repository_json = escape_script_json(&serde_json::to_string(REPOSITORY_URL)?);
    let version_json = escape_script_json(&serde_json::to_string(SEQCHECK_VERSION)?);
    let report_input_json = escape_script_json(&serde_json::to_string(
        &report_input_path.display().to_string(),
    )?);
    let report_output_json = escape_script_json(&serde_json::to_string(
        &html_output_path.display().to_string(),
    )?);
    let app_js = APP_JS.replace("</script", "<\\/script");

    Ok(TEMPLATE_HTML
        .replace("__SEQCHECK_STYLE__", STYLE_CSS)
        .replace("__SEQCHECK_REPORT_JSON__", &report_json)
        .replace("__SEQSPEC_LIB_JSON__", &lib_json)
        .replace("__SEQCHECK_GENERATED_AT__", &generated_at_json)
        .replace("__SEQCHECK_REPOSITORY__", &repository_json)
        .replace("__SEQCHECK_VERSION__", &version_json)
        .replace("__SEQCHECK_REPORT_INPUT__", &report_input_json)
        .replace("__SEQCHECK_REPORT_OUTPUT__", &report_output_json)
        .replace("__SEQCHECK_APP_JS__", &app_js))
}

pub fn load_optional_seqspec_lib_data(
    report_path: &Path,
    report: &Report,
    spec_override: Option<&Path>,
) -> Option<SeqspecLibData> {
    let spec_path = resolve_spec_path(report_path, report, spec_override)?;
    let assay = seqspec::utils::load_spec_path(&spec_path).ok()?;
    build_seqspec_lib_data(&assay, &report.meta.modality).ok()
}

pub fn build_seqspec_lib_data(spec: &Assay, modality: &str) -> Result<SeqspecLibData> {
    let libspec = spec
        .get_libspec(modality)
        .with_context(|| format!("modality '{}' is not present in the library_spec", modality))?;

    let mut region_nodes = Vec::new();
    let mut regions = Vec::new();
    let mut spans = HashMap::new();
    let mut total_bp = 0;
    for child in &libspec.regions {
        total_bp = collect_region_layout(
            child,
            0,
            total_bp,
            None,
            Vec::new(),
            Vec::new(),
            &mut region_nodes,
            &mut regions,
            &mut spans,
        );
    }

    let reads = spec
        .get_seqspec(modality)
        .into_iter()
        .filter_map(|read| {
            let primer_span = spans.get(&read.primer_id).copied()?;
            let read_len = read.max_len.max(0);
            let (start, end) = if read.strand == "neg" {
                (primer_span.start - read_len, primer_span.start)
            } else {
                (primer_span.end, primer_span.end + read_len)
            };

            Some(SeqspecLibRead {
                read_id: read.read_id,
                name: read.name.clone(),
                label: read.name,
                min_len: read.min_len,
                max_len: read.max_len,
                strand: read.strand,
                start,
                end,
                primer_id: read.primer_id,
                files: read.files.iter().map(file_view).collect(),
            })
        })
        .collect();

    Ok(SeqspecLibData {
        assay_id: spec.assay_id.clone(),
        assay_name: spec.name.clone(),
        modality: modality.to_string(),
        library_region_id: libspec.region_id,
        seqspec_version: spec.seqspec_version.clone(),
        sequence_protocols: seq_protocol_rows(spec.sequence_protocol.as_ref(), modality),
        sequence_kits: seq_kit_rows(spec.sequence_kit.as_ref(), modality),
        library_protocols: lib_protocol_rows(spec.library_protocol.as_ref(), modality),
        library_kits: lib_kit_rows(spec.library_kit.as_ref(), modality),
        region_nodes,
        regions,
        reads,
        total_bp,
    })
}

pub fn severity_rank(assessment_type: &AssessmentType) -> usize {
    match assessment_type {
        AssessmentType::Error => 0,
        AssessmentType::Warning => 1,
        AssessmentType::Interpretation => 2,
        AssessmentType::Pass => 3,
    }
}

pub fn severity_label(assessment_type: &AssessmentType) -> &'static str {
    match assessment_type {
        AssessmentType::Interpretation => "info",
        AssessmentType::Error => "error",
        AssessmentType::Warning => "warning",
        AssessmentType::Pass => "pass",
    }
}

pub fn generated_timestamp_utc() -> String {
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    format_unix_timestamp_utc(seconds)
}

fn resolve_spec_path(
    report_path: &Path,
    report: &Report,
    spec_override: Option<&Path>,
) -> Option<PathBuf> {
    if let Some(path) = spec_override {
        return Some(path.to_path_buf());
    }

    let spec_source = report.meta.spec.clone();
    if spec_source.is_empty() {
        return None;
    }
    if seqspec::utils::is_remote_source(&spec_source) {
        return None;
    }

    let spec_path = PathBuf::from(spec_source);
    if spec_path.is_absolute() || spec_path.exists() {
        return Some(spec_path);
    }

    report_path.parent().map(|parent| parent.join(spec_path))
}

#[allow(clippy::too_many_arguments)]
fn collect_region_layout(
    region: &Region,
    depth: usize,
    bp_start: i64,
    parent_region_id: Option<String>,
    path_region_ids: Vec<String>,
    path_names: Vec<String>,
    region_nodes: &mut Vec<SeqspecLibRegion>,
    leaves: &mut Vec<SeqspecLibRegion>,
    spans: &mut HashMap<String, RegionSpan>,
) -> i64 {
    let mut region_path_ids = path_region_ids;
    region_path_ids.push(region.region_id.clone());
    let mut region_path_names = path_names;
    region_path_names.push(region.name.clone());
    let end = if region.regions.is_empty() {
        let len = region.max_len.max(0);
        let node = SeqspecLibRegion {
            region_id: region.region_id.clone(),
            region_type: region.region_type.display(),
            name: region.name.clone(),
            sequence_type: region.sequence_type.clone(),
            sequence: if region.sequence.is_empty() {
                region.get_sequence()
            } else {
                region.sequence.clone()
            },
            min_len: region.min_len,
            max_len: region.max_len,
            len,
            bp_start,
            bp_end: bp_start + len,
            depth,
            parent_region_id,
            path_region_ids: region_path_ids,
            path_names: region_path_names,
            is_leaf: true,
            child_region_ids: Vec::new(),
            onlist: onlist_view(region.onlist.clone()),
        };
        region_nodes.push(node.clone());
        leaves.push(node);
        bp_start + len
    } else {
        let mut current = bp_start;
        for child in &region.regions {
            current = collect_region_layout(
                child,
                depth + 1,
                current,
                Some(region.region_id.clone()),
                region_path_ids.clone(),
                region_path_names.clone(),
                region_nodes,
                leaves,
                spans,
            );
        }
        region_nodes.push(SeqspecLibRegion {
            region_id: region.region_id.clone(),
            region_type: region.region_type.display(),
            name: region.name.clone(),
            sequence_type: region.sequence_type.clone(),
            sequence: if region.sequence.is_empty() {
                region.get_sequence()
            } else {
                region.sequence.clone()
            },
            min_len: region.min_len,
            max_len: region.max_len,
            len: current - bp_start,
            bp_start,
            bp_end: current,
            depth,
            parent_region_id,
            path_region_ids: region_path_ids,
            path_names: region_path_names,
            is_leaf: false,
            child_region_ids: region
                .regions
                .iter()
                .map(|child| child.region_id.clone())
                .collect(),
            onlist: onlist_view(region.onlist.clone()),
        });
        current
    };

    spans.insert(
        region.region_id.clone(),
        RegionSpan {
            start: bp_start,
            end,
        },
    );

    end
}

fn onlist_view(onlist: Option<SeqspecOnlist>) -> Option<SeqspecLibOnlist> {
    onlist.map(|onlist| SeqspecLibOnlist {
        file_id: onlist.file_id,
        filename: onlist.filename,
        filetype: onlist.filetype,
        filesize: onlist.filesize,
        url: onlist.url,
        urltype: onlist.urltype,
        md5: onlist.md5,
    })
}

fn file_view(file: &SeqspecFile) -> SeqspecLibFile {
    SeqspecLibFile {
        file_id: file.file_id.clone(),
        filename: file.filename.clone(),
        filetype: file.filetype.clone(),
        filesize: file.filesize,
        url: file.url.clone(),
        urltype: file.urltype.clone(),
        md5: file.md5.clone(),
    }
}

fn seq_protocol_rows(
    entries: Option<&Vec<SeqProtocol>>,
    modality: &str,
) -> Vec<SeqspecMetadataRow> {
    entries
        .map(|entries| {
            entries
                .iter()
                .filter(|entry| entry.modality == modality)
                .map(|entry| SeqspecMetadataRow {
                    protocol_id: Some(entry.protocol_id.clone()),
                    kit_id: None,
                    name: Some(entry.name.clone()),
                })
                .collect()
        })
        .unwrap_or_default()
}

fn seq_kit_rows(entries: Option<&Vec<SeqKit>>, modality: &str) -> Vec<SeqspecMetadataRow> {
    entries
        .map(|entries| {
            entries
                .iter()
                .filter(|entry| entry.modality == modality)
                .map(|entry| SeqspecMetadataRow {
                    protocol_id: None,
                    kit_id: Some(entry.kit_id.clone()),
                    name: entry.name.clone(),
                })
                .collect()
        })
        .unwrap_or_default()
}

fn lib_protocol_rows(
    entries: Option<&Vec<LibProtocol>>,
    modality: &str,
) -> Vec<SeqspecMetadataRow> {
    entries
        .map(|entries| {
            entries
                .iter()
                .filter(|entry| entry.modality == modality)
                .map(|entry| SeqspecMetadataRow {
                    protocol_id: Some(entry.protocol_id.clone()),
                    kit_id: None,
                    name: Some(entry.name.clone()),
                })
                .collect()
        })
        .unwrap_or_default()
}

fn lib_kit_rows(entries: Option<&Vec<LibKit>>, modality: &str) -> Vec<SeqspecMetadataRow> {
    entries
        .map(|entries| {
            entries
                .iter()
                .filter(|entry| entry.modality == modality)
                .map(|entry| SeqspecMetadataRow {
                    protocol_id: None,
                    kit_id: Some(entry.kit_id.clone()),
                    name: entry.name.clone(),
                })
                .collect()
        })
        .unwrap_or_default()
}

fn escape_script_json(json: &str) -> String {
    json.replace("</", "<\\/")
}

fn format_unix_timestamp_utc(seconds_since_epoch: u64) -> String {
    let days = (seconds_since_epoch / 86_400) as i64;
    let seconds_of_day = (seconds_since_epoch % 86_400) as u32;
    let hour = seconds_of_day / 3_600;
    let minute = (seconds_of_day % 3_600) / 60;
    let second = seconds_of_day % 60;
    let (year, month, day) = civil_from_days(days);

    format!("{year:04}-{month:02}-{day:02} {hour:02}:{minute:02}:{second:02} UTC")
}

fn civil_from_days(days_since_epoch: i64) -> (i64, u32, u32) {
    let z = days_since_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let day_of_era = z - era * 146_097;
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let year = year_of_era + era * 400;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_prime = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * month_prime + 2) / 5 + 1;
    let month = month_prime + if month_prime < 10 { 3 } else { -9 };
    let year = year + if month <= 2 { 1 } else { 0 };

    (year, month as u32, day as u32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::{
        Assessment, AtomicResult, MetricData, MetricDataKind, MetricItem, ReportMeta,
    };
    use seqspec::read::Read;
    use seqspec::region::Region;
    use serde_json::json;
    use std::fs;

    fn sample_report() -> Report {
        Report {
            report_schema_version: "0.1.0".to_string(),
            meta: ReportMeta {
                command: "check".to_string(),
                spec: "tests/fixtures/synthetic/spec.yaml".to_string(),
                modality: "rna".to_string(),
                requested_reads: 100,
            },
            results: vec![
                AtomicResult {
                    check: "onlist".to_string(),
                    files: vec!["synthetic_R1.fastq".to_string()],
                    reads: vec!["synthetic_R1".to_string()],
                    regions: vec!["barcode".to_string()],
                    expected: vec![MetricItem {
                        id: "e1".to_string(),
                        name: "onlist_source".to_string(),
                        description: "Whitelist source".to_string(),
                        data: MetricData {
                            kind: MetricDataKind::Scalar,
                            value: json!("tests/fixtures/synthetic/onlists/synthetic_barcodes.txt"),
                            unit: None,
                        },
                    }],
                    observed: vec![MetricItem {
                        id: "o1".to_string(),
                        name: "fetch_load_status".to_string(),
                        description: "Load status".to_string(),
                        data: MetricData {
                            kind: MetricDataKind::Scalar,
                            value: json!("ok"),
                            unit: None,
                        },
                    }],
                    assessment: vec![Assessment {
                        assessment_type: AssessmentType::Error,
                        code: "missing_onlist_resource".to_string(),
                        description: "The whitelist source could not be loaded.".to_string(),
                        expected_ids: vec!["e1".to_string()],
                        observed_ids: vec!["o1".to_string()],
                    }],
                },
                AtomicResult {
                    check: "length".to_string(),
                    files: vec!["synthetic_R1.fastq".to_string()],
                    reads: vec!["synthetic_R1".to_string()],
                    regions: Vec::new(),
                    expected: vec![MetricItem {
                        id: "e1".to_string(),
                        name: "expected_max_len".to_string(),
                        description: "Maximum length".to_string(),
                        data: MetricData {
                            kind: MetricDataKind::Scalar,
                            value: json!(12),
                            unit: Some("bp".to_string()),
                        },
                    }],
                    observed: vec![MetricItem {
                        id: "o1".to_string(),
                        name: "out_of_range_fraction".to_string(),
                        description: "Fraction out of range".to_string(),
                        data: MetricData {
                            kind: MetricDataKind::Scalar,
                            value: json!(0.2),
                            unit: Some("fraction".to_string()),
                        },
                    }],
                    assessment: vec![Assessment {
                        assessment_type: AssessmentType::Warning,
                        code: "length_out_of_range".to_string(),
                        description: "Observed read lengths fall outside the seqspec range."
                            .to_string(),
                        expected_ids: vec!["e1".to_string()],
                        observed_ids: vec!["o1".to_string()],
                    }],
                },
                AtomicResult {
                    check: "random".to_string(),
                    files: vec!["synthetic_R1.fastq".to_string()],
                    reads: vec!["synthetic_R1".to_string()],
                    regions: vec!["umi".to_string()],
                    expected: vec![MetricItem {
                        id: "e1".to_string(),
                        name: "max_entropy_bits".to_string(),
                        description: "Max entropy".to_string(),
                        data: MetricData {
                            kind: MetricDataKind::Scalar,
                            value: json!(4.0),
                            unit: Some("bits".to_string()),
                        },
                    }],
                    observed: vec![MetricItem {
                        id: "o1".to_string(),
                        name: "top_sequences".to_string(),
                        description: "Top sequences".to_string(),
                        data: MetricData {
                            kind: MetricDataKind::Records,
                            value: json!([{ "sequence": "AA", "count": 2 }]),
                            unit: None,
                        },
                    }],
                    assessment: vec![Assessment {
                        assessment_type: AssessmentType::Interpretation,
                        code: "random_sequence_entropy".to_string(),
                        description: "Observed sequence entropy is 0.5 of the theoretical maximum."
                            .to_string(),
                        expected_ids: vec!["e1".to_string()],
                        observed_ids: vec!["o1".to_string()],
                    }],
                },
                AtomicResult {
                    check: "coverage".to_string(),
                    files: vec!["synthetic_R1.fastq".to_string()],
                    reads: vec!["synthetic_R1".to_string()],
                    regions: Vec::new(),
                    expected: vec![MetricItem {
                        id: "e1".to_string(),
                        name: "expected_regions".to_string(),
                        description: "Expected regions".to_string(),
                        data: MetricData {
                            kind: MetricDataKind::Records,
                            value: json!([{ "region_id": "barcode", "start": 0, "stop": 4 }]),
                            unit: None,
                        },
                    }],
                    observed: vec![MetricItem {
                        id: "o1".to_string(),
                        name: "full_read_coverage_fraction".to_string(),
                        description: "Coverage fraction".to_string(),
                        data: MetricData {
                            kind: MetricDataKind::Scalar,
                            value: json!(1.0),
                            unit: Some("fraction".to_string()),
                        },
                    }],
                    assessment: vec![Assessment {
                        assessment_type: AssessmentType::Pass,
                        code: "full_read_coverage_complete".to_string(),
                        description: "All sampled reads fully cover the projected read geometry."
                            .to_string(),
                        expected_ids: vec!["e1".to_string()],
                        observed_ids: vec!["o1".to_string()],
                    }],
                },
            ],
        }
    }

    fn nested_spec() -> Assay {
        let fixed_a = Region::new(
            "fixed_a".into(),
            "linker".into(),
            "fixed a".into(),
            "fixed".into(),
            "AAA".into(),
            3,
            3,
            None,
            vec![],
        );
        let fixed_t = Region::new(
            "fixed_t".into(),
            "linker".into(),
            "fixed t".into(),
            "fixed".into(),
            "T".into(),
            1,
            1,
            None,
            vec![],
        );
        let joined_block = Region::new(
            "joined_block".into(),
            "linker".into(),
            "joined block".into(),
            "joined".into(),
            "AAAT".into(),
            4,
            4,
            None,
            vec![fixed_a, fixed_t],
        );
        let umi = Region::new(
            "umi".into(),
            "umi".into(),
            "umi".into(),
            "random".into(),
            "XX".into(),
            2,
            2,
            None,
            vec![],
        );
        let libspec = Region::new(
            "rna".into(),
            "rna".into(),
            "rna".into(),
            "joined".into(),
            "AAATXX".into(),
            6,
            6,
            None,
            vec![joined_block, umi],
        );
        let read = Read::new(
            "rna_R1".into(),
            "Read 1".into(),
            "rna".into(),
            "joined_block".into(),
            2,
            2,
            "pos".into(),
            vec![],
        );

        Assay::new(
            "nested-assay".into(),
            "Nested Assay".into(),
            "".into(),
            "2026-03-24".into(),
            "nested regions".into(),
            vec!["rna".into()],
            "".into(),
            vec![read],
            vec![libspec],
            None,
            None,
            None,
            None,
            Some("0.4.0".into()),
        )
    }

    #[test]
    fn test_render_report_html_embeds_payloads() {
        let html = render_report_html(
            &sample_report(),
            None,
            "2026-03-24 01:02:03 UTC",
            Path::new("/tmp/report.json"),
            Path::new("/tmp/report.html"),
        )
        .unwrap();

        assert!(html.contains("id=\"seqcheck-report-data\""));
        assert!(html.contains("\"report_schema_version\": \"0.1.0\""));
        assert!(html.contains("window.SEQCHECK_GENERATED_AT = \"2026-03-24 01:02:03 UTC\";"));
        assert!(html
            .contains("window.SEQCHECK_REPOSITORY = \"https://github.com/sbooeshaghi/seqcheck\";"));
        assert!(html.contains("window.SEQCHECK_VERSION = \"0.2.0\";"));
        assert!(html.contains("window.SEQCHECK_REPORT_INPUT = \"/tmp/report.json\";"));
        assert!(html.contains("window.SEQCHECK_REPORT_OUTPUT = \"/tmp/report.html\";"));
        assert!(html.contains("id=\"seqspec-lib-data\""));
        assert!(html.contains(">null</script>") || html.contains(">\nnull\n</script>"));
        assert!(html.contains("meta-section-head\">Run<"));
        assert!(!html.contains("meta-section-head\">Inputs<"));
    }

    #[test]
    fn test_load_report_json_rejects_invalid_assessment_references() {
        let mut invalid = sample_report();
        invalid.results[0].assessment[0].observed_ids = vec!["missing".to_string()];

        let path = temp_path("bad-report.json");
        fs::write(&path, serde_json::to_string_pretty(&invalid).unwrap()).unwrap();

        let error = load_report_json(&path).unwrap_err().to_string();
        assert!(error.contains("report JSON failed validation"));

        let _ = fs::remove_file(path);
    }

    #[test]
    fn test_build_seqspec_lib_data_orders_regions_and_maps_reads() {
        let spec =
            seqspec::utils::load_spec_path(Path::new("tests/fixtures/bad_geometry/spec.yaml"))
                .unwrap();
        let lib_data = build_seqspec_lib_data(&spec, "rna").unwrap();

        assert_eq!(lib_data.assay_id, "bad-geometry");
        assert_eq!(lib_data.library_region_id, "rna");
        assert_eq!(lib_data.total_bp, 22);
        assert_eq!(lib_data.sequence_protocols.len(), 0);
        assert_eq!(lib_data.regions[0].region_id, "primer1");
        assert_eq!(lib_data.regions[1].region_id, "barcode");
        assert_eq!(lib_data.regions[5].region_id, "primer2");
        assert_eq!(lib_data.regions[5].bp_start, 18);
        assert_eq!(lib_data.regions[1].path_region_ids, vec!["barcode"]);
        assert_eq!(lib_data.regions[1].path_names, vec!["Barcode"]);
        assert_eq!(
            lib_data.regions[1].onlist.as_ref().unwrap().url,
            "onlists/bad_barcodes.txt"
        );

        assert_eq!(lib_data.reads.len(), 2);
        assert_eq!(lib_data.reads[0].read_id, "bad_R1");
        assert_eq!(lib_data.reads[0].name, "Bad Read 1");
        assert_eq!(lib_data.reads[0].min_len, 6);
        assert_eq!(lib_data.reads[0].start, 4);
        assert_eq!(lib_data.reads[0].end, 10);
        assert_eq!(lib_data.reads[0].files.len(), 1);
        assert_eq!(lib_data.reads[0].files[0].file_id, "bad_R1.fastq");
        assert_eq!(lib_data.reads[1].read_id, "bad_R2");
        assert_eq!(lib_data.reads[1].start, 4);
        assert_eq!(lib_data.reads[1].end, 18);
    }

    #[test]
    fn test_build_seqspec_lib_data_omits_unresolved_primer_reads() {
        let mut spec =
            seqspec::utils::load_spec_path(Path::new("tests/fixtures/bad_geometry/spec.yaml"))
                .unwrap();
        spec.sequence_spec[0].primer_id = "missing_primer".to_string();

        let lib_data = build_seqspec_lib_data(&spec, "rna").unwrap();
        assert_eq!(lib_data.reads.len(), 1);
        assert_eq!(lib_data.reads[0].read_id, "bad_R2");
    }

    #[test]
    fn test_build_seqspec_lib_data_keeps_nested_regions() {
        let spec = nested_spec();
        let lib_data = build_seqspec_lib_data(&spec, "rna").unwrap();

        assert_eq!(lib_data.total_bp, 6);
        assert_eq!(lib_data.regions.len(), 3);
        assert_eq!(lib_data.regions[0].region_id, "fixed_a");
        assert_eq!(lib_data.regions[1].region_id, "fixed_t");
        assert_eq!(lib_data.regions[2].region_id, "umi");

        let parent = lib_data
            .region_nodes
            .iter()
            .find(|region| region.region_id == "joined_block")
            .unwrap();
        assert!(!parent.is_leaf);
        assert_eq!(parent.bp_start, 0);
        assert_eq!(parent.bp_end, 4);
        assert_eq!(parent.child_region_ids, vec!["fixed_a", "fixed_t"]);
        assert_eq!(parent.path_region_ids, vec!["joined_block"]);

        assert_eq!(lib_data.reads.len(), 1);
        assert_eq!(lib_data.reads[0].primer_id, "joined_block");
        assert_eq!(lib_data.reads[0].start, 4);
        assert_eq!(lib_data.reads[0].end, 6);
    }

    #[test]
    fn test_severity_helpers_match_ui_contract() {
        assert_eq!(severity_rank(&AssessmentType::Error), 0);
        assert_eq!(severity_rank(&AssessmentType::Warning), 1);
        assert_eq!(severity_rank(&AssessmentType::Interpretation), 2);
        assert_eq!(severity_rank(&AssessmentType::Pass), 3);

        assert_eq!(severity_label(&AssessmentType::Error), "error");
        assert_eq!(severity_label(&AssessmentType::Warning), "warning");
        assert_eq!(severity_label(&AssessmentType::Interpretation), "info");
        assert_eq!(severity_label(&AssessmentType::Pass), "pass");
    }

    #[test]
    fn test_generated_timestamp_utc_format() {
        let rendered = format_unix_timestamp_utc(0);
        assert_eq!(rendered, "1970-01-01 00:00:00 UTC");
    }

    fn temp_path(name: &str) -> PathBuf {
        let mut path = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        path.push(format!("seqcheck-{name}-{nanos}"));
        path
    }
}
