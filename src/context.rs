use crate::auth::RemoteAccess;
use anyhow::{anyhow, bail, Context, Result};
use flate2::read::GzDecoder;
use seqspec::assay::Assay;
use seqspec::file::File;
use seqspec::models::region_type::RegionTypeValue;
use seqspec::onlist::Onlist;
use seqspec::read::Read;
use seqspec::region::{Region, RegionCoordinate};
use serde::Serialize;
use std::collections::{BTreeMap, HashSet};
use std::io::{BufRead, BufReader, Read as IoRead};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Serialize)]
pub struct ExpectedRegion {
    pub region_id: String,
    pub name: String,
    pub region_type: RegionTypeValue,
    pub sequence_type: String,
    pub start: usize,
    pub stop: usize,
    pub expected_sequence: String,
}

#[derive(Debug, Clone)]
pub struct LoadedInputs {
    pub spec_source: String,
    pub spec: Assay,
    pub input_check: InputCheck,
    pub inputs: Vec<ResolvedInput>,
    pub remote_access: RemoteAccess,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct ExpectedFile {
    pub read_id: String,
    pub file_id: String,
    pub filename: String,
    pub url_basename: String,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct MatchedInput {
    pub input_path: String,
    pub read_id: String,
    pub file_id: String,
    pub matched_by: String,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct InputCheck {
    pub expected_files: Vec<ExpectedFile>,
    pub supplied_inputs: Vec<String>,
    pub matched_inputs: Vec<MatchedInput>,
    pub missing_expected_files: Vec<ExpectedFile>,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerClassificationKind {
    FixedScannable,
    GhostPrimer,
    NonScannablePrimer,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct PrimerClassification {
    pub kind: PrimerClassificationKind,
    pub scannable: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

#[derive(Debug, Clone)]
pub struct ResolvedInput {
    pub input_source: String,
    pub input_path: PathBuf,
    pub read: Read,
    pub matched_file: Option<File>,
    pub matched_by: String,
    pub coordinates: Vec<RegionCoordinate>,
    pub primer_region: Region,
    pub primer_classification: PrimerClassification,
    pub spec_base: Option<PathBuf>,
    pub remote_access: RemoteAccess,
}

#[derive(Debug, Clone)]
pub struct LoadedOnlist {
    pub source: String,
    pub entries: HashSet<String>,
}

#[derive(Debug, Clone)]
struct Candidate {
    rank: usize,
    read: Read,
    matched_file: Option<File>,
    matched_by: &'static str,
}

impl ResolvedInput {
    pub fn file_id(&self) -> String {
        self.matched_file
            .as_ref()
            .map(|file| file.file_id.clone())
            .unwrap_or_else(|| self.read.read_id.clone())
    }

    pub fn expected_regions(&self) -> Vec<ExpectedRegion> {
        self.coordinates
            .iter()
            .map(|region| ExpectedRegion {
                region_id: region.region.region_id.clone(),
                name: region.region.name.clone(),
                region_type: region.region.region_type.clone(),
                sequence_type: region.region.sequence_type.clone(),
                start: usize::try_from(region.start).unwrap_or_default(),
                stop: usize::try_from(region.stop).unwrap_or_default(),
                expected_sequence: region.region.sequence.clone(),
            })
            .collect()
    }

    pub fn expected_stop(&self) -> usize {
        self.coordinates
            .iter()
            .map(|region| usize::try_from(region.stop).unwrap_or_default())
            .max()
            .unwrap_or_default()
    }

    pub fn match_key(&self) -> (String, Option<String>) {
        (
            self.read.read_id.clone(),
            self.matched_file.as_ref().map(|file| file.file_id.clone()),
        )
    }
}

pub fn load_spec(spec_source: &str, remote_access: &RemoteAccess) -> Result<Assay> {
    Ok(normalize_spec_version(load_spec_raw(
        spec_source,
        remote_access,
    )?))
}

pub fn load_spec_raw(spec_source: &str, remote_access: &RemoteAccess) -> Result<Assay> {
    let spec = if seqspec::utils::is_remote_source(spec_source) {
        remote_access.with_reader(spec_source, |mut reader| {
            let mut data = Vec::new();
            reader.read_to_end(&mut data)?;
            seqspec::utils::load_spec_bytes(&data)
        })?
    } else {
        let spec_path = Path::new(spec_source);
        if !spec_path.exists() {
            bail!("spec file does not exist: {}", spec_source);
        }
        seqspec::utils::load_spec_path(spec_path)?
    };
    Ok(spec)
}

fn normalize_spec_version(spec: Assay) -> Assay {
    let version = spec
        .seqspec_version
        .clone()
        .unwrap_or_else(|| "0.0.0".to_string());

    match version.as_str() {
        "0.0.0" | "0.1.0" | "0.1.1" | "0.2.0" | "0.3.0" | "0.4.0" => {
            seqspec::seqspec_upgrade::seqspec_upgrade(spec, &version)
        }
        _ => spec,
    }
}

pub fn load_resolved_inputs(
    spec_source: &str,
    modality: &str,
    fastqs: &[String],
    _n_reads: usize,
    auth_profile: Option<&str>,
) -> Result<LoadedInputs> {
    let remote_access = RemoteAccess::load(auth_profile)?;
    let spec = load_spec(spec_source, &remote_access)?;
    let spec_base = seqspec::utils::spec_base_from_source(spec_source);
    let expected_files = expected_files_for_modality(&spec, modality)?;

    let inputs = fastqs
        .iter()
        .map(|fastq| resolve_input(&spec, modality, spec_base.as_deref(), fastq, &remote_access))
        .collect::<Result<Vec<_>>>()?;

    ensure_unique_matches(&inputs)?;

    let input_check = build_input_check(expected_files, fastqs, &inputs);

    Ok(LoadedInputs {
        spec_source: spec_source.to_string(),
        spec,
        input_check,
        inputs,
        remote_access,
    })
}

pub fn filter_regions_by_id(
    input: &ResolvedInput,
    region_id: &str,
) -> Result<Vec<RegionCoordinate>> {
    let filtered: Vec<RegionCoordinate> = input
        .coordinates
        .iter()
        .filter(|region| region.region.region_id == region_id)
        .cloned()
        .collect();

    if filtered.is_empty() {
        bail!(
            "region '{}' is not present on read '{}' for {}",
            region_id,
            input.read.read_id,
            input.input_source
        );
    }

    Ok(filtered)
}

pub fn filter_regions_by_sequence_type(
    input: &ResolvedInput,
    sequence_type: &str,
) -> Vec<RegionCoordinate> {
    input
        .coordinates
        .iter()
        .filter(|region| region.region.sequence_type == sequence_type)
        .cloned()
        .collect()
}

pub fn load_onlist(
    spec_base: Option<&Path>,
    region: &Region,
    remote_access: &RemoteAccess,
) -> Result<LoadedOnlist> {
    let onlist = region
        .onlist
        .as_ref()
        .with_context(|| format!("region '{}' has no onlist", region.region_id))?;
    let source = onlist_source(spec_base, onlist)?;

    let entries = match onlist.urltype.as_str() {
        "local" => {
            let lines = seqspec::utils::read_local_list(Path::new(&source))
                .map_err(|err| anyhow!("failed to read onlist '{}': {}", source, err))?;
            normalize_onlist_text(region, &lines.join("\n"))?
        }
        "http" | "https" | "ftp" => read_remote_onlist_entries(remote_access, region, &source)
            .with_context(|| format!("failed to stream remote onlist '{}'", source))?,
        other => bail!(
            "unsupported onlist urltype '{}' for region '{}'",
            other,
            region.region_id
        ),
    };

    Ok(LoadedOnlist { source, entries })
}

fn resolve_input(
    spec: &Assay,
    modality: &str,
    spec_base: Option<&Path>,
    input_source: &str,
    remote_access: &RemoteAccess,
) -> Result<ResolvedInput> {
    let input_path = if seqspec::utils::is_remote_source(input_source) {
        PathBuf::from(input_source_basename(input_source)?)
    } else {
        let path = PathBuf::from(input_source);
        if !path.exists() {
            bail!("FASTQ does not exist: {}", path.display());
        }
        path
    };

    let basename = input_source_basename(input_source)?;
    let candidate = resolve_candidate(spec, modality, &basename)?;
    let primer_region = resolve_primer_region(spec, modality, &candidate.read)?;
    let primer_classification = classify_primer_region(&primer_region);
    let (_, regions) =
        seqspec::utils::map_read_id_to_regions(spec, modality, &candidate.read.read_id)
            .map_err(|err| anyhow!("failed to map read '{}': {}", candidate.read.read_id, err))?;
    let coordinates = seqspec::utils::itx_read(
        seqspec::utils::project_regions_to_coordinates(regions),
        0,
        candidate.read.max_len,
    );

    Ok(ResolvedInput {
        input_source: input_source.to_string(),
        input_path,
        read: candidate.read,
        matched_file: candidate.matched_file,
        matched_by: candidate.matched_by.to_string(),
        coordinates,
        primer_region,
        primer_classification,
        spec_base: spec_base.map(Path::to_path_buf),
        remote_access: remote_access.clone(),
    })
}

fn input_source_basename(input_source: &str) -> Result<String> {
    if seqspec::utils::is_remote_source(input_source) {
        let trimmed = input_source
            .split_once('?')
            .map(|(head, _)| head)
            .unwrap_or(input_source)
            .split_once('#')
            .map(|(head, _)| head)
            .unwrap_or(input_source);
        let basename = trimmed.rsplit('/').next().unwrap_or_default().to_string();
        if basename.is_empty() {
            bail!("remote FASTQ URL has no basename: {}", input_source);
        }
        Ok(basename)
    } else {
        Path::new(input_source)
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| name.to_string())
            .ok_or_else(|| anyhow!("invalid FASTQ path: {}", input_source))
    }
}

fn resolve_candidate(spec: &Assay, modality: &str, basename: &str) -> Result<Candidate> {
    let reads = spec.get_seqspec(modality);
    if reads.is_empty() {
        bail!("modality '{}' is not present in the seqspec", modality);
    }

    let mut candidates = Vec::new();

    for read in reads {
        for file in &read.files {
            if file.file_id == basename {
                candidates.push(Candidate {
                    rank: 0,
                    read: read.clone(),
                    matched_file: Some(file.clone()),
                    matched_by: "file_id",
                });
            }
            if file.filename == basename {
                candidates.push(Candidate {
                    rank: 1,
                    read: read.clone(),
                    matched_file: Some(file.clone()),
                    matched_by: "filename",
                });
            }
            let url_basename = Path::new(&file.url)
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or_default();
            if !url_basename.is_empty() && url_basename == basename {
                candidates.push(Candidate {
                    rank: 2,
                    read: read.clone(),
                    matched_file: Some(file.clone()),
                    matched_by: "url_basename",
                });
            }
        }

        if read.read_id == basename {
            candidates.push(Candidate {
                rank: 3,
                read: read.clone(),
                matched_file: None,
                matched_by: "read_id",
            });
        }
    }

    candidates.sort_by(|left, right| {
        left.rank
            .cmp(&right.rank)
            .then_with(|| left.read.read_id.cmp(&right.read.read_id))
    });

    let Some(best) = candidates.first().cloned() else {
        bail!(
            "could not match '{}' to any read in modality '{}'; tried file_id, filename, url basename, then read_id",
            basename,
            modality
        );
    };

    if let Some(conflict) = candidates.iter().skip(1).find(|candidate| {
        candidate.rank == best.rank && candidate.read.read_id != best.read.read_id
    }) {
        bail!(
            "FASTQ '{}' matched multiple reads at the same priority: '{}' and '{}'",
            basename,
            best.read.read_id,
            conflict.read.read_id
        );
    }

    Ok(best)
}

fn onlist_source(spec_base: Option<&Path>, onlist: &Onlist) -> Result<String> {
    if onlist.urltype == "local" {
        let relative = PathBuf::from(
            seqspec::utils::local_onlist_locator(onlist).map_err(|err| anyhow!(err))?,
        );
        let resolved = if relative.is_absolute() {
            relative
        } else if let Some(base) = spec_base {
            base.join(relative)
        } else {
            bail!(
                "cannot resolve local onlist '{}' without a local seqspec source",
                onlist.filename
            );
        };
        Ok(resolved.to_string_lossy().to_string())
    } else {
        Ok(onlist.url.clone())
    }
}

fn expected_files_for_modality(spec: &Assay, modality: &str) -> Result<Vec<ExpectedFile>> {
    let reads = spec.get_seqspec(modality);
    if reads.is_empty() {
        bail!("modality '{}' is not present in the seqspec", modality);
    }

    Ok(reads
        .into_iter()
        .flat_map(|read| {
            read.files.into_iter().map(move |file| ExpectedFile {
                read_id: read.read_id.clone(),
                file_id: file.file_id,
                filename: file.filename.clone(),
                url_basename: Path::new(&file.url)
                    .file_name()
                    .and_then(|name| name.to_str())
                    .unwrap_or_default()
                    .to_string(),
            })
        })
        .collect())
}

fn build_input_check(
    expected_files: Vec<ExpectedFile>,
    fastqs: &[String],
    inputs: &[ResolvedInput],
) -> InputCheck {
    let matched_inputs = inputs
        .iter()
        .map(|input| MatchedInput {
            input_path: input.input_source.clone(),
            read_id: input.read.read_id.clone(),
            file_id: input.file_id(),
            matched_by: input.matched_by.clone(),
        })
        .collect::<Vec<_>>();

    let matched_keys = inputs
        .iter()
        .filter_map(|input| {
            input
                .matched_file
                .as_ref()
                .map(|file| (input.read.read_id.clone(), file.file_id.clone()))
        })
        .collect::<HashSet<_>>();

    let missing_expected_files = expected_files
        .iter()
        .filter(|expected| {
            !matched_keys.contains(&(expected.read_id.clone(), expected.file_id.clone()))
        })
        .cloned()
        .collect();

    InputCheck {
        expected_files,
        supplied_inputs: fastqs.to_vec(),
        matched_inputs,
        missing_expected_files,
    }
}

fn ensure_unique_matches(inputs: &[ResolvedInput]) -> Result<()> {
    let mut seen = BTreeMap::<(String, Option<String>), String>::new();

    for input in inputs {
        let key = input.match_key();
        if let Some(previous_path) = seen.insert(key.clone(), input.input_source.clone()) {
            let file_label = key
                .1
                .as_deref()
                .map(|file_id| format!("file '{}'", file_id))
                .unwrap_or_else(|| "read-only match".to_string());
            bail!(
                "FASTQ '{}' and '{}' both resolved to read '{}' ({})",
                previous_path,
                input.input_source,
                key.0,
                file_label
            );
        }
    }

    Ok(())
}

fn resolve_primer_region(spec: &Assay, modality: &str, read: &Read) -> Result<Region> {
    let libspec = spec
        .get_libspec(modality)
        .ok_or_else(|| anyhow!("modality '{}' is not present in the library_spec", modality))?;

    libspec
        .get_region_by_id(&read.primer_id)
        .into_iter()
        .next()
        .ok_or_else(|| {
            anyhow!(
                "read '{}' primer_id '{}' does not resolve to a region in modality '{}'",
                read.read_id,
                read.primer_id,
                modality
            )
        })
}

fn classify_primer_region(region: &Region) -> PrimerClassification {
    let primer_len = usize::try_from(region.max_len.max(0)).unwrap_or_default();
    if primer_len == 0 {
        return PrimerClassification {
            kind: PrimerClassificationKind::GhostPrimer,
            scannable: false,
            reason: Some("zero-length primer anchor".to_string()),
        };
    }

    if region.sequence_type != "fixed" {
        return PrimerClassification {
            kind: PrimerClassificationKind::NonScannablePrimer,
            scannable: false,
            reason: Some(format!(
                "sequence_type '{}' is not fixed",
                region.sequence_type
            )),
        };
    }

    if region.sequence.is_empty() {
        return PrimerClassification {
            kind: PrimerClassificationKind::NonScannablePrimer,
            scannable: false,
            reason: Some("fixed primer sequence is empty".to_string()),
        };
    }

    if !sequence_is_concrete_dna(&region.sequence) {
        return PrimerClassification {
            kind: PrimerClassificationKind::NonScannablePrimer,
            scannable: false,
            reason: Some("fixed primer sequence contains non-ACGT characters".to_string()),
        };
    }

    PrimerClassification {
        kind: PrimerClassificationKind::FixedScannable,
        scannable: true,
        reason: None,
    }
}

fn sequence_is_concrete_dna(sequence: &str) -> bool {
    sequence
        .chars()
        .all(|base| matches!(base.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T'))
}

fn normalize_onlist_reader<R>(reader: R) -> Result<HashSet<String>>
where
    R: IoRead,
{
    let mut entries = HashSet::new();
    for line in BufReader::new(reader).lines() {
        let token = line?
            .split('\t')
            .next()
            .unwrap_or_default()
            .trim()
            .to_string();
        if !token.is_empty() {
            entries.insert(token);
        }
    }
    Ok(entries)
}

fn normalize_onlist_text(region: &Region, text: &str) -> Result<HashSet<String>> {
    if let Some(entries) = try_parse_delimited_onlist(region, text)? {
        return Ok(entries);
    }
    normalize_onlist_reader(std::io::Cursor::new(text))
}

fn try_parse_delimited_onlist(region: &Region, text: &str) -> Result<Option<HashSet<String>>> {
    let candidate_lines = text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect::<Vec<_>>();
    if candidate_lines.is_empty() {
        return Ok(None);
    }

    let mut best_delimiter = None;
    let mut best_column_rows = 0usize;
    let mut best_width = 0usize;
    for delimiter in [',', '\t'] {
        let rows_with_columns = candidate_lines
            .iter()
            .filter(|line| line.contains(delimiter))
            .count();
        let width = candidate_lines
            .iter()
            .map(|line| line.split(delimiter).count())
            .max()
            .unwrap_or_default();
        if rows_with_columns > best_column_rows
            || (rows_with_columns == best_column_rows && width > best_width)
        {
            best_delimiter = Some(delimiter);
            best_column_rows = rows_with_columns;
            best_width = width;
        }
    }

    let Some(delimiter) = best_delimiter else {
        return Ok(None);
    };
    if best_column_rows == 0 || best_width < 2 {
        return Ok(None);
    }

    let rows = candidate_lines
        .iter()
        .map(|line| {
            line.split(delimiter)
                .map(normalize_cell)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let has_header = detect_header(region, &rows);
    let header = has_header.then(|| rows[0].clone());
    let data_rows = if has_header { &rows[1..] } else { &rows[..] };
    if data_rows.is_empty() {
        return Ok(None);
    }

    let mut column_match_counts = vec![0usize; best_width];
    for row in data_rows {
        for (idx, value) in row.iter().enumerate() {
            if field_matches_region(value, region) {
                column_match_counts[idx] += 1;
            }
        }
    }

    let candidate_columns = column_match_counts
        .iter()
        .enumerate()
        .filter_map(|(idx, count)| if *count > 0 { Some(idx) } else { None })
        .collect::<Vec<_>>();
    if candidate_columns.is_empty() {
        return Ok(None);
    }

    let selected_columns = select_onlist_columns(region, header.as_deref(), &column_match_counts);
    let selected_columns = if selected_columns.is_empty() {
        candidate_columns
    } else {
        selected_columns
    };

    let mut entries = HashSet::new();
    for row in data_rows {
        for idx in &selected_columns {
            if let Some(value) = row.get(*idx) {
                if field_matches_region(value, region) {
                    entries.insert(value.clone());
                }
            }
        }
    }

    if entries.is_empty() {
        return Ok(None);
    }

    Ok(Some(entries))
}

fn normalize_cell(value: &str) -> String {
    value.trim().trim_matches('"').to_string()
}

fn field_matches_region(value: &str, region: &Region) -> bool {
    if value.is_empty() || !sequence_is_concrete_dna(value) {
        return false;
    }
    let len = value.len();
    let min_len = usize::try_from(region.min_len).unwrap_or_default();
    let max_len = usize::try_from(region.max_len).unwrap_or_default();
    if max_len == 0 {
        return len > 0;
    }
    len >= min_len && len <= max_len
}

fn detect_header(region: &Region, rows: &[Vec<String>]) -> bool {
    if rows.len() < 2 {
        return false;
    }
    let first_matches = rows[0]
        .iter()
        .filter(|value| field_matches_region(value, region))
        .count();
    let subsequent_matches = rows[1..]
        .iter()
        .take(5)
        .map(|row| {
            row.iter()
                .filter(|value| field_matches_region(value, region))
                .count()
        })
        .sum::<usize>();

    first_matches == 0 && subsequent_matches > 0
}

fn select_onlist_columns(
    region: &Region,
    header: Option<&[String]>,
    column_match_counts: &[usize],
) -> Vec<usize> {
    let candidate_columns = column_match_counts
        .iter()
        .enumerate()
        .filter_map(|(idx, count)| if *count > 0 { Some(idx) } else { None })
        .collect::<Vec<_>>();
    if candidate_columns.is_empty() {
        return Vec::new();
    }

    if let Some(header) = header {
        let mut best_score = i32::MIN;
        let mut best_columns = Vec::new();
        for idx in &candidate_columns {
            let score =
                region_header_score(region, header.get(*idx).map(|s| s.as_str()).unwrap_or(""))
                    * 1000
                    + i32::try_from(column_match_counts[*idx]).unwrap_or_default();
            if score > best_score {
                best_score = score;
                best_columns.clear();
                best_columns.push(*idx);
            } else if score == best_score {
                best_columns.push(*idx);
            }
        }
        if best_score
            > i32::try_from(column_match_counts[*candidate_columns.first().unwrap()])
                .unwrap_or_default()
        {
            return best_columns;
        }
    }

    let max_matches = candidate_columns
        .iter()
        .map(|idx| column_match_counts[*idx])
        .max()
        .unwrap_or_default();
    candidate_columns
        .into_iter()
        .filter(|idx| column_match_counts[*idx] == max_matches)
        .collect()
}

fn region_header_score(region: &Region, header: &str) -> i32 {
    let header = header.to_ascii_lowercase();
    let mut score = 0;

    if region.region_type.is_index5() {
        if header.contains("i5") || header.contains("index5") || header.contains("index 5") {
            score += 8;
        }
        if header.contains("index2") {
            score += 6;
        }
        if header.contains("i7") {
            score -= 8;
        }
    }
    if region.region_type.is_index7() {
        if header.contains("i7") || header.contains("index7") || header.contains("index 7") {
            score += 8;
        }
        if header.contains("index(") {
            score += 2;
        }
        if header.contains("i5") || header.contains("index2") {
            score -= 8;
        }
    }
    if region.region_type.is_cell_barcode() {
        if header.contains("barcode") {
            score += 8;
        }
        if header.contains("cell") || header.contains("cb") {
            score += 4;
        }
    }
    if region.region_type.is_molecule_barcode() && header.contains("umi") {
        score += 8;
    }
    if region.region_type.has_term("RGN:measure:guide")
        && (header.contains("spacer")
            || header.contains("guide")
            || header.contains("grna")
            || header.contains("sgrna"))
    {
        score += 8;
    }

    for token in tokenize_region_metadata(region) {
        if token.len() > 1 && header.contains(&token) {
            score += 2;
        }
    }

    score
}

fn tokenize_region_metadata(region: &Region) -> Vec<String> {
    format!(
        "{} {} {}",
        region.region_id, region.name, region.region_type
    )
    .split(|c: char| !c.is_ascii_alphanumeric())
    .filter(|token| !token.is_empty())
    .map(|token| token.to_ascii_lowercase())
    .collect()
}

fn read_remote_onlist_entries(
    remote_access: &RemoteAccess,
    region: &Region,
    url: &str,
) -> Result<HashSet<String>> {
    remote_access.with_reader(url, |mut reader| {
        let mut data = Vec::new();
        reader.read_to_end(&mut data)?;
        let text = if url.ends_with(".gz") {
            let mut decoder = GzDecoder::new(&data[..]);
            let mut text = String::new();
            decoder.read_to_string(&mut text)?;
            text
        } else {
            String::from_utf8(data)?
        };
        normalize_onlist_text(region, &text)
    })
}

#[cfg(test)]
fn normalize_onlist_tokens(lines: Vec<String>) -> HashSet<String> {
    lines
        .into_iter()
        .filter_map(|line| {
            let token = line
                .split('\t')
                .next()
                .unwrap_or_default()
                .trim()
                .to_string();
            if token.is_empty() {
                None
            } else {
                Some(token)
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use seqspec::assay::Assay;
    use seqspec::file::File;
    use seqspec::onlist::Onlist;
    use seqspec::read::Read;
    use seqspec::region::Region;
    use std::io::{Read as IoReadTrait, Write};
    use std::net::TcpListener;
    use std::thread;

    fn sample_file(file_id: &str, filename: &str, url: &str) -> File {
        File::new(
            file_id.to_string(),
            filename.to_string(),
            "fastq".to_string(),
            0,
            url.to_string(),
            "local".to_string(),
            String::new(),
        )
    }

    fn sample_assay() -> Assay {
        let read = Read::new(
            "rna_R1".to_string(),
            "RNA Read 1".to_string(),
            "rna".to_string(),
            "primer".to_string(),
            8,
            8,
            "pos".to_string(),
            vec![sample_file(
                "R1.fastq.gz",
                "R1.fastq.gz",
                "fastqs/R1.fastq.gz",
            )],
        );
        let library = Region::new(
            "rna".to_string(),
            "rna".to_string(),
            "RNA".to_string(),
            "joined".to_string(),
            "AAAANNNNTT".to_string(),
            10,
            10,
            None,
            vec![
                Region::new(
                    "primer".to_string(),
                    "truseq_read1".to_string(),
                    "Primer".to_string(),
                    "fixed".to_string(),
                    "AAAA".to_string(),
                    4,
                    4,
                    None,
                    vec![],
                ),
                Region::new(
                    "barcode".to_string(),
                    "barcode".to_string(),
                    "Barcode".to_string(),
                    "onlist".to_string(),
                    "NNNN".to_string(),
                    4,
                    4,
                    Some(Onlist::new(
                        "barcodes.txt".to_string(),
                        "barcodes.txt".to_string(),
                        "txt".to_string(),
                        0,
                        "barcodes.txt".to_string(),
                        "local".to_string(),
                        String::new(),
                    )),
                    vec![],
                ),
                Region::new(
                    "tail".to_string(),
                    "linker".to_string(),
                    "Tail".to_string(),
                    "fixed".to_string(),
                    "TT".to_string(),
                    2,
                    2,
                    None,
                    vec![],
                ),
            ],
        );

        Assay::new(
            "test".to_string(),
            "test".to_string(),
            String::new(),
            "2026-03-19".to_string(),
            String::new(),
            vec!["rna".to_string()],
            String::new(),
            vec![read],
            vec![library],
            None,
            None,
            None,
            None,
            Some("0.4.0".to_string()),
        )
    }

    fn primer_region(assay: &Assay, region_id: &str) -> Region {
        assay
            .get_libspec("rna")
            .unwrap()
            .get_region_by_id(region_id)
            .into_iter()
            .next()
            .unwrap()
    }

    #[test]
    fn test_resolve_candidate_matches_by_file_id() {
        let candidate = resolve_candidate(&sample_assay(), "rna", "R1.fastq.gz").unwrap();
        assert_eq!(candidate.read.read_id, "rna_R1");
        assert_eq!(candidate.matched_by, "file_id");
    }

    #[test]
    fn test_resolve_candidate_matches_by_read_id() {
        let candidate = resolve_candidate(&sample_assay(), "rna", "rna_R1").unwrap();
        assert_eq!(candidate.read.read_id, "rna_R1");
        assert_eq!(candidate.matched_by, "read_id");
    }

    #[test]
    fn test_load_onlist_reads_first_column() {
        let root = std::env::temp_dir().join(format!(
            "seqcheck-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&root).unwrap();
        std::fs::write(root.join("barcodes.txt"), "AAAA\t1\nCCCC\t2\n").unwrap();

        let library = sample_assay().get_libspec("rna").unwrap();
        let barcode = library.get_region_by_id("barcode").pop().unwrap();
        let access = RemoteAccess::anonymous();
        let loaded = load_onlist(Some(&root), &barcode, &access).unwrap();

        assert!(loaded.entries.contains("AAAA"));
        assert!(loaded.entries.contains("CCCC"));
        assert_eq!(loaded.entries.len(), 2);

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn test_onlist_source_prefers_url_when_present_for_local_onlists() {
        let root = PathBuf::from("/tmp/spec-root");
        let onlist = Onlist::new(
            "ol".to_string(),
            "barcodes.txt".to_string(),
            "txt".to_string(),
            0,
            "nested/barcodes.txt".to_string(),
            "local".to_string(),
            String::new(),
        );

        let source = onlist_source(Some(&root), &onlist).unwrap();
        assert_eq!(source, "/tmp/spec-root/nested/barcodes.txt");
    }

    #[test]
    fn test_onlist_source_errors_when_local_url_is_empty() {
        let root = PathBuf::from("/tmp/spec-root");
        let onlist = Onlist::new(
            "ol".to_string(),
            "barcodes.txt".to_string(),
            "txt".to_string(),
            0,
            String::new(),
            "local".to_string(),
            String::new(),
        );

        let error = onlist_source(Some(&root), &onlist).unwrap_err();
        assert_eq!(
            error.to_string(),
            "local onlist 'barcodes.txt' has empty url"
        );
    }

    #[test]
    fn test_normalize_onlist_reader_matches_line_tokenizer() {
        let text = "AAAA\t1\nCCCC\t2\n\nGGGG\n";
        let observed = normalize_onlist_reader(text.as_bytes()).unwrap();
        let expected = normalize_onlist_tokens(vec![
            "AAAA\t1".to_string(),
            "CCCC\t2".to_string(),
            String::new(),
            "GGGG".to_string(),
        ]);

        assert_eq!(observed, expected);
    }

    #[test]
    fn test_load_spec_upgrades_0_3_0_to_current_version() {
        let mut assay = sample_assay();
        assay.seqspec_version = Some("0.3.0".to_string());

        let root = std::env::temp_dir().join(format!(
            "seqcheck-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&root).unwrap();
        let path = root.join("spec.yaml");
        std::fs::write(&path, assay.to_bytes().unwrap()).unwrap();

        let loaded = load_spec(path.to_str().unwrap(), &RemoteAccess::anonymous()).unwrap();
        assert_eq!(loaded.seqspec_version.as_deref(), Some("0.5.0"));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn test_normalize_spec_upgrades_0_4_0_region_types() {
        let loaded = normalize_spec_version(sample_assay());

        assert_eq!(loaded.seqspec_version.as_deref(), Some("0.5.0"));
        assert!(primer_region(&loaded, "barcode")
            .region_type
            .has_term("RGN:partition:cell"));
    }

    #[test]
    fn test_normalize_onlist_reader_handles_gzip_stream() {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"AAAA\t1\nCCCC\t2\n").unwrap();
        let compressed = encoder.finish().unwrap();

        let observed = normalize_onlist_reader(GzDecoder::new(&compressed[..])).unwrap();

        assert!(observed.contains("AAAA"));
        assert!(observed.contains("CCCC"));
        assert_eq!(observed.len(), 2);
    }

    #[test]
    fn test_normalize_onlist_text_parses_illumina_dual_index_csv_for_index7() {
        let region = Region::new_with_region_type_value(
            "idx7".to_string(),
            RegionTypeValue::from(vec![
                "RGN:partition:sample".to_string(),
                "RGN:technical:index7".to_string(),
            ]),
            "Index 7".to_string(),
            "onlist".to_string(),
            String::new(),
            10,
            10,
            Some(Onlist::new(
                "ol".to_string(),
                "illumina.csv.gz".to_string(),
                "csv.gz".to_string(),
                0,
                "illumina.csv.gz".to_string(),
                "local".to_string(),
                String::new(),
            )),
            vec![],
        );
        let text = "\
# comment\n\
index_name,index(i7),index2_workflow_a(i5),index2_workflow_b(i5)\n\
SI-A1,CCTGTCAGGG,AGTGTTACCT,AGGTAACACT\n\
SI-A2,GTGGATCAAA,GCCAACCCTG,CAGGGTTGGC\n";

        let observed = normalize_onlist_text(&region, text).unwrap();

        assert!(observed.contains("CCTGTCAGGG"));
        assert!(observed.contains("GTGGATCAAA"));
        assert_eq!(observed.len(), 2);
    }

    #[test]
    fn test_region_header_score_uses_ontology_semantics() {
        let cell_barcode = Region::new_with_region_type_value(
            "cell_id".to_string(),
            RegionTypeValue::from(vec!["RGN:partition:cell".to_string()]),
            "Cell identifier".to_string(),
            "onlist".to_string(),
            String::new(),
            16,
            16,
            None,
            vec![],
        );
        let guide = Region::new_with_region_type_value(
            "feature".to_string(),
            RegionTypeValue::from(vec![
                "RGN:measure:guide".to_string(),
                "RGN:classify:perturbation".to_string(),
            ]),
            "Perturbation feature".to_string(),
            "onlist".to_string(),
            String::new(),
            20,
            20,
            None,
            vec![],
        );

        assert!(region_header_score(&cell_barcode, "cell_barcode") >= 8);
        assert!(region_header_score(&guide, "sgRNA spacer") >= 8);
    }

    #[test]
    fn test_normalize_onlist_text_parses_illumina_dual_index_csv_for_index5() {
        let region = Region::new(
            "idx5".to_string(),
            "index5".to_string(),
            "Index 5".to_string(),
            "onlist".to_string(),
            String::new(),
            10,
            10,
            Some(Onlist::new(
                "ol".to_string(),
                "illumina.csv.gz".to_string(),
                "csv.gz".to_string(),
                0,
                "illumina.csv.gz".to_string(),
                "local".to_string(),
                String::new(),
            )),
            vec![],
        );
        let text = "\
# comment\n\
index_name,index(i7),index2_workflow_a(i5),index2_workflow_b(i5)\n\
SI-A1,CCTGTCAGGG,AGTGTTACCT,AGGTAACACT\n\
SI-A2,GTGGATCAAA,GCCAACCCTG,CAGGGTTGGC\n";

        let observed = normalize_onlist_text(&region, text).unwrap();

        assert!(observed.contains("AGTGTTACCT"));
        assert!(observed.contains("AGGTAACACT"));
        assert!(observed.contains("GCCAACCCTG"));
        assert!(observed.contains("CAGGGTTGGC"));
        assert_eq!(observed.len(), 4);
    }

    #[test]
    fn test_read_remote_onlist_entries_streams_gzip_over_http() {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"AAAA\t1\nCCCC\t2\n").unwrap();
        let body = encoder.finish().unwrap();

        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let addr = listener.local_addr().unwrap();
        let server = thread::spawn(move || {
            let (mut stream, _) = listener.accept().unwrap();
            let mut buffer = [0_u8; 4096];
            let _ = stream.read(&mut buffer);
            let response = format!(
                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            );
            stream.write_all(response.as_bytes()).unwrap();
            stream.write_all(&body).unwrap();
        });

        let access = RemoteAccess::anonymous();
        let region = Region::new(
            "barcode".to_string(),
            "barcode".to_string(),
            "Barcode".to_string(),
            "onlist".to_string(),
            String::new(),
            4,
            4,
            None,
            vec![],
        );
        let observed = read_remote_onlist_entries(
            &access,
            &region,
            &format!("http://{}/barcodes.txt.gz", addr),
        )
        .unwrap();

        server.join().unwrap();

        assert!(observed.contains("AAAA"));
        assert!(observed.contains("CCCC"));
        assert_eq!(observed.len(), 2);
    }

    #[test]
    fn test_load_resolved_inputs_accepts_remote_spec_and_fastq() {
        let _guard = crate::auth::test_env_lock()
            .lock()
            .unwrap_or_else(|err| err.into_inner());
        std::env::remove_var("SEQCHECK_AUTH_CONFIG");

        let spec_bytes = sample_assay().to_bytes().unwrap();
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"@r1\nAAAACCCC\n+\nFFFFFFFF\n").unwrap();
        let fastq_bytes = encoder.finish().unwrap();

        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let addr = listener.local_addr().unwrap();
        let server = thread::spawn(move || {
            for _ in 0..2 {
                let (mut stream, _) = listener.accept().unwrap();
                let mut buffer = [0_u8; 4096];
                let bytes_read = stream.read(&mut buffer).unwrap();
                let request = String::from_utf8_lossy(&buffer[..bytes_read]);
                let (path, body) = if request.starts_with("GET /spec.yaml") {
                    ("/spec.yaml", &spec_bytes)
                } else {
                    ("/R1.fastq.gz", &fastq_bytes)
                };
                assert!(request.contains(path));
                let response = format!(
                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                    body.len()
                );
                stream.write_all(response.as_bytes()).unwrap();
                stream.write_all(body).unwrap();
            }
        });

        let spec_url = format!("http://{}/spec.yaml", addr);
        let fastq_url = format!("http://{}/R1.fastq.gz", addr);
        let loaded =
            load_resolved_inputs(&spec_url, "rna", std::slice::from_ref(&fastq_url), 1, None)
                .unwrap();

        assert_eq!(loaded.spec_source, spec_url);
        assert_eq!(loaded.input_check.supplied_inputs, vec![fastq_url.clone()]);
        assert_eq!(loaded.inputs.len(), 1);
        assert_eq!(loaded.inputs[0].input_source, fastq_url);
        assert!(loaded.inputs[0].spec_base.is_none());
        let sampled =
            crate::scan::scan_fastq_records(&loaded.inputs[0], 1, |_, _, _| Ok(())).unwrap();
        assert_eq!(sampled, 1);

        server.join().unwrap();
    }

    #[test]
    fn test_classify_primer_region_fixed_scannable() {
        let classification = classify_primer_region(&primer_region(&sample_assay(), "primer"));
        assert_eq!(
            classification.kind,
            PrimerClassificationKind::FixedScannable
        );
        assert!(classification.scannable);
        assert!(classification.reason.is_none());
    }

    #[test]
    fn test_classify_primer_region_ghost_primer() {
        let region = Region::new(
            "ghost".to_string(),
            "truseq_read1".to_string(),
            "Ghost Primer".to_string(),
            "fixed".to_string(),
            String::new(),
            0,
            0,
            None,
            vec![],
        );

        let classification = classify_primer_region(&region);

        assert_eq!(classification.kind, PrimerClassificationKind::GhostPrimer);
        assert!(!classification.scannable);
        assert_eq!(
            classification.reason.as_deref(),
            Some("zero-length primer anchor")
        );
    }

    #[test]
    fn test_classify_primer_region_non_scannable() {
        let region = Region::new(
            "bad_primer".to_string(),
            "barcode".to_string(),
            "Bad Primer".to_string(),
            "random".to_string(),
            "XX".to_string(),
            2,
            2,
            None,
            vec![],
        );

        let classification = classify_primer_region(&region);

        assert_eq!(
            classification.kind,
            PrimerClassificationKind::NonScannablePrimer
        );
        assert!(!classification.scannable);
        assert_eq!(
            classification.reason.as_deref(),
            Some("sequence_type 'random' is not fixed")
        );
    }

    #[test]
    fn test_filter_regions_by_id_finds_expected_region() {
        let assay = sample_assay();
        let spec_base = PathBuf::from(".");
        let (_, regions) = seqspec::utils::map_read_id_to_regions(&assay, "rna", "rna_R1").unwrap();
        let coordinates = seqspec::utils::itx_read(
            seqspec::utils::project_regions_to_coordinates(regions),
            0,
            8,
        );
        let input = ResolvedInput {
            input_source: "R1.fastq.gz".to_string(),
            input_path: PathBuf::from("R1.fastq.gz"),
            read: assay.get_read("rna_R1").unwrap(),
            matched_file: None,
            matched_by: "read_id".to_string(),
            coordinates,
            primer_region: primer_region(&assay, "primer"),
            primer_classification: classify_primer_region(&primer_region(&assay, "primer")),
            spec_base: Some(spec_base),
            remote_access: RemoteAccess::anonymous(),
        };

        let barcode = filter_regions_by_id(&input, "barcode").unwrap();
        assert_eq!(barcode.len(), 1);
        assert_eq!(barcode[0].region.region_id, "barcode");
    }
}
