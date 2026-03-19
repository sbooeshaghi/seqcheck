use anyhow::{anyhow, bail, Context, Result};
use flate2::read::GzDecoder;
use seqspec::assay::Assay;
use seqspec::file::File;
use seqspec::onlist::Onlist;
use seqspec::read::Read;
use seqspec::region::{Region, RegionCoordinate};
use serde::Serialize;
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Read as IoRead};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

#[derive(Debug, Clone, Serialize)]
pub struct ExpectedRegion {
    pub region_id: String,
    pub name: String,
    pub region_type: String,
    pub sequence_type: String,
    pub start: usize,
    pub stop: usize,
    pub expected_sequence: String,
}

#[derive(Debug, Clone)]
pub struct ResolvedInput {
    pub input_path: PathBuf,
    pub read: Read,
    pub matched_file: Option<File>,
    pub matched_by: String,
    pub coordinates: Vec<RegionCoordinate>,
    pub spec_base: PathBuf,
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
}

pub fn load_spec(spec_path: &Path) -> Result<Assay> {
    if !spec_path.exists() {
        bail!("spec file does not exist: {}", spec_path.display());
    }

    Ok(seqspec::utils::load_spec(&spec_path.to_path_buf()))
}

pub fn load_resolved_inputs(
    spec_path: &Path,
    modality: &str,
    fastqs: &[PathBuf],
) -> Result<(Assay, Vec<ResolvedInput>)> {
    let spec = load_spec(spec_path)?;
    let spec_base = spec_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf();

    let inputs = fastqs
        .iter()
        .map(|fastq| resolve_input(&spec, modality, &spec_base, fastq))
        .collect::<Result<Vec<_>>>()?;

    Ok((spec, inputs))
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
            input.input_path.display()
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

pub fn load_onlist(spec_base: &Path, region: &Region) -> Result<LoadedOnlist> {
    let onlist = region
        .onlist
        .as_ref()
        .with_context(|| format!("region '{}' has no onlist", region.region_id))?;
    let source = onlist_source(spec_base, onlist);

    let entries = match onlist.urltype.as_str() {
        "local" => {
            let lines = seqspec::utils::read_local_list(Path::new(&source))
                .map_err(|err| anyhow!("failed to read onlist '{}': {}", source, err))?;
            normalize_onlist_lines(lines)
        }
        "http" | "https" | "ftp" => read_remote_onlist_entries(&source)
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
    spec_base: &Path,
    input_path: &Path,
) -> Result<ResolvedInput> {
    if !input_path.exists() {
        bail!("FASTQ does not exist: {}", input_path.display());
    }

    let basename = input_path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| anyhow!("invalid FASTQ path: {}", input_path.display()))?;

    let candidate = resolve_candidate(spec, modality, basename)?;
    let (_, regions) =
        seqspec::utils::map_read_id_to_regions(spec, modality, &candidate.read.read_id)
            .map_err(|err| anyhow!("failed to map read '{}': {}", candidate.read.read_id, err))?;
    let coordinates = seqspec::utils::itx_read(
        seqspec::utils::project_regions_to_coordinates(regions),
        0,
        candidate.read.max_len,
    );

    Ok(ResolvedInput {
        input_path: input_path.to_path_buf(),
        read: candidate.read,
        matched_file: candidate.matched_file,
        matched_by: candidate.matched_by.to_string(),
        coordinates,
        spec_base: spec_base.to_path_buf(),
    })
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

fn onlist_source(spec_base: &Path, onlist: &Onlist) -> String {
    if onlist.urltype == "local" {
        let relative = if onlist.url.is_empty() {
            PathBuf::from(&onlist.filename)
        } else {
            PathBuf::from(&onlist.url)
        };
        let resolved = if relative.is_absolute() {
            relative
        } else {
            spec_base.join(relative)
        };
        resolved.to_string_lossy().to_string()
    } else {
        onlist.url.clone()
    }
}

fn normalize_onlist_lines(lines: Vec<String>) -> HashSet<String> {
    normalize_onlist_reader(std::io::Cursor::new(lines.join("\n"))).unwrap_or_default()
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

fn read_remote_onlist_entries(url: &str) -> Result<HashSet<String>> {
    let mut child = Command::new("curl")
        .arg("-fsSL")
        .arg(url)
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .with_context(|| format!("failed to spawn curl for '{}'", url))?;

    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| anyhow!("curl did not provide stdout for '{}'", url))?;

    let entries = if url.ends_with(".gz") {
        normalize_onlist_reader(GzDecoder::new(stdout))?
    } else {
        normalize_onlist_reader(stdout)?
    };

    let status = child
        .wait()
        .with_context(|| format!("failed to wait for curl while reading '{}'", url))?;
    if !status.success() {
        bail!("curl exited with status {} while reading '{}'", status, url);
    }

    Ok(entries)
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
    use std::io::Write;
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
        let loaded = load_onlist(&root, &barcode).unwrap();

        assert!(loaded.entries.contains("AAAA"));
        assert!(loaded.entries.contains("CCCC"));
        assert_eq!(loaded.entries.len(), 2);

        std::fs::remove_dir_all(root).unwrap();
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
    fn test_read_remote_onlist_entries_streams_gzip_over_http() {
        if Command::new("curl")
            .arg("--version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_err()
        {
            return;
        }

        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"AAAA\t1\nCCCC\t2\n").unwrap();
        let body = encoder.finish().unwrap();

        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let addr = listener.local_addr().unwrap();
        let server = thread::spawn(move || {
            let (mut stream, _) = listener.accept().unwrap();
            let response = format!(
                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            );
            stream.write_all(response.as_bytes()).unwrap();
            stream.write_all(&body).unwrap();
        });

        let observed =
            read_remote_onlist_entries(&format!("http://{}/barcodes.txt.gz", addr)).unwrap();

        server.join().unwrap();

        assert!(observed.contains("AAAA"));
        assert!(observed.contains("CCCC"));
        assert_eq!(observed.len(), 2);
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
            input_path: PathBuf::from("R1.fastq.gz"),
            read: assay.get_read("rna_R1").unwrap(),
            matched_file: None,
            matched_by: "read_id".to_string(),
            coordinates,
            spec_base,
        };

        let barcode = filter_regions_by_id(&input, "barcode").unwrap();
        assert_eq!(barcode.len(), 1);
        assert_eq!(barcode[0].region.region_id, "barcode");
    }
}
