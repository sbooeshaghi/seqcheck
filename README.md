# seqcheck

`seqcheck` checks whether observed FASTQ reads are consistent with a `seqspec` file.

This matters because a `seqspec` can pass `seqspec check` and still describe the wrong read geometry. The spec may be valid YAML, and the internal region tree may be structurally consistent, but a human may still encode the assay logic incorrectly. That can place a barcode, UMI, or feature in the wrong read, or even duplicate the same logical element across multiple reads. `seqcheck` inspects the reads directly and reports those mismatches.

## Installation

For local development:

```bash
cargo test
```

For a standalone install:

```bash
cargo install --path .
```

## Scope

`seqspec check` validates the `seqspec` file itself.

`seqcheck` validates the reads that the `seqspec` claims to describe.

The current Rust implementation is exact-only. It does not do fuzzy barcode matching, approximate fixed-sequence matching, or FastQC-style quality modules. Each invocation works on one modality and one or more explicit FASTQ paths.

## Commands

`seqcheck` has these subcommands:

- `check`: run the core read-validation checks in one report
- `length`: compare observed read lengths to the read length range in the spec
- `coverage`: report how often each expected region is fully covered by the read
- `fixed`: check exact matches for fixed-sequence regions
- `onlist`: check exact onlist membership for onlist regions
- `primer`: classify the primer anchor and scan reads for exact primer or reverse-complement hits
- `random`: summarize whole-sequence entropy for `sequence_type=random` regions by grouping exact observed sequences and comparing that entropy to the theoretical DNA maximum for the region length
- `cut`: extract an exact region slice from each covered read
- `hist`: count exact region slices across reads
- `version`: print the `seqcheck` version and the `seqspec` file version

All metric commands use the same core arguments:

```bash
seqcheck <command> -s SPEC -m MODALITY [-n N] [--format text|json] FASTQ...
```

- `-s, --spec`: path to the `seqspec` YAML file
- `-m, --modality`: modality to inspect
- `-n, --n-reads`: number of reads to inspect per FASTQ. `0` means all reads
- `--format text|json`: human-readable text or structured JSON output
- `--auth-profile`: optional auth profile for remote resources declared in the seqspec
- `FASTQ...`: one or more FASTQ files to inspect

FASTQ paths are resolved against the spec in this order: `file_id`, `filename`, `url` basename, then `read_id`.

`seqcheck check` runs the core QC surface in one report:

- `input_check`
- `length`
- `coverage`
- `primer`
- `fixed`
- `onlist`
- `random`

It does not include `cut` or `hist`, which remain explicit drill-down tools.

Each successful metric report emits an atomic `input_check` result first. It records the expected FASTQ inventory from the seqspec for that modality, the supplied inputs, the resolved matches, and any expected files that were not supplied. Missing expected files are warnings, not hard failures, so subset runs remain valid.

## Remote Auth

Some seqspec files point at remote whitelist or reference resources that require HTTP auth. `seqcheck` now resolves remote auth by host and applies it through one shared fetch path.

Auth profiles live in `auth.toml` and point to env vars, not raw secrets on the command line. `seqcheck` looks for the config in this order:

1. `SEQCHECK_AUTH_CONFIG`
2. `$XDG_CONFIG_HOME/seqcheck/auth.toml`
3. `$HOME/.config/seqcheck/auth.toml`

Example config:

```toml
[profiles.igvf]
hosts = ["api.data.igvf.org", "data.igvf.org"]
kind = "basic"
username_env = "IGVF_ACCESS_KEY_ID"
password_env = "IGVF_ACCESS_KEY_SECRET"
```

Or initialize the same profile directly:

```bash
seqcheck auth init \
  --profile igvf \
  --host api.data.igvf.org \
  --host data.igvf.org \
  --kind basic \
  --username-env IGVF_ACCESS_KEY_ID \
  --password-env IGVF_ACCESS_KEY_SECRET
```

Example usage:

```bash
export IGVF_ACCESS_KEY_ID=...
export IGVF_ACCESS_KEY_SECRET=...

seqcheck onlist \
  --auth-profile igvf \
  -s spec.yaml \
  -m crispr \
  sample_R1.fastq.gz
```

If `--auth-profile` is omitted, `seqcheck` matches profiles by URL host. If a named profile is supplied, the host must still match that profile.

There is also a hidden inspector for local debugging:

```bash
seqcheck auth path
seqcheck auth list
seqcheck auth resolve https://api.data.igvf.org/reference-files/...
```

## Check Catalog

The table below summarizes the read-based checks in `seqcheck`.

| Command | Inputs | Check | Reports |
| --- | --- | --- | --- |
| `length` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Compare each observed read length to the matched read's `min_len` and `max_len` in the seqspec. | Sampled count, expected min/max length, observed min/max length, out-of-range count, out-of-range fraction. |
| `coverage` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Project expected regions onto each read and ask whether the observed read is long enough to cover each region. Also flag logical regions that appear in multiple reads. | Expected region coordinates, per-region covered count/fraction, full-read coverage count/fraction, overlap/geometry warnings. |
| `fixed` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Extract each `sequence_type=fixed` slice and compare it exactly to the expected sequence. Matching is strand-aware by default, and the command also scans the full read for shifted exact hits of the same motif. | Sampled count, covered count/fraction, short-read count, exact-match count/fraction, orientation match counts, offset histogram, absent count, multi-hit count, top non-matching sequences. |
| `onlist` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Extract each `sequence_type=onlist` slice and test exact membership in the loaded whitelist. | Sampled count, covered count/fraction, short-read count, exact-onlist count/fraction, offlist count, top offlist sequences, onlist source, onlist size. |
| `primer` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Classify the `primer_id` anchor for each read, then scan the full read for exact primer-sequence and reverse-complement hits. | Primer classification, primer sequence/length, start-hit and internal-hit counts/fractions for forward and reverse-complement matches, absent count/fraction, top hit positions. |
| `random` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Extract each `sequence_type=random` slice, group exact observed sequences, and compute Shannon entropy over that sequence distribution. | Sampled count, covered count/fraction, short-read count, unique-sequence count, entropy in bits, maximum possible entropy, entropy fraction, top observed sequences. |
| `cut` | `spec`, `modality`, one or more FASTQs, required `region_id`, optional `n_reads` | Extract the exact sequence slice for one named region from each covered read. | Sampled count, covered count, short-read count, extracted sequence per covered read. |
| `hist` | `spec`, `modality`, one or more FASTQs, required `region_id`, optional `n_reads` | Extract the exact sequence slice for one named region and count exact sequence occurrences across reads. | Sampled count, covered count, short-read count, exact histogram of observed region sequences. |

`version` is a utility command, not a read check. It reports the `seqcheck` version, the seqspec file version, and the assay id.

## Examples

Check read lengths with the current DOGMA fixture in the sibling `seqspec` repo:

```bash
seqcheck length \
  --format json \
  -s ../seqspec/tests/fixtures/spec.yaml \
  -m rna \
  -n 100 \
  ../seqspec/tests/fixtures/fastqs/rna_R1_SRR18677638.fastq.gz \
  ../seqspec/tests/fixtures/fastqs/rna_R2_SRR18677638.fastq.gz
```

Run the core checks in one pass over the command surface:

```bash
seqcheck check \
  --format json \
  -s tests/fixtures/synthetic/spec.yaml \
  -m rna \
  -n 0 \
  tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
```

Check exact onlist membership for the local index read in the upgraded `10xv3` example:

```bash
seqcheck onlist \
  --format json \
  -s examples/10xv3/spec.yaml \
  -m rna \
  -n 100 \
  examples/10xv3/fastqs/I1.fastq.gz
```

Inspect a synthetic barcode slice and its histogram:

```bash
seqcheck cut \
  -s tests/fixtures/synthetic/spec.yaml \
  -m rna \
  -r barcode \
  -n 0 \
  tests/fixtures/synthetic/fastqs/synthetic_R1.fastq

seqcheck hist \
  -s tests/fixtures/synthetic/spec.yaml \
  -m rna \
  -r barcode \
  -n 0 \
  tests/fixtures/synthetic/fastqs/synthetic_R1.fastq

seqcheck random \
  --format json \
  -s tests/fixtures/synthetic/spec.yaml \
  -m rna \
  -n 0 \
  tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
```

## Output

Each read-check command now emits the same top-level JSON object:

- `report_schema_version`
- `meta`
- `results`

`meta` contains:

- `command`
- `spec`
- `modality`
- `requested_reads`

`results` is a flat list of atomic result objects. A single command can emit many results. Each result has:

- `check`
- `files`
- `reads`
- `regions`
- `expected`
- `observed`
- `assessment`

`files`, `reads`, and `regions` are seqspec ids only. They point back to entities already defined in the spec.

`expected` and `observed` both use the same metric-item shape:

- `id`
- `name`
- `description`
- `data`

`data` is shallow and typed:

- `kind`: `scalar`, `series`, or `records`
- `value`
- `unit` when it matters

`assessment` contains structured interpretation:

- `type`: `pass`, `warning`, `error`, or `interpretation`
- `code`
- `description`
- `expected_ids`
- `observed_ids`

This means the JSON is the source of truth. Stdout is rendered from the same atomic results, so a human and an AI agent see the same underlying report.

For example, `coverage` emits:

- one `input_check` result per invocation
- one file-level `coverage` result per matched file/read
- one region-level `coverage` result per matched file/read/region
- cross-file `coverage` results when the same logical region or curated biological region type appears in multiple reads

Abbreviated JSON example:

```json
{
  "report_schema_version": "0.1.0",
  "meta": {
    "command": "coverage",
    "spec": "../seqspec/tests/fixtures/spec.yaml",
    "modality": "rna",
    "requested_reads": 100
  },
  "results": [
    {
      "check": "input_check",
      "files": ["rna_R1_SRR18677638.fastq.gz", "rna_R2_SRR18677638.fastq.gz"],
      "reads": ["rna_R1", "rna_R2"],
      "regions": [],
      "expected": [...],
      "observed": [...],
      "assessment": [...]
    },
    {
      "check": "coverage",
      "files": ["rna_R1_SRR18677638.fastq.gz"],
      "reads": ["rna_R1"],
      "regions": ["rna_umi"],
      "expected": [...],
      "observed": [...],
      "assessment": [...]
    }
  ]
}
```

Abbreviated stdout for the same schema:

```text
seqcheck coverage
spec: ../seqspec/tests/fixtures/spec.yaml
modality: rna
requested_reads: 100

check: input_check
files: rna_R1_SRR18677638.fastq.gz, rna_R2_SRR18677638.fastq.gz
reads: rna_R1, rna_R2
regions: (none)
assessment:
  - [pass] all_expected_files_matched: All expected modality files were supplied and matched uniquely.

check: coverage
files: rna_R1_SRR18677638.fastq.gz
reads: rna_R1
regions: rna_umi
assessment:
  - [pass] region_coverage_complete: All sampled reads cover region 'rna_umi'.
expected:
  - expected_region:
    - name=UMI region_id=rna_umi region_type=umi sequence_type=random start=16 stop=28
observed:
  - sampled_count: 100 count
  - covered_count: 100 count
  - covered_fraction: 1 fraction
```

`onlist` reads remote whitelist files directly from their URL. It does not require you to stage the file locally first. If a whitelist cannot be loaded, `seqcheck onlist` now emits a structured error result instead of failing the whole report.

## Fixtures

- `examples/10xv3/spec.yaml` is a current-format local example used to test index-read behavior.
- `tests/fixtures/synthetic/` contains a small fully local assay used for golden CLI tests.
- `tests/fixtures/primer_cases/` covers fixed, ghost, and non-scannable primer definitions.
- `tests/fixtures/neg_fixed/` covers reverse-complement fixed matching on a negative-strand read.
- `examples/dogmaseq-lll/spec.yaml` is legacy prototype material. It is not the source-of-truth test fixture for the current Rust implementation.

## Notes

The presentation in [docs/Example of seqspec Failing in the CRISPR-pipeline.pptx](docs/Example%20of%20seqspec%20Failing%20in%20the%20CRISPR-pipeline.pptx) shows the motivating failure mode for this tool: a spec can be internally valid but still miss the effective read geometry that the data reveals. The `coverage` command surfaces that by flagging logical regions that appear in multiple reads, and `fixed` plus `onlist` check whether the observed read content matches the encoded positions.

The `primer` command adds another direct geometry check. Sequencing usually starts at the primer, so the primer sequence itself should usually not appear inside the read. If primer hits are common, or if the primer definition is not a concrete fixed sequence, that is a strong signal that the seqspec is anchoring the read incorrectly.
