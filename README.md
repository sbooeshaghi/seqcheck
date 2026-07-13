# seqcheck

`seqcheck` validates sequencing reads against a `seqspec`.

This matters because a `seqspec` can be internally valid and still encode the wrong read geometry. A barcode, UMI, primer anchor, or feature can be placed in the wrong read, shifted to the wrong coordinates, or backed by the wrong whitelist. `seqcheck` reads the FASTQs directly and reports those mismatches.

## What `seqcheck` checks

`seqspec check` validates the `seqspec` file itself.

`seqcheck` validates the reads that the `seqspec` claims to describe.

The current Rust implementation is exact-only:

- no fuzzy barcode matching
- no approximate fixed-sequence matching
- no FastQC-style quality modules

Each run works on one modality and one or more explicit FASTQ paths.

## Installation

For local development:

```bash
cargo test
```

For a standalone install:

```bash
cargo install --path .
```

## Quick start

All read-check commands use the same core interface:

```bash
seqcheck <command> -s SPEC -m MODALITY [-n N] [--format text|json] FASTQ...
```

Arguments:

- `-s, --spec`: path to the `seqspec` YAML file
- `-m, --modality`: modality to inspect
- `-n, --n-reads`: number of reads to inspect per FASTQ. `0` means all reads
- `--format text|json`: human-readable text or structured JSON output
- `--auth-profile`: optional auth profile for remote resources declared in the seqspec
- `FASTQ...`: one or more FASTQ files to inspect

FASTQ paths are resolved against the spec in this order:

1. `file_id`
2. `filename`
3. `url` basename
4. `read_id`

Run the core checks in one report:

```bash
seqcheck check \
  --format json \
  -o report.json \
  -s tests/fixtures/synthetic/spec.yaml \
  -m rna \
  -n 0 \
  tests/fixtures/synthetic/fastqs/synthetic_R1.fastq
```

Render that JSON as a self-contained HTML report:

```bash
seqcheck report -i report.json -o report.html
```

## IGVF Audit Runner

The repo includes a standalone audit helper for IGVF seqspec-backed FASTQ datasets:

```bash
python scripts/igvf_audit.py --output-root tmp/igvf_audit --limit 10 --public-only
```

The runner expects:

- a usable `seqcheck` binary, typically `target/debug/seqcheck` during local development
- a usable `seqspec` command, or a sibling `../seqspec` checkout that the script can call directly

Optional overrides:

```bash
export SEQCHECK_BIN=/path/to/seqcheck
export SEQSPEC_BIN=/path/to/seqspec
```

The runner:

- enumerates released, validated IGVF seqspec configuration files through the portal API
- queries `seqspec` directly on the remote seqspec URL to get the file version, modalities, and expected FASTQ contracts
- runs `seqcheck` directly on the remote seqspec URL and remote FASTQ URLs
- runs `seqcheck check` per modality and stores one JSON report per run
- writes flat `runs`, `diagnostics`, `metrics`, `failures`, and `lab_summary` catalogs
- records seqcheck, seqspec, and FastQC identities in a study manifest
- keys cached reports by tool, input, sampling, and report-schema identity
- validates that reports and flat catalogs reconcile before returning success

Portal seqspec YAMLs are often old tagged `v0.3.x` files. The current remote audit path relies on the toolchain itself to handle those:

- `seqspec version`, `seqspec info`, and `seqspec file` work on remote seqspec URLs
- `seqcheck` upgrades old seqspecs internally before validating reads

The runner no longer downloads or stages local spec or FASTQ inputs.

For controlled-access FASTQs, define a matching `igvf` auth profile in both `seqspec` and `seqcheck`, then export the required credential env vars for those profiles.

If the credentials are not ready, the runner will still audit public data, but controlled-access FASTQs will be skipped.

Example:

```bash
export IGVF_ACCESS_KEY_ID=...
export IGVF_ACCESS_KEY_SECRET=...
```

The runner passes `--auth-profile igvf` to `seqspec` and `seqcheck` only when the local profile exists and its configured env vars are set.

## Study FASTQ Sampler

The study sampler creates bounded prefix or seeded reservoir samples without
changing FASTQ records. Reservoir sampling reads the complete source stream but
retains only the requested records in memory.

```bash
python scripts/sample_fastq.py \
  --input R1.fastq.gz \
  --input R2.fastq.gz \
  --output-root tmp/sample-seed-17 \
  --method reservoir \
  --n-reads 10000 \
  --seed 17 \
  --synchronize-mates
```

`--synchronize-mates` checks normalized read names in lockstep and uses one
shared reservoir for every input. The output manifest records source metadata,
sampling settings, sampler identity, record counts, and output checksums.

## Commands

`seqcheck` has these subcommands:

- `check`: run the core checks in one report
- `length`: compare observed read lengths to the read length range in the spec
- `coverage`: report how often each expected region is fully covered by the read
- `fixed`: check exact matches for fixed-sequence regions
- `onlist`: check exact onlist membership for onlist regions
- `primer`: classify the primer anchor and scan reads for exact primer or reverse-complement hits
- `random`: compute whole-sequence entropy for `sequence_type=random` regions
- `cut`: extract an exact region slice from each covered read
- `hist`: count exact region slices across reads
- `report`: render a self-contained HTML QC report from `seqcheck` JSON
- `version`: print the `seqcheck` version and the seqspec file version

`seqcheck check` runs:

- `input_check`
- `length`
- `coverage`
- `primer`
- `fixed`
- `onlist`
- `random`

It does not run `cut` or `hist`. Those remain drill-down tools.

## HTML Report

`seqcheck report` takes an existing JSON report and writes one offline HTML file:

```bash
seqcheck report -i report.json -o report.html
```

Arguments:

- `-i, --input`: path to a `seqcheck` JSON report
- `-o, --output`: path to the output HTML report
- `--spec`: optional seqspec YAML path to use for the library diagram. By default `seqcheck report` uses `meta.spec` from the JSON report

The HTML report is built from the same atomic result objects used in the JSON output. It includes:

- a header with report metadata
- a run section with the `seqcheck check` and `seqcheck report` commands
- the seqspec library structure diagram when the spec is available
- severity counts
- a flat result list sorted by severity, with expandable expected and observed detail

## Objects and Check Scope

`seqcheck` works on a small set of objects. Each check attaches to the narrowest object that makes sense. These scopes also appear in the JSON report through the `files`, `reads`, and `regions` id lists on each result.

| Object | Meaning in `seqcheck` | Checks on that object |
| --- | --- | --- |
| input set | The supplied FASTQs resolved against the seqspec for one modality | `input_check` |
| resolved read | One matched seqspec file/read pair and its read geometry | `length`, file-level `coverage`, `primer` |
| projected region | One seqspec region projected onto one resolved read | region-level `coverage`, `cut`, `hist` |
| fixed region | A projected region with `sequence_type=fixed` | `fixed` |
| onlist region | A projected region with `sequence_type=onlist` | `onlist` |
| random region | A projected region with `sequence_type=random` | `random` |
| primer anchor | The region referenced by `primer_id` for a resolved read | `primer` |
| shared region across reads | The same logical region, or the same curated biological region type, visible in more than one read | cross-file `coverage` overlap results |

## Check Catalog

The table below summarizes the read-based checks.

| Command | Inputs | Check | Reports |
| --- | --- | --- | --- |
| `length` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Compare each observed read length to the matched read's `min_len` and `max_len` in the seqspec. | Sampled count, expected min/max length, observed min/max length, out-of-range count, out-of-range fraction |
| `coverage` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Project expected regions onto each read and ask whether the observed read is long enough to cover each region. Also emit overlap interpretations when the same logical region appears in multiple reads. | Expected region coordinates, per-region covered count/fraction, full-read coverage count/fraction, overlap assessments |
| `fixed` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Extract each `sequence_type=fixed` slice and compare it exactly to the expected sequence. Matching is strand-aware by default. The command also scans the full read for shifted exact hits of the same motif. | Sampled count, covered count/fraction, short-read count, exact-match count/fraction, orientation counts, offset histogram, absent count, multi-hit count, top non-matching sequences |
| `onlist` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Extract each `sequence_type=onlist` slice and test exact membership in the loaded whitelist. | Sampled count, covered count/fraction, short-read count, exact-onlist count/fraction, offlist count, top offlist sequences, onlist source, onlist size |
| `primer` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Classify the `primer_id` anchor for each read, then scan the full read for exact primer-sequence and reverse-complement hits. | Primer classification, primer sequence/length, start-hit and internal-hit counts/fractions for forward and reverse-complement matches, absent count/fraction, top hit positions |
| `random` | `spec`, `modality`, one or more FASTQs, optional `n_reads` | Extract each `sequence_type=random` slice, group exact observed sequences, and compute Shannon entropy over that sequence distribution. | Sampled count, covered count/fraction, short-read count, unique-sequence count, entropy in bits, maximum possible entropy, entropy fraction, top observed sequences |
| `cut` | `spec`, `modality`, one or more FASTQs, required `region_id`, optional `n_reads` | Extract the exact sequence slice for one named region from each covered read. | Sampled count, covered count, short-read count, extracted sequence per covered read |
| `hist` | `spec`, `modality`, one or more FASTQs, required `region_id`, optional `n_reads` | Extract the exact sequence slice for one named region and count exact sequence occurrences across reads. | Sampled count, covered count, short-read count, exact histogram of observed region sequences |

`version` is a utility command, not a read check. It reports the `seqcheck` version, the seqspec file version, and the assay id.

## Check Outcomes

Not every check uses all four outcomes. The table below gives the usual meaning of each outcome for each check.

| Check | What it checks | `Pass` means | `Warning` means | `Error` means | `Interpretation` means |
| --- | --- | --- | --- | --- | --- |
| `input_check` | Matches the supplied FASTQs to seqspec file ids and read ids. Also checks whether primer anchors are scannable. | The expected files matched uniquely. | A read points to a primer anchor that cannot be scanned as sequence. | The supplied FASTQs do not match the seqspec. | A ghost primer anchor is allowed for projection, but there is no sequence motif to scan. |
| `length` | Compares observed read lengths to the seqspec `min_len` and `max_len`. | The sampled reads stay within the declared range. | Some reads fall outside the declared range. | Usually not emitted here. | Usually not emitted here. |
| `coverage` | Checks whether the reads cover the projected read geometry and named regions. | The read covers the full geometry or named region. | Coverage is incomplete. | Usually not emitted here. | A shared region type appears in multiple reads, or another geometry pattern is being reported. |
| `primer` | Checks the primer anchor named by `primer_id`. | A scannable primer is present as expected. | The primer is not scannable or the expected motif is absent. | Usually not emitted here. | The primer anchor is informative, but it cannot be treated as a sequence motif. |
| `fixed` | Checks fixed regions such as linkers or adapters. | The fixed region matches exactly in covered reads. | The exact match fraction is low or absent. | Usually not emitted here. | The region matches only part of the time, often in one dominant orientation. |
| `onlist` | Checks barcodes, indices, or other onlist-backed regions against their whitelist. | The covered sequences are onlist. | Some covered sequences are offlist. | The onlist resource could not be loaded. | Usually not emitted here. |
| `random` | Summarizes the sequence distribution of random regions such as UMI, barcode, cDNA, or scaffold sequence. | Usually not the main output of this check. | Usually not emitted here. | Usually not emitted here. | Reports entropy and the most frequent observed sequences. |

## Report Format

Every read-check command emits the same top-level JSON object:

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

`expected` and `observed` use the same metric-item shape:

- `id`
- `name`
- `description`
- `data`

`data` is shallow and typed:

- `kind`: `scalar`, `series`, or `records`
- `value`
- `unit` when it matters

`assessment` carries the interpretation:

- `type`: `pass`, `warning`, `error`, or `interpretation`
- `code`
- `description`
- `expected_ids`
- `observed_ids`

This makes the JSON the source of truth. Stdout is rendered from the same result objects, so the human-readable report and the machine-readable report describe the same checks.

Each successful report emits an `input_check` result first. It records:

- the expected FASTQ inventory from the seqspec for that modality
- the supplied inputs
- the resolved matches
- any expected files that were not supplied

Missing expected files are warnings, not hard failures, so subset runs remain valid.

`coverage` emits three kinds of results:

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

## Remote Auth

Some seqspec files point at remote whitelist or reference resources that require HTTP auth. `seqcheck` resolves remote auth by host and applies it through one shared fetch path.

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

You can also initialize the same profile directly:

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

Hidden auth helpers for local debugging:

```bash
seqcheck auth path
seqcheck auth list
seqcheck auth resolve https://api.data.igvf.org/reference-files/...
```

`onlist` reads remote whitelist files directly from their URL. It does not require a local staged copy. If a whitelist cannot be loaded, `seqcheck onlist` emits a structured error result instead of failing the whole report.

## Examples

Check read lengths with the DOGMA fixture in the sibling `seqspec` repo:

```bash
seqcheck length \
  --format json \
  -s ../seqspec/tests/fixtures/spec.yaml \
  -m rna \
  -n 100 \
  ../seqspec/tests/fixtures/fastqs/rna_R1_SRR18677638.fastq.gz \
  ../seqspec/tests/fixtures/fastqs/rna_R2_SRR18677638.fastq.gz
```

Run the core checks in one pass:

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

## Fixtures

- `examples/10xv3/spec.yaml` is a current-format local example used to test index-read behavior
- `tests/fixtures/synthetic/` contains a small local assay used for golden CLI tests
- `tests/fixtures/primer_cases/` covers fixed, ghost, and non-scannable primer definitions
- `tests/fixtures/neg_fixed/` covers reverse-complement fixed matching on a negative-strand read
- `examples/dogmaseq-lll/spec.yaml` is legacy prototype material and is not the source-of-truth fixture for the current Rust implementation

## Notes

The presentation in [docs/Example of seqspec Failing in the CRISPR-pipeline.pptx](docs/Example%20of%20seqspec%20Failing%20in%20the%20CRISPR-pipeline.pptx) shows the motivating failure mode for this tool: a spec can be internally valid and still miss the effective read geometry that the data reveals.

The `coverage` command surfaces one common case. A logical region can appear in multiple reads. That can reflect intended paired-end overlap, or it can reflect a read-geometry mismatch when the observed reads extend farther than the spec implies.

The `primer` command surfaces another case. Sequencing usually starts at the primer, so the primer sequence itself should usually not appear inside the read. If primer hits are common, or if the primer definition is not a concrete fixed sequence, that is strong evidence that the seqspec is anchoring the read incorrectly.

## Development Note

Codex was used to help build out parts of the Rust `seqspec` library and the `seqcheck` CLI, including implementation, testing, and documentation.
