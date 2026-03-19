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

- `length`: compare observed read lengths to the read length range in the spec
- `coverage`: report how often each expected region is fully covered by the read
- `fixed`: check exact matches for fixed-sequence regions
- `onlist`: check exact onlist membership for onlist regions
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
- `FASTQ...`: one or more FASTQ files to inspect

FASTQ paths are resolved against the spec in this order: `file_id`, `filename`, `url` basename, then `read_id`.

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
```

## Output

Each metric command emits the same JSON envelope:

- `spec`
- `modality`
- `command`
- `n_reads`
- `files`

Each file entry contains:

- `input_path`
- `read_id`
- `file_id`
- `matched_by`
- `results`

The `results` payload depends on the command:

- `length`: expected and observed read length ranges
- `coverage`: expected coordinates and coverage fractions
- `coverage` also emits warnings when the same `region_id` or biological `region_type` is assigned to multiple reads in one invocation. This is the main guardrail for the failure mode shown in the CRISPR presentation under `docs/`.
- `fixed`: per-region fixed-sequence match counts and top mismatches
- `onlist`: per-region onlist match counts and top offlist sequences
- `cut`: extracted region sequences
- `hist`: counts of extracted region sequences

## Fixtures

- `examples/10xv3/spec.yaml` is a current-format local example used to test index-read behavior.
- `tests/fixtures/synthetic/` contains a small fully local assay used for golden CLI tests.
- `examples/dogmaseq-lll/spec.yaml` is legacy prototype material. It is not the source-of-truth test fixture for the current Rust implementation.

## Notes

The presentation in [docs/Example of seqspec Failing in the CRISPR-pipeline.pptx](docs/Example%20of%20seqspec%20Failing%20in%20the%20CRISPR-pipeline.pptx) shows the motivating failure mode for this tool: a spec can be internally valid but still place biological elements in the wrong reads because of human geometry encoding choices. The `coverage` command now surfaces this explicitly with duplicate-assignment warnings, and `fixed` plus `onlist` check whether the observed read content matches the encoded positions.
