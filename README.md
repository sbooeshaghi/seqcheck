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
- `FASTQ...`: one or more FASTQ files to inspect

FASTQ paths are resolved against the spec in this order: `file_id`, `filename`, `url` basename, then `read_id`.

Each successful metric report also includes an `input_check` block. It records the expected FASTQ inventory from the seqspec for that modality, the supplied inputs, the resolved matches, and any expected files that were not supplied. Missing expected files are warnings, not hard failures, so subset runs remain valid.

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

Each metric command emits the same JSON envelope:

- `spec`
- `modality`
- `command`
- `n_reads`
- `input_check`
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
- `coverage` also emits warnings when the same `region_id` or biological `region_type` appears in multiple reads in one invocation. This can reflect intended overlapping paired-end reads, or a seqspec/read-geometry mismatch when the observed reads extend farther than expected. This is the main guardrail for the failure mode shown in the CRISPR presentation under `docs/`.
- `fixed`: per-region fixed-sequence match counts, strand-aware orientation counts, shifted-hit offset histograms, and top mismatches
- `onlist`: per-region onlist match counts and top offlist sequences
- `onlist` reads remote whitelist files directly from their URL. It does not require you to stage the file locally first.
- `primer`: primer classification, exact primer-hit fractions, and dominant hit positions
- `random`: per-region exact-sequence counts, Shannon entropy in bits over the observed sequence distribution, and entropy as a fraction of the theoretical `2 * region_length` DNA maximum
- `cut`: extracted region sequences
- `hist`: counts of extracted region sequences

## Fixtures

- `examples/10xv3/spec.yaml` is a current-format local example used to test index-read behavior.
- `tests/fixtures/synthetic/` contains a small fully local assay used for golden CLI tests.
- `tests/fixtures/primer_cases/` covers fixed, ghost, and non-scannable primer definitions.
- `tests/fixtures/neg_fixed/` covers reverse-complement fixed matching on a negative-strand read.
- `examples/dogmaseq-lll/spec.yaml` is legacy prototype material. It is not the source-of-truth test fixture for the current Rust implementation.

## Notes

The presentation in [docs/Example of seqspec Failing in the CRISPR-pipeline.pptx](docs/Example%20of%20seqspec%20Failing%20in%20the%20CRISPR-pipeline.pptx) shows the motivating failure mode for this tool: a spec can be internally valid but still miss the effective read geometry that the data reveals. The `coverage` command surfaces that by flagging logical regions that appear in multiple reads, and `fixed` plus `onlist` check whether the observed read content matches the encoded positions.

The `primer` command adds another direct geometry check. Sequencing usually starts at the primer, so the primer sequence itself should usually not appear inside the read. If primer hits are common, or if the primer definition is not a concrete fixed sequence, that is a strong signal that the seqspec is anchoring the read incorrectly.
