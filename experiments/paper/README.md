# Paper Experiment Workspace

This directory separates paper runs from source fixtures and prior exploratory
audits. Git tracks protocols and analysis code. It ignores generated `runs/` and
large `artifacts/`.

Each run directory must contain its own manifests, tables, validation records,
and normalized specs. Do not seed a final run by copying an entire prior output
directory because that can retain stale specs. Raw FASTQ samples belong in a
separately versioned artifact bundle with checksums.

The cohort review and freeze procedure is defined in
`docs/COHORT_REVIEW_PROTOCOL.md`.

Run the Phase 0 acceptance gate before generating study results:

```bash
python3 scripts/run_phase0.py \
  --case-manifest experiments/paper/protocol/phase0_cases.json \
  --output-root experiments/paper/runs/phase0-acceptance-<date> \
  --seqcheck-bin target/release/seqcheck \
  --seqspec-bin ../seqspec/target/release/seqspec \
  --fastqc-command fastqc
```

The command refuses to overwrite a prior run. It must report `valid: true` in
`validation/phase0.json`. The tracked case manifest uses a true 0.3 local spec
and pinned 0.5 HTTPS inputs. Raw and normalized results must match after replacing
only the materialized paths of declared specs and resources; metric values and
assessments are never excluded from the parity check.

The Experiment 1 condition matrix is fixed in
`protocol/sampling_calibration.json`. After the reviewed cohort is frozen,
select one or two nonredundant FASTQs from each calibration configuration:

```bash
python3 scripts/build_sampling_cases.py \
  --cohort-manifest experiments/paper/runs/<cohort-run>/freeze/manifests/cohort_frozen.json \
  --sampling-protocol experiments/paper/protocol/sampling_calibration.json \
  --seqspec-bin ../seqspec/target/release/seqspec \
  --output-root experiments/paper/runs/<run-id>/case-selection
```

The primary FASTQ has the largest indexed `RGN:measure:*` span. A secondary
FASTQ is retained only when it adds a new measure, partition, or technical
ontology term. The selector verifies the frozen cohort and spec hashes, records
declared compressed sizes, and estimates the two source traversals needed for
bounded sampling and the complete-stream reference.

Create all prefix and reservoir samples for one selected FASTQ in one complete
source traversal:

```bash
python3 scripts/sample_fastq_matrix.py \
  --input reads.fastq.gz \
  --output-root experiments/paper/runs/<run-id>/samples/<fastq-id> \
  --n-reads 1000 \
  --n-reads 10000 \
  --n-reads 100000 \
  --seed 17 \
  --seed 29 \
  --seed 43 \
  --include-prefix
```

Repeat `--input` and add `--synchronize-mates` only when the read files form a
paired set that must retain matching names. The sampler consumes every source
to compute its content hash, writes deterministic gzip files, records all 12
sample conditions per FASTQ, and refuses to overwrite an existing output root.
The complete-stream seqcheck report is the thirteenth condition and runs from
the source rather than a retained copy.

Run the complete calibration after case selection:

```bash
python3 scripts/run_sampling_calibration.py \
  --case-selection-manifest experiments/paper/runs/<run-id>/case-selection/manifests/case_selection.json \
  --sampling-protocol experiments/paper/protocol/sampling_calibration.json \
  --seqcheck-bin target/release/seqcheck \
  --output-root experiments/paper/runs/<run-id>/sampling-calibration
```

The runner creates each bounded sample matrix in one source traversal, runs the
complete-stream reference, measures every planned seqcheck replicate, and
exports canonical reports, raw metrics, scalar errors, and performance rows. It
refuses to overwrite a run or accept changed input hashes. The run is usable
only when `validation/sampling_calibration.json` reports `valid: true`.

Analyze a valid calibration and freeze the sampling policy:

```bash
python3 scripts/analyze_sampling_calibration.py \
  --study-manifest experiments/paper/runs/<run-id>/sampling-calibration/manifests/study.json \
  --analysis-protocol experiments/paper/protocol/sampling_analysis.json \
  --output-root experiments/paper/runs/<run-id>/sampling-analysis
```

The analyzer collapses region-level values within each FASTQ and clusters
multiple FASTQs from the same configuration before computing accuracy limits or
bootstrap intervals. `validation/sampling_analysis.json` must report
`valid: true`. The generated policy is ready for later experiments only when
`policy/sampling_policy.json` also reports `frozen: true`; insufficient endpoint
coverage remains a valid analysis but does not pass the policy-freeze gate.
