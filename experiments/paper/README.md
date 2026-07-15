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
