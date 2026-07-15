# Cohort Review Protocol

The paper cohort must contain five verified configurations from each of six
assay families. Automated metadata and modality rules create the review pool,
but they do not approve a baseline configuration.

## Review Unit

One row is a seqspec configuration proposed for one assay family. The same
configuration can appear in more than one family when portal metadata is
ambiguous. At most one row for a configuration can enter the final cohort.

Reviewers compare four records:

1. The normalized seqspec linked in `normalized_spec_path`
2. The assay protocol linked in the reviewer protocol field
3. The portal-linked FASTQ accessions in `fastq_accessions`
4. The seqspec-declared FASTQ accessions in `expected_fastq_accessions`

## Independent Decisions

Two reviewers inspect each row independently. The coordinator gives each
reviewer a copy of `cohort_reviews.template.csv` with the other reviewer's
columns blank. Reviewers must not see the other decision before both copies are
locked.

Each reviewer records their name, decision, rationale, protocol URL, and date.
Allowed decisions are:

- `include`: the protocol, modality, read structure, and exact FASTQ set support
  this row as an unmodified baseline for the proposed family.
- `exclude`: the row is not an unmodified baseline for the proposed family.
- `inconclusive`: the available protocol or metadata cannot resolve the row.

The coordinator merges the two locked copies without changing candidate fields.
If the decisions disagree, or either decision is `inconclusive`, an adjudicator
records `include` or `exclude` plus a rationale. An included row sets
`final_family` to its proposed `family_id`. An excluded row leaves
`final_family` blank. Reviewers and adjudicators leave `split` blank.

## Freeze Gate

Run `scripts/freeze_cohort.py` after the merged review table is complete. The
command refuses to freeze unless:

- Every candidate row has two distinct, complete reviews.
- Every disagreement or inconclusive decision has a final adjudication.
- Included specs use seqspec 0.5.0 and match the proposed modality.
- Included specs pass local structural and resource-aware checks.
- Spec-declared FASTQs exactly match portal-linked FASTQs.
- Configurations and normalized structure/FASTQ keys are unique.
- Exactly five rows enter each family.

The command assigns two calibration and three evaluation configurations per
family by a seeded hash only after all checks pass. It writes 12 calibration and
18 evaluation rows and refuses to overwrite an existing frozen cohort.

```bash
uv run python scripts/freeze_cohort.py \
  --candidate-manifest experiments/paper/runs/<run>/manifests/cohort_candidates.json \
  --reviews experiments/paper/runs/<run>/tables/cohort_reviews.csv \
  --family-rules docs/cohort_family_rules.json \
  --output-root experiments/paper/runs/<run>/freeze
```

The current feature-tag source has only four configurations with exact FASTQ
reconciliation. Do not weaken the gate to fill the fifth slot. Add another
public source or approve a separately versioned correction that preserves the
original record and its diff.
