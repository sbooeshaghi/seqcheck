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
reviewer a separate package made by `scripts/manage_cohort_reviews.py`.
Reviewers must not see the other package or decision before both returned sheets
are locked.

Prepare both packages from the authoritative candidate manifest and table:

```bash
uv run python scripts/manage_cohort_reviews.py prepare \
  --candidate-manifest experiments/paper/runs/<run>/manifests/cohort_candidates.json \
  --output-root experiments/paper/runs/<run>/review/packages
```

Give `review/packages/reviewer_1/` to reviewer 1 and
`review/packages/reviewer_2/` to reviewer 2. Each directory contains a review
sheet and a package manifest. The two packages have different content-addressed
package identifiers but the same candidate evidence.

Each reviewer records their name, decision, rationale, protocol URL, and date.
These five generic columns are the only editable fields in a reviewer sheet.
The reviewer must use the same identity on every row and must not edit the
package identifier, slot, or copied candidate evidence. The merge command
requires a complete decision for every candidate; it rejects partial sheets.
Allowed decisions are:

- `include`: the protocol, modality, read structure, and exact FASTQ set support
  this row as an unmodified baseline for the proposed family.
- `exclude`: the row is not an unmodified baseline for the proposed family.
- `inconclusive`: the available protocol or metadata cannot resolve the row.

Merge the two locked sheets with their original package manifests:

```bash
uv run python scripts/manage_cohort_reviews.py merge \
  --candidate-manifest experiments/paper/runs/<run>/manifests/cohort_candidates.json \
  --reviewer-1-package experiments/paper/runs/<run>/review/packages/reviewer_1/review_package.json \
  --reviewer-1-sheet experiments/paper/runs/<run>/review/packages/reviewer_1/cohort_review.csv \
  --reviewer-2-package experiments/paper/runs/<run>/review/packages/reviewer_2/review_package.json \
  --reviewer-2-sheet experiments/paper/runs/<run>/review/packages/reviewer_2/cohort_review.csv \
  --output-root experiments/paper/runs/<run>/review/merged
```

The merge command verifies the exact candidate keys and every copied evidence
value, checks that the reviewers are distinct, and reconstructs candidate fields
from the authoritative candidate table rather than either reviewer sheet. It
writes the combined review table plus a manifest that hashes both locked inputs.
It never overwrites an existing package or merged review.

Preserve the merged output unchanged. Create the adjudication table from that
output and edit only `adjudication_decision`, `adjudication_rationale`, and
`final_family`. If the decisions disagree, or either decision is
`inconclusive`, an adjudicator records `include` or `exclude` plus a rationale.
An included row sets `final_family` to its proposed `family_id`. An excluded row
leaves `final_family` blank. Reviewers and adjudicators leave `split` blank. The
freeze command independently rejects any changed stable candidate evidence.

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
  --reviews experiments/paper/runs/<run>/review/adjudication/cohort_reviews.csv \
  --family-rules docs/cohort_family_rules.json \
  --output-root experiments/paper/runs/<run>/freeze
```

The current feature-tag source has only four configurations with exact FASTQ
reconciliation. Do not weaken the gate to fill the fifth slot. Add another
public source or approve a separately versioned correction that preserves the
original record and its diff.
