# Seqcheck Paper Experimental Plan

Status: draft protocol
Date: July 13, 2026

## Study Question

Seqcheck tests whether sequencing reads agree with the assay structure declared in
a seqspec file. For example, a seqspec can be valid YAML and still put a cell
barcode in the wrong read or assign the wrong coordinates to a linker. Generic
FASTQ quality control cannot test that contract because it does not know the
assay structure.

The study asks:

> Can a machine-readable assay specification support accurate, localized, and
> scalable checks that observed reads implement the declared measurement?

The study will test four claims:

1. Seqcheck detects and localizes known read-versus-spec inconsistencies.
2. Seqcheck complements rather than replaces generic FASTQ quality control.
3. Region ontology terms make equivalent regions comparable across independently
   named seqspec files and support better role-specific interpretation.
4. Specification-aware QC is practical at consortium scale using bounded read
   samples rather than complete FASTQ files.

The study will not claim that every warning is an experimental defect. A warning
is a measured inconsistency until a frozen policy or expert adjudication gives it
a stronger interpretation.

## Outcome Rules

Each experiment distinguishes three kinds of result:

| Result | Meaning |
| --- | --- |
| Completion check | Confirms that the planned observations exist and their denominators reconcile |
| Expected result | States the directional hypothesis before evaluation data are examined |
| Scientific target | Gives a numerical result that would support the intended paper claim |

An experiment is complete when its completion checks pass, even when its expected
result or scientific target is not met. A missed target changes the claim; it does
not justify changing the target after looking at the evaluation data.

Thresholds and applicability rules may change during calibration. Once the
evaluation cohort is unlocked, any change to a tool, policy, or endpoint creates
a new study run. The complete evaluation must then be repeated under the new run
identifier. Results from different locked runs will not be pooled.

Primary rates use fixed denominators:

| Metric | Denominator |
| --- | --- |
| Sensitivity | Applicable perturbed configuration-operator pairs |
| Specificity | Verified unmodified configuration-check pairs |
| Localization accuracy | Detected perturbations for which file, read, or region localization applies |
| End-to-end localization | All applicable perturbations, including missed detections |
| Mapping precision | Manually adjudicated ontology mappings with sampling weights |
| Audit completion | Eligible configuration-modality runs in the frozen portal manifest |
| Confirmation precision | Reviewed high-confidence candidates, with conservative and evaluable versions reported |

Calculate each rate within configuration first and then macro-average across
configurations. Confidence intervals resample configurations, not individual
mutations or regions.

## Experimental Units

The primary independent unit is a seqspec configuration and its matched FASTQ
set for one modality. Multiple perturbations derived from the same configuration
are paired observations, not independent replicates.

The controlled benchmark will contain 30 public configurations, with five from
each of six assay families:

1. Droplet-partitioned RNA measurement
2. Combinatorial-indexing RNA sequencing
3. Single-cell or single-nucleus ATAC sequencing
4. Joint RNA and chromatin multiome sequencing
5. Perturbation or guide-capture sequencing
6. Protein, feature, or sample-tag barcoding

Configurations must have a loadable seqspec, complete expected FASTQ metadata,
and accessible referenced resources. Where possible, each family will include
at least two submitting laboratories. Configurations with identical normalized
library structure and identical FASTQ sets will be deduplicated.

Two reviewers must verify each benchmark pair against the available assay
protocol, expected read layout, and linked FASTQ metadata before it can serve as
an unmodified control. Unresolved pairs will be excluded before the calibration
and evaluation split. Cohort selection will not use seqcheck warning counts or
controlled perturbation performance.

Build the review pool with `scripts/build_cohort_candidates.py` and the versioned
rules in `docs/cohort_family_rules.json`. Selection has two stages. First, exact
portal assay titles and assay terms produce a lab-balanced metadata proposal
pool. The builder then normalizes every proposal and uses the seqspec modality
to produce a second lab-balanced review pool. For example, a `Perturb-seq` title
can describe either a guide read or its companion RNA read, so only a `crispr`
or `guide` modality can enter the perturbation review pool. A normalization
failure cannot enter the review pool because its modality is unknown, but it
remains in the proposal table with the failure reason.

Both stages balance submitting laboratories first and portal-linked FASTQ counts
within each laboratory second. FASTQ count is input metadata, not a validation
result. This stratification keeps uncommon index-read layouts represented
without favoring a specification because its declared FASTQs happen to match.

Structural checks, resource checks, FASTQ reconciliation, and seqcheck results
do not affect which modality-qualified proposals enter the review pool. They are
diagnostic fields for review and later exclusion. This rule avoids selecting
specifications because they happen to pass the tools being evaluated. The
builder writes both `cohort_proposals.csv` and `cohort_candidates.csv` so every
proposal, exclusion, and selected row can be reconciled.

The frozen portal snapshot contains 18 released and validated configurations
titled `scRNA-seq`, but all 18 point to controlled FASTQs. Public droplet RNA
proposals therefore include RNA configurations from `10x multiome` and `10x
multiome with MULTI-seq`. This family describes the RNA read measurement, not an
RNA-only experiment. The actual `rna` modality and human protocol review still
have to agree.

The feature-tag source has a separate known issue. In the current snapshot, 105
of 109 public configurations link an i2 FASTQ that the seqspec does not declare;
only four reconcile exactly. Exact reconciliation must not be weakened to fill
the family quota. Reviewers must either find another public source or preserve
the original spec and approve an explicit, versioned correction before a fifth
configuration can serve as a baseline. The uncorrected record remains an audit
observation.

The proposed `IGVFFI7663RZZB` correction is registered as immutable review
evidence. Both independent reviewers must approve that exact correction; an
adjudicator cannot override either reviewer for a corrected baseline. The freeze
gate revalidates the corrected seqspec and recomputes its FASTQ, structure, and
deduplication evidence while preserving the original candidate fields.

The builder normalizes proposals to seqspec 0.5.0, records schema/structural
validation separately from network-dependent resource validation, fingerprints
the library/read structure without FASTQ bindings, and compares spec-declared
FASTQs with portal-linked FASTQs. Resource validation runs only after a spec
passes the internal checks, so an unavailable endpoint cannot hide a stable
structural failure. Volatile resource outcomes and retry counts are recorded but
do not change the candidate selection identifier. The builder writes a blank
combined-review template and marks the cohort as unfrozen. A rerun never
overwrites an existing review sheet and reports when that sheet belongs to a
different selection identifier.

Reviewers follow `docs/COHORT_REVIEW_PROTOCOL.md`.
`scripts/manage_cohort_reviews.py` creates two independent, content-addressed
review packages and merges only complete sheets from distinct reviewers. The
merge verifies copied evidence against the candidate table and reconstructs the
combined table from the authoritative candidate rows. After both reviews and
any adjudication are locked, `scripts/freeze_cohort.py` reconciles the final
review table to the candidate manifest and assigns the seeded 12/18 split. The
command writes no frozen cohort when any review, baseline check, family quota,
or deduplication rule fails.

Two configurations from each family will form a 12-configuration calibration
set. The remaining 18 configurations will be locked before thresholds are chosen
and used only for evaluation. Existing repository fixtures will be used for
software tests, not as paper benchmark observations.

Thirty configurations is the minimum cohort, not a fixed maximum. Before opening
evaluation results, use the frozen applicability matrix to confirm that each
primary perturbation has at least 12 applicable evaluation configurations and
that a configuration-level simulation gives a projected 95% interval half-width
of at most 0.15 for the primary sensitivity estimate. Add verified configurations
before evaluation if either check fails. Assay-family results with only three
evaluation configurations will be descriptive rather than inferential.

## Common Sampling Design

The default audit sample is 10,000 records from each FASTQ, not 10,000 records
across the complete configuration. This keeps index, barcode, and biological
reads represented independently.

The study will use three sample sizes:

| Sample | Records per FASTQ | Purpose |
| --- | ---: | --- |
| Screen | 10,000 | Default controlled benchmark and IGVF audit |
| Escalate | 100,000 | Rare events, unstable estimates, and borderline cases |
| Complete stream | All records | Sampling validation on a small selected cohort only |

The current `--n-reads` behavior reads the first records from each FASTQ. This is
efficient for remote compressed files, but file order could bias the result. The
sampling experiment will compare FASTQ prefixes with seeded reservoir samples
from the complete stream. Reservoir sampling reads the complete input but does
not retain it on disk.

The controlled benchmark will persist only the bounded public samples needed to
reproduce the paper. Perturbed reads will be generated into temporary files or
streams and removed after their reports are written. The full IGVF audit will
stream remote FASTQs and retain reports, metrics, and manifests rather than read
data.

Every sample record must include:

- Source configuration, modality, FASTQ accession, read id, and URL
- Sampling method, requested count, observed count, and seed
- Source content identifier or checksum when one is available
- Retrieval time and access status
- Seqspec, seqcheck, and sampler commit or version

Paired FASTQ samples used for downstream analyses must preserve read names and
select the same names from every mate.

## Phase 0: Make Runs Reproducible

These tasks must be complete before generating paper results.

| Task | Implemented evidence | Completion criterion |
| --- | --- | --- |
| Normalize audit specs to seqspec 0.5.0 | `igvf_audit.py` and `run_phase0.py` record raw and normalized versions; the acceptance gate upgrades a true 0.3 fixture | A test verifies old portal specs normalize to 0.5.0 |
| Record run provenance | Study manifests record exact input, tool, script, invocation, runtime, sampling, and report identities | Each run records tool commits, command, sampling method, seed, UTC time, and report schema |
| Export raw metrics | The audit and Phase 0 runners flatten every expected and observed metric with scope and ontology terms | One row per expected or observed metric is written with result scope and ontology terms |
| Prevent stale cache reuse | Cache identities include executable or script hashes, input identity, sampling policy, and report schema | Cache keys include tool versions, input identity, sampling settings, and report schema |
| Add a study sampling harness | `sample_fastq.py` handles individual probes; `sample_fastq_matrix.py` creates every planned prefix and reservoir condition in one complete traversal with synchronized-mate support and a shared manifest | Tested scripts write seeded reservoir samples, preserve full FASTQ records, synchronize mates when requested, and record stable content identities |
| Freeze external tools | Audit and Phase 0 manifests capture the FastQC executable, version, hash, and limits configuration | FastQC version and configuration are captured in the study manifest |
| Create a clean output contract | New run roots contain immutable inputs, reports, tables, manifests, and fail-closed reconciliation records | Every table reconciles to a manifest entry and is generated in a new run directory |

Required validation for this phase:

```text
cargo test
cargo clippy --all-targets --all-features -- -D warnings
python3 -m unittest discover -s tests
python3 scripts/run_phase0.py --case-manifest experiments/paper/protocol/phase0_cases.json ...
```

The integrated runner materializes and hashes every spec, normalizes each case
to 0.5.0, runs raw and normalized reports, exports every expected and observed
metric, and checks all row denominators. It also repeats a seeded reservoir
sample, changes the seed, and probes cache identities after changing the tool or
sampling policy. The gate compares complete result objects. It canonicalizes
only materialized paths for declared specs and resources so local path changes
do not masquerade as measurement changes.

### Phase 0 Goal and Check

| Item | Definition |
| --- | --- |
| Measurable goal | Produce one self-describing report bundle from a legacy 0.3 spec, a current 0.5 spec, local FASTQs, and remote FASTQs |
| Expected result | Existing check metrics remain unchanged apart from intentional ontology grouping and added provenance |
| Scientific target | Not applicable; this is an engineering gate |
| Completion check | All tests pass; both old specs normalize to 0.5.0; metric-table row counts equal the metrics in source reports; identical sampling seeds produce identical sample checksums; changed tool or sampling settings invalidate the cache |
| If unmet | Fix the instrumentation and repeat Phase 0 before selecting the paper cohort |

## Experiment 1: Sampling Accuracy and Runtime

### Question

How many reads are needed for stable seqcheck metrics, and does FASTQ-prefix
sampling differ materially from uniform sampling over the complete file?

### Goal and Check

| Item | Definition |
| --- | --- |
| Measurable goal | Choose the smallest sample that estimates common seqcheck fractions accurately and determine whether prefix sampling is biased |
| Primary endpoint | Absolute error from the complete-stream metric and paired prefix-minus-reservoir error |
| Expected result | A 10,000-read sample is adequate for common events, while rare-event and entropy estimates sometimes require 100,000 reads |
| Scientific target | For common fraction metrics at 10,000 reads, median absolute error is at most 0.01 and the 95th percentile is at most 0.05; the 95% confidence interval for mean prefix-minus-reservoir error lies inside -0.01 to 0.01; peak local memory at 100,000 reads is below 1 GiB |
| Completion check | Every selected FASTQ has a complete-stream reference, all planned prefix samples, three reservoir seeds per size, and reconciled metric and performance rows |
| If unmet | Use 100,000 reads for the affected metric or replace prefix sampling in the main audit; report the failed 10,000-read result |

### Cohort

Use 12 benchmark configurations, two from each assay family. Select one or two
FASTQs from each configuration to cover biological inserts, barcode or index
reads, and fixed technical regions. Cap the complete-stream cohort at 24 FASTQs.
Select files with manageable public transfer sizes so complete streaming is
feasible without retaining full files.

The selection is deterministic. `build_sampling_cases.py` chooses the FASTQ
with the largest indexed `RGN:measure:*` span, then retains one additional FASTQ
only when it adds an ontology role or term. This favors a biological measurement
read plus a complementary partition or technical read without paying to stream
redundant files. The selection manifest records each declared compressed size
and the minimum transfer for the bounded-sample and complete-reference streams.

### Procedure

For each selected FASTQ, run seqcheck at 1,000, 10,000, and 100,000 reads using:

- Prefix sampling
- Three reservoir samples with fixed seeds 17, 29, and 43, produced with all
  prefix samples during one complete input stream
- The complete stream as the reference

The machine-readable condition definition is
`experiments/paper/protocol/sampling_calibration.json`. It defines 12 bounded
sample conditions and one complete-stream reference per FASTQ. The sample
matrix manifest must contain all 12 bounded conditions before seqcheck runs.

Execute the frozen matrix with `scripts/run_sampling_calibration.py`. The runner
verifies the case-selection and protocol hashes, creates each bounded matrix in
one source traversal, runs the complete reference, and writes reconciled report,
metric, scalar-error, and performance tables. A run advances to analysis only
when `validation/sampling_calibration.json` reports `valid: true`.

Run each local performance condition three times after one warm-up run. Measure
wall time, peak resident memory, records processed, compressed bytes read when
available, and report size. Report local computation and remote transfer
separately.

Compare these observed metrics with the complete-stream result:

- Out-of-range read fraction
- Full-read and region coverage fractions
- Fixed-sequence exact-match fraction
- Onlist fraction
- Primer hit and absence fractions
- Unique-sequence count and normalized entropy

Raw unique-sequence count is not an absolute-error endpoint because it scales
with the number of sampled reads. Compare it between prefix and reservoir
samples at the same sample size and show it as a richness curve. Compare
normalized sequence entropy directly with the complete-stream value.

### Analysis

For fraction metrics, report signed error, absolute error, and 95% limits across
FASTQs. For entropy, report absolute normalized-entropy error. Compare prefix and
reservoir samples at the same sample size with paired differences. Report the
probability of observing at least one event as a function of event prevalence and
sample size.

Aggregate region instances within each FASTQ, then cluster FASTQs from the same
seqspec configuration before computing percentiles or bootstrap intervals. This
keeps a configuration with two selected read files from counting as two
independent experiments.

The default sample size will be chosen before the evaluation benchmark. The
provisional rule is to retain 10,000 reads when common fraction metrics have a
median absolute error at most 0.01, the 95th percentile error is at most 0.05,
and prefix sampling shows no material systematic shift. Metrics that fail these
criteria will use 100,000 reads or remain quantitative interpretation outputs.
The frozen endpoint list, aggregation rules, bootstrap settings, and thresholds
are defined in `experiments/paper/protocol/sampling_analysis.json`.

### Outputs

- Sampling-error table by check, metric, assay family, and read role
- Prefix-versus-reservoir paired comparison
- Runtime and memory scaling curves
- A frozen default and escalation rule for later experiments

## Experiment 2: Controlled Read-Versus-Spec Perturbations

### Question

Does seqcheck detect the inconsistency it was designed to detect, and does it
localize the result to the correct file, read, and region?

### Goal and Check

| Item | Definition |
| --- | --- |
| Measurable goal | Measure defect detection, clean-sample specificity, and file/read/region localization on known perturbations |
| Primary endpoint | Macro-averaged sensitivity, specificity, and localization accuracy across configurations |
| Expected result | Deterministic contract errors are detected and localized reliably; quantitative signals increase monotonically with the perturbed read fraction |
| Scientific target | At least 95% sensitivity for applicable deterministic perturbations, at least 95% specificity on verified unmodified controls, at least 95% correct localization, at least 90% sensitivity for applicable read-level perturbations present in 1% of reads, and a monotonic metric response in at least 90% of applicable event-fraction series |
| Completion check | Every applicable configuration-operator pair has a clean report, perturbation manifest, expected target metric, three seeded replicates where stochastic, and a classified outcome; every inapplicable pair has a recorded reason |
| If unmet | Preserve the result; inspect failures by operator and assay. Fixes made after evaluation starts require a new run and complete evaluation rerun |

### Perturbation Set

Each operator has an applicability rule. For example, a whitelist perturbation
is applied only to a region with `sequence_type: onlist`. Inapplicable cases are
recorded rather than counted as failures.

The operator definitions, stochastic fractions, seeds, expected process
outcomes, and target metrics are frozen in
`experiments/paper/protocol/perturbations.json`. Before any FASTQ or spec is
modified, `scripts/build_perturbation_cases.py` emits one content-addressed
applicability record per configuration-operator pair. An applicable record
names the exact file, read, region, coordinate interval, sequence, and resource
identity used to generate the condition. This separates a missing biological
feature from a failed detection and prevents the executor from selecting a
different target after results are seen.

Each spec mutation must also declare whether `seqspec check` is expected to pass.
The primary read-versus-spec benchmark uses internally valid-but-wrong specs.
Every generated spec is checked before seqcheck runs, and a mutation is excluded
from that primary analysis if its observed validity does not match its manifest.
S10 is an expected `seqspec check` failure because the declared local resource
is deliberately absent. It tests structured missing-resource handling and is
reported outside the internally-valid-but-wrong primary subset.

| Id | Perturbation | Target signal |
| --- | --- | --- |
| S01 | Supply an unexpected or ambiguously named FASTQ | Input matching result |
| S02 | Swap two FASTQ contents while retaining expected file names | Read-scoped structural metrics |
| S03 | Omit one expected FASTQ from the invocation | Missing expected file result |
| S04 | Change declared read length bounds | Out-of-range length fraction |
| S05 | Shift a testable region boundary by 1, 4, or 8 bases | Fixed, onlist, or coverage metric |
| S06 | Flip the declared strand of a testable region | Orientation and exact-match metrics |
| S07 | Point a read to the wrong primer anchor | Primer classification and hit location |
| S08 | Replace or reverse-complement a fixed motif | Fixed exact-match and offset metrics |
| S09 | Replace an onlist with an incompatible onlist | Exact onlist fraction |
| S10 | Make an onlist unavailable | Structured missing-resource error |
| D01 | Truncate a fraction of reads | Length and region coverage fractions |
| D02 | Corrupt a fixed motif in a fraction of reads | Fixed exact-match fraction |
| D03 | Replace onlist sequences with offlist sequences | Offlist fraction and top sequences |
| D04 | Insert, remove, or reverse a primer motif | Primer hit position and orientation |

Read-level perturbations will use event fractions of 0.01%, 0.1%, 1%, 5%, and
10% when the selected sample size can represent the event. Each stochastic
condition will use three fixed seeds. Perturbations will modify one property at a
time.

### Calibration and Frozen Evaluation

Run all applicable operators on the 12 calibration configurations. Use this set
only to fix metric thresholds, applicability rules, and binary detection rules.
Write these rules to a versioned policy file before running the 18 evaluation
configurations.

`scripts/analyze_perturbation_calibration.py` implements this boundary. It pairs
each condition with its clean report and restricts evidence to the declared
target scope. Numeric endpoints require at least six independent applicable
configurations and 80% directional consistency. Their thresholds use the larger
of a unit-specific floor and half the lower empirical decile of calibration
effects. Stochastic thresholds are fit at the predeclared 1% event fraction
after collapsing the three seeds within each configuration. Assessment rules
must appear in at least 80% of seeds within a configuration and in at least 80%
of independent configurations. These constants are frozen in
`experiments/paper/protocol/perturbation_analysis.json`.

After the policy is frozen, `scripts/sample_policy_cases.py` applies the frozen
sampling method, record count, and seed to the evaluation split. Evaluation
conditions are materialized from that content-addressed sample bundle.
`scripts/evaluate_perturbation_execution.py` then applies the frozen detection
policy without selecting a new metric, assessment code, process pattern, or
threshold. It rejects calibration samples and evaluation variants that were not
represented in the policy.

For every evaluation condition, retain the raw metric response even when a binary
threshold is not defensible. A perturbation is localized correctly only when the
reported file, read, and region include the mutated target at every applicable
scope.

### Primary Endpoints

- Sensitivity by perturbation class
- Specificity on matched unmodified samples
- Correct file, read, and region localization
- Metric effect size relative to the paired clean sample
- Detection probability as a function of event fraction and sample size

Compute macro-averages across configurations. Use configuration-level cluster
bootstrap intervals so repeated perturbations from one source do not inflate the
sample size. The confidence level, resample count, and random seed are frozen in
`experiments/paper/protocol/perturbation_analysis.json`. An incomplete endpoint
or missed target remains a valid result and does not trigger policy refitting.

### Outputs

- One perturbation manifest row per generated condition
- One tool-run row per condition
- Detection and localization matrix
- Sensitivity curves over corruption fraction
- Paired clean-versus-perturbed metric table

## Experiment 3: Complementarity With Existing QC

### Question

Do seqspec check, FastQC, and seqcheck detect different and complementary classes
of problems?

### Goal and Check

| Item | Definition |
| --- | --- |
| Measurable goal | Compare the three tools on paired structural, quality-only, and composition perturbations without treating unrelated baseline warnings as detections |
| Primary endpoint | Baseline-corrected detection rate by perturbation class and whether a tool localizes the affected seqspec entity |
| Expected result | Seqcheck has higher detection and localization for read-versus-spec defects; FastQC detects quality-only defects that leave sequences unchanged; some composition defects are visible to both tools |
| Scientific target | On internally valid-but-wrong structural defects, seqcheck exceeds FastQC by at least 25 percentage points and localizes at least 95% of its detections; `seqspec check` rejects 100% of schema-invalid controls; FastQC detects at least 90% of quality-only conditions at the frozen severity; unchanged sequences produce no seqcheck metric or assessment differences for quality-only conditions |
| Completion check | Every condition has matched clean outputs from all three tools, parsed module or assessment results, tool versions, and a predeclared expected scope |
| If unmet | Report the overlap directly and narrow the complementarity claim; do not add post hoc perturbations to improve the heatmap |

Run all three tools on the same bounded clean and perturbed samples:

| Tool | Tested contract |
| --- | --- |
| `seqspec check` | Internal validity of the seqspec document |
| FastQC | Generic sequence composition and base-quality patterns |
| `seqcheck` | Agreement between reads and the declared assay structure |

Use the structural perturbations from Experiment 2. Add four quality-only FASTQ
perturbations that change quality strings but leave names and nucleotide
sequences unchanged:

- Lower quality scores across every cycle
- Lower quality scores in selected cycles
- Lower quality scores in a controlled fraction of reads
- Create a two-component mixture of high- and low-quality reads

These are strict FastQC-positive and seqcheck-negative controls because seqcheck
does not currently consume quality strings. Also add adapter contamination,
sequence duplication, and GC-composition shifts as exploratory overlap
conditions. Those sequence-changing conditions may alter seqcheck entropy or
top-sequence metrics, so they are not negative controls. All sequence-changing
conditions must preserve read length and every declared fixed, onlist, and primer
interval. Composition changes will be restricted to eligible measurement
payloads.

Add two schema-invalid seqspec controls, one missing a required field and one
violating a region constraint. `seqspec check` should reject both. These controls
are reported separately from the internally valid-but-wrong structural set.

For each tool, a detection is a new or worsened predeclared signal relative to
the matched clean sample. Existing unrelated warnings on a clean assay do not
count as detection of a perturbation.

The primary output is a perturbation-by-tool heatmap. The intended conclusion is
complementarity. Seqcheck is not expected to detect generic base-quality defects,
and FastQC is not expected to identify the seqspec region whose coordinates were
changed.

## Experiment 4: Region Ontology Utility

### Question

Does the region ontology preserve existing behavior while making region-level QC
more comparable and more useful across assays?

### Goal and Check

| Item | Definition |
| --- | --- |
| Measurable goal | Measure upgrade parity, mapping accuracy, and the incremental QC value of ontology terms beyond `sequence_type` |
| Primary endpoint | Unexpected query changes, weighted manual mapping precision, reviewer agreement, and false-positive rate at 90% sensitivity for ontology-stratified versus `sequence_type`-only interpretation |
| Expected result | Deterministic legacy mappings preserve tool behavior, curated mappings are precise, and semantic grouping reduces false positives caused by pooling regions with different nominal roles |
| Scientific target | Zero unexpected query changes after upgrade; at least 95% weighted mapping precision; Cohen's kappa of at least 0.80; at least a 20% relative false-positive-rate reduction for model 3 versus model 2; and an upper 95% configuration-bootstrap confidence bound below zero for `FPR(model 3) - FPR(model 2)` |
| Completion check | Every eligible region has an extraction and upgrade row; all 300 validation regions have two blinded reviews and an adjudication; every model comparison uses the frozen calibration and evaluation split |
| If unmet | Preserve ambiguous labels as unknown, revise erroneous mappings only in a new ontology version, and limit the paper claim to the ontology benefits supported by the data |

### Mapping Survey

For every eligible IGVF seqspec, extract one row per region before upgrade with:

- Original `region_type`, region name, region id, and `sequence_type`
- Assay, modality, read, and parent region context
- Upgraded ontology term list
- Mapping source and whether the result is unknown

Report:

- Number and frequency of original labels
- Fraction mapped deterministically
- Fraction mapped to `RGN:unknown:unclassified`
- Number of canonical terms used
- Frequency and examples of multi-term regions
- Number of specifications whose query results change unexpectedly after upgrade

### Manual Mapping Validation

Two reviewers will independently annotate a stratified sample of 300 regions.
The sample will overrepresent rare labels, unknown mappings, and multi-role
regions, while retaining sampling weights for population estimates. Reviewers
will see the original region and assay context but not the registry mapping. At
least one reviewer must not have participated in drafting the original mapping.

Report weighted mapping accuracy, term-level precision, inter-reviewer agreement,
and adjudicated disagreement categories. Mapping an ambiguous region to unknown
is considered correct when the available context does not support a more precise
term.

### Analytical Utility

For the controlled benchmark, compare reference models that group a metric by:

1. No region semantics
2. `sequence_type` only
3. `sequence_type` plus ontology term
4. Assay family plus ontology term

Use role-appropriate metrics:

- Entropy for molecule partitions versus transcript measurements
- Onlist fraction for cell, sample, feature, and guide-related regions
- Exact-match fraction for linker, primer, adapter, and capture regions

Fit thresholds on the calibration configurations and evaluate on the locked
configurations. Compare false-positive rate at fixed sensitivity, calibration,
and area under the precision-recall curve. The primary ontology comparison is
model 3 versus model 2 because it isolates the value added beyond
`sequence_type`.

Only ontology terms represented by at least three independent configurations in
both calibration and evaluation will receive a term-specific model comparison.
Less frequent terms will remain in the mapping survey and descriptive results.

If ontology grouping does not improve held-out interpretation, report the mapping
and query benefits without claiming a QC accuracy improvement.

## Experiment 5: Current IGVF Audit

### Question

What can specification-aware QC measure at current consortium scale?

### Goal and Check

| Item | Definition |
| --- | --- |
| Measurable goal | Produce a complete, current, and internally reconciled audit of every eligible released IGVF seqspec configuration using the frozen sampling policy |
| Primary endpoint | Completion rate, categorized failure rate, sampled-record count, and distributions of raw QC metrics by assay and ontology term |
| Expected result | Completion is close to or better than the previous 96.67% audit, failures are fully classified, and the audit yields quantitative candidate inconsistencies without requiring complete FASTQ processing |
| Scientific target | At least 95% of eligible configuration-modality runs complete after frozen retries; zero unclassified runner exceptions; exact reconciliation among manifests, reports, assessments, and metrics; and identical canonical report hashes for a repeated 100-run subset |
| Completion check | Every portal object has one eligibility outcome, every eligible modality has one completed or classified-failure outcome, and all audit validation checks below pass |
| If unmet | Separate external access or transport loss from tool failure, repair tool failures under a new run identifier, and narrow the reported population to the documented eligible subset |

Create a frozen portal manifest containing every released, validated seqspec
configuration considered for the audit. Record the retrieval time, portal query,
configuration accession, linked FASTQ accessions, access status, and resource
URLs. Run in a new output directory using fixed seqspec and seqcheck commits.

Define eligibility separately for public and controlled-access data using the
credentials available at manifest freeze. Report completion by access class. Do
not retain controlled reads in the benchmark artifact bundle.

Use the sample size selected in Experiment 1, provisionally 10,000 reads per
FASTQ. Escalate selected findings to 100,000 reads according to the frozen rule.
Do not process complete files for the main audit.

Report separate denominators for:

- Portal configurations found
- Eligible configurations
- Modalities
- Expected and resolved FASTQ files
- Sampled records
- Completed reports
- Portal, transport, authentication, and tool failures
- Specification validation failures
- Read-versus-spec metrics and assessments

Summarize quantitative metrics by assay family, `sequence_type`, and ontology
term. Report warnings as candidate inconsistencies unless they meet a frozen
benchmark policy. Do not use laboratory rankings as a primary result.

The following checks must pass before analysis:

- Every completed manifest entry resolves to exactly one report.
- No report exists without a manifest entry.
- All report cache keys match the frozen run configuration.
- Run, failure, assessment, and metric table denominators reconcile.
- The raw and normalized seqspec versions are recorded.
- A rerun of a fixed subset produces identical reports apart from timestamps.

## Experiment 6: Adjudication and Downstream Consequences

### Question

Which audit findings represent confirmed specification or data problems, and do
corrections change a useful downstream result?

### Goal and Check

| Item | Definition |
| --- | --- |
| Measurable goal | Estimate confirmation precision for audit candidates and test whether selected confirmed corrections change a predeclared downstream endpoint |
| Primary endpoint | Adjudicated confirmation precision by assessment class, reviewer agreement, pass-control agreement, and paired before-versus-after endpoint change |
| Expected result | High-confidence candidates are confirmed more often than pass controls, and selected corrections improve the endpoint directly affected by the incorrect specification or resource |
| Scientific target | At least 70% confirmation precision among predeclared high-confidence candidates, Cohen's kappa of at least 0.70, at least 90% of pass controls judged consistent, at least three confirmed cases suitable for presentation, and at least two cases with a measurable downstream change in the predeclared direction |
| Completion check | At least 50 candidates and 20 blinded pass controls receive two reviews and final status; every case study has immutable before and after inputs, commands, and results |
| If unmet | Report adjudication precision without claiming confirmed real-world utility; omit downstream-rescue claims that lack a confirmed paired result |

Select at least 50 candidate findings across checks, assay families, ontology
terms, event sizes, and submitting groups. Include all rare hard errors and a
stratified sample of warnings and interpretations. Add 20 randomly selected pass
cases, matched by assay family, and blind reviewers to whether each case was a
candidate or pass control. Two reviewers will classify each case as:

- Confirmed specification problem
- Confirmed read or resource problem
- Intended assay design
- Seqcheck problem
- Inconclusive

Define the high-confidence candidate class from the frozen controlled-benchmark
policy before drawing the review sample. Report conservative confirmation
precision with inconclusive cases counted as unconfirmed and evaluable precision
with inconclusive cases shown separately.

Review packets will hide the original pass, warning, error, or interpretation
label and the high-confidence classification. Reviewers will receive the relevant
spec context, observed metric, and linked assay information needed to assess the
case.

When possible, ask the submitting group to confirm the interpretation. Preserve
the original finding, reviewer rationale, submitter response, proposed correction,
and final status in an adjudication table.

Select three to five confirmed cases for detailed presentation. At least two
cases should compare downstream results before and after correcting the seqspec or
resource. Candidate endpoints include:

- Valid cell-barcode fraction
- Reads assigned to a whitelist
- Recovered cell count
- Molecules or genes detected per cell
- Guide or perturbation assignment rate
- Usable reads after extraction or trimming

Choose case studies using a frozen rule based on defect class, assay coverage,
data access, and availability of a relevant pipeline, not on the size of the
after-correction result. Record one primary endpoint and its expected direction
before running each corrected analysis.

Complete-file processing is allowed only for these selected downstream cases.
The primary adjudication result is confirmation precision by assessment class.
Inconclusive cases remain in the denominator but are reported separately.

## Data and Output Contract

Paper runs will use a new `experiments/paper/` workspace. Generated outputs will
not be mixed with repository fixtures or the previous audit.

```text
experiments/paper/
  protocol/
    endpoints.yaml
    policy.yaml
  manifests/
    study.json
    cohort.csv
    samples.csv
    perturbations.jsonl
  reports/
    seqcheck/
    fastqc/
    seqspec_check/
  tables/
    tool_runs.csv.gz
    assessments.csv.gz
    metrics.csv.gz
    ontology_regions.csv.gz
    reviews.csv
    adjudication.csv
  validation/
    reconciliation.json
  analysis/
  figures/
```

Each derived table must be reproducible from manifests and immutable reports.
Analysis scripts must not read ad hoc files outside this layout. Large temporary
FASTQ data will live outside the repository and can be deleted after report
generation. Bounded benchmark samples will be retained in a separately versioned
artifact bundle with checksums rather than committed to Git.

## Execution Order and Gates

| Step | Measurable goal | Expected result | Gate to next step |
| --- | --- | --- | --- |
| 1. Instrumentation | Produce versioned reports, manifests, raw metric tables, deterministic samples, and cache keys | Existing metrics remain stable while provenance becomes complete | Phase 0 completion check passes with zero unresolved validation errors |
| 2. Cohort freeze | Select 30 verified, deduplicated configurations with five per assay family | The cohort spans the six families and is not selected for favorable seqcheck results | All 30 have two completed reviews, immutable source records, and a recorded 12-calibration/18-evaluation assignment |
| 3. Sampling calibration | Complete Experiment 1 on at most 24 selected FASTQs | A bounded sample gives accurate common metrics and defines explicit escalation cases | Sample-size and sampling-method rules are written and hashed before perturbation evaluation |
| 4. Detection calibration | Run every applicable perturbation on 12 calibration configurations | Operators produce their target signals and yield a usable frozen policy | Applicability matrix, severities, thresholds, and endpoint code are versioned; no evaluation result has been opened |
| 5. Controlled evaluation | Run Experiments 2 and 3 once on 18 locked configurations | Seqcheck meets structural detection targets and shows complementary coverage with FastQC | Every planned condition has reconciled outputs; results are accepted whether targets pass or fail |
| 6. Ontology evaluation | Complete the mapping survey, 300-region review, and held-out model comparison | Mapping is precise and ontology grouping improves role-specific interpretation | Review rows, adjudications, upgrade parity checks, and model outputs reconcile exactly |
| 7. Consortium audit | Audit the frozen eligible IGVF snapshot with the frozen tools and sampling rule | At least 95% of eligible runs complete and all failures are classified | Audit completion check passes, including the repeated 100-run canonical-hash check |
| 8. Adjudication | Blind-review at least 50 candidates and 20 pass controls and run selected corrections | High-confidence candidates confirm more often than controls and at least two corrections change the target endpoint | All reviews are resolved or marked inconclusive; immutable paired case-study outputs exist |
| 9. Analysis freeze | Regenerate all tables and figures from a clean study manifest | A second clean generation reproduces every canonical result | Canonical table and figure-data hashes match; software and artifact manifests are archived |

Any tool bug discovered during locked evaluation must be fixed and documented.
The affected experiment is then rerun in full under a new run identifier. Results
from before and after the fix must not be combined.

## Planned Paper Results

| Result | Primary figure or table |
| --- | --- |
| Seqcheck measurement model and check scopes | Workflow figure and check catalog |
| Sampling accuracy and computational cost | Convergence and runtime figure |
| Controlled detection and localization | Sensitivity curves and defect matrix |
| Complementarity with FastQC and seqspec check | Perturbation-by-tool heatmap |
| Ontology mapping and analytical utility | Mapping flow and stratified distributions |
| Current IGVF audit | Cohort flow and assay-stratified metric summary |
| Confirmed findings and downstream corrections | Three to five case-study panels |

The paper is ready to draft only after the controlled evaluation, clean IGVF
audit, ontology validation, and adjudication tables are frozen. The existing IGVF
audit remains useful for planning but is not a substitute for those results.
