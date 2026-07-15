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
  --cohort-split calibration \
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

Before generating Experiment 2 conditions, build the applicability inventory
from the content-addressed case selection:

```bash
python3 scripts/build_perturbation_cases.py \
  --case-selection-manifest experiments/paper/runs/<run-id>/case-selection/manifests/case_selection.json \
  --perturbation-protocol experiments/paper/protocol/perturbations.json \
  --seqspec-bin ../seqspec/target/release/seqspec \
  --output-root experiments/paper/runs/<run-id>/perturbation-inventory
```

The builder runs `seqspec check`, resolves read and region coordinates through
the seqspec CLI, and writes one row for every configuration-operator pair. Each
applicable row records the exact FASTQ, read, region, coordinates, sequence, and
resource identity needed by the perturbation executor. Each inapplicable row
records a reason. The inventory is usable only when
`validation/perturbation_cases.json` reports `valid: true` and all selected
variants are declared in `protocol/perturbations.json`.

Materialize conditions from that inventory and the retained sampling-study
matrices only after the sampling policy is frozen:

```bash
python3 scripts/materialize_perturbations.py \
  --perturbation-inventory experiments/paper/runs/<run-id>/perturbation-inventory/manifests/perturbation_inventory.json \
  --sampling-study experiments/paper/runs/<run-id>/sampling-calibration/manifests/study.json \
  --sampling-policy experiments/paper/runs/<run-id>/sampling-analysis/policy/sampling_policy.json \
  --perturbation-protocol experiments/paper/protocol/perturbations.json \
  --seqspec-bin ../seqspec/target/release/seqspec \
  --yq-bin "$(command -v yq)" \
  --output-root experiments/paper/runs/<run-id>/perturbations
```

The materializer reuses the frozen prefix or reservoir sample without another
source traversal. It hashes every generated spec, FASTQ, selected-record list,
and bundled local resource. A stochastic fraction is generated only when it
maps to an integer record count. `validation/materialization.json` must report
`valid: true` before seqcheck execution begins.

Run seqcheck once on every materialized condition and export the raw evidence:

```bash
python3 scripts/run_perturbation_calibration.py \
  --materialization-manifest experiments/paper/runs/<run-id>/perturbations/manifests/materialization.json \
  --execution-protocol experiments/paper/protocol/perturbation_execution.json \
  --seqcheck-bin target/release/seqcheck \
  --output-root experiments/paper/runs/<run-id>/perturbation-execution
```

The runner enforces the declared process outcome for each condition. Successful
conditions must write a valid JSON report; expected command failures retain
stdout and stderr but no synthetic report. Raw metrics, assessments, ontology
terms, command lines, process outcomes, runtime, and memory are written before
any detection threshold is fit. `validation/perturbation_execution.json` must
report `valid: true` before calibration analysis.

Fit the predeclared detection policy on calibration configurations only:

```bash
python3 scripts/analyze_perturbation_calibration.py \
  --execution-manifest experiments/paper/runs/<run-id>/perturbation-execution/manifests/execution.json \
  --analysis-protocol experiments/paper/protocol/perturbation_analysis.json \
  --output-root experiments/paper/runs/<run-id>/perturbation-analysis
```

The analyzer pairs each condition with its configuration's clean report and
retains only metrics at the declared file, read, and region target. Seeded
stochastic repeats are collapsed within a configuration. Stochastic policy
endpoints use the predeclared 1% anchor, while deterministic endpoints use all
applicable calibration conditions. A numeric endpoint requires at least six
independent configurations and 80% directional consistency. Its threshold is
the larger of the unit floor and half the lower empirical decile of absolute
effects. Assessment evidence must reach 80% support both across seeds within a
configuration and across configurations. Expected process failures must match a
declared stderr pattern; file localization also requires the affected input name
in stderr.

The output includes paired metric effects, endpoint candidates, condition-level
detection and localization calls, directional monotonicity checks, and a
content-addressed detection policy. `validation/perturbation_analysis.json` may
be valid while the policy remains unfrozen. The policy is usable on the locked
evaluation set only when `policy/detection_policy.json` reports `frozen: true`.
Scientific target checks in the validation file describe calibration behavior
only; final paper estimates come from the 18 evaluation configurations.

Do not inspect the evaluation split until both policies are frozen. Then select
and sample the evaluation FASTQs with the exact sampling rule chosen during
calibration:

```bash
python3 scripts/build_sampling_cases.py \
  --cohort-manifest experiments/paper/runs/<cohort-run>/freeze/manifests/cohort_frozen.json \
  --sampling-protocol experiments/paper/protocol/sampling_calibration.json \
  --seqspec-bin ../seqspec/target/release/seqspec \
  --cohort-split evaluation \
  --output-root experiments/paper/runs/<run-id>/evaluation-case-selection

python3 scripts/sample_policy_cases.py \
  --case-selection-manifest experiments/paper/runs/<run-id>/evaluation-case-selection/manifests/case_selection.json \
  --sampling-study experiments/paper/runs/<run-id>/sampling-calibration/manifests/study.json \
  --sampling-policy experiments/paper/runs/<run-id>/sampling-analysis/policy/sampling_policy.json \
  --sampler-script scripts/sample_fastq.py \
  --output-root experiments/paper/runs/<run-id>/evaluation-samples
```

Build the evaluation applicability inventory as above, but use the evaluation
case-selection manifest. Materialize from the immutable sample bundle rather
than the calibration study, then execute every condition:

```bash
python3 scripts/materialize_perturbations.py \
  --perturbation-inventory experiments/paper/runs/<run-id>/evaluation-perturbation-inventory/manifests/perturbation_inventory.json \
  --sample-bundle experiments/paper/runs/<run-id>/evaluation-samples/manifests/sample_bundle.json \
  --sampling-policy experiments/paper/runs/<run-id>/sampling-analysis/policy/sampling_policy.json \
  --perturbation-protocol experiments/paper/protocol/perturbations.json \
  --seqspec-bin ../seqspec/target/release/seqspec \
  --yq-bin "$(command -v yq)" \
  --output-root experiments/paper/runs/<run-id>/evaluation-perturbations

python3 scripts/run_perturbation_calibration.py \
  --materialization-manifest experiments/paper/runs/<run-id>/evaluation-perturbations/manifests/materialization.json \
  --execution-protocol experiments/paper/protocol/perturbation_execution.json \
  --seqcheck-bin target/release/seqcheck \
  --output-root experiments/paper/runs/<run-id>/evaluation-execution
```

Apply the frozen detection rules without fitting a new endpoint or threshold:

```bash
python3 scripts/evaluate_perturbation_execution.py \
  --execution-manifest experiments/paper/runs/<run-id>/evaluation-execution/manifests/execution.json \
  --detection-policy experiments/paper/runs/<run-id>/perturbation-analysis/policy/detection_policy.json \
  --analysis-protocol experiments/paper/protocol/perturbation_analysis.json \
  --output-root experiments/paper/runs/<run-id>/perturbation-evaluation
```

The evaluator rejects calibration samples, a policy derived from the evaluation
execution, changed policy content, and evaluation variants absent from the
frozen policy. It writes the five predeclared endpoints with deterministic
configuration-cluster bootstrap intervals. A structurally valid evaluation is
retained even when an endpoint is incomplete or a scientific target is missed.
Conditions deliberately expected to fail `seqspec check`, such as S10, remain in
the raw call and operator tables but are excluded from primary sensitivity and
localization estimates.
