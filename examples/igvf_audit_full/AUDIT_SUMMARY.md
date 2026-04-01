# IGVF `seqcheck` audit

Date: March 30, 2026

I ran `seqcheck` on all `10,134` seqspec files on the IGVF portal against their sequencing reads. This was possible because the portal already has seqspec files attached to submitted datasets, so we can check the reads against the declared assay design directly rather than waiting for something strange to show up downstream. Each audit used one seqspec file and `10,000` reads from each linked FASTQ file. I ran the full audit with `8` workers on my MacBook M2, and the whole run finished in `9,740.97` seconds (`2h 42m 21s`).

The main result is that this kind of audit is practical at full IGVF scale. It also catches real problems. Some seqspecs are malformed, some point to missing resources, and some complete but still produce hard `seqcheck` errors. At the same time, a large fraction of the output currently lands in the interpretation layer. These are measured quantities without policy attached to them yet (e.g. how do we interpret location distribution of a fixed region?). I think the next step is to get help to decide what those policies should be for different assays and `seqspec` regions.

A summary of the audit is below:

- Number of seqspecs checked: `10,134`
- Number of FASTQ files indexed from the portal: `49,921`
- Completed audits: `9,797`
- Failed seqspecs: `337`
- Completion rate: `96.67%`

If you want to work on a spec, the main files from this audit are here in this folder. [failures.csv](failures.csv) has seqspec- and portal-level failures. [runs.csv](runs.csv) has one row per completed run. [diagnostics.csv.gz](diagnostics.csv.gz) has the row-level `seqcheck` assessments. [lab_summary.csv](lab_summary.csv) has the lab totals. [reports.tar.gz](reports.tar.gz) has the full per-run JSON reports.

## `seqspec` and portal audit

The table below summarizes the main seqspec- and portal-level issues from [failures.csv](failures.csv), plus the `24` completed runs that still had a hard `seqcheck` error.

| Outcome | Count | Percent | Description | Examples |
| --- | ---: | ---: | --- | --- |
| No linked FASTQs | 144 | `1.42%` | These seqspecs do not have linked FASTQ files, so they cannot be checked in this FASTQ-only pass. | `IGVFFI0259MQIB`, `IGVFFI0580VOOX`, `IGVFFI0632NLGY`, `IGVFFI0711ZNGR`, `IGVFFI0805GKFP`, `IGVFFI0932SZUE`, `IGVFFI1005FOWV`, `IGVFFI1020QDUF` |
| Malformed seqspec YAML | 84 | `0.83%` | These are real malformed seqspec files on the portal. | `IGVFFI0110WYBL`, `IGVFFI0226HJUA`, `IGVFFI0721PRUY`, `IGVFFI0744UOFD`, `IGVFFI0829KXQE`, `IGVFFI0974LZEY`, `IGVFFI1210FJUE`, `IGVFFI1333FBJB` |
| Missing FASTQ metadata | 24 | `0.24%` | The seqspec points to FASTQ names or paths that do not resolve through IGVF FASTQ metadata. | `IGVFFI0212ZPUD`, `IGVFFI0282ROMN`, `IGVFFI1017KMXC`, `IGVFFI1471UJSK`, `IGVFFI2418ERHU`, `IGVFFI2429KBGN`, `IGVFFI2442UPCW`, `IGVFFI2791HITN` |
| Missing or inaccessible seqspec URL | 42 | `0.41%` | `seqspec version`, `info`, or `file` could not read the remote seqspec URL. This combines `21` version failures, `14` info failures, and `7` file failures. | `IGVFFI0464XYHE`, `IGVFFI3091GPCZ`, `IGVFFI3096OJLP`, `IGVFFI3102WZNV`, `IGVFFI3104GUKT`, `IGVFFI3105FROG`, `IGVFFI3105NWFB`, `IGVFFI3107CPXS` |
| Missing onlist resource | 24 | `0.24%` | The run completed, but `seqcheck` could not load an onlist resource referenced by the seqspec. | `IGVFFI0459TBIF`, `IGVFFI0933UFCY`, `IGVFFI1308UNYG`, `IGVFFI2644YTDQ`, `IGVFFI3212WSNE`, `IGVFFI3463NWOK`, `IGVFFI3733NIWD`, `IGVFFI4331AVET` |
| Transport failures while reading seqspec URLs | 27 | `0.27%` | These were connection reset failures while reading remote seqspec files. This looks like tool/runtime follow-up rather than bad portal content. | `IGVFFI0055XEHV`, `IGVFFI1163NTEB`, `IGVFFI1742ECFU`, `IGVFFI2071AEQQ`, `IGVFFI3100MHYG`, `IGVFFI3105GAIN`, `IGVFFI3105WHXP`, `IGVFFI3113BANQ` |
| Runner exceptions | 16 | `0.16%` | All `16` had the same message: `Expecting value: line 1 column 1 (char 0)`. | `IGVFFI1426OVDM`, `IGVFFI2574ZRKS`, `IGVFFI2982KSVM`, `IGVFFI3146GSKQ`, `IGVFFI4106VMGU`, `IGVFFI4661RHDG`, `IGVFFI5208SSVS`, `IGVFFI7257YBUF` |

The clearest portal-side problems from this pass are:

- `84` malformed seqspec YAMLs
- `24` seqspecs with FASTQ references that do not resolve through IGVF metadata
- `24` completed runs with missing onlist resources

## `seqcheck` audit

`seqspec check` asks whether the seqspec is valid as a file. `seqcheck` asks whether the reads are consistent with the seqspec they are supposed to describe. In that sense, this audit is a check on the seqspec with respect to the reads. To see the checks run in this audit please review the `seqcheck` [README](../../README.md#check-outcomes) under `Check Outcomes`. 

In total there were `246,697` assessments (e.g. `Pass`, `Warning`, `Interpretation`, `Error`) across the completed runs.

| Outcome | Count | Percent | Description | Examples |
| --- | ---: | ---: | --- | --- |
| `Pass` | `154,951` | `62.81%` | The seqspec and the sampled reads agree for that check. | `IGVFFI0001OUZN` contributes `16` pass assessments. |
| `Warning` | `39,671` | `16.08%` | The check found something that is still usable, but likely worth review. | `IGVFFI0001OUZN` contributes `6` warnings, including `offlist_sequences_detected`. |
| `Interpretation` | `52,051` | `21.10%` | The check is describing a pattern in the reads or in the seqspec that may or may not reflect a problem. | `IGVFFI0001OUZN` contributes `10` interpretations, including `shared_region_type_visible_in_multiple_reads`. |
| `Error` | `24` | `0.01%` | The check found something that should be treated as a hard failure. | `IGVFFI0459TBIF` contributes one error, `missing_onlist_resource`. |

Runs with at least one error: `24`

All `24` errors had the same code: `missing_onlist_resource`.

## Defining QC policy

One important thing that came out of this audit is that the `Interpretation` layer is helpful but underspecified. Many of the interpretation outputs are measured quantities that could become warnings or errors once we decide what the right rules are. For example, `IGVFFI0001OUZN` has three `random_sequence_entropy` interpretations:

- `umi` at `0.5258` of the theoretical maximum
- `poly_A` at `0.1555`
- `cdna` at `0.0442`

Those are real numbers computed from the reads and right now they are reported, but there is no assay-specific rule that says when an entropy value is acceptable, when it should raise a warning, and when it should be treated as a hard failure. The answer is likely different for a UMI, a polyA segment, and a cDNA payload, so this cannot be one universal cutoff.

The same accession also has three `fixed_exact_match_partial` interpretations:

- `rna-linker1` exact match fraction `0.9906`
- `rna-linker2` exact match fraction `0.9584`
- `rna-linker3` exact match fraction `0.9146`

Again, these are real and useful quantities, but we have not yet said what the thresholds should be. For a given fixed region we need a policy that says what exact-match fraction is expected, what fraction is suspicious, and what fraction means the read structure or chemistry is likely wrong.

The same file also has a `shared_region_type_visible_in_multiple_reads` interpretation for `barcode`. In this case the barcode shows up multiple reads (like the CRISPR assays). That could likely reflect intended overlap but how that overlap is supposed to be handled downstream still remains underspecified. So here too, we need a policy that says when overlap is expected and how downstream tooling should handle it (e.g. `seqspec index --no-overlap`)

More generally, when I say we need policy rules, I mean we need to define things like:

- which checks should remain descriptive and which should become actionable
- assay-specific cutoffs for entropy, exact-match fractions, offlist fractions, and related metrics
- whether thresholds should depend on region type, for example `umi` versus `cdna` versus `barcode`
- whether some patterns should be warnings in one assay and errors in another
- what metadata should accompany a failure so a lab knows exactly what to fix

I also think there are likely other informative checks that are not in `seqcheck` yet and should be added. There may be additional diagnostics around overlap handling, barcode structure, guide structure, orientation consistency, length distributions, or other assay-specific properties that would also be useful to compute and benchmark.

It would be helpful to have labs work on describing the kinds of checks they think are useful to be made in their reads so we can define what those policy rules should look like, and which additional checks are worth adding. In practice that means deciding on cutoffs, thresholds, assay-specific expectations, and what should count as a true failure rather than just an interpretation.

## Lab breakdown

This table comes from [lab_summary.csv](lab_summary.csv).

### How to read this table

- `# completed seqspecs` is the number of completed audit runs for that lab.
- `# failed seqspecs` is the number of unique seqspec accessions from that lab that failed in any way. This includes seqspecs that did not produce a completed run row, and seqspecs that completed but still had a hard `seqcheck` error such as `missing_onlist_resource`.
- `seqcheck passes`, `seqcheck warnings`, `seqcheck errors`, and `seqcheck interpretations` are totals summed across that lab's completed runs.
- `seqcheck errors` here is the number of error assessments across completed runs, not the number of failed seqspecs.

Examples:

- `IGVFFI0001OUZN` contributes one completed run to `Jason Buenrostro, Broad`, with `16` pass, `6` warning, and `10` interpretation assessments.
- `IGVFFI0110WYBL` contributes one failed seqspec to `Jay Shendure, UW` because the YAML is malformed.
- `IGVFFI0459TBIF` contributes one failed seqspec to `Charles Gersbach, Duke` because it completed with the hard error `missing_onlist_resource`.

| Lab | # completed seqspecs | # failed seqspecs | seqcheck passes | seqcheck warnings | seqcheck errors | seqcheck interpretations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `(unknown)` | 0 | 16 | 0 | 0 | 0 | 0 |
| `Ali Mortazavi, UCI` | 1940 | 14 | 27296 | 9673 | 0 | 6219 |
| `Ansuman Satpathy, Stanford` | 619 | 7 | 10796 | 1850 | 0 | 2164 |
| `Charles Gersbach, Duke` | 210 | 56 | 2458 | 434 | 24 | 678 |
| `Gary Hon, UT` | 233 | 1 | 3548 | 1871 | 0 | 1564 |
| `Hao Wu, UPenn` | 50 | 2 | 400 | 100 | 0 | 150 |
| `Harinder Singh, University of Pittsburgh` | 53 | 1 | 1026 | 369 | 0 | 286 |
| `Jason Buenrostro, Broad` | 2552 | 25 | 36974 | 14733 | 0 | 23786 |
| `Jay Shendure, UW` | 139 | 78 | 1903 | 196 | 0 | 698 |
| `Jesse Engreitz, Stanford` | 116 | 0 | 2339 | 844 | 0 | 787 |
| `Jimmie Ye, UCSF` | 569 | 2 | 10469 | 1365 | 0 | 2845 |
| `Kathrin Plath, UCLA` | 10 | 0 | 120 | 10 | 0 | 80 |
| `Leif Ludwig, MDC` | 298 | 2 | 3393 | 499 | 0 | 843 |
| `Maya Kasowski, Stanford` | 32 | 144 | 508 | 140 | 0 | 184 |
| `Nadav Ahituv, UCSF` | 1 | 0 | 13 | 2 | 0 | 2 |
| `Ryan Corces, Gladstone Institute UCSF` | 2907 | 13 | 52637 | 7220 | 0 | 11520 |
| `Thomas Quertermous, Stanford` | 68 | 0 | 1071 | 365 | 0 | 245 |

## Main files

- Run catalog: [runs.csv](runs.csv)
- Failure catalog: [failures.csv](failures.csv)
- Diagnostic catalog: [diagnostics.csv.gz](diagnostics.csv.gz)
- Lab summary: [lab_summary.csv](lab_summary.csv)
- Report archive: [reports.tar.gz](reports.tar.gz)

This analysis was performed with Codex.
