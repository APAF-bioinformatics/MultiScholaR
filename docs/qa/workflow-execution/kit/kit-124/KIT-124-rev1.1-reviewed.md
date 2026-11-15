---
tags: ["ticket", "kit", "ruv", "statistics", "all-omics"]
id: "KIT-124"
title: "RUV in every omic: convention, reference match and calibration, up to the user's decision"
description: "For each of the package's RUV implementations: is it used as the field uses it, does the code match the reference packages, and is the test after it calibrated against the field's alternatives. Evidence and a decision table only; no package behaviour changes."
status: in_progress
priority: "P1"
dependencies: ["KIT-123"]
time_estimate: "40h"
date-created: "2026-11-16"
date-updated: "2026-11-16"
needs_planning: false
ticket_revision: "1.1"
plan: "00-KIT-PLAN.md"
files_to_modify: []
---

# KIT-124: RUV in every omic: convention, reference match and calibration, up to the user's decision

**Authority.** User, 2026-11-16: "what do you need to test across everything to make sure that the implementation of RUV is best in class and matches the current convention for how it is used"; "this sounds like this needs a comprehensive ticket"; after the audit, on applying its repairs and on comparators: "those decisions should follow best practices and best practices only".

**Revision 1.1 (2026-11-16).** Revision 1.0 was audited the same day (`docs/qa/workflow-execution/kit/kit-124/audit-20261116.md`; the audited text is kept as `KIT-124-rev1.0-audited.md`). This revision applies every finding; the remediation record is `audit-20261116-remediation.md`. No run of this ticket preceded it.

**Started 2026-11-16** (user: "go."): W0 and W1. The simulations wait for the independent read.

## Bounded outcome

One decision table for the user, one row per omic, backed by a convention table, reference-match tests and a pre-specified simulation study. **This ticket changes no package behaviour.** What the user approves at the gate is implemented in a later ticket, authored then.

Takes over parked rows 10 and 11, KIT-123's former R3 and KIT-114 (b).

## Standards this ticket is held to

Each is re-read at its source in W0 before it is relied on.

- **Simulation studies:** planned and reported as aims, data-generating mechanisms, targets, methods and performance measures, with Monte Carlo standard errors and the number of repetitions set by the precision needed (Morris, White, Crowther 2019, Statistics in Medicine).
- **Benchmarking:** the field's current methods are in the comparison; measures are fixed in advance; simulated and real data are both used (Weber et al. 2019, Genome Biology).
- **Error-rate robustness:** an empirical rate within 0.5 to 1.5 times the nominal level (Bradley 1978, the liberal criterion).
- **Pre-specification:** rules, cells, arms and seeds are fixed here before any run. Runs made before this revision are exploratory.

## The five implementations today

| | Reference usage (source per cell fixed in W0) | Proteomics classic (DIA with IQ, LFQ, TMT) | DIA_LIMPA | Transcriptomics | Lipid, metab |
|---|---|---|---|---|---|
| Controls | Known unchanged features (RUV-III-C paper: spike-ins); RUVSeq's example: all but the most significant quarter | 1 to 20% least significant | 1 to 30% least significant | 5 to 30% least significant | 1 to 20% least significant |
| Replicate structure | Technical replicates | Technical replicates if declared, else the group | The same | The group | The group |
| k | No rule; small, by diagnostics | Elbow, no upper limit (KIT-123 R1 bounds it) | Bounded, 2 to 5 | Bounded, 3 or 4 | Elbow, no upper limit |
| Use in the test | Patterns in the model; corrected values for display (RUVSeq, limma) | Subtract, then test | Covariates (default) or correction | Covariates | Subtract, then test |
| Level after correction | Kept | Shifted (KIT-123 R2) | Shifted in the QC matrix and in correction mode (KIT-123 R2) | Kept (limma centres its covariates; checked) | Same call as classic: expected, not measured |
| Code | | `findBestNegCtrlPercentage`, `ruvIII_C_Varying` | `findBestNegCtrlPercentageCapped`, `computeRuvCovariates`, `dpcDE` | `trxRuvGrid`, `RUVSeq::RUVs`, edgeR | `runLipidPerAssayRuvOptimization`, `func_metab_norm.R`, `ruvIII_C_Varying` |

Every omic also chooses its controls and its k by the same objective: the largest gap between "all features" and "controls" in canonical correlation with the group.

### Exploratory evidence (2026-11-16; classic kernels, one simulator, 20 seeds; not confirming)

`kit-114/b_ruv_modes*.tsv`. Achieved false discovery rate at a nominal 5%.

| Design | Unwanted factor, 20% true effects | The same, 40% | No unwanted factor |
|---|---|---|---|
| No RUV | 0.02 to 0.09 | 0.07 | 0.045 to 0.050 |
| Today: sweep's controls, subtract, then test | 0.07 to 0.48 | 0.22 to 0.30 | 0.065 to 0.084 |
| True null controls, subtract, then test | 0.09 to 0.28 | 0.14 to 0.21 | 0.09 to 0.17 |
| Sweep's controls, patterns in the model | 0.08 to 0.42 | 0.14 to 0.27 | 0.045 to 0.055 |
| True null controls, patterns in the model | 0.04 to 0.07 | 0.05 to 0.06 | 0.050 to 0.057 |
| Controls outside the 25% most significant, patterns in the model | 0.06 to 0.11 | 0.10 to 0.11 | 0.049 to 0.057 |

Read as hypotheses for the study below: subtracting before testing inflates the rate by itself; a control set chosen with the result is biased when an unwanted factor exists; RUV costs power when there is nothing to remove (0.72 against 0.90 for today's design at 50% missing).

## Scope

In: the RUV step and the test after it in proteomics classic, DIA_LIMPA (both modes), transcriptomics, lipidomics and metabolomics. The peptide-level method is placed or scoped out by W0.

Out: limpa's quantification and precision weights, between-sample normalisation, the q-value method (KIT-114 (a), decided), enrichment, real data before the real-data gate, and every change to package behaviour.

## W0. Convention, inventory and consumers (reading only)

1. **Convention table.** For each implementation and each choice: what the code does (file and line), what the method's paper says in its Methods and what its authors' own code does, and whether they agree. Each cell of "Reference usage" cites its source and says "example" where it is one. The references taken through summaries on 2026-11-16 (the RUV-III-C paper's Methods; the third-party `runRUV_III_PRPS`) and the four standards above are re-read at source.
2. **Lipid and metab convention:** pooled-QC, internal-standard and replicate-based methods (hRUV among them), and where RUV with data-chosen controls stands.
3. **Inventory:** per omic, whether a known batch can be put in the test's model; which designs carry technical replicates or pooled QC samples; which route uses the peptide-level method.
4. **Consumers.** One table: each producer of a corrected matrix or of the W design, what it writes, and every reader (figures, heatmaps, exports, report, enrichment's input, integration workbook, session export, study parameters).

Output: `docs/qa/workflow-execution/kit/kit-124/convention.md`. Exit: no cell unsourced or unmarked.

## W1. Reference match (tests only)

The same input through the package's function and through the reference package, results equal: classic and lipid/metab against `RUVIIIC::RUVIII_C_Varying`; DIA_LIMPA against `ruv::RUVIII`; transcriptomics against `RUVSeq::RUVs` and edgeR. Existing tests are listed, missing ones added. Exit: one passing reference-match test per implementation.

## Simulation Setup

### Aims

1. Is the test after each RUV design calibrated?
2. Does the design gain over no RUV where unwanted variation exists, and do no harm where it does not?
3. How does it stand against the field's alternatives?

### Data-generating mechanisms

| Simulator | Data | Runs the package's |
|---|---|---|
| S1 | log2 protein matrix (from `kit-114/sim/b_ruv_modes.R`) | classic sweep, `ruvIII_C_Varying`, limma kernel |
| S2 | precursor-level data from `limpa::simProteinDataSet` with an unwanted factor added | DPC quantification, `findBestNegCtrlPercentageCapped`, `computeRuvCovariates`, `dpcDE`, both modes |
| S3 | negative-binomial counts | `trxRuvGrid`, RUVs, the edgeR test |
| S4 | lipid and metab matrices, several assays | their sweeps, `ruvIII_C_Varying`, their DA |

**S4's model and comparators are fixed after W0, in a revision of this ticket, before S4 runs.** It includes drift along injection order.

**Primary cells: a centre and one setting moved at a time (16 per simulator).**

| Setting | Centre | Moved to |
|---|---|---|
| Per group | 6 | 3; 12 |
| True effects | 20% | 10%; 40% |
| Unwanted factor (sd of loadings) | 0.5 | none; 1.0 |
| Its correlation with group | 0 | 0.3; 0.6 |
| Number of factors | 1 | 3 |
| Groups | 2 | 3 |
| Technical replicates | none | duplicates |
| Missing (where the omic has it) | 25% | 0; 50% |
| Noise | normal | t with 4 degrees of freedom |
| Design | groups only | a known, unbalanced batch |

### Targets

The log2 group difference of each feature and whether it is zero.

### Methods (arms)

| # | Arm | Purpose |
|---|---|---|
| 1 | No RUV | baseline |
| 2 | The true factor in the model | ceiling |
| 3 | Today's design, with KIT-123's fixes | status quo |
| 4 | Today's controls, patterns in the model | isolates the test form |
| 5 to 7 | Wide first-pass controls, patterns in the model: outside the 25% most significant; first-pass q above 0.5; outside the top (1 - pi0) share | outcome-based wide sets |
| 8 | All features as controls, patterns in the model | no selection by the result |
| 9 | 50 known nulls as controls, patterns in the model (stands for a supplied list: spike-ins, housekeeping, internal standards) | the papers' usage |
| 10 | Arm 9, subtract then test | isolates the test form with good controls |
| 11 to 13 | Arm 8 with k fixed at 1; at 2; estimated and allowed to be 0 | k without the group; "no RUV needed" as a result |
| 14, 15 | `ruv::RUV4` with today's controls; with arm 9's | the RUV authors' own test |
| 16 | `ruv::RUVinv` with arm 9's controls | the same, no k |
| 17 | limma with surrogate variables (`sva`) | the most used alternative |
| 18 | The known batch in the model (batch cell only) | the standard for a known batch |

- Where a cell has technical replicates, arms 3 to 13 run with them as the replicate structure and again with the group.
- S2 runs arms 3, 4 and 10 as DIA_LIMPA's correction and covariates modes with `dpcDE`. S3 uses RUVs, the vignette's RUVg design and `svaseq` for arms 3, 14 to 17.
- `sva` is installed in a scratch library that is used by the simulators only; the pinned libraries are not touched. Its version is recorded.

### Performance measures (each with its Monte Carlo standard error)

- share of true nulls with p under 0.05, and under 0.01;
- mean false discovery proportion among features the package would call (its q-value under 0.05);
- share of runs with a false discovery proportion above 0.10;
- power;
- bias and root mean square error of the log2 differences of the true effects;
- share of runs that refuse or fail.

### Repetitions and seeds

- A pilot of 50 runs at the centre cell of each simulator gives the between-run standard deviations.
- Repetitions per primary cell: the smallest number that brings the Monte Carlo standard error to 0.005 for the null rate and 0.0075 for the paired difference in false discovery proportion; at least 200, at most 1,000.
- Seeds start at 100001. The exploratory runs used 9001 to 10300.
- Every result file has a sidecar: commit, library path, R and package versions, command, seed range.

### Null split

One group's biological units are split at random into two false groups, technical replicates kept together, at least 3 units per false group; the whole RUV and test sequence is run. Measures: the share of splits with any call, and the share of features with p under 0.05. Run in each simulator at the centre cell, and once end to end through the worker operations on a synthetic corpus. Where a group is too small, group labels are permuted instead and the result is marked as the weaker check.

### Compute

At most 60 core-hours per simulator. If the repetitions needed do not fit, arms are dropped from the end of this order: 7, 6, 12, 15, and the cells left at 200 repetitions are reported as such.

### Independent read

Before W3 runs, this ticket and the simulators are read by a reviewer that did not write them. Its findings are answered in a revision.

## W3. Decision rule (fixed here)

For a design in a primary cell, with one-sided 97.5% Monte Carlo bounds:

1. **Calibration:** the null rate at 0.05 is at most 0.075.
2. **False discoveries:** its false discovery proportion exceeds no-RUV's, run for run, by at most 0.02.
3. **Benefit, where an unwanted factor exists:** its power is not below no-RUV's by more than 0.02.
4. **No harm, where none exists:** its power is not below no-RUV's by more than 0.02, and rule 1 holds.
5. **Robustness:** it refuses or fails in at most 1% of runs.
6. **Null split:** the share of splits with any call is at most 0.075.

- **Pass:** the bound is inside the limit in every primary cell.
- **Fail:** in some cell the whole interval lies outside a limit.
- **Inconclusive:** an interval straddles a limit at 1,000 repetitions. It is reported as that and is not rerun.

A design is recommended for an omic only if it passes. Among passing designs the one closest to the reference usage is preferred, then the one with the higher power. If none passes for an omic, the recommendation is RUV off by default there.

Output: `kit-124/outcome.md`: per omic, each arm against each rule in each primary cell, the null split, and the recommendation.

## Gate: the user's decision

One table, one row per omic: keep, change (to which design) or off by default. A change to DIA_LIMPA's control rule departs from the limpa script, which KIT-041 brought in "with its behaviour kept exactly"; the table names it as such. Implementation is authored after this gate, as its own ticket, from the consumer table of W0.

## Acceptance Criteria

- [ ] W0: the convention, inventory and consumer tables, every cell sourced.
- [ ] W1: a passing reference-match test per implementation.
- [ ] S1 to S3 and the null split proven by pilots; S4 specified in a revision after W0, then proven.
- [ ] The independent read answered.
- [ ] W3: every arm against every rule in every primary cell, with Monte Carlo standard errors and sidecars; the outcome with a recommendation per omic.
- [ ] The user's decision recorded per omic.

## Regression Checks

- [ ] `git diff` over `R/` is empty for this ticket.
- [ ] The W1 test files and their dependents pass, one file per process; `R CMD check` at its baseline.

## Backward Compatibility

- [ ] No behaviour, default, stored state or output changes.

## Rollback Note

- [ ] The W1 test commits revert alone; evidence files are additive.

## Constraints

R stays at 4.6.1 and the pinned libraries are not installed into. No real-data values are read, aggregates and hashes only. No file under `tools/profiling/` changes. No `roxygen2::document()`. Commit locally, never push, no attribution lines.
