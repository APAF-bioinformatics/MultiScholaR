---
tags: ["ticket", "kit", "ruv", "statistics", "all-omics"]
id: "KIT-124"
title: "RUV in every omic: convention, reference match, calibration, and one design"
description: "Establish for each of the package's RUV implementations that it is used as the field uses it, that the code matches the reference packages, and that the differential analysis after it is calibrated; then bring every omic to one design."
status: pending
priority: "P1"
dependencies: ["KIT-123"]
time_estimate: "24h"
date-created: "2026-11-16"
date-updated: "2026-11-16"
needs_planning: false
ticket_revision: "1.0"
plan: "00-KIT-PLAN.md"
files_to_modify: []
---

# KIT-124: RUV in every omic: convention, reference match, calibration, and one design

**Authority.** User, 2026-11-16: "what do you need to test across everything to make sure that the implementation of RUV is best in class and matches the current convention for how it is used", then "this sounds like this needs a comprehensive ticket". Earlier the same day: "do all that is recommended, as long as it follows DS and current proteomics best practices and conventions (feel free to research that as well or ref original paper/implementation)".

**Takes over** parked rows 10 and 11, the control-set and test-model question of KIT-123 (its R3), and KIT-114 (b). KIT-123 keeps its two bug fixes (bounded k, levels kept) and lands them first.

**Nothing here changes a default without the user's decision on measured evidence** (the gate after W3).

## Why

RUV needs three choices: which features are the negative controls, how many unwanted patterns (k) to remove, and whether the patterns are subtracted before testing or given to the test as covariates. The package has five implementations that make them differently.

| | Recommended practice | Proteomics classic (DIA with IQ, LFQ, TMT) | DIA_LIMPA | Transcriptomics | Lipid, metab |
|---|---|---|---|---|---|
| Controls | Known unchanged features, or a large empirical set | 1 to 20% least significant | 1 to 30% least significant | 5 to 30% least significant | 1 to 20% least significant |
| Replicate structure | Technical replicates | Technical replicates if declared, else the group | The same | The group | The group |
| k | Small | Elbow, no upper limit (KIT-123 R1 bounds it) | Bounded, 2 to 5 | Bounded, 3 or 4 | Elbow, no upper limit |
| Use in the test | Patterns in the model; corrected values for display | Subtract, then test | Covariates (default) or correction | Covariates | Subtract, then test |
| Level after correction | Kept | Shifted (KIT-123 R2 restores it) | Shifted in the QC matrix and in correction mode (KIT-123 R2) | Not affected | Same call as classic: expected, not measured |
| Code | | `findBestNegCtrlPercentage`, `ruvIII_C_Varying` | `findBestNegCtrlPercentageCapped`, `computeRuvCovariates`, `dpcDE` | `trxRuvGrid`, `RUVSeq::RUVs`, edgeR | `runLipidPerAssayRuvOptimization`, `func_metab_norm.R`, `ruvIII_C_Varying` |

Also in the package and to be placed in the table by W0: the peptide-level method (`func_pept_s4_norm_methods.R`).

### Measured so far (classic proteomics, one simulator; `kit-114/b_ruv_modes*.tsv`)

Achieved false discovery rate at a nominal 5%; 1,000 proteins, one unwanted factor independent of group, group as the replicate structure, 20 seeds per cell.

| Design | 20% true effects, six cells | 40% true effects, 6 per group |
|---|---|---|
| No RUV | 0.02 to 0.09 | 0.07 |
| Today: sweep's controls, subtract, then test | 0.07 to 0.48 | 0.22 to 0.30 |
| True null controls, subtract, then test | 0.09 to 0.28 | 0.14 to 0.21 |
| Sweep's controls, patterns in the model | 0.08 to 0.42 | 0.14 to 0.27 |
| True null controls, patterns in the model | 0.04 to 0.07 | 0.05 to 0.06 |
| Controls outside the 25% most significant, patterns in the model | 0.06 to 0.11 | 0.10 to 0.11 |

- **Two causes.** Subtracting before testing inflates the rate even with true controls. A small set of least-significant features is a biased control set and inflates it even with the patterns in the model.
- **A fixed 25% cut is not enough.** With 40% true effects the wide set contains true effects and the rate is 0.10 to 0.11. The cut has to follow the data.
- **Not measured:** DIA_LIMPA on its own optimiser and `dpcDE`, transcriptomics, lipid, metab. The "patterns in the model" rows used the limma kernel with the classic sweep's controls.

### References read (2026-11-16)

- RUVSeq vignette and Risso et al. 2014: normalised counts "should be used only for exploration"; the DE design "includes both the covariates of interest ... and the factors of unwanted variation"; empirical controls are "all but the top 5,000 genes"; a small control set was "problematic".
- limma `removeBatchEffect`: "not intended to be used prior to linear modelling ... include the batch factors in the linear model".
- Nygaard, Rødland, Hovig 2016; Li et al. 2023: adjusting first and then testing as if nothing had been estimated exaggerates confidence, "essentially independent of sample size".
- Poulos et al. 2020 (RUV-III-C): technical replicates across instruments, spike-in yeast peptides as controls, k chosen by an outcome check.
- RUV-III-PRPS implementation (`runRUV_III_PRPS`): W from the controls centred on their means. CRAN `ruv::RUVIII` 0.9.7.1 and RUVIIIC do not centre.
- **Not yet read:** the lipidomics and metabolomics convention (pooled-QC and internal-standard methods, hRUV).

## Scope

In: the RUV step and the test after it, in proteomics classic, DIA_LIMPA (both modes), transcriptomics, lipidomics, metabolomics and the peptide-level method.

Out: limpa's quantification and precision weights, between-sample normalisation, the q-value method (KIT-114 (a), decided), enrichment, real data before the real-data gate.

## Work, in order

### W0. Convention read and inventory (no code)

- For each implementation, one row per choice: what the code does (file and line), what the method's paper and reference code do, and whether they agree.
- Read the lipid and metab convention and state where RUV with data-chosen controls stands in it, and what the package already offers beside it (internal standards).
- Output: `docs/qa/workflow-execution/kit/kit-124/convention.md`. Exit: every cell of the table above is backed by a code reference and a source, or marked unknown.

### W1. Reference match

- The same input through the package's function and through the reference package, results equal: classic against `RUVIIIC::RUVIII_C_Varying`; DIA_LIMPA against `ruv::RUVIII`; transcriptomics against `RUVSeq::RUVs` and edgeR; lipid and metab against `RUVIIIC`.
- Tests that exist are listed, not rewritten; missing ones are added. Exit: one passing reference-match test per implementation.

### W2. Simulators and the null split (no package change)

Four simulators with known truth, each writing one row per run and arm, under `kit-124/sim/`:

| Simulator | Data | Runs the package's |
|---|---|---|
| S1 | log2 protein matrix (exists: `kit-114/sim/b_ruv_modes.R`) | classic sweep, `ruvIII_C_Varying`, limma kernel |
| S2 | precursor-level data with limpa's missingness model | DPC quantification, `findBestNegCtrlPercentageCapped`, `computeRuvCovariates`, `dpcDE`, both modes |
| S3 | negative-binomial counts | `trxRuvGrid`, RUVs, the edgeR test |
| S4 | log2 matrices of a few hundred features, several assays | the lipid and metab sweeps, `ruvIII_C_Varying`, their DA |

- **Arms in every simulator:** no RUV; today's design; today's controls with the patterns in the model; a data-following wide control set with the patterns in the model; true null controls with the patterns in the model (the ceiling); and the correct-then-test form of the last two.
- **Settings varied:** 10, 20, 40% true effects; unwanted factor absent, moderate, strong; 3, 6, 12 per group; two groups and three; with and without technical replicates; 0, 25, 50% missing where the omic has missing values.
- **The wide control set, candidates fixed here before any run:** (a) every feature outside the 25% most significant; (b) every feature whose first-pass q-value is above 0.5; (c) every feature outside the top (1 - pi0) share, pi0 estimated from the first-pass p-values.
- **Null split:** one group's replicates split at random into two false groups and the whole RUV and test sequence run. Every call is false. Run on each simulator and, at the real-data gate, on real samples (aggregates only).
- Each simulator: a pilot of 2 seeds to prove it completes, then one full run of 20 seeds per cell. One correction pass for an instrumentation fault.

### W3. Measurement and the decision table

**Rule, fixed here before any run of W3.** A design passes in a simulator if, in every cell:

1. its mean achieved FDR is no higher than the larger of the no-RUV rate and the nominal 5%, plus two standard errors of the paired difference;
2. its mean power is no lower than no-RUV's minus two standard errors;
3. it refuses or fails in no run;
4. in the null split it calls no more features than no-RUV does, plus two standard errors.

Why the anchor is "the larger of no-RUV and nominal": KIT-123's rule anchored on no-RUV alone, and in one cell no-RUV sat at 0.024, where even true controls with the patterns in the model (0.059) did not pass. Stated here so that the change of rule is visible: the six cells already run would pass for the wide set under this rule and the 40% cells would be at its edge.

Output: `kit-124/outcome.md` with, per omic, today's design and each candidate against the rule, and a recommendation: keep, change (to which design), or switch RUV off by default.

### Gate: the user's decision

One table, one row per omic. Each change of a default, of a control rule or of what the test reads is the user's. A change to DIA_LIMPA's control rule departs from the limpa script, which KIT-041 brought in "with its behaviour kept exactly": it is named as such in the table.

### W4. Implementation (only what the gate approves), one mechanism at a time

Candidate mechanisms, each with a failing-first test and its own commit:

- the control-set rule, shared by the omics that take it;
- the test reads the uncorrected values with the patterns in the model (classic proteomics, lipid, metab), the corrected values kept for figures and exports and labelled so in the GUI and the report;
- k bounded in lipid and metab, as KIT-123 R1 does for classic;
- the level kept after correction in lipid, metab and the peptide method, as KIT-123 R2 does for proteomics;
- provenance: the control rule, the number of controls, k and the mode recorded with every RUV step.

### W5. Re-record and exit

One re-record per omic touched, each difference attributed to one mechanism; journeys exact against the new records; full suite in both modes, `R CMD check`, e2e lanes; NEWS under "Scientific changes". At the real-data gate: per omic, the chosen k, the number of controls, the range of corrected values and the null split, as aggregates.

## Scientific contract

- No default, control rule or test input changes before the gate. W0 to W3 change no package code; they add tests (W1) and evidence.
- A simulator result is evidence about a design, not about a dataset. The real-data pass at the end is the check that the simulators' assumptions were not the deciding factor.
- Rules, arms and settings are those written here. A change to them after a run is a new revision of this ticket with the reason, not an edit.

## Acceptance criteria

- [ ] W0: the convention table, every cell sourced.
- [ ] W1: a passing reference-match test per implementation.
- [ ] W2: four simulators and the null split, each proven by a pilot.
- [ ] W3: every arm in every cell against the rule; the outcome with a recommendation per omic.
- [ ] Gate: the user's decision recorded per omic.
- [ ] W4, W5: each approved mechanism landed with a failing-first test; records re-made; journeys exact; no regression in the full suite; real-data aggregates at the gate.

## Constraints

R stays at 4.6.1. No real-data values are read, aggregates and hashes only. `tools/profiling/*` changes only as a harness revision with a protocol amendment. No `roxygen2::document()`. Commit locally, never push, no attribution lines.

## Rollback

W0 to W3 add files only. Each W4 mechanism is one commit and reverts alone.
