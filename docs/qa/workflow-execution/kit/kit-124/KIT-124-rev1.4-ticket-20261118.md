---
tags: ["ticket", "kit", "ruv", "statistics", "all-omics"]
id: "KIT-124"
title: "RUV in every omic: convention, reference match and calibration, up to the user's decision"
description: "For each of the package's RUV implementations: convention, reference match and calibration against the field's alternatives. Evidence and a decision table only; no package behaviour changes."
status: in_progress
priority: "P1"
dependencies: ["KIT-123"]
time_estimate: "40h"
date-created: "2026-11-16"
date-updated: "2026-11-16"
needs_planning: false
ticket_revision: "1.4"
plan: "00-KIT-PLAN.md"
files_to_modify: []
---

# KIT-124: RUV in every omic: convention, reference match and calibration, up to the user's decision

**Authority.** User, 2026-11-16: the validation question and "this sounds like this needs a comprehensive ticket"; on open choices, "those decisions should follow best practices and best practices only"; "go."; "keep continuing agentically till this is done."

**Revisions.** 1.0 audited (`kit-124/audit-20261116.md`); 1.1 applied the audit and was read independently (`independent-review-20261116.md`: not ready to run). **1.2 (this text)** applies that review (R-01 to R-20) and the convention read (W0). Earlier texts are kept in `docs/qa/workflow-execution/kit/kit-124/`. No confirming run preceded this revision.

## Bounded outcome

One decision table for the user, one row per omic. **No package behaviour changes here.** Implementation is authored after the gate.

## Done

- **W0** (`kit-124/convention.md`, `c5576693`): the authority for what each implementation does and what its sources say. Its main results:
  - No source supports choosing a narrow set of least significant features sized by the canonical-correlation gap. Data-chosen controls in the literature are wide and fixed in advance (RUVSeq: all but the top 5,000 genes; RUV-4: q above 0.5); the RUV-III-C paper used peptides known not to change and peptides seen in every sample.
  - **The sources disagree on the test.** The RUV-III and RUV-III-C papers test on corrected values, with technical replicates. RUV-2, RUV-4, RUVSeq and limma put the factors in the model. Both are reference usages; neither is "the" convention.
  - The reference packages `ruv::RUVIII` and `RUVIIIC` shift feature levels; the RUV-III-PRPS authors' code centres. Lipid, metab, the peptide method and the transcriptomics RUVIII export still shift.
  - k: every source asks for judgement; published values run from 1 to 20.
  - Lipid and metab: the field's default is pooled-QC regression on injection order and internal standards per class; no guideline read recommends RUV with data-chosen controls and the group as replicates.
- **W1** (`e5ecd785` on branch `kit-123`): each implementation equals its reference package (proteomics up to the constant per protein that KIT-123 R2 removes).

## Standards

Morris, White, Crowther 2019 (aims, mechanisms, targets, methods, measures, Monte Carlo error; failures are not random). Weber et al. 2019 (current methods in the comparison; simulations shown to reflect real data). Bradley 1978, through secondary sources: the liberal band 0.025 to 0.075 at 0.05. Nygaard et al. 2016; Li et al. 2023. Efron 2007 and Leek and Storey 2008 on dependence and the run-to-run spread of the false discovery proportion.

**Claim limit (R-18).** The simulators generate RUV's own linear model. A fail under it is strong evidence against a design; a pass is weak evidence for it. The decision table says so, and the real-data gate later compares the simulators' aggregates (missing share, variance, factor strength) with real ones.

## Simulation Setup

### Simulators

Each has a parameter table in its header and, with its pilot, is committed before any confirming run.

| | Data | Runs the package's |
|---|---|---|
| S1 | log2 protein matrix, 1,000 proteins | classic sweep, `ruvIII_C_Varying`, the limma kernel |
| S2 | `limpa::simProteinDataSet`, 1,000 proteins; the factor added to the complete precursor matrix, missingness then redrawn with the returned DPC coefficients (R-11) | DPC quantification, `findBestNegCtrlPercentageCapped`, `computeRuvCovariates`, `dpcDE`, both modes |
| S3 | negative-binomial counts, 5,000 genes | `trxRuvGrid`, RUVs, the edgeR test |
| S4 | lipid and metab matrices, 300 features, with pooled-QC injections, drift along injection order and a batch step | their sweeps, `ruvIII_C_Varying`, their DA, the package's internal-standard normalisation |

Common to all: rows permuted (R-19); effects of 1 log2 with sign as the cell says; the factor's loadings normal, its strength stated as the ratio of its variance to the median residual variance (centre: 3); the factor's correlation with group fixed by construction; the test's settings those of each omic's default path, pinned by a test (R-16); seed = 100000 + 10000 x cell + run, one stream per arm.

### Cells

| | Cell | |
|---|---|---|
| C0 | centre: 6 per group, 20% effects, one factor, uncorrelated, 25% missing where the omic has it, 2 groups | core |
| C1, C2 | 12 per group; 10% effects | core |
| C3 | no unwanted factor | core |
| C4 | no true effects (global null), with the factor | core |
| C5 | factor correlated with group, 0.3 | core |
| C6 | technical duplicates (4 units per group, 2 injections each) | core |
| C7 | 3 groups | core |
| S1 to S10 | 3 per group; 40% effects; correlation 0.6; 3 factors; 0% and 50% missing; heavy tails; a known unbalanced batch; 80% of effects up; loadings tied to abundance | stress |
| S11 to S14 | 3 per group with correlation 0.3; 3 per group with 3 factors; no factor with 3 per group; no factor with 50% missing | stress |

**1,000 runs per core cell and 500 per stress cell, for every arm** (no pilot-dependent number, R-08). A pilot of 20 runs per cell proves completion and gives the cost. Ceiling: 150 core-hours per simulator; if a simulator's cost exceeds it, stress cells drop to 250 runs and that is reported.

### Arms

| # | Arm | k | Test on |
|---|---|---|---|
| 1 | No RUV | | |
| 2 | The true factor in the model (ceiling) | | uncorrected |
| 3 | Today's design (with KIT-123's fixes) | today's rule | as today |
| 4 | Today's controls and k, factors in the model | today's | uncorrected |
| 5 | Controls outside the 25% most significant, factors in the model | today's | uncorrected |
| 6 | Controls with first-pass q above 0.5, factors in the model | today's | uncorrected |
| 7 | All features as controls, factors in the model | E | uncorrected |
| 8 | Oracle controls (50 true nulls), factors in the model | E | uncorrected |
| 9 | Degraded supplied set (50 controls, 10 of them true effects), factors in the model | E | uncorrected |
| 10 | Arm 8, subtract then test | E | corrected |
| 11 | Arm 7 with k fixed at 1 | 1 | uncorrected |
| 12 | `ruv::RUV4`, oracle controls, `p.ebayes` | E | |
| 13 | `ruv::RUVinv`, oracle controls, `p.ebayes` | none | |
| 14 | limma with surrogate variables (`sva`) | E | uncorrected |
| 15 | The papers' RUV-III usage: technical replicates as the replicate structure, all features as controls, subtract then test (cell C6 only) | E | corrected |
| 16 | Arm 15 with the factors in the model (cell C6 only) | E | uncorrected |
| 17 | The known batch in the model (batch cell only) | | uncorrected |

- **E** is one estimator, the same for every arm that uses it, and it may return 0 (then the arm is no RUV): `sva::num.sv(method = "be", B = 20)` on the data with the group model, fixed seed (R-05).
- Arms that need complete data (12 to 14) are scored on the features complete in every sample; every other arm is scored on the full set and again on that subset (R-13).
- S2 runs arms 3 and 4 as DIA_LIMPA's correction and covariates modes. S3 uses RUVs for the package arms, RUVg for arm 5, RUVr as an added arm and `svaseq` for arm 14. S4 adds: internal-standard normalisation alone; pooled-QC injections as the replicate structure with all features as controls; LOESS on injection order fitted to the pooled QCs (the field's default, computed in the simulator).
- `sva` 3.60.0 is in a scratch library (`kit/scratch-k124-lib`) used by the simulators only.

### Measures (each with its Monte Carlo standard error over runs)

Share of true nulls with p under 0.05 and under 0.01; mean false discovery proportion among the package's calls (q under 0.05); share of runs with a false discovery proportion above 0.10; power; bias and root mean square error of the log2 differences of true effects; designed refusals and crashes, counted apart. A refused run is scored as its no-RUV result (R-07).

### Null check (R-09)

Cell C4 is the global null at full group size. Measure: the share of runs with any call, beside arm 1's. A split of one group's units into false groups is run once end to end through the worker operations on a synthetic corpus, as a smoke test only. Neither tests control selection in the presence of true effects, nor a factor tied to group; cells C5 and S3 do.

## W3. Decision rule (fixed here; applied to every arm, arm 1 included)

One-sided 97.5% bounds; exact binomial bounds for shares.

| # | Rule | Limit |
|---|---|---|
| 1 | Null rate at 0.05, and at 0.01 | upper bound at most 0.075, and 0.015 |
| 2 | Mean false discovery proportion | upper bound at most 0.075 |
| 3 | Share of runs with false discovery proportion above 0.10 | upper bound at most the true-factor arm's share plus 0.05 |
| 4 | Gain, in cells with a factor | lower bound of the power gain over no RUV above 0, or of the reduction in rule 3's share above 0 |
| 5 | No harm, in cells without a factor | upper bound of the power loss against no RUV at most 0.02 |
| 6 | Crashes | upper bound at most 0.01 |
| 7 | Global null (C4): share of runs with any call | upper bound at most 0.075 |

- **Rule 5 is a stated judgement.** Fitting two extra factors at 6 per group costs about 0.02 power by arithmetic, so a design that always removes two cannot pass it; one that can return k = 0 can. The measured losses are in the table for the user to overrule.
- A cell where the true-factor arm fails rules 1 or 2 is void and reported.
- **Pass:** every rule in every core cell. **Fail:** some bound in a core cell outside its limit with the whole interval outside. **Inconclusive:** otherwise. Stress cells never fail a design: they give its conditions of use (for example a minimum group size).

**Recommendation per omic.**
- A design that passes is eligible. Among eligible designs the order of preference is: the method paper's own usage where the data have what it needs (arms 15, 16 with technical replicates; arm 8 or 9 where controls are supplied); then all features as controls with the factors in the model (arm 7); then the others, by the lower share in rule 3. Power is not a tie-break (R-15).
- If no RUV design passes and arm 1 passes rules 1, 2, 6 and 7: RUV off by default.
- If neither passes: reported as such, with the cells and the measures.

Output: `kit-124/outcome.md`.

## Gate: the user's decision

One row per omic: keep, change (to which design, under which conditions) or off by default. A change to DIA_LIMPA's control rule departs from the limpa script; the table names it.

## Acceptance Criteria

- [x] W0, W1.
- [ ] S1 to S4 committed with parameter tables and pilots.
- [ ] Every arm against every rule in every cell, with standard errors and sidecars (commit, libraries, versions, command, seeds).
- [ ] `outcome.md` with a recommendation per omic and the claim limit.
- [ ] The user's decision recorded.

## Regression Checks

- [ ] `git diff` over `R/` empty for this ticket; W1's files pass.

## Backward Compatibility

- [ ] No behaviour, default, stored state or output changes.

## Rollback Note

- [ ] Evidence and tests are additive.

## Constraints

R 4.6.1; the pinned libraries are not installed into. No real-data values. Nothing under `tools/profiling/` changes. No `roxygen2::document()`. Commit locally, never push, no attribution lines.

## Revision 1.3 (2026-11-16): decisions from the S1 and S3 pilots, before any confirming run

The pilots (20 runs per cell, runs 9001 to 9020; `kit-124/s1-pilot/`, `s3-pilot/`) proved completion and cost. Decisions taken from their instrumentation facts only, not from which design looked better:

- **S3, cell C6 (technical duplicates) is reported, not core.** The transcriptomics path fits library replicates of one biological unit as independent samples (no handling of the design's technical-replicate column), so even the true-factor arm has a null rate of 0.18. By the rule the cell is void. The cause is the package's, not RUV's: parked row 23.
- **`ruv::RUVinv`'s own stop at 3 per group** ("fourth-smallest eigenvalue ... not currently supported") is a method refusal, not a crash.
- **Refused runs** carry their data set's no-RUV result (as written); arm 4's refusals ("k is not below the controls observed in every sample") are the package's rule and count as designed refusals.
- **S3 has no 50%-missing cell;** S5 and S14 are low-count cells (library size divided by 10).
- **E on data with missing values** is computed on the features complete in every sample.
- **Bias** is sign-aligned (negative = shrunk towards zero).
- **Cost:** S1 about 75 to 149 core-hours, S3 about 70 to 135, by machine load. Confirming runs start when the measured suite and journeys are off the machine.
- A void cell (the true-factor arm failing rules 1 or 2: S12 in S1 by the pilot) is reported and judges no design.

### Revision 1.3, additions from the S4 pilot (2026-11-16, before any confirming run)

Pilot: 20 runs per cell, runs 9001 to 9020, 22 cells, no run died (`kit-124/s4-pilot/`). Decisions from instrumentation facts only:

- **S4, cell C6 (duplicate injections) is reported, not core**, as in S3: the lipid and metab test has no replicate block, so the true-factor arm's null rate is 0.137. The cause is the package's (the technical-replicate column is read and not used); it joins parked row 23 as an omics-parity defect.
- **S4's missing share is the detection limit's target, not the share among tested features.** The package's sweep and test refuse a feature with no value in some group, so the simulator keeps features observed at least once per group. The realised share (C0 16%, S6 28%) is recorded per run in `s4-timing.tsv` and is the figure reported; the 25% and 50% labels name the detection limit. The limit is not moved to hit the label (infeasible in the cells without unwanted variation).
- **`ruv::RUVinv` with 5 internal standards** is a method refusal (needs as many controls as residual dimensions); its `invvar` eigenvalue errors at 3 per group are the same refusal, not a crash of a package design.
- **Drift is balanced against group by construction;** chance confounding of drift with group is outside the claim.
- **The silent empty result of `runLipidsDA` on a kernel error** is recorded as a defect (parked row 24), not fixed here.
- **Cost:** about 105 core-hours, under the ceiling.

### Revision 1.3, additions from the S2 pilot (2026-11-16, before any confirming run)

Pilot: 20 runs per core cell, 5 per stress cell, runs 9001 on (`kit-124/s2-pilot/`). Decisions from instrumentation facts only:

- **Engine: fastlimpadpc**, which `resolveLimpaDpcEngine` selects when it is installed (a default run on this machine). The two engines agreed to 0.0013 log2 on a probe. The claim is stated for that engine; limpa itself is not rerun.
- **S5 (no missing values): `limpa::dpcDE` fails in some runs** (constant DPC standard errors, NA variance-trend weights, "No residual degrees of freedom"). It is the route's behaviour on complete data, counted as crashes of the route in a stress cell, and parked as row 25. S5 is not a core cell.
- **The correlation of the factor with group is set in expectation in S2** (realised value per run in `s2-data.tsv`); S1, S3 and S4 state their own construction. The difference is reported with the results.
- **Arm 16 repeats arm 7 in C6** on this route (declared technical replicates are always the replicate structure); reported once.
- **Arm 6's refusals** ("RUV-III needs more than k control proteins") are designed refusals and carry the no-RUV result.
- **Cost:** 196 core-hours projected under load, about half unloaded (C0: 22 s per run). The ceiling of 150 core-hours stands. The confirming run is made on an otherwise idle machine, core cells first; if the ceiling is reached, the stress cells not finished are reported as not run and nothing is retried.

### Confirming studies: declaration before the run (2026-11-16)

- **Purpose:** the rev 1.2 rules on every arm, per simulator. **Tree:** pinned worktree `kit/wt-k124-sim` at `0b6fcb6d` (KIT-123 R1 and R2 are in it, so "today" is the package after the two bug repairs). **Libraries:** `kit/scratch-k124-lib` then `browser-smoke.y5ordS/library`; one BLAS/OpenMP thread.
- **Runs:** from run 1; 1,000 per core cell, 500 per stress cell; seeds by each script's rule. Output `kit-124/s{1,2,3,4}-confirm/`.
- **Order and cores:** 6 cores per simulator, two simulators at a time (S1 then S2; S3 then S4), core cells first, then all cells (`kit/scripts-k124/confirm.sh`).
- **Ceiling:** 150 core-hours per simulator = 25 h of wall at 6 cores over both phases, by `timeout`. Cells not finished at the ceiling are reported as not run.
- **Attempts:** one. A pause for a measured journey (stop, then the same command; finished cells are kept by each script's restart rule, proved by a two-cell lifecycle test, `kit/scripts-k124/lifecycle/`) is not a new attempt: seeds and runs are the same.
- **Decision:** the rule table per arm and cell, then `kit-124/outcome.md` and the user's gate.

## Revision 1.4 (2026-11-17, plan 4.7 audit M-04, M-11): how the gate is presented. No rule and no run changes.

- **Rule 4 in C4** (aggregation note A12) is put to the user with both readings side by side, as the first line of the table; neither is applied silently.
- **From verdict to option:** the rule proposed in KIT-129 revision 1.1 is shown with the table for the user's confirmation.
- **DIA_LIMPA (M-11):** S2 ran on the `fastlimpadpc` engine; the default engine is `limpa`. The outcome cites the evidence that the two give the same results (KIT-053's exactness record) or states the limit of the row.

## Gate: the user's decision (2026-11-18)

- Option set: tiered, every design offered in every omic with its measured figure and a note beside it ("tiered option above is perfect"; "didnt i say we offer all as options? with tooltips?"): Recommended / Valid with conditions (technical replicates; supplied controls; first-pass null controls) / Legacy, warned (today's design) / None. The choice is written into the report's methods text.
- Pre-selected default in all four omics: pseudo-replicates (groups), all features as controls, factors in the model (arm 7), k estimated or user-chosen, not fixed at 1 ("okay." after the per-omic table and the literature comparison; recorded as agreement).
- The paper states where each design follows its method paper and where it departs (factors in the model against the RUV-III papers' corrected-then-tested usage), with the simulation table.
- Open: k estimation (`sva::num.sv` is not a dependency); conditions of use from the stress cells still running (classic proteomics, DIA_LIMPA).
Evidence: `kit-124/outcome.md`, `kit-124/design-review-20261118.md`; dashboard https://claude.ai/artifact/WRjGxjQEczpETtMdVM3Ycc.

**User, 2026-11-18 (platform and budget).** "run all simulators on github, free up this laptop": every simulator moves to GitHub-hosted runners (branch `sim-k124`); result tables may be uploaded as artifacts for all four simulators. The local run is stopped when S1's lane ends; S2's extra-hard cells are not run on the laptop. Asked whether the unfinished extra-hard cells run to completion on GitHub, one attempt each, with the 6-hour job limit as the only ceiling: "yes." This replaces the per-simulator wall-time ceiling (a laptop protection) for the cells not finished locally; runner results sit beside laptop results only after the equivalence proof under the fixed CI tolerance rule.
