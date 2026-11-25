# KIT-255 outcome (gate G0 only)

Ticket `tickets/kit/KIT-255.md` **r1.2** (the plan audit's Re-check 1, RC1-01 to RC1-06, applied; E-4 ruled by the
coordinator). Branch `kit-255-g0`: scripts only, from **v0.9-staging 7327607f**, never merged; nothing under `R/`.
**OD-1, OD-2 and OD-3 are the owner's and open.** Nothing here assumes an answer; R1 and every code change wait for
OD-1.

Labels: OBSERVED = a runner's output (run id given). READ = by reading the code at the commit named. INFERRED =
follows from what was read or observed, not executed itself.

## Step 0 (READ, before G0)

**Branch point.** `git rev-parse origin/v0.9-staging` = 7327607f. `git diff --stat 29065062 7327607f` and
`git diff --stat e0747fc1 7327607f`: `docs/methods-decision-ledger.md` only. Every R and test line r1.1 cites (read at
29065062) therefore reads as r1.1 states. Spot checks at 7327607f: `omicDaPreflight` `R/mod_omic_da.R:398`; metab
`preflight_first = TRUE` `R/utils_omic_gui_config.R:2716` and its validator `:2723`; the alignment message
`R/func_omics_design_shared.R:142`; the per-assay column removal `R/func_lipid_s4_correlation_methods.R:393`, metab
`R/func_metab_s4_qc_methods.R:468`; `samples_to_remove_unique` lipid `:406`, metab `:481`.

**R2's scope, counted again.** The ticket's grep over `tests/testthat` without helpers: **26 files** (27 with
`helper-module-ci-metab-design.R`), as r1.1. The five files added by name (`test-gui-entry-lipidomics.R`,
`test-gui-entry-metabolomics.R`, `test-metab-03c-norm-module-seams.R`, `test-integration-input-lipid-metab.R`,
`test-kit207-summary-head-guard.R`) exist and none is in the grep's list: **31**, plus the new file = **32**, as r1.1.

**Fixtures that might remove a sample in one assay only.**
- Module-CI (`helper-module-ci-lipid-norm.R:308-335`; metab `helper-module-ci-metab-norm.R:387`): one pairs table
  given to every assay, so a sample fails in every assay or in none (READ). They cannot make the defect.
- `lipidDeclarationsObject` (two assays, seeds 31 and 32) and `metabDeclarationsObject` (LCMS_Pos and GCMS, seeds 31
  and 32): dense synthetic values; whether one of their 24 samples falls below 0.5 in one assay only is not decidable by
  reading (correlations are computed). Not observed in G0 (out of its declared items); recorded by R1 on base, as the
  ticket says.

## G0 (declared before the request)

- **Purpose:** observe on base, before any code, phase 0 (the managed harness's lifecycle check) and o1 to o5 for
  lipidomics and metabolomics, with K1, K2 and the control K0 (ticket section 10, r1.2). Observation only.
- **Binding:** v0.9-staging 7327607f plus `tools/kit255/g0.R` (this branch); R 4.6.1, reference BLAS, ubuntu-22.04.
- **Ceiling:** one job, 75 job-minutes (5.00 core-hours); script timeout 60 minutes. Expected about 35.
- **Predictions:** in the script (`k255Predictions`) and written to `predictions.tsv` before the package is loaded;
  o4's ("no metabolomics number changes") is **predicted by reading**.
- **What the result decides:** stop rule (0) (a prediction differs: report, the ticket is revised before R1); (i) (K2
  does not fail in exactly LCMS_Neg: instrument defect, the construction is fixed and the job resent, at most three);
  (ii) (o4 shows any metabolomics difference: the claim "no metabolomics number changes" is withdrawn and OD-1/OD-2 are
  put again); (iv) (a managed row ending with the head unchanged and the failure in the log is an instrument defect
  until phase 0 has passed). E7 is counted only if phase 0 passes.
- **Cited, not repeated:** KIT-242 G0 C4b (run 45594264982): the kernel's split on an emptied column (B_6 out of the
  design and A2, in A1), the lipid declaration-route refusal ("Assay 'A1' has sample column(s) missing from design:
  B_6") and the metabolomics declaration-route completion (tables for A1 and A2, 40 rows each).

## G0: result (OBSERVED)

Runs (`runs.tsv`): **45606617982** at 6fb5078d (staging 7327607f + script) and **45609166200** at fd9da180 (instrument
rerun 1 of 3, VAL-54). Tables: `g0-45606617982/`, `g0-45609166200/` (`predictions.tsv` written before the package was
loaded; `observations.tsv`; `conditions.tsv`; `provenance.txt`). R 4.6.1, reference BLAS, ubuntu-22.04.

**The instrument defects of the first run, and the fix.** Two recorders produced no result; nothing was judged on them:
- o3: the recorder searched the worker's raw handoff, which holds references only (`checkpoint`, `reference`, `native`,
  ...), and found nothing. Fix: the same worker run is accepted on a controller (`omicsWorkerAction`, as
  `omicsWorkerJourney` does) and the contrast's value and the manifest acceptance are read.
- o4, declaration arm: saving the transient object read from the store into a fresh owner stopped ("invalid hydration
  reference"). Fix: the fresh owner gets a plain object built from the same assay tables and design
  (`as.data.frame`), with `kit172Object`'s args (the metab DA reads no other arg here; `ruv_uncorrected_assays` is
  absent on both sides).
- No prediction, threshold, object or comparison changed (VAL-54 (i)). Every other observation of the rerun is
  identical to the first run's (`diff` of the two tables without `seconds`, o3 and o4).

**Observed against predicted** (both omics unless stated; values from 45609166200):

| Item | Predicted | Observed | Holds |
|---|---|---|---|
| Phase 0, managed run and Apply(0.5) on K2 | complete; head `<prefix>_correlation_filtered` | complete (7 actions delivered, no error in the log); head `<prefix>_correlation_filtered` | yes |
| Phase 0, managed design = in-process design | TRUE | TRUE; **PHASE 0 PASSED** (lipid, metab) | yes |
| Phase 0, B_6 after log2 and quantile | out of the design and LCMS_Neg, in LCMS_Pos, both routes | design 11, LCMS_Pos 12 with B_6, LCMS_Neg 11 without, in process and managed | yes |
| Stop rule (i), K2 | fails in exactly LCMS_Neg | B_6's best within-group r: design-state values LCMS_Pos 0.975839, LCMS_Neg -0.014881; after log2 and quantile 0.981240 and 0.012826 | yes (not triggered) |
| o1 (K2, declaration) | design 11 without B_6; LCMS_Pos 12 with; LCMS_Neg 11 without | as predicted | yes |
| K1, A2's between-sample step | warning "Assay 'A2': At least one sample column contains only NA values. Skipping normalization."; A2 equals its log2 input | warning present; A2 equal to its log2 input (non-IS rows, tolerance 0) TRUE; A1 FALSE (normalised) | yes |
| K1, correlation state | design 11 without B_6; A1 12 with; A2 11 without (as C4b) | as predicted | yes |
| o2 (K1) tab in process, controller | `status=error; stage=preflight; errorMessage=Assay 'A1' has sample column(s) missing from design: B_6` | exactly that, lipid and metab; one `[error]` notification with the message | yes |
| o2 (K2) | the same with `'LCMS_Pos'` | exactly that, lipid and metab; one `[error]` notification | yes |
| o3 (K2) lipid, worker | carries the message; no assay table | worker run complete; contrast value: no assay table, errors "Assay 'LCMS_Pos' has sample column(s) missing from design: B_6" (listed twice in the value); the manifest acceptance stops with that message | yes |
| o3 (K2) metab, worker | complete; tables for both assays; no error | complete; LCMS_Pos 520 rows, LCMS_Neg 520 rows; no error; manifest accepted | yes |
| o4 (K1, K2) metab, `runMetabolitesDA` (predicted by reading) | identical per assay and contrast | `identical` TRUE for A1, A2 (K1) and LCMS_Pos, LCMS_Neg (K2) | yes |
| o4 (K1, K2) metab, declaration route (predicted by reading) | identical | `identical` TRUE for every assay, K1 and K2 | yes |
| o5 (K1, K2) | writes; A1 / LCMS_Pos feature table has a column B_6; samples 11; every matrix 11 | writes (6 files); feature table columns `feature_id, display_name, assay, is_internal_standard, <id>, <annotation>, B_6`; samples 11 (no B_6); every matrix 11 sample columns | yes |
| K0 (control) | nothing removed; both DAs complete with both assays | design 12, both assays 12; both DAs: LCMS_Pos and LCMS_Neg, 520 rows each, no error | yes |

**Verdict: G0 passed.** Phase 0 passed for both omics (E7 may be counted). Stop rules (0), (i), (ii) and (iv) did not
fire. o4 confirms the claim "no metabolomics number changes under (a)" on both metabolomics routes and both objects
(stop rule (ii) not triggered).

**Inferred (not executed):** o2 shows that today the metabolomics tab without a worker refuses exactly as lipidomics
does, while o3 shows that the metabolomics worker route completes on the same object. This is the route split of ticket
2.4, now observed on the tab's controller and on the worker. o5 shows that the stray column reaches the integration
input's feature table as an annotation column. It reaches neither the matrices nor the samples table.

## Not checked

- E1 to E10 and G1 to G6 (R1 and R2 wait for OD-1).
- The lipid and metab declarations fixtures (`lipidDeclarationsObject`, `metabDeclarationsObject`): whether they
  remove a sample in one assay only; R1 records this on base.
- The record's per-step `dimensions` and the report's figures (ticket 13).
- Replay of an older record (E8, OD-2).
- RUV-applied objects: G0 ran RUV skipped only. The `ruv_uncorrected_assays` path was read, not run.
- Windows, the Electron build and real data.

## Core-hours

G0: two jobs, 13.6 + 13.5 = 27.1 job-minutes = **1.81 core-hours**: 0.91 against G0's expected 2.33 (ceiling 5.00),
plus 0.90 from the instrument budget (ceiling 15.00, one of three reruns used).
