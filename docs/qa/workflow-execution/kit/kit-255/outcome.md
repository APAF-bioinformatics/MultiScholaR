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
