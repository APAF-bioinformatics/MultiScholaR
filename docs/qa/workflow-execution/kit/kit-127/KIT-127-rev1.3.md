---
tags: ["ticket", "kit", "ruv", "lipidomics", "metabolomics"]
id: "KIT-127"
title: "Lipidomics and metabolomics RUV: levels kept, and the factors in the model"
description: "The two repairs proteomics has: each feature keeps its mean after RUV-III (KIT-123 R2), and the DA reads uncorrected values with the factors as covariates."
status: in_progress
priority: "P1"
dependencies: ["KIT-124"]
time_estimate: "12h"
date-created: "2026-11-17"
date-updated: "2026-11-18"
needs_planning: false
ticket_revision: "1.3"
plan: "00-KIT-PLAN.md"
files_to_modify: ["R/func_omics_ruv_factors.R", "R/func_lipid_da_model.R", "R/func_metab_da_orchestration.R", "DESCRIPTION", "man/runLipidsDA.Rd", "man/runMetabolitesDA.Rd", "tests/testthat/test-kit127-omics-ruv-factors.R", "NEWS.md"]
---

# KIT-127: Lipidomics and metabolomics RUV: levels kept, and the factors in the model

**Authority.** User, 2026-11-17: the road to the submission (plan 4.7): "ticket series looks great. LETS GO."; RUV: "calibrated default, other valid options selectable, subtract-then-test as a warned legacy switch": "yes"; the distributable is Electron around a bundled R, library and Python, Linux and Windows; 31 December 2026 is hard and the paper claims integration and Windows; no comparison with other tools before submission; GitHub Actions on the public repository (push of `performance-v2` allowed, "yes."). User, 2026-11-18 (KIT-124's gate): arm 7 the pre-selected default in all four omics, k estimated or user-chosen and not fixed at 1; every design an option; today's design the warned legacy option.

## Bounded outcome

Parity with proteomics: (1) the level shift after the correction is removed (the shared kernel fix of KIT-123 applied to the lipid and metab methods); (2) the DA uses the factors in the model; the corrected matrix stays for figures and export. Internal standards as a user-supplied control list where the gate chooses it.

## Exit

Failing-first tests; KIT-124's S4 arm reproduced through the package's entry point; no re-record here (KIT-129).

## Revision 1.1 (2026-11-17, plan 4.7 audit M-01)

Depends on KIT-124 only (KIT-125 runs beside it). The refusal rule of KIT-126 revision 1.1 applies here.

## Revision 1.2 (2026-11-17, KIT-060 D15)

The corrected matrix stays a product of the RUV step (for integration and figures) after the DA reads the factors in the model.

## Revision 1.3 (2026-11-18): authored in full (W0), after the gate

This text is the ticket. The ticket's own file (`tickets/kit/KIT-127.md`, not under version control) lives in the main checkout, which the authoring session could not write to; this copy is to be put there.

Read for this: AGENTS.md and the four convention files; KIT-124 (the gate of 2026-11-18), `kit-124/outcome.md`, `kit-124/design-review-20261118.md`, `kit-123/outcome.md`, KIT-128 (its branch `kit-128`, not merged), KIT-126 and KIT-129 (stubs), and the simulator `kit-124/sim/s4_lipid_metab.R` at `0b6fcb6d` (arms 3, 4, 7, 8, Qs).

### The work item

- **Outcome:** a kernel, the same for both omics, that (a) computes the RUV-III factors W per assay for the two gate designs and lets the differential test read the uncorrected matrix with W in the limma model, and (b) returns the corrected matrix with each feature's level kept. Today's path stays as it is and is what runs by default.
- **Not here:** GUI wiring, the workflow declarations and their step maps, the default, any re-record (all KIT-129); technical replicates as the RUV replicate structure (below, D3).
- **Exit:** the tests below pass; the legacy files are unchanged; the affected test files pass one at a time.

### What the simulator measured (the reference)

Per assay, on the log2 matrix of the study samples (`s4_lipid_metab.R`: `modelArm`, `ruvW`, `testWithW`, `daKernel`):

1. Y = the features observed in every sample (injections by rows); M = the replicate matrix of the groups.
2. Controls: every such feature (arm 7), or the internal standards, which are rows beside the features and are not tested (arm 8).
3. W = `Y[, ctl]` centred on its column means, times `ac' (ac ac')^-1`, `ac` the first k rows of `ruv::RUVIII(Y, M, ctl, k, return.info = TRUE)$fullalpha` on the controls (the formula of `computeRuvCovariates`, KIT-123 R2); centred; named W1..Wk.
4. The test: `runTestsContrastsLipidDA` on the uncorrected matrix, missing values kept, with the design plus the W columns and the formula `~ 0 + group + W1 + ... + Wk`, trend and robust, no TREAT. A feature whose observed samples cannot estimate every coefficient of that model is left out (it would stop the fit for the whole assay).
5. Refusals, scored as no RUV for that assay: k not estimable; k above the controls observed in every sample; no residual degree of freedom left (k above samples minus replicate sets minus 1). k = 0 is no RUV.
6. k was E, `sva::num.sv(method = "be")`, per assay. Arm 11 (k = 1) lost its power in lipid and metab (0.07 against 0.43) and with three factors.

### The kernel (new file `R/func_omics_ruv_factors.R`, internal functions)

| Function | Does |
|---|---|
| `runOmicsRuvFactors(theObject, ruv_design, ruv_number_k, ruv_k_source, ruv_control_ids, ruv_grouping_variable)` | Steps 1 to 3 and 5 per assay; then the corrected object. Returns `list(corrected, factors, parameters)`. |
| `boundOmicsRuvSweepK(theObject, ruv_results)` | Today's estimate of k per assay (the automatic sweep's `best_k`), limited to `calculateAdaptiveMaxK(samples)`; returns k with the sweep's value and the limit. |
| `applyOmicsRuvFactorsToDaInput(expr_matrix, design_matrix, formula_string, sample_id_col, assay_factors)` | Step 4's input: the supplied controls removed, the W columns on the design, the formula extended, the features that cannot be estimated removed and named. |

The DA entry points `runLipidsDA` and `runMetabolitesDA` gain `ruv_factors = NULL`. `NULL` is today's behaviour, unchanged. With the kernel's `factors` and the **uncorrected** object they run step 4 per assay through the unchanged fit and contrast halves.

**Declared parameters** (the kernel's arguments; KIT-129 lifts them into the `<prefix>_ruv` and `<prefix>_da_contrast` declarations with its harness revision):

| Parameter | Values | Note |
|---|---|---|
| `ruv_design` | `factors_all_features` (arm 7), `factors_supplied_controls` (arm 8) | No default. Today's design is not a value: it is today's functions, untouched. |
| `ruv_number_k` | a whole number of 0 or more, or one per assay by name | No default; never fixed at 1. 0 is no RUV. |
| `ruv_k_source` | `user`, `sweep_bounded` | No default. Recorded with k per assay; a label, not a switch. |
| `ruv_control_ids` | feature ids per assay | Required for `factors_supplied_controls`, refused otherwise. |
| `ruv_grouping_variable` | a design column | The replicate sets (groups: pseudo-replicates, as today). Default: the object's group column. |

Recorded per assay in `factors`: status (`ok`, `none`, `refused`) and reason, W, k, the k asked for and its source, the control ids and their count, the number of features observed in every sample, the ids not tested, the design and the grouping variable.

**The corrected matrix (b).** `ruvIII_C_Varying` as it is (`RUVIIIC::RUVIII_C_Varying`, missing values kept) with the same replicate sets, k and controls (all features, or the supplied ids), then `restoreRuvFeatureLevels` per assay against the uncorrected values. An assay with no RUV (k = 0 or refused) is returned as it came. An assay the correction drops or leaves unchanged stops the run (class `omics_ruv_correction_failed`).

**Refusal (KIT-126 revision 1.1).** An assay whose factors cannot be computed gets no RUV, with a warning of class `omics_ruv_factors_refused` that names the assay and the reason; never a silent fallback. A wrong argument stops (class `omics_ruv_factors_invalid`).

### Decisions taken at authoring (for the user to overrule)

- **D1. The workflow declarations are not touched.** Adding a parameter to `<prefix>_ruv` or `<prefix>_da_contrast` changes what every run records and the harness step maps under `tools/profiling/`, which this ticket may not edit and KIT-129 revises once. "Declared parameters" are therefore the kernel's arguments above.
- **D2. Supplied controls are not tested.** Arm 8 measured standards beside the features. In the package they are rows of the assay, so the DA leaves them out, as the ITSD normalisation removes its standards. They stay in the corrected matrix.
- **D3. Technical replicates as the RUV replicate structure (arms 15 and 16; KIT-128 D3) are not in this cut.** Their only S4 evidence is cell C6, "reported, not core", measured with the unblocked test and to be re-measured after KIT-128. The replicate sets here are the groups.
- **D4. k.** `sva` is not added. Today's estimate is kept as the automatic sweep's k, bounded by the sample-size limit. Neither a user's k nor the bounded sweep's k was an arm of S4 with all features as controls: arm 7 used E. This is the open point of the gate.
- **D5. The corrected matrix is RUV-III-C's own correction, not Y minus W alpha.** W exists only for the features observed in every sample; RUV-III-C corrects features with missing values too, and it is what arms 10, 15 and Qs called. The corrected matrix and W share controls, k and replicate sets, not arithmetic.

### (b) and the stored references

(b) is implemented behind the new path only: `R/func_lipid_s4_ruv_methods.R` and `R/func_metab_s4_ruv_methods.R` are not edited, so today's correction is bit-identical. Which stored references (b) would move if it were applied to today's path is measured by a probe (a scratch copy with `restoreRuvFeatureLevels` in both methods, the candidate test files one at a time) and listed under "Status".

### Where this meets KIT-128 (branch `kit-128`, not merged)

- `runLipidsDA` and `runMetabolitesDA`: KIT-128 passes `group_id` and `technical_replicate_id` to the fit half and reports failed assays; this ticket changes the `data`, `design_matrix` and `formula_string` lines of the same call and adds one block before it. Different lines; both man pages are edited by both.
- The fit halves (`fitLipidDaAssayModel`, `fitMetabDaAssayModel`) are not edited here. After the merge a design with declared technical replicates and `ruv_factors` fits the block on the design with the W columns (`duplicateCorrelation` with W in the model). That combination is not an arm of KIT-124.
- `NEWS.md`: both add under "Scientific changes".

### Tests (failing first), `tests/testthat/test-kit127-omics-ruv-factors.R`

On a seeded fixture (two assays, two groups of six, a sample-wise factor, missing values under a detection limit, five standards per assay):

1. Arm 7 through `runOmicsRuvFactors` and `runLipidsDA(ruv_factors =)` equals the simulator's computation (its `ruvW`, `testWithW` and `daKernel`, copied into the test): W, the tested features, logFC, p and q.
2. Arm 8 likewise, with the standards as rows of the assays.
3. Metabolomics gives the same W and tables as lipidomics (`runMetabolitesDA`).
4. The corrected matrix keeps each feature's mean and equals today's correction up to a constant per feature.
5. Today's path: `ruvIII_C_Varying` still shifts levels; `runLipidsDA` and `runMetabolitesDA` without `ruv_factors` equal the unchanged kernel on the same matrix.
6. Refusals, k = 0, the k record, the bounded sweep's k, wrong arguments.

### Verification budget

Targeted file first, then the lipid and metab DA and RUV files that call the edited functions, one R process at a time under `nice -n 19` (the machine is shared). One correction pass. No full suite, no journey, no re-record (KIT-129).

## Acceptance Criteria

- [ ] Arms 7 and 8 reproduced through the package's entry points (tests 1 to 3).
- [ ] The corrected matrix of the new path keeps each feature's level (test 4).
- [ ] Today's path bit-identical: the two legacy method files unchanged; test 5; the affected existing files pass.
- [ ] k has no default, is recorded with its source, and the bounded estimate is available (test 6).
- [ ] The stored references (b) would move, listed.
- [ ] NEWS under "Scientific changes".

## Regression Checks

- [ ] `git diff performance-v2 -- R/func_lipid_s4_ruv_methods.R R/func_metab_s4_ruv_methods.R R/utils_workflow_omics_declarations.R` empty.
- [ ] The lipid and metab DA and RUV test files pass one at a time.

## Backward Compatibility

- [ ] No default, stored state, oracle or journey changes. `runLipidsDA` and `runMetabolitesDA` gain a last argument that defaults to today's behaviour.

## Rollback Note

- [ ] Additive: one new file, one optional argument in two functions.

## Constraints

R 4.6.1; the pinned libraries are not installed into. No edits to `tests/testthat/helper-workflow-checkpoint-fixtures.R`, stored oracles or `tools/profiling/`. No `roxygen2::document()`. Branch `kit-127` from `performance-v2`; commit locally, never push, no attribution lines.

## Status

(filled at close)
