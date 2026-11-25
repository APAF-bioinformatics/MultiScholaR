# KIT-262, step 0: base, re-read lines, and the owner's answers

Ticket: KIT-262 r1.1 (`tickets/kit/KIT-262.md`, git-ignored; SHA-256 of the file as read: 873d551847b11333...). Plan register: AUDIT-KIT-262-plan-20261126. Scope of this implementation pass: step 0, M1 and M2 and the skip rule's verdict. M3 is not started (it needs KIT-242 Part B on staging and the coordinator's go). No package code, no NEWS, no version change.

## Base

Branch `kit-262`, from `origin/v0.9-staging` at **15d86bd4e3bb2438d7f198b18742439176b49079** (merge train 18). KIT-190 (train 15), KIT-217 (train 15) and KIT-204 are on it. At the time of branching `origin/v0.9-staging` had moved to 176f9d3f (a ledger row only; no file under `R/` changed).

## The owner's answers (VAL-302, 2026-11-26, "all as recommended."), as taken here

| Item | Choice taken | Effect on this pass |
|---|---|---|
| O-1 | (ii): C2 (proteomics' rule) is reported for coverage and k in M2 only, not tested in M3 and not recommended; C1 is tested at one-sided 0.05 | `K262_O1 = ii`: the curve hook is used (C2 needs each percentage's curve); C2's usable set, winner and k are written beside C1's in M1 and M2 |
| O-2 | (a): "no worse" on cell averages with VAL-185's margins; changed-data-set differences and the present rule's baselines beside the verdict | no effect before M3 |
| O-3 | yes to both: (1) the two 24-sample stress cells, (2) the M3 skip rule | the stress-cell choice (largest C1-not-P share among the six declared 24-sample cells) and the skip rule are computed from M2 |
| O-4 | (a): 18 cells of 4 samples (2 groups of 2), reported only: k and residual degrees of freedom; no error rates at 2 x 2 | M2 has 198 cells; cells 181 to 198 are excluded from the gain rule and from the stress-cell selection |

The ticket text is not edited.

## Cited lines re-read on this base

The sweep and its twin are byte-identical to `origin/kit-190-ruv-control-selection` (git diff of the three files: empty); the cancor-once test file is also identical (it carries KIT-190's trace print).

| What | Ticket's line (kit-190 6b465f1e) | On this base |
|---|---|---|
| Lipid sweep: percentages `seq(percentage_min, percentage_max, by = 1)` | :82 | `R/func_lipid_norm_ruv_optimization_helpers.R:82` |
| R1 / R2 (error statuses) | :93-128, :168-188 | :93-128, :169-172 |
| R3 (fewer than 5 controls) | :203, :210-225 | :203-220 |
| R4 `invalid_cancor_plot` | :232-252 | :247 |
| R5 `isTRUE(best_k > complete_controls)` | :257 | :257 (metabolomics twin :582) |
| R6 `error_scoring` | :273-293 | :288 (twin :613) |
| winner: composite, then lowest k, highest separation, lowest percentage | :330-337 | :330 onwards |
| `omicsRuvCompleteRows`, `omicsRuvCompleteControlCount` | `R/func_metab_norm.R:264-282` | :264, :274 |
| `findBestKElbow`, `calculateSeparationScore`, `calculateCompositeScore`, `calculateAdaptiveMaxK` | `R/func_prot_norm_optimization_helpers.R:64-91`, :219-244, :270-288 | :64, :172, :219, :270 |
| proteomics' truncated curve before the elbow | :901-928 (kit-190) | `R/func_prot_norm_optimization_helpers.R:905` (bounded curve), :928 (R5) |
| the cap after the sweep | `R/utils_workflow_omics_declarations.R:541-570` | `omicsRuvKBound` at :736; `boundOmicsRuvSweepK` at `R/func_omics_ruv_factors.R:109` |
| the filter the reference step calls | (named in the declaration) | `lipidomicsDeclarationBinding$filter_intensity` (`R/utils_workflow_lipid_declarations.R`) = `lipidIntensityFiltering(..., lipids_min_values_per_group, lipids_min_groups)`; the log2 step = `logTransformAssays(theObject, offset = 0)`; the sweep = `binding$optimise_ruv` = `runLipidPerAssayRuvOptimization` |

Contradictions with the ticket's text found: none, apart from line numbers that moved in `omicsRuvKBound`'s file (above).

## What the sweep's trace does not carry (resolved before any run)

`optimization_results` holds status, k, separation, composite and the error reason; the number of complete controls is only in the reason of R5 rows. M1 recomputes it for every row with the package's own helpers (control selection, `omicsRuvCompleteRows`, `omicsRuvCompleteControlCount`). In M2 it is read from R5 reasons only (the ANOVA per percentage would double the cost); for the other rows it is not recorded, and nothing needs it: a row with status `ok` has k <= c, so it is admissible under every rule.
