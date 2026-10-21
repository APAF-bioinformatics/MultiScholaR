# KIT-053: DIA_LIMPA differential abundance with dpcDE and W covariates, enrichment, report and record; the journey measured

Status 2026-10-21: **STATUS_PLACEHOLDER** (ticket revision 1.1). Implemented at `CODE_SHA`; install `kit_053` = `core/browser-smoke.INSTALL_ID`; counted: one single journey member on the 300-sample DIA-NN corpus (authorised in advance by the user on 2026-10-21 with the afk authorisation). The `da_limpa` step-map unit is recorded in the plan's decision record for the user's review.

## One mechanism: the DA kernel's limpa branch, and the route's closure

- **DA branch** (`differentialAbundanceAnalysisHelper`, `R/func_prot_da_model_methods.R`): an object carrying the quantified EList (the DPC quantification's results, reattached from its owned value) is fitted with `limpa::dpcDE` on the EList subset to the design's samples, the design built from the user's formula plus the `W1..Wk` columns present in the design (the recorded `formula_used`), `plot = FALSE`, and KIT-035's resolved technical-replicate block passed through (`dpcDE` forwards it to `voomaLmFitWithImputation`); `convertDpcDAToStandardFormat` converts the fit (contrasts.fit, eBayes with the trend/robust flags, topTable), with the q-values now through the limma path's own `adjustProtDaResultTable` (qvalue, BH fallback, the fallback recorded as a q-value warning as on the limma path); the converter's tables are normalised to the limma tables' shape (ids as row names; the shared code adds the id and comparison columns and derives the log intensity). TREAT does not apply on this branch and is recorded (`treat_applied = FALSE`, the requested cutoff kept). Every other object takes the limma path unchanged. The helper returns `da_model` (`model`, `formula_used`, `w_terms`, `k`, `block_reason`, `uses_duplicate_correlation`, `treat_applied`); the workbook wrapper `daAnalysisWrapperFunction` is left as it was (no caller in `R/`; the deletion list of KIT-048).
- **Declaration, identities and manifest** (`R/utils_workflow_dia_da.R`): `diaContrastOperation(route)` carries the route factor; `diaActiveRoute()` gives `limpa` on the `DIA_LIMPA` type and `iq` otherwise (every DIA, LFQ and TMT identity unchanged); the base-fit and full-result identities read the same; the contrast product carries `model`; the manifest's `model` is the kernel's account (`dpcDE` or `limma`, the formula used, the W terms, k, the block) with the limma base fit's factors beside it and the route.
- **Record and replay** (`R/utils_workflow_analysis_record.R`): `analyses$differential_abundance$model` names the model; the report overlay gains the route's values (the DPC engine and version, the fit and band, the normalisation method, the peptide-count thresholds, the covariates' k, percentage, controls and replicate structure, the DA model, the formula used, k, the block reason); the replay drives the route's operations (the worker handlers for the peptide and protein stages, the correlation from the head state when no normalisation checkpoint exists, `protein_ruv_covariates` committed as the runner commits it, then the DA on the record's parameters).
- **Report** (`inst/reports/proteomics/DIANN_limpa_report.rmd`): on a record naming the DPC quantification, the method text is the route's (no intensity, sample or replicate filters; the log2 matrix and the peptide-level accession cleanup; DPC-Quant on the named engine with its fit and band, Li et al. 2025; the peptide-count filter; no protein missing-value filter; the normalisation as run, none by default; RUV-III as covariates of the model, never a correction) and the DA text names dpcDE with the W terms and no TREAT; the five passages the ticket names are no longer printed on the route; the DIA-NN template and the limpa template's non-route text are untouched except the DPC citation.
- **Normalisation tab on the route**: the RUV covariates are the tab's own action after the correlation sub-step (`apply_ruv_covariates`, rendered on `DIA_LIMPA` only, with its summary), the limpa script's order (correlation, then RUV); the KIT-041 hook on the correlation click is removed. The step list has `correlation` then `ruv_covariates`.
- **Harness**: `kitJourneySteps("DIA_LIMPA")` = import → design → qvalue → rollup → peptide_log2 → peptide_accession → dpc_quant → protein_normalisation_limpa → peptide_count → correlation → ruv_covariates → session_export → da; the driver's route import sets the entry point's radio before the pinned import actions (the pinned helper untouched), the route's tabs and actions, the normalisation pass-through before the correlation, the covariates click; a ticket may bind the tab's normalisation settings (`norm_method = none`, `ruv_mode = covariates`), part of the case binding; `kit-001/step-map-limpa.json` (the route's rows and the shared steps) hashed as `step_map_limpa_sha256`; ticket `053` bound as a whole-journey single member on the route with oracles `limpa_<corpus>` (`KIT_053_CORPUS=screen` for the lifecycle run); the route's oracles recorded by `kit-053/oracle-record.R` and bound as `limpa_<corpus>_<step>`.
- **Reproduced defect repaired**: `filterMinNumPeptidesPerProtein` subset the quantified EList by `genes$protein.id`, absent on limpa 1.4.x (KIT-047 found the same column missing for the quantification), so the EList was emptied and the DA had nothing to fit; the ids are the row names of `E` when the column is absent, and the standard errors and observation counts are subset with it. KIT-041's screen check compared the emptied EList on both paths (equal, both empty); its matrix comparisons stand.

## Evidence

SCREEN_EVIDENCE

## Counted member

COUNTED_EVIDENCE

## Deviations and claim limits

- The `da_limpa` unit was recorded without the user's confirmation row (afk; decision record).
- The route's RUV covariates are a separate tab action after the correlation filter (the ticket's "through the managed step" reads as the in-process whole-matrix declaration, the `normalization-ruv` pattern; KIT-024 refuses worker dispatch of a whole-matrix in-process declaration).
- `daAnalysisWrapperFunction` is neither reduced nor removed (no caller in `R/`; `test-prot-07d` keeps its test); its own dpcDE call (no block, the old q-value call) is the workbook's and is on KIT-048's deletion list.
- The script-equivalence run uses the route's replicate rule (KIT-035); on the screen design (no technical replicates) it coincides with the script's `group`.
- The q-value fallback on the dpcDE branch: on small protein sets `qvalue()` cannot estimate pi0 and BH is used, recorded as a q-value warning (the limma path's rule); on the 2,500-protein fixture qvalue() succeeds.

## Rollback

`git revert CODE_SHA` (and the closure commits): the `DIA_LIMPA` route then stops at DA on the previous build (its DA kernel has no limpa branch and the route factor is the constant); DIA route stores and identities are unaffected.

## Hand-off

The `DIA_LIMPA` route is complete on the owner. KIT-044 takes both DIA-NN routes for its capability fallbacks; KIT-048 the deletion list (the workbook DA wrapper's dpcDE branch, the `use_limpa` flags); KIT-050 the Limpa populations and oracles for CASE-PERF-JOURNEY-002 on a non-synthetic corpus and the `fast_limpa_dpc` publication requirement.
