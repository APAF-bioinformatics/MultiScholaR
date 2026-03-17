# Specification Report: Metabolomics & Lipidomics Refactor Fixes
# Trace ID: tr_20260318_063210_

## 1. Goal
Address reviewer feedback by enhancing the Imputation UI with diagnostics, fixing S4 implementation anti-patterns, and robustifying filtering logic.

## 2. Architecture Changes
### 2.1. GUI Enhancements (`mod_*_norm_impute.R`)
- **Sub-tab Architecture**: Add `tabsetPanel` inside the Imputation tab with "Controls" and "Diagnostics" sub-tabs.
- **Limpa Diagnostics**:
    - Display DPC plot (ggplot2) in the Diagnostics sub-tab.
    - Save DPC plot to `[analysis_dir]/QC/Imputation/`.
- **MissForest Diagnostics**:
    - Display OOB error summary table.
    - Save summary to `[analysis_dir]/QC/Imputation/`.

### 2.2. S4 Logic Refinement (`func_metab_s4_objects.R`, `func_lipid_s4_objects.R`)
- **TryCatch Refactor**: Replace `<<-` super-assignment with direct `tryCatch` return assignment.
- **MissForest I/O**: Refactor `missForest` calls to handle `imputed_matrix_t` assignments safely within the return value of `tryCatch`.
- **Filtering Logic**: 
    - Inject a check to compare `min_samples_per_group` against `design_matrix` group sizes.
    - Emit a `warning()` if the threshold is logically impossible for any group.

### 2.3. Imputation Improvements
- **Half-Min Fix**: Add `length(pos_vals) == 0` check to prevent `min(numeric(0))` warnings.
- **Checkpoint Alignment**: Ensure `mod_*_norm_impute.R` calls `.capture_checkpoint` with ID `cp02_imputed`.

## 3. Implementation Plan
- [ ] Refactor S4 methods in `R/func_metab_s4_objects.R` and `R/func_lipid_s4_objects.R` (logic + anti-patterns).
- [ ] Update `mod_metab_norm_impute.R` and `mod_lipid_norm_impute.R` UI (sub-tabs).
- [ ] Update `mod_metab_norm_impute.R` and `mod_lipid_norm_impute.R` Server (plotting, saving, diagnostics).
- [ ] Add strict group size validation to filtering S4 methods.
- [ ] Verify fixes with existing and new test cases.

<!-- APAF Bioinformatics | specification_report_refactor_fixes.md | Approved | 2026-03-18 -->
