# Todo: Metabolomics & Lipidomics Refactor Fixes

## Phase 1: S4 Logic Cleanup
- [ ] Refactor `metaboliteMissingValueImputationLimpa` (Remove `<<-`, use `tryCatch` return).
- [ ] Refactor `metaboliteMissingValueImputationMissForest` (Safe variable assignment).
- [ ] Refactor `lipidMissingValueImputationLimpa` (Remove `<<-`).
- [ ] Refactor `lipidMissingValueImputationMissForest` (Safe variable assignment).
- [ ] Implement group-size validation in `metaboliteIntensityFiltering` and `lipidIntensityFiltering`.

## Phase 2: Imputation GUI Enhancements
- [ ] Update `mod_metab_norm_impute.R` UI with Diagnostics sub-tab.
- [ ] Update `mod_metab_norm_impute.R` Server (DPC/OOB diagnostics + saving).
- [ ] Update `mod_lipid_norm_impute.R` UI with Diagnostics sub-tab.
- [ ] Update `mod_lipid_norm_impute.R` Server (DPC/OOB diagnostics + saving).
- [ ] Fix `half_min` warning logic in both servers.
- [ ] Ensure `cp02_imputed` capture in both servers.

## Phase 3: Validation
- [ ] Verify unit tests.
- [ ] Manual check of checkpoint generation.

<!-- APAF Bioinformatics | todo.md | Approved | 2026-03-18 -->
