# Specification Report: Metabolomics and Lipidomics GUI Refactor

## 1. Goal
Refactor the monolithic `mod_metab_norm.R` and `mod_lipid_norm.R` files into a modular architecture matching Proteomics 2.0. Implement a `testthat` harness using checkpoint-based verification.

## 2. Architecture
### 2.1. Metabolomics Normalization (`mod_metab_norm.R`)
- **Orchestrator**: `mod_metab_norm.R` will become a lightweight coordinator.
- **Sub-modules**:
    - `mod_metab_norm_impute.R`: New module for missing value imputation (KNN, Min, Zero).
    - `mod_metab_norm_itsd.R`: Extracted logic for Internal Standard normalization.
    - `mod_metab_norm_scale.R`: Extracted logic for Log2 transform and between-sample normalization (Loess, Quantile, etc.).

### 2.2. Lipidomics Normalization (`mod_lipid_norm.R`)
- **Orchestrator**: `mod_lipid_norm.R`.
- **Sub-modules**:
    - `mod_lipid_norm_impute.R`
    - `mod_lipid_norm_itsd.R`
    - `mod_lipid_norm_scale.R`

### 2.3. Testing Harness
- Implement `test-metab-*.R` and `test-lipid-*.R` using `testthat`.
- Use `.capture_checkpoint()` to generate gold-standard RDS fixtures.
- Since interactive capture is not possible for the agent, a mock script `scripts/generate_norm_checkpoints.R` will be created to simulate the reactive environment and save the fixtures.

## 3. Implementation Plan
### Phase 1: Setup & Mocking
- [ ] Create `scripts/generate_norm_checkpoints.R` to generate baseline RDS fixtures for Metabolomics and Lipidomics using the `mock` dataset.
- [ ] Verify fixtures are generated in `tests/fixtures/`.

### Phase 2: Modularization (Metabolomics)
- [ ] Extract Imputation logic to `mod_metab_norm_impute.R`.
- [ ] Extract ITSD logic to `mod_metab_norm_itsd.R`.
- [ ] Extract Scaling/Normalization logic to `mod_metab_norm_scale.R`.
- [ ] Update `mod_metab_norm.R` to orchestrate these sub-modules.

### Phase 3: Modularization (Lipidomics)
- [ ] Extract Imputation logic to `mod_lipid_norm_impute.R`.
- [ ] Extract ITSD logic to `mod_lipid_norm_itsd.R`.
- [ ] Extract Scaling/Normalization logic to `mod_lipid_norm_scale.R`.
- [ ] Update `mod_lipid_norm.R` to orchestrate these sub-modules.

### Phase 4: Verification
- [ ] Run `devtools::test(filter="metab")`.
- [ ] Run `devtools::test(filter="lipid")`.
- [ ] Verify UI consistency in the Shiny application.

## 4. Constraints
- Zero tolerance for imperative loops (use `purrr`).
- Maintain APAF watermarks.
- Follow existing S4 object patterns.

<!-- APAF Bioinformatics | specification_report.md | Approved | 2026-03-17 -->
