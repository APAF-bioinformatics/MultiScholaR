# Demonolithing Parity Audit

Date: `2026-04-26`  
Branch: `janitor`  
Repo: `/home/doktersmol/Documents/MultiScholaR`

## Executive Summary

This closeout reached the repo-wide coverage target and preserved the previously cleared demonolithing parity bundle gates.

- Package coverage closeout: `63,105 / 70,076 = 90.0522%`
- Margin over `90.0%` gate: `36` covered lines
- Bundle parity status from the last completed FA-C4 refresh: `233 / 233` bundles passed
- Low-coverage bundles: `0`
- Comparison gaps: `0`
- Regressions vs baseline: `0`

The parity conclusion here is evidence-layered rather than a single fresh monolithic rerun. Bundle parity is carried from the completed run `coverage-evidence-20260423Tglobal-fa-c4-refresh-bp`; package coverage closeout is established by an exact line-union over the final base probe `coverage-target-probe-plus-20260425c` plus `38` focused follow-up artifacts.

## Evidence Base

Primary evidence used for this report:

- Bundle parity report: `.refactor-fidelity-audit/reports/latest-coverage-evidence.md`
- Bundle manifest: `.refactor-fidelity-audit/reports/latest-bundles.json`
- Dated summary artifact: `.refactor-fidelity-audit/reports/demonolithing-parity-audit-20260426.json`

Method notes:

- Bundle parity is not inferred from coverage alone.
- Package coverage is the exact union of covered executable lines across `39` captured artifacts.
- The final `90.0522%` figure is therefore a deterministic threshold closeout, even though it is not a single all-in-one rerun.
- No fresh integrated closeout was rerun after the last coverage wave; this report combines the last completed parity evidence run with the final package-coverage union.

## Headline Findings

| Layer | Result | Status |
| --- | ---: | --- |
| Bundle parity | `233 / 233` pass | Cleared |
| Package coverage | `63,105 / 70,076` | Cleared |
| Package coverage percent | `90.0522%` | Cleared |
| Margin over `90.0%` | `36` lines | Cleared |

## Per-Omic Breakdown

### Proteomics

Manifest bundle context:

- Bundle count: `92`
- Bundle status in completed FA-C4 evidence: no open parity failures surfaced in the proteomics slice

Coverage status:

- Covered: `23,866 / 26,234`
- Coverage: `90.9735%`

What materially moved:

- Proteomics was the largest absolute reservoir and remained the largest covered slice at closeout.
- High-leverage focused/shared tests were added or hardened around DA rendering, DA model stats, FASTA processing, QC coordination, peptide support, peptide replicate filters, S4 normalization methods, and peptide/protein workflow seams.
- The proteomics volcano helper branch issue was also normalized by qualifying `stringr::fixed()` inside `generateProtDAVolcanoStatic()`.

Representative coverage drivers:

- `tests/testthat/test-prot-da-render-direct-shared.R`
- `tests/testthat/test-prot-07d-da-model-stats-shared.R`
- `tests/testthat/test-prot-qc-coordinator-direct-shared.R`
- `tests/testthat/test-prot-qc-correlation-direct-shared.R`
- `tests/testthat/test-prot-s4-norm-methods-direct-shared.R`
- `tests/testthat/test-prot-fasta-processing-direct-shared.R`

Residual concentration after closeout:

- `R/func_prot_qc_peptide_methods.R` `26 / 196` covered, `170` uncovered
- `R/func_prot_da_model_methods.R` `308 / 418` covered, `110` uncovered
- `R/func_prot_da_heatmap.R` `98 / 185` covered, `87` uncovered
- `R/mod_prot_qc_protein_intensity.R` `218 / 301` covered, `83` uncovered
- `R/mod_prot_qc_peptide_intensity.R` `233 / 315` covered, `82` uncovered

Assessment:

- Proteomics is audit-acceptable for the current gate and is over `90%` within its own file slice.
- Further movement beyond the current repo target would still likely come substantially from proteomics because of its total surface area.

### Metabolomics

Manifest bundle context:

- Bundle count: `42`
- Bundle status in completed FA-C4 evidence: no open parity failures surfaced in the metabolomics slice

Coverage status:

- Covered: `14,071 / 15,120`
- Coverage: `93.0622%`

What materially moved:

- Metabolomics ended as the strongest omics-specific slice by package coverage.
- Focused coverage work landed around import metrics/organisms, DA core engine behavior, volcano/static plotting, S4 DA methods, S4 plotting, standalone normalization, QC ITSD behavior, and QC helper wrappers.

Representative coverage drivers:

- `tests/testthat/test-metab-da-core-engine-direct-shared.R`
- `tests/testthat/test-metab-da-volcano-static-direct-shared.R`
- `tests/testthat/test-metab-s4-da-methods-direct-shared.R`
- `tests/testthat/test-metab-s4-plotting-direct-shared.R`
- `tests/testthat/test-metab-import-metrics-organisms-direct-shared.R`
- `tests/testthat/test-metab-norm-standalone-shared.R`

Residual concentration after closeout:

- `R/mod_metab_norm_server_helpers.R` `1,799 / 1,977` covered, `178` uncovered
- `R/func_metab_s4_norm_methods.R` `1,539 / 1,644` covered, `105` uncovered
- `R/func_metab_s4_da_methods.R` `2 / 97` covered, `95` uncovered
- `R/func_metab_norm_support_helpers.R` `48 / 125` covered, `77` uncovered
- `R/mod_metab_qc_itsd.R` `219 / 287` covered, `68` uncovered

Assessment:

- Metabolomics is in the best shape of the three main omics families.
- Most remaining metabolomics work is deeper normalization/server-helper coverage rather than broad parity risk.

### Lipidomics

Manifest bundle context:

- Bundle count: `64`
- Bundle status in completed FA-C4 evidence: no open parity failures surfaced in the lipidomics slice

Coverage status:

- Covered: `13,430 / 15,142`
- Coverage: `88.6937%`

What materially moved:

- Lipidomics closed several previously cold reservoirs but remains under `90%` within its own slice even though the repo-wide target is cleared.
- The biggest gains came from DA export, S4 plotting methods, S4 DA plot methods, normalization/RUV helpers, and characterization coverage around import, QC, normalization, and design seams.

Representative coverage drivers:

- `tests/testthat/test-lipid-da-export-shared.R`
- `tests/testthat/test-lipid-12b-norm-ruv-helpers-shared.R`
- `tests/testthat/test-lipid-10-s4-plotting-methods.R`
- `tests/testthat/test-lipid-07-s4-da-plot-methods.R`
- `tests/testthat/test-lipid-12a-norm-module-characterization.R`

Residual concentration after closeout:

- `R/func_lipid_da_plotting.R` `151 / 432` covered, `281` uncovered
- `R/func_lipid_s4_plotting_methods.R` `327 / 504` covered, `177` uncovered
- `R/func_lipid_norm_support_helpers.R` `48 / 212` covered, `164` uncovered
- `R/mod_lipid_import_server_helpers.R` `585 / 745` covered, `160` uncovered
- `R/func_lipid_s4_da_plot_methods.R` `87 / 208` covered, `121` uncovered

Assessment:

- Lipidomics is the main remaining omics-specific headroom reservoir if a later push targets `91-95%`.
- Current evidence does not show a parity failure, but coverage density is visibly thinner here than in proteomics or metabolomics.

### Integration / Multi-omics

Manifest bundle context:

- Dedicated manifest bundle count in the current `latest-bundles.json`: `0`
- This slice is represented primarily through shared and workflow-oriented coverage rather than a dedicated parity-bundle family in the current closeout manifest

Coverage status:

- Covered: `1,616 / 2,044`
- Coverage: `79.0607%`

What materially moved:

- The multi-omics slice improved through shared enrichment coverage, STRING orchestration/retrieval coverage, and branch-hardening around shared compare suites.
- This area remains the clearest package-level laggard after the repo cleared `90%`.

Representative coverage drivers:

- `tests/testthat/test-multiomics-enrich-shared.R`
- `tests/testthat/test-workflow-server-shared.R`
- `tests/testthat/test-workflow-ui-shared.R`

Residual concentration after closeout:

- `R/func_multiomics_enrich.R` `1,616 / 2,012` covered, `396` uncovered
- `R/func_multiomics_mofa.R` `0 / 32` covered, `32` uncovered

Assessment:

- Integration/multi-omics is the single largest non-general obstacle to pushing significantly beyond `90%`.
- The remaining work is real function-surface coverage, not cleanup noise.

### Transcriptomics

Manifest bundle context:

- Dedicated transcriptomics bundle count in the current `latest-bundles.json`: `0`

Coverage status:

- No transcriptomics-specific extracted file slice was present in the exact-union category map used for this closeout.

Assessment:

- Transcriptomics did not present as a distinct demonolithing parity or coverage driver in this final wave.
- Any transcriptomics stability benefit in this closeout is indirect, via shared/general helper coverage rather than transcriptomics-specific test expansion.

### Shared / General Infrastructure

Manifest bundle context:

- Bundle count: `35`
- Bundle status in completed FA-C4 evidence: no open parity failures surfaced in the shared/general slice

Coverage status:

- Covered: `10,122 / 11,536`
- Coverage: `87.7427%`

What materially moved:

- The final repo-wide push was finished here.
- The decisive late lifts came from shared enrichment coverage, shared file-management/reporting coverage, `loadDependencies()` coverage, general plotting helpers, general PCA/RLE helpers, GitHub management helpers, volcano/Glimma helpers, and file-management configuration/path contracts.

Representative coverage drivers:

- `tests/testthat/test-general-enrichment-shared.R`
- `tests/testthat/test-general-filemgmt-report-helpers-shared.R`
- `tests/testthat/test-general-filemgmt-load-dependencies-shared.R`
- `tests/testthat/test-general-filemgmt-config-direct-shared.R`
- `tests/testthat/test-general-pca-rle-direct-shared.R`
- `tests/testthat/test-general-volcano-glimma-direct-shared.R`
- `tests/testthat/test-github-management-direct-shared.R`

Residual concentration after closeout:

- `R/func_general_enrichment.R` `1,948 / 2,331` covered, `383` uncovered
- `R/func_peptide_qc_imputation.R` `322 / 509` covered, `187` uncovered
- `R/func_general_plotting_volcano_glimma_helpers.R` `210 / 372` covered, `162` uncovered
- `R/func_general_filemgmt_config_helpers.R` `1,196 / 1,318` covered, `122` uncovered
- `R/func_general_helpers.R` `535 / 642` covered, `107` uncovered

Assessment:

- Shared/general infrastructure is where the repo-wide threshold was ultimately closed.
- It is also where a future `91-95%` push would still find meaningful, but not especially cheap, headroom.

## Runtime and Tooling Adjustments Made During Closeout

Non-test changes in the final audit wave were small and targeted:

- `R/func_prot_da_volcano_static.R`
  - qualified `stringr::fixed()` inside the fuzzy-match branch used by `generateProtDAVolcanoStatic()`
- `DESCRIPTION`
  - corrected the local collate ordering between `func_prot_s4_objects.R` and `func_prot_s4_grid.R`
- `tools/refactor/fidelity_audit.py`
  - added `compat` as a recognized shared coverage test family
  - defaulted bundle manifest loading to `.refactor-fidelity-audit/reports/latest-bundles.json` when no manifest path is passed
  - taught the tool to discover shared compare suites by marker comment
  - allowed `surface::manifest::*_legacy` helper bundles to inherit a passing measurement from the matching `*_active` sibling when the same shared suite proved the active bundle
  - expanded inventory/entity and contract-search indexing to better line up shared compare tests with moved helpers
- `tests/testthat/helpers-scoped-mocked-bindings.R`
  - added a shared test helper for safe binding rebinding/restore in namespace-scoped coverage tests

## Residual Risk and Merge Readiness

Merge-readiness conclusion:

- The repo-wide package coverage gate is cleared.
- The demonolithing parity bundle evidence remains cleared from the last completed FA-C4 refresh.
- No open bundle-level comparison gaps or regressions are outstanding in the completed evidence set.

Residual risk concentration:

- The remaining technical debt is coverage density, not a known parity regression.
- The main remaining headroom reservoirs are `func_multiomics_enrich`, `func_general_enrichment`, and several lipid/proteomics plotting/helper files.
- If a later program targets `95%`, the likely worklist starts in:
  - `R/func_multiomics_enrich.R`
  - `R/func_general_enrichment.R`
  - `R/func_lipid_da_plotting.R`
  - `R/func_peptide_qc_imputation.R`
  - `R/mod_metab_norm_server_helpers.R`

## Final Judgment

As of `2026-04-26`, this branch is at a done stage for the current audit objective:

- Demonolithing parity evidence: cleared at the bundle layer
- Repo-wide package coverage: cleared above `90%`
- Remaining work, if desired, is an optional next-phase depth push rather than a blocker for this closeout
