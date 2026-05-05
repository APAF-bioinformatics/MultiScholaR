# E2E-008 Metabolomics LC/GC Browser Matrix

E2E-008 covers the metabolomics GUI workflow across LC, GC, and combined LC+GC
assay layouts. The test proves that metabolomics can be loaded through the
browser, assigned to a shared design, finalized through QC, normalized, exported
as a filtered session, reloaded into differential abundance, analyzed, copied to
publication outputs, and rendered through the summary/report path without assay
state collapse.

The intent is wiring and data-integrity fidelity. The fixture data are synthetic
and deliberately small so failures identify GUI, state, export, or report
regressions rather than biological interpretation.

## Covered Lanes

| Lane | Fixture inputs | Assays | Contract |
|---|---|---|---|
| `metab_lc` | `metab_lc/lcms_pos_features.tsv`, `metab_lc/lcms_neg_features.tsv` | `LCMS_Pos`, `LCMS_Neg` | Two LC assays remain partitioned through import, normalization export, DA, and report artifacts |
| `metab_gc` | `metab_gc/gcms_features.tsv` | `GCMS` | Single GC assay follows the same end-to-end path without requiring a second assay |
| `metab_combined` | `metab_combined/combined_lcms_features.tsv`, `metab_combined/combined_gcms_features.tsv` | `LCMS_Pos`, `GCMS` | Combined LC+GC state preserves assay identity, file prefixes, feature counts, and report metadata |

All three lanes use the canonical `WT` versus `KO` design with the friendly
contrast `KO_vs_WT` and the raw limma contrast `groupKO-groupWT`.

## Browser Contract

The tracked scenario lives in
`tests/testthat/test-e2e-metabolomics-lc-gc.R`. Each lane drives the production
app through these touch points:

- launches `run_app(test_mode = TRUE)` through the shared `shinytest2` harness;
- selects metabolomics and creates an isolated temporary project;
- uploads one or two custom assay files through browser-safe file inputs;
- configures assay names, vendor format, custom metabolite ID mapping, sample
  column pattern, and name sanitization through the real import module;
- waits for setup/import to complete and validates the app state digest;
- completes the metabolomics design builder from `tests/testdata/e2e/manifest.json`;
- asserts design files including `design_matrix.tab`, `contrast_strings.tab`,
  `assay_manifest.txt`, `column_mapping.json`, and per-assay `data_cln_*.tab`;
- finalizes QC through the production QC handoff;
- runs the normalization export path with no between-sample normalization, no
  RUV, and skipped correlation filtering;
- validates `metab_filtered_session_data_latest.rds` and
  `metab_filtered_session_summary.txt`;
- reloads the filtered session in differential abundance;
- runs metabolomics DA and asserts assay-specific output tables;
- saves workflow arguments, copies publication artifacts, generates the report,
  and downloads the rendered report;
- checks browser console errors and visible Shiny error notifications.

## Required Evidence

The E2E test fails unless every lane proves the following invariants:

- import state reports `detected_format = "custom"`;
- imported assay count and assay names match the manifest;
- imported sample count matches the lane manifest;
- design artifacts are written under the metabolomics project source directory;
- `assay_manifest.txt` exactly reflects the expected assay set;
- every expected assay has a corresponding cleaned assay table;
- normalized session export records `omic_type = "metabolomics"`;
- normalized session export contains `current_s4_object`;
- normalized session export preserves `assay_names` and `feature_counts`;
- normalization records `normalization_complete = TRUE`;
- skipped correlation filtering records `correlation_filtering_complete = TRUE`;
- normalization records `normalization_method = "none"`;
- RUV records `ruv_mode = "skip"`;
- ITSD normalization records `itsd_applied = FALSE`;
- DA output includes TSV and Excel files matching
  `de_<assay-prefix>_metabolites_<contrast>_long_annot.*`;
- DA output preserves both `comparison = groupKO-groupWT` and
  `friendly_name = KO_vs_WT`;
- expected assay prefixes are present: `posmode`, `negmode`, or sanitized assay
  names such as `gcms`;
- `study_parameters.txt` records metabolomics, the friendly contrast, and the
  raw limma contrast;
- `results_summary/metabolomics/Publication_tables/DA_results_metabolomics.xlsx`
  exists and is non-empty;
- `session_state_*.RDS` is exported from the summary path.

## Production Regressions Caught

E2E-008 exposed and locked down several real wiring defects:

- metabolomics and lipidomics workflow gate observers could downgrade completed
  steps back to pending after a late upstream observer fired; the gate observers
  are now monotonic;
- metabolomics skipped correlation filtering dispatched the next state without
  marking `correlation_filtering_complete` or retaining the filtered S4 object;
- custom metabolomics imports with no annotation column crashed DA because
  `annotation_id_column = NA` produced a zero-length data-frame column;
- study parameter serialization dropped friendly contrast names even though the
  design builder and DA outputs retained them.

These are production fixes, not test-only bypasses. Unit and characterization
tests cover the helper-level behavior, while the browser matrix proves the
helpers remain wired into the GUI.

## Fixture Contract

The custom metabolomics E2E fixtures use a deliberately narrow tabular contract:

- `Feature.Name` is the metabolite identifier column;
- sample columns are `WT_1`, `WT_2`, `WT_3`, `KO_1`, `KO_2`, and `KO_3`;
- the design file maps those samples to `WT` and `KO`;
- LC fixtures split positive and negative mode into separate assay uploads;
- GC fixtures exercise the single-assay branch;
- combined fixtures exercise mixed assay names in one metabolomics session.

`tests/testthat/helper-e2e-fixture-integrity.R` validates the `custom` fixture
shape and the optional `assay2_file` manifest field so browser tests fail early
when fixture drift invalidates the scenario.

## Focused Runner

Run the full ticket suite with:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(e2e-manifest|e2e-fixture-integrity|e2e-browser-harness|workflow-server-shared|metab-da-orchestration-direct-shared|metab-03c-norm-module-seams|metab-03d-norm-module-characterization|general-filemgmt-config-direct-shared|general-filemgmt-config-roundtrip-contracts|metab-05-summary-module-characterization|e2e-metabolomics-lc-gc)$",
  reporter = "summary"
)
```

For the browser-only acceptance gate:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(e2e-metabolomics-lc-gc)$",
  reporter = "summary"
)
```

## CI Notes

This matrix is appropriate for the nightly/full GUI job and for release-candidate
promotion. It is heavier than a minimal PR smoke because it launches three full
browser sessions and renders reports. If CI runtime becomes a constraint, keep
`metab_gc` as the single-assay smoke lane and retain all three lanes for the full
matrix; do not drop the combined lane from release gating because it is the
assay-state collision sentinel.
