# E2E-009 Lipidomics Multi-Assay Browser Matrix

E2E-009 covers the canonical lipidomics GUI workflow as a true multi-assay
session. The test proves that LipidSearch data can be loaded through the browser,
partitioned into named assays, assigned to one shared design, finalized through
QC, normalized, exported as a filtered session, reloaded into differential
abundance, analyzed per assay and combined, copied to publication outputs, and
rendered through the summary/report path.

The target is app wiring and data-integrity fidelity rather than biological
interpretation. The fixture is synthetic, deterministic, and deliberately small
so failures isolate GUI, state, assay-partitioning, export, DA, or report
regressions.

## Covered Lanes

| Scenario | Fixture inputs | Assays | Contract |
|---|---|---|---|
| Canonical full workflow | `lipidsearch_lcms_pos.txt`, `lipidsearch_lcms_neg.txt` | `LCMS_Pos`, `LCMS_Neg` | Two lipid assays remain partitioned through import, design, QC, normalization export, DA reload, DA output, and report artifacts |
| GCMS import smoke | `lipidsearch_lcms_pos.txt`, `lipidsearch_gcms.txt` | `LCMS_Pos`, `GCMS` | A non-LC assay name imports through the same browser contract without positive/negative-mode assumptions |

Both scenarios use the canonical six-sample `WT` versus `KO` design with the
friendly contrast `KO_vs_WT` and the raw limma contrast `groupKO-groupWT`.

## Browser Contract

The tracked scenario lives in
`tests/testthat/test-e2e-lipidomics-canonical.R`. The full lane drives the
production app through these touch points:

- launches `run_app(test_mode = TRUE)` through the shared `shinytest2` harness;
- selects lipidomics and creates an isolated temporary project;
- uploads primary and secondary LipidSearch files through browser-safe file inputs;
- configures assay names, vendor format, and name sanitization through the real import module;
- waits for setup/import to complete and validates the app state digest;
- completes the lipidomics design builder from `tests/testdata/e2e/manifest.json`;
- asserts design files including `design_matrix.tab`, `contrast_strings.tab`,
  `assay_manifest.txt`, `column_mapping.json`, and per-assay `data_cln_*.tab`;
- finalizes QC through the production S4 finalization handoff;
- runs normalization with ITSD disabled, no between-sample normalization, no
  RUV, and skipped correlation filtering;
- validates `lipid_filtered_session_data_latest.rds` and
  `lipid_filtered_session_summary.txt`;
- reloads the filtered session in differential abundance;
- runs lipidomics DA and asserts per-assay output tables plus combined/per-assay
  dropdown availability;
- saves workflow arguments, copies publication artifacts, generates the report,
  and downloads the rendered report;
- checks browser console errors and visible Shiny error notifications.

## Required Evidence

The E2E test fails unless the lane proves the following invariants:

- import state reports `detected_format = "lipidsearch"`;
- imported assay count and assay names match the manifest;
- imported sample count matches the lane manifest;
- `column_mapping.json` records `LipidName`, `LipidClass`, and all six sample columns;
- `assay_manifest.txt` exactly reflects `LCMS_Pos` and `LCMS_Neg`;
- every expected assay has a corresponding cleaned assay table;
- normalized session export records `omic_type = "lipidomics"`;
- normalized session export contains `current_s4_object`;
- normalized session export preserves `assay_names` and `feature_counts`;
- normalization records `normalization_complete = TRUE`;
- skipped correlation filtering records `correlation_filtering_complete = TRUE`;
- normalization records `normalization_method = "none"`;
- RUV records `ruv_mode = "skip"`;
- ITSD normalization records `itsd_applied = FALSE`;
- DA output includes TSV and Excel files matching
  `de_<assay-prefix>_lipids_<contrast>_long_annot.*`;
- DA output preserves both `comparison = groupKO-groupWT` and
  `friendly_name = KO_vs_WT`;
- DA output contains both assay identities in the `assay` column;
- expected assay prefixes are present: `posmode` and `negmode`;
- `study_parameters.txt` records lipidomics, assay names, the friendly contrast,
  and the raw limma contrast;
- `results_summary/lipidomics/Publication_tables/DA_results_lipidomics.xlsx`
  exists and is non-empty;
- `session_state_*.RDS` is exported from the summary path.

## Production Regressions Caught

E2E-009 locks down two lipidomics-specific wiring risks:

- the second lipid assay import path previously always used the MS-DIAL reader,
  even when the primary file was detected as LipidSearch; assay two now routes
  through the detected/vendor reader and forwards the active column mapping;
- skipped lipid correlation filtering previously advanced the workflow without
  marking `correlation_filtering_complete` or retaining the current S4 object as
  the filtered handoff; the exported DA session now records the same completion
  contract as the explicit correlation-filter path.

These are production fixes, not test-only bypasses. Unit and characterization
tests cover helper-level behavior, while the browser lane proves those helpers
remain wired into the GUI.

## Fixture Contract

The canonical lipidomics fixtures use a narrow LipidSearch-style wide contract:

- `LipidName` is the lipid identifier column;
- `LipidClass` is the annotation column;
- quantitative sample columns are `WT_1`, `WT_2`, `WT_3`, `KO_1`, `KO_2`, and `KO_3`;
- the design file maps those samples to `WT` and `KO`;
- LCMS positive and negative mode are imported as separate browser uploads;
- a GCMS-named secondary fixture exercises non-LC assay naming at import time.

`tests/testthat/helper-e2e-fixture-integrity.R` validates the LipidSearch fixture
shape and the optional `assay2_file` manifest field so browser tests fail early
when fixture drift invalidates the scenario.

## Focused Runner

Run the full ticket slice with:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(e2e-manifest|e2e-fixture-integrity|e2e-browser-harness|workflow-server-shared|lipid-00b-import-module-seams|lipid-00c-import-module-characterization|lipid-12-norm-module-contracts|lipid-12a-norm-module-characterization|lipid-13a-design-module-characterization|lipid-13b-design-wrapper-characterization|lipid-14a-summary-module-characterization|general-filemgmt-config-direct-shared|general-filemgmt-config-roundtrip-contracts|general-filemgmt-s4-contracts|e2e-lipidomics-canonical)$",
  reporter = "summary"
)
```

For the browser-only acceptance gate:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(e2e-lipidomics-canonical)$",
  reporter = "summary"
)
```

## CI Notes

This matrix belongs in the full GUI/nightly/release-candidate job. It is heavier
than a minimal PR smoke because it launches a full browser workflow, writes DA
tables and plots, and renders a report. If CI runtime becomes constrained, keep
the GCMS import smoke as an early wiring sentinel and reserve the full LCMS
positive plus negative workflow for release gating.
