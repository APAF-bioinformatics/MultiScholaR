# E2E-KIT101 Transcriptomics Screen Browser Lane

E2E-KIT101 covers the transcriptomics GUI workflow as a manifest lane, in the form of
the proteomics, metabolomics and lipidomics lanes. The test proves that a count table
can be loaded through the browser, assigned to a design, filtered and finalised through
QC, normalised (TMM/logCPM, the RUV grid and RUVs), analysed with the edgeR
quasi-likelihood DA, exported as the RUVIII integration matrix, copied to publication
outputs, rendered through the summary/report path and exported as a session.

The target is app wiring and data-integrity fidelity rather than biological
interpretation. The fixture is synthetic, deterministic and small, so failures isolate
GUI, state, normalisation, DA, export or report regressions. The KIT-097 journey
(`tests/testthat/test-e2e-trx-journey.R`, the `trx_edge_12` corpus) remains the
transcriptomics proof on the KIT-040 reference corpus; this lane is the manifest-driven
lane that the impact routing and the CI runner select.

## Covered Lane

| Scenario | Fixture inputs | Assays | Contract |
|---|---|---|---|
| Screen full workflow | `trx_screen/counts.tsv` | `counts` | One count assay through import, design, QC, normalisation (RUVs), DA (edgeR QL), the RUVIII export, the publication copy, the report and the session export |

The design is the fixture's twelve runs (three biological samples per group, two
technical replicates each) in the groups `RPMI` and `Sera`, with the friendly contrast
`Sera_vs_RPMI`.

## Browser Contract

The tracked scenario lives in `tests/testthat/test-e2e-transcriptomics-screen.R`. The
lane drives the production app through these touch points:

- launches `run_app(test_mode = TRUE)` through the shared `shinytest2` harness;
- selects transcriptomics and creates an isolated temporary project;
- uploads the count table through the browser-safe file input (`assay1_file_std`) and
  processes the import;
- completes the design builder from the lane's `design.tsv` (factors, metadata per group,
  the manifest's contrast) and saves it;
- applies the expression filter and finalises QC;
- runs the normalisation (`trx-norm-run`) and the DA (`trx-da-run`);
- runs the RUVIII export and waits for `M = group`;
- saves the workflow arguments, copies the publication artifacts, seeds the lane's report
  template stub, generates and downloads the report, and exports the session state;
- checks for visible Shiny error notifications.

## Required Evidence

The E2E test fails unless the lane proves the following invariants:

- `analysis_record.json` records the `trx_expression_filter`, `trx_qc_finalize` and
  `trx_ruvs` steps;
- the DA's source state is `trx_ruvs` and its model is `edgeR_QL`;
- the RUVIII export records `ruviii_m_source = "group"`;
- every file the manifest's `expected_exports` names exists:
  `DA_results_transcriptomics.xlsx`, `transcriptomics_normalized_data.tsv` and
  `ruviii_before_and_after.png`;
- the downloaded report is non-empty and a `session_state_*.RDS` was exported;
- the lane completes within 1,800 s.

## Fixture Contract

`tests/testdata/e2e/trx_screen/` is written by `tests/testdata/e2e/generate_fixtures.R`
through the KIT-040 corpus tool (`tools/profiling/kit_omics_corpus.R`,
`kitOmicsRecipe("trx_screen", 3, 300)`, seed 20261106). It holds `counts.tsv` (`GeneID`,
`GeneName` and one column per run, in a seeded permutation of the design order),
`design_matrix.tab`, `contrast_strings.tab`, the kit `design.tsv`, `truth.json` and the
corpus `manifest.json`. The lane's manifest entry records the recipe and the sha256 of
every file; `test-e2e-manifest.R` checks them. The files are byte-identical to the module
CI fixture `tests/testdata/mci/fx/transcriptomics/trx_screen/`.

## Focused Runner

For the browser-only acceptance gate (a Chromium-based browser is required; set
`CHROMOTE_CHROME` when it is not on the path):

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(e2e-transcriptomics-screen)$",
  reporter = "summary"
)
```

Or through the CI runner: `Rscript tools/ci/run-e2e-ci.R --lane trx_screen`.

## CI Notes

The lane belongs in the full GUI/nightly/release-candidate job beside the other omics'
lanes. The app process loads MultiScholaR from the installed library, so the lane runs
against the install.
