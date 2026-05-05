# E2E-007 Proteomics Enrichment Backend Matrix

E2E-007 covers the proteomics enrichment GUI branch after differential abundance
analysis. It proves that both supported enrichment paths can be selected from the
browser, executed without live external services, exported, copied into
publication artifacts, and serialized into report-facing workflow parameters.

The goal is wiring fidelity, not biological interpretation. The deterministic
terms are synthetic and intentionally small so the lane can run headlessly in CI
and fail for GUI/state/export regressions rather than network availability.

## Browser Contract

The tracked scenario lives in
`tests/testthat/test-e2e-proteomics-enrichment.R` and runs two full browser lanes:

- `gprofiler2` with taxon `9606`, which is treated as a supported organism;
- `clusterprofiler` with taxon `999999`, which forces the unsupported-organism
  fallback route.

Each lane drives the production app through these handoffs:

- launches `run_app(test_mode = TRUE)` through the shared browser harness;
- selects proteomics and creates an isolated temporary project;
- uploads the canonical DIA report fixture and FASTA through import;
- completes the design builder from the fixture manifest;
- runs DIA peptide/protein QC through the canonical browser helper;
- exports the normalized filtered session;
- reloads that filtered session in differential abundance analysis;
- runs DA and waits for the enrichment tab to become available;
- selects the enrichment backend tab and organism taxon explicitly;
- runs enrichment with deterministic thresholds;
- downloads the enrichment ZIP archive;
- runs the proteomics summary/report workflow;
- asserts publication-copy artifacts and saved study parameters;
- checks browser console errors and Shiny error notifications.

## Required Evidence

The test fails unless the selected backend is observable in all of the following
places:

- the enrichment status output contains `Method: gprofiler2` or
  `Method: clusterprofiler`;
- the pathway-enrichment TSV files under
  `results/proteomics/pathway_enrichment` contain an `analysis_method` column
  matching the selected backend;
- the downloaded enrichment ZIP contains the expected backend-specific TSV:
  `gprofiler2_results.tsv` for `gprofiler2` or
  `clusterProfileR_results.tsv` for `clusterprofiler`;
- the downloaded ZIP does not contain the opposite backend TSV, which catches
  silent branch collapse;
- `enrichment_analysis_summary.txt` records the selected backend;
- `results_summary/proteomics/Publication_tables/Pathway_enrichment_results_proteomics.xlsx`
  exists and carries the selected backend into the workbook data;
- `study_parameters.txt` includes `[Enrichment Analysis UI Parameters]`;
- `study_parameters.txt` records `database_source = <backend>`;
- `study_parameters.txt` records the taxon used to choose that backend;
- browser console errors and visible Shiny error notifications are absent.

## Deterministic Backend Behavior

The browser app runs enrichment in a separate R process, so test-process mocks are
not sufficient for required E2E coverage. E2E-007 adds app-side deterministic
test-mode seams in `R/mod_prot_enrich_server_helpers.R`:

- deterministic UniProt annotations are generated when no annotation table is
  available in `test_mode`;
- deterministic `gprofiler2` results are represented as a list with a
  `result` data frame and gprofiler-like metadata;
- deterministic `clusterprofiler` results are represented by the lightweight S4
  class `ProtEnrichDeterministicResult` with a `result` slot;
- both deterministic branches write the same pathway-enrichment TSV files that
  downstream publication/report code expects;
- both branches populate `EnrichmentResults` with backend-specific data, plots,
  summaries, and method labels.

Production behavior still delegates to `processEnrichments()`. The deterministic
process is selected only through `resolveProtEnrichProcessEnrichmentsFn()` when
`is_test_mode()` is true.

## Regression Surface

E2E-007 specifically guards these touch points:

- the proteomics workflow must mount the enrichment module under the real
  `enrichment_analysis` namespace;
- browser automation must drive the same `organism_taxid`,
  `enrichment_method_tabs`, contrast, cutoff, run, and download controls used by
  users;
- DA output must provide selectable contrasts to enrichment after filtered-session
  reload;
- supported-organism routing must preserve the `gprofiler2` method label;
- unsupported-organism routing must preserve the `clusterprofiler` method label;
- result table aggregation must accept both gprofiler-style lists and S4 objects
  with a `result` slot;
- enrichment exports must remain backend-specific and must not include both branch
  files for one selected backend;
- publication table generation must carry the backend label rather than stripping
  branch metadata;
- summary/report parameter serialization must include the selected taxon,
  thresholds, and backend.

## Focused Runner

Run the ticket suite with:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(e2e-browser-harness|e2e-enrichment-doubles|prot-11-enrichment-module-contracts|e2e-proteomics-enrichment)$",
  reporter = "summary"
)
```

After harness changes, also run the shared proteomics browser regressions:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(e2e-proteomics-dia|e2e-proteomics-dia-limpa)$",
  reporter = "summary"
)
```
