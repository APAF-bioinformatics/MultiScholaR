# KIT-119: enrichment stores no plot objects; input identities; id stripping; exact Reactome organism (2026-11-15)

Scope: the executive calls of 2026-11-15, rows 8 and 9 ([../parked-20261113.md](../parked-20261113.md)). Implemented in a separate worktree and brought onto `performance-v2` as four commits. The full suite, `R CMD check` and the e2e lanes are KIT-118's exit run, which includes these commits ([../kit-118/exit/](../kit-118/exit/)).

| | Change | Commit |
|---|---|---|
| 1 | A run stores no plot objects. One builder (`buildEnrichmentGgplot`, `buildEnrichmentPlot`, `resolveEnrichmentPlotly`; `R/func_general_enrichment_engine_helpers.R`) draws a plot from the stored result. The run writes its HTML and PNG with it and the tab renders with it. The jitter and the label placement have a fixed seed. The S4 class is unchanged; its two plot slots stay empty on the real paths. | `2c7c07de` |
| 2 | A107-06: the enrichment step's outcomes record `annotation_sha256` and `da_inputs_sha256` (one per contrast). | `1473c495` |
| 3 | The g:Profiler custom background ids are stripped as the query ids are (`stripEnrichmentIdSuffix`). | `d3ed1ed3` |
| 4 | `runReactomeEnrichment` matches the organism name in full, ignoring case, and stops with a classed error listing the supported names on an unknown or ambiguous one. | `47e45199` |

## Stored object size (one contrast, both directions)

| Terms per direction | | In memory | Serialized |
|---|---|---|---|
| 20 | with plots | 34.4 MB | 62.7 MB |
| 20 | without | 0.02 MB | 0.005 MB |
| 500 | with plots | 35.1 MB | 65.8 MB |
| 500 | without | 0.29 MB | 0.08 MB |

This is a fixture measurement, not a journey: no journey drives the enrichment tab.

## Choices made differently from the plan

- **Plot inputs as an attribute, not columns.** The clusterProfiler plot needs each term's source and name and the q cutoff. They are carried as the attribute `enrichment_plot_inputs` on the stored `enrichResult`, because the tab copies `@result` whole into its results table and downloads, where extra columns would have shown. The attribute survives the state save and restore (tested).
- **Hashes of the tables as read, not of the files.** `da_inputs_sha256` is the digest of each DA table as the step read it (before the organism filter), by the existing `diaAnnotationDigest`. It avoids duplicating the reader's file-name rule. `sha256sum` of the TSV on disk will not equal the recorded value. The replay does not check the hashes.

## Tests (per file, before -> after, no failure in any)

`test-general-enrichment-shared.R` 104 -> 136; `test-multiomics-enrich-shared.R` 91 -> 108; `test-workflow-class-declarations.R` 181 -> 188; `test-workflow-analysis-record.R` 292 -> 299; unchanged counts: `test-omics-enrich-core.R` 241, `test-omic-enrich-module.R` 1070, `test-module-ci-prot-enrich.R` 79, `test-e2e-enrichment-doubles.R` 34, `test-workflow-dia-products.R` 189, `test-module-ci-prot-da.R` 250, `test-omic-design-module.R` 160, `test-omic-summary-module.R` 1703, `test-omic-gui-config.R` 1906, `test-gui-entry-proteomics.R` 193.

Re-frozen in `test-omic-enrich-module.R`, each with its attribution comment:

- Item 1: 12 values in the two live-kernel cases (the stored state, the bus value, the returned object, the store-entry count, the two plot outputs' text, and the KIT-113 captures' plot digests). The deterministic-kernel cases did not move.
- Item 2: the seven KIT-113 state captures, `config$outcomes` only.

## Left as it is

- The plot is rebuilt on every render of the tab; that cost was not measured.
- A plot that fails to build shows "Plot error: ..." in the tab.
- The clusterProfiler path strips neither query nor background ids; they agree with each other.
- `man/EnrichmentResults-class.Rd` still says the class stores plots. `kit-001/step-map.json`'s enrichment row cites line numbers that have moved (harness-bound file).
- The integration workbook passes an organism Reactome does not have; it now stops (parked row 16).
