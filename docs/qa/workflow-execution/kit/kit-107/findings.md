# KIT-107 findings

Append only. Ids are F107-NN. A fix that does not change the methodology is made here, with a test that fails first. Anything else is parked.

| Id | Where | What | Status |
|---|---|---|---|
| F107-01 | `mod_prot_enrich_analysis_helpers.R:537`, `:387` | A second Run on the same tab fails. The first run overwrites the tab's `da_results_data` with `manifest$results`, so it is no longer a manifest, and the second Run falls into the legacy branch: "[error 10] Error in enrichment analysis: Missing contrasts_tbl or design_matrix...". No state or product is written. | Confirmed by W1's probe (`w1/defect-probes.txt`). **Assigned to W3**, fixed with a test that fails first. Not frozen in W1. |
| F107-02 | `mod_prot_enrich_observer_helpers.R:202`, `mod_prot_enrich_state_helpers.R:590` | Enrichment results are never restored after an Open. The tab starts with `all_enrichment_results = list()`, and the restore returns early on anything non-NULL. On an owner holding `enrichment_completed`, a fresh tab shows "No gprofiler2 results available." | Confirmed by W1's probe. **Assigned to W3**, fixed with a test that fails first. W1 freezes only the function-level restore. |
| F107-03 | `createDAResultsForEnrichment` and the run | With a DA results table missing (deleted before the Run), the run completes on a NULL DA table and saves `enrichment_completed` and a product, with only a warning. A product is written for an input that is not there. | **Assigned to W3 (lead decision):** refuse the run with an error notification, writing no state or product. This is a guard, not a method change. The test fails first. W1's frozen `refusals/table_removed` capture is re-frozen once, with attribution. |
| (W1 note) | `analysis_method_display` | A non-numeric taxon makes the method display output error. In the app the error shows in that output, without a crash. W1 records it as today's behaviour. | Recorded, not changed. |
| (P-07) | the organism filter | With a mapping found, the run fails with "object 'clean_acc' not found". | Already parked (P-07); frozen in W1 as today's behaviour. |
