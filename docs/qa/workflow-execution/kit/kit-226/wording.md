# KIT-226: new user-facing texts (drafts for the owner's wording pass, VAL-41)

Each text below is a draft written by the implementer. The owner sets the final wording; a change of any of them moves
the stored references named.

| # | Where | Before | Draft now | Stored references it moves |
|---|---|---|---|---|
| 1 | g:Profiler results table, archive `gprofiler2_results.tsv`: new column name (`R/utils_workflow_enrichment.R`, builder) | (none) | `p_adjust_method` | W1 table leaves (section 5 (i)), `prot_dia.json` 6 leaves |
| 2 | Value of that column when the answer states no correction (`gprofilerAdjustmentMethod`, `R/func_general_enrichment_engine_helpers.R`) | (none) | `not recorded` | same as 1 (every test-mode table reads it) |
| 3 | Value when stated | (none) | copied verbatim from g:Profiler (not reworded by design, E-1); a live table is expected to read `gSCS` or `g_SCS` for the default and `fdr` or a normalised spelling for BH, whichever the service echoes (not observed) | none of its own |
| 4 | Published `<contrast>_<dir>_enrichment_results.tsv` and the enrichment workbook: new column | (none) | `p_adjust_method` (same values as 2 and 3) | W1 `kernel_gprofiler` `run` `files` leaves of the two pathway TSVs |
| 5 | g:Profiler figure, y-axis title (`buildEnrichmentGgplot`) | `-log10(p-value)` | `-log10(adjusted p-value)` | W1 `kernel_gprofiler` `run` plot leaf (or leaves, `:1842`, `:1847`) |
| 6 | g:Profiler figure, hover text line | `P-value: 0.0012` | `Adjusted p-value: 0.0012` | as 5 |
| 7 | Tab archive `enrichment_analysis_summary.txt`, g:Profiler line (`R/mod_omic_enrich_panels.R`) | `Q-value Cutoff: 0.05` | `Adjusted p-value cutoff: 0.05` | W1 `gprofiler` `download` summary leaf; `test-enrichment-significance.R:497` |
| 8 | NEWS entry (`NEWS.md`, Fixes) | (none) | as written on the branch | none |

**For the wording pass (audit K226-I05):**
- On the g:Profiler route one number now has three names on the tab: the input "Significance cutoff"
  (`R/mod_omic_enrich.R:100`), the status "Q-value cutoff: 0.050" and waiting step "3. Set Q-value cutoff"
  (`R/mod_omic_enrich_panels.R:57`, `:77`, parked by VAL-128), and the archive's "Adjusted p-value cutoff".
- "not recorded" (row 2) shows in every test-mode table, including a demo with test features on; the auditor suggests
  "not stated by g:Profiler".
- The method column names the correction, not that g:Profiler corrects each data source separately (the report says
  so; the table and files do not).
- The NEWS line now says that files written or downloaded before the fix keep the old labels until the enrichment is
  run again.

Kept as they are (column names are clusterProfiler's, data rather than prose): `p.adjust`, `Count`. Not touched: the
status and waiting texts that say "Q-value cutoff" on both routes (`R/mod_omic_enrich_panels.R:57`, `:77`; parked by
VAL-128 for the wording pass).
