# KIT-249 wording drafts for the owner's pass (VAL-41)

Drafts only: `NEWS.md` on the branch still carries the 82470db1 text until the owner words it. Source:
AUDIT-KIT-249-impl-20261125 I05, I06, I10; ticket section 15a.

| # | Where | Today (82470db1) | Draft |
|---|---|---|---|
| W1 | DA tab notice, `R/mod_prot_da_handlers_run_analysis.R:15-16` | "Some contrasts used the existing BH fallback after q-value estimation failed." | name the comparisons, as `main`'s modal did; drop "existing" (owner question, I06 (a)) |
| W2a | NEWS, KIT-249 entry | "the volcano plot was empty" | "the volcano plot showed no protein and the interactive volcano was replaced by a 'not available' note" |
| W2b | NEWS, KIT-249 entry | "This happens when no p-value in the comparison reaches 0.95, which is likely only when a few tens of proteins or fewer are unchanged, for example in a small protein table or in a comparison in which nearly every protein changes strongly." | "This happens when no p-value in the comparison reaches 0.95. That is likely when only a few tens of proteins or fewer are unchanged, for example in a small protein table or in a comparison in which nearly every protein changes strongly; it can also happen in a larger table when the unchanged proteins' p-values are not spread evenly between 0 and 1 (for example after an uncorrected shift between the groups)." |
| W3 | NEWS, KIT-249 entry, OD-1 (a) | "Analyses recorded before this change now replay with the Benjamini-Hochberg values their record names. Saved results of such comparisons change when the analysis is replayed" | "Analyses recorded before this change now replay with the Benjamini-Hochberg values their record names, so for those comparisons a replay differs from the results saved before this change (significant counts, plots, enrichment and report numbers included)" |
| W4 | `diaInteractiveVolcanoUnavailableText`, `R/utils_workflow_dia_da.R:1184-1188` | "(the FDR could not be estimated, for example on too few proteins)" | on this route the note now appears only when no row has a finite FDR and log2 fold change; the owner words it (I10) |
| W5 (optional, I06 (b)) | NEWS, and the report beside the FDR sentence | none | "Benjamini-Hochberg values are more conservative than Storey's q-values would be here: they assume every protein could be unchanged, which is far from true exactly when this happens." |
