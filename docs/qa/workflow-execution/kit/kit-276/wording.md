# KIT-276 drafts for the owner's wording pass (VAL-41)

The tests assert these as drafts and are edited with the owner's final words. Source: `protNormRuvSkipTexts()` in
`R/mod_prot_norm_workflow_helpers.R`.

1. "Skip RUV" after a covariates or correction run (summary text and notice):
   "RUV-III skipped. The earlier RUV-III result is set aside: the differential analysis will use the correlation-filtered
   table without RUV factors or correction."
   With a differential analysis or enrichment made on the RUV table, a second sentence:
   " Run the differential analysis again."
2. The correlation filter applied again after a RUV step (RUV summary):
   "The correlation filter was applied again, so the table no longer carries the RUV-III result. Compute the RUV-III
   covariates again to use them."
3. The busy refusal of the Skip while the managed runner has a step active:
   "A step is still running. Press \"Compute RUV-III Covariates\" again when it has finished."
4. NEWS entry (`NEWS.md`, `### Fixes`, first bullet, marked DRAFT): "DIA with limpa: what the normalisation tab, the analysis
   record, the study parameters, the session file and the report say about RUV now agrees with what the differential analysis
   uses ..." (full text in `NEWS.md`).
