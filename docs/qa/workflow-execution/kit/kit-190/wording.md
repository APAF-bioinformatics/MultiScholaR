# KIT-190: wording drafts for the owner's pass (VAL-41)

DRAFTS. The owner sets the final wording. None of these changes a number.

1. **Automatic k reduced** (`omicsRuvKBound`, `R/utils_workflow_omics_declarations.R`; a classed warning
   `multischolar_ruv_k_notice`, the record's `ruv_k_notice_per_assay`, and since I01 a log line and a warning
   notification in the tab on both routes): the audit's O3 text,
   "{assay}: the automatic search found {k_found} unwanted factors. With {n} samples the application removes at most
   {cap} (the limit proteomics uses), so {cap} were removed."
   The earlier draft ("… while leaving enough samples to judge the differences between groups …") is withdrawn: the
   cap does not look at the groups (I05), so it must not ship.
2. **Manual k above the cap** (same function and delivery; the chosen k is kept):
   "{assay}: you chose to remove {k} unwanted factors. With {n} samples the application removes at most {cap} in
   Automatic mode (the limit proteomics uses); your choice was kept."
   The earlier draft ("… at most {cap} are advised, to leave enough samples to judge …") is withdrawn for the same
   reason.
2a. **The tab's RUV summary where the cap reduced k** (`omicNormRuvSummary`, `R/mod_omic_norm_panels.R`), first line:
   "Factors removed: {k_removed} (the search found {k_found}; at most {cap} with {n} samples)". Where the cap did not
   reduce k the summary is unchanged ("Best k: {k}").
2b. **The session-export summary** (`omicsRuvExportSummaryLine`, used by `R/utils_workflow_saved_exports.R`), where the
   cap reduced k: "{assay}: k={k_removed} (found {k_found}, at most {cap} with {n} samples), %=…, controls=…";
   otherwise unchanged.
3. **Refusals of an unusable grouping** (7.2, unchanged since the first repair): "Control features cannot be chosen for
   {assay}: the design has no column '{grouping}', so the features cannot be compared between groups. Choose a grouping
   column that is in the design." and "Control features cannot be chosen for {assay}: after leaving out the pooled
   samples, fewer than two groups of '{grouping}' remain." with "Choose a grouping column with at least two groups that
   are not pools."
4. **The RUV design panel's legacy note** (7.5): as in `R/utils_omic_gui_config.R`.
