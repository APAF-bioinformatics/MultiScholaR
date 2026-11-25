# KIT-190: wording drafts for the owner's pass (VAL-41)

DRAFTS. The owner sets the final wording. None of these changes a number.

1. **Automatic k reduced** (`omicsRuvKBound`, `R/utils_workflow_omics_declarations.R`; a classed warning
   `multischolar_ruv_k_notice` and the record's `ruv_k_notice_per_assay`):
   "{assay}: the automatic search found {k_found} unwanted factors. With {n} samples, at most {cap} can be removed while
   leaving enough samples to judge the differences between groups, so {cap} were removed."
2. **Manual k above the cap** (same function; warning only, the chosen k is kept):
   "{assay}: you chose to remove {k} unwanted factors. With {n} samples, at most {cap} are advised, to leave enough
   samples to judge the differences between groups. Your choice was kept."
3. **Refusals of an unusable grouping** (7.2, unchanged since the first repair): "Control features cannot be chosen for
   {assay}: the design has no column '{grouping}', so the features cannot be compared between groups. Choose a grouping
   column that is in the design." and "Control features cannot be chosen for {assay}: after leaving out the pooled
   samples, fewer than two groups of '{grouping}' remain." with "Choose a grouping column with at least two groups that
   are not pools."
4. **The RUV design panel's legacy note** (7.5): as in `R/utils_omic_gui_config.R`.
