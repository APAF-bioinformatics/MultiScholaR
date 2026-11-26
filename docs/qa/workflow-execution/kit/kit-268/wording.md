# KIT-268 wording: drafts for the owner's pass (VAL-41, VAL-312, VAL-314)

The sentence (VAL-314, as the wording draft; one text in the report templates and in the record's `fdr_fallback_reason`):

> For these comparisons the q-value column holds Benjamini-Hochberg adjusted p-values, which equal Storey's q-values with the proportion of unchanged features (pi0) set to 1; because pi0 could not be estimated for these comparisons, the values are larger than Storey's q-values would be by a factor of 1/pi0 of unknown size (1 or more), so they are more conservative than, and not on the same scale as, the q-values of the other comparisons.

Where it is read: the six reports' differential-abundance paragraph, appended after "..., where q-values could not be computed and Benjamini-Hochberg adjusted p-values were used instead." (only when a comparison fell back); `analysis_record.json`, `fdr_fallback_reason` of each fallback comparison, after the existing text. The record's Word and Excel outputs are unchanged (VAL-314 (2)). The record uses the identical sentence, so the two agree; its "the other comparisons" reads against the record's other entries.

The notice (type warning, persistent; one text for proteomics, lipidomics, metabolomics; comparisons grouped by assay):

> Storey's q-values could not be computed for these comparisons, so Benjamini-Hochberg adjusted p-values are shown in the q-value column instead: A_vs_B; C_vs_D in LCMS_Neg.

Code: `daQvalueFallbackSentence`, `daQvalueFallbackText` (R/func_omics_da_shared.R); the six templates hold their own copy of the sentence.
The NEWS sentence is under "Fixes" of the single unreleased heading of NEWS.md.
