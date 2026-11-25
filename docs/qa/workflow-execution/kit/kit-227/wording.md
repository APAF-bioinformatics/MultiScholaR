# KIT-227: the RUV sentences the templates print

**Draft; the owner sets the final wording.** The owner accepted drafts W1 to W8 of the ticket (section 7) on
2026-11-25 as the basis for the tests (OD-1); the final wording comes in a later wording pass. The templates carry
the sentences below exactly. `{...}` is a recorded value; a value written as 'N/A' (or missing) counts as not
recorded.

## W1. Proteomics (LFQ, TMT, DIA-NN, DIA with limpa without the limpa quantification), METHOD DETAILS

When RUV was applied, the paragraph opens with the two unchanged sentences ("To remove batch effects ... between
different samples."), then:

- Automatic: "For this study, an automatic optimisation process was run, which identified that using {pct}% of
  features as empirical negative controls ({n} proteins) provided the best separation of {groups}. Based on this, {k}
  unwanted components were selected for removal." (Byte-identical to the earlier text for the column `group`.)
- Manual: "For this study, no automatic optimisation was run: the analyst chose to use {pct}% of features as
  empirical negative controls ({n} proteins) and to remove {k} unwanted components." A recorded column other than
  `group` adds " The controls were chosen by comparing {groups}."
- Manual, percentage not recorded: "For this study, no automatic optimisation was run: the analyst set the share of
  features used as empirical negative controls (not recorded; {n} proteins were used) and chose to remove {k} unwanted
  components."
- Applied, mode not recorded: "For this study, {pct}% of features were used as empirical negative controls ({n}
  proteins) and {k} unwanted components were removed; whether these values came from an automatic optimisation or
  were set by the analyst is not recorded."
- Not applied (replaces the paragraph): "RUV-III batch correction was not applied. The analysis proceeded directly
  with normalised data for differential abundance analysis."
- Not recorded (replaces the paragraph): "Whether RUV-III batch correction was applied is not recorded for this
  analysis."

`{groups}`: "biological groups" for `group`; "the sample groups given by the design column '{g}'" for another
column; "the sample groups (the design column used is not recorded)" when unrecorded.

**Added by the implementer (not in W1; for the owner's pass):** the cases R03 requires that the drafts do not cover.
- Automatic, percentage not recorded: "For this study, an automatic optimisation process was run, which selected the
  share of features used as empirical negative controls (not recorded; {n} proteins were used) that provided the best
  separation of {groups}. Based on this, {k} unwanted components were selected for removal."
- Mode and percentage not recorded: "For this study, a share of features (not recorded; {n} proteins were used)
  served as empirical negative controls and {k} unwanted components were removed; whether these values came from an
  automatic optimisation or were set by the analyst is not recorded."
- Count not recorded: "({n} proteins)" reads "(number of proteins not recorded)", and "(not recorded; {n} proteins
  were used)" reads "(not recorded, nor the number of proteins used)". k not recorded: "{k}" reads "an unrecorded
  number of".

## W2. The count sentence (LFQ, TMT, DIA-NN)

"Following RUV-III batch correction and subsequent missing value cleanup, {N} proteins were retained for
differential abundance analysis." Printed only when RUV was applied.

## W3. Results, Figure 1 caption (LFQ, TMT, DIA-NN, DIA with limpa without the limpa quantification, lipidomics, metabolomics)

Applied and not applied: as before. Not recorded: "Figure 1 shows four quality-control views at each stage of
processing drawn for this analysis. Whether RUV-III correction was applied is not recorded for this analysis."
followed by four bullets without panel letters: in proteomics the earlier caption's four bullets; in lipidomics and
metabolomics (whose captions have no such list) the four views of the template's own Quality Control method list
("Principal component analysis (PCA) to visualize sample clustering", "Density boxplots to assess distribution
consistency", "Relative log-expression (RLE) plots to detect systematic biases", "Pearson correlation matrices to
evaluate replicate consistency").

## W4. DIA with limpa (limpa quantification), Results, Figure 1 caption

- Covariates: as before.
- Correction: right column "after RUV-III correction (the corrected matrix entered the differential abundance
  model)".
- RUV not run: "Figure 1 shows four quality-control views at each stage of processing drawn for this analysis.
  RUV-III correction was not applied to this dataset." followed by the four bullets without panel letters.

## W5. Lipidomics / metabolomics, METHOD DETAILS, per-assay RUV paragraph

- `group`: unchanged ("RUV correction was applied independently to each assay using the 'group' variable to preserve
  biological signal. {per-assay sentences}").
- Another column: "RUV correction was applied independently to each assay. Samples with the same value of the design
  column '{g}' were treated as replicates, so the differences between those values were kept and other differences
  between samples could be removed as unwanted variation. {per-assay sentences}"
- Unrecorded: "RUV correction was applied independently to each assay (the design column that defined the replicate
  samples is not recorded). {per-assay sentences}"
- Manual (record `ruv_mode` 'manual'): appends " The number of unwanted components and the share of features used as
  controls were set by the analyst; no automatic optimisation was run."

## W6. Lipidomics / metabolomics, the older single-value form

W1's Automatic, Manual and mode-not-recorded texts (and the implementer's additions above), with "lipids" or
"metabolites" and the template's spelling "optimization"; the Automatic text keeps the template's own form: "For this
study, an automatic optimization process identified that using {pct}% of features as empirical negative controls ({n}
lipids) provided the best separation of {groups}. Based on this, {k} unwanted components were selected for removal."
The mode is the record's only (the file never carries it), so without a record this form prints the
mode-not-recorded text.

## W7. Lipidomics / metabolomics, not applied or not recorded

- Not applied: "RUV-III batch correction was not applied. The analysis proceeded directly with normalized data for
  differential abundance analysis."
- Not recorded: "Whether RUV-III batch correction was applied is not recorded for this analysis."

## W8. Transcriptomics, Methods, "Unwanted variation (RUVs)", first paragraph

- Automatic: unchanged.
- Manual: "Unwanted variation was estimated with RUVs (RUVSeq) on inferred empirical negative control genes: the least
  significant genes of a preliminary group-only quasi-likelihood GLM, taken at a chosen percentage of the filtered
  genes. They are inferred from the data, not known in advance. The selection mode was manual: the analyst set k = {k}
  and {pct} % of the genes, giving {n} inferred empirical negative controls; no grid of percentages was searched." The
  second paragraph is unchanged.
