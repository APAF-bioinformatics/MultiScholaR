# KIT-272 wording items for the owner's pass (VAL-41)

Drafts only; the substance is decided (VAL-332). The tests hold the sentences in file-local constants, so a later
wording is a text edit of the templates and the constants together (AUDIT-KIT-272-impl-20261126 I03, I08).

## Drafts to confirm

| # | Where | Draft |
|---|---|---|
| D1 | Results, three-stage bullet (panels d-f) and its two variants (c-d; not recorded, no letters) | "Principal component scores by group ..." (ticket section 5, A') |
| D2 | PCA bullet, three variants | "samples plotted on the first two principal components of the log2 protein values (centred, not scaled), computed at each stage from the tenth of the proteins with the largest ratio of standard deviation to mean of their log2 values" |
| D3 | Methods phrase | "box plots of the principal component scores" |
| D4 | Placeholder label (`kinds`) | "principal component scores by group" |
| D5 | Lipidomics and metabolomics bullet (Methods and RUV-not-recorded) | "Box plots of the samples' scores on the first two principal components of every lipid (metabolite) with at least two values, a separate analysis from the PCA panels, which use the tenth of the lipids (metabolites) with the largest ratio of standard deviation to mean" |

## Added by the implementation, not in the ticket's section 5 item 5

- The lipid/metab drafts add "of the lipids" / "of the metabolites" after "the tenth" (the ticket has "the tenth with the
  largest ratio"). It clarifies the sentence and does not alter its substance.

## Not changed here; for the owner

- The lipid/metab Figure 1 caption keeps the word "density": "density boxplots of the samples' scores on the first two
  principal components" (`inst/reports/lipidomics/lipidomics_report.rmd:796`,
  `inst/reports/metabolomics/metabolomics_report.rmd:783`). It is imprecise, not false: no density is estimated. Suggested
  wording: "box plots of the samples' scores on the first two principal components (a separate analysis from the
  principal component analysis panels)". Not pinned by any test.
- Optional: the "Color by" bracket after "group" in D1 (ticket section 5).
- Optional: "within 1.5 times the interquartile range of the box" can be read as 1.5 times the box's length; ggplot2's
  rule is "no further than 1.5 times the interquartile range beyond the box".
- The older rows this replaces: `kit-236/wording.md:36` (P18, "Density boxplots to assess distribution consistency") and
  `kit-259/wording.md` W8 ("the reports call it 'Density plots'") are stale since KIT-272 (VAL-332); not edited here.
