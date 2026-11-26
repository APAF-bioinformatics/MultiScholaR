# KIT-266 user-facing text, as drafts for the wording pass (VAL-41)

The substance of each is decided (VAL-316); the words are the owner's to set. `\n` is a line break inside an image title.
Where the code differs from the ticket's draft this is said.

| Id | Where | Text |
|---|---|---|
| D1a | Figure 1 image title, left column (`diaLimpaCompositeColumnLabels`) | `DPC-Quant estimates\n(log2)` |
| D1b | centre, no between-sample method | `After the peptide count filter\n(no between-sample\nnormalisation applied)` |
| D1c | centre, a method (cyclic loess, quantile, scale; any other name as recorded) | `After <name> normalisation\nand the peptide count filter`. The cyclic loess form has a first line of 32 characters, the ticket's limit is 30; not rendered by R3 (its lane uses method none) |
| D1d | right column, covariates mode | `RUV-III-adjusted\n(for quality control only;\nnot the data analysed)` |
| D1e | right column, correction mode | `RUV-III-corrected\n(the data analysed)` |
| D1f | the limpa `QC_metrics` chunk's `fig.cap` (limpa template only) | `Figure 1: QC metrics across normalisation stages. On DIA with limpa the left and centre columns are before the sample-correlation filter and the right column is after it.` Observation: with a two-column figure there is no right column; the sentence is the ticket's draft |
| D2a | Results, three stages, covariates | `Figure 1 shows four quality-control views at three stages of processing: DPC-Quant's protein estimates (a complete log2 matrix; model-based where peptides were not detected) (left column: panels a, d, g, j); after the peptide count filter, with no between-sample normalisation applied, before the sample-correlation filter (centre column: b, e, h, k); and for the RUV-III-adjusted matrix after the sample-correlation filter (drawn for quality control only; the differential abundance model carried the RUV-III factors as covariates) (right column: c, f, i, l).` With a recorded method: `after the <name> normalisation and the peptide count filter` replaces `after the peptide count filter, with no between-sample normalisation applied` |
| D2b | Results, three stages, correction | the same with `and after RUV-III correction (the corrected matrix entered the differential abundance model), applied after the sample-correlation filter (right column: c, f, i, l).` |
| D2c | Results, two stages (RUV skipped, no covariates step, a record without the factor count, k = 0) | `Figure 1 shows four quality-control views at two stages of processing: DPC-Quant's protein estimates (a complete log2 matrix; model-based where peptides were not detected) (left column: panels a, c, e, g) and after the peptide count filter, with no between-sample normalisation applied, before the sample-correlation filter (right column: b, d, f, h). RUV-III correction was not applied to this dataset.` The bullets d to f (the density panels) are KIT-272's and untouched |
| D3 | note in a left cell whose figure is from another table (`writeDiaLimpaPreNormalisationNotes`) | `DPC-Quant estimates <PCA, RLE, PCA density, Pearson correlation>: could not be drawn.` then `the figures saved in this session for this column come from a later table, or no saved table from before the between-sample normalisation step is available.` |
| D4 | KIT-259's placeholder sentence, limpa template (`QC_metrics_placeholders`) | stage words: `of DPC-Quant's protein estimates`, `after the peptide count filter and any between-sample normalisation`, `of the RUV-III-adjusted matrix after the sample-correlation filter`; panel m: `RUV-III canonical correlation` (also for `ruv_covariates_cancor`) |
| D5 | the covariates step's summary text (`finishProtNormRuvCovariates`) | `... The quantified matrix is unchanged; the adjusted matrix is shown in the RUV-III-adjusted QC figures only.` |
| D6 | NEWS entry | under `### Fixes`, as the ticket's section 13 draft |
