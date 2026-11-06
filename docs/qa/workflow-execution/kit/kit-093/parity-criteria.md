# KIT-093 parity criteria (committed before the first comparison)

## The reference

- The fingerprint stores `docs/qa/workflow-execution/kit/kit-040/reference-<corpus>.rds` of the KIT-040 addendum (2026-11-06), 95 products each. Before comparing, each test checks the store's sha256 against the value KIT-040 recorded ([addendum-093-compare.txt](../kit-040/addendum-093-compare.txt)): `trx_edge_12` `f460b376008219c667fdca32bd1184813ff1f5b8ba8c53d31f45c63b728cb12c`, `trx_edge_16` `ac118c0ada3411ad5510efdd385e0ae08a4b6b3cecfb6c9a4ba436c75a7ecf59`, `trx_screen` `dd058d7bd344e606899f95dac67786ef8a60df473ffede1e71918bdca0d66840`.
- The inputs are regenerated in the test by KIT-040's generator (`tools/profiling/kit_omics_corpus.R`, `trx_screen` with 6, 3 and 4 biological samples per group), and each regenerated `manifest.json` must have the store's `corpus_manifest_sha256` before any kernel runs. They are read as the notebook reads them (rmd:121, 136-140): `vroom::vroom(show_col_types = FALSE)`, `Run` as character, `group` a factor with levels `RPMI`, `Sera`.

## What "identical" means here

- A kernel or renderer product is identical to its reference product when the sha256 of `serialize(x, NULL, version = 3)` equals the store's fingerprint. This is KIT-040's binding of every reference product. Equal serialisations imply `identical()`; the converse is not assumed, so a hash mismatch is a failure, never a pass.
- One expectation per reference product, named by the reference product.
- On a mismatch, the diagnosis (not the test) loads the full product from `$P/kit/trx-reference/reference-<corpus>.rds` and reports the cells and the size of the difference.

## Mapping (kit-040/product-map.tsv)

Products are compared under the names the product map gives. Where the reference product is a shape the notebook built from several objects, the test builds that shape from the kernel products with the driver's own expression (`tools/profiling/trx_notebook_reference.R`), not a new one:

| Reference product | Built from |
|---|---|
| `grid_members` | `grid_members` (the successful members) and `grid_failures` (status `failed`, the other fields `NA`), one row per percentage, the driver's `rbind` |
| `grid_cancor_data`, `grid_control_genes` | each member's `cancor_data` and `control_genes`, named `pct_<percentage>` |
| `best_res` | `best[c("percentage", "n_controls", "best_k", "separation_score", "k_penalty", "composite_score", "ctl_idx", "control_genes")]` |
| `dge_filtered_samples` | `norm_factors_filtered`, which holds edgeR's samples table (`group`, `lib.size`, `norm.factors`) |
| `pca_<stage>_variance` | `pca_variance`: the driver's rule, `var` of each PC column of the plotted scores over their sum |
| `pca_<stage>_axis_labels` | the renderer's PCA panel, `unlist(ggplot2::get_labs(p)[c("x", "y")])` |
| `rle_<stage>_box` | the driver's column subset of `layer_data(p, 1)` of the renderer's RLE panel |
| `rle_<stage>_whiskers` | `rle_whiskers`: the driver's `tapply` rule on the plotted RLE data |
| `norm_grid_layers`, `ruviii_grid_layers` | the renderer's grid: every layer of every panel in reading order (the driver's flattening of the patchwork) |
| `summary_text` | the renderer's lines and grid table written as the notebook writes them (`writeLines`, then `write.table(append = TRUE)`), read back, the `Timestamp: ` line dropped |
| stable CSV/TSV outputs | `trxDeTables()` written with the notebook's writers (`write.csv(row.names = FALSE)`, `vroom::vroom_write`), sha256 against the store's `stable_outputs` |

Figures are compared by `layer_data()` of every layer and by their plotted data frames; plot objects and PNG/PDF/xlsx/HTML bytes are never compared.

## Declared divergences (none changes a notebook product)

1. `grid_failures` is an added product (the notebook drops a failed member silently).
2. Refusals where the notebook would fabricate or crash: the design table, the count policy, the preflight classes, `trx_ruviii_replicates_invalid`.
3. Designs with other than two groups are refused (`trx_design_groups_unsupported`): the notebook's `coef = 2` ranking defines only a two-group test. *Superseded by the KIT-093 addendum (revision 1.2, 2026-11-06): any number of groups, the ranking over every group coefficient (with two groups the notebook's `coef = 2`); the comparisons and their criteria above are unchanged.*
4. `trxGlimmaInputs()$anno` is the notebook's `glimma_anno` reindexed to `qlf`'s feature order, as the ticket's renderer table requires. The notebook passes it in `top_genes` order, while `glimmaVolcano` binds it to the table by position in `lrt`'s order. Its reference is the derived `glimma_anno_aligned`; the notebook's own `glimma_anno` stays frozen and is not a renderer product.
5. Products added to a kernel's output that are notebook objects, so that the next kernel or the renderer can be fed without recomputation: `ctl_mask` (trxRuviii, rmd:707-708), `pdata` (trxRuvs, `pData(set_ruvs)`, rmd:414) and `pca_axis_labels` (trxPcaStats, rmd:474, 483-484).

## Tolerances

None. Every comparison is exact. A non-identical result from an unchanged call is a defect, never a tolerance.

## Not compared

The summary file's timestamp line; the technical-replicate RUVIII product, which has no notebook reference (its oracle is the direct `ruv::RUVIII` call of the acceptance criterion).
