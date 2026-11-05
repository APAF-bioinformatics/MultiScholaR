# KIT-093: the count class, design validation, the notebook stages as kernels and the publication renderers

Ticket revision 1.1. Implemented 2026-11-06 on `performance-v2`. In-process tests only, no counted members.

## What was done

- **KIT-040 addendum (first step, approved).** The driver `tools/profiling/trx_notebook_reference.R` also freezes the publication figures' `layer_data()` and plotted data, the workbook's `Optimal_Parameters`, the Glimma annotation and the summary text (25 added products, 95 in all). Bindings were committed before the counted renders (`8ad5f71e`), and the three references plus a repeat were re-rendered under KIT-040's rules (`b248c28a`). On each corpus the 70 KIT-040 products are 70 of 70 `identical()` before and after. Record: kit-040/outcome.md, "KIT-040 addendum (2026-11-06, for KIT-093)".
- **Class** (`R/func_trx_s4_class.R`): `TranscriptomicsAssayData` in the lipid/metab shape. Validity runs `trxAssayInputs()` (count table, design checks, count policy). `createTranscriptomicsAssayData(count_data, design_matrix, group_levels, roles, args)` stores `args$group_levels`. The declaration entry is in `workflowClassDeclarations()`, and the owner test's expected set gains the class.
- **Design** (`R/func_trx_design.R`, `validateTrxDesign()`): labels attach by run ID, with no positional fallback. Columns are renamed to `Run`/`group`, `group` is a factor with the declared levels, and `count_index` is returned. Seven classed refusals sit under `trx_design_invalid`, including `trx_design_groups_unsupported` for anything other than two groups.
- **Count policy** (`R/func_trx_counts.R`, `validateTrxCounts()`): integer-valued doubles pass. The refusals are `trx_counts_na`, `trx_counts_nonfinite`, `trx_counts_negative`, `trx_counts_fractional` and `trx_counts_zero_library`, all under `trx_counts_invalid`. The policy runs in the class validity and at every kernel that takes counts.
- **Kernels.** `trxExpressionFilter`, `trxTmmLogCpm`, `trxRuvGrid` (with `trxSelectGridWinner`, `grid_failures` and the CCA seam `trxCancor`), `trxRuvs`, `trxDaContrast` with `trxParseContrastFile`, `trxAnnotateDe`, `trxRuviii` (`m_source` group or technical_replicate), and `trxPcaStats`/`trxRleStats`. The preflight refusals are under `trx_error`.
- **Renderers** (`R/func_trx_publication.R`): `trxCcaPanels`, `trxNormStageGrid`, `trxVolcano`, `trxSignificanceBar`, `trxDeTables`, `trxGlimmaInputs`, `trxRuviiiGrid` and `trxSummaryText`. They are the notebook's functions with qualified calls. They return figures with `layers`/`data` and write no files.
- **DESCRIPTION.** Suggests gains edgeR, RUVSeq, EDASeq and BiocGenerics, and the nine files enter `Collate`. `ruv` and `limma` stay in Imports. NAMESPACE is unchanged (nothing exported).

## Evidence

- **Environment** ([environment.txt](environment.txt)): every item matches KIT-040: R 4.6.1, Bioconductor 3.23, reference BLAS/LAPACK 3.12.0, edgeR 4.10.1, limma 3.68.4, ruv 0.9.7.1, RUVSeq 1.46.0, EDASeq 2.46.0, Biobase 2.72.0, BiocGenerics 0.58.1 and the plotting packages. No thread variable is set. [parity-criteria.md](parity-criteria.md) was committed with it (`b352e3c4`) before the first comparison.
- **Test first** ([test-first.txt](test-first.txt)): before the functions existed, every new test errored, and the owner test failed on the class set.
- **Parity, kernels** (`test-trx-kernels.R`, 341 pass): on `trx_screen`, `trx_edge_12` (k cap 3) and `trx_edge_16` (k cap 4), every kernel product's `serialize` sha256 equals its KIT-040 fingerprint, one expectation per reference product. The test also asserts that the mapped products plus the listed notebook intermediates cover the store exactly. The first comparison passed; no correction pass was used for parity.
- **Parity, renderers** (`test-trx-publication.R`, 174 pass): on the three corpora, every renderer product matches: the figures' `layer_data()` per layer and plotted data, the grids' panels in reading order and their annotation, the PCA axis labels, the RLE box statistics, the workbook sheets, the aligned Glimma annotation and the summary text without its timestamp. The DE CSV and the top-genes TSV, written with the notebook's writers, equal the stable output hashes.
- **Determinism:** two consecutive runs on the screen corpus are `identical()`. With the count columns reversed and the ID columns moved, every `trx_edge_12` product still matches the reference. A `group_id` of `condition` also matches every product.
- **Tie and failure:** a constructed tie selects the first member. A forced failure of the 30% member (CCA seam mocked) appears in `grid_failures` as `kit093_forced_failure`. The other five members remain, and the winner and its controls match the reference. When every member fails, the result is `trx_grid_no_member`.
- **Refusals:** every class in the design, counts and preflight tables, plus `trx_ruviii_replicates_invalid`, is raised by `expect_error(class =)`. KIT-040's fixtures (regenerated, sha256 equal to bindings.json) cover NA, negative and fractional values. In-test variants cover NaN, ±Inf, the zero library, a fractional matrix at a kernel entry, an empty filter, fewer than 10 genes, a W confounded with group, unparseable and non-estimable contrasts, and the absent, NA, group-spanning (`replicates`) and singleton technical-replicate columns. KIT-040's whole-double fixture passes.
- **RUVs counts:** `ruvs_counts` is `identical()` to the filtered counts on all three corpora.
- **Technical-replicate RUVIII:** the product is `identical()` to the direct `ruv::RUVIII(Y = t(log_cpm_tmm), M = ruv::replicate.matrix(tech_rep_group), ctl, k, return.info = TRUE)$newY`. It differs from the group default, and its `m_source`/`m_label` are set.
- **Labels:** `control_source` is `inferred_least_significant` on the grid, RUVs and both RUVIII products. `top_table` has `FDR`, equal to `p.adjust(PValue, "BH")`, and no q-value column. The volcano's y axis reads FDR (BH).
- **Targeted files** ([tests.txt](tests.txt), one process per file, kit_040 library, source tree): `test-trx-design.R` 67, `test-trx-kernels.R` 341, `test-trx-publication.R` 174, `test-workflow-omics-owner.R` 74 and `test-workflow-class-declarations.R` 181 pass. All have 0 fail and 0 error, and the trx files show **SKIP 0**. The dependents' test `test-workflow-checkpoint-dia-state.R` has 174 pass and its 4 pre-existing quarantine skips. Each trx file takes 6-38 s against the 900 s ceiling.
- **Dependents** ([dependents.txt](dependents.txt)): `R/utils_workflow_dia_state.R:229,256` now accepts one more assay class (additive). It is covered by the owner and dia-state tests.
- **R CMD check** ([check.log](check.log), `--no-tests --no-manual`, HEAD `6862f662`, kit_040 library): 3 WARNINGs and 3 NOTEs. The items are identical to KIT-040's check, with no new global-variable or undeclared-import item.
- **Install** ([installed-test.txt](installed-test.txt)): `kit_093` = `$P/core/browser-smoke.hBkN8K` (source `6862f662`), installed through the kit runner from kit_040's library. The fifteen KIT-040 packages were carried (RUVSeq, EDASeq and the rest). One run of `test-trx-kernels.R` against the installed build, `MultiScholaR` loaded from that library: 341 pass, 0 fail, 0 error, 0 skip.

## Recorded divergences and decisions

- **SummarizedExperiment:** the class does not inherit it (`conventions/R-bioinformatics.md:28-33`). As for lipid and metab, the kit's declaration is the container contract (TRX-005 decision 1).
- **Temporary attach:** `trxRuvs()` attaches BiocGenerics for the single `RUVSeq::RUVs` call (`withr::with_package`). RUVSeq's SeqExpressionSet method calls `counts()` unqualified and works only when BiocGenerics is attached, which the notebook gets through RUVSeq's Depends. The call is unchanged, and the search path is restored afterwards. This diverges from the conventions' "no library() in package functions".
- **Glimma annotation:** `trxGlimmaInputs()$anno` is aligned to `qlf`'s feature order, as the ticket's renderer table says. The notebook passes it in `top_genes` order, and `glimmaVolcano` binds it to its table by position, so the notebook's Glimma tooltips are misaligned. Declared in parity-criteria.md (item 4).
- **Single surviving gene:** `drop = FALSE` in the filter subsets keeps one surviving gene as a matrix, where the notebook would fail. Results are identical otherwise.
- **Added kernel products:** `ctl_mask`, `pdata` and `pca_axis_labels` are notebook objects. `norm_factors_filtered` holds edgeR's samples table (its `norm.factors` column is the factors), because the reference product is that table.
- **Volcano label:** the y label and subtitle say "FDR (BH)" (ticket, KIT097-F03). Labels are not compared products. The summary text keeps the notebook's wording, because it is a compared product.
- **Open user decision:** more than two groups is refused (`trx_design_groups_unsupported`) and not implemented.

## Limits

- Synthetic corpora only; real data at the final gate. The technical-replicate product has no notebook reference; its oracle is the direct call. With technical replicates in the corpus, the notebook's edgeR treats them as independent samples (pseudo-replication), and parity says nothing about that design.
- Under `R CMD check` the trx tests skip, because `docs/` and `tools/` are not shipped.

## Hand-offs

- **KIT-095:** the kernels and renderers named above, the edgeR fit (`fit`) with `method = "edgeR_QL_robust"`, and `m_source` for the RUVIII export step. The filter's `min_sample_size`/`median_lib_size` are for the partitioned variant. `grid_failures` and `trxGlimmaInputs()` cover the Glimma write and its recorded outcome.
- **KIT-096:** `validateTrxCounts()` and `validateTrxDesign()` with their classes, and `args$group_levels`.
- **KIT-097:** the FDR (BH) rule.

Rollback: `git revert` of this ticket's commits (and of the addendum's, which restores KIT-040's 70-product stores). The full KIT-040 references are kept as `$P/kit/trx-reference/reference-<corpus>.kit040-v1.rds`. Nothing was installed into any library by this ticket beyond the kit install below.
