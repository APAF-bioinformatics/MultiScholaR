# KIT-093 remediation, 2026-11-06

**Mode:** `audit_and_repair`, on the user's "fix all recommendations" (2026-11-06). The audit report ([KIT-093-AUDIT.md](KIT-093-AUDIT.md)) is unchanged. The only edits are `tickets/kit/KIT-093.md` (revision 1.0 → 1.1) and this file. No code, other ticket, index or plan was touched, and nothing was installed.

**Validation:** `next-openai-ticket.sh KIT-093` passes (`validation.valid: true`, `openai_validation.valid: true`, status consistent, 10 acceptance items).

## Frontmatter changes, for the index sync

| Field | Before | After |
|---|---|---|
| `title` | "Transcriptomics re-cut, part 2: the count class, design validation and the notebook stages as package kernels, identical to the reference" | "Transcriptomics re-cut, part 2: the count class, design validation, the notebook stages as package kernels and publication renderers, identical to the reference" |
| `time_estimate` | 8h | 12h (the renderers are added, KIT095-F02) |
| `ticket_revision` | 1.0 | 1.1 |
| `files_to_modify` | — | adds `R/func_trx_counts.R`, `R/func_trx_publication.R`, `tests/testthat/test-trx-publication.R` and `tests/testthat/test-workflow-omics-owner.R` |
| `dependencies` | KIT-040 | unchanged (KIT-040 also carries the install, user decision 1) |

The index title currently reads "…count class, design validation, kernels identical to the reference", so it needs a sync.

## Closure per finding

| ID | Severity | Repair in KIT-093 rev 1.1, or route | State |
|---|---|---|---|
| KIT093-F01 | blocking | User decision 1 applied. KIT-040 installs RUVSeq and EDASeq into the KIT library: isolated, no upgrades, versions recorded. KIT-093's "Environment" section makes that library a precondition: every test run uses it, `kit-093/environment.txt` must match KIT-040's record, and the ticket stops if the library is absent. **Routed:** the install and the version record belong to KIT-040 (KIT040-F01). | closed in the ticket; execution waits for KIT-040 |
| KIT093-F02 | should-fix | `ruv` and `limma` stay in Imports (`DESCRIPTION:522`, `NAMESPACE:816`). Only edgeR, RUVSeq, EDASeq and BiocGenerics go into Suggests; Biobase is already there. `Biobase::pData` and `BiocGenerics::counts` are guarded and qualified. A backward-compatibility item checks that `NAMESPACE:816` is unchanged. | closed |
| KIT093-F03 | should-fix | The filter section defines the pre-pass as `lib.size × TMM norm.factors` on the unfiltered, design-ordered counts. Then come the unchanged `filterByExpr`, the zero-variance filter (`keep_var`) and the second TMM pass (`norm_factors_filtered`). Each is a named product. | closed |
| KIT093-F04 | should-fix | The kernel table now covers: the preliminary non-robust QL ranking (`prelim_order`), the control sets, raw logCPM, contrast parsing (`trxParseContrastFile`), the annotation join with `DisplaySymbol` and `top_de`, Glimma inputs and the summary text (renderers). Every computational chunk at rmd:115-779 maps to a kernel or renderer. Glimma's HTML write and its explicit success/failure record are **routed** to KIT-095's publication operation. | closed (Glimma write routed to KIT-095) |
| KIT093-F05 | should-fix | The kernel and renderer tables name every product. Parity is one expectation per product. Figures are compared by `layer_data()`, never by objects or file hashes. **Routed to KIT-040** (KIT040-F04): freeze every product named in KIT-093's tables. | closed in the ticket; KIT-040 must freeze the list |
| KIT093-F06 | should-fix | Preflight classes added: `trx_filter_empty`, `trx_grid_too_few_features`, `trx_grid_no_member`, `trx_da_rank_deficient`, `trx_contrast_unparseable` and `trx_contrast_not_estimable` (the zero library is under the count policy). The winner is the first maximum of the unrounded composite score. Failed members are dropped as in the notebook and recorded in the added product `grid_failures`, a declared divergence. Acceptance criteria cover the tie and the forced failure. | closed |
| KIT093-F07 | should-fix | `validateTrxCounts()` refuses NA/NaN, ±Inf, negative, fractional and zero-library values, each with its own class. Integer-valued doubles pass, and counts are stored as double. The check runs in validity and at every kernel entry, before any statistics. `round()` is kept verbatim. `ruvs_counts` must be `identical` to the filtered counts. | closed |
| KIT093-F08 | should-fix | `validateTrxDesign()` returns the design in design-file order with `count_index`. Levels are declared and stored in `args$group_levels`. The model frame column is always `group`. Two acceptance tests: shuffled columns, and a non-`group` `group_id`. **Routed to KIT-040** (KIT040-F02): the corpus uses the notebook's `RPMI`/`Sera`. | closed in the ticket; KIT-040 holds the levels |
| KIT093-F09 | should-fix | `m_source = "technical_replicate"` reads the column named by `technical_replicate_id` (`tech_rep_group` in the corpus), never `replicates`. Three classed refusals are under `trx_ruviii_replicates_invalid`. The oracle is a direct `ruv::RUVIII` call, and the product must differ from the group default. `m_label` is set, and the pseudo-replication claim limit is stated. User decision 4: KIT-095 places RUVIII as a product-only export after DA with `M` in its parameters. **Routed to KIT-040** (KIT040-F07): two technical replicates per biological sample. | closed in the ticket; corpus routed to KIT-040 |
| KIT093-F10 | should-fix | `test-workflow-omics-owner.R` is in `files_to_modify`, and its expected set is the only intended observable change outside the new files. The declaration entry is written out. KIT-093 owns it and KIT-095 only reads it. The DESCRIPTION `Collate` addition is named. **KIT-095's** `files_to_modify` still lists `utils_workflow_class_declarations.R`; that ticket's own repair should narrow it to a read. | closed in KIT-093; KIT-095 note routed |
| KIT093-F11 | should-fix | The tests read the reference from `kit-040/` by repository path and check its sha256 first. The trx test files must show SKIP 0 on the source tree with the KIT library (skips are allowed only under `R CMD check`). One installed-build run of `test-trx-kernels.R` is required. | closed |
| KIT093-F12 | should-fix | The SummarizedExperiment divergence is recorded with its reason (the kit's declaration as the container contract, as for lipid/metab). Validity plus `createTranscriptomicsAssayData()` are required. The limits of about 50 lines and 5–7 parameters are restated. | closed |
| KIT093-F13 | minor | No RNG and no seed parameter are stated, with the reason. The two-run determinism test is an acceptance item. `environment.txt` records BLAS/LAPACK, threads and versions against KIT-040's. | closed |
| KIT093-F14 | minor | Tolerances live in `kit-093/parity-criteria.md`, committed before the first comparison. They are admissible only for a declared divergence; a non-identical result from an unchanged call is a defect. PCA is pinned to `stats::prcomp`, never `pcaWithMissing`. | closed |
| KIT093-F15 | minor | `control_source = "inferred_least_significant"` on the grid, RUVs and RUVIII products, with an acceptance item. RUVIII also carries `m_label`. | closed |
| KIT093-F16 | minor | A kernel interface table (kernel, file, notebook lines, products) and a renderer table were added. There is a 900 s ceiling per test file and one attempt plus one correction pass. The estimate is 12 h. `needs_planning` stays false because the interface is now in the ticket. | closed |

## Routed-in findings (from the other audits)

| ID | Repair in KIT-093 rev 1.1 | State |
|---|---|---|
| KIT096-F07 | The declared level order is stored in `args$group_levels`, and kernels build the factor from it. More than two groups is refused (`trx_design_groups_unsupported`), because the notebook's `coef = 2` ranking defines only a two-group test. KIT-096's design tab writes the order. | closed; extending to three or more groups is a user decision (below) |
| KIT096-F08 / KIT099-F07 / KIT095-F03 | The filter section defines: the pre-pass as TMM-effective library size on the unfiltered counts; the minimum group size and the median library size, recorded as `filterByExpr.default` computes them, for KIT-095's partitioned variant only; the zero-variance filter; and the second TMM pass. The whole kernel's mask always comes from the unchanged edgeR call. | closed (KIT-099's reporting of the pre-pass phase stays with KIT-099) |
| KIT095-F02 | KIT-093 owns the publication renderers in `R/func_trx_publication.R`: the CCA panels and workbook data; the norm-stage PCA/RLE grid; the volcano and the significance bar plot; the DE table and the top genes; Glimma inputs; the RUVIII before/after grid; the summary text. They return figure objects with `layer_data()` and data frames, and write no files. KIT-095's publication operations write them. The notebook's functions are moved, not copied from package helpers, and shared helpers are not modified (decision 3; display reuse by parameter is KIT-097's). | closed |
| KIT097-F03 | `top_table` keeps edgeR's BH `FDR` under that name. No Storey q-value is computed or substituted. Renderers label it "FDR (BH)". An acceptance item checks that there is no q-value column. | closed in KIT-093; the tab labels stay with KIT-097 |
| KIT096-F05 | `validateTrxCounts()` is callable on its own, before any statistics, with one class per fault under `trx_counts_invalid`. It is handed to KIT-096's import tab. | closed in KIT-093; display stays with KIT-096 |

## Counts

| Group | Total | Closed in KIT-093 | Closed, part routed | Routed elsewhere |
|---|---|---|---|---|
| KIT-093's own findings (F01–F16) | 16 | 10 | 6 | — |
| Routed-in findings | 7 | 7 | — | — |

For the six own findings closed with a part routed:
- KIT-040 freezes the products and provides the levels, the technical-replicate structure and the install (F01, F05, F08, F09).
- KIT-095 writes Glimma and narrows its file list (F04, F10).

Of the routed-in findings, KIT096-F08, KIT099-F07 and KIT095-F03 are one repair.

## Not closed here, or needing a decision

1. **More than two groups.** KIT-093 refuses them, and the notebook defines no ranking for them. Allowing them (for example by an F-test over all group coefficients in the preliminary fit) changes the science, so it needs a user decision, and then a kernel and reference addendum.
2. **KIT-040 must carry the routed items** (install, frozen product list, `RPMI`/`Sera`, `tech_rep_group`). Until it does, KIT-093's parity acceptance cannot run.
3. **KIT-095 ownership.** KIT-095's `files_to_modify` should treat `R/utils_workflow_class_declarations.R` as read-only for the TRX entry. That is outside this file's write scope.
