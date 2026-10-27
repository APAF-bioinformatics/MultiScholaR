# KIT-064: static audit repairs (roar 2026-10-27), exactness-gated

Status 2026-10-28: **completed** (user decision 2026-10-28, "keep the repair": the `drop = FALSE` in `plotPcaHelper` stays; the LFQ and TMT fixture PCA product changes from a degenerate NaN plot to a PCA, and KIT-065 re-bases the fixture records on its build). Install `kit_064` = `core/browser-smoke.tdDG9e`. Code 1ab2318b (item A 3880ddb0; items B1-B6 4ddb9a3b, a5c715a7, f5905a28, 1ed229c2, cdbe9f22, 8b6c03d9; item C 83f4ec6e + 7765d9b7; the BUG-18 revert 1ab2318b). Base build `kit_045` = `core/browser-smoke.JGgjxW`. No counted member. The finding register is [AUDIT-ALIGNMENT-20261028.md](AUDIT-ALIGNMENT-20261028.md); this outcome closes its "Action" column.

## Closure list against the register

| Group | Register action | Done | Where | Proof |
|---|---|---|---|---|
| BUG-02 | item A | yes | `saveListOfFunctionalEnrichmentHeatmaps`: `base_path = results_dir`; its guard `&&` | the enrichment shared test's `savePlot` seam now takes `base_path` and asserts it |
| BUG-03 | item A | yes | `initializeProtDesignImportedWorkflowState`: the type is the config's, else detected (`DIA` or `TMT`); the `NULL` comparison is gone | `test-prot-04a`: a protein-level import without `workflow_type` resolves to `TMT` (aborted before) |
| BUG-07 | item A | yes | the STRING-DB table filters on `input$q_cutoff %||% 0.05` | `test-prot-11`: explicit cutoff and the default both filter |
| BUG-08 | item A | yes | the five render helpers read their arguments through `protEnrichRenderValue()` inside the render expression; the four registration sites pass accessors (`function() enrichmentData$...`, `function() input$...`, `rawContrastNameFn` itself) | `test-prot-11`: a re-evaluating fake renderer sees the placeholder first and the populated table after the source changes (a value captured at registration could not); GUI reproduction is KIT-065's |
| BUG-15 (652) | item A | yes | `plotVolcano` facets only on a non-NULL, non-NA, non-empty formula | `test-general-volcano-glimma`: `""` and `NA` give an unfaceted ggplot (aborted before) |
| BUG-14 | item B1 | plotting and DA helpers | `plotPcaHelper`/`plotPcaListHelper` (4 subsets), `runTests` (2), `runTestsContrasts` (1) with `drop = FALSE`; the two kernel sites (`func_prot_s4_missingness.R`, `func_prot_limpa_qc_helpers.R`) **not changed** | replay |
| BUG-16 | item B1 | yes | `seq_len` in `runTests`, the cnetplot helper (8 sites, `seq.int(n + 1, length.out = )` for the tail), the composite generator | replay; the enrichment and norm tests |
| BUG-04/05 | item B1 | yes | `union()` for the dormant pair; scalar guards for `weights` (3) and `sample_rows_list` | `test-prot-07d`; the replay's DA tables (every caller passes `NA`) |
| BUG-18/19 | item B1 | 19 yes; **18 reverted** | the `prcomp` fallback drops constant columns. The threshold `max(2L, floor(0.05 * n))` was applied, then reverted at 1ab2318b: the LFQ and TMT fixture records differed from KIT-048's in exactly one product, the normalisation step's `dia_normalization_pca` (the PCA scores are recorded science; the fixtures have 5 samples); the other B1 edits stay | the records before and after the revert (below) |
| BUG-23 | item B2 | proteomics and shared | `vapply` in `func_prot_da_volcano_glimma.R` and the STRING-DB reporting helper; lipid/metab sites carried | tests |
| BUG-15 (class) | item B2 | non-kernel sites | `&&`/`||` in the revigo, enrichment plotting, pathway and volcano helpers (7 sites); the phospho sites **not changed** (their operand may be a vector; a `&&` would abort on it) | tests |
| BUG-13 | item B2 | yes | `dplyr::all_of()` at the three phospho sites and in `removeEmptyRows` | replay (`removeEmptyRows` sits under the kit's filters) |
| BUG-22 (devices) | item B3 | yes | twelve device sites close on error (`tryCatch(..., finally = dev.off())`; `on.exit` in `saveListOfPdfs`); `dev.size(units =)` | tests; no output change |
| BUG-17 (enrichment) | item B4 | yes | the barplot's `-log10` limit clamped at `.Machine$double.xmin` | tests |
| BUG-11 (revigo) | item B | **not a defect** | `plotPathwayCluster` maps literal shapes with `scale_shape_identity()` | read |
| BUG-24 | item B5 | enrichment and annotation | `requireOptionalPackage()` at the entry of 15 functions (clusterProfiler, GSEABase, tm, UniProt.ws, AnnotationDbi, GO.db, mixOmics, KEGGREST, tidytext); the kit backends and `limpa` stay behind their probes | tests (the guards pass where the packages are installed) |
| BUG-12 | item B6 | yes | the shadowed `plotVolcano` generic and the two dead unexported generics removed; `plotVolcano.Rd` documents the plain function | `NAMESPACE` unchanged (checked); R CMD check |
| BUG-25 | item B | partial | the three grouped returns: the phosphosite helper already ends in `ungroup`; the lipid/metab ones carried; the stringr empty pattern is not at the reported line (not located; recorded) | read |
| BUG-20 | item B | **not applied** | the percentage messages print `NaN`/`Inf` without aborting; the kernel scores untouched (recorded as a cost/benefit decision, not a repair) | read |
| BUG-01, 22 (pins), 24 (backends) | item C | yes | the fixture's fork mode refuses on a non-forking R with a message; the `RAYON_NUM_THREADS` pin, `set.seed(42)` (3), `options(multischolar.test_mode)` and the two `testthat` modes carry a comment naming the decision | read |
| BUG-06 | item B3 | yes | `units =` | read |
| BUG-09 | item D | **not done** (the user's word) | | |
| BUG-10, 17 (lipid), 21 (`.GlobalEnv`), 16/23 (lipid/metab) | stage 3 | carried to KIT-030/031 | | |
| BUG-21 (`<<-`), 11 (PCA), 12 (`createGridQC*`) | not defects | | | read |

## Exactness gates

- KIT-028 replay (in process, the screen origin, 16 states) after item A on a clean worktree of 3880ddb0: **exact, 0 mismatches, 142 s**. After items B and C on 83f4ec6e: **exact, 0 mismatches, 173 s**. After the BUG-18 revert on 1ab2318b: **exact, 0 mismatches, 174 s**. (The two remaining item C commits are comments only.)
- The three screen records re-recorded on the first KIT-064 build (`core/browser-smoke.r4eG1H`, 7765d9b7) and compared state by state, field by field (class, rows, components, support, summary, products) with KIT-048's records on `kit_048`: DIA_LIMPA exact on all 13 states; LFQ and TMT identical on every state except the `correlation` state's `dia_normalization_pca` product (its content hash; the Pearson pairs and selection products identical), traced to the BUG-18 threshold and reverted. Re-recorded on the final build (`core/browser-smoke.tdDG9e`, 1ab2318b): DIA_LIMPA not re-recorded (unchanged code path); LFQ and TMT **still differ in that one product**, so the threshold was not the cause. The product was hydrated from both stores and diffed leaf by leaf: KIT-048's PCA product on the LFQ fixture is degenerate (`data` 1 × 11, x label `PC1 (NaN% of top 10% CV)`, limits ±1), KIT-064's is a real PCA (`data` 6 × 11, `PC1 (100% of top 10% CV)`, fitted limits). The cause is item B1's `drop = FALSE` in `plotPcaHelper`: on the fixtures exactly one feature passes the CV filter, the former `data_abundant[which(cvs >= cv_threshold), ]` collapsed it to a vector and the PCA fell into its degenerate branch; that is roar BUG-14's predicted failure, observed. On the screen and 300-sample corpora many features pass, nothing collapses, and the DIA replay stayed exact. **User decision (2026-10-28): keep the repair.** The records made on the KIT-064 builds are filed beside this outcome (`oracles-lfq_fixture.rds`, `oracles-tmt_fixture.rds` on tdDG9e; `oracles-limpa_screen.rds` on r4eG1H, the limpa path unchanged by the later revert); they are evidence here, and KIT-065 records its own on its build. The replay JSONs of the three gates are filed as `replay-check-kit064-*.json`.
- Targeted tests on the tree after item A (five files): `test-prot-11` 1258/0 (+5 expectations), `test-module-ci-prot-enrich` 92/1, `test-prot-04a` 194/4 (+2), `test-general-enrichment-shared` 104/2 (+1), `test-general-volcano-glimma` 9/0 (+4), all pass/skip at or above the checkpoint baseline. After items B and C, the 33 test files that name a changed function (the dependents grep): 2,300 pass, 0 fail, 0 errors, 35 skips.
- `R CMD check --no-tests --no-manual` at 7765d9b7: 3 WARNINGs, 3 NOTEs, the baseline; at 1ab2318b: 3 WARNINGs, 3 NOTEs, the baseline. Install `kit_064` = `core/browser-smoke.tdDG9e` (1ab2318b); the first install r4eG1H (7765d9b7) is the one the product difference was found on.
- `NAMESPACE`: no diff against ee7806a9. No kernel body changed except the two allowed sites gated by the replay (`removeEmptyRows`'s selection, the DA helpers' `drop = FALSE`).

## Deviations recorded

- Item C landed in two commits: the first stopped at the kegg helper's duplicated `set.seed` line (an edit-ordering slip), the second completed the notes.
- The replay was run after item A and after items B and C together, not after each of the six B commits; each B commit is a separate revert unit all the same.
- BUG-20 (zero-denominator guards at message sites) and the phospho `&` sites were not changed, with the reasons in the table.
