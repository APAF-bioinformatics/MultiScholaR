# KIT-221 outcome

Base: `origin/v0.9-staging` 176f9d3ff. Tests-only commit 41c69770f (nothing under `R/`); fix commit 4654fd600.
The sections below the R1 table were added after the coordinator's rulings of 2026-11-26; R1 stands unchanged.

## R1 (run 45713440865, request kit221-r1-41c69770, one test job and one script job)

Observed (tables in `r1/`): `test-kit221-da-several-comparisons.R` 43 passed, 38 failed, 0 errors, 0 skipped (predicted
42 / 39); `test-kit218-da-wrapper-tail.R` (G4) 80 passed, 0 failed, 2 skipped (G1 and M2 skip because the fix is
present); `g1.R` exit 1 after 2 s.

| Test | Predicted passed / failed | Observed | Assertion that failed |
|---|---|---|---|
| W1 | 2 / 8 | as predicted | first: `sort(unique(rows$comparison))` is `B_vs_A` only; then 60 not 120 rows (twice), `num_sig_da_molecules` and the first-go table `B_vs_A` only, 1 volcano row not 2, no `log2FC:C_vs_A` column, long table `B_vs_A` only. No error |
| W1b | 0 / 1 | as predicted | `expect_null(run$error)`; M1 = "number of items to replace is not a multiple of replacement length" (the ticket's reading) |
| W2 | 8 / 2 | as predicted | comparisons set; C_vs_A rows (0 against 60); B_vs_A identical |
| W3 | 8 / 5 | as predicted | `length(contrasts_results$results)` 1 not 2; comparisons; 60 not 120 rows; set; C_vs_A rows |
| W4 | 10 / 10 | as predicted | as W1 plus set and C_vs_A rows |
| W5 | 9 / 5 | 10 passed / 4 failed (the `length` assertion passes; see below) | comparisons `B_vs_A` only; `nrow` n not 2n; set; C_vs_A rows |
| W7 | 5 / 5 | as predicted | set; B_vs_A rows; `log2FC`, `raw_pvalue`, `fdr_qvalue` of the rows labelled C_vs_A differ from the one-comparison call; mean `log2FC` of P11..P16 under `C_vs_A` = 1.0322 (B - A's effect: the factor premise holds; stop rule 6 does not apply) |
| W8 | 0 / 1 | as predicted | `expect_null(run$error)` |
| W9 | 0 / 2 | as predicted | both calls stop with "incorrect number of dimensions", not of class `prot_da_contrasts_invalid` |

**Missed prediction (W5, one expectation).** I predicted that `length(multi$value$contrasts_results$results)` is 1
on base. It is 2: `differentialAbundanceAnalysis()` keeps every comparison in `contrasts_results` (ticket 2.1: only
the derived tables are first-only), which the ticket states and I mis-applied to W5 (W3's assertion is right: the
sibling's `contrasts_results` holds one). The row therefore passes on base; the other four W5 expectations fail as
predicted. Per the ticket (a missed prediction: stop and report) no further run was sent and no fix was committed.

**Instrument defect (g1.R).** `ProteinQuantitativeData` is not exported, so `g1.R` (run through the installed package)
stops at the first call: "could not find function". No digest was produced, no value judged. The mechanics fix is
`ns$ProteinQuantitativeData(...)` in `make_object`; not sent (the stop above), so none of the three instrument reruns is used.

## Not checked

Everything after R1: the fix's behaviour, G1 and G3 digests, R2 and R3, the check comparison, NEWS, ledger row.

## Records added 2026-11-26 (coordinator's rulings)

**Dated correction of R1's prediction.** Predicted 42 passed / 39 failed for the new file, observed 43 / 38: one
expectation, W5's `length(contrasts_results$results) == 2`, which passes on base because `differentialAbundanceAnalysis()`
keeps every comparison in `contrasts_results` (ticket 2.1). Cause: a mis-read of 2.1. No expectation was changed after
the result (the test now carries a comment marking this assertion a guard that passes on base and on head); the
ticket's correction pass is not used. R1's result stands and is not re-sent.

**The fix** (commit 4654fd600) was drafted before R1's verdict and committed after it; it was not changed in the
light of the result beyond what ticket 3.1 specifies (R1 showed no behaviour the draft had to answer).

**Instrument rerun 1 of 3** (g1.R, the namespace accessor for `ProteinQuantitativeData`): run 45715339444, script job
alone on e80e558d (R/ = base), predictions in `r1-predictions.md`; about 1 core-hour. Head digests: one script job
attached to R2.

**R2 scope (45 files, the ticket priced about 40; accepted by the coordinator).** The ticket's grep at 176f9d3ff,
helpers dropped, gives 43 files (the new file included); plus `test-gui-entry-proteomics.R` and
`test-workflow-operations.R` (not matched by the grep): 45.
- test-dia-limpa-ruv-rules.R
- test-enrichment-significance.R
- test-gui-entry-proteomics.R
- test-kit128-no-biological-replication-prot.R
- test-kit129-run-path.R
- test-kit174-ruv-grouping.R
- test-kit200-prot-da-no-value-proteins.R
- test-kit201-treat-moderation.R
- test-kit204-prot-da-not-tested.R
- test-kit205-enrichment-background.R
- test-kit208-gprofiler-identifiers.R
- test-kit212-da-current.R
- test-kit215-not-tested-lists.R
- test-kit218-da-wrapper-tail.R
- test-kit220-da-threshold-consumers.R
- test-kit221-da-several-comparisons.R
- test-kit243-summary-current-outputs.R
- test-kit248-correlation-skip-record.R
- test-kit249-prot-da-qvalue-fallback.R
- test-lipid-07-s4-da-plot-methods.R
- test-lipid-08-s4-da-analysis-methods.R
- test-metab-02k-da-numsig-barplot-characterization.R
- test-metab-02p-da-analysis-characterization.R
- test-metab-s4-da-methods-direct-shared.R
- test-module-ci-prot-da.R
- test-omic-da-module.R
- test-omics-da-list-routing.R
- test-prot-07b-da-handlers-compat.R
- test-prot-07c-da-results-characterization.R
- test-prot-07c-da-results-output-shared.R
- test-prot-07-da-analysis.R
- test-prot-07d-da-model-stats-shared.R
- test-prot-annotation-uniprot-empty.R
- test-prot-da-compat-shared.R
- test-prot-design-log2.R
- test-protein-da-golden-master.R
- test-prot-ruv-factors-in-model.R
- test-prot-s4-norm-methods-direct-shared.R
- test-reference-mapping-enrichment.R
- test-workflow-analysis-record.R
- test-workflow-figure-paths.R
- test-workflow-memory-mechanism-member.R
- test-workflow-operations.R
- test-workflow-saved-da.R
- test-workflow-saved-installed.R

**Correction of 2026-11-26 (implementation audit I01):** the ticket's grep at 176f9d3ff gives 43 test files without the
new file (44 with it), not 43 with it: dropping 'helpers' by a name match also dropped
`test-metab-02q-da-analysis-helper-characterization.R`, which calls the metabolomics `differentialAbundanceAnalysisHelper`.
The declared scope was 46 files; R2 ran 45. The omitted file tests `R/func_metab_s4_da_methods.R`, which this change does
not touch (by reading unaffected; 19 passed on staging, merge train h). It runs in the merge train.

R2 price and predictions: `r2-predictions.md` (expected 36.0, ceiling 55.2 core-hours with the head script job).

**Stored values the change could move (declared by file and line; searched `tests/` whole, `tests/testdata/`
included, `data/benchmark/` excluded by path).** The change alters outputs only for a call with more than one
comparison, and the message of a non-data-frame input. Searches: every call of `differentialAbundanceAnalysis(`,
`protein_deAnalysisWrapperFunction(` and `deAnalysisWrapperFunction(` under `tests/`; the pinned texts
'incorrect number of dimensions', 'number of items to replace', 'subscript out of bounds'; the names of the returned
fields (`da_proteins_long`, `significant_rows`, `contrasts_results_table`, `num_sig_da_molecules`) in non-R files of
`tests/`. Result: none moved. Every proteomics call passes one comparison (e.g. `test-kit129-run-path.R:313, :570,
:619, :706, :898`, `test-prot-07-da-analysis.R:139`, `test-prot-07d-da-model-stats-shared.R:698-745`,
`test-workflow-saved-da.R:54`); the calls with two comparisons (`test-kit200-prot-da-no-value-proteins.R:743`,
`test-kit249-prot-da-qvalue-fallback.R:342`) use the list method on lipidomics and metabolomics assays, not changed;
the three sibling tests (`test-prot-da-compat-shared.R`, `test-protein-da-golden-master.R:16`,
`tests/test_limpa_connection.R:20`) pass a one-row table or NULL. The pinned texts occur only in
`test-general-enrichment-shared.R:175` (unrelated). The field names occur in `tests/testdata/gui-entry/prot_*.json`
(application captures, one comparison per kernel call) and in `tests/testdata/kit227/base/*.rmd` and
`tests/testdata/kit236/base/*.rmd` (report text); none holds a value of a multi-comparison call. No stored value,
digest or pinned text is predicted to move; none is edited.

**Register notes for the coordinator (ticket text minors, not applied by me).** R02: in 2.3 replace "No reported
number is wrong: each comparison's moderated statistics" by "For comparisons given as text, no reported number is
wrong: each comparison's moderated statistics". R03: front matter `description` "By reading at v0.9-staging e6791275:"
-> "By reading at v0.9-staging 0cf0e440 (re-read at 9a0bf9f0 and d048898d; no change under R/, inst/, tests/):";
2's heading "(SHOWED at e6791275 unless marked)" -> "(SHOWED at 0cf0e440, unchanged at d048898d, unless marked)";
8 Files "`R/func_prot_da_compat.R:115`" -> "`R/func_prot_da_compat.R:114-115` (3.1 item 2)". R04: the tests-only
commit carried (A)'s rows (OD-1 was answered (A) before it), so its text is not needed. Step 0 also notes that the
ticket's sibling check is placed right after `contrasts_tbl` is read (before line 114, as 3.1 item 2 requires).

**A4 (full suite): by the checkpoint after merge (coordinator's ruling, 2026-11-26).** R3 is not run by this ticket.

**Instrument rerun 1 result.** Run 45715339444: script exit 0 after 5 s; `g1-base.tsv` (93 digests: limma 19, limpa 20,
ruv 17, kernel 19, alias 18) is kept in this directory. Prediction (exit 0, five cases, at least ten rows each) met.

## R2 (run 45716909041, base e80e558d, head 53ac0f87; R/ at head = 4654fd600; observed)

Request kit221-r2-53ac0f87: 45 files, 3 shards per side, check on both sides, head script. Tables in `r2/`.

| Item | Predicted | Observed |
|---|---|---|
| new file, head | 124 passed, 0 failed (W5's guard included), verdict `fixed` | 124 passed, 0 failed, 0 errors, 0 skipped; `fixed` (base 43/38, as R1) |
| the other 44 files | `same` | `same` (44 of 44, per-file counts and masked failure messages) |
| failures present on both sides | unchanged | `test-dia-limpa-ruv-rules.R` 1, `test-gui-entry-proteomics.R` 101, `test-kit129-run-path.R` 6, `test-kit220-da-threshold-consumers.R` 2, `test-omic-da-module.R` 81, `test-workflow-saved-da.R` 196 on both sides (stored references recorded on the laptop; the ticket cites 171 for the last, the run gives 196 on both sides, verdict `same`) |
| `test-kit218-da-wrapper-tail.R` (G4) | 80 passed / 2 skipped both | 80 / 2 both |
| R CMD check | identical, nothing only on head | base and head: 0 ERROR, 3 WARNING, 3 NOTE; 6 findings identical; only on head: none; only on base: none |
| G1 and G3 | every digest equal, base against head | `g1-base.tsv` (run 45715339444, 93 rows) and `g1-head.tsv` equal in every case, field and digest (cut -f1-3 diff empty); head row count 93 |

No missed prediction on R2.

Head observations (R2 artifact, `logs/test-kit221-da-several-comparisons.log`): W1 with two comparisons gives
`names(da_proteins_long)` = Protein.Ids, fdr_qvalue, raw_pvalue, log2FC, comparison, numerator, denominator,
log2norm.1.left to log2norm.4.left, log2norm.1.right to log2norm.4.right, raw.1.left to raw.4.left, raw.1.right to
raw.4.right, gene_names (the per-sample columns numbered, as NEWS says); W7 on head: mean log2FC of P11..P16 under
B_vs_A = 1.0322, under C_vs_A = 0.0764 (on base 1.0322 sat under the label C_vs_A). Of the 93 digests, 20 are of the
four mocked pre-DA plots (empty ggplots, the same digest in every case); 73 are of DA outputs.

**Citation slips (implementation audit I08, 2026-11-26; recorded, nothing else edited).** (1) The header of
`tests/testthat/test-kit221-da-several-comparisons.R` cites `kit218MockPlots :86-91`; the helper is at
`tests/testthat/test-kit218-da-wrapper-tail.R:79-83`. It is a test comment and is left as it is (no test edit). (2) The R1
price in the core-hours table below (2.0 expected, 6.67 ceiling) is the implementer's; the ticket's R1 row is 1.0 expected,
4.0 ceiling. (3) `r2-predictions.md` repeats the stale "171" failing expectations of `test-workflow-saved-da.R`; the run
gives 196 on both sides (verdict `same` held, as the R2 table above says).

**NEWS refusal sentence (implementation audit I04, 2026-11-26): a DRAFT for the wording pass (VAL-41).** The clause
"(before: 'incorrect number of dimensions')" was true only for a vector or a list. It now reads, in `NEWS.md` on the merge
train: "(before: an error, 'incorrect number of dimensions', for a vector or a list; a one-column character matrix was
read, by every row in `differentialAbundanceAnalysis()` and by its first value in `protein_deAnalysisWrapperFunction()`)".
A one-column character matrix, accepted before, is now refused, as VAL-193 decided. The owner sets the final wording.

## Core-hours (job minutes x 4 / 60, from the jobs' start and end times)

| Run | Job minutes | Core-hours | Expected | Ceiling |
|---|---|---|---|---|
| R1, 45713440865 | 17.5 | 1.17 | 2.0 | 6.67 |
| g1 base, 45715339444 (instrument rerun 1 of 3) | 10.4 | 0.69 | 1.0 | 2.67 |
| R2, 45716909041 | 351.9 | 23.46 | 36.0 | 55.2 |
| Total | 379.8 | 25.3 | | ticket total 473.33 (R3 excluded by the ruling) |

## Acceptance

| # | Criterion | Result |
|---|---|---|
| A0 | OD-1 and OD-4 answered | met: OD-1 (A) by VAL-193 (1); OD-4 to KIT-230 (parked); no change to 3, 3.1 or 7 |
| A1 | tests-only commit: failing-first rows fail at the stated values; G1, G3, G4 pass; M1 recorded; head: all W rows pass | met with the recorded correction (W5's `length` assertion is a guard); M1 = "number of items to replace is not a multiple of replacement length"; head 124/0 |
| A2 | G1 and G3 digests equal base against head | met (93 digests; the base side from the script job on the tests-only side, same R/ as base) |
| A3 | R2: every file `same` except the new file; nothing only on head in check | met |
| A4 | full suite base against head | by the checkpoint after merge (coordinator's ruling, 2026-11-26) |
| A5 | NEWS, ledger row; `NAMESPACE`, `man/`, `DESCRIPTION` unchanged; "For the change report" | NEWS entry and ledger row `VAL-357` added; `git diff` of NAMESPACE, man/, DESCRIPTION is empty; section below |

## For the change report

**What was wrong and who was affected.** Only people who script the package and pass several comparisons in one call.
`differentialAbundanceAnalysis()` fitted every comparison, then built its significance table, wide and long tables,
both count tables and the volcano and p-value plots from the first comparison only (only `contrasts_results` kept
all); `protein_deAnalysisWrapperFunction()` (and its alias `deAnalysisWrapperFunction()`) given a base data frame of
several rows analysed only its first row. No message in either case. For comparisons given as text no reported
number was wrong; comparisons were missing. With a comparisons column stored as a factor whose levels were not in row
order, the S4 method labelled one comparison's statistics with another's name (R1, W7: the effect proteins' mean
log2FC under the label C_vs_A was 1.03). Not affected: the application, replay and the managed worker (one
comparison per kernel call, `R/utils_workflow_dia_da.R:535-543`); the shipped workbooks' DA call (archived by KIT-234).

**Code references.** Base 176f9d3ff: `R/func_prot_da_model_methods.R:295` (`contrasts_tbl[, 1]`), `:430`, `:433`
(`results[[1]]`, `names(...)[1]`), `:477-480` and `:509-514` (one-table nested lists); `R/func_prot_da_compat.R:115`
(`contrasts_tbl[, 1][[1]]`). Fix 4654fd600: `R/func_prot_da_model_methods.R` check and `as.character(contrasts_tbl[[1]])`
at `:295-302`, `contrast_tables` (`purrr::imap`) at `:443-461`, `list(contrast_tables)` at the two calls; 
`R/func_prot_da_compat.R:32-38` (check) and `:121` (`as.character(contrasts_tbl[[1]])`).

**In earlier releases.** The same pattern is in every release read: `main` 0.4.1.2 (`main:R/func_prot_da.R:4093`
`contrasts_tbl[, 1]`, `:4119` `results[[1]]`, `:4168-4171`, `:4200-4203`, the active copy; the same in
`main:R/func_prot_s4_objects.R:780`, `:806`; the sibling `main:R/da_analysis_function_wrapper.R:112`), the published
v0.3.6.3 (`R/func_prot_de.R:3278`, `:3303`), v0.5.0-rc1 (`R/func_prot_da_model_methods.R:356`, `:382`;
`R/func_prot_da_compat.R:112`) and v0.9.0 (`:295`, `:427`; compat `:115`). In each, the application's DA step called the
method once per comparison with a one-row text data frame (`main:R/mod_prot_da.R:1153-1168`, `:1316`) and the DIA
workbooks called the sibling once per comparison with a one-row tibble slice
(`main:inst/workbooks/proteomics/DIA_workflow_experienced.rmd:1107-1109`): users of the application and the workbooks
were not affected; scripted users of every release were, in the same way (by reading), the factor case included
(`main:R/func_prot_da.R:2595`, `:2822`). The archived workbooks' fixed-name summary files, which keep only the last
comparison, are KIT-230's, parked with the workbooks (VAL-193 (2), VAL-209). A multi-comparison result written with
`outputDaAnalysisResults()` now fills its 'all' files with every comparison (before: the first).

**What changed in results.** A scripted call with several comparisons now returns every comparison in every table and
plot; each comparison's `log2FC`, `raw_pvalue`, `fdr_qvalue` equal those of a call with it alone (identical on the
limma, limpa-with-a-fake-`dpcDE` and RUV-factors routes, W2, W4, W5). With several comparisons the long table's
per-sample columns are numbered by position. A tibble now works in the S4 method (it stopped before). A
non-data-frame stops with `prot_da_contrasts_invalid` (before: "incorrect number of dimensions"). Comparisons are read
as text, closing the factor mislabelling. One-comparison outputs are unchanged (G1, G3).

**Evidence and limits.** R1 45713440865, g1 base 45715339444, R2 45716909041. One synthetic design (60 proteins, three
groups of four; the RUV route 12 runs relabelled A, B, C with one factor); the limpa route uses a fake `dpcDE`. Not
tested: digit-named groups with several comparisons (the group-name mapping is applied per table), real data, the full
suite (A4, after merge). `outputDaResultsAllContrasts()` still expects one comparison per list element (ticket F3).

## Not checked

The full suite; digit-named groups with several comparisons; the limpa route with the real `dpcDE`; designs with
covariates beyond the RUV route (identity to rounding is stated in NEWS, tested only bitwise on the three routes
here); real data; Windows.
