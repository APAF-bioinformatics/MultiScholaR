# KIT-204 outcome: proteomics DA rows without a p-value are counted as "Not tested"

Ticket `tickets/kit/KIT-204.md` revision 1.2 (step 0 applied re-check 1's R01 and R02 and re-read every citation at
`v0.9-staging` dd769671). Branch `kit-204` from dd769671 (KIT-200 merged, VAL-211). No R was run on the laptop; every
run below is a GitHub-hosted `run-tests` run (R 4.6.1, KIT-134 snapshot library, reference BLAS, ubuntu-22.04).

## Commits

| Commit | Content |
|---|---|
| 9813e4f7 | tests only: `tests/testthat/test-kit204-prot-da-not-tested.R` (N1 to N9) |
| d77ac7b7 | tests only: two instrument repairs after R1 (below) |
| b5bf4b5c | the fix (7.1), `man/` pages, NEWS, ledger row `VAL-323` |
| 29848a2f | the one correction after R2: the three barplot filters leave out "Not significant" and "Not tested" |
| (this commit) | records |
| 89ce2a1f | merge with staging 7327607f (KIT-220 conflicts as the audit register's section 6) |
| d1bdcfe3 | merge with staging abfc99ac (ledger only resolved by hand; see "Second merge") |
| d792847f | R2d declared before sending |

## Runs (`runs.tsv`)

R1 45574481570, R1b 45576574709, R2 45578443378, R2c 45585667100, R3 45592061384. Core-hours 115.4 in all
(0.75 + 0.69 + 14.50 + 14.55 + 84.89; job minutes x 4 / 60, set-up included) against an expected 85.8 and a ceiling of
262.7. Above the expected figure because the correction's R2 was spent (14.55) and R2/R3 took longer than KIT-186's
and KIT-192's measured figures.

## What was observed

### R1 (instrument defect, no result judged on the two items)

Observed: N1, N2b, N3 (i), N5, N6 and N7 failed at their stated values. Two failures were in the test script's mechanics:

- N2's B-A expected value was NA. The test computed "not significant" as `!(fdr_qvalue < 0.05)` over tested rows, which
  is NA on a row with a finite p-value and a missing q-value; the package's `case_when` counts such a row as
  "Not significant" (its last branch). B-A on this fixture has such rows (see "Recorded, not changed" 1). Repair: the
  expected value counts a missing q-value as not significant, as the package does. No expected value or threshold
  changed.
- N9 errored before any assertion (`dynGet("da_q_val_thresh")` not found): the test did not record on the object the
  arguments `produceDIADaPublication` records (`R/utils_workflow_dia_da.R:586-590`). Repair: the test records them as
  that function does. Because B-A on this fixture has no significant protein (same cause), N9's manifest also gets the
  contrast C-B, which carries the group-B effect, so that N3 (iii)'s precondition (a NumSig figure is drawn) can hold.

Instrument rerun 1 of 3 (VAL-54 conditions i to iv: no expected value, threshold, declared comparison or stop rule
changed; recorded here and in `runs.tsv`; ceiling unchanged).

### R1b (tests-only d77ac7b7): failing first, as stated

`r1b-failures.txt`, `r1b-tests-only.log`. 50 passed, 31 failed, 0 errors.

| Test | Failed at |
|---|---|
| N1 | Not significant 3 (expected 1); Not tested absent (expected 2); three statuses; counts `c(1L, 2L, 1L)` |
| N2 | C-A, both count tables: Not significant 60 (expected 59), Not tested absent (expected 1). B-A passes (no untested row), as the ticket states |
| N2b | B-A and C-A, both tables: Not significant 61 (expected 60 and 59), Not tested absent (expected 1 and 2); summary line absent for B-A and C-A |
| N3 (i) | no "Not tested" bar in C-B's every-status barplot (preconditions held: C-B has a significant protein, P5 untested) |
| N5 | line absent: manifest C-A (n 1), manifest B-A (n 0), list C-A (n 1). The list without `raw_pvalue` passes (no line) |
| N6 | line absent (the computed five lines otherwise equal) |
| N7 | `total` = -1 (the seeded v1 product served); `not_tested` absent |
| N9 | no "Not tested" row for C-A in `da_proteins_num_sig_da_molecules.tab` |

Guards passing on the tests-only commit: N3 (ii), N3 (iii), N4 (a), N4 (b), N8; N9's report-sentence guard. Stop rule 1
(N1 or N2 passes on base) did not apply.

### R2 (fix b5bf4b5c against base dd769671): a result caused by the change

`r2-compare.tsv`, `r2-compare-failures.txt`, `r2-check-compare.txt`. Scope (section 9; the grep at the branch point,
helpers dropped, gave 24 files; added: the GUI-entry proteomics file holding the five lanes, `test-lipid-08-...`,
`test-metab-02p-...`, `test-metab-02q-...` (the metabolomics S4 DA counterparts), KIT-200's file and the new file; 30
in all): `test-gui-entry-proteomics.R test-omic-da-module.R test-omic-enrich-module.R test-workflow-saved-da.R
test-general-filemgmt-report-helpers-shared.R test-lipid-02ab-da-module-characterization.R
test-lipid-02b-da-module-contracts.R test-lipid-07-s4-da-plot-methods.R
test-metab-02aa-da-display-filter-characterization.R test-metab-02k-da-numsig-barplot-characterization.R
test-module-ci-prot-da.R test-omic-summary-module.R test-prot-07-da-analysis.R test-prot-07b-da-handlers-compat.R
test-prot-07c-da-results-characterization.R test-prot-07c-da-results-output-shared.R
test-prot-07d-da-model-stats-shared.R test-prot-da-compat-shared.R test-prot-da-render-direct-shared.R
test-prot-summary-report-templates.R test-workflow-dia-products.R test-workflow-figure-paths.R
test-workflow-memory-mechanism-member.R test-workflow-omics-application-reads.R test-workflow-operations.R
test-lipid-08-s4-da-analysis-methods.R test-metab-02p-da-analysis-characterization.R
test-metab-02q-da-analysis-helper-characterization.R test-kit200-prot-da-no-value-proteins.R
test-kit204-prot-da-not-tested.R`.

Observed: the 13 `prot_main` summary values moved, each by the one line (as declared); the five GUI-entry lanes,
`test-omic-enrich-module.R` and `test-workflow-saved-da.R` same; check identical (0 ERROR, 3 WARNING, 3 NOTE on both;
none only on head). Not declared: `test-prot-07c-da-results-output-shared.R:347` failed on head. Its synthetic count
table uses the labels "Up" and "Down" (`:159`), not the package's; base's filter `status != "Not significant"` drew a
NumSig figure from them, the fix's `status %in% c("Significant and Up", "Significant and Down")` (7.1's text) drew
none. That test holds no stored reference, so stop rule 2 does not apply; it is a result caused by the change, and the
ticket's one correction (10) was spent on it.

Correction (29848a2f): the three filters (C4 guard and barplot, C6, C7) leave out `c("Not significant", "Not tested")`
instead of keeping the two significant labels. On every table the package makes the two are the same rows; on any other
table the figure is base's except that "Not tested" is left out, which is outcome 3's "those figures are unchanged". All
files are in 7.1 and no stored reference moved beyond the declared, so R3 was not to be repeated (10).

### R2c (29848a2f against dd769671): green

`r2c-compare.tsv`, `r2c-compare-failures.txt`, `r2c-check-compare.txt`, `r2c-kit204-head.log`. 28 same;
`test-omic-da-module.R` new-failure = exactly the 13 `prot_main` summary values ("Not tested (no value in a compared
group): 0" inserted before "Fold-change cutoff"); the new file 81 passed, 0 failed (only on head); check identical on
both sides, nothing only on head. Stop rules 2 and 3 not met (N4 passes on base, R1b, and head, R2c).

### R3 (full suite, 29848a2f against dd769671): green

`r3-full-suite-compare.tsv`, `r3-full-suite-compare-failures.txt`. 482 files: 480 same; `test-omic-da-module.R`
new-failure, exactly the 13 `prot_main` summary values (13 distinct steps, each the one line); the new file only on
head. 38 files fail identically on base and head (stored references recorded on the laptop, KIT-156; not this ticket's),
and `test-omic-da-module.R` fails on both sides with exactly these 13 more on head (53 → 66): of R3's failure
blocks, 965 are byte-identical on both sides, 13 are only on head and none only on base. Longest job 56.2 minutes (ceiling 70).

## Departure from ticket 6.4, by the coordinator

The implementation audit (AUDIT-KIT-204-impl-20261125, I02) found that the merge with staging conflicts in three files
under `R/` and one Rd page (KIT-220), so ticket 6.4 would ask for R3b on the merge. The coordinator ruled on
2026-11-25 that the next full-suite checkpoint on staging stands in for R3b (as for KIT-218's A4 and KIT-222's A5). R2b
on the merged tree is run (below, when recorded).

## R2b on the merge with staging (declared before sending)

Merge 89ce2a1f = `kit-204` (fdbb6751) + `v0.9-staging` 7327607f, the five conflicts resolved as
AUDIT-KIT-204-impl-20261125 section 6 states. R2b: base 7327607f (staging's head), head 89ce2a1f, `check: true`, 3
shards per side, ceiling 60 minutes per test job and 90 per check (36.0 core-hours, within the ticket's R2b row). Files:
R2's 30 (above) plus `test-kit220-da-threshold-consumers.R`, `test-kit218-da-wrapper-tail.R`, `test-kit212-da-current.R`,
`test-kit215-not-tested-lists.R` (KIT-249 and KIT-201 are not on staging). Predictions (the register's "the train must
hold"):
- `test-kit204-prot-da-not-tested.R`: only on head, 81 passed, 0 failed;
- `test-omic-da-module.R`: exactly the 13 `prot_main` `values` expectations newly failing, each by the one line, nothing
  else;
- `test-kit220-da-threshold-consumers.R`: `same`;
- the five GUI-entry lanes (`test-gui-entry-proteomics.R`): `same`;
- every other file `same`;
- R CMD check: findings identical, none only on head.
Decides: whether the merged tree holds 6.2's verdicts. R3b: the coordinator's ruling (the next full-suite checkpoint).

### R2b observed (45605914511, fd962ed9 against 7327607f): green, every prediction held

`r2b-merge-compare.tsv`, `r2b-merge-compare-failures.txt`, `r2b-merge-check-compare.txt`. 34 files: 32 `same`
(`test-kit220-da-threshold-consumers.R`, `test-kit218-...`, `test-kit212-...`, `test-kit215-...` and
`test-gui-entry-proteomics.R` among them); `test-omic-da-module.R` new-failure = exactly the 13 `prot_main` summary
values (13 distinct steps, each the one line; 211 → 198 passed); `test-kit204-prot-da-not-tested.R` only on head, 81
passed, 0 failed. R CMD check: 0 ERROR, 3 WARNING, 3 NOTE on both sides, none only on head. 248.4 job-minutes (16.6
core-hours); total for the ticket 132.0 core-hours of 262.7.

## The owner's decision VAL-278 (2026-11-26)

The 13 `prot_main` stored summary values of `tests/testthat/test-omic-da-module.R` (`:3485` to `:3973`; unchanged
on staging abfc99ac) are not edited and wait for the single end-of-project re-record; KIT-204 may merge after its
re-merge with staging and one comparison run (below). The ledger row `VAL-323` cites VAL-278.

## Second merge with staging (abfc99ac) and R2d (declared before sending)

**The merge.** d1bdcfe3 = `kit-204` 6e4f5994 + `origin/v0.9-staging` abfc99ac (KIT-239, 214, 226, 241, 224, 251, 252,
250, 244, 243, 201, 240, 245, 256, 259 and the full-suite checkpoint since 7327607f). Files changed on both sides
since 7327607f: `NEWS.md`, `R/func_prot_da_model_methods.R`, `R/utils_workflow_dia_da.R` (auto-merged) and
`docs/methods-decision-ledger.md` (the one textual conflict: staging's rows through VAL-279 kept, then KIT-204's row
`VAL-323`, now citing VAL-278 in its source column and its stored-references sentence). NEWS: KIT-204's entry
under the one unreleased heading, unchanged. No conflict under `R/`, `tests/`, `inst/` or `man/`.

**One shared function body, kept as merged (the coordinator's ruling of 2026-11-26).** In
`R/func_prot_da_model_methods.R`, method `differentialAbundanceAnalysisHelper` (`setMethod` at `:85`), git auto-merged
three hunks about 95 lines apart (merged-tree lines):
- KIT-201 `:414-416`: `if (is.null(ruv_model)) da_model <- c(da_model, workflowTreatModerationRecord())`, after the
  contrasts are run;
- KIT-204 `:512-518`: the first count call passes `p_value_column = !!rlang::sym(raw_pvalue_column)` after KIT-220's
  `q_val_thresh`;
- KIT-204 `:604-625`: the `"Not tested"` first `case_when` branch and the two barplot filters.
Reading (INFERRED, not run): KIT-201's line writes only `da_model`, which is read once, at `:462`
(`return_list$da_model`); KIT-204's hunks read and write only the count tables and their barplots. Neither reads a
variable the other writes. Guarded by the predictions below on `test-kit201-treat-moderation.R`,
`test-prot-ruv-factors-in-model.R` and the new file. `R/utils_workflow_dia_da.R`: KIT-243 changes
`runDIAReferenceDaAnalysis` and adds three functions; KIT-204 changes only `resolveDIAResultSelection`: no shared body.

**Scope, 37 files.** The ticket's grep on the merged tree (helpers dropped) finds 28 files: the 27 it found on
6e4f5994 plus `test-kit243-summary-current-outputs.R`. With the ticket's named files and R2b's additions, the scope is
R2b's 34 files plus `test-kit243-summary-current-outputs.R`, `test-kit201-treat-moderation.R` (KIT-201 changed the
shared method; named by the register's section 6) and `test-prot-ruv-factors-in-model.R` (the coordinator's guard).
KIT-249 is not on staging (its file is not run). `test-prot-norm-core.R` is not in scope, so only
`MULTISCHOLAR_TRX_MANIFEST_NOTE=1`.

**R2d.** Base abfc99ac (staging), head = the merged `kit-204` with this declaration, `check: true`, 3 shards per side,
`extra_packages: ["sva"]` (as the checkpoint, whose counts for `test-prot-ruv-factors-in-model.R` ran with `sva`),
reference BLAS. Ceiling 600 job-minutes (40.0 core-hours); ticket total at most 172.0 of 262.7. Decides: whether the
merged tree is ready for a merge train.

**Predictions.** Staging's counts are the checkpoint's head counts (`checkpoint-20261126/compare.tsv`, run 45644467847
on 828dd105) or merge train 10's (`merge-train-20261126c/compare.tsv`, run 45650328815 on 96e72da2) where the file is
there; since 828dd105 staging changed nothing in these files' code paths (KIT-259: normalisation figures; docs), so
the counts carry over (INFERRED).

| File | Base (staging) | Head | Verdict |
|---|---|---|---|
| `test-kit204-prot-da-not-tested.R` | absent | 81/0 | only-on-head, passing in full: the Not tested counts unaffected by KIT-201's moderation (INFERRED: `raw_pvalue` is missing only for untested rows under either moderation) |
| `test-omic-da-module.R` | 1289/80 | 1288/81 | new-failure, see below |
| `test-kit201-treat-moderation.R` | 89/0 | 89/0 | same (the TREAT moderation record still written) |
| `test-prot-ruv-factors-in-model.R` | 312/0 | 312/0 | same |
| `test-workflow-saved-da.R` | 776/175 | 776/175 | same, messages identical |
| `test-gui-entry-proteomics.R` | 168/25 | 168/25 | same |
| `test-omic-enrich-module.R` | 1011/59 | 1011/59 | same |
| `test-omic-summary-module.R` | 1675/28 | 1675/28 | same |
| `test-kit243-summary-current-outputs.R` | 246/0 | 246/0 | same |
| `test-kit220-...`, `kit218` (2 skipped), `kit212`, `kit215` | 79, 80, 309, 112 passed, 0 failed | as base | same |
| every other file (24) | the checkpoint's head count | as base | same |
| R CMD check | 0 E, 3 W, 3 N | as base | findings identical, none only on head |

`test-omic-da-module.R`, the 13 `prot_main` summary values, checked one by one against staging's failure list
(KIT-201's `compare-failures.txt` of R2b 45607796537 and R3 45592673648; the checkpoint reproduces their count,
`prot_main` 189 passed): KIT-201 already fails 12 of the 13 `values` expectations on staging, at element 8 (the `<ui>`
digest): `dendrogram_off`, `download`, `save_heatmap`, `select`, `table_down`, `table_limit`, `tree_dynamic`,
`tree_height`, `tree_k`, `tree_none`, `volcano_labels`, `volcano_ticked`. On head each of these 12 still fails, its
message gaining KIT-204's element-5 difference (the one line before "Fold-change cutoff"): message-changed, no count
change. The 13th, `volcano_unticked values`, passes on staging and newly fails on head by exactly the one line. So:
head 1288 passed and 81 failed (`(b) in process: prot_main` 189 → 188 passed); nothing else in the file moves.
The step attribution is read from the failure records, where each expectation's label follows its message.

**Stop rule.** Any other difference is a finding: recorded and reported, not fixed or re-run. An unfavourable valid
result is not re-run; an instrument defect with no result may be fixed and re-sent on the brief's four conditions, at
most three times.

### R2d observed (45666442643, d792847f against abfc99ac): ready for a merge train; one base-side finding

`r2d-merge-compare.tsv`, `r2d-merge-compare-failures.txt`, `r2d-merge-check-compare.txt`. Every head count is as
predicted.

| File | Predicted | Observed | Held |
|---|---|---|---|
| `test-kit204-prot-da-not-tested.R` | only on head, 81/0 | only on head, 81/0 | yes |
| `test-omic-da-module.R` | 1289/80 → 1288/81; 12 values message-changed, `volcano_unticked values` newly failing | 1289/80 → 1288/81; `prot_main` 189 → 188 passed; the 12 values gain exactly the one line beside KIT-201's `<ui>` difference, `volcano_unticked values` fails by exactly the one line; every other block identical | yes (the tool labels the file `message-changed`, not `new-failure`, because the `prot_main` test already fails on base) |
| `test-kit201-treat-moderation.R` | same 89/0 | same 89/0 | yes |
| `test-prot-ruv-factors-in-model.R` | same 312/0 | same 312/0 | yes |
| `test-workflow-saved-da.R`, `test-gui-entry-proteomics.R`, `test-omic-enrich-module.R`, `test-omic-summary-module.R` | same at 776/175, 168/25, 1011/59, 1675/28 | same at those counts | yes |
| `test-kit243-...`, `kit220`, `kit218`, `kit212`, `kit215` | same 246, 79, 80, 309, 112 passed | same at those counts | yes |
| `test-workflow-figure-paths.R` | same 181/0 | base 179/2, head 181/0 (`fixed`) | **no: a finding on the base side** |
| the other 23 files | same at the checkpoint's counts | same at those counts | yes |
| R CMD check | identical | 0 E, 3 W, 3 N both sides; 6 findings identical, none only on head or only on base | yes |

**Finding (base side only; recorded, not re-run).** On base (staging abfc99ac), two expectations of
`test-workflow-figure-paths.R` ("the figure check names the missing files and, on Windows, the over-long path") failed
with `expect_no_match(..., "259")`. The messages quote the R session's temporary folder, `/tmp/RtmpZBB2x8/file13648562592ac`,
whose random name contains "259". The over-long-path warning itself was absent. Head passed 181/0, as did staging in
the checkpoint and in train 10. Reading (INFERRED): the test is flaky whenever a random temporary name contains
"259". KIT-204 does not touch this file or the code it tests. A test-mechanics defect on staging for the coordinator;
not this ticket's to fix.

Cost: 252.5 job-minutes (16.8 core-hours) of the 600 ceiling; longest job 45.7 minutes. Ticket total 148.8 core-hours
of 262.7.

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A1 | met | R1b: each failing-first test fails at its stated value by an assertion; guards N3 (ii), N3 (iii), N4, N8 pass. R2c, R3: the new file 81 passed, 0 failed on head |
| A2 | met (after the correction) | R2c `r2c-compare.tsv`: 6.2's verdicts only (the 13 values); no compare-error; the only new-failure is the declared one. R2 had one undeclared new failure, corrected (above) |
| A3 | met | R2 and R2c: check findings identical base and head, none only on head |
| A4 | met | R3: 480 same, the declared 13 values, the new file only on head |
| A5 | met | NEWS under the one unreleased heading (Fixes); ledger row `VAL-323`; `DESCRIPTION`, lock files and `NAMESPACE` unchanged; `man/countStatDaGenes.Rd` and `man/printCountDaGenesTable.Rd` hand-edited only |
| A6 | met | "For the change report" below |

Plan section 11 row (E5, from R3): "| KIT-204 | `tests/testthat/test-omic-da-module.R`, `omicDaW1Frozen[["prot_main"]]`:
the `values` of 13 steps (`:3485` to `:3973`), each gaining the line 'Not tested (no value in a compared group): 0'
before 'Fold-change cutoff' (13 failing expectations); no written file and no GUI-entry reference | on merge | R3
45592061384 |".

## Deviations from the ticket's text (recorded)

1. `countStatDaGenes` does not add a first `case_when` branch `is.na({{ p_value_column }})`: with the default NULL that
   expression is `is.na(NULL)` and would break the unchanged path. It keeps today's `case_when` and then sets
   "Not tested" on the rows whose p-value is missing, only when a column is given. Same statuses, same counts; the
   NULL path is today's line for line.
2. The three barplot filters: `!status %in% c("Not significant", "Not tested")` instead of 7.1's
   `status %in% c("Significant and Up", "Significant and Down")` (the R2 correction above).
3. N3 reads the drawn statuses from each plot's data, asserting that `ggplot2::layer_data` has one row per data row
   (layer data hold positions, not labels).
4. N8 calls `printCountDaGenesTable` in the nested shape with the lipid/metab callers' other arguments
   (`formula_string = NA`, no p-value column). Their literal shape (a data frame per element) gives the shared helper
   an unnamed list whose comparisons cannot be named (INFERRED; the same on base and head, not this ticket's).
5. N9's manifest has a third contrast, C-B (R1's instrument repair).

## Recorded, not changed (two lines each)

1. **Lost q-value fallback (pre-existing, outside this ticket).** In `runTestsContrasts`, when `qvalue()` fails, the
   fallback `q_values_all[valid_p_idx] <- p.adjust(...)` is written inside the `tryCatch` error handler
   (`R/func_prot_da_model_stats.R:680-686` and `:803-809` at dd769671, the topTreat and topTable branches), a local copy, so
   every q-value stays NA and nothing is called significant; `qvalue_failed` is set with `<<-` and is recorded. On
   KIT-200's 60-protein fixture B-A (six proteins with a +1 effect, residual sd 0.3, 4 v 4) R1 and R1b show every
   tested B-A row with a missing q-value and 0 significant, while C-A and C-B are adjusted (INFERRED from the outputs
   and the code; the run's messages were silenced). `adjustProtDaResultTable` (`:379`) assigns the fallback correctly.
   Recommended: a ticket of its own (the DA tab is affected: N6's tab route also shows 0 significant on B-A).
2. The legacy proteomics wrappers' and the lipid/metab legacy count tables still say "Not significant" (S3, KD-B; F3,
   F7), as decided.

## For the change report

**What was wrong, who was affected.** In proteomics differential abundance a protein that could not be tested in a
comparison (no value in a compared group, or since KIT-200 no value in any analysed sample) has no p-value and no
q-value. The DA tab's summary counted it only in "Total genes", and the tables of significant counts (the stored
`num_sig_da_molecules_first_go` and `num_sig_da_molecules`, the published `*_num_sig_da_molecules.tab`/`.xlsx`, and
`outputDaAnalysisResults`' files) called it "Not significant", stating a test that did not happen. Every proteomics
analysis with such a protein in a comparison (application and the `differentialAbundanceAnalysis()` method; DIA, LFQ,
TMT, limma and limpa routes) was affected. Lipidomics and metabolomics already show the line (KIT-172).

**Code at base dd769671 and at the fix 29848a2f.**

| Place | Base | Fix |
|---|---|---|
| Display counts | `R/utils_workflow_dia_da.R:1290` (token v1), `:1299`, `:1304` | `:1294` (v2), `:1303` (`not_tested`), `:1308`, `:1316` |
| Tab summary | `R/utils_workflow_dia_da_views.R:287-290`, `:305-312` | `:293` (manifest), `:308-317` (list) |
| Count function | `R/func_prot_da_results_summary.R:13-42` | `:17-56` (`p_value_column`, `:29`, `:39`, `:47`); helper `:71`; printer `:178`, filter `:198` |
| S4 method | `R/func_prot_da_model_methods.R:508-512`, `:597-615` | `:508-514` (p column passed), `:602` (Not tested), `:616`, `:619` (filters) |
| Publication figure | `R/func_prot_da_results_methods.R:360-361` | `:359-361` |

Legacy (`main` 436cdbe8): the same counting in `R/func_prot_da.R:3037-3066` (`TRUE ~ "Not significant"` at
`:3051`) and `:4292-4301` (`:4297`); the tab's summary without the line in `R/mod_prot_da.R:1899-1906`.

**What changed in results.** For a comparison with n untested proteins (n > 0): the count tables gain a
"Not tested" row of n and their "Not significant" count falls by n; significant, up and down counts are unchanged; the
every-status barplot draws a "Not tested" bar; the significant-only figures and the report's sentences are unchanged.
The tab's summary always shows "Not tested (no value in a compared group): n", also when n is 0. For an analysis with
no untested protein every table and file is unchanged (R2c, R3: no written file's digest and no GUI-entry reference
moved); only the summary gains its ": 0" line (13 frozen display values, plan section 11). Analyses saved by an
earlier version keep their written count tables until the analysis is run again.

**Evidence.** R1b 45576574709 (failing first); R2c 45585667100 and R3 45592061384 (base against head). Synthetic
data only: KIT-200's 60-protein, three-group fixture (seed 200001) with one protein missing in group C and one
missing everywhere. Limits: no real data; Windows not run; the GUI-entry lanes hold no untested protein (no digest
moved), so the published `.tab` change is shown by N9 on fixtures only.

**On the merged tree (d792847f, staging abfc99ac merged).** This section still holds. The behaviour is unchanged and
R2d 45666442643 reproduces the new file's 81/0 and the same 13 summary values. The same places sit at different
lines on the merged tree:
- display counts: `R/utils_workflow_dia_da.R:1604` (v2), `:1613`, `:1618`, `:1626`;
- tab summary: `R/utils_workflow_dia_da_views.R:302`, `:318`;
- count function: `R/func_prot_da_results_summary.R:27`, `:39`, `:47`, filter `:202`;
- S4 method: `R/func_prot_da_model_methods.R:512-518`, `:607`, `:621`, `:624`;
- publication figure: `R/func_prot_da_results_methods.R:366-368`.

KIT-201's TREAT moderation now runs before these counts on the legacy route. It changes q-values but not which rows
lack a p-value, so the "Not tested" counts are the same under either moderation (INFERRED; R2d's new-file pass is
consistent with it).

## Not checked

- Real data; Windows.
- The rendered tab in a browser (the summary is checked as the text `diaDaSummaryStatsView` returns).
- The limpa (dpcDE) route's count tables with an untested protein (by reading, the same code; no limpa fixture here).
- Why `qvalue()` fails on B-A of this fixture (its messages were silenced by the test); the cause is inferred.
- The merge with staging after KIT-220 (merged first, VAL-230): `printCountDaGenesTable` has KIT-220's
  `q_val_thresh` and then this ticket's `p_value_column` (`countStatDaGenes` and `countStatDaGenesHelper` are not
  touched by KIT-220); three files under `R/` and one Rd page conflict, so the merger's rule (ticket 6.4) asks for
  R2b and R3b on the merge.
- The second merge (abfc99ac): no full suite on the merged tree. R2d covered the 37 files only; the shared method
  `differentialAbundanceAnalysisHelper` was judged by reading plus R2d's guard files.
