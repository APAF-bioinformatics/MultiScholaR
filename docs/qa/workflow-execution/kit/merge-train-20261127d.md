# Merge train 38 of 2026-11-27: KIT-219 (Copy and the Report present only the enrichment files the current enrichment made)

**Not stacked.** Branched from `origin/v0.9-staging` at `e4a6d4aea` (fetched before building; train 37 has landed: KIT-242 VAL-385, KIT-260 VAL-386, VAL-387). Branch `merge-train-20261127d`,
worktree `/home/doktersmol/Projects/APAF/MultiScholaR-train38`, base `e4a6d4aea`.

## Train

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-219 | `origin/kit-219` 3743790a (holds staging a67c4161, merge base with staging) | 21ed5db2b | `VAL-388` (the ticket's row); the coordinator numbers it at landing (next free VAL-388). No stored-reference row is needed: the audit found that no stored value moves (no stored capture, digest or reference changes; the GUI-entry lanes store the report as `<present>`) | AUDIT-KIT-219-impl-20261126 (ready after rulings and record fixes, which are done) |

Head (bound by the run): `21ed5db2b`, the merge commit; later commits are under `docs/` only.

### Conflicts

`git merge --no-ff origin/kit-219` conflicted in `docs/methods-decision-ledger.md` only (one hunk at the end of the VAL rows). `NEWS.md` merged without conflict; `inst/reports/lipidomics/lipidomics_report.rmd` and `metabolomics_report.rmd` auto-merged;
nothing under `R/`, `tests/`, `inst/` or `man/` conflicted. Ledger: the branch's copy of VAL-376 to VAL-383 is older (staging has since filled their status cells), so staging's rows (VAL-376 to VAL-387) are kept whole and the member's row follows as `VAL-388` (line 2560).
(A slip while resolving: my first resolution script aborted on a wrong line index and I committed the merge with the markers still in; I noticed at once, fixed the file and amended that unpushed merge commit, `21ed5db2b`; nothing was pushed with markers.)

### The merged tree against the base

`git diff --name-only e4a6d4aea 21ed5db2b` outside `docs/qa/` is exactly the member's own change, 12 files, none extra: `R/func_general_filemgmt_results_summary_helpers.R` (+38 -4), `R/utils_workflow_summary.R` (+126 -4), `man/copyToResultsSummary.Rd` (+9 -2),
`inst/reports/lipidomics/lipidomics_report.rmd` (+12 -7), `metabolomics_report.rmd` (+12 -7), `proteomics/DIANN_limpa_report.rmd` (+6 -6), `DIANN_report.rmd` (+3 -3), `LFQ_report.rmd` (+3 -3), `TMT_report.rmd` (+3 -3),
`tests/testthat/test-kit219-enrichment-files-current.R` (new, 619 lines), `NEWS.md` (+13), `docs/methods-decision-ledger.md` (+1 row); plus 36 files under `docs/qa/workflow-execution/kit/kit-219/`. The insertion and deletion counts of every file equal the member's own against its merge base (a67c4161).
Staging's changes since the member's merge base (84 commits, 42 files under R, tests, inst, man; of the files KIT-219 touches only the two report templates changed after a67c4161): the lipidomics template hunks at `:429` (+3 lines), `:543` (+6) and `:552` (KIT-242/KIT-260 `FormatMethodDetails`) and the metabolomics template hunks at `:417` (+3), `:532` (+6), `:541`
sit beside KIT-219's hunks at `:1630-1654` and `:1586-1610` (the enrichment-folder block); both sets of hunks are present (the merged file differs from the base by KIT-219's five hunks only). KIT-273's Copy branch (`writeDIARuvDeliveryTables` at `R/func_general_filemgmt_results_summary_helpers.R:457`) and KIT-216's and KIT-266's changes to the proteomics templates and `R/utils_workflow_summary.R`
are in the member's merge base and intact beside KIT-219's hunks (`@@ -12`, `-22`, `-30`, `-365`, `-719`, `-721`, `-948` of the helpers file; `-226`, `-430`, `-434`, `-450`, `-452`, `-455`, `-527` of `utils_workflow_summary.R`); the merged diff of both R files against the base has only those hunks.

## Declaration (committed and pushed before the request)

**Purpose.** The merge check of KIT-219 on staging e4a6d4aea (which holds KIT-242 and KIT-260): base against head, R CMD check on both sides. One run.

**Scope: 51 files** (list below, in the request order; `merge-train-20261127d-predictions.tsv`). The ticket's formula as the member's records state it (`predictions.md`, "Run R2"; ticket section 10): `git grep -lE "copyToResultsSummary|copyWorkflowSummaryResults|renderWorkflowSummaryReport|pathway_dir|_enrichment_results|Enrichment_Plots|Pathway_enrichment_results|diaEnrichmentProduct" -- tests/testthat`,
helpers dropped (40 files on the merged tree), plus the four GUI-entry lane files, `test-e2e-proteomics-enrichment.R`, `test-prot-summary-report-templates.R`, `test-go-propagation.R`, `test-kit212-da-current.R`, `test-kit216-prot-session-after-enrichment.R`, the new file and the eight template-naming tests = **50, the same list as the member's `r2-predictions.tsv`**.
**My count equals the member's: KIT-242 and KIT-260 added no file that the formula matches** (their new files do not contain these names). Plus `test-omics-module-globals.R` for continuity (the diff touches no `R/` file matching `^(mod|func)_(lipid|metab|omic|omics)...R$`, so the rule does not require it) = **51**.

**Settings** (the form of train 36's request): 5 shards per side; `file_timeout_minutes` 45; `file_timeouts` 90 for `test-kit266-limpa-combined-qc-figure.R`, 120 for the four GUI-entry files (`test-kit247-design-resave.R` and `test-prot-norm-core.R` are not in scope); `env` `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`;
`check` true on both sides; `blas` reference; `os` ubuntu-22.04; `extra_packages` none; no scripts. Files are dealt to shards by position; the order below puts the heavy files in different shards.

Planned minutes per shard (test minutes from train 37's committed seconds; the new file 121 s from the member's head records; job minutes = test minutes + about 10 setup):

| Shard | Test min | Job min | Files |
|---|---|---|---|
| 1 | 29.8 | 39.8 | kit266-limpa-combined-qc-figure, general-filemgmt-report-helpers-shared, general-filemgmt-path-contracts, module-ci-prot-summary, lipid-14a-summary-module-characterization, omics-module-globals, metab-05-summary-module-characterization, e2e-proteomics-enrichment, module-ci-lipid-summary, module-ci-sentinels, kit226-gprofiler-columns |
| 2 | 29.5 | 39.5 | gui-entry-proteomics, omic-summary-module, kit205-enrichment-background, enrichment-significance, workflow-dia-products, omic-gui-config, module-ci-metab-summary, omic-design-module, general-enrichment-shared, module-ci-prot-enrich |
| 3 | 28.8 | 38.8 | kit212-da-current, kit257-norm-own-files, kit271-limpa-table-name, workflow-analysis-record, kit215-not-tested-lists, kit219-enrichment-files-current, kit208-gprofiler-identifiers, module-ci-metab-da, multiomics-enrich-shared, prot-11a-enrichment-module-characterization |
| 4 | 28.5 | 38.5 | kit216-prot-session-after-enrichment, kit273-limpa-ruv-delivery, gui-entry-metabolomics, kit207-summary-head-guard, kit213-study-parameters-labels, omics-enrich-core, gui-entry-transcriptomics, workflow-gui-diagnostic, lipid-15-summary-module-seams, prot-summary-report-templates |
| 5 | 28.5 | 38.5 | kit243-summary-current-outputs, omic-enrich-module, kit259b-second-correction, omics-summary-core, workflow-omics-figures, gui-entry-lipidomics, go-propagation, workflow-memory-mechanism-member, prot-11-enrichment-module-contracts, prot-12-summary-module-contracts |

Largest shard 29.8 test minutes, about 40 job minutes of 300 (13 %), far under the 75 % bound (225). **Ceiling 1,400 job-minutes = 93.3 core-hours** (core-hours = job-minutes x 4 / 60). **Expected about 425 job-minutes (28.3 core-hours)**: 2 x (145 + 5 x 10) = 390, checks 2 x 16 = 32, plan and results about 3.

### Predictions (passed/failed/errors/skipped)

Base counts come from train 37's full-suite head table (`merge-train-20261126aa-tests-head.tsv`; staging e4a6d4aea holds that tree), not from older scoped trains. A base-side count that differs from it is REPORTED, not explained away.

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit219-enrichment-files-current.R | absent | 157/0/0/0 | only-on-head |
| test-omic-summary-module.R | 1675/28 | 1675/28 | same |
| test-omics-summary-core.R | 676/25/0/2 | same | same |
| test-omic-enrich-module.R | 1011/59 | same (a one-side `stringdb_plot` digest difference in its "(b) in process: refusals" is the known flake: recorded, not repaired) | same |
| test-gui-entry-proteomics.R 92/101; test-gui-entry-lipidomics.R 19/16; test-gui-entry-metabolomics.R 57/48; test-gui-entry-transcriptomics.R 16/13 | as stated | as stated, messages identical | same |
| test-kit216-prot-session-after-enrichment.R | 152/0/0/1 | same | same |
| test-kit212-da-current.R 309/0; test-kit243-summary-current-outputs.R 246/0 | as stated | same | same |
| test-omics-module-globals.R | 14/0 | 14/0 | same |
| the other files | train 37's table | same | same |
| R CMD check | 0 ERROR, 3 WARNING, 3 NOTE | identical | nothing only on head |

Totals expected: 51 files; 1 only-on-head, 50 same; no new-failure, no fixed, no message-changed, no count-changed. Any difference outside this declaration is a RESULT: reported, not repaired.

### Scope (51 files)

- test-e2e-proteomics-enrichment.R
- test-enrichment-significance.R
- test-general-enrichment-shared.R
- test-general-filemgmt-path-contracts.R
- test-general-filemgmt-report-helpers-shared.R
- test-go-propagation.R
- test-gui-entry-lipidomics.R
- test-gui-entry-metabolomics.R
- test-gui-entry-proteomics.R
- test-gui-entry-transcriptomics.R
- test-kit205-enrichment-background.R
- test-kit207-summary-head-guard.R
- test-kit208-gprofiler-identifiers.R
- test-kit212-da-current.R
- test-kit213-study-parameters-labels.R
- test-kit215-not-tested-lists.R
- test-kit216-prot-session-after-enrichment.R
- test-kit219-enrichment-files-current.R
- test-kit226-gprofiler-columns.R
- test-kit243-summary-current-outputs.R
- test-kit257-norm-own-files.R
- test-kit259b-second-correction.R
- test-kit266-limpa-combined-qc-figure.R
- test-kit271-limpa-table-name.R
- test-kit273-limpa-ruv-delivery.R
- test-lipid-14a-summary-module-characterization.R
- test-lipid-15-summary-module-seams.R
- test-metab-05-summary-module-characterization.R
- test-module-ci-lipid-summary.R
- test-module-ci-metab-da.R
- test-module-ci-metab-summary.R
- test-module-ci-prot-enrich.R
- test-module-ci-prot-summary.R
- test-module-ci-sentinels.R
- test-multiomics-enrich-shared.R
- test-omic-design-module.R
- test-omic-enrich-module.R
- test-omic-gui-config.R
- test-omic-summary-module.R
- test-omics-enrich-core.R
- test-omics-module-globals.R
- test-omics-summary-core.R
- test-prot-11-enrichment-module-contracts.R
- test-prot-11a-enrichment-module-characterization.R
- test-prot-12-summary-module-contracts.R
- test-prot-summary-report-templates.R
- test-workflow-analysis-record.R
- test-workflow-dia-products.R
- test-workflow-gui-diagnostic.R
- test-workflow-memory-mechanism-member.R
- test-workflow-omics-figures.R

### Decision and attempts

One run. A shard that times out or an instrument failure with no result: handed back, not re-sent. A 502/403 from the API is not a run failure (re-attach once). An HTTP 500/502 at dispatch with no run created: started once by hand.

## Results (run 45842971860, request train38-kit219-21ed5db2, head 21ed5db2b, base e4a6d4aea)

Run success; every file ran to completion (no timeout, no crash, no withheld file); the watch did not need a re-attach. Tables beside this file: `merge-train-20261127d-compare.tsv`, `-check-compare.txt`, `-tests-head.tsv`, `-tests-base.tsv`.
Base versus head, 51 files: same 50, only-on-head 1, new-failure 0, count-changed 0, message-changed 0, fixed 0. Tests head: passed 10225, failed 310, errors 3, skipped 45; base: 10068, 310, 3, 45 (the difference 157 is the new file).
**Every prediction held: TRAIN GREEN for KIT-219.**

| File | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| test-kit219-enrichment-files-current.R | absent / 157/0 | absent | 157/0/0/0 (226 s) | only-on-head: HELD |
| test-omic-summary-module.R, test-omics-summary-core.R, test-omic-enrich-module.R | 1675/28, 676/25/0/2, 1011/59 | as predicted | same | same: HELD (the `stringdb_plot` flake did not appear: 0 mentions in `compare-failures.txt`) |
| the four GUI-entry files | 92/101, 19/16, 57/48, 16/13, messages identical | as predicted | same | same: HELD |
| test-kit216-prot-session-after-enrichment.R 152/0/0/1; test-kit212-da-current.R 309/0; test-kit243-summary-current-outputs.R 246/0; test-omics-module-globals.R 14/0 | same | as predicted | same | same: HELD |
| the other files | same | equal to train 37's table | same | same: HELD |
| R CMD check | identical 0/3/3 | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | identical, 6 findings on both, none only on head: HELD |

No base-side count differed from train 37's table (all 50 existing files compared by script; the only mismatch the script printed is the new file's "absent" against its 0/0/0/0 placeholder row).

**Budget** (jobs API, summed, skipped jobs 0): **469.55 job-minutes = 31.30 core-hours** (core-hours = job-minutes x 4 / 60). Ceiling 1,400 (93.3), expected about 425 (28.3): 10 % above the expectation, 34 % of the ceiling. Longest shard 56.82 job-minutes (`tests base 1`, with `test-kit266-limpa-combined-qc-figure.R`, which took 2,536 s on head, limit 90 min): 19 % of 300. Check head 13.38, base 9.68; plan 0.42; results 0.53.

Decision: KIT-219 is green on the merged tree.

## What the coordinator resolves at landing

- Ledger: `VAL-388` (`docs/methods-decision-ledger.md:2560`) takes its number (next free VAL-388). Status cells of the rows it cites: VAL-179 (SD-5) and VAL-187: fill "implemented by KIT-219 (VAL-n), merged by train 38" where they wait for this merge. No stored-reference row (the audit found no stored value moves).
- Placeholders in the member's records and NEWS: `docs/qa/workflow-execution/kit/kit-219/outcome.md:50` (row A5) mentions the ledger row as `VAL-<merge>` (a record mention, to read as the row's number); `NEWS.md` has none. (My first line here said none; the grep output showed that one hit.)

## Merge record (coordinator)

- Landed on `v0.9-staging` on 2026-11-27 on top of e4a6d4ae, the train's own base: the merge applied without conflict.
- Ledger: the ticket's row is VAL-388; status cells of VAL-179 (SD-5) and VAL-187 filled. No stored-reference row: no stored value moves.
- The staged change outside `docs/` equals the member's own change (two R files, the hand-edited Rd file, six report templates, the new test file, NEWS).
- Run 45842971860 on the merge head 21ed5db2: 51 files and the check; 50 same; the new file 157/0 on head only; no base-side count differed from merge train 37's table; check identical; 469.55 job-minutes (31.30 core-hours) of 1,400.
- The builder reports one slip that never reached the remote: its first merge commit held conflict markers in the ledger; it was corrected before any push.
- Open with the owner (not part of this train): the sentence the proteomics report prints when an enrichment is no longer current.
