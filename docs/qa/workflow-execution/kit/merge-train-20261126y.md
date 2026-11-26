# Merge train 32 of 2026-11-26: KIT-268 alone (the notice, report sentence and record field naming the comparisons that fell back to Benjamini-Hochberg), STACKED on merge train 31

**This train IS stacked on merge train 31** (KIT-272, branch `origin/merge-train-20261126x`, out and not yet landed). It
exists because KIT-268 and KIT-272 both edit `tests/testthat/test-kit236-report-sentences.R` and the report templates, and
the combination must be run.

## Base and head

- **Train 31 base**: branch tip `8105c2620` (its merge head `f211b6675`, KIT-272, plus docs-only commits: the register's
  record fixes, its declaration, the render-binding merge). `git diff --stat f211b6675 8105c2620 -- R tests inst man` is
  empty: the tip is code-identical to the merge head train 31's record names.
- **Staging**: `origin/v0.9-staging` at `5ab4e01bf` is ahead of train 31's base `ec27f69e0`. It brings merge train 28
  (KIT-278: `test-omics-module-globals.R` allows `mod_omic_norm.R`; KIT-269: `formatConfigList` lists vector settings in
  full, the writer's failure line; ledger VAL-360, VAL-361) and its records. Nothing newer than `5ab4e01bf` existed when I
  built.
- **B** (the comparison's base side) = `92205ca9393c404954d0c72fad43550687d85470`: `git merge --no-ff origin/v0.9-staging`
  on train 31's tip. One conflict, `docs/methods-decision-ledger.md` only (staging's rows VAL-360 and VAL-361 first, train
  31's row stays `VAL-<train31-1>`). NEWS merged by itself. Nothing under `R/`, `tests/`, `inst/` or `man/` conflicted.
  So KIT-272's change is on BOTH sides of the comparison.
- **Head (code-bearing)** = `c0f4b3f7681172cd1172d9a3d3310fd30e99cff3`: `git merge --no-ff origin/kit-268` (`49d1a5c83`) on B.
  The docs-only commit after it (`6383dcb12`, record fixes) and this file change no code, template or test.

## Train

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-268 (VAL-312, VAL-314, VAL-351) | `origin/kit-268` 49d1a5c83 | c0f4b3f76 | `VAL-<train32-1>` (the method row), `VAL-<train32-2>` (the stored-reference row); the coordinator numbers both at landing | AUDIT-KIT-268-impl-20261126 (closure check: READY FOR A MERGE TRAIN; C02 applied below) |

### Head differs from B by exactly KIT-268's change

`git diff --name-only 92205ca93 c0f4b3f76` lists, outside `docs/qa/`: `NEWS.md`, `docs/methods-decision-ledger.md`,
`R/func_omics_da_shared.R`, `R/mod_omic_da.R`, `R/mod_prot_da_handlers_run_analysis.R`, `R/utils_workflow_analysis_record.R`,
the six templates (`inst/reports/lipidomics/lipidomics_report.rmd`, `inst/reports/metabolomics/metabolomics_report.rmd`,
`inst/reports/proteomics/{DIANN_limpa,DIANN,LFQ,TMT}_report.rmd`), `tests/testthat/test-kit236-report-sentences.R`,
`tests/testthat/test-kit268-bh-fallback-notice.R`, `tests/testthat/test-omic-da-module.R`; and 169 files under
`docs/qa/workflow-execution/kit/kit-268/`. That is the closure check's list (C.2) and nothing else.
`git diff --name-only origin/kit-268 c0f4b3f76 -- R inst man tests` lists only what comes from staging's trains and train 31
(KIT-278, KIT-269, KIT-272, KIT-255 and earlier), plus the merged test-kit236 file.

### Conflicts met and their resolutions

- **`tests/testthat/test-kit236-report-sentences.R`** (expected): KIT-255's block (on staging at the end of the file) and
  KIT-268's T4 block both append at the same place. Both kept, KIT-255's first, then KIT-268's; git's alignment put KIT-255's
  loop body and KIT-268's header side by side, and the shared closing `    }` / `})` was moved so that each block is whole.
  KIT-272's edits (`k272OldBullet`, `k272NewBullet`, `k272Substitute`, in the file's earlier part) merged by themselves and
  are untouched. Merged file, top-level order of the added/relevant parts: `k272OldBullet`, `k272NewBullet`,
  `k272Substitute` (`:498-508`) ... `test_that("KIT-236 T22 (I09) ...")` (`:723`) ->
  `test_that("KIT-255 (I01, I02): the comment counts each sample set missing once ...")` (`:736`) ->
  `k268Existing`, `k268Added` (`:753-754`) -> `test_that("KIT-268 T4: ...")` (`:759`, last, ends the file). No duplicated
  top-level assignment (`grep -oE '^name <-' | sort | uniq -d` is empty); no name clash between `k255`-side locals,
  `k268*` and `k272*`.
- **`NEWS.md`**: train 31's first Fixes bullet (KIT-272, Figure 1 panels d to f) kept first, then KIT-268's bullet; the
  following bullets unchanged.
- **`docs/methods-decision-ledger.md`**: every row kept: staging's through VAL-361, then `VAL-<train31-1>` (KIT-272), then
  `VAL-<train32-1>` (KIT-268's method row, renamed from `VAL-<merge>`), then `VAL-<train32-2>` (the stored-reference row, the
  closure check's exact text, added in the docs commit `6383dcb12`). The same commit applies C02 to
  `VAL-<train32-1>`: the proteomics notice "driven by a test on the DIA (iq) route only (the W3 fixture)", the label-free,
  TMT and limpa routes established by reading the shared handler; and the report sentence of lipidomics and metabolomics, like
  their record, waits for KIT-242.
- No other file under `R/`, `tests/`, `inst/` or `man/` conflicted.

## Declaration (committed and pushed before any run is sent)

**Purpose.** The merge check of KIT-268 on the stacked tree: base B (staging 5ab4e01b + train 31) against head (B +
KIT-268), the closure check's scope with KIT-272's files added, R CMD check on both sides.

**Scope recount: 44 files.** The ticket's grep with `warnDaQvalueFallback` added
(`git grep -lE "qvalue_warnings|fdr_method|fdr_fallback_reason|workflowRecordDaFdrMethod|completeOmicDaRun|da_server_run_analysis_handler|warn_bh_fallback|ReportFdrSentence|Storey|Benjamini-Hochberg adjusted|warnDaQvalueFallback" -- 'tests/testthat/test-*.R' ':!data/benchmark'`)
gives **36** files on the merged tree (35 on staging plus the new KIT-268 file; the tests train 25, 27, 28 and KIT-272 added
do not match). Added: `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R`, `test-omic-summary-module.R` (the GUI
entries and the summary module, R2's 39), `test-omics-module-globals.R`, `test-kit196-da-no-value-features.R` (I03), and
KIT-272's `test-kit272-prot-box-panel-text.R`, `test-kit227-report-ruv-sentences.R`, `test-kit235-lipid-metab-report-parser.R`.
36 + 8 = 44; none of the eight is in the grep's 36. All 44 exist on the head.

**Run settings, copied from KIT-268's R2 request (`kit268-r2-e3660e02.json`)**: `tools_ref` empty; `extra_packages`
`sva`, `limpa`, `mixOmics`; `blas: reference`; `file_timeout_minutes` 45; `file_timeouts` 120 for the two GUI-entry files;
`check: true`; no scripts. Changed from R2: `shards` 6 (as R2; the list order is a longest-first deal so the shards balance:
six round-robin shards of 8, 8, 7, 7, 7, 7 files, 14.0 to 14.5 test-minutes each at the committed timings) and **`env`
`MULTISCHOLAR_TRX_MANIFEST_NOTE=1`** (R2's `env` was empty; transcriptomics corpus files are in scope through
`test-kit243-summary-current-outputs.R` and `test-omic-summary-module.R`, the standing rule since train 28).

**Timing arithmetic.** Committed per-file seconds (largest seen in the committed tables) sum to 5,172 s = 86.2
test-minutes per side. Per shard about 14.5 test-minutes plus library setup (about 9.5 minutes) and compile: about 26-30
job-minutes; the longest single file is 670 s (`test-omic-da-module.R`). Even at twice the committed timings no job reaches
39 minutes, so none passes 75. The estimate rests on R2's measurement: R2 measured 282.0 job-minutes
for 39 files, 6 shards per side and both checks; adding 5 files and KIT-272's three (about 20 test-minutes per side) gives
about 310. **Expected 310 job-minutes (20.7 core-hours); ceiling 420 job-minutes = 28.0 core-hours** (core-hours = job-minutes
x 4 / 60).

### Predictions

Base side = B, which already holds KIT-272, so KIT-272's changes are on both sides. Base counts are the newest committed
figures (the closure check's table; train 31's and KIT-272's records for the three KIT-272 files). **Where the run's own base
side differs, the run's base side is the reference**, and I say so here in advance.

| File | Base on B (p/f/e/s) | Head | Verdict |
|---|---|---|---|
| `test-kit268-bh-fallback-notice.R` | absent | 30/0/0/0 | only-on-head |
| `test-kit236-report-sentences.R` | 297/0/0/0 (staging's 293 incl. KIT-255 + KIT-272's 4) | 321/0/0/0 (+24: KIT-268's T4) | count-changed |
| `test-omic-da-module.R` | 1288/81/0/0 | 1288/81/0/0, identical failures; A112-03 passes on both sides | same |
| `test-gui-entry-lipidomics.R` | 19/16/0/0 | 19/16/0/0; at step `da` only the failure message also shows the added warning | message-changed |
| `test-gui-entry-metabolomics.R` | 57/48/0/0 | same | same (a `message-changed` on metab_gc was once seen: the counts are the claim) |
| `test-omic-enrich-module.R` | 1011/59/0/0 | same; a one-side `stringdb_plot` digest difference in "(b) in process: refusals" is the known flake I06: recorded, no re-run | same |
| `test-omics-module-globals.R` | 14/0/0/0 | same | same |
| `test-kit196-da-no-value-features.R` | 540/0/0/0 | same | same |
| `test-omic-summary-module.R` | 1675/28/0/0 | same | same |
| `test-workflow-analysis-record.R` | 351/0/0/0 | same | same |
| `test-workflow-saved-da.R` | 755/196/0/0 | same | same |
| `test-workflow-operation-registry.R` | 1107/2/0/0 | same | same |
| `test-kit220-da-threshold-consumers.R` | 71/2/0/0 | same | same |
| `test-kit127-omics-ruv-factors.R`, `test-dia-limpa-ruv-rules.R` | 425/1/0/0, 333/1/0/0 | same | same |
| `test-tech-reps-limma.R` | 0/0/1/0 | same (an error on both) | same |
| `test-kit272-prot-box-panel-text.R` | 373/0/0/0 | 373/0/0/0 | same (KIT-272 on both sides) |
| `test-kit227-report-ruv-sentences.R` | 186/0/0/0 | 186/0/0/0 | same |
| `test-kit235-lipid-metab-report-parser.R` | 58/0/0/0 | 58/0/0/0 | same |
| the other files of the 44 (`test-kit243` 246, `test-kit249` 144, `test-kit174` 77, `test-go-propagation` 321/0/0/1, `test-enrichment-significance` 181, `test-omic-gui-config` 1911, `test-module-ci-lipid-da` 140, `test-module-ci-metab-da` 138, `test-module-ci-prot-da` 250/0/0/1, `test-ruv-control-anova-grouping` 139, `test-workflow-omics-da-contrasts` 138, `test-workflow-prot-declarations` 303/0/0/2, `test-prot-*` and `test-lipid-*`, `test-metab-*` characterization files at the checkpoint's counts) | as committed in the newest table | as base | same |
| R CMD check | 0 ERRORs, 3 WARNINGs, 3 NOTEs (train 31's figures) | identical findings, nothing only on head | identical |

Totals expected: 44 files; 1 `only-on-head`, 1 `count-changed`, 1 `message-changed`, 41 `same`; no `new-failure`, no
`fixed`.

**Masked comparison for the lipid GUI-entry file.** From `compare-failures.txt`, for each of the 16 failing steps compare
the base and head messages after masking temporary and checkout paths; for the step `da`, additionally delete from the head
message the one added notification (type `warning`, "Storey's q-values could not be computed for these comparisons, so
Benjamini-Hochberg adjusted p-values are shown in the q-value column instead: KO_vs_WT in LCMS_Pos; KO_vs_WT in
LCMS_Neg."). Predicted: after that deletion all 16 messages are identical, the other 15 steps' messages are identical
without any deletion, and `metab_*` files have no non-identical message.

### Decision and attempts

- Every prediction holds: TRAIN GREEN for KIT-268 on the stacked tree; landing needs the coordinator's ledger numbers and
  status cells.
- Any other verdict is a result: recorded as it is, not re-sent, handed back. A mechanics failure with no result: handed
  back with the evidence; no re-send without the coordinator's word. A 500/502 at dispatch with no run created: started once
  by hand (`gh workflow run run-tests.yml --ref run-tests -f request=<file>`), recorded.
- If train 31 itself is not green, this run still stands as an observation of the stacked tree; the coordinator decides
  the landing.

## What the coordinator resolves at landing

- Land train 31 first, then this train (it is stacked on it; its base commits are train 31's). Then: the numbers of
  `VAL-<train31-1>` (KIT-272), `VAL-<train32-1>` (KIT-268) and `VAL-<train32-2>` (the stored-reference row), in that order.
- Status cells: VAL-312, VAL-314 (3) and VAL-351 (KIT-268 implemented/merged by train 32).
- Placeholders: `docs/qa/workflow-execution/kit/kit-268/outcome.md:42` (a record's mention of the ledger row `VAL-<merge>`), the
  ledger's `VAL-<train32-1>` and `VAL-<train32-2>` rows, and this record's `VAL-<train32-n>` mentions (train 31's
  `kit-272/outcome.md:65` is train 31's to resolve).
- `tickets/kit/KIT-268.md` status line still reads "R2 next" (closure check C03, coordinator bookkeeping).

## Results (run 45741255942, request `train32-compare-c0f4b3f7`, base B `92205ca93`, head `c0f4b3f76`)

The run is green (tests success, check success). Tables: `merge-train-20261126y-compare.tsv` and
`merge-train-20261126y-check-compare.txt` (copies of the run's `compare.tsv` and `check-compare.txt`). One run, no re-send.
`meta.txt` confirms the request as sent: 44 files, 6 shards, `check` true, `env` `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`.

**Verdicts: 44 files: 41 `same`, 1 `count-changed`, 1 `message-changed`, 1 `only-on-head`; no `new-failure`, no `fixed`. One
base count was missed (below), identically on both sides. R CMD check identical.**

| File | Predicted (base / head) | Observed (base / head) | Verdict |
|---|---|---|---|
| `test-kit268-bh-fallback-notice.R` | absent / 30/0/0/0 | absent / 30/0/0/0 | only-on-head: HELD |
| `test-kit236-report-sentences.R` | 297 / 321 (0 failed) | 297/0/0/0 / 321/0/0/0 | count-changed: HELD (+24 = T4) |
| `test-omic-da-module.R` | 1288/81 same | 1288/81/0/0 both sides, identical failures | same: HELD |
| `test-gui-entry-lipidomics.R` | 19/16 both, message-changed at `da` only | 19/16/0/0 both | message-changed: HELD (masked comparison below) |
| `test-gui-entry-metabolomics.R` | 57/48 same | 57/48/0/0 both | same: HELD |
| `test-omic-enrich-module.R` | 1011/59 same | 1011/59/0/0 both | same: HELD (no `stringdb_plot` difference appeared; I06 not seen) |
| `test-omics-module-globals.R` | 14/0 same | 14/0/0/0 both | same: HELD |
| `test-kit196-da-no-value-features.R` | 540/0 same | 540/0/0/0 both | same: HELD |
| `test-omic-summary-module.R` | 1675/28 same | 1675/28/0/0 both | same: HELD |
| `test-workflow-analysis-record.R` | 351/0 same | 351/0/0/0 both | same: HELD |
| `test-workflow-saved-da.R` | 755/196 same | 755/196/0/0 both | same: HELD |
| `test-kit272-prot-box-panel-text.R` | 373/0 both | 373/0/0/0 both | same: HELD |
| `test-kit227-report-ruv-sentences.R` | 186/0 both | 186/0/0/0 both | same: HELD |
| `test-kit235-lipid-metab-report-parser.R` | 58/0 both | 58/0/0/0 both | same: HELD |
| `test-kit220-da-threshold-consumers.R` | 71/2/0/0 both | **79/0/0/0 both** | same, but the BASE COUNT WAS MISSED: MISSED prediction (identical on both sides; fewer failures than the committed 71/2) |
| the other 29 files | as base counts of the closure check | as predicted, each `same` | HELD |
| R CMD check | 0/3/3 identical | base 0 ERRORs, 3 WARNINGs, 3 NOTEs; head the same; 6 findings identical, 0 only on head, 0 only on base | identical: HELD |

**The one miss: `test-kit220-da-threshold-consumers.R` reads 79 passed, 0 failed on BOTH sides, against the declared
71/2/0/0.** Observed fact: the counts are identical on B and on the head, so KIT-268 changes nothing in this file. Inferred,
not shown: the committed 71/2 comes from requests without `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`; the file uses the
transcriptomics declaration reference (`trxDeclReference`, `:204`), whose corpus-manifest digest check the variable turns into
a note (the same effect train 32 saw in train 28's three files, `merge-train-20261126u.md`). My declaration took the base
counts from tables that were run without the variable while setting the variable in this request; I did not catch that the
two diverge for this file. A prediction miss, recorded; no re-run.

**Masked comparison, `test-gui-entry-lipidomics.R`** (`compare-failures.txt`, one block, base side and head side of the 16
failing steps; performed offline with Python on the fetched file, no R). Method: split each side into its 16 step blocks
(`lipid_canonical at step '<id>'`), delete from the head's `da` block the 12 lines of the added notification's difference
(from "`actual$notifications` is length 4" through "`expected$notifications[[4]]` is absent": the length 4 against 3, the
type `warning` against `message` at position 3, the message "Storey's q-values could not be computed for these comparisons,
so Benjamini-Hochberg adjusted p-values are shown in the q-value column instead: KO_vs_WT in LCMS_Pos; KO_vs_WT in
LCMS_Neg." against "DA results saved to disk (tables, volcano plots, heatmaps)", and the extra fourth notification), then
compare block by block. Result: both sides have the same 16 steps; **all 16 blocks are identical after that one deletion**
(15 identical without any deletion; the `da` block's remaining lines, the state-digest differences, are identical). The
only other difference is the position of the `da` block in the concatenated message (first on the head, third on the base),
not its content; the text before the first block is identical. The head side has 8 mentions of `notifications` against 0 on
the base, all inside the deleted chunk. The notification text is the one stored in the closure check's text for
`VAL-<train32-2>`: it matches exactly. `metab_*` of `test-gui-entry-metabolomics.R`: verdict `same`, no non-identical
message.

**Cost (jobs API).** 276.07 job-minutes = 18.40 core-hours, against the declared expected 310 and ceiling 420 (28.0). Longest
job 24.7 minutes (tests base 3), far under 75. Shard timings are in the run; no job passed 25 minutes.

**Decision.** Every declared prediction held except the one base count of `test-kit220-da-threshold-consumers.R`
(identical on both sides, favourable). By the declaration's rule that is a result, recorded, not re-run. KIT-268 on the
stacked tree (staging 5ab4e01b + train 31 + KIT-268): no regression; the stacked combination of KIT-272 and KIT-268 in
`test-kit236-report-sentences.R` passes 321/0; the lipid GUI-entry oracle moves at step `da` by the one added warning only.
