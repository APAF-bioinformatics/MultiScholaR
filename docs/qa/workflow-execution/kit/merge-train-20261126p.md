# Merge train 23 of 2026-11-26: KIT-259 (second and third corrections)

**This train is NOT stacked.** It is branched from `origin/v0.9-staging` at `176f9d3ffa100f0480fd5100a32c002c5a5cdbe5`
(fetched before building; staging had not moved by code since: 176f9d3f is the ledger commit VAL-329 on top of train 18).
It does not contain trains 19 to 22 (`merge-train-20261126l`, `m`, `n`, `o`: other agents', still out). **KIT-271,
KIT-266 and KIT-272 depend on this landing** (KIT-266 needs `writeProtNormCompositePlaceholders` and the list on
KIT-257's removal list; KIT-272 needs the four templates' chunks).

## Train

Base: `176f9d3ffa100f0480fd5100a32c002c5a5cdbe5` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-259 (second and third corrections) | `origin/kit-259b` 2b0b379d | 460253ba | `VAL-330` (replaces `VAL-<merge>`) | AUDIT-KIT-259-impl-20261126, closure checks 1 and 2 (verdict: ready after listed test-or-records repairs; applied) |

Head (bound by the run): `460253ba5d2f7b7a2ec9defe6b5cb78fc12f9f9a` (`git merge --no-ff kit-259b`, which already
contains staging 176f9d3f; the ledger placeholder replaced in the same commit). This file and the predictions are in the
declaration commit after it, under `docs/` only; the records commit after the run likewise.

### Conflicts met and their resolutions

- None. kit-259b already contains staging 176f9d3f, so the merge was clean, ledger included.
- **`docs/methods-decision-ledger.md`**: staging's rows kept through VAL-329; KIT-259's row, `VAL-<merge>`, is now
  `VAL-330` (the only occurrence; the branch's records name no placeholder).
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`. Code files changed against the base:
  `R/func_general_filemgmt_results_summary_helpers.R`, `R/func_lipid_s4_correlation_methods.R`,
  `R/func_metab_s4_qc_methods.R`, `R/mod_prot_norm_workflow_helpers.R`, `R/utils_workflow_dia_norm.R`,
  `R/utils_workflow_dia_products.R`, `R/utils_workflow_dia_sessions.R`, the four proteomics templates under
  `inst/reports/proteomics/`, `tests/testthat/test-kit259-pearson-histogram-figure-failure.R` (the declared B-E1 edit)
  and the new `tests/testthat/test-kit259b-second-correction.R`. No `man/` or `NAMESPACE` change.

### Path lengths

Checked before sending: the longest path new since the base is 89 bytes with the `MultiScholaR/` prefix
(`docs/qa/workflow-execution/kit/kit-259/r2c3-45701181156/compare-failures.txt`), under the 100-byte limit; this
train's own files are shorter.

## Declaration (written before any run is sent)

**Purpose.** KIT-259's merge check on the merged tree (comparison with R CMD check), and the render job of closure
check C.4 (the proteomics reports from the four GUI lanes).

**Scope recount.** The comparison of run 45701181156 had 83 files (the second correction's 77 plus the third pass's
six). Recounted on THIS base with the ticket's two greps (R2's list: `plotPearson`, `rebuildWorkflowDiaQcPlot`,
`buildProtNormQcPlotProduct`, `buildProtNormCorrelationPlot`, `runDIAReferenceNormalization`,
`publishDIAReferenceNormalization`, `generateProtNormPostNormalizationQcArtifacts`,
`generateProtNormPreNormalizationQcArtifacts`, `dispatchProtNormNormalization`; A3's `generateLipidQcPlots`,
`generateMetabQcPlots`, `(lipid|metab)_qc_figures`), the four named files, and the third pass's additions
(`drawProtNormQcFigure`, `removeProtNormOwnFiles`, the proteomics templates, `collectDIANativeGuiFiles`,
`restoreDIANativeGuiFiles`, `copyToResultsSummary`, the combined figure's file names): **all 83 files of the old list
still qualify, and two files that landed since the old base 82a5df02 qualify as well**:
`test-kit204-prot-da-not-tested.R` (KIT-204, train 15; reads `DIANN_report.rmd`) and
`test-kit261-norm-run-button.R` (KIT-261, train 16; names `lipid_qc_figures`). **Scope = 85 test files** (83 + 2),
each of the two as `same` expected. The four proteomics templates are byte-identical between 82a5df02 and this base; staging's code changes since 82a5df02
(trains 15 to 18) touch `R/func_metab_s4_qc_methods.R` (KIT-217, other methods) and none of this branch's other R files
(merge: no conflict). `test-prot-norm-core.R` is in scope: no
`MULTISCHOLAR_KIT105_W3_DUMP`, as in run 45701181156.

**Run settings, copied from run 45701181156's request (`kit259c-r2-c6e8814d.json`):** `env` empty; no extra packages;
`blas: reference`; `shards: 5`; `file_timeout_minutes: 45`; `file_timeouts` 120 for `test-prot-norm-core.R` and
`test-gui-entry-proteomics.R`; `check: true`; no scripts; `tools_ref` empty; test list = that request's 83 plus the two
files above. Base = 176f9d3f, head = 460253ba. The render is a SEPARATE request (as `kit259c-render-c6e8814d.json`: one
scripts job, `docs/qa/workflow-execution/kit/kit-259/render3.R`, `NOT_CRAN=true`, `out` `kit259-render`, `out_pattern`
`\.(txt|tsv)$`, script timeout 55, no tests, no check, `blas: reference`); `render3.R` is used unchanged: it compares
the head's templates with `git show 82a5df02d:inst/reports/proteomics/...`, which equal this base's templates, so its
V0 comparison is with the base's text.

**Ceilings.** Comparison 600 job-minutes = 40.0 core-hours (the branch's own comparison cost 490.1 for 83 files).
Render 60 job-minutes = 4.0 core-hours (the branch's cost 23.2). Total ceiling 660 job-minutes = 44.0 core-hours.

### Predictions

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit259b-second-correction.R | absent | 90 passed / 0 failed | only-on-head |
| test-prot-norm-core.R | 1548 passed / 306 failed (train 18's tables; train 17's `merge-train-20261126i`) | 1547 / 307 | message-changed: exactly ONE new failure, "the normalisation action and the correlation action give the recorded output through the moved functions", the pinned call count at `:187-189` (`resolveProtNormCompositeFigureInputs` 1 against 2); left unedited by the coordinator's ruling of 2026-11-26, waits for the end-of-project re-record (VAL-282 / VAL-301 precedent). ANY other change in that file is a failed prediction |
| test-kit259-pearson-histogram-figure-failure.R | 85/0 | 85/0 | same |
| test-metab-01k-s4-plot-pearson-characterization.R | 17/0 | 17/0 | same |
| test-kit204-prot-da-not-tested.R, test-kit261-norm-run-button.R | read from run | same as base | same |
| every other file (79) | read from run | same as base, failure lists identical | same |
| R CMD check | | | identical findings (0 E / 3 W / 3 N as the branch's last run), nothing only on head |
| Render V0 | | the head's report text identical to the base's on prot_lfq, prot_tmt, prot_dia, prot_dia_limpa; drawings 3, 3, 3, 1 | |
| Render V1 | | exactly one added sentence on LFQ, TMT and DIA, naming panels a and h; none on DIA with limpa | |
| Render V2 | | a list with no combined figure gives no sentence | |

Totals expected: 85 files, 83 `same`, 1 `only-on-head`, 1 `message-changed` (test-prot-norm-core.R), no `new-failure`
verdict and no `count-changed` outside it; render `checks.tsv` 27 of
27 as predicted. (The new KIT-259b file is the one absent on base.) Expected about 500
job-minutes for the comparison.

### Decision and attempts

- Every prediction holds: TRAIN GREEN; landing needs the coordinator's ledger ruling below.
- Any other verdict is a result: recorded as it is, not re-sent, reported.
- One comparison run and one render run. A re-send only for an instrument defect with no result judged (VAL-54, at most
  three). An HTTP 500/502 at dispatch with no run created: started once by hand with the README's command, recorded.

## What the coordinator resolves at landing

- The ledger: `VAL-330` is KIT-259's row (`<id>`/`<t>` placeholders, if any remain in its text, fill at landing); trains
  19 to 22's rows are reconciled as they land.
- The stored-value row for `test-prot-norm-core.R:188` (the skip run's count of `resolveProtNormCompositeFigureInputs`,
  1 to 2; no output changes), as VAL-282 and VAL-301.
- `NEWS.md`: KIT-259's bullets (merged without conflict here).

## Result

Two runs, both dispatched without error (no HTTP 500/502; the API's 403 rate limit hit only the watcher, not the
runs). Requests `train23-compare-460253ba` (run 45716988827) and `train23-render-460253ba` (run 45716996784), request
commit 9d5ca71d. Head 460253ba vs base 176f9d3f (staging). Tables: `merge-train-20261126p/` (comparison) and
`merge-train-20261126p-render/` (render), same directory as this file.

**Job-minutes.** Comparison 580.9 (tests head 62.8/44.5/60.4/62.0/54.8, tests base 53.3/55.6/48.8/53.0/54.0, check head
15.1, check base 15.6, plan 0.3, results 0.7) against the ceiling of 600 = 38.73 core-hours. Render 23.7 (scripts 22.9;
the script itself 792 s) against 60 = 1.58 core-hours. **Total 604.6 job-minutes = 40.31 core-hours** (ceilings 660 =
44.0).

| Item | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| Scope | 85 files | | 85 | HELD (83 + `test-kit204-prot-da-not-tested.R`, `test-kit261-norm-run-button.R`) |
| test-kit259b-second-correction.R | only-on-head, 90/0 | absent | 90/0, 0 skipped | HELD |
| test-prot-norm-core.R | message-changed, 1548/306 to 1547/307, one new failure at `:187-189` | 1548/306 | 1547/307 | HELD (below) |
| test-kit259-pearson-histogram-figure-failure.R | same, 85/0 | 85/0 | 85/0 | HELD |
| test-metab-01k-s4-plot-pearson-characterization.R | same, 17/0 | 17/0 | 17/0 | HELD |
| the other 80 files incl. kit204 and kit261 | same | same | same | HELD |
| R CMD check | identical | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | HELD: 6 findings identical, only on head none, only on base none |
| Render V0 | text identical to base's, drawings 3, 3, 3, 1 | | | HELD, all four lanes |
| Render V1 | one added sentence (a and h) on LFQ, TMT, DIA; none on DIA with limpa | | | HELD |
| Render V2 | no sentence | | | HELD |

Comparison totals: 85 files, 83 `same`, 1 `only-on-head`, 1 `message-changed`, no `new-failure`, no `count-changed`.

**The one change in `test-prot-norm-core.R`**, from `compare-failures.txt` (head):
the test "the normalisation action and the correlation action give the recorded output through the moved functions",
kind failure: "Expected `skipped$calls` to be identical to `c(...)`. `actual[2:7]`: 0 0 1 2 1 1; `expected[2:7]`: 0 0 1 1 1 1"
(the pinned call count of `resolveProtNormCompositeFigureInputs`, `:187-189`, 2 observed against 1 pinned). The other
failures of that block (the `correlation`, `normalization_manual` and `normalization_skip` digests) are identical on
both sides. No other change in the file.

**Render** (`merge-train-20261126p-render/checks.tsv`): 27 of 27 checks as predicted. V0 text equal on prot_lfq,
prot_tmt, prot_dia, prot_dia_limpa; drawings 3, 3, 3, 1; no placeholder sentence. V1 one sentence, letters a and h, on
LFQ, TMT, DIA, the rest of the text equal to V0; none on DIA with limpa. V2 no sentence on any lane.

**Decision.** Every prediction held. TRAIN GREEN, with the one pinned-count message change that is the coordinator's
stored-value row (the end-of-project re-record list). The pre-existing failures in the other files (for example
`test-gui-entry-proteomics.R` 92/101, `test-omic-summary-module.R` 1675/28, `test-workflow-saved-da.R` 755/196) are
identical on both sides.
