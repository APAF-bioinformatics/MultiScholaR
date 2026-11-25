# Merge train 13 of 2026-11-26: KIT-254, stacked on train 12

**This train is stacked on train 12.** Its base is train 12's tip, `origin/merge-train-20261126e` at
`eb6c40006e5a57e8e15f84d555391898223a3b8c`. That tip carries both of train 12's merges, KIT-249 (758d23ee) and KIT-236
(55413b59), and its declaration commit. Train 12 is itself stacked on train 11. Train 11 has landed (staging 0a774e0b);
train 12 has not, and its run is out.

**This train lands only if train 12 lands with its code unchanged.** If train 12 changes or does not land, this train
must be rebuilt.

## Train

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-254 | `origin/kit-254` bb015ef8 | a32edeb1 | `VAL-<train13-1>` (the fix), `VAL-<train13-2>` (the replay statements); `VAL-<train13-c>` reserved for the coordinator's correction of VAL-247 (not written here) | AUDIT-KIT-254-impl-20261126, closure check |

- The run binds head `a32edeb1`; the full hash is in the request.
- This file, the predictions table and the records commit are under `docs/qa/` only.
- The rows are numbered by the coordinator at landing.

### Conflicts met and their resolutions

**`R/mod_omic_norm.R`.** This is the one authorised code resolution, M2 against KIT-257, which is on the base.
- KIT-257's `omicNormOwnFiles` is kept whole.
- Then a closing `}` and a blank line.
- Then KIT-254's `omicNormCorrelationStart` and `clearOmicNormCorrelation`, unchanged. The shared `}` that follows
  closes `clearOmicNormCorrelation`.
- The merged file differs from KIT-254's side only by KIT-257's own edits: `own_files`, `remove =`, the two
  `*_unremoved` texts, `removed_then_stopped`, the removal in `omicNormApplyRunStart`, and `omicNormOwnFiles`.

```r
        sprintf("%s_ruv_corrected_%s.png", rep(safe, each = 4L), c("pca", "density", "rle", "correlation")))
}

#' KIT-254: before a correlation Apply or Skip at a correlation-step head ...
omicNormCorrelationStart <- function(manager, normData, addLog) {
...
clearOmicNormCorrelation <- function(normData) {
    ...
    invisible(normData)
}
```

**`NEWS.md`.** An append: the base's entries (KIT-259, KIT-257, KIT-249, KIT-236 and the rest) first, then KIT-254's,
under the one unreleased heading.

**`docs/methods-decision-ledger.md`.** An append. The base's rows are kept as they are: VAL-273 to VAL-280,
`VAL-<train11-1>`, `VAL-<train12-1>` and `VAL-<train12-2>`. KIT-254's two rows follow with clean status cells. Its
`VAL-<merge>` and `VAL-<merge-2>` placeholders became `VAL-<train13-1>` and `VAL-<train13-2>`, in the ledger and in
`kit-254/outcome.md`. VAL-247's correction is not made here.

**No other conflict.** `DESCRIPTION`, `NAMESPACE` and `man/` are unchanged.

### Path lengths

Checked by script before sending: 36 paths are new since the base. The longest is 76 bytes with the `MultiScholaR/`
prefix. None exceeds 100.

## Declaration (written before the run is sent)

**Purpose.** KIT-254's merge check on the merged tree (closure check, "The train's scope").

**Scope.** 62 files, base and head, R CMD check on both sides:
- R2c's 58 files (`kit254-r2c-904a5514.json`): the ticket's grep, the three frame files, the new file,
  `test-kit214-sample-columns.R` and nine replay files;
- the tenth replay file, `test-kit201-treat-moderation.R`;
- `test-kit243-summary-current-outputs.R` and `test-kit252-managed-failure-text.R`, which never ran with KIT-254;
- `test-kit257-norm-own-files.R`, because KIT-257 is in the base.

**Run settings.**
- Environment on both sides: `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1`
  (`test-prot-norm-core.R` is in scope).
- No extra packages, as KIT-254's runs.
- 4 shards per side, with the files listed longest first by the checkpoint's seconds.
- Per-file kill limit of 120 minutes for `test-prot-norm-core.R`. It is not a threshold.

**Ceiling.** 660 job-minutes: KIT-254's R2 ceiling (8 + 2 jobs), ticket section 10. The ticket's E-2 and R3 row make
the train's run its merge check.

**Base values.** Per file, in `merge-train-20261126f-predictions.tsv`.
- Train 12's declared head values where its scope holds the file.
- Otherwise train 11's declared values (checkpoint 45644467847 on 828dd105, plus KIT-257's file).
- In the scope here: `test-kit257-norm-own-files.R` 133/0 and `test-prot-norm-core.R` 1548/306, with KIT-257's four
  moved labels.
- `test-workflow-saved-da.R`, which KIT-249 moves to 755/196, is not in this scope.
- Every file has a value on record; none is read from the run alone.

### Predictions (closure check)

- `test-kit254-omics-correlation-skip.R`: only on head, 379 passed, 0 failed.
- `test-omics-norm-state-guards.R`: `count-changed`, 647/0 to 651/0 (the T1/T8 set-up edit; A4).
- Every other file `same` at the base's counts, including:
  - `test-prot-norm-core.R` 1548/306, with failures identical;
  - `test-workflow-analysis-record.R` 351/0;
  - `test-kit243-summary-current-outputs.R` 246/0;
  - `test-kit252-managed-failure-text.R` 72/0;
  - `test-kit257-norm-own-files.R` 133/0;
  - `test-kit201-treat-moderation.R` 89/0;
  - `test-workflow-omics-binding-steps.R` 58/1.
- No replay consumer gains a correlation statement or a `workflow_replay_correlation_skip` warning: per-file warning
  counts are equal base against head in every file but the new one.
- R CMD check: findings identical on both sides, nothing only on head, and the same portable-file-names path count on
  both sides.

### Decision and attempts

- Every prediction holds: TRAIN GREEN, conditional on train 12 landing unchanged.
- Any other verdict is a result: it is reported with its files and counts, not repaired here.
- One run. A dispatch failure (HTTP 500 or 502, no run created) is started once by hand with the README's command, and
  recorded.

## Result

Run 45669839790 bound head a32edeb1 against base eb6c4000, train 12's tip (request `merge-train-20261126f-a32edeb1`,
run-tests commit d16a4dd5). It was dispatched normally, sent once, and concluded success. It used 570.7 job-minutes of
the 660 ceiling, 28.8 of them for the check. The tables are in `merge-train-20261126f/`: `compare.tsv`, `tests.tsv`,
`tests-base.tsv`, `check-compare.txt`, `check-summary.txt`, `meta.txt`.

Of the 62 files, 60 are `same` and 2 are not. Both were predicted:

| File | Base p/f/e/s | Head p/f/e/s | Verdict | Prediction |
|---|---|---|---|---|
| test-kit254-omics-correlation-skip.R | absent | 379/0/0/0 | only-on-head | 379/0: held |
| test-omics-norm-state-guards.R | 647/0/0/0 | 651/0/0/0 | count-changed | 647 to 651: held |

- **Counts.** Every file's base and head counts equal the declared expectations in
  `merge-train-20261126f-predictions.tsv`, with no exception. For example: `test-prot-norm-core.R` 1548/306,
  `test-workflow-analysis-record.R` 351/0, the kit243 file 246/0, the kit252 file 72/0, the kit257 file 133/0, the
  kit201 file 89/0, `test-workflow-omics-binding-steps.R` 58/1.
- **Replay statements.** The per-file warning counts are identical base against head in all 61 files present on both
  sides, so no replay consumer gained a correlation statement or a `workflow_replay_correlation_skip` warning. Only
  the new file's log carries the statements.
- **R CMD check.** 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; all 6 findings identical, none only on head. The
  portable-file-names NOTE lists 238 paths on each side.

**Verdict:** every prediction held. TRAIN GREEN.

**Train 12.** It landed while this run was out, as bbbfdc1a (KIT-249 VAL-295, KIT-236 VAL-297, VAL-296). eb6c4000 is
an ancestor of staging, and `git diff eb6c4000 bbbfdc1a` outside `docs/` is empty. So train 12 landed with its code
unchanged, and the condition this train was built on holds.

**At landing.** The train lands by a merge. Its ledger still carries `VAL-<train11-1>`, `VAL-<train12-1>` and
`VAL-<train12-2>` in their pre-landing forms, which staging numbers VAL-281, VAL-295 and VAL-297. Staging's rows win
that conflict. The coordinator numbers `VAL-<train13-1>` and `VAL-<train13-2>` and writes `VAL-<train13-c>`, the
correction of VAL-247.

**Not checked:**
- the trees after the run's head (docs only);
- files outside the 62;
- the owner's wording items, C2 and C3.
