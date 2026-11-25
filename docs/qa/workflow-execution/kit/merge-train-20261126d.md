# Merge train 11 of 2026-11-26: KIT-257, stacked on train 10

**This train is stacked on train 10.** Its base is train 10's tip, `origin/merge-train-20261126c` at
`53e9955408dbcb30906a6c200756056a9895b542`: staging 828dd105 plus KIT-259, merge 4baea49f, with train 10's declaration.
Train 10 had not landed when this train was built; its run was still out.

**This train lands only if train 10 lands with its code unchanged.** If train 10 changes or does not land, this train's
evidence does not hold, and the train must be rebuilt.

Staging was abe2199c when this train was built: 828dd105 plus the checkpoint's docs. Those docs reach this branch
through KIT-257's own merge of abe2199c (ebd3c039), as additions; they did not conflict.

## Train

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-257 | `origin/kit-257` ebd3c039 | 1292afc4 | `VAL-<train11-1>`; `VAL-<train11-1b>` reserved for the coordinator's stored-value row (not written here) | AUDIT-KIT-257-impl-20261126, closure check |

The run binds head `1292afc4` (full hash in the request). The only commits after it are this file and the records
commit, both under `docs/qa/`. The rows are numbered by the coordinator at landing.

The merge commit's message summarises the change as "on every proteomics route". It is wider than that: every omic,
per the NEWS entry. The message was left unamended.

### Conflicts met and their resolutions

**`R/mod_prot_norm_workflow_helpers.R`.** This was the one authorised code resolution, at the KIT-259 hunk (register
closure check, "KIT-259's hunk"). The conflict block was replaced by the register's exact text:
- KIT-257's `fail` handler is kept.
- KIT-259's `withProtNormFigureNotice` wrapper and comment are kept.
- The two lines after the block were already KIT-259's.
- No other line was edited.

The merged file differs from KIT-257's side only by KIT-259's three edits: the comment and wrapper in `pass_through`,
this hunk, and `publish`'s wrapper.

Resolved text:

```r
      if (!is.null(runner)) dispatchProtNormNormalization(runner, settings, workflowData, normData, experimentPaths,
        omicType, experimentLabel, generatePostNormalizationQcFn, generateRuvCorrectedQcFn, getRuvGroupingVariableFn,
        finish, fail, preQc = generatePreNormalizationQcFn) else
        # KIT-259 (B2): after the completion notice, one warning names the figures that could not be drawn.
        withProtNormFigureNotice(finish(runDIAReferenceNormalization(settings, workflowData, normData, experimentPaths,
          omicType, experimentLabel, generatePostNormalizationQcFn, generateRuvCorrectedQcFn,
          getRuvGroupingVariableFn, pre_qc = generatePreNormalizationQcFn)), showNotificationFn)
```

**`docs/methods-decision-ledger.md`.** It is an append:
- The base's rows are kept through VAL-272 and `VAL-<train10-1>`.
- KIT-257's row follows as `VAL-<train11-1>`, with a clean status that names `VAL-<train11-1b>` for the two stored
  values.
- Its `VAL-<merge>` in `kit-257/outcome.md` was replaced to match.

**`NEWS.md`.** It merged without conflict. KIT-259's entries come first, then KIT-257's, under the one unreleased
heading.

**No other conflict.** `DESCRIPTION`, `NAMESPACE` and `man/` are unchanged. M2 (KIT-254) does not arise, because
KIT-254 is not in the base.

### Path lengths

Checked by script before sending:
- 29 paths are new since the base, most of them the checkpoint's docs.
- The longest is 87 bytes with the `MultiScholaR/` prefix.
- None exceeds 100.

## Declaration (written before the run is sent)

**Purpose.** KIT-257's merge check on the merged tree, from its register's closure check, "The train: scope and
predictions".

**Scope.** 81 files, base and head, with R CMD check on both sides:
- R2's 55 files (`kit-257/r2-files.txt`);
- the 4 declared dependents and the 3 from I10: `test-kit129-run-path.R`, `test-prot-design-log2.R`,
  `test-workflow-native-session.R`, `test-workflow-saved-da.R`, `test-trx-da-without-ruv.R`,
  `test-kit196-da-no-value-features.R`, `test-module-ci-trx-da.R`;
- `test-kit243-summary-current-outputs.R`, `test-kit244-limpa-pass-through.R`, `test-kit252-managed-failure-text.R`;
- M1's 15 DA test files;
- `test-kit259-pearson-histogram-figure-failure.R`, because KIT-259 is in the base.

**Run settings.**
- Environment `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1` on both sides.
- No extra packages, as KIT-257's R2' and train 10.
- 6 shards per side.
- Per-file kill limits of 120 minutes for `test-prot-norm-core.R` and the four GUI-entry files. These are not
  thresholds.

**Ceiling.** 740 job-minutes. This is KIT-257's merge-check ceiling: its R2 row, 8 × 70 + 2 × 90, ticket section 10.
The ticket's E-5 makes the train's run its merge check.

**Base values.** These are the base's own counts from the latest run on record for 828dd105: the full-suite checkpoint
45644467847, head side. Train 10 moves one file in this scope: `test-kit259-pearson-histogram-figure-failure.R`, new
at train 10's head, 85/0/0 as train 10 declares (its run is out). Train 10 predicts `same` for
`test-metab-01k-s4-plot-pearson-characterization.R` and `test-workflow-saved-exports.R`; neither is in this scope.

The checkpoint installed `sva`, `limpa` and `mixOmics`; this run installs none. Skip counts in files that need those
packages may therefore read lower than the checkpoint's. The verdict, base against head in this run, is what is
predicted. Every file's expected counts are in `merge-train-20261126d-predictions.tsv`.

### Predictions (register closure check)

| File | Verdict | Expected |
|---|---|---|
| `test-kit257-norm-own-files.R` | only on head | 133 passed, 0 failed, 0 skipped |
| `test-prot-norm-core.R` | `message-changed` | 1548/306 on both sides, changing at exactly four labels: the `files` of W1 (a) and W2 for `dia_normalization_refused` and `dia_normalization_retry`, bbfad929 → 53f7545d and b51063f4. These are the stored values at `:767` and `:801`, declared under `VAL-<train11-1b>`. Nothing else changes. |

Every other file is `same` at the base's counts, including:
- `test-workflow-saved-da.R` 776/175, messages identical;
- `test-kit243-summary-current-outputs.R` 246/0;
- `test-kit244-limpa-pass-through.R` 204/0;
- `test-kit252-managed-failure-text.R` 72/0;
- `test-kit259-pearson-histogram-figure-failure.R` 85/0;
- M1's 15 files;
- the GUI-entry files.

R CMD check is predicted to give identical findings on both sides, with nothing only on head. The portable-file-names
NOTE should list the same count on both sides.

### Decision and attempts

- If every prediction holds, the train is GREEN, conditional on train 10 landing unchanged.
- Any other verdict is a result. It is reported with its files and counts, and is not repaired here.
- One run is sent. One re-send is allowed, only for an instrument defect where no result was judged.

## Result

(filled after the run)
