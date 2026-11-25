# Merge train 9 of 2026-11-26: KIT-240, KIT-245, KIT-256

**This train is STACKED on merge train 8.** It was branched from `origin/merge-train-20261126a` at its tip
`b24fc15984329100926bda0eac77f073e160e34e` (train 8's declaration commit; train 8's code head aceb7190: KIT-244,
KIT-243, KIT-201), whose own run was still out when this train was built. Staging was 42846c0e (trains 6 and 7
landed). **Train 9 lands only if train 8 lands unchanged.** If train 8's check fails, or train 8 is rebuilt or changed
in anything but its records and ledger numbers, this train is rebuilt and its check repeated.

## Train

Base: `b24fc15984329100926bda0eac77f073e160e34e` (train 8's tip, fetched before building).

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-240 | `origin/kit-240` acf726ea | e1dc2154 | `VAL-<train9-1>` | AUDIT-KIT-240-impl-20261125, closure check (2026-11-26) C4 |
| 2 | KIT-245 | `origin/kit-245` 7b376683 | 75c5ddf6 | `VAL-<train9-2>` (Part S), `VAL-<train9-3>` (product) | AUDIT-KIT-245-impl-20261125, closure check (f) |
| 3 | KIT-256 | `origin/kit-256` c3d1b946 | 700853e7 | `VAL-<train9-4>` | AUDIT-KIT-256-impl-20261126, closure check |

Head (bound by the run): `700853e7` (full hash in the request). This file and the records commit after the run are
under `docs/qa/` only.

### Conflicts met and their resolutions

- **`docs/methods-decision-ledger.md`, all three merges.** The base's rows kept in order (staging's through VAL-263,
  then train 8's `VAL-<train8-1>` to `VAL-<train8-4>`), then train 9's rows in train order.
  - KIT-240: its row `VAL-<merge>` written `VAL-<train9-1>`; the same placeholder in the status cells it edits
    (INT-003, VAL-221) and in `kit-240/outcome.md` (commit table) replaced likewise.
  - KIT-245: the Part S row `VAL-<train9-2>`, the product row `VAL-<train9-3>`; in both, "correction <runs>" filled
    with R1c 45617558634, W1'-head 45617907710, W1'-fix 45619494004, R2c 45626340616 (closure I06). "the train's
    run <id>" is left to fill after the run.
  - KIT-256: `VAL-<train9-4>`; also in `kit-256/outcome.md` (commit table, A5 row). Its line 144 (an instruction
    naming the placeholder) left as written.
- **`NEWS.md`, merges 2 and 3.** KIT-245's two Fixes entries after the base's (KIT-244's), KIT-256's after KIT-245's;
  every entry kept, under the one unreleased heading. NEWS merged automatically in merge 1.
- **No conflict under `R/`, `tests/`, `inst/`, `man/` or metadata.** `DESCRIPTION` and `NAMESPACE` unchanged by the
  train.

### Owed at landing, outside this branch

- Train 8's placeholders first (`VAL-<train8-1>` to `-4`, `-2b`), then train 9's `VAL-<train9-1>` to `-4`.
- KIT-240: the release-checklist entries and KIT-228's partner note (closure C2, not merge conditions).
- KIT-245: "the train's run <id>" in its two rows; N02 is a separate ticket.
- KIT-256: its row's status cell with R3 45631885056 (N2) and VAL-249's status cell, as in `kit-256/outcome.md`.

### Path lengths

Checked by script before sending: 105 paths are new since the base; none exceeds 100 bytes with the `MultiScholaR/`
prefix; the longest is 89 (`docs/qa/workflow-execution/kit/kit-240/b/checks-45608074480-139184126330.txt`).

## Declaration (written before the run is sent)

**Purpose.** Each member's merge check on the merged tree, once, on the train.

**Scope.** 58 files, deduplicated, base and head, R CMD check on both sides, and one script:
- KIT-240: the 27 files of its merged run I+T (`k240-merged-bded5eb2`, run 45608060716);
- KIT-245: the 13 files of R2c (`kit245-r2c-91b269b7`), as in closure (f);
- KIT-256: the 17 regression files of ticket section 10 and the new file (`kit256-r2-57238dcd`).
- Script on the head: `docs/qa/workflow-execution/kit/kit-240/default-identity.R`, arguments the train's base (full
  hash) and `kit240-identity`, as KIT-240's I+T requested it. kit-138's gate-off script is not in the request.

**Run settings.** `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` on both sides; no extra packages (none of the three members'
merge checks used any); reference BLAS; 6 shards per side; per-file kill limits of 120 minutes for the four GUI-entry
files (not thresholds).

**Ceiling.** 1,260 job-minutes: KIT-240's merge-check (I+T) ceiling 24.00 core-hours = 360; KIT-245's R2c 300;
KIT-256's R2 600.

**Base values.** The base is train 8's tip, not staging. Of the 58 files, 21 are also in train 8's scope; of those
train 8 declares a move only for `test-gui-entry-proteomics.R` (189/4 to 168/25 on its head, KIT-201; under train 8's
environment with `sva` and `MULTISCHOLAR_KIT105_W3_DUMP=1`, which this run does not set). Every other overlap is
declared `same` by train 8, so it is predicted here at the members' last measured values. Where a member measured
without the trx flag and the flag can change a count, the prediction is `same` with the count read from the run.

### Predictions, KIT-240 (closure C4)

- `test-kit240-integration-default.R` only on head, 92/0, 1 skipped (T9).
- `test-integration-gui.R` 401 to 405, `test-integration-model.R` 263 to 265, `count-changed`, 0 failed.
- The other 24 `same`: `-record` 167, `-results` 258, `-outputs` 183, `-partitioned-alignment` 253,
  `-removed-factors` 520 (3 skipped), the eight `*-mofa.R` files skipped only, `test-workflow-checkpoint-harness.R`
  628/0/1 error identical, `test-prot-05b-norm-module-contracts.R` 440/0 (22 skipped), the rest at I+T's values.
- Identity script: every check TRUE; per mode settings note 1, Results note 2, mark 21, Enrichment insertion 8
  (0 on base).

### Predictions, KIT-245 (closure (f))

- `test-kit245-github-push-safety.R` only on head 87/0, `test-kit245-github-panel.R` only on head 42/0; S6 run, not
  skipped.
- `test-github-management-direct-shared.R` 21/0 and `test-kit234-workbooks-withdrawn.R` 29/0 on both sides.
- `test-omic-summary-module.R` 1675/28 on both, the 28 failures identical (train 8 declares it `same`).
- `test-module-ci-trx-summary.R` `same`, count read from the run (the register's 83/2 was without the trx flag;
  train 4, with the flag, read 86/2). `test-optional-packages.R` 11/1, `test-omic-gui-config.R` 1911/0,
  `test-module-ui-wrapper-direct-shared.R` 45/0 on both; the four quarantined summary files 0/0.

### Predictions, KIT-256 (closure check)

- `test-kit256-tech-rep-selection-cleared.R` only on head, 20/0/0.
- The 17 `same`, failure messages identical: `test-gui-entry-transcriptomics.R` 16/13; `test-module-ci-trx-design.R`
  66/0/1 error; `test-workflow-ruv-pearson-diagnostic.R` 301/0/2 errors/5 skipped;
  `test-workflow-gui-memory-lifetimes.R` 81/0 (11 skipped); `test-workflow-omics-residency.R` 97/0; the rest at R2's
  values (`test-kit222-replicate-keys.R` 456/0, `test-prot-04-design.R` 488/0 (118 skipped), and so on).
- `test-gui-entry-proteomics.R`: `same`, count read from the run (expected near train 8's head 168/25, which carries
  KIT-201's moves; not derivable here because train 8's environment differs and its run is out).

### Predictions, R CMD check

- Findings identical on both sides, nothing only on head; the portable-file-names NOTE lists the same number of paths
  on both sides.

### Decision and attempts

- Every prediction holds: TRAIN GREEN, landing after train 8 lands unchanged.
- Otherwise the member is named with its files and counts; it is not repaired here.
- One run; one re-send only for an instrument defect with no result judged.

## Result

Run 45636060029 (request `merge-train-20261126b-700853e7`, run-tests commit c590695d), head 700853e7 against train 8's
tip b24fc159, sent once; conclusion success; 354.8 job-minutes of the 1,260 ceiling (tests 315.0, check 26.7, script
11.8, plan and results 1.2). Tables in `merge-train-20261126b/` (`compare.tsv`, `tests.tsv`, `tests-base.tsv`,
`check-compare.txt`, `check-summary.txt`, `meta.txt`, `identity-console.log`).

Observed, 58 files: 52 `same`, 6 not `same`:

| file | base p/f/e/s | head p/f/e/s | verdict | ticket, prediction |
|---|---|---|---|---|
| test-kit240-integration-default.R | absent | 92/0/0/1 | only-on-head | KIT-240, as declared |
| test-integration-gui.R | 401/0/0/0 | 405/0/0/0 | count-changed | KIT-240, as declared |
| test-integration-model.R | 263/0/0/0 | 265/0/0/0 | count-changed | KIT-240, as declared |
| test-kit245-github-push-safety.R | absent | 87/0/0/0 | only-on-head | KIT-245, as declared (S6 ran: 0 skipped) |
| test-kit245-github-panel.R | absent | 42/0/0/0 | only-on-head | KIT-245, as declared |
| test-kit256-tech-rep-selection-cleared.R | absent | 20/0/0/0 | only-on-head | KIT-256, as declared |

Counts read from the run, as declared:
- `test-gui-entry-proteomics.R` 168/25 on both sides, `same` (equal to train 8's declared head, KIT-201's moves).
- `test-module-ci-trx-summary.R` 87/1 on both sides, `same`.

Every other file is `same` at the predicted values (for example `test-omic-summary-module.R` 1675/28,
`test-gui-entry-transcriptomics.R` 16/13, `test-workflow-checkpoint-harness.R` 628/0/1 error,
`test-workflow-ruv-pearson-diagnostic.R` 301/0/2 errors/5 skipped, the eight MOFA files skipped only).

Identity script (KIT-240): exit 0, "all checks pass: TRUE"; in both modes settings note 1, Results note 2, mark 21,
Enrichment insertion 8 (0 on base); head FALSE = base unset; head unset is not the unavailable panel.

R CMD check: 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; all 6 findings identical, none only on head or only on base.

Verdicts:
- KIT-240: holds exactly.
- KIT-245: holds exactly.
- KIT-256: holds exactly.
- Check: identical.

Train verdict: TRAIN GREEN, on condition that train 8 lands unchanged (this run's base is train 8's tip b24fc159,
whose own run was out when this train was built).

Not checked: train 8's own result; any file outside the 58; KIT-240's MOFA lane and bundle lane (stand, not repeated);
the full suite (each member's A-criterion at the next checkpoint).
