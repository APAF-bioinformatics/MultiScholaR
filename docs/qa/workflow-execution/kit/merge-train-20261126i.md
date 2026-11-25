# Merge train 16 of 2026-11-26: KIT-261, stacked on train 15

**This train is stacked on train 15.** Its base is train 15's tip, `origin/merge-train-20261126h` at
`069c2832b0eae98d023cdd7e58855748c8e9942a`. That tip already holds all three of train 15's merges, KIT-190 (588c960c),
KIT-204 (ab5ba116) and KIT-217 (4308c5ae), and its declaration commit. Train 15 is itself stacked on train 14, which
has landed (staging d7cfe160). Train 15 had not landed when this train was built.

**This train lands only if train 15 lands with its code unchanged.** If train 15 changes or does not land, this train
must be rebuilt.

## Train

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-261 | `origin/kit-261` e7623fb1 | 5ff6cfcf | `VAL-<train16-1>` | AUDIT-KIT-261-impl-20261126 (ready; I01 scope, I02 ledger) |

The run binds head `5ff6cfcf` (full hash in the request). The only later commits are this file, the predictions table
and the records commit, all under `docs/qa/`. The row is numbered by the coordinator at landing.

### Conflicts met and their resolutions

**`docs/methods-decision-ledger.md`, the only conflict.** It was resolved as I02 prescribes:
- Every row of the train's base is kept in its order. These are the base's forms of VAL-282 and of the train 13 to 15
  rows (`VAL-<train13-1>` … `VAL-<train15-3>`).
- KIT-261's row follows as `VAL-<train16-1>`.
- The status cell carries I02's text, with the train's run to fill.

KIT-261's branch also carried staging's numbered VAL-298 to VAL-301, which came in through its merge of d7cfe160. They
were not carried here, because the base holds those rows in its own (placeholder) forms. The landing reconciles them
with staging.

Staging's in-place correction of VAL-247 (train 13's landing) did not conflict, and it came in with KIT-261's merge.
`kit-261/outcome.md`'s `VAL-<merge>` was replaced.

**Merged without conflict:**
- `NEWS.md`: KIT-261's entry is under `### Fixes` of the one unreleased heading, after the base's entries.
- `R/mod_omic_norm.R`.
- The staging docs KIT-261 brought: the train 13 and 14 records, as additions.

There was no conflict under `R/`, `tests/` or `inst/`. `DESCRIPTION`, `NAMESPACE` and `man/` are unchanged.

### Path lengths

Checked by script before sending: 25 paths are new since the base. The longest is 83 bytes with the `MultiScholaR/`
prefix. None exceeds 100.

## Declaration (written before the run is sent)

**Purpose:** KIT-261's merge check on the merged tree (AUDIT-KIT-261-impl I01).

**Scope: 44 files, base and head, R CMD check on both sides.**
- The ticket's grep (10.3), recounted at this train's base 069c2832, gives 35 files. These are R2's 34 plus
  `test-kit254-omics-correlation-skip.R`. No file train 15 added matches the grep.
- The 8 helper-driven files: `test-gui-entry-{lipidomics,metabolomics,transcriptomics}.R`,
  `test-workflow-omics-worker.R`, `test-workflow-omics-residency.R` and `test-module-ci-trx-{norm,da,summary}.R`.
- `test-kit261-norm-run-button.R`.

**Run settings:**
- `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1` on both sides.
- Extra packages `limpa`, `mixOmics` and `sva`, as train 15, so that the base counts declared by train 15 apply.
- 4 shards per side, with files listed longest first.
- Per-file limits of 120 minutes for `test-prot-norm-core.R` and the three GUI-entry files. These are kill limits, not
  thresholds.

**Ceiling: 820 job-minutes.** This is KIT-261's R2 ceiling (8 + 2 jobs, ticket section 11), which serves as its merge
check.

**Base values.** `merge-train-20261126i-predictions.tsv` gives each file's base value from the first available source:
1. train 15's declared head values, which apply to 32 of the files;
2. otherwise trains 14 to 11's declared values;
3. otherwise the checkpoint 45644467847 on 828dd105.

Every file has a value on record. Train 15 moves three files in this scope: `test-gui-entry-lipidomics.R` 19/16,
`test-gui-entry-metabolomics.R` 57/48 and `test-omic-da-module.R` 1288/81. Its declared heads are used as this base.

### Predictions (AUDIT-KIT-261-impl I01)

- `test-kit261-norm-run-button.R`: only on head, 105 passed, 0 failed.
- Every other file `same`, with identical failure messages, at the base's counts in the table. For example:
  - `test-prot-norm-core.R` 1548/306;
  - `test-kit254-omics-correlation-skip.R` 379/0;
  - `test-kit257-norm-own-files.R` 133/0;
  - `test-omics-norm-state-guards.R` 651/0;
  - `test-gui-entry-transcriptomics.R` 16/13.
- R CMD check: findings identical on both sides, nothing only on head. The portable-file-names NOTE has the same count
  on both sides.

### Decision and attempts

- Every prediction holds: TRAIN GREEN, conditional on train 15 landing unchanged.
- Any other verdict is a result. It is reported with its files and counts, and not repaired here.
- One run. A dispatch failure (HTTP 500 or 502, no run created) is started once by hand and recorded.

## Result

(filled after the run)
