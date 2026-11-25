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

Run 45688253564 (request `merge-train-20261126i-5ff6cfcf`, run-tests commit 01b040e5): head 5ff6cfcf against base
069c2832 (train 15's tip). Sent once, dispatched normally, conclusion success. It used 430.2 job-minutes of the 820
ceiling (check 29.1). Tables are in `merge-train-20261126i/`: `compare.tsv`, `tests.tsv`, `tests-base.tsv`,
`check-compare.txt`, `check-summary.txt` and `meta.txt`.

Of the 44 files, 42 are `same` and 2 are not.

| File | Base p/f/e/s | Head p/f/e/s | Verdict | Prediction |
|---|---|---|---|---|
| test-kit261-norm-run-button.R | absent | 105/0/0/0 | only-on-head | 105/0: held |
| test-gui-entry-metabolomics.R | 57/48/0/0 | 57/48/0/0 | message-changed | `same` predicted: **not held** (observation 1) |

All other files are `same` base against head. Their counts equal the declared base values, with one exception:
`test-kit254-omics-correlation-skip.R` (observation 2).

**Observation 1: `test-gui-entry-metabolomics.R`, metab_gc lane, base side only.** This is an instrument artifact,
not KIT-261. The evidence is OBSERVED; the cause is INFERRED.
- 8 base failure records differ from the head's. Line numbers are masked. Both sides have 48 failures, and the other 40
  records are identical.
- Each of the 8 differs only in lines that name one artifact,
  `scripts/metabolomics/export_complete_dia-export-f897e88f6.json`, as present in actual and absent in expected.
- The lane's normaliser masks only hexadecimal runs of 10 or more characters
  (`tests/testthat/helper-gui-entry.R:60`). The expected value holds `…-export-<HEX>.json`.
- On the base side, this run's random export id had only nine characters, `f897e88f6`, so it escaped the mask.
- On the head side the id was masked, and no unmasked export id appears in the head's failures.
- Inferred: the export id is drawn at random. When it has fewer than 10 hexadecimal characters, the normaliser leaves it
  visible. That is a test-harness instability, independent of the code under test.
- KIT-261 changes neither the export nor its id. Its change is the run button's state, `R/mod_omic_norm.R`.
- The counts are equal, 57/48, on both sides.

**Observation 2: `test-kit254-omics-correlation-skip.R` reads 376/3 on both sides.** The declared base value was 379/0,
from train 15's declaration ("same; base train 13 declared head"). The verdict is `same`, so this is not KIT-261.
- The 3 failures are KIT-254's E7 lipid_ruv assertions, identical on both sides:
  - "a Skip after an Apply commits the run's result", once;
  - "a re-applied filter reads the run's result", twice: `kit254Runs`/`runs` against `kit254HiRuns`/`kit254LoRuns`.
- E7 runs at the constructed threshold 0.960024, which was chosen from the lipid-with-RUV fixture's correlations.
- Train 15's KIT-190 changes how RUV control features are chosen. Inferred: that moves those correlations on train 15's
  tree.
- Train 15's own run (45681065939) was still in progress when this was read. If it measured this file, it should see
  the same move on its head: 376/3, against its prediction of 379/0. That is train 15's finding, not this train's.

**R CMD check:** 0 ERRORs, 3 WARNINGs and 3 NOTEs on both sides. All 6 findings are identical, and none is only on
head. The portable-file-names NOTE lists 238 paths on each side.

**Verdict.** KIT-261's predictions held: its new file passed 105/0, nothing else moved base against head that its change
could reach, and check is identical. Two observations were not predicted, and neither is attributable to KIT-261:
- a base-side random-id masking artifact in `test-gui-entry-metabolomics.R`;
- a base-side count in `test-kit254-omics-correlation-skip.R`, carried from train 15.

Per the decision rule they are reported, not repaired. The coordinator rules whether the train is green.

**Landing.** The train lands only after train 15 lands with its code unchanged. Staging is 4d588ba3 (ledger rows through
VAL-313). At landing:
- the ledger keeps staging's numbered rows over this branch's placeholder forms;
- `VAL-<train16-1>` is numbered.

**Not checked:**
- the trees after the run's head (docs only);
- files outside the 44;
- train 15's own result;
- whether a 9-character export id recurs (no re-run: a valid result is not re-run).
