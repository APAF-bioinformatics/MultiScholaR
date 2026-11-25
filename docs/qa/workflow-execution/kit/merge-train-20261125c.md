# Merge train 3 of 2026-11-25: KIT-212, KIT-215, KIT-174 (a), KIT-225

Coordinator's instruction after train 2 landed (staging dd769671, KIT-200 = VAL-211): same
procedure as trains 1 and 2 (`merge-train-20261125.md`, `merge-train-20261125b.md`).

## Train

Base: `dd769671687aa8cdf5d7858ce3e6f804a3cd0fa6` (`origin/v0.9-staging`).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-212 | `origin/kit-212` 650d5adf | 742b76b8 (records 5f76dcf0: closure check C6) | `VAL-<train3-1>` | AUDIT-KIT-212-impl-20261125, closure check |
| 2 | KIT-215 | `origin/kit-215` c5c29c68 | 29d0ec04 | `VAL-<train3-2>` | AUDIT-KIT-215-impl-20261125, closure check; failing-first 45525963187 |
| 3 | KIT-174 (a) | `origin/kit-174` 65656b00 | 3d80d17e | `VAL-<train3-3>` | AUDIT-KIT-174a-impl-20261125 (I11) |
| 4 | KIT-225 | `origin/kit-225` 600dcc6e | 8bce32aa | `VAL-<train3-4>` | AUDIT-KIT-225-impl-20261125 (merge readiness) |

Head (bound by the run): `8bce32aa` (full hash in the request). This file and the records
commit after the run are under `docs/qa/` only.

Conflicts met: the ledger in all four merges (staging's rows through VAL-211 kept, each ticket's
row after them, written `VAL-<train3-n>` for the coordinator to number); `NEWS.md` in merges 2, 3
and 4, each the same adjacency at the top of "### Fixes" under the one unreleased heading,
resolved by keeping both entries (earlier train members first). No code conflict. `DESCRIPTION`
and `NAMESPACE` unchanged. KIT-174's `VAL-<merge>` replaced in its three places (ledger row,
UVR-021's status, `parked-20261113.md` row 20); each ticket's `outcome.md` placeholder likewise.

## Declaration (written before the run is sent)

- Purpose: the check on the merged tree that each member's closure requires, done once on the
  train.
- Scope: the union of the four merge-check scopes plus KIT-200's file, 108 files, base and
  head, R CMD check on both sides:
  - KIT-212: R4's 60 files (`kit-212/r4-compare.tsv`) and `test-integration-removed-factors.R`;
  - KIT-215: R2's 26 files (named in the closure check C5) and `test-kit212-da-current.R`.
    `test-kit220-da-threshold-consumers.R` is left out: it is on KIT-220's branch, which is not
    in this train;
  - KIT-174 (a): the files of its register's I11;
  - KIT-225: the 30 files of `kit225-r2-17643d77.json`;
  - `test-kit200-prot-da-no-value-proteins.R` (KIT-200, now in staging).
- Environment, the same on both sides: `extra_packages` sva, limpa, mixOmics;
  `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`; per-file limits of 120 minutes for the four GUI-entry
  files and `test-prot-norm-core.R` (kill limits, not thresholds); 10 shards per side.
- Ceiling: 2,730 job-minutes, the sum of the four members' merge-check ceilings (KIT-212 R2b
  1,080; KIT-215 R2b 540; KIT-174 R2-type check 34.0 core-hours = 510; KIT-225 R2 600).
- Predictions (each register; judged file by file):
  - KIT-212: `test-kit212-da-current.R` only on head, 309/0; `test-omics-summary-core.R`
    676/25 on both sides, `message-changed` acceptable only if label-only (E16, closure C2);
    `test-omic-summary-module.R` `same` (1675/28); `test-kit207-summary-head-guard.R` `same`.
  - KIT-215: `test-kit215-not-tested-lists.R` only on head, 112/0; its other R2 files `same`.
  - KIT-174 (a): `test-kit174-ruv-grouping.R` only on head, 77/0 (A1 45522632448); its I11
    files `same`.
  - KIT-225: `test-kit225-prot-norm-reset-source.R` only on head, 56/0/0 skipped;
    `test-prot-05b-norm-module-contracts.R` `count-changed` 439 to 440 passed, 0 failed;
    `test-prot-norm-core.R` 1548/306 on both sides; `test-kit129-run-path.R` 711/0 on both
    sides (trx flag set).
  - KIT-200's file `same`.
  - Every other file `same`. R CMD check: findings identical on both sides, or fewer on head.
- Decision: all predictions hold: TRAIN GREEN. Any other difference is attributed to a member
  where possible and reported with its files and counts; the train is not rebuilt here.
- Attempts: one run; one re-send only for an instrument defect with no result judged.

## Result

Run 45573711299 (request `merge-train-20261125c-8bce32aa`, run-tests commit cde483c6), head
8bce32aa against base dd769671, sent once; conclusion success; 810.5 job-minutes of the 2,730
ceiling (check 28.1). Tables in `merge-train-20261125c/` (`compare.tsv`, `tests.tsv`,
`tests-base.tsv`, `check-compare.txt`, `check-summary.txt`, `meta.txt`).

Observed, 108 files: 102 `same`, 6 not `same`:

| file | base p/f/e/s | head p/f/e/s | verdict | member |
|---|---|---|---|---|
| test-kit212-da-current.R | absent | 309/0/0/0 | only-on-head | KIT-212, predicted |
| test-omics-summary-core.R | 676/25/0/2 | 676/25/0/2 | message-changed | KIT-212, predicted (label-only, below) |
| test-kit215-not-tested-lists.R | absent | 112/0/0/0 | only-on-head | KIT-215, predicted |
| test-kit174-ruv-grouping.R | absent | 77/0/0/0 | only-on-head | KIT-174 (a), predicted |
| test-kit225-prot-norm-reset-source.R | absent | 56/0/0/0 | only-on-head | KIT-225, predicted |
| test-prot-05b-norm-module-contracts.R | 439/0/0/22 | 440/0/0/22 | count-changed | KIT-225, predicted |

Predicted `same` and observed `same`: `test-omic-summary-module.R` 1675/28,
`test-kit207-summary-head-guard.R` 292/0, `test-prot-norm-core.R` 1548/306,
`test-kit129-run-path.R` 711/0 (1 skipped), `test-kit200-prot-da-no-value-proteins.R` 538/0,
`test-workflow-omics-binding-steps.R` 58/1, the four GUI-entry files (proteomics 189/4,
lipidomics 35/0, metabolomics 105/0, transcriptomics 16/13).

`test-omics-summary-core.R`, label-only (KIT-212 E16): 25 failure records on each side. After
masking temporary paths and replacing the expectation label
`omicSummaryW1Frozen[[case]][[step]][[field]]` with `omicsSummaryCoreExpected(case, step, field)`,
the two lists are identical. Two raw records differ, as in KIT-212's R4.

R CMD check: 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides, all 6 findings identical, none only on
head.

Verdicts: KIT-212, KIT-215, KIT-174 (a) and KIT-225 each meet their predictions; no other file
differs. TRAIN GREEN.

When the result was read, staging had moved to 3eeaf100: ledger rows VAL-212 to VAL-225 only
(`git diff dd769671 3eeaf100` touches only `docs/methods-decision-ledger.md`). The train lands by
a merge; the rows stay `VAL-<train3-1>` to `VAL-<train3-4>` for the coordinator to number.

Not checked: the trees after the run's head (docs only); files outside the 108 (no member's full
suite on the train; KIT-225's R3 is left to the coordinator, as its `outcome.md` states); the merge
with staging 3eeaf100.
