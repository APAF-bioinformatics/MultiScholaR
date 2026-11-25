# Merge train 4 of 2026-11-25: KIT-220, KIT-227, KIT-235, KIT-234, KIT-222, KIT-218

Coordinator's instruction of 2026-11-25: same procedure as trains 1 to 3 (`merge-train-20261125.md`,
`merge-train-20261125b.md`, `merge-train-20261125c.md`). Train 4 is stacked on train 3: it lands only after train 3
lands. If train 3's check fails or train 3 is rebuilt, this train is rebuilt.

## Train

Stacked base: `5174a798627cbe425cfc620b273625ff3dcfe11f` = train 3's head `origin/merge-train-20261125c`
c719b6c6 merged with `origin/v0.9-staging` 3eeaf100 (staging's ledger rows VAL-212 to VAL-225 kept in order,
train 3's placeholder rows `VAL-<train3-1>` to `VAL-<train3-4>` after them, unchanged; no other file differs).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-220 | `origin/kit-220` a662f3e1 | 50239c31 | `VAL-<train4-1>` | AUDIT-KIT-220-impl-20261125, closure check C.7 |
| 2 | KIT-227 | `origin/kit-227` 5e7efe5f | eb0db68c | none (owner decisions are VAL-202, on staging) | AUDIT-KIT-227-impl-20261125, section 6, I11 |
| 3 | KIT-235 | `origin/kit-235` 36683696 | 540191a9 | none | AUDIT-KIT-235-impl-20261125, section 4 |
| 4 | KIT-234 | `origin/kit-234` 84f8ca4a | c8b6a458 | `VAL-<train4-2>` | AUDIT-KIT-234-impl-20261125, I07 |
| 5 | KIT-222 | `origin/kit-222` 75501687 | 11d55342 | `VAL-<train4-3>` | AUDIT-KIT-222-impl-20261125, closure check (4) |
| 6 | KIT-218 | `origin/kit-218` dae1220f | b7ee701f | `VAL-<train4-4>` | AUDIT-KIT-218-impl-20261125, closure check "Merge" |

Head (bound by the run): `b7ee701ffcc104e3636e5cf444b10def6f98a0da`. This file and the records commit after the run
are under `docs/qa/` only.

Conflicts met and how each was resolved:

- Staging into train 3: the ledger only (above).
- KIT-220: the ledger (its row after train 3's placeholders, written `VAL-<train4-1>`); `NEWS.md` (KIT-220's four
  entries after train 3's, both kept); and the one authorised code conflict, `R/utils_workflow_trx_declarations.R` in
  `computeWorkflowTrxDaPublicationAction` (about lines 904 to 918). Kept: KIT-212's `files <- trxDaPublicationNames(name)`
  with `files$table`, `files$volcano` and `files$bar`, TOGETHER WITH KIT-220's `fdr_thresh = qvalue` on all three calls
  (`trxAnnotateDe`, `trxVolcano`, `trxSignificanceBar`), exactly the text of AUDIT-KIT-220 C.7. KIT-220's
  `qvalue <- manifest$parameters$qvalue` line, the `trxSummaryText(..., fdr_thresh = qvalue)` call and
  `trxVolcanoWidget` merged without conflict. `VAL-<merge>` replaced in the ledger row and in `kit-220/outcome.md`
  (its staging-merge line).
- KIT-227: `NEWS.md` only (its entry after KIT-220's). The seven report templates merged automatically. The branch
  carries no `VAL-<merge>`.
- KIT-235: no conflict (NEWS merged automatically, its entry after KIT-227's).
- KIT-234: the ledger only (`VAL-<train4-2>`; also in `kit-234/outcome.md`, its commit table).
- KIT-222: the ledger only (`VAL-<train4-3>`; also in `kit-222/outcome.md`, its commit table).
  `R/func_omics_design_shared.R` at the train's head has blob `651fc1b6` (observed), so C1 is not re-run.
- KIT-218: the ledger (`VAL-<train4-4>`) and `NEWS.md` (both entries kept, KIT-218's after the train's earlier entries). Edits owed
  at merge applied: `KIT-<new>` becomes KIT-249 in its three places (twice in `kit-218/outcome.md`, once in the ledger
  row); `VAL-<merge>` becomes `VAL-<train4-4>` in `kit-218/outcome.md` (three places); `runs.tsv` row A1 "(W1, W3, W4
  prints)" becomes "(printed in W1, W2 and W3)".
- No test file conflicted. No other code conflict. `DESCRIPTION` (0.9.0), `NAMESPACE` unchanged by the train's merges;
  every NEWS entry is under the one unreleased heading, "### Fixes". No NEWS entry for KIT-222 (owner's decision).

## Declaration (written before the run is sent)

- Purpose: the check on the merged tree that each member's closure requires, done once on the train.
- Base: the stacked base 5174a798 (train 3's head with staging 3eeaf100). Head: b7ee701f.
- Scope: the union of the six members' merge-check lists, 165 files, base and head, R CMD check on both sides:
  - KIT-220: the 66 files of `kit-220/r2b-compare.tsv` (the new file included), `test-gui-entry-proteomics.R`,
    `test-gui-entry-transcriptomics.R`, `test-omic-gui-config.R`, `test-kit215-not-tested-lists.R`,
    `test-kit212-da-current.R` (70 distinct);
  - KIT-227: the 27 files of `kit-227/r2-45508506125/compare.tsv` (the new file included),
    `test-kit220-da-threshold-consumers.R`, `test-module-ci-trx-summary.R` (28 distinct);
  - KIT-235: the 15 files of `kit-235/r2-45577137662/compare.tsv` (the new file included);
  - KIT-234: the 14 files of R3's request (`kit234-r3-merged-51d53171`, run-tests c21faa52; the new file included);
  - KIT-222: R2's 72 files (`kit-222/r2/compare.tsv`) and the 20 added in run 45574235919
    (`kit-222/r3a/compare.tsv`), 92 distinct;
  - KIT-218: the 11 files of `kit-218/r2/compare.tsv` (the new file included).
  Partner files named in the registers whose tickets are not in this train (`test-kit201-treat-moderation.R`,
  `test-kit204-prot-da-not-tested.R`, `test-kit214-sample-columns.R`, `test-kit224-ruv-feature-levels.R`,
  `test-kit217-intensity-filter-missing.R`) are not on the train's head and are left out.
- Scripts on the head: `docs/qa/workflow-execution/kit/kit-234/w6.R` (KIT-234 I07 (3); timeout 30 minutes, as R3).
  No other member's register asks for a script at the train (KIT-222's C1 is not re-run, blob above).
- Environment, the same on both sides: `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`, as the earlier trains. No
  `extra_packages`: no member's register asks for limpa, mixOmics or sva, and every member's own merge-check evidence
  ran without them. Reference BLAS, ubuntu-22.04. 12 shards per side (files ordered by recorded duration, dealt so the
  shards balance). Per-file kill limits (not thresholds): 120 minutes for the four GUI-entry files, 90 for
  `test-workflow-saved-da.R`, 45 otherwise.
- Ceiling: 3,105 job-minutes, the sum of the members' merge-check ceilings (KIT-220 R2 540; KIT-227 R2 540;
  KIT-235 R2 540; KIT-234 R3 345; KIT-222 R2 840; KIT-218 R2 300). Expected about 1,100.
- Predictions (each from its register; judged file by file, base against head):
  - KIT-220: `test-kit220-da-threshold-consumers.R` only on head, 79/0. `test-workflow-saved-da.R`
    `message-changed`, 780/171 on both sides, the messages differing at exactly the 37 places of C.5 (the 19 W1
    leaves of `savedDaW1Frozen` and their 18 W2 twins), nothing else. Every other file of its 70 `same`, including
    the re-frozen files (GUI-entry captures, `test-omic-da-module.R`, `test-omic-enrich-module.R`,
    `test-omics-da-core.R`), `test-omic-gui-config.R` and the partner files `test-kit215-not-tested-lists.R` and
    `test-kit212-da-current.R`.
  - KIT-227: `test-kit227-report-ruv-sentences.R` only on head, 166/0. Its other 26 files `same`.
    `test-module-ci-trx-summary.R` `same`, and on the head's log MCI-032.1 counted, 0 skips, the KIT-040 skip
    message absent (condition (iv)).
  - KIT-235: `test-kit235-lipid-metab-report-parser.R` only on head, 50/0 (in its R2 the base carried the test file;
    here the base lacks it). Its other 14 files `same` (`test-omic-summary-module.R` 1675/28 and
    `test-omics-summary-core.R` 676/25 with the same failure lists).
  - KIT-234: `test-enrichment-significance.R` 191 to 181 passed, `test-go-propagation.R` 331 to 321 (1 skipped on
    both), `test-integration-removed-factors.R` 526 to 520 (3 skipped on both), `test-integration-workbook.R` 53 to
    43, each `count-changed` with 0 failed; `test-integration-workbook-mofa.R` 3 skipped on both sides;
    `test-kit234-workbooks-withdrawn.R` only on head, 29/0, 0 skipped; `test-workflow-saved-da.R` as KIT-220 above;
    the other files `same`. W6: exit 0, `""` and 0 matching paths.
  - KIT-222: `test-kit222-replicate-keys.R` only on head, 456/0. Its other 91 files `same`, the pre-existing failures
    identical on both sides.
  - KIT-218: `test-kit218-da-wrapper-tail.R` only on head, 80/0 with G1 and M2 skipped (2). Its other 10 files
    `same` (`test-prot-ruv-factors-in-model.R` its one known failure identical).
  - Every other file `same`. R CMD check: findings identical on both sides, or fewer on head.
- Decision: all predictions hold: TRAIN GREEN (it lands after train 3, and only as built on train 3's c719b6c6).
  Any other difference is reported with its files and counts and attributed to a member only where the evidence
  shows it; the train is not rebuilt here.
- Attempts: one run; one re-send only for an instrument defect with no result judged.

## Result

(filled after the run)
