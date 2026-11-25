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

Run 45586240983 (request `merge-train-20261125d-b7ee701f`, run-tests commit ca111465), head b7ee701f against the
stacked base 5174a798, sent once; conclusion success (every job success); 811.0 job-minutes of the 3,105 ceiling
(tests 768.5, check 29.7, scripts 11.9, plan and results 0.8). Tables in `merge-train-20261125d/` (`compare.tsv`,
`tests.tsv`, `tests-base.tsv`, `check-compare.txt`, `check-summary.txt`, `meta.txt`, `w6-console.log`).

Observed, 165 files: 154 `same`, 11 not `same`:

| file | base p/f/e/s | head p/f/e/s | verdict | ticket, prediction |
|---|---|---|---|---|
| test-kit220-da-threshold-consumers.R | absent | 79/0/0/0 | only-on-head | KIT-220, as declared |
| test-workflow-saved-da.R | 780/171/0/0 | 780/171/0/0 | message-changed | KIT-220, as declared (below) |
| test-kit227-report-ruv-sentences.R | absent | 166/0/0/0 | only-on-head | KIT-227, as declared |
| test-kit235-lipid-metab-report-parser.R | absent | 50/0/0/0 | only-on-head | KIT-235, as declared |
| test-enrichment-significance.R | 191/0/0/0 | 181/0/0/0 | count-changed | KIT-234, as declared |
| test-go-propagation.R | 331/0/0/1 | 321/0/0/1 | count-changed | KIT-234, as declared |
| test-integration-removed-factors.R | 526/0/0/3 | 520/0/0/3 | count-changed | KIT-234, as declared |
| test-integration-workbook.R | 53/0/0/0 | 43/0/0/0 | count-changed | KIT-234, as declared |
| test-kit234-workbooks-withdrawn.R | absent | 29/0/0/0 | only-on-head | KIT-234, as declared |
| test-kit222-replicate-keys.R | absent | 456/0/0/0 | only-on-head | KIT-222, as declared |
| test-kit218-da-wrapper-tail.R | absent | 80/0/0/2 | only-on-head | KIT-218, as declared (G1, M2 skipped) |

- `test-workflow-saved-da.R`: of the 171 failure messages on each side, 134 are identical and 37 differ (37 only on
  base, 37 only on head, matched by test and leaf). The 37 are exactly the 19 W1 leaves of `savedDaW1Frozen`
  (dia_qvalue_only products, dependencies, manifest, head_manifest, display, files; dia_refused products,
  dependencies, head_manifest, display, files; `files` of dia_one_contrast, dia_two_contrasts, dia_repeat,
  dia_reopened_rerun, dia_non_head, dia_changed_arguments, dia_retry, limpa_one_contrast) and their 18 W2 twins
  (all but dia_non_head).
- `test-module-ci-trx-summary.R` `same` (86/2/0/0 on both sides); head log: MCI-032.1 counted, skip 0, no
  "KIT-040" in the log (KIT-227 condition (iv)).
- Among the `same` files: `test-integration-workbook-mofa.R` 3 skipped on both; `test-omic-gui-config.R` 1911/0;
  `test-kit212-da-current.R` 309/0; `test-kit215-not-tested-lists.R` 112/0; `test-gui-entry-proteomics.R` 189/4;
  `test-gui-entry-transcriptomics.R` 16/13; `test-gui-entry-lipidomics.R` 35/0; `test-gui-entry-metabolomics.R`
  105/0; `test-omic-summary-module.R` 1675/28; `test-omics-summary-core.R` 676/25 (2 skipped);
  `test-module-ci-trx-design.R` 66/0/1 error; `test-prot-ruv-factors-in-model.R` 306/1;
  `test-integration-input-transcriptomics.R` 78/0 on both sides (KIT-222's register named 3 pre-existing failures
  here, on its own run's base; on this train both sides pass: still `same`).
- W6 (KIT-234): exit 0; `system.file("workbooks")` = `""`; the built tarball holds 8,120 paths, 0 matching
  `^MultiScholaR/(inst/workbooks|archive)(/|$)`.
- R CMD check: 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; 5 findings identical; ONE finding differs: the NOTE
  "checking for portable file names" lists 246 paths on head against 238 on base. The 8 added paths are all
  KIT-235's committed G0 records over 100 bytes (`docs/qa/workflow-execution/kit/kit-235/g0-45573071526/g0b/`
  `{lip,met}-c17-src-analysis_record.json`, `{lip,met}-c17-src-study_parameters.txt`,
  `{lip,met}-c18-sr-study_parameters.txt`, `{lip,met}-c19-sr-analysis_record.json`, 102 or 103 bytes with the
  `MultiScholaR/` prefix), committed in 96e3c813. KIT-235's own R2 did not show it because both its base and head
  already held those files. No other path changed; no code, Rd or namespace finding differs.

Verdicts against the predictions:

- KIT-220: holds exactly (new file 79/0; saved-da message-changed at exactly the 37; every other file `same`).
- KIT-227: holds exactly (new file 166/0; its files `same`; condition (iv) met).
- KIT-235: its tests hold exactly (new file 50/0; its 14 files `same`). NOT as declared on the check: its records add
  8 paths over 100 bytes to the portable-file-names NOTE (a records defect: the implementer brief's 100-byte rule).
- KIT-234: holds exactly (four count changes, new file 29/0, W6 exit 0, saved-da as KIT-220).
- KIT-222: holds exactly (new file 456/0; its 91 other files `same`). C1 not re-run: blob `651fc1b6` at the head.
- KIT-218: holds exactly (new file 80/0, 2 skipped).
- Check: not identical, one NOTE's path list grown by KIT-235's 8 record paths (above).

Train verdict: every test prediction holds; the check prediction does not hold, by one NOTE text attributable to
KIT-235's record paths only. Not TRAIN GREEN as declared; the coordinator decides (the record files can be renamed
or moved under 100 bytes, a docs-only change that no test reads; the run is not repeated here).

Train 3 landed on staging (16112917; KIT-212 = VAL-226, KIT-215 = VAL-227, KIT-174 (a) = VAL-228, KIT-225 =
VAL-229) while the run was in progress. `git diff --stat 5174a798 origin/v0.9-staging`: 13 files, all
`docs/methods-decision-ledger.md` or under `docs/qa/workflow-execution/kit/` (the numbered rows, four outcome files,
`parked-20261113.md`, train 3's record files); no other file differs. The train's `VAL-<train3-n>` and
`VAL-<train4-k>` placeholders are left for the coordinator at landing.

Not checked: the trees after the run's head (docs only); any file outside the 165; MOFA lane (not required for
KIT-234 alone); the landing merge with staging 16112917 (ledger rows only, per the coordinator).

## Correction (coordinator's ruling, the train's one correction pass)

Cause of the check difference in run 45586240983: 8 of KIT-235's G0 record files (committed in 96e3c813) are over
100 bytes with the `MultiScholaR/` prefix. KIT-235's own R2 check (45577137662) could not show this: its base
f53db55a already held the files, so the NOTE was identical on both sides there. For KIT-235's records: note added in
`kit-235/outcome.md`.

Correction, commit e84c5571 (docs only): the 8 files renamed with `git mv`, contents unchanged (100 % renames), by
one rule for these 8 names only: `analysis_record` becomes `record`, `study_parameters` becomes `params`. References
updated: `kit-235/g0-45573071526/g0b/files.tsv` (8 rows, column 4) and `kit-235/outcome.md:59`. Left verbatim, as run
outputs: `kit-235/g0-45573071526/console.log`, `g0c/checks.txt`, `kit-235/r2-45577137662/check-summary.txt`, and this
train's `merge-train-20261125d/` tables. Script check on the new head: 255 paths absent from the base 5174a798, none
over 100 bytes with the prefix; the longest is
`docs/qa/workflow-execution/kit/kit-235/g0-45573071526/g0b/met-c04-src-design_matrix.tab`, 100 bytes.
`git diff --name-only b7ee701f e84c5571`: 18 paths, all under `docs/qa/`, so the test verdicts of run 45586240983
stand for this head.

### Declaration of the second run (written before it is sent)

- Purpose: repeat R CMD check only, after the renames.
- Base 5174a798627cbe425cfc620b273625ff3dcfe11f; head e84c5571c0752af9519975b59ab994e8d5562128.
- Tests: none (`"tests": []`, which the request format allows); `check: true`; no scripts; env and BLAS as the first
  run.
- Prediction: check findings identical on both sides, including the "portable file names" NOTE with 238 paths on
  both sides; 0 ERRORs, 3 WARNINGs, 3 NOTEs on each.
- Ceiling: 60 job-minutes (the first run's check took 29.7 for both sides), charged to the train's ceiling (811.0 of
  3,105 used before it).
- Decision: holds: the check prediction is met and the train is green. Does not hold: stop and report; no further run.
