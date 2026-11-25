# KIT-246: outcome (first part: tests-only commit and R1)

Ticket `tickets/kit/KIT-246.md` r1.2 (step 0: AUDIT-KIT-246-plan-20261125 Re-check 1, R1-01 to R1-05, applied
with the register's texts). Branch `kit-246` from `v0.9-staging` at 3eeaf100.

**State: R1 done, as expected. Stopped by instruction before the fix commit.** The fix commit and R2 wait for
KIT-239's merge (E5 reads its `help_scaling`) and the owner's answer to item A (the wording; recommended option (a),
ticket section 14).

## What was run

| Step | Commit | Run | What it decides |
|---|---|---|---|
| R1 | 2d81e890 (tests only: one `test_that` block appended to `tests/testthat/test-integration-record.R`; nothing under `R/`) | 45584885011 (request `kit246-r1-2d81e890-2`), base none, 1 shard, file timeout 45 min | stop rule 2 (every row of the failing-first list fails at its base value) |

Declared before the run (ticket sections 10 and 11): purpose failing-first; binds 2d81e890; ceiling 45 job-minutes
(3.00 core-hours), expected 24; E1 fails by assertion (the paragraph contains `OLD`, "weak factors" and "no practical
difference"), G1 passes, every existing test passes.

An earlier request, `kit246-r1-2d81e890` (run 45584833302), was pushed with the request file still blank (the
edit had not been written before the commit). The run's `plan` job refused it after 7 s; no test job started and no
result was produced. It was resent under a new name with every field stated and nothing else changed (`runs.tsv`).

## Observed (R1, run 45584885011; `r1-tests.tsv`, `r1-failures.txt`, `r1-tests-only.log`, `r1-meta.txt`)

- `test-integration-record.R`: 139 passed, 6 failed, 0 errors, 0 skipped, 0 warnings, 667.4 s. Test job 23.2 min
  (library set-up 717 s), plan and results jobs 0.8 min.
- All six failures are in the new block, at lines 434, 435 and 436, each twice (the paragraph as built and as read
  back from JSON), each `actual: TRUE`, `expected: FALSE`:
  - 434: `grepl("weak factors", methods, fixed = TRUE)`;
  - 435: `grepl("no practical difference", methods, fixed = TRUE)`;
  - 436: `grepl(old, methods, fixed = TRUE)` (`old` = the current sentence written out).
- No failure or error at line 433 (the paragraph contains the text table's `methods_scaling_feature_zscore`) or in
  G1 (lines 438-447: one parameter named `scaling`; with its value replaced by `encodeWorkflowRecordValue("none")`
  the paragraph contains `methods_scaling_none` and "without per-feature scaling." and not the z-score sentence).
- No failure in any of the file's other blocks.

## Inferred

- E1 fails at exactly its stated base value: today every record with the default scaling prints the old sentence,
  with "weak factors" and "no practical difference", from the built record and from its JSON alike. Stop rule 2 is
  not triggered.
- G1 is a valid guard on base: the encoded replacement (R1-01) decodes, and the scaling clause follows the recorded
  value.
- The passes not separately counted per block: 139 passed = the existing blocks plus E1's two "contains" checks and
  G1's four; no failure falls outside E1, so every existing test passes (A1).

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met | R1: E1 fails by assertion at its base value (6 failures, lines 434-436, both renderings); G1 passes; no other failure in the file |
| A2 | open | fix commit with the owner's words, then R2 (waits for KIT-239's merge and item A) |
| A3 | open | next checkpoint's full suite |
| A4 | open | NEWS, INT and VAL ledger rows, "For the change report" (with the fix) |

## Budget

R1: 24.0 job-minutes = 1.60 core-hours (expected 1.60; ceiling 3.00). The refused request: plan job only, about
0.1 job-minutes. Ticket total so far 1.61 of 53.33 core-hours.

## Not checked

- Per-block pass counts (the runner reports per file); the attribution of the 139 passes rests on the failures'
  line numbers.
- The fix, E2, E5, R2, `R CMD check` (not yet due).
- KIT-239's L4/L5 in this file (R1 branched from 3eeaf100, before KIT-239's merge; R1-05).
- The owner's words (item A).
