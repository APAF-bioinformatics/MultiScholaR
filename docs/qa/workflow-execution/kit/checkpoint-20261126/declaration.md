# Full-suite checkpoint on staging 828dd105 (2026-11-26)

This is the coordinator's checkpoint after trains 8 and 9 landed. It records the state of the whole test suite and of
R CMD check on the landed tree. Each test file runs in its own process, as in the programme's full-suite runs.

The branch `checkpoint-20261126` is docs only. The coordinator merges it.

## Binding

- **Head:** `828dd1051d2acab075d6d04f887b7c9e9ccaa48e` (`origin/v0.9-staging` after trains 8 and 9; ledger rows through
  VAL-272).
- **Reference:** `dd769671687aa8cdf5d7858ce3e6f804a3cd0fa6`. This is the staging tree of the most recent full-suite
  results on record:
  - KIT-201 R3, run 45592673648: base dd769671, head c25e3e1d, 1,242.4 job-minutes, 20 shards per side, `sva`, no env;
  - KIT-220 R3: base dd769671, trx switch set.

  It is also train 3's base, so trains 3 to 9 (`merge-train-20261125c` to `-g`, `merge-train-20261126a` and `-b`)
  account for every change from it to the head.
- **How the reference is used:** it runs again in this run as `base_commit`, under the same environment as the head.
  The file-by-file comparison (`compare.tsv`, `compare-failures.txt`) therefore has one runner image, one library and
  one environment on both sides. The recorded R3 base tables are a cross-check: base side against run 45592673648's
  base, for files outside the environment's reach.
- **Request form:** KIT-201 R3's (`kit201-r3-c25e3e1d.json`), with these settings:
  - `tests: "all"`, 20 shards per side;
  - the same file timeouts (120 minutes for the four GUI-entry files and `test-prot-norm-core.R`);
  - extra packages `sva`, `limpa` and `mixOmics`, so that KIT-225's E3 runs as in train 3;
  - env `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1` on both sides;
  - `check: true`.
- **Scripts job:** KIT-240's `docs/qa/workflow-execution/kit/kit-240/default-identity.R`, as its `outcome.md` release
  checklist asks. It runs against base `67aed57f`, staging before train 9 (without KIT-240), with a 60-minute limit, as
  in train 9. It runs on the head install in the same run, so no second run is needed.

## Scope

Every `tests/testthat/test-*.R` file on the head: **505 files**. The reference tree has 481. 24 files were added since,
and none was removed (`git diff --name-status dd769671 828dd105 -- tests/testthat`). 36 existing test files were
modified, all by the merged tickets.

## Tickets on this tree, and the post-merge full suites it stands in for

| Ticket | On 828dd105? | Merged by | Stands in for |
|---|---|---|---|
| KIT-201 | yes | train 8 (VAL-266, VAL-267) | R3b |
| KIT-243 | yes | train 8 (VAL-265) | R3 |
| KIT-218 | yes | train 4 (VAL-233) | A4 |
| KIT-222 | yes | train 4 (VAL-232) | A5 |
| KIT-240 | yes | train 9 (VAL-269) | A7, plus the identity script |
| KIT-245 | yes | train 9 (VAL-270, VAL-271) | its full suite |
| KIT-256 | yes | train 9 (VAL-272) | its full suite |
| KIT-204 | **no** | — | not served |
| KIT-190 | **no** (ledger decisions only) | — | not served |
| KIT-236 | **no** (branch `kit-236`) | — | not served |

## Ceiling

1,553.0 job-minutes. This is the reference run's 1,242.4 plus 25 %.

## Prediction (written before the run)

`predictions.tsv`, beside this file, lists 51 files: every file a merged ticket declared as differing from its train's
base in trains 3 to 9. For each it gives:
- the declaring ticket and its VAL row;
- the base value (first train's base);
- the expected head value (last train's head);
- the expected verdict.

The 51 are of four kinds:
- **24 new files:** only on head, with no failure.
- **Count changes declared by tickets:**
  - KIT-225 (VAL-229);
  - KIT-234 (VAL-231): integration-removed-factors, go-propagation, enrichment-significance, integration-workbook;
  - KIT-239 (VAL-255) and KIT-240 (VAL-269): integration gui, record, results, model;
  - KIT-224 (VAL-206): integration-input-transcriptomics, kit127, lipid-12c, metab-03a.
- **Stored references that already fail on runners and wait for the end-of-project re-record**, moving only by the
  declared lines:
  - VAL-181 and VAL-189: the declarations and binding-steps snapshots;
  - VAL-259: KIT-214's lines;
  - VAL-206 and VAL-213: KIT-224's lines and digests;
  - VAL-226: KIT-212's E16 label;
  - VAL-230: KIT-220's W1 values;
  - VAL-236: the gui-entry-transcriptomics export leaf;
  - VAL-246: KIT-244's re-pin of `test-prot-norm-core.R:45`;
  - VAL-266 and VAL-268: KIT-201's moved references (D5) and KIT-243's saved-da `files` values.
- **One declared fix:** `test-prot-ruv-factors-in-model.R` (KIT-201).

Every other file of the 505 is predicted `same`.

Known flaky leaves, named in advance:
- KIT-201 I13's three `refusals` leaves in `test-omic-enrich-module.R`. If one fails on base, that file reads 22 to
  59, not 19 to 59.
- Plan section 11's `gprofiler_plot`.

R CMD check is predicted to give findings identical on both sides, matching the trains: 0 ERRORs, 3 WARNINGs and 3
NOTEs, with 238 paths in the portable-file-names NOTE (train 8). On the head this NOTE lists no path new since dd769671.

The KIT-240 identity script is predicted to give its verdict file all TRUE, with exit 0 (train 9's result).

## Decision rule

- A file that differs from the reference, or fails, and is **not** in `predictions.tsv` is a **finding**. It is
  recorded with its test names and messages, and is not fixed or re-run here.
- A predicted file whose difference goes beyond its declared kind is also a finding.
- An unfavourable valid result is never re-run.
- An instrument defect with no result judged may be re-sent under VAL-54, at most three times, and is recorded.
- **CLEAN** means every difference is explained by `predictions.tsv`.

## Result

**Verdict: CLEAN.** Every difference is explained by `predictions.tsv`, and there are no findings.

Run 45644467847 (request `checkpoint-20261126-828dd105`, run-tests commit cabf3d8f) was sent once and finished with
conclusion success. Its 45 jobs used 1,353.9 job-minutes of the 1,553.0 ceiling: check 30.8, scripts 11.1. Tables
beside this file: `compare.tsv`, `tests.tsv`, `tests-base.tsv`, `check-compare.txt`, `check-summary.txt`, `meta.txt`;
KIT-240's identity outputs are under `id240/`.

### Verdicts, 505 files

| Verdict | Files | Predicted |
|---|---|---|
| same | 454 | 454 |
| only-on-head | 24 | 24: the new files, every one passing (failed 0) |
| count-changed | 13 | 13 |
| message-changed | 9 | 9 |
| new-failure | 4 | 4: KIT-201's moved references |
| fixed | 1 | 1: `test-prot-ruv-factors-in-model.R` |

All 51 predicted files differ, and no other file does. Each predicted file's base and head counts equal the
predicted values exactly: the first train's base and the last train's head (passed, failed, skipped), as listed in
`predictions.tsv`.

Selected rows:

| File | Base | Head |
|---|---|---|
| test-workflow-saved-da.R | 780/171 | 776/175 |
| test-gui-entry-proteomics.R | 189/4 | 168/25 |
| test-omic-da-module.R | 1316/53 | 1289/80 |
| test-omic-enrich-module.R | 1051/19 | 1011/59 |
| test-omics-enrich-core.R | 232/9/1 | 221/20/1 |
| test-integration-gui.R | 397 | 405 |
| test-integration-results.R | 169 | 258 |
| test-integration-record.R | 133 | 167 |
| test-integration-model.R | 263 | 265 |
| test-workflow-trx-declarations.R | 309/33 | 313/31 |
| test-module-ci-trx-summary.R | 86/2 | 87/1 |
| test-prot-norm-core.R | 1548/306 | 1548/306 |

In `test-omic-enrich-module.R` the I13 leaf passed on base, so the move is 19 to 59, not 22.

### Message level

Every failing file on the head was compared with the latest train run that held it, by failure records with temporary
paths masked.
- The head's failure records are identical to train 8's head (45626605316, tree aceb7190) in every failing file that
  train ran. That covers the 13 message-changed and new-failure files except the declarations snapshot,
  `test-kit127-omics-ruv-factors.R`, and ten pre-existing failing `same` files.
- `test-workflow-omics-declarations.R` is identical to train 6's head (45608483144).
- `test-trx-publication.R` (a `same` file) is identical to train 6's head.
- So no message moved beyond what the trains observed and declared.
- Of the 20 files that fail on both sides with `same`, nine have no train record to compare. They are `same` against
  the reference, so their failures are the reference's own, recorded before the trains:
  `test-general-filemgmt-load-dependencies-shared.R`, `test-module-ci-trx-design.R`, `-norm.R`, `-qc.R`,
  `test-optional-packages.R`, `test-trx-design.R`, `test-workflow-omics-owner.R`, `test-workflow-perf101-diagnostic.R`
  and `test-workflow-sleep-guard.R`.

### R CMD check

- 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; all 6 findings identical, none only on head. This matches the trains.
- The portable-file-names NOTE lists 238 paths on each side, the same as train 8. No path new since dd769671 is in it.

### KIT-240 identity script

`default-identity.R` 67aed57f against 828dd105 exited 0 after 27 s. The verdict file holds 23 checks, all TRUE ("all
checks pass: TRUE"):
- (i): head with the option unset equals base with it TRUE, after the declared substitutions (counts 1, 2 and 21 in
  base and 0 in head; the Enrichment mark inserted 8 times on head);
- (ii): head FALSE equals base unset; the module is the unavailable panel;
- (iii): head renders the same in two processes.

This is the release-checklist item from KIT-240's `outcome.md` that this run could carry.

### Tickets served

On this tree, the checkpoint stands in for these post-merge full suites:
- KIT-201 R3b;
- KIT-243 R3;
- KIT-218 A4;
- KIT-222 A5;
- KIT-240 A7;
- KIT-245;
- KIT-256.

KIT-204, KIT-190 and KIT-236 are not merged and are not served.

### Findings

None.

### Not checked

- The tree after this record (docs only).
- The other KIT-240 release-checklist items (the systemd observation, the enrichment defaults): those are the owner's.
- The e2e browser lanes, which skip in run-tests.
