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

(filled after the run)
