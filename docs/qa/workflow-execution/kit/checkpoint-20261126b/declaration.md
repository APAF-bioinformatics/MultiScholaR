# Full-suite checkpoint on staging 15d86bd4 (2026-11-26), second of the day

Written before the run is sent. The branch `checkpoint-20261126b` is docs only. The coordinator merges it.

## Purpose

This run stands in for the full-suite acceptance rows of the tickets merged since the previous checkpoint
(`checkpoint-20261126`, run 45644467847 on 828dd105): KIT-259, KIT-257, KIT-249, KIT-236, KIT-254, KIT-248, KIT-190,
KIT-204 (its R3b), KIT-217, KIT-261, KIT-246 (its A3) and KIT-263, that is merge trains 10 to 18
(`merge-train-20261126c` to `-k`). It records the whole suite and R CMD check on the landed tree. Each test file runs
in its own process.

## Binding

- **Head:** `15d86bd4e3bb2438d7f198b18742439176b49079` (`origin/v0.9-staging`, fetched; it had not moved since the
  task was written, so no newer head applies).
- **Previous checkpoint's head, used as the base:** `828dd1051d2acab075d6d04f887b7c9e9ccaa48e`. It is the tree of
  run 45644467847's head side. It runs again in this run as `base_commit`, under the same environment, library and
  runner image as the head. The file-by-file comparison is therefore one runner image on both sides, and it is read a
  second way: the base side of this run against run 45644467847's head table
  (`checkpoint-20261126/tests.tsv`, 505 files) as a cross-check for files outside the reach of the environment.
- **Request form:** the previous checkpoint's, field for field (`checkpoint-20261126-828dd105.json`):
  - `tests: "all"`, 20 shards; file timeouts 45 minutes by default and 120 for `test-gui-entry-proteomics.R`,
    `-lipidomics.R`, `-metabolomics.R`, `-transcriptomics.R` and `test-prot-norm-core.R`;
  - extra packages `sva`, `limpa`, `mixOmics`;
  - env `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1` on both sides;
  - `check: true`; reference BLAS; ubuntu-22.04.
- **Scripts job: none.** The previous checkpoint's only script was KIT-240's `kit-240/default-identity.R`, requested
  as KIT-240's release-checklist item. That script is being retired by KIT-267 and is not requested. The previous
  checkpoint's rule requires no other identity script, so none is requested; kit-138's identity logic is not
  requested either.

## Scope

Every `tests/testthat/test-*.R` file on the head: **516 files**. The previous checkpoint ran 505. 11 files were added
since (`git diff --name-status 828dd105 15d86bd4 -- tests/testthat`), and none was removed.

## Ceiling

1,692.4 job-minutes = 112.8 core-hours (job-minutes x 4 / 60). This is the previous checkpoint's 1,353.9 job-minutes
plus 25 %. The previous figure included 11.1 job-minutes of the scripts job, which this run does not send.

## Predicted differences (written before the run)

`predictions.tsv`, beside this file, lists 27 files: every file the trains 10 to 18 declared or observed as moving
against its base, with the base being the previous checkpoint's head. Each row names its source train and section.
The head value is the latest observed value in the train chain (train 15's run 45681065939 for the files it held;
the later trains 16, 17 and 18 moved none of them except `test-integration-record.R` by train 17). The 27 are:

- 11 new files, only on head, every one passing: the KIT-204, KIT-217, KIT-190 (`test-ruv-control-anova-grouping.R`,
  139/0 with sva), KIT-236, KIT-249, KIT-248, KIT-254 (379/0, the edited file), KIT-257, KIT-259, KIT-261 and KIT-263
  files.
- 7 `count-changed`: `test-kit227-report-ruv-sentences.R` 166 to 162; `test-kit235-lipid-metab-report-parser.R`
  50 to 54; `test-omics-norm-state-guards.R` 647 to 651; `test-workflow-partitions.R` 439 to 440;
  `test-workflow-operation-registry.R` 1095/2 to 1107/2; `test-metab-qc-shared.R` 20 to 23;
  `test-integration-record.R` 167 to 200.
- 5 `message-changed`: `test-prot-norm-core.R` (1548/306 on both sides, 18 labels changed in message);
  `test-workflow-saved-da.R` 776/175 to 755/196; `test-workflow-omics-binding-steps.R` 58/1 to 57/2;
  `test-workflow-omics-declarations.R` 45/23 to 43/25; `test-omic-da-module.R` 1289/80 to 1288/81.
- 4 `new-failure`: `test-omic-qc-module.R` 160/0 to 156/4; `test-gui-entry-proteomics.R` 168/25 to 92/101;
  `test-gui-entry-lipidomics.R` 35/0 to 19/16; `test-gui-entry-metabolomics.R` 105/0 to 57/48.

Every other file of the 516, **489**, is predicted `same`, at the previous checkpoint's counts. This includes the
pre-existing failing files, `test-kit212-da-current.R`, the other integration files and the three trx module-CI
files.

Known flaky leaves, named in advance. A flake is recorded as such and is never re-run away:
- `tests/testthat/test-workflow-figure-paths.R`: a failure whose raw message quotes a random temporary folder whose
  name contains "259".
- `tests/testthat/test-gui-entry-metabolomics.R`: the random 9-hexadecimal-character export id escapes the GUI-entry
  mask (`tests/testthat/helper-gui-entry.R:60`, which masks runs of 10 or more); train 16 saw it on a base side. It
  would show as a changed message with the same 48 failures.
- Plan section 11's `gprofiler_plot` and KIT-201 I13's `refusals` leaves in `test-omic-enrich-module.R` (the
  previous checkpoint's flakes).

Known stored-reference failures that move only by their declared lines (end-of-project re-record) are inside the
27 above: the gui-entry, saved-da, declarations, binding-steps and prot-norm-core files.

R CMD check is predicted to give identical findings on both sides, matching trains 10 to 18 and the previous
checkpoint: 0 ERRORs, 3 WARNINGs, 3 NOTEs, with the same count of paths (238 at the previous checkpoint) in the
portable-file-names NOTE on both sides. The head adds paths under `docs/qa/` and `tests/`; the NOTE is read for any path
new since 828dd105.

## Decision rule

As the previous checkpoint's:
- A file that differs from the base side or fails, and is not in `predictions.tsv`, is UNPREDICTED. It is listed
  with its failure text, not interpreted, not fixed and not re-run.
- A predicted file whose difference goes beyond its declared kind or count is UNPREDICTED.
- An unfavourable valid result is never re-run.
- An instrument defect with no result judged may be re-sent under VAL-54, at most three times, and is recorded.
- Nothing is fixed here.

## Result

Run 45712202804 (request `checkpoint-20261126b-15d86bd4`, run-tests commit 9ab38c67, dispatched normally by the
push) was sent once and finished with conclusion success. Package commit 15d86bd4, base 828dd105, 516 files, 20
shards per side, no scripts job. 45 jobs used 1,325.4 job-minutes (88.4 core-hours) of the 1,692.4 ceiling (112.8
core-hours); the check jobs 26.9. Tables beside this file: `compare.tsv`, `tests.tsv`, `tests-base.tsv`,
`check-compare.txt`, `check-summary.txt`, `meta.txt`. The one `gh run watch` ended early on an HTTP 403 (rate limit,
not a run failure); the run was read afterwards by status only.

Cross-check: the base side of this run equals run 45644467847's head table in all 505 files (counts, failed, errors,
skipped), so the environment is unchanged.

### Verdicts, 516 files

| Verdict | Files | Predicted |
|---|---|---|
| same | 488 | 489 |
| only-on-head | 11 | 11 |
| count-changed | 7 | 7 |
| message-changed | 5 | 5 |
| new-failure | 5 | 4 |

All 27 predicted files differ, each at exactly its predicted head count (see `predictions.tsv`). One further file
differs: UNPREDICTED.

### Unpredicted difference

`test-omics-module-globals.R`: base 14/0/0/0, head 13/1/0/0, `new-failure`. Failing test: "KIT-030: the
super-assignments in the lipid and metab code are the listed local closures". Message (`unpredicted-omics-module-globals.txt`):

    Expected `setdiff(files_hit, allowed_files)` to be identical to `character(0)`.
    `actual`:   "mod_omic_norm.R"
    `expected`:  (empty)

Not interpreted, not fixed, not re-run. No train scope held this file.

### Flakes

Neither named flake appeared: `test-workflow-figure-paths.R` is `same` 181/0; no 9-hexadecimal-character export id
appears in the failures of `test-gui-entry-metabolomics.R` on either side (its counts, 57/48, are as predicted).

### R CMD check

0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; all 6 findings identical, none only on head or only on base. This is
the previous checkpoint's result. The portable-file-names NOTE lists the same paths as the previous checkpoint's
(240 by a line count of both `check-summary.txt` files, equal).
