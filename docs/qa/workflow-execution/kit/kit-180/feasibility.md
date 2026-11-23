# KIT-180 feasibility gates F1 to F4

## F1: can S1 drive the lanes without editing a helper file? Met.

Observed in R0 (run 45120111092, tests-only commit 09f242bb) and again in R3 (run 45124757845, fix eea8775b):

- S1 ran the three lanes in arms B and B-ruv (R0) and B, Q, B-ruv, Q-ruv (R3), every lane and arm `ok`
  (`r0/qc-effect.tsv`, `r3/qc-effect.tsv`, measure `status`).
- Mechanism: `pkgload::load_all(".", export_all = TRUE, helpers = FALSE)`, then
  `testthat::source_test_helpers("tests/testthat", env = globalenv())`, and the wrapped functions defined in the
  global environment (the pattern of `kit-101/diagnosis/f19_pipeline.R`). This differs in form from the ticket's
  wording (`helpers = TRUE` and the helpers' environment); the effect is the same and no helper file is edited.
- Each protein state's distinct protein ids and rows were read; the DA table and the RUV product were read from the
  owner; the step texts show the defaults applied (`r3/step-texts.tsv`: "Min replicates per group: 2", "Min groups
  required: 2", "Intensity cutoff percentile: 1.0%"). The provenance values read through
  `step_provenance[[state]]$values` came back empty (no `provenance:` rows); the step texts are the evidence.
- `f1.tsv` (both runs): a git worktree of the commit; a fresh `Rscript` process that starts with its working
  directory set to the worktree's root resolves `testthat::test_path()` inside that tree; `guiEntryExpectLane("prot_tmt")`
  inside `testthat::test_that()` returned TRUE with `MULTISCHOLAR_GUI_ENTRY_RECORD_DIR`, `..._ORACLE_DIR` and
  `..._FREEZE=1` set, wrote the record and the frozen file to the redirected directories, and left the checkout's
  oracle unchanged. The child process sets its working directory with `setwd()` as its first line rather than
  being started there; the effect is the same.

## F2: does a label-free and a TMT record replay on base? Not answered (instrumentation failed twice).

- R0 (run 45120111092): `replay-check.R` stopped on both lanes at its own check `length(stores) == 1L`: the session
  owner's store is not under the lane's root, where the script looked (`r0/replay.tsv`).
- R0b, the gate's one correction pass (run 45122598592, run branch `run-k180-base` c9ca1ed5 = 09f242bb with the
  corrected script, `R/` equal to base): the store root read from the owner; both lanes stopped with
  `simpleError invalid 'file' argument` after 160 s, before any replay result was written (`r0b/replay.tsv`).
  Not diagnosed further: the correction budget is spent.
- So F2's question is open: no GUI-made label-free or TMT record was replayed. This is an instrumentation failure,
  not the "not supported" outcome of the stop rule. T5 (an LFQ record made in process without QC steps, replayed
  into a copy holding the import and design states) passed on the tests-only commit (R1, run 45120113965) and on
  head (R2, run 45124754518); TMT replay is tested by nothing.

## F3: does `test-e2e-proteomics-tmt-lfq.R` run on run-tests? No.

R2 (run 45124754518): 3 skipped on base and head (no browser). The rewritten browser lane is shipped untested by a
run; route behaviour is evidenced by T1, T4 (the control) and T6 (the shared helper's dispatch, fake driver, passing
on head). The distributable bundle lane runs `proteomics-tmt-lfq` in Chrome at the merger's version bump; its outcome
is to be added to `runs.tsv` then.

## F4: reference rows for the three lanes. Met, before R4.

- Source `tickets/kit/evidence-run-44953187101/failures-base.txt`, SHA-256
  1b918062f8c94a3c10859a965f5f1120ea3c692ae70181039846c66a43730bde (checked with `sha256sum` on 2026-11-23; byte-equal
  to `failures.txt`), lines 1917-2431: 12 blocks, 4 per lane (`summary_save`, `summary_copy`, `summary_report`,
  `summary_export`).
- Each block names lane and step and lists only `artifacts` keys under `integration/proteomics/` present on one side
  only (head `.parquet`, frozen `.tsv`, three files), no `states` difference: class `integration_artifacts` (KIT-159
  r1.5 G1(a)), 0 rows of any other class.
- Committed: `g1a-reference-protein-lanes.tsv` (lane, step, class, differing fields and keys, the oracle's source
  SHA, source path, line range and SHA-256) and the twelve blocks as `run-44953187101-protein-lanes-failures.txt`
  (SHA-256 e9d5e7a38954c88d03fa678b7e7019b66039c9bb4df2f048e414a90ee95ebaa5).
