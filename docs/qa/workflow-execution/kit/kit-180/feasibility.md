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

## F2: does a label-free and a TMT record replay on base? Yes, in order; the replayed states are not identical().

- R0 (run 45120111092): `replay-check.R` stopped on both lanes at its own check `length(stores) == 1L` (the session
  owner's store is not under the lane's root). Instrument failure, no result (`r0/replay.tsv`).
- R0b, the gate's correction pass (run 45122598592): both lanes stopped with `invalid 'file' argument`, no location,
  no result (`r0b/replay.tsv`).
- R0c, instrument run 1 under VAL-54 (run 45147720852): the script now records where a lane stops:
  `file.exists(sources)` with `sources` NULL (no recorded input carries a source path). No result (`r0c/replay.tsv`).
- R0d, instrument run 2 under VAL-54 (run 45149789668, run branch `run-k180-base` c7251d9c, `R/` equal to base
  a60c4e9a), the result (`r0d/replay.tsv`). On `prot_lfq` and `prot_tmt` the GUI lane's record (import, design,
  the four protein QC steps, normalisation with RUV skipped, correlation filter; valid) replayed into a copy cut to
  `initial` and `protein_s4_initial`, with the inputs checked by SHA-256 (`maxquant_report` or `pd_tmt_report`,
  `fasta`): the replayed state history equals the record's (8 states) on both lanes. **None of the 6 later states is
  `identical()` to the first run's, on either lane.** Not diagnosed: no comparison by digest or by slot was made, so
  whether the difference is in the data or only in object attributes (for example `@args`) is not known. Observed on
  base, before this ticket's change; the replay reads no chain rule or tab status (2.5), so this ticket does not
  change it. A result, so not re-run.
- What follows for T5: F2's stop rule ("not supported") does not apply, because the records do replay. T5 (an LFQ
  record made in process without QC steps) passed on the tests-only commit (R1) and on head (R2, R5), with every
  replayed state `identical()`. The ticket's claim "a replayed record repeats the steps it recorded" is shown (step
  order); "every replayed state identical" is shown for T5's in-process record only, not for GUI-made records.

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
