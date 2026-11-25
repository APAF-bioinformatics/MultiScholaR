# KIT-266 gate G0: declaration (written before the run)

**Purpose.** Observe the BASE (no product code, no test expectation, no stored file touched) for KIT-266 (r1.3, section 11): why the DIA-with-limpa route builds no combined QC figure, and the memory cost of building one in a fresh process.

**Bound commit.** Staging `894b9cc82` (the owner's decisions VAL-314 to VAL-318). This branch's only difference from it is the files under `docs/qa/workflow-execution/kit/kit-266/g0/`; `R/`, `inst/`, `tests/` are byte-identical (check: `git diff 894b9cc82 HEAD --stat`). The reuse condition of section 11 applies to any later branch point.

**Files.** `g0.R` (driver; the lane `guiEntryRunLane("prot_dia_limpa")` observed through traces of `guiEntryArtifacts`, `runProtNormTabEntryWorkflow`, `copyToResultsSummary`, `generateProtNormPreNormalizationQcArtifacts`; all trace output written to absolute paths under `kit266-g0/`), `memory-build.R` (one fresh process per figure). Outputs: `controls.txt`, `events.log`, `events.tsv`, `steps.tsv`, `predictions.tsv`, `oracle-compare.tsv`, `memory-*.txt`, `memory-readings.txt`.

**Positive controls (before any reading is trusted).** Each trace is exercised on a stand-in call and must record an event and a log line; the memory measure is exercised by a fresh process that touches a 0.224 GiB vector (pass: rise 0.18 to 0.60 GiB). A control that records nothing means the instrument is blind: the lane is not run, stop and report.

**Predictions (section 11).**
- P1 at step `norm` only `post_norm_*` exist (of the pre/post/ruv/composite figures).
- P2 `pre_norm_*` first exist at step `ruv_covariates`.
- P3 no composite figure file at any step.
- P4 the final `pre_normalization` source is exactly the `correlation_filtered` state's entity.
- P5 `copyToResultsSummary`'s returned list names "Composite QC (PDF)" and "Composite QC (PNG)" as missing (traced return value).
- P6 the tab entry (`runProtNormTabEntryWorkflow`) fires at step `ruv_covariates` with `current_state == "correlation_filtered"`.
- P7 the Reset fallback (`normalization_input_reference` NULL; `revertProtNormStateManagerToPreNormalization`) selects `correlation_filtered`.
- Per-step record against `tests/testdata/gui-entry/prot_dia_limpa.json`: INFORMATION ONLY (KIT-271 moves 13 lines of that oracle).

**Memory rule.** In a fresh R process with the package loaded, the three-column figure (12 stage images and `ruv_covariates_cancor.png`, the existing builder `generateProtNormCompositeQcFigure` with a two-line local resolver) and, in its own fresh process, the two-column figure (`ruvMode = "skip"`). Rise = VmHWM after the build minus VmRSS just before it (`/proc/self/status`, as KIT-122); every `Vm*` line printed verbatim. Declared prediction: a rise of 0.3 GiB or more for the three-column figure. This run only REPORTS the figure; the coordinator rules (0.3 GiB or more moves the managed limpa build to the worker, section 4 item 10).

**Stop rules.** Any of P1 to P7 false: stop rule 1 (record and report, no interpretation). A blind control: stop and report. Instrument defects (driver mechanics, no observation) may be fixed and re-sent under the brief's four conditions, at most three; an observed-false prediction is a result and is never re-run.

**Budget.** One scripts job (4 cores): expected 25, ceiling 45 job-minutes (script timeout 40); core-hours = job-minutes x 4 / 60 (expected 1.67, ceiling 3.00).

## Instrument defect 1 (run 45711141374, 2026-11-26), recorded under VAL-54
The run stopped at the positive controls, by design: the artifacts trace fired and wrote (events.log: `type=capture; n_qc=2; qc=post_norm_pca.png,pre_norm_pca.png`) but the step id read `NA`, because `parent.frame()` inside a trace's tracer does not return the caller. Mechanics only: no lane reading was taken. Fix: the capture frame is found by scanning `sys.frames()` for the frame holding `step` and `refusal_before`. No prediction, threshold or stop rule changed. Records: `run1-45711141374/`.

## Result record (run 45713081362, commit 6dc54955; files under run2-45713081362/)
All five positive controls PASS; P1 to P7 as predicted (predictions.tsv). Memory: three-column rise 0.914 GiB, two-column 0.539 GiB (memory-readings.txt). Caveat on the informational per-step comparison: P7 moves the head at the `session_export` capture, so the comparison is valid up to `ruv_covariates` (all 16 steps identical to the oracle) and the later steps differ partly because of the observation itself.
