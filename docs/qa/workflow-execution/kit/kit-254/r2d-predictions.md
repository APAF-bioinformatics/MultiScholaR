# KIT-254 statement's correction run R2d: predictions, written before the request

Base db727407 (kit-254 with staging 828dd105 merged; the statement code of 904a5514; the tests of this correction:
E7 at t_hi 0.960024, E17 both ways). Head 267c9507 (the one-place fix in `workflowReplayDesignSamples`). Scope (the
coordinator's ruling): the files that reach `workflowReplayDesignSamples`, `correlation_statement` or
`workflow_replay_correlation_skip` by grep on the merged tree (only the new file), `test-workflow-analysis-record.R`,
`test-omics-norm-state-guards.R`. No R CMD check (the merge train repeats the full scope and check). Ceiling 150
job-minutes.

- New file, head: 0 failures. E7's restored sample assertions pass: Apply(0.960024) on lipid with RUV removes **A_1
  only** (PREDICTED from G0b's best partners: A_1 0.957297 below, B_4 0.962751 and every other sample above). E17 "pos"
  names only LCMS_Pos and E17 "both" names LCMS_Pos and LCMS_Neg, lipid and metab.
- New file, base: fails exactly E17 "pos" for lipid and metab (2 failures: base names "LCMS_Pos, LCMS_Neg"). E17
  "both" passes on base (its wrong rule happens to be right when the sample fails in every assay). E7's restored rows
  behave as on head (base has the main fix and the statement; the one-place fix touches only the assay detail).
- `test-workflow-analysis-record.R`, `test-omics-norm-state-guards.R`: `same` (651/0 for the guards file).
- If E7's precondition shows another sample than A_1, that is a result (recorded, not adjusted).
