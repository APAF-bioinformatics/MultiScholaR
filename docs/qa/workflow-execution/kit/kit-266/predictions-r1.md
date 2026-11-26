# R1 declaration (before the run)

Purpose: the tests-only commit shows every new E row failing by assertion on the code before the change, and every G row
passing. Binding: the commit that adds the tests and the declared edits (nothing under R/ or inst/); base `none`. Files run:
`test-kit266-limpa-combined-qc-figure.R`, `test-kit259b-second-correction.R`, `test-kit227-report-ruv-sentences.R`. One job,
ceiling 45 job-minutes (3.00 core-hours), expected 25 (1.67). Decision: E rows failing at the stated values and G rows passing
=> A1 met and R2 may be sent; an E row that passes, or a G row that fails, on this code is a stop (rule 2); an error (not a
failure) in a test body, or a crashed file, is an instrument defect if no result was produced (VAL-54), else a result.

## New file, by row (head = the old code)

Failing by assertion: E1 (in process), E1 (managed), E2, E3, E4, E5, E6, E7, E8 (in process), E8 (managed), E8b, E9, E9b,
E10, E11, E12 (a), (c), E12 (b), E13, E13 (caption; its classic-template expectations pass), E14, E15 x4 (apply/skip,
in process/managed), E16, E17, E18, E19 (fails on `is.function(saveWorkflowDiaLimpaCompositeFigure)`), E20, E21.
Values: no composite file (FALSE x3), 0 builder calls, 0 actions of the new kinds, `resolveDiaLimpaCompositeFigureInputs` and
`buildDiaLimpaCombinedQcFigure` absent, the old covariates summary text, the old limpa Results text (three stages for k = 0,
"at each stage of processing" for a record without `ruvcov_k`), the old caption.
Passing: G1, G2, G2b, G3, G4, G5, G6, and the three classic-template expectations of E13 (caption).
Errors: none predicted. Skips: none predicted (limpa, mixOmics, qs2 and duckdb are in the runner library).
An E row that "passes on base" is a stop: none is predicted to. E2, E3, E13 and E12 (b) also skip after their failing
assertion (testthat reports the failure and the skip).

## test-kit259b-second-correction.R (declared edits)

5 base-only failures: C01-G1 (2 expectations: the list and the composite exist), the DIANN_limpa iteration of (ii)-E3 (3
expectations: limpa's stage words and the `ruv_covariates_cancor` panel m). Every other test of the file: as on staging.

## test-kit227-report-ruv-sentences.R (declared edit)

1 base-only failure: E7 `skipped` (the two-stage limpa text of O-4). Every other test of the file: as on staging.

## Not run here

`test-prot-norm-core.R` (W3 limpa block `:1640-1672`, A111-01 `:1776-1815`) and `test-prot-norm-worker-figures.R` (`:136`,
`:137`, `:144-145`) carry the declared edits whose values are in the ticket (section 8 item 6); they are evidenced at R2
(base side). Predictions for them: base-only failures on exactly those assertions (`test-prot-norm-core.R` 2: the W3 limpa
operations lists and the `composite_files` set-aside count the head makes, the A111-01 length and operations list;
`test-prot-norm-worker-figures.R` 3); a base-only failure elsewhere is a stop.

## R1 result and the instrument re-send (VAL-54, re-send 1 of at most 3 across R1 and M1)

Run 45737408677 (27.8 job-minutes): new file 44 passed / 101 failed / 0 errors / 6 skipped; kit259b 5 failed; kit227 1 failed, as
predicted above. Every E row failed at the stated value by assertion; G1 to G6 and the classic expectations passed. The six
skips are the five rows that skip after a failed assertion (E2, E3, E13, E12 (b), E19) and E11. One defect of the mechanics: E11's
export failed in the test's own helper ("This does not look like a tar archive", the tab's export function on a limpa
project), so it did not fail for its stated reason (nothing to carry). Fix (no expected value, threshold or comparison
changes): E11 uses the native export as KIT-271's test runs it. Re-send: the new file only, to observe E11 (and all rows once
more) on the same code; prediction: as above, E11 failing at `file.exists(restored)` or at the byte comparison, not at the
export. Ceiling 45 job-minutes, expected 15.
