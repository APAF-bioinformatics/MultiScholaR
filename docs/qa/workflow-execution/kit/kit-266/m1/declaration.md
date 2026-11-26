# M1 declaration (before the run)

Purpose: the controller process's memory rise at the delivery of each limpa_composite action (the two-column build at the
pass-through, the three-column build at the covariates step) is under 0.3 GiB (A5a, KIT-122's figure applied to the moved
event). Binding: the fix commit on kit-266 plus this script. One scripts job, ceiling 60 job-minutes, expected 30.
Prediction: rise under 0.3 GiB at both points; the reset control passes (VmHWM within 5 percent of VmRSS after the reset);
the path control passes (`writeDIACompositeFigureFiles` TRUE at both, no fallback logged, composite PNG 1500 x 3300 after the
pass-through and 1950 x 3900 or 1950 x 3300 after the covariates step); the positive control reads 0.3 GiB or more (G0:
0.539 and 0.914). Decision: a rise of 0.3 GiB or more at either point is a stop (rule 7); a failed control is an instrument
defect (VAL-54, re-send 2 of 3 in total across R1 and M1).

## Result of the first send (run 45743367152, 13.4 job-minutes)

Reset control passed; positive control 0.540 (two columns) and 0.914 GiB (three). The PATH control failed: `writeDIACompositeFigureFiles`
returned FALSE at both build points and the fallback was logged, so the controller built the figure itself (rises 0.558 and 0.762
GiB are the fallback's, not a measurement of the worker branch). Stop rule 7: a failed path control is an instrument defect if the
cause is the instrument; the cause is not yet known (this script loads the package with `load_all` while the child runs the
installed one; or the sibling function's save/root does not reach the controller's store). Before re-sending M1, a head check of
the new test file (R1c, below) shows whether the in-process double finds the worker's files saved (E16 to E18, E20, E21).
