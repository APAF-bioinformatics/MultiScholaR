# M1 declaration (before the run)

Purpose: the controller process's memory rise at the delivery of each limpa_composite action (the two-column build at the
pass-through, the three-column build at the covariates step) is under 0.3 GiB (A5a, KIT-122's figure applied to the moved
event). Binding: the fix commit on kit-266 plus this script. One scripts job, ceiling 60 job-minutes, expected 30.
Prediction: rise under 0.3 GiB at both points; the reset control passes (VmHWM within 5 percent of VmRSS after the reset);
the path control passes (`writeDIACompositeFigureFiles` TRUE at both, no fallback logged, composite PNG 1500 x 3300 after the
pass-through and 1950 x 3900 or 1950 x 3300 after the covariates step); the positive control reads 0.3 GiB or more (G0:
0.539 and 0.914). Decision: a rise of 0.3 GiB or more at either point is a stop (rule 7); a failed control is an instrument
defect (VAL-54, re-send 2 of 3 in total across R1 and M1).
