# KIT-151 d, Q1 results: which RUV output feeds integration, scored once

Declaration `ruv-declaration.md` (9bce6039, pushed alone). Lifecycle 44785248037 (fc5e8042; strength rule chose C1:
the ideal arm recovered S in 10 of 10 data sets for both candidates at both N, `ruv-lifecycle-candidates.tsv`).
Counted run 44786899460 (92ff4512): 12 jobs, 720 of 720 trainings ok, 0 failed, 0 at the 10,000-round limit (most
2,451). Job wall time 7,902 s: **8.8 core-hours** (ceiling 25). Processors: AMD EPYC 7763, 9V45, 9V74; Intel Xeon
Platinum 8370C. sva 3.60.0; RUVIII_C_Varying against ruv::RUVIII largest difference 1.7e-8. Tables: `ruv-rows.tsv`,
`ruv-summary.tsv`, `ruv-pairs.tsv`, `ruv-check.tsv`; scorer output `ruv-score.md` (unchanged). No correction pass.

Design (as declared): two omics (1,000 and 500 features), group effect G (weight 1), planted shared within-group factor
S (weight 1.5, 20% of features in both), one unwanted factor per omic U (weight 2, 40%); `main`: U unrelated to the
group; `conf`: U correlated 0.6 with the group. 30 data sets per N and scenario.

**The calibrated default's count gave k = 2 in both omics in all 120 data sets** (none refused), so the calibrated
arm is the same model as RUV-III k = 2 in every data set.

## Results (S recovered = matched |r| >= 0.9, of 30; medians)

| N | scenario | arm | S recovered | G recovered | group factor VE % (mean of omics) | R2 vs group | share of S left | share of G left |
|---|---|---|---|---|---|---|---|---|
| 20 | main | raw | 30 | 30 | 7.8 | 0.994 | 1 | 1 |
| 20 | main | ideal | 30 | 30 | 15.4 | 0.995 | 1 | 1 |
| 20 | main | k = 1 | 29 | 30 | 16.2 | 0.995 | 0.985 | 0.999 |
| 20 | main | k = 2 = calibrated | **0** | 30 | 23.3 | 0.996 | **0.007** | 0.998 |
| 20 | main | k = 3 | 0 | 30 | 24.9 | 0.996 | 0.006 | 0.997 |
| 20 | conf | raw | 30 | **6** | 5.4 | **0.76** | 1 | 1 |
| 20 | conf | ideal | 30 | 30 | 15.4 | 0.995 | 1 | 1 |
| 20 | conf | k = 1 | 16 | 19 | 17.6 | 0.86 | 0.966 | 0.999 |
| 20 | conf | k = 2 = calibrated | **0** | 30 | 25.4 | 0.996 | 0.007 | 0.999 |
| 20 | conf | k = 3 | 0 | 30 | 27.1 | 0.997 | 0.006 | 0.997 |
| 60 | main | raw / ideal | 30 / 30 | 30 / 30 | 6.8 / 13.3 | 0.996 / 0.996 | 1 / 1 | 1 / 1 |
| 60 | main | k = 1 | 30 | 30 | 13.6 | 0.996 | 0.991 | 0.999 |
| 60 | main | k = 2 = calibrated / k = 3 | 0 / 0 | 30 / 30 | 18.8 / 19.3 | 0.997 / 0.997 | 0.004 / 0.004 | 0.998 / 0.996 |
| 60 | conf | raw / ideal | 30 / 30 | 30 / 30 | 9.4 / 13.3 | 0.983 / 0.996 | 1 / 1 | 1 / 1 |
| 60 | conf | k = 1 | 28 | 28 | 14.1 | 0.987 | 0.981 | 0.999 |
| 60 | conf | k = 2 = calibrated / k = 3 | 0 / 0 | 30 / 30 | 19.6 / 20.0 | 0.997 / 0.997 | 0.004 / 0.004 | 0.998 / 0.996 |

Paired against ideal (mean [95% seed-bootstrap interval]):

| N | scenario | arm | abs r S | abs r G | group VE (points) | R2 | -log10 p |
|---|---|---|---|---|---|---|---|
| 20 | main | raw | -0.001 no practical difference | -0.000 | -7.5 [-7.6, -7.5] | -0.000 | -0.18 [-0.49, +0.12] |
| 20 | main | calibrated | **-0.988 [-0.990, -0.986] worse** | +0.001 | **+7.9 [+7.9, +8.0]** | +0.002 [+0.001, +0.004] | +1.5 [+0.8, +2.3] |
| 20 | conf | raw | -0.004 no practical difference | **-0.161 [-0.191, -0.134] worse** | -10.5 | **-0.289 [-0.336, -0.246]** | -15.2 |
| 20 | conf | k = 1 | -0.106 [-0.142, -0.074] worse | -0.095 [-0.129, -0.063] worse | +3.5 | -0.171 | -10.6 |
| 20 | conf | calibrated | **-0.987 [-0.989, -0.985] worse** | +0.001 no practical difference | +9.9 [+9.7, +10.1] | +0.003 | +2.1 |
| 60 | main | calibrated | -0.991 worse | +0.000 | +5.5 [+5.5, +5.6] | +0.001 | +2.1 |
| 60 | conf | raw | -0.000 | -0.006 no practical difference | -4.1 | -0.013 [-0.014, -0.011] | -17.3 [-18.8, -15.8] |
| 60 | conf | calibrated | -0.991 worse | +0.000 | +6.2 | +0.001 | +2.8 |

(Full table, against raw too: `ruv-pairs.tsv`.)

## What is shown

- With the groups as replicate sets and every feature a control, the package's calibrated default (k estimated at 2
  in every data set here, matching the two within-group axes per omic: U and S) **removed the shared within-group
  factor in every data set** (S recovered 0 of 30 in all four N x scenario cells; 0.4 to 0.7% of its signal left in the
  corrected data), while keeping the group effect (99.8% left). k = 1 took the unwanted factor and kept S (29 to 30 of
  30 unrelated; 16 and 28 of 30 when confounded). The ideal correction found S in 30 of 30 everywhere, so the design
  can show removal (baseline met).
- Unconfounded (`main`), the uncorrected input was as good as the ideal one for S and for the group factor's |r|
  and R^2; only the group factor's variance explained is lower (7.8% against 15.4% at N = 20), because the unwanted
  variation stays in the denominator.
- Confounded (`conf`, U correlated 0.6 with the group), at N = 20 the uncorrected input distorts the group factor
  (recovered in 6 of 30; R^2 0.76 against 0.995); the calibrated correction restores it to the ideal level (|r| and R^2
  equal to ideal) but still removes S. k = 1 does neither reliably. At N = 60 the distortion of the uncorrected input
  is small (|r| -0.006, R^2 -0.013) but the group p-value is weaker by 17 orders.
- Against the ideal input the corrected input raises the group factor's variance explained by 5.5 to 9.9 percentage
  points (x1.4 to x1.6) and R^2 by at most 0.003; this is the variance of S (and of the remaining noise) being taken
  out of the denominator, not group signal added.

## What is not shown

One simulated design: S weaker than each omic's unwanted factor and stronger than the group effect, both within-group
axes counted by the estimate; Gaussian, linear, complete data, two groups, two omics, one training seed. If S were
stronger than U, k = 1 would remove S instead; if the count stopped at 1, the calibrated arm would behave like k = 1.
Real data sets, technical-replicate designs, supplied control lists and the other omics' kernels (transcriptomics'
`ruv::RUVIII`, which agreed to 1.7e-8 here on complete data) were not run. The package itself was not installed: the
count is `sva::num.sv` with the package's arguments and the correction `RUVIIIC::RUVIII_C_Varying`, as the
declaration states.

Methods sentence the result supports: "In simulations, the default RUV-III design (groups as replicate sets, every
feature a negative control, number of factors estimated) removed a planted biological factor shared across omics
that varied within groups (recovered in 0 of 30 data sets against 30 of 30 with the true unwanted factor removed) and
raised the group factor's share of variance explained; it restored the group factor when the unwanted factor was
correlated with the groups. Integration on corrected inputs therefore cannot be expected to recover within-group
structure shared across omics."
