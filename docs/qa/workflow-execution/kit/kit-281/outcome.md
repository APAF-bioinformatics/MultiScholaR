# KIT-281 outcome (gate scored once, run 45770099879, package commit 3b9e6288; scorer output files unaltered in `results/`)

## Decision, as the scorer prints it
**`KIT-281 score: outcome Calibrated`** (`outcome.tsv`: outcome Calibrated, counted TRUE, bootstrap_B 10000, range_source g3-table; consequence line of the scorer: "Part C: the tab default becomes 1 (only after KIT-280 is on staging)"). Nothing is interpreted beyond the declared rule (declaration-counted.md).

## Deciding quantities (the twelve settings pooled, R = 1,260 replicate units = 4 cells x 315; bootstrap limits one-sided 97.5 %, B = 10,000; jackknife-normal beside, reported only)
| Item | Estimate | Bootstrap upper (decides) | Limit | State | Jackknife upper (reported) |
|---|---|---|---|---|---|
| rho_g (alpha_B1(H_B1 without H_B) / alpha_B(H_B)) | 1.0535 | 1.0716 | 1.10 | pass | 1.0718 |
| G1a FDP_BH(B1) | 0.04123 | 0.04170 | 0.05 | pass | 0.04171 |
| G1b FDP_BH(B1) - FDP_BH(B) | 0.0000555 | 0.000132 | 0.0025 | pass | 0.000133 |
| G2 BH power(B) - BH power(B1) | -0.004722 | -0.004512 | 0.005 | pass | -0.004396 |
F5' (B1.m6_l_H): estimate -0.00333, two-sided 95 % interval [-0.01002, 0.00351]: does not fire. Resamples dropped (not finite): 0 for every row. F1 0, F2 0, determinism 5/5, completeness 15,120 of 15,120 rows, `bootstrap_B` 10000, `range_declared` g3-table (`checks.tsv`, all pass).

## q-method as reported
0 metrics rows record a q-method other than qvalue (`qmethod-fallback.tsv`, checks `qmethod` 0): reported, the decision reads no q-value.

## Reported only (decide nothing)
- Per design, rho_g [two-sided 95 % limits], flagged (ell = 1) and other (ell = 0) proteins: A 1.0686 [1.0429, 1.0945], flagged 1.7498, other 0.9943; B 1.0311 [0.9884, 1.0747], flagged 1.4806, other 0.8863; C 1.0229 [0.9712, 1.0754], flagged 1.5463, other 0.7991; D 1.0281 [0.9839, 1.0720], flagged 1.5344, other 0.8067. Pooled: 1.0535 [1.0353, 1.0716], flagged 1.5994 [1.5474, 1.6520], other 0.9396 [0.9209, 0.9587].
- Tested added null proteins per masked data set (pooled 72.92; A 173.2, B 43.2, C 29.9, D 45.4).
- Phi pooled: B 0.9995, B1 0.9948. Positive control (Phi_S, reported, decides and stops nothing): 1.5162, lower limit 1.4718 against 1.20: passes.
- Per setting (labels at one-sided 0.05/12): rho_g estimates 0.989 to 1.087; 11 labelled "not determined", 1 "within margin" (D at vd25) (`labels-settings.tsv`).
Every claim sentence carries the per-design range: pooled over the twelve settings (the prespecified analysis) the ratio was 1.05 (upper limit 1.07; margin 1.10); per design the ratio was 1.02 to 1.07.

## Cost
Scoring run 13.0 min created to updated (script 38 s) = 0.87 core-hours (ceiling 2.67). Running total 21.73 + 0.87 = 22.60 of 74.67.

## Placeholder for the ledger row (decision text not drafted here)
`VAL-<merge>`: KIT-281 gate scored once, outcome as above; figures in `results/`; to be written by the coordinator.

## Not checked
No per-design decision (the gate decides on the pool only); no real data, random missingness, intensity-dependent variance, unequal groups, other percentiles; Part C not started (needs KIT-280 on staging); the transfer of the decision to the DIA route (peptide filter runs first); that the design effect of the added proteins in B, C and D equals Ag's (assumed). Nothing re-run or re-scored.
