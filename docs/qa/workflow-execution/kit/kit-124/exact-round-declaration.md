# KIT-124 exactness round: declaration (2026-11-18, written before any runner job)

Authority: the user's decision M16 of 2026-11-18 (`design-review-20261118.md` section 6 on `performance-v2`): four changes to the recommended RUV design are measured, each as its own arm and all together, before any of them is put in the package. Simulations run on GitHub-hosted runners (M13). Nothing under `R/`, `tests/`, `DESCRIPTION` or `tools/` changes; the arms are computed in the simulators from limma, sva and the package's own kernels.

## Purpose

To give, for each change, the figures on which the user decides whether it goes into the package: its effect on calibration by the ticket's rules 1 to 7 (`KIT-124-rev1.4-ticket-20261118.md`, W3), and whether the specific defect it addresses is reduced.

## Arms

Every arm is an existing arm with the named change and nothing else; it runs on that arm's random stream, on the same data sets and seeds as the confirming study. Arm numbers 17 to 24 were taken.

| Arm | Is | Change | Simulators |
|---|---|---|---|
| 25 | arm 7 | (a) each comparison's standard error exact per feature: the comparison is fitted as a coefficient of a reparameterised model matrix (`limma::contrastAsCoef`) instead of being formed by `limma::contrasts.fit` from one correlation matrix for all features | S1, S2, S4 |
| 26 | arm 7 | (b) the classic proteomics test applies the variance trend and robust settings it sets: `limma::treat(fit, lfc = 0, trend = TRUE, robust = TRUE)` instead of `treat(fit, lfc = 0)` | S1 |
| 27 | arm 1 (no RUV) | (b) | S1 |
| 28 | arm 3 (today's design) | (b) | S1 |
| 29 | arm 7 | (c) the estimated count E computed after the biological units are removed (`sva::num.sv` with the unit as the model), not only the groups | S1, S2, S4; cells C6 and C8 |
| 30 | arm 16 (technical replicates as the replicate structure) | (c) | S1, S4; cells C6 and C8 |
| 31 | arm 7 | (d) E from 200 permutations instead of 20, the same seed | S1, S2, S3, S4 |
| 32 | arm 7 | (e) all together: (a), (b) and (d) where each applies, and (c) in C6 and C8 | S1, S2, S4 |
| 33 | arm 16 | (e) all together, on the technical-replicate design | S1, S4; cells C6 and C8 |

Where a change has no arm, and why:

| Change | Simulator | Reason |
|---|---|---|
| (a) | S3 | The package's test is edgeR's quasi-likelihood F-test of the contrast (`glmQLFTest(fit, contrast = )`), which fits the contrast itself per gene. No standard error is formed from a common correlation matrix; nothing is approximate |
| (a) | S2 | Has an arm, for a different cause than S1 and S4: the quantified matrix has no missing values, but `limpa::dpcDE` fits with per-observation weights, and with W in the model `contrasts.fit` is approximate under weights as it is under missing values |
| (b) | S2, S3, S4 | No TREAT step: S2's `dpcDE` branch and S4's kernel (a cutoff of 0 is passed on as NA) report the p-values of `eBayes(trend, robust)`; S3 is edgeR |
| (b), two variants | S1 | "treat() given the settings" and "no TREAT step at a cutoff of 0" give the same p-values (shown per data set to 1e-12 in the local check), so one arm, not two. The arm is computed as treat() given the settings |
| (c) | S3 | On the tree where technical replicates are handled (kit-128) the libraries of a unit are summed before anything else, so E already sees one sample per biological unit. C6 and C8 are still run there with arms 1, 2, 7 and 31 |
| (c), arm 30 | S2 | On the DIA route arm 7 already takes the declared technical replicates as the replicate structure (arm 16 repeats arm 7), so arm 29 is the variant of both |
| (e) | S3 | Would be arm 31; not run under a second number |

Arm 33 was not in the request: it is the all-together design on a declared technical-replicate design, added because that is the design such a data set is given by default (M6). It can be ignored.

With E = 0 an arm is the no-RUV result, as arm 7: arm 1's result for arms 29 to 31, and the no-RUV fit under the arm's own test for arms 25, 26, 32, 33. In S1 a refusal of arms 25, 26, 28, 32, 33 is scored as the no-RUV fit under the arm's own test; elsewhere as arm 1's result, as before.

## Cells and runs

The eight core cells C0 to C7, and one new cell: **C8, the technical-duplicate cell with no unwanted factor** (data generation as C6 with the number of factors set to 0; seed index 22 in S1, S3 and S4, 8 in S2). C8 is reported and judges no design; it is a no-factor cell (rule 5 applies, not rule 4). 1,000 runs per cell, runs 1 to 1,000, seeds as the confirming study (100000 + 10000 x cell index + run).

## Package trees

| Branch | Package tree | Simulators and cells | Why |
|---|---|---|---|
| `sim-k124-exact` (cut from `sim-k124` at `0baeb076`) | `0b6fcb6d`, the tree of the confirming study | S1 and S2: C0 to C8. S3 and S4: C0 to C5 and C7 | The laptop's confirm rows exist for this tree, so the gate can be read. Classic proteomics and DIA handle technical replicates on this tree (the test blocks on the unit) |
| `sim-k129-exact` (cut from `sim-k129-r0a` at `3db3a2ed`) | `ace1f4de`, kit-128 | S3 and S4: C6 and C8 | On `0b6fcb6d` the technical-duplicate cell is void in transcriptomics and lipid/metab (duplicates tested as independent samples); kit-128 repairs it, and run r0a (44203261310) measured C6 on it |

Simulator files: the commit of each branch that carries this declaration. `aggregate.R` on both branches is the version on `performance-v2` (sha256 `c5c6275f...`) with three constants changed so that C8 is a no-factor, reported cell.

## Gate

Before any new arm is scored:

- Tree `0b6fcb6d`, cells C0 to C7: arms 1 and 2 of the runner equal the laptop's confirm rows of the same seeds under the fixed rule (`packaging/verify/compare-results.R`: relative 1e-6, absolute floor 1e-9, discrete values exact), read on the values the study's rules use (null rates, false discovery proportion, power, k, counts), as M19 decides; bias and rmse are reported.
- Tree kit-128, cell C6: arms 1, 2 and 7 (S4 also arm 16) equal run r0a's rows under the same rule, read the same way. No laptop rows exist for that tree.
- C8 is new on both trees: there is nothing to compare it with, and that is stated with its figures.
- S1 arm 3 is run so that arm 28 has its base on the same machine; it is compared with the laptop's arm 3 and the result reported, not gated (its platform dependence is recorded in `gha-proof-20261118/README.md`).

A cell that does not meet the gate is not scored and is reported.

## Attempts and limits

One attempt per cell. A job that fails for an infrastructure reason may be rerun once as the same command (`gh run rerun <run> --failed`); seeds are per run, so a rerun recomputes the same rows. A job that reaches its time limit is reported, not relaunched with other settings. Per-job limit: 300 minutes for the simulator (the workflow's maximum; the job's own limit is 345 minutes). A cell may be split into run ranges (the workflow's `chunk`) where one job would not finish in the limit; the rows are the same. Sizing from the lifecycle run is added below before the launch; if the round exceeds about 120 core-hours it is not launched without the user's word.

## Reading, fixed in advance

Rules 1 to 7 are the ticket's, computed by `aggregate.R` with their limits unchanged. "No worse" means: in every judged core cell, on every rule, the arm's result (pass, inconclusive, fail) is not worse than its base arm's in the same cell and subset. No other threshold is introduced.

| Change | Adopted if | Also reported | Rejected if |
|---|---|---|---|
| (a) arm 25 against arm 7 | No worse than arm 7 (nothing to tune) | The change in the type I error at 0.05 and 0.01 and in the false discovery proportion, per cell | Worse than arm 7 on a rule in a judged core cell |
| (b) arm 26 against 7, 27 against 1, 28 against 3 | No worse than the base arm | The same differences per cell; the restated figures of arms 1 and 3 | Worse than the base arm on a rule in a judged core cell |
| (c) arm 29 against 7, 30 against 16 | In C8 the mean E and the share of data sets with E above 0 are lower than the base arm's; in C6 no worse than the base arm (where C6 is a judged cell: S1, S2); in C8 rule 5 no worse than the base arm's | E's distribution in C8 and C6; calibration and power against arms 1 and 2 in both cells | E in C8 not lower than the base arm's, or worse than the base arm on a rule in C6 (S1, S2), or on rule 5 in C8 |
| (d) arm 31 against arm 7 | No worse than arm 7. The seed-dependence of the on/off verdict cannot be measured here (one seed per data set) and is not part of the reading | E's distribution per cell, the share of data sets with E above 0 in C3, calibration | Worse than arm 7 on a rule in a judged core cell |
| (e) arm 32 against arm 7 (33 against 16) | No worse than the base arm | The same figures as above | Worse than the base arm on a rule in a judged core cell |

Anything between adopted and rejected (for example a rule that moves from pass to inconclusive within the Monte Carlo error) is reported as such and is the user's decision. C6 and C8 in S3 and S4 are reported cells and judge no design.

## Local check before the runner (laptop, R 4.6.1, one R process, niced)

| Check | Result |
|---|---|
| Existing arms unchanged, tree `0b6fcb6d`, runs 1 to 3: S1 arms 1, 2, 3, 7, 16 in C0, C3, C6; S2 arms 1, 2, 7, 16 in C0, C3, C6; S3 arms 1, 2, 7 in C0, C3; S4 arms 1, 2, 7 in C0, C3, C7 | Equal to the laptop's confirm rows under the fixed rule; largest difference 0 |
| Existing arms unchanged, kit-128 tree, C6, runs 1 to 3, against run r0a | S3 arms 1, 2, 7 and S4 arms 1, 2, 7, 16: within the rule (largest relative difference 1.1e-9). S4 arm 15 (not used here) differs from the runner's rows in one of three runs; r0a's own script gives the same rows as the new one on the laptop, so this is laptop against runner, not the edit |
| The rebuilt S1 kernel with neither change | `identical()` to the package's kernel, with and without W, in every data set checked (`K124_EXACT_CHECK=1`) |
| treat() given the settings against no TREAT step, S1 | p-values equal to 1e-12 in every data set checked |
| The reparameterised fit without W, S2 and S4 | Estimates and p-values equal to arm 1's to 1e-8 in every data set checked |
| New arms | Every row `ok`; scored by `aggregate.R`; C8 scored under rule 5 as a reported cell |
| Deliberate failure | An unknown cell stops the simulator with exit status 1 and writes nothing |

## Sizing and launch

To be completed from the lifecycle run before the launch.
