# KIT-124: RUV in every omic, confirming studies (core cells complete, 2026-11-18; stress cells running)

**All four simulators' core cells are scored** (S2 on 2026-11-18); the stress cells (conditions of use, never a verdict) are still running and are added when scored. Nothing here is the user's decision yet; the gate table at the end is filled in when S2 is scored. Rules: the ticket's W3 table (revision 1.2), applied by `aggregate.R` (`0bb67eb2`) to the pinned tree `0b6fcb6d` (KIT-123 R1 and R2 in). Scores of the core phases: scratchpad `k124-agg/s{1,3,4}-core/` (copied here with the final aggregation).

## The binding question before any per-omic row: rule 4 in the global-null cell (A12)

Rule 4 ("gain over no RUV in cells with a factor") cannot be met literally in C4, where nothing is differentially abundant: there is no power to gain, and the rule-3 reduction has nothing to reduce. Read literally, **every design fails C4 in S1 and S3, including the true-factor arm** (in S4, C4 is void; see below). Revision 1.4 puts both readings to the user.

Under the second reading, "inconclusive" means "as good as the true-factor arm": every inconclusive verdict in S1 and S3 comes only from C4, rules 2 and 7, where the true-factor arm itself sits at the limit (S1 arm 2 any-call 0.061, upper bound 0.078; S3 arm 2 0.075 with edgeR's quasi-likelihood test). In S3 no design can reach "pass" in C4 however many runs are made. S1's global-null excess sits in proteins with missing values (full 0.061 against complete 0.022), not in RUV.

| Reading | Effect |
|---|---|
| Literal: rule 4 applies in C4 | No design passes in any omic except arms 15/16 (technical replicates, S1 only). |
| C4 judged by rules 1, 2, 3, 6, 7 only (rule 4 undecidable there) | The verdicts in the column "without rule 4 in C4" below. |

## Per simulator (primary subset; failing core cells with the rules they fail)

### S1, classic proteomics (complete; 8 core cells, 202,000 rows)

| Arm | Design | Verdict (literal) | Without rule 4 in C4 | Failing core cells | Share FDP > 0.10 (mean over core) | Power C0 |
|---|---|---|---|---|---|---|
| 1 | no RUV | fail | fail | rule 4 in every cell with a factor | 0.006 | 0.49 |
| 2 | true factor in the model (ceiling) | fail | inconclusive | C4(4) | 0.028 | 0.91 |
| 3 | today (sweep, subtract, test) | fail | fail | C0, C2, C3, C4, C5, C7 (rules 1-3) | 0.417 | 0.92 |
| 4 | today's controls and k, W in the model | fail | fail | C5(1,2,3) | 0.140 | 0.89 |
| 5, 6 | data-driven controls, W in the model | fail | fail | C5(2,3) | 0.057 | 0.91 |
| 7 | all proteins as controls, k = E, W in the model | fail | fail | C1(3) | 0.051 | 0.91 |
| 11 | arm 7 with k = 1 | fail | fail | C1(3) | 0.052 | 0.91 |
| 8 | oracle null controls, W in the model (diagnostic) | fail | inconclusive | C4(4) | 0.035 | 0.91 |
| 14 | `sva` (complete proteins only; compare with arm 7 on the complete subset: power 0.95) | fail | inconclusive | C4(4) | 0.113 | 0.95 |
| 15 | technical replicates as RUV structure, subtract (run in C6 only) | pass (C6 only) | pass (C6 only) | none | 0.012 | n/a |
| 16 | the same, W in the model (run in C6 only) | pass (C6 only) | pass (C6 only) | none | 0.032 | n/a |

The "failing core cells" column lists the cells named by the aggregation for the primary subset; C4(4) applies to every arm under the literal reading (some rows omit it).

Plain reading: today's design is confirmed miscalibrated. All proteins as controls with the factors in the model (arm 7) misses one bound in one cell (C1, 12 per group: the share of runs above 10% FDP is 0.115 against a limit of 0.057; its mean FDP, 0.066, is within rule 2; it is genuinely liberal there, null rate 0.058). No RUV is conservative (null rate 0.014 to 0.015 in cells with a factor, except 0.034 with the factor tied to the groups, C5, and 0.028 with technical duplicates, C6; corrected 2026-11-18 from "0.010 to 0.016" against `scored/s1-core/k124-cells.tsv`) and loses 42 percentage points of power (0.91 to 0.49). **Arms 15/16 "pass" only because they were run in one cell (C6):** in C6 arm 7 does as well (mean FDP 0.054, power 0.896; arm 16 0.054, 0.896; arm 15 0.047, 0.871). Arm 7 fails only because it was also judged in C1 and C4. Reading "every core cell" as "every core cell where the arm runs" is aggregation note A14 (new; to the user): whether a design untested in C1, C4 and C5 can be eligible is the user's call. The order "15 > 16" is arm-number order; the ticket gives no order inside that tier.

### S3, transcriptomics (complete; 7 core cells judged, C6 reported)

| Arm | Design | Verdict (literal) | Without rule 4 in C4 | Failing core cells | Share FDP > 0.10 | Power C0 |
|---|---|---|---|---|---|---|
| 1 | no RUV | fail | fail | rule 4 in every cell with a factor | 0.000 | 0.005 |
| 2 | true factor in the model (ceiling) | fail | inconclusive | C4(4) | 0.011 | 0.77 |
| 3 | today (`trxRuvGrid`, RUVs, test) | fail | fail | C3(5), C5(1,2,3) | 0.143 | 0.73 |
| 7 | RUVs, all genes, groups as replicate sets, k = E | fail | inconclusive | C4(4) | 0.011 | 0.77 |
| 11 | arm 7 with k = 1 | fail | inconclusive | C4(4) | 0.011 | 0.77 |
| R | RUVr (residuals) | fail | fail | C5(1,2,3) | 0.154 | 0.77 |
| 14, 14E | `svaseq` | fail | fail | C3(5) | 0.011 | 0.77 |
| 5, 6, 8, 9, 10 | other control sets | fail | fail | see `k124-arms.tsv` | | |

Plain reading: RUVs with all genes as controls and the groups as replicate sets matches the true-factor ceiling in every core cell; under the second reading of rule 4 it is the only practical design not failing anything (arm 11 needs no k estimate). Today's design and RUVr fail when the factor is tied to the groups (C5: FDP 0.30 and 0.52). No RUV is conservative but finds almost nothing when the unwanted factor is strong. C6 (library replicates) is reported only: every RUV arm has a mean FDP of 0.27 to 0.35 there, the package's replicate handling (parked row 23, repaired by KIT-128, not in this tree).

### S4, lipidomics and metabolomics (complete; C4 void by the scoring unit)

C4 is void by rule 2: the true-factor arm's mean FDP in the global null is 0.130 (lower bound 0.109, limit 0.075), while its per-test null rate is 0.047. **This is an artefact of the scoring unit, not of the ceiling arm:** S4 scores one FDP per run over both assays (and all contrasts in C7), while each assay controls its FDR on its own; under the global null a run's FDP is 1 when either assay calls anything, about 1 - (1 - a)^2. The implied per-assay any-call rate (0.067) matches the single-family arms of S1 (0.061) and S3 (0.075). S4 has therefore no global-null check; a per-assay score needs per-assay rows the outputs do not hold (a new run, outside the budget). The pooling may also shape rule 3's tail in low-power S4 cells; its limit is pooled the same way.

| Arm | Design | Verdict | Failing core cells | Share FDP > 0.10 | Power C0 |
|---|---|---|---|---|---|
| 1 | no correction | fail | rule 4 in every cell with a factor | 0.066 | 0.003 |
| 2 | true drift subtracted, factor in the model (ceiling) | pass | none | 0.143 | 0.50 |
| 3 | today (sweep, `ruvIII_C_Varying`, test) | fail | rules 1-3 in every judged cell | 0.915 | 0.42 |
| 7 | all features as controls, k = E, W in the model | fail | C1(3), C5(3) | 0.206 | 0.43 |
| 8 | internal standards as controls, W in the model | fail | C5(3), C7(3) | 0.180 | 0.37 |
| 11 | arm 7 with k = 1 | fail | C3(3), C5(3) | 0.116 | 0.07 |
| Qs | pooled QC as the replicate set, subtract | fail (inconclusive on complete features) | C1(3), C5(3) | 0.140 | 0.26 |
| Qm | the same, W in the model | fail | C0, C1, C5, C7 (3) | 0.217 | 0.42 |
| I, L | ITSD normalisation alone; QC-RLSC | fail | C3(5) | 0.066, 0.071 | 0.002, 0.030 |
| 14 | `sva` (complete features only; arm 7 on the complete subset: share 0.291, power 0.69) | fail | C2, C5, C7 | 0.310 | 0.70 |
| 13 | `RUVinv` with the standards | not applicable | refused in 1,000 of 1,000 runs in every core cell (scored as arm 1 by the refusal rule; counted among failures by the aggregation) | | |

Plain reading: no design passes. Today's design is the worst measured (more than 90% of runs above 10% FDP). The closest are all features or the standards as controls with the factor in the model, and pooled QC as the replicate set: each misses rule 3 (the FDP tail) in two cells. No correction is calibrated but has almost no power under drift.

### S2, DIA_LIMPA (complete; 8 core cells, 112,000 rows; scored 2026-11-18 03:30)

| Arm | Design | Verdict (literal) | Without rule 4 in C4 | Failing core cells | Share FDP > 0.10 | Power C0 |
|---|---|---|---|---|---|---|
| 1 | no RUV | fail | fail | rules 3, 4 in every cell with a factor; C5 rules 1-3 (null rate 0.116, mean FDP 0.19) | 0.153 | 0.56 |
| 2 | true factor in the model (ceiling) | inconclusive | pass | none | 0.013 | 0.91 |
| 3 | the route's correction mode (sweep, `computeRuvCovariates`, test on the corrected EList) | fail | fail | rules 1-3 in every cell | 0.643 | 0.91 |
| 4 | the route's covariates mode, today's default (the sweep's controls and k, W in the model) | fail | fail | rules 1-3 in every cell but C3 | 0.570 | 0.90 |
| 5, 6 | data-driven controls, covariates | fail | fail | rules 1-3 in most cells | 0.260, 0.328 | 0.92, 0.88 |
| 7 | all proteins as controls, k = E, covariates | inconclusive | **pass** | none | 0.018 | 0.91 |
| 11 | arm 7 with k = 1 | inconclusive | **pass** | none | 0.018 | 0.91 |
| 8 | oracle null controls, covariates (diagnostic) | inconclusive | inconclusive | none | 0.032 | 0.91 |
| 9, 10 | degraded supplied controls; oracle controls subtracted | fail | fail | rule 3 / rules 1-3 | 0.160, 0.275 | 0.90, 0.93 |
| 12, 13, 14 | `RUV4`, `RUVinv`, `sva` | fail | fail | C3(5), C6(1-3) and others | 0.15-0.16 | 0.82-0.88 |
| 15 | technical replicates as the structure, subtract (C6 only) | inconclusive | inconclusive | none | 0.075 | n/a |

Plain reading: **one design is valid in every core cell.** All proteins as controls with the factors as covariates (arm 7, and arm 11 with k = 1, which needs no estimate of k) matches the true-factor ceiling cell by cell (null rate 0.050 to 0.052, mean FDP 0.044 to 0.055, power within 0.01 of the ceiling) and is the only practical design with no failure under the second reading of rule 4; literally it is inconclusive only through the C4 rule-4 point that no design, the ceiling included, can meet. Today's default (covariates mode with the sweep's controls) and correction mode are miscalibrated in nearly every cell (more than half the runs above 10% FDP). Unlike the other omics, **no RUV also fails here**: the unwanted factor inflates the tail (15% of runs above 10% FDP) and, when tied to the groups, the null rate (0.116). The aggregation's branch: "neither passes" literally; "eligible: 7 > 11" without rule 4 in C4.

**Engine (M-11):** S2 ran on `fastlimpadpc`, which the route selects when installed. The two engines agreed to 0.0013 log2 on the S2 pilot's probe; no separate exactness record of the two engines exists (KIT-053's record covers the route, not engine identity). The DIA_LIMPA row is therefore claimed for `fastlimpadpc`; for `limpa` it is expected to transfer (RUV acts on the quantified matrix, which the engines agree on to 0.0013 log2) but is not measured. Arm 6's 1,396 refusals ("RUV-III needs more than k control proteins") are designed refusals carrying the no-RUV result.

## What the rule gives, literally (ticket W3, "Recommendation per omic")

| Omic | Literal branch | With rule 4 undecidable in C4 |
|---|---|---|
| Classic proteomics | Without technical replicates: RUV off by default (no RUV design passes; arm 1 passes rules 1, 2, 6, 7). With them: arms 15/16 eligible only under A14 (see S1) | Same (arm 7 still fails C1 rule 3) |
| Transcriptomics | RUV off by default | Arms 7/11 inconclusive (no failure): selectable with a note under KIT-129's verdict-to-option rule |
| Lipidomics, metabolomics | RUV off by default | Same (C4 is void anyway) |
| DIA_LIMPA | Neither passes (no RUV fails rules 1-3 too): reported as such | Arm 7 (then 11) eligible: all proteins as controls, factors as covariates |

## For the gate (to be completed with S2)

- Both readings of rule 4 in C4, side by side (above).
- The verdict-to-option rule proposed in KIT-129 revision 1.1, for confirmation: a passing design is selectable; the default per omic is the user's choice among passing ones; a failing design is not offered except subtract-then-test as the warned legacy switch; an inconclusive or not-run design is selectable with a note; diagnostic arms (oracle controls) are never options.
- Read literally, that rule leaves most omics with no RUV option at all, and "RUV off by default" means accepting the power losses above (transcriptomics and lipid/metab near zero under strong unwanted variation). This is the user's call; the measures are here to overrule the rule's literal reading if the user chooses to.
- The `sva` choice (k from `sva::num.sv` in the arms that take k = E; `sva` is not a dependency; IR-M6).
- Stress cells: conditions of use (minimum group size, high missingness), added when scored.

## Independent audit of the scoring and this draft (2026-11-17, read-only)

`aggregate.R` at `0bb67eb2` (byte-identical to the script that scored the core phases) implements the W3 rules faithfully; every figure recomputed from the raw outputs matched (S1 C0, C1, C4, C5, C6; S3 C4, C5; S4 C4), and no per-omic verdict changes. Applied to this draft: S4's C4 void explained by the scoring unit (F1); arms 15/16 "pass" in the one cell they run, A14 (F2); S1's literal branch is RUV off by default (F3); "inconclusive" = as good as the true-factor arm (F4); arm 14 against complete-subset values (F5); C1 described correctly (F6); the void named by rule 2 (F7); the failing-cells column's omissions noted (F8); arm 13 in S4 is not applicable (F10); wording (F13). Recorded, not changed: rule 1 uses run-level normal bounds rather than a pooled binomial (A1; chosen because features within a run are dependent), and rule 4's FDP-tail reduction uses normal bounds on paired 0/1 differences (anti-conservative with few discordant pairs; no RUV design's rule-4 pass depends on it) (F9); shares carry exact bounds but no MCSE column (F11); refused rows with no measures would drop silently (latent, not in any core cell) (F12). Limit: FDP was not recomputed from truth tables (the simulators write per-run results only).

## New arms and the technical-duplicate re-measurement, scored 2026-11-18

Two runs on GitHub-hosted runners, scored with `aggregate.R` (unchanged, sha256 `c5c6275f`). Scoring only; nothing was rerun. Tables, sidecars, job tables and the comparison files: `scored/newarms-gha/` and `scored/r0a-gha/`. The per-run tables are not committed (5 to 12 MB per simulator): their sha256 are in each folder's `not-committed.sha256` and in the sidecars, and they can be downloaded again from the runs' artifacts (digests in `jobs.tsv`).

| Run | Content | Jobs | Attempts | Rows |
|---|---|---|---|---|
| 44200924581 ("newarms"; workflow commit `8cdcc5d6`, package tree `0b6fcb6d`, simulator files identical to `35e5a25a`) | Arms 1, 2, 22, 23, 24 (S4: 1, 2, 22, 23), C0 to C7, 1,000 runs | 32 simulator jobs and the plan job, all succeeded | 1 for every job; no rerun | 1,000 per arm, cell and subset; every row `ok`; no refusal, no crash |
| 44203261310 ("r0a"; workflow commit `3db3a2ed`, package tree `ace1f4de`, kit-128) | C6 only: S3 arms 1, 2, 7; S4 arms 1, 2, 7, 15, 16; 1,000 runs | 2 simulator jobs and the plan job, all succeeded | 1; no rerun | 1,000 per arm and subset; every row `ok` |

Arm 22 is arm 8's design (known nulls; in S4 the internal standards) with k from the package's separation rule on those controls. Arm 23 is arm 5's design (first-pass nulls) with k = E, the estimated count. Arm 24 (S1 to S3) is arm 5's design with the separation rule computed on arm 5's own control set. Cells: C0 centre (6 per group, 20% effects, one factor); C1 12 per group; C2 10% effects; C3 no unwanted factor; C4 no true effects; C5 factor tied to the groups (0.3); C6 technical duplicates; C7 three groups.

### The equivalence gate (declared before the launch)

Arms 1 and 2 of the runner against the laptop's confirm rows for the same seeds, under `packaging/verify/compare-results.R` (relative 1e-6, absolute floor 1e-9, discrete exact). Files: `scored/newarms-gha/compare/`.

| Simulator | Rows compared | Discrete differences | Rows outside the rule | Gate |
|---|---|---|---|---|
| S1 | 32,000 (8 cells, 1,000 runs, 2 arms, 2 subsets) | 0 | 0 (largest relative difference 2.2e-10) | pass, all cells |
| S2 | 16,000 | 0 | 0 (3.1e-8) | pass, all cells |
| S3 | 16,000 | 0 | 8 rows (9 values): `bias` in C0 (1 row), C1 (2), C5 (1), C6 (1), C7 (3), and `rmse` in one of the C7 rows | **fails in C0, C1, C5, C6, C7**; passes in C2, C3, C4 |
| S4 | 32,000 | 0 | 0 (identical) | pass, all cells |

In S3 the null rates, the false discovery proportion, power, k, the feature counts and the number of calls are identical in all 16,000 rows. The values outside the rule are the bias of the estimated effects (relative differences 1.0e-6 to 6.7e-5, absolute 3.6e-9 to 1.9e-7) and one rmse (relative 6.4e-6); they are listed in `compare/gate-rows-outside.tsv`. No rule of the ticket reads bias or rmse. The gate as declared is nevertheless not met in those five cells, so **arms 22, 23 and 24 are not scored in S3 C0, C1, C5, C6 and C7**; S3 was aggregated from C2, C3 and C4 only. Whether a difference confined to bias and rmse may be accepted is the user's decision; the data are downloaded and one `aggregate.R` call scores the five cells if it is. The 20-run proof had not shown this (bias and rmse within 5.6e-8 there).

Arms 1 and 2 scored from the runner's rows give the same figures as the committed laptop tables in every cell scored (86 rows, every measure identical).

### The new arms beside the arm each differs from by one choice

Primary subset (`full`). Arms 8 and 5 are read from the laptop tables (`scored/s*-core`); arms 22 to 24 from the runner. Standard errors are Monte Carlo errors over 1,000 runs per arm; the arms share their data sets, so a difference between two arms is known more precisely than the two errors suggest, and no paired error was computed. Mean k has no standard error in the tables. The limit of rule 3 is the ceiling arm's share plus 0.05. Every table: `scored/newarms-gha/new-arms-beside-reference.tsv`.

**S1, classic proteomics**

| Cell | Arm | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Cell verdict (rules not passed) |
|---|---|---|---|---|---|---|---|---|
| C0 | 8 (laptop) | 1.00 | 0.050 (0.0003) | 0.051 (0.0007) | 0.025 [0.016, 0.037] | 0.063 | 0.912 (0.0008) | pass |
| C0 | **22** | 1.50 | 0.051 (0.0003) | 0.052 (0.0007) | 0.019 [0.011, 0.030] | 0.063 | 0.903 (0.0009) | pass |
| C0 | 5 (laptop) | 1.41 | 0.051 (0.0003) | 0.052 (0.0006) | 0.014 [0.008, 0.023] | 0.063 | 0.908 (0.0009) | pass |
| C0 | **23** | 1.00 | 0.050 (0.0003) | 0.051 (0.0006) | 0.014 [0.008, 0.023] | 0.063 | 0.914 (0.0008) | pass |
| C0 | **24** | 1.59 | 0.051 (0.0003) | 0.052 (0.0006) | 0.012 [0.006, 0.021] | 0.063 | 0.906 (0.0009) | pass |
| C1 | 8 (laptop) | 1.00 | 0.050 (0.0003) | 0.051 (0.0007) | 0.020 [0.012, 0.031] | 0.057 | 0.977 (0.0004) | pass |
| C1 | **22** | 2.96 | 0.051 (0.0003) | 0.052 (0.0007) | 0.026 [0.017, 0.038] | 0.057 | 0.974 (0.0004) | pass |
| C1 | 5 (laptop) | 2.96 | 0.051 (0.0003) | 0.052 (0.0006) | 0.014 [0.008, 0.023] | 0.057 | 0.975 (0.0004) | pass |
| C1 | **23** | 1.00 | 0.050 (0.0003) | 0.050 (0.0006) | 0.010 [0.005, 0.018] | 0.057 | 0.977 (0.0004) | pass |
| C1 | **24** | 2.97 | 0.051 (0.0003) | 0.052 (0.0006) | 0.014 [0.008, 0.023] | 0.057 | 0.975 (0.0004) | pass |
| C2 | 8 (laptop) | 1.00 | 0.050 (0.0003) | 0.051 (0.0009) | 0.049 [0.036, 0.064] | 0.106 | 0.868 (0.0013) | pass |
| C2 | **22** | 1.36 | 0.050 (0.0003) | 0.052 (0.0009) | 0.058 [0.044, 0.074] | 0.106 | 0.858 (0.0015) | pass |
| C2 | 5 (laptop) | 1.28 | 0.050 (0.0003) | 0.053 (0.0009) | 0.054 [0.041, 0.070] | 0.106 | 0.865 (0.0014) | pass |
| C2 | **23** | 1.00 | 0.050 (0.0003) | 0.051 (0.0009) | 0.051 [0.038, 0.067] | 0.106 | 0.871 (0.0013) | pass |
| C2 | **24** | 1.47 | 0.050 (0.0003) | 0.053 (0.0009) | 0.058 [0.044, 0.074] | 0.106 | 0.862 (0.0014) | pass |
| C3 | 8 (laptop) | 0.12 | 0.050 (0.0003) | 0.052 (0.0007) | 0.016 [0.009, 0.026] | 0.064 | 0.922 (0.0008) | pass |
| C3 | **22** | 1.95 | 0.052 (0.0003) | 0.057 (0.0007) | 0.033 [0.023, 0.046] | 0.064 | 0.900 (0.0010) | fail: rule 5 fail |
| C3 | 5 (laptop) | 1.99 | 0.052 (0.0003) | 0.056 (0.0007) | 0.032 [0.022, 0.045] | 0.064 | 0.906 (0.0009) | pass |
| C3 | **23** | 0.12 | 0.050 (0.0003) | 0.052 (0.0007) | 0.014 [0.008, 0.023] | 0.064 | 0.922 (0.0008) | pass |
| C3 | **24** | 1.99 | 0.052 (0.0003) | 0.056 (0.0007) | 0.032 [0.022, 0.045] | 0.064 | 0.906 (0.0009) | pass |
| C4 | 8 (laptop) | 1.00 | 0.050 (0.0003) | 0.062 (0.0076) | 0.062 [0.048, 0.079] | 0.111 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |
| C4 | **22** | 1.26 | 0.050 (0.0003) | 0.062 (0.0076) | 0.062 [0.048, 0.079] | 0.111 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |
| C4 | 5 (laptop) | 1.22 | 0.050 (0.0003) | 0.062 (0.0076) | 0.062 [0.048, 0.079] | 0.111 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |
| C4 | **23** | 1.00 | 0.050 (0.0003) | 0.061 (0.0076) | 0.061 [0.047, 0.078] | 0.111 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |
| C4 | **24** | 1.35 | 0.050 (0.0003) | 0.066 (0.0079) | 0.066 [0.051, 0.083] | 0.111 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |
| C5 | 8 (laptop) | 1.00 | 0.051 (0.0004) | 0.053 (0.0008) | 0.029 [0.020, 0.041] | 0.064 | 0.898 (0.0009) | pass |
| C5 | **22** | 1.48 | 0.051 (0.0003) | 0.053 (0.0007) | 0.031 [0.021, 0.044] | 0.064 | 0.889 (0.0010) | pass |
| C5 | 5 (laptop) | 1.57 | 0.065 (0.0004) | 0.078 (0.0010) | 0.205 [0.180, 0.231] | 0.064 | 0.900 (0.0009) | fail: rule 2 fail; rule 3 fail |
| C5 | **23** | 1.00 | 0.065 (0.0004) | 0.077 (0.0010) | 0.204 [0.179, 0.230] | 0.064 | 0.908 (0.0008) | fail: rule 2 fail; rule 3 fail |
| C5 | **24** | 1.57 | 0.065 (0.0004) | 0.078 (0.0010) | 0.198 [0.174, 0.224] | 0.064 | 0.900 (0.0009) | fail: rule 2 fail; rule 3 fail |
| C6 | 8 (laptop) | 1.01 | 0.048 (0.0003) | 0.049 (0.0007) | 0.021 [0.013, 0.032] | 0.059 | 0.894 (0.0009) | pass |
| C6 | **22** | 2.89 | 0.054 (0.0003) | 0.060 (0.0007) | 0.050 [0.037, 0.065] | 0.059 | 0.885 (0.0011) | inconclusive: rule 3 inconclusive |
| C6 | 5 (laptop) | 2.98 | 0.048 (0.0003) | 0.048 (0.0006) | 0.007 [0.003, 0.014] | 0.059 | 0.880 (0.0010) | pass |
| C6 | **23** | 1.01 | 0.049 (0.0003) | 0.050 (0.0006) | 0.013 [0.007, 0.022] | 0.059 | 0.896 (0.0009) | pass |
| C6 | **24** | 2.98 | 0.048 (0.0003) | 0.048 (0.0006) | 0.007 [0.003, 0.014] | 0.059 | 0.880 (0.0010) | pass |
| C7 | 8 (laptop) | 1.00 | 0.050 (0.0003) | 0.051 (0.0009) | 0.059 [0.045, 0.075] | 0.100 | 0.889 (0.0012) | pass |
| C7 | **22** | 2.97 | 0.051 (0.0003) | 0.053 (0.0009) | 0.061 [0.047, 0.078] | 0.100 | 0.866 (0.0015) | pass |
| C7 | 5 (laptop) | 2.92 | 0.051 (0.0003) | 0.055 (0.0009) | 0.065 [0.051, 0.082] | 0.100 | 0.876 (0.0014) | pass |
| C7 | **23** | 1.00 | 0.050 (0.0003) | 0.051 (0.0009) | 0.050 [0.037, 0.065] | 0.100 | 0.892 (0.0012) | pass |
| C7 | **24** | 2.98 | 0.052 (0.0003) | 0.055 (0.0009) | 0.066 [0.051, 0.083] | 0.100 | 0.875 (0.0014) | pass |

**S2, DIA_LIMPA**

| Cell | Arm | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Cell verdict (rules not passed) |
|---|---|---|---|---|---|---|---|---|
| C0 | 8 (laptop) | 1.00 | 0.050 (0.0003) | 0.051 (0.0007) | 0.022 [0.014, 0.033] | 0.051 | 0.905 (0.0012) | pass |
| C0 | **22** | 1.61 | 0.051 (0.0003) | 0.053 (0.0007) | 0.025 [0.016, 0.037] | 0.051 | 0.903 (0.0012) | pass |
| C0 | 5 (laptop) | 1.78 | 0.069 (0.0009) | 0.094 (0.0022) | 0.286 [0.258, 0.315] | 0.051 | 0.916 (0.0008) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C0 | **23** | 1.00 | 0.068 (0.0009) | 0.090 (0.0022) | 0.261 [0.234, 0.289] | 0.051 | 0.916 (0.0008) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C0 | **24** | 1.68 | 0.069 (0.0009) | 0.093 (0.0022) | 0.282 [0.254, 0.311] | 0.051 | 0.916 (0.0008) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C1 | 8 (laptop) | 1.00 | 0.049 (0.0003) | 0.048 (0.0006) | 0.016 [0.009, 0.026] | 0.052 | 0.991 (0.0002) | pass |
| C1 | **22** | 3.00 | 0.051 (0.0003) | 0.052 (0.0006) | 0.019 [0.011, 0.030] | 0.052 | 0.989 (0.0003) | pass |
| C1 | 5 (laptop) | 3.00 | 0.062 (0.0006) | 0.076 (0.0015) | 0.182 [0.159, 0.207] | 0.052 | 0.991 (0.0002) | fail: rule 1 inconclusive; rule 2 inconclusive; rule 3 fail |
| C1 | **23** | 1.00 | 0.060 (0.0006) | 0.071 (0.0014) | 0.151 [0.129, 0.175] | 0.052 | 0.991 (0.0002) | fail: rule 3 fail |
| C1 | **24** | 2.88 | 0.062 (0.0006) | 0.076 (0.0014) | 0.180 [0.157, 0.205] | 0.052 | 0.991 (0.0002) | fail: rule 1 inconclusive; rule 2 inconclusive; rule 3 fail |
| C2 | 8 (laptop) | 1.00 | 0.050 (0.0003) | 0.052 (0.0009) | 0.058 [0.044, 0.074] | 0.071 | 0.865 (0.0021) | inconclusive: rule 3 inconclusive |
| C2 | **22** | 1.54 | 0.051 (0.0003) | 0.054 (0.0010) | 0.061 [0.047, 0.078] | 0.071 | 0.863 (0.0021) | inconclusive: rule 3 inconclusive |
| C2 | 5 (laptop) | 1.69 | 0.070 (0.0009) | 0.111 (0.0030) | 0.381 [0.351, 0.412] | 0.071 | 0.880 (0.0014) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C2 | **23** | 1.00 | 0.069 (0.0009) | 0.108 (0.0030) | 0.359 [0.329, 0.390] | 0.071 | 0.880 (0.0014) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C2 | **24** | 1.67 | 0.070 (0.0009) | 0.111 (0.0030) | 0.383 [0.353, 0.414] | 0.071 | 0.880 (0.0014) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C3 | 8 (laptop) | 0.17 | 0.049 (0.0002) | 0.049 (0.0005) | 0.002 [0.000, 0.007] | 0.053 | 0.932 (0.0006) | pass |
| C3 | **22** | 2.00 | 0.052 (0.0003) | 0.053 (0.0006) | 0.008 [0.003, 0.016] | 0.053 | 0.925 (0.0007) | pass |
| C3 | 5 (laptop) | 2.00 | 0.052 (0.0003) | 0.054 (0.0006) | 0.012 [0.006, 0.021] | 0.053 | 0.931 (0.0006) | pass |
| C3 | **23** | 0.17 | 0.050 (0.0003) | 0.049 (0.0005) | 0.003 [0.001, 0.009] | 0.053 | 0.932 (0.0006) | pass |
| C3 | **24** | 2.00 | 0.052 (0.0003) | 0.054 (0.0006) | 0.012 [0.006, 0.021] | 0.053 | 0.931 (0.0006) | pass |
| C4 | 8 (laptop) | 1.00 | 0.050 (0.0003) | 0.054 (0.0072) | 0.054 [0.041, 0.070] | 0.101 | n/a | inconclusive: rule 4 inconclusive |
| C4 | **22** | 1.50 | 0.051 (0.0003) | 0.059 (0.0075) | 0.059 [0.045, 0.075] | 0.101 | n/a | inconclusive: rule 4 inconclusive; rule 7 inconclusive |
| C4 | 5 (laptop) | 1.68 | 0.070 (0.0008) | 0.216 (0.0130) | 0.216 [0.191, 0.243] | 0.101 | n/a | fail: rule 1 fail; rule 2 fail; rule 3 fail; rule 4 fail; rule 7 fail |
| C4 | **23** | 1.00 | 0.069 (0.0008) | 0.209 (0.0129) | 0.209 [0.184, 0.236] | 0.101 | n/a | fail: rule 1 fail; rule 2 fail; rule 3 fail; rule 4 fail; rule 7 fail |
| C4 | **24** | 1.64 | 0.070 (0.0008) | 0.218 (0.0131) | 0.218 [0.193, 0.245] | 0.101 | n/a | fail: rule 1 fail; rule 2 fail; rule 3 fail; rule 4 fail; rule 7 fail |
| C5 | 8 (laptop) | 1.00 | 0.051 (0.0004) | 0.054 (0.0008) | 0.049 [0.036, 0.064] | 0.053 | 0.878 (0.0024) | inconclusive: rule 3 inconclusive |
| C5 | **22** | 1.65 | 0.052 (0.0004) | 0.056 (0.0008) | 0.052 [0.039, 0.068] | 0.053 | 0.875 (0.0024) | inconclusive: rule 3 inconclusive |
| C5 | 5 (laptop) | 1.83 | 0.092 (0.0015) | 0.149 (0.0036) | 0.528 [0.497, 0.559] | 0.053 | 0.902 (0.0013) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C5 | **23** | 1.00 | 0.091 (0.0015) | 0.146 (0.0036) | 0.513 [0.482, 0.544] | 0.053 | 0.902 (0.0013) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C5 | **24** | 1.73 | 0.092 (0.0015) | 0.149 (0.0036) | 0.527 [0.496, 0.558] | 0.053 | 0.902 (0.0013) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C6 | 8 (laptop) | 1.01 | 0.052 (0.0003) | 0.056 (0.0008) | 0.040 [0.029, 0.054] | 0.071 | 0.758 (0.0016) | pass |
| C6 | **22** | 2.99 | 0.049 (0.0003) | 0.050 (0.0007) | 0.026 [0.017, 0.038] | 0.071 | 0.731 (0.0017) | pass |
| C6 | 5 (laptop) | 3.00 | 0.060 (0.0006) | 0.073 (0.0014) | 0.174 [0.151, 0.199] | 0.071 | 0.762 (0.0013) | fail: rule 2 inconclusive; rule 3 fail |
| C6 | **23** | 1.01 | 0.063 (0.0006) | 0.081 (0.0015) | 0.232 [0.206, 0.259] | 0.071 | 0.777 (0.0013) | fail: rule 1 inconclusive; rule 2 fail; rule 3 fail |
| C6 | **24** | 2.96 | 0.060 (0.0006) | 0.073 (0.0015) | 0.174 [0.151, 0.199] | 0.071 | 0.763 (0.0013) | fail: rule 2 inconclusive; rule 3 fail |
| C7 | 8 (laptop) | 1.00 | 0.049 (0.0003) | 0.050 (0.0006) | 0.012 [0.006, 0.021] | 0.053 | 0.920 (0.0010) | pass |
| C7 | **22** | 2.99 | 0.052 (0.0003) | 0.054 (0.0006) | 0.016 [0.009, 0.026] | 0.053 | 0.913 (0.0010) | pass |
| C7 | 5 (laptop) | 3.00 | 0.070 (0.0008) | 0.093 (0.0020) | 0.299 [0.271, 0.328] | 0.053 | 0.926 (0.0007) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C7 | **23** | 1.00 | 0.067 (0.0008) | 0.087 (0.0019) | 0.267 [0.240, 0.296] | 0.053 | 0.927 (0.0007) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C7 | **24** | 2.90 | 0.070 (0.0008) | 0.093 (0.0019) | 0.299 [0.271, 0.328] | 0.053 | 0.926 (0.0007) | fail: rule 1 fail; rule 2 fail; rule 3 fail |

**S3, transcriptomics (C2, C3, C4 only: the other cells did not pass the gate)**

| Cell | Arm | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Cell verdict (rules not passed) |
|---|---|---|---|---|---|---|---|---|
| C2 | 8 (laptop) | 1.00 | 0.056 (0.0002) | 0.056 (0.0005) | 0.021 [0.013, 0.032] | 0.051 | 0.696 (0.0008) | pass |
| C2 | **22** | 2.99 | 0.054 (0.0002) | 0.053 (0.0006) | 0.022 [0.014, 0.033] | 0.051 | 0.533 (0.0046) | pass |
| C2 | 5 (laptop) | 2.97 | 0.055 (0.0001) | 0.054 (0.0004) | 0.001 [0.000, 0.006] | 0.051 | 0.637 (0.0030) | pass |
| C2 | **23** | 1.00 | 0.056 (0.0001) | 0.055 (0.0004) | 0.000 [0.000, 0.004] | 0.051 | 0.709 (0.0008) | pass |
| C2 | **24** | 2.96 | 0.055 (0.0001) | 0.053 (0.0004) | 0.001 [0.000, 0.006] | 0.051 | 0.638 (0.0030) | pass |
| C3 | 8 (laptop) | 1.35 | 0.052 (0.0001) | 0.043 (0.0003) | 0.000 [0.000, 0.004] | 0.050 | 0.696 (0.0045) | fail: rule 5 fail |
| C3 | **22** | 1.99 | 0.052 (0.0001) | 0.043 (0.0003) | 0.000 [0.000, 0.004] | 0.050 | 0.662 (0.0039) | fail: rule 5 fail |
| C3 | 5 (laptop) | 2.01 | 0.053 (0.0001) | 0.044 (0.0002) | 0.000 [0.000, 0.004] | 0.050 | 0.770 (0.0005) | pass |
| C3 | **23** | 1.35 | 0.053 (0.0001) | 0.044 (0.0002) | 0.000 [0.000, 0.004] | 0.050 | 0.774 (0.0006) | pass |
| C3 | **24** | 2.08 | 0.053 (0.0001) | 0.044 (0.0002) | 0.000 [0.000, 0.004] | 0.050 | 0.770 (0.0006) | pass |
| C4 | 8 (laptop) | 1.00 | 0.056 (0.0002) | 0.085 (0.0088) | 0.085 [0.068, 0.104] | 0.125 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |
| C4 | **22** | 2.00 | 0.055 (0.0002) | 0.078 (0.0085) | 0.078 [0.062, 0.096] | 0.125 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |
| C4 | 5 (laptop) | 1.49 | 0.055 (0.0001) | 0.077 (0.0084) | 0.077 [0.061, 0.095] | 0.125 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |
| C4 | **23** | 1.00 | 0.056 (0.0001) | 0.082 (0.0087) | 0.082 [0.066, 0.101] | 0.125 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |
| C4 | **24** | 1.63 | 0.055 (0.0001) | 0.073 (0.0082) | 0.073 [0.058, 0.091] | 0.125 | n/a | fail: rule 2 inconclusive; rule 4 fail; rule 7 inconclusive |

**S4, lipidomics and metabolomics (no arm 24; C4 and C6 void on this tree)**

| Cell | Arm | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Cell verdict (rules not passed) |
|---|---|---|---|---|---|---|---|---|
| C0 | 8 (laptop) | 1.95 | 0.044 (0.0005) | 0.050 (0.0016) | 0.146 [0.125, 0.169] | 0.184 | 0.366 (0.0035) | pass |
| C0 | **22** | 4.07 | 0.050 (0.0006) | 0.077 (0.0028) | 0.317 [0.288, 0.347] | 0.184 | 0.262 (0.0038) | fail: rule 2 inconclusive; rule 3 fail |
| C0 | 5 (laptop) | 4.24 | 0.060 (0.0006) | 0.097 (0.0018) | 0.452 [0.421, 0.483] | 0.184 | 0.502 (0.0027) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C0 | **23** | 1.95 | 0.044 (0.0005) | 0.052 (0.0014) | 0.135 [0.114, 0.158] | 0.184 | 0.437 (0.0033) | pass |
| C1 | 8 (laptop) | 2.00 | 0.042 (0.0005) | 0.041 (0.0012) | 0.066 [0.051, 0.083] | 0.125 | 0.685 (0.0023) | pass |
| C1 | **22** | 4.06 | 0.043 (0.0005) | 0.047 (0.0013) | 0.100 [0.082, 0.120] | 0.125 | 0.644 (0.0025) | pass |
| C1 | 5 (laptop) | 4.15 | 0.048 (0.0005) | 0.054 (0.0012) | 0.105 [0.087, 0.126] | 0.125 | 0.712 (0.0020) | inconclusive: rule 3 inconclusive |
| C1 | **23** | 2.00 | 0.041 (0.0004) | 0.041 (0.0010) | 0.038 [0.027, 0.052] | 0.125 | 0.715 (0.0021) | pass |
| C2 | 8 (laptop) | 1.96 | 0.044 (0.0005) | 0.064 (0.0030) | 0.281 [0.253, 0.310] | 0.305 | 0.266 (0.0041) | inconclusive: rule 3 inconclusive |
| C2 | **22** | 4.25 | 0.050 (0.0006) | 0.118 (0.0054) | 0.415 [0.384, 0.446] | 0.305 | 0.170 (0.0038) | fail: rule 2 fail; rule 3 fail |
| C2 | 5 (laptop) | 5.35 | 0.072 (0.0007) | 0.188 (0.0036) | 0.752 [0.724, 0.778] | 0.305 | 0.449 (0.0038) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C2 | **23** | 1.96 | 0.044 (0.0004) | 0.062 (0.0025) | 0.276 [0.248, 0.305] | 0.305 | 0.342 (0.0040) | pass |
| C3 | 8 (laptop) | 0.14 | 0.049 (0.0004) | 0.058 (0.0013) | 0.149 [0.127, 0.173] | 0.187 | 0.573 (0.0027) | pass |
| C3 | **22** | 2.03 | 0.053 (0.0005) | 0.071 (0.0015) | 0.243 [0.217, 0.271] | 0.187 | 0.491 (0.0037) | fail: rule 3 fail; rule 5 fail |
| C3 | 5 (laptop) | 2.06 | 0.059 (0.0005) | 0.083 (0.0014) | 0.334 [0.305, 0.364] | 0.187 | 0.620 (0.0025) | fail: rule 1 inconclusive; rule 2 fail; rule 3 fail |
| C3 | **23** | 0.14 | 0.049 (0.0004) | 0.059 (0.0012) | 0.153 [0.131, 0.177] | 0.187 | 0.582 (0.0026) | pass |
| C4 (void) | 8 (laptop) | 1.95 | 0.044 (0.0005) | 0.141 (0.0110) | 0.141 [0.120, 0.164] | 0.180 | n/a | fail: rule 2 fail; rule 4 fail; rule 7 fail |
| C4 (void) | **22** | 2.73 | 0.035 (0.0006) | 0.125 (0.0105) | 0.125 [0.105, 0.147] | 0.180 | n/a | fail: rule 2 fail; rule 4 fail; rule 7 fail |
| C4 (void) | 5 (laptop) | 5.83 | 0.079 (0.0012) | 0.737 (0.0139) | 0.737 [0.709, 0.764] | 0.180 | n/a | fail: rule 1 fail; rule 2 fail; rule 3 fail; rule 4 fail; rule 7 fail |
| C4 (void) | **23** | 1.95 | 0.044 (0.0004) | 0.140 (0.0110) | 0.140 [0.119, 0.163] | 0.180 | n/a | fail: rule 2 fail; rule 4 fail; rule 7 fail |
| C5 | 8 (laptop) | 1.95 | 0.049 (0.0008) | 0.068 (0.0026) | 0.245 [0.219, 0.273] | 0.194 | 0.276 (0.0034) | fail: rule 3 fail |
| C5 | **22** | 4.10 | 0.051 (0.0006) | 0.088 (0.0035) | 0.346 [0.317, 0.376] | 0.194 | 0.185 (0.0034) | fail: rule 2 fail; rule 3 fail |
| C5 | 5 (laptop) | 4.24 | 0.111 (0.0010) | 0.239 (0.0032) | 0.926 [0.908, 0.941] | 0.194 | 0.471 (0.0026) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C5 | **23** | 1.95 | 0.093 (0.0009) | 0.174 (0.0031) | 0.758 [0.730, 0.784] | 0.194 | 0.432 (0.0030) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C6 (void) | 8 (laptop) | 1.99 | 0.113 (0.0008) | 0.223 (0.0025) | 0.941 [0.925, 0.955] | 1.041 | 0.582 (0.0028) | fail: rule 1 fail; rule 2 fail |
| C6 (void) | **22** | 4.08 | 0.109 (0.0008) | 0.221 (0.0028) | 0.910 [0.891, 0.927] | 1.041 | 0.511 (0.0033) | fail: rule 1 fail; rule 2 fail |
| C6 (void) | 5 (laptop) | 4.19 | 0.180 (0.0010) | 0.384 (0.0024) | 1.000 [0.996, 1.000] | 1.041 | 0.695 (0.0022) | fail: rule 1 fail; rule 2 fail |
| C6 (void) | **23** | 1.99 | 0.126 (0.0007) | 0.252 (0.0023) | 0.978 [0.967, 0.986] | 1.041 | 0.635 (0.0024) | fail: rule 1 fail; rule 2 fail |
| C7 | 8 (laptop) | 1.99 | 0.046 (0.0004) | 0.066 (0.0018) | 0.193 [0.169, 0.219] | 0.163 | 0.325 (0.0029) | fail: rule 3 fail |
| C7 | **22** | 4.28 | 0.048 (0.0004) | 0.081 (0.0022) | 0.318 [0.289, 0.348] | 0.163 | 0.232 (0.0029) | fail: rule 2 fail; rule 3 fail |
| C7 | 5 (laptop) | 4.75 | 0.062 (0.0004) | 0.115 (0.0019) | 0.576 [0.545, 0.607] | 0.163 | 0.441 (0.0024) | fail: rule 1 fail; rule 2 fail; rule 3 fail |
| C7 | **23** | 1.99 | 0.051 (0.0004) | 0.072 (0.0017) | 0.223 [0.198, 0.250] | 0.163 | 0.405 (0.0025) | fail: rule 2 inconclusive; rule 3 fail |

Plain reading.

- **Supplied controls, separation rule against estimated count (22 against 8).** The separation rule picks more factors than there are: 1.3 to 3.0 where the truth is 1 (the count: 1.00), about 2 where there is no factor (the count: 0.12 to 0.17, and 1.35 in transcriptomics), 4.1 to 4.3 in lipid/metab where the truth is 2 (the count: 1.95). Classic proteomics: a small cost in both directions, mean FDP higher by 0 to 0.005 in seven cells and by 0.011 with technical duplicates (C6, 0.060 against 0.049, rule 3 inconclusive), power lower by 0.003 to 0.023; with no factor (C3) the power loss against no RUV is 0.023 and rule 5 fails where arm 8 passed. DIA_LIMPA: practically no difference (mean FDP within 0.005, power within 0.007, except C6: FDP 0.050 against 0.056, power 0.731 against 0.758); the same verdicts. Transcriptomics, three cells only: calibration unchanged, power clearly lower at 10% effects (C2: 0.533 against 0.696) and lower with no factor (0.662 against 0.696). Lipid/metab: worse on both counts in every judged cell (centre cell: mean FDP 0.077 against 0.050, share above 0.10 0.317 against 0.146, power 0.262 against 0.366), and with no factor it fails rules 3 and 5 where arm 8 passed. In short: the separation rule never helps the supplied-controls design; it costs little in the two proteomics routes, power in transcriptomics, and calibration and power in lipid/metab.
- **First-pass nulls, estimated count against separation rule (23 against 5).** Classic proteomics: no material change (mean FDP within 0.004); the failure with the factor tied to the groups stays (C5: mean FDP 0.077 against 0.078, share above 0.10 0.204 against 0.205); power is higher by 0.002 to 0.016. DIA_LIMPA: **not rescued.** Mean FDP 0.071 to 0.146 in the cells with a factor and true effects (arm 5: 0.073 to 0.149), 15% to 51% of data sets above 0.10, rule 3 failed in each of them as for arm 5; with technical duplicates it is slightly worse (0.081 against 0.073). The poor figure there comes from the control set, not from k. Lipid/metab: **largely rescued where the factor is not tied to the groups**: mean FDP 0.052 against 0.097 (C0), 0.041 against 0.054 (C1), 0.062 against 0.188 (C2), 0.059 against 0.083 (C3); these four cells now pass every rule; no refusals (arm 5 refused 28 runs in C2). Power is lower than arm 5's (0.437 against 0.502 in C0), but arm 5's power was measured at about twice the nominal FDP; the ceiling is 0.498. The same holds against arm 8 with the same k: the first-pass set gives more power than the standards (0.437 against 0.366) at the same mean FDP (0.052 against 0.050). It still fails with the factor tied to the groups (C5: mean FDP 0.174, type I 0.093; arm 5 0.239 and 0.111) and misses rule 3 with three groups (C7: share 0.223 against a limit of 0.163; mean FDP 0.072). Transcriptomics: in the three cells scored the power loss is gone (10% effects, C2: 0.709 against 0.637, the ceiling's 0.709; no factor: 0.774 against 0.770) with calibration unchanged. The centre cell, where arm 5's loss was reported (0.60 against 0.77), is among the cells not scored.
- **Arm 24 against arm 5.** Not materially different in any cell scored: mean FDP within 0.004, power within 0.003, the same rule results in every cell. The two differ only in the k chosen in some runs (mean k within 0.2).
- Lipid/metab C4 (void by the scoring unit) and C6 (void on this tree, before the technical-replicate repair) judge nothing; their rows are in the table for completeness.

### The technical-duplicate cell after the repair (run r0a, kit-128 tree)

"Before" is the committed laptop scoring on tree `0b6fcb6d` (`scored/s3-core`, `scored/s4-core`); "after" is the runner on kit-128's tree with the simulators of `sim-k129-r0a`. No laptop rows exist on that tree, so this run has **no equivalence comparison**; its only link to the laptop is the tiny local proof recorded in KIT-129 (R0). In S3 the two libraries of a unit are now summed before the analysis (8 samples instead of 16), arms 15 and 16 cannot be formed, and arm 2's factor is the mean of a unit's two library scores, an approximate ceiling.

| Cell, subset | Arm (tree) | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Cell verdict (rules not passed) |
|---|---|---|---|---|---|---|---|---|
| S3 C6, full (void) | 1, before | 0.00 | 0.056 (0.0003) | 0.070 (0.0005) | 0.025 [0.016, 0.037] | 1.050 | 0.299 (0.0008) | fail: rule 4 fail |
| S3 C6, full | 1, after | 0.00 | 0.015 (0.0003) | 0.000 (0.0001) | 0.000 [0.000, 0.004] | 0.150 | 0.002 (0.0005) | fail: rule 4 fail |
| S3 C6, full (void) | 2, before | 1.00 | 0.179 (0.0003) | 0.275 (0.0006) | 1.000 [0.996, 1.000] | 1.050 | 0.840 (0.0004) | fail: rule 1 fail; rule 2 fail |
| S3 C6, full | 2, after | 1.00 | 0.055 (0.0005) | 0.047 (0.0013) | 0.100 [0.082, 0.120] | 0.150 | 0.372 (0.0034) | pass |
| S3 C6, full (void) | 7, before | 1.00 | 0.179 (0.0003) | 0.275 (0.0006) | 1.000 [0.996, 1.000] | 1.050 | 0.840 (0.0004) | fail: rule 1 fail; rule 2 fail |
| S3 C6, full | 7, after | 1.07 | 0.055 (0.0005) | 0.050 (0.0012) | 0.104 [0.086, 0.125] | 0.150 | 0.421 (0.0027) | pass |
| S3 C6, full (void) | 15, before | 1.00 | 0.176 (0.0002) | 0.271 (0.0005) | 1.000 [0.996, 1.000] | 1.050 | 0.834 (0.0004) | fail: rule 1 fail; rule 2 fail |
| S3 C6, full (void) | 16, before | 1.00 | 0.178 (0.0003) | 0.274 (0.0006) | 1.000 [0.996, 1.000] | 1.050 | 0.840 (0.0004) | fail: rule 1 fail; rule 2 fail |
| S4 C6, full (void) | 1, before | n/a | 0.027 (0.0004) | 0.035 (0.0034) | 0.129 [0.109, 0.151] | 1.041 | 0.061 (0.0015) | fail: rule 4 fail |
| S4 C6, full | 1, after | n/a | 0.018 (0.0003) | 0.019 (0.0031) | 0.060 [0.046, 0.077] | 0.223 | 0.030 (0.0011) | fail: rule 4 fail |
| S4 C6, full (void) | 2, before | n/a | 0.135 (0.0007) | 0.271 (0.0023) | 0.991 [0.983, 0.996] | 1.041 | 0.697 (0.0023) | fail: rule 1 fail; rule 2 fail |
| S4 C6, full | 2, after | n/a | 0.047 (0.0005) | 0.055 (0.0015) | 0.173 [0.150, 0.198] | 0.223 | 0.410 (0.0031) | pass |
| S4 C6, full (void) | 7, before | 1.99 | 0.135 (0.0008) | 0.272 (0.0024) | 0.987 [0.978, 0.993] | 1.041 | 0.635 (0.0025) | fail: rule 1 fail; rule 2 fail |
| S4 C6, full | 7, after | 1.99 | 0.049 (0.0005) | 0.058 (0.0017) | 0.202 [0.178, 0.228] | 0.223 | 0.354 (0.0030) | inconclusive: rule 3 inconclusive |
| S4 C6, full (void) | 15, before | 1.99 | 0.158 (0.0009) | 0.336 (0.0025) | 0.998 [0.993, 1.000] | 1.041 | 0.681 (0.0025) | fail: rule 1 fail; rule 2 fail |
| S4 C6, full | 15, after | 1.99 | 0.044 (0.0005) | 0.051 (0.0018) | 0.169 [0.146, 0.194] | 0.223 | 0.293 (0.0035) | pass |
| S4 C6, full (void) | 16, before | 1.99 | 0.134 (0.0008) | 0.271 (0.0024) | 0.989 [0.980, 0.994] | 1.041 | 0.631 (0.0025) | fail: rule 1 fail; rule 2 fail |
| S4 C6, full | 16, after | 1.99 | 0.049 (0.0005) | 0.058 (0.0017) | 0.187 [0.163, 0.213] | 0.223 | 0.356 (0.0030) | pass |
| S4 C6, complete (void) | 1, before | n/a | 0.030 (0.0005) | 0.044 (0.0029) | 0.199 [0.175, 0.225] | 1.044 | 0.198 (0.0035) | fail: rule 4 fail |
| S4 C6, complete | 1, after | n/a | 0.018 (0.0004) | 0.025 (0.0029) | 0.089 [0.072, 0.108] | 0.361 | 0.098 (0.0029) | fail: rule 4 fail |
| S4 C6, complete (void) | 2, before | n/a | 0.164 (0.0011) | 0.337 (0.0030) | 0.994 [0.987, 0.998] | 1.044 | 0.864 (0.0024) | fail: rule 1 fail; rule 2 fail |
| S4 C6, complete | 2, after | n/a | 0.054 (0.0007) | 0.074 (0.0021) | 0.311 [0.282, 0.341] | 0.361 | 0.579 (0.0043) | inconclusive: rule 2 inconclusive |
| S4 C6, complete (void) | 7, before | 1.99 | 0.172 (0.0012) | 0.351 (0.0030) | 0.994 [0.987, 0.998] | 1.044 | 0.853 (0.0025) | fail: rule 1 fail; rule 2 fail |
| S4 C6, complete | 7, after | 1.99 | 0.057 (0.0007) | 0.079 (0.0023) | 0.322 [0.293, 0.352] | 0.361 | 0.544 (0.0045) | inconclusive: rule 2 inconclusive |
| S4 C6, complete (void) | 15, before | 1.99 | 0.176 (0.0012) | 0.361 (0.0029) | 0.998 [0.993, 1.000] | 1.044 | 0.857 (0.0025) | fail: rule 1 fail; rule 2 fail |
| S4 C6, complete | 15, after | 1.99 | 0.053 (0.0007) | 0.072 (0.0023) | 0.288 [0.260, 0.317] | 0.361 | 0.521 (0.0045) | inconclusive: rule 2 inconclusive |
| S4 C6, complete (void) | 16, before | 1.99 | 0.171 (0.0012) | 0.347 (0.0030) | 0.993 [0.986, 0.997] | 1.044 | 0.848 (0.0025) | fail: rule 1 fail; rule 2 fail |
| S4 C6, complete | 16, after | 1.99 | 0.057 (0.0007) | 0.079 (0.0023) | 0.331 [0.302, 0.361] | 0.361 | 0.547 (0.0044) | inconclusive: rule 2 inconclusive; rule 3 inconclusive |

Plain reading: **the cell is no longer void in either simulator.** The ceiling arm now meets rules 1 and 2 on the full set (transcriptomics: type I 0.055, upper bound 0.056, at 0.01 0.0118, mean FDP 0.047; lipid/metab: 0.047, 0.0099, 0.055). On the complete-feature subset of lipid/metab the ceiling's mean FDP is 0.074 with an upper bound of 0.078, so rule 2 is inconclusive there: flagged, not void. All features as controls (arm 7) is calibrated in both (mean FDP 0.050 and 0.058); in lipid/metab its share above 0.10 (0.202, upper bound 0.228) sits at the limit of rule 3 (0.223), inconclusive. Technical replicates as the structure in lipid/metab: subtracting (arm 15) mean FDP 0.051 and power 0.293; in the model (arm 16) 0.058 and 0.356, the same as arm 7 (0.354); both pass in this one cell. Power is far below the old figures because those counted each duplicate as an independent sample. In transcriptomics arm 7's power (0.421) is above the ceiling arm's (0.372): the ceiling there is approximate, as noted. The cell stays "reported, not core" in both simulators, so these figures judge no design.
