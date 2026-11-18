# KIT-124: RUV in every omic, confirming studies (core cells complete, 2026-11-18; extra-hard cells complete and scored 2026-11-18)

**Corrected 2026-11-18 (final scoring of the extra-hard cells).** The title read "stress cells running" and the paragraph below says they "are still running": every extra-hard cell that is defined has now run and is scored in the last section, "Extra-hard scenarios, final, and conditions of use". They change no verdict. One qualifier to the per-simulator readings below, which are core-cell statements: no correction is conservative (S1, S3) or calibrated (S4) in the core cells, where the factor is tied to the groups at 0.3 at most; at 0.6 its mean false discovery proportion is 0.593 (S1), 0.450 (S3) and 0.642 (S4), and 0.533 in S2 (last section).

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
- Stress cells: conditions of use (minimum group size, high missingness), added when scored. (Added 2026-11-18: the last section.)

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
| S3 | 16,000 | 0 | 8 rows (9 values): `bias` in C0 (1 row), C1 (2), C5 (1), C6 (1), C7 (3), and `rmse` in one of the C7 rows | as declared: **fails in C0, C1, C5, C6, C7**, passes in C2, C3, C4; as read under M19: pass, with the exception stated below |
| S4 | 32,000 | 0 | 0 (identical) | pass, all cells |

In S3 the null rates, the false discovery proportion, power, k, the feature counts and the number of calls are identical in all 16,000 rows. The values outside the rule are the bias of the estimated effects (relative differences 1.0e-6 to 6.7e-5, absolute 3.6e-9 to 1.9e-7) and one rmse (relative 6.4e-6); they are listed in `compare/gate-rows-outside.tsv`. No rule of the ticket reads bias or rmse. The gate as first declared (every value within the rule) is therefore not met in those five cells, and at the first scoring arms 22, 23 and 24 were held back in S3 C0, C1, C5, C6 and C7 (commit `38047a96`). **The user accepted the exception on 2026-11-18 (decision M19, `design-review-20261118.md` section 6): the gate is read on the values the rules use, identical in all 16,000 rows; the 8 rows outside the tolerance on bias and rmse are accepted and reported.** The five cells were then scored from the same downloaded data, and every figure of the cells scored before is unchanged. `compare/gate-summary.tsv` carries both readings per cell. The 20-run proof had not shown the difference (bias and rmse within 5.6e-8 there).

| Cell | Run | Arm | Measure | Laptop | Runner | Relative difference |
|---|---|---|---|---|---|---|
| C0 | 564 | 2 | bias | 0.0103163943 | 0.0103163618 | 3.2e-6 |
| C1 | 151 | 1 | bias | -0.0009244668 | -0.0009244741 | 7.9e-6 |
| C1 | 977 | 1 | bias | -0.0031892499 | -0.0031892567 | 2.1e-6 |
| C5 | 614 | 1 | bias | -0.0035918684 | -0.0035918720 | 1.0e-6 |
| C6 | 323 | 1 | bias | 0.0061616440 | 0.0061616503 | 1.0e-6 |
| C7 | 218 | 1 | bias | -0.0026642314 | -0.0026640533 | 6.7e-5 |
| C7 | 218 | 1 | rmse | 0.5565873986 | 0.5565838332 | 6.4e-6 |
| C7 | 218 | 2 | bias | -0.0072835920 | -0.0072834015 | 2.6e-5 |
| C7 | 823 | 2 | bias | -0.0044673143 | -0.0044674502 | 3.0e-5 |

Arms 1 and 2 scored from the runner's rows give the same figures as the committed laptop tables in every cell (96 rows, every measure identical).

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

**S3, transcriptomics (all eight cells; C6 void on this tree)**

| Cell | Arm | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Cell verdict (rules not passed) |
|---|---|---|---|---|---|---|---|---|
| C0 | 8 (laptop) | 1.00 | 0.056 (0.0002) | 0.048 (0.0003) | 0.002 [0.000, 0.007] | 0.050 | 0.763 (0.0006) | pass |
| C0 | **22** | 3.00 | 0.054 (0.0002) | 0.047 (0.0005) | 0.005 [0.002, 0.012] | 0.050 | 0.612 (0.0045) | pass |
| C0 | 5 (laptop) | 2.92 | 0.055 (0.0001) | 0.047 (0.0005) | 0.004 [0.001, 0.010] | 0.050 | 0.601 (0.0066) | pass |
| C0 | **23** | 1.00 | 0.056 (0.0001) | 0.048 (0.0002) | 0.000 [0.000, 0.004] | 0.050 | 0.774 (0.0005) | pass |
| C0 | **24** | 2.96 | 0.055 (0.0001) | 0.047 (0.0004) | 0.003 [0.001, 0.009] | 0.050 | 0.600 (0.0065) | pass |
| C1 | 8 (laptop) | 1.00 | 0.056 (0.0002) | 0.049 (0.0003) | 0.002 [0.000, 0.007] | 0.050 | 0.954 (0.0002) | pass |
| C1 | **22** | 3.04 | 0.055 (0.0002) | 0.048 (0.0003) | 0.000 [0.000, 0.004] | 0.050 | 0.934 (0.0007) | pass |
| C1 | 5 (laptop) | 3.04 | 0.057 (0.0001) | 0.049 (0.0005) | 0.008 [0.003, 0.016] | 0.050 | 0.816 (0.0083) | pass |
| C1 | **23** | 1.00 | 0.057 (0.0001) | 0.050 (0.0002) | 0.000 [0.000, 0.004] | 0.050 | 0.956 (0.0002) | pass |
| C1 | **24** | 3.00 | 0.057 (0.0001) | 0.050 (0.0005) | 0.009 [0.004, 0.017] | 0.050 | 0.827 (0.0079) | pass |
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
| C5 | 8 (laptop) | 1.00 | 0.055 (0.0002) | 0.048 (0.0004) | 0.003 [0.001, 0.009] | 0.050 | 0.721 (0.0007) | pass |
| C5 | **22** | 2.99 | 0.054 (0.0002) | 0.046 (0.0005) | 0.008 [0.003, 0.016] | 0.050 | 0.580 (0.0045) | pass |
| C5 | 5 (laptop) | 2.93 | 0.058 (0.0002) | 0.064 (0.0044) | 0.125 [0.105, 0.147] | 0.050 | 0.107 (0.0059) | fail: rule 3 fail |
| C5 | **23** | 1.00 | 0.071 (0.0002) | 0.071 (0.0004) | 0.022 [0.014, 0.033] | 0.050 | 0.741 (0.0006) | fail: rule 1 fail |
| C5 | **24** | 2.54 | 0.061 (0.0003) | 0.064 (0.0035) | 0.100 [0.082, 0.120] | 0.050 | 0.255 (0.0097) | fail: rule 3 fail |
| C6 (void) | 8 (laptop) | 1.00 | 0.178 (0.0004) | 0.274 (0.0008) | 1.000 [0.996, 1.000] | 1.050 | 0.832 (0.0005) | fail: rule 1 fail; rule 2 fail |
| C6 (void) | **22** | 3.01 | 0.208 (0.0007) | 0.351 (0.0013) | 0.998 [0.993, 1.000] | 1.050 | 0.719 (0.0048) | fail: rule 1 fail; rule 2 fail |
| C6 (void) | 5 (laptop) | 2.99 | 0.214 (0.0008) | 0.350 (0.0015) | 0.993 [0.986, 0.997] | 1.050 | 0.779 (0.0048) | fail: rule 1 fail; rule 2 fail |
| C6 (void) | **23** | 1.00 | 0.179 (0.0003) | 0.275 (0.0006) | 1.000 [0.996, 1.000] | 1.050 | 0.840 (0.0005) | fail: rule 1 fail; rule 2 fail |
| C6 (void) | **24** | 3.06 | 0.215 (0.0008) | 0.353 (0.0015) | 0.991 [0.983, 0.996] | 1.050 | 0.778 (0.0048) | fail: rule 1 fail; rule 2 fail |
| C7 | 8 (laptop) | 1.00 | 0.057 (0.0002) | 0.051 (0.0004) | 0.003 [0.001, 0.009] | 0.050 | 0.752 (0.0006) | pass |
| C7 | **22** | 3.02 | 0.056 (0.0002) | 0.050 (0.0003) | 0.000 [0.000, 0.004] | 0.050 | 0.666 (0.0029) | pass |
| C7 | 5 (laptop) | 2.98 | 0.055 (0.0001) | 0.057 (0.0065) | 0.086 [0.069, 0.105] | 0.050 | 0.023 (0.0039) | fail: rule 3 fail; rule 4 fail |
| C7 | **23** | 1.00 | 0.058 (0.0001) | 0.051 (0.0003) | 0.000 [0.000, 0.004] | 0.050 | 0.763 (0.0006) | pass |
| C7 | **24** | 2.99 | 0.055 (0.0001) | 0.064 (0.0068) | 0.091 [0.074, 0.111] | 0.050 | 0.049 (0.0057) | fail: rule 2 inconclusive; rule 3 fail; rule 4 fail |

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

- **Supplied controls, separation rule against estimated count (22 against 8).** The separation rule picks more factors than there are: 1.3 to 3.0 where the truth is 1 (the count: 1.00), about 2 where there is no factor (the count: 0.12 to 0.17, and 1.35 in transcriptomics), 4.1 to 4.3 in lipid/metab where the truth is 2 (the count: 1.95). Classic proteomics: a small cost in both directions, mean FDP higher by 0 to 0.005 in seven cells and by 0.011 with technical duplicates (C6, 0.060 against 0.049, rule 3 inconclusive), power lower by 0.003 to 0.023; with no factor (C3) the power loss against no RUV is 0.023 and rule 5 fails where arm 8 passed. DIA_LIMPA: practically no difference (mean FDP within 0.005, power within 0.007, except C6: FDP 0.050 against 0.056, power 0.731 against 0.758); the same verdicts. Transcriptomics: k is 3.0 where the truth is 1; calibration is unchanged (mean FDP 0.043 to 0.053 against 0.043 to 0.056), power is lower in every cell with true effects: centre 0.612 (SE 0.0045) against 0.763 (0.0006), 12 per group 0.934 against 0.954, 10% effects 0.533 against 0.696, factor tied to the groups 0.580 against 0.721, three groups 0.666 against 0.752, no factor 0.662 against 0.696 (both fail rule 5 there). Lipid/metab: worse on both counts in every judged cell (centre cell: mean FDP 0.077 against 0.050, share above 0.10 0.317 against 0.146, power 0.262 against 0.366), and with no factor it fails rules 3 and 5 where arm 8 passed. In short: the separation rule never helps the supplied-controls design; it costs little in the two proteomics routes, 0.02 to 0.16 of power in transcriptomics, and calibration and power in lipid/metab.
- **First-pass nulls, estimated count against separation rule (23 against 5).** Classic proteomics: no material change (mean FDP within 0.004); the failure with the factor tied to the groups stays (C5: mean FDP 0.077 against 0.078, share above 0.10 0.204 against 0.205); power is higher by 0.002 to 0.016. DIA_LIMPA: **not rescued.** Mean FDP 0.071 to 0.146 in the cells with a factor and true effects (arm 5: 0.073 to 0.149), 15% to 51% of data sets above 0.10, rule 3 failed in each of them as for arm 5; with technical duplicates it is slightly worse (0.081 against 0.073). The poor figure there comes from the control set, not from k. Lipid/metab: **largely rescued where the factor is not tied to the groups**: mean FDP 0.052 against 0.097 (C0), 0.041 against 0.054 (C1), 0.062 against 0.188 (C2), 0.059 against 0.083 (C3); these four cells now pass every rule; no refusals (arm 5 refused 28 runs in C2). Power is lower than arm 5's (0.437 against 0.502 in C0), but arm 5's power was measured at about twice the nominal FDP; the ceiling is 0.498. The same holds against arm 8 with the same k: the first-pass set gives more power than the standards (0.437 against 0.366) at the same mean FDP (0.052 against 0.050). It still fails with the factor tied to the groups (C5: mean FDP 0.174, type I 0.093; arm 5 0.239 and 0.111) and misses rule 3 with three groups (C7: share 0.223 against a limit of 0.163; mean FDP 0.072). Transcriptomics: **the power loss is gone.** Centre cell: power 0.774 (SE 0.0005) against arm 5's 0.601 (0.0066), equal to the ceiling's 0.774, at the same mean FDP (0.048 against 0.047). 12 per group 0.956 against 0.816; 10% effects 0.709 against 0.637; three groups 0.763 against 0.023 (arm 5 there also fails rule 3); no factor 0.774 against 0.770; each at the ceiling's power or within 0.008 of it, with mean FDP 0.044 to 0.055. With the factor tied to the groups (C5) power is 0.741 against arm 5's 0.107, but the design is liberal there: type I 0.071 (SE 0.0002; arm 5 0.058) and mean FDP 0.071 (0.0004; the ceiling 0.047), so rule 1 fails. Its share of data sets above 0.10 is 0.022 there, against arm 5's 0.125. Arm 5 failed that cell too (rule 3), in a different way: few calls, and a tail of bad data sets.
- **Arm 24 against arm 5.** Not materially different in classic proteomics and DIA_LIMPA: mean FDP within 0.004, power within 0.003, the same rule results in every cell, mean k within 0.2. In transcriptomics the same holds in six cells (power within 0.011, which is inside the two arms' errors where it is largest, 12 per group: 0.827, SE 0.0079, against 0.816, SE 0.0083), but not with the factor tied to the groups (C5: power 0.255, SE 0.0097, against 0.107, SE 0.0059; mean k 2.54 against 2.93; both fail rule 3) nor quite with three groups (0.049, SE 0.0057, against 0.023, SE 0.0039; both fail). Both keep arm 5's power loss in the centre cell (0.600 against 0.601). Arm 24 is arm 5 with a slightly different k in some runs; it removes none of arm 5's failures.
- Lipid/metab C4 (void by the scoring unit) and C6 of lipid/metab and of transcriptomics (void on this tree, before the technical-replicate repair) judge nothing; their rows are in the table for completeness.

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

## Extra-hard scenarios, final, and conditions of use (scored 2026-11-18)

Every extra-hard cell that is defined has now run, 500 runs per arm, and is scored. Scoring and reading only: nothing was rerun and no code changed. One scoring per simulator, covering its main cells and all its extra-hard cells: `scored/final-20261118/s1` to `s4` (`aggregate.R` unchanged, sha256 `c5c6275f`, run at commit `19409bb6`). The extra-hard cells state conditions of use; they never change a design's verdict (ticket W3), and the verdict columns of the four final tables equal the committed `scored/s*-core` tables.

### Where each cell ran

| Simulator | Extra-hard cells on the laptop (tree `0b6fcb6d`) | On hosted runners (run 44201244793, label `hard`, package tree `0b6fcb6d`) | Not run |
|---|---|---|---|
| S1 classic proteomics | S1 to S8 | S9 to S14 | none |
| S2 DIA_LIMPA | none | S1 to S14 | none |
| S3 transcriptomics | S1 to S5, S7 to S13 | S14 | S6: not defined for this simulator (revision 1.3: no 50%-missing cell) |
| S4 lipidomics and metabolomics | S1 to S14 | none | none |

- 55 extra-hard cells, each complete: runs 1 to 500 for every arm the cell has, every row with a status. No cell exists on both platforms and no cell mixes rows of the two. Per cell and platform (arms, rows, ok, refused, crashed): `scored/final-20261118/inventory.tsv`.
- The laptop's S1 phase ended at its ceiling (exit 124, 2026-11-18 09:38) while working on S9; the simulator writes a cell only when it is complete, so no partial S9 rows exist and the runner's complete S9 is used. S3's give-back phase ended at its ceiling (08:30) before S14; the runner's S14 is used. S2's extra-hard cells never ran on the laptop.
- The hosted-runner run: 21 simulator jobs and the plan job, all succeeded, attempt 1 for each, no rerun (`scored/final-20261118/jobs.tsv`, with the artifact digests). The 21 artifacts hold result tables, per-cell restart files, the simulators' meta and log files and the per-job library check and console output; no input data. The per-run tables are not committed (sha256 in `not-committed.sha256` and in the sidecars).
- Refusals and crashes in the option-set designs: none, except **S2 S5 (DIA_LIMPA, no missing values): 252 of 500 runs crash in every arm**, no correction and the ceiling included (the route's `dpcDE` failure on complete data, parked row 25); the figures of that cell are over the 248 runs that completed and rule 6 fails for every arm. In S4 arm 5 refuses 412, 414 and 426 of 500 runs in the 3-per-group cells with a factor (S1, S11, S12), 20 in S13 and 2 in S6, and arm 3 refuses 12 in S6 (refused runs are scored as no correction, by the refusal rule). The other refusals are in arms 4, 6 and 13 (`k124-cells.tsv`).

### Equivalence between the platforms

No extra-hard cell exists on both platforms, so arms 1 and 2 could not be compared for these cells. **The runner's extra-hard cells rest on the proof of 2026-11-18** (`gha-proof-20261118/README.md`): arms 1, 2, 7, 8 in core cells within the fixed rule, and all arms in core cells and one extra-hard cell per simulator within the rule except arm 3; decision M19 fixes how the gate is read. Two checks that could be made here were made: the main-cell rows of the four final tables are byte-identical to `scored/s1-core` to `s4-core` (688 rows), and the extra-hard rows already in `scored/s34-stress-partial` (S3: 9 cells, 118 rows; S4: 11 cells, 365 rows) are byte-identical in the final tables, which read the same files.

**† Legacy (arm 3) on the runner.** Arm 3's per-run values differ between machines (the proof and the bounded BLAS test; accepted and stated by the user). Every arm-3 figure marked † below (S1 S9 to S14, S2 S1 to S14, S3 S14) is the runner's and would not be reproduced run for run on the laptop; in the proof's 20 runs the per-cell means stayed close (mean false discovery proportion 0.0369 on the laptop against 0.0351 on the runner in S1 C6; 0.4134 against 0.4140 in S4 C0).

### The figures

Full set of features, 500 runs per arm (S2 S5: 248). Type I is the share of true nulls with p below 0.05; standard errors are Monte Carlo errors over runs; the share of data sets above 0.10 carries exact one-sided 97.5% bounds; the rule 3 limit is the ceiling's share plus 0.05. "Rules not passed" lists rules 1 to 6 as `aggregate.R` judged them (rule 1 also reads the 0.01 level, in `k124-cells.tsv`); "none" fails rule 4 by construction wherever a factor exists (no gain over itself). A **void** scenario is one where the ceiling (true factor known) itself fails rule 1 or 2: it judges no design. No technical-replicate design (arms 15, 16) is defined in any extra-hard cell. Arms share their data sets, so a difference between two arms is known more precisely than the two errors suggest; no paired error was computed. Source: `scored/final-20261118/option-arms-extra-hard.tsv` (with arm 11, k = 1, and arm 17, the known batch in the model).

**S1, classic proteomics (S1 to S8 laptop, S9 to S14 runner)**

| Scenario | Design | Scored runs | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Rules not passed |
|---|---|---|---|---|---|---|---|---|---|
| S1 3 per group | 2 ceiling | 500 | 1.00 | 0.050 (0.0004) | 0.050 (0.0010) | 0.018 [0.008, 0.034] | 0.068 | 0.661 (0.0033) | within every limit |
| S1 3 per group | **7 recommended** | 500 | 1.00 | 0.051 (0.0004) | 0.051 (0.0011) | 0.022 [0.011, 0.039] | 0.068 | 0.660 (0.0033) | within every limit |
| S1 3 per group | 8 known controls | 500 | 1.00 | 0.050 (0.0005) | 0.051 (0.0011) | 0.028 [0.015, 0.047] | 0.068 | 0.644 (0.0036) | within every limit |
| S1 3 per group | 5 first-pass nulls | 500 | 1.35 | 0.052 (0.0005) | 0.055 (0.0012) | 0.042 [0.026, 0.063] | 0.068 | 0.621 (0.0046) | within every limit |
| S1 3 per group | 3 legacy | 500 | 1.35 | 0.188 (0.0052) | 0.327 (0.0098) | 0.928 [0.902, 0.949] | 0.068 | 0.853 (0.0054) | 1 fail; 2 fail; 3 fail |
| S1 3 per group | 1 none | 500 | 0.00 | 0.010 (0.0002) | 0.002 (0.0020) | 0.002 [0.000, 0.011] | 0.068 | 0.000 (0.0000) | 4 fail |
| S2 40% of features changed | 2 ceiling | 500 | 1.00 | 0.050 (0.0005) | 0.051 (0.0007) | 0.000 [0.000, 0.007] | 0.050 | 0.948 (0.0006) | within every limit |
| S2 40% of features changed | **7 recommended** | 500 | 1.00 | 0.054 (0.0006) | 0.057 (0.0008) | 0.020 [0.010, 0.036] | 0.050 | 0.948 (0.0006) | within every limit |
| S2 40% of features changed | 8 known controls | 500 | 1.00 | 0.050 (0.0005) | 0.051 (0.0007) | 0.006 [0.001, 0.017] | 0.050 | 0.947 (0.0006) | within every limit |
| S2 40% of features changed | 5 first-pass nulls | 500 | 1.51 | 0.054 (0.0005) | 0.056 (0.0008) | 0.016 [0.007, 0.031] | 0.050 | 0.944 (0.0007) | within every limit |
| S2 40% of features changed | 3 legacy | 500 | 1.51 | 0.072 (0.0008) | 0.082 (0.0013) | 0.244 [0.207, 0.284] | 0.050 | 0.947 (0.0035) | 1 fail; 2 fail; 3 fail |
| S2 40% of features changed | 1 none | 500 | 0.00 | 0.014 (0.0003) | 0.007 (0.0003) | 0.000 [0.000, 0.007] | 0.050 | 0.623 (0.0018) | 4 fail |
| S3 factor tied to groups, 0.6 | 2 ceiling | 500 | 1.00 | 0.050 (0.0004) | 0.050 (0.0009) | 0.008 [0.002, 0.020] | 0.058 | 0.822 (0.0017) | within every limit |
| S3 factor tied to groups, 0.6 | **7 recommended** | 500 | 1.00 | 0.057 (0.0005) | 0.062 (0.0012) | 0.076 [0.054, 0.103] | 0.058 | 0.840 (0.0017) | 3 inconclusive |
| S3 factor tied to groups, 0.6 | 8 known controls | 500 | 1.00 | 0.055 (0.0006) | 0.058 (0.0013) | 0.070 [0.049, 0.096] | 0.058 | 0.834 (0.0018) | 3 inconclusive |
| S3 factor tied to groups, 0.6 | 5 first-pass nulls | 500 | 1.75 | 0.146 (0.0013) | 0.273 (0.0035) | 0.990 [0.977, 0.997] | 0.058 | 0.877 (0.0014) | 1 fail; 2 fail; 3 fail |
| S3 factor tied to groups, 0.6 | 3 legacy | 500 | 1.75 | 0.437 (0.0056) | 0.684 (0.0065) | 0.980 [0.964, 0.990] | 0.058 | 0.829 (0.0065) | 1 fail; 2 fail; 3 fail |
| S3 factor tied to groups, 0.6 | 1 none | 500 | 0.00 | 0.258 (0.0008) | 0.593 (0.0025) | 1.000 [0.993, 1.000] | 0.058 | 0.762 (0.0019) | 1 fail; 2 fail; 3 fail; 4 fail |
| S4 three factors | 2 ceiling | 500 | 3.00 | 0.049 (0.0004) | 0.048 (0.0009) | 0.010 [0.003, 0.023] | 0.060 | 0.863 (0.0015) | within every limit |
| S4 three factors | **7 recommended** | 500 | 3.00 | 0.054 (0.0005) | 0.057 (0.0011) | 0.054 [0.036, 0.078] | 0.060 | 0.862 (0.0016) | 3 inconclusive |
| S4 three factors | 8 known controls | 500 | 3.00 | 0.048 (0.0004) | 0.048 (0.0009) | 0.018 [0.008, 0.034] | 0.060 | 0.852 (0.0016) | within every limit |
| S4 three factors | 5 first-pass nulls | 500 | 1.20 | 0.009 (0.0002) | 0.003 (0.0011) | 0.004 [0.000, 0.014] | 0.060 | 0.071 (0.0063) | within every limit |
| S4 three factors | 3 legacy | 500 | 1.20 | 0.021 (0.0008) | 0.016 (0.0017) | 0.022 [0.011, 0.039] | 0.060 | 0.143 (0.0089) | within every limit |
| S4 three factors | 1 none | 500 | 0.00 | 0.005 (0.0001) | 0.000 (0.0000) | 0.000 [0.000, 0.007] | 0.060 | 0.000 (0.0000) | 4 fail |
| S5 no missing values | 2 ceiling | 500 | 1.00 | 0.050 (0.0003) | 0.050 (0.0007) | 0.002 [0.000, 0.011] | 0.052 | 0.949 (0.0007) | within every limit |
| S5 no missing values | **7 recommended** | 500 | 1.01 | 0.050 (0.0003) | 0.051 (0.0007) | 0.002 [0.000, 0.011] | 0.052 | 0.948 (0.0007) | within every limit |
| S5 no missing values | 8 known controls | 500 | 1.01 | 0.049 (0.0004) | 0.049 (0.0008) | 0.008 [0.002, 0.020] | 0.052 | 0.948 (0.0008) | within every limit |
| S5 no missing values | 5 first-pass nulls | 500 | 1.30 | 0.050 (0.0003) | 0.050 (0.0008) | 0.000 [0.000, 0.007] | 0.052 | 0.949 (0.0008) | within every limit |
| S5 no missing values | 3 legacy | 500 | 1.30 | 0.066 (0.0006) | 0.082 (0.0013) | 0.250 [0.213, 0.290] | 0.052 | 0.957 (0.0008) | 1 fail; 2 fail; 3 fail |
| S5 no missing values | 1 none | 500 | 0.00 | 0.010 (0.0002) | 0.005 (0.0003) | 0.000 [0.000, 0.007] | 0.052 | 0.603 (0.0019) | 4 fail |
| S6 50% missing | 2 ceiling | 500 | 1.00 | 0.049 (0.0005) | 0.051 (0.0012) | 0.054 [0.036, 0.078] | 0.104 | 0.887 (0.0017) | within every limit |
| S6 50% missing | **7 recommended** | 500 | 1.01 | 0.056 (0.0007) | 0.062 (0.0016) | 0.124 [0.096, 0.156] | 0.104 | 0.886 (0.0017) | 3 inconclusive |
| S6 50% missing | 8 known controls | 500 | 1.01 | 0.049 (0.0005) | 0.050 (0.0012) | 0.044 [0.028, 0.066] | 0.104 | 0.883 (0.0018) | within every limit |
| S6 50% missing | 5 first-pass nulls | 500 | 1.35 | 0.049 (0.0005) | 0.053 (0.0012) | 0.048 [0.031, 0.071] | 0.104 | 0.876 (0.0019) | within every limit |
| S6 50% missing | 3 legacy | 500 | 1.35 | 0.071 (0.0010) | 0.099 (0.0023) | 0.462 [0.418, 0.507] | 0.104 | 0.858 (0.0065) | 1 fail; 2 fail; 3 fail |
| S6 50% missing | 1 none | 500 | 0.00 | 0.016 (0.0003) | 0.008 (0.0006) | 0.000 [0.000, 0.007] | 0.104 | 0.424 (0.0035) | 4 fail |
| S7 heavy tails | 2 ceiling | 500 | 1.00 | 0.042 (0.0004) | 0.035 (0.0008) | 0.002 [0.000, 0.011] | 0.052 | 0.922 (0.0010) | within every limit |
| S7 heavy tails | **7 recommended** | 500 | 1.01 | 0.046 (0.0005) | 0.041 (0.0010) | 0.018 [0.008, 0.034] | 0.052 | 0.919 (0.0017) | within every limit |
| S7 heavy tails | 8 known controls | 500 | 1.01 | 0.042 (0.0004) | 0.035 (0.0009) | 0.004 [0.000, 0.014] | 0.052 | 0.919 (0.0016) | within every limit |
| S7 heavy tails | 5 first-pass nulls | 500 | 1.37 | 0.043 (0.0004) | 0.037 (0.0008) | 0.002 [0.000, 0.011] | 0.052 | 0.914 (0.0017) | within every limit |
| S7 heavy tails | 3 legacy | 500 | 1.37 | 0.067 (0.0008) | 0.084 (0.0019) | 0.290 [0.251, 0.332] | 0.052 | 0.927 (0.0026) | 2 fail; 3 fail |
| S7 heavy tails | 1 none | 500 | 0.00 | 0.010 (0.0002) | 0.003 (0.0003) | 0.000 [0.000, 0.007] | 0.052 | 0.534 (0.0024) | 4 fail |
| S8 known unbalanced batch | 2 ceiling | 500 | 1.00 | 0.051 (0.0004) | 0.055 (0.0009) | 0.020 [0.010, 0.036] | 0.070 | 0.872 (0.0014) | within every limit |
| S8 known unbalanced batch | **7 recommended** | 500 | 2.00 | 0.055 (0.0005) | 0.059 (0.0011) | 0.056 [0.038, 0.080] | 0.070 | 0.876 (0.0014) | 3 inconclusive |
| S8 known unbalanced batch | 8 known controls | 500 | 2.00 | 0.050 (0.0005) | 0.053 (0.0010) | 0.038 [0.023, 0.059] | 0.070 | 0.869 (0.0014) | within every limit |
| S8 known unbalanced batch | 5 first-pass nulls | 500 | 1.63 | 0.058 (0.0010) | 0.068 (0.0019) | 0.214 [0.179, 0.253] | 0.070 | 0.759 (0.0079) | 3 fail |
| S8 known unbalanced batch | 3 legacy | 500 | 1.63 | 0.183 (0.0051) | 0.346 (0.0106) | 0.728 [0.687, 0.767] | 0.070 | 0.755 (0.0092) | 1 fail; 2 fail; 3 fail |
| S8 known unbalanced batch | 1 none | 500 | 0.00 | 0.020 (0.0003) | 0.009 (0.0011) | 0.010 [0.003, 0.023] | 0.070 | 0.144 (0.0036) | 4 fail |
| S9 80% of effects up | 2 ceiling | 500 | 1.00 | 0.051 (0.0004) | 0.052 (0.0009) | 0.012 [0.004, 0.026] | 0.062 | 0.920 (0.0011) | within every limit |
| S9 80% of effects up | **7 recommended** | 500 | 1.00 | 0.053 (0.0004) | 0.056 (0.0009) | 0.018 [0.008, 0.034] | 0.062 | 0.920 (0.0011) | within every limit |
| S9 80% of effects up | 8 known controls | 500 | 1.00 | 0.050 (0.0004) | 0.051 (0.0009) | 0.018 [0.008, 0.034] | 0.062 | 0.918 (0.0011) | within every limit |
| S9 80% of effects up | 5 first-pass nulls | 500 | 1.10 | 0.051 (0.0004) | 0.053 (0.0009) | 0.016 [0.007, 0.031] | 0.062 | 0.919 (0.0011) | within every limit |
| S9 80% of effects up | 3 legacy † | 500 | 1.10 | 0.067 (0.0007) | 0.086 (0.0015) | 0.288 [0.249, 0.330] | 0.062 | 0.922 (0.0032) | 1 fail; 2 fail; 3 fail |
| S9 80% of effects up | 1 none | 500 | 0.00 | 0.015 (0.0002) | 0.007 (0.0004) | 0.000 [0.000, 0.007] | 0.062 | 0.516 (0.0025) | 4 fail |
| S10 loadings tied to abundance | 2 ceiling | 500 | 1.00 | 0.050 (0.0004) | 0.050 (0.0009) | 0.012 [0.004, 0.026] | 0.062 | 0.913 (0.0011) | within every limit |
| S10 loadings tied to abundance | **7 recommended** | 500 | 1.00 | 0.052 (0.0004) | 0.053 (0.0009) | 0.014 [0.006, 0.029] | 0.062 | 0.914 (0.0011) | within every limit |
| S10 loadings tied to abundance | 8 known controls | 500 | 1.00 | 0.050 (0.0004) | 0.050 (0.0009) | 0.012 [0.004, 0.026] | 0.062 | 0.912 (0.0012) | within every limit |
| S10 loadings tied to abundance | 5 first-pass nulls | 500 | 1.85 | 0.051 (0.0004) | 0.052 (0.0009) | 0.012 [0.004, 0.026] | 0.062 | 0.902 (0.0013) | within every limit |
| S10 loadings tied to abundance | 3 legacy † | 500 | 1.85 | 0.079 (0.0009) | 0.112 (0.0021) | 0.596 [0.552, 0.639] | 0.062 | 0.912 (0.0042) | 1 fail; 2 fail; 3 fail |
| S10 loadings tied to abundance | 1 none | 500 | 0.00 | 0.015 (0.0002) | 0.007 (0.0004) | 0.000 [0.000, 0.007] | 0.062 | 0.525 (0.0026) | 4 fail |
| S11 3 per group, factor tied 0.3 | 2 ceiling | 500 | 1.00 | 0.051 (0.0004) | 0.051 (0.0011) | 0.036 [0.021, 0.056] | 0.086 | 0.614 (0.0036) | within every limit |
| S11 3 per group, factor tied 0.3 | **7 recommended** | 500 | 1.00 | 0.052 (0.0004) | 0.055 (0.0011) | 0.032 [0.018, 0.051] | 0.086 | 0.628 (0.0036) | within every limit |
| S11 3 per group, factor tied 0.3 | 8 known controls | 500 | 1.00 | 0.052 (0.0005) | 0.054 (0.0012) | 0.044 [0.028, 0.066] | 0.086 | 0.608 (0.0042) | within every limit |
| S11 3 per group, factor tied 0.3 | 5 first-pass nulls | 500 | 1.51 | 0.061 (0.0005) | 0.067 (0.0014) | 0.128 [0.100, 0.160] | 0.086 | 0.601 (0.0044) | 3 fail |
| S11 3 per group, factor tied 0.3 | 3 legacy † | 500 | 1.51 | 0.299 (0.0060) | 0.521 (0.0085) | 0.964 [0.944, 0.979] | 0.086 | 0.852 (0.0069) | 1 fail; 2 fail; 3 fail |
| S11 3 per group, factor tied 0.3 | 1 none | 500 | 0.00 | 0.017 (0.0004) | 0.004 (0.0024) | 0.006 [0.001, 0.017] | 0.086 | 0.000 (0.0000) | 4 fail |
| S12 3 per group, three factors (**void**) | 2 ceiling | 500 | 3.00 | 0.078 (0.0009) | 0.127 (0.0049) | 0.632 [0.588, 0.674] | 0.682 | 0.382 (0.0128) | 1 fail; 2 fail |
| S12 3 per group, three factors (**void**) | **7 recommended** | 500 | 3.00 | 0.079 (0.0010) | 0.119 (0.0046) | 0.572 [0.527, 0.616] | 0.682 | 0.359 (0.0131) | 1 fail; 2 fail |
| S12 3 per group, three factors (**void**) | 8 known controls | 500 | 3.00 | 0.074 (0.0009) | 0.124 (0.0056) | 0.606 [0.562, 0.649] | 0.682 | 0.341 (0.0124) | 1 fail; 2 fail |
| S12 3 per group, three factors (**void**) | 5 first-pass nulls | 500 | 1.19 | 0.007 (0.0002) | 0.006 (0.0035) | 0.006 [0.001, 0.017] | 0.682 | 0.000 (0.0000) | 4 inconclusive |
| S12 3 per group, three factors (**void**) | 3 legacy † | 500 | 1.19 | 0.056 (0.0040) | 0.106 (0.0100) | 0.192 [0.158, 0.229] | 0.682 | 0.103 (0.0100) | 1 fail; 2 fail |
| S12 3 per group, three factors (**void**) | 1 none | 500 | 0.00 | 0.002 (0.0001) | 0.006 (0.0035) | 0.006 [0.001, 0.017] | 0.682 | 0.000 (0.0000) | 4 fail |
| S13 no factor, 3 per group | 2 ceiling | 500 | 0.00 | 0.050 (0.0004) | 0.051 (0.0010) | 0.026 [0.014, 0.044] | 0.076 | 0.746 (0.0026) | within every limit |
| S13 no factor, 3 per group | **7 recommended** | 500 | 0.13 | 0.051 (0.0005) | 0.053 (0.0011) | 0.034 [0.020, 0.054] | 0.076 | 0.741 (0.0029) | within every limit |
| S13 no factor, 3 per group | 8 known controls | 500 | 0.13 | 0.051 (0.0004) | 0.053 (0.0011) | 0.040 [0.025, 0.061] | 0.076 | 0.739 (0.0030) | within every limit |
| S13 no factor, 3 per group | 5 first-pass nulls | 500 | 1.86 | 0.057 (0.0005) | 0.067 (0.0013) | 0.144 [0.114, 0.178] | 0.076 | 0.617 (0.0044) | 3 fail; 5 fail |
| S13 no factor, 3 per group | 3 legacy † | 500 | 1.86 | 0.281 (0.0037) | 0.518 (0.0071) | 0.980 [0.964, 0.990] | 0.076 | 0.928 (0.0035) | 1 fail; 2 fail; 3 fail |
| S13 no factor, 3 per group | 1 none | 500 | 0.00 | 0.050 (0.0004) | 0.051 (0.0010) | 0.026 [0.014, 0.044] | 0.076 | 0.746 (0.0026) | within every limit |
| S14 no factor, 50% missing | 2 ceiling | 500 | 0.00 | 0.051 (0.0005) | 0.052 (0.0011) | 0.038 [0.023, 0.059] | 0.088 | 0.908 (0.0015) | within every limit |
| S14 no factor, 50% missing | **7 recommended** | 500 | 0.13 | 0.051 (0.0005) | 0.053 (0.0011) | 0.036 [0.021, 0.056] | 0.088 | 0.906 (0.0016) | within every limit |
| S14 no factor, 50% missing | 8 known controls | 500 | 0.13 | 0.051 (0.0005) | 0.052 (0.0011) | 0.040 [0.025, 0.061] | 0.088 | 0.906 (0.0016) | within every limit |
| S14 no factor, 50% missing | 5 first-pass nulls | 500 | 1.97 | 0.054 (0.0005) | 0.057 (0.0013) | 0.068 [0.048, 0.094] | 0.088 | 0.870 (0.0020) | 3 inconclusive; 5 fail |
| S14 no factor, 50% missing | 3 legacy † | 500 | 1.97 | 0.066 (0.0011) | 0.089 (0.0022) | 0.376 [0.333, 0.420] | 0.088 | 0.822 (0.0083) | 1 fail; 2 fail; 3 fail; 5 fail |
| S14 no factor, 50% missing | 1 none | 500 | 0.00 | 0.051 (0.0005) | 0.052 (0.0011) | 0.038 [0.023, 0.059] | 0.088 | 0.908 (0.0015) | within every limit |

**S2, DIA_LIMPA (all fourteen on the runner)**

| Scenario | Design | Scored runs | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Rules not passed |
|---|---|---|---|---|---|---|---|---|---|
| S1 3 per group | 2 ceiling | 500 | 1.00 | 0.048 (0.0004) | 0.045 (0.0012) | 0.020 [0.010, 0.036] | 0.070 | 0.541 (0.0083) | within every limit |
| S1 3 per group | **7 recommended** | 500 | 1.00 | 0.058 (0.0011) | 0.070 (0.0031) | 0.132 [0.104, 0.165] | 0.070 | 0.602 (0.0048) | 2 inconclusive; 3 fail |
| S1 3 per group | 8 known controls | 500 | 1.00 | 0.059 (0.0014) | 0.073 (0.0035) | 0.156 [0.125, 0.191] | 0.070 | 0.592 (0.0051) | 1 inconclusive; 2 inconclusive; 3 fail |
| S1 3 per group | 5 first-pass nulls | 500 | 1.73 | 0.093 (0.0027) | 0.158 (0.0067) | 0.476 [0.431, 0.521] | 0.070 | 0.654 (0.0026) | 1 fail; 2 fail; 3 fail |
| S1 3 per group | 3 legacy † | 500 | 1.73 | 0.286 (0.0063) | 0.505 (0.0086) | 0.986 [0.971, 0.994] | 0.070 | 0.826 (0.0029) | 1 fail; 2 fail; 3 fail |
| S1 3 per group | 1 none | 500 | 0.00 | 0.060 (0.0042) | 0.083 (0.0086) | 0.178 [0.145, 0.214] | 0.070 | 0.153 (0.0102) | 1 inconclusive; 2 inconclusive; 3 fail; 4 fail |
| S2 40% of features changed | 2 ceiling | 500 | 1.00 | 0.049 (0.0004) | 0.049 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.941 (0.0014) | within every limit |
| S2 40% of features changed | **7 recommended** | 500 | 1.00 | 0.050 (0.0004) | 0.051 (0.0006) | 0.002 [0.000, 0.011] | 0.050 | 0.943 (0.0012) | within every limit |
| S2 40% of features changed | 8 known controls | 500 | 1.00 | 0.050 (0.0005) | 0.050 (0.0007) | 0.006 [0.001, 0.017] | 0.050 | 0.942 (0.0012) | within every limit |
| S2 40% of features changed | 5 first-pass nulls | 500 | 1.85 | 0.083 (0.0018) | 0.100 (0.0026) | 0.390 [0.347, 0.434] | 0.050 | 0.948 (0.0009) | 1 fail; 2 fail; 3 fail |
| S2 40% of features changed | 3 legacy † | 500 | 1.85 | 0.207 (0.0067) | 0.248 (0.0071) | 0.758 [0.718, 0.795] | 0.050 | 0.938 (0.0015) | 1 fail; 2 fail; 3 fail |
| S2 40% of features changed | 1 none | 500 | 0.00 | 0.053 (0.0038) | 0.070 (0.0056) | 0.184 [0.151, 0.221] | 0.050 | 0.692 (0.0057) | 1 inconclusive; 2 inconclusive; 3 fail; 4 fail |
| S3 factor tied to groups, 0.6 | 2 ceiling | 500 | 1.00 | 0.048 (0.0003) | 0.047 (0.0008) | 0.006 [0.001, 0.017] | 0.056 | 0.741 (0.0076) | within every limit |
| S3 factor tied to groups, 0.6 | **7 recommended** | 500 | 1.00 | 0.056 (0.0005) | 0.065 (0.0013) | 0.084 [0.061, 0.112] | 0.056 | 0.786 (0.0055) | 3 fail |
| S3 factor tied to groups, 0.6 | 8 known controls | 500 | 1.00 | 0.059 (0.0009) | 0.073 (0.0021) | 0.180 [0.147, 0.217] | 0.056 | 0.782 (0.0056) | 2 inconclusive; 3 fail |
| S3 factor tied to groups, 0.6 | 5 first-pass nulls | 500 | 1.86 | 0.154 (0.0026) | 0.301 (0.0063) | 0.894 [0.864, 0.920] | 0.056 | 0.858 (0.0028) | 1 fail; 2 fail; 3 fail |
| S3 factor tied to groups, 0.6 | 3 legacy † | 500 | 1.86 | 0.469 (0.0069) | 0.677 (0.0067) | 0.994 [0.983, 0.999] | 0.056 | 0.865 (0.0021) | 1 fail; 2 fail; 3 fail |
| S3 factor tied to groups, 0.6 | 1 none | 500 | 0.00 | 0.324 (0.0082) | 0.533 (0.0125) | 0.836 [0.801, 0.867] | 0.056 | 0.745 (0.0034) | 1 fail; 2 fail; 3 fail; 4 fail |
| S4 three factors | 2 ceiling | 500 | 3.00 | 0.049 (0.0004) | 0.049 (0.0008) | 0.004 [0.000, 0.014] | 0.054 | 0.824 (0.0053) | within every limit |
| S4 three factors | **7 recommended** | 500 | 2.63 | 0.066 (0.0018) | 0.088 (0.0042) | 0.188 [0.155, 0.225] | 0.054 | 0.867 (0.0024) | 1 fail; 2 fail; 3 fail |
| S4 three factors | 8 known controls | 500 | 2.63 | 0.066 (0.0019) | 0.088 (0.0042) | 0.224 [0.188, 0.263] | 0.054 | 0.855 (0.0026) | 1 fail; 2 fail; 3 fail |
| S4 three factors | 5 first-pass nulls | 500 | 1.86 | 0.094 (0.0033) | 0.152 (0.0073) | 0.432 [0.388, 0.477] | 0.054 | 0.821 (0.0037) | 1 fail; 2 fail; 3 fail |
| S4 three factors | 3 legacy † | 500 | 1.86 | 0.173 (0.0048) | 0.316 (0.0096) | 0.788 [0.750, 0.823] | 0.054 | 0.843 (0.0036) | 1 fail; 2 fail; 3 fail |
| S4 three factors | 1 none | 500 | 0.00 | 0.055 (0.0027) | 0.073 (0.0058) | 0.228 [0.192, 0.267] | 0.054 | 0.455 (0.0078) | 2 inconclusive; 3 fail; 4 fail |
| S5 no missing values | 2 ceiling | 248 | 1.00 | 0.049 (0.0005) | 0.049 (0.0010) | 0.000 [0.000, 0.015] | 0.050 | 0.997 (0.0009) | 6 fail |
| S5 no missing values | **7 recommended** | 248 | 1.00 | 0.051 (0.0005) | 0.050 (0.0010) | 0.004 [0.000, 0.022] | 0.050 | 0.997 (0.0007) | 6 fail |
| S5 no missing values | 8 known controls | 248 | 1.00 | 0.050 (0.0006) | 0.050 (0.0012) | 0.012 [0.003, 0.035] | 0.050 | 0.997 (0.0008) | 6 fail |
| S5 no missing values | 5 first-pass nulls | 248 | 1.71 | 0.068 (0.0016) | 0.088 (0.0037) | 0.262 [0.208, 0.321] | 0.050 | 0.997 (0.0005) | 1 fail; 2 fail; 3 fail; 6 fail |
| S5 no missing values | 3 legacy † | 248 | 1.71 | 0.228 (0.0107) | 0.367 (0.0158) | 0.819 [0.765, 0.864] | 0.050 | 0.980 (0.0019) | 1 fail; 2 fail; 3 fail; 6 fail |
| S5 no missing values | 1 none | 248 | 0.00 | 0.044 (0.0059) | 0.066 (0.0107) | 0.141 [0.100, 0.191] | 0.050 | 0.676 (0.0086) | 1 inconclusive; 2 inconclusive; 3 fail; 4 fail; 6 fail |
| S6 50% missing | 2 ceiling | 500 | 1.00 | 0.063 (0.0005) | 0.068 (0.0013) | 0.116 [0.089, 0.147] | 0.166 | 0.514 (0.0031) | within every limit |
| S6 50% missing | **7 recommended** | 500 | 1.00 | 0.066 (0.0005) | 0.071 (0.0014) | 0.138 [0.109, 0.171] | 0.166 | 0.524 (0.0025) | 1 inconclusive; 3 inconclusive |
| S6 50% missing | 8 known controls | 500 | 1.00 | 0.066 (0.0006) | 0.073 (0.0016) | 0.162 [0.131, 0.197] | 0.166 | 0.517 (0.0025) | 1 inconclusive; 2 inconclusive; 3 inconclusive |
| S6 50% missing | 5 first-pass nulls | 500 | 1.84 | 0.081 (0.0012) | 0.109 (0.0033) | 0.404 [0.361, 0.448] | 0.166 | 0.544 (0.0021) | 1 fail; 2 fail; 3 fail |
| S6 50% missing | 3 legacy † | 500 | 1.84 | 0.164 (0.0040) | 0.307 (0.0091) | 0.862 [0.829, 0.891] | 0.166 | 0.568 (0.0032) | 1 fail; 2 fail; 3 fail |
| S6 50% missing | 1 none | 500 | 0.00 | 0.095 (0.0026) | 0.136 (0.0064) | 0.428 [0.384, 0.473] | 0.166 | 0.318 (0.0043) | 1 fail; 2 fail; 3 fail; 4 fail |
| S7 heavy tails | 2 ceiling | 500 | 1.00 | 0.050 (0.0004) | 0.048 (0.0007) | 0.002 [0.000, 0.011] | 0.052 | 0.925 (0.0015) | within every limit |
| S7 heavy tails | **7 recommended** | 500 | 1.00 | 0.052 (0.0004) | 0.051 (0.0009) | 0.016 [0.007, 0.031] | 0.052 | 0.928 (0.0013) | within every limit |
| S7 heavy tails | 8 known controls | 500 | 1.00 | 0.052 (0.0005) | 0.051 (0.0011) | 0.040 [0.025, 0.061] | 0.052 | 0.926 (0.0013) | 3 inconclusive |
| S7 heavy tails | 5 first-pass nulls | 500 | 1.76 | 0.072 (0.0014) | 0.096 (0.0033) | 0.306 [0.266, 0.348] | 0.052 | 0.934 (0.0010) | 1 fail; 2 fail; 3 fail |
| S7 heavy tails | 3 legacy † | 500 | 1.76 | 0.212 (0.0064) | 0.359 (0.0105) | 0.814 [0.777, 0.847] | 0.052 | 0.925 (0.0020) | 1 fail; 2 fail; 3 fail |
| S7 heavy tails | 1 none | 500 | 0.00 | 0.047 (0.0036) | 0.068 (0.0072) | 0.146 [0.116, 0.180] | 0.052 | 0.602 (0.0067) | 2 inconclusive; 3 fail; 4 fail |
| S8 known unbalanced batch | 2 ceiling | 500 | 1.00 | 0.048 (0.0003) | 0.049 (0.0008) | 0.002 [0.000, 0.011] | 0.052 | 0.861 (0.0025) | within every limit |
| S8 known unbalanced batch | **7 recommended** | 500 | 1.98 | 0.052 (0.0010) | 0.055 (0.0019) | 0.026 [0.014, 0.044] | 0.052 | 0.870 (0.0020) | within every limit |
| S8 known unbalanced batch | 8 known controls | 500 | 1.98 | 0.052 (0.0011) | 0.057 (0.0021) | 0.058 [0.039, 0.082] | 0.052 | 0.861 (0.0021) | 3 inconclusive |
| S8 known unbalanced batch | 5 first-pass nulls | 500 | 1.78 | 0.100 (0.0024) | 0.164 (0.0054) | 0.726 [0.685, 0.765] | 0.052 | 0.830 (0.0057) | 1 fail; 2 fail; 3 fail |
| S8 known unbalanced batch | 3 legacy † | 500 | 1.78 | 0.314 (0.0058) | 0.547 (0.0086) | 0.938 [0.913, 0.957] | 0.052 | 0.827 (0.0041) | 1 fail; 2 fail; 3 fail |
| S8 known unbalanced batch | 1 none | 500 | 0.00 | 0.043 (0.0025) | 0.040 (0.0050) | 0.090 [0.066, 0.119] | 0.052 | 0.220 (0.0074) | 3 fail; 4 fail |
| S9 80% of effects up | 2 ceiling | 500 | 1.00 | 0.048 (0.0004) | 0.047 (0.0008) | 0.004 [0.000, 0.014] | 0.054 | 0.900 (0.0022) | within every limit |
| S9 80% of effects up | **7 recommended** | 500 | 1.00 | 0.050 (0.0004) | 0.050 (0.0008) | 0.008 [0.002, 0.020] | 0.054 | 0.904 (0.0018) | within every limit |
| S9 80% of effects up | 8 known controls | 500 | 1.00 | 0.050 (0.0005) | 0.051 (0.0010) | 0.030 [0.017, 0.049] | 0.054 | 0.902 (0.0018) | within every limit |
| S9 80% of effects up | 5 first-pass nulls | 500 | 1.52 | 0.073 (0.0014) | 0.101 (0.0035) | 0.330 [0.289, 0.373] | 0.054 | 0.912 (0.0012) | 1 fail; 2 fail; 3 fail |
| S9 80% of effects up | 3 legacy † | 500 | 1.52 | 0.207 (0.0063) | 0.351 (0.0106) | 0.792 [0.754, 0.827] | 0.054 | 0.911 (0.0020) | 1 fail; 2 fail; 3 fail |
| S9 80% of effects up | 1 none | 500 | 0.00 | 0.058 (0.0041) | 0.088 (0.0080) | 0.192 [0.158, 0.229] | 0.054 | 0.561 (0.0074) | 1 inconclusive; 2 inconclusive; 3 fail; 4 fail |
| S10 loadings tied to abundance | 2 ceiling | 500 | 1.00 | 0.048 (0.0003) | 0.048 (0.0008) | 0.002 [0.000, 0.011] | 0.052 | 0.902 (0.0023) | within every limit |
| S10 loadings tied to abundance | **7 recommended** | 500 | 1.00 | 0.050 (0.0004) | 0.052 (0.0009) | 0.010 [0.003, 0.023] | 0.052 | 0.907 (0.0018) | within every limit |
| S10 loadings tied to abundance | 8 known controls | 500 | 1.00 | 0.050 (0.0004) | 0.052 (0.0010) | 0.026 [0.014, 0.044] | 0.052 | 0.905 (0.0017) | within every limit |
| S10 loadings tied to abundance | 5 first-pass nulls | 500 | 1.87 | 0.072 (0.0013) | 0.100 (0.0032) | 0.344 [0.302, 0.387] | 0.052 | 0.916 (0.0012) | 1 fail; 2 fail; 3 fail |
| S10 loadings tied to abundance | 3 legacy † | 500 | 1.87 | 0.209 (0.0065) | 0.353 (0.0105) | 0.828 [0.792, 0.860] | 0.052 | 0.909 (0.0022) | 1 fail; 2 fail; 3 fail |
| S10 loadings tied to abundance | 1 none | 500 | 0.00 | 0.059 (0.0042) | 0.088 (0.0081) | 0.202 [0.168, 0.240] | 0.052 | 0.578 (0.0073) | 1 inconclusive; 2 inconclusive; 3 fail; 4 fail |
| S11 3 per group, factor tied 0.3 | 2 ceiling | 500 | 1.00 | 0.048 (0.0004) | 0.046 (0.0024) | 0.034 [0.020, 0.054] | 0.084 | 0.482 (0.0095) | within every limit |
| S11 3 per group, factor tied 0.3 | **7 recommended** | 500 | 1.00 | 0.061 (0.0014) | 0.078 (0.0036) | 0.164 [0.133, 0.199] | 0.084 | 0.560 (0.0062) | 1 inconclusive; 2 inconclusive; 3 fail |
| S11 3 per group, factor tied 0.3 | 8 known controls | 500 | 1.00 | 0.061 (0.0016) | 0.082 (0.0041) | 0.180 [0.147, 0.217] | 0.084 | 0.548 (0.0065) | 1 inconclusive; 2 inconclusive; 3 fail |
| S11 3 per group, factor tied 0.3 | 5 first-pass nulls | 500 | 1.75 | 0.113 (0.0033) | 0.210 (0.0080) | 0.600 [0.556, 0.643] | 0.084 | 0.644 (0.0027) | 1 fail; 2 fail; 3 fail |
| S11 3 per group, factor tied 0.3 | 3 legacy † | 500 | 1.75 | 0.338 (0.0070) | 0.567 (0.0085) | 0.982 [0.966, 0.992] | 0.084 | 0.824 (0.0030) | 1 fail; 2 fail; 3 fail |
| S11 3 per group, factor tied 0.3 | 1 none | 500 | 0.00 | 0.093 (0.0055) | 0.142 (0.0112) | 0.270 [0.232, 0.311] | 0.084 | 0.206 (0.0111) | 1 fail; 2 fail; 3 fail; 4 fail |
| S12 3 per group, three factors | 2 ceiling | 500 | 3.00 | 0.046 (0.0005) | 0.053 (0.0058) | 0.114 [0.087, 0.145] | 0.164 | 0.173 (0.0092) | within every limit |
| S12 3 per group, three factors | **7 recommended** | 500 | 1.39 | 0.092 (0.0035) | 0.164 (0.0086) | 0.448 [0.404, 0.493] | 0.164 | 0.390 (0.0087) | 1 fail; 2 fail; 3 fail |
| S12 3 per group, three factors | 8 known controls | 500 | 1.39 | 0.092 (0.0035) | 0.165 (0.0086) | 0.448 [0.404, 0.493] | 0.164 | 0.378 (0.0086) | 1 fail; 2 fail; 3 fail |
| S12 3 per group, three factors | 5 first-pass nulls | 500 | 1.79 | 0.122 (0.0038) | 0.236 (0.0091) | 0.648 [0.604, 0.690] | 0.164 | 0.494 (0.0088) | 1 fail; 2 fail; 3 fail |
| S12 3 per group, three factors | 3 legacy † | 500 | 1.79 | 0.283 (0.0060) | 0.513 (0.0100) | 0.900 [0.870, 0.925] | 0.164 | 0.725 (0.0086) | 1 fail; 2 fail; 3 fail |
| S12 3 per group, three factors | 1 none | 500 | 0.00 | 0.057 (0.0027) | 0.078 (0.0067) | 0.226 [0.190, 0.265] | 0.164 | 0.100 (0.0066) | 1 inconclusive; 2 inconclusive; 3 fail; 4 fail |
| S13 no factor, 3 per group | 2 ceiling | 500 | 0.00 | 0.047 (0.0004) | 0.047 (0.0009) | 0.008 [0.002, 0.020] | 0.058 | 0.700 (0.0017) | within every limit |
| S13 no factor, 3 per group | **7 recommended** | 500 | 0.14 | 0.048 (0.0004) | 0.048 (0.0009) | 0.010 [0.003, 0.023] | 0.058 | 0.702 (0.0018) | within every limit |
| S13 no factor, 3 per group | 8 known controls | 500 | 0.14 | 0.048 (0.0004) | 0.047 (0.0009) | 0.010 [0.003, 0.023] | 0.058 | 0.700 (0.0018) | within every limit |
| S13 no factor, 3 per group | 5 first-pass nulls | 500 | 1.99 | 0.053 (0.0004) | 0.056 (0.0011) | 0.036 [0.021, 0.056] | 0.058 | 0.702 (0.0020) | within every limit |
| S13 no factor, 3 per group | 3 legacy † | 500 | 1.99 | 0.167 (0.0013) | 0.326 (0.0033) | 0.994 [0.983, 0.999] | 0.058 | 0.854 (0.0026) | 1 fail; 2 fail; 3 fail |
| S13 no factor, 3 per group | 1 none | 500 | 0.00 | 0.047 (0.0004) | 0.047 (0.0009) | 0.008 [0.002, 0.020] | 0.058 | 0.700 (0.0017) | within every limit |
| S14 no factor, 50% missing (**void**) | 2 ceiling | 500 | 0.00 | 0.064 (0.0004) | 0.074 (0.0011) | 0.138 [0.109, 0.171] | 0.188 | 0.551 (0.0017) | 1 fail; 2 inconclusive |
| S14 no factor, 50% missing (**void**) | **7 recommended** | 500 | 0.19 | 0.065 (0.0004) | 0.075 (0.0011) | 0.150 [0.120, 0.184] | 0.188 | 0.553 (0.0017) | 1 fail; 2 inconclusive |
| S14 no factor, 50% missing (**void**) | 8 known controls | 500 | 0.19 | 0.065 (0.0004) | 0.075 (0.0011) | 0.154 [0.123, 0.189] | 0.188 | 0.551 (0.0017) | 1 fail; 2 inconclusive; 3 inconclusive |
| S14 no factor, 50% missing (**void**) | 5 first-pass nulls | 500 | 2.00 | 0.069 (0.0005) | 0.079 (0.0012) | 0.202 [0.168, 0.240] | 0.188 | 0.566 (0.0017) | 1 fail; 2 fail; 3 inconclusive |
| S14 no factor, 50% missing (**void**) | 3 legacy † | 500 | 2.00 | 0.080 (0.0007) | 0.110 (0.0019) | 0.568 [0.523, 0.612] | 0.188 | 0.594 (0.0022) | 1 fail; 2 fail; 3 fail |
| S14 no factor, 50% missing (**void**) | 1 none | 500 | 0.00 | 0.064 (0.0004) | 0.074 (0.0011) | 0.138 [0.109, 0.171] | 0.188 | 0.551 (0.0017) | 1 fail; 2 inconclusive |

**S3, transcriptomics (S14 on the runner, the rest laptop; S6 not defined)**

| Scenario | Design | Scored runs | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Rules not passed |
|---|---|---|---|---|---|---|---|---|---|
| S1 3 per group | 2 ceiling | 500 | 1.00 | 0.053 (0.0002) | 0.046 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.369 (0.0013) | within every limit |
| S1 3 per group | **7 recommended** | 500 | 1.00 | 0.053 (0.0002) | 0.046 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.369 (0.0013) | within every limit |
| S1 3 per group | 8 known controls | 500 | 1.00 | 0.054 (0.0003) | 0.047 (0.0008) | 0.010 [0.003, 0.023] | 0.050 | 0.345 (0.0017) | within every limit |
| S1 3 per group | 5 first-pass nulls | 500 | 2.66 | 0.052 (0.0002) | 0.044 (0.0012) | 0.014 [0.006, 0.029] | 0.050 | 0.250 (0.0060) | within every limit |
| S1 3 per group | 3 legacy | 500 | 2.66 | 0.052 (0.0003) | 0.046 (0.0008) | 0.010 [0.003, 0.023] | 0.050 | 0.306 (0.0033) | within every limit |
| S1 3 per group | 1 none | 500 | 0.00 | 0.004 (0.0000) | 0.000 (0.0000) | 0.000 [0.000, 0.007] | 0.050 | 0.000 (0.0000) | 4 fail |
| S2 40% of features changed | 2 ceiling | 500 | 1.00 | 0.056 (0.0002) | 0.035 (0.0002) | 0.000 [0.000, 0.007] | 0.050 | 0.835 (0.0005) | within every limit |
| S2 40% of features changed | **7 recommended** | 500 | 1.00 | 0.056 (0.0002) | 0.035 (0.0002) | 0.000 [0.000, 0.007] | 0.050 | 0.835 (0.0005) | within every limit |
| S2 40% of features changed | 8 known controls | 500 | 1.00 | 0.056 (0.0003) | 0.035 (0.0003) | 0.000 [0.000, 0.007] | 0.050 | 0.827 (0.0005) | within every limit |
| S2 40% of features changed | 5 first-pass nulls | 500 | 2.78 | 0.053 (0.0002) | 0.035 (0.0075) | 0.036 [0.021, 0.056] | 0.050 | 0.092 (0.0117) | 3 inconclusive; 4 inconclusive |
| S2 40% of features changed | 3 legacy | 500 | 2.78 | 0.057 (0.0003) | 0.036 (0.0003) | 0.000 [0.000, 0.007] | 0.050 | 0.793 (0.0025) | within every limit |
| S2 40% of features changed | 1 none | 500 | 0.00 | 0.010 (0.0001) | 0.001 (0.0001) | 0.000 [0.000, 0.007] | 0.050 | 0.114 (0.0012) | 4 fail |
| S3 factor tied to groups, 0.6 | 2 ceiling | 500 | 1.00 | 0.055 (0.0002) | 0.048 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.553 (0.0012) | within every limit |
| S3 factor tied to groups, 0.6 | **7 recommended** | 500 | 1.00 | 0.058 (0.0002) | 0.052 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.575 (0.0011) | within every limit |
| S3 factor tied to groups, 0.6 | 8 known controls | 500 | 1.00 | 0.055 (0.0003) | 0.048 (0.0006) | 0.004 [0.000, 0.014] | 0.050 | 0.538 (0.0016) | within every limit |
| S3 factor tied to groups, 0.6 | 5 first-pass nulls | 500 | 3.00 | 0.055 (0.0002) | 0.062 (0.0101) | 0.076 [0.054, 0.103] | 0.050 | 0.002 (0.0012) | 2 inconclusive; 3 fail |
| S3 factor tied to groups, 0.6 | 3 legacy | 500 | 3.00 | 0.509 (0.0022) | 0.724 (0.0016) | 1.000 [0.993, 1.000] | 0.050 | 0.676 (0.0010) | 1 fail; 2 fail; 3 fail |
| S3 factor tied to groups, 0.6 | 1 none | 500 | 0.00 | 0.337 (0.0010) | 0.450 (0.0030) | 1.000 [0.993, 1.000] | 0.050 | 0.483 (0.0012) | 1 fail; 2 fail; 3 fail; 4 fail |
| S4 three factors | 2 ceiling | 500 | 3.00 | 0.056 (0.0002) | 0.048 (0.0004) | 0.000 [0.000, 0.007] | 0.050 | 0.763 (0.0008) | within every limit |
| S4 three factors | **7 recommended** | 500 | 3.00 | 0.056 (0.0002) | 0.049 (0.0004) | 0.000 [0.000, 0.007] | 0.050 | 0.762 (0.0008) | within every limit |
| S4 three factors | 8 known controls | 500 | 3.00 | 0.056 (0.0005) | 0.048 (0.0008) | 0.012 [0.004, 0.026] | 0.050 | 0.718 (0.0012) | within every limit |
| S4 three factors | 5 first-pass nulls | 500 | 1.55 | 0.017 (0.0009) | 0.009 (0.0008) | 0.000 [0.000, 0.007] | 0.050 | 0.180 (0.0129) | within every limit |
| S4 three factors | 3 legacy | 500 | 1.55 | 0.016 (0.0008) | 0.008 (0.0008) | 0.002 [0.000, 0.011] | 0.050 | 0.160 (0.0121) | within every limit |
| S4 three factors | 1 none | 500 | 0.00 | 0.005 (0.0002) | 0.000 (0.0000) | 0.000 [0.000, 0.007] | 0.050 | 0.000 (0.0000) | 4 fail |
| S5 low counts (library size / 10) | 2 ceiling | 500 | 1.00 | 0.056 (0.0004) | 0.049 (0.0007) | 0.000 [0.000, 0.007] | 0.050 | 0.754 (0.0013) | within every limit |
| S5 low counts (library size / 10) | **7 recommended** | 500 | 1.00 | 0.055 (0.0003) | 0.047 (0.0006) | 0.000 [0.000, 0.007] | 0.050 | 0.753 (0.0013) | within every limit |
| S5 low counts (library size / 10) | 8 known controls | 500 | 1.00 | 0.054 (0.0005) | 0.047 (0.0008) | 0.012 [0.004, 0.026] | 0.050 | 0.735 (0.0014) | within every limit |
| S5 low counts (library size / 10) | 5 first-pass nulls | 500 | 3.00 | 0.053 (0.0003) | 0.051 (0.0026) | 0.034 [0.020, 0.054] | 0.050 | 0.440 (0.0120) | 3 inconclusive |
| S5 low counts (library size / 10) | 3 legacy | 500 | 3.00 | 0.056 (0.0004) | 0.050 (0.0007) | 0.002 [0.000, 0.011] | 0.050 | 0.688 (0.0028) | within every limit |
| S5 low counts (library size / 10) | 1 none | 500 | 0.00 | 0.009 (0.0003) | 0.001 (0.0004) | 0.002 [0.000, 0.011] | 0.050 | 0.002 (0.0002) | 4 fail |
| S7 heavy tails (**void**) | 2 ceiling | 500 | 1.00 | 0.085 (0.0003) | 0.022 (0.0004) | 0.000 [0.000, 0.007] | 0.050 | 0.352 (0.0013) | 1 fail |
| S7 heavy tails (**void**) | **7 recommended** | 500 | 1.00 | 0.085 (0.0003) | 0.022 (0.0004) | 0.000 [0.000, 0.007] | 0.050 | 0.352 (0.0013) | 1 fail |
| S7 heavy tails (**void**) | 8 known controls | 500 | 1.00 | 0.086 (0.0004) | 0.022 (0.0005) | 0.002 [0.000, 0.011] | 0.050 | 0.316 (0.0017) | 1 fail |
| S7 heavy tails (**void**) | 5 first-pass nulls | 500 | 2.97 | 0.069 (0.0005) | 0.020 (0.0026) | 0.022 [0.011, 0.039] | 0.050 | 0.142 (0.0046) | within every limit |
| S7 heavy tails (**void**) | 3 legacy | 500 | 2.97 | 0.064 (0.0003) | 0.018 (0.0006) | 0.000 [0.000, 0.007] | 0.050 | 0.204 (0.0033) | within every limit |
| S7 heavy tails (**void**) | 1 none | 500 | 0.00 | 0.064 (0.0003) | 0.000 (0.0000) | 0.000 [0.000, 0.007] | 0.050 | 0.000 (0.0000) | 4 fail |
| S8 known unbalanced batch | 2 ceiling | 500 | 1.00 | 0.056 (0.0002) | 0.048 (0.0004) | 0.000 [0.000, 0.007] | 0.050 | 0.725 (0.0009) | within every limit |
| S8 known unbalanced batch | **7 recommended** | 500 | 1.00 | 0.056 (0.0002) | 0.049 (0.0004) | 0.000 [0.000, 0.007] | 0.050 | 0.728 (0.0009) | within every limit |
| S8 known unbalanced batch | 8 known controls | 500 | 1.00 | 0.056 (0.0003) | 0.049 (0.0006) | 0.004 [0.000, 0.014] | 0.050 | 0.712 (0.0011) | within every limit |
| S8 known unbalanced batch | 5 first-pass nulls | 500 | 3.00 | 0.058 (0.0002) | 0.067 (0.0060) | 0.184 [0.151, 0.221] | 0.050 | 0.042 (0.0038) | 2 inconclusive; 3 fail; 4 fail |
| S8 known unbalanced batch | 3 legacy | 500 | 3.00 | 0.241 (0.0011) | 0.433 (0.0020) | 1.000 [0.993, 1.000] | 0.050 | 0.692 (0.0016) | 1 fail; 2 fail; 3 fail |
| S8 known unbalanced batch | 1 none | 500 | 0.00 | 0.036 (0.0001) | 0.004 (0.0003) | 0.000 [0.000, 0.007] | 0.050 | 0.083 (0.0017) | 4 fail |
| S9 80% of effects up (**void**) | 2 ceiling | 500 | 1.00 | 0.075 (0.0004) | 0.074 (0.0006) | 0.040 [0.025, 0.061] | 0.090 | 0.726 (0.0009) | 1 fail; 2 inconclusive |
| S9 80% of effects up (**void**) | **7 recommended** | 500 | 1.00 | 0.075 (0.0004) | 0.074 (0.0006) | 0.038 [0.023, 0.059] | 0.090 | 0.726 (0.0009) | 1 fail; 2 inconclusive |
| S9 80% of effects up (**void**) | 8 known controls | 500 | 1.00 | 0.075 (0.0005) | 0.075 (0.0009) | 0.094 [0.070, 0.123] | 0.090 | 0.711 (0.0010) | 1 fail; 2 inconclusive; 3 inconclusive |
| S9 80% of effects up (**void**) | 5 first-pass nulls | 500 | 2.67 | 0.063 (0.0003) | 0.058 (0.0017) | 0.014 [0.006, 0.029] | 0.090 | 0.473 (0.0109) | within every limit |
| S9 80% of effects up (**void**) | 3 legacy | 500 | 2.67 | 0.074 (0.0004) | 0.074 (0.0007) | 0.044 [0.028, 0.066] | 0.090 | 0.678 (0.0028) | 1 fail |
| S9 80% of effects up (**void**) | 1 none | 500 | 0.00 | 0.014 (0.0002) | 0.000 (0.0000) | 0.000 [0.000, 0.007] | 0.090 | 0.001 (0.0002) | 4 fail |
| S10 loadings tied to abundance | 2 ceiling | 500 | 1.00 | 0.056 (0.0002) | 0.048 (0.0004) | 0.000 [0.000, 0.007] | 0.050 | 0.774 (0.0007) | within every limit |
| S10 loadings tied to abundance | **7 recommended** | 500 | 1.00 | 0.056 (0.0002) | 0.048 (0.0004) | 0.000 [0.000, 0.007] | 0.050 | 0.774 (0.0007) | within every limit |
| S10 loadings tied to abundance | 8 known controls | 500 | 1.00 | 0.055 (0.0003) | 0.048 (0.0006) | 0.006 [0.001, 0.017] | 0.050 | 0.764 (0.0008) | within every limit |
| S10 loadings tied to abundance | 5 first-pass nulls | 500 | 3.00 | 0.055 (0.0002) | 0.049 (0.0008) | 0.006 [0.001, 0.017] | 0.050 | 0.623 (0.0089) | within every limit |
| S10 loadings tied to abundance | 3 legacy | 500 | 3.00 | 0.055 (0.0002) | 0.047 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.727 (0.0024) | within every limit |
| S10 loadings tied to abundance | 1 none | 500 | 0.00 | 0.010 (0.0001) | 0.000 (0.0002) | 0.000 [0.000, 0.007] | 0.050 | 0.026 (0.0011) | 4 fail |
| S11 3 per group, factor tied 0.3 | 2 ceiling | 500 | 1.00 | 0.053 (0.0002) | 0.047 (0.0006) | 0.000 [0.000, 0.007] | 0.050 | 0.315 (0.0013) | within every limit |
| S11 3 per group, factor tied 0.3 | **7 recommended** | 500 | 1.00 | 0.055 (0.0002) | 0.049 (0.0006) | 0.000 [0.000, 0.007] | 0.050 | 0.324 (0.0013) | within every limit |
| S11 3 per group, factor tied 0.3 | 8 known controls | 500 | 1.00 | 0.054 (0.0004) | 0.049 (0.0009) | 0.012 [0.004, 0.026] | 0.050 | 0.289 (0.0020) | within every limit |
| S11 3 per group, factor tied 0.3 | 5 first-pass nulls | 500 | 2.54 | 0.061 (0.0003) | 0.065 (0.0029) | 0.080 [0.058, 0.107] | 0.050 | 0.230 (0.0056) | 3 fail |
| S11 3 per group, factor tied 0.3 | 3 legacy | 500 | 2.54 | 0.108 (0.0018) | 0.186 (0.0048) | 0.728 [0.687, 0.767] | 0.050 | 0.334 (0.0035) | 1 fail; 2 fail; 3 fail |
| S11 3 per group, factor tied 0.3 | 1 none | 500 | 0.00 | 0.013 (0.0002) | 0.000 (0.0000) | 0.000 [0.000, 0.007] | 0.050 | 0.000 (0.0000) | 4 fail |
| S12 3 per group, three factors | 2 ceiling | 500 | 3.00 | 0.051 (0.0002) | 0.044 (0.0006) | 0.000 [0.000, 0.007] | 0.050 | 0.341 (0.0017) | within every limit |
| S12 3 per group, three factors | **7 recommended** | 500 | 2.00 | 0.005 (0.0001) | 0.000 (0.0000) | 0.000 [0.000, 0.007] | 0.050 | 0.000 (0.0000) | 4 fail |
| S12 3 per group, three factors | 8 known controls | 500 | 2.00 | 0.005 (0.0001) | 0.000 (0.0000) | 0.000 [0.000, 0.007] | 0.050 | 0.000 (0.0000) | 4 fail |
| S12 3 per group, three factors | 5 first-pass nulls | 500 | 1.95 | 0.021 (0.0010) | 0.016 (0.0010) | 0.000 [0.000, 0.007] | 0.050 | 0.123 (0.0073) | within every limit |
| S12 3 per group, three factors | 3 legacy | 500 | 1.95 | 0.021 (0.0011) | 0.017 (0.0013) | 0.022 [0.011, 0.039] | 0.050 | 0.094 (0.0063) | within every limit |
| S12 3 per group, three factors | 1 none | 500 | 0.00 | 0.003 (0.0001) | 0.000 (0.0000) | 0.000 [0.000, 0.007] | 0.050 | 0.000 (0.0000) | 4 fail |
| S13 no factor, 3 per group | 2 ceiling | 500 | 0.00 | 0.053 (0.0002) | 0.044 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.386 (0.0011) | within every limit |
| S13 no factor, 3 per group | **7 recommended** | 500 | 0.57 | 0.052 (0.0002) | 0.044 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.380 (0.0012) | within every limit |
| S13 no factor, 3 per group | 8 known controls | 500 | 0.57 | 0.052 (0.0002) | 0.044 (0.0013) | 0.018 [0.008, 0.034] | 0.050 | 0.310 (0.0054) | 5 fail |
| S13 no factor, 3 per group | 5 first-pass nulls | 500 | 2.04 | 0.052 (0.0002) | 0.044 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.359 (0.0014) | 5 fail |
| S13 no factor, 3 per group | 3 legacy | 500 | 2.04 | 0.052 (0.0002) | 0.043 (0.0006) | 0.000 [0.000, 0.007] | 0.050 | 0.324 (0.0031) | 5 fail |
| S13 no factor, 3 per group | 1 none | 500 | 0.00 | 0.053 (0.0002) | 0.044 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.386 (0.0011) | within every limit |
| S14 no factor, low counts | 2 ceiling | 500 | 0.00 | 0.052 (0.0002) | 0.043 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.788 (0.0011) | within every limit |
| S14 no factor, low counts | **7 recommended** | 500 | 2.95 | 0.054 (0.0003) | 0.046 (0.0005) | 0.002 [0.000, 0.011] | 0.050 | 0.784 (0.0012) | within every limit |
| S14 no factor, low counts | 8 known controls | 500 | 2.95 | 0.051 (0.0002) | 0.042 (0.0008) | 0.010 [0.003, 0.023] | 0.050 | 0.569 (0.0091) | 5 fail |
| S14 no factor, low counts | 5 first-pass nulls | 500 | 2.07 | 0.053 (0.0002) | 0.044 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.773 (0.0013) | within every limit |
| S14 no factor, low counts | 3 legacy † | 500 | 2.07 | 0.052 (0.0002) | 0.044 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.749 (0.0030) | 5 fail |
| S14 no factor, low counts | 1 none | 500 | 0.00 | 0.052 (0.0002) | 0.043 (0.0005) | 0.000 [0.000, 0.007] | 0.050 | 0.788 (0.0011) | within every limit |

**S4, lipidomics and metabolomics (all on the laptop)**

| Scenario | Design | Scored runs | Mean k | Type I at 0.05 (SE) | Mean FDP (SE) | Share FDP > 0.10 [97.5% bounds] | Rule 3 limit | Power (SE) | Rules not passed |
|---|---|---|---|---|---|---|---|---|---|
| S1 3 per group (**void**) | 2 ceiling | 500 | n/a | 0.055 (0.0007) | 0.126 (0.0052) | 0.554 [0.509, 0.598] | 0.604 | 0.153 (0.0038) | 1 inconclusive; 2 fail |
| S1 3 per group (**void**) | **7 recommended** | 500 | 1.47 | 0.040 (0.0010) | 0.101 (0.0069) | 0.370 [0.328, 0.414] | 0.604 | 0.084 (0.0040) | 2 fail |
| S1 3 per group (**void**) | 8 known controls | 500 | 1.47 | 0.040 (0.0011) | 0.098 (0.0075) | 0.336 [0.295, 0.379] | 0.604 | 0.067 (0.0035) | 2 fail |
| S1 3 per group (**void**) | 5 first-pass nulls | 500 | 0.78 | 0.027 (0.0010) | 0.070 (0.0071) | 0.228 [0.192, 0.267] | 0.604 | 0.045 (0.0035) | 2 inconclusive |
| S1 3 per group (**void**) | 3 legacy | 500 | 1.85 | 0.265 (0.0096) | 0.546 (0.0131) | 0.830 [0.794, 0.862] | 0.604 | 0.285 (0.0090) | 1 fail; 2 fail; 3 fail |
| S1 3 per group (**void**) | 1 none | 500 | n/a | 0.014 (0.0005) | 0.021 (0.0062) | 0.022 [0.011, 0.039] | 0.604 | 0.000 (0.0001) | 4 fail |
| S2 40% of features changed | 2 ceiling | 500 | n/a | 0.046 (0.0007) | 0.046 (0.0012) | 0.038 [0.023, 0.059] | 0.088 | 0.634 (0.0029) | within every limit |
| S2 40% of features changed | **7 recommended** | 500 | 1.97 | 0.053 (0.0009) | 0.055 (0.0015) | 0.106 [0.080, 0.136] | 0.088 | 0.563 (0.0035) | 3 inconclusive |
| S2 40% of features changed | 8 known controls | 500 | 1.97 | 0.043 (0.0008) | 0.041 (0.0013) | 0.032 [0.018, 0.051] | 0.088 | 0.508 (0.0042) | within every limit |
| S2 40% of features changed | 5 first-pass nulls | 500 | 3.94 | 0.058 (0.0009) | 0.066 (0.0015) | 0.172 [0.140, 0.208] | 0.088 | 0.571 (0.0031) | 1 inconclusive; 3 fail |
| S2 40% of features changed | 3 legacy | 500 | 2.72 | 0.203 (0.0057) | 0.318 (0.0074) | 0.866 [0.833, 0.895] | 0.088 | 0.427 (0.0077) | 1 fail; 2 fail; 3 fail |
| S2 40% of features changed | 1 none | 500 | n/a | 0.011 (0.0004) | 0.005 (0.0023) | 0.014 [0.006, 0.029] | 0.088 | 0.004 (0.0004) | 4 fail |
| S3 factor tied to groups, 0.6 (ceiling inconclusive) | 2 ceiling | 500 | n/a | 0.049 (0.0007) | 0.072 (0.0036) | 0.312 [0.272, 0.355] | 0.362 | 0.253 (0.0043) | 2 inconclusive |
| S3 factor tied to groups, 0.6 (ceiling inconclusive) | **7 recommended** | 500 | 1.94 | 0.071 (0.0026) | 0.130 (0.0080) | 0.446 [0.402, 0.491] | 0.362 | 0.173 (0.0050) | 1 fail; 2 fail; 3 fail |
| S3 factor tied to groups, 0.6 (ceiling inconclusive) | 8 known controls | 500 | 1.94 | 0.063 (0.0027) | 0.120 (0.0090) | 0.344 [0.302, 0.387] | 0.362 | 0.106 (0.0048) | 1 fail; 2 fail; 3 inconclusive |
| S3 factor tied to groups, 0.6 (ceiling inconclusive) | 5 first-pass nulls | 500 | 3.08 | 0.219 (0.0023) | 0.537 (0.0047) | 1.000 [0.993, 1.000] | 0.362 | 0.450 (0.0042) | 1 fail; 2 fail; 3 fail |
| S3 factor tied to groups, 0.6 (ceiling inconclusive) | 3 legacy | 500 | 2.16 | 0.294 (0.0046) | 0.633 (0.0062) | 0.998 [0.989, 1.000] | 0.362 | 0.495 (0.0063) | 1 fail; 2 fail; 3 fail |
| S3 factor tied to groups, 0.6 (ceiling inconclusive) | 1 none | 500 | n/a | 0.252 (0.0014) | 0.642 (0.0039) | 1.000 [0.993, 1.000] | 0.362 | 0.338 (0.0039) | 1 fail; 2 fail; 3 fail; 4 fail |
| S4 three factors | 2 ceiling | 500 | n/a | 0.047 (0.0007) | 0.055 (0.0019) | 0.130 [0.102, 0.163] | 0.180 | 0.505 (0.0044) | within every limit |
| S4 three factors | **7 recommended** | 500 | 1.93 | 0.036 (0.0006) | 0.037 (0.0021) | 0.110 [0.084, 0.141] | 0.180 | 0.311 (0.0049) | within every limit |
| S4 three factors | 8 known controls | 500 | 1.93 | 0.031 (0.0006) | 0.030 (0.0023) | 0.098 [0.073, 0.127] | 0.180 | 0.227 (0.0049) | within every limit |
| S4 three factors | 5 first-pass nulls | 500 | 4.74 | 0.060 (0.0008) | 0.100 (0.0026) | 0.460 [0.416, 0.505] | 0.180 | 0.491 (0.0039) | 1 fail; 2 fail; 3 fail |
| S4 three factors | 3 legacy | 500 | 3.38 | 0.179 (0.0042) | 0.456 (0.0089) | 0.942 [0.918, 0.961] | 0.180 | 0.415 (0.0080) | 1 fail; 2 fail; 3 fail |
| S4 three factors | 1 none | 500 | n/a | 0.010 (0.0003) | 0.010 (0.0042) | 0.012 [0.004, 0.026] | 0.180 | 0.002 (0.0003) | 4 fail |
| S5 no missing values | 2 ceiling | 500 | n/a | 0.050 (0.0006) | 0.059 (0.0017) | 0.144 [0.114, 0.178] | 0.194 | 0.719 (0.0035) | within every limit |
| S5 no missing values | **7 recommended** | 500 | 1.98 | 0.053 (0.0007) | 0.066 (0.0018) | 0.182 [0.149, 0.219] | 0.194 | 0.704 (0.0039) | 3 inconclusive |
| S5 no missing values | 8 known controls | 500 | 1.98 | 0.048 (0.0008) | 0.063 (0.0023) | 0.188 [0.155, 0.225] | 0.194 | 0.608 (0.0055) | 3 inconclusive |
| S5 no missing values | 5 first-pass nulls | 500 | 4.02 | 0.062 (0.0007) | 0.091 (0.0020) | 0.394 [0.351, 0.438] | 0.194 | 0.768 (0.0033) | 1 fail; 2 fail; 3 fail |
| S5 no missing values | 3 legacy | 500 | 3.27 | 0.165 (0.0041) | 0.381 (0.0083) | 0.932 [0.906, 0.952] | 0.194 | 0.652 (0.0081) | 1 fail; 2 fail; 3 fail |
| S5 no missing values | 1 none | 500 | n/a | 0.003 (0.0002) | 0.002 (0.0020) | 0.002 [0.000, 0.011] | 0.194 | 0.003 (0.0004) | 4 fail |
| S6 50% missing by label (realised 28%) | 2 ceiling | 500 | n/a | 0.045 (0.0006) | 0.056 (0.0023) | 0.186 [0.153, 0.223] | 0.236 | 0.350 (0.0041) | within every limit |
| S6 50% missing by label (realised 28%) | **7 recommended** | 500 | 1.92 | 0.045 (0.0007) | 0.066 (0.0029) | 0.276 [0.237, 0.317] | 0.236 | 0.272 (0.0042) | 3 fail |
| S6 50% missing by label (realised 28%) | 8 known controls | 500 | 1.92 | 0.040 (0.0007) | 0.054 (0.0029) | 0.212 [0.177, 0.250] | 0.236 | 0.243 (0.0042) | 3 inconclusive |
| S6 50% missing by label (realised 28%) | 5 first-pass nulls | 500 | 4.77 | 0.064 (0.0010) | 0.124 (0.0037) | 0.570 [0.525, 0.614] | 0.236 | 0.342 (0.0035) | 1 fail; 2 fail; 3 fail |
| S6 50% missing by label (realised 28%) | 3 legacy | 500 | 3.46 | 0.231 (0.0050) | 0.565 (0.0076) | 0.976 [0.958, 0.988] | 0.236 | 0.336 (0.0066) | 1 fail; 2 fail; 3 fail |
| S6 50% missing by label (realised 28%) | 1 none | 500 | n/a | 0.018 (0.0005) | 0.012 (0.0044) | 0.018 [0.008, 0.034] | 0.236 | 0.002 (0.0003) | 4 fail |
| S7 heavy tails | 2 ceiling | 500 | n/a | 0.042 (0.0006) | 0.041 (0.0015) | 0.054 [0.036, 0.078] | 0.104 | 0.621 (0.0037) | within every limit |
| S7 heavy tails | **7 recommended** | 500 | 1.97 | 0.046 (0.0008) | 0.046 (0.0019) | 0.108 [0.082, 0.139] | 0.104 | 0.543 (0.0048) | 3 inconclusive |
| S7 heavy tails | 8 known controls | 500 | 1.97 | 0.040 (0.0007) | 0.038 (0.0018) | 0.076 [0.054, 0.103] | 0.104 | 0.483 (0.0054) | within every limit |
| S7 heavy tails | 5 first-pass nulls | 500 | 4.37 | 0.058 (0.0008) | 0.080 (0.0024) | 0.310 [0.270, 0.353] | 0.104 | 0.568 (0.0036) | 1 inconclusive; 2 inconclusive; 3 fail |
| S7 heavy tails | 3 legacy | 500 | 3.27 | 0.190 (0.0047) | 0.453 (0.0088) | 0.948 [0.925, 0.966] | 0.104 | 0.478 (0.0083) | 1 fail; 2 fail; 3 fail |
| S7 heavy tails | 1 none | 500 | n/a | 0.011 (0.0004) | 0.009 (0.0036) | 0.014 [0.006, 0.029] | 0.104 | 0.003 (0.0005) | 4 fail |
| S8 known unbalanced batch | 2 ceiling | 500 | n/a | 0.052 (0.0007) | 0.062 (0.0020) | 0.184 [0.151, 0.221] | 0.234 | 0.551 (0.0039) | within every limit |
| S8 known unbalanced batch | **7 recommended** | 500 | 1.00 | 0.094 (0.0008) | 0.139 (0.0032) | 0.680 [0.637, 0.721] | 0.234 | 0.379 (0.0037) | 1 fail; 2 fail; 3 fail |
| S8 known unbalanced batch | 8 known controls | 500 | 1.00 | 0.096 (0.0010) | 0.155 (0.0036) | 0.718 [0.676, 0.757] | 0.234 | 0.361 (0.0040) | 1 fail; 2 fail; 3 fail |
| S8 known unbalanced batch | 5 first-pass nulls | 500 | 3.87 | 0.103 (0.0010) | 0.204 (0.0036) | 0.908 [0.879, 0.932] | 0.234 | 0.461 (0.0036) | 1 fail; 2 fail; 3 fail |
| S8 known unbalanced batch | 3 legacy | 500 | 2.25 | 0.195 (0.0032) | 0.460 (0.0065) | 0.998 [0.989, 1.000] | 0.234 | 0.507 (0.0059) | 1 fail; 2 fail; 3 fail |
| S8 known unbalanced batch | 1 none | 500 | n/a | 0.048 (0.0006) | 0.036 (0.0046) | 0.132 [0.104, 0.165] | 0.234 | 0.047 (0.0020) | 4 fail |
| S9 80% of effects up | 2 ceiling | 500 | n/a | 0.046 (0.0006) | 0.051 (0.0018) | 0.112 [0.086, 0.143] | 0.162 | 0.520 (0.0040) | within every limit |
| S9 80% of effects up | **7 recommended** | 500 | 1.96 | 0.049 (0.0007) | 0.059 (0.0022) | 0.178 [0.145, 0.214] | 0.162 | 0.464 (0.0046) | 3 inconclusive |
| S9 80% of effects up | 8 known controls | 500 | 1.96 | 0.042 (0.0007) | 0.051 (0.0022) | 0.146 [0.116, 0.180] | 0.162 | 0.391 (0.0052) | 3 inconclusive |
| S9 80% of effects up | 5 first-pass nulls | 500 | 4.16 | 0.059 (0.0008) | 0.096 (0.0026) | 0.434 [0.390, 0.479] | 0.162 | 0.522 (0.0041) | 1 fail; 2 fail; 3 fail |
| S9 80% of effects up | 3 legacy | 500 | 2.49 | 0.138 (0.0041) | 0.356 (0.0099) | 0.828 [0.792, 0.860] | 0.162 | 0.352 (0.0081) | 1 fail; 2 fail; 3 fail |
| S9 80% of effects up | 1 none | 500 | n/a | 0.011 (0.0004) | 0.015 (0.0053) | 0.016 [0.007, 0.031] | 0.162 | 0.003 (0.0003) | 4 fail |
| S10 loadings tied to abundance | 2 ceiling | 500 | n/a | 0.046 (0.0006) | 0.054 (0.0019) | 0.148 [0.118, 0.182] | 0.198 | 0.499 (0.0044) | within every limit |
| S10 loadings tied to abundance | **7 recommended** | 500 | 1.68 | 0.041 (0.0007) | 0.049 (0.0023) | 0.158 [0.127, 0.193] | 0.198 | 0.370 (0.0069) | within every limit |
| S10 loadings tied to abundance | 8 known controls | 500 | 1.68 | 0.038 (0.0007) | 0.045 (0.0026) | 0.140 [0.111, 0.174] | 0.198 | 0.309 (0.0065) | within every limit |
| S10 loadings tied to abundance | 5 first-pass nulls | 500 | 4.13 | 0.061 (0.0008) | 0.095 (0.0025) | 0.430 [0.386, 0.475] | 0.198 | 0.514 (0.0040) | 1 fail; 2 fail; 3 fail |
| S10 loadings tied to abundance | 3 legacy | 500 | 3.40 | 0.179 (0.0045) | 0.430 (0.0088) | 0.940 [0.915, 0.959] | 0.198 | 0.440 (0.0076) | 1 fail; 2 fail; 3 fail |
| S10 loadings tied to abundance | 1 none | 500 | n/a | 0.014 (0.0004) | 0.013 (0.0048) | 0.016 [0.007, 0.031] | 0.198 | 0.007 (0.0006) | 4 fail |
| S11 3 per group, factor tied 0.3 (**void**) | 2 ceiling | 500 | n/a | 0.053 (0.0007) | 0.133 (0.0067) | 0.518 [0.473, 0.563] | 0.568 | 0.125 (0.0035) | 2 fail |
| S11 3 per group, factor tied 0.3 (**void**) | **7 recommended** | 500 | 1.48 | 0.053 (0.0014) | 0.148 (0.0087) | 0.466 [0.422, 0.511] | 0.568 | 0.083 (0.0036) | 1 inconclusive; 2 fail |
| S11 3 per group, factor tied 0.3 (**void**) | 8 known controls | 500 | 1.48 | 0.050 (0.0014) | 0.146 (0.0098) | 0.408 [0.365, 0.453] | 0.568 | 0.061 (0.0032) | 1 inconclusive; 2 fail |
| S11 3 per group, factor tied 0.3 (**void**) | 5 first-pass nulls | 500 | 0.84 | 0.043 (0.0018) | 0.144 (0.0094) | 0.376 [0.333, 0.420] | 0.568 | 0.065 (0.0042) | 1 inconclusive; 2 fail |
| S11 3 per group, factor tied 0.3 (**void**) | 3 legacy | 500 | 1.82 | 0.268 (0.0084) | 0.571 (0.0124) | 0.854 [0.820, 0.884] | 0.568 | 0.300 (0.0092) | 1 fail; 2 fail; 3 fail |
| S11 3 per group, factor tied 0.3 (**void**) | 1 none | 500 | n/a | 0.016 (0.0005) | 0.013 (0.0050) | 0.014 [0.006, 0.029] | 0.568 | 0.000 (0.0001) | 4 fail |
| S12 3 per group, three factors (**void**) | 2 ceiling | 500 | n/a | 0.055 (0.0007) | 0.133 (0.0059) | 0.542 [0.497, 0.586] | 0.592 | 0.148 (0.0040) | 1 inconclusive; 2 fail |
| S12 3 per group, three factors (**void**) | **7 recommended** | 500 | 1.35 | 0.027 (0.0007) | 0.079 (0.0085) | 0.212 [0.177, 0.250] | 0.592 | 0.035 (0.0025) | 2 inconclusive |
| S12 3 per group, three factors (**void**) | 8 known controls | 500 | 1.35 | 0.025 (0.0008) | 0.077 (0.0092) | 0.172 [0.140, 0.208] | 0.592 | 0.024 (0.0019) | 2 inconclusive |
| S12 3 per group, three factors (**void**) | 5 first-pass nulls | 500 | 0.75 | 0.020 (0.0008) | 0.052 (0.0062) | 0.158 [0.127, 0.193] | 0.592 | 0.026 (0.0023) | within every limit |
| S12 3 per group, three factors (**void**) | 3 legacy | 500 | 1.96 | 0.271 (0.0100) | 0.549 (0.0128) | 0.828 [0.792, 0.860] | 0.592 | 0.292 (0.0092) | 1 fail; 2 fail; 3 fail |
| S12 3 per group, three factors (**void**) | 1 none | 500 | n/a | 0.010 (0.0004) | 0.010 (0.0045) | 0.010 [0.003, 0.023] | 0.592 | 0.000 (0.0000) | 4 fail |
| S13 no factor, 3 per group (**void**) | 2 ceiling | 500 | n/a | 0.053 (0.0007) | 0.094 (0.0048) | 0.418 [0.374, 0.463] | 0.468 | 0.161 (0.0044) | 2 fail |
| S13 no factor, 3 per group (**void**) | **7 recommended** | 500 | 0.13 | 0.055 (0.0007) | 0.105 (0.0047) | 0.466 [0.422, 0.511] | 0.468 | 0.176 (0.0046) | 2 fail; 3 inconclusive |
| S13 no factor, 3 per group (**void**) | 8 known controls | 500 | 0.13 | 0.054 (0.0007) | 0.097 (0.0047) | 0.446 [0.402, 0.491] | 0.468 | 0.162 (0.0044) | 2 fail; 3 inconclusive |
| S13 no factor, 3 per group (**void**) | 5 first-pass nulls | 500 | 2.03 | 0.074 (0.0009) | 0.170 (0.0042) | 0.754 [0.714, 0.791] | 0.468 | 0.302 (0.0051) | 1 fail; 2 fail; 3 fail |
| S13 no factor, 3 per group (**void**) | 3 legacy | 500 | 1.69 | 0.167 (0.0028) | 0.405 (0.0069) | 0.970 [0.951, 0.983] | 0.468 | 0.454 (0.0065) | 1 fail; 2 fail; 3 fail |
| S13 no factor, 3 per group (**void**) | 1 none | 500 | n/a | 0.053 (0.0007) | 0.094 (0.0048) | 0.418 [0.374, 0.463] | 0.468 | 0.161 (0.0044) | 2 fail |
| S14 no factor, 50% missing by label (realised 16%) | 2 ceiling | 500 | n/a | 0.047 (0.0006) | 0.054 (0.0019) | 0.152 [0.122, 0.187] | 0.202 | 0.474 (0.0038) | within every limit |
| S14 no factor, 50% missing by label (realised 16%) | **7 recommended** | 500 | 0.15 | 0.049 (0.0006) | 0.056 (0.0021) | 0.162 [0.131, 0.197] | 0.202 | 0.478 (0.0038) | within every limit |
| S14 no factor, 50% missing by label (realised 16%) | 8 known controls | 500 | 0.15 | 0.048 (0.0006) | 0.055 (0.0020) | 0.158 [0.127, 0.193] | 0.202 | 0.471 (0.0039) | within every limit |
| S14 no factor, 50% missing by label (realised 16%) | 5 first-pass nulls | 500 | 2.12 | 0.059 (0.0007) | 0.084 (0.0023) | 0.344 [0.302, 0.387] | 0.202 | 0.523 (0.0035) | 1 inconclusive; 2 fail; 3 fail |
| S14 no factor, 50% missing by label (realised 16%) | 3 legacy | 500 | 1.88 | 0.117 (0.0031) | 0.272 (0.0091) | 0.770 [0.731, 0.806] | 0.202 | 0.419 (0.0055) | 1 fail; 2 fail; 3 fail; 5 fail |
| S14 no factor, 50% missing by label (realised 16%) | 1 none | 500 | n/a | 0.047 (0.0006) | 0.054 (0.0019) | 0.152 [0.122, 0.187] | 0.202 | 0.474 (0.0038) | within every limit |

### Where the ceiling itself fails (these scenarios cannot fault a design)

| Omic | Scenario | Ceiling: type I at 0.05 (SE), at 0.01, mean FDP (SE) | Reading |
|---|---|---|---|
| Classic proteomics | S12, 3 per group with three factors (runner) | 0.078 (0.0009), 0.027, 0.127 (0.0049) | Void, as the pilot had shown |
| DIA_LIMPA | S14, no factor, 50% missing (runner) | 0.064 (0.0004), 0.016, 0.074 (0.0011) | Void by rule 1. There is no factor here, so the ceiling is the no-correction result: the route itself is liberal at 50% missing. With a factor and 50% missing (S6) the ceiling passes but sits high (0.063, mean FDP 0.068, SE 0.0013) |
| Transcriptomics | S7, heavy tails | 0.085 (0.0003), 0.010, 0.022 (0.0004) | Void by rule 1 (the mean FDP is low) |
| Transcriptomics | S9, 80% of effects up | 0.075 (0.0004), 0.018, 0.074 (0.0006) | Void by rule 1 |
| Lipid/metab | S1, S11, S12: 3 per group with a factor | 0.053 to 0.055, 0.014 to 0.015, 0.126 to 0.133 (0.0052 to 0.0067) | Void by rule 2 |
| Lipid/metab | S13: 3 per group, no unwanted variation | 0.053 (0.0007), 0.012, 0.094 (0.0048) | Void by rule 2; the ceiling is the no-correction result here. Whether this is the two-assay scoring unit described for C4 above or the test at 3 per group was not separated |
| Lipid/metab | S3, factor tied to groups 0.6 | 0.049 (0.0007), 0.011, 0.072 (0.0036) | Not void: rule 2 inconclusive, flagged |

In DIA_LIMPA with no missing values (S5) the ceiling is calibrated on the 248 runs that completed, and half the runs crash in every arm: the scenario is a condition of the route, not of any RUV design.

### Conditions of use per omic (measured; mean FDP with its standard error unless stated)

"Within every limit" means rules 1 to 3 pass and rule 4 (or 5), where it applies, is not failed. "Tail" is rule 3: the share of data sets above 0.10 against the ceiling's share plus 0.05. Void scenarios are left out of each row and listed in the last column.

**Classic proteomics** (13 scenarios judged; S9 to S14 on the runner)

| Design | Within every limit | Inconclusive on the tail only (mean FDP within its limit) | Outside a limit | Void |
|---|---|---|---|---|
| Recommended (7) | 3 per group (0.051, SE 0.0011; power 0.660 against the ceiling's 0.661); 40% changed (0.057); no missing values; heavy tails; 80% up (0.056); abundance-tied loadings (0.053); 3 per group with the factor tied 0.3 (0.055); no factor at 3 per group and at 50% missing (0.053, 0.053) | Factor tied 0.6 (0.062, SE 0.0012; share 0.076 [0.054, 0.103] against 0.058); three factors (0.057; 0.054 against 0.060); 50% missing (0.062, SE 0.0016; 0.124 against 0.104); known unbalanced batch (0.059; 0.056 against 0.070) | none | S12 |
| Known controls (8) | every scenario but one (0.035 to 0.054) | Factor tied 0.6 (0.058, SE 0.0013; share 0.070 against 0.058) | none | S12 |
| First-pass nulls (5) | 3 per group; 40% changed; no and 50% missing; heavy tails; 80% up; abundance-tied loadings (0.037 to 0.056). Three factors: calibrated (0.003) with power 0.071 against the ceiling's 0.863 (mean k 1.20) | none | Factor tied 0.6 (type I 0.146, mean FDP 0.273, 99% of data sets above 0.10); known batch (tail 0.214 against 0.070; power 0.759 against 0.872); 3 per group with the factor tied 0.3 (tail 0.128 against 0.086); no factor at 3 per group (tail 0.144 against 0.076; power 0.617 against 0.746 without correction); no factor with 50% missing (power 0.870 against 0.908) | S12 |
| Legacy (3) | Three factors only (0.016; power 0.143) | none | Every other scenario: mean FDP 0.082 to 0.684 (laptop cells), 0.086 to 0.521 † (runner cells) | S12 |
| None (1) | Calibrated in 12 scenarios (mean FDP at most 0.009 with a factor); power 0.000 at 3 per group with a factor and with three factors, 0.144 to 0.623 in the other scenarios with a factor | none | Factor tied 0.6: type I 0.258 (SE 0.0008), mean FDP 0.593, every data set above 0.10 | S12 |

**DIA_LIMPA** (13 scenarios judged; all on the runner; engine `fastlimpadpc`)

| Design | Within every limit | Inconclusive | Outside a limit | Void |
|---|---|---|---|---|
| Recommended (7) | 40% changed (0.051, SE 0.0006); heavy tails (0.051); known unbalanced batch (0.055, SE 0.0019; mean k 1.98); 80% up (0.050); abundance-tied loadings (0.052); no factor at 3 per group (0.048); no missing values on the 248 runs that completed (0.050) | 50% missing (type I 0.066, 0.0152 at 0.01; mean FDP 0.071, SE 0.0014; the ceiling 0.068) | 3 per group (tail 0.132 [0.104, 0.165] against 0.070; mean FDP 0.070, SE 0.0031, against the ceiling's 0.045); factor tied 0.6 (tail 0.084 [0.061, 0.112] against 0.056; mean FDP 0.065, SE 0.0013, within its limit); three factors (type I 0.066, mean FDP 0.088, SE 0.0042, tail 0.188; mean k 2.63 where 3 exist); 3 per group with the factor tied 0.3 (tail 0.164 against 0.084; mean FDP 0.078, SE 0.0036); 3 per group with three factors (mean FDP 0.164, SE 0.0086; mean k 1.39) | S14 |
| Known controls (8) | 40% changed; 80% up; abundance-tied loadings; no factor at 3 per group (0.047 to 0.052); no missing values on the 248 runs that completed (0.050) | Heavy tails and known batch (tail only: 0.040 and 0.058 against 0.052); 50% missing (0.073) | The same five scenarios as the recommended design, each slightly further out: 3 per group (tail 0.156); factor tied 0.6 (0.180); three factors (mean FDP 0.088); 3 per group with the factor tied (0.180); 3 per group with three factors (mean FDP 0.165) | S14 |
| First-pass nulls (5) | No factor at 3 per group (0.056) | none | Every scenario with a factor: mean FDP 0.088 to 0.301, 26% to 89% of data sets above 0.10 | S14 |
| Legacy (3) † | none | none | Every scenario: mean FDP 0.248 to 0.677 † | S14 |
| None (1) | No factor at 3 per group (0.047) | none | The tail in every scenario with a factor (9% to 84% of data sets above 0.10; mean FDP 0.040 to 0.142); factor tied 0.6: type I 0.324 (SE 0.0082), mean FDP 0.533 | S14 |

**Transcriptomics** (11 scenarios judged; S14 on the runner; no 50%-missing scenario; S5 and S14 are low-count scenarios)

| Design | Within every limit | Inconclusive | Outside a limit | Void |
|---|---|---|---|---|
| Recommended (7) | 3 per group (0.046, SE 0.0005; power 0.369, equal to the ceiling); 40% changed (0.035; power 0.835); factor tied 0.6 (0.052, SE 0.0005; power 0.575 against 0.553); three factors (0.049; power 0.762 against 0.763); low counts (0.047); known batch (0.049); abundance-tied loadings (0.048); 3 per group with the factor tied 0.3 (0.049); no factor at 3 per group (0.044; power 0.380 against 0.386 without correction); no factor with low counts (0.046, SE 0.0005; power 0.784 against 0.788; mean k 2.95 where no factor exists) | none | 3 per group with three factors: no calls at all (type I 0.005, power 0.000 against the ceiling's 0.341, SE 0.0017; mean k 2.00): calibrated, no gain over no correction (rule 4) | S7, S9 |
| Known controls (8) | The same scenarios with a factor as the recommended design, power 0.008 to 0.044 lower (0.035 to 0.049) | none | 3 per group with three factors (power 0.000); no factor: power lost against no correction, 0.310 against 0.386 at 3 per group and 0.569 (SE 0.0091) against 0.788 with low counts (rule 5) | S7, S9 |
| First-pass nulls (5) | 3 per group (power 0.250 against 0.369); three factors (0.180 against 0.763); abundance-tied loadings (0.623 against 0.774); 3 per group with three factors (0.123); no factor with low counts (0.773 against 0.788) | 40% changed (power 0.092 against 0.835); low counts (power 0.440 against 0.754) | Factor tied 0.6 (tail 0.076 against 0.050; power 0.002); known batch (tail 0.184; power 0.042, below no correction's 0.083); 3 per group with the factor tied 0.3 (tail 0.080 against 0.050; mean FDP 0.065, SE 0.0029); no factor at 3 per group (power 0.359 against 0.386) | S7, S9 |
| Legacy (3) | Calibrated (0.008 to 0.050) in 3 per group, 40% changed, three factors, low counts, abundance-tied loadings, 3 per group with three factors; power below the ceiling in each (0.094 to 0.793) | none | Factor tied 0.6 (mean FDP 0.724); known batch (0.433); 3 per group with the factor tied 0.3 (0.186); no factor: power 0.324 against 0.386 at 3 per group and 0.749 † against 0.788 with low counts | S7, S9 |
| None (1) | Calibrated in 10 scenarios; power at most 0.114 in every scenario with a factor other than the one at right | none | Factor tied 0.6: type I 0.337 (SE 0.0010), mean FDP 0.450, every data set above 0.10 | S7, S9 |

**Lipidomics and metabolomics** (10 scenarios judged; all on the laptop; "50% missing" is the detection limit's label, realised 28% in S6 and 16% in S14)

| Design | Within every limit | Inconclusive on the tail only | Outside a limit | Void |
|---|---|---|---|---|
| Recommended (7) | Three factors (0.037, SE 0.0021; power 0.311 against the ceiling's 0.505); abundance-tied loadings (0.049; power 0.370 against 0.499); no factor with 50% missing (0.056; power 0.478 against 0.474) | 40% changed (0.055, SE 0.0015; tail 0.106 against 0.088); no missing values (0.066, SE 0.0018; 0.182 against 0.194); heavy tails (0.046; 0.108 against 0.104); 80% up (0.059; 0.178 against 0.162) | Factor tied 0.6 (type I 0.071, SE 0.0026; mean FDP 0.130, SE 0.0080; the ceiling 0.072, itself inconclusive); known unbalanced batch (type I 0.094, mean FDP 0.139, SE 0.0032, 68% of data sets above 0.10; mean k 1.00; the ceiling 0.062); 50% missing (tail only: 0.276 [0.237, 0.317] against 0.236, the lower bound just outside; mean FDP 0.066, SE 0.0029, within its limit) | S1, S11, S12, S13 |
| Internal standards as controls (8) | 40% changed (0.041); three factors (0.030; power 0.227); heavy tails (0.038); abundance-tied loadings (0.045); no factor with 50% missing (0.055). Power below the recommended design's in every scenario with a factor | No missing values (0.063); 50% missing (0.054); 80% up (0.051) | Factor tied 0.6 (mean FDP 0.120, SE 0.0090); known unbalanced batch (0.155, SE 0.0036) | S1, S11, S12, S13 |
| First-pass nulls (5) | none | none | The tail in every scenario (17% to 100% of data sets above 0.10); mean FDP 0.066 (40% changed) to 0.537 (factor tied 0.6) | S1, S11, S12, S13 |
| Legacy (3) | none | none | Every scenario: mean FDP 0.272 to 0.633 | S1, S11, S12, S13 |
| None (1) | Calibrated in 9 scenarios (mean FDP at most 0.054); power at most 0.007 under drift, 0.047 with the known batch | none | Factor tied 0.6: type I 0.252 (SE 0.0014), mean FDP 0.642, every data set above 0.10 | S1, S11, S12, S13 |

In the known-batch scenario of lipid/metab, putting the known batch in the model without RUV (arm 17) gives a mean FDP of 0.242 (SE 0.0042); in DIA_LIMPA 0.064 (SE 0.0067; tail 0.134 against 0.052); in classic proteomics 0.010 and in transcriptomics 0.048, both within every limit.

### Plain reading

- **The recommended design never left a limit in classic proteomics, and in transcriptomics only at 3 per group with three factors, where it makes no calls.** In classic proteomics four scenarios are inconclusive on the tail with the mean FDP at 0.057 to 0.062.
- **It stops holding in DIA_LIMPA** at 3 per group (the tail; mean FDP 0.070 to 0.078 with and without a tied factor), with three factors (0.088; the estimated count returns 2.63 on average) and with both (0.164), and on the tail with the factor tied 0.6 (mean FDP 0.065). **In lipid/metab** it stops with the factor tied 0.6 (0.130) and with a known unbalanced batch (0.139; the estimated count returns 1), and misses the tail narrowly at the higher missing share.
- **The number of factors is the weak point where there are several or few samples.** k = 1 loses all power with three factors (arm 11: power 0.003 in classic proteomics, 0.000 in transcriptomics); the estimated count returns 3.00 with three factors at 6 per group in classic proteomics and transcriptomics, is low in DIA_LIMPA (2.63 of 3; 1.39 of 3 at 3 per group) and in transcriptomics at 3 per group (2.00 of 3, no calls), and returns 2.95 where no factor exists in the low-count transcriptomics scenario without cost there (power 0.784 against 0.788).
- **No correction is not safe when the unwanted factor follows the groups strongly:** at 0.6 its mean FDP is 0.593 (classic proteomics), 0.533 (DIA_LIMPA), 0.450 (transcriptomics) and 0.642 (lipid/metab). The recommended design gives 0.062, 0.065, 0.052 and 0.130 in the same scenarios.
- **Known controls** behave as the recommended design does, with a little less power; in transcriptomics they lose power when there is no factor (0.569 against 0.788 with low counts). **First-pass nulls** hold in classic proteomics except under a tied factor, a known batch or no factor, and lose nearly all power there with three factors (0.071); they fail every DIA_LIMPA scenario with a factor and every lipid/metab scenario, and in transcriptomics they are calibrated on the mean but lose most of their power in most scenarios. These first-pass figures are arm 5's (the separation rule for k); the estimated-count variant (arm 23) was not run in the extra-hard cells.
- **Legacy** is outside the limits in nearly every scenario of classic proteomics, DIA_LIMPA and lipid/metab (mean FDP 0.082 to 0.684), and in transcriptomics where the factor is tied to the groups or a batch is present (0.186 to 0.724). Runner figures marked †.
- **Not measurable here:** 3 per group in lipid/metab (the ceiling fails in all four scenarios), 50% missing without a factor in DIA_LIMPA (the route is liberal on its own), heavy tails and one-directional effects in transcriptomics, and 3 per group with three factors in classic proteomics.
