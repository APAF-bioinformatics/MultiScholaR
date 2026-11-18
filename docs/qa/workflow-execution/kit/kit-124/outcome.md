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
