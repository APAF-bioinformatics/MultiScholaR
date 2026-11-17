# KIT-124: RUV in every omic, confirming studies (DRAFT, 2026-11-17)

**Draft.** S1, S3 and S4 core cells are scored; S2 (DIA_LIMPA) core cells are running (4 of 8 at 16:00) and every simulator's stress cells follow. Nothing here is the user's decision yet; the gate table at the end is filled in when S2 is scored. Rules: the ticket's W3 table (revision 1.2), applied by `aggregate.R` (`0bb67eb2`) to the pinned tree `0b6fcb6d` (KIT-123 R1 and R2 in). Scores of the core phases: scratchpad `k124-agg/s{1,3,4}-core/` (copied here with the final aggregation).

## The binding question before any per-omic row: rule 4 in the global-null cell (A12)

Rule 4 ("gain over no RUV in cells with a factor") cannot be met literally in C4, where nothing is differentially abundant: there is no power to gain, and the rule-3 reduction has nothing to reduce. Read literally, **every design fails C4, including the true-factor arm**. Revision 1.4 puts both readings to the user:

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
| 14 | `sva` (complete proteins) | fail | inconclusive | C4(4) | 0.113 | 0.95 |
| 15 | technical replicates as RUV structure, subtract (C6 only) | pass | pass | none | 0.012 | n/a |
| 16 | the same, W in the model (C6 only) | pass | pass | none | 0.032 | n/a |

Plain reading: today's design is confirmed miscalibrated. All proteins as controls with the factors in the model misses one bound in one cell (C1, 12 per group with 10% effects: the share of runs above 10% FDP is 0.115 against a limit of 0.057; its mean FDP, 0.066, is within rule 2). No RUV is calibrated but loses about 40% of the power. Where technical replicates exist, arms 15/16 pass.

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

Plain reading: RUVs with all genes as controls and the groups as replicate sets matches the true-factor ceiling in every core cell; under the second reading of rule 4 it is the only practical design not failing anything (arm 11 needs no k estimate). Today's design and RUVr fail when the factor is tied to the groups (C5: FDP 0.30 and 0.52). No RUV is calibrated but finds almost nothing when the unwanted factor is strong. C6 (library replicates) is reported only: every RUV arm has a mean FDP of 0.27 there, the package's replicate handling (parked row 23, repaired by KIT-128, not in this tree).

### S4, lipidomics and metabolomics (complete; C4 void)

C4 is void: the true-factor arm calls something in 13% of global-null runs (limit 0.075), so the cell judges no design (rule as written).

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
| 14 | `sva` | fail | C2, C5, C7 | 0.310 | 0.70 |

Plain reading: no design passes. Today's design is the worst measured (more than 90% of runs above 10% FDP). The closest are all features or the standards as controls with the factor in the model, and pooled QC as the replicate set: each misses rule 3 (the FDP tail) in two cells. No correction is calibrated but has almost no power under drift.

### S2, DIA_LIMPA

Pending (core cells 4 of 8 at 16:00 on 2026-11-17). Engine `fastlimpadpc`; the claim for the default `limpa` engine is stated with KIT-053's exactness record or limited (revision 1.4, M-11).

## What the rule gives, literally (ticket W3, "Recommendation per omic")

| Omic | Literal branch | With rule 4 undecidable in C4 |
|---|---|---|
| Classic proteomics | Eligible only with technical replicates (15 > 16); otherwise no design passes, and arm 1 fails rule 4: "reported as such" | Same (arm 7 still fails C1 rule 3) |
| Transcriptomics | RUV off by default | Arms 7/11 inconclusive (no failure): selectable with a note under KIT-129's verdict-to-option rule |
| Lipidomics, metabolomics | RUV off by default | Same (C4 is void anyway) |
| DIA_LIMPA | pending | pending |

## For the gate (to be completed with S2)

- Both readings of rule 4 in C4, side by side (above).
- The verdict-to-option rule proposed in KIT-129 revision 1.1, for confirmation: a passing design is selectable; the default per omic is the user's choice among passing ones; a failing design is not offered except subtract-then-test as the warned legacy switch; an inconclusive or not-run design is selectable with a note; diagnostic arms (oracle controls) are never options.
- Read literally, that rule leaves most omics with no RUV option at all, and "RUV off by default" means accepting the power losses above (transcriptomics and lipid/metab near zero under strong unwanted variation). This is the user's call; the measures are here to overrule the rule's literal reading if the user chooses to.
- The `sva` choice (k from `sva::num.sv` in the arms that take k = E; `sva` is not a dependency; IR-M6).
- Stress cells: conditions of use (minimum group size, high missingness), added when scored.
