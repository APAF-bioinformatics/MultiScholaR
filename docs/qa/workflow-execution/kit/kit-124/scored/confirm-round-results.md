# KIT-124 confirmation round: results (Part 1 and Part 2, 2026-11-19)

Declaration: `confirm-round-declaration.md` (decisions M29, M34, M38). Package under test: `kit-129-integration` 44b5f2195e355bb6e6efe1dae21905e30b81c26c, installed in every job. Savings L1 and L2 on (proven identical, declaration section 7). Pass rules 1 to 7 unchanged: `aggregate.R` of `sim-k124-exact` 2ba5cf87, sha256 `089d852c1f3b08d70096b1218ef26550f393cffd213e4bdd9980453075fd2a66`, named in `.github/k124-sim/score.json` and checked by the scoring job before any result was read. Scored once on a hosted runner (workflow `KIT-124 scoring`, run 44480610832); tables in `confirm-gha/` (layout of `final-20261118/`: `s1`..`s4/k124-*.tsv`, `inventory.tsv`, `jobs.tsv`, `not-committed.sha256`, `gate/`). Raw tables are not committed (sha256 in `not-committed.sha256`; Part 1's in `confirm-p1-tables.sha256`).

## Runs

| Part | Run | Label | Jobs | Outcome |
|---|---|---|---|---|
| Binding lifecycle | 44434645817 | confirm-lc4 | 6 | as declared; gate passed |
| 1 identity | 44436553563 | confirm-p1 | 40 | every expected identity held except S1 legacy (C6): added to Part 2 by the declared rule |
| 2 core | 44443281402 | confirm-p2c | 72 | 70 succeeded; **S3 C1 and S3 C7 reached the 300-minute limit** (nothing uploaded; one attempt; a limit is not an infrastructure failure, so not rerun) |
| 2 extra-hard | 44443334468 | confirm-p2x | 4 | all succeeded |
| 2 transcriptomics C1, C7 (amendment, section 12) | 44596252326 | confirm-p2t | 8 (ranges of 250 runs) | all succeeded; longest job 80 minutes; 33.8 core-hours; gate passed |

The two cells that timed out (transcriptomics C1, 12 per group; C7, three groups) were rerun in run ranges after the user's approval of 2026-11-20 (declaration section 12; nothing else changed). Transcriptomics was then scored once more with the same workflow and `aggregate.R` (run 44612000970, label `confirm-gha-s3`, from the top-level tables because C1 and C7 come in ranges), **only because those two cells were missing**; classic proteomics, DIA and lipid/metab were not re-scored. The seven earlier transcriptomics cells' rows of `k124-cells.tsv` are identical in the two scores (60 rows). `confirm-gha/s3/` holds the complete score; `confirm-gha/s3-rescore/` its inventory, jobs and raw-table hashes; `confirm-gha/gate/confirm-p2t/` the gate of C1 and C7 (arms 1 and 2 agree with the committed reference in every rule-bearing value).

## Gate (arms 1 and 2 against the committed reference; M19 values)

Every gated cell agrees on every rule-bearing value in both parts (S1 to S4, C0 to C8 where run; S4 S1, S3, S8 against the laptop rows). Reported, not rule-bearing: S3 C6 bias outside the rule (relative 1.8e-6). S3's cell S3 has no reference (declared). The gate jobs ended `failure` because of the **S1 legacy comparison**: arm 3 differs from the `exact-gha` rows in C6 and C8 (same seeds and code; calls, null rates, fdp and power differ), while C0 to C5 and C7 agree. Part 1 had shown C6 only. This is the machine dependence of the legacy correction on the technical-duplicate designs (declaration 7b); legacy figures for C6 and C8 depend on the runner machine.

## Results per omic and option (primary subset `full`; 1,000 runs per core cell, 500 per extra-hard cell)

Columns: verdict over the judged core cells (rules 1 to 7; the C4 reading without rule 4 in brackets where it differs); type I at 0.05 in C0 (and the largest over the judged core cells); mean FDP at the 5 % target in C0 (association 0), C5 (association 0.3) and the extra-hard S3 (0.6) where run; the worst core cell with true effects (the app's figure); power in C0; refused rows. "Was" = the earlier figure from the arm the app quotes now (`final-20261118`, `newarms-gha`, `r0a-gha`) and, where it exists, the exactness round's final-method arm.

### Classic proteomics (S1)

| Option (shipped arm) | Verdict | Type I C0 (max) | FDP C0 / C5 | Worst core FDP | Power C0 | Refused | Was (arm, source): worst FDP, power, verdict |
|---|---|---|---|---|---|---|---|
| Recommended (34) | identity with arm 32 (Part 1): fail | 0.053 (0.058 C1) | 0.055 / 0.058 | 0.067 C1 | 0.916 | 0 | app: arm 7, 0.066 C1, 0.91, fail |
| Technical replicates (38) | **pass** (C6) | C6 0.052, C8 0.050 | C6 0.055, C8 0.052 | 0.055 C6 | C6 0.897, C8 0.904 | 0 | app: arm 16, 0.054 C6, 0.90, pass; arm 33 identical |
| Supplied list (35) | fail (inconclusive without rule 4 in C4) | 0.051 (0.051 C5) | 0.052 / 0.054 | 0.054 C5 | 0.914 | 0 | app: arm 8, 0.053 C5, 0.91, fail (inconclusive) |
| Contaminated list (36) | fail | 0.063 (0.083 C1) | 0.076 / 0.078 | 0.121 C1 | 0.915 | 0 | arm 9: 0.119 C1, 0.913, fail (app note: 11.9 %) |
| Supplied list, separation rule (41) | fail (C3 rule 5, C4 rule 4) | 0.051 (0.055 C6) | 0.053 / 0.054 | 0.060 C6 | 0.906 | 0 | app: arm 22, 0.060 C6, 0.90, fail |
| First-pass nulls (39) | fail (C4 rule 4; C5 rules 2, 3) | 0.051 (0.066 C5) | 0.052 / 0.079 | 0.079 C5 | 0.916 | 0 | app: arm 23, 0.077 C5, 0.91, fail |
| First-pass nulls, separation rule (40) | fail (C4; C5 rules 2, 3) | 0.051 (0.066 C5) | 0.053 / 0.079 | 0.079 C5 | 0.911 | 0 | app: arm 5, 0.078 C5, 0.91, fail |
| Legacy (3) | fail | 0.072 (0.194 C5) | 0.098 / 0.385 | 0.385 C5 | 0.920 | 0 | app: arm 3, 0.385 C5, 0.92, fail (C0 to C5, C7 identical; C6, C8 machine-dependent) |
| No RUV (37, M30's test) | fail (rule 4) | 0.014 (0.050 C3) | 0.007 / 0.033 | 0.052 C3 | 0.493 | 0 | app: arm 1, 0.052 C3, 0.49; arm 27 identical |

### DIA with limpa (S2)

| Option | Verdict | Type I C0 (max) | FDP C0 / C5 | Worst core FDP | Power C0 | Refused | Was |
|---|---|---|---|---|---|---|---|
| Recommended / technical replicates (34) | identity with arm 32 (Part 1) | 0.051 (0.052 C6) | 0.052 / 0.053 | 0.056 C6 | 0.910 | 0 | app: arm 7, 0.055 C6, 0.91 |
| Supplied list (35) | inconclusive (C4, C5) | 0.050 (0.052 C6) | 0.052 / 0.055 | 0.056 C6 | 0.907 | 0 | app: arm 8, 0.056 C6, 0.91, inconclusive |
| Contaminated list (36) | fail | 0.060 (0.068 C1) | 0.072 / 0.074 | 0.090 C1 | 0.905 | 0 | arm 9: 0.089 C1 (app note 8.9 %), fail |
| Supplied list, separation rule (41; now the rule on the list) | inconclusive (C2, C4, C5) | 0.051 (0.053 C5) | 0.054 / 0.057 | 0.057 C5 | 0.904 | 0 | app: arm 22, 0.056 C5, 0.90, inconclusive |
| First-pass nulls (39) | fail | 0.068 (0.091 C5) | 0.092 / 0.147 | 0.147 C5 | 0.918 | 0 | app: arm 23, 0.146 C5, 0.92, fail |
| First-pass nulls, separation rule (40) | fail | 0.070 (0.093 C5) | 0.095 / 0.151 | 0.151 C5 | 0.917 | 0 | app: arm 5, 0.149 C5, 0.92, fail |
| Legacy, factors in the test (4) | fail | 0.170 (0.265 C5) | 0.290 / 0.435 | 0.435 C5 | 0.897 | 0 | app: arm 4, 0.435, 0.90 (identical) |
| Legacy, corrected (3) | fail | 0.194 (0.287 C5) | 0.331 / 0.466 | 0.466 C5 | 0.911 | 0 | app: arm 3, 0.466, 0.91 (identical) |
| No RUV (37) | identity with arm 1 (Part 1) | 0.047 (0.116 C5) | 0.065 / 0.191 | 0.191 C5 | 0.558 | 0 | app: arm 1, 0.191 C5, 0.56 |

### Transcriptomics (S3; every core cell, C0 to C8, and the extra-hard S3)

| Option | Verdict (failing cells) | Type I C0 (max) | FDP C0 / C5 / S3 | Worst core FDP | Power C0 | C1: FDP, power | C7: FDP, power | Was (app's arm): worst FDP, power, C7 |
|---|---|---|---|---|---|---|---|---|
| Recommended (34) | identity with arm 31 (Part 1); extra-hard S3: pass | 0.056 (0.058 C7) | 0.048 / 0.048 / **0.052** | 0.055 C2 | 0.774 | 0.049, 0.956 | 0.051, 0.763 | arm 7: 0.055 C2, 0.77; S3 0.052 |
| Supplied list (35) | fail (C3 rule 5, C4 rule 4) | 0.056 (0.057 C7) | 0.048 / 0.048 / — | 0.056 C2 | 0.763 | 0.049, 0.954 | 0.051, 0.752 | arm 8: 0.056 C2, 0.76 (same cells) |
| Contaminated list (36) | fail (C0, C1, C2, C3, C4, C5, C7) | 0.064 (0.073 C1) | 0.061 / 0.065 / — | 0.076 C1 | 0.757 | 0.076, 0.951 | 0.064, 0.748 | arm 9: 0.076 C1 (app 7.6 %), same failing cells |
| Supplied list, separation rule (41) | fail (C3, C4) | 0.054 (0.056 C7) | 0.047 / 0.046 / — | 0.053 C2 | 0.612 | 0.048, 0.934 | 0.050, 0.666 | arm 22: 0.053 C2, 0.61 (identical) |
| First-pass nulls (39) | fail (C4; C5 rule 1) | 0.056 (0.071 C5) | 0.048 / 0.071 / — | 0.071 C5 | 0.774 | 0.050, 0.956 | 0.051, 0.763 | arm 23: 0.071 C5, 0.77 (identical) |
| First-pass nulls, separation rule (40) | fail (C4; C5 rule 3; C7 rules 3, 4) | 0.055 (0.061 C5) | 0.047 / 0.064 / — | 0.064 C5 | 0.600 | 0.050, 0.827 | 0.064, **0.049** | app: arm 5, 0.064 C5, 0.60, C7 power **0.023** ("2 %"); arm 24 identical to the shipped arm |
| Legacy (3) | fail (C3, C4, C5) | 0.055 (0.179 C5) | 0.047 / 0.299 / 0.724 | 0.299 C5 | 0.729 | 0.051, 0.948 | 0.050, 0.731 | arm 3: identical |
| No RUV (37) | identity with arm 1 (Part 1) | 0.010 (0.053 C3) | 0.000 / 0.003 / 0.450 | 0.044 C3 | 0.005 | | | arm 1: identical |

No transcriptomics verdict changes: every option fails exactly where its earlier arm failed (C1 and C7 now included). One quoted number changes: with three groups (C7) first-pass nulls with the separation rule find 4.9 % of the real differences on the shipped code (the app's note says 2 %, from arm 5, whose k came from the sweep; the shipped rule is arm 24's, the separation rule on the first pass's own controls).

### Lipidomics and metabolomics (S4; one table for both)

| Option | Verdict | Type I C0 (max) | FDP C0 / C5 / S3 | Worst core FDP | Power C0 | Refused | Was |
|---|---|---|---|---|---|---|---|
| Recommended (34) | fail (C0 rule 3; C1, C5 rules 2, 3; C2, C7 rule 3) | 0.054 (0.060 C1) | 0.072 / **0.085** / **0.148** | 0.085 C5 | 0.474 | 0 | app: arm 7, 0.072 C5, S3 0.130, 0.43, fail (C1, C5 rule 3); exact arm 32: 0.088 C5, 0.467 |
| Recommended, extra-hard | S1 0.111 (rule 2 fail), S3 0.148, S8 0.180 | | | | | | app note: S1 10 %, S3 13 %, S8 14 % (arm 7) |
| Technical replicates (38) | not judged (C6, C8 reported) | C6 0.049 | C6 0.063, C8 0.065 | | C6 0.341, C8 0.472 | 0 | app: arm 16 (r0a), C6 0.058, 0.36 |
| Supplied list, standards (35) | fail (C0, C5, C7 rule 3) | 0.049 (0.053 C5) | 0.063 / 0.079 | 0.079 C5 | 0.399 | 0 | app: arm 8, 0.068 C5, 0.37, fail (C5, C7 rule 3) |
| Contaminated standards (36) | fail | 0.121 (0.185 C1) | 0.288 / 0.290 | 0.399 C2 | 0.426 | 0 | arm 9: 0.371 C1 (app note 37 %) |
| Supplied, separation rule (41) | fail | 0.049 (0.053 C3) | 0.081 / 0.095 | 0.095 C5 | 0.264 | **2,327 of 9,000 rows** | app: arm 22, 0.118 C2, 0.26, fail |
| First-pass nulls (39) | fail (C5 rules 1 to 3; C7 rules 2, 3) | 0.050 (0.102 C5) | 0.063 / **0.203** | 0.203 C5 | 0.475 | 0 | app: arm 23, 0.174 C5, 0.44 |
| First-pass nulls, separation rule (40) | fail | 0.065 (0.120 C5) | 0.108 / 0.261 | 0.261 C5 | 0.527 | 510 of 9,000 | app: arm 5, 0.239 C5, 0.50 |
| Legacy (3) | fail | 0.193 (0.207 C1) | 0.470 / 0.472 | 0.557 C2 | 0.423 | 0 | app: arm 3, 0.556 C2, 0.42 |
| No RUV (37) | identity with arm 1 | 0.012 (0.048 C3) | 0.005 / 0.043 | 0.059 C7 | 0.003 | 0 | app: arm 1, 0.059 C7, 0.003 |

**Refusals in lipid/metab (scored as no RUV, R-07).** Arm 41 (standards with the separation rule): 2,729 assay refusals in 2,327 rows, all "k = 5 is not below the 5 controls observed in every sample" (decision M32: 5 standards per assay; the separation rule returned k = 5). Arm 40: 472 "k leaves no residual degree of freedom" and 84 "k above the limit of the biological units" (C6, C8; M28). No refusal in any other shipped arm or omic; no crash anywhere.

## What changed

- **Verdicts:** no option in any omic changes verdict between the earlier figure and the shipped code (transcriptomics complete over C0 to C8 since the amendment run). Classic technical replicates still pass.
- **Figures that move (shipped code against the figure the app quotes):** lipid/metab recommended worst core FDP 0.072 -> 0.085 (C5) and power 0.43 -> 0.47, extra-hard S3 0.130 -> 0.148, S1 0.101 -> 0.111, S8 0.139 -> 0.180 (M29's exact standard error and M30's 200 permutations, as the exactness round showed); lipid/metab supplied 0.068 -> 0.079; first-pass 0.174 -> 0.203; first-pass with the separation rule 0.239 -> 0.261; contaminated standards 0.371 -> 0.399; standards with the separation rule 0.118 -> 0.095 with 26 % of rows refused by M32 (the option then runs without RUV in a quarter of data sets); lipid/metab technical replicates C6 0.058 -> 0.063. Classic, DIA and transcriptomics figures move by at most 0.004, except transcriptomics first-pass nulls with the separation rule at three groups: power 0.023 -> 0.049 (the app's note: 2 %).
- **Machine-dependent:** the legacy corrected design in classic proteomics (C6, C8) and lipid/metab (section 7 of the declaration): the same seed gives different calls on different runner machines. DIA's bias and rmse differ at 1e-15 between machines (not rule-bearing).
