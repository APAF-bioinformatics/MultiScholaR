# KIT-124 confirmation round: declaration, two parts (decision M34; rewritten 2026-11-19 before any job of the round)

Authority: decisions M29 ("every option's figure is re-measured on the shipped code in a final confirmation round") and M34 (user, 2026-11-19: "(1) for every option, the shipped package code is shown to give the same numbers as the simulator arm it replaces on the same seeds, on a modest number of runs per cell, under the fixed equivalence rule read on the rule-bearing values (M19); where that holds, the earlier measurements stand for the shipped code; (2) options that are not identical are re-measured in full on the shipped code: lipid/metab, the technical-replicate designs once their kernel is repaired, and classic proteomics 'No RUV' under the new test (M30)"), with M28, M30 to M33 (`design-review-20261118.md` section 6 on `performance-v2`). Simulations run on GitHub-hosted runners (M13); no compute on the laptop. Nothing under `R/`, `tests/`, `DESCRIPTION` or `tools/` changes on this branch; the package is installed in each job from the bound commit.

Every choice below that is not carried over from an earlier declaration (`exact-round-declaration.md`, the confirming study, `final-20261118`, the first version of this file) is marked **[for the user to approve]**. The first version of this file (cf572809, one round of every option at full size) is replaced by this one; its lifecycle run 44341242746 stays on record (section 8).

## 1. Package under test

**BINDING LINE (the one line the final commit changes): `package_commit = tools_ref = 2bbb1f6d29e4f4a1767ebb461c08858883370680` (`kit-129-integration`, 2026-11-19, fifth pass item 5).** [for the user to approve: this commit, or the `kit-129-integration` tip at launch; another commit is merged here with the `ours` strategy (the tree of this branch does not change), the hash replaced in this line and in the three launch files, and section 7's lifecycle repeated on it]

It carries M28 to M33, the technical-replicate design and first-pass nulls in every omic, the run-path wiring, classic "No RUV" and a count of 0 as the declared state without factors with M30's test, supplied lists on the classic and lipid/metab run paths, and the separation rule on the run path where it is defined. Made reachable on this branch by 1880f7df (merge, `ours`). Installed as before: params `package: "installed"`, the library action installs MultiScholaR from that checkout with the two internal packages built in the job from `DEPS_READ_TOKEN` (never printed, cached or uploaded); the job checks the checkout is the bound commit and that the namespace holds every entry point the arms call.

## 2. The shipped arms (34 to 41) and how each calls the package

Each calls the installed package as the run path does (the compute halves the tab and the managed worker call), not the work tree; each runs on its own random stream; the count's seed is each simulator's E seed. Arms 1 to 33 are unchanged.

| Arm | S1 classic proteomics | S2 DIA with limpa | S3 transcriptomics | S4 lipid/metab |
|---|---|---|---|---|
| 34 recommended (groups, all features, estimated count) | `computeProtRuvModelFactors` with the run path's replicate sets (`protRuvModelReplicateSets`); count 0 = `declareProtRuvNoFactors`; `differentialAbundanceAnalysis` | `diaRuvCovariatesOperation("declared")` kernel, all proteins; in C6/C8 this is also the technical-replicate option (one computation on this route) | `trxChooseRuvDesign` + `trxRuvFactors` + `trxDaFit`/`trxDaApply` | `runOmicsRuvFactors` + `runLipidsDA`, same through the metab classes (must agree) |
| 35 supplied list, estimated count | arm 8's oracle set | as S1 | arm 8's 50 genes | the standards |
| 36 contaminated list, estimated count | arm 9's set (not C4) | as S1 (not C4) | arm 9's set | 2 of 10 standards changing |
| 37 No RUV | `declareProtRuvNoFactors(ruv_design = "none")`: M30's test | DA on the quantified object | design `none` | no factors (k reported NA, as arm 1) |
| 38 technical replicates (C6, C8) | units as replicate sets (`technical_replicates`) | (= 34) | unavailable (summed at the design step) | `ruv_grouping_variable = "extract"` |
| 39 first-pass nulls, estimated count | kernel's own first pass; units on C6/C8 | kernel's own (dpcDE F test) | kernel's own | kernel's own |
| 40 first-pass nulls, separation rule | the sweep's k by `resolveProtRuvParameters` (tab defaults 1-20 %, on arm 3's stream) | kernel's rule (the sweep's k) | `trxSeparationRuleK` | `omicsRuvSeparationRuleK` |
| 41 supplied list, separation rule | `protRuvSeparationRuleK` on the list | kernel's rule (the sweep's k, whatever the list) | `trxSeparationRuleK` on arm 8's genes | `omicsRuvSeparationRuleK` on the standards |

Legacy designs are the simulators' existing arms 3 (S1, S3, S4) and 3, 4 (S2), which call the package's own legacy functions; they run on the installed package. A user's typed k is not simulated (no rule to simulate); its figure is the estimated count's (KIT-129 section 1a). **[for the user to approve]**

Refusals (carried over, R-07): a classed refusal of the package (`ruv_k_refused`, `ruv_factor_count_unavailable`, `prot_ruv_factors_unavailable`, `omics_ruv_no_residual_df`, the block-correlation and no-replication classes; S3 `trx_error`; S2 a run the kernel leaves without RUV other than by a count of 0; S4 an assay the kernel refuses) is `refused` and scored as arm 1. Any other error is `crashed`, counts in rule 6, and fails the job after its upload (`check-arms.R`).

## 3. Options per omic: shipped arm, the simulator arm it replaces, identity expected?

"Identity" = agreement under `packaging/verify/compare-results.R` of the bound commit (relative 1e-6, absolute floor 1e-9, discrete values exact) on the M19 rule-bearing values (status, k, n_features, n_true, n_calls, null_p05, null_p01, fdp, power), per cell, run and subset. The counterpart is the arm the option's figure in `omicGuiRuvTable*` came from, or, where the final method was measured, its exact-round variant.

| Omic | Option (app id) | Shipped ~ simulator | Identity expected? | Why (not) | Part |
|---|---|---|---|---|---|
| Classic | `groups_all_features` | 34 ~ 32 | yes | arm 32 is the final method; run path proven (KIT-129 `runpath`) | 1 |
| Classic | `technical_replicates` | 38 ~ 33 | yes | arm 33 is the final method on the units; proven (`techrep`, `runpath`) | 1, and 2 (M34 names it) |
| Classic | `supplied_controls` | 35 ~ 8 | no | arm 8: E with 20 permutations, `contrasts.fit` standard error, treat without trend/robust | 2 |
| Classic | supplied, contaminated | 36 ~ 9 | no | as 35 | 2 |
| Classic | supplied, separation rule | 41 ~ 22 | k yes, values no | same rule on the same list; arm 22 has the earlier test | 2 |
| Classic | `first_pass_nulls` | 39 ~ 23 | no | arm 23: 20 permutations, earlier test, groups not units on C6/C8 | 2 |
| Classic | first-pass, separation rule | 40 ~ 5 | k yes, values no | same sweep, same stream; arm 5 has the earlier test | 2 |
| Classic | `legacy` | 3 (installed) ~ 3 (`exact-gha` rows) | yes | legacy code unchanged (M3); compared with the earlier runner rows (section 5) | 1 |
| Classic | `none` | 37 ~ 27 | yes | arm 27 is arm 1 under M30's test | 1, and 2 (M34 names it) |
| DIA | `groups_all_features` / `technical_replicates` | 34 ~ 32 | yes | proven (`runpath`, `m30b`) | 1 |
| DIA | `supplied_controls` | 35 ~ 8 | no | 20 permutations, `contrasts.fit` standard error | 2 |
| DIA | supplied, contaminated | 36 ~ 9 | no | as 35 | 2 |
| DIA | supplied, separation rule | 41 ~ 22 | no | **the kernel's rule is today's sweep's k, not arm 22's rule on the list** (see 9.4) | 2 |
| DIA | `first_pass_nulls` | 39 ~ 23 | no | 20 permutations, standard error; with 3 groups (C7) the first pass is an F test, arm 23's a B-A test | 2 |
| DIA | first-pass, separation rule | 40 ~ 5 | k yes, values no | standard error | 2 |
| DIA | `legacy` (factors in test), `legacy_correction` | 4, 3 (installed) ~ earlier rows | yes | code unchanged; the earlier rows are laptop rows only (section 5) | 2 |
| DIA | `none` | 37 ~ 1 | yes | DIA's no RUV unchanged | 1 |
| Trx | `groups_all_features` | 34 ~ 31 | yes | arm 31 is the final method (edgeR fits each comparison exactly) | 1 |
| Trx | `supplied_controls` | 35 ~ 8 | no (by construction) | arm 8's k is E with 20 permutations; identical only in runs where the two counts agree | 2 |
| Trx | supplied, contaminated | 36 ~ 9 | no (by construction) | as 35 | 2 |
| Trx | supplied, separation rule | 41 ~ 22 | yes | run path = arm 22 (`separation` proof) | 1; C6, C8 in 2 (no earlier rows on the summed tree) |
| Trx | `first_pass_nulls` | 39 ~ 23 | no (by construction) | 20 permutations | 2 |
| Trx | first-pass, separation rule | 40 ~ 24 | yes | run path = arm 24 | 1; C6, C8 in 2 |
| Trx | `legacy` | 3 (installed) ~ earlier rows | yes | code unchanged; earlier rows laptop only | 2 |
| Trx | `none` | 37 ~ 1 | yes | unchanged | 1 |
| Lipid/metab | `groups_all_features` | 34 ~ 32 | subset `complete` yes, `full` no | see below | 2 |
| Lipid/metab | `technical_replicates` | 38 ~ 33 | `complete` yes, `full` no | as 34 | 2 |
| Lipid/metab | `supplied_controls`, contaminated | 35 ~ 8, 36 ~ 9 | no | 20 permutations, `contrasts.fit` standard error | 2 |
| Lipid/metab | supplied, separation rule | 41 ~ 22 | k yes on `complete`; values no | earlier test | 2 |
| Lipid/metab | `first_pass_nulls` | 39 ~ 23 | no | 20 permutations, earlier test | 2 |
| Lipid/metab | first-pass, separation rule | 40 ~ 5 | k yes on `complete`; values no | earlier test | 2 |
| Lipid/metab | `legacy` | 3 (installed) ~ earlier rows | no | the level repair (U8, KIT-123 R2) changes legacy's corrected matrix; the earlier rows predate it | 2 |
| Lipid/metab | `none` | 37 ~ 1 | yes | unchanged | 1 (and in 2 as arm 1) |

**What differs in lipid/metab, exactly.** S4 forms its own W on the features observed in every study **and pooled-QC** injection. The package's kernel forms W on the features observed in every sample **of the object it holds**, and the object holds the study injections only (as S4 builds it; the shipped kernels take the samples of the object's design, and the pooled-QC injections are not among the groups tested). On the subset `complete` (the data restricted to the features complete in study and QC injections) both rules select the same features, so the shipped arm can equal the simulator's (lifecycle: 34~32, 38~33, 40~5, 41~22 identical there). On the subset `full` they select different features: the difference is a property of **the comparison** (the simulator uses information the shipped object does not hold), not a defect of the shipped method; the shipped arm's own figure on `full` (the primary subset) is the one that ships. Giving the package the QC injections would not be "the same information": in the object they would be a further group in the replicate sets and the test. So lipid/metab is re-measured in Part 2 (as M34 says), and Part 1's identity for it is read on `complete` only.

## 4. Part 1, identity

- **Cells and runs [for the user to approve].** Runs 1 to 50 of every core cell C0 to C8 (S3 C6 and C8 on the summed tree) and of the extra-hard cells the option notes quote: S3 for transcriptomics (recommended 5.2 % and legacy 72 % "at a correlation of 0.6"), S1, S3 and S8 for lipid/metab ("10 % at 3 per group", "13 % at 0.6", "14 % with an unbalanced batch"); classic and DIA notes quote none. Why 50: the identity claim is per run and exact, so one disagreeing run decides; 50 runs per cell, 450 to 600 per simulator and 2,000 data sets over the round, cover the spread of k and the refusal paths of each cell far better than the proofs' one run per cell, at about a twentieth of the full cost. A pair that holds in all of them is taken as identical by construction, which the proofs on runs 1 of C0 and C6 already argue.
- **Arms**: arms 1 and 2, every shipped arm, every counterpart, and the legacy arms (`.github/k124-sim/params-confirm-part1.json`).
- **Reading**: the workflow's compare job (`compare-arms.R`, then the bound commit's `compare-results.R`), each subset apart. A pair whose identity is expected and holds in every run of every cell: the earlier rows of its simulator arm stand for the shipped code (`exact-gha` 44233215700, `exact-k128-gha` 44233219582, `newarms-gha` 44200924581). A pair "not expected" is reported (its differences are the stated cause or not) and is in Part 2 regardless of what the 50 runs show (agreement there is per-run coincidence, e.g. two counts that happen to agree).
- **Rule written now for an unexpected failure**: any pair declared "identity expected" that fails on any rule-bearing value in any run of any cell, or whose shipped arm crashes, is added to Part 2 with every cell of its simulator at 1,000 runs (and the quoted extra-hard cells at 500), with no other change; the failing rows and the cause found are reported before Part 2 launches. **[for the user to approve]**
- **Legacy against earlier rows.** S1 arm 3 (installed) against the `exact-gha` runner rows of the same seeds under the same rule, as the arms 1 and 2 gate (section 6). S2, S3 and S4 legacy rows exist on the laptop only, so they go to Part 2 (cheap: their sweep or grid is shared or under 4 s per run). **[for the user to approve]**

## 5. Part 2, full re-measurement on the shipped code

Fixed list (launch files `params-confirm-part2-core.json`, `params-confirm-part2-extra.json`); Part 1 can only add to it (rule above):

| Simulator | Arms (plus 1 and 2) | Cells |
|---|---|---|
| S1 classic | 35, 36, 37 (M34), 38 (M34; C6, C8), 39, 40, 41 | C0 to C8, 1,000 runs |
| S2 DIA | 35, 36, 39, 40, 41; legacy 3, 4 | C0 to C8, 1,000 runs |
| S3 trx | 35, 36, 39; 40, 41 (needed in C6, C8; run in every cell at about 4 s per run); legacy 3 | C0 to C8, 1,000 runs |
| S3 trx | 34, legacy 3 | S3, 500 runs |
| S4 lipid/metab | 34 to 41; legacy 3 | C0 to C8, 1,000 runs |
| S4 lipid/metab | 34 | S1, S3, S8, 500 runs |

Pass rules 1 to 7 **unchanged** through `aggregate.R` as on `sim-k124-exact` (sha256 `089d852c...`): rule 1 null rate at 0.05 and 0.01 (upper bounds 0.075, 0.015), rule 2 mean false discovery proportion (0.075), rule 3 share above 0.10 against arm 2's share + 0.05, rule 4 the reduction against arm 1 in cells with a factor, rule 5 power loss against arm 1 without one (0.02), rule 6 crash share (0.01), rule 7 global null C4 (0.075); void cells, the primary subset (`full`), C8 reported, extra-hard cells listing failing rules without failing a design: as before. `aggregate.R`'s recommendation table does not know arms 34 to 41 and is not read; the per-arm and per-cell verdicts are. The figure per option is read as before (mean FDP in the worst core cell with true effects; power in C0).

## 6. Gate (carried over, M19)

Before a cell of either part is scored, its arms 1 and 2 equal the earlier hosted-runner rows of the same seeds under the same rule (references as in section 4; for the extra-hard cells the runner rows of `final-20261118` 44201244793 where they exist). A cell that fails the gate is not scored and is reported. **Not yet built for this round:** the gate and the S1 legacy comparison need a runner step that downloads the reference run's artifacts by run id and runs `compare-laptop.R` and `compare-results.R` (the exact round did this on the laptop, now excluded). The `exact-gha` artifacts expire about 2026-12-18. **[for the user to approve: a `gate` job added to the workflow before Part 1]**

## 7. Lifecycle run (the instrument check; 2026-11-19)

Run **44419300892** (branch `sim-k124-confirm`, commit c1cae181, label `confirm-lc2`), package installed from 2bbb1f6d in every job (meta lines confirm), S1 C6, S2 C0, S3 C0, S4 C6, runs 1 and 2, 2 cores, every shipped arm beside its counterpart. All six jobs succeeded; every shipped arm `ok` in every run (no refusal, no crash: arm 38 now runs in S1 and S4); the compare job compared every declared pair. The failure path (a crashed shipped arm is uploaded, then fails the job) was shown on this instrument by 44341242746 and not repeated. Tables (not committed): sha256 `50c17531...` (s1), `b6fe491b...` (s2), `6bf4473c...` (s3), `5153fe45...` (s4); artifacts kept 30 days.

Per pair (rule-bearing values; largest relative difference where not identical):

| Pair | S1 C6 (full / complete) | S2 C0 | S3 C0 | S4 C6 (full / complete) |
|---|---|---|---|---|
| 34 ~ 32 / 31 | identical / identical | identical | identical | differ (fdp 0.12, power 0.13, null_p05 0.25) / identical |
| 35 ~ 8 | null_p05 0.057 / identical | fdp 0.13, null_p01 0.33 (k equal) | identical | differ (k: counts differ) / differ |
| 36 ~ 9 | fdp 0.016, n_calls / fdp 0.016 | fdp 0.14 | identical | differ / differ |
| 37 ~ 27 / 1 | identical / identical | identical | identical | identical / identical |
| 38 ~ 33 | identical / identical | (not formed) | (not formed) | differ (fdp 0.43) / identical |
| 39 ~ 23 | null_p05 0.029 / identical | fdp 0.005 (k equal) | identical | differ / differ (k) |
| 40 ~ 5 / 24 | k equal; fdp 0.008, power 0.008 / null_p05 0.18 | k equal; fdp 0.005 | identical | k equal; fdp 1, power 0.32 / identical |
| 41 ~ 22 | k equal; fdp 0.15 / identical | **k differs**; fdp 0.18 | identical | k equal; power 0.05 / identical |

Every expected identity held; every difference is of the declared kind (S3's 35, 36, 39 agreed because the two counts happened to agree in both runs). Seconds per run inside a worker, one core (all arms of the table):

| Simulator, cell | Shipped arms | Counterparts and others | Whole run (job's simulator seconds x cores / runs) |
|---|---|---|---|
| S1 C6 | 34-36, 38, 39: 2.1-2.2; 37: 1.0; 40: 13.2 (its sweep); 41: 1.6 | 3: 13.2; 32: 1.1; others under 1 | about 48 s (other cells: the sweep is 23-27 s, about 70 s) |
| S2 C0 | 34-36, 39: 2.9-3.7; 37: 0.4; 40, 41: 15.8 each (each runs the sweep) | 32: 1.3; others under 0.2 (shared quantification and sweep about 18 s) | about 65 s (C6: 34-36 9 s, 40 36 s, about 150 s, from 44341242746) |
| S3 C0 | 34-36, 39: 16.3-17.9 (each its own count, 200 permutations); 37: 0.6; 40: 2.2; 41: 1.6 | 31: 8.1; 3: 2.4; others 0.6-1.5 | about 100 s |
| S4 C6 | 34-36, 38-40: 2.2-2.7 (both classes); 37: 0.7; 41: 1.3 | 3: 3.6; others under 1 | about 60 s (both subsets) |

## 8. Budget (from the lifecycle; core-hours are simulator compute; one cell measured per simulator, so a cell may differ by about half either way)

| | Runs | Core-hours | Jobs (4 cores, limit 300 min) | Wall time at 20 concurrent jobs |
|---|---|---|---|---|
| Part 1 | S1 450, S2 450, S3 500, S4 600 | S1 6-9, S2 10-11, S3 14, S4 10: **about 40 to 45** | 40 (one per simulator and cell) | **about 1.5 h** |
| Part 2 core | 9,000 per simulator | S1 70-105, S2 about 195 (C6, C8 at about 140 s per run), S3 about 175, S4 about 100: **about 540 to 575** | 63 (S2 in ranges of 400 runs, S3 of 500) | **about 16 to 18 h** |
| Part 2 extra-hard | S3 500, S4 1,500 | **about 13** | 4 | about 1 h |

Each job also spends about 6 to 8 minutes building the library (not counted above: about 5 job-hours in Part 1, 8 in Part 2). One attempt per job; one infrastructure rerun (`gh run rerun <run> --failed`, the same rows). A job that reaches its limit or a crash of a shipped arm is reported, not relaunched with other settings.

**Part 2 is about 4.5 times the carried-over limit of about 120 core-hours [the user's decision].** Levers that keep every number identical (each checked in Part 1 by k identity), none implemented:

| Lever | Saving, about | Note |
|---|---|---|
| L1 One estimated count per data set shared by the shipped arms (the package's own `estimateRuvFactorCount` with the same seed; each kernel otherwise repeats it) | 85 (S3), 10 (others) | k identity per run already shown (prove `count`, lifecycle) |
| L2 One DIA sweep per data set for arms 40 and 41 (the kernel's sweep call does not depend on the control rule) | 60 (S2) | |
| L3 Part 2 reuses Part 1's runs 1 to 50 (same seeds, commit, code) | 5 % | needs a first-run entry in the workflow |
| L4 Drop S1 37 and 38 from Part 2 if Part 1 shows identity (arms 27, 33 have 1,000-run rows) | 5 | M34 names them for Part 2; they are cheap |

With L1 and L2 Part 2 is about 400 to 430 core-hours.

## 9. Points for the user

1. Bind 2bbb1f6d, or the tip at launch (section 1).
2. Part 1's cells and 50 runs; the extra-hard cells quoted by the notes only (section 4).
3. The rule adding an unexpected failure to Part 2 (section 4).
4. **DIA, supplied list with the separation rule:** the app's figure (0.056, `omicGuiRuvTableDiaLimpa`) cites arm 22 (the separation rule on the list); the shipped kernel uses today's sweep k whatever the list (lifecycle: k differs in both runs). Part 2's arm 41 measures what ships; the figure source must change, or the kernel follow arm 22.
5. Legacy: S1 by identity against `exact-gha`; S2, S3, S4 re-measured (section 4); lipid/metab legacy is not expected to be identical (level repair).
6. Lipid/metab identity read on subset `complete`; `full` re-measured (section 3).
7. A user's typed k not simulated (section 2).
8. The gate job (section 6) must be added before Part 1.
9. Part 2's size and the levers L1 to L4 (section 8).

## 10. Launch, once approved

1. If the binding changes: merge the commit with `git merge -s ours <hash>`, replace the hash in section 1 and the three launch files, and repeat the lifecycle (`params.json` as for 44419300892, new label).
2. Add the gate job (item 8) and its own lifecycle check.
3. Part 1: `cp .github/k124-sim/params-confirm-part1.json .github/k124-sim/params.json`, commit, push `sim-k124-confirm` (the push of `params.json` starts the run). Read the compare artifact `k124-confirm-p1-compare` and the gate; apply section 4's rule; record the outcome here.
4. Part 2: the same with `params-confirm-part2-core.json`, then `params-confirm-part2-extra.json` (each push one run), plus any pairs Part 1 added. Assemble (`assemble.sh <run> <label>`) and score with `aggregate.R` unchanged.

`params.json` holds the lifecycle's parameters; a push that leaves it alone starts nothing.
