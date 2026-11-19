# KIT-124 confirmation round: declaration, two parts (decision M34; rewritten 2026-11-19 before any job of the round)

Authority: decisions M29 ("every option's figure is re-measured on the shipped code in a final confirmation round") and M34 (user, 2026-11-19: "(1) for every option, the shipped package code is shown to give the same numbers as the simulator arm it replaces on the same seeds, on a modest number of runs per cell, under the fixed equivalence rule read on the rule-bearing values (M19); where that holds, the earlier measurements stand for the shipped code; (2) options that are not identical are re-measured in full on the shipped code: lipid/metab, the technical-replicate designs once their kernel is repaired, and classic proteomics 'No RUV' under the new test (M30)"), with M28, M30 to M33 (`design-review-20261118.md` section 6 on `performance-v2`). Simulations run on GitHub-hosted runners (M13); no compute on the laptop. Nothing under `R/`, `tests/`, `DESCRIPTION` or `tools/` changes on this branch; the package is installed in each job from the bound commit.

**Approved (decision M38, user, 2026-11-19, recorded on `performance-v2`):** the round as declared here, Part 1 then Part 2, **with savings L1 and L2** (section 8); a pair expected to be identical that fails joins Part 2 (section 4); lipid/metab identity read on the subset `complete`; a typed k not simulated; the carried-over limit of about 120 core-hours lifted for this round. Not yet decided: the binding commit (section 1). **The binding and the launch wait on the main session: nothing of this round is launched from this branch until it sends the commit.** The marks "[for the user to approve]" below are kept as the record of what M38 approved.

Every choice below that is not carried over from an earlier declaration (`exact-round-declaration.md`, the confirming study, `final-20261118`, the first version of this file) is marked **[for the user to approve]**. The first version of this file (cf572809, one round of every option at full size) is replaced by this one; its lifecycle run 44341242746 stays on record (section 8).

## 1. Package under test

**BINDING LINE (open; the one line the final commit changes): `package_commit = tools_ref = 44b5f2195e355bb6e6efe1dae21905e30b81c26c`** (the `kit-129-integration` tip after the KIT-129 sixth pass, sent by the main session; M38 rejected 2bbb1f6d. The sixth pass makes DIA's supplied list with the separation rule follow the rule on the list, so DIA's pair 41~22 is then expected identical in k.) Binding: merge it here with the `ours` strategy (the tree of this branch does not change), write its hash in this line and in the three launch files, and repeat section 7's lifecycle on it (section 10).

The lifecycles below ran on 2bbb1f6d (`kit-129-integration`, fifth pass item 5), made reachable by 1880f7df (merge, `ours`); it carries

M28 to M33, the technical-replicate design and first-pass nulls in every omic, the run-path wiring, classic "No RUV" and a count of 0 as the declared state without factors with M30's test, supplied lists on the classic and lipid/metab run paths, and the separation rule on the run path where it is defined. Made reachable on this branch by 1880f7df (merge, `ours`). Installed as before: params `package: "installed"`, the library action installs MultiScholaR from that checkout with the two internal packages built in the job from `DEPS_READ_TOKEN` (never printed, cached or uploaded); the job checks the checkout is the bound commit and that the namespace holds every entry point the arms call.

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
| DIA | supplied, separation rule | 41 ~ 22 | k yes, values no | from the bound commit's code (sixth pass item 3): the kernel runs `diaRuvSeparationRuleK` on the list, S2's `separationK` (proven `diaseparation`, C0 run 1), so k equals arm 22's; arm 22 tests with `contrasts.fit`'s standard error, the shipped DA with M29's exact one, so the values differ (as 35~8). Written before `confirm-lc4` was read. At 2bbb1f6d the kernel used the sweep's k (lifecycle 44419300892: k differed) | 2 |
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

## 6. Gate (carried over, M19; built 2026-11-19, decision M38)

Before a cell of either part is scored, its arms 1 and 2 equal the earlier rows of the same seeds under the same rule (`compare-results.R` of the bound commit, relative 1e-6, absolute floor 1e-9, discrete values exact, read on the M19 values). A cell that fails the gate is not scored and is reported.

- **The job.** `gate` in `.github/workflows/k124-sim.yml` runs after the simulators of every round with the package installed: per simulator `.github/k124-sim/gate.R` takes the run's arms 1 and 2 (and in S1 arm 3, the legacy comparison of section 4) and the committed reference rows of the same keys through `compare-laptop.R` and `compare-results.R`, and writes `<sim>-gate-summary.tsv` and `s1-legacy-summary.tsv` (artifact `k124-<label>-gate`). A rule-bearing value outside the rule, or a run key the reference does not hold, fails the job after its tables are uploaded; a cell with no reference row at all is reported "no reference" (a warning).
- **The committed reference** (`scored/confirm-reference/`, so that the round does not depend on artifacts that expire about 2026-12-18): the rows of arms 1 and 2 (and S1 arm 3), every column as written by the simulators, extracted by a fixed text filter (`reference-extract.sh`, `reference-build.sh`, `reference-spec.tsv`: rows selected and sorted, nothing computed):

| File | Rows | Source |
|---|---|---|
| `s1-runs.tsv.gz` | 54,000 (arms 1, 2, 3; C0 to C8; runs 1 to 1,000; full and complete) | `exact-gha` 44233215700, 18 artifacts |
| `s2-runs.tsv.gz` | 18,000 (arms 1, 2; C0 to C8) | `exact-gha`, 9 artifacts |
| `s3-runs.tsv.gz` | 18,000 (C0 to C5, C7; C6, C8 on the summed tree) | `exact-gha`, 7 artifacts; `exact-k128-gha` 44233219582, 2 |
| `s4-runs.tsv.gz` | 36,000 (C0 to C5, C7; C6, C8 on the kit-128 tree; both subsets) | `exact-gha`, 7; `exact-k128-gha`, 2 |
| `s4-runs.laptop-extra.tsv.gz` | 6,000 (arms 1, 2; S1, S3, S8; runs 1 to 500; both subsets) | the laptop confirm run, tree 0b6fcb6d (2026-11-18 04:04 to 06:16): the only rows of those cells |

  `sources.tsv` gives each source table's sha256 and the rows taken, `artifacts.tsv` the artifacts' digests and expiry, `tables.sha256` the decompressed files. **Verified once against the downloaded artifacts** (lifecycle 44423530060, params `verify_reference: true`): every listed artifact present with the recorded table sha256, the rows re-extracted on the runner byte-identical to the committed copy, the laptop rows equal to their recorded sha256 (they cannot be re-downloaded). A missing artifact stops the step with its name.
- **Extra-hard cells without reference:** S3's cell S3 (its earlier rows are on the laptop in a binary file only, which cannot be read without R there): reported "no reference", its gate not applied **[for the user]**. S4's S1, S3, S8 are gated against laptop rows (reference BLAS against the runners' OpenBLAS, as M19 accepted for the exact round).

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


**Second lifecycle (savings and gate), run 44423530060** (commit bd57e073, label `confirm-lc3`): the cells, runs, arms and package of `confirm-lc2` with `share: "count,sweep"` and `verify_reference: true`. All seven jobs succeeded. The reference verification passed (above). The gate: S1, S2, S3, S4 arms 1 and 2 and S1 arm 3 agree with the committed reference in every rule-bearing value (status 0).

Identity of the savings, every row of `confirm-lc3` against `confirm-lc2` by cell, run, arm and subset (the same seeds and package; no sharing in lc2):

| Simulator | Shipped arms 34 to 41 | Every other arm |
|---|---|---|
| S1 C6 | byte-identical (76 rows of the table) | byte-identical |
| S3 C0 | byte-identical | byte-identical |
| S4 C6 | byte-identical | byte-identical but arm 3 (legacy) run 2, `full`: null_p05 0.300 against 0.295, null_p01, bias, rmse; arm 3 calls neither shared function, so this is the legacy lipid sweep differing between two runner machines |
| S2 C0 | every rule-bearing value identical; bias or rmse differ in the 15th significant digit (relative about 1e-15) | the same 1e-15 differences in arms 1, 2, 4, 5, 8, 9, 22, 23, 32, which share nothing: the DPC quantification is not bit-reproducible across runner machines |

So L1 and L2 change no value: where the table is reproducible across machines it is byte-identical, and in S2 the only differences are the cross-machine rounding every arm shows. Two findings for the reading: S2's bias and rmse are not bit-reproducible across machines (rule-bearing values are; M19 reads those), and S4's legacy arm 3 can differ between machines in a rule-bearing value (null_p05) for the same seed, so lipid/metab legacy figures carry that machine dependence (it is re-measured in Part 2, not read for identity).

Seconds per run with the savings (one core; the first arm that needs the count or the sweep pays for it, the others take it):

| Simulator, cell | Shipped arms | Whole run |
|---|---|---|
| S1 C6 | 34: 2.2 (pays the count); 35, 36, 39: 1.2-1.3; 37: 1.1; 38: 2.1 (units: its own count); 40: 13.4 (sweep); 41: 1.4 | about 51 s |
| S2 C0 | 34: 5.2 (count); 35, 36, 37, 39: 0.8-1.2; 40: 31.0 (the sweep, on a slower machine than lc2's 15.8); 41: 1.0 | about 84 s on this machine (lc2's machine about 65 s without savings) |
| S3 C0 | 34: 19.5 (count); 35, 36: 1.3; 39: 2.0; 40: 2.8; 41: 1.5 | about 72 s (was about 100) |
| S4 C6 | 34-36, 38: 1.7-1.8; 39: 1.5; 40: 2.3; 41: 1.3 | about 60 s |

## 7b. Binding lifecycle and Part 1 (2026-11-19, on 44b5f2195e355bb6e6efe1dae21905e30b81c26c)

**Binding lifecycle `confirm-lc4`, run 44434645817** (lc3's cells, runs, arms and savings, with the reference verified): every job succeeded; every pair as section 3 declares, DIA 41~22 now equal in k (values differ, as written before the run); the gate passed (status 0).

**Part 1 `confirm-p1`, run 44436553563** (runs 1 to 50; 40 simulator jobs, all succeeded; 49 minutes wall; about 42 core-hours of simulator compute: S1 8.0, S2 15.3, S3 12.9, S4 6.1). Every shipped arm `ok` in every row, but S4's 40 and 41 (and some 34 to 39 assays) `refused` in 473 of 8,600 rows, by the kernel's limits (below). Tables not committed; sha256 in `confirm-p1-tables.sha256`.

Identity, per pair and subset, read once (rule-bearing values; "identical" = every run of the cell agrees):

| Pair | Subset | Cells | Identical cells | Largest relative difference | Declared | As declared? |
|---|---|---|---|---|---|---|
| S1 34~32 | full, complete | 9 | 9 | 0 | yes | yes |
| S1 37~27 | full, complete | 9 | 9 | 0 | yes | yes |
| S1 38~33 | full, complete | 2 (C6, C8) | 2 | 0 | yes | yes |
| S1 35~8, 36~9, 39~23 | both | 9 (36: 8) | 0 | up to 1 | no | yes |
| S1 40~5, 41~22 | both | 9 | 0 (k equal in all 9) | up to 1 | k yes, values no | yes |
| S2 34~32, 37~1 | full | 9 | 9 | 0 | yes | yes |
| S2 35~8, 36~9, 39~23 | full | 9 (36: 8) | 0 | up to 1 | no | yes |
| S2 40~5, 41~22 | full | 9 | 0 (k equal in all 9) | 0.5, 1 | k yes, values no | yes |
| S3 34~31, 37~1, 40~24, 41~22 | full | 10 (C0 to C8, S3) | 10 | 0 | yes | yes |
| S3 35~8, 36~9, 39~23 | full | 10 | 8 (not C3, C8) | 0.81, 1, 0.14 | no (by construction) | yes |
| S4 34~32, 38~33, 37~1 | complete | 12 (38: 2) | all | 0 | yes | yes |
| S4 37~1 | full | 12 | 12 | 0 | yes | yes |
| S4 34~32, 38~33 | full | 12, 2 | 1, 1 | 1 | no | yes |
| S4 35~8, 36~9, 39~23 | both | 12 | 6 on complete, 0 on full | 1 | no | yes |
| S4 40~5, 41~22 | complete | 12 | 11, 3 | 0.78, 1 | k yes on complete, values no | **no, in k**: see below |

**Not as declared (S4 40~5 and 41~22, k on `complete`).** The separation rule gives the same k on both sides; the shipped kernel then refuses an assay where that k is not below the controls observed in every sample (decision M32; 5 standards per assay for 41) or leaves no residual degree of freedom, and reports k = 0 for it (scored as arm 1), where arms 5 and 22 accepted k equal to the controls. It follows from M32, not from a call error; both pairs were in Part 2 already, so nothing is added.

**Gate.** Arms 1 and 2 agree with the committed reference in every gated cell of S1 to S4 (C0 to C8; S4 S1, S3, S8). S3's cell S3: no reference (declared). **The S1 legacy comparison failed in C6:** 21 of the 100 rows of arm 3 (runs 1 to 50, both subsets) differ from `exact-gha` (same k; calls, null rates, fdp and power differ; e.g. run 1); C0 to C5, C7, C8 agree. In `confirm-lc4` runs 1 and 2 of C6 had agreed, and in this run run 1 differs: the same code and seed give different legacy results on different runner machines in C6, as S4's legacy arm 3 did in section 7 (RUV-III-C's correction on the technical-duplicate design is not reproducible across machines). The gate job therefore ended `failure` (by design, after the upload).

**Section 4's rule applied:** the expected-identical pair S1 legacy (arm 3) failed, so arm 3 joins Part 2 for S1 at 1,000 runs in every core cell (`params-confirm-part2-core.json`, S1 arms now 1, 2, 3, 35 to 41). No other expected-identical pair failed; no shipped arm crashed. Note for the reading: S1 legacy C6 (and S4 legacy) figures are machine-dependent; Part 2 measures them once, on whichever machines run them.

## 8. Budget (with L1 and L2; from Part 1's measured seconds)

| | Runs | Core-hours, about | Jobs (4 cores, limit 300 min) | Wall time at 20 concurrent jobs |
|---|---|---|---|---|
| Part 1 (done) | 2,000 | 42 measured | 40 | 49 min |
| Part 2 core (`params-confirm-part2-core.json`, final list: S1 1, 2, 3, 35 to 41; S2 1, 2, 3, 4, 35, 36, 39, 40, 41; S3 1, 2, 3, 35, 36, 39, 40, 41; S4 1, 2, 3, 34 to 41) | 9,000 per simulator | S1 about 150 (about 60 s per run: the legacy sweep and arm 40's sweep), S2 about 250 (about 100 s: quantification, legacy sweep, arm 40's sweep; C8 about 145 s), S3 about 110 (about 44 s; arm 35 pays the count), S4 about 90 (about 37 s): **about 600** | 72 (S1 in ranges of 500 runs, S2 of 300) | **about 12 to 14 h** |
| Part 2 extra-hard (`params-confirm-part2-extra.json`) | S3 500, S4 1,500 | **about 14** | 4 | about 1 h |

The core estimate counts the whole job time per run on 4 cores (load imbalance at a job's end included), so it is an upper reading; by arm seconds alone it is about 470. Each job also builds the library (6 to 8 minutes) and every run has a gate job. One attempt per job; one infrastructure rerun. The limit of about 120 core-hours is lifted (M38).

## 9. Decisions (M38) and what remains open

Approved by M38: sections 2 to 8 as written, with L1 and L2; the rule of section 4; lipid/metab identity on `complete`; a typed k not simulated; the limit lifted. Open:

1. **The binding commit** (section 1): the main session sends it. Until then nothing is launched.
2. **DIA, supplied list with the separation rule:** at 2bbb1f6d the kernel used today's sweep k whatever the list (k differed from arm 22's in both lifecycle runs); the sixth pass makes it follow the rule on the list. On the bound commit the pair 41~22 is then expected identical in k (values still differ by the standard error, M29), and section 3's row reads "k yes, values no"; the binding lifecycle shows it.
3. S3's extra-hard cell S3 has no gate reference (section 6).
4. The cross-machine findings of section 7 (S2 bias and rmse at 1e-15; S4 legacy arm 3 at a rule-bearing value).

## 10. Launch sequence (after the main session sends the binding commit)

```
cd /home/doktersmol/Projects/APAF/MultiScholaR-sim-k124-confirm
git fetch origin kit-129-integration sim-k124-confirm && git status        # clean, on sim-k124-confirm
H=<the binding commit, 40 hex digits>
git merge -s ours --no-edit -m "KIT-124 confirmation round: kit-129-integration $H made reachable (ours; tree unchanged)" "$H"
git diff HEAD^1 HEAD --stat                                                 # empty: the tree did not change
sed -i "s/BINDING-COMMIT-SENT-BY-THE-MAIN-SESSION/$H/g" .github/k124-sim/params-confirm-*.json   # the launch files carry this placeholder, which the workflow refuses
sed -i "s/package_commit = tools_ref = _\{40\}/package_commit = tools_ref = $H/" docs/qa/workflow-execution/kit/kit-124/scored/confirm-round-declaration.md
# binding lifecycle: lc3's cells, runs, arms, savings and verification on the bound commit
T=$(mktemp); jq --arg h "$H" '.label="confirm-lc4" | .package_commit=$h | .tools_ref=$h' .github/k124-sim/params.json > "$T" && cp "$T" .github/k124-sim/params.json
git commit -am "KIT-124 confirmation round: bound to kit-129-integration $H; the binding lifecycle (confirm-lc4)" && git push origin sim-k124-confirm
#   read: every job success; k124-confirm-lc4-compare (expected identities as section 3, DIA 41~22 equal in k); k124-confirm-lc4-gate status 0
# Part 1
cp .github/k124-sim/params-confirm-part1.json .github/k124-sim/params.json
git commit -am "KIT-124 confirmation round: Part 1 (confirm-p1) on $H" && git push origin sim-k124-confirm
#   read k124-confirm-p1-compare and k124-confirm-p1-gate; apply section 4's rule; record the outcome here; add any failing pair to
#   params-confirm-part2-core.json (its arm to the simulator's arms) before Part 2
# Part 2 (two pushes, one run each)
cp .github/k124-sim/params-confirm-part2-core.json .github/k124-sim/params.json && git commit -am "KIT-124 confirmation round: Part 2 core (confirm-p2c) on $H" && git push origin sim-k124-confirm
cp .github/k124-sim/params-confirm-part2-extra.json .github/k124-sim/params.json && git commit -am "KIT-124 confirmation round: Part 2 extra-hard (confirm-p2x) on $H" && git push origin sim-k124-confirm
#   assemble.sh <run> <label> <out root>, then aggregate.R unchanged; the gate artifacts decide which cells are scored
```

`params.json` holds the last lifecycle's parameters (`confirm-lc3`); a push that leaves it alone starts nothing.

## 11. Part 2 and the scoring (2026-11-19; results in `confirm-round-results.md`, tables in `confirm-gha/`)

Part 2 core run 44443281402 (`confirm-p2c`, 72 jobs) and extra-hard run 44443334468 (`confirm-p2x`, 4 jobs) on the bound commit, with L1 and L2, launched as listed in section 8 (decision M38). Two jobs, S3 C1 and S3 C7, reached the 300-minute limit and uploaded nothing; per section 8 they are reported and not relaunched (a relaunch in run ranges is a settings change: the user's decision). Scored once by the scoring workflow (run 44480610832) with `aggregate.R` 089d852c (named in `score.json` before any result was read), one call per simulator. Gate: every gated cell agrees on every rule-bearing value; the S1 legacy comparison differs in C6 and C8 (machine-dependent, as in 7b). Nothing was re-set or adjusted after reading.
