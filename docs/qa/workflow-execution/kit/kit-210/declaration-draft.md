# KIT-210 declaration: DRAFT (not the declaration)

**Status: draft, written 2026-11-24 against v0.9-staging 53713252 with no run of any kind; revised the same day for the
auditor's Re-check 1 (R1-01 to R1-16, `tickets/kit/audits/AUDIT-KIT-210-plan-20261124.md`) and the coordinator's
rulings C-1 to C-5.** It is not `declaration.md` and binds nothing. Every measured figure is a named blank `[L1: ...]`
that the blind lifecycle run L1 fills. The script hashes below are **provisional**: the scripts are untested (no R on
the laptop) and will change at L1 under VAL-54 if L1 shows a mechanics defect, and again if the bound commit changes a
call's signature (section 9). The declaration is committed alone, after L1 and the auditor's check of this draft and
of L1, before any counted run (KIT-210 10, A5). Design: KIT-210 revision 1.4 (`tickets/kit/KIT-210.md`), the owner's
decisions VAL-152 and VAL-156.

## 1. Bound code

- Package: the bound commit of E1, `[E1: the first v0.9-staging commit with KIT-190, KIT-196, KIT-200, KIT-201 and
  KIT-202, and the VAL-142 alignment and KIT-188's ruling if either lands before and changes code on the filter, the
  RUV-III-C path or the legacy DA (the coordinator's C-1)]`. Not available on 2026-11-24.
- Every request's `package_commit`: this ticket's branch `kit-210-ruviiic-remedy` rebased or merged onto the bound
  commit, adding only `docs/qa/workflow-execution/kit/kit-210/` `[E1: commit]`.
- Scripts (SHA-256, provisional; recomputed in the declaration):
  - `sim-k210.R` 4e520db7a41404fd24ba3356ec50d193f82fb08fb682872ab9cfd509c12403fe
  - `ruviiic-ref210.R` 61b86b48c67c8d9b2dd36416242f44a085eefd921e430c2cac7e09c9b4f89add
  - `score-k210.R` 25dcd8bc2fb8a45aedd2a795f412e1c5858ad9096467358707d99ebdd7d9cacf
  - derived from `kit-185/sim-k185.R` 5ffc96049efbcba0c7c465c78201465aa7a6cc13490b20afa6e705524aa234ce, `ruviiic-ref.R`
    019ee0a6a039f3bf1b842dbff6543e964dda08ce393908a933ade1ea40e39130, `score-k185.R`
    6847852cc8f4032fb41b96d0e89fa14b9c5d8dec2d45e0a27f22935d85e62554 (the merged files; R1-15)
- Runner: R 4.6.1, the snapshot library of the bound commit, `extra_packages: ["sva"]`, reference BLAS,
  `OMP_NUM_THREADS=1`, one RUV-III-C OpenMP thread, 4 workers, one stream per replicate (KIT-185's).

## 2. Generator, cells, masks, filter, k, controls (KIT-185's, unchanged; KIT-210 6.1)

As `kit-185/declaration.md` section 2, with three differences: seeds
`set.seed(210000000 + 10000 c + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")`
(c as KIT-185; r 1 to 400 counted, 401 to 403 lifecycle); the filter, normalisation and DA are the bound commit's
(E10 states whether the 13 non-P2 cells' inputs equal KIT-185's); no presence variant (KIT-185's G2 decided it is not
taken). Part A: k planted, every feature a potential control, replicate sets as the route resolves them,
normalisation none (cyclic loess only in the sensitivity variant).

## 3. The arms (KIT-210 5, 6.2; `sim-k210.R` fitCandidate and helpers)

Per replicate and mask (`none` included): S, C3, C4a, C4b, C6, C7, O, N; lipidomics and metabolomics also the
reported C7r5. Guard g = max(20, 3k) on C3, C4a, C4b.

| Arm | Exactly as coded |
|---|---|
| S | `ruvIII_C_Varying` S4 method of the cell's class, then the route's post-processing (`fitShipped`, KIT-185's) |
| C3 | `ruviiicRef(directions = "complete", guard = g)` per assay; an assay with fewer than g complete features: the arm is refused (`candidate_too_few_complete`); post-processing as S's |
| C4a | `ruviiicRef(directions = "subset", guard = g)` per pattern of observed runs; |F(S)| < g leaves the target NA and is counted (`guard_targets`), not a refusal |
| C4b | per assay: W from the features observed in every run, all controls, by the declared route's own function (proteomics `computeRuvCovariates` on the complete subset with S's replicate sets; lipid/metab `computeOmicsRuvAssayFactors` with S's replicate matrix); fewer than g complete features: refused; W not formed (that function's own refusal, or k or fewer complete controls): refused (`candidate_w_refused`); then `limma::removeBatchEffect(y, covariates = W, design = M)`; a feature with fewer than k residual dimensions: NA; post-processing as S's |
| C6 | each feature centred on its observed mean, the S4 method on the centred matrix, the mean added back, then the route's post-processing on the uncentred values (the proteomics filter is value-dependent: `R/func_prot_qc_filtering_helpers.R:229-245`) |
| C7 | `ruviiicRef(directions = "nipals", ncomp)` with `imputeNipalsMatrix`; **ncomp = the route's own count when it chooses k**: proteomics 5 (`R/func_prot_norm_optimization_helpers.R:679`), lipidomics/metabolomics the object's stored `ruvCancor` argument, otherwise 2 (`R/func_metab_norm.R:298`; the simulated objects store none: 2). **Features with no observed value are left out of the NIPALS input; they add nothing to the cross-product, as in the zero-filled form** (the package's NIPALS stops on such a column, `R/func_general_nipals.R:22-24`; R1-01); their number is reported (`nipals_left_out`) |
| C7r5 (reported) | as C7 with ncomp = 5; computed for lipidomics/metabolomics only (on proteomics it is C7) |
| O, N | oracle, none (KIT-185's) |
| C (L1 only, E1) | KIT-185's unguarded complete-feature comparator, run only with `K210_CANDIDATES=off` |

Per data set and assay (rows `corr = "assay"`): realised missing share after the filter, kappa-hat (A-8: Q-hat with the
co-missing term, P over all runs, delta-hat scaled by n / n_c at the part-A k), its perturbation norm and eigengap, the
sine of the largest principal angle between S's directions (zero-filled cross-product) and C3's, the number of complete
features.

## 4. Estimands, feature sets, primary endpoint, rules (KIT-210 6.3 to 6.5; `score-k210.R`)

As KIT-210 r1.4 sections 6.3 to 6.5, implemented as follows:

- F_m(c): the features of F_m^filt finite in c's output on m, in S's on m and in S's on `none`. Rows of c and of S, O, N,
  S0, O0, N0 on F_m(c) carry `fset = c`; S's own classifications use `fset = "C3"` (KIT-185's F_m with C3 for C). When C3
  does not return, S's set falls back to S's own, as KIT-185's did when C did not return (R1-16). The DA rows
  (F_m^filt) carry `fset = "C3"` for every arm.
- Delta, D, complete-data departure and excess features lost as 6.3. Excess features lost is lost(c) - lost(S) on
  F_m^filt per replicate, not oracle-referenced; edges 0.02 / 0.10: **the power band's material and serious edges
  (VAL-67), by VAL-152 (D6)** (A-1; R1-16).
- The candidates' DA also runs on `none` (forced by 6.2 and 6.3: the complete-data departure of FDP and power).
- J per route class: legacy B, C*, FDP, power; declared B, C*, excess lost. "refused" when the route or the candidate's
  guard refused in more than 5 % of the unit's replicates; "not determined" above 5 % not completed (KIT-185's rule).
- Safe rates per (candidate, route class, cell, mechanism) by KIT-185's rule; **the more severe mechanism is the one
  with the lower safe rate; on a tie, the one with the lower realised share at the safe rate** (R1-11; C-5). The same
  rule picks the H-rate threshold from KIT-185's safe-rate table.
- Rules 1 to 4 as 6.4. Rule 2 (d) fails when any deciding cell scored for S has no safe rate for the candidate (R1-10).
  The complete-data departure family is 56 and uses B, C*, FDP and power on both route classes (as 6.5; R1-16).
- **Across-candidate adjustment (6.5; R1-08):** J, the safe rates and the more severe mechanism are built twice, at
  KIT-185's limits and at the K = 6 limits (one-sided 1 - 0.05/6 for "no material effect", 0.05/(456 x 6) for
  "material" and "serious"); rule 1, rule 2 (a) to (d) and rule 3 (a) of the robust verdict read the second set. T1, T2
  and T4a read S's robust rule 1; T3 reads the candidates' robust rule 1 and robust verdict. C7r5 is reported and not in K.

## 5. H-rate and H-kappa (KIT-210 5, 6.6)

- **The hold or warning acts per data set (decided by the owner, VAL-156):** a data set is admitted when every one of its
  assays meets the condition; any assay over the threshold holds or warns the whole data set; the data set's kappa-hat
  is the largest over its assays. (The estimands are per data set; a per-assay hold would need per-assay B, C*, FDP and
  power and is not measured.)
- A data set whose kappa-hat cannot be computed (fewer than k + 1 complete features) is held (R1-16).
- H-rate threshold per cell: KIT-185's `results/safe-rate.tsv`, the realised share at the safe rate of the more severe
  mechanism (section 4's rule, ties included).
- kappa*: fitted on replicates 1-200 of the deciding units, thresholds the 2 % to 100 % quantiles of the data sets'
  kappa-hat, margin 1 - 0.05/K; validated on 201-400 by 6.4's H rule (units with at least 30 admitted data sets each,
  the others pooled per route class); a unit is cell, variant and mask, so the sensitivity variants are never pooled
  (R1-09); out of design on P0lo/P0hi, the sensitivity variants and the run-aligned masks; leave-one-shape-out reported.

## 6. Lifecycle L1 (KIT-210 6.9; VAL-11, VAL-150)

Replicates 401 to 403 of every part-A cell, part B (401-402), the sensitivity panel and the run-aligned masks, in seven
requests on the bound commit (section 8).

**What leaves each job (the whole list):** per script `console.log` (start and finish lines, unit counts; no value) and
these files under the request's `out` directory, each matched by an exact `out_pattern`:

| File | Content |
|---|---|
| `status-<job>.tsv` | per unit, mask and arm (S, each candidate, O, N, mask, assay): completed, or crashed with its error class. **No seconds** (a refusal returns in about 0 s; R1-02) |
| `hashes-<job>.tsv` | per unit: SHA-256 of its full rows (the rows stay on the runner) |
| `checks-<job>.tsv` | pass or fail of E1', E2', E3', E4', E4'', E6, E7, E8, G1, G2 (a check on a masked arm writes pass or fail only; checks on `none` and on constructed cases keep their value) |
| `timing-<job>.tsv` | per unit: seconds per component summed over masks (filter, S, each candidate, kappa-hat, DA, estimands), warnings, DA errors |
| `provenance-<job>.json` | commit, versions, seeds, hashes of the scripts and of the files above |
| `check-e1.tsv` | E1 (KIT-185's arms, sim-k210.R with the candidates off against sim-k185.R): pass or fail per cell, row counts, the names of differing columns |
| `check-e10.tsv` | E10 **with the candidates on** (the counted code path, KIT-185's seed base, replicates 1-3 of the 14 cells; C-4, R1-03): filter counts; B, C*, A1 of S, O, N, S0, O0, N0 and of C3 against KIT-185's C where C3 returned; the DA columns, deciding for lipid/metab, listed as expected (KIT-201) for proteomics (R1-04); pass or fail and the names of differing columns |
| `check-e1pp.tsv` | E1'': pass or fail, the names of differing columns |
| `check-scale.tsv` | the scorer (with R1-08 and R1-09 applied, R1-07) on synthetic rows of the counted run's size: rows, seconds, peak memory, completed or the error |

Nothing else: no refusal count, guard count, features lost, |F(S)|, kappa-hat, admission column or estimand value of
a masked arm or of a candidate's `none` row. The full rows and every intermediate file are written under
`k210-private/` (outside every `out`), read only by the comparison modes in the same job, and never uploaded. The scorer
is not run on L1 rows; its scale test runs on synthetic rows. G1 uses only S's `none` rows. E4' runs in every cell
whose replicate sets are the groups, the D5 cells included, not P5 (R1-12; C-5).

Gates G1 to G3 as KIT-210 6.9. L1's measured figures (named blanks):

- `[L1: seconds per replicate per cell and per candidate]` (from `timing-*.tsv`)
- `[L1: arms completed / crashed per candidate]` (from `status-*.tsv`)
- `[L1: E1, E1', E1'', E2', E3', E4', E4'', E6, E7, E8, E10, G1, G2 pass or fail]`
- `[L1: scorer scale test: rows, seconds, peak memory]`
- `[L1: G3: measured seconds against section 11 per cell; job splits]`

## 7. Budget (KIT-210 11; to be restated from L1's measured seconds)

- L1: **about 27 expected / about 45 ceiling** core-hours (the coordinator's C-3, from the auditor's RC7: indicative,
  scaled from KIT-185's timing files; seven requests of at most 270 minutes of script timeouts each).
- Ticket total (r1.4): about 831 expected / about 1,527 ceiling with C6, C7 and the reported C7r5.
- Counted jobs, splits and timeouts `[L1: from measured seconds, G3]`.

## 8. The blind lane (written, not sent)

Branch `lane-k210-l1` (from `run-tests` 5ef1a49d), files under `.github/run-tests/drafts/` (not `requests/`, so a push
of that branch or a merge of it starts no run-tests run). Seven request drafts, `package_commit` "BOUND_COMMIT" (the
plan job would refuse them as they stand); script timeouts per request (minutes): 240, 255, 270, 195, 195, 175, 240,
because the scripts job stops at 300 minutes and its upload has no `if: always()` (R1-05; C-2):

| Request | Content |
|---|---|
| `k210-l1-1` | lifecycle: proteomics part A with the E7 fault; the D5 cells; the run-aligned masks |
| `k210-l1-2` | lifecycle: lipid/metab part A; the sensitivity panel; part B |
| `k210-l1-3` | E1 on the proteomics part-A units; E1''; E10 (candidates on) |
| `k210-l1-4` | E1 on the lipid/metab part-A units and the sensitivity panel |
| `k210-l1-5` | E1 on part B |
| `k210-l1-6` | E1 on the D5 cells and the run-aligned masks |
| `k210-l1-7` | the scorer's scale test alone |

My check of the workflow (`run-tests.yml`, `run-scripts.sh`), what leaves a scripts job:
(1) `$RUNNER_TEMP/scripts-out` is the only uploaded directory (artifact `scripts`, 30 days), then `results`;
(2) it holds `blas.txt`, `cpu.txt`, `library.tsv`, `r-version.txt`, `scripts.tsv` and per script `console.log`
and `out/`; (3) `out/` is a copy of the script's declared `out` directory with every file not matching `out_pattern`
deleted before the upload; a script with no `out` uploads no files; (4) the job log (visible on the run page) shows the
last 20 lines of each console; (5) `scan.R` withholds token-shaped and build-output files (it is not a blindness
filter). So blindness rests on: each script writing value-free lines to its console; the private rows under
`k210-private/` (never an `out`); each script's `out_pattern` naming only its own files. The KIT-184 failure (R02:
the scorer writing every table into an uploaded folder and a step summary) cannot recur here: the scorer's comparison
and scale modes write one table each to the `out` and everything else under `k210-private/`, and run-tests' step summary
carries only each script's index, path, exit status and seconds (`.github/run-tests/assemble.R:259-261`). A push of
`lane-k210-l1` starts the repository's Module CI (not run-tests); that run is cancelled at once and recorded.

## 9. Merge-sensitive calls (every package call goes through one wrapper; the bound commit's behaviour is what L1 sees)

| Wrapper in `sim-k210.R` | Package call | Pending change that can alter it | What L1 checks |
|---|---|---|---|
| `filterArm` (proteomics) | `updateMissingValueParameters`, `removeRowsWithMissingValuesPercent` | VAL-142 alignment (filter and test on threshold-level values), before the release commit; KIT-187 already merged | E10 (filter counts), G1 |
| `filterArm` (lipid/metab) | `lipidIntensityFiltering`, `metaboliteIntensityFiltering` | **KIT-202** (intensity-filter samples; on the bound commit's list, C-1). Inferred no-op for the simulated tables (their only non-sample columns are the character `feature_id` and `name`; the filter counts numeric columns, `R/func_lipid_qc.R:47`) | E10 |
| `buildArm` | `normaliseBetweenSamples` (sensitivity "loess" only) | none known | G2 |
| `fitShipped`, `fitC6` | `ruvIII_C_Varying` S4 methods (proteomics, lipid, metab), `removeRowsWithMissingValuesPercent` | KIT-188's ruling if it changes the RUV-III-C path (E1's list) | E2', E3', E10 |
| `postProcess` | `restoreRuvFeatureLevels`, `removeRowsWithMissingValuesPercent` | VAL-142 (the filter) | E2' |
| `daTable` (proteomics) | `differentialAbundanceAnalysis` | **KIT-201**: TREAT takes the moderated test's trend and robust (`R/func_prot_da_model_stats.R`); the default (option `MultiScholaR.treat_moderation_rule` unset) is the new rule, so the wrapper needs no change; changes FDP and power at cutoff 0. **KIT-200** (a protein with no value in the legacy DA; on the list, C-1): inferred no-op for part A, which passes only rows with a finite value; part B's route DA less certain | E10 lists proteomics DA differences as expected |
| `daTable` (lipid/metab) | `runLipidsDA`, `runMetabolitesDA` (`da_lipids_long`, `da_metabolites_long`, `failed_assays`) | **KIT-196** (features with no value; `fitLipidDaAssayModel`, `applyLipidDaContrasts`, `runLipidsDA`, metab the same) and **KIT-201** (TREAT in `applyLipidDaContrasts`); the long table's columns `P.Value`, `fdr_qvalue`, `assay`, the id column must stay as read | G2, E10 (deciding for lipid/metab) |
| `c4bFactors` | `computeRuvCovariates`; `computeOmicsRuvAssayFactors` | KIT-190 touches `R/func_prot_norm_optimization_helpers.R` (control selection, `getNegCtrlProtAnovaHelper`), not these functions | E4', E4'' |
| `fitC4b` | `limma::removeBatchEffect` | moves with the lock's limma (R1-14) | E4'' |
| `declaredW` (E4') | `computeProtRuvModelFactors`, `runOmicsRuvFactors` | KIT-188's ruling (if any) | E4' |
| `nipalsRuns` | `imputeNipalsMatrix` | none known | E2' (C7 = S on `none`) |
| `replicateSets` | `resolveProtRuvReplicateBlocks`, `getRuvIIIReplicateMatrixHelper` | none known | G1 |
| `protObject`, `omicObject` | the S4 constructors `ProteinQuantitativeData`, `LipidomicsAssayData`, `MetaboliteAssayData` | a slot change at the bound commit breaks them (R1-14) | G2 (every unit completes) |
| part B `candidateB` (declared route, lipid/metab) | `restoreRuvFeatureLevels` in place of the route's `restoreOmicsRuvAssayLevels` | the same arithmetic by its documentation; reported part only (R1-14) | none (reported) |
| part B `routeB` | `diaNormalizationOperation`, `omicsRuvOperation` kernels, the bindings' `extract_best_k`, `extract_controls` | **KIT-190** (the legacy lipid/metab control selection; `omicsRuvOperation`); KIT-188 | part B status (reported) |
| `c7Components` | reads `object@args$ruvCancor$num_components_to_impute` | KIT-190 adds code to `R/func_metab_norm.R`; line 298's rule must be re-read at the bound commit | declaration re-read |

## 10. Provisional hashes and open items before the declaration

- Provisional SHA-256 at the commit that adds this revision: `sim-k210.R` 4e520db7a41404fd24ba3356ec50d193f82fb08fb682872ab9cfd509c12403fe, `ruviiic-ref210.R` 61b86b48c67c8d9b2dd36416242f44a085eefd921e430c2cac7e09c9b4f89add,
  `score-k210.R` 25dcd8bc2fb8a45aedd2a795f412e1c5858ad9096467358707d99ebdd7d9cacf.
- E1's bound commit; the rebase of this branch onto it; the scripts' final hashes.
- The auditor's check of this draft and of L1's outputs.
- Poulos et al. 2020 not read: a limit (KIT-210 14).
