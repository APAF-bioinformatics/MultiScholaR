# KIT-210 declaration: DRAFT (not the declaration)

**Status: draft, written 2026-11-24 against v0.9-staging 53713252 with no run of any kind.** It is not `declaration.md`
and binds nothing. Every measured figure is a named blank `[L1: ...]` that the blind lifecycle run L1 fills. The
script hashes below are **provisional**: the scripts are untested (no R on the laptop) and will change at L1 under
VAL-54 if L1 shows a mechanics defect, and again if the bound commit changes a call's signature (section 9). The
declaration is committed alone, after L1 and the auditor's check of this draft and of L1, before any counted run
(KIT-210 10, A5). Design: KIT-210 revision 1.3 (`tickets/kit/KIT-210.md`), owner's decisions VAL-152.

## 1. Bound code

- Package: the bound commit of E1, `[E1: the first v0.9-staging commit with KIT-190, KIT-196, KIT-201 and the
  pre-release filter / RUV-path / legacy-DA changes (today the VAL-142 alignment; KIT-188's ruling if it changes code on
  these paths)]`. Not available on 2026-11-24 (none of KIT-190, 196, 201 is merged on 53713252).
- Every request's `package_commit`: this ticket's branch `kit-210-ruviiic-remedy` rebased or merged onto the bound
  commit, adding only `docs/qa/workflow-execution/kit/kit-210/` `[E1: commit]`.
- Scripts (SHA-256, provisional, at the draft commit; recomputed in the declaration):
  - `sim-k210.R` 52f3cbc1887ebc49aa166236693f583fba618badd91131c24016ade8ec5149b8 (provisional)
  - `ruviiic-ref210.R` d6868348de8fea0aa4d25abe32d5be80cc02847423eec5207e9282454acbb7f3 (provisional)
  - `score-k210.R` 3bb349b022f5c82bcbd5aa8af1079b8cd57625b82711e7954f1c19bed034825c (provisional)
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
| C7 | `ruviiicRef(directions = "nipals", ncomp)` with `imputeNipalsMatrix`; **ncomp = the route's own count when it chooses k**: proteomics 5 (`R/func_prot_norm_optimization_helpers.R:679`), lipidomics/metabolomics the object's stored `ruvCancor` argument, otherwise 2 (`R/func_metab_norm.R:298`; the simulated objects store none: 2) |
| C7r5 (reported) | as C7 with ncomp = 5; computed for lipidomics/metabolomics only (on proteomics it is C7) |
| O, N | oracle, none (KIT-185's) |
| C (L1 only, E1) | KIT-185's unguarded complete-feature comparator, run only with `K210_CANDIDATES=off` |

Per data set and assay (rows `corr = "assay"`): realised missing share after the filter, kappa-hat (A-8: Q-hat with the
co-missing term, P over all runs, delta-hat scaled by n / n_c at the part-A k), its perturbation norm and eigengap, the
sine of the largest principal angle between S's directions (zero-filled cross-product) and C3's, the number of complete
features.

## 4. Estimands, feature sets, primary endpoint, rules (KIT-210 6.3 to 6.5; `score-k210.R`)

As KIT-210 r1.3 sections 6.3 to 6.5, implemented as follows (the choices below were not fixed by the ticket; each is
named for the auditor, **A-n**):

- F_m(c): the features of F_m^filt finite in c's output on m, in S's on m and in S's on `none`. Rows of c and of S, O, N,
  S0, O0, N0 on F_m(c) carry `fset = c`; S's own classifications use `fset = "C3"` (KIT-185's F_m with C3 for C). The DA
  rows (F_m^filt) carry `fset = "C3"` for every arm.
- Delta, D, complete-data departure and excess features lost as 6.3. **A-1**: excess features lost is lost(c) - lost(S)
  on F_m^filt per replicate, not oracle-referenced; its material edge is 0.02 (D6) and **its serious edge 0.10** (the
  power band's, which the ticket does not state).
- **A-2**: the candidates' DA also runs on `none` (needed for the complete-data departure of FDP and power).
- J per route class: legacy B, C*, FDP, power; declared B, C*, excess lost. "refused" when the route or the candidate's
  guard refused in more than 5 % of the unit's replicates; "not determined" above 5 % not completed (KIT-185's rule).
- Safe rates per (candidate, route class, cell, mechanism) by KIT-185's rule; the more severe mechanism = the lower safe
  rate (**A-3**: ties take the first mechanism listed, which gives the same rate).
- Rules 1 to 4 as 6.4; families as 6.5 (456; departure 56; robustness at K = 6: S, C3, C4a, C4b, C6, C7; C7r5 is
  reported and not in K).
- Terminal outcome per route class as 6.9; S's "safe under stated conditions" taken at the K-adjusted limits (robust),
  as 6.5 asks of any verdict the owner acts on (**A-4**).

## 5. H-rate and H-kappa (KIT-210 5, 6.6)

- **A-5**: a lipidomics data set (two assays) is admitted when every assay meets the condition; the data set's
  kappa-hat is the largest over its assays (the correction is per assay, the estimands are per data set).
- H-rate threshold per cell: KIT-185's `results/safe-rate.tsv`, the realised share at the safe rate of the more severe
  mechanism (**A-6**: the mechanism with the lower safe rate in that table).
- kappa*: fitted on replicates 1-200 of the deciding units, thresholds the 2 % to 100 % quantiles of the data sets'
  kappa-hat, margin 1 - 0.05/K; validated on 201-400 by 6.4's H rule (units with at least 30 admitted data sets each,
  the others pooled per route class); out of design on P0lo/P0hi, the sensitivity variants and the run-aligned masks;
  leave-one-shape-out reported.

## 6. Lifecycle L1 (KIT-210 6.9; VAL-11, VAL-150)

Replicates 401 to 403 of every part-A cell, part B, the sensitivity panel and the run-aligned masks, in three requests
on the bound commit (the lane in section 8).

**What leaves each job (the whole list):** per script `console.log` (start and finish lines, unit counts; no value) and
these files under the request's `out` directory, each matched by an exact `out_pattern`:

| File | Content |
|---|---|
| `status-<job>.tsv` | per unit, mask and arm (S, each candidate, O, N, mask, assay): completed, or crashed with its error class; seconds |
| `hashes-<job>.tsv` | per unit: SHA-256 of its full rows (the rows stay on the runner) |
| `checks-<job>.tsv` | pass or fail of E1', E2', E3', E4', E4'', E6, E7, E8, G1, G2 (a check on a masked arm writes pass or fail only; checks on `none` and on constructed cases keep their value) |
| `timing-<job>.tsv` | per unit: seconds per component (filter, S, each candidate, kappa-hat, DA, estimands), warnings, DA errors |
| `provenance-<job>.json` | commit, versions, seeds, hashes of the scripts and of the files above |
| `check-e1.tsv`, `check-e10.tsv` | E1 and E10: pass or fail per cell, row counts, the names of differing columns (no value) |
| `check-e1pp.tsv` | E1'': pass or fail, the names of differing columns |
| `check-scale.tsv` | the scorer on synthetic rows of the counted run's size: rows, seconds, peak memory, completed or the error |

Nothing else: no refusal count, guard count, features lost, |F(S)|, kappa-hat, admission column or estimand value of
a masked arm or of a candidate's `none` row. The full rows and every intermediate file are written under
`k210-private/` (outside every `out`), read only by the comparison modes in the same job, and never uploaded. The scorer
is not run on L1 rows; its scale test runs on synthetic rows. G1 uses only S's `none` rows.

Gates G1 to G3 as KIT-210 6.9. L1's measured figures (named blanks):

- `[L1: seconds per replicate per cell and per candidate]` (from `timing-*.tsv`)
- `[L1: arms completed / crashed per candidate]` (from `status-*.tsv`)
- `[L1: E1, E1', E1'', E2', E3', E4', E4'', E6, E7, E8, E10, G1, G2 pass or fail]`
- `[L1: scorer scale test: rows, seconds, peak memory]`
- `[L1: G3: measured seconds against section 11 per cell; job splits]`

**Note for the auditor (A-7):** the coordinator's list of L1 outputs names "formation and refusal counts". KIT-210 6.9
forbids refusal and guard counts of masked arms in L1, so the draft writes only completed / crashed per arm; a candidate
that formed is "completed". If refusal counts are wanted, 6.9 must change first.

## 7. Budget (KIT-210 11; to be restated from L1's measured seconds)

- Decided design about 811 expected / 1,500 ceiling; with C7r5 (r1.3) about 819 / 1,514.
- L1 row of 11: 15.3 / 32.0. **Indicative (not measured):** the candidates roughly double KIT-185's lifecycle units
  (KIT-185 L1: 6.2 core-hours in three jobs), and E1 (both sim scripts on every L1 unit), E10 (replicates 1-3 of 14
  cells, KIT-185's arms) and the scale test add about 10; about 27 expected, inside the 32.0 ceiling but above the 15.3
  expected `[L1: measured]`.
- Counted jobs, splits and timeouts `[L1: from measured seconds, G3]`.

## 8. The blind lane (written, not sent)

Branch `lane-k210-l1` (from `run-tests` 5ef1a49d; commit 3f8f988f), files under `.github/run-tests/drafts/` (not `requests/`, so a push of that
branch or a merge of it starts nothing). Three request drafts, `package_commit` "BOUND_COMMIT" (the plan job would
refuse them as they stand):

- `k210-l1a`: proteomics part A (P0, P0r, P1, P2, P3, P4, P5, P0lo, P0hi; with `K210_FAULT=1`), E1 on the same units,
  E1'' and the scorer scale test.
- `k210-l1b`: lipidomics and metabolomics part A (L1, L2, M1, M2), the sensitivity panel (P0, M2), part B, E10 (the 14
  part-A cells, replicates 1-3, KIT-185's seeds and arms).
- `k210-l1c`: the D5 cells (P4b, L1b, M1b) and the run-aligned masks (P0, P1, L2, M2), E1 on those units.

My check of the workflow (`run-tests.yml`, `run-scripts.sh` at run-tests 9db9fed2), what leaves a scripts job:
(1) `$RUNNER_TEMP/scripts-out` is the only uploaded directory (artifact `scripts`, 30 days), then `results`;
(2) it holds `blas.txt`, `cpu.txt`, `library.tsv`, `r-version.txt`, `scripts.tsv` and per script `console.log`
and `out/`; (3) `out/` is a copy of the script's declared `out` directory with every file not matching `out_pattern`
deleted before the upload; a script with no `out` uploads no files; (4) the job log (visible on the run page) shows the
last 20 lines of each console; (5) `scan.R` withholds token-shaped and build-output files (it is not a blindness
filter). So blindness rests on: each script writing value-free lines to its console; the private rows under
`k210-private/` (never an `out`); each script's `out_pattern` naming only its own files. The KIT-184 failure (R02:
the scorer writing every table into an uploaded folder and a step summary) cannot recur here: the scorer's comparison
and scale modes write one table each to the `out` and everything else under `k210-private/`, and run-tests' step summary
carries only each script's index, path, exit status and seconds (`.github/run-tests/assemble.R:259-261`).

## 9. Merge-sensitive calls (every package call goes through one wrapper; the bound commit's behaviour is what L1 sees)

| Wrapper in `sim-k210.R` | Package call | Pending change that can alter it | What L1 checks |
|---|---|---|---|
| `filterArm` (proteomics) | `updateMissingValueParameters`, `removeRowsWithMissingValuesPercent` | VAL-142 alignment (filter and test on threshold-level values), before the release commit; KIT-187 already merged | E10 (filter counts), G1 |
| `filterArm` (lipid/metab) | `lipidIntensityFiltering`, `metaboliteIntensityFiltering` | none known among KIT-190/196/201; KIT-202 (intensity-filter samples) is not merged and not in E1's list: the coordinator should say whether it is pre-release | E10 |
| `buildArm` | `normaliseBetweenSamples` (sensitivity "loess" only) | none known | G2 |
| `fitShipped`, `fitC6` | `ruvIII_C_Varying` S4 methods (proteomics, lipid, metab), `removeRowsWithMissingValuesPercent` | KIT-188's ruling if it changes the RUV-III-C path (E1's list) | E2', E3', E1 |
| `postProcess` | `restoreRuvFeatureLevels`, `removeRowsWithMissingValuesPercent` | VAL-142 (the filter) | E2' |
| `daTable` (proteomics) | `differentialAbundanceAnalysis` | **KIT-201**: TREAT takes the moderated test's trend and robust (`R/func_prot_da_model_stats.R`); the default (option `MultiScholaR.treat_moderation_rule` unset) is the new rule, so the wrapper needs no change; changes FDP and power at cutoff 0 | E10 lists FDP/power differences as expected |
| `daTable` (lipid/metab) | `runLipidsDA`, `runMetabolitesDA` (`da_lipids_long`, `da_metabolites_long`, `failed_assays`) | **KIT-196** (features with no value; `fitLipidDaAssayModel`, `applyLipidDaContrasts`, `runLipidsDA`, metab the same) and **KIT-201** (TREAT in `applyLipidDaContrasts`); the long table's columns `P.Value`, `fdr_qvalue`, `assay`, the id column must stay as read | G2 (no DA stop), status |
| `c4bFactors` | `computeRuvCovariates`; `computeOmicsRuvAssayFactors` | KIT-190 touches `R/func_prot_norm_optimization_helpers.R` (control selection, `getNegCtrlProtAnovaHelper`), not these functions | E4', E4'' |
| `declaredW` (E4') | `computeProtRuvModelFactors`, `runOmicsRuvFactors` | KIT-188's ruling (if any) | E4' |
| `nipalsRuns` | `imputeNipalsMatrix` | none known | E2' (C7 = S on `none`) |
| `replicateSets` | `resolveProtRuvReplicateBlocks`, `getRuvIIIReplicateMatrixHelper` | none known | G1 |
| part B `routeB` | `diaNormalizationOperation`, `omicsRuvOperation` kernels, the bindings' `extract_best_k`, `extract_controls` | **KIT-190** (the legacy lipid/metab control selection; `omicsRuvOperation`); KIT-188 | part B status (reported) |
| `c7Components` | reads `object@args$ruvCancor$num_components_to_impute` | KIT-190 adds code to `R/func_metab_norm.R`; line 298's rule must be re-read at the bound commit | declaration re-read |

## 10. Open items before the declaration

- E1's bound commit; the rebase of this branch onto it; the scripts' final hashes.
- The auditor's check of this draft and of L1's outputs (the coordinator's instruction).
- A-1 to A-7 above.
- Poulos et al. 2020 not read: a limit (KIT-210 14).
