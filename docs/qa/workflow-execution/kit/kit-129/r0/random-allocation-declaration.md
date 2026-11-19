# KIT-129 R0 part B, supplement: the unrelated cells with the factor allocated at random (declared 2026-11-19, before any of its data exist)

Authority: decision M37 (user, 2026-11-19), with M36, in `docs/qa/workflow-execution/kit/kit-124/design-review-20261118.md`
section 6 (origin/performance-v2 `74bdf127`).

## Purpose

The thresholds of the check that an unwanted factor follows the groups (now presented as a detection limit, M36) were
validated as firing in about 5% (5% point) or about 10% (10% point) of "unrelated" data sets (M24, M27). In the
classic proteomics, transcriptomics and lipid/metab simulators those data sets had an unwanted factor with EXACTLY zero
association with the groups in the sample. Real unblocked experiments have chance imbalance. This supplement measures,
with the factor allocated at random (association zero in expectation only): (a) how often the estimated association
exceeds the validated limit, and (b) whether false discoveries stay at target there.

**Descriptive.** No pass mark, no verdict. It changes no threshold, no verdict and no earlier record, whatever it shows.

## How the association was built before (sim-k129-r0b at `5a1ee379`)

| Simulator | Construction | Lines |
|---|---|---|
| S1 classic proteomics | a normal draw residualised on the group model (`lm.fit(basis, e)`, basis = `model.matrix(~ group)`) and scaled; `w = rho g + sqrt(1 - rho^2) e`: exact | `s1_classic.R` 274-277 |
| S3 transcriptomics | a normal draw residualised on the group model (`qr.resid`), orthonormalised; `W1 = rho u + sqrt(1 - rho^2) z1`: exact | `s3_trx.R` 213-216 |
| S4 lipid/metab | `residuals(lm(rnorm(n) ~ group))`; `f = rho z + sqrt(1 - rho^2) e`: exact. The drift's run order is the most balanced of 2,000 random orders (line 282) | `s4_lipid_metab.R` 288-289 |
| S2 DIA_LIMPA (not rerun) | plain normal noise plus `rho` times the group score: rho in expectation only | `s2_limpa.R` 68, 305-306 |

## The variant (commit `c2bd327e`)

New cells in S1, S3 and S4, the existing cells unchanged:

| Cell | Index (seed) | Counterpart | Class |
|---|---|---|---|
| C0R | 22 | C0 (centre, 6 per group) | core |
| C4R | 23 | C4 (no true effects) | core |
| S1R | 24 | S1 (3 per group, factor uncorrelated: the unrelated counterpart of S11) | stress |

The factor's scores are the same normal draw as in the counterpart, standardised (sample mean 0, variance 1, so the
strength is unchanged) but NOT made orthogonal to the group model. Everything else (effects, sample sizes, missingness,
loadings, controls, arms, streams) is the counterpart's. In S4 the drift's run order stays chosen as in every cell
(balanced): only the sample-wise factor is randomised. C3 has no factor and is not rerun. The new cells are never part
of `all`. Existing cells: the code path is textually unchanged (each change is behind the cell's `random` flag).

**Mode.** The association mode runs no test, so it cannot give false-discovery measures. The run therefore uses the
study mode with K124_ARMS, and a random-allocation cell also writes `s<i>-assoc.tsv` (columns of `association.R`) with
the statistic of the W each arm actually tested with, and the realised true association. One run gives both.

## Lifecycle (done before this declaration; runs outside the declared range)

- Run 44419142081 (commit `120c689c`): C0 and C0R, runs 1 to 4, each simulator, the arms below. All 6 jobs succeeded
  on their first attempt (about 3.7 core-hours).
- Existing cell unchanged: C0 rows of arms 1 and 2, runs 1 to 4, against the earlier newarms run 44200924581 (same
  seeds): S3 and S4 identical in every column but seconds. S1: every count, null rate and FDP identical; `bias` and
  `rmse` differ in the 15th significant digit. The S1 code path for C0 is textually the same and the runner's R,
  namespaces and BLAS library are the same, so this is read as OpenBLAS choosing different kernels per CPU model.
- C0R writes `s<i>-assoc.tsv` with status, k, group_association and true_association for every arm and run (S4 per
  assay).
- Scoring lifecycle: run 44421002832 ran `random-allocation.R` on that run (purpose `lifecycle`, C0R, runs 1 to 4)
  and wrote every table below.

## Declaration

- **Simulators:** S1, S3, S4. **Cells:** C0R, C4R, S1R. **Runs:** 500 per cell, runs 3001 to 3500 (`first_run` 3001).
  Seeds: seed = 100000 + 10000 x cell index + run, so 323001-323500 (C0R), 333001-333500 (C4R), 343001-343500 (S1R).
  No cell with index 22 to 24 existed before this supplement, so these seeds have never been used. Arm streams as in
  the study (seed x 100 + arm).
- **Arms:** the designs with a validated threshold, plus arm 1 (no correction) and arm 2 (ceiling: the true factor in
  the model). S1: 1, 2, 5, 7, 8. S3: 1, 2, 7. S4: 1, 2, 7, 8.
- **Thresholds, used as they are (fires = strictly above):**

| Simulator | Design | Threshold | Set at | Source |
|---|---|---|---|---|
| S1 | 5 | 0.16376377 | 5% point | threshold-followup-S1.tsv (M27) |
| S1 | 7 | 0.058077853 | 5% point | threshold-followup-S1.tsv (M27) |
| S1 | 8 | 0.070068989 | 5% point | threshold-followup-S1.tsv (M27) |
| S3 | 7 | 0.038965628 | 10% point | threshold-by-design-S3.tsv (M24) |
| S4 | 7 | 0.088611623 | 10% point | threshold-by-design-S4.tsv (M24) |
| S4 | 8 | 0.19494193 | 5% point | threshold-followup-S4.tsv (M27) |

- **Quantities, fixed now** (per simulator, cell, and C0R+C4R together; per design):
  - the TRUE in-sample association (one per run): quartiles, median, 90th percentile, and the share of runs above the threshold;
  - the ESTIMATED largest association (rows with a W): quartiles, median, 90th percentile;
  - the share of non-crashed rows above the threshold (rows without a W cannot fire and stay in the denominator;
    S4 over both assays, as in the calibration);
  - type I at 0.05 and 0.01 (mean `null_p05`, `null_p01`) and the mean false discovery proportion at the 5% target
    (mean `fdp`), `full` subset, crashed rows left out, refused rows scored as the run's no-RUV result (as aggregate.R),
    for the design and for arms 1 and 2 on the same runs: overall; split by whether the run exceeded the design's
    threshold (S4: in either assay); and by thirds of the true association (cut at its 1/3 and 2/3 quantiles).
- **Scoring:** `random-allocation.R`, sha256 `281b70d15cf65dcf0cf4f758138a2ae1e4e001ab053fe807cc203543327033c3`,
  run once with purpose `score` by `.github/workflows/k129-score.yml` (sha256
  `b469591e5f93cd8e596d7df2b7eaa5b441a75e509f4fd3c5d8f7fdee8563b6dc`) on a hosted runner. It refuses if its results are
  already committed (read once). Raw tables are recorded by sha256 and not committed. Only aggregates are committed.
- **Simulators** at this declaration: `s1_classic.R` `88a5a426...`, `s3_trx.R` `f50ae13b...`, `s4_lipid_metab.R`
  `6a26a3e0...` (sha256).
- **Budget:** at most 40 core-hours for the whole supplement, the lifecycle included (job duration x 4). One attempt
  per job. One rerun of a job that fails for infrastructure, recorded. An unfavourable valid result is reported as it is.

## Order, shown by separate pushes of `sim-k129-r0b`

1. `c2bd327e` the variant, workflows and scoring script; `120c689c` lifecycle launch; `d3538ff4` scoring lifecycle.
2. This file, pushed alone.
3. Then the launch params (label `r0b-random`, runs 3001 to 3500, the three cells, the arms above).
4. After the run, the scoring params (purpose `score`), then the results.
