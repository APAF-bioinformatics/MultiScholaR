# KIT-211 lifecycle L1 (blind)

Run 45713532325 (request `kit211-l1-15829c8f`), package commit 15829c8fc8a505b0b13800a5cd1c8b105c72049d (the scripts and
declaration as committed; `git diff --stat 894b9cc8 15829c8f -- R/ src/ inst/ tests/ DESCRIPTION NAMESPACE` is empty), R 4.6.1,
reference BLAS, runner AMD EPYC 9V74 (4 cores). Job 23.4 min (script 832 s) = 1.56 core-hours against 2.00 expected and 3.67
ceiling. Run once; every check passed at the first attempt; predictions (declaration.md) all held.

Blind: the L1 rows and scores stayed in the runner's temporary directory. Uploaded: `l1/checks.tsv`, `l1/timing.tsv`,
`l1/hashes.tsv`, `l1/console.txt`, plus `l1/g0.txt` and `l1/e0.txt` (made locally, present in the checkout).

| Check | Declared expectation | Observed |
|---|---|---|
| E0 | pass | pass (local, `l1/e0.txt`) |
| E1 | 120/120; DA flag expected FALSE | 120 of 120; DA-dependent flag FALSE (as predicted: KIT-201 moves p-values) |
| E2, E3, E7 | 60/60 | 60 of 60 each |
| E9 | 60/60, warnings 0 | 60 of 60, warnings counted 0 |
| E10 | 60/60, deviation 0 | 60 of 60, maximum absolute deviation 0 |
| E8 | 2/2 | 2 of 2 |
| E4 | 56 sub-checks | 56 of 56 (11 scenarios: outcomes K1 to K5, N0 twice, Stop (F4 and F1), a withheld candidate (F3, F5), the join with an incomplete range, empty resamples, hand estimates and jackknife SEs) |
| E5 | identical | two runs exit 0; 11 files each; identical |
| E6 | killed and continued | killed with 1 of 4 complete; continued with 3; rows and metrics identical to uninterrupted; each replicate once with six masks |
| G2 | planned <= 150 | deciding-cell job 41.0 min, Ag 38.4 min (both below 150; timeouts 120 and 139 stand; no split) |

Timing (`l1/timing.tsv`): 16.4 s per replicate-mask all arms together (deciding cells, 4 workers), 15.4 s for Ag; the planning model
assumed 24.6 s expected, so the counted jobs are expected to cost less than the ticket's figures. No result, count, rate or
threshold left L1 beyond units counts and the identity's deviation. One departure from the ticket's wording, declared
beforehand (declaration.md, resolution 3): L1's two ranges are two directories (`rows/r9101`, `rows/r9102`) each holding every cell.
