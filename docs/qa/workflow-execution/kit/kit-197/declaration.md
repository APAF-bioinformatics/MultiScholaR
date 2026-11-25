# KIT-197 declaration (committed alone, after L1 and before any counted job)

Ticket `tickets/kit/KIT-197.md` revision 1.2 (the plan audit and its re-check applied; the owner's answers D0 to D5,
methods ledger VAL-121). Everything below is fixed before S_a is sent. Between S_a and S_b nothing declared changes
(ticket 6, R1-03): the scripts, seeds, ranges, estimands, margins, rules and the owner's answers are those committed
here, and S_b runs the hashes below. A defect found in S_a (a VAL-54 instrument defect, or the SC′ re-scoring) is fixed
under VAL-54's conditions in a dated section committed with its diff before S_b, and may not change an estimand, a
margin, a rule, a seed or the deciding range.

## Bindings

| Item | Value |
|---|---|
| Package code | `v0.9-staging` at 19997323dd9a0efac5830696c7d3ac617560dae6 (the branch point; `git diff --stat e6c42d1d 19997323 -- R/ src/ DESCRIPTION NAMESPACE inst/config` lists only five `R/*integration*` files, as at 701458be) |
| Scripts commit | 5a702b7b6484e97e03ac515b2dd2c8b0b130c80e (L1 ran on it, run 45316713443, all checks pass (lifecycle.md)) |
| `sim-k197.R` SHA-256 | 71001de9b7b1c547e1e27fd92fb1b57a067121e336fa09e134f3317f1fa7c599 |
| `score-k197.R` SHA-256 | f546947985617c5a6c0e79edd2d2f4a062db2ed0bea922b9927af6a1a9e9e426 |
| `l1-k197.R` SHA-256 (lifecycle only) | 50d55e99c5eb99a4534a4ea93c39aca9fbe4738ca6a783e086b69cd78294da5e |
| S2's script `kit-187/filter-effect.R` | SHA-256 1147a5aa480861581210f1a088c06731e046b20e15d8cdcbe6775a5f58906c9c (4a58a3be, unchanged at 6416dbd6 and 19997323) |
| Every request | `package_commit` = this declaration's commit; `base_commit` "none"; `tests` []; `blas` "reference"; `extra_packages` []; `out_pattern` `"\\.(tsv|gz|json|txt|md)$"` |

E0 (local text check, `e0-k197.py`, no R): the verbatim blocks of `sim-k197.R` equal `filter-effect.R` at 6416dbd6
line for line, except line 215 (`daTable` also keeps `logFC = tab$logFC, t = tab$t, ave = tab$AveExpr`). Block hashes
(SHA-256 of the block's lines as copied, newline-terminated):

| `filter-effect.R` lines | What | SHA-256 of the copy |
|---|---|---|
| 23-234 (215 changed) | KIT-185's generator, masks, object, DA call | 95b0e2e984d2c1cb314a51fc7ee2d43ba4b2dba085f06db8870aa5bbaa39c9d6 |
| 237 | `%\|\|%` | d98bc0062d8d7b73b772c17e1677237a9b7dd464b0337a8735bc115e576842d8 |
| 244-246 | `MIN_REPS`, `MIN_GROUPS`, `PERCENTILE` | a8e791a9206f191703d42cb03f3aee4acf69e511ee678578ba2afc823e3aa13a |
| 248-249 | `k187Abort` | 7e8b32362b94aacf3ead8d9c745615cf85230f9e9da79f592fd3d3810dfc0949 |
| 263-270 | `filterKept` | dd5c8597eb2a34fcb218ba6991322d0f35bd8ab9abea620315aec3103414982e |
| 276-286 | `perGroupKept` | 7ae2657d940ad68963484bce6490dca10ca2eaee6d95b7cd04151249a3a8097f |

## Design

- Cells A (2 groups of 3), B (4 of 3), C (5 of 3), D (3 of 4), indices 201 to 204; 3,000 proteins; strength 0.
- Masks vd5, vd10, vd25 (value-dependent, S2's) and mcar5, mcar10, mcar25 (the replicate's own `V < rate`); the six
  masks of a replicate share its values.
- Seeds: `set.seed(185000000 + 10000 × idx + r, "Mersenne-Twister", "Inversion", "Rejection")`, idx = 200 + i, so
  187000000 + 10000 i + r.
- Ranges: S_a r = 1 to 100 (S2's draws; gate G1 and the descriptive decomposition; no label); S_b r = 101 to 200 (the
  new draws; they alone decide); L1 r = 9001, 9002 (and 9001 to 9004 for E6; never scored as evidence).
- H = the package's filter at the bound commit (`filterKept`, minimums 2 and 2, percentile 1); O = the old rule
  in-script (`perGroupKept` with (m, g) = (3, 2) for A, B, C and (2, 3) for D); Δ = H \ O; Ref = O with no NA in
  `control` or `treated`; the whole list H.
- DA on `masked[H, ]` (S2's call); a protein is tested when its p-value and `logFC` are finite, called when q < 0.05;
  t_r = the largest p called in (cell, mask, r) (0 with no call); ℓ = 1 when `control` or `treated` holds a finite
  value at or below the filter's threshold (1st percentile, type 7, of every finite value of the masked matrix);
  s = +1 (NA in `treated` only), −1 (NA in `control` only), 0 otherwise; strata one_treated, one_control, both,
  other_only, none.

## Estimands (ticket 5.3, as `score-k197.R` computes them)

All ratios are ratios of sums over the replicates of a setting, or of the pool. "Tested" everywhere.

- R1: α_S = Σ null p < 0.05 / Σ null tested (S = Δ, O, Ref); **ρ_O = α_Δ / α_O** (deciding for B); ρ_Ref; κ_S = α_S / 0.05;
  the same at p < 0.01; c_S = Σ null called / Σ_r n0_S,r t_r.
- R2: m6_NA = Σ s × logFC / count over tested null proteins of Δ with s ≠ 0 and ℓ = 0; m6_ℓ the same with ℓ = 1.
- R3: FDR_S = V_S / (V_S + T_S); FDR_S^cal = Σ_r n0_S,r t_r / (Σ_r n0_S,r t_r + T_S); **E_comp = FDR_Δ^cal − FDR_O^cal**
  (deciding for P); E_cal = (FDR_Δ − FDR_Δ^cal) − (FDR_O − FDR_O^cal); E_total = FDR_Δ − FDR_O. With
  f(π0, pow, t) = π0 t / (π0 t + (1 − π0) pow), π0_S = Σ n0 / Σ (n0 + n1), pow_S = T_S / Σ n1, t̄_S = Σ n0 t / Σ n0:
  power first: power part f(π0_O, pow_Δ, t̄_O) − f(π0_O, pow_O, t̄_O), null-share part f(π0_Δ, pow_Δ, t̄_Δ) −
  f(π0_O, pow_Δ, t̄_O); null share first: null-share part f(π0_Δ, pow_O, t̄_Δ) − f(π0_O, pow_O, t̄_O), power part
  f(π0_Δ, pow_Δ, t̄_Δ) − f(π0_Δ, pow_O, t̄_Δ). Reported beside: revision 1.0's Ref-rate counterfactual (Δ's tested
  nulls called at Ref's own null call rate of the replicate; E_null = FDR_Δ − FDR_Δ^Ref, E_rest = FDR_Δ^Ref − FDR_Ref),
  FDR_Δ − FDR_H, S2's mean-of-ratios `added_fdp` (mean over replicates with an added set), π0, power, t̄.
- R4: M5_S = mean (logFC − β) over tested planted proteins of S (Δ, O), M5_Δ − M5_O, by direction, attenuation
  A_S = mean sign(β)(logFC − β), and the missing-side-signed mean over planted proteins of Δ with s ≠ 0.
- R5: ρ_O, E_cal, E_comp, m6_NA, m6_ℓ, A_Δ under vd minus under mcar, same cell and rate, same resamples; and for
  the pools.
- R6: D = calls with s × logFC > 0 / calls, over null proteins of Δ with s ≠ 0; all, ℓ = 0, ℓ = 1.
- R7: within Δ by stratum (and by ℓ): ρ_O, tested nulls per replicate, mean residual variance of the null proteins
  (pooled within-group variance of the observed values); the same variance for O's nulls by stratum.
- Channel share (printed for a label bias or both): (n_p05,Δ,ℓ=1 − α_O n_Δ,ℓ=1) / (n_p05,Δ − α_O n_Δ) over tested
  null proteins.

## Uncertainty

- Bootstrap: B = 10,000; `set.seed(197000001, "Mersenne-Twister", "Inversion", "Rejection")` once; for cells A, B, C,
  D in that order `sample.int(R, B × R, replace = TRUE)` filled by row into a B × R index matrix; one index set per
  cell, shared by its six masks and all statistics (R1-01); a pool sums its settings under their cells' common
  resamples. Percentile limits, `quantile(type = 7)`.
- Jackknife: delete one replicate of one cell with all its masks; units are the cell's R replicates for a setting and
  all 4R for a pool; SE = sqrt((n − 1)/n Σ (θ_(i) − θ̄)²). Jackknife-normal limits (estimate ± z SE at the same levels)
  give a second label, printed; the decision uses the bootstrap limits.
- Every table: R, estimate, bootstrap SE, jackknife SE, 2.5 % and 97.5 % limits, count of non-finite resamples.

## Labels, decision and outcome (ticket 5.6)

| Component | Present | Absent | Otherwise |
|---|---|---|---|
| B | lower limit of ρ_O above 1 | 95 % upper limit of ρ_O at or below 1.10 | not determined |
| P | lower limit of E_comp above 0 | 95 % upper limit of E_comp at or below 0.005 | not determined |

- Lower limit for "present": the 0.05/12 quantile per setting; the 0.025 quantile for the pool. Upper limit for
  "absent": the 0.95 quantile. Both met: **absent**, printed "non-zero, within the margin".
- Label: bias (B present, P absent); power (P present, B absent); both; neither (both absent); not determined.
- Settings labelled: the 12 vd and the 12 mcar settings. The **deciding settings** are the ten vd settings other
  than D at vd5 and vd10.
- **Primary decision**: ρ_O and E_comp pooled over the ten deciding settings (ratios of sums; each cell's
  replicates resampled once for all its masks), labelled at one-sided 0.025 for present and 95 % for absent. The
  matching mcar reading is the pool of the same ten (cell, rate) pairs under mcar, labelled by the same rules
  (R1-04).
- Attribution of each present component (pool and each vd setting): value-dependent if that component is absent under
  the matching mcar reading, not specific to value dependence if present there, otherwise not determined.
- F1: any protein in O and not in H (L1, S_a, S_b): stop. F2: G1 fails: stop, S_b not sent. F3: E2's first half fails
  in L1: stop. F4 (S_b only): the two-sided Bonferroni interval (quantiles 0.05/48 and 1 − 0.05/48) of m6_NA in any of
  the 24 settings excludes 0: no label is assigned, the owner is told.
- G1 (S_a, before any estimand): per vd setting, the means over r = 1 to 100 of |O| (S2 base `kept`), |H| (head
  `kept`), |Δ| (`added`), head `calls`, `power`, `fdp`, `added_share_called`, `added_fdp` (and the count of replicates
  for the last two) equal `kit-187/s2/cells.tsv` to a relative 1e-9.

| Outcome | Condition (the first that holds) |
|---|---|
| Stop | F1, F2, F3 (first half) or F4 |
| O5 | pooled label neither |
| O4 | pooled label not determined (no replicates are added under this declaration) |
| O1 | pooled label bias |
| O2 | pooled label power |
| O3 | pooled label both |

## Sentence (ticket 5.7), as `score-k197.R` fills it

The template of the outcome with the pooled estimates (vd pool). "More | not more" unchanged proteins in O2 and O3:
"more" when the 95 % interval of π0_Δ − π0_O lies above 0, else "not more". "Does | does not | also appears"
clauses: the matching mcar pool's state of B (O1) or P (O2); O3 prints both. The power and null-share parts are
printed in the power-first order with the other order beside. The R4 clause (attenuation of changed proteins in Δ and
O) and the deciding settings whose own label differs from the pooled one are added; the closing clause of 5.7 ends
every sentence. O5 and Stop: no sentence, no NEWS clause until the owner decides.

## Budget and attempts

Expected 10.67 core-hours, ceiling 25.33 (ticket 10; core-hours = job minutes × 4 / 60). G2 from L1's seconds per
replicate-mask. L1 once (repeated only after a script change, inside S′); S_a once; S_b once; scoring once per job.
An unfavourable or not-determined result is recorded, never re-run.
