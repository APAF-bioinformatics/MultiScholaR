# KIT-190: declaration of the bounded-k measurement (VAL-182 (2), (3)), revision 6

**Revision 6 (2026-11-25): the scripts, committed for the auditor's script check before the lifecycle.** No design
change. Files and SHA-256 at this commit:
- `effect.R` (the `bound` mode: phases `lc`, `lc-stage1`, `stage1`, `stage2`; the generator argument; arms; checks;
  instrumentation): 6cd59c4c6120dd12a930ca08cae9b332a72b9d1c7bc6b35422d78210d4d955ff
- `reestimate.R` (`estimate`, `plan`, `selftest`; section 8): b171924cae9493d3a6b9ff080906ac27d5b3a1f7488579269b5026fd53053226
- `decide.R` (sections 5 and 6: the tests, the intersection-union decision, the reported items; `selftest` for item
  (vi)); a third script, named here for the first time: d908af20d51cc58d03f16e29fcf2c501b973837bb921e5cac2985bb7480cf01c
- `size2stage.py` (unchanged from revision 5): 07cf59df0bfc528062cd49e9644fcf2bb184a21cb46cb1ab105c19178ec18d27

Points for the script check, where the scripts read the text more narrowly than it could be read:
1. **E1-T1b, lipid only:** 20 seeds (172101 to 172120), 40 journeys; "21 seeds" in section 5 counted the metab seed,
   which this lipid-only measurement excludes. T1b is reported, not deciding (VAL-185).
2. **Check (i)** compares k, percentage, number of controls and calls with the committed rows; the control identities
   themselves were not committed in `e/` and `ed/`, so they cannot be compared.
3. **E1 runs in stage 2 only** (it is reported, has a fixed size, and plays no part in the re-estimation).
4. **Secondary differences** (k, power, fold-change error) are reported with paired t intervals, not sign-flip
   intervals (descriptive only).
5. **The lifecycle's runs:** `effect.R bound lc` (R = 2, all cells and E1's two T11 fixtures; checks (i) to (iv),
   (vii)); `effect.R bound lc` with `K190_FAULT=1`, N-r only (check (v)); `decide.R selftest` (vi);
   `reestimate.R selftest` (viii); `effect.R bound lc-stage1` twice, the second with `K190_PLANT=1` (ix), each with its
   own `out`. In the lifecycle's stage-1 runs the unbounded old-cell journeys are reused from the committed rows, as in
   stage 1.


**Revision 5 (2026-11-25):** the audit's "Declaration audit, revision 4" (`tickets/kit/audits/AUDIT-KIT-190-rule-d-20261125.md`)
applied with its exact texts: R4-01 (skewness check with Johnson's statistic; owner's word pending, section 5), R4-02
(stage 1's console and provenance; item (ix)), R4-03 (R_c computed by `reestimate.R`; stage 2 always runs), R4-04
(`size2stage.py` reproduces the table; the appendix's ceiling). Revision 4 is commit 2c78be95 (SHA-256 516447008754...93e0).

**Revision 4 (2026-11-25):** B04 decided by the owner (methods ledger VAL-185): option C, with one blinded re-estimation
of the number of data sets (section 8); ceiling 100 core-hours; 8 deciding cells, E1-T1b reported; margins 0.10 (null
cells) and 0.05 (effect cells). Every former "[B04 OPEN]" value is filled from option C below. Revision 3 is commit
d7c5ef98 (option C added; SHA-256 cd1acd5b...0c55); revision 2 is ae8e85af.

Revision 2 (2026-11-25) applied the declaration audit (`tickets/kit/audits/AUDIT-KIT-190-rule-d-20261125.md`,
"Declaration audit") B01 to B03 and B05 to B10 with its exact texts. **B04 is OPEN**: the margins, R per cell, whether
E1-T1b decides, and the ceiling were the owner's (decided in revision 4); the options put to the owner are in the
appendix. Revision 1 is commit fe1eeb49 (SHA-256 31f453e7ca78...cd97cd).

Written before any code or run of this measurement, to be audited before the lifecycle run is sent. Authority: methods
ledger VAL-182 (2026-11-25); audit `tickets/kit/audits/AUDIT-KIT-190-rule-d-20261125.md`, "The smallest further
measurement" and K190-M11. Nothing here is run until the auditor has checked it.

## 1. Question and decision rule (fixed now)

Does rule (d) with k bounded (D+B) give no worse error control than the first repair with k bounded (F+B), in every
cell? **If yes in every cell, D+B ships; if not in any cell, F+B ships** (VAL-182 (3)). The bound's value is decided by
the owner with this result (VAL-182 (1)); this measurement uses the cap 2 that VAL-182 names. The result applies to 12
samples (cap 2) and to a cap applied after the sweep. A bound that truncates the curve before the elbow, or the adaptive
cap at larger designs (3 to 5), is not measured here; the ticket that implements the bound reports, on these seeds, how
many journeys' percentage or k differ from this emulation.

## 2. Arms (trees fixed by commit)

| Arm | Control rule | k | Tree | Source of the figures |
|---|---|---|---|---|
| B | legacy (row order) | Automatic, unbounded | d4e2eb2c | reused from `e/`, `ed/` (identical there) in the six old cells and E1; run in N-w, E-w |
| F | first repair (q, table order) | Automatic, unbounded | d6d67efe | reused from `e/`; run in N-w, E-w |
| D | rule (d) | Automatic, unbounded | a1c12d28 | reused from `ed/`; run in N-w, E-w |
| **F+B** | first repair | min(k_F, 2), Manual at F's chosen percentage | d6d67efe | run in every cell |
| **D+B** | rule (d) | min(k_D, 2), Manual at D's chosen percentage | a1c12d28 | run in every cell |
| B+B (reference) | legacy | min(k_B, 2), Manual at B's chosen percentage | d4e2eb2c | run in N-w and E-w only: the group-blind reference the audit's model (K190-M11) points to |

- **Bounded arms through the existing Manual mode** (`kit172LegacyRuv("manual")` with `manual_percentage` = the unbounded
  arm's chosen percentage for that journey and `manual_k` = min(k, 2)): no package code. Manual uses the same control
  rule as Automatic in each tree (lipid and metab, D8 (b)). This emulates a cap applied after the sweep; a cap that
  truncates the curve before the elbow (and could move the chosen percentage) is measured in the bound's own runs.
- In the six old cells and E1, the unbounded k and percentage come from the committed rows (`e/*/journeys.tsv` for F,
  `ed/*/journeys.tsv` for D and B), matched by cell, class, replicate and data set. A journey whose unbounded arm
  dropped the assay (no percentage) has no bounded counterpart; it is counted and reported, and its pair is excluded
  from the paired comparison (with the count stated).
- In N-w and E-w each journey runs B, F and D unbounded first, then B+B, F+B and D+B on the same data.
- **Class: lipid only.** Lipid and metab gave identical rows in all 1,292 journeys of `e/` and `ed/` (audit K190-M01).
  A metab check of 5 replicates per cell for F+B and D+B is added (inside the ceiling) and must equal the lipid rows;
  if it does not, that is reported and the lipid rows still decide.

## 3. Cells, R and seeds

| Cell | Effects | Factor's sample pattern | Row order | R | Seed |
|---|---|---|---|---|---|
| N-r, N-s, E-r, E-s, E-p, E-g | as `effect.R` | balanced within groups (unchanged) | as `effect.R` | at least 100 (N-r, N-s) and 200 (E-r, E-s, E-p, E-g); final R by section 8 | `190000000 + 10000 × cell + r`, cells 1 to 6 (unchanged) |
| **N-w** | none | random (below) | as generated | at least 250; final R by section 8 | cell 7 |
| **E-w** | 10 % ±1 log2 | random (below) | as generated | at least 200; final R by section 8 | cell 8 |
| E1 | T1b's 42 journeys (lipid seeds 172101 to 172120, with and without missing cells; seed 172101 metab excluded as metab); T11 concentration and F172 (lipid) | as KIT-172 | as KIT-172 | 1 each | KIT-172's |

- **The generator change, exactly.** One added argument to `e2Values` (the declared copy of `kit172Values`),
  `factor_pattern = c("balanced", "random")`. `"balanced"` (default) keeps the line
  `logs <- logs + outer(loadings, rep(c(-1, -1, -1, 1, 1, 1), 2L))` and the random stream unchanged, so cells 1 to 6
  reproduce the committed rows exactly (checked in the lifecycle). `"random"` replaces the pattern by
  `w <- stats::rnorm(12)`, drawn immediately after `loadings` (so the stream differs only from that point, only in the new
  cells), then `w <- (w - mean(w)) / stats::sd(w) * sqrt(12 / 11)`: centred, with the balanced pattern's sample
  variance exactly (the pattern ±1 has mean 0 and sample variance 12/11). So the
  factor's strength (loadings N(0, 0.5²) on log2) is unchanged; only its alignment with the groups is left to chance.
  Its correlation with the group indicator over 12 samples has standard deviation 0.30 (r² ~ Beta(1/2, 5)); mean |r|
  0.25, median 0.22, |r| > 0.5 in about 10 % of data sets; the realised
  correlation is recorded per journey and reported (quartiles; and the endpoints split at |r| > 0.5, descriptively).
- **R and power.** Option C (VAL-185) gives each cell's minimum R: N-r, N-s 100; E-r, E-s, E-p, E-g 200; N-w 250;
  E-w 200. The final R per cell is set once, by section 8's fixed formula, from the first 50 seeds of the capped arms;
  it is never below option C's. For replicates beyond 50 in the balanced cells, the unbounded F and D journeys are run,
  not reused. **Power, stated before the run.** At option C's R, with arms truly equal, the probability that D+B passes
  every decisive cell is 0.80 (paired SD 0.20, N-w discordance 0.25) to 0.86 (discordance 0.12), and 0.93 to 0.995 at
  SD 0.10; section 8 raises R to keep it at 0.80 or more if the observed spread is larger, within the ceiling. A true
  excess equal to the margin passes its cell with probability at most α (section 8 gives the two-stage procedure's
  simulated size).

## 4. Endpoints

| Cells | Primary (decides) | Secondary (reported, decides nothing) |
|---|---|---|
| N-r, N-s, N-w | **family-wise error**: the data set has any call at q < 0.05 (yes/no per journey) | calls per data set; k; mean \|log2 FC − 0\| |
| E-r, E-s, E-p, E-g, E-w (deciding); E1-T1b (reported, VAL-185) | **false discovery proportion** per journey (0 with no calls) | power over all planted features; k; mean \|log2 FC − truth\| over tested features; calls |
| E1-T11 (two fixtures) | none (n = 1 each): false calls, k reported | |

Also recorded per journey (instrumentation, VAL-182; turns the audit's M05 and M06 into shown or not shown): the k
diagnostic's two curves at the chosen percentage (`All` and `Control` canonical correlations for K = 1 to the maximum),
the separation score, the counts of features with residual dimensions (observed runs minus observed replicate sets)
below, at and above k, and the number of features RUV-III removed (all NA after correction). For journeys whose unbounded
arm is reused from `e/` or `ed/`, the script recomputes the k diagnostic's curves and separation score by one `ruvCancor`
call at the recorded percentage on the regenerated RUV input, and the residual-dimension counts at both the unbounded k
and the bounded k; the features removed by RUV-III are counted for every journey the script runs.

## 5. "No worse": the exact paired comparison

Per cell c, pairs are the same data set under D+B and F+B.

- **Null cells (FWER).** b = data sets with a call under D+B only, c′ = under F+B only, n = pairs; Δ = p10 − p01. Test
  of H0: Δ ≥ 0.10 against H1: Δ < 0.10 by Tango's score statistic (Tango 1998, Stat Med 17:891), Z = (b − c′ − n·0.10) /
  √(n(2q̂ + 0.10·0.90)), with q̂ the restricted maximum-likelihood estimate of p01 under Δ = 0.10. No worse if Z < −z_{1−α},
  α as in 'Multiplicity'. This test is asymptotic. Beside it is reported the Berger–Boos exact unconditional p-value
  (statistic b − c′, nuisance p10 + p01 restricted to its 1 − α/10 Clopper–Pearson interval, plus α/10); if the two
  disagree, that is reported to the owner, and the Tango verdict stands as declared. (Margin 0.10, VAL-185.)
- **Effect cells (FDP).** d_i = FDP(D+B)_i − FDP(F+B)_i. Test of H0: mean(d) ≥ 0.05 against H1: mean(d) < 0.05 by a
  one-sided sign-flip test of the studentised mean of e_i = d_i − 0.05: statistic t = mean(e)/(sd(e)/√n), reference
  distribution of t over 100,000 random sign vectors (seed 19001) plus the observed one. p = the share of sign vectors
  with t* ≤ t. No worse if p ≤ α. The test is asymptotically valid without symmetry and exact under symmetry; the Monte
  Carlo error of p is reported. The paired t bound is reported beside it. **E1-T1b:** the unit is the seed (21 seeds); e
  is the mean of a seed's two journeys. (Margin 0.05, VAL-185.) E1-T1b is
  reported with this test, not deciding (VAL-185).
  **Skewness check, fixed now (owner's word pending: R4-01; the run is not sent until it is given).** In each effect
  cell, on all R_c pairs, the sample skewness g of e is computed (third central moment over the cubed SD, both with
  divisor n). If |g| ≤ 1, the cell is decided as above. If |g| > 1, it is decided by Johnson's skewness-corrected
  statistic (Johnson 1978, JASA 73:536), t_J = [ē + μ̂₃/(6s²n) + μ̂₃ē²/(3s⁴)] / (s/√n), with μ̂₃ the third central moment
  (divisor n) and s the SD (divisor n − 1). The cell is no worse if t_J < −z_{0.95}. The sign-flip result and g are
  reported beside it. Simulated size with the re-estimation (audit, revision 4): normal 0.054; exponential (g ≈ 2)
  0.053–0.058; zero-inflated (80 % zeros, g ≈ 4) 0.053, where the studentised test alone gives 0.099. g is computed in
  the final analysis, on all R_c seeds, from the paired differences only; it is not computed in stage 1 and does not
  feed the re-estimation (`reestimate.R` computes only the SD and the discordant share).
- **Multiplicity.** The decision is an intersection–union test: D+B ships only if every decisive cell rejects its H0. Each
  cell is tested at one-sided α = 0.05; the probability of shipping D+B when it is worse by at least the margin in any
  cell is then at most 0.05 (Berger 1982, Technometrics 24:295). Per-cell confidence bounds are reported at the same
  level, labelled non-simultaneous.
- **Direction.** Only "D+B worse" is tested (one-sided); D+B better is reported but changes nothing (the rule is
  "no worse").
- **Margins (VAL-185).** FDP 0.05 (the route's nominal FDR); null-cell family-wise rate 0.10 (about 5 data sets in 50,
  against earlier differences of 0.40 to 0.52). Fixed before the run and not changed after.
- Power is not part of "no worse" (VAL-182 (3) speaks of the rule being no worse; the error endpoints are the owner's
  concern in X168). Power and k differences are reported with paired 95 % intervals (sign-flip), unadjusted, labelled
  descriptive.

## 6. What follows

| Outcome | Follows |
|---|---|
| D+B shown no worse in every decisive cell | rule (d) with the bound ships (VAL-182 (3)); the bound's value goes to the owner with the k and FDP-by-k rows |
| D+B not shown no worse in any decisive cell (its test does not reject at α, including an inconclusive cell) | the first repair with the bound ships (row-order dependence stays, stated in NEWS) |
| Mixed: not shown no worse only in N-w or E-w | as the row above: F+B ships. Reported to KIT-188's default ruling with the B+B reference (VAL-182's carried risk) |
| Both bounded rules far worse than B+B in N-w or E-w (descriptive: difference of rates or of mean FDP above 0.10) | reported to the owner and KIT-188 as the measured form of K190-M11; changes nothing here |
| A journey dropped or refused in one arm only | the pair is excluded, counted and reported; if more than 5 pairs in a cell are excluded, that cell's verdict is reported with the count and the owner is asked before shipping |

## 7. Lifecycle run, stop rules, budget

- **Script.** `docs/qa/workflow-execution/kit/kit-190/effect.R` extended (a `bound` mode; the generator argument; the
  arms list; the instrumentation) and committed, its SHA-256 recorded in `runs.tsv` and in `provenance.json`, before
  the lifecycle; the counted run uses that same hash (any change after the lifecycle is an instrument repair under
  VAL-54 and is recorded).
- **Before the lifecycle:** the extended script is checked against this declaration by the auditor (a short closure
  check) before the lifecycle is sent.
- **Lifecycle (R = 2 per cell, all arms, E1 two T11 fixtures).** Must show: (i) cells 1 to 6 under `"balanced"`
  reproduce the committed `e/`/`ed/` values of the first two replicates exactly (k, percentage, controls); (ii) every
  bounded journey ran with k = min(k_unbounded, 2) and the unbounded percentage, read back from its record; (iii) the
  bounded controls equal each rule's oracle at that percentage (F2's oracle, rule-specific); (iv) the instrumentation
  columns are filled; (v) a fault run (one bounded journey given k + 1) stops with a classed error after writing;
  (vi) the decision code run on two constructed inputs: identical arms (must pass) and a shift of +0.20 in every pair
  (must fail), in both the null and the effect test; (vii) the lifecycle's journey times projected to R as declared; if
  the projection exceeds the ceiling, the counted run is split by cell, never reduced in R. Ceiling 30 job-minutes. A failed lifecycle is an instrument defect (VAL-54, at most three further runs).
- **Counted run.** Job 1 (old cells, E1): trees d6d67efe and a1c12d28 at once, 2 workers each; job 2 (N-w, E-w, metab
  check): trees d4e2eb2c, d6d67efe and a1c12d28, unbounded arms first and bounded after, 1 worker per tree or as the
  lifecycle's timing allows. Further jobs split by cell if a job would pass 300 minutes. Stage 1 (seeds 1 to 50 of the
  capped arms in every cell; with the unbounded arms in N-w and E-w) runs first, then the re-estimation (section 8), then
  stage 2 (the remaining seeds). About 48 core-hours expected at option C's R; **ceiling 100 core-hours** in all,
  lifecycle, stage 1, re-estimation and stage 2 included (VAL-185).
- **Stop rules.** A falsifier ((ii) or (iii) above failing in any counted journey) stops the run with a classed error
  after writing what it has; reported, not re-run except under VAL-54 when the run's own log shows a mechanics failure
  with no result produced.
- **An unfavourable valid result is reported and not re-run**, whichever arm it disfavours.
- **Blindness.** The script and this declaration are committed and hashed before the lifecycle; the counted run binds
  the commit the lifecycle passed on. No figure of the counted run is looked at before all cells are written.


## Appendix (B04): the two options put to the owner

Common to both: intersection–union at one-sided α = 0.05 per cell (B03); E1-T1b and the two T11 fixtures reported, not
decisive (audit B04 (d)); decisive cells N-r, N-s, N-w (null) and E-r, E-s, E-p, E-g, E-w (effect): 8 cells; null-cell
FWER margin 0.10; B+B in N-w and E-w only.

| | Option A (auditor's preference) | Option B (cheaper) |
|---|---|---|
| FDP margin | 0.05 | 0.10 |
| R: N-r, N-s | 50, 50 | 50, 50 |
| R: E-r, E-s, E-p, E-g | 100 each (replicates 51 to 100 run unbounded F and D as well) | 50 each |
| R: N-w, E-w | 150, 150 | 150, 150 |
| E1 (lipid) | T1b 42 journeys + 2 T11 fixtures, reported | the same |
| Journey-seconds (estimate) | about 68,800 s | about 50,600 s |
| Expected core-hours (lifecycle included) | about 26 | about 20 |
| Ceiling | 35 | 30 |
| P(D+B passes all 8 cells), arms equal, paired SD of FDP difference 0.20; N-w discordance 0.12 / 0.25 | 0.32 / 0.27 | 0.75 / 0.61 |
| same, paired SD 0.10 | 0.84 / 0.69 | 0.84 / 0.69 |
| P(pass all), true differences as the first (unbounded) measurement saw (FDP +0.23 or more; null any-call +0.42 or more) | below 1e-9 | below 1e-9 |
| A true FDP excess equal to the margin | passes its cell with probability ≤ 0.05 | the same |

**How these were computed** (stdlib Python, `random` and `statistics.NormalDist`, no R; script in the implementer's
scratchpad `k190py/power.py`, not committed; indicative):
- Null cells: Tango's statistic with q̂ by bisection on the restricted score, simulated (20,000 data sets; 2,000 to 5,000
  for the "observed difference" rows) with p10 = p01 = discordance/2 for equal arms. N-r and N-s use discordance 0.02
  (bounded k at most 2: VAL-182 records no null data set with k of 2 or less with a call), giving 0.939 each at R = 50;
  N-w uses 0.12 or 0.25 (the audit's model, K190-M11), giving 0.952 or 0.785 at R = 150. These reproduce the audit's
  B01 figures (0.94 / 0.58 / 0.39 at R = 50).
- Effect cells: normal approximation to the studentised test, P = Φ((margin − Δ)·√n / SD − z_0.95), with the audit's
  paired-SD estimates 0.10 and 0.20.
- Joint probability: the product over the 8 cells (independent data sets).
- Journey-seconds: bounded Manual 17 s; unbounded B 15 s, F 25 s, D 32 s (`ed/*/provenance.json`, audit's cost
  section); core-hours = journey-seconds / 3,600 × 1.25 (set-up, compile, idle) plus about 2 for the lifecycle.

### Option C (added at the coordinator's request, 2026-11-25)

FDP margin 0.05 and null-cell margin 0.10 as in A; 8 deciding cells; T1b reported. R per cell: the smallest multiples
of 50 (searched over 50 to 400 per cell type, cheapest combination kept) with P(pass all | arms equal) ≥ 0.80 at paired
SD 0.20 and N-w discordance 0.25.

| | Option C |
|---|---|
| R: N-r, N-s | 100, 100 |
| R: E-r, E-s, E-p, E-g | 200 each |
| R: N-w, E-w | 250, 200 |
| Journey-seconds (estimate) | about 132,700 s |
| Expected core-hours (lifecycle included) | about 48 |
| Ceiling | 100 (VAL-185) |
| P(pass all), arms equal, SD 0.20, N-w discordance 0.25 / 0.12 | 0.80 / 0.86 |
| same, SD 0.10 | 0.93 / 0.995 |
| P(pass all), first measurement's differences | below 1e-9 |

Computed as for A and B (null cells simulated with 4,000 data sets per R; N-r and N-s at discordance 0.02 need R = 100
because 0.939² at R = 50 alone caps the joint probability near 0.88).

**Reuse.** The six balanced cells keep their generator and seed rule (`190000000 + 10000 × cell + r`), so replicates 1
to 50 reuse the committed unbounded rows of `e/` and `ed/` (B, F, D) and only their bounded arms are run; replicates 51
onwards are new seeds and run F and D unbounded plus F+B and D+B. The costs above count it so (bounded only for
r ≤ 50; unbounded F and D plus both bounded arms for r > 50; all six arms in N-w and E-w).

### Observed spread in the first measurement (committed rows only; stdlib Python on `e/*/journeys.tsv` and `ed/*/journeys.tsv`, lipid)

| Cell | SD of FDP(D) − FDP(F), unbounded, n = 50 | SD of FDP(F) − FDP(B), unbounded | pairs with both k ≤ 2 (SD) |
|---|---|---|---|
| E-r | 0.417 | 0.313 | 4 (0.035) |
| E-s | 0.354 | 0.240 | 1 |
| E-p | 0.406 | 0.514 | 6 (0.312) |
| E-g | 0.418 | 0.337 | 1 |

| Cell | discordance, any call, D vs F unbounded (D only / F only) | discordance, F vs B unbounded |
|---|---|---|
| N-r | 0.64 (28 / 4) | 0.04 |
| N-s | 0.58 (25 / 4) | 0.02 |

These are for the **unbounded** arms, whose k differ widely (D 5 to 6 on average); the bounded arms' spread is not
known. The pairs with both k ≤ 2 are too few to estimate it. Unbounded, the observed paired SD (0.35 to 0.42) is above
both assumptions (0.10, 0.20); the null discordance between two low-k arms (F against B: 0.02 to 0.04) is near the 0.02
assumed for N-r and N-s.

## 8. One blinded re-estimation of R (VAL-185)

**Stage 1.** Seeds 1 to 50 of the capped arms (F+B, D+B) in every cell, in N-w and E-w together with their unbounded
arms (B, F, D, B+B), lipid. The stage-1 job writes every journey row, side log and intermediate file to a directory
outside `out`, which is discarded with the runner. Its console prints only fixed lines: the number of journeys per cell
and arm, the counts of the mechanics falsifiers (the k and percentage read-back and the oracle), seconds, and the
`reestimate.tsv` table. `out` holds only `reestimate.tsv` and `provenance.json`. `provenance.json`'s fields are commits,
script hashes, seeds, journey counts, falsifier counts, and the median seconds per journey of each arm type; no field
derives from calls, FDP, power, fold changes or a per-arm rate.

**The re-estimation script** (`docs/qa/workflow-execution/kit/kit-190/reestimate.R`, committed and hashed with
`effect.R` before the lifecycle, checked by the auditor). It reads the stage-1 journey rows of the two capped arms on the
runner and writes one file, `reestimate.tsv`, with exactly four columns: `cell`, `endpoint`, `n_pairs`, `value`, where
`endpoint` is `fdp_difference_sd` (effect cells: the sample SD of d_i = FDP(D+B)_i − FDP(F+B)_i, from `stats::sd`, which
is invariant to a shift of every d_i) or `any_call_discordant_share` (null cells: (b + c′) / n, which is invariant to
exchanging b and c′). It computes no mean difference, test statistic or per-arm rate; it prints only this table.
- **Enforced** by the script itself: before writing it stops unless `names(out)` is exactly the four columns and
  `endpoint` takes only those two values; no other object is saved; the console shows only the table.
- **Checked** in the lifecycle (item (viii)): the script run on planted inputs must give byte-identical files for (a) an
  effect cell whose d_i are all shifted by +0.30 and the same cell unshifted, and (b) a null cell with (b, c′) = (10, 0),
  (0, 10) and (5, 5) on 50 pairs; and the auditor reads the script against this section before the lifecycle.

**The formula, fixed now.** α = 0.05 one-sided, z = z_0.95 = 1.6449, per-cell target t = 0.80^(1/8) = 0.97250 (so the
product over the 8 deciding cells is 0.80 when every cell meets it).
- Effect cell c with stage-1 SD ŝ_c: P_c(R) = Φ(0.05·√R / ŝ_c − z).
- Null cell c with stage-1 discordant share δ̂_c: P_c(R) = Φ(0.10·√R / √δ̂_c − z) (the normal approximation to Tango's
  test when the arms are equal: b − c′ has variance n·δ).
- R_c = the smallest multiple of 50 that is at least option C's R for that cell and gives P_c(R_c) ≥ t, at most 1,000
  per cell. If ŝ_c or δ̂_c is 0, R_c = option C's R.
- **Ceiling.** Projected cost = spent so far (lifecycle and stage 1, from the runs) + stage 2's journeys (seeds 1 to R_c,
  see "Stage 2") × the stage-1 median seconds per journey of each arm type × 1.25. If it exceeds 100 core-hours, t is
  lowered for all cells together in steps of 0.0005 until it fits; R_c never goes below option C's R. The achieved joint
  probability Π P_c(R_c) under equality is reported with the final R. (Option C's own R costs about 48; it cannot bind.)
- The final R per cell is written to `runs.tsv` and `outcome.md` before stage 2 is sent, and is not changed after.
- R_c is computed by `reestimate.R`, from `reestimate.tsv`, the core-hours recorded in `runs.tsv` for the lifecycle and
  stage 1, and the stage-1 median seconds in `provenance.json`. Its output (R_c per cell, the target t used, the achieved
  Π P_c(R_c)) is committed before stage 2 is sent. If P_c(1,000) < t, then R_c = 1,000. Stage 2 runs at these R_c
  whatever stage 1 showed: no early stop, and no change of cells, arms, margins or tests.

**Stage 2.** Seeds 1 to R_c in every cell. Seeds 1 to 50 are run again (their rows were not returned): they are
deterministic, and stage 2 checks that re-running `reestimate.R` on its own seeds 1 to 50 reproduces the stage-1
`reestimate.tsv` byte for byte (a falsifier if not). **The final tests use all R_c seeds.** The stage-1 re-run costs
about 8 core-hours, inside the ceiling.

**Why the error rate is kept, and its limits.**
- R is set from a nuisance parameter only (the spread of the paired differences; the share of discordant pairs), never
  from the observed difference. For normally distributed paired differences the sample SD is independent of the sample
  mean, so a rule that uses only the SD cannot steer R toward a favourable mean; the remaining effect is the use of the
  stage-1 SD in the final statistic, which the simulation below finds negligible.
- The FDP differences are bounded and may be skewed, so the independence is only approximate; and the studentised test
  itself is anticonservative in the lower tail when the differences are right-skewed, at any fixed R (B02 accepted it as
  asymptotically valid).
- **Paired binary endpoint: the known limit.** The discordant share δ = p10 + p01 and the difference Δ = p10 − p01 are not
  separable: δ ≥ |Δ|, and at a given p01 a larger Δ means a larger δ. So δ̂ carries information on Δ, and the rule's R
  depends weakly on the difference. Its effect on the size is shown only by the simulation below; I know of no general
  result that bounds it, and none is claimed.

**Simulated size of the whole two-stage procedure at the margin** (`size2stage.py`, committed beside this file; stdlib
Python; run `python3 -I size2stage.py`, which prints these rows with their seeds; floors are option C's R (null 100 for
N-r and N-s, 250 for N-w; effect 200); the formula above, per-cell maximum 1,000; the final test is Tango's statistic in
null cells and, in effect cells, the studentised t with the R4-01 skewness check (Johnson's statistic when |g| > 1), the t
reference standing in for the sign-flip test; one cell at a time):

| Endpoint | Setting | Simulated cells | Size (nominal 0.05) | Mean final R |
|---|---|---|---|---|
| null (Tango, margin 0.10), floor 100 | discordance 0.10 / 0.25 | 6,000 each | 0.059 / 0.046 | 163 / 350 |
| null, floor 250 | discordance 0.12 / 0.25 / 0.40 | 6,000 each | 0.049 / 0.049 / 0.051 | 254 / 353 / 539 |
| effect (margin 0.05), normal | SD 0.20 / 0.40 | 20,000 each | 0.051 / 0.050 | 232 / 821 |
| effect, right-skewed (centred exponential, skewness 2) | SD 0.20 / 0.40 | 20,000 each | 0.057 / 0.055 | 245 / 769 |
| effect, zero-inflated (80 % zeros, else exponential; skewness about 4) | SD 0.20 / 0.40 | 20,000 each | 0.052 / 0.057 | 267 / 672 |
| audit re-simulation (revision 4), 5,000 cells each | null, discordance 0.10 to 0.40 | | 0.047 to 0.055 (floors 100, 250) | |
| | effect with the R4-01 rule: normal / exponential / zero-inflated | | 0.054 / 0.053 to 0.058 / 0.053 (studentised alone: 0.099 zero-inflated) | |

Monte Carlo SE about 0.0015 (20,000) and 0.003 (6,000). At the margin the discordance cannot be below 0.10 (δ ≥ |Δ|), so
lower discordances are not size settings. Read: with the R4-01 rule in place the whole two-stage procedure keeps its size
within 0.046 to 0.059 in every setting simulated; the largest, 0.059, is the asymptotic Tango test at discordance 0.10
with the floor of 100.

**Lifecycle additions.** Item (viii) above (planted inputs to `reestimate.R`); and (ix): on the lifecycle's own stage-1 run: the artifact holds the infrastructure's files
(`scripts.tsv`, `library.tsv`, `cpu.txt`, `blas.txt`, `r-version.txt`, `console.log`) and, under `out/`, exactly
`reestimate.tsv` and `provenance.json`; `console.log` contains only the fixed lines above; and a second stage-1
lifecycle run, with the D+B FDP rows of one cell shifted by +0.30 and the any-call flags of one null cell swapped
between arms (a planted switch of the bound mode, lifecycle only), gives an artifact identical to the first after
removing the seconds fields. The lifecycle
stays at R = 2 per cell, after the auditor's script check.
