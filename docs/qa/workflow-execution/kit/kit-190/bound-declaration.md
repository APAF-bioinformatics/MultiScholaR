# KIT-190: declaration of the bounded-k measurement (VAL-182 (2), (3)), revision 2

Revision 2 (2026-11-25) applies the declaration audit (`tickets/kit/audits/AUDIT-KIT-190-rule-d-20261125.md`,
"Declaration audit") B01 to B03 and B05 to B10 with its exact texts. **B04 is OPEN**: the margins, R per cell, whether
E1-T1b decides, and the ceiling are the owner's; every such value below is marked **[B04 OPEN]**, and the two options
put to the owner are worked out in the appendix. Revision 1 is commit fe1eeb49 (SHA-256 31f453e7ca78...cd97cd).

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
| N-r, N-s, E-r, E-s, E-p, E-g | as `effect.R` | balanced within groups (unchanged) | as `effect.R` | **[B04 OPEN]** (50 in the null cells; 50 or 100 in the effect cells) | `190000000 + 10000 × cell + r`, cells 1 to 6 (unchanged) |
| **N-w** | none | random (below) | as generated | **[B04 OPEN]** (150 in both options) | cell 7 |
| **E-w** | 10 % ±1 log2 | random (below) | as generated | **[B04 OPEN]** (150 in both options) | cell 8 |
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
- **R and power [B04 OPEN].** R per cell is the owner's choice between the appendix's options A and B (for replicates
  beyond 50 in the balanced cells, the unbounded F and D journeys are run, not reused). **Power, stated before the run.**
  With arms truly equal, the probability that D+B passes every decisive cell is **[B04 OPEN: the chosen option's figure
  from the appendix]**; a true FDP excess of 0.05 passes with probability at most α per cell.

## 4. Endpoints

| Cells | Primary (decides) | Secondary (reported, decides nothing) |
|---|---|---|
| N-r, N-s, N-w | **family-wise error**: the data set has any call at q < 0.05 (yes/no per journey) | calls per data set; k; mean \|log2 FC − 0\| |
| E-r, E-s, E-p, E-g, E-w; E1-T1b **[B04 OPEN: decisive or reported; reported in both options]** | **false discovery proportion** per journey (0 with no calls) | power over all planted features; k; mean \|log2 FC − truth\| over tested features; calls |
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
  disagree, that is reported to the owner, and the Tango verdict stands as declared. (Margin 0.10: **[B04 OPEN]**, 0.10
  in both options.)
- **Effect cells (FDP).** d_i = FDP(D+B)_i − FDP(F+B)_i. Test of H0: mean(d) ≥ 0.05 against H1: mean(d) < 0.05 by a
  one-sided sign-flip test of the studentised mean of e_i = d_i − 0.05: statistic t = mean(e)/(sd(e)/√n), reference
  distribution of t over 100,000 random sign vectors (seed 19001) plus the observed one. p = the share of sign vectors
  with t* ≤ t. No worse if p ≤ α. The test is asymptotically valid without symmetry and exact under symmetry; the Monte
  Carlo error of p is reported. The paired t bound is reported beside it. **E1-T1b:** the unit is the seed (21 seeds); e
  is the mean of a seed's two journeys. (Margin 0.05: **[B04 OPEN]**, 0.05 in option A, 0.10 in option B; the text's
  "0.05" is replaced by the chosen margin throughout this bullet.)
- **Multiplicity.** The decision is an intersection–union test: D+B ships only if every decisive cell rejects its H0. Each
  cell is tested at one-sided α = 0.05; the probability of shipping D+B when it is worse by at least the margin in any
  cell is then at most 0.05 (Berger 1982, Technometrics 24:295). Per-cell confidence bounds are reported at the same
  level, labelled non-simultaneous.
- **Direction.** Only "D+B worse" is tested (one-sided); D+B better is reported but changes nothing (the rule is
  "no worse").
- **Margins [B04 OPEN].** The owner's choice (appendix): FDP margin 0.05 (twice the precision cost) or 0.10 (twice the
  nominal FDR, well under the earlier differences of +0.23 to +0.39); null-cell FWER margin 0.10 (about 5 data sets in
  50, against earlier differences of 0.40 to 0.52). Fixed before the run and not changed after.
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
  lifecycle's timing allows. Further jobs split by cell if a job would pass 300 minutes. Expected core-hours and
  **ceiling: [B04 OPEN]** (appendix: option A about 26 expected, ceiling 35; option B about 20 expected, ceiling 30;
  lifecycle included).
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
