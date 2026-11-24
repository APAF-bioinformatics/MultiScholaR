# KIT-190: declaration of the bounded-k measurement (VAL-182 (2), (3))

Written before any code or run of this measurement, to be audited before the lifecycle run is sent. Authority: methods
ledger VAL-182 (2026-11-25); audit `tickets/kit/audits/AUDIT-KIT-190-rule-d-20261125.md`, "The smallest further
measurement" and K190-M11. Nothing here is run until the auditor has checked it.

## 1. Question and decision rule (fixed now)

Does rule (d) with k bounded (D+B) give no worse error control than the first repair with k bounded (F+B), in every
cell? **If yes in every cell, D+B ships; if not in any cell, F+B ships** (VAL-182 (3)). The bound's value is decided by
the owner with this result (VAL-182 (1)); this measurement uses the cap 2 that VAL-182 names.

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
| N-r, N-s, E-r, E-s, E-p, E-g | as `effect.R` | balanced within groups (unchanged) | as `effect.R` | 50 | `190000000 + 10000 × cell + r`, cells 1 to 6 (unchanged) |
| **N-w** | none | random (below) | as generated | 50 | cell 7 |
| **E-w** | 10 % ±1 log2 | random (below) | as generated | 50 | cell 8 |
| E1 | T1b's 42 journeys (lipid seeds 172101 to 172120, with and without missing cells; seed 172101 metab excluded as metab); T11 concentration and F172 (lipid) | as KIT-172 | as KIT-172 | 1 each | KIT-172's |

- **The generator change, exactly.** One added argument to `e2Values` (the declared copy of `kit172Values`),
  `factor_pattern = c("balanced", "random")`. `"balanced"` (default) keeps the line
  `logs <- logs + outer(loadings, rep(c(-1, -1, -1, 1, 1, 1), 2L))` and the random stream unchanged, so cells 1 to 6
  reproduce the committed rows exactly (checked in the lifecycle). `"random"` replaces the pattern by
  `w <- stats::rnorm(12)`, drawn immediately after `loadings` (so the stream differs only from that point, only in the new
  cells), then `w <- (w - mean(w)) / stats::sd(w) * sqrt(12 / 11)`: centred, with the balanced pattern's sample
  variance exactly (the pattern ±1 has mean 0 and sample variance 12/11). So the
  factor's strength (loadings N(0, 0.5²) on log2) is unchanged; only its alignment with the groups is left to chance.
  Its correlation with the group indicator over 12 samples has standard deviation about 1/sqrt(11) = 0.30; the realised
  correlation is recorded per journey and reported (quartiles; and the endpoints split at |r| > 0.5, descriptively).
- **Precision at R = 50.** Paired binary differences (null cells): with discordance near 0.2 to 0.4, the exact
  conditional interval's half-width on the difference of rates is about 0.10 to 0.14. Paired FDP differences (effect
  cells): per-journey FDP SD about 0.2 to 0.3 (from `ed/`), paired-difference SD about 0.2, so the standard error of the
  mean difference is about 0.03 and the Bonferroni-adjusted one-sided bound sits about 0.08 above the mean.

## 4. Endpoints

| Cells | Primary (decides) | Secondary (reported, decides nothing) |
|---|---|---|
| N-r, N-s, N-w | **family-wise error**: the data set has any call at q < 0.05 (yes/no per journey) | calls per data set; k; mean \|log2 FC − 0\| |
| E-r, E-s, E-p, E-g, E-w, E1-T1b | **false discovery proportion** per journey (0 with no calls) | power over all planted features; k; mean \|log2 FC − truth\| over tested features; calls |
| E1-T11 (two fixtures) | none (n = 1 each): false calls, k reported | |

Also recorded per journey (instrumentation, VAL-182; turns the audit's M05 and M06 into shown or not shown): the k
diagnostic's two curves at the chosen percentage (`All` and `Control` canonical correlations for K = 1 to the maximum),
the separation score, the counts of features with residual dimensions (observed runs minus observed replicate sets)
below, at and above k, and the number of features RUV-III removed (all NA after correction).

## 5. "No worse": the exact paired comparison

Per cell c, pairs are the same data set under D+B and F+B.

- **Null cells (FWER).** b = data sets with a call under D+B only, c' = under F+B only, n = pairs. Statistic: the paired
  difference of rates Δ = (b − c') / n. Interval: exact conditional (the discordant pairs' split b ~ Binomial(b + c',
  θ), Clopper-Pearson upper bound θ_U at one-sided level α/m), giving Δ_U = (b + c') (2 θ_U − 1) / n. **No worse if
  Δ_U < 0.10.**
- **Effect cells (FDP).** d_i = FDP(D+B)_i − FDP(F+B)_i. Statistic: the mean of d_i. Interval: one-sided upper bound by
  the exact sign-flip permutation distribution of the mean (all 2^n sign patterns are not enumerable at n = 50, so
  100,000 random sign patterns with a fixed seed, 19001; Monte Carlo error on the bound reported). **No worse if the
  upper bound < 0.05.**
- **Multiplicity.** m = 9 decisive cells (3 null, 6 effect including E1-T1b); one-sided α = 0.05 / 9 per cell
  (Bonferroni), so the decision holds with family-wise confidence 95 %.
- **Direction.** Only "D+B worse" is tested (one-sided); D+B better is reported but changes nothing (the rule is
  "no worse").
- **Margins and their justification.** FDP margin 0.05: the nominal false discovery rate the route targets is 0.05, so
  a mean excess smaller than the target itself is taken as no material difference; the measured spread (section 3) lets
  a true difference of 0 pass with high probability. FWER margin 0.10: 5 data sets in 50; a smaller margin cannot be
  resolved at R = 50 (exact half-width 0.10 to 0.14), and VAL-182's own figures differ by 22 to 26 data sets in 50
  between the unbounded rules, so 0.10 separates material differences from noise at this size. Both margins are
  proposals for the auditor and the owner; they are fixed before the run and not changed after.
- Power is not part of "no worse" (VAL-182 (3) speaks of the rule being no worse; the error endpoints are the owner's
  concern in X168). Power and k differences are reported with paired 95 % intervals (sign-flip), unadjusted, labelled
  descriptive.

## 6. What follows

| Outcome | Follows |
|---|---|
| D+B no worse in all 9 cells | rule (d) with the bound ships (VAL-182 (3)); the bound's value goes to the owner with the k and FDP-by-k rows |
| D+B worse in any cell | the first repair with the bound ships (row-order dependence stays, stated in NEWS) |
| Mixed: worse only in N-w or E-w | as "worse in any cell": F+B ships. Reported to KIT-188's default ruling with the B+B reference (VAL-182's carried risk) |
| Both bounded rules far worse than B+B in N-w or E-w (descriptive: difference of rates or of mean FDP above 0.10) | reported to the owner and KIT-188 as the measured form of K190-M11; changes nothing here |
| A journey dropped or refused in one arm only | the pair is excluded, counted and reported; if more than 5 pairs in a cell are excluded, that cell's verdict is reported with the count and the owner is asked before shipping |

## 7. Lifecycle run, stop rules, budget

- **Script.** `docs/qa/workflow-execution/kit/kit-190/effect.R` extended (a `bound` mode; the generator argument; the
  arms list; the instrumentation) and committed, its SHA-256 recorded in `runs.tsv` and in `provenance.json`, before
  the lifecycle; the counted run uses that same hash (any change after the lifecycle is an instrument repair under
  VAL-54 and is recorded).
- **Lifecycle (R = 2 per cell, all arms, E1 two T11 fixtures).** Must show: (i) cells 1 to 6 under `"balanced"`
  reproduce the committed `e/`/`ed/` values of the first two replicates exactly (k, percentage, controls); (ii) every
  bounded journey ran with k = min(k_unbounded, 2) and the unbounded percentage, read back from its record; (iii) the
  bounded controls equal each rule's oracle at that percentage (F2's oracle, rule-specific); (iv) the instrumentation
  columns are filled; (v) a fault run (one bounded journey given k + 1) stops with a classed error after writing.
  Ceiling 30 job-minutes. A failed lifecycle is an instrument defect (VAL-54, at most three further runs).
- **Counted run.** Two jobs split by cell (old cells and E1; N-w, E-w with the metab check), each two trees at once,
  2 workers per tree. Expected 9 to 10 core-hours (audit's estimate: about 692 bounded journeys at about 17 s plus 100
  data sets × 6 arms in the new cells); **ceiling 20 core-hours** in all, lifecycle included.
- **Stop rules.** A falsifier ((ii) or (iii) above failing in any counted journey) stops the run with a classed error
  after writing what it has; reported, not re-run except under VAL-54 when the run's own log shows a mechanics failure
  with no result produced.
- **An unfavourable valid result is reported and not re-run**, whichever arm it disfavours.
- **Blindness.** The script and this declaration are committed and hashed before the lifecycle; the counted run binds
  the commit the lifecycle passed on. No figure of the counted run is looked at before all cells are written.
