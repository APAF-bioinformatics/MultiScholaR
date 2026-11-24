# T-COV rank re-run (KIT-159 r1.13, VAL-69): declaration, written before the request

**Purpose.** Re-run the four rank cells of T-COV (`pauc_diff`, class `rank`) with the re-pinned
generator of 10.2 (r1.13): C4-144 (form replicate, g = 4), C4-146 (replicate, g = 6), C4e-148
(block, g = 5), C4e-150 (block, g = 14). No other cell is re-run; every other row of
`coverage-table.tsv` keeps the classification it received in T2.

**Generator and arms (r1.13).** 4,840 features, 40 changed; the group effects scaled per form and g
by `valRankEffect` (grid 0.05 to 1 by 0.05; 50 seeded calibration simulations per grid value; the
value whose plain arm's mean normalised pAUC at FPR 0.01 lies in [0.4, 0.8] and is nearest 0.6;
no value in the band refuses the cell). Plain arm: moderated limma p; second arm: the same fit's
unmoderated t p (a different ranking, not a monotone transform). The chosen effect is recorded in
the table (`effect`).

**Rules applied (4.4, r1.13).** Approximate run-level band: coverage below level − 0.02 fails;
`conservative` when the simulated misses are below a quarter of the nominal miss rate. A cell is
`degenerate` (not offered) when more than 1 % of its intervals have zero width (`zero_width`
printed) or its target lies within one Monte Carlo SE of a bound of [−1, 1].

**Sizes and binding.** 2,000 simulations and 20,000 target simulations per cell; the ranges of T2
(`t2r-plan.tsv`, request R: 1, 2, 5, 5 ranges); 4 concurrent chunks; seeds as T2 (by cell and
index only). Bound to the harness commit named in the request (the implementation commit after the
T-RANK tests pass); the declaration is committed before the request.

**Ceiling.** One `run-tests` scripts request, script limit 280 minutes (`timeout_minutes`),
at most about 20 core-hours, inside T2's remaining ceiling (54 − 18.59 = 35.41). Expected wall:
the measured T2 process time of these cells (23,994 s) / f_sim 2.19 + calibration and set-up,
about 3.5 hours.

**What the result decides.** For each rank regime (form x g x level): offered (pass or
conservative), not offered (fail), or not offered (degenerate). `USAGE.md`'s table and the H2 row
are updated from this run only. It is evidence: one correction pass under the ticket's rule; an
unfavourable valid result is recorded, not re-run. KIT-162 is told the result before its
declaration (VAL-69).
