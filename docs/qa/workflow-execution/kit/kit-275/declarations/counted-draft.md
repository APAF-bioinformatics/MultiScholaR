# KIT-275 counted tiers: DRAFT declaration (`declaration.md` is committed alone after L1 and G2, before any counted job)

Not sent, not final. The text of sections 4 to 7 of ticket r1.2 as accepted is the declaration; this draft lists what the final file carries beyond it.

- Bound commit: `<40 digits>`; tier V (if it binds later): `<40 digits>` with its own G0 and L1-V.
- SHA-256 of `sim-k275.R`, `score-k275.R`, `fixture-k275.R`, `selftest/make-selftest.py`, `selftest/check-selftest.py`, `selftest/run-selftest.R`, `e0-diff.py`, this declaration's draft: `<hashes>`.
- Arms, versions, tiers, designs, R (200 for tier D; 100 for K, F, S, V, M, T, E; 20 for ID), seeds (generator 275M, planting 276M, bootstrap 277M with t by tier), endpoints and labels, margins (M_abs 0.05, M_act 0.50, M_rho 0.95, M_match 0.03, M_any 0.10, M_alpha 0.0125), levels, the outcome table, the falsifiers F1 to F4: as sections 4 to 7.
- Jobs (re-planned by G2 from L1's seconds; the ticket's figures are 38 jobs, every plan at most 85.31 minutes): `sim-k275.R <tiers> <design> <first> <last> <out>` with `K275_JOB`, `K275_DECLARATION`, `OMP_NUM_THREADS=1`, `timeout_minutes` as section 3 of `lifecycle.md`; no `K275_SMOKE`, no `K275_LIFECYCLE`.
  D: `D D3|D4|D6 1 50 | 51 100 | 101 150 | 151 200` (12 jobs, 115 min); K: `K ... 1 50 | 51 100` (6, 128); F and S: `F,S ... 1 50 | 51 100` (6, 77); M: `M D6 1 100` (110); ID: `ID D4 1 20` (17); T: `T ... 1 100` (3, 96); E: `E ... 1 100` (3, 92); V (waits): `V ... 1 50 | 51 100` (6, 98) with `K275_TIER_V_BOUND`.
- Order: tier D first; its PC1 and NC1 are read before the other tiers are scored; F1 or F2 stops the rest. Scoring once (`score-k275.R <rows> <out>` with `K275_G1=<G1 table>`), 40 job-minutes ceiling.
- Running-ceiling rule: before every request `runs.tsv` gets used + this ceiling + the ceilings still planned, at most 350.73 core-hours.
