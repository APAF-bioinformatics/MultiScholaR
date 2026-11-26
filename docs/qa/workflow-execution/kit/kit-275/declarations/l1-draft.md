# KIT-275 L1 (the blind lifecycle): DRAFT declaration, not sent

To be completed and sent only after the script audit and the coordinator's word. Placeholders in angle brackets.

- Purpose: prove the instrument through completion and failure (VAL-11): every arm of every tier at 3 data sets per design (r = 9101, 9102, 9103) and E9's second run of the first data set of each (tier, design); checks E0 to E12; the scorer blind on L1's own rows. Writes no failure rate, correlation, standard-error ratio, any call or rate.
- Binding: the bound commit `<40 digits of the first staging commit holding KIT-255 with I01; today ec27f69e0>`; scripts' commit `<kit-275 tip>` merged with that commit if it moved (G0 re-read).
- Request (one): `base_commit` none; `extra_packages` `["sva"]`; `fetch_refs` `["kit-188-ruv-default-calibration"]`; scripts:
  1. `selftest/run-selftest.R` (E0, E4), limit 15;
  2. `fixture-k275.R D3,D4,D6 1 20 k275-g1` is G1 (its own request, before L1);
  3. `sim-k275.R all D3,D4,D6 9101 9103 k275-l1` with `K275_LIFECYCLE=1 K275_FAULT=1 K275_JOB=l1 K275_DECLARATION=<declaration path> OMP_NUM_THREADS=1`, `timeout_minutes` 107 (130 if tier V is bound with the others; from section 3 of `lifecycle.md`, re-derived from the smoke's timing only as information).
- Ceiling: 117 job-minutes = 7.80 core-hours (140 = 9.33 with V); ledger row written before sending (`runs.tsv`).
- Writes only (`k275-l1/`): `checks-l1.tsv`, `timing-l1.tsv`, `hashes-l1.tsv`, `scorer-l1.tsv`, console lines, `g0.txt`. Row files stay under `tempdir()`.
- Decision: every check passes and every arm class has a completed journey at every design, then G2 re-plans the counted jobs from `timing-l1.tsv` and recomputes the ceiling against 350.73; else instrument defects are repaired inside the reserve under VAL-54's four conditions (at most three), results are never repaired.
