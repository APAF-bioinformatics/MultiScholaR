# KIT-275 L1 (the blind lifecycle): declaration (see `declarations/l1-declaration.md` for the sent text)

Requests, in this order: **G1 first** (its own request: `fixture-k275.R D3,D4,D6 1 20 k275-g1`, then `score-k275.R --g1 k275-g1/g1.tsv k275-g1`, ceiling 30 job-minutes = 2.00 core-hours; not blind; run once; F3 stops everything), then **L1**:

1. `selftest/run-selftest.R` (parse, E0, E4), `timeout_minutes` 15;
2. `sim-k275.R all D3,D4,D6 9101 9103 k275-l1` with `K275_LIFECYCLE=1 K275_FAULT=1 K275_JOB=l1 K275_DECLARATION=<path> OMP_NUM_THREADS=1`, `timeout_minutes` 114; no smoke flag, no planted guard.

Upload pattern `\.(tsv|txt|json)$`. Ceiling 114 + 15 + 10 = 139 job-minutes = 9.27 core-hours. Writes only checks, timings, hashes, statuses and the blind scorer's lines; row files stay under `tempdir()`.
E7's fault and refusal sit on the first unit that E9 does not repeat (data set 9102).
