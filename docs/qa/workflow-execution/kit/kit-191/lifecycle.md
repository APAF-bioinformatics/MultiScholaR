# KIT-191 L1: the permutation test's time and memory

Run 45285106025 (run-tests, scripts job, ubuntu-22.04, reference BLAS), package 53290941, script `lifecycle-k191.R`;
tables in `l1/`. Each case in a fresh process: inputs built (random values, seed 20261191), `VmRSS` read immediately
before the call, the call timed, `VmHWM` read after. A measurement of the instrument's cost, not a result about the method.

## Observed

| Case | G | n | sets (sizes) | B | seconds | increase (MB) |
|---|---|---|---|---|---|---|
| a | 4,000 | 12 | 150 (5-50) | 1,000 | 1.1 | 127 |
| b1 | 4,000 | 12 | 5,000 (5-500) | 1,000 | 35.0 | 317 |
| b2 | 4,000 | 12 | 5,000 (5-500) | 10,000 | 259.1 | 393 |
| c1 | 20,000 | 500 | 2,000 (5-500) | 1,000 | 25.6 | 473 |
| c2 | 20,000 | 500 | 2,000 (5-500) | 10,000 | 232.0 | 561 |
| d | 20,000 | 4 | 20 (5) | 10,000 | 21.3 | 312 |

Every case completed with every p-value finite. Largest increase 561 MB: the stop rule (1 GB) did not fire. The
declared invalid call (`permutations = 0`) was refused (`integration_enrichment_parameters_invalid`).

## Coefficients (twice the seconds per unit of the term each case is dominated by)

- a (n G B, cases c1, c2): 2 x 2.5603e-9 = 5.12e-9
- c (G B, case d): 2 x 1.0632e-7 = 2.13e-7 (written by the script)
- s (sum of set sizes x B, cases b1, b2): 2 x 2.7417e-8 = 5.48e-8

Instrument defect, no re-run: the script wrote `NA` for a and s because the product of its integer columns overflowed
(500 x 20,000 x 10,000 > 2^31). The measured seconds are in `l1/lifecycle.tsv`; a and s were computed from them by the
script's own declared rule in doubles (c1: 25.603 / 1e10; c2: 231.972 / 1e11; b1: 34.956 / 1.274991e9; b2: 259.149 /
1.274991e10; the larger of each pair). No measurement was changed; the script now multiplies in doubles (2f54469e).
The coefficients are in `INTEGRATION_ENRICHMENT_PERMUTATION_COST` with the run id, and U12 asserts them.

## Inferred

The prediction over-states every case (attributing a whole case's time to its dominant term, then doubling): case a
predicts 1.3 s (measured 1.1), b2 710 s (259), c2 about 840 s (232). The study's size (case a) is far below the tab's
120 s limit; at Gene Ontology sizes with 10,000 reorderings the tab refuses and a script runs after the message.
The time of the 2 minutes the tab would wait was not measured in the app (ticket 14).
