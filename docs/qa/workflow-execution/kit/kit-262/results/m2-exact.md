# Exact-binomial counterpart of the declared gain rule (description only; the declared rule and its verdict stand)

Nominal level per cell: 0.05 / 180 = 0.000278. Smallest count of 40 whose upper tail P(X >= x | p = 0.05) is at most that: **9** (P = 0.000130).
Declared Wilson rule (>= 7 of 40): a 17, b 19, e 1 cells. Exact test (>= 9 of 40): a 14, b 16, e 0 cells.

Operating characteristic of the declared rule (score test at R = 40): per-cell size at a true share of 0.05 is P(X >= 7) = 0.0034 against a nominal 0.000278;
family-wise over 180 cells at a true share of 0.05 everywhere: 1 - (1 - 0.0034)^180 = 0.46 (the ticket's own figure 0.46); at 0.03: 0.030; at 0.02: 0.0024.
The rule therefore controls error near a true share of 0.03, not 0.05. The extra gain of C2 over C1 (e) rests on one cell (82) at 8 of 40.

Cells with a under the exact test: 45, 63, 79, 80, 97, 113, 114, 131, 132, 147, 159, 161, 165, 179
