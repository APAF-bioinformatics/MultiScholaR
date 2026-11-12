# KIT-120 memory finding: the bisect (2026-11-13)

The DIA full journey stopped after `design` (members `9xxd`, about five minutes each), application process (`controller_anon`, GiB) from each member's `memory-trace.tsv`. Builds `kit_924` to `kit_930` were installed for this diagnosis only; `kit_928`, `kit_929` and `kit_930` are from the side branches `kit120-diag-a/b/c` (never merged).

| Member | Build | Note | Import peak | Design minimum | Design end |
|---|---|---|---|---|---|
| 922d | kit_922, before the ticket | run first | 2.32 | 1.32 | 2.15 |
| 922e | kit_922 | run last | 2.32 | 1.32 | 2.20 |
| 926d | kit_926, W1 |  | 2.32 | 1.33 | 1.60 |
| 927d | kit_927, W3 |  | 2.32 | 1.32 | 2.20 |
| 924d | kit_924, W2 |  | 2.33 | 1.99 | 2.71 |
| 928d | kit_928, W3 + W2 config only | side branch | 2.31 | 2.06 | 2.66 |
| 929d | kit_929, W3 + W2 tab code only | side branch | 2.48 | 2.10 | 2.90 |
| 925d | kit_925, W4 (switch) |  | 2.30 | 2.05 | 2.77 |
| 923d | kit_923, W5 (exit) |  | 2.30 | 2.05 | 2.65 |
| 930d | kit_930, exit + gc(full = TRUE) after design save | side branch, diagnosis only | 2.30 | 1.87 | 1.87 |

- Releasing builds fall to about 1.32 GiB about 48 s into `design` (a collection while the design builder works); the others never fall.
- Running the before-ticket build last (922e) still releases: not an order effect.
- W1 and W3 release; W2 does not, and **each half of W2 alone stops it** (928d, 929d).
- **930d: one `gc(full = TRUE)` after the design save frees 0.74 GiB (2.61 to 1.87)**, ending the step lower than the before-ticket build (2.20). The memory was garbage, collectable at any time; nothing holds it.
