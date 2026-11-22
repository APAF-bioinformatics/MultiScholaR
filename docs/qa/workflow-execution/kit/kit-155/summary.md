# KIT-155 measurement summary

Rows: 84; instrumentation invalid: 0

## Stages (median of repeats; range in brackets when 3 repeats)

| Size | Build | Stage | Outcome | R heap peak MiB (range) | Process peak MiB (range) | Wall s (range) | n |
|---|---|---|---|---|---|---|---|
| S | before | write | completed | 100.4 (0.0) | 702.9 (0.3) | 0.2 (0.0) | 3 |
| S | before | import | completed | 356.2 (0.0) | 987.2 (3.4) | 2.3 (0.2) | 3 |
| S | before | align | completed | 261.8 (0.0) | 877.6 (0.2) | 1.0 (0.1) | 3 |
| S | before | train | completed | 261.8 (0.0) | 962.0 (1.8) | 26.9 (1.0) | 3 |
| S | after | write | completed | 100.5 (0.0) | 702.8 (0.1) | 0.2 (0.0) | 3 |
| S | after | import | completed | 356.2 (0.0) | 984.8 (14.1) | 2.2 (0.2) | 3 |
| S | after | align | completed | 255.4 (0.0) | 871.3 (0.2) | 0.9 (0.1) | 3 |
| S | after | train | completed | 255.3 (0.0) | 957.3 (0.8) | 27.4 (0.2) | 3 |
| M | before | write | completed | 110.4 (0.0) | 933.8 (1.3) | 2.9 (0.1) | 3 |
| M | before | import | completed | 471.9 (0.0) | 1161.7 (4.0) | 23.9 (1.6) | 3 |
| M | before | align | completed | 413.0 (0.0) | 1182.2 (7.2) | 16.3 (0.7) | 3 |
| M | before | train | completed | 541.8 (0.0) | 1212.9 (5.9) | 217.4 (5.5) | 3 |
| M | after | write | completed | 110.4 (0.0) | 935.6 (1.9) | 2.8 (0.1) | 3 |
| M | after | import | completed | 471.9 (0.0) | 1161.6 (13.3) | 23.4 (1.8) | 3 |
| M | after | align | completed | 335.1 (0.0) | 1074.6 (18.0) | 17.2 (0.7) | 3 |
| M | after | train | completed | 463.9 (0.0) | 1210.9 (5.3) | 218.4 (0.8) | 3 |
| L | before | write | completed | 126.3 (0.0) | 1766.3 (4.6) | 13.8 (4.2) | 3 |
| L | before | import | completed | 589.2 (0.0) | 1567.4 (6.6) | 79.2 (33.5) | 3 |
| L | before | align | completed | 746.1 (0.0) | 1560.7 (12.1) | 59.2 (9.0) | 3 |
| L | before | train | completed | 848.2 (0.0) | 1558.9 (15.2) | 652.9 (160.2) | 3 |
| L | after | write | completed | 126.3 (0.0) | 1768.7 (2.8) | 15.4 (4.6) | 4 |
| L | after | import | completed | 589.2 (0.0) | 1571.3 (11.1) | 103.2 (25.1) | 3 |
| L | after | align | completed | 456.8 (1.2) | 1176.9 (7.5) | 82.1 (7.9) | 3 |
| L | after | train | completed | 558.9 (1.2) | 1265.8 (4.7) | 833.3 (161.8) | 3 |

## Per mechanism (median; tolerance max(64 MiB or 1 s, 5 %, 2 x noise floor))

| Size | Mechanism | Stage | Metric | Before | After | Change | Noise | Tolerance | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| S | KIT-155: alignment pages by bytes | align | R heap peak (MiB) | 261.8 | 255.4 | -6.4 | 0.0  | 64.0 | no measurable change |
| S | KIT-155: alignment pages by bytes | align | process peak (MiB) | 877.6 | 871.3 | -6.3 | 0.2  | 64.0 | no measurable change |
| S | KIT-155: alignment pages by bytes | align | wall (s) | 1.0 | 0.9 | -0.1 | 0.1  | 1.0 | no measurable change |
| S | KIT-155: alignment pages by bytes | train | R heap peak (MiB) | 261.8 | 255.3 | -6.5 | 0.0  | 64.0 | no measurable change |
| S | KIT-155: alignment pages by bytes | train | process peak (MiB) | 962.0 | 957.3 | -4.7 | 1.8  | 64.0 | no measurable change |
| S | KIT-155: alignment pages by bytes | train | wall (s) | 26.9 | 27.4 | 0.5 | 1.0  | 2.0 | no measurable change |
| M | KIT-155: alignment pages by bytes | align | R heap peak (MiB) | 413.0 | 335.1 | -77.9 | 0.0  | 64.0 | reduced |
| M | KIT-155: alignment pages by bytes | align | process peak (MiB) | 1182.2 | 1074.6 | -107.6 | 18.0  | 64.0 | reduced |
| M | KIT-155: alignment pages by bytes | align | wall (s) | 16.3 | 17.2 | 0.9 | 0.7  | 1.4 | no measurable change |
| M | KIT-155: alignment pages by bytes | train | R heap peak (MiB) | 541.8 | 463.9 | -77.9 | 0.0  | 64.0 | reduced |
| M | KIT-155: alignment pages by bytes | train | process peak (MiB) | 1212.9 | 1210.9 | -2.0 | 5.9  | 64.0 | no measurable change |
| M | KIT-155: alignment pages by bytes | train | wall (s) | 217.4 | 218.4 | 1.0 | 5.5  | 11.0 | no measurable change |
| L | KIT-155: alignment pages by bytes | align | R heap peak (MiB) | 746.1 | 456.8 | -289.3 | 1.2  | 64.0 | reduced |
| L | KIT-155: alignment pages by bytes | align | process peak (MiB) | 1560.7 | 1176.9 | -383.8 | 12.1  | 78.0 | reduced |
| L | KIT-155: alignment pages by bytes | align | wall (s) | 59.2 | 82.1 | 22.9 | 9.0  | 18.0 | increased |
| L | KIT-155: alignment pages by bytes | train | R heap peak (MiB) | 848.2 | 558.9 | -289.3 | 1.2  | 64.0 | reduced |
| L | KIT-155: alignment pages by bytes | train | process peak (MiB) | 1558.9 | 1265.8 | -293.1 | 15.2  | 77.9 | reduced |
| L | KIT-155: alignment pages by bytes | train | wall (s) | 652.9 | 833.3 | 180.4 | 161.8  | 323.6 | no measurable change |

## Claim 1

| Claim | Build | Size | Verdict | Detail |
|---|---|---|---|---|
| 1 as declared | before | M | fail | heap peak 541.8 MiB vs bound 259.0 MiB; heap peak 541.8 MiB vs bound 259.0 MiB; heap peak 541.8 MiB vs bound 259.0 MiB |
| 1 restated (INT-27) | before | M | fail | overhead 261.8 MiB (before train at S, n=3); growth 280.0, 280.0, 280.0 MiB vs bound 259.0 MiB + 2 x noise 0.0 |
| 1 as declared | before | L | fail | heap peak 848.2 MiB vs bound 263.6 MiB; heap peak 848.2 MiB vs bound 263.6 MiB; heap peak 848.2 MiB vs bound 263.6 MiB |
| 1 restated (INT-27) | before | L | fail | overhead 261.8 MiB (before train at S, n=3); growth 586.4, 586.4, 586.4 MiB vs bound 263.6 MiB + 2 x noise 0.0 |
| 1 as declared | after | M | fail | heap peak 463.9 MiB vs bound 259.0 MiB; heap peak 463.9 MiB vs bound 259.0 MiB; heap peak 463.9 MiB vs bound 259.0 MiB |
| 1 restated (INT-27) | after | M | pass | overhead 255.3 MiB (after train at S, n=3); growth 208.6, 208.6, 208.6 MiB vs bound 259.0 MiB + 2 x noise 0.0 |
| 1 as declared | after | L | fail | heap peak 558.9 MiB vs bound 263.6 MiB; heap peak 558.8 MiB vs bound 263.6 MiB; heap peak 560.0 MiB vs bound 263.6 MiB |
| 1 restated (INT-27) | after | L | fail | overhead 255.3 MiB (after train at S, n=3); growth 303.6, 303.5, 304.7 MiB vs bound 263.6 MiB + 2 x noise 1.2 |

## Python against the dense floor (reported, not claimed)

| Size | Build | Case | Python peak MiB | Dense floor MiB | Largest selected MiB | Process peak MiB | Iterations |
|---|---|---|---|---|---|---|---|
| S | before | capped | 62.2 | 0.6 | 0.2 | 962.7 | 2046 |
| S | after | capped | 61.7 | 0.6 | 0.2 | 957.3 | 2046 |
| M | before | capped | 83.5 | 6.1 | 1.5 | 1212.9 | 3916 |
| M | after | capped | 83.5 | 6.1 | 1.5 | 1210.9 | 3916 |
| L | before | capped | 121.4 | 15.3 | 3.8 | 1549.1 | 4361 |
| L | after | capped | 121.8 | 15.3 | 3.8 | 1263.5 | 4361 |

## Model sanity

- S train: the same factors (rounded to 1e-8) and iterations in every build
- M train: the same factors (rounded to 1e-8) and iterations in every build
- L train: DIFFERENT: after [('18f39bdcd7d0a574', 4786), ('5bf078ac9cb00f1f', 4361)]; before [('18f39bdcd7d0a574', 4786), ('5bf078ac9cb00f1f', 4361)]
