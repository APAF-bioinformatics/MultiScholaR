# KIT-155 measurement summary

Rows: 72; instrumentation invalid: 0

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
| L | before | write | completed | 126.3 | 1765.1 | 13.8 | 2 |
| L | before | import | completed | 589.2 | 1568.1 | 77.7 | 2 |
| L | before | align | completed | 746.1 | 1560.0 | 58.5 | 2 |
| L | before | train | completed | 848.2 | 1561.6 | 649.7 | 2 |
| L | after | write | completed | 126.3 | 1767.4 | 13.2 | 2 |
| L | after | import | completed | 589.2 | 1577.2 | 80.0 | 1 |
| L | after | align | completed | 457.9 | 1171.2 | 74.5 | 1 |
| L | after | train | completed | 560.0 | 1268.2 | 682.4 | 1 |

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
| L | KIT-155: alignment pages by bytes | align | R heap peak (MiB) | 746.1 | 457.9 | -288.2 | 0.0 noise of M | 64.0 | reduced |
| L | KIT-155: alignment pages by bytes | align | process peak (MiB) | 1560.0 | 1171.2 | -388.8 | 18.0 noise of M | 78.0 | reduced |
| L | KIT-155: alignment pages by bytes | align | wall (s) | 58.5 | 74.5 | 16.0 | 0.7 noise of M | 2.9 | increased |
| L | KIT-155: alignment pages by bytes | train | R heap peak (MiB) | 848.2 | 560.0 | -288.2 | 0.0 noise of M | 64.0 | reduced |
| L | KIT-155: alignment pages by bytes | train | process peak (MiB) | 1561.6 | 1268.2 | -293.4 | 5.9 noise of M | 78.1 | reduced |
| L | KIT-155: alignment pages by bytes | train | wall (s) | 649.7 | 682.4 | 32.7 | 5.5 noise of M | 32.5 | increased |

## Claim 1

| Claim | Build | Size | Verdict | Detail |
|---|---|---|---|---|
| 1 as declared | before | M | fail | heap peak 541.8 MiB vs bound 259.0 MiB; heap peak 541.8 MiB vs bound 259.0 MiB; heap peak 541.8 MiB vs bound 259.0 MiB |
| 1 restated (INT-27) | before | M | fail | overhead 261.8 MiB (before train at S, n=3); growth 280.0, 280.0, 280.0 MiB vs bound 259.0 MiB + 2 x noise 0.0 |
| 1 as declared | before | L | fail | heap peak 848.2 MiB vs bound 263.6 MiB; heap peak 848.2 MiB vs bound 263.6 MiB |
| 1 restated (INT-27) | before | L | fail | overhead 261.8 MiB (before train at S, n=3); growth 586.4, 586.4 MiB vs bound 263.6 MiB + 2 x noise 0.0 |
| 1 as declared | after | M | fail | heap peak 463.9 MiB vs bound 259.0 MiB; heap peak 463.9 MiB vs bound 259.0 MiB; heap peak 463.9 MiB vs bound 259.0 MiB |
| 1 restated (INT-27) | after | M | pass | overhead 255.3 MiB (after train at S, n=3); growth 208.6, 208.6, 208.6 MiB vs bound 259.0 MiB + 2 x noise 0.0 |
| 1 as declared | after | L | fail | heap peak 560.0 MiB vs bound 263.6 MiB |
| 1 restated (INT-27) | after | L | fail | overhead 255.3 MiB (after train at S, n=3); growth 304.7 MiB vs bound 263.6 MiB + 2 x noise 0.0 |

## Python against the dense floor (reported, not claimed)

| Size | Build | Case | Python peak MiB | Dense floor MiB | Largest selected MiB | Process peak MiB | Iterations |
|---|---|---|---|---|---|---|---|
| S | before | capped | 62.2 | 0.6 | 0.2 | 962.7 | 2046 |
| S | after | capped | 61.7 | 0.6 | 0.2 | 957.3 | 2046 |
| M | before | capped | 83.5 | 6.1 | 1.5 | 1212.9 | 3916 |
| M | after | capped | 83.5 | 6.1 | 1.5 | 1210.9 | 3916 |
| L | before | capped | 121.4 | 15.3 | 3.8 | 1564.3 | 4786 |
| L | after | capped | 121.4 | 15.3 | 3.8 | 1268.2 | 4786 |

## Model sanity

- S train: the same factors (rounded to 1e-8) and iterations in every build
- M train: the same factors (rounded to 1e-8) and iterations in every build
- L train: the same factors (rounded to 1e-8) and iterations in every build
