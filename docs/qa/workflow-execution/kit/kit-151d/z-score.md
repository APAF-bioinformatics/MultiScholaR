# KIT-151 d, Q2: z-score and scale_views

Cells: 1110, ok 1110, failed 0 (failed: ). CPU models: AMD EPYC 7763 64-Core Processor (300 cells), AMD EPYC 9V45 96-Core Processor (318 cells), AMD EPYC 9V74 80-Core Processor (196 cells), Intel(R) Xeon(R) 6973P-C (96 cells), Intel(R) Xeon(R) Platinum 8370C CPU @ 2.80GHz (26 cells), INTEL(R) XEON(R) PLATINUM 8573C (174 cells).
Fits at the round limit (10,000): 0; most rounds: 2366.

| N | cell | data sets | worst-of-three median (10th pct) | F1 | F2 | F3 | F2 >= 0.9 | iterations median (max) |
|---|---|---|---|---|---|---|---|---|
| 20 | m0_ldown_znone_sv0 | 30 | 0.966 (0.901) | 0.996 | 0.972 | 0.987 | 27 | 778.5 (1571) |
| 20 | m0_ldown_znone_sv1 | 30 | 0.964 (0.821) | 0.996 | 0.965 | 0.987 | 26 | 966 (1421) |
| 20 | m0_lup_znone_sv0 | 30 | 0.963 (0.306) | 0.996 | 0.967 | 0.987 | 23 | 273.5 (501) |
| 20 | m0_lup_znone_sv1 | 30 | 0.964 (0.821) | 0.996 | 0.965 | 0.987 | 26 | 966 (1421) |
| 20 | m0_none_znone_sv0 | 30 | 0.959 (0.324) | 0.996 | 0.963 | 0.987 | 21 | 631 (1141) |
| 20 | m0_none_znone_sv1 | 30 | 0.964 (0.821) | 0.996 | 0.965 | 0.987 | 26 | 966 (1421) |
| 20 | m0_orig_zall_sv0 | 30 | 0.969 (0.926) | 0.996 | 0.97 | 0.987 | 29 | 793.5 (1126) |
| 20 | m0_orig_zall_sv1 | 30 | 0.969 (0.926) | 0.996 | 0.97 | 0.987 | 29 | 786 (1116) |
| 20 | m0_orig_zlarge_sv1 | 30 | 0.968 (0.939) | 0.996 | 0.972 | 0.986 | 29 | 873.5 (1256) |
| 20 | m0_orig_znone_sv0 | 30 | 0.966 (0.943) | 0.996 | 0.971 | 0.987 | 28 | 491 (856) |
| 20 | m0_orig_znone_sv1 | 30 | 0.964 (0.821) | 0.996 | 0.965 | 0.987 | 26 | 966 (1421) |
| 20 | m0_orig_zsmall_sv1 | 30 | 0.97 (0.937) | 0.996 | 0.973 | 0.987 | 29 | 948.5 (1676) |
| 20 | m1_ldown_znone_sv0 | 30 | 0.596 (0.114) | 0.996 | 0.596 | 0.991 | 12 | 1021 (2036) |
| 20 | m1_ldown_znone_sv1 | 30 | 0.569 (0.264) | 0.996 | 0.569 | 0.991 | 12 | 1006 (1926) |
| 20 | m1_lup_znone_sv0 | 30 | 0.567 (0.203) | 0.996 | 0.567 | 0.991 | 11 | 451 (716) |
| 20 | m1_lup_znone_sv1 | 30 | 0.569 (0.264) | 0.996 | 0.569 | 0.991 | 12 | 1006 (1926) |
| 20 | m1_none_znone_sv0 | 30 | 0.622 (0.218) | 0.996 | 0.622 | 0.991 | 13 | 848.5 (1756) |
| 20 | m1_none_znone_sv1 | 30 | 0.569 (0.264) | 0.996 | 0.569 | 0.991 | 12 | 1006 (1926) |
| 20 | m1_orig_zall_sv0 | 30 | 0.958 (0.877) | 0.996 | 0.96 | 0.991 | 26 | 1043.5 (1806) |
| 20 | m1_orig_zall_sv1 | 30 | 0.958 (0.877) | 0.996 | 0.96 | 0.991 | 26 | 1056 (1796) |
| 20 | m1_orig_zlarge_sv1 | 30 | 0.821 (0.245) | 0.996 | 0.821 | 0.991 | 15 | 1023.5 (1836) |
| 20 | m1_orig_znone_sv0 | 30 | 0.758 (0.226) | 0.997 | 0.758 | 0.991 | 14 | 733.5 (1276) |
| 20 | m1_orig_znone_sv1 | 30 | 0.569 (0.264) | 0.996 | 0.569 | 0.991 | 12 | 1006 (1926) |
| 20 | m1_orig_zsmall_sv1 | 30 | 0.965 (0.923) | 0.996 | 0.967 | 0.992 | 28 | 1076 (2366) |
| 60 | m0_orig_zall_sv1 | 30 | 0.984 (0.976) | 0.998 | 0.984 | 0.996 | 30 | 1028.5 (1771) |
| 60 | m0_orig_znone_sv1 | 30 | 0.985 (0.978) | 0.998 | 0.985 | 0.996 | 30 | 1081 (1651) |
| 60 | m0_orig_zsmall_sv1 | 30 | 0.985 (0.976) | 0.998 | 0.985 | 0.996 | 30 | 1068.5 (1951) |
| 60 | m1_ldown_znone_sv0 | 30 | 0.979 (0.969) | 0.999 | 0.979 | 0.997 | 30 | 1073.5 (1751) |
| 60 | m1_ldown_znone_sv1 | 30 | 0.98 (0.969) | 0.999 | 0.98 | 0.997 | 30 | 1196 (1751) |
| 60 | m1_lup_znone_sv0 | 30 | 0.979 (0.97) | 0.999 | 0.979 | 0.997 | 30 | 491 (1216) |
| 60 | m1_lup_znone_sv1 | 30 | 0.98 (0.969) | 0.999 | 0.98 | 0.997 | 30 | 1196 (1751) |
| 60 | m1_none_znone_sv0 | 30 | 0.979 (0.969) | 0.999 | 0.979 | 0.997 | 30 | 1016 (1941) |
| 60 | m1_none_znone_sv1 | 30 | 0.98 (0.969) | 0.999 | 0.98 | 0.997 | 30 | 1196 (1751) |
| 60 | m1_orig_zall_sv1 | 30 | 0.979 (0.967) | 0.999 | 0.979 | 0.997 | 30 | 1206 (2131) |
| 60 | m1_orig_znone_sv0 | 30 | 0.979 (0.969) | 0.999 | 0.979 | 0.997 | 30 | 783.5 (1666) |
| 60 | m1_orig_znone_sv1 | 30 | 0.98 (0.969) | 0.999 | 0.98 | 0.997 | 30 | 1196 (1751) |
| 60 | m1_orig_zsmall_sv1 | 30 | 0.979 (0.967) | 0.999 | 0.979 | 0.997 | 30 | 1113.5 (1881) |

## Declared comparisons (b minus a; mean [95% seed-bootstrap interval], reading with margin 0.02)

N = 20

| comparison | worst-of-three | F2 | F3 | F1 |
|---|---|---|---|---|
| Z1 z-score all vs none, missing on | +0.301 [+0.200, +0.406] better (26/4) | +0.301 [+0.200, +0.406] better (26/4) | -0.001 [-0.001, -0.000] no practical difference (11/19) | -0.000 [-0.000, +0.000] no practical difference (16/14) |
| Z2 z-score all vs none, missing off | +0.058 [-0.001, +0.132] inconclusive (17/13) | +0.057 [-0.002, +0.131] inconclusive (15/15) | +0.000 [-0.000, +0.001] no practical difference (16/14) | -0.000 [-0.000, +0.000] no practical difference (13/17) |
| Z4 z-score small data sets vs none, missing on | +0.317 [+0.214, +0.422] better (29/1) | +0.317 [+0.214, +0.422] better (29/1) | -0.000 [-0.000, +0.000] no practical difference (15/15) | +0.000 [-0.000, +0.000] no practical difference (17/13) |
| Z5 z-score large data sets vs none, missing on | +0.044 [-0.023, +0.123] inconclusive (15/15) | +0.044 [-0.023, +0.123] inconclusive (14/16) | -0.000 [-0.001, +0.000] no practical difference (11/19) | -0.000 [-0.000, +0.000] no practical difference (11/19) |
| Z6 z-score all vs small only, missing on | -0.016 [-0.050, +0.004] inconclusive (8/22) | -0.016 [-0.050, +0.004] inconclusive (8/22) | -0.000 [-0.001, +0.000] no practical difference (11/19) | -0.000 [-0.000, +0.000] no practical difference (11/19) |
| Z7 z-score small data sets vs none, missing off | +0.061 [+0.002, +0.135] better (20/10) | +0.061 [+0.002, +0.135] better (18/12) | +0.000 [-0.000, +0.001] no practical difference (21/9) | +0.000 [-0.000, +0.000] no practical difference (17/13) |
| Z8 z-score large data sets vs none, missing off | +0.058 [-0.002, +0.132] inconclusive (15/15) | +0.057 [-0.002, +0.132] inconclusive (13/17) | +0.000 [-0.000, +0.001] no practical difference (14/16) | -0.000 [-0.000, -0.000] no practical difference (7/23) |
| Z9 z-score all vs none, scale_views off, missing on | +0.239 [+0.143, +0.342] better (24/6) | +0.240 [+0.143, +0.342] better (24/6) | -0.000 [-0.001, +0.000] no practical difference (9/21) | -0.000 [-0.000, +0.000] no practical difference (13/17) |
| Z10 z-score all: scale_views on vs off, missing on | -0.000 [-0.001, -0.000] no practical difference (6/24) | -0.000 [-0.001, -0.000] no practical difference (5/25) | -0.000 [-0.000, +0.000] no practical difference (16/14) | -0.000 [-0.000, +0.000] no practical difference (10/20) |
| Z11 z-score all vs none, scale_views off, missing off | +0.020 [-0.004, +0.065] inconclusive (15/15) | +0.019 [-0.005, +0.064] inconclusive (14/16) | +0.000 [-0.000, +0.001] no practical difference (16/14) | -0.000 [-0.000, +0.000] no practical difference (16/14) |
| S scale_views on vs off, scales none, missing on | -0.014 [-0.084, +0.046] inconclusive (13/17) | -0.014 [-0.084, +0.046] inconclusive (13/17) | +0.000 [-0.000, +0.001] no practical difference (18/12) | -0.000 [-0.000, +0.000] no practical difference (16/14) |
| S scale_views on vs off, scales orig, missing on | -0.062 [-0.145, +0.010] inconclusive (13/17) | -0.062 [-0.145, +0.010] inconclusive (13/17) | +0.000 [-0.000, +0.001] no practical difference (15/15) | -0.000 [-0.000, +0.000] no practical difference (16/14) |
| S scale_views on vs off, scales lup, missing on | +0.031 [-0.079, +0.143] inconclusive (15/15) | +0.031 [-0.079, +0.143] inconclusive (16/14) | +0.000 [-0.000, +0.001] no practical difference (14/16) | -0.000 [-0.000, +0.000] no practical difference (14/16) |
| S scale_views on vs off, scales ldown, missing on | +0.028 [-0.083, +0.140] inconclusive (12/18) | +0.028 [-0.083, +0.140] inconclusive (12/18) | -0.000 [-0.000, +0.000] no practical difference (11/19) | -0.000 [-0.000, +0.000] no practical difference (12/18) |
| S scale_views on vs off, scales none, missing off | +0.107 [+0.025, +0.203] better (18/12) | +0.107 [+0.026, +0.204] better (20/10) | -0.000 [-0.000, +0.000] no practical difference (13/17) | -0.000 [-0.000, +0.000] no practical difference (14/16) |
| S scale_views on vs off, scales orig, missing off | -0.038 [-0.101, +0.006] inconclusive (13/17) | -0.038 [-0.101, +0.006] inconclusive (16/14) | -0.000 [-0.000, +0.000] no practical difference (17/13) | +0.000 [-0.000, +0.000] no practical difference (20/10) |
| S scale_views on vs off, scales lup, missing off | +0.084 [-0.013, +0.195] inconclusive (19/11) | +0.084 [-0.013, +0.195] inconclusive (20/10) | -0.000 [-0.000, +0.000] no practical difference (16/14) | +0.000 [-0.000, +0.000] no practical difference (17/13) |
| S scale_views on vs off, scales ldown, missing off | -0.017 [-0.079, +0.045] inconclusive (16/14) | -0.017 [-0.079, +0.045] inconclusive (19/11) | -0.000 [-0.000, +0.000] no practical difference (10/20) | +0.000 [-0.000, +0.000] no practical difference (17/13) |
| Z3 z-score gain, missing on minus missing off | +0.243 [+0.124, +0.358] better (25/5) | +0.244 [+0.125, +0.359] better (25/5) | -0.001 [-0.002, +0.000] no practical difference (12/18) | +0.000 [-0.000, +0.000] no practical difference (15/15) |
| SX scale_views effect at scales orig minus at none, missing on | -0.048 [-0.112, +0.003] inconclusive (13/17) | -0.048 [-0.112, +0.003] inconclusive (13/17) | +0.000 [-0.000, +0.000] no practical difference (11/19) | -0.000 [-0.000, +0.000] no practical difference (15/15) |
| SX scale_views effect at scales lup minus at none, missing on | +0.045 [-0.061, +0.158] inconclusive (16/14) | +0.045 [-0.061, +0.158] inconclusive (17/13) | -0.000 [-0.000, +0.000] no practical difference (10/20) | -0.000 [-0.000, +0.000] no practical difference (15/15) |
| SX scale_views effect at scales ldown minus at none, missing on | +0.042 [-0.072, +0.165] inconclusive (14/16) | +0.042 [-0.072, +0.165] inconclusive (14/16) | -0.000 [-0.001, +0.000] no practical difference (9/21) | -0.000 [-0.000, +0.000] no practical difference (15/15) |
| SX scale_views effect at scales orig minus at none, missing off | -0.145 [-0.249, -0.054] worse (14/16) | -0.145 [-0.250, -0.054] worse (14/16) | +0.000 [-0.000, +0.000] no practical difference (18/12) | +0.000 [-0.000, +0.000] no practical difference (18/12) |
| SX scale_views effect at scales lup minus at none, missing off | -0.022 [-0.120, +0.075] inconclusive (20/10) | -0.023 [-0.120, +0.075] inconclusive (18/12) | +0.000 [-0.000, +0.000] no practical difference (19/11) | +0.000 [-0.000, +0.000] no practical difference (18/12) |
| SX scale_views effect at scales ldown minus at none, missing off | -0.124 [-0.228, -0.033] worse (12/18) | -0.125 [-0.229, -0.033] worse (11/19) | +0.000 [-0.000, +0.001] no practical difference (16/14) | +0.000 [-0.000, +0.000] no practical difference (18/12) |

N = 60

| comparison | worst-of-three | F2 | F3 | F1 |
|---|---|---|---|---|
| Z1 z-score all vs none, missing on | -0.000 [-0.001, +0.001] no practical difference (15/15) | -0.000 [-0.001, +0.001] no practical difference (15/15) | -0.000 [-0.001, -0.000] no practical difference (3/27) | +0.000 [-0.000, +0.000] no practical difference (17/13) |
| Z2 z-score all vs none, missing off | -0.001 [-0.002, -0.000] no practical difference (14/16) | -0.001 [-0.002, -0.000] no practical difference (13/17) | +0.000 [+0.000, +0.000] no practical difference (16/14) | +0.000 [-0.000, +0.000] no practical difference (15/15) |
| Z4 z-score small data sets vs none, missing on | -0.000 [-0.001, +0.001] no practical difference (16/14) | -0.000 [-0.001, +0.001] no practical difference (16/14) | -0.000 [-0.000, +0.000] no practical difference (16/14) | +0.000 [-0.000, +0.000] no practical difference (15/15) |
| Z6 z-score all vs small only, missing on | +0.000 [-0.000, +0.000] no practical difference (16/14) | +0.000 [-0.000, +0.000] no practical difference (16/14) | -0.000 [-0.001, -0.000] no practical difference (2/28) | +0.000 [-0.000, +0.000] no practical difference (17/13) |
| Z7 z-score small data sets vs none, missing off | -0.001 [-0.001, -0.000] no practical difference (15/15) | -0.001 [-0.001, -0.000] no practical difference (14/16) | +0.000 [-0.000, +0.000] no practical difference (18/12) | +0.000 [+0.000, +0.000] no practical difference (18/12) |
| S scale_views on vs off, scales none, missing on | -0.000 [-0.000, +0.000] no practical difference (18/12) | -0.000 [-0.000, +0.000] no practical difference (18/12) | +0.000 [-0.000, +0.000] no practical difference (12/18) | +0.000 [-0.000, +0.000] no practical difference (20/10) |
| S scale_views on vs off, scales orig, missing on | -0.000 [-0.000, +0.000] no practical difference (13/17) | -0.000 [-0.000, +0.000] no practical difference (13/17) | -0.000 [-0.000, +0.000] no practical difference (15/15) | -0.000 [-0.000, +0.000] no practical difference (11/19) |
| S scale_views on vs off, scales lup, missing on | +0.000 [+0.000, +0.000] no practical difference (21/9) | +0.000 [+0.000, +0.000] no practical difference (21/9) | -0.000 [-0.000, +0.000] no practical difference (15/15) | -0.000 [-0.000, +0.000] no practical difference (15/15) |
| S scale_views on vs off, scales ldown, missing on | -0.000 [-0.000, +0.000] no practical difference (18/12) | -0.000 [-0.000, +0.000] no practical difference (18/12) | +0.000 [-0.000, +0.000] no practical difference (16/14) | +0.000 [-0.000, +0.000] no practical difference (14/16) |
| Z3 z-score gain, missing on minus missing off | +0.001 [+0.000, +0.001] no practical difference (20/10) | +0.001 [+0.000, +0.001] no practical difference (20/10) | -0.000 [-0.001, -0.000] no practical difference (8/22) | +0.000 [-0.000, +0.000] no practical difference (17/13) |
| SX scale_views effect at scales orig minus at none, missing on | -0.000 [-0.000, +0.000] no practical difference (8/22) | -0.000 [-0.000, +0.000] no practical difference (8/22) | -0.000 [-0.000, +0.000] no practical difference (10/20) | -0.000 [-0.000, +0.000] no practical difference (11/19) |
| SX scale_views effect at scales lup minus at none, missing on | +0.000 [+0.000, +0.000] no practical difference (18/12) | +0.000 [+0.000, +0.000] no practical difference (18/12) | -0.000 [-0.000, +0.000] no practical difference (11/19) | -0.000 [-0.000, -0.000] no practical difference (12/18) |
| SX scale_views effect at scales ldown minus at none, missing on | -0.000 [-0.000, +0.000] no practical difference (14/16) | -0.000 [-0.000, +0.000] no practical difference (14/16) | -0.000 [-0.000, +0.000] no practical difference (16/14) | -0.000 [-0.000, +0.000] no practical difference (13/17) |

scale_views on, scales none against each other pattern (abs r of F1 to F3; expected near 0): N 20 m1 m1_orig_znone_sv1 max    0 (100% within 1e-6); N 20 m1 m1_lup_znone_sv1 max    0 (100% within 1e-6); N 20 m1 m1_ldown_znone_sv1 max    0 (100% within 1e-6); N 20 m0 m0_orig_znone_sv1 max    0 (100% within 1e-6); N 20 m0 m0_lup_znone_sv1 max    0 (100% within 1e-6); N 20 m0 m0_ldown_znone_sv1 max    0 (100% within 1e-6); N 60 m1 m1_orig_znone_sv1 max    0 (100% within 1e-6); N 60 m1 m1_lup_znone_sv1 max    0 (100% within 1e-6); N 60 m1 m1_ldown_znone_sv1 max    0 (100% within 1e-6)

Tie to k151b (instrumentation, not a gate): 156 of 180 abs r values agree within 1e-6 (same CPU model: 42 of 42); iterations equal in 156 of 180; largest difference 0.316.
