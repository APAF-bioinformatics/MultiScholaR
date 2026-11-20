# KIT-151 d, Q1: which RUV output feeds integration

Cells: 720, ok 720, failed 0 (failed: ). CPU models: AMD EPYC 7763 64-Core Processor (420 cells), AMD EPYC 9V45 96-Core Processor (120 cells), AMD EPYC 9V74 80-Core Processor (120 cells), Intel(R) Xeon(R) Platinum 8370C CPU @ 2.80GHz (60 cells). RUVIII_C_Varying against ruv::RUVIII, largest abs difference: 1.66e-08.
Fits at the round limit (10,000): 0; most rounds: 2451.

## Per arm (median over data sets; S and G recovered = matched abs r >= 0.9)

| N | scenario | arm | S recovered | median abs r S | mean abs r S [95%] | G recovered | group factor VE % | R2 vs group | -log10 p | S kept | G kept | U kept | k (calibrated, A; B) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | main | raw | 30/30 | 0.999 | 0.998 [0.997, 0.998] | 30/30 | 7.78 | 0.994 | 20.9 |    1 |    1 |    1 |  |
| 20 | main | ideal | 30/30 | 0.999 | 0.998 [0.998, 0.999] | 30/30 | 15.4 | 0.995 | 21.2 |    1 |    1 | 0.0929 |  |
| 20 | main | ruv_k1 | 29/30 | 0.986 | 0.975 [0.963, 0.984] | 30/30 | 16.2 | 0.995 | 21.3 | 0.985 | 0.999 | 0.05 |  |
| 20 | main | ruv_k2 | 0/30 | 0.00993 | 0.0105 [0.00867, 0.0124] | 30/30 | 23.3 | 0.996 | 22.1 | 0.00688 | 0.998 | 0.00289 |  |
| 20 | main | ruv_k3 | 0/30 | 0.0122 | 0.012 [0.01, 0.014] | 30/30 | 24.9 | 0.996 | 22.3 | 0.00641 | 0.997 | 0.00268 |  |
| 20 | main | calibrated | 0/30 | 0.00993 | 0.0105 [0.00867, 0.0124] | 30/30 | 23.3 | 0.996 | 22.1 | 0.00688 | 0.998 | 0.00289 | 2:30; 2:30 |
| 20 | conf | raw | 30/30 | 0.997 | 0.994 [0.991, 0.997] | 6/30 | 5.38 | 0.76 | 6.25 |    1 |    1 |    1 |  |
| 20 | conf | ideal | 30/30 | 0.999 | 0.998 [0.998, 0.999] | 30/30 | 15.4 | 0.995 | 21.2 |    1 |    1 | 0.116 |  |
| 20 | conf | ruv_k1 | 16/30 | 0.91 | 0.892 [0.856, 0.924] | 19/30 | 17.6 | 0.862 | 8.45 | 0.966 | 0.999 | 0.0811 |  |
| 20 | conf | ruv_k2 | 0/30 | 0.0106 | 0.0111 [0.00904, 0.0133] | 30/30 | 25.4 | 0.996 | 22.8 | 0.00688 | 0.999 | 0.00444 |  |
| 20 | conf | ruv_k3 | 0/30 | 0.0107 | 0.0113 [0.00928, 0.0135] | 30/30 | 27.1 | 0.997 | 23.1 | 0.00644 | 0.997 | 0.00417 |  |
| 20 | conf | calibrated | 0/30 | 0.0106 | 0.0111 [0.00904, 0.0133] | 30/30 | 25.4 | 0.996 | 22.8 | 0.00688 | 0.999 | 0.00444 | 2:30; 2:30 |
| 60 | main | raw | 30/30 | 0.999 | 0.999 [0.999, 0.999] | 30/30 | 6.78 | 0.996 | 70.3 |    1 |    1 |    1 |  |
| 60 | main | ideal | 30/30 | 0.999 | 0.999 [0.999, 0.999] | 30/30 | 13.3 | 0.996 | 70.9 |    1 |    1 | 0.0526 |  |
| 60 | main | ruv_k1 | 30/30 | 0.994 | 0.991 [0.988, 0.994] | 30/30 | 13.6 | 0.996 |   71 | 0.991 | 0.999 | 0.0219 |  |
| 60 | main | ruv_k2 | 0/30 | 0.00793 | 0.00799 [0.00639, 0.00963] | 30/30 | 18.8 | 0.997 | 72.5 | 0.00408 | 0.998 | 0.00164 |  |
| 60 | main | ruv_k3 | 0/30 | 0.00535 | 0.0062 [0.00481, 0.00774] | 30/30 | 19.3 | 0.997 | 72.6 | 0.00399 | 0.996 | 0.00163 |  |
| 60 | main | calibrated | 0/30 | 0.00793 | 0.00799 [0.00639, 0.00963] | 30/30 | 18.8 | 0.997 | 72.5 | 0.00408 | 0.998 | 0.00164 | 2:30; 2:30 |
| 60 | conf | raw | 30/30 | 0.999 | 0.999 [0.999, 0.999] | 30/30 | 9.35 | 0.983 | 52.4 |    1 |    1 |    1 |  |
| 60 | conf | ideal | 30/30 | 0.999 | 0.999 [0.999, 0.999] | 30/30 | 13.3 | 0.996 | 70.9 |    1 |    1 | 0.0658 |  |
| 60 | conf | ruv_k1 | 28/30 | 0.977 | 0.962 [0.943, 0.977] | 28/30 | 14.1 | 0.987 | 56.1 | 0.981 | 0.999 | 0.038 |  |
| 60 | conf | ruv_k2 | 0/30 | 0.00791 | 0.00792 [0.00633, 0.00956] | 30/30 | 19.6 | 0.997 | 73.1 | 0.00408 | 0.998 | 0.00255 |  |
| 60 | conf | ruv_k3 | 0/30 | 0.00586 | 0.00694 [0.00557, 0.0084] | 30/30 |   20 | 0.997 | 73.3 | 0.004 | 0.996 | 0.00255 |  |
| 60 | conf | calibrated | 0/30 | 0.00791 | 0.00792 [0.00633, 0.00956] | 30/30 | 19.6 | 0.997 | 73.1 | 0.00408 | 0.998 | 0.00255 | 2:30; 2:30 |

## Paired against ideal (b minus ideal; mean [95% seed-bootstrap interval], reading)

| N | scenario | arm | abs r S | abs r G | group VE (points) | R2 | -log10 p | S kept |
|---|---|---|---|---|---|---|---|---|
| 20 | main | raw | -0.001 [-0.001, -0.000] no practical difference | -0.000 [-0.000, +0.000] no practical difference | -7.542 [-7.602, -7.480] lower | -0.000 [-0.001, +0.000] inconclusive | -0.175 [-0.492, +0.120] inconclusive | +0.000 [+0.000, +0.000] inconclusive |
| 20 | main | ruv_k1 | -0.024 [-0.035, -0.014] worse | +0.000 [-0.000, +0.001] no practical difference | +0.826 [+0.744, +0.920] higher | +0.000 [-0.001, +0.002] inconclusive | +0.317 [-0.129, +0.900] inconclusive | -0.022 [-0.029, -0.015] lower |
| 20 | main | ruv_k2 | -0.988 [-0.990, -0.986] worse | +0.001 [+0.001, +0.002] no practical difference | +7.944 [+7.880, +8.005] higher | +0.002 [+0.001, +0.004] higher | +1.511 [+0.833, +2.313] higher | -0.993 [-0.994, -0.993] lower |
| 20 | main | ruv_k3 | -0.986 [-0.988, -0.984] worse | +0.001 [+0.001, +0.002] no practical difference | +9.557 [+9.481, +9.630] higher | +0.003 [+0.002, +0.004] higher | +1.807 [+1.207, +2.502] higher | -0.993 [-0.994, -0.993] lower |
| 20 | main | calibrated | -0.988 [-0.990, -0.986] worse | +0.001 [+0.001, +0.002] no practical difference | +7.944 [+7.880, +8.005] higher | +0.002 [+0.001, +0.004] higher | +1.511 [+0.833, +2.313] higher | -0.993 [-0.994, -0.993] lower |
| 20 | conf | raw | -0.004 [-0.007, -0.002] no practical difference | -0.161 [-0.191, -0.134] worse | -10.542 [-11.127, -9.990] lower | -0.289 [-0.336, -0.246] lower | -15.152 [-15.946, -14.359] lower | +0.000 [+0.000, +0.000] inconclusive |
| 20 | conf | ruv_k1 | -0.106 [-0.142, -0.074] worse | -0.095 [-0.129, -0.063] worse | +3.490 [+2.366, +4.759] higher | -0.171 [-0.231, -0.116] lower | -10.622 [-12.792, -8.306] lower | -0.044 [-0.056, -0.032] lower |
| 20 | conf | ruv_k2 | -0.987 [-0.989, -0.985] worse | +0.001 [+0.001, +0.002] no practical difference | +9.900 [+9.743, +10.058] higher | +0.003 [+0.002, +0.004] higher | +2.133 [+1.416, +2.942] higher | -0.993 [-0.994, -0.993] lower |
| 20 | conf | ruv_k3 | -0.987 [-0.989, -0.985] worse | +0.002 [+0.001, +0.002] no practical difference | +11.580 [+11.408, +11.751] higher | +0.003 [+0.002, +0.005] higher | +2.367 [+1.726, +3.097] higher | -0.993 [-0.994, -0.993] lower |
| 20 | conf | calibrated | -0.987 [-0.989, -0.985] worse | +0.001 [+0.001, +0.002] no practical difference | +9.900 [+9.743, +10.058] higher | +0.003 [+0.002, +0.004] higher | +2.133 [+1.416, +2.942] higher | -0.993 [-0.994, -0.993] lower |
| 60 | main | raw | -0.000 [-0.000, -0.000] no practical difference | -0.000 [-0.000, -0.000] no practical difference | -6.566 [-6.595, -6.538] lower | -0.000 [-0.000, -0.000] lower | -0.344 [-0.619, -0.090] lower | +0.000 [+0.000, +0.000] inconclusive |
| 60 | main | ruv_k1 | -0.008 [-0.011, -0.005] no practical difference | +0.000 [-0.000, +0.000] no practical difference | +0.234 [+0.202, +0.270] higher | +0.000 [-0.000, +0.000] inconclusive | +0.200 [-0.118, +0.536] inconclusive | -0.011 [-0.015, -0.007] lower |
| 60 | main | ruv_k2 | -0.991 [-0.993, -0.989] worse | +0.000 [+0.000, +0.001] no practical difference | +5.538 [+5.511, +5.566] higher | +0.001 [+0.000, +0.001] higher | +2.115 [+1.265, +3.107] higher | -0.996 [-0.996, -0.996] lower |
| 60 | main | ruv_k3 | -0.993 [-0.994, -0.991] worse | +0.000 [+0.000, +0.001] no practical difference | +5.965 [+5.933, +6.000] higher | +0.001 [+0.000, +0.001] higher | +2.426 [+1.526, +3.479] higher | -0.996 [-0.996, -0.996] lower |
| 60 | main | calibrated | -0.991 [-0.993, -0.989] worse | +0.000 [+0.000, +0.001] no practical difference | +5.538 [+5.511, +5.566] higher | +0.001 [+0.000, +0.001] higher | +2.115 [+1.265, +3.107] higher | -0.996 [-0.996, -0.996] lower |
| 60 | conf | raw | -0.000 [-0.000, +0.000] no practical difference | -0.006 [-0.007, -0.005] no practical difference | -4.100 [-4.251, -3.951] lower | -0.013 [-0.014, -0.011] lower | -17.296 [-18.803, -15.777] lower | +0.000 [+0.000, +0.000] inconclusive |
| 60 | conf | ruv_k1 | -0.037 [-0.056, -0.022] worse | -0.029 [-0.052, -0.011] worse | +1.293 [+0.836, +1.880] higher | -0.053 [-0.094, -0.021] lower | -18.512 [-25.608, -11.814] lower | -0.023 [-0.031, -0.016] lower |
| 60 | conf | ruv_k2 | -0.991 [-0.993, -0.989] worse | +0.000 [+0.000, +0.001] no practical difference | +6.220 [+6.156, +6.286] higher | +0.001 [+0.001, +0.001] higher | +2.825 [+1.936, +3.821] higher | -0.996 [-0.996, -0.996] lower |
| 60 | conf | ruv_k3 | -0.992 [-0.993, -0.991] worse | +0.001 [+0.000, +0.001] no practical difference | +6.657 [+6.590, +6.726] higher | +0.001 [+0.001, +0.001] higher | +3.101 [+2.160, +4.161] higher | -0.996 [-0.996, -0.996] lower |
| 60 | conf | calibrated | -0.991 [-0.993, -0.989] worse | +0.000 [+0.000, +0.001] no practical difference | +6.220 [+6.156, +6.286] higher | +0.001 [+0.001, +0.001] higher | +2.825 [+1.936, +3.821] higher | -0.996 [-0.996, -0.996] lower |

## Declared readings
- N = 20, main, raw: keeps the within-group factor (S recovered 30 of 30; ideal 30).
- N = 20, main, ruv_k1: keeps the within-group factor (S recovered 29 of 30; ideal 30).
- N = 20, main, ruv_k2: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 20, main, ruv_k3: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 20, main, calibrated: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 20, conf, raw: keeps the within-group factor (S recovered 30 of 30; ideal 30).
- N = 20, conf, ruv_k1: partly removes the within-group factor (S recovered 16 of 30; ideal 30).
- N = 20, conf, ruv_k2: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 20, conf, ruv_k3: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 20, conf, calibrated: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 60, main, raw: keeps the within-group factor (S recovered 30 of 30; ideal 30).
- N = 60, main, ruv_k1: keeps the within-group factor (S recovered 30 of 30; ideal 30).
- N = 60, main, ruv_k2: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 60, main, ruv_k3: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 60, main, calibrated: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 60, conf, raw: keeps the within-group factor (S recovered 30 of 30; ideal 30).
- N = 60, conf, ruv_k1: keeps the within-group factor (S recovered 28 of 30; ideal 30).
- N = 60, conf, ruv_k2: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 60, conf, ruv_k3: removes the within-group factor (S recovered 0 of 30; ideal 30).
- N = 60, conf, calibrated: removes the within-group factor (S recovered 0 of 30; ideal 30).
