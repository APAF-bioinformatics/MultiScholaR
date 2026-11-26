# KIT-211 results (the declared labels only; the owner edits the wording, VAL-41)

Scoring run 45730910292 (`score-k211.R`, no environment variables), on the rows of the ten counted jobs at package commit
656b08c65cd3cac19af0aab3725178b7b9ffb7fd. Checks read before the outcome (`results/checks.tsv`): complete 36,000 of 36,000 rows; ranges incomplete 0;
F1 0; F2 0; t_r monotone failures 0; BH fallbacks 0; determinism (S-A1, cell A r 201 to 205) 5 of 5 identical; bootstrap_B 10,000 (also in `outcome.tsv`).
R = 200 for every cell (A, B, C, D, Av, Ag).

**Outcome: K2** (the positive control passes; classes A = C, B = E, B′ = E; no candidate withheld). Expected before the run: K2 with A E or C.

## Primary, positive control, classes (declared rules; labels.tsv)

Positive control: S's Φ is 1.4897; its one-sided 97.5 % lower limit is 1.4301, above M = 1.20: passes.
The classes use one-sided 99.17 % limits (0.025/3); the decision table below shows two-sided 95 % intervals, so read the limits here beside each class (DR-05):

| Arm | class | primary state (one-sided 99.17 %) | phi_lower | phi (estimate) | phi_upper | power state | dpow (S minus X) [one-sided limits] | guard state | dfdr (X minus S) [one-sided limits] |
|---|---|---|---|---|---|---|---|---|---|
| A | C | removes | 0.7636 | 0.8405 | 0.9183 | not determined | 0.0049 [0.0047, 0.0052] | not raised | 0.000002 [-0.000109, 0.000108] |
| B | E | removes | 0.9872 | 1.0181 | 1.0490 | no material cost | -0.0037 [-0.0039, -0.0035] | not raised | 0.000070 [-0.000012, 0.000152] |
| B′ | E | removes | 1.1042 | 1.1416 | 1.1791 | no material cost | -0.0032 [-0.0034, -0.0030] | not raised | 0.000117 [0.000035, 0.000197] |

Jackknife-normal labels equal the bootstrap labels for every candidate (`primary_differs_jack` FALSE for all three).
Per-setting labels of Φ at one-sided 0.05/72 are in `results/labels-settings.tsv` (reported, not deciding). Deciding settings whose own label differs from the pooled one (22 setting-candidate labels): A vd5 A (excess remains), A vd5 B (not determined), A vd5 B′ (excess remains), A vd10 A (not determined), A vd10 B (not determined), A vd10 B′ (excess remains), A vd25 A (not determined), A vd25 B (not determined), A vd25 B′ (not determined), B vd5 B′ (excess remains), B vd10 A (not determined), B vd10 B (excess remains), B vd10 B′ (excess remains), B vd25 A (not determined), B vd25 B (excess remains), B vd25 B′ (excess remains), C vd10 B′ (excess remains), C vd25 A (not determined), C vd25 B (excess remains), C vd25 B′ (excess remains), D vd25 B (not determined), D vd25 B′ (not determined).

## Falsifiers (pooled over the 24 settings; per-setting Bonferroni intervals reported, not deciding)

| Falsifier | statistic | pooled estimate | two-sided limits at 0.05/4 | fires | settings whose Bonferroni (24) interval excludes 0 |
|---|---|---|---|---|---|
| F3 | A.m6_l_H | -0.00481 | [-0.01812, 0.00771] | FALSE | 1 |
| F4 | S.m6_na_added | 0.00022 | [-0.00088, 0.00130] | FALSE | 0 |
| F5_B | B.m6_l_H | -0.00391 | [-0.01184, 0.00383] | FALSE | 0 |
| F5_Bp | Bp.m6_l_H | -0.00522 | [-0.01818, 0.00707] | FALSE | 0 |

## The decision table (the generated `results/decision-table.md`, verbatim below its title)


Outcome: **K2**. Positive control (S's lower 97.5 % limit of Phi above 1.20): passes.
Replicates: A 200, B 200, C 200, D 200, Av 200, Ag 200. B = 10000.
No candidate withheld.

## Subgroup (the ten deciding settings, vd; ratios of sums; limits two-sided 95 %)

| Arm | Phi (95 % limits) | Class | rho_l against O | rho_l0 (C0) | share of S's excess removed | PhiH |
|---|---|---|---|---|---|---|
| S (shipped) | 1.490 [1.430, 1.548] | PC passes | 2.600 [2.475, 2.727] |  |  | 1.549 [1.495, 1.601] |
| A | 0.840 [0.779, 0.904] | removes / C | 1.467 [1.361, 1.579] |  | 1.326 [1.191, 1.468] | 0.854 [0.797, 0.911] |
| B | 1.018 [0.993, 1.043] | removes / E | 1.777 [1.712, 1.844] |  | 0.963 [0.917, 1.016] | 1.051 [1.027, 1.075] |
| B' | 1.142 [1.111, 1.173] | removes / E | 1.993 [1.919, 2.071] |  | 0.711 [0.660, 0.761] | 1.176 [1.147, 1.206] |
| N (no filter) |  |  | 1.746 [1.699, 1.794] |  |  |  |
| C0 (complete data) |  |  | NA [NA, NA] | 1.220 [1.192, 1.247] |  |  |

## Near the threshold

| Arm | alpha, true level within 1 log2 above the threshold | ratio to C0 | power there | attenuation there |
|---|---|---|---|---|
| S (shipped) | 0.0505 [0.0492, 0.0519] | 1.012 [0.984, 1.041] | 0.1959 [0.1928, 0.1990] | -0.0447 [-0.0497, -0.0396] |
| A | 0.0466 [0.0453, 0.0480] | 0.934 [0.907, 0.963] | 0.1311 [0.1286, 0.1337] | -0.1440 [-0.1496, -0.1385] |
| B | 0.0526 [0.0513, 0.0540] | 1.054 [1.027, 1.083] | 0.2329 [0.2296, 0.2363] | -0.0301 [-0.0349, -0.0253] |
| B' | 0.0525 [0.0512, 0.0538] | 1.051 [1.023, 1.080] | 0.2309 [0.2275, 0.2342] | -0.0310 [-0.0357, -0.0261] |
| N (no filter) | 0.0537 [0.0527, 0.0548] | 1.076 [1.052, 1.102] | 0.2970 [0.2932, 0.3008] | -0.0281 [-0.0323, -0.0241] |
| C0 (complete data) | 0.0499 [0.0489, 0.0509] | 1 | 0.8318 [0.8286, 0.8351] | -0.0007 [-0.0031, 0.0017] |

## Whole list

| Arm | kept | tested | called | power | power change against S (S minus X) | false share | false share change (X minus S) |
|---|---|---|---|---|---|---|---|
| S (shipped) | 2687.3 | 2654.8 | 430.90 | 0.6816 [0.6800, 0.6831] |  | 0.0510 [0.0503, 0.0517] |  |
| A | 2687.3 | 2649.9 | 427.78 | 0.6766 [0.6751, 0.6782] | 0.0049 [0.0048, 0.0051] | 0.0510 [0.0503, 0.0517] | 0.0000 [-0.0001, 0.0001] |
| B | 2711.6 | 2673.3 | 433.26 | 0.6852 [0.6837, 0.6868] | -0.0037 [-0.0038, -0.0035] | 0.0511 [0.0503, 0.0518] | 0.0001 [0.0000, 0.0001] |
| B' | 2707.2 | 2669.2 | 432.97 | 0.6847 [0.6832, 0.6863] | -0.0032 [-0.0033, -0.0030] | 0.0511 [0.0504, 0.0518] | 0.0001 [0.0000, 0.0002] |
| N (no filter) | 3000.0 | 2811.6 | 445.52 | 0.7033 [0.7017, 0.7048] | -0.0217 [-0.0222, -0.0213] | 0.0529 [0.0522, 0.0536] | 0.0019 [0.0017, 0.0021] |
| C0 (complete data) | 3000.0 | 3000.0 | 525.35 | 0.8315 [0.8302, 0.8329] | -0.1500 [-0.1510, -0.1490] | 0.0503 [0.0495, 0.0511] | -0.0007 [-0.0011, -0.0002] |

## Fold change

| Arm | m6_l on the added set | m6_l on the kept set | attenuation, whole list | attenuation, added set |
|---|---|---|---|---|
| S (shipped) | 0.1169 [0.0991, 0.1347] | 0.0502 [0.0320, 0.0679] | -0.0074 [-0.0084, -0.0065] | -0.0192 [-0.0215, -0.0171] |
| A | -0.0008 [-0.0161, 0.0147] | -0.0067 [-0.0218, 0.0084] | -0.0101 [-0.0111, -0.0092] | -0.0325 [-0.0348, -0.0303] |
| B | 0.0232 [0.0124, 0.0340] | -0.0029 [-0.0141, 0.0083] | -0.0072 [-0.0081, -0.0063] | -0.0170 [-0.0193, -0.0149] |
| B' | 0.0317 [0.0187, 0.0447] | -0.0012 [-0.0147, 0.0122] | -0.0072 [-0.0081, -0.0062] | -0.0169 [-0.0191, -0.0147] |
| N (no filter) | 0.0219 [0.0125, 0.0315] | 0.0014 [-0.0084, 0.0111] | -0.0079 [-0.0088, -0.0069] | -0.0180 [-0.0201, -0.0160] |
| C0 (complete data) | NA [NA, NA] | NA [NA, NA] | 0.0003 [-0.0006, 0.0012] | NA [NA, NA] |

## What changes in users' results

| Arm | share of S's calls also called | new calls per data set | proteins gained per data set | proteins lost per data set | kept but not tested in the contrast |
|---|---|---|---|---|---|
| A | 0.9913 [0.9911, 0.9916] | 0.62 [0.57, 0.66] | 0.0 [0.0, 0.0] | 0.0 [0.0, 0.0] | 37.5 [37.2, 37.7] |
| B | 0.9987 [0.9986, 0.9988] | 2.90 [2.82, 2.99] | 24.3 [24.1, 24.4] | 0.0 [0.0, 0.0] | 38.3 [38.0, 38.6] |
| B' | 0.9989 [0.9988, 0.9990] | 2.55 [2.47, 2.64] | 19.9 [19.7, 20.1] | 0.0 [0.0, 0.0] | 38.0 [37.7, 38.3] |
| N (no filter) | 0.9965 [0.9963, 0.9967] | 16.13 [15.86, 16.41] | 312.7 [312.2, 313.1] | 0.0 [0.0, 0.0] | 188.4 [187.9, 188.9] |

## Intensity-dependent variance (Av, reported only, no label)

| stat | vd5 | vd10 | vd25 |
|---|---|---|---|
| S.rho_l | 1.7659 [1.5281, 2.0180] | 1.4994 [1.1825, 1.8387] | 2.3437 [1.3270, 3.4847] |
| S.pow | 0.7001 [0.6968, 0.7034] | 0.6614 [0.6580, 0.6646] | 0.5413 [0.5380, 0.5447] |
| S.fdr_whole | 0.0491 [0.0474, 0.0507] | 0.0504 [0.0486, 0.0522] | 0.0495 [0.0477, 0.0513] |
| S.alpha_band | 0.0348 [0.0306, 0.0391] | 0.0338 [0.0297, 0.0381] | 0.0432 [0.0351, 0.0519] |
| A.rho_l | 1.8816 [1.6347, 2.1303] | 1.9663 [1.6153, 2.3325] | 3.0321 [1.9366, 4.1523] |
| A.pow | 0.7011 [0.6978, 0.7044] | 0.6623 [0.6590, 0.6656] | 0.5416 [0.5382, 0.5449] |
| A.fdr_whole | 0.0497 [0.0480, 0.0513] | 0.0507 [0.0490, 0.0525] | 0.0496 [0.0478, 0.0514] |
| A.alpha_band | 0.0382 [0.0341, 0.0426] | 0.0371 [0.0329, 0.0415] | 0.0432 [0.0351, 0.0518] |
| B.rho_l | 1.8100 [1.6806, 1.9433] | 1.8943 [1.7124, 2.0846] | 2.6195 [2.1562, 3.1018] |
| B.pow | 0.7006 [0.6972, 0.7040] | 0.6623 [0.6590, 0.6656] | 0.5412 [0.5378, 0.5445] |
| B.fdr_whole | 0.0495 [0.0478, 0.0511] | 0.0508 [0.0490, 0.0526] | 0.0496 [0.0478, 0.0514] |
| B.alpha_band | 0.0505 [0.0462, 0.0547] | 0.0485 [0.0436, 0.0535] | 0.0566 [0.0484, 0.0651] |
| Bp.rho_l | 2.0381 [1.8823, 2.2017] | 2.0162 [1.8081, 2.2360] | 2.7496 [2.2560, 3.2613] |
| Bp.pow | 0.7007 [0.6974, 0.7042] | 0.6623 [0.6590, 0.6656] | 0.5412 [0.5378, 0.5446] |
| Bp.fdr_whole | 0.0494 [0.0478, 0.0510] | 0.0508 [0.0490, 0.0526] | 0.0497 [0.0479, 0.0515] |
| Bp.alpha_band | 0.0490 [0.0447, 0.0534] | 0.0478 [0.0429, 0.0528] | 0.0568 [0.0486, 0.0653] |
| N.rho_l | 1.7290 [1.6345, 1.8251] | 1.9512 [1.8283, 2.0785] | 2.6143 [2.4506, 2.7865] |
| N.alpha_band | 0.0537 [0.0503, 0.0571] | 0.0546 [0.0513, 0.0579] | 0.0628 [0.0595, 0.0661] |
| C0.alpha_band | 0.0499 [0.0469, 0.0529] | 0.0498 [0.0469, 0.0528] | 0.0507 [0.0487, 0.0528] |
| Phi.S | 1.0213 [0.8985, 1.1484] | 0.7684 [0.6105, 0.9356] | 0.8965 [0.5102, 1.3275] |
| Phi.A | 1.0882 [0.9440, 1.2364] | 1.0077 [0.8316, 1.1972] | 1.1598 [0.7463, 1.5808] |
| Phi.B | 1.0468 [0.9958, 1.0981] | 0.9708 [0.8938, 1.0497] | 1.0020 [0.8373, 1.1673] |
| Phi.Bp | 1.1788 [1.1118, 1.2466] | 1.0333 [0.9449, 1.1222] | 1.0517 [0.8785, 1.2277] |

## Groups required (Ag, reported only, no label): 2 groups against 1 under B and B' (VAL-204 (2), (7); VAL-315)

| stat | vd5 | vd10 | vd25 |
|---|---|---|---|
| S.pow | 0.8309 [0.8292, 0.8326] | 0.7449 [0.7429, 0.7470] | 0.5397 [0.5370, 0.5422] |
| S.fdr_whole | 0.0502 [0.0487, 0.0517] | 0.0506 [0.0490, 0.0522] | 0.0502 [0.0484, 0.0520] |
| S.kept_per_r | 2806.9350 [2805.8099, 2808.0351] | 2627.1150 [2625.6399, 2628.6700] | 2037.7150 [2035.2397, 2040.1401] |
| S.pow_down | 0.6925 [0.6894, 0.6956] | 0.5570 [0.5534, 0.5606] | 0.2954 [0.2917, 0.2990] |
| S.powbh | 0.8305 [0.8289, 0.8322] | 0.7446 [0.7426, 0.7467] | 0.5394 [0.5367, 0.5420] |
| S.powbh_down | 0.6922 [0.6891, 0.6954] | 0.5568 [0.5532, 0.5604] | 0.2952 [0.2916, 0.2988] |
| S.onoff_removed | 1.0000 [1.0000, 1.0000] | 0.9748 [0.9726, 0.9770] | 0.9146 [0.9121, 0.9171] |
| S.onoff_powbh | 0.0000 [0.0000, 0.0000] | 0.0251 [0.0229, 0.0273] | 0.0853 [0.0828, 0.0878] |
| A.pow | 0.8309 [0.8292, 0.8326] | 0.7450 [0.7430, 0.7470] | 0.5397 [0.5370, 0.5422] |
| A.fdr_whole | 0.0504 [0.0489, 0.0519] | 0.0508 [0.0492, 0.0524] | 0.0502 [0.0484, 0.0520] |
| A.kept_per_r | 2806.9350 [2805.8099, 2808.0351] | 2627.1150 [2625.6399, 2628.6700] | 2037.7150 [2035.2397, 2040.1401] |
| A.pow_down | 0.6925 [0.6894, 0.6956] | 0.5570 [0.5534, 0.5606] | 0.2954 [0.2917, 0.2990] |
| A.powbh | 0.8305 [0.8289, 0.8322] | 0.7446 [0.7426, 0.7467] | 0.5394 [0.5367, 0.5420] |
| A.powbh_down | 0.6922 [0.6891, 0.6953] | 0.5568 [0.5532, 0.5604] | 0.2952 [0.2916, 0.2988] |
| A.onoff_removed | 1.0000 [1.0000, 1.0000] | 0.9748 [0.9726, 0.9770] | 0.9146 [0.9121, 0.9171] |
| A.onoff_powbh | 0.0000 [0.0000, 0.0000] | 0.0251 [0.0229, 0.0273] | 0.0853 [0.0828, 0.0878] |
| B.pow | 0.8733 [0.8718, 0.8749] | 0.7724 [0.7705, 0.7745] | 0.5481 [0.5455, 0.5507] |
| B.fdr_whole | 0.0502 [0.0487, 0.0516] | 0.0506 [0.0490, 0.0522] | 0.0499 [0.0482, 0.0517] |
| B.kept_per_r | 2844.5200 [2843.3050, 2845.7100] | 2650.8900 [2649.3500, 2652.5101] | 2044.4700 [2041.9800, 2046.8951] |
| B.pow_down | 0.7705 [0.7675, 0.7736] | 0.6059 [0.6023, 0.6095] | 0.3090 [0.3054, 0.3126] |
| B.powbh | 0.8729 [0.8713, 0.8745] | 0.7722 [0.7702, 0.7742] | 0.5478 [0.5451, 0.5503] |
| B.powbh_down | 0.7701 [0.7671, 0.7732] | 0.6057 [0.6021, 0.6092] | 0.3088 [0.3051, 0.3124] |
| B.onoff_removed | 0.8230 [0.8153, 0.8308] | 0.8493 [0.8443, 0.8542] | 0.8941 [0.8914, 0.8967] |
| B.onoff_powbh | 0.1769 [0.1691, 0.1846] | 0.1505 [0.1457, 0.1555] | 0.1058 [0.1031, 0.1084] |
| Bp.pow | 0.8686 [0.8670, 0.8702] | 0.7708 [0.7688, 0.7728] | 0.5478 [0.5452, 0.5504] |
| Bp.fdr_whole | 0.0501 [0.0486, 0.0515] | 0.0506 [0.0490, 0.0522] | 0.0499 [0.0481, 0.0517] |
| Bp.kept_per_r | 2836.5250 [2835.3350, 2837.6900] | 2647.7150 [2646.1900, 2649.3151] | 2044.0000 [2041.5000, 2046.4401] |
| Bp.pow_down | 0.7618 [0.7588, 0.7649] | 0.6029 [0.5993, 0.6065] | 0.3086 [0.3049, 0.3122] |
| Bp.powbh | 0.8682 [0.8665, 0.8698] | 0.7705 [0.7685, 0.7725] | 0.5475 [0.5449, 0.5501] |
| Bp.powbh_down | 0.7614 [0.7584, 0.7645] | 0.6027 [0.5991, 0.6062] | 0.3084 [0.3047, 0.3120] |
| Bp.onoff_removed | 0.8484 [0.8414, 0.8554] | 0.8589 [0.8540, 0.8637] | 0.8947 [0.8920, 0.8974] |
| Bp.onoff_powbh | 0.1515 [0.1445, 0.1585] | 0.1409 [0.1362, 0.1458] | 0.1051 [0.1025, 0.1078] |
| B1.pow | 0.9474 [0.9464, 0.9485] | 0.8883 [0.8867, 0.8899] | 0.7206 [0.7184, 0.7229] |
| B1.fdr_whole | 0.0502 [0.0487, 0.0516] | 0.0508 [0.0493, 0.0522] | 0.0498 [0.0483, 0.0513] |
| B1.kept_per_r | 2976.6250 [2976.0050, 2977.2400] | 2920.5550 [2919.5550, 2921.5700] | 2633.1800 [2631.4100, 2634.9900] |
| B1.pow_down | 0.9025 [0.9005, 0.9045] | 0.7969 [0.7941, 0.7998] | 0.5303 [0.5265, 0.5342] |
| B1.powbh | 0.9468 [0.9457, 0.9479] | 0.8877 [0.8861, 0.8893] | 0.7199 [0.7177, 0.7222] |
| B1.powbh_down | 0.9018 [0.8997, 0.9038] | 0.7962 [0.7933, 0.7991] | 0.5297 [0.5259, 0.5337] |
| B1.onoff_removed | 0.1107 [0.1049, 0.1165] | 0.1538 [0.1495, 0.1580] | 0.2648 [0.2614, 0.2682] |
| B1.onoff_powbh | 0.5252 [0.5158, 0.5344] | 0.4655 [0.4591, 0.4716] | 0.3580 [0.3538, 0.3624] |
| Bp1.pow | 0.9423 [0.9412, 0.9434] | 0.8864 [0.8848, 0.8881] | 0.7203 [0.7181, 0.7226] |
| Bp1.fdr_whole | 0.0502 [0.0487, 0.0516] | 0.0509 [0.0494, 0.0523] | 0.0498 [0.0483, 0.0513] |
| Bp1.kept_per_r | 2961.3300 [2960.6550, 2961.9900] | 2911.2950 [2910.2800, 2912.3300] | 2630.6000 [2628.8300, 2632.4000] |
| Bp1.pow_down | 0.8931 [0.8910, 0.8952] | 0.7934 [0.7905, 0.7964] | 0.5298 [0.5260, 0.5338] |
| Bp1.powbh | 0.9417 [0.9406, 0.9428] | 0.8858 [0.8842, 0.8875] | 0.7196 [0.7174, 0.7219] |
| Bp1.powbh_down | 0.8924 [0.8903, 0.8946] | 0.7928 [0.7899, 0.7957] | 0.5292 [0.5254, 0.5332] |
| Bp1.onoff_removed | 0.1439 [0.1372, 0.1507] | 0.1654 [0.1609, 0.1698] | 0.2657 [0.2622, 0.2691] |
| Bp1.onoff_powbh | 0.4954 [0.4857, 0.5052] | 0.4548 [0.4484, 0.4610] | 0.3572 [0.3530, 0.3616] |
| d1.B.pow | 0.0741 [0.0727, 0.0755] | 0.1159 [0.1142, 0.1176] | 0.1725 [0.1703, 0.1748] |
| d1.Bp.pow | 0.0737 [0.0723, 0.0751] | 0.1157 [0.1140, 0.1173] | 0.1725 [0.1703, 0.1748] |
| d1.B.fdr_whole | -0.0000 [-0.0004, 0.0004] | 0.0002 [-0.0003, 0.0007] | -0.0001 [-0.0010, 0.0007] |
| d1.Bp.fdr_whole | 0.0001 [-0.0003, 0.0005] | 0.0003 [-0.0003, 0.0008] | -0.0001 [-0.0009, 0.0008] |
| d1.B.powbh_down | 0.1317 [0.1291, 0.1344] | 0.1905 [0.1875, 0.1936] | 0.2210 [0.2174, 0.2246] |
| d1.Bp.powbh_down | 0.1310 [0.1284, 0.1336] | 0.1901 [0.1871, 0.1932] | 0.2209 [0.2173, 0.2245] |
| d1.B.onoff_powbh | 0.3482 [0.3384, 0.3579] | 0.3150 [0.3090, 0.3210] | 0.2522 [0.2484, 0.2561] |
| d1.Bp.onoff_powbh | 0.3439 [0.3342, 0.3535] | 0.3138 [0.3079, 0.3198] | 0.2521 [0.2483, 0.2560] |

## Fixed rows

| Column | Content |
|---|---|
| Routes affected | DIA without limpa, LFQ, TMT; the legacy limma test (with and without TREAT) and the declared designs; not DIA with limpa; strict mode: A changes nothing (strict mode keeps only proteins with no value at or below the threshold); B and B' change it to "no NA" |
| What users are told and see | A: the lowest values of the table (at most `percentile / 100` of them) are left out of the test; the test treats them as absent although they are known to be low (left-censored at the threshold, not missing at random), so fold changes of changed proteins near the threshold are biased toward 0; the share left out of the test is the percentile itself (1 % at the default, up to 10 % at the tab's maximum), and its costs near the threshold grow with it; only the default is measured. B: the percentile control no longer does anything at this step. B': the percentile becomes a label-blind low-intensity rule as in lipidomics and metabolomics, with its proportion (0.5); B' is the rule lipidomics and metabolomics ship after KIT-217, with proteomics' minimums. At the default percentile B' removes only proteins with at least half of all samples below the table's 1st percentile, so it keeps nearly B's proteins. The choice between B and B' is whether the percentile control keeps a meaning, not a measured difference. S: today's help text |
| Flag column (SD8; decided with the direction) | whether results tables mark proteins kept with a value at or below the threshold; how many per data set: the kept counts of the first table of 'Whole list' and the added-set counts per arm |
| Code change it needs | A: NA at the DA input on both DA entry points, a record field, a replay marker; B: the presence rule in one helper, the control, record and texts, a replay marker; B': B plus the share rule and its proportion parameter |
| Measured where | the simulated design only: one generator whose variance does not depend on intensity (in real label-free and DIA data the lowest intensities are the noisiest, which changes both what A discards and how large the subgroup's excess is under B); TREAT with the trend and robust settings of the bound commit; equal groups of 3 or 4; no normalisation or RUV between the filter and the test; the default percentile and minimums; strict mode not run. Which direction is better on real data is not measured. One reported cell (Av) with intensity-dependent variance shows the direction of that effect in one design |


## Groups required (VAL-204 (2), (7); VAL-315): the ticket's 6.9 text, verbatim, beside Ag's figures above (DR-04)

Groups required (VAL-204 (2), (7); VAL-315) | the cell Ag (18.3): B and B′ at 2 groups required against 1, paired (power on strong decreases by the shipped rule and by BH, the on/off share removed, whole-list power and false share), with KIT-217's cell Le beside (B′ at (2, 2) and (2, 1) correspond to its Aobs_g2 and Aobs_any: BH power 0.712 against 0.819 at vd25). Not comparable in level: Le is lipid (500 features in two assays, median SD 0.5), 2 x 6, the 10th percentile, R = 100; Ag is proteomics (3,000 features, median SD 0.27), 2 x 3, the 1st percentile, R = 200; only the direction and size of each cell's paired difference are read side by side. In 2 x 3 a protein kept at 1 group with no value in a compared group is still not tested there (KIT-204), so 1 group helps only proteins left with one value in the lower group. Not run: S's and A's filter (values at or below τ counted missing) at 1 group.

Ag's figures at vd25 (BH power on the 300 strong decreases): B 0.3088 at 2 groups, 0.5297 at 1 group (paired change 0.2210 [0.2174, 0.2246]); B′ 0.3084 and 0.5292 (0.2209 [0.2173, 0.2245]). The same at vd5 and vd10: `d1.B.powbh_down` and `d1.Bp.powbh_down` in the table above.

## Per-setting statement beside every pooled label (RA-01)

Pooled over the ten deciding settings (the prespecified analysis), B and B′ are within the 1.20 margin (1.02; 1.14); in every individual setting the ratio is above 1 (B 1.04 to 1.38, B′ 1.11 to 1.48), with an excess remaining (one-sided 0.05/72) in 3 settings for B and 7 for B′. B′'s pooled limits exclude 1 ("non-zero, within the margin"). The pooled values are lower than every setting's because the arms' subgroups are distributed differently across the settings (shares in `outcome.md`). Per-setting point estimates of Φ (S, A, B, B′): A vd5 1.31, 1.35, 1.11, 1.39; A vd10 0.90, 1.38, 1.18, 1.46; A vd25 1.36, 1.58, 1.28, 1.48; B vd5 1.62, 0.90, 1.05, 1.25; B vd10 1.91, 0.94, 1.17, 1.29; B vd25 2.20, 0.91, 1.38, 1.40; C vd5 1.44, 0.86, 1.04, 1.12; C vd10 1.68, 0.70, 1.07, 1.11; C vd25 2.10, 0.96, 1.35, 1.35; D vd25 1.99, 0.74, 1.18, 1.19. A is above 1 in the 2 x 3 cell, where S's own excess is not resolved. A's costs: whole-list power change 0.0049 (limits around the 0.005 margin, not determined); attenuation near the threshold -0.144 against -0.045 shipped (about three times as large). Never write "B leaves no excess beyond the unfiltered analysis".

The 1.22 figure's comparator: on complete data with no filter, proteins holding a value at or below the threshold against proteins holding no such value (not against proteins kept under both rules).

## Draft Methods sentence (the auditor's proposal RA5; a draft for the owner, VAL-41; the clause naming the direction is added after the choice)

In simulation (3,000 proteins; two to five groups of three or three groups of four; 5, 10 and 25 % of values missing with probability decreasing with intensity; 200 data sets per setting, ten settings pooled as prespecified), unchanged proteins that the filter newly kept and that held a value at or below its intensity threshold reached p < 0.05 2.6 times as often as proteins kept under both the earlier and the current rule, against 1.75 times for the same comparison without a filter; counting such values as missing in the filter while the test used them thus multiplied the excess by 1.49 (one-sided 97.5 % lower limit 1.43). Counting them as present kept this factor within a prespecified margin of 1.20 (1.02; 1.14 with an added group-blind low-intensity rule) without loss of power or change in the false discovery proportion, although in each individual setting it remained above 1 (1.04 to 1.38; 1.11 to 1.48); leaving them out of the test as well lowered it to 0.84, at a whole-list power loss of 0.005 and with three-fold stronger attenuation of fold changes near the threshold. Even in complete data without a filter, proteins holding such a value reached p < 0.05 1.22 times as often as other proteins. Simulated data only: one generator whose variance does not depend on intensity, equal group sizes, no normalisation or RUV between filter and test.

## Proposed corrections of KIT-197's merged texts (PROPOSALS only: NEWS.md, `kit-197/` and the other ledger rows are not edited here; the coordinator applies them at landing after the owner's wording pass)

KIT-197's NEWS clause (KIT-187's bullet), ledger VAL-149 and the Methods draft in `kit-197/results.md` attribute the 2.5-fold figure to the filter and test inconsistency alone. The figure stands (2.60 here; KIT-197 2.54); the attribution changes. Pooled over ten settings; in the 2 x 3 design not resolved. The split: against the unfiltered analysis on the same masked data Φ_S = 1.49 [1.43, 1.55]; the same ratio without a filter is 1.75; on complete data, proteins holding such a value reach p < 0.05 1.22 times as often as proteins holding no such value. Proposed replacement for the attribution: "unchanged proteins holding a value at or below the threshold reach p < 0.05 about 1.5 times as often under the shipped filter as the corresponding newly kept proteins in the same data analysed without a filter (pooled over ten settings; not resolved in the 2 x 3 design); the ratio against proteins kept under both rules is about 2.6, and about 1.75 without a filter; on complete data such proteins reach p < 0.05 about 1.2 times as often as proteins holding no such value".
