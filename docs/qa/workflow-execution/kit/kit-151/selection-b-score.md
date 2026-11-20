# KIT-151 selection evidence, run b: score

Run counted; seeds 1..30 (30); cells ok 1020 of 1020; failed 0; missing 0; cells stopped at the round limit (10000): 0. CPU models: {'Intel(R) Xeon(R) Platinum 8370C CPU @ 2.80GHz': 2, 'AMD EPYC 7763 64-Core Processor': 11, 'AMD EPYC 9V74 80-Core Processor': 10, 'INTEL(R) XEON(R) PLATINUM 8573C': 5, 'AMD EPYC 9V45 96-Core Processor': 2}.

## prev: matched |r|, median / 10th percentile / worst seed

| arm | F1 | F2 | F3 | worst of three | iterations median (max); at limit |
|---|---|---|---|---|---|
| none | 0.987 / 0.873 / 0.542 | 0.405 / 0.186 / 0.081 | 0.963 / 0.623 / 0.583 | 0.405 / 0.186 / 0.081 | 371 (736); 0 |
| scale_views | 0.990 / 0.874 / 0.542 | 0.378 / 0.148 / 0.078 | 0.974 / 0.623 / 0.586 | 0.378 / 0.148 / 0.078 | 528 (836); 0 |
| raw250 | 0.994 / 0.993 / 0.990 | 0.958 / 0.919 / 0.887 | 0.988 / 0.977 / 0.945 | 0.956 / 0.919 / 0.887 | 908 (1746); 0 |
| raw1000 | 0.997 / 0.995 / 0.989 | 0.950 / 0.897 / 0.527 | 0.994 / 0.967 / 0.940 | 0.946 / 0.897 / 0.527 | 776 (1911); 0 |
| raw2000 | 0.997 / 0.990 / 0.940 | 0.939 / 0.771 / 0.232 | 0.994 / 0.942 / 0.754 | 0.938 / 0.723 / 0.232 | 906 (1306); 0 |
| random1000 | 0.979 / 0.969 / 0.941 | 0.953 / 0.905 / 0.309 | 0.938 / 0.884 / 0.801 | 0.918 / 0.873 / 0.309 | 788 (1241); 0 |
| trend1000 | 0.997 / 0.995 / 0.988 | 0.949 / 0.898 / 0.838 | 0.994 / 0.967 / 0.946 | 0.948 / 0.898 / 0.838 | 736 (1276); 0 |
| trendratio1000 | 0.997 / 0.995 / 0.989 | 0.949 / 0.904 / 0.835 | 0.994 / 0.967 / 0.939 | 0.947 / 0.904 / 0.835 | 776 (1891); 0 |
| raw1000_z | 0.998 / 0.996 / 0.990 | 0.960 / 0.918 / 0.852 | 0.994 / 0.969 / 0.949 | 0.958 / 0.918 / 0.852 | 751 (2011); 0 |

### prev: paired differences (b minus a), worst-of-three per seed: mean [95% CI] reading

| group | a | b | mean [CI] | b better / worse | worst seed diff | reading |
|---|---|---|---|---|---|---|
| selecting at all | scale_views | raw250 | +0.509 [+0.418, +0.595] | 30 / 0 | +0.011 | better |
| selecting at all | scale_views | raw1000 | +0.486 [+0.397, +0.571] | 30 / 0 | +0.003 | better |
| selecting at all | scale_views | raw2000 | +0.424 [+0.312, +0.524] | 28 / 2 | -0.506 | better |
| selecting at all | scale_views | random1000 | +0.456 [+0.365, +0.542] | 29 / 1 | -0.066 | better |
| selecting at all | scale_views | trend1000 | +0.498 [+0.409, +0.582] | 30 / 0 | +0.008 | better |
| selecting at all | scale_views | trendratio1000 | +0.498 [+0.408, +0.581] | 30 / 0 | +0.007 | better |
| statistic | raw1000 | random1000 | -0.030 [-0.079, +0.015] | 11 / 19 | -0.544 | inconclusive |
| statistic | raw1000 | trend1000 | +0.012 [-0.004, +0.040] | 16 / 14 | -0.053 | inconclusive |
| statistic | raw1000 | trendratio1000 | +0.012 [-0.003, +0.039] | 16 / 14 | -0.026 | inconclusive |
| statistic | random1000 | trend1000 | +0.042 [+0.013, +0.086] | 20 / 10 | -0.042 | better |
| statistic | random1000 | trendratio1000 | +0.042 [+0.013, +0.084] | 20 / 10 | -0.041 | better |
| statistic | trend1000 | trendratio1000 | +0.000 [-0.002, +0.003] | 17 / 13 | -0.015 | no practical difference |
| number | raw1000 | raw250 | +0.023 [+0.005, +0.054] | 21 / 9 | -0.010 | better |
| number | raw1000 | raw2000 | -0.061 [-0.134, +0.000] | 9 / 21 | -0.622 | inconclusive |
| z-score | raw1000 | raw1000_z | +0.024 [+0.008, +0.050] | 25 / 5 | -0.008 | better |
| scale_views | none | scale_views | -0.024 [-0.080, +0.008] | 15 / 15 | -0.760 | inconclusive |

## trend: matched |r|, median / 10th percentile / worst seed

| arm | F1 | F2 | F3 | worst of three | iterations median (max); at limit |
|---|---|---|---|---|---|
| none | 0.998 / 0.879 / 0.854 | 0.657 / 0.248 / 0.112 | 0.989 / 0.639 / 0.555 | 0.576 / 0.248 / 0.112 | 296 (671); 0 |
| scale_views | 0.998 / 0.872 / 0.845 | 0.322 / 0.180 / 0.124 | 0.991 / 0.626 / 0.550 | 0.322 / 0.180 / 0.124 | 664 (1166); 0 |
| raw250 | 0.924 / 0.857 / 0.669 | 0.971 / 0.947 / 0.890 | 0.548 / 0.298 / 0.205 | 0.548 / 0.298 / 0.205 | 941 (2781); 0 |
| raw1000 | 0.989 / 0.982 / 0.980 | 0.969 / 0.918 / 0.156 | 0.922 / 0.849 / 0.821 | 0.899 / 0.829 / 0.156 | 781 (1406); 0 |
| raw2000 | 0.994 / 0.990 / 0.986 | 0.952 / 0.276 / 0.184 | 0.972 / 0.938 / 0.919 | 0.935 / 0.276 / 0.184 | 996 (1601); 0 |
| random1000 | 0.990 / 0.981 / 0.973 | 0.969 / 0.386 / 0.240 | 0.967 / 0.931 / 0.865 | 0.945 / 0.386 / 0.240 | 861 (1206); 0 |
| trend1000 | 0.996 / 0.993 / 0.991 | 0.965 / 0.822 / 0.335 | 0.987 / 0.969 / 0.949 | 0.964 / 0.822 / 0.335 | 936 (1491); 0 |
| trendratio1000 | 0.999 / 0.997 / 0.997 | 0.962 / 0.298 / 0.211 | 0.995 / 0.977 / 0.941 | 0.960 / 0.298 / 0.211 | 1014 (1701); 0 |
| raw1000_z | 0.989 / 0.982 / 0.978 | 0.969 / 0.943 / 0.870 | 0.923 / 0.842 / 0.681 | 0.907 / 0.842 / 0.681 | 786 (1076); 0 |

### trend: paired differences (b minus a), worst-of-three per seed: mean [95% CI] reading

| group | a | b | mean [CI] | b better / worse | worst seed diff | reading |
|---|---|---|---|---|---|---|
| selecting at all | scale_views | raw250 | +0.080 [-0.042, +0.196] | 20 / 10 | -0.747 | inconclusive |
| selecting at all | scale_views | raw1000 | +0.418 [+0.299, +0.528] | 23 / 7 | -0.309 | better |
| selecting at all | scale_views | raw2000 | +0.294 [+0.148, +0.433] | 23 / 7 | -0.556 | better |
| selecting at all | scale_views | random1000 | +0.391 [+0.280, +0.499] | 27 / 3 | -0.310 | better |
| selecting at all | scale_views | trend1000 | +0.462 [+0.349, +0.569] | 28 / 2 | -0.264 | better |
| selecting at all | scale_views | trendratio1000 | +0.354 [+0.221, +0.481] | 24 / 6 | -0.439 | better |
| statistic | raw1000 | random1000 | -0.027 [-0.142, +0.089] | 22 / 8 | -0.704 | inconclusive |
| statistic | raw1000 | trend1000 | +0.044 [-0.035, +0.122] | 26 / 4 | -0.583 | inconclusive |
| statistic | raw1000 | trendratio1000 | -0.065 [-0.178, +0.047] | 21 / 9 | -0.713 | inconclusive |
| statistic | random1000 | trend1000 | +0.071 [-0.007, +0.154] | 22 / 8 | -0.594 | inconclusive |
| statistic | random1000 | trendratio1000 | -0.037 [-0.155, +0.077] | 18 / 12 | -0.768 | inconclusive |
| statistic | trend1000 | trendratio1000 | -0.109 [-0.202, -0.029] | 15 / 15 | -0.764 | worse |
| number | raw1000 | raw250 | -0.338 [-0.423, -0.239] | 1 / 29 | -0.740 | worse |
| number | raw1000 | raw2000 | -0.124 [-0.234, -0.025] | 19 / 11 | -0.760 | worse |
| z-score | raw1000 | raw1000_z | +0.040 [-0.014, +0.114] | 13 / 17 | -0.213 | inconclusive |
| scale_views | none | scale_views | -0.149 [-0.232, -0.077] | 7 / 23 | -0.682 | worse |

## full: matched |r|, median / 10th percentile / worst seed

| arm | F1 | F2 | F3 | worst of three | iterations median (max); at limit |
|---|---|---|---|---|---|
| none | 0.995 / 0.881 / 0.621 | 0.347 / 0.211 / 0.166 | 0.989 / 0.646 / 0.565 | 0.347 / 0.211 / 0.166 | 504 (836); 0 |
| scale_views | 0.995 / 0.871 / 0.624 | 0.313 / 0.193 / 0.104 | 0.988 / 0.646 / 0.562 | 0.313 / 0.193 / 0.104 | 1024 (1706); 0 |
| raw250 | 0.938 / 0.837 / 0.792 | 0.931 / 0.342 / 0.159 | 0.522 / 0.248 / 0.137 | 0.444 / 0.238 / 0.137 | 1644 (3486); 0 |
| raw1000 | 0.991 / 0.987 / 0.984 | 0.486 / 0.226 / 0.177 | 0.971 / 0.940 / 0.922 | 0.486 / 0.226 / 0.177 | 1191 (1791); 0 |
| raw2000 | 0.995 / 0.991 / 0.987 | 0.480 / 0.201 / 0.118 | 0.987 / 0.964 / 0.936 | 0.480 / 0.201 / 0.118 | 1291 (2061); 0 |
| random1000 | 0.985 / 0.976 / 0.967 | 0.702 / 0.221 / 0.124 | 0.953 / 0.886 / 0.773 | 0.702 / 0.221 / 0.124 | 1088 (1941); 0 |
| trend1000 | 0.995 / 0.992 / 0.990 | 0.721 / 0.218 / 0.066 | 0.986 / 0.970 / 0.933 | 0.721 / 0.218 / 0.066 | 1368 (2191); 0 |
| trendratio1000 | 0.998 / 0.997 / 0.996 | 0.618 / 0.306 / 0.148 | 0.994 / 0.976 / 0.949 | 0.618 / 0.306 / 0.148 | 1221 (1671); 0 |
| raw1000_z | 0.990 / 0.987 / 0.985 | 0.959 / 0.919 / 0.399 | 0.969 / 0.948 / 0.942 | 0.951 / 0.919 / 0.399 | 1276 (1966); 0 |
| raw250_obs50 | 0.975 / 0.957 / 0.938 | 0.933 / 0.486 / 0.303 | 0.939 / 0.855 / 0.851 | 0.884 / 0.486 / 0.303 | 1228 (2046); 0 |
| raw1000_obs50 | 0.994 / 0.991 / 0.990 | 0.916 / 0.161 / 0.113 | 0.986 / 0.972 / 0.966 | 0.916 / 0.161 / 0.113 | 1131 (2091); 0 |
| raw2000_obs50 | 0.996 / 0.993 / 0.990 | 0.493 / 0.276 / 0.166 | 0.991 / 0.968 / 0.938 | 0.493 / 0.276 / 0.166 | 1101 (1611); 0 |
| random1000_obs50 | 0.985 / 0.980 / 0.974 | 0.674 / 0.172 / 0.067 | 0.958 / 0.894 / 0.415 | 0.603 / 0.172 / 0.067 | 954 (1781); 0 |
| trend1000_obs50 | 0.996 / 0.995 / 0.994 | 0.510 / 0.268 / 0.084 | 0.991 / 0.975 / 0.964 | 0.510 / 0.268 / 0.084 | 1016 (1926); 0 |
| trendratio1000_obs50 | 0.998 / 0.997 / 0.996 | 0.540 / 0.280 / 0.220 | 0.994 / 0.977 / 0.936 | 0.540 / 0.280 / 0.220 | 1008 (1716); 0 |
| raw1000_z_obs50 | 0.994 / 0.992 / 0.991 | 0.961 / 0.867 / 0.312 | 0.985 / 0.973 / 0.954 | 0.961 / 0.867 / 0.312 | 1118 (2241); 0 |

### full: paired differences (b minus a), worst-of-three per seed: mean [95% CI] reading

| group | a | b | mean [CI] | b better / worse | worst seed diff | reading |
|---|---|---|---|---|---|---|
| selecting at all | scale_views | raw250 | +0.056 [-0.040, +0.149] | 18 / 12 | -0.586 | inconclusive |
| selecting at all | scale_views | raw1000 | +0.204 [+0.082, +0.325] | 20 / 10 | -0.556 | better |
| selecting at all | scale_views | raw2000 | +0.166 [+0.018, +0.310] | 22 / 8 | -0.742 | better |
| selecting at all | scale_views | random1000 | +0.265 [+0.130, +0.396] | 21 / 9 | -0.660 | better |
| selecting at all | scale_views | trend1000 | +0.286 [+0.169, +0.403] | 23 / 7 | -0.326 | better |
| selecting at all | scale_views | trendratio1000 | +0.278 [+0.162, +0.395] | 21 / 9 | -0.208 | better |
| statistic | raw1000 | random1000 | +0.061 [-0.052, +0.177] | 14 / 16 | -0.581 | inconclusive |
| statistic | raw1000 | trend1000 | +0.083 [-0.004, +0.179] | 16 / 14 | -0.272 | inconclusive |
| statistic | raw1000 | trendratio1000 | +0.074 [-0.033, +0.189] | 16 / 14 | -0.529 | inconclusive |
| statistic | random1000 | trend1000 | +0.022 [-0.094, +0.135] | 20 / 10 | -0.770 | inconclusive |
| statistic | random1000 | trendratio1000 | +0.013 [-0.086, +0.108] | 19 / 11 | -0.748 | inconclusive |
| statistic | trend1000 | trendratio1000 | -0.008 [-0.110, +0.103] | 8 / 22 | -0.575 | inconclusive |
| number | raw1000 | raw250 | -0.148 [-0.262, -0.035] | 9 / 21 | -0.749 | worse |
| number | raw1000 | raw2000 | -0.037 [-0.151, +0.074] | 14 / 16 | -0.695 | inconclusive |
| z-score | raw1000 | raw1000_z | +0.347 [+0.250, +0.444] | 28 / 2 | -0.002 | better |
| z-score | raw1000_obs50 | raw1000_z_obs50 | +0.234 [+0.138, +0.339] | 26 / 4 | -0.026 | better |
| scale_views | none | scale_views | -0.086 [-0.166, -0.017] | 12 / 18 | -0.638 | worse |
| eligibility | raw250 | raw250_obs50 | +0.358 [+0.267, +0.446] | 27 / 3 | -0.226 | better |
| eligibility | raw1000 | raw1000_obs50 | +0.097 [+0.000, +0.199] | 17 / 13 | -0.438 | better |
| eligibility | raw2000 | raw2000_obs50 | +0.007 [-0.075, +0.083] | 17 / 13 | -0.625 | inconclusive |
| eligibility | random1000 | random1000_obs50 | -0.033 [-0.127, +0.061] | 12 / 18 | -0.637 | inconclusive |
| eligibility | trend1000 | trend1000_obs50 | -0.066 [-0.154, +0.015] | 14 / 16 | -0.638 | inconclusive |
| eligibility | trendratio1000 | trendratio1000_obs50 | -0.034 [-0.126, +0.050] | 15 / 15 | -0.652 | inconclusive |
| eligibility | raw1000_z | raw1000_z_obs50 | -0.016 [-0.051, +0.005] | 15 / 15 | -0.474 | inconclusive |

## Tie to the follow-up (prev, seeds 1..10, against n20strong slow)

```
{
  "compared": 150,
  "identical_1e-9": 90,
  "within_1e-6": 97,
  "old_stopped_by_convergence": 120,
  "of_those_within_1e-6": 94,
  "max_abs_diff_old_converged": 0.7879589455525546
}
```

