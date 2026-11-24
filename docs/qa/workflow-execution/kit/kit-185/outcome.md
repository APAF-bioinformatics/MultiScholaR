# KIT-185 outcome: RUV-III-C and missing values, measured

Ticket revision 1.2; declaration `declaration.md` (commit 432f703b, sections 1 to 6; dated sections 7, 8 and 10 after
part A's results, section 9 before the part-B re-run; section 11, P2's addendum, committed alone at 6779d56e before P2
ran). Package code: `v0.9-staging` at ffdfdfb76a9eacdc681f2098f75ac574e2cf8af4; P2 at 6779d56e, which contains
`v0.9-staging` 6416dbd6 (KIT-187's count filter; merged into this branch at 070d34e1; between the two, the RUV-III-C S4
method, its helpers, the proteomics DA and the factor count are unchanged). Runs and core-hours: `runs.tsv`. Rows:
`rows/`; the jobs' checks `checks/`; timing and provenance `jobs/`; scored tables `results/` (final scoring run 45333019858).

## Answer

**Outcome O3** (the ticket's rule, 5.7) on all 14 deciding cells: the 456 declared classifications, 45 "serious", 53
"material", 312 "no material effect", 18 "not determined", 28 "refused by the route" (re-score run 45333019858). P2 ran
after KIT-187 (VAL-71, VAL-108; declaration section 11) and its 36 classifications are included; the other 420 are
unchanged from the first scoring (every non-P2 row of `results/cells.tsv` identical). P2 (two groups of three): at
random missingness B (0.088, 0.21, 0.24 at 0.5, 2, 10 %) and power (0.045, 0.13, 0.057) material from 0.5 %,
attributable to the zero-filling; the false discovery proportion material at random 10 % (+0.030; the comparator
serious, +0.215); value-dependent FDP serious from 10 % (+0.086, +0.21, +0.30 at 10, 25, 50 %; realised 4.6, 8.6,
13.2 %), the comparator higher (S − C −0.049 to −0.022, intervals excluding 0), so none of it is the zero-filling;
value-dependent power not determined at 2 %; safe rates 0 at random, 0.5 % value-dependent. **No remedy is decided here** (VAL-56 keeps it for its own decision).

Observed, in plain words:
1. **At random (MCAR) missingness, the zero-filling itself damages the correction.** From 0.5 % (P2), 2 % (P0, P0r,
   P1, L2, M2, P5) or 10 % (P3, P4) the shipped correction leaves 0.09 to 0.40 of the planted unwanted variance (B) that it removed
   without missing values, and the legacy test loses up to 0.45 of power beyond the oracle's own loss (L2 MCAR 10 %).
   The comparator with directions from complete features leaves none of it (C's Delta between -0.08 and +0.05; classified
   "no material effect" in every such case): 41 of the 98 material or serious deciding classifications are "attributable
   to the zero-filling" (35 at random, 6 value-dependent: P5 VD 25 and 50 % B, P0 and P0r VD 25 % power, M1b VD 2 % B, M2 VD 2 % power).
   The comparator is not free of cost at random missingness: its false discovery proportion is serious at 10 % in P0
   (+0.031), P0r (+0.035), P1 (+0.093), P2 (+0.215), L2 (+0.042) and M2 (+0.082) and at 2 % in P2 (+0.040), where the
   shipped correction's is not (−0.07 to −0.16: it leaves unwanted variance and calls less; P2 10 %: +0.030, material),
   and at P4 10 % it is higher than the shipped correction's (+0.128 against +0.077). Directions from complete features
   remove the zero-filling's damage to B and power, but not the rise in false discoveries that comes with correcting
   features with missing values (inferred: the degrees of freedom the correction uses are not counted by the legacy
   test, VAL-53 (1)).
2. **At value-dependent (left-censored) missingness, the zero-filling does little to B, but the legacy test's false
   discovery proportion rises, mostly not because of the zero-filling.** From 10 % (P0, P0r, P1, P2, P3, L1, L2) or 25 %
   (P4), 2 % in M1, M1b and L1b: +0.016 to +0.34 beyond the oracle's change (P1 VD 50 %: +0.34; P0 VD 25 %: +0.066).
   The comparator is serious wherever S is serious (26 of 26), and in the proteomics cells it raises the false discovery
   proportion as much or more (S - C from -0.049, P2 VD 10 %, to +0.001). In four lipid and metabolomics cells S exceeds the
   comparator, with 95 % intervals that exclude 0: L1 VD 25 % (+0.020), M1 VD 10 % (+0.017; 211 assessable replicates),
   M1b batch-block 10 % (+0.022) and 25 % (+0.041; 16 assessable replicates). There, 8 to 26 % of the rise is the
   zero-filling, on top of a serious rise that both share. Where few replicates have max(20, 3k) complete features (M2
   VD 50 % 2 of 400, M1b batch-block 25 % 16, M1 VD 10 % 211), the comparator's own classification says little about
   directions from complete features. Inferred: the shared part follows from correcting and testing features with
   missing values (5.9's falsification case).
3. **The group effect (C\*)** is removed materially in the lipid and metabolomics cells at value-dependent rates (L1,
   L1b, M1 from 10 %; L2, M2 from 25 %; M1b from 2 %, serious, 0.12), at random in L1 (2 %), M2 (10 %) and P4 (10 %,
   0.047), and in M1b batch-block (material at 10 %, serious at 25 %); the comparator removes about as much (within
   0.015) except in M1b (batch-block 10 %: 0.083 against 0.064; 25 %: 0.245 against 0.156; VD 2 %: 0.119 against 0.097,
   "partly").
4. **Batch-block missingness** (D5): P4b and L1b no material effect up to 10 %; P4b FDP serious at 25 %; M1b FDP serious from 10 %, C\* material at
   10 % and serious at 25 %, power material at 25 %.
5. **Refused by the route** (the package's guard before RUV-III-C, k above the controls observed in every sample,
   `R/func_lipid_s4_ruv_methods.R:678-682`): L1 VD 50 % (34 % of replicates), L1
   MCAR 10 % (100 %), M1 VD 25 % (7 %), VD 50 % (94 %), MCAR 2 % (12 %), MCAR 10 % (100 %), M1b VD 10 % (16 %).

## Classifications (deciding cells; family 456, Bonferroni t on R - 1 df at 1 - 0.05/456; R = 400)

**B, unwanted variance left** (. no material effect, M material, S serious, ? not determined, R refused by the route, - not in the cell's design)

| cell | vd0.1 | vd0.5 | vd2 | vd10 | vd25 | vd50 | mcar0.5 | mcar2 | mcar10 | bb2 | bb10 | bb25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P0 | . | . | . | . | . | . | . | M | S | - | - | - |
| P0r | . | . | . | . | . | . | . | M | S | - | - | - |
| P1 | . | . | . | . | . | . | . | M | S | - | - | - |
| P2 | . | . | . | . | . | . | M | M | M | - | - | - |
| P3 | . | . | . | . | . | . | . | . | M | - | - | - |
| P4 | . | . | . | . | . | . | . | . | M | - | - | - |
| P5 | . | . | . | ? | M | M | ? | M | M | - | - | - |
| L1 | . | . | . | . | . | R | . | . | R | - | - | - |
| L2 | . | . | . | . | . | . | . | M | S | - | - | - |
| M1 | . | . | . | . | R | R | . | R | R | - | - | - |
| M2 | . | . | . | . | . | . | ? | M | S | - | - | - |
| P4b | - | - | . | . | - | - | - | - | - | . | . | . |
| L1b | - | - | . | . | - | - | - | - | - | . | . | . |
| M1b | - | - | M | R | - | - | - | - | - | . | . | . |

**C*, group effect removed** (. no material effect, M material, S serious, ? not determined, R refused by the route, - not in the cell's design)

| cell | vd0.1 | vd0.5 | vd2 | vd10 | vd25 | vd50 | mcar0.5 | mcar2 | mcar10 | bb2 | bb10 | bb25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P0 | . | . | . | . | . | . | . | . | . | - | - | - |
| P0r | . | . | . | . | . | . | . | . | . | - | - | - |
| P1 | . | . | . | . | . | . | . | . | . | - | - | - |
| P2 | . | . | . | . | . | . | . | . | . | - | - | - |
| P3 | . | . | . | . | . | . | . | . | . | - | - | - |
| P4 | . | . | . | . | . | . | . | . | M | - | - | - |
| P5 | . | . | . | . | . | . | . | . | . | - | - | - |
| L1 | . | . | . | M | M | R | . | M | R | - | - | - |
| L2 | . | . | . | . | M | M | . | . | . | - | - | - |
| M1 | . | . | . | M | R | R | ? | R | R | - | - | - |
| M2 | . | . | . | ? | M | S | . | . | M | - | - | - |
| P4b | - | - | . | . | - | - | - | - | - | . | . | . |
| L1b | - | - | . | M | - | - | - | - | - | . | . | . |
| M1b | - | - | S | R | - | - | - | - | - | . | M | S |

**false discovery proportion** (. no material effect, M material, S serious, ? not determined, R refused by the route, - not in the cell's design)

| cell | vd0.1 | vd0.5 | vd2 | vd10 | vd25 | vd50 | mcar0.5 | mcar2 | mcar10 | bb2 | bb10 | bb25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P0 | . | . | . | M | S | S | . | . | . | - | - | - |
| P0r | . | . | . | M | S | S | . | . | . | - | - | - |
| P1 | . | . | . | S | S | S | . | . | . | - | - | - |
| P2 | . | . | . | S | S | S | . | . | M | - | - | - |
| P3 | . | . | . | M | S | S | . | . | . | - | - | - |
| P4 | . | . | . | . | M | S | . | . | S | - | - | - |
| P5 | . | . | . | . | . | . | . | . | . | - | - | - |
| L1 | . | . | ? | S | S | R | ? | S | R | - | - | - |
| L2 | . | . | . | M | S | S | . | . | . | - | - | - |
| M1 | . | ? | S | S | R | R | S | R | R | - | - | - |
| M2 | . | . | . | ? | S | S | . | . | . | - | - | - |
| P4b | - | - | . | ? | - | - | - | - | - | . | . | S |
| L1b | - | - | M | S | - | - | - | - | - | . | . | ? |
| M1b | - | - | S | R | - | - | - | - | - | ? | S | S |

**power lost** (. no material effect, M material, S serious, ? not determined, R refused by the route, - not in the cell's design)

| cell | vd0.1 | vd0.5 | vd2 | vd10 | vd25 | vd50 | mcar0.5 | mcar2 | mcar10 | bb2 | bb10 | bb25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P0 | . | . | . | M | M | . | . | M | M | - | - | - |
| P0r | . | . | . | M | M | . | . | M | M | - | - | - |
| P1 | . | . | . | . | . | . | . | M | S | - | - | - |
| P2 | . | . | ? | . | . | . | M | M | M | - | - | - |
| P3 | . | . | . | . | . | . | . | . | S | - | - | - |
| P4 | . | . | . | . | . | . | . | . | . | - | - | - |
| P5 | . | . | ? | M | S | S | ? | M | M | - | - | - |
| L1 | . | . | . | . | . | R | . | . | R | - | - | - |
| L2 | . | . | ? | . | . | . | M | S | S | - | - | - |
| M1 | . | . | . | . | R | R | . | R | R | - | - | - |
| M2 | . | . | M | ? | ? | M | M | S | S | - | - | - |
| P4b | - | - | . | . | - | - | - | - | - | . | . | . |
| L1b | - | - | . | . | - | - | - | - | - | . | . | . |
| M1b | - | - | . | R | - | - | - | - | - | . | . | M |

Full values (mean, sd, Monte Carlo error, limits, resolvable edges, C's classification, attribution with intervals and
the attributable share): `results/cells.tsv`, `results/attribution.tsv`.

## Safe rate per cell and mechanism (5.7; nominal rate, realised share after the filter beside it)

| cell | mechanism | safe rate (nominal) | realised share at it | highest measured |
|---|---|---|---|---|
| L1 | vd | 0.5 % | 0.20 % | 50 % |
| L1 | mcar | 0 % | 0.00 % | 10 % |
| L1b | vd | 0 % | 0.00 % | 10 % |
| L1b | bb | 10 % | 9.96 % | 25 % |
| L2 | vd | 0.5 % | 0.21 % | 50 % |
| L2 | mcar | 0 % | 0.00 % | 10 % |
| M1 | vd | 0.1 % | 0.04 % | 50 % |
| M1 | mcar | 0 % | 0.00 % | 10 % |
| M1b | vd | 0 % | 0.00 % | 10 % |
| M1b | bb | 0 % | 0.00 % | 25 % |
| M2 | vd | 0.5 % | 0.21 % | 50 % |
| M2 | mcar | 0 % | 0.00 % | 10 % |
| P0 | vd | 2 % | 1.55 % | 50 % |
| P0 | mcar | 0.5 % | 0.50 % | 10 % |
| P0r | vd | 2 % | 1.53 % | 50 % |
| P0r | mcar | 0.5 % | 0.50 % | 10 % |
| P1 | vd | 2 % | 1.41 % | 50 % |
| P1 | mcar | 0.5 % | 0.50 % | 10 % |
| P2 | vd | 0.5 % | 0.31 % | 50 % |
| P2 | mcar | 0 % | 0.00 % | 10 % |
| P3 | vd | 2 % | 1.62 % | 50 % |
| P3 | mcar | 2 % | 2.00 % | 10 % |
| P4 | vd | 10 % | 8.81 % | 50 % |
| P4 | mcar | 2 % | 2.00 % | 10 % |
| P4b | vd | 2 % | 1.74 % | 10 % |
| P4b | bb | 10 % | 8.48 % | 25 % |
| P5 | vd | 0.5 % | 0.38 % | 50 % |
| P5 | mcar | 0 % | 0.00 % | 10 % |

## Onset (lowest rate classified material and serious)

| cell | mechanism | estimand | lowest material | lowest serious | realised share at material onset |
|---|---|---|---|---|---|
| L1 | vd | Cstar | 10 % | - | 9.41 % |
| L1 | vd | fdp | 10 % | 10 % | 9.41 % |
| L1 | mcar | Cstar | 2 % | - | 2.01 % |
| L1 | mcar | fdp | 2 % | 2 % | 2.01 % |
| L1b | vd | Cstar | 10 % | - | 9.42 % |
| L1b | vd | fdp | 2 % | 10 % | 1.09 % |
| L2 | vd | Cstar | 25 % | - | 24.84 % |
| L2 | vd | fdp | 10 % | 25 % | 9.05 % |
| L2 | mcar | B | 2 % | 10 % | 2.03 % |
| L2 | mcar | power | 0.5 % | 2 % | 0.50 % |
| M1 | vd | Cstar | 10 % | - | 9.49 % |
| M1 | vd | fdp | 2 % | 2 % | 1.07 % |
| M1 | mcar | fdp | 0.5 % | 0.5 % | 0.50 % |
| M1b | vd | B | 2 % | - | 1.13 % |
| M1b | vd | Cstar | 2 % | 2 % | 1.13 % |
| M1b | vd | fdp | 2 % | 2 % | 1.13 % |
| M1b | bb | Cstar | 10 % | 25 % | 9.84 % |
| M1b | bb | fdp | 10 % | 10 % | 9.84 % |
| M1b | bb | power | 25 % | - | 25.06 % |
| M2 | vd | Cstar | 25 % | 50 % | 24.92 % |
| M2 | vd | fdp | 25 % | 25 % | 24.92 % |
| M2 | vd | power | 2 % | - | 1.11 % |
| M2 | mcar | B | 2 % | 10 % | 2.03 % |
| M2 | mcar | Cstar | 10 % | - | 10.08 % |
| M2 | mcar | power | 0.5 % | 2 % | 0.51 % |
| P0 | vd | fdp | 10 % | 25 % | 7.46 % |
| P0 | vd | power | 10 % | - | 7.46 % |
| P0 | mcar | B | 2 % | 10 % | 2.00 % |
| P0 | mcar | power | 2 % | - | 2.00 % |
| P0r | vd | fdp | 10 % | 25 % | 7.44 % |
| P0r | vd | power | 10 % | - | 7.44 % |
| P0r | mcar | B | 2 % | 10 % | 2.00 % |
| P0r | mcar | power | 2 % | - | 2.00 % |
| P1 | vd | fdp | 10 % | 10 % | 6.25 % |
| P1 | mcar | B | 2 % | 10 % | 2.00 % |
| P1 | mcar | power | 2 % | 10 % | 2.00 % |
| P2 | vd | fdp | 10 % | 10 % | 4.60 % |
| P2 | mcar | B | 0.5 % | - | 0.49 % |
| P2 | mcar | fdp | 10 % | - | 8.30 % |
| P2 | mcar | power | 0.5 % | - | 0.49 % |
| P3 | vd | fdp | 10 % | 25 % | 8.07 % |
| P3 | mcar | B | 10 % | - | 10.00 % |
| P3 | mcar | power | 10 % | 10 % | 10.00 % |
| P4 | vd | fdp | 25 % | 50 % | 21.67 % |
| P4 | mcar | B | 10 % | - | 10.00 % |
| P4 | mcar | Cstar | 10 % | - | 10.00 % |
| P4 | mcar | fdp | 10 % | 10 % | 10.00 % |
| P4b | bb | fdp | 25 % | 25 % | 20.24 % |
| P5 | vd | B | 25 % | - | 18.34 % |
| P5 | vd | power | 10 % | 25 % | 7.89 % |
| P5 | mcar | B | 2 % | - | 2.00 % |
| P5 | mcar | power | 2 % | - | 2.00 % |

## Reported, not classified

- **P0lo / P0hi** (`results/reported.tsv`): level dependence as predicted at random missingness. B at MCAR 2 % is 0.027
  at level 14, 0.16 at 20 (P0) and 0.29 at 26; at value-dependent rates B stays below 0.03 at all three levels.
- **Run-aligned masks** (D5, reported): B 0.14 to 0.24 in P0 and P1 at 2 % and 10 % (inferred: the zero-filling direction
  follows the run with most zeros); 0.02 to 0.03 in L2, M2; power loss 0.04 to 0.13.
- **Part B** (`results/partb.tsv`, R = 100). The declared designs' DA contrasts were measured for lipidomics and
  metabolomics only: on L2 and M2, FDP moved by at most +0.006 and power by at most 0.026 up to VD 10 % (n = 100 each);
  at VD 50 % L2's power loss is 0.17 on the 62 of 100 replicates that returned a table, but in 52 of them one assay's
  DA failed and its planted features count as not called; on the 10 replicates where both assays completed it is
  0.11 ± 0.06 (mean ± 1.96 Monte Carlo errors, n = 10); M2 has no value (see the note below). The proteomics part-B DA contrasts were first lost to a script defect (instrument events) and were re-run under
  VAL-98 (1) (six jobs at df60efa9, 100 replicates, the reference arm complete in every replicate); the values, oracle-
  referenced as part A's, mean ± 1.96 Monte Carlo errors:

  | route | cell | mask | B | C* | FDP | power lost |
  |---|---|---|---|---|---|---|
  | declared | P0 | VD 2 % | +0.012 ± 0.002 | +0.001 ± 0.000 | +0.000 ± 0.001 | +0.001 ± 0.001 |
  | declared | P0 | VD 10 % | +0.006 ± 0.002 | +0.004 ± 0.000 | +0.001 ± 0.002 | +0.011 ± 0.002 |
  | declared | P0 | VD 50 % | -0.018 ± 0.003 | +0.017 ± 0.002 | +0.009 ± 0.004 | +0.090 ± 0.005 |
  | legacy | P0 | VD 2 % | -0.000 ± 0.042 | +0.002 ± 0.000 | -0.001 ± 0.007 | +0.014 ± 0.013 |
  | legacy | P0 | VD 10 % | +0.030 ± 0.063 | +0.008 ± 0.001 | +0.000 ± 0.012 | +0.050 ± 0.020 |
  | legacy | P0 | VD 50 % | +0.123 ± 0.216 | +0.021 ± 0.003 | +0.065 ± 0.026 | +0.105 ± 0.036 |
  | legacy | P1 | VD 2 % | -0.038 ± 0.050 | +0.002 ± 0.000 | +0.013 ± 0.015 | -0.017 ± 0.033 |
  | legacy | P1 | VD 10 % | +0.015 ± 0.065 | +0.009 ± 0.001 | +0.009 ± 0.023 | +0.044 ± 0.045 |
  | legacy | P4b | VD 2 % | -0.193 ± 0.067 | +0.002 ± 0.000 | +0.013 ± 0.005 | -0.001 ± 0.003 |
  | legacy | P4b | VD 10 % | -0.266 ± 0.070 | +0.006 ± 0.001 | +0.022 ± 0.006 | +0.013 ± 0.005 |

  The declared proteomics design (W from complete features in the model) moves FDP by at most +0.001 and power by at
  most 0.011 up to VD 10 % (at VD 50 %: FDP +0.009 ± 0.004, power 0.090 ± 0.005). The legacy route with its own sweep and k is noisier and, at VD 50 % in
  P0, raises FDP by +0.065 and loses 0.105 of power; in P4b its FDP rises by +0.013 to +0.022 at VD 2 and 10 %. Part B's B, C*
  and A1 are complete where the route corrected. The legacy route with its own sweep is noisy: its k and controls change
  with the mask, and B is unstable where the sweep picks few controls. The legacy lipidomics and metabolomics method
  returns an assay uncorrected (with a log warning; the user sees no refusal) at L2 VD 2, 10 and 50 % in 6, 78 and 100 %
  of replicates, M2 VD 10 and 50 % in 42 and 100 %, L1b VD 10 % in 11 %, and M1b VD 2 and 10 % in 2 and 34 %.
  `sim-k185.R:856-864` records this as `refused: assay_not_corrected`; section 7.3 would count it as `not_applied`.
- **Part B, declared lipid and metabolomics DA at VD 50 % (observed; owner VAL-116 (6)).** Counted from the part-B rows
  (`rows/b-*.tsv.gz`; every part-B cell and route checked for `da_failed`): only the declared route at VD 50 % is
  affected. L2: in 52 of 100 replicates one assay's DA stopped (`da_failed` = 1), and in 38 more no DA table was returned
  at all; M2: no DA table in 100 of 100. In the rows the message of a whole-table failure was not recorded (the
  `tryCatch` around the route's DA kept none). No message was observed here. What was observed: every assay whose DA
  failed is an assay for which the declared route estimated k = 0 (no factor; L2: 52 replicates with one such assay,
  all 52 failed on it; 38 with both, no table; M2: k = 0 in 100 of 100, no table), and every assay with k >= 1
  completed. With no factor the DA takes the plain limma path (`R/func_lipid_da_model.R:513-528`), where KIT-196
  diagnoses the stop `prior.weights contain NA values` on a feature with no value in any sample; that it is this stop
  is inferred, not shown (KIT-196 2.6). Reported beside it: at VD 50 % the calibrated design (`num_sv_be`) estimated
  k = 0 for at least one assay in 90 of 100 L2 and all 100 M2 replicates, although two factors are planted; at VD 2
  and 10 % it estimated 2. No other part-B cell (legacy or declared, L2, M2, L1b, M1b, any mask) has a failed assay.
  The L2 and M2 VD 50 % declared DA figures therefore do not stand as measurements (L2: 10 complete replicates, 52 with
  one assay, 38 with none; M2: none); the replicates that returned are not a random subset. These cells are re-run after KIT-196 merges
  only if their figures are cited. Part A is not affected (G2; no deciding DA failed).
- **Reported-only summary** (`results/reported-extra.tsv`, `results/qvalue-methods.tsv`; `summ-k185.R`, VAL-98 (2)):
  in the features whose observed runs leave exactly k residual dimensions (the closed-form branch) the shipped
  correction's FDP rises by +0.29 to +0.52 beyond the oracle at VD 25 and 50 % (P0, P1, L2, M2), and by +0.06 to +0.30
  in the features with more than k; the comparator's rise in the same strata is as large (within 0.03; M2 VD 25 %
  +0.38 against +0.29), as expected, since the closed-form branch takes no zero-filled value (ticket 2.1): this part is
  not the zero-filling; the stratum figures are per-replicate proportions over few calls (at VD 25 %, P0 4.4, P1 62,
  L2 3.3, M2 1.2 calls per replicate), not shares of all false calls (those are 5 to 25 %, below); the oracle's own FDP
  change is at most +0.08 over all filtered features (M2 VD 50 %) and at most +0.10 within the closed-form stratum
  (L2 VD 50 %), both already subtracted; the shipped correction loses 2 to 3 % of the filtered proteins at VD 25-50 %
  in P0 and P1, 8 to 41 % of the features in L2 and M2 (the comparator loses the same features). The calls per
  replicate and the in-stratum oracle change are counted from the committed rows (`rows/a-*`), not from a results table.
- **Sensitivity** (`results/sensitivity.tsv`, P0 and M2, R = 200): with strength 1 (weaker factors) B at VD 2 % rises
  to 0.33 to 0.39 and MCAR 2 % to 0.48 to 0.62 (larger than any deciding value; the material edge is crossed at VD 2 %
  where part A is "no material effect"); strength 10 leaves B near 0. Against part A's MCAR 2 % B (P0 0.16, M2 0.19),
  k = 1 gives 0.06 and 0.09, k = 4 gives 0.04 and 0.03 (each a change larger than the 0.05 edge), cyclic loess 0.17 and
  0.19 (no change). FDP and power move by more than their edges in some variants (P0 k = 4 VD 10 % FDP +0.059 against
  +0.016; M2 strength 1 power at 2 % 0.25 to 0.35).
  Beside the D2 reading: **at weaker unwanted variation the zero-filling matters at value-dependent rates too.**

## D2: reading of KIT-160 to KIT-166 (O3)

As O2's mapping, plus: where a data set's realised share is read at a rate classified "serious" in its matched cell and
mechanism, the RUV-III-C-corrected results of that route are reported as measured but not as evidence that the
correction works, and a remedy is decided before any RUV-based claim on that route. Each data set is placed by its
ticket's `ruv-input.tsv` (11.1). Safe rates are low: 0 to 2 % at random in every cell, 0.1 to 2 % value-dependent in
most (0 in L1b and M1b; 10 % only in P4), so (inferred) real data sets will mostly be read above the safe rate. The declared designs' DA is not
exposed to the zero-filling by construction (E6). Part B shows its change is small up to VD 10 % (lipidomics and
metabolomics; proteomics after the VAL-98 (1) re-run); their corrected matrices (figures, export, the
integration input) carry the limitation measured here. KIT-161's PXD000279 is read against P2.

How D2 applies to these results, as the owner decided (VAL-98 (3)):
1. **Mechanism.** The mechanism with the lower safe rate in the matched cell. On a tie (P3: 2 % and 2 %; M1b: 0 and 0)
   VAL-98 (3)'s "the lower applies" does not choose, and the two mechanisms do not read the same above the safe rate
   (`results/cells.tsv`): P3 at 10 % is "serious" at random (power; B material) but "material" value-dependent (FDP;
   realised 8.1 %); M1b at 2 % is "serious" value-dependent (C\*, FDP; B material; realised 1.13 %) but "not
   determined" batch-block (FDP; realised 1.3 %). P2 is not a tie: the rule's mechanism (at random, safe rate 0) is at
   most "material" up to 10 %, while value-dependent missingness is "serious" from 10 % (realised 4.6 %). **D2's reading
   for P3, M1b and P2 is awaiting the owner's ruling** (close audit C02); KIT-161 (P3, P2) and KIT-166 (M1b) are not
   read under this item until it is given.
2. **Rates.** A data set's realised share (11.1 item 2, after the route's filter) is compared with the cells' realised
   shares after the filter (the "realised share" columns), not with the nominal rates.
3. **Coverage.** Where the matched cell and mechanism were measured only up to 10 % (at random in L1, L2, M1, M2, P0,
   P0r, P1, P2, P4 and P5; value-dependent in P4b and L1b), a data set whose realised share is above 10 % is "not covered",
   and its RUV-III-C result stays held (VAL-56).
4. For proteomics, part B's DA contrasts are those of the re-run (VAL-98 (1), the Part B bullet above).
In practice, a data set on an affected route is read as evidence of correction only at or below the safe rate of its
cell's D2 mechanism (0 to 2 %). **No remedy is decided here**: VAL-56 keeps the remedy for its own decision on this
measurement.

## Acceptance criteria

| # | Evidence | Met |
|---|---|---|
| A1 | `git diff --stat ffdfdfb7..HEAD`: only `docs/qa/workflow-execution/kit/kit-185/` and the ledger row | met |
| A2 | `lifecycle.md`; L1 E1 to E4, E6, E7 pass after one instrument run (E3 mechanics); E1-neg failed on P2 at the old filter (VAL-71) and passed after KIT-187 (P2 lifecycle 45314701113: 1.85 to 2.54); G1, G2 pass; G3 no re-split | met |
| A3 | declaration 432f703b committed alone after L1 and VAL-67, VAL-71, before S; script hashes match the files S ran (section 1); section 7 dated, scoring only; section 9 (df60efa9, with the two-line fix) before the part-B re-run; section 10 (ff73dbff) before the summary ran; section 11 committed alone (6779d56e) before P2's lifecycle and jobs | met |
| A4 | every job's replicates written (P5 jobs 1 and 3 timed out; their missing replicates 117-133 and 380-400 re-sent once, VAL-54, run 45268841462); refused counts in cells.tsv; E5 identical; E2 470 of 470 recorded pass (P2's 20 included; 10 not run, P5 timed-out jobs); E8 pass | met |
| A5 | `results/` | met. The deciding estimands (B, C*, FDP, power), A1 and the null rate for all 14 cells, with mean, sd, Monte Carlo error, R and limits (`cells.tsv`); the estimands the pinned scoring does not summarise (features lost by residual dimension and planted status, B2, the incomplete-feature recovery, every canonical correlation, the comp/inc and residual-dimension strata, the oracle's own change, the uncorrected arm's departure, the arms' mechanics, the q-value methods and their S/S0 mismatches) in the reported-only summary `reported-extra.tsv` and `qvalue-methods.tsv` (`summ-k185.R`, VAL-98 (2)). |
| A6 | this file | met |
| A7 | ledger row `VAL-140` on the branch | met |
| A8 | `runs.tsv`: 443.3 core-hours in all, lifecycle, P2, the part-B re-run and scoring included (expected about 514, ceiling 966) | met |
| A9 | `git diff ffdfdfb7 147adde8` (current `v0.9-staging`; the release commit is not yet made) over the files of 2.2 and the filters: three files changed, all KIT-172's log step: `R/utils_workflow_omics_declarations.R` (251-320), `R/func_lipid_s4_normalization_methods.R` (the log transform, 470-610) and `R/func_lipid_qc.R` (`computeLipidFilteringMetrics` passes `logged`); it does not touch the RUV-III-C call or its inputs (the data enter after the log step) | met at 147adde8; to recheck at the release commit |

## Observed and inferred

Observed: every number above, from the committed rows through `score-k185.R`, except the sentences marked inferred. Inferred: the cause of the value-dependent
FDP rise (the correction of incomplete features, which the comparator shares; the residual dimension strata in the rows
can test this); the generality beyond this generator (13). Not from the pinned scoring (the implementation audit's
count from the rows): the closed-form features (exactly k residual dimensions) carry about 5 % (P0 VD 25 %), 21 % (L2 VD
25 %) and 25 % (P1 VD 50 %) of S's excess false calls; the rest sit in features with more than k residual dimensions, so
5.8's closed-form expectation is not the main source.

## Limits

- One synthetic generator: planted factors of total strength 3 (1 and 10 in the sensitivity panel), orthogonal to the
  groups except P0r; normal residuals; logistic left-censoring or uniform missingness; batch-aligned missingness only
  in the D5 cells and run-aligned only in the reported masks; no peptide-to-protein roll-up.
- Levels about 20 in the deciding cells (14 and 26 reported); a data set whose median log2 value is outside 14 to 26 is
  not covered.
- k fixed at the truth in the deciding part; real use estimates it (part B and the sensitivity panel, reported).
- Normalisation off except one sensitivity variant; the DA is the legacy route's test.
- The measurement shows what the correction does to these data, not what it did to any real data set, nor which remedy
  is best in the package.

## Instrument events

- L1: E3 parameters; part-B comparator (assay split, control order) — one instrument run (L1fix).
- S: P5 ran at 298.5 s per replicate; L1 measured 306 s against the 204.8 s estimate (+49 %), above G3's 40 %, and G3
  was applied to an inferred 260 s instead (a deviation from G3 as written, disclosed in `lifecycle.md` before S). Jobs 1
  and 3 hit their 175-minute timeouts with 116 and 112 of 133 replicates written. The 38 unwritten replicates (117-133,
  380-400) were re-sent once under VAL-54 (run 45268841462, same commit and seeds), in place of 5.7's "did not complete"
  rule, which would have made P5 "not determined". **A recorded departure, accepted by the owner: VAL-98 (4).** The
  timed-out jobs' rows were assembled outside R from their uploaded part files (text identical, verified by the
  implementation audit). Those jobs wrote no checks, so E2 did not run for P5 replicates 20-100 and 280-360 (every 20th);
  "E2 470 of 470" counts the checks that exist.
- Part B (found by the implementation audit; re-run, see the end of this item): the reference DA arm on the `none` state restricted to F_filt
  (`sim-k185.R:994`) failed in most proteomics replicates and the `tryCatch` kept no error. The P0, P1 and P4b part-B FDP
  and power contrasts are missing, or rest on 7 to 20 of 100 replicates. Reported part only; no classification changes.
  Inferred, by reading only (no run): `subsetState` (`sim-k185.R:897-909`) subsets `protein_quant_table` by a logical vector and keeps
  the surviving row names (for example 1, 2, 5, ...), which are no longer the automatic ones; the DA then calls
  `tibble::column_to_rownames(theObject@protein_quant_table, ...)` (`R/func_prot_da_model_methods.R:287`), which stops
  when the data frame already has row names (the likely message: "`.data` must be a data frame without row names"). The
  replicates that returned a table are those where the subset left automatic row names. A fixed job needs
  `rownames(state@protein_quant_table) <- NULL` after the subset in `subsetState`, a new hash in a dated declaration
  section, and the six proteomics part-B jobs re-run (b-p0, b-p1, b-p4b, two each). Done under VAL-98 (1): declaration section 9,
  lifecycle 45291628697 (the reference arm returned a table in every cell), six jobs at df60efa9 (45.7 core-hours), the
  reference DA complete in all 100 replicates of each cell.
- Scoring: E5 and E8 compared tables by inferred column type; fixed (declaration section 7), re-scored once; every
  scored table identical between the two scoring runs (owner told 2026-11-24; ticket section 10).
- Job splits used `round()` boundaries (1-133, 134-267, 268-400) where the declaration's section 6 wrote 1-134 etc.;
  the same replicates.

## Not checked

- The release commit (A9 is stated at 147adde8).
- The run path's own declared lipid/metab DA kernel (E6 on L2 mirrors its input rule).
- Real data; other generators, residual distributions, normalisers beyond cyclic loess.
- The message of the declared lipid and metabolite part-B DA failures (not recorded); that they are KIT-196's stop is
  inferred from the k = 0 coincidence and the code path, not shown.
- D2's reading for P3, M1b and P2 (close audit C02): awaiting the owner's ruling.

## For the change report

No behaviour changed in this ticket; it measured a limitation of the legacy pipeline that a later remedy may change.

What it measured. On the classic proteomics routes (DIA without limpa, LFQ, TMT) and on the lipidomics and metabolomics
routes, MultiScholaR passes matrices that still contain missing values to RUV-III-C (`R/func_prot_s4_norm_methods.R:248-254`,
`R/func_lipid_s4_ruv_methods.R:691-697` at ffdfdfb7). RUV-III-C (RUVIIIC 1.0.19, `src/RUVIIIC_Varying.cpp:80-85`) sets
those missing values to 0 when it estimates the unwanted directions. DIA with limpa and transcriptomics do not use
RUV-III-C.

Who is affected. Inferred from the simulation (real data were not measured; see Limits): users of these routes whose data have more missing values than the matched shape's safe rate (0 to
2 % in most shapes; 10 % only in the 50-run TMT-like shape) get RUV-III-C-corrected values, and on the legacy test
calls, affected as tabulated above. The declared designs' test (W from complete features in the model) does not see
the zero-filling by construction; their corrected matrices (figures, export, the integration input) do.

What the measurement shows (outcome O3, VAL-140). At random missingness from 0.5 % (two groups of three) or 2 % (most shapes), the zero-filling itself
leaves 0.09 to 0.40 of the planted unwanted variation in the data and costs the legacy test up to 0.45 of its power beyond what the missing values cost the truth;
directions from complete features remove this. At value-dependent missingness from 10 % (from 2 % in three shapes), the legacy test's false
discovery proportion rises by +0.016 to +0.34 beyond the truth's own change, and the complete-feature version rises
as much, so that rise is mostly not the zero-filling. The complete-feature version is not a ready remedy: at 10 %
random missingness its own false discovery proportion is seriously raised in six shapes. With the declared design at
value-dependent 50 %, small lipid and metabolite designs estimated no factor (k = 0) in 90 and 100 of 100 data sets,
and their DA did not complete for those assays (cause inferred, see Part B).

Evidence. Rows `rows/`, scored tables `results/` (final scoring run 45333019858), runs and commits `runs.tsv`; the
declaration `declaration.md` (432f703b, with dated sections 7 to 11).

Limits. One synthetic generator (log2 levels about 20; two planted unwanted factors of total strength 3; k at the
truth; no normalisation; the legacy route's limma test); 14 study shapes, 400 data sets each; no real data; other
designs, k, normalisers and tests not measured. No remedy is decided (VAL-56).


## Merge record (coordinator)

Merged into `v0.9-staging` on 2026-11-24 from af6bad3f (base 3d4a585d) after the closure audit's re-check 1 (register `AUDIT-KIT-185-close-20261124`): C01 and C03 to C07 closed; R1-01's two phrase edits applied in the merge commit; the ledger row is VAL-140. Nothing ships with this merge: no code, no NEWS entry, no version change. Open with the owner: C02 (how D2 reads P3, M1b and P2; KIT-161 and KIT-166 are not read under that item until it is given) and the remedy (VAL-56).

## Note added 2026-11-24 (coordinator): the D2 ruling

Where this file says D2's reading for P3, M1b and P2 is "awaiting the owner's ruling", the ruling has since been given (methods ledger VAL-143): each study shape is read at the more severe of the mechanisms measured for it. The earlier sentences are left as written; this note supersedes them. KIT-161 (P3, P2) and KIT-166 (M1b) are read under VAL-143. The remedy is not decided here: a separate ticket measures the candidates first (VAL-144).
