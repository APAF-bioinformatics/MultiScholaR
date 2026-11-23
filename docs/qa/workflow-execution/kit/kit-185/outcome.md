# KIT-185 outcome: RUV-III-C and missing values, measured (P2 held for KIT-187)

Ticket revision 1.2; declaration `declaration.md` (commit 432f703b, sections 1 to 6; section 7 dated, after the results).
Package code: `v0.9-staging` at ffdfdfb76a9eacdc681f2098f75ac574e2cf8af4. Runs and core-hours: `runs.tsv`. Rows:
`rows/`; the jobs' checks `checks/`; timing and provenance `jobs/`; scored tables `results/` (re-score run 45283782973).

## Answer

**Outcome O3** (the ticket's rule, 5.7) on the 13 deciding cells run: 420 of the 456 declared classifications, 42
"serious", 46 "material", 287 "no material effect", 17 "not determined", 28 "refused by the route". **P2 is held**
for KIT-187 (VAL-71): its 36 classifications, its onset and safe rate, and therefore the outcome's final form wait for
it. P2 can only add classifications; O3 is already the most severe outcome, so the overall outcome stays O3 whatever P2
shows, but P2's own reading (KIT-161, PXD000279) waits.

Observed, in plain words:
1. **At random (MCAR) missingness, the zero-filling itself damages the correction.** From 2 % (P0, P0r, P1, L2, M2,
   P5) or 10 % (P3, P4) the shipped correction leaves 0.12 to 0.40 of the planted unwanted variance (B) that it removed
   without missing values, and the legacy test loses up to 0.45 of power beyond the oracle's own loss (L2 MCAR 10 %).
   The comparator with directions from complete features leaves none of it (C's Delta between -0.06 and +0.05; classified
   "no material effect" in every such case): 35 of the 88 material or serious deciding classifications are "attributable
   to the zero-filling" (29 at random, 6 value-dependent: P5 VD 25 and 50 % B, P0 and P0r VD 25 % power, M1b VD 2 % B, M2 VD 2 % power).
2. **At value-dependent (left-censored) missingness, the zero-filling does little to B, but the legacy test's false
   discovery proportion rises.** From 10 % (P0, P0r, P1, P3, L1, L2) or 25 % (P4), 2 % in M1, M1b and L1b: +0.016 to
   +0.34 beyond the oracle's change (P1 VD 50 %: +0.34; P0 VD 25 %: +0.066). The comparator raises it as much or more
   (C serious wherever S is serious), so this is **not** the zero-filling: it follows from correcting and testing
   features with missing values (attribution "not determined" or "not attributable"; section 5.9's falsification case).
3. **The group effect (C\*)** is removed materially only in the lipid and metabolomics cells from VD 10 % (L1, L1b,
   M1, M2) and in M1b from 2 % (serious, 0.12); the comparator removes about as much (not the zero-filling, except M1b
   VD 2 %, "partly").
4. **Batch-block missingness** (D5): P4b and L1b no material effect up to 10 %; P4b FDP serious at 25 %; M1b C\* and
   FDP serious from 10 %.
5. **Refused by the route** (RUV-III-C's own rule, k above the complete controls): L1 VD 50 % (34 % of replicates), L1
   MCAR 10 % (100 %), M1 VD 25 % (7 %), VD 50 % (94 %), MCAR 2 % (12 %), MCAR 10 % (100 %), M1b VD 10 % (16 %).

## Classifications (deciding cells; family 456, Bonferroni t on R - 1 df at 1 - 0.05/456; R = 400)

**B, unwanted variance left** (. no material effect, M material, S serious, ? not determined, R refused by the route, - not in the cell's design)

| cell | vd0.1 | vd0.5 | vd2 | vd10 | vd25 | vd50 | mcar0.5 | mcar2 | mcar10 | bb2 | bb10 | bb25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P0 | . | . | . | . | . | . | . | M | S | - | - | - |
| P0r | . | . | . | . | . | . | . | M | S | - | - | - |
| P1 | . | . | . | . | . | . | . | M | S | - | - | - |
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
- **Run-aligned masks** (D5, reported): B 0.14 to 0.24 in P0 and P1 at 2 % and 10 % (the zero-filling direction follows
  the run with most zeros); 0.02 to 0.03 in L2, M2; power loss 0.04 to 0.13.
- **Part B** (`results/partb.tsv`, R = 100): the declared designs (P0, L2, M2), which test the uncorrected values with
  W from complete features, move FDP by at most +0.006 and power by at most 0.03 up to VD 10 % (L2 VD 50 % power 0.17).
  The legacy route with its own sweep is noisy (its k and controls change with the mask; B unstable where the sweep
  picks few controls) and refuses assays at L2 and M2 VD 10 and 50 % ("assay_not_corrected" in 22 to 92 % of arms).
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
most (10 % only in P4), so real data sets will mostly be read above the safe rate. The declared designs' DA is not
exposed by construction (E6) and part B shows its change is small; their corrected matrices (figures, export, the
integration input) carry the limitation measured here. KIT-161's PXD000279 (P2's shape) waits for KIT-187 and P2.

## Acceptance criteria

| # | Evidence | Met |
|---|---|---|
| A1 | `git diff --stat ffdfdfb7..HEAD`: only `docs/qa/workflow-execution/kit/kit-185/` and the ledger row | met |
| A2 | `lifecycle.md`; L1 E1 to E4, E6, E7 pass after one instrument run (E3 mechanics); E1-neg fails on P2 (VAL-71, P2 held); G1, G2 pass; G3 no re-split | met except P2 (held) |
| A3 | declaration 432f703b committed alone after L1 and VAL-67, VAL-71, before S; script hashes match the files S ran (section 1); section 7 dated, scoring only | met |
| A4 | every job's replicates written (P5 jobs 1 and 3 timed out; their missing replicates 117-133 and 380-400 re-sent once, VAL-54, run 45268841462); refused counts in cells.tsv; E5 identical; E2 450 of 450 pass; E8 pass | met |
| A5 | `results/` | met for 13 cells; P2 waits |
| A6 | this file | met |
| A7 | ledger row `VAL-<merge>` on the branch | met |
| A8 | `runs.tsv`: 375.8 core-hours in all, lifecycle and scoring included (expected about 514 with P2, ceiling 966) | met |
| A9 | `git diff ffdfdfb7 147adde8` (current `v0.9-staging`; the release commit is not yet made) over the files of 2.2 and the filters: one change, `R/utils_workflow_omics_declarations.R` (KIT-172's log step, lines 251-320); it does not touch the RUV-III-C call or its inputs (the data enter after the log step) | met at 147adde8; to recheck at the release commit |

## Observed and inferred

Observed: every number above, from the committed rows through `score-k185.R`. Inferred: the cause of the value-dependent
FDP rise (the correction of incomplete features, which the comparator shares; the residual dimension strata in the rows
can test this); the generality beyond this generator (13).

## Instrument events

- L1: E3 parameters; part-B comparator (assay split, control order) — one instrument run (L1fix).
- S: P5 planned 302 s per replicate against 260 inferred at G3; jobs 1 and 3 hit their 175-minute timeouts with 116 and
  112 of 133 replicates written; the 38 replicates left were re-sent once (VAL-54). G3's inference understated P5.
- Scoring: E5 and E8 compared tables by inferred column type; fixed (declaration section 7), re-scored once; every
  scored table identical between the two scoring runs.
- Job splits used `round()` boundaries (1-133, 134-267, 268-400) where the declaration's section 6 wrote 1-134 etc.;
  the same replicates.

## Not checked

- P2 (held for KIT-187).
- The release commit (A9 is stated at 147adde8).
- The run path's own declared lipid/metab DA kernel (E6 on L2 mirrors its input rule).
- Real data; other generators, residual distributions, normalisers beyond cyclic loess.
