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
   discovery proportion rises, mostly not because of the zero-filling.** From 10 % (P0, P0r, P1, P3, L1, L2) or 25 %
   (P4), 2 % in M1, M1b and L1b: +0.016 to +0.34 beyond the oracle's change (P1 VD 50 %: +0.34; P0 VD 25 %: +0.066).
   The comparator is serious wherever S is serious (26 of 26), and in the proteomics cells it raises the false discovery
   proportion as much or more (S - C from -0.023 to +0.001). In four lipid and metabolomics cells S exceeds the
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
- **Run-aligned masks** (D5, reported): B 0.14 to 0.24 in P0 and P1 at 2 % and 10 % (inferred: the zero-filling direction
  follows the run with most zeros); 0.02 to 0.03 in L2, M2; power loss 0.04 to 0.13.
- **Part B** (`results/partb.tsv`, R = 100). The declared designs' DA contrasts were measured for lipidomics and
  metabolomics only: on L2 and M2, FDP moved by at most +0.006 and power by at most 0.026 up to VD 10 % (n = 100 each);
  at VD 50 % L2's power loss is 0.17 on 62 of 100 replicates, and M2 has no value. **For proteomics the part-B DA
  contrasts were not measured** (instrument defect, see Instrument events). The reference arm (the route's DA on the
  `none` state restricted to F_filt, `sim-k185.R:993-995`) returned no table in most replicates and the error was not
  recorded, so `partb.tsv` holds P0 declared n = 20, 0, 0 (VD 2, 10, 50 %), P0 legacy 11, 0, 0, P1 legacy 0, 0 and P4b
  legacy 7, 0. These few replicates are not a random subset and are not read. The declared proteomics design's DA is not
  exposed to the zero-filling by construction (E6); its change with missing values is not measured here. Part B's B, C*
  and A1 are complete where the route corrected. The legacy route with its own sweep is noisy: its k and controls change
  with the mask, and B is unstable where the sweep picks few controls. The legacy lipidomics and metabolomics method
  returns an assay uncorrected (with a log warning; the user sees no refusal) at L2 VD 2, 10 and 50 % in 6, 78 and 100 %
  of replicates, M2 VD 10 and 50 % in 42 and 100 %, L1b VD 10 % in 11 %, and M1b VD 2 and 10 % in 2 and 34 %.
  `sim-k185.R:856-864` records this as `refused: assay_not_corrected`; section 7.3 would count it as `not_applied`.
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
exposed to the zero-filling by construction (E6). Part B shows its change is small for lipidomics and metabolomics up
to VD 10 %, and did not measure it for proteomics (instrument defect), so for proteomics D2's "larger of part A's and
part B's Delta_tot" uses part A alone; their corrected matrices (figures, export, the
integration input) carry the limitation measured here. KIT-161's PXD000279 (P2's shape) waits for KIT-187 and P2.

How D2 applies to these results. (1) In 11 of the 13 cells the mechanism with the lower safe rate (at random in nine
cells, value-dependent in P4b and L1b) was measured up to 10 %, so a data set whose realised share exceeds 10 % is "not
covered" under D2 and its RUV-III-C result stays held (VAL-56). (2) In P3 and M1b both mechanisms have the same safe
rate; which one applies is the owner's decision. (3) A data set's realised share (11.1 item 2, after the route's filter)
is compared with the cells' realised shares (the "realised share" columns), not with the nominal rates, subject to the
owner's confirmation. (4) For proteomics, part B gave no DA contrast (instrument defect), so the reading uses part A
alone. In practice, a data set on an affected route is read as evidence of correction only at or below the safe rate
of its cell's D2 mechanism (0 to 2 %). Open for the owner: the mechanism when safe rates tie (P3, M1b); realised against
nominal rates; data sets above 10 % on the at-random mechanism ("not covered" as D2 is written).

## Acceptance criteria

| # | Evidence | Met |
|---|---|---|
| A1 | `git diff --stat ffdfdfb7..HEAD`: only `docs/qa/workflow-execution/kit/kit-185/` and the ledger row | met |
| A2 | `lifecycle.md`; L1 E1 to E4, E6, E7 pass after one instrument run (E3 mechanics); E1-neg fails on P2 (VAL-71, P2 held); G1, G2 pass; G3 no re-split | met except P2 (held) |
| A3 | declaration 432f703b committed alone after L1 and VAL-67, VAL-71, before S; script hashes match the files S ran (section 1); section 7 dated, scoring only | met |
| A4 | every job's replicates written (P5 jobs 1 and 3 timed out; their missing replicates 117-133 and 380-400 re-sent once, VAL-54, run 45268841462); refused counts in cells.tsv; E5 identical; E2 450 of 450 recorded pass (10 not run, P5 timed-out jobs); E8 pass | met |
| A5 | `results/` | partly met. The deciding estimands (B, C*, FDP, power), A1 and the null rate for 13 cells, with mean, sd, Monte Carlo error, R and limits. Not summarised by the pinned scoring (`score-k185.R:113` covers six estimands): (f) features lost, B2, the incomplete-feature recovery, the full canonical correlations, the comp/inc and residual-dimension strata, the oracle's own change and the uncorrected arm's departure (5.5), the per-arm per-run shares and residual-dimension counts (5.2), and the q-value-method mismatches (5.6: 733 BH and 19 mixed DA calls in part A, in M1, M1b, M2 and L2). Their per-replicate values are in `rows/`. P2 waits. |
| A6 | this file | met |
| A7 | ledger row `VAL-<merge>` on the branch | met |
| A8 | `runs.tsv`: 375.8 core-hours in all, lifecycle and scoring included (expected about 514 with P2, ceiling 966) | met |
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
  rule, which would have made P5 "not determined". **This is a departure awaiting the owner's ratification.** The
  timed-out jobs' rows were assembled outside R from their uploaded part files (text identical, verified by the
  implementation audit). Those jobs wrote no checks, so E2 did not run for P5 replicates 20-100 and 280-360 (every 20th);
  "E2 450 of 450" counts the checks that exist.
- Part B (found by the implementation audit): the reference DA arm on the `none` state restricted to F_filt
  (`sim-k185.R:994`) failed in most proteomics replicates and the `tryCatch` kept no error. The P0, P1 and P4b part-B FDP
  and power contrasts are missing, or rest on 7 to 20 of 100 replicates. Reported part only; no classification changes.
  Inferred, by reading only (no run): `subsetState` (`sim-k185.R:897-909`) subsets `protein_quant_table` by a logical vector and keeps
  the surviving row names (for example 1, 2, 5, ...), which are no longer the automatic ones; the DA then calls
  `tibble::column_to_rownames(theObject@protein_quant_table, ...)` (`R/func_prot_da_model_methods.R:287`), which stops
  when the data frame already has row names (the likely message: "`.data` must be a data frame without row names"). The
  replicates that returned a table are those where the subset left automatic row names. A fixed job needs
  `rownames(state@protein_quant_table) <- NULL` after the subset in `subsetState`, a new hash in a dated declaration
  section, and the six proteomics part-B jobs re-run (b-p0, b-p1, b-p4b, two each); the owner decides (VAL-54).
- Scoring: E5 and E8 compared tables by inferred column type; fixed (declaration section 7), re-scored once; every
  scored table identical between the two scoring runs (owner to be told; ticket section 10).
- Job splits used `round()` boundaries (1-133, 134-267, 268-400) where the declaration's section 6 wrote 1-134 etc.;
  the same replicates.

## Not checked

- P2 (held for KIT-187).
- Not summarised by the pinned scoring (rows only): (f) features lost, B2, the incomplete-feature recovery, the full
  canonical correlations, the comp/inc and residual-dimension strata, the oracle's own change and the uncorrected arm's
  departure, the per-arm per-run shares and residual-dimension counts, the q-value-method mismatches.
- Part B's proteomics DA contrasts (instrument defect).
- The release commit (A9 is stated at 147adde8).
- The run path's own declared lipid/metab DA kernel (E6 on L2 mirrors its input rule).
- Real data; other generators, residual distributions, normalisers beyond cyclic loess.
