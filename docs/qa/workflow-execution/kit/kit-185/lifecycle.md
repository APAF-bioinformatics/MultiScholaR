# KIT-185 lifecycle (L1), blind

Ticket revision 1.2. Branch `kit-185-ruviiic-missing`, base `v0.9-staging` at ffdfdfb76a9eacdc681f2098f75ac574e2cf8af4.
Runs, commits and core-hours: `runs.tsv`. Returned files: `l1/` (per job: blind rows, `checks-*.tsv`, `timing-*.tsv`,
`provenance-*.json`). The masked arms' estimands were not written (only completion, status, timing, mechanics counts and
the SHA-256 of the full rows, column `rows_sha256` of the timing tables); the `none` arm is written in full (G1).

## Observed

| Check / gate | Result | Where |
|---|---|---|
| E1 (a) raw and post-processed, every cell, replicate and mask | pass everywhere (largest 2.4e-7 at P4b bb2 raw, tolerance 1e-6 above 20 runs; at most 1.4e-11 at 20 runs or fewer) | `l1/*/checks-*.tsv` |
| E1 (b), one replicate per cell, none / VD 10 / MCAR 2 (BB 10 in the D5 cells) | pass in every cell (largest 3.2e-8, M1b VD 10); features in one replicate set excluded: 10 (L1b BB 10), 1 (L2, M2 VD 10) | same |
| E1-neg | pass in every cell **except P2** (C equals S to 3e-14 on every masked arm, all three replicates) | `l1/l1a-1/checks-l1a.tsv` |
| E2 | pass everywhere | same |
| E3 | P0 failed at L1a in the mechanics (the parameter list carried `ruv_mode` twice); after the fix (L1fix) identical; L2 3.6e-15, M2 3.6e-15 after row-centring (tolerance 1e-10) | `l1/fix-1/`, `l1/l1b-1/` |
| E4 | B 0.25, rho 1, C* 0.10 to 1e-12 | every checks file |
| E6 | P0 identical p-values (differentialAbundanceAnalysis on the declared state, its table replaced by C's); L2 identical (the declared DA's input rule mirrored, `R/utils_workflow_omics_declarations.R:801`; the run path's kernel needs a store manager and was not called) | `l1/fix-1/`, `l1/l1b-1/` |
| E7 | refusal recorded `refused` with `prot_ruv_insufficient_controls` / `omics_ruv_insufficient_controls` in every cell; the forced error recorded `crashed` (7 rows, P0 401 VD 0.1, C) and all 27 units completed | `l1/l1a-1/` |
| G1 | pass in every part-A cell: rho_min 0.926 (M1b) to 0.9998, B(S, none) on F^comp from -0.21 to -0.004 | checks `G1` |
| G2 | pass: no lipidomics or metabolomics DA failed or stopped, shipped filter and presence rule alike (L1, L2, M1, M2, L1b, M1b); the presence rule is therefore not taken | checks `G2` |
| Refusals by the route seen | L1 MCAR 10 %, M1 VD 50 % and MCAR 10 % (all three replicates); part B legacy `assay_not_corrected` at L2 and M2 VD 10 % (1 of 2) and VD 50 % (2 of 2), M1b one arm | rows |

Instrument defects (VAL-54, one run, L1fix): E3's parameters; the part-B comparator handed an empty assay group (L2, L1b)
and resolved the lipid route's controls by name instead of by position. No expected value, threshold or stop rule changed.

## G3, time per replicate (lifecycle, which adds the reference Z and E1's second compiled call to every mask)

| Cell | measured s / rep | estimate (10) | S s | C s |
|---|---|---|---|---|
| P0, P0r, P0lo, P0hi | 155-180 | 203.2 | 27 | 21 |
| P1 | 223 | 404.3 | 36 | 21 |
| P2 | 151 | 202.2 | 25 | 20 |
| P3 | 351 | 414.4 | 51 | 30 |
| P4 | 487 | 488.8 | 75 | 53 |
| P5 | 306 | 204.8 | 30 | 22 |
| L1, L2, M1, M2 | 48, 26, 123, 12 | 84.2, 53.3, 238.5, 36.1 | | |
| P4b, L1b, M1b | 360, 21, 68 | 279.6, 48.4, 136.5 | 82 | 26 |
| part B P0; L2 + M2 | 561; 59 | 647.4; 178.1 | | |
| part B P1, P4b, L1b + M1b | 400, 400, 61 | 780.5, 695.3, 171.1 | | |
| run-aligned P0 + L2 + M2; P1 | 60; 61 | 74.9; 102.0 | | |
| sensitivity P0 (5 variants); M2 | 441; 20 | 365.8; 66.7 | | |

C is below S in every cell (the plan assumed C = 2 S). Inferred: P5 is 49 % above its estimate with the lifecycle-only work
included (Z about 22 s and a second reference fit about 22 s per replicate); without it about 260 s, 27 % above. P4b 29 %,
sensitivity P0 21 %. Within the 40 % rule on the S-run configuration; P5 at 260 s plans 143.9 min per job in its 3 jobs,
under its 175-minute timeout. No re-split is needed and the ceilings are unchanged.

## The finding that stops the work before the declaration: P2 receives no missing value

Observed (mask rows of `l1/l1a-1/rows/a-P2-401-403-l1.tsv.gz`): after the route's filter, the realised missing share of
P2 is 0 on every mask; the filter removes 78 (VD 0.1 %) to 2,533 (VD 50 %) of 3,000 proteins, every incomplete one.

Read in the code: `updateMissingValueParameters` rounds the groupwise cutoff to three decimals
(`R/func_general_parameter_helpers.R:258`: (3 - 2) / 3 becomes 33.333 %), and the helper removes a group whose missing
share is **greater** than the cutoff (`R/func_prot_qc_filtering_helpers.R:238`, `perc_missing_per_group >
groupwise_percentage_cutoff`): one missing value of 3 is 33.3333 % > 33.333 %, so the group fails, and with two groups
`min_groups = 2` allows no failing group (cutoff 0 %). Inferred: with three runs per group the shipped defaults keep only
proteins observed in every run, so RUV-III-C never sees a missing value in this design (a rounding effect: the intent of
`min_reps_per_group = 2` is to allow one missing value of three). Groups of 4, 6, 8 and 10 are not affected (exact cutoffs).

Consequence for the declared measurement: E1-neg cannot pass on P2 (there is no difference to see), and P2's deciding
contrasts would measure only the loss of proteins to the filter (referenced out of the DA contrasts by the oracle on the
same rows), not the zero-filling. The ticket treats a failed E1-neg as an instrument defect, but this is not a mechanics
failure and no fix of the script can change it without changing the route's filter. By the brief's rule (unsure whether
instrument defect or result: it is a result), the work stops here and the owner decides.

## Later lifecycles (part-B re-run and P2)

L1pb (45291628697, df60efa9, VAL-98 (1)): blind, replicate 401 of P0, P1, P4b part B; the reference DA arm returned a
table in every cell and route (`qmethod` recorded). Outputs `l1/l1pb-1/` (committed after the close audit, C06 (1);
their SHA-256 equal those in `provenance-l1pb.json`). L1p2 (45314701113, 6779d56e, section 11): blind, P2 replicates
401-403; E1-neg 2.54, 1.85, 2.30; G1, E1 (a), E1 (b), E2, E7 pass. Outputs `l1/l1p2-1/`.

## Not checked

- The S-run timing itself (L1 measures with the lifecycle's extra work).
- The declared lipid/metab DA through the run path's own contrast kernel (E6 on L2 mirrors its input rule).
- Part B with more than two replicates; the refusal rates at R = 400.
