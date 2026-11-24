# KIT-188 lifecycle (L1 and L1p), blind

Ticket revision 1.6. Base: `v0.9-staging` at 98df468b (before KIT-190: the positional selection, so on this base the LA
request runs what the ticket calls LP). Scripts at 049c51b4 (L1) and ef66c8f8 (the L1 rerun and L1p); `git diff
98df468b ef66c8f8 --stat` lists only this directory. Runs and costs: `runs.tsv`. Every estimand column of the written
rows is NA (blind; checked on `l1/l/rows/l-L12-401-403-l1.tsv.gz`: 0 of 723 rows carry `fdp`, `null_p05`, `power`,
`s5`, `err_max`, `calls`); each unit's unblinded SHA-256 is in its `timing-*.tsv`.

## Observed

| Item | L1 (L12, L64, M12, M168; 401-403) | L1p (P0, P2: 401-403; P1, P16, P3, P4b: 401) |
|---|---|---|
| Units | 36 of 36 after one instrument rerun (below) | 262 of 262 |
| E1 (declared operations and variants, journeys that ran) | 120 of 120 pass (`l1/*/checks-*.tsv`) | 142 of 142 pass |
| E2 (Skip, gate off, against design `none`, gate on: identical DA) | 5 of 5 (L) and 2 of 2 (M) pass | **fails in all 10 data sets**: largest difference over logFC, t, p, q 0.029 to 0.39 |
| E3 (defaults read from the configuration and tab builders) | pass (lipid, metab) | pass (proteomics) |
| E4 (estimand functions, constructed case) | pass | pass |
| E6 (oracle construction, at most 1e-12) | 135 of 135 pass (147 with the rerun) | 262 of 262 pass |
| E7 (refusal; deliberate crash recorded, job continues) | refusal pass (L64, M12, M168, L12 rerun); crash pass (M job, L12 rerun; in the L job it read FALSE only because the three L12 units stopped) | refusal pass in all six; crash pass (P0) |
| **E9** (shipped selection) | **385 of 385 pass**: every legacy Automatic and Manual journey's controls are the first round(N_v x pct / 100) rows among the rows with at least two values, in table order, in every assay; the per-percentage counts match for all 20 percentages; the control method called directly gives the positional set at every percentage 1 to 20; the RUV input rows keep the written order | **142 of 142 pass**: LA's controls equal the helper's ANOVA rule on `group` with `qvalue` (and differ from the BH rule) |
| E8 (reported, never re-run) | 90 of 90 within 1e-8 (12 runs) or 1e-6 (above 20) | 46 of 46 |
| E5r (reported jobs' LA and OR rows equal the deciding jobs') | 7 of 7 pass where both units ran without the deliberate fault; L12 401-403 (units stopped) and the fault units L64 401, M168 401 differ by construction | not applicable |

**G1, route crashes (a finding for the owner, before S).** The shipped lipidomics and metabolomics DA stops with
`prior.weights contain NA values` on many masked data sets, for every arm, **including the oracle (no RUV, normalisation
`none`)**, so the stop is in the DA step on sparse data, not in RUV. Journeys stopped (deciding jobs): L12 VD 20 % 30 of
72; M12 VD 20 % 57 of 72, VD 5 % 6 of 72; M168 VD 20 % 13 of 48; mask BA: L64 3 of 48, M168 13 of 48; none with mask
`none`, none in L64's VD masks; the reported and M-unit jobs the same pattern. At more than 5 % of a condition's data
sets these conditions classify as "no result" (5.7), and the oracle reference is missing for S5, F1 and F2 there.
Proteomics: no journey stopped. Other statuses (counted, as 5.7 states): legacy Automatic `dropped` and `uncorrected`
under VD 20 % and BA; legacy Manual `refused` (`omics_ruv_insufficient_controls`) in L64 VD 20 %; LCm `not_run` when LA
gives no k or controls in an assay (24 in M, 16 in L12, 8 in L64) and `refused` (`ruv_k_above_limit`).

**Instrument defect (VAL-54), one rerun.** L1-L's three L12 reported units stopped with `subscript out of bounds` (no
result). Fix at ef66c8f8: M-df scores only the features the truth covers, by-name lookups guarded, rows written as each
unit completes, crashed journeys recorded as G1, not E1. Rerun 45293899542 (1 of 3 allowed): 3 of 3 units, every check
passes. No expected value, edge, comparison or stop rule changed.

**Time.** Core-seconds per data set (sequential, one core per child): L12 deciding 1,166 (18 conditions), reported
1,529, munit 960; L64 deciding 2,072, reported 3,012, ba 616; M12 deciding 705 (18 conditions), reported 797, munit 547;
M168 deciding 2,577, reported 3,669, ba 767. Part P: deciding P0 2,485 (12 base conditions; 1,203 for N2, E2), P1 6,859,
P16 7,989, P2 2,088, P3 9,203, P4b 6,501; LM with OR 447 to 1,782; P4b BA 1,951. L1's include the lifecycle-only checks
(E2, E9's direct sweep) and journeys shortened by the DA stops.

## Inferred (not observed)

- G3 (Part L, indicative; the formal G3 uses L1r's repaired LA): measured over estimated per data set L12 deciding 1.25,
  L64 1.33, M12 1.51, M168 1.38; reported L64 1.08, M12 1.43, M168 1.10, L12 1.37; munit L12 1.40, M12 1.60; ba L64 1.18,
  M168 1.23. Re-split at the same R (plan at most 170 minutes per job) gives about 1,220 core-hours expected and 2,260
  ceiling with D6, D8 and the 25 % reserve, above the accepted 978.9 / 1,809.3 (VAL-84): by G3 the owner decides.
- Part P (D4b, VAL-84 (8)): recomputed on the measured times at R = 400 deciding (D6 included), LM with OR and P4b BA at
  R = 100, the same job rule: about 4,531 core-hours expected (plus 21 for L1p as run) and 8,422 ceiling with the 25 %
  reserve, 386 S jobs; inside the indicative bracket (3,356 to 6,557; 6,236 to 12,131). P1 is above both indicative
  scalings (808 against 604 / 483 deciding core-hours). Measured on one data set in four of six cells.
- E2 (proteomics): the legacy Skip DA and the declared `none` DA are different tests (18.1: the legacy `treat` without
  trend or robust), so RUV off has two figures in proteomics; which one 5.3's NO measures is for the owner.

## Not checked

- Where in limma or the package `prior.weights contain NA values` arises (the DA step; not traced, no R here).
- LP's job type (identical to LA on this base) and job D; the closed-form branch recorder (`cfBranch`) ran only in L1p
  (proteomics) and not on lipid or metabolomics data; E10 and the repaired LA's time (L1r, after KIT-190).
- `score-k188.R` has not run.
