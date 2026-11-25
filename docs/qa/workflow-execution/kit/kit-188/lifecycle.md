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

## Small lifecycle check of the owed script changes (r1.13; runs 45334020028, 45334023492)

At aee81861 = staging ea93d948 (KIT-177 merged; KIT-190, KIT-196, KIT-201 not merged) plus the scripts at bce0d986.
Observed: both jobs ran to completion (2 of 2 and 24 of 24 units). Exercised: the per-assay `no_da` path through a
deliberate failure (E7_noda pass) and through real ones, so journeys that L1 recorded as `ok` while one assay's DA had
failed are now `no_da` for that family (L12 at VD 20 %: LA, NO, OR in N1 and E2); the crash path (E7_crash pass, L and P);
the refusal path (P2); E4b (tie arithmetic, constructed); the tie count on 23 of 23 proteomics legacy journeys, the replica
equal to the package's selection in every one; the design-step transform taken from the commit (`log2(x); 0 missing` in all
12 checks and in the provenance); the Skip/`none` gate through its failing branch (largest difference 0.20; expected before
KIT-201); K172c in both classes, all five arms `ok` (values blinded). E1, E2 (lipid), E3, E4, E6, E8, E9 pass.
Not exercised: the gate's passing branch and `K188_NON=1` (need KIT-201); E10 (needs KIT-190); LP; the reported, M-unit
and BA job types of Part L (unchanged code paths since L1, apart from `no_da`); metabolomics cells (only through K172c).

## Script changes for revisions 1.16 to 1.21 (written, not run; no R was available)

Ticket revision 1.21 (sections 35 to 40). Scripts as written on `kit-188-ruv-default-calibration`, started from origin tip
e2d18a79 and then merged with origin/v0.9-staging 176f9d3ff (no conflict). Nothing was run: no lifecycle run, no request file,
`runs.tsv` untouched. The checks made here are parse-level only (bracket balance by a Python scan of both files; every package
function the new code calls was found at origin/v0.9-staging 176f9d3ff by `git grep`). **Nothing below has executed, and the
R-marked predictions (35.6 r1.20) cannot be written until KIT-242, KIT-260 and KIT-255 are on staging.**

Line numbers are those of `sim-k188.R` (S) and `score-k188.R` (C) at the commit that carries this section.

| Requirement | Ticket | Where |
|---|---|---|
| `Name` annotation column holding every row's name, made the annotation column; class column and slot regex unchanged | 35.2 (1), 5.2 | S:345-380 (`omicObject`, `annotation_column` S:363) |
| E12 per assay (a) no `IS_` row left, (b) record count = standards at `qc_complete` (0 when absent), (c) a non-standard value changed; `itsd_max_change`, class `itsd_not_applied:<assay>:<letters>` | 35.2 (2), 7.4 E12 | S:483-502 (`checkItsdActed`), called S:659; status/class S:660-668 |
| E12 format-default pattern assertion, once per job, defaults read from the bound commit, through `stringr::str_detect` | 35.2 (2), R16-P16, 40.3 (8) | S:1873-1890, called S:2124 |
| E13 per assay, matched by feature id and run; `norm_max_change`; class `norm_not_applied:<assay>`; not for OR or method `none` | 35.2 (2), 7.4 E13, R16-P17 | S:504-519, called S:678-685 |
| `itsd_ok` and `norm_ok` per-assay columns, never overwritten; `status` stays `instrument`, a later error is appended to `class` | 35.2 (2), R16-P07 | S:554-620, S:748-762 |
| E14 per-assay `ruv_status` (`not_reached`, `failed`, `refused`, `dropped`, `uncorrected` by name at tolerance 0, `corrected`; declared designs: the factors' status), `ruv_changed`, `ruv_k_found`, `ruv_n_controls`, `ruv_no_k_reason` (status counts and the first `error_reason`) | 35.2 (3), R16-P05, R16-P06 | S:520-530 (`sameByName`), S:554-620 (`journeyAssayColumns`) |
| Comparison moved out of the `status == "ok"` guard; journey status keeps its precedence | 35.2 (3) | S:776-781 |
| Product fields `ruv_k_sweep`, `ruv_k_used`, `ruv_k_cap`, `ruv_k_samples`, `ruv_k_notice`; record fields equal the product's (`ruv_record_match`), `ruv_record_applied`; the product's `ruv_status` beside the script's (`ruv_status_product`, `ruv_status_agree`) | 35.2 (3) r1.20 | S:593-613, record read S:697-698 |
| E14 count (a status per assay for every journey that reached RUV) and the script-against-product agreement, an `identical()`-strictness difference passing with a note for KIT-242 | 7.4 E14, R20-P14, 40.3 (7) | S:1566-1575 |
| DA outcome per assay: `da_status` ok / failed / no_da / not_reached, `da_class`, `da_message` | 35.2 (3) | S:582-584 (`journeyAssayColumns`), S:725-732 (manifest errors) |
| `corr_removed_one_assay` (assay:sample listed where the correlation filter left a sample out of one assay only) | 35.2 (3) r1.21 R20-P11, 14 | S:531-545 (`oneAssayRemovals`), S:702 |
| The two different `best_k`: the trace's is the sweep's k (`k_sweep`), the product's is the k used (`k`); Part M, LCm, E8, M-amp, M-df take the k used | 40.3 (2), 5.8 R20-P04 | S:790-801 (`legacyChoice`); users unchanged (`choiceOf`) |
| M-warn's first signal is `k_sweep > k_cap`; the other four signals; sensitivity and specificity with Clopper-Pearson | 5.8, R20-P04 | C:411-440; `logfc_max` recorded S:1005 |
| KIT-217's two settings passed explicitly from E3's read values | 40.3 (3), 7.4 E3 | S:642-644; read S:396 (`routeDefaults`) |
| E3 records the cyclic-loess minimum's state, does not stop on it; marker tokens; the import pattern as KIT-260 ships it | 7.4 E3, R20-P06 | S:433-452, `ITSD_DEFAULT_PATTERN` S:432 |
| `quiet()` records warning and classed-message classes per journey and assay (`warning_classes`) | 35.2 (2) R20-P08, 40.3 (4) | S:101-122 |
| L24: cell definition (c = 11, 12 per group, no N2 or E2), deciding job type only, R = 100 | 5.2, 38.5, 38.3 | S:78-80, `jobApplies` S:1491-1500 |
| Plants `K188_PLANT`: `itsd`, `itsdpos`, `emptysample`, `ruvskip`, `ruvdrop`, `daerror`, `normplant`, `exitplant` as job `l1c` (one child each) | 35.2 (5), 35.6 r1.18, r1.21 | S:1892-1925 (list, predictions), `runPlantUnit` S:1963-2005, units S:2082-2100; `itsd`/`itsdpos` in `omicObject` S:345-380; `normplant` S:676; `daerror` S:709 |
| `exitplant`: the child `quit(status = 3)` before writing rows; the parent's unit row checked | 35.6 R20-P10, 7.4 E15 | S:2002, parent check in `runChild` |
| Plants act on the generator's input, the request's argument, the matrix E13 reads, or the child's exit, never on the route's state | 40.3 (6) | S:1892-1896 (comment); `normplant` S:676 copies only the comparison matrix |
| E7: an invalid typed k of -1 (a stop by rule); the r1.20 plant kept as an observation | 7.4 E7, R20-P07 | S:1849-1867 |
| R-marked predictions: read from a TSV written into lifecycle.md (`K188_PREDICTIONS`), compared, and the run stops after writing its outputs when one is not met; a run sent without the file stops at the start | 35.6 r1.20, 40.3 (8) | S:1930-1960, S:2086-2093, S:2276-2283 |
| LP on the bound commit with `MultiScholaR.ruv_control_grouping = "dropped"` and `MultiScholaR.ruv_k_bound = "none"`; LAs with `MultiScholaR.ruv_feature_levels = "shifted"`; set inside the journey, restored on exit, asserted at each journey's start (a leak stops the unit); each token asserted equal to the marker table's `earlier` | 5.3, 34, 40.3 (5), R20-P02 | S:137-160, S:622-630, jobs S:1478 and S:1510 |
| L1r check: LP's controls equal E9's positional rule (with LP's markers set for the direct control call) and LP's k equals the uncapped sweep's k | 40.3 (10), 5.3 | S:1105-1180 (`checkE9L`), called S:1562 |
| E10 (rule (d), the bound, the record fields, K172c's k at most 2) | 7.4 E10 r1.20 | S:1268-1315, called S:1563 |
| E11 (LAs against LA: same RUV input, controls and k; per-feature constant at most 1e-9; logFC 1e-10) | 34, 7.4 E11 | S:1316-1345, called S:1619 |
| Scorer: stop on `itsd_ok` / `norm_ok` FALSE or `instrument` rows | 38.3 (i), 35.2 | C:61-65 |
| Scorer: unit rows expanded to every arm and condition of their job as `crashed`; declared range (deciding 1 to 400, others and L24 1 to 100) checked, a gap stops | 38.3 (ii), 35.2, R16-P02 | C:66-100 |
| Scorer: gate status (crashed, refused, no_da, dropped, unit crash), DA-table flag, journeys without a DA table out of every estimand including S3 | 38.2 (2), (3), 38.3 (iii), (iv) | C:107-146, C:196-215 |
| Scorer: best and worst case, `sensitive` flag | 38.2 (4), 38.3 (vi) | C:160-175, C:199-203 |
| Scorer: paired comparisons on the data sets where the arms compared and the oracle returned a DA table; `n_paired_dropped` | 38.2 (5), 38.3 (vii) | C:148-152, C:207, C:218-250, K1's efficiency clause C:319-336 |
| Scorer: phi with both denominators and its one-sided 97.5 % Clopper-Pearson limit; a stop at RUV is not an event; GA's event is `refused` factors | 38.2 (6), 40.1 draft, 40.3 (9) | C:269-286, C:303-318 |
| Scorer: K1's phi clause (limit below 0.05 in every deciding condition); K2 states GA's phi; VAL-248 label; the outcome written per denominator until the owner confirms | 38.2 (7), (8), 38.9 (1) | C:303-350 |
| Scorer: marking of conditions where LA equalled RUV off in more than 5 % of data sets (`lan0_share`, `lan0_mark`) | 38.2 (9), 38.3 (ix) | C:269-286 |
| Scorer: classifications at R = 100 labelled as unable to meet S5 (the upper limit at zero events is computed from n, not typed) | 38.3, R20-P13 | C:255-263 |
| Scorer: L24 and the reported arms enter no family, verdict or outcome | 38.5 | C:25-28 (`OMIC_GROUP`), C:216 and the verdict set C:295 |
| Scorer: cap-binding share, `k_sweep`, k used and cap, share of data sets with an assay not corrected, share with `da_status` failed | 5.5, 35.2 | C:287-292, phi.tsv columns |
| Scorer: per-assay RUV and DA status counts (a second file, `statuses-assay.tsv`, so `statuses.tsv` keeps one count per journey) | 35.2 (`score-k188.R`) | C:352-358 |
| Scorer: LP-LA (M-repair) and LAs-LA tables paired per data set with intervals | 5.8, 34, R20-P02 | C:388-410 |

## What waits for KIT-242, KIT-260 and KIT-255 (on the bound commit, not yet on staging at 176f9d3ff)

- **KIT-260**: `ITSD_DEFAULT_PATTERN` (S:432) is the string 5.2 states for VAL-241 OD-1 (c), cited from the ticket; at 176f9d3ff the
  import still pre-fills `w/o MS2:`, so E3 and the E12 pattern assertion read as a route change there. E12's pattern check reads
  `getLipidomicsColumnDefaults` and `getMetabolomicsColumnDefaults` at run time, so it needs no edit. The per-assay record of the
  rows used as standards (VAL-241 OD-2 (b)) is only recorded by name: `itsd_record_fields` lists the names found in
  `@args$ITSDNormalization` and in the analysis record's ITSD step; E12 (b) still compares the existing `itsd_counts_per_assay`.
- **KIT-242 (A and B)**: the predictions marked **R** in 35.6 (r1.20) are not written here and must be: the per-assay ITSD and
  normalisation records, their notice classes and report sentence (`itsd`, `itsdpos`, `emptysample`), the trace status string and
  the per-assay `ruv_applied` and reason fields (`ruvskip`, `ruvdrop`), whether "fewer than 5 controls at every percentage" is a
  genuine no k (Q6.2 (c)) and not a stop (`ruvskip`, `ruvdrop`), what `omicsCyclicLoessMinimumActive()` is at the bound commit,
  whether the minimum is evaluated with the empty column (`emptysample`), and the class of E7's two stops. Re-read at the bound
  commit, write them as a TSV (`key` = `<plant>/<field>`, `expected`; fields as `flattenJourney` names them, e.g.
  `ruv_no_k_reason.LCMS_Pos`; `re:` for a Perl regular expression, `range:a:b`, `NA`; keys `E7_refusal/status`, `E7_refusal/class`,
  `E7_invalid_k/status`, `E7_invalid_k/class` for E7) into this file, and pass its absolute path in `K188_PREDICTIONS`. The L1c job
  stops at its start without it, and stops at its end if a written prediction is not met (S:2276).
  The field names the script reads from the product and the record (`k_sweep`, `k_cap`, `k_samples`, `k_notice`, `best_k`,
  `ruv_status`; record `k_sweep_per_assay`, `k_cap_per_assay`, `k_per_assay`, `ruv_status_per_assay`) are those of 176f9d3ff; re-read
  them at the bound commit (a difference stops L1r, 7.4 E10, E14).
- **KIT-255**: `corr_removed_one_assay` is computed from the matrices before and after the correlation step (S:531-545); whether
  KIT-255 leaves such a sample in the design as missing in one assay (so it shows as an all-missing column) or removes it is read
  at the bound commit, and the function's reading of "removed" (absent column, or no finite value) is checked against it.
- Nothing else in the scripts depends on the three tickets.

## For the author of the declaration and the auditor

- Parse level only: the new code has never been parsed by R. Expect instrument defects in a first lifecycle run (VAL-54 allows
  fixing script mechanics and re-sending, at most three times per item).
- The scorer needs `K188_SCORE_R` only for E16 (a constructed row set); a counted run leaves it unset.
