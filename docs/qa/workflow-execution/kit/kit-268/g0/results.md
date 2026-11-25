# KIT-268 G0 results (run 45708195790, base 4d588ba3 + docs; tables verbatim under run-45708195790/scripts/5/out/)

Positive control (tied_pvalues lipid through completeOmicDaRun): fallback comparisons LCMS_Pos, LCMS_Neg, GCMS x groupKO-groupWT,
groupRES-groupWT (6 keys, each an exact match on `contrasts`; rule result = friendly name). Non-empty: the instrument is live.

## O-3 table: stored `da` steps that would gain a notification (OBSERVED: results$qvalue_warnings non-empty at the lane's one completeOmicDaRun call)
| file | step | line | lane | assay | comparisons that fell back |
|---|---|---|---|---|---|
| tests/testdata/gui-entry/lipid_canonical.json | da | 2402 (key 2401) | lipid_canonical | LCMS_Pos | groupKO-groupWT (KO_vs_WT) |
| tests/testdata/gui-entry/lipid_canonical.json | da | 2402 | lipid_canonical | LCMS_Neg | groupKO-groupWT (KO_vs_WT) |
| tests/testdata/gui-entry/metab_lc.json | da | 2392 | metab_lc | none | none (call observed, 0 keys) |
| tests/testdata/gui-entry/metab_gc.json | da | 2316 | metab_gc | none | none (call observed, 0 keys) |
| tests/testdata/gui-entry/metab_combined.json | da | 2392 | metab_combined | none | none (call observed, 0 keys) |

So one stored step moves (lipid_canonical, `da`, one added warning note); the three metabolomics oracles do not.
`test-omic-summary-module.R`: no call of completeOmicDaRun by any lipid or metab capture (4 calls, all transcriptomics, 0 keys): nothing moves.
`test-omic-da-module.R`: 15 calls (lipid_main, lipid_no_friendly, metab_main, metab_no_friendly, the managed lipid/metab/trx, the W3/W4 volcano tests, trx), all 0 keys: no W1 capture moves.

## Mapping observation
Real route (lipid_canonical): keys are `groupKO-groupWT` under each assay; the lane's contrasts table has contrasts `groupKO-groupWT`, friendly_names `KO_vs_WT`: exact match on `contrasts` on both assays; the ticket's rule gives KO_vs_WT (in friendly_names). Control (stub): the same, 6 of 6. No key carried an "=friendly" suffix. Stop rule 3 not triggered on the routes observed (the proteomics mapping is T2b's).

## Limits of the instrument
`step` is NA in every row (the harness's loop variable was not found in the frames): steps are inferred from the lane having one traced call each (the `da` step; each lane's stored `da` step holds the completion note). Lane liveness: each of the four lanes shows exactly one call. The stub's tied_pvalues keys are the same assay -> comparison shape as the real kernel's.
Test counts in `*-tests.tsv` are context: the two GUI-entry files 35/0 and 35 x 3/0; da-module and summary show their known stored-reference failures (not judged here).
