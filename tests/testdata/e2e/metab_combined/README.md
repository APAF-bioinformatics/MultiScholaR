# metab_combined — Multi-Assay Metabolomics Fixture

Minimal synthetic seed data for the combined multi-assay metabolomics workflow lane.
Exercises the `createMetaboliteAssayData()` path with three simultaneous assays.

## Files

| File | Description |
|------|-------------|
| `seed_lcms_pos.csv` | MS-DIAL export, positive ion mode (8 features) |
| `seed_lcms_neg.csv` | MS-DIAL export, negative ion mode (8 features) |
| `seed_gcms.csv` | MS-DIAL export, GC-MS format (8 features) |
| `design.tsv` | Sample-to-group mapping (ctrl vs treat, 2 replicates each) |

These files are identical to their counterparts in `metab_lc/` and `metab_gc/`, combined here to exercise the multi-assay import path (`assayList` with LCMS_Pos, LCMS_Neg, GCMS).

## Expected Assay Names

When imported, assays should be named: `LCMS_Pos`, `LCMS_Neg`, `GCMS`.

## ITSD Features (per assay)

- LCMS_Pos: `ITSD_Caffeine_d9`, `ITSD_Leucine_d10`
- LCMS_Neg: `ITSD_Succinic_acid_d4`, `ITSD_Glutamic_acid_d5`
- GCMS: `ITSD_Ribitol_d4`, `ITSD_Myristic_acid_d27`

All ITSDs match the `_d[0-9]+$` regex from `getLipidomicsColumnDefaults("msdial")`.

## Expected Behavior

- Each assay independently has ~half features up in ctrl, half up in treat
- ITSD features stable across groups for normalization testing
- Combined object exercises cross-assay DA orchestration via `assay_names` iteration
