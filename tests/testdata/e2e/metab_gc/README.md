# metab_gc — GC-MS Metabolomics Fixture

Minimal synthetic seed data for the GC-MS (gas chromatography) metabolomics workflow lane.

## Files

| File | Description |
|------|-------------|
| `seed_gcms.csv` | MS-DIAL export, GC-MS format (8 features, 4 samples) |
| `design.tsv` | Sample-to-group mapping (ctrl vs treat, 2 replicates each) |

## Format

GC-MS seed file uses MS-DIAL alignment export CSV format without the m/z column (GC-MS uses RI not m/z for identification):

- `Alignment ID` — numeric feature identifier
- `Average Rt(min)` — retention time
- `Name` — metabolite name or annotation
- `ctrl_1`, `ctrl_2`, `treat_1`, `treat_2` — sample peak intensities

Note: No `Average Mz` or `Adduct type` columns — GC-MS features lack this information in MS-DIAL exports.

## ITSD Features

Internal standards are identifiable by the `_d[0-9]+$` suffix pattern:

- `ITSD_Ribitol_d4`
- `ITSD_Myristic_acid_d27`

ITSD intensities are roughly equal across all 4 samples to support normalization testing.

## Expected Behavior

- Features 1–3: ~2x higher in ctrl group (Alanine, Valine, Glucose)
- Features 4–6: ~2x higher in treat group (Palmitic_acid, Stearic_acid, Metab_GC_001)
- Features 7–8: ITSD (stable across groups)
- Expected contrast: `treat_vs_ctrl`
