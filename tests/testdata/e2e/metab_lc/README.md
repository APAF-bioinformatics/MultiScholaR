# metab_lc — LC-MS Metabolomics Fixture

Minimal synthetic seed data for the LC-MS (liquid chromatography) metabolomics workflow lane.

## Files

| File | Description |
|------|-------------|
| `seed_lcms_pos.csv` | MS-DIAL export, positive ion mode (8 features, 4 samples) |
| `seed_lcms_neg.csv` | MS-DIAL export, negative ion mode (8 features, 4 samples) |
| `design.tsv` | Sample-to-group mapping (ctrl vs treat, 2 replicates each) |

## Format

Both seed files use MS-DIAL alignment export CSV format:

- `Alignment ID` — numeric feature identifier
- `Average Rt(min)` — retention time
- `Average Mz` — precursor m/z
- `Name` — metabolite name or annotation
- `Adduct type` — ion adduct
- `ctrl_1`, `ctrl_2`, `treat_1`, `treat_2` — sample peak intensities

## ITSD Features

Internal standards are identifiable by the `_d[0-9]+$` suffix pattern (matching `importLipidMSDIALData()` is_pattern):

- Positive mode: `ITSD_Caffeine_d9`, `ITSD_Leucine_d10`
- Negative mode: `ITSD_Succinic_acid_d4`, `ITSD_Glutamic_acid_d5`

ITSD intensities are roughly equal across all 4 samples to support normalization testing.

## Expected Behavior

- Features 1–3 in each file: ~2x higher in ctrl group
- Features 4–6 in each file: ~2x higher in treat group
- Features 7–8 in each file: ITSD (stable across groups)
- Expected contrast: `treat_vs_ctrl`
