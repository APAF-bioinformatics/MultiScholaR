# Lipidomics Canonical E2E Fixture

Synthetic seed data for the lipidomics canonical workflow lane.

## Files

| File | Description |
|------|-------------|
| `seed_msdial.csv` | MS-DIAL height matrix export (10 lipids × 4 samples) |
| `design.tsv` | Sample-to-group mapping (2 groups, 2 replicates each) |
| `lipidsearch_results.txt` | LipidSearch format fixture (pre-existing) |

## seed_msdial.csv

MS-DIAL wide-format export with the following structure:

- **Annotation columns**: `Alignment ID`, `Average Rt(min)`, `Average Mz`, `Name`, `Adduct type`, `Ontology`
- **Sample columns**: `ctrl_1`, `ctrl_2`, `treat_1`, `treat_2`
- **10 lipid species** across 4 classes (PC, PE, TG, SM)
- **Missing values**: rows 3, 4, 7, 10 contain zeros to exercise QC filtering

### Lipid species

| Class | Species | Adduct |
|-------|---------|--------|
| PC | PC 34:1, PC 36:2 | [M+H]+ |
| PE | PE 36:2, PE 38:3 | [M+H]+ |
| TG | TG 52:3, TG 54:3, TG 56:6 | [M+NH4]+ |
| SM | SM 34:1, SM 36:1, SM 42:2 | [M+H]+ |

Lipid names follow `CLASS XX:Y` nomenclature (total carbons : total double bonds) for
correct class detection by `func_lipid_import_detection.R`.

## design.tsv

| Sample | Group |
|--------|-------|
| ctrl_1 | ctrl |
| ctrl_2 | ctrl |
| treat_1 | treat |
| treat_2 | treat |

## Usage

```r
msdial_result <- importLipidMSDIALData(
  "tests/testdata/e2e/lipid_canonical/seed_msdial.csv"
)
design <- read.delim("tests/testdata/e2e/lipid_canonical/design.tsv")
```
