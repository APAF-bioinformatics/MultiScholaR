# Lipidomics Canonical E2E Fixture

Synthetic LipidSearch seed data for the lipidomics canonical workflow lane.
The lane is intentionally multi-assay because lipidomics state is represented
as a named assay list throughout import, design, QC, normalization, DA, and
summary/report generation.

## Files

| File | Description |
|---|---|
| `lipidsearch_lcms_pos.txt` | Primary LCMS positive-mode LipidSearch assay |
| `lipidsearch_lcms_neg.txt` | Secondary LCMS negative-mode LipidSearch assay |
| `lipidsearch_gcms.txt` | GCMS-named LipidSearch smoke fixture for non-LC assay naming |
| `design.tsv` | Six-sample `WT` versus `KO` design used by all assays |

## Import Shape

All three data files are tab-separated LipidSearch-style wide matrices with:

- identifier and annotation columns: `LipidName`, `LipidClass`, `FattyAcid`, `IonType`, `BaseRt`, `CalcMz`, and `Grade`;
- six quantitative sample columns: `WT_1`, `WT_2`, `WT_3`, `KO_1`, `KO_2`, and `KO_3`;
- two strong synthetic group effects so DA output is deterministic enough for E2E artifact assertions.

The canonical manifest lane imports `lipidsearch_lcms_pos.txt` as `LCMS_Pos`
and `lipidsearch_lcms_neg.txt` as `LCMS_Neg`. The GCMS fixture is used by the
E2E-009 import smoke to prove the browser lane does not assume positive/negative
LC-only assay names.
