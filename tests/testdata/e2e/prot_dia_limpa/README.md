# prot_dia_limpa E2E Fixture

## Data Shape

| Property | Value |
|---|---|
| File | seed_report.tsv |
| Format | DIA-NN report.tsv (long format) with MNAR missing values |
| Proteins | 6 (P00001–P00006) |
| Samples | 6 (WT_1–3, KO_1–3) |
| Precursors per protein | 2 |
| Total rows | 72 |
| Missing values (NA) | ~22 rows (~30%) |

## Groups

| Sample | Group |
|---|---|
| WT_1, WT_2, WT_3 | WT (control) |
| KO_1, KO_2, KO_3 | KO (treatment) |

## Missing Value Pattern (MNAR)

Missing values are concentrated in low-abundance proteins in specific conditions,
simulating "Missing Not At Random" (MNAR) data typical of DIA proteomics:

| Protein | Missing In | Reason |
|---|---|---|
| P00005 | WT_2, WT_3, all WT_3 precursors | Low abundance in WT, near LOD |
| P00006 | All WT samples, KO_3 | MNAR: present in KO, absent in WT |
| P00004 | WT_1 (1 precursor), WT_2 (both), KO_3 (1 precursor) | Sporadic low detection |

## Expected Significant Proteins (KO_vs_WT)

| Protein | Gene | Direction | Approx log2FC |
|---|---|---|---|
| P00001 | GENE1 | UP in KO | ~1.05 |
| P00002 | GENE2 | UP in KO | ~1.02 |
| P00005 | GENE5 | UP in KO | ~1.5 (MNAR, requires imputation) |
| P00006 | GENE6 | UP in KO | ~Inf (MNAR, requires imputation) |

## Import Function

`importDIANNData()` in `R/func_prot_import_readers.R`, followed by limpa imputation.

Required columns consumed: `Protein.Group`, `Stripped.Sequence`, `Run`, `Precursor.Quantity`, `Q.Value`, `PG.Q.Value`
