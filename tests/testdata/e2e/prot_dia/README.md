# prot_dia E2E Fixture

## Data Shape

| Property | Value |
|---|---|
| File | report.tsv |
| Format | DIA-NN report.tsv (long format) |
| Proteins | 5 (P00001–P00005) |
| Samples | 6 (WT_1–3, KO_1–3) |
| Precursors per protein | 2 |
| Total rows | 60 |
| FASTA | proteins.fasta (P00001–P00006 UniProt-style headers) |

## Groups

| Sample | Group |
|---|---|
| WT_1, WT_2, WT_3 | WT (control) |
| KO_1, KO_2, KO_3 | KO (treatment) |

## Expected Significant Proteins (KO_vs_WT)

| Protein | Gene | Direction | Approx log2FC |
|---|---|---|---|
| P00001 | GENE1 | UP in KO | ~1.05 |
| P00002 | GENE2 | UP in KO | ~1.02 |

P00003–P00005 are not expected to be significant (<0.1 log2FC between groups).

## Import Function

`importDIANNData()` in `R/func_prot_import_readers.R`

Required columns consumed include `Protein.Group`, `Stripped.Sequence`, `Run`,
`Precursor.Quantity`, `Precursor.Normalised`, `Q.Value`, `Global.Q.Value`,
`PG.Q.Value`, `Global.PG.Q.Value`, and `Proteotypic`. The synthetic fixture uses
the same valid `0.001` confidence value for run-specific and global gates and
marks every invented precursor as proteotypic.
