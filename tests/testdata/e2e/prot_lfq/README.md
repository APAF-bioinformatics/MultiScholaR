# prot_lfq E2E Fixture

## Data Shape

| Property | Value |
|---|---|
| File | seed_combined_protein.tsv |
| Format | FragPipe combined_protein.tsv (wide format) |
| Proteins | 5 (P00001–P00005) |
| Samples | 4 (WT_1, WT_2, KO_1, KO_2) |
| Total rows | 5 |

## Groups

| Sample | Group |
|---|---|
| WT_1, WT_2 | WT (control) |
| KO_1, KO_2 | KO (treatment) |

## Column Format

Intensity columns follow FragPipe combined_protein.tsv convention:
`{sample_name} MaxLFQ Intensity`

The `importFragPipeData()` function detects columns ending with `MaxLFQ Intensity`,
strips the suffix to extract sample names, and pivots to long format.

## Expected Significant Proteins (KO_vs_WT)

| Protein | Gene | Direction | Approx log2FC |
|---|---|---|---|
| P00001 | GENE1 | UP in KO | ~1.05 |
| P00002 | GENE2 | UP in KO | ~1.01 |

P00003–P00005 are not expected to be significant (<0.1 log2FC between groups).

## Import Function

`importFragPipeData()` in `R/func_prot_import_readers.R`

Required columns: `Protein ID` (used as protein identifier), `{sample} MaxLFQ Intensity` columns
