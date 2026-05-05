# prot_lfq E2E Fixture

## Data Shape

| Property | Value |
|---|---|
| Files | proteinGroups.txt; seed_combined_protein.tsv |
| Formats | MaxQuant proteinGroups.txt; FragPipe combined protein |
| Proteins | 5 (P00001–P00005) |
| Samples | 6 (WT_1–3, KO_1–3) |
| Total rows | 5 |
| FASTA | proteins.fasta (P00001–P00006 UniProt-style headers) |

## Groups

| Sample | Group |
|---|---|
| WT_1, WT_2, WT_3 | WT (control) |
| KO_1, KO_2, KO_3 | KO (treatment) |

## Column Format

MaxQuant intensity columns follow proteinGroups convention:
`LFQ.intensity.{sample_name}`

FragPipe intensity columns follow combined-protein convention:
`{sample_name} MaxLFQ Intensity`

The `importMaxQuantData()` function detects `LFQ.intensity.*` columns,
extracts sample names, filters contaminants/reverse hits, and pivots to long format.

## Expected Significant Proteins (KO_vs_WT)

| Protein | Gene | Direction | Approx log2FC |
|---|---|---|---|
| P00001 | GENE1 | UP in KO | ~1.05 |
| P00002 | GENE2 | UP in KO | ~1.01 |

P00003–P00005 are not expected to be significant (<0.1 log2FC between groups).

## Import Function

`importMaxQuantData()` and `importFragPipeData()` in `R/func_prot_import_readers.R`

Required columns: `Protein.IDs`, `Gene.names`, `LFQ.intensity.*` columns
for MaxQuant; `Protein ID`, `Gene`, and `* MaxLFQ Intensity` columns for FragPipe.
