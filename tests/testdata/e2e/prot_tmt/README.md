# prot_tmt E2E Fixture

## Data Shape

| Property | Value |
|---|---|
| File | seed_proteins.tsv |
| Format | Proteome Discoverer TMT protein export (wide format) |
| Proteins | 5 (P00001–P00005) |
| Samples | 4 (WT_1, WT_2, KO_1, KO_2) |
| TMT Channels | 126=WT_1, 127N=WT_2, 128N=KO_1, 129N=KO_2 |

## Groups

| Run ID (after import) | Sample | Group | TMT Channel |
|---|---|---|---|
| 126_WT_1 | WT_1 | WT | 126 |
| 127N_WT_2 | WT_2 | WT | 127N |
| 128N_KO_1 | KO_1 | KO | 128N |
| 129N_KO_2 | KO_2 | KO | 129N |

## Column Format

Abundance columns follow Proteome Discoverer export convention:
`Abundance: F{fraction}: {channel}, {sample_name}`

The `importProteomeDiscovererTMTData()` function renames these to `{channel}_{sample_name}`
(e.g. `Abundance: F1: 126, WT_1` → `126_WT_1`) then pivots to long format.
The `design.tsv` sample column must use these post-import Run IDs.

## Expected Significant Proteins (KO_vs_WT)

| Protein | Gene | Direction | Approx log2FC |
|---|---|---|---|
| P00001 | GENE1 | UP in KO | ~1.05 |
| P00002 | GENE2 | UP in KO | ~1.01 |

P00003–P00005 are not expected to be significant (<0.1 log2FC between groups).

## Import Function

`importProteomeDiscovererTMTData()` in `R/func_prot_import_tmt.R`

Required columns: `Accession` (renamed to `Protein.Ids`), `Abundance: F[n]: [channel], [sample]` columns
