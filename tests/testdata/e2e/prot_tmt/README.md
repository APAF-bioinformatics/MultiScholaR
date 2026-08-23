# prot_tmt E2E Fixture

## Data Shape

| Property | Value |
|---|---|
| File | pd_tmt_peptides.tsv |
| Format | Proteome Discoverer TMT peptide export (wide format) |
| Proteins | 5 (P00001–P00005) |
| Samples | 6 (WT_1–3, KO_1–3) |
| TMT Channels | 126=WT_1, 127N=WT_2, 127C=WT_3, 128N=KO_1, 128C=KO_2, 129N=KO_3 |
| FASTA | proteins.fasta (P00001–P00006 UniProt-style headers) |

## Groups

| Run ID (after import) | Sample | Group | TMT Channel |
|---|---|---|---|
| 126_WT_1 | WT_1 | WT | 126 |
| 127N_WT_2 | WT_2 | WT | 127N |
| 127C_WT_3 | WT_3 | WT | 127C |
| 128N_KO_1 | KO_1 | KO | 128N |
| 128C_KO_2 | KO_2 | KO | 128C |
| 129N_KO_3 | KO_3 | KO | 129N |

## Column Format

Abundance columns follow Proteome Discoverer export convention:
`Abundance: F{fraction}: {channel}, {sample}`

The `importProteomeDiscovererTMTData()` function reads peptide rows keyed by
`Annotated.Sequence` and `Master.Protein.Accessions`, then pivots abundance
channels to long-format run IDs.

## Raw Fixture Effects (KO_vs_WT)

| Protein | Gene | Direction | Approx log2FC |
|---|---|---|---|
| P00001 | GENE1 | Up in KO | ~1.03 |
| P00002 | GENE2 | Down in KO | ~-0.90 |

P00003-P00005 have raw absolute log2 fold changes below 0.1. This five-protein
fixture validates parser identity, workflow routing, and current transform
parity. It does not establish differential-abundance sensitivity or power.
The current default cyclic-loess workflow produces no BH-significant proteins
from this intentionally tiny fixture; that result is frozen separately as a
current-behavior oracle.

## Import Function

`importProteomeDiscovererTMTData()` in `R/func_prot_import_tmt.R`

Required columns: `Annotated.Sequence`, `Master.Protein.Accessions`, `Abundance.*` columns
