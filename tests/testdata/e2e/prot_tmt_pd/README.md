# prot_tmt_pd E2E Fixture

Synthetic Proteome Discoverer TMT export for the `prot_tmt` lane (KIT-101 revision 2.2, F-07). Written by
`tests/testdata/e2e/generate_fixtures.R` from the lane's `fixture_recipe` in `manifest.json`, which also records each
file's sha256. The earlier six-sample fixture stays in `prot_tmt/`, where the native-session tests and the KIT-021
lifecycle members read it.

## Data Shape

| Property | Value |
|---|---|
| File | `pd_tmt_peptides.tsv`: one row per protein (`Annotated.Sequence`, `Master.Protein.Accessions`), wide abundances |
| Proteins | 400 (P00001–P00400) |
| Samples | 10 in one TMT10 plex: WT_1–WT_5 on 126–128C, KO_1–KO_5 on 129N–131 |
| Abundance columns | `Abundance: F1: <channel>, <sample>`; the importer names each run `<channel>_<sample>` |
| Missing | about 2% of cells left empty, as Proteome Discoverer writes them |
| FASTA | `proteins.fasta`: the 400-protein universe all five proteomics lanes share |
| Design | `design.tsv` (`sample`, `group`); the design builder matches `WT_1` to the run `126_WT_1` (or its sanitised form) |

## Signal

Protein abundances come from the KIT-021 protein corpus recipe (`tools/profiling/kit_protein_corpus.R`,
`kitProteinRecipe("tmt")`, `kitProteinMatrix`, seed 20261110):
- log2 abundance about N(20, 2.2), with residual SD 0.4;
- P00001–P00040 (10%) carry a KO effect drawn from N(0, 1.2).
