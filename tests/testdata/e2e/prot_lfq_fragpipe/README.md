# prot_lfq_fragpipe E2E Fixture

Synthetic FragPipe `combined_protein.tsv` for the `prot_lfq_fragpipe` lane (KIT-101 revision 2.2, F-07). Written by
`tests/testdata/e2e/generate_fixtures.R` from the lane's `fixture_recipe` in `manifest.json`, which also records each
file's sha256. The earlier six-sample table, `prot_lfq/seed_combined_protein.tsv`, is kept for the KIT-021 lifecycle
members.

## Data Shape

| Property | Value |
|---|---|
| File | `combined_protein.tsv`: `Protein`, `Protein ID`, `Entry Name`, `Gene`, `Length`, `Organism`, `Protein Description`, `Protein Existence`, `Coverage`, `Protein Probability`, `<sample> MaxLFQ Intensity` |
| Proteins | 400 (P00001–P00400) |
| Samples | 10 (WT_1–WT_5, KO_1–KO_5) |
| Missing | about 8% of cells, written as 0 (read as NA since F-18; see `prot_lfq_maxquant/README.md`) |
| FASTA | `proteins.fasta`: the 400-protein universe all five proteomics lanes share |
| Design | `design.tsv` (`sample`, `group`) |

## Signal

Protein abundances come from the KIT-021 protein corpus recipe (`kitProteinRecipe("lfq")`, `kitProteinMatrix`, seed
20261112). P00001–P00040 (10%) carry a KO effect drawn from N(0, 1.2). As in `prot_lfq_maxquant/`, about 8% of
cells are written as 0.
