# prot_lfq_maxquant E2E Fixture

Synthetic MaxQuant `proteinGroups.txt` for the `prot_lfq` lane (KIT-101 revision 2.2, F-07). Written by
`tests/testdata/e2e/generate_fixtures.R` from the lane's `fixture_recipe` in `manifest.json`, which also records each
file's sha256. The earlier six-sample fixture stays in `prot_lfq/`, where the native-session tests and the KIT-021
lifecycle members read it.

## Data Shape

| Property | Value |
|---|---|
| File | `proteinGroups.txt`: dotted MaxQuant headers (`Protein.IDs`, `Majority.protein.IDs`, `Protein.names`, `Gene.names`, `LFQ.intensity.<sample>`, `Potential.contaminant`, `Reverse`) |
| Proteins | 400 (P00001–P00400), and 5 contaminant (`CON__`) and 5 decoy (`REV__`) rows that the reader drops; one of each leads the table, as in the KIT-021 corpus |
| Samples | 10 (WT_1–WT_5, KO_1–KO_5) |
| Missing | about 8% of cells (the KIT-021 recipe), written as 0, as MaxQuant writes "not quantified" |
| FASTA | `proteins.fasta`: the 400-protein universe all five proteomics lanes share |
| Design | `design.tsv` (`sample`, `group`) |

## Signal

Protein abundances come from the KIT-021 protein corpus recipe (`kitProteinRecipe("lfq")`, `kitProteinMatrix`, seed
20261111). P00001–P00040 (10%) carry a KO effect drawn from N(0, 1.2).

Revision 2.2 phase 1 wrote no missing values: the import then kept a 0 as a quantity, log2(0 + 1), and the
correlation filter dropped samples. Since F-18 (8a2a3581) the MaxQuant and FragPipe readers read a 0 as NA. The
recipe's missing values are therefore back, written as 0 (phase 2a).
