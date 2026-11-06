# prot_dia_limpa E2E Fixture

Synthetic DIA-NN main report for the DIA_LIMPA lane (KIT-101 revision 2.2, F-07). It has the same layout as
`prot_dia/report.tsv`, from its own seed. Written by `tests/testdata/e2e/generate_fixtures.R` from the lane's
`fixture_recipe` in `manifest.json`, which also records each file's sha256.

## Data Shape

| Property | Value |
|---|---|
| File | `report.tsv`: DIA-NN main report, long format, one row per precursor and run |
| Proteins | 400 (P00001–P00400, genes GENE1–GENE400); 1 to 4 observed peptides each; charge 2, some also charge 3 |
| Runs | 8 (WT_1–WT_4, KO_1–KO_4) |
| Rows | 8,886 (a missing precursor has no row, as DIA-NN writes it) |
| FASTA | `proteins.fasta`: the 400-protein universe all five proteomics lanes share |
| Design | `design.tsv` (`sample`, `group`) |

The columns and the confidence-filter content are those described in `prot_dia/README.md`.

## Signal and missing values

Protein abundances come from the KIT-021 protein corpus recipe (`kitProteinMatrix`, seed 20261109):
- P00001–P00040 carry a KO effect drawn from N(0, 1.2);
- about 8% of protein-run cells are absent;
- a precursor is more often absent at low intensity (an MNAR component, `plogis(17 - log2 quantity)`).

This gives limpa's detection-probability model missing values to work with.

## Not used by the lane

`seed_report.tsv` is an earlier hand-made table. It is kept, and recorded under `retained_fixtures` in the manifest.
