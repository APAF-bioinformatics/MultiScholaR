# prot_dia E2E Fixture

Synthetic DIA-NN main report for the canonical DIA lane (KIT-101 revision 2.2, F-07). Written by
`tests/testdata/e2e/generate_fixtures.R` from the lane's `fixture_recipe` in `manifest.json`, which also records each
file's sha256.

## Data Shape

| Property | Value |
|---|---|
| File | `report.tsv`: DIA-NN main report, long format, one row per precursor and run |
| Proteins | 400 (P00001–P00400, genes GENE1–GENE400); 1 to 4 observed peptides each; charge 2, some also charge 3 |
| Runs | 8 (WT_1–WT_4, KO_1–KO_4) |
| Rows | 8,989 (a missing precursor has no row, as DIA-NN writes it) |
| FASTA | `proteins.fasta`: the 400-protein universe all five proteomics lanes share; each observed peptide is in its protein's sequence |
| Design | `design.tsv` (`sample`, `group`) |

Columns: `Run`, `Protein.Group`, `Protein.Ids`, `Protein.Names`, `Genes`, `Modified.Sequence` (cysteines carry
`(UniMod:4)`), `Stripped.Sequence`, `Precursor.Id`, `Precursor.Charge`, `Q.Value`, `Global.Q.Value`, `PG.Q.Value`,
`Global.PG.Q.Value`, `Proteotypic`, `Precursor.Quantity`, `Precursor.Normalised`. These are the columns the import
(`importDIANNData`), the Q-value filter (`srlQvalueProteotypicPeptideCleanHelper`, `diann_main_report` mode) and the
precursor rollup (`rollUpPrecursorToPeptideHelper`) read.

The confidence filter removes rows: about 2% of rows fail the run q-value, about 1% of precursors fail the global
q-value, four unchanged proteins fail the global protein-group q-value (0.03), and about 3% of peptides are not
proteotypic.

## Signal

Protein abundances come from the KIT-021 protein corpus recipe (`tools/profiling/kit_protein_corpus.R`,
`kitProteinMatrix`, seed 20261108):
- log2 abundance about N(20, 2.2), with residual SD 0.4;
- P00001–P00040 (10%) carry a KO effect drawn from N(0, 1.2);
- about 8% of protein-run cells are absent.

Each precursor adds a fixed offset and run noise (SD 0.15).

`Precursor.Quantity` is raw and keeps the recipe's per-sample loading. `Precursor.Normalised` removes that loading as
DIA-NN's cross-run normalisation would (KIT-101 F-19): each run is centred on its median log2 quantity, with a small
residual run error (SD 0.1 log2). Before F-19 the normalised column kept the loading, which at seed 20261109 differed
by -0.30 log2 between KO and WT.

## Not used by the lane

`seed_report.tsv` is an earlier hand-made table. It is kept, and recorded under `retained_fixtures` in the manifest.
