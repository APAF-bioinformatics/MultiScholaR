# E2E-005 Proteomics TMT and LFQ Browser Matrix

E2E-005 covers the protein-level proteomics branches that intentionally diverge
from the DIA peptide-to-protein path. The matrix proves that TMT and LFQ can be
loaded through the GUI, assigned to a design, advanced through the protein-level
bypass semantics, normalized, exported as a filtered session, reloaded in DA,
analyzed, rendered through the correct report template, and exported as a summary
session without falling back to the DIA report route.

## Covered Lanes

| Lane | Import entry point | Fixture | Workflow | Report |
|---|---|---|---|---|
| `prot_tmt` | Proteome Discoverer TMT | `prot_tmt/pd_tmt_peptides.tsv` | `TMT` | `TMT_report.rmd` |
| `prot_lfq` | MaxQuant LFQ | `prot_lfq/proteinGroups.txt` | `LFQ` | `LFQ_report.rmd` |
| `prot_lfq_fragpipe` | FragPipe LFQ | `prot_lfq/seed_combined_protein.tsv` | `LFQ` | `LFQ_report.rmd` |

## Runtime Assertions

The browser test asserts each branch at the same high-value handoff points:

- setup/import completes and sets `design_matrix` to `pending`;
- the state digest reports the expected workflow type (`TMT` or `LFQ`);
- design assignment maps sanitized imported run names back to manifest samples;
- protein-level workflows mark `quality_control` complete and expose
  normalization without running the DIA peptide QC sequence;
- normalization exports `filtered_session_data_latest.rds` and
  `filtered_session_summary.txt`;
- the filtered session records workflow type, R6 current state, complete state
  map, state history, protein count, and sample count;
- DA reload restores the filtered session and completes differential abundance;
- session summary writes `study_parameters.txt`, renders the branch-specific
  report, and exports `session_state_*.RDS`;
- non-DIA branches do not create or route through `DIANN_report.rmd`;
- browser console errors and Shiny error notifications fail the test.

## Fixture Notes

The TMT fixture uses peptide-export naming from Proteome Discoverer:
`Abundance: F{n}: {channel}, {sample}` and
`Master.Protein.Accessions`. The importer also tolerates the simplified
`Abundance.{channel}` form so direct importer tests cover both real and reduced
exports.

The LFQ fixture directory contains both supported LFQ entry surfaces. MaxQuant is
represented by `LFQ.intensity.{sample}` columns in `proteinGroups.txt`; FragPipe
is represented by `{sample} MaxLFQ Intensity` columns in
`seed_combined_protein.tsv`. Both resolve to a long protein table with
`Protein.Ids`, `Run`, and `Intensity` so downstream design, normalization, DA,
and report modules receive the same contract.

## Test Entrypoint

Run the focused ticket suite with:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(e2e-manifest|e2e-fixture-integrity|e2e-browser-harness|prot-import-tmt-direct-shared|prot-import-maxquant-direct-shared|prot-01c-import-module-contracts|e2e-proteomics-tmt-lfq)$",
  reporter = "summary"
)
```
