# Protein QC Seam Map

## Goal

Apply exact-source extraction checkpoints for `R/func_prot_qc.R` while keeping
live protein QC and rollup behavior frozen behind the current public APIs.

## Current Position In The Flow

- target: `R/func_prot_qc.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `archival closeout verified on 2026-04-13; R/func_prot_qc.R remains a reconciled breadcrumb stub and the focused gate stays green`
- next step: `Treat this handover as archival for R/func_prot_qc.R and continue bucket 4 under a fresh classification/handover pass for the remaining QC and rollup files.`

## Existing Safety Net

- `tests/testthat/test-prot-02-qc-filtering.R`
- `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
- `tests/testthat/test-prot-03-rollup.R`

## Notes

- This bucket previously had no active handover; this file now records the
  first protein-level QC stabilization stop point.
- Wave 1 manifest:
  `tools/refactor/manifest-qc-protein-wave1.yml`
  now verifies and applies cleanly into
  `R/func_prot_qc_filtering_helpers.R`,
  covering:
  `removeEmptyRows()`,
  `removeProteinsWithOnlyOneReplicateHelper()`,
  `removeRowsWithMissingValues()`,
  and
  `removeRowsWithMissingValuesPercentHelper()`.
- The live collate artifact was emitted at
  `tools/refactor/collate-qc-protein-wave1.txt`.
- `DESCRIPTION` `Collate:` now includes
  `func_prot_qc_filtering_helpers.R`
  immediately after
  `func_prot_qc.R`.
- The focused gate reran green after the live apply checkpoint:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
- Wave 2 manifest:
  `tools/refactor/manifest-qc-protein-wave2.yml`
  now verifies, stages, and applies cleanly into staged
  `tools/refactor/staging/wave2_proteomics_qc_protein_correlation_helpers/R/func_prot_qc_correlation_helpers.R`,
  and live
  `R/func_prot_qc_correlation_helpers.R`,
  covering:
  `getPairsOfSamplesTable()`,
  `calulatePearsonCorrelation()`,
  `calculatePearsonCorrelationMatrix()`,
  `calculatePearsonCorrelationOptimized()`,
  `calulatePearsonCorrelationForSamplePairsHelper()`,
  `filterSamplesByPeptideCorrelationThreshold()`,
  `findSamplesPairBelowPeptideCorrelationThreshold()`,
  and
  `filterSamplesByProteinCorrelationThresholdHelper()`.
- The staged wave-2 collate artifact now exists at
  `tools/refactor/staging/wave2_proteomics_qc_protein_correlation_helpers/tools/refactor/collate-qc-protein-wave2.txt`.
- The live wave-2 collate artifact now exists at
  `tools/refactor/collate-qc-protein-wave2.txt`.
- `DESCRIPTION` `Collate:` now also includes
  `func_prot_qc_correlation_helpers.R`
  immediately after
  `func_prot_qc_filtering_helpers.R`.
- The focused gate reran green after the live wave-2 apply checkpoint:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
- `tests/testthat/test-prot-02-qc-filtering.R` and
  `tests/testthat/test-prot-03-rollup.R`
  still skip snapshot-validity assertions when the cp02/cp03 fixtures are only
  Git LFS pointers and the binary artifacts are absent.
- Wave 3 manifest:
  `tools/refactor/manifest-qc-protein-wave3.yml`
  now verifies, stages, and applies cleanly into staged
  `tools/refactor/staging/wave3_proteomics_qc_protein_support_helpers/R/func_prot_qc_support_helpers.R`,
  and live
  `R/func_prot_qc_support_helpers.R`,
  covering:
  `avgReplicateProteinIntensity()`,
  `calculatePercentMissingPeptidePerReplicate()`,
  `calculatePercentMissingProteinPerReplicate()`,
  `calculatePercentMissingPerProtein()`,
  `calculateMissingValuesPerProteinFishersTest()`,
  `getRowsToKeepList()`,
  `averageValuesFromReplicates()`,
  and
  `proteinTechRepCorrelationHelper()`.
- The staged wave-3 collate artifact now exists at
  `tools/refactor/staging/wave3_proteomics_qc_protein_support_helpers/tools/refactor/collate-qc-protein-wave3.txt`.
- The live wave-3 collate artifact now exists at
  `tools/refactor/collate-qc-protein-wave3.txt`.
- `DESCRIPTION` `Collate:` now also includes
  `func_prot_qc_support_helpers.R`
  immediately after
  `func_prot_qc_correlation_helpers.R`.
- The focused gate reran green after the live wave-3 apply checkpoint:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
- After the live wave-3 apply, `R/func_prot_qc.R` is down to `1255` lines with
  `8` remaining top-level functions and still classifies as
  `direct-extraction-ready`.
- Wave 4 manifest:
  `tools/refactor/manifest-qc-protein-wave4.yml`
  now verifies and stages cleanly into
  `tools/refactor/staging/wave4_proteomics_qc_protein_reporting_shell/R/func_prot_qc_reporting_helpers.R`
  and
  `tools/refactor/staging/wave4_proteomics_qc_protein_reporting_shell/R/func_prot_qc_replicate_helpers.R`,
  covering:
  `checkProteinNAPercentages()`,
  `getProteinNARecommendations()`,
  `validatePostImputationProteinData()`,
  `getSamplesCorrelationMatrix()`,
  `updateProteinFiltering()`,
  and
  `removeProteinWithOnlyOneReplicate()`.
- The staged wave-4 collate artifact now exists at
  `tools/refactor/staging/wave4_proteomics_qc_protein_reporting_shell/tools/refactor/collate-qc-protein-wave4.txt`.
- The focused gate reran green after the staged wave-4 checkpoint:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
- Wave 4 now also applies live via
  `tools/refactor/manifest-qc-protein-wave4.yml`
  into
  `R/func_prot_qc_reporting_helpers.R`
  and
  `R/func_prot_qc_replicate_helpers.R`,
  covering the remaining protein-QC reporting and replicate shell:
  `checkProteinNAPercentages()`,
  `getProteinNARecommendations()`,
  `validatePostImputationProteinData()`,
  `getSamplesCorrelationMatrix()`,
  `updateProteinFiltering()`,
  and
  `removeProteinWithOnlyOneReplicate()`.
- The live wave-4 collate artifact now exists at
  `tools/refactor/collate-qc-protein-wave4.txt`.
- `DESCRIPTION` `Collate:` now also includes
  `func_prot_qc_reporting_helpers.R`
  and
  `func_prot_qc_replicate_helpers.R`
  immediately after
  `func_prot_qc_support_helpers.R`.
- The focused gate reran green after the live wave-4 apply checkpoint:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
- `tests/testthat/test-prot-02-qc-filtering.R` and
  `tests/testthat/test-prot-03-rollup.R`
  still skip snapshot-validity assertions when the cp02/cp03 fixtures are only
  Git LFS pointers and the binary artifacts are absent.
- Live `R/func_prot_qc.R` is now a `34`-line breadcrumb stub with `0`
  remaining top-level functions for this target file, and its stale extraction
  TODO inventory has been removed so the file only documents the helper-file
  ownership that remains after the applied waves.
- `R/func_prot_qc.R` no longer needs additional stabilization work; bucket 4
  remains open for the remaining QC and rollup files.
- On 2026-04-13, the focused gate reran green again for:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
  while the cp02/cp03 snapshot-validity assertions continued to skip when the
  Git LFS binary fixtures were absent and only pointer files were present.
