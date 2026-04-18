# Peptide QC Seam Map

## Goal

Stage exact-source extraction for func_prot_qc_peptide.R while keeping live peptide-QC behavior frozen behind the current public APIs.

## Current Position In The Flow

- target: `R/func_prot_qc_peptide.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `wave 4 from manifest-qc-peptide-wave4.yml is now applied live into R/func_prot_qc_peptide_methods.R`
- next step: `Manual target is complete; keep this handover as the archival seam record and continue bucket 4 on the remaining QC and rollup files when needed.`

## Existing Safety Net

- `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
- `tests/testthat/test-prot-02-qc-filtering.R`
- `tests/testthat/test-prot-03-rollup.R`

## Notes

Manual target bucket 0 for proteomics peptide QC stabilization.

- Wave 1 manifest:
  `tools/refactor/manifest-qc-peptide-wave1.yml`
  now verifies and applies cleanly against the current source tree, moving:
  `check_case_collision_columns()`,
  `peptideIntensityFilteringHelper()`,
  and
  `removePeptidesWithMissingValuesPercentHelper()`
  into
  `R/func_prot_qc_peptide_group_filters.R`.
- `DESCRIPTION` `Collate:` now includes
  `func_prot_qc_peptide_group_filters.R`
  immediately after
  `func_prot_qc_peptide.R`
  so the shared
  `resolvePeptideQcColumnName()`
  seam remains available before the extracted helpers load.
- Focused gate rerun now passes:
  - `tests/testthat/test-prot-02-qc-filtering.R`
    skips the snapshot-validity check when `cp02_qc_filtered_peptide.rds` is only a Git LFS pointer and still exercises the live helper coverage
  - `tests/testthat/test-prot-03-rollup.R`
    skips the snapshot-validity check when `cp03_rolled_up_protein.rds` is only a Git LFS pointer and still exercises the live rollup helper coverage
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
    still passes against the live code path after the earlier column-resolution seam
- `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  already passes after routing both staged wave 1 helpers through the new
  `resolvePeptideQcColumnName()` seam, which preserves bare-symbol defaults
  without forcing missing objects during column-name resolution.
- Wave 2 manifest:
  `tools/refactor/manifest-qc-peptide-wave2.yml`
  now verifies, stages cleanly into
  `tools/refactor/staging/wave2_proteomics_qc_peptide_replicate_filters/R/func_prot_qc_peptide_replicate_filters.R`,
  and applies live into
  `R/func_prot_qc_peptide_replicate_filters.R`,
  covering:
  `removePeptidesWithOnlyOneReplicateHelper()`,
  `filterMinNumPeptidesPerProteinHelper()`,
  `filterMinNumPeptidesPerSampleHelper()`,
  and
  `srlQvalueProteotypicPeptideCleanHelper()`.
- `DESCRIPTION` `Collate:` now includes
  `func_prot_qc_peptide_replicate_filters.R`
  immediately after
  `func_prot_qc_peptide_group_filters.R`
  and the apply-time collate artifact was regenerated at
  `tools/refactor/collate-qc-peptide-wave2.txt`.
- The focused gate reruns green after the live apply checkpoint:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
- After the live apply, `R/func_prot_qc_peptide.R` is down to `897` lines with
  `6` remaining top-level functions and still classifies as
  `direct-extraction-ready`.
- Wave 3 manifest:
  `tools/refactor/manifest-qc-peptide-wave3.yml`
  now verifies and stages cleanly into
  `tools/refactor/staging/wave3_proteomics_qc_peptide_support_helpers/R/func_prot_qc_peptide_support.R`,
  covering:
  `checkPeptideNAPercentages()`,
  `removePeptidesOnlyInHek293()`,
  `compareTwoPeptideDataObjects()`,
  `summarisePeptideObject()`,
  and
  `calculatePeptidePearsonCorrelation()`.
- The staged wave-3 collate artifact was emitted at
  `tools/refactor/staging/wave3_proteomics_qc_peptide_support_helpers/tools/refactor/collate-qc-peptide-wave3.txt`.
- Wave 3 now also applies live into
  `R/func_prot_qc_peptide_support.R`,
  covering the remaining non-method helper cluster without hand-rewriting the
  extracted bodies.
- `DESCRIPTION` `Collate:` now also includes
  `func_prot_qc_peptide_support.R`
  immediately after
  `func_prot_qc_peptide_replicate_filters.R`
  and the apply-time collate artifact was regenerated at
  `tools/refactor/collate-qc-peptide-wave3.txt`.
- After the live apply, `R/func_prot_qc_peptide.R` is down to `579` lines and
  retains only the shared `resolvePeptideQcColumnName()` seam plus the
  remaining S4 wrapper methods.
- The focused gate reruns green after the live apply checkpoint:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
- The next safe stop point is staging one final exact-source wave for the
  remaining S4 wrapper methods before deciding whether the manual target can be
  archived.
- Wave 4 manifest:
  `tools/refactor/manifest-qc-peptide-wave4.yml`
  now verifies and stages cleanly into
  `tools/refactor/staging/wave4_proteomics_qc_peptide_s4_methods/R/func_prot_qc_peptide_methods.R`,
  covering the remaining S4 wrapper methods:
  `peptideIntensityFiltering()`,
  `removePeptidesWithMissingValuesPercent()`,
  `removePeptidesWithOnlyOneReplicate()`,
  `filterMinNumPeptidesPerProtein()`,
  `filterMinNumPeptidesPerSample()`,
  and
  `srlQvalueProteotypicPeptideClean()`.
- The staged wave-4 collate artifact was emitted at
  `tools/refactor/staging/wave4_proteomics_qc_peptide_s4_methods/tools/refactor/collate-qc-peptide-wave4.txt`.
- The focused gate reran green after the staging checkpoint:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
- Wave 4 now also applies live into
  `R/func_prot_qc_peptide_methods.R`,
  covering the remaining S4 wrapper-method cluster without hand-rewriting the
  extracted bodies.
- `DESCRIPTION` `Collate:` now also includes
  `func_prot_qc_peptide_methods.R`
  immediately after
  `func_prot_qc_peptide_support.R`
  and the apply-time collate artifact was regenerated at
  `tools/refactor/collate-qc-peptide-wave4.txt`.
- The focused gate reran green after the live apply checkpoint:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-03-rollup.R`
- `R/func_prot_qc_peptide.R` is now down to `245` lines and retains only the
  shared `resolvePeptideQcColumnName()` seam.
- Manual bucket 0 target is complete; this handover is now archival.
