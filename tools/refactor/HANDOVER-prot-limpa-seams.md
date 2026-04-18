# Protein Limpa Seam Map

## Goal

Reconcile `R/func_prot_limpa.R` as a breadcrumb stub while keeping the
live proteomics limpa behavior frozen behind the current public APIs.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_prot_limpa.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `breadcrumb_stub_reconciled_after_limpa_qc_wave3_apply`
- next step: `No further seams remain in R/func_prot_limpa.R; if later limpa cleanup is needed, classify the extracted live owner files separately rather than reopening the breadcrumb stub.`

## Existing Safety Net

- `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`

## Notes

Manual target bucket 0 for protein limpa stabilization.

- This manual target previously had no active handover; this file now records
  the first protein-limpa stabilization stop point.
- April 16, 2026 checkpoint drafted, staged, reviewed, and applied
  `tools/refactor/manifest-prot-limpa-wave1.yml`, moving the exact
  `proteinMissingValueImputationLimpa` `PeptideQuantitativeData` DPC-Quant
  method block out of `R/func_prot_limpa.R` into the live
  `R/func_prot_s4_missingness.R` target, preserving the existing missingness
  helper symbols there, and passing
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-prot-limpa-wave1.yml`.
- `tools/refactor/collate-prot-limpa-wave1.txt` records that the existing
  `DESCRIPTION` `Collate:` order already loaded `func_prot_s4_missingness.R`
  before `func_prot_limpa.R`, so this checkpoint did not require a load-order
  change.
- Live `R/func_prot_limpa.R` is now a `941`-line helper-only breadcrumb with
  `2` top-level expressions, keeping only `generateLimpaQCPlots()` and
  `convertDpcDAToStandardFormat()`.
- The stale duplicate `proteinMissingValueImputationLimpa`
  `ProteinQuantitativeData` method block was removed from
  `R/func_prot_limpa.R`, so both proteomics limpa S4 methods now live in
  `R/func_prot_s4_missingness.R`, which measures `650` lines after the applied
  wave.
- The focused design gate reran green with `1572` passes and the same expected
  Git LFS snapshot skip:
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
- April 16, 2026 checkpoint drafted, staged, reviewed, and applied
  `tools/refactor/manifest-prot-limpa-wave2.yml`, moving the exact
  `convertDpcDAToStandardFormat()` helper block out of
  `R/func_prot_limpa.R` into the live
  `R/func_prot_limpa_da_helpers.R` target, passing both
  `Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-prot-limpa-wave2.yml`
  and
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-prot-limpa-wave2.yml`.
- `tools/refactor/collate-prot-limpa-wave2.txt` recorded the new helper target,
  and `DESCRIPTION` `Collate:` now loads
  `func_prot_limpa_da_helpers.R` immediately before `func_prot_limpa.R`.
- Live `R/func_prot_limpa.R` now measures `812` lines and is a QC-only
  breadcrumb with `1` top-level expression, while the new
  `R/func_prot_limpa_da_helpers.R` leaf measures `129` lines.
- The focused design gate reran green again with `1572` passes and the same
  expected Git LFS snapshot skip:
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
- `R/func_prot_limpa.R` still needs more stabilization work because
  `generateLimpaQCPlots()` remains the last standalone helper surface; the next
  checkpoint should continue with a bounded exact-source helper wave rather than
  more live S4 seam cleanup.
- April 16, 2026 checkpoint drafted and staged
  `tools/refactor/manifest-prot-limpa-wave3.yml`, resolving the exact
  `generateLimpaQCPlots()` symbol from `R/func_prot_limpa.R:33-811` into staged
  target
  `tools/refactor/staging/wave3_proteomics_prot_limpa_qc_helper/R/func_prot_limpa_qc_helpers.R`
  without rewriting live `R/` sources, and passing
  `Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-prot-limpa-wave3.yml`.
- `tools/refactor/collate-prot-limpa-wave3.txt` now records the future helper
  target `R/func_prot_limpa_qc_helpers.R`, which should load before
  `func_prot_limpa.R` when the staged wave is reviewed and applied.
- The staged QC helper leaf measures `780` lines while live
  `R/func_prot_limpa.R` still measures `812` lines until the wave is applied,
  so the public wrapper identity remains in place for this stop point.
- The focused design gate reran green again with `1572` passes and the same
  expected Git LFS snapshot skip:
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
- April 16, 2026 checkpoint reviewed and applied
  `tools/refactor/manifest-prot-limpa-wave3.yml`, moving the exact
  `generateLimpaQCPlots()` helper body out of `R/func_prot_limpa.R` into the
  live `R/func_prot_limpa_qc_helpers.R` target, passing both
  `Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-prot-limpa-wave3.yml`
  and
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-prot-limpa-wave3.yml`.
- `tools/refactor/collate-prot-limpa-wave3.txt` now reflects the live QC helper
  target, and `DESCRIPTION` `Collate:` now loads
  `func_prot_limpa_qc_helpers.R` immediately before `func_prot_limpa.R`.
- Live `R/func_prot_limpa.R` now measures `33` lines with `0` top-level
  expressions, while the new `R/func_prot_limpa_qc_helpers.R` leaf measures
  `780` lines and owns `generateLimpaQCPlots()`.
- The focused design gate reran green again with `1572` passes and the same
  expected Git LFS snapshot skip:
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
- April 16, 2026 checkpoint reconciled the stale breadcrumb shell in
  `R/func_prot_limpa.R`, updating the header to point only at the live helper
  owners in `R/func_prot_s4_missingness.R`,
  `R/func_prot_limpa_da_helpers.R`, and
  `R/func_prot_limpa_qc_helpers.R`, and adding an explicit intentionally-empty
  breadcrumb stub marker with no executable code.
- The focused design gate reran green again with `1572` passes and the same
  expected Git LFS snapshot skip:
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
- Live `R/func_prot_limpa.R` is now a `32`-line breadcrumb stub with `0`
  top-level expressions and no stale ownership claims.
- `R/func_prot_limpa.R` no longer needs additional stabilization work.
  Any later limpa cleanup should move to the extracted live owner files under
  a fresh classification and handover pass.
