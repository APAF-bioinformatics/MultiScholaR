# Peptide S4 Seam Map

## Goal

Freeze peptide constructor and matrix-building behavior, then stage exact extraction waves for the remaining peptide S4 methods without changing the public S4 surface.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_pept_s4_objects.R`
- classification: `direct-extraction-ready`
- next step: `Stage the adjacent peptide plotPca QC block into func_pept_s4_qc_methods.R now that the first peptide QC helper is live and ordered before the wrapper.`
- latest completed checkpoint: `peptide_s4_plot_rle_wave_applied`
- latest completed checkpoint summary:
  `completed one bounded stabilization checkpoint by adding a direct peptide plotRle characterization to tests/testthat/test-prot-04-design.R, drafting/staging/reviewing/applying tools/refactor/manifest-pept-s4-wave3.yml, materializing R/func_pept_s4_qc_methods.R with the exact plotRle PeptideQuantitativeData method block, removing that block from R/func_pept_s4_objects.R, updating DESCRIPTION Collate so func_pept_s4_qc_methods.R now loads between func_pept_s4_norm_methods.R and func_pept_s4_objects.R, and rerunning the focused design gate green with 1511 passes and the same expected Git LFS snapshot skip.`

## Existing Safety Net

- `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`

## Notes

Manual target bucket 0 for `R/func_pept_s4_objects.R`.
After the wave 3 apply checkpoint, `R/func_pept_s4_objects.R` measures `1516` lines while `R/func_pept_s4_core.R` remains `205` lines, `R/func_pept_s4_norm_methods.R` remains `153` lines, and the new `R/func_pept_s4_qc_methods.R` helper carries `36` exact-source lines.
`tools/refactor/collate-pept-s4-wave1.txt` records the new helper target and the new `DESCRIPTION` collate entry now loads `func_pept_s4_core.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave2.txt` records the normalization helper target and the updated `DESCRIPTION` collate entry now loads `func_pept_s4_norm_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave3.txt` records the first peptide QC helper target and the updated `DESCRIPTION` collate entry now loads `func_pept_s4_qc_methods.R` before `func_pept_s4_objects.R`.
The next bounded checkpoint should stay exact-source and move one adjacent peptide QC or imputation block into the live peptide QC helper without reopening the class-definition or normalization surfaces.
