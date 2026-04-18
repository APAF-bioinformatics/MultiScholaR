# Peptide S4 Seam Map

## Goal

Freeze peptide constructor and matrix-building behavior, then stage exact extraction waves for the remaining peptide S4 methods without changing the public S4 surface.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_pept_s4_objects.R`
- classification: `direct-extraction-ready`
- next step: `No additional live method blocks remain in func_pept_s4_objects.R; treat the file as a stabilized breadcrumb/generics shell and shift later work to a different backlog target if needed.`
- latest completed checkpoint: `peptide_s4_plot_density_ggplot_wave_applied`
- latest completed checkpoint summary:
  `completed one bounded stabilization checkpoint by reviewing and applying tools/refactor/manifest-pept-s4-wave14.yml, moving the exact ggplot2::ggplot plotDensity method block into the live R/func_pept_s4_qc_methods.R target, removing that block from R/func_pept_s4_objects.R, preserving the existing DESCRIPTION collate ordering for the already-live QC helper, passing tools/refactor/check_wave_apply.R, and rerunning the focused design gate green with 1572 passes and the same expected Git LFS snapshot skip.`

## Existing Safety Net

- `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`

## Notes

Manual target bucket 0 for `R/func_pept_s4_objects.R`.
The first fourteen reviewed peptide S4 waves are now live, with `R/func_pept_s4_core.R` preserving the public `PeptideQuantitativeData` class shell plus `PeptideQuantitativeDataDiann` and `calcPeptideMatrix`, `R/func_pept_s4_norm_methods.R` holding the reviewed normalization, negative-control, and RUV method cluster, `R/func_pept_s4_qc_methods.R` carrying the reviewed `plotRle`, `plotPcaDispatch`, both `plotPca` methods, both `plotDensity` methods, `plotPearson`, and `pearsonCorForSamplePairs` blocks, `R/func_pept_s4_accession_methods.R` holding the reviewed `chooseBestProteinAccession` method, and `R/func_pept_s4_missingness.R` carrying the reviewed `peptideMissingValueImputationLimpa` method.
After the wave 14 apply checkpoint, `R/func_pept_s4_objects.R` measures `162` lines while `R/func_pept_s4_core.R` remains `205` lines, `R/func_pept_s4_norm_methods.R` remains `851` lines, `R/func_pept_s4_qc_methods.R` now measures `387` lines, `R/func_pept_s4_accession_methods.R` remains `110` exact-source lines, and `R/func_pept_s4_missingness.R` remains `210` lines.
`tools/refactor/collate-pept-s4-wave1.txt` records the new helper target and the new `DESCRIPTION` collate entry now loads `func_pept_s4_core.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave2.txt` records the normalization helper target and the updated `DESCRIPTION` collate entry now loads `func_pept_s4_norm_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave3.txt` records the first peptide QC helper target and the updated `DESCRIPTION` collate entry now loads `func_pept_s4_qc_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave4.txt` records the unchanged peptide QC helper target while confirming the existing `DESCRIPTION` collate entry still loads `func_pept_s4_qc_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave5.txt` records the unchanged peptide QC helper target while confirming the existing `DESCRIPTION` collate entry still loads `func_pept_s4_qc_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave6.txt` records the unchanged peptide QC helper target while confirming the existing `DESCRIPTION` collate entry still loads `func_pept_s4_qc_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave7.txt` records the new accession helper target and the updated `DESCRIPTION` collate entry now loads `func_pept_s4_accession_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave8.txt` records the unchanged normalization helper target while confirming the existing `DESCRIPTION` collate entry still loads `func_pept_s4_norm_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave9.txt` records the unchanged normalization helper target while confirming the existing `DESCRIPTION` collate entry still loads `func_pept_s4_norm_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave10.txt` records the unchanged normalization helper target while confirming the existing `DESCRIPTION` collate entry still loads `func_pept_s4_norm_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave11.txt` records the unchanged normalization helper target while confirming the existing `DESCRIPTION` collate entry still loads `func_pept_s4_norm_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave12.txt` records the unchanged normalization helper target while confirming the existing `DESCRIPTION` collate entry still loads `func_pept_s4_norm_methods.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave13.txt` records the new missingness helper target and the updated `DESCRIPTION` collate entry now loads `func_pept_s4_missingness.R` before `func_pept_s4_objects.R`.
`tools/refactor/collate-pept-s4-wave14.txt` records the unchanged peptide QC helper target while confirming the existing `DESCRIPTION` collate entry still loads `func_pept_s4_qc_methods.R` before `func_pept_s4_objects.R`.
April 16, 2026 wave 10 applied `wave10_proteomics_peptide_s4_ruv_cancor_fast`, moving the exact `ruvCancorFast` `PeptideQuantitativeData` method block into `R/func_pept_s4_norm_methods.R`, removing it from `R/func_pept_s4_objects.R`, passing the manifest verify/apply checks, and rerunning the focused design gate green with `1572` passes plus the same expected Git LFS snapshot skip.
April 16, 2026 wave 11 applied `wave11_proteomics_peptide_s4_ruv_cancor`, moving the exact `ruvCancor` `PeptideQuantitativeData` method block into `R/func_pept_s4_norm_methods.R`, removing it from `R/func_pept_s4_objects.R`, passing the manifest verify/apply checks, and rerunning the focused design gate green with `1572` passes plus the same expected Git LFS snapshot skip.
April 16, 2026 wave 12 applied `wave12_proteomics_peptide_s4_ruviii_c_varying`, moving the exact `ruvIII_C_Varying` `PeptideQuantitativeData` method block into `R/func_pept_s4_norm_methods.R`, removing it from `R/func_pept_s4_objects.R`, passing the manifest verify/apply checks, and rerunning the focused design gate green with `1572` passes plus the same expected Git LFS snapshot skip.
April 16, 2026 wave 13 applied `wave13_proteomics_peptide_s4_missingness_limpa`, moving the exact `peptideMissingValueImputationLimpa` `PeptideQuantitativeData` method block into `R/func_pept_s4_missingness.R`, removing that block from `R/func_pept_s4_objects.R`, recording the new missingness-helper collate target in `tools/refactor/collate-pept-s4-wave13.txt`, updating `DESCRIPTION` so `func_pept_s4_missingness.R` loads before `func_pept_s4_objects.R`, and rerunning the focused design gate green with `1572` passes plus the same expected Git LFS snapshot skip.
April 16, 2026 wave 14 staged `wave14_proteomics_peptide_s4_plot_density_ggplot`, drafting
`tools/refactor/manifest-pept-s4-wave14.yml`, verifying it against current
sources, and staging the exact `plotDensity` `ggplot2::ggplot` method block
into
`tools/refactor/staging/wave14_proteomics_peptide_s4_plot_density_ggplot/R/func_pept_s4_qc_methods.R`
with `tools/refactor/collate-pept-s4-wave14.txt` recording the unchanged QC
helper collate target, while keeping live `R/func_pept_s4_objects.R`,
`R/func_pept_s4_qc_methods.R`, and `DESCRIPTION` unchanged; the focused design
gate stayed green with `1572` passes plus the same expected Git LFS snapshot
skip.
April 16, 2026 wave 14 applied `wave14_proteomics_peptide_s4_plot_density_ggplot`, moving the exact `plotDensity` `ggplot2::ggplot` method block into `R/func_pept_s4_qc_methods.R`, removing that block from `R/func_pept_s4_objects.R`, preserving the existing `DESCRIPTION` collate ordering for the already-live QC helper target, passing the manifest verify/apply checks, and rerunning the focused design gate green with `1572` passes plus the same expected Git LFS snapshot skip.
No additional live method blocks remain in `R/func_pept_s4_objects.R`; treat the file as a stabilized breadcrumb/generics shell for this target and move later work to a different backlog target if further bucket work is needed.
