# Proteomics S4 Seam Map

## Goal

Freeze constructor and validity behavior, then stage exact extraction waves for proteomics S4 families without changing the public S4 surface.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_prot_s4_objects.R`
- classification: `direct-extraction-ready`
- next step: `No additional live method blocks remain in func_prot_s4_objects.R; treat the file as a stabilized breadcrumb/generics shell and shift later work to a different backlog target if needed.`
- latest completed checkpoint: `filter_min_num_peptides_per_protein_wave_applied`
- latest completed checkpoint summary:
  `completed one bounded stabilization checkpoint by drafting, staging, reviewing, and applying tools/refactor/manifest-prot-s4-wave38.yml, moving the exact filterMinNumPeptidesPerProtein ProteinQuantitativeData method block into the live R/func_prot_s4_qc_methods.R target, removing that block from R/func_prot_s4_objects.R, preserving the existing DESCRIPTION collate ordering for the already-live QC methods helper, passing tools/refactor/check_wave_apply.R, and rerunning the focused design gate green with 1499 passes and the same expected Git LFS snapshot skip.`

## Existing Safety Net

- `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`

## Notes

Bucket 8 target: func_prot_s4_objects.R
The first thirty-eight reviewed proteomics S4 waves are now live, with func_prot_s4_core.R preserving the public ProteinQuantitativeData class shell plus setProteinData, func_prot_s4_missingness.R holding the extracted missing-value/imputation method block plus the reviewed preservePeptideNaValues method and helper, func_prot_s4_da_methods.R carrying both reviewed DA entrypoint methods, func_prot_s4_da_results.R holding the reviewed DA result-export method, func_prot_s4_norm_controls.R holding the reviewed low-CV negative-control selector plus getNegCtrlProtAnova, func_prot_s4_norm_methods.R now holding the reviewed normaliseBetweenSamples, ruvCancor, and ruvIII_C_Varying methods, func_prot_s4_replicates.R holding the reviewed averageTechReps plus getRuvIIIReplicateMatrix methods, func_prot_s4_qc_methods.R now holding the reviewed removeProteinsWithOnlyOneReplicate, removeRowsWithMissingValuesPercent, filterSamplesByProteinCorrelationThreshold, cleanDesignMatrix, proteinIntensityFiltering, plotRle, plotRleList, savePlotRleList, plotPca, plotPcaList, plotPcaBox, plotDensityList, savePlotDensityList, filterMinNumPeptidesPerProtein, getPcaMatrix, proteinTechRepCorrelation, plotPearson, and pearsonCorForSamplePairs blocks, func_prot_s4_accession_methods.R holding the reviewed chooseBestProteinAccession plus chooseBestProteinAccessionSumDuplicates methods, and func_prot_s4_grid.R now holding the reviewed GridPlotData class block, the reviewed InitialiseGrid method, and the reviewed createGridQC method.
After the wave 38 apply checkpoint, func_prot_s4_objects.R measures 288 lines while the live func_prot_s4_qc_methods.R target carries 1287 exact-source lines after adding filterMinNumPeptidesPerProtein; the existing DESCRIPTION Collate order already loaded the QC methods helper before func_prot_s4_objects.R, so no additional load-order sync was required.
No additional live method blocks remain in func_prot_s4_objects.R; treat the file as a stabilized breadcrumb/generics shell for this target and move later work to a different backlog target if further bucket work is needed.
