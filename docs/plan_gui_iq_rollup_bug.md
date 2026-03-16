# Plan: Fix GUI IQ Rollup Bug

## Problem 
In the R/Shiny GUI proteomics workflow, the module `mod_prot_qc_protein_rollup.R` performs peptide-to-protein rollup using the `iq` package. It writes the `peptide_s4` long data to a TSV, runs `iq::process_long_format()`, and then reads the result back with `vroom::vroom()`. It then tries to create a `ProteinQuantitativeData` S4 object. 

The S4 object's validity check (`func_prot_s4_objects.R`) throws the error:  
`"Samples in protein data and design matrix must be the same"`  
when processing the neurolincs parquet dataset.

This happens for two possible reasons:
1. **Samples dropped by `iq`**: `iq` filters out peptides not meeting Q-value thresholds (`filter_double_less = c("Q.Value" = "0.01", ...)`). If a sample has no remaining valid proteins after IQ processing/filtering, it is dropped from the wide-format output. However, the `mod_prot_qc_protein_rollup.R` script passes the *original* `peptide_s4@design_matrix` directly to `ProteinQuantitativeData()`.
2. **Column Name Mangling**: `vroom::vroom()` defaults to `.name_repair = "unique"`. If the sample names are numeric or contain special characters, `vroom` or `iq` might alter the column names, leading to a mismatch with the exact strings in the design matrix `Run` column.

## Proposed Changes

Modify `R/mod_prot_qc_protein_rollup.R` to:

1. Read the `iq` output file with `.name_repair = "minimal"` to prevent unintended formatting changes by `vroom` to sample names.
2. Interrogate the `protein_log2_quant` columns to identify exactly which samples survived the `iq` rollup.
3. Filter the `peptide_s4@design_matrix` to keep only the samples present in the `iq` output, and match their names exactly.
4. Add robust logging to report if any samples were completely dropped by the IQ algorithm.

### Specific Code Changes (`R/mod_prot_qc_protein_rollup.R`)

```r
        # Read IQ output - use .name_repair = "minimal" to preserve sample IDs
        protein_log2_quant <- vroom::vroom(iq_output_file, .name_repair = "minimal")
        
        logger::log_info("Protein Processing: Aligning design matrix with IQ output")
        
        # Identify samples in the IQ output
        iq_samples <- setdiff(colnames(protein_log2_quant), "Protein.Ids")
        
        # Identify if any samples were dropped
        original_samples <- peptide_s4@design_matrix |> dplyr::pull("Run")
        dropped_samples <- setdiff(original_samples, iq_samples)
        
        if (length(dropped_samples) > 0) {
          msg <- sprintf("IQ rollup dropped %d samples entirely due to lack of valid protein data: %s", 
                         length(dropped_samples), paste(dropped_samples, collapse = ", "))
          logger::log_info(msg)
          shiny::showNotification(msg, type = "warning", duration = 10)
        }
        
        # Filter design matrix to match the surviving samples
        aligned_design_matrix <- peptide_s4@design_matrix |>
          dplyr::filter(!!sym("Run") %in% iq_samples)
        
        logger::log_info("Protein Processing: Creating ProteinQuantitativeData S4 object")
        
        # Create ProteinQuantitativeData S4 object
        protein_obj <- ProteinQuantitativeData(
          protein_quant_table = protein_log2_quant,
          protein_id_column = "Protein.Ids",
          protein_id_table = protein_log2_quant |> dplyr::distinct(Protein.Ids),
          design_matrix = aligned_design_matrix,  # Use aligned design matrix
          sample_id = "Run",
          group_id = "group",
          technical_replicate_id = "replicates",
          args = peptide_s4@args
        )
```

## Verification
- We can load the neurolincs parquet test dataset.
- Launch the GUI and trace the steps to the protein rollup stage.
- Run the actual `iq` rollup with the dataset.
- Ensure the `ProteinQuantitativeData()` object is created successfully without the validity check exception.
