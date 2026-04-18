#'@export
setMethod( f = "ruvIII_C_Varying"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL) {
             
             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id
             
             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)
             k <- checkParamsObjectFunctionSimplify( theObject, "ruv_number_k", NULL)
             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)
             
             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
             theObject <- updateParamInObject(theObject, "ruv_number_k")
             theObject <- updateParamInObject(theObject, "ctrl")

             # Create matrix from peptide_data (exactly like protein version)
             current_quant_column <- theObject@norm_quantity_column
             
             # Create unique peptide identifier for matrix rownames
             peptide_unique_id <- paste(peptide_data[[theObject@protein_id_column]], 
                                       peptide_data[[theObject@peptide_sequence_column]], 
                                       sep = "%")
             
             # Create wide format data frame then convert to matrix (exactly like protein version)
             temp_peptide_wide <- peptide_data |>
               mutate(peptide_id = peptide_unique_id) |>
               select(peptide_id, !!sym(sample_id), !!sym(current_quant_column)) |>
               pivot_wider(names_from = !!sym(sample_id), 
                          values_from = !!sym(current_quant_column),
                          values_fn = mean)

             # Convert to matrix exactly like protein version: column_to_rownames() then as.matrix()
             normalised_frozen_peptide_matrix_filt <- temp_peptide_wide |>
               column_to_rownames("peptide_id") |>
               as.matrix()
             
             Y <-  t( normalised_frozen_peptide_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id))])
             
             M <- getRuvIIIReplicateMatrixHelper( design_matrix
                                                  , !!sym(sample_id)
                                                  , !!sym(ruv_grouping_variable))
             
             cln_mat <- RUVIII_C_Varying( k = ruv_number_k
                                          , Y = Y
                                          , M = M
                                          , toCorrect = colnames(Y)
                                          , potentialControls = names( ctrl[which(ctrl)] ) )
             
             # Remove samples with no values
             cln_mat_2 <- cln_mat[rowSums(is.na(cln_mat) | is.nan(cln_mat)) != ncol(cln_mat),]
             
             # Remove proteins with no values
             cln_mat_3 <- t(cln_mat_2)
             cln_mat_4 <- cln_mat_3[rowSums(is.na(cln_mat_3) | is.nan(cln_mat_3)) != ncol(cln_mat_3),]
             
             # Convert back to data frame (exactly like protein version)
             ruv_normalised_results_cln <- cln_mat_4 |>
               as.data.frame() |>
               rownames_to_column("peptide_id")

             # Update peptide_data by joining with RUV-corrected values and replacing the quant column
             updated_peptide_data <- peptide_data |>
               mutate(peptide_id = peptide_unique_id) |>
               select(-!!sym(current_quant_column)) |>
               left_join(ruv_normalised_results_cln |> 
                        pivot_longer(cols = -peptide_id, 
                                    names_to = sample_id, 
                                    values_to = current_quant_column),
                        by = c("peptide_id", sample_id)) |>
               select(-peptide_id)

             # Update both slots
             theObject@peptide_data <- updated_peptide_data
             theObject@peptide_matrix <- cln_mat_4
             
             theObject <- cleanDesignMatrixPeptide(theObject)
             
             return( theObject )
           })

