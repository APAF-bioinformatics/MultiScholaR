#'@export
setMethod(f = "chooseBestProteinAccession"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject, delim=NULL, seqinr_obj=NULL
                              , seqinr_accession_column=NULL
                              , replace_zero_with_na = NULL # Note: this param is kept for signature consistency but not used at peptide level
                              , aggregation_method = NULL) {

            peptide_data <- theObject@peptide_data
            protein_id_column <- theObject@protein_id_column
            is_logged <- theObject@is_logged_data
            verbose <- TRUE # Assuming verbose is desired, can be made a parameter

            # --- Parameter Handling ---
            delim <- checkParamsObjectFunctionSimplify(theObject, "delim",  default_value =  ";")
            seqinr_obj <- checkParamsObjectFunctionSimplify(theObject, "seqinr_obj",  default_value = NULL)
            seqinr_accession_column <- checkParamsObjectFunctionSimplify(theObject
                                                                       , "seqinr_accession_column"
                                                                       , default_value = "uniprot_acc")
            aggregation_method <- checkParamsObjectFunctionSimplify(theObject
                                                                  , "aggregation_method"
                                                                  , default_value = "mean")

            if (!aggregation_method %in% c("sum", "mean", "median")) {
              stop("aggregation_method must be one of: 'sum', 'mean', 'median'")
            }

            theObject <- updateParamInObject(theObject, "delim")
            theObject <- updateParamInObject(theObject, "seqinr_obj")
            theObject <- updateParamInObject(theObject, "seqinr_accession_column")
            theObject <- updateParamInObject(theObject, "aggregation_method")

            if (verbose) {
              log_info("Choosing best protein accession at the peptide level...")
              log_info("Aggregation method for duplicate peptides will be: '{aggregation_method}'")
            }

            # --- Robustly handle list input for seqinr_obj ---
            if (is.list(seqinr_obj) && !is.data.frame(seqinr_obj) && "aa_seq_tbl_final" %in% names(seqinr_obj)) {
              if (verbose) log_info("Extracting data frame from seqinr_obj list...")
              seqinr_obj <- seqinr_obj$aa_seq_tbl_final
            }


            # --- Map Old to New IDs ---
            accession_mapping <- chooseBestProteinAccessionHelper(input_tbl = peptide_data,
                                                                  acc_detail_tab = seqinr_obj,
                                                                  accessions_column = !!sym(protein_id_column),
                                                                  row_id_column = seqinr_accession_column,
                                                                  group_id = !!sym(protein_id_column),
                                                                  delim = delim)

            # --- Update Protein IDs ---
            updated_peptide_data <- peptide_data |>
              dplyr::left_join(accession_mapping |> dplyr::select(!!sym(protein_id_column), !!sym(seqinr_accession_column)),
                               by = setNames(protein_id_column, protein_id_column)) |>
              dplyr::select(-!!sym(protein_id_column)) |>
              dplyr::rename(!!sym(protein_id_column) := !!sym(seqinr_accession_column))

            # --- Aggregate Duplicates ---
            if (verbose) {
              log_info("Aggregating peptide quantities for newly resolved protein groups...")
            }

            grouping_cols <- c(theObject@sample_id, protein_id_column, theObject@peptide_sequence_column)
            quant_cols <- c(theObject@raw_quantity_column, theObject@norm_quantity_column)
            meta_cols <- setdiff(colnames(updated_peptide_data), c(grouping_cols, quant_cols))

            aggregated_data <- updated_peptide_data |>
              dplyr::group_by(across(all_of(grouping_cols))) |>
              dplyr::summarise(
                # Aggregate quantitative columns based on the chosen method
                across(all_of(quant_cols), function(x) {
                  # Remove NAs for calculation
                  x_clean <- x[!is.na(x)]
                  if (length(x_clean) == 0) return(NA_real_)

                  # Perform aggregation
                  if (aggregation_method == "mean") {
                    mean(x_clean)
                  } else if (aggregation_method == "median") {
                    median(x_clean)
                  } else if (aggregation_method == "sum") {
                    if (is_logged) {
                      log2(sum(2^x_clean)) # Sum on linear scale, then log2 transform back
                    } else {
                      sum(x_clean)
                    }
                  }
                }),
                # Keep the first value for all metadata columns
                across(all_of(meta_cols), dplyr::first),
                .groups = "drop"
              )

            if (verbose) {
              log_info("Aggregation complete. Before: {nrow(updated_peptide_data)} rows. After: {nrow(aggregated_data)} rows.")
            }

            theObject@peptide_data <- aggregated_data

            # --- Recalculate Matrix ---
            if (verbose) {
              log_info("Recalculating peptide matrix with updated protein accessions...")
            }
            theObject <- calcPeptideMatrix(theObject)

            return(theObject)
          })
