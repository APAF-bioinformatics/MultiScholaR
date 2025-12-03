# ============================================================================
# func_prot_s4_objects.R
# ============================================================================
# Purpose: Proteomics S4 class definitions and methods
# 
# This file contains S4 class definitions and methods specific to proteomics,
# including PeptideQuantitativeData, ProteinQuantitativeData classes,
# their constructors, and proteomics-specific S4 methods.
#
# Consolidated from:
# - peptideVsSamplesS4Objects.R (PeptideQuantitativeData class + 23 methods)
# - proteinVsSamplesS4Objects.R (ProteinQuantitativeData method)
# - protein_de_analysis_wrapper.R (DE analysis methods)
#
# Dependencies:
# - methods package
# - func_general_s4_generics.R (for generic definitions)
# ============================================================================

# ==========================================
# ProteinQuantitativeData S4 Class
# ==========================================

#' ProteinQuantitativeData S4 Class
#'
#' @description
#' An S4 class to store and manage protein-level quantitative data
#' from proteomics experiments along with experimental design metadata.
#'
#' @slot protein_quant_table A data.frame containing protein quantification data
#' @slot protein_id_column Character string naming the protein ID column
#' @slot design_matrix A data.frame containing experimental design
#' @slot protein_id_table A data.frame containing protein ID information
#' @slot sample_id Character string naming the sample ID column
#' @slot group_id Character string naming the group column
#' @slot technical_replicate_id Character string naming the replicate column
#' @slot args A list of additional arguments
#'
#' @exportClass ProteinQuantitativeData
ProteinQuantitativeData <- setClass("ProteinQuantitativeData"

         , slots = c(

                      # Protein vs Sample quantitative data

                      protein_quant_table = "data.frame"

                      , protein_id_column = "character"



                      # Design Matrix Information

                      , design_matrix = "data.frame"

                      , protein_id_table = "data.frame"

                      , sample_id="character"

                      , group_id="character"

                      , technical_replicate_id="character"

                      , args = "list")



         , prototype = list(

           # Protein vs Sample quantitative data

           protein_id_column = "Protein.Ids"

          , protein_id_table = data.frame()

           # Design Matrix Information

           , sample_id="Sample_id"

           , group_id="group"

           , technical_replicate_id="replicates"

           , args = NULL

           )



         , validity = function(object) {

           if( !is.data.frame(object@protein_quant_table) ) {

             stop("protein_quant_table must be a data.frame")

           }

           if( !is.character(object@protein_id_column) ) {

             stop("protein_id_column must be a character")

           }

           if( !is.data.frame(object@design_matrix) ) {

             stop("design_matrix must be a data.frame")

           }

           if( !is.character(object@sample_id) ) {

             stop("sample_id must be a character")

           }

           if( !is.character(object@group_id) ) {

             stop("group_id must be a character")

           }

           if( !is.character(object@technical_replicate_id) ) {

             stop("technical_replicate_id must be a character")

           }



            if( ! object@protein_id_column %in% colnames(object@protein_quant_table) ) {

                stop("Protein ID column must be in the protein data table")

            }



           if( ! object@sample_id %in% colnames(object@design_matrix) ) {

             stop("Sample ID column must be in the design matrix")

           }





           #Need to check the rows names in design matrix and the column names of the data table

           samples_in_protein_quant_table <- setdiff(colnames( object@protein_quant_table), object@protein_id_column)

           samples_in_design_matrix <- object@design_matrix |> dplyr::pull( !! sym( object@sample_id ) )



           if( length( which( sort(samples_in_protein_quant_table) != sort(samples_in_design_matrix) )) > 0 ) {

             stop("Samples in protein data and design matrix must be the same" )

           }



         }



)



# Initialize method to ensure design_matrix's sample_id column is character

setMethod("initialize", "ProteinQuantitativeData",

  function(.Object, ...) {

    # Capture all arguments passed to the constructor

    args_list <- list(...)



    # If design_matrix and sample_id are provided in the arguments

    if (!is.null(args_list$design_matrix) && !is.null(args_list$sample_id)) {

      if (args_list$sample_id %in% names(args_list$design_matrix)) {

        args_list$design_matrix[[args_list$sample_id]] <- as.character(args_list$design_matrix[[args_list$sample_id]])

        logger::log_debug(

          "Initialize ProteinQuantitativeData: Converted column '{args_list$sample_id}' in design_matrix to character."

        )

      } else {

        # This case should ideally be caught by validity, but good to be aware

        logger::log_warn(

          "Initialize ProteinQuantitativeData: Sample ID column '{args_list$sample_id}' not found in the provided design_matrix during initialization."

        )

      }

    }

    

    # Reconstruct the call to the default initializer with potentially modified args_list

    # This uses do.call to pass arguments correctly after modification

    .Object <- do.call(callNextMethod, c(list(.Object), args_list))

    

    # Explicitly return the object

    return(.Object)

  }

)

# ==========================================
# Content from proteinVsSamplesS4Objects.R
# ==========================================
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Protein-Level Missing Value Imputation using limpa Package
#'
#' This function applies limpa's DPC-based missing value imputation directly to 
#' protein-level quantification data. This is useful when you already have protein
#' quantification data with missing values that need to be handled.
#'
#' @param theObject A ProteinQuantitativeData object with protein-level data
#' @param dpc_results DPC results to use. If NULL, will estimate using dpc_slope
#' @param dpc_slope Default DPC slope to use if no DPC results available (default: 0.8)
#' @param quantified_protein_column Name for the column containing quantified protein values
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param chunk When verbose=TRUE, how often to output progress information (default: 1000)
#'
#' @details
#' This method treats each protein as a separate "feature" and applies DPC-based
#' imputation using limpa's dpcImpute function. This is appropriate when you have
#' protein-level data with missing values that follow intensity-dependent patterns.
#'
#' @return Updated ProteinQuantitativeData object with imputed protein values
#'
#' @export
setMethod(f="proteinMissingValueImputationLimpa"
          , signature="ProteinQuantitativeData"
          , definition = function(theObject, 
                                  dpc_results = NULL,
                                  dpc_slope = 0.8,
                                  quantified_protein_column = NULL,
                                  verbose = TRUE,
                                  chunk = 1000) {
            
            # Load required packages
            if (!requireNamespace("limpa", quietly = TRUE)) {
              stop("limpa package is required but not installed. Please install it using: BiocManager::install('limpa')")
            }
            
            # Parameter validation and defaults
            quantified_protein_column <- if (is.null(quantified_protein_column)) {
              "Protein.Imputed.Limpa"
            } else {
              quantified_protein_column
            }
            
            # Extract data from protein object
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            sample_id_column <- theObject@sample_id
            design_matrix <- theObject@design_matrix
            
            if (verbose) {
              log_info("Starting limpa-based protein-level missing value imputation...")
            }
            
            # Convert to matrix format (proteins as rows, samples as columns)
            # First, identify sample columns (exclude protein ID and other metadata)
            sample_columns <- setdiff(colnames(protein_quant_table), 
                                     c(protein_id_column, "description", "gene_name", 
                                       "protein_name", "organism", "length"))
            
            if (verbose) {
              log_info("Converting protein data to matrix format...")
              log_info("Found {length(sample_columns)} sample columns")
            }
            
            # Create protein matrix
            protein_matrix <- protein_quant_table |>
              dplyr::select(all_of(c(protein_id_column, sample_columns))) |>
              tibble::column_to_rownames(protein_id_column) |>
              as.matrix()

            if (verbose) {
              log_info("Protein matrix dimensions: {nrow(protein_matrix)} proteins x {ncol(protein_matrix)} samples")
              log_info("Missing value percentage: {round(100 * mean(is.na(protein_matrix)), 1)}%")
            }
            
            # Check if we need log2 transformation
            # Assume if max value > 50, data is not log2 transformed
            max_val <- max(protein_matrix, na.rm = TRUE)
            needs_log_transform <- max_val > 50
            
            if (needs_log_transform) {
              if (verbose) {
                log_info("Converting to log2 scale for limpa (max value: {round(max_val, 2)})...")
              }
              protein_matrix <- log2(protein_matrix + 1)
            } else {
              if (verbose) {
                log_info("Data appears to be log2-scale already (max value: {round(max_val, 2)})")
              }
            }
            
            # Handle infinite or NaN values
            if (any(is.infinite(protein_matrix) | is.nan(protein_matrix), na.rm = TRUE)) {
              if (verbose) {
                log_warn("Infinite or NaN values detected. Replacing with NA...")
              }
              protein_matrix[is.infinite(protein_matrix) | is.nan(protein_matrix)] <- NA
            }
            
            # Get or estimate DPC parameters
            dpc_params <- NULL
            if (!is.null(dpc_results)) {
              dpc_params <- dpc_results
              if (verbose) {
                log_info("Using provided DPC results")
              }
            } else {
              if (verbose) {
                log_info("Estimating DPC from protein data...")
              }
            tryCatch({
                dpcfit <- limpa::dpc(protein_matrix)
                dpc_params <- dpcfit
                if (verbose) {
                  log_info("DPC parameters estimated:")
                  log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                  log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
                }
            }, error = function(e) {
                if (verbose) {
                  log_warn("DPC estimation failed, using default slope: {dpc_slope}")
                }
                dpc_params <- NULL
              })
            }
            
            # Apply missing value imputation
            if (verbose) {
              log_info("Applying protein-level missing value imputation...")
            }
            
              tryCatch({
              # Apply dpcImpute to protein matrix
              if (!is.null(dpc_params)) {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc = dpc_params, verbose = verbose, chunk = chunk)
            } else {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc.slope = dpc_slope, verbose = verbose, chunk = chunk)
              }
              
              if (verbose) {
                log_info("Protein-level imputation completed successfully")
                log_info("No missing values remaining: {!any(is.na(imputed_result$E))}")
              }
              
              # Extract imputed matrix
              imputed_matrix <- imputed_result$E
              
              # Transform back to original scale if necessary
              if (needs_log_transform) {
                if (verbose) {
                  log_info("Converting back from log2 scale...")
                }
                imputed_matrix <- 2^imputed_matrix - 1
                # Ensure no negative values
                imputed_matrix[imputed_matrix < 0] <- 0
              }
              
              # Convert back to long format
              if (verbose) {
                log_info("Converting imputed data back to original format...")
              }
              
              imputed_long <- imputed_matrix |>
                      as.data.frame() |>
                tibble::rownames_to_column(protein_id_column) |>
                tidyr::pivot_longer(cols = -all_of(protein_id_column), 
                                   names_to = sample_id_column, 
                                   values_to = quantified_protein_column)
              
              # Merge with original protein data
              updated_protein_data <- protein_quant_table |>
                dplyr::left_join(imputed_long, by = c(protein_id_column, sample_id_column))
              
              # Update the object
              theObject@protein_quant_table <- updated_protein_data
              
              # Store DPC results for future reference
              if (is.null(theObject@args)) {
                theObject@args <- list()
              }
              
              theObject@args$limpa_protein_imputation_results <- list(
                dpc_parameters_used = if (!is.null(dpc_params)) {
                  if (is.list(dpc_params) && !is.null(dpc_params$dpc)) {
                    dpc_params$dpc  # Extract parameters from DPC object
                  } else if (is.numeric(dpc_params)) {
                    dpc_params  # Already numeric parameters  
             } else {
                    c(NA, dpc_slope)
                  }
                } else {
                  c(NA, dpc_slope)
                },
                dpc_object_used = if (is.list(dpc_params) && !is.null(dpc_params$dpc)) dpc_params else NULL,
                quantified_protein_column = quantified_protein_column,
                missing_percentage_before = round(100 * mean(is.na(protein_matrix)), 1),
                missing_percentage_after = 0,  # DPC imputation produces complete data
                imputation_method = "limpa_dpc_protein_imputation",
                total_proteins_imputed = nrow(imputed_matrix)
              )
              
              if (verbose) {
                log_info("limpa protein-level imputation completed successfully!")
                log_info("New imputed column: {quantified_protein_column}")
                log_info("DPC results stored in object@args$limpa_protein_imputation_results")
              }
            
            return(theObject)
              
              }, error = function(e) {
              log_error(paste("Error during limpa protein imputation:", e$message))
              stop(paste("limpa protein imputation failed:", e$message))
            })
          })

# ==========================================
# Content from protein_de_analysis_wrapper.R
# ==========================================
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Differential Expression Analysis for Proteomics
#' 
#' This file contains modular functions for proteomics differential expression analysis,
#' breaking up the monolithic deAnalysisWrapperFunction into focused components.
#' 
#' Functions:
#' - differentialExpressionAnalysis: Main S4 generic for DE analysis
#' - differentialExpressionAnalysisHelper: Core limma-based analysis
#' - generateVolcanoPlotGlimma: Interactive volcano plot generation
#' - generateDEHeatmap: Heatmap visualization with clustering
#' 
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#' @export
setMethod( f ="differentialExpressionAnalysis"
           , signature = "ProteinQuantitativeData"
           , definition=function( theObject
                                  , contrasts_tbl = NULL
                                  , formula_string = NULL
                                  , group_id = NULL
                                  , de_q_val_thresh = NULL
                                  , treat_lfc_cutoff = NULL
                                  , eBayes_trend = NULL
                                  , eBayes_robust = NULL
                                  , args_group_pattern = NULL
                                  , args_row_id = NULL
                                  , qvalue_column = "fdr_qvalue"
                                  , raw_pvalue_column = "raw_pvalue" ) {

  # IMMEDIATE ERROR CATCH - Check if we even get here
  message("*** ENTERING differentialExpressionAnalysis METHOD ***")
  message(sprintf("*** METHOD SIGNATURE MATCHED: ProteinQuantitativeData ***"))
  
  # Try to catch the index error immediately
  tryCatch({
    message("*** Testing parameter access ***")
    if (!is.null(contrasts_tbl)) {
      test <- contrasts_tbl[[1]]
      message("*** Parameter access successful ***")
    }
  }, error = function(e) {
    message(sprintf("*** IMMEDIATE ERROR: %s ***", e$message))
    message(sprintf("*** ERROR CLASS: %s ***", class(e)))
    stop(e)
  })

  message("--- Entering differentialExpressionAnalysis ---")
  message(sprintf("   differentialExpressionAnalysis: theObject class = %s", class(theObject)))
  message(sprintf("   differentialExpressionAnalysis: contrasts_tbl provided = %s", !is.null(contrasts_tbl)))
  if(!is.null(contrasts_tbl)) {
    message(sprintf("   differentialExpressionAnalysis: contrasts_tbl dims = %d x %d", nrow(contrasts_tbl), ncol(contrasts_tbl)))
    message(sprintf("   differentialExpressionAnalysis: contrasts_tbl content = %s", paste(contrasts_tbl[[1]], collapse=", ")))
  }

  # Wrap the helper function call in tryCatch to get better error info
  message("   differentialExpressionAnalysis: About to call differentialExpressionAnalysisHelper...")
  
  results_list <- tryCatch({
    differentialExpressionAnalysisHelper(  theObject
                                          , contrasts_tbl = contrasts_tbl
                                          , formula_string = formula_string
                                          , group_id = group_id
                                          , de_q_val_thresh = de_q_val_thresh
                                          , treat_lfc_cutoff = treat_lfc_cutoff
                                          , eBayes_trend = eBayes_trend
                                          , eBayes_robust = eBayes_robust
                                          , args_group_pattern = args_group_pattern
                                          , args_row_id = args_row_id
                                          , qvalue_column = qvalue_column
                                          , raw_pvalue_column = raw_pvalue_column
                                          )
  }, error = function(e) {
    # CRITICAL FIX: Use paste() for logger calls in error handlers to avoid interpolation bug
    message(paste("   differentialExpressionAnalysis ERROR in helper function:", e$message))
    message(paste("   differentialExpressionAnalysis ERROR call stack:", capture.output(traceback())))
    stop(e)
  })

  message("   differentialExpressionAnalysis: Helper function completed successfully!")
  message("--- Exiting differentialExpressionAnalysis ---")
  return(results_list)

})

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
setMethod( f ="differentialExpressionAnalysisHelper"
           , signature = "ProteinQuantitativeData"
           , definition=function( theObject
                                  , contrasts_tbl = NULL
                                  , formula_string = NULL
                                  , group_id = NULL
                                  , de_q_val_thresh = NULL
                                  , treat_lfc_cutoff = NULL
                                  , eBayes_trend = NULL
                                  , eBayes_robust = NULL
                                  , args_group_pattern = NULL
                                  , args_row_id = NULL
                                  , qvalue_column = "fdr_qvalue"
                                  , raw_pvalue_column = "raw_pvalue" ) {

  message("--- Entering differentialExpressionAnalysisHelper ---")
  
  # IMMEDIATE PARAMETER VALIDATION TO CATCH INDEX ERROR
  message("   PARAMETER VALIDATION Step: Checking all input parameters...")
  message(sprintf("      theObject is NULL: %s", is.null(theObject)))
  message(sprintf("      contrasts_tbl is NULL: %s", is.null(contrasts_tbl)))
  if (!is.null(contrasts_tbl)) {
    message(sprintf("      contrasts_tbl class: %s", class(contrasts_tbl)))
    message(sprintf("      contrasts_tbl type: %s", typeof(contrasts_tbl)))
    message("      contrasts_tbl print:")
    print(contrasts_tbl)
    message("      Trying to access contrasts_tbl[[1]]...")
    tryCatch({
      first_col <- contrasts_tbl[[1]]
      message(sprintf("      SUCCESS: contrasts_tbl[[1]] class: %s", class(first_col)))
      message(sprintf("      SUCCESS: contrasts_tbl[[1]] length: %d", length(first_col)))
      message(sprintf("      SUCCESS: contrasts_tbl[[1]] content: %s", paste(first_col, collapse=", ")))
    }, error = function(e) {
      message(sprintf("      ERROR accessing contrasts_tbl[[1]]: %s", e$message))
    })
  }
  
  message("   DEBUG66: Initial parameter inspection")
  message(sprintf("      DEBUG66: theObject class = %s", class(theObject)))
  message(sprintf("      DEBUG66: contrasts_tbl param is.null = %s", is.null(contrasts_tbl)))
  if (!is.null(contrasts_tbl)) {
    message(sprintf("      DEBUG66: contrasts_tbl param class = %s", class(contrasts_tbl)))
    message("      DEBUG66: contrasts_tbl param structure:")
    str(contrasts_tbl)
  }

  message("   differentialExpressionAnalysisHelper Step: Extracting parameters...")
  # Extract parameters from S4 object with fallbacks
  contrasts_tbl <- checkParamsObjectFunctionSimplify( theObject, "contrasts_tbl", contrasts_tbl)
  formula_string <- checkParamsObjectFunctionSimplify( theObject, "formula_string", "~ 0 + group")
  group_id <- checkParamsObjectFunctionSimplify( theObject, "group_id", "group")
  de_q_val_thresh <- checkParamsObjectFunctionSimplify( theObject, "de_q_val_thresh", 0.05)
  treat_lfc_cutoff <- checkParamsObjectFunctionSimplify( theObject, "treat_lfc_cutoff", 0)
  eBayes_trend <- checkParamsObjectFunctionSimplify( theObject, "eBayes_trend", TRUE)
  eBayes_robust <- checkParamsObjectFunctionSimplify( theObject, "eBayes_robust", TRUE)
  args_group_pattern <- checkParamsObjectFunctionSimplify( theObject, "args_group_pattern", "(\\d+)")
  args_row_id <- checkParamsObjectFunctionSimplify( theObject, "args_row_id", "uniprot_acc")
  
  message("   DEBUG66: After parameter extraction")
  message(sprintf("      DEBUG66: contrasts_tbl class = %s", class(contrasts_tbl)))
  message("      DEBUG66: contrasts_tbl structure:")
  str(contrasts_tbl)

  message(sprintf("   differentialExpressionAnalysisHelper: formula_string = %s", formula_string))
  message(sprintf("   differentialExpressionAnalysisHelper: group_id = %s", group_id))
  message(sprintf("   differentialExpressionAnalysisHelper: de_q_val_thresh = %f", de_q_val_thresh))

  message("   differentialExpressionAnalysisHelper Step: Processing group names...")
  # Handle group names that start with numbers (same pattern as original wrapper)
  design_matrix <- theObject@design_matrix
  group_col <- design_matrix[[group_id]]
  message(sprintf("   differentialExpressionAnalysisHelper: group_col length = %d", length(group_col)))
  message(sprintf("   differentialExpressionAnalysisHelper: group_col content = %s", paste(head(group_col, 10), collapse=", ")))
  
  # Check if any group names start with numbers and create mapping
  starts_with_number <- grepl("^[0-9]", group_col)
  message(sprintf("   differentialExpressionAnalysisHelper: any start with number? %s", any(starts_with_number)))
  
  if(any(starts_with_number)) {
    message("   differentialExpressionAnalysisHelper Step: Fixing group names that start with numbers...")
    
    original_groups <- unique(group_col)
    message(sprintf("   differentialExpressionAnalysisHelper: original_groups = %s", paste(original_groups, collapse=", ")))
    message(sprintf("   differentialExpressionAnalysisHelper: About to process %d original groups with purrr::map_chr", length(original_groups)))
    
    tryCatch({
      message("      DEBUG66: Before safe_groups purrr::map_chr")
      message(sprintf("      DEBUG66: original_groups class = %s", class(original_groups)))
      message(sprintf("      DEBUG66: original_groups length = %d", length(original_groups)))
      message("      DEBUG66: original_groups content:")
      print(original_groups)
      
      safe_groups <- purrr::map_chr(original_groups, \(x) {
        message(sprintf("         DEBUG66: Processing group item: '%s' (class: %s)", x, class(x)))
        result <- if(grepl("^[0-9]", x)) paste0("grp_", x) else x
        message(sprintf("         DEBUG66: Result for '%s' -> '%s'", x, result))
        result
      })
      message("   differentialExpressionAnalysisHelper: safe_groups processing SUCCESS")
      message("      DEBUG66: safe_groups result:")
      print(safe_groups)
    }, error = function(e) {
      message(sprintf("   differentialExpressionAnalysisHelper ERROR in safe_groups purrr::map_chr: %s", e$message))
      message("      DEBUG66: Error details:")
      message(sprintf("      DEBUG66: Error class: %s", class(e)))
      print(str(e))
      stop(e)
    })
    
    group_mapping <- setNames(original_groups, safe_groups)
    message("   differentialExpressionAnalysisHelper: group_mapping created")
    
    # Update design matrix with safe names
    message("   differentialExpressionAnalysisHelper: About to update design matrix with purrr::map_chr")
    tryCatch({
      design_matrix[[group_id]] <- purrr::map_chr(group_col, \(x) {
        if(grepl("^[0-9]", x)) paste0("grp_", x) else x
      })
      message("   differentialExpressionAnalysisHelper: design matrix update SUCCESS")
    }, error = function(e) {
      message(sprintf("   differentialExpressionAnalysisHelper ERROR in design matrix purrr::map_chr: %s", e$message))
      stop(e)
    })
    
    # Update contrasts table if it exists
    if(!is.null(contrasts_tbl)) {
      message("   differentialExpressionAnalysisHelper: About to update contrasts table with purrr::map_chr")
      message(sprintf("   differentialExpressionAnalysisHelper: contrasts_tbl[[1]] = %s", paste(contrasts_tbl[[1]], collapse=", ")))
      
      message("      DEBUG66: Contrasts table inspection before purrr::map_chr")
      message(sprintf("      DEBUG66: contrasts_tbl class = %s", class(contrasts_tbl)))
      message(sprintf("      DEBUG66: contrasts_tbl dim = %d x %d", nrow(contrasts_tbl), ncol(contrasts_tbl)))
      message("      DEBUG66: contrasts_tbl structure:")
      str(contrasts_tbl)
      message(sprintf("      DEBUG66: contrasts_tbl[[1]] class = %s", class(contrasts_tbl[[1]])))
      message(sprintf("      DEBUG66: contrasts_tbl[[1]] length = %d", length(contrasts_tbl[[1]])))
      message("      DEBUG66: group_mapping content:")
      print(group_mapping)
      
      tryCatch({
        contrasts_tbl[[1]] <- purrr::map_chr(contrasts_tbl[[1]], \(x) {
          message(sprintf("         DEBUG66: Processing contrast: '%s' (class: %s)", x, class(x)))
          message(sprintf("         DEBUG66: group_mapping names: %s", paste(names(group_mapping), collapse=", ")))
          result <- x
          for(orig in names(group_mapping)) {
            message(sprintf("            DEBUG66: Checking gsub: '%s' -> '%s'", group_mapping[orig], orig))
            result <- gsub(group_mapping[orig], orig, result, fixed = TRUE)
          }
          message(sprintf("         DEBUG66: Final result: '%s' -> '%s'", x, result))
          result
        })
        message("   differentialExpressionAnalysisHelper: contrasts table update SUCCESS")
        message("      DEBUG66: Updated contrasts_tbl[[1]]:")
        print(contrasts_tbl[[1]])
      }, error = function(e) {
        message(sprintf("   differentialExpressionAnalysisHelper ERROR in contrasts table purrr::map_chr: %s", e$message))
        message("      DEBUG66: Error details:")
        message(sprintf("      DEBUG66: Error class: %s", class(e)))
        print(str(e))
        stop(e)
      })
    }
    
    theObject@design_matrix <- design_matrix
  }

  message("   differentialExpressionAnalysisHelper Step: Updating S4 object parameters...")
  # Update object with validated parameters
  theObject <- updateParamInObject(theObject, "contrasts_tbl")
  theObject <- updateParamInObject(theObject, "formula_string")
  theObject <- updateParamInObject(theObject, "group_id")
  theObject <- updateParamInObject(theObject, "de_q_val_thresh")
  theObject <- updateParamInObject(theObject, "treat_lfc_cutoff")
  theObject <- updateParamInObject(theObject, "eBayes_trend")
  theObject <- updateParamInObject(theObject, "eBayes_robust")
  theObject <- updateParamInObject(theObject, "args_group_pattern")
  theObject <- updateParamInObject(theObject, "args_row_id")

  return_list <- list()
  return_list$theObject <- theObject

  message("   differentialExpressionAnalysisHelper Step: Generating QC plots...")

  # Generate QC plots (RLE, PCA, density)
  rle_plot <- plotRle(theObject = theObject, theObject@group_id) +
    theme(axis.text.x = element_text(size = 13)) +
    theme(axis.text.y = element_text(size = 13)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12)) +
    xlab("Samples")

  return_list$rle_plot <- rle_plot

  pca_plot <- plotPca( theObject
                      , grouping_variable = theObject@group_id
                      , label_column = ""
                      , title = ""
                      , font_size = 8) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))

  return_list$pca_plot <- pca_plot

  pca_plot_with_labels <- plotPca( theObject
                                  , grouping_variable = theObject@group_id
                                  , label_column = theObject@sample_id
                                  , title = ""
                                  , font_size = 8) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))

  return_list$pca_plot_with_labels <- pca_plot_with_labels

  # Count the number of values
  return_list$plot_num_of_values <- plotNumOfValuesNoLog(theObject@protein_quant_table)

  message("   differentialExpressionAnalysisHelper Step: Running limma contrasts analysis...")

  # Prepare design matrix for limma
  message("   DEBUG66: About to prepare design matrix rownames")
  message(paste("      DEBUG66: theObject@sample_id =", theObject@sample_id))
  message(paste("      DEBUG66: design_matrix dims before =", nrow(theObject@design_matrix), "x", ncol(theObject@design_matrix)))
  message("      DEBUG66: design_matrix column names:")
  print(colnames(theObject@design_matrix))
  
  rownames( theObject@design_matrix ) <- theObject@design_matrix |> dplyr::pull( one_of(theObject@sample_id ))
  
  message("      DEBUG66: design_matrix rownames set successfully")

  # Run the core limma analysis using existing function
  protein_quant_matrix <- as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column))
  contrast_strings_to_use <- contrasts_tbl[, 1]
  message(paste("   differentialExpressionAnalysisHelper: About to call runTestsContrasts with", length(contrast_strings_to_use), "contrasts"))
  message(paste("   differentialExpressionAnalysisHelper: protein_quant_matrix dims =", nrow(protein_quant_matrix), "x", ncol(protein_quant_matrix)))
  
  contrasts_results <- runTestsContrasts(protein_quant_matrix,
                                         contrast_strings = contrast_strings_to_use,
                                         design_matrix = theObject@design_matrix,
                                         formula_string = formula_string,
                                         weights = NA,
                                         treat_lfc_cutoff = as.double(treat_lfc_cutoff),
                                         eBayes_trend = as.logical(eBayes_trend),
                                         eBayes_robust = as.logical(eBayes_robust))

  message("   differentialExpressionAnalysisHelper Step: runTestsContrasts completed successfully!")

  # The result from runTestsContrasts is a LIST with a 'results'
  # element, which ITSELF is a list of tables. For a single contrast run, we 
  # need to extract the first table from the nested list.
  
  # Check if the expected structure exists
  if (is.null(contrasts_results) || is.null(contrasts_results$results) || length(contrasts_results$results) == 0) {
      stop("Error: DE analysis function did not return results in the expected format.")
  }
  
  # Extract the results table (it's the first element in the nested list)
  de_results_table <- contrasts_results$results[[1]]
  
  # Get the name of the contrast from the list element name
  contrast_name <- names(contrasts_results$results)[1]
  
  # Ensure it's a data frame before proceeding
  if (!is.data.frame(de_results_table)) {
    stop("Error: DE analysis results are not in the expected data frame format.")
  }

  # CRITICAL FIX 3.0: The 'topTreat' table needs a 'comparison' column added to it,
  # containing the name of the contrast. It also needs the protein IDs from rownames.
  de_results_table <- de_results_table |>
    tibble::rownames_to_column(var = args_row_id) |>
    dplyr::mutate(comparison = contrast_name)

  # Map back to original group names in results if needed
  if(exists("group_mapping")) {
    contrasts_results_table <- de_results_table |>
      dplyr::mutate(comparison = purrr::map_chr(comparison, \(x) {
        result <- x
        for(safe_name in names(group_mapping)) {
          result <- gsub(safe_name, group_mapping[safe_name], result, fixed = TRUE)
        }
        result
      }))
  } else {
    contrasts_results_table <- de_results_table
  }

  return_list$contrasts_results <- contrasts_results
  return_list$contrasts_results_table <- contrasts_results_table
  
  # Extract qvalue warnings if present
  if (!is.null(contrasts_results$qvalue_warnings) && length(contrasts_results$qvalue_warnings) > 0) {
    return_list$qvalue_warnings <- contrasts_results$qvalue_warnings
    message(sprintf("   differentialExpressionAnalysisHelper Step: qvalue() failed for %d contrast(s)", length(contrasts_results$qvalue_warnings)))
  }

  message("   differentialExpressionAnalysisHelper Step: Preparing data for visualization...")
  message(paste("   DEBUG66: contrast_name for list naming =", contrast_name))

  # Prepare data for volcano plots
  # CRITICAL FIX: getSignificantData expects list_of_de_tables to be a list where each element 
  # is itself a list of data.frames. Since we're processing one contrast at a time, we need to 
  # wrap the single data.frame in an additional list layer.
  # CRITICAL FIX 2: The list element name MUST contain "=" delimiter because countStatDeGenesHelper
  # expects to split on "=" to separate comparison name from expression. Use contrast_name which
  # contains the full format like "H4_vs_WT=groupH4-groupWT"
  nested_list <- list(contrasts_results_table)
  names(nested_list) <- contrast_name  # Use actual contrast name with "=" delimiter
  significant_rows <- getSignificantData(list_of_de_tables = list(nested_list),
                                         list_of_descriptions = list("RUV applied"),
                                         row_id = !!sym(args_row_id),
                                         p_value_column = !!sym(raw_pvalue_column),
                                         q_value_column = !!sym(qvalue_column),
                                         fdr_value_column = fdr_value_bh_adjustment,
                                         log_q_value_column = lqm,
                                         log_fc_column = logFC,
                                         comparison_column = "comparison",
                                         expression_column = "log_intensity",
                                         facet_column = analysis_type,
                                         q_val_thresh = de_q_val_thresh) |>
    dplyr::rename(log2FC = "logFC")

  return_list$significant_rows <- significant_rows

  # Generate volcano plots
  volplot_plot <- plotVolcano(significant_rows,
                              log_q_value_column = lqm,
                              log_fc_column = log2FC,
                              q_val_thresh = de_q_val_thresh,
                              formula_string = "analysis_type ~ comparison")

  return_list$volplot_plot <- volplot_plot

  # Count significant molecules
  # CRITICAL FIX: Same as above - wrap in additional list layer and use contrast_name
  nested_list_for_count <- list(contrasts_results_table)
  names(nested_list_for_count) <- contrast_name  # Use actual contrast name with "=" delimiter
  num_sig_de_molecules <- printCountDeGenesTable(list_of_de_tables = list(nested_list_for_count),
                                                 list_of_descriptions = list("RUV applied"),
                                                 formula_string = "analysis_type ~ comparison")

  return_list$num_sig_de_molecules_first_go <- num_sig_de_molecules

  # P-values distribution plot
  pvalhist <- printPValuesDistribution(significant_rows,
                                       p_value_column = !!sym(raw_pvalue_column),
                                       formula_string = "analysis_type ~ comparison")

  return_list$pvalhist <- pvalhist

  message("   differentialExpressionAnalysisHelper Step: Creating results tables...")

  # Create wide format output
  counts_table_to_use <- theObject@protein_quant_table

  norm_counts <- counts_table_to_use |>
    as.data.frame() |>
    column_to_rownames(args_row_id) |>
    set_colnames(paste0(colnames(counts_table_to_use[-1]), ".log2norm")) |>
    rownames_to_column(args_row_id)

  return_list$norm_counts <- norm_counts

  de_proteins_wide <- significant_rows |>
    dplyr::filter(analysis_type == "RUV applied") |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    pivot_wider(id_cols = c(!!sym(args_row_id)),
                names_from = c(comparison),
                names_sep = ":",
                values_from = c(log2FC, !!sym(qvalue_column), !!sym(raw_pvalue_column))) |>
    left_join(counts_table_to_use, by = join_by( !!sym(args_row_id) == !!sym(theObject@protein_id_column))) |>
    left_join(theObject@protein_id_table, by = join_by( !!sym(args_row_id) == !!sym(theObject@protein_id_column))) |>
    dplyr::arrange(across(matches(paste0("!!sym(", qvalue_column, ")")))) |>
    distinct()

  return_list$de_proteins_wide <- de_proteins_wide

  # Create long format output
  de_proteins_long <- createDeResultsLongFormat( lfc_qval_tbl = significant_rows |>
                                                   dplyr::filter(analysis_type == "RUV applied"),
                                                 norm_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
                                                 raw_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
                                                 row_id = args_row_id,
                                                 sample_id = theObject@sample_id,
                                                 group_id = group_id,
                                                 group_pattern = args_group_pattern,
                                                 design_matrix_norm = theObject@design_matrix,
                                                 design_matrix_raw = theObject@design_matrix,
                                                 protein_id_table = theObject@protein_id_table)

  return_list$de_proteins_long <- de_proteins_long

  # Static volcano plots with gene names
  static_volcano_plot_data <- de_proteins_long |>
    mutate( lqm = -log10(!!sym(qvalue_column))) |>
    dplyr::mutate(label = case_when( !!sym(qvalue_column) < de_q_val_thresh ~ "Significant",
                                     TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when( !!sym(qvalue_column) < de_q_val_thresh ~ "purple",
                                      TRUE ~ "black")) |>
    dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))

  list_of_volcano_plots <- static_volcano_plot_data %>%
    group_by( comparison) %>%
    nest() %>%
    ungroup() %>%
    mutate( title = paste( comparison)) %>%
    mutate( plot = purrr::map2( data, title, \(x,y) { plotOneVolcanoNoVerticalLines(x, y
                                                                                     , log_q_value_column = lqm
                                                                                     , log_fc_column = log2FC) } ) )

  return_list$list_of_volcano_plots <- list_of_volcano_plots

  # Additional summary statistics
  num_sig_de_molecules <- significant_rows %>%
    dplyr::mutate(status = case_when(!!sym(qvalue_column) >= de_q_val_thresh ~ "Not significant",
                                     log2FC >= 0 & !!sym(qvalue_column) < de_q_val_thresh ~ "Significant and Up",
                                     log2FC < 0 & !!sym(qvalue_column) < de_q_val_thresh ~ "Significant and Down",
                                     TRUE ~ "Not significant")) %>%
    group_by( comparison, status) %>%
    summarise(counts = n()) %>%
    ungroup()

  return_list$num_sig_de_molecules <- num_sig_de_molecules

  # Create barplots if significant results exist
  if (num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") |>
      nrow() > 0 ) {

    num_sig_de_genes_barplot_only_significant <- num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~ comparison)

    return_list$num_sig_de_genes_barplot_only_significant <- num_sig_de_genes_barplot_only_significant

    num_sig_de_genes_barplot_with_not_significant <- num_sig_de_molecules %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~ comparison)

    return_list$num_sig_de_genes_barplot_with_not_significant <- num_sig_de_genes_barplot_with_not_significant
  }

  message("--- Exiting differentialExpressionAnalysisHelper ---")
  return(return_list)

})

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Generate Interactive Volcano Plot using Glimma
#' 
#' @export
generateVolcanoPlotGlimma <- function( de_results_list
                                       , selected_contrast = NULL
                                       , uniprot_tbl = NULL
                                       , args_row_id = "uniprot_acc"
                                       , fdr_column = "fdr_qvalue"
                                       , raw_p_value_column = "raw_pvalue"
                                       , log2fc_column = "log2FC"
                                       , de_q_val_thresh = 0.05
                                       , uniprot_id_column = "Entry"
                                       , gene_names_column = "gene_names"
                                       , display_columns = c( "best_uniprot_acc" )) {

  message("--- Entering generateVolcanoPlotGlimma ---")
  message(paste("   generateVolcanoPlotGlimma Arg: selected_contrast =", selected_contrast))

  if (is.null(de_results_list) || is.null(de_results_list$de_proteins_long)) {
    message("   generateVolcanoPlotGlimma: No DE results available")
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    message("   generateVolcanoPlotGlimma: No contrast selected")
    return(NULL)
  }

  # Get the contrast-specific results
  de_proteins_long <- de_results_list$de_proteins_long
  contrasts_results <- de_results_list$contrasts_results

  # CRITICAL FIX: Extract comparison name from selected_contrast if it contains "="
  # The selected_contrast might be the full name like "GA_Elevated.minus.GA_Control=groupGA_Elevated-groupGA_Control"
  # But the comparison column in de_proteins_long contains only "GA_Elevated.minus.GA_Control"
  comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
  message(paste("   generateVolcanoPlotGlimma: Searching for comparison =", comparison_to_search))
  message(paste("   generateVolcanoPlotGlimma: Original selected_contrast =", selected_contrast))

  # Filter for selected contrast using the extracted comparison name
  contrast_data <- de_proteins_long |>
    dplyr::filter(comparison == comparison_to_search)

  if (nrow(contrast_data) == 0) {
    message(paste("   generateVolcanoPlotGlimma: No data found for contrast", comparison_to_search))
    # DEBUG: Show what comparisons are actually available
    available_comparisons <- unique(de_proteins_long$comparison)
    message(paste("   generateVolcanoPlotGlimma: Available comparisons:", paste(available_comparisons, collapse = ", ")))
    return(NULL)
  }

  message(paste("   generateVolcanoPlotGlimma: Found", nrow(contrast_data), "rows for contrast", comparison_to_search))
  message("   generateVolcanoPlotGlimma Step: Preparing volcano plot data...")

  # Prepare volcano plot table
  volcano_plot_tab <- contrast_data |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":") |> purrr::map_chr(1))

  # Add UniProt annotations if available
  if (!is.null(uniprot_tbl)) {
    volcano_plot_tab <- volcano_plot_tab |>
      left_join(uniprot_tbl, by = join_by( best_uniprot_acc == !!sym( uniprot_id_column) )) |>
      dplyr::rename( UNIPROT_GENENAME = !!sym(gene_names_column) ) |>
      mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){
        tryCatch({
          if(is.na(x) || is.null(x) || x == "") {
            ""
          } else {
            split_result <- str_split(x, " |:")[[1]]
            if(length(split_result) > 0) split_result[1] else ""
          }
        }, error = function(e) "")
      })) |>
      dplyr::mutate(gene_name = UNIPROT_GENENAME )
  } else {
    volcano_plot_tab <- volcano_plot_tab |>
      dplyr::mutate(gene_name = best_uniprot_acc)
  }

  # Add visualization columns
  volcano_plot_tab <- volcano_plot_tab |>
    mutate( lqm = -log10(!!sym(fdr_column))) |>
    dplyr::mutate(label = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "Not sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC < 1",
                                    TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "orange",
                                     abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                     abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "blue",
                                     TRUE ~ "black")) |>
    dplyr::mutate(analysis_type = comparison) |>
    dplyr::select( best_uniprot_acc, lqm, !!sym(fdr_column), !!sym(raw_p_value_column), !!sym(log2fc_column), comparison, label, colour, gene_name
                   , any_of( display_columns )) |>
    dplyr::mutate( my_alpha = case_when ( gene_name != "" ~ 1
                                          , TRUE ~ 0.5))

  message("   generateVolcanoPlotGlimma Step: Generating interactive plot widget...")

  # Find the coefficient index for this contrast
  # CRITICAL FIX: Use grep for reliable pattern matching instead of stringr
  coef_names <- colnames(contrasts_results$fit.eb$coefficients)
  message(paste("   generateVolcanoPlotGlimma: Available coefficient names:", paste(coef_names, collapse = ", ")))
  message(paste("   generateVolcanoPlotGlimma: Looking for pattern:", paste0("^", comparison_to_search, "=")))
  
  # Use grep to find coefficient that starts with friendly name followed by "="
  coef_index <- grep(paste0("^", comparison_to_search, "="), coef_names)
  message(paste("   generateVolcanoPlotGlimma: Pattern match found indices:", paste(coef_index, collapse = ", ")))
  
  # Fallback: try exact match with the friendly name (for backwards compatibility)
  if (length(coef_index) == 0) {
    message("   generateVolcanoPlotGlimma: Pattern match failed, trying exact match with friendly name")
    coef_index <- which(coef_names == comparison_to_search)
  }
  
  # Last resort: try exact match with selected_contrast
  if (length(coef_index) == 0) {
    message("   generateVolcanoPlotGlimma: Exact friendly name match failed, trying selected_contrast")
    coef_index <- which(coef_names == selected_contrast)
  }

  if (length(coef_index) == 0) {
    message(paste("   generateVolcanoPlotGlimma: FINAL FAILURE - No coefficient found for:", selected_contrast))
    message(paste("   generateVolcanoPlotGlimma: Also tried:", comparison_to_search))
    message("   generateVolcanoPlotGlimma: Available coefficient patterns:")
    for (i in seq_along(coef_names)) {
      message(sprintf("     [%d] %s", i, coef_names[i]))
    }
    message(paste("   generateVolcanoPlotGlimma: Pattern attempted:", paste0("^", comparison_to_search, "=")))
    return(NULL)
  }

  message(paste("   generateVolcanoPlotGlimma: Using coefficient index", coef_index, "for contrast", coef_names[coef_index]))

  # Prepare counts matrix and groups
  counts_mat <- de_results_list$theObject@protein_quant_table |>
    column_to_rownames(de_results_list$theObject@protein_id_column) |>
    as.matrix()

  this_design_matrix <- de_results_list$theObject@design_matrix
  rownames( this_design_matrix ) <- this_design_matrix[, de_results_list$theObject@sample_id]
  this_groups <- this_design_matrix[colnames( counts_mat), "group"]

  # Generate the Glimma widget
  glimma_widget <- getGlimmaVolcanoProteomicsWidget( contrasts_results$fit.eb
                                                    , coef = coef_index
                                                    , volcano_plot_tab = volcano_plot_tab
                                                    , uniprot_column = best_uniprot_acc
                                                    , gene_name_column = gene_name
                                                    , display_columns = display_columns
                                                    , counts_tbl = counts_mat
                                                    , groups = this_groups )

  message("--- Exiting generateVolcanoPlotGlimma ---")
  return(glimma_widget)
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Generate DE Results Heatmap with Advanced Clustering
#' 
#' @export
generateDEHeatmap <- function( de_results_list
                              , selected_contrast = NULL
                              , top_n_genes = 50
                              , clustering_method = "ward.D2"
                              , distance_method = "euclidean"
                              , cluster_rows = TRUE
                              , cluster_cols = TRUE
                              , scale_data = "row"
                              , color_scheme = "RdBu"
                              , show_gene_names = FALSE
                              , de_q_val_thresh = 0.05
                              , qvalue_column = "fdr_qvalue"
                              , log2fc_column = "log2FC") {

  message("--- Entering generateDEHeatmap ---")
  message(sprintf("   generateDEHeatmap Arg: selected_contrast = %s", selected_contrast))
  message(sprintf("   generateDEHeatmap Arg: top_n_genes = %d", top_n_genes))

  if (is.null(de_results_list) || is.null(de_results_list$de_proteins_long)) {
    message("   generateDEHeatmap: No DE results available")
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    message("   generateDEHeatmap: No contrast selected")
    return(NULL)
  }

  message("   generateDEHeatmap Step: Filtering data for selected contrast...")

  # Get the contrast-specific results
  de_proteins_long <- de_results_list$de_proteins_long
  
  # Filter for selected contrast and significant results
  contrast_data <- de_proteins_long |>
    dplyr::filter(comparison == selected_contrast) |>
    dplyr::filter(!!sym(qvalue_column) < de_q_val_thresh)

  if (nrow(contrast_data) == 0) {
    message(sprintf("   generateDEHeatmap: No significant genes found for contrast %s", selected_contrast))
    return(NULL)
  }

  message(sprintf("   generateDEHeatmap Step: Found %d significant genes", nrow(contrast_data)))

  # Order by absolute fold change and take top N
  top_genes_data <- contrast_data |>
    dplyr::arrange(desc(abs(!!sym(log2fc_column)))) |>
    head(top_n_genes)

  message(sprintf("   generateDEHeatmap Step: Selected top %d genes for heatmap", nrow(top_genes_data)))

  # Extract expression matrix for these genes
  protein_ids <- top_genes_data |> dplyr::pull(!!sym(names(top_genes_data)[1])) # First column should be protein ID

  # Get expression data from the S4 object
  expr_matrix <- de_results_list$theObject@protein_quant_table |>
    dplyr::filter(!!sym(de_results_list$theObject@protein_id_column) %in% protein_ids) |>
    column_to_rownames(de_results_list$theObject@protein_id_column) |>
    as.matrix()

  if (nrow(expr_matrix) == 0) {
    message("   generateDEHeatmap: No expression data found for selected genes")
    return(NULL)
  }

  message(sprintf("   generateDEHeatmap Step: Expression matrix dims = %d rows, %d cols", nrow(expr_matrix), ncol(expr_matrix)))

  # Apply scaling
  if (scale_data == "row") {
    message("   generateDEHeatmap Step: Applying row scaling...")
    expr_matrix <- t(scale(t(expr_matrix)))
  } else if (scale_data == "column") {
    message("   generateDEHeatmap Step: Applying column scaling...")
    expr_matrix <- scale(expr_matrix)
  } else if (scale_data == "both") {
    message("   generateDEHeatmap Step: Applying both row and column scaling...")
    expr_matrix <- scale(t(scale(t(expr_matrix))))
  } else {
    message("   generateDEHeatmap Step: No scaling applied")
  }

  message("   generateDEHeatmap Step: Setting up clustering parameters...")

  # Set up clustering
  row_clust <- NULL
  col_clust <- NULL

  if (cluster_rows || cluster_cols) {
    # Calculate distance matrices
    if (distance_method %in% c("pearson", "spearman")) {
      if (cluster_rows) {
        if (distance_method == "pearson") {
          row_dist <- as.dist(1 - cor(t(expr_matrix), method = "pearson", use = "complete.obs"))
        } else {
          row_dist <- as.dist(1 - cor(t(expr_matrix), method = "spearman", use = "complete.obs"))
        }
        row_clust <- hclust(row_dist, method = clustering_method)
      }
      
      if (cluster_cols) {
        if (distance_method == "pearson") {
          col_dist <- as.dist(1 - cor(expr_matrix, method = "pearson", use = "complete.obs"))
        } else {
          col_dist <- as.dist(1 - cor(expr_matrix, method = "spearman", use = "complete.obs"))
        }
        col_clust <- hclust(col_dist, method = clustering_method)
      }
    } else {
      # Standard distance metrics
      if (cluster_rows) {
        row_dist <- dist(expr_matrix, method = distance_method)
        row_clust <- hclust(row_dist, method = clustering_method)
      }
      
      if (cluster_cols) {
        col_dist <- dist(t(expr_matrix), method = distance_method)
        col_clust <- hclust(col_dist, method = clustering_method)
      }
    }
  }

  message("   generateDEHeatmap Step: Setting up color scheme...")

  # Choose color palette
  colors <- switch(color_scheme,
    "RdBu" = colorRampPalette(c("red", "white", "blue"))(100),
    "RdYlBu" = colorRampPalette(c("red", "yellow", "blue"))(100),
    "coolwarm" = colorRampPalette(c("blue", "white", "red"))(100),
    "viridis" = viridis::viridis(100),
    "plasma" = viridis::plasma(100),
    "inferno" = viridis::inferno(100),
    colorRampPalette(c("red", "white", "blue"))(100)
  )

  message("   generateDEHeatmap Step: Generating heatmap...")

  # Generate the heatmap
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    message("   generateDEHeatmap: Using ComplexHeatmap for advanced visualization")
    
    heatmap_plot <- ComplexHeatmap::Heatmap(
      expr_matrix,
      name = "Expression",
      col = colors,
      cluster_rows = if (cluster_rows) row_clust else FALSE,
      cluster_columns = if (cluster_cols) col_clust else FALSE,
      show_row_names = show_gene_names,
      show_column_names = TRUE,
      show_row_dend = cluster_rows,
      show_column_dend = cluster_cols,
      column_title = paste("Heatmap - Top", nrow(expr_matrix), "genes\nContrast:", selected_contrast),
      heatmap_legend_param = list(title = "Scaled\nExpression")
    )
  } else {
    message("   generateDEHeatmap: Using base R heatmap")
    
    # Create a plot function that returns the heatmap
    heatmap_plot <- function() {
      heatmap(
        expr_matrix,
        main = paste("Heatmap - Top", nrow(expr_matrix), "genes\nContrast:", selected_contrast),
        col = colors,
        Rowv = if (cluster_rows) as.dendrogram(row_clust) else NA,
        Colv = if (cluster_cols) as.dendrogram(col_clust) else NA,
        labRow = if (show_gene_names) rownames(expr_matrix) else rep("", nrow(expr_matrix)),
        cexRow = if (show_gene_names) 0.8 else 0.1
      )
    }
  }

  message("--- Exiting generateDEHeatmap ---")
  return(heatmap_plot)
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
setMethod(f = "outputDeResultsAllContrasts",
          signature = "ProteinQuantitativeData",
          definition = function(theObject,
                               de_results_list_all_contrasts = NULL,
                               uniprot_tbl = NULL,
                               de_output_dir = NULL,
                               publication_graphs_dir = NULL,
                               file_prefix = "de_proteins",
                               args_row_id = NULL,
                               gene_names_column = "gene_names",
                               uniprot_id_column = "Entry") {

  message("--- Entering outputDeResultsAllContrasts ---")
  message(sprintf("   outputDeResultsAllContrasts: de_output_dir = %s", de_output_dir))
  message(sprintf("   outputDeResultsAllContrasts: file_prefix = %s", file_prefix))
  
  # Extract parameters from S4 object with fallbacks
  uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", uniprot_tbl)
  de_output_dir <- checkParamsObjectFunctionSimplify(theObject, "de_output_dir", de_output_dir)
  publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", publication_graphs_dir)
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", args_row_id)
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", gene_names_column)
  uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", uniprot_id_column)
  
  # Ensure output directory exists (CRITICAL FIX!)
  if (!dir.exists(de_output_dir)) {
    dir.create(de_output_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("   outputDeResultsAllContrasts: Created output directory: %s", de_output_dir))
  }
  
  message(sprintf("   outputDeResultsAllContrasts: Processing %d contrasts", length(de_results_list_all_contrasts)))
  
  # Write results for each contrast with UNIQUE filenames
  for (contrast_name in names(de_results_list_all_contrasts)) {
    message(sprintf("   outputDeResultsAllContrasts: Writing files for contrast: %s", contrast_name))
    
    contrast_result <- de_results_list_all_contrasts[[contrast_name]]
    
    if (!is.null(contrast_result$de_proteins_long)) {
      # Clean contrast name for safe filename
      safe_contrast_name <- gsub("[^A-Za-z0-9_-]", "_", contrast_name)
      
             # Create annotated version (FIXED: proper conditional logic)
       de_proteins_long_annot <- contrast_result$de_proteins_long |>
         mutate(uniprot_acc_cleaned = str_split(!!sym(args_row_id), "-") |>
                  purrr::map_chr(1))
       
       # Add UniProt annotations if available
       if (!is.null(uniprot_tbl) && nrow(uniprot_tbl) > 0) {
         de_proteins_long_annot <- de_proteins_long_annot |>
           left_join(uniprot_tbl, by = join_by(uniprot_acc_cleaned == !!sym(uniprot_id_column))) |>
           dplyr::select(-uniprot_acc_cleaned)
       } else {
         de_proteins_long_annot <- de_proteins_long_annot |>
           dplyr::select(-uniprot_acc_cleaned)
       }
       
       # Add gene_name column with proper conditional logic
       if ("gene_names" %in% names(de_proteins_long_annot)) {
         de_proteins_long_annot <- de_proteins_long_annot |>
           mutate(gene_name = purrr::map_chr(gene_names, \(x){
             tryCatch({
               if(is.na(x) || is.null(x) || x == "") {
                 ""
               } else {
                 split_result <- str_split(x, " |:")[[1]]
                 if(length(split_result) > 0) split_result[1] else ""
               }
             }, error = function(e) "")
           }))
       } else if (gene_names_column %in% names(de_proteins_long_annot)) {
         de_proteins_long_annot <- de_proteins_long_annot |>
           mutate(gene_name = purrr::map_chr(.data[[gene_names_column]], \(x){
             tryCatch({
               if(is.na(x) || is.null(x) || x == "") {
                 ""
               } else {
                 split_result <- str_split(x, " |:")[[1]]
                 if(length(split_result) > 0) split_result[1] else ""
               }
             }, error = function(e) "")
           }))
       } else {
         de_proteins_long_annot <- de_proteins_long_annot |>
           mutate(gene_name = "")
       }
       
       # Relocate gene_name column
       de_proteins_long_annot <- de_proteins_long_annot |>
         relocate(gene_name, .after = !!sym(args_row_id))
      
      # CRITICAL FIX: Use contrast-specific filename!
      contrast_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.tsv")
      output_path <- file.path(de_output_dir, contrast_filename)
      
      message(sprintf("   outputDeResultsAllContrasts: Writing %s", output_path))
      
      # Write the file
      vroom::vroom_write(de_proteins_long_annot, output_path)
      
      # Verify file was written
      if (file.exists(output_path)) {
        file_size <- file.size(output_path)
        message(sprintf("   outputDeResultsAllContrasts: SUCCESS - File written: %s (%d bytes)", 
                       contrast_filename, file_size))
      } else {
        message(sprintf("   outputDeResultsAllContrasts: ERROR - File NOT written: %s", output_path))
      }
      
      # Also write Excel version
      excel_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.xlsx")
      excel_path <- file.path(de_output_dir, excel_filename)
      writexl::write_xlsx(de_proteins_long_annot, excel_path)
      
      message(sprintf("   outputDeResultsAllContrasts: Also wrote Excel: %s", excel_filename))
      
      # âœ… NEW: Generate volcano plot for this contrast
      message(sprintf("   outputDeResultsAllContrasts: Generating volcano plot for contrast: %s", contrast_name))
      
      # Extract parameters - try direct access first, then use checkParams
      de_q_val_thresh <- if (!is.null(theObject@args$outputDeResultsAllContrasts$de_q_val_thresh)) {
        theObject@args$outputDeResultsAllContrasts$de_q_val_thresh
      } else {
        checkParamsObjectFunctionSimplify(theObject, "de_q_val_thresh", 0.05)
      }
      
      fdr_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$fdr_column)) {
        theObject@args$outputDeResultsAllContrasts$fdr_column
      } else {
        checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
      }
      
      log2fc_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$log2fc_column)) {
        theObject@args$outputDeResultsAllContrasts$log2fc_column
      } else {
        checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
      }
      
      # Prepare data for volcano plot (same logic as in differentialExpressionAnalysisHelper)
      volcano_data <- contrast_result$de_proteins_long |>
        mutate(lqm = -log10(!!sym(fdr_column))) |>
        dplyr::mutate(label = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "Significant",
                                       TRUE ~ "Not sig.")) |>
        dplyr::mutate(colour = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                        TRUE ~ "black")) |>
        dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))
      
      # Generate the volcano plot
      volcano_plot <- plotOneVolcanoNoVerticalLines(
        volcano_data, 
        contrast_name,
        log_q_value_column = lqm,
        log_fc_column = !!sym(log2fc_column)
      )
      
      # Create Volcano_Plots directory if it doesn't exist
      volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
      if (!dir.exists(volcano_dir)) {
        dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
        message(sprintf("   outputDeResultsAllContrasts: Created volcano plots directory: %s", volcano_dir))
      }
      
      # Save volcano plot with contrast-specific filename
      volcano_png_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".png"))
      volcano_pdf_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".pdf"))
      
      # Save as PNG
      tryCatch({
        ggplot2::ggsave(volcano_png_path, volcano_plot, width = 7, height = 7, dpi = 300)
        message(sprintf("   outputDeResultsAllContrasts: Saved volcano plot PNG: %s", basename(volcano_png_path)))
      }, error = function(e) {
        message(sprintf("   outputDeResultsAllContrasts: ERROR saving PNG: %s", e$message))
      })
      
      # Save as PDF
      tryCatch({
        ggplot2::ggsave(volcano_pdf_path, volcano_plot, width = 7, height = 7)
        message(sprintf("   outputDeResultsAllContrasts: Saved volcano plot PDF: %s", basename(volcano_pdf_path)))
      }, error = function(e) {
        message(sprintf("   outputDeResultsAllContrasts: ERROR saving PDF: %s", e$message))
      })
    }
  }
  
  # âœ… NEW: Generate combined multi-page PDF with all volcano plots
  message("   outputDeResultsAllContrasts: Creating combined volcano plots PDF...")
  
  volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
  combined_pdf_path <- file.path(volcano_dir, "all_volcano_plots_combined.pdf")
  
  tryCatch({
    # Re-extract parameters for combined PDF generation (in case they were modified)
    de_q_val_thresh <- if (!is.null(theObject@args$outputDeResultsAllContrasts$de_q_val_thresh)) {
      theObject@args$outputDeResultsAllContrasts$de_q_val_thresh
    } else {
      0.05
    }
    
    fdr_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$fdr_column)) {
      theObject@args$outputDeResultsAllContrasts$fdr_column
    } else {
      "fdr_qvalue"
    }
    
    log2fc_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$log2fc_column)) {
      theObject@args$outputDeResultsAllContrasts$log2fc_column
    } else {
      "log2FC"
    }
    
    # Collect all volcano plots
    all_volcano_plots <- list()
    
    for (contrast_name in names(de_results_list_all_contrasts)) {
      contrast_result <- de_results_list_all_contrasts[[contrast_name]]
      
      if (!is.null(contrast_result$de_proteins_long)) {
        # Recreate the volcano plot
        volcano_data <- contrast_result$de_proteins_long |>
          mutate(lqm = -log10(!!sym(fdr_column))) |>
          dplyr::mutate(label = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "Significant",
                                         TRUE ~ "Not sig.")) |>
          dplyr::mutate(colour = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                          TRUE ~ "black")) |>
          dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))
        
        volcano_plot <- plotOneVolcanoNoVerticalLines(
          volcano_data, 
          contrast_name,
          log_q_value_column = lqm,
          log_fc_column = !!sym(log2fc_column)
        )
        
        all_volcano_plots[[contrast_name]] <- volcano_plot
      }
    }
    
    # Generate multi-page PDF
    if (length(all_volcano_plots) > 0) {
      pdf(file = combined_pdf_path, width = 7, height = 7, onefile = TRUE)
      purrr::walk(all_volcano_plots, print)
      invisible(dev.off())
      message(sprintf("   outputDeResultsAllContrasts: Created combined PDF with %d volcano plots: %s", 
                     length(all_volcano_plots), basename(combined_pdf_path)))
    }
    
  }, error = function(e) {
    message(sprintf("   outputDeResultsAllContrasts: ERROR creating combined PDF: %s", e$message))
  })
  
  message("--- Exiting outputDeResultsAllContrasts ---")
  return(TRUE)
}) 

# ==========================================
# Missing S4 Methods recovered from GUI branch
# ==========================================

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Get Low Coefficient of Variation Proteins
#' 
#' Identifies proteins with the lowest coefficient of variation for use as 
#' negative controls in RUV normalization.
#'
#' @param theObject A ProteinQuantitativeData object
#' @param percentage_as_neg_ctrl Percentage of proteins to use as negative controls (default: 10)
#' @param num_neg_ctrl Number of negative control proteins (overrides percentage if set)
#' @return A logical vector indicating which proteins are selected as negative controls
#' @export
setMethod( f = "getLowCoefficientOfVariationProteins"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject
                                  , percentage_as_neg_ctrl = NULL
                                  , num_neg_ctrl = NULL) {

             percentage_as_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject, "percentage_as_neg_ctrl", 10)
             num_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject
                                                                , "num_neg_ctrl"
                                                                , round(nrow( theObject@protein_quant_table) * percentage_as_neg_ctrl / 100, 0))

             theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
             theObject <- updateParamInObject(theObject, "num_neg_ctrl")

  list_of_control_genes <- theObject@protein_quant_table |>
    column_to_rownames(theObject@protein_id_column) |>
    t() |>
    as.data.frame() |>
    summarise( across(everything(), ~sd(.)/mean(.))) |>
    t() |>
    as.data.frame() |>
    dplyr::rename( coefficient_of_variation = "V1") |>
    tibble::rownames_to_column(theObject@protein_id_column) |>
    arrange( coefficient_of_variation) |>
    head(num_neg_ctrl)

  control_gene_index_helper <- theObject@protein_quant_table |>
    dplyr::select(theObject@protein_id_column) |>
    mutate( index = row_number()) |>
    left_join( list_of_control_genes, by = theObject@protein_id_column)  |>
    mutate( is_selected = case_when( is.na(coefficient_of_variation) ~ FALSE
                                     , TRUE ~  TRUE ) ) |>
    arrange( index) |>
    dplyr::select( !!sym(theObject@protein_id_column), is_selected) |>
    column_to_rownames(theObject@protein_id_column) |>
    t()

  control_gene_index <- control_gene_index_helper[1,]

  control_gene_index

})


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Average Technical Replicates
#' 
#' Averages protein intensities across technical replicates, collapsing the
#' design matrix accordingly.
#'
#' @param theObject A ProteinQuantitativeData object
#' @param design_matrix_columns Additional columns to keep in the design matrix
#' @return The modified ProteinQuantitativeData object with averaged technical replicates
#' @rdname averageTechReps
#' @export
setMethod( f = "averageTechReps"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, design_matrix_columns=c()  ) {

             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             theObject@protein_quant_table <- protein_quant_table |>
               pivot_longer( cols = !matches( protein_id_column)
                             , names_to = sample_id
                             , values_to = "Log2.Protein.Imputed") |>
               left_join( design_matrix
                          , by = join_by( !!sym(sample_id) == !!sym(sample_id))) |>
               group_by( !!sym(protein_id_column), !!sym(replicate_group_column) )  |>
               summarise( Log2.Protein.Imputed = mean( Log2.Protein.Imputed, na.rm = TRUE)) |>
               ungroup() |>
               pivot_wider( names_from = !!sym(replicate_group_column)
                            , values_from = Log2.Protein.Imputed)

              theObject@sample_id <- theObject@technical_replicate_id

              theObject@design_matrix <- design_matrix |>
                dplyr::select(-!!sym( sample_id)) |>
                dplyr::select(all_of( unique( c( replicate_group_column,  group_id,  design_matrix_columns) ))) |>
                distinct()

              theObject@sample_id <- replicate_group_column
              theObject@technical_replicate_id <- NA_character_

              theObject <- cleanDesignMatrix(theObject)

              theObject

           }) 

           ##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeProteinsWithOnlyOneReplicate"
           , def=function( theObject, core_utilisation = NULL, grouping_variable = NULL) {
             standardGeneric("removeProteinsWithOnlyOneReplicate")
           }
           , signature=c("theObject", "core_utilisation", "grouping_variable"))

#'@export
setMethod(f="removeProteinsWithOnlyOneReplicate"
          , definition=function( theObject, core_utilisation = NULL, grouping_variable = NULL) {
            protein_quant_table <- theObject@protein_quant_table
            samples_id_tbl <- theObject@design_matrix
            sample_id_tbl_sample_id_column <- theObject@sample_id
            # replicate_group_column <- theObject@technical_replicate_id
            protein_id_column <- theObject@protein_id_column

            input_table_sample_id_column <- theObject@sample_id
            quantity_column <- "log_values"

            grouping_variable <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                                              , "grouping_variable"
                                                                              , NULL)

            core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                                   , "core_utilisation"
                                                                   , NA)

            theObject <- updateParamInObject(theObject, "grouping_variable")
            theObject <- updateParamInObject(theObject, "core_utilisation")

            data_long_cln <- protein_quant_table  |>
              pivot_longer( cols=!matches(protein_id_column)
                            , names_to = input_table_sample_id_column
                            , values_to = quantity_column)

            protein_quant_table <- removeProteinsWithOnlyOneReplicateHelper( input_table = data_long_cln
                                                                , samples_id_tbl = samples_id_tbl
                                                                , input_table_sample_id_column = !!sym( input_table_sample_id_column )
                                                                , sample_id_tbl_sample_id_column = !!sym( sample_id_tbl_sample_id_column)
                                                                , replicate_group_column = !!sym(grouping_variable)
                                                                , protein_id_column = !!sym( protein_id_column)
                                                                , quantity_column = !!sym( quantity_column)
                                                                , core_utilisation = core_utilisation)


            theObject@protein_quant_table <- protein_quant_table |>
              pivot_wider( id_cols = !!sym( protein_id_column)
                           , names_from = !!sym( input_table_sample_id_column)
                           , values_from = !!sym( quantity_column) )

            theObject <- cleanDesignMatrix(theObject)

            updated_object <- theObject

            return(updated_object)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric(name="removeRowsWithMissingValuesPercent"
           , def=function( theObject
                           , ruv_grouping_variable = NULL
                           , groupwise_percentage_cutoff = NULL
                           , max_groups_percentage_cutoff = NULL
                           , proteins_intensity_cutoff_percentile = NULL ) {
             standardGeneric("removeRowsWithMissingValuesPercent")
           }
           , signature=c("theObject"
                         , "ruv_grouping_variable"
                         , "groupwise_percentage_cutoff"
                         , "max_groups_percentage_cutoff"
                         , "proteins_intensity_cutoff_percentile" ))

#'@export
setMethod( f = "removeRowsWithMissingValuesPercent"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject
                                  , ruv_grouping_variable = NULL
                                  , groupwise_percentage_cutoff = NULL
                                  , max_groups_percentage_cutoff = NULL
                                  , proteins_intensity_cutoff_percentile = NULL) {

             print("--- Entering removeRowsWithMissingValuesPercent S4 Method ---")
             
             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             print("   removeRowsWithMissingValuesPercent: Extracting input arguments...")
             print(sprintf("      Arg: protein_id_column = %s", protein_id_column))
             print(sprintf("      Arg: sample_id = %s", sample_id))
             print(sprintf("      Arg: group_id = %s", group_id))
             print(sprintf("      Data State (protein_quant_table): Dims = %d rows, %d cols", nrow(protein_quant_table), ncol(protein_quant_table)))
             print(sprintf("      Data State (design_matrix): Dims = %d rows, %d cols", nrow(design_matrix), ncol(design_matrix)))
             print(head(design_matrix))

             # print(groupwise_percentage_cutoff)
             # print(min_protein_intensity_threshold )

             print("   removeRowsWithMissingValuesPercent: Resolving parameters with checkParamsObjectFunctionSimplify...")
             ruv_grouping_variable <- checkParamsObjectFunctionSimplify(theObject
                                                                        , "ruv_grouping_variable"
                                                                        , NULL)
             print(sprintf("      Resolved ruv_grouping_variable = %s", ifelse(is.null(ruv_grouping_variable), "NULL", ruv_grouping_variable)))
             
             groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify(theObject
                                                                              , "groupwise_percentage_cutoff"
                                                                              , 50)
             print(sprintf("      Resolved groupwise_percentage_cutoff = %g", groupwise_percentage_cutoff))
             
             max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify(theObject
                                                                               , "max_groups_percentage_cutoff"
                                                                               , 50)
             print(sprintf("      Resolved max_groups_percentage_cutoff = %g", max_groups_percentage_cutoff))
             
             proteins_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify(theObject
                                                                                   , "proteins_intensity_cutoff_percentile"
                                                                                   , 1)
             print(sprintf("      Resolved proteins_intensity_cutoff_percentile = %g", proteins_intensity_cutoff_percentile))

             print("   removeRowsWithMissingValuesPercent: Updating parameters in S4 object...")
             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
             theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "proteins_intensity_cutoff_percentile")

             print("   removeRowsWithMissingValuesPercent: About to call helper function...")
             print(sprintf("      Helper Args: cols = %s", protein_id_column))
             print(sprintf("      Helper Args: sample_id = %s", sample_id))
             print(sprintf("      Helper Args: row_id = %s", protein_id_column))
             print(sprintf("      Helper Args: grouping_variable = %s", ruv_grouping_variable))
             print(sprintf("      Helper Args: groupwise_percentage_cutoff = %g", groupwise_percentage_cutoff))
             print(sprintf("      Helper Args: max_groups_percentage_cutoff = %g", max_groups_percentage_cutoff))
             print(sprintf("      Helper Args: proteins_intensity_cutoff_percentile = %g", proteins_intensity_cutoff_percentile))

             theObject@protein_quant_table <- removeRowsWithMissingValuesPercentHelper( protein_quant_table
                                                                           , cols= protein_id_column
                                                                           , design_matrix = design_matrix
                                                                           , sample_id = !!sym(sample_id)
                                                                           , row_id = !!sym(protein_id_column)
                                                                           , grouping_variable = !!sym(ruv_grouping_variable)
                                                                           , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                                                           , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                                                           , proteins_intensity_cutoff_percentile = proteins_intensity_cutoff_percentile
                                                                           , temporary_abundance_column = "Log_Abundance")

             print(sprintf("   removeRowsWithMissingValuesPercent: Helper function returned. New dims = %d rows, %d cols", 
                           nrow(theObject@protein_quant_table), ncol(theObject@protein_quant_table)))

             print("   removeRowsWithMissingValuesPercent: Cleaning design matrix...")
             theObject <- cleanDesignMatrix(theObject)

             print("--- Exiting removeRowsWithMissingValuesPercent S4 Method ---")
             return(theObject)

           })

#'@title Choose Best Protein Accession
#'@export
setGeneric(name="chooseBestProteinAccession"
           , def=function(theObject, delim=NULL, seqinr_obj=NULL, seqinr_accession_column=NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
             standardGeneric("chooseBestProteinAccession")
           }
           , signature=c("theObject", "delim", "seqinr_obj", "seqinr_accession_column"))

#'@rdname chooseBestProteinAccession
#'@export
#'@param theObject The object of class ProteinQuantitativeData
#'@param delim The delimiter used to split the protein accessions
#'@param seqinr_obj The object of class Seqinr::seqinr
#'@param seqinr_accession_column The column in the seqinr object that contains the protein accessions
#'@param replace_zero_with_na Replace zero values with NA
#'@param aggregation_method Method to aggregate protein values: "sum", "mean", or "median" (default: "sum")
setMethod(f = "chooseBestProteinAccession"
          , signature="ProteinQuantitativeData"
          , definition=function(theObject, delim=NULL, seqinr_obj=NULL
                              , seqinr_accession_column=NULL
                              , replace_zero_with_na = NULL
                              , aggregation_method = NULL) {
            
            cat("\n\n╔═══════════════════════════════════════════════════════════════════════════╗\n")
            cat("║           ENTERING chooseBestProteinAccession (DEBUG66)               ║\n")
            cat("╚═══════════════════════════════════════════════════════════════════════════╝\n\n")

            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            
            cat(">>> STEP 1: PARAMETER INPUTS (before checkParamsObjectFunctionSimplify) <<<\n")
            cat(sprintf("   Input parameter 'delim' (user provided): %s\n", ifelse(is.null(delim), "NULL", delim)))
            cat(sprintf("   Input parameter 'seqinr_obj' is NULL: %s\n", is.null(seqinr_obj)))
            cat(sprintf("   Input parameter 'seqinr_accession_column': %s\n", ifelse(is.null(seqinr_accession_column), "NULL", seqinr_accession_column)))
            cat(sprintf("   Input parameter 'replace_zero_with_na': %s\n", ifelse(is.null(replace_zero_with_na), "NULL", replace_zero_with_na)))
            cat(sprintf("   Input parameter 'aggregation_method': %s\n", ifelse(is.null(aggregation_method), "NULL", aggregation_method)))
            cat(sprintf("   Protein ID column from S4: %s\n", protein_id_column))
            cat(sprintf("   Protein table dimensions: %d rows x %d cols\n", nrow(protein_quant_table), ncol(protein_quant_table)))
            cat(sprintf("   First 5 Protein IDs: %s\n", paste(head(protein_quant_table[[protein_id_column]], 5), collapse = ", ")))
            cat("\n")

            delim <- checkParamsObjectFunctionSimplify(theObject, "delim",  default_value =  " |;|:|\\|")
            seqinr_obj <- checkParamsObjectFunctionSimplify(theObject, "seqinr_obj",  default_value = NULL)
            seqinr_accession_column <- checkParamsObjectFunctionSimplify(theObject
                                                                       , "seqinr_accession_column"
                                                                       , default_value = NULL)
            replace_zero_with_na <- checkParamsObjectFunctionSimplify(theObject
                                                                    , "replace_zero_with_na"
                                                                    , default_value = FALSE)
            aggregation_method <- checkParamsObjectFunctionSimplify(theObject
                                                                  , "aggregation_method"
                                                                  , default_value = "sum")
            
            cat(">>> STEP 2: PARAMETERS AFTER checkParamsObjectFunctionSimplify <<<\n")
            cat(sprintf("   Resolved 'delim': '%s'\n", delim))
            cat(sprintf("   Resolved 'seqinr_obj' is NULL: %s\n", is.null(seqinr_obj)))
            cat(sprintf("   Resolved 'seqinr_accession_column': %s\n", ifelse(is.null(seqinr_accession_column), "NULL", seqinr_accession_column)))
            cat(sprintf("   Resolved 'replace_zero_with_na': %s\n", replace_zero_with_na))
            cat(sprintf("   Resolved 'aggregation_method': %s\n", aggregation_method))
            cat("\n")

            if (!aggregation_method %in% c("sum", "mean", "median")) {
              stop("aggregation_method must be one of: 'sum', 'mean', 'median'")
            }

            theObject <- updateParamInObject(theObject, "delim")
            theObject <- updateParamInObject(theObject, "seqinr_obj")
            theObject <- updateParamInObject(theObject, "seqinr_accession_column")
            theObject <- updateParamInObject(theObject, "replace_zero_with_na")
            theObject <- updateParamInObject(theObject, "aggregation_method")

            evidence_tbl_cleaned <- protein_quant_table |>
              distinct() |>
              mutate(row_id = row_number() -1)
            
            cat(">>> STEP 3: CALLING chooseBestProteinAccessionHelper <<<\n")
            cat(sprintf("   Passing 'delim' to helper: '%s'\n", delim))
            cat(sprintf("   First 5 Protein IDs to process: %s\n", paste(head(evidence_tbl_cleaned[[protein_id_column]], 5), collapse = ", ")))
            cat(sprintf("   Number of rows to process: %d\n", nrow(evidence_tbl_cleaned)))
            cat("\n")

            accession_gene_name_tbl <- chooseBestProteinAccessionHelper(input_tbl = evidence_tbl_cleaned,
                                                                      acc_detail_tab = seqinr_obj,
                                                                      accessions_column = !!sym(protein_id_column),
                                                                      row_id_column = seqinr_accession_column,
                                                                      group_id = row_id,
                                                                      delim = delim)
            
            cat(">>> STEP 4: RETURNED FROM chooseBestProteinAccessionHelper <<<\n")
            cat(sprintf("   Number of rows returned: %d\n", nrow(accession_gene_name_tbl)))
            if (nrow(accession_gene_name_tbl) > 0) {
              cat(sprintf("   First 5 accessions: %s\n", paste(head(accession_gene_name_tbl[[seqinr_accession_column]], 5), collapse = ", ")))
            }
            cat("\n")

            protein_log2_quant_cln <- evidence_tbl_cleaned |>
              left_join(accession_gene_name_tbl |>
                         dplyr::distinct(row_id, !!sym(as.character(seqinr_accession_column)))
                       , by = join_by(row_id)) |>
              mutate(!!sym(theObject@protein_id_column) := !!sym(as.character(seqinr_accession_column))) |>
              dplyr::select(-row_id, -!!sym(as.character(seqinr_accession_column)))

            protein_id_table <- evidence_tbl_cleaned |>
              left_join(accession_gene_name_tbl |>
                         dplyr::distinct(row_id, !!sym(as.character(seqinr_accession_column)))
                       , by = join_by(row_id)) |>
              distinct(uniprot_acc, !!sym(protein_id_column)) |>
              mutate(!!sym(paste0(protein_id_column, "_list")) := !!sym(protein_id_column)) |>
              mutate(!!sym(protein_id_column) := !!sym("uniprot_acc")) |>
              distinct(!!sym(protein_id_column), !!sym(paste0(protein_id_column, "_list"))) |>
              group_by(!!sym(protein_id_column)) |>
              summarise(!!sym(paste0(protein_id_column, "_list")) := paste(!!sym(paste0(protein_id_column, "_list")), collapse = ";")) |>
              ungroup() |>
              mutate(!!sym(paste0(protein_id_column, "_list")) := purrr::map_chr(!!sym(paste0(protein_id_column, "_list"))
                                                                                , \(x){ paste(unique(sort(str_split(x, ";")[[1]])), collapse=";") }))

            summed_data <- protein_log2_quant_cln |>
              mutate(!!sym(protein_id_column) := purrr::map_chr(!!sym(protein_id_column), \(x){ str_split(x, delim)[[1]][1] })) |>
              pivot_longer(
                cols = !matches(protein_id_column),
                names_to = "sample_id",
                values_to = "temporary_values_choose_accession"
              ) |>
              group_by(!!sym(protein_id_column), sample_id) |>
              summarise(
                is_na = sum(is.na(temporary_values_choose_accession)),
                temporary_values_choose_accession = case_when(
                  all(is.na(temporary_values_choose_accession)) ~ NA_real_,
                  aggregation_method == "sum" ~ sum(temporary_values_choose_accession, na.rm = TRUE),
                  aggregation_method == "mean" ~ mean(temporary_values_choose_accession, na.rm = TRUE),
                  aggregation_method == "median" ~ median(temporary_values_choose_accession, na.rm = TRUE)
                ),
                num_values = n()
              ) |>
              ungroup() |>
              pivot_wider(
                id_cols = !!sym(protein_id_column),
                names_from = sample_id,
                values_from = temporary_values_choose_accession,
                values_fill = NA_real_
              )

            if(replace_zero_with_na == TRUE) {
              summed_data[is.na(summed_data)] <- NA
            }
            
            cat(">>> STEP 5: CALLING rankProteinAccessionHelper <<<\n")
            cat(sprintf("   Passing 'delim' to helper: '%s'\n", delim))
            cat(sprintf("   Number of protein_id_table rows to rank: %d\n", nrow(protein_id_table)))
            if (nrow(protein_id_table) > 0) {
              list_col_name <- paste0(protein_id_column, "_list")
              cat(sprintf("   First 5 protein ID lists: %s\n", paste(head(protein_id_table[[list_col_name]], 5), collapse = ", ")))
            }
            cat("\n")

            protein_id_table <- rankProteinAccessionHelper(input_tbl = protein_id_table,
                                                         acc_detail_tab = seqinr_obj,
                                                         accessions_column = !!sym(paste0(protein_id_column, "_list")),
                                                         row_id_column = seqinr_accession_column,
                                                         group_id = !!sym(protein_id_column),
                                                         delim = delim) |>
              dplyr::rename(!!sym(paste0(protein_id_column, "_list")) := seqinr_accession_column) |>
              dplyr::select(-num_gene_names, -gene_names, -is_unique)
            
            cat(">>> STEP 6: RETURNED FROM rankProteinAccessionHelper <<<\n")
            cat(sprintf("   Number of rows returned: %d\n", nrow(protein_id_table)))
            cat("\n")

            theObject@protein_id_table <- protein_id_table
            theObject@protein_quant_table <- summed_data[, colnames(protein_quant_table)]
            
            cat(">>> STEP 7: FINAL OUTPUT <<<\n")
            cat(sprintf("   Final protein_quant_table dimensions: %d rows x %d cols\n", nrow(theObject@protein_quant_table), ncol(theObject@protein_quant_table)))
            cat(sprintf("   Final protein_id_table dimensions: %d rows x %d cols\n", nrow(theObject@protein_id_table), ncol(theObject@protein_id_table)))
            cat(sprintf("   First 5 final Protein IDs: %s\n", paste(head(theObject@protein_quant_table[[protein_id_column]], 5), collapse = ", ")))
            cat("\n")
            cat("╔═══════════════════════════════════════════════════════════════════════════╗\n")
            cat("║           EXITING chooseBestProteinAccession (DEBUG66)                ║\n")
            cat("╚═══════════════════════════════════════════════════════════════════════════╝\n\n")

            return(theObject)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="chooseBestProteinAccessionSumDuplicates"
           , def=function( theObject, delim, quant_columns_pattern, islogged ) {
             standardGeneric("chooseBestProteinAccessionSumDuplicates")
           }
           , signature=c("theObject", "delim", "quant_columns_pattern", "islogged" ))

#'@export
setMethod( f = "chooseBestProteinAccessionSumDuplicates"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, delim=";", quant_columns_pattern = "\\d+", islogged = TRUE ) {

             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column

             protein_log2_quant_cln <- protein_quant_table |>
               mutate( !!sym(protein_id_column) := str_split_i(!!sym( protein_id_column), delim, 1 ) ) |>
               group_by( !!sym(protein_id_column) ) |>
               summarise ( across( matches(quant_columns_pattern)
                                   , \(x){ if(islogged==TRUE) {
                                                log2(sum(2^x, na.rm = TRUE))
                                           } else {
                                                sum(x, na.rm = TRUE)
                                           }
                                         } )) |>
               ungroup()

             theObject@protein_quant_table <- protein_log2_quant_cln

             theObject

           })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="filterSamplesByProteinCorrelationThreshold"
           , def=function( theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = NULL ) {
             standardGeneric("filterSamplesByProteinCorrelationThreshold")
           }
           , signature=c("theObject", "pearson_correlation_per_pair", "min_pearson_correlation_threshold" ))

#'@export
setMethod( f = "filterSamplesByProteinCorrelationThreshold"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = NULL  ) {

             pearson_correlation_per_pair <- checkParamsObjectFunctionSimplify( theObject
                                                                           , "pearson_correlation_per_pair"
                                                                           , default_value = NULL)
             min_pearson_correlation_threshold <- checkParamsObjectFunctionSimplify( theObject
                                                                                , "min_pearson_correlation_threshold"
                                                                                , default_value = 0.75)

             theObject <- updateParamInObject(theObject, "pearson_correlation_per_pair")
             theObject <- updateParamInObject(theObject, "min_pearson_correlation_threshold")

             filtered_table <- filterSamplesByProteinCorrelationThresholdHelper (
               pearson_correlation_per_pair
               , protein_intensity_table = theObject@protein_quant_table
               , min_pearson_correlation_threshold = min_pearson_correlation_threshold
               , filename_column_x = !!sym( paste0( theObject@sample_id, ".x") )
               , filename_column_y = !!sym( paste0( theObject@sample_id, ".y") )
               , protein_id_column = theObject@protein_id_column
               , correlation_column = pearson_correlation )

             theObject@protein_quant_table <- filtered_table

             theObject <- cleanDesignMatrix(theObject)

             theObject
             })

#Format the design matrix so that only metadata for samples in the protein data are retained, and also
# sort the sample IDs in the same order as the data matrix

#'@export
setMethod( f ="cleanDesignMatrix"
           , signature = "ProteinQuantitativeData"
           , definition=function( theObject ) {

            samples_id_vector <- setdiff(colnames(theObject@protein_quant_table), theObject@sample_id )

             theObject@design_matrix <- data.frame( temp_sample_id = samples_id_vector )  |>
               inner_join( theObject@design_matrix
                          , by = join_by ( temp_sample_id == !!sym(theObject@sample_id)) ) |>
               dplyr::rename( !!sym(theObject@sample_id) := "temp_sample_id" ) |>
               dplyr::filter( !!sym( theObject@sample_id) %in% samples_id_vector )


             return(theObject)
           })
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod( f="proteinIntensityFiltering"
           , signature="ProteinQuantitativeData"
           , definition = function( theObject
                                    , proteins_intensity_cutoff_percentile = NULL
                                    , proteins_proportion_of_samples_below_cutoff = NULL
                                    , core_utilisation = NULL) {
             protein_quant_table <- theObject@protein_quant_table

             proteins_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject
                                                                               , "proteins_intensity_cutoff_percentile"
                                                                               , NULL)
             proteins_proportion_of_samples_below_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                               , "proteins_proportion_of_samples_below_cutoff"
                                                                               , NULL)
             core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                                    , "core_utilisation"
                                                                    , NA)

             theObject <- updateParamInObject(theObject, "proteins_intensity_cutoff_percentile")
             theObject <- updateParamInObject(theObject, "proteins_proportion_of_samples_below_cutoff")
             theObject <- updateParamInObject(theObject, "core_utilisation")


             data_long_cln <- protein_quant_table  |>
               pivot_longer( cols=!matches(theObject@protein_id_column)
                             , names_to = theObject@sample_id
                             , values_to = "log_values")  |>
               mutate( temp = "")

             min_peptide_intensity_threshold <- ceiling( quantile( data_long_cln$log_values, na.rm=TRUE, probs = c(proteins_intensity_cutoff_percentile) ))[1]

             peptide_normalised_pif_cln <- peptideIntensityFilteringHelper( data_long_cln
                                                                      , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                                                      , proteins_proportion_of_samples_below_cutoff = proteins_proportion_of_samples_below_cutoff
                                                                      , protein_id_column = !!sym( theObject@protein_id_column)
                                                                      , peptide_sequence_column = temp
                                                                      , peptide_quantity_column = log_values
                                                                      , core_utilisation = core_utilisation)


             theObject@protein_quant_table <- peptide_normalised_pif_cln |>
               dplyr::select( -temp) |>
               pivot_wider( id_cols = theObject@protein_id_column , names_from = !!sym(theObject@sample_id), values_from = log_values)

             theObject <- cleanDesignMatrix(theObject)

             updated_object <- theObject

          return(updated_object)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## MISSING METHODS COPIED FROM proteins4.R - Added during golem audit fix
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod( f ="setProteinData"
           , signature = "ProteinQuantitativeData"
            , definition=function( theObject, protein_quant_table, protein_id_column ) {
              theObject@protein_quant_table <- protein_quant_table
              theObject@protein_id_column <- protein_id_column

              return(theObject)
            })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod(f="plotRle"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            design_matrix <- as.data.frame(design_matrix)

            if(!is.null(sample_label)) {
              if ( sample_label %in% colnames(design_matrix)) {
                rownames( design_matrix) <- design_matrix[,sample_label]
                colnames( frozen_protein_matrix ) <- design_matrix[,sample_label]

              } } else {
                rownames( design_matrix) <- design_matrix[,sample_id]
              }

            # print( design_matrix)

            rowinfo_vector <- NA
            if( !is.na(grouping_variable)){
              rowinfo_vector <-  design_matrix[colnames(frozen_protein_matrix), grouping_variable]
            }

            print(rownames( design_matrix))
            print(colnames( frozen_protein_matrix))
            print(rowinfo_vector)
              rle_plot_before_cyclic_loess <- plotRleHelper( t(frozen_protein_matrix)
                                                       , rowinfo = rowinfo_vector
                                                       , yaxis_limit = yaxis_limit)

            return( rle_plot_before_cyclic_loess)

          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setMethod(f="plotRleList"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, list_of_columns, yaxis_limit = c()) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            design_matrix <- as.data.frame(design_matrix)
            rownames( design_matrix) <- design_matrix[,sample_id]

            # print( design_matrix)

            runOneRle <- function( column_name) {
              rowinfo_vector <- NA

              if ( column_name %in% colnames(design_matrix) ) {
                rowinfo_vector <- design_matrix[colnames(frozen_protein_matrix), column_name]
              }

              rle_plot_before_cyclic_loess <- plotRleHelper( t(frozen_protein_matrix)
                                                             , rowinfo = rowinfo_vector
                                                             , yaxis_limit = yaxis_limit)

              return( rle_plot_before_cyclic_loess)
            }

            list_of_rle_plots <- purrr::map( list_of_columns, runOneRle)

            names(list_of_rle_plots) <- list_of_columns

            return( list_of_rle_plots)

          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
savePlotRleList <- function( input_list, prefix = "RLE", suffix = c("png", "pdf"), output_dir ) {

  list_of_filenames <- expand_grid( column=names(input_list), suffix=suffix)  |>
    mutate( filename= paste0( "RLE", "_", column , ".", suffix))  |>
    left_join( tibble( column =names( input_list)
                       ,  plots= input_list )
               , by=join_by(column ))


  purrr::walk2( list_of_filenames$plots
                , list_of_filenames$filename
                , \(.x, .y){
                  ggsave( plot=.x, filename= file.path(output_dir, .y))
                } )

  list_of_filenames

}



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod(f="plotPca"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size=8, cv_percentile = 0.90) {
            # Defensive checks
            if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
              stop("grouping_variable must be a single character string")
            }
            
            if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
              stop("shape_variable must be NULL or a single character string")
            }
            
            if (!grouping_variable %in% colnames(theObject@design_matrix)) {
              stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
            }
            
            if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
              stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
            }

            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix_pca <- frozen_protein_matrix
            frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

            if(is.na(label_column) || label_column == "") {
              label_column <- ""
            }

            required_cols <- c(sample_id, grouping_variable)
            if (!is.null(shape_variable)) {
              required_cols <- c(required_cols, shape_variable)
            }
            missing_cols <- setdiff(required_cols, colnames(design_matrix))
            if (length(missing_cols) > 0) {
              stop(sprintf("Missing columns in design matrix: %s", paste(missing_cols, collapse = ", ")))
            }

            tryCatch({
              pca_plot <- plotPcaHelper(frozen_protein_matrix_pca,
                                           design_matrix,
                                           sample_id_column = sample_id,
                                           grouping_variable = grouping_variable,
                                           shape_variable = shape_variable,
                                           label_column = label_column,
                                           title = title,
                                           geom.text.size = font_size)
              return(pca_plot)
            }, error = function(e) {
              stop(sprintf("Error in plotPcaHelper: %s", e$message))
            })
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setMethod(f="plotPcaList"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, grouping_variables_list, label_column, title, font_size=8, cv_percentile = 0.90) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix_pca <- frozen_protein_matrix
            frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

            if( is.na(label_column) || label_column == "") {
              label_column <- ""
            }

            pca_plots_list <- plotPcaListHelper( frozen_protein_matrix_pca
                                                 , design_matrix
                                                 , sample_id_column =  sample_id
                                                 , grouping_variables_list = grouping_variables_list
                                                 , label_column =  label_column
                                                 , title = title
                                                 , geom.text.size = font_size )

            return( pca_plots_list)
          })


#' @export
savePlotPcaList <- function( input_list, prefix = "PCA", suffix = c("png", "pdf"), output_dir ) {

  list_of_filenames <- expand_grid( column=names(input_list), suffix=suffix)  |>
    mutate( filename= paste0( "RLE", "_", column , ".", suffix))  |>
    left_join( tibble( column =names( input_list)
                       ,  plots= input_list )
               , by=join_by(column ))


  purrr::walk2( list_of_filenames$plots
                , list_of_filenames$filename
                , \(.x, .y){
                  ggsave( plot=.x, filename= file.path(output_dir, .y))
                } )

  list_of_filenames

}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setMethod(f="getPcaMatrix"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id


            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix_pca <- frozen_protein_matrix
            frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA


            pca_mixomics_before_cyclic_loess <- mixOmics::pca(t(as.matrix(frozen_protein_matrix_pca)))$variates$X |>
              as.data.frame()    |>
              rownames_to_column(sample_id)  |>
              left_join(design_matrix, by = sample_id  )


            return( pca_mixomics_before_cyclic_loess)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setMethod( f = "proteinTechRepCorrelation"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject,  tech_rep_num_column = NULL, tech_rep_remove_regex = NULL ) {
             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             sample_id <- theObject@sample_id
             tech_rep_column <- theObject@technical_replicate_id

             tech_rep_num_column <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_num_column", NULL)
             tech_rep_remove_regex <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", NULL)

             theObject <- updateParamInObject(theObject, "tech_rep_num_column")
             theObject <- updateParamInObject(theObject, "tech_rep_remove_regex")

             frozen_protein_matrix <- protein_quant_table |>
               column_to_rownames(protein_id_column) |>
               as.matrix()

             frozen_protein_matrix_pca <- frozen_protein_matrix
             frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

             protein_matrix_tech_rep <-proteinTechRepCorrelationHelper( design_matrix, frozen_protein_matrix_pca
                                                                        , protein_id_column = protein_id_column
                                                                        , sample_id_column=sample_id
                                                                        , tech_rep_column = tech_rep_column
                                                                        , tech_rep_num_column = tech_rep_num_column
                                                                        , tech_rep_remove_regex = tech_rep_remove_regex )

             return( protein_matrix_tech_rep )
           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Plot Pearson Correlation
#' @param theObject is an object of the type ProteinQuantitativeData
#' @param tech_rep_remove_regex samples containing this string are removed from correlation analysis (e.g. if you have lots of pooled sample and want to remove them)
#' @param correlation_group is the group where every pair of samples are compared
#' @export
setMethod(f="plotPearson",
          signature="ProteinQuantitativeData",
          definition=function(theObject, tech_rep_remove_regex = "pool", correlation_group = NA) {

            correlation_group_to_use <- correlation_group

            if( is.na( correlation_group)) {
              correlation_group_to_use <- theObject@technical_replicate_id
            }

            correlation_vec <- pearsonCorForSamplePairs(theObject
                                                        , tech_rep_remove_regex
                                                        , correlation_group = correlation_group_to_use)

            pearson_plot <- correlation_vec |>
              ggplot(aes(pearson_correlation)) +
              geom_histogram(breaks = seq(min(round(correlation_vec$pearson_correlation - 0.5, 2), na.rm = TRUE), 1, 0.001)) +
              scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4)) +
              xlab("Pearson Correlation") +
              ylab("Counts") +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank())

            return(pearson_plot)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create empty QC Grid
#' @export
setClass("GridPlotData",
         slots = list(
           pca_plots = "list",
           density_plots = "list",
           rle_plots = "list",
           pearson_plots = "list",
           cancor_plots = "list",
           pca_titles = "list",
           density_titles = "list",
           rle_titles = "list",
           pearson_titles = "list",
           cancor_titles = "list"
         ))

#' @export
setGeneric("InitialiseGrid", function(dummy = NULL) {
  standardGeneric("InitialiseGrid")
})

#' @export
setMethod("InitialiseGrid", 
          signature(dummy = "ANY"),
          function(dummy = NULL) {
            new("GridPlotData",
                pca_plots = list(),
                density_plots = list(),
                rle_plots = list(),
                pearson_plots = list(),
                cancor_plots = list(),
                pca_titles = list(),
                density_titles = list(),
                rle_titles = list(),
                pearson_titles = list(),
                cancor_titles = list())
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create a QC composite figure

#' @export
setGeneric(name = "createGridQC",
           def = function(theObject, pca_titles = NULL, density_titles = NULL, rle_titles = NULL, pearson_titles = NULL, cancor_titles = NULL, ncol = 3, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
             standardGeneric("createGridQC")
           },
           signature = c("theObject"))

#' @export
setMethod(f = "createGridQC",
          signature = "GridPlotData",
          definition = function(theObject, pca_titles = NULL, density_titles = NULL, rle_titles = NULL, pearson_titles = NULL, cancor_titles = NULL, ncol = 3, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
            
            # Use stored titles if not provided as parameters
            pca_titles <- if(is.null(pca_titles)) theObject@pca_titles else pca_titles
            density_titles <- if(is.null(density_titles)) theObject@density_titles else density_titles
            rle_titles <- if(is.null(rle_titles)) theObject@rle_titles else rle_titles
            pearson_titles <- if(is.null(pearson_titles)) theObject@pearson_titles else pearson_titles
            cancor_titles <- if(is.null(cancor_titles)) theObject@cancor_titles else cancor_titles
            
            createLabelPlot <- function(title) {
              # Option 1: Use xlim to expand the plot area and position text at left edge
              ggplot() + 
                annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
                xlim(0, 1) +  # Explicitly set the x limits
                theme_void() +
                theme(
                  plot.margin = margin(5, 5, 5, 5),
                  panel.background = element_blank()
                )
            }
            
            # Create basic plots without titles
            createPcaPlot <- function(plot) {
              plot +
                xlim(-40, 45) + ylim(-30, 25) +
                theme(text = element_text(size = 15),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank())
            }
            
            createDensityPlot <- function(plot) {
              # For all plots, just apply the theme without adding title
              if (inherits(plot, "patchwork")) {
                plot & 
                  theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    text = element_text(size = 15)
                  )
              } else {
                plot +
                  theme(text = element_text(size = 15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank())
              }
            }
            
            createRlePlot <- function(plot) {
              plot +
                theme(text = element_text(size = 15),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
            }
            
            createPearsonPlot <- function(plot) {
              plot +
                theme(text = element_text(size = 15))
            }
            
            createCancorPlot <- function(plot) {
              plot +
                theme(text = element_text(size = 15),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank())
            }
            
            # Create plots without titles - FIXED ORDER
            # Define the correct order: before_cyclic_loess, before_ruvIIIc, after_ruvIIIc
            plot_order <- c("pca_plot_before_cyclic_loess_group", "pca_plot_before_ruvIIIc_group", "pca_plot_after_ruvIIIc_group")
            density_order <- c("density_plot_before_cyclic_loess_group", "density_plot_before_ruvIIIc_group", "density_plot_after_ruvIIIc_group")
            rle_order <- c("rle_plot_before_cyclic_loess_group", "rle_plot_before_ruvIIIc_group", "rle_plot_after_ruvIIIc_group")
            pearson_order <- c("pearson_correlation_pair_before_cyclic_loess", "pearson_correlation_pair_before_ruvIIIc", "pearson_correlation_pair_after_ruvIIIc_group")
            cancor_order <- c("cancor_plot_before_cyclic_loess", "cancor_plot_before_ruvIIIc", "cancor_plot_after_ruvIIIc")
            
            # Extract plots in the correct order using lapply
            created_pca_plots <- lapply(plot_order, function(name) {
              if (!is.null(theObject@pca_plots[[name]])) createPcaPlot(theObject@pca_plots[[name]]) else NULL
            })
            created_pca_plots <- created_pca_plots[!sapply(created_pca_plots, is.null)]
            
            created_density_plots <- lapply(density_order, function(name) {
              if (!is.null(theObject@density_plots[[name]])) createDensityPlot(theObject@density_plots[[name]]) else NULL
            })
            created_density_plots <- created_density_plots[!sapply(created_density_plots, is.null)]
            
            created_rle_plots <- lapply(rle_order, function(name) {
              if (!is.null(theObject@rle_plots[[name]])) createRlePlot(theObject@rle_plots[[name]]) else NULL
            })
            created_rle_plots <- created_rle_plots[!sapply(created_rle_plots, is.null)]
            
            created_pearson_plots <- lapply(pearson_order, function(name) {
              if (!is.null(theObject@pearson_plots[[name]])) createPearsonPlot(theObject@pearson_plots[[name]]) else NULL
            })
            created_pearson_plots <- created_pearson_plots[!sapply(created_pearson_plots, is.null)]
            
            created_cancor_plots <- lapply(cancor_order, function(name) {
              if (!is.null(theObject@cancor_plots[[name]])) createCancorPlot(theObject@cancor_plots[[name]]) else NULL
            })
            # Don't filter out NULL plots to maintain column alignment
            
            # Create label plots
            pca_labels <- lapply(pca_titles, createLabelPlot)
            density_labels <- lapply(density_titles, createLabelPlot)
            rle_labels <- lapply(rle_titles, createLabelPlot)
            pearson_labels <- lapply(pearson_titles, createLabelPlot)
            cancor_labels <- lapply(cancor_titles, createLabelPlot)
            
            # Combine with labels above each row - modified to keep legends with their plots
            plot_sections <- list()
            height_values <- c()
            
            # Add PCA plots if they exist
            if(length(theObject@pca_plots) > 0) {
              plot_sections <- append(plot_sections, list(
                wrap_plots(pca_labels, ncol = ncol),
                wrap_plots(created_pca_plots, ncol = ncol)
              ))
              height_values <- c(height_values, 0.1, 1)
            }
            
            # Add Density plots if they exist
            if(length(theObject@density_plots) > 0) {
              plot_sections <- append(plot_sections, list(
                wrap_plots(density_labels, ncol = ncol),
                wrap_plots(created_density_plots, ncol = ncol)
              ))
              height_values <- c(height_values, 0.1, 1)
            }
            
            # Add RLE plots if they exist
            if(length(theObject@rle_plots) > 0) {
              plot_sections <- append(plot_sections, list(
                wrap_plots(rle_labels, ncol = ncol),
                wrap_plots(created_rle_plots, ncol = ncol)
              ))
              height_values <- c(height_values, 0.1, 1)
            }
            
            # Add Pearson plots if they exist
            if(length(theObject@pearson_plots) > 0) {
              plot_sections <- append(plot_sections, list(
                wrap_plots(pearson_labels, ncol = ncol),
                wrap_plots(created_pearson_plots, ncol = ncol)
              ))
              height_values <- c(height_values, 0.1, 1)
            }
            
            # Add Cancor plots if they exist (check for any non-NULL plots)
            if(length(theObject@cancor_plots) > 0 && any(!sapply(created_cancor_plots, is.null))) {
              # Replace NULL plots with empty plots to maintain column alignment
              cancor_plots_aligned <- lapply(created_cancor_plots, function(plot) {
                if (is.null(plot)) {
                  ggplot() + theme_void()  # Empty plot for NULL positions
                } else {
                  plot
                }
              })
              
              plot_sections <- append(plot_sections, list(
                wrap_plots(cancor_labels, ncol = ncol),
                wrap_plots(cancor_plots_aligned, ncol = ncol)
              ))
              height_values <- c(height_values, 0.1, 1)
            }
            
            # Create combined plot from sections
            combined_plot <- wrap_plots(plot_sections, ncol = 1) +
              plot_layout(heights = height_values)

            if (!is.null(save_path)) {
              # Calculate dynamic width based on number of columns
              plot_width <- 4 + (ncol * 3)  # Base width + 3 units per column
              plot_height <- 4 + (length(height_values) * 2)  # Base height + 2 units per row
              
              sapply(c("png", "pdf", "svg"), function(ext) {
                ggsave(
                  plot = combined_plot,
                  filename = file.path(save_path, paste0(file_name, ".", ext)),
                  width = plot_width,
                  height = plot_height
                )
              })
              message(paste("Plots saved in", save_path))
            }
            
            return(combined_plot)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@title Normalise between Arrays
#'@export
#'@param theObject Object of class ProteinQuantitativeData
#'@param normalisation_method Method to use for normalisation. Options are cyclicloess, quantile, scale, none
setMethod(f="normaliseBetweenSamples"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject,  normalisation_method= NULL) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            normalisation_method <- checkParamsObjectFunctionSimplify( theObject
                                                                       , "normalisation_method"
                                                                       , "cyclicloess")

            theObject <- updateParamInObject(theObject, "normalisation_method")

            frozen_protein_matrix <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            frozen_protein_matrix[!is.finite(frozen_protein_matrix)] <- NA

            normalised_frozen_protein_matrix <- frozen_protein_matrix

            print(paste0("normalisation_method = ", normalisation_method))

            switch( normalisation_method
                    , cyclicloess = {
                      normalised_frozen_protein_matrix <- normalizeCyclicLoess( frozen_protein_matrix )
                    }
                    , quantile = {
                      normalised_frozen_protein_matrix <- normalizeQuantiles( frozen_protein_matrix  )
                    }
                    , scale = {
                      normalised_frozen_protein_matrix <- normalizeMedianAbsValues( frozen_protein_matrix  )
                    }
                    , none = {
                      normalised_frozen_protein_matrix <- frozen_protein_matrix
                    }
            )

            normalised_frozen_protein_matrix[!is.finite(normalised_frozen_protein_matrix)] <- NA

            # normalised_frozen_protein_matrix_filt <- as.data.frame( normalised_frozen_protein_matrix ) |>
            #   dplyr::filter( if_all( everything(), \(x) { !is.na(x) } ) ) |>
            #   as.matrix()

            theObject@protein_quant_table <- normalised_frozen_protein_matrix |>
                      as.data.frame() |>
                      rownames_to_column(protein_id_column)

            theObject <- cleanDesignMatrix(theObject)

            updated_object <- theObject

            return(updated_object)

          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Pearson Correlation for Sample Pairs
#' @param theObject is an object of the type ProteinQuantitativeData
#' @param tech_rep_remove_regex samples containing this string are removed from correlation analysis (e.g. if you have lots of pooled sample and want to remove them)
#' @param correlation_group is the group where every pair of samples are compared
#'@export
setMethod(f="pearsonCorForSamplePairs"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject, tech_rep_remove_regex = NULL, correlation_group = NA ) {
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            replicate_group_column <- theObject@technical_replicate_id
            if(!is.na( correlation_group )) {
              replicate_group_column <- correlation_group
            }

            tech_rep_remove_regex <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", "pool")
            theObject <- updateParamInObject(theObject, "tech_rep_remove_regex")

            frozen_mat_pca_long <- protein_quant_table |>
              pivot_longer( cols=!matches(protein_id_column)
                            , values_to = "Protein.normalised"
                            , names_to = sample_id) |>
              left_join( design_matrix
                         , by = join_by( !!sym(sample_id) == !!sym(sample_id))) |>
              mutate( temp = "")

            # Detect available cores and configure parallel processing
            num_available_cores <- parallel::detectCores()
            num_workers <- min(max(1, num_available_cores - 1), 8)  # Use n-1 cores, max 8
            message(sprintf("*** PEARSON: Detected %d CPU cores, using %d workers for parallel processing ***", 
                            num_available_cores, num_workers))
            
            # Calculate expected pair count for user information
            sample_count <- nrow(design_matrix |> dplyr::select(!!sym(sample_id), !!sym(replicate_group_column)))
            estimated_pairs <- choose(sample_count, 2)
            message(sprintf("*** PEARSON: Processing %d samples (~%d pairwise correlations) ***", 
                            sample_count, estimated_pairs))
            message("*** PEARSON: Starting pairwise correlation calculations... ***")

            correlation_results_before_cyclic_loess <- calulatePearsonCorrelationForSamplePairsHelper( design_matrix |>
                                                                                                         dplyr::select( !!sym(sample_id), !!sym(replicate_group_column) )
                                                                                                       , run_id_column = sample_id
                                                                                                       , replicate_group_column = replicate_group_column
                                                                                                       , frozen_mat_pca_long
                                                                                                       , num_of_cores = num_workers
                                                                                                       , sample_id_column = sample_id
                                                                                                       , protein_id_column = protein_id_column
                                                                                                       , peptide_sequence_column = "temp"
                                                                                                       , peptide_normalised_column = "Protein.normalised")
            
            message("*** PEARSON: Correlation calculations complete ***")

            correlation_vec_before_cyclic_loess <- correlation_results_before_cyclic_loess |>
              dplyr::filter( !str_detect(!!sym(replicate_group_column), tech_rep_remove_regex )  )

           return( correlation_vec_before_cyclic_loess)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="getNegCtrlProtAnova"
           , def=function( theObject
                           , ruv_grouping_variable  = NULL
                           , percentage_as_neg_ctrl  = NULL
                           , num_neg_ctrl  = NULL
                           , ruv_qval_cutoff = NULL
                           , ruv_fdr_method = NULL ) {
             standardGeneric("getNegCtrlProtAnova")
           }
           , signature=c("theObject", "ruv_grouping_variable", "num_neg_ctrl", "ruv_qval_cutoff", "ruv_fdr_method"))

#'@export
setMethod(f="getNegCtrlProtAnova"
          , signature="ProteinQuantitativeData"
          , definition=function( theObject
                                 , ruv_grouping_variable = NULL
                                 , percentage_as_neg_ctrl = NULL
                                 , num_neg_ctrl = NULL
                                 , ruv_qval_cutoff = NULL
                                 , ruv_fdr_method = NULL ) {

            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            design_matrix <- theObject@design_matrix
            group_id <- theObject@group_id
            sample_id <- theObject@sample_id

            normalised_frozen_protein_matrix_filt <- protein_quant_table |>
              column_to_rownames(protein_id_column) |>
              as.matrix()

            ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", "replicates")
            percentage_as_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject, "percentage_as_neg_ctrl", 10)
            num_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject
                                                               , "num_neg_ctrl"
                                                               , round(nrow( theObject@protein_quant_table) * percentage_as_neg_ctrl / 100, 0))
            ruv_qval_cutoff <- checkParamsObjectFunctionSimplify( theObject, "ruv_qval_cutoff", 0.05)
            ruv_fdr_method <- checkParamsObjectFunctionSimplify( theObject, "ruv_fdr_method", "BH")

            theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
            theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
            theObject <- updateParamInObject(theObject, "num_neg_ctrl")
            theObject <- updateParamInObject(theObject, "ruv_qval_cutoff")
            theObject <- updateParamInObject(theObject, "ruv_fdr_method")

            control_genes_index <- getNegCtrlProtAnovaHelper( normalised_frozen_protein_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id)) ]
                                                        , design_matrix = design_matrix |>
                                                          column_to_rownames(sample_id) |>
                                                          dplyr::select( -!!sym(group_id))
                                                        , grouping_variable = ruv_grouping_variable
                                                        , percentage_as_neg_ctrl = percentage_as_neg_ctrl
                                                        , num_neg_ctrl = num_neg_ctrl
                                                        , ruv_qval_cutoff = ruv_qval_cutoff
                                                        , ruv_fdr_method = ruv_fdr_method )

            return(control_genes_index)
          })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod( f = "ruvCancor"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL) {
             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id

             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)
             num_components_to_impute <- checkParamsObjectFunctionSimplify( theObject, "num_components_to_impute", 2)
             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)

             theObject <- updateParamInObject(theObject, "ctrl")
             theObject <- updateParamInObject(theObject, "num_components_to_impute")
             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")

             if(! ruv_grouping_variable %in% colnames(design_matrix)) {
               stop( paste0("The 'ruv_grouping_variable = "
                            , ruv_grouping_variable
                            , "' is not a column in the design matrix.") )
             }

             if( is.na(num_components_to_impute) || num_components_to_impute < 1) {
               stop(paste0("The num_components_to_impute = ", num_components_to_impute, " value is invalid."))
             }

             if( length( ctrl) < 5 ) {
               stop(paste0( "The number of negative control molecules entered is less than 5. Please check the 'ctl' parameter."))
             }

             normalised_frozen_protein_matrix_filt <- protein_quant_table |>
               column_to_rownames(protein_id_column) |>
               as.matrix()

             Y <-  t( normalised_frozen_protein_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id))])
             if( length(which( is.na(normalised_frozen_protein_matrix_filt) )) > 0 ) {
               Y <- impute.nipals( t( normalised_frozen_protein_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id))])
                                   , ncomp=num_components_to_impute)
             }

             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                dplyr::pull(!!sym(ruv_grouping_variable)),
                                              ctl = ctrl)
             cancorplot_r2


           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setGeneric(name="getRuvIIIReplicateMatrix"
           , def=function( theObject,  ruv_grouping_variable = NULL) {
             standardGeneric("getRuvIIIReplicateMatrix")
           }
           , signature=c("theObject", "ruv_grouping_variable"))

#'@export
setMethod( f = "getRuvIIIReplicateMatrix"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, ruv_grouping_variable = NULL) {
             protein_quant_table <- theObject@protein_quant_table
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)

             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")

             ruvIII_replicates_matrix <- getRuvIIIReplicateMatrixHelper( design_matrix
                                                                   , !!sym(sample_id)
                                                                   , !!sym(ruv_grouping_variable))
             return( ruvIII_replicates_matrix)
           })


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setMethod( f = "ruvIII_C_Varying"
           , signature="ProteinQuantitativeData"
           , definition=function( theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL) {
             protein_quant_table <- theObject@protein_quant_table
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

             normalised_frozen_protein_matrix_filt <- protein_quant_table |>
               column_to_rownames(protein_id_column) |>
               as.matrix()

             Y <-  t( normalised_frozen_protein_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id))])

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

             ruv_normalised_results_cln <- cln_mat_4 |>
               as.data.frame() |>
               rownames_to_column(protein_id_column)

             theObject@protein_quant_table <- ruv_normalised_results_cln

             theObject <- cleanDesignMatrix(theObject)

             return( theObject )
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setGeneric(name="preservePeptideNaValues"
           , def=function( peptide_obj, protein_obj)  {
             standardGeneric("preservePeptideNaValues")
           }
           , signature=c("peptide_obj", "protein_obj" ))

#'@export
setMethod( f = "preservePeptideNaValues"
           , signature=c( "PeptideQuantitativeData", "ProteinQuantitativeData" )
           , definition= function( peptide_obj, protein_obj) {
             preservePeptideNaValuesHelper( peptide_obj, protein_obj)
           })

#'@export
preservePeptideNaValuesHelper <- function( peptide_obj, protein_obj) {

  sample_id_column <- peptide_obj@sample_id
  protein_id_column <- peptide_obj@protein_id_column

  check_peptide_value <- peptide_obj@peptide_data |>
    group_by( !!sym( sample_id_column), !!sym(protein_id_column) ) |>
    summarise( Peptide.Normalised = sum( Peptide.Normalised, na.rm=TRUE)
               , is_na = sum( is.na(Peptide.Normalised ))
               , num_values = n() ) |>
    mutate( Peptide.Normalised = if_else( is_na == num_values, NA_real_, Peptide.Normalised)) |>
    ungroup() |>
    arrange( !!sym( sample_id_column)) |>
    pivot_wider( id_cols = !!sym(protein_id_column)
                 , names_from = !!sym(sample_id_column)
                 , values_from = Peptide.Normalised
                 , values_fill = NA_real_)

  check_peptide_value_cln <- check_peptide_value[rownames(protein_obj@protein_quant_table)
                                                 , colnames(  protein_obj@protein_quant_table)]

  if( length( which (rownames(protein_obj@protein_quant_table) ==  rownames(check_peptide_value_cln))) != nrow(check_peptide_value_cln) ) {
    stop("The rows in the protein object and the peptide object do not match")
  }

  if( length( which( colnames( protein_obj@protein_quant_table) == colnames(check_peptide_value_cln) )) != ncol(check_peptide_value_cln) ) {
    stop("The columns in the protein object and the peptide object do not match")
  }

  protein_obj@protein_quant_table [is.na(check_peptide_value_cln)] <- NA

  protein_obj
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod(f="plotPcaBox"
          , signature="gg"
          , definition=function(theObject, grouping_variable, title = "", font_size = 8) {
            # For gg class objects, create a copy and change its class to ggplot
            gg_obj <- theObject
            class(gg_obj) <- "ggplot"
            
            # Then call the ggplot method
            plotPcaBox(gg_obj, grouping_variable, title, font_size)
          })

#'@export
setMethod(f="plotPcaBox"
          , signature="ggplot"
          , definition=function(theObject, grouping_variable, title = "", font_size = 8) {
            # First try to get data directly from the ggplot object's data element
            if (!is.null(theObject$data) && is.data.frame(theObject$data)) {
              pca_data <- as_tibble(theObject$data)
            } else {
              # Fall back to other extraction methods
              pca_data <- as_tibble(ggplot_build(theObject)$data[[1]])
              
              # If the data doesn't have PC1/PC2, try to extract from the plot's environment
              if (!("PC1" %in% colnames(pca_data) && "PC2" %in% colnames(pca_data))) {
                # Try to get the data from the plot's environment
                if (exists("data", envir = environment(theObject$mapping$x))) {
                  pca_data <- as_tibble(get("data", envir = environment(theObject$mapping$x)))
                } else {
                  stop("Could not extract PCA data from the ggplot object")
                }
              }
            }
            
            # Check if grouping variable exists in the data
            if (!grouping_variable %in% colnames(pca_data)) {
              stop(sprintf("grouping_variable '%s' not found in the data", grouping_variable))
            }
            
            # Create PC1 boxplot
            pc1_box <- ggplot(pca_data, aes(x = !!sym(grouping_variable), y = PC1, fill = !!sym(grouping_variable))) +
              geom_boxplot(notch = TRUE) +
              theme_bw() +
              labs(title = title,
                   x = "",
                   y = "PC1") +
              theme(
                legend.position = "none",
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                text = element_text(size = font_size),
                plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )
            
            # Add explicit fill scale to support >6 discrete levels
            categorical_colors <- getCategoricalColourPalette()
            pc1_box <- pc1_box + scale_fill_manual(values = categorical_colors)
            
            # Create PC2 boxplot
            pc2_box <- ggplot(pca_data, aes(x = !!sym(grouping_variable), y = PC2, fill = !!sym(grouping_variable))) +
              geom_boxplot(notch = TRUE) +
              theme_bw() +
              labs(x = "",
                   y = "PC2") +
              theme(
                legend.position = "none",
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                text = element_text(size = font_size),
                plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )
            
            # Add explicit fill scale to support >6 discrete levels
            pc2_box <- pc2_box + scale_fill_manual(values = categorical_colors)
            
            # Combine plots with minimal spacing
            combined_plot <- pc1_box / pc2_box + 
              plot_layout(heights = c(1, 1)) +
              plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
            
            return(combined_plot)
          }) 

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod(f="plotDensityList"
          , signature="ProteinQuantitativeData"
          , definition=function(theObject, grouping_variables_list, title = "", font_size = 8) {
            
            # Create a list of density plots for each grouping variable
            density_plots_list <- purrr::map(grouping_variables_list, function(group_var) {
              tryCatch({
                plotDensity(theObject, 
                           grouping_variable = group_var,
                           title = title,
                           font_size = font_size)
              }, error = function(e) {
                warning(sprintf("Error creating density plot for %s: %s", group_var, e$message))
              return(NULL)
                  })
                              })
            
            # Name the list elements with the grouping variables
            names(density_plots_list) <- grouping_variables_list
            
            # Remove any NULL elements (failed plots)
            density_plots_list <- density_plots_list[!sapply(density_plots_list, is.null)]
            
            return(density_plots_list)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
savePlotDensityList <- function(input_list, prefix = "Density", suffix = c("png", "pdf"), output_dir) {
  
  list_of_filenames <- expand_grid(column = names(input_list), suffix = suffix) |>
    mutate(filename = paste0(prefix, "_", column, ".", suffix)) |>
    left_join(tibble(column = names(input_list),
              plots = input_list),
              by = join_by(column))
  
  purrr::walk2(list_of_filenames$plots,
               list_of_filenames$filename,
               \(.x, .y) {
                 ggsave(plot = .x, filename = file.path(output_dir, .y))
               })
  
  list_of_filenames
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
