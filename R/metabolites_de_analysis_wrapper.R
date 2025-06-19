#' MetabolomicsDifferentialAbundanceResults S4 Class
#'
#' @description
#' S4 class to store essential results from metabolomics differential abundance analysis.
#' This class contains the original data object, fitted model, and results table.
#'
#' @slot theObject The original MetaboliteAssayData object used for analysis
#' @slot fit.eb The fitted eBayes model from limma analysis
#' @slot contrasts_results_table Data frame with differential abundance statistics
#'
#' @export
setClass("MetabolomicsDifferentialAbundanceResults"
         , slots = c(
             theObject = "MetaboliteAssayData"
           , fit.eb = "MArrayLM"
           , contrasts_results_table = "list"
           , num_sig_diff_exp_bar_plot = "list"
           , num_sig_diff_table = "data.frame"
           , volcano_plot = "list"
           , interactive_volcano_plot = "list"
           , p_value_dist_plot = "list"
           , results_table_long = "data.frame"
           , results_table_wide = "data.frame"
         ),
         prototype = list(
             theObject = NULL
           , fit.eb = NULL
           , contrasts_results_table = list()
           , num_sig_diff_exp_bar_plot =  list()
           , num_sig_diff_table =  data.frame()
           , volcano_plot =  list()
           , interactive_volcano_plot = list()
           , p_value_dist_plot =  list()
           , results_table_long = data.frame()
           , results_table_wide = data.frame()
         )
)

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric( name ="differentialAbundanceAnalysis"
            , def=function(objectsList
                           , contrasts_tbl = NULL
                           , formula_string = NULL
                           , de_q_val_thresh = NULL
                           , treat_lfc_cutoff = NULL
                           , eBayes_trend = NULL
                           , eBayes_robust = NULL
                           , args_group_pattern = NULL) {
              standardGeneric("differentialAbundanceAnalysis")})



#'@export
setMethod( f ="differentialAbundanceAnalysis"
           , signature = "list"
           , definition=function( objectsList
                                  , contrasts_tbl = NULL
                                  , formula_string = NULL
                                  , de_q_val_thresh = NULL
                                  , treat_lfc_cutoff = NULL
                                  , eBayes_trend = NULL
                                  , eBayes_robust = NULL
                                  , args_group_pattern = NULL ) {

             # Validate that all objects in the list are MetaboliteAssayData
             if (!all(purrr::map_lgl(objectsList, ~inherits(.x, "MetaboliteAssayData")))) {
               stop("All objects in objectsList must be of class MetaboliteAssayData")
             }

             # Run DE analysis and explicitly set names
             results_list <- purrr::map(    objectsList
                                            , \( obj) {
                                              differentialAbundanceAnalysisHelper(  obj
                                                                                    , contrasts_tbl = contrasts_tbl
                                                                                    , formula_string = formula_string
                                                                                    , de_q_val_thresh = de_q_val_thresh
                                                                                    , treat_lfc_cutoff = treat_lfc_cutoff
                                                                                    , eBayes_trend = eBayes_trend
                                                                                    , eBayes_robust = eBayes_robust
                                                                                    , args_group_pattern = args_group_pattern
                                              )
                                            })

             # Set names if the input list had names
             if (!is.null(names(objectsList))) {
               names(results_list) <- names(objectsList)
             }

             return(results_list)

           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
setGeneric( name ="differentialAbundanceAnalysisHelper"
            , def=function(theObject
                           , contrasts_tbl = NULL
                           , formula_string = NULL
                           , de_q_val_thresh = NULL
                           , treat_lfc_cutoff = NULL
                           , eBayes_trend = NULL
                           , eBayes_robust = NULL
                           , args_group_pattern = NULL) {
              standardGeneric("differentialAbundanceAnalysisHelper")
            })

#'@export
setMethod( f ="differentialAbundanceAnalysisHelper"
           , signature = "MetaboliteAssayData"
           , definition=function( theObject
                                  , contrasts_tbl = NULL
                                  , formula_string = NULL
                                  , de_q_val_thresh = NULL
                                  , treat_lfc_cutoff = NULL
                                  , eBayes_trend = NULL
                                  , eBayes_robust = NULL
                                  , args_group_pattern = NULL) {

  message("--- Entering differentialAbundanceAnalysisHelper ---")

  contrasts_tbl <- checkParamsObjectFunctionSimplify( theObject, "contrasts_tbl", NULL)
  formula_string <- checkParamsObjectFunctionSimplify( theObject, "formula_string", " ~ 0 + group")
  de_q_val_thresh <- checkParamsObjectFunctionSimplify( theObject, "de_q_val_thresh", 0.05)
  treat_lfc_cutoff <- checkParamsObjectFunctionSimplify( theObject, "treat_lfc_cutoff", 0)
  eBayes_trend <- checkParamsObjectFunctionSimplify( theObject, "eBayes_trend", TRUE)
  eBayes_robust <- checkParamsObjectFunctionSimplify( theObject, "eBayes_robust", TRUE)
  args_group_pattern <- checkParamsObjectFunctionSimplify( theObject, "args_group_pattern", "(\\d+)")

  #message("   differentialAbundanceAnalysisHelper Arg: contrasts_tbl = ")
  #print(utils::str(contrasts_tbl))

  print(formula_string)

  # Add preprocessing for group names that start with numbers
  design_matrix <- theObject@design_matrix
  group_col <- design_matrix[["group"]]

  # Check if any group names start with numbers and create mapping
  starts_with_number <- grepl("^[0-9]", group_col)
  if(any(starts_with_number)) {
    original_groups <- unique(group_col)
    safe_groups <- purrr::map_chr(original_groups, \(x) {
      if(grepl("^[0-9]", x)) paste0("grp_", x) else x
    })
    group_mapping <- setNames(original_groups, safe_groups)

    # Update design matrix with safe names
    design_matrix[["group"]] <- purrr::map_chr(group_col, \(x) {
      if(grepl("^[0-9]", x)) paste0("grp_", x) else x
    })

    # Update contrasts table if it exists
    if(!is.null(contrasts_tbl)) {
      contrasts_tbl[[1]] <- purrr::map_chr(contrasts_tbl[[1]], \(x) {
        for(orig in names(group_mapping)) {
          x <- gsub(group_mapping[orig], orig, x, fixed = TRUE)
        }
        x
      })
    }

    theObject@design_matrix <- design_matrix
  }

  theObject <- updateParamInObject(theObject, "contrasts_tbl")
  theObject <- updateParamInObject(theObject, "formula_string")
  theObject <- updateParamInObject(theObject, "de_q_val_thresh")
  theObject <- updateParamInObject(theObject, "treat_lfc_cutoff")
  theObject <- updateParamInObject(theObject, "eBayes_trend")
  theObject <- updateParamInObject(theObject, "eBayes_robust")
  theObject <- updateParamInObject(theObject, "args_group_pattern")

  return_list <- list()
  return_list$theObject <- theObject

  message("--- Entering differentialAbundanceAnalysis ---")
  message(sprintf("   differentialAbundanceAnalysis: theObject class = %s", class(theObject)))

  # Debug: Print parameter values before conversion
  message(sprintf("   eBayes_trend value: %s (class: %s)", eBayes_trend, class(eBayes_trend)))
  message(sprintf("   eBayes_robust value: %s (class: %s)", eBayes_robust, class(eBayes_robust)))

  # Ensure logical values - handle character strings that might represent logical values
  if (is.character(eBayes_trend)) {
    eBayes_trend <- as.logical(toupper(eBayes_trend) %in% c("TRUE", "T", "1", "YES"))
    message(sprintf("   eBayes_trend converted from character '%s' to logical: %s", eBayes_trend, eBayes_trend))
  } else {
    eBayes_trend <- as.logical(eBayes_trend)
  }

  if (is.character(eBayes_robust)) {
    eBayes_robust <- as.logical(toupper(eBayes_robust) %in% c("TRUE", "T", "1", "YES"))
    message(sprintf("   eBayes_robust converted from character '%s' to logical: %s", eBayes_robust, eBayes_robust))
  } else {
    eBayes_robust <- as.logical(eBayes_robust)
  }

  message(sprintf("   eBayes_trend after conversion: %s (class: %s)", eBayes_trend, class(eBayes_trend)))
  message(sprintf("   eBayes_robust after conversion: %s (class: %s)", eBayes_robust, class(eBayes_robust)))

  ## Compare the different experimental groups and obtain lists of differentially expressed metabolites

  rownames( theObject@design_matrix ) <- theObject@design_matrix |> dplyr::pull( one_of(theObject@sample_id ))

  # Prepare data matrix for DE analysis
  data_matrix <- NA

  matrix_data <- as.matrix(theObject@metabolite_data[, -1]) # Exclude Name column
  colnames(matrix_data) <- colnames(theObject@metabolite_data)[-1]
  rownames(matrix_data) <- theObject@metabolite_data$Name
  data_matrix <- matrix_data
  message("   differentialAbundanceAnalysisHelper Step: Calling runTestsContrasts...")
  contrasts_results <- runTestsContrasts(
    data_matrix,
    contrast_strings = contrasts_tbl$contrasts,
    design_matrix = theObject@design_matrix,
    formula_string = formula_string,
    treat_lfc_cutoff = treat_lfc_cutoff,
    eBayes_trend = eBayes_trend,
    eBayes_robust = eBayes_robust
  )
  message("   differentialAbundanceAnalysisHelper Step: runTestsContrasts completed.")

 #  # Combine all contrast results into a single data frame
 # # message("   Data State (contrasts_results$results) Structure:")
 # # utils::str(contrasts_results$results)
 # # message("   Data State (contrasts_results$results) Head:")
 # # print(head(contrasts_results$results))
 #  contrasts_results_table <-  dplyr::bind_rows(contrasts_results$results, .id = "comparison")

  # Map back to original group names in results if needed
  if(exists("group_mapping")) {
    contrasts_results$results <- contrasts_results$results |>
      purrr::map( \(results_table){
        results_table |>
          dplyr::mutate(comparison = purrr::map_chr(comparison, \(x) {
            result <- x
            for(safe_name in names(group_mapping)) {
              result <- gsub(safe_name, group_mapping[safe_name], result, fixed = TRUE)
            }
            result
          }))

      })
  }


  return_list$fit.eb <- contrasts_results$fit.eb
  return_list$contrasts_results_table <- contrasts_results$results |>
    purrr::map( \(result_table) {
      result_table |>
        rownames_to_column(var = theObject@metabolite_id_column)
    })

  # Create and return the S4 object
  result_object <- new("MetabolomicsDifferentialAbundanceResults",
                       theObject = return_list$theObject,
                       fit.eb = return_list$fit.eb,
                       contrasts_results_table = return_list$contrasts_results_table
  )
  message("--- Exiting differentialAbundanceAnalysisHelper ---")
  return(result_object)
})


# Helper function to get counts table
getCountsTable <- function(obj) {
  if (inherits(obj, "MetaboliteAssayData")) {
    message(sprintf("   Getting counts table for object of class: %s", class(obj)[1]))
    message(sprintf("   Returning metabolite_data with dimensions: %d rows, %d cols",
                    nrow(obj@metabolite_data), ncol(obj@metabolite_data)))
    obj@metabolite_data
  } else if (inherits(obj, "ProteinQuantitativeData")) {
    message(sprintf("   Returning protein_quant_table with dimensions: %d rows, %d cols",
                    nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)))
    obj@protein_quant_table
  } else {
    message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
    stop("Unsupported object type")
  }
}
