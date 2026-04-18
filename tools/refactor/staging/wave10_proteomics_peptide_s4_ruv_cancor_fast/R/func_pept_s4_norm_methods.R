#' @title Fast version of ruvCancor for optimization
#' @name ruvCancorFast,PeptideQuantitativeData-method
#' @description This function provides a lightweight alternative to ruvCancor that skips
#' the expensive DPC imputation step, making it suitable for use during
#' optimization loops where speed is critical.
#' @param theObject A PeptideQuantitativeData object
#' @param ctrl Control features for RUV
#' @param num_components_to_impute Number of components to impute
#' @param ruv_grouping_variable Grouping variable for RUV
#' @param simple_imputation_method Method for simple missing value handling.
#'   Options: "none" (default), "mean", "median", "min"
#' @export
setMethod( f = "ruvCancorFast"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ctrl= NULL, num_components_to_impute=NULL, 
                                  ruv_grouping_variable = NULL, simple_imputation_method = "none") {
             
             peptide_matrix <- theObject@peptide_matrix
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
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
             
             # Skip expensive DPC imputation - use simple approach for optimization
             peptide_matrix_working <- log2(peptide_matrix + 1)  # Add small constant to avoid log(0)
             
             # Handle missing values with simple method if requested
             if (simple_imputation_method != "none" && anyNA(peptide_matrix_working)) {
               if (simple_imputation_method == "mean") {
                 peptide_matrix_working <- apply(peptide_matrix_working, 2, function(x) {
                   x[is.na(x)] <- mean(x, na.rm = TRUE)
                   return(x)
                 })
               } else if (simple_imputation_method == "median") {
                 peptide_matrix_working <- apply(peptide_matrix_working, 2, function(x) {
                   x[is.na(x)] <- median(x, na.rm = TRUE)
                   return(x)
                 })
               } else if (simple_imputation_method == "min") {
                 peptide_matrix_working <- apply(peptide_matrix_working, 2, function(x) {
                   x[is.na(x)] <- min(x, na.rm = TRUE)
                   return(x)
                 })
               }
             }
             
             # Remove or winsorize extreme values (lightweight version)
             peptide_matrix_clean <- apply(peptide_matrix_working, 2, function(x) {
               q <- quantile(x, probs = c(0.01, 0.99), na.rm=TRUE)
               x[x < q[1]] <- q[1]
               x[x > q[2]] <- q[2]
               return(x)
             })
             
             # Use the matrix directly without expensive DPC imputation
             Y <- t(peptide_matrix_clean[, design_matrix |> dplyr::pull(!!sym(sample_id))])
             
             # Generate canonical correlation plot (this is what we actually need for optimization)
             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                dplyr::pull(!!sym(ruv_grouping_variable)),
                                              ctl = ctrl)
             
             return(cancorplot_r2)
           })

