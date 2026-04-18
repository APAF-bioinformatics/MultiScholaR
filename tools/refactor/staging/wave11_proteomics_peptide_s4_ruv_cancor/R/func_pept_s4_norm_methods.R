##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' RUV Canonical Correlation Analysis for Peptide Data
#' 
#' Performs Remove Unwanted Variation (RUV) using canonical correlation analysis
#' with DPC imputation for missing values.
#' 
#' @param theObject A PeptideQuantitativeData object
#' @param ctrl Vector of control gene indices for RUV
#' @param num_components_to_impute Number of principal components for imputation
#' @param ruv_grouping_variable Column name in design matrix for RUV grouping
#' 
#' @export
#' @title RUV Canonical Correlation for PeptideQuantitativeData
#' @name ruvCancor,PeptideQuantitativeData-method
#' @export
setMethod( f = "ruvCancor"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL) {
             
             
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

             # Remove or winsorize extreme values
             peptide_matrix_clean <- apply(log2(peptide_matrix), 2, function(x) {
               q <- quantile(x, probs = c(0.01, 0.99), na.rm=TRUE)
               x[x < q[1]] <- q[1]
               x[x > q[2]] <- q[2]
               return(x)
             })
             
             
             dpcfit <- dpc(peptide_matrix_clean)
             
             peptide_matrix_complete <- dpcImpute(peptide_matrix_clean, dpc=dpcfit)$E
             
             
             Y <-  t( peptide_matrix_complete[,design_matrix |> dplyr::pull(!!sym(sample_id))]) 

             print("steps")
             
             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                dplyr::pull(!!sym(ruv_grouping_variable)),
                                              ctl = ctrl)
             cancorplot_r2
             
             
           })

