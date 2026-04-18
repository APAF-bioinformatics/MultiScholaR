# ----------------------------------------------------------------------------
# updateRuvParameters
# ----------------------------------------------------------------------------
##################################################################################################################
#' @export
updateRuvParameters <- function(config_list, best_k, control_genes_index, percentage_as_neg_ctrl) {
  config_list$ruvParameters$best_k <- best_k
  config_list$ruvParameters$num_neg_ctrl <- length(control_genes_index)
  config_list$ruvParameters$percentage_as_neg_ctrl <- percentage_as_neg_ctrl
  
  # Print the number of negative controls (as in the original code)
  config_list$ruvParameters$num_neg_ctrl
  
  # Return the updated config list
  return(config_list)
}

# ----------------------------------------------------------------------------
# getRuvIIIReplicateMatrixHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'Converts a design matrix to a biological replicate matrix for use with ruvIII.
#'@param design_matrix The design matrix with the sample ID in one column and the experimental group in another column
#'@param sample_id_column The name of the column with the sample ID, tidyverse style input.
#'@param grouping_variable The name of the column with the experimental group, tidyverse style input.
#'@param temp_column The name of the temporary column that indicates which samples are biological replicates of the same experimental group.
#'@return A numeric matrix with rows as samples, columns as experimental group, and a value of 1 for samples within the same experimental group represented by the same column, and a value of zero otherwise.
#'@export
getRuvIIIReplicateMatrixHelper <- function(design_matrix, sample_id_column, grouping_variable, temp_column = is_replicate_temp) {

  ruvIII_replicates_matrix <- design_matrix |>
    dplyr::select({ { sample_id_column } }, { { grouping_variable } }) |>
    mutate({ { temp_column } } := 1) |>
    pivot_wider(id_cols = as_string( as_name(enquo(sample_id_column ))),
                names_from = { { grouping_variable } },
                values_from = { { temp_column } },
                values_fill = 0) |>
    column_to_rownames(as_string(as_name(enquo(sample_id_column)))) |>
    as.matrix()

  ruvIII_replicates_matrix
}

# ----------------------------------------------------------------------------
# getNegCtrlProtAnovaHelper
# ----------------------------------------------------------------------------
#' Identify negative control proteins for use in removal of unwanted variation, using an ANOVA test.
#' @param data_matrix A matrix containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The row ID are the protein accessions. The data is preferably median-scaled with missing values imputed.
#' @param design_matrix A data frame with the design matrix. Matches sample IDs to group IDs.
#' @param grouping_variable The name of the column with the experimental group, as a string.
#' @param num_neg_ctrl The number of negative control genes to select. Typically the number of genes with the highest q-value (e.g. least statistically significant). Default is 100
#' @param ruv_qval_cutoff The FDR threshold. No proteins with q-values lower than this value are included in the list of negative control proteins. This means the number of negative control proteins could be less than the number specified in \code{num_neg_ctrl} when some were excluded by this threshold.
#' @param ruv_fdr_method The FDR calculation method, default is "qvalue". The other option is "BH"
#' @param exclude_pool_samples Logical. If TRUE (default), automatically exclude samples from groups containing "Pool" in their name from ANOVA calculation. Pool/QC samples are excluded from negative control selection but remain in RUV-III correction.
#' @return A boolean vector which indicates which row in the input data matrix is a control gene. The row is included if the value is TRUE. The names of each element is the row ID / protein accessions of the input data matrix.
#'@export
getNegCtrlProtAnovaHelper <- function(data_matrix
                                , design_matrix
                                , grouping_variable = "group"
                                , percentage_as_neg_ctrl = 10
                                , num_neg_ctrl = round( nrow(data_matrix)*percentage_as_neg_ctrl/100, 0)
                                , ruv_qval_cutoff = 0.05
                                , ruv_fdr_method = "qvalue"
                                , exclude_pool_samples = TRUE) {

  message("--- DEBUG66: Entering getNegCtrlProtAnovaHelper ---")
  message(sprintf("   DEBUG66 Arg: grouping_variable = %s", grouping_variable))
  message(sprintf("   DEBUG66 Arg: percentage_as_neg_ctrl = %s", ifelse(is.null(percentage_as_neg_ctrl), "NULL", percentage_as_neg_ctrl)))
  message(sprintf("   DEBUG66 Arg: num_neg_ctrl = %s", ifelse(is.null(num_neg_ctrl), "NULL", num_neg_ctrl)))
  message(sprintf("   DEBUG66 Arg: exclude_pool_samples = %s", exclude_pool_samples))
  message(sprintf("   DEBUG66 Data State: data_matrix dims = %d x %d", nrow(data_matrix), ncol(data_matrix)))
  message(sprintf("   DEBUG66 Data State: design_matrix dims = %d x %d", nrow(design_matrix), ncol(design_matrix)))

  ## Both percentage_as_neg_ctrl and num_neg_ctrl is missing, and number of proteins >= 50 use only 10 percent of the proteins as negative control by default
  if((is.null(percentage_as_neg_ctrl) ||
     is.na(percentage_as_neg_ctrl) ) &&
     (is.null(num_neg_ctrl) ||
      is.na(num_neg_ctrl) ) &&
     nrow(data_matrix) >= 50 ) {
    num_neg_ctrl <- round( nrow(data_matrix)*10/100, 0)
    warnings( paste0( getFunctionName(), ": Using 10% of proteins from the input matrix as negative controls by default.\n"))
    message("   DEBUG66: Defaulting to 10% negative controls")
  } else if (!is.null(percentage_as_neg_ctrl) &
             !is.na(percentage_as_neg_ctrl)) {
    num_neg_ctrl <- round( nrow(data_matrix)*percentage_as_neg_ctrl/100, 0)
  } else if(!is.null(num_neg_ctrl) &
       !is.na(num_neg_ctrl)) {
    num_neg_ctrl <- as.integer(num_neg_ctrl)
  } else {
    stop(paste0( getFunctionName(), ": Please provide either percentage_as_neg_ctrl or num_neg_ctrl.\n"))
  }
  
  message(sprintf("   DEBUG66: Target num_neg_ctrl = %d", num_neg_ctrl))

  # --- Pool/QC Sample Handling ---
  # Pool/QC samples should not influence negative control selection but can help
  # estimate unwanted variation patterns in RUV-III correction
  data_matrix_for_anova <- data_matrix
  design_matrix_for_anova <- design_matrix
  
  if (exclude_pool_samples) {
    message("   DEBUG66: Checking for Pool/QC samples to exclude...")
    # Detect Pool/QC groups (case-insensitive detection)
    all_groups <- unique(design_matrix[[grouping_variable]])
    is_pool_group <- grepl("pool", all_groups, ignore.case = TRUE)
    pool_group_names <- all_groups[is_pool_group]
    
    if (length(pool_group_names) > 0) {
      message(sprintf("   DEBUG66: Found pool groups: %s", paste(pool_group_names, collapse=", ")))
      # Identify samples that belong to Pool groups
      samples_in_design <- rownames(design_matrix)
      pool_samples <- samples_in_design[design_matrix[[grouping_variable]] %in% pool_group_names]
      
      # Filter to samples present in both design matrix and data matrix columns
      pool_samples_in_data <- intersect(pool_samples, colnames(data_matrix))
      
      if (length(pool_samples_in_data) > 0) {
        # Log Pool/QC exclusion details
        message(sprintf("*** NEG CTRL ANOVA: Detected %d Pool/QC group(s): %s ***", 
                       length(pool_group_names), paste(pool_group_names, collapse = ", ")))
        message(sprintf("*** NEG CTRL ANOVA: Excluding %d Pool samples from ANOVA ***", 
                       length(pool_samples_in_data)))
        message(sprintf("   DEBUG66: Pool samples: %s", paste(pool_samples_in_data, collapse=", ")))
        
        # Create filtered matrices excluding Pool samples
        non_pool_samples <- setdiff(colnames(data_matrix), pool_samples_in_data)
        data_matrix_for_anova <- data_matrix[, non_pool_samples, drop = FALSE]
        design_matrix_for_anova <- design_matrix[non_pool_samples, , drop = FALSE]
        
        message(sprintf("*** NEG CTRL ANOVA: Using %d non-Pool samples for negative control selection ***", 
                       ncol(data_matrix_for_anova)))
        
        # Validate sufficient samples remain for ANOVA
        remaining_groups <- unique(design_matrix_for_anova[[grouping_variable]])
        if (length(remaining_groups) < 2) {
          stop(paste("*** NEG CTRL ANOVA ERROR: After excluding Pool samples, only", 
                    length(remaining_groups), 
                    "group(s) remain. Need at least 2 groups for ANOVA. ***"))
        }
        
        # Check samples per group
        samples_per_group <- table(design_matrix_for_anova[[grouping_variable]])
        groups_with_single_sample <- sum(samples_per_group < 2)
        if (groups_with_single_sample > 0) {
          message(sprintf("*** NEG CTRL ANOVA WARNING: %d group(s) have only 1 sample after Pool exclusion ***", 
                         groups_with_single_sample))
        }
      } else {
        message("*** NEG CTRL ANOVA: Pool groups detected but no Pool samples found in data matrix ***")
      }
    } else {
      message("*** NEG CTRL ANOVA: No Pool/QC groups detected - using all samples for ANOVA ***")
    }
  }

  ## Inspired by matANOVA function from PhosR package: http://www.bioconductor.org/packages/release/bioc/html/PhosR.html

  message("   DEBUG66: Preparing for ANOVA...")
  grps <- design_matrix_for_anova[colnames(data_matrix_for_anova), grouping_variable]
  message(sprintf("   DEBUG66: ANOVA groups vector length: %d", length(grps)))

  ps <- apply(data_matrix_for_anova, 1, function(x) {
       if( length( unique( grps[!is.na(x)] )  ) > 1 ) {
         summary(stats::aov(as.numeric(x) ~ grps))[[1]][["Pr(>F)"]][1]
       } else {
          return(NA_real_)
       }
    })

  message(sprintf("   DEBUG66: ANOVA p-values calculated. NA count: %d", sum(is.na(ps))))
  ps[is.na(ps)] <- 1

  aov <- c()

  if ( ruv_fdr_method == "qvalue") {
    message("   DEBUG66: Calculating q-values...")
    aov <- qvalue(unlist(ps))$qvalues
  } else if ( ruv_fdr_method == "BH") {
    message("   DEBUG66: Calculating BH adjusted p-values...")
    aov <- qvalue(unlist(ps), pi0=1)$qvalues
  } else {
    error( paste( "Input FDR method", ruv_fdr_method, "not valid") )
  }

  message(sprintf("   DEBUG66: q-values calculated. Range: %.4f - %.4f", min(aov, na.rm=TRUE), max(aov, na.rm=TRUE)))

  filtered_list <- aov[aov > ruv_qval_cutoff]
  message(sprintf("   DEBUG66: Proteins with q > %.4f: %d", ruv_qval_cutoff, length(filtered_list)))

  list_size <- ifelse(num_neg_ctrl > length(filtered_list), length(filtered_list), num_neg_ctrl)
  message(sprintf("   DEBUG66: Selecting top %d proteins", list_size))

  control_genes <- names(sort(filtered_list, decreasing = TRUE)[1:list_size])

  #nrow(data_matrix) - length(control_genes)
  control_genes_index <- rownames(data_matrix) %in% control_genes
  names(control_genes_index) <- rownames(data_matrix)

  message(sprintf("   DEBUG66: Final control genes selected: %d", sum(control_genes_index)))
  message("--- DEBUG66: Exiting getNegCtrlProtAnovaHelper ---")

  return(control_genes_index)

}

# ----------------------------------------------------------------------------
# extractRuvResults
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
extractRuvResults <- function(results_list) {

  extracted <- purrr::map(results_list, \(x){ x$results })

  names(extracted) <- names(results_list)

  return(extracted)
}

# ----------------------------------------------------------------------------
# scaleCenterAndFillMissing
# ----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#' @export
#' @title Scale, Center and Fill Missing Values
#'@param input_matrix Samples are columns, rows are proteins or peptides
scaleCenterAndFillMissing <- function( input_matrix) {

  input_matrix_scaled <- scale(input_matrix, center = TRUE, scale = TRUE)


  min_data_point <- min(input_matrix_scaled, na.rm=TRUE)


  input_matrix_scaled_fill_missing <-  input_matrix_scaled
  input_matrix_scaled_fill_missing[is.na(input_matrix_scaled_fill_missing)] <- min_data_point*2

  input_matrix_scaled_fill_missing
}

