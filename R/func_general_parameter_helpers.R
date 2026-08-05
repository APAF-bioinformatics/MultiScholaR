# ----------------------------------------------------------------------------
# setArgsDefault
# ----------------------------------------------------------------------------
#' @export
setArgsDefault <- function(args, value_name, as_func, default_val=NA ) {

  if(isArgumentDefined(args, value_name))
  {
    args<-parseType(args,
                    c(value_name)
                    ,as_func)
  }else {
    logwarn(paste0( value_name, " is undefined, default value set to ", as.character(default_val), "."))
    args[[ value_name ]] <- default_val
  }

  return(args)
}

# ----------------------------------------------------------------------------
# checkParamsObjectFunctionSimplify
# ----------------------------------------------------------------------------
#' Check the parameters in the arguments list and the function parameters to see what param applies
#' @export
checkParamsObjectFunctionSimplify <- function(theObject, param_name_string, default_value = NULL) {

  function_name <- getFunctionNameSecondLevel()

  # print(function_name)
  param_value <- dynGet(param_name_string)

  # Fix: Safely access nested list to avoid index errors
  object_value <- NULL
  if (!is.null(theObject@args) && 
      !is.null(theObject@args[[function_name]]) && 
      is.list(theObject@args[[function_name]])) {
    object_value <- theObject@args[[function_name]][[param_name_string]]
  }

  # print(paste0("param_value = ", param_value))

  error <- paste0(function_name,  paste0(": '", param_name_string, "' is not defined.\n") )

  if( !is.null(param_value) ) {
    return( param_value)
  } else if( !is.null(object_value) ) {
    # print("use object value")
    return( object_value)
  } else if( !is.null(default_value) ) {
    return( default_value)
  } else {
    stop( error )
  }

}

# ----------------------------------------------------------------------------
# checkParamsObjectFunctionSimplifyAcceptNull
# ----------------------------------------------------------------------------
#' Check the parameters in the arguments list and the function parameters to see what param applies
#' @export
checkParamsObjectFunctionSimplifyAcceptNull <- function(theObject, param_name_string, default_value = NULL) {

  function_name <- getFunctionNameSecondLevel()

  # print(function_name)
  param_value <- dynGet(param_name_string)
  object_value <- theObject@args[[function_name]][[param_name_string]]

  # print(paste0("param_value = ", param_value))

  error <- paste0(function_name, ": '", param_name_string, "' is not defined.\n")

  if( !is.null(param_value) ) {
    return( param_value)
  } else if( !is.null(object_value) ) {
    return( object_value)
  } else if ( !is.null(default_value) ) {
    return( default_value)
  } else {
    warning(error)
    return( NULL )
  }

}

# ----------------------------------------------------------------------------
# updateParamInObject
# ----------------------------------------------------------------------------
#' Update the parameter in the object
#'@export
updateParamInObject <- function(theObject, param_name_string) {

  function_name <- getFunctionNameSecondLevel()

  theObject@args[[function_name]][[param_name_string]] <- dynGet(param_name_string)

  theObject

}

# ----------------------------------------------------------------------------
# updateMissingValueParameters
# ----------------------------------------------------------------------------
#' Update Missing Value Parameters in Configuration List and S4 Object
#'
#' @description
#' Automatically calculates and updates the missing value filtering parameters in the configuration list
#' and S4 object @args based on the experimental design matrix. The function ensures at least a specified 
#' number of groups have sufficient quantifiable values for analysis.
#'
#' @param theObject An S4 object containing the experimental data and design matrix.
#' @param min_reps_per_group Integer specifying the minimum number of replicates required in each passing group.
#'                          If a group has fewer total replicates than this value, the minimum is adjusted.
#' @param min_groups Integer specifying the minimum number of groups required to have sufficient
#'                  quantifiable values. Default is 2.
#' @param config_list_name The name of the global config list variable (defaults to "config_list").
#' @param env The environment where the global config list resides (defaults to .GlobalEnv).
#'
#' @return Updated S4 object with synchronized @args and global config_list
#'
#' @details
#' The function calculates:
#' - groupwise_percentage_cutoff: Based on minimum required replicates per group
#' - max_groups_percentage_cutoff: Based on minimum required groups
#' 
#' Both the S4 object's @args slot and the global config_list are updated to maintain synchronization.
#'
#' @examples
#' \dontrun{
#' protein_log2_quant_cln <- updateMissingValueParameters(
#'   theObject = protein_log2_quant_cln, 
#'   min_reps_per_group = 2, 
#'   min_groups = 2
#' )
#' }
#'
#' @export
updateMissingValueParameters <- function(theObject, 
                                       min_reps_per_group = 2, 
                                       min_groups = 2,
                                       function_name = "removeRowsWithMissingValuesPercent",
                                       grouping_variable = "group",
                                       config_list_name = "config_list",
                                       env = .GlobalEnv) {
    
    # --- Input Validation ---
    if (!isS4(theObject)) {
        stop("'theObject' must be an S4 object.")
    }
    if (!"design_matrix" %in% methods::slotNames(theObject)) {
        stop("'theObject' must have a '@design_matrix' slot.")
    }
    if (!"args" %in% methods::slotNames(theObject)) {
        stop("'theObject' must have an '@args' slot.")
    }
    if (!exists(config_list_name, envir = env)) {
        stop("Global config list '", config_list_name, "' not found in the specified environment.")
    }
    
    # Extract design matrix from S4 object
    design_matrix <- theObject@design_matrix
    
    # Retrieve the global config list
    config_list <- get(config_list_name, envir = env)
    
    # Get number of replicates per group
    reps_per_group_tbl <- design_matrix |>
        dplyr::group_by(!!rlang::sym(grouping_variable)) |>
        dplyr::summarise(n_reps = dplyr::n(), .groups = "drop")
    
    # Get total number of groups
    total_groups <- nrow(reps_per_group_tbl)
    
    if (min_groups > total_groups) {
        stop("min_groups cannot be larger than total number of groups")
    }
    
    # Calculate percentage missing allowed for each group
    group_thresholds <- reps_per_group_tbl |>
        dplyr::mutate(
            adjusted_min_reps = pmin(n_reps, min_reps_per_group),
            max_missing = n_reps - adjusted_min_reps,
            missing_percent = round((max_missing / n_reps) * 100, 3)
        )
    
    # Use a consistent percentage threshold across all groups
    # Take the maximum percentage to ensure all groups meet minimum requirements
    groupwise_cutoff <- max(group_thresholds$missing_percent)
    
    # Calculate maximum failing groups allowed
    max_failing_groups <- total_groups - min_groups
    max_groups_cutoff <- round((max_failing_groups / total_groups) * 100, 3)
    
    # Update global config_list
    if (is.null(config_list[[function_name]])) {
        config_list[[function_name]] <- list()
    }
    config_list[[function_name]]$groupwise_percentage_cutoff <- groupwise_cutoff
    config_list[[function_name]]$max_groups_percentage_cutoff <- max_groups_cutoff
    
    # Function specific parameters
    if (function_name == "removeRowsWithMissingValuesPercent") {
        config_list[[function_name]]$proteins_intensity_cutoff_percentile <- 1
    } else if (function_name == "peptideIntensityFiltering") {
        config_list[[function_name]]$peptides_intensity_cutoff_percentile <- 1
    }
    
    # Assign updated config_list back to global environment
    assign(config_list_name, config_list, envir = env)
    
    # Update S4 object's internal args to match
    if (is.null(theObject@args[[function_name]])) {
        theObject@args[[function_name]] <- list()
    }
    theObject@args[[function_name]]$groupwise_percentage_cutoff <- groupwise_cutoff
    theObject@args[[function_name]]$max_groups_percentage_cutoff <- max_groups_cutoff
    
    # Create informative message
    basic_msg <- sprintf(
        "Updated missing value parameters in both config_list and S4 object @args for function '%s':
    - Requiring at least %d replicates in each passing group (varies by group size)
    - Requiring at least %d out of %d groups to pass (%.3f%% failing groups allowed)
    - groupwise_percentage_cutoff set to %.3f%%
    - max_groups_percentage_cutoff set to %.3f%%",
        function_name,
        min_reps_per_group,
        min_groups,
        total_groups,
        max_groups_cutoff,
        groupwise_cutoff,
        max_groups_cutoff
    )
    
    # Add details for each group
    group_detail_strings <- group_thresholds |>
        dplyr::mutate(
            detail = sprintf("    Group %s: %d out of %d replicates required (%.3f%% missing allowed)",
                             !!rlang::sym(grouping_variable), adjusted_min_reps, n_reps, missing_percent)
        ) |>
        dplyr::pull(detail)
        
    group_details <- paste(group_detail_strings, collapse = "\n")
    
    # Print the message
    message(paste(basic_msg, "\n\nGroup details:", group_details, sep = "\n"))
    message("[OK] S4 object @args and global config_list are now synchronized")
    
    return(theObject)
}

