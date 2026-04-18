#' @export
setMethod(
  f = "plotDensityList",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variables_list, title = "", font_size = 8) {
    # Create a list of density plots for each grouping variable
    density_plots_list <- purrr::map(grouping_variables_list, function(group_var) {
      tryCatch(
        {
          plotDensity(theObject,
            grouping_variable = group_var,
            title = title,
            font_size = font_size
          )
        },
        error = function(e) {
          warning(sprintf("Error creating density plot for %s: %s", group_var, e$message))
          return(NULL)
        }
      )
    })

    # Name the list elements with the grouping variables
    names(density_plots_list) <- grouping_variables_list

    # Remove any NULL elements (failed plots)
    density_plots_list <- density_plots_list[!sapply(density_plots_list, is.null)]

    return(density_plots_list)
  }
)

