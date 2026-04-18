#' plot number of significant differentially expressed metabolites
#' @export
setMethod(
    f = "plotNumSigDiffExpBarPlot",
    signature = "list",
    definition = function(objectsList) {
        return_object_list <- purrr::imap(
            objectsList,
            function(object, idx) {
                ## Count the number of up or down significant differentially expressed metabolites.
                # The contrasts_results_table is already a list of data frames (one per contrast)
                # So we don't need to wrap it in another list
                num_sig_de_molecules_first_go <- printCountDaGenesTable(
                    list_of_da_tables = object@contrasts_results_table,
                    list_of_descriptions = names(object@contrasts_results_table),
                    formula_string = NA
                )


                object@num_sig_diff_exp_bar_plot <- num_sig_de_molecules_first_go$plot

                object@num_sig_diff_table <- num_sig_de_molecules_first_go$table

                object
            }
        )

        return_object_list
    }
)

