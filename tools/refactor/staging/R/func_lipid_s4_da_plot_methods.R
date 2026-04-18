#' plot number of significant differentially expressed lipids
#' @export
setMethod(
    f = "plotNumSigDiffExpBarPlot",
    signature = "list",
    definition = function(objectsList) {
        return_object_list <- purrr::imap(
            objectsList,
            function(object, idx) {
                ## Count the number of up or down significant differentially expressed lipids.
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

#' Plot static volcano plot (without gene names)
#' @export
setMethod(
    f = "plotVolcanoS4",
    signature = "list",
    definition = function(objectsList,
                          da_q_val_thresh = 0.05,
                          qvalue_column = "fdr_qvalue",
                          log2fc_column = "logFC") {
        return_object_list <- purrr::imap(
            objectsList,
            function(object, idx) {
                volcano_colors <- c(
                    "Significant Up" = "red",
                    "Significant Down" = "blue",
                    "Not significant" = "grey"
                )

                ## Plot static volcano plot
                static_volcano_plot_data <- object@contrasts_results_table |>
                    bind_rows(.id = "comparison") |>
                    mutate(lqm = -log10(!!sym(qvalue_column))) |>
                    dplyr::mutate(label = case_when(
                        !!sym(qvalue_column) < da_q_val_thresh & !!sym(log2fc_column) > 0 ~ "Significant Up",
                        !!sym(qvalue_column) < da_q_val_thresh & !!sym(log2fc_column) < 0 ~ "Significant Down",
                        TRUE ~ "Not significant"
                    )) |>
                    # This is the key change: ensure the 'label' column is a factor with all possible levels.
                    # This makes the scales identical across all plots, allowing patchwork to merge them.
                    dplyr::mutate(label = factor(label, levels = names(volcano_colors))) |>
                    dplyr::mutate(colour = case_when(
                        !!sym(qvalue_column) < da_q_val_thresh & !!sym(log2fc_column) < 0 ~ "blue",
                        !!sym(qvalue_column) < da_q_val_thresh & !!sym(log2fc_column) > 0 ~ "red",
                        TRUE ~ "grey"
                    )) |>
                    dplyr::mutate(colour = factor(colour, levels = c("blue", "grey", "red"))) |>
                    dplyr::mutate(title = str_split_i(comparison, "=", 1))

                list_of_volcano_plots_tbl <- static_volcano_plot_data |>
                    group_by(comparison, title) |>
                    nest() |>
                    mutate(plot = purrr::map2(data, title, \(x, y) {
                        plotOneVolcanoNoVerticalLines(x,
                            paste0(idx, " - ", y),
                            log_q_value_column = "lqm",
                            log_fc_column = log2fc_column
                        ) +
                            scale_color_manual(
                                values = volcano_colors,
                                name = "Significance",
                                limits = names(volcano_colors)
                            )
                    }))

                # THE FIX: Extract the 'plot' column to get a LIST of plots,
                # and correctly name the list elements.
                plots_list <- list_of_volcano_plots_tbl$plot
                names(plots_list) <- list_of_volcano_plots_tbl$comparison

                # Assign the LIST of plots to the slot, not the whole table
                object@volcano_plot <- plots_list
                return(object)
            }
        )

        return(return_object_list)
    }
)

