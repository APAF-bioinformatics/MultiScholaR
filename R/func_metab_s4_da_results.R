# Helper function to get counts table
getCountsTable <- function(obj) {
    if (inherits(obj, "MetaboliteAssayData")) {
        message(sprintf("   Getting counts table for object of class: %s", class(obj)[1]))
        message(sprintf(
            "   Returning metabolite_data with dimensions: %d rows, %d cols",
            nrow(obj@metabolite_data), ncol(obj@metabolite_data)
        ))
        obj@metabolite_data
    } else if (inherits(obj, "ProteinQuantitativeData")) {
        message(sprintf(
            "   Returning protein_quant_table with dimensions: %d rows, %d cols",
            nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)
        ))
        obj@protein_quant_table
    } else {
        message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
        stop("Unsupported object type")
    }
}

# ==========================================
# Content from metabolite_da_analysis_wrapper.R
# ==========================================
#' MetabolomicsDifferentialAbundanceResults S4 Class
#'
#' @description
#' S4 class to store essential results from metabolomics differential abundance analysis.
#' This class contains the original data object, fitted model, and results table.
#'
#' @slot theObject The original MetaboliteAssayData object used for analysis
#' @slot fit.eb The fitted eBayes model from limma analysis (ANY type to allow limma MArrayLM when available)
#' @slot contrasts_results_table Data frame with differential abundance statistics
#'
#' @export
setClass("MetabolomicsDifferentialAbundanceResults",
    slots = c(
        theObject = "MetaboliteAssayData",
        fit.eb = "ANY",
        contrasts_results_table = "list",
        num_sig_diff_exp_bar_plot = "list",
        num_sig_diff_table = "data.frame",
        volcano_plot = "list",
        interactive_volcano_plot = "list",
        p_value_dist_plot = "list",
        results_table_long = "data.frame",
        results_table_wide = "data.frame"
    ),
    prototype = list(
        theObject = NULL,
        fit.eb = NULL,
        contrasts_results_table = list(),
        num_sig_diff_exp_bar_plot = list(),
        num_sig_diff_table = data.frame(),
        volcano_plot = list(),
        interactive_volcano_plot = list(),
        p_value_dist_plot = list(),
        results_table_long = data.frame(),
        results_table_wide = data.frame()
    )
)

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

#' Plot static volcano plot (without gene names)
#' @export
setMethod(
    f = "plotVolcanoS4",
    signature = "list",
    definition = function(
      objectsList,
      da_q_val_thresh = 0.05,
      qvalue_column = "fdr_qvalue",
      log2fc_column = "logFC"
    ) {
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

#' @export
setMethod(
    f = "getDaResultsWideFormat",
    signature = "list",
    definition = function(
      objectsList,
      qvalue_column = "fdr_qvalue",
      raw_pvalue_column = "raw_pvalue",
      log2fc_column = "logFC"
    ) {
        return_object_list <- purrr::map(objectsList, function(object) {
            # Correctly access the metabolite data from the nested 'theObject' slot.
            # This defensively handles cases where the slot might hold a list of
            # data frames (correct) or a single data frame (incorrect but handled).
            counts_data_slot <- object@theObject@metabolite_data
            counts_table_to_use <- if (is.list(counts_data_slot) && !is.data.frame(counts_data_slot)) {
                counts_data_slot[[1]]
            } else {
                counts_data_slot
            }

            id_col_name <- object@theObject@metabolite_id_column

            # Bind the list of data frames into a single tidy data frame
            tidy_results <- object@contrasts_results_table |>
                dplyr::bind_rows(.id = "comparison") |>
                dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

            # Pivot the tidy data frame to a wide format using the correct ID column
            wide_results <- tidy_results |>
                tidyr::pivot_wider(
                    id_cols = c(!!sym(id_col_name)),
                    names_from = c(comparision_short),
                    names_sep = ":",
                    values_from = c(!!sym(log2fc_column), !!sym(qvalue_column), !!sym(raw_pvalue_column))
                ) |>
                # Join with original counts using the correct ID column.
                dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
                dplyr::arrange(dplyr::across(matches(qvalue_column))) |>
                dplyr::distinct()

            # Assign to the correct slot and return the object
            object@results_table_wide <- wide_results
            return(object)
        })

        return(return_object_list)
    }
)

# Get the differential expression results in wide format
#' @export
setMethod(
    f = "getDaResultsLongFormat",
    signature = "list",
    definition = function(objectsList) {
        return_object_list <- purrr::map(objectsList, function(object) {
            # Correctly access the metabolite data from the nested 'theObject' slot.
            # This defensively handles cases where the slot might hold a list of
            # data frames (correct) or a single data frame (incorrect but handled).
            counts_data_slot <- object@theObject@metabolite_data
            counts_table_to_use <- if (is.list(counts_data_slot) && !is.data.frame(counts_data_slot)) {
                counts_data_slot[[1]]
            } else {
                counts_data_slot
            }

            id_col_name <- object@theObject@metabolite_id_column

            # Bind the list of data frames into a single tidy data frame
            tidy_results <- object@contrasts_results_table |>
                dplyr::bind_rows(.id = "comparison") |>
                dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

            # Pivot the tidy data frame to a wide format using the correct ID column
            long_results <- tidy_results |>
                # Join with original counts using the correct ID column.
                dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
                dplyr::distinct()

            print(head(long_results))

            # Assign to the correct slot and return the object
            object@results_table_long <- long_results
            return(object)
        })

        return(return_object_list)
    }
)

## Create proteomics interactive volcano plot
#' @export
setMethod(
    f = "plotInteractiveVolcano",
    signature = "list",
    definition =
        function(objectsList, anno_list = NULL) {
            list_of_objects <- purrr::imap(
                objectsList,
                \(de_output_object, idx) {
                    # This helper function now correctly uses its own argument 'r_obj'
                    # and includes error handling for the qvalue calculation.
                    updateWithQvalue <- function(r_obj) {
                        r_obj_output <- r_obj
                        if (is.null(r_obj$p.value)) {
                            return(r_obj_output)
                        }

                        for (coef in seq_len(ncol(r_obj$p.value))) {
                            p_values <- r_obj$p.value[, coef]
                            if (any(is.na(p_values))) {
                                warning(sprintf("NA p-values found in coefficient %d. These will result in NA q-values.", coef), call. = FALSE)
                            }
                            q_result <- tryCatch(
                                {
                                    qvalue::qvalue(p_values)$qvalues
                                },
                                error = function(e) {
                                    warning(sprintf("qvalue calculation failed for coefficient %d: %s. Returning original p-values as q-values.", coef, e$message), call. = FALSE)
                                    p_values # Fallback to p-values on error
                                }
                            )
                            r_obj_output$p.value[, coef] <- q_result
                        }
                        r_obj_output
                    }

                    my_fit_eb <- updateWithQvalue(de_output_object@fit.eb)

                    counts_matrix <- de_output_object@theObject@metabolite_data[[1]] |>
                        column_to_rownames(de_output_object@theObject@metabolite_id_column) |>
                        as.matrix()

                    # Defensive measure: ensure my_fit_eb components have rownames.
                    # `MArrayLM` objects don't have a `rownames<-` method, so we edit the components.
                    if (!is.null(my_fit_eb$coefficients) && is.null(rownames(my_fit_eb$coefficients))) {
                        if (nrow(my_fit_eb$coefficients) == nrow(counts_matrix)) {
                            warning("`my_fit_eb` components were missing rownames. Assigning them from `counts_matrix`.", call. = FALSE)
                            feature_names <- rownames(counts_matrix)
                            rownames(my_fit_eb$coefficients) <- feature_names
                            rownames(my_fit_eb$p.value) <- feature_names
                            if (!is.null(my_fit_eb$t)) {
                                rownames(my_fit_eb$t) <- feature_names
                            }
                            if (!is.null(my_fit_eb$stdev.unscaled)) {
                                rownames(my_fit_eb$stdev.unscaled) <- feature_names
                            }
                            if (!is.null(my_fit_eb$genes)) {
                                rownames(my_fit_eb$genes) <- feature_names
                            }
                        } else {
                            stop("`my_fit_eb` components are missing rownames, and their row count does not match `counts_matrix`. Cannot proceed.")
                        }
                    }

                    groups <- data.frame(Run = colnames(counts_matrix)) |>
                        left_join(de_output_object@theObject@design_matrix,
                            by = join_by(!!sym(de_output_object@theObject@sample_id) == !!sym(de_output_object@theObject@sample_id))
                        ) |>
                        dplyr::pull(genotype_group)


                    list_of_glimma_objs <- purrr::map(
                        seq_len(ncol(de_output_object@fit.eb$p.value)),
                        \(idxb) {
                            id_col_name <- de_output_object@theObject@metabolite_id_column

                            # More robust way to create the base annotation table
                            base_anno_tbl <- data.frame(id = rownames(my_fit_eb))
                            colnames(base_anno_tbl) <- id_col_name


                            # If a user-provided anno_tbl exists, join it.
                            anno_tbl_joined <- if (!is.null(anno_list)) {
                                print(paste("idx=", idx))
                                anno_tbl <- anno_list[[idx]]
                                # Defensive check: ensure anno_tbl is a data.frame and has the join column
                                if (!is.data.frame(anno_tbl) || !id_col_name %in% colnames(anno_tbl)) {
                                    warning(sprintf("Provided 'anno_tbl' is not a data.frame or is missing the join column '%s'. Ignoring 'anno_tbl'.", id_col_name))
                                    base_anno_tbl
                                } else {
                                    dplyr::left_join(
                                        base_anno_tbl |>
                                            mutate(!!sym(id_col_name) := purrr::map_chr(!!sym(id_col_name), as.character)),
                                        anno_tbl |>
                                            mutate(!!sym(id_col_name) := purrr::map_chr(!!sym(id_col_name), as.character)),
                                        by = id_col_name
                                    )
                                }
                            } else {
                                base_anno_tbl
                            }

                            # Glimma requires a data.frame with rownames set to the feature IDs.
                            anno_df_for_glimma <- as.data.frame(anno_tbl_joined)
                            rownames(anno_df_for_glimma) <- anno_df_for_glimma[[id_col_name]]


                            Glimma::glimmaVolcano(my_fit_eb,
                                coef = idxb,
                                anno = anno_df_for_glimma,
                                counts = counts_matrix,
                                groups = groups,
                                display.columns = if (!is.null(anno_tbl)) colnames(anno_tbl) else NULL,
                                status = decideTests(my_fit_eb, adjust.method = "none"),
                                p.adj.method = "none",
                                transform.counts = "none"
                            )
                        }
                    )

                    names(list_of_glimma_objs) <- colnames(my_fit_eb$coefficients)
                    de_output_object@interactive_volcano_plot <- list_of_glimma_objs

                    de_output_object
                }
            )

            list_of_objects
        }
)

