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

