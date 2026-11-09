# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-106 W5 (D5): proteomics' DA views, moved unchanged from the render handlers (mod_prot_da_render_volcano.R,
# mod_prot_da_render_heatmap.R, mod_prot_da_render_table.R, deleted) and the init handler's two outputs
# (mod_prot_da_handlers_init.R:318-336) at fe97aac1. Each view is a render's body, given the inputs and the tab's
# values it read (as promises, so each is read where the body read it); the shared frame (mod_omic_da.R) registers
# them with their render kinds kept. Kernels are called by name, so a test's namespace mock reaches them.
# Recorded, not fixed here: the cluster summary prints "NULL" after its no-cluster text (F106-03) and the download
# writes a placeholder (F106-04, P-26).

#' The contrasts text (mod_prot_da_handlers_init.R:319-326 at fe97aac1).
#' @param texts The group's `texts`.
#' @noRd
diaDaContrastsView <- function(contrastsAvailable, texts) {
    if (!is.null(contrastsAvailable)) {
        paste(contrastsAvailable, collapse = "\n")
    } else {
        texts[["contrasts_none"]]
    }
}

#' The status text (mod_prot_da_handlers_init.R:328-335 at fe97aac1).
#' @noRd
diaDaStatusView <- function(analysisComplete, contrastsAvailable, texts) {
    if (analysisComplete) {
        diaDaStatusText(length(contrastsAvailable))
    } else {
        texts[["status_waiting"]]
    }
}

#' The interactive volcano: the DA-time page, or why it is absent (mod_prot_da_render_volcano.R:17-28 at fe97aac1).
#' @noRd
diaDaVolcanoGlimmaView <- function(daResultsList, selectedContrast, daQValThresh, s4Object, outputDir) {
    shiny::req(selectedContrast, daResultsList)
    generateProtDAVolcanoPlotGlimma(
        da_results_list = daResultsList,
        selected_contrast = selectedContrast,
        da_q_val_thresh = daQValThresh,
        args_row_id = diaGuiProteinKey(s4Object),
        output_dir = outputDir
    ) %||% shiny::p(diaInteractiveVolcanoUnavailableText(selectedContrast)) # KIT-101 F-08: say why it is absent
}

#' The static volcano with its label inputs (mod_prot_da_render_volcano.R:32-44 at fe97aac1).
#' @noRd
diaDaStaticVolcanoView <- function(daResultsList, selectedContrast, daQValThresh, lfcThreshold, showLabels, labelTopN) {
    shiny::req(selectedContrast, daResultsList)
    generateProtDAVolcanoStatic(
        da_results_list = daResultsList,
        selected_contrast = selectedContrast,
        da_q_val_thresh = daQValThresh,
        lfc_threshold = lfcThreshold,
        show_labels = if (is.null(showLabels)) TRUE else showLabels,
        n_labels = if (is.null(labelTopN)) 10 else labelTopN
    )
}

#' The heatmap (mod_prot_da_render_heatmap.R:31-97 at fe97aac1): its clusters and drawing kept on the tab's values for
#' the cluster summary and Save Heatmap.
#' @param daData The tab's reactive values.
#' @noRd
diaDaHeatmapView <- function(daData, selectedContrast, topN, clusterMethod, distanceMethod, clustering, scaling, colorScheme,
    showLabels, daQValThresh, treeCutMethod, nClusters, cutHeight, minClusterSize, showDendro) {
    shiny::req(selectedContrast, daData$da_results_list)

    tryCatch(
        {
            # Create a compatible data structure for the plotting function
            plot_data_structure <- if (isDIAResultManifest(daData$da_results_list)) daData$da_results_list else list(
                da_proteins_long = daData$da_results_list$da_proteins_long,
                theObject = daData$da_results_list$theObject
            )

            # The heatmap contrast input should now match the comparison column in da_proteins_long
            cat(sprintf("   HEATMAP: Looking for contrast = %s\n", selectedContrast))

            heatmap_result <- generateProtDAHeatmap(
                da_results_list = plot_data_structure,
                selected_contrast = selectedContrast,
                top_n_genes = topN,
                clustering_method = clusterMethod,
                distance_method = distanceMethod,
                cluster_rows = clustering %in% c("both", "row"),
                cluster_cols = clustering %in% c("both", "column"),
                scale_data = scaling,
                color_scheme = colorScheme,
                show_gene_names = showLabels,
                da_q_val_thresh = daQValThresh,
                tree_cut_method = treeCutMethod,
                n_clusters = nClusters,
                cut_height = cutHeight,
                min_cluster_size = minClusterSize,
                show_dendrogram = if (is.null(showDendro)) TRUE else isTRUE(showDendro)
            )

            if (!is.null(heatmap_result)) {
                # Check if it's a list (new format) or just the plot (old format/fallback)
                if (is.list(heatmap_result) && "plot" %in% names(heatmap_result)) {
                    # Store cluster info and plot object
                    daData$current_row_clusters <- heatmap_result$row_clusters
                    daData$current_col_clusters <- heatmap_result$col_clusters
                    daData$current_heatmap_plot <- heatmap_result$plot

                    # Return plot for rendering
                    heatmap_result$plot
                } else {
                    # Fallback for old return format
                    heatmap_result
                }
            } else {
                # No significant genes found
                plot(1, 1,
                    type = "n", axes = FALSE, xlab = "", ylab = "",
                    main = paste("No significant genes found for contrast:", selectedContrast)
                )
                text(1, 1, "Adjust significance thresholds\nor select different contrast", cex = 1.2)
            }
        },
        error = function(e) {
            cat(paste("*** ERROR in heatmap generation:", e$message, "\n"))
            plot(1, 1,
                type = "n", axes = FALSE, xlab = "", ylab = "",
                main = "Error generating heatmap"
            )
            text(1, 1, paste("Error:", e$message), cex = 1.0)
        }
    )
}

#' The tree cut's cluster summary, printed (mod_prot_da_render_heatmap.R:104-130 at fe97aac1; KIT-111 F106-03: the
#' no-clusters text returns invisibly, so `renderPrint` prints no `NULL` after it).
#' @noRd
diaDaClusterSummaryView <- function(treeCutMethod, rowClusters) {
    shiny::req(treeCutMethod != "none")

    if (is.null(rowClusters)) {
        cat("No clusters defined. Enable clustering and tree cutting on the heatmap.")
        return(invisible(NULL))
    }

    clusters <- rowClusters
    n_clusters <- length(unique(clusters))

    cat(sprintf("Total Clusters: %d\n", n_clusters))
    cat("----------------------------------\n")

    unique_sorted_clusters <- sort(unique(clusters))
    purrr::walk(unique_sorted_clusters, function(i) {
        members <- names(clusters)[clusters == i]
        cat(sprintf("\nCluster %d (%d proteins):\n", i, length(members)))

        # If we have gene names, try to show them too
        # This requires looking up the proteins in the da_results list
        cat(paste(head(members, 20), collapse = ", "))
        if (length(members) > 20) {
            cat(sprintf(", ... and %d more", length(members) - 20))
        }
        cat("\n")
    })
}

#' The results table (mod_prot_da_render_table.R:6-99 at fe97aac1). KIT-024: on the native session it is paged
#' server-side from the store (the significance selection applied in the read); `pager` holds the reader and the
#' selection's row count for the table's filter.
#' @noRd
diaDaResultsTableView <- function(daResultsList, tableContrast, tableSignificance, daQValThresh, lfcCutoff, tableMaxRows, pager) {
    shiny::req(tableContrast, daResultsList)
    pager$read <- NULL

    tryCatch(
        {
            if (isDIAResultManifest(daResultsList)) {
                manifest <- daResultsList
                contrast <- tableContrast
                significance <- tableSignificance %||% "all"
                qvalue <- daQValThresh
                lfc <- lfcCutoff
                page_length <- as.integer(tableMaxRows %||% 25L)
                # The owner is resolved here, in the render; the page requests run outside any reactive domain.
                manager <- getDIAActiveManager()
                first <- readDIAResultPage(manifest, contrast, significance, qvalue, lfc, limit = max(1L, min(page_length, 10000L)),
                    manager = manager)
                pager$rows <- as.integer(if (identical(significance, "all")) first$full_population_rows else first$counts[[significance]])
                pager$read <- function(offset, limit) readDIAResultPage(manifest, contrast, significance, qvalue, lfc,
                    offset = offset, limit = limit, manager = manager)$data
                return(DT::datatable(first$data, options = diaDisplayTableOptions(page_length, scrollX = TRUE,
                    dom = "Bfrtip", buttons = c("copy", "csv", "excel")), extensions = "Buttons", rownames = FALSE) |>
                    DT::formatRound(columns = c("log2FC", "raw_pvalue", "fdr_qvalue"), digits = 4))
            }
            # Get DA results from the new format
            if (!is.null(daResultsList$da_proteins_long)) {
                # Debug: Show available contrasts in data vs selected contrast
                available_contrasts <- unique(daResultsList$da_proteins_long$comparison)
                cat(sprintf("   DA TABLE: Available contrasts in data: %s\n", paste(available_contrasts, collapse = ", ")))
                cat(sprintf("   DA TABLE: Selected contrast from UI: %s\n", tableContrast))

                # Filter for selected contrast (the table contrast now matches the comparison column)
                current_results <- daResultsList$da_proteins_long |>
                    dplyr::filter(comparison == tableContrast)

                cat(sprintf("   DA TABLE: Found %d rows for contrast %s\n", nrow(current_results), tableContrast))

                if (nrow(current_results) > 0) {
                    # Filter based on significance selection
                    if (tableSignificance == "significant") {
                        current_results <- current_results[current_results$fdr_qvalue < daQValThresh, ]
                    } else if (tableSignificance == "up") {
                        current_results <- current_results[current_results$fdr_qvalue < daQValThresh &
                            current_results$log2FC > lfcCutoff, ]
                    } else if (tableSignificance == "down") {
                        current_results <- current_results[current_results$fdr_qvalue < daQValThresh &
                            current_results$log2FC < -lfcCutoff, ]
                    }

                    # Limit rows
                    if (nrow(current_results) > tableMaxRows) {
                        current_results <- current_results[1:tableMaxRows, ]
                    }

                    # Select relevant columns for display - use the correct protein ID column
                    if (!is.null(daResultsList$theObject)) {
                        protein_id_column <- daResultsList$theObject@protein_id_column
                        cat(sprintf("   DA TABLE: Using protein ID column = %s\n", protein_id_column))
                    } else {
                        # Fallback: try to find protein ID column in the data
                        possible_protein_cols <- c("Protein.Ids", "uniprot_acc", "protein_id", "Protein_ID")
                        protein_id_column <- intersect(possible_protein_cols, names(current_results))[1]
                        cat(sprintf("   DA TABLE: Using fallback protein ID column = %s\n", protein_id_column))
                    }

                    display_columns <- c(protein_id_column, "log2FC", "raw_pvalue", "fdr_qvalue")
                    current_results <- current_results |>
                        dplyr::select(any_of(display_columns))

                    DT::datatable(
                        current_results,
                        options = list(
                            pageLength = 25,
                            scrollX = TRUE,
                            dom = "Bfrtip",
                            buttons = c("copy", "csv", "excel")
                        ),
                        extensions = "Buttons"
                    ) |>
                        DT::formatRound(columns = c("log2FC", "raw_pvalue", "fdr_qvalue"), digits = 4)
                } else {
                    # No results for this contrast
                    DT::datatable(data.frame(Message = "No results available for selected contrast"))
                }
            } else {
                # No DA results available
                DT::datatable(data.frame(Message = "No DE analysis results available"))
            }
        },
        error = function(e) {
            cat(paste("*** ERROR in DA results table:", e$message, "\n"))
            DT::datatable(data.frame(Message = paste("Error:", e$message)))
        }
    )
}

#' The summary statistics (mod_prot_da_render_table.R:103-149 at fe97aac1).
#' @noRd
diaDaSummaryStatsView <- function(daResultsList, tableContrast, daQValThresh, lfcCutoff) {
    shiny::req(tableContrast, daResultsList)

    tryCatch(
        {
            if (isDIAResultManifest(daResultsList)) {
                counts <- readDIAResultPage(daResultsList, tableContrast,
                    qvalue = daQValThresh, lfc = lfcCutoff, limit = 1L)$counts
                return(paste(sprintf("Total genes: %d", counts$total),
                    sprintf("Significant (q < %.3f): %d", daQValThresh, counts$significant),
                    sprintf("Up-regulated: %d", counts$up), sprintf("Down-regulated: %d", counts$down),
                    sprintf("Fold-change cutoff: %.2f", lfcCutoff), sep = "\n"))
            }
            if (!is.null(daResultsList$da_proteins_long)) {
                # Filter for selected contrast
                current_results <- daResultsList$da_proteins_long |>
                    dplyr::filter(comparison == tableContrast)

                if (nrow(current_results) > 0) {
                    total_genes <- nrow(current_results)
                    significant <- sum(current_results$fdr_qvalue < daQValThresh, na.rm = TRUE)
                    up_reg <- sum(current_results$fdr_qvalue < daQValThresh &
                        current_results$log2FC > lfcCutoff, na.rm = TRUE)
                    down_reg <- sum(current_results$fdr_qvalue < daQValThresh &
                        current_results$log2FC < -lfcCutoff, na.rm = TRUE)

                    paste(
                        sprintf("Total genes: %d", total_genes),
                        sprintf("Significant (q < %.3f): %d", daQValThresh, significant),
                        sprintf("Up-regulated: %d", up_reg),
                        sprintf("Down-regulated: %d", down_reg),
                        sprintf("Fold-change cutoff: %.2f", lfcCutoff),
                        sep = "\n"
                    )
                } else {
                    "No results available for selected contrast"
                }
            } else {
                "No DE analysis results available"
            }
        },
        error = function(e) {
            paste("Error calculating statistics:", e$message)
        }
    )
}

#' The download's content (mod_prot_da_render_table.R:156-159 at fe97aac1): a placeholder (F106-04, P-26).
#' @noRd
diaDaDownloadContent <- function(file) {
    # Placeholder for creating downloadable zip file
    writeLines("Differential expression results would be packaged here", file)
}
