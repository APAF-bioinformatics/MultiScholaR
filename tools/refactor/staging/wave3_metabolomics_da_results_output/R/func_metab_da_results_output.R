# ----------------------------------------------------------------------------
# outputMetabDaResultsAllContrasts
# ----------------------------------------------------------------------------
#' Output all metabolomics DE results to disk
#'
#' @description Writes DE results tables, volcano plots, and heatmaps to disk
#'   for all contrasts in a metabolomics DA analysis. Outputs are split by
#'   assay mode (posmode/negmode) and contrast, matching the proteomics workflow.
#'
#' @details Output filenames follow the pattern:
#'   - `de_posmoda_metabolites_{contrast}_long_annot.xlsx`
#'   - `de_negmoda_metabolites_{contrast}_long_annot.xlsx`
#'
#' @param da_results_list Results list from `runMetabolitesDA()`.
#' @param da_output_dir Directory for DE result tables.
#' @param publication_graphs_dir Directory for publication-quality figures.
#' @param da_q_val_thresh Q-value threshold for significance (default 0.05).
#' @param lfc_threshold Log fold-change threshold for volcano plot lines.
#' @param heatmap_top_n Number of top metabolites for heatmap (default 50).
#' @param heatmap_clustering Clustering option: "both", "row", "column", "none".
#' @param heatmap_color_scheme Color scheme for heatmap.
#'
#' @return TRUE if successful, FALSE otherwise.
#'
#' @importFrom vroom vroom_write
#' @importFrom writexl write_xlsx
#' @importFrom ggplot2 ggsave
#' @importFrom grDevices pdf dev.off png
#' @importFrom dplyr filter group_by summarise n mutate select
#' @importFrom purrr walk map
#' @importFrom logger log_info log_error log_warn
#' @export
outputMetabDaResultsAllContrasts <- function(
  da_results_list,
  da_output_dir,
  publication_graphs_dir,
  da_q_val_thresh = 0.05,
  lfc_threshold = 1,
  heatmap_top_n = 50,
  heatmap_clustering = "both",
  heatmap_color_scheme = "RdBu"
) {
    logger::log_info("--- Entering outputMetabDaResultsAllContrasts ---")
    logger::log_info(sprintf("   da_output_dir = %s", da_output_dir))
    logger::log_info(sprintf("   publication_graphs_dir = %s", publication_graphs_dir))

    # Validate inputs
    if (is.null(da_results_list) || is.null(da_results_list$da_metabolites_long)) {
        logger::log_error("   No DE results available")
        return(FALSE)
    }

    # Normalize paths
    if (!is.null(da_output_dir)) {
        da_output_dir <- gsub("//+", "/", da_output_dir)
        da_output_dir <- normalizePath(da_output_dir, winslash = "/", mustWork = FALSE)
    }

    if (!is.null(publication_graphs_dir)) {
        publication_graphs_dir <- gsub("//+", "/", publication_graphs_dir)
        publication_graphs_dir <- normalizePath(publication_graphs_dir, winslash = "/", mustWork = FALSE)
    }

    # Create output directories
    if (!dir.exists(da_output_dir)) {
        dir.create(da_output_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info(sprintf("   Created da_output_dir: %s", da_output_dir))
    }

    volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
    if (!dir.exists(volcano_dir)) {
        dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info(sprintf("   Created Volcano_Plots directory: %s", volcano_dir))
    }

    heatmap_dir <- file.path(publication_graphs_dir, "Heatmaps")
    if (!dir.exists(heatmap_dir)) {
        dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info(sprintf("   Created Heatmaps directory: %s", heatmap_dir))
    }

    numsigde_dir <- file.path(publication_graphs_dir, "NumSigDeMolecules")
    if (!dir.exists(numsigde_dir)) {
        dir.create(numsigde_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info(sprintf("   Created NumSigDeMolecules directory: %s", numsigde_dir))
    }

    # Get DE results with sample intensity columns (already included from runMetabolitesDA)
    da_metabolites_long <- da_results_list$da_metabolites_long

    # Get unique assays and contrasts
    assays <- unique(da_metabolites_long$assay)
    contrasts <- unique(da_metabolites_long$comparison)

    logger::log_info(sprintf("   Found %d assays: %s", length(assays), paste(assays, collapse = ", ")))
    logger::log_info(sprintf("   Found %d contrasts: %s", length(contrasts), paste(contrasts, collapse = ", ")))

    # Map assay names to mode prefixes
    # LCMS_Pos -> posmode, LCMS_Neg -> negmode
    get_mode_prefix <- function(assay_name) {
        if (grepl("pos", assay_name, ignore.case = TRUE)) {
            return("posmode")
        } else if (grepl("neg", assay_name, ignore.case = TRUE)) {
            return("negmode")
        } else {
            # Fallback: use sanitized assay name
            return(gsub("[^A-Za-z0-9]", "", tolower(assay_name)))
        }
    }

    # Storage for combined PDFs
    all_volcano_plots <- list()
    all_heatmap_plots <- list()
    all_numsig_tables <- list()

    # Process each ASSAY x CONTRAST combination
    for (assay_name in assays) {
        mode_prefix <- get_mode_prefix(assay_name)
        logger::log_info(sprintf("   Processing assay: %s (mode: %s)", assay_name, mode_prefix))

        for (contrast_name in contrasts) {
            logger::log_info(sprintf("      Processing contrast: %s", contrast_name))

            # Filter data for this assay + contrast
            assay_contrast_data <- da_metabolites_long |>
                dplyr::filter(assay == assay_name & comparison == contrast_name)

            if (nrow(assay_contrast_data) == 0) {
                logger::log_warn(sprintf("      No data for %s / %s, skipping", assay_name, contrast_name))
                next
            }

            # =====================================================================
            # 1. Write DE results tables (TSV and Excel)
            # Filename format: de_{mode}_metabolites_{contrast}_long_annot.xlsx
            # =====================================================================
            tryCatch(
                {
                    # Reorder columns: ID, name first, then stats, then intensity columns
                    priority_cols <- c(
                        "metabolite_id", "metabolite_name",
                        "logFC", "raw_pvalue", "fdr_qvalue", "significant",
                        "comparison", "friendly_name", "numerator", "denominator"
                    )
                    priority_cols <- intersect(priority_cols, colnames(assay_contrast_data))

                    # Get sample intensity columns
                    intensity_cols <- grep("^intensity\\.", colnames(assay_contrast_data), value = TRUE)

                    # Get any remaining columns (exclude assay since we're splitting by it)
                    other_cols <- setdiff(
                        colnames(assay_contrast_data),
                        c(priority_cols, intensity_cols, "assay")
                    )

                    # Final column order
                    final_col_order <- c(priority_cols, other_cols, sort(intensity_cols))
                    output_data <- assay_contrast_data[, final_col_order, drop = FALSE]

                    # Build filename: de_{mode}_metabolites_{contrast}_long_annot
                    file_base <- paste0("de_", mode_prefix, "_metabolites_", contrast_name, "_long_annot")

                    # Write TSV
                    tsv_path <- file.path(da_output_dir, paste0(file_base, ".tsv"))
                    vroom::vroom_write(output_data, tsv_path)
                    logger::log_info(sprintf(
                        "      Wrote TSV: %s (%d rows, %d cols)",
                        basename(tsv_path), nrow(output_data), ncol(output_data)
                    ))

                    # Write Excel
                    xlsx_path <- file.path(da_output_dir, paste0(file_base, ".xlsx"))
                    writexl::write_xlsx(output_data, xlsx_path)
                    logger::log_info(sprintf("      Wrote Excel: %s", basename(xlsx_path)))
                },
                error = function(e) {
                    logger::log_error(sprintf(
                        "      Error writing tables for %s / %s: %s",
                        assay_name, contrast_name, e$message
                    ))
                }
            )

            # =================================================================
            # 2. Generate and save volcano plots (per assay, per contrast)
            # =================================================================
            tryCatch(
                {
                    # Build a filtered results list for just this assay
                    assay_filtered_results <- da_results_list
                    assay_filtered_results$da_metabolites_long <- assay_contrast_data

                    volcano_plot <- generateMetabDAVolcanoStatic(
                        da_results_list = assay_filtered_results,
                        selected_contrast = contrast_name,
                        selected_assay = assay_name,
                        da_q_val_thresh = da_q_val_thresh,
                        lfc_threshold = lfc_threshold,
                        show_labels = TRUE,
                        n_labels = 15
                    )

                    if (!is.null(volcano_plot)) {
                        # Unique key for combined PDF
                        plot_key <- paste0(mode_prefix, "_", contrast_name)
                        all_volcano_plots[[plot_key]] <- volcano_plot

                        # Filename: {mode}_{contrast}_volcano
                        volcano_base <- paste0(mode_prefix, "_", contrast_name)

                        # Save PNG
                        volcano_png <- file.path(volcano_dir, paste0(volcano_base, ".png"))
                        ggplot2::ggsave(volcano_png, volcano_plot, width = 8, height = 7, dpi = 300)
                        logger::log_info(sprintf("      Saved volcano PNG: %s", basename(volcano_png)))

                        # Save PDF
                        volcano_pdf <- file.path(volcano_dir, paste0(volcano_base, ".pdf"))
                        ggplot2::ggsave(volcano_pdf, volcano_plot, width = 8, height = 7)
                        logger::log_info(sprintf("      Saved volcano PDF: %s", basename(volcano_pdf)))
                    }
                },
                error = function(e) {
                    logger::log_error(sprintf(
                        "      Error generating volcano for %s / %s: %s",
                        assay_name, contrast_name, e$message
                    ))
                }
            )

            # =================================================================
            # 3. Generate and save heatmaps (per assay, per contrast)
            # =================================================================
            tryCatch(
                {
                    # Build a filtered results list for just this assay
                    assay_filtered_results <- da_results_list
                    assay_filtered_results$da_metabolites_long <- assay_contrast_data

                    heatmap_obj <- generateMetabDAHeatmap(
                        da_results_list = assay_filtered_results,
                        selected_contrast = contrast_name,
                        selected_assay = assay_name,
                        top_n = heatmap_top_n,
                        clustering_method = "ward.D2",
                        distance_method = "euclidean",
                        cluster_rows = heatmap_clustering %in% c("both", "row"),
                        cluster_cols = heatmap_clustering %in% c("both", "column"),
                        scale_data = "row",
                        color_scheme = heatmap_color_scheme,
                        show_metabolite_names = FALSE,
                        da_q_val_thresh = da_q_val_thresh
                    )

                    if (!is.null(heatmap_obj)) {
                        # Unique key for combined PDF
                        plot_key <- paste0(mode_prefix, "_", contrast_name)
                        all_heatmap_plots[[plot_key]] <- heatmap_obj

                        # Filename: {mode}_{contrast}_heatmap
                        heatmap_base <- paste0(mode_prefix, "_", contrast_name)

                        # Save PNG
                        heatmap_png <- file.path(heatmap_dir, paste0(heatmap_base, "_heatmap.png"))
                        grDevices::png(heatmap_png, width = 10, height = 8, units = "in", res = 300)
                        ComplexHeatmap::draw(heatmap_obj)
                        grDevices::dev.off()
                        logger::log_info(sprintf("      Saved heatmap PNG: %s", basename(heatmap_png)))

                        # Save PDF
                        heatmap_pdf <- file.path(heatmap_dir, paste0(heatmap_base, "_heatmap.pdf"))
                        grDevices::pdf(heatmap_pdf, width = 10, height = 8)
                        ComplexHeatmap::draw(heatmap_obj)
                        grDevices::dev.off()
                        logger::log_info(sprintf("      Saved heatmap PDF: %s", basename(heatmap_pdf)))
                    } else {
                        logger::log_warn(sprintf(
                            "      No significant metabolites for heatmap: %s / %s",
                            assay_name, contrast_name
                        ))
                    }
                },
                error = function(e) {
                    logger::log_error(sprintf(
                        "      Error generating heatmap for %s / %s: %s",
                        assay_name, contrast_name, e$message
                    ))
                }
            )

            # =================================================================
            # 4. Calculate NumSigDeMolecules for this assay/contrast
            # =================================================================
            tryCatch(
                {
                    sig_summary <- assay_contrast_data |>
                        dplyr::summarise(
                            total = dplyr::n(),
                            significant = sum(significant != "NS", na.rm = TRUE),
                            up_regulated = sum(significant == "Up", na.rm = TRUE),
                            down_regulated = sum(significant == "Down", na.rm = TRUE),
                            .groups = "drop"
                        ) |>
                        dplyr::mutate(
                            assay = assay_name,
                            mode = mode_prefix,
                            contrast = contrast_name,
                            q_threshold = da_q_val_thresh
                        )

                    table_key <- paste0(mode_prefix, "_", contrast_name)
                    all_numsig_tables[[table_key]] <- sig_summary
                },
                error = function(e) {
                    logger::log_error(sprintf(
                        "      Error calculating NumSigDE for %s / %s: %s",
                        assay_name, contrast_name, e$message
                    ))
                }
            )
        } # End contrast loop
    } # End assay loop

    # =========================================================================
    # 5. Create combined volcano plots PDF
    # =========================================================================
    if (length(all_volcano_plots) > 0) {
        tryCatch(
            {
                combined_volcano_pdf <- file.path(volcano_dir, "all_volcano_plots_combined.pdf")
                grDevices::pdf(combined_volcano_pdf, width = 8, height = 7, onefile = TRUE)
                purrr::walk(all_volcano_plots, print)
                grDevices::dev.off()
                logger::log_info(sprintf(
                    "   Created combined volcano PDF: %d plots",
                    length(all_volcano_plots)
                ))
            },
            error = function(e) {
                logger::log_error(sprintf("   Error creating combined volcano PDF: %s", e$message))
            }
        )
    }

    # =========================================================================
    # 6. Create combined heatmaps PDF
    # =========================================================================
    if (length(all_heatmap_plots) > 0) {
        tryCatch(
            {
                combined_heatmap_pdf <- file.path(heatmap_dir, "all_heatmaps_combined.pdf")
                grDevices::pdf(combined_heatmap_pdf, width = 10, height = 8, onefile = TRUE)
                purrr::walk(all_heatmap_plots, function(hm) {
                    ComplexHeatmap::draw(hm)
                })
                grDevices::dev.off()
                logger::log_info(sprintf(
                    "   Created combined heatmap PDF: %d plots",
                    length(all_heatmap_plots)
                ))
            },
            error = function(e) {
                logger::log_error(sprintf("   Error creating combined heatmap PDF: %s", e$message))
            }
        )
    }

    # =========================================================================
    # 7. Write NumSigDeMolecules summary table
    # =========================================================================
    if (length(all_numsig_tables) > 0) {
        tryCatch(
            {
                combined_numsig <- dplyr::bind_rows(all_numsig_tables)

                # Write TSV
                numsig_tsv <- file.path(numsigde_dir, "metabolites_num_sig_de_molecules.tab")
                vroom::vroom_write(combined_numsig, numsig_tsv)
                logger::log_info(sprintf("   Wrote NumSigDE table: %s", basename(numsig_tsv)))

                # Write Excel
                numsig_xlsx <- file.path(numsigde_dir, "metabolites_num_sig_de_molecules.xlsx")
                writexl::write_xlsx(combined_numsig, numsig_xlsx)
                logger::log_info(sprintf("   Wrote NumSigDE Excel: %s", basename(numsig_xlsx)))

                # Create bar plot
                numsig_plot <- ggplot2::ggplot(
                    combined_numsig,
                    ggplot2::aes(x = contrast, y = significant, fill = mode)
                ) +
                    ggplot2::geom_bar(stat = "identity", position = "dodge") +
                    ggplot2::labs(
                        title = "Number of Significant DE Metabolites",
                        x = "Contrast",
                        y = "Number of Significant Metabolites",
                        fill = "Assay"
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(
                        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
                    )

                # Save bar plot
                numsig_png <- file.path(numsigde_dir, "metabolites_num_sig_barplot.png")
                ggplot2::ggsave(numsig_png, numsig_plot, width = 10, height = 6, dpi = 300)
                logger::log_info(sprintf("   Saved NumSigDE barplot: %s", basename(numsig_png)))
            },
            error = function(e) {
                logger::log_error(sprintf("   Error writing NumSigDE summary: %s", e$message))
            }
        )
    }

    logger::log_info("--- Exiting outputMetabDaResultsAllContrasts ---")
    return(TRUE)
}

