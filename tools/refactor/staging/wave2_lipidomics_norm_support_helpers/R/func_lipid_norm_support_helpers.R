#' Build ITSD Selection Table for DT Display
#'
#' @title Build ITSD selection table for DT display
#' @description Creates a data frame suitable for DT::datatable display with
#'              pre-selection of candidate internal standards based on pattern matching.
#'
#' @param assay_data Data frame for a single assay
#' @param lipid_id_col Column name for lipid IDs
#' @param annotation_cols Character vector of columns to search for IS patterns
#' @param is_patterns Named list of regex patterns for IS detection
#'
#' @return Data frame with columns: feature_id, annotation, mean_intensity, cv_percent, is_candidate
#'
#' @importFrom dplyr select mutate across all_of any_of filter arrange desc
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect str_to_lower
#' @importFrom stats sd
#' @noRd
buildLipidItsdSelectionTable <- function(
    assay_data
    , lipid_id_col
    , annotation_cols = NULL
    , is_patterns = list(
        itsd = "(?i)itsd|istd|internal.?standard"
        , rt = "(?i)^rt_|_rt$|retention.?time"
        , isd = "(?i)^is_|_is$|^isd_|_isd$"
        , labeled = "(?i)_d[0-9]+$|_13c|_15n|deuterated|labeled"
    )
) {
    # --- Input validation ---
    stopifnot(
        is.data.frame(assay_data)
        , is.character(lipid_id_col)
        , lipid_id_col %in% colnames(assay_data)
    )

    # --- Identify sample columns (numeric) vs metadata columns ---
    numeric_cols <- names(assay_data)[sapply(assay_data, is.numeric)]
    sample_cols <- setdiff(numeric_cols, lipid_id_col)

    if (length(sample_cols) == 0) {
        warning("No sample columns found in assay_data")
        return(data.frame(
            feature_id = character(0)
            , annotation = character(0)
            , mean_intensity = numeric(0)
            , cv_percent = numeric(0)
            , is_candidate = logical(0)
        ))
    }

    # --- Determine annotation columns to search ---
    if (is.null(annotation_cols)) {
        # Default: search lipid ID column and common annotation columns
        potential_anno_cols <- c(
            lipid_id_col
            , "lipid"
            , "lipid_identification"
            , "annotation"
            , "compound_name"
            , "name"
        )
        annotation_cols <- intersect(potential_anno_cols, colnames(assay_data))
    }

    if (length(annotation_cols) == 0) {
        annotation_cols <- lipid_id_col
    }

    # --- Calculate summary statistics per feature ---
    feature_stats <- assay_data |>
        dplyr::mutate(
            feature_id = as.character(.data[[lipid_id_col]])
            , mean_intensity = rowMeans(
                dplyr::across(dplyr::all_of(sample_cols))
                , na.rm = TRUE
            )
            , sd_intensity = apply(
                dplyr::across(dplyr::all_of(sample_cols))
                , 1
                , sd
                , na.rm = TRUE
            )
        ) |>
        dplyr::mutate(
            cv_percent = ifelse(
                mean_intensity > 0
                , (sd_intensity / mean_intensity) * 100
                , NA_real_
            )
        )

    # --- Create annotation string for pattern matching ---
    if (length(annotation_cols) > 1) {
        feature_stats$annotation <- apply(
            feature_stats[, annotation_cols, drop = FALSE]
            , 1
            , \(x) paste(na.omit(x), collapse = " | ")
        )
    } else {
        feature_stats$annotation <- as.character(feature_stats[[annotation_cols[1]]])
    }

    # --- Detect IS candidates using patterns ---
    combined_pattern <- paste(unlist(is_patterns), collapse = "|")

    feature_stats$is_candidate <- stringr::str_detect(
        stringr::str_to_lower(feature_stats$annotation)
        , combined_pattern
    )

    # --- Select and arrange output columns ---
    result <- feature_stats |>
        dplyr::select(
            feature_id
            , annotation
            , mean_intensity
            , cv_percent
            , is_candidate
        ) |>
        dplyr::arrange(dplyr::desc(is_candidate), cv_percent)

    return(result)
}

#' Generate QC Plots for Lipidomics Data and Save to Disk
#'
#' @title Generate QC plots for lipidomics data and save to disk
#' @description Generates PCA, RLE, Density, and Correlation plots for each assay
#'              in a LipidomicsAssayData object and saves them as PNG files.
#'
#' @param theObject LipidomicsAssayData S4 object
#' @param experiment_paths List of experiment directory paths (must contain lipid_qc_dir)
#' @param stage One of "post_filter", "post_norm", "ruv_corrected"
#' @param grouping_variable Column name from design matrix for plot coloring
#' @param shape_variable Column name from design matrix for plot shapes (optional)
#' @param font_size Font size for plots
#'
#' @return Invisible list of file paths saved (named by assay and plot type)
#'
#' @importFrom ggplot2 ggsave
#' @importFrom purrr imap set_names
#' @importFrom logger log_info log_warn log_error
#' @noRd
generateLipidQcPlots <- function(
    theObject
    , experiment_paths
    , stage = c("post_filter", "post_norm", "ruv_corrected")
    , grouping_variable = NULL
    , shape_variable = NULL
    , font_size = 8
) {
    # --- Input validation ---
    stage <- match.arg(stage)

    stopifnot(
        inherits(theObject, "LipidomicsAssayData")
        , !is.null(experiment_paths$lipid_qc_dir)
    )

    qc_dir <- experiment_paths$lipid_qc_dir
    if (!dir.exists(qc_dir)) {
        dir.create(qc_dir, recursive = TRUE)
    }

    # --- Get assay names ---
    assay_list <- theObject@lipid_data
    assay_names <- names(assay_list)

    if (length(assay_names) == 0) {
        logger::log_warn("No assays found in LipidomicsAssayData object")
        return(invisible(list()))
    }

    # --- Determine grouping variable ---
    if (is.null(grouping_variable)) {
        grouping_variable <- theObject@group_id
    }

    # --- Stage prefix mapping ---
    stage_prefix <- switch(
        stage
        , "post_filter" = "pre_norm"
        , "post_norm" = "post_norm"
        , "ruv_corrected" = "ruv_corrected"
    )

    # --- Generate plots for each assay ---
    saved_paths <- list()

    for (assay_name in assay_names) {
        safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
        logger::log_info(paste("Generating QC plots for assay:", assay_name))

        # --- PCA Plot ---
        tryCatch({
            pca_plots <- plotPca(
                theObject
                , grouping_variable = grouping_variable
                , label_column = ""
                , shape_variable = shape_variable
                , title = ""
                , font_size = font_size
            )

            # Select plot for this assay
            pca_plot <- if (is.list(pca_plots) && assay_name %in% names(pca_plots)) {
                pca_plots[[assay_name]]
            } else if (inherits(pca_plots, "ggplot")) {
                pca_plots
            } else if (is.list(pca_plots) && length(pca_plots) > 0) {
                pca_plots[[1]]
            } else {
                NULL
            }

            if (!is.null(pca_plot)) {
                pca_path <- file.path(qc_dir, sprintf("%s_%s_pca.png", safe_name, stage_prefix))
                ggplot2::ggsave(pca_path, pca_plot, width = 8, height = 6, dpi = 150)
                saved_paths[[paste0(assay_name, "_pca")]] <- pca_path
            }
        }, error = function(e) {
            logger::log_warn(paste("Error generating PCA for", assay_name, ":", e$message))
        })

        # --- RLE Plot ---
        tryCatch({
            rle_plots <- plotRle(
                theObject
                , group = grouping_variable
                , yaxis_limit = c(-6, 6)
            )

            rle_plot <- if (is.list(rle_plots) && assay_name %in% names(rle_plots)) {
                rle_plots[[assay_name]]
            } else if (inherits(rle_plots, "ggplot")) {
                rle_plots
            } else if (is.list(rle_plots) && length(rle_plots) > 0) {
                rle_plots[[1]]
            } else {
                NULL
            }

            if (!is.null(rle_plot)) {
                rle_path <- file.path(qc_dir, sprintf("%s_%s_rle.png", safe_name, stage_prefix))
                ggplot2::ggsave(rle_path, rle_plot, width = 10, height = 6, dpi = 150)
                saved_paths[[paste0(assay_name, "_rle")]] <- rle_path
            }
        }, error = function(e) {
            logger::log_warn(paste("Error generating RLE for", assay_name, ":", e$message))
        })

        # --- Density Plot ---
        tryCatch({
            # Generate PCA first for density plot
            pca_for_density <- plotPca(
                theObject
                , grouping_variable = grouping_variable
                , label_column = ""
                , shape_variable = shape_variable
                , title = ""
                , font_size = font_size
            )

            density_plots <- plotDensity(
                theObject
                , grouping_variable = grouping_variable
                , title = ""
                , font_size = font_size
            )

            density_plot <- if (is.list(density_plots) && assay_name %in% names(density_plots)) {
                density_plots[[assay_name]]
            } else if (inherits(density_plots, "ggplot")) {
                density_plots
            } else if (is.list(density_plots) && length(density_plots) > 0) {
                density_plots[[1]]
            } else {
                NULL
            }

            if (!is.null(density_plot)) {
                density_path <- file.path(qc_dir, sprintf("%s_%s_density.png", safe_name, stage_prefix))
                ggplot2::ggsave(density_path, density_plot, width = 8, height = 6, dpi = 150)
                saved_paths[[paste0(assay_name, "_density")]] <- density_path
            }
        }, error = function(e) {
            logger::log_warn(paste("Error generating Density for", assay_name, ":", e$message))
        })

        # --- Pearson Correlation Plot ---
        tryCatch({
            corr_plots <- plotPearson(
                theObject
                , correlation_group = grouping_variable
            )

            corr_plot <- if (is.list(corr_plots) && assay_name %in% names(corr_plots)) {
                corr_plots[[assay_name]]
            } else if (inherits(corr_plots, "ggplot")) {
                corr_plots
            } else if (is.list(corr_plots) && length(corr_plots) > 0) {
                corr_plots[[1]]
            } else {
                NULL
            }

            if (!is.null(corr_plot)) {
                corr_path <- file.path(qc_dir, sprintf("%s_%s_correlation.png", safe_name, stage_prefix))
                ggplot2::ggsave(corr_path, corr_plot, width = 10, height = 8, dpi = 150)
                saved_paths[[paste0(assay_name, "_correlation")]] <- corr_path
            }
        }, error = function(e) {
            logger::log_warn(paste("Error generating Correlation for", assay_name, ":", e$message))
        })
    }

    logger::log_info(paste("Generated", length(saved_paths), "QC plots for stage:", stage))
    return(invisible(saved_paths))
}

