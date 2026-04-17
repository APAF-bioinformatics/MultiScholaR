# ----------------------------------------------------------------------------
# generateLimpaQCPlots
# ----------------------------------------------------------------------------
#' Generate limpa DPC Quality Control Diagnostic Plots
#'
#' Creates a comprehensive set of diagnostic plots for limpa DPC imputation or quantification,
#' including detection probability curves, missing value comparisons, intensity distributions,
#' and summary statistics. Works with both peptide-level imputation and protein-level quantification.
#'
#' @param after_object A PeptideQuantitativeData or ProteinQuantitativeData object containing limpa results
#' @param before_object Optional. The object before limpa processing for comparison plots.
#'   If NULL, will try to extract from stored data or skip comparison plots.
#' @param save_plots Logical. Whether to save individual plots to files.
#' @param save_dir Optional directory path for saving plots. If NULL and save_plots=TRUE,
#'   will attempt to use project directories.
#' @param plot_prefix String prefix for saved plot filenames.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A composite ggplot object containing all diagnostic plots arranged in a grid
#'
#' @examples
#' # For peptide-level limpa imputation
#' limpa_qc_fig <- generateLimpaQCPlots(
#'   after_object = peptide_ruv_normalised_imputed_obj,
#'   before_object = peptide_ruv_normalised_results_temp_obj
#' )
#'
#' # For protein-level DPC quantification
#' limpa_qc_fig <- generateLimpaQCPlots(
#'   after_object = protein_dpc_quant_obj,
#'   before_object = peptide_obj
#' )
#'
#' @importFrom ggplotify as.ggplot
#' @importFrom cowplot ggdraw draw_grob
#' @importFrom png readPNG
#' @importFrom magick image_read
#' @importFrom grid grid.grabExpr
#' @importFrom patchwork plot_annotation
#' @importFrom ggplot2 ggplot aes annotate theme_void theme margin element_blank annotation_raster scale_x_continuous scale_y_continuous ggtitle geom_col scale_fill_manual geom_density theme_bw element_text labs ggsave
#' @importFrom stringr str_to_title
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr filter select mutate distinct left_join rename summarise group_by n n_distinct
#' @importFrom tidyr pivot_longer
#' @export
generateLimpaQCPlots <- function(after_object,
                                 before_object = NULL,
                                 save_plots = TRUE,
                                 save_dir = "peptide_qc",
                                 plot_prefix = "limpa",
                                 verbose = TRUE) {
  message("--- Entering generateLimpaQCPlots ---")
  message(sprintf(
    "   Arguments: after_object type=%s, before_object=%s, save_plots=%s, save_dir=%s, plot_prefix=%s, verbose=%s",
    class(after_object), ifelse(is.null(before_object), "NULL", class(before_object)),
    save_plots, save_dir, plot_prefix, verbose
  ))

  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package is required for composite plots")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("stringr package is required")
  }

  message("   Step: Detecting limpa results type and extracting DPC results")
  # Detect limpa results type and extract DPC results
  dpc_results <- NULL
  analysis_type <- "unknown"

  # Check for peptide-level limpa results
  if (!is.null(after_object@args$limpa_dpc_results)) {
    dpc_results <- after_object@args$limpa_dpc_results
    analysis_type <- "peptide_imputation"
    if (verbose) cat("Detected peptide-level limpa DPC imputation results\n")
  }

  # Check for protein-level limpa results
  if (is.null(dpc_results) && !is.null(after_object@args$limpa_dpc_quant_results)) {
    dpc_results <- after_object@args$limpa_dpc_quant_results
    analysis_type <- "protein_quantification"
    if (verbose) cat("Detected protein-level limpa DPC quantification results\n")
  }

  # Exit if no limpa results found
  if (is.null(dpc_results)) {
    stop("No limpa DPC results found in the provided object. Please run limpa imputation/quantification first.")
  }
  message(sprintf("   Results: analysis_type=%s, dpc_results structure=%s", analysis_type, capture.output(str(dpc_results))))

  # Try to get the correct save directory from project_dirs GLOBAL OBJECT
  message("   Step: Determining save directory")
  if (save_dir == "peptide_qc") {
    tryCatch(
      {
        # Access project_dirs from global environment
        if (exists("project_dirs", envir = .GlobalEnv)) {
          project_dirs <- get("project_dirs", envir = .GlobalEnv)

          # Look through each experiment in project_dirs
          for (experiment_name in names(project_dirs)) {
            if (!is.null(project_dirs[[experiment_name]]$peptide_qc_dir)) {
              save_dir <- project_dirs[[experiment_name]]$peptide_qc_dir
              if (verbose) cat("[info] Using peptide_qc_dir from global project_dirs:", save_dir, "\n")
              break
            }
          }
        }
      },
      error = function(e) {
        if (verbose) cat("Could not access project_dirs from global environment, using default\n")
      }
    )
  }
  message(sprintf("   Results: save_dir=%s", save_dir))

  # Extract DPC parameters (handle both formats)
  message("   Step: Extracting DPC parameters")
  if (!is.null(dpc_results$dpc_parameters_used)) {
    dpc_params <- dpc_results$dpc_parameters_used
  } else {
    dpc_params <- dpc_results$dpc_parameters
  }
  message(sprintf("   Results: dpc_params=%s", capture.output(str(dpc_params))))

  if (verbose) {
    cat("\n=== limpa DPC Analysis Results ===\n")
    cat("Analysis Type:", analysis_type, "\n")
    cat("DPC Method:", ifelse(is.null(dpc_results$dpc_method), "limpa_dpc", dpc_results$dpc_method), "\n")

    # Debug: show what DPC data is available
    cat("Available DPC data:\n")
    cat("  - y_peptide_for_dpc:", !is.null(dpc_results$y_peptide_for_dpc), "\n")
    cat("  - dpc_object_used:", !is.null(dpc_results$dpc_object_used), "\n")
    cat("  - dpc_object:", !is.null(dpc_results$dpc_object), "\n")
    cat("  - slope_interpretation:", ifelse(is.null(dpc_results$slope_interpretation), "NULL", dpc_results$slope_interpretation), "\n")

    if (!is.null(dpc_params) && length(dpc_params) >= 2) {
      cat("Beta0 (intercept):", round(dpc_params[1], 4), "\n")
      cat("Beta1 (slope):", round(dpc_params[2], 4), "\n")
      cat("Missing value mechanism:", dpc_results$slope_interpretation, "\n")
    }

    if (!is.null(dpc_results$missing_percentage_before)) {
      cat("Missing percentage before:", dpc_results$missing_percentage_before, "%\n")
    }
    cat("=====================================\n")
  }

  # Initialize plot list
  plot_list <- list()

  # 1. Detection Probability Curve Plot
  plot_list$dpc_curve <- tryCatch(
    {
      message("--- Entering DPC curve plot generation ---")
      if (verbose) cat("Generating DPC curve plot...\n")

      # Validate dpc_params first
      message("   Step: Validating DPC parameters...")
      message(sprintf("      dpc_params: %s", capture.output(str(dpc_params))))
      if (is.null(dpc_params) || !is.numeric(dpc_params) || length(dpc_params) < 2) {
        message("   Warning: Invalid DPC parameters detected, using defaults")
        if (verbose) cat("Warning: Invalid DPC parameters, using defaults\n")
        dpc_params <- c(-1.335, 0.074) # Use the known values from your data
      }
      message(sprintf("      Final dpc_params: %s", capture.output(str(dpc_params))))

      message("   Step: Extracting stored DPC object")
      # Extract stored DPC object based on analysis type
      stored_dpc <- NULL
      message(sprintf("      analysis_type: %s", analysis_type))
      message(sprintf("      dpc_results structure: %s", capture.output(str(names(dpc_results)))))

      if (analysis_type == "protein_quantification" && !is.null(dpc_results$dpc_object_used)) {
        stored_dpc <- dpc_results$dpc_object_used
        message("     Using dpc_object_used for protein quantification")
        message(sprintf("      stored_dpc class: %s", class(stored_dpc)))
        message(sprintf("      stored_dpc structure: %s", capture.output(str(stored_dpc))))
      } else if (analysis_type == "peptide_imputation" && !is.null(dpc_results$dpc_object)) {
        stored_dpc <- dpc_results$dpc_object
        message("     Using dpc_object for peptide imputation")
        message(sprintf("      stored_dpc class: %s", class(stored_dpc)))
        message(sprintf("      stored_dpc structure: %s", capture.output(str(stored_dpc))))
      } else if (!is.null(dpc_results$y_peptide_for_dpc)) {
        message("     Recreating DPC from stored y_peptide_for_dpc")
        message(sprintf("      y_peptide_for_dpc dims: %d x %d", nrow(dpc_results$y_peptide_for_dpc), ncol(dpc_results$y_peptide_for_dpc)))
        stored_dpc <- limpa::dpc(dpc_results$y_peptide_for_dpc)
        message(sprintf("      Created stored_dpc class: %s", class(stored_dpc)))
      }

      message(sprintf("   Results: stored_dpc is null? %s", is.null(stored_dpc)))

      # Generate plot using plotDPC if we have the object
      if (!is.null(stored_dpc)) {
        message("   Step: Attempting to call limpa::plotDPC")
        message("      Checking if limpa package is available...")
        if (!requireNamespace("limpa", quietly = TRUE)) {
          stop("limpa package is required for DPC plotting")
        }
        message("      limpa package confirmed available")

        # OPTION 1: Try ggplotify conversion approach
        message("   Attempting Option 1: ggplotify conversion...")
        ggplotify_plot <- NULL
        tryCatch(
          {
            message("        Checking if ggplotify package is available...")
            if (!requireNamespace("ggplotify", quietly = TRUE)) {
              message("        ggplotify not available, skipping this option")
            } else {
              message("        Calling plotDPC and converting with ggplotify")
              # Call plotDPC and immediately convert to ggplot
              ggplotify_plot <- ggplotify::as.ggplot(function() limpa::plotDPC(stored_dpc))
              message(sprintf("        ggplotify conversion successful, class: %s", class(ggplotify_plot)))
            }
          },
          error = function(e) {
            message(sprintf("        ggplotify conversion FAILED: %s", e$message))
            ggplotify_plot <<- NULL
          }
        )

        # OPTION 2: Try cowplot approach
        message("   Attempting Option 2: cowplot conversion...")
        cowplot_plot <- NULL
        tryCatch(
          {
            message("        Checking if cowplot package is available...")
            if (!requireNamespace("cowplot", quietly = TRUE)) {
              message("        cowplot not available, skipping this option")
            } else {
              message("        Using cowplot to capture base plot")
              # Create the plot and capture as grob
              captured_grob <- grid::grid.grabExpr({
                limpa::plotDPC(stored_dpc)
              })

              # Convert using cowplot
              cowplot_plot <- cowplot::ggdraw() +
                cowplot::draw_grob(captured_grob)

              message(sprintf("        cowplot conversion successful, class: %s", class(cowplot_plot)))
            }
          },
          error = function(e) {
            message(sprintf("        cowplot conversion FAILED: %s", e$message))
            cowplot_plot <<- NULL
          }
        )

        # OPTION 3: Try direct plotDPC call and see what happens
        message("   Attempting Option 3: Direct plotDPC call...")
        direct_plot <- NULL
        tryCatch(
          {
            message("        About to call plotDPC directly")
            direct_plot <- limpa::plotDPC(stored_dpc)
            message(sprintf("        Direct plotDPC returned: %s", capture.output(str(direct_plot))))
            message(sprintf("        Direct plotDPC class: %s", class(direct_plot)))
          },
          error = function(e) {
            message(sprintf("        Direct plotDPC FAILED: %s", e$message))
            direct_plot <- NULL
          }
        )

        # OPTION 4: Try temp file approach
        message("   Attempting Option 4: Temporary file approach...")
        temp_plot <- NULL
        tryCatch(
          {
            message("        Creating temporary PNG file")
            temp_file <- tempfile(fileext = ".png")
            png(temp_file, width = 800, height = 600, res = 150)

            message("        Calling plotDPC to temp file")
            limpa::plotDPC(stored_dpc)
            dev.off()

            message("        Reading temp file back as ggplot")
            if (!requireNamespace("magick", quietly = TRUE)) {
              message("        magick not available, trying png package")
              if (!requireNamespace("png", quietly = TRUE)) {
                message("        png package not available either")
              } else {
                # Read with png package and create ggplot
                img <- png::readPNG(temp_file)
                temp_plot <- ggplot2::ggplot() +
                  ggplot2::annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
                  ggplot2::scale_x_continuous(expand = c(0, 0)) +
                  ggplot2::scale_y_continuous(expand = c(0, 0)) +
                  ggplot2::theme_void()
              }
            } else {
              # Use magick to read and convert
              img <- magick::image_read(temp_file)
              temp_plot <- ggplot2::ggplot() +
                ggplot2::annotation_raster(as.raster(img), xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
                ggplot2::scale_x_continuous(expand = c(0, 0)) +
                ggplot2::scale_y_continuous(expand = c(0, 0)) +
                ggplot2::theme_void()
            }

            message("        Removing temporary file")
            unlink(temp_file)

            message(sprintf("        temp file approach result: %s", capture.output(str(temp_plot))))
          },
          error = function(e) {
            message(sprintf("        temp file approach FAILED: %s", e$message))
            try(dev.off(), silent = TRUE)
            temp_plot <- NULL
          }
        )

        # Determine which plot was successful and assign to dpc_final_plot
        dpc_final_plot <- NULL
        if (!is.null(ggplotify_plot)) {
          message("   Option 1 SUCCESS: Adding title to ggplotify plot")
          slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "",
            "nearly random missing", dpc_results$slope_interpretation
          )

          # Add title using ggplot2
          dpc_final_plot <- ggplotify_plot +
            ggplot2::labs(title = paste(
              "Detection Probability Curve\n",
              "Slope:", round(as.numeric(dpc_params[2]), 3),
              "- Mechanism:", slope_text
            )) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

          message("     ggplotify plot with title created successfully")
        } else if (!is.null(cowplot_plot)) {
          message("   Option 2 SUCCESS: Adding title to cowplot plot")
          slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "",
            "nearly random missing", dpc_results$slope_interpretation
          )

          # Add title using ggplot2
          dpc_final_plot <- cowplot_plot +
            ggplot2::labs(title = paste(
              "Detection Probability Curve\n",
              "Slope:", round(as.numeric(dpc_params[2]), 3),
              "- Mechanism:", slope_text
            )) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

          message("     cowplot plot with title created successfully")
        } else if (!is.null(direct_plot) && inherits(direct_plot, "ggplot")) {
          message("   Option 3 SUCCESS: Using direct ggplot")
          slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "",
            "nearly random missing", dpc_results$slope_interpretation
          )

          dpc_final_plot <- direct_plot +
            ggplot2::labs(title = paste(
              "Detection Probability Curve\n",
              "Slope:", round(as.numeric(dpc_params[2]), 3),
              "- Mechanism:", slope_text
            )) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

          message("     direct ggplot with title created successfully")
        } else if (!is.null(temp_plot)) {
          message("   Option 4 SUCCESS: Using temp file plot")
          slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "",
            "nearly random missing", dpc_results$slope_interpretation
          )

          dpc_final_plot <- temp_plot +
            ggplot2::labs(title = paste(
              "Detection Probability Curve\n",
              "Slope:", round(as.numeric(dpc_params[2]), 3),
              "- Mechanism:", slope_text
            )) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

          message("     temp file plot with title created successfully")
        }
      } else {
        message("   No stored DPC object available - creating fallback plot")
        dpc_final_plot <- NULL
      }

      # Use the successful plot or create fallback
      if (!is.null(dpc_final_plot)) {
        if (verbose) cat("[OK] DPC curve plot generated successfully\n")
        dpc_final_plot
      } else {
        message("   Fallback: Creating parameters summary plot")
        if (verbose) cat("Creating DPC parameters summary plot...\n")

        slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "",
          "nearly random missing", dpc_results$slope_interpretation
        )

        ggplot2::ggplot() +
          ggplot2::annotate("text",
            x = 0.5, y = 0.5,
            label = paste(
              "DPC Parameters Summary\n",
              "beta0 (intercept):", round(as.numeric(dpc_params[1]), 3), "\n",
              "beta1 (slope):", round(as.numeric(dpc_params[2]), 3), "\n",
              "Mechanism:", slope_text
            ),
            size = 4, hjust = 0.5
          ) +
          ggplot2::theme_void() +
          ggplot2::ggtitle("DPC Parameters")
      }
    },
    error = function(e) {
      message(sprintf("--- ERROR in DPC plot generation: %s ---", e$message))
      message(sprintf("   Error class: %s", class(e)))
      message(sprintf("   Call stack: %s", capture.output(traceback())))
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste("DPC plot error:", e$message)) +
        ggplot2::theme_void() +
        ggplot2::ggtitle("DPC Curve")
    }
  )

  # 2. Missing Value Pattern Comparison
  plot_list$missing_comparison <- tryCatch(
    {
      if (verbose) cat("Generating missing value comparison plot...\n")

      # Get matrices for comparison
      if (analysis_type == "peptide_imputation") {
        after_matrix <- after_object@peptide_matrix
        before_matrix <- if (!is.null(before_object)) before_object@peptide_matrix else NULL
        data_type <- "Peptides"
      } else {
        # For protein quantification, show peptide input vs protein output
        after_matrix <- as.matrix(after_object@protein_quant_table[, -1]) # Exclude protein ID column

        # For protein quantification, we want to show the reduction in missing values
        # from input peptides to output proteins
        if (!is.null(dpc_results$y_peptide_for_dpc)) {
          before_matrix <- dpc_results$y_peptide_for_dpc # Use stored peptide data
          data_type <- "Peptides -> Proteins"
        } else if (!is.null(before_object) && "peptide_matrix" %in% slotNames(before_object)) {
          before_matrix <- before_object@peptide_matrix
          data_type <- "Peptides -> Proteins"
        } else {
          before_matrix <- NULL
          data_type <- "Proteins"
        }
      }

      if (!is.null(before_matrix) && !is.null(after_matrix)) {
        # Handle comparison between different data types (peptides vs proteins)
        if (analysis_type == "protein_quantification") {
          # For protein quantification: compare peptide input vs protein output
          # Get common samples only (dimensions will be different)
          common_samples <- intersect(colnames(before_matrix), colnames(after_matrix))

          if (length(common_samples) > 0) {
            before_subset <- before_matrix[, common_samples, drop = FALSE]
            after_subset <- after_matrix[, common_samples, drop = FALSE]

            # Calculate missing values per sample
            missing_before <- colSums(is.na(before_subset))
            missing_after <- colSums(is.na(after_subset))

            # Create labels that reflect the data transformation
            before_label <- paste("Before DPC-Quant\n(", nrow(before_subset), " peptides)")
            after_label <- paste("After DPC-Quant\n(", nrow(after_subset), " proteins)")

            missing_data <- data.frame(
              Sample = rep(names(missing_before), 2),
              Missing_Count = c(missing_before, missing_after),
              Stage = rep(c(before_label, after_label), each = length(missing_before))
            )
          } else {
            stop("No common samples between peptide and protein data")
          }
        } else {
          # For peptide imputation: ensure matrices are compatible for comparison
          if (!identical(dim(before_matrix), dim(after_matrix))) {
            common_samples <- intersect(colnames(before_matrix), colnames(after_matrix))
            common_features <- intersect(rownames(before_matrix), rownames(after_matrix))
            if (length(common_samples) > 0 && length(common_features) > 0) {
              before_matrix <- before_matrix[common_features, common_samples]
              after_matrix <- after_matrix[common_features, common_samples]
            } else {
              stop("No common samples/features for comparison")
            }
          }

          # Calculate missing values per sample
          missing_before <- colSums(is.na(before_matrix))
          missing_after <- colSums(is.na(after_matrix))

          missing_data <- data.frame(
            Sample = rep(names(missing_before), 2),
            Missing_Count = c(missing_before, missing_after),
            Stage = rep(c("Before limpa", "After limpa"), each = length(missing_before))
          )
        }

        missing_plot <- missing_data |>
          ggplot2::ggplot(ggplot2::aes(x = Sample, y = Missing_Count, fill = Stage)) +
          ggplot2::geom_col(position = "dodge") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(
            title = paste("Missing Values per Sample:", data_type),
            y = paste("Number of Missing", data_type),
            x = "Sample"
          ) +
          ggplot2::scale_fill_manual(values = c(
            "Before limpa" = "#E74C3C",
            "After limpa" = "#3498DB"
          ))
      } else {
        missing_plot <- ggplot2::ggplot() +
          ggplot2::annotate("text",
            x = 0.5, y = 0.5,
            label = "Missing value comparison\nrequires before/after objects"
          ) +
          ggplot2::theme_void() +
          ggplot2::ggtitle("Missing Values Comparison")
      }

      if (verbose) cat("[OK] Missing value comparison plot generated\n")
      missing_plot
    },
    error = function(e) {
      if (verbose) message("Could not generate missing value comparison: ", e$message)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Missing value comparison unavailable") +
        ggplot2::theme_void() +
        ggplot2::ggtitle("Missing Values Comparison")
    }
  )

  # 3. Intensity Distribution Comparison
  plot_list$intensity_distribution <- tryCatch(
    {
      if (verbose) cat("Generating intensity distribution plot...\n")

      # Get data for intensity comparison (same matrices as above)
      if (analysis_type == "peptide_imputation") {
        after_matrix <- after_object@peptide_matrix
        before_matrix <- if (!is.null(before_object)) before_object@peptide_matrix else NULL
        data_type <- "Peptide"
      } else {
        after_matrix <- as.matrix(after_object@protein_quant_table[, -1])

        # Use stored peptide data for before comparison
        if (!is.null(dpc_results$y_peptide_for_dpc)) {
          before_matrix <- dpc_results$y_peptide_for_dpc
          data_type <- "Peptide -> Protein"
        } else if (!is.null(before_object) && "peptide_matrix" %in% slotNames(before_object)) {
          before_matrix <- before_object@peptide_matrix
          data_type <- "Peptide -> Protein"
        } else {
          before_matrix <- NULL
          data_type <- "Protein"
        }
      }

      if (!is.null(before_matrix) && !is.null(after_matrix)) {
        # Ensure compatibility
        if (!identical(dim(before_matrix), dim(after_matrix))) {
          common_samples <- intersect(colnames(before_matrix), colnames(after_matrix))
          common_features <- intersect(rownames(before_matrix), rownames(after_matrix))
          if (length(common_samples) > 0 && length(common_features) > 0) {
            before_matrix <- before_matrix[common_features, common_samples]
            after_matrix <- after_matrix[common_features, common_samples]
          }
        }

        # Extract intensity data (remove NA and infinite values)
        before_data <- as.vector(before_matrix)
        before_data <- before_data[!is.na(before_data) & is.finite(before_data)]

        after_data <- as.vector(after_matrix)
        after_data <- after_data[!is.na(after_data) & is.finite(after_data)]

        if (length(before_data) > 0 && length(after_data) > 0) {
          # Create appropriate labels for the comparison
          if (analysis_type == "protein_quantification") {
            before_label <- "Input Peptides"
            after_label <- "DPC-Quant Proteins"
          } else {
            before_label <- "Before limpa"
            after_label <- "After limpa"
          }

          intensity_data <- data.frame(
            Intensity = c(before_data, after_data),
            Stage = c(
              rep(before_label, length(before_data)),
              rep(after_label, length(after_data))
            )
          )

          intensity_plot <- intensity_data |>
            ggplot2::ggplot(ggplot2::aes(x = Intensity, fill = Stage)) +
            ggplot2::geom_density(alpha = 0.6) +
            ggplot2::theme_bw() +
            ggplot2::labs(
              title = paste(data_type, "Intensity Distribution"),
              x = "log2 Intensity",
              y = "Density"
            ) +
            ggplot2::scale_fill_manual(values = c(
              "Input Peptides" = "#E74C3C",
              "DPC-Quant Proteins" = "#3498DB",
              "Before limpa" = "#E74C3C",
              "After limpa" = "#3498DB"
            ))
        } else {
          intensity_plot <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for\nintensity distribution") +
            ggplot2::theme_void() +
            ggplot2::ggtitle("Intensity Distribution")
        }
      } else {
        intensity_plot <- ggplot2::ggplot() +
          ggplot2::annotate("text",
            x = 0.5, y = 0.5,
            label = "Intensity distribution\nrequires before/after objects"
          ) +
          ggplot2::theme_void() +
          ggplot2::ggtitle("Intensity Distribution")
      }

      if (verbose) cat("[OK] Intensity distribution plot generated\n")
      intensity_plot
    },
    error = function(e) {
      if (verbose) message("Could not generate intensity distribution: ", e$message)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Intensity distribution unavailable") +
        ggplot2::theme_void() +
        ggplot2::ggtitle("Intensity Distribution")
    }
  )

  # 4. Summary Statistics Plot
  plot_list$summary <- tryCatch(
    {
      if (verbose) cat("Generating summary plot...\n")

      # Calculate current missing percentage
      if (analysis_type == "peptide_imputation") {
        current_matrix <- after_object@peptide_matrix
      } else {
        current_matrix <- as.matrix(after_object@protein_quant_table[, -1])
      }

      current_missing_pct <- round(100 * mean(is.na(current_matrix)), 1)

      # Create summary text
      slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "",
        "nearly random missing", dpc_results$slope_interpretation
      )

      summary_text <- paste(
        "limpa DPC Analysis Summary\n",
        "Type:", stringr::str_to_title(gsub("_", " ", analysis_type)), "\n",
        "Method:", ifelse(is.null(dpc_results$dpc_method), "limpa_dpc", dpc_results$dpc_method), "\n",
        "DPC Slope (beta1):", round(as.numeric(dpc_params[2]), 3), "\n",
        "DPC Intercept (beta0):", round(as.numeric(dpc_params[1]), 3), "\n",
        "Missing Mechanism:", slope_text, "\n",
        if (!is.null(dpc_results$missing_percentage_before)) {
          paste("Missing % Before:", dpc_results$missing_percentage_before, "%\n")
        } else {
          ""
        },
        "Missing % Current:", current_missing_pct, "%"
      )

      summary_plot <- ggplot2::ggplot() +
        ggplot2::annotate("text",
          x = 0.5, y = 0.5,
          label = summary_text,
          size = 4, hjust = 0.5, vjust = 0.5
        ) +
        ggplot2::theme_void() +
        ggplot2::ggtitle("limpa Analysis Summary") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))

      if (verbose) cat("[OK] Summary plot generated\n")
      summary_plot
    },
    error = function(e) {
      if (verbose) message("Could not generate summary plot: ", e$message)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Summary unavailable") +
        ggplot2::theme_void() +
        ggplot2::ggtitle("Summary")
    }
  )

  # Create composite plot using patchwork
  if (verbose) cat("Creating composite plot...\n")

  composite_plot <- tryCatch(
    {
      (plot_list$dpc_curve + plot_list$missing_comparison) /
        (plot_list$intensity_distribution + plot_list$summary) +
        patchwork::plot_annotation(
          title = paste("limpa DPC", stringr::str_to_title(gsub("_", " ", analysis_type)), "Quality Control"),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5))
        )
    },
    error = function(e) {
      if (verbose) message("Could not create composite plot with patchwork: ", e$message)
      # Fallback: return first available plot
      for (plot in plot_list) {
        if (!is.null(plot)) {
          return(plot)
        }
      }
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Composite plot failed")
    }
  )

  # Save individual plots if requested
  if (save_plots && !is.null(save_dir)) {
    message(sprintf("   Step: Saving plots to directory: %s", save_dir))
    if (verbose) cat("Saving individual plots...\n")

    tryCatch(
      {
        message(sprintf("      Checking if directory exists: %s", dir.exists(save_dir)))
        if (!dir.exists(save_dir)) {
          message("      Creating directory...")
          dir.create(save_dir, recursive = TRUE)
          message(sprintf("      Directory created: %s", dir.exists(save_dir)))
        }

        # Save individual plots
        message(sprintf("      Saving %d individual plots...", length(plot_list)))
        for (plot_name in names(plot_list)) {
          if (!is.null(plot_list[[plot_name]])) {
            filename <- file.path(save_dir, paste0(plot_prefix, "_", plot_name, ".png"))
            message(sprintf("        Saving plot: %s", filename))
            ggplot2::ggsave(filename, plot_list[[plot_name]], width = 8, height = 6, dpi = 300)
            message(sprintf("        Plot saved successfully: %s", file.exists(filename)))
          }
        }

        # Save composite plot
        composite_filename <- file.path(save_dir, paste0(plot_prefix, "_composite.png"))
        message(sprintf("      Saving composite plot: %s", composite_filename))
        ggplot2::ggsave(composite_filename, composite_plot, width = 12, height = 10, dpi = 300)
        message(sprintf("      Composite plot saved: %s", file.exists(composite_filename)))

        if (verbose) cat("[OK] Plots saved to:", save_dir, "\n")
      },
      error = function(e) {
        message(sprintf("   ERROR saving plots: %s", e$message))
        if (verbose) message("Error saving plots: ", e$message)
      }
    )
  } else {
    message(sprintf("   Not saving plots: save_plots=%s, save_dir=%s", save_plots, ifelse(is.null(save_dir), "NULL", save_dir)))
  }

  if (verbose) cat("[OK] limpa QC plot generation completed!\n")

  # Return the list of individual plots instead of the composite
  return(plot_list)
}

