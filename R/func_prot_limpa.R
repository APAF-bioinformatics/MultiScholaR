# ============================================================================
# func_prot_limpa.R
# ============================================================================
# Purpose: Protein-level limpa (DPC) imputation helper functions
# 
# This file contains standalone helper functions for limpa-based protein
# quantification and imputation workflows:
# - QC diagnostic plot generation for limpa analysis
# - Format conversion for dpcDE results to standard DE format
#
# Note: S4 methods for proteinMissingValueImputationLimpa remain in
# limpa_functions.R due to R package loading order requirements.
#
# Source: Extracted from R/limpa_functions.R
#
# Dependencies:
# - limpa, ggplot2, patchwork, dplyr
# ============================================================================

# Functions will be extracted here by the extraction tool from:
# - limpa_functions.R



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
#' @export
generateLimpaQCPlots <- function(after_object,
                                before_object = NULL,
                                save_plots = TRUE,
                                save_dir = "peptide_qc",
                                plot_prefix = "limpa",
                                verbose = TRUE) {
  
  message("--- Entering generateLimpaQCPlots ---")
  message(sprintf("   Arguments: after_object type=%s, before_object=%s, save_plots=%s, save_dir=%s, plot_prefix=%s, verbose=%s",
                  class(after_object), ifelse(is.null(before_object), "NULL", class(before_object)),
                  save_plots, save_dir, plot_prefix, verbose))
  
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
    tryCatch({
      # Access project_dirs from global environment
      if (exists("project_dirs", envir = .GlobalEnv)) {
        project_dirs <- get("project_dirs", envir = .GlobalEnv)
        
        # Look through each experiment in project_dirs
        for (experiment_name in names(project_dirs)) {
          if (!is.null(project_dirs[[experiment_name]]$peptide_qc_dir)) {
            save_dir <- project_dirs[[experiment_name]]$peptide_qc_dir
            if (verbose) cat("âœ“ Using peptide_qc_dir from global project_dirs:", save_dir, "\n")
            break
          }
        }
      }
    }, error = function(e) {
      if (verbose) cat("Could not access project_dirs from global environment, using default\n")
    })
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
  plot_list$dpc_curve <- tryCatch({
    message("--- Entering DPC curve plot generation ---")
    if (verbose) cat("Generating DPC curve plot...\n")
    
    # Validate dpc_params first
    message("   Step: Validating DPC parameters...")
    message(sprintf("      dpc_params: %s", capture.output(str(dpc_params))))
    if (is.null(dpc_params) || !is.numeric(dpc_params) || length(dpc_params) < 2) {
      message("   Warning: Invalid DPC parameters detected, using defaults")
      if (verbose) cat("Warning: Invalid DPC parameters, using defaults\n")
      dpc_params <- c(-1.335, 0.074)  # Use the known values from your data
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
      tryCatch({
        message("        Checking if ggplotify package is available...")
        if (!requireNamespace("ggplotify", quietly = TRUE)) {
          message("        ggplotify not available, skipping this option")
        } else {
          message("        Calling plotDPC and converting with ggplotify")
          # Call plotDPC and immediately convert to ggplot
          ggplotify_plot <- ggplotify::as.ggplot(function() limpa::plotDPC(stored_dpc))
          message(sprintf("        ggplotify conversion successful, class: %s", class(ggplotify_plot)))
        }
      }, error = function(e) {
        message(sprintf("        ggplotify conversion FAILED: %s", e$message))
        ggplotify_plot <<- NULL
      })
      
      # OPTION 2: Try cowplot approach
      message("   Attempting Option 2: cowplot conversion...")
      cowplot_plot <- NULL
      tryCatch({
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
      }, error = function(e) {
        message(sprintf("        cowplot conversion FAILED: %s", e$message))
        cowplot_plot <<- NULL
      })
      
      # OPTION 3: Try direct plotDPC call and see what happens
      message("   Attempting Option 3: Direct plotDPC call...")
      direct_plot <- NULL
      tryCatch({
        message("        About to call plotDPC directly")
        direct_plot <- limpa::plotDPC(stored_dpc)
        message(sprintf("        Direct plotDPC returned: %s", capture.output(str(direct_plot))))
        message(sprintf("        Direct plotDPC class: %s", class(direct_plot)))
      }, error = function(e) {
        message(sprintf("        Direct plotDPC FAILED: %s", e$message))
        direct_plot <- NULL
      })
      
      # OPTION 4: Try temp file approach  
      message("   Attempting Option 4: Temporary file approach...")
      temp_plot <- NULL
      tryCatch({
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
              ggplot2::scale_x_continuous(expand = c(0,0)) +
              ggplot2::scale_y_continuous(expand = c(0,0)) +
              ggplot2::theme_void()
          }
        } else {
          # Use magick to read and convert
          img <- magick::image_read(temp_file)
          temp_plot <- ggplot2::ggplot() + 
            ggplot2::annotation_raster(as.raster(img), xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
            ggplot2::scale_x_continuous(expand = c(0,0)) +
            ggplot2::scale_y_continuous(expand = c(0,0)) +
            ggplot2::theme_void()
        }
        
        message("        Removing temporary file")
        unlink(temp_file)
        
        message(sprintf("        temp file approach result: %s", capture.output(str(temp_plot))))
      }, error = function(e) {
        message(sprintf("        temp file approach FAILED: %s", e$message))
        try(dev.off(), silent = TRUE)
        temp_plot <- NULL
      })
      
      # Determine which plot was successful and assign to dpc_final_plot
      dpc_final_plot <- NULL
      if (!is.null(ggplotify_plot)) {
        message("   Option 1 SUCCESS: Adding title to ggplotify plot")
        slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "", 
                             "nearly random missing", dpc_results$slope_interpretation)
        
        # Add title using ggplot2 
        dpc_final_plot <- ggplotify_plot +
          ggplot2::labs(title = paste("Detection Probability Curve\n", 
                                      "Slope:", round(as.numeric(dpc_params[2]), 3),
                                      "- Mechanism:", slope_text)) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        
        message("     ggplotify plot with title created successfully")
      } else if (!is.null(cowplot_plot)) {
        message("   Option 2 SUCCESS: Adding title to cowplot plot")
        slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "", 
                             "nearly random missing", dpc_results$slope_interpretation)
        
        # Add title using ggplot2 
        dpc_final_plot <- cowplot_plot +
          ggplot2::labs(title = paste("Detection Probability Curve\n", 
                                      "Slope:", round(as.numeric(dpc_params[2]), 3),
                                      "- Mechanism:", slope_text)) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        
        message("     cowplot plot with title created successfully")
      } else if (!is.null(direct_plot) && inherits(direct_plot, "ggplot")) {
        message("   Option 3 SUCCESS: Using direct ggplot")
        slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "", 
                             "nearly random missing", dpc_results$slope_interpretation)
        
        dpc_final_plot <- direct_plot +
          ggplot2::labs(title = paste("Detection Probability Curve\n", 
                                      "Slope:", round(as.numeric(dpc_params[2]), 3),
                                      "- Mechanism:", slope_text)) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        
        message("     direct ggplot with title created successfully")
      } else if (!is.null(temp_plot)) {
        message("   Option 4 SUCCESS: Using temp file plot")
        slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "", 
                             "nearly random missing", dpc_results$slope_interpretation)
        
        dpc_final_plot <- temp_plot +
           ggplot2::labs(title = paste("Detection Probability Curve\n", 
                                       "Slope:", round(as.numeric(dpc_params[2]), 3),
                                       "- Mechanism:", slope_text)) +
           ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        
        message("     temp file plot with title created successfully")
      }
      
    } else {
      message("   No stored DPC object available - creating fallback plot")
      dpc_final_plot <- NULL
    }
    
    # Use the successful plot or create fallback
    if (!is.null(dpc_final_plot)) {
      if (verbose) cat("âœ“ DPC curve plot generated successfully\n")
      dpc_final_plot
    } else {
      message("   Fallback: Creating parameters summary plot")
      if (verbose) cat("Creating DPC parameters summary plot...\n")
      
      slope_text <- ifelse(is.null(dpc_results$slope_interpretation) || dpc_results$slope_interpretation == "", 
                           "nearly random missing", dpc_results$slope_interpretation)
      
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, 
                         label = paste("DPC Parameters Summary\n",
                                       "Î²0 (intercept):", round(as.numeric(dpc_params[1]), 3), "\n",
                                       "Î²1 (slope):", round(as.numeric(dpc_params[2]), 3), "\n",
                                       "Mechanism:", slope_text),
                         size = 4, hjust = 0.5) +
        ggplot2::theme_void() +
        ggplot2::ggtitle("DPC Parameters")
    }
    
  }, error = function(e) {
    message(sprintf("--- ERROR in DPC plot generation: %s ---", e$message))
    message(sprintf("   Error class: %s", class(e)))
    message(sprintf("   Call stack: %s", capture.output(traceback())))
    ggplot2::ggplot() + 
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste("DPC plot error:", e$message)) +
      ggplot2::theme_void() +
      ggplot2::ggtitle("DPC Curve")
  })
  
  # 2. Missing Value Pattern Comparison
  plot_list$missing_comparison <- tryCatch({
    if (verbose) cat("Generating missing value comparison plot...\n")
    
         # Get matrices for comparison
     if (analysis_type == "peptide_imputation") {
       after_matrix <- after_object@peptide_matrix
       before_matrix <- if (!is.null(before_object)) before_object@peptide_matrix else NULL
       data_type <- "Peptides"
     } else {
       # For protein quantification, show peptide input vs protein output
       after_matrix <- as.matrix(after_object@protein_quant_table[, -1])  # Exclude protein ID column
       
       # For protein quantification, we want to show the reduction in missing values
       # from input peptides to output proteins
       if (!is.null(dpc_results$y_peptide_for_dpc)) {
         before_matrix <- dpc_results$y_peptide_for_dpc  # Use stored peptide data
         data_type <- "Peptides â†’ Proteins"
       } else if (!is.null(before_object) && "peptide_matrix" %in% slotNames(before_object)) {
         before_matrix <- before_object@peptide_matrix
         data_type <- "Peptides â†’ Proteins"
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
        ggplot2::scale_fill_manual(values = c("Before limpa" = "#E74C3C", 
                                             "After limpa" = "#3498DB"))
      
    } else {
      missing_plot <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, 
                         label = "Missing value comparison\nrequires before/after objects") +
        ggplot2::theme_void() +
        ggplot2::ggtitle("Missing Values Comparison")
    }
    
    if (verbose) cat("âœ“ Missing value comparison plot generated\n")
    missing_plot
    
  }, error = function(e) {
    if (verbose) message("Could not generate missing value comparison: ", e$message)
    ggplot2::ggplot() + 
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Missing value comparison unavailable") +
      ggplot2::theme_void() +
      ggplot2::ggtitle("Missing Values Comparison")
  })
  
  # 3. Intensity Distribution Comparison
  plot_list$intensity_distribution <- tryCatch({
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
         data_type <- "Peptideâ†’Protein"
       } else if (!is.null(before_object) && "peptide_matrix" %in% slotNames(before_object)) {
         before_matrix <- before_object@peptide_matrix
         data_type <- "Peptideâ†’Protein"
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
           Stage = c(rep(before_label, length(before_data)),
                     rep(after_label, length(after_data)))
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
           ggplot2::scale_fill_manual(values = c("Input Peptides" = "#E74C3C",
                                                "DPC-Quant Proteins" = "#3498DB",
                                                "Before limpa" = "#E74C3C", 
                                                "After limpa" = "#3498DB"))
      } else {
        intensity_plot <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for\nintensity distribution") +
          ggplot2::theme_void() +
          ggplot2::ggtitle("Intensity Distribution")
      }
    } else {
      intensity_plot <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, 
                         label = "Intensity distribution\nrequires before/after objects") +
        ggplot2::theme_void() +
        ggplot2::ggtitle("Intensity Distribution")
    }
    
    if (verbose) cat("âœ“ Intensity distribution plot generated\n")
    intensity_plot
    
  }, error = function(e) {
    if (verbose) message("Could not generate intensity distribution: ", e$message)
    ggplot2::ggplot() + 
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Intensity distribution unavailable") +
      ggplot2::theme_void() +
      ggplot2::ggtitle("Intensity Distribution")
  })
  
  # 4. Summary Statistics Plot
  plot_list$summary <- tryCatch({
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
                         "nearly random missing", dpc_results$slope_interpretation)
     
     summary_text <- paste(
       "limpa DPC Analysis Summary\n",
       "Type:", stringr::str_to_title(gsub("_", " ", analysis_type)), "\n",
       "Method:", ifelse(is.null(dpc_results$dpc_method), "limpa_dpc", dpc_results$dpc_method), "\n",
       "DPC Slope (Î²1):", round(as.numeric(dpc_params[2]), 3), "\n",
       "DPC Intercept (Î²0):", round(as.numeric(dpc_params[1]), 3), "\n",
       "Missing Mechanism:", slope_text, "\n",
       if (!is.null(dpc_results$missing_percentage_before)) {
         paste("Missing % Before:", dpc_results$missing_percentage_before, "%\n")
       } else "",
       "Missing % Current:", current_missing_pct, "%"
     )
    
    summary_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = summary_text,
                       size = 4, hjust = 0.5, vjust = 0.5) +
      ggplot2::theme_void() +
      ggplot2::ggtitle("limpa Analysis Summary") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
    
    if (verbose) cat("âœ“ Summary plot generated\n")
    summary_plot
    
  }, error = function(e) {
    if (verbose) message("Could not generate summary plot: ", e$message)
    ggplot2::ggplot() + 
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Summary unavailable") +
      ggplot2::theme_void() +
      ggplot2::ggtitle("Summary")
  })
  
  # Create composite plot using patchwork
  if (verbose) cat("Creating composite plot...\n")
  
  composite_plot <- tryCatch({
    (plot_list$dpc_curve + plot_list$missing_comparison) /
    (plot_list$intensity_distribution + plot_list$summary) +
    patchwork::plot_annotation(
      title = paste("limpa DPC", stringr::str_to_title(gsub("_", " ", analysis_type)), "Quality Control"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5))
    )
  }, error = function(e) {
    if (verbose) message("Could not create composite plot with patchwork: ", e$message)
    # Fallback: return first available plot
    for (plot in plot_list) {
      if (!is.null(plot)) return(plot)
    }
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Composite plot failed")
  })
  
  # Save individual plots if requested
  if (save_plots && !is.null(save_dir)) {
    message(sprintf("   Step: Saving plots to directory: %s", save_dir))
    if (verbose) cat("Saving individual plots...\n")
    
    tryCatch({
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
      
      if (verbose) cat("âœ“ Plots saved to:", save_dir, "\n")
      
    }, error = function(e) {
      message(sprintf("   ERROR saving plots: %s", e$message))
      if (verbose) message("Error saving plots: ", e$message)
    })
  } else {
    message(sprintf("   Not saving plots: save_plots=%s, save_dir=%s", save_plots, ifelse(is.null(save_dir), "NULL", save_dir)))
  }
  
  if (verbose) cat("âœ“ limpa QC plot generation completed!\n")
  
  # Return the list of individual plots instead of the composite
  return(plot_list)
}


# ----------------------------------------------------------------------------
# convertDpcDEToStandardFormat
# ----------------------------------------------------------------------------
#' Convert limpa dpcDE Results to Standard Format
#'
#' This function converts the output from limpa::dpcDE() to the same format as
#' runTestsContrasts() to ensure compatibility with existing DE analysis workflows.
#'
#' @param dpc_fit MArrayLM object from limpa::dpcDE()
#' @param contrast_strings Character vector of contrast strings in "name=expression" format
#' @param design_matrix Design matrix used in the analysis
#' @param eBayes_trend Logical, whether empirical Bayes trend was used
#' @param eBayes_robust Logical, whether robust empirical Bayes was used
#'
#' @return List with same structure as runTestsContrasts output
#' @export
convertDpcDEToStandardFormat <- function(dpc_fit, 
                                          contrast_strings, 
                                          design_matrix,
                                          eBayes_trend = TRUE,
                                          eBayes_robust = TRUE) {
  
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma package is required for convertDpcDEToStandardFormat")
  }
  
  cat("   convertDpcDEToStandardFormat: Converting dpcDE results to standard format\n")
  cat("   convertDpcDEToStandardFormat: Processing", length(contrast_strings), "contrasts\n")
  
  # Create contrast matrix from contrast strings
  # Extract just the contrast expressions (after the "=" sign)
  contrast_expressions <- sapply(contrast_strings, function(x) {
    parts <- strsplit(x, "=")[[1]]
    if (length(parts) == 2) {
      return(parts[2])
    } else {
      return(x)  # Fallback if no "=" found
    }
  })
  
  # Create contrast matrix using limma
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_expressions, 
                                          levels = colnames(design_matrix))
  
  cat("   convertDpcDEToStandardFormat: Contrast matrix dims:", nrow(contrast_matrix), "x", ncol(contrast_matrix), "\n")
  
  # Apply contrasts to the dpcDE fit
  contrast_fit <- limma::contrasts.fit(dpc_fit, contrast_matrix)
  
  # Apply empirical Bayes
  eb_fit <- limma::eBayes(contrast_fit, trend = eBayes_trend, robust = eBayes_robust)
  
  # Extract results for each contrast
  results_list <- list()
  
  for (i in seq_along(contrast_strings)) {
    contrast_name <- names(contrast_strings)[i]
    if (is.null(contrast_name) || contrast_name == "") {
      # Extract friendly name from contrast string if no name provided
      parts <- strsplit(contrast_strings[i], "=")[[1]]
      if (length(parts) == 2) {
        contrast_name <- parts[1]
      } else {
        contrast_name <- paste0("Contrast_", i)
      }
    }
    
    # Get results for this contrast using limma::topTable
    contrast_results <- limma::topTable(
      eb_fit, 
      coef = i, 
      number = Inf, 
      sort.by = "none",  # Keep original order
      adjust.method = "BH"
    )
    
    # Rename columns to match runTestsContrasts output format, and add qvalue
    if (!"P.Value" %in% colnames(contrast_results)) {
        stop("P.Value column not found in topTable results for dpcDE.")
    }

    contrast_results <- contrast_results |>
      dplyr::mutate(
        # Ensure the qvalue library is used for fdr_qvalue, consistent with runTestsContrasts
        fdr_qvalue = qvalue::qvalue(P.Value)$qvalues,
        # Keep the BH adjustment in its own column for full compatibility
        fdr_value_bh_adjustment = adj.P.Val
      ) |>
      dplyr::rename(
        raw_pvalue = P.Value
      ) |>
      dplyr::mutate(
        comparison = contrast_name,
        uniprot_acc = rownames(contrast_results),
        log_intensity = AveExpr
      ) |>
      dplyr::select(
        uniprot_acc, 
        comparison, 
        logFC, 
        log_intensity,
        raw_pvalue, 
        fdr_qvalue, 
        fdr_value_bh_adjustment,
        everything(),
        -adj.P.Val # remove original adj.P.Val to avoid confusion
      )
    
    # Store in results list (matching runTestsContrasts structure)
    # The name of the list element MUST be the full contrast string for downstream parsing
    full_contrast_string <- contrast_strings[i]
    results_list[[full_contrast_string]] <- contrast_results
    
    cat("   convertDpcDEToStandardFormat: Processed contrast", contrast_name, "with", nrow(contrast_results), "proteins\n")
  }
  
  # Return in same format as runTestsContrasts, which is a list of tables
  return_object <- list(
    results = results_list, # Return the list of data frames
    fit.eb = eb_fit,
    dpc_method_used = TRUE
  )
  
  cat("   convertDpcDEToStandardFormat: Conversion completed successfully\n")
  
  return(return_object)
}


# ==========================================
# S4 Methods from limpa_functions.R
# ==========================================
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Protein-Level Missing Value Imputation using limpa Package
#'
#' This function applies limpa's DPC-based missing value imputation directly to 
#' protein-level quantification data. This is useful when you already have protein
#' quantification data with missing values that need to be handled.
#'
#' @param theObject A ProteinQuantitativeData object with protein-level data
#' @param dpc_results DPC results to use. If NULL, will estimate using dpc_slope
#' @param dpc_slope Default DPC slope to use if no DPC results available (default: 0.8)
#' @param quantified_protein_column Name for the column containing quantified protein values
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param chunk When verbose=TRUE, how often to output progress information (default: 1000)
#'
#' @details
#' This method treats each protein as a separate "feature" and applies DPC-based
#' imputation using limpa's dpcImpute function. This is appropriate when you have
#' protein-level data with missing values that follow intensity-dependent patterns.
#'
#' @return Updated ProteinQuantitativeData object with imputed protein values
#'
#' @importFrom limpa dpc dpcImpute
#' @export
setMethod(f="proteinMissingValueImputationLimpa"
          , signature="ProteinQuantitativeData"
          , definition = function(theObject, 
                                  dpc_results = NULL,
                                  dpc_slope = 0.8,
                                  quantified_protein_column = NULL,
                                  verbose = TRUE,
                                  chunk = 1000) {
            
            # Load required packages
            if (!requireNamespace("limpa", quietly = TRUE)) {
              stop("limpa package is required but not installed. Please install it using: BiocManager::install('limpa')")
            }
            
            # Parameter validation and defaults
            quantified_protein_column <- if (is.null(quantified_protein_column)) {
              "Protein.Imputed.Limpa"
            } else {
              quantified_protein_column
            }
            
            # Extract data from protein object
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            sample_id_column <- theObject@sample_id
            design_matrix <- theObject@design_matrix
            
            if (verbose) {
              log_info("Starting limpa-based protein-level missing value imputation...")
            }
            
            # Convert to matrix format (proteins as rows, samples as columns)
            # First, identify sample columns (exclude protein ID and other metadata)
            sample_columns <- setdiff(colnames(protein_quant_table), 
                                     c(protein_id_column, "description", "gene_name", 
                                       "protein_name", "organism", "length"))
            
            if (verbose) {
              log_info("Converting protein data to matrix format...")
              log_info("Found {length(sample_columns)} sample columns")
            }
            
            # Create protein matrix
            protein_matrix <- protein_quant_table |>
              dplyr::select(all_of(c(protein_id_column, sample_columns))) |>
              tibble::column_to_rownames(protein_id_column) |>
              as.matrix()
            
            if (verbose) {
              log_info("Protein matrix dimensions: {nrow(protein_matrix)} proteins x {ncol(protein_matrix)} samples")
              log_info("Missing value percentage: {round(100 * mean(is.na(protein_matrix)), 1)}%")
            }
            
            # Check if we need log2 transformation
            # Assume if max value > 50, data is not log2 transformed
            max_val <- max(protein_matrix, na.rm = TRUE)
            needs_log_transform <- max_val > 50
            
            if (needs_log_transform) {
              if (verbose) {
                log_info("Converting to log2 scale for limpa (max value: {round(max_val, 2)})...")
              }
              protein_matrix <- log2(protein_matrix + 1)
            } else {
              if (verbose) {
                log_info("Data appears to be log2-scale already (max value: {round(max_val, 2)})")
              }
            }
            
            # Handle infinite or NaN values
            if (any(is.infinite(protein_matrix) | is.nan(protein_matrix), na.rm = TRUE)) {
              if (verbose) {
                log_warn("Infinite or NaN values detected. Replacing with NA...")
              }
              protein_matrix[is.infinite(protein_matrix) | is.nan(protein_matrix)] <- NA
            }
            
            # Get or estimate DPC parameters
            dpc_params <- NULL
            if (!is.null(dpc_results)) {
              dpc_params <- dpc_results
              if (verbose) {
                log_info("Using provided DPC results")
              }
            } else {
              if (verbose) {
                log_info("Estimating DPC from protein data...")
              }
              tryCatch({
                dpcfit <- limpa::dpc(protein_matrix)
                dpc_params <- dpcfit
                if (verbose) {
                  log_info("DPC parameters estimated:")
                  log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                  log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
                }
              }, error = function(e) {
                if (verbose) {
                  log_warn("DPC estimation failed, using default slope: {dpc_slope}")
                }
                dpc_params <- NULL
              })
            }
            
            # Apply missing value imputation
            if (verbose) {
              log_info("Applying protein-level missing value imputation...")
            }
            
            tryCatch({
              # Apply dpcImpute to protein matrix
              if (!is.null(dpc_params)) {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc = dpc_params, verbose = verbose, chunk = chunk)
              } else {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc.slope = dpc_slope, verbose = verbose, chunk = chunk)
              }
              
              if (verbose) {
                log_info("Protein-level imputation completed successfully")
                log_info("No missing values remaining: {!any(is.na(imputed_result$E))}")
              }
              
              # Extract imputed matrix
              imputed_matrix <- imputed_result$E
              
              # Transform back to original scale if necessary
              if (needs_log_transform) {
                if (verbose) {
                  log_info("Converting back from log2 scale...")
                }
                imputed_matrix <- 2^imputed_matrix - 1
                # Ensure no negative values
                imputed_matrix[imputed_matrix < 0] <- 0
              }
              
              # Convert back to long format
              if (verbose) {
                log_info("Converting imputed data back to original format...")
              }
              
              imputed_long <- imputed_matrix |>
                as.data.frame() |>
                tibble::rownames_to_column(protein_id_column) |>
                tidyr::pivot_longer(cols = -all_of(protein_id_column), 
                                   names_to = sample_id_column, 
                                   values_to = quantified_protein_column)
              
              # Merge with original protein data
              updated_protein_data <- protein_quant_table |>
                dplyr::left_join(imputed_long, by = c(protein_id_column, sample_id_column))
              
              # Update the object
              theObject@protein_quant_table <- updated_protein_data
              
              # Store DPC results for future reference
              if (is.null(theObject@args)) {
                theObject@args <- list()
              }
              
              theObject@args$limpa_protein_imputation_results <- list(
                dpc_parameters_used = if (!is.null(dpc_params)) {
                  if (is.list(dpc_params) && !is.null(dpc_params$dpc)) {
                    dpc_params$dpc  # Extract parameters from DPC object
                  } else if (is.numeric(dpc_params)) {
                    dpc_params  # Already numeric parameters  
                  } else {
                    c(NA, dpc_slope)
                  }
                } else {
                  c(NA, dpc_slope)
                },
                dpc_object_used = if (is.list(dpc_params) && !is.null(dpc_params$dpc)) dpc_params else NULL,
                quantified_protein_column = quantified_protein_column,
                missing_percentage_before = round(100 * mean(is.na(protein_matrix)), 1),
                missing_percentage_after = 0,  # DPC imputation produces complete data
                imputation_method = "limpa_dpc_protein_imputation",
                total_proteins_imputed = nrow(imputed_matrix)
              )
              
              if (verbose) {
                log_info("limpa protein-level imputation completed successfully!")
                log_info("New imputed column: {quantified_protein_column}")
                log_info("DPC results stored in object@args$limpa_protein_imputation_results")
              }
              
              return(theObject)
              
            }, error = function(e) {
              log_error(paste("Error during limpa protein imputation:", e$message))
              stop(paste("limpa protein imputation failed:", e$message))
            })
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Protein Quantification using limpa Package DPC-Quant Method
#'
#' This function uses the limpa package's Detection Probability Curve (DPC) approach
#' combined with DPC-Quant for sophisticated protein-level quantification from peptide data.
#' This method leverages the DPC estimated at the peptide level to perform robust
#' protein quantification that accounts for missing value mechanisms.
#'
#' @param theObject A PeptideQuantitativeData object with peptide-level data
#' @param dpc_results DPC results from peptide-level analysis. If NULL, will use
#'   results stored in theObject@args$limpa_dpc_results. If those don't exist,
#'   will estimate DPC using default slope.
#' @param dpc_slope Default DPC slope to use if no DPC results available (default: 0.8)
#' @param quantified_protein_column Name for the new column containing quantified protein values.
#'   Default is "Protein.Quantified.Limpa"
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param chunk When verbose=TRUE, how often to output progress information (default: 1000)
#'
#' @details
#' The DPC-Quant method represents missing values probabilistically using the Detection
#' Probability Curve and returns maximum posterior estimates for all protein log2-expression
#' values. Unlike simple imputation, this method provides:
#' 
#' 1. Protein-level quantification from peptide data using hierarchical modeling
#' 2. Uncertainty quantification via standard errors for each protein estimate
#' 3. No missing values in final protein summaries
#' 4. Proper handling of peptides with sparse observations
#'
#' The process follows these steps:
#' 1. Extract or estimate Detection Probability Curve from peptide data
#' 2. Apply dpcQuant() to obtain protein-level summaries
#' 3. Extract quantified protein matrix and uncertainty estimates
#' 4. Create ProteinQuantitativeData object with results
#'
#' @return ProteinQuantitativeData object with quantified protein values and metadata
#'
#' @importFrom limpa dpcQuant
#' @export
setMethod(f="proteinMissingValueImputationLimpa"
          , signature="PeptideQuantitativeData"
          , definition = function(theObject, 
                                  dpc_results = NULL,
                                  dpc_slope = 0.8,
                                  quantified_protein_column = NULL,
                                  verbose = TRUE,
                                  chunk = 1000) {
            
            # Load required packages
            if (!requireNamespace("limpa", quietly = TRUE)) {
              stop("limpa package is required but not installed. Please install it using: BiocManager::install('limpa')")
            }
            
            # Parameter validation and defaults
            quantified_protein_column <- if (is.null(quantified_protein_column)) {
              "Protein.Quantified.Limpa"
            } else {
              quantified_protein_column
            }
            
            dpc_slope <- checkParamsObjectFunctionSimplify(theObject, "dpc_slope", dpc_slope)
            verbose <- checkParamsObjectFunctionSimplify(theObject, "verbose", verbose)
            chunk <- checkParamsObjectFunctionSimplify(theObject, "chunk", chunk)
            
            # Extract data from peptide object
            peptide_data <- theObject@peptide_data
            peptide_matrix <- theObject@peptide_matrix
            protein_id_column <- theObject@protein_id_column
            peptide_sequence_column <- theObject@peptide_sequence_column
            sample_id_column <- theObject@sample_id
            design_matrix <- theObject@design_matrix
            
            if (verbose) {
              log_info("Starting limpa DPC-Quant protein quantification...")
              log_info("Peptide data dimensions: {nrow(peptide_matrix)} peptides x {ncol(peptide_matrix)} samples")
              log_info("Missing value percentage in peptides: {round(100 * mean(is.na(peptide_matrix)), 1)}%")
            }
            
            # Check if peptide matrix exists and is properly formatted
            if (is.null(peptide_matrix) || length(peptide_matrix) == 0) {
              stop("Peptide matrix not found. Please run calcPeptideMatrix() first.")
            }
            
            # Prepare peptide data for limpa (peptides as rows, samples as columns)
            y_peptide <- peptide_matrix
            
            # Ensure data is log2 transformed for limpa
            if (!theObject@is_logged_data) {
              if (verbose) {
                log_info("Converting to log2 scale for limpa DPC-Quant...")
              }
              y_peptide <- log2(y_peptide + 1)
            } else {
              if (verbose) {
                log_info("Using existing log2 transformed peptide data")
              }
            }
            
            # Handle infinite or NaN values
            if (any(is.infinite(y_peptide) | is.nan(y_peptide), na.rm = TRUE)) {
              if (verbose) {
                log_warn("Infinite or NaN values detected. Replacing with NA...")
              }
              y_peptide[is.infinite(y_peptide) | is.nan(y_peptide)] <- NA
            }
            
            # NEW: Remove rows that are all NA, which cause limpa optimization to fail
            all_na_rows <- rowSums(is.na(y_peptide)) == ncol(y_peptide)
            if (any(all_na_rows)) {
              num_removed <- sum(all_na_rows)
              if (verbose) {
                log_warn("Found and removed {num_removed} peptides with no observations across all samples.")
              }
              y_peptide <- y_peptide[!all_na_rows, ]
            }
            
            # Get or estimate DPC parameters
            dpc_params <- NULL
            slope_interpretation <- "nearly random missing"  # Default interpretation
            if (!is.null(dpc_results)) {
              # Use provided DPC results
              dpc_params <- dpc_results
              if (verbose) {
                log_info("Using provided DPC results")
                if (is.numeric(dpc_results) && length(dpc_results) == 2) {
                  log_info("  beta0 (intercept): {round(dpc_results[1], 4)}")
                  log_info("  beta1 (slope): {round(dpc_results[2], 4)}")
                }
              }
            } else if (!is.null(theObject@args$limpa_dpc_results)) {
              # Prefer full DPC object if available, otherwise use parameters
              if (!is.null(theObject@args$limpa_dpc_results$dpc_object)) {
                dpc_params <- theObject@args$limpa_dpc_results$dpc_object
                # Extract slope interpretation if available
                if (!is.null(theObject@args$limpa_dpc_results$slope_interpretation)) {
                  slope_interpretation <- theObject@args$limpa_dpc_results$slope_interpretation
                }
                if (verbose) {
                  log_info("Using full DPC object from previous peptide analysis:")
                  log_info("  beta0 (intercept): {round(dpc_params$dpc[1], 4)}")
                  log_info("  beta1 (slope): {round(dpc_params$dpc[2], 4)}")
                  log_info("  Interpretation: {slope_interpretation}")
                }
              } else if (!is.null(theObject@args$limpa_dpc_results$dpc_parameters)) {
                dpc_params <- theObject@args$limpa_dpc_results$dpc_parameters
                # Extract slope interpretation if available
                if (!is.null(theObject@args$limpa_dpc_results$slope_interpretation)) {
                  slope_interpretation <- theObject@args$limpa_dpc_results$slope_interpretation
                }
                if (verbose) {
                  log_info("Using DPC parameters from previous peptide analysis:")
                  log_info("  beta0 (intercept): {round(dpc_params[1], 4)}")
                  log_info("  beta1 (slope): {round(dpc_params[2], 4)}")
                  log_info("  Interpretation: {slope_interpretation}")
                }
              }
            } else {
              # Estimate DPC parameters on-the-fly if not available
              if (verbose) {
                log_info("No DPC results found. Estimating DPC parameters from peptide data...")
              }
              tryCatch({
                # Estimate DPC using the peptide data
                dpcfit <- limpa::dpc(y_peptide)
                dpc_params <- dpcfit  # Use the full DPC object
                
                # Interpret the slope
                slope_interpretation <- if (dpcfit$dpc[2] < 0.3) {
                  "nearly random missing"
                } else if (dpcfit$dpc[2] < 0.7) {
                  "moderate intensity-dependent missing" 
                } else if (dpcfit$dpc[2] < 1.2) {
                  "strong intensity-dependent missing"
                } else {
                  "very strong intensity-dependent missing (approaching left-censoring)"
                }
                
                if (verbose) {
                  log_info("DPC parameters estimated on-the-fly:")
                  log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                  log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
                  log_info("  Interpretation: {slope_interpretation}")
                }
              }, error = function(e) {
                if (verbose) {
                  log_warn("DPC estimation failed: {e$message}. Using default slope: {dpc_slope}")
                }
                dpc_params <<- NULL  # Will let dpcQuant estimate with default slope
                slope_interpretation <<- "unable to determine (DPC estimation failed)"
              })
            }
            
            # Create protein ID mapping from peptide data
            # Extract protein IDs from rownames (format: "ProteinID%PeptideSequence")
            rownames_split <- strsplit(rownames(y_peptide), "%")
            protein_ids <- sapply(rownames_split, function(x) x[1])
            
            # Calculate peptide and peptidoform counts per protein for future filtering
            # This ensures we carry this critical info to the protein level
            peptide_summary <- peptide_data |>
              dplyr::group_by(!!sym(protein_id_column)) |>
              dplyr::summarise(
                peptide_count = dplyr::n_distinct(!!sym(peptide_sequence_column)),
                peptidoform_count = dplyr::n(), # Total peptidoforms is just the number of rows per protein
                .groups = "drop"
              )
            
            if (verbose) {
              unique_proteins <- nrow(peptide_summary)
              log_info("Found {unique_proteins} unique proteins from {nrow(y_peptide)} peptidoforms")
              
              peptide_count_summary <- peptide_summary |>
                dplyr::summarise(
                  min_peptides = min(peptide_count),
                  max_peptides = max(peptide_count),
                  median_peptides = median(peptide_count),
                  min_peptidoforms = min(peptidoform_count),
                  max_peptidoforms = max(peptidoform_count),
                  median_peptidoforms = median(peptidoform_count),
                  proteins_with_2plus_peptides = sum(peptide_count >= 2)
                )
              log_info("Unique peptide counts per protein - min: {peptide_count_summary$min_peptides}, max: {peptide_count_summary$max_peptides}, median: {round(peptide_count_summary$median_peptides, 1)}")
              log_info("Total peptidoform counts per protein - min: {peptide_count_summary$min_peptidoforms}, max: {peptide_count_summary$max_peptidoforms}, median: {round(peptide_count_summary$median_peptidoforms, 1)}")
              log_info("Proteins with >=2 unique peptides: {peptide_count_summary$proteins_with_2plus_peptides} / {unique_proteins}")
            }
            
            # Apply DPC-Quant for protein quantification
            if (verbose) {
              log_info("Applying DPC-Quant for protein quantification...")
            }
            
            tryCatch({
              # Create EList-like object for limpa
              y_elist <- list(
                E = y_peptide,
                genes = data.frame(
                  protein.id = protein_ids,
                  peptide.id = rownames(y_peptide),
                  stringsAsFactors = FALSE
                )
              )
              class(y_elist) <- "EList"
              
              # Apply dpcQuant
              if (!is.null(dpc_params)) {
                protein_quant_result <- limpa::dpcQuant(
                  y = y_elist,
                  protein.id = "protein.id",
                  dpc = dpc_params,
                  verbose = verbose,
                  chunk = chunk
                )
              } else {
                protein_quant_result <- limpa::dpcQuant(
                  y = y_elist,
                  protein.id = "protein.id",
                  dpc.slope = dpc_slope,
                  verbose = verbose,
                  chunk = chunk
                )
              }
              
              if (verbose) {
                log_info("DPC-Quant completed successfully")
                log_info("Quantified proteins: {nrow(protein_quant_result$E)}")
                log_info("Samples: {ncol(protein_quant_result$E)}")
                log_info("No missing values in protein quantification: {!any(is.na(protein_quant_result$E))}")
              }
              
              # Extract quantified protein matrix
              protein_matrix <- protein_quant_result$E
              rownames(protein_matrix) <- protein_quant_result$genes$protein.id
              colnames(protein_matrix) <- colnames(y_peptide)
              
              # Extract standard errors and observation counts
              standard_errors <- protein_quant_result$other$standard.error
              n_observations <- protein_quant_result$other$n.observations
              
              # Transform back to original scale if necessary
              if (!theObject@is_logged_data) {
                if (verbose) {
                  log_info("Converting back from log2 scale...")
                }
                protein_matrix <- 2^protein_matrix - 1
                # Ensure no negative values
                protein_matrix[protein_matrix < 0] <- 0
              }
              
              # Convert protein matrix to wide format data frame (samples as columns)
              if (verbose) {
                log_info("Converting protein data to wide format...")
              }
              
              protein_wide <- protein_matrix |>
                as.data.frame() |>
                tibble::rownames_to_column(protein_id_column)
              
              # --- NEW: Synchronize the peptide summary with the actual quantified proteins ---
              # dpcQuant may not be able to quantify all proteins present in the input.
              # We must ensure the summary table only contains proteins that are in the final matrix.
              final_protein_ids <- protein_wide[[protein_id_column]]
              peptide_summary_synced <- peptide_summary |>
                dplyr::filter(!!sym(protein_id_column) %in% final_protein_ids)
              
              if (verbose) {
                log_info("Synchronized peptide summary table with quantified proteins. Kept {nrow(peptide_summary_synced)} of {nrow(peptide_summary)} entries.")
              }
              # --- END NEW ---
              
              # Create ProteinQuantitativeData object
              if (verbose) {
                log_info("Creating ProteinQuantitativeData object...")
              }
              
              protein_obj <- ProteinQuantitativeData(
                protein_quant_table = protein_wide,
                protein_id_column = protein_id_column,
                protein_id_table = protein_wide |> dplyr::distinct(!!sym(protein_id_column)),
                design_matrix = design_matrix,
                sample_id = sample_id_column,
                group_id = theObject@group_id,
                technical_replicate_id = theObject@technical_replicate_id,
                args = theObject@args
              )
              
              # Store DPC-Quant results in args for future reference
              if (is.null(protein_obj@args)) {
                protein_obj@args <- list()
              }
              
              protein_obj@args$limpa_dpc_quant_results <- list(
                dpc_parameters_used = if (!is.null(dpc_params)) {
                  if (is.list(dpc_params) && !is.null(dpc_params$dpc)) {
                    dpc_params$dpc  # Extract parameters from DPC object
                  } else if (is.numeric(dpc_params)) {
                    dpc_params  # Already numeric parameters
                  } else {
                    c(NA, dpc_slope)
                  }
                } else {
                  c(NA, dpc_slope)
                },
                dpc_object_used = if (is.list(dpc_params) && !is.null(dpc_params$dpc)) dpc_params else NULL,
                slope_interpretation = slope_interpretation,  # Store slope interpretation for QC plots
                y_peptide_for_dpc = y_peptide,  # Store peptide data for DPC plotting
                quantified_elist = protein_quant_result, # Store the entire EList object
                standard_errors = standard_errors,
                n_observations = n_observations,
                peptide_counts_per_protein = peptide_summary_synced, # Use the new SYNCED summary table
                missing_percentage_peptides = round(100 * mean(is.na(y_peptide)), 1),
                missing_percentage_proteins = 0,  # DPC-Quant produces complete protein data
                quantification_method = "limpa_dpc_quant",
                total_proteins_quantified = nrow(protein_matrix),
                total_peptides_used = nrow(y_peptide),
                dpc_method = "limpa_dpc_quant"  # Identify the method type
              )
              
              if (verbose) {
                log_info("limpa DPC-Quant protein quantification completed successfully!")
                log_info("Protein data stored in wide format with {ncol(protein_matrix)} sample columns")
                log_info("Proteins quantified: {nrow(protein_matrix)}")
                log_info("DPC-Quant results stored in object@args$limpa_dpc_quant_results")
              }
              
              return(protein_obj)
              
            }, error = function(e) {
              log_error(paste("Error during limpa DPC-Quant:", e$message))
              stop(paste("limpa DPC-Quant failed:", e$message))
            })
          })
