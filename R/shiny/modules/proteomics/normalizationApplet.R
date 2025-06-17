#' @title normalizationAppletModule
#'
#' @description A Shiny module for the Normalization step of the proteomics
#' workflow. Automatically generates pre-normalization QC plots and handles
#' normalization methods and RUV-III batch correction.
#'
#' @name normalizationAppletModule
NULL

#' @rdname normalizationAppletModule
#' @export
#' @import shiny
#' @import shinydashboard
normalizationAppletUI <- function(id) {
  ns <- NS(id)
  
  # This is handled in the main workflow UI now
  # Return empty div as placeholder
  shiny::div()
}

#' @rdname normalizationAppletModule 
#' @export
normalizationAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Initialize reactive values for normalization state
    norm_data <- reactiveValues(
      pre_norm_qc_generated = FALSE,
      normalization_complete = FALSE,
      ruv_complete = FALSE,
      
      # QC plot objects - 3 columns: post-filtering, post-normalization, ruv-corrected
      qc_plots = list(
        post_filtering = list(
          pca = NULL,
          density = NULL, 
          rle = NULL,
          correlation = NULL
        ),
        post_normalization = list(
          pca = NULL,
          density = NULL,
          rle = NULL, 
          correlation = NULL
        ),
        ruv_corrected = list(
          pca = NULL,
          density = NULL,
          rle = NULL,
          correlation = NULL
        )
      ),
      
      # Normalization results
      normalized_protein_obj = NULL,
      ruv_normalized_obj = NULL,
      best_k = NULL,
      control_genes_index = NULL
    )
    
    # Helper function to get current plot aesthetics
    getPlotAesthetics <- function() {
      list(
        color_var = if(is.null(input$color_variable) || input$color_variable == "") "group" else input$color_variable,
        shape_var = if(is.null(input$shape_variable) || input$shape_variable == "") "group" else input$shape_variable
      )
    }
    
    # Helper function to generate pre-normalization QC plots
    generatePreNormalizationQc <- function() {
      message("Generating pre-normalization QC plots...")
      
      # Get current S4 object using CORRECT R6 methods
      shiny::req(workflow_data$state_manager)
      current_state <- workflow_data$state_manager$current_state
      current_s4 <- workflow_data$state_manager$getState(current_state)
      
      if (is.null(current_s4)) {
        stop("No S4 object available for QC plot generation")
      }
      
      # Get current plot aesthetics
      aesthetics <- getPlotAesthetics()
      
      # Generate QC plots and store in norm_data$qc_plots$post_filtering
      
      # PCA plot
      norm_data$qc_plots$post_filtering$pca <- plotPca(
        current_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var, 
        title = "Pre-Normalization PCA",
        font_size = 8
      )
      
      # RLE plot  
      norm_data$qc_plots$post_filtering$rle <- plotRle(
        current_s4,
        group = aesthetics$color_var,
        yaxis_limit = c(-6, 6)
      )
      
      # Density plot
      norm_data$qc_plots$post_filtering$density <- plotDensity(
        norm_data$qc_plots$post_filtering$pca,
        grouping_variable = aesthetics$color_var
      )
      
      # Correlation plot
      norm_data$qc_plots$post_filtering$correlation <- plotPearson(
        current_s4,
        tech_rep_remove_regex = "pool",
        correlation_group = aesthetics$color_var
      )
      
      message("Pre-normalization QC plots generated successfully")
    }
    
    # Helper function to generate post-normalization QC plots
    generatePostNormalizationQc <- function(normalized_s4) {
      aesthetics <- getPlotAesthetics()
      
      norm_data$qc_plots$post_normalization$pca <- plotPca(
        normalized_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
        title = "Post-Normalization PCA", 
        font_size = 8
      )
      
      norm_data$qc_plots$post_normalization$rle <- plotRle(
        normalized_s4,
        group = aesthetics$color_var,
        yaxis_limit = c(-6, 6)
      )
      
      norm_data$qc_plots$post_normalization$density <- plotDensity(
        norm_data$qc_plots$post_normalization$pca,
        grouping_variable = aesthetics$color_var
      )
      
      norm_data$qc_plots$post_normalization$correlation <- plotPearson(
        normalized_s4,
        tech_rep_remove_regex = "pool",
        correlation_group = aesthetics$color_var
      )
    }
    
    # Helper function to generate RUV-corrected QC plots
    generateRuvCorrectedQc <- function(ruv_corrected_s4) {
      aesthetics <- getPlotAesthetics()
      
      norm_data$qc_plots$ruv_corrected$pca <- plotPca(
        ruv_corrected_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
        title = "RUV-Corrected PCA",
        font_size = 8
      )
      
      norm_data$qc_plots$ruv_corrected$rle <- plotRle(
        ruv_corrected_s4,
        group = aesthetics$color_var, 
        yaxis_limit = c(-6, 6)
      )
      
      norm_data$qc_plots$ruv_corrected$density <- plotDensity(
        norm_data$qc_plots$ruv_corrected$pca,
        grouping_variable = aesthetics$color_var
      )
      
      norm_data$qc_plots$ruv_corrected$correlation <- plotPearson(
        ruv_corrected_s4,
        tech_rep_remove_regex = "pool",
        correlation_group = aesthetics$color_var
      )
    }
    
    # Update plot aesthetic choices based on design matrix
    observe({
      if (!is.null(workflow_data$design_matrix)) {
        design_cols <- colnames(workflow_data$design_matrix)
        
        # Filter to common experimental variables
        available_vars <- intersect(design_cols, c("group", "factor1", "factor2", "technical_replicate_id", "sample_id"))
        
        if (length(available_vars) > 0) {
          # Update color variable choices
          shiny::updateSelectInput(session, "color_variable",
            choices = available_vars,
            selected = if("group" %in% available_vars) "group" else available_vars[1]
          )
          
          # Update shape variable choices  
          shiny::updateSelectInput(session, "shape_variable",
            choices = available_vars,
            selected = if("group" %in% available_vars) "group" else available_vars[1]
          )
        }
      }
    })
    
    # Auto-trigger pre-normalization QC when normalization tab is accessed
    # This observes when the filtering steps are complete
    observe({
      # Check if protein filtering is complete using CORRECT R6 pattern
      if (!is.null(workflow_data$state_manager) && 
          workflow_data$state_manager$current_state == "protein_replicate_filtered" &&
          !norm_data$pre_norm_qc_generated) {
        
        message("Auto-triggering pre-normalization QC (chunk 24)")
        
        tryCatch({
          # Generate pre-normalization QC plots (chunk 24 equivalent)
          generatePreNormalizationQc()
          norm_data$pre_norm_qc_generated <- TRUE
          
        }, error = function(e) {
          message(paste("Error generating pre-normalization QC:", e$message))
          shiny::showNotification(
            paste("Error generating pre-normalization QC:", e$message),
            type = "error",
            duration = 10
          )
        })
      }
    })
    
    # Regenerate plots when aesthetic variables change
    observeEvent(c(input$color_variable, input$shape_variable), {
      # Only regenerate if we have plots already
      if (norm_data$pre_norm_qc_generated) {
        message("Regenerating pre-normalization QC plots with new aesthetics...")
        tryCatch({
          generatePreNormalizationQc()
        }, error = function(e) {
          message(paste("Error regenerating pre-normalization QC:", e$message))
        })
      }
      
      # Also regenerate post-normalization plots if they exist
      if (norm_data$normalization_complete && !is.null(norm_data$normalized_protein_obj)) {
        message("Regenerating post-normalization QC plots with new aesthetics...")
        tryCatch({
          generatePostNormalizationQc(norm_data$normalized_protein_obj)
        }, error = function(e) {
          message(paste("Error regenerating post-normalization QC:", e$message))
        })
      }
      
      # Also regenerate RUV-corrected plots if they exist
      if (norm_data$ruv_complete && !is.null(norm_data$ruv_normalized_obj)) {
        message("Regenerating RUV-corrected QC plots with new aesthetics...")
        tryCatch({
          generateRuvCorrectedQc(norm_data$ruv_normalized_obj)
        }, error = function(e) {
          message(paste("Error regenerating RUV-corrected QC:", e$message))
        })
      }
    })
    
    # Normalization button logic
    observeEvent(input$run_normalization, {
      message("Starting normalization workflow...")
      
      tryCatch({
        # Get current S4 object using CORRECT R6 methods
        shiny::req(workflow_data$state_manager)
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        
        if (is.null(current_s4)) {
          stop("No S4 object available for normalization")
        }
        
        # Update progress
        shiny::withProgress(message = "Running normalization...", value = 0, {
          
          # Step 1: Between-samples normalization (chunk 25 equivalent)
          shiny::incProgress(0.2, detail = "Normalizing between samples...")
          
          # Update config with user selection
          if (exists("config_list", envir = .GlobalEnv)) {
            config_list$normaliseBetweenSamples$method <- input$norm_method
            assign("config_list", config_list, envir = .GlobalEnv)
          }
          
          normalized_s4 <- normaliseBetweenSamples(current_s4, normalisation_method = input$norm_method)
          norm_data$normalized_protein_obj <- normalized_s4
          
          # Generate post-normalization QC plots
          shiny::incProgress(0.2, detail = "Generating post-normalization QC plots...")
          generatePostNormalizationQc(normalized_s4)
          norm_data$normalization_complete <- TRUE
          
          # Step 2: RUV-III correction (chunks 26-27 equivalent)
          shiny::incProgress(0.2, detail = "Determining RUV parameters...")
          
          # Get RUV parameters based on mode
          if (input$ruv_mode == "automatic") {
            # Use default parameters (placeholder for automatic tuning)
            percentage_as_neg_ctrl <- 5
            ruv_k <- 3  # Default, will be optimized
          } else {
            percentage_as_neg_ctrl <- input$ruv_percentage
            ruv_k <- if(is.null(input$ruv_k) || is.na(input$ruv_k)) 3 else input$ruv_k
          }
          
          # Get negative control proteins
          control_genes_index <- getNegCtrlProtAnova(
            normalized_s4,
            ruv_grouping_variable = getPlotAesthetics()$color_var,
            percentage_as_neg_ctrl = percentage_as_neg_ctrl,
            ruv_qval_cutoff = 0.05,
            ruv_fdr_method = "BH"
          )
          norm_data$control_genes_index <- control_genes_index
          
          # Determine best k using canonical correlation
          if (input$ruv_mode == "automatic") {
            cancor_result <- ruvCancor(
              normalized_s4,
              ctrl = control_genes_index,
              num_components_to_impute = 2,
              ruv_grouping_variable = getPlotAesthetics()$color_var
            )
            # Find best k (simplified - take the one with highest canonical correlation)
            best_k <- which.max(cancor_result) + 1  # Add 1 since index starts at 0
            norm_data$best_k <- best_k
          } else {
            best_k <- ruv_k
            norm_data$best_k <- best_k
          }
          
          # Apply RUV-III correction
          shiny::incProgress(0.2, detail = "Applying RUV-III batch correction...")
          
          ruv_corrected_s4 <- ruvIII_C_Varying(
            normalized_s4,
            ruv_grouping_variable = getPlotAesthetics()$color_var,
            ruv_number_k = best_k,
            ctrl = control_genes_index
          )
          
          # Clean up any proteins with excessive missing values after RUV
          ruv_corrected_s4_clean <- removeRowsWithMissingValuesPercent(
            ruv_corrected_s4,
            ruv_grouping_variable = getPlotAesthetics()$color_var,
            groupwise_percentage_cutoff = 60,
            max_groups_percentage_cutoff = 60,
            proteins_intensity_cutoff_percentile = 1
          )
          
          norm_data$ruv_normalized_obj <- ruv_corrected_s4_clean
          
          # Generate RUV-corrected QC plots
          shiny::incProgress(0.2, detail = "Generating RUV-corrected QC plots...")
          generateRuvCorrectedQc(ruv_corrected_s4_clean)
          norm_data$ruv_complete <- TRUE
          
          # Update workflow data with final normalized object
          workflow_data$ruv_normalised_for_de_analysis_obj <- ruv_corrected_s4_clean
          
          # Update tab status to enable differential expression
          workflow_data$tab_status$normalization <- "complete"
          workflow_data$tab_status$differential_expression <- "pending"
          
        })
        
        shiny::showNotification(
          "Normalization and RUV correction completed successfully!",
          type = "success",
          duration = 5
        )
        
        message("Normalization workflow completed successfully")
        
      }, error = function(e) {
        message(paste("Error in normalization workflow:", e$message))
        shiny::showNotification(
          paste("Error in normalization:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
    
    # Reset normalization button logic
    observeEvent(input$reset_normalization, {
      message("Resetting normalization...")
      
      # Reset normalization state
      norm_data$normalization_complete <- FALSE
      norm_data$ruv_complete <- FALSE
      norm_data$normalized_protein_obj <- NULL
      norm_data$ruv_normalized_obj <- NULL
      norm_data$best_k <- NULL
      norm_data$control_genes_index <- NULL
      
      # Clear post-normalization and RUV plots
      norm_data$qc_plots$post_normalization <- list(pca = NULL, density = NULL, rle = NULL, correlation = NULL)
      norm_data$qc_plots$ruv_corrected <- list(pca = NULL, density = NULL, rle = NULL, correlation = NULL)
      
      # Reset workflow data
      workflow_data$ruv_normalised_for_de_analysis_obj <- NULL
      workflow_data$tab_status$normalization <- "pending"
      workflow_data$tab_status$differential_expression <- "disabled"
      
      shiny::showNotification(
        "Normalization has been reset to pre-normalization state",
        type = "warning",
        duration = 3
      )
      
      message("Normalization reset completed")
    })
    
    # Render QC plots - Post-filtering column
    output$pca_post_filtering <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_filtering$pca)) {
        norm_data$qc_plots$post_filtering$pca
      } else {
        plot.new()
        text(0.5, 0.5, "Pre-normalization QC not yet generated", cex = 1.2)
      }
    })
    
    output$density_post_filtering <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_filtering$density)) {
        norm_data$qc_plots$post_filtering$density
      } else {
        plot.new()
        text(0.5, 0.5, "Pre-normalization QC not yet generated", cex = 1.2)
      }
    })
    
    output$rle_post_filtering <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_filtering$rle)) {
        norm_data$qc_plots$post_filtering$rle
      } else {
        plot.new()
        text(0.5, 0.5, "Pre-normalization QC not yet generated", cex = 1.2)
      }
    })
    
    output$correlation_post_filtering <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_filtering$correlation)) {
        norm_data$qc_plots$post_filtering$correlation
      } else {
        plot.new()
        text(0.5, 0.5, "Pre-normalization QC not yet generated", cex = 1.2)
      }
    })
    
    # Render QC plots - Post-normalization column
    output$pca_post_normalization <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_normalization$pca)) {
        norm_data$qc_plots$post_normalization$pca
      } else {
        plot.new()
        text(0.5, 0.5, "Run normalization to generate plots", cex = 1.2)
      }
    })
    
    output$density_post_normalization <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_normalization$density)) {
        norm_data$qc_plots$post_normalization$density
      } else {
        plot.new()
        text(0.5, 0.5, "Run normalization to generate plots", cex = 1.2)
      }
    })
    
    output$rle_post_normalization <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_normalization$rle)) {
        norm_data$qc_plots$post_normalization$rle
      } else {
        plot.new()
        text(0.5, 0.5, "Run normalization to generate plots", cex = 1.2)
      }
    })
    
    output$correlation_post_normalization <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_normalization$correlation)) {
        norm_data$qc_plots$post_normalization$correlation
      } else {
        plot.new()
        text(0.5, 0.5, "Run normalization to generate plots", cex = 1.2)
      }
    })
    
    # Render QC plots - RUV-corrected column
    output$pca_ruv_corrected <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$ruv_corrected$pca)) {
        norm_data$qc_plots$ruv_corrected$pca
      } else {
        plot.new()
        text(0.5, 0.5, "Run RUV correction to generate plots", cex = 1.2)
      }
    })
    
    output$density_ruv_corrected <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$ruv_corrected$density)) {
        norm_data$qc_plots$ruv_corrected$density
      } else {
        plot.new()
        text(0.5, 0.5, "Run RUV correction to generate plots", cex = 1.2)
      }
    })
    
    output$rle_ruv_corrected <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$ruv_corrected$rle)) {
        norm_data$qc_plots$ruv_corrected$rle
      } else {
        plot.new()
        text(0.5, 0.5, "Run RUV correction to generate plots", cex = 1.2)
      }
    })
    
    output$correlation_ruv_corrected <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$ruv_corrected$correlation)) {
        norm_data$qc_plots$ruv_corrected$correlation
      } else {
        plot.new()
        text(0.5, 0.5, "Run RUV correction to generate plots", cex = 1.2)
      }
    })
    
    # Return normalization data for potential use by parent module
    return(norm_data)
  })
} 