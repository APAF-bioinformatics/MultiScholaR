# ============================================================================
# mod_metab_norm.R
# ============================================================================
# Purpose: Metabolomics normalization Shiny module
#
# This module provides a multi-step normalization pipeline with visualization
# of before/after effects for metabolomics data.
# ============================================================================

#' @title Metabolomics Normalization Module
#' @description A Shiny module for multi-step normalization of metabolomics data.
#'              Includes ITSD normalization, log2 transformation, cyclic loess,
#'              and optional RUV-III-C.
#' @name mod_metab_norm
NULL

#' @rdname mod_metab_norm
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags checkboxInput numericInput plotOutput conditionalPanel
#' @importFrom shinyjqui jqui_resizable
mod_metab_norm_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::h3("Metabolite Normalization Pipeline")
                    
                    # Progress indicator
                    , shiny::uiOutput(ns("normalization_progress"))
                    
                    , shiny::hr()
                    
                    , shiny::fluidRow(
                        # Left column: Normalization steps
                        shiny::column(4
                            # Step 1: ITSD Normalization
                            , shiny::wellPanel(
                                style = "border-left: 4px solid #3498db;"
                                , shiny::h4(
                                    shiny::icon("flask")
                                    , " Step 1: ITSD Normalization"
                                )
                                , shiny::p("Normalize using internal standards (optional).")
                                , shiny::checkboxInput(
                                    ns("apply_itsd")
                                    , "Apply ITSD normalization"
                                    , value = FALSE
                                )
                                , shiny::conditionalPanel(
                                    condition = paste0("input['", ns("apply_itsd"), "']")
                                    , shiny::selectInput(
                                        ns("itsd_method")
                                        , "Normalization method"
                                        , choices = c(
                                            "Median of IS" = "median"
                                            , "Mean of IS" = "mean"
                                            , "Single IS" = "single"
                                        )
                                    )
                                )
                                , shiny::actionButton(
                                    ns("run_itsd")
                                    , "Apply ITSD"
                                    , class = "btn-primary"
                                    , width = "100%"
                                )
                            )
                            
                            # Step 2: Log2 Transform
                            , shiny::wellPanel(
                                style = "border-left: 4px solid #2ecc71;"
                                , shiny::h4(
                                    shiny::icon("superscript")
                                    , " Step 2: Log2 Transform"
                                )
                                , shiny::p("Apply log2 transformation to stabilize variance.")
                                , shiny::numericInput(
                                    ns("log_offset")
                                    , "Offset for zeros"
                                    , value = 1
                                    , min = 0
                                    , step = 0.1
                                )
                                , shiny::actionButton(
                                    ns("run_log2")
                                    , "Apply Log2"
                                    , class = "btn-primary"
                                    , width = "100%"
                                )
                            )
                            
                            # Step 3: Cyclic Loess
                            , shiny::wellPanel(
                                style = "border-left: 4px solid #9b59b6;"
                                , shiny::h4(
                                    shiny::icon("wave-square")
                                    , " Step 3: Cyclic Loess"
                                )
                                , shiny::p("Between-sample normalization using cyclic loess.")
                                , shiny::numericInput(
                                    ns("loess_span")
                                    , "Span parameter"
                                    , value = 0.4
                                    , min = 0.1
                                    , max = 1.0
                                    , step = 0.1
                                )
                                , shiny::actionButton(
                                    ns("run_loess")
                                    , "Apply Cyclic Loess"
                                    , class = "btn-primary"
                                    , width = "100%"
                                )
                            )
                            
                            # Step 4: RUV-III-C (optional)
                            , shiny::wellPanel(
                                style = "border-left: 4px solid #e74c3c;"
                                , shiny::h4(
                                    shiny::icon("filter")
                                    , " Step 4: RUV-III-C (Optional)"
                                )
                                , shiny::p("Remove unwanted variation using negative controls.")
                                , shiny::checkboxInput(
                                    ns("apply_ruv")
                                    , "Apply RUV-III-C"
                                    , value = FALSE
                                )
                                , shiny::conditionalPanel(
                                    condition = paste0("input['", ns("apply_ruv"), "']")
                                    , shiny::numericInput(
                                        ns("ruv_k")
                                        , "Number of factors (k)"
                                        , value = 2
                                        , min = 1
                                        , max = 10
                                        , step = 1
                                    )
                                    , shiny::numericInput(
                                        ns("ruv_pval")
                                        , "P-value threshold for neg controls"
                                        , value = 0.5
                                        , min = 0.1
                                        , max = 1.0
                                        , step = 0.1
                                    )
                                )
                                , shiny::actionButton(
                                    ns("run_ruv")
                                    , "Apply RUV-III-C"
                                    , class = "btn-primary"
                                    , width = "100%"
                                )
                            )
                            
                            , shiny::hr()
                            
                            # Finalize
                            , shiny::actionButton(
                                ns("finalize_norm")
                                , "Finalize Normalization"
                                , class = "btn-success"
                                , width = "100%"
                                , icon = shiny::icon("check")
                            )
                        )
                        
                        # Right column: Visualization
                        , shiny::column(8
                            , shiny::h4("Normalization Results")
                            , shiny::verbatimTextOutput(ns("norm_log"))
                            
                            , shiny::hr()
                            
                            # Visualization tabs
                            , shiny::tabsetPanel(
                                id = ns("viz_tabs")
                                , shiny::tabPanel(
                                    "RLE Plot"
                                    , shiny::br()
                                    , shinyjqui::jqui_resizable(
                                        shiny::plotOutput(ns("rle_plot"), height = "400px")
                                    )
                                )
                                , shiny::tabPanel(
                                    "PCA"
                                    , shiny::br()
                                    , shinyjqui::jqui_resizable(
                                        shiny::plotOutput(ns("pca_plot"), height = "400px")
                                    )
                                )
                                , shiny::tabPanel(
                                    "Density"
                                    , shiny::br()
                                    , shinyjqui::jqui_resizable(
                                        shiny::plotOutput(ns("density_plot"), height = "400px")
                                    )
                                )
                                , shiny::tabPanel(
                                    "Boxplot"
                                    , shiny::br()
                                    , shinyjqui::jqui_resizable(
                                        shiny::plotOutput(ns("boxplot"), height = "400px")
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
}

#' @rdname mod_metab_norm
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText renderPlot showNotification removeNotification tags
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_density geom_point labs theme_minimal theme facet_wrap scale_color_brewer coord_flip
#' @importFrom logger log_info log_error log_warn
mod_metab_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        # Local reactive values
        local_data <- shiny::reactiveValues(
            normalization_log = character(0)
            , current_step = 0
            , steps_completed = c(itsd = FALSE, log2 = FALSE, loess = FALSE, ruv = FALSE)
        )
        
        # Normalization progress indicator
        output$normalization_progress <- shiny::renderUI({
            steps <- c(
                "ITSD" = local_data$steps_completed["itsd"]
                , "Log2" = local_data$steps_completed["log2"]
                , "Loess" = local_data$steps_completed["loess"]
                , "RUV" = local_data$steps_completed["ruv"]
            )
            
            step_badges <- lapply(names(steps), function(step_name) {
                completed <- steps[[step_name]]
                badge_class <- if (completed) "badge-success" else "badge-secondary"
                icon_name <- if (completed) "check" else "circle"
                
                shiny::tags$span(
                    class = paste("badge", badge_class)
                    , style = "margin-right: 10px; padding: 8px 12px;"
                    , shiny::icon(icon_name)
                    , step_name
                )
            })
            
            shiny::tags$div(
                style = "text-align: center; margin-bottom: 15px;"
                , step_badges
            )
        })
        
        # Helper function to add to log
        add_log <- function(message) {
            timestamp <- format(Sys.time(), "%H:%M:%S")
            local_data$normalization_log <- c(
                local_data$normalization_log
                , sprintf("[%s] %s", timestamp, message)
            )
        }
        
        # Render normalization log
        output$norm_log <- shiny::renderText({
            if (length(local_data$normalization_log) == 0) {
                return("Normalization log will appear here as you apply steps...")
            }
            paste(local_data$normalization_log, collapse = "\n")
        })
        
        # Step 1: ITSD Normalization
        shiny::observeEvent(input$run_itsd, {
            shiny::req(workflow_data$state_manager)
            
            if (!input$apply_itsd) {
                add_log("ITSD normalization skipped")
                local_data$steps_completed["itsd"] <- TRUE
                return()
            }
            
            shiny::showNotification("Applying ITSD normalization...", id = "norm_working", duration = NULL)
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                # Apply ITSD normalization
                normalized_s4 <- normaliseUntransformedData(
                    theObject = current_s4
                    , method = input$itsd_method
                )
                
                workflow_data$state_manager$saveState(
                    state_name = "metab_itsd_norm"
                    , s4_data_object = normalized_s4
                    , config_object = workflow_data$config_list
                    , description = sprintf("ITSD normalization applied (method: %s)", input$itsd_method)
                )
                
                add_log(sprintf("ITSD normalization applied (method: %s)", input$itsd_method))
                local_data$steps_completed["itsd"] <- TRUE
                
                shiny::removeNotification("norm_working")
                shiny::showNotification("ITSD normalization complete", type = "message")
                
            }, error = function(e) {
                add_log(sprintf("ERROR in ITSD: %s", e$message))
                logger::log_error(paste("ITSD normalization error:", e$message))
                shiny::removeNotification("norm_working")
                shiny::showNotification(paste("Error:", e$message), type = "error")
            })
        })
        
        # Step 2: Log2 Transform
        shiny::observeEvent(input$run_log2, {
            shiny::req(workflow_data$state_manager)
            
            shiny::showNotification("Applying log2 transformation...", id = "norm_working", duration = NULL)
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                # Apply log2 transformation
                transformed_s4 <- logTransformAssays(
                    theObject = current_s4
                    , offset = input$log_offset
                )
                
                workflow_data$state_manager$saveState(
                    state_name = "metab_log2"
                    , s4_data_object = transformed_s4
                    , config_object = workflow_data$config_list
                    , description = sprintf("Log2 transformation applied (offset: %g)", input$log_offset)
                )
                
                add_log(sprintf("Log2 transformation applied (offset: %g)", input$log_offset))
                local_data$steps_completed["log2"] <- TRUE
                
                shiny::removeNotification("norm_working")
                shiny::showNotification("Log2 transformation complete", type = "message")
                
            }, error = function(e) {
                add_log(sprintf("ERROR in Log2: %s", e$message))
                logger::log_error(paste("Log2 transformation error:", e$message))
                shiny::removeNotification("norm_working")
                shiny::showNotification(paste("Error:", e$message), type = "error")
            })
        })
        
        # Step 3: Cyclic Loess
        shiny::observeEvent(input$run_loess, {
            shiny::req(workflow_data$state_manager)
            
            shiny::showNotification("Applying cyclic loess normalization...", id = "norm_working", duration = NULL)
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                # Apply cyclic loess
                normalized_s4 <- normaliseBetweenSamples(
                    theObject = current_s4
                    , method = "cyclicloess"
                    , span = input$loess_span
                )
                
                workflow_data$state_manager$saveState(
                    state_name = "metab_cyclicloess"
                    , s4_data_object = normalized_s4
                    , config_object = workflow_data$config_list
                    , description = sprintf("Cyclic loess normalization applied (span: %g)", input$loess_span)
                )
                
                add_log(sprintf("Cyclic loess normalization applied (span: %g)", input$loess_span))
                local_data$steps_completed["loess"] <- TRUE
                
                shiny::removeNotification("norm_working")
                shiny::showNotification("Cyclic loess normalization complete", type = "message")
                
            }, error = function(e) {
                add_log(sprintf("ERROR in Cyclic Loess: %s", e$message))
                logger::log_error(paste("Cyclic loess error:", e$message))
                shiny::removeNotification("norm_working")
                shiny::showNotification(paste("Error:", e$message), type = "error")
            })
        })
        
        # Step 4: RUV-III-C
        shiny::observeEvent(input$run_ruv, {
            shiny::req(workflow_data$state_manager)
            
            if (!input$apply_ruv) {
                add_log("RUV-III-C skipped")
                local_data$steps_completed["ruv"] <- TRUE
                return()
            }
            
            shiny::showNotification("Applying RUV-III-C...", id = "norm_working", duration = NULL)
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                # Apply RUV-III-C
                normalized_s4 <- ruvIII_C_Varying(
                    theObject = current_s4
                    , k = input$ruv_k
                    , pval_threshold = input$ruv_pval
                )
                
                workflow_data$state_manager$saveState(
                    state_name = "metab_ruv"
                    , s4_data_object = normalized_s4
                    , config_object = workflow_data$config_list
                    , description = sprintf("RUV-III-C applied (k: %d, pval: %g)", input$ruv_k, input$ruv_pval)
                )
                
                add_log(sprintf("RUV-III-C applied (k: %d, pval: %g)", input$ruv_k, input$ruv_pval))
                local_data$steps_completed["ruv"] <- TRUE
                
                shiny::removeNotification("norm_working")
                shiny::showNotification("RUV-III-C complete", type = "message")
                
            }, error = function(e) {
                add_log(sprintf("ERROR in RUV-III-C: %s", e$message))
                logger::log_error(paste("RUV-III-C error:", e$message))
                shiny::removeNotification("norm_working")
                shiny::showNotification(paste("Error:", e$message), type = "error")
            })
        })
        
        # Finalize normalization
        shiny::observeEvent(input$finalize_norm, {
            shiny::req(workflow_data$state_manager)
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                workflow_data$state_manager$saveState(
                    state_name = "metab_norm_complete"
                    , s4_data_object = current_s4
                    , config_object = workflow_data$config_list
                    , description = "Normalization pipeline complete"
                )
                
                workflow_data$tab_status$normalization <- "complete"
                
                add_log("Normalization pipeline finalized")
                
                shiny::showNotification("Normalization complete! Proceed to Differential Analysis.", type = "message")
                
            }, error = function(e) {
                logger::log_error(paste("Error finalizing normalization:", e$message))
                shiny::showNotification(paste("Error:", e$message), type = "error")
            })
        })
        
        # Visualization: RLE Plot
        output$rle_plot <- shiny::renderPlot({
            shiny::req(workflow_data$state_manager)
            
            current_s4 <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) NULL)
            
            if (is.null(current_s4) || !inherits(current_s4, "MetaboliteAssayData")) {
                return(NULL)
            }
            
            # Get first assay for visualization
            assay_list <- current_s4@metabolite_data
            if (length(assay_list) == 0) return(NULL)
            
            assay_data <- assay_list[[1]]
            quant_info <- getMetaboliteQuantData(assay_data)
            
            if (ncol(quant_info$quant_data) == 0) return(NULL)
            
            # Calculate RLE
            log_data <- log2(quant_info$quant_data + 1)
            medians <- apply(log_data, 1, median, na.rm = TRUE)
            rle_data <- sweep(log_data, 1, medians)
            
            # Convert to long format
            rle_long <- tidyr::pivot_longer(
                as.data.frame(rle_data)
                , cols = everything()
                , names_to = "Sample"
                , values_to = "RLE"
            )
            
            ggplot2::ggplot(rle_long, ggplot2::aes(x = Sample, y = RLE)) +
                ggplot2::geom_boxplot(fill = "#3498db", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                ggplot2::labs(
                    title = "Relative Log Expression (RLE) Plot"
                    , x = "Sample"
                    , y = "RLE"
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
        })
        
        # Visualization: PCA Plot
        output$pca_plot <- shiny::renderPlot({
            shiny::req(workflow_data$state_manager)
            
            current_s4 <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) NULL)
            
            if (is.null(current_s4) || !inherits(current_s4, "MetaboliteAssayData")) {
                return(NULL)
            }
            
            assay_list <- current_s4@metabolite_data
            if (length(assay_list) == 0) return(NULL)
            
            assay_data <- assay_list[[1]]
            quant_info <- getMetaboliteQuantData(assay_data)
            
            if (ncol(quant_info$quant_data) < 3) return(NULL)
            
            # Prepare data for PCA
            log_data <- log2(quant_info$quant_data + 1)
            log_data <- log_data[complete.cases(log_data), ]
            
            if (nrow(log_data) < 10) return(NULL)
            
            # Run PCA
            pca_result <- prcomp(t(log_data), scale. = TRUE)
            pca_df <- as.data.frame(pca_result$x[, 1:2])
            pca_df$Sample <- rownames(pca_df)
            
            # Get group info if available
            dm <- current_s4@design_matrix
            group_col <- current_s4@group_id
            sample_col <- current_s4@sample_id
            
            if (nrow(dm) > 0 && group_col %in% names(dm) && sample_col %in% names(dm)) {
                pca_df <- merge(pca_df, dm[, c(sample_col, group_col)], by.x = "Sample", by.y = sample_col, all.x = TRUE)
                names(pca_df)[names(pca_df) == group_col] <- "Group"
            } else {
                pca_df$Group <- "Sample"
            }
            
            var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
            
            ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
                ggplot2::geom_point(size = 3) +
                ggplot2::labs(
                    title = "PCA Plot"
                    , x = sprintf("PC1 (%.1f%%)", var_explained[1])
                    , y = sprintf("PC2 (%.1f%%)", var_explained[2])
                ) +
                ggplot2::theme_minimal() +
                ggplot2::scale_color_brewer(palette = "Set1")
        })
        
        # Visualization: Density Plot
        output$density_plot <- shiny::renderPlot({
            shiny::req(workflow_data$state_manager)
            
            current_s4 <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) NULL)
            
            if (is.null(current_s4) || !inherits(current_s4, "MetaboliteAssayData")) {
                return(NULL)
            }
            
            assay_list <- current_s4@metabolite_data
            if (length(assay_list) == 0) return(NULL)
            
            assay_data <- assay_list[[1]]
            quant_info <- getMetaboliteQuantData(assay_data)
            
            if (ncol(quant_info$quant_data) == 0) return(NULL)
            
            # Convert to long format
            log_data <- log2(quant_info$quant_data + 1)
            long_data <- tidyr::pivot_longer(
                as.data.frame(log_data)
                , cols = everything()
                , names_to = "Sample"
                , values_to = "Intensity"
            )
            
            ggplot2::ggplot(long_data, ggplot2::aes(x = Intensity, color = Sample)) +
                ggplot2::geom_density(alpha = 0.5) +
                ggplot2::labs(
                    title = "Intensity Distribution (log2)"
                    , x = "Log2 Intensity"
                    , y = "Density"
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(legend.position = "none")
        })
        
        # Visualization: Boxplot
        output$boxplot <- shiny::renderPlot({
            shiny::req(workflow_data$state_manager)
            
            current_s4 <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) NULL)
            
            if (is.null(current_s4) || !inherits(current_s4, "MetaboliteAssayData")) {
                return(NULL)
            }
            
            assay_list <- current_s4@metabolite_data
            if (length(assay_list) == 0) return(NULL)
            
            assay_data <- assay_list[[1]]
            quant_info <- getMetaboliteQuantData(assay_data)
            
            if (ncol(quant_info$quant_data) == 0) return(NULL)
            
            log_data <- log2(quant_info$quant_data + 1)
            long_data <- tidyr::pivot_longer(
                as.data.frame(log_data)
                , cols = everything()
                , names_to = "Sample"
                , values_to = "Intensity"
            )
            
            ggplot2::ggplot(long_data, ggplot2::aes(x = Sample, y = Intensity)) +
                ggplot2::geom_boxplot(fill = "#2ecc71", alpha = 0.7) +
                ggplot2::labs(
                    title = "Sample Intensity Distribution (log2)"
                    , x = "Sample"
                    , y = "Log2 Intensity"
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
        })
    })
}

