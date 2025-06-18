#' @title differentialExpressionAppletModule
#'
#' @description A Shiny module for the Differential Expression step of the proteomics
#' workflow. Handles limma-based statistical analysis and visualization of results
#' across user-defined contrasts.
#'
#' @name differentialExpressionAppletModule
NULL

#' @rdname differentialExpressionAppletModule
#' @export
#' @import shiny
#' @import shinydashboard
differentialExpressionAppletUI <- function(id) {
  ns <- NS(id)
  
  shiny::fluidRow(
    shiny::column(3,
      shiny::wellPanel(
        shiny::h4("DE Analysis Settings"),
        
        # Formula input (will be populated from S4 @args)
        shiny::textAreaInput(
          ns("formula_string"),
          "Model Formula:",
          value = "~ 0 + group",
          height = "60px"
        ),
        shiny::helpText("Formula fetched from S4 object @args"),
        
        shiny::hr(),
        
        # Analysis parameters
        shiny::h5("Analysis Parameters"),
        shiny::numericInput(
          ns("de_q_val_thresh"),
          "Q-value Threshold:",
          value = 0.05,
          min = 0.001,
          max = 0.2,
          step = 0.005
        ),
        shiny::numericInput(
          ns("treat_lfc_cutoff"),
          "Log Fold-Change Cutoff:",
          value = 0,
          min = 0,
          max = 2,
          step = 0.1
        ),
        
        shiny::hr(),
        
        # Available contrasts display
        shiny::h5("Available Contrasts"),
        shiny::verbatimTextOutput(ns("contrasts_display")),
        shiny::helpText("Contrasts defined in design matrix"),
        
        shiny::br(),
        
        # Main action button
        shiny::actionButton(
          ns("run_de_analysis"),
          "Run DE Analysis",
          class = "btn-primary",
          width = "100%",
          icon = shiny::icon("play")
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Download button
        shiny::downloadButton(
          ns("download_de_results"),
          "Download All Results",
          class = "btn-success",
          width = "100%"
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Status display
        shiny::h5("Analysis Status"),
        shiny::verbatimTextOutput(ns("de_status"))
      )
    ),
    
    # Main results panel with 3 tabs
    shiny::column(9,
      shiny::tabsetPanel(
        id = ns("de_results_tabs"),
        
        # Tab 1: Volcano Plot (Interactive with Glimma)
        shiny::tabPanel(
          "Volcano Plot",
          icon = shiny::icon("chart-scatter"),
          shiny::br(),
          
          # Contrast selector for volcano plot
          shiny::fluidRow(
            shiny::column(6,
              shiny::selectInput(
                ns("volcano_contrast"),
                "Select Contrast:",
                choices = NULL,
                width = "100%"
              )
            ),
            shiny::column(6,
              shiny::checkboxInput(
                ns("volcano_interactive"),
                "Interactive Plot (Glimma)",
                value = TRUE
              )
            )
          ),
          
          shiny::hr(),
          
          # Volcano plot output (will switch between static and interactive)
          shiny::conditionalPanel(
            condition = "input.volcano_interactive == true",
            ns = ns,
            shiny::div(
              id = ns("volcano_glimma_container"),
              style = "height: 600px;",
              shiny::htmlOutput(ns("volcano_glimma"))
            )
          ),
          
          shiny::conditionalPanel(
            condition = "input.volcano_interactive == false", 
            ns = ns,
            shinyjqui::jqui_resizable(
              shiny::plotOutput(ns("volcano_static"), height = "600px")
            )
          )
        ),
        
        # Tab 2: Heatmap (Interactive)
        shiny::tabPanel(
          "Heatmap",
          icon = shiny::icon("th"),
          shiny::br(),
          
          # Contrast selector and heatmap options
          shiny::fluidRow(
            shiny::column(3,
              shiny::selectInput(
                ns("heatmap_contrast"),
                "Select Contrast:",
                choices = NULL,
                width = "100%"
              )
            ),
            shiny::column(3,
              shiny::numericInput(
                ns("heatmap_top_n"),
                "Top N Genes:",
                value = 50,
                min = 10,
                max = 500,
                step = 10
              )
            ),
            shiny::column(3,
              shiny::selectInput(
                ns("heatmap_clustering"),
                "Apply Clustering:",
                choices = list(
                  "Both rows & columns" = "both",
                  "Rows only" = "row", 
                  "Columns only" = "column",
                  "None" = "none"
                ),
                selected = "both"
              )
            ),
            shiny::column(3,
              shiny::selectInput(
                ns("heatmap_scaling"),
                "Data Scaling:",
                choices = list(
                  "Row (gene) scaling" = "row",
                  "Column (sample) scaling" = "column",
                  "Both" = "both",
                  "None" = "none"
                ),
                selected = "row"
              )
            )
          ),
          
          # Advanced clustering controls (collapsible)
          shiny::conditionalPanel(
            condition = "input.heatmap_clustering != 'none'",
            ns = ns,
            shiny::wellPanel(
              shiny::h5("Advanced Clustering Options"),
              
              shiny::fluidRow(
                shiny::column(3,
                  shiny::selectInput(
                    ns("heatmap_cluster_method"),
                    "Clustering Method:",
                    choices = list(
                      "Ward (minimum variance)" = "ward.D2",
                      "Ward (original)" = "ward.D",
                      "Complete linkage" = "complete",
                      "Single linkage" = "single", 
                      "Average linkage" = "average",
                      "McQuitty (WPGMA)" = "mcquitty",
                      "Median (WPGMC)" = "median",
                      "Centroid (UPGMC)" = "centroid"
                    ),
                    selected = "ward.D2"
                  )
                ),
                shiny::column(3,
                  shiny::selectInput(
                    ns("heatmap_distance_method"),
                    "Distance Metric:",
                    choices = list(
                      "Euclidean" = "euclidean",
                      "Manhattan" = "manhattan",
                      "Pearson correlation" = "pearson",
                      "Spearman correlation" = "spearman",
                      "Maximum" = "maximum",
                      "Canberra" = "canberra",
                      "Binary" = "binary",
                      "Minkowski" = "minkowski"
                    ),
                    selected = "euclidean"
                  )
                ),
                shiny::column(3,
                  shiny::selectInput(
                    ns("heatmap_tree_cut_method"),
                    "Tree Cutting:",
                    choices = list(
                      "Number of clusters" = "k_clusters",
                      "Height cutoff" = "height_cutoff",
                      "Dynamic tree cutting" = "dynamic",
                      "No cutting" = "none"
                    ),
                    selected = "k_clusters"
                  )
                ),
                shiny::column(3,
                  # Dynamic UI for tree cutting parameters
                  shiny::conditionalPanel(
                    condition = "input.heatmap_tree_cut_method == 'k_clusters'",
                    ns = ns,
                    shiny::numericInput(
                      ns("heatmap_n_clusters"),
                      "Number of Clusters:",
                      value = 4,
                      min = 2,
                      max = 20,
                      step = 1
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = "input.heatmap_tree_cut_method == 'height_cutoff'",
                    ns = ns,
                    shiny::numericInput(
                      ns("heatmap_cut_height"),
                      "Cut Height:",
                      value = 0.5,
                      min = 0.1,
                      max = 2.0,
                      step = 0.1
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = "input.heatmap_tree_cut_method == 'dynamic'",
                    ns = ns,
                    shiny::numericInput(
                      ns("heatmap_min_cluster_size"),
                      "Min Cluster Size:",
                      value = 3,
                      min = 2,
                      max = 20,
                      step = 1
                    )
                  )
                )
              ),
              
              # Additional heatmap display options
              shiny::fluidRow(
                shiny::column(4,
                  shiny::checkboxInput(
                    ns("heatmap_show_dendro"),
                    "Show Dendrogram",
                    value = TRUE
                  )
                ),
                shiny::column(4,
                  shiny::checkboxInput(
                    ns("heatmap_show_labels"),
                    "Show Gene Labels",
                    value = FALSE
                  )
                ),
                shiny::column(4,
                  shiny::selectInput(
                    ns("heatmap_color_scheme"),
                    "Color Scheme:",
                    choices = list(
                      "Red-Blue" = "RdBu",
                      "Red-Yellow-Blue" = "RdYlBu", 
                      "Blue-White-Red" = "coolwarm",
                      "Viridis" = "viridis",
                      "Plasma" = "plasma",
                      "Inferno" = "inferno"
                    ),
                    selected = "RdBu"
                  )
                )
              )
            )
          ),
          
          shiny::hr(),
          
          # Interactive heatmap output
          shiny::div(
            id = ns("heatmap_container"),
            style = "height: 650px;",
            shinyjqui::jqui_resizable(
              shiny::plotOutput(ns("heatmap_plot"), height = "600px")
            )
          ),
          
          # Cluster summary (if tree cutting is applied)
          shiny::conditionalPanel(
            condition = "input.heatmap_tree_cut_method != 'none'",
            ns = ns,
            shiny::wellPanel(
              shiny::h5("Cluster Information"),
              shiny::verbatimTextOutput(ns("cluster_summary"))
            )
          )
        ),
        
        # Tab 3: DE Results Table (Interactive)
        shiny::tabPanel(
          "DE Results Table",
          icon = shiny::icon("table"),
          shiny::br(),
          
          # Contrast selector and table options
          shiny::fluidRow(
            shiny::column(4,
              shiny::selectInput(
                ns("table_contrast"),
                "Select Contrast:",
                choices = NULL,
                width = "100%"
              )
            ),
            shiny::column(4,
              shiny::selectInput(
                ns("table_significance"),
                "Show:",
                choices = list(
                  "All Results" = "all",
                  "Significant Only" = "significant",
                  "Up-regulated" = "up",
                  "Down-regulated" = "down"
                ),
                selected = "significant"
              )
            ),
            shiny::column(4,
              shiny::numericInput(
                ns("table_max_rows"),
                "Max Rows:",
                value = 1000,
                min = 100,
                max = 10000,
                step = 100
              )
            )
          ),
          
          shiny::hr(),
          
          # Summary statistics
          shiny::fluidRow(
            shiny::column(12,
              shiny::wellPanel(
                shiny::h5("Summary Statistics"),
                shiny::verbatimTextOutput(ns("de_summary_stats"))
              )
            )
          ),
          
          # Interactive results table
          DT::DTOutput(ns("de_results_table"))
        )
      )
    )
  )
}

#' @rdname differentialExpressionAppletModule 
#' @export
differentialExpressionAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    message("=== DIFFERENTIAL EXPRESSION MODULE SERVER STARTED ===")
    message(sprintf("Module ID: %s", id))
    
    # Initialize reactive values for DE state
    de_data <- reactiveValues(
      de_results_list = NULL,
      contrasts_available = NULL,
      analysis_complete = FALSE,
      current_s4_object = NULL,
      formula_from_s4 = NULL,
      current_row_clusters = NULL,
      current_col_clusters = NULL
    )
    
    # Fetch formula and contrasts from S4 object when module loads
    observe({
      # Get current S4 object from state manager
      if (!is.null(workflow_data$state_manager)) {
        current_state <- workflow_data$state_manager$current_state
        
        # Should be getting the normalized protein object
        if (current_state %in% c("protein_replicate_filtered", "normalized", "ruv_corrected")) {
          current_s4 <- workflow_data$state_manager$getState(current_state)
          
          if (!is.null(current_s4)) {
            de_data$current_s4_object <- current_s4
            
            # Extract formula from S4 @args
            if ("deAnalysisParameters" %in% names(current_s4@args)) {
              if ("formula_string" %in% names(current_s4@args$deAnalysisParameters)) {
                formula_from_s4 <- current_s4@args$deAnalysisParameters$formula_string
                de_data$formula_from_s4 <- formula_from_s4
                
                # Update UI with formula from S4
                shiny::updateTextAreaInput(
                  session,
                  "formula_string",
                  value = formula_from_s4
                )
              }
            }
            
            # Get contrasts from design matrix or S4 args
            if (!is.null(current_s4@design_matrix)) {
              # For now, extract unique groups to create basic contrasts
              # This should ideally come from the design matrix applet
              if (exists("contrasts_tbl", envir = .GlobalEnv)) {
                contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
                de_data$contrasts_available <- contrasts_tbl$comparison
              } else {
                # Create basic contrasts from groups
                groups <- unique(current_s4@design_matrix$group)
                if (length(groups) >= 2) {
                  # Create pairwise contrasts (placeholder)
                  basic_contrasts <- combn(groups, 2, function(x) paste(x[1], "-", x[2]), simplify = TRUE)
                  de_data$contrasts_available <- basic_contrasts
                }
              }
            }
          }
        }
      }
    })
    
    # Update contrast dropdowns when contrasts become available
    observe({
      if (!is.null(de_data$contrasts_available)) {
        contrast_choices <- setNames(de_data$contrasts_available, de_data$contrasts_available)
        
        # Update all contrast selectors
        shiny::updateSelectInput(session, "volcano_contrast", choices = contrast_choices)
        shiny::updateSelectInput(session, "heatmap_contrast", choices = contrast_choices)
        shiny::updateSelectInput(session, "table_contrast", choices = contrast_choices)
      }
    })
    
    # Display available contrasts
    output$contrasts_display <- shiny::renderText({
      if (!is.null(de_data$contrasts_available)) {
        paste(de_data$contrasts_available, collapse = "\n")
      } else {
        "No contrasts available.\nComplete normalization first."
      }
    })
    
    # Display analysis status
    output$de_status <- shiny::renderText({
      if (de_data$analysis_complete) {
        paste(
          "✅ Analysis Complete\n",
          sprintf("Contrasts analyzed: %d\n", length(de_data$contrasts_available)),
          "Results available in all tabs"
        )
      } else {
        "⏳ Waiting for analysis...\nClick 'Run DE Analysis' to start"
      }
    })
    
    # Main DE Analysis Button Logic
    observeEvent(input$run_de_analysis, {
      message("=== STARTING DIFFERENTIAL EXPRESSION ANALYSIS ===")
      
      shiny::req(de_data$current_s4_object, de_data$contrasts_available)
      
      shiny::showNotification("Running differential expression analysis...", id = "de_working", duration = NULL)
      
      tryCatch({
        shiny::withProgress(message = "Running DE analysis...", value = 0, {
          
          # Step 1: Prepare data for limma analysis
          shiny::incProgress(0.2, detail = "Preparing data for limma...")
          
          # This is where the actual limma-based DE analysis would happen
          # For now, create placeholder structure
          
          # Step 2: Run limma analysis for each contrast
          shiny::incProgress(0.4, detail = "Running limma analysis...")
          
          # Placeholder: Create mock DE results for each contrast
          mock_de_results <- list()
          for (contrast in de_data$contrasts_available) {
            # This would be replaced with actual deAnalysisWrapperFunction call
            mock_de_results[[contrast]] <- data.frame(
              Protein.Ids = paste0("PROTEIN_", 1:100),
              logFC = rnorm(100, 0, 2),
              AveExpr = rnorm(100, 10, 3),
              t = rnorm(100, 0, 3),
              P.Value = runif(100, 0, 1),
              adj.P.Val = runif(100, 0, 1),
              B = rnorm(100, 0, 2),
              stringsAsFactors = FALSE
            )
          }
          
          # Step 3: Generate plots for all contrasts
          shiny::incProgress(0.6, detail = "Generating plots...")
          
          # Placeholder for plot generation
          # This would create volcano plots, heatmaps, etc. for each contrast
          
          # Step 4: Save results
          shiny::incProgress(0.8, detail = "Saving results...")
          
          de_data$de_results_list <- mock_de_results
          de_data$analysis_complete <- TRUE
          
          # Update workflow data
          workflow_data$de_analysis_results_list <- mock_de_results
          workflow_data$tab_status$differential_expression <- "complete"
          workflow_data$tab_status$enrichment_analysis <- "pending"
          
          shiny::incProgress(1.0, detail = "Complete!")
        })
        
        shiny::showNotification(
          "Differential expression analysis completed successfully!",
          type = "success",
          duration = 5
        )
        
        message("=== DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED ===")
        
      }, error = function(e) {
        message(paste("*** ERROR in DE analysis:", e$message))
        shiny::showNotification(
          paste("Error in DE analysis:", e$message),
          type = "error",
          duration = 10
        )
      })
      
      shiny::removeNotification("de_working")
    })
    
    # Volcano Plot Rendering
    output$volcano_glimma <- shiny::renderUI({
      shiny::req(input$volcano_contrast, de_data$de_results_list)
      
      # Placeholder for interactive Glimma volcano plot
      shiny::div(
        style = "text-align: center; padding: 50px;",
        shiny::h4("Interactive Volcano Plot"),
        shiny::p(paste("Contrast:", input$volcano_contrast)),
        shiny::p("Glimma interactive plot would render here"),
        shiny::p("(Placeholder - will be replaced with actual glimmaVolcano output)")
      )
    })
    
    output$volcano_static <- shiny::renderPlot({
      shiny::req(input$volcano_contrast, de_data$de_results_list)
      
      # Placeholder for static volcano plot
      plot(1:10, 1:10, main = paste("Static Volcano Plot -", input$volcano_contrast),
           xlab = "Log2 Fold Change", ylab = "-Log10 P-value")
    })
    
    # Heatmap Rendering
    output$heatmap_plot <- shiny::renderPlot({
      shiny::req(input$heatmap_contrast, de_data$de_results_list)
      
      # Get the current results for the selected contrast
      current_results <- de_data$de_results_list[[input$heatmap_contrast]]
      
      if (!is.null(current_results)) {
        # Filter for significant results and get top N
        significant_results <- current_results[current_results$adj.P.Val < input$de_q_val_thresh, ]
        
        if (nrow(significant_results) > 0) {
          # Order by absolute fold change and take top N
          significant_results <- significant_results[order(abs(significant_results$logFC), decreasing = TRUE), ]
          top_genes <- head(significant_results, input$heatmap_top_n)
          
          # Create mock expression matrix for these genes
          # In real implementation, this would extract from the S4 object
          n_samples <- 12  # Mock number of samples
          expr_matrix <- matrix(
            rnorm(nrow(top_genes) * n_samples, mean = 0, sd = 2),
            nrow = nrow(top_genes),
            ncol = n_samples,
            dimnames = list(
              top_genes$Protein.Ids,
              paste0("Sample_", 1:n_samples)
            )
          )
          
          # Apply scaling
          if (input$heatmap_scaling == "row") {
            expr_matrix <- t(scale(t(expr_matrix)))
          } else if (input$heatmap_scaling == "column") {
            expr_matrix <- scale(expr_matrix)
          } else if (input$heatmap_scaling == "both") {
            expr_matrix <- scale(t(scale(t(expr_matrix))))
          }
          
          # Set up clustering parameters
          if (input$heatmap_clustering != "none") {
            # Calculate distance matrix based on user selection
            if (input$heatmap_distance_method %in% c("pearson", "spearman")) {
              # For correlation-based distances
              if (input$heatmap_distance_method == "pearson") {
                row_dist <- as.dist(1 - cor(t(expr_matrix), method = "pearson"))
                col_dist <- as.dist(1 - cor(expr_matrix, method = "pearson"))
              } else {
                row_dist <- as.dist(1 - cor(t(expr_matrix), method = "spearman"))
                col_dist <- as.dist(1 - cor(expr_matrix, method = "spearman"))
              }
            } else {
              # For standard distance metrics
              row_dist <- dist(expr_matrix, method = input$heatmap_distance_method)
              col_dist <- dist(t(expr_matrix), method = input$heatmap_distance_method)
            }
            
            # Perform hierarchical clustering
            row_clust <- hclust(row_dist, method = input$heatmap_cluster_method)
            col_clust <- hclust(col_dist, method = input$heatmap_cluster_method)
            
            # Apply tree cutting if requested
            if (input$heatmap_tree_cut_method != "none") {
              if (input$heatmap_tree_cut_method == "k_clusters") {
                row_clusters <- cutree(row_clust, k = input$heatmap_n_clusters)
                col_clusters <- cutree(col_clust, k = input$heatmap_n_clusters)
              } else if (input$heatmap_tree_cut_method == "height_cutoff") {
                row_clusters <- cutree(row_clust, h = input$heatmap_cut_height)
                col_clusters <- cutree(col_clust, h = input$heatmap_cut_height)
              } else if (input$heatmap_tree_cut_method == "dynamic") {
                # For dynamic tree cutting, use a simplified approach
                # In real implementation, would use dynamicTreeCut package
                row_clusters <- cutree(row_clust, k = max(2, min(input$heatmap_min_cluster_size, nrow(expr_matrix)/3)))
                col_clusters <- cutree(col_clust, k = max(2, min(input$heatmap_min_cluster_size, ncol(expr_matrix)/3)))
              }
              
              # Store cluster information for summary
              de_data$current_row_clusters <- row_clusters
              de_data$current_col_clusters <- col_clusters
            }
          }
          
          # Create the heatmap
          # Set up clustering parameters for heatmap function
          if (input$heatmap_clustering == "both") {
            cluster_rows <- if (input$heatmap_clustering != "none") row_clust else FALSE
            cluster_cols <- if (input$heatmap_clustering != "none") col_clust else FALSE
          } else if (input$heatmap_clustering == "row") {
            cluster_rows <- if (input$heatmap_clustering != "none") row_clust else FALSE
            cluster_cols <- FALSE
          } else if (input$heatmap_clustering == "column") {
            cluster_rows <- FALSE
            cluster_cols <- if (input$heatmap_clustering != "none") col_clust else FALSE
          } else {
            cluster_rows <- FALSE
            cluster_cols <- FALSE
          }
          
          # Choose color palette
          colors <- switch(input$heatmap_color_scheme,
            "RdBu" = colorRampPalette(c("red", "white", "blue"))(100),
            "RdYlBu" = colorRampPalette(c("red", "yellow", "blue"))(100),
            "coolwarm" = colorRampPalette(c("blue", "white", "red"))(100),
            "viridis" = viridis::viridis(100),
            "plasma" = viridis::plasma(100),
            "inferno" = viridis::inferno(100),
            colorRampPalette(c("red", "white", "blue"))(100)
          )
          
          # Generate the heatmap
          if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
            # Use ComplexHeatmap if available for more control
            ComplexHeatmap::Heatmap(
              expr_matrix,
              name = "Expression",
              col = colors,
              cluster_rows = cluster_rows,
              cluster_columns = cluster_cols,
              show_row_names = input$heatmap_show_labels,
              show_column_names = TRUE,
              show_row_dend = input$heatmap_show_dendro && input$heatmap_clustering %in% c("both", "row"),
              show_column_dend = input$heatmap_show_dendro && input$heatmap_clustering %in% c("both", "column"),
              column_title = paste("Heatmap - Top", input$heatmap_top_n, "genes\nContrast:", input$heatmap_contrast),
              heatmap_legend_param = list(title = "Scaled\nExpression")
            )
          } else {
            # Fallback to base heatmap
            heatmap(
              expr_matrix,
              main = paste("Heatmap - Top", input$heatmap_top_n, "genes\nContrast:", input$heatmap_contrast),
              col = colors,
              Rowv = if (input$heatmap_clustering %in% c("both", "row")) as.dendrogram(row_clust) else NA,
              Colv = if (input$heatmap_clustering %in% c("both", "column")) as.dendrogram(col_clust) else NA,
              labRow = if (input$heatmap_show_labels) rownames(expr_matrix) else rep("", nrow(expr_matrix)),
              cexRow = if (input$heatmap_show_labels) 0.8 else 0.1
            )
          }
        } else {
          # No significant genes found
          plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
               main = paste("No significant genes found for contrast:", input$heatmap_contrast))
          text(1, 1, "Adjust significance thresholds\nor select different contrast", cex = 1.2)
        }
      }
    })
    
    # Cluster Summary Output
    output$cluster_summary <- shiny::renderText({
      if (!is.null(de_data$current_row_clusters) && !is.null(de_data$current_col_clusters)) {
        row_cluster_summary <- table(de_data$current_row_clusters)
        col_cluster_summary <- table(de_data$current_col_clusters)
        
        paste0(
          "Gene Clusters:\n",
          paste(sprintf("Cluster %s: %d genes", names(row_cluster_summary), row_cluster_summary), collapse = "\n"),
          "\n\nSample Clusters:\n",
          paste(sprintf("Cluster %s: %d samples", names(col_cluster_summary), col_cluster_summary), collapse = "\n"),
          "\n\nClustering Method: ", input$heatmap_cluster_method,
          "\nDistance Metric: ", input$heatmap_distance_method,
          "\nTree Cutting: ", input$heatmap_tree_cut_method
        )
      } else {
        "No cluster information available.\nEnable tree cutting to see cluster details."
      }
    })
    
    # DE Results Table Rendering
    output$de_results_table <- DT::renderDT({
      shiny::req(input$table_contrast, de_data$de_results_list)
      
      current_results <- de_data$de_results_list[[input$table_contrast]]
      
      if (!is.null(current_results)) {
        # Filter based on significance selection
        if (input$table_significance == "significant") {
          current_results <- current_results[current_results$adj.P.Val < input$de_q_val_thresh, ]
        } else if (input$table_significance == "up") {
          current_results <- current_results[current_results$adj.P.Val < input$de_q_val_thresh & 
                                           current_results$logFC > input$treat_lfc_cutoff, ]
        } else if (input$table_significance == "down") {
          current_results <- current_results[current_results$adj.P.Val < input$de_q_val_thresh & 
                                           current_results$logFC < -input$treat_lfc_cutoff, ]
        }
        
        # Limit rows
        if (nrow(current_results) > input$table_max_rows) {
          current_results <- current_results[1:input$table_max_rows, ]
        }
        
        DT::datatable(
          current_results,
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel')
          ),
          extensions = 'Buttons'
        ) |>
        DT::formatRound(columns = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"), digits = 4)
      }
    })
    
    # Summary statistics
    output$de_summary_stats <- shiny::renderText({
      shiny::req(input$table_contrast, de_data$de_results_list)
      
      current_results <- de_data$de_results_list[[input$table_contrast]]
      
      if (!is.null(current_results)) {
        total_genes <- nrow(current_results)
        significant <- sum(current_results$adj.P.Val < input$de_q_val_thresh, na.rm = TRUE)
        up_reg <- sum(current_results$adj.P.Val < input$de_q_val_thresh & 
                     current_results$logFC > input$treat_lfc_cutoff, na.rm = TRUE)
        down_reg <- sum(current_results$adj.P.Val < input$de_q_val_thresh & 
                       current_results$logFC < -input$treat_lfc_cutoff, na.rm = TRUE)
        
        paste(
          sprintf("Total genes: %d", total_genes),
          sprintf("Significant (q < %.3f): %d", input$de_q_val_thresh, significant),
          sprintf("Up-regulated: %d", up_reg),
          sprintf("Down-regulated: %d", down_reg),
          sprintf("Fold-change cutoff: %.2f", input$treat_lfc_cutoff),
          sep = "\n"
        )
      }
    })
    
    # Download handler for results
    output$download_de_results <- shiny::downloadHandler(
      filename = function() {
        paste0("DE_results_", Sys.Date(), ".zip")
      },
      content = function(file) {
        # Placeholder for creating downloadable zip file
        writeLines("Differential expression results would be packaged here", file)
      }
    )
    
    # Return DE data for potential use by parent module
    return(de_data)
  })
} 