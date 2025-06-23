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
differentialExpressionAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    cat("--- Entering differentialExpressionAppletServer ---\n")
    cat(sprintf("   differentialExpressionAppletServer Arg: id = %s\n", id))
    cat(sprintf("   differentialExpressionAppletServer Arg: workflow_data is NULL = %s\n", is.null(workflow_data)))
    cat(sprintf("   differentialExpressionAppletServer Arg: selected_tab is NULL = %s\n", is.null(selected_tab)))
    
    cat("=== DIFFERENTIAL EXPRESSION MODULE SERVER STARTED ===\n")
    cat(sprintf("Module ID: %s\n", id))
    
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
    
    cat("   differentialExpressionAppletServer Step: Reactive values initialized\n")
    
    # CRITICAL FIX: Add tab selection observer like normalization module
    if (!is.null(selected_tab)) {
      cat("   differentialExpressionAppletServer Step: Setting up tab selection observer\n")
      observeEvent(selected_tab(), {
        
        cat("--- Entering tab selection observer ---\n")
        cat(sprintf("   tab_observer Step: Selected tab = %s\n", selected_tab()))
        
        # Only trigger if DE tab is selected
        if (!is.null(selected_tab()) && selected_tab() == "de") {
          
          cat("=== DE TAB CLICKED ===\n")
          cat(sprintf("   DE TAB Step: workflow_data$state_manager is NULL = %s\n", is.null(workflow_data$state_manager)))
          
          if (!is.null(workflow_data$state_manager)) {
            current_state <- workflow_data$state_manager$current_state
            
            # Define valid states where this tab can be active
            valid_states_for_de_tab <- c("correlation_filtered")
            
            cat(sprintf("   DE TAB Step: Current state = '%s'\n", current_state))
            cat(sprintf("   DE TAB Step: Valid states for DE = %s\n", paste(valid_states_for_de_tab, collapse = ", ")))
            
            # Auto-trigger only fires if we've completed correlation filtering
            if (current_state == "correlation_filtered") {
              
              cat("*** AUTO-TRIGGERING DE INITIALIZATION (correlation-filtered state found) ***\n")
              
              tryCatch({
                # Initialize DE analysis setup
                cat("   DE TAB Step: Getting S4 object from state manager...\n")
                current_s4 <- workflow_data$state_manager$getState(current_state)
                
                if (!is.null(current_s4)) {
                  cat(sprintf("   DE TAB Step: S4 object retrieved, class = %s\n", class(current_s4)))
                  de_data$current_s4_object <- current_s4
                  
                  # Check for contrasts_tbl in global environment FIRST
                  cat("   DE TAB Step: Checking for contrasts_tbl in global environment...\n")
                  cat(sprintf("   DE TAB Step: contrasts_tbl exists in global env: %s\n", exists("contrasts_tbl", envir = .GlobalEnv)))
                  if (exists("contrasts_tbl", envir = .GlobalEnv)) {
                    contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
                    cat("   DE TAB Step: Found contrasts_tbl in global environment\n")
                    cat("   DE TAB Step: contrasts_tbl structure:\n")
                    str(contrasts_tbl)
                    cat("   DE TAB Step: contrasts_tbl content:\n")
                    print(contrasts_tbl)
                    
                    # Validate contrasts_tbl has content
                    if (is.null(contrasts_tbl) || nrow(contrasts_tbl) == 0) {
                      cat("   DE TAB Step: contrasts_tbl exists but is empty, falling back to auto-generation\n")
                      de_data$contrasts_available <- NULL
                    } else if ("comparison" %in% names(contrasts_tbl)) {
                      de_data$contrasts_available <- contrasts_tbl$comparison
                      cat(sprintf("   DE TAB Step: Set contrasts from comparison column: %s\n", paste(de_data$contrasts_available, collapse = ", ")))
                    } else if ("contrasts" %in% names(contrasts_tbl)) {
                      de_data$contrasts_available <- contrasts_tbl$contrasts
                      cat(sprintf("   DE TAB Step: Set contrasts from contrasts column: %s\n", paste(de_data$contrasts_available, collapse = ", ")))
                    } else {
                      cat("   DE TAB Step: contrasts_tbl found but no recognized column names\n")
                      cat("   DE TAB Step: Available column names:\n")
                      print(names(contrasts_tbl))
                      # Try first column if it has content
                      if (ncol(contrasts_tbl) > 0) {
                        de_data$contrasts_available <- contrasts_tbl[[1]]
                        cat(sprintf("   DE TAB Step: Using first column: %s\n", paste(de_data$contrasts_available, collapse = ", ")))
                      } else {
                        cat("   DE TAB Step: contrasts_tbl has no usable content\n")
                        de_data$contrasts_available <- NULL
                      }
                    }
                  } else {
                    cat("   DE TAB Step: No contrasts_tbl in global environment. Will attempt auto-generation.\n")
                    de_data$contrasts_available <- NULL
                  }
                  
                  # Only auto-generate contrasts if user-specified ones weren't found or were empty
                  if (is.null(de_data$contrasts_available) || length(de_data$contrasts_available) == 0) {
                    cat("   DE TAB Step: No valid user-specified contrasts found, creating basic contrasts\n")
                    shiny::showNotification("No user-specified contrasts found. Auto-generating all possible contrasts.", type = "warning", duration = 5)
                    
                    if (!is.null(current_s4@design_matrix)) {
                      groups <- unique(current_s4@design_matrix$group)
                      if (length(groups) >= 2) {
                        group_prefixed <- paste0("group", groups)
                        basic_contrasts <- combn(group_prefixed, 2, function(x) paste0(x[1], "-", x[2]), simplify = TRUE)
                        de_data$contrasts_available <- basic_contrasts
                        cat(sprintf("   DE TAB Step: Created basic contrasts: %s\n", paste(basic_contrasts, collapse = ", ")))
                      }
                    }
                  }
                  
                  # Extract formula from S4 @args
                  cat("   DE TAB Step: Checking for formula in S4 @args...\n")
                  if ("deAnalysisParameters" %in% names(current_s4@args)) {
                    if ("formula_string" %in% names(current_s4@args$deAnalysisParameters)) {
                      formula_from_s4 <- current_s4@args$deAnalysisParameters$formula_string
                      de_data$formula_from_s4 <- formula_from_s4
                      cat(sprintf("   DE TAB Step: Formula from S4 = %s\n", formula_from_s4))
                      
                      # Update UI with formula from S4
                      shiny::updateTextAreaInput(
                        session,
                        "formula_string",
                        value = formula_from_s4
                      )
                    } else {
                      cat("   DE TAB Step: No formula_string in deAnalysisParameters\n")
                    }
                  } else {
                    cat("   DE TAB Step: No deAnalysisParameters in S4 @args\n")
                  }
                  
                  # CRITICAL FIX: Validate and fix contrast format to match design matrix columns
                  if (!is.null(de_data$contrasts_available) && !is.null(current_s4@design_matrix)) {
                    cat("   DE TAB Step: Validating contrast format against design matrix...\n")
                    groups <- unique(current_s4@design_matrix$group)
                    cat(sprintf("   DE TAB Step: Available groups: %s\n", paste(groups, collapse = ", ")))
                    
                    # Check if contrasts need "group" prefix (for formula ~ 0 + group)
                    current_contrasts <- de_data$contrasts_available
                    cat(sprintf("   DE TAB Step: Current contrasts before validation: %s\n", paste(current_contrasts, collapse = ", ")))
                    
                    # If contrasts don't start with "group" but should (based on formula), fix them
                    if (any(grepl("~ 0 \\+ group", input$formula_string))) {
                      cat("   DE TAB Step: Formula uses '~ 0 + group', checking if contrasts need group prefix...\n")
                      
                      # Check if any contrast references raw group names without prefix
                      needs_prefix <- any(sapply(groups, function(g) any(grepl(paste0("\\b", g, "\\b"), current_contrasts))))
                      
                      if (needs_prefix) {
                        cat("   DE TAB Step: Contrasts need group prefix, fixing...\n")
                        fixed_contrasts <- current_contrasts
                        for (group in groups) {
                          # Replace "GA_Control" with "groupGA_Control" etc.
                          fixed_contrasts <- gsub(paste0("\\b", group, "\\b"), paste0("group", group), fixed_contrasts)
                        }
                        de_data$contrasts_available <- fixed_contrasts
                        cat(sprintf("   DE TAB Step: Fixed contrasts: %s\n", paste(fixed_contrasts, collapse = ", ")))
                      } else {
                        cat("   DE TAB Step: Contrasts already have correct format\n")
                      }
                    }
                  }
                  
                } else {
                  cat("   DE TAB Step: S4 object is NULL\n")
                }
                
                cat("*** DE INITIALIZATION COMPLETED SUCCESSFULLY ***\n")
                
              }, error = function(e) {
                cat(paste("*** ERROR in DE initialization:", e$message, "\n"))
                shiny::showNotification(
                  paste("Error initializing DE analysis:", e$message),
                  type = "error",
                  duration = 10
                )
              })
              
            } else {
              cat(sprintf("*** State '%s' is not valid for DE analysis. User needs to complete normalization and correlation filtering. ***\n", current_state))
              shiny::showNotification(
                "Please complete the normalization and correlation filtering steps before accessing differential expression analysis.",
                type = "warning",
                duration = 5
              )
            }
          } else {
            cat("*** workflow_data$state_manager is NULL - cannot check state ***\n")
          }
        } else {
          cat(sprintf("   tab_observer Step: Tab '%s' is not DE tab, ignoring\n", selected_tab()))
        }
        
        cat("--- Exiting tab selection observer ---\n")
      }, ignoreInit = TRUE)
    } else {
      cat("   differentialExpressionAppletServer Step: No selected_tab parameter provided - tab selection observer NOT set up\n")
    }
    
    # BACKUP: Fetch formula and contrasts from S4 object when state updates (fallback method)
    observeEvent(workflow_data$state_update_trigger, {
      cat("--- Entering state update trigger observer ---\n")
      cat("\n\n=== DE TAB TRIGGERED VIA STATE UPDATE: Checking for S4 object and contrasts ===\n")
      
      # Get current S4 object from state manager
      if (!is.null(workflow_data$state_manager)) {
        current_state <- workflow_data$state_manager$current_state
        cat(sprintf("   DE TAB Step: Current state = %s\n", current_state))
        
        # Should be getting the correlation-filtered protein object (final state before DE)
        if (current_state == "correlation_filtered") {
          cat("   DE TAB Step: State is valid for DE analysis (correlation-filtered state found)\n")
          current_s4 <- workflow_data$state_manager$getState(current_state)
          
          if (!is.null(current_s4)) {
            cat(sprintf("   DE TAB Step: S4 object retrieved, class = %s\n", class(current_s4)))
            de_data$current_s4_object <- current_s4
            
            # Extract formula from S4 @args
            if ("deAnalysisParameters" %in% names(current_s4@args)) {
              if ("formula_string" %in% names(current_s4@args$deAnalysisParameters)) {
                formula_from_s4 <- current_s4@args$deAnalysisParameters$formula_string
                de_data$formula_from_s4 <- formula_from_s4
                cat(sprintf("   DE TAB Step: Formula from S4 = %s\n", formula_from_s4))
                
                # Update UI with formula from S4
                shiny::updateTextAreaInput(
                  session,
                  "formula_string",
                  value = formula_from_s4
                )
              } else {
                cat("   DE TAB Step: No formula_string in deAnalysisParameters\n")
              }
            } else {
              cat("   DE TAB Step: No deAnalysisParameters in S4 @args\n")
            }
            
            # Get contrasts from design matrix or S4 args
            if (!is.null(current_s4@design_matrix)) {
              cat("   DE TAB Step: Design matrix found in S4 object\n")
              cat(sprintf("   DE TAB Step: Design matrix dims = %d rows, %d cols\n", nrow(current_s4@design_matrix), ncol(current_s4@design_matrix)))
              
              # Check for contrasts_tbl in global environment
              cat(sprintf("   DE TAB Step: contrasts_tbl exists in global env: %s\n", exists("contrasts_tbl", envir = .GlobalEnv)))
              if (exists("contrasts_tbl", envir = .GlobalEnv)) {
                contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
                cat("   DE TAB Step: Found contrasts_tbl in global environment\n")
                cat("   DE TAB Step: contrasts_tbl structure:\n")
                str(contrasts_tbl)
                cat("   DE TAB Step: contrasts_tbl content:\n")
                print(contrasts_tbl)
                
                # Validate contrasts_tbl has content
                if (is.null(contrasts_tbl) || nrow(contrasts_tbl) == 0) {
                  cat("   DE TAB Step: contrasts_tbl exists but is empty, falling back to auto-generation\n")
                  de_data$contrasts_available <- NULL
                } else if ("comparison" %in% names(contrasts_tbl)) {
                  de_data$contrasts_available <- contrasts_tbl$comparison
                  cat(sprintf("   DE TAB Step: Set contrasts from comparison column: %s\n", paste(de_data$contrasts_available, collapse = ", ")))
                } else if ("contrasts" %in% names(contrasts_tbl)) {
                  de_data$contrasts_available <- contrasts_tbl$contrasts
                  cat(sprintf("   DE TAB Step: Set contrasts from contrasts column: %s\n", paste(de_data$contrasts_available, collapse = ", ")))
                } else {
                  cat("   DE TAB Step: contrasts_tbl found but no recognized column names\n")
                  cat("   DE TAB Step: Available column names:\n")
                  print(names(contrasts_tbl))
                  # Try first column if it has content
                  if (ncol(contrasts_tbl) > 0) {
                    de_data$contrasts_available <- contrasts_tbl[[1]]
                    cat(sprintf("   DE TAB Step: Using first column: %s\n", paste(de_data$contrasts_available, collapse = ", ")))
                  } else {
                    cat("   DE TAB Step: contrasts_tbl has no usable content\n")
                    de_data$contrasts_available <- NULL
                  }
                }
              } else {
                cat("   DE TAB Step: No contrasts_tbl in global environment. Will attempt auto-generation.\n")
                de_data$contrasts_available <- NULL
              }
              
              # Only auto-generate contrasts if user-specified ones weren't found or were empty
              if (is.null(de_data$contrasts_available) || length(de_data$contrasts_available) == 0) {
                cat("   DE TAB Step: No valid user-specified contrasts found, creating basic contrasts\n")
                shiny::showNotification("No user-specified contrasts found. Auto-generating all possible contrasts.", type = "warning", duration = 5)
                
                if (!is.null(current_s4@design_matrix)) {
                  groups <- unique(current_s4@design_matrix$group)
                  if (length(groups) >= 2) {
                    group_prefixed <- paste0("group", groups)
                    basic_contrasts <- combn(group_prefixed, 2, function(x) paste0(x[1], "-", x[2]), simplify = TRUE)
                    de_data$contrasts_available <- basic_contrasts
                    cat(sprintf("   DE TAB Step: Created basic contrasts: %s\n", paste(basic_contrasts, collapse = ", ")))
                  }
                }
              }
            } else {
              cat("   DE TAB Step: No design matrix found in S4 object\n")
            }
          } else {
            cat("   DE TAB Step: S4 object is NULL\n")
          }
        } else {
          cat(sprintf("   DE TAB Step: State '%s' not valid for DE analysis (expecting: correlation_filtered)\n", current_state))
        }
      } else {
        cat("   DE TAB Step: workflow_data$state_manager is NULL\n")
      }
      
      cat("=== DE TAB: Contrast detection complete ===\n")
      cat("--- Exiting state update trigger observer ---\n")
    }, ignoreInit = TRUE)
    
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
        "No contrasts available.\nComplete normalization and\ncorrelation filtering first."
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
      cat("=== STARTING DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
      
      shiny::req(de_data$current_s4_object, de_data$contrasts_available)
      
      shiny::showNotification("Running differential expression analysis...", id = "de_working", duration = NULL)
      
      tryCatch({
        shiny::withProgress(message = "Running DE analysis...", value = 0, {
          
          # Step 1: Prepare contrasts table for analysis
          shiny::incProgress(0.2, detail = "Preparing contrasts for analysis...")
          
          # Get contrasts from global environment or create basic ones
          contrasts_tbl <- NULL
          if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
            cat("   DE ANALYSIS Step: Using existing contrasts_tbl from global environment\n")
            cat("   DE ANALYSIS Step: contrasts_tbl structure:\n")
            str(contrasts_tbl)
          } else {
            # Create contrasts table that matches the original format
            # The original expects a data frame with contrasts in the first column
            cat("   DE ANALYSIS Step: Creating contrasts_tbl from de_data$contrasts_available\n")
            cat("   DE ANALYSIS Step: Available contrasts:\n")
            print(de_data$contrasts_available)
            
            contrasts_tbl <- data.frame(
              contrasts = de_data$contrasts_available,
              stringsAsFactors = FALSE
            )
            cat("   DE ANALYSIS Step: Created contrasts_tbl:\n")
            str(contrasts_tbl)
            print(contrasts_tbl)
          }
          
          # DEBUG: Check what will be passed to runTestsContrasts - FIX: Get ALL contrasts, not just first
          contrast_strings_to_use <- contrasts_tbl[, 1]  # Get all rows from first column, not just first element
          cat("   DE ANALYSIS Step: contrast_strings_to_use (what goes to runTestsContrasts):\n")
          print(contrast_strings_to_use)
          cat("   DE ANALYSIS Step: Length of contrast_strings_to_use:\n")
          print(length(contrast_strings_to_use))
          
          # Step 2: Run actual differential expression analysis
          shiny::incProgress(0.4, detail = "Running limma analysis...")
          
          # CRITICAL FIX: Follow the original RMarkdown approach - run DE analysis for ONE contrast at a time
          cat("   DE ANALYSIS Step: Running DE analysis for each contrast separately (like original workflow)\n")
          
          # Use original contrast names as keys to match UI dropdown expectations
          # Don't transform them - keep them exactly as they appear in contrasts_tbl
          contrast_names <- contrasts_tbl$contrasts
          
          cat(sprintf("   DE ANALYSIS Step: Generated contrast names: %s\n", paste(contrast_names, collapse = ", ")))
          
          # Run DE analysis for each contrast separately using purrr::map (like original)
          de_results_list <- seq_len(nrow(contrasts_tbl)) |>
            purrr::set_names(contrast_names) |>
            purrr::map(\(contrast_idx) {
              cat(sprintf("   DE ANALYSIS Step: Processing contrast %d: %s\n", contrast_idx, contrasts_tbl$contrasts[contrast_idx]))
              
              # Get single contrast (one row) like original workflow
              single_contrast_tbl <- contrasts_tbl |> dplyr::slice(contrast_idx)
              cat(sprintf("   DE ANALYSIS Step: Single contrast table for index %d:\n", contrast_idx))
              cat("   DE ANALYSIS Step: Structure of single_contrast_tbl:\n")
              str(single_contrast_tbl)
              cat("   DE ANALYSIS Step: Content of single_contrast_tbl:\n")
              print(single_contrast_tbl)
              cat(sprintf("   DE ANALYSIS Step: Contrast string being passed: '%s'\n", single_contrast_tbl$contrasts[1]))
              
              # Use the modular DE analysis function for this single contrast
              cat("   DE ANALYSIS Step: About to call differentialExpressionAnalysis...\n")
              cat(sprintf("   DE ANALYSIS Step: S4 object class: %s\n", class(de_data$current_s4_object)))
              cat(sprintf("   DE ANALYSIS Step: Formula: %s\n", input$formula_string))
              
              # CRITICAL: Ensure S4 object has the required parameters in @args
              # The differentialExpressionAnalysis function expects parameters in @args$differentialExpressionAnalysis
              if (is.null(de_data$current_s4_object@args$differentialExpressionAnalysis)) {
                cat("   DE ANALYSIS Step: Adding differentialExpressionAnalysis parameters to S4 @args\n")
                de_data$current_s4_object@args$differentialExpressionAnalysis <- list(
                  contrasts_tbl = single_contrast_tbl,
                  formula_string = input$formula_string,
                  de_q_val_thresh = input$de_q_val_thresh,
                  treat_lfc_cutoff = input$treat_lfc_cutoff,
                  eBayes_trend = TRUE,
                  eBayes_robust = TRUE,
                  args_group_pattern = "(\\d+)",
                  args_row_id = de_data$current_s4_object@protein_id_column,
                  qvalue_column = "fdr_qvalue",
                  raw_pvalue_column = "raw_pvalue"
                )
              }
              
              # Also add for the helper function
              if (is.null(de_data$current_s4_object@args$differentialExpressionAnalysisHelper)) {
                cat("   DE ANALYSIS Step: Adding differentialExpressionAnalysisHelper parameters to S4 @args\n")
                de_data$current_s4_object@args$differentialExpressionAnalysisHelper <- list(
                  contrasts_tbl = single_contrast_tbl,
                  formula_string = input$formula_string,
                  de_q_val_thresh = input$de_q_val_thresh,
                  treat_lfc_cutoff = input$treat_lfc_cutoff,
                  eBayes_trend = TRUE,
                  eBayes_robust = TRUE,
                  args_group_pattern = "(\\d+)",
                  args_row_id = de_data$current_s4_object@protein_id_column,
                  qvalue_column = "fdr_qvalue",
                  raw_pvalue_column = "raw_pvalue"
                )
              }
              
              result <- tryCatch({
                # DEBUG: Check if function exists
                cat("   DE ANALYSIS Step: Checking if differentialExpressionAnalysis function exists...\n")
                if (!exists("differentialExpressionAnalysis")) {
                  cat("   DE ANALYSIS Step: ERROR - differentialExpressionAnalysis function not found!\n")
                  cat("   DE ANALYSIS Step: Attempting to use MultiScholaR namespace...\n")
                  
                  # Try with explicit namespace
                  MultiScholaR::differentialExpressionAnalysis(
                    theObject = de_data$current_s4_object,
                    contrasts_tbl = single_contrast_tbl,
                    formula_string = input$formula_string,
                    de_q_val_thresh = input$de_q_val_thresh,
                    treat_lfc_cutoff = input$treat_lfc_cutoff,
                    qvalue_column = "fdr_qvalue",
                    raw_pvalue_column = "raw_pvalue"
                  )
                } else {
                  cat("   DE ANALYSIS Step: differentialExpressionAnalysis function found\n")
                  
                  # Verify S4 method dispatch
                  cat(sprintf("   DE ANALYSIS Step: S4 object class: %s\n", class(de_data$current_s4_object)))
                  cat(sprintf("   DE ANALYSIS Step: Available methods for differentialExpressionAnalysis:\n"))
                  print(showMethods("differentialExpressionAnalysis", printTo = FALSE))
                  
                  differentialExpressionAnalysis(
                    theObject = de_data$current_s4_object,
                    contrasts_tbl = single_contrast_tbl,
                    formula_string = input$formula_string,
                    de_q_val_thresh = input$de_q_val_thresh,
                    treat_lfc_cutoff = input$treat_lfc_cutoff,
                    qvalue_column = "fdr_qvalue",
                    raw_pvalue_column = "raw_pvalue"
                  )
                }
              }, error = function(e) {
                cat(paste("*** ERROR in differentialExpressionAnalysis call:", e$message, "\n"))
                cat(paste("*** ERROR object: \n"))
                print(str(e))
                cat("*** Attempting detailed error diagnosis:\n")
                
                # Check parameter types
                cat(sprintf("   - de_data$current_s4_object class: %s\n", class(de_data$current_s4_object)))
                cat(sprintf("   - single_contrast_tbl class: %s\n", class(single_contrast_tbl)))
                cat(sprintf("   - input$formula_string: %s\n", input$formula_string))
                cat(sprintf("   - input$de_q_val_thresh: %s\n", input$de_q_val_thresh))
                cat(sprintf("   - input$treat_lfc_cutoff: %s\n", input$treat_lfc_cutoff))
                
                stop(e)
              })
              
              cat(sprintf("   DE ANALYSIS Step: Successfully completed contrast %d\n", contrast_idx))
              result
            })
          
          cat(sprintf("   DE ANALYSIS Step: Completed DE analysis for %d contrasts\n", length(de_results_list)))
          cat(sprintf("   DE ANALYSIS Step: DE results list names: %s\n", paste(names(de_results_list), collapse = ", ")))
          
          # Step 3: Combine results into expected format for UI components
          shiny::incProgress(0.8, detail = "Processing results...")
          
          # Update UI dropdowns to match what's in the data's comparison column (before = sign)
          cat("   DE ANALYSIS Step: Updating UI dropdowns to match data comparison column...\n")
          # Extract just the comparison part (before =) to match what's in de_proteins_long$comparison
          comparison_names <- stringr::str_extract(contrast_names, "^[^=]+")
          contrast_choices <- setNames(comparison_names, comparison_names)
          shiny::updateSelectInput(session, "volcano_contrast", choices = contrast_choices)
          shiny::updateSelectInput(session, "heatmap_contrast", choices = contrast_choices)
          shiny::updateSelectInput(session, "table_contrast", choices = contrast_choices)
          cat(sprintf("   DE ANALYSIS Step: Updated UI dropdowns with comparison names: %s\n", paste(comparison_names, collapse = ", ")))
          cat(sprintf("   DE ANALYSIS Step: (These match the comparison column in de_proteins_long)\n"))
          
          # Combine all de_proteins_long results into a single dataframe
          cat("   DE ANALYSIS Step: Combining de_proteins_long results from all contrasts...\n")
          combined_de_proteins_long <- NULL
          
          for(contrast_name in names(de_results_list)) {
            if (!is.null(de_results_list[[contrast_name]]$de_proteins_long)) {
              contrast_results <- de_results_list[[contrast_name]]$de_proteins_long
              cat(sprintf("   DE ANALYSIS Step: Adding %d rows from contrast %s\n", nrow(contrast_results), contrast_name))
              
              if (is.null(combined_de_proteins_long)) {
                combined_de_proteins_long <- contrast_results
              } else {
                combined_de_proteins_long <- dplyr::bind_rows(combined_de_proteins_long, contrast_results)
              }
            }
          }
          
          cat(sprintf("   DE ANALYSIS Step: Combined de_proteins_long has %d total rows\n", nrow(combined_de_proteins_long)))
          
          # Create the expected structure for UI components
          combined_results <- list(
            de_proteins_long = combined_de_proteins_long,
            individual_contrasts = de_results_list  # Keep individual results for other purposes
          )
          
          # Add other shared elements from the first contrast result (they should be the same)
          if (length(de_results_list) > 0) {
            first_result <- de_results_list[[1]]
            if (!is.null(first_result$theObject)) {
              combined_results$theObject <- first_result$theObject
            }
          }
          
          de_data$de_results_list <- combined_results
          de_data$analysis_complete <- TRUE
          
          # Update workflow data
          workflow_data$de_analysis_results_list <- de_results_list
          workflow_data$tab_status$differential_expression <- "complete"
          workflow_data$tab_status$enrichment_analysis <- "pending"
          
          shiny::incProgress(1.0, detail = "Complete!")
        })
        
        shiny::showNotification(
          "Differential expression analysis completed successfully!",
          type = "success",
          duration = 5
        )
        
        cat("=== DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED ===\n")
        
      }, error = function(e) {
        cat(paste("*** ERROR in DE analysis:", e$message, "\n"))
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
      
      tryCatch({
        # Get UniProt annotations if available
        uniprot_tbl <- NULL
        if (exists("uniprot_dat_cln", envir = .GlobalEnv)) {
          uniprot_tbl <- get("uniprot_dat_cln", envir = .GlobalEnv)
        }
        
        # Generate interactive Glimma volcano plot using new function
        # Find the correct contrast results using functional programming
        contrast_keys <- names(de_data$de_results_list$individual_contrasts)
        
        # Extract comparison parts and find match using vectorized operations
        comparison_parts <- stringr::str_extract(contrast_keys, "^[^=]+")
        matching_index <- which(comparison_parts == input$volcano_contrast)
        
        if (length(matching_index) > 0) {
          full_contrast_name <- contrast_keys[matching_index[1]]
          selected_contrast_results <- de_data$de_results_list$individual_contrasts[[full_contrast_name]]
          cat(sprintf("   VOLCANO: Found matching contrast results for %s\n", full_contrast_name))
        } else {
          selected_contrast_results <- NULL
          full_contrast_name <- NULL
        }
        
        if(is.null(selected_contrast_results)) {
          cat(sprintf("   VOLCANO: No contrast results found for %s\n", input$volcano_contrast))
        } else {
          # Create a compatible data structure for the plotting function
          plot_data_structure <- list(
            de_proteins_long = de_data$de_results_list$de_proteins_long,
            contrasts_results = selected_contrast_results$contrasts_results,
            theObject = de_data$de_results_list$theObject
          )
          
          # The volcano contrast input should now match the comparison column in de_proteins_long
          cat(sprintf("   VOLCANO: Looking for contrast = %s\n", input$volcano_contrast))
          cat(sprintf("   VOLCANO: Using full contrast name for coefficients = %s\n", full_contrast_name))
          
          # Fix: Use the correct protein ID column from the S4 object
          protein_id_column <- de_data$de_results_list$theObject@protein_id_column
          cat(sprintf("   VOLCANO: Using protein ID column = %s\n", protein_id_column))
          
          glimma_widget <- generateVolcanoPlotGlimma(
            de_results_list = plot_data_structure,
            selected_contrast = full_contrast_name,  # Use full name for coefficient lookup
            uniprot_tbl = uniprot_tbl,
            de_q_val_thresh = input$de_q_val_thresh,
            args_row_id = protein_id_column  # Pass the correct column name
          )
        }
        
        if (!is.null(glimma_widget)) {
          # Return the Glimma widget
          glimma_widget
        } else {
          shiny::div(
            style = "text-align: center; padding: 20px;",
            shiny::h5("No data available for selected contrast"),
            shiny::p("Please check contrast selection and ensure analysis has completed.")
          )
        }
        
      }, error = function(e) {
        cat(paste("*** ERROR in volcano plot generation:", e$message, "\n"))
        shiny::div(
          style = "text-align: center; padding: 20px;",
          shiny::h5("Error generating volcano plot"),
          shiny::p(paste("Error:", e$message))
        )
      })
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
      
      tryCatch({
        # Generate heatmap using new function
        # Create a compatible data structure for the plotting function
        plot_data_structure <- list(
          de_proteins_long = de_data$de_results_list$de_proteins_long,
          theObject = de_data$de_results_list$theObject
        )
        
        # The heatmap contrast input should now match the comparison column in de_proteins_long
        cat(sprintf("   HEATMAP: Looking for contrast = %s\n", input$heatmap_contrast))
        
        heatmap_plot <- generateDEHeatmap(
          de_results_list = plot_data_structure,
          selected_contrast = input$heatmap_contrast,
          top_n_genes = input$heatmap_top_n,
          clustering_method = input$heatmap_cluster_method,
          distance_method = input$heatmap_distance_method,
          cluster_rows = input$heatmap_clustering %in% c("both", "row"),
          cluster_cols = input$heatmap_clustering %in% c("both", "column"),
          scale_data = input$heatmap_scaling,
          color_scheme = input$heatmap_color_scheme,
          show_gene_names = input$heatmap_show_labels,
          de_q_val_thresh = input$de_q_val_thresh
        )
        
        if (!is.null(heatmap_plot)) {
          if (is.function(heatmap_plot)) {
            # For base R heatmap functions
            heatmap_plot()
          } else {
            # For ComplexHeatmap objects
            heatmap_plot
          }
        } else {
          # No significant genes found
          plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
               main = paste("No significant genes found for contrast:", input$heatmap_contrast))
          text(1, 1, "Adjust significance thresholds\nor select different contrast", cex = 1.2)
        }
        
      }, error = function(e) {
        cat(paste("*** ERROR in heatmap generation:", e$message, "\n"))
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
             main = "Error generating heatmap")
        text(1, 1, paste("Error:", e$message), cex = 1.2)
      })
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
      
      tryCatch({
        # Get DE results from the new format
        if (!is.null(de_data$de_results_list$de_proteins_long)) {
          # Debug: Show available contrasts in data vs selected contrast
          available_contrasts <- unique(de_data$de_results_list$de_proteins_long$comparison)
          cat(sprintf("   DE TABLE: Available contrasts in data: %s\n", paste(available_contrasts, collapse = ", ")))
          cat(sprintf("   DE TABLE: Selected contrast from UI: %s\n", input$table_contrast))
          
          # Filter for selected contrast (input$table_contrast now matches comparison column)
          current_results <- de_data$de_results_list$de_proteins_long |>
            dplyr::filter(comparison == input$table_contrast)
          
          cat(sprintf("   DE TABLE: Found %d rows for contrast %s\n", nrow(current_results), input$table_contrast))
          
          if (nrow(current_results) > 0) {
            # Filter based on significance selection
            if (input$table_significance == "significant") {
              current_results <- current_results[current_results$fdr_qvalue < input$de_q_val_thresh, ]
            } else if (input$table_significance == "up") {
              current_results <- current_results[current_results$fdr_qvalue < input$de_q_val_thresh & 
                                               current_results$log2FC > input$treat_lfc_cutoff, ]
            } else if (input$table_significance == "down") {
              current_results <- current_results[current_results$fdr_qvalue < input$de_q_val_thresh & 
                                               current_results$log2FC < -input$treat_lfc_cutoff, ]
            }
            
            # Limit rows
            if (nrow(current_results) > input$table_max_rows) {
              current_results <- current_results[1:input$table_max_rows, ]
            }
            
            # Select relevant columns for display - use the correct protein ID column
            if (!is.null(de_data$de_results_list$theObject)) {
              protein_id_column <- de_data$de_results_list$theObject@protein_id_column
              cat(sprintf("   DE TABLE: Using protein ID column = %s\n", protein_id_column))
            } else {
              # Fallback: try to find protein ID column in the data
              possible_protein_cols <- c("Protein.Ids", "uniprot_acc", "protein_id", "Protein_ID")
              protein_id_column <- intersect(possible_protein_cols, names(current_results))[1]
              cat(sprintf("   DE TABLE: Using fallback protein ID column = %s\n", protein_id_column))
            }
            
            display_columns <- c(protein_id_column, "log2FC", "raw_pvalue", "fdr_qvalue")
            current_results <- current_results |>
              dplyr::select(any_of(display_columns))
            
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
            DT::formatRound(columns = c("log2FC", "raw_pvalue", "fdr_qvalue"), digits = 4)
          } else {
            # No results for this contrast
            DT::datatable(data.frame(Message = "No results available for selected contrast"))
          }
        } else {
          # No DE results available
          DT::datatable(data.frame(Message = "No DE analysis results available"))
        }
        
      }, error = function(e) {
        cat(paste("*** ERROR in DE results table:", e$message, "\n"))
        DT::datatable(data.frame(Message = paste("Error:", e$message)))
      })
    })
    
    # Summary statistics
    output$de_summary_stats <- shiny::renderText({
      shiny::req(input$table_contrast, de_data$de_results_list)
      
      tryCatch({
        if (!is.null(de_data$de_results_list$de_proteins_long)) {
          # Filter for selected contrast
          current_results <- de_data$de_results_list$de_proteins_long |>
            dplyr::filter(comparison == input$table_contrast)
          
          if (nrow(current_results) > 0) {
            total_genes <- nrow(current_results)
            significant <- sum(current_results$fdr_qvalue < input$de_q_val_thresh, na.rm = TRUE)
            up_reg <- sum(current_results$fdr_qvalue < input$de_q_val_thresh & 
                         current_results$log2FC > input$treat_lfc_cutoff, na.rm = TRUE)
            down_reg <- sum(current_results$fdr_qvalue < input$de_q_val_thresh & 
                           current_results$log2FC < -input$treat_lfc_cutoff, na.rm = TRUE)
            
            paste(
              sprintf("Total genes: %d", total_genes),
              sprintf("Significant (q < %.3f): %d", input$de_q_val_thresh, significant),
              sprintf("Up-regulated: %d", up_reg),
              sprintf("Down-regulated: %d", down_reg),
              sprintf("Fold-change cutoff: %.2f", input$treat_lfc_cutoff),
              sep = "\n"
            )
          } else {
            "No results available for selected contrast"
          }
        } else {
          "No DE analysis results available"
        }
        
      }, error = function(e) {
        paste("Error calculating statistics:", e$message)
      })
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
