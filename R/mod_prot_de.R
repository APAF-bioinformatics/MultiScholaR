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
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 h5 textAreaInput helpText hr numericInput actionButton icon verbatimTextOutput br downloadButton
mod_prot_de_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::wellPanel(
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
        
        # Load filtered session button
        shiny::actionButton(
          ns("load_filtered_session"),
          "Load Filtered Session",
          class = "btn-info",
          width = "100%",
          icon = shiny::icon("upload")
        ),
        shiny::helpText("Load previously exported filtered data for DE analysis"),
        
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
          icon = shiny::icon("chart-area"),
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
  )
  )
}

#' @rdname differentialExpressionAppletModule 
#' @export
mod_prot_de_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns  # Define namespace function for module
    
    cat("--- Entering mod_prot_de_server ---\n")
    cat(sprintf("   mod_prot_de_server Arg: id = %s\n", id))
    cat(sprintf("   mod_prot_de_server Arg: workflow_data is NULL = %s\n", is.null(workflow_data)))
    cat(sprintf("   mod_prot_de_server Arg: selected_tab is NULL = %s\n", is.null(selected_tab)))
    
    cat("=== DIFFERENTIAL EXPRESSION MODULE SERVER STARTED ===\n")
    cat(sprintf("Module ID: %s\n", id))
    
    # Initialize reactive values for DE state
    de_data <- shiny::reactiveValues(
      de_results_list = NULL,
      contrasts_available = NULL,
      analysis_complete = FALSE,
      current_s4_object = NULL,
      formula_from_s4 = NULL,
      current_row_clusters = NULL,
      current_col_clusters = NULL
    )
    
    cat("   mod_prot_de_server Step: Reactive values initialized\n")
    
    # CRITICAL FIX: Add tab selection observer like normalization module
    if (!is.null(selected_tab)) {
      cat("   mod_prot_de_server Step: Setting up tab selection observer\n")
      shiny::observeEvent(selected_tab(), {
        
        cat("--- Entering tab selection observer ---\n")
        cat(sprintf("   tab_observer Step: Selected tab = %s\n", selected_tab()))
        
        # Only trigger if DE tab is selected
        if (!is.null(selected_tab()) && selected_tab() == "de") {
          
          cat("=== DE TAB CLICKED ===\n")
          cat(sprintf("   DE TAB Step: workflow_data$state_manager is NULL = %s\n", is.null(workflow_data$state_manager)))
          
          if (!is.null(workflow_data$state_manager)) {
            current_state <- workflow_data$state_manager$current_state
            
            # Define valid states where this tab can be active
            valid_states_for_de_tab <- c("normalized", "ruv_corrected", "correlation_filtered")
            
            cat(sprintf("   DE TAB Step: Current state = '%s'\n", current_state))
            cat(sprintf("   DE TAB Step: Valid states for DE = %s\n", paste(valid_states_for_de_tab, collapse = ", ")))
            
            # Auto-trigger only fires if we've completed correlation filtering or at least normalization
            if (current_state %in% valid_states_for_de_tab) {
              
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
              cat(sprintf("*** State '%s' is not valid for DE analysis. User needs to complete normalization (with or without RUV) and correlation filtering. ***\n", current_state))
              shiny::showNotification(
                "Please complete the normalization (with or without RUV) and correlation filtering steps before accessing differential expression analysis.",
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
      cat("   mod_prot_de_server Step: No selected_tab parameter provided - tab selection observer NOT set up\n")
    }
    
    # BACKUP: Fetch formula and contrasts from S4 object when state updates (fallback method)
    shiny::observeEvent(workflow_data$state_update_trigger, {
      cat("--- Entering state update trigger observer ---\n")
      cat("\n\n=== DE TAB TRIGGERED VIA STATE UPDATE: Checking for S4 object and contrasts ===\n")
      
      # Get current S4 object from state manager
      if (!is.null(workflow_data$state_manager)) {
        current_state <- workflow_data$state_manager$current_state
        cat(sprintf("   DE TAB Step: Current state = %s\n", current_state))
        
        # Should be getting the correlation-filtered protein object (final state before DE)
        # Also accept normalized or ruv_corrected states
        if (current_state %in% c("normalized", "ruv_corrected", "correlation_filtered")) {
          cat(sprintf("   DE TAB Step: State is valid for DE analysis (%s state found)\n", current_state))
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
          cat(sprintf("   DE TAB Step: State '%s' not valid for DE analysis (expecting: normalized, ruv_corrected, or correlation_filtered)\n", current_state))
        }
      } else {
        cat("   DE TAB Step: workflow_data$state_manager is NULL\n")
      }
      
      cat("=== DE TAB: Contrast detection complete ===\n")
      cat("--- Exiting state update trigger observer ---\n")
    }, ignoreInit = TRUE)
    
    # Update contrast dropdowns when contrasts become available
    shiny::observe({
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
    
    # Load Filtered Session Button Logic
    shiny::observeEvent(input$load_filtered_session, {
      cat("=== LOAD FILTERED SESSION BUTTON CLICKED ===\n")
      
      tryCatch({
        # Check for latest session file
        source_dir <- experiment_paths$source_dir
        if (is.null(source_dir) || !dir.exists(source_dir)) {
          stop("Could not find the source directory to load session data.")
        }
        
        latest_session_file <- file.path(source_dir, "filtered_session_data_latest.rds")
        
        if (!file.exists(latest_session_file)) {
          stop("No exported session data found. Please export session data from the Normalization tab first.")
        }
        
        shiny::withProgress(message = "Loading filtered session data...", value = 0, {
          
          # Step 1: Load session data
          shiny::incProgress(0.3, detail = "Reading session file...")
          session_data <- readRDS(latest_session_file)
          
          cat("*** LOAD: Session data loaded successfully ***\n")
          cat(sprintf("*** LOAD: Export timestamp: %s ***\n", session_data$export_timestamp))
          cat(sprintf("*** LOAD: R6 current state: %s ***\n", session_data$r6_current_state_name))
          cat(sprintf("*** LOAD: Protein count: %d ***\n", session_data$final_protein_count))
          cat(sprintf("*** LOAD: Sample count: %d ***\n", session_data$final_sample_count))
          
          # Step 2: Restore R6 state manager completely
          shiny::incProgress(0.3, detail = "Restoring R6 state...")
          
          # Restore the complete R6 state structure
          workflow_data$state_manager$states <- session_data$r6_complete_states
          workflow_data$state_manager$state_history <- session_data$r6_state_history
          workflow_data$state_manager$current_state <- session_data$r6_current_state_name
          
          cat("*** LOAD: R6 state manager completely restored ***\n")
          cat(sprintf("*** LOAD: Restored %d states ***\n", length(session_data$r6_complete_states)))
          cat(sprintf("*** LOAD: State history: %s ***\n", paste(unlist(session_data$r6_state_history), collapse = " -> ")))
          cat(sprintf("*** LOAD: Current state set to: %s ***\n", session_data$r6_current_state_name))
          
          # Step 3: Restore workflow data
          shiny::incProgress(0.2, detail = "Restoring workflow data...")
          
          workflow_data$contrasts_tbl <- session_data$contrasts_tbl
          workflow_data$design_matrix <- session_data$design_matrix
          workflow_data$config_list <- session_data$config_list
          
          # Set global contrasts_tbl for DE analysis
          if (!is.null(session_data$contrasts_tbl)) {
            assign("contrasts_tbl", session_data$contrasts_tbl, envir = .GlobalEnv)
            cat("*** LOAD: Restored contrasts_tbl to global environment ***\n")
          }
          
          # Set global config_list
          if (!is.null(session_data$config_list)) {
            assign("config_list", session_data$config_list, envir = .GlobalEnv)
            cat("*** LOAD: Restored config_list to global environment ***\n")
          }
          
          # Step 4: Update DE module state
          shiny::incProgress(0.15, detail = "Updating DE module...")
          
          de_data$current_s4_object <- session_data$current_s4_object
          
          # Set contrasts available from loaded session
          if (!is.null(session_data$contrasts_tbl)) {
            if ("contrasts" %in% names(session_data$contrasts_tbl)) {
              de_data$contrasts_available <- session_data$contrasts_tbl$contrasts
            } else {
              de_data$contrasts_available <- session_data$contrasts_tbl[, 1]
            }
          }
          
          # Extract formula from S4 object
          if (!is.null(session_data$current_s4_object) && "deAnalysisParameters" %in% names(session_data$current_s4_object@args)) {
            if ("formula_string" %in% names(session_data$current_s4_object@args$deAnalysisParameters)) {
              formula_from_s4 <- session_data$current_s4_object@args$deAnalysisParameters$formula_string
              de_data$formula_from_s4 <- formula_from_s4
              
              # Update UI with formula from S4
              shiny::updateTextAreaInput(
                session,
                "formula_string",
                value = formula_from_s4
              )
            }
          }
          
          # ✅ NEW: Load UniProt annotations if available
          shiny::incProgress(0.05, detail = "Loading UniProt annotations...")
          
          cat("*** LOAD: Checking for UniProt annotations ***\n")
          tryCatch({
            # Try to load uniprot_dat_cln.RDS from scripts directory
            scripts_uniprot_path <- file.path(source_dir, "uniprot_dat_cln.RDS")
            
            if (file.exists(scripts_uniprot_path)) {
              cat(sprintf("*** LOAD: Found uniprot_dat_cln.RDS at %s ***\n", scripts_uniprot_path))
              
              uniprot_dat_cln <- readRDS(scripts_uniprot_path)
              
              # Store in workflow data and global environment
              workflow_data$uniprot_dat_cln <- uniprot_dat_cln
              assign("uniprot_dat_cln", uniprot_dat_cln, envir = .GlobalEnv)
              
              cat(sprintf("*** LOAD: Successfully loaded %d UniProt annotations ***\n", nrow(uniprot_dat_cln)))
              log_info(paste("Loaded UniProt annotations from", scripts_uniprot_path))
              
              # Show notification to user
              shiny::showNotification(
                sprintf("UniProt annotations loaded: %d protein annotations available for enrichment analysis", nrow(uniprot_dat_cln)),
                type = "message",
                duration = 5
              )
              
            } else {
              cat(sprintf("*** LOAD: No uniprot_dat_cln.RDS found at %s ***\n", scripts_uniprot_path))
              log_info("No UniProt annotations file found - enrichment analysis may have limited functionality")
              
              shiny::showNotification(
                "No UniProt annotations found in session. Enrichment analysis may be limited.",
                type = "warning",
                duration = 5
              )
            }
            
          }, error = function(e) {
            cat(sprintf("*** LOAD: Error loading UniProt annotations: %s ***\n", e$message))
            log_warn(paste("Error loading UniProt annotations:", e$message))
            
            shiny::showNotification(
              paste("Warning: Could not load UniProt annotations:", e$message),
              type = "warning",
              duration = 8
            )
          })
          
          # ✅ NEW: Restore FASTA metadata for complete audit trail
          cat("*** LOAD: Restoring FASTA metadata ***\n")
          if (!is.null(session_data$fasta_metadata)) {
            workflow_data$fasta_metadata <- session_data$fasta_metadata
            cat(sprintf("*** LOAD: FASTA metadata restored (format: %s, sequences: %d) ***\n", 
                       session_data$fasta_metadata$fasta_format,
                       session_data$fasta_metadata$num_sequences))
          } else {
            cat("*** LOAD: No FASTA metadata in session data ***\n")
          }
          
          # ✅ NEW: Restore accession cleanup results
          cat("*** LOAD: Restoring accession cleanup results ***\n")
          if (!is.null(session_data$accession_cleanup_results)) {
            workflow_data$accession_cleanup_results <- session_data$accession_cleanup_results
            cat(sprintf("*** LOAD: Accession cleanup results restored (applied: %s, method: %s) ***\n",
                       session_data$accession_cleanup_results$cleanup_applied,
                       session_data$accession_cleanup_results$aggregation_method))
          } else {
            cat("*** LOAD: No accession cleanup results in session data ***\n")
          }
          
          # ✅ NEW: Restore complete RUV optimization results  
          cat("*** LOAD: Restoring complete RUV optimization results ***\n")
          if (!is.null(session_data$ruv_optimization_result)) {
            workflow_data$ruv_optimization_result <- session_data$ruv_optimization_result
            cat(sprintf("*** LOAD: RUV optimization results restored (k: %d, percentage: %.1f%%) ***\n",
                       session_data$ruv_optimization_result$best_k,
                       session_data$ruv_optimization_result$best_percentage))
          } else {
            cat("*** LOAD: No complete RUV optimization results in session data ***\n")
          }
          
          # ✅ NEW: Restore QC parameters
          cat("*** LOAD: Restoring QC parameters ***\n")
          if (!is.null(session_data$qc_params)) {
            workflow_data$qc_params <- session_data$qc_params
            param_count <- length(unlist(session_data$qc_params, recursive = FALSE))
            cat(sprintf("*** LOAD: QC parameters restored (%d parameter groups) ***\n", param_count))
          } else {
            cat("*** LOAD: No QC parameters in session data ***\n")
          }
          
          # Update tab status
          workflow_data$tab_status$normalization <- "complete"
          workflow_data$tab_status$differential_expression <- "pending"
          
          # Trigger state update for other modules
          workflow_data$state_update_trigger <- Sys.time()
          
        })
        
        # Create summary message
        summary_msg <- sprintf(
          "Session loaded successfully!\n\nData Summary:\n- Proteins: %d\n- Samples: %d\n- Contrasts: %d\n- State: %s\n- Export time: %s\n\nReady for differential expression analysis.",
          session_data$final_protein_count,
          session_data$final_sample_count,
          ifelse(is.null(session_data$contrasts_tbl), 0, nrow(session_data$contrasts_tbl)),
          session_data$r6_current_state_name,
          format(session_data$export_timestamp, "%Y-%m-%d %H:%M:%S")
        )
        
        shiny::showNotification(
          summary_msg,
          type = "message",
          duration = 10
        )
        
        cat("=== LOAD FILTERED SESSION COMPLETED SUCCESSFULLY ===\n")
        
      }, error = function(e) {
        # CRITICAL FIX: Use paste() for logger calls in error handlers to avoid interpolation bug
        error_msg <- paste("Error loading session:", e$message)
        cat(paste("***", error_msg, "\n"))
        log_error(error_msg) 
        
        shiny::showNotification(
          error_msg,
          type = "error",
          duration = 10
        )
      })
    })
    
    # Main DE Analysis Button Logic
    shiny::observeEvent(input$run_de_analysis, {
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
          
          # CRITICAL FIX: Use the correct column for contrast strings
          # The downstream functions expect "comparison=expression" format (from full_format column)
          # NOT just the raw contrast expression (from contrasts column)
          if ("full_format" %in% names(contrasts_tbl)) {
            contrast_strings_to_use <- contrasts_tbl$full_format  # Use full_format column
            cat("   DE ANALYSIS Step: Using full_format column for contrast strings\n")
          } else {
            cat("   DE ANALYSIS Step: No full_format column found, auto-generating from raw contrasts\n")
            # Auto-generate full_format column from raw contrasts
            raw_contrasts <- contrasts_tbl[, 1]
            
            # Generate friendly names and full format
            full_format_strings <- sapply(raw_contrasts, function(contrast_string) {
              # Remove "group" prefixes if present for friendly name
              clean_string <- gsub("^group", "", contrast_string)
              clean_string <- gsub("-group", "-", clean_string)
              
              # Create friendly name by replacing - with _vs_
              friendly_name <- gsub("-", "_vs_", clean_string)
              
              # Create full format: friendly_name=original_contrast_string
              paste0(friendly_name, "=", contrast_string)
            })
            
            contrast_strings_to_use <- full_format_strings
            cat("   DE ANALYSIS Step: Auto-generated full_format strings:\n")
            print(contrast_strings_to_use)
          }
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
              
              # CRITICAL FIX: Ensure the single contrast table uses the correct format
              # Create a modified version that has the full_format in the first column
              if ("full_format" %in% names(single_contrast_tbl)) {
                # Create a temporary table with full_format as the first column for runTestsContrasts
                single_contrast_for_analysis <- data.frame(
                  contrasts = single_contrast_tbl$full_format,
                  stringsAsFactors = FALSE
                )
                cat("   DE ANALYSIS Step: Using full_format for individual contrast analysis\n")
              } else {
                cat("   DE ANALYSIS Step: No full_format available, auto-generating from raw contrast\n")
                # Auto-generate full_format for this single contrast
                raw_contrast <- single_contrast_tbl$contrasts[1]
                
                # Remove "group" prefixes if present for friendly name
                clean_string <- gsub("^group", "", raw_contrast)
                clean_string <- gsub("-group", "-", clean_string)
                
                # Create friendly name by replacing - with _vs_
                friendly_name <- gsub("-", "_vs_", clean_string)
                
                # Create full format: friendly_name=original_contrast_string
                full_format_string <- paste0(friendly_name, "=", raw_contrast)
                
                single_contrast_for_analysis <- data.frame(
                  contrasts = full_format_string,
                  stringsAsFactors = FALSE
                )
                cat(sprintf("   DE ANALYSIS Step: Auto-generated full format: %s\n", full_format_string))
              }
              cat(sprintf("   DE ANALYSIS Step: Single contrast table for index %d:\n", contrast_idx))
              cat("   DE ANALYSIS Step: Structure of single_contrast_for_analysis:\n")
              str(single_contrast_for_analysis)
              cat("   DE ANALYSIS Step: Content of single_contrast_for_analysis:\n")
              print(single_contrast_for_analysis)
              cat(sprintf("   DE ANALYSIS Step: Contrast string being passed: '%s'\n", single_contrast_for_analysis$contrasts[1]))
              
              # Use the modular DE analysis function for this single contrast
              cat("   DE ANALYSIS Step: About to call differentialExpressionAnalysis...\n")
              cat(sprintf("   DE ANALYSIS Step: S4 object class: %s\n", class(de_data$current_s4_object)))
              cat(sprintf("   DE ANALYSIS Step: Formula: %s\n", input$formula_string))
              
              # CRITICAL: Ensure S4 object has the required parameters in @args
              # The differentialExpressionAnalysis function expects parameters in @args$differentialExpressionAnalysis
              if (is.null(de_data$current_s4_object@args$differentialExpressionAnalysis)) {
                cat("   DE ANALYSIS Step: Adding differentialExpressionAnalysis parameters to S4 @args\n")
                de_data$current_s4_object@args$differentialExpressionAnalysis <- list(
                  contrasts_tbl = single_contrast_for_analysis,
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
              
              # ✅ NEW: Store UI parameters in @args for session summary
              cat("   DE ANALYSIS Step: Storing UI parameters in S4 @args\n")
              if (is.null(de_data$current_s4_object@args$deAnalysisUI)) {
                de_data$current_s4_object@args$deAnalysisUI <- list()
              }
              
              de_data$current_s4_object@args$deAnalysisUI <- list(
                q_value_threshold = input$de_q_val_thresh,
                log_fold_change_cutoff = input$treat_lfc_cutoff,
                formula_string = input$formula_string,
                timestamp = Sys.time()
              )
              
              # ✅ NEW: Also store UI parameters in workflow_data for sessionSummary
              workflow_data$de_ui_params <- list(
                q_value_threshold = input$de_q_val_thresh,
                log_fold_change_cutoff = input$treat_lfc_cutoff,
                treat_enabled = input$treat_lfc_cutoff > 0,
                formula_string = input$formula_string,
                timestamp = Sys.time()
              )
              cat("   DE ANALYSIS Step: Stored UI parameters in workflow_data for sessionSummary\n")
              
              # ✅ NEW: Update R6 state manager with UI parameters
              cat("   DE ANALYSIS Step: Updating R6 state with DE UI parameters\n")
              tryCatch({
                # Find the current data state and update it
                current_data_states <- c("correlation_filtered", "normalized", "ruv_corrected", "protein_replicate_filtered")
                available_states <- workflow_data$state_manager$getHistory()
                current_data_state <- purrr::detect(current_data_states, ~ .x %in% available_states)
                
                if (!is.null(current_data_state)) {
                  # CRITICAL FIX: Use the correct 'saveState' method instead of non-existent 'updateState'
                  # Retrieve existing config and description to preserve them
                  existing_state_info <- workflow_data$state_manager$states[[current_data_state]]
                  
                  workflow_data$state_manager$saveState(
                    state_name = current_data_state,
                    s4_data_object = de_data$current_s4_object,
                    config_object = existing_state_info$config,
                    description = existing_state_info$description
                  )
                  cat(sprintf("   DE ANALYSIS Step: Updated state '%s' with DE UI parameters\n", current_data_state))
                }
              }, error = function(e) {
                # CRITICAL FIX: Use paste() for logger calls in error handlers
                cat(paste("   DE ANALYSIS Step: Warning - could not update state with UI parameters:", e$message, "\n"))
              })
              
              # Also add for the helper function
              if (is.null(de_data$current_s4_object@args$differentialExpressionAnalysisHelper)) {
                cat("   DE ANALYSIS Step: Adding differentialExpressionAnalysisHelper parameters to S4 @args\n")
                de_data$current_s4_object@args$differentialExpressionAnalysisHelper <- list(
                  contrasts_tbl = single_contrast_for_analysis,
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
                    contrasts_tbl = single_contrast_for_analysis,
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
                  print(methods::showMethods("differentialExpressionAnalysis", printTo = FALSE))
                  
                  differentialExpressionAnalysis(
                    theObject = de_data$current_s4_object,
                    contrasts_tbl = single_contrast_for_analysis,
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
          
          # Update UI dropdowns to use friendly names that match the data's comparison column
          cat("   DE ANALYSIS Step: Updating UI dropdowns to match data comparison column...\n")
          
          # Use friendly names from contrasts_tbl since they match de_proteins_long$comparison
          if ("friendly_names" %in% names(contrasts_tbl)) {
            friendly_names <- contrasts_tbl$friendly_names
            contrast_choices <- setNames(friendly_names, friendly_names)
            cat(sprintf("   DE ANALYSIS Step: Using friendly_names from contrasts_tbl: %s\n", paste(friendly_names, collapse = ", ")))
          } else {
            # Fallback: extract friendly names from full_format (part before =)
            full_format_strings <- contrasts_tbl$full_format
            friendly_names <- stringr::str_extract(full_format_strings, "^[^=]+")
            contrast_choices <- setNames(friendly_names, friendly_names)
            cat(sprintf("   DE ANALYSIS Step: Extracted friendly names from full_format: %s\n", paste(friendly_names, collapse = ", ")))
          }
          
          shiny::updateSelectInput(session, "volcano_contrast", choices = contrast_choices)
          shiny::updateSelectInput(session, "heatmap_contrast", choices = contrast_choices)
          shiny::updateSelectInput(session, "table_contrast", choices = contrast_choices)
          cat(sprintf("   DE ANALYSIS Step: Updated UI dropdowns with friendly names: %s\n", paste(names(contrast_choices), collapse = ", ")))
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
          
          # Check for qvalue() failures and show prominent warning notification
          all_qvalue_warnings <- list()
          for(contrast_name in names(de_results_list)) {
            if (!is.null(de_results_list[[contrast_name]]$qvalue_warnings) && 
                length(de_results_list[[contrast_name]]$qvalue_warnings) > 0) {
              all_qvalue_warnings <- c(all_qvalue_warnings, de_results_list[[contrast_name]]$qvalue_warnings)
            }
          }
          
          if (length(all_qvalue_warnings) > 0) {
            # Map contrast names to friendly names for user-friendly message
            failed_contrasts <- names(all_qvalue_warnings)
            friendly_failed_names <- c()
            if (exists("contrasts_tbl", envir = .GlobalEnv)) {
              contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
              for (failed_contrast in failed_contrasts) {
                # Extract the part before = if it exists
                contrast_base <- stringr::str_extract(failed_contrast, "^[^=]+")
                if (is.na(contrast_base)) contrast_base <- failed_contrast
                
                # Try to find matching friendly name
                if ("friendly_names" %in% names(contrasts_tbl) && "contrasts" %in% names(contrasts_tbl)) {
                  match_idx <- which(contrasts_tbl$contrasts == contrast_base)
                  if (length(match_idx) > 0) {
                    friendly_failed_names <- c(friendly_failed_names, contrasts_tbl$friendly_names[match_idx[1]])
                  } else {
                    friendly_failed_names <- c(friendly_failed_names, contrast_base)
                  }
                } else {
                  friendly_failed_names <- c(friendly_failed_names, contrast_base)
                }
              }
            } else {
              friendly_failed_names <- failed_contrasts
            }
            
            # Create prominent warning modal dialog
            warning_title <- paste0("⚠️ IMPORTANT: Statistical Analysis Warning")
            warning_body <- paste0(
              "<div style='font-size: 14px; line-height: 1.6;'>",
              "<p><strong>The q-value calculation failed for ", length(all_qvalue_warnings), " contrast(s):</strong></p>",
              "<p style='margin-left: 20px;'>", paste(friendly_failed_names, collapse = ", "), "</p>",
              "<p>The analysis used <strong>Benjamini-Hochberg FDR correction (p.adjust)</strong> instead.</p>",
              "<hr>",
              "<p><strong>⚠️ INTERPRET RESULTS WITH CAUTION:</strong></p>",
              "<ul style='margin-left: 20px;'>",
              "<li>The p-value distribution may be problematic</li>",
              "<li>Results may be less reliable than normal</li>",
              "<li>Consider reviewing your experimental design and data quality</li>",
              "<li>Check diagnostic messages in the console for details</li>",
              "</ul>",
              "</div>"
            )
            
            # Show prominent modal dialog
            shiny::showModal(
              shiny::modalDialog(
                title = shiny::tags$div(shiny::tags$strong(warning_title), style = "color: #d9534f; font-size: 18px;"),
                shiny::HTML(warning_body),
                size = "l",  # Large modal
                easyClose = FALSE,  # User must click button to close
                footer = shiny::tagList(
                  shiny::actionButton(ns("acknowledge_qvalue_warning"), 
                                     "I Understand - Continue", 
                                     class = "btn-warning",
                                     style = "font-weight: bold;")
                )
              )
            )
            
            # Handle acknowledgment button
            shiny::observeEvent(input$acknowledge_qvalue_warning, {
              shiny::removeModal()
            }, once = TRUE)
            
            cat(sprintf("   DE ANALYSIS Step: ⚠️ qvalue() failed for %d contrast(s) - user notification shown\n", length(all_qvalue_warnings)))
          }
          
          # Update workflow data
          workflow_data$de_analysis_results_list <- de_results_list
          workflow_data$tab_status$differential_expression <- "complete"
          workflow_data$tab_status$enrichment_analysis <- "pending"
          
          # ✅ FIXED: Write DE results to disk using new S4 method
          shiny::incProgress(0.9, detail = "Writing results to disk...")
          
          cat("   DE ANALYSIS Step: Writing DE results to disk for enrichment analysis...\n")
          
          tryCatch({
            # Get UniProt annotations for output
            uniprot_tbl <- NULL
            if (exists("uniprot_dat_cln", envir = .GlobalEnv)) {
              uniprot_tbl <- get("uniprot_dat_cln", envir = .GlobalEnv)
              cat("   DE ANALYSIS Step: Found uniprot_dat_cln for annotations\n")
            } else {
              cat("   DE ANALYSIS Step: No uniprot_dat_cln found - proceeding without annotations\n")
            }
            
            # Get the S4 object from the first contrast (they should all be the same)
            first_result <- de_results_list[[1]]
            s4_object <- first_result$theObject
            
            # Use the properly constructed de_output_dir from experiment_paths
            de_output_dir <- experiment_paths$de_output_dir
            cat(sprintf("   DE ANALYSIS Step: Using de_output_dir from experiment_paths: %s\n", de_output_dir))
            
            # Update S4 object parameters for file writing
            s4_object@args$outputDeResultsAllContrasts <- list(
              uniprot_tbl = uniprot_tbl,
              de_output_dir = de_output_dir,
              publication_graphs_dir = experiment_paths$publication_graphs_dir,
              file_prefix = "de_proteins",
              args_row_id = s4_object@protein_id_column,
              gene_names_column = "gene_names",
              uniprot_id_column = "Entry",
              de_q_val_thresh = input$de_q_val_thresh,
              fdr_column = "fdr_qvalue",
              log2fc_column = "log2FC"
            )
            
            # Call the new S4 method that writes all contrasts properly
            cat(sprintf("   DE ANALYSIS Step: Calling outputDeResultsAllContrasts for %d contrasts\n", length(de_results_list)))
            
            success <- outputDeResultsAllContrasts(
              theObject = s4_object,
              de_results_list_all_contrasts = de_results_list,
              uniprot_tbl = uniprot_tbl,
              de_output_dir = de_output_dir,
              publication_graphs_dir = experiment_paths$publication_graphs_dir,
              file_prefix = "de_proteins",
              args_row_id = s4_object@protein_id_column,
              gene_names_column = "gene_names",
              uniprot_id_column = "Entry"
            )
            
            if (success) {
              cat("   DE ANALYSIS Step: All DE results written to disk successfully\n")
            } else {
              cat("   DE ANALYSIS Step: ERROR - outputDeResultsAllContrasts returned FALSE\n")
            }
            
            # ✅ REMOVED: Per-contrast outputDeAnalysisResults calls - not needed!
            # outputDeResultsAllContrasts already handles everything including volcano plots
            cat("   DE ANALYSIS Step: All output handled by outputDeResultsAllContrasts\n")
            
          }, error = function(e) {
            cat(sprintf("   DE ANALYSIS Step: Error writing results to disk: %s\n", e$message))
            # Don't fail the entire analysis for file writing issues
            # NOTE: Removed log_warn call to avoid logger interpolation bug in error handlers
            message(paste("Could not write DE results to disk:", e$message))
          })
          
          shiny::incProgress(1.0, detail = "Complete!")
        })
        
        shiny::showNotification(
          "Differential expression analysis completed successfully!",
          type = "message",
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
        # Initialize glimma_widget to prevent "object not found" errors
        glimma_widget <- NULL
        
        # Find the correct contrast results using friendly name mapping
        contrast_keys <- names(de_data$de_results_list$individual_contrasts)
        
        # Map friendly names from contrasts_tbl to the actual stored contrast keys
        selected_contrast_results <- NULL
        full_contrast_name <- NULL
        
        if (exists("contrasts_tbl", envir = .GlobalEnv)) {
          contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
          
          # Find the row in contrasts_tbl that matches the selected friendly name
          if ("friendly_names" %in% names(contrasts_tbl)) {
            matching_row <- which(contrasts_tbl$friendly_names == input$volcano_contrast)
            if (length(matching_row) > 0) {
              # Get the raw contrast name that corresponds to this friendly name
              raw_contrast <- contrasts_tbl$contrasts[matching_row[1]]
              
              # Find the stored result that matches this raw contrast
              matching_key <- which(stringr::str_detect(contrast_keys, raw_contrast))
              if (length(matching_key) > 0) {
                full_contrast_name <- contrast_keys[matching_key[1]]
                selected_contrast_results <- de_data$de_results_list$individual_contrasts[[full_contrast_name]]
                cat(sprintf("   VOLCANO: Found matching contrast results for %s\n", full_contrast_name))
                cat(sprintf("   VOLCANO: Mapped friendly name '%s' to stored key '%s'\n", input$volcano_contrast, full_contrast_name))
              }
            }
          }
        }
        
        if (is.null(selected_contrast_results)) {
          cat(sprintf("   VOLCANO: No contrast results found for %s\n", input$volcano_contrast))
          cat(sprintf("   VOLCANO: Available stored keys: %s\n", paste(contrast_keys, collapse = ", ")))
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
          cat(sprintf("   VOLCANO: Using stored contrast key for results = %s\n", full_contrast_name))
          
          # Fix: Use the correct protein ID column from the S4 object
          protein_id_column <- de_data$de_results_list$theObject@protein_id_column
          cat(sprintf("   VOLCANO: Using protein ID column = %s\n", protein_id_column))
          
          # CRITICAL FIX: Pass the friendly name to generateVolcanoPlotGlimma since that's what the data contains
          glimma_widget <- tryCatch({
            generateVolcanoPlotGlimma(
              de_results_list = plot_data_structure,
              selected_contrast = input$volcano_contrast,  # Use friendly name that matches de_proteins_long$comparison
              uniprot_tbl = uniprot_tbl,
              de_q_val_thresh = input$de_q_val_thresh,
              args_row_id = protein_id_column  # Pass the correct column name
            )
          }, error = function(e) {
            cat(sprintf("   VOLCANO: Error in generateVolcanoPlotGlimma: %s\n", e$message))
            NULL
          })
          
          cat(sprintf("   VOLCANO: glimma_widget created successfully: %s\n", !is.null(glimma_widget)))
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

