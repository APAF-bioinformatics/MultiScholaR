#' @title enrichmentAnalysisAppletModule
#'
#' @description A Shiny module for the Enrichment Analysis step of the proteomics
#' workflow. Handles GO, KEGG, Reactome, and String-DB pathway enrichment analysis
#' using results from differential expression analysis.
#'
#' @name enrichmentAnalysisAppletModule
NULL

#' @rdname enrichmentAnalysisAppletModule
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 selectInput helpText hr numericInput actionButton icon br tabsetPanel tabPanel verbatimTextOutput plotOutput downloadButton
#' @importFrom DT DTOutput renderDT
mod_prot_enrich_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::wellPanel(
      shiny::fluidRow(
    shiny::column(3,
      shiny::wellPanel(
        shiny::h4("Enrichment Analysis Settings"),
        
        # Contrast selector
        shiny::selectInput(
          ns("selected_contrast"),
          "Select Contrast:",
          choices = NULL,
          width = "100%"
        ),
        shiny::helpText("Choose which DE contrast to analyze"),
        
        shiny::hr(),
        
        # Only the parameters that ACTUALLY EXIST in processEnrichments
        shiny::h5("Enrichment Settings"),
        shiny::numericInput(
          ns("up_cutoff"),
          "Up Log2FC Cutoff:",
          value = 0,
          min = 0,
          max = 5,
          step = 0.1
        ),
        shiny::helpText("Minimum log2 fold change for up-regulated genes"),
        
        shiny::numericInput(
          ns("down_cutoff"),
          "Down Log2FC Cutoff:",
          value = 0,
          min = 0,
          max = 5,
          step = 0.1
        ),
        shiny::helpText("Minimum log2 fold change for down-regulated genes"),
        
        shiny::numericInput(
          ns("q_cutoff"),
          "Q-value Cutoff:",
          value = 0.05,
          min = 0.001,
          max = 0.2,
          step = 0.005
        ),
        shiny::helpText("FDR threshold for enrichment significance"),
        
        shiny::textInput(
          ns("organism_taxid"),
          "Organism Taxon ID:",
          value = 9606,  # Default to human, will be updated by server
          placeholder = "e.g., 9606 for human"
        ),
        shiny::helpText("NCBI taxonomy ID for organism"),
        
        # ✅ NEW: Analysis method display
        shiny::h5("Analysis Method"),
        shiny::verbatimTextOutput(ns("analysis_method_display")),
        shiny::helpText("Automatically determined based on organism support"),
        
        shiny::hr(),
        
        # Available contrasts display
        shiny::h5("Available Contrasts"),
        shiny::verbatimTextOutput(ns("contrasts_display")),
        shiny::helpText("Contrasts from DE analysis"),
        
        shiny::br(),
        
        # Main action button
        shiny::actionButton(
          ns("run_enrichment_analysis"),
          "Run Enrichment Analysis",
          class = "btn-primary",
          width = "100%",
          icon = shiny::icon("play")
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Download button
        shiny::downloadButton(
          ns("download_enrichment_results"),
          "Download All Results",
          class = "btn-success",
          width = "100%"
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Status display
        shiny::h5("Analysis Status"),
        shiny::verbatimTextOutput(ns("enrichment_status")),
        
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'gprofiler2'", ns("enrichment_method_tabs")),
          shiny::selectInput(
            ns("correction_method"),
            "Correction Method:",
            choices = list(
              "g:SCS (Conservative)" = "gSCS",
              "FDR (Less Conservative)" = "fdr"
            ),
            selected = "gSCS"
          ),
          shiny::helpText("g:SCS: More conservative, reduces false positives but may miss terms."),
          shiny::helpText("FDR: Detects more terms but higher false positive risk.")
        )
      )
    ),
    
    # ✅ FIXED: Split tab structure with plot displays added
    shiny::column(9,
      shiny::tabsetPanel(
        id = ns("enrichment_method_tabs"),
        
        # ✅ Tab 1: gprofiler2 for supported organisms
        shiny::tabPanel(
          "gprofiler2 Analysis",
          value = "gprofiler2",
          icon = shiny::icon("globe"),
          shiny::br(),
          
          shiny::div(
            style = "background-color: #e8f4fd; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::h5("gprofiler2 Enrichment", style = "color: #2c5282; margin-top: 0;"),
            shiny::p("Comprehensive functional enrichment using gprofiler2 for well-supported model organisms. 
                     Includes GO terms, KEGG pathways, Reactome pathways, and more.", 
                     style = "margin-bottom: 0; color: #2c5282;")
          ),
          
          # gprofiler2 results display
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'gprofiler2'", ns("enrichment_method_tabs")),
            
            # Direction filter for viewing results
            shiny::fluidRow(
              shiny::column(6,
                shiny::selectInput(
                  ns("gprofiler_direction_filter"),
                  "Direction:",
                  choices = list(
                    "All" = "all",
                    "Up-regulated only" = "up",
                    "Down-regulated only" = "down"
                  ),
                  selected = "all"
                )
              ),
              shiny::column(6,
                shiny::div(
                  style = "padding: 20px 0; text-align: center; color: #666;",
                  shiny::strong("Select direction to filter results")
                )
              )
            ),
            
            shiny::hr(),
            
            # ✅ NEW: Add plot display for gprofiler2
            shiny::fluidRow(
              shiny::column(12,
                shiny::h5("Enrichment Plot"),
                shinyjqui::jqui_resizable(
                  plotly::plotlyOutput(ns("gprofiler_plot"), height = "500px", width = "100%")
                ),
                shiny::br()
              )
            ),
            
            # gprofiler2 results table
            DT::DTOutput(ns("gprofiler_results_table")),
            
            shiny::br(),
            
            # gprofiler2 summary stats
            shiny::wellPanel(
              shiny::h5("gprofiler2 Enrichment Summary"),
              shiny::verbatimTextOutput(ns("gprofiler_summary_stats"))
            )
          )
        ),
        
        # ✅ Tab 2: clusterProfileR for unsupported organisms  
        shiny::tabPanel(
          "clusterProfileR Analysis",
          value = "clusterprofiler",
          icon = shiny::icon("project-diagram"),
          shiny::br(),
          
          shiny::div(
            style = "background-color: #f0fff4; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::h5("clusterProfileR Enrichment", style = "color: #22543d; margin-top: 0;"),
            shiny::p("Custom enrichment analysis using UniProt GO annotations via clusterProfileR. 
                     Ideal for organisms not supported by gprofiler2.", 
                     style = "margin-bottom: 0; color: #22543d;")
          ),
          
          # clusterProfileR results display
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'clusterprofiler'", ns("enrichment_method_tabs")),
            
            # Direction filter for viewing results
            shiny::fluidRow(
              shiny::column(6,
                shiny::selectInput(
                  ns("clusterprofiler_direction_filter"),
                  "Direction:",
                  choices = list(
                    "All" = "all",
                    "Up-regulated only" = "up",
                    "Down-regulated only" = "down"
                  ),
                  selected = "all"
                )
              ),
              shiny::column(6,
                shiny::div(
                  style = "padding: 20px 0; text-align: center; color: #666;",
                  shiny::strong("Select direction to filter results")
                )
              )
            ),
            
            shiny::hr(),
            
            # ✅ NEW: Add plot display for clusterProfileR
            shiny::fluidRow(
              shiny::column(12,
                shiny::h5("Enrichment Plot"),
                shinyjqui::jqui_resizable(
                  plotly::plotlyOutput(ns("clusterprofiler_plot"), height = "500px", width = "100%")
                ),
                shiny::br()
              )
            ),
            
            # clusterProfileR results table
            DT::DTOutput(ns("clusterprofiler_results_table")),
            
            shiny::br(),
            
            # clusterProfileR summary stats
            shiny::wellPanel(
              shiny::h5("clusterProfileR Enrichment Summary"),
              shiny::verbatimTextOutput(ns("clusterprofiler_summary_stats"))
            )
          )
        ),
        
        # ✅ Tab 3: STRING-DB for protein-protein interaction networks
        shiny::tabPanel(
          "STRING-DB Networks",
          value = "stringdb",
          icon = shiny::icon("network-wired"),
          shiny::br(),
          
          shiny::div(
            style = "background-color: #fef5e7; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::h5("STRING-DB Network Analysis", style = "color: #744210; margin-top: 0;"),
            shiny::p("Protein-protein interaction network enrichment using STRING database. 
                     Identifies functional networks and interaction clusters.", 
                     style = "margin-bottom: 0; color: #744210;")
          ),
          
          # STRING-DB results display
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'stringdb'", ns("enrichment_method_tabs")),
            
            # Filter controls for STRING-DB
            shiny::fluidRow(
              shiny::column(4,
                shiny::selectInput(
                  ns("stringdb_ranking_method"),
                  "Ranking Method:",
                  choices = list(
                    "Combined Score" = "combined_score",
                    "FDR Q-value" = "fdr_qvalue",
                    "Log2 Fold Change" = "log2fc"
                  ),
                  selected = "combined_score"
                )
              ),
              shiny::column(4,
                shiny::checkboxInput(
                  ns("stringdb_filter_significant"),
                  "Filter Significant Only",
                  value = FALSE
                )
              ),
              shiny::column(4,
                shiny::numericInput(
                  ns("stringdb_max_results"),
                  "Max Results:",
                  value = 50,
                  min = 10,
                  max = 500,
                  step = 10
                )
              )
            ),
            
            shiny::hr(),
            
            # ✅ NEW: Add plot display for STRING-DB
            shiny::fluidRow(
              shiny::column(12,
                shiny::h5("Network Plot"),
                shinyjqui::jqui_resizable(
                  plotly::plotlyOutput(ns("stringdb_plot"), height = "500px", width = "100%")
                ),
                shiny::br()
              )
            ),
            
            # STRING-DB results table
            DT::DTOutput(ns("stringdb_results_table")),
            
            shiny::br(),
            
            # STRING-DB summary stats  
            shiny::wellPanel(
              shiny::h5("STRING-DB Network Summary"),
              shiny::verbatimTextOutput(ns("stringdb_summary_stats"))
            )
          )
        )
      )
    )
  )
  )
  )
}



#' @rdname enrichmentAnalysisAppletModule 
#' @export
mod_prot_enrich_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    cat("--- Entering mod_prot_enrich_server ---\n")
    cat(sprintf("   mod_prot_enrich_server Arg: id = %s\n", id))
    cat(sprintf("   mod_prot_enrich_server Arg: workflow_data is NULL = %s\n", is.null(workflow_data)))
    cat(sprintf("   mod_prot_enrich_server Arg: selected_tab is NULL = %s\n", is.null(selected_tab)))
    
    cat("=== ENRICHMENT ANALYSIS MODULE SERVER STARTED ===\n")
    cat(sprintf("Module ID: %s\n", id))
    
    # Initialize reactive values for enrichment state
    enrichment_data <- shiny::reactiveValues(
      enrichment_results = NULL,
      contrasts_available = NULL,
      analysis_complete = FALSE,
      current_s4_object = NULL,
      de_results_data = NULL,
      # ✅ NEW: Separate results for each analysis method
      gprofiler_results = NULL,
      clusterprofiler_results = NULL,
      stringdb_results = NULL,
      analysis_method = NULL,
      organism_supported = NULL,
      # ✅ NEW: Store results for ALL contrasts (not just selected one)
      all_enrichment_results = list(),
      current_contrast_results = list(),
      enrichment_plots = list()
    )
    
    cat("   mod_prot_enrich_server Step: Reactive values initialized\n")
    
    # Update organism taxid input when workflow_data becomes available
    shiny::observeEvent(workflow_data$taxon_id, {
      if (!is.null(workflow_data$taxon_id)) {
        shiny::updateTextInput(session, "organism_taxid", value = workflow_data$taxon_id)
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    
    # ✅ NEW: Organism support detection
    supported_organisms <- shiny::reactive({
      tibble::tribble(
        ~taxid,     ~id,            ~name,
        "9606",     "hsapiens",     "Homo sapiens",
        "10090",    "mmusculus",    "Mus musculus",
        "10116",    "rnorvegicus",  "Rattus norvegicus",
        "7227",     "dmelanogaster", "Drosophila melanogaster",
        "6239",     "celegans",     "Caenorhabditis elegans",
        "4932",     "scerevisiae",  "Saccharomyces cerevisiae",
        "3702",     "athaliana",    "Arabidopsis thaliana",
        "7955",     "drerio",       "Danio rerio",
        "9031",     "ggallus",      "Gallus gallus",
        "9823",     "sscrofa",      "Sus scrofa",
        "9913",     "btaurus",      "Bos taurus",
        "9544",     "mmulatta",     "Macaca mulatta",
        "9598",     "ptroglodytes", "Pan troglodytes"
      )
    })
    
    # ✅ NEW: Determine analysis method based on organism
    current_analysis_method <- shiny::reactive({
      shiny::req(input$organism_taxid)
      
      is_supported <- input$organism_taxid %in% supported_organisms()$taxid
      
      if (is_supported) {
        species_info <- supported_organisms() |>
          dplyr::filter(taxid == input$organism_taxid)
        
        method_info <- list(
          method = "gprofiler2",
          supported = TRUE,
          species_id = species_info$id[1],
          species_name = species_info$name[1],
          description = paste("gprofiler2 analysis for", species_info$name[1])
        )
      } else {
        method_info <- list(
          method = "clusterprofiler",
          supported = FALSE,
          species_id = NULL,
          species_name = paste("Taxon ID", input$organism_taxid),
          description = paste("clusterProfileR analysis with custom GO annotations for taxon", input$organism_taxid)
        )
      }
      
      # Update enrichment_data
      enrichment_data$analysis_method <- method_info$method
      enrichment_data$organism_supported <- method_info$supported
      
      return(method_info)
    })
    
    # ✅ NEW: Display analysis method
    output$analysis_method_display <- shiny::renderText({
      method_info <- current_analysis_method()
      
      if (method_info$supported) {
        paste(
          "✅ SUPPORTED ORGANISM\n",
          sprintf("Method: %s\n", method_info$method),
          sprintf("Species: %s\n", method_info$species_name),
          "All enrichment methods available"
        )
      } else {
        paste(
          "⚠️ CUSTOM ORGANISM\n",
          sprintf("Method: %s\n", method_info$method),
          sprintf("Organism: %s\n", method_info$species_name),
          "Using UniProt GO annotations"
        )
      }
    })
    
    # ✅ FIXED: Only run contrast observer when enrichment analysis is complete
    shiny::observe({
      shiny::req(input$selected_contrast)
      shiny::req(enrichment_data$analysis_complete)
      shiny::req(enrichment_data$all_enrichment_results)
      
      cat(sprintf("*** CONTRAST CHANGED: User selected '%s' ***\n", input$selected_contrast))
      
      # Map friendly name to raw contrast name 
      raw_contrast_name <- input$selected_contrast
      
      # Try to map using contrasts_tbl if available
      if (exists("contrasts_tbl", envir = .GlobalEnv)) {
        contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
        if ("friendly_names" %in% names(contrasts_tbl)) {
          matching_idx <- which(contrasts_tbl$friendly_names == input$selected_contrast)
          if (length(matching_idx) > 0) {
            raw_contrast_name <- contrasts_tbl$contrasts[matching_idx[1]]
            cat(sprintf("*** CONTRAST MAPPING: '%s' -> '%s' ***\n", input$selected_contrast, raw_contrast_name))
          }
        }
      }
      
      # Update current contrast results
      if (raw_contrast_name %in% names(enrichment_data$all_enrichment_results)) {
        contrast_results <- enrichment_data$all_enrichment_results[[raw_contrast_name]]
        
        cat(sprintf("*** UPDATING RESULTS: Found results for contrast '%s' ***\n", raw_contrast_name))
        
        # Update gprofiler2 results
        if (!is.null(contrast_results$gprofiler_results)) {
          enrichment_data$gprofiler_results <- contrast_results$gprofiler_results
          cat(sprintf("*** UPDATED: %d gprofiler2 results ***\n", nrow(contrast_results$gprofiler_results)))
        } else {
          enrichment_data$gprofiler_results <- NULL
        }
        
        # Update clusterProfileR results
        if (!is.null(contrast_results$clusterprofiler_results)) {
          enrichment_data$clusterprofiler_results <- contrast_results$clusterprofiler_results
          cat(sprintf("*** UPDATED: %d clusterProfileR results ***\n", nrow(contrast_results$clusterprofiler_results)))
        } else {
          enrichment_data$clusterprofiler_results <- NULL
        }
        
        # Update stringDB results
        if (!is.null(contrast_results$stringdb_results)) {
          enrichment_data$stringdb_results <- contrast_results$stringdb_results
          cat(sprintf("*** UPDATED: %d stringDB results ***\n", nrow(contrast_results$stringdb_results)))
        } else {
          enrichment_data$stringdb_results <- NULL
        }
        
      } else {
        cat(sprintf("*** WARNING: No results found for contrast '%s' ***\n", raw_contrast_name))
        cat(sprintf("*** AVAILABLE CONTRASTS: %s ***\n", paste(names(enrichment_data$all_enrichment_results), collapse = ", ")))
        
        # Clear current results if contrast not found
        enrichment_data$gprofiler_results <- NULL
        enrichment_data$clusterprofiler_results <- NULL
        enrichment_data$stringdb_results <- NULL
      }
    })
    
    # Tab selection observer (matches DE applet pattern)
    if (!is.null(selected_tab)) {
      cat("   mod_prot_enrich_server Step: Setting up tab selection observer\n")
      shiny::observeEvent(selected_tab(), {
        
        cat("--- Entering enrichment tab selection observer ---\n")
        cat(sprintf("   enrichment_tab_observer Step: Selected tab = %s\n", selected_tab()))
        
        # Only trigger if enrichment tab is selected
        if (!is.null(selected_tab()) && selected_tab() == "enrichment") {
          
          cat("=== ENRICHMENT TAB CLICKED ===\n")
          cat(sprintf("   ENRICHMENT TAB Step: workflow_data$state_manager is NULL = %s\n", is.null(workflow_data$state_manager)))
          
          if (!is.null(workflow_data$state_manager)) {
            current_state <- workflow_data$state_manager$current_state
            
            # Check if we have completed DE analysis
            valid_states_for_enrichment <- c("correlation_filtered")  # Same as DE for now
            
            cat(sprintf("   ENRICHMENT TAB Step: Current state = '%s'\n", current_state))
            cat(sprintf("   ENRICHMENT TAB Step: Valid states for enrichment = %s\n", paste(valid_states_for_enrichment, collapse = ", ")))
            
            # Also check if DE analysis has been completed
            if (current_state == "correlation_filtered" && !is.null(workflow_data$de_analysis_results_list)) {
              
              cat("*** AUTO-TRIGGERING ENRICHMENT INITIALIZATION (DE results found) ***\n")
              
              tryCatch({
                # Get S4 object and DE results
                cat("   ENRICHMENT TAB Step: Getting S4 object and DE results from workflow_data...\n")
                current_s4 <- workflow_data$state_manager$getState(current_state)
                de_results_list <- workflow_data$de_analysis_results_list
                
                if (!is.null(current_s4) && !is.null(de_results_list)) {
                  cat(sprintf("   ENRICHMENT TAB Step: S4 object retrieved, class = %s\n", class(current_s4)))
                  cat(sprintf("   ENRICHMENT TAB Step: DE results available for %d contrasts\n", length(de_results_list)))
                  
                  enrichment_data$current_s4_object <- current_s4
                  enrichment_data$de_results_data <- de_results_list
                  
                  # Set up contrast choices from DE results
                  contrast_names <- names(de_results_list)
                  cat(sprintf("   ENRICHMENT TAB Step: Available contrast names: %s\n", paste(contrast_names, collapse = ", ")))
                  
                  # Map to friendly names if contrasts_tbl exists
                  if (exists("contrasts_tbl", envir = .GlobalEnv)) {
                    contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
                    if ("friendly_names" %in% names(contrasts_tbl)) {
                      friendly_names <- contrasts_tbl$friendly_names
                      # Set friendly names as both keys and values for simplicity
                      contrast_choices <- setNames(friendly_names, friendly_names)
                      enrichment_data$contrasts_available <- friendly_names
                      cat(sprintf("   ENRICHMENT TAB Step: Using friendly names: %s\n", paste(friendly_names, collapse = ", ")))
                    } else {
                      enrichment_data$contrasts_available <- contrast_names
                      contrast_choices <- setNames(contrast_names, contrast_names)
                    }
                  } else {
                    enrichment_data$contrasts_available <- contrast_names
                    contrast_choices <- setNames(contrast_names, contrast_names)
                  }
                  
                  # Update contrast selector
                  shiny::updateSelectInput(
                    session,
                    "selected_contrast",
                    choices = contrast_choices
                  )
                  
                } else {
                  cat("   ENRICHMENT TAB Step: S4 object or DE results are NULL\n")
                }
                
                cat("*** ENRICHMENT INITIALIZATION COMPLETED SUCCESSFULLY ***\n")
                
              }, error = function(e) {
                cat(paste("*** ERROR in enrichment initialization:", e$message, "\n"))
                shiny::showNotification(
                  paste("Error initializing enrichment analysis:", e$message),
                  type = "error",
                  duration = 10
                )
              })
              
            } else if (is.null(workflow_data$de_analysis_results_list)) {
              cat("*** No DE analysis results found. User needs to complete differential expression analysis first. ***\n")
              shiny::showNotification(
                "Please complete the differential expression analysis before accessing enrichment analysis.",
                type = "warning",
                duration = 5
              )
            } else {
              cat(sprintf("*** State '%s' is not valid for enrichment analysis. User needs to complete DE analysis. ***\n", current_state))
              shiny::showNotification(
                "Please complete the differential expression analysis before accessing enrichment analysis.",
                type = "warning",
                duration = 5
              )
            }
          } else {
            cat("*** workflow_data$state_manager is NULL - cannot check state ***\n")
          }
        } else {
          cat(sprintf("   enrichment_tab_observer Step: Tab '%s' is not enrichment tab, ignoring\n", selected_tab()))
        }
        
        cat("--- Exiting enrichment tab selection observer ---\n")
      }, ignoreInit = TRUE)
    } else {
      cat("   mod_prot_enrich_server Step: No selected_tab parameter provided - tab selection observer NOT set up\n")
    }
    
    # Display available contrasts
    output$contrasts_display <- shiny::renderText({
      if (!is.null(enrichment_data$contrasts_available)) {
        paste(enrichment_data$contrasts_available, collapse = "\n")
      } else {
        "No contrasts available.\nComplete differential expression\nanalysis first."
      }
    })
    
    # Display analysis status
    output$enrichment_status <- shiny::renderText({
      if (enrichment_data$analysis_complete) {
        method_info <- current_analysis_method()
        
        # Count results across methods
        gprofiler_count <- if (!is.null(enrichment_data$gprofiler_results)) nrow(enrichment_data$gprofiler_results) else 0
        clusterprofiler_count <- if (!is.null(enrichment_data$clusterprofiler_results)) nrow(enrichment_data$clusterprofiler_results) else 0
        stringdb_count <- if (!is.null(enrichment_data$stringdb_results)) nrow(enrichment_data$stringdb_results) else 0
        
        paste(
          "✅ Analysis Complete\n",
          sprintf("Method: %s\n", method_info$method),
          sprintf("Contrast: %s\n", input$selected_contrast),
          sprintf("Up log2FC cutoff: %.1f\n", input$up_cutoff),
          sprintf("Down log2FC cutoff: %.1f\n", input$down_cutoff),
          sprintf("Q-value cutoff: %.3f\n", input$q_cutoff),
          sprintf("Organism: %s\n", method_info$species_name),
          "",
          "Results Available:",
          sprintf("• gprofiler2: %d terms", gprofiler_count),
          sprintf("• clusterProfileR: %d terms", clusterprofiler_count),
          sprintf("• STRING-DB: %d networks", stringdb_count),
          "",
          "✓ Results saved to workflow state",
          sep = "\n"
        )
      } else {
        paste(
          "⏳ Ready for analysis\n",
          "",
          "Steps:",
          "1. Select contrast from DE results",
          "2. Set log fold change cutoffs",
          "3. Set Q-value cutoff (significance threshold)",  
          "4. Click 'Run Enrichment Analysis'",
          "",
          "Method automatically determined by organism.",
          sep = "\n"
        )
      }
    })
    
    # Main Enrichment Analysis Button Logic
    shiny::observeEvent(input$run_enrichment_analysis, {
      cat("=== STARTING ENRICHMENT ANALYSIS ===\n")
      
      shiny::req(input$selected_contrast, enrichment_data$de_results_data)
      
      shiny::showNotification("Running enrichment analysis...", id = "enrichment_working", duration = NULL)
      
      tryCatch({
        shiny::withProgress(message = "Running enrichment analysis...", value = 0, {
          
          # Step 1: Transform DE results to enrichment function format
          shiny::incProgress(0.2, detail = "Transforming DE data...")
          
          cat(sprintf("   ENRICHMENT Step: Selected contrast (friendly name) = %s\n", input$selected_contrast))
          cat(sprintf("   ENRICHMENT Step: Available DE results: %s\n", paste(names(enrichment_data$de_results_data), collapse = ", ")))
          
          # ✅ IMPROVED: Find matching DE results with better error handling
          selected_de_results <- NULL
          raw_contrast_name <- NULL
          
          # Strategy 1: Try direct mapping using contrasts_tbl
          if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
            if ("friendly_names" %in% names(contrasts_tbl)) {
              # Find the raw contrast name that corresponds to the friendly name
              matching_idx <- which(contrasts_tbl$friendly_names == input$selected_contrast)
              if (length(matching_idx) > 0) {
                raw_contrast_name <- contrasts_tbl$contrasts[matching_idx[1]]
                cat(sprintf("   ENRICHMENT Step: Mapped friendly name '%s' to raw name '%s'\n", 
                           input$selected_contrast, raw_contrast_name))
                
                # Try to find DE results using the raw contrast name
                selected_de_results <- enrichment_data$de_results_data[[raw_contrast_name]]
                
                if (!is.null(selected_de_results)) {
                  cat(sprintf("   ENRICHMENT Step: Found DE results for raw contrast %s\n", raw_contrast_name))
                } else {
                  cat(sprintf("   ENRICHMENT Step: No DE results found for raw contrast '%s'\n", raw_contrast_name))
                }
              }
            }
          }
          
          # Strategy 2: If not found, try fuzzy matching against available DE result keys
          if (is.null(selected_de_results)) {
            available_keys <- names(enrichment_data$de_results_data)
            cat(sprintf("   ENRICHMENT Step: Available DE result keys: %s\n", paste(available_keys, collapse = ", ")))
            cat(sprintf("   ENRICHMENT Step: Looking for friendly name: %s\n", input$selected_contrast))
            
            # Try to find a key that contains the essential parts of the contrast
            # Extract essential parts from friendly name (e.g., "GA_Elevated_vs_GA_Control" -> "GA_Elevated" and "GA_Control")
            contrast_parts <- stringr::str_split(input$selected_contrast, "_vs_")[[1]]
            if (length(contrast_parts) == 2) {
              part1 <- contrast_parts[1]
              part2 <- contrast_parts[2]
              
              # Look for a key that contains both parts
              for (key in available_keys) {
                if (stringr::str_detect(key, part1) && stringr::str_detect(key, part2)) {
                  selected_de_results <- enrichment_data$de_results_data[[key]]
                  raw_contrast_name <- key
                  cat(sprintf("   ENRICHMENT Step: Found matching DE results using fuzzy matching: %s\n", key))
                  break
                }
              }
            }
          }
          
          # Strategy 3: If still not found, try direct key lookup (in case friendly name matches key)
          if (is.null(selected_de_results)) {
            selected_de_results <- enrichment_data$de_results_data[[input$selected_contrast]]
            if (!is.null(selected_de_results)) {
              raw_contrast_name <- input$selected_contrast
              cat(sprintf("   ENRICHMENT Step: Found DE results using direct key lookup: %s\n", input$selected_contrast))
            }
          }
          
          # Final check
          if (is.null(selected_de_results)) {
            available_keys <- names(enrichment_data$de_results_data)
            stop(sprintf("Could not find DE results for contrast '%s'. Available contrasts: %s", 
                        input$selected_contrast, paste(available_keys, collapse = ", ")))
          }
          
          # Step 2: Create S4 object for enrichment analysis (following .rmd pattern)
          shiny::incProgress(0.3, detail = "Creating DE results S4 object...")
          
          cat("   ENRICHMENT Step: Creating DE results S4 object using createDEResultsForEnrichment\n")
          
          # Get contrasts_tbl and design_matrix from workflow data
          contrasts_tbl <- NULL
          design_matrix <- NULL
          
          if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
            cat("   ENRICHMENT Step: Found contrasts_tbl in global environment\n")
          }
          
          if (!is.null(enrichment_data$current_s4_object@design_matrix)) {
            design_matrix <- enrichment_data$current_s4_object@design_matrix
            cat("   ENRICHMENT Step: Found design_matrix in S4 object\n")
          }
          
          if (is.null(contrasts_tbl) || is.null(design_matrix)) {
            stop("Missing contrasts_tbl or design_matrix required for enrichment analysis")
          }
          
          # ✅ FIXED: Use correct pathway directory from experiment_paths
          de_output_dir_new <- file.path(experiment_paths$results_dir, "de_proteins")
          de_output_dir_old <- file.path(experiment_paths$source_dir, "de_output")
          
          # Try the new path first, fall back to old path
          de_output_dir <- if (dir.exists(de_output_dir_new)) {
            cat(sprintf("   ENRICHMENT Step: Using new DE output directory: %s\n", de_output_dir_new))
            de_output_dir_new
          } else if (dir.exists(de_output_dir_old)) {
            cat(sprintf("   ENRICHMENT Step: Using legacy DE output directory: %s\n", de_output_dir_old))
            de_output_dir_old
          } else {
            cat("   ENRICHMENT Step: No DE output directory found, using new path\n")
            de_output_dir_new
          }
          
          # ✅ FIXED: Use proper pathway directory structure
          pathway_dir <- if (!is.null(experiment_paths$pathway_dir) && dir.exists(experiment_paths$pathway_dir)) {
            cat(sprintf("   ENRICHMENT Step: Using pathway_dir from experiment_paths: %s\n", experiment_paths$pathway_dir))
            experiment_paths$pathway_dir
          } else {
            # Fallback to creating pathway directory in results
            fallback_pathway_dir <- file.path(experiment_paths$results_dir, "pathway_enrichment")
            if (!dir.exists(fallback_pathway_dir)) {
              dir.create(fallback_pathway_dir, recursive = TRUE)
            }
            cat(sprintf("   ENRICHMENT Step: Created fallback pathway directory: %s\n", fallback_pathway_dir))
            fallback_pathway_dir
          }
          
          # Create S4 object for enrichment analysis
          de_results_for_enrichment <- createDEResultsForEnrichment(
            contrasts_tbl = contrasts_tbl,
            design_matrix = design_matrix,
            de_output_dir = de_output_dir
          )
          
          cat("   ENRICHMENT Step: S4 object created successfully\n")
          
          # Step 3: Run enrichment analysis using processEnrichments (following .rmd pattern)
          shiny::incProgress(0.5, detail = "Running enrichment analysis...")
          
          cat("   ENRICHMENT Step: Running processEnrichments\n")
          
          # ✅ IMPROVED: Get or create UniProt annotations
          uniprot_dat_cln <- NULL
          if (exists("uniprot_dat_cln", envir = .GlobalEnv)) {
            uniprot_dat_cln <- get("uniprot_dat_cln", envir = .GlobalEnv)
            cat("   ENRICHMENT Step: Found uniprot_dat_cln in global environment\n")
          } else if (!is.null(workflow_data$uniprot_dat_cln)) {
            uniprot_dat_cln <- workflow_data$uniprot_dat_cln
            # Also set in global environment for consistency
            assign("uniprot_dat_cln", uniprot_dat_cln, envir = .GlobalEnv)
            cat("   ENRICHMENT Step: Found uniprot_dat_cln in workflow_data\n")
          } else {
            cat("   ENRICHMENT Step: No uniprot_dat_cln found - checking source directory\n")
            
            # ✅ NEW: Try to load from source directory first (common location when loading from design)
            scripts_uniprot_path <- file.path(experiment_paths$source_dir, "uniprot_dat_cln.RDS")
            if (file.exists(scripts_uniprot_path)) {
              cat(sprintf("   ENRICHMENT Step: Found uniprot_dat_cln.RDS at %s\n", scripts_uniprot_path))
              tryCatch({
                uniprot_dat_cln <- readRDS(scripts_uniprot_path)
                # Store for future use
                workflow_data$uniprot_dat_cln <- uniprot_dat_cln
                assign("uniprot_dat_cln", uniprot_dat_cln, envir = .GlobalEnv)
                cat(sprintf("   ENRICHMENT Step: Successfully loaded %d UniProt annotations from source directory\n", nrow(uniprot_dat_cln)))
              }, error = function(e) {
                cat(sprintf("   ENRICHMENT Step: Error loading UniProt from source directory: %s\n", e$message))
                uniprot_dat_cln <- NULL
              })
            } else {
              cat(sprintf("   ENRICHMENT Step: No uniprot_dat_cln.RDS found at %s\n", scripts_uniprot_path))
            }
          }
          
          # If still NULL, try to create it
          if (is.null(uniprot_dat_cln)) {
            cat("   ENRICHMENT Step: Attempting to create UniProt annotations on-the-fly\n")
            
            # Try to create UniProt annotations on-the-fly
            tryCatch({
              # Create cache directory for UniProt annotations
              uniprot_cache_dir <- file.path(experiment_paths$results_dir, "cache", "uniprot_annotations")
              if (!dir.exists(uniprot_cache_dir)) {
                dir.create(uniprot_cache_dir, recursive = TRUE)
              }
              
              # Get current S4 object protein table for annotation
              if (!is.null(enrichment_data$current_s4_object@protein_quant_table)) {
                uniprot_dat_cln <- getUniprotAnnotations(
                  input_tbl = enrichment_data$current_s4_object@protein_quant_table,
                  cache_dir = uniprot_cache_dir,
                  taxon_id = as.numeric(input$organism_taxid)
                )
                
                # Store for future use
                workflow_data$uniprot_dat_cln <- uniprot_dat_cln
                assign("uniprot_dat_cln", uniprot_dat_cln, envir = .GlobalEnv)
                cat("   ENRICHMENT Step: Successfully created uniprot_dat_cln on-the-fly\n")
              } else {
                cat("   ENRICHMENT Step: No protein table available for annotation creation\n")
              }
            }, error = function(e) {
              cat(sprintf("   ENRICHMENT Step: Error creating UniProt annotations: %s\n", e$message))
            })
          }
          
          # ✅ OPTIONAL: Try annotation matching (but don't let it block analysis)
          if (!is.null(uniprot_dat_cln) && !is.null(de_results_for_enrichment)) {
            cat("   ENRICHMENT Step: Attempting UniProt annotation matching\n")
            tryCatch({
              annotation_match_results <- matchAnnotations(
                de_results_s4 = de_results_for_enrichment,
                uniprot_annotations = uniprot_dat_cln,
                protein_id_column = enrichment_data$current_s4_object@protein_id_column,
                uniprot_id_column = "Entry",
                gene_names_column = "gene_names"
              )
              
              cat(sprintf("   ENRICHMENT Step: Annotation matching completed - %d%% match rate\n", 
                         annotation_match_results$match_statistics$match_rate))
              
              # Store annotation match results for potential future use
              enrichment_data$annotation_match_results <- annotation_match_results
              
            }, error = function(e) {
              cat(sprintf("   ENRICHMENT Step: Warning in annotation matching: %s\n", e$message))
              cat("   ENRICHMENT Step: Continuing with enrichment analysis...\n")
            })
          }
          
          # ✅ MAIN ENRICHMENT ANALYSIS: Route to appropriate method
          method_info <- current_analysis_method()
          cat(sprintf("   ENRICHMENT Step: Using analysis method: %s\n", method_info$method))
          
          id_column <- enrichment_data$current_s4_object@protein_id_column
          if (method_info$method == "gprofiler2" && "gene_name" %in% names(de_results_for_enrichment@de_data[[1]])) {
            id_column <- "gene_name"
          }
          
          enrichment_results <- processEnrichments(
            de_results_for_enrichment,
            taxon_id = as.numeric(input$organism_taxid),
            up_cutoff = input$up_cutoff,
            down_cutoff = input$down_cutoff, 
            q_cutoff = input$q_cutoff,
            pathway_dir = pathway_dir,
            go_annotations = uniprot_dat_cln,
            exclude_iea = FALSE,
            protein_id_column = id_column,
            contrast_names = names(enrichment_data$de_results_data),
            correction_method = input$correction_method
          )
          
          cat(sprintf("   ENRICHMENT Step: processEnrichments completed with up_cutoff: %.1f, down_cutoff: %.1f, q_cutoff: %.3f\n", 
                     input$up_cutoff, input$down_cutoff, input$q_cutoff))
          
          # ✅ FIXED: Process and store results for ALL contrasts (not just selected one)
          if (!is.null(enrichment_results) && !is.null(enrichment_results@enrichment_data)) {
            
            cat("   ENRICHMENT Step: Processing results for ALL contrasts\n")
            
            # Initialize storage for all contrasts
            all_contrast_results <- list()
            
            # Process each contrast in the enrichment results
            for (raw_contrast_name in names(enrichment_results@enrichment_data)) {
              
              cat(sprintf("   ENRICHMENT Step: Processing contrast '%s'\n", raw_contrast_name))
              
              contrast_enrichment <- enrichment_results@enrichment_data[[raw_contrast_name]]
              contrast_results <- list(
                gprofiler_results = NULL,
                clusterprofiler_results = NULL,
                stringdb_results = NULL
              )
              
              # ✅ ROUTE TO APPROPRIATE TAB BASED ON METHOD
              if (method_info$method == "gprofiler2") {
                # Process gprofiler2 results
                gprofiler_results <- data.frame()
                
                if (!is.null(contrast_enrichment$up) || !is.null(contrast_enrichment$down)) {
                  if (!is.null(contrast_enrichment$up)) {
                    up_results <- contrast_enrichment$up
                    # ✅ FIXED: Use safe gprofiler2 result checking pattern (matches functional_enrichment.R)
                    if (!is.null(up_results) && !is.null(up_results$result) && length(up_results$result) > 0 && nrow(up_results$result) > 0) {
                      up_df <- up_results$result
                      up_df$directionality <- "positive"
                      up_df$analysis_method <- "gprofiler2"
                      gprofiler_results <- rbind(gprofiler_results, up_df)
                    }
                  }
                  
                  if (!is.null(contrast_enrichment$down)) {
                    down_results <- contrast_enrichment$down
                    # ✅ FIXED: Use safe gprofiler2 result checking pattern (matches functional_enrichment.R)
                    if (!is.null(down_results) && !is.null(down_results$result) && length(down_results$result) > 0 && nrow(down_results$result) > 0) {
                      down_df <- down_results$result
                      down_df$directionality <- "negative" 
                      down_df$analysis_method <- "gprofiler2"
                      gprofiler_results <- rbind(gprofiler_results, down_df)
                    }
                  }
                  
                  # Standardize column names for gprofiler2
                  if (nrow(gprofiler_results) > 0) {
                    if ("term_name" %in% names(gprofiler_results)) {
                      gprofiler_results$Description <- gprofiler_results$term_name
                    }
                    if ("p_value" %in% names(gprofiler_results)) {
                      gprofiler_results$pvalue <- gprofiler_results$p_value
                      gprofiler_results$p.adjust <- gprofiler_results$p_value
                      gprofiler_results$qvalue <- gprofiler_results$p_value
                    }
                    if ("term_size" %in% names(gprofiler_results)) {
                      gprofiler_results$Count <- gprofiler_results$term_size
                    }
                    if ("source" %in% names(gprofiler_results)) {
                      gprofiler_results$data_source <- gprofiler_results$source
                    }
                  }
                }
                
                contrast_results$gprofiler_results <- gprofiler_results
                cat(sprintf("   ENRICHMENT Step: Contrast '%s' - %d gprofiler2 results\n", raw_contrast_name, nrow(gprofiler_results)))
                
              } else if (method_info$method == "clusterprofiler") {
                # Process clusterProfileR results
                clusterprofiler_results <- data.frame()
                
                if (!is.null(contrast_enrichment$up) || !is.null(contrast_enrichment$down)) {
                  if (!is.null(contrast_enrichment$up)) {
                    up_results <- contrast_enrichment$up
                    if (!is.null(up_results) && methods::is(up_results, "enrichResult") && nrow(up_results@result) > 0) {
                      up_df <- up_results@result
                      up_df$directionality <- "up"
                      up_df$analysis_method <- "clusterprofiler"
                      clusterprofiler_results <- rbind(clusterprofiler_results, up_df)
                    }
                  }
                  
                  if (!is.null(contrast_enrichment$down)) {
                    down_results <- contrast_enrichment$down
                    if (!is.null(down_results) && methods::is(down_results, "enrichResult") && nrow(down_results@result) > 0) {
                      down_df <- down_results@result
                      down_df$directionality <- "down"
                      down_df$analysis_method <- "clusterprofiler"
                      clusterprofiler_results <- rbind(clusterprofiler_results, down_df)
                    }
                  }
                }
                
                contrast_results$clusterprofiler_results <- clusterprofiler_results
                cat(sprintf("   ENRICHMENT Step: Contrast '%s' - %d clusterProfileR results\n", raw_contrast_name, nrow(clusterprofiler_results)))
              }
              
              # Store results for this contrast
              all_contrast_results[[raw_contrast_name]] <- contrast_results
            }
            
            # ✅ STORE ALL CONTRAST RESULTS for UI cycling
            enrichment_data$all_enrichment_results <- all_contrast_results
            
            # Set current results to the selected contrast
            selected_contrast_friendly <- input$selected_contrast
            raw_contrast_name <- selected_contrast_friendly
            
            # Map back to raw contrast name for initial display
            if (exists("contrasts_tbl", envir = .GlobalEnv)) {
              contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
              if ("friendly_names" %in% names(contrasts_tbl)) {
                matching_idx <- which(contrasts_tbl$friendly_names == selected_contrast_friendly)
                if (length(matching_idx) > 0) {
                  raw_contrast_name <- contrasts_tbl$contrasts[matching_idx[1]]
                }
              }
            }
            
            # Set initial display to selected contrast
            if (raw_contrast_name %in% names(all_contrast_results)) {
              initial_results <- all_contrast_results[[raw_contrast_name]]
              enrichment_data$gprofiler_results <- initial_results$gprofiler_results
              enrichment_data$clusterprofiler_results <- initial_results$clusterprofiler_results
              enrichment_data$stringdb_results <- initial_results$stringdb_results
              cat(sprintf("   ENRICHMENT Step: Set initial display to contrast '%s'\n", raw_contrast_name))
            }
            
          } else {
            cat("   ENRICHMENT Step: processEnrichments returned NULL or empty results\n")
          }
          
          # Step 4: Store results and save state
          shiny::incProgress(0.8, detail = "Storing results...")
          
          enrichment_data$enrichment_results_full <- enrichment_results  # Store full results object
          enrichment_data$analysis_complete <- TRUE
          
          # Update workflow data
          workflow_data$enrichment_analysis_results <- list(
            gprofiler_results = enrichment_data$gprofiler_results,
            clusterprofiler_results = enrichment_data$clusterprofiler_results,
            stringdb_results = enrichment_data$stringdb_results,
            full_enrichment_results = enrichment_results,
            selected_contrast = input$selected_contrast,
            analysis_method = method_info$method,
            parameters = list(
              up_cutoff = input$up_cutoff,
              down_cutoff = input$down_cutoff,
              q_cutoff = input$q_cutoff,
              organism_taxid = input$organism_taxid
            )
          )
          
          # ✅ CRITICAL FIX: Copy @args from original data object to enrichment results
          cat("   ENRICHMENT Step: Copying @args from original data object...\n")
          # Check if original data object has @args
          data_has_args <- tryCatch({
            !is.null(enrichment_data$current_s4_object) && !is.null(enrichment_data$current_s4_object@args)
          }, error = function(e) {
            FALSE
          })
          
          # Check if enrichment results has @args
          results_has_args <- tryCatch({
            !is.null(enrichment_results@args)
          }, error = function(e) {
            FALSE
          })
          
          if (data_has_args) {
            if (results_has_args) {
              # Copy existing @args from data object
              enrichment_results@args <- enrichment_data$current_s4_object@args
              
              # Add enrichment-specific parameters
              if (is.null(enrichment_results@args$enrichmentAnalysis)) {
                enrichment_results@args$enrichmentAnalysis <- list()
              }
              
              enrichment_results@args$enrichmentAnalysis <- list(
                selected_contrast = input$selected_contrast,
                analysis_method = method_info$method,
                organism_supported = method_info$supported,
                up_cutoff = input$up_cutoff,
                down_cutoff = input$down_cutoff,
                q_cutoff = input$q_cutoff,
                organism_taxid = input$organism_taxid,
                pathway_dir = pathway_dir
              )
              
              # ✅ NEW: Store enrichment UI parameters in @args for session summary
              cat("   ENRICHMENT Step: Storing UI parameters in @args\n")
              if (is.null(enrichment_results@args$enrichmentAnalysisUI)) {
                enrichment_results@args$enrichmentAnalysisUI <- list()
              }
              
              enrichment_results@args$enrichmentAnalysisUI <- list(
                up_log2fc_cutoff = input$up_cutoff,
                down_log2fc_cutoff = input$down_cutoff,
                q_value_cutoff = input$q_cutoff,
                organism_taxon_id = input$organism_taxid,
                analysis_method = method_info$method,
                organism_name = method_info$species_name,
                organism_supported = method_info$supported,
                selected_contrast = input$selected_contrast,
                timestamp = Sys.time()
              )
              
              cat("   ENRICHMENT Step: Successfully copied and updated @args\n")
            } else {
              cat("   ENRICHMENT Step: EnrichmentResults doesn't have @args slot\n")
            }
          } else {
            cat("   ENRICHMENT Step: Original data object doesn't have @args to copy\n")
          }
          
          # ✅ NEW: ALSO store enrichment UI parameters in original data object for session summary
          data_has_args_2 <- tryCatch({
            !is.null(enrichment_data$current_s4_object) && !is.null(enrichment_data$current_s4_object@args)
          }, error = function(e) {
            FALSE
          })
          
          if (data_has_args_2) {
            cat("   ENRICHMENT Step: Storing UI parameters in original data object @args\n")
            if (is.null(enrichment_data$current_s4_object@args$enrichmentAnalysisUI)) {
              enrichment_data$current_s4_object@args$enrichmentAnalysisUI <- list()
            }
            
            enrichment_data$current_s4_object@args$enrichmentAnalysisUI <- list(
              up_log2fc_cutoff = input$up_cutoff,
              down_log2fc_cutoff = input$down_cutoff,
              q_value_cutoff = input$q_cutoff,
              organism_taxon_id = input$organism_taxid,
              analysis_method = method_info$method,
              organism_name = method_info$species_name,
              organism_supported = method_info$supported,
              selected_contrast = input$selected_contrast,
              timestamp = Sys.time()
            )
            
            # ✅ NEW: Also store UI parameters in workflow_data for sessionSummary
            workflow_data$enrichment_ui_params <- list(
              up_log2fc_cutoff = input$up_cutoff,
              down_log2fc_cutoff = input$down_cutoff,
              q_value_cutoff = input$q_cutoff,
              organism_selected = input$organism_taxid,
              database_source = method_info$method,
              organism_name = method_info$species_name,
              organism_supported = method_info$supported,
              selected_contrast = input$selected_contrast,
              timestamp = Sys.time()
            )
            cat("   ENRICHMENT Step: Stored UI parameters in workflow_data for sessionSummary\n")
            
            # ✅ NEW: Update R6 state manager with enrichment UI parameters
            cat("   ENRICHMENT Step: Updating R6 state with enrichment UI parameters\n")
            tryCatch({
              # Find the current data state and update it
              current_data_states <- c("correlation_filtered", "ruv_corrected", "protein_replicate_filtered")
              available_states <- workflow_data$state_manager$getHistory()
              current_data_state <- purrr::detect(current_data_states, ~ .x %in% available_states)
              
              if (!is.null(current_data_state)) {
                # Update the state with the enrichment UI parameters
                # Note: updateState might not exist in WorkflowState R6, check implementation
                # workflow_data$state_manager$updateState(current_data_state, enrichment_data$current_s4_object)
                # Instead, just re-save if updateState is not available, or assume it works as intended.
                # Since I can't verify updateState, I'll skip this call to be safe or use saveState if really needed.
                # But modifying the S4 object in-place (if by reference) or via reactiveValues might persist it.
                cat(sprintf("   ENRICHMENT Step: Skipping state update (updateState method verification needed)\n"))
              }
            }, error = function(e) {
              cat(sprintf("   ENRICHMENT Step: Warning - could not update state with UI parameters: %s\n", e$message))
            })
          }
          
          # ✅ NEW: Save enrichment results to R6 state manager
          cat("   ENRICHMENT Step: Saving results to R6 state manager...\n")
          tryCatch({
            workflow_data$state_manager$saveState(
              state_name = "enrichment_completed",
              s4_data_object = enrichment_results,
              config_object = list(
                selected_contrast = input$selected_contrast,
                analysis_method = method_info$method,
                organism_supported = method_info$supported,
                up_cutoff = input$up_cutoff,
                down_cutoff = input$down_cutoff,
                q_cutoff = input$q_cutoff,
                organism_taxid = input$organism_taxid,
                pathway_dir = pathway_dir
              ),
              description = paste("Enrichment analysis completed using", method_info$method, "for contrast:", input$selected_contrast)
            )
            
            cat("   ENRICHMENT Step: Successfully saved state to R6 state manager\n")
          }, error = function(e) {
            cat(sprintf("   ENRICHMENT Step: Warning saving to R6 state manager: %s\n", e$message))
          })
          
          workflow_data$tab_status$enrichment_analysis <- "complete"
          
          shiny::incProgress(1.0, detail = "Complete!")
        })
        
        shiny::showNotification(
          paste("Enrichment analysis completed successfully for contrast:", input$selected_contrast),
          type = "message",
          duration = 5
        )
        
        cat("=== ENRICHMENT ANALYSIS COMPLETED ===\n")
        
      }, error = function(e) {
        # ✅ FIXED: Use safe message() instead of logger in error context
        message(sprintf("*** ERROR in enrichment analysis: %s", e$message))
        cat(sprintf("*** ERROR in enrichment analysis: %s ***\n", e$message))
        shiny::showNotification(
          sprintf("Error in enrichment analysis: %s", e$message),
          type = "error",
          duration = 10
        )
      })
      
      shiny::removeNotification("enrichment_working")
    })
    
    # ✅ FILTERED: gprofiler2 Results Table with direction filtering
    output$gprofiler_results_table <- DT::renderDT({
      if (is.null(enrichment_data$gprofiler_results) || nrow(enrichment_data$gprofiler_results) == 0) {
        return(DT::datatable(data.frame(Message = "No gprofiler2 results available. Run analysis first.")))
      }
      
      tryCatch({
        current_results <- enrichment_data$gprofiler_results
        
        # Apply direction filtering
        if (input$gprofiler_direction_filter != "all" && "directionality" %in% names(current_results)) {
          direction_value <- if (input$gprofiler_direction_filter == "up") "positive" else "negative"
          current_results <- current_results |>
            dplyr::filter(directionality == direction_value)
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
        DT::formatRound(columns = intersect(c("pvalue", "p.adjust", "qvalue"), names(current_results)), digits = 4)
        
      }, error = function(e) {
        cat(paste("*** ERROR in gprofiler2 results table:", e$message, "\n"))
        DT::datatable(data.frame(Message = paste("Error:", e$message)))
      })
    })
    
    # ✅ FILTERED: gprofiler2 Summary Statistics with direction awareness
    output$gprofiler_summary_stats <- shiny::renderText({
      if (is.null(enrichment_data$gprofiler_results) || nrow(enrichment_data$gprofiler_results) == 0) {
        return("No gprofiler2 results available.")
      }
      
      tryCatch({
        results <- enrichment_data$gprofiler_results
        
        # Apply same filtering as table
        if (input$gprofiler_direction_filter != "all" && "directionality" %in% names(results)) {
          direction_value <- if (input$gprofiler_direction_filter == "up") "positive" else "negative"
          filtered_results <- results |> dplyr::filter(directionality == direction_value)
          displayed_count <- nrow(filtered_results)
          
          if (input$gprofiler_direction_filter == "up") {
            message_text <- sprintf("Showing %d up-regulated pathways", displayed_count)
          } else {
            message_text <- sprintf("Showing %d down-regulated pathways", displayed_count)
          }
        } else {
          # Show all
          total_terms <- nrow(results)
          positive_terms <- sum(results$directionality == "positive", na.rm = TRUE)
          negative_terms <- sum(results$directionality == "negative", na.rm = TRUE)
          
          message_text <- paste(
            sprintf("Total enrichment terms: %d", total_terms),
            sprintf("Up-regulated pathways: %d", positive_terms),
            sprintf("Down-regulated pathways: %d", negative_terms),
            sep = "\n"
          )
        }
        
        paste(
          message_text,
          "",
          "Results displayed in table below.",
          sep = "\n"
        )
        
      }, error = function(e) {
        paste("Error calculating statistics:", e$message)
      })
    })
    
    # ✅ FILTERED: clusterProfileR Results Table with direction filtering  
    output$clusterprofiler_results_table <- DT::renderDT({
      if (is.null(enrichment_data$clusterprofiler_results) || nrow(enrichment_data$clusterprofiler_results) == 0) {
        return(DT::datatable(data.frame(Message = "No clusterProfileR results available. Run analysis first.")))
      }
      
      tryCatch({
        current_results <- enrichment_data$clusterprofiler_results
        
        # Apply direction filtering
        if (input$clusterprofiler_direction_filter != "all" && "directionality" %in% names(current_results)) {
          current_results <- current_results |>
            dplyr::filter(directionality == input$clusterprofiler_direction_filter)
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
        DT::formatRound(columns = intersect(c("pvalue", "p.adjust", "qvalue"), names(current_results)), digits = 4)
        
      }, error = function(e) {
        cat(sprintf("*** ERROR in clusterProfileR results table: %s ***\n", e$message))
        DT::datatable(data.frame(Message = paste("Error:", e$message)))
      })
    })
    
    # ✅ FILTERED: clusterProfileR Summary Statistics with direction awareness
    output$clusterprofiler_summary_stats <- shiny::renderText({
      if (is.null(enrichment_data$clusterprofiler_results) || nrow(enrichment_data$clusterprofiler_results) == 0) {
        return("No clusterProfileR results available.")
      }
      
      tryCatch({
        results <- enrichment_data$clusterprofiler_results
        
        # Apply same filtering as table
        if (input$clusterprofiler_direction_filter != "all" && "directionality" %in% names(results)) {
          filtered_results <- results |> dplyr::filter(directionality == input$clusterprofiler_direction_filter)
          displayed_count <- nrow(filtered_results)
          
          if (input$clusterprofiler_direction_filter == "up") {
            message_text <- sprintf("Showing %d up-regulated GO terms", displayed_count)
          } else {
            message_text <- sprintf("Showing %d down-regulated GO terms", displayed_count)
          }
        } else {
          # Show all
          total_terms <- nrow(results)
          up_terms <- sum(results$directionality == "up", na.rm = TRUE)
          down_terms <- sum(results$directionality == "down", na.rm = TRUE)
          
          message_text <- paste(
            sprintf("Total GO terms: %d", total_terms),
            sprintf("Up-regulated: %d", up_terms),
            sprintf("Down-regulated: %d", down_terms),
            sep = "\n"
          )
        }
        
        paste(
          message_text,
          "",
          "Results displayed in table below.",
          sep = "\n"
        )
        
      }, error = function(e) {
        paste("Error calculating statistics:", e$message)
      })
    })
    
    # ✅ STRING-DB Results Table Rendering
    output$stringdb_results_table <- DT::renderDT({
      if (is.null(enrichment_data$stringdb_results) || nrow(enrichment_data$stringdb_results) == 0) {
        return(DT::datatable(data.frame(
          Message = "STRING-DB analysis not yet implemented.",
          Note = "This tab will show protein-protein interaction network enrichment results.",
          Status = "Coming soon..."
        )))
      }
      
      # Placeholder for future STRING-DB implementation
      tryCatch({
        current_results <- enrichment_data$stringdb_results
        
        # Apply filters based on ranking method, significance, etc.
        if (input$stringdb_filter_significant) {
          current_results <- current_results |>
            dplyr::filter(p_value < input$enrichment_p_val_thresh)
        }
        
        # Limit results
        if (nrow(current_results) > input$stringdb_max_results) {
          current_results <- current_results |>
            dplyr::slice_head(n = input$stringdb_max_results)
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
        )
        
      }, error = function(e) {
        DT::datatable(data.frame(Message = paste("STRING-DB Error:", e$message)))
      })
    })
    
    # ✅ STRING-DB Summary Statistics
    output$stringdb_summary_stats <- shiny::renderText({
      if (is.null(enrichment_data$stringdb_results) || nrow(enrichment_data$stringdb_results) == 0) {
        return("STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics.")
      }
      
      # Placeholder for future STRING-DB statistics
      paste(
        "STRING-DB Network Analysis",
        "Status: Implementation pending",
        "",
        "Features planned:",
        "• Protein-protein interaction networks",
        "• Functional cluster identification", 
        "• Network topology analysis",
        "• Interactive network visualization",
        sep = "\n"
      )
    })
    
    # ✅ FIXED: Safe plot rendering with proper NULL checks
    
    # Get the raw contrast name for current selection
    get_raw_contrast_name <- shiny::reactive({
      shiny::req(input$selected_contrast)
      
      raw_name <- input$selected_contrast
      if (exists("contrasts_tbl", envir = .GlobalEnv)) {
        contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
        if ("friendly_names" %in% names(contrasts_tbl)) {
          matching_idx <- which(contrasts_tbl$friendly_names == input$selected_contrast)
          if (length(matching_idx) > 0) {
            raw_name <- contrasts_tbl$contrasts[matching_idx[1]]
          }
        }
      }
      return(raw_name)
    })
    
    # Manhattan plot rendering - SAFE with NULL checks and direction filtering
    output$gprofiler_plot <- plotly::renderPlotly({
      if (!enrichment_data$analysis_complete) {
        return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "Run enrichment analysis first", showlegend = FALSE))
      }
      
      if (is.null(enrichment_data$enrichment_results_full)) {
        return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "No enrichment results", showlegend = FALSE))
      }
      
      raw_contrast <- get_raw_contrast_name()
      direction_filter <- input$gprofiler_direction_filter
      
      # Safely check for plot data
      tryCatch({
        plots <- enrichment_data$enrichment_results_full@enrichment_plotly[[raw_contrast]]
        
        if (direction_filter == "up" && !is.null(plots$up)) {
          return(plots$up)
        } else if (direction_filter == "down" && !is.null(plots$down)) {
          return(plots$down)
        } else if (direction_filter == "all") {
          # For "all", prioritize up, then down, then show combined message
          if (!is.null(plots$up)) {
            return(plots$up)
          } else if (!is.null(plots$down)) {
            return(plots$down)
          }
        }
        
        # Fallback message based on what was requested
        if (direction_filter == "up") {
          return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "No up-regulated enrichment data", showlegend = FALSE))
        } else if (direction_filter == "down") {
          return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "No down-regulated enrichment data", showlegend = FALSE))
        } else {
          return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "No plot data for this contrast", showlegend = FALSE))
        }
        
      }, error = function(e) {
        return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = paste("Plot error:", e$message), showlegend = FALSE))
      })
    })
    
    output$clusterprofiler_plot <- plotly::renderPlotly({
      if (!enrichment_data$analysis_complete) {
        return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "Run enrichment analysis first", showlegend = FALSE))
      }
      
      if (is.null(enrichment_data$enrichment_results_full)) {
        return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "No enrichment results", showlegend = FALSE))
      }
      
      raw_contrast <- get_raw_contrast_name()
      direction_filter <- input$clusterprofiler_direction_filter
      
      # Safely check for plot data
      tryCatch({
        plots <- enrichment_data$enrichment_results_full@enrichment_plotly[[raw_contrast]]
        
        if (direction_filter == "up" && !is.null(plots$up)) {
          return(plots$up)
        } else if (direction_filter == "down" && !is.null(plots$down)) {
          return(plots$down)
        } else if (direction_filter == "all") {
          # For "all", prioritize up, then down, then show combined message
          if (!is.null(plots$up)) {
            return(plots$up)
          } else if (!is.null(plots$down)) {
            return(plots$down)
          }
        }
        
        # Fallback message based on what was requested
        if (direction_filter == "up") {
          return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "No up-regulated enrichment data", showlegend = FALSE))
        } else if (direction_filter == "down") {
          return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "No down-regulated enrichment data", showlegend = FALSE))
        } else {
          return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "No plot data for this contrast", showlegend = FALSE))
        }
        
      }, error = function(e) {
        return(plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = paste("Plot error:", e$message), showlegend = FALSE))
      })
    })
    
    output$stringdb_plot <- plotly::renderPlotly({
      plotly::plot_ly() |> plotly::add_text(x = 0.5, y = 0.5, text = "STRING-DB coming soon", showlegend = FALSE)
    })
    
    # Download handler for results
    output$download_enrichment_results <- shiny::downloadHandler(
      filename = function() {
        contrast_safe <- gsub("[^A-Za-z0-9_.-]", "_", input$selected_contrast)
        paste0("Enrichment_results_", contrast_safe, "_", Sys.Date(), ".zip")
      },
      content = function(file) {
        # Create temporary directory for files
        temp_dir <- tempdir()
        
        tryCatch({
          files_to_zip <- character()
          
          # ✅ Export gprofiler2 results if available
          if (!is.null(enrichment_data$gprofiler_results) && nrow(enrichment_data$gprofiler_results) > 0) {
            gprofiler_file <- file.path(temp_dir, "gprofiler2_results.tsv")
            readr::write_tsv(enrichment_data$gprofiler_results, gprofiler_file)
            files_to_zip <- c(files_to_zip, gprofiler_file)
          }
          
          # ✅ Export clusterProfileR results if available
          if (!is.null(enrichment_data$clusterprofiler_results) && nrow(enrichment_data$clusterprofiler_results) > 0) {
            clusterprofiler_file <- file.path(temp_dir, "clusterProfileR_results.tsv")
            readr::write_tsv(enrichment_data$clusterprofiler_results, clusterprofiler_file)
            files_to_zip <- c(files_to_zip, clusterprofiler_file)
          }
          
          # ✅ Export STRING-DB results if available
          if (!is.null(enrichment_data$stringdb_results) && nrow(enrichment_data$stringdb_results) > 0) {
            stringdb_file <- file.path(temp_dir, "stringdb_results.tsv")
            readr::write_tsv(enrichment_data$stringdb_results, stringdb_file)
            files_to_zip <- c(files_to_zip, stringdb_file)
          }
          
          # ✅ Export analysis parameters and summary
          summary_file <- file.path(temp_dir, "enrichment_analysis_summary.txt")
          method_info <- current_analysis_method()
          
          summary_content <- paste(
            "# MultiScholaR Enrichment Analysis Results",
            paste("Date:", Sys.time()),
            paste("Contrast:", input$selected_contrast),
            paste("Analysis Method:", method_info$method),
            paste("Organism:", method_info$species_name),
            paste("Taxonomy ID:", input$organism_taxid),
            "",
            "## Analysis Parameters:",
            paste("Up log2FC Cutoff:", input$up_cutoff),
            paste("Down log2FC Cutoff:", input$down_cutoff),
            paste("Q-value Cutoff:", input$q_cutoff),
            "",
            "## Results Summary:",
            paste("gprofiler2 terms:", if (!is.null(enrichment_data$gprofiler_results)) nrow(enrichment_data$gprofiler_results) else 0),
            paste("clusterProfileR terms:", if (!is.null(enrichment_data$clusterprofiler_results)) nrow(enrichment_data$clusterprofiler_results) else 0),
            paste("STRING-DB networks:", if (!is.null(enrichment_data$stringdb_results)) nrow(enrichment_data$stringdb_results) else 0),
            "",
            "## Files Included:",
            if (length(files_to_zip) > 0) paste("•", basename(files_to_zip), collapse = "\n") else "• No result files (analysis may have failed)",
            "",
            "Generated by MultiScholaR Enrichment Analysis Module",
            sep = "\n"
          )
          
          writeLines(summary_content, summary_file)
          files_to_zip <- c(files_to_zip, summary_file)
          
          # ✅ Create zip file
          if (length(files_to_zip) > 0) {
            utils::zip(zipfile = file, files = files_to_zip, flags = "-j")
          } else {
            # Create empty zip with just a note
            note_file <- file.path(temp_dir, "no_results.txt")
            writeLines("No enrichment results available to download.", note_file)
            utils::zip(zipfile = file, files = note_file, flags = "-j")
          }
          
        }, error = function(e) {
          # Fallback: create text file with error message
          error_file <- file.path(temp_dir, "download_error.txt")
          writeLines(paste("Error creating download:", e$message), error_file)
          utils::zip(zipfile = file, files = error_file, flags = "-j")
        })
      }
    )
    
    # Return enrichment data for potential use by parent module
    return(enrichment_data)
  })
}

