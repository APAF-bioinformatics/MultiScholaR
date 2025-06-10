# MultiScholaR Interactive Workflow Application
# Main entry point for the Shiny app

# Load required packages with explicit namespace handling
library(shiny)
library(shinyjs)
library(shinydashboard)
library(logger)

# Handle namespace conflicts explicitly
library(DT, warn.conflicts = FALSE)

# Source helper functions to get access to loadDependencies
if (file.exists("../../helper_functions.R")) {
  source("../../helper_functions.R")
} else if (exists("loadDependencies", mode = "function")) {
  # Function already available (e.g., from package namespace)
} else {
  warning("loadDependencies function not found. Dependencies may not be properly loaded.")
}

# Source shiny_applets.R for RunApplet function
# Try multiple possible locations
shiny_applets_paths <- c(
  "../../shiny_applets.R",
  "../shiny_applets.R",
  "R/shiny_applets.R",
  file.path(dirname(getwd()), "shiny_applets.R")
)

shiny_applets_loaded <- FALSE
for (path in shiny_applets_paths) {
  if (file.exists(path)) {
    source(path)
    message(paste("Loaded shiny_applets.R from:", path))
    shiny_applets_loaded <- TRUE
    break
  }
}

if (!shiny_applets_loaded && !exists("RunApplet")) {
  warning("shiny_applets.R not found and RunApplet function not available.")
}

# Ensure all dependencies are loaded
tryCatch({
  if (exists("loadDependencies", mode = "function")) {
    loadDependencies(verbose = FALSE)
  }
}, error = function(e) {
  logger::log_warn("Failed to load some dependencies: {e$message}")
})

# Source all module files
sourceDir <- function(path, recursive = TRUE) {
  if (!dir.exists(path)) {
    logger::log_warn("Directory does not exist: {path}")
    return(invisible(NULL))
  }
  files <- list.files(path, pattern = "\\.R$", full.names = TRUE, recursive = recursive)
  lapply(files, source)
}

# Source common modules
sourceDir("modules/common")
sourceDir("ui/common")
sourceDir("server/common")

# Pre-load proteomics modules (as they're ready)
sourceDir("modules/proteomics")
sourceDir("ui/proteomics")
sourceDir("server/proteomics")

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "MultiScholaR Workflow"),
  
  dashboardSidebar(
    sidebarMenu(
      id = "main_menu",
      menuItem("Home", tabName = "home", icon = shiny::icon("home")),
      
      # Dynamic menu items will be added based on omics selection
      uiOutput("dynamic_menu")
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    
    # CSS for better styling - use shiny:: prefix
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .omic-selection-box {
          border: 2px solid #3c8dbc;
          border-radius: 10px;
          padding: 20px;
          margin: 10px;
          background-color: white;
          cursor: pointer;
          transition: all 0.3s;
        }
        .omic-selection-box:hover {
          transform: scale(1.05);
          box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }
        .omic-selection-box.selected {
          background-color: #d1ecf1;
          border-color: #bee5eb;
        }
      "))
    ),
    
    shinydashboard::tabItems(
      # Home tab with omics selection
      shinydashboard::tabItem(
        tabName = "home",
        fluidRow(
          column(12,
            h2("Welcome to MultiScholaR Interactive Workflow"),
            p("Select the type(s) of omics data you want to analyze:"),
            br()
          )
        ),
        
        fluidRow(
          column(4,
            shiny::div(
              id = "proteomics_box",
              class = "omic-selection-box",
              onclick = "toggleOmicSelection('proteomics')",
              h3(shiny::icon("dna"), "Proteomics"),
              p("Analyze protein abundance data from DIA-NN or similar tools")
            )
          ),
          column(4,
            shiny::div(
              id = "metabolomics_box",
              class = "omic-selection-box",
              onclick = "toggleOmicSelection('metabolomics')",
              h3(shiny::icon("flask"), "Metabolomics"),
              p("Process and analyze metabolite profiling data")
            )
          ),
          column(4,
            shiny::div(
              id = "transcriptomics_box",
              class = "omic-selection-box",
              onclick = "toggleOmicSelection('transcriptomics')",
              h3(shiny::icon("chart-line"), "Transcriptomics"),
              p("Analyze gene expression data from RNA-seq")
            )
          )
        ),
        
        fluidRow(
          column(4,
            shiny::div(
              id = "lipidomics_box",
              class = "omic-selection-box",
              onclick = "toggleOmicSelection('lipidomics')",
              h3(shiny::icon("oil-can"), "Lipidomics"),
              p("Comprehensive lipid profiling and analysis")
            )
          ),
          column(4,
            shiny::div(
              id = "integration_box",
              class = "omic-selection-box disabled",
              h3(shiny::icon("project-diagram"), "Integration"),
              p("Multi-omics integration (available when 2+ omics selected)")
            )
          )
        ),
        
        br(),
        
        fluidRow(
          column(12,
            actionButton(
              "start_analysis",
              "Start Analysis",
              class = "btn-primary btn-lg",
              width = "100%",
              icon = shiny::icon("play")
            )
          )
        ),
        
        # Hidden inputs to track selections
        shinyjs::hidden(
          shiny::div(
            id = "selection_tracker",
            checkboxInput("proteomics_selected", "Proteomics", FALSE),
            checkboxInput("metabolomics_selected", "Metabolomics", FALSE),
            checkboxInput("transcriptomics_selected", "Transcriptomics", FALSE),
            checkboxInput("lipidomics_selected", "Lipidomics", FALSE),
            checkboxInput("integration_selected", "Integration", FALSE)
          )
        )
      ),
      
      # Dynamic tab items will be added here
      uiOutput("dynamic_tabs")
    ),
    
    # JavaScript for handling omics selection - use shiny:: prefix
    shiny::tags$script(shiny::HTML("
      var selectedOmics = [];
      
      function toggleOmicSelection(omic) {
        if (omic === 'integration') return; // Integration is auto-enabled
        
        var box = document.getElementById(omic + '_box');
        var checkbox = document.getElementById(omic + '_selected');
        
        if (selectedOmics.includes(omic)) {
          selectedOmics = selectedOmics.filter(item => item !== omic);
          box.classList.remove('selected');
        } else {
          selectedOmics.push(omic);
          box.classList.add('selected');
        }
        
        // Update hidden checkboxes
        checkbox.checked = selectedOmics.includes(omic);
        checkbox.dispatchEvent(new Event('change'));
        
        // Enable/disable integration based on selection count
        var integrationBox = document.getElementById('integration_box');
        if (selectedOmics.length >= 2) {
          integrationBox.classList.remove('disabled');
          integrationBox.onclick = function() { toggleOmicSelection('integration'); };
        } else {
          integrationBox.classList.add('disabled');
          integrationBox.classList.remove('selected');
          integrationBox.onclick = null;
          selectedOmics = selectedOmics.filter(item => item !== 'integration');
          document.getElementById('integration_selected').checked = false;
          document.getElementById('integration_selected').dispatchEvent(new Event('change'));
        }
        
        Shiny.setInputValue('selected_omics', selectedOmics);
      }
    "))
  )
)

# Define server
server <- function(input, output, session) {
  # Initialize reactive values
  values <- reactiveValues(
    selected_omics = character(),
    workflow_state = list(),
    current_tab = "home",
    project_dirs = NULL,
    experiment_label = NULL
  )
  
  # Track omics selection
  observe({
    values$selected_omics <- input$selected_omics
  })
  
  # Handle start analysis button
  observeEvent(input$start_analysis, {
    if (length(values$selected_omics) == 0) {
      showNotification(
        "Please select at least one omics type to analyze",
        type = "warning"
      )
      return()
    }
    
    # Show modal to get experiment label
    showModal(modalDialog(
      title = "Experiment Setup",
      textInput("experiment_label", 
                "Experiment Label:", 
                value = paste0("analysis_", format(Sys.Date(), "%Y%m%d")),
                placeholder = "e.g., workshop_data"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_setup", "Continue", class = "btn-primary")
      )
    ))
  })
  
  # Handle experiment setup confirmation
  observeEvent(input$confirm_setup, {
    req(input$experiment_label)
    
    values$experiment_label <- input$experiment_label
    removeModal()
    
    # Log the selection
    logger::log_info("Starting analysis for: {paste(values$selected_omics, collapse = ', ')}")
    logger::log_info("Experiment label: {values$experiment_label}")
    
    # Initialize project directories
    tryCatch({
      # Don't actually create directories - just set up the structure needed for the workflow
      values$project_dirs <- list()
      
      for (omic in values$selected_omics) {
        key <- paste0(omic, "_", values$experiment_label)
        values$project_dirs[[key]] <- list(
          base_dir = getwd(),
          data_dir = file.path(tempdir(), "multischolar", values$experiment_label, "data", omic),
          results_dir = file.path(tempdir(), "multischolar", values$experiment_label, "results", omic),
          source_dir = file.path(tempdir(), "multischolar", values$experiment_label, "source", omic),
          # These are used by various functions but won't actually be created
          results_summary_dir = file.path(tempdir(), "multischolar", values$experiment_label, "results_summary", omic),
          qc_dir = file.path(tempdir(), "multischolar", values$experiment_label, "results", omic, "qc"),
          de_output_dir = file.path(tempdir(), "multischolar", values$experiment_label, "results", omic, "de_output")
        )
      }
      
      logger::log_info("Created temporary workflow structure with keys: {paste(names(values$project_dirs), collapse = ', ')}")
      
    }, error = function(e) {
      logger::log_error("Failed to create workflow structure: {e$message}")
      showNotification("Failed to create workflow structure", type = "error")
      return()
    })
    
    # Load omics-specific modules if not already loaded
    for (omic in values$selected_omics) {
      if (omic != "proteomics") { # proteomics already loaded
        tryCatch({
          sourceDir(paste0("modules/", omic))
          sourceDir(paste0("ui/", omic))
          sourceDir(paste0("server/", omic))
          logger::log_info("Loaded modules for {omic}")
        }, error = function(e) {
          logger::log_warn("Could not load modules for {omic}: {e$message}")
        })
      }
    }
    
    # Generate dynamic UI elements
    output$dynamic_menu <- renderUI({
      menu_items <- lapply(values$selected_omics, function(omic) {
        menuItem(
          text = tools::toTitleCase(omic),
          tabName = paste0(omic, "_tab"),
          icon = switch(omic,
            proteomics = shiny::icon("dna"),
            metabolomics = shiny::icon("flask"),
            transcriptomics = shiny::icon("chart-line"),
            lipidomics = shiny::icon("oil-can"),
            integration = shiny::icon("project-diagram")
          )
        )
      })
      do.call(tagList, menu_items)
    })
    
    # Generate dynamic tab content
    output$dynamic_tabs <- renderUI({
      message(sprintf("--- Entering dynamic_tabs renderUI ---"))
      message(sprintf("   Selected omics: %s", paste(values$selected_omics, collapse = ", ")))
      
      tab_items <- lapply(values$selected_omics, function(omic) {
        message(sprintf("   Processing omic: %s", omic))
        
        shinydashboard::tabItem(
          tabName = paste0(omic, "_tab"),
          shiny::h2(tools::toTitleCase(omic), " Workflow"),
          shiny::h4("Experiment: ", values$experiment_label),
          shiny::br(),
          
          # Load the appropriate UI module
          tryCatch({
            message(sprintf("   [%s] Starting UI module load", omic))
            
            # Call the UI function for this omics type
            ui_function_name <- paste0(omic, "WorkflowUi")
            message(sprintf("   [%s] Looking for function: %s", omic, ui_function_name))
            
            if (exists(ui_function_name)) {
              message(sprintf("   [%s] Function exists check: TRUE", omic))
              ui_function <- get(ui_function_name)
              message(sprintf("   [%s] Function retrieved. Type: %s, Class: %s", 
                            omic, typeof(ui_function), class(ui_function)))
              
              # Check if it's actually a function
              if (!is.function(ui_function)) {
                message(sprintf("   [%s] ERROR: Object is not a function!", omic))
                stop(sprintf("%s is not a function", ui_function_name))
              }
              
              # Try to inspect the function
              message(sprintf("   [%s] Function formals:", omic))
              print(formals(ui_function))
              
              # Call the function
              module_id <- paste0(omic, "_workflow")
              message(sprintf("   [%s] Calling function with id: %s", omic, module_id))
              
              result <- ui_function(module_id)
              
              message(sprintf("   [%s] Function call completed. Result type: %s, class: %s", 
                            omic, typeof(result), class(result)))
              message(sprintf("   [%s] Result structure:", omic))
              utils::str(result, max.level = 2)
              
              return(result)
            } else {
              message(sprintf("   [%s] Function exists check: FALSE", omic))
              shiny::div(
                class = "alert alert-warning",
                shiny::h4("Module not yet implemented"),
                shiny::p(paste("The", omic, "workflow module is under development."))
              )
            }
          }, error = function(e) {
            message(sprintf("   [%s] ERROR in tryCatch: %s", omic, e$message))
            message(sprintf("   [%s] Error call stack:", omic))
            print(sys.calls())
            
            shiny::div(
              class = "alert alert-danger",
              shiny::h4("Error loading module"),
              shiny::p("Error details:", e$message)
            )
          })
        )
      })
      
      message(sprintf("   Total tab items created: %d", length(tab_items)))
      result <- do.call(shiny::tagList, tab_items)
      message(sprintf("--- Exiting dynamic_tabs renderUI ---"))
      return(result)
    })
    
    # Initialize server modules for each selected omics
    for (omic in values$selected_omics) {
      message(sprintf("--- Initializing server module for: %s ---", omic))
      
      tryCatch({
        # Call the server function for this omics type
        server_function_name <- paste0(omic, "WorkflowServer")
        message(sprintf("   Looking for server function: %s", server_function_name))
        
        if (exists(server_function_name)) {
          message(sprintf("   Server function exists: TRUE"))
          server_function <- get(server_function_name)
          message(sprintf("   Server function type: %s, class: %s", 
                        typeof(server_function), class(server_function)))
          
          # Check if it's actually a function
          if (!is.function(server_function)) {
            message(sprintf("   ERROR: Server object is not a function!"))
            stop(sprintf("%s is not a function", server_function_name))
          }
          
          # Pass the correct parameters
          message(sprintf("   Calling server function with parameters:"))
          message(sprintf("     id: %s_workflow", omic))
          message(sprintf("     project_dirs: %s", typeof(values$project_dirs)))
          message(sprintf("     omic_type: %s", omic))
          message(sprintf("     experiment_label: %s", values$experiment_label))
          
          workflow_data <- server_function(
            id = paste0(omic, "_workflow"),
            project_dirs = values$project_dirs,
            omic_type = omic,
            experiment_label = values$experiment_label
          )
          
          message(sprintf("   Server function returned. Type: %s", typeof(workflow_data)))
          
          # Store the workflow data
          values$workflow_state[[omic]] <- workflow_data
          message(sprintf("   Workflow data stored for %s", omic))
        } else {
          message(sprintf("   Server function exists: FALSE"))
        }
      }, error = function(e) {
        message(sprintf("   ERROR initializing server for %s: %s", omic, e$message))
        logger::log_error("Could not initialize server for {omic}: {e$message}")
      })
    }
    message("--- Server module initialization complete ---")
    
    # Switch to first omics tab
    shinydashboard::updateTabItems(session, "main_menu", paste0(values$selected_omics[1], "_tab"))
  })
}

# Run the app
shinyApp(ui = ui, server = server) 