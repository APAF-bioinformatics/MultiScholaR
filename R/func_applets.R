# ----------------------------------------------------------------------------
# Applet Launcher Functions
# ----------------------------------------------------------------------------
# These functions launch standalone Shiny applets from R Markdown workflows.
# They are separate from the main MultiScholaR GUI (run_app).
# ----------------------------------------------------------------------------


# ============================================================================
# Workflow Stepper UI Component
# ============================================================================

#' @title Render Workflow Stepper UI
#' @description Creates a modern, professional workflow progress stepper component
#'              that displays step completion status with visual indicators.
#' @param steps List of step definitions. Each step should have:
#'   \itemize{
#'     \item \code{name}: Display name for the step
#'     \item \code{key}: Key matching the tab_status list entry
#'     \item \code{icon}: FontAwesome icon name (without "fa-" prefix)
#'   }
#' @param tab_status Named list of step statuses. Values should be one of:
#'   "complete", "pending", or "disabled".
#' @return A \code{shiny::tags$div} containing the stepper HTML
#' @examples
#' \dontrun{
#' steps <- list(
#'   list(name = "Import", key = "setup_import", icon = "file-import"),
#'   list(name = "Design", key = "design_matrix", icon = "th")
#' )
#' tab_status <- list(setup_import = "complete", design_matrix = "pending")
#' render_workflow_stepper(steps, tab_status)
#' }
#' @export
#' @importFrom shiny tags icon
render_workflow_stepper <- function(steps, tab_status) {
    # Determine status for each step
    step_states <- vapply(steps, function(step) {
        status <- tab_status[[step$key]]
        if (is.null(status) || status == "disabled") {
            return("pending")
        } else if (status == "complete") {
            return("completed")
        } else {
            return("current")
        }
    }, character(1))
    
    # Find the first non-completed step to mark as current
    # Only one step should be "current" - the first pending/active one
    found_current <- FALSE
    for (i in seq_along(step_states)) {
        if (step_states[i] == "current" && !found_current) {
            found_current <- TRUE
        } else if (step_states[i] == "current" && found_current) {
            step_states[i] <- "pending"
        } else if (step_states[i] == "pending" && !found_current) {
            step_states[i] <- "current"
            found_current <- TRUE
        }
    }
    
    # If no current step found (all complete), leave as is
    # If no steps are pending/current, all are complete
    
    # Build stepper elements
    step_elements <- lapply(seq_along(steps), function(i) {
        step <- steps[[i]]
        state <- step_states[i]
        
        # Determine step class based on state
        step_class <- paste0("stepper-step step-", state)
        circle_class <- paste0("stepper-circle ", state)
        
        # Circle content: checkmark for completed, icon for current/pending
        circle_content <- if (state == "completed") {
            shiny::icon("check")
        } else {
            shiny::icon(step$icon)
        }
        
        shiny::tags$div(
            class = step_class
            , shiny::tags$div(
                class = circle_class
                , circle_content
            )
            , shiny::tags$div(
                class = "stepper-label"
                , step$name
            )
        )
    })
    
    # Wrap in container
    shiny::tags$div(
        class = "workflow-stepper"
        , step_elements
    )
}

#' @title Run an applet
#' @description Launch a standalone Shiny applet from R Markdown workflows.
#' Currently supports the design matrix builder applet.
#' @param applet_type The type of applet to run (e.g., "designMatrix").
#' @param omic_type The omics context for the applet (e.g., "proteomics", "metabolomics").
#' @param experiment_label The label used when setting up directories (e.g., "workshop_data"). 
#'   This is used to find the correct paths within the project_dirs_object.
#' @param project_dirs_object_name The name of the list object in the parent frame that holds 
#'   the directory structures (typically the output of setupDirectories). Defaults to "project_dirs".
#' @param force Logical; if TRUE, skips user confirmation for overwriting existing files.
#' @return Invisibly returns the result from the applet (e.g., design matrix, contrasts).
#' @export
RunApplet <- function(applet_type, omic_type, experiment_label, 
                      project_dirs_object_name = "project_dirs", force = FALSE) {
  # Load required packages
  require(shiny)
  require(DT)
  require(gtools) # For mixedsort
  require(dplyr)  # For data manipulation
  require(logger) # For logging

  # --- Input Validation ---
  valid_applet_types <- c("designMatrix") # Add more as they are developed
  if (!applet_type %in% valid_applet_types) {
    stop(
      "Invalid applet_type specified. Valid options are: ",
      paste(valid_applet_types, collapse = ", ")
    )
  }

  valid_omic_types <- c("proteomics", "metabolomics", "lipidomics", "transcriptomics", "integration")
  if (!omic_type %in% valid_omic_types) {
    stop(
      "Invalid omic_type specified. Valid options are: ",
      paste(valid_omic_types, collapse = ", ")
    )
  }

  if (!is.character(experiment_label) || length(experiment_label) != 1 || experiment_label == "") {
    stop("'experiment_label' must be a single non-empty character string (e.g., 'workshop_data').")
  }
  
  if (!is.character(project_dirs_object_name) || length(project_dirs_object_name) != 1 || project_dirs_object_name == "") {
    stop("'project_dirs_object_name' must be a single non-empty character string.")
  }

  # --- Derive source_dir from project_dirs object and experiment_label ---
  if (!exists(project_dirs_object_name, envir = parent.frame())) {
    stop(paste0("Error: Project directories object ", sQuote(project_dirs_object_name), 
                " not found in the calling environment. Please run setupDirectories() first."))
  }
  project_dirs <- get(project_dirs_object_name, envir = parent.frame())
  if (!is.list(project_dirs)) {
    stop(paste0("Error: ", sQuote(project_dirs_object_name), " is not a list as expected."))
  }

  # The key in project_dirs combines omic_type and experiment_label
  list_key <- paste(omic_type, experiment_label, sep = "_")
  if (!list_key %in% names(project_dirs)) {
    stop(paste0("Error: No directory information found for ", sQuote(list_key), " in ", 
                sQuote(project_dirs_object_name), ". Ensure omic_type and experiment_label are correct."))
  }
  
  # Get the paths for the current omic type
  current_omic_paths <- project_dirs[[list_key]]
  if (!is.list(current_omic_paths) || !"source_dir" %in% names(current_omic_paths)) {
    stop(paste0("Error: The entry for ", sQuote(list_key), " does not contain 'source_dir'."))
  }
  
  source_dir <- current_omic_paths$source_dir
  output_dir <- current_omic_paths$output_dir
  
  # --- Run the appropriate applet ---
  if (applet_type == "designMatrix") {
    result <- .runDesignMatrixApplet(
      source_dir = source_dir,
      output_dir = output_dir,
      omic_type = omic_type,
      force = force
    )
  }
  
  invisible(result)
}


#' @title Run Design Matrix Builder Applet (internal)
#' @description Internal function that launches the design matrix builder as a standalone Shiny app.
#' @param source_dir Path to source directory containing sample files.
#' @param output_dir Path to output directory for saving results.
#' @param omic_type The omics type (proteomics, metabolomics, etc.)
#' @param force Logical; if TRUE, skips confirmation for overwriting.
#' @return List containing design_matrix, data_cln, contrasts_tbl, config_list
#' @keywords internal
.runDesignMatrixApplet <- function(source_dir, output_dir, omic_type, force = FALSE) {
  
  # Check for existing output
  output_file <- file.path(output_dir, "design_matrix_results.rds")
  if (file.exists(output_file) && !force) {
    response <- readline(prompt = paste0(
      "Design matrix results already exist at:\n  ", output_file, 
      "\nOverwrite? (y/n): "
    ))
    if (!tolower(response) %in% c("y", "yes")) {
      message("Applet cancelled. Existing results preserved.")
      return(invisible(NULL))
    }
  }
  
  # Read sample files from source directory
  sample_files <- list.files(source_dir, pattern = "\\.(tsv|csv|txt)$", full.names = TRUE)
  if (length(sample_files) == 0) {
    stop(paste0("No sample files (tsv/csv/txt) found in: ", source_dir))
  }
  
  # Read the first file to get column structure (assuming all files have same structure)
  # This is used to build the initial data table
  first_file <- sample_files[1]
  if (grepl("\\.tsv$", first_file)) {
    sample_data <- utils::read.delim(first_file, stringsAsFactors = FALSE, nrows = 5)
  } else {
    sample_data <- utils::read.csv(first_file, stringsAsFactors = FALSE, nrows = 5)
  }
  
  # Build initial data table with file info
  data_tbl <- data.frame(
    File = basename(sample_files),
    FullPath = sample_files,
    SampleName = gsub("\\.(tsv|csv|txt)$", "", basename(sample_files)),
    stringsAsFactors = FALSE
  )
  
  # Default config for standalone mode
  config_list <- list(
    omic_type = omic_type,
    source_dir = source_dir,
    output_dir = output_dir,
    formula = "~ 0 + Condition"
  )
  

  # Default column mapping for standalone mode
  column_mapping <- list(
    sample_id_col = "SampleName",
    run_col = "File",
    intensity_col = NULL
  )
  
  # Create reactive wrappers for module
  data_tbl_reactive <- shiny::reactiveVal(data_tbl)
  config_list_reactive <- shiny::reactiveVal(config_list)
  column_mapping_reactive <- shiny::reactiveVal(column_mapping)
  
  # Result storage
  applet_result <- NULL
  
  # Build and run the Shiny app
  ui <- mod_prot_design_builder_ui("designBuilder")
  
  server <- function(input, output, session) {
    # Call the module server
    result <- mod_prot_design_builder_server(
      id = "designBuilder",
      data_tbl = data_tbl_reactive,
      config_list = config_list_reactive,
      column_mapping = column_mapping_reactive
    )
    
    # Observe the result and stop app when saved
    shiny::observe({
      res <- result()
      if (!is.null(res)) {
        # Save to file
        saveRDS(res, output_file)
        message("Design matrix saved to: ", output_file)
        
        # Store result and stop app
        applet_result <<- res
        shiny::stopApp(res)
      }
    })
  }
  
  # Run the app (blocking)
  message("Launching Design Matrix Builder...")
  message("Source directory: ", source_dir)
  message("Output will be saved to: ", output_file)
  
  result <- shiny::runApp(
    shiny::shinyApp(ui = ui, server = server),
    launch.browser = TRUE
  )
  
  return(result)
}

