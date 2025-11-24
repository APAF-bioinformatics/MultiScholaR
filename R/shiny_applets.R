#' @title Run an applet
#' @param applet_type The type of applet to run (e.g., "designMatrix").
#' @param omic_type The omics context for the applet (e.g., "proteomics", "metabolomics").
#' @param experiment_label The label used when setting up directories (e.g., "workshop_data"). This is used to find the correct paths within the project_dirs_object.
#' @param project_dirs_object_name The name of the list object in the parent frame that holds the directory structures (typically the output of setupDirectories). Defaults to "project_dirs".
#' @param force Logical; if TRUE, skips user confirmation for overwriting existing files (specific to some applets).
#' @export
RunApplet <- function(applet_type, omic_type, experiment_label, project_dirs_object_name = "project_dirs", force = FALSE) {
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
    stop(paste0("Error: Project directories object ", sQuote(project_dirs_object_name), " not found in the calling environment. Please run setupDirectories() first."))
  }
  project_dirs <- get(project_dirs_object_name, envir = parent.frame())
  if (!is.list(project_dirs)) {
    stop(paste0("Error: ", sQuote(project_dirs_object_name), " is not a list as expected."))
  }

  # The key in project_dirs combines omic_type and experiment_label
  list_key <- paste(omic_type, experiment_label, sep = "_")
  if (!list_key %in% names(project_dirs)) {
    stop(paste0("Error: No directory information found for ", sQuote(list_key), " in ", sQuote(project_dirs_object_name), ". Ensure omic_type and experiment_label are correct."))
  }
  
  # Get the paths for the current omic type
  current_omic_paths <- project_dirs[[list_key]]
  if (!is.list(current_omic_paths) || !"source_dir" %in% names(current_omic_paths)) {
      stop(paste0("Error: 'source_dir' not found for ", sQuote(list_key), ". Check the structure of ", sQuote(project_dirs_object_name), "."))
  }
  
  # The 'source_dir' is where config files should be written.
  derived_source_dir <- current_omic_paths$source_dir

  # Validate the derived_source_dir
  if (!is.character(derived_source_dir) || length(derived_source_dir) != 1 || !dir.exists(derived_source_dir)) {
    stop(
      paste0("Error: Derived 'source_dir' (", derived_source_dir, ") for ", sQuote(list_key), " is not a single character string representing an existing directory. ",
             "Please check the output of setupDirectories() and ensure the path is correct.")
    )
  }
  
  # For simplicity within the applet logic, assign it to a variable named source_dir.
  source_dir <- derived_source_dir

  # --- Omic Type Specific Logic ---
  switch(omic_type,
    "proteomics" = {
      # Load the required module using a reliable method
      module_path <- system.file("shiny", "modules", "proteomics", "designMatrixBuilderModule.R", package = "MultiScholaR")

      if (!file.exists(module_path)) {
        stop("Could not find designMatrixBuilderModule.R. Please ensure the package is properly installed and the file exists in 'inst/shiny/modules/proteomics/'.")
      }
      
      source(module_path)
      
      if (applet_type == "designMatrix") {
        log_info(paste("Launching Design Matrix applet for PROTEOMICS. Using source_dir:", source_dir))

        # Check for data_tbl and config_list in parent frame
        if (!exists("data_tbl", envir = parent.frame())) {
          stop("Proteomics 'data_tbl' not found in the current environment for designMatrix applet.")
        }
        if (!exists("config_list", envir = parent.frame())) {
          stop("Proteomics 'config_list' not found in the current environment for designMatrix applet.")
        }

        # Get data_tbl and config_list from parent frame
        data_tbl <- get("data_tbl", envir = parent.frame())
        config_list <- get("config_list", envir = parent.frame())

        # Check for existing files
        existing_files <- list.files(source_dir, pattern = "^(design_matrix|data_cln).*\\.tab$")
        if (length(existing_files) > 0 && !force) {
          message("Found existing design/cleaned data files in ", source_dir, " for this proteomics experiment:")
          message(paste(" -", existing_files, collapse = "\n"))
          message("\nTo overwrite, relaunch with force = TRUE.")
          return(invisible(NULL))
        }

        # --- Shiny App Definition (using the new module) ---

        ui <- fluidPage(
          title = "Design Matrix Builder Applet",
          designMatrixBuilderUI("builder")
        )

        server <- function(input, output, session) {
          
          module_results_rv <- designMatrixBuilderServer(
            "builder",
            data_tbl = reactive(data_tbl),
            config_list = reactive(config_list)
          )

          observeEvent(module_results_rv(), {
            results <- module_results_rv()
            req(results)

            log_info("Standalone Applet: Saving results and closing.")
            
            showModal(shiny::modalDialog(
              title = "Saving Data",
              p("Saving results and closing applet..."),
              footer = NULL, easyClose = FALSE
            ))

            # Assign results back to the parent environment
            assign("design_matrix", results$design_matrix, envir = parent.frame())
            assign("data_cln", results$data_cln, envir = parent.frame())
            assign("config_list", results$config_list, envir = parent.frame())
            if (!is.null(results$contrasts_tbl)) {
              assign("contrasts_tbl", results$contrasts_tbl, envir = parent.frame())
            }

            # Also assign results to global environment (for R Markdown compatibility)
            assign("design_matrix", results$design_matrix, envir = .GlobalEnv)
            assign("data_cln", results$data_cln, envir = .GlobalEnv)
            assign("config_list", results$config_list, envir = .GlobalEnv)
            if (!is.null(results$contrasts_tbl)) {
              assign("contrasts_tbl", results$contrasts_tbl, envir = .GlobalEnv)
            }
            log_info("Results assigned to both parent frame and global environment.")

            # Write files to disk
            tryCatch({
              design_matrix_path <- file.path(source_dir, "design_matrix.tab")
              write.table(results$design_matrix, file = design_matrix_path, sep = "\t", row.names = FALSE, quote = FALSE)

              data_cln_path <- file.path(source_dir, "data_cln.tab")
              write.table(results$data_cln, file = data_cln_path, sep = "\t", row.names = FALSE, quote = FALSE)

              if (!is.null(results$contrasts_tbl) && nrow(results$contrasts_tbl) > 0) {
                contrast_path <- file.path(source_dir, "contrast_strings.tab")
                writeLines(results$contrasts_tbl$contrasts, contrast_path)
              }
            }, error = function(e) {
                log_error(paste("Error saving files in standalone applet:", e$message))
            })
            
            removeModal()
            stopApp(results)
          })
        }

        # Run the app
        result <- shiny::runApp(
          shiny::shinyApp(ui, server),
          launch.browser = TRUE
        )

        return(invisible(result))
      } else {
        stop(paste("Applet type", applet_type, "not implemented for proteomics yet."))
      }
    },
    { # Default for the omic_type switch
      log_error(paste("Unsupported omic_type provided to RunApplet:", omic_type))
      stop(paste("Unsupported omic_type:", omic_type))
    }
  )
  invisible(NULL)
}
