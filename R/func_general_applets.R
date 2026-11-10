# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
    result <- runDesignMatrixAppletFrom(
      source_dir = source_dir,
      output_dir = output_dir,
      omic_type = omic_type,
      force = force,
      caller_env = parent.frame()
    )
  }
  
  invisible(result)
}


#' @title Run Design Matrix Builder Applet (internal)
#' @description Internal function that launches the proteomics design matrix builder as a
#'   standalone Shiny app for the notebook workbooks. The applet works on the caller's long
#'   peptide table (`data_tbl`, with its `Run` column) and `config_list`, as the workbooks'
#'   data-management chunk defines them, and holds the table in its own session owner over a
#'   temporary store, which is closed and removed when the applet ends (KIT-048: the builder has
#'   one path, the owner's). On save it writes `design_matrix.tab`, `data_cln.tab` and, when
#'   contrasts were defined, `contrast_strings.tab` to `source_dir`, and
#'   `design_matrix_results.rds` to `output_dir`, and assigns `design_matrix`, `data_cln`,
#'   `config_list` and `contrasts_tbl` in the caller's environment, where the workbooks' later
#'   chunks read them.
#' @param source_dir Path to the project's source directory, where the design files are written.
#' @param output_dir Path to output directory for saving results.
#' @param omic_type The omics type. Only `"proteomics"` has a standalone applet.
#' @param force Logical; if TRUE, skips confirmation for overwriting.
#' @return List containing design_matrix, data_cln, contrasts_tbl, config_list
#' @keywords internal
.runDesignMatrixApplet <- function(source_dir, output_dir, omic_type, force = FALSE) {
  runDesignMatrixAppletFrom(source_dir, output_dir, omic_type, force, caller_env = parent.frame())
}

#' The design matrix applet (`.runDesignMatrixApplet()`), run for a given calling environment:
#' `RunApplet()` passes its own caller's.
#' @param source_dir,output_dir,omic_type,force As `.runDesignMatrixApplet()`.
#' @param caller_env The environment holding `data_tbl` and `config_list`, which receives the
#'   saved design.
#' @return As `.runDesignMatrixApplet()`.
#' @noRd
runDesignMatrixAppletFrom <- function(source_dir, output_dir, omic_type, force, caller_env) {
  if (!identical(omic_type, "proteomics")) {
    cli::cli_abort(
      c(
        "the design matrix applet has no standalone builder for {.val {omic_type}}.",
        i = "It serves proteomics only; build other designs in the app's design step."
      ),
      class = "design_applet_omic_unsupported"
    )
  }
  inputs <- readDesignAppletInputs(caller_env)

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

  # The applet's own session store, removed with its owner when the applet ends.
  store_root <- tempfile("design-applet-store-")
  dir.create(store_root)
  applet <- new.env(parent = emptyenv())
  on.exit(closeDesignAppletStore(applet$owner, store_root), add = TRUE)

  # Run the app (blocking)
  message("Launching Design Matrix Builder...")
  message("Design files will be written to: ", source_dir)
  message("Output will be saved to: ", output_file)

  result <- shiny::runApp(
    shiny::shinyApp(
      ui = mod_prot_design_builder_ui("designBuilder"),
      server = designAppletServer(inputs, source_dir, output_file, store_root, applet)
    ),
    launch.browser = TRUE
  )

  # The workbooks' later chunks read these from the environment the applet was run in.
  for (name in names(result)) {
    if (!is.null(result[[name]])) assign(name, result[[name]], envir = caller_env)
  }

  return(result)
}

#' The standalone design applet's inputs, as the workbooks define them before the design chunk:
#' the long peptide table `data_tbl` (a `Run` column names each row's run) and `config_list`.
#' @param caller_env The environment the applet was called from.
#' @return A list: `data_tbl`, `config_list` and the builder's `column_mapping` (the run column
#'   and, when present, DIA-NN's precursor quantity column, as the DIA-NN import maps them).
#' @noRd
readDesignAppletInputs <- function(caller_env) {
  data_tbl <- get0("data_tbl", envir = caller_env)
  if (!is.data.frame(data_tbl)) {
    cli::cli_abort(
      c(
        "the design matrix applet found no {.code data_tbl} table where it was called.",
        i = "Read the search report into {.code data_tbl} first (the workbook's data chunk)."
      ),
      class = "design_applet_input_missing"
    )
  }
  if (!"Run" %in% names(data_tbl)) {
    cli::cli_abort(
      "the design matrix applet needs a {.code Run} column in {.code data_tbl}.",
      class = "design_applet_input_invalid"
    )
  }
  config_list <- get0("config_list", envir = caller_env)
  if (!is.list(config_list)) {
    cli::cli_abort(
      c(
        "the design matrix applet found no {.code config_list} where it was called.",
        i = "Read the configuration first (the workbook's {.code readConfigFile()} call)."
      ),
      class = "design_applet_input_missing"
    )
  }
  quantity_col <- intersect(c("Precursor.Normalised", "Precursor.Quantity"), names(data_tbl))
  column_mapping <- list(run_col = "Run")
  if (length(quantity_col) > 0) column_mapping$quantity_col <- quantity_col[[1]]
  list(data_tbl = data_tbl, config_list = config_list, column_mapping = column_mapping)
}

#' The standalone design applet's server. It creates the session's owner over the applet's
#' temporary store as the app does (`createDIASessionOwner`, bound to the DIA route as a DIA-NN
#' import binds it), hands the table to the owner as the import does, and runs the design tab's
#' builder on it. A saved design is written out and ends the app.
#' @param inputs The applet's inputs (`readDesignAppletInputs()`).
#' @param source_dir,output_file Where the design files and the result file are written.
#' @param store_root The directory holding the applet's temporary store.
#' @param applet An environment; `owner` is set to the session's owner, for the applet's cleanup.
#' @return A Shiny server function.
#' @noRd
designAppletServer <- function(inputs, source_dir, output_file, store_root, applet) {
  function(input, output, session) {
    workflow_data <- shiny::reactiveValues(config_list = inputs$config_list)
    store_paths <- list(source_dir = store_root)
    shiny::isolate({
      workflow_data$state_manager <- createDIASessionOwner(workflow_data, store_paths, session)
      applet$owner <- initializeDIANativeImport(
        workflow_data, store_paths, session,
        format = "diann", data_type = "peptide"
      )
      withDIAStateScope(applet$owner, function() {
        applyProtImportResultToWorkflow(
          workflowData = workflow_data,
          dataImportResult = list(
            data = inputs$data_tbl,
            data_type = "peptide",
            column_mapping = inputs$column_mapping
          ),
          format = "diann",
          fastaPath = NULL
        )
      })
    })

    saved_design <- registerProtDesignBuilderModule(workflow_data, moduleId = "designBuilder")

    # Write the saved design out and stop the app with it
    shiny::observeEvent(getDIAReactiveReference(saved_design), {
      withDIAStateScope(applet$owner, function() {
        saved <- saved_design()
        shiny::req(saved)
        shiny::stopApp(writeDesignAppletResult(saved, applet$owner, source_dir, output_file))
      })
    }, ignoreNULL = TRUE)
  }
}

#' Write a saved design as the workbooks read it and return it with its table in memory.
#' @param saved The builder's saved result (`design_matrix`, `contrasts_tbl`, `config_list` and
#'   the cleaned table's `data_reference`).
#' @param owner The applet's session owner.
#' @param source_dir,output_file Where the design files and the result file are written.
#' @return A list: `design_matrix`, `data_cln`, `contrasts_tbl` and `config_list`.
#' @noRd
writeDesignAppletResult <- function(saved, owner, source_dir, output_file) {
  result <- list(
    design_matrix = saved$design_matrix,
    data_cln = resolveDIAOwnedValue(saved$data_reference, owner),
    contrasts_tbl = saved$contrasts_tbl,
    config_list = saved$config_list
  )
  utils::write.table(
    result$design_matrix,
    file = file.path(source_dir, "design_matrix.tab"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  writeDIADesignTable(owner, saved$data_reference, file.path(source_dir, "data_cln.tab"))
  if (!is.null(result$contrasts_tbl) && nrow(result$contrasts_tbl) > 0) {
    writeLines(result$contrasts_tbl$contrasts, file.path(source_dir, "contrast_strings.tab"))
  }
  saveRDS(result, output_file)
  message("Design matrix saved to: ", source_dir)
  result
}

#' Close the standalone design applet's owner (the session's end has usually closed it already)
#' and remove its temporary store.
#' @param owner The applet's session owner, or NULL when no session started.
#' @param store_root The directory holding the applet's temporary store.
#' @return NULL, invisibly.
#' @noRd
closeDesignAppletStore <- function(owner, store_root) {
  tryCatch(closeDIAStateAdapter(owner), finally = unlink(store_root, recursive = TRUE))
  invisible(NULL)
}
