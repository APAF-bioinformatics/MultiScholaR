#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinyjs
#' @import shinydashboard
#' @import logger
#' @import shinyFiles
#' @noRd
app_server <- function(input, output, session) {
  # Initialize the logger for Shiny
  # Note: setup_shiny_logger should be exported or available in the package
  if (exists("setup_shiny_logger")) {
    setup_shiny_logger()
  } else {
    logger::log_info("Custom logger setup not found, using default")
  }

  # Initialize reactive values
  values <- reactiveValues(
    selected_omics = character(),
    workflow_state = list(),
    current_tab = "home",
    project_dirs = NULL,
    experiment_label = NULL,
    project_base_dir = NULL,
    initialized_omics = character() # Track which modules have been spun up
  )

  # Log Terminal UI Rendering
  output$log_terminal_ui <- shiny::renderUI({
    shinyjqui::jqui_resizable(
      shiny::div(
        id = "log-terminal-wrapper",
        shiny::div(
          id = "log-terminal-panel",
          shiny::div(
            id = "log-terminal-header",
            onclick = "Shiny.setInputValue('toggle_log', Math.random())",
            shiny::h5(shiny::icon("terminal"), "Log Output")
          ),
          shinyjs::hidden(
            shiny::div(
              id = "log-content-wrapper",
              shiny::verbatimTextOutput("log_terminal_content")
            )
          )
        )
      )
    )
  })

  # Render log content
  output$log_terminal_content <- shiny::renderText({
    if (exists("log_messages")) {
      log_messages()
    } else {
      "Log messages not available"
    }
  })

  # Toggle log visibility
  observeEvent(input$toggle_log, {
    shinyjs::toggle(id = "log-content-wrapper", anim = TRUE)
  })

  # Set up volumes for shinyFiles at the top level
  logger::log_info("--- MultiScholaR App Starting ---")
  volumes <- NULL
  if (requireNamespace("shinyFiles", quietly = TRUE)) {
    message("   app_server: shinyFiles package is available")
    volumes <- shinyFiles::getVolumes()() # Call the function to get the actual volumes list
    message(sprintf(
      "   app_server: Created volumes. Type = %s, class = %s",
      typeof(volumes), paste(class(volumes), collapse = ", ")
    ))

    # Inspect volumes
    tryCatch(
      {
        message(sprintf("   app_server: volumes length = %d", length(volumes)))
        if (length(volumes) > 0) {
          message(sprintf(
            "   app_server: Volume names: %s",
            paste(names(volumes), collapse = ", ")
          ))
        }
      },
      error = function(e) {
        message(sprintf("   app_server ERROR inspecting volumes: %s", e$message))
      }
    )
  } else {
    message("   app_server: shinyFiles package NOT available")
  }

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

    # Show modal to get experiment label and output directory
    showModal(modalDialog(
      title = "Experiment Setup",
      size = "m",
      shiny::fluidRow(
        shiny::column(
          12,
          shiny::textInput("experiment_label",
            "Experiment Label:",
            value = paste0("analysis_", format(Sys.Date(), "%Y%m%d")),
            placeholder = "e.g., workshop_data",
            width = "100%"
          ),
          shiny::br(),
          shiny::h4("Project Output Directory"),
          shiny::p("Select where to create the project directories:"),
          shiny::fluidRow(
            shiny::column(
              9,
              shiny::textInput("project_dir",
                "Directory Path:",
                value = normalizePath("~", winslash = "/"),
                placeholder = "Enter directory path",
                width = "100%"
              )
            ),
            shiny::column(
              3,
              if (requireNamespace("shinyFiles", quietly = TRUE)) {
                shinyFiles::shinyDirButton("browse_dir",
                  "Browse...",
                  "Select Project Output Directory",
                  icon = shiny::icon("folder-open"),
                  style = "width: 100%; margin-top: 25px;"
                )
              } else {
                shiny::actionButton("browse_dir_fallback",
                  "Browse...",
                  icon = shiny::icon("folder-open"),
                  width = "100%",
                  style = "margin-top: 25px;"
                )
              }
            )
          ),
          shiny::tags$small(
            shiny::tags$i(
              "A new folder will be created here with your experiment label"
            )
          )
        )
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_setup", "Continue", class = "btn-primary")
      )
    ))
  })


  # Set up shinyFiles directory browser if available
  if (requireNamespace("shinyFiles", quietly = TRUE) && !is.null(volumes)) {
    shinyFiles::shinyDirChoose(input, "browse_dir",
      roots = volumes,
      session = session
    )

    observeEvent(input$browse_dir, {
      if (!is.null(input$browse_dir) && !is.integer(input$browse_dir)) {
        path <- shinyFiles::parseDirPath(volumes, input$browse_dir)
        if (length(path) > 0) {
          updateTextInput(session, "project_dir", value = path)
        }
      }
    })
  } else {
    # Fallback for when shinyFiles is not available
    observeEvent(input$browse_dir_fallback, {
      showNotification(
        "Please install the 'shinyFiles' package for directory browsing, or enter the path manually.",
        type = "message",
        duration = 5
      )
    })
  }

  # Handle experiment setup confirmation
  observeEvent(input$confirm_setup, {
    req(input$experiment_label)
    req(input$project_dir)

    # Validate directory
    if (!dir.exists(input$project_dir)) {
      showNotification(
        paste("Base directory does not exist:", input$project_dir),
        type = "error",
        duration = NULL
      )
      return()
    }

    values$experiment_label <- input$experiment_label
    values$project_base_dir <- file.path(input$project_dir, input$experiment_label)

    # Check if the final experiment directory already exists
    if (dir.exists(values$project_base_dir)) {
      # Show modal asking user what to do
      showModal(modalDialog(
        title = "Project Already Exists",
        paste("A project named '", values$experiment_label, "' already exists at this location."),
        "What would you like to do?",
        footer = tagList(
          modalButton("Cancel"),
          actionButton("reuse_project", "Use Existing Project"),
          actionButton("overwrite_project", "Overwrite Project", class = "btn-danger")
        )
      ))
    } else {
      # Directory doesn't exist, so proceed with creation
      removeModal()
      # Using a reactiveVal to trigger the setup process
      setup_action("create_new")
    }
  })

  # Reactive value to control the setup process
  setup_action <- reactiveVal(NULL)

  # Observer for "Use Existing"
  observeEvent(input$reuse_project, {
    removeModal()
    setup_action("reuse_existing")
  })

  # Observer for "Overwrite"
  observeEvent(input$overwrite_project, {
    removeModal()
    setup_action("create_new") # Same action as creating new, but 'force' will handle it
  })

  # Centralized observer to run setupDirectories
  observeEvent(setup_action(), {
    action <- setup_action()
    req(action)

    experiment_dir <- values$project_base_dir

    tryCatch(
      {
        # Call setupDirectories with the appropriate arguments

        # For both "create_new" and "overwrite", we pre-stage scripts and use force=TRUE
        if (action == "create_new") {
          logger::log_info("Creating new project directories...")
          # Pre-stage script templates
          for (omic in values$selected_omics) {
            # Try installed package path first (inst/workbooks/{omic}/)
            source_template_dir <- system.file("workbooks", omic, package = "MultiScholaR")

            # Fallback if package not installed or path differs (dev mode)
            if (source_template_dir == "" || !dir.exists(source_template_dir)) {
              source_template_dir <- file.path("inst", "workbooks", omic)
            }

            dest_template_dir <- file.path(experiment_dir, "scripts", omic)

            if (dir.exists(source_template_dir)) {
              # Copy all files from the workbooks folder to scripts destination
              fs::dir_copy(source_template_dir, dest_template_dir, overwrite = TRUE)
              logger::log_info(paste("Copied workbook templates for", omic, "to", dest_template_dir))
            } else {
              logger::log_warn(paste("No workbook templates found for", omic, "at", source_template_dir))
            }
          }

          if (exists("setupDirectories")) {
            project_dirs_list <- setupDirectories(
              base_dir = experiment_dir,
              omic_types = values$selected_omics,
              label = NULL,
              force = TRUE # This will handle overwriting if the button was clicked
            )
          } else {
            # Fallback if function not found (dev mode issue?)
            logger::log_error("setupDirectories function not found")
            stop("setupDirectories function not found")
          }
        } else if (action == "reuse_existing") {
          logger::log_info("Reusing existing project directories...")
          if (exists("setupDirectories")) {
            project_dirs_list <- setupDirectories(
              base_dir = experiment_dir,
              omic_types = values$selected_omics,
              label = NULL,
              force = FALSE,
              reuse_existing = TRUE # Use the new parameter
            )
          } else {
            stop("setupDirectories function not found")
          }
        } else {
          return() # Do nothing if action is not recognized
        }

        values$project_dirs <- project_dirs_list

        logger::log_info(paste("Project setup complete. Action:", action))

        # Common logic to proceed after setup
        proceed_with_analysis()
      },
      error = function(e) {
        logger::log_error(paste("Failed to set up project directories:", e$message))
        showNotification(paste("Failed to set up project directories:", e$message), type = "error", duration = NULL)
      }
    )

    # Reset the action so it can be triggered again
    setup_action(NULL)
  })

  # Helper function to continue the app initialization after directory setup
  proceed_with_analysis <- function() {
    # Modules are auto-loaded by package/load_all, so no need to sourceDir

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
          tryCatch(
            {
              message(sprintf("   [%s] Starting UI module load", omic))

              # Call the UI function for this omics type
              # Refactored to use mod_{omic}_ui convention
              ui_function_name <- paste0("mod_", omic, "_ui")
              message(sprintf("   [%s] Looking for function: %s", omic, ui_function_name))

              if (exists(ui_function_name)) {
                message(sprintf("   [%s] Function exists check: TRUE", omic))
                ui_function <- get(ui_function_name)
                message(sprintf(
                  "   [%s] Function retrieved. Type: %s, Class: %s",
                  omic, typeof(ui_function), class(ui_function)
                ))

                # Check if it's actually a function
                if (!is.function(ui_function)) {
                  message(sprintf("   [%s] ERROR: Object is not a function!"))
                  stop(sprintf("%s is not a function", ui_function_name))
                }

                # Try to inspect the function
                message(sprintf("   [%s] Function formals:", omic))
                print(formals(ui_function))

                # Call the function
                # Using standard ID
                module_id <- paste0(omic, "_workflow")
                message(sprintf("   [%s] Calling function with id: %s", omic, module_id))

                result <- ui_function(module_id)

                message(sprintf(
                  "   [%s] Function call completed. Result type: %s, class: %s",
                  omic, typeof(result), class(result)
                ))

                return(result)
              } else {
                message(sprintf("   [%s] Function exists check: FALSE", omic))
                shiny::div(
                  class = "alert alert-warning",
                  shiny::h4("Module not yet implemented"),
                  shiny::p(paste("The", omic, "workflow module is under development."))
                )
              }
            },
            error = function(e) {
              message(sprintf("   [%s] ERROR in tryCatch: %s", omic, e$message))
              message(sprintf("   [%s] Error call stack:", omic))
              print(sys.calls())

              shiny::div(
                class = "alert alert-danger",
                shiny::h4("Error loading module"),
                shiny::p("Error details:", e$message)
              )
            }
          )
        )
      })

      message(sprintf("   Total tab items created: %d", length(tab_items)))
      result <- do.call(shiny::tagList, tab_items)
      message(sprintf("--- Exiting dynamic_tabs renderUI ---"))
      return(result)
    })

    # Switch to first omics tab if staying on "home"
    if (values$current_tab == "home" && length(values$selected_omics) > 0) {
      shinydashboard::updateTabItems(session, "main_menu", paste0(values$selected_omics[1], "_tab"))
    }
  }

  # =========================================================================
  # REACTIVE WATCHDOG: Server Module Initialization
  # =========================================================================
  # This observer watches for changes in selected_omics. If a new omic is
  # selected (that hasn't been initialized yet), it spins up the server module.
  # This ensures that even late-added omics get their server logic running.
  observe({
    # Only proceed if the project is set up
    req(values$project_dirs)
    req(values$selected_omics)

    # Identify which selected omics haven't been initialized
    new_omics <- setdiff(values$selected_omics, values$initialized_omics)

    if (length(new_omics) > 0) {
      logger::log_info(sprintf("Watchdog detected new omics to initialize: %s", paste(new_omics, collapse = ", ")))

      # --- 1. Perform Filesystem Setup for New Omics ---
      # We need to ensure folders exist and templates are copied for these new omics
      experiment_dir <- values$project_base_dir

      # Copy templates for new omics
      for (omic in new_omics) {
        # Try installed package path first (inst/workbooks/{omic}/)
        source_template_dir <- system.file("workbooks", omic, package = "MultiScholaR")

        # Fallback if package not installed or path differs (dev mode)
        if (source_template_dir == "" || !dir.exists(source_template_dir)) {
          source_template_dir <- file.path("inst", "workbooks", omic)
        }

        dest_template_dir <- file.path(experiment_dir, "scripts", omic)

        if (dir.exists(source_template_dir)) {
          # Copy all files from the workbooks folder to scripts destination
          # Ensure destination parent exists
          if (!dir.exists(dirname(dest_template_dir))) {
            dir.create(dirname(dest_template_dir), recursive = TRUE)
          }
          fs::dir_copy(source_template_dir, dest_template_dir, overwrite = TRUE)
          logger::log_info(paste("Copied workbook templates for", omic, "to", dest_template_dir))
        } else {
          logger::log_warn(paste("No workbook templates found for", omic, "at", source_template_dir))
        }
      }

      # Update directory structure (Additive)
      if (exists("setupDirectories")) {
        # We pass the full list of selected omics. setupDirectories should be robust enough
        # to just ensure folders exist. We use reuse_existing=TRUE and force=FALSE.
        new_project_dirs <- setupDirectories(
          base_dir = experiment_dir,
          omic_types = values$selected_omics,
          label = NULL,
          force = FALSE,
          reuse_existing = TRUE
        )
        values$project_dirs <- new_project_dirs
        logger::log_info("Updated project directories for new omics.")
      } else {
        logger::log_error("setupDirectories function not found during dynamic setup")
      }

      for (omic in new_omics) {
        message(sprintf("--- Initializing server module for: %s ---", omic))

        tryCatch(
          {
            # Call the server function for this omics type
            server_function_name <- paste0("mod_", omic, "_server")
            message(sprintf("   Looking for server function: %s", server_function_name))

            if (exists(server_function_name)) {
              server_function <- get(server_function_name)

              # Check if it's actually a function
              if (!is.function(server_function)) {
                stop(sprintf("%s is not a function", server_function_name))
              }

              message(sprintf("   Calling server function for %s", omic))

              # Initialize the Module
              workflow_data <- server_function(
                id = paste0(omic, "_workflow"),
                project_dirs = values$project_dirs,
                omic_type = omic,
                experiment_label = values$experiment_label,
                volumes = volumes
              )

              # Store the workflow data
              values$workflow_state[[omic]] <- workflow_data
              message(sprintf("   Workflow data stored for %s", omic))

              # Mark as initialized immediately to prevent re-entry
              values$initialized_omics <- c(values$initialized_omics, omic)
            } else {
              message(sprintf("   Server function exists: FALSE for %s", omic))
              # Mark as initialized to stop trying
              values$initialized_omics <- c(values$initialized_omics, omic)
            }
          },
          error = function(e) {
            message(sprintf("   ERROR initializing server for %s: %s", omic, e$message))
            logger::log_error(paste("Could not initialize server for", omic, ":", e$message))
            # Even on error, we might want to mark as initialized to prevent infinite loops,
            # or we might want to retry? For now, let's NOT mark it so it might try again if fixed?
            # Actually, better to mark it to avoid crashing the observer repeatedly.
            values$initialized_omics <- c(values$initialized_omics, omic)
          }
        )
      }
    }
  })
}
