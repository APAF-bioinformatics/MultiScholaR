# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

localNamespaceBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localNamespaceBinding(
      env = env,
      name = name,
      value = bindings[[name]],
      .local_envir = .local_envir
    )
  }
}

renderSharedOutput <- function(value) {
  paste(capture.output(print(value)), collapse = "\n")
}

test_that("app server preserves setup flow, dynamic tabs, and watchdog module wiring", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$modals <- list()
  captured$notifications <- list()
  captured$setup_calls <- list()
  captured$server_calls <- list()
  captured$update_tab_calls <- list()

  base_dir <- tempfile("app-shell-base-")
  dir.create(base_dir, recursive = TRUE)
  withr::defer(unlink(base_dir, recursive = TRUE, force = TRUE))

  localNamespaceBindings(
    package_ns,
    list(
      setup_shiny_logger = function() {
        captured$logger_initialized <- TRUE
        invisible(NULL)
      },
      log_messages = function() "log line 1\nlog line 2",
      setupDirectories = function(base_dir, omic_types, label = NULL, force = FALSE, reuse_existing = FALSE) {
        captured$setup_calls[[length(captured$setup_calls) + 1L]] <- list(
          base_dir = base_dir,
          omic_types = omic_types,
          label = label,
          force = force,
          reuse_existing = reuse_existing
        )

        result <- list()
        for (omic in omic_types) {
          result[[omic]] <- list(
            base_dir = file.path(base_dir, omic),
            source_dir = file.path(base_dir, "scripts", omic),
            results_dir = file.path(base_dir, "results", omic)
          )
        }
        result
      },
      mod_proteomics_ui = function(id) shiny::div(class = "dynamic-ui", `data-id` = id, "proteomics ui"),
      mod_lipidomics_ui = function(id) shiny::div(class = "dynamic-ui", `data-id` = id, "lipidomics ui"),
      mod_proteomics_server = function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
        captured$server_calls[[length(captured$server_calls) + 1L]] <- list(
          id = id,
          project_dirs = project_dirs,
          omic_type = omic_type,
          experiment_label = experiment_label,
          volumes = volumes
        )
        list(kind = "proteomics", id = id)
      },
      mod_lipidomics_server = function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
        captured$server_calls[[length(captured$server_calls) + 1L]] <- list(
          id = id,
          project_dirs = project_dirs,
          omic_type = omic_type,
          experiment_label = experiment_label,
          volumes = volumes
        )
        list(kind = "lipidomics", id = id)
      }
    )
  )

  local_mocked_bindings(
    showModal = function(ui, ...) {
      captured$modals[[length(captured$modals) + 1L]] <- renderSharedOutput(ui)
      invisible(NULL)
    },
    removeModal = function(...) {
      captured$modal_removed <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = "default", duration = 5, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(
        ui = paste(ui, collapse = ""),
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    .package = "shiny"
  )

  local_mocked_bindings(
    getVolumes = function() {
      function() c(Home = base_dir)
    },
    shinyDirChoose = function(...) invisible(NULL),
    parseDirPath = function(...) base_dir,
    .package = "shinyFiles"
  )

  local_mocked_bindings(
    updateTabItems = function(session, inputId, selected, ...) {
      captured$update_tab_calls[[length(captured$update_tab_calls) + 1L]] <- list(
        inputId = inputId,
        selected = selected
      )
      invisible(NULL)
    },
    .package = "shinydashboard"
  )

  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    app_server,
    {
      session$flushReact()

      expect_true(isTRUE(captured$logger_initialized))
      expect_match(output$log_terminal_content, "log line 1", fixed = TRUE)

      session$setInputs(selected_omics = c("proteomics", "lipidomics"))
      session$setInputs(start_analysis = 1)
      session$flushReact()

      session$setInputs(experiment_label = "demo_project", project_dir = base_dir)
      session$setInputs(confirm_setup = 1)
      session$flushReact()

      expect_equal(captured$setup_calls[[1]]$omic_types, c("proteomics", "lipidomics"))
      expect_true(captured$setup_calls[[1]]$force)
      expect_false(captured$setup_calls[[1]]$reuse_existing)

      dynamic_menu <- renderSharedOutput(output$dynamic_menu)
      expect_match(dynamic_menu, "Proteomics", fixed = TRUE)
      expect_match(dynamic_menu, "Lipidomics", fixed = TRUE)

      dynamic_tabs <- renderSharedOutput(output$dynamic_tabs)
      expect_match(dynamic_tabs, "proteomics ui", fixed = TRUE)
      expect_match(dynamic_tabs, "lipidomics ui", fixed = TRUE)

      expect_equal(captured$update_tab_calls[[1]]$inputId, "main_menu")
      expect_equal(captured$update_tab_calls[[1]]$selected, "proteomics_tab")

      session$flushReact()
      expect_equal(length(captured$server_calls), 2L)
      expect_equal(captured$server_calls[[1]]$id, "proteomics_workflow")
      expect_equal(captured$server_calls[[1]]$omic_type, "proteomics")
      expect_equal(captured$server_calls[[2]]$id, "lipidomics_workflow")
      expect_equal(captured$server_calls[[2]]$omic_type, "lipidomics")

      expect_equal(sort(values$initialized_omics), c("lipidomics", "proteomics"))
      expect_equal(values$workflow_state$proteomics$kind, "proteomics")
      expect_equal(values$workflow_state$lipidomics$kind, "lipidomics")
      expect_equal(values$experiment_label, "demo_project")
      expect_true(!is.null(values$project_dirs$proteomics))
      expect_true(!is.null(values$project_dirs$lipidomics))
    }
  )
})
