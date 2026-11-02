# Extracted from test-metab-01q-qc-intensity-module-characterization.R:1400

# prequel ----------------------------------------------------------------------
library(testthat)
source("helpers-scoped-mocked-bindings.R")
register_binding_teardown(
  asNamespace("shiny"),
  c(
    "moduleServer",
    "reactiveVal",
    "observeEvent",
    "renderText",
    "renderUI",
    "renderPlot",
    "req",
    "showNotification",
    "removeNotification"
  ),
  .local_envir = environment()
)
register_binding_teardown(
  asNamespace("logger"),
  c("log_info", "log_error", "log_warn"),
  .local_envir = environment()
)
register_binding_teardown(
  asNamespace("MultiScholaR"),
  c("metaboliteIntensityFiltering", "updateMetaboliteFiltering"),
  .local_envir = environment()
)
repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")

      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}
loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "mod_metab_qc_intensity_ui_helpers.R"),
    file.path(repo_root, "R", "mod_metab_qc_intensity.R")
  ),
  symbols = c(
    "buildMetabIntensityAssayTabsUi",
    "renderMetabIntensityFilterPlot",
    "buildMetabIntensityFilterStats",
    "updateMetabIntensityFilterQcPlot",
    "buildMetabIntensityFilterSummary",
    "reportMetabIntensityFilterSuccess",
    "reportMetabIntensityFilterRevertSuccess",
    "mod_metab_qc_intensity_server"
  ),
  env = environment()
)
if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id"
    )
  )
}
.metab_intensity_test_env <- environment()
makeMetabIntensityAssayData <- function(metabolite_data, metabolite_id_column = "metabolite_id") {
  sample_columns <- unique(unlist(lapply(
    metabolite_data,
    function(assay) setdiff(names(assay), metabolite_id_column)
  )))
  sample_columns <- sort(sample_columns)

  args <- list(
    Class = "MetaboliteAssayData",
    metabolite_data = metabolite_data,
    metabolite_id_column = metabolite_id_column
  )

  slots <- methods::slotNames("MetaboliteAssayData")
  if ("annotation_id_column" %in% slots) {
    args$annotation_id_column <- "annotation_id"
  }
  if ("database_identifier_type" %in% slots) {
    args$database_identifier_type <- "test_id"
  }
  if ("internal_standard_regex" %in% slots) {
    args$internal_standard_regex <- NA_character_
  }
  if ("design_matrix" %in% slots) {
    args$design_matrix <- data.frame(
      Sample_ID = sample_columns,
      group = "A",
      stringsAsFactors = FALSE
    )
  }
  if ("sample_id" %in% slots) {
    args$sample_id <- "Sample_ID"
  }
  if ("group_id" %in% slots) {
    args$group_id <- "group"
  }
  if ("technical_replicate_id" %in% slots) {
    args$technical_replicate_id <- NA_character_
  }
  if ("args" %in% slots) {
    args$args <- list()
  }

  do.call(methods::new, args)
}
metabIntensityHasExtractedSeams <- function() {
  all(vapply(
    c(
      "buildMetabIntensityAssayTabsUi",
      "renderMetabIntensityFilterPlot",
      "buildMetabIntensityFilterStats",
      "updateMetabIntensityFilterQcPlot",
      "buildMetabIntensityFilterSummary",
      "reportMetabIntensityFilterSuccess",
      "reportMetabIntensityFilterRevertSuccess"
    ),
    exists,
    logical(1),
    envir = .metab_intensity_test_env,
    inherits = FALSE
  ))
}

# test -------------------------------------------------------------------------
skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )
captured <- new.env(parent = emptyenv())
captured$notifications <- list()
captured$reactive_values <- list()
server_env <- environment(mod_metab_qc_intensity_server)
binding_names <- c(
    "buildMetabIntensityFilterStats",
    "updateMetabIntensityFilterQcPlot",
    "reportMetabIntensityFilterSuccess",
    "metaboliteIntensityFiltering"
  )
had_bindings <- stats::setNames(
    vapply(
      binding_names,
      function(name) exists(name, envir = server_env, inherits = FALSE),
      logical(1)
    ),
    binding_names
  )
old_bindings <- vector("list", length(binding_names))
names(old_bindings) <- binding_names
for (name in binding_names) {
    if (had_bindings[[name]]) {
      old_bindings[[name]] <- get(name, envir = server_env, inherits = FALSE)
    }
  }
on.exit({
    for (name in binding_names) {
      if (had_bindings[[name]]) {
        assign(name, old_bindings[[name]], envir = server_env)
      } else if (exists(name, envir = server_env, inherits = FALSE)) {
        rm(list = name, envir = server_env)
      }
    }
  }, add = TRUE)
current_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c("m1", "m2"),
        Sample1 = c(10, 20),
        Sample2 = c(5, 25),
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        metabolite_id = c("m3"),
        Sample1 = 7,
        Sample2 = 9,
        stringsAsFactors = FALSE
      )
    )
  )
filtered_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "m2",
        Sample1 = 20,
        Sample2 = 25,
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        metabolite_id = "m3",
        Sample1 = 7,
        Sample2 = 9,
        stringsAsFactors = FALSE
      )
    )
  )
expected_stats <- data.frame(
    Assay = c("Plasma", "Serum"),
    Original = c(2, 1),
    Filtered = c(1, 1),
    Removed = c(1, 0),
    Percent_Retained = c(50, 100),
    stringsAsFactors = FALSE
  )
assign(
    "metaboliteIntensityFiltering",
    function(
      theObject,
      metabolites_intensity_cutoff_percentile,
      metabolites_proportion_of_samples_below_cutoff
    ) {
      captured$filter_call <- list(
        the_object = theObject,
        intensity_cutoff_percentile = metabolites_intensity_cutoff_percentile,
        proportion_below_cutoff = metabolites_proportion_of_samples_below_cutoff
      )
      filtered_s4
    },
    envir = server_env
  )
assign(
    "updateMetabIntensityFilterQcPlot",
    function(filteredS4, omicType, setFilterPlotFn, ...) {
      captured$qc_plot_helper_call <- list(
        filtered_s4 = filteredS4,
        omic_type = omicType,
        set_filter_plot_fn = setFilterPlotFn
      )
      setFilterPlotFn("plot-token")
      invisible("plot-token")
    },
    envir = server_env
  )
assign(
    "reportMetabIntensityFilterSuccess",
    function(
      stats,
      intensityCutoffPercentile,
      proportionBelowCutoff,
      savedStateName,
      output,
      ...
    ) {
      captured$success_report_call <- list(
        stats = stats,
        intensity_cutoff_percentile = intensityCutoffPercentile,
        proportion_below_cutoff = proportionBelowCutoff,
        state_name = savedStateName,
        output = output
      )
      output$filter_results <- "success-report-token"
      invisible("success-report-token")
    },
    envir = server_env
  )
assign(
    "buildMetabIntensityFilterStats",
    function(currentS4, filteredS4) {
      captured$stats_call <- list(
        current_s4 = currentS4,
        filtered_s4 = filteredS4
      )
      expected_stats
    },
    envir = server_env
  )
scoped_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          apply_filter = TRUE,
          revert_filter = FALSE,
          intensity_cutoff_percentile = 10,
          proportion_below_cutoff = 0.5
        ),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      store <- new.env(parent = emptyenv())
      store$set_count <- 0L
      store$last_value <- value
      captured$reactive_values[[length(captured$reactive_values) + 1L]] <- store

      function(new_value) {
        if (missing(new_value)) {
          store$last_value
        } else {
          store$set_count <- store$set_count + 1L
          store$last_value <- new_value
          invisible(NULL)
        }
      }
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderUI = function(expr) substitute(expr),
    renderPlot = function(expr) substitute(expr),
    req = function(...) list(...)[[1]],
    showNotification = function(message, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- c(
        list(message = message),
        list(...)
      )
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    .env = asNamespace("shiny")
  )
scoped_mocked_bindings(
    log_info = function(message) invisible(NULL),
    log_error = function(message) {
      captured$error_log_message <- message
      invisible(NULL)
    },
    log_warn = function(message) {
      captured$warn_log_message <- message
      invisible(NULL)
    },
    .env = asNamespace("logger")
  )
workflow_data <- omicsTestOwner(current_s4, "metab_raw_data_s4", config = list(config = "value"))
testthat::local_mocked_bindings(
    metaboliteIntensityFiltering = get("metaboliteIntensityFiltering", envir = server_env, inherits = FALSE),
    .package = "MultiScholaR"
  )
mod_metab_qc_intensity_server(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )
expect_identical(captured$module_id, "intensity")
expect_equal(captured$filter_call$the_object@metabolite_data, current_s4@metabolite_data, ignore_attr = TRUE)
expect_identical(captured$filter_call$intensity_cutoff_percentile, 10)
expect_identical(captured$filter_call$proportion_below_cutoff, 0.5)
expect_equal(captured$stats_call$current_s4@metabolite_data, current_s4@metabolite_data, ignore_attr = TRUE)
