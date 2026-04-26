# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

multiScholaRNamespace <- function() {
  asNamespace("MultiScholaR")
}

getMultiScholaRBinding <- function(name) {
  get(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

if (!methods::isClass("LipidomicsAssayData")) {
  methods::setClass(
    "LipidomicsAssayData",
    slots = c(
      lipid_data = "list",
      design_matrix = "data.frame",
      group_id = "character",
      lipid_id_column = "character",
      sample_id = "character",
      annotation_id_column = "character",
      database_identifier_type = "character",
      internal_standard_regex = "character",
      technical_replicate_id = "character",
      args = "list"
    ),
    prototype = list(
      lipid_data = list(),
      design_matrix = data.frame(),
      group_id = "group",
      lipid_id_column = "lipid_id",
      sample_id = "Sample_ID",
      annotation_id_column = "annotation_id",
      database_identifier_type = "test_id",
      internal_standard_regex = NA_character_,
      technical_replicate_id = NA_character_,
      args = list()
    )
  )
}

makeLipidQcItsdCharacterizationData <- function() {
  sample_columns <- c("S1", "S2")
  lipid_data <- list(
    Positive = data.frame(
      lipid_id = c("IS_anchor", "Lipid_A"),
      annotation_id = c("internal_standard", "lipid_a"),
      S1 = c(100, 50),
      S2 = c(110, 55),
      stringsAsFactors = FALSE
    )
  )

  args <- list(
    Class = "LipidomicsAssayData",
    lipid_data = lipid_data,
    design_matrix = data.frame(
      Sample_ID = sample_columns,
      group = c("A", "A"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Sample_ID",
    group_id = "group",
    lipid_id_column = "lipid_id"
  )

  slots <- methods::slotNames("LipidomicsAssayData")
  if ("annotation_id_column" %in% slots) {
    args$annotation_id_column <- "annotation_id"
  }
  if ("database_identifier_type" %in% slots) {
    args$database_identifier_type <- "test_id"
  }
  if ("internal_standard_regex" %in% slots) {
    args$internal_standard_regex <- "^IS_"
  }
  if ("technical_replicate_id" %in% slots) {
    args$technical_replicate_id <- NA_character_
  }
  if ("args" %in% slots) {
    args$args <- list(label = "lipid_qc_itsd_fixture")
  }

  do.call(methods::new, args)
}

makeLipidQcItsdWorkflowHarness <- function(current_s4 = makeLipidQcItsdCharacterizationData()) {
  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- current_s4
  state_manager$getState <- function() state_manager$current_state

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager

  list(
    workflow_data = workflow_data,
    state_manager = state_manager,
    current_s4 = current_s4
  )
}

test_that("mod_lipid_qc_itsd_server preserves public internal-standard analysis behavior", {
  skip_if_not(
    exists("mod_lipid_qc_itsd_server", envir = multiScholaRNamespace(), inherits = FALSE),
    "mod_lipid_qc_itsd_server is unavailable in this ref."
  )

  server <- getMultiScholaRBinding("mod_lipid_qc_itsd_server")
  harness <- makeLipidQcItsdWorkflowHarness()
  notifications <- new.env(parent = emptyenv())
  notifications$shown <- list()
  notifications$removed <- character()

  local_mocked_bindings(
    showNotification = function(ui, id = NULL, duration = 5, type = "default", ...) {
      notifications$shown[[length(notifications$shown) + 1L]] <- list(
        ui = as.character(ui),
        id = id,
        duration = duration,
        type = type
      )
      invisible(if (is.null(id)) "notification" else id)
    },
    removeNotification = function(id, ...) {
      notifications$removed <- c(notifications$removed, id)
      invisible(NULL)
    },
    .package = "shiny"
  )
  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    server,
    args = list(
      workflow_data = harness$workflow_data,
      omic_type = "lipidomics",
      experiment_label = "Lipidomics"
    ),
    {
      session$setInputs(is_pattern = "^IS_")
      session$setInputs(analyze_is = 1)
      session$flushReact()

      expect_match(output$is_results, "Internal Standard Analysis Complete", fixed = TRUE)
      expect_match(output$is_results, "Pattern used: ^IS_", fixed = TRUE)
      expect_match(output$is_results, "Total IS detected: 1", fixed = TRUE)
      expect_match(output$is_results, "Positive: 1 internal standards", fixed = TRUE)
      expect_match(htmltools::renderTags(output$is_summary)$html, "Positive: 1 IS", fixed = TRUE)
      expect_match(htmltools::renderTags(output$is_viz_tabs)$html, "CV Distribution", fixed = TRUE)
    }
  )

  expect_identical(notifications$removed, "is_analysis_working")
  expect_true(any(vapply(
    notifications$shown,
    function(entry) identical(entry$ui, "Found 1 internal standards"),
    logical(1)
  )))
})
