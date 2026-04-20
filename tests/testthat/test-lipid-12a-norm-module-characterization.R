library(testthat)
library(shiny)

if (!methods::isClass("LipidomicsAssayData")) {
  methods::setClass(
    "LipidomicsAssayData",
    slots = c(
      lipid_data = "list",
      lipid_id_column = "character",
      annotation_id_column = "character",
      database_identifier_type = "character",
      internal_standard_regex = "character",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character",
      technical_replicate_id = "character",
      args = "list"
    ),
    prototype = list(
      lipid_data = list(),
      lipid_id_column = "database_identifier",
      annotation_id_column = "lipid_identification",
      database_identifier_type = "Unknown",
      internal_standard_regex = NA_character_,
      design_matrix = data.frame(),
      sample_id = "Sample_ID",
      group_id = "group",
      technical_replicate_id = NA_character_,
      args = list()
    )
  )
}

makeLipidNormCharacterizationData <- function() {
  methods::new(
    "LipidomicsAssayData",
    lipid_data = list(
      `Positive Mode` = data.frame(
        database_identifier = c("L1", "L2"),
        lipid_identification = c("A1", "A2"),
        Sample1 = c(10, 20),
        Sample2 = c(30, 40),
        stringsAsFactors = FALSE
      )
    ),
    lipid_id_column = "database_identifier",
    annotation_id_column = "lipid_identification",
    database_identifier_type = "Mock",
    internal_standard_regex = "",
    design_matrix = data.frame(
      Sample_ID = c("Sample1", "Sample2"),
      group = c("A", "B"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Sample_ID",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = list()
  )
}

makeLipidNormCharacterizationHarness <- function(current_s4 = makeLipidNormCharacterizationData()) {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- list(
    getState = function(...) current_s4
  )
  workflow_data$design_matrix <- data.frame(
    Sample_ID = c("Sample1", "Sample2"),
    group = c("A", "B"),
    batch = c("B1", "B2"),
    stringsAsFactors = FALSE
  )

  list(
    workflow_data = workflow_data,
    experiment_paths = list(lipid_qc_dir = tempdir()),
    session = shiny::MockShinySession$new()
  )
}

launchLipidNormModule <- function(harness, selected_tab = NULL) {
  shiny::withReactiveDomain(harness$session, {
    mod_lipid_norm_server(
      id = "norm",
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "lipidomics",
      experiment_label = "Lipidomics",
      selected_tab = selected_tab
    )
  })

  harness$session$getReturned()
}

renderLipidNormOutput <- function(session, output_id) {
  htmltools::renderTags(session$getOutput(paste0("norm-", output_id)))$html
}

test_that("mod_lipid_norm_server preserves startup log placeholder behavior", {
  harness <- makeLipidNormCharacterizationHarness()

  local_mocked_bindings(
    generateLipidQcPlots = function(...) invisible(NULL),
    .env = asNamespace("MultiScholaR")
  )

  launchLipidNormModule(harness)

  norm_log_html <- renderLipidNormOutput(harness$session, "norm_log")

  expect_match(
    norm_log_html,
    "Normalization log will appear here as you apply steps...",
    fixed = TRUE
  )
})

test_that("mod_lipid_norm_server preserves selected-tab pre-QC error logging behavior", {
  harness <- makeLipidNormCharacterizationHarness()
  selected_tab_value <- shiny::reactiveVal("overview")

  local_mocked_bindings(
    generateLipidQcPlots = function(...) stop("pre qc boom"),
    .env = asNamespace("MultiScholaR")
  )

  launchLipidNormModule(
    harness,
    selected_tab = function() selected_tab_value()
  )

  selected_tab_value("norm")
  shiny:::flushReact()

  norm_log_html <- renderLipidNormOutput(harness$session, "norm_log")

  expect_match(norm_log_html, "Error generating Pre-QC: pre qc boom", fixed = TRUE)
})
