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

makeLipidQcS4CharacterizationData <- function(label = "qc_s4_fixture") {
  lipid_data <- list(
    Positive = data.frame(
      lipid_id = c(paste0(label, "_l1"), paste0(label, "_l2")),
      annotation_id = c("ann_a", "ann_b"),
      S1 = c(10, NA_real_),
      S2 = c(20, 30),
      stringsAsFactors = FALSE
    )
  )
  sample_columns <- c("S1", "S2")

  args <- list(
    Class = "LipidomicsAssayData",
    lipid_data = lipid_data,
    design_matrix = data.frame(
      Sample_ID = sample_columns,
      group = c("A", "B"),
      batch = c("B1", "B2"),
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
    args$internal_standard_regex <- NA_character_
  }
  if ("technical_replicate_id" %in% slots) {
    args$technical_replicate_id <- NA_character_
  }
  if ("args" %in% slots) {
    args$args <- list(label = label)
  }

  do.call(methods::new, args)
}

# KIT-072: the harness's bus holds a real lipidomics owner (the modules' only state manager) with the data committed
# under the import, duplicates and intensity states; the finalise commits the head's record under the QC state.
makeLipidQcS4WorkflowHarness <- function(
  current_s4 = makeLipidQcS4CharacterizationData(),
  history = c("import_complete", "duplicates_resolved", "intensity_filtered"),
  env = parent.frame()
) {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$config_list <- list(qc = "complete")
  workflow_data$tab_status <- list(
    quality_control = "pending",
    normalization = "locked"
  )
  omicsTestOwner(current_s4, history[[1L]], workflow = workflow_data, env = env)
  for (state in history[-1L]) workflow_data$state_manager$saveState(state, current_s4, workflow_data$config_list, state)

  list(
    workflow_data = workflow_data,
    state_manager = workflow_data$state_manager,
    current_s4 = current_s4
  )
}

test_that("mod_lipid_qc_s4_server preserves public finalize behavior", {
  harness <- makeLipidQcS4WorkflowHarness()
  mod_lipid_qc_s4_server <- getMultiScholaRBinding("mod_lipid_qc_s4_server")

  local_mocked_bindings(
    showNotification = function(...) invisible(NULL),
    .package = "shiny"
  )
  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )
  local_mocked_bindings(
    updateLipidFiltering = function(...) grid::nullGrob(),
    .env = multiScholaRNamespace()
  )

  testServer(
    mod_lipid_qc_s4_server,
    args = list(
      workflow_data = harness$workflow_data,
      omic_type = "lipidomics",
      experiment_label = "Lipidomics"
    ),
    {
      session$setInputs(finalize_qc = 1)
    }
  )

  expect_identical(harness$workflow_data$tab_status$quality_control, "complete")
  expect_identical(harness$state_manager$current_state, "lipid_qc_complete")
  expect_equal(harness$state_manager$getState()@lipid_data, harness$current_s4@lipid_data, ignore_attr = TRUE)
})
