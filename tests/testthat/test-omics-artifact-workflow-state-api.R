workflowStateRoundTripCases <- function() {
  list(
    proteomics = list(
      workflow_type = "DIA",
      object = module_ci_prot_norm_object()
    ),
    metabolomics = list(
      workflow_type = "metabolomics_standard",
      object = module_ci_metab_norm_object()
    ),
    lipidomics = list(
      workflow_type = "lipidomics_standard",
      object = module_ci_lipid_norm_object()
    )
  )
}

test_that("WorkflowState exposes compatible state and metadata methods", {
  manager <- WorkflowState$new(audit_enabled = FALSE)
  expect_identical(manager$getCurrentStateName(), "initial")
  expect_identical(manager$getWorkflowType(), "LFQ")
  expect_false(manager$isAuditEnabled())
  expect_true(manager$hasState("initial"))
  expect_false(manager$hasState("missing"))

  object <- list(marker = "payload")
  expect_identical(
    manager$saveState("filtered", object, list(threshold = 0.05), "filtered"),
    "filtered"
  )
  expect_identical(manager$getState(), object)
  expect_identical(manager$getStateConfig(), list(threshold = 0.05))
  expect_false("data" %in% names(manager$getStateMetadata()))
  expect_identical(manager$getCurrentStateName(), "filtered")
})

test_that("active lineage and immutable transition events remain separate", {
  manager <- WorkflowState$new()
  manager$saveState("first", list(value = 1L), list(), "first")
  manager$saveState("second", list(value = 2L), list(), "second")
  eventsBeforeRevert <- manager$getEvents()

  manager$revertToState("first")

  expect_identical(manager$getHistory(), c("initial", "first"))
  expect_true(manager$hasState("second"))
  expect_identical(manager$getState("second"), list(value = 2L))
  expect_identical(manager$revertToState("second"), list(value = 2L))
  expect_identical(manager$getHistory(), c("initial", "first", "second"))
  expect_identical(
    vapply(eventsBeforeRevert, `[[`, character(1), "event_type"),
    c("initialized", "state_saved", "state_saved")
  )
  expect_identical(
    tail(manager$getEvents(), 1L)[[1L]]$event_type,
    "state_reverted"
  )
})

test_that("WorkflowState observers receive transitions and can unsubscribe", {
  manager <- WorkflowState$new()
  observed <- list()
  unsubscribe <- manager$observeTransitions(function(event) {
    observed[[length(observed) + 1L]] <<- event
  })

  manager$saveState("first", list(), list(), "first")
  manager$setWorkflowType("DIA")
  unsubscribe()
  unsubscribe()
  manager$saveState("second", list(), list(), "second")

  expect_identical(
    vapply(observed, `[[`, character(1), "event_type"),
    c("state_saved", "workflow_type_set")
  )
  expect_identical(vapply(observed, `[[`, integer(1), "revision"), 2:3)
})

test_that("versioned WorkflowState manifests round trip every implemented omic", {
  for (omicName in names(workflowStateRoundTripCases())) {
    case <- workflowStateRoundTripCases()[[omicName]]
    manager <- WorkflowState$new()
    manager$setWorkflowType(case$workflow_type)
    manager$saveState(
      "imported",
      case$object,
      list(omic = omicName, stage = "import"),
      "imported"
    )
    auditRecord <- list(
      record_id = paste0(omicName, "-record"),
      stage_id = "normalization"
    )
    manager$saveState(
      "normalized",
      case$object,
      list(omic = omicName, method = "none"),
      "normalized",
      audit_metadata = auditRecord
    )
    manager$revertToState("imported")
    manifest <- manager$exportState()

    restored <- WorkflowState$new(audit_enabled = FALSE)
    restored$restoreState(manifest)

    expect_identical(restored$getCurrentStateName(), manager$getCurrentStateName())
    expect_identical(restored$getHistory(), manager$getHistory())
    expect_identical(restored$getWorkflowType(), manager$getWorkflowType())
    expect_identical(restored$isAuditEnabled(), manager$isAuditEnabled())
    expect_identical(restored$getAuditRecords(), manager$getAuditRecords())
    expect_identical(restored$getState("imported"), manager$getState("imported"))
    expect_true(methods::validObject(restored$getState("imported")))
    expect_identical(
      restored$getState("normalized"),
      manager$getState("normalized")
    )
    expect_identical(
      restored$getStateMetadata("normalized"),
      manager$getStateMetadata("normalized")
    )
    expect_identical(
      head(restored$getEvents(), length(manifest$events)),
      manifest$events
    )
    expect_identical(
      tail(restored$getEvents(), 1L)[[1L]]$event_type,
      "state_restored"
    )
  }
})

test_that("manifest validation rejects invalid input before state mutation", {
  manager <- WorkflowState$new()
  manager$saveState("ready", list(value = 1L), list(), "ready")
  baseline <- manager$exportState()

  invalidManifests <- list(
    future = within(baseline, schema_version <- 999L),
    backend = within(baseline, backend <- "artifact"),
    current = within(baseline, current_state <- "missing"),
    lineage = within(baseline, active_lineage <- "initial"),
    states = within(baseline, states[["ready"]] <- list(data = list())),
    events = within(baseline, events[[1L]]$revision <- 99L),
    eventHead = within(baseline, events[[length(events)]]$state_name <- "initial"),
    missingEvents = within(baseline, events <- list()),
    embeddedVersion = within(baseline, schema_version <- 1.5)
  )
  expectedMessages <- c(
    "future WorkflowState schema",
    "backend is invalid",
    "current state is inconsistent",
    "active lineage is inconsistent",
    "malformed state record",
    "malformed transition events",
    "transition head is inconsistent",
    "transition events are missing",
    "schema version is missing"
  )

  for (index in seq_along(invalidManifests)) {
    expect_error(manager$restoreState(invalidManifests[[index]]), expectedMessages[[index]])
    expect_identical(manager$exportState(), baseline)
  }
  expect_error(manager$restoreState(list()), "schema version is missing")
  expect_error(
    manager$restoreState(
      within(baseline, schema_version <- 1.5),
      schema_version = 1L
    ),
    "schema version is inconsistent"
  )
  expect_identical(manager$exportState(), baseline)
})

test_that("audit record identifiers cannot be reused with changed metadata", {
  manager <- WorkflowState$new()
  manager$saveState(
    "first",
    list(value = 1L),
    list(),
    "first",
    audit_metadata = list(record_id = "audit-1", threshold = 0.05)
  )
  baseline <- manager$exportState()

  expect_error(
    manager$saveState(
      "second",
      list(value = 2L),
      list(),
      "second",
      audit_metadata = list(record_id = "audit-1", threshold = 0.1)
    ),
    "audit record IDs must be immutable"
  )
  expect_identical(manager$exportState(), baseline)

  inconsistent <- baseline
  inconsistent$states$first$audit_metadata$threshold <- 0.1
  expect_error(
    manager$restoreState(inconsistent),
    "state audits are inconsistent"
  )
  expect_identical(manager$exportState(), baseline)
})

test_that("legacy session lists restore only through the explicit legacy version", {
  object <- module_ci_prot_norm_object()
  legacy <- list(
    r6_current_state_name = "filtered",
    r6_complete_states = list(filtered = object),
    r6_state_history = "filtered",
    current_s4_object = object,
    workflow_type = "DIA"
  )
  manager <- WorkflowState$new()

  expect_error(manager$restoreState(legacy), "schema version is missing")
  manager$restoreState(legacy, schema_version = 0L)

  expect_identical(manager$getCurrentStateName(), "filtered")
  expect_identical(manager$getHistory(), "filtered")
  expect_identical(manager$getState(), object)
  expect_identical(manager$getWorkflowType(), "DIA")
  expect_identical(manager$getEvents()[[1L]]$event_type, "state_restored")
})

test_that("legacy four-argument state-manager doubles remain compatible", {
  manager <- new.env(parent = emptyenv())
  manager$states <- list(initial = list(value = 1L))
  manager$state_history <- "initial"
  manager$current_state <- "initial"
  manager$workflow_type <- "DIA"
  manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    manager$states[[state_name]] <- s4_data_object
    manager$state_history <- c(manager$state_history, state_name)
    manager$current_state <- state_name
  }

  manager$saveState("filtered", list(value = 2L), list(), "filtered")
  snapshot <- workflowStateLegacySnapshot(manager)

  expect_identical(snapshot$r6_current_state_name, "filtered")
  expect_identical(snapshot$r6_state_history, c("initial", "filtered"))
  expect_identical(workflowStateType(manager), "DIA")
})

test_that("production callers do not access WorkflowState fields directly", {
  sourceFiles <- list.files(
    testthat::test_path("..", "..", "R"),
    pattern = "[.]R$",
    full.names = TRUE
  )
  ownerFiles <- c("utils_workflow_state.R", "utils_workflow_state_compat.R")
  sourceFiles <- sourceFiles[!basename(sourceFiles) %in% ownerFiles]
  sourceLines <- unlist(lapply(sourceFiles, readLines, warn = FALSE), use.names = FALSE)
  directAccess <- paste0(
    "(state_manager|stateManager)[$]",
    "(states|current_state|state_history|audit_records|audit_enabled|workflow_type)\\b"
  )

  expect_false(any(grepl(directAccess, sourceLines, perl = TRUE)))
})
