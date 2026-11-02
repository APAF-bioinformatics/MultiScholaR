# KIT-075: the lipid and metab tabs hold references. After every step of each omic's journey, run under testServer on
# the real module server and session owner, no LipidomicsAssayData / MetaboliteAssayData is reachable from the
# session (every module's outputs, observers and reactive values) or the bus. The journey must reach every step's state.

# Item 1: the owned-value helpers take the owner explicitly. Without it they find only the proteomics session owner,
# so an omics value either passed through unsaved or landed in the proteomics owner's store.
test_that("an omics value is saved to the omics owner passed, with no proteomics owner in the session", {
  owner <- omicsTestOwner(omicsOwnerLipidObject(), "lipid_raw_data_s4")$state_manager
  table <- data.frame(sample = c("S1", "S2"), value = c(1.5, 2.5), stringsAsFactors = FALSE)
  shiny::withReactiveDomain(shiny::MockShinySession$new(), {
    expect_null(getWorkflowActiveOwner())
    holder <- new.env(parent = emptyenv())
    setWorkflowOwnedPayload(holder, "data_cln", table, manager = owner)
    expect_true(isSavedWorkflowReference(holder$data_cln))
    expect_identical(holder$data_cln$owner, owner$.context$store$identity)
    expect_equal(as.data.frame(getWorkflowOwnedPayload(holder, "data_cln", manager = owner)), table)
    expect_identical(getWorkflowOwnedInput(holder, "data_cln", manager = owner), holder$data_cln)
    value <- savedReactiveVal(table, "design_builder_table", manager = owner)
    expect_true(isSavedWorkflowReference(shiny::isolate(getWorkflowReactiveReference(value))))
    expect_equal(as.data.frame(shiny::isolate(value())), table)
  })
})

test_that("with a proteomics owner in the session, the omics value lands in the omics owner's store", {
  owner <- omicsTestOwner(omicsOwnerLipidObject(), "lipid_raw_data_s4")$state_manager
  session <- shiny::MockShinySession$new()
  proteomics <- new.env(parent = emptyenv())
  proteomics$config_list <- list()
  proteomics_root <- withr::local_tempdir()
  proteomics$state_manager <- createWorkflowSessionOwner(proteomics, list(source_dir = proteomics_root), session)
  withr::defer(closeWorkflowStateOwner(proteomics$state_manager))
  expect_identical(session$userData$dia_capacity_manager, proteomics$state_manager)
  table <- data.frame(sample = c("S1", "S2"), value = c(1.5, 2.5), stringsAsFactors = FALSE)
  shiny::withReactiveDomain(session, {
    holder <- new.env(parent = emptyenv())
    setWorkflowOwnedPayload(holder, "data_cln", table, manager = owner)
    expect_identical(holder$data_cln$owner, owner$.context$store$identity)
    expect_false(identical(holder$data_cln$owner, proteomics$state_manager$.context$store$identity))
  })
})

test_that("the walk finds a held object and a held assay table, through closures, reactive values and attributes", {
  object <- omicsOwnerLipidObject()
  table <- data.frame(id = 1:2, S1 = 1:2, S2 = 3:4)
  held <- new.env(parent = emptyenv())
  held$kept <- local({ captured <- object; function() captured })
  expect_identical(as.character(omicsResidencyWalk(list(root = held))), "root$kept<closure>$captured")
  values <- shiny::reactiveValues(assay = table)
  expect_true(length(attr(omicsResidencyWalk(list(root = values), table_columns = c("S1", "S2")), "tables")) == 1L)
  expect_length(attr(omicsResidencyWalk(list(root = list(meta = data.frame(id = 1:2, S1 = 1:2))), table_columns = c("S1", "S2")), "tables"), 0L)
  expect_length(omicsResidencyWalk(list(root = list(attr_free = 1:3))), 0L)
})

expected_states <- function(prefix) c(import = "initial", design = paste0(prefix, "_raw_data_s4"),
  intensity = paste0(prefix, "_intensity_filtered"), duplicates = paste0(prefix, "_duplicates_resolved"),
  qc_finalize = paste0(prefix, "_qc_complete"), normalisation = paste0(prefix, "_ruv_corrected"),
  correlation = paste0(prefix, "_correlation_filtered"), session_export = paste0(prefix, "_correlation_filtered"),
  da_load = paste0(prefix, "_correlation_filtered"), da = paste0(prefix, "_correlation_filtered"),
  summary = paste0(prefix, "_correlation_filtered"))

for (omic in c("lipidomics", "metabolomics")) {
  test_that(paste0(omic, ": no assay object is reachable after any journey step"), {
    hits <- suppressMessages(suppressWarnings(omicsResidencyJourney(omic)))
    prefix <- if (identical(omic, "lipidomics")) "lipid" else "metab"
    expect_identical(unlist(attr(hits, "states")), expected_states(prefix))
    # A hit is reported by its holder: the path after the last output, observer or reactive-values hop.
    holder <- function(path) sub("^.*(\\$private\\$outs\\[\\[|\\.values\\$private\\$map|\\$userData\\$)", "\\1", path)
    for (step in names(hits)) expect(length(hits[[step]]) == 0L,
      sprintf("%s: assay objects reachable at\n  %s", step, paste(unique(holder(hits[[step]])), collapse = "\n  ")))
    # The import and design tables are held by reference too (the import's own file-selection preview is outside this
    # walk, which drives the import through its processing helper; recorded in the outcome).
    tables <- attr(hits, "tables")
    for (step in names(hits)) {
      expect(length(tables[[step]]) == 0L, sprintf("%s: assay tables reachable at\n  %s", step,
        paste(unique(holder(tables[[step]])), collapse = "\n  ")))
      expect(all(attr(hits, "bus")[[step]] %in% c("reference", "empty")), sprintf("%s: the bus holds %s", step,
        paste(names(attr(hits, "bus")[[step]]), attr(hits, "bus")[[step]], collapse = ", ")))
    }
  })
}
