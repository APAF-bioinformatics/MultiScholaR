# KIT-070: the on-disk state owner for lipidomics and metabolomics. Item 1: the owner layer accepts the two omics
# (identity, route, session registration per omic, the metadata reader's class list, the analysis record's omic);
# a real owner in a temporary project commits, reads and reverts the declared classes.
library(testthat)

omicsOwnerSession <- function() {
  session <- new.env(parent = emptyenv())
  session$userData <- new.env(parent = emptyenv())
  session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
  session
}

omicsOwnerDesign <- function() {
  data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
}

omicsOwnerLipidObject <- function(n = 6L) {
  assay <- function(prefix) {
    data.frame(lipid_id = paste0(prefix, "_L", seq_len(n)), lipid_name = paste("Lipid", seq_len(n)),
      S1 = seq_len(n) + 0.5, S2 = seq_len(n) + 1.5, S3 = seq_len(n) + 2.5, S4 = seq_len(n) + 3.5,
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  methods::new("LipidomicsAssayData", lipid_data = list(pos = assay("pos"), neg = assay("neg")),
    design_matrix = omicsOwnerDesign(), sample_id = "Run", group_id = "group",
    lipid_id_column = "lipid_id", annotation_id_column = "lipid_name")
}

omicsOwnerMetabObject <- function(n = 5L) {
  assay <- function(prefix) {
    data.frame(metabolite_id = paste0(prefix, "_M", seq_len(n)), metabolite_name = paste("Met", seq_len(n)),
      S1 = seq_len(n) * 2, S2 = seq_len(n) * 3, S3 = seq_len(n) * 4, S4 = seq_len(n) * 5,
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  methods::new("MetaboliteAssayData", metabolite_data = list(LCMS_Pos = assay("pos"), GCMS = assay("gc")),
    design_matrix = omicsOwnerDesign(), sample_id = "Run", group_id = "group",
    metabolite_id_column = "metabolite_id", annotation_id_column = "metabolite_name")
}

test_that("KIT-070: the metadata reader's class list is every declared assay class", {
  expect_setequal(workflowDeclaredClassNames("assay"),
    c("PeptideQuantitativeData", "ProteinQuantitativeData", "LipidomicsAssayData", "MetaboliteAssayData"))
  expect_true(all(c("LipidomicsDifferentialAbundanceResults", "EnrichmentResults") %in% workflowDeclaredClassNames("result")))
  expect_identical(workflowRouteForType("LIPIDOMICS"), "lipidomics")
  expect_identical(workflowRouteForType("METABOLOMICS"), "metabolomics")
  expect_identical(workflowNativeStoreLabel(NULL, "LIPIDOMICS"), "LIPIDOMICS-native")
  expect_identical(workflowNativeStoreLabel("msdial", "METABOLOMICS"), "MSDIAL-native")
})

test_that("KIT-070: a lipidomics session owner carries its omic, leaves the proteomics slot alone and binds at the design step", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "source"); dir.create(source_dir)
  session <- omicsOwnerSession()
  workflow <- list2env(list(state_update_trigger = NULL))
  manager <- createDIASessionOwner(workflow, list(source_dir = source_dir), session, omic_type = "lipidomics")
  withr::defer(closeDIAStateAdapter(workflow$state_manager %||% manager))
  workflow$state_manager <- manager
  expect_s3_class(manager, "DIAStateAdapter")
  expect_true(isDIASessionOwnerPending(manager))
  expect_identical(manager$.context$store$identity$project_id, "LIPIDOMICS-pending")
  expect_identical(manager$.context$omic_type, "lipidomics")
  expect_null(session$userData$dia_capacity_manager)
  expect_true(usesDIAReferenceStates(manager))
  # The import records the type on the pending owner; the design step binds it.
  manager$setWorkflowType("LIPIDOMICS")
  expect_identical(manager$workflow_type, "LIPIDOMICS")
  expect_error(manager$setWorkflowType("lipidomics_standard"), class = "saved_entity_scope")
  placeholder_root <- dirname(manager$.context$store$directory)
  bound <- shiny::withReactiveDomain(session, bindWorkflowOwnerTypeAtDesign(workflow, "LIPIDOMICS", "lipidsearch"))
  expect_false(identical(bound, manager))
  expect_identical(workflow$state_manager, bound)
  expect_false(isDIASessionOwnerPending(bound))
  expect_false(dir.exists(placeholder_root))
  expect_identical(bound$.context$store$identity$workflow_id, "LIPIDOMICS")
  expect_identical(bound$.context$store$identity$project_id, "LIPIDSEARCH-native")
  expect_identical(bound$.context$omic_type, "lipidomics")
  expect_null(session$userData$dia_capacity_manager)
  # The same type again passes (a second design pass); another type is refused.
  expect_identical(shiny::withReactiveDomain(session, bindWorkflowOwnerTypeAtDesign(workflow, "LIPIDOMICS", "lipidsearch")), bound)
  expect_error(shiny::withReactiveDomain(session, bindWorkflowOwnerTypeAtDesign(workflow, "METABOLOMICS")), class = "saved_entity_scope")
})

omicsOwnerRoundTrip <- function(omic_type, workflow_type, object, filtered, assay_slot, states) {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "source"); dir.create(source_dir)
  session <- omicsOwnerSession()
  workflow <- list2env(list(state_update_trigger = NULL, config_list = list(threshold = 0.5)))
  workflow$state_manager <- createDIASessionOwner(workflow, list(source_dir = source_dir), session, omic_type = omic_type)
  withr::defer(closeDIAStateAdapter(workflow$state_manager))
  shiny::withReactiveDomain(session, bindWorkflowOwnerTypeAtDesign(workflow, workflow_type))
  owner <- workflow$state_manager
  # Every module save site calls the owner's saveState: a reference commit of the declared class.
  owner$saveState(states[[1]], object, workflow$config_list, "the design's object")
  record <- readWorkflowStateMetadata(owner)
  expect_true(isDIAComponentCheckpoint(record$data))
  expect_identical(record$data$class, class(object)[[1]])
  expect_identical(sort(names(record$data$components[[assay_slot]])), sort(names(methods::slot(object, assay_slot))))
  expect_identical(record$config, workflow$config_list)
  hydrated <- owner$getState()
  expect_s4_class(hydrated, class(object)[[1]])
  for (name in names(methods::slot(object, assay_slot))) {
    expect_equal(as.data.frame(methods::slot(hydrated, assay_slot)[[name]]), methods::slot(object, assay_slot)[[name]], ignore_attr = TRUE)
  }
  expect_equal(as.data.frame(hydrated@design_matrix), object@design_matrix, ignore_attr = TRUE)
  expect_identical(hydrated@sample_id, object@sample_id)
  owner$saveState(states[[2]], filtered, workflow$config_list, "a filtered object")
  expect_identical(owner$current_state, states[[2]])
  expect_identical(as.character(unlist(owner$getHistory())), c("initial", states))
  expect_identical(nrow(methods::slot(owner$getState(), assay_slot)[[1]]), nrow(methods::slot(filtered, assay_slot)[[1]]))
  # Revert by name: the current state and the history move back; the reverted object is the first commit's.
  reverted <- owner$revertToState(states[[1]])
  expect_identical(owner$current_state, states[[1]])
  expect_identical(as.character(unlist(owner$getHistory())), c("initial", states[[1]]))
  expect_identical(nrow(methods::slot(reverted, assay_slot)[[1]]), nrow(methods::slot(object, assay_slot)[[1]]))
  # The analysis record of this owner carries the omic and validates.
  analysis <- buildWorkflowAnalysisRecord(owner)
  expect_identical(analysis$omic_type, omic_type)
  expect_true(validateWorkflowAnalysisRecord(analysis)$valid)
  invisible(owner)
}

test_that("KIT-070: a lipidomics owner commits, reads and reverts LipidomicsAssayData states", {
  object <- omicsOwnerLipidObject()
  filtered <- object
  filtered@lipid_data <- lapply(filtered@lipid_data, function(df) df[1:3, , drop = FALSE])
  omicsOwnerRoundTrip("lipidomics", "LIPIDOMICS", object, filtered, "lipid_data",
    c("lipid_raw_data_s4", "lipid_intensity_filtered"))
})

test_that("KIT-070: a metabolomics owner commits, reads and reverts MetaboliteAssayData states", {
  object <- omicsOwnerMetabObject()
  filtered <- object
  filtered@metabolite_data <- lapply(filtered@metabolite_data, function(df) df[1:2, , drop = FALSE])
  omicsOwnerRoundTrip("metabolomics", "METABOLOMICS", object, filtered, "metabolite_data",
    c("metab_raw_data_s4", "metab_intensity_filtered"))
})

# ---------------------------------------------------------------------------
# Item 2: both module servers hold the session owner from session start; no in-memory manager remains.
# ---------------------------------------------------------------------------
test_that("KIT-070: the lipidomics and metabolomics module servers create the session owner at start", {
  skip_if_not(workflowOwnerCapable())
  for (omic in c("lipidomics", "metabolomics")) {
    root <- withr::local_tempdir()
    source_dir <- file.path(root, "source"); dir.create(source_dir)
    dirs <- stats::setNames(list(list(source_dir = source_dir)), omic)
    server <- get(paste0("mod_", omic, "_server"), envir = asNamespace("MultiScholaR"), inherits = FALSE)
    shiny::testServer(server, args = list(project_dirs = dirs, omic_type = omic, experiment_label = "Owner test", volumes = NULL), {
      owner <- shiny::isolate(workflow_data$state_manager)
      expect_s3_class(owner, "DIAStateAdapter")
      expect_true(isDIASessionOwnerPending(owner))
      expect_identical(owner$.context$store$identity$project_id, paste0(toupper(omic), "-pending"))
      expect_identical(owner$.context$omic_type, omic)
      expect_identical(owner$.context$experiment_paths$source_dir, source_dir)
      expect_null(session$userData$dia_capacity_manager)
      expect_identical(session$userData$omics_workflow_data[[omic]], workflow_data)
      closeDIAStateAdapter(owner)
    })
  }
})

test_that("KIT-070: no lipidomics or metabolomics module constructs an in-memory state manager", {
  r_dir <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = TRUE)
  files <- list.files(r_dir, "^mod_(lipid|metab)", full.names = TRUE)
  expect_gt(length(files), 20L)
  hits <- unlist(lapply(files, function(f) {
    lines <- readLines(f, warn = FALSE)
    code <- sub("#.*$", "", lines)
    idx <- grep("WorkflowState\\$new|workflowStateClass|stateManagerFactory", code)
    if (length(idx)) paste0(basename(f), ":", idx) else character()
  }))
  expect_identical(hits, character())
})
