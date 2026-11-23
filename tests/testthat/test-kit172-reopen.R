# KIT-172 T9: a project saved by the base commit (offset 1; tests/testdata/kit-172/base-projects.tar.gz, written by
# docs/qa/workflow-execution/kit/kit-172/t9-base-projects.R in a worktree of that commit, its SHA in
# tests/testdata/kit-172/base-projects.sha) reopens, replays to identical states, and the normalisation tab shows the
# recorded-offset sentence; a new project records offset 0 and shows none.
library(testthat)

kit172RecordedSentence <- function(manager, prefix) {
  html <- NULL
  module <- function(id) shiny::moduleServer(id, function(input, output, session) {
    output$sentence <- omicNormRecordedOffsetUi(shiny::reactiveValues(state_manager = manager, tab_status = list()),
      shiny::reactiveValues(post_log2_obj = NULL), prefix)
  })
  shiny::testServer(module, {
    value <- output$sentence
    html <<- if (is.null(value) || !nzchar(paste(as.character(value$html), collapse = ""))) NULL else as.character(value$html)
  })
  html
}

kit172Paths <- function(root, class) {
  qc <- if (identical(class, "lipid")) "lipid_qc_dir" else "metabolite_qc_dir"
  stats::setNames(list(file.path(root, "source"), file.path(root, "qc"), file.path(root, "da"), file.path(root, "pub")),
    c("source_dir", qc, "da_output_dir", "publication_graphs_dir"))
}

for (class in c("lipid", "metab")) {
  test_that(paste("KIT-172 T9:", class, "a project saved by the base commit reopens, replays exactly and shows its recorded offset"), {
    skip_if_not(workflowOwnerCapable())
    archive <- testthat::test_path("..", "testdata", "kit-172", "base-projects.tar.gz")
    expect_true(file.exists(archive))
    unpacked <- withr::local_tempdir()
    utils::untar(archive, exdir = unpacked)
    root <- file.path(unpacked, class)
    prefix <- kit172Binding(class)$prefix
    source_dir <- file.path(root, "source")
    store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    identity <- readSavedStoreIdentity(file.path(store_root, ".workflow-entities", "store.json"))
    owner <- newDIAStateAdapter(newWorkflowEntityStore(store_root, identity$project_id, identity$workflow_id), reopen = TRUE)
    expect_identical(omicNormRecordedLog2Offset(owner, prefix), 1)
    html <- kit172RecordedSentence(owner, prefix)
    expect_true(!is.null(html) && grepl(sprintf(omicGuiLog2Texts()[["recorded_offset"]], "1"), html, fixed = TRUE))
    closeDIAStateAdapter(owner)
    # The replay on a copy cut to the design state: every state identical to the one the base commit wrote.
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
    harness <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
    copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
    prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(prefix, "_raw_data_s4")), api)
    inputs <- c(stats::setNames(file.path(source_dir, "import_sources", "LCMS_Pos_export.txt"), paste0(prefix, "_assay_1")),
      design_file = file.path(source_dir, "design_matrix.tab"))
    record <- jsonlite::read_json(file.path(root, "analysis_record.json"), simplifyVector = FALSE)
    replay <- kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, kit172Paths(root, class),
      withr::local_tempdir(), api))
    expect_true(paste0(prefix, "_log2") %in% replay$states)
    expect_identical(as.character(replay$mismatches), character())
  })

  test_that(paste("KIT-172 T9:", class, "a new project records offset 0 and shows no recorded-offset sentence"), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    values <- kit172Values(172001L)$values[c(1:40, 516:520), ]
    journey <- kit172Journey(file.path(root, "j"), kit172Object(values, class), class, da = FALSE)
    withr::defer(closeDIAStateAdapter(journey$owner))
    expect_identical(omicNormRecordedLog2Offset(journey$owner, journey$prefix), 0)
    expect_null(kit172RecordedSentence(journey$owner, journey$prefix))
  })
}
