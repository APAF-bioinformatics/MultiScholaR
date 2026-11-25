# KIT-258 (VAL-254): lipidomics and metabolomics duplicate resolution keeps, among the rows of one duplicated id, the
# row measured (not NA) in the most of the averaged columns (the design's samples, KIT-214), then the highest mean of
# its measured values, then the first in the table (today's tie rule, Q-C1). Tests E1 to E4, G1 to G3, T1 (a), T1 (b)
# and T2 of the ticket's section 10. E5 (zeros, OD-1) and E6 (the 12-sample example, OD-2) are written after the
# owner's answers and are not here. Every "not measured" value below is NA; no duplicated row holds a zero.
# Base values (KIT-214's rule, the highest mean alone) are the ticket's: E1 to E4 keep row A (E3 row Q).

# ---- Helper section: the tables -------------------------------------------------------------------------------------

# Six design samples, two groups of three.
k258Samples <- sprintf("S%d", 1:6)
k258Design <- function() {
  groups <- rep(c("A", "B"), each = 3L)
  data.frame(Run = k258Samples, group = groups, tech_rep_group = paste(groups, c(1:3, 1:3), sep = "_"),
    stringsAsFactors = FALSE)
}

# Rows from a list of (feature_id, name, six values); `name` is text, unique per row (the row's annotation).
k258Table <- function(...) {
  rows <- list(...)
  values <- do.call(rbind, lapply(rows, function(row) row$values))
  colnames(values) <- k258Samples
  cbind(data.frame(feature_id = vapply(rows, function(row) row$id, ""), name = vapply(rows, function(row) row$name, ""),
    stringsAsFactors = FALSE), as.data.frame(values))
}
k258Row <- function(id, name, values) list(id = id, name = name, values = values)

# E1: id D, row A measured in 2 of 6 (mean 210,000), row B in 6 of 6 (mean 150,000); id E once.
k258RowA <- function() k258Row("D", "A", c(200000, 220000, NA, NA, NA, NA))
k258RowB <- function() k258Row("D", "B", c(140000, 145000, 150000, 150000, 155000, 160000))
k258RowC <- function() k258Row("E", "C", c(1000, 1100, 1200, 1300, 1400, 1500))
k258E1Table <- function() k258Table(k258RowA(), k258RowB(), k258RowC())
# E2: E1 with B listed first.
k258E2Table <- function() k258Table(k258RowB(), k258RowA(), k258RowC())
# E3: id F, three rows: P 6/6 mean 100; Q 2/6 mean 500; R 6/6 mean 120; id G once.
k258E3Table <- function() k258Table(
  k258Row("F", "P", c(90, 95, 100, 100, 105, 110)),
  k258Row("F", "Q", c(NA, 450, NA, 550, NA, NA)),
  k258Row("F", "R", c(110, 115, 120, 120, 125, 130)),
  k258Row("G", "S", c(10, 11, 12, 13, 14, 15)))
# G1: equal completeness, U 6/6 at 120 and V 6/6 at 150 (the mean decides).
k258G1Table <- function() k258Table(
  k258Row("H", "U", c(110, 115, 120, 120, 125, 130)),
  k258Row("H", "V", c(140, 145, 150, 150, 155, 160)),
  k258Row("J", "Y", c(10, 11, 12, 13, 14, 15)))
# G3: equal completeness and equal mean (each mean exactly 100): W listed first.
k258G3Table <- function() k258Table(
  k258Row("K", "W", c(90, 95, 100, 100, 105, 110)),
  k258Row("K", "X", c(100, 100, 100, 100, 100, 100)),
  k258Row("L", "Z", c(10, 11, 12, 13, 14, 15)))

# The ticket's cases on the declared step: table, duplicated id, the name of the row kept now and under KIT-214.
k258Cases <- list(
  E1 = list(table = k258E1Table, id = "D", head = "B", base = "A"),
  E2 = list(table = k258E2Table, id = "D", head = "B", base = "A"),
  E3 = list(table = k258E3Table, id = "F", head = "R", base = "Q"),
  G1 = list(table = k258G1Table, id = "H", head = "V", base = "V"),
  G3 = list(table = k258G3Table, id = "K", head = "W", base = "W"))

k258Build <- function(omic, table) {
  if (identical(omic, "lipidomics")) {
    createLipidomicsAssayData(lipid_data = list(LCMS_Pos = table), design_matrix = k258Design(),
      lipid_id_column = "feature_id", annotation_id_column = "name", sample_id = "Run", group_id = "group",
      technical_replicate_id = "tech_rep_group")
  } else {
    createMetaboliteAssayData(metabolite_data = list(LCMS_Pos = table), design_matrix = k258Design(),
      metabolite_id_column = "feature_id", annotation_id_column = "name", sample_id = "Run", group_id = "group",
      technical_replicate_id = "tech_rep_group")
  }
}

# ---- Helper section: reading an object, running the step (after KIT-214's test file) --------------------------------

k258Lipid <- function(object) methods::is(object, "LipidomicsAssayData")
k258Binding <- function(object) if (k258Lipid(object)) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
k258Assays <- function(object) methods::slot(object, k258Binding(object)$data_slot)

# The object with a text column `.row` (the row's position) in each assay: neither numeric nor a design sample.
k258WithRow <- function(object) {
  assays <- lapply(k258Assays(object), function(table) {
    table$.row <- as.character(seq_len(nrow(table)))
    table
  })
  methods::slot(object, k258Binding(object)$data_slot, check = FALSE) <- assays
  object
}

# The kept row of `id` in a resolved table: its name (annotation), or how many rows the id has when not one.
k258KeptName <- function(table, id) {
  rows <- table[as.character(table$feature_id) == id, , drop = FALSE]
  if (nrow(rows) != 1L) return(sprintf("%d rows of id %s", nrow(rows), id))
  rows$name
}

# The kept row of `id` in the first assay of a resolved object, read by `.row` and by its name.
k258Kept <- function(object, id) {
  table <- k258Assays(object)[[1L]]
  rows <- table[as.character(table$feature_id) == id, , drop = FALSE]
  if (nrow(rows) != 1L) return(sprintf("%d rows of id %s", nrow(rows), id))
  c(row = rows$.row, name = rows$name)
}

# The declared step on its own owner in its own store root (KIT-214's `k214Step`).
k258Step <- function(object, variant = c("whole", "partitioned")) {
  variant <- match.arg(variant)
  binding <- k258Binding(object)
  root <- withr::local_tempdir()
  fixture <- if (k258Lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
  withr::defer(closeDIAStateAdapter(fixture$owner))
  step <- omicsDuplicateResolutionStep(binding)
  kit172Quiet(runWorkflowOmicsReferenceOperation(fixture$owner, workflowOperationEntry(step$operation)$declaration,
    c(step$parameters, list(variant = variant)), step$state_name, fixture$workflow$config_list, step$description))
  list(object = suppressWarnings(fixture$owner$getState(step$state_name)),
    variant = fixture$owner$.context$snapshot$step_provenance[[step$state_name]]$variant,
    record = buildWorkflowAnalysisRecord(fixture$owner), state = step$state_name)
}

# A rule marker of a recorded step: its source and decoded value, NULL when the step does not list it.
k258Marker <- function(record, state, marker) {
  step <- Filter(function(s) identical(s$state, state), record$steps)
  if (!length(step)) return(NULL)
  rows <- Filter(function(p) identical(p$name, marker), step[[1L]]$parameters)
  if (!length(rows)) return(NULL)
  list(source = rows[[1L]]$source, value = decodeWorkflowRecordValue(rows[[1L]]$value))
}

# The two markers as a step computed now records them (KIT-214's samples marker, this ticket's rank marker), and the
# option under which a step computes with KIT-214's rule.
k258Design0 <- list(source = "fixed", value = "design")
k258RankNow <- list(source = "fixed", value = "measured_then_mean")
k258KIT214 <- list(MultiScholaR.duplicate_resolution_rank = "mean")

# ---- E1 to E3, G1, G3: the declared step, both omics, both variants ---------------------------------------------------

for (case in names(k258Cases)) {
  for (omic in c("lipidomics", "metabolomics")) {
    for (variant in c("whole", "partitioned")) {
      local({
        case <- case
        omic <- omic
        variant <- variant
        spec <- k258Cases[[case]]
        test_that(sprintf("KIT-258 %s (%s, %s): duplicate resolution keeps row %s", case, omic, variant, spec$head), {
          skip_if_not(workflowOwnerCapable())
          table <- spec$table()
          object <- k258WithRow(k258Build(omic, table))
          run <- k258Step(object, variant)
          expect_identical(run$variant, variant)
          expected_row <- as.character(which(table$name == spec$head))
          expect_identical(k258Kept(run$object, spec$id), c(row = expected_row, name = spec$head))
          # One row per id.
          expect_identical(nrow(k258Assays(run$object)[[1L]]), length(unique(table$feature_id)))
        })
      })
    }
  }
}

# ---- E4: the exported kernels and the S4 methods on E1's table (scripted route) -------------------------------------

test_that("KIT-258 E4: the exported kernels keep the most completely measured row", {
  table <- k258E1Table()
  for (kernel in c("resolveDuplicateFeaturesByIntensity", "resolveLipidDuplicateFeaturesByIntensity")) {
    resolved <- kit172Quiet(get(kernel, asNamespace("MultiScholaR"))(table, "feature_id", k258Samples))
    expect_identical(k258KeptName(resolved, "D"), "B", info = kernel)
    expect_identical(k258KeptName(resolved, "E"), "C", info = kernel)
  }
})

for (omic in c("lipidomics", "metabolomics")) {
  local({
    omic <- omic
    test_that(sprintf("KIT-258 E4 (%s): the S4 method resolveDuplicateFeatures keeps the most completely measured row",
      omic), {
      table <- k258E1Table()
      # The S4 methods require their descriptive-name column (`lipid`, `metabolite`); no internal-standard regex.
      table[[if (identical(omic, "lipidomics")) "lipid" else "metabolite"]] <- table$name
      object <- k258Build(omic, table)
      resolved <- kit172Quiet(resolveDuplicateFeatures(object))
      expect_identical(k258KeptName(k258Assays(resolved)[[1L]], "D"), "B")
      expect_identical(k258KeptName(k258Assays(resolved)[[1L]], "E"), "C")
    })
  })
}

# ---- G2: no duplicated id, returned unchanged (passes on base and head) ----------------------------------------------

for (omic in c("lipidomics", "metabolomics")) {
  local({
    omic <- omic
    test_that(sprintf("KIT-258 G2 (%s): a table without a duplicated id is returned unchanged (numeric and text id)",
      omic), {
      text_id <- k258E1Table()[c(1L, 3L), , drop = FALSE]
      numeric_id <- text_id
      numeric_id$feature_id <- c(1000L, 2000L)
      for (table in list(text = text_id, numeric = numeric_id)) {
        object <- k258Build(omic, table)
        before <- k258Assays(object)[[1L]]
        resolved <- kit172Quiet(k258Binding(object)$resolve_duplicates(object)$object)
        expect_identical(k258Assays(resolved)[[1L]], before, info = class(table$feature_id))
      }
    })
  })
}

# ---- T1: records made under KIT-214's rule and now, replayed ----------------------------------------------------------
# As KIT-214's T6 (iv) (`k214Replay`): the owner's record written, the owner closed, a copy cut to the design state and
# the record replayed on it; the replay compared with the origin.

k258Paths <- function(root, lipid) {
  qc <- if (lipid) "lipid_qc_dir" else "metabolite_qc_dir"
  stats::setNames(list(file.path(root, "source"), file.path(root, "qc"), file.path(root, "da"), file.path(root, "pub")),
    c("source_dir", qc, "da_output_dir", "publication_graphs_dir"))
}

# `options` are set while the step is computed (KIT-214's rule: rank "mean"); the step runs as the module runs it.
k258ReplayOwner <- function(root, object, options = list()) {
  fixture <- if (k258Lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
  withr::with_options(options, kit172Quiet(runWorkflowOmicsReferenceDuplicateResolution(k258Binding(object),
    fixture$owner, fixture$workflow$config_list)))
  fixture
}

k258Replay <- function(fixture, root, object, id) {
  prefix <- k258Binding(object)$prefix
  state <- paste0(prefix, "_duplicates_resolved")
  kept <- k258KeptName(k258Assays(suppressWarnings(fixture$owner$getState(state)))[[1L]], id)
  json <- file.path(root, "analysis_record.json")
  writeWorkflowAnalysisRecord(fixture$owner, json)
  record <- jsonlite::read_json(json, simplifyVector = FALSE)
  closeDIAStateAdapter(fixture$owner)
  source_dir <- file.path(root, "source")
  store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE,
    full.names = TRUE)[[1L]]))
  copy <- tempfile("k258-copy-")
  dir.create(copy)
  stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
  harness <- new.env(parent = globalenv())
  sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
  prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(prefix, "_raw_data_s4")), api)
  inputs <- c(stats::setNames(file.path(source_dir, "import_sources", "LCMS_Pos_export.txt"), paste0(prefix, "_assay_1")),
    design_file = file.path(source_dir, "design_matrix.tab"))
  replay <- kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, k258Paths(root, k258Lipid(object)),
    tempfile("k258-scratch-"), api))
  identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
  replayed <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
  replayed_kept <- k258KeptName(k258Assays(suppressWarnings(replayed$getState(state)))[[1L]], id)
  replayed_record <- buildWorkflowAnalysisRecord(replayed)
  closeDIAStateAdapter(replayed)
  list(record = record, replayed_record = replayed_record, mismatches = as.character(replay$mismatches), kept = kept,
    replayed_kept = replayed_kept, state = state)
}

for (omic in c("lipidomics", "metabolomics")) {
  local({
    omic <- omic
    test_that(sprintf("KIT-258 T1 (a) (%s): a record made under KIT-214's rule replays its row, without the rank marker",
      omic), {
      skip_if_not(workflowOwnerCapable())
      object <- k258Build(omic, k258E1Table())
      root <- withr::local_tempdir()
      a <- k258Replay(k258ReplayOwner(root, object, k258KIT214), root, object, "D")
      expect_identical(a$kept, "A")
      expect_identical(a$replayed_kept, "A")
      expect_identical(a$mismatches, character())
      expect_null(k258Marker(a$record, a$state, "duplicate_resolution_rank"))
      expect_null(k258Marker(a$replayed_record, a$state, "duplicate_resolution_rank"))
      expect_identical(k258Marker(a$record, a$state, "duplicate_resolution_samples"), k258Design0)
    })

    test_that(sprintf("KIT-258 T1 (b) (%s): a record made now lists the rank marker and replays row B", omic), {
      skip_if_not(workflowOwnerCapable())
      object <- k258Build(omic, k258E1Table())
      root <- withr::local_tempdir()
      b <- k258Replay(k258ReplayOwner(root, object), root, object, "D")
      expect_identical(b$kept, "B")
      expect_identical(b$replayed_kept, "B")
      expect_identical(b$mismatches, character())
      expect_identical(k258Marker(b$record, b$state, "duplicate_resolution_rank"), k258RankNow)
      expect_identical(k258Marker(b$replayed_record, b$state, "duplicate_resolution_rank"), k258RankNow)
      expect_identical(k258Marker(b$record, b$state, "duplicate_resolution_samples"), k258Design0)
    })
  })
}

# ---- T2: a managed step records both markers (as KIT-214's T6 (v)) ---------------------------------------------------

test_that("KIT-258 T2: a managed duplicate step records the samples marker and the rank marker", {
  skip_if_not(workflowOwnerCapable())
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
    assertWorkflowController = function(...) invisible(TRUE), .package = "MultiScholaR")
  for (lipid in c(TRUE, FALSE)) {
    binding <- if (lipid) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
    root <- withr::local_tempdir()
    fixture <- if (lipid) lipidDeclarationsOwner(root) else metabDeclarationsOwner(root)
    withr::defer(closeDIAStateAdapter(fixture$owner))
    step <- omicsDuplicateResolutionStep(binding)
    worker <- omicsWorkerRun(fixture, step$operation, step$parameters)
    expect_identical(worker$status, "complete", info = paste(binding$prefix, worker$error))
    kit172Quiet(completeWorkflowAction(list(manager = fixture$owner, accepting = TRUE), step$operation, worker$result,
      context = list(state_name = step$state_name, config = list(), description = step$description),
      values = step$parameters))
    record <- buildWorkflowAnalysisRecord(fixture$owner)
    expect_identical(k258Marker(record, step$state_name, "duplicate_resolution_samples"), k258Design0,
      info = binding$prefix)
    expect_identical(k258Marker(record, step$state_name, "duplicate_resolution_rank"), k258RankNow,
      info = binding$prefix)
  }
})
