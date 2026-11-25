# KIT-217 (VAL-173 (2); the owner's decisions of 2026-11-25): the lipidomics and metabolomics intensity filter counts
# present (finite) values per group. A feature is kept when at least g groups each hold at least min(n_i, m) values
# and, as before, its share of samples below the percentile threshold is under the proportion. Defaults m = 2, g = 1;
# groups named like "pool" are not counted; a design without a group column, and a design whose groups are all pool
# groups, are refused; records made before replay with the rule they were made under (marker
# `intensity_filter_missing`). Tests T1 to T11 and guard G-C1 of the ticket's section 5.6 and 5.7 (T-M1 to T-M3 edit
# test-metab-qc-shared.R, test-metab-02t-metabolite-intensity-filtering-characterization.R and
# test-metab-s4-qc-shared.R). Every kept set below was derived by hand in the ticket (2.4, 5.6). Values: `+` = 1000,
# `-` = 10, `.` = missing. On the commit before the fix each test fails at the value the ticket states.

# ---- Fixtures --------------------------------------------------------------------------------------------------------

k217Value <- c(`+` = 1000, `-` = 10, `.` = NA_real_)

# A feature table from patterns: one string per feature, one character per sample (`|` separates groups, for reading).
k217Table <- function(patterns, samples) {
  rows <- lapply(patterns, function(p) unname(k217Value[strsplit(gsub("|", "", p, fixed = TRUE), "")[[1L]]]))
  values <- do.call(rbind, unname(rows))
  stopifnot(ncol(values) == length(samples))
  colnames(values) <- samples
  cbind(data.frame(feature = names(patterns), annotation = paste0("a_", names(patterns)), stringsAsFactors = FALSE),
    as.data.frame(values, check.names = FALSE))
}

# The ten-feature table of ticket 2.4 (2 x 6).
k217Ten <- c(F1 = "++++++|++++++", F2 = "------|++++++", F3 = "......|......", F4 = "++++++|......",
  F5 = "+++...|+++...", F6 = "+.....|+.....", F7 = "++....|++++++", F8 = "++....|......", F9 = "---...|+++...",
  F10 = "-----+|-----+")
k217Samples <- c(paste0("A", 1:6), paste0("B", 1:6))
k217Groups <- rep(c("A", "B"), each = 6L)

k217Design <- function(samples, groups = NULL) {
  design <- data.frame(Run = samples, stringsAsFactors = FALSE)
  if (!is.null(groups)) design$group <- groups
  design$tech_rep_group <- samples
  design
}

k217Object <- function(omic, table, design) {
  if (identical(omic, "lipid")) {
    methods::new("LipidomicsAssayData", lipid_data = list(LCMS_Pos = table), design_matrix = design,
      lipid_id_column = "feature", annotation_id_column = "annotation", sample_id = "Run", group_id = "group",
      technical_replicate_id = "tech_rep_group", database_identifier_type = "Unknown", internal_standard_regex = "^IS_",
      args = list())
  } else {
    methods::new("MetaboliteAssayData", metabolite_data = list(LCMS_Pos = table), design_matrix = design,
      metabolite_id_column = "feature", annotation_id_column = "annotation", sample_id = "Run", group_id = "group",
      technical_replicate_id = "tech_rep_group", database_identifier_type = "Unknown", internal_standard_regex = "^IS_",
      args = list())
  }
}

k217TenObject <- function(omic) k217Object(omic, k217Table(k217Ten, k217Samples), k217Design(k217Samples, k217Groups))

k217Lipid <- function(object) methods::is(object, "LipidomicsAssayData")
k217Binding <- function(object) if (k217Lipid(object)) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
# The two new declared names, written out (the commit before the fix does not declare them).
k217Names <- function(object) {
  if (k217Lipid(object)) {
    list(percentile = "lipids_intensity_cutoff_percentile", proportion = "lipids_proportion_of_samples_below_cutoff",
      m = "lipids_min_values_per_group", g = "lipids_min_groups")
  } else {
    list(percentile = "metabolites_intensity_cutoff_percentile", proportion = "metabolites_proportion_of_samples_below_cutoff",
      m = "metabolites_min_values_per_group", g = "metabolites_min_groups")
  }
}
k217Kept <- function(object) {
  as.character(methods::slot(object, k217Binding(object)$data_slot)[[1L]][["feature"]])
}

# The exported helper of one omic on a table with the id column "feature".
k217Helper <- function(omic, table, threshold, proportion, ...) {
  if (identical(omic, "lipid")) {
    lipidIntensityFilteringHelper(assay_table = table, min_lipid_intensity_threshold = threshold,
      lipids_proportion_of_samples_below_cutoff = proportion, lipid_id_column = "feature", ...)
  } else {
    metaboliteIntensityFilteringHelper(assay_table = table, min_metabolite_intensity_threshold = threshold,
      metabolites_proportion_of_samples_below_cutoff = proportion, metabolite_id_column = "feature", ...)
  }
}

# On the commit before the fix a call with an argument it does not have errors; the test then fails at is.null().
k217Try <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) NULL)

# The S4 method by its generic's name; `m` and `g` are passed only when given.
k217Method <- function(object, percentile, proportion, m = NULL, g = NULL) {
  names <- k217Names(object)
  args <- list(object)
  args[[names$percentile]] <- percentile
  args[[names$proportion]] <- proportion
  if (!is.null(m)) args[[names$m]] <- m
  if (!is.null(g)) args[[names$g]] <- g
  suppressWarnings(suppressMessages(do.call(if (k217Lipid(object)) "lipidIntensityFiltering" else
    "metaboliteIntensityFiltering", args)))
}

k217Config <- list(qc = list(mode = "kit217"))

# The declared step through the wrapper the modules use, with the step built by omicsIntensityFilterStep and the two
# settings set by name when given (the commit before the fix does not declare them, so it drops them).
k217RunStep <- function(owner, object, percentile, proportion, variant = "partitioned", m = NULL, g = NULL,
                        config = k217Config) {
  names <- k217Names(object)
  step <- omicsIntensityFilterStep(k217Binding(object), percentile, proportion)
  if (!is.null(m)) step$parameters[[names$m]] <- m
  if (!is.null(g)) step$parameters[[names$g]] <- g
  step$parameters$variant <- variant
  kit172Quiet(runWorkflowOmicsReferenceStep(owner, step, config))
}

# A bare owner holding the object (refusals and kept sets; KIT-202 T5's construction).
k217Owner <- function(object, env = parent.frame()) {
  omicsTestOwner(object, "raw", config = k217Config, env = env)$state_manager
}

# The declared step on its own owner; returns the kept ids of the filtered state, the step's value and the record.
k217Step <- function(object, percentile, proportion, variant = "partitioned", m = NULL, g = NULL) {
  owner <- k217Owner(object)
  out <- k217RunStep(owner, object, percentile, proportion, variant, m, g)
  state <- paste0(k217Binding(object)$prefix, "_intensity_filtered")
  list(kept = k217Kept(owner$getState(state)), value = out$value, record = buildWorkflowAnalysisRecord(owner),
    state = state, variant = owner$.context$snapshot$step_provenance[[state]]$variant)
}

k217Param <- function(record, state, name) {
  step <- Filter(function(s) identical(s$state, state), record$steps)
  if (!length(step)) return(NULL)
  rows <- Filter(function(p) identical(p$name, name), step[[1L]]$parameters)
  if (!length(rows)) return(NULL)
  list(source = rows[[1L]]$source, value = decodeWorkflowRecordValue(rows[[1L]]$value))
}

k217Earlier <- list(MultiScholaR.intensity_filter_missing = "not_counted")

# The kept sets of ticket 5.6 (T2 on the commit before the fix: today's rule).
k217HeadDefault <- c("F1", "F4", "F5", "F7", "F8", "F9")
k217HeadG2 <- c("F1", "F5", "F7", "F9")
k217Base <- c("F1", "F3", "F4", "F5", "F6", "F7", "F8", "F9")

# ---- T1: the exported helpers ----------------------------------------------------------------------------------------

for (omic in c("lipid", "metab")) {
  test_that(sprintf("KIT-217 T1 (%s): the helper counts present values per group, g = 1 by default and selectable", omic), {
    table <- k217Table(k217Ten, k217Samples)
    call <- function(...) k217Try(k217Helper(omic, table, 100, 0.5, sample_columns = k217Samples,
      group_labels = k217Groups, min_values_per_group = 2L, ...))
    defaults <- k217Try(k217Helper(omic, table, 100, 0.5, sample_columns = k217Samples, group_labels = k217Groups))
    expect_false(is.null(defaults))
    if (!is.null(defaults)) expect_identical(defaults$feature, k217HeadDefault)
    at_two <- call(min_groups = 2L)
    expect_false(is.null(at_two))
    if (!is.null(at_two)) expect_identical(at_two$feature, k217HeadG2)
    at_one <- call(min_groups = 1L)
    expect_false(is.null(at_one))
    if (!is.null(at_one)) expect_identical(at_one$feature, k217HeadDefault)
  })
}

# ---- T2: the S4 methods ----------------------------------------------------------------------------------------------

for (omic in c("lipid", "metab")) {
  test_that(sprintf("KIT-217 T2 (%s): the S4 method keeps by the per-group count at the defaults, g selectable", omic), {
    object <- k217TenObject(omic)
    # tau = 529.75: type 7 over the 66 finite values (19 tens, 47 thousands); h = 65 x 0.285 = 18.525.
    values <- as.matrix(k217Table(k217Ten, k217Samples)[, k217Samples])
    expect_equal(unname(stats::quantile(as.vector(values), 0.285, na.rm = TRUE, type = 7)), 529.75)
    expect_identical(k217Kept(k217Method(object, 28.5, 0.5)), k217HeadDefault)
    at_two <- k217Try(k217Method(object, 28.5, 0.5, g = 2L))
    expect_false(is.null(at_two))
    if (!is.null(at_two)) expect_identical(k217Kept(at_two), k217HeadG2)
  })
}

# ---- T3: four groups, g from 1 to 4 ----------------------------------------------------------------------------------

test_that("KIT-217 T3: four groups of 3, g from 1 to 4 (metabolite helper)", {
  samples <- paste0(rep(c("A", "B", "C", "D"), each = 3L), 1:3)
  groups <- rep(c("A", "B", "C", "D"), each = 3L)
  table <- k217Table(c(G1 = "+++|+++|...|...", G2 = "+++|...|...|...", G3 = "++.|++.|++.|++.",
    G4 = "+..|+..|+..|+.."), samples)
  expected <- list(`1` = c("G1", "G2", "G3"), `2` = c("G1", "G3"), `3` = "G3", `4` = "G3")
  defaults <- k217Try(k217Helper("metab", table, 100, 0.5, sample_columns = samples, group_labels = groups,
    min_values_per_group = 2L))
  expect_false(is.null(defaults))
  if (!is.null(defaults)) expect_identical(defaults$feature, expected[["1"]])
  for (g in 1:4) {
    out <- k217Try(k217Helper("metab", table, 100, 0.5, sample_columns = samples, group_labels = groups,
      min_values_per_group = 2L, min_groups = g))
    expect_false(is.null(out), info = g)
    if (!is.null(out)) expect_identical(out$feature, expected[[as.character(g)]], info = g)
  }
})

# ---- T4: unequal groups, exact per group -----------------------------------------------------------------------------

for (omic in c("lipid", "metab")) {
  test_that(sprintf("KIT-217 T4 (%s): groups of 3 and 6 with m = 4 need 3 and 4 values", omic), {
    samples <- c(paste0("A", 1:3), paste0("B", 1:6))
    groups <- rep(c("A", "B"), c(3L, 6L))
    table <- k217Table(c(U1 = "+++|++++..", U2 = "++.|++++++", U3 = "++.|+++..."), samples)
    for (case in list(list(g = 2L, kept = "U1"), list(g = 1L, kept = c("U1", "U2")))) {
      out <- k217Try(k217Helper(omic, table, 100, 0.5, sample_columns = samples, group_labels = groups,
        min_values_per_group = 4L, min_groups = case$g))
      expect_false(is.null(out), info = case$g)
      if (!is.null(out)) expect_identical(out$feature, case$kept, info = case$g)
    }
  })
}

# ---- T5: the declared step, both variants ----------------------------------------------------------------------------

for (omic in c("lipid", "metab")) {
  for (variant in c("partitioned", "whole")) {
    test_that(sprintf("KIT-217 T5 (%s, %s): the declared step at its defaults keeps by the per-group count", omic,
      variant), {
      skip_if_not(workflowOwnerCapable())
      object <- k217TenObject(omic)
      names <- k217Names(object)
      run <- k217Step(object, 28.5, 0.5, variant)
      expect_identical(run$kept, k217HeadDefault)
      expect_identical(run$variant, variant)
      m <- k217Param(run$record, run$state, names$m)
      g <- k217Param(run$record, run$state, names$g)
      expect_false(is.null(m))
      expect_false(is.null(g))
      if (!is.null(m)) expect_equal(m$value, 2)
      if (!is.null(g)) expect_equal(g$value, 1)
      expect_equal(unname(run$value$removed_missing[["LCMS_Pos"]]), 2)
      expect_equal(unname(run$value$removed_low[["LCMS_Pos"]]), 2)
    })
  }
}

# ---- T6: refusals of the settings ------------------------------------------------------------------------------------

k217Invalid <- list(list(m = 0L), list(m = 1.5), list(g = 0L), list(g = 3L))

for (omic in c("lipid", "metab")) {
  test_that(sprintf("KIT-217 T6 (%s): values per group and groups required out of range are refused", omic), {
    object <- k217TenObject(omic)
    for (case in k217Invalid) {
      info <- paste(names(case), unlist(case))
      expect_error(k217Method(object, 28.5, 0.5, m = case$m, g = case$g),
        class = "omics_intensity_filter_counts_invalid", info = info)
    }
    skip_if_not(workflowOwnerCapable())
    for (case in k217Invalid) {
      info <- paste(names(case), unlist(case))
      owner <- k217Owner(object)
      expect_error(k217RunStep(owner, object, 28.5, 0.5, m = case$m, g = case$g),
        class = "omics_intensity_filter_counts_invalid", info = info)
    }
  })
}

# ---- T7: replay ------------------------------------------------------------------------------------------------------

test_that("KIT-217 T7 (i): the replay rule of the intensity filter follows the step's marker", {
  markers <- get0("workflowReplayRuleMarkers", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  expect_true(!is.null(markers[["intensity_filter_missing"]]))
  if (!is.null(markers[["intensity_filter_missing"]])) {
    entry <- function(source) list(name = "intensity_filter_missing", source = source,
      value = encodeWorkflowRecordValue("group_counts"), provenance = "fixed")
    for (operation in c("lipid_intensity_filter", "metab_intensity_filter")) {
      step <- function(parameters) list(state = "s", operation = operation, parameters = parameters)
      expect_identical(workflowReplayRule(step(list(entry("fixed"))), "intensity_filter_missing"), "group_counts",
        info = operation)
      expect_identical(workflowReplayRule(step(list()), "intensity_filter_missing"), "not_counted", info = operation)
      expect_identical(workflowReplayRule(step(list(entry("registry"))), "intensity_filter_missing"), "not_counted",
        info = operation)
    }
  }
})

for (omic in c("lipid", "metab")) {
  test_that(sprintf("KIT-217 T7 (ii) (%s): under the earlier rule the head gives today's kept set, unrefused", omic), {
    withr::local_options(k217Earlier)
    object <- k217TenObject(omic)
    expect_identical(k217Kept(k217Method(object, 28.5, 0.5)), k217Base)
    no_column <- k217Object(omic, k217Table(k217Ten, k217Samples), k217Design(k217Samples))
    expect_identical(k217Kept(k217Method(no_column, 28.5, 0.5)), k217Base)
    skip_if_not(workflowOwnerCapable())
    for (variant in c("partitioned", "whole")) {
      expect_identical(k217Step(object, 28.5, 0.5, variant)$kept, k217Base, info = variant)
    }
  })

  test_that(sprintf("KIT-217 T7 (iii) (%s): a step made now records the marker", omic), {
    skip_if_not(workflowOwnerCapable())
    run <- k217Step(k217TenObject(omic), 28.5, 0.5)
    expect_identical(k217Param(run$record, run$state, "intensity_filter_missing"),
      list(source = "fixed", value = "group_counts"))
  })
}

# T7 (iv) and (v), T11: KIT-202 T6 (iv)'s construction (tests/testthat/test-kit202-intensity-filter-samples.R,
# k202ReplayOwner and k202Replay): the owner's record written, the owner closed, a copy cut to the design state and the
# record replayed on it; the replay compared with the origin.
k217Paths <- function(root, lipid) {
  qc <- if (lipid) "lipid_qc_dir" else "metabolite_qc_dir"
  stats::setNames(list(file.path(root, "source"), file.path(root, "qc"), file.path(root, "da"), file.path(root, "pub")),
    c("source_dir", qc, "da_output_dir", "publication_graphs_dir"))
}

# One owner in its own store root; `earlier = TRUE` computes the step under the earlier rule, so its provenance has
# no marker, as a build before this change wrote it; `g` sets groups required by name.
k217ReplayOwner <- function(root, object, earlier = FALSE, g = NULL) {
  fixture <- if (k217Lipid(object)) lipidDeclarationsOwner(root, object, k217Config) else
    metabDeclarationsOwner(root, object, k217Config)
  run <- function() k217RunStep(fixture$owner, object, 28.5, 0.5, g = g, config = fixture$workflow$config_list)
  if (earlier) withr::with_options(k217Earlier, run()) else run()
  fixture
}

k217Replay <- function(fixture, root, object) {
  prefix <- k217Binding(object)$prefix
  state <- paste0(prefix, "_intensity_filtered")
  kept <- k217Kept(fixture$owner$getState(state))
  json <- file.path(root, "analysis_record.json")
  writeWorkflowAnalysisRecord(fixture$owner, json)
  record <- jsonlite::read_json(json, simplifyVector = FALSE)
  closeDIAStateAdapter(fixture$owner)
  source_dir <- file.path(root, "source")
  store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE,
    full.names = TRUE)[[1L]]))
  copy <- tempfile("k217-copy-")
  dir.create(copy)
  stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
  harness <- new.env(parent = globalenv())
  sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
  prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(prefix, "_raw_data_s4")), api)
  inputs <- c(stats::setNames(file.path(source_dir, "import_sources", "LCMS_Pos_export.txt"), paste0(prefix, "_assay_1")),
    design_file = file.path(source_dir, "design_matrix.tab"))
  replay <- kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, k217Paths(root, k217Lipid(object)),
    tempfile("k217-scratch-"), api))
  identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
  replayed <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
  replayed_record <- buildWorkflowAnalysisRecord(replayed)
  closeDIAStateAdapter(replayed)
  list(record = record, replayed_record = replayed_record, mismatches = as.character(replay$mismatches), kept = kept,
    state = state)
}

k217Marker <- function(record, state) k217Param(record, state, "intensity_filter_missing")

for (omic in c("lipid", "metab")) {
  test_that(sprintf("KIT-217 T7 (iv) and T11 (%s): records made before and after the change replay, each with its own rule",
    omic), {
    skip_if_not(workflowOwnerCapable())
    object <- k217TenObject(omic)
    names <- k217Names(object)
    root_a <- withr::local_tempdir()
    root_b <- withr::local_tempdir()
    a <- k217Replay(k217ReplayOwner(root_a, object, earlier = TRUE), root_a, object)
    b <- k217Replay(k217ReplayOwner(root_b, object), root_b, object)
    # Non-vacuity: the two rules give different filtered states on this table (F3 and F6).
    expect_identical(setdiff(a$kept, b$kept), c("F3", "F6"))
    expect_identical(a$mismatches, character())
    expect_identical(b$mismatches, character())
    expect_null(k217Marker(a$record, a$state))
    expect_null(k217Marker(a$replayed_record, a$state))
    expect_identical(k217Marker(b$record, b$state), list(source = "fixed", value = "group_counts"))
    expect_identical(k217Marker(b$replayed_record, b$state), list(source = "fixed", value = "group_counts"))
    # T11 (SD-17): the record made under the earlier rule shows the two settings as not applicable, not as missing.
    for (name in c(names$m, names$g)) {
      expect_identical(k217Param(a$record, a$state, name)$value, "not applicable: earlier rule", info = name)
      expect_false(paste0("declared parameter not recorded: ", a$state, "/", name) %in%
        validateWorkflowAnalysisRecord(a$record)$incomplete, info = name)
    }
  })

  test_that(sprintf("KIT-217 T7 (v) (%s): a record made at g = 2 replays at g = 2", omic), {
    skip_if_not(workflowOwnerCapable())
    object <- k217TenObject(omic)
    root <- withr::local_tempdir()
    replay <- k217Replay(k217ReplayOwner(root, object, g = 2L), root, object)
    recorded <- k217Param(replay$record, replay$state, k217Names(object)$g)
    expect_false(is.null(recorded))
    if (!is.null(recorded)) expect_equal(recorded$value, 2)
    expect_identical(replay$kept, k217HeadG2)
    expect_identical(replay$mismatches, character())
  })
}

# T7 (vi): a managed step, as KIT-202 T6 (v): the worker's entry point (omicsWorkerRun), completed with
# completeWorkflowAction.
test_that("KIT-217 T7 (vi): a managed intensity filter step records the marker", {
  skip_if_not(workflowOwnerCapable())
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
    assertWorkflowController = function(...) invisible(TRUE), .package = "MultiScholaR")
  for (lipid in c(TRUE, FALSE)) {
    binding <- if (lipid) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
    root <- withr::local_tempdir()
    fixture <- if (lipid) lipidDeclarationsOwner(root) else metabDeclarationsOwner(root)
    withr::defer(closeDIAStateAdapter(fixture$owner))
    step <- omicsIntensityFilterStep(binding, 10, 0.5)
    worker <- omicsWorkerRun(fixture, step$operation, step$parameters)
    expect_identical(worker$status, "complete", info = paste(binding$prefix, worker$error))
    kit172Quiet(completeWorkflowAction(list(manager = fixture$owner, accepting = TRUE), step$operation, worker$result,
      context = list(state_name = step$state_name, config = list(), description = step$description),
      values = step$parameters))
    expect_identical(k217Marker(buildWorkflowAnalysisRecord(fixture$owner), step$state_name),
      list(source = "fixed", value = "group_counts"), info = binding$prefix)
  }
})

# ---- T8: the helpers without groups ----------------------------------------------------------------------------------

for (omic in c("lipid", "metab")) {
  test_that(sprintf("KIT-217 T8 (%s): without group labels only features with no value are removed, with a warning",
    omic), {
    table <- k217Table(k217Ten, k217Samples)
    expect_warning(out <- k217Helper(omic, table, 100, 0.5, sample_columns = k217Samples),
      class = "omics_intensity_filter_groups_guessed")
    expect_identical(out$feature, c("F1", "F4", "F5", "F6", "F7", "F8", "F9"))
    ignored <- k217Try(k217Helper(omic, table, 100, 0.5, sample_columns = k217Samples, min_values_per_group = 1L,
      min_groups = 2L))
    expect_false(is.null(ignored))
    if (!is.null(ignored)) expect_identical(ignored$feature, c("F1", "F4", "F5", "F6", "F7", "F8", "F9"))
  })
}

# ---- T9: the tab -----------------------------------------------------------------------------------------------------

for (omic in c("lipidomics", "metabolomics")) {
  test_that(sprintf("KIT-217 T9 (%s): the tab's two settings default to 2 and 1 and reach the step", omic), {
    config <- omicGuiConfig(omic)
    html <- as.character(htmltools::renderTags(mod_omic_qc_intensity_ui("x", config))$html)
    input_tag <- function(id) regmatches(html, regexpr(sprintf("<input[^>]*id=\"x-%s\"[^>]*>", id), html))
    for (case in list(c("min_values_per_group", "2"), c("min_groups", "1"))) {
      tag <- input_tag(case[[1L]])
      expect_length(tag, 1L)
      if (length(tag)) expect_match(tag, sprintf("value=\"%s\"", case[[2L]]), fixed = TRUE, info = case[[1L]])
    }
    binding <- config$binding$binding
    captured <- new.env(parent = emptyenv())
    captured$steps <- list()
    testthat::local_mocked_bindings(
      omicQcHeadMetadata = function(...) invisible(NULL),
      assertOmicQcHeadNotNormalised = function(...) invisible(NULL),
      runWorkflowOmicsQcAction = function(binding, workflow_data, step, ...) {
        captured$steps[[length(captured$steps) + 1L]] <- step
        invisible(NULL)
      },
      .package = "MultiScholaR")
    workflow <- shiny::reactiveValues(state_manager = "owner")
    names <- if (identical(omic, "lipidomics")) c("lipids_min_values_per_group", "lipids_min_groups") else
      c("metabolites_min_values_per_group", "metabolites_min_groups")
    suppressMessages(shiny::testServer(mod_omic_qc_intensity_server, args = list(workflow_data = workflow,
      omic_type = omic, experiment_label = "kit217", config = config), {
      session$setInputs(intensity_cutoff_percentile = 10, proportion_below_cutoff = 0.5, min_values_per_group = 2,
        min_groups = 1)
      session$setInputs(apply_filter = 1L)
      session$setInputs(min_groups = 2)
      session$setInputs(apply_filter = 2L)
    }))
    expect_length(captured$steps, 2L)
    if (length(captured$steps) == 2L) {
      expect_equal(captured$steps[[1L]]$parameters[[names[[1L]]]], 2)
      expect_equal(captured$steps[[1L]]$parameters[[names[[2L]]]], 1)
      expect_equal(captured$steps[[2L]]$parameters[[names[[1L]]]], 2)
      expect_equal(captured$steps[[2L]]$parameters[[names[[2L]]]], 2)
    }
  })
}

# ---- T10: groups, pools and the design's refusals --------------------------------------------------------------------

test_that("KIT-217 T10 (i): pool groups are not counted (S4 method; g through the declared step)", {
  samples <- c(paste0("A", 1:3), paste0("B", 1:3), paste0("P", 1:3))
  groups <- rep(c("A", "B", "Pool"), each = 3L)
  table <- k217Table(c(P1 = "...|...|+++", P2 = "+++|...|+++", P3 = "++.|++.|..."), samples)
  object <- k217Object("metab", table, k217Design(samples, groups))
  expect_identical(k217Kept(k217Method(object, 10, 0.5)), c("P2", "P3"))
  skip_if_not(workflowOwnerCapable())
  expect_identical(k217Step(object, 10, 0.5, g = 2L)$kept, "P3")
  expect_error(k217RunStep(k217Owner(object), object, 10, 0.5, g = 3L), class = "omics_intensity_filter_counts_invalid")
})

for (omic in c("lipid", "metab")) {
  test_that(sprintf("KIT-217 T10 (ii) (%s): no group column, or only pool groups, is refused", omic), {
    no_column <- k217Object(omic, k217Table(k217Ten, k217Samples), k217Design(k217Samples))
    samples <- c(paste0("P", 1:3), paste0("Q", 1:3))
    pools <- k217Object(omic, k217Table(c(P1 = "+++|+++", P2 = "++.|..."), samples),
      k217Design(samples, rep(c("Pool", "QC_pool"), each = 3L)))
    expect_error(k217Method(no_column, 28.5, 0.5), class = "omics_intensity_filter_counts_invalid")
    expect_error(k217Method(pools, 10, 0.5), class = "omics_intensity_filter_counts_invalid")
    message <- tryCatch({ k217Method(pools, 10, 0.5); "" }, error = function(e) conditionMessage(e))
    expect_match(message, "Give at least one group of samples a name without", fixed = TRUE)
    skip_if_not(workflowOwnerCapable())
    expect_error(k217RunStep(k217Owner(no_column), no_column, 28.5, 0.5),
      class = "omics_intensity_filter_counts_invalid")
    expect_error(k217RunStep(k217Owner(pools), pools, 10, 0.5), class = "omics_intensity_filter_counts_invalid")
  })

  test_that(sprintf("KIT-217 T10 (iii) (%s): an object with no assays keeps its warning before any group check", omic), {
    design <- data.frame(Run = "S1", stringsAsFactors = FALSE)
    object <- if (identical(omic, "lipid")) methods::new("LipidomicsAssayData", lipid_data = list(), design_matrix = design,
      sample_id = "Run", lipid_id_column = "feature") else methods::new("MetaboliteAssayData", metabolite_data = list(),
      design_matrix = design, sample_id = "Run", metabolite_id_column = "feature")
    names <- k217Names(object)
    args <- stats::setNames(list(object, 50, 0.6), c("theObject", names$percentile, names$proportion))
    expect_warning(suppressMessages(do.call(if (identical(omic, "lipid")) "lipidIntensityFiltering" else
      "metaboliteIntensityFiltering", args)), "no assays")
  })

  test_that(sprintf("KIT-217 T10 (iv) (%s): a design with one counted group is filtered at g = 1, refused at g = 2", omic), {
    designs <- list(
      one = list(samples = paste0("A", 1:3), groups = rep("A", 3L),
        patterns = c(Q1 = "+++", Q2 = "+..", Q3 = "...")),
      with_pool = list(samples = c(paste0("A", 1:3), paste0("P", 1:3)), groups = rep(c("A", "Pool"), each = 3L),
        patterns = c(Q1 = "+++|...", Q2 = "+..|+++", Q3 = "...|+++")))
    for (name in names(designs)) {
      d <- designs[[name]]
      object <- k217Object(omic, k217Table(d$patterns, d$samples), k217Design(d$samples, d$groups))
      expect_identical(k217Kept(k217Method(object, 10, 0.5)), "Q1", info = name)
      expect_error(k217Method(object, 10, 0.5, g = 2L), class = "omics_intensity_filter_counts_invalid", info = name)
    }
    skip_if_not(workflowOwnerCapable())
    for (name in names(designs)) {
      d <- designs[[name]]
      object <- k217Object(omic, k217Table(d$patterns, d$samples), k217Design(d$samples, d$groups))
      expect_identical(k217Step(object, 10, 0.5)$kept, "Q1", info = name)
      expect_error(k217RunStep(k217Owner(object), object, 10, 0.5, g = 2L),
        class = "omics_intensity_filter_counts_invalid", info = name)
    }
  })

  test_that(sprintf("KIT-217 T10 (v) (%s): a sample with a missing or blank group is in no counted group, not refused",
    omic), {
    samples <- c(paste0("A", 1:3), paste0("B", 1:3), "s7", "s8")
    groups <- c(rep(c("A", "B"), each = 3L), NA, "")
    table <- k217Table(c(R1 = "...|...|++", R2 = "++.|...|++"), samples)
    object <- k217Object(omic, table, k217Design(samples, groups))
    expect_identical(k217Kept(k217Method(object, 10, 0.5)), "R2")
    out <- k217Try(k217Helper(omic, table, 1000, 0.5, sample_columns = samples, group_labels = groups))
    expect_false(is.null(out))
    if (!is.null(out)) expect_identical(out$feature, "R2")
  })
}

# ---- Guard G-C1: complete data unchanged (passes before and after the fix) -------------------------------------------

# The application's import path, copied from KIT-202's tests (test-kit202-intensity-filter-samples.R at 185ebf75,
# k202RunInput and k202Lane, themselves from KIT-198's G0 script); here it returns the design object.
k217RunInput <- function(omic, vendor_format, assay1_file, assay1_name, assay2_file = NULL, assay2_name = NULL,
                         custom = NULL, sanitize = TRUE, design = NULL) {
  config <- omicGuiConfig(omic)
  import <- config$import
  preview <- config$import$preview(assay1_file, import = import)
  headers <- names(preview$assay_data)
  if (identical(vendor_format, "custom")) {
    id_col <- resolveOmicsImportColumn(headers, vendor_format, NULL, custom$id)
    annotation_col <- resolveOmicsImportColumn(headers, vendor_format, NULL, custom$annotation)
    sample_pattern <- custom$samples
    is_pattern <- custom$is_pattern
  } else {
    id_col <- preview$updates$id$selected
    annotation_col <- preview$updates$annotation$selected
    sample_pattern <- ""
    is_pattern <- preview$updates$is_pattern %||% ""
  }
  exclude <- if (identical(omic, "lipidomics")) OMICS_IMPORT_NORMALISED_COLUMNS
  sample_columns <- resolveOmicsImportSampleColumns(preview$assay_data, vendor_format, sample_pattern,
    preview$import_result, exclude)
  source_dir <- file.path(tempfile("kit217-import-"), "source")
  dir.create(source_dir, recursive = TRUE)
  request <- list(assay1_name = assay1_name, assay1_data = preview$assay_data, assay1_file = assay1_file,
    assay2_file = assay2_file, assay2_name = assay2_name, vendor_format = vendor_format,
    detected_format = preview$detected_format, id_col = id_col, annotation_col = annotation_col,
    sample_columns = sample_columns, sanitize_names = sanitize, is_pattern = is_pattern,
    experiment_paths = list(source_dir = source_dir))
  imported <- importOmicsTables(request, config)
  data_cln <- buildOmicsDesignSaveDataList(imported$assays, design$sample, imported$column_mapping, import$id_mapping_key)
  design_matrix <- data.frame(Run = design$sample, group = design$group,
    tech_rep_group = paste(design$group, ave(design$sample, design$group, FUN = seq_along), sep = "_"),
    stringsAsFactors = FALSE)
  config$design$build_object(list(data_cln = data_cln, design_matrix = design_matrix, config_list = list()),
    imported$column_mapping)
}

k217Lane <- function(lane_id) {
  lane <- guiEntryLane(lane_id)
  dir <- normalizePath(testthat::test_path("..", "testdata", "e2e", lane$fixture_dir), mustWork = TRUE)
  assays <- unlist(lane$assays)
  design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  second <- !is.null(lane$assay2_file) && length(assays) >= 2L
  custom <- list(id = "Feature.Name", annotation = "", samples = "^(WT|KO)_", is_pattern = "^ITSD|_d[0-9]+$")
  suppressWarnings(suppressMessages(k217RunInput(lane$omic_type, lane$import_tool, file.path(dir, lane$seed_file),
    assays[[1L]], assay2_file = if (second) file.path(dir, lane$assay2_file), assay2_name = if (second) assays[[2L]],
    custom = if (identical(lane$import_tool, "custom")) custom, sanitize = FALSE, design = design)))
}

# Today's rule written from the ticket (L-1 to L-3), not from package code: the type-7 percentile of the design
# samples' values of the assay; kept when the share of values strictly below it is under the proportion.
k217Oracle <- function(table, samples, percentile, proportion) {
  values <- as.matrix(as.data.frame(table)[, samples, drop = FALSE])
  threshold <- unname(stats::quantile(as.vector(values), probs = percentile / 100, na.rm = TRUE, type = 7))
  rowSums(values < threshold, na.rm = TRUE) / length(samples) < proportion
}

# Groups required, passed only where the method has the formal (after the fix); before it the call is today's.
k217HasCounts <- function(object) {
  k217Names(object)$g %in% names(formals(if (k217Lipid(object)) lipidIntensityFiltering else metaboliteIntensityFiltering))
}

# G-C1's objects: the four GUI-entry lanes (complete tables, ticket 5.11) at g = 1 and g = 2; the module-CI QC objects
# at g = 1 only (each of their assays holds one missing value, L_missing / M_missing in group WT, which a per-group
# count at g = 2 removes by design: they are not complete inputs; recorded in the outcome).
k217Complete <- list(lipid_canonical = 1:2, metab_lc = 1:2, metab_gc = 1:2, metab_combined = 1:2,
  module_ci_lipid = 1L, module_ci_metab = 1L)

for (name in names(k217Complete)) {
  test_that(sprintf("KIT-217 G-C1 (%s): complete data keep today's features at every setting", name), {
    object <- switch(name, module_ci_lipid = module_ci_lipid_qc_object(), module_ci_metab = module_ci_metab_qc_object(),
      k217Lane(name))
    binding <- k217Binding(object)
    id <- methods::slot(object, binding$id_slot)
    samples <- as.character(object@design_matrix[[object@sample_id]])
    for (percentile in c(1, 10)) for (proportion in c(0.2, 0.5)) for (g in k217Complete[[name]]) {
      info <- sprintf("%s at (%s, %s, g = %s)", name, percentile, proportion, g)
      shipped <- k217Method(object, percentile, proportion, g = if (k217HasCounts(object)) g)
      earlier <- withr::with_options(k217Earlier, k217Method(object, percentile, proportion))
      for (assay in names(methods::slot(object, binding$data_slot))) {
        table <- methods::slot(object, binding$data_slot)[[assay]]
        keep <- k217Oracle(table, intersect(samples, names(table)), percentile, proportion)
        filtered <- methods::slot(shipped, binding$data_slot)[[assay]]
        expect_identical(as.character(filtered[[id]]), as.character(table[[id]][keep]), info = paste(info, assay))
        expect_identical(filtered, methods::slot(earlier, binding$data_slot)[[assay]], info = paste(info, assay))
      }
    }
  })
}
