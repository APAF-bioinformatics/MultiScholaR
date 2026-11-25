# KIT-202 (VAL-118 (1), VAL-130): the lipidomics and metabolomics intensity filter counts only the design's samples,
# in the percentile that sets the threshold and in the share of samples below it; the id column is never counted.
# Tests T1 to T7 of the ticket's section 9 (T8 edits tests/testthat/test-metab-qc-shared.R). Values not shown by
# KIT-198's gate G0 (run 45314460795) were computed by hand from the fixtures (ticket section 9). Quantiles are R's
# default type 7.

# ---- Helper section: the application's import path ------------------------------------------------------------------
# Copied from KIT-198's G0 script, docs/qa/workflow-execution/kit/kit-198/g0-numeric-columns.R at commit b0e913b8
# (branch kit-198-intensity-filter-wording), function `run_input`: the configuration's preview, the panel's column
# resolvers, `importOmicsTables`, `buildOmicsDesignSaveDataList` with the input's design, then the configuration's
# `design$build_object`. Here it returns the design object.
k202RunInput <- function(omic, vendor_format, assay1_file, assay1_name, assay2_file = NULL, assay2_name = NULL,
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
  source_dir <- file.path(tempfile("kit202-import-"), "source")
  dir.create(source_dir, recursive = TRUE)
  request <- list(assay1_name = assay1_name, assay1_data = preview$assay_data, assay1_file = assay1_file,
    assay2_file = assay2_file, assay2_name = assay2_name, vendor_format = vendor_format,
    detected_format = preview$detected_format, id_col = id_col, annotation_col = annotation_col,
    sample_columns = sample_columns, sanitize_names = sanitize, is_pattern = is_pattern,
    experiment_paths = list(source_dir = source_dir))
  imported <- importOmicsTables(request, config)
  if (is.null(design)) {
    id_and_annotation <- c(imported$column_mapping[[import$id_mapping_key]], imported$column_mapping$annotation_col)
    samples <- setdiff(imported$sample_columns, id_and_annotation)
    design <- data.frame(sample = samples, group = sub("_.*$", "", samples), stringsAsFactors = FALSE)
  }
  data_cln <- buildOmicsDesignSaveDataList(imported$assays, design$sample, imported$column_mapping, import$id_mapping_key)
  design_matrix <- data.frame(Run = design$sample, group = design$group,
    tech_rep_group = paste(design$group, ave(design$sample, design$group, FUN = seq_along), sep = "_"),
    stringsAsFactors = FALSE)
  config$design$build_object(list(data_cln = data_cln, design_matrix = design_matrix, config_list = list()),
    imported$column_mapping)
}

k202E2e <- function(...) normalizePath(testthat::test_path("..", "testdata", "e2e", ...), mustWork = TRUE)

# G0's MS-DIAL inputs: the metabolite seed through the metabolomics MS-DIAL import; module_ci_lipid_msdial_assay() with
# `Peak ID` = 1:3 through the lipid MS-DIAL import. Panel defaults; every resolved sample assigned.
k202MetabMsdial <- function() {
  k202RunInput("metabolomics", "msdial", k202E2e("metab_lc", "seed_lcms_pos.csv"), "LCMS_Pos")
}

k202LipidMsdial <- function() {
  assay <- module_ci_lipid_msdial_assay()
  assay[["Peak ID"]] <- seq_len(nrow(assay))
  path <- file.path(tempfile("kit202-lipid-"), "msdial_lipid_pos.csv")
  dir.create(dirname(path), recursive = TRUE)
  utils::write.csv(assay, path, row.names = FALSE)
  k202RunInput("lipidomics", "msdial", path, "LCMS_Pos")
}

# G0's GUI-entry lanes with their manifest import settings and design (custom lanes: id Feature.Name, no annotation,
# samples "^(WT|KO)_", internal standards "^ITSD|_d[0-9]+$"; sanitising off).
k202Lane <- function(lane_id) {
  lane <- guiEntryLane(lane_id)
  dir <- k202E2e(lane$fixture_dir)
  assays <- unlist(lane$assays)
  design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  second <- !is.null(lane$assay2_file) && length(assays) >= 2L
  custom <- list(id = "Feature.Name", annotation = "", samples = "^(WT|KO)_", is_pattern = "^ITSD|_d[0-9]+$")
  k202RunInput(lane$omic_type, lane$import_tool, file.path(dir, lane$seed_file), assays[[1L]],
    assay2_file = if (second) file.path(dir, lane$assay2_file), assay2_name = if (second) assays[[2L]],
    custom = if (identical(lane$import_tool, "custom")) custom, sanitize = FALSE, design = design)
}

# T7: the metabolomics workbook's shape: createMetaboliteAssayData (metabolomics_workflow_starter.rmd:245) on the seed
# with every column kept (`Alignment ID`, `Average Rt(min)`, `Average Mz` numeric) and the 4-sample design.
k202WorkbookObject <- function() {
  seed <- utils::read.csv(k202E2e("metab_lc", "seed_lcms_pos.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  design <- data.frame(Run = c("ctrl_1", "ctrl_2", "treat_1", "treat_2"), group = c("ctrl", "ctrl", "treat", "treat"),
    stringsAsFactors = FALSE)
  createMetaboliteAssayData(metabolite_data = list(LCMS_Pos = seed), design_matrix = design,
    metabolite_id_column = "Alignment ID", annotation_id_column = "Name", sample_id = "Run", group_id = "group")
}

k202Cache <- new.env(parent = emptyenv())
k202Object <- function(name) {
  if (is.null(k202Cache[[name]])) {
    k202Cache[[name]] <- suppressWarnings(suppressMessages(switch(name,
      metab_msdial = k202MetabMsdial(),
      lipid_msdial = k202LipidMsdial(),
      workbook = k202WorkbookObject(),
      module_ci_lipid = module_ci_lipid_qc_object(),
      module_ci_metab = module_ci_metab_qc_object(),
      k202Lane(name))))
  }
  k202Cache[[name]]
}

# ---- Helper section: reading an object, running a route -------------------------------------------------------------

k202Lipid <- function(object) methods::is(object, "LipidomicsAssayData")
k202Binding <- function(object) if (k202Lipid(object)) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
k202Assays <- function(object) methods::slot(object, k202Binding(object)$data_slot)
k202IdColumn <- function(object) methods::slot(object, k202Binding(object)$id_slot)
k202Samples <- function(object) as.character(object@design_matrix[[object@sample_id]])

# The kept ids of a one-assay object (the MS-DIAL objects: integer ids), sorted.
k202Kept <- function(object) sort(as.integer(k202Assays(object)[[1L]][[k202IdColumn(object)]]))

# The S4 method with its messages kept: list(object, thresholds by assay as the method's message prints them).
k202Method <- function(object, percentile, proportion) {
  messages <- character()
  result <- withCallingHandlers(suppressWarnings(if (k202Lipid(object)) {
    lipidIntensityFiltering(object, lipids_intensity_cutoff_percentile = percentile,
      lipids_proportion_of_samples_below_cutoff = proportion)
  } else {
    metaboliteIntensityFiltering(object, metabolites_intensity_cutoff_percentile = percentile,
      metabolites_proportion_of_samples_below_cutoff = proportion)
  }), message = function(m) {
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  lines <- grep("^Calculated minimum intensity threshold for assay '", messages, value = TRUE)
  thresholds <- stats::setNames(as.numeric(sub(".*': ", "", trimws(lines))),
    sub("^Calculated minimum intensity threshold for assay '([^']*)'.*$", "\\1", trimws(lines)))
  list(object = result, thresholds = thresholds)
}

# The rule written from the ticket (2.3, 5.3), not from package code: the type-7 percentile of the design samples'
# values of the assay (missing removed); a feature's share is its values strictly below it (missing not counted) over
# the number of design samples; kept when the share is below the proportion.
k202Oracle <- function(table, samples, percentile, proportion) {
  values <- as.matrix(as.data.frame(table)[, samples, drop = FALSE])
  threshold <- unname(stats::quantile(as.vector(values), probs = percentile / 100, na.rm = TRUE, type = 7))
  share <- rowSums(values < threshold, na.rm = TRUE) / length(samples)
  list(threshold = threshold, keep = share < proportion)
}

# The declared step on its own owner in its own store root (the identity holds neither the variant nor the rule marker:
# a second run on one owner returns the first run's saved product, ticket section 9 T2). Returns the kept ids, the
# variant the provenance names and the analysis record; the owner is closed and the root removed on return.
k202Step <- function(object, percentile, proportion, variant = c("partitioned", "whole")) {
  variant <- match.arg(variant)
  binding <- k202Binding(object)
  root <- withr::local_tempdir()
  fixture <- if (k202Lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
  withr::defer(closeDIAStateAdapter(fixture$owner))
  if (identical(variant, "whole")) kit089WholeIntensity(env = environment())
  kit172Quiet(runWorkflowOmicsReferenceIntensityFilter(binding, fixture$owner, percentile, proportion,
    fixture$workflow$config_list))
  state <- paste0(binding$prefix, "_intensity_filtered")
  list(kept = k202Kept(fixture$owner$getState(state)),
    variant = fixture$owner$.context$snapshot$step_provenance[[state]]$variant,
    record = buildWorkflowAnalysisRecord(fixture$owner), state = state)
}

k202Marker <- function(record, state) {
  step <- Filter(function(s) identical(s$state, state), record$steps)
  if (!length(step)) return(NULL)
  rows <- Filter(function(p) identical(p$name, "intensity_filter_samples"), step[[1L]]$parameters)
  if (!length(rows)) return(NULL)
  list(source = rows[[1L]]$source, value = decodeWorkflowRecordValue(rows[[1L]]$value))
}

k202Earlier <- list(MultiScholaR.intensity_filter_samples = "numeric_columns")

# T1's settings and values. Head: the design's samples; base: every numeric column (G0 shown; lipid kept set
# hand-computed, ticket section 9 T1 (b)).
k202T1 <- list(
  metab_msdial = list(proportion = 0.5, head = list(threshold = 16410, kept = c(1L, 2L, 3L, 4L, 7L, 8L)),
    base = list(threshold = 4.9, kept = 1:8)),
  lipid_msdial = list(proportion = 0.2, head = list(threshold = 100.9, kept = c(2L, 3L)),
    base = list(threshold = 2.2, kept = 3L)))

# ---- T1: the S4 methods on the application's import path -------------------------------------------------------------

for (name in names(k202T1)) {
  test_that(sprintf("KIT-202 T1 (%s): the S4 method counts only the design's samples", name), {
    expected <- k202T1[[name]]
    object <- k202Object(name)
    out <- k202Method(object, 10, expected$proportion)
    expect_equal(unname(out$thresholds[["LCMS_Pos"]]), expected$head$threshold)
    expect_identical(k202Kept(out$object), expected$head$kept)
  })
}

test_that("KIT-202 T1 (b) guard: at 0.5 both rules keep the three lipids", {
  out <- k202Method(k202Object("lipid_msdial"), 10, 0.5)
  expect_identical(k202Kept(out$object), 1:3)
})

# ---- T2: the declared step, both variants, each on its own owner -----------------------------------------------------

for (name in names(k202T1)) {
  for (variant in c("partitioned", "whole")) {
    test_that(sprintf("KIT-202 T2 (%s, %s): the declared step counts only the design's samples", name, variant), {
      skip_if_not(workflowOwnerCapable())
      expected <- k202T1[[name]]
      run <- k202Step(k202Object(name), 10, expected$proportion, variant)
      expect_identical(run$variant, variant)
      expect_identical(run$kept, expected$head$kept)
    })
  }
}

# ---- T3: the exported helpers ----------------------------------------------------------------------------------------

k202Helper <- function(object, threshold, proportion, ...) {
  table <- k202Assays(object)[[1L]]
  if (k202Lipid(object)) {
    lipidIntensityFilteringHelper(assay_table = table, min_lipid_intensity_threshold = threshold,
      lipids_proportion_of_samples_below_cutoff = proportion, lipid_id_column = k202IdColumn(object), ...)
  } else {
    metaboliteIntensityFilteringHelper(assay_table = table, min_metabolite_intensity_threshold = threshold,
      metabolites_proportion_of_samples_below_cutoff = proportion, metabolite_id_column = k202IdColumn(object), ...)
  }
}

k202HelperIds <- function(table, object) sort(as.integer(table[[k202IdColumn(object)]]))

for (name in names(k202T1)) {
  test_that(sprintf("KIT-202 T3 (a) (%s): the helper counts the sample columns it is given", name), {
    expected <- k202T1[[name]]
    object <- k202Object(name)
    result <- tryCatch(k202Helper(object, expected$head$threshold, expected$proportion,
      sample_columns = k202Samples(object)), error = function(e) NULL)
    expect_false(is.null(result))
    if (!is.null(result)) expect_identical(k202HelperIds(result, object), expected$head$kept)
  })
}

# T3 (b), D1 (d): without `sample_columns` every numeric column except the id is counted, with a classed warning.
# Metabolite at 4.9 and 0.2: no intensity is below 4.9, all 8 kept (base: ids 1 to 4 counted, share 1/5, removed).
# Lipid at 2.2 and 0.2: all 3 kept (base: ids 1 and 2 counted, share 1/5, removed).
k202T3b <- list(metab_msdial = list(threshold = 4.9, kept = 1:8), lipid_msdial = list(threshold = 2.2, kept = 1:3))

for (name in names(k202T3b)) {
  test_that(sprintf("KIT-202 T3 (b) (%s): without sample_columns the id column is not counted, with a warning", name), {
    expected <- k202T3b[[name]]
    object <- k202Object(name)
    expect_warning(result <- k202Helper(object, expected$threshold, 0.2),
      class = "omics_intensity_filter_samples_guessed")
    expect_identical(k202HelperIds(result, object), expected$kept)
  })
}

# ---- T4: unaffected tables (passes on base and head) -----------------------------------------------------------------

k202Unaffected <- c("lipid_canonical", "metab_lc", "metab_gc", "metab_combined", "module_ci_lipid", "module_ci_metab")

for (name in k202Unaffected) {
  test_that(sprintf("KIT-202 T4 (%s): an object whose numeric columns are its samples keeps the oracle's features", name), {
    object <- k202Object(name)
    id <- k202IdColumn(object)
    for (percentile in c(1, 10)) {
      for (proportion in c(0.2, 0.5)) {
        info <- sprintf("%s at (%s, %s)", name, percentile, proportion)
        out <- k202Method(object, percentile, proportion)
        for (assay in names(k202Assays(object))) {
          table <- k202Assays(object)[[assay]]
          oracle <- k202Oracle(table, k202Samples(object), percentile, proportion)
          # The method's message prints the threshold to 7 significant digits.
          expect_equal(unname(out$thresholds[[assay]]), oracle$threshold, tolerance = 1e-6, info = paste(info, assay))
          expect_identical(as.character(k202Assays(out$object)[[assay]][[id]]), as.character(table[[id]][oracle$keep]),
            info = paste(info, assay))
        }
      }
    }
  })
}

# ---- T5: a design sample absent or not numeric is refused (D3) -------------------------------------------------------

# A module-CI object with its first design sample removed from (absent) or turned to text in (text) its first assay;
# validity is bypassed with `methods::slot<-(check = FALSE)`.
k202Invalid <- function(name, kind = c("absent", "text")) {
  kind <- match.arg(kind)
  object <- k202Object(name)
  assays <- k202Assays(object)
  sample <- k202Samples(object)[[1L]]
  if (identical(kind, "absent")) assays[[1L]][[sample]] <- NULL else assays[[1L]][[sample]] <- as.character(assays[[1L]][[sample]])
  methods::slot(object, k202Binding(object)$data_slot, check = FALSE) <- assays
  object
}

for (name in c("module_ci_lipid", "module_ci_metab")) {
  for (kind in c("absent", "text")) {
    test_that(sprintf("KIT-202 T5 (a) (%s, %s): the S4 method refuses the design sample", name, kind), {
      expect_error(suppressMessages(k202Method(k202Invalid(name, kind), 10, 0.5)),
        class = "omics_intensity_filter_samples_invalid")
    })
  }
  test_that(sprintf("KIT-202 T5 (a) (%s, text): the partitioned step refuses the design sample", name), {
    skip_if_not(workflowOwnerCapable())
    object <- k202Invalid(name, "text")
    owner <- omicsTestOwner(object, "raw", config = list(qc = list(mode = "kit202")))$state_manager
    expect_error(kit172Quiet(runWorkflowOmicsReferenceIntensityFilter(k202Binding(object), owner, 10, 0.5,
      list(qc = list(mode = "kit202")))), class = "omics_intensity_filter_samples_invalid")
  })
}

test_that("KIT-202 T5 (b) guard: the application's objects are not refused", {
  for (name in c(names(k202T1), k202Unaffected)) {
    object <- k202Object(name)
    expect_no_error(k202Method(object, 10, 0.5))
  }
})

# ---- T6: the replay rule (D2 (a)) ------------------------------------------------------------------------------------

test_that("KIT-202 T6 (i): the replay rule of the intensity filter follows the step's marker", {
  markers <- get0("workflowReplayRuleMarkers", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  expect_true(!is.null(markers[["intensity_filter_samples"]]))
  if (!is.null(markers[["intensity_filter_samples"]])) {
    entry <- function(source) list(name = "intensity_filter_samples", source = source,
      value = encodeWorkflowRecordValue("design"), provenance = "fixed")
    for (operation in c("lipid_intensity_filter", "metab_intensity_filter")) {
      step <- function(parameters) list(state = "s", operation = operation, parameters = parameters)
      expect_identical(workflowReplayRule(step(list(entry("fixed"))), "intensity_filter_samples"), "design",
        info = operation)
      expect_identical(workflowReplayRule(step(list()), "intensity_filter_samples"), "numeric_columns", info = operation)
      expect_identical(workflowReplayRule(step(list(entry("registry"))), "intensity_filter_samples"), "numeric_columns",
        info = operation)
    }
  }
})

for (name in names(k202T1)) {
  test_that(sprintf("KIT-202 T6 (ii) (%s): under the earlier rule the head gives the base's values", name), {
    withr::local_options(k202Earlier)
    expected <- k202T1[[name]]
    out <- k202Method(k202Object(name), 10, expected$proportion)
    expect_equal(unname(out$thresholds[["LCMS_Pos"]]), expected$base$threshold)
    expect_identical(k202Kept(out$object), expected$base$kept)
    skip_if_not(workflowOwnerCapable())
    for (variant in c("partitioned", "whole")) {
      # New owners in new store roots (k202Step), never one holding a product of the other rule.
      run <- k202Step(k202Object(name), 10, expected$proportion, variant)
      expect_identical(run$variant, variant)
      expect_identical(run$kept, expected$base$kept, info = variant)
    }
  })

  test_that(sprintf("KIT-202 T6 (iii) (%s): a step made now records the marker", name), {
    skip_if_not(workflowOwnerCapable())
    run <- k202Step(k202Object(name), 10, k202T1[[name]]$proportion)
    expect_identical(k202Marker(run$record, run$state), list(source = "fixed", value = "design"))
  })
}

# T6 (iv): as KIT-190's T7 (iv) (branch kit-190-ruv-control-selection, tests/testthat/test-ruv-control-anova-grouping.R,
# k190Replay) and KIT-187's T8 (iv) (tests/testthat/test-prot-missing-value-counts.R): the owner's record written, the
# owner closed, a copy cut to the design state and the record replayed on it; the replay compared with the origin.
k202Paths <- function(root, lipid) {
  qc <- if (lipid) "lipid_qc_dir" else "metabolite_qc_dir"
  stats::setNames(list(file.path(root, "source"), file.path(root, "qc"), file.path(root, "da"), file.path(root, "pub")),
    c("source_dir", qc, "da_output_dir", "publication_graphs_dir"))
}

# One owner per rule in its own store root on an affected object; `earlier = TRUE` computes the step under the earlier
# rule, so its provenance has no marker, as a build before this change wrote it.
k202ReplayOwner <- function(root, object, proportion, earlier) {
  fixture <- if (k202Lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
  run <- function() kit172Quiet(runWorkflowOmicsReferenceIntensityFilter(k202Binding(object), fixture$owner, 10,
    proportion, fixture$workflow$config_list))
  if (earlier) withr::with_options(k202Earlier, run()) else run()
  fixture
}

k202Replay <- function(fixture, root, object) {
  prefix <- k202Binding(object)$prefix
  state <- paste0(prefix, "_intensity_filtered")
  kept <- k202Kept(fixture$owner$getState(state))
  json <- file.path(root, "analysis_record.json")
  writeWorkflowAnalysisRecord(fixture$owner, json)
  record <- jsonlite::read_json(json, simplifyVector = FALSE)
  closeDIAStateAdapter(fixture$owner)
  source_dir <- file.path(root, "source")
  store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE,
    full.names = TRUE)[[1L]]))
  copy <- tempfile("k202-copy-")
  dir.create(copy)
  stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
  harness <- new.env(parent = globalenv())
  sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
  prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(prefix, "_raw_data_s4")), api)
  inputs <- c(stats::setNames(file.path(source_dir, "import_sources", "LCMS_Pos_export.txt"), paste0(prefix, "_assay_1")),
    design_file = file.path(source_dir, "design_matrix.tab"))
  replay <- kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, k202Paths(root, k202Lipid(object)),
    tempfile("k202-scratch-"), api))
  identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
  replayed <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
  replayed_record <- buildWorkflowAnalysisRecord(replayed)
  closeDIAStateAdapter(replayed)
  list(record = record, replayed_record = replayed_record, mismatches = as.character(replay$mismatches), kept = kept,
    state = state)
}

for (name in names(k202T1)) {
  test_that(sprintf("KIT-202 T6 (iv) (%s): records made before and after the change replay, each with its own rule", name), {
    skip_if_not(workflowOwnerCapable())
    object <- k202Object(name)
    proportion <- k202T1[[name]]$proportion
    root_a <- withr::local_tempdir()
    root_b <- withr::local_tempdir()
    a <- k202Replay(k202ReplayOwner(root_a, object, proportion, earlier = TRUE), root_a, object)
    b <- k202Replay(k202ReplayOwner(root_b, object, proportion, earlier = FALSE), root_b, object)
    # Non-vacuity: the two rules give different filtered states on this table.
    expect_false(identical(a$kept, b$kept))
    expect_identical(a$mismatches, character())
    expect_identical(b$mismatches, character())
    expect_null(k202Marker(a$record, a$state))
    expect_null(k202Marker(a$replayed_record, a$state))
    expect_identical(k202Marker(b$record, b$state), list(source = "fixed", value = "design"))
    expect_identical(k202Marker(b$replayed_record, b$state), list(source = "fixed", value = "design"))
  })
}

# T6 (v): a managed step, as KIT-190's T7 (vi): the worker's entry point (omicsWorkerRun, as
# tests/testthat/test-kit129-run-path.R), completed with completeWorkflowAction.
test_that("KIT-202 T6 (v): a managed intensity filter step records the marker", {
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
    expect_identical(k202Marker(buildWorkflowAnalysisRecord(fixture$owner), step$state_name),
      list(source = "fixed", value = "design"), info = binding$prefix)
  }
})

# ---- T7: the metabolomics workbook's shape ---------------------------------------------------------------------------

test_that("KIT-202 T7: retention time, m/z and the id kept by the workbook are not counted", {
  # Head: threshold 16,410 over the 32 sample values; Metab_001 and Metab_002 are below in 2 of 4 and removed at 0.25.
  # Base (hand-computed): 56 counted values, threshold 2.725; features 1 and 2 have id and retention time below
  # (share 2/7), removed; features 3 and 4 one of 7.
  out <- k202Method(k202Object("workbook"), 10, 0.25)
  expect_identical(k202Kept(out$object), c(1L, 2L, 3L, 4L, 7L, 8L))
})
