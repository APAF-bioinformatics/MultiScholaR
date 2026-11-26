# KIT-260 (VAL-241: OD-1 (c), OD-2 (b), OD-3 (a), OD-4 (a)): the internal-standard pattern proposed for MS-DIAL files.
# Ticket section 7, rows E1 to E9 and G1 to G4 under the decided options. Objects are the ticket's section 4: two
# assays, 2 groups x 6 samples (KIT-172's design), 40 ordinary rows ("Metabolite 001" ...) and 20 named rows, id
# `Alignment ID` integers 1 to 60 (named rows 41 to 60), name column `Metabolite name` (metab) or `Name` (lipid), the
# pattern taken from the reader's default at the commit under test. Case k draws with
# set.seed(260000 + k, "Mersenne-Twister", "Inversion", "Rejection") before any package call. Expected texts are written
# out, so that a commit without them fails on the value.

KIT260_HEAD <- "^IS_|_d[0-9]+$|ISTD|^IS-|\\(d[0-9]+\\)"
KIT260_BASE <- "^IS_|_d[0-9]+$|ISTD|IS-|w/o MS2:"
KIT260_NAMED <- c("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine",
  "w/o MS2: Alanine", "w/o MS2: Serine", "w/o MS2: Valine", "w/o MS2: Citric acid", "w/o MS2: Lactic acid",
  "CIS-ACONITIC ACID", "HIS-LEU", "cis-Aconitic acid", "PC 15:0_18:1(d7)", "Unknown",
  "w/o MS2:Glucose", "no MS2: Glucose", "low score: Glucose", "w/o MS2:IS_Caffeine", "Caffeine-d9", NA)
# The 9 rows the base default used and the decided default does not: the `w/o MS2:` names and the capitalised names.
KIT260_RETURNED <- c("w/o MS2: Alanine", "w/o MS2: Serine", "w/o MS2: Valine", "w/o MS2: Citric acid", "w/o MS2: Lactic acid",
  "CIS-ACONITIC ACID", "HIS-LEU", "w/o MS2:Glucose", "w/o MS2:IS_Caffeine")
KIT260_STANDARDS_HEAD <- c("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine", "PC 15:0_18:1(d7)")
KIT260_NOTE <- paste("This project was imported with an internal-standard pattern that also matches features MS-DIAL",
  "annotated without MS/MS, or names containing 'IS-' in capitals; such features were used as internal standards",
  "and removed (assay A1, 9 rows; assay A2, 9 rows).")

kit260Seed <- function(k) set.seed(260000L + k, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
kit260Runs <- function() kit172Design()$Run

#' Two assays of log2 values: row levels N(20, 1.5^2), one loading per sample N(0, 0.5^2), noise N(0, 0.25^2).
kit260Logs <- function(n) {
  runs <- kit260Runs()
  loading <- stats::rnorm(length(runs), 0, 0.5)
  out <- list()
  for (a in c("A1", "A2")) {
    level <- stats::rnorm(n, 20, 1.5)
    noise <- matrix(stats::rnorm(n * length(runs), 0, 0.25), n)
    out[[a]] <- outer(level, rep(1, length(runs))) + outer(rep(1, n), loading) + noise
  }
  out
}

kit260Table <- function(y, names, name_col) {
  table <- data.frame(seq_len(nrow(y)), names, check.names = FALSE, stringsAsFactors = FALSE)
  names(table) <- c("Alignment ID", name_col)
  cbind(table, stats::setNames(as.data.frame(2^unname(y)), kit260Runs()))
}

#' The section-4 object of `class`, with `pattern` (the reader's default at this commit unless given) and `names`.
kit260Object <- function(class, k, pattern = NULL, names = c(sprintf("Metabolite %03d", 1:40), KIT260_NAMED)) {
  kit260Seed(k)
  y <- kit260Logs(length(names))
  args <- list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
  if (identical(class, "lipid")) {
    pattern <- pattern %||% getLipidomicsColumnDefaults("msdial")$is_pattern
    methods::new("LipidomicsAssayData", lipid_data = lapply(y, kit260Table, names = names, name_col = "Name"),
      design_matrix = kit172Design(), lipid_id_column = "Alignment ID", annotation_id_column = "Name", sample_id = "Run",
      group_id = "group", technical_replicate_id = "tech_rep_group", database_identifier_type = "Unknown",
      internal_standard_regex = pattern, args = args)
  } else {
    pattern <- pattern %||% getMetabolomicsColumnDefaults("msdial")$is_pattern
    methods::new("MetaboliteAssayData", metabolite_data = lapply(y, kit260Table, names = names, name_col = "Metabolite name"),
      design_matrix = kit172Design(), metabolite_id_column = "Alignment ID", annotation_id_column = "Metabolite name",
      sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group", database_identifier_type = "Unknown",
      internal_standard_regex = pattern, args = args)
  }
}

kit260NameColumn <- function(class) if (identical(class, "lipid")) "Name" else "Metabolite name"

#' An owner at the QC-complete state of `object`.
kit260Owner <- function(root, object, class) {
  fixture <- kit172Quiet(kit172Owner(root, object, class))
  prefix <- kit172Binding(class)$prefix
  kit172Quiet(commitWorkflowOmicsStoreStep(fixture$owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"),
    fixture$workflow$config_list, "QC complete"))
  fixture
}

#' An assay of a state as a data frame.
kit260Assay <- function(owner, binding, state, assay) as.data.frame(methods::slot(owner$getState(state), binding$data_slot)[[assay]],
  check.names = FALSE)

kit260Step <- function(owner, state) {
  record <- kit172Quiet(buildWorkflowAnalysisRecord(owner))
  for (step in record$steps) if (identical(step$state, state)) return(step)
  NULL
}

# --- E1, E2, G1, G2: the strings ---------------------------------------------------------------------------------------

test_that("KIT-260 E1: neither MS-DIAL default matches MS-DIAL's 'w/o MS2:' mark, with or without the space", {
  for (pattern in c(getMetabolomicsColumnDefaults("msdial")$is_pattern, getLipidomicsColumnDefaults("msdial")$is_pattern)) {
    expect_identical(stringr::str_detect(c("w/o MS2: Glucose", "w/o MS2:Glucose"), pattern), c(FALSE, FALSE))
  }
})

test_that("KIT-260 E2: ordinary names in capitals holding 'IS-' are not matched, on the default and on the custom preview", {
  for (pattern in c(getMetabolomicsColumnDefaults("msdial")$is_pattern, getLipidomicsColumnDefaults("msdial")$is_pattern)) {
    expect_identical(stringr::str_detect(c("CIS-ACONITIC ACID", "HIS-LEU"), pattern), c(FALSE, FALSE))
  }
  paths <- module_ci_metab_import_paths()
  file <- module_ci_metab_write_table(module_ci_metab_custom_assay(), file.path(paths$raw_dir, "custom.tsv"))
  preview <- suppressWarnings(previewOmicsImport(file, omicGuiConfig("metabolomics")$import))
  expect_identical(stringr::str_detect(c("CIS-ACONITIC ACID", "HIS-LEU"), preview$updates$is_pattern), c(FALSE, FALSE))
})

test_that("KIT-260 G1: the intended standards are matched by the default; the parenthesised label too (OD-1 (c))", {
  for (pattern in c(getMetabolomicsColumnDefaults("msdial")$is_pattern, getLipidomicsColumnDefaults("msdial")$is_pattern)) {
    # Guard: matched before and after the change.
    expect_true(all(stringr::str_detect(c("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine"), pattern)))
    # OD-1 (c): MS-DIAL's own form for a labelled lipid standard (not matched before).
    expect_true(stringr::str_detect("PC 15:0_18:1(d7)", pattern))
  }
})

test_that("KIT-260 G2: every other format's default is unchanged; the import proposes the MS-DIAL default as before", {
  expected <- list(progenesis = "^IS_|ISTD", xcms = "^IS_|ISTD", compound_discoverer = "^IS_|ISTD|Internal", custom = NA_character_)
  for (format in names(expected)) {
    expect_identical(getMetabolomicsColumnDefaults(format)$is_pattern, expected[[format]], info = format)
    expect_identical(getLipidomicsColumnDefaults(format)$is_pattern, expected[[format]], info = format)
  }
  expect_identical(getLipidomicsColumnDefaults("lipidsearch")$is_pattern, "^IS_|ISTD|Internal|_IS")
  lipid_paths <- module_ci_lipid_import_paths()
  lipid_files <- c(
    lipidsearch = module_ci_lipid_write_table(module_ci_lipid_lipidsearch_assay(), file.path(lipid_paths$raw_dir, "ls.txt")),
    msdial = module_ci_lipid_write_table(module_ci_lipid_msdial_assay(), file.path(lipid_paths$raw_dir, "msdial_positive.tsv")),
    custom = module_ci_lipid_write_table(module_ci_lipid_custom_assay(), file.path(lipid_paths$raw_dir, "custom.tsv")),
    malformed = module_ci_lipid_write_malformed(file.path(lipid_paths$raw_dir, "malformed.tsv")))
  for (name in names(lipid_files)) {
    preview <- suppressWarnings(previewOmicsImport(lipid_files[[name]], omicGuiConfig("lipidomics")$import))
    want <- if (identical(name, "lipidsearch")) "^IS_|ISTD|Internal|_IS" else getLipidomicsColumnDefaults("msdial")$is_pattern
    expect_identical(preview$updates$is_pattern, want, info = name)
  }
  metab_paths <- module_ci_metab_import_paths()
  metab_files <- c(
    msdial = module_ci_metab_write_table(module_ci_metab_msdial_assay(), file.path(metab_paths$raw_dir, "msdial_positive.tsv")),
    custom = module_ci_metab_write_table(module_ci_metab_custom_assay(), file.path(metab_paths$raw_dir, "custom.tsv")))
  for (name in names(metab_files)) {
    preview <- suppressWarnings(previewOmicsImport(metab_files[[name]], omicGuiConfig("metabolomics")$import))
    expect_identical(preview$updates$is_pattern, getMetabolomicsColumnDefaults("msdial")$is_pattern, info = name)
  }
})

test_that("KIT-260 E9: the custom-layout preview proposes the decided pattern (metabolomics and lipidomics)", {
  metab_paths <- module_ci_metab_import_paths()
  metab <- module_ci_metab_write_table(module_ci_metab_custom_assay(), file.path(metab_paths$raw_dir, "custom.tsv"))
  expect_identical(suppressWarnings(previewOmicsImport(metab, omicGuiConfig("metabolomics")$import))$updates$is_pattern, KIT260_HEAD)
  lipid_paths <- module_ci_lipid_import_paths()
  lipid <- module_ci_lipid_write_table(module_ci_lipid_custom_assay(), file.path(lipid_paths$raw_dir, "custom.tsv"))
  expect_identical(suppressWarnings(previewOmicsImport(lipid, omicGuiConfig("lipidomics")$import))$updates$is_pattern, KIT260_HEAD)
})

# --- E3, E4: pattern-mode ITSD on the section-4 objects ------------------------------------------------------------------

for (class in c("metab", "lipid")) {
  test_that(sprintf("KIT-260 %s: pattern-mode ITSD uses the 5 labelled standards; the 9 other rows stay and are normalised (%s)",
    if (identical(class, "metab")) "E3" else "E4", class), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    fixture <- kit260Owner(root, kit260Object(class, 1L), class)
    withr::defer(closeDIAStateAdapter(fixture$owner))
    kit172Quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, fixture$owner, "median", NULL, fixture$workflow$config_list))
    after <- fixture$owner$getState(paste0(prefix, "_itsd_norm"))
    expect_identical(after@args$ITSDNormalization$itsd_counts_per_assay, list(A1 = 5L, A2 = 5L))
    name_col <- kit260NameColumn(class)
    for (assay in c("A1", "A2")) {
      before_table <- kit260Assay(fixture$owner, binding, paste0(prefix, "_qc_complete"), assay)
      after_table <- kit260Assay(fixture$owner, binding, paste0(prefix, "_itsd_norm"), assay)
      kept <- after_table[[name_col]]
      expect_true(all(KIT260_RETURNED %in% kept), info = assay)
      expect_false(any(KIT260_STANDARDS_HEAD %in% kept), info = assay)
      # The missing-annotation row is dropped in pattern mode on both sides (P12 (a); not changed here).
      expect_identical(nrow(after_table), 54L, info = assay)
      # The returned rows are scaled by the same per-sample factor as an ordinary row.
      ratio <- function(id) {
        b <- unlist(before_table[before_table[["Alignment ID"]] == id, kit260Runs()])
        a <- unlist(after_table[after_table[["Alignment ID"]] == id, kit260Runs()])
        a / b
      }
      ordinary <- ratio(1L)
      for (id in before_table[["Alignment ID"]][before_table[[name_col]] %in% KIT260_RETURNED]) {
        expect_equal(ratio(id), ordinary, tolerance = 1e-12, info = paste(assay, id))
      }
    }
  })
}

# --- E5, G-type integration naming (OD-4 (a)) ----------------------------------------------------------------------------

test_that("KIT-260 E5: the integration input, ITSD unticked, marks exactly the decided default's rows (OD-4 (a))", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  binding <- kit172Binding("metab")
  fixture <- kit260Owner(root, kit260Object("metab", 1L), "metab")
  withr::defer(closeDIAStateAdapter(fixture$owner))
  config <- fixture$workflow$config_list
  kit172Quiet({
    runWorkflowOmicsReferenceLog2(binding, fixture$owner, 0, config)
    runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, fixture$owner, "cyclicloess", config)
  })
  source <- kit172Quiet(withWorkflowStateScope(fixture$owner, function() omicsIntegrationInputSource(binding, fixture$owner)))
  for (view in names(source$features)) {
    features <- source$features[[view]]
    expect_identical(sort(as.character(features[["Metabolite name"]][features$is_internal_standard %in% TRUE])),
      sort(KIT260_STANDARDS_HEAD), info = view)
  }
  # OD-4 (a): the marked rows named per data set, "display name (identifier)".
  views <- names(source$features)
  object <- methods::new("IntegrationViewsData", views = stats::setNames(lapply(views, function(v) data.frame(feature_id = character())), views),
    features = source$features, inputs = list(metabolomics = list(omic = "metabolomics", views = views, manifest = source$manifest)))
  named <- integrationInternalStandardNames(object, views)
  expect_identical(named[[views[[1L]]]], c("IS_Caffeine (41)", "Caffeine_d9 (42)", "Caffeine-ISTD (43)", "IS-Leucine (44)",
    "PC 15:0_18:1(d7) (53)"))
  expect_length(integrationInternalStandardLegacy(object, views), 0L)
  preview <- list(parameters = list(missing_samples = "drop_incomplete", technical_replicates = "average", internal_standards = "drop"),
    aligned = list(view_summary = data.frame(view = views, sparse_dropped = 0L), shared_samples = 12L), all_samples = 12L,
    replicates = data.frame(view = views, runs = 12L, biological_samples = 12L, repeated = 0L, standards = 5L),
    standard_names = named, standard_legacy = list())
  facts <- integrationAlignmentFacts(preview)
  expect_true(sprintf("%s: internal standards marked: IS_Caffeine (41); Caffeine_d9 (42); Caffeine-ISTD (43); IS-Leucine (44); PC 15:0_18:1(d7) (53).",
    views[[1L]]) %in% facts)
})

# --- E6, E7: what the user is shown (OD-2 (b)) ---------------------------------------------------------------------------

#' A runner that runs the managed actions in this process (test-kit172-downstream.R's, copied).
kit260InProcessRunner <- function(manager) {
  list(progress = function(...) list(inc = function(...) invisible(NULL), close = function() invisible(NULL)),
    run = function(actions, done = function() NULL, failed = function(error) stop(error)) {
      outcome <- tryCatch({
        for (action in actions) {
          if (is.function(action$before)) action$before()
          if (is.null(action$context)) {
            action$deliver(list(display = list()))
          } else {
            step <- list(operation = action$operation, parameters = action$parameters, state_name = action$context$state_name,
              description = action$context$description)
            result <- runWorkflowOmicsReferenceStep(manager, step, action$context$config)
            action$deliver(list(stateReference = result$reference, value = result$value))
          }
        }
        NULL
      }, error = function(e) e)
      if (inherits(outcome, "error")) failed(outcome) else done()
    })
}

#' One normalisation run of the tab (route `in_process` or `managed`) with the log and notifications recorded.
kit260TabRun <- function(fixture, class, route) {
  binding <- kit172Binding(class)
  binding$qc_plots <- function(...) invisible(NULL)
  normData <- new.env(parent = emptyenv())
  log <- character()
  notes <- list()
  addLog <- function(entry) log <<- c(log, entry)
  testthat::local_mocked_bindings(showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
    type = "default", session = NULL) {
    notes[[length(notes) + 1L]] <<- list(type = type, message = paste(as.character(ui), collapse = " "))
    invisible(NULL)
  }, .package = "shiny")
  settings <- list(apply_itsd = TRUE, itsd_aggregation = "median", itsd_feature_ids = NULL, log_offset = 0,
    norm_method = "cyclicloess", ruv_mode = "skip", ruv_parameters = list(), grouping_variable = "group", shape_variable = "group",
    assay_names = NULL, qc_dir = NULL, experiment_paths = list())
  kit172Quiet(if (identical(route, "in_process")) {
    runWorkflowOmicsInProcessNormalisation(binding, fixture$owner, settings, normData, fixture$workflow$config_list, addLog,
      advance = function(detail) NULL, notify = function(...) shiny::showNotification(...), finish = function() NULL,
      itsd_table = NULL)
  } else {
    runWorkflowOmicsManagedNormalisation(binding, kit260InProcessRunner(fixture$owner), fixture$owner, settings, normData,
      fixture$workflow$config_list, addLog, finish = function() NULL, failed = function(error) addLog(paste("ERROR:", conditionMessage(error))))
  })
  list(log = log, notes = notes)
}

test_that("KIT-260 E6: the tab names the rows used per assay (log and notification, both routes) and the record lists them", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  owners <- list()
  withr::defer(for (owner in owners) closeDIAStateAdapter(owner))
  rows <- paste("IS_Caffeine (Alignment ID 41); Caffeine_d9 (Alignment ID 42); Caffeine-ISTD (Alignment ID 43);",
    "IS-Leucine (Alignment ID 44); PC 15:0_18:1(d7) (Alignment ID 53)")
  for (route in c("in_process", "managed")) {
    fixture <- kit260Owner(file.path(root, route), kit260Object("metab", 1L), "metab")
    owners[[length(owners) + 1L]] <- fixture$owner
    run <- kit260TabRun(fixture, "metab", route)
    expect_false(any(startsWith(run$log, "ERROR")), info = paste(route, paste(run$log, collapse = " | ")))
    for (assay in c("A1", "A2")) {
      line <- sprintf("ITSD, assay %s: 5 rows used as internal standards, then removed: %s", assay, rows)
      expect_identical(sum(run$log == line), 1L, info = paste(route, assay))
    }
    text <- paste("Internal-standard normalisation used 5 (A1) and 5 (A2) rows as internal standards and removed them;",
      "the log names them.")
    expect_identical(sum(vapply(run$notes, function(note) identical(note$message, text), logical(1))), 1L, info = route)
    step <- kit260Step(fixture$owner, "metab_itsd_norm")
    expect_false(is.null(step), info = route)
    expect_identical(unlist(step$outcomes$itsd_features_per_assay$A1), strsplit(rows, "; ", fixed = TRUE)[[1L]], info = route)
    expect_null(step$outcomes$itsd_note)
  }
})

test_that("KIT-260 E7: the pattern box's help text says where it is matched, that it is case-sensitive, and what a match does", {
  for (omic in c("metabolomics", "lipidomics")) {
    html <- paste(as.character(shiny::tagList(omicImportAuxInputs(shiny::NS("x"), omicGuiConfig(omic)$import))), collapse = " ")
    expect_true(grepl("case-sensitive", html, fixed = TRUE), info = omic)
    expect_true(grepl("annotation column", html, fixed = TRUE), info = omic)
    expect_true(grepl("is then removed from the data", html, fixed = TRUE), info = omic)
  }
})

# --- E8: a project whose stored pattern is the 0.9.0 default (OD-3 (a)) -------------------------------------------------

test_that("KIT-260 E8: the stored 0.9.0 pattern governs; the record and the report carry the note only when such rows matched", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  binding <- kit172Binding("metab")
  plain <- c(sprintf("Metabolite %03d", 1:40), ifelse(KIT260_NAMED %in% KIT260_RETURNED, sprintf("Plain %02d", seq_along(KIT260_NAMED)),
    KIT260_NAMED))
  cases <- list(legacy = kit260Object("metab", 8L, pattern = KIT260_BASE),
    twin = kit260Object("metab", 8L, pattern = KIT260_BASE, names = plain))
  for (case in names(cases)) {
    fixture <- kit260Owner(file.path(root, case), cases[[case]], "metab")
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit172Quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config))
    itsd <- owner$getState("metab_itsd_norm")
    # (i) The stored pattern governs: 13 rows used on the legacy object (as the base commit), 4 on the twin.
    expect_identical(itsd@args$ITSDNormalization$itsd_counts_per_assay,
      if (identical(case, "legacy")) list(A1 = 13L, A2 = 13L) else list(A1 = 4L, A2 = 4L), info = case)
    step <- kit260Step(owner, "metab_itsd_norm")
    # (ii) The note: in the record (the record viewer prints every outcome), and in the study parameters and the report.
    if (identical(case, "legacy")) expect_identical(step$outcomes$itsd_note, KIT260_NOTE) else expect_null(step$outcomes$itsd_note)
    dir <- file.path(root, paste0("sp-", case))
    dir.create(dir)
    kit172Quiet(utils::capture.output(createWorkflowArgsFromConfig("kit260", source_dir_path = dir, final_s4_object = itsd,
      workflow_data = list(state_manager = owner))))
    lines <- readLines(file.path(dir, "study_parameters.txt"), warn = FALSE)
    note_line <- paste("* Note:", KIT260_NOTE)
    expect_identical(sum(trimws(lines) == note_line), if (identical(case, "legacy")) 1L else 0L, info = case)
    if (identical(case, "legacy")) {
      template <- system.file("reports", "metabolomics", "metabolomics_report.rmd", package = "MultiScholaR")
      env <- new.env(parent = asNamespace("MultiScholaR"))
      env$`%||%` <- function(a, b) if (!is.null(a)) a else b
      env$d66_log <- function(...) invisible(NULL)
      report <- readLines(template, warn = FALSE, encoding = "UTF-8")
      for (label in c("Hierarchical Parser", "Sample Analysis Functions", "Method Section")) {
        start <- grep(paste0("^```\\{r ", label, "[,}]"), report)
        end <- start + which(report[-seq_len(start)] == "```")[[1L]]
        chunk <- tempfile(fileext = ".Rmd")
        writeLines(report[start:end], chunk, useBytes = TRUE)
        code <- tempfile(fileext = ".R")
        suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
        for (e in parse(code, keep.source = FALSE, encoding = "UTF-8")) suppressWarnings(eval(e, env))
      }
      numbers <- env$ExtractStudyParametersNew(env$ReadStudyParamsNew(file.path(dir, "study_parameters.txt")))
      expect_identical(numbers$itsd_note, KIT260_NOTE)
      method_text <- paste(suppressWarnings(env$FormatMethodDetails(numbers)), collapse = " ")
      expect_true(grepl(KIT260_NOTE, method_text, fixed = TRUE))
      # (i) Replayed with today's default: the stored pattern still governs, every state as recorded.
      json <- file.path(root, "legacy-record.json")
      kit172Quiet(writeWorkflowAnalysisRecord(owner, json))
      history <- unlist(owner$getHistory())
      closeDIAStateAdapter(owner)
      source_dir <- fixture$source_dir
      copy <- withr::local_tempdir()
      stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
      store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
      copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
      record <- jsonlite::read_json(json, simplifyVector = FALSE)
      harness <- new.env(parent = globalenv())
      sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
      api <- function(name) get(name, asNamespace("MultiScholaR"))
      prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", "metab_raw_data_s4"), api)
      inputs <- c(stats::setNames(file.path(source_dir, "import_sources", c("A1_export.txt", "A2_export.txt")), c("metab_assay_1", "metab_assay_2")),
        design_file = file.path(source_dir, "design_matrix.tab"))
      replay <- kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, fixture$paths, withr::local_tempdir(), api))
      expect_identical(replay$states, history)
      expect_identical(as.character(replay$mismatches), character())
    } else {
      closeDIAStateAdapter(owner)
    }
  }
})

# --- G3, G4 ----------------------------------------------------------------------------------------------------------------

test_that("KIT-260 G3: ITSD with hand-picked standards is the arithmetic it was (median, mean-centred, standards removed)", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  binding <- kit172Binding("metab")
  fixture <- kit260Owner(root, kit260Object("metab", 1L), "metab")
  withr::defer(closeDIAStateAdapter(fixture$owner))
  ids <- list(A1 = as.character(41:44), A2 = as.character(41:44))
  kit172Quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, fixture$owner, "median", ids, fixture$workflow$config_list))
  for (assay in c("A1", "A2")) {
    before <- kit260Assay(fixture$owner, binding, "metab_qc_complete", assay)
    after <- kit260Assay(fixture$owner, binding, "metab_itsd_norm", assay)
    runs <- kit260Runs()
    standard <- apply(as.matrix(before[before[["Alignment ID"]] %in% 41:44, runs]), 2L, stats::median, na.rm = TRUE)
    factor <- mean(standard) / standard
    kept <- before[!before[["Alignment ID"]] %in% 41:44, , drop = FALSE]
    expected <- sweep(as.matrix(kept[, runs]), 2L, factor, `*`)
    at <- match(kept[["Alignment ID"]], after[["Alignment ID"]])
    expect_false(anyNA(at), info = assay)
    expect_equal(unname(as.matrix(after[at, runs])), unname(expected), tolerance = 1e-12, info = assay)
  }
})

test_that("KIT-260 G4: KIT-242's C14 object, 'w/o MS2:' rows kept as features by the default", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  binding <- kit172Binding("metab")
  names <- sprintf("Metabolite %03d", 1:40)
  names[1:5] <- sprintf("w/o MS2: Feature %d", 1:5)
  names[6:7] <- c("PC 15:0_18:1(d7)", "LPC 18:1(d7)")
  object <- kit260Object("metab", 14L, names = names)
  fixture <- kit260Owner(root, object, "metab")
  withr::defer(closeDIAStateAdapter(fixture$owner))
  kit172Quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, fixture$owner, "median", NULL, fixture$workflow$config_list))
  after <- kit260Assay(fixture$owner, binding, "metab_itsd_norm", "A1")
  expect_true(all(sprintf("w/o MS2: Feature %d", 1:5) %in% after[["Metabolite name"]]))
  expect_false(any(c("PC 15:0_18:1(d7)", "LPC 18:1(d7)") %in% after[["Metabolite name"]]))
})

# --- Correction pass (AUDIT-KIT-260-impl-20261126): I02 (a), I05, I10 ---------------------------------------------------

#' Run pattern-mode ITSD with the product as a build before KIT-260 wrote it (the kernel's rows-used fields absent),
#' then the record, the study-parameters file and the report as an older project reopened at this version builds them.
kit260OlderItsd <- function(root, class, object) {
  binding <- kit172Binding(class)
  fixture <- kit260Owner(root, object, class)
  local({
    testthat::local_mocked_bindings(omicsItsdRowsUsed = function(...) list(), .package = "MultiScholaR")
    kit172Quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, fixture$owner, "median", NULL, fixture$workflow$config_list))
  })
  fixture
}

kit260StudyLines <- function(root, owner, state) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  kit172Quiet(utils::capture.output(createWorkflowArgsFromConfig("kit260", source_dir_path = root, final_s4_object = owner$getState(state),
    workflow_data = list(state_manager = owner))))
  readLines(file.path(root, "study_parameters.txt"), warn = FALSE)
}

#' The report's parser and Methods paragraph on a study-parameters file (the template's own chunks).
kit260ReportNote <- function(class, file) {
  spec <- if (identical(class, "lipid")) list(dir = "lipidomics", file = "lipidomics_report.rmd", extract = "ExtractLipidStudyParameters") else
    list(dir = "metabolomics", file = "metabolomics_report.rmd", extract = "ExtractStudyParametersNew")
  env <- new.env(parent = asNamespace("MultiScholaR"))
  env$`%||%` <- function(a, b) if (!is.null(a)) a else b
  env$d66_log <- function(...) invisible(NULL)
  report <- readLines(system.file("reports", spec$dir, spec$file, package = "MultiScholaR"), warn = FALSE, encoding = "UTF-8")
  for (label in c("Hierarchical Parser", "Sample Analysis Functions", "Method Section")) {
    start <- grep(paste0("^```\\{r ", label, "[,}]"), report)
    end <- start + which(report[-seq_len(start)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(report[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    for (e in parse(code, keep.source = FALSE, encoding = "UTF-8")) suppressWarnings(eval(e, env))
  }
  numbers <- env[[spec$extract]](env$ReadStudyParamsNew(file))
  list(note = numbers$itsd_note, method = paste(suppressWarnings(env$FormatMethodDetails(numbers)), collapse = " "))
}

kit260Plain <- function() c(sprintf("Metabolite %03d", 1:40), ifelse(KIT260_NAMED %in% KIT260_RETURNED,
  sprintf("Plain %02d", seq_along(KIT260_NAMED)), KIT260_NAMED))

for (class in c("metab", "lipid")) {
  test_that(sprintf("KIT-260 I02 (record): an older project's record carries the note, read from its stored states (%s)", class), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    prefix <- kit172Binding(class)$prefix
    for (case in c("legacy", "twin")) {
      object <- if (identical(case, "legacy")) kit260Object(class, 8L, pattern = KIT260_BASE) else
        kit260Object(class, 8L, pattern = KIT260_BASE, names = kit260Plain())
      # The ITSD product is as an older build wrote it: no rows-used fields.
      fixture <- kit260OlderItsd(file.path(root, case), class, object)
      step <- kit260Step(fixture$owner, paste0(prefix, "_itsd_norm"))
      if (identical(case, "legacy")) {
        expect_identical(step$outcomes$itsd_note, KIT260_NOTE, info = class)
        expect_length(unlist(step$outcomes$itsd_features_per_assay$A1), 13L)
      } else {
        expect_null(step$outcomes$itsd_note)
      }
      closeDIAStateAdapter(fixture$owner)
    }
  })

  test_that(sprintf("KIT-260 I02 (study parameters): an older project's study-parameters file carries the note (%s)", class), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    prefix <- kit172Binding(class)$prefix
    for (case in c("legacy", "twin")) {
      object <- if (identical(case, "legacy")) kit260Object(class, 8L, pattern = KIT260_BASE) else
        kit260Object(class, 8L, pattern = KIT260_BASE, names = kit260Plain())
      fixture <- kit260OlderItsd(file.path(root, case), class, object)
      lines <- kit260StudyLines(file.path(root, paste0("sp-", case)), fixture$owner, paste0(prefix, "_itsd_norm"))
      expect_identical(sum(trimws(lines) == paste("* Note:", KIT260_NOTE)), if (identical(case, "legacy")) 1L else 0L,
        info = paste(class, case))
      closeDIAStateAdapter(fixture$owner)
    }
  })

  test_that(sprintf("KIT-260 I02 (report): an older project's report carries the note (%s)", class), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    prefix <- kit172Binding(class)$prefix
    fixture <- kit260OlderItsd(file.path(root, "legacy"), class, kit260Object(class, 8L, pattern = KIT260_BASE))
    withr::defer(closeDIAStateAdapter(fixture$owner))
    dir <- file.path(root, "sp")
    kit260StudyLines(dir, fixture$owner, paste0(prefix, "_itsd_norm"))
    report <- kit260ReportNote(class, file.path(dir, "study_parameters.txt"))
    expect_identical(report$note, KIT260_NOTE, info = class)
    expect_true(grepl(KIT260_NOTE, report$method, fixed = TRUE), info = class)
  })
}

test_that("KIT-260 I10: the lipid study parameters and lipid report carry the note on a project run at this version", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  binding <- kit172Binding("lipid")
  fixture <- kit260Owner(root, kit260Object("lipid", 8L, pattern = KIT260_BASE), "lipid")
  withr::defer(closeDIAStateAdapter(fixture$owner))
  kit172Quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, fixture$owner, "median", NULL, fixture$workflow$config_list))
  dir <- file.path(root, "sp")
  lines <- kit260StudyLines(dir, fixture$owner, "lipid_itsd_norm")
  expect_identical(sum(trimws(lines) == paste("* Note:", KIT260_NOTE)), 1L)
  report <- kit260ReportNote("lipid", file.path(dir, "study_parameters.txt"))
  expect_identical(report$note, KIT260_NOTE)
  expect_true(grepl(KIT260_NOTE, report$method, fixed = TRUE))
})

test_that("KIT-260 I10: the integration tab states the note when marked rows were matched through the removed alternatives", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  binding <- kit172Binding("metab")
  fixture <- kit260Owner(root, kit260Object("metab", 8L, pattern = KIT260_BASE), "metab")
  withr::defer(closeDIAStateAdapter(fixture$owner))
  config <- fixture$workflow$config_list
  kit172Quiet({
    runWorkflowOmicsReferenceLog2(binding, fixture$owner, 0, config)
    runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, fixture$owner, "cyclicloess", config)
  })
  source <- kit172Quiet(withWorkflowStateScope(fixture$owner, function() omicsIntegrationInputSource(binding, fixture$owner)))
  views <- names(source$features)
  object <- methods::new("IntegrationViewsData", views = stats::setNames(lapply(views, function(v) data.frame(feature_id = character())), views),
    features = source$features, inputs = list(metabolomics = list(omic = "metabolomics", views = views, manifest = source$manifest)))
  legacy <- integrationInternalStandardLegacy(object, views)
  expect_identical(lengths(legacy), stats::setNames(c(9L, 9L), views))
  preview <- list(parameters = list(missing_samples = "drop_incomplete", technical_replicates = "average", internal_standards = "drop"),
    aligned = list(view_summary = data.frame(view = views, sparse_dropped = 0L), shared_samples = 12L), all_samples = 12L,
    replicates = data.frame(view = views, runs = 12L, biological_samples = 12L, repeated = 0L, standards = 13L),
    standard_names = integrationInternalStandardNames(object, views), standard_legacy = legacy)
  note <- paste("This project was imported with an internal-standard pattern that also matches features MS-DIAL annotated",
    "without MS/MS, or names containing 'IS-' in capitals; such features were used as internal standards and removed",
    sprintf("(data set %s, 9 rows; data set %s, 9 rows).", views[[1L]], views[[2L]]))
  expect_true(note %in% integrationAlignmentFacts(preview))
})

test_that("KIT-260 I10: the integration model and its record name the marked rows per data set (stand-in trainer)", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  integration_dir <- file.path(root, "integration")
  input <- function(omic, views, samples, standards = list()) {
    matrices <- lapply(views, function(m) data.frame(feature_id = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE))
    features <- lapply(stats::setNames(names(views), names(views)), function(view) {
      table <- data.frame(feature_id = rownames(views[[view]]), display_name = tolower(rownames(views[[view]])), stringsAsFactors = FALSE)
      if (!is.null(standards[[view]])) table$is_internal_standard <- standards[[view]]
      table
    })
    manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2", correction = NULL,
      source_state = list(name = "kit260-test"), sample_id_column = "Run", group_column = "group",
      technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
    writeIntegrationInput(file.path(integration_dir, omic), manifest, samples, matrices, features)
  }
  samples <- function(runs, biological) data.frame(Run = runs, group = sub("_.*$", "", biological),
    replicates = as.integer(sub("^.*_", "", biological)), tech_rep_group = biological, stringsAsFactors = FALSE)
  biological <- c("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
  p <- rbind(P1 = c(10, 11, 12, 13, 14, 15), P2 = c(5, 7, 6, 9, 8, 10), P3 = c(1, 1.5, 2, 2.5, 3.5, 3))
  colnames(p) <- paste0("p", seq_along(biological))
  input("proteomics", list(proteomics = p), samples(colnames(p), biological))
  l <- rbind(L1 = c(20, 22, 21, 23, 25, 27), L2 = c(3, 4, 4.5, 5, 6, 7), ITSD1 = rep(15, 6))
  colnames(l) <- paste0("l", seq_along(biological))
  input("lipidomics", list(lipidomics_LCMS_Pos = l), samples(colnames(l), biological),
    standards = list(lipidomics_LCMS_Pos = c(FALSE, FALSE, TRUE)))
  manager <- newIntegrationStateOwner(file.path(root, "store"))
  withr::defer(closeWorkflowStateOwner(manager))
  scan <- integrationScanProject(integration_dir)
  integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
  trainer <- function(views, options) {
    samples <- colnames(views[[1L]])
    k <- options$num_factors
    labels <- paste0("Factor", seq_len(k))
    list(factors = matrix(seq_len(length(samples) * k) / 10, length(samples), k, dimnames = list(samples, labels)),
      weights = lapply(views, function(v) matrix(0.5, nrow(v), k, dimnames = list(rownames(v), labels))),
      r2_per_factor = matrix(1, k, length(views), dimnames = list(labels, names(views))),
      r2_total = stats::setNames(seq_along(views) + 0.25, names(views)), elbo = c(-20, -10, NaN), iterations = 2L,
      options = options, versions = list(MOFA2 = "stand-in"), messages = character(),
      model_file = list(bytes = as.raw(1:3), sha256 = strrep("a", 64L), size = 3))
  }
  run <- suppressWarnings(runIntegrationModel(manager, list(views = sort(scan$views$view)), probe = function() list(available = TRUE,
    reason = NA_character_), trainer = trainer))
  expect_identical(run$value$internal_standards_marked, list(lipidomics_LCMS_Pos = "itsd1 (ITSD1)"))
  outcomes <- integrationModelRecordOutcomes(manager$.context, list(reference = run$reference))
  expect_identical(outcomes$internal_standards_marked, list(lipidomics_LCMS_Pos = "itsd1 (ITSD1)"))
})

test_that("KIT-260 I05: a record from a build before KIT-260 whose ITSD used standards replays with identical values", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  binding <- kit172Binding("metab")
  # The ITSD step as a build before KIT-260 made it (VAL-276): under the earlier token of the replay marker KIT-242 and
  # KIT-260 share (`MultiScholaR.itsd_outcome` = "applied_flag"), so the step records no marker and its product carries
  # none of the rows-used fields.
  fixture <- kit260Owner(root, kit260Object("metab", 8L, pattern = KIT260_BASE), "metab")
  owner <- fixture$owner
  config <- fixture$workflow$config_list
  withr::with_options(list(MultiScholaR.itsd_outcome = "applied_flag"),
    kit172Quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config)))
  # The gate itself (AUDIT-KIT-260-impl CC03): under the earlier token the origin's ITSD product holds only `final`.
  origin_spec <- readWorkflowStateMetadata(owner, owner$.context$snapshot$states$metab_itsd_norm)$data
  origin_product <- hydrateSavedWorkflowEntity(owner$.context$store, origin_spec$products$metab_itsd_normalisation_science)$value
  expect_identical(names(origin_product), "final")
  kit172Quiet({
    runWorkflowOmicsReferenceLog2(binding, owner, 0, config)
    runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config)
  })
  json <- file.path(root, "older-record.json")
  local({
    # The record as the older build wrote it: the ITSD step's outcome without KIT-260's fields.
    testthat::local_mocked_bindings(omicsItsdRowsFromStates = function(...) list(), .package = "MultiScholaR")
    kit172Quiet(writeWorkflowAnalysisRecord(owner, json))
  })
  history <- unlist(owner$getHistory())
  closeDIAStateAdapter(owner)
  record <- jsonlite::read_json(json, simplifyVector = FALSE)
  itsd_step <- Filter(function(step) identical(step$state, "metab_itsd_norm"), record$steps)[[1L]]
  expect_null(itsd_step$outcomes$itsd_features_per_assay)
  source_dir <- fixture$source_dir
  copy <- withr::local_tempdir()
  stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
  store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
  copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
  harness <- new.env(parent = globalenv())
  sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", "metab_raw_data_s4"), api)
  inputs <- c(stats::setNames(file.path(source_dir, "import_sources", c("A1_export.txt", "A2_export.txt")), c("metab_assay_1", "metab_assay_2")),
    design_file = file.path(source_dir, "design_matrix.tab"))
  # The replay runs this version's kernel (the rows-used fields are added to the ITSD product).
  replay <- kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, fixture$paths, withr::local_tempdir(), api))
  expect_identical(replay$states, history)
  # What the comparison does with the additive fields: every state and every product judged by its sha256.
  expect_identical(as.character(replay$mismatches), character())
  # The values: each replayed state's assays identical to the origin's.
  origin <- newDIAStateAdapter(newWorkflowEntityStore(store_root, readSavedStoreIdentity(file.path(store_root, ".workflow-entities",
    "store.json"))$project_id, readSavedStoreIdentity(file.path(store_root, ".workflow-entities", "store.json"))$workflow_id), reopen = TRUE)
  withr::defer(closeDIAStateAdapter(origin))
  identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
  again <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
  withr::defer(closeDIAStateAdapter(again))
  for (state in c("metab_itsd_norm", "metab_log2", "metab_normalized")) {
    a <- methods::slot(origin$getState(state), binding$data_slot)
    b <- methods::slot(again$getState(state), binding$data_slot)
    for (assay in names(a)) expect_identical(as.data.frame(b[[assay]]), as.data.frame(a[[assay]]), info = paste(state, assay))
  }
})

test_that("KIT-260 I10: E6's managed arm through the worker's own delivery (the tab's observers, omicsManagedChainMocks)", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  fixture <- kit172Quiet(kit172Owner(file.path(root, "tab"), kit260Object("metab", 1L), "metab"))
  owner <- fixture$owner
  withr::defer(closeDIAStateAdapter(owner))
  probe <- new.env()
  notes <- new.env()
  notes$updates <- list(); notes$notes <- list(); notes$modals <- integer(); notes$choices <- list()
  guiEntryInstallBrowserDouble(notes)
  omicsManagedChainMocks(fixture, probe)
  kit172Quiet(utils::capture.output(shiny::testServer(omicsTabServer(fixture, "metabolomics", probe, TRUE), {
    session$userData$probe_owner <- owner
    inputs <- omicsTabInputs()
    inputs[c("apply_itsd", "itsd_aggregation", "norm_method", "ruv_mode", "log_offset")] <- list(TRUE, "median", "cyclicloess", "skip", 0)
    do.call(session$setInputs, inputs)
    commitWorkflowOmicsStoreStep(owner, "metab_qc_finalize", "metab_qc_complete", list(), "QC processing complete - ready for normalization")
    session$setInputs(run_normalization = 1)
    for (i in seq_len(10L)) {
      session$elapse(200)
      session$flushReact()
    }
  })))
  expect_true("metab_itsd_normalisation" %in% probe$operations)
  rows <- paste("IS_Caffeine (Alignment ID 41); Caffeine_d9 (Alignment ID 42); Caffeine-ISTD (Alignment ID 43);",
    "IS-Leucine (Alignment ID 44); PC 15:0_18:1(d7) (Alignment ID 53)")
  for (assay in c("A1", "A2")) {
    line <- sprintf("ITSD, assay %s: 5 rows used as internal standards, then removed: %s", assay, rows)
    expect_identical(sum(probe$log == line), 1L, info = paste(assay, paste(probe$log, collapse = " | ")))
  }
  text <- "Internal-standard normalisation used 5 (A1) and 5 (A2) rows as internal standards and removed them; the log names them."
  expect_identical(sum(vapply(notes$notes, function(note) identical(note$message, text), logical(1))), 1L)
})

test_that("KIT-260 VAL-276 (mirror of I05): a record made at this version replays exactly and its ITSD product carries the rows used", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  binding <- kit172Binding("metab")
  fixture <- kit260Owner(root, kit260Object("metab", 8L), "metab")
  owner <- fixture$owner
  config <- fixture$workflow$config_list
  kit172Quiet({
    runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config)
    runWorkflowOmicsReferenceLog2(binding, owner, 0, config)
    runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config)
  })
  json <- file.path(root, "record.json")
  kit172Quiet(writeWorkflowAnalysisRecord(owner, json))
  history <- unlist(owner$getHistory())
  closeDIAStateAdapter(owner)
  record <- jsonlite::read_json(json, simplifyVector = FALSE)
  itsd_step <- Filter(function(step) identical(step$state, "metab_itsd_norm"), record$steps)[[1L]]
  marker <- Filter(function(p) identical(p$name, "itsd_outcome"), itsd_step$parameters)
  expect_length(marker, 1L)
  expect_identical(length(unlist(itsd_step$outcomes$itsd_features_per_assay$A1)), 5L)
  source_dir <- fixture$source_dir
  copy <- withr::local_tempdir()
  stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
  store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
  copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
  harness <- new.env(parent = globalenv())
  sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", "metab_raw_data_s4"), api)
  inputs <- c(stats::setNames(file.path(source_dir, "import_sources", c("A1_export.txt", "A2_export.txt")), c("metab_assay_1", "metab_assay_2")),
    design_file = file.path(source_dir, "design_matrix.tab"))
  replay <- kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, fixture$paths, withr::local_tempdir(), api))
  expect_identical(replay$states, history)
  expect_identical(as.character(replay$mismatches), character())
  # The replayed ITSD product carries the rows used, as the origin's does.
  identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
  again <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
  withr::defer(closeDIAStateAdapter(again))
  spec <- readWorkflowStateMetadata(again, again$.context$snapshot$states$metab_itsd_norm)$data
  product <- hydrateSavedWorkflowEntity(again$.context$store, spec$products$metab_itsd_normalisation_science)$value
  expect_identical(product$itsd$rows_used$A1$annotation, KIT260_STANDARDS_HEAD)
})
