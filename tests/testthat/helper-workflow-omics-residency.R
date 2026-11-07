# KIT-075: the residency walk. Everything reachable from the roots (environments, lists, attributes, closures, R6 and
# reactive internals), stopping at namespaces, attached packages and the global environment; lazy and active bindings
# are not forced. Returns the paths at which an object of one of `classes` was reached.

omicsResidencyWalk <- function(roots, classes = c("LipidomicsAssayData", "MetaboliteAssayData", "TranscriptomicsAssayData"), max_nodes = 500000L,
    table_columns = character()) {
  seen <- new.env(hash = TRUE, parent = emptyenv())
  hits <- character()
  tables <- character()
  nodes <- 0L
  skip_attributes <- c("names", "class", "row.names", "dim", "dimnames", "srcref", "srcfile", "wholeSrcref", "levels", "tsp")
  stop_env <- function(env) {
    identical(env, globalenv()) || identical(env, emptyenv()) || identical(env, baseenv()) || isNamespace(env) ||
      startsWith(environmentName(env), "package:")
  }
  visit <- function(x, path) {
    nodes <<- nodes + 1L
    if (nodes > max_nodes) stop("residency walk exceeded ", max_nodes, " nodes at ", path)
    if (is.null(x)) return(invisible())
    # Rendered grid structures (grobs, gtables, units) hold drawing data, never an object or a closure.
    if (inherits(x, c("grob", "gtable", "gTree", "gList", "unit", "gpar", "viewport"))) return(invisible())
    if (isS4(x) && any(vapply(classes, function(cls) methods::is(x, cls), logical(1)))) {
      hits <<- c(hits, path)
      return(invisible())
    }
    # A plain table carrying the corpus's sample columns is an assay table in memory.
    if (is.data.frame(x) && length(table_columns) && sum(names(x) %in% table_columns) >= 2L) tables <<- c(tables, path)
    if (is.environment(x)) {
      if (stop_env(x)) return(invisible())
      key <- rlang::obj_address(x)
      if (!is.null(seen[[key]])) return(invisible())
      seen[[key]] <- TRUE
      names <- ls(x, all.names = TRUE, sorted = FALSE)
      lazy <- tryCatch(rlang::env_binding_are_lazy(x, names), error = function(e) rep(TRUE, length(names)))
      for (i in seq_along(names)) {
        name <- names[[i]]
        if (isTRUE(lazy[[i]]) || bindingIsActive(name, x)) next
        # A missing-argument binding cannot be read; it holds nothing.
        value <- tryCatch({
          value <- get0(name, envir = x, inherits = FALSE)
          force(is.null(value)) # touching the empty symbol errors
          value
        }, error = function(e) NULL)
        visit(value, paste0(path, "$", name))
      }
    } else if (is.function(x)) {
      env <- environment(x)
      if (!is.null(env)) visit(env, paste0(path, "<closure>"))
    } else if (is.list(x) && !is.data.frame(x)) {
      # Unclassed, so no class's `[[` method runs (reactivevalues refuses positional indexing).
      plain <- unclass(x)
      element_names <- names(plain)
      for (i in seq_along(plain)) visit(plain[[i]], paste0(path, "[[", if (length(element_names) && nzchar(element_names[[i]])) element_names[[i]] else i, "]]"))
    } else if (isS4(x) && !is.environment(x)) {
      for (slot in methods::slotNames(x)) visit(methods::slot(x, slot), paste0(path, "@", slot))
    }
    attrs <- attributes(x)
    for (name in setdiff(names(attrs), skip_attributes)) visit(attrs[[name]], paste0(path, "<attr:", name, ">"))
    invisible()
  }
  for (name in names(roots)) visit(roots[[name]], name)
  structure(unique(hits), tables = unique(tables))
}

#' The default value of every input a module UI renders (numeric/text value, checkbox state, the checked radio, the
#' selected option), by input id, so a `testServer` journey sets what the browser would send before an action.
omicsUiInputDefaults <- function(ui) {
  out <- list()
  attr_of <- function(tag, name) tag$attribs[[name]]
  visit <- function(tag) {
    if (inherits(tag, "shiny.tag.list") || (is.list(tag) && !inherits(tag, "shiny.tag"))) {
      for (child in tag) visit(child)
      return(invisible())
    }
    if (!inherits(tag, "shiny.tag")) return(invisible())
    id <- attr_of(tag, "id")
    if (identical(tag$name, "input")) {
      type <- attr_of(tag, "type") %||% "text"
      if (identical(type, "checkbox") && !is.null(id)) {
        out[[id]] <<- !is.null(attr_of(tag, "checked"))
      } else if (identical(type, "radio")) {
        if (!is.null(attr_of(tag, "checked"))) out[[attr_of(tag, "name")]] <<- attr_of(tag, "value")
      } else if (!is.null(id) && !is.null(attr_of(tag, "data-from"))) {
        # sliderInput: the value is `data-from` (and `data-to` for a range).
        to <- attr_of(tag, "data-to")
        out[[id]] <<- as.numeric(c(attr_of(tag, "data-from"), if (identical(attr_of(tag, "data-type"), "double")) to))
      } else if (!is.null(id) && !is.null(attr_of(tag, "value"))) {
        value <- attr_of(tag, "value")
        out[[id]] <<- if (identical(type, "number")) as.numeric(value) else value
      }
    } else if (identical(tag$name, "select") && !is.null(id)) {
      # selectInput renders its options as one HTML string.
      html <- paste(vapply(tag$children, function(child) paste(as.character(child), collapse = ""), character(1)), collapse = "")
      options <- regmatches(html, gregexpr("<option value=\"[^\"]*\"( selected)?", html))[[1]]
      chosen <- c(grep(" selected$", options, value = TRUE), options)
      if (length(chosen)) out[[id]] <<- sub("^<option value=\"([^\"]*)\".*$", "\\1", chosen[[1]])
    }
    for (child in tag$children) visit(child)
    invisible()
  }
  visit(ui)
  out
}

#' KIT-075: one omic's journey under `testServer` on the real omic module server and its session owner. A small screen
#' corpus is written by `tools/profiling/kit_omics_corpus.R`. The import runs through the tab's processing helper, as
#' the KIT-031/032 recorders call it. Every later step goes through the tabs' own inputs: the kit driver's input ids
#' (`runKitOmicsStep`), with each tab's UI defaults set first, since `testServer` renders no UI. After each step the
#' residency walk runs from the session and the bus.
#' @param omic `"lipidomics"` or `"metabolomics"`.
#' @return A named list, step to the walk's hit paths; attribute `states`, step to the head's state.
omicsResidencyJourney <- function(omic) {
  short <- if (identical(omic, "lipidomics")) "lipid" else "metab"
  repo <- normalizePath(testthat::test_path("..", ".."))
  corpus_tool <- file.path(repo, "tools", "profiling", "kit_omics_corpus.R")
  testthat::skip_if_not(file.exists(corpus_tool), "the corpus generator is in the source tree only")
  testthat::skip_if_not(workflowOwnerCapable())
  corpus_env <- new.env(parent = globalenv())
  sys.source(corpus_tool, corpus_env)
  root <- withr::local_tempdir(.local_envir = parent.frame())
  corpus_dir <- file.path(root, "corpus")
  kind <- paste0(short, "_screen")
  recipe <- corpus_env$kitOmicsRecipe(kind, biological_per_group = 3L, features = 60L)
  writer <- if (identical(omic, "lipidomics")) corpus_env$writeKitLipidCorpus else corpus_env$writeKitMetabCorpus
  invisible(utils::capture.output(writer(corpus_dir, recipe)))
  assays <- if (identical(omic, "lipidomics")) {
    c(LCMS_Pos = file.path(corpus_dir, "lipidsearch_lcms_pos.txt"), LCMS_Neg = file.path(corpus_dir, "lipidsearch_lcms_neg.txt"))
  } else {
    c(LCMS_Pos = file.path(corpus_dir, "msdial_lcms_pos.txt"), GCMS = file.path(corpus_dir, "msdial_gcms.txt"))
  }
  design_in <- utils::read.delim(file.path(corpus_dir, "design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true", .local_envir = parent.frame())
  project_dirs <- suppressMessages(setupDirectories(base_dir = file.path(root, "project", "kit075"), omic_types = omic, label = NULL,
    force = TRUE, reuse_existing = TRUE))
  # The app's summary and report helpers find the project directories in the global environment.
  old_dirs <- get0("project_dirs", envir = globalenv(), inherits = FALSE)
  assign("project_dirs", project_dirs, envir = globalenv())
  withr::defer(if (is.null(old_dirs)) rm("project_dirs", envir = globalenv()) else assign("project_dirs", old_dirs, envir = globalenv()),
    envir = parent.frame())
  server <- if (identical(omic, "lipidomics")) mod_lipidomics_server else mod_metabolomics_server
  ui <- function(name) get(paste0("mod_", short, "_", name, "_ui"), envir = asNamespace("MultiScholaR"))
  tabs <- paste0(omic, "_tabs")
  quiet <- function(...) NULL
  acc <- new.env(parent = emptyenv())
  acc$out <- list()
  acc$states <- list()
  acc$tables <- list()
  acc$bus <- list()
  shiny::testServer(server, args = list(id = paste0(omic, "_workflow"), project_dirs = project_dirs, omic_type = omic,
    experiment_label = "kit075"), {
    paths <- project_dirs[[omic]]
    set <- function(values) do.call(session$setInputs, values)
    defaults <- function(ui_name, id) set(omicsUiInputDefaults(ui(ui_name)(id)))
    walk <- function(step) {
      walked <- omicsResidencyWalk(list(session = session, workflow_data = workflow_data), table_columns = design_in$sample)
      acc$tables[[step]] <- attr(walked, "tables")
      acc$out[[step]] <- as.character(walked)
      acc$states[[step]] <- shiny::isolate(workflow_data$state_manager$current_state)
      acc$bus[[step]] <- vapply(c("data_tbl", "data_cln"), function(field) {
        value <- shiny::isolate(workflow_data[[field]])
        if (is.null(value)) "empty" else if (isSavedWorkflowReference(value)) "reference" else "resident"
      }, character(1))
    }
    # Import, through the processing the tab's process_import calls (KIT-102 W5: the shared feature-table panel's, the
    # core import then the frame's completion).
    import_config <- omicGuiConfig(omic)
    preview <- omicImportFeatureTablePreview(assays[[1]], import_config)
    vendor <- if (identical(omic, "lipidomics")) "lipidsearch" else "msdial"
    exclude <- if (identical(omic, "lipidomics")) OMICS_IMPORT_NORMALISED_COLUMNS
    processOmicImportFeatureTable(workflow_data, list(assay1_name = names(assays)[[1]], assay1_data = preview$assay_data,
      assay1_file = assays[[1]], assay2_file = assays[[2]], assay2_name = names(assays)[[2]], vendor_format = vendor,
      detected_format = preview$detected_format, id_col = preview$updates$id$selected,
      annotation_col = preview$updates$annotation$selected,
      sample_columns = resolveOmicsImportSampleColumns(preview$assay_data, vendor, "", preview$import_result, exclude),
      is_pattern = preview$updates$is_pattern %||% "", sanitize_names = FALSE, experiment_paths = paths), import_config,
      complete = function(...) completeOmicImport(workflow_data, import_config, notify = FALSE),
      showNotification = quiet, removeNotification = quiet, logInfo = quiet)
    rm(preview)
    session$flushReact()
    walk("import")
    # Design, through the builder's own inputs (the driver's actions; KIT-032's recorder for the metab builder).
    set(stats::setNames(list("design"), tabs))
    builder <- function(name) paste0("design-builder-", name)
    groups <- unique(design_in$group)
    for (i in seq_along(groups)) {
      set(stats::setNames(list(groups[[i]]), builder("new_factor")))
      set(stats::setNames(list(i), builder("add_factor")))
    }
    for (i in seq_along(groups)) {
      set(stats::setNames(list(design_in$sample[design_in$group == groups[[i]]], groups[[i]], "", "", 1L),
        builder(c("selected_runs", "factor1_select", "factor2_select", "factor3_select", "replicate_start"))))
      set(stats::setNames(list(i), builder("assign_metadata")))
    }
    clicks <- 0L
    for (group in groups) for (replicate in unique(design_in$replicates[design_in$group == group])) {
      samples <- design_in$sample[design_in$group == group & design_in$replicates == replicate]
      if (length(samples) < 2L) next
      clicks <- clicks + 1L
      set(stats::setNames(list(samples, "manual", as.integer(replicate)),
        builder(c("tech_rep_samples", "tech_rep_assignment_mode", "manual_replicate_number"))))
      set(stats::setNames(list(clicks), builder("assign_tech_reps")))
    }
    set(stats::setNames(list(groups[[2]], groups[[1]], "~ 0 + group"), builder(c("contrast_group1", "contrast_group2", "formula_string"))))
    set(stats::setNames(list(1L), builder("add_contrast")))
    set(stats::setNames(list(1L), builder("save_results")))
    walk("design")
    # QC: the intensity filter, the duplicates, the finalise.
    set(stats::setNames(list("qc"), tabs))
    defaults(paste0("qc_intensity"), "qc-intensity")
    set(list(`qc-intensity-apply_filter` = 1L))
    walk("intensity")
    defaults("qc_duplicates", "qc-duplicates")
    set(list(`qc-duplicates-detect_duplicates` = 1L))
    set(list(`qc-duplicates-resolve_duplicates` = 1L))
    walk("duplicates")
    defaults("qc_s4", "qc-s4_finalize")
    set(list(`qc-s4_finalize-finalize_qc` = 1L))
    walk("qc_finalize")
    # Normalisation: the one tab action (ITSD, log2, between-sample, RUV-III automatic), the correlation filter, the export.
    set(stats::setNames(list("norm"), tabs))
    defaults("norm", "norm")
    set(list(`norm-apply_itsd` = TRUE, `norm-norm_method` = "cyclicloess", `norm-ruv_mode` = "automatic"))
    set(list(`norm-run_normalization` = 1L))
    walk("normalisation")
    set(list(`norm-min_pearson_correlation_threshold` = 0.5))
    set(list(`norm-apply_correlation_filter` = 1L))
    walk("correlation")
    set(list(`norm-export_session` = 1L))
    walk("session_export")
    # DA: load the exported session, run.
    set(stats::setNames(list("de"), tabs))
    defaults("da", "de")
    set(list(`de-load_filtered_session` = 1L))
    walk("da_load")
    set(list(`de-formula_string` = "~ 0 + group", `de-da_q_val_thresh` = 0.05, `de-treat_lfc_cutoff` = 0))
    set(list(`de-run_da_analysis` = 1L))
    walk("da")
    # Summary: the publication copy.
    set(stats::setNames(list("summary"), tabs))
    defaults("summary", "summary")
    set(list(`summary-experiment_label` = "kit075", `summary-description` = "KIT-075 residency walk"))
    set(list(`summary-copy_to_publication` = 1L))
    walk("summary")
  })
  structure(acc$out, states = acc$states, tables = acc$tables, bus = acc$bus)
}

#' KIT-096: the transcriptomics journey under `testServer` on the real module server and its session owner: the import
#' through the tab's processing helper, the design through the builder's own inputs, the QC tab's filter and finalise
#' through theirs. A KIT-040 trx_screen corpus (3 biological samples per group, 300 genes). The walk runs after each step.
#' @return As `omicsResidencyJourney`.
trxResidencyJourney <- function() {
  repo <- normalizePath(testthat::test_path("..", ".."))
  corpus_tool <- file.path(repo, "tools", "profiling", "kit_omics_corpus.R")
  testthat::skip_if_not(file.exists(corpus_tool), "the corpus generator is in the source tree only")
  testthat::skip_if_not(workflowOwnerCapable())
  testthat::skip_if_not_installed("edgeR")
  corpus_env <- new.env(parent = globalenv())
  sys.source(corpus_tool, corpus_env)
  root <- withr::local_tempdir(.local_envir = parent.frame())
  corpus_dir <- file.path(root, "corpus")
  invisible(utils::capture.output(corpus_env$writeKitTrxCorpus(corpus_dir, corpus_env$kitOmicsRecipe("trx_screen", 3L, 300L))))
  design_in <- utils::read.delim(file.path(corpus_dir, "design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  counts_file <- file.path(corpus_dir, "counts.tsv")
  withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true", .local_envir = parent.frame())
  project_dirs <- suppressMessages(setupDirectories(base_dir = file.path(root, "project", "kit096"), omic_types = "transcriptomics",
    label = NULL, force = TRUE, reuse_existing = TRUE))
  acc <- new.env(parent = emptyenv())
  acc$out <- list()
  acc$states <- list()
  acc$tables <- list()
  acc$bus <- list()
  shiny::testServer(mod_transcriptomics_server, args = list(id = "transcriptomics_workflow", project_dirs = project_dirs,
    omic_type = "transcriptomics", experiment_label = "kit096"), {
    set <- function(values) do.call(session$setInputs, values)
    walk <- function(step) {
      walked <- omicsResidencyWalk(list(session = session, workflow_data = workflow_data), table_columns = design_in$sample)
      acc$tables[[step]] <- attr(walked, "tables")
      acc$out[[step]] <- as.character(walked)
      acc$states[[step]] <- shiny::isolate(workflow_data$state_manager$current_state)
      acc$bus[[step]] <- vapply(c("data_tbl", "data_cln"), function(field) {
        value <- shiny::isolate(workflow_data[[field]])
        if (is.null(value)) "empty" else if (isSavedWorkflowReference(value)) "reference" else "resident"
      }, character(1))
    }
    quiet <- function(...) NULL
    # KIT-102 W5: the shared feature-table panel's preview and processing, then the frame's completion.
    import_config <- omicGuiConfig("transcriptomics")
    preview <- omicImportFeatureTablePreview(counts_file, import_config)
    processOmicImportFeatureTable(workflow_data, list(assay1_name = "counts", assay1_data = preview$assay_data,
      assay1_file = counts_file, vendor_format = "trx_counts", detected_format = preview$detected_format,
      id_col = preview$updates$id$selected, annotation_col = preview$updates$annotation$selected,
      experiment_paths = project_dirs$transcriptomics), import_config,
      complete = function(...) completeOmicImport(workflow_data, import_config, notify = FALSE),
      showNotification = quiet, removeNotification = quiet, logInfo = quiet)
    rm(preview)
    session$flushReact()
    walk("import")
    set(list(transcriptomics_tabs = "design"))
    builder <- function(name) paste0("design-builder-", name)
    groups <- unique(design_in$group)
    for (i in seq_along(groups)) {
      set(stats::setNames(list(groups[[i]]), builder("new_factor")))
      set(stats::setNames(list(i), builder("add_factor")))
    }
    for (i in seq_along(groups)) {
      set(stats::setNames(list(design_in$sample[design_in$group == groups[[i]]], groups[[i]], "", "", 1L),
        builder(c("selected_runs", "factor1_select", "factor2_select", "factor3_select", "replicate_start"))))
      set(stats::setNames(list(i), builder("assign_metadata")))
    }
    clicks <- 0L
    for (set_name in unique(design_in$tech_rep_group)) {
      samples <- design_in$sample[design_in$tech_rep_group == set_name]
      clicks <- clicks + 1L
      set(stats::setNames(list(samples, "manual", as.integer(design_in$replicates[match(samples[[1L]], design_in$sample)])),
        builder(c("tech_rep_samples", "tech_rep_assignment_mode", "manual_replicate_number"))))
      set(stats::setNames(list(clicks), builder("assign_tech_reps")))
    }
    set(stats::setNames(list(groups[[2]], groups[[1]], "~ 0 + group"), builder(c("contrast_group1", "contrast_group2", "formula_string"))))
    set(stats::setNames(list(1L), builder("add_contrast")))
    set(stats::setNames(list(1L), builder("save_results")))
    walk("design")
    set(list(transcriptomics_tabs = "qc"))
    set(list(`qc-filter-apply_filter` = 1L))
    walk("filter")
    set(list(`qc-s4_finalize-finalize_qc` = 1L))
    walk("qc_finalize")
  })
  structure(acc$out, states = acc$states, tables = acc$tables, bus = acc$bus)
}
