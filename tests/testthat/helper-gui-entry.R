# KIT-101: the entry-point behaviour suite. Each lane drives one omic's public module server (`mod_<omic>_server`)
# under `shiny::testServer` from import to the summary and report, through the namespaced inputs only, with the import's
# fileInput fallback (test mode). Under this driver there is no managed application, so every step runs in process
# through the same core wrappers; the e2e browser lanes remain the managed proof.
#
# Browser double: `testServer` has no browser, so an `update*Input()` / `updateTabsetPanel()` call reaches no client.
# The lane records each such update and applies its value to the session's inputs after the current flush, as the
# browser would send it back. Notifications and modals are recorded. No MultiScholaR internal is mocked; the only other
# doubles are the e2e enrichment doubles (`install_enrichment_doubles()`, gprofiler2/clusterProfiler) and the offline
# UniProt cache below, and every lane runs with the HTTP(S) proxy pointed at a closed port, so a network call fails.
#
# Captures after every step: the head state, the state history and `tab_status`; a content digest of every committed
# state read through the owner (`getState`); every output's name with its text value, plot or image alt text, or table
# digest; the project's artifact paths, with the sha256 of the deterministic tables (volatile-field policy below); the
# notifications the step raised. A refusal step also records that the head and history are unchanged.
#
# Volatile-field policy (after helper-module-ci-oracles.R's table checksum: tables are compared as written TSV): text
# is masked for the lane's temporary root, timestamps, durations, hexadecimal (10 or more characters) and UUID
# identifiers, shiny's random tabset ids and process ids;
# sha256 is taken for .tsv/.tab/.csv/.fasta artifacts only (records, logs, RDS/qs2 stores, figures, workbooks and
# reports carry times, paths or encoder bytes and are recorded by path); store internals (`.multischolar`,
# `.workflow*`, duckdb, qs2) are not listed: the store is recorded as its committed entries, counted per masked entity
# path (KIT-101 F-13); numbers in state digests are rounded to 10 significant digits.
#
# Oracles: tests/testdata/gui-entry/<lane>.json, written with MULTISCHOLAR_GUI_ENTRY_FREEZE=1 and the source SHA.
# MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=<dir> writes each run's record there (the determinism runs).

guiEntryLaneIds <- function() c("prot_dia", "prot_dia_limpa", "prot_tmt", "prot_lfq", "prot_lfq_fragpipe", "lipid_canonical",
  "metab_lc", "metab_gc", "metab_combined", "trx_screen")

#' The frozen oracles (MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR redirects them: the lifecycle proof's provisional oracles).
guiEntryOracleDir <- function() {
  dir <- Sys.getenv("MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR")
  if (nzchar(dir)) dir else testthat::test_path("..", "testdata", "gui-entry")
}

guiEntryFixtureDir <- function(lane) normalizePath(testthat::test_path("..", "testdata", "e2e", lane$fixture_dir), mustWork = TRUE)

#' The lane's manifest entry (tests/testdata/e2e/manifest.json), read without the e2e helper's validation.
guiEntryLane <- function(lane_id) {
  raw <- jsonlite::read_json(testthat::test_path("..", "testdata", "e2e", "manifest.json"), simplifyVector = FALSE)
  lanes <- stats::setNames(raw$lanes, vapply(raw$lanes, `[[`, "", "lane_id"))
  if (is.null(lanes[[lane_id]])) stop("no e2e manifest lane: ", lane_id)
  lanes[[lane_id]]
}

# --------------------------------------------------------------------------------------------------------------------
# Normalisation and digests.
# --------------------------------------------------------------------------------------------------------------------

guiEntryMask <- function(text, roots = character()) {
  text <- paste(as.character(text), collapse = "\n")
  for (root in unique(roots[nzchar(roots)])) text <- gsub(root, "<ROOT>", text, fixed = TRUE)
  package_root <- system.file(package = "MultiScholaR")
  if (nzchar(package_root)) text <- gsub(package_root, "<PKG>", text, fixed = TRUE)
  text <- gsub(normalizePath(tempdir(), mustWork = FALSE), "<TMP>", text, fixed = TRUE)
  text <- gsub("/tmp/Rtmp[A-Za-z0-9]+", "<TMP>", text)
  text <- gsub("[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}", "<UUID>", text)
  # Hexadecimal identifiers of 10 or more characters: tempfile suffixes in file names, htmlwidget and plotly ids.
  text <- gsub("(?<![0-9a-zA-Z])[0-9a-f]{10,}(?![0-9a-zA-Z])", "<HEX>", text, perl = TRUE)
  # shiny's random tabset ids (tabsetPanel / navlistPanel inside rendered UI).
  text <- gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"<ID>\"", text)
  text <- gsub("tab-[0-9]+-([0-9]+)", "tab-<ID>-\\1", text)
  text <- gsub("(?<![0-9])[0-9]{4}-[0-9]{2}-[0-9]{2}[ T_][0-9]{2}[:-]?[0-9]{2}([:-]?[0-9]{2}(\\.[0-9]+)?)?( ?[A-Z]{3,4}\\b|[+-][0-9]{2}:?[0-9]{2}|Z)?",
    "<TIME>", text, perl = TRUE)
  text <- gsub("(?<![0-9])[0-9]{8}_[0-9]{6}(?![0-9])", "<TIME>", text, perl = TRUE)
  text <- gsub("(?<![0-9])[0-9]{2}:[0-9]{2}:[0-9]{2}(\\.[0-9]+)?(?![0-9])", "<TIME>", text, perl = TRUE)
  text <- gsub("(?<![0-9])[0-9]{4}-[0-9]{2}-[0-9]{2}(?![0-9])", "<DATE>", text, perl = TRUE)
  text <- gsub("[0-9]+(\\.[0-9]+)? ?(secs|seconds|sec|s|ms|mins|minutes)\\b", "<DURATION>", text, perl = TRUE)
  text <- gsub("(pid|PID|process) ?[0-9]+", "\\1 <PID>", text)
  text
}

#' A content digest that ignores environments and functions and rounds numbers, so that it depends on the data only.
guiEntryCanonical <- function(x, depth = 0L) {
  if (depth > 12L) return("<deep>")
  if (is.null(x)) return(NULL)
  if (is.environment(x) || is.function(x) || inherits(x, "R6")) return(paste0("<", class(x)[[1L]], ">"))
  if (isS4(x)) {
    slots <- methods::slotNames(x)
    return(c(list(`<class>` = class(x)[[1L]]), stats::setNames(lapply(slots, function(s) guiEntryCanonical(methods::slot(x, s), depth + 1L)), slots)))
  }
  if (is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    return(list(`<names>` = names(x), `<rows>` = nrow(x), `<cols>` = lapply(x, function(col) guiEntryCanonical(col, depth + 1L))))
  }
  if (is.factor(x)) return(list(`<levels>` = levels(x), `<codes>` = as.integer(x)))
  if (is.list(x)) return(lapply(unclass(x), function(el) guiEntryCanonical(el, depth + 1L)))
  if (is.numeric(x)) {
    out <- signif(as.numeric(x), 10L)
    if (!is.null(dim(x))) out <- list(`<dim>` = dim(x), `<dimnames>` = dimnames(x), `<values>` = out)
    else if (!is.null(names(x))) out <- list(`<names>` = names(x), `<values>` = out)
    return(out)
  }
  if (inherits(x, c("POSIXt", "Date", "difftime"))) return("<time>")
  if (is.character(x) || is.logical(x)) {
    if (!is.null(dim(x))) return(list(`<dim>` = dim(x), `<values>` = as.vector(x)))
    return(x)
  }
  paste0("<", class(x)[[1L]], ">")
}

guiEntryDigest <- function(x) digest::digest(guiEntryCanonical(x), algo = "sha256")

# --------------------------------------------------------------------------------------------------------------------
# Captures.
# --------------------------------------------------------------------------------------------------------------------

#' Every output of the session by name: text values, plot/image alt text, table digests; download outputs by name only.
guiEntryOutputs <- function(root, roots) {
  private <- root$.__enclos_env__$private
  names <- sort(names(private$outs))
  downloads <- private$file_generators$keys()
  out <- list()
  for (name in names) {
    # A download output renders its file on read; it is recorded by name only.
    if (name %in% downloads || any(endsWith(downloads, paste0("-", name)))) { out[[name]] <- "<download>"; next }
    value <- tryCatch(root$getOutput(name), error = function(e) {
      if (inherits(e, c("shiny.silent.error", "shiny.output.cancel"))) "<req>" else paste("<error>", guiEntryMask(conditionMessage(e), roots))
    })
    out[[name]] <- guiEntryOutputValue(value, roots)
  }
  out
}

guiEntryOutputValue <- function(value, roots) {
  if (is.null(value)) return("<NULL>")
  if (is.character(value)) {
    text <- guiEntryMask(paste(value, collapse = "\n"), roots)
    # Long text (rendered tables, logs) is held as its digest.
    return(if (nchar(text) > 1000L) paste("<text>", digest::digest(text, algo = "sha256", serialize = FALSE)) else text)
  }
  if (is.list(value) && !is.null(value$src)) {
    kind <- if (startsWith(as.character(value$src), "data:image")) "plot" else "image"
    return(paste0("<", kind, "> alt=", guiEntryMask(value$alt %||% "", roots)))
  }
  if (is.list(value) && !is.null(value$html)) {
    # renderUI: the HTML with the masks; the input ids it holds are what the browser would offer.
    return(paste("<ui>", digest::digest(guiEntryMask(as.character(value$html), roots), algo = "sha256")))
  }
  if (is.list(value) && !is.null(value$x)) {
    data <- value$x$data %||% value$x
    return(paste("<widget>", digest::digest(guiEntryCanonical(data), algo = "sha256")))
  }
  if (is.list(value) || is.data.frame(value)) return(paste("<table>", guiEntryDigest(value)))
  paste0("<", class(value)[[1L]], ">")
}

#' The project's artifacts: relative paths; sha256 of deterministic tables; the owner's store by its committed entries.
#' KIT-101 revision 2.2 (F-13): the store's raw file count varied by one to three between runs (files the store stages
#' and releases); the store is now captured as its committed entries only, an entity directory with its `commit.json`,
#' counted per masked entity path (store and request ids masked), so the capture is a deterministic list.
guiEntryArtifacts <- function(base_dir, roots) {
  files <- sort(list.files(base_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE))
  store <- grepl("(^|/)\\.multischolar(/|$)|(^|/)\\.workflow|\\.duckdb|\\.qs2$|\\.wal$", files)
  commits <- grep("(^|/)\\.workflow-entities/entities/[^/]+/commit\\.json$", files[store], value = TRUE)
  entries <- table(vapply(dirname(commits), guiEntryMask, "", roots = roots, USE.NAMES = FALSE))
  out <- list(`<store entries>` = as.list(stats::setNames(as.integer(entries), names(entries))))
  for (file in files[!store]) {
    key <- guiEntryMask(file, roots)
    out[[key]] <- if (grepl("\\.(tsv|tab|csv|fasta)$", file, ignore.case = TRUE)) {
      text <- readLines(file.path(base_dir, file), warn = FALSE)
      digest::digest(guiEntryMask(text, roots), algo = "sha256", serialize = FALSE)
    } else "<present>"
  }
  out
}

guiEntryHistory <- function(manager) {
  if (is.null(manager)) return(character())
  tryCatch(as.character(unlist(manager$.context$snapshot$state_history)), error = function(e) character())
}

guiEntryStateDigests <- function(manager, cache) {
  history <- guiEntryHistory(manager)
  out <- list()
  for (name in unique(history)) {
    out[[name]] <- tryCatch(guiEntryDigest(withWorkflowStateScope(manager, function() manager$getState(name))),
      error = function(e) paste("<error>", conditionMessage(e)))
  }
  out
}

# --------------------------------------------------------------------------------------------------------------------
# Doubles.
# --------------------------------------------------------------------------------------------------------------------

#' The offline UniProt double (KIT-101 W2.2): the design step's `getUniprotAnnotationsFull` reads
#' `<results_dir>/cache/uniprot_annotations/uniprot_annotations.RDS` before it would go to the network
#' (`getUniprotAnnotations`, func_prot_annotation_uniprot_helpers.R); the committed cache for the fixture's accessions is
#' placed there before the lane starts.
guiEntryInstallUniprotCache <- function(experiment_paths) {
  cache <- readRDS(testthat::test_path("..", "testdata", "gui-entry", "uniprot", "uniprot_annotations.RDS"))
  target <- file.path(experiment_paths$results_dir, "cache", "uniprot_annotations")
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, file.path(target, "uniprot_annotations.RDS"))
  invisible(target)
}

#' The browser double and the notification log, scoped to `env`.
guiEntryInstallBrowserDouble <- function(log, env = parent.frame()) {
  queue_update <- function(session, input_id, value, choices = NULL) {
    if (is.null(session)) return(invisible(NULL))
    if (!is.null(choices)) log$choices[[session$ns(input_id)]] <- as.character(unlist(choices, use.names = FALSE))
    if (is.null(value)) return(invisible(NULL))
    log$updates[[length(log$updates) + 1L]] <- list(id = session$ns(input_id), value = value)
    invisible(NULL)
  }
  notify <- function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL, type = "default", session = NULL) {
    log$notes[[length(log$notes) + 1L]] <- list(type = type, message = paste(as.character(ui), collapse = " "))
    invisible(id %||% "notification")
  }
  quiet <- function(...) invisible(NULL)
  bindings <- list(
    updateSelectInput = function(session = shiny::getDefaultReactiveDomain(), inputId, label = NULL, choices = NULL, selected = NULL) queue_update(session, inputId, selected, choices),
    updateSelectizeInput = function(session = shiny::getDefaultReactiveDomain(), inputId, label = NULL, choices = NULL, selected = NULL, options = list(), server = FALSE) queue_update(session, inputId, selected, choices),
    updateRadioButtons = function(session = shiny::getDefaultReactiveDomain(), inputId, label = NULL, choices = NULL, selected = NULL, inline = FALSE, choiceNames = NULL, choiceValues = NULL) queue_update(session, inputId, selected, choices %||% choiceValues),
    updateTextInput = function(session = shiny::getDefaultReactiveDomain(), inputId, label = NULL, value = NULL, placeholder = NULL) queue_update(session, inputId, value),
    updateTextAreaInput = function(session = shiny::getDefaultReactiveDomain(), inputId, label = NULL, value = NULL, placeholder = NULL) queue_update(session, inputId, value),
    updateNumericInput = function(session = shiny::getDefaultReactiveDomain(), inputId, label = NULL, value = NULL, min = NULL, max = NULL, step = NULL) queue_update(session, inputId, value),
    updateCheckboxInput = function(session = shiny::getDefaultReactiveDomain(), inputId, label = NULL, value = NULL) queue_update(session, inputId, value),
    updateSliderInput = function(session = shiny::getDefaultReactiveDomain(), inputId, label = NULL, value = NULL, min = NULL, max = NULL, step = NULL, timeFormat = NULL, timezone = NULL) queue_update(session, inputId, value),
    updateTabsetPanel = function(session = shiny::getDefaultReactiveDomain(), inputId, selected = NULL) queue_update(session, inputId, selected),
    showNotification = notify, removeNotification = quiet,
    showModal = function(ui, session = shiny::getDefaultReactiveDomain()) { log$modals <- c(log$modals, 1L); invisible(NULL) },
    removeModal = quiet)
  do.call(testthat::local_mocked_bindings, c(bindings, list(.package = "shiny", .env = env)))
  # The package's own imports of these functions (importFrom(shiny, ...)) are bound in its imports environment.
  imported <- intersect(names(bindings), unlist(getNamespaceImports("MultiScholaR")$shiny, use.names = FALSE))
  if (length(imported)) do.call(testthat::local_mocked_bindings, c(bindings[imported], list(.package = "MultiScholaR", .env = env)))
  invisible(log)
}

# --------------------------------------------------------------------------------------------------------------------
# The lane runner.
# --------------------------------------------------------------------------------------------------------------------

guiEntryServerSpec <- function(omic) {
  switch(omic,
    proteomics = list(server = mod_proteomics_server, id = "proteomics_workflow", tabs = "workflow_tabs"),
    lipidomics = list(server = mod_lipidomics_server, id = "lipidomics_workflow", tabs = "lipidomics_tabs"),
    metabolomics = list(server = mod_metabolomics_server, id = "metabolomics_workflow", tabs = "metabolomics_tabs"),
    transcriptomics = list(server = mod_transcriptomics_server, id = "transcriptomics_workflow", tabs = "transcriptomics_tabs"))
}

#' Run one lane. Returns the record: lane, source SHA, the steps' captures in order, the seconds taken.
#' @param fault Optional `list(step = <id>, fn = function() ...)` run in place of that step's actions (the failure check).
guiEntryRunLane <- function(lane_id, fault = NULL) {
  if (!workflowOwnerCapable()) stop("the session owner is not available (qs2/duckdb): the entry-point lane cannot run")
  started <- Sys.time()
  lane <- guiEntryLane(lane_id)
  omic <- lane$omic_type
  spec <- guiEntryServerSpec(omic)
  steps <- guiEntrySteps(lane)
  root <- withr::local_tempdir(pattern = paste0("gui-entry-", lane_id, "-"))
  root <- normalizePath(root)
  withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true", http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
    HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "")
  withr::local_options(warn = 1)
  withr::local_seed(20261107L)
  project_dirs <- suppressMessages(setupDirectories(base_dir = file.path(root, "project", lane_id), omic_types = omic, label = NULL,
    force = TRUE, reuse_existing = TRUE))
  paths <- project_dirs[[omic]]
  old_dirs <- get0("project_dirs", envir = globalenv(), inherits = FALSE)
  assign("project_dirs", project_dirs, envir = globalenv())
  withr::defer(if (is.null(old_dirs)) rm("project_dirs", envir = globalenv()) else assign("project_dirs", old_dirs, envir = globalenv()))
  if (identical(omic, "proteomics")) {
    guiEntryInstallUniprotCache(paths)
    install_enrichment_doubles()
  }
  log <- new.env(parent = emptyenv())
  log$updates <- list(); log$notes <- list(); log$modals <- integer(); log$choices <- list()
  guiEntryInstallBrowserDouble(log)
  inputs <- list(dir = guiEntryFixtureDir(lane), scratch = file.path(root, "inputs"), root = root, lane = lane, paths = paths)
  dir.create(inputs$scratch, showWarnings = FALSE)
  acc <- new.env(parent = emptyenv())
  acc$steps <- list()
  roots <- c(root, paths$base_dir)
  suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(spec$server, args = list(id = spec$id, project_dirs = project_dirs,
    omic_type = omic, experiment_label = lane_id), {
    rootSession <- session$rootScope()
    prefix <- paste0(spec$id, "-")
    apply_updates <- function() {
      for (i in seq_len(50L)) {
        if (!length(log$updates)) return(invisible())
        pending <- log$updates; log$updates <- list()
        values <- stats::setNames(lapply(pending, `[[`, "value"), vapply(pending, `[[`, "", "id"))
        values <- values[!duplicated(names(values), fromLast = TRUE)]
        do.call(rootSession$setInputs, values)
      }
    }
    api <- list(
      set = function(values) {
        names(values) <- paste0(prefix, names(values))
        do.call(rootSession$setInputs, values)
        apply_updates()
      },
      defaults = function(ui_fn, id) {
        values <- omicsUiInputDefaults(ui_fn(id))
        if (length(values)) { names(values) <- paste0(prefix, names(values)); do.call(rootSession$setInputs, values); apply_updates() }
      },
      upload = function(id, path) {
        copy <- file.path(inputs$scratch, basename(path))
        file.copy(path, copy, overwrite = TRUE)
        value <- data.frame(name = basename(path), size = file.size(copy), type = "", datapath = copy, stringsAsFactors = FALSE)
        values <- stats::setNames(list(value), paste0(prefix, id))
        do.call(rootSession$setInputs, values)
        apply_updates()
      },
      input = function(id) shiny::isolate(rootSession$input[[paste0(prefix, id)]]),
      choices = function(id) log$choices[[paste0(prefix, id)]] %||% character(),
      modals = function() log$modals,
      workflow_data = workflow_data, session = session, inputs = inputs, tabs = spec$tabs)
    manager <- function() shiny::isolate(workflow_data$state_manager)
    capture <- function(step, refusal_before = NULL) {
      session$flushReact()
      head <- tryCatch(shiny::isolate(manager()$current_state), error = function(e) NULL) %||% "<none>"
      history <- guiEntryHistory(manager())
      record <- list(
        head = head, history = history,
        tab_status = lapply(shiny::isolate(workflow_data$tab_status), as.character),
        states = guiEntryStateDigests(manager()),
        outputs = guiEntryOutputs(rootSession, roots),
        artifacts = guiEntryArtifacts(paths$base_dir, roots),
        notifications = lapply(log$notes, function(note) list(type = note$type, message = guiEntryMask(note$message, roots))))
      if (!is.null(refusal_before)) {
        record$refusal_check <- list(head_unchanged = identical(refusal_before$head, head), history_unchanged = identical(refusal_before$history, history),
          refused = length(log$notes) > 0L)
      }
      log$notes <- list()
      acc$steps[[step$id]] <- c(list(step = step$id, tab = step$tab, kind = step$kind), record)
    }
    for (step in steps) {
      before <- if (!identical(step$kind, "action")) list(head = tryCatch(shiny::isolate(manager()$current_state), error = function(e) NULL) %||% "<none>",
        history = guiEntryHistory(manager())) else NULL
      if (!is.null(step$tab)) api$set(stats::setNames(list(step$tab), spec$tabs))
      if (!is.null(fault) && identical(fault$step, step$id)) fault$fn() else step$run(api)
      capture(step, before)
    }
    shiny::isolate(closeDIAStateAdapter(workflow_data$state_manager))
  }))))
  record <- list(lane = lane_id, omic = omic, steps = acc$steps, seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1))
  record_dir <- Sys.getenv("MULTISCHOLAR_GUI_ENTRY_RECORD_DIR")
  if (nzchar(record_dir)) {
    dir.create(record_dir, recursive = TRUE, showWarnings = FALSE)
    guiEntryWriteJson(record, file.path(record_dir, paste0(lane_id, ".json")))
  }
  record
}

guiEntryWriteJson <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, null = "null", na = "string", digits = NA, pretty = TRUE)
}

#' The record as an oracle holds it (JSON round trip), without the timing.
guiEntryComparable <- function(record) {
  x <- jsonlite::fromJSON(jsonlite::toJSON(record$steps, auto_unbox = TRUE, null = "null", na = "string", digits = NA), simplifyVector = FALSE)
  x
}

guiEntryOraclePath <- function(lane_id) file.path(guiEntryOracleDir(), paste0(lane_id, ".json"))

guiEntrySourceSha <- function() {
  repo <- testthat::test_path("..", "..")
  sha <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)), error = function(e) character())
  if (length(sha) && nzchar(sha[[1L]])) sha[[1L]] else Sys.getenv("MULTISCHOLAR_GUI_ENTRY_SOURCE_SHA", "unknown")
}

#' Run a lane and compare it with its frozen oracle, step by step (or freeze it under MULTISCHOLAR_GUI_ENTRY_FREEZE=1).
#' A missing oracle fails; there is no skip.
guiEntryExpectLane <- function(lane_id, ceiling = 900) {
  record <- guiEntryRunLane(lane_id)
  testthat::expect_lt(record$seconds, ceiling)
  steps <- guiEntryComparable(record)
  testthat::expect_identical(names(steps), vapply(guiEntrySteps(guiEntryLane(lane_id)), `[[`, "", "id"))
  # A refusal leaves the head and the history as they were, with its message; a silent refusal raises none.
  for (step in steps) if (!identical(step$kind, "action")) {
    check <- step$refusal_check
    testthat::expect_true(isTRUE(check$head_unchanged) && isTRUE(check$history_unchanged), info = paste(lane_id, step$step, "head and history unchanged"))
    testthat::expect_identical(isTRUE(check$refused), identical(step$kind, "refusal"), info = paste(lane_id, step$step, "refusal message"))
  }
  path <- guiEntryOraclePath(lane_id)
  if (identical(Sys.getenv("MULTISCHOLAR_GUI_ENTRY_FREEZE"), "1")) {
    guiEntryWriteJson(list(lane = lane_id, source_sha = guiEntrySourceSha(), load_mode = "load_all", r_version = R.version.string,
      frozen = format(Sys.Date()), steps = steps), path)
    testthat::succeed(paste("froze", lane_id))
    return(invisible(record))
  }
  if (!file.exists(path)) {
    testthat::fail(paste("no frozen oracle for lane", lane_id, "at", path))
    return(invisible(record))
  }
  oracle <- jsonlite::read_json(path, simplifyVector = FALSE)
  testthat::expect_identical(names(steps), names(oracle$steps), info = paste(lane_id, "step order"))
  for (id in names(oracle$steps)) {
    testthat::expect_identical(steps[[id]], oracle$steps[[id]], info = paste0(lane_id, " at step '", id, "' (oracle ", oracle$source_sha, ")"))
  }
  invisible(record)
}

# --------------------------------------------------------------------------------------------------------------------
# The lanes' steps. Each step: id, the tab it opens (the tabset value), whether it is a refusal, and its actions through
# the namespaced inputs (ids relative to the omic's module).
# --------------------------------------------------------------------------------------------------------------------

guiEntrySteps <- function(lane) {
  switch(lane$omic_type,
    transcriptomics = guiEntryTrxSteps(lane),
    lipidomics = guiEntryOmicsSteps(lane, "lipid"),
    metabolomics = guiEntryOmicsSteps(lane, "metab"),
    proteomics = guiEntryProtSteps(lane))
}

guiEntryStep <- function(id, tab, run, refusal = c("none", "message", "silent")) {
  refusal <- match.arg(refusal)
  list(id = id, tab = tab, run = run, kind = switch(refusal, none = "action", message = "refusal", silent = "silent_refusal"))
}

#' The design builder through its own inputs (the KIT-075/096 journeys): factors, metadata per group, technical
#' replicate sets when the design has them, one contrast, save.
guiEntryDesignActions <- function(api, design, contrast, prefix = "design-builder-", tech_rep = TRUE) {
  builder <- function(name) paste0(prefix, name)
  groups <- unique(design$group)
  for (i in seq_along(groups)) {
    api$set(stats::setNames(list(groups[[i]]), builder("new_factor")))
    api$set(stats::setNames(list(i), builder("add_factor")))
  }
  for (i in seq_along(groups)) {
    api$set(stats::setNames(list(design$sample[design$group == groups[[i]]], groups[[i]], "", "", 1L),
      builder(c("selected_runs", "factor1_select", "factor2_select", "factor3_select", "replicate_start"))))
    api$set(stats::setNames(list(i), builder("assign_metadata")))
  }
  if (tech_rep && !is.null(design$tech_rep_group)) {
    clicks <- 0L
    for (set_name in unique(design$tech_rep_group)) {
      samples <- design$sample[design$tech_rep_group == set_name]
      if (length(samples) < 2L) next
      clicks <- clicks + 1L
      api$set(stats::setNames(list(samples, "manual", as.integer(design$replicates[match(samples[[1L]], design$sample)])),
        builder(c("tech_rep_samples", "tech_rep_assignment_mode", "manual_replicate_number"))))
      api$set(stats::setNames(list(clicks), builder("assign_tech_reps")))
    }
  }
  api$set(stats::setNames(list(contrast[[1L]], contrast[[2L]], "~ 0 + group"), builder(c("contrast_group1", "contrast_group2", "formula_string"))))
  api$set(stats::setNames(list(1L), builder("add_contrast")))
  api$set(stats::setNames(list(1L), builder("save_results")))
}

guiEntryReadDesign <- function(api) {
  utils::read.delim(file.path(api$inputs$dir, "design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
}

guiEntryTrxSteps <- function(lane) {
  profile <- trxOmicProfile()
  ui <- function(fn) function(id) fn(id, profile = profile)
  list(
    guiEntryStep("import_refusal", "import", refusal = "message", function(api) {
      # A count table with one negative count (KIT-040's counts_negative case) is refused with no state.
      lines <- readLines(file.path(api$inputs$dir, "counts.tsv"))
      cells <- strsplit(lines[[2L]], "\t", fixed = TRUE)[[1L]]
      cells[[3L]] <- "-1"
      lines[[2L]] <- paste(cells, collapse = "\t")
      bad <- file.path(api$inputs$scratch, "negative", "counts.tsv")
      dir.create(dirname(bad), showWarnings = FALSE)
      writeLines(lines, bad)
      api$defaults(ui(mod_lipid_import_ui), "import")
      api$upload("import-assay1_file_std", bad)
      api$set(list(`import-process_import` = 1L))
    }),
    guiEntryStep("import", "import", function(api) {
      api$upload("import-assay1_file_std", file.path(api$inputs$dir, "counts.tsv"))
      api$set(list(`import-process_import` = 2L))
    }),
    guiEntryStep("design_refusal", "design", refusal = "message", function(api) {
      api$set(list(`design-builder-save_results` = 1L))
    }),
    guiEntryStep("design", "design", function(api) {
      design <- guiEntryReadDesign(api)
      guiEntryDesignActions(api, design, c("Sera", "RPMI"))
    }),
    guiEntryStep("qc_filter", "qc", function(api) {
      api$defaults(mod_trx_qc_filter_ui, "qc-filter")
      api$set(list(`qc-filter-apply_filter` = 1L))
    }),
    guiEntryStep("qc_finalize", "qc", function(api) {
      api$defaults(ui(mod_lipid_qc_s4_ui), "qc-s4_finalize")
      api$set(list(`qc-s4_finalize-finalize_qc` = 1L))
    }),
    guiEntryStep("da_refusal", "de", refusal = "message", function(api) {
      # DA before the normalisation: the head has no RUVs factors.
      api$defaults(ui(mod_lipid_da_ui), "de")
      api$set(list(`de-run_da_analysis` = 1L))
    }),
    guiEntryStep("norm_refusal", "norm", refusal = "message", function(api) {
      # A manual k above the cap (MCI-030.2), refused with its class.
      api$defaults(ui(mod_lipid_norm_ui), "norm")
      api$set(list(`norm-ruv_mode` = "manual", `norm-ruv_k` = 50L))
      api$set(list(`norm-run_normalization` = 1L))
    }),
    guiEntryStep("norm", "norm", function(api) {
      api$defaults(ui(mod_lipid_norm_ui), "norm")
      api$set(list(`norm-run_normalization` = 2L))
    }),
    guiEntryStep("da", "de", function(api) {
      api$defaults(ui(mod_lipid_da_ui), "de")
      api$set(list(`de-run_da_analysis` = 2L))
    }),
    guiEntryStep("summary_refusal", "summary", refusal = "silent", function(api) {
      # The report before the publication copy: the button is gated (req) and says nothing.
      api$defaults(ui(mod_lipid_summary_ui), "summary")
      api$set(list(`summary-experiment_label` = "trx_screen"))
      api$set(list(`summary-generate_report` = 1L))
    }),
    guiEntryStep("summary_ruviii", "summary", function(api) {
      api$set(list(`summary-ruviii_m_source` = "group"))
      api$set(list(`summary-run_ruviii_export` = 1L))
    }),
    guiEntryStep("summary_save", "summary", function(api) {
      api$set(list(`summary-experiment_label` = "trx_screen", `summary-description` = "KIT-101 entry-point lane"))
      api$set(list(`summary-save_workflow_args` = 1L))
    }),
    guiEntryStep("summary_copy", "summary", function(api) api$set(list(`summary-copy_to_publication` = 1L))),
    guiEntryStep("summary_report", "summary", function(api) api$set(list(`summary-generate_report` = 2L))),
    guiEntryStep("summary_export", "summary", function(api) api$set(list(`summary-export_session_state` = 1L))))
}

#' The lipidomics and metabolomics lanes: the e2e lane's import inputs (helper-e2e-browser.R e2e_upload_lane_inputs),
#' then the KIT-075 residency journey's steps through the tabs' own inputs (intensity, duplicates, ITSD analysis and
#' finalise; normalisation, correlation filter, session export; DA load and run; the summary).
guiEntryOmicsSteps <- function(lane, short) {
  ui <- function(name) get(paste0("mod_", short, "_", name, "_ui"), envir = asNamespace("MultiScholaR"))
  assays <- unlist(lane$assays)
  import_inputs <- function(api) {
    api$defaults(ui("import"), "import")
    api$set(list(`import-assay1_name` = assays[[1L]], `import-vendor_format` = lane$import_tool))
    api$upload("import-assay1_file_std", file.path(api$inputs$dir, lane$seed_file))
    if (!is.null(lane$assay2_file) && length(assays) >= 2L) {
      api$set(list(`import-assay2_name` = assays[[2L]]))
      api$upload("import-assay2_file_std", file.path(api$inputs$dir, lane$assay2_file))
    }
    if (identical(lane$import_tool, "custom")) {
      api$set(list(`import-metabolite_id_col_custom` = "Feature.Name", `import-annotation_col_custom` = "",
        `import-sample_cols_pattern` = "^(WT|KO)_", `import-is_pattern` = "^ITSD|_d[0-9]+$"))
    }
    api$set(list(`import-sanitize_names` = FALSE))
  }
  list(
    guiEntryStep("import_refusal", "import", refusal = "message", function(api) {
      # The lane's files with an empty primary assay name.
      import_inputs(api)
      api$set(list(`import-assay1_name` = ""))
      api$set(list(`import-process_import` = 1L))
    }),
    guiEntryStep("import", "import", function(api) {
      import_inputs(api)
      api$set(list(`import-process_import` = 2L))
    }),
    guiEntryStep("design_refusal", "design", refusal = "message", function(api) {
      api$set(list(`design-builder-save_results` = 1L))
    }),
    guiEntryStep("design", "design", function(api) guiEntryDesignActions(api, guiEntryReadDesign(api), c("KO", "WT"))),
    guiEntryStep("qc_intensity", "qc", function(api) {
      api$defaults(ui("qc_intensity"), "qc-intensity")
      api$set(list(`qc-intensity-apply_filter` = 1L))
    }),
    guiEntryStep("qc_duplicates", "qc", function(api) {
      api$defaults(ui("qc_duplicates"), "qc-duplicates")
      api$set(list(`qc-duplicates-detect_duplicates` = 1L))
      api$set(list(`qc-duplicates-resolve_duplicates` = 1L))
    }),
    guiEntryStep("qc_itsd_refusal", "qc", refusal = "message", function(api) {
      # The internal-standard analysis on assays with no standard matching the pattern.
      api$defaults(ui("qc_itsd"), "qc-itsd")
      api$set(list(`qc-itsd-analyze_is` = 1L))
    }),
    guiEntryStep("qc_finalize", "qc", function(api) {
      api$defaults(ui("qc_s4"), "qc-s4_finalize")
      api$set(list(`qc-s4_finalize-finalize_qc` = 1L))
    }),
    guiEntryStep("da_refusal", "de", refusal = "message", function(api) {
      # DA before any filtered session was exported or loaded.
      api$defaults(ui("da"), "de")
      api$set(list(`de-load_filtered_session` = 1L))
    }),
    guiEntryStep("norm_refusal", "norm", refusal = "message", function(api) {
      # The session export before the normalisation.
      api$defaults(ui("norm"), "norm")
      api$set(list(`norm-export_session` = 1L))
    }),
    guiEntryStep("norm", "norm", function(api) {
      api$defaults(ui("norm"), "norm")
      api$set(list(`norm-apply_itsd` = TRUE, `norm-norm_method` = "cyclicloess", `norm-ruv_mode` = "skip"))
      api$set(list(`norm-run_normalization` = 1L))
    }),
    guiEntryStep("correlation", "norm", function(api) {
      api$set(list(`norm-min_pearson_correlation_threshold` = 0.5))
      api$set(list(`norm-apply_correlation_filter` = 1L))
    }),
    guiEntryStep("session_export", "norm", function(api) api$set(list(`norm-export_session` = 2L))),
    guiEntryStep("da_load", "de", function(api) api$set(list(`de-load_filtered_session` = 2L))),
    guiEntryStep("da", "de", function(api) {
      api$set(list(`de-formula_string` = "~ 0 + group", `de-da_q_val_thresh` = 0.05, `de-treat_lfc_cutoff` = 0))
      api$set(list(`de-run_da_analysis` = 1L))
    }),
    guiEntryStep("summary_refusal", "summary", refusal = "silent", function(api) {
      # The report before the publication copy: the button is gated (req) and says nothing.
      api$defaults(ui("summary"), "summary")
      api$set(list(`summary-experiment_label` = lane$lane_id))
      api$set(list(`summary-generate_report` = 1L))
    }),
    guiEntryStep("summary_save", "summary", function(api) {
      api$set(list(`summary-experiment_label` = lane$lane_id, `summary-description` = "KIT-101 entry-point lane"))
      api$set(list(`summary-save_workflow_args` = 1L))
    }),
    guiEntryStep("summary_copy", "summary", function(api) api$set(list(`summary-copy_to_publication` = 1L))),
    guiEntryStep("summary_report", "summary", function(api) api$set(list(`summary-generate_report` = 2L))),
    guiEntryStep("summary_export", "summary", function(api) api$set(list(`summary-export_session_state` = 1L))))
}

#' The action buttons a sub-tab UI renders (`<button id=... class="... action-button">`).
guiEntryActionIds <- function(ui) {
  html <- as.character(htmltools::renderTags(ui)$html)
  tags <- regmatches(html, gregexpr("<button[^>]*>", html))[[1L]]
  tags <- grep("action-button", tags, value = TRUE, fixed = TRUE)
  sub("^.* id=\"([^\"]+)\".*$", "\\1", tags)
}

#' The small e2e fixtures need the e2e lanes' permissive QC settings (helper-e2e-browser.R e2e_run_prot_dia_qc); each
#' is set only where the sub-tab renders that input.
guiEntryProtQcSettings <- function() {
  list(use_strict_mode = FALSE, min_reps_per_group = 1, min_groups = 1, intensity_cutoff_percentile = 0.1,
    proteins_intensity_cutoff_percentile = 0.1, min_peptides_per_protein = 1, min_peptidoforms_per_protein = 1,
    min_peptides_per_sample = 1, proportion_missing_values = 1.0, parallel_cores = 1)
}

#' The proteomics lanes: the e2e lane's import inputs, the design builder, every QC sub-tab of the route's spec tables
#' (`getProtQcPeptideModuleSpecs` / `getProtQcProteinModuleSpecs`, the production configuration) in order, the
#' normalisation (RUV skipped on these six-sample fixtures, as the e2e lanes run it) with the correlation filter (and the
#' limpa route's RUV covariates), the session export, the DA, the enrichment (prot_dia, through the doubles) and the
#' summary.
guiEntryProtSteps <- function(lane) {
  route <- if (isTRUE(lane$use_limpa)) "DIA_LIMPA" else lane$workflow_type
  peptide_level <- route %in% c("DIA", "DIA_LIMPA")
  qc_specs <- c(
    if (peptide_level) lapply(getProtQcPeptideModuleSpecs(route), function(spec) c(spec, list(group = "peptide_qc", tabs = "peptide_filter_tabs"))),
    lapply(getProtQcProteinModuleSpecs(route), function(spec) c(spec, list(group = "protein_qc", tabs = "protein_filter_tabs"))))
  qc_tabset <- if (peptide_level) "quality_control-qc_tabs_lfq" else "quality_control-qc_tabs_tmt"
  settings <- guiEntryProtQcSettings()
  qc_step <- function(spec) {
    force(spec)
    guiEntryStep(paste0("qc_", sub("_qc$", "", spec$group), "_", spec$id), "qc", function(api) {
      sub_id <- paste0("quality_control-", spec$group, "-", spec$id)
      ui_fn <- get(spec$uiFn, envir = asNamespace("MultiScholaR"))
      api$set(stats::setNames(list(if (identical(spec$group, "peptide_qc")) "Peptide QC" else "Protein QC"), qc_tabset))
      api$set(stats::setNames(list(spec$title), paste0("quality_control-", spec$group, "-", spec$tabs)))
      api$defaults(ui_fn, sub_id)
      rendered <- names(omicsUiInputDefaults(ui_fn(sub_id)))
      overrides <- settings[paste0(sub_id, "-", names(settings)) %in% rendered]
      if (length(overrides)) api$set(stats::setNames(overrides, paste0(sub_id, "-", names(overrides))))
      actions <- grep("-apply_", guiEntryActionIds(ui_fn(sub_id)), value = TRUE)
      if (!length(actions)) stop("no apply action in ", spec$uiFn)
      api$set(stats::setNames(list(1L), actions[[1L]]))
    })
  }
  import_inputs <- function(api) {
    api$defaults(mod_prot_import_ui, "setup_import")
    api$defaults(function(id) buildProtImportFormatSpecificOptionsUi(lane$import_tool, shiny::NS(id)), "setup_import")
    api$upload("setup_import-search_results_standard", file.path(api$inputs$dir, lane$seed_file))
    api$upload("setup_import-fasta_file_standard", file.path(api$inputs$dir, lane$fasta_file))
    api$set(list(`setup_import-format_override` = lane$import_tool))
    if (peptide_level) api$set(list(`setup_import-dia_route` = route))
  }
  norm_steps <- list(
    guiEntryStep("norm_refusal", "normalization", refusal = "message", function(api) {
      # The filtered-session export before the normalisation.
      api$defaults(mod_prot_norm_ui, "normalization")
      api$set(list(`normalization-export_filtered_session` = 1L))
    }),
    guiEntryStep("norm", "normalization", function(api) {
      api$defaults(mod_prot_norm_ui, "normalization")
      if (!identical(route, "DIA_LIMPA")) api$set(list(`normalization-ruv_mode` = "skip"))
      api$set(list(`normalization-run_normalization` = 1L))
    }),
    guiEntryStep("correlation", "normalization", function(api) {
      api$set(list(`normalization-norm_qc_tabs` = "Correlation Filtering", `normalization-min_pearson_correlation_threshold` = 0.5))
      api$set(list(`normalization-apply_correlation_filter` = 1L))
    }),
    if (identical(route, "DIA_LIMPA")) guiEntryStep("ruv_covariates", "normalization", function(api) {
      api$set(list(`normalization-apply_ruv_covariates` = 1L))
    }),
    guiEntryStep("session_export", "normalization", function(api) api$set(list(`normalization-export_filtered_session` = 2L))))
  c(list(
    guiEntryStep("import_refusal", "setup", refusal = "silent", function(api) {
      # The process button with the format chosen and no search results file: gated (req), no message.
      api$defaults(mod_prot_import_ui, "setup_import")
      api$set(list(`setup_import-format_override` = lane$import_tool))
      api$set(list(`setup_import-process_data` = 1L))
    }),
    guiEntryStep("import", "setup", function(api) {
      import_inputs(api)
      api$set(list(`setup_import-process_data` = 2L))
    }),
    guiEntryStep("design_refusal", "design", refusal = "message", function(api) {
      api$set(list(`design_matrix-builder-save_results` = 1L))
    }),
    guiEntryStep("design", "design", function(api) {
      design <- guiEntryReadDesign(api)
      design$sample <- e2e_match_available_values(design$sample, api$choices("design_matrix-builder-selected_runs"))
      guiEntryDesignActions(api, design, c("KO", "WT"), prefix = "design_matrix-builder-", tech_rep = FALSE)
    })),
    lapply(qc_specs, qc_step),
    list(guiEntryStep("da_refusal", "da", refusal = "message", function(api) {
      # The filtered session loaded before any was exported.
      api$defaults(mod_prot_da_ui, "differential_expression")
      api$set(list(`differential_expression-load_filtered_session` = 1L))
    })),
    Filter(Negate(is.null), norm_steps),
    list(
    guiEntryStep("da_load", "da", function(api) api$set(list(`differential_expression-load_filtered_session` = 2L))),
    guiEntryStep("da", "da", function(api) {
      api$set(list(`differential_expression-formula_string` = "~ 0 + group", `differential_expression-da_q_val_thresh` = 0.05,
        `differential_expression-treat_lfc_cutoff` = 0))
      api$set(list(`differential_expression-run_da_analysis` = 1L))
      if (length(api$modals())) api$set(list(`differential_expression-acknowledge_qvalue_warning` = 1L))
    })),
    if (identical(lane$lane_id, "prot_dia")) list(
      guiEntryStep("enrichment", "enrichment", function(api) {
        api$defaults(mod_prot_enrich_ui, "enrichment_analysis")
        api$set(list(`enrichment_analysis-organism_taxid` = "9606", `enrichment_analysis-enrichment_method_tabs` = "gprofiler2"))
        api$set(list(`enrichment_analysis-selected_contrast` = lane$expected_contrasts[[1L]], `enrichment_analysis-up_cutoff` = 0,
          `enrichment_analysis-down_cutoff` = 0, `enrichment_analysis-q_cutoff` = 0.05))
        api$set(list(`enrichment_analysis-run_enrichment_analysis` = 1L))
      })),
    list(
    guiEntryStep("summary_refusal", "session_summary", refusal = "silent", function(api) {
      api$defaults(mod_prot_summary_ui, "session_summary")
      api$set(list(`session_summary-experiment_label` = lane$lane_id))
      api$set(list(`session_summary-generate_report` = 1L))
    }),
    guiEntryStep("summary_save", "session_summary", function(api) {
      api$set(list(`session_summary-experiment_label` = lane$lane_id, `session_summary-description` = "KIT-101 entry-point lane"))
      api$set(list(`session_summary-save_workflow_args` = 1L))
    }),
    guiEntryStep("summary_copy", "session_summary", function(api) api$set(list(`session_summary-copy_to_publication` = 1L))),
    guiEntryStep("summary_report", "session_summary", function(api) api$set(list(`session_summary-generate_report` = 2L))),
    guiEntryStep("summary_export", "session_summary", function(api) api$set(list(`session_summary-export_session_state` = 1L)))))
}
