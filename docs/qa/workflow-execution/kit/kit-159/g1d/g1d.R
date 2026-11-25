# KIT-159 G1D (ticket section 2.1): the prot_dia_limpa DA table at the oracle's tree against run 44953187101's base.
#
# Usage (run-tests scripts job, from the script commit's root):
#   Rscript docs/qa/workflow-execution/kit/kit-159/g1d/g1d.R <base commit> <head commit> <out dir>
# Internal (one child process per side, started by this script with its working directory at that side's tree):
#   Rscript --vanilla g1d.R child <side> <private dir>
#
# Per side (base = the oracle's tree, head = run 44953187101's base), one after the other: a detached worktree of the
# commit under $RUNNER_TEMP/g1d/<side>, the workflow's native step for a package tree, then a child process that loads
# the tree as the test jobs do (pkgload::load_all, the test helpers sourced as testthat::test_file sources them) and runs
# the lane's existing entry point guiEntryRunLane("prot_dia_limpa") with the lane's own fault hook at
# summary_refusal. The hook copies the DA table to $RUNNER_TEMP/g1d-private/<side>/ (outside <out dir>) with its digest
# computed as the lane computes it, then runs that step's own actions, so every step runs as the lane defines it. The
# lane's record is written to <out dir>/records/<side>/ (removed by the request's out_pattern before upload).
#
# Only three aggregate tables are written for upload: <out dir>/status.tsv (mechanism checks), g1d-captures.tsv
# (reproduction checks and per-step capture counts) and g1d-compare.tsv (the DA tables compared: what the lane's mask
# leaves visible, i.e. headers, row counts, key order, non-numeric values, missing values and the masked lines that
# differ (KIT-159-C13 item 1); per numeric column the largest absolute and relative difference; keys on one side only;
# calls changed). No feature-level row. A base side that cannot load or run the lane on this library is recorded
# `incompatible` (step and error class only) and the comparison is not made (KIT-159-C13 item 2).
# Exit status 0 unless a check is `mechanism` or `stop`.

LANE <- "prot_dia_limpa"
DA_KEY <- "results/proteomics/da_proteins/da_proteins_groupKO-groupWT_long_annot.tsv"
EXPECTED_DIGEST <- c(
  base = "80cfedaae858a2e93b87185583156c17b106b157ac02e3eacfc068edf7679b47",
  head = "8a4c2dd50ad521daf9b407f9d3beb2ba5bff2b3bd393b3a48fdf13e7bcd112b9")
RULE_COLUMNS <- c("log2FC", "raw_pvalue", "fdr_qvalue")
RULE_TOLERANCE <- 1e-9
QVALUE <- 0.05
CHILD_TIMEOUT_SECONDS <- 3600
# The upload scanner's patterns (.github/run-tests/scan.R, strict): a detail that matches one is not written.
UPLOAD_PATTERNS <- c(
  "\\.rs:[0-9]", "src/rust|\\.cargo|Cargo\\.(toml|lock)", "rustc|cargo |panicked at",
  "Compiling [a-z_]+ v[0-9]|error\\[E[0-9]+\\]|warning: unused",
  "(^|[ /])(gcc|g\\+\\+|clang|cc1|ld)( |:)|-fPIC|-std=gnu", "make(\\[[0-9]+\\])?: \\*\\*\\*",
  "\\.(c|cc|cpp|h|hpp|o):[0-9]+", "fast_limpa_dpc|fast_duckdb|internal-build-logs",
  "gh[pousr]_[A-Za-z0-9]{20,}|github_pat_|Bearer ")

oneLine <- function(text) gsub("[\t\r\n]+", " ", paste(as.character(text), collapse = " "))

cleanDetail <- function(text) {
  text <- substr(oneLine(text), 1L, 300L)
  hit <- vapply(UPLOAD_PATTERNS, function(pattern) grepl(pattern, text, useBytes = TRUE), logical(1L))
  if (any(hit)) "<detail not written: it matches an upload pattern>" else text
}

formatSeconds <- function(seconds) if (length(seconds) != 1L || is.na(seconds)) "" else sprintf("%.1f", seconds)

formatValue <- function(value) {
  if (is.null(value) || !length(value)) return("")
  value <- value[[1L]]
  if (is.na(value)) return("NA")
  if (is.double(value)) return(sprintf("%.17g", value))
  as.character(value)
}

# Evaluate `expr`; return list(ok, value or error, seconds).
attempt <- function(expr) {
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(list(ok = TRUE, value = expr), error = function(e) list(ok = FALSE, error = e))
  result$seconds <- proc.time()[["elapsed"]] - started
  result
}

errorDetail <- function(error) paste0(class(error)[[1L]], ": ", conditionMessage(error))

# --------------------------------------------------------------------------------------------------------------------
# Child: one side, working directory at that side's tree root.
# --------------------------------------------------------------------------------------------------------------------

childMain <- function(side, private_dir) {
  stages_path <- file.path(private_dir, "child-stages.tsv")
  writeStage <- function(stage, outcome, detail, seconds) {
    cat(paste(stage, outcome, oneLine(detail), formatSeconds(seconds), sep = "\t"), "\n", sep = "", file = stages_path,
      append = TRUE)
  }
  # KIT-159-C13 item 2: on the base side (the 28d5a4be tree on the 76b57c0a library) an error of the load or of a lane
  # step is a result, `incompatible`, written with the step and the error's class only (no message text).
  failureOutcome <- function(stage, error, seconds, step = NA_character_) {
    if (identical(side, "base")) {
      writeStage(stage, "incompatible", sprintf("step=%s; class=%s", if (is.na(step)) stage else step,
        class(error)[[1L]]), seconds)
    } else {
      writeStage(stage, "mechanism", if (is.na(step)) errorDetail(error) else paste0("step=", step, "; ",
        errorDetail(error)), seconds)
    }
  }
  # As test_that() runs a test (testthat edition 3): width 80, no colour or unicode, C collation and language.
  testthat::local_reproducible_output()

  loaded <- attempt({
    suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
    suppressPackageStartupMessages(library(testthat))
    # testthat::test_file() sources the helpers again into the test environment (a child of the global environment),
    # with TESTTHAT set and the working directory at the test directory; the test files call those copies.
    helper_env <- new.env(parent = globalenv())
    test_dir <- normalizePath("tests/testthat")
    withr::with_dir(test_dir, withr::with_envvar(c(TESTTHAT = "true"),
      suppressMessages(testthat::source_test_helpers(test_dir, env = helper_env))))
    helper_env
  })
  if (!loaded$ok) {
    failureOutcome("load", loaded$error, loaded$seconds)
    return(1L)
  }
  writeStage("load", "ok", "", loaded$seconds)
  helper_env <- loaded$value

  run_lane <- get("guiEntryRunLane", envir = helper_env, mode = "function")
  mask <- get("guiEntryMask", envir = helper_env, mode = "function")
  copy_path <- file.path(private_dir, "da-table.tsv")
  digest_path <- file.path(private_dir, "da-table-digest.txt")

  # The lane's fault hook at summary_refusal (helper-gui-entry.R guiEntryRunLane: `fault$fn()` in place of the step's
  # actions). The DA table is copied as step `da` left it, its digest computed as guiEntryArtifacts computes it, and
  # the step's own actions are then run, so the later steps run as the lane defines them.
  fault_fn <- function() {
    caller <- parent.frame()
    step <- get("step", envir = caller, mode = "list")
    api <- get("api", envir = caller, mode = "list")
    if (!identical(step$id, "summary_refusal")) stop("the fault hook ran at an unexpected step")
    roots <- NULL
    for (frame in rev(sys.frames())) {
      if (exists("roots", envir = frame, inherits = FALSE) && exists("acc", envir = frame, inherits = FALSE)) {
        roots <- get("roots", envir = frame, inherits = FALSE)
        break
      }
    }
    if (is.null(roots)) stop("the lane's masking roots were not found")
    base_dir <- get("project_dirs", envir = globalenv())[["proteomics"]][["base_dir"]]
    source_path <- file.path(base_dir, DA_KEY)
    if (!file.exists(source_path)) stop("the DA table is absent after step da")
    if (!file.copy(source_path, copy_path, overwrite = TRUE)) stop("the DA table could not be copied")
    text <- readLines(source_path, warn = FALSE)
    masked <- mask(text, roots)
    writeLines(digest::digest(masked, algo = "sha256", serialize = FALSE), digest_path)
    # The masked text itself (what the digest sees), private, for the line comparison in the parent.
    writeLines(masked, file.path(private_dir, "da-table-masked.txt"))
    step$run(api)
  }

  # KIT-159-C13 item 2: on the base side (the 28d5a4be tree on the 76b57c0a library) an error of the load or of a lane
  # step is a result, `incompatible`, written with the step and the error's class only (no message text).
  ids <- tryCatch(names(jsonlite::read_json(file.path("tests", "testdata", "gui-entry", paste0(LANE, ".json")),
    simplifyVector = FALSE)$steps), error = function(e) character())
  failingStep <- function() {
    completed <- character()
    for (frame in rev(sys.frames())) {
      if (exists("roots", envir = frame, inherits = FALSE) && exists("acc", envir = frame, inherits = FALSE)) {
        completed <- names(get("acc", envir = frame, inherits = FALSE)$steps) %||% character()
        break
      }
    }
    following <- setdiff(ids, completed)
    if (length(following)) following[[1L]] else paste0("after ", if (length(completed)) completed[[length(completed)]] else "start")
  }
  failed_step <- NA_character_

  # An expectation inside a lane step (summary_report) is recorded and the lane continues, as inside test_that().
  failures <- 0L
  lane <- attempt(withCallingHandlers(
    run_lane(LANE, fault = list(step = "summary_refusal", fn = fault_fn)),
    expectation = function(condition) {
      if (inherits(condition, c("expectation_failure", "expectation_error"))) failures <<- failures + 1L
      for (name in c("muffle_expectation", "continue_test")) {
        restart <- findRestart(name)
        if (!is.null(restart)) invokeRestart(restart)
      }
    },
    error = function(condition) {
      # The last error signalled names the step the lane stopped at (errors a step catches itself come earlier).
      if (!inherits(condition, c("expectation", "shiny.silent.error"))) failed_step <<- tryCatch(failingStep(),
        error = function(e) NA_character_)
    }))
  if (!lane$ok) {
    failureOutcome("lane", lane$error, lane$seconds, failed_step)
    return(1L)
  }
  writeStage("lane", "ok", sprintf("lane record %s s; %d failed expectations inside the lane's steps",
    formatSeconds(as.numeric(lane$value$seconds)), failures), lane$seconds)
  0L
}

# --------------------------------------------------------------------------------------------------------------------
# Parent.
# --------------------------------------------------------------------------------------------------------------------

parentMain <- function(base_commit, head_commit, out_arg, script_path) {
  root <- normalizePath(".")
  out <- file.path(root, out_arg)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out)
  status_path <- file.path(out, "status.tsv")
  compare_path <- file.path(out, "g1d-compare.tsv")
  captures_path <- file.path(out, "g1d-captures.tsv")
  writeLines(paste("stage", "side", "outcome", "detail", "seconds", sep = "\t"), status_path)
  writeLines(paste("quantity", "column", "value", sep = "\t"), compare_path)
  writeLines(paste("quantity", "step", "field", "value", sep = "\t"), captures_path)

  failed <- FALSE
  status <- function(stage, side, outcome, detail = "", seconds = NA_real_) {
    # `incompatible` (the base side cannot run on this library) is a result, not a failure of the script.
    if (outcome %in% c("mechanism", "stop")) failed <<- TRUE
    cat(paste(stage, side, outcome, cleanDetail(detail), formatSeconds(seconds), sep = "\t"), "\n", sep = "",
      file = status_path, append = TRUE)
    cat(sprintf("%s %s: %s\n", side, stage, outcome))
  }
  # A line that would match an upload pattern (a column name, say) is written with its text replaced, so that the
  # table is not withheld.
  safeLine <- function(fields) {
    line <- paste(vapply(fields, oneLine, ""), collapse = "\t")
    hit <- vapply(UPLOAD_PATTERNS, function(pattern) grepl(pattern, line, useBytes = TRUE), logical(1L))
    if (any(hit)) paste(c(fields[[1L]], rep("<not written: matches an upload pattern>", length(fields) - 1L)),
      collapse = "\t") else line
  }
  emitCompare <- function(quantity, column, value) {
    cat(safeLine(c(quantity, column, formatValue(value))), "\n", sep = "", file = compare_path, append = TRUE)
  }
  emitCapture <- function(quantity, step, field, value) {
    cat(safeLine(c(quantity, step, field, formatValue(value))), "\n", sep = "", file = captures_path, append = TRUE)
  }

  commits <- c(base = base_commit, head = head_commit)
  if (!all(grepl("^[0-9a-f]{40}$", commits))) {
    status("arguments", "both", "mechanism", "the two commits must be full 40-digit hashes")
    return(1L)
  }
  temp <- Sys.getenv("RUNNER_TEMP")
  if (!nzchar(temp)) temp <- tempdir()
  trees_dir <- file.path(temp, "g1d")
  private_root <- file.path(temp, "g1d-private")
  dir.create(trees_dir, recursive = TRUE, showWarnings = FALSE)
  rscript <- file.path(R.home("bin"), "Rscript")

  sides <- list()
  for (side in names(commits)) {
    sides[[side]] <- runSide(side, commits[[side]], root, trees_dir, file.path(private_root, side),
      file.path(out, "records", side), rscript, script_path, status)
  }

  # Reproduction checks and per-step capture counts (from the two records).
  records <- lapply(sides, `[[`, "record")
  for (side in names(commits)) {
    observed <- if (is.null(records[[side]])) "<absent>" else records[[side]]$steps$da$artifacts[[DA_KEY]] %||% "<absent>"
    emitCapture("da_digest", "da", side, observed)
    emitCapture("da_digest_expected", "da", side, EXPECTED_DIGEST[[side]])
    emitCapture("reproduced", "da", side, identical(observed, EXPECTED_DIGEST[[side]]))
  }
  oracle <- sides$base$oracle
  emitCapture("oracle_file_da_digest", "da", "base",
    if (is.null(oracle)) "<absent>" else oracle$steps$da$artifacts[[DA_KEY]] %||% "<absent>")
  if (!is.null(records$base) && !is.null(records$head)) {
    steps <- unique(c(names(records$base$steps), names(records$head$steps)))
    for (step in steps) {
      for (field in c("states", "outputs", "artifacts")) {
        a <- records$base$steps[[step]][[field]]
        b <- records$head$steps[[step]][[field]]
        keys_a <- names(a) %||% character()
        keys_b <- names(b) %||% character()
        both <- intersect(keys_a, keys_b)
        differ <- sum(!vapply(both, function(key) identical(a[[key]], b[[key]]), logical(1L)))
        emitCapture("keys_base", step, field, length(keys_a))
        emitCapture("keys_head", step, field, length(keys_b))
        emitCapture("keys_both", step, field, length(both))
        emitCapture("keys_differ", step, field, differ)
      }
    }
  }

  # The comparison of the two DA tables (from the private copies).
  if (isTRUE(sides$base$incompatible) || isTRUE(sides$head$incompatible)) {
    status("compare", "both", "not_made", "a side is incompatible with this library (KIT-159-C13 item 2)")
  } else if (isTRUE(sides$base$copy_ok) && isTRUE(sides$head$copy_ok)) {
    compared <- attempt(compareTables(sides$base$copy_path, sides$head$copy_path, status, emitCompare,
      file.path(private_root, "base", "da-table-masked.txt"), file.path(private_root, "head", "da-table-masked.txt")))
    if (!compared$ok) status("compare", "both", "mechanism", errorDetail(compared$error), compared$seconds)
    else if (isTRUE(compared$value)) status("compare", "both", "ok", "", compared$seconds)
  } else {
    status("compare", "both", "stop", "a DA table copy is missing or does not match its record")
  }
  if (failed) 1L else 0L
}

runSide <- function(side, commit, root, trees_dir, private_dir, record_dir, rscript, script_path, status) {
  result <- list(record = NULL, oracle = NULL, copy_ok = FALSE, copy_path = file.path(private_dir, "da-table.tsv"))
  dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(record_dir, recursive = TRUE, showWarnings = FALSE)
  tree <- file.path(trees_dir, side)
  git_log <- file.path(private_dir, "git.log")

  made <- attempt(processx::run("git", c("-C", root, "worktree", "add", "--detach", tree, commit),
    error_on_status = FALSE, stdout = git_log, stderr = "2>&1"))
  if (!made$ok) { status("worktree", side, "mechanism", errorDetail(made$error), made$seconds); return(result) }
  if (made$value$status != 0L) {
    status("worktree", side, "mechanism", sprintf("git worktree add exited %s", made$value$status), made$seconds)
    return(result)
  }
  status("worktree", side, "ok", "", made$seconds)

  head <- attempt(processx::run("git", c("-C", tree, "rev-parse", "HEAD"), error_on_status = FALSE))
  observed <- if (head$ok) trimws(head$value$stdout) else ""
  if (!identical(observed, commit)) {
    status("head_commit", side, "mechanism", if (head$ok) "HEAD differs from the commit" else errorDetail(head$error),
      head$seconds)
    return(result)
  }
  status("head_commit", side, "ok", "", head$seconds)

  # The workflow's native-code step for a package tree (run-tests.yml, job scripts, step "Install package_commit").
  if (file.exists(file.path(tree, "src", "build_dia_rsum.R")) &&
      !file.exists(file.path(tree, "inst", "duckdb", "linux_amd64", "manifest.rds"))) {
    native <- attempt(processx::run(rscript, c("--vanilla", "src/build_dia_rsum.R"), wd = tree, error_on_status = FALSE,
      stdout = file.path(private_dir, "native.log"), stderr = "2>&1"))
    if (!native$ok) { status("native", side, "mechanism", errorDetail(native$error), native$seconds); return(result) }
    if (native$value$status != 0L || !file.exists(file.path(tree, "inst", "duckdb", "linux_amd64", "manifest.rds"))) {
      status("native", side, "mechanism", sprintf("src/build_dia_rsum.R exited %s", native$value$status), native$seconds)
      return(result)
    }
    status("native", side, "ok", "", native$seconds)
  } else {
    status("native", side, "ok", "not needed (no build script, or the manifest is present)")
  }

  # The child: its own process, working directory at the tree's root, the lane's record directory set.
  child <- attempt(processx::run(rscript, c("--vanilla", script_path, "child", side, private_dir), wd = tree,
    env = c("current", MULTISCHOLAR_GUI_ENTRY_RECORD_DIR = record_dir), error_on_status = FALSE,
    stdout = file.path(private_dir, "child.log"), stderr = "2>&1", timeout = CHILD_TIMEOUT_SECONDS, cleanup_tree = TRUE))
  stages_path <- file.path(private_dir, "child-stages.tsv")
  incompatible <- FALSE
  if (file.exists(stages_path)) {
    for (line in readLines(stages_path, warn = FALSE)) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      fields <- c(fields, rep("", max(0L, 4L - length(fields))))[1:4]
      if (identical(fields[[2L]], "incompatible")) incompatible <- TRUE
      status(fields[[1L]], side, fields[[2L]], fields[[3L]], suppressWarnings(as.numeric(fields[[4L]])))
    }
  }
  if (!child$ok) { status("child", side, "mechanism", errorDetail(child$error), child$seconds); return(result) }
  if (incompatible && !isTRUE(child$value$timeout)) {
    status("child", side, "incompatible", sprintf("the child exited %s", child$value$status), child$seconds)
    result$incompatible <- TRUE
    return(result)
  }
  if (isTRUE(child$value$timeout)) {
    status("child", side, "mechanism", sprintf("timeout after %d s", CHILD_TIMEOUT_SECONDS), child$seconds)
    return(result)
  }
  if (!identical(as.integer(child$value$status), 0L)) {
    status("child", side, "mechanism", sprintf("the child exited %s", child$value$status), child$seconds)
    return(result)
  }
  status("child", side, "ok", "", child$seconds)

  # The record: present, with the oracle's step ids in order.
  oracle_path <- file.path(tree, "tests", "testdata", "gui-entry", paste0(LANE, ".json"))
  record_path <- file.path(record_dir, paste0(LANE, ".json"))
  read <- attempt(list(oracle = jsonlite::read_json(oracle_path, simplifyVector = FALSE),
    record = jsonlite::read_json(record_path, simplifyVector = FALSE)))
  if (!read$ok) { status("record", side, "mechanism", errorDetail(read$error), read$seconds); return(result) }
  result$oracle <- read$value$oracle
  result$record <- read$value$record
  oracle_ids <- names(read$value$oracle$steps)
  record_ids <- names(read$value$record$steps)
  if (length(oracle_ids) != 24L || !identical(record_ids, oracle_ids)) {
    status("record", side, "stop", sprintf("the record has %d step ids, the oracle %d; equal in order: %s",
      length(record_ids), length(oracle_ids), identical(record_ids, oracle_ids)))
    return(result)
  }
  status("record", side, "ok", sprintf("%d step ids, as the oracle's", length(record_ids)))

  # The copy: present, its digest equal to the record's digest of the DA table at step da.
  digest_path <- file.path(private_dir, "da-table-digest.txt")
  recorded <- read$value$record$steps$da$artifacts[[DA_KEY]] %||% "<absent>"
  copied <- if (file.exists(digest_path)) readLines(digest_path, warn = FALSE)[[1L]] else "<absent>"
  if (!file.exists(result$copy_path) || !identical(copied, recorded)) {
    status("copy", side, "stop", sprintf("copy present: %s; its digest equals the record's at da: %s",
      file.exists(result$copy_path), identical(copied, recorded)))
    return(result)
  }
  status("copy", side, "ok", "the copy's digest equals the record's at da")
  result$copy_ok <- TRUE
  result
}

# The two DA tables by key; returns TRUE when compared, FALSE when a check stopped it.
compareTables <- function(base_path, head_path, status, emit, masked_base_path, masked_head_path) {
  read <- function(path) data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE,
    integer64 = "double")
  tables <- list(base = read(base_path), head = read(head_path))
  # KIT-159-C13 item 1: what the oracle's mask leaves visible. Each side's header (column names in order; design
  # labels, no feature row) and its row count.
  headers_identical <- identical(names(tables$base), names(tables$head))
  for (side in names(tables)) {
    emit(paste0("header_", side), "", paste(names(tables[[side]]), collapse = " | "))
    emit(paste0("columns_", side), "", ncol(tables[[side]]))
  }
  emit("header_identical", "", headers_identical)
  # Missing values per column per side (every column of either side).
  na_counts_equal <- TRUE
  for (column in union(names(tables$base), names(tables$head))) {
    counts <- vapply(tables, function(table) if (column %in% names(table)) sum(is.na(table[[column]])) else NA_integer_,
      integer(1L))
    emit("na_count_base", column, counts[["base"]])
    emit("na_count_head", column, counts[["head"]])
    if (!identical(counts[["base"]], counts[["head"]])) na_counts_equal <- FALSE
  }
  emit("na_counts_equal", "", na_counts_equal)
  required <- c("Protein.Ids", RULE_COLUMNS)
  ok <- TRUE
  for (side in names(tables)) {
    missing <- setdiff(required, names(tables[[side]]))
    if (length(missing)) {
      status("columns", side, "stop", paste("missing:", paste(missing, collapse = ", ")))
      ok <- FALSE
    } else {
      status("columns", side, "ok", sprintf("%d columns", ncol(tables[[side]])))
    }
  }
  if (!ok) return(FALSE)
  use_comparison <- all(vapply(tables, function(table) "comparison" %in% names(table), logical(1L)))
  keyOf <- function(table) {
    if (use_comparison) paste(table$Protein.Ids, table$comparison, sep = "\t") else as.character(table$Protein.Ids)
  }
  keys <- lapply(tables, keyOf)
  for (side in names(tables)) {
    duplicated_keys <- sum(duplicated(keys[[side]]))
    if (duplicated_keys > 0L) {
      status("keys", side, "stop", sprintf("%d duplicated keys", duplicated_keys))
      ok <- FALSE
    } else {
      status("keys", side, "ok", sprintf("%d rows, keys unique", nrow(tables[[side]])))
    }
  }
  if (!ok) return(FALSE)

  common <- intersect(keys$base, keys$head)
  index_base <- match(common, keys$base)
  index_head <- match(common, keys$head)
  only_base <- sum(!keys$base %in% keys$head)
  only_head <- sum(!keys$head %in% keys$base)
  emit("key", "", if (use_comparison) "Protein.Ids+comparison" else "Protein.Ids")
  emit("rows_base", "", nrow(tables$base))
  emit("rows_head", "", nrow(tables$head))
  emit("keys_both", "", length(common))
  emit("keys_only_base", "", only_base)
  emit("keys_only_head", "", only_head)
  key_order_identical <- identical(keys$base, keys$head)
  emit("key_order_identical", "", key_order_identical)
  for (column in setdiff(names(tables$base), names(tables$head))) emit("column_only_base", column, 1L)
  for (column in setdiff(names(tables$head), names(tables$base))) emit("column_only_head", column, 1L)

  max_rel <- list()
  one_side <- list()
  numeric_columns <- character()
  non_numeric_differ <- 0L
  types_differ <- 0L
  for (column in intersect(names(tables$base), names(tables$head))) {
    a <- tables$base[[column]][index_base]
    b <- tables$head[[column]][index_head]
    if (is.logical(a) && all(is.na(a)) && is.numeric(b)) a <- as.numeric(a)
    if (is.logical(b) && all(is.na(b)) && is.numeric(a)) b <- as.numeric(b)
    if (is.numeric(a) && is.numeric(b)) {
      a <- as.numeric(a)
      b <- as.numeric(b)
      both <- !is.na(a) & !is.na(b)
      missing_one <- xor(is.na(a), is.na(b))
      x <- a[both]
      y <- b[both]
      same <- x == y
      absolute <- abs(x - y)
      absolute[same] <- 0
      absolute[is.nan(absolute)] <- Inf
      relative <- absolute / pmax(abs(x), abs(y))
      relative[same] <- 0
      relative[is.na(relative)] <- Inf
      max_rel[[column]] <- if (length(relative)) max(relative) else 0
      one_side[[column]] <- sum(missing_one)
      numeric_columns <- c(numeric_columns, column)
      emit("rows_compared", column, sum(both))
      emit("max_abs_diff", column, if (length(absolute)) max(absolute) else 0)
      emit("max_rel_diff", column, max_rel[[column]])
      emit("na_one_side", column, sum(missing_one))
      emit("values_differ", column, sum(!same) + sum(missing_one))
    } else if (is.numeric(a) != is.numeric(b)) {
      emit("column_type_differs", column, paste(class(a)[[1L]], class(b)[[1L]], sep = "/"))
      types_differ <- types_differ + 1L
    } else {
      x <- as.character(a)
      y <- as.character(b)
      differ <- xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
      emit("values_differ", column, sum(differ))
      non_numeric_differ <- non_numeric_differ + sum(differ)
    }
  }
  emit("non_numeric_values_differ", "", non_numeric_differ)
  emit("column_types_differ", "", types_differ)

  # The masked texts the two digests are taken of, line by line in file order: lines that differ, and how many of
  # those differ only inside fields of columns numeric on both sides (same field count, same column names there).
  masked <- list(base = readLines(masked_base_path, warn = FALSE), head = readLines(masked_head_path, warn = FALSE))
  splitFields <- function(line) {
    fields <- strsplit(paste0(line, "\tEND"), "\t", fixed = TRUE)[[1L]]
    fields[-length(fields)]
  }
  shared <- min(length(masked$base), length(masked$head))
  lines_differ <- abs(length(masked$base) - length(masked$head))
  numeric_only <- 0L
  for (i in seq_len(shared)) {
    if (identical(masked$base[[i]], masked$head[[i]])) next
    lines_differ <- lines_differ + 1L
    if (i == 1L) next
    x <- splitFields(masked$base[[i]])
    y <- splitFields(masked$head[[i]])
    if (length(x) != length(y)) next
    positions <- which(x != y)
    names_x <- names(tables$base)[positions]
    names_y <- names(tables$head)[positions]
    if (length(positions) && !anyNA(names_x) && identical(names_x, names_y) && all(names_x %in% numeric_columns)) {
      numeric_only <- numeric_only + 1L
    }
  }
  emit("masked_lines_base", "", length(masked$base))
  emit("masked_lines_head", "", length(masked$head))
  emit("masked_lines_differ", "", lines_differ)
  emit("masked_lines_differ_numeric_only", "", numeric_only)
  emit("masked_lines_differ_other", "", lines_differ - numeric_only)

  # Calls at the route's defaults: up when fdr_qvalue < 0.05 and log2FC > 0, down when fdr_qvalue < 0.05 and
  # log2FC < 0, otherwise not called.
  callOf <- function(table) {
    q <- suppressWarnings(as.numeric(table$fdr_qvalue))
    l <- suppressWarnings(as.numeric(table$log2FC))
    called <- !is.na(q) & q < QVALUE & !is.na(l)
    ifelse(called & l > 0, "up", ifelse(called & l < 0, "down", "none"))
  }
  calls <- lapply(tables, callOf)
  calls_changed <- sum(calls$base[index_base] != calls$head[index_head])
  emit("calls_changed", "", calls_changed)
  for (side in names(tables)) {
    emit(paste0("called_", side), "", sum(calls[[side]] != "none"))
    emit(paste0("called_up_", side), "", sum(calls[[side]] == "up"))
    emit(paste0("called_down_", side), "", sum(calls[[side]] == "down"))
  }

  # The decision rule's quantities (section 2.1); applied to the result only when both reproduction checks hold.
  rule_rel <- if (all(RULE_COLUMNS %in% names(max_rel))) max(unlist(max_rel[RULE_COLUMNS])) else NA_real_
  rule_one_side <- if (all(RULE_COLUMNS %in% names(one_side))) sum(unlist(one_side[RULE_COLUMNS])) else NA_integer_
  # KIT-159-C13: the structural items must be equal as well (headers, rows, key order, non-numeric values, column
  # types, missing values, masked lines differing only inside numeric fields).
  structure_equal <- headers_identical && nrow(tables$base) == nrow(tables$head) && key_order_identical &&
    non_numeric_differ == 0L && types_differ == 0L && na_counts_equal && lines_differ == numeric_only
  rule_met <- !is.na(rule_rel) && rule_rel <= RULE_TOLERANCE && !is.na(rule_one_side) && rule_one_side == 0L &&
    calls_changed == 0L && only_base == 0L && only_head == 0L && structure_equal
  emit("rule_max_rel_diff", paste(RULE_COLUMNS, collapse = "+"), rule_rel)
  emit("rule_na_one_side", paste(RULE_COLUMNS, collapse = "+"), rule_one_side)
  emit("rule_structure_equal", "", structure_equal)
  emit("rule_met", "", rule_met)
  TRUE
}

# --------------------------------------------------------------------------------------------------------------------
# Entry.
# --------------------------------------------------------------------------------------------------------------------

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) == 3L && identical(arguments[[1L]], "child")) {
  code <- tryCatch(childMain(arguments[[2L]], normalizePath(arguments[[3L]])), error = function(e) {
    cat(paste("child", "mechanism", oneLine(errorDetail(e)), "", sep = "\t"), "\n", sep = "",
      file = file.path(arguments[[3L]], "child-stages.tsv"), append = TRUE)
    1L
  })
  quit(save = "no", status = code)
}
if (length(arguments) != 3L) stop("usage: Rscript g1d.R <base commit> <head commit> <out dir>")
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]))
code <- tryCatch(parentMain(arguments[[1L]], arguments[[2L]], arguments[[3L]], script_path), error = function(e) {
  out <- file.path(normalizePath("."), arguments[[3L]])
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  cat(paste("script", "both", "mechanism", cleanDetail(errorDetail(e)), "", sep = "\t"), "\n", sep = "",
    file = file.path(out, "status.tsv"), append = TRUE)
  cat("script: mechanism\n")
  1L
})
quit(save = "no", status = code)
