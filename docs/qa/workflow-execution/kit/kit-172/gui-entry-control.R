# KIT-172 section 7(c): the GUI-entry control on live runs and the re-freeze of the lipid and metabolomics oracles.
#
# Usage (run-tests scripts lane, from the package root of the head; one request; `out` = control-out):
#   Rscript docs/qa/workflow-execution/kit/kit-172/gui-entry-control.R <base sha> <out dir> [--base-records <dir>]
#
# Three passes, each a child Rscript process under pkgload::load_all (the oracles' load mode) with the tests' helpers
# in an environment under the global one: `base` in a git worktree of the base commit; `head1`, the head with the
# normalisation tab's defaults carrying `norm-log_offset = 1` (a script-local wrapper of `omicsUiInputDefaults`; no
# helper file is edited); `head0`, the head at its defaults under MULTISCHOLAR_GUI_ENTRY_FREEZE=1 (guiEntryExpectLane:
# the lane's record and its re-frozen oracle). Each child's working directory is its own tree's root, and its records
# go to <out>/<pass>/ (MULTISCHOLAR_GUI_ENTRY_RECORD_DIR, absolute). The passes run side by side (one process each,
# the runner's 4 cores); `lipid_canonical` first in each pass, then the mechanism is checked (records written, the base
# child's source SHA, the `norm_log` capture at `norm` of head1 carrying "(offset: 1 )" and of head0 not), then the
# three metabolomics lanes. With --base-records the base records are read from <dir> and the base pass is not run.
#
# Out: status.tsv (step, lane, outcome: pass, stop with its reason, mechanism with the error), gui-entry-control.tsv
# (compareLane's rows, packaging/verify/run-lanes.R:104-135, with the comparison, the class, the phase and whether the
# step is before `norm`; first line `# cpu: <model>`), gui-entry-digests.tsv (per lane and step the digest of the whole
# capture: frozen, base, head1, head0), declared-keys.tsv (the declared keys' expected and observed values at every
# step), k156-gui-entry-failures.tsv (the reference rows used, copied), the three passes' records and refrozen/<lane>.json.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L, grepl("^[0-9a-f]{40}$", args[[1L]]))
base_sha <- args[[1L]]
dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[[2L]])
base_records <- if (length(args) >= 4L && identical(args[[3L]], "--base-records")) normalizePath(args[[4L]]) else NULL
head_root <- normalizePath(".")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(out, "status.tsv")
cat("step\tlane\toutcome\n", file = status_file)
status <- function(step, lane, outcome) cat(sprintf("%s\t%s\t%s\n", step, lane, gsub("[\t\n]", " ", outcome)), file = status_file, append = TRUE)
finish <- function() quit(save = "no", status = 0L)

lanes <- c("lipid_canonical", "metab_lc", "metab_gc", "metab_combined")
omic_of <- function(lane) if (startsWith(lane, "lipid")) "lipidomics" else "metabolomics"
cpu <- tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]), error = function(e) "unknown")
status("start", "", sprintf("head %s, base %s, cpu %s, base records %s", paste(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), collapse = ""),
  base_sha, cpu, base_records %||% "none"))

# --- compareLane, extracted as a definition (the file is never sourced) ---------------------------------------------
lanes_file <- parse(file.path(head_root, "packaging", "verify", "run-lanes.R"), keep.source = FALSE)
for (expression in lanes_file) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) && is.name(expression[[2L]]) &&
    as.character(expression[[2L]]) %in% c("compareLane", "truncate")) {
    eval(expression, globalenv())
  }
}
stopifnot(is.function(compareLane), is.function(truncate))

# --- The reference rows (step 1) ---------------------------------------------------------------------------------
reference <- utils::read.delim(file.path(head_root, "docs", "qa", "workflow-execution", "kit", "kit-172", "k156-gui-entry-failures.tsv"),
  stringsAsFactors = FALSE)
file.copy(file.path(head_root, "docs", "qa", "workflow-execution", "kit", "kit-172", "k156-gui-entry-failures.tsv"), out, overwrite = TRUE)

# --- The lane fixtures hold no missing or zero cell in a sample column (step 2's not-tested count is 0) --------------
manifest <- jsonlite::read_json(file.path(head_root, "tests", "testdata", "e2e", "manifest.json"), simplifyVector = FALSE)$lanes
for (lane in lanes) {
  entry <- Filter(function(x) identical(x$lane_id, lane), manifest)[[1L]]
  for (file in unlist(entry[c("seed_file", "assay2_file")])) {
    table <- utils::read.delim(file.path(head_root, "tests", "testdata", "e2e", entry$fixture_dir, file), colClasses = "character",
      check.names = FALSE, na.strings = character())
    samples <- grep("^(WT|KO)_", names(table), value = TRUE)
    cells <- unlist(table[samples], use.names = FALSE)
    bad <- sum(!nzchar(trimws(cells)) | cells %in% c("NA", "NaN") | suppressWarnings(as.numeric(cells)) %in% 0)
    if (!length(samples) || bad > 0L) {
      status("fixtures", lane, sprintf("stop: %s has %d sample columns and %d missing or zero sample cells", file, length(samples), bad))
      finish()
    }
  }
  status("fixtures", lane, "pass")
}

# --- The frozen oracles, copied before head0 re-freezes them -------------------------------------------------------
frozen_dir <- file.path(out, "frozen-before")
dir.create(frozen_dir, showWarnings = FALSE)
for (lane in lanes) file.copy(file.path(head_root, "tests", "testdata", "gui-entry", paste0(lane, ".json")), frozen_dir, overwrite = TRUE)

# --- The base worktree --------------------------------------------------------------------------------------------
scratch <- Sys.getenv("RUNNER_TEMP", tempdir())
worktree <- file.path(scratch, "kit172-base")
if (is.null(base_records) && system2("git", c("worktree", "add", "--detach", shQuote(worktree), base_sha)) != 0L) {
  status("worktree", "", "mechanism: git worktree add failed")
  finish()
}

# --- The child ----------------------------------------------------------------------------------------------------
child <- file.path(scratch, "kit172-control-child.R")
writeLines(c(
  'args <- commandArgs(trailingOnly = TRUE); pass <- args[[1L]]; record_dir <- args[[2L]]; lanes <- strsplit(args[[3L]], ",")[[1L]]',
  '# VAL-55: sort as testthat does (testthat 3.3.2 local_reproducible_output: withr::local_collate("C"), i.e. LC_COLLATE C).',
  'invisible(Sys.setlocale("LC_COLLATE", "C"))',
  'Sys.setenv(MULTISCHOLAR_GUI_ENTRY_RECORD_DIR = record_dir, NOT_CRAN = "true", MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR = "")',
  'if (identical(pass, "head0")) Sys.setenv(MULTISCHOLAR_GUI_ENTRY_FREEZE = "1") else Sys.unsetenv("MULTISCHOLAR_GUI_ENTRY_FREEZE")',
  'suppressMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, attach_testthat = FALSE, quiet = TRUE))',
  'env <- new.env(parent = globalenv())',
  'suppressMessages(testthat::source_test_helpers("tests/testthat", env = env))',
  'dir.create(record_dir, recursive = TRUE, showWarnings = FALSE)',
  'writeLines(env$guiEntrySourceSha(), file.path(record_dir, "source-sha.txt"))',
  'if (identical(pass, "head1")) {',
  '  packaged <- env$omicsUiInputDefaults',
  '  env$omicsUiInputDefaults <- function(...) {',
  '    values <- packaged(...)',
  '    at <- which(endsWith(as.character(names(values)), "norm-norm_method"))',
  '    if (length(at)) values[[sub("norm_method$", "log_offset", names(values)[[at[[1L]]]])]] <- 1',
  '    values',
  '  }',
  '}',
  'if (!identical(pass, "base")) {',
  '  texts <- get("omicGuiLog2Texts", envir = asNamespace("MultiScholaR"))()',
  '  # The capture of the sentence as a renderUI output, through the lanes own output reader.',
  '  rendered <- NULL',
  '  shiny::testServer(function(id) shiny::moduleServer(id, function(input, output, session) {',
  '    output$x <- shiny::renderUI(shiny::helpText(sprintf(texts[["recorded_offset"]], "1")))',
  '  }), { rendered <<- env$guiEntryOutputValue(output$x, character()) })',
  '  writeLines(rendered, file.path(record_dir, "expected-recorded-offset.txt"))',
  '}',
  'for (lane in lanes) {',
  '  started <- Sys.time()',
  '  if (identical(pass, "head0")) env$guiEntryExpectLane(lane) else env$guiEntryRunLane(lane)',
  '  cat(sprintf("%s\\t%.1f\\n", lane, as.numeric(difftime(Sys.time(), started, units = "mins"))), file = file.path(record_dir, "minutes.tsv"), append = TRUE)',
  '}'), child)

# Native code compiled once per tree before the passes start side by side (two load_all in one tree at once would race).
compile <- file.path(scratch, "kit172-compile.R")
writeLines('suppressMessages(pkgload::load_all(".", quiet = TRUE))', compile)
for (tree in c(if (is.null(base_records)) worktree, head_root)) {
  code <- system2("bash", c("-c", shQuote(sprintf("cd %s && Rscript %s > %s 2>&1", shQuote(tree), shQuote(compile),
    shQuote(file.path(out, paste0("compile-", basename(tree), ".log")))))))
  status("compile", basename(tree), if (code == 0L) "pass" else sprintf("exit %d", code))
}

runPasses <- function(passes, lane_set) {
  jobs <- lapply(passes, function(pass) {
    tree <- if (identical(pass, "base")) worktree else head_root
    log <- file.path(out, paste0(pass, "-", gsub(",", "-", lane_set), ".log"))
    list(pass = pass, command = sprintf("cd %s && Rscript %s %s %s %s > %s 2>&1", shQuote(tree), shQuote(child), pass,
      shQuote(file.path(out, pass)), lane_set, shQuote(log)))
  })
  codes <- parallel::mclapply(jobs, function(job) system2("bash", c("-c", shQuote(job$command))), mc.cores = length(jobs),
    mc.preschedule = FALSE)
  stats::setNames(unlist(codes), passes)
}
readRecord <- function(dir, lane) {
  path <- file.path(dir, paste0(lane, ".json"))
  if (!file.exists(path)) return(NULL)
  jsonlite::read_json(path, simplifyVector = FALSE)$steps
}
passes <- if (is.null(base_records)) c("base", "head1", "head0") else c("head1", "head0")
if (!is.null(base_records)) {
  dir.create(file.path(out, "base"), showWarnings = FALSE)
  file.copy(list.files(base_records, full.names = TRUE), file.path(out, "base"), overwrite = TRUE)
}

# --- Round 1: lipid_canonical in every pass, then the mechanism check -----------------------------------------------
codes <- runPasses(passes, "lipid_canonical")
records <- lapply(c(base = "base", head1 = "head1", head0 = "head0"), function(pass) readRecord(file.path(out, pass), "lipid_canonical"))
problems <- character()
for (pass in passes) if (codes[[pass]] != 0L) problems <- c(problems, sprintf("%s exit %d", pass, codes[[pass]]))
for (pass in names(records)) if (is.null(records[[pass]])) problems <- c(problems, sprintf("%s wrote no record", pass))
base_child_sha <- tryCatch(readLines(file.path(out, "base", "source-sha.txt"), warn = FALSE)[[1L]], error = function(e) "<none>")
status("mechanism", "", sprintf("base child source sha %s", base_child_sha))
if (!identical(base_child_sha, base_sha)) problems <- c(problems, sprintf("base child sha %s, not %s", base_child_sha, base_sha))
log_at <- function(steps) steps$norm$outputs[["lipidomics_workflow-norm-norm_log"]] %||% "<absent>"
if (!length(problems)) {
  b <- log_at(records$base); h1 <- log_at(records$head1); h0 <- log_at(records$head0)
  digest_form <- function(x) startsWith(x, "<text> ")
  head1_ok <- if (digest_form(h1)) identical(h1, b) else grepl("(offset: 1 )", h1, fixed = TRUE)
  head0_ok <- if (digest_form(h0)) !identical(h0, h1) else !grepl("(offset:", h0, fixed = TRUE)
  if (!head1_ok) problems <- c(problems, "head1's norm_log at norm does not carry offset 1")
  if (!head0_ok) problems <- c(problems, "head0's norm_log at norm still names an offset")
}
if (length(problems)) {
  status("mechanism", "lipid_canonical", paste("mechanism:", paste(problems, collapse = "; ")))
  finish()
}
status("mechanism", "lipid_canonical", "pass")

# --- Round 2: the metabolomics lanes in every pass -------------------------------------------------------------------
codes <- runPasses(passes, "metab_lc,metab_gc,metab_combined")
for (pass in passes) if (codes[[pass]] != 0L) {
  status("passes", pass, sprintf("mechanism: exit %d", codes[[pass]]))
  finish()
}

# --- The comparisons ------------------------------------------------------------------------------------------------
expected_ui <- readLines(file.path(out, "head1", "expected-recorded-offset.txt"), warn = FALSE)[[1L]]
not_tested <- "Not tested (no value in a compared group): 0"
rows <- list(); digests <- list(); declared <- list(); stops <- character()
classify_frozen <- function(lane, row, metab_range) {
  if (metab_range) return("metab_fixture_4")
  integration <- identical(row$field, "artifacts") && (startsWith(row$key, "integration/") && (row$runner == "<absent>" || row$oracle == "<absent>") ||
    identical(row$key, "<store entries>"))
  if (integration) "integration_artifacts" else "unexplained"
}
for (lane in lanes) {
  omic <- omic_of(lane)
  frozen <- jsonlite::read_json(file.path(frozen_dir, paste0(lane, ".json")), simplifyVector = FALSE)$steps
  base <- readRecord(file.path(out, "base"), lane)
  head1 <- readRecord(file.path(out, "head1"), lane)
  head0 <- readRecord(file.path(out, "head0"), lane)
  if (is.null(base) || is.null(head1) || is.null(head0)) {
    status("records", lane, "mechanism: a pass wrote no record")
    finish()
  }
  steps <- names(base)
  norm_at <- match("norm", steps)
  da_at <- match("da", steps)
  before_norm <- function(step) match(step, steps) < norm_at
  for (step in steps) digests[[length(digests) + 1L]] <- data.frame(lane = lane, step = step,
    frozen = digest::digest(frozen[[step]], algo = "sha256"), base = digest::digest(base[[step]], algo = "sha256"),
    head1 = digest::digest(head1[[step]], algo = "sha256"), head0 = digest::digest(head0[[step]], algo = "sha256"), stringsAsFactors = FALSE)
  # The declared keys: expected and observed at every step.
  offset_key <- paste0(omic, "_workflow-norm-log2_recorded_offset")
  summary_key <- paste0(omic, "_workflow-de-da_summary_stats")
  for (step in steps) {
    i <- match(step, steps)
    base_summary <- base[[step]]$outputs[[summary_key]]
    expected <- c(if (i < norm_at) "<NULL>" else expected_ui,
      if (i >= da_at && is.character(base_summary) && !startsWith(base_summary, "<")) paste0(base_summary, "\n", not_tested) else base_summary %||% "<absent>")
    observed <- c(head1[[step]]$outputs[[offset_key]] %||% "<absent>", head1[[step]]$outputs[[summary_key]] %||% "<absent>")
    declared[[length(declared) + 1L]] <- data.frame(lane = lane, step = step, key = c(paste0("outputs/", offset_key), paste0("outputs/", summary_key)),
      expected = expected, observed_head1 = observed, equal = expected == observed, stringsAsFactors = FALSE)
  }
  declared_keys <- c(offset_key, summary_key)
  # Step 1: base against the frozen file.
  r1 <- compareLane(lane, base, frozen)
  if (!is.null(r1)) {
    metab_range <- startsWith(lane, "metab_") & vapply(r1$step, function(s) !before_norm(s), logical(1))
    r1$class <- vapply(seq_len(nrow(r1)), function(i) classify_frozen(lane, r1[i, ], metab_range[[i]]), "")
    r1$comparison <- "frozen-base"
    rows[[length(rows) + 1L]] <- r1
    if (any(r1$class == "unexplained")) stops <- c(stops, sprintf("%s: base differs from the frozen file outside the classes", lane))
  }
  differing <- sort(unique(if (is.null(r1)) character() else r1$step))
  expected_steps <- sort(reference$step[reference$lane == lane])
  if (!identical(differing, expected_steps)) {
    stops <- c(stops, sprintf("%s: base's differing steps (%s) are not the reference list (%s)", lane, paste(differing, collapse = ","),
      paste(expected_steps, collapse = ",")))
  }
  # Step 2: head at offset 1 against base.
  r2 <- compareLane(lane, head1, base)
  if (!is.null(r2)) {
    r2$class <- ifelse(r2$field == "outputs" & r2$key %in% declared_keys, "declared_change (KIT-172)", "unexplained")
    r2$comparison <- "base-head1"
    rows[[length(rows) + 1L]] <- r2
    if (any(r2$class == "unexplained")) stops <- c(stops, sprintf("%s: head at offset 1 differs from base outside the declared keys", lane))
  }
  # Step 3: head at the default against head at offset 1.
  r3 <- compareLane(lane, head0, head1)
  if (!is.null(r3)) {
    r3$class <- ifelse(vapply(r3$step, before_norm, logical(1)), "unexplained", "plain_log2 (KIT-172)")
    r3$comparison <- "head1-head0"
    rows[[length(rows) + 1L]] <- r3
    if (any(r3$class == "unexplained")) stops <- c(stops, sprintf("%s: the default differs from offset 1 before norm", lane))
  }
  # Base against head at the default: the class of the same key's base-head1 difference, if any, plus plain_log2.
  r4 <- compareLane(lane, head0, base)
  if (!is.null(r4)) {
    declared_here <- if (is.null(r2)) character() else paste(r2$step, r2$field, r2$key)[r2$class == "declared_change (KIT-172)"]
    from_norm <- !vapply(r4$step, before_norm, logical(1))
    is_declared <- paste(r4$step, r4$field, r4$key) %in% declared_here
    r4$class <- ifelse(is_declared & from_norm, "declared_change (KIT-172); plain_log2 (KIT-172)",
      ifelse(is_declared, "declared_change (KIT-172)", ifelse(from_norm, "plain_log2 (KIT-172)", "unexplained")))
    r4$comparison <- "base-head0"
    rows[[length(rows) + 1L]] <- r4
  }
}
control <- do.call(rbind, rows)
lane_steps <- lapply(stats::setNames(lanes, lanes), function(lane) names(readRecord(file.path(out, "base"), lane)))
control$phase <- vapply(seq_len(nrow(control)), function(i) {
  steps <- lane_steps[[control$lane[[i]]]]
  if (match(control$step[[i]], steps) < match("norm", steps)) "before_norm" else "from_norm"
}, "")
control <- control[, c("comparison", "lane", "step", "phase", "field", "key", "kind", "runner", "oracle", "class")]
writeLines(sprintf("# cpu: %s; changed step per lane: norm; runner = the later side, oracle = the earlier side of the comparison", cpu),
  file.path(out, "gui-entry-control.tsv"))
suppressWarnings(utils::write.table(control, file.path(out, "gui-entry-control.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
  append = TRUE))
utils::write.table(do.call(rbind, digests), file.path(out, "gui-entry-digests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
declared <- do.call(rbind, declared)
utils::write.table(declared, file.path(out, "declared-keys.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(declared$equal)) stops <- c(stops, sprintf("%d declared-key values differ from their expected value", sum(!declared$equal)))

# --- The re-frozen files ------------------------------------------------------------------------------------------
dir.create(file.path(out, "refrozen"), showWarnings = FALSE)
for (lane in lanes) file.copy(file.path(head_root, "tests", "testdata", "gui-entry", paste0(lane, ".json")), file.path(out, "refrozen"), overwrite = TRUE)
if (length(stops)) {
  for (reason in stops) status("control", "", paste("stop:", reason))
} else {
  status("control", "", "pass: steps 1 to 3")
}
finish()
