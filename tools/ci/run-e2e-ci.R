#!/usr/bin/env Rscript

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_file_arg[[1L]]), mustWork = FALSE)
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
  if (dir.exists(file.path(repo_root, "tools", "ci"))) {
    setwd(repo_root)
  }
}

default_args <- list(
  lane = character(),
  filter = character(),
  reporter = "summary",
  artifact_dir = file.path("tests", "testthat", "_e2e_artifacts"),
  manifest = file.path("tests", "testdata", "e2e", "manifest.json"),
  dry_run = FALSE
)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

usage <- function() {
  cat(paste(
    "Usage: Rscript tools/ci/run-e2e-ci.R [options]",
    "",
    "Options:",
    "  --lane <lane_id>      May be repeated; comma-separated accepted",
    "  --filter <regex>      May be repeated; comma-separated accepted",
    "  --reporter <testthat_reporter>",
    "  --artifact-dir <path>",
    "  --manifest <tests/testdata/e2e/manifest.json>",
    "  --dry-run <true|false>  Write the run manifest without executing tests",
    "  --help",
    sep = "\n"
  ))
}

split_values <- function(x) {
  x <- unlist(strsplit(as.character(x), "[,\n\r]+"), use.names = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}

parse_bool <- function(x) {
  isTRUE(tolower(as.character(x)) %in% c("1", "true", "yes", "y"))
}

parse_args <- function(argv) {
  args <- default_args
  idx <- 1L
  while (idx <= length(argv)) {
    token <- argv[[idx]]
    if (identical(token, "--help")) {
      usage()
      quit(status = 0L)
    }
    if (!startsWith(token, "--")) {
      stop(sprintf("Unexpected positional argument: %s", token), call. = FALSE)
    }

    key_value <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", key_value[[1L]], fixed = TRUE)
    if (!key %in% names(args)) {
      stop(sprintf("Unknown option: %s", token), call. = FALSE)
    }
    if (length(key_value) == 2L) {
      value <- key_value[[2L]]
    } else {
      idx <- idx + 1L
      if (idx > length(argv)) {
        stop(sprintf("Missing value for option: %s", token), call. = FALSE)
      }
      value <- argv[[idx]]
    }
    if (key %in% c("lane", "filter")) {
      args[[key]] <- c(args[[key]], split_values(value))
    } else if (identical(key, "dry_run")) {
      args[[key]] <- parse_bool(value)
    } else {
      args[[key]] <- value
    }
    idx <- idx + 1L
  }
  args
}

utc_now <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

read_lanes <- function(path) {
  raw <- jsonlite::read_json(path, simplifyVector = FALSE)
  lanes <- raw$lanes
  names(lanes) <- vapply(lanes, `[[`, character(1), "lane_id")
  lanes
}

write_manifest <- function(path, payload) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$artifact_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(MULTISCHOLAR_E2E_ARTIFACT_DIR = args$artifact_dir)

lanes <- read_lanes(args$manifest)
selected_lanes <- list()
if (length(args$lane) > 0L) {
  missing <- setdiff(args$lane, names(lanes))
  if (length(missing) > 0L) {
    stop(sprintf("Unknown E2E lane(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  selected_lanes <- lanes[args$lane]
}

lane_filters <- unique(vapply(selected_lanes, function(lane) lane$test_filter %||% NA_character_, character(1)))
lane_filters <- lane_filters[!is.na(lane_filters) & nzchar(lane_filters)]
filters <- unique(c(args$filter, lane_filters))
if (length(filters) == 0L) {
  stop("No E2E lanes or filters selected", call. = FALSE)
}
test_filter <- sprintf("(%s)", paste(filters, collapse = "|"))

run_started_at <- utc_now()
test_error <- NULL
tryCatch(
  {
    if (!isTRUE(args$dry_run)) {
      Sys.setenv(NOT_CRAN = "true")
      pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
      testthat::test_dir("tests/testthat", filter = test_filter, reporter = args$reporter)
    }
  },
  error = function(err) {
    test_error <<- conditionMessage(err)
  }
)

result <- if (is.null(test_error)) "passed" else "failed"
payload <- list(
  schema_version = "1.0.0",
  generated_at = utc_now(),
  started_at = run_started_at,
  finished_at = utc_now(),
  invocation = args,
  result = result,
  failure_reason = test_error,
  test_filter = test_filter,
  lanes = selected_lanes,
  filters = as.list(filters),
  artifact_dir = args$artifact_dir,
  run_artifacts = as.list(c(
    file.path(args$artifact_dir, "e2e-run-manifest.json")
  ))
)
write_manifest(file.path(args$artifact_dir, "e2e-run-manifest.json"), payload)
message(sprintf(
  "E2E-CI result=%s lanes=%d filters=%s artifact_dir=%s",
  result,
  length(selected_lanes),
  paste(filters, collapse = ","),
  args$artifact_dir
))
if (!is.null(test_error)) {
  quit(status = 1L)
}
