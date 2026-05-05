#!/usr/bin/env Rscript

default_args <- list(
  artifact_dir = file.path("tests", "testthat", "_module_ci_artifacts"),
  input = NULL,
  output_json = NULL,
  output_md = NULL
)

usage <- function() {
  cat(
    paste(
      "Usage: Rscript tools/ci/module-ci-scorecard.R [options]",
      "",
      "Options:",
      "  --artifact-dir <path>",
      "  --input <module-ci-run-manifest.json>",
      "  --output-json <path>",
      "  --output-md <path>",
      "  --help",
      sep = "\n"
    )
  )
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
      args[[key]] <- key_value[[2L]]
    } else {
      idx <- idx + 1L
      if (idx > length(argv)) {
        stop(sprintf("Missing value for option: %s", token), call. = FALSE)
      }
      args[[key]] <- argv[[idx]]
    }
    idx <- idx + 1L
  }
  args
}

as_chr <- function(x) {
  as.character(unlist(x, use.names = FALSE))
}

utc_now <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

scenario_status <- function(scenario, run_result) {
  status <- scenario$result
  if (is.null(status) || !nzchar(status)) {
    status <- run_result
  }
  if (identical(status, "skipped")) {
    return("skipped")
  }
  if (identical(status, "passed")) {
    return("passed")
  }
  if (identical(status, "failed")) {
    return("failed")
  }
  "unknown"
}

scenario_failure_reason <- function(scenario, status, run_failure_reason) {
  reason <- scenario$failure_reason
  if (is.null(reason) || length(reason) == 0L || !nzchar(reason)) {
    reason <- run_failure_reason
  }
  if ((is.null(reason) || length(reason) == 0L || !nzchar(reason)) && identical(status, "skipped")) {
    reason <- "skipped; not counted as passed"
  }
  if (is.null(reason) || length(reason) == 0L) {
    return("")
  }
  as.character(reason)
}

scorecard_rows <- function(manifest) {
  lapply(manifest$scenarios, function(scenario) {
    status <- scenario_status(scenario, manifest$result)
    list(
      scenario_id = scenario$scenario_id,
      ticket_id = scenario$ticket_id,
      omic = scenario$omic,
      module = scenario$module,
      runtime = scenario$runtime,
      ci_lane = scenario$ci_lane,
      result = status,
      artifacts = as.list(unique(c(as_chr(scenario$artifacts), as_chr(scenario$run_artifacts)))),
      failure_reason = scenario_failure_reason(scenario, status, manifest$failure_reason)
    )
  })
}

count_statuses <- function(rows) {
  statuses <- c("passed", "failed", "skipped", "unknown")
  counts <- stats::setNames(as.list(rep(0L, length(statuses))), statuses)
  observed <- table(vapply(rows, `[[`, character(1), "result"))
  for (status in names(observed)) {
    counts[[status]] <- as.integer(observed[[status]])
  }
  counts
}

overall_result <- function(counts) {
  if (counts$failed > 0L || counts$unknown > 0L) {
    return("failed")
  }
  if (counts$skipped > 0L) {
    return("skipped")
  }
  "passed"
}

markdown_escape <- function(x) {
  x <- gsub("\\|", "\\\\|", x)
  gsub("\n", " ", x, fixed = TRUE)
}

write_markdown <- function(path, scorecard) {
  rows <- scorecard$scenarios
  lines <- c(
    "# Module CI Scorecard",
    "",
    sprintf("- Generated: `%s`", scorecard$generated_at),
    sprintf("- Result: `%s`", scorecard$result),
    sprintf("- Passed: `%s`", scorecard$counts$passed),
    sprintf("- Failed: `%s`", scorecard$counts$failed),
    sprintf("- Skipped: `%s`", scorecard$counts$skipped),
    sprintf("- Unknown: `%s`", scorecard$counts$unknown),
    "",
    "| Scenario | Omic | Module | Runtime | Result | Failure reason |",
    "| --- | --- | --- | --- | --- | --- |"
  )

  for (row in rows) {
    lines <- c(lines, sprintf(
      "| `%s` | `%s` | `%s` | `%s` | `%s` | %s |",
      markdown_escape(row$scenario_id),
      markdown_escape(row$omic),
      markdown_escape(row$module),
      markdown_escape(row$runtime),
      markdown_escape(row$result),
      markdown_escape(row$failure_reason)
    ))
  }

  writeLines(lines, path, useBytes = TRUE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$input)) {
  args$input <- file.path(args$artifact_dir, "module-ci-run-manifest.json")
}
if (is.null(args$output_json)) {
  args$output_json <- file.path(args$artifact_dir, "module-ci-scorecard.json")
}
if (is.null(args$output_md)) {
  args$output_md <- file.path(args$artifact_dir, "module-ci-scorecard.md")
}

if (!file.exists(args$input)) {
  stop(sprintf("Run manifest not found: %s", args$input), call. = FALSE)
}

manifest <- jsonlite::read_json(args$input, simplifyVector = FALSE)
rows <- scorecard_rows(manifest)
counts <- count_statuses(rows)
scorecard <- list(
  schema_version = "1.0.0",
  generated_at = utc_now(),
  source_manifest = args$input,
  artifact_dir = args$artifact_dir,
  result = overall_result(counts),
  counts = counts,
  scenarios = rows,
  artifact_retention = manifest$artifact_retention
)

dir.create(dirname(args$output_json), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(scorecard, args$output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
write_markdown(args$output_md, scorecard)

message(sprintf(
  "module-CI scorecard result=%s passed=%d failed=%d skipped=%d unknown=%d output=%s",
  scorecard$result,
  counts$passed,
  counts$failed,
  counts$skipped,
  counts$unknown,
  args$output_json
))

if (identical(scorecard$result, "failed")) {
  quit(status = 1L)
}
