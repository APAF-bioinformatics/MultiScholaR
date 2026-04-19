#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(argv) {
  parsed <- list(
    project_root = NULL,
    test_files = character()
  )

  index <- 1
  while (index <= length(argv)) {
    flag <- argv[[index]]
    if (flag == "--project-root") {
      index <- index + 1
      parsed$project_root <- argv[[index]]
    } else if (flag == "--test-file") {
      index <- index + 1
      parsed$test_files <- c(parsed$test_files, argv[[index]])
    } else {
      stop(sprintf("Unknown argument: %s", flag), call. = FALSE)
    }
    index <- index + 1
  }

  if (is.null(parsed$project_root) || !nzchar(parsed$project_root)) {
    stop("Missing required --project-root", call. = FALSE)
  }

  parsed
}

read_package_name <- function(root) {
  desc_path <- file.path(root, "DESCRIPTION")
  if (!file.exists(desc_path)) {
    stop("Could not find DESCRIPTION under the project root.", call. = FALSE)
  }

  desc <- read.dcf(desc_path)
  if (!"Package" %in% colnames(desc)) {
    stop("DESCRIPTION is missing the Package field.", call. = FALSE)
  }

  package_name <- trimws(desc[1, "Package"])
  if (!nzchar(package_name)) {
    stop("DESCRIPTION Package field is empty.", call. = FALSE)
  }

  package_name
}

normalize_relpath <- function(path, root) {
  normalized_root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (startsWith(normalized_path, paste0(normalized_root, "/"))) {
    return(sub(paste0("^", normalized_root, "/"), "", normalized_path))
  }
  normalized_path
}

summarize_signal_lines <- function(lines, label, limit = 10L) {
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & lines != "NULL"]
  if (!length(lines)) {
    return(list())
  }

  summary_lines <- sprintf("%s: %s", label, utils::head(lines, limit))
  if (length(lines) > limit) {
    summary_lines <- c(summary_lines, sprintf("%s: ... %d additional line(s)", label, length(lines) - limit))
  }

  as.list(summary_lines)
}

expectation_counts <- function(results) {
  counts <- c(success = 0L, failure = 0L, skip = 0L, warning = 0L)
  for (item in results) {
    item_class <- class(item)
    if ("expectation_skip" %in% item_class) {
      counts[["skip"]] <- counts[["skip"]] + 1L
    } else if ("expectation_warning" %in% item_class) {
      counts[["warning"]] <- counts[["warning"]] + 1L
    } else if ("expectation_success" %in% item_class) {
      counts[["success"]] <- counts[["success"]] + 1L
    } else {
      counts[["failure"]] <- counts[["failure"]] + 1L
    }
  }
  counts
}

test_status <- function(counts) {
  if (counts[["failure"]] > 0L) {
    return("failed")
  }
  if (counts[["skip"]] > 0L && counts[["success"]] == 0L) {
    return("skipped")
  }
  "passed"
}

serialize_test_results <- function(results) {
  lapply(results, function(item) {
    counts <- expectation_counts(item$results)
    list(
      test_name = item$test,
      status = test_status(counts),
      expectation_counts = as.list(counts),
      messages = unname(vapply(
        item$results,
        function(result) {
          result$message %||% ""
        },
        character(1)
      ))
    )
  })
}

emit_json <- function(payload) {
  cat(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"), "\n", sep = "")
}

fixture_path <- Sys.getenv("FIDELITY_COVERAGE_FIXTURE_JSON", "")
if (nzchar(fixture_path)) {
  cat(paste(readLines(fixture_path, warn = FALSE), collapse = "\n"), "\n", sep = "")
  quit(save = "no", status = 0)
}

parsed <- parse_args(args)
project_root <- normalizePath(parsed$project_root, winslash = "/", mustWork = TRUE)

selected_test_files <- parsed$test_files
if (!length(selected_test_files)) {
  test_root <- file.path(project_root, "tests", "testthat")
  if (dir.exists(test_root)) {
    selected_test_files <- list.files(
      test_root,
      pattern = "^test-.*\\.R$",
      full.names = TRUE
    )
  }
}

selected_test_files <- normalizePath(selected_test_files, winslash = "/", mustWork = FALSE)
selected_test_files <- selected_test_files[file.exists(selected_test_files)]
selected_test_files <- unique(selected_test_files)

if (Sys.getenv("FIDELITY_COVERAGE_FORCE_UNAVAILABLE", "") %in% c("1", "true", "TRUE")) {
  emit_json(list(
    status = "tool_unavailable",
    tool_status = "forced_unavailable",
    project_root = project_root,
    test_files = vapply(selected_test_files, normalize_relpath, character(1), root = project_root),
    files = list(),
    line_percent = NULL,
    lines_total = 0,
    lines_covered = 0,
    notes = list("Coverage tool availability was disabled via FIDELITY_COVERAGE_FORCE_UNAVAILABLE.")
  ))
  quit(save = "no", status = 0)
}

if (!requireNamespace("covr", quietly = TRUE)) {
  emit_json(list(
    status = "tool_unavailable",
    tool_status = "covr_missing",
    project_root = project_root,
    test_files = vapply(selected_test_files, normalize_relpath, character(1), root = project_root),
    files = list(),
    line_percent = NULL,
    lines_total = 0,
    lines_covered = 0,
    notes = list("covr is not installed in the active R library.")
  ))
  quit(save = "no", status = 0)
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  emit_json(list(
    status = "tool_unavailable",
    tool_status = "pkgload_missing",
    project_root = project_root,
    test_files = vapply(selected_test_files, normalize_relpath, character(1), root = project_root),
    files = list(),
    line_percent = NULL,
    lines_total = 0,
    lines_covered = 0,
    notes = list("pkgload is not installed in the active R library.")
  ))
  quit(save = "no", status = 0)
}

if (!requireNamespace("testthat", quietly = TRUE)) {
  emit_json(list(
    status = "tool_unavailable",
    tool_status = "testthat_missing",
    project_root = project_root,
    test_files = vapply(selected_test_files, normalize_relpath, character(1), root = project_root),
    files = list(),
    line_percent = NULL,
    lines_total = 0,
    lines_covered = 0,
    notes = list("testthat is not installed in the active R library.")
  ))
  quit(save = "no", status = 0)
}

source_root <- file.path(project_root, "R")
source_files <- character()
if (dir.exists(source_root)) {
  source_files <- list.files(
    source_root,
    pattern = "\\.R$",
    full.names = TRUE
  )
}

if (!length(source_files)) {
  emit_json(list(
    status = "runner_error",
    tool_status = "no_source_files",
    project_root = project_root,
    test_files = vapply(selected_test_files, normalize_relpath, character(1), root = project_root),
    files = list(),
    line_percent = NULL,
    lines_total = 0,
    lines_covered = 0,
    notes = list("No R source files were found under the project root.")
  ))
  quit(save = "no", status = 0)
}

if (!length(selected_test_files)) {
  emit_json(list(
    status = "runner_error",
    tool_status = "no_test_files",
    project_root = project_root,
    test_files = list(),
    files = list(),
    line_percent = NULL,
    lines_total = 0,
    lines_covered = 0,
    notes = list("No testthat files were selected for coverage capture.")
  ))
  quit(save = "no", status = 0)
}

package_name <- read_package_name(project_root)
project_source_prefix <- paste0(normalizePath(source_root, winslash = "/", mustWork = FALSE), "/")

result <- tryCatch(
  {
    warnings <- character()
    messages <- character()
    output <- character()
    reporter_results <- list()

    with_handlers <- function(expr) {
      withCallingHandlers(
        expr,
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          messages <<- c(messages, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
    }

    output <- with_handlers(capture.output({
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(project_root)

      old_device <- getOption("device")
      options(device = function(...) grDevices::pdf(file = tempfile(fileext = ".pdf")))
      on.exit(options(device = old_device), add = TRUE)

      pkgload::load_all(project_root, quiet = TRUE, export_all = TRUE)
      namespace_env <- if (isNamespaceLoaded(package_name)) {
        asNamespace(package_name)
      } else {
        pkgload::ns_env(package_name)
      }

      covr:::trace_environment(namespace_env)
      on.exit({
        covr:::reset_traces()
        covr:::clear_counters()
      }, add = TRUE)

      old_covr_env <- Sys.getenv("R_COVR", unset = NA_character_)
      Sys.setenv(R_COVR = "true")
      on.exit({
        if (is.na(old_covr_env)) {
          Sys.unsetenv("R_COVR")
        } else {
          Sys.setenv(R_COVR = old_covr_env)
        }
      }, add = TRUE)

      invisible(lapply(selected_test_files, function(test_file) {
        file_results <- testthat::test_file(test_file, reporter = testthat::ListReporter$new())
        reporter_results <<- c(reporter_results, file_results)
        NULL
      }))
    }, type = "output"))

    coverage <- covr:::as_coverage(covr:::.counters)
    coverage_df <- as.data.frame(coverage)

    filename_col <- intersect(c("filename", "file"), names(coverage_df))
    line_col <- intersect(c("line", "srcref"), names(coverage_df))
    value_col <- intersect(c("value", "count"), names(coverage_df))
    first_line_col <- intersect(c("first_line"), names(coverage_df))
    last_line_col <- intersect(c("last_line"), names(coverage_df))
    has_range_columns <- length(first_line_col) && length(last_line_col)
    if (!length(filename_col) || !length(value_col) || (!length(line_col) && !has_range_columns)) {
      stop("covr output schema is missing filename, value, and line-range columns.", call. = FALSE)
    }

    filename_col <- filename_col[[1]]
    value_col <- value_col[[1]]
    if (length(line_col)) {
      line_col <- line_col[[1]]
    } else {
      first_line_col <- first_line_col[[1]]
      last_line_col <- last_line_col[[1]]
    }

    if (length(line_col)) {
      coverage_df <- coverage_df[!is.na(coverage_df[[line_col]]), , drop = FALSE]
    } else {
      coverage_df <- coverage_df[
        !is.na(coverage_df[[first_line_col]]) & !is.na(coverage_df[[last_line_col]]),
        ,
        drop = FALSE
      ]
    }
    coverage_df[[filename_col]] <- normalizePath(coverage_df[[filename_col]], winslash = "/", mustWork = FALSE)
    coverage_df <- coverage_df[startsWith(coverage_df[[filename_col]], project_source_prefix), , drop = FALSE]
    grouped <- split(coverage_df, coverage_df[[filename_col]])

    files <- lapply(names(grouped), function(filename) {
      rows <- grouped[[filename]]
      if (length(line_col)) {
        all_lines <- sort(unique(as.integer(rows[[line_col]])))
        covered_lines <- sort(unique(as.integer(rows[[line_col]][rows[[value_col]] > 0])))
      } else {
        row_ranges <- Map(
          function(first_line, last_line) {
            seq.int(as.integer(first_line), as.integer(last_line))
          },
          rows[[first_line_col]],
          rows[[last_line_col]]
        )
        all_lines <- sort(unique(unlist(row_ranges, use.names = FALSE)))
        covered_lines <- sort(unique(unlist(
          row_ranges[rows[[value_col]] > 0],
          use.names = FALSE
        )))
      }
      uncovered_lines <- sort(setdiff(all_lines, covered_lines))
      lines_total <- length(all_lines)
      lines_covered <- length(covered_lines)
      line_percent <- if (lines_total > 0) (lines_covered / lines_total) * 100 else NULL

      list(
        file_path = normalize_relpath(filename, project_root),
        line_percent = line_percent,
        lines_total = lines_total,
        lines_covered = lines_covered,
        covered_lines = covered_lines,
        uncovered_lines = uncovered_lines
      )
    })

    total_lines <- sum(vapply(files, function(item) item$lines_total, numeric(1)))
    covered_lines <- sum(vapply(files, function(item) item$lines_covered, numeric(1)))
    line_percent <- if (total_lines > 0) (covered_lines / total_lines) * 100 else NULL
    tests <- serialize_test_results(reporter_results)
    test_statuses <- unname(vapply(tests, function(item) item$status, character(1)))
    failure_count <- sum(test_statuses == "failed")
    skip_count <- sum(test_statuses == "skipped")
    failed_test_names <- unname(vapply(
      tests[test_statuses == "failed"],
      function(item) item$test_name,
      character(1)
    ))
    capture_status <- if (failure_count > 0L) "test_failures" else "completed"
    tool_status <- if (failure_count > 0L) "testthat_failures" else "ok"

    list(
      status = capture_status,
      tool_status = tool_status,
      project_root = project_root,
      test_files = vapply(selected_test_files, normalize_relpath, character(1), root = project_root),
      test_count = length(tests),
      failure_count = failure_count,
      skip_count = skip_count,
      failed_tests = as.list(failed_test_names),
      tests = tests,
      files = files,
      line_percent = line_percent,
      lines_total = total_lines,
      lines_covered = covered_lines,
      notes = c(
        list("Coverage captured with pkgload::load_all + covr namespace tracing + testthat::test_file."),
        if (failure_count > 0L) {
          c(
            list(sprintf("Selected coverage test subset had %d failing test(s).", failure_count)),
            summarize_signal_lines(failed_test_names, "failed_test")
          )
        } else {
          list()
        },
        summarize_signal_lines(warnings, "warning"),
        summarize_signal_lines(messages, "message"),
        summarize_signal_lines(output, "output")
      )
    )
  },
  error = function(err) {
    list(
      status = "runner_error",
      tool_status = "covr_error",
      project_root = project_root,
      test_files = vapply(selected_test_files, normalize_relpath, character(1), root = project_root),
      files = list(),
      line_percent = NULL,
      lines_total = 0,
      lines_covered = 0,
      notes = list(conditionMessage(err))
    )
  }
)

emit_json(result)
