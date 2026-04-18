#!/usr/bin/env Rscript

resolve_refactor_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }

  frames <- sys.frames()
  for (idx in rev(seq_along(frames))) {
    ofile <- frames[[idx]]$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }

  normalizePath("tools/refactor", mustWork = TRUE)
}

refactor_dir <- resolve_refactor_dir()
source(file.path(refactor_dir, "fidelity_inventory.R"))

parse_args <- function(args) {
  opts <- list(
    case_file = NULL,
    case_id = NULL,
    repo_root = ".",
    side = NULL,
    pretty = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (identical(arg, "--case-file")) {
      i <- i + 1L
      opts$case_file <- args[[i]]
    } else if (identical(arg, "--case-id")) {
      i <- i + 1L
      opts$case_id <- args[[i]]
    } else if (identical(arg, "--repo-root")) {
      i <- i + 1L
      opts$repo_root <- args[[i]]
    } else if (identical(arg, "--side")) {
      i <- i + 1L
      opts$side <- args[[i]]
    } else if (identical(arg, "--pretty")) {
      opts$pretty <- TRUE
    } else {
      abort("Unknown argument: ", arg)
    }

    i <- i + 1L
  }

  if (is.null(opts$case_file) || is.null(opts$case_id) || is.null(opts$side)) {
    abort("--case-file, --case-id, and --side are required")
  }

  if (!opts$side %in% c("baseline", "target")) {
    abort("--side must be 'baseline' or 'target'")
  }

  opts
}

read_case_file <- function(path) {
  payload <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  cases <- payload$cases
  if (is.null(cases) || !length(cases)) {
    abort("Behavior case file does not contain any cases: ", path)
  }
  cases
}

find_case <- function(cases, case_id) {
  matches <- Filter(function(case) identical(case$id, case_id), cases)
  if (!length(matches)) {
    abort("Behavior case not found: ", case_id)
  }
  if (length(matches) != 1L) {
    abort("Behavior case id is not unique: ", case_id)
  }
  matches[[1]]
}

case_files_for_side <- function(case, side) {
  field <- if (identical(side, "baseline")) "baseline_files" else "target_files"
  files <- case[[field]]
  if (is.null(files) || !length(files)) {
    files <- case$files
  }

  if (is.null(files)) {
    return(character())
  }

  unname(unlist(files))
}

execute_code_lines <- function(lines, env) {
  if (is.null(lines)) {
    return(invisible(NULL))
  }

  if (is.character(lines) && length(lines) == 1L) {
    lines <- list(lines)
  }

  for (line in unname(unlist(lines))) {
    eval(parse(text = line), envir = env)
  }
}

serialize_value <- function(value) {
  jsonlite::serializeJSON(value, digits = NA)
}

run_case <- function(case, repo_root, side) {
  env <- globalenv()
  files <- case_files_for_side(case, side)

  warnings <- character()
  messages <- character()
  output <- character()
  error_info <- NULL
  raw_result <- NULL
  normalized_value <- NULL

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

  output <- tryCatch(
    with_handlers(capture.output({
      for (relative_path in files) {
        absolute_path <- file.path(repo_root, relative_path)
        if (!file.exists(absolute_path)) {
          abort("Behavior replay source file does not exist for ", side, ": ", relative_path)
        }
        sys.source(absolute_path, envir = env)
      }

      execute_code_lines(case$setup_code, env)
      raw_result <- eval(parse(text = case$call_expr), envir = env)
      assign("result", raw_result, envir = env)
      normalized_value <- eval(parse(text = case$normalize_expr %||% "result"), envir = env)
    }, type = "output")),
    error = function(e) {
      error_info <<- list(
        message = conditionMessage(e),
        class = unname(class(e))
      )
      character()
    }
  )

  normalized_serialized <- NULL
  raw_serialized <- NULL
  if (is.null(error_info)) {
    raw_serialized <- serialize_value(raw_result)
    normalized_serialized <- serialize_value(normalized_value)
  }

  list(
    case_id = case$id,
    family = case$family %||% "pure_helper",
    entity_key = case$entity_key %||% NA_character_,
    side = side,
    status = if (is.null(error_info)) "ok" else "error",
    raw_result_json = raw_serialized,
    normalized_result_json = normalized_serialized,
    warnings = warnings,
    messages = messages,
    output = output,
    error = error_info
  )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_root <- normalizePath(opts$repo_root, mustWork = TRUE)
  case_file <- normalizePath(opts$case_file, mustWork = TRUE)
  cases <- read_case_file(case_file)
  case <- find_case(cases, opts$case_id)
  result <- run_case(case, repo_root, opts$side)

  cat(
    jsonlite::toJSON(
      result,
      auto_unbox = TRUE,
      pretty = isTRUE(opts$pretty),
      null = "null",
      na = "null"
    )
  )
  cat("\n")
}

tryCatch(
  main(),
  error = function(e) {
    writeLines(conditionMessage(e), con = stderr())
    quit(status = 1)
  }
)
