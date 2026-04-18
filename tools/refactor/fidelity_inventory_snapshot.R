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

source(file.path(resolve_refactor_dir(), "fidelity_inventory.R"))

parse_args <- function(args) {
  opts <- list(
    repo_root = ".",
    pretty = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (identical(arg, "--repo-root")) {
      i <- i + 1L
      opts$repo_root <- args[[i]]
    } else if (identical(arg, "--pretty")) {
      opts$pretty <- TRUE
    } else {
      abort("Unknown argument: ", arg)
    }

    i <- i + 1L
  }

  opts
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  snapshot <- inventory_snapshot(opts$repo_root)
  cat(
    jsonlite::toJSON(
      snapshot,
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
