#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

parse_args <- function(args) {
  opts <- list(
    repo_root = ".",
    root = "R",
    manifest = NULL,
    output = NULL,
    ideal_min = 150L,
    ideal_max = 500L,
    acceptable_max = 800L,
    soft_max = 1000L
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg == "--repo-root") {
      i <- i + 1
      opts$repo_root <- args[[i]]
    } else if (arg == "--root") {
      i <- i + 1
      opts$root <- args[[i]]
    } else if (arg == "--manifest") {
      i <- i + 1
      opts$manifest <- args[[i]]
    } else if (arg == "--output") {
      i <- i + 1
      opts$output <- args[[i]]
    } else if (arg == "--ideal-min") {
      i <- i + 1
      opts$ideal_min <- as.integer(args[[i]])
    } else if (arg == "--ideal-max") {
      i <- i + 1
      opts$ideal_max <- as.integer(args[[i]])
    } else if (arg == "--acceptable-max") {
      i <- i + 1
      opts$acceptable_max <- as.integer(args[[i]])
    } else if (arg == "--soft-max") {
      i <- i + 1
      opts$soft_max <- as.integer(args[[i]])
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }

    i <- i + 1
  }

  opts
}

rel_path <- function(path, repo_root) {
  sub(
    paste0("^", normalizePath(repo_root, winslash = "/", mustWork = TRUE), "/?"),
    "",
    normalizePath(path, winslash = "/", mustWork = TRUE)
  )
}

line_count <- function(path) {
  length(readLines(path, warn = FALSE))
}

classify_budget <- function(loc, ideal_min, ideal_max, acceptable_max, soft_max) {
  if (loc < ideal_min) {
    return("small")
  }
  if (loc <= ideal_max) {
    return("ideal")
  }
  if (loc <= acceptable_max) {
    return("acceptable")
  }
  if (loc <= soft_max) {
    return("soft_cap")
  }
  "oversized"
}

collect_files <- function(repo_root, root, manifest_path = NULL) {
  root_path <- file.path(repo_root, root)
  if (!dir.exists(root_path)) {
    stop("Root does not exist: ", root, call. = FALSE)
  }

  files <- list.files(root_path, pattern = "\\.R$", full.names = TRUE)

  if (is.null(manifest_path)) {
    return(files)
  }

  manifest <- yaml::read_yaml(file.path(repo_root, manifest_path))
  targets <- unique(vapply(Filter(function(entry) {
    !is.null(entry$target)
  }, manifest$entries), function(entry) entry$target, character(1)))
  targets <- file.path(repo_root, targets)
  targets[file.exists(targets)]
}

build_table <- function(files, repo_root, opts) {
  if (!length(files)) {
    return(data.frame())
  }

  df <- data.frame(
    file = vapply(files, rel_path, character(1), repo_root = repo_root),
    loc = vapply(files, line_count, integer(1)),
    stringsAsFactors = FALSE
  )

  df$budget <- vapply(
    df$loc,
    classify_budget,
    character(1),
    ideal_min = opts$ideal_min,
    ideal_max = opts$ideal_max,
    acceptable_max = opts$acceptable_max,
    soft_max = opts$soft_max
  )

  df <- df[order(df$loc, decreasing = TRUE), , drop = FALSE]
  rownames(df) <- NULL
  df
}

format_markdown <- function(df, opts, title) {
  if (!nrow(df)) {
    return(c(
      sprintf("# %s", title),
      "",
      "No matching `.R` files were found."
    ))
  }

  levels <- c("oversized", "soft_cap", "acceptable", "ideal", "small")
  counts <- table(factor(df$budget, levels = levels))

  lines <- c(
    sprintf("# %s", title),
    "",
    "Budget bands:",
    sprintf("- `ideal`: %d-%d LOC", opts$ideal_min, opts$ideal_max),
    sprintf("- `acceptable`: %d-%d LOC", opts$ideal_max + 1L, opts$acceptable_max),
    sprintf("- `soft_cap`: %d-%d LOC", opts$acceptable_max + 1L, opts$soft_max),
    sprintf("- `oversized`: > %d LOC", opts$soft_max),
    sprintf("- `small`: < %d LOC", opts$ideal_min),
    "",
    "## Summary",
    ""
  )

  for (name in levels) {
    lines <- c(lines, sprintf("- `%s`: %d", name, counts[[name]]))
  }

  lines <- c(lines, "", "## Largest Files", "")

  for (i in seq_len(min(40L, nrow(df)))) {
    row <- df[i, , drop = FALSE]
    lines <- c(lines, sprintf("- `%s` %d LOC [%s](%s)", row$budget, row$loc, row$file, file.path(normalizePath(".", mustWork = TRUE), row$file)))
  }

  lines
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_root <- normalizePath(opts$repo_root, mustWork = TRUE)
  files <- collect_files(repo_root, opts$root, opts$manifest)
  df <- build_table(files, repo_root, opts)

  title <- if (is.null(opts$manifest)) {
    "File Size Audit"
  } else {
    sprintf("File Size Audit: %s", basename(opts$manifest))
  }

  report <- format_markdown(df, opts, title)

  if (!is.null(opts$output)) {
    output_path <- file.path(repo_root, opts$output)
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(report, con = output_path)
  } else {
    cat(paste(report, collapse = "\n"), "\n")
  }
}

main()
