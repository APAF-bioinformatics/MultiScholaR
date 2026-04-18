#!/usr/bin/env Rscript

parse_args <- function(args) {
  opts <- list(
    repo_root = ".",
    output = NULL,
    verbose = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg == "--repo-root") {
      i <- i + 1
      opts$repo_root <- args[[i]]
    } else if (arg == "--output") {
      i <- i + 1
      opts$output <- args[[i]]
    } else if (arg == "--verbose") {
      opts$verbose <- TRUE
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }

    i <- i + 1
  }

  opts
}

collect_scan_files <- function(repo_root) {
  search_roots <- c("R", "tests", "dev", "inst", "docs", "Workbooks")
  exts <- c(".R", ".r", ".Rmd", ".rmd", ".qmd", ".md")
  files <- character()

  for (root in search_roots) {
    abs_root <- file.path(repo_root, root)
    if (!dir.exists(abs_root)) {
      next
    }

    root_files <- list.files(abs_root, recursive = TRUE, full.names = TRUE)
    root_files <- root_files[file.info(root_files)$isdir %in% FALSE]
    keep <- Reduce(`|`, lapply(exts, function(ext) endsWith(root_files, ext)))
    files <- c(files, root_files[keep])
  }

  files <- unique(files)

  excluded <- grepl("(^|/)tools/refactor/staging/", files) |
    grepl("(^|/)renv/", files) |
    grepl("(^|/)docs/.+_files/", files)

  files[!excluded]
}

rel_path <- function(path, repo_root) {
  sub(paste0("^", normalizePath(repo_root, winslash = "/", mustWork = TRUE), "/?"), "", normalizePath(path, winslash = "/", mustWork = TRUE))
}

classify_scope <- function(path) {
  if (grepl("^R/", path)) return("package")
  if (grepl("^tests/", path)) return("tests")
  if (grepl("^dev/", path)) return("dev")
  if (grepl("^inst/", path)) return("inst")
  if (grepl("^docs/", path)) return("docs")
  if (grepl("^Workbooks/", path)) return("workbooks")
  "other"
}

scan_patterns <- function() {
  list(
    list(
      id = "source_r_file",
      severity = "high",
      description = "Direct source() of R/ file",
      regex = "source\\s*\\(\\s*['\"]R/[^'\"]+\\.R['\"]"
    ),
    list(
      id = "parse_r_file",
      severity = "high",
      description = "Direct parse(file=...) of R/ file",
      regex = "parse\\s*\\(.*file\\s*=\\s*['\"]R/[^'\"]+\\.R['\"]"
    ),
    list(
      id = "readlines_r_file",
      severity = "medium",
      description = "Direct readLines() of R/ file",
      regex = "readLines\\s*\\(\\s*['\"]R/[^'\"]+\\.R['\"]"
    ),
    list(
      id = "file_exists_r_file",
      severity = "medium",
      description = "Direct file.exists() check against R/ file",
      regex = "file\\.exists\\s*\\(\\s*['\"]R/[^'\"]+\\.R['\"]"
    ),
    list(
      id = "literal_r_file",
      severity = "low",
      description = "Literal R/ filename reference",
      regex = "['\"]R/[^'\"]+\\.R['\"]"
    )
  )
}

find_matches <- function(file, repo_root) {
  lines <- readLines(file, warn = FALSE)
  rel <- rel_path(file, repo_root)
  scope <- classify_scope(rel)
  patterns <- scan_patterns()
  rows <- list()

  for (i in seq_along(lines)) {
    line <- lines[[i]]

    for (pat in patterns) {
      if (!grepl(pat$regex, line, perl = TRUE)) {
        next
      }

      match <- regmatches(line, regexpr(pat$regex, line, perl = TRUE))
      rows[[length(rows) + 1L]] <- data.frame(
        scope = scope,
        file = rel,
        line = i,
        severity = pat$severity,
        pattern = pat$id,
        description = pat$description,
        snippet = trimws(line),
        match = if (length(match)) match[[1]] else "",
        stringsAsFactors = FALSE
      )
      break
    }
  }

  if (!length(rows)) {
    return(NULL)
  }

  do.call(rbind, rows)
}

format_markdown <- function(df, repo_root) {
  if (!nrow(df)) {
    return(c(
      "# Refactor Coupling Audit",
      "",
      "No filename-coupled references to `R/*.R` were found in the scanned text files."
    ))
  }

  sev_order <- c(high = 1L, medium = 2L, low = 3L)
  df <- df[order(unname(sev_order[df$severity]), df$scope, df$file, df$line), , drop = FALSE]

  summary_counts <- aggregate(
    rep(1L, nrow(df)),
    by = list(severity = df$severity, scope = df$scope),
    FUN = sum
  )
  names(summary_counts)[[3]] <- "count"
  summary_counts <- summary_counts[order(unname(sev_order[summary_counts$severity]), summary_counts$scope), , drop = FALSE]

  lines <- c(
    "# Refactor Coupling Audit",
    "",
    "This report flags filename-coupled references to `R/*.R` that can break wave-based splitting.",
    "",
    "## Summary",
    ""
  )

  for (i in seq_len(nrow(summary_counts))) {
    row <- summary_counts[i, , drop = FALSE]
    lines <- c(lines, sprintf("- `%s` / `%s`: %d", row$severity, row$scope, row$count))
  }

  lines <- c(
    lines,
    "",
    "## Findings",
    ""
  )

  for (i in seq_len(nrow(df))) {
    row <- df[i, , drop = FALSE]
    lines <- c(
      lines,
      sprintf(
        "- `%s` `%s` [%s:%d](%s:%d): %s",
        row$severity,
        row$scope,
        row$file,
        row$line,
        file.path(repo_root, row$file),
        row$line,
        row$snippet
      )
    )
  }

  lines
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_root <- normalizePath(opts$repo_root, mustWork = TRUE)
  files <- collect_scan_files(repo_root)

  findings <- lapply(files, find_matches, repo_root = repo_root)
  findings <- Filter(Negate(is.null), findings)
  df <- if (length(findings)) do.call(rbind, findings) else data.frame()
  report <- format_markdown(df, repo_root)

  if (!is.null(opts$output)) {
    output_path <- file.path(repo_root, opts$output)
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(report, con = output_path)
    if (isTRUE(opts$verbose)) {
      message("Wrote audit report: ", opts$output)
    }
  } else {
    cat(paste(report, collapse = "\n"), "\n")
  }
}

main()
