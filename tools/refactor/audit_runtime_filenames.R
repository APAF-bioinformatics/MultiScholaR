#!/usr/bin/env Rscript

parse_args <- function(args) {
  opts <- list(
    repo_root = ".",
    root = "R",
    output = NULL,
    check = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--check") {
      opts$check <- TRUE
      i <- i + 1L
      next
    }
    if (!arg %in% c("--repo-root", "--root", "--output")) {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for ", arg, call. = FALSE)
    }
    key <- gsub("-", "_", sub("^--", "", arg))
    opts[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }

  opts
}

classify_filename <- function(filename) {
  core <- c("app_ui.R", "app_server.R", "run_app.R", "launch_app.R", "zzz.R")
  function_domains <- c(
    "general", "github", "lipid", "metab", "multiomics", "omics",
    "pept", "phospho", "prot"
  )
  module_domains <- c(
    "lipid", "lipidomics", "metab", "metabolomics", "prot", "proteomics"
  )

  if (filename %in% core) {
    return(c(family = "golem_entrypoint", violation = ""))
  }

  function_match <- regexec(
    "^func_([a-z][a-z0-9]*)_([a-z0-9]+(?:_[a-z0-9]+)*)[.]R$",
    filename,
    perl = TRUE
  )
  function_parts <- regmatches(filename, function_match)[[1L]]
  if (length(function_parts)) {
    if (!function_parts[[2L]] %in% function_domains) {
      return(c(
        family = "function",
        violation = paste0("unrecognized function domain: ", function_parts[[2L]])
      ))
    }
    return(c(family = "function", violation = ""))
  }

  module_match <- regexec(
    "^mod_([a-z][a-z0-9]*)(?:_[a-z0-9]+)*[.]R$",
    filename,
    perl = TRUE
  )
  module_parts <- regmatches(filename, module_match)[[1L]]
  if (length(module_parts)) {
    if (!module_parts[[2L]] %in% module_domains) {
      return(c(
        family = "module",
        violation = paste0("unrecognized module domain: ", module_parts[[2L]])
      ))
    }
    return(c(family = "module", violation = ""))
  }

  if (grepl("^utils_[a-z0-9]+(?:_[a-z0-9]+)*[.]R$", filename, perl = TRUE)) {
    return(c(family = "utility", violation = ""))
  }

  c(family = "unclassified", violation = "does not match a runtime filename form")
}

read_collate <- function(description_path) {
  description <- read.dcf(description_path)
  if (!"Collate" %in% colnames(description)) {
    return(character())
  }
  entries <- strsplit(description[[1L, "Collate"]], "[[:space:]]+")[[1L]]
  entries <- gsub("^'|'$", "", entries)
  entries[nzchar(entries)]
}

format_report <- function(inventory, collate, actual) {
  violations <- inventory[nzchar(inventory$violation), , drop = FALSE]
  duplicate_collate <- unique(collate[duplicated(collate)])
  missing_from_collate <- setdiff(actual, collate)
  missing_from_runtime <- setdiff(collate, actual)
  family_counts <- sort(table(inventory$family), decreasing = TRUE)

  lines <- c(
    "# Runtime Filename Audit",
    "",
    "The package uses its established `func_` prefix as the local equivalent of",
    "`golem::add_fct()`. Runtime files must use lower snake case and one of the",
    "documented Golem entrypoint, function, module, or utility forms.",
    "",
    "## Summary",
    "",
    sprintf("- Runtime `.R` files: `%d`", length(actual)),
    sprintf("- Conforming filenames: `%d`", sum(!nzchar(inventory$violation))),
    sprintf("- Filename violations: `%d`", nrow(violations)),
    sprintf("- `Collate` entries: `%d`", length(collate)),
    sprintf("- Duplicate `Collate` entries: `%d`", length(duplicate_collate)),
    sprintf("- Runtime files absent from `Collate`: `%d`", length(missing_from_collate)),
    sprintf("- `Collate` entries absent from runtime: `%d`", length(missing_from_runtime)),
    "",
    "## Families",
    ""
  )

  for (family in names(family_counts)) {
    lines <- c(lines, sprintf("- `%s`: `%d`", family, family_counts[[family]]))
  }

  lines <- c(lines, "", "## Violations", "")
  if (!nrow(violations) && !length(duplicate_collate) &&
      !length(missing_from_collate) && !length(missing_from_runtime)) {
    lines <- c(lines, "None.")
  } else {
    for (i in seq_len(nrow(violations))) {
      lines <- c(
        lines,
        sprintf("- `%s`: %s", violations$filename[[i]], violations$violation[[i]])
      )
    }
    for (filename in duplicate_collate) {
      lines <- c(lines, sprintf("- `%s`: duplicate `Collate` entry", filename))
    }
    for (filename in missing_from_collate) {
      lines <- c(lines, sprintf("- `%s`: absent from `Collate`", filename))
    }
    for (filename in missing_from_runtime) {
      lines <- c(lines, sprintf("- `%s`: `Collate` entry has no runtime file", filename))
    }
  }

  list(
    lines = lines,
    failure_count = nrow(violations) + length(duplicate_collate) +
      length(missing_from_collate) + length(missing_from_runtime)
  )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_root <- normalizePath(opts$repo_root, mustWork = TRUE)
  runtime_root <- file.path(repo_root, opts$root)
  if (!dir.exists(runtime_root)) {
    stop("Runtime root does not exist: ", opts$root, call. = FALSE)
  }

  actual <- sort(basename(list.files(runtime_root, pattern = "[.]R$", full.names = TRUE)))
  classified <- lapply(actual, classify_filename)
  inventory <- data.frame(
    filename = actual,
    family = vapply(classified, `[[`, character(1L), "family"),
    violation = vapply(classified, `[[`, character(1L), "violation"),
    stringsAsFactors = FALSE
  )
  collate <- read_collate(file.path(repo_root, "DESCRIPTION"))
  report <- format_report(inventory, collate, actual)

  if (is.null(opts$output)) {
    cat(paste(report$lines, collapse = "\n"), "\n")
  } else {
    output <- if (grepl("^/", opts$output)) {
      opts$output
    } else {
      file.path(repo_root, opts$output)
    }
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    writeLines(report$lines, output)
  }

  if (isTRUE(opts$check) && report$failure_count > 0L) {
    quit(status = 1L)
  }
}

main()
