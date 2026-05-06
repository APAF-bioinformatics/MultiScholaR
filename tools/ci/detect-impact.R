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
  changed_files = character(),
  changed_files_file = NULL,
  base_ref = NULL,
  head_ref = NULL,
  map = file.path("tools", "ci", "impact-map.json"),
  output = NULL,
  summary = NULL,
  github_output = NULL,
  artifact_dir = NULL,
  print_json = FALSE
)

impact_levels <- c("none", "targeted", "omic", "broad", "full")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

usage <- function() {
  cat(paste(
    "Usage: Rscript tools/ci/detect-impact.R [options]",
    "",
    "Options:",
    "  --changed-files <comma_or_newline_separated_paths>",
    "  --changed-files-file <path>",
    "  --base-ref <git_ref_or_sha>",
    "  --head-ref <git_ref_or_sha>",
    "  --map <impact-map.json>",
    "  --output <impact.json>",
    "  --summary <impact.md>",
    "  --github-output <GITHUB_OUTPUT file>",
    "  --artifact-dir <dir>",
    "  --print-json <true|false>",
    "  --help",
    sep = "\n"
  ))
}

parse_bool <- function(x) {
  isTRUE(tolower(as.character(x)) %in% c("1", "true", "yes", "y"))
}

split_paths <- function(x) {
  x <- unlist(strsplit(as.character(x), "[,\n\r]+"), use.names = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}

parse_args <- function(argv) {
  args <- default_args
  if (length(argv) == 0L) {
    return(args)
  }

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

    if (identical(key, "changed_files")) {
      args$changed_files <- c(args$changed_files, split_paths(value))
    } else if (identical(key, "print_json")) {
      args$print_json <- parse_bool(value)
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

compact_json <- function(x) {
  jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", pretty = FALSE)
}

read_changed_file <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(character())
  }
  if (!file.exists(path)) {
    stop(sprintf("changed-files file not found: %s", path), call. = FALSE)
  }
  split_paths(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

git_changed_files <- function(base_ref, head_ref) {
  if (is.null(base_ref) || is.null(head_ref) || !nzchar(base_ref) || !nzchar(head_ref)) {
    return(list(files = character(), error = NULL))
  }
  if (grepl("^0+$", base_ref)) {
    return(list(files = character(), error = "push base ref is all zeros"))
  }

  output <- tryCatch(
    system2("git", c("diff", "--name-status", base_ref, head_ref), stdout = TRUE, stderr = TRUE),
    error = function(err) {
      structure(conditionMessage(err), status = 1L)
    }
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    return(list(files = character(), error = paste(output, collapse = "\n")))
  }

  files <- character()
  for (line in output) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) >= 2L) {
      files <- c(files, parts[-1L])
    }
  }
  list(files = unique(files[nzchar(files)]), error = NULL)
}

as_list <- function(x) {
  if (is.null(x)) list() else x
}

as_chr <- function(x) {
  as.character(unlist(x, use.names = FALSE))
}

empty_payload <- function(args, errors = character()) {
  list(
    schema_version = "1.0.0",
    resolver_version = "1.0.0",
    generated_at = utc_now(),
    map_path = args$map,
    map_schema_version = NA_character_,
    impact_level = "none",
    changed_files = as.list(character()),
    matched_rules = list(),
    no_op_files = as.list(character()),
    unmatched_files = as.list(character()),
    suppressed_no_op_files = as.list(character()),
    module_matrix = list(include = list()),
    browser_matrix = list(include = list()),
    e2e_matrix = list(include = list()),
    module_count = 0L,
    browser_count = 0L,
    e2e_count = 0L,
    run_cross_omic = FALSE,
    run_report_export = FALSE,
    impact_summary = "No changed files selected CI.",
    errors = as.list(errors)
  )
}

broad_fallback_rules <- function(reason) {
  list(list(
    id = "builtin-broad-fallback",
    owner = "qa",
    impact_level = "broad",
    reason = reason,
    module_ci = list(
      list(name = "foundation", omic = "all", module = "foundation", runtime = "unit-contract"),
      list(name = "cross-module-integrity", omic = "all", module = "cross_module_integrity", runtime = "unit-contract"),
      list(name = "proteomics", omic = "proteomics", module = "all", runtime = "unit-contract"),
      list(name = "metabolomics", omic = "metabolomics", module = "all", runtime = "unit-contract"),
      list(name = "lipidomics", omic = "lipidomics", module = "all", runtime = "unit-contract")
    ),
    e2e = list(
      list(name = "proteomics-dia", lane = "prot_dia", filter = "^e2e-proteomics-dia$"),
      list(name = "metabolomics-lc-gc", lane = "metab_lc,metab_gc,metab_combined", filter = "^e2e-metabolomics-lc-gc$"),
      list(name = "lipidomics-canonical", lane = "lipid_canonical", filter = "^e2e-lipidomics-canonical$")
    ),
    browser = list(list(name = "harness", filter = "^(module-ci-harness|e2e-browser-harness)$")),
    run_cross_omic = TRUE,
    run_report_export = TRUE
  ))
}

validate_rule <- function(rule, index) {
  required <- c("id", "paths", "impact_level", "reason")
  missing <- setdiff(required, names(rule))
  if (length(missing) > 0L) {
    stop(sprintf("Rule %d missing fields: %s", index, paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (!rule$impact_level %in% impact_levels) {
    stop(sprintf("Rule %s has unsupported impact_level: %s", rule$id, rule$impact_level), call. = FALSE)
  }
  for (pattern in as_chr(rule$paths)) {
    tryCatch(grepl(pattern, ""), error = function(err) {
      stop(sprintf("Rule %s has invalid regex '%s': %s", rule$id, pattern, conditionMessage(err)), call. = FALSE)
    })
  }
}

read_impact_map <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("impact map not found: %s", path), call. = FALSE)
  }
  impact_map <- jsonlite::read_json(path, simplifyVector = FALSE)
  if (is.null(impact_map$schema_version) || is.null(impact_map$rules) || length(impact_map$rules) == 0L) {
    stop("impact map must contain schema_version and non-empty rules", call. = FALSE)
  }
  for (idx in seq_along(impact_map$rules)) {
    validate_rule(impact_map$rules[[idx]], idx)
  }
  impact_map
}

rule_matches_path <- function(rule, path) {
  any(vapply(as_chr(rule$paths), function(pattern) grepl(pattern, path, perl = TRUE), logical(1)))
}

matrix_key <- function(row, fields) {
  paste(vapply(fields, function(field) as.character(row[[field]] %||% ""), character(1)), collapse = "\r")
}

dedupe_rows <- function(rows, fields) {
  if (length(rows) == 0L) {
    return(list())
  }
  keys <- vapply(rows, matrix_key, character(1), fields = fields)
  rows[!duplicated(keys)]
}

add_artifact_dirs <- function(rows, prefix) {
  if (length(rows) == 0L) {
    return(rows)
  }
  lapply(seq_along(rows), function(idx) {
    row <- rows[[idx]]
    safe_name <- gsub("[^A-Za-z0-9._-]+", "-", row$name %||% paste0("row-", idx))
    row$artifact_dir <- file.path(prefix, safe_name)
    row
  })
}

level_max <- function(levels) {
  levels <- levels[levels %in% impact_levels]
  if (length(levels) == 0L) {
    return("none")
  }
  impact_levels[[max(match(levels, impact_levels))]]
}

resolve_impact <- function(changed_files, impact_map, args, seed_errors = character()) {
  changed_files <- sort(unique(gsub("^\\./", "", changed_files)))
  payload <- empty_payload(args, errors = seed_errors)
  payload$changed_files <- as.list(changed_files)
  payload$map_schema_version <- impact_map$schema_version

  if (length(changed_files) == 0L && length(seed_errors) == 0L) {
    return(payload)
  }

  matched <- list()
  matched_rule_ids <- character()
  no_op_files <- character()
  unmatched_files <- character()
  module_rows <- list()
  browser_rows <- list()
  e2e_rows <- list()
  levels <- character()
  run_cross_omic <- FALSE
  run_report_export <- FALSE

  rules <- impact_map$rules
  if (length(seed_errors) > 0L) {
    rules <- c(broad_fallback_rules(paste(seed_errors, collapse = "; ")), rules)
  }

  for (path in changed_files) {
    path_rules <- rules[vapply(rules, rule_matches_path, logical(1), path = path)]
    specific_rules <- path_rules[!grepl("(catchall|unknown-r-source-broad)$", vapply(path_rules, `[[`, character(1), "id"))]
    if (length(specific_rules) > 0L) {
      path_rules <- specific_rules
    }
    if (length(path_rules) == 0L) {
      unmatched_files <- c(unmatched_files, path)
      path_rules <- broad_fallback_rules(sprintf("Unmatched risky path: %s", path))
    }

    path_non_noop <- FALSE
    for (rule in path_rules) {
      rule_id <- rule$id
      matched_rule_ids <- c(matched_rule_ids, rule_id)
      if (!isTRUE(rule$no_op)) {
        path_non_noop <- TRUE
      }
      levels <- c(levels, rule$impact_level)
      module_rows <- c(module_rows, as_list(rule$module_ci))
      browser_rows <- c(browser_rows, as_list(rule$browser))
      e2e_rows <- c(e2e_rows, as_list(rule$e2e))
      run_cross_omic <- isTRUE(run_cross_omic || isTRUE(rule$run_cross_omic))
      run_report_export <- isTRUE(run_report_export || isTRUE(rule$run_report_export))
      matched <- c(matched, list(list(
        id = rule_id,
        owner = rule$owner %||% NA_character_,
        path = path,
        impact_level = rule$impact_level,
        no_op = isTRUE(rule$no_op),
        reason = rule$reason
      )))
    }

    if (!path_non_noop) {
      no_op_files <- c(no_op_files, path)
    }
  }

  module_rows <- dedupe_rows(module_rows, c("omic", "module", "runtime"))
  browser_rows <- dedupe_rows(browser_rows, c("filter"))
  e2e_rows <- dedupe_rows(e2e_rows, c("filter"))

  module_rows <- add_artifact_dirs(module_rows, file.path("tests", "testthat", "_module_ci_artifacts", "impact", "module"))
  browser_rows <- add_artifact_dirs(browser_rows, file.path("tests", "testthat", "_module_ci_artifacts", "impact", "browser"))
  e2e_rows <- add_artifact_dirs(e2e_rows, file.path("tests", "testthat", "_e2e_artifacts", "impact"))

  impact_level <- level_max(levels)
  if (length(seed_errors) > 0L) {
    impact_level <- level_max(c(impact_level, "broad"))
  }

  payload$impact_level <- impact_level
  payload$matched_rules <- matched
  payload$no_op_files <- as.list(sort(unique(no_op_files)))
  payload$suppressed_no_op_files <- as.list(sort(unique(intersect(no_op_files, changed_files))))
  payload$unmatched_files <- as.list(sort(unique(unmatched_files)))
  payload$module_matrix <- list(include = module_rows)
  payload$browser_matrix <- list(include = browser_rows)
  payload$e2e_matrix <- list(include = e2e_rows)
  payload$module_count <- length(module_rows)
  payload$browser_count <- length(browser_rows)
  payload$e2e_count <- length(e2e_rows)
  payload$run_cross_omic <- run_cross_omic
  payload$run_report_export <- run_report_export
  payload$impact_summary <- sprintf(
    "impact=%s changed=%d modules=%d browser=%d e2e=%d rules=%s",
    impact_level,
    length(changed_files),
    length(module_rows),
    length(browser_rows),
    length(e2e_rows),
    paste(sort(unique(matched_rule_ids)), collapse = ",")
  )
  payload
}

write_summary <- function(path, payload) {
  if (is.null(path) || !nzchar(path)) {
    return(invisible(NULL))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  changed_arg <- paste(unlist(payload$changed_files), collapse = ",")
  if (!nzchar(changed_arg)) {
    changed_arg <- "<none>"
  }
  module_commands <- if (payload$module_count == 0L) {
    "- none"
  } else {
    vapply(payload$module_matrix$include, function(row) {
      sprintf(
        "- `Rscript tools/ci/run-module-ci.R --omic %s --module %s --runtime %s --reporter summary --artifact-dir %s`",
        row$omic,
        row$module,
        row$runtime,
        row$artifact_dir
      )
    }, character(1))
  }
  e2e_commands <- if (payload$e2e_count == 0L) {
    "- none"
  } else {
    vapply(payload$e2e_matrix$include, function(row) {
      sprintf(
        "- `Rscript tools/ci/run-e2e-ci.R --lane %s --filter %s --reporter summary --artifact-dir %s`",
        row$lane,
        row$filter,
        row$artifact_dir
      )
    }, character(1))
  }
  lines <- c(
    "# CI Impact Routing Summary",
    "",
    sprintf("- Impact level: `%s`", payload$impact_level),
    sprintf("- Changed files: `%d`", length(payload$changed_files)),
    sprintf("- Module rows: `%d`", payload$module_count),
    sprintf("- Browser rows: `%d`", payload$browser_count),
    sprintf("- E2E rows: `%d`", payload$e2e_count),
    sprintf("- Cross-omic: `%s`", payload$run_cross_omic),
    sprintf("- Report/export: `%s`", payload$run_report_export),
    sprintf("- Resolver version: `%s`", payload$resolver_version),
    sprintf("- Map schema version: `%s`", payload$map_schema_version),
    "",
    "## Changed Files",
    if (length(payload$changed_files) == 0L) "- none" else paste0("- `", unlist(payload$changed_files), "`"),
    "",
    "## No-op Files",
    if (length(payload$no_op_files) == 0L) "- none" else paste0("- `", unlist(payload$no_op_files), "`"),
    "",
    "## Suppressed No-op Files",
    if (length(payload$suppressed_no_op_files) == 0L) "- none" else paste0("- `", unlist(payload$suppressed_no_op_files), "`"),
    "",
    "## Matched Rules",
    if (length(payload$matched_rules) == 0L) "- none" else vapply(payload$matched_rules, function(rule) {
      sprintf("- `%s` on `%s`: %s", rule$id, rule$path, rule$reason)
    }, character(1)),
    "",
    "## Selected Module CI",
    if (payload$module_count == 0L) "- none" else vapply(payload$module_matrix$include, function(row) {
      sprintf("- `%s`: `%s/%s/%s` -> `%s`", row$name, row$omic, row$module, row$runtime, row$artifact_dir)
    }, character(1)),
    "",
    "## Selected Browser",
    if (payload$browser_count == 0L) "- none" else vapply(payload$browser_matrix$include, function(row) {
      sprintf("- `%s`: `%s` -> `%s`", row$name, row$filter, row$artifact_dir)
    }, character(1)),
    "",
    "## Selected E2E",
    if (payload$e2e_count == 0L) "- none" else vapply(payload$e2e_matrix$include, function(row) {
      sprintf("- `%s`: lanes `%s`, filter `%s` -> `%s`", row$name, row$lane, row$filter, row$artifact_dir)
    }, character(1)),
    "",
    "## Unmatched Files",
    if (length(payload$unmatched_files) == 0L) "- none" else paste0("- `", unlist(payload$unmatched_files), "`"),
    "",
    "## Resolver Errors",
    if (length(payload$errors) == 0L) "- none" else paste0("- `", unlist(payload$errors), "`"),
    "",
    "## Local Reproduction",
    sprintf(
      "- `Rscript tools/ci/detect-impact.R --changed-files %s --output tests/testthat/_module_ci_artifacts/impact/local/impact.json --summary tests/testthat/_module_ci_artifacts/impact/local/impact.md`",
      changed_arg
    ),
    module_commands,
    e2e_commands,
    if (isTRUE(payload$run_cross_omic)) {
      "- `Rscript -e 'Sys.setenv(NOT_CRAN = \"true\"); pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE); testthat::test_dir(\"tests/testthat\", filter = \"^(e2e-cross-omic|workflow-server-shared|module-ci-sentinels)$\", reporter = \"summary\")'`"
    } else {
      "- cross-omic gate not selected"
    },
    if (isTRUE(payload$run_report_export)) {
      "- `Rscript -e 'Sys.setenv(NOT_CRAN = \"true\"); pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE); testthat::test_dir(\"tests/testthat\", filter = \"^(general-filemgmt-report-helpers-shared|general-filemgmt-export-contracts|e2e-proteomics-enrichment)$\", reporter = \"summary\")'`"
    } else {
      "- report/export gate not selected"
    }
  )
  writeLines(lines, path)
  invisible(path)
}

write_github_output <- function(path, payload) {
  if (is.null(path) || !nzchar(path)) {
    return(invisible(NULL))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- c(
    sprintf("impact_level=%s", payload$impact_level),
    sprintf("module_count=%s", payload$module_count),
    sprintf("browser_count=%s", payload$browser_count),
    sprintf("e2e_count=%s", payload$e2e_count),
    sprintf("run_cross_omic=%s", tolower(as.character(payload$run_cross_omic))),
    sprintf("run_report_export=%s", tolower(as.character(payload$run_report_export))),
    sprintf("module_matrix=%s", compact_json(payload$module_matrix)),
    sprintf("browser_matrix=%s", compact_json(payload$browser_matrix)),
    sprintf("e2e_matrix=%s", compact_json(payload$e2e_matrix)),
    sprintf("impact_summary=%s", payload$impact_summary)
  )
  write(lines, file = path, append = TRUE)
  invisible(path)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- parse_args(argv)
  if (!is.null(args$artifact_dir) && nzchar(args$artifact_dir)) {
    dir.create(args$artifact_dir, recursive = TRUE, showWarnings = FALSE)
    args$output <- args$output %||% file.path(args$artifact_dir, "impact.json")
    args$summary <- args$summary %||% file.path(args$artifact_dir, "impact.md")
  }

  changed <- c(args$changed_files, read_changed_file(args$changed_files_file))
  git_result <- git_changed_files(args$base_ref, args$head_ref)
  changed <- c(changed, git_result$files)
  seed_errors <- character()
  if (!is.null(git_result$error)) {
    seed_errors <- c(seed_errors, git_result$error)
  }
  if (length(seed_errors) > 0L && length(changed) == 0L) {
    changed <- "<resolver-error>"
  }

  payload <- tryCatch(
    {
      impact_map <- read_impact_map(args$map)
      resolve_impact(changed, impact_map, args, seed_errors = seed_errors)
    },
    error = function(err) {
      fallback_map <- list(schema_version = NA_character_, rules = broad_fallback_rules(conditionMessage(err)))
      resolve_impact(
        unique(c(changed, "<resolver-error>")),
        fallback_map,
        args,
        seed_errors = c(seed_errors, conditionMessage(err))
      )
    }
  )

  if (!is.null(args$output) && nzchar(args$output)) {
    dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(payload, args$output, auto_unbox = TRUE, pretty = TRUE, null = "null")
  }
  write_summary(args$summary, payload)
  write_github_output(args$github_output, payload)
  if (isTRUE(args$print_json) || is.null(args$output)) {
    cat(compact_json(payload), "\n")
  }
  if (payload$impact_level %in% impact_levels) {
    quit(status = 0L)
  }
  quit(status = 1L)
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
