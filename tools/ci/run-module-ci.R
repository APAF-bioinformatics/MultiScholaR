#!/usr/bin/env Rscript

default_args <- list(
  omic = NULL,
  module = NULL,
  runtime = "unit-contract",
  scenario = NULL,
  reporter = "summary",
  artifact_dir = file.path("tests", "testthat", "_module_ci_artifacts")
)

usage <- function() {
  cat(
    paste(
      "Usage: Rscript tools/ci/run-module-ci.R [options]",
      "",
      "Options:",
      "  --omic <all|proteomics|metabolomics|lipidomics>",
      "  --module <foundation|fixtures|import|design|qc|qc_peptide|qc_protein|normalization|differential_abundance|enrichment|summary_report|browser|e2e>",
      "  --runtime <unit-contract|module-browser|module-artifact|workflow-e2e|release-full>",
      "  --scenario <scenario_id_or_pack_id>",
      "  --reporter <testthat_reporter>",
      "  --artifact-dir <path>",
      "  --help",
      sep = "\n"
    )
  )
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

read_json_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  jsonlite::read_json(path, simplifyVector = FALSE)
}

as_chr <- function(x) {
  as.character(unlist(x, use.names = FALSE))
}

utc_now <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

matches_filter <- function(value, filter_value, include_all_value = FALSE) {
  if (is.null(filter_value) || identical(filter_value, "all")) {
    return(TRUE)
  }
  if (isTRUE(include_all_value) && identical(value, "all")) {
    return(TRUE)
  }
  identical(value, filter_value)
}

manifest_scenarios <- function(manifest_path) {
  manifest <- read_json_if_exists(manifest_path)
  if (is.null(manifest) || is.null(manifest$scenarios)) {
    return(list())
  }

  lapply(manifest$scenarios, function(scenario) {
    list(
      scenario_id = scenario$scenario_id,
      ticket_id = scenario$ticket_id,
      item_id = scenario$item_id,
      omic = scenario$omic,
      module = scenario$module_family,
      runtime = scenario$runtime,
      ci_lane = scenario$ci_lane,
      fixture = scenario$fixture$path,
      artifacts = lapply(scenario$expected_artifacts, function(artifact) artifact$path),
      source = "module_ci_manifest"
    )
  })
}

fixture_pack_scenarios <- function(pack_path) {
  catalog <- read_json_if_exists(pack_path)
  if (is.null(catalog) || is.null(catalog$packs)) {
    return(list())
  }

  lapply(catalog$packs, function(pack) {
    list(
      scenario_id = pack$pack_id,
      ticket_id = pack$ticket_id,
      item_id = paste(as_chr(pack$item_ids), collapse = ","),
      omic = pack$omic,
      module = "fixtures",
      runtime = "unit-contract",
      ci_lane = pack$ci_lane,
      fixture = pack$data_path,
      artifacts = as.list(c(pack$data_path, pack$design_path, pack$oracle_path)),
      source = "module_ci_fixture_pack",
      fixture_class = pack$fixture_class,
      import_format = pack$import_format,
      push_safe = pack$push_safe
    )
  })
}

select_scenarios <- function(scenarios, args) {
  Filter(function(scenario) {
    matches_filter(scenario$scenario_id, args$scenario) &&
      matches_filter(scenario$omic, args$omic, include_all_value = TRUE) &&
      matches_filter(scenario$module, args$module) &&
      matches_filter(scenario$runtime, args$runtime)
  }, scenarios)
}

tests_for_selection <- function(scenarios) {
  tests <- character()
  sources <- unique(vapply(scenarios, `[[`, character(1), "source"))
  modules <- unique(vapply(scenarios, `[[`, character(1), "module"))

  if ("module_ci_manifest" %in% sources) {
    tests <- c(tests, "module-ci-manifest")
  }
  if ("foundation" %in% modules) {
    tests <- c(tests, "module-ci-harness")
  }
  if ("import" %in% modules) {
    tests <- c(tests, "module-ci-oracles")
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "import") && identical(scenario$omic, "proteomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-prot-import")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "import") && identical(scenario$omic, "metabolomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-metab-import")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "import") && identical(scenario$omic, "lipidomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-lipid-import")
    }
  }
  if ("design" %in% modules) {
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "design") && identical(scenario$omic, "proteomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-prot-design")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "design") && identical(scenario$omic, "metabolomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-metab-design")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "design") && identical(scenario$omic, "lipidomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-lipid-design")
    }
  }
  if ("qc_peptide" %in% modules) {
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "qc_peptide") && identical(scenario$omic, "proteomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-prot-peptide-qc")
    }
  }
  if ("qc_protein" %in% modules) {
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "qc_protein") && identical(scenario$omic, "proteomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-prot-protein-qc")
    }
  }
  if ("qc" %in% modules) {
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "qc") && identical(scenario$omic, "metabolomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-metab-qc")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "qc") && identical(scenario$omic, "lipidomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-lipid-qc")
    }
  }
  if ("normalization" %in% modules) {
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "normalization") && identical(scenario$omic, "proteomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-prot-norm")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "normalization") && identical(scenario$omic, "metabolomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-metab-norm")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "normalization") && identical(scenario$omic, "lipidomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-lipid-norm")
    }
  }
  if ("differential_abundance" %in% modules) {
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "differential_abundance") && identical(scenario$omic, "proteomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-prot-da")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "differential_abundance") && identical(scenario$omic, "metabolomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-metab-da")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "differential_abundance") && identical(scenario$omic, "lipidomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-lipid-da")
    }
  }
  if ("enrichment" %in% modules) {
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "enrichment") && identical(scenario$omic, "proteomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-prot-enrich")
    }
  }
  if ("summary_report" %in% modules) {
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "summary_report") && identical(scenario$omic, "proteomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-prot-summary")
    }
    if (any(vapply(scenarios, function(scenario) {
      identical(scenario$module, "summary_report") && identical(scenario$omic, "metabolomics")
    }, logical(1)))) {
      tests <- c(tests, "module-ci-metab-summary")
    }
  }
  if ("fixtures" %in% modules) {
    tests <- c(tests, "module-ci-fixture-integrity", "module-ci-oracles")
  }

  unique(tests)
}

write_run_manifest <- function(path, payload) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, path = path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$artifact_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(MULTISCHOLAR_MODULE_CI_ARTIFACT_DIR = args$artifact_dir)

all_scenarios <- c(
  manifest_scenarios(file.path("tests", "testdata", "module_ci", "manifest.json")),
  fixture_pack_scenarios(file.path("tests", "testdata", "module_ci", "fixture_packs.json"))
)
selected <- select_scenarios(all_scenarios, args)
if (length(selected) == 0L) {
  stop("No module-CI scenarios matched the requested filters", call. = FALSE)
}

test_names <- tests_for_selection(selected)
if (length(test_names) == 0L) {
  stop("No module-CI test files mapped to the selected scenarios", call. = FALSE)
}

test_filter <- sprintf("^(%s)$", paste(test_names, collapse = "|"))
run_started_at <- utc_now()
test_error <- NULL

tryCatch(
  {
    Sys.setenv(NOT_CRAN = "true")
    pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
    testthat::test_dir(
      "tests/testthat",
      filter = test_filter,
      reporter = args$reporter
    )
  },
  error = function(err) {
    test_error <<- conditionMessage(err)
  }
)

result <- if (is.null(test_error)) "passed" else "failed"
selected <- lapply(selected, function(scenario) {
  scenario$result <- result
  scenario$failure_reason <- test_error
  scenario$run_artifacts <- as.list(c(
    file.path(args$artifact_dir, "module-ci-run-manifest.json"),
    file.path(args$artifact_dir, "module-ci-scorecard.json"),
    file.path(args$artifact_dir, "module-ci-scorecard.md")
  ))
  scenario
})

payload <- list(
  schema_version = "1.0.0",
  generated_at = utc_now(),
  started_at = run_started_at,
  finished_at = utc_now(),
  invocation = args,
  result = result,
  failure_reason = test_error,
  test_filter = test_filter,
  test_names = as.list(test_names),
  artifact_dir = args$artifact_dir,
  scenarios = selected,
  artifact_retention = list(
    scorecards = as.list(c("module-ci-run-manifest.json", "module-ci-scorecard.json", "module-ci-scorecard.md")),
    module_ci_artifacts = "tests/testthat/_module_ci_artifacts/**",
    e2e_artifacts = "tests/testthat/_e2e_artifacts/**",
    allowed_file_types = as.list(c(
      "screenshots",
      "browser_logs",
      "dom_snapshots",
      "state_digests",
      "session_exports",
      "reports",
      "publication_tables",
      "scorecards"
    ))
  )
)

write_run_manifest(file.path(args$artifact_dir, "module-ci-run-manifest.json"), payload)
message(sprintf(
  "module-CI result=%s scenarios=%d tests=%s artifact_dir=%s",
  result,
  length(selected),
  paste(test_names, collapse = ","),
  args$artifact_dir
))

if (!is.null(test_error)) {
  quit(status = 1L)
}
