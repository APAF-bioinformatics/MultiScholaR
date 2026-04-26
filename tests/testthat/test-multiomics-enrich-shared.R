# fidelity-coverage-compare: shared
library(testthat)

printStringDbFunctionalEnrichmentBarGraph <- get(
  "printStringDbFunctionalEnrichmentBarGraph",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

localNamespaceBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localNamespaceBinding(
      env = env,
      name = name,
      value = bindings[[name]],
      .local_envir = .local_envir
    )
  }
}

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

localFakeCheckmate <- function(env = parent.frame()) {
  old_lib_paths <- .libPaths()
  fake_lib <- file.path(tempdir(), "multischolar-fake-checkmate-lib")
  installed_pkg <- file.path(fake_lib, "checkmate")
  pkg_dir <- file.path(tempdir(), "multischolar-fake-checkmate-src")

  if ("checkmate" %in% loadedNamespaces()) {
    try(unloadNamespace("checkmate"), silent = TRUE)
  }

  if (!dir.exists(installed_pkg)) {
    unlink(pkg_dir, recursive = TRUE)
    dir.create(file.path(pkg_dir, "R"), recursive = TRUE, showWarnings = FALSE)
    dir.create(fake_lib, recursive = TRUE, showWarnings = FALSE)

    writeLines(
      c(
        "Package: checkmate",
        "Version: 0.0.1",
        "Title: Fake checkmate test double",
        "Description: Minimal namespace used by MultiScholaR coverage tests.",
        "License: MIT",
        "Encoding: UTF-8"
      ),
      file.path(pkg_dir, "DESCRIPTION")
    )
    writeLines(
      "export(assertDataFrame,assertString,assertChoice,assertNumber,assert,checkCharacter,checkFactor,assertNumeric,assertList)",
      file.path(pkg_dir, "NAMESPACE")
    )
    writeLines(
      c(
        "assertDataFrame <- function(x, min.rows = 0, .var.name = 'x') {",
        "  if (!is.data.frame(x) || nrow(x) < min.rows) {",
        "    stop(sprintf('%s must be a data frame with at least %d row(s).', .var.name, min.rows), call. = FALSE)",
        "  }",
        "  invisible(TRUE)",
        "}",
        "",
        "assertString <- function(x, min.chars = 0, .var.name = 'x') {",
        "  if (!is.character(x) || length(x) != 1 || is.na(x) || nchar(x) < min.chars) {",
        "    stop(sprintf('%s must be a non-missing scalar string.', .var.name), call. = FALSE)",
        "  }",
        "  invisible(TRUE)",
        "}",
        "",
        "assertChoice <- function(x, choices, .var.name = 'x') {",
        "  if (length(x) != 1 || is.na(x) || !(x %in% choices)) {",
        "    stop(sprintf('%s must be one of the allowed choices.', .var.name), call. = FALSE)",
        "  }",
        "  invisible(TRUE)",
        "}",
        "",
        "assertNumber <- function(x, lower = -Inf, upper = Inf, .var.name = 'x') {",
        "  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < lower || x > upper) {",
        "    stop(sprintf('%s must be a scalar number in range.', .var.name), call. = FALSE)",
        "  }",
        "  invisible(TRUE)",
        "}",
        "",
        "checkCharacter <- function(x, any.missing = FALSE, min.len = 0, ...) {",
        "  is.character(x) && length(x) >= min.len && (any.missing || !any(is.na(x)))",
        "}",
        "",
        "checkFactor <- function(x, any.missing = FALSE, min.len = 0, ...) {",
        "  is.factor(x) && length(x) >= min.len && (any.missing || !any(is.na(x)))",
        "}",
        "",
        "assert <- function(..., .var.name = 'x') {",
        "  checks <- list(...)",
        "  ok <- any(vapply(checks, function(item) identical(item, TRUE), logical(1)))",
        "  if (!ok) {",
        "    stop(sprintf('%s failed validation.', .var.name), call. = FALSE)",
        "  }",
        "  invisible(TRUE)",
        "}",
        "",
        "assertNumeric <- function(x, any.missing = FALSE, min.len = 0, .var.name = 'x') {",
        "  if (!is.numeric(x) || length(x) < min.len || (!any.missing && any(is.na(x)))) {",
        "    stop(sprintf('%s must be numeric.', .var.name), call. = FALSE)",
        "  }",
        "  invisible(TRUE)",
        "}",
        "",
        "assertList <- function(x, min.len = 0, .var.name = 'x') {",
        "  if (!is.list(x) || length(x) < min.len) {",
        "    stop(sprintf('%s must be a list.', .var.name), call. = FALSE)",
        "  }",
        "  invisible(TRUE)",
        "}"
      ),
      file.path(pkg_dir, "R", "checkmate.R")
    )
    utils::install.packages(pkg_dir, lib = fake_lib, repos = NULL, type = "source", quiet = TRUE)
  }

  .libPaths(c(fake_lib, old_lib_paths))

  withr::defer({
    if ("checkmate" %in% loadedNamespaces()) {
      try(unloadNamespace("checkmate"), silent = TRUE)
    }
    .libPaths(old_lib_paths)
  }, envir = env)

  invisible(fake_lib)
}

localFakeBiocEnrichmentStack <- function(env = parent.frame()) {
  old_lib_paths <- .libPaths()
  old_options <- options(
    multischolar.fake_clusterProfiler.enricher = NULL,
    multischolar.fake_clusterProfiler.GSEA = NULL,
    multischolar.fake_clusterProfiler.enricher_result = NULL,
    multischolar.fake_clusterProfiler.gsea_result = NULL,
    multischolar.fake_keggrest.keggList = NULL
  )

  fake_lib <- file.path(tempdir(), "multischolar-fake-bioc-enrichment-lib")
  dir.create(fake_lib, recursive = TRUE, showWarnings = FALSE)

  for (pkg in c("clusterProfiler", "KEGGREST", "BiocParallel")) {
    if (pkg %in% loadedNamespaces()) {
      try(unloadNamespace(pkg), silent = TRUE)
    }
  }

  install_fake_pkg <- function(pkg_name, description_lines, namespace_lines, r_lines) {
    installed_pkg <- file.path(fake_lib, pkg_name)
    pkg_dir <- file.path(tempdir(), paste0("multischolar-fake-", tolower(pkg_name), "-src"))

    if (!dir.exists(installed_pkg)) {
      unlink(pkg_dir, recursive = TRUE)
      dir.create(file.path(pkg_dir, "R"), recursive = TRUE, showWarnings = FALSE)

      writeLines(description_lines, file.path(pkg_dir, "DESCRIPTION"))
      writeLines(namespace_lines, file.path(pkg_dir, "NAMESPACE"))
      writeLines(r_lines, file.path(pkg_dir, "R", paste0(pkg_name, ".R")))

      utils::install.packages(pkg_dir, lib = fake_lib, repos = NULL, type = "source", quiet = TRUE)
    }
  }

  install_fake_pkg(
    "clusterProfiler",
    c(
      "Package: clusterProfiler",
      "Version: 0.0.1",
      "Title: Fake clusterProfiler test double",
      "Description: Minimal namespace used by MultiScholaR enrichment coverage tests.",
      "License: MIT",
      "Encoding: UTF-8",
      "Depends: methods"
    ),
    c(
      "export(enricher,GSEA)",
      "exportClasses(enrichResult,gseaResult)"
    ),
    c(
      "methods::setClass('enrichResult', slots = c(result = 'data.frame'))",
      "methods::setClass('gseaResult', slots = c(result = 'data.frame'))",
      "",
      ".call_handler <- function(handler, args) {",
      "  formal_names <- names(formals(handler))",
      "  if (is.null(formal_names) || '...' %in% formal_names) {",
      "    return(do.call(handler, args))",
      "  }",
      "  do.call(handler, args[intersect(names(args), formal_names)])",
      "}",
      "",
      "enricher <- function(gene, universe = NULL, TERM2GENE = NULL, TERM2NAME = NULL, ...) {",
      "  handler <- getOption('multischolar.fake_clusterProfiler.enricher')",
      "  args <- list(gene = gene, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, ...)",
      "  if (is.function(handler)) {",
      "    return(.call_handler(handler, args))",
      "  }",
      "  result <- getOption('multischolar.fake_clusterProfiler.enricher_result')",
      "  if (is.null(result)) {",
      "    result <- data.frame()",
      "  }",
      "  methods::new('enrichResult', result = result)",
      "}",
      "",
      "GSEA <- function(geneList, TERM2GENE = NULL, TERM2NAME = NULL, ...) {",
      "  handler <- getOption('multischolar.fake_clusterProfiler.GSEA')",
      "  args <- list(geneList = geneList, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, ...)",
      "  if (is.function(handler)) {",
      "    return(.call_handler(handler, args))",
      "  }",
      "  result <- getOption('multischolar.fake_clusterProfiler.gsea_result')",
      "  if (is.null(result)) {",
      "    result <- data.frame()",
      "  }",
      "  methods::new('gseaResult', result = result)",
      "}"
    )
  )

  install_fake_pkg(
    "KEGGREST",
    c(
      "Package: KEGGREST",
      "Version: 0.0.1",
      "Title: Fake KEGGREST test double",
      "Description: Minimal namespace used by MultiScholaR KEGG coverage tests.",
      "License: MIT",
      "Encoding: UTF-8"
    ),
    "export(keggList)",
    c(
      "keggList <- function(database) {",
      "  handler <- getOption('multischolar.fake_keggrest.keggList')",
      "  if (is.function(handler)) {",
      "    return(handler(database))",
      "  }",
      "  stats::setNames(character(), character())",
      "}"
    )
  )

  install_fake_pkg(
    "BiocParallel",
    c(
      "Package: BiocParallel",
      "Version: 0.0.1",
      "Title: Fake BiocParallel test double",
      "Description: Minimal namespace used by MultiScholaR enrichment coverage tests.",
      "License: MIT",
      "Encoding: UTF-8"
    ),
    "export(SerialParam)",
    c(
      "SerialParam <- function(...) {",
      "  structure(list(...), class = 'SerialParam')",
      "}"
    )
  )

  .libPaths(c(fake_lib, old_lib_paths))

  withr::defer({
    options(old_options)
    for (pkg in c('clusterProfiler', 'KEGGREST', 'BiocParallel')) {
      if (pkg %in% loadedNamespaces()) {
        try(unloadNamespace(pkg), silent = TRUE)
      }
    }
    .libPaths(old_lib_paths)
  }, envir = env)

  invisible(fake_lib)
}

test_that("STRING submission helper preserves payload construction and failure handling", {
  localFakeCheckmate()

  request_log <- new.env(parent = emptyenv())
  request_log$bodies <- list()
  request_log$post_mode <- "success"
  request_log$json_mode <- "success"

  local_mocked_bindings(
    p_load = function(...) TRUE,
    .package = "pacman"
  )
  local_mocked_bindings(
    POST = function(url, body, encode) {
      request_log$bodies[[length(request_log$bodies) + 1L]] <- list(
        url = url,
        body = body,
        encode = encode
      )
      if (identical(request_log$post_mode, "error")) {
        stop("network down")
      }
      if (identical(request_log$post_mode, "http_error")) {
        return(list(status = 500L, text = "{\"status\":\"error\",\"message\":\"bad\"}"))
      }
      list(status = 200L, text = "[{\"job_id\":\"job-123\",\"status\":\"submitted\"}]")
    },
    content = function(x, as = "parsed", encoding = NULL) {
      if (identical(as, "text")) {
        return(x$text)
      }
      stop("unexpected content mode")
    },
    http_error = function(x) isTRUE(x$status >= 400L),
    status_code = function(x) x$status,
    .package = "httr"
  )
  local_mocked_bindings(
    fromJSON = function(text, ...) {
      if (identical(request_log$json_mode, "error")) {
        stop("invalid json")
      }
      list(list(job_id = "job-123", status = "submitted"))
    },
    .package = "jsonlite"
  )

  input_tbl <- data.frame(
    protein = c("P1", "P2", NA, "P4"),
    score = c(3.4, NA, 5.0, -2.1),
    stringsAsFactors = FALSE
  )

  result <- submitStringDBEnrichment(
    input_data_frame = input_tbl,
    identifier_column_name = "protein",
    value_column_name = "score",
    caller_identity = "demo-run",
    api_key = "token-1",
    species = 9606
  )

  expect_identical(result$job_id, "job-123")
  expect_identical(result$api_key, "token-1")
  expect_identical(request_log$bodies[[1]]$encode, "form")
  expect_identical(request_log$bodies[[1]]$body$species, "9606")
  expect_identical(
    request_log$bodies[[1]]$body$identifiers,
    paste(c("P1\t3.4", "P4\t-2.1"), collapse = "\n")
  )

  request_log$post_mode <- "http_error"
  http_error_result <- submitStringDBEnrichment(
    input_data_frame = input_tbl,
    identifier_column_name = "protein",
    value_column_name = "score",
    caller_identity = "demo-run",
    api_key = "token-2"
  )
  expect_null(http_error_result$job_id)
  expect_identical(http_error_result$submission_response$status, "error")

  request_log$post_mode <- "error"
  post_error_result <- submitStringDBEnrichment(
    input_data_frame = input_tbl,
    identifier_column_name = "protein",
    value_column_name = "score",
    caller_identity = "demo-run",
    api_key = "token-3"
  )
  expect_identical(post_error_result$job_id, "job-123")
  expect_identical(post_error_result$submission_response$status, "submitted")

  request_log$post_mode <- "success"
  request_log$json_mode <- "error"
  parse_error_result <- submitStringDBEnrichment(
    input_data_frame = input_tbl,
    identifier_column_name = "protein",
    value_column_name = "score",
    caller_identity = "demo-run",
    api_key = "token-4"
  )
  expect_null(parse_error_result$job_id)
  expect_null(parse_error_result$submission_response$message)
})

test_that("STRING download helpers preserve graph and table retrieval behavior", {
  localFakeCheckmate()

  state <- new.env(parent = emptyenv())
  state$mode <- "graph_ok"

  local_mocked_bindings(
    GET = function(url, query = NULL) {
      switch(
        state$mode,
        graph_ok = list(status = 200L, raw = charToRaw("PNGDATA"), text = ""),
        graph_http_error = list(status = 404L, raw = raw(), text = "not found"),
        graph_empty = list(status = 200L, raw = raw(), text = ""),
        table_ok = list(
          status = 200L,
          text = paste("category\ttermDescription", "KEGG\tPathway A", sep = "\n")
        ),
        table_http_error = list(status = 500L, text = "server error"),
        stop("request boom")
      )
    },
    content = function(x, as = "parsed", encoding = NULL) {
      if (identical(as, "raw")) {
        return(x$raw)
      }
      if (identical(as, "text")) {
        return(x$text)
      }
      stop("unexpected content mode")
    },
    http_error = function(x) isTRUE(x$status >= 400L),
    status_code = function(x) x$status,
    .package = "httr"
  )

  graph <- downloadStringDBGraph("https://string.test/graph")
  expect_equal(graph, charToRaw("PNGDATA"))

  state$mode <- "graph_http_error"
  expect_null(downloadStringDBGraph("https://string.test/graph"))

  state$mode <- "graph_empty"
  expect_null(downloadStringDBGraph("https://string.test/graph"))

  state$mode <- "table_ok"
  table_result <- downloadStringDBResultsFile("https://string.test/results")
  expect_identical(table_result$category[[1]], "KEGG")

  state$mode <- "table_http_error"
  expect_null(downloadStringDBResultsFile("https://string.test/results"))

  local_mocked_bindings(
    read_tsv = function(...) stop("bad tsv"),
    .package = "readr"
  )
  state$mode <- "table_ok"
  expect_null(downloadStringDBResultsFile("https://string.test/results"))
})

test_that("STRING retrieval helper polls until completion and reports terminal failures", {
  localFakeCheckmate()

  package_ns <- asNamespace("MultiScholaR")
  polling <- new.env(parent = emptyenv())
  polling$count <- 0L
  polling$mode <- "success"
  base_ns <- asNamespace("base")

  localNamespaceBinding(
    base_ns,
    "Sys.sleep",
    function(time) invisible(time)
  )
  localNamespaceBindings(
    package_ns,
    list(
      downloadStringDBResultsFile = function(download_url) {
        expect_identical(download_url, "https://string.test/download")
        data.frame(termDescription = "Pathway A", stringsAsFactors = FALSE)
      },
      downloadStringDBGraph = function(graph_url) {
        expect_identical(graph_url, "https://string.test/graph")
        charToRaw("PNG")
      }
    )
  )
  local_mocked_bindings(
    GET = function(url, query = NULL) {
      polling$count <- polling$count + 1L
      if (identical(polling$mode, "http_error")) {
        return(list(status = 500L, text = "{\"status\":\"error\",\"message\":\"fail\"}"))
      }
      list(status = 200L, text = sprintf("poll-%d", polling$count))
    },
    content = function(x, as = "parsed", encoding = NULL) x$text,
    http_error = function(x) isTRUE(x$status >= 400L),
    status_code = function(x) x$status,
    .package = "httr"
  )
  local_mocked_bindings(
    fromJSON = function(text, simplifyDataFrame = TRUE) {
      if (identical(polling$mode, "success")) {
        if (polling$count == 1L) {
          return(data.frame(
            status = "running",
            message = "queued",
            stringsAsFactors = FALSE
          ))
        }
        return(data.frame(
          status = "success",
          message = "done",
          page_url = "https://string.test/page",
          download_url = "https://string.test/download",
          graph_url = "https://string.test/graph",
          stringsAsFactors = FALSE
        ))
      }
      data.frame(status = "error", message = "failed", stringsAsFactors = FALSE)
    },
    .package = "jsonlite"
  )

  results <- retrieveStringDBEnrichmentResults(
    submission_info = list(job_id = "job-1", api_key = "token"),
    polling_interval_seconds = 1,
    max_polling_attempts = 3
  )

  expect_identical(results$page_url, "https://string.test/page")
  expect_identical(results$download_url, "https://string.test/download")
  expect_equal(results$graph_image_content, charToRaw("PNG"))
  expect_identical(results$enrichment_data$termDescription[[1]], "Pathway A")

  polling$count <- 0L
  polling$mode <- "api_error"
  expect_null(
    retrieveStringDBEnrichmentResults(
      submission_info = list(job_id = "job-1", api_key = "token"),
      polling_interval_seconds = 1,
      max_polling_attempts = 2
    )
  )

  polling$count <- 0L
  polling$mode <- "http_error"
  expect_null(
    retrieveStringDBEnrichmentResults(
      submission_info = list(job_id = "job-1", api_key = "token"),
      polling_interval_seconds = 1,
      max_polling_attempts = 1
    )
  )
})

test_that("species discovery helpers standardize API results and limit wrapper output", {
  species_fun <- makeFunctionWithOverrides(
    getStringDbSpecies,
    list(syms = rlang::syms)
  )

  local_mocked_bindings(
    GET = function(url, query = NULL) list(status = 200L, text = "species"),
    content = function(x, as = "parsed", encoding = NULL) x$text,
    http_error = function(x) FALSE,
    status_code = function(x) x$status,
    .package = "httr"
  )
  local_mocked_bindings(
    fromJSON = function(text, flatten = TRUE) {
      data.frame(
        taxon_id = c("10090", "9606", "559292"),
        species_name = c("Mus musculus", "Homo sapiens", "Saccharomyces cerevisiae"),
        short_name = c("mouse", "human", "yeast"),
        stringsAsFactors = FALSE
      )
    },
    .package = "jsonlite"
  )

  filtered <- species_fun(search_term = "sapiens")
  expect_true(all(c("species_id", "official_name", "compact_name") %in% names(filtered)))
  expect_true("Homo sapiens" %in% filtered$official_name)

  package_ns <- asNamespace("MultiScholaR")
  localNamespaceBinding(
    package_ns,
    "getStringDbSpecies",
    function(search_term = NULL, api_key = NULL) {
      tibble::tibble(
        species_id = c("1", "2", "3"),
        official_name = c("alpha species", "beta species", "gamma species")
      )
    }
  )
  truncated <- searchStringDbSpecies("species", show_top_n = 2)
  expect_equal(nrow(truncated), 2L)

  localNamespaceBinding(
    package_ns,
    "getStringDbSpecies",
    function(search_term = NULL, api_key = NULL) tibble::tibble()
  )
  expect_equal(nrow(searchStringDbSpecies("missing")), 0L)
  expect_error(searchStringDbSpecies(""), "Please provide a species name")
})

test_that("rank-based STRING enrichment helpers preserve scoring, file output, and S4 ranking flows", {
  package_ns <- asNamespace("MultiScholaR")
  recorded <- new.env(parent = emptyenv())
  recorded$submit <- list()
  recorded$mofa <- list()

  localNamespaceBindings(
    package_ns,
    list(
      submitStringDBEnrichment = function(input_data_frame, ...) {
        recorded$submit[[length(recorded$submit) + 1L]] <- input_data_frame
        list(job_id = "job-1", api_key = "token", submission_response = list(status = "submitted"))
      },
      retrieveStringDBEnrichmentResults = function(submission_info, ...) {
        list(
          enrichment_data = data.frame(
            termDescription = "Pathway A",
            enrichmentScore = 1.2,
            falseDiscoveryRate = 0.02,
            genesMapped = 2L,
            stringsAsFactors = FALSE
          ),
          page_url = "https://string.test/page",
          download_url = "https://string.test/download",
          graph_url = "https://string.test/graph",
          graph_image_content = charToRaw("PNG")
        )
      },
      runOneStringDbRankEnrichmentMofa = function(input_table,
                                                  identifier_column_name,
                                                  value_column_name,
                                                  result_label,
                                                  results_dir,
                                                  ...) {
        recorded$mofa <- list(
          input_table = input_table,
          identifier_column_name = identifier_column_name,
          value_column_name = value_column_name,
          result_label = result_label,
          results_dir = results_dir
        )
        data.frame(termDescription = "MOFA Pathway", stringsAsFactors = FALSE)
      }
    )
  )

  fixture_dir <- tempfile("string-rank-")
  dir.create(fixture_dir, recursive = TRUE)
  withr::defer(unlink(fixture_dir, recursive = TRUE, force = TRUE))

  ranked_input <- data.frame(
    Protein.Ids = c("P1:iso1", "P2:iso2", "P3"),
    log2FC = c(2, -1.5, 0.5),
    fdr_qvalue = c(0.01, 0.02, 0.5),
    stringsAsFactors = FALSE
  )

  one_result <- runOneStringDbRankEnrichment(
    input_table = ranked_input,
    result_label = "contrast_a",
    pathway_dir = fixture_dir,
    api_key = "token"
  )

  expect_equal(nrow(one_result), 1L)
  expect_equal(recorded$submit[[1]]$protein_id, c("P1", "P3", "P2"))
  expect_equal(
    round(recorded$submit[[1]]$score, 6),
    round(c(2, 0.30103, -1.69897), 6)
  )
  expect_true(file.exists(file.path(fixture_dir, "string_db", "contrast_a_string_enrichment_results.tab")))
  expect_true(file.exists(file.path(fixture_dir, "string_db", "contrast_a_string_enrichment_graph.png")))

  mofa_dir <- tempfile("string-mofa-")
  dir.create(mofa_dir, recursive = TRUE)
  withr::defer(unlink(mofa_dir, recursive = TRUE, force = TRUE))

  mofa_result <- suppressWarnings(
    runOneStringDbRankEnrichmentMofa(
      input_table = data.frame(protein_id = c("P1", "P2"), score = c(2, -1), stringsAsFactors = FALSE),
      result_label = "factor1",
      results_dir = mofa_dir,
      api_key = "token"
    )
  )
  expect_equal(nrow(mofa_result), 1L)

  da_obj <- methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrast = "unused"),
    da_data = list(
      contrast_a = data.frame(
        Protein.Ids = c("P1:iso1", "P2:iso2", "P3"),
        fdr_qvalue = c(0.01, 0.20, 0.03),
        log2FC = c(2.5, -1, -3),
        stringsAsFactors = FALSE
      ),
      contrast_b = data.frame(
        Protein.Ids = c("P4", "P5"),
        fdr_qvalue = c(0.04, 0.20),
        log2FC = c(1.5, -0.5),
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE)
  )

  da_result <- runStringDbEnrichmentFromDAResults(
    da_results_for_enrichment = da_obj,
    contrast_name = "contrast_a",
    ranking_method = "combined_score",
    filter_significant = TRUE,
    fdr_threshold = 0.05,
    results_dir = fixture_dir,
    api_key = "token"
  )
  expect_equal(nrow(da_result), 1L)
  expect_equal(recorded$mofa$identifier_column_name, "protein_id")
  expect_equal(recorded$mofa$value_column_name, "score")
  expect_equal(recorded$mofa$result_label, "contrast_a_combined_score")
  expect_equal(recorded$mofa$input_table$protein_id, c("P1", "P3"))

  localNamespaceBinding(
    package_ns,
    "runStringDbEnrichmentFromDAResults",
    function(da_results_for_enrichment, contrast_name, ranking_method, result_label, ...) {
      if (identical(contrast_name, "contrast_b")) {
        stop("contrast failed")
      }
      data.frame(
        contrast = contrast_name,
        ranking_method = ranking_method,
        result_label = result_label,
        stringsAsFactors = FALSE
      )
    }
  )

  multi_result <- runStringDbEnrichmentFromDAResultsMultiple(
    da_results_for_enrichment = da_obj,
    contrast_names = c("contrast_a", "contrast_b"),
    ranking_method = "none",
    api_key = "token"
  )
  expect_identical(multi_result$contrast_a$contrast[[1]], "contrast_a")
  expect_null(multi_result$contrast_b)
})

test_that("STRING contrast aggregation and plotting helpers preserve cache and filtered plot structure", {
  package_ns <- asNamespace("MultiScholaR")
  call_log <- new.env(parent = emptyenv())
  call_log$labels <- character()

  localNamespaceBinding(
    package_ns,
    "runOneStringDbRankEnrichment",
    function(input_table, result_label, pathway_dir, ...) {
      call_log$labels <- c(call_log$labels, result_label)
      data.frame(
        termID = c("GO:1", "GO:2"),
        termDescription = c("Cell cycle response term", "Metabolism term"),
        category = c("GO Process", "KEGG"),
        enrichmentScore = c(3.2, 2.4),
        falseDiscoveryRate = c(0.01, 0.02),
        direction = c("top", "bottom"),
        genesMapped = c(4L, 3L),
        stringsAsFactors = FALSE
      )
    }
  )

  fixture_dir <- tempfile("string-aggregate-")
  pathway_dir <- file.path(fixture_dir, "pathway")
  dir.create(pathway_dir, recursive = TRUE)
  withr::defer(unlink(fixture_dir, recursive = TRUE, force = TRUE))

  project_dirs <- list(
    proteomics_demo = list(pathway_dir = pathway_dir)
  )
  da_list <- list(
    `alpha=treated_control` = list(
      da_proteins_long = data.frame(
        Protein.Ids = c("P1", "P2"),
        log2FC = c(2, -1),
        fdr_qvalue = c(0.01, 0.02),
        stringsAsFactors = FALSE
      )
    ),
    `beta=treated_control` = list(
      de_proteins_long = data.frame(
        Protein.Ids = c("P3", "P4"),
        log2FC = c(1.5, -1.2),
        fdr_qvalue = c(0.03, 0.04),
        stringsAsFactors = FALSE
      )
    )
  )

  combined <- runStringDbEnrichmentAllContrasts(
    da_analysis_results_list = da_list,
    project_dirs = project_dirs,
    omic_type = "proteomics",
    experiment_label = "demo",
    force_refresh = TRUE
  )

  expect_equal(sort(unique(combined$comparison)), c("alpha=treated (M-C)", "beta=treated (M-C)"))
  expect_equal(call_log$labels, c("alpha", "beta"))
  expect_true(file.exists(file.path(pathway_dir, "string_db", "all_enrichment_results.rds")))

  localNamespaceBinding(
    package_ns,
    "runOneStringDbRankEnrichment",
    function(...) stop("cache should be used")
  )
  cached <- runStringDbEnrichmentAllContrasts(
    da_analysis_results_list = da_list,
    project_dirs = project_dirs,
    omic_type = "proteomics",
    experiment_label = "demo",
    force_refresh = FALSE
  )
  expect_equal(nrow(cached), nrow(combined))

  localNamespaceBinding(
    package_ns,
    "printStringDbFunctionalEnrichmentBarGraph",
    function(input_table, word_limit = 10) {
      ggplot2::ggplot(input_table, ggplot2::aes(enrichmentScore, termDescription)) +
        ggplot2::geom_point()
    }
  )

  plot_results <- plotStringDbEnrichmentResults(
    project_dirs = project_dirs,
    omic_type = "proteomics",
    experiment_label = "demo",
    enrichment_results = data.frame(
      comparison = c("A", "A", "A", "A"),
      category = c("GO Process", "KEGG", "InterPro", "Custom"),
      termID = c("GO:1", "K1", "IPR:1", "C1"),
      termDescription = c("Cell cycle response term", "Pathway signal term", "Domain feature term", "Custom bucket term"),
      enrichmentScore = c(4, 3, 2, 1),
      falseDiscoveryRate = c(0.01, 0.02, 0.03, 0.04),
      direction = c("top", "bottom", "both ends", "top"),
      genesMapped = c(4L, 3L, 2L, 1L),
      stringsAsFactors = FALSE
    ),
    top_n_terms = 3,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_true(is.list(plot_results$plots))
  expect_s3_class(plot_results$plots$go_terms, "ggplot")
  expect_s3_class(plot_results$plots$pathways, "ggplot")
  expect_s3_class(plot_results$plots$protein_domains, "ggplot")
  expect_s3_class(plot_results$plots$others, "ggplot")

  plot_fun <- makeFunctionWithOverrides(
    printStringDbFunctionalEnrichmentBarGraph,
    list(
      reorder_within = function(x, by, within) x,
      scale_y_reordered = function(...) ggplot2::scale_y_discrete()
    )
  )
  direct_plot <- plot_fun(
    data.frame(
      comparison = "A",
      category = "GO Process",
      termDescription = "A very long description for truncation checks",
      enrichmentScore = 3.5,
      falseDiscoveryRate = 0.01,
      direction = "top",
      genesMapped = 4L,
      stringsAsFactors = FALSE
    ),
    word_limit = 3
  )
  expect_s3_class(direct_plot, "ggplot")
})

test_that("metabolomics orchestration helpers preserve combined assay and mapped-name behavior", {
  package_ns <- asNamespace("MultiScholaR")
  results_dir <- tempfile("metab-enrichment-")
  data_dir <- tempfile("metab-data-")
  dir.create(results_dir, recursive = TRUE)
  dir.create(data_dir, recursive = TRUE)
  withr::defer({
    unlink(results_dir, recursive = TRUE, force = TRUE)
    unlink(data_dir, recursive = TRUE, force = TRUE)
  })

  metabolomics_obj <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      data.frame(
        metabolite = c("met_lc_a", "met_lc_b"),
        database_identifier = c("CHEBI:1", "CHEBI:2"),
        metabolite_identification = c("met_lc_a", "met_lc_b"),
        S1 = c(10, 20),
        stringsAsFactors = FALSE
      ),
      data.frame(
        metabolite = c("met_gc_a", "met_gc_b"),
        database_identifier = c("CHEBI:3", "ITSD:99"),
        metabolite_identification = c("met_gc_a", "ITSD control"),
        S1 = c(30, 40),
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(Sample_ID = "S1", group = "A", stringsAsFactors = FALSE)
  )

  weights <- data.frame(
    view = c("metabolome_lc", "metabolome_lc", "metabolome_gc"),
    factor = c("Factor1", "Factor1", "Factor1"),
    feature = c("met_lc_a_metabolome_lc", "met_lc_b_metabolome_lc", "met_gc_a_metabolome_gc"),
    value = c(2.5, -1.5, 1.2),
    stringsAsFactors = FALSE
  )
  mapping_table <- data.frame(
    KEGG = c("C00001", "C00003"),
    ChEBI = c("1", "3"),
    stringsAsFactors = FALSE
  )

  localNamespaceBindings(
    package_ns,
    list(
      getProjectPaths = function(omic_type, experiment_label) {
        list(
          data_dir = data_dir,
          integration_enrichment_plots_dir = results_dir
        )
      },
      runKeggEnrichment = function(ranked_list, mapping_table, project_dirs, omic_type, experiment_label, assay_name, ...) {
        data.frame(
          termDescription = paste("KEGG", assay_name),
          enrichmentScore = 2.1,
          falseDiscoveryRate = 0.01,
          genesMapped = 1L,
          mappedIDs = names(ranked_list)[1],
          comparison = assay_name,
          stringsAsFactors = FALSE
        )
      },
      runReactomeEnrichment = function(ranked_list, mapping_table, project_dirs, omic_type, experiment_label, assay_name, ...) {
        data.frame(
          termDescription = paste("Reactome", assay_name),
          enrichmentScore = 1.8,
          falseDiscoveryRate = 0.02,
          genesMapped = 1L,
          mappedIDs = names(ranked_list)[1],
          comparison = assay_name,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  assay_results <- runMetabolomicsEnrichmentAnalysis(
    weights = weights,
    metabolomics_obj = metabolomics_obj,
    mapping_table = mapping_table,
    project_dirs = list(),
    omic_type = "metabolomics",
    experiment_label = "demo",
    assay_name = "metabolome_lc"
  )
  expect_equal(sort(unique(assay_results$category)), c("KEGG", "Reactome"))
  expect_true(file.exists(file.path(results_dir, "metabolome_lc_enrichment_results.tab")))

  localNamespaceBinding(
    package_ns,
    "runMetabolomicsEnrichmentAnalysis",
    function(weights, metabolomics_obj, mapping_table, project_dirs, omic_type, experiment_label, assay_name, ...) {
      if (identical(assay_name, "metabolome_lc")) {
        return(data.frame(
          termDescription = "LC Pathway",
          enrichmentScore = 2.5,
          falseDiscoveryRate = 0.01,
          genesMapped = 2L,
          mappedIDs = "CHEBI:1/cpd:C00001",
          comparison = "metabolome_lc",
          category = "KEGG",
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        termDescription = "GC Pathway",
        enrichmentScore = 1.5,
        falseDiscoveryRate = 0.03,
        genesMapped = 1L,
        mappedIDs = "CHEBI:3",
        comparison = "metabolome_gc",
        category = "Reactome",
        stringsAsFactors = FALSE
      )
    }
  )

  combined_results <- runMetabolomicsPathwayEnrichment(
    weights = weights,
    metabolomics_obj = metabolomics_obj,
    mapping_table = mapping_table,
    project_dirs = list(),
    omic_type = "metabolomics",
    experiment_label = "demo"
  )

  expect_equal(sort(unique(combined_results$assay)), c("GC-MS", "LC-MS"))
  expect_equal(unique(combined_results$comparison), "Metabolome")
  expect_match(combined_results$mappedNames[[1]], "met_lc_a")
  expect_true(file.exists(file.path(results_dir, "combined_metabolomics_enrichment_results.tab")))
})

test_that("KEGG and Reactome helpers preserve early-return branches when inputs cannot reach enrichment", {
  package_ns <- asNamespace("MultiScholaR")
  data_dir <- tempfile("metab-pathway-data-")
  dir.create(data_dir, recursive = TRUE)
  withr::defer(unlink(data_dir, recursive = TRUE, force = TRUE))

  localNamespaceBinding(
    package_ns,
    "getProjectPaths",
    function(omic_type, experiment_label) list(data_dir = data_dir)
  )

  mapping_table <- data.frame(
    KEGG = c("C00001", "C00002"),
    ChEBI = c("7", "8"),
    stringsAsFactors = FALSE
  )

  kegg_empty <- runKeggEnrichment(
    ranked_list = c(`CHEBI:1` = 2, `CHEBI:2` = -1),
    mapping_table = mapping_table,
    project_dirs = list(),
    omic_type = "metabolomics",
    experiment_label = "demo",
    assay_name = "metabolome_lc"
  )
  expect_equal(nrow(kegg_empty), 0L)

  readr::write_tsv(
    tibble::tibble(
      chebi_id = c("CHEBI:7", "CHEBI:8"),
      reactome_id = c("R-HSA-1", "R-HSA-2"),
      reactome_term = c("Path A", "Path B"),
      organism = c("Homo sapiens", "Homo sapiens")
    ),
    file.path(data_dir, "chebi_to_reactome.tab")
  )

  reactome_empty <- runReactomeEnrichment(
    ranked_list = c(`CHEBI:1` = 2, `CHEBI:2` = -1),
    mapping_table = mapping_table,
    project_dirs = list(),
    omic_type = "metabolomics",
    experiment_label = "demo",
    assay_name = "metabolome_lc",
    reactome_organism = "Homo sapiens"
  )
  expect_equal(nrow(reactome_empty), 0L)
})

test_that("KEGG and Reactome helpers preserve mapped enrichment formatting on success and fallback branches", {
  localFakeBiocEnrichmentStack()
  package_ns <- asNamespace("MultiScholaR")
  data_dir <- tempfile("metab-pathway-rich-")
  dir.create(data_dir, recursive = TRUE)
  withr::defer(unlink(data_dir, recursive = TRUE, force = TRUE))

  localNamespaceBinding(
    package_ns,
    "getProjectPaths",
    function(omic_type, experiment_label) list(data_dir = data_dir)
  )

  readr::write_tsv(
    tibble::tibble(
      pathway = c("path:kpn00010", "path:kpn00010", "path:kpn00020"),
      compound = c("cpd:C00001", "cpd:C00002", "cpd:C00003")
    ),
    file.path(data_dir, "species_specific_pathway_to_compound_tbl.tab")
  )
  readr::write_tsv(
    tibble::tibble(
      chebi_id = c("CHEBI:1", "CHEBI:2", "CHEBI:1", "CHEBI:3"),
      reactome_id = c("R-HSA-1", "R-HSA-1", "R-MUS-2", "R-HSA-2"),
      reactome_term = c("Path A", "Path A", "Mouse Path", "Path B"),
      organism = c("Homo sapiens", "Homo sapiens", "Mus musculus", "Homo sapiens")
    ),
    file.path(data_dir, "chebi_to_reactome.tab")
  )

  mapping_table <- data.frame(
    KEGG = c("C00001", "C00002", "C00003"),
    ChEBI = c("1", "2", "3"),
    stringsAsFactors = FALSE
  )

  options(
    multischolar.fake_keggrest.keggList = function(database) {
      stats::setNames(
        c("Glycolysis / Gluconeogenesis", "Citrate cycle"),
        c("path:map00010", "path:map00020")
      )
    },
    multischolar.fake_clusterProfiler.GSEA = function(geneList, TERM2GENE, TERM2NAME, ...) {
      if (all(grepl("^R-", TERM2GENE$term))) {
        stop("force enricher fallback")
      }

      methods::new(
        "gseaResult",
        result = data.frame(
          Description = "Glycolysis / Gluconeogenesis",
          NES = 1.8,
          p.adjust = 0.01,
          setSize = 2L,
          core_enrichment = paste(utils::head(names(geneList), 2L), collapse = "/"),
          stringsAsFactors = FALSE
        )
      )
    },
    multischolar.fake_clusterProfiler.enricher = function(gene, universe, TERM2GENE, TERM2NAME, ...) {
      methods::new(
        "enrichResult",
        result = data.frame(
          Description = TERM2NAME$name[[1]],
          p.adjust = 0.02,
          Count = length(gene),
          geneID = paste(gene, collapse = "/"),
          stringsAsFactors = FALSE
        )
      )
    }
  )

  kegg_results <- runKeggEnrichment(
    ranked_list = c(`CHEBI:1` = 2.5, `CHEBI:2` = 1.2, `CHEBI:3` = -0.7),
    mapping_table = mapping_table,
    project_dirs = list(),
    omic_type = "metabolomics",
    experiment_label = "demo",
    assay_name = "assay_lc",
    kegg_species_code = "kpn"
  )
  expect_equal(nrow(kegg_results), 1L)
  expect_identical(kegg_results$termDescription[[1]], "Glycolysis / Gluconeogenesis")
  expect_identical(kegg_results$comparison[[1]], "assay_lc")
  expect_match(kegg_results$mappedIDs[[1]], "cpd:C00001", fixed = TRUE)
  expect_equal(kegg_results$genesMapped[[1]], 2L)

  reactome_results <- runReactomeEnrichment(
    ranked_list = c(`CHEBI:1` = 1.5, `CHEBI:2` = 0.4, `CHEBI:3` = -0.7),
    mapping_table = mapping_table,
    project_dirs = list(),
    omic_type = "metabolomics",
    experiment_label = "demo",
    assay_name = "assay_gc",
    reactome_organism = "sapiens"
  )
  expect_equal(nrow(reactome_results), 1L)
  expect_identical(reactome_results$termDescription[[1]], "Path A")
  expect_identical(reactome_results$comparison[[1]], "assay_gc")
  expect_identical(reactome_results$genesMapped[[1]], 2L)
  expect_match(reactome_results$mappedIDs[[1]], "CHEBI:1/CHEBI:2", fixed = TRUE)
})

test_that("STRING validation and discovery helpers preserve current fallback branches", {
  localFakeCheckmate()

  local_mocked_bindings(
    p_load = function(...) TRUE,
    .package = "pacman"
  )

  local_mocked_bindings(
    POST = function(url, body, encode) {
      list(status = 200L, text = "{\"status\":\"error\",\"message\":\"bad request\"}")
    },
    content = function(x, as = "parsed", encoding = NULL) {
      if (!is.null(x$text)) {
        return(x$text)
      }
      "{}"
    },
    http_error = function(x) isTRUE(x$status >= 400L),
    status_code = function(x) if (is.null(x$status)) 500L else x$status,
    .package = "httr"
  )
  local_mocked_bindings(
    fromJSON = function(text, ...) {
      list(status = "error", message = "bad request")
    },
    .package = "jsonlite"
  )

  expect_error(
    submitStringDBEnrichment(
      input_data_frame = data.frame(protein = "P1", stringsAsFactors = FALSE),
      identifier_column_name = "protein",
      value_column_name = "protein",
      caller_identity = "demo-run",
      api_key = "token"
    ),
    "different"
  )

  expect_error(
    submitStringDBEnrichment(
      input_data_frame = data.frame(
        protein = c(NA_character_, NA_character_),
        score = c(NA_real_, NA_real_),
        stringsAsFactors = FALSE
      ),
      identifier_column_name = "protein",
      value_column_name = "score",
      caller_identity = "demo-run",
      api_key = "token"
    ),
    "No valid identifier"
  )

  error_status <- submitStringDBEnrichment(
    input_data_frame = data.frame(
      protein = c("P1", "P2"),
      score = c(2, -1),
      stringsAsFactors = FALSE
    ),
    identifier_column_name = "protein",
    value_column_name = "score",
    caller_identity = "demo-run",
    api_key = "token"
  )
  expect_null(error_status$job_id)
  expect_identical(error_status$submission_response$status, "error")

  species_missing_httr <- makeFunctionWithOverrides(
    getStringDbSpecies,
    list(
      requireNamespace = function(pkg, quietly = TRUE) {
        pkg != "httr"
      }
    )
  )
  expect_error(species_missing_httr(), "Package 'httr' is required")

  local_mocked_bindings(
    GET = function(url, query = NULL) list(status = 200L, text = "species"),
    content = function(x, as = "parsed", encoding = NULL) x$text,
    http_error = function(x) FALSE,
    status_code = function(x) x$status,
    .package = "httr"
  )
  local_mocked_bindings(
    fromJSON = function(text, flatten = TRUE) {
      data.frame(
        id_code = c("2", "1"),
        species_label = c("Beta species", "Alpha species"),
        short_name = c("beta", "alpha"),
        stringsAsFactors = FALSE
      )
    },
    .package = "jsonlite"
  )

  species_fun <- makeFunctionWithOverrides(
    getStringDbSpecies,
    list(syms = rlang::syms)
  )
  species <- species_fun(search_term = "species")
  expect_equal(species$species_id, c("1", "2"))
  expect_true(all(c("species_id", "compact_name") %in% names(species)))
})

test_that("STRING orchestration helpers preserve validation and empty-result branches", {
  da_obj <- methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrast = "contrast_a"),
    da_data = list(
      contrast_a = data.frame(
        Protein.Ids = c("P1", "P2"),
        fdr_qvalue = c(0.2, 0.3),
        log2FC = c(1, -1),
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE)
  )

  expect_error(
    runStringDbEnrichmentFromDAResults(
      da_results_for_enrichment = da_obj,
      ranking_method = "invalid"
    ),
    "Invalid ranking_method"
  )

  expect_error(
    runStringDbEnrichmentFromDAResults(
      da_results_for_enrichment = da_obj,
      filter_significant = TRUE,
      fdr_threshold = 0.05
    ),
    "No proteins remain after filtering"
  )

  expect_error(
    runStringDbEnrichmentFromDAResultsMultiple(
      da_results_for_enrichment = da_obj,
      contrast_names = "missing_contrast"
    ),
    "Invalid contrast names"
  )

  expect_error(
    runStringDbEnrichmentAllContrasts(
      da_analysis_results_list = NULL,
      project_dirs = list(),
      omic_type = "proteomics",
      experiment_label = "demo",
      de_analysis_results_list = NULL
    ),
    "Either da_analysis_results_list or de_analysis_results_list"
  )

  expect_error(
    runStringDbEnrichmentAllContrasts(
      da_analysis_results_list = list(),
      project_dirs = list(),
      omic_type = "proteomics",
      experiment_label = "demo"
    ),
    "non-empty list"
  )

  project_dirs <- list(proteomics_demo = list(pathway_dir = tempdir()))
  expect_null(
    plotStringDbEnrichmentResults(
      project_dirs = project_dirs,
      omic_type = "proteomics",
      experiment_label = "demo",
      enrichment_results = data.frame(),
      save_plots = FALSE,
      print_plots = FALSE
    )
  )
  expect_null(
    suppressWarnings(
      plotStringDbEnrichmentResults(
        project_dirs = project_dirs,
        omic_type = "proteomics",
        experiment_label = "demo",
        enrichment_results = data.frame(
          comparison = "A",
          category = "GO Process",
          stringsAsFactors = FALSE
        ),
        save_plots = FALSE,
        print_plots = FALSE
      )
    )
  )
})
