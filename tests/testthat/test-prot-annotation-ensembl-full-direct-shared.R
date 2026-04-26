# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

localFakeGprofiler2 <- function(env = parent.frame()) {
  old_lib_paths <- .libPaths()
  old_options <- options(multischolar.fake_gprofiler2.gconvert = NULL)

  fake_lib <- file.path(tempdir(), "multischolar-fake-gprofiler2-lib")
  fake_pkg <- file.path(tempdir(), "multischolar-fake-gprofiler2-src")
  installed_pkg <- file.path(fake_lib, "gprofiler2")

  if ("gprofiler2" %in% loadedNamespaces()) {
    try(unloadNamespace("gprofiler2"), silent = TRUE)
  }

  if (!dir.exists(installed_pkg)) {
    unlink(fake_pkg, recursive = TRUE)
    dir.create(file.path(fake_pkg, "R"), recursive = TRUE, showWarnings = FALSE)
    dir.create(fake_lib, recursive = TRUE, showWarnings = FALSE)

    writeLines(
      c(
        "Package: gprofiler2",
        "Version: 0.0.1",
        "Title: Fake gprofiler2 test double",
        "Description: Minimal namespace used by MultiScholaR annotation coverage tests.",
        "License: MIT",
        "Encoding: UTF-8"
      ),
      file.path(fake_pkg, "DESCRIPTION")
    )
    writeLines("export(gconvert)", file.path(fake_pkg, "NAMESPACE"))
    writeLines(
      c(
        "gconvert <- function(query, organism, target, mthreshold = Inf, filter_na = FALSE) {",
        "  handler <- getOption('multischolar.fake_gprofiler2.gconvert')",
        "  if (is.function(handler)) {",
        "    return(handler(query = query, organism = organism, target = target, mthreshold = mthreshold, filter_na = filter_na))",
        "  }",
        "  data.frame(input = character(), target = character(), name = character(), stringsAsFactors = FALSE)",
        "}"
      ),
      file.path(fake_pkg, "R", "gprofiler2.R")
    )
    utils::install.packages(fake_pkg, lib = fake_lib, repos = NULL, type = "source", quiet = TRUE)
  }

  .libPaths(c(fake_lib, old_lib_paths))

  withr::defer({
    options(old_options)
    if ("gprofiler2" %in% loadedNamespaces()) {
      try(unloadNamespace("gprofiler2"), silent = TRUE)
    }
    .libPaths(old_lib_paths)
  }, envir = env)

  invisible(fake_lib)
}

extractProteinIdFromHeader <- get(
  ".extractProteinIdFromHeader",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
detectEnsemblIds <- get(
  "detectEnsemblIds",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
taxonIdToGprofilerOrganism <- get(
  "taxonIdToGprofilerOrganism",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
convertEnsemblToUniprot <- get(
  "convertEnsemblToUniprot",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
getUniprotAnnotationsFull <- get(
  "getUniprotAnnotationsFull",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

test_that("annotation conversion helpers preserve header parsing, detection, taxon mapping, and gprofiler2 fallbacks", {
  expect_identical(extractProteinIdFromHeader(">sp|P11111|PROTA_HUMAN desc"), "P11111")
  expect_identical(extractProteinIdFromHeader(">generic|Q22222|desc"), "Q22222")
  expect_identical(extractProteinIdFromHeader(">ENSEMBL:ENSP000003 extra"), "ENSP000003")
  expect_identical(extractProteinIdFromHeader(">P33333 trailing words"), ">P33333")

  empty_detection <- detectEnsemblIds(c(NA_character_, ""))
  mixed_detection <- detectEnsemblIds(c("ENSP000001", "P11111"))
  ensembl_detection <- detectEnsemblIds(c("ENSP000001", "ENSP000002.1", "ENSP000003"))

  expect_false(empty_detection$is_ensembl)
  expect_identical(empty_detection$detection_rate, 0)
  expect_false(mixed_detection$is_ensembl)
  expect_true(ensembl_detection$is_ensembl)
  expect_identical(ensembl_detection$ensembl_prefix, "ENSP")

  expect_identical(taxonIdToGprofilerOrganism(9606), "hsapiens")
  expect_error(
    taxonIdToGprofilerOrganism(999999),
    "Taxon ID 999999 not supported",
    fixed = TRUE
  )

  missing_pkg_convert <- makeFunctionWithOverrides(
    convertEnsemblToUniprot,
    list(requireNamespace = function(...) FALSE)
  )
  expect_error(
    missing_pkg_convert("ENSP000001.1", "hsapiens"),
    "Package 'gprofiler2' is required",
    fixed = TRUE
  )

  localFakeGprofiler2()
  options(
    multischolar.fake_gprofiler2.gconvert = function(query, organism, target, mthreshold, filter_na) {
      data.frame(
        input = c("ENSP000001", "ENSP000001", "ENSP000002"),
        target = c("P11111", "P11111_ALT", "Q22222"),
        name = c("Protein 1", "Protein 1 alt", "Protein 2"),
        stringsAsFactors = FALSE
      )
    }
  )

  converted <- convertEnsemblToUniprot(
    c("ENSP000001.1", "ENSP000002.2", "ENSP000003.3"),
    "hsapiens"
  )

  expect_identical(converted$original_id, c("ENSP000001.1", "ENSP000002.2", "ENSP000003.3"))
  expect_identical(converted$converted_id, c("P11111", "Q22222", NA_character_))
  expect_identical(converted$conversion_status, c("success", "success", "failed"))
  expect_identical(converted$n_targets, c(2, 1, 0))

  options(multischolar.fake_gprofiler2.gconvert = function(...) stop("gconvert boom"))
  errored <- suppressWarnings(convertEnsemblToUniprot("ENSP999999.1", "hsapiens"))

  expect_identical(errored$original_id, "ENSP999999.1")
  expect_true(all(is.na(errored$converted_id)))
  expect_identical(errored$conversion_status, "error")
  expect_identical(errored$n_targets, 0)
})

test_that("getUniprotAnnotationsFull preserves direct UniProt lookup, validation, and metadata attachment", {
  expect_error(
    getUniprotAnnotationsFull(NULL, "Protein.Group", tempfile("annot-cache-")),
    "data_tbl must be a non-null data frame",
    fixed = TRUE
  )
  expect_error(
    getUniprotAnnotationsFull(data.frame(Protein.Group = character()), "Protein.Group", tempfile("annot-cache-")),
    "data_tbl cannot be empty",
    fixed = TRUE
  )
  expect_error(
    getUniprotAnnotationsFull(data.frame(other = "P1"), "Protein.Group", tempfile("annot-cache-")),
    "Protein ID column 'Protein.Group' not found",
    fixed = TRUE
  )

  captured_input <- NULL
  direct_helper <- makeFunctionWithOverrides(
    getUniprotAnnotationsFull,
    list(
      .extractProteinIdFromHeader = extractProteinIdFromHeader,
      normalizeUniprotAccession = function(string, remove_isoform = TRUE) {
        sub("-.*$", "", string)
      },
      getUniprotAnnotations = function(input_tbl, cache_dir, taxon_id, progress_callback = NULL) {
        captured_input <<- input_tbl
        data.frame(
          Entry = input_tbl$Protein.Ids,
          gene_name = paste0("GENE_", seq_len(nrow(input_tbl))),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  cache_dir <- tempfile("direct-uniprot-cache-")
  result <- direct_helper(
    data_tbl = data.frame(
      Protein.Group = c(
        ">sp|P11111|PROTA_HUMAN desc",
        ">generic|Q22222|desc",
        ">TREMBL:A2A5Y0 isoform",
        "",
        NA_character_
      ),
      stringsAsFactors = FALSE
    ),
    protein_id_column = "Protein.Group",
    cache_dir = cache_dir,
    taxon_id = 9606
  )

  expect_true(dir.exists(cache_dir))
  expect_identical(captured_input$Protein.Ids, c("P11111", "Q22222", "A2A5Y0"))
  expect_false(attr(result, "ensembl_conversion_applied"))
  expect_identical(attr(result, "annotation_source"), "DATA_OPTIMIZED")
  expect_identical(attr(result, "original_protein_groups"), 5L)
  expect_identical(attr(result, "unique_proteins_processed"), 3L)
  expect_identical(attr(result, "taxon_id"), 9606)
  expect_true(inherits(attr(result, "processing_timestamp"), "POSIXct"))
})

test_that("getUniprotAnnotationsFull preserves Ensembl conversion merge behavior and minimal fallback output", {
  success_helper <- makeFunctionWithOverrides(
    getUniprotAnnotationsFull,
    list(
      .extractProteinIdFromHeader = extractProteinIdFromHeader,
      normalizeUniprotAccession = function(string, remove_isoform = TRUE) {
        sub("\\.\\d+$", "", string)
      },
      detectEnsemblIds = function(protein_ids) {
        list(is_ensembl = TRUE, ensembl_prefix = "ENSP", detection_rate = 1)
      },
      taxonIdToGprofilerOrganism = function(taxon_id) "hsapiens",
      convertEnsemblToUniprot = function(ensembl_ids, organism_code) {
        data.frame(
          original_id = ensembl_ids,
          converted_id = c("P11111", NA_character_),
          conversion_status = c("success", "failed"),
          n_targets = c(1, 0),
          stringsAsFactors = FALSE
        )
      },
      getUniprotAnnotations = function(input_tbl, cache_dir, taxon_id, progress_callback = NULL) {
        data.frame(
          Entry = c("P11111", "P99999"),
          annotation = c("Alpha", "Beta"),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  success_result <- suppressWarnings(success_helper(
    data_tbl = data.frame(
      Protein.Ids = c("ENSP000001.1", "ENSP000002.2"),
      stringsAsFactors = FALSE
    ),
    protein_id_column = "Protein.Ids",
    cache_dir = tempfile("ensembl-cache-success-"),
    taxon_id = 9606
  ))

  expect_true(attr(success_result, "ensembl_conversion_applied"))
  expect_equal(attr(success_result, "ensembl_conversion_success_rate"), 0.5)
  expect_true("Original_Ensembl_ID" %in% names(success_result))
  expect_identical(
    success_result$Original_Ensembl_ID[success_result$Entry == "P11111"],
    "ENSP000001"
  )
  expect_identical(
    success_result$Original_Ensembl_ID[success_result$Entry == "P99999"],
    "P99999"
  )

  fallback_helper <- makeFunctionWithOverrides(
    getUniprotAnnotationsFull,
    list(
      .extractProteinIdFromHeader = extractProteinIdFromHeader,
      normalizeUniprotAccession = function(string, remove_isoform = TRUE) {
        sub("\\.\\d+$", "", string)
      },
      detectEnsemblIds = function(protein_ids) {
        list(is_ensembl = TRUE, ensembl_prefix = "ENSP", detection_rate = 1)
      },
      taxonIdToGprofilerOrganism = function(taxon_id) {
        stop("taxonomy boom")
      },
      getUniprotAnnotations = function(input_tbl, cache_dir, taxon_id, progress_callback = NULL) {
        stop("annotation boom")
      }
    )
  )

  fallback_result <- suppressWarnings(fallback_helper(
    data_tbl = data.frame(
      Protein.Ids = c("ENSP000001.1", "ENSP000002.2"),
      stringsAsFactors = FALSE
    ),
    protein_id_column = "Protein.Ids",
    cache_dir = tempfile("ensembl-cache-fallback-"),
    taxon_id = 9606
  ))

  expect_identical(fallback_result$Protein.Ids, c("ENSP000001", "ENSP000002"))
  expect_true(all(is.na(fallback_result$UNIPROT_GENENAME)))
})
