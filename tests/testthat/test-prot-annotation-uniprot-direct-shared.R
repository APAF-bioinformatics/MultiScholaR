# fidelity-coverage-compare: shared
library(testthat)

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
    localNamespaceBinding(env, name, bindings[[name]], .local_envir = .local_envir)
  }
}

test_that("UniProt annotation helpers preserve cache, download, and empty-table behavior", {
  package_ns <- asNamespace("MultiScholaR")
  input_tbl <- data.frame(
    Protein.Ids = c("P12345", "Q8N158-2", "bad-id"),
    stringsAsFactors = FALSE
  )

  download_dir <- tempfile("uniprot-download-")
  dir.create(download_dir)
  output_path <- file.path(download_dir, "uniprot.tsv")
  cache_dir <- tempfile("uniprot-cache-")
  dir.create(cache_dir)

  response_text <- paste(
    c(
      "Entry\tGene Names\tProtein names\tProtein existence\tAnnotation score\tGene Ontology IDs",
      "P12345\tGENE1\tProtein 1\t1\t5\tGO:1",
      "Q8N158\tGENE2\tProtein 2\t2\t7\tGO:2"
    ),
    collapse = "\n"
  )
  request_log <- new.env(parent = emptyenv())
  request_log$query <- NULL
  request_log$status <- 200L

  local_mocked_bindings(
    GET = function(url, query, ...) {
      request_log$query <- query$query
      list(status = request_log$status, text = response_text)
    },
    status_code = function(x) x$status,
    content = function(x, as = "parsed", encoding = NULL) x$text,
    timeout = function(seconds) seconds,
    .package = "httr"
  )

  downloaded <- directUniprotDownload(
    input_tbl = input_tbl,
    output_path = output_path,
    taxon_id = 9606,
    batch_size = 10,
    timeout = 60,
    api_delay = 0
  )

  request_log$status <- 500L
  failed_download <- directUniprotDownload(
    input_tbl = input_tbl,
    output_path = tempfile(fileext = ".tsv"),
    taxon_id = 9606,
    batch_size = 10,
    timeout = 60,
    api_delay = 0
  )

  standardized <- standardizeUniprotColumns(data.frame(
    Protein.existence = "1",
    Protein.names = "Protein 1",
    Gene.Names.primary = "GENE1",
    Annotation.score = "5",
    stringsAsFactors = FALSE
  ))
  empty_tbl <- createEmptyUniprotTable()

  localNamespaceBindings(
    package_ns,
    list(
      directUniprotDownload = function(...) {
        data.frame(
          Entry = "P12345",
          Gene.Ontology.IDs = "GO:1",
          Protein.names = "Protein 1",
          Protein.existence = "1",
          Annotation.score = "5",
          stringsAsFactors = FALSE
        )
      },
      uniprotGoIdToTerm = function(.data, ...) {
        transform(.data, Go.Term = "Term 1")
      }
    )
  )

  annotations_downloaded <- getUniprotAnnotations(
    input_tbl = data.frame(Protein.Ids = "P12345", stringsAsFactors = FALSE),
    cache_dir = cache_dir,
    taxon_id = 9606,
    force_download = TRUE,
    api_delay = 0
  )

  localNamespaceBinding(
    package_ns,
    "directUniprotDownload",
    function(...) stop("directUniprotDownload should not run on cache hit")
  )
  annotations_cached <- getUniprotAnnotations(
    input_tbl = data.frame(Protein.Ids = "P12345", stringsAsFactors = FALSE),
    cache_dir = cache_dir,
    taxon_id = 9606,
    force_download = FALSE,
    api_delay = 0
  )

  failure_cache_dir <- tempfile("uniprot-cache-empty-")
  dir.create(failure_cache_dir)
  localNamespaceBinding(package_ns, "directUniprotDownload", function(...) NULL)
  annotations_empty <- suppressWarnings(getUniprotAnnotations(
    input_tbl = data.frame(Protein.Ids = "P12345", stringsAsFactors = FALSE),
    cache_dir = failure_cache_dir,
    taxon_id = 9606,
    force_download = TRUE,
    api_delay = 0
  ))

  legacy_dir <- tempfile("legacy-uniprot-")
  dir.create(legacy_dir)
  legacy_tbl <- data.frame(
    Protein.Ids = "P12345",
    UNIPROT_GENENAME = "GENE1",
    stringsAsFactors = FALSE
  )
  saveRDS(legacy_tbl, file.path(legacy_dir, "uniprot_dat.rds"))
  legacy_annotations <- getUniProtAnnotation(
    input_table = data.frame(Protein.Ids = "P12345", stringsAsFactors = FALSE),
    output_dir = legacy_dir
  )

  expect_equal(nrow(downloaded), 2L)
  expect_true(file.exists(output_path))
  expect_true(grepl("organism_id:9606", request_log$query, fixed = TRUE))
  expect_true(grepl("P12345", request_log$query, fixed = TRUE))
  expect_true(grepl("Q8N158-2", request_log$query, fixed = TRUE))
  expect_false(grepl("bad-id", request_log$query, fixed = TRUE))
  expect_true("From" %in% names(downloaded))
  expect_null(failed_download)

  expect_identical(standardized$Protein_existence[[1]], "1")
  expect_identical(standardized$Protein_names[[1]], "Protein 1")
  expect_identical(standardized$gene_names[[1]], "GENE1")
  expect_identical(standardized$annotation_score[[1]], 5)
  expect_identical(colnames(empty_tbl), c("Entry", "From", "gene_names", "Protein_existence", "Protein_names", "annotation_score"))

  expect_true(file.exists(file.path(cache_dir, "uniprot_annotations.RDS")))
  expect_identical(annotations_downloaded$Protein_names[[1]], "Protein 1")
  expect_identical(annotations_cached$Protein_names[[1]], "Protein 1")
  expect_identical(nrow(annotations_empty), 0L)
  expect_identical(legacy_annotations, legacy_tbl)
})
