# DIA-025: bounded ownership and exact-payload checks, not a GUI journey.

test_that("repair checks keep the explicitly installed namespace authoritative", {
  library <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")
  skip_if(!nzchar(library), "installed check service only")
  expect_identical(normalizePath(getNamespaceInfo(asNamespace("MultiScholaR"), "path")),
    normalizePath(file.path(library, "MultiScholaR")))
})

diaRepairReference <- function() {
  path <- Sys.getenv("MULTISCHOLAR_DIA008_B0_LIBRARY")
  skip_if(!nzchar(path), "requires the frozen installed predecessor")
  db <- new.env(parent = baseenv())
  lazyLoad(file.path(path, "MultiScholaR", "R", "MultiScholaR"), envir = db)
  ref <- new.env(parent = asNamespace("MultiScholaR"))
  for (name in c(".peptideQcCanonicalise", ".peptideQcDigest", ".peptideQcIdentityKey",
      "generateProtDAVolcanoPlotGlimma", "loadProtNormImageAsPlot")) {
    fn <- get(name, db)
    environment(fn) <- ref
    ref[[name]] <- fn
  }
  ref
}

diaRepairMeasure <- function(fn) {
  log <- withr::local_tempfile()
  initial_heap <- gc(reset = TRUE)
  started <- proc.time()[["elapsed"]]
  utils::Rprofmem(log, threshold = 10000)
  on.exit(utils::Rprofmem(NULL), add = TRUE)
  value <- fn()
  utils::Rprofmem(NULL)
  elapsed <- proc.time()[["elapsed"]] - started
  heap <- gc()
  stopifnot(file.size(log) < 64 * 1024^2)
  lines <- readLines(log)
  allocations <- as.numeric(sub(" .*", "", grep("^[0-9]+ :", lines, value = TRUE)))
  value_hash <- digest::digest(value, algo = "sha256")
  list(value_hash = value_hash, observation = list(elapsed = elapsed,
    initial_vector_heap_mib = unname(initial_heap["Vcells", 2L]),
    vector_heap_high_water_mib = unname(heap["Vcells", 6L]),
    allocations_over_10kb = sum(allocations), largest_allocation = max(c(0, allocations))))
}

diaRepairRecord <- function(name, observation) {
  receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
  if (nzchar(receipt)) jsonlite::write_json(observation,
    file.path(dirname(receipt), paste0("repair-", name, ".json")),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
}

test_that("step-owned workers terminate but saved S4 and borrowed pools remain usable", {
    skip("quarantined: Q-22 (KIT-048) the legacy-mode manager path these tests drove is retired; the owner is the one path")
  object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
  paths <- list(source_dir = withr::local_tempdir(), protein_qc_dir = withr::local_tempdir())
  workflow <- new.env(parent = emptyenv())
  workflow$state_manager <- WorkflowState$new()
  workflow$state_manager$saveState("raw", object, list(), "raw")
  owned <- list()
  create <- function(n) {
    pool <- multidplyr::new_cluster(n)
    owned[[length(owned) + 1L]] <<- pool
    pool
  }
  on.exit(for (pool in owned) stopDIAWorkerCluster(pool), add = TRUE)
  # Keep this lifecycle fixture small; real numeric parallel dispatch is checked below.
  result <- runProteinReplicateFilterApplyStep(workflow, paths, "group", 2L,
    newClusterFn = create, isTestModeFn = function() FALSE,
    removeProteinsWithOnlyOneReplicateFn = function(theObject, core_utilisation, grouping_variable) {
      expect_true(all(vapply(core_utilisation, \(worker) worker$is_alive(), logical(1))))
      theObject
    })
  expect_false(any(vapply(owned[[1L]], \(worker) worker$is_alive(), logical(1))))
  spec <- newDIAWorkerSpec(2L)
  expect_identical(result$filteredS4@args$removeProteinsWithOnlyOneReplicate$core_utilisation, spec)
  expect_identical(result$filteredS4@protein_quant_table, object@protein_quant_table)
  session_file <- withr::local_tempfile()
  saveRDS(workflow$state_manager, session_file)
  restored <- readRDS(session_file)$getState()
  expect_identical(restored, result$filteredS4)
  baseline <- removeProteinsWithOnlyOneReplicate(object, NA, "group")
  actual <- removeProteinsWithOnlyOneReplicate(restored, grouping_variable = "group")
  expect_identical(actual@protein_quant_table, baseline@protein_quant_table)
  expect_identical(actual@design_matrix, baseline@design_matrix)
  expect_identical(actual@args$removeProteinsWithOnlyOneReplicate$core_utilisation, spec)
  expect_error(runProteinReplicateFilterApplyStep(workflow, paths, "group", 2L,
    newClusterFn = create, isTestModeFn = function() FALSE,
    removeProteinsWithOnlyOneReplicateFn = function(...) stop("scientific failure")), "scientific failure")
  expect_false(any(vapply(owned[[2L]], \(worker) worker$is_alive(), logical(1))))
  expect_error(runProteinReplicateFilterApplyStep(workflow, paths, "group", 2L,
    newClusterFn = create, isTestModeFn = function() FALSE,
    removeProteinsWithOnlyOneReplicateFn = function(theObject, ...) theObject,
    writeTsvFn = function(...) stop("publication failure")), "publication failure")
  expect_false(any(vapply(owned[[3L]], \(worker) worker$is_alive(), logical(1))))
  testthat::with_mocked_bindings({
    expect_error(removeProteinsWithOnlyOneReplicate(restored, grouping_variable = "group"),
      "restored calculation failure")
  }, startDIAWorkerCluster = function(spec) create(spec$workers),
    removeProteinsWithOnlyOneReplicateHelper = function(...) stop("restored calculation failure"))
  expect_false(any(vapply(owned[[4L]], \(worker) worker$is_alive(), logical(1))))
  borrowed <- startDIAWorkerCluster(spec)
  on.exit(stopDIAWorkerCluster(borrowed), add = TRUE)
  borrowed_result <- removeProteinsWithOnlyOneReplicate(object, borrowed, "group")
  expect_identical(borrowed_result@protein_quant_table, baseline@protein_quant_table)
  expect_true(all(vapply(borrowed, \(worker) worker$is_alive(), logical(1))))
  workflow$state_manager$revertToState("raw")
  expect_identical(workflow$state_manager$getState(), object)
})

test_that("audit keys preserve frozen full-column encodings and digests", {
  reference <- diaRepairReference()
  set.seed(25025)
  n <- 40003L
  latin <- iconv("caf\u00e9", to = "latin1")
  data <- data.frame(id = sample(c("a\u001fb", "<NA>", NA, "", "\u03b1", latin), n, TRUE),
    value = rep(c(NA, NaN, Inf, -Inf, -0, 1e-100, 1e100, pi), length.out = n),
    integer = seq_len(n), flag = rep(c(TRUE, FALSE, NA), length.out = n))
  attr(data, "provenance") <- list(source = "synthetic", factor = factor("x", c("x", "unused")))
  variants <- list(data, data[n:1L, ], data[0, ], data[, FALSE],
    transform(data, group = factor(id)), transform(data, date = as.Date("2026-01-01")))
  for (value in variants) {
    expect_identical(.peptideQcIdentityKey(value), reference$.peptideQcIdentityKey(value))
    expect_identical(.peptideQcCanonicalise(value), reference$.peptideQcCanonicalise(value))
    expect_identical(.peptideQcDigest(value), reference$.peptideQcDigest(value))
  }
  wide <- as.data.frame(rep(data, 4L), check.names = TRUE)
  before <- diaRepairMeasure(function() reference$.peptideQcIdentityKey(wide))
  after <- diaRepairMeasure(function() .peptideQcIdentityKey(wide))
  expect_identical(after$value_hash, before$value_hash)
  diaRepairRecord("audit-keys", list(rows = n, columns = ncol(wide),
    before = before$observation, after = after$observation,
    meaning = "isolated key construction, not complete Apply or process peak"))
})

test_that("composite image plots retain the raster but not the decoded image frame", {
  reference <- diaRepairReference()
  path <- withr::local_tempfile(fileext = ".png")
  png::writePNG(array(seq(0, 1, length.out = 256L * 256L * 3L), c(256L, 256L, 3L)), path)
  before <- reference$loadProtNormImageAsPlot(path)
  after <- loadProtNormImageAsPlot(path)
  before_retains_image <- exists("img", before$plot_env, inherits = FALSE)
  expect_true(before_retains_image || identical(before$plot_env, baseenv()))
  expect_identical(after$plot_env, baseenv())
  expect_identical(after$layers[[1L]]$geom_params$grob$raster,
    before$layers[[1L]]$geom_params$grob$raster)
  outputs <- lapply(list(before, after), function(plot) {
    output <- withr::local_tempfile(fileext = ".png")
    ggplot2::ggsave(output, plot, width = 3, height = 3, dpi = 72)
    png::readPNG(output)
  })
  expect_identical(outputs[[1L]], outputs[[2L]])
  saved <- withr::local_tempfile()
  saveRDS(after, saved)
  expect_silent(ggplot2::ggplotGrob(readRDS(saved)))
  diaRepairRecord("plot-owners", list(predecessor_retains_decoded_frame = before_retains_image,
    decoded_image_bytes = if (before_retains_image)
      as.numeric(object.size(get("img", before$plot_env))) else 0,
    candidate_retains_decoded_frame = FALSE, pixels_identical = TRUE))
})

test_that("removing discarded Glimma preparation preserves the complete widget payload", {
  reference <- diaRepairReference()
  object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
  # Enlarge only this bounded view fixture, not the frozen GUI workload.
  rows <- rep(seq_len(nrow(object@protein_quant_table)), length.out = 2048L)
  object@protein_quant_table <- object@protein_quant_table[rows, ]
  object@protein_quant_table$Protein.Group <- paste0("P", seq_along(rows))
  ids <- object@protein_quant_table$Protein.Group
  data <- data.frame(Protein.Ids = ids, comparison = "T_vs_C",
    log2FC = rep(c(-2, 0, 2), length.out = length(ids)),
    fdr_qvalue = rep(c(0, 0.01, 0.8), length.out = length(ids)), raw_pvalue = 0.01,
    custom.display = factor(rep(c("a", NA), length.out = length(ids))))
  object@protein_quant_table <- dplyr::rename(object@protein_quant_table, Protein.Ids = Protein.Group)
  object@protein_id_column <- "Protein.Ids"
  input <- list(da_proteins_long = data, theObject = object)
  # Warm shared namespaces; measurements are isolated allocation evidence only.
  invisible(generateProtDAVolcanoPlotGlimma(input, "T_vs_C", args_row_id = "Protein.Ids"))
  expected <- reference$generateProtDAVolcanoPlotGlimma(
    input, "T_vs_C", args_row_id = "Protein.Ids", display_columns = "custom.display")
  actual <- generateProtDAVolcanoPlotGlimma(
    input, "T_vs_C", args_row_id = "Protein.Ids", display_columns = "custom.display")
  expect_identical(actual$x, expected$x)
  expect_identical(actual$prepend, expected$prepend)
  expect_identical(actual$dependencies, expected$dependencies)
  rm(expected, actual)
  before <- diaRepairMeasure(function() reference$generateProtDAVolcanoPlotGlimma(
    input, "T_vs_C", args_row_id = "Protein.Ids", display_columns = "custom.display")$x)
  after <- diaRepairMeasure(function() generateProtDAVolcanoPlotGlimma(
    input, "T_vs_C", args_row_id = "Protein.Ids", display_columns = "custom.display")$x)
  expect_identical(after$value_hash, before$value_hash)
  diaRepairRecord("glimma-preparation", list(rows = length(ids),
    before = before$observation, after = after$observation, widget_payload_identical = TRUE))
})


test_that("the unrelated density S4 dispatch failure is checked against the frozen predecessor", {
  reference <- Sys.getenv("MULTISCHOLAR_DIA008_B0_LIBRARY")
  skip_if(!nzchar(reference), "frozen installed predecessor required")
  path <- withr::local_tempfile()
  saveRDS(makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")), path)
  check <- function(library, file) {
    .libPaths(unique(c(library, .libPaths())))
    ns <- loadNamespace("MultiScholaR", lib.loc = library)
    tryCatch({ ns$plotDensity(readRDS(file), grouping_variable = "group"); NULL },
      error = function(error) conditionMessage(error))
  }
  before <- callr::r(check, list(reference, path))
  after <- callr::r(check, list(Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY"), path))
  expect_identical(after, before)
  diaRepairRecord("existing-density-dispatch", list(before = before, after = after,
    changed_by_memory_repair = FALSE))
})
