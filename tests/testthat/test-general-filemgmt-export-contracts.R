# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  funOverride <- fun
  environment(funOverride) <- list2env(replacements, parent = environment(fun))
  funOverride
}

withCleanGlobalObjects <- function(objectNames, code) {
  hadExisting <- vapply(
    objectNames,
    function(name) exists(name, envir = .GlobalEnv, inherits = FALSE),
    logical(1)
  )
  oldValues <- lapply(seq_along(objectNames), function(i) {
    if (hadExisting[[i]]) {
      get(objectNames[[i]], envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
  })
  names(oldValues) <- objectNames

  on.exit({
    for (name in rev(objectNames)) {
      if (hadExisting[[name]]) {
        assign(name, oldValues[[name]], envir = .GlobalEnv)
      } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  force(code)
}

test_that("savePlot saves an RDS and uses cairo_pdf for PDF output", {
  fixtureDir <- tempfile("general-filemgmt-export-single-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  plotObject <- structure(list(id = "single-plot"), class = c("gg", "ggplot"))
  ggsaveCalls <- list()

  local_mocked_bindings(
    ggsave = function(filename, plot, device, width, height, ...) {
      ggsaveCalls[[length(ggsaveCalls) + 1]] <<- list(
        filename = filename,
        plot = plot,
        device = device,
        width = width,
        height = height,
        dots = list(...)
      )
      invisible(filename)
    },
    .env = environment(savePlot)
  )

  savePlot(
    plot = plotObject,
    base_path = fixtureDir,
    plot_name = "summary_plot",
    formats = c("pdf", "png"),
    width = 4,
    height = 5,
    bg = "white"
  )

  expect_true(file.exists(file.path(fixtureDir, "summary_plot.rds")))
  expect_identical(readRDS(file.path(fixtureDir, "summary_plot.rds")), plotObject)
  expect_length(ggsaveCalls, 2)
  expect_identical(
    vapply(ggsaveCalls, `[[`, character(1), "filename"),
    c(
      file.path(fixtureDir, "summary_plot.pdf"),
      file.path(fixtureDir, "summary_plot.png")
    )
  )
  expect_true(identical(ggsaveCalls[[1]]$device, grDevices::cairo_pdf))
  expect_identical(ggsaveCalls[[2]]$device, "png")
  expect_equal(ggsaveCalls[[1]]$width, 4)
  expect_equal(ggsaveCalls[[1]]$height, 5)
  expect_identical(ggsaveCalls[[1]]$dots$bg, "white")
  expect_identical(ggsaveCalls[[2]]$plot, plotObject)
})

test_that("savePlot writes each gg plot from a named list and skips non-gg members", {
  fixtureDir <- tempfile("general-filemgmt-export-list-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  plotList <- list(
    alpha = structure(list(id = "alpha"), class = c("gg", "ggplot")),
    notes = "skip-me",
    beta = structure(list(id = "beta"), class = c("gg", "ggplot"))
  )
  ggsaveCalls <- list()

  local_mocked_bindings(
    ggsave = function(filename, plot, device, width, height, ...) {
      ggsaveCalls[[length(ggsaveCalls) + 1]] <<- list(
        filename = filename,
        plot = plot,
        device = device,
        width = width,
        height = height,
        dots = list(...)
      )
      invisible(filename)
    },
    .env = environment(savePlot)
  )

  savePlot(
    plot = plotList,
    base_path = fixtureDir,
    plot_name = "bundle",
    formats = "png",
    width = 6,
    height = 3
  )

  expect_true(file.exists(file.path(fixtureDir, "bundle.rds")))
  expect_identical(readRDS(file.path(fixtureDir, "bundle.rds")), plotList)
  expect_length(ggsaveCalls, 2)
  expect_identical(
    vapply(ggsaveCalls, `[[`, character(1), "filename"),
    c(
      file.path(fixtureDir, "bundle_alpha.png"),
      file.path(fixtureDir, "bundle_beta.png")
    )
  )
  expect_true(all(vapply(ggsaveCalls, function(call) identical(call$device, "png"), logical(1))))
  expect_identical(ggsaveCalls[[1]]$plot, plotList$alpha)
  expect_identical(ggsaveCalls[[2]]$plot, plotList$beta)
})

test_that("save_plot delegates to savePlot with the same arguments", {
  recordedCall <- NULL
  savePlotAlias <- makeFunctionWithOverrides(
    save_plot,
    list(
      savePlot = function(plot, base_path, plot_name, formats, width, height, ...) {
        recordedCall <<- list(
          plot = plot,
          base_path = base_path,
          plot_name = plot_name,
          formats = formats,
          width = width,
          height = height,
          dots = list(...)
        )
        invisible(NULL)
      }
    )
  )

  plotObject <- structure(list(id = "delegated-plot"), class = c("gg", "ggplot"))

  savePlotAlias(
    plot = plotObject,
    base_path = "/tmp/export-path",
    plot_name = "delegated",
    formats = "svg",
    width = 9,
    height = 4,
    dpi = 300
  )

  expect_identical(recordedCall$plot, plotObject)
  expect_identical(recordedCall$base_path, "/tmp/export-path")
  expect_identical(recordedCall$plot_name, "delegated")
  expect_identical(recordedCall$formats, "svg")
  expect_equal(recordedCall$width, 9)
  expect_equal(recordedCall$height, 4)
  expect_identical(recordedCall$dots$dpi, 300)
})

test_that("write_results writes tabular output under results_dir/protein_qc", {
  fixtureDir <- tempfile("general-filemgmt-write-results-")
  proteinQcDir <- file.path(fixtureDir, "protein_qc")
  dir.create(proteinQcDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  exportData <- data.frame(
    protein = c("P1", "P2"),
    score = c(1.5, 2.5),
    stringsAsFactors = FALSE
  )

  withCleanGlobalObjects("results_dir", {
    assign("results_dir", fixtureDir, envir = .GlobalEnv)

    write_results(exportData, "protein-qc.tsv")

    writtenPath <- file.path(proteinQcDir, "protein-qc.tsv")
    expect_true(file.exists(writtenPath))
    expect_identical(
      utils::read.delim(writtenPath, stringsAsFactors = FALSE),
      exportData
    )
  })
})
