# KIT-241: the proteomics enrichment "Download All Results" archive is written without an
# external zip program. E rows make the external program unreachable (R_ZIPCMD names a path that
# does not exist; utils::zip reads R_ZIPCMD first), so they falsify on a runner that has zip.
# G1 and E4 are guards that pass on base by design; the failing-first rule applies to E1 to E3.

kit241GprofilerTable <- function() {
  data.frame(
    term_id = "GO:0006915", term_name = "apoptotic process", source = "GO:BP",
    p_value = 0.001, intersection_size = 3L, stringsAsFactors = FALSE
  )
}

kit241WriteArchive <- function(archive, tempDir, ...) {
  writeProtEnrichResultsDownloadArchive(
    file = archive,
    selectedContrast = "B_vs_A",
    methodInfo = list(method = "gprofiler2", species_name = "Homo sapiens"),
    organismTaxid = "9606",
    upCutoff = 0.5,
    downCutoff = -0.5,
    qCutoff = 0.05,
    tempDir = tempDir,
    ...
  )
}

kit241Listing <- function(archive) sort(utils::unzip(archive, list = TRUE)$Name)

kit241NoZipProgram <- function(env = parent.frame()) {
  withr::local_envvar(R_ZIPCMD = file.path(withr::local_tempdir(.local_envir = env), "no-zip-program"),
    .local_envir = env)
}

test_that("E1: the archive is written with no zip program (one g:Profiler table)", {
  kit241NoZipProgram()
  archive <- file.path(withr::local_tempdir(), "archive.zip")
  result <- kit241WriteArchive(archive, withr::local_tempdir(), gprofilerResults = kit241GprofilerTable())

  expect_true(file.exists(archive))
  if (file.exists(archive)) {
    expect_identical(kit241Listing(archive), c("enrichment_analysis_summary.txt", "gprofiler2_results.tsv"))
  }
  expect_identical(result$status, "ok")
})

test_that("E2: the archive is written with no zip program (no result tables)", {
  kit241NoZipProgram()
  archive <- file.path(withr::local_tempdir(), "archive.zip")
  result <- kit241WriteArchive(archive, withr::local_tempdir())

  expect_true(file.exists(archive))
  if (file.exists(archive)) {
    expect_identical(kit241Listing(archive), "enrichment_analysis_summary.txt")
  }
  expect_identical(result$status, "ok")
})

test_that("E3: the error note is archived with no zip program", {
  kit241NoZipProgram()
  archive <- file.path(withr::local_tempdir(), "archive.zip")
  result <- kit241WriteArchive(
    archive, withr::local_tempdir(),
    gprofilerResults = kit241GprofilerTable(),
    writeTsvFn = function(...) stop("write failed")
  )

  expect_true(file.exists(archive))
  if (file.exists(archive)) {
    expect_identical(kit241Listing(archive), "download_error.txt")
  }
  expect_identical(result$status, "error")
  expect_identical(result$error, "write failed")
})

test_that("E4 (guard): an archiver that fails is never returned as status ok", {
  archive <- file.path(withr::local_tempdir(), "archive.zip")
  expect_error(
    kit241WriteArchive(
      archive, withr::local_tempdir(),
      gprofilerResults = kit241GprofilerTable(),
      zipFn = function(zipfile, files, flags) stop("archive not written")
    ),
    "archive not written"
  )
})

test_that("G1 (guard): members are stored under their base names and the summary is what was written", {
  archive <- file.path(withr::local_tempdir(), "archive.zip")
  tempDir <- withr::local_tempdir()
  result <- kit241WriteArchive(
    archive, tempDir,
    gprofilerResults = kit241GprofilerTable(),
    clusterprofilerResults = data.frame(ID = "GO:0006915", Description = "apoptotic process",
      p.adjust = 0.002, stringsAsFactors = FALSE)
  )

  expect_identical(result$status, "ok")
  expect_true(file.exists(archive))
  if (file.exists(archive)) {
    expect_identical(kit241Listing(archive), sort(basename(result$files)))
    extracted <- withr::local_tempdir()
    utils::unzip(archive, files = "enrichment_analysis_summary.txt", exdir = extracted)
    expect_identical(
      readLines(file.path(extracted, "enrichment_analysis_summary.txt")),
      readLines(file.path(tempDir, "enrichment_analysis_summary.txt"))
    )
  }
})
