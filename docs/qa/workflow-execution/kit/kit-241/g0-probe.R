# KIT-241 G0 and W1 probe (ticket section 8): the proteomics enrichment download writer
# (writeProtEnrichResultsDownloadArchive) under the installed bundle's R, in the environment the
# application gives R. Run with the installed Rscript --vanilla; no test framework, no network.
#
# Usage: Rscript --vanilla g0-probe.R <output folder> <host label> <G0|W1>
#
# Prints one line per observation, "KIT241<TAB><gate><TAB><host><TAB><key><TAB><value>":
#   (1) the environment: R, R_ZIPCMD, Sys.which("zip"), PATH, .libPaths(), the zip package, tempdir();
#   (2) case A: the writer with its own default zipFn;
#   (3) case B: the same with zipFn = zip::zipr (the fix, exercised on the code under test);
#   (4) Windows only, case C: case B under a folder with non-ASCII letters (P09).
# Per case: the returned status, warnings, error, whether the archive exists, its size, its listing
# with R's unzip or the error, and whether each extracted member equals the file the writer wrote.
# Each existing archive and the written members are copied to <output folder>/<case>/ (R removes
# its session temporary folder on exit) and listed in <output folder>/cases.tsv for the Windows
# step's own extractors. A failure inside a case is an observation; a failure outside the cases
# (for example the package does not load) stops the probe before its last line
# "<gate> END <host>", which the calling step requires (instrument fault, VAL-54).

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
outDir <- args[[1L]]
host <- args[[2L]]
gate <- args[[3L]]
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

flat <- \(x) gsub("[\t\r\n]+", " | ", paste(x, collapse = " ; "))
say <- \(key, value) cat("KIT241", gate, host, key, paste0(flat(value), "\n"), sep = "\t")

# (1) The environment.
isWindows <- identical(.Platform$OS.type, "windows")
say("R", R.version.string)
say("os", paste(Sys.info()[["sysname"]], Sys.info()[["release"]]))
zipCmd <- Sys.getenv("R_ZIPCMD", unset = NA)
say("R_ZIPCMD", if (is.na(zipCmd)) "unset" else paste0("set: '", zipCmd, "'"))
say("Sys.which(zip)", if (nzchar(Sys.which("zip"))) Sys.which("zip") else "<none>")
if (!isWindows) say("file.exists(/usr/bin/zip)", file.exists("/usr/bin/zip"))
say("PATH", Sys.getenv("PATH"))
say(".libPaths()", .libPaths())
say("packageVersion(zip)", tryCatch(as.character(utils::packageVersion("zip")), error = conditionMessage))
say("find.package(zip)", tryCatch(find.package("zip"), error = conditionMessage))
say("packageVersion(MultiScholaR)", as.character(utils::packageVersion("MultiScholaR")))
say("MultiScholaR SourceCommit",
    utils::packageDescription("MultiScholaR", fields = "Config/MultiScholaR/SourceCommit"))
say("tempdir()", tempdir())
say("nchar(tempdir())", nchar(tempdir()))
if (isWindows) {
  say("R_RTOOLS45_PATH", Sys.getenv("R_RTOOLS45_PATH", unset = "<unset>"))
  say("RTOOLS45_HOME", Sys.getenv("RTOOLS45_HOME", unset = "<unset>"))
}

writer <- utils::getFromNamespace("writeProtEnrichResultsDownloadArchive", "MultiScholaR")
say("writer default zipFn", deparse(formals(writer)$zipFn))

# One-row synthetic g:Profiler table (ticket section 4).
gprofilerTable <- data.frame(
  term_id = "GO:0006915", term_name = "apoptotic process", source = "GO:BP",
  p_value = 0.001, intersection_size = 3L, stringsAsFactors = FALSE
)
expectedMembers <- c("enrichment_analysis_summary.txt", "gprofiler2_results.tsv")
casesFile <- file.path(outDir, "cases.tsv")
casesRows <- "case\tarchive\twritten_dir\tmembers"

runCase <- function(label, zipFn, parent) {
  caseDir <- file.path(parent, paste0("kit241-case-", label))
  unlink(caseDir, recursive = TRUE)
  dir.create(caseDir, recursive = TRUE, showWarnings = FALSE)
  say(paste0(label, ".folder"), caseDir)
  say(paste0(label, ".nchar(folder)"), nchar(caseDir))
  archive <- file.path(caseDir, "Enrichment_results_B_vs_A.zip")
  callArgs <- list(
    file = archive, selectedContrast = "B_vs_A",
    methodInfo = list(method = "gprofiler2", species_name = "Homo sapiens"),
    organismTaxid = "9606", upCutoff = 0.5, downCutoff = -0.5, qCutoff = 0.05,
    gprofilerResults = gprofilerTable, tempDir = caseDir
  )
  if (!is.null(zipFn)) callArgs$zipFn <- zipFn
  caught <- character()
  raised <- NA_character_
  result <- withCallingHandlers(
    tryCatch(do.call(writer, callArgs), error = \(e) {
      raised <<- conditionMessage(e)
      NULL
    }),
    warning = \(w) {
      caught <<- c(caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  say(paste0(label, ".status"), if (is.null(result)) "<no value: an error left the writer>" else result$status)
  say(paste0(label, ".warnings"), if (length(caught)) caught else "<none>")
  say(paste0(label, ".error"), if (!is.na(raised)) paste("raised:", raised) else
    if (!is.null(result$error)) paste("returned:", result$error) else "<none>")
  exists <- file.exists(archive)
  say(paste0(label, ".archive_exists"), exists)
  say(paste0(label, ".archive_size"), if (exists) file.size(archive) else NA)
  listing <- if (exists) tryCatch(sort(utils::unzip(archive, list = TRUE)$Name),
    error = \(e) paste("unzip error:", conditionMessage(e))) else "<no archive>"
  say(paste0(label, ".listing"), listing)
  say(paste0(label, ".listing_is_expected"), identical(listing, sort(expectedMembers)))
  written <- if (!is.null(result)) result$files else file.path(caseDir, expectedMembers)
  if (exists) {
    extracted <- file.path(caseDir, "r-unzip")
    equal <- tryCatch({
      utils::unzip(archive, exdir = extracted)
      same <- vapply(written, \(f) {
        target <- file.path(extracted, basename(f))
        file.exists(target) && identical(unname(tools::md5sum(f)), unname(tools::md5sum(target)))
      }, logical(1))
      if (all(same)) "members equal" else paste("differ:", paste(basename(written)[!same], collapse = ", "))
    }, error = \(e) paste("error:", conditionMessage(e)))
    say(paste0(label, ".R_unzip_extract"), equal)
    keep <- file.path(outDir, label)
    dir.create(file.path(keep, "written"), recursive = TRUE, showWarnings = FALSE)
    file.copy(archive, file.path(keep, "download.zip"), overwrite = TRUE)
    file.copy(written, file.path(keep, "written"), overwrite = TRUE)
    casesRows <<- c(casesRows, paste(label, file.path(keep, "download.zip"), file.path(keep, "written"),
      paste(basename(written), collapse = ";"), sep = "\t"))
  }
  invisible(NULL)
}

zipr <- \(zipfile, files, flags) zip::zipr(zipfile = zipfile, files = files)
runCase("A", NULL, tempdir())
runCase("B", zipr, tempdir())
if (isWindows) {
  nonAscii <- file.path(tempdir(), "Ü-数据")
  dir.create(nonAscii, showWarnings = FALSE)
  runCase("C", zipr, nonAscii)
}

con <- file(casesFile, open = "w", encoding = "UTF-8")
writeLines(casesRows, con)
close(con)
say("cases.tsv", casesFile)
cat(gate, " END ", host, "\n", sep = "")
