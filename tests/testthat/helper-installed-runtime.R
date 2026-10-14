# KIT-034: skip conditions for tests that need the installed package or a browser, replacing the Q-04 and Q-14
# quarantine entries. Under `devtools::load_all` no installed `MultiScholaR` binary exists for a fresh R process or a
# worker to load; these tests run against the installed proof library (set `MULTISCHOLAR_INSTALLED_RUNTIME=1`
# with the library on `R_LIBS`) and are recorded in `docs/qa/workflow-execution/kit/kit-034/`.
installedRuntimeAvailable <- function() {
  identical(Sys.getenv("MULTISCHOLAR_INSTALLED_RUNTIME"), "1") &&
    nzchar(system.file("R", "MultiScholaR.rdb", package = "MultiScholaR"))
}

skipIfNoInstalledRuntime <- function() {
  testthat::skip_if_not(installedRuntimeAvailable(),
    "needs the installed package (MULTISCHOLAR_INSTALLED_RUNTIME=1 with the proof library on R_LIBS); see docs/qa/workflow-execution/kit/kit-034/")
}

browserRuntimeAvailable <- function() {
  installedRuntimeAvailable() && requireNamespace("chromote", quietly = TRUE) &&
    !is.null(tryCatch(chromote::find_chrome(), error = function(error) NULL))
}

skipIfNoBrowserRuntime <- function() {
  testthat::skip_if_not(browserRuntimeAvailable(),
    "needs the installed package and a Chrome binary (CHROMOTE_CHROME or on PATH); see docs/qa/workflow-execution/kit/kit-034/")
}
