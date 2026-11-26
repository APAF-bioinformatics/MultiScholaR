# Entry point for the run-tests `scripts` mechanism (Rscript from the package root): E16, then its verdict as the exit status.
options(warnPartialMatchDollar = TRUE) # the guard of the scripts it runs; each Rscript sets its own
here <- normalizePath("docs/qa/workflow-execution/kit/kit-188/e16")
out <- file.path(getwd(), "k188-smoke-e16")
quit(save = "no", status = system2("bash", c(shQuote(file.path(here, "run-e16.sh")), shQuote(here), shQuote(out))))
