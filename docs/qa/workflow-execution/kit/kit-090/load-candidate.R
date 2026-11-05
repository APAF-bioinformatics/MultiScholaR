# load-candidate.R: source() this from your own R scripts. It loads the installed package and puts the three functions from
# func_general_nipals.R (in this folder) into the package namespace in place of the installed ones, so the package and the
# tests call yours. It never shows the installed versions. Usage in a script: source("load-candidate.R"); load_candidate()
suppressMessages(library(MultiScholaR))
load_candidate <- function(file = file.path(getwd(), "func_general_nipals.R")) {
  ns <- asNamespace("MultiScholaR")
  candidate <- new.env(parent = ns)
  sys.source(file, envir = candidate)
  for (name in c("nipalsComponents", "imputeNipalsMatrix", "pcaWithMissing")) {
    if (!exists(name, envir = candidate, inherits = FALSE)) stop("func_general_nipals.R does not define ", name, call. = FALSE)
    fn <- get(name, envir = candidate)
    environment(fn) <- ns
    unlockBinding(name, ns); assign(name, fn, envir = ns); lockBinding(name, ns)
  }
  invisible(TRUE)
}
