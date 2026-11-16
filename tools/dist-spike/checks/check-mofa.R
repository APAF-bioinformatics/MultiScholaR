# KIT-130 piece 3: the standalone CPython in the folder, with mofapy2, trains a toy model through
# MOFA2::run_mofa(use_basilisk = FALSE); nothing is downloaded at run time.
#
# Usage: <bundle>/.../Rscript --vanilla check-mofa.R <bundle directory> <online|offline>
# "offline" asserts that this process cannot reach the network (the caller has cut it), which is
# what makes a successful run proof that nothing was fetched.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "common.R"))
bundle <- expectedBundle()
expect_offline <- identical(commandArgs(trailingOnly = TRUE)[2L], "offline")

reachable <- canReachNetwork()
cat(sprintf("network reachable from R: %s (expected offline: %s)\n", reachable, expect_offline))
cat(sprintf("python on the PATH: '%s'\n", paste(Sys.which(c("python", "python3")), collapse = "' '")))
print(Sys.getenv(c("RETICULATE_PYTHON", "PYTHONNOUSERSITE", "PYTHONPATH", "PYTHONHOME")))
if (expect_offline) {
    stopifnot(
        "the network is cut for this process" = !reachable,
        "no other Python on the PATH" = length(foreignOnPath(c("python", "python3"), bundle)) == 0L
    )
}

suppressPackageStartupMessages(library(MOFA2))
example_env <- new.env()
load(system.file("extdata", "test_data.RData", package = "MOFA2"), envir = example_env)
mofa_object <- MOFA2::create_mofa(example_env$dt)
training_options <- MOFA2::get_default_training_options(mofa_object)
training_options$seed <- 130L
training_options$maxiter <- 200L
mofa_object <- MOFA2::prepare_mofa(mofa_object, training_options = training_options)

model_file <- file.path(tempdir(), "kit-130-toy-model.hdf5")
messages <- character()
trained <- withCallingHandlers(
    MOFA2::run_mofa(mofa_object, outfile = model_file, use_basilisk = FALSE),
    message = \(m) messages <<- c(messages, conditionMessage(m))
)

python <- reticulate::py_config()
print(python)
module_file <- \(name) reticulate::import(name)[["__file__"]]
python_paths <- c(
    python = python$python,
    libpython = python$libpython,
    mofapy2 = module_file("mofapy2"),
    numpy = module_file("numpy"),
    h5py = module_file("h5py")
)
print(python_paths)
python_reachable <- tryCatch(
    {
        reticulate::py_run_string("import socket; socket.create_connection(('pypi.org', 80), 10).close()")
        TRUE
    },
    error = \(e) FALSE
)
cat(sprintf("network reachable from Python: %s\n", python_reachable))

factors <- MOFA2::get_factors(trained)
cat(sprintf(
    "trained: status '%s', %d factors, %d groups, model file %d bytes, mofapy2 %s\n",
    trained@status, ncol(factors[[1L]]), length(factors), file.size(model_file),
    reticulate::import("importlib.metadata")$version("mofapy2")
))

stopifnot(
    "the model is trained" = trained@status == "trained",
    "at least one factor was learnt" = ncol(factors[[1L]]) >= 1L,
    "factor values are finite" = all(is.finite(unlist(factors))),
    "Python, its library and its modules are inside the bundle" = all(insideBundle(python_paths, bundle)),
    "basilisk was not used" = !any(grepl("using basilisk", messages, fixed = TRUE))
)
if (expect_offline) {
    stopifnot("the network is cut for Python too" = !python_reachable)
}

checkPassed(if (expect_offline) "Python for MOFA2 (offline)" else "Python for MOFA2 (online)")
