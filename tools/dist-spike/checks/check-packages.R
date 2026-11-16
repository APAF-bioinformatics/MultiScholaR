# KIT-130 piece 2: limpa, RUVIIIC and MOFA2 load from the bundled library; one call each on a
# toy input.
#
# Usage: <bundle>/.../Rscript --vanilla check-packages.R <bundle directory>

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "common.R"))
bundle <- expectedBundle()

for (pkg in c("limpa", "RUVIIIC", "MOFA2")) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    cat(sprintf("%s %s from %s\n", pkg, as.character(utils::packageVersion(pkg)), find.package(pkg)))
}

set.seed(130L)

# limpa: fit the detection probability curve on a simulated peptide set.
peptides <- limpa::simProteinDataSet(n.peptides = 100, n.groups = 2, samples.per.group = 3)
dpc_fit <- limpa::dpc(peptides)
cat(sprintf("limpa::dpc: beta0 = %.3f, beta1 = %.3f\n", dpc_fit$dpc[[1L]], dpc_fit$dpc[[2L]]))
stopifnot(length(dpc_fit$dpc) == 2L, all(is.finite(dpc_fit$dpc)))

# RUVIIIC: remove two factors from 12 samples in 4 replicate sets, with missing values
# (the compiled C++ path).
n_samples <- 12L
n_features <- 60L
abundance <- matrix(
    stats::rnorm(n_samples * n_features, mean = 20),
    n_samples,
    n_features,
    dimnames = list(paste0("s", seq_len(n_samples)), paste0("p", seq_len(n_features)))
)
abundance[sample(length(abundance), 40L)] <- NA
controls <- colnames(abundance)[1:20]
abundance[, controls] <- stats::rnorm(n_samples * length(controls), mean = 20)
replicate_design <- stats::model.matrix(~ 0 + factor(rep(seq_len(4L), each = 3L)))
rownames(replicate_design) <- rownames(abundance)
corrected <- RUVIIIC::RUVIII_C(
    k = 2L,
    Y = abundance,
    M = replicate_design,
    toCorrect = colnames(abundance),
    controls = controls,
    progress = FALSE
)
cat(sprintf("RUVIIIC::RUVIII_C: %d x %d corrected, %d finite\n", nrow(corrected), ncol(corrected), sum(is.finite(corrected))))
stopifnot(identical(dim(corrected), dim(abundance)), sum(is.finite(corrected)) > 0L)

# MOFA2: build and prepare a model object from the package's own example data (the R side and
# its HDF5 stack; training through Python is piece 3).
example_env <- new.env()
load(system.file("extdata", "test_data.RData", package = "MOFA2"), envir = example_env)
mofa_object <- MOFA2::prepare_mofa(MOFA2::create_mofa(example_env$dt))
cat(sprintf("MOFA2::prepare_mofa: status '%s', %d views, %d groups\n", mofa_object@status,
    length(MOFA2::views_names(mofa_object)), length(MOFA2::groups_names(mofa_object))))
stopifnot(methods::is(mofa_object, "MOFA"), mofa_object@status == "untrained")
hdf5_file <- tempfile(fileext = ".h5")
rhdf5::h5createFile(hdf5_file)
rhdf5::h5write(matrix(as.numeric(1:6), 2L), hdf5_file, "toy")
stopifnot(identical(dim(rhdf5::h5read(hdf5_file, "toy")), c(2L, 3L)))

# Every namespace these calls pulled in came from the bundle.
namespace_paths <- vapply(loadedNamespaces(), \(ns) getNamespaceInfo(ns, "path"), character(1L))
outside <- namespace_paths[!insideBundle(namespace_paths, bundle)]
cat(sprintf("namespaces loaded: %d, outside the bundle: %d\n", length(namespace_paths), length(outside)))
if (length(outside) > 0L) {
    print(outside)
    stop("namespaces were loaded from outside the bundle")
}

checkPassed("limpa, RUVIIIC, MOFA2")
