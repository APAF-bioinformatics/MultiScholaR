# KIT-130 packaging spike, build time only: write the pip requirements for the bundled Python
# from the pins the installed MOFA2 carries for its own basilisk environment. run_mofa() checks
# the mofapy2 version it finds against this pin, so the bundle must match it.
#
# Usage: Rscript mofapy2-requirements.R <output directory>
# Writes requirements-pinned.txt (mofapy2 and MOFA2's dependency pins) and
# requirements-minimal.txt (the mofapy2 pin alone, dependencies left to pip).

out_dir <- commandArgs(trailingOnly = TRUE)[[1L]]
mofa_ns <- asNamespace("MOFA2")
mofapy2_version <- get0(".mofapy2_version", envir = mofa_ns, inherits = FALSE)
dependency_pins <- get0(".mofapy2_dependencies", envir = mofa_ns, inherits = FALSE)
if (is.null(mofapy2_version)) {
    stop("MOFA2 ", as.character(utils::packageVersion("MOFA2")), " carries no mofapy2 version pin")
}

mofapy2_pin <- paste0("mofapy2==", mofapy2_version)
# basilisk pins use conda syntax (name=version); pip needs name==version. Python itself is the
# standalone build, not a pip requirement.
pip_pins <- sub("^([A-Za-z0-9_.-]+)=+(.+)$", "\\1==\\2", as.character(dependency_pins))
pip_pins <- pip_pins[!grepl("^python==", pip_pins)]

cat(sprintf("MOFA2 %s pins: %s\n", as.character(utils::packageVersion("MOFA2")),
    paste(c(mofapy2_pin, as.character(dependency_pins)), collapse = ", ")))
writeLines(c(mofapy2_pin, pip_pins), file.path(out_dir, "requirements-pinned.txt"))
writeLines(mofapy2_pin, file.path(out_dir, "requirements-minimal.txt"))
