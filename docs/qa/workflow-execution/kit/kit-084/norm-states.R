# KIT-084: the log2 operation pairs' population and oracle, read from a record's stored ITSD and log2 states. The records
# bind one state per journey step, and the normalisation step's is the RUV state, so these two were never bound.
# Usage: KIT_INSTALLED library on R_LIBS; Rscript docs/qa/workflow-execution/kit/kit-084/norm-states.R <oracle rds> <prefix> <output rds>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE)
repo <- normalizePath(".")
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
recorded <- readRDS(args[[1L]])
prefix <- args[[2L]]
states <- paste0(prefix, c("_itsd_norm", "_log2"))
head <- harness$kitReadHead(recorded$store, api)
stopifnot(identical(head$head_sha256, recorded$head_sha256), all(states %in% names(head$states)))
fingerprints <- harness$kitFingerprintStore(recorded$store, api)
stopifnot(identical(fingerprints$head_sha256, recorded$head_sha256))
saveRDS(lapply(stats::setNames(states, c("itsd", "log2")), function(state) list(state = state, entity_id = head$states[[state]],
    rows = head$rows[[state]], fingerprint = fingerprints$states[[state]])), args[[3L]])
cat("written", args[[3L]], paste(states, unlist(head$rows[states])), "\n")
