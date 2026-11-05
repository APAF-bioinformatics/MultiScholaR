# KIT-092 (protocol v1.7): re-fingerprint one in-process oracle record from its store under the content form.
# Usage, from the repository root, with the recording build's library first on R_LIBS:
#   Rscript docs/qa/workflow-execution/kit/kit-092/refingerprint.R <oracles-*.rds> <v1.6 workflow_kit.R>
# The store is only read. The record's head must still be the store's head (else the oracle is listed, not guessed).
# Each population is recomputed twice: under the v1.6 harness, where it must equal the recorded fingerprint exactly (the
# recomputation is faithful), and under the current harness, where class, rows, components, support and summary must be
# unchanged and only products may differ. The record is then rewritten with the v1.7 fingerprints (the old file stays in
# git history) and every changed product is printed.
args <- commandArgs(TRUE)
oracle_file <- args[[1L]]
old_harness_file <- args[[2L]]
repo <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
recorded <- readRDS(oracle_file)
label <- basename(oracle_file)
library_dir <- dirname(recorded$library)
if (!dir.exists(recorded$store)) {
    cat("ORACLE", label, "UNAVAILABLE store missing:", recorded$store, "\n")
    quit(status = 0L)
}
stopifnot(identical(normalizePath(find.package("MultiScholaR")), normalizePath(recorded$library)))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
current <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), current)
old <- new.env(parent = globalenv()); sys.source(old_harness_file, old)
head <- current$kitReadHead(recorded$store, api)
if (!identical(head$head_sha256, recorded$head_sha256)) {
    cat("ORACLE", label, "UNAVAILABLE store head moved on:", head$head_sha256, "recorded", recorded$head_sha256, "\n")
    quit(status = 0L)
}
fingerprints <- function(harness) {
    store <- harness$kitFingerprintStore(recorded$store, api)
    stopifnot(identical(store$head_sha256, recorded$head_sha256))
    attached <- current$kitStepProducts(recorded$populations, store$products, rendered = FALSE)
    lapply(stats::setNames(names(recorded$populations), names(recorded$populations)), function(step) {
        state <- store$states[[recorded$populations[[step]]$state]]
        if (!is.null(state)) state$products <- c(state$products, attached[[step]])
        state
    })
}
started <- Sys.time()
before <- fingerprints(old)
after <- fingerprints(current)
from_sha256 <- digest::digest(file = oracle_file, algo = "sha256")
changed_products <- 0L
for (step in names(recorded$populations)) {
    was <- recorded$populations[[step]]$fingerprint
    if (!identical(before[[step]], was)) stop(label, " ", step, ": the v1.6 recomputation differs from the record")
    now <- after[[step]]
    science <- setdiff(union(names(was), names(now)), "products")
    if (!identical(was[science], now[science])) stop(label, " ", step, ": a state field changed under the content form")
    old_names <- names(was$products) %||% character()
    new_names <- names(now$products) %||% character()
    same <- vapply(old_names, function(name) name %in% new_names && identical(was$products[[name]], now$products[[name]]), logical(1))
    for (name in old_names[!same]) cat(sprintf("  %-14s %-12s was %s\n", label, step, name))
    for (name in setdiff(new_names, old_names[same])) cat(sprintf("  %-14s %-12s now %s\n", label, step, name))
    changed_products <- changed_products + sum(!same)
    recorded$populations[[step]]$fingerprint <- now
}
recorded$canonical_form <- list(protocol = "1.7", ticket = "KIT-092", refingerprinted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    from_sha256 = from_sha256, library = library_dir,
    harness = "tools/profiling/workflow_kit.R (kitFingerprintStore, kitStepProducts; nested tables by pages)",
    note = "every state field reproduced unchanged; products re-fingerprinted; the v1.6 file is in git history")
saveRDS(recorded, oracle_file, compress = FALSE, version = 3)
cat("ORACLE", label, "store", basename(recorded$store), "head", substr(recorded$head_sha256, 1, 12), "| steps", length(recorded$populations),
    "| v1.6 recomputation equal to the record | state fields unchanged | products changed", changed_products,
    "| from", substr(from_sha256, 1, 12), "to", substr(digest::digest(file = oracle_file, algo = "sha256"), 1, 12),
    "|", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "s\n")
