# content-compare.R <proof root> <journey A> <journey B>: every state and every product of two journey stores, compared
# by content. The harness's canonical form (tools/profiling/workflow_kit.R, unchanged) hashes a table nested inside a
# product by its file bytes; here, in this script's copy only, a nested table is hashed by its pages as a state's own
# components already are. A partitioned write gives the same content in different bytes. Synthetic corpora.
args <- commandArgs(TRUE); P <- args[[1L]]
suppressMessages(library(MultiScholaR)); api <- function(name) get(name, asNamespace("MultiScholaR"))
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
original <- h$kitCanonical
h$kitCanonical <- function(value, store = NULL, api = NULL, depth = 0L, memo = NULL) {
  if (is.list(value) && !is.data.frame(value) && is.list(value$payload) && identical(value$payload$codec, "parquet_duckdb") &&
      is.character(value$entity_id) && !is.null(store) && !is.null(api)) {
    key <- paste0("table:", value$entity_id)
    if (!is.null(memo) && !is.null(memo[[key]])) return(memo[[key]])
    out <- list(.ref = "table", content = h$kitFingerprintTable(store, value, api)$sha256)
    if (!is.null(memo)) memo[[key]] <- out
    return(out)
  }
  original(value, store, api, depth, memo)
}
environment(h$kitCanonical) <- h
print_store <- function(label) {
  sd <- h$kitFindStore(file.path(P, "kit/runtime", label, "candidate-1", "project"))
  fp <- h$kitFingerprintStore(sd, api)
  head <- h$kitReadHead(sd, api)
  list(states = fp$states, history = head$history,
    products = sort(vapply(fp$products, function(p) paste0(p$product_type %||% p$kind, ":", p$sha256), "")))
}
a <- print_store(args[[2L]]); b <- print_store(args[[3L]])
cat("history identical", identical(a$history, b$history), "states", length(a$history), "\n")
bad <- 0L
for (s in a$history) {
  x <- a$states[[s]]; y <- b$states[[s]]
  same <- c(class = identical(x$class, y$class), rows = identical(x$rows, y$rows), components = identical(x$components, y$components),
    support = identical(x$support, y$support), summary = identical(x$summary_sha256, y$summary_sha256), products = identical(x$products, y$products))
  if (!all(same)) bad <- bad + 1L
  cat(sprintf("STATE %-30s %s\n", s, paste(names(same), same, sep = "=", collapse = " ")))
}
cat("PRODUCTS", args[[2L]], length(a$products), args[[3L]], length(b$products), "identical multiset", identical(unname(a$products), unname(b$products)), "\n")
only_a <- setdiff(a$products, b$products); only_b <- setdiff(b$products, a$products)
if (length(only_a)) cat("  only in", args[[2L]], ":", substr(only_a, 1, 60), "\n")
if (length(only_b)) cat("  only in", args[[3L]], ":", substr(only_b, 1, 60), "\n")
# Rendered figure pages are reported, never judged: their bytes differ between members of one build (KIT-078, as the harness's
# member reproduction does with `^rendered_asset:`).
judged <- function(x) unname(x[!grepl("^rendered_asset:", x)])
cat("JUDGED PRODUCTS", length(judged(a$products)), length(judged(b$products)), "identical", identical(judged(a$products), judged(b$products)),
  "| rendered pages", sum(grepl("^rendered_asset:", a$products)), sum(grepl("^rendered_asset:", b$products)), "\n")
cat("RESULT", args[[2L]], "vs", args[[3L]], if (bad == 0L && identical(judged(a$products), judged(b$products))) "content identical" else "DIFFERENT", "\n")
