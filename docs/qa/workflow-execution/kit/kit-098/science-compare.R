# KIT-098 (plan section 5): the cause for `trx_*_science` in kitOperationIdentityProducts, shown on the two recordings. Their
# upstream states (design to RUVs) come from the same inputs and settings; for each such state the state's components and
# summary, the attached `trx.*:<sha12>` products and the `trx_*_science` products are compared between trx98a and trx98b.
# Usage: R_LIBS=<kit_098 library> Rscript docs/qa/workflow-execution/kit/kit-098/science-compare.R
repo <- normalizePath(".")
a <- readRDS(file.path(repo, "docs/qa/workflow-execution/kit/kit-098/oracles-trx98a_trx_screen.rds"))
b <- readRDS(file.path(repo, "docs/qa/workflow-execution/kit/kit-098/oracles-trx98b_trx_screen.rds"))
for (step in c("design", "expression_filter", "qc_finalize", "normalisation", "da")) {
    x <- a$populations[[step]]$fingerprint; y <- b$populations[[step]]$fingerprint
    science <- grep("_science$", names(x$products), value = TRUE)
    attached <- grep("^trx[.]", names(x$products), value = TRUE)
    same <- function(n) identical(x$products[[n]]$sha256, y$products[[n]]$sha256)
    cat(sprintf("%-18s components %s summary %s | attached trx.* products equal %d of %d | *_science equal %d of %d (%s)\n", step,
        identical(x$components, y$components), identical(x$summary_sha256, y$summary_sha256),
        sum(vapply(attached, function(n) n %in% names(y$products) && same(n), logical(1))), length(attached),
        sum(vapply(science, same, logical(1))), length(science), paste(science, collapse = ",")))
}
