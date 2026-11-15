# KIT-118: two in-process oracle records (before and after a mechanism) compared step by step: the state's class and
# rows, each component and support entry, each product by name, and the summary hash. Usage (anywhere):
#   Rscript compare-records.R <before.rds> <after.rds> <out.tsv>
args <- commandArgs(TRUE)
`%||%` <- function(p, q) if (is.null(p)) q else p
a <- readRDS(args[[1L]]); b <- readRDS(args[[2L]])
digest_of <- function(x) digest::digest(x, algo = "sha256")
rows <- list()
add <- function(step, field, same, note = "") rows[[length(rows) + 1L]] <<- data.frame(step = step, field = field,
    identical = same, note = note)
for (step in union(names(a$populations), names(b$populations))) {
    x <- a$populations[[step]]$fingerprint; y <- b$populations[[step]]$fingerprint
    if (is.null(x) || is.null(y)) { add(step, "present", is.null(x) && is.null(y)); next }
    add(step, "class", identical(x$class, y$class))
    add(step, "rows", identical(x$rows, y$rows), paste(x$rows %||% NA, "->", y$rows %||% NA))
    for (part in c("components", "support", "products")) for (name in union(names(x[[part]]), names(y[[part]]))) {
        u <- x[[part]][[name]]; v <- y[[part]][[name]]
        note <- if (is.list(u) && is.list(v) && !is.null(u$rows) && !is.null(v$rows)) paste(u$rows, "->", v$rows) else ""
        add(step, paste0(part, "/", name), identical(digest_of(u), digest_of(v)), note)
    }
    add(step, "summary_sha256", identical(x$summary_sha256, y$summary_sha256))
}
`%||%` <- function(p, q) if (is.null(p)) q else p
out <- do.call(rbind, rows)
utils::write.table(out, args[[3L]], sep = "\t", quote = FALSE, row.names = FALSE)
cat("STEPS", length(unique(out$step)), "COMPARED", nrow(out), "DIFFERENT", sum(!out$identical), "\n")
print(out[!out$identical, ], row.names = FALSE)
