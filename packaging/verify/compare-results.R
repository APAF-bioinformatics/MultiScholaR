# KIT-134 revision 2.0.4: compare two saved results under the fixed comparison rule.
#
# Usage: Rscript --vanilla packaging/verify/compare-results.R <a.rds> <b.rds> <out.tsv> [name ...]
#
# <a.rds> and <b.rds> hold named lists, or a KIT-040 reference (a list with $products). The named elements compared are
# those given, else every name both hold (a name in only one is reported). The rule, fixed before any run was read:
#   doubles agree elementwise when |a - b| <= 1e-6 * max(|a|, |b|) or |a - b| <= 1e-9; the positions and kinds of NA,
#   NaN and Inf exactly;
#   everything else exactly: integers, logicals, characters, factors and their levels, names, dimensions, dimnames,
#   row names, classes, and list structure.
# One row per leaf: path, kind, compared elements, max relative and absolute difference (doubles), elements outside
# the rule or differing, verdict (agree, outside, differ). Prints a summary; the exit status is 0 (a measurement).

REL <- 1e-6
ABS <- 1e-9

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 3L) stop("usage: Rscript --vanilla compare-results.R <a.rds> <b.rds> <out.tsv> [name ...]")
unwrap <- function(x) if (is.list(x) && !is.null(x$products) && is.list(x$products)) x$products else x
a_all <- unwrap(readRDS(arguments[[1L]]))
b_all <- unwrap(readRDS(arguments[[2L]]))
wanted <- if (length(arguments) > 3L) arguments[-(1:3)] else intersect(names(a_all), names(b_all))

rows <- list()
add <- function(path, kind, n = NA_integer_, max_rel = NA_real_, max_abs = NA_real_, outside = NA_integer_, verdict, note = "") {
    rows[[length(rows) + 1L]] <<- data.frame(path = path, kind = kind, n = n, max_rel = max_rel, max_abs = max_abs,
        outside = outside, verdict = verdict, note = note, stringsAsFactors = FALSE)
}

compare <- function(a, b, path) {
    if (!identical(class(a), class(b))) {
        return(add(path, "class", verdict = "differ", note = paste(paste(class(a), collapse = "/"), "vs", paste(class(b), collapse = "/"))))
    }
    for (attribute in c("names", "dim", "dimnames", "levels")) {
        if (!identical(attr(a, attribute), attr(b, attribute))) add(paste0(path, "@", attribute), "attribute", verdict = "differ")
    }
    if (is.data.frame(a) && !identical(attr(a, "row.names"), attr(b, "row.names"))) add(paste0(path, "@row.names"), "attribute", verdict = "differ")
    if (is.list(a)) {
        if (length(a) != length(b)) return(add(path, "list", verdict = "differ", note = sprintf("length %d vs %d", length(a), length(b))))
        labels <- if (is.null(names(a))) paste0("[[", seq_along(a), "]]") else paste0("$", names(a))
        for (i in seq_along(a)) compare(a[[i]], b[[i]], paste0(path, labels[[i]]))
        return(invisible())
    }
    if (length(a) != length(b)) return(add(path, typeof(a), verdict = "differ", note = sprintf("length %d vs %d", length(a), length(b))))
    if (is.double(a) && !is.factor(a)) {
        a <- as.vector(a); b <- as.vector(b)
        kinds <- identical(is.na(a), is.na(b)) && identical(is.nan(a), is.nan(b)) &&
            identical(a[is.infinite(a)], b[is.infinite(b)]) && identical(is.infinite(a), is.infinite(b))
        finite <- is.finite(a) & is.finite(b)
        d <- abs(a[finite] - b[finite])
        m <- pmax(abs(a[finite]), abs(b[finite]))
        ok <- d <= REL * m | d <= ABS
        relative <- ifelse(m > 0, d / m, 0)
        return(add(path, "double", n = sum(finite), max_rel = if (length(d)) max(relative) else 0,
            max_abs = if (length(d)) max(d) else 0, outside = sum(!ok),
            verdict = if (kinds && all(ok)) "agree" else "outside", note = if (kinds) "" else "NA/NaN/Inf positions differ"))
    }
    differing <- if (is.atomic(a)) sum(!(a == b | (is.na(a) & is.na(b))), na.rm = TRUE) + sum(xor(is.na(a), is.na(b))) else NA_integer_
    add(path, typeof(a), n = length(a), outside = differing, verdict = if (identical(a, b)) "agree" else "differ")
}

for (name in wanted) {
    if (!name %in% names(a_all) || !name %in% names(b_all)) {
        add(name, "missing", verdict = "differ", note = sprintf("in a: %s; in b: %s", name %in% names(a_all), name %in% names(b_all)))
        next
    }
    compare(a_all[[name]], b_all[[name]], name)
}
result <- do.call(rbind, rows)
utils::write.table(result, arguments[[3L]], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
options(width = 200)
print(result[result$verdict != "agree" | result$kind == "double", ], row.names = FALSE)
cat(sprintf("\n%d leaves compared: %d agree, %d outside the rule, %d differ (discrete or structure)\n",
    nrow(result), sum(result$verdict == "agree"), sum(result$verdict == "outside"), sum(result$verdict == "differ")))
doubles <- result[result$kind == "double", ]
if (nrow(doubles)) cat(sprintf("largest relative difference over the doubles: %.3g (%s)\n",
    max(doubles$max_rel), doubles$path[which.max(doubles$max_rel)]))
