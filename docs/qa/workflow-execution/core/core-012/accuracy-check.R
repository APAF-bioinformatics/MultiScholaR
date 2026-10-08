# CORE-012 step 1 (in-process, not counted): ordered long-double sum vs DuckDB kahan_sum vs exact rational sum.
# Usage: R_LIBS=<control library>:<reference library> Rscript accuracy-check.R <input parquet> <result json>
args <- commandArgs(TRUE)
suppressPackageStartupMessages({ library(gmp); library(Rmpfr) })
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
DBI::dbExecute(con, "SET threads = 8"); DBI::dbExecute(con, "SET preserve_insertion_order = true")
groups <- DBI::dbGetQuery(con, sprintf(paste(
    "SELECT list(\"Precursor.Quantity\" ORDER BY \".multischolar_order\") AS lr,",
    "list(\"Precursor.Normalised\" ORDER BY \".multischolar_order\") AS ln,",
    "kahan_sum(\"Precursor.Quantity\") AS kr, kahan_sum(\"Precursor.Normalised\") AS kn",
    "FROM read_parquet('%s') GROUP BY Run, \"Protein.Group\", \"Stripped.Sequence\"",
    "ORDER BY Run, \"Protein.Group\", \"Stripped.Sequence\""), args[[1L]]))
# The kernel's semantics: skip missing, long-double accumulation in source order, NA if nothing observed.
ordered <- function(lists) vapply(lists, function(v) { v <- v[!is.na(v)]; if (length(v)) sum(v) else NA_real_ }, numeric(1))
exactBy <- function(lists) {
    clean <- lapply(lists, \(v) v[!is.na(v)])
    total <- as.bigq(rep(0, length(clean)))
    for (j in seq_len(max(lengths(clean)))) total <- total + as.bigq(vapply(clean, \(v) if (length(v) >= j) v[[j]] else 0, numeric(1)))
    list(total = total, observed = lengths(clean) > 0L)
}
deviation <- function(approx, exact) {
    out <- rep(0, length(approx)); live <- exact$observed & !is.na(approx)
    e <- mpfr(exact$total[live], precBits = 256L); a <- mpfr(as.bigq(approx[live]), precBits = 256L)
    nz <- as.numeric(abs(e)) > 0
    out[live][nz] <- as.numeric(abs(a[nz] - e[nz]) / abs(e[nz]))
    out[live][!nz] <- as.numeric(abs(a[!nz]))
    out
}
result <- list(input = args[[1L]], groups = nrow(groups),
    group_sizes = as.list(table(pmin(lengths(groups$lr), 6L))))
for (column in c("raw", "normalised")) {
    lists <- if (column == "raw") groups$lr else groups$ln
    kahan <- if (column == "raw") groups$kr else groups$kn
    ord <- ordered(lists)
    exact <- exactBy(lists)
    d_ord <- deviation(ord, exact); d_kahan <- deviation(kahan, exact)
    same <- (is.na(ord) & is.na(kahan)) | (!is.na(ord) & !is.na(kahan) & ord == kahan)
    finite <- ord[is.finite(ord)]; finite_k <- kahan[is.finite(kahan)]
    cut <- lapply(c(0.5, 1, 2.5, 5, 10, 25, 50, 75, 90, 95, 99), function(p) {
        t_o <- unname(ceiling(stats::quantile(finite, p / 100)))
        t_k <- unname(ceiling(stats::quantile(finite_k, p / 100)))
        list(percentile = p, threshold_ordered = t_o, threshold_kahan = t_k,
            identical_threshold = identical(t_o, t_k),
            identical_retention = identical(ord >= t_o, kahan >= t_k))
    })
    result[[column]] <- list(bit_identical_groups = sum(same), differing_groups = sum(!same),
        missing_ordered = sum(is.na(ord)), missing_kahan = sum(is.na(kahan)),
        max_relative_deviation_ordered = max(d_ord), max_relative_deviation_kahan = max(d_kahan),
        within_tolerance_ordered = all(d_ord <= 2^-50), within_tolerance_kahan = all(d_kahan <= 2^-50),
        cutoffs = cut)
}
jsonlite::write_json(result, args[[2L]], auto_unbox = TRUE, pretty = TRUE, digits = NA)
for (column in c("raw", "normalised")) with(result[[column]], cat(sprintf(
    "%s: identical %d / %d groups; max rel dev ordered %.3g kahan %.3g; within 2^-50: ordered %s kahan %s; cutoffs identical %d/%d\n",
    column, bit_identical_groups, result$groups, max_relative_deviation_ordered, max_relative_deviation_kahan,
    within_tolerance_ordered, within_tolerance_kahan,
    sum(vapply(cutoffs, \(x) x$identical_threshold && x$identical_retention, logical(1))), length(cutoffs))))
