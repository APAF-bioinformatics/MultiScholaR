# CORE-012 determinism check (in-process, not counted): does DuckDB kahan_sum depend on run or thread count?
# Usage: R_LIBS=<control library> Rscript determinism-check.R <input parquet> <result json>
args <- commandArgs(TRUE)
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
DBI::dbExecute(con, "SET preserve_insertion_order = true")
from <- sprintf("FROM read_parquet(%s) GROUP BY Run, \"Protein.Group\", \"Stripped.Sequence\" ORDER BY Run, \"Protein.Group\", \"Stripped.Sequence\"",
    DBI::dbQuoteString(con, args[[1L]]))
kahan <- function(threads) {
    DBI::dbExecute(con, sprintf("SET threads = %d", threads))
    DBI::dbGetQuery(con, paste("SELECT kahan_sum(\"Precursor.Quantity\") AS r, kahan_sum(\"Precursor.Normalised\") AS n", from))
}
DBI::dbExecute(con, "SET threads = 8")
lists <- DBI::dbGetQuery(con, paste("SELECT list(\"Precursor.Quantity\" ORDER BY \".multischolar_order\") AS lr,",
    "list(\"Precursor.Normalised\" ORDER BY \".multischolar_order\") AS ln", from))
# The C kernel's semantics: skip missing, long-double accumulation in source order.
ordered <- function(l) vapply(l, function(v) { v <- v[!is.na(v)]; if (length(v)) sum(v) else NA_real_ }, numeric(1))
ord <- list(r = ordered(lists$lr), n = ordered(lists$ln))
differ <- function(a, b) sum(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b))
eight <- lapply(1:4, function(i) kahan(8L))
one <- kahan(1L)
result <- list(input = args[[1L]], groups = nrow(one), duckdb = as.character(packageVersion("duckdb")))
for (col in c("r", "n")) {
    name <- if (col == "r") "raw" else "normalised"
    diff_one_eight <- xor(is.na(one[[col]]), is.na(eight[[1]][[col]])) |
        (!is.na(one[[col]]) & one[[col]] != eight[[1]][[col]])
    diff_eight_ordered <- !is.na(ord[[col]]) & eight[[1]][[col]] != ord[[col]]
    result[[name]] <- list(
        repeats_at_8_threads_differing = vapply(2:4, \(i) differ(eight[[i]][[col]], eight[[1]][[col]]), numeric(1)),
        one_vs_eight_threads_differing = sum(diff_one_eight),
        one_thread_vs_ordered_differing = differ(one[[col]], ord[[col]]),
        eight_threads_vs_ordered_differing = differ(eight[[1]][[col]], ord[[col]]),
        same_groups_one_eight_and_eight_ordered = identical(which(diff_one_eight), which(diff_eight_ordered)))
}
jsonlite::write_json(result, args[[2L]], auto_unbox = TRUE, pretty = TRUE)
for (name in c("raw", "normalised")) with(result[[name]], cat(sprintf(
    "%s: repeats@8 differing %s; 1 vs 8 threads %d; 1 thread vs ordered %d; 8 threads vs ordered %d; same groups %s\n",
    name, paste(repeats_at_8_threads_differing, collapse = "/"), one_vs_eight_threads_differing,
    one_thread_vs_ordered_differing, eight_threads_vs_ordered_differing, same_groups_one_eight_and_eight_ordered)))
