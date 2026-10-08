# CORE grouping-memory diagnostic (2026-10-08, user-approved; not a counted run, no performance claim).
# Runs the P14 arm F grouping shape on one input in a fresh process and reports peak resident memory.
# Usage: R_LIBS=<control library> Rscript grouping-memory-probe.R <variant> <input parquet> <output dir>
args <- commandArgs(TRUE); variant <- args[[1L]]; input <- args[[2L]]; out <- args[[3L]]
hwm <- function() as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE))) * 1024
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
DBI::dbExecute(con, "SET threads = 8"); DBI::dbExecute(con, "SET preserve_insertion_order = true")
DBI::dbExecute(con, "SET memory_limit = '12GiB'"); DBI::dbExecute(con, "SET temp_directory = ''")
DBI::dbExecute(con, sprintf("CREATE TEMP VIEW src AS SELECT Run, \"Protein.Group\", \"Stripped.Sequence\", \"Modified.Sequence\", \"Precursor.Quantity\" AS raw, \"Precursor.Normalised\" AS norm, \".multischolar_order\" AS o FROM read_parquet('%s')", input))
integer_keys <- variant %in% c("int_keys", "lower_bound", "int_streaming")
if (integer_keys) {
  for (k in c("Run", "Protein.Group", "Stripped.Sequence"))
    DBI::dbExecute(con, sprintf("CREATE TEMP TABLE \"d_%s\" AS SELECT v, (row_number() OVER (ORDER BY v))::INTEGER AS c FROM (SELECT DISTINCT \"%s\" AS v FROM src)", k, k))
  DBI::dbExecute(con, "CREATE TEMP TABLE prepared AS SELECT a.c AS k1, b.c AS k2, d.c AS k3, s.\"Modified.Sequence\" AS m, s.raw, s.norm, s.o FROM src s JOIN \"d_Run\" a ON s.Run = a.v JOIN \"d_Protein.Group\" b ON s.\"Protein.Group\" = b.v JOIN \"d_Stripped.Sequence\" d ON s.\"Stripped.Sequence\" = d.v")
} else DBI::dbExecute(con, "CREATE TEMP TABLE prepared AS SELECT Run AS k1, \"Protein.Group\" AS k2, \"Stripped.Sequence\" AS k3, \"Modified.Sequence\" AS m, raw, norm, o FROM src")
ordered <- !variant %in% c("streaming_sums", "lower_bound", "int_streaming", "kahan")
distinct_m <- variant %in% c("count_distinct", "lower_bound")
sums <- if (ordered) "list_sum(list(raw ORDER BY o)) AS r, list_sum(list(norm ORDER BY o)) AS n" else
  if (variant == "kahan") "kahan_sum(raw) AS r, kahan_sum(norm) AS n" else "sum(raw) AS r, sum(norm) AS n"
forms <- if (distinct_m) "count(DISTINCT m) AS f" else "list_unique(list(m)) AS f"
query <- sprintf("SELECT k1, k2, k3, %s, %s, count(raw), count(norm), count(*) FROM prepared GROUP BY k1, k2, k3", sums, forms)
before <- hwm(); gc()
t <- system.time(DBI::dbExecute(con, sprintf("COPY (%s ORDER BY k1, k2, k3) TO '%s' (FORMAT PARQUET)", query, file.path(out, paste0(variant, ".parquet")))))[["elapsed"]]
groups <- DBI::dbGetQuery(con, sprintf("SELECT count(*) n FROM read_parquet('%s')", file.path(out, paste0(variant, ".parquet"))))$n
cat(sprintf("RESULT variant=%s rows_in=%s groups=%s hwm_before_query=%.0fMB hwm_after_query=%.0fMB query_delta=%.0fMB seconds=%.2f\n",
  variant, DBI::dbGetQuery(con, "SELECT count(*) FROM prepared")[[1]], groups, before / 2^20, hwm() / 2^20, (hwm() - before) / 2^20, t))
