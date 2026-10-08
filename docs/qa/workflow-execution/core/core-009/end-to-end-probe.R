# CORE-009 end-to-end key-encoding check (2026-10-08, in-process; not a counted run, no performance claim).
# Unlike diagnostics/grouping-memory-probe.R, dictionary build, encoding join and label decode are inside the measured span.
# Usage: R_LIBS=<control library> Rscript end-to-end-probe.R <strings|int_inline> <input parquet> <output dir>
args <- commandArgs(TRUE); variant <- args[[1L]]; input <- args[[2L]]; out <- args[[3L]]
hwm <- function() as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE))) * 1024
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
for (s in c("SET threads = 8", "SET preserve_insertion_order = true", "SET memory_limit = '12GiB'", "SET temp_directory = ''")) DBI::dbExecute(con, s)
DBI::dbExecute(con, sprintf("CREATE TEMP VIEW src AS SELECT Run, \"Protein.Group\", \"Stripped.Sequence\", \"Modified.Sequence\" AS m, \"Precursor.Quantity\" AS raw, \"Precursor.Normalised\" AS norm, \".multischolar_order\" AS o FROM read_parquet('%s')", input))
gc(); before <- hwm()
t <- system.time({
  if (variant == "int_inline") {
    for (k in c("Run", "Protein.Group", "Stripped.Sequence"))
      DBI::dbExecute(con, sprintf("CREATE TEMP TABLE \"d_%s\" AS SELECT v, (row_number() OVER (ORDER BY v))::INTEGER AS c FROM (SELECT DISTINCT \"%s\" AS v FROM src)", k, k))
    DBI::dbExecute(con, "CREATE TEMP VIEW prepared AS SELECT a.c AS k1, b.c AS k2, d.c AS k3, s.m, s.raw, s.norm, s.o FROM src s JOIN \"d_Run\" a ON s.Run = a.v JOIN \"d_Protein.Group\" b ON s.\"Protein.Group\" = b.v JOIN \"d_Stripped.Sequence\" d ON s.\"Stripped.Sequence\" = d.v")
    agg <- "SELECT k1, k2, k3, list_sum(list(raw ORDER BY o)) AS r, list_sum(list(norm ORDER BY o)) AS n, list_unique(list(m)) AS f, count(raw), count(norm), count(*) c FROM prepared GROUP BY k1, k2, k3"
    q <- sprintf("SELECT a.v AS Run, b.v AS pg, d.v AS ss, g.* EXCLUDE (k1, k2, k3) FROM (%s) g JOIN \"d_Run\" a ON g.k1 = a.c JOIN \"d_Protein.Group\" b ON g.k2 = b.c JOIN \"d_Stripped.Sequence\" d ON g.k3 = d.c ORDER BY g.k1, g.k2, g.k3", agg)
  } else {
    q <- "SELECT Run, \"Protein.Group\", \"Stripped.Sequence\", list_sum(list(raw ORDER BY o)) AS r, list_sum(list(norm ORDER BY o)) AS n, list_unique(list(m)) AS f, count(raw), count(norm), count(*) c FROM src GROUP BY ALL ORDER BY 1, 2, 3"
  }
  DBI::dbExecute(con, sprintf("COPY (%s) TO '%s' (FORMAT PARQUET)", q, file.path(out, paste0(variant, ".parquet"))))
})[["elapsed"]]
cat(sprintf("RESULT variant=%s delta=%.0fMB seconds=%.2f\n", variant, (hwm() - before) / 2^20, t))
