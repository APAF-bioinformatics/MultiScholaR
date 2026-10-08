# CORE-011 premise probe (2026-10-08, in-process; not a counted run, no performance claim).
# Can the P14 grouping shape spill under a DuckDB memory_limit below its working set?
# Usage: R_LIBS=<control library> Rscript spill-probe.R <ordered|plain> <memory_limit> <input parquet> <output dir>
args <- commandArgs(TRUE); shape <- args[[1L]]; limit <- args[[2L]]; input <- args[[3L]]; out <- args[[4L]]
hwm <- function() as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE))) * 1024
tmp <- file.path(out, paste0("spill-", shape, "-", limit)); unlink(tmp, recursive = TRUE); dir.create(tmp)
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
for (s in c("SET threads = 8", "SET preserve_insertion_order = true", sprintf("SET memory_limit = '%s'", limit),
  sprintf("SET temp_directory = '%s'", tmp), "SET max_temp_directory_size = '20GiB'")) DBI::dbExecute(con, s)
DBI::dbExecute(con, sprintf("CREATE TEMP VIEW src AS SELECT * FROM read_parquet('%s')", input))
sums <- if (shape == "ordered") paste("list_sum(list(\"Precursor.Quantity\" ORDER BY \".multischolar_order\")) AS r,",
    "list_sum(list(\"Precursor.Normalised\" ORDER BY \".multischolar_order\")) AS n, list_unique(list(\"Modified.Sequence\")) AS f") else
    "sum(\"Precursor.Quantity\") AS r, sum(\"Precursor.Normalised\") AS n"
q <- sprintf("SELECT Run, \"Protein.Group\", \"Stripped.Sequence\", %s, count(*) c FROM src GROUP BY ALL ORDER BY 1, 2, 3", sums)
gc(); before <- hwm()
res <- tryCatch({ t <- system.time(DBI::dbExecute(con, sprintf("COPY (%s) TO '%s' (FORMAT PARQUET)", q,
    file.path(out, paste0("out-", shape, "-", limit, ".parquet")))))[["elapsed"]]; sprintf("ok seconds=%.2f", t) },
  error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 160)))
cat(sprintf("RESULT shape=%s limit=%s %s rss_growth=%.0fMB\n", shape, limit, res, (hwm() - before) / 2^20))
