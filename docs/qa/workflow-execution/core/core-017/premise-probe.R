# CORE-017 premise probe (2026-10-08, in-process; not a counted run, no performance claim).
# Whole-input duplicate check and admission grouping statistics versus the same queries per batch of 3.
# Usage: R_LIBS=<control library> Rscript premise-probe.R <global|per_batch> <input parquet>
args <- commandArgs(TRUE); variant <- args[[1L]]; input <- args[[2L]]
hwm <- function() as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE))) * 1024
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
for (s in c("SET threads = 8", "SET preserve_insertion_order = true", "SET memory_limit = '12GiB'", "SET temp_directory = ''")) DBI::dbExecute(con, s)
DBI::dbExecute(con, sprintf("CREATE TEMP VIEW qc_input AS SELECT * FROM read_parquet('%s')", input))
b <- function(f) paste(sprintf("coalesce(octet_length(encode(CAST(\"%s\" AS VARCHAR))),0)", f), collapse = " + ")
dup <- function(w) sprintf("SELECT 1 FROM qc_input WHERE %s GROUP BY Run, \"Protein.Group\", \"Precursor.Id\" HAVING count(*) > 1 LIMIT 1", w)
stats <- function(w) sprintf("SELECT Run AS sample, count(*) AS g, sum(n) AS r, max(n) AS mg, max(vb + kb) AS mgb, sum(vb + kb) AS sb FROM (SELECT Run, \"Protein.Group\", \"Stripped.Sequence\", count(*) AS n, sum(32 + %s) AS vb, max(256 + %s) AS kb FROM qc_input WHERE %s GROUP BY ALL) GROUP BY 1", b("Modified.Sequence"), b(c("Run", "Protein.Group", "Stripped.Sequence")), w)
gc(); before <- hwm()
t <- system.time({
  if (variant == "global") { d <- nrow(DBI::dbGetQuery(con, dup("TRUE"))); s <- DBI::dbGetQuery(con, stats("TRUE")) } else {
    samples <- DBI::dbGetQuery(con, "SELECT DISTINCT Run FROM qc_input ORDER BY 1")$Run
    parts <- split(samples, ceiling(seq_along(samples) / ceiling(length(samples) / 3)))
    d <- 0L; s <- NULL
    for (p in parts) { w <- sprintf("Run IN (%s)", paste(DBI::dbQuoteString(con, p), collapse = ","))
      d <- d + nrow(DBI::dbGetQuery(con, dup(w))); s <- rbind(s, DBI::dbGetQuery(con, stats(w))) }
  }
})[["elapsed"]]
s <- s[order(s$sample), ]
cat(sprintf("RESULT variant=%s delta=%.0fMB seconds=%.2f dup=%d samples=%d groups=%.0f state=%.0f digest=%s\n", variant, (hwm() - before) / 2^20, t, d, nrow(s), sum(s$g), sum(s$sb), digest::digest(unname(as.list(s)))))
