.libPaths(c("/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/core/browser-smoke.UgUYyi/library", .libPaths()))
suppressMessages(library(fastduckdb))
dir <- commandArgs(TRUE)[[1L]]; path <- file.path(dir, "stage.parquet"); n <- 4000000L
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE)); DBI::dbExecute(con, "SET threads = 8"); DBI::dbExecute(con, "SET preserve_insertion_order = true")
DBI::dbExecute(con, sprintf("COPY (SELECT printf('run_%%04d', i // 40000 + 1) AS Run, printf('PG_%%06d', (i // 4) %% 10000) AS \"Protein.Group\",
  printf('PEPTIDE%%07dK', i %% 400000) AS \"Stripped.Sequence\", (i %% 9)::INTEGER AS peptidoform_count, (i %% 7)::INTEGER AS identification_peptide_count,
  (i %% 11)::INTEGER AS identification_peptidoform_count, (i %% 5)::INTEGER AS peptides_for_protein_count, (i %% 13)::INTEGER AS peptidoforms_for_protein_count,
  printf('P%%05d', (i // 4) %% 10000) AS \"Protein.Ids\", random() * 1e7 AS raw_observed, random() * 1e7 AS normalised_observed, (i %% 3 + 1)::INTEGER AS source_count,
  row_number() OVER (ORDER BY i) AS \".dia_output_order\" FROM range(%d) t(i) ORDER BY i) TO '%s' (FORMAT PARQUET, COMPRESSION SNAPPY, ROW_GROUP_SIZE 10240)", n, path))
rg <- DBI::dbGetQuery(con, sprintf("SELECT count(DISTINCT row_group_id) g FROM parquet_metadata('%s')", path))$g
cat(sprintf("file %.0f MB, %d rows, %d row groups\n", file.size(path) / 1e6, n, rg))
current <- function() { fastduckdb::clear_parquet_cache(); rows <- 0
  for (start in seq.int(1L, n, by = 10000L)) { end <- min(n, start + 9999L)
    part <- fastduckdb::read_parquet_fast(path, n_threads = 8L, lazy_strings = FALSE, filter = list(.dia_output_order = c(start, end)))
    rows <- rows + nrow(part) }
  fastduckdb::clear_parquet_cache(); rows }
aligned <- function() { fastduckdb::clear_parquet_cache(); rows <- 0; carry <- NULL
  for (start in seq.int(1L, n, by = 10240L)) { end <- min(n, start + 10239L)
    part <- fastduckdb::read_parquet_fast(path, n_threads = 8L, lazy_strings = FALSE, filter = list(.dia_output_order = c(start, end)))
    if (!is.null(carry)) part <- rbind(carry, part)
    while (nrow(part) >= 10000L) { batch <- part[seq_len(10000L), , drop = FALSE]; rows <- rows + nrow(batch); part <- part[-seq_len(10000L), , drop = FALSE] }
    carry <- part }
  if (!is.null(carry)) rows <- rows + nrow(carry)
  fastduckdb::clear_parquet_cache(); rows }
readonly <- function() { fastduckdb::clear_parquet_cache(); rows <- 0
  for (start in seq.int(1L, n, by = 10240L)) { end <- min(n, start + 10239L)
    rows <- rows + nrow(fastduckdb::read_parquet_fast(path, n_threads = 8L, lazy_strings = FALSE,
      filter = list(.dia_output_order = c(start, end)))) }
  fastduckdb::clear_parquet_cache(); rows }
for (rep in 1:3) for (arm in c("current", "readonly")) {
  t <- system.time(r <- get(arm)())[["elapsed"]]; stopifnot(r == n)
  cat(sprintf("rep %d %-8s %.2f s\n", rep, arm, t)) }
