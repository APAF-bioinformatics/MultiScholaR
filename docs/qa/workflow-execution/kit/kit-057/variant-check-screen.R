# KIT-057: each production variant against the whole variant, in process on a copy of a DIA_LIMPA project (the
# accepted screen lifecycle member's, or the counted member's for the full corpus): the same source state, both
# computes, the outputs compared field by field (matrix, tables by page hash, config, args) and the process's peak
# RSS (VmHWM) read after each. Usage: [KIT053_INSTALLED=1] R_LIBS=<library> Rscript .../variant-check-screen.R <project copy> <step> <label>
if (identical(Sys.getenv("KIT053_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
ns <- asNamespace("MultiScholaR"); api <- function(n) get(n, ns)
for (f in ls(ns)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
a <- commandArgs(TRUE); project <- a[[1]]; step <- a[[2]]; label <- a[[3]]
# Memory is a process high-water mark, so each variant runs in its own process: `run:<variant>` computes one variant and
# saves its record beside the JSON; `compare` reads the records and writes the JSON. Without a fourth argument both
# variants run in this process (exactness only; the second variant's peak is not its own).
mode <- if (length(a) >= 4L) a[[4]] else "both"
record_path <- function(variant) file.path("docs/qa/workflow-execution/kit/kit-057", paste0("variant-check-", step, "-", label, "-", variant, ".rds"))
harness <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", harness)
store <- harness$kitFindStore(project); root <- dirname(store)
id <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(root, id$project_id, id$workflow_id), reopen = TRUE)
hwm <- function() as.numeric(strsplit(grep("VmHWM", readLines("/proc/self/status"), value = TRUE), "[[:space:]]+")[[1]][2]) / 2^20
table_sha <- function(reference) { sc <- newWorkflowTableReadScope(manager$.context$store, reference, reader = "fastduckdb"); on.exit(releaseWorkflowTableReadScope(sc))
  n <- sc$metadata$rows; hashes <- character(); if (n > 0L) for (off in seq.int(0L, n - 1L, by = 10000L)) { d <- readSavedWorkflowTable(manager$.context$store, reference, sc$metadata$columns, offset = off, limit = min(10000L, n - off), max_columns = max(64L, length(sc$metadata$columns)), scope = sc)$data; attr(d, "dia_source_order") <- NULL; attr(d, "row.names") <- .set_row_names(nrow(d)); hashes <- c(hashes, digest::digest(d, algo = "sha256")) }
  list(rows = n, columns = sc$metadata$columns, sha256 = digest::digest(paste(hashes, collapse = "\n"), algo = "sha256", serialize = FALSE)) }
spec <- switch(step,
  peptide_log2 = list(source = "precursor_rollup", variants = c("whole", "streamed"), compute = function(ref, variant) computeDIAReferencePeptideLog2Matrix(manager, ref, list(), variant = variant)),
  peptide_accession = list(source = "log2_matrix", variants = c("whole", "partitioned"), compute = function(ref, variant) {
    fasta_path <- list.files(project, pattern = "^aa_seq_tbl_final[.]RDS$", recursive = TRUE, full.names = TRUE)[[1]]
    sink(tempfile()); on.exit(sink())
    computeDIAReferencePeptideAccessionCleanup(manager, ref, list(delimiter = ";", aggregation_method = "mean", fasta_table = fasta_path), variant = variant, now_fn = function() as.POSIXct("2026-10-23", tz = "UTC")) }),
  stop("unknown step"))
out <- list(schema = "kit.variant_check.v1", ticket = "KIT-057", step = step, label = label, date = format(Sys.Date()),
  head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), source_state = spec$source, results = list())
withDIAStateScope(manager, function() {
  ref <- getDIAStateReference(manager, spec$source)
  values <- list()
  run_variants <- if (startsWith(mode, "run:")) sub("^run:", "", mode) else if (identical(mode, "compare")) character() else spec$variants
  if (identical(mode, "compare")) for (variant in spec$variants) values[[variant]] <- readRDS(record_path(variant))
  for (variant in run_variants) {
    gc(); t0 <- Sys.time(); before <- hwm()
    value <- spec$compute(ref, variant)
    seconds <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1); peak <- hwm()
    cp <- value$checkpoint
    matrix_ref <- cp$components$peptide_matrix
    m <- if (!is.null(matrix_ref)) { h <- hydrateSavedWorkflowEntity(manager$.context$store, matrix_ref); if (is.list(h) && !is.null(h$value)) h$value else h } else NULL
    rec <- list(seconds = seconds, process_hwm_gib = round(peak, 2), hwm_before_gib = round(before, 2),
      state_name = value$state_name, config = value$config, execution = value$execution[setdiff(names(value$execution), c("batch_sizes"))],
      tables = lapply(cp$components[intersect(names(cp$components), c("peptide_data", "design_matrix"))], table_sha),
      matrix = if (!is.null(m)) list(dim = dim(m), sha256 = digest::digest(m, algo = "sha256"), row_head = head(rownames(m), 3), col_head = head(colnames(m), 3)) else NULL,
      args = cp$slots$args, slots = cp$slots[setdiff(names(cp$slots), "args")])
    values[[variant]] <- rec
    if (startsWith(mode, "run:")) saveRDS(rec, record_path(variant))
    cat(sprintf("%-11s %6.1f s | HWM %5.2f GiB | matrix %s | table %s\n", variant, seconds, peak, paste(rec$matrix$dim, collapse = "x"), substr(rec$tables$peptide_data$sha256, 1, 12)))
  }
  if (startsWith(mode, "run:")) return(invisible(NULL))
  w <- values$whole; s <- values[[spec$variants[[2]]]]
  comparisons <- list(matrix_identical = identical(w$matrix$sha256, s$matrix$sha256), matrix_dimnames_identical = identical(w$matrix[c("dim","row_head","col_head")], s$matrix[c("dim","row_head","col_head")]),
    peptide_data_identical = identical(w$tables$peptide_data$sha256, s$tables$peptide_data$sha256), design_identical = identical(w$tables$design_matrix$sha256, s$tables$design_matrix$sha256),
    config_identical = identical(w$config, s$config), args_identical = identical(w$args, s$args), slots_identical = identical(w$slots, s$slots), state_identical = identical(w$state_name, s$state_name))
  if (!comparisons$args_identical) { diffs <- names(w$args)[!vapply(names(w$args), function(n) identical(w$args[[n]], s$args[[n]]), logical(1))]; comparisons$args_differing <- diffs }
  if (!comparisons$slots_identical) comparisons$slots_differing <- names(w$slots)[!vapply(names(w$slots), function(n) identical(w$slots[[n]], s$slots[[n]]), logical(1))]
  out$results <<- lapply(values, function(v) v[setdiff(names(v), c("args", "slots"))]); out$comparisons <<- comparisons
  cat("EXACT:", all(unlist(comparisons[c("matrix_identical","peptide_data_identical","design_identical","config_identical","args_identical","slots_identical")])), "|", paste(names(comparisons), vapply(comparisons, function(x) paste(x, collapse = ","), ""), sep = "=", collapse = " "), "\n")
})
if (startsWith(mode, "run:")) quit(save = "no")
jsonlite::write_json(out, file.path("docs/qa/workflow-execution/kit/kit-057", paste0("variant-check-", step, "-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
