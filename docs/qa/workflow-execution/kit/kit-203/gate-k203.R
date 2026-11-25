# KIT-203 gates G0 (base) and G1 (head): KIT-189's lists of cell 1 (B1, null, L = 25, parent terms), replicates 1 to
# 400, through the shipped `processEnrichments()`. Evidence code under docs/qa: no package code. Runs on a GitHub-hosted
# runner (run-tests `scripts`) with MultiScholaR installed from the bound commit, from the package root.
#
# Usage:
#   Rscript docs/qa/workflow-execution/kit/kit-203/gate-k203.R <mode> <out dir> [<first> <last>]
#   mode  g0: the base commit; per list the written up TSV's `significant` count against KIT-189's A-Q `calls`, the
#             stored result's `p.adjust < 0.05` count against A-BH; also KIT-173's null draw 5 (KIT-189's E0) is
#             drawn, run (expected 156 tested, 120 q-value calls, 0 BH calls) and its 126 accessions written.
#         g1: the head commit; per list the default run's TSV `significant` count against A-BHf `calls` (and
#             `p_adjust_family`, `family_size` against the formula and m), then the run with
#             `significance_rule = "qvalue"`: its count against A-Q.
#   first, last  replicate range (default 1 to 400)
# Environment: K203_WORKERS (default 4).
#
# Before anything is compared: R and the versions of clusterProfiler, enrichit, qvalue, GO.db and AnnotationDbi are
# checked against KIT-189's provenance (kit-189/prov/provenance-s-1-1-1000.json); a difference is written to
# env-check.tsv and the script stops with exit status 3 (an environment mismatch, not a result).
# KIT-189's generator is parsed from kit-189/sim-k189.R: only the definitions named in GENERATOR below are evaluated;
# the SHA-256 of their deparsed text is recorded. Each list's `list_sha` is checked against the committed row.
# The first five lists run one after another with their timing logged; the rest on K203_WORKERS forked workers.
# Writes to <out dir>: env-check.tsv, generator.tsv, <mode>.tsv (one row per list), <mode>-summary.tsv, and under g0
# e0-accessions.tsv and e0.tsv.

args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(2L, 4L), args[[1L]] %in% c("g0", "g1"))
mode <- args[[1L]]
out_dir <- args[[2L]]
first_rep <- if (length(args) == 4L) as.integer(args[[3L]]) else 1L
last_rep <- if (length(args) == 4L) as.integer(args[[4L]]) else 400L
workers <- as.integer(Sys.getenv("K203_WORKERS", "4"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
elapsed <- function() proc.time()[["elapsed"]]
t_start <- elapsed()
k189 <- file.path("docs", "qa", "workflow-execution", "kit", "kit-189")
writeTsv <- function(x, name) utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")

suppressPackageStartupMessages(library(MultiScholaR))

# ---- environment against KIT-189's provenance ----------------------------------------------------------------------
prov <- jsonlite::read_json(file.path(k189, "prov", "provenance-s-1-1-1000.json"))
pkgs <- c("clusterProfiler", "enrichit", "qvalue", "GO.db", "AnnotationDbi")
env <- data.frame(item = c("R", pkgs),
    expected = c(prov$r_version, vapply(pkgs, function(p) prov$versions[[p]], "")),
    installed = c(R.version.string, vapply(pkgs, function(p) as.character(utils::packageVersion(p)), "")),
    stringsAsFactors = FALSE)
env$same <- env$expected == env$installed
writeTsv(env, "env-check.tsv")
print(env)
if (!all(env$same)) {
    writeLines("ENVIRONMENT MISMATCH: no comparison made")
    quit(status = 3L)
}

# ---- KIT-189's generator, parsed -----------------------------------------------------------------------------------
GENERATOR <- c("CELLS", "STRATA", "VALUES", "seedKit", "sortC", "b2Entries", "loadBackground", "drawReplicate")
exprs <- parse(file.path(k189, "sim-k189.R"), keep.source = FALSE)
picked <- Filter(function(e) is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) &&
    as.character(e[[2L]]) %in% GENERATOR, as.list(exprs))
stopifnot(identical(sort(vapply(picked, function(e) as.character(e[[2L]]), "")), sort(GENERATOR)))
g <- new.env(parent = globalenv())
g$h <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = g$h)
g$h$goPropFixtureDir <- function(...) file.path("tests", "testdata", "go-propagation", ...)
g$build <- get("buildUniprotGoGeneSets", asNamespace("MultiScholaR"))
g$elapsed <- elapsed
for (e in picked) eval(e, g)
generator <- data.frame(
    script_sha256 = digest::digest(file = file.path(k189, "sim-k189.R"), algo = "sha256"),
    definitions = paste(GENERATOR, collapse = ","),
    deparsed_sha256 = digest::digest(paste(unlist(lapply(picked, deparse, width.cutoff = 500L)), collapse = "\n"),
        algo = "sha256", serialize = FALSE),
    stringsAsFactors = FALSE)
writeTsv(generator, "generator.tsv")
print(generator)
h <- g$h
# The kernel's plot files are not written (no comparison depends on them), as in KIT-189's simulation.
restore_widget <- h$goPropReplaceBinding("MultiScholaR", "saveSelfContainedWidget", function(...) invisible(NULL))
restore_ggsave <- h$goPropReplaceBinding("ggplot2", "ggsave", function(...) invisible(NULL))
g$BG <- list(B1 = g$loadBackground("B1"))
bg <- g$BG$B1
writeLines(sprintf("set-up %.1f s; B1 universe %d proteins", elapsed() - t_start, length(bg$U)))

# ---- the committed rows ---------------------------------------------------------------------------------------------
rows <- utils::read.delim(gzfile(file.path(k189, "rows", "s-1-1-1000.tsv.gz")), colClasses = "character",
    na.strings = "NA", quote = "")
rows <- rows[rows$propagation == "is_a_part_of_regulates", , drop = FALSE]
rowOf <- function(r, arm) {
    x <- rows[as.integer(rows$replicate) == r & rows$arm == arm, , drop = FALSE]
    stopifnot(nrow(x) == 1L)
    x
}

# ---- one call through the shipped route (KIT-189 runRoute, without go_propagation: the shipped default) -------------
runRoute <- function(list_ids, extra = list()) {
    table <- data.frame(Protein.IDs = bg$entries, fdr_qvalue = ifelse(bg$entries %in% list_ids, 0.01, 0.5),
        log2FC = ifelse(bg$entries %in% list_ids, 1, 0), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "effect"),
        da_data = list(effect = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    pathway_dir <- tempfile("k203-pw-")
    dir.create(pathway_dir)
    on.exit(unlink(pathway_dir, recursive = TRUE), add = TRUE)
    warnings <- character()
    results <- NULL
    t0 <- elapsed()
    status <- tryCatch({
        utils::capture.output(results <- withCallingHandlers(suppressMessages(do.call(processEnrichments, c(list(
            da_results = da, taxon_id = "999999", q_cutoff = 0.05, pathway_dir = pathway_dir,
            go_annotations = bg$annotations, exclude_iea = FALSE, protein_id_column = "Protein.IDs",
            contrast_names = "effect"), extra))), warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }))
        "ok"
    }, error = function(e) paste("error:", conditionMessage(e)))
    if (identical(status, "ok") && any(grepl("Error processing up-regulated genes", warnings, fixed = TRUE))) {
        status <- "swallowed"
    }
    tsv_path <- file.path(pathway_dir, "effect_up_enrichment_results.tsv")
    tsv <- if (file.exists(tsv_path)) utils::read.delim(tsv_path, colClasses = "character", na.strings = character(),
        quote = "") else NULL
    list(status = status, up = if (is.null(results)) NULL else results@enrichment_data$effect$up, tsv = tsv,
        seconds = elapsed() - t0)
}
countTrue <- function(x) sum(x %in% "TRUE")
bhf <- function(up) {
    p <- up@result$pvalue
    stats::p.adjust(c(p, rep(1, length(up@geneSets) - length(p))), "BH")[seq_along(p)]
}

# ---- one list ------------------------------------------------------------------------------------------------------
runList <- function(r) {
    t0 <- elapsed()
    tryCatch({
        draw <- g$drawReplicate(as.list(g$CELLS[g$CELLS$idx == 1L, , drop = FALSE]), r)
        q_row <- rowOf(r, "Q")
        sha <- digest::digest(draw$list, algo = "sha256")
        base <- data.frame(replicate = r, list_sha_same = identical(sha, q_row$list_sha),
            m_row = as.integer(q_row$m), m_kept_row = as.integer(q_row$m_kept),
            q_calls_row = as.integer(q_row$calls), bh_calls_row = as.integer(rowOf(r, "BH")$calls),
            bhf_calls_row = as.integer(rowOf(r, "BHf")$calls), stringsAsFactors = FALSE)
        if (identical(mode, "g0")) {
            rc <- runRoute(draw$list)
            up <- rc$up
            ok <- identical(rc$status, "ok") && !is.null(up) && !is.null(rc$tsv)
            x <- data.frame(base, status = rc$status, m = if (ok) length(up@geneSets) else NA_integer_,
                nrow_tsv = if (ok) nrow(rc$tsv) else NA_integer_,
                q_calls_tsv = if (ok) countTrue(rc$tsv$significant) else NA_integer_,
                bh_calls = if (ok) sum(up@result$p.adjust < 0.05) else NA_integer_,
                bhf_calls_formula = if (ok) sum(bhf(up) < 0.05) else NA_integer_,
                seconds = round(rc$seconds, 2), stringsAsFactors = FALSE)
            x$match <- x$list_sha_same %in% TRUE & x$status == "ok" & x$m %in% x$m_row & x$nrow_tsv %in% x$m_kept_row &
                x$q_calls_tsv %in% x$q_calls_row & x$bh_calls %in% x$bh_calls_row &
                x$bhf_calls_formula %in% x$bhf_calls_row
        } else {
            rd <- runRoute(draw$list)
            rq <- runRoute(draw$list, list(significance_rule = "qvalue"))
            ok_d <- identical(rd$status, "ok") && !is.null(rd$up) && !is.null(rd$tsv)
            ok_q <- identical(rq$status, "ok") && !is.null(rq$up) && !is.null(rq$tsv)
            family_same <- ok_d && identical(rd$up@result$p_adjust_family, bhf(rd$up)) &&
                all(rd$up@result$family_size == length(rd$up@geneSets))
            x <- data.frame(base, status_default = rd$status, status_qvalue = rq$status,
                m = if (ok_d) length(rd$up@geneSets) else NA_integer_,
                nrow_tsv = if (ok_d) nrow(rd$tsv) else NA_integer_,
                bhf_calls_tsv = if (ok_d) countTrue(rd$tsv$significant) else NA_integer_,
                rule_column = if (ok_d) paste(unique(rd$tsv$significance_rule), collapse = ",") else NA_character_,
                family_same = family_same,
                q_calls_tsv = if (ok_q) countTrue(rq$tsv$significant) else NA_integer_,
                qvalue_columns_absent = ok_q && !any(c("p_adjust_family", "family_size", "significance_rule") %in%
                    c(names(rq$tsv), names(rq$up@result))),
                seconds = round(rd$seconds + rq$seconds, 2), stringsAsFactors = FALSE)
            x$match <- x$list_sha_same %in% TRUE & x$status_default == "ok" & x$status_qvalue == "ok" &
                x$m %in% x$m_row & x$nrow_tsv %in% x$m_kept_row & x$bhf_calls_tsv %in% x$bhf_calls_row &
                x$rule_column %in% "bh_whole_family" & x$family_same %in% TRUE & x$q_calls_tsv %in% x$q_calls_row &
                x$qvalue_columns_absent %in% TRUE
        }
        x$wall <- round(elapsed() - t0, 2)
        x
    }, error = function(e) data.frame(replicate = r, instrument_error = conditionMessage(e), stringsAsFactors = FALSE))
}

# ---- run -----------------------------------------------------------------------------------------------------------
reps <- first_rep:last_rep
lead <- utils::head(reps, 5L)
t_run <- elapsed()
lead_rows <- lapply(lead, function(r) {
    x <- runList(r)
    writeLines(sprintf("list %d: %s (%.1f s)", r, if ("match" %in% names(x)) paste("match", x$match) else
        paste("instrument error:", x$instrument_error), elapsed() - t_run))
    x
})
rest <- setdiff(reps, lead)
rest_rows <- if (length(rest)) parallel::mclapply(rest, runList, mc.cores = workers, mc.preschedule = TRUE) else list()
all_rows <- c(lead_rows, rest_rows)
bad <- !vapply(all_rows, function(x) is.data.frame(x) && "match" %in% names(x), logical(1))
if (any(bad)) {
    errors <- vapply(all_rows[bad], function(x) if (is.data.frame(x)) x$instrument_error else
        paste(as.character(x), collapse = " "), "")
    writeTsv(data.frame(replicate = reps[bad], error = gsub("[\t\r\n]+", " ", errors)), sprintf("%s-errors.tsv", mode))
}
table_out <- do.call(rbind, all_rows[!bad])
writeTsv(table_out, sprintf("%s.tsv", mode))
writeLines(sprintf("%d of %d lists compared in %.1f s; %d instrument errors", nrow(table_out), length(reps),
    elapsed() - t_run, sum(bad)))

summary <- if (identical(mode, "g0")) {
    data.frame(lists = length(reps), compared = nrow(table_out), instrument_errors = sum(bad),
        list_sha_same = sum(table_out$list_sha_same), status_ok = sum(table_out$status == "ok"),
        m_values = paste(sort(unique(table_out$m)), collapse = ","),
        m_kept_range = paste(range(table_out$nrow_tsv), collapse = "-"),
        q_lists_with_call = sum(table_out$q_calls_tsv > 0), q_lists_with_call_row = sum(table_out$q_calls_row > 0),
        q_calls_total = sum(table_out$q_calls_tsv), bh_lists_with_call = sum(table_out$bh_calls > 0),
        bh_calls_total = sum(table_out$bh_calls), bhf_lists_with_call = sum(table_out$bhf_calls_formula > 0),
        bhf_calls_total = sum(table_out$bhf_calls_formula), lists_matching = sum(table_out$match),
        verdict = if (!any(bad) && all(table_out$match)) "match on every list" else "MISMATCH or incomplete",
        stringsAsFactors = FALSE)
} else {
    data.frame(lists = length(reps), compared = nrow(table_out), instrument_errors = sum(bad),
        list_sha_same = sum(table_out$list_sha_same), m_values = paste(sort(unique(table_out$m)), collapse = ","),
        default_lists_with_call = sum(table_out$bhf_calls_tsv > 0),
        default_lists_with_call_row = sum(table_out$bhf_calls_row > 0),
        default_calls_total = sum(table_out$bhf_calls_tsv), family_same = sum(table_out$family_same),
        qvalue_lists_with_call = sum(table_out$q_calls_tsv > 0), qvalue_calls_total = sum(table_out$q_calls_tsv),
        lists_matching = sum(table_out$match),
        verdict = if (!any(bad) && all(table_out$match)) "match on every list" else "MISMATCH or incomplete",
        stringsAsFactors = FALSE)
}
writeTsv(summary, sprintf("%s-summary.tsv", mode))
print(t(summary))

# ---- g0: KIT-173's null draw 5 (KIT-189 E0, sim-k189.R:477-499), its accessions written for the fixture ------------
if (identical(mode, "g0")) {
    entries173 <- sort(bg$raw$Entry)
    prop <- bg$sets$is_a_part_of_regulates
    planted173 <- lapply(c("GO:0032790", "GO:0006396", "GO:0044282"), function(term)
        sort(unique(prop$term2gene[[2L]][prop$term2gene[[1L]] == term])))
    null_size <- as.integer(stats::median(2L * (lengths(planted173) %/% 2L)))
    g$seedKit(20261123L + 100L * 4L + 5L)
    list173 <- sort(sample(entries173, null_size))
    writeTsv(data.frame(entry = list173), "e0-accessions.tsv")
    rc <- runRoute(list173, list(go_propagation = "none"))
    e0 <- data.frame(list_size = length(list173), status = rc$status,
        tested = if (is.null(rc$tsv)) NA_integer_ else nrow(rc$tsv),
        q_calls_tsv = if (is.null(rc$tsv)) NA_integer_ else countTrue(rc$tsv$significant),
        bh_calls = if (is.null(rc$up)) NA_integer_ else sum(rc$up@result$p.adjust < 0.05),
        bhf_calls_formula = if (is.null(rc$up)) NA_integer_ else sum(bhf(rc$up) < 0.05),
        m = if (is.null(rc$up)) NA_integer_ else length(rc$up@geneSets),
        accessions_sha256 = digest::digest(file = file.path(out_dir, "e0-accessions.tsv"), algo = "sha256"),
        stringsAsFactors = FALSE)
    e0$as_kit189 <- e0$list_size == 126L && e0$tested %in% 156L && e0$q_calls_tsv %in% 120L && e0$bh_calls %in% 0L
    writeTsv(e0, "e0.tsv")
    print(e0)
}
restore_widget()
restore_ggsave()
writeLines(sprintf("done in %.1f s", elapsed() - t_start))
if (any(bad)) quit(status = 1L)
