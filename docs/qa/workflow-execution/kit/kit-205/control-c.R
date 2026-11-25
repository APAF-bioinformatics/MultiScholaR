# KIT-205 control C (ticket 9, A2; stop rule 2). Evidence code under docs/qa: no package code. Runs on a GitHub-hosted
# runner (run-tests `scripts`) from the package root of the commit under test (head), after head is installed.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-205/control-c.R <out dir> <base commit>
#
# The base commit's tree is installed into a private library. Each side runs in its own Rscript process
# (`--side`): base with that library, head with the installed head. Both build the same DA tables (seeded, on
# KIT-173's committed E. coli fixture, as tests/testthat/test-kit205-enrichment-background.R builds them) and run
# processEnrichments() on them, the clusterProfiler route (taxon 999999; direct annotations, and the shipped parent
# terms) with the real enricher, the g:Profiler route (taxon 9606) with `perform_enrichment` replaced by a recorder (no
# service call). Plot files are not written (their writers are stubbed); the results TSV is.
#   (a) tables with no untested row: head(table) must equal base(table) byte for byte;
#   (b) tables with untested rows: head(table) must equal base(table[tested, ]) byte for byte, where tested is a finite
#       fdr_qvalue and, where the column exists, a finite raw_pvalue (the same base-R subset the head code takes).
# Compared per case, contrast and direction: clusterProfiler, the TSV's SHA-256, `identical()` of `@result`, of the
# plot inputs and of `@universe`; g:Profiler, `identical()` of every argument the recorder received (the list's rows
# and the custom background). For (b) the change against base on the whole table is also given: the largest
# |log10 p| change over the terms both runs tested, and the terms only one run tested.
# Also writes the exact hypergeometric rows of ticket 2.5 (stats::phyper) to control-c-phyper.tsv.
# Writes <out dir>/control-c.tsv, control-c-phyper.tsv, control-c-verdict.txt and the side logs. Exit status 0 (a
# difference is a result: the verdict names it).

args <- commandArgs(trailingOnly = TRUE)

# --- the tables (as the test file's k205Ids, k205Untested, k205Table) ------------------------------------------------
k205Ids <- function(h) {
    annotations <- h$goPropConvert(h$goPropReadTsv("uniprot-ecoli-83333.tsv.gz"))
    sets <- MultiScholaR:::buildUniprotGoGeneSets(annotations, "none")$term2gene
    pairs <- unique(data.frame(term = as.character(sets[[1L]]), gene = as.character(sets[[2L]]),
        stringsAsFactors = FALSE))
    sizes <- table(pairs$term)
    term <- sort(names(sizes)[sizes >= 10L & sizes <= 40L], method = "radix")[[1L]]
    in_term <- sort(unique(pairs$gene[pairs$term == term]), method = "radix")
    pool <- sort(setdiff(unique(pairs$gene), in_term), method = "radix")
    pool_pairs <- pairs[pairs$gene %in% pool, , drop = FALSE]
    counts <- table(pool_pairs$term)
    big <- names(counts)[order(-as.integer(counts), names(counts), method = "radix")][[1L]]
    up <- sort(unique(pool_pairs$gene[pool_pairs$term == big]), method = "radix")[seq_len(15L)]
    rest <- setdiff(pool, up)
    down <- rest[seq(1L, by = 5L, length.out = 10L)]
    fill <- setdiff(rest[seq(3L, by = 5L, length.out = 120L)], down)[seq_len(85L)]
    list(annotations = annotations, term = term, members = in_term[seq_len(10L)], up = up, down = down, fill = fill,
        others = c(up, down, fill))
}

k205Table <- function(ids, untested = character(), raw_p = TRUE, seed = 205L) {
    accessions <- c(ids$members, ids$others)
    set.seed(seed)
    p <- stats::runif(length(accessions))
    planted <- c(ids$up, ids$down)
    p[match(planted, accessions)] <- seq(1e-6, 1e-5, length.out = length(planted))
    log2fc <- ifelse(accessions %in% ids$up, 1, ifelse(accessions %in% ids$down, -1, 0))
    p[accessions %in% untested] <- NA_real_
    log2fc[accessions %in% untested] <- NA_real_
    q <- rep(NA_real_, length(p))
    q[is.finite(p)] <- stats::p.adjust(p[is.finite(p)], method = "BH")
    table <- data.frame(Protein.IDs = accessions, raw_pvalue = p, fdr_qvalue = q, log2FC = log2fc,
        stringsAsFactors = FALSE)
    if (!raw_p) table$raw_pvalue <- NULL
    table
}

testedRows <- function(table) {
    tested <- is.finite(table$fdr_qvalue)
    if ("raw_pvalue" %in% names(table)) tested <- tested & is.finite(table$raw_pvalue)
    if (all(tested)) table else table[tested, , drop = FALSE]
}

#' The cases: name, kind (a or b), route, propagation and the tables
cases <- function(ids) {
    in_term <- c(ids$members[seq_len(6L)], ids$fill[c(5L, 19L, 33L, 47L, 61L, 75L)])
    spread <- ids$fill[c(2L, 16L, 30L, 44L, 58L, 72L, 7L, 21L, 35L, 49L, 63L, 77L)]
    no_p <- k205Table(ids, seed = 207L)
    no_p_rows <- match(ids$fill[c(9L, 23L)], no_p$Protein.IDs)
    no_p$raw_pvalue[no_p_rows] <- NA_real_
    no_q <- k205Table(ids, raw_p = FALSE, seed = 208L)
    no_q$fdr_qvalue[match(ids$fill[c(11L, 25L)], no_q$Protein.IDs)] <- NA_real_
    list(
        a_cluster_none = list(kind = "a", route = "clusterprofiler", propagation = "none",
            tables = list(c1 = k205Table(ids), c2 = k205Table(ids, seed = 206L))),
        a_cluster_shipped = list(kind = "a", route = "clusterprofiler", propagation = "is_a_part_of_regulates",
            tables = list(c1 = k205Table(ids))),
        a_cluster_k189_form = list(kind = "a", route = "clusterprofiler", propagation = "none",
            tables = list(c1 = k205Table(ids, raw_p = FALSE))),
        a_gprofiler = list(kind = "a", route = "gprofiler", propagation = NA_character_,
            tables = list(c1 = k205Table(ids), c2 = k205Table(ids, seed = 206L))),
        b_cluster_none = list(kind = "b", route = "clusterprofiler", propagation = "none",
            tables = list(c1 = k205Table(ids, in_term), c2 = k205Table(ids, spread, seed = 206L))),
        b_cluster_shipped = list(kind = "b", route = "clusterprofiler", propagation = "is_a_part_of_regulates",
            tables = list(c1 = k205Table(ids, in_term))),
        b_cluster_finite_q_no_p = list(kind = "b", route = "clusterprofiler", propagation = "none",
            tables = list(c1 = no_p)),
        b_cluster_k189_form = list(kind = "b", route = "clusterprofiler", propagation = "none",
            tables = list(c1 = no_q)),
        b_gprofiler = list(kind = "b", route = "gprofiler", propagation = NA_character_,
            tables = list(c1 = k205Table(ids, in_term), c2 = k205Table(ids, spread, seed = 206L))))
}

# --- one side ---------------------------------------------------------------------------------------------------------
runSide <- function(side, file) {
    suppressPackageStartupMessages(library(MultiScholaR))
    Sys.setenv(MULTISCHOLAR_TEST_MODE = "false")
    h <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
    restores <- list(h$goPropReplaceBinding("MultiScholaR", "saveSelfContainedWidget", function(...) invisible(NULL)),
        h$goPropReplaceBinding("ggplot2", "ggsave", function(...) invisible(NULL)))
    recorded <- new.env()
    recorded$calls <- list()
    restores <- c(restores, list(h$goPropReplaceBinding("MultiScholaR", "perform_enrichment", function(...) {
        recorded$calls[[length(recorded$calls) + 1L]] <- list(...)
        NULL
    })))
    ids <- k205Ids(h)
    da <- function(tables) methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = names(tables)),
        da_data = tables, design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    run <- function(case, tables) {
        dir <- tempfile("k205-c-")
        dir.create(dir)
        recorded$calls <- list()
        results <- NULL
        if (identical(case$route, "clusterprofiler")) {
            utils::capture.output(results <- suppressMessages(suppressWarnings(processEnrichments(da(tables), "999999",
                q_cutoff = 0.05, pathway_dir = dir, go_annotations = ids$annotations, exclude_iea = FALSE,
                protein_id_column = "Protein.IDs", contrast_names = names(tables), go_propagation = case$propagation))))
            out <- list()
            for (contrast in names(tables)) for (direction in c("up", "down")) {
                r <- results@enrichment_data[[contrast]][[direction]]
                tsv <- file.path(dir, paste0(contrast, "_", direction, "_enrichment_results.tsv"))
                out[[paste0(contrast, "/", direction)]] <- if (is.null(r)) NULL else list(result = r@result,
                    plot_inputs = attr(r, "enrichment_plot_inputs", exact = TRUE), universe = r@universe,
                    tsv_sha256 = if (file.exists(tsv)) digest::digest(file = tsv, algo = "sha256") else NA_character_)
            }
            out
        } else {
            utils::capture.output(suppressMessages(suppressWarnings(processEnrichments(da(tables), "9606",
                q_cutoff = 0.05, pathway_dir = dir, exclude_iea = FALSE, protein_id_column = "Protein.IDs",
                contrast_names = names(tables)))))
            calls <- recorded$calls
            keys <- vapply(seq_along(calls), function(i) paste0(rep(names(tables), each = 2L)[[i]], "/",
                if (all(calls[[i]]$data_subset$log2FC > 0)) "up" else "down"), "")
            stats::setNames(calls, keys)
        }
    }
    store <- list()
    for (name in names(cases(ids))) {
        case <- cases(ids)[[name]]
        store[[name]] <- list(full = run(case, case$tables))
        if (identical(side, "base") && identical(case$kind, "b")) {
            store[[name]]$reduced <- run(case, lapply(case$tables, testedRows))
        }
        store[[name]]$left_out <- vapply(case$tables, function(t) nrow(t) - nrow(testedRows(t)), integer(1))
        cat(sprintf("side %s: case %s done\n", side, name))
    }
    for (restore in restores) restore()
    saveRDS(list(store = store, version = as.character(utils::packageVersion("MultiScholaR")),
        path = find.package("MultiScholaR"), formals = names(formals(processEnrichments)),
        tested_rule = exists("enrichmentTestedRows", envir = asNamespace("MultiScholaR"))), file)
}

if (length(args) && identical(args[[1L]], "--side")) {
    # The side's functions see the package's internals (the namespace's parents end at the global environment).
    suppressPackageStartupMessages(library(MultiScholaR))
    environment(runSide) <- asNamespace("MultiScholaR")
    runSide(args[[2L]], args[[3L]])
    quit(status = 0L)
}

# --- main ---------------------------------------------------------------------------------------------------------------
stopifnot(length(args) == 2L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
base <- args[[2L]]
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k205-control-c")
src <- file.path(work, "base-src")
lib <- file.path(work, "base-lib")
for (d in c(src, lib)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
rscript <- file.path(R.home("bin"), "Rscript")
archive_exit <- system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src)))
install_exit <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load",
    paste0("--library=", shQuote(lib)), shQuote(src)), stdout = file.path(out, "base-install.log.txt"),
    stderr = file.path(out, "base-install.log.txt"))
sides <- c(base = file.path(work, "base.rds"), head = file.path(work, "head.rds"))
exits <- c(base = system2(rscript, c(shQuote(script), "--side", "base", shQuote(sides[["base"]])),
        stdout = file.path(out, "side-base.log.txt"), stderr = file.path(out, "side-base.log.txt"),
        env = paste0("R_LIBS=", lib)),
    head = system2(rscript, c(shQuote(script), "--side", "head", shQuote(sides[["head"]])),
        stdout = file.path(out, "side-head.log.txt"), stderr = file.path(out, "side-head.log.txt")))
if (!all(file.exists(sides))) {
    writeLines(c(sprintf("archive %d install %d side exits %s", archive_exit, install_exit,
        paste(names(exits), exits, collapse = ", ")), "VERDICT: INSTRUMENT (a side wrote no result)"),
        file.path(out, "control-c-verdict.txt"))
    quit(status = 0L)
}
b <- readRDS(sides[["base"]])
h <- readRDS(sides[["head"]])

logp <- function(entry) if (is.null(entry$result) || !nrow(entry$result)) stats::setNames(numeric(), character()) else
    stats::setNames(log10(entry$result$pvalue), entry$result$ID)
rows <- list()
for (name in names(h$store)) {
    kind <- substr(name, 1L, 1L)
    route <- if (grepl("gprofiler", name)) "gprofiler" else "clusterprofiler"
    head_full <- h$store[[name]]$full
    against <- if (identical(kind, "a")) b$store[[name]]$full else b$store[[name]]$reduced
    base_full <- b$store[[name]]$full
    for (key in sort(union(names(head_full), names(against)))) {
        x <- head_full[[key]]
        y <- against[[key]]
        contrast <- sub("/.*$", "", key)
        equal <- if (identical(route, "gprofiler")) c(arguments = identical(x, y), custom_bg = identical(x$custom_bg,
            y$custom_bg), list_rows = identical(x$data_subset, y$data_subset)) else c(result = identical(x$result, y$result),
            plot_inputs = identical(x$plot_inputs, y$plot_inputs), universe = identical(x$universe, y$universe),
            tsv = identical(x$tsv_sha256, y$tsv_sha256))
        change <- NA_real_
        only_head <- NA_integer_
        only_base <- NA_integer_
        if (identical(kind, "b") && identical(route, "clusterprofiler")) {
            ph <- logp(x)
            pb <- logp(base_full[[key]])
            both <- intersect(names(ph), names(pb))
            change <- if (length(both)) max(abs(ph[both] - pb[both])) else NA_real_
            only_head <- length(setdiff(names(ph), names(pb)))
            only_base <- length(setdiff(names(pb), names(ph)))
        }
        rows[[length(rows) + 1L]] <- data.frame(case = name, kind = kind, route = route, contrast = contrast,
            direction = sub("^.*/", "", key), rows_left_out = h$store[[name]]$left_out[[contrast]],
            compared = if (identical(kind, "a")) "head(table) = base(table)" else "head(table) = base(table[tested, ])",
            checks = paste(names(equal), collapse = ","), all_equal = all(equal),
            not_equal = paste(names(equal)[!equal], collapse = ","),
            max_abs_log10p_change_vs_base_whole_table = signif(change, 6), terms_only_head = only_head,
            terms_only_base = only_base, stringsAsFactors = FALSE)
    }
}
table <- do.call(rbind, rows)
utils::write.table(table, file.path(out, "control-c.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
phyper <- data.frame(universe = c("today: 3,000 annotated rows, 300 untested", "untested spread like the rest (5 in the term)",
    "none of the 300 in the term", "20 of the term's 50 untested", "2 % untested (60), none in the term"),
    N = c(3000, 2700, 2700, 2700, 2940), K = c(50, 45, 50, 30, 50), n = 150)
phyper$expected <- phyper$n * phyper$K / phyper$N
phyper$p_x_ge_8 <- stats::phyper(7, phyper$K, phyper$N - phyper$K, phyper$n, lower.tail = FALSE)
utils::write.table(phyper, file.path(out, "control-c-phyper.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
expected_rows <- sum(vapply(h$store, function(s) length(s$full), integer(1)))
ok <- nrow(table) > 0L && all(table$all_equal) && nrow(table) == expected_rows &&
    isTRUE(h$tested_rule) && !isTRUE(b$tested_rule) && all(table$rows_left_out[table$kind == "b"] > 0L)
verdict <- c(sprintf("base %s (library %s, enrichmentTestedRows: %s); head %s (library %s, enrichmentTestedRows: %s)",
        base, b$path, b$tested_rule, system2("git", c("rev-parse", "HEAD"), stdout = TRUE), h$path, h$tested_rule),
    sprintf("archive %d, install %d, side exits %s", archive_exit, install_exit, paste(names(exits), exits, collapse = ", ")),
    sprintf("rows compared %d (head keys %d); all equal: %d; not equal: %d", nrow(table), expected_rows,
        sum(table$all_equal), sum(!table$all_equal)),
    if (ok) "VERDICT: control C holds ((a) head = base; (b) head = base on the table without its untested rows)" else
        "VERDICT: STOP (a difference, a missing row or a side that did not run the code under test)")
writeLines(verdict, file.path(out, "control-c-verdict.txt"))
cat(verdict, sep = "\n")
print(table)
print(phyper)
