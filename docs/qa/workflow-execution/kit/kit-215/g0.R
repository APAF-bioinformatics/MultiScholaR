# KIT-215 gate G0 (ticket 9 A1, 10 G0; stop rule 3). Evidence code under docs/qa: no package code. Runs on a
# GitHub-hosted runner (run-tests `scripts`) from the package root of the base commit (staging after KIT-205 and
# KIT-208), after it is installed. No service is called (proxies point at a closed port).
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-215/g0.R <out dir>
#
# (1) versions.tsv: R, enrichit, clusterProfiler, gson.
# (2) engine.tsv: the four empty-result conditions the plan audit decoded from enrichit 0.2.0 (ticket 2.1), each one
#     call of the installed `clusterProfiler::enricher` on a synthetic TERM2GENE (T1: 12 genes, T2: 3 genes, T3: 30
#     genes; the universe is every gene): (i) n = 0 -> NULL; (ii) m' = 0 -> NULL; (iii) a tested term holds a list
#     gene -> a result with at least one row; (iv) m = 0 -> NULL. `pass` compares the observed with the expected.
# (3) g0.tsv: per fixture family, comparison and direction: the list size, n (list genes in the route's TERM2GENE),
#     the annotated universe, m (terms of 10 to 500 universe genes), m' (those holding a list gene), whether the
#     installed `enricher` returns NULL for the list, and the class VAL-170 gives it (no_list, not_tested, refused_c3,
#     refused_c4, tested). Families: w1_module (the W1/W113 module fixture's table and annotation,
#     helper-module-ci-prot-enrich.R, as written: the captured runs fit their own DA and call a fake `enricher`);
#     t7_complete (KIT-173's T7 project, helper-go-propagation.R `goPropRunProject`); k175_incomplete (KIT-175's
#     `k175Incomplete`: the first 25 annotation rows removed); k208_lacks_up (KIT-208's T7 annotation: the p_term
#     proteins removed). The project families run the shipped GO route (taxon 999999, the shipped parent terms) and
#     record the DA tables the kernel receives; the lists are formed from them as `processEnrichments` forms them.
# Exit status 0; a difference is a result (engine.tsv `pass`).

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out_dir <- args[[1L]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeTsv <- function(x, name) utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
Sys.setenv(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9", HTTP_PROXY = "http://127.0.0.1:9",
    HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "", MULTISCHOLAR_TEST_MODE = "false")
options(multischolar.test_mode = FALSE)
suppressPackageStartupMessages(library(MultiScholaR))
ns <- asNamespace("MultiScholaR")
h <- new.env(parent = ns)
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
sys.source(file.path("tests", "testthat", "helper-module-ci-prot-enrich.R"), envir = h)
version_of <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")

# ---- (1) versions ---------------------------------------------------------------------------------------------------
versions <- data.frame(item = c("R", "enrichit", "clusterProfiler", "gson", "GO.db"),
    version = c(R.version.string, vapply(c("enrichit", "clusterProfiler", "gson", "GO.db"), version_of, "")),
    stringsAsFactors = FALSE)
writeTsv(versions, "versions.tsv")
print(versions)

# ---- the counts, as the ticket's 3 item 1 states them ---------------------------------------------------------------
counts <- function(genes, term2gene, background) {
    t2g_genes <- as.character(term2gene[[2L]])
    t2g_terms <- as.character(term2gene[[1L]])
    annotated <- unique(t2g_genes)
    universe_ids <- intersect(annotated, as.character(background))
    keep <- t2g_genes %in% universe_ids
    sets <- lapply(split(t2g_genes[keep], t2g_terms[keep]), unique)
    tested <- sets[lengths(sets) >= 10L & lengths(sets) <= 500L]
    genes <- unique(as.character(genes))
    list(n = sum(genes %in% annotated), universe = length(universe_ids), m = length(tested),
        m_prime = sum(vapply(tested, function(set) any(genes %in% set), NA)))
}
enrich <- function(genes, term2gene, background, term2name = NA) {
    suppressMessages(suppressWarnings(clusterProfiler::enricher(gene = genes, universe = background,
        TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 0.05, pAdjustMethod = "BH")))
}

# ---- (2) the engine conditions --------------------------------------------------------------------------------------
synthetic <- data.frame(TERM = c(rep("T1", 12L), rep("T2", 3L), rep("T3", 30L)),
    GENE = c(sprintf("G%02d", 1:12), "L1", "L2", "L3", sprintf("B%02d", 1:30)), stringsAsFactors = FALSE)
everything <- unique(synthetic$GENE)
engine_cases <- list(
    i_n0 = list(genes = c("Z1", "Z2", "Z3"), t2g = synthetic, expected = "NULL"),
    ii_mprime0 = list(genes = c("L1", "L2", "L3"), t2g = synthetic, expected = "NULL"),
    iii_result = list(genes = sprintf("G%02d", 1:5), t2g = synthetic, expected = "result"),
    iv_m0 = list(genes = c("L1", "L2", "L3"), t2g = synthetic[synthetic$TERM == "T2", , drop = FALSE],
        expected = "NULL"))
engine <- do.call(rbind, lapply(names(engine_cases), function(name) {
    case <- engine_cases[[name]]
    set.seed(215L)
    result <- enrich(case$genes, case$t2g, everything)
    observed <- if (is.null(result)) "NULL" else if (nrow(result@result) > 0L) "result" else "result_zero_rows"
    k <- counts(case$genes, case$t2g, everything)
    data.frame(case = name, list_size = length(case$genes), n = k$n, universe = k$universe, m = k$m,
        m_prime = k$m_prime, expected = case$expected, observed = observed,
        rows = if (is.null(result)) NA_integer_ else nrow(result@result), class = paste(class(result), collapse = "/"),
        pass = identical(observed, case$expected), stringsAsFactors = FALSE)
}))
writeTsv(engine, "engine.tsv")
print(engine)

# ---- (3) the fixture families ---------------------------------------------------------------------------------------
classify <- function(size, k, is_null) {
    if (size == 0L) return("no_list")
    if (!is_null) return("tested")
    if (k$universe == 0L) return("refused_c3")
    if (k$m == 0L) return("refused_c4")
    if (k$n == 0L) return("not_tested")
    "tested_no_term"
}
family_rows <- function(family, tables, annotations, id_column, q_cutoff = 0.05, up_cutoff = 0, down_cutoff = 0) {
    sets <- MultiScholaR:::buildUniprotGoGeneSets(annotations, MultiScholaR:::goPropagationShipped)
    do.call(rbind, lapply(names(tables), function(comparison) {
        table <- MultiScholaR:::enrichmentTestedRows(tables[[comparison]])$rows
        significant <- table[!is.na(table$fdr_qvalue) & table$fdr_qvalue < q_cutoff, , drop = FALSE]
        lists <- list(up = significant[[id_column]][significant$log2FC > up_cutoff],
            down = significant[[id_column]][significant$log2FC < -down_cutoff])
        background <- unique(table[[id_column]])
        do.call(rbind, lapply(names(lists), function(direction) {
            genes <- lists[[direction]]
            k <- counts(genes, sets$term2gene, background)
            set.seed(215L)
            is_null <- if (length(genes)) is.null(enrich(genes, sets$term2gene, background, sets$term2name)) else NA
            data.frame(family = family, comparison = comparison, direction = direction, list_size = length(genes),
                n = k$n, universe = k$universe, m = k$m, m_prime = k$m_prime, enricher_null = is_null,
                val170 = classify(length(genes), k, isTRUE(is_null)), stringsAsFactors = FALSE)
        }))
    }))
}

# The W1/W113 module fixture as written (6 proteins; the GO case's taxon 999999 annotation).
module_table <- h$module_ci_prot_enrich_da_table()
module_annotations <- h$module_ci_prot_enrich_annotations(organism = "custom")
rows <- list(family_rows("w1_module", list(`groupB-groupA` = module_table), module_annotations, "uniprot_acc"))
rows[[1L]]$run_on_base <- "not run (captured runs use a fake enricher)"

# The project families: the DA tables the kernel receives, recorded by a wrapper around the shipped kernel.
project_family <- function(family, annotations) {
    seen <- new.env(parent = emptyenv())
    kernel <- function(...) {
        arguments <- list(...)
        seen$tables <- arguments$da_results@da_data
        seen$id_column <- arguments$protein_id_column
        seen$q <- arguments$q_cutoff
        seen$up <- arguments$up_cutoff
        seen$down <- arguments$down_cutoff
        MultiScholaR:::processEnrichments(...)
    }
    root <- tempfile(paste0("kit215-g0-", family, "-"))
    dir.create(root)
    run <- tryCatch({
        utils::capture.output(suppressMessages(suppressWarnings(h$goPropRunProject(root, annotations,
            process_fn = kernel))))
        "completed"
    }, error = function(e) paste("error:", conditionMessage(e)))
    message(family, ": ", run)
    if (is.null(seen$tables)) {
        return(data.frame(family = family, comparison = NA, direction = NA, list_size = NA, n = NA, universe = NA,
            m = NA, m_prime = NA, enricher_null = NA, val170 = NA, run_on_base = run, stringsAsFactors = FALSE))
    }
    out <- family_rows(family, seen$tables, annotations, seen$id_column, seen$q, seen$up, seen$down)
    out$run_on_base <- run
    out
}
full <- h$goPropAnnotate(h$goPropReadTsv("t7-annotation-ecoli.tsv"), tempfile("kit215-g0-cache-"))
roles <- h$goPropReadTsv("t7-fixture-entries.tsv")
rows[[2L]] <- project_family("t7_complete", full)
rows[[3L]] <- project_family("k175_incomplete", full[-seq_len(25L), , drop = FALSE])
rows[[4L]] <- project_family("k208_lacks_up", full[!full$Entry %in% roles$entry[roles$role == "p_term"], ,
    drop = FALSE])
g0 <- do.call(rbind, rows)
writeTsv(g0, "g0.tsv")
print(g0)
writeLines(c(sprintf("engine conditions (i) to (iv) as decoded: %s", all(engine$pass)),
    sprintf("enrichit %s, clusterProfiler %s", version_of("enrichit"), version_of("clusterProfiler"))),
    file.path(out_dir, "g0-verdict.txt"))
