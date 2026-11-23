# KIT-173 section 9.1 (R2e): the effect of GO parent-term propagation on proteomics enrichment results, descriptive.
#
# Usage, from the package root on a run-tests runner, after the fix commit is installed:
#   Rscript docs/qa/workflow-execution/kit/kit-173/effect-measurement.R <out dir>
#
# The committed E. coli K-12 table (tests/testdata/go-propagation/uniprot-ecoli-83333.tsv.gz, UniProt 2026_04, reviewed)
# through `processEnrichments` (taxon 999999, q 0.05, up cutoff 0) with go_propagation "none", "is_a_part_of" and
# "is_a_part_of_regulates", on identical lists for the three values. Arms (rule fixed in the ticket; terms drawn in F1,
# effect-planted-terms.tsv): 1 specific GO:0032790, 2 mid-level GO:0006396, 3 broad GO:0044282 (sizes under the shipped
# value), 20 draws each: half the term's proteins plus as many drawn from the other fixture proteins; 4 null, 20 random
# lists of the median planted size. Draw d of arm a uses seed 20261123 + 100 a + d. The DA table covers every fixture
# protein: the list at q 0.01 and log2FC 1, the others at q 0.5 and log2FC 0. The kernel's plot files are not written
# (the save functions are replaced by no-ops; no result depends on them). No pass rule (D8): the tables go to the owner.
# Writes effect-draws.tsv (one row per arm, draw and value), effect.tsv (per arm and value) and effect-sets.tsv.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "k173-effect-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(MultiScholaR))
h <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
ns <- asNamespace("MultiScholaR")
build <- get("buildUniprotGoGeneSets", ns)
restore_widget <- h$goPropReplaceBinding("MultiScholaR", "saveSelfContainedWidget", function(...) invisible(NULL))
restore_ggsave <- h$goPropReplaceBinding("ggplot2", "ggsave", function(...) invisible(NULL))

raw <- h$goPropReadTsv("uniprot-ecoli-83333.tsv.gz")
annotations <- h$goPropConvert(raw)
entries <- sort(raw$Entry)
values <- c("none", "is_a_part_of", "is_a_part_of_regulates")
arms <- data.frame(arm = 1:4, name = c("specific", "mid", "broad", "null"),
    term = c("GO:0032790", "GO:0006396", "GO:0044282", NA), stringsAsFactors = FALSE)
q <- 0.05
roots <- c("GO:0008150", "GO:0003674", "GO:0005575")

# The gene sets per value (the builder `processEnrichments` calls), their sizes among the fixture's proteins.
sets <- lapply(stats::setNames(values, values), function(value) build(annotations, value))
universe <- lapply(sets, function(s) intersect(entries, unique(s$term2gene[[2L]])))
sizes <- lapply(values, function(value) table(sets[[value]]$term2gene[[1L]][sets[[value]]$term2gene[[2L]] %in%
    universe[[value]]]))
names(sizes) <- values
in_range <- lapply(sizes, function(x) names(x)[x >= 10 & x <= 500])
aspect_of <- AnnotationDbi::Ontology(GO.db::GOTERM)
set_rows <- do.call(rbind, lapply(values, function(value) data.frame(value = value,
    background_N = length(universe[[value]]), pairs = nrow(sets[[value]]$term2gene),
    terms = length(sizes[[value]]), terms_in_range = length(in_range[[value]]),
    in_range_BP = sum(aspect_of[in_range[[value]]] == "BP", na.rm = TRUE),
    in_range_MF = sum(aspect_of[in_range[[value]]] == "MF", na.rm = TRUE),
    in_range_CC = sum(aspect_of[in_range[[value]]] == "CC", na.rm = TRUE),
    entering_at_10 = length(setdiff(in_range[[value]], names(sizes$none)[sizes$none >= 10])),
    leaving_at_500 = length(intersect(names(sizes[[value]])[sizes[[value]] > 500], in_range$none)),
    root_only_proteins = if (identical(value, "none")) NA_integer_ else
        length(intersect(attr(sets[[value]], "root_only_proteins"), entries)),
    go_db_version = attr(sets[[value]], "go_db_version") %||% NA,
    go_source_date = attr(sets[[value]], "go_source_date") %||% NA, stringsAsFactors = FALSE)))
utils::write.table(set_rows, file.path(out, "effect-sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(set_rows)

ancestors <- function(term) setdiff(GO.db::GOBPANCESTOR[[term]], "all")
offspring <- function(term) GO.db::GOBPOFFSPRING[[term]]
planted_members <- lapply(arms$term[1:3], function(term) sort(unique(sets$is_a_part_of_regulates$term2gene[[2L]][
    sets$is_a_part_of_regulates$term2gene[[1L]] == term])))
planted_sizes <- 2L * (lengths(planted_members) %/% 2L)
null_size <- as.integer(stats::median(planted_sizes))

drawList <- function(arm, draw) {
    set.seed(20261123L + 100L * arm + draw)
    if (arm == 4L) return(sort(sample(entries, null_size)))
    members <- planted_members[[arm]]
    half <- sample(members, length(members) %/% 2L)
    sort(c(half, sample(setdiff(entries, members), length(half))))
}

runOne <- function(list_ids, value) {
    table <- data.frame(Protein.IDs = entries, fdr_qvalue = ifelse(entries %in% list_ids, 0.01, 0.5),
        log2FC = ifelse(entries %in% list_ids, 1, 0), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "effect"),
        da_data = list(effect = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    results <- NULL
    utils::capture.output(results <- suppressMessages(suppressWarnings(processEnrichments(da, "999999", q_cutoff = q,
        pathway_dir = tempdir(), go_annotations = annotations, exclude_iea = FALSE, protein_id_column = "Protein.IDs",
        contrast_names = "effect", go_propagation = value))))
    up <- results@enrichment_data$effect$up
    if (is.null(up)) NULL else up@result
}

rows <- list()
for (a in arms$arm) {
    term <- arms$term[[a]]
    related <- if (is.na(term)) NULL else list(ancestors = ancestors(term), offspring = offspring(term))
    for (d in 1:20) {
        list_ids <- drawList(a, d)
        for (value in values) {
            started <- Sys.time()
            res <- runOne(list_ids, value)
            tested <- if (is.null(res)) character() else as.character(res$ID)
            set_size <- if (is.null(res)) integer() else as.integer(sub("/.*$", "", res$BgRatio))
            sig_q <- if (is.null(res)) character() else tested[!is.na(res$qvalue) & res$qvalue < q]
            sig_bh <- if (is.null(res)) character() else tested[res$p.adjust < q]
            kind <- function(ids) if (is.na(term)) c(planted = 0L, ancestor = 0L, descendant = 0L,
                other = length(ids)) else c(planted = sum(ids == term), ancestor = sum(ids %in% related$ancestors),
                descendant = sum(ids %in% related$offspring),
                other = sum(!ids %in% c(term, related$ancestors, related$offspring)))
            kq <- kind(sig_q)
            rows[[length(rows) + 1L]] <- data.frame(arm = arms$name[[a]], planted_term = term %||% NA, draw = d,
                value = value, list_size = length(list_ids),
                background_N = if (is.null(res)) NA_integer_ else as.integer(unique(sub("^.*/", "", res$BgRatio)))[[1L]],
                tested = length(tested), tested_BP = sum(aspect_of[tested] == "BP", na.rm = TRUE),
                tested_MF = sum(aspect_of[tested] == "MF", na.rm = TRUE),
                tested_CC = sum(aspect_of[tested] == "CC", na.rm = TRUE), tested_roots = sum(tested %in% roots),
                size_q1 = if (length(set_size)) stats::quantile(set_size, 0.25) else NA,
                size_median = if (length(set_size)) stats::median(set_size) else NA,
                size_q3 = if (length(set_size)) stats::quantile(set_size, 0.75) else NA,
                significant_q = length(sig_q), significant_bh = length(sig_bh),
                qvalue_all_na = !is.null(res) && all(is.na(res$qvalue)),
                sig_q_planted = kq[["planted"]], sig_q_ancestor = kq[["ancestor"]],
                sig_q_descendant = kq[["descendant"]], sig_q_other = kq[["other"]],
                planted_tested = !is.na(term) && term %in% tested, planted_recovered_q = !is.na(term) && term %in% sig_q,
                planted_recovered_bh = !is.na(term) && term %in% sig_bh,
                seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
                stringsAsFactors = FALSE)
        }
    }
    message(sprintf("arm %s done", arms$name[[a]]))
}
draws <- do.call(rbind, rows)
utils::write.table(draws, file.path(out, "effect-draws.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summary <- do.call(rbind, lapply(split(draws, list(draws$arm, draws$value), drop = TRUE), function(x) data.frame(
    arm = x$arm[[1L]], value = x$value[[1L]], draws = nrow(x), list_size = stats::median(x$list_size),
    background_N = stats::median(x$background_N, na.rm = TRUE), tested_median = stats::median(x$tested),
    tested_BP = stats::median(x$tested_BP), tested_MF = stats::median(x$tested_MF), tested_CC = stats::median(x$tested_CC),
    size_median = stats::median(x$size_median, na.rm = TRUE),
    significant_q_median = stats::median(x$significant_q), significant_q_iqr = stats::IQR(x$significant_q),
    significant_bh_median = stats::median(x$significant_bh), draws_qvalue_all_na = sum(x$qvalue_all_na),
    draws_any_significant_q = sum(x$significant_q > 0), draws_any_significant_bh = sum(x$significant_bh > 0),
    sig_q_planted = sum(x$sig_q_planted), sig_q_ancestor = sum(x$sig_q_ancestor),
    sig_q_descendant = sum(x$sig_q_descendant), sig_q_other = sum(x$sig_q_other),
    planted_tested = sum(x$planted_tested), planted_recovered_q = sum(x$planted_recovered_q),
    planted_recovered_bh = sum(x$planted_recovered_bh), seconds_median = stats::median(x$seconds),
    stringsAsFactors = FALSE)))
summary <- summary[order(match(summary$arm, arms$name), match(summary$value, values)), ]
utils::write.table(summary, file.path(out, "effect.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE)
restore_widget()
restore_ggsave()
