# KIT-173 T5 evidence (section 5.2): agreement of the route's yeast GO sets with g:Profiler's propagated sets.
#
# Usage, from the package root on a run-tests runner, after the fix commit is installed:
#   Rscript docs/qa/workflow-execution/kit/kit-173/t5-agreement.R <out dir>
#
# The same chain and measures as the test T5 (helper-go-propagation.R): the committed yeast UniProt table through
# `uniprotGoIdToTerm` and `processEnrichments`, TERM2GENE captured as passed to `clusterProfiler::enricher`, for the
# three values; the scope is set by the shipped value. Writes:
#   t5-agreement.tsv   per value and term: aspect, sizes, containment C, precision P, Jaccard
#   t5-summary.tsv     per value: terms, median and IQR of C and P; the rule's three conditions for the shipped value
#   t5-two-relation.tsv  the two-relation option, descriptive: C over the terms with no descendant reached only through a
#                      regulates edge, the other BP terms' C separately, P over all terms
#   t5-lowest.tsv      the 20 terms with the lowest C and the 20 with the lowest P (shipped value)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "k173-t5-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(MultiScholaR))
h <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
restore_widget <- h$goPropReplaceBinding("MultiScholaR", "saveSelfContainedWidget", function(...) invisible(NULL))
restore_ggsave <- h$goPropReplaceBinding("ggplot2", "ggsave", function(...) invisible(NULL))

raw <- h$goPropReadTsv("uniprot-yeast-559292.tsv.gz")
annotations <- h$goPropConvert(raw)
gmt <- h$goPropReadGmt()
values <- c("none", "is_a_part_of", "is_a_part_of_regulates")
sets <- lapply(stats::setNames(values, values), function(value) h$goPropYeastTerm2Gene(annotations, value, tempdir()))
agreement <- h$goPropAgreement(sets, raw, gmt, "is_a_part_of_regulates")
agreement$jaccard <- agreement$both / (agreement$ours + agreement$reference - agreement$both)
utils::write.table(agreement, file.path(out, "t5-agreement.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

oln <- h$goPropOln(raw)
summary <- do.call(rbind, lapply(values, function(value) {
    x <- agreement[agreement$value == value, ]
    data.frame(value = value, terms = nrow(x), C_median = stats::median(x$C), C_q1 = stats::quantile(x$C, 0.25),
        C_q3 = stats::quantile(x$C, 0.75), P_median = stats::median(x$P, na.rm = TRUE),
        P_q1 = stats::quantile(x$P, 0.25, na.rm = TRUE), P_q3 = stats::quantile(x$P, 0.75, na.rm = TRUE),
        P_undefined = sum(is.na(x$P)), jaccard_median = stats::median(x$jaccard), stringsAsFactors = FALSE)
}))
shipped <- summary[summary$value == "is_a_part_of_regulates", ]
none <- summary[summary$value == "none", ]
summary$rule_C_ge_0.90 <- summary$value == "is_a_part_of_regulates" & shipped$C_median >= 0.90
summary$rule_P_ge_0.90 <- summary$value == "is_a_part_of_regulates" & shipped$P_median >= 0.90
summary$rule_C_gain_ge_0.20 <- summary$value == "is_a_part_of_regulates" & (shipped$C_median - none$C_median) >= 0.20
utils::write.table(summary, file.path(out, "t5-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE)
writeLines(c(sprintf("uniprot_entries\t%d", nrow(raw)), sprintf("gmt_genes\t%d", length(unique(gmt$GENE))),
    sprintf("gmt_genes_unmapped\t%d", sum(!unique(gmt$GENE) %in% oln$GENE)),
    sprintf("gmt_terms\t%d", length(unique(gmt$TERM))),
    sprintf("gmt_terms_not_in_go_db\t%d", sum(!unique(gmt$TERM) %in% names(AnnotationDbi::Ontology(GO.db::GOTERM))))),
    file.path(out, "t5-mapping.tsv"))

# Two-relation option: BP terms gaining descendants through regulates edges (GO.db ancestor table against the
# recursive-SQL two-relation closure), descriptive.
bp <- AnnotationDbi::toTable(GO.db::GOBPANCESTOR)
bp <- paste(bp[[1L]], bp[[2L]])[bp[[2L]] != "all"]
sql <- DBI::dbGetQuery(GO.db::GO_dbconn(), paste("WITH RECURSIVE anc(c, a) AS (SELECT _id, _parent_id FROM go_bp_parents",
    "WHERE relationship_type IN ('isa', 'part of') UNION SELECT anc.c, p._parent_id FROM anc JOIN go_bp_parents p",
    "ON p._id = anc.a WHERE p.relationship_type IN ('isa', 'part of')) SELECT t1.go_id AS term, t2.go_id AS ancestor",
    "FROM anc JOIN go_term t1 ON t1._id = anc.c JOIN go_term t2 ON t2._id = anc.a WHERE t2.go_id <> 'all'"))
only_five <- setdiff(bp, paste(sql$term, sql$ancestor))
gaining <- unique(sub("^\\S+ ", "", only_five))
two <- agreement[agreement$value == "is_a_part_of", ]
two_rows <- data.frame(measure = c("bp_terms_gaining_descendants_through_regulates", "terms_in_scope",
    "terms_without_regulates_only_descendant", "C_median_those_terms", "terms_with_regulates_only_descendant",
    "C_median_those_terms", "P_median_all_terms"),
    value = c(length(gaining), nrow(two), sum(!two$term %in% gaining), stats::median(two$C[!two$term %in% gaining]),
        sum(two$term %in% gaining), stats::median(two$C[two$term %in% gaining]), stats::median(two$P, na.rm = TRUE)))
utils::write.table(two_rows, file.path(out, "t5-two-relation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(two_rows)

s <- agreement[agreement$value == "is_a_part_of_regulates", ]
lowest <- rbind(cbind(by = "C", utils::head(s[order(s$C, s$term), ], 20L)),
    cbind(by = "P", utils::head(s[order(s$P, s$term), ], 20L)))
utils::write.table(lowest, file.path(out, "t5-lowest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
restore_widget()
restore_ggsave()
