# KIT-191 U14 (ticket 10): every figure the package shows about the four tests, regenerated from KIT-183's committed
# tables (docs/qa/workflow-execution/kit/kit-183/s/error-rates.tsv, power.tsv, familywise.tsv; run 45225854502) and
# compared with the package: each row of INTEGRATION_ENRICHMENT_TEST_EVIDENCE (value_text equal), and each figure of the
# notes, the few-samples sentences (the strength and limits sentences read the table). Kinds: metabolites =
# metabolomics_LCMS; genes and proteins = proteomics and transcriptomics. Tests: camera, camera_floor, permutation
# (the table's pcgse_perm), gsea. Rounding as the scorer: sprintf("%.1f", 100 * x) for error rates; power as a share
# "%.2f" (or "%.3f" when the largest is below 0.1), as a percentage "%.0f" of 100 x (or "%.1f" below 0.1; 0 as "0").
#
# Usage (from the installed head's package root): Rscript evidence-k191.R <out dir>
# Writes <out dir>/evidence.tsv, texts.tsv, verdict.txt.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
out_dir <- arguments[[1L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressMessages(library(MultiScholaR))
ns <- asNamespace("MultiScholaR")
s <- "docs/qa/workflow-execution/kit/kit-183/s"
er <- utils::read.delim(file.path(s, "error-rates.tsv"), stringsAsFactors = FALSE, colClasses = c(rho = "character"))
pw <- utils::read.delim(file.path(s, "power.tsv"), stringsAsFactors = FALSE)
fw <- utils::read.delim(file.path(s, "familywise.tsv"), stringsAsFactors = FALSE)
views <- list(metabolites = "metabolomics_LCMS", genes_proteins = c("proteomics", "transcriptomics"))
tests <- c(camera = "camera", camera_floor = "camera_floor", permutation = "pcgse_perm", gsea = "gsea")

P <- function(x) sprintf("%.1f", 100 * x)
range_text <- function(v, f = P) { a <- f(min(v)); b <- f(max(v)); if (a == b) a else paste(a, "to", b) }
power_share <- function(v) { f <- if (max(v) >= 0.1) "%.2f" else "%.3f"; range_text(v, function(x) sprintf(f, x)) }
power_pct <- function(x) if (x == 0) "0" else if (x >= 0.1) sprintf("%.0f", 100 * x) else sprintf("%.1f", 100 * x)
cells <- function(n, kind, factors, test, classes, rho = NULL) {
    rows <- er[er$n == n & er$view %in% views[[kind]] & er$kind %in% factors & er$test == tests[[test]] & er$class %in% classes, ]
    if (!is.null(rho)) rows <- rows[rows$rho == rho, ]
    stopifnot(nrow(rows) > 0L)
    rows
}
largest <- function(rows) { r <- rows[which.max(rows$T), ]; sprintf("95 %% interval %s to %s %%", P(r$lower), P(r$upper)) }
power <- function(n, kind, test, sizes) pw$power[pw$n == n & pw$view %in% views[[kind]] & pw$test == tests[[test]] &
    pw$size %in% sizes]

rows <- list()
add <- function(kind, test, n, measure, text) rows[[length(rows) + 1L]] <<- data.frame(kind = kind, test = test, samples = n,
    measure = measure, regenerated = text, stringsAsFactors = FALSE)
for (n in c(12L, 6L, 20L)) for (kind in names(views)) for (test in names(tests)) {
    base <- if (test == "gsea") c("C0", "C2") else c("C0", "C1", "C2")
    A <- cells(n, kind, "F_ind", test, base)
    text <- sprintf("%s %% (largest: %s)", range_text(A$T), largest(A))
    if (test == "gsea") text <- sprintf("%s; terms whose features vary together: %s %%", text,
        range_text(cells(n, kind, "F_ind", test, "C1")$T))
    add(kind, test, n, "A", text)
    B <- cells(n, kind, "F_ind", test, "C4")
    add(kind, test, n, "B", sprintf("%s %% (largest: %s)", range_text(B$T), largest(B)))
    if (n == 12L) {
        C <- cells(n, kind, "F_main", test, base)
        text <- paste(range_text(C$T), "%")
        if (test == "gsea") text <- sprintf("%s; terms whose features vary together: %s %%", text,
            range_text(cells(n, kind, "F_main", test, "C1")$T))
        add(kind, test, n, "C", text)
        add(kind, test, n, "D", paste(range_text(cells(n, kind, "F_main", test, "C1", "0.3")$T), "%"))
        E <- cells(n, kind, c("F_weak", "F_none"), test, c("C4", "C4b", "C1out"))
        e <- E[which.max(E$T), ]
        add(kind, test, n, "E", sprintf("%s %% (95 %% interval %s to %s %%)", P(e$T), P(e$lower), P(e$upper)))
        familywise <- fw$familywise[fw$n == n & fw$view %in% views[[kind]] & fw$kind %in% c("F_ind", "F_main") &
            fw$test == tests[[test]]]
        add(kind, test, n, "F", paste(range_text(familywise), "%"))
        if (kind == "metabolites") add(kind, test, n, "C3", paste(range_text(cells(n, kind, "F_ind", test, "C3")$T), "%"))
    }
    if (kind == "metabolites") add(kind, test, n, "I", power_share(power(n, kind, test, 5L))) else {
        add(kind, test, n, "G", power_share(power(n, kind, test, 10L)))
        add(kind, test, n, "H", power_share(power(n, kind, test, 30L)))
    }
}
regenerated <- do.call(rbind, rows)
package <- ns$INTEGRATION_ENRICHMENT_TEST_EVIDENCE
key <- function(x) paste(x$kind, x$test, x$samples, x$measure)
regenerated$package <- package$value_text[match(key(regenerated), key(package))]
regenerated$equal <- !is.na(regenerated$package) & regenerated$package == regenerated$regenerated
utils::write.table(regenerated, file.path(out_dir, "evidence.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# The notes and few-samples sentences: each figure regenerated and looked for in the package's text.
texts <- ns$integrationGuiTexts()
AC <- function(kind, test) {
    base <- if (test == "gsea") c("C0", "C2") else c("C0", "C1", "C2")
    range_text(c(cells(12L, kind, "F_ind", test, base)$T, cells(12L, kind, "F_main", test, base)$T))
}
Bmax <- function(n, kind, test) P(max(cells(n, kind, "F_ind", test, "C4")$T))
Brange <- function(kind, test) range_text(cells(12L, kind, "F_ind", test, "C4")$T)
Arange <- function(kind, test) range_text(cells(12L, kind, "F_ind", test, if (test == "gsea") c("C0", "C2") else
    c("C0", "C1", "C2"))$T)
Crange <- function(kind, test) range_text(cells(12L, kind, "F_main", test, c("C0", "C1", "C2"))$T)
E <- function(kind, test) { x <- cells(12L, kind, c("F_weak", "F_none"), test, c("C4", "C4b", "C1out")); P(max(x$T)) }
D <- function(kind, test) range_text(cells(12L, kind, "F_main", test, "C1", "0.3")$T)
pp <- function(v) { a <- power_pct(min(v)); b <- power_pct(max(v)); if (a == b) a else paste(a, "to", b) }
checks <- list(
    c("enrich_test_camera_note_genes_proteins", sprintf("in %s %% of cases", AC("genes_proteins", "camera"))),
    c("enrich_test_camera_note_genes_proteins", sprintf("in up to %s %%", Bmax(12L, "genes_proteins", "camera"))),
    c("enrich_test_camera_note_genes_proteins", sprintf("in up to %s %% on factors", E("genes_proteins", "camera"))),
    c("enrich_test_camera_note_genes_proteins", sprintf("found in %s %% of cases", pp(power(12L, "genes_proteins", "camera", 10L)))),
    c("enrich_test_camera_note_genes_proteins", sprintf("30 such features in %s %%", pp(power(12L, "genes_proteins", "camera", 30L)))),
    c("enrich_test_camera_note_metabolites", sprintf("in %s %% of cases", AC("metabolites", "camera"))),
    c("enrich_test_camera_note_metabolites", sprintf("up to %s %% even when", Bmax(12L, "metabolites", "camera"))),
    c("enrich_test_camera_note_metabolites", sprintf("in up to %s %% on factors", E("metabolites", "camera"))),
    c("enrich_test_camera_note_metabolites", sprintf("found in %s %% of cases", pp(power(12L, "metabolites", "camera", 5L)))),
    c("enrich_test_camera_floor_note_genes_proteins", sprintf("in %s %% of cases against", Arange("genes_proteins", "camera_floor"))),
    c("enrich_test_camera_floor_note_genes_proteins", sprintf("in up to %s %% (the default test: up to %s %%)",
        Bmax(12L, "genes_proteins", "camera_floor"), Bmax(12L, "genes_proteins", "camera"))),
    c("enrich_test_camera_floor_note_genes_proteins", sprintf("found in %s %% of cases (the default test: %s %%; 30 features: %s %%)",
        pp(power(12L, "genes_proteins", "camera_floor", 10L)), pp(power(12L, "genes_proteins", "camera", 10L)),
        pp(power(12L, "genes_proteins", "camera_floor", 30L)))),
    c("enrich_test_camera_floor_note_genes_proteins", sprintf("it found at most %s %%",
        power_pct(max(power(6L, "genes_proteins", "camera_floor", c(10L, 30L)))))),
    c("enrich_test_camera_floor_note_genes_proteins", sprintf("called in up to %s %%, as with", E("genes_proteins", "camera_floor"))),
    c("enrich_test_camera_floor_note_metabolites", sprintf("in %s %% of cases against", Arange("metabolites", "camera_floor"))),
    c("enrich_test_camera_floor_note_metabolites", sprintf("(the default test: %s %%)", Arange("metabolites", "camera"))),
    c("enrich_test_camera_floor_note_metabolites", sprintf("found in %s %% of cases (the default test: %s %%)",
        pp(power(12L, "metabolites", "camera_floor", 5L)), pp(power(12L, "metabolites", "camera", 5L)))),
    c("enrich_test_camera_floor_note_metabolites", sprintf("called in up to %s %%, as with", E("metabolites", "camera_floor"))),
    c("enrich_test_permutation_note_genes_proteins", sprintf("in %s %% of cases against", Arange("genes_proteins", "permutation"))),
    c("enrich_test_permutation_note_genes_proteins", sprintf("(the most variable features: %s %%)", Brange("genes_proteins", "permutation"))),
    c("enrich_test_permutation_note_genes_proteins", sprintf("and in %s %% on the main factor", Crange("genes_proteins", "permutation"))),
    c("enrich_test_permutation_note_genes_proteins", sprintf("found in %s %% of cases (the default test: %s %%), of 30 in %s %%",
        pp(power(12L, "genes_proteins", "permutation", 10L)), pp(power(12L, "genes_proteins", "camera", 10L)),
        pp(power(12L, "genes_proteins", "permutation", 30L)))),
    c("enrich_test_permutation_note_genes_proteins", sprintf("it found %s %% of such terms (10 to 100 features)",
        pp(power(6L, "genes_proteins", "permutation", c(10L, 30L, 100L))))),
    c("enrich_test_permutation_note_genes_proteins", sprintf("called in up to %s %%.", E("genes_proteins", "permutation"))),
    c("enrich_test_permutation_note_metabolites", sprintf("in %s %% of cases against", Arange("metabolites", "permutation"))),
    c("enrich_test_permutation_note_metabolites", sprintf("and in %s %% on the main factor", Crange("metabolites", "permutation"))),
    c("enrich_test_permutation_note_metabolites", sprintf("found in %s %% of cases (the default test: %s %%)",
        pp(power(12L, "metabolites", "permutation", 5L)), pp(power(12L, "metabolites", "camera", 5L)))),
    c("enrich_test_permutation_note_metabolites", sprintf("called in up to %s %% (the default test: %s %%)",
        E("metabolites", "permutation"), E("metabolites", "camera"))),
    c("enrich_test_gsea_note_genes_proteins", sprintf("in %s %% of cases (terms of 10 to 100", D("genes_proteins", "gsea"))),
    c("enrich_test_gsea_note_genes_proteins", sprintf("found in %s %% of cases, a figure", pp(power(12L, "genes_proteins", "gsea", 10L)))),
    c("enrich_test_gsea_note_metabolites", sprintf("in %s %% of cases (terms of 5 to 10", D("metabolites", "gsea"))),
    c("enrich_test_gsea_note_metabolites", sprintf("found in %s %% of cases (in result", pp(power(12L, "metabolites", "gsea", 5L)))))
for (kind in names(views)) for (test in names(tests)) {
    key6 <- paste0("enrich_test_few_samples_", test, "_", kind)
    base <- if (test == "gsea") c("C0", "C2") else c("C0", "C1", "C2")
    checks[[length(checks) + 1L]] <- c(key6, sprintf("in up to %s %% of cases", P(max(cells(6L, kind, "F_ind", test, base)$T))))
    if (test == "gsea") checks[[length(checks) + 1L]] <- c(key6, sprintf("(%s %% when the term's features vary together",
        P(max(cells(6L, kind, "F_ind", test, "C1")$T))))
    checks[[length(checks) + 1L]] <- c(key6, sprintf("most variable features: up to %s %%", Bmax(6L, kind, test)))
    checks[[length(checks) + 1L]] <- c(key6, sprintf("were found in %s %%.", pp(power(6L, kind, test,
        if (kind == "metabolites") 5L else 30L))))
}
found <- data.frame(key = vapply(checks, `[[`, character(1), 1L), figure = vapply(checks, `[[`, character(1), 2L),
    stringsAsFactors = FALSE)
found$text_present <- found$key %in% names(texts)
found$found <- vapply(seq_len(nrow(found)), function(i) {
    text <- texts[[found$key[[i]]]]
    if (grepl("^enrich_test_few_samples_", found$key[[i]])) text <- sprintf(text, 6L)
    grepl(found$figure[[i]], text, fixed = TRUE)
}, logical(1))
utils::write.table(found, file.path(out_dir, "texts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
verdict <- c(evidence_rows = nrow(regenerated), evidence_package_rows = nrow(package),
    evidence_equal = sum(regenerated$equal), text_figures = nrow(found), text_figures_found = sum(found$found))
writeLines(sprintf("%s\t%d", names(verdict), verdict), file.path(out_dir, "verdict.txt"))
print(verdict)
print(regenerated[!regenerated$equal, ])
print(found[!found$found, ])
stopifnot(nrow(regenerated) == nrow(package), all(regenerated$equal), all(found$found))
