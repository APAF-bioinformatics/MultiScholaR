# usage: Rscript compare-content.R <before_root> <after_root>
# Pairs files by path with random ids masked; RDS compared structurally with timestamps, paths and environments
# masked; xlsx sheet by sheet; text with timestamp lines masked; other files by byte equality.
args <- commandArgs(TRUE); b <- args[1]; a <- args[2]
norm <- function(p) gsub("[0-9a-f]{10,}", "ID", gsub("([A-Za-z0-9_.-]*?-)[0-9a-f]{6,}(/|$)", "\\1ID\\2", p))
lst <- function(root) { f <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE); setNames(f, norm(substring(f, nchar(root) + 2L))) }
B <- lst(b); A <- lst(a)
scrub <- function(x) {
    if (inherits(x, "POSIXt") || inherits(x, "Date")) return("<time>")
    if (is.environment(x)) return("<env>")
    if (is.function(x)) return("<fn>")
    if (is.character(x)) { x <- gsub("/home/\\S+|/tmp/\\S+", "<path>", x); x <- gsub("\\d{4}-\\d{2}-\\d{2}[ T]?\\d{0,2}:?\\d{0,2}:?\\d{0,2}", "<ts>", x); return(x) }
    if (isS4(x)) { for (s in methods::slotNames(x)) { v <- tryCatch(methods::slot(x, s), error = function(e) NULL); if (!is.null(v)) attr(x, paste0(".s.", s)) <- scrub(v) }; return(x) }
    if (is.list(x)) { x[] <- lapply(x, scrub); return(x) }
    x
}
cmp <- function(pb, pa) {
    ext <- tolower(tools::file_ext(pb))
    if (ext %in% c("rds")) { ok <- isTRUE(all.equal(scrub(readRDS(pb)), scrub(readRDS(pa)), check.attributes = TRUE)); return(if (ok) "identical (timestamps/paths masked)" else paste("DIFFERS:", paste(head(all.equal(scrub(readRDS(pb)), scrub(readRDS(pa))), 3), collapse = " | "))) }
    if (ext == "xlsx") { sb <- readxl::excel_sheets(pb); sa <- readxl::excel_sheets(pa); if (!identical(sb, sa)) return("DIFFERS: sheets"); for (s in sb) { if (!isTRUE(all.equal(as.data.frame(readxl::read_excel(pb, s)), as.data.frame(readxl::read_excel(pa, s))))) return(paste("DIFFERS: sheet", s)) }; return("identical (cell contents)") }
    tb <- tryCatch(readLines(pb, warn = FALSE), error = function(e) NULL); ta <- tryCatch(readLines(pa, warn = FALSE), error = function(e) NULL)
    if (!is.null(tb) && !is.null(ta) && ext %in% c("txt","tab","tsv","csv","json","ini","md")) { m <- function(t) gsub("\\d{4}-\\d{2}-\\d{2}[ T]?\\d{0,2}:?\\d{0,2}:?\\d{0,2}(\\.\\d+)?( [A-Z]{2,4})?", "<ts>", gsub("/home/\\S+|/tmp/\\S+", "<path>", t)); return(if (identical(m(tb), m(ta))) "identical (timestamps/paths masked)" else "DIFFERS") }
    if (identical(unname(tools::md5sum(pb)), unname(tools::md5sum(pa)))) "identical (bytes)" else "DIFFERS (bytes)"
}
keys <- intersect(names(B), names(A)); tally <- table(character())
res <- vapply(keys, function(k) cmp(B[[k]], A[[k]]), character(1))
byte_same <- vapply(keys, function(k) identical(unname(tools::md5sum(B[[k]])), unname(tools::md5sum(A[[k]]))), logical(1))
cat(sprintf("paired %d; byte-identical %d; content-identical after masking %d; differing %d; only-before %d; only-after %d\n", length(keys), sum(byte_same), sum(!byte_same & !grepl("^DIFFERS", res)), sum(grepl("^DIFFERS", res)), length(setdiff(names(B), names(A))), length(setdiff(names(A), names(B)))))
for (k in keys[!byte_same]) cat(sprintf("  %-95s %s\n", k, res[[k]]))
for (k in setdiff(names(B), names(A))) cat("  ONLY-BEFORE", k, "\n"); for (k in setdiff(names(A), names(B))) cat("  ONLY-AFTER", k, "\n")
