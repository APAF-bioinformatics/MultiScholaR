# KIT-108 W5 (D8 content checks): results_summary files identical W4 vs W5, xlsx sheet by sheet; report text identical
# with the timestamp masked. Args: <w4 work dir> <w5 work dir> <out tsv>
a <- commandArgs(TRUE)
w4 <- normalizePath(a[1]); w5 <- normalizePath(a[2]); out <- a[3]
summaryDir <- function(work) {
    d <- list.dirs(work, recursive = TRUE)
    d[basename(dirname(d)) == "results_summary" & !grepl("_backup_[0-9]{8}_[0-9]{6}$", d)]
}
s4 <- summaryDir(w4); s5 <- summaryDir(w5)
stopifnot(length(s4) == 1L, length(s5) == 1L, basename(s4) == basename(s5))
rel <- function(d) sort(list.files(d, recursive = TRUE, all.files = TRUE))
f4 <- rel(s4); f5 <- rel(s5)
docText <- function(path, root) {
    txt <- system2(rmarkdown::pandoc_exec(), c("-t", "plain", "--wrap=none", shQuote(path)), stdout = TRUE)
    # The replica's own work directory (the W4 or W5 copy) appears in the path warnings; it is masked too.
    txt <- gsub(root, "<WORK>", txt, fixed = TRUE)
    # The report's generation timestamp (date and time) is masked.
    gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}[ T][0-9]{2}:[0-9]{2}(:[0-9]{2})?", "<TIMESTAMP>", txt)
}
hashOne <- function(path, root) {
    lower <- tolower(path)
    if (grepl("[.]xlsx$", lower)) {
        sheets <- openxlsx::getSheetNames(path)
        h <- vapply(sheets, function(s) digest::digest(openxlsx::read.xlsx(path, sheet = s), algo = "sha256"), "")
        return(data.frame(unit = paste0("sheet:", sheets), sha256 = unname(h), method = "sha256_xlsx_sheet"))
    }
    if (grepl("[.]docx$", lower)) {
        txt <- docText(path, root)
        return(data.frame(unit = "text", sha256 = digest::digest(txt, algo = "sha256"), method = "sha256_pandoc_plain_timestamp_masked",
            lines = length(txt), masked = sum(grepl("<TIMESTAMP>", txt))))
    }
    data.frame(unit = "file", sha256 = digest::digest(file = path, algo = "sha256"), method = "sha256")
}
rows <- list()
for (f in union(f4, f5)) {
    h4 <- if (f %in% f4) hashOne(file.path(s4, f), w4) else NULL
    h5 <- if (f %in% f5) hashOne(file.path(s5, f), w5) else NULL
    if (is.null(h4) || is.null(h5)) {
        rows[[length(rows) + 1L]] <- data.frame(file = f, unit = NA, method = NA, w4 = !is.null(h4), w5 = !is.null(h5), identical = FALSE)
        next
    }
    units <- union(h4$unit, h5$unit)
    rows[[length(rows) + 1L]] <- data.frame(file = f, unit = units, method = h4$method[1],
        w4 = h4$sha256[match(units, h4$unit)], w5 = h5$sha256[match(units, h5$unit)],
        identical = h4$sha256[match(units, h4$unit)] %in% h5$sha256[match(units, h5$unit)] &
            !is.na(match(units, h4$unit)) & !is.na(match(units, h5$unit)))
}
res <- do.call(rbind, lapply(rows, function(r) { r$w4 <- as.character(r$w4); r$w5 <- as.character(r$w5); r }))
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d files, %d units, %d identical, %d differ; sheets %d; report units %d\n", basename(out), length(union(f4, f5)), nrow(res),
    sum(res$identical), sum(!res$identical), sum(res$method == "sha256_xlsx_sheet", na.rm = TRUE), sum(res$method == "sha256_pandoc_plain_timestamp_masked", na.rm = TRUE)))
if (any(!res$identical)) print(res[!res$identical, c("file", "unit", "method")])
