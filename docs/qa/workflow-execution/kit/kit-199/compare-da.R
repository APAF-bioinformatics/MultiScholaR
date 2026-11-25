# KIT-199 script S3: compare two da-runner.tsv files (ticket section 6.2; R4, and under S6b R4's second script).
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-199/compare-da.R <a.tsv> <b.tsv> <out dir>
#
# Both files hold the DA table's keys (`Protein.Ids`, `comparison`) and its 20 numeric columns as written (control.R).
# Compared exactly: the headers, the keys, their order, the missing values and the calls at the app's q threshold
# (R/func_prot_da_model_methods.R:597-608: q < 0.05 and log2FC >= 0 is up, q < 0.05 and log2FC < 0 is down).
# Per numeric column: values compared, values whose text differs, the largest absolute and relative difference, and
# the values outside SWD-083's rule (packaging/verify/compare-results.R: |a - b| <= 1e-6 * max(|a|, |b|) or
# |a - b| <= 1e-9). Verdict per column and overall: `identical` (every text equal), `agree` (within the rule, keys,
# order, missing values and calls exact) or `outside`.
# Writes <out dir>/compare-da.tsv (quantity, column, value). A measurement: exit status 0 whenever the comparison is
# made; 1 when it cannot be made (a file unreadable, a header without the keys, duplicated keys).

REL <- 1e-6
ABS <- 1e-9
Q_THRESHOLD <- 0.05
KEYS <- c("Protein.Ids", "comparison")

readTable <- function(path) {
    lines <- readLines(path, warn = FALSE)
    fields <- lapply(strsplit(paste0(lines, "\tEND"), "\t", fixed = TRUE), \(row) row[-length(row)])
    header <- fields[[1L]]
    if (any(lengths(fields) != length(header))) stop(basename(path), ": rows with another field count")
    if (!all(KEYS %in% header)) stop(basename(path), ": the key columns are missing")
    body <- do.call(rbind, fields[-1L])
    colnames(body) <- header
    body
}

callsOf <- function(table) {
    q <- suppressWarnings(as.numeric(table[, "fdr_qvalue"]))
    l <- suppressWarnings(as.numeric(table[, "log2FC"]))
    called <- !is.na(q) & q < Q_THRESHOLD & !is.na(l)
    ifelse(called & l >= 0, "up", ifelse(called & l < 0, "down", "none"))
}

compareColumn <- function(x_text, y_text) {
    x <- suppressWarnings(as.numeric(x_text))
    y <- suppressWarnings(as.numeric(y_text))
    na_same <- identical(is.na(x), is.na(y))
    both <- !is.na(x) & !is.na(y)
    d <- abs(x[both] - y[both])
    m <- pmax(abs(x[both]), abs(y[both]))
    relative <- ifelse(m > 0, d / m, 0)
    ok <- d <= REL * m | d <= ABS
    text_differ <- sum(x_text != y_text)
    list(compared = sum(both), text_differ = text_differ, max_abs = if (length(d)) max(d) else 0,
        max_rel = if (length(d)) max(relative) else 0, outside = sum(!ok) + sum(xor(is.na(x), is.na(y))),
        na_same = na_same,
        verdict = if (text_differ == 0L) "identical" else if (na_same && all(ok)) "agree" else "outside")
}

compareMain <- function(a_path, b_path, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    rows <- list()
    emit <- function(quantity, column, value) {
        if (is.double(value)) value <- sprintf("%.17g", value)
        rows[[length(rows) + 1L]] <<- data.frame(quantity = quantity, column = column, value = as.character(value),
            stringsAsFactors = FALSE)
    }
    sha_a <- digest::digest(file = a_path, algo = "sha256")
    sha_b <- digest::digest(file = b_path, algo = "sha256")
    emit("sha256_a", "", sha_a)
    emit("sha256_b", "", sha_b)
    emit("bytes_identical", "", identical(sha_a, sha_b))
    a <- readTable(a_path)
    b <- readTable(b_path)
    header_identical <- identical(colnames(a), colnames(b))
    emit("header_identical", "", header_identical)
    emit("rows_a", "", nrow(a))
    emit("rows_b", "", nrow(b))
    key_a <- paste(a[, "Protein.Ids"], a[, "comparison"], sep = "\t")
    key_b <- paste(b[, "Protein.Ids"], b[, "comparison"], sep = "\t")
    if (anyDuplicated(key_a) || anyDuplicated(key_b)) stop("duplicated keys")
    common <- intersect(key_a, key_b)
    ia <- match(common, key_a)
    ib <- match(common, key_b)
    emit("keys_only_a", "", sum(!key_a %in% key_b))
    emit("keys_only_b", "", sum(!key_b %in% key_a))
    key_order_identical <- identical(key_a, key_b)
    emit("key_order_identical", "", key_order_identical)

    calls_a <- callsOf(a)
    calls_b <- callsOf(b)
    calls_changed <- sum(calls_a[ia] != calls_b[ib])
    for (side in c("a", "b")) {
        calls <- if (side == "a") calls_a else calls_b
        emit(paste0("called_up_", side), "", sum(calls == "up"))
        emit(paste0("called_down_", side), "", sum(calls == "down"))
    }
    emit("calls_changed", "", calls_changed)

    numeric_columns <- setdiff(intersect(colnames(a), colnames(b)), KEYS)
    verdicts <- character()
    na_all_same <- TRUE
    for (column in numeric_columns) {
        result <- compareColumn(a[ia, column], b[ib, column])
        emit("values_compared", column, result$compared)
        emit("values_text_differ", column, result$text_differ)
        emit("max_abs_diff", column, result$max_abs)
        emit("max_rel_diff", column, result$max_rel)
        emit("outside_rule", column, result$outside)
        emit("missing_values_same", column, result$na_same)
        emit("verdict", column, result$verdict)
        verdicts <- c(verdicts, result$verdict)
        na_all_same <- na_all_same && result$na_same
    }
    emit("numeric_columns", "", length(numeric_columns))
    exact_parts <- header_identical && key_order_identical && nrow(a) == nrow(b) && length(common) == nrow(a) &&
        calls_changed == 0L && na_all_same
    overall <- if (!exact_parts || any(verdicts == "outside")) "outside" else if (all(verdicts == "identical"))
        "identical" else "agree"
    emit("verdict_overall", "", overall)
    table <- do.call(rbind, rows)
    utils::write.table(table, file.path(out_dir, "compare-da.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    print(table, row.names = FALSE)
    0L
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("usage: Rscript compare-da.R <a.tsv> <b.tsv> <out dir>")
code <- tryCatch(compareMain(arguments[[1L]], arguments[[2L]], arguments[[3L]]), error = \(e) {
    cat("compare-da: the comparison could not be made:", conditionMessage(e), "\n")
    1L
})
quit(save = "no", status = code)
