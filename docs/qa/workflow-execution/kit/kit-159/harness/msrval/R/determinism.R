# Determinism: comparing two jobs' tables, and fingerprints of tables that never leave a job
# (KIT-159 section 7.11).

#' Fingerprint of a table that may not leave the job (D8)
#'
#' Rows, columns, per numeric column the sum, minimum, maximum and the 0.01, 0.5 and 0.99 quantiles
#' (17 significant digits), and the SHA-256 of the table with numbers rounded to 10 significant
#' digits, so a rounding difference and a real one can be told apart.
#'
#' @param table A data frame.
#' @return A list: `rows`, `cols`, `numeric` (one row per numeric column), `sha256_10`.
#' @export
valFingerprint <- function(table) {
    checkmate::assert_data_frame(table)
    num <- names(table)[vapply(table, is.numeric, logical(1))]
    stats17 <- do.call(rbind, lapply(num, function(col) {
        x <- table[[col]][is.finite(table[[col]])]
        q <- if (length(x)) stats::quantile(x, c(0.01, 0.5, 0.99), names = FALSE, type = 7) else
            rep(NA_real_, 3L)
        data.frame(column = col, sum = signif(sum(x), 17), min = signif(min(x, Inf), 17),
                   max = signif(max(x, -Inf), 17), q01 = signif(q[[1L]], 17),
                   q50 = signif(q[[2L]], 17), q99 = signif(q[[3L]], 17), stringsAsFactors = FALSE)
    }))
    rounded <- table
    for (col in num) rounded[[col]] <- ifelse(is.finite(table[[col]]),
                                              sprintf("%.9e", table[[col]]), "NA")
    text <- utils::capture.output(utils::write.table(rounded, sep = "\t", quote = FALSE,
                                                     row.names = FALSE))
    list(rows = nrow(table), cols = ncol(table),
         numeric = stats17 %||% data.frame(column = character(), sum = numeric(), min = numeric(),
                                           max = numeric(), q01 = numeric(), q50 = numeric(),
                                           q99 = numeric()),
         sha256_10 = digest::digest(paste(text, collapse = "\n"), algo = "sha256",
                                    serialize = FALSE))
}

#' Compare one table of two jobs aligned by key columns
#' @noRd
valCompareTable <- function(pa, pb, keys) {
    a <- as.data.frame(data.table::fread(pa, sep = "\t"))
    b <- as.data.frame(data.table::fread(pb, sep = "\t"))
    if (!identical(names(a), names(b))) {
        return(list(result = "differs", max_abs_diff = NA_real_, max_rel_diff = NA_real_,
                    cause = "columns differ"))
    }
    keyOf <- function(x) {
        if (length(keys)) do.call(paste, c(x[keys], sep = "\r")) else as.character(seq_len(nrow(x)))
    }
    ka <- keyOf(a)
    kb <- keyOf(b)
    common <- intersect(ka, kb)
    one_side <- length(setdiff(ka, kb)) + length(setdiff(kb, ka))
    ia <- match(common, ka)
    ib <- match(common, kb)
    num <- setdiff(names(a)[vapply(a, is.numeric, logical(1))], keys)
    disc <- setdiff(names(a), c(num, keys))
    disc_mismatch <- sum(vapply(disc, function(col) {
        sum(as.character(a[[col]][ia]) != as.character(b[[col]][ib]), na.rm = TRUE)
    }, numeric(1)))
    abs_d <- rel_d <- 0
    for (col in num) {
        x <- a[[col]][ia]
        y <- b[[col]][ib]
        ok <- is.finite(x) & is.finite(y)
        d <- abs(x[ok] - y[ok])
        den <- pmax(abs(x[ok]), abs(y[ok]))
        r <- ifelse(den == 0, 0, d / den)
        abs_d <- max(abs_d, d, 0)
        rel_d <- max(rel_d, r, 0)
        if (any(is.finite(x) != is.finite(y))) disc_mismatch <- disc_mismatch + 1
    }
    within <- one_side == 0 && disc_mismatch == 0 && rel_d <= 1e-12
    list(result = if (within) "within_rounding" else "differs", max_abs_diff = abs_d,
         max_rel_diff = rel_d,
         cause = sprintf("rows on one side %d; discrete mismatches %d", one_side, disc_mismatch))
}

#' Compare the result tables of two jobs (G8, section 7.11)
#'
#' Equal SHA-256: `identical`. Otherwise, aligned by key columns: rows on one side only, discrete
#' mismatches, the largest absolute and relative numeric differences; `within_rounding` when every
#' relative difference is at most 1e-12 and discrete columns agree; else `differs`.
#'
#' @param dir_a,dir_b The two jobs' output folders.
#' @param keys Named list: table file = key columns.
#' @return A data frame, one row per table (the rows of `determinism.tsv`).
#' @export
valCompareRuns <- function(dir_a, dir_b, keys = list()) {
    files <- union(list.files(dir_a, pattern = "\\.tsv$", recursive = TRUE),
                   list.files(dir_b, pattern = "\\.tsv$", recursive = TRUE))
    do.call(rbind, lapply(files, function(f) {
        pa <- file.path(dir_a, f)
        pb <- file.path(dir_b, f)
        sa <- if (file.exists(pa)) valSha256(pa) else NA_character_
        sb <- if (file.exists(pb)) valSha256(pb) else NA_character_
        cmp <- if (is.na(sa) || is.na(sb)) {
            list(result = "differs", max_abs_diff = NA_real_, max_rel_diff = NA_real_,
                 cause = "present in one job only")
        } else if (identical(sa, sb)) {
            list(result = "identical", max_abs_diff = 0, max_rel_diff = 0, cause = "")
        } else {
            valCompareTable(pa, pb, keys[[f]] %||% character())
        }
        data.frame(table = f, sha256_job1 = sa, sha256_job2 = sb,
                   identical = identical(cmp$result, "identical"), result = cmp$result,
                   max_abs_diff = cmp$max_abs_diff, max_rel_diff = cmp$max_rel_diff,
                   cause = cmp$cause, stringsAsFactors = FALSE)
    }))
}
