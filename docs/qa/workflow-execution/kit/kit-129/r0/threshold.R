# KIT-129 R0, part B: the threshold of the factor-tied-to-groups warning and its acceptance, per omic (ticket
# section 3, "Declaration, before any run"). Reads the association tables the simulators write in their association
# mode (s<n>-assoc.tsv, columns in association.R). Base R only.
#
#   Rscript threshold.R set   <S1|S2|S3|S4> <result dir> <assoc.tsv> [<assoc.tsv> ...]
#   Rscript threshold.R check <S1|S2|S3|S4> <result dir> <assoc.tsv> [<assoc.tsv> ...]
#
# Two separate invocations, and the split is enforced, not promised:
#   set    takes the setting set only. Before any statistic is read, the cell and run columns alone are read and the
#          invocation stops if a row lies outside the setting set's runs (so a check-set file handed to it is refused
#          unread). It stops if a threshold already exists in the result dir (one attempt).
#   check  takes the check set only (the same refusal, the other way round), reads the threshold `set` wrote, never
#          computes one, and stops if a result already exists (one attempt).
#
# The declaration (quoted from the ticket):
#   "Setting set: runs 1 to 1000 of C0, C3, C4 (correlation 0) and C5 (0.3), runs 1 to 500 of stress S3 (0.6) and S11
#    (3 per group, 0.3), for arms 3 or 4, 5, 7 and 8, per simulator. Check set: runs 2001 to 2500 of the same cells,
#    never looked at while setting. Threshold: the value on the setting set that fires in at most 10% of the
#    correlation-0 runs. Accepted, per omic, when on the check set it fires in at most 10% of correlation-0 runs and at
#    least 80% of the runs at 0.3 or more. One attempt."
# Readings made here (each is in the report for the user; none was chosen after a result):
#   - One threshold per omic, over the four arms and the three correlation-0 cells together (in S4 over both assays):
#     the declaration names one value and accepts "per omic". The per-cell and per-arm shares are written beside it.
#   - "Fires" = the largest factor association is above the threshold (strictly). The threshold is the smallest
#     observed value with at most 10% of the correlation-0 rows above it.
#   - A row without a W (the arm was no RUV on that data set, E = 0, or refused) cannot fire and stays in the
#     denominator: the package shows no warning without a W. The shares among rows with a W are written too.
#   - A crashed row is left out of both counts and reported.
#   - A set must be complete as declared (every cell, run and arm once), or there is no result.
# K129_MECHANICS=1 (environment): a mechanics check on a partial set. Completeness is not required, the files are
# prefixed "mechanics-" and every line says it is not a result. It cannot write or read a real threshold.

DECLARED <- list(
    set = list(C0 = 1:1000, C3 = 1:1000, C4 = 1:1000, C5 = 1:1000, S3 = 1:500, S11 = 1:500),
    check = list(C0 = 2001:2500, C3 = 2001:2500, C4 = 2001:2500, C5 = 2001:2500, S3 = 2001:2500, S11 = 2001:2500))
ZERO <- c("C0", "C3", "C4")                               # correlation 0
TIED <- c("C5", "S3", "S11")                              # correlation 0.3 or more
ARMS <- list(S1 = c("4", "5", "7", "8"), S2 = c("4", "5", "7", "8"), S3 = c("3", "5", "7", "8"), S4 = c("4", "5", "7", "8"))
OMIC <- c(S1 = "classic proteomics", S2 = "DIA_LIMPA", S3 = "transcriptomics", S4 = "lipidomics and metabolomics")
FIRE_LIMIT <- 0.10
DETECT_LIMIT <- 0.80

args <- commandArgs(TRUE)
if (length(args) < 4L || !args[[1L]] %in% c("set", "check") || !args[[2L]] %in% names(OMIC)) {
    stop("usage: Rscript threshold.R <set|check> <S1|S2|S3|S4> <result dir> <assoc.tsv> [<assoc.tsv> ...]")
}
phase <- args[[1L]]; simulator <- args[[2L]]; result_dir <- args[[3L]]; files <- args[-(1:3)]
mechanics <- nzchar(Sys.getenv("K129_MECHANICS", ""))
prefix <- if (mechanics) "mechanics-" else ""
label <- if (mechanics) "MECHANICS CHECK ON A PARTIAL SET, NOT A RESULT: " else ""
out <- function(name) file.path(result_dir, paste0(prefix, name, "-", simulator, ".tsv"))
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# The rows of one phase's set and no other. The cell and run columns are read first, alone.
readSet <- function(files, phase) {
    declared <- DECLARED[[phase]]
    set_name <- c(set = "setting", check = "check")[[phase]]
    for (file in files) {
        header <- names(utils::read.delim(file, nrows = 1L, check.names = FALSE))
        keys <- utils::read.delim(file, colClasses = ifelse(header %in% c("cell", "run"), NA, "NULL"), stringsAsFactors = FALSE)
        inside <- vapply(seq_len(nrow(keys)), function(i) keys$cell[[i]] %in% names(declared) && keys$run[[i]] %in% declared[[keys$cell[[i]]]], logical(1L))
        if (!all(inside)) {
            stop(sprintf("%s holds %d row(s) outside the %s set (e.g. %s run %s); only its cell and run columns were read. Refused.",
                file, sum(!inside), set_name, keys$cell[which(!inside)[[1L]]], keys$run[which(!inside)[[1L]]]))
        }
    }
    rows <- do.call(rbind, lapply(files, utils::read.delim, stringsAsFactors = FALSE,
        colClasses = c(arm = "character", assay = "character", factor_associations = "character")))
    stopifnot(all(rows$simulator == simulator), all(rows$arm %in% ARMS[[simulator]]))
    key <- paste(rows$cell, rows$run, rows$arm, rows$assay)
    if (anyDuplicated(key)) stop("a cell, run and arm occurs more than once: ", key[duplicated(key)][[1L]])
    assays <- unique(rows$assay)
    wanted <- unlist(lapply(names(declared), function(cell) as.vector(outer(as.vector(outer(paste(cell, declared[[cell]]),
        ARMS[[simulator]], paste)), assays, paste))))
    missing <- setdiff(wanted, key)
    if (length(missing) && !mechanics) {
        stop(sprintf("the %s set is not complete as declared: %d of %d rows are missing (e.g. %s). No result.", set_name,
            length(missing), length(wanted), missing[[1L]]))
    }
    rows$class <- ifelse(rows$cell %in% ZERO, "correlation 0", "correlation 0.3 or more")
    rows
}

# Shares of rows that fire, by scope. Crashed rows are counted and left out.
shares <- function(rows, threshold) {
    scopes <- list(class = rows$class, cell = paste(rows$class, rows$cell, sep = " | "),
        arm = paste(rows$class, paste("arm", rows$arm), sep = " | "),
        `cell and arm` = paste(rows$class, rows$cell, paste("arm", rows$arm), sep = " | "))
    do.call(rbind, lapply(names(scopes), function(scope) do.call(rbind, lapply(split(rows, scopes[[scope]]), function(part) {
        used <- part[part$status != "crashed", , drop = FALSE]
        with_w <- used[used$status == "ok", , drop = FALSE]
        fired <- sum(with_w$group_association > threshold)
        data.frame(scope = scope, class = part$class[[1L]], cell = if (scope %in% c("cell", "cell and arm")) part$cell[[1L]] else "all",
            arm = if (scope %in% c("arm", "cell and arm")) part$arm[[1L]] else "all", rows = nrow(used), rows_with_W = nrow(with_w),
            crashed = sum(part$status == "crashed"), fired = fired, share_fired = round(fired / nrow(used), 4),
            share_fired_of_rows_with_W = if (nrow(with_w)) round(fired / nrow(with_w), 4) else NA_real_,
            association_median = if (nrow(with_w)) round(stats::median(with_w$group_association), 4) else NA_real_,
            association_p90 = if (nrow(with_w)) round(stats::quantile(with_w$group_association, 0.9, names = FALSE), 4) else NA_real_,
            true_association_median = round(stats::median(used$true_association, na.rm = TRUE), 4), stringsAsFactors = FALSE)
    }))))
}
writeTsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (phase == "set") {
    if (file.exists(out("threshold"))) stop(out("threshold"), " exists: the threshold is set once (one attempt)")
    rows <- readSet(files, "set")
    zero <- rows[rows$class == "correlation 0" & rows$status != "crashed", , drop = FALSE]
    values <- sort(zero$group_association[zero$status == "ok"], decreasing = TRUE)
    allowed <- floor(FIRE_LIMIT * nrow(zero))
    # the smallest observed value with at most `allowed` rows above it; 0 when fewer rows than that have a W at all
    threshold <- if (length(values) > allowed) values[[allowed + 1L]] else 0
    writeTsv(data.frame(simulator = simulator, omic = OMIC[[simulator]], threshold = signif(threshold, 8L),
        rule = sprintf("%sfires when the largest factor association is above the threshold; at most %d of %d correlation-0 rows above it",
            label, allowed, nrow(zero)),
        correlation_0_rows = nrow(zero), correlation_0_rows_with_W = length(values), set_files = paste(basename(files), collapse = ","),
        stringsAsFactors = FALSE), out("threshold"))
    writeTsv(shares(rows, threshold), out("setting"))
    cat(sprintf("%s%s (%s): threshold = %s on the setting set (%d rows; fires in %.1f%% of the correlation-0 rows)\n", label,
        simulator, OMIC[[simulator]], format(signif(threshold, 4L)), nrow(rows), 100 * sum(values > threshold) / nrow(zero)))
} else {
    if (!file.exists(out("threshold"))) stop(out("threshold"), " does not exist: run `set` on the setting set first")
    if (file.exists(out("result"))) stop(out("result"), " exists: the check set is read once (one attempt)")
    threshold <- utils::read.delim(out("threshold"), stringsAsFactors = FALSE)$threshold
    rows <- readSet(files, "check")
    table <- shares(rows, threshold)
    by_class <- table[table$scope == "class", ]
    zero <- by_class$share_fired[by_class$class == "correlation 0"]
    tied <- by_class$share_fired[by_class$class == "correlation 0.3 or more"]
    accepted <- length(zero) == 1L && length(tied) == 1L && zero <= FIRE_LIMIT && tied >= DETECT_LIMIT
    verdict <- if (accepted) sprintf("threshold = %s, accepted", format(signif(threshold, 4L))) else "no verdict"
    line <- sprintf("%s%s (%s): %s (check set: fires in %.1f%% of correlation-0 rows, limit %.0f%%; in %.1f%% of rows at 0.3 or more, needed %.0f%%; threshold %s from the setting set)",
        label, simulator, OMIC[[simulator]], verdict, 100 * zero, 100 * FIRE_LIMIT, 100 * tied, 100 * DETECT_LIMIT,
        format(signif(threshold, 4L)))
    writeTsv(data.frame(simulator = simulator, omic = OMIC[[simulator]], threshold = threshold, share_fired_correlation_0 = zero,
        share_fired_correlation_0.3_or_more = tied, accepted = accepted, verdict = paste0(label, verdict), stringsAsFactors = FALSE,
        check.names = FALSE), out("result"))
    writeTsv(table, out("check"))
    cat(line, "\n", sep = "")
}
