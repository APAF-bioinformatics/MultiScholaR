# KIT-129 R0, part B: the amended rule of the factor-tied-to-groups warning, one threshold per design (arm) and omic
# (decision M23, 2026-11-18; the amendment is recorded in kit-129/r0/amendment-20261118.md on kit-129-r0b-record,
# committed before this script computed anything). A sibling of threshold.R, which stays as it was and still scores
# the original pooled rule. Base R only.
#
#   Rscript threshold-by-design.R set   <S1|S2|S3|S4> <result dir> <assoc.tsv> [<assoc.tsv> ...]
#   Rscript threshold-by-design.R check <S1|S2|S3|S4> <result dir> <assoc.tsv> [<assoc.tsv> ...]
#
# The same separation and refusals as threshold.R:
#   set    takes the setting set only: the cell and run columns alone are read first and the invocation stops if a row
#          lies outside the setting set's runs. It stops if the per-design thresholds already exist (one attempt).
#   check  takes the check set only (the same refusal, the other way round), reads the thresholds `set` wrote, never
#          computes one, and stops if a result already exists (one attempt).
#
# The amended rule: per design and omic, the threshold is the value on the setting set that fires in at most 10% of
# that design's correlation-0 runs; accepted, per design and omic, when on the check set it fires in at most 10% of
# that design's correlation-0 runs and at least 80% of its runs at 0.3 or more. One attempt; otherwise "no verdict".
# Conventions as in threshold.R: fires = strictly above; the threshold is the smallest observed value with at most
# 10% of the design's correlation-0 rows above it (0 when fewer rows than that have a W); a row without a W cannot
# fire and stays in the denominator; a crashed row is left out of both counts and reported; cells classed by their
# declared correlation; S4 over both assays; a set must be complete as declared, or there is no result.
# Files: threshold-by-design-S<n>.tsv, setting-by-design-S<n>.tsv (set); result-by-design-S<n>.tsv,
# check-by-design-S<n>.tsv (check).

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
    stop("usage: Rscript threshold-by-design.R <set|check> <S1|S2|S3|S4> <result dir> <assoc.tsv> [<assoc.tsv> ...]")
}
phase <- args[[1L]]; simulator <- args[[2L]]; result_dir <- args[[3L]]; files <- args[-(1:3)]
out <- function(name) file.path(result_dir, paste0(name, "-by-design-", simulator, ".tsv"))
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
    if (length(missing)) {
        stop(sprintf("the %s set is not complete as declared: %d of %d rows are missing (e.g. %s). No result.", set_name,
            length(missing), length(wanted), missing[[1L]]))
    }
    rows$class <- ifelse(rows$cell %in% ZERO, "correlation 0", "correlation 0.3 or more")
    rows
}

# Shares of one design's rows that fire under its threshold, by class and by cell.
sharesOfArm <- function(rows, arm, threshold) {
    rows <- rows[rows$arm == arm, , drop = FALSE]
    scopes <- list(class = rows$class, cell = paste(rows$class, rows$cell, sep = " | "))
    do.call(rbind, lapply(names(scopes), function(scope) do.call(rbind, lapply(split(rows, scopes[[scope]]), function(part) {
        used <- part[part$status != "crashed", , drop = FALSE]
        with_w <- used[used$status == "ok", , drop = FALSE]
        fired <- sum(with_w$group_association > threshold)
        data.frame(arm = arm, threshold = threshold, scope = scope, class = part$class[[1L]],
            cell = if (scope == "cell") part$cell[[1L]] else "all", rows = nrow(used), rows_with_W = nrow(with_w),
            crashed = sum(part$status == "crashed"), fired = fired, share_fired = round(fired / nrow(used), 4),
            share_fired_of_rows_with_W = if (nrow(with_w)) round(fired / nrow(with_w), 4) else NA_real_,
            association_median = if (nrow(with_w)) round(stats::median(with_w$group_association), 4) else NA_real_,
            association_p90 = if (nrow(with_w)) round(stats::quantile(with_w$group_association, 0.9, names = FALSE), 4) else NA_real_,
            true_association_median = round(stats::median(used$true_association, na.rm = TRUE), 4), stringsAsFactors = FALSE)
    }))))
}
writeTsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (phase == "set") {
    if (file.exists(out("threshold"))) stop(out("threshold"), " exists: the thresholds are set once (one attempt)")
    rows <- readSet(files, "set")
    thresholds <- do.call(rbind, lapply(ARMS[[simulator]], function(arm) {
        zero <- rows[rows$arm == arm & rows$class == "correlation 0" & rows$status != "crashed", , drop = FALSE]
        values <- sort(zero$group_association[zero$status == "ok"], decreasing = TRUE)
        allowed <- floor(FIRE_LIMIT * nrow(zero))
        threshold <- if (length(values) > allowed) values[[allowed + 1L]] else 0
        data.frame(simulator = simulator, omic = OMIC[[simulator]], arm = arm, threshold = signif(threshold, 8L),
            rule = sprintf("fires when the largest factor association is above the threshold; at most %d of %d correlation-0 rows of the design above it",
                allowed, nrow(zero)),
            correlation_0_rows = nrow(zero), correlation_0_rows_with_W = length(values),
            set_files = paste(basename(files), collapse = ","), stringsAsFactors = FALSE)
    }))
    writeTsv(thresholds, out("threshold"))
    table <- do.call(rbind, lapply(seq_len(nrow(thresholds)), function(i) sharesOfArm(rows, thresholds$arm[[i]], thresholds$threshold[[i]])))
    writeTsv(table, out("setting"))
    by_class <- table[table$scope == "class", ]
    for (i in seq_len(nrow(thresholds))) {
        part <- by_class[by_class$arm == thresholds$arm[[i]], ]
        cat(sprintf("%s (%s) arm %s: threshold = %s on the setting set (fires in %.1f%% of its correlation-0 rows, %.1f%% of its rows at 0.3 or more)\n",
            simulator, OMIC[[simulator]], thresholds$arm[[i]], format(signif(thresholds$threshold[[i]], 4L)),
            100 * part$share_fired[part$class == "correlation 0"], 100 * part$share_fired[part$class == "correlation 0.3 or more"]))
    }
} else {
    if (!file.exists(out("threshold"))) stop(out("threshold"), " does not exist: run `set` on the setting set first")
    if (file.exists(out("result"))) stop(out("result"), " exists: the check set is read once (one attempt)")
    thresholds <- utils::read.delim(out("threshold"), stringsAsFactors = FALSE, colClasses = c(arm = "character"))
    stopifnot(identical(thresholds$arm, ARMS[[simulator]]))
    rows <- readSet(files, "check")
    table <- do.call(rbind, lapply(seq_len(nrow(thresholds)), function(i) sharesOfArm(rows, thresholds$arm[[i]], thresholds$threshold[[i]])))
    by_class <- table[table$scope == "class", ]
    result <- do.call(rbind, lapply(seq_len(nrow(thresholds)), function(i) {
        part <- by_class[by_class$arm == thresholds$arm[[i]], ]
        zero <- part$share_fired[part$class == "correlation 0"]
        tied <- part$share_fired[part$class == "correlation 0.3 or more"]
        accepted <- length(zero) == 1L && length(tied) == 1L && zero <= FIRE_LIMIT && tied >= DETECT_LIMIT
        verdict <- if (accepted) sprintf("threshold = %s, accepted", format(signif(thresholds$threshold[[i]], 4L))) else "no verdict"
        cat(sprintf("%s (%s) arm %s: %s (check set: fires in %.1f%% of its correlation-0 rows, limit %.0f%%; in %.1f%% of its rows at 0.3 or more, needed %.0f%%; threshold %s from the setting set)\n",
            simulator, OMIC[[simulator]], thresholds$arm[[i]], verdict, 100 * zero, 100 * FIRE_LIMIT, 100 * tied, 100 * DETECT_LIMIT,
            format(signif(thresholds$threshold[[i]], 4L))))
        data.frame(simulator = simulator, omic = OMIC[[simulator]], arm = thresholds$arm[[i]], threshold = thresholds$threshold[[i]],
            share_fired_correlation_0 = zero, share_fired_correlation_0.3_or_more = tied, accepted = accepted, verdict = verdict,
            stringsAsFactors = FALSE, check.names = FALSE)
    }))
    writeTsv(result, out("result"))
    writeTsv(table, out("check"))
}
