# KIT-125 (ticket revision 1.4): what the four per-omic scripts share. Evidence only: nothing under R/, tests/ or
# tools/ is changed or sourced apart from pkgload::load_all(".").
#
# A per-omic script is run from the root of the pinned work tree the KIT-124 simulations used (0b6fcb6d), with R_LIBS
# as KIT-124's resume-20261118.sh sets it and one BLAS/OpenMP thread:
#   Rscript <this folder>/k125_<omic>.R <out dir> <data dir> <splits> [first split = 1]
# <data dir> is the omic's folder of public files outside the repository; <out dir> receives
#   runs.tsv    one row per split type, split and arm: status, k, n_features, n_calls, null_p05, null_p01, seconds
#   splits.tsv  the units of each false group (public sample names)
#   events.tsv  the message of every refusal and crash
#   meta.txt    commit, libraries, versions, command, seed rule, sha256 of the input files, the pool
# Environment, for the lifecycle test only: K125_FAIL_AT=<split number> makes the recommended design's arm stop with
# an error in that split, to prove that a failure is recorded and the run goes on.
#
# The null split. A pool of biological units of ONE condition. Per split, 12 units are drawn from the pool and put
# in two false groups of 6 (the group size of the simulations' centre cell). Every p under 0.05 and every call is
# false. Split types:
#   random  12 units drawn from the pool without regard to anything; the first 6 drawn are group A.
#   batch   (pools with two recorded batches only) group A = 4 units of batch 1 + 2 of batch 2, group B = 2 + 4:
#           the correlation of batch with the false group is 1/3 by layout (the simulations' "known unbalanced
#           batch" layout and close to cell C5's 0.3).
# Seeds. seed(split) = 5,000,000 + 100,000 * omic index + 10,000 * type index (random 0, batch 1) + split number.
# The split's units are drawn under set.seed(seed); arm a runs under set.seed(seed * 100 + a); the estimate of k
# uses the stream its simulator used. Splits of one invocation must be distinct (drawn without replacement): the
# script stops before any arm if two are the same partition.
# Statuses, as in the simulators: ok; refused = a designed refusal, scored as the same split's no-RUV result;
# crashed = any other error, measures NA.
# Measures per row: n_features = features with a finite p; null_p05, null_p01 = share of them with p under 0.05,
# 0.01; n_calls = features the omic's rule would call (q under 0.05). "Any call" = n_calls > 0.

OMIC_INDEX <- c(classic = 1L, dia_limpa = 2L, trx = 3L, metab = 4L)
TYPE_INDEX <- c(random = 0L, batch = 1L)
RUN_COLUMNS <- c("omic", "split_type", "split", "seed", "arm", "design", "status", "k", "n_features", "n_calls",
    "null_p05", "null_p01", "seconds")

k125Args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 3L) stop("usage: Rscript k125_<omic>.R <out dir> <data dir> <splits> [first split = 1]")
    list(out_dir = normalizePath(args[[1L]], mustWork = FALSE), data_dir = normalizePath(args[[2L]], mustWork = TRUE),
        splits = as.integer(args[[3L]]), first = if (length(args) >= 4L) as.integer(args[[4L]]) else 1L,
        fail_at = suppressWarnings(as.integer(Sys.getenv("K125_FAIL_AT", ""))))
}

splitSeed <- function(omic, type, split) 5000000L + 100000L * OMIC_INDEX[[omic]] + 10000L * TYPE_INDEX[[type]] + split

# The units of one split. `pool` = unit ids; `batch` = their batch (two levels), for type "batch".
drawSplit <- function(pool, batch, type, seed) {
    set.seed(seed)
    if (type == "random") {
        units <- sample(pool, 12L)
        return(list(A = units[1:6], B = units[7:12]))
    }
    levels <- sort(unique(batch))
    stopifnot(length(levels) == 2L)
    first <- sample(pool[batch == levels[[1L]]], 6L)
    second <- sample(pool[batch == levels[[2L]]], 6L)
    list(A = c(first[1:4], second[1:2]), B = c(first[5:6], second[3:6]))
}
splitKey <- function(s) paste(sort(c(paste(sort(s$A), collapse = ","), paste(sort(s$B), collapse = ","))), collapse = " | ")

nullMeasures <- function(p, q) {
    tested <- is.finite(p)
    data.frame(n_features = sum(tested), n_calls = sum(tested & is.finite(q) & q < 0.05),
        null_p05 = mean(p[tested] < 0.05), null_p01 = mean(p[tested] < 0.01))
}
NA_MEASURES <- data.frame(n_features = NA_integer_, n_calls = NA_integer_, null_p05 = NA_real_, null_p01 = NA_real_)

appendTsv <- function(x, path) {
    new <- !file.exists(path)
    utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = new, append = !new, na = "NA")
}

# The study. `prepare(split)` returns list(arms = named list of functions, each returning list(p, q, k) or
# list(same_as_no_ruv = TRUE)); names are "<arm number>" and `designs` maps them to labels. `is_refusal(e)` tells a
# designed refusal from a crash. `pool`, `batch` as drawSplit; `types` the split types this pool allows.
k125Run <- function(omic, args, pool, batch, types, designs, prepare, is_refusal, inputs, packages, notes = character()) {
    dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
    path <- function(name) file.path(args$out_dir, name)
    if (file.exists(path("runs.tsv"))) stop("the out dir already holds runs.tsv: one attempt per out dir")
    started <- Sys.time()
    numbers <- args$first - 1L + seq_len(args$splits)
    plan <- do.call(rbind, lapply(types, function(type) data.frame(type = type, split = numbers, stringsAsFactors = FALSE)))
    plan$seed <- mapply(function(type, split) splitSeed(omic, type, split), plan$type, plan$split)
    drawn <- lapply(seq_len(nrow(plan)), function(i) drawSplit(pool, batch, plan$type[[i]], plan$seed[[i]]))
    keys <- paste(plan$type, vapply(drawn, splitKey, ""))
    if (anyDuplicated(keys)) stop("two splits of this invocation are the same partition")
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
    git <- function(...) tryCatch(paste(system2("git", c(...), stdout = TRUE, stderr = FALSE), collapse = " "), error = function(e) NA_character_)
    sha <- function(file) unname(tools::sha256sum(file))
    meta <- c(
        sprintf("KIT-125 null split, omic %s (ticket revision 1.4)", omic),
        sprintf("package tree (working directory): %s; commit %s; changes under R/, tests/, src/: %s", getwd(), git("rev-parse", "HEAD"),
            nzchar(git("status", "--porcelain", "--untracked-files=no", "--", "R", "tests", "src"))),
        sprintf("script: %s; sha256 %s; common sha256 %s", script, sha(script), sha(file.path(dirname(script), "k125_common.R"))),
        sprintf("command: Rscript %s", paste(c(script, commandArgs(TRUE)), collapse = " ")),
        sprintf("K125_FAIL_AT: %s", if (is.na(args$fail_at)) "not set" else args$fail_at),
        sprintf("R_LIBS: %s", Sys.getenv("R_LIBS")), sprintf("library paths: %s", paste(.libPaths(), collapse = ":")),
        sprintf("threads: OMP_NUM_THREADS=%s OPENBLAS_NUM_THREADS=%s; one R process", Sys.getenv("OMP_NUM_THREADS"), Sys.getenv("OPENBLAS_NUM_THREADS")),
        sprintf("%s; %s; BLAS %s", R.version.string, R.version$platform, extSoftVersion()[["BLAS"]]),
        sprintf("packages: %s", paste(vapply(packages, function(p) paste(p, as.character(utils::packageVersion(p))), ""), collapse = "; ")),
        sprintf("input files (sha256): %s", paste(sprintf("%s %s", basename(inputs), vapply(inputs, sha, "")), collapse = "; ")),
        sprintf("pool: %d units%s", length(pool), if (is.null(batch)) "; no batch recorded" else
            sprintf("; batches %s", paste(names(table(batch)), as.integer(table(batch)), sep = " = ", collapse = ", "))),
        sprintf("split types: %s; splits %d to %d of each", paste(types, collapse = ", "), min(numbers), max(numbers)),
        "seed rule: seed = 5000000 + 100000 * omic index + 10000 * type index + split; units under set.seed(seed); arm a under set.seed(seed * 100 + a)",
        notes, sprintf("started: %s", format(started, "%Y-%m-%d %H:%M:%S %Z")))
    writeLines(meta, path("meta.txt"))
    for (i in seq_len(nrow(plan))) {
        type <- plan$type[[i]]; split <- plan$split[[i]]; seed <- plan$seed[[i]]; units <- drawn[[i]]
        appendTsv(data.frame(omic = omic, split_type = type, split = split, seed = seed, group_A = paste(units$A, collapse = ","),
            group_B = paste(units$B, collapse = ",")), path("splits.tsv"))
        event <- function(arm, status, e) appendTsv(data.frame(omic = omic, split_type = type, split = split, arm = arm, status = status,
            message = gsub("[\t\r\n]+", " ", conditionMessage(e))), path("events.tsv"))
        row <- function(arm, status, k, m, seconds) data.frame(omic = omic, split_type = type, split = split, seed = seed, arm = arm,
            design = designs[[arm]], status = status, k = k, m, seconds = round(seconds, 3), stringsAsFactors = FALSE)[RUN_COLUMNS]
        t0 <- proc.time()[["elapsed"]]
        prepared <- tryCatch(prepare(units, seed), error = function(e) e)
        if (inherits(prepared, "error")) {
            event("all", "crashed", prepared)
            appendTsv(do.call(rbind, lapply(names(designs), function(arm) row(arm, "crashed", NA, NA_MEASURES, NA_real_))), path("runs.tsv"))
            next
        }
        prepare_seconds <- proc.time()[["elapsed"]] - t0
        rows <- list(); reference <- NULL
        for (arm in names(designs)) {
            set.seed(seed * 100L + as.integer(arm))
            t0 <- proc.time()[["elapsed"]]
            result <- tryCatch({
                if (arm == "7" && !is.na(args$fail_at) && split == args$fail_at) stop("K125_FAIL_AT: the deliberate failure of the lifecycle test")
                prepared$arms[[arm]]()
            }, error = function(e) e)
            seconds <- proc.time()[["elapsed"]] - t0
            if (inherits(result, "error")) {
                refused <- is_refusal(result) && !is.null(reference) && arm != "1"
                event(arm, if (refused) "refused" else "crashed", result)
                rows[[arm]] <- row(arm, if (refused) "refused" else "crashed", NA, if (refused) reference else NA_MEASURES, seconds)
                next
            }
            if (isTRUE(result$same_as_no_ruv)) {
                rows[[arm]] <- row(arm, "ok", 0L, reference, seconds)
                next
            }
            m <- nullMeasures(result$p, result$q)
            if (arm == "1") reference <- m
            rows[[arm]] <- row(arm, "ok", result$k, m, seconds)
        }
        rows <- do.call(rbind, rows)
        rows$seconds[rows$arm == "1"] <- rows$seconds[rows$arm == "1"] + round(prepare_seconds, 3)   # the shared preparation is in arm 1's time
        appendTsv(rows, path("runs.tsv"))
        cat(sprintf("%s %s %s split %d done (%.1f s)\n", format(Sys.time(), "%H:%M:%S"), omic, type, split, sum(rows$seconds, na.rm = TRUE)))
    }
    cat(sprintf("finished: %s; wall %.0f s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), as.numeric(difftime(Sys.time(), started, units = "secs"))),
        file = path("meta.txt"), append = TRUE)
}
