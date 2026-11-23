# KIT-189 E7 (7.4; re-check R01): the continuation path of sim-k189.R, proven in the lifecycle run L1.
#
# Usage, from the package root on a run-tests runner, after sim-k189.R's part l1 in the same job:
#   Rscript docs/qa/workflow-execution/kit/kit-189/e7-k189.R <l1 out dir> <out dir>
#
# sim-k189.R on B1 null L = 126 (cell 2), replicates 9011 to 9050, is run under the operating system's `timeout` with
# a limit from L1's logged start-up, B1 set-up and median replicate wall in cell 2 (about half the run's length), so it
# is killed with some replicates complete and some not; the replicates without a complete part then run as a
# continuation from the list built here; then the same 40 replicates run uninterrupted. Every row goes to k189-blind/
# (never uploaded: the lifecycle is blind); only the check result is written to <out dir>/checks-e7.tsv.

args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
l1_dir <- args[[1L]]
out_dir <- args[[2L]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sim <- file.path("docs", "qa", "workflow-execution", "kit", "kit-189", "sim-k189.R")
blind <- file.path("k189-blind", "e7")
dir.create(blind, showWarnings = FALSE, recursive = TRUE)
reps <- 9011:9050
clean <- function(x) substr(gsub("[\"\t\r\n]+", " ", x), 1L, 600L)

result <- tryCatch({
    timing <- utils::read.delim(file.path(l1_dir, "lifecycle-timing.tsv"))
    limit <- ceiling(timing$startup_seconds + timing$setup_B1 + 5 * timing$cell2_rep_wall_median)
    run <- function(label, list_file = NULL, limit_seconds = NULL) {
        dir <- file.path(blind, label)
        cmd_args <- c(sim, label, "2", "9011", "9050", dir, list_file)
        env <- c(paste0("K189_JOB=", label))
        t0 <- proc.time()[["elapsed"]]
        status <- if (is.null(limit_seconds)) {
            system2("env", c(env, "Rscript", cmd_args), stdout = file.path(blind, paste0(label, ".log")),
                stderr = file.path(blind, paste0(label, ".log")))
        } else {
            system2("timeout", c("--kill-after=10", paste0(limit_seconds, "s"), "env", env, "Rscript", cmd_args),
                stdout = file.path(blind, paste0(label, ".log")), stderr = file.path(blind, paste0(label, ".log")))
        }
        list(dir = dir, status = status, seconds = proc.time()[["elapsed"]] - t0)
    }
    killed <- run("e7k", limit_seconds = limit)
    rows_k <- file.path(killed$dir, "rows")
    complete_k <- list.files(rows_k, pattern = "^part-2-[0-9]+\\.tsv$")
    tmp_k <- list.files(rows_k, pattern = "\\.tmp$")
    gz_k <- list.files(rows_k, pattern = "\\.tsv\\.gz$")
    done_k <- as.integer(sub("^part-2-([0-9]+)\\.tsv$", "\\1", complete_k))
    missing <- setdiff(reps, done_k)
    list_file <- file.path(blind, "e7-continuation-list.txt")
    writeLines(as.character(missing), list_file)
    cont <- run("e7c", list_file = list_file)
    whole <- run("e7u")
    # the union: the killed run's complete parts (only names ending .tsv: a .tmp part is never read) and the
    # continuation's assembled rows
    readRows <- function(files) do.call(rbind, lapply(files, function(f) utils::read.delim(f, colClasses = "character", quote = "",
        na.strings = character())))
    union_rows <- rbind(readRows(file.path(rows_k, complete_k)),
        readRows(list.files(file.path(cont$dir, "rows"), pattern = "\\.tsv\\.gz$", full.names = TRUE)))
    whole_rows <- readRows(list.files(file.path(whole$dir, "rows"), pattern = "\\.tsv\\.gz$", full.names = TRUE))
    key <- function(x) x[order(as.integer(x$replicate), x$propagation, x$arm), setdiff(names(x), c("part", "seconds"))]
    per_rep <- table(factor(union_rows$replicate[union_rows$arm == "Q" & union_rows$propagation == "none"], levels = reps))
    a <- key(union_rows)
    b <- key(whole_rows)
    rownames(a) <- NULL
    rownames(b) <- NULL
    pass <- killed$status %in% c(124L, 137L) && length(complete_k) >= 1L && length(missing) >= 1L &&
        cont$status == 0L && whole$status == 0L && all(per_rep == 1L) && nrow(union_rows) == 14L * length(reps) &&
        identical(a, b)
    data.frame(id = "E7", pass = pass, detail = clean(sprintf(paste(
        "limit %d s; killed run exit %d after %.0f s: %d complete parts, %d .tmp parts left (not read), %d assembled files;",
        "continuation of %d replicates exit %d (%.0f s); uninterrupted exit %d (%.0f s); union holds each of the 40",
        "replicates once: %s; union rows %d; rows equal to the uninterrupted run's without part and seconds: %s")),
        limit, killed$status, killed$seconds, length(complete_k), length(tmp_k), length(gz_k), length(missing), cont$status,
        cont$seconds, whole$status, whole$seconds, all(per_rep == 1L), nrow(union_rows), identical(a, b))),
        stringsAsFactors = FALSE)
}, error = function(e) data.frame(id = "E7", pass = NA, detail = clean(paste("instrument:", conditionMessage(e))),
    stringsAsFactors = FALSE))
utils::write.table(result, file.path(out_dir, "checks-e7.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(sprintf("check E7: %s", result$pass))
writeLines(result$detail)
