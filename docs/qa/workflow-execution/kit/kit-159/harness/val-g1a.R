# KIT-159 G1(a) entry script (section 2; KIT-199 12.3, D7 / VAL-120 (2); VAL-133 (2)). From the
# package root, in a `run-tests` scripts request at the commit under test (full git history):
#
#   Rscript docs/qa/workflow-execution/kit/kit-159/harness/val-g1a.R <out dir>
#       [--lanes <l1,l2,...>] [--parts 1|1,2] [--reference <g1a-reference.tsv>]
#       [--da-link <evidence folder>] [--not-refrozen <lane=KIT-n,...>]
#
# Writes <out>/environment.tsv first (the job's CPU model, image, the ten D7 keys, the lockfile's
# SHA-256), then part (1) (every lane twice, installed against load_all, each in a child process:
# <out>/part1/g1a-part1.tsv and the lane records), then with --parts 1,2 part (2) against the
# oracles of HEAD through each lane's re-freeze chain (<out>/g1a-part2.tsv; the reference rows are
# required; prot_dia_limpa's DA table by value against the --da-link folder, default KIT-199's),
# and last <out>/gates.md. A refusal prints one template line with its class and exits 3.

harness_dir <- "docs/qa/workflow-execution/kit/kit-159/harness"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: val-g1a.R <out dir> [options]", call. = FALSE)
out <- args[[1L]]
opts <- args[-1L]
optValue <- function(name, default = NULL) {
    i <- match(name, opts)
    if (is.na(i) || i == length(opts)) default else opts[[i + 1L]]
}
tmp_root <- Sys.getenv("RUNNER_TEMP", tempdir())
val_lib <- file.path(tmp_root, "val-lib")
if (!requireNamespace("msrval", lib.loc = val_lib, quietly = TRUE)) {
    dir.create(val_lib, recursive = TRUE, showWarnings = FALSE)
    status <- system2(file.path(R.home("bin"), "R"),
                      c("CMD", "INSTALL", paste0("--library=", val_lib),
                        file.path(harness_dir, "msrval")),
                      stdout = FALSE, stderr = FALSE)
    if (!identical(as.integer(status), 0L)) stop("msrval did not install", call. = FALSE)
}
.libPaths(c(val_lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

refused <- function(e) {
    msrval::valStatus(phase = "g1a", status = "refused", class = class(e)[[1L]])
    quit(save = "no", status = 3L)
}
refuse <- function(cls) {
    refused(structure(class = c(cls, "val_error", "error", "condition"),
                      list(message = cls, call = NULL)))
}

dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
at <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
env <- msrval::valJobEnvironment(".")
msrval::valWriteJobEnvironment(env, file.path(out, "environment.tsv"))

all_lanes <- c("prot_dia", "prot_dia_limpa", "prot_lfq", "prot_lfq_fragpipe", "prot_tmt",
               "trx_screen", "lipid_canonical", "metab_lc", "metab_gc", "metab_combined")
lanes <- strsplit(optValue("--lanes", paste(all_lanes, collapse = ",")), ",")[[1L]]
parts <- strsplit(optValue("--parts", "1,2"), ",")[[1L]]
reference <- optValue("--reference")
if (!all(lanes %in% all_lanes) || !all(parts %in% c("1", "2"))) refuse("val_spec_invalid")
if ("2" %in% parts && (is.null(reference) || !file.exists(reference))) refuse("val_spec_invalid")

tryCatch({
    part1 <- NULL
    part2 <- NULL
    da <- NULL
    p1_dir <- file.path(out, "part1")
    if ("1" %in% parts || "2" %in% parts) {
        part1 <- msrval::valG1aPart1(".", lanes, p1_dir)
        msrval::valStatus(phase = "g1a", status = "part1",
                          count = sum(part1$status == "identical"))
    }
    if ("2" %in% parts) {
        read_steps <- function(lane) {
            p <- file.path(p1_dir, "installed", paste0(lane, ".json"))
            if (!file.exists(p)) return(NULL)
            jsonlite::read_json(p, simplifyVector = FALSE)$steps
        }
        live <- stats::setNames(lapply(lanes, read_steps), lanes)
        if ("prot_dia_limpa" %in% lanes) {
            da <- msrval::valG1aDaInputs(".", optValue("--da-link",
                                                       "docs/qa/workflow-execution/kit/kit-199"),
                                         env, file.path(p1_dir, "da-live.tsv"))
        }
        nr <- optValue("--not-refrozen")
        not_refrozen <- if (is.null(nr)) NULL else {
            kv <- strsplit(strsplit(nr, ",")[[1L]], "=", fixed = TRUE)
            stats::setNames(vapply(kv, `[[`, "", 2L), vapply(kv, `[[`, "", 1L))
        }
        part2 <- msrval::valG1aPart2(".", at, lanes, live, msrval::valReadReferenceRows(reference),
                                     msrval::valLaneComparator("."), da = da,
                                     not_refrozen = not_refrozen)
        utils::write.table(part2, file.path(out, "g1a-part2.tsv"), sep = "\t", quote = FALSE,
                           row.names = FALSE, na = "")
    }
    msrval::valG1aGatesMd(file.path(out, "gates.md"), at, env, part1, part2, da)
}, val_error = refused)
msrval::valStatus(phase = "g1a", status = "completed", count = length(lanes))
