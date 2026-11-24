# KIT-159 G1(r) script (section 2, G1(r); re-audit 2, R2.6 item 1). From the package root, in a
# `run-tests` scripts request at the commit under test (full git history; no product code runs):
#
#   Rscript docs/qa/workflow-execution/kit/kit-159/harness/val-g1r.R <out dir>
#
# Reads the committed source blocks (kit-159/g1a-reference-source.txt, checked against the SHA-256
# in g1a-reference-source-prov.txt), classes every difference with each lane's step order at the
# reference point 12d68ca7 and KIT-199's declared rows, checks the counts, and compares G1(r)'s own
# rows with every reference file a re-freezing ticket committed from the same source (found by the
# source's SHA-256, `git grep`). Writes <out>/g1a-reference.tsv (G1(r)'s own rows), g1r-diffs.tsv
# (lane, step, field, key, class; no values), g1r-checks.tsv, g1r-adoption.tsv,
# g1r-disagreements.tsv and g1r-status.txt (`built`, or `stop` with the reason). Exit 0 either way;
# a block the reader cannot read or a difference without a class writes `stop` with its class.

harness_dir <- "docs/qa/workflow-execution/kit/kit-159/harness"
kit_dir <- "docs/qa/workflow-execution/kit"
source_sha <- "1b918062f8c94a3c10859a965f5f1120ea3c692ae70181039846c66a43730bde"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: val-g1r.R <out dir>", call. = FALSE)
out <- args[[1L]]
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
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)

writeTsv <- function(x, name) {
    utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE,
                       na = "")
}
status <- function(state, reason) {
    writeLines(c(paste("state", state, sep = "\t"), paste("reason", reason, sep = "\t"),
                 paste("commit", system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                       sep = "\t")),
               file.path(out, "g1r-status.txt"))
    cat(sprintf("[val-g1r] %s: %s\n", state, reason))
}
gitShow <- function(at, path) {
    system2("git", c("show", paste0(at, ":", path)), stdout = TRUE)
}

checks <- data.frame(check = character(), value = character(), expected = character(),
                     pass = logical(), stringsAsFactors = FALSE)
check <- function(name, value, expected) {
    checks[nrow(checks) + 1L, ] <<- list(name, as.character(value), as.character(expected),
                                         identical(as.character(value), as.character(expected)))
}

source_file <- file.path(kit_dir, "kit-159", "g1a-reference-source.txt")
prov <- msrval::valReadTextTable(paste(c("key\tvalue", readLines(
    file.path(kit_dir, "kit-159", "g1a-reference-source-prov.txt"))), collapse = "\n"))
prov <- stats::setNames(prov$value, prov$key)
check("source_sha256_file", digest::digest(file = source_file, algo = "sha256"), prov[["sha256"]])
check("source_sha256_run_file", prov[["source_sha256"]], source_sha)

result <- tryCatch({
    blocks <- msrval::valReadG1rBlocks(source_file)
    files <- vapply(blocks, `[[`, "", "file")
    pairs <- paste(vapply(blocks, `[[`, "", "lane"), vapply(blocks, `[[`, "", "step"))
    check("blocks", length(blocks), 69L)
    per_file <- table(files)
    for (f in c("lipidomics", "metabolomics", "proteomics", "transcriptomics")) {
        check(paste0("blocks_", f), per_file[[paste0("test-gui-entry-", f, ".R")]],
              c(lipidomics = 4L, metabolomics = 30L, proteomics = 22L,
                transcriptomics = 13L)[[f]])
    }
    check("lane_step_pairs", length(unique(pairs)), 69L)
    diffs <- msrval::valG1rDiffs(blocks)
    check("blocks_without_difference",
          sum(table(factor(paste(diffs$lane, diffs$step), levels = unique(pairs))) == 0L), 0L)
    lanes <- unique(diffs$lane)
    since <- "12d68ca78442a50b46f52e33525dd85cf86df5e7"
    lane_steps <- stats::setNames(lapply(lanes, function(lane) {
        text <- gitShow(since, sprintf("tests/testdata/gui-entry/%s.json", lane))
        names(jsonlite::fromJSON(paste(text, collapse = "\n"), simplifyVector = FALSE)$steps)
    }), lanes)
    declared_path <- file.path(kit_dir, "kit-199", "declared-keys.tsv")
    declared <- msrval::valReadTextTable(declared_path)
    classed <- msrval::valClassifyLiveDiffs(diffs, lane_steps, declared = declared)
    writeTsv(classed[, c("lane", "step", "field", "key", "class")], "g1r-diffs.tsv")
    if (anyNA(classed$class)) {
        writeTsv(classed[is.na(classed$class), c("lane", "step", "field", "key")],
                 "g1r-unclassed.tsv")
        stop(structure(class = c("val_g1r_no_class", "val_error", "error", "condition"),
                       list(message = "a difference without a class of G1(a)", call = NULL)))
    }
    dec <- classed[classed$class %in% "declared_change (KIT-199)", , drop = FALSE]
    check("declared_rows_values", msrval:::valDeclaredValuesOk(dec, declared), TRUE)
    rows <- msrval::valG1rRows(classed, source = paste(
        "docs/qa/workflow-execution/kit/kit-159/g1a-reference-source.txt (the 69 GUI-entry",
        "blocks of run 44953187101's failures.txt)"), source_sha256 = source_sha)
    writeTsv(rows, "g1a-reference.tsv")
    check("rows", nrow(rows), 73L)
    per_class <- table(rows$class)
    check("rows_per_class", paste(sprintf("%s=%d", names(per_class), as.integer(per_class)),
                                  collapse = ","),
          paste("declared_change (KIT-199)=6", "integration_artifacts=24",
                "kit128_technical_replicates=13", "metab_fixture_4=30", sep = ","))
    for (lane in lanes) {
        check(paste0("counts_", lane),
              paste(sprintf("%s=%d", names(msrval:::valStepClassCounts(
                  classed[classed$lane == lane, ])),
                  msrval:::valStepClassCounts(classed[classed$lane == lane, ])), collapse = ","),
              paste(sprintf("%s=%d", names(msrval:::valClassCounts(rows$class[rows$lane == lane])),
                            msrval:::valClassCounts(rows$class[rows$lane == lane])),
                    collapse = ","))
    }
    found <- system2("git", c("grep", "-l", source_sha, "--", paste0(kit_dir, "/")),
                     stdout = TRUE)
    earlier <- found[grepl("[.]tsv$", found) & !startsWith(found, file.path(kit_dir, "kit-159"))]
    adoption <- data.frame(file = character(), groups = integer(), disagreements = integer(),
                           stringsAsFactors = FALSE)
    dis_all <- NULL
    for (f in earlier) {
        e <- msrval::valReadReferenceRows(f)
        dis <- msrval::valReferenceDisagreements(rows, e)
        adoption[nrow(adoption) + 1L, ] <- list(f, length(unique(paste(e$lane, e$step,
                                                                       e$class))), nrow(dis))
        if (nrow(dis)) dis_all <- rbind(dis_all, cbind(file = f, dis))
    }
    writeTsv(adoption, "g1r-adoption.tsv")
    if (!is.null(dis_all)) writeTsv(dis_all, "g1r-disagreements.tsv")
    list(adoption = adoption, dis = dis_all)
}, val_error = function(e) {
    status("stop", paste(class(e)[[1L]], conditionMessage(e)))
    NULL
})
writeTsv(checks, "g1r-checks.tsv")
if (!is.null(result)) {
    if (!all(checks$pass)) {
        status("stop", paste("checks failed:", paste(checks$check[!checks$pass], collapse = ", ")))
    } else if (!is.null(result$dis)) {
        status("stop", paste("a ticket's committed rows disagree with G1(r)'s own rows:",
                             paste(unique(result$dis$file), collapse = ", ")))
    } else {
        status("built", paste("adopted unchanged:", paste(result$adoption$file, collapse = ", ")))
    }
}
