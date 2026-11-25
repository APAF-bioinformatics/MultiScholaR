# Widening in phase score (T-SCORE-2), the results note, the hand-off folder (revision 1.8 item 1,
# C21), the ticket coverage files (revision 1.8 item 11) and the shell pre-push check (5.2).

test_that("T-SCORE-2: phase score widens a run-level interval bound by bound to the requests' within-run rows", {
    run <- valNewInterval(-0.01, 0.02, "jackknife", n_units = 4, conf = 0.975, form = "replicate")
    within <- data.frame(arm = "legacy", comparator = "baseline", metric = "mae",
                         scope = "within_run", ci_low = c(-0.03, -0.005), ci_high = c(0.01, 0.015),
                         request = c("r1", "r2"))
    w <- valWidenFromRequests(run, within, arm = "legacy", comparator = "baseline", metric = "mae")
    expect_identical(c(w$ci_low, w$ci_high), c(-0.03, 0.02))
    expect_true("widened_to_within_run" %in% w$flags)
    run8 <- valNewInterval(-0.01, 0.02, "jackknife", n_units = 8, conf = 0.975, form = "replicate")
    w8 <- valWidenFromRequests(run8, within, "legacy", "baseline", "mae", widen = FALSE)
    expect_identical(c(w8$ci_low, w8$ci_high), c(-0.01, 0.02))
    expect_identical(w8$within_run_reported, c(-0.03, 0.015))
    expect_error(valWidenFromRequests(run, within, "legacy", "baseline", "mae", widen = FALSE),
                 class = "val_spec_invalid")
})

test_that("the results note prints the measured coverage beside every verdict and the attribution", {
    cmp <- data.frame(member = c("legacy", "calibrated"), metric = "mae", arm = c("legacy", "calibrated"),
                      comparator = "baseline", outcome = c("no_worse", "inconclusive"),
                      conf = 0.975, regime = "C4", coverage = c(0.972, 0.968),
                      mc_se = c(0.0037, 0.0039), statistic = "median absolute error",
                      family = "ruv_designs", m = 2L)
    att <- data.frame(role = "counts", kind = "data", licence = "CC BY 4.0",
                      redistribution = "permitted", citation = "GSE1, doi:10.0/x")
    path <- tempfile(fileext = ".md")
    valWriteResults(path, cmp, att)
    txt <- readLines(path)
    expect_true(any(grepl("legacy.*no_worse.*measured coverage 0.972 \\(Monte Carlo SE 0.0037\\)", txt)))
    # VAL-137 (2): each component at the one-sided 95 % Clopper-Pearson upper bound of its misses
    # (56 and 64 of 2,000), not at its point estimate (0.028 + 0.032 = 0.06)
    bound <- sum(stats::qbeta(0.95, c(56, 64) + 1, 2000 - c(56, 64)))
    expect_true(any(grepl(paste0("ruv_designs: m = 2, nominal family level 0.95, family error ",
                                 "bound from measured misses ", format(round(bound, 4))),
                          txt, fixed = TRUE)))
    expect_true(any(grepl("^## Sources and attribution", txt)))
    expect_true(any(grepl("GSE1, doi:10.0/x", txt)))
})

handoffDir <- function(ev) {
    h <- file.path(ev, "ticket/handoff/lipid")
    dir.create(h, recursive = TRUE)
    writeLines("PAR1 samples", file.path(h, "samples.parquet"))
    writeLines("PAR1 matrix", file.path(h, "matrix_lipid.parquet"))
    files <- list(samples.parquet = sha256File(file.path(h, "samples.parquet")),
                  matrix_lipid.parquet = sha256File(file.path(h, "matrix_lipid.parquet")))
    jsonlite::write_json(list(schema = "multischolar.integration_input.v2", files = files),
                         file.path(h, "integration_input.json"), auto_unbox = TRUE)
    root <- file.path(ev, "ticket/handoff")
    writeLines(c("sample\tlipid", "s1\tL1"), file.path(root, "sample-match.tsv"))
    writeLines(c("species\tset", "human\tup"), file.path(root, "string-requests.tsv"))
    valWriteHandoffSha(root)
    root
}

handoffOutputs <- function() {
    list(list(path = "ticket/handoff/", handoff = TRUE, key = "folder",
              columns = list(`sample-match.tsv` = c("sample", "lipid"),
                             `string-requests.tsv` = c("species", "set")),
              lineage = c("open", "tool")))
}

test_that("the hand-off folder: admitted only counted with every data input permitted, and checked", {
    dir <- inputRepo()
    man <- licenceManifest(dir)[c(1, 3), ]
    ev <- evidenceDir(man)
    root <- handoffDir(ev)
    expect_true(valCheckEvidence(ev, man, outputs = handoffOutputs(), phase = "counted"))
    expect_error(valCheckEvidence(ev, man, outputs = handoffOutputs(), phase = "lifecycle"),
                 "admitted only in the counted stage", class = "val_publish_forbidden",
                 fixed = TRUE)
    man_closed <- licenceManifest(dir)
    ev2 <- evidenceDir(man_closed)
    handoffDir(ev2)
    expect_error(valCheckEvidence(ev2, man_closed, outputs = handoffOutputs(), phase = "counted"),
                 "a data input of the job is not permitted", class = "val_publish_forbidden",
                 fixed = TRUE)
    writeLines("changed", file.path(root, "lipid/matrix_lipid.parquet"))
    valWriteHandoffSha(root)
    expect_error(valCheckEvidence(ev, man, outputs = handoffOutputs(), phase = "counted"),
                 "matrix_lipid.parquet differs from its manifest", class = "val_publish_forbidden",
                 fixed = TRUE)
    ev3 <- evidenceDir(man)
    r3 <- handoffDir(ev3)
    writeLines("x", file.path(r3, "notes.txt"))
    valWriteHandoffSha(r3)
    expect_error(valCheckEvidence(ev3, man, outputs = handoffOutputs(), phase = "counted"),
                 "notes.txt is not a hand-off file", class = "val_publish_forbidden",
                 fixed = TRUE)
    ev4 <- evidenceDir(man)
    r4 <- handoffDir(ev4)
    unlink(file.path(r4, "lipid/samples.parquet"))
    valWriteHandoffSha(r4)
    expect_error(valCheckEvidence(ev4, man, outputs = handoffOutputs(), phase = "counted"),
                 "samples.parquet named by a manifest is absent", class = "val_publish_forbidden",
                 fixed = TRUE)
    ev5 <- evidenceDir(man)
    r5 <- handoffDir(ev5)
    writeLines(c("sample\tother", "s1\tL1"), file.path(r5, "sample-match.tsv"))
    valWriteHandoffSha(r5)
    expect_error(valCheckEvidence(ev5, man, outputs = handoffOutputs(), phase = "counted"),
                 "sample-match.tsv has other columns", class = "val_publish_forbidden",
                 fixed = TRUE)
    ev6 <- evidenceDir(man)
    r6 <- handoffDir(ev6)
    writeLines("late", file.path(r6, "string-requests.tsv"))
    expect_error(valCheckEvidence(ev6, man, outputs = handoffOutputs(), phase = "counted"),
                 "string-requests.tsv differs from handoff.sha256", class = "val_publish_forbidden",
                 fixed = TRUE)
    two <- c(handoffOutputs(), list(list(path = "ticket/extra.tsv", handoff = TRUE, key = "none",
                                         columns = list(list(name = "n", type = "count")),
                                         labels = list(), max_rows = 1L, lineage = "open")))
    expect_error(valCheckEvidence(ev3, man, outputs = two, phase = "counted"),
                 class = "val_spec_invalid")
})

test_that("ticket-coverage.tsv has its fixed columns; ticket-coverage-cells.tsv is admitted", {
    dir <- inputRepo()
    man <- licenceManifest(dir)
    ev <- evidenceDir(man)
    tc <- emptyTable("ticket-coverage.tsv")
    utils::write.table(tc, file.path(ev, "ticket-coverage.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    cells <- emptyTable("ticket-coverage-cells.tsv")
    utils::write.table(cells, file.path(ev, "ticket-coverage-cells.tsv"), sep = "\t",
                       quote = FALSE, row.names = FALSE)
    expect_true(valCheckEvidence(ev, man))
    tc$extra <- 1
    utils::write.table(tc, file.path(ev, "ticket-coverage.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    expect_error(valCheckEvidence(ev, man), "ticket-coverage.tsv has other columns",
                 class = "val_publish_forbidden", fixed = TRUE)
})

test_that("check-evidence.sh: the pre-push allow-list and hash check without R", {
    harness <- Sys.getenv("VAL_HARNESS_DIR")
    expect_true(nzchar(harness))
    script <- file.path(harness, "check-evidence.sh")
    dir <- inputRepo()
    man <- licenceManifest(dir)[c(1, 3), ]
    ev <- evidenceDir(man)
    valWriteOutputsSha(ev)
    run <- function(d) suppressWarnings(system2("bash", c(script, d), stdout = TRUE, stderr = TRUE))
    ok <- run(ev)
    expect_null(attr(ok, "status"))
    writeLines("x", file.path(ev, "notes.txt"))
    valWriteOutputsSha(ev)
    expect_false(is.null(attr(run(ev), "status")))
    unlink(file.path(ev, "notes.txt"))
    valWriteOutputsSha(ev)
    write("tampered", file.path(ev, "cells.tsv"), append = TRUE)
    expect_false(is.null(attr(run(ev), "status")))
    ev2 <- evidenceDir(man)
    cells <- emptyTable("cells.tsv")
    cells$extra <- 1
    utils::write.table(cells, file.path(ev2, "cells.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    valWriteOutputsSha(ev2)
    expect_false(is.null(attr(run(ev2), "status")))
})
