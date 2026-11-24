# The licence guard (KIT-159 5.2, VAL-47; section 10.3: T-LIC-1 to T-LIC-11, T-LOG-1, T-LOG-2).

test_that("T-LIC-1 to T-LIC-4: unlisted file, extra column, per-feature class, too many rows", {
    dir <- inputRepo()
    man <- licenceManifest(dir)
    ev <- evidenceDir(man)
    expect_true(valCheckEvidence(ev, man))
    writeLines("x", file.path(ev, "notes.txt"))
    expect_error(valCheckEvidence(ev, man), "notes.txt is not on the allow-list",
                 class = "val_publish_forbidden", fixed = TRUE)
    ev2 <- evidenceDir(man)
    cells <- emptyTable("cells.tsv")
    cells$extra <- 1
    writeTsv(cells, file.path(ev2, "cells.tsv"))
    expect_error(valCheckEvidence(ev2, man), "cells.tsv has other columns",
                 class = "val_publish_forbidden", fixed = TRUE)
    ev3 <- evidenceDir(man)
    met <- emptyTable("metrics.tsv", 2L)
    met$class <- c("changed", "P69905")
    writeTsv(met, file.path(ev3, "metrics.tsv"))
    expect_error(valCheckEvidence(ev3, man, classes = c("changed", "null")),
                 "metrics.tsv has an undeclared class", class = "val_publish_forbidden",
                 fixed = TRUE)
    met$class <- c("changed", "null")
    writeTsv(met, file.path(ev3, "metrics.tsv"))
    expect_true(valCheckEvidence(ev3, man, classes = c("changed", "null")))
    expect_error(valCheckEvidence(ev3, man, classes = c("changed", "null"),
                                  row_bounds = list(metrics.tsv = 1L)),
                 "metrics.tsv has more rows than its bound", class = "val_publish_forbidden",
                 fixed = TRUE)
})

test_that("T-LIC-5 to T-LIC-8: copies, archives and the size of scored.tsv.gz", {
    dir <- inputRepo()
    man <- licenceManifest(dir)
    ev <- evidenceDir(man)
    dir.create(file.path(ev, "inputs"))
    file.copy(file.path(dir, "data/closed.tsv"), file.path(ev, "inputs/closed.tsv"))
    expect_error(valCheckEvidence(ev, man), "closed.tsv copies a not-permitted input",
                 class = "val_publish_forbidden", fixed = TRUE)
    ev2 <- evidenceDir(man)
    dir.create(file.path(ev2, "inputs"))
    file.copy(file.path(dir, "data/open.tsv"), file.path(ev2, "inputs/open.tsv"))
    expect_true(valCheckEvidence(ev2, man))
    ev3 <- evidenceDir(man)
    writeLines("x", file.path(ev3, "bundle.zip"))
    expect_error(valCheckEvidence(ev3, man), "bundle.zip is an archive or binary",
                 class = "val_publish_forbidden", fixed = TRUE)
    ev4 <- evidenceDir(man[man$role != "closed", ])
    big <- file.path(ev4, "scored.tsv.gz")
    con <- gzfile(big, "w")
    withr::with_seed(1, writeLines(c(paste(SCORED_COLUMNS, collapse = "\t"),
                                     replicate(90000, paste(sample(c(0:9, letters, LETTERS), 120,
                                                                   TRUE), collapse = ""))), con))
    close(con)
    expect_gt(file.size(big), 5 * 1024^2)
    expect_error(valCheckEvidence(ev4, man[man$role != "closed", ], scored_lineage = "open"),
                 "scored.tsv.gz is over 5 MB", class = "val_publish_forbidden",
                 fixed = TRUE)
})

test_that("T-LIC-9: declared outputs, keyed lineage, and path placeholders", {
    dir <- inputRepo()
    man <- licenceManifest(dir)
    outputs <- list(
        list(path = "ticket/pca/{set}-{arm}-scores.tsv", placeholders = list(set = "s1"),
             key = "none", columns = list(list(name = "arm", type = "label"),
                                          list(name = "pc1", type = "statistic")),
             labels = list(arm = c("a1", "a2")), max_rows = 10L, lineage = "closed"),
        list(path = "ticket/keyed.tsv", key = "protein", columns = list(
            list(name = "protein", type = "key"), list(name = "value", type = "statistic")),
            max_rows = 10L, lineage = c("open", "tool"))
    )
    arms <- c("a1", "a2")
    expect_setequal(valExpandOutputPaths(outputs[[1]], arms = arms),
                    c("ticket/pca/s1-a1-scores.tsv", "ticket/pca/s1-a2-scores.tsv"))
    ev <- evidenceDir(man)
    writeTsv(data.frame(arm = "a1", pc1 = 0.5), file.path(ev, "ticket/pca/s1-a1-scores.tsv"))
    writeTsv(data.frame(protein = "P1", value = 1), file.path(ev, "ticket/keyed.tsv"))
    expect_true(valCheckEvidence(ev, man, outputs = outputs, arms = arms))
    writeTsv(data.frame(arm = "a1", pc1 = 0.5, extra = 1),
             file.path(ev, "ticket/pca/s1-a2-scores.tsv"))
    expect_error(valCheckEvidence(ev, man, outputs = outputs, arms = arms),
                 "columns differ from the declaration", class = "val_publish_forbidden",
                 fixed = TRUE)
    ev2 <- evidenceDir(man)
    writeTsv(data.frame(arm = "a1", pc1 = 0.5), file.path(ev2, "ticket/pca/s2-a1-scores.tsv"))
    expect_error(valCheckEvidence(ev2, man, outputs = outputs, arms = arms),
                 "s2-a1-scores.tsv is not on the allow-list", class = "val_publish_forbidden",
                 fixed = TRUE)
    closed_keyed <- outputs
    closed_keyed[[2]]$lineage <- "closed"
    ev3 <- evidenceDir(man)
    writeTsv(data.frame(protein = "P1", value = 1), file.path(ev3, "ticket/keyed.tsv"))
    expect_error(valCheckEvidence(ev3, man, outputs = closed_keyed, arms = arms),
                 "keyed.tsv is a keyed table", class = "val_publish_forbidden",
                 fixed = TRUE)
    undeclared <- outputs[1]
    undeclared[[1]]$path <- "ticket/pca/{variant}-{arm}.tsv"
    expect_error(valExpandOutputPaths(undeclared[[1]], arms = arms), class = "val_spec_invalid")
    unused <- outputs[1]
    unused[[1]]$placeholders <- list(set = "s1", other = "x")
    expect_error(valExpandOutputPaths(unused[[1]], arms = arms), class = "val_spec_invalid")
    key_value <- outputs[1]
    key_value[[1]]$placeholders <- list(set = "P69905")
    ev4 <- evidenceDir(man)
    writeTsv(data.frame(arm = "a1", pc1 = 0.5), file.path(ev4, "ticket/pca/P69905-a1-scores.tsv"))
    expect_error(valCheckEvidence(ev4, man, outputs = key_value, arms = arms,
                                  universe = c("P69905", "P68871")),
                 "names a feature key", class = "val_publish_forbidden",
                 fixed = TRUE)
})

test_that("T-LIC-10: the lifecycle keeps scored.tsv.gz and primary call counts private", {
    dir <- inputRepo()
    man <- licenceManifest(dir)[c(1, 3), ]
    ev <- evidenceDir(man)
    con <- gzfile(file.path(ev, "scored.tsv.gz"), "w")
    writeLines(paste(SCORED_COLUMNS, collapse = "\t"), con)
    close(con)
    expect_true(valCheckEvidence(ev, man, scored_lineage = "open", phase = "counted"))
    expect_error(valCheckEvidence(ev, man, scored_lineage = "open", phase = "lifecycle"),
                 "scored.tsv.gz never leaves a lifecycle job", class = "val_publish_forbidden",
                 fixed = TRUE)
    ev2 <- evidenceDir(man)
    met <- emptyTable("metrics.tsv")
    met$unit <- "u01"
    met$metric <- "keys_called"
    writeTsv(met, file.path(ev2, "metrics.tsv"))
    expect_error(valCheckEvidence(ev2, man, phase = "lifecycle", primary_units = "u01"),
                 "holds a call count of a primary unit", class = "val_publish_forbidden",
                 fixed = TRUE)
})

test_that("T-LIC-11: result tables of a not-permitted input, label keys, attribution", {
    dir <- inputRepo()
    man <- licenceManifest(dir)
    out <- list(list(path = "ticket/summary.tsv", key = "none",
                     columns = list(list(name = "group", type = "label"),
                                    list(name = "n", type = "count")),
                     labels = list(group = c("A", "B", "P69905")), max_rows = 10L,
                     lineage = "closed"))
    ev <- evidenceDir(man)
    writeTsv(data.frame(group = c("A", "B"), n = c(3, 4)), file.path(ev, "ticket/summary.tsv"))
    expect_true(valCheckEvidence(ev, man, outputs = out, universe = c("P69905", "P68871")))
    writeTsv(data.frame(group = c("A", "P69905"), n = c(3, 4)),
             file.path(ev, "ticket/summary.tsv"))
    expect_error(valCheckEvidence(ev, man, outputs = out, universe = c("P69905", "P68871")),
                 "a label is a feature key", class = "val_publish_forbidden",
                 fixed = TRUE)
    ev2 <- evidenceDir(man)
    writeTsv(man[1, c("role", "kind", "licence", "redistribution", "citation")],
             file.path(ev2, "attribution.tsv"))
    expect_error(valCheckEvidence(ev2, man), "attribution.tsv does not list every input",
                 class = "val_publish_forbidden", fixed = TRUE)
    expect_true(valPublishable(man, "open"))
    expect_true(valPublishable(man, c("open", "tool")))
    expect_false(valPublishable(man, c("open", "closed")))
})

test_that("T-LOG-1 and T-LOG-2: a child's text, system2 output and fd 2 reach only the private log", {
    private <- withr::local_tempdir()
    script <- tempfile(fileext = ".R")
    writeLines(c(
        "writeLines('SENTINEL-FD2-7731', file('/dev/stderr'))",
        "cat('SENTINEL-OUT-7731\\n')",
        "message('SENTINEL-MSG-7731')",
        "system2('echo', 'SENTINEL-ECHO-7731')",
        "msrval::valStatus(phase = 'lifecycle', status = 'completed', count = 3L)"
    ), script)
    lines <- valRunAsChild(script, character(), private)
    expect_false(any(grepl("SENTINEL", lines)))
    expect_true(any(grepl("phase=lifecycle", lines)))
    log_text <- unlist(lapply(list.files(private, recursive = TRUE, full.names = TRUE), readLines))
    for (s in c("OUT", "MSG", "ECHO", "FD2")) {
        expect_true(any(grepl(paste0("SENTINEL-", s, "-7731"), log_text)))
    }
    expect_true(valPrivateLogs(licenceManifest(inputRepo())))
    expect_false(valPrivateLogs(licenceManifest(inputRepo())[c(1, 3), ]))
})
