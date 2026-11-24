# The entry script itself (re-check E2): `val-run.R` lifecycle, then counted, then score, on the
# toy ticket in a temporary git repository with committed declarations, so that its JSON,
# retrieval and check path runs; through completion (a refused unit inside the counted run) and
# through failures (a missing declaration; a declaration changed after its commit), each refused
# with its class and exit 3.

EVR <- "docs/qa/workflow-execution/kit/kit-998"

# Write and commit one declaration stage of the toy ticket alone; returns the commit.
toyStage <- function(repo, stage, release, ticket_code) {
    sd <- file.path(EVR, stage)
    dir.create(file.path(repo, sd), recursive = TRUE, showWarnings = FALSE)
    file.copy(ticket_code, file.path(repo, sd, "ticket.R"), overwrite = TRUE)
    writeLines("# Toy declaration (msrval entry-script test)", file.path(repo, sd, "declaration.md"))
    writeManifest(manifestRow(repo, "corpus", "open.tsv", commit = TRUE),
                  file.path(repo, sd, "inputs.tsv"))
    utils::write.table(repairUnits(5L), file.path(repo, sd, "units-d1.tsv"), sep = "\t",
                       quote = FALSE, row.names = FALSE)
    sha <- function(p) sha256File(file.path(repo, sd, p))
    decl <- repairDeclaration(arms = c("good", "null", "refusing"))
    decl$datasets[[1L]]$units <- list(table = "units-d1.tsv")
    decl <- c(list(schema = "msrval.declaration.v1", revision = "1.0", date = "2026-11-24",
                   release_commit = release,
                   harness_tree = gitDo(repo, "rev-parse", paste0("HEAD:", HARNESS_PATH)),
                   environment = list(tools_ref = strrep("1", 40),
                                      library_sha256 = valLibrarySha(valLibraryTable()),
                                      extra_packages = list()),
                   inputs = list(path = "inputs.tsv", sha256 = sha("inputs.tsv")),
                   code = list(list(path = "ticket.R", sha256 = sha("ticket.R"))),
                   attempts = list(lifecycle = 1L, counted = 1L, correction = 1L)),
              decl)
    decl$stage <- stage
    jsonlite::write_json(decl, file.path(repo, sd, "declaration.json"), auto_unbox = TRUE,
                         null = "null", pretty = TRUE)
    gitDo(repo, "add", "--", sd)
    gitDo(repo, "commit", "-q", "-m", paste("declaration", stage))
    gitDo(repo, "rev-parse", "HEAD")
}

# Move a request's run folder into the evidence folder and commit it (what the builder does after
# fetching an artifact); val-out/ is then removed.
commitRun <- function(repo, label) {
    file.copy(file.path(repo, "val-out", "runs", label), file.path(repo, EVR, "runs"),
              recursive = TRUE)
    unlink(file.path(repo, "val-out"), recursive = TRUE)
    gitDo(repo, "add", "--", file.path(EVR, "runs"))
    gitDo(repo, "commit", "-q", "-m", paste("run", label))
}

test_that("T-RUN-E2: val-run.R runs lifecycle, counted and score, and refuses with its class", {
    script <- file.path(Sys.getenv("VAL_HARNESS_DIR"), "val-run.R")
    skip_if(!file.exists(script), "the harness folder is not known (run through val-tests.R)")
    rt <- childRunnerTemp()
    skip_if(is.na(rt), "msrval is not installed in a val-lib (run through val-tests.R)")
    repo <- gitRepo()
    inputRepo(repo)
    gitDo(repo, "add", "-A")
    gitDo(repo, "commit", "-q", "-m", "data")
    release <- gitDo(repo, "rev-parse", "HEAD")
    code <- writeRepairTicket(withr::local_tempdir())
    toyStage(repo, "lifecycle", release, code)
    dir.create(file.path(repo, EVR, "runs"), recursive = TRUE, showWarnings = FALSE)

    # Lifecycle: through completion, blind.
    life <- runEntryScript(script, c("lifecycle", EVR, "lifecycle", "life1", "--dataset", "d1",
                                     "--chunks", "1"), repo, rt)
    expect_identical(life$status, 0L, info = paste(utils::tail(life$lines, 20L), collapse = "\n"))
    expect_true(any(grepl("phase=lifecycle status=completed", statusLine(life$lines))))
    lp <- file.path(repo, "val-out", "runs", "life1", "provenance.json")
    expect_true(file.exists(lp))
    expect_true(valCheckProvenance(lp))
    prov <- jsonlite::read_json(lp)
    expect_true(prov$blinded)
    expect_match(prov$code$declaration$commit, "^[0-9a-f]{40}$")
    expect_identical(prov$code$release_commit, release)
    commitRun(repo, "life1")

    # Counted: committed alone after the lifecycle's provenance; a unit refused by one arm.
    toyStage(repo, "counted", release, code)
    cnt <- runEntryScript(script, c("counted", EVR, "counted", "c1", "--dataset", "d1",
                                    "--chunks", "2"), repo, rt)
    expect_identical(cnt$status, 0L, info = paste(utils::tail(cnt$lines, 20L), collapse = "\n"))
    expect_true(any(grepl("phase=counted status=completed", statusLine(cnt$lines))))
    cells <- utils::read.delim(file.path(repo, "val-out", "runs", "c1", "cells.tsv"),
                               stringsAsFactors = FALSE)
    expect_identical(cells$exit[cells$unit == "u04" & cells$arm == "refusing"], "refused")
    expect_true(all(cells$exit[cells$arm != "refusing"] == "completed"))
    expect_true(file.exists(file.path(repo, "val-out", "runs", "c1", "results.md")))
    commitRun(repo, "c1")

    # Score: the committed counted run, once.
    sc <- runEntryScript(script, c("score", EVR, "counted", "s1", "--dataset", "d1",
                                   "--runs", "c1"), repo, rt)
    expect_identical(sc$status, 0L, info = paste(utils::tail(sc$lines, 20L), collapse = "\n"))
    expect_true(any(grepl("phase=score status=completed", statusLine(sc$lines))))
    expect_true(length(list.files(file.path(repo, "val-out"), "^comparisons\\.tsv$",
                                  recursive = TRUE)) >= 1L)

    # A missing declaration: refused with its class, exit 3, nothing run.
    unlink(file.path(repo, "val-out"), recursive = TRUE)
    miss <- runEntryScript(script, c("counted", EVR, "nostage", "x1", "--dataset", "d1"), repo, rt)
    expect_identical(miss$status, 3L)
    expect_true(any(grepl("status=refused class=val_declaration_missing",
                          statusLine(miss$lines))))
    expect_false(dir.exists(file.path(repo, "val-out", "runs", "x1")))

    # A declaration changed after its commit: refused with its class, exit 3.
    tk <- file.path(repo, EVR, "counted", "ticket.R")
    cat("\n# changed\n", file = tk, append = TRUE)
    gitDo(repo, "commit", "-q", "-a", "-m", "change the ticket code")
    chg <- runEntryScript(script, c("counted", EVR, "counted", "c2", "--dataset", "d1"), repo, rt)
    expect_identical(chg$status, 3L)
    expect_true(any(grepl("status=refused class=val_declaration_changed",
                          statusLine(chg$lines))))
    expect_false(dir.exists(file.path(repo, "val-out", "runs", "c2")))
})
