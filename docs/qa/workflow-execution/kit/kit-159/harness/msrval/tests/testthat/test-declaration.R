# Declarations in two stages and their check (KIT-159 5.3; section 10.3: T-DECL-1 to T-DECL-17,
# T-ENV-1).

EV <- "docs/qa/workflow-execution/kit/kit-999"

baseLibrary <- function() {
    data.frame(package = c("limma", "shiny", "testthat"), version = c("3.66.0", "1.14.0", "3.3.2"))
}

# Write and commit one declaration stage alone; returns the commit.
declStage <- function(repo, stage, release, code = list(prepare.R = "prepareInputs <- function() 1",
                                                       scoring.R = "score <- function() 1"),
                      changes = list(), preparation = NULL, extra_files = list(),
                      design = "sample\tgroup\nA1\tA") {
    sd <- file.path(EV, stage)
    files <- list()
    files[[file.path(sd, "declaration.md")]] <- "# Declaration"
    files[[file.path(sd, "inputs.tsv")]] <- "dataset\trole"
    files[[file.path(sd, "design.tsv")]] <- design
    for (p in names(code)) files[[file.path(sd, p)]] <- code[[p]]
    for (p in names(extra_files)) files[[p]] <- extra_files[[p]]
    for (p in names(files)) {
        dir.create(dirname(file.path(repo, p)), recursive = TRUE, showWarnings = FALSE)
        writeLines(files[[p]], file.path(repo, p))
    }
    sha <- function(p) sha256File(file.path(repo, sd, p))
    decl <- list(
        schema = "msrval.declaration.v1", ticket = "KIT-999", stage = stage, revision = "1.0",
        date = "2026-11-23", release_commit = release,
        harness_tree = gitDo(repo, "rev-parse", paste0("HEAD:", HARNESS_PATH)),
        environment = list(tools_ref = strrep("1", 40), library_sha256 = valLibrarySha(baseLibrary()),
                           extra_packages = list("lintr", "cyclocomp")),
        inputs = list(path = "inputs.tsv", sha256 = sha("inputs.tsv")),
        preparation = preparation,
        code = lapply(c(names(code), "design.tsv"), function(p) list(path = p, sha256 = sha(p))),
        changes_from_lifecycle = changes,
        datasets = list(), arms = list(),
        metrics = list(list(id = "mae", `function` = "valMedianAbsError", class = "accuracy",
                            margin = 0.05, conf = 0.95)),
        blind = list(seed = 159101L), ceilings = list(job_minutes = 240L),
        attempts = list(lifecycle = 1L, counted = 1L, correction = 1L)
    )
    jsonlite::write_json(decl, file.path(repo, sd, "declaration.json"), auto_unbox = TRUE,
                         null = "null", pretty = TRUE)
    gitDo(repo, "add", "-A")
    gitDo(repo, "commit", "-q", "-m", paste("declaration", stage))
    gitDo(repo, "rev-parse", "HEAD")
}

lifecycleProvenance <- function(repo, label = "life1") {
    p <- list(schema = "msrval.provenance.v1", stage = "lifecycle", phase = "lifecycle",
              blinded = TRUE, label = label)
    path <- file.path(EV, "runs", label, "provenance.json")
    dir.create(file.path(repo, dirname(path)), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(p, file.path(repo, path), auto_unbox = TRUE)
    gitDo(repo, "add", "-A")
    gitDo(repo, "commit", "-q", "-m", "lifecycle provenance")
    gitDo(repo, "rev-parse", "HEAD")
}

check <- function(repo, stage = "lifecycle", ...) {
    valCheckDeclaration(repo, EV, stage, at = "HEAD", ...)
}

test_that("T-DECL-1, T-DECL-2: a missing declaration; a declaration committed with another file", {
    repo <- gitRepo()
    release <- gitDo(repo, "rev-parse", "HEAD")
    expect_error(check(repo), class = "val_declaration_missing")
    declStage(repo, "lifecycle", release,
              extra_files = stats::setNames(list("notes"), file.path(EV, "lifecycle/notes.md")))
    expect_error(check(repo), class = "val_declaration_not_committed_alone")
})

test_that("T-DECL-3, T-DECL-4: a changed design file or scoring file is refused", {
    repo <- gitRepo()
    release <- gitDo(repo, "rev-parse", "HEAD")
    declStage(repo, "lifecycle", release)
    expect_true(check(repo)$ok)
    gitCommitFiles(repo, stats::setNames(list("sample\tgroup\nA1\tB"),
                                         file.path(EV, "lifecycle/design.tsv")))
    expect_error(check(repo), class = "val_declaration_changed")
    repo2 <- gitRepo()
    declStage(repo2, "lifecycle", gitDo(repo2, "rev-parse", "HEAD"))
    gitCommitFiles(repo2, stats::setNames(list("score <- function() 2"),
                                          file.path(EV, "lifecycle/scoring.R")))
    expect_error(check(repo2), class = "val_declaration_changed")
})

test_that("T-DECL-5 to T-DECL-7: evidence beside the harness is accepted; harness and product changes are refused", {
    repo <- gitRepo()
    release <- gitDo(repo, "rev-parse", "HEAD")
    declStage(repo, "lifecycle", release)
    gitCommitFiles(repo, stats::setNames(list("gate notes"),
                                         "docs/qa/workflow-execution/kit/kit-159/gates.md"))
    expect_true(check(repo)$ok)
    gitCommitFiles(repo, stats::setNames(list("harness v2"), file.path(HARNESS_PATH, "VERSION")))
    expect_error(check(repo), class = "val_harness_mismatch")
    repo2 <- gitRepo()
    declStage(repo2, "lifecycle", gitDo(repo2, "rev-parse", "HEAD"))
    gitCommitFiles(repo2, list(`R/f.R` = "f <- function() 2"))
    expect_error(check(repo2), class = "val_product_tree_mismatch")
})

test_that("T-DECL-8 and T-ENV-1: the library table against the declared hash and extras", {
    repo <- gitRepo()
    declStage(repo, "lifecycle", gitDo(repo, "rev-parse", "HEAD"))
    lib <- baseLibrary()
    expect_true(check(repo, library_table = lib)$ok)
    with_all <- rbind(lib, data.frame(package = c("lintr", "cyclocomp"), version = c("3.2", "1.1")))
    expect_true(check(repo, library_table = with_all)$ok)
    with_part <- rbind(lib, data.frame(package = "lintr", version = "3.2"))
    expect_error(check(repo, library_table = with_part), class = "val_environment_mismatch")
    other <- lib
    other$version[1] <- "3.67.0"
    expect_error(check(repo, library_table = other), class = "val_environment_mismatch")
    env <- valCheckEnvironment(with_all, valLibrarySha(lib), c("lintr", "cyclocomp"))
    expect_identical(env$extras, "whole_closure")
    expect_identical(valCheckEnvironment(lib, valLibrarySha(lib), "lintr")$extras, "none")
})

test_that("T-DECL-9 to T-DECL-11 and T-DECL-16: the two-stage order and the listed changes", {
    repo <- gitRepo()
    release <- gitDo(repo, "rev-parse", "HEAD")
    declStage(repo, "lifecycle", release)
    declStage(repo, "counted", release)
    expect_error(check(repo, "counted"), class = "val_declaration_after_run")
    repo2 <- gitRepo()
    release2 <- gitDo(repo2, "rev-parse", "HEAD")
    declStage(repo2, "lifecycle", release2)
    lifecycleProvenance(repo2)
    declStage(repo2, "counted", release2)
    expect_true(check(repo2, "counted")$ok)
    repo3 <- gitRepo()
    release3 <- gitDo(repo3, "rev-parse", "HEAD")
    declStage(repo3, "lifecycle", release3)
    lifecycleProvenance(repo3)
    declStage(repo3, "counted", release3, code = list(prepare.R = "prepareInputs <- function() 1",
                                                      scoring.R = "score <- function() 3"))
    expect_error(check(repo3, "counted"), class = "val_declaration_changed")
    repo4 <- gitRepo()
    release4 <- gitDo(repo4, "rev-parse", "HEAD")
    declStage(repo4, "lifecycle", release4)
    lifecycleProvenance(repo4)
    declStage(repo4, "counted", release4,
              code = list(prepare.R = "prepareInputs <- function() 1",
                          scoring.R = "score <- function() 3"),
              changes = list(list(path = "scoring.R", reason = "a fix decided in the lifecycle")))
    res <- check(repo4, "counted")
    expect_true(res$ok)
    expect_true("scoring.R" %in% res$declaration_diff$path)
    repo5 <- gitRepo()
    release5 <- gitDo(repo5, "rev-parse", "HEAD")
    declStage(repo5, "lifecycle", release5)
    lifecycleProvenance(repo5)
    declStage(repo5, "counted", release5, code = list(prepare.R = "prepareInputs <- function() 1",
                                                      scoring.R = "score <- function() 1",
                                                      extra.R = "x <- 1"))
    expect_error(check(repo5, "counted"), class = "val_declaration_changed")
})

test_that("T-DECL-12: a correction at another commit than the counted run's is refused", {
    repo <- gitRepo()
    release <- gitDo(repo, "rev-parse", "HEAD")
    declStage(repo, "lifecycle", release)
    lifecycleProvenance(repo)
    counted_commit <- declStage(repo, "counted", release)
    prov <- list(schema = "msrval.provenance.v1", stage = "counted", phase = "counted",
                 code = list(package_commit = counted_commit))
    dir.create(file.path(repo, EV, "runs/c1"), recursive = TRUE)
    jsonlite::write_json(prov, file.path(repo, EV, "runs/c1/provenance.json"), auto_unbox = TRUE)
    gitDo(repo, "add", "-A")
    gitDo(repo, "commit", "-q", "-m", "counted evidence")
    expect_error(valCheckCorrection(repo, EV, "c1", at = "HEAD"), class = "val_declaration_changed")
    expect_true(valCheckCorrection(repo, EV, "c1", at = counted_commit))
})

test_that("T-DECL-13: an amended run with its listed replacements is accepted and stamped", {
    repo <- gitRepo()
    release <- gitDo(repo, "rev-parse", "HEAD")
    declStage(repo, "lifecycle", release)
    gitCommitFiles(repo, stats::setNames(list("score <- function() 9"),
                                         file.path(EV, "lifecycle/scoring.R")))
    amend <- list(label = "made after the results were seen",
                  code = list(list(path = "scoring.R",
                                   sha256 = sha256File(file.path(repo, EV, "lifecycle/scoring.R")))))
    jsonlite::write_json(amend, file.path(repo, EV, "lifecycle/declaration-amendments.json"),
                         auto_unbox = TRUE)
    gitDo(repo, "add", "-A")
    gitDo(repo, "commit", "-q", "-m", "amendment")
    res <- check(repo, phase = "amended")
    expect_identical(res$amended, "made after the results were seen")
    gitCommitFiles(repo, stats::setNames(list("prepareInputs <- function() 2"),
                                         file.path(EV, "lifecycle/prepare.R")))
    expect_error(check(repo, phase = "amended"), class = "val_declaration_changed")
})

test_that("T-DECL-14: a changed preparation file is refused", {
    repo <- gitRepo()
    release <- gitDo(repo, "rev-parse", "HEAD")
    prep_commit <- gitCommitFiles(repo, stats::setNames(list("a\t1"),
                                                        file.path(EV, "rollup/R0/m1.tsv")))
    prep <- list(commit = prep_commit,
                 files = list(list(path = "rollup/R0/m1.tsv",
                                   sha256 = sha256File(file.path(repo, EV, "rollup/R0/m1.tsv")))))
    declStage(repo, "lifecycle", release, preparation = prep)
    expect_true(check(repo)$ok)
    gitCommitFiles(repo, stats::setNames(list("a\t2"), file.path(EV, "rollup/R0/m1.tsv")))
    expect_error(check(repo), class = "val_declaration_changed")
})

test_that("T-DECL-15: annotation 'none' on identifiers more than 1 % UniProt-like is refused", {
    ids_refseq <- sprintf("WP_%09d.1", 1:200)
    expect_true(valCheckAnnotationExemption(ids_refseq, list(expected = "none", reason = "RefSeq")))
    ids_mixed <- c(ids_refseq, c("P69905", "Q9Y6K9", "O15530"))
    expect_error(valCheckAnnotationExemption(ids_mixed, list(expected = "none", reason = "RefSeq")),
                 class = "val_spec_invalid")
})

test_that("T-DECL-17: the evidence check refuses two corrections of one counted run", {
    runs <- data.frame(run_id = c("1", "2", "3"), label = c("c1", "c1", "c1"),
                       stage = "counted", phase = c("counted", "correction", "correction"))
    expect_error(valCheckRuns(runs), class = "val_spec_invalid")
    expect_true(valCheckRuns(runs[1:2, ]))
})
