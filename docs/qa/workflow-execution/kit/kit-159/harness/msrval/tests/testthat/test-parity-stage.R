# The parity stage inside the bundle jobs (KIT-159 7.9, G13): the request file, the bundle's commit
# check (KIT-176 A02, PKG_PATHS by reference) and the workflow's guard on the one admitted upload.

parityRequest <- function(dir, ...) {
    req <- list(
        schema = "msrval.parity_request.v1", platform = "linux",
        routes = list("prot_lfq", "metab_lc"),
        inputs = list(list(route = "prot_lfq", role = "search_results", path = "gen/lfq.tsv",
                           sha256 = strrep("a", 64L))),
        release_commit = strrep("1", 40L), harness_commit = strrep("2", 40L),
        harness_tree = strrep("3", 40L)
    )
    over <- list(...)
    req[names(over)] <- over
    path <- file.path(dir, "req.json")
    jsonlite::write_json(req, path, auto_unbox = TRUE)
    path
}

test_that("T-PAR-1: a parity request is read and checked", {
    dir <- withr::local_tempdir()
    req <- valParityRequest(parityRequest(dir))
    expect_identical(req$platform, "linux")
    expect_identical(unlist(req$routes), c("prot_lfq", "metab_lc"))
    two <- valParityRequest(parityRequest(dir, routes = list("prot_dia_two_species")))
    expect_identical(unlist(two$routes), "prot_dia_two_species")
    expect_error(valParityRequest(parityRequest(dir, routes = list("nope"))),
                 class = "val_spec_invalid")
    expect_error(valParityRequest(parityRequest(dir, platform = "windows",
                                                routes = list("prot_dia"))),
                 class = "val_spec_invalid")
    expect_error(valParityRequest(parityRequest(dir, release_commit = "abc")),
                 class = "val_spec_invalid")
    bad_hash <- list(list(route = "prot_lfq", role = "x", path = "p", sha256 = "zz"))
    expect_error(valParityRequest(parityRequest(dir, inputs = bad_hash)),
                 class = "val_spec_invalid")
})

test_that("T-PAR-2: the bundle's commit is the release commit, or equal over PKG_PATHS", {
    repo <- gitRepo()
    first <- gitDo(repo, "rev-parse", "HEAD")
    second <- gitCommitFiles(repo, list(NOTES.txt = "notes"), "notes")
    expect_identical(valCheckBundleCommit(list(package_commit = first), first, "R", repo),
                     "equal")
    expect_identical(valCheckBundleCommit(list(package_commit = first), second, "R", repo),
                     "same_package_paths")
    third <- gitCommitFiles(repo, list(`R/f.R` = "f <- function() 2"), "change R")
    expect_error(valCheckBundleCommit(list(package_commit = first), third, "R", repo),
                 class = "val_product_tree_mismatch")
})

test_that("T-PAR-3: distributable.yml has stage parity and one guarded upload per platform job", {
    root <- driverPkgRoot()
    skip_if(is.na(root), "VAL_HARNESS_DIR not set")
    wf <- file.path(root, ".github", "workflows", "distributable.yml")
    expect_true(valCheckParityWorkflow(wf))
    text <- readLines(wf)
    dir <- withr::local_tempdir()
    unguarded <- file.path(dir, "unguarded.yml")
    writeLines(sub("steps.parity-check.outcome == 'success'", "always()", text, fixed = TRUE),
               unguarded)
    expect_error(valCheckParityWorkflow(unguarded), class = "val_spec_invalid")
    extra <- file.path(dir, "extra.yml")
    at <- grep("name: parity-results", text, fixed = TRUE)[[1L]]
    writeLines(append(text, c("      - uses: actions/upload-artifact@v4",
                              "        with:", "          name: leak", "          path: x"),
                      after = at - 2L), extra)
    expect_error(valCheckParityWorkflow(extra), class = "val_spec_invalid")
})

test_that("T-PAR-4: the request's copied L1 values must equal L1's committed evidence (5 (e))", {
    repo <- gitRepo()
    l1 <- list(prot_lfq = list(digest = strrep("d", 64L),
                               fingerprint = list(sha256_10 = strrep("e", 64L))))
    path <- "docs/qa/workflow-execution/kit/kit-159/lifecycle/l1-digests.json"
    dir.create(file.path(repo, dirname(path)), recursive = TRUE)
    jsonlite::write_json(l1, file.path(repo, path), auto_unbox = TRUE)
    commit <- gitCommitFiles(repo, list(), "l1 evidence")
    req <- list(release_commit = commit, l1_path = path, l1 = l1)
    expect_true(valCheckParityL1(req, repo, at = commit))
    req$l1$prot_lfq$digest <- strrep("f", 64L)
    expect_error(valCheckParityL1(req, repo, at = commit), class = "val_spec_invalid")
})
