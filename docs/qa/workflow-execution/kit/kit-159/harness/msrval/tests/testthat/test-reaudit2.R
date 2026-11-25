# Re-audit 2 of KIT-159 (register AUDIT-KIT-159-impl-recheck-20261124.md, section "Re-audit 2"):
# R2-01 (a missing DA result on both sides is not agreement) and R2-02 (the family level applies
# to every metric a member names).

reKit <- function(da = function(root, api) "d1") {
    list(prepare = function(root, pre, api) root,
         held = function(root, api) c("initial", "design"),
         fingerprint = function(root, api) list(states = list()),
         mismatches = function(origin, replay) character(), da = da)
}

reRecord <- function() {
    list(steps = list(list(operation = "import", state = "initial"),
                      list(operation = "design_builder", state = "design"),
                      list(operation = "normalisation", state = "norm")))
}

test_that("R2-01: no DA result on either side is unsupported, never identical", {
    o <- list(states = list())
    none <- valReplayCompare(reKit(), o, o, NA_character_, NA_character_)
    expect_identical(none$status, "unsupported")
    expect_identical(none$class, "val_replay_no_da_result")
    expect_false(none$identical)
    lost <- valReplayCompare(reKit(), o, o, "d1", NA_character_)
    expect_identical(lost$status, "differs")
    expect_identical(lost$differences, "da_results")
    same <- valReplayCompare(reKit(), o, o, "d1", "d1")
    expect_identical(same$status, "identical")
    expect_identical(same$da_fingerprint_origin, "d1")
    expect_identical(same$da_fingerprint_replay, "d1")
})

test_that("R2-01: through the store replay and into the cell's replay block", {
    origin <- withr::local_tempdir()
    writeLines("x", file.path(origin, "marker"))
    ok <- function(record, root, inputs, paths, scratch) invisible(NULL)
    res <- valReplayStore(reRecord(), origin, character(), list(), withr::local_tempdir(),
                          reKit(da = function(root, api) NA_character_),
                          api = function(name) NULL, replay = ok)
    expect_identical(res$status, "unsupported")
    expect_identical(res$class, "val_replay_no_da_result")
    good <- valReplayStore(reRecord(), origin, character(), list(), withr::local_tempdir(),
                           reKit(), api = function(name) NULL, replay = ok)
    expect_identical(good$status, "identical")
    b <- msrval:::valReplayBlock(good, TRUE, TRUE)
    expect_identical(b$da_fingerprint_origin, "d1")
    expect_identical(b$da_fingerprint_replay, "d1")
})

test_that("R2-02: a metric a family member names is at 1 - 0.05/m even without primary_family", {
    d <- repairDeclaration(family = TRUE, conf = 0.975)
    d$metrics[[2L]]$primary_family <- NULL
    expect_true(valCheckDeclarationSchema(d, repairUnits()))
    d$metrics[[2L]]$conf <- 0.95
    expect_error(valCheckDeclarationSchema(d, repairUnits()), class = "val_spec_invalid",
                 regexp = "mae2")
})
