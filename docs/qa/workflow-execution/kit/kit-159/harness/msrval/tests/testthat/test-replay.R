# Replay (KIT-159 7.3: valReplay; statuses as G5: identical, differs, unsupported), with the
# record builder and the replay as functions so that the comparison is tested without a store.

daRows <- function(shift = 0) {
    data.frame(feature_id = c("P1", "P2"), assay = "all", contrast = "KO_vs_WT",
               log2fc = c(1, -1) + shift, p = c(0.01, 0.2), s_own = c(0.02, 0.2),
               call_statistic = "qvalue", called_own = c(TRUE, FALSE), mean_expression = c(10, 11))
}

test_that("a replay giving the same DA table is identical", {
    res <- valReplay(manager = NULL, original_da = daRows(),
                     build = function(manager) list(steps = list()),
                     replay = function(record) daRows())
    expect_identical(res$status, "identical")
    expect_true(is.na(res$class))
})

test_that("a replay giving other values differs, with the largest difference", {
    res <- valReplay(manager = NULL, original_da = daRows(),
                     build = function(manager) list(), replay = function(record) daRows(1e-3))
    expect_identical(res$status, "differs")
    expect_equal(res$max_abs_diff, 1e-3, tolerance = 1e-9)
    expect_identical(res$rows_only_one_side, 0L)
})

test_that("a replay the package refuses is unsupported, with its class", {
    refusing <- function(record) {
        cli::cli_abort("cannot replay", class = c("workflow_replay_invalid", "rlang_error"))
    }
    res <- valReplay(manager = NULL, original_da = daRows(), build = function(manager) list(),
                     replay = refusing)
    expect_identical(res$status, "unsupported")
    expect_identical(res$class, "workflow_replay_invalid")
})

# ---------------------------------------------------------------------------------------------
# Replay of declared cells on the store (G5; audit I07, re-check E5): the package's own store
# preparation and fingerprints, the prepared store holding only the import and design states, a
# differing replay reported as such, the result in provenance.
# ---------------------------------------------------------------------------------------------

fpState <- function(sha) {
    list(components = list(m = list(sha256 = sha)), support = list(), products = list(),
         summary_sha256 = "s", rows = 10L)
}

replayRecord <- function() {
    list(steps = list(list(state = "initial", operation = "import"),
                      list(state = "design", operation = "design_builder"),
                      list(state = "norm", operation = "normalization-ruv")))
}

fakeKit <- function(da = function(root, api) "d1") {
    list(prepare = function(root, pre, api) root,
         held = function(root, api) c("initial", "design"),
         fingerprint = function(root, api) list(states = list()),
         mismatches = function(origin, replay) character(), da = da)
}

test_that("T-REP-1: the store comparison names a differing state and a differing DA result", {
    root <- driverPkgRoot()
    skip_if(is.na(root), "the package checkout is not known (run through val-tests.R)")
    kit <- valReplayKit(root)
    expect_true(all(msrval:::REPLAY_KIT_FUNCTIONS %in% attr(kit, "helpers")[["function"]]))
    o <- list(states = list(import = fpState("a"), norm = fpState("b")))
    same <- valReplayCompare(kit, o, o, "d1", "d1")
    expect_identical(same$status, "identical")
    expect_true(same$identical)
    r <- o
    r$states$norm <- fpState("c")
    d <- valReplayCompare(kit, o, r, "d1", "d2")
    expect_identical(d$status, "differs")
    expect_false(d$identical)
    expect_false(d$da_identical)
    expect_setequal(d$differences, c("norm/components/m", "da_results"))
})

test_that("T-REP-2: a refused replay is unsupported with its class; a differing one is not swallowed", {
    origin <- withr::local_tempdir()
    writeLines("x", file.path(origin, "marker"))
    refusing <- function(...) {
        cli::cli_abort("no", class = c("workflow_replay_unsupported", "rlang_error"))
    }
    res <- valReplayStore(replayRecord(), origin, character(), list(), withr::local_tempdir(),
                          fakeKit(), api = function(name) NULL, replay = refusing)
    expect_identical(res$status, "unsupported")
    expect_identical(res$class, "workflow_replay_unsupported")
    expect_identical(res$preconditions, c("initial", "design"))
    expect_identical(readLines(file.path(origin, "marker")), "x")
    ok <- function(record, root, inputs, paths, scratch) invisible(NULL)
    moved <- fakeKit(da = function(root, api) if (grepl("/store$", root)) "d2" else "d1")
    res2 <- valReplayStore(replayRecord(), origin, character(), list(), withr::local_tempdir(),
                           moved, api = function(name) NULL, replay = ok)
    expect_identical(res2$status, "differs")
    expect_identical(res2$differences, "da_results")
    expect_identical(res2$prepared_states, c("initial", "design"))
    expect_true(file.exists(file.path(origin, "marker")))
})

test_that("T-REP-3: the replay block of a cell: declared, not declared, not run", {
    decl <- list(replay = list(list(arms = list("wf"), units = list("u01"))))
    expect_true(valReplayDeclared(decl, "u01", "wf"))
    expect_false(valReplayDeclared(decl, "u02", "wf"))
    expect_false(valReplayDeclared(decl, "u01", "null"))
    expect_false(valReplayDeclared(list(), "u01", "wf"))
    expect_true(valReplayDeclared(list(replay = list(list(arms = list("wf")))), "u09", "wf"))
    b <- msrval:::valReplayBlock(list(status = "differs", class = NA_character_, identical = FALSE,
                                      da_identical = FALSE, preconditions = "initial",
                                      prepared_states = "initial",
                                      differences = "da_results"), TRUE, TRUE)
    expect_identical(b$status, "differs")
    expect_identical(unlist(b$differences), "da_results")
    expect_identical(msrval:::valReplayBlock(NULL, FALSE, TRUE)$status, "not_declared")
    expect_identical(msrval:::valReplayBlock(NULL, TRUE, FALSE)$status, "not_run")
})

test_that("T-REP-4: a declared replay cell runs on metab_combined; its result is written and checked", {
    root <- driverPkgRoot()
    skip_if(is.na(root), "the package checkout is not known (run through val-tests.R)")
    dir <- withr::local_tempdir()
    code <- replayTicket(dir, "metab_combined", root)
    out <- valRunPhase("lifecycle", replayDeclaration("metab_combined"), "d1", repairUnits(1L),
                       "rep", code, dir, file.path(dir, "w"), file.path(dir, "o"), chunks = 1L,
                       arms = "wf", manifest = repairManifest(), pkg_root = root)
    prov <- jsonlite::read_json(file.path(out, "provenance.json"))
    expect_true(valCheckProvenance(file.path(out, "provenance.json")))
    cell <- prov$cells[[1L]]
    expect_identical(cell$exit, "completed")
    rp <- cell$replay
    expect_false(identical(rp$status, "not_built"))
    expect_identical(rp$status, "identical", info = paste(unlist(rp$differences), rp$class))
    pre <- unlist(rp$preconditions)
    expect_true(length(pre) >= 2L)
    # The prepared store held only the record's import and design states before the replay.
    expect_setequal(unlist(rp$prepared_states), pre)
    # E4: the extracted helpers and the per-unit truths are in provenance.
    fns <- vapply(prov$code$helpers_extracted, function(h) h[["function"]], character(1))
    expect_true(all(c("guiEntryRunLane", "kitPrepareReplayProject") %in% fns))
    expect_true(all(vapply(prov$code$helpers_extracted, function(h) {
        grepl("^[0-9a-f]{64}$", h$sha256_of_source)
    }, logical(1))))
    pu <- prov$truth$per_unit
    expect_identical(length(pu), 1L)
    expect_identical(pu[[1L]]$unit, "u01")
    expect_identical(pu[[1L]]$truth_id, "unit")
    expect_match(pu[[1L]]$sha256, "^[0-9a-f]{64}$")
})
