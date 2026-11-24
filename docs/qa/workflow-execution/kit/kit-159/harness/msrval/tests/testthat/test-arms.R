# The kernel identity gate (KIT-159 7.5; T-KER-1, T-KER-2), the post-cell hook's named arguments
# (7.3, revision 1.8 item 2, C22) and prepareInputs' result tables (7.8, revision 1.8 item 3).

daTable <- function(shift = 0) {
    data.frame(feature_id = paste0("k", 1:5), contrast = "B_vs_A",
               log2fc = c(1, -0.5, 0.2, 0, 2) + shift, p = c(0.01, 0.2, 0.5, 0.9, 0.001))
}

test_that("T-KER-1: a kernel equal to the workflow arm passes on seeded identity units and a refusal unit", {
    units <- sprintf("u%02d", 1:10)
    workflow <- function(u) {
        if (u == "u07") return(list(refused = TRUE, step = "norm", class = "too_few"))
        list(refused = FALSE, da = daTable())
    }
    res <- valKernelIdentity(units, kernel = workflow, workflow = workflow, identity_seed = 159102L,
                             refusal_units = "u07")
    expect_true(res$passed)
    expect_identical(length(res$units), 4L)
    expect_identical(res$units, sort(unique(res$units)))
    expect_true("u07" %in% res$units)
    again <- valKernelIdentity(units, kernel = workflow, workflow = workflow,
                               identity_seed = 159102L, refusal_units = "u07")
    expect_identical(again$units, res$units)
    none <- valKernelIdentity(units, kernel = workflow, workflow = workflow,
                              identity_seed = 159102L)
    expect_identical(none$note, "refusal_path_not_exercised")
})

test_that("T-KER-2: a deliberately different kernel is refused", {
    units <- sprintf("u%02d", 1:10)
    workflow <- function(u) list(refused = FALSE, da = daTable())
    kernel <- function(u) list(refused = FALSE, da = daTable(shift = 1e-6))
    expect_error(valKernelIdentity(units, kernel, workflow, identity_seed = 159102L),
                 class = "val_kernel_identity_failed")
    wrong_step <- function(u) list(refused = TRUE, step = "da", class = "too_few")
    refusing <- function(u) list(refused = TRUE, step = "norm", class = "too_few")
    expect_error(valKernelIdentity(units, wrong_step, refusing, identity_seed = 159102L),
                 class = "val_kernel_identity_failed")
})

test_that("the hook receives the named arguments only through formals that name them (C22)", {
    got <- new.env()
    named_hook <- function(manager, results, summary_dir, tables, unit, out_dir) {
        got$unit <- unit
        got$out_dir <- out_dir
        got$test_mode <- getOption("multischolar.test_mode")
        got$env <- Sys.getenv("MULTISCHOLAR_TEST_MODE")
        TRUE
    }
    valCallHook(named_hook, manager = NULL, results = NULL, summary_dir = NULL, tables = list(),
                unit = "u01", arm = "shipped", design_path = "d.tsv", out_dir = "/abs/val-out")
    expect_identical(got$unit, "u01")
    expect_identical(got$out_dir, "/abs/val-out")
    expect_false(got$test_mode)
    expect_identical(got$env, "false")
    dots_hook <- function(manager, results, summary_dir, tables, ...) {
        got$dots <- names(list(...))
        TRUE
    }
    valCallHook(dots_hook, manager = NULL, results = NULL, summary_dir = NULL, tables = list(),
                unit = "u01", arm = "shipped", design_path = "d.tsv", out_dir = "/abs/val-out")
    expect_length(got$dots, 0L)
})

test_that("prepareInputs' result tables are copied, checked, and an undeclared one refused", {
    outputs <- list(list(path = "ticket/feasibility/{arm}.tsv", key = "none",
                         columns = list(list(name = "group", type = "label"),
                                        list(name = "n", type = "count")),
                         labels = list(group = c("A", "B")), max_rows = 10L, lineage = "open"))
    run_dir <- withr::local_tempdir()
    tables <- list(`ticket/feasibility/a1.tsv` = data.frame(group = c("A", "B"), n = c(3, 4)))
    paths <- valCopyResultTables(tables, outputs, run_dir, arms = "a1")
    expect_true(file.exists(file.path(run_dir, "ticket/feasibility/a1.tsv")))
    expect_identical(paths, "ticket/feasibility/a1.tsv")
    bad <- list(`ticket/other.tsv` = data.frame(group = "A", n = 1))
    expect_error(valCopyResultTables(bad, outputs, withr::local_tempdir(), arms = "a1"),
                 "ticket/other.tsv is not a declared output", class = "val_publish_forbidden",
                 fixed = TRUE)
    wrong <- list(`ticket/feasibility/a1.tsv` = data.frame(group = "Z", n = 1))
    expect_error(valCopyResultTables(wrong, outputs, withr::local_tempdir(), arms = "a1"),
                 "a label outside the declared vocabulary", class = "val_publish_forbidden",
                 fixed = TRUE)
})
