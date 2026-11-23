# Units, chunks and processes (KIT-159 7.6; section 10.3: T-CHUNK-1 to T-CHUNK-7, T-RUN-1,
# T-RNG-1). Cells run in child `Rscript` processes; the toy runner's behaviour is set per cell.

toyRunner <- function(counter_dir) {
    path <- tempfile(fileext = ".R")
    writeLines(c(
        "runCell <- function(cell) {",
        sprintf("    cat('x', file = file.path('%s', paste0(cell$run, '.count')), append = TRUE)",
                counter_dir),
        "    switch(cell$behaviour,",
        "        ok = list(value = 1),",
        "        error = stop('toy failure'),",
        "        refuse = msrval::valArmRefused(step = 'norm', class = 'toy_refusal'),",
        "        busy = { repeat { x <- sum(stats::runif(1e4)) } },",
        "        kill = tools::pskill(Sys.getpid(), tools::SIGKILL),",
        "        random = list(value = stats::runif(1)),",
        "        sleep = { Sys.sleep(2); list(value = 1) })",
        "}"
    ), path)
    path
}

toyCells <- function(behaviour, run = NULL) {
    n <- length(behaviour)
    data.frame(cell = sprintf("c%02d", seq_len(n)), unit = sprintf("u%02d", seq_len(n)),
               run = run %||% sprintf("r%02d", seq_len(n)), arm = "a1", behaviour = behaviour,
               seed = 159000L + seq_len(n), stringsAsFactors = FALSE)
}

runToy <- function(cells, k = 1L, cell_minutes = 1, counter_dir = withr::local_tempdir()) {
    valRunCells(cells, runner = toyRunner(counter_dir), k = k, cell_minutes = cell_minutes,
                work_dir = withr::local_tempdir(.local_envir = parent.frame()))
}

test_that("T-CHUNK-1 and T-CHUNK-2: an erroring cell and a refusing cell keep the others", {
    res <- runToy(toyCells(c("ok", "error", "refuse", "ok")))
    expect_identical(res$exit, c("completed", "did_not_complete", "refused", "completed"))
    expect_identical(res$stop_reason[2], "error")
    expect_identical(res$step[3], "norm")
    expect_identical(res$class[3], "toy_refusal")
})

test_that("T-CHUNK-3 and T-CHUNK-5: a cell over its own limit times out and later cells complete", {
    res <- runToy(toyCells(c("ok", "busy", "ok")), cell_minutes = 0.05)
    expect_identical(res$exit, c("completed", "did_not_complete", "completed"))
    expect_identical(res$stop_reason[2], "timeout")
})

test_that("T-CHUNK-4 and T-CHUNK-6: a killed chunk; its unstarted cells complete in a fresh chunk", {
    res <- runToy(toyCells(c("ok", "kill", "ok", "ok")))
    expect_identical(res$exit, c("completed", "did_not_complete", "completed", "completed"))
    expect_identical(res$stop_reason[2], "memory")
    expect_true(all(res$fresh_chunk[3:4]))
    expect_false(res$fresh_chunk[1])
})

test_that("T-CHUNK-7: four chunks run at the same time", {
    res <- runToy(toyCells(rep("sleep", 4)), k = 4L)
    chunks <- attr(res, "chunks")
    expect_identical(nrow(chunks), 4L)
    expect_lt(max(chunks$started), min(chunks$ended))
    expect_true(all(chunks$k == 4L))
})

test_that("T-RUN-1: two unit rows naming one run execute it once and are scored twice", {
    counter <- withr::local_tempdir()
    cells <- toyCells(c("ok", "ok"), run = c("r01", "r01"))
    res <- runToy(cells, counter_dir = counter)
    expect_identical(nrow(res), 2L)
    expect_identical(res$unit, c("u01", "u02"))
    expect_identical(nchar(readLines(file.path(counter, "r01.count"), warn = FALSE)), 1L)
})

test_that("T-RNG-1: a cell that draws random numbers is flagged", {
    res <- runToy(toyCells(c("ok", "random")))
    expect_identical(res$random_state_changed, c(FALSE, TRUE))
    expect_true(all(is.finite(res$vmhwm_mib)))
})

test_that("valUnitsDropReplicate deletes replicate r from every group; a group without it stays whole", {
    design <- data.frame(sample = c("A1", "A2", "A3", "B1", "B2"), group = c("A", "A", "A", "B", "B"),
                         replicate = c(1L, 2L, 3L, 1L, 2L))
    d3 <- valUnitsDropReplicate(design, 3L)
    expect_identical(d3$design$sample, c("A1", "A2", "B1", "B2"))
    expect_true(d3$unequal_deletion)
    d1 <- valUnitsDropReplicate(design, 1L)
    expect_identical(d1$design$sample, c("A2", "A3", "B2"))
    expect_false(d1$unequal_deletion)
})
