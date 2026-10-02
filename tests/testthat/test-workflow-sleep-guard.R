test_that("measurement sleep inhibition ends on success and failure", {
    runner <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), runner)
    check <- function(fail) {
        root <- withr::local_tempdir()
        process <- runner$checkpointInhibitSleep(root)
        on.exit({ process$kill_tree(); process$wait(2000) }, add = TRUE)
        record <- jsonlite::read_json(file.path(root, "sleep-inhibitor.json"))
        expect_true(process$is_alive())
        expect_identical(record$mode, "block")
        if (fail) stop("injected failure after inhibition")
        process
    }
    process <- check(FALSE)
    expect_false(process$is_alive())
    expect_error(check(TRUE), "injected failure after inhibition")
    entries <- jsonlite::fromJSON(processx::run("systemd-inhibit",
        c("--list", "--json=short"))$stdout, simplifyVector = FALSE)
    expect_false(any(vapply(entries, \(entry)
        identical(entry$who, paste0("MultiScholaR-proof-", Sys.getpid())), logical(1))))
})
