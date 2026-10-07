# Fixture-only checks: never invoke a diagnostic member or load a scientific package.
perf101TestCode <- function() {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools/profiling/workflow_perf101.R"), environment)
    environment
}

test_that("logical comparisons ignore batch partitioning and reject scientific differences", {
    code <- perf101TestCode()
    data <- data.frame(Run = rep(c("a", "b"), length.out = 20003L),
        amount = c(NA_real_, 0, seq_len(20001L)),
        annotation = factor(rep("x", 20003L), levels = c("x", "unused")))
    make_reader <- function(value, partition) {
        parts <- split(seq_len(nrow(value)), ceiling(seq_len(nrow(value)) / partition))
        function(offset, limit) {
            wanted <- seq.int(offset + 1L, offset + limit)
            # Simulate physical partitions crossing the requested logical page.
            rows <- unlist(lapply(parts, function(part) part[part %in% wanted]), use.names = FALSE)
            list(data = value[rows, , drop = FALSE], source_order = rows,
                source_row_names = as.character(rows))
        }
    }
    left <- make_reader(data, 10000L)
    expect_length(code$comparePerf101Pages(left, make_reader(data, 2880L), nrow(data)), 3L)
    for (kind in c("quantity", "row_order", "type", "factor_levels", "missing")) {
        changed <- data
        if (kind == "quantity") changed$amount[10001L] <- -1
        if (kind == "row_order") changed <- changed[c(2L, 1L, 3:nrow(changed)), ]
        if (kind == "type") changed$amount <- as.integer(changed$amount)
        if (kind == "factor_levels") levels(changed$annotation) <- c("x", "different")
        if (kind == "missing") changed$amount[1L] <- 0
        expect_error(code$comparePerf101Pages(left, make_reader(changed, 2880L),
            nrow(data)), "logical page differs", info = kind)
    }
    expect_identical(code$comparePerf101Pages(left, rows = 0L), character())
})

test_that("native binding works without a correction and refuses ambiguity", {
    code <- perf101TestCode()
    root <- withr::local_tempdir()
    put <- function(member, status) {
        folder <- file.path(root, member)
        dir.create(folder, showWarnings = FALSE)
        saveRDS(list(run_dir = folder), file.path(root, paste0(member, "-attempt.rds")))
        saveRDS(list(member = "D0", status = status), file.path(folder, "perf101-result.rds"))
    }
    expect_error(code$bindPerf101Control(root), "exactly one")
    put("D0", "complete")
    control <- code$bindPerf101Control(root)
    expect_identical(control$directory, file.path(root, "D0"))
    expect_identical(control$receipt_sha256,
        digest::digest(file = control$receipt_path, algo = "sha256"))
    put("D0R", "complete")
    expect_error(code$bindPerf101Control(root), "exactly one")
    put("D0", "failed")
    expect_identical(code$bindPerf101Control(root)$directory, file.path(root, "D0R"))
})

test_that("request source snapshots are self-contained and detect missing or changed bytes", {
    root <- normalizePath(testthat::test_path("..", ".."))
    shell <- readLines(file.path(root, "tools/profiling/run_workflow_checkpoint_proof.sh"))
    functions <- paste(shell[seq_len(grep("^proof_script_dir=", shell) - 1L)], collapse = "\n")
    temporary <- withr::local_tempdir()
    snapshot <- file.path(temporary, "source")
    script <- file.path(temporary, "snapshot.sh")
    writeLines(c(functions, 'checkpointSnapshotPerf101Sources "$1" "$2"'), script)
    expect_identical(processx::run("bash", c(script, root, snapshot))$status, 0L)
    driver <- new.env(parent = baseenv())
    sys.source(file.path(snapshot, "tools/profiling/workflow_checkpoint_dia_driver.R"), driver)
    hashes <- driver$bindPerf101Sources(snapshot)
    expect_length(hashes, 9L)
    expect_identical(unname(hashes[["tools/profiling/workflow_perf101.R"]]),
        digest::digest(file = file.path(root, "tools/profiling/workflow_perf101.R"), algo = "sha256"))
    frozen <- new.env(parent = baseenv())
    sys.source(file.path(snapshot, "tools/profiling/workflow_perf101.R"), frozen)
    expect_true(is.function(frozen$comparePerf101Pages))
    expect_identical(processx::run("bash", c(script, root, snapshot),
        error_on_status = FALSE)$status, 2L)
    # Missing or modified frozen bytes cannot masquerade as the original manifest.
    writeLines("changed", file.path(snapshot, "R/func_prot_qc_peptide_methods.R"))
    expect_false(identical(hashes, driver$bindPerf101Sources(snapshot)))
    expect_identical(unname(hashes[["R/func_prot_qc_peptide_methods.R"]]),
        digest::digest(file = file.path(root, "R/func_prot_qc_peptide_methods.R"), algo = "sha256"))
    unlink(file.path(snapshot, "R/func_prot_qc_peptide_methods.R"))
    expect_error(driver$bindPerf101Sources(snapshot), "incomplete")
})
