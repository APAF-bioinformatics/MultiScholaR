# Tracked generator for docs/qa/workflow-execution/core/bindings.json (CORE-004; audit CORE-004-P03).
# Run from the repository root after any harness, protocol or build change, then commit both files.
# Usage: Rscript bindings-generator.R core_004=<install run dir> [core_005=<install run dir> ...]
args <- commandArgs(TRUE)
stopifnot(length(args) >= 1L, all(grepl("^core_[0-9]{3}=.+", args)))
installs <- stats::setNames(normalizePath(sub("^[^=]+=", "", args), mustWork = TRUE), sub("=.*$", "", args))
# Recorded outcome per CORE build; a control must be `accepted` (coreControl).
build_outcomes <- list(
    core_004 = list(outcome = "harness", note = "lifecycle build; not a control (harness ticket, no mechanism under test)"),
    core_005 = list(outcome = "pending", note = "accepted or reverted by CORE-009's control pair on this build (CORE-006, CORE-007 and CORE-012 ran no pairs) against the PERF-108 full_p14 receipts (peak <= 1.00, elapsed <= 1.05)"))
stopifnot(all(names(installs) %in% names(build_outcomes)))
repo <- normalizePath(".")
stopifnot(file.exists(file.path(repo, "DESCRIPTION")))
path <- file.path(repo, "docs/qa/workflow-execution/core/bindings.json")
b <- jsonlite::read_json(path)
closeout <- readRDS(file.path(repo, "docs/qa/workflow-execution/perf-108/closeout/bindings.rds"))
probe <- function(library) {
    out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(paste0(
        '.libPaths(c("', library, '", .libPaths())); ns <- asNamespace("MultiScholaR"); ',
        'cat(get("workflowExecutionBuildIdentity", ns)(), get("diaRollupNativeInfo", ns)()$manifest$binary_sha256, ',
        'normalizePath(getNamespaceInfo(ns, "path")), sep = "\\n")'))), stdout = TRUE)
    out <- tail(out, 3L)
    stopifnot(identical(out[[3L]], normalizePath(file.path(library, "MultiScholaR"))))
    list(build = out[[1L]], exact_sum = out[[2L]])
}
control <- b$control
control_probe <- probe(control$library)
stopifnot(identical(control_probe$build, control$build), identical(control_probe$exact_sum, control$exact_sum_sha256))
b$owner <- "CORE-001; extended by CORE-004 and later CORE tickets"
b$builds <- list(perf108 = list(ticket = "PERF-108", library = control$library, build = control$build,
    exact_sum = control$exact_sum_sha256))
b$outcomes <- list(list(ticket = "PERF-108", outcome = "accepted", build = "perf108"))
for (key in names(installs)) {
    install <- installs[[key]]
    library <- file.path(install, "library")
    new <- probe(library)
    ticket <- sub("^core_", "CORE-", key)
    b$builds[[key]] <- list(ticket = ticket, library = library, build = new$build, exact_sum = new$exact_sum,
        install_member = basename(install),
        source_commit = trimws(readLines(file.path(install, "source-commit.txt"))[[1L]]),
        runtime_patch_bytes = file.size(file.path(install, "runtime.patch")),
        untracked_runtime_files = file.exists(file.path(install, "untracked-runtime.sha256")))
    b$outcomes[[length(b$outcomes) + 1L]] <- c(list(ticket = ticket, build = key), build_outcomes[[key]])
}
oracle <- closeout$p14_oracle
pagehash <- function(x) digest::digest(paste(x, collapse = "\n"), algo = "sha256", serialize = FALSE)
stopifnot(identical(pagehash(oracle$screen_pages), b$oracles$p14$screen_pages_sha256),
    identical(pagehash(oracle$full_pages), b$oracles$p14$full_pages_sha256),
    identical(digest::digest(file = oracle$receipt, algo = "sha256"), oracle$receipt_sha256))
encode <- function(value) as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = NA))
b$oracles$p14$screen_page_hashes <- as.list(oracle$screen_pages)
b$oracles$p14$full_page_hashes <- as.list(oracle$full_pages)
b$oracles$p14$full_summary_json <- encode(oracle$full_summary)
reused <- closeout$reused_rollup
stage <- readRDS(reused$path)
stopifnot(identical(digest::digest(file = reused$path, algo = "sha256"), reused$sha256))
b$control_receipts <- list(full_p14 = list(source = reused$path, sha256 = reused$sha256,
    note = "PERF-108 ordered-rollup `current` members: P14 on the jEmE2F build with every experimental option off; the route the CORE control runs. Reference only: CORE runs control pairs under the CORE protocol and harness.",
    members = lapply(stats::setNames(c("p14-current-1", "p14-current-2"), c("p14-current-1", "p14-current-2")), function(id) {
        op <- stage$cases[[id]]$operation
        stopifnot(identical(op$config$build, control$build), identical(op$status, "complete"),
            isTRUE(op$measurement_valid), isTRUE(stage$cases[[id]]$verifier_completion$accepted))
        list(build = op$config$build, peak_bytes = op$peak_bytes, elapsed_seconds = op$elapsed_seconds,
            swap_bytes = op$swap_bytes, verified = TRUE)
    })))
files <- c("tools/profiling/workflow_core.R", "tools/profiling/workflow_perf106_repair.R",
    "tools/profiling/workflow_perf106.R", "tools/profiling/workflow_perf108.R",
    "tools/profiling/run_workflow_checkpoint_proof.R", "tools/profiling/run_workflow_checkpoint_proof.sh",
    "tools/profiling/workflow_checkpoint_dia_driver.R", "tests/testthat/test-workflow-execution-adapters.R",
    "tests/testthat/helper-workflow-bulk-rollup.R", "tests/testthat/helper-workflow-checkpoint-fixtures.R")
b$harness <- list(root = "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/core/runtime",
    coordinator = "tools/profiling/workflow_perf106_repair.R",
    protocol = "docs/qa/workflow-execution/core/protocol.md",
    protocol_sha256 = digest::digest(file = file.path(repo, "docs/qa/workflow-execution/core/protocol.md"), algo = "sha256"),
    plan_sha256_at_protocol_freeze = digest::digest(file = file.path(repo, "tickets/openai/06-CORE-INFRA-PLAN.md"), algo = "sha256"),
    files = as.list(stats::setNames(vapply(file.path(repo, files), \(f) digest::digest(file = f, algo = "sha256"),
        character(1), USE.NAMES = FALSE), files)),
    tickets = structure(list(), names = character()))
b$authority$plan_sha256_note <- "CORE-001 value; superseded by harness.plan_sha256_at_protocol_freeze (CORE-004)"
jsonlite::write_json(b, path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
for (key in names(installs)) cat(key, "build", b$builds[[key]]$build, "exact_sum", b$builds[[key]]$exact_sum, "\n")
