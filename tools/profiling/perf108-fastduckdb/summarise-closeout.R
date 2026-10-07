# Summarise retained receipts; never run science or grant default activation.
directory <- file.path(here::here(), "docs/qa/workflow-execution/perf-108/closeout")
bindings <- readRDS(file.path(directory, "bindings.rds"))
path <- file.path(bindings$root, "closure-filters/stage.rds")
receipt <- readRDS(path)
members <- lapply(names(receipt$cases), function(id) {
    item <- receipt$cases[[id]]
    op <- item$operation
    if (is.null(op)) return(NULL)
    data.frame(id = id, status = op$status, peak_gib = op$peak_bytes / 1024^3,
        elapsed_seconds = op$elapsed_seconds, sampled_anon_peak_gib = op$sampled_anon_bytes / 1024^3,
        sampled_file_peak_gib = op$sampled_file_bytes / 1024^3, swap_bytes = op$swap_bytes,
        measurement_valid = isTRUE(op$measurement_valid),
        max_monitor_gap_seconds = op$execution$monitor$max_observation_interval,
        independently_verified = isTRUE(item$verifier_completion$accepted))
})
utils::write.csv(do.call(rbind, members), file.path(directory, "members.csv"), row.names = FALSE)
decision_path <- file.path(dirname(path), "decision.rds")
decision <- if (file.exists(decision_path)) readRDS(decision_path) else list(status = "pending")
evidence <- function(path) list(path = path, sha256 = digest::digest(file = path, algo = "sha256"))
checks <- lapply(c("ZQ", "PC"), function(slot) {
    unit <- jsonlite::fromJSON(readLines(file.path(directory, paste0(slot, ".log")), n = 1L))$unit
    run <- file.path(dirname(bindings$root), sub("[.]service$", "", sub("^multischolar-dia003-", "", unit)))
    terminal <- jsonlite::read_json(file.path(run, "terminal-interval.json"))
    value <- readRDS(file.path(run, "perf108-result.rds"))
    list(slot = slot, complete = isTRUE(terminal$service_succeeded) && identical(value$status, "complete"),
        peak_gib = as.numeric(terminal$raw_properties$MemoryPeak) / 1024^3,
        elapsed_seconds = terminal$elapsed_seconds,
        result = evidence(file.path(run, "perf108-result.rds")),
        terminal = evidence(file.path(run, "terminal-interval.json")), run = run)
})
jsonlite::write_json(list(owner = "PERF-108", release_status = "held_pending_final_reconciliation",
    default_activation = FALSE, decision = decision, checks = checks,
    stage = evidence(path), bindings = evidence(file.path(directory, "bindings.rds")),
    reused_rollup = bindings$reused_rollup,
    measurement_basis = "whole managed operation through publication and teardown; kernel peak includes file cache",
    package_check_limitation = "--no-install source check; generic installed tests are separately recorded in ZQ"),
    file.path(directory, "result.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
