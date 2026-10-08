# CORE-004 harness: case expansion, bindings, control selection and the frozen decision rule.
# Synthetic receipts only; no systemd, no installed library.

coreHarness <- function() {
    env <- new.env(parent = globalenv())
    root <- if (exists("workflow_test_root")) workflow_test_root else here::here()
    sys.source(file.path(root, "tools/profiling/workflow_core.R"), env)
    env
}

coreSyntheticBindings <- function() {
    input <- function(rows) list(store = "/store", head_sha256 = strrep("a", 64), rows = rows,
        reference = list(entity_id = strrep("b", 64)), parameters = list())
    build <- function(key) list(library = paste0("/lib/", key), build = strrep(substr(key, 1, 1), 64),
        exact_sum = paste0("sum-", key))
    list(root = "/root", protocol_sha256 = strrep("c", 64),
        builds = list(perf108 = build("p"), core_004 = build("h"), core_005 = build("f"), core_006 = build("s")),
        outcomes = list(list(ticket = "PERF-108", outcome = "accepted", build = "perf108")),
        tickets = list(`006` = list(control = "core_005", candidate = "core_006"),
            `010` = list(control = "core_005", candidate = "core_006"),
            `017` = list(control = "core_005", candidate = "core_006"),
            `010-control` = list(control = "perf108", candidate = "core_005")),
        screen = list(p14 = input(960000L)), full = list(p14 = input(24000000L), primary = input(6225000L)))
}

coreSyntheticReceipt <- function(cases, peak = c(control = 7e9, candidate = 6e9),
                                 elapsed = c(control = 100, candidate = 100), edit = NULL) {
    receipt <- list(cases = list())
    for (case in cases) {
        operation <- list(status = "complete", measurement_valid = TRUE, lifecycle_passed = TRUE,
            swap_bytes = 0, peak_bytes = peak[[case$role]], elapsed_seconds = elapsed[[case$role]],
            phase_seconds = list(staged_result_read_summary_write = elapsed[[case$role]] / 4),
            phase_peak_bytes = list(relational_validation = peak[[case$role]] / 2,
                working_memory_admission = peak[[case$role]] / 3,
                relational_aggregate_order_summary_write = peak[[case$role]],
                ordered_result_stage = peak[[case$role]]),
            phase_status = list(relational_validation = "complete", working_memory_admission = "complete",
                relational_aggregate_order_summary_write = "complete", ordered_result_stage = "complete"),
            accounting = list(cpu_seconds = elapsed[[case$role]]),
            execution = list(cleanup = list(worker_stopped = TRUE, monitor_stopped = TRUE),
                monitor = list(max_observation_interval = 0.05), peak_bytes = peak[[case$role]],
                elapsed_seconds = elapsed[[case$role]]))
        receipt$cases[[case$id]] <- list(passed = TRUE, operation = operation,
            verifier_completion = list(accepted = TRUE))
    }
    if (is.function(edit)) receipt <- edit(receipt)
    receipt
}

test_that("CORE stages expand to fixed cases with per-role builds and bound digests", {
    h <- coreHarness()
    b <- coreSyntheticBindings()
    life <- h$coreCases(b, "core-lifecycle-complete")
    expect_length(life, 1L)
    expect_identical(life[[1]]$selector, "core:lifecycle:complete")
    expect_identical(life[[1]]$library, "/lib/h")
    expect_identical(life[[1]]$input$rows, 960000L)
    expect_false(life[[1]]$expect_failure)
    fail <- h$coreCases(b, "core-lifecycle-fail")[[1]]
    expect_identical(fail$fault, "after_aggregate")
    expect_true(fail$expect_failure)
    expect_false(identical(fail$binding, life[[1]]$binding))
    pair <- h$coreCases(b, "core-006")
    expect_identical(vapply(pair, `[[`, character(1), "id"),
        c("control-1", "candidate-1", "candidate-2", "control-2"))
    expect_identical(vapply(pair, `[[`, character(1), "library"), c("/lib/f", "/lib/s", "/lib/s", "/lib/f"))
    expect_identical(pair[[2]]$selector, "core:006:candidate:1")
    expect_identical(pair[[1]]$input$rows, 24000000L)
    expect_identical(length(unique(vapply(pair, `[[`, character(1), "binding"))), 4L)
    expect_identical(h$coreCases(b, "core-006"), pair)
    changed <- b
    changed$protocol_sha256 <- strrep("d", 64)
    expect_false(identical(h$coreCases(changed, "core-006")[[1]]$binding, pair[[1]]$binding))
    expect_error(h$coreCases(b, "core-007"), "no bound control")
    expect_error(h$coreCases(b, "core-011-overflow"), "CORE-011")
    expect_error(h$coreCases(b, "core-015"))
    expect_identical(h$coreSelectorStage("core:006:candidate:2"), list(stage = "core-006", id = "candidate-2"))
    expect_identical(h$coreSelectorStage("core:lifecycle:fail"), list(stage = "core-lifecycle-fail", id = "member"))
    expect_error(h$coreSelectorStage("core:015:candidate:1"))
})

test_that("the control is the most recent accepted build, never a reverted one", {
    h <- coreHarness()
    outcomes <- list(list(ticket = "PERF-108", outcome = "accepted", build = "perf108"),
        list(ticket = "CORE-005", outcome = "accepted", build = "core_005"),
        list(ticket = "CORE-006", outcome = "reverted", build = "core_006"),
        list(ticket = "CORE-007", outcome = "pending", build = "core_007"))
    expect_identical(h$coreControl(outcomes), "core_005")
    expect_identical(h$coreControl(outcomes[1]), "perf108")
    expect_error(h$coreControl(outcomes[3:4]), "no accepted build")
})

test_that("coreDecision accepts only when every eligibility and comparison rule passes", {
    h <- coreHarness()
    b <- coreSyntheticBindings()
    cases <- h$coreCases(b, "core-006")
    accepted <- h$coreDecision(coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 7e9),
        elapsed = c(control = 100, candidate = 90)), cases, "006")
    expect_identical(accepted$status, "accept")
    expect_null(accepted$failing_rule)
    expect_false(accepted$promotion_authorized)
    invalid_monitor <- h$coreDecision(coreSyntheticReceipt(cases, edit = function(r) {
        r$cases[["candidate-2"]]$operation$execution$monitor$max_observation_interval <- 0.26
        r
    }), cases, "006")
    expect_identical(invalid_monitor$status, "invalid")
    expect_identical(invalid_monitor$failing_rule, "candidate-2:measurement_valid")
    # Protocol v5: a busy host may stretch the monitor gap to 250 ms; the kernel peak is unaffected.
    busy <- h$coreDecision(coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 7e9),
        elapsed = c(control = 100, candidate = 90), edit = function(r) {
            r$cases[["control-2"]]$operation$execution$monitor$max_observation_interval <- 0.131
            r$cases[["control-2"]]$operation$measurement_valid <- FALSE
            r
        }), cases, "006")
    expect_identical(busy$status, "accept")
    expect_equal(busy$time_pair_ratios$wall, c(0.9, 0.9))
    lost <- h$coreDecision(coreSyntheticReceipt(cases, edit = function(r) {
        r$cases[["candidate-1"]]$operation$execution$monitor_lost_after_result <- TRUE
        r
    }), cases, "006")
    expect_identical(lost$failing_rule, "candidate-1:measurement_valid")
    failed_control <- h$coreDecision(coreSyntheticReceipt(cases, edit = function(r) {
        r$cases[["control-2"]]$operation$status <- "memory_guard"
        r
    }), cases, "006")
    expect_identical(failed_control$status, "invalid")
    expect_identical(failed_control$failing_rule, "control-2:status_complete")
    unverified <- h$coreDecision(coreSyntheticReceipt(cases, edit = function(r) {
        r$cases[["control-1"]]$verifier_completion$accepted <- FALSE
        r
    }), cases, "006")
    expect_identical(unverified$failing_rule, "control-1:verifier_accepted")
    swap <- h$coreDecision(coreSyntheticReceipt(cases, edit = function(r) {
        r$cases[["candidate-1"]]$operation$swap_bytes <- 4096
        r
    }), cases, "006")
    expect_identical(swap$failing_rule, "candidate-1:zero_swap")
    slow <- h$coreDecision(coreSyntheticReceipt(cases, elapsed = c(control = 100, candidate = 104),
        peak = c(control = 7e9, candidate = 7e9)), cases, "006")
    expect_identical(slow$status, "reject")
    expect_identical(slow$failing_rule, "primary_metric_worse:phase:staged_result_read_summary_write")
    # Time is reported, not judged (protocol v5): slower CPU and wall time still pass.
    contended <- h$coreDecision(coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 7e9),
        elapsed = c(control = 100, candidate = 90), edit = function(r) {
            for (id in c("candidate-1", "candidate-2")) {
                r$cases[[id]]$operation$elapsed_seconds <- 140
                r$cases[[id]]$operation$accounting$cpu_seconds <- 130
            }
            r
        }), cases, "006")
    expect_identical(contended$status, "accept")
    expect_equal(contended$time_pair_ratios, list(cpu = c(1.3, 1.3), wall = c(1.4, 1.4)))
    heavier <- h$coreDecision(coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 7.1e9),
        elapsed = c(control = 100, candidate = 90)), cases, "006")
    expect_identical(heavier$failing_rule, "guard_breached:peak")
    missing <- h$coreDecision(coreSyntheticReceipt(cases, edit = function(r) {
        for (id in names(r$cases)) r$cases[[id]]$operation$phase_seconds <- NULL
        r
    }), cases, "006")
    expect_identical(missing$failing_rule, "primary_metric_unavailable:phase:staged_result_read_summary_write")
})

test_that("a memory gain with slower time is accepted and the time cost reported", {
    h <- coreHarness()
    cases <- h$coreCases(coreSyntheticBindings(), "core-010-control")
    gain <- h$coreDecision(coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 4e9),
        elapsed = c(control = 100, candidate = 120)), cases, "010-control")
    expect_identical(gain$status, "accept")
    expect_equal(gain$ratios$peak, 4 / 7)
    expect_equal(gain$time_pair_ratios$cpu, c(1.2, 1.2))
    worse <- h$coreDecision(coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 7.2e9),
        elapsed = c(control = 100, candidate = 120)), cases, "010-control")
    expect_identical(worse$status, "reject")
    expect_identical(worse$failing_rule, "primary_metric_worse:peak")
})

test_that("lifecycle decisions never relabel a failure as success", {
    h <- coreHarness()
    b <- coreSyntheticBindings()
    complete <- h$coreCases(b, "core-lifecycle-complete")
    expect_identical(h$coreLifecycleDecision(coreSyntheticReceipt(complete), complete)$status,
        "lifecycle_complete_observed")
    fail <- h$coreCases(b, "core-lifecycle-fail")
    failed <- coreSyntheticReceipt(fail, edit = function(r) {
        r$cases$member$operation$status <- "worker_failed"
        r$cases$member$operation$phase_status$ordered_result_stage <- "failed"
        r
    })
    decision <- h$coreLifecycleDecision(failed, fail)
    expect_identical(decision$status, "lifecycle_failure_observed")
    expect_true(decision$expect_failure)
    expect_identical(h$coreLifecycleDecision(coreSyntheticReceipt(fail), fail)$failing_rule, "failure_recorded")
    uncleaned <- failed
    uncleaned$cases$member$operation$execution$cleanup$monitor_stopped <- FALSE
    expect_identical(h$coreLifecycleDecision(uncleaned, fail)$failing_rule, "cleanup")
    untimed <- failed
    untimed$cases$member$operation$phase_status$ordered_result_stage <- "complete"
    expect_identical(h$coreLifecycleDecision(untimed, fail)$failing_rule, "phase_timing")
    unobserved <- coreSyntheticReceipt(complete, edit = function(r) {
        r$cases$member$operation$phase_peak_bytes$working_memory_admission <- NA_real_
        r
    })
    expect_identical(h$coreLifecycleDecision(unobserved, complete)$failing_rule, "phase_timing")
    broken <- coreSyntheticReceipt(complete, edit = function(r) {
        r$cases$member$operation$peak_bytes <- 16 * 1024^3
        r
    })
    expect_identical(h$coreLifecycleDecision(broken, complete)$failing_rule, "peak_below_early_stop")
})

test_that("every guard holds per matched pair, so one breaching pair cannot be averaged away", {
    h <- coreHarness()
    expect_error(h$coreCases(coreSyntheticBindings(), "core-007"), "no bound control")
    cases <- h$coreCases(coreSyntheticBindings(), "core-017")
    # Checks-phase medians improve, but pair 2's whole-run peak alone is 1.02.
    uneven <- coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 7e9), edit = function(r) {
        for (id in c("candidate-1", "candidate-2")) r$cases[[id]]$operation$phase_peak_bytes$relational_validation <- 2e9
        r$cases[["candidate-2"]]$operation$peak_bytes <- 7.14e9
        r
    })
    decision <- h$coreDecision(uneven, cases, "017")
    expect_identical(decision$status, "reject")
    expect_identical(decision$failing_rule, "guard_breached:peak")
    expect_equal(decision$pair_ratios$peak, c(1.00, 1.02))
})

test_that("a counted member is refused below 12 GiB or with another managed unit active", {
    h <- coreHarness()
    quiet <- function(...) utils::modifyList(list(available = 13 * 1024^3, units = character(),
        profile = "performance", cpu_busy = 0.01), list(...))
    expect_null(h$coreSetupRefusal(quiet()))
    expect_null(h$coreSetupRefusal(quiet(profile = NA_character_)))
    expect_match(h$coreSetupRefusal(quiet(available = 11.9 * 1024^3)), "below 12 GiB")
    expect_match(h$coreSetupRefusal(quiet(available = 20 * 1024^3, units = "multischolarwfx0-controller.service")),
        "1 other managed unit")
    expect_match(h$coreSetupRefusal(quiet(available = NA_real_)), "unobserved")
    # Protocol v5: the power profile and other host work are recorded, never a refusal.
    expect_null(h$coreSetupRefusal(quiet(profile = "low-power", cpu_busy = 0.21)))
})

test_that("phase events and 50 ms samples give phase seconds and sampled phase peaks", {
    h <- coreHarness()
    events <- data.frame(uptime = c(10, 11, 12, 12.5, 13, 14, 20, 21), id = c(1L, 2L, 2L, 3L, 3L, 1L, 4L, 4L),
        phase = c("outer", "check", "check", "check", "check", "outer", "late", "late"),
        event = c("start", "start", "complete", "start", "complete", "complete", "start", "failed"))
    trace <- data.frame(uptime = c(10.5, 11.5, 12.25, 12.75, 13.5, 15), current = c(1, 5, 9, 3, 2, 8) * 1e9,
        anon = c(1, 4, 6, 2, 1, 7) * 1e9)
    m <- h$corePhaseMeasures(events, trace)
    expect_equal(m$phase_seconds$check, 1.5)
    expect_equal(m$phase_seconds$outer, 4)
    expect_equal(m$phase_peak_bytes$check, 5e9)
    expect_equal(m$phase_peak_bytes$outer, 9e9)
    expect_equal(m$phase_peak_anon_bytes$check, 4e9)
    expect_equal(m$phase_peak_anon_bytes$outer, 6e9)
    expect_true(is.na(m$phase_peak_bytes$late))
    expect_identical(unlist(m$phase_status), c(outer = "complete", check = "complete", late = "failed"))
    open <- h$corePhaseMeasures(events[1:7, ], trace)
    expect_identical(open$phase_status$late, "open")
    op <- list(phase_peak_bytes = list(relational_validation = 4e9, working_memory_admission = 6e9))
    expect_equal(h$coreMetric(op, "phase_peak:checks"), 6e9)
    expect_equal(h$coreMetric(op, "phase_peak:relational_validation"), 4e9)
    expect_true(is.na(h$coreMetric(list(), "phase_peak:checks")))
})

test_that("CORE-017 expands, measures the checks phase and keeps the per-pair guards", {
    h <- coreHarness()
    b <- coreSyntheticBindings()
    cases <- h$coreCases(b, "core-017")
    expect_identical(vapply(cases, `[[`, character(1), "selector"),
        c("core:017:control:1", "core:017:candidate:1", "core:017:candidate:2", "core:017:control:2"))
    expect_true(all(vapply(cases, `[[`, logical(1), "phase_events")))
    expect_true(h$coreCases(b, "core-lifecycle-fail")[[1]]$phase_events)
    expect_identical(h$coreSelectorStage("core:017:candidate:2"), list(stage = "core-017", id = "candidate-2"))
    lower <- coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 7e9), edit = function(r) {
        for (id in c("candidate-1", "candidate-2")) {
            r$cases[[id]]$operation$phase_peak_bytes$relational_validation <- 2e9
            r$cases[[id]]$operation$phase_peak_bytes$working_memory_admission <- 1.5e9
        }
        r
    })
    accepted <- h$coreDecision(lower, cases, "017")
    expect_identical(accepted$status, "accept")
    expect_equal(accepted$ratios$`phase_peak:checks`, 2 / 3.5)
    same <- h$coreDecision(coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 7e9)), cases, "017")
    expect_identical(same$status, "accept")
    heavier <- h$coreDecision(coreSyntheticReceipt(cases, peak = c(control = 7e9, candidate = 7.2e9)), cases, "017")
    expect_identical(heavier$failing_rule, "primary_metric_worse:phase_peak:checks")
    missing <- h$coreDecision(coreSyntheticReceipt(cases, edit = function(r) {
        r$cases[["control-1"]]$operation$phase_peak_bytes <- NULL
        r
    }), cases, "017")
    expect_identical(missing$failing_rule, "primary_metric_unavailable:phase_peak:checks")
})

test_that("the worker phase hook flushes boundaries on completion and on the forced failure", {
    h <- coreHarness()
    fixture <- paste("p <- getOption(\"multischolar.dia_rollup_phase\");",
        "p(\"relational_validation\", function() Sys.sleep(0.05));",
        "p(\"relational_aggregate_order_summary_write\", function() p(\"ordered_result_stage\", function() 1))")
    for (fault in c("none", "after_aggregate")) {
        dir <- withr::local_tempdir()
        request <- file.path(dir, "request.rds")
        saveRDS(list(), request)
        expression <- sub("MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))", fixture,
            h$corePhaseHookExpression(fault), fixed = TRUE)
        status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
            c("--vanilla", "-e", shQuote(expression), shQuote(request)), stdout = FALSE, stderr = FALSE))
        expect_identical(status == 0L, fault == "none")
        events <- utils::read.delim(file.path(dir, "phase-events.tsv"), stringsAsFactors = FALSE)
        m <- h$corePhaseMeasures(events, data.frame(uptime = numeric(), current = numeric()))
        expect_identical(m$phase_status$relational_validation, "complete")
        expect_identical(m$phase_status$ordered_result_stage, if (fault == "none") "complete" else "failed")
        expect_identical(m$phase_status$relational_aggregate_order_summary_write,
            if (fault == "none") "complete" else "failed")
        expect_gte(m$phase_seconds$relational_validation, 0.04)
    }
    expect_error(h$corePhaseHookExpression("cancel"))
})

test_that("CORE-010 runs a side-by-side control stage, then a fixed budget with a placeholder", {
    h <- coreHarness()
    b <- coreSyntheticBindings()
    control <- h$coreCases(b, "core-010-control")
    expect_identical(vapply(control, `[[`, character(1), "library"), c("/lib/p", "/lib/f", "/lib/f", "/lib/p"))
    expect_identical(control[[2]]$selector, "core:010-control:candidate:1")
    expect_identical(h$coreSelectorStage("core:010-control:control:2"), list(stage = "core-010-control", id = "control-2"))
    expect_error(h$coreCases(b, "core-010"), "budget is not fixed")
    # CORE-005 qualification: whole-run peak not worse per pair, time reported.
    qualified <- h$coreDecision(coreSyntheticReceipt(control, peak = c(control = 7.1e9, candidate = 7.0e9),
        elapsed = c(control = 100, candidate = 130)), control, "010-control")
    expect_identical(qualified$status, "accept")
    expect_equal(qualified$time_pair_ratios$wall, c(1.3, 1.3))
    reference <- coreSyntheticReceipt(control, edit = function(r) {
        for (id in names(r$cases)) r$cases[[id]]$operation$sampled_anon_bytes <- 6.5 * 1024^3
        r
    })
    budget <- h$coreBudgetBytes(reference, control)
    expect_equal(budget, 5.75 * 1024^3)
    # Protocol v7: with CORE-005 reverted, B comes from the PERF-108 build's members.
    split <- reference
    for (id in c("control-1", "control-2")) split$cases[[id]]$operation$sampled_anon_bytes <- 6.72 * 1024^3
    expect_equal(h$coreBudgetBytes(split, control, role = "control"), 6 * 1024^3)
    expect_equal(h$coreBallastBytes(16 * 1024^3, budget), 16 * 1024^3 - 7.75 * 1024^3)
    expect_null(h$coreBallastBytes(7 * 1024^3, budget))
    b$tickets$`010`$budget_bytes <- budget
    stage <- h$coreCases(b, "core-010")
    expect_identical(vapply(stage, `[[`, character(1), "id"), c("control-1", "candidate-1", "candidate-2"))
    expect_true(all(vapply(stage, \(x) x$budget_bytes == budget, logical(1))))
    run <- function(control_status, candidate_status = "complete") coreSyntheticReceipt(stage, edit = function(r) {
        r$cases[["control-1"]]$operation$status <- control_status
        for (id in c("candidate-1", "candidate-2")) {
            r$cases[[id]]$operation$status <- candidate_status
            r$cases[[id]]$operation$admission <- list(batches = 9L)
        }
        r
    })
    fits <- h$coreBudgetDecision(run("host_headroom"), stage, list(receipt = reference, cases = control))
    expect_identical(fits$status, "fits_budget")
    expect_identical(unname(fits$batches), c(9L, 9L))
    expect_equal(unname(fits$time_ratios$wall), c(1, 1))
    expect_identical(h$coreBudgetDecision(run("complete"), stage)$status, "budget_not_binding")
    over <- h$coreBudgetDecision(run("host_headroom", "worker_failed"), stage)
    expect_identical(over$status, "exceeds_budget")
    expect_identical(over$failing_rule, "candidate-1:worker_failed")
    odd <- h$coreBudgetDecision(run("controller_failed"), stage)
    expect_identical(odd$status, "invalid")
    expect_identical(odd$failing_rule, "control-1:unexpected_status:controller_failed")
})
