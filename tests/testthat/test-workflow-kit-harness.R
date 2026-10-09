# KIT-001 harness: step map, case expansion, selectors, phase measures and the lifecycle decision.
# Synthetic receipts only; no systemd, no installed library.

kitHarness <- function() {
    env <- new.env(parent = globalenv())
    root <- if (exists("workflow_test_root")) workflow_test_root else here::here()
    sys.source(file.path(root, "tools/profiling/workflow_kit.R"), env)
    env
}

kitSyntheticBindings <- function() {
    file <- function(name) list(path = paste0("/corpus/", name), bytes = 1, sha256 = strrep("a", 64))
    corpus <- function(samples) list(profile = "calibration", samples = samples,
        files = list(`report.tsv` = file("report.tsv"), `design.tsv` = file("design.tsv"), `proteins.fasta` = file("proteins.fasta")))
    build <- function(key) list(library = paste0("/lib/", key), build = strrep(substr(key, 1, 1), 64), exact_sum = paste0("sum-", key))
    list(root = "/root", protocol_sha256 = strrep("c", 64), step_map_sha256 = strrep("d", 64),
        builds = list(core_018 = build("n"), kit_006 = build("k")), control = "core_018",
        outcomes = list(list(ticket = "CORE-018", outcome = "control", build = "core_018")),
        corpora = list(screen = corpus(8L), full = corpus(300L)),
        journey = list(member = list(programme = "kit-journey-v1", deadline = 5400),
            settings = list(), ceiling_seconds = 3600),
        lifecycle = list(ceiling_seconds = 300L, fault = list(name = "guard_after_intensity", after_step = "intensity", failed_step = "sample")))
}

kitSyntheticReceipt <- function(cases, h, edit = NULL) {
    receipt <- list(cases = list())
    steps <- h$kitStepIds()
    for (case in cases) {
        operation <- list(status = "complete", swap_bytes = 0, peak_bytes = 6e9, elapsed_seconds = 200,
            max_sample_interval = 0.06, cleanup = list(slice_stopped = TRUE), resident_end_bytes = 1.2e9,
            populations_complete = TRUE, head_unchanged = NA,
            phase_status = as.list(stats::setNames(rep("complete", length(steps)), steps)),
            phase_peak_anon_bytes = as.list(stats::setNames(rep(3e9, length(steps)), steps)))
        receipt$cases[[case$id]] <- list(passed = TRUE, operation = operation, verifier_completion = list(accepted = TRUE))
    }
    if (is.function(edit)) receipt <- edit(receipt)
    receipt
}

test_that("the journey step map lists every DIA step in production order", {
    h <- kitHarness()
    ids <- h$kitStepIds()
    expect_identical(ids, c("import", "design", "qvalue", "protein_peptide", "rollup", "intensity", "sample", "replicate",
        "imputation", "iq", "accession", "protein_intensity", "duplicates", "protein_replicate", "normalisation_ruv",
        "correlation", "session_export", "da"))
    steps <- h$kitJourneySteps()
    expect_identical(steps[[1]]$state, "initial")
    expect_identical(steps[[2]]$state, "raw_data_s4")
    expect_identical(steps[[7]]$state, "sample_filtered")
    expect_null(steps[[17]]$state)
})

test_that("KIT stages expand to fixed cases with bound corpora, builds and digests", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    life <- h$kitCases(b, "kit-lifecycle-complete")
    expect_length(life, 1L)
    expect_identical(life[[1]]$selector, "kit:lifecycle:complete")
    expect_identical(life[[1]]$library, "/lib/n")
    expect_identical(life[[1]]$corpus, "screen")
    expect_false(life[[1]]$expect_failure)
    fail <- h$kitCases(b, "kit-lifecycle-fail")[[1]]
    expect_identical(fail$fault, "guard_after_intensity")
    expect_true(fail$expect_failure)
    expect_false(identical(fail$binding, life[[1]]$binding))
    changed <- b
    changed$protocol_sha256 <- strrep("e", 64)
    expect_false(identical(h$kitCases(changed, "kit-lifecycle-complete")[[1]]$binding, life[[1]]$binding))
    changed <- b
    changed$step_map_sha256 <- strrep("f", 64)
    expect_false(identical(h$kitCases(changed, "kit-lifecycle-complete")[[1]]$binding, life[[1]]$binding))
    expect_error(h$kitCases(b, "kit-journey"), "KIT-005")
    expect_error(h$kitCases(b, "kit-006"), "KIT-005")
    b$journey$pair <- list(control = "core_018", candidate = "kit_006")
    pair <- h$kitCases(b, "kit-journey")
    expect_identical(vapply(pair, `[[`, "", "library"), c("/lib/n", "/lib/k"))
    expect_identical(vapply(pair, `[[`, "", "corpus"), c("full", "full"))
    expect_identical(pair[[2]]$selector, "kit:journey:candidate:1")
    expect_identical(h$kitSelectorStage("kit:lifecycle:fail"), list(stage = "kit-lifecycle-fail", id = "member"))
    expect_identical(h$kitSelectorStage("kit:journey:control:1"), list(stage = "kit-journey", id = "control-1"))
    expect_error(h$kitSelectorStage("core:lifecycle:fail"))
    expect_identical(h$kitMemberSeconds("kit-lifecycle-complete"), 300L)
    expect_identical(h$kitMemberSeconds("kit-journey", b), 3600L)
    b$journey$ceiling_seconds <- 6000
    expect_error(h$kitMemberSeconds("kit-journey", b), "90 min")
    expect_identical(h$kitControl(b), "core_018")
    b$outcomes <- c(b$outcomes, list(list(ticket = "KIT-006", outcome = "accepted", build = "kit_006")))
    expect_identical(h$kitControl(b), "kit_006")
})

test_that("step marks and samples give per-step anonymous peaks and the resident end", {
    h <- kitHarness()
    phases <- data.frame(uptime = c(10, 14, 14, 20, 20, 25, 31),
        step = c("import", "import", "design", "design", "qvalue", "qvalue", "journey"),
        event = c("start", "complete", "start", "complete", "start", "failed", "teardown"), stringsAsFactors = FALSE)
    trace <- data.frame(uptime = c(11, 13, 15, 19, 21, 24, 26.5, 28, 30), current = c(2, 5, 4, 7, 3, 6, 4, 3, 2) * 1e9,
        anon = c(1, 4, 3, 6, 2, 5, 3, 1.5, 1.2) * 1e9)
    m <- h$kitPhaseMeasures(phases, trace)
    expect_equal(m$phase_seconds$import, 4)
    expect_equal(m$phase_peak_anon_bytes$design, 6e9)
    expect_equal(m$phase_peak_bytes$import, 5e9)
    expect_identical(unlist(m$phase_status), c(import = "complete", design = "complete", qvalue = "failed"))
    expect_equal(m$resident_end_bytes, 1.2e9)
    open <- h$kitPhaseMeasures(phases[1:5, ], trace)
    expect_identical(open$phase_status$qvalue, "open")
    expect_true(is.na(open$resident_end_bytes))
})

test_that("lifecycle decisions never relabel a failure as success", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    complete <- h$kitCases(b, "kit-lifecycle-complete")
    expect_identical(h$kitLifecycleDecision(kitSyntheticReceipt(complete, h), complete, b)$status, "lifecycle_complete_observed")
    unobserved <- kitSyntheticReceipt(complete, h, edit = function(r) { r$cases$member$operation$phase_peak_anon_bytes$da <- NA_real_; r })
    expect_identical(h$kitLifecycleDecision(unobserved, complete, b)$failing_rule, "phase_timing")
    no_end <- kitSyntheticReceipt(complete, h, edit = function(r) { r$cases$member$operation$resident_end_bytes <- NA_real_; r })
    expect_identical(h$kitLifecycleDecision(no_end, complete, b)$failing_rule, "resident_end_unobserved")
    gap <- kitSyntheticReceipt(complete, h, edit = function(r) { r$cases$member$operation$max_sample_interval <- 0.3; r })
    expect_identical(h$kitLifecycleDecision(gap, complete, b)$failing_rule, "measurement_valid")
    slow <- kitSyntheticReceipt(complete, h, edit = function(r) { r$cases$member$operation$elapsed_seconds <- 301; r })
    expect_identical(h$kitLifecycleDecision(slow, complete, b)$failing_rule, "elapsed_within_ceiling")
    partial <- kitSyntheticReceipt(complete, h, edit = function(r) { r$cases$member$operation$populations_complete <- FALSE; r })
    expect_identical(h$kitLifecycleDecision(partial, complete, b)$failing_rule, "populations_incomplete")
    fail <- h$kitCases(b, "kit-lifecycle-fail")
    failed <- kitSyntheticReceipt(fail, h, edit = function(r) {
        r$cases$member$operation$status <- "journey_failed"
        r$cases$member$operation$phase_status$sample <- "failed"
        r$cases$member$operation$guard <- list(reason = "injected_guard_stop")
        r$cases$member$operation$head_unchanged <- TRUE
        r
    })
    decision <- h$kitLifecycleDecision(failed, fail, b)
    expect_identical(decision$status, "lifecycle_failure_observed")
    expect_true(decision$expect_failure)
    expect_identical(h$kitLifecycleDecision(kitSyntheticReceipt(fail, h), fail, b)$failing_rule, "failure_recorded")
    moved <- failed
    moved$cases$member$operation$head_unchanged <- FALSE
    expect_identical(h$kitLifecycleDecision(moved, fail, b)$failing_rule, "head_unchanged")
    wrong_step <- failed
    wrong_step$cases$member$operation$phase_status$sample <- "complete"
    expect_identical(h$kitLifecycleDecision(wrong_step, fail, b)$failing_rule, "phase_timing")
})

test_that("a counted member is refused below 12 GiB or with another managed unit active", {
    h <- kitHarness()
    quiet <- function(...) utils::modifyList(list(available = 13 * 1024^3, units = character(), profile = "balanced", cpu_busy = 0.1), list(...))
    expect_null(h$kitSetupRefusal(quiet()))
    expect_match(h$kitSetupRefusal(quiet(available = 11.9 * 1024^3)), "below 12 GiB")
    expect_match(h$kitSetupRefusal(quiet(units = "multischolarwfx0-controller.service")), "1 other managed unit")
    expect_null(h$kitSetupRefusal(quiet(profile = "low-power", cpu_busy = 0.4)))
})

test_that("path masking keeps summaries comparable across members", {
    h <- kitHarness()
    masked <- h$kitMaskPaths(list(a = "/tmp/run-1/project/iq_output.txt", b = 3, c = list(d = "/x/y")))
    expect_identical(masked$a, "<path>/iq_output.txt")
    expect_identical(masked$c$d, "<path>/y")
    expect_identical(masked$b, 3)
    deep <- "/leaf"
    for (i in seq_len(5000L)) deep <- list(deep)
    expect_identical(h$kitMaskPaths(deep), deep)
    plot <- structure(list(path = "/x/y"), class = "gg")
    expect_identical(h$kitMaskPaths(list(plot = plot))$plot$path, "/x/y")
    frame <- data.frame(path = "/a/b", stringsAsFactors = FALSE)
    expect_identical(h$kitMaskPaths(frame)$path, "<path>/b")
})
