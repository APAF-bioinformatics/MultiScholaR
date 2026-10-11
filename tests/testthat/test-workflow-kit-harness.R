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
        lifecycle = list(ceiling_seconds = 300L, fault = list(name = "guard_after_intensity", after_step = "intensity", failed_step = "sample")),
        operation_lifecycle = list(ceiling_seconds = 300L, fault = list(name = "after_aggregate", operation = "precursor_rollup",
            population = "screen_p14", failed_phase = "ordered_result_stage", error_pattern = "CORE-004 forced failure")),
        populations = list(screen_p14 = list(store = "/store/screen", head_sha256 = strrep("1", 64), entity_id = strrep("2", 64), rows = 960000L),
            full_p14 = list(store = "/store/full", head_sha256 = strrep("3", 64), entity_id = strrep("4", 64), rows = 24000000L),
            full_primary = list(store = "/store/primary", head_sha256 = strrep("5", 64), entity_id = strrep("6", 64), rows = 6225000L,
                parameters = list(peptides_per_sample_cutoff = 3000L, inclusion_list = list("run_0076", "run_0077")))),
        tickets = list(`006` = list(control = "core_018", candidate = "kit_006",
            pairs = list(list(operation = "precursor_rollup", population = "full_p14"),
                list(operation = "minimum_peptides_per_sample", population = "full_primary")))))
}

kitSyntheticOperationReceipt <- function(cases, h, edit = NULL) {
    receipt <- list(cases = list())
    phases <- h$kitRequiredOperationPhases
    for (case in cases) {
        operation <- list(status = "complete", swap_bytes = 0, peak_bytes = 5e9, elapsed_seconds = 120,
            max_sample_interval = 0.06, cleanup = list(slice_stopped = TRUE), head_unchanged = NA,
            phase_status = as.list(stats::setNames(rep("complete", length(phases)), phases)))
        receipt$cases[[case$id]] <- list(passed = TRUE, operation = operation, verifier_completion = list(accepted = TRUE))
    }
    if (is.function(edit)) receipt <- edit(receipt)
    receipt
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
    b$journey$pair <- list(control = "core_018", candidate = "kit_006")
    pair <- h$kitCases(b, "kit-journey")
    expect_identical(vapply(pair, `[[`, "", "library"), c("/lib/n", "/lib/k"))
    expect_identical(vapply(pair, `[[`, "", "corpus"), c("full", "full"))
    expect_identical(pair[[2]]$selector, "kit:journey:candidate:1")
    expect_true(all(vapply(pair, `[[`, logical(1), "record")))
    # KIT-005 baseline stage: two control repetitions and the production-default member (no journey record).
    b$journey$pair <- list(control = "core_018", candidate = NULL, reps = 2L, production_default = TRUE)
    baseline <- h$kitCases(b, "kit-journey")
    expect_identical(vapply(baseline, `[[`, "", "id"), c("control-1", "control-2", "default-1"))
    expect_identical(vapply(baseline, `[[`, "", "selector"), c("kit:journey:control:1", "kit:journey:control:2", "kit:journey:default:1"))
    expect_identical(vapply(baseline, `[[`, logical(1), "record"), c(TRUE, TRUE, FALSE))
    expect_false(identical(baseline[[1]]$binding, baseline[[2]]$binding))
    expect_false(identical(baseline[[1]]$binding, baseline[[3]]$binding))
    # Operation members (KIT-005): the lifecycle pair on screen_p14 and a ticket's control/candidate pair.
    op <- h$kitCases(b, "kit-lifecycle-operation-complete")[[1]]
    expect_identical(op$kind, "operation")
    expect_identical(op$operation, "precursor_rollup")
    expect_identical(op$input$entity_id, strrep("2", 64))
    expect_identical(op$library, "/lib/n")
    expect_false(op$expect_failure)
    lifecycle_build <- b
    lifecycle_build$operation_lifecycle$build <- "kit_006"
    expect_identical(h$kitCases(lifecycle_build, "kit-lifecycle-operation-complete")[[1]]$library, "/lib/k")
    op_fail <- h$kitCases(b, "kit-lifecycle-operation-fail")[[1]]
    expect_identical(op_fail$fault, "after_aggregate")
    expect_true(op_fail$expect_failure)
    expect_identical(op_fail$selector, "kit:lifecycle:operation:fail")
    ticket <- h$kitCases(b, "kit-006")
    expect_identical(vapply(ticket, `[[`, "", "id"), c("control-1", "candidate-1", "control-2", "candidate-2"))
    expect_identical(vapply(ticket, `[[`, "", "library"), c("/lib/n", "/lib/k", "/lib/n", "/lib/k"))
    expect_identical(ticket[[2]]$input$store, "/store/full")
    expect_identical(ticket[[2]]$selector, "kit:006:candidate:1")
    expect_identical(ticket[[3]]$operation, "minimum_peptides_per_sample")
    expect_identical(ticket[[4]]$input$store, "/store/primary")
    expect_identical(ticket[[4]]$parameters$peptides_per_sample_cutoff, 3000L)
    # JSON arrays arrive as lists of scalars; the pinned action validator needs atomic vectors (KIT-006 attempt 1).
    expect_identical(ticket[[4]]$parameters$inclusion_list, c("run_0076", "run_0077"))
    single <- b
    single$tickets[["007"]] <- list(operation = "precursor_rollup", population = "screen_p14", control = "core_018", candidate = "kit_006")
    expect_identical(vapply(h$kitCases(single, "kit-007"), `[[`, "", "id"), c("control-1", "candidate-1"))
    expect_error(h$kitCases(b, "kit-007"), "KIT-007")
    expect_identical(h$kitSelectorStage("kit:lifecycle:fail"), list(stage = "kit-lifecycle-fail", id = "member"))
    expect_identical(h$kitSelectorStage("kit:lifecycle:operation:fail"), list(stage = "kit-lifecycle-operation-fail", id = "member"))
    expect_identical(h$kitSelectorStage("kit:journey:control:1"), list(stage = "kit-journey", id = "control-1"))
    expect_identical(h$kitSelectorStage("kit:journey:default:1"), list(stage = "kit-journey", id = "default-1"))
    expect_identical(h$kitSelectorStage("kit:006:candidate:2"), list(stage = "kit-006", id = "candidate-2"))
    expect_error(h$kitSelectorStage("core:lifecycle:fail"))
    expect_error(h$kitSelectorStage("kit:journey:default:2"))
    expect_identical(h$kitMemberSeconds("kit-lifecycle-complete"), 300L)
    expect_identical(h$kitMemberSeconds("kit-lifecycle-operation-fail"), 300L)
    expect_identical(h$kitMemberSeconds("kit-006"), 600L)
    expect_identical(h$kitMemberSeconds("kit-journey", b), 3600L)
    b$journey$ceiling_seconds <- 15000
    expect_error(h$kitMemberSeconds("kit-journey", b), "4 h")
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
    # Journey members tolerate isolated gaps up to 500 ms (protocol v1.3); operation members keep 250 ms.
    gap <- kitSyntheticReceipt(complete, h, edit = function(r) { r$cases$member$operation$max_sample_interval <- 0.3; r })
    expect_identical(h$kitLifecycleDecision(gap, complete, b)$status, "lifecycle_complete_observed")
    wide <- kitSyntheticReceipt(complete, h, edit = function(r) { r$cases$member$operation$max_sample_interval <- 0.6; r })
    expect_identical(h$kitLifecycleDecision(wide, complete, b)$failing_rule, "measurement_valid")
    expect_false(h$kitMeasurementValid(gap$cases$member$operation, "operation"))
    expect_true(h$kitMeasurementValid(gap$cases$member$operation, "journey"))
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

test_that("operation lifecycle decisions need every P14 phase or the mid-step failure with the head unchanged", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    complete <- h$kitCases(b, "kit-lifecycle-operation-complete")
    expect_identical(h$kitOperationLifecycleDecision(kitSyntheticOperationReceipt(complete, h), complete, b)$status,
        "lifecycle_complete_observed")
    missing_phase <- kitSyntheticOperationReceipt(complete, h, edit = function(r) { r$cases$member$operation$phase_status$ordered_result_stage <- "open"; r })
    expect_identical(h$kitOperationLifecycleDecision(missing_phase, complete, b)$failing_rule, "phase_timing")
    slow <- kitSyntheticOperationReceipt(complete, h, edit = function(r) { r$cases$member$operation$elapsed_seconds <- 301; r })
    expect_identical(h$kitOperationLifecycleDecision(slow, complete, b)$failing_rule, "elapsed_within_ceiling")
    fail <- h$kitCases(b, "kit-lifecycle-operation-fail")
    failed <- kitSyntheticOperationReceipt(fail, h, edit = function(r) {
        r$cases$member$operation$status <- "worker_failed"
        r$cases$member$operation$phase_status$ordered_result_stage <- "failed"
        r$cases$member$operation$head_unchanged <- TRUE
        r
    })
    expect_identical(h$kitOperationLifecycleDecision(failed, fail, b)$status, "lifecycle_failure_observed")
    expect_identical(h$kitOperationLifecycleDecision(kitSyntheticOperationReceipt(fail, h), fail, b)$failing_rule, "failure_recorded")
    moved <- failed
    moved$cases$member$operation$head_unchanged <- FALSE
    expect_identical(h$kitOperationLifecycleDecision(moved, fail, b)$failing_rule, "head_unchanged")
    between <- failed
    between$cases$member$operation$phase_status$ordered_result_stage <- "complete"
    expect_identical(h$kitOperationLifecycleDecision(between, fail, b)$failing_rule, "phase_timing")
})

test_that("the baseline journey decision needs both control members eligible and the second reproduced", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$journey$pair <- list(control = "core_018", candidate = NULL, reps = 2L, production_default = TRUE)
    cases <- h$kitCases(b, "kit-journey")
    receipt <- kitSyntheticReceipt(cases, h, edit = function(r) {
        r$cases[["control-2"]]$verifier_completion$reproduced <- TRUE
        r$cases[["control-2"]]$operation$sampled_anon_bytes <- 4.2e9
        r$cases[["control-1"]]$operation$sampled_anon_bytes <- 4e9
        r$cases[["default-1"]]$operation$populations_complete <- FALSE
        r
    })
    decision <- h$kitJourneyDecision(receipt, cases, b)
    expect_identical(decision$status, "baseline_recorded")
    expect_equal(decision$noise$peak_anon, 1.05)
    expect_true(decision$production_default$recorded)
    expect_false(decision$promotion_authorized)
    not_reproduced <- receipt
    not_reproduced$cases[["control-2"]]$verifier_completion$reproduced <- FALSE
    decision <- h$kitJourneyDecision(not_reproduced, cases, b)
    expect_identical(decision$status, "invalid")
    expect_identical(unname(decision$failing_rule), "not_reproduced")
    expect_null(decision$noise)
    partial <- receipt
    partial$cases[["control-1"]]$operation$populations_complete <- FALSE
    expect_identical(unname(h$kitJourneyDecision(partial, cases, b)$failing_rule), "populations_incomplete")
    slow_default <- receipt
    slow_default$cases[["default-1"]]$operation$elapsed_seconds <- 3601
    decision <- h$kitJourneyDecision(slow_default, cases, b)
    expect_identical(decision$status, "baseline_recorded")
    expect_false(decision$production_default$recorded)
})

test_that("canonical fingerprints drop run-specific ids, times and plot objects and reproduce across stores", {
    h <- kitHarness()
    reference <- function(id, sha, codec = "parquet_duckdb") list(entity_id = id, kind = "checkpoint",
        payload = list(codec = codec, sha256 = sha, relative_path = paste0("entities/", id, "/payload")))
    a <- list(normalized = reference(strrep("a", 64), strrep("1", 64)), timestamp = Sys.time(), path = "/run-1/x.tsv",
        plot = structure(list(), class = "gg"), table = data.frame(v = 1:3), nested = list(mtime = 1, value = 2))
    b <- list(normalized = reference(strrep("b", 64), strrep("1", 64)), timestamp = Sys.time() + 100, path = "/run-2/x.tsv",
        plot = structure(list(a = 1), class = "gg"), table = data.frame(v = 1:3), nested = list(mtime = 2, value = 2))
    expect_identical(h$kitCanonical(a), h$kitCanonical(b))
    expect_identical(h$kitCanonical(a)$normalized$sha256, strrep("1", 64))
    expect_identical(h$kitCanonical(a)$plot, "<gg>")
    expect_null(h$kitCanonical(a)$timestamp)
    expect_false(identical(h$kitCanonical(a), h$kitCanonical(c(a, list(extra = 1)))))
    fingerprint <- function(sha) list(class = "ProteinQuantitativeData", rows = 10L, summary_sha256 = "s",
        components = list(protein_quant_table = list(kind = "table", sha256 = sha)), support = list(), products = list())
    first <- list(step = list(fingerprint = fingerprint("x")), other = list(fingerprint = NULL))
    expect_length(h$kitReproduceMismatches(first, list(step = list(fingerprint = fingerprint("x")))), 0L)
    expect_identical(h$kitReproduceMismatches(first, list(step = list(fingerprint = fingerprint("y")))), "step/components/protein_quant_table")
    expect_identical(h$kitReproduceMismatches(first, list()), "step/missing")
    expect_true(h$kitStateFingerprintMatches(c(fingerprint("x"), list(extra = 1)), fingerprint("x")))
    expect_false(h$kitStateFingerprintMatches(fingerprint("x"), fingerprint("y")))
    files <- function(pdf, png, stamp = "20261010_184024", id = "3005b574e4ee6") data.frame(
        path = c("a.pdf", paste0("filtering_qc/", stamp, "/b.png"), paste0("filtering_qc/", stamp, "/.workflow-export/qvalue-", strrep("a", 40), "/journal.json"),
            paste0("filtered_session_", stamp, "_native-export-", id, ".msdia.tar"), "filtered_session_summary.txt", "fasta_metadata.RDS"),
        bytes = c(10, 20, 5, 100, 7, 9), sha256 = c(pdf, png, "j", "t", "s", "f"), stringsAsFactors = FALSE)
    same <- files("p1", "q", "20261010_200235", "29442991abd56")
    same$sha256[4:5] <- c("t2", "s2")
    expect_length(h$kitFileMismatches(files("p1", "q"), same), 0L)
    expect_identical(h$kitFileMismatches(files("p1", "q1"), files("p2", "q2")), "filtering_qc/<time>/b.png")
    expect_identical(h$kitFileMismatches(files("p1", "q"), transform(files("p1", "q"), sha256 = replace(sha256, 6, "f2"))), "fasta_metadata.RDS")
    expect_null(h$kitCanonical(list(processing_timestamp = Sys.time(), x = 1))$processing_timestamp)
})

test_that("a journey population is identified by its per-step checkpoint descriptor and materialised as a snapshot store", {
    h <- kitHarness()
    store <- file.path(withr::local_tempdir(), "store")
    dir.create(file.path(store, "entities", "dia-state-aaa"), recursive = TRUE)
    dir.create(file.path(store, "entities", "dia-state-bbb"))
    writeLines('{"schema":1}', file.path(store, "store.json"))
    saveRDS(list(entity_id = strrep("1", 64), kind = "checkpoint"), file.path(store, "entities", "dia-state-aaa", "descriptor.rds"))
    saveRDS(list(entity_id = strrep("2", 64), kind = "checkpoint"), file.path(store, "entities", "dia-state-bbb", "descriptor.rds"))
    # The live head has moved on to the second checkpoint.
    file.copy(file.path(store, "entities", "dia-state-bbb", "descriptor.rds"), file.path(store, "dia-state-head.rds"))
    step_head <- digest::digest(file = file.path(store, "entities", "dia-state-aaa", "descriptor.rds"), algo = "sha256")
    expect_identical(h$kitFindHeadEntity(store, step_head), "dia-state-aaa")
    expect_null(h$kitFindHeadEntity(store, strrep("0", 64)))
    journey <- list(store = store, head_sha256 = step_head, head_entity = "dia-state-aaa", entity_id = strrep("5", 64), state = "imputed")
    expect_null(h$kitPopulationInvalid(journey))
    # Without the head entity the stale live head is refused, as the audit found (PH-01).
    expect_identical(h$kitPopulationInvalid(journey[setdiff(names(journey), "head_entity")]), "head_hash_mismatch")
    expect_identical(h$kitPopulationInvalid(c(journey[setdiff(names(journey), "head_sha256")], list(head_sha256 = strrep("0", 64)))), "head_hash_mismatch")
    expect_identical(h$kitPopulationInvalid(list(store = file.path(store, "missing"))), "store_missing")
    live <- list(store = store, head_sha256 = digest::digest(file = file.path(store, "dia-state-head.rds"), algo = "sha256"))
    expect_null(h$kitPopulationInvalid(live))
    expect_identical(h$kitPopulationStore(live, file.path(store, "..", "unused")), store)
    snapshot <- h$kitPopulationStore(journey, file.path(dirname(store), "member", "population-store"))
    expect_identical(basename(snapshot), ".workflow-entities")
    expect_true(dir.exists(file.path(snapshot, "entities", "dia-state-bbb")))
    expect_identical(digest::digest(file = file.path(snapshot, "dia-state-head.rds"), algo = "sha256"), step_head)
    expect_identical(readRDS(file.path(snapshot, "dia-state-head.rds"))$entity_id, strrep("1", 64))
    expect_true(file.exists(file.path(snapshot, "store.json")))
    # The evidence store is untouched.
    expect_identical(digest::digest(file = file.path(store, "dia-state-head.rds"), algo = "sha256"), live$head_sha256)
    expect_error(h$kitPopulationStore(journey, dirname(snapshot)), "dir.exists")
    b <- kitSyntheticBindings()
    b$populations$journey_imputation <- c(journey, list(rows = 10L))
    b$tickets[["009"]] <- list(operation = "minimum_peptides_per_sample", population = "journey_imputation", control = "core_018", candidate = "kit_006")
    case <- h$kitCases(b, "kit-009")[[1L]]
    expect_identical(case$input$head_entity, "dia-state-aaa")
    expect_identical(case$input$state, "imputed")
})

test_that("the production-default member is recorded when it completes or ends at a protocol guard with its step peaks", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$journey$pair <- list(control = "core_018", candidate = NULL, reps = 2L, production_default = TRUE)
    cases <- h$kitCases(b, "kit-journey")
    receipt <- kitSyntheticReceipt(cases, h, edit = function(r) {
        r$cases[["control-2"]]$verifier_completion$reproduced <- TRUE
        r$cases[["default-1"]]$operation$populations_complete <- FALSE
        r
    })
    expect_true(h$kitJourneyDecision(receipt, cases, b)$production_default$recorded)
    expect_identical(h$kitJourneyDecision(receipt, cases, b)$production_default$outcome, "complete")
    guarded <- receipt
    guarded$cases[["default-1"]]$operation$status <- "memory_guard"
    guarded$cases[["default-1"]]$operation$outer_guard <- "memory_guard"
    guarded$cases[["default-1"]]$operation$peak_bytes <- 16.6e9
    guarded$cases[["default-1"]]$operation$resident_end_bytes <- NA_real_
    guarded$cases[["default-1"]]$operation$phase_status$imputation <- "failed"
    guarded$cases[["default-1"]]$operation$phase_status[10:18] <- NULL
    guarded$cases[["default-1"]]$verifier_completion <- list(accepted = FALSE)
    decision <- h$kitJourneyDecision(guarded, cases, b)
    expect_identical(decision$status, "baseline_recorded")
    expect_true(decision$production_default$recorded)
    expect_identical(decision$production_default$outcome, "memory_guard")
    expect_null(decision$production_default$failing_rule)
    expect_identical(h$kitDefaultOutcome(guarded$cases[["default-1"]]$operation), "memory_guard")
    gap <- guarded
    gap$cases[["default-1"]]$operation$max_sample_interval <- 0.6
    expect_identical(h$kitJourneyDecision(gap, cases, b)$production_default$failing_rule, "measurement_valid")
    swapped <- guarded
    swapped$cases[["default-1"]]$operation$swap_bytes <- 4096
    expect_false(h$kitJourneyDecision(swapped, cases, b)$production_default$recorded)
    crashed <- guarded
    crashed$cases[["default-1"]]$operation$status <- "journey_failed"
    crashed$cases[["default-1"]]$operation$outer_guard <- NULL
    expect_identical(h$kitJourneyDecision(crashed, cases, b)$production_default$failing_rule, "status_complete_or_guard")
    no_peaks <- guarded
    no_peaks$cases[["default-1"]]$operation$phase_peak_anon_bytes$import <- NA_real_
    expect_identical(h$kitJourneyDecision(no_peaks, cases, b)$production_default$failing_rule, "phase_peaks_recorded")
})

test_that("the operation pair decision needs every member eligible and each candidate not worse than 1.02 x its control", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    cases <- h$kitCases(b, "kit-006")
    receipt <- kitSyntheticOperationReceipt(cases, h, edit = function(r) {
        for (id in names(r$cases)) {
            r$cases[[id]]$operation$sampled_anon_bytes <- 4e9
            r$cases[[id]]$operation$science_execution <- list(admission = list(batches = 3L))
        }
        r$cases[["candidate-1"]]$operation$peak_bytes <- 5e9 * 1.01
        r$cases[["candidate-2"]]$operation$peak_bytes <- 5e9 * 0.9
        r$cases[["candidate-2"]]$operation$science_execution$admission$batches <- 5L
        # A control build's fixed-three policy records count_batch_sizes only (KIT-006 correction 1).
        r$cases[["control-2"]]$operation$science_execution <- list(admission = list(count_batch_sizes = c(100L, 100L, 100L)))
        r
    })
    decision <- h$kitOperationDecision(receipt, cases, b)
    expect_identical(decision$status, "accepted")
    expect_null(decision$failing_rule)
    expect_equal(decision$pairs[["1"]]$peak_ratio, 1.01)
    expect_false(decision$pairs[["1"]]$reduction)
    expect_true(decision$pairs[["2"]]$reduction)
    expect_identical(decision$pairs[["2"]]$batches$candidate, 5L)
    expect_identical(decision$pairs[["2"]]$batches$control, 3L)
    expect_identical(decision$pairs[["2"]]$operation, "minimum_peptides_per_sample")
    expect_false(decision$promotion_authorized)
    worse <- receipt
    worse$cases[["candidate-1"]]$operation$peak_bytes <- 5e9 * 1.03
    rejected <- h$kitOperationDecision(worse, cases, b)
    expect_identical(rejected$status, "rejected")
    expect_identical(rejected$failing_rule, "not_worse")
    anon <- receipt
    anon$cases[["candidate-2"]]$operation$sampled_anon_bytes <- 4e9 * 1.05
    expect_identical(h$kitOperationDecision(anon, cases, b)$status, "rejected")
    unverified <- receipt
    unverified$cases[["control-2"]]$verifier_completion$accepted <- FALSE
    invalid <- h$kitOperationDecision(unverified, cases, b)
    expect_identical(invalid$status, "invalid")
    expect_identical(unname(invalid$failing_rule), "verifier_accepted")
    gap <- receipt
    gap$cases[["candidate-1"]]$operation$max_sample_interval <- 0.3
    expect_identical(unname(h$kitOperationDecision(gap, cases, b)$failing_rule), "measurement_valid")
})

test_that("the sample-filter verifier fingerprints a result exactly as the CORE native-oracle recipe did", {
    skip_if_not_installed("fastduckdb")
    h <- kitHarness()
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small", if (exists("workflow_test_root")) workflow_test_root else here::here()))
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, basename(root), "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "input")
    value <- computeDIAReferenceSampleFilter(manager, getDIAStateReference(manager),
        list(peptides_per_sample_cutoff = 3L, inclusion_list = NULL), output_mode = "physical")
    actual <- h$kitFingerprintFilterResult(manager, value$checkpoint, api)
    # The CORE recipe (tools/profiling/workflow_perf106.R, verifyPerf106Result), applied here verbatim.
    reference <- value$checkpoint$components$peptide_data
    scope <- newWorkflowTableReadScope(store, reference, reader = "fastduckdb")
    hashes <- character()
    for (offset in seq.int(0L, max(0L, scope$metadata$rows - 1L), by = 10000L)) {
        page <- readSavedWorkflowTable(store, reference, scope$metadata$columns, offset = offset,
            limit = min(10000L, scope$metadata$rows - offset), scope = scope)$data
        hashes[[length(hashes) + 1L]] <- digest::digest(page, algo = "sha256")
    }
    rows <- scope$metadata$rows
    columns <- scope$metadata$columns
    releaseWorkflowTableReadScope(scope)
    design <- materializeDIATable(manager$.context, value$checkpoint$components$design_matrix)
    support <- lapply(value$checkpoint$support$filterMinNumPeptidesPerSample, function(r) materializeDIATable(manager$.context, r))
    expect_identical(actual$page_hashes, hashes)
    expect_identical(actual$rows, rows)
    expect_identical(actual$samples, nrow(design))
    expect_identical(actual$columns, columns)
    expect_identical(actual$design_sha256, digest::digest(design, algo = "sha256"))
    expect_identical(actual$summary_sha256, digest::digest(value$checkpoint$slots$args$filterMinNumPeptidesPerSample$filter_summary, algo = "sha256"))
    expect_identical(actual$support_sha256, digest::digest(support, algo = "sha256"))
    expect_gt(length(hashes), 0L)
})

test_that("a ticket journey pair runs both builds as users run them and is decided per completed step against a control that may guard", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$builds$kit_011 <- list(library = "/lib/e", build = strrep("e", 64), exact_sum = "sum-e")
    b$tickets[["011"]] <- list(kind = "journey", control = "core_018", candidate = "kit_011", record = FALSE)
    cases <- h$kitCases(b, "kit-journey-011")
    expect_identical(vapply(cases, `[[`, "", "id"), c("control-1", "candidate-1"))
    expect_identical(vapply(cases, `[[`, "", "library"), c("/lib/n", "/lib/e"))
    expect_false(cases[[1]]$record)
    expect_identical(cases[[2]]$selector, "kit:journey-011:candidate:1")
    expect_identical(h$kitSelectorStage("kit:journey-011:candidate:1"), list(stage = "kit-journey-011", id = "candidate-1"))
    expect_identical(h$kitMemberSeconds("kit-journey-011", b), 3600L)
    expect_error(h$kitCases(b, "kit-journey-012"), "journey pair")
    steps <- h$kitStepIds()
    receipt <- kitSyntheticReceipt(cases, h, edit = function(r) {
        # The control (in-memory default) guards in imputation after eight steps; the candidate completes lower.
        control <- r$cases[["control-1"]]$operation
        control$status <- "memory_guard"; control$outer_guard <- "memory_guard"; control$peak_bytes <- 16.6e9
        control$sampled_anon_bytes <- 6e9
        control$resident_end_bytes <- NA_real_; control$populations_complete <- FALSE
        control$phase_status[steps[10:18]] <- NULL; control$phase_status$imputation <- "failed"
        control$phase_peak_anon_bytes <- as.list(stats::setNames(c(2.7, 3.5, 6.4, 6.6, 9.7, 9.3, 9.3, 9.3, 14.3) * 1e9, steps[1:9]))
        r$cases[["control-1"]]$operation <- control
        r$cases[["control-1"]]$verifier_completion <- list(accepted = TRUE)
        candidate <- r$cases[["candidate-1"]]$operation
        candidate$sampled_anon_bytes <- 5.7e9
        # Import and design do not differ between the paths; the data steps do.
        candidate$phase_peak_anon_bytes <- as.list(stats::setNames(c(2.7e9, 3.5e9, rep(3e9, length(steps) - 2L)), steps))
        r$cases[["candidate-1"]]$operation <- candidate
        r
    })
    decision <- h$kitJourneyTicketDecision(receipt, cases, b)
    expect_identical(decision$status, "accepted")
    expect_identical(decision$control_outcome, "memory_guard")
    expect_true(decision$control_peak_is_lower_bound)
    expect_identical(decision$steps_compared, steps[1:8])
    # The control's anonymous peak at the guard (6e9 in the synthetic receipt) is a lower bound; 5.7 / 6 = 0.95.
    expect_equal(decision$whole_slice_anon_ratio, 5.7 / 6, tolerance = 1e-6)
    expect_equal(decision$step_ratios$rollup, 3 / 9.7, tolerance = 1e-6)
    expect_false(decision$promotion_authorized)
    not_reduced <- receipt
    not_reduced$cases[["candidate-1"]]$operation$sampled_anon_bytes <- 5.8e9
    rejected <- h$kitJourneyTicketDecision(not_reduced, cases, b)
    expect_identical(rejected$status, "rejected")
    expect_identical(rejected$failing_rule, "required_reduction")
    worse_step <- receipt
    worse_step$cases[["candidate-1"]]$operation$phase_peak_anon_bytes$import <- 2.7e9 * 1.03
    expect_identical(h$kitJourneyTicketDecision(worse_step, cases, b)$failing_rule, "step_not_worse")
    guarded <- receipt
    guarded$cases[["candidate-1"]]$operation$status <- "memory_guard"
    expect_identical(h$kitJourneyTicketDecision(guarded, cases, b)$status, "invalid")
    unverified <- receipt
    unverified$cases[["candidate-1"]]$verifier_completion$accepted <- FALSE
    expect_identical(unname(h$kitJourneyTicketDecision(unverified, cases, b)$failing_rule), "verifier_accepted")
    complete_control <- receipt
    complete_control$cases[["control-1"]] <- receipt$cases[["candidate-1"]]
    complete_control$cases[["control-1"]]$operation$sampled_anon_bytes <- 6e9
    complete_control$cases[["control-1"]]$operation$phase_peak_anon_bytes <- as.list(stats::setNames(c(2.7e9, 3.5e9, rep(4e9, length(steps) - 2L)), steps))
    both <- h$kitJourneyTicketDecision(complete_control, cases, b)
    expect_identical(both$status, "accepted")
    expect_equal(both$whole_slice_anon_ratio, 0.95)
    expect_equal(both$resident_end_ratio, 1)
    oracles <- h$kitOracleJourneyPopulations(b)
    expect_identical(names(oracles), steps)
    expect_null(oracles$import$fingerprint)
})
