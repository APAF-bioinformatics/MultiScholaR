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
                list(operation = "minimum_peptides_per_sample", population = "full_primary"))),
            `008` = list(kind = "gui", build = "kit_006", corpus = "screen", reference = "kit-lifecycle-complete",
                operations = c("precursor_rollup", "minimum_peptides_per_sample"), ceiling_seconds = 600L)))
}

# KIT-092 (protocol v1.7): a single journey is refused at binding without an oracle for every step it runs that leaves a
# state, so a synthetic single journey binds one per such step.
kitSyntheticOracles <- function(b, h, route, prefix) {
    for (step in h$kitJourneySteps(route)) if (!is.null(step$state) && !identical(step$state, "initial"))
        b$oracles[[paste0(prefix, "_", step$id)]] <- list(class = "Synthetic", rows = 1L,
            components = list(table = list(kind = "table", sha256 = "x")))
    b
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
        "correlation", "session_export", "da", "report")) # KIT-065: the report closes every route
    steps <- h$kitJourneySteps()
    expect_identical(steps[[1]]$state, "initial")
    expect_identical(steps[[2]]$state, "raw_data_s4")
    expect_identical(steps[[7]]$state, "sample_filtered")
    expect_null(steps[[17]]$state)
})

test_that("KIT-021: the protein route's step list, a whole-journey single member with a suffixed key and in-process oracles", {
    h <- kitHarness()
    # The DIA list is unchanged; the protein route has no peptide stage and no IQ rollup and its design commits protein_s4_initial.
    expect_identical(h$kitStepIds("DIA"), h$kitStepIds())
    protein <- h$kitStepIds("protein")
    expect_identical(protein, c("import", "design", "accession", "protein_intensity", "duplicates", "protein_replicate",
        "normalisation_ruv", "correlation", "session_export", "da", "report"))
    expect_identical(h$kitJourneySteps("protein")[[2]]$state, "protein_s4_initial")
    expect_error(h$kitStepIds("limpa"), "route")
    b <- kitSyntheticBindings()
    b$builds$kit_021 <- list(library = "/lib/p", build = strrep("p", 64), exact_sum = "sum-p")
    b$corpora$lfq_300 <- b$corpora$full
    b$corpora$lfq_300$report_file <- "proteinGroups.txt"; b$corpora$lfq_300$format <- "maxquant"; b$corpora$lfq_300$route <- "protein"
    b$harness$step_map_protein_sha256 <- strrep("m", 64)
    b$tickets[["021"]] <- list(kind = "single_journey", build = "kit_021", corpus = "lfq_300", route = "protein", whole_journey = TRUE,
        oracles = "lfq_300", ceiling_seconds = 14400L)
    b$tickets[["021t"]] <- list(kind = "single_journey", build = "kit_021", corpus = "lfq_300", route = "protein", whole_journey = TRUE,
        oracles = "lfq_300", ceiling_seconds = 1800L)
    b <- kitSyntheticOracles(b, h, "protein", "lfq_300")
    cases <- h$kitCases(b, "kit-journey-021")
    expect_length(cases, 1L)
    case <- cases[[1L]]
    expect_identical(case$route, "protein")
    expect_identical(case$oracles, "lfq_300")
    expect_null(case$stop_after)
    expect_identical(case$selector, "kit:journey-021:candidate:1")
    expect_identical(h$kitSelectorStage("kit:journey-021t:candidate:1"), list(stage = "kit-journey-021t", id = "candidate-1"))
    expect_identical(h$kitCases(b, "kit-journey-021t")[[1L]]$selector, "kit:journey-021t:candidate:1")
    # Each key binds its own ceiling; the protein step map is part of the binding; a protein route without the map refuses.
    expect_identical(h$kitMemberSeconds("kit-journey-021", b), 14400L)
    expect_identical(h$kitMemberSeconds("kit-journey-021t", b), 1800L)
    other <- b; other$harness$step_map_protein_sha256 <- strrep("n", 64)
    expect_false(identical(h$kitCases(other, "kit-journey-021")[[1L]]$binding, case$binding))
    none <- b; none$harness$step_map_protein_sha256 <- NULL
    expect_error(h$kitCases(none, "kit-journey-021"), "protein step map")
    # A single journey without a stop step must say so; the DIA stop-after form is unchanged.
    bare <- b; bare$tickets[["021"]]$whole_journey <- NULL
    expect_error(h$kitCases(bare, "kit-journey-021"), "build or stop_after")
    # The oracles are looked up under the bound prefix for the route's steps.
    b$oracles$lfq_300_design <- list(class = "ProteinQuantitativeData", components = list(protein_quant_table = "x"))
    oracles <- h$kitOracleJourneyPopulations(b, "protein", "lfq_300")
    expect_identical(names(oracles), protein)
    expect_identical(oracles$design$fingerprint$class, "ProteinQuantitativeData")
    expect_null(oracles$import$fingerprint)
    # The decision takes the route's whole step list when there is no stop step.
    receipt <- kitSyntheticReceipt(cases, h, edit = function(r) {
        op <- r$cases[["candidate-1"]]$operation
        op$status <- "complete"; op$populations_complete <- TRUE
        op$phase_status <- stats::setNames(as.list(rep("complete", length(protein))), protein)
        op$phase_peak_anon_bytes <- stats::setNames(as.list(rep(2e9, length(protein))), protein)
        op$phase_seconds <- stats::setNames(as.list(rep(10, length(protein))), protein)
        op$sampled_anon_bytes <- 2e9; op$resident_end_bytes <- 1.5e9
        r$cases[["candidate-1"]]$operation <- op
        r$cases[["candidate-1"]]$verifier_completion <- list(accepted = TRUE)
        r
    })
    decision <- h$kitSingleJourneyDecision(receipt, cases, b)
    expect_identical(decision$status, "accepted")
    expect_identical(decision$steps, protein)
    expect_identical(decision$ticket, "021")
    expect_null(decision$stop_after)
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
    lifecycle_build$lifecycle$build <- "kit_006"
    expect_identical(h$kitCases(lifecycle_build, "kit-lifecycle-operation-complete")[[1]]$library, "/lib/k")
    expect_identical(h$kitCases(lifecycle_build, "kit-lifecycle-complete")[[1]]$library, "/lib/k")
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

test_that("KIT-068: the lipidomics and metabolomics routes, their in-memory lifecycle members and decisions", {
    h <- kitHarness()
    ids <- c("import", "design", "intensity", "duplicates", "qc_finalize", "normalisation", "correlation", "session_export", "da", "report")
    expect_identical(h$kitStepIds("lipidomics"), ids)
    expect_identical(h$kitStepIds("metabolomics"), ids)
    lipid <- h$kitJourneySteps("lipidomics")
    expect_null(lipid[[1]]$state)
    expect_identical(lipid[[2]]$state, "lipid_raw_data_s4")
    expect_identical(lipid[[5]]$state, "lipid_qc_complete")
    expect_identical(lipid[[6]]$state, "lipid_ruv_corrected")
    expect_identical(h$kitJourneySteps("lipidomics", "skip")[[6]]$state, "lipid_normalized")
    expect_identical(lipid[[9]]$state, "lipid_correlation_filtered")
    expect_identical(vapply(lipid, `[[`, "", "kind"), c("store", "store", "omics_qc", "omics_qc", "omics_qc", "omics_norm", "omics_norm", "store", "da", "report"))
    metab <- h$kitJourneySteps("metabolomics")
    expect_identical(metab[[2]]$state, "metab_raw_data_s4")
    expect_identical(metab[[7]]$state, "metab_correlation_filtered")
    expect_identical(h$kitStepIds("DIA"), h$kitStepIds())
    expect_error(h$kitStepIds("lipid"), "route")
    # The omics lifecycle stages: one in-memory member per route on the bound build, the route's screen corpus, the
    # route's step map in the binding, the fault carried on the case.
    b <- kitSyntheticBindings()
    file <- function(name) list(path = paste0("/corpus/", name), bytes = 1, sha256 = strrep("a", 64))
    b$corpora$lipid_screen <- list(profile = "lipid_screen", samples = 24L, route = "lipidomics", format = "lipidsearch", vendor = "lipidsearch",
        report_file = "lipidsearch_lcms_pos.txt", assay_files = list(LCMS_Pos = "lipidsearch_lcms_pos.txt", LCMS_Neg = "lipidsearch_lcms_neg.txt"),
        files = list(`lipidsearch_lcms_pos.txt` = file("lipidsearch_lcms_pos.txt"), `lipidsearch_lcms_neg.txt` = file("lipidsearch_lcms_neg.txt"), `design.tsv` = file("design.tsv")))
    b$corpora$metab_screen <- list(profile = "metab_screen", samples = 24L, route = "metabolomics", format = "msdial", vendor = "msdial",
        report_file = "msdial_lcms_pos.txt", assay_files = list(LCMS_Pos = "msdial_lcms_pos.txt", GCMS = "msdial_gcms.txt"),
        files = list(`msdial_lcms_pos.txt` = file("msdial_lcms_pos.txt"), `msdial_gcms.txt` = file("msdial_gcms.txt"), `design.tsv` = file("design.tsv")))
    b$harness <- list(step_map_lipid_sha256 = strrep("l", 64), step_map_metab_sha256 = strrep("m", 64))
    b$lifecycle_omics <- list(ceiling_seconds = 1800L, build = "kit_006",
        fault = list(name = "guard_after_qc_finalize", after_step = "qc_finalize", failed_step = "normalisation"),
        normalisation = list(apply_itsd = TRUE, norm_method = "cyclicloess", ruv_mode = "automatic", correlation_threshold = 0.5))
    # KIT-092 (protocol v1.7): in memory only where declared; the KIT-068 builds had no store.
    b$lifecycle_omics$in_memory <- TRUE
    complete <- h$kitCases(b, "kit-lifecycle-lipid")
    expect_length(complete, 1L)
    case <- complete[[1L]]
    expect_identical(case$id, "member")
    expect_identical(case$route, "lipidomics")
    expect_identical(case$corpus, "lipid_screen")
    expect_identical(case$library, "/lib/k")
    expect_false(case$expect_failure)
    expect_false(case$record)
    expect_null(case$oracles)
    expect_identical(case$fault_spec$failed_step, "normalisation")
    expect_identical(case$normalisation$ruv_mode, "automatic")
    expect_identical(case$selector, "kit:lifecycle:lipid:complete")
    expect_identical(h$kitSelectorStage("kit:lifecycle:lipid:complete"), list(stage = "kit-lifecycle-lipid", id = "member"))
    expect_identical(h$kitSelectorStage("kit:lifecycle:metab:fail"), list(stage = "kit-lifecycle-metab-fail", id = "member"))
    expect_identical(h$kitMemberSeconds("kit-lifecycle-lipid", b), 1800L)
    expect_identical(h$kitMemberSeconds("kit-lifecycle-metab-fail", b), 1800L)
    expect_identical(h$kitMemberSeconds("kit-lifecycle-complete", b), 300L)
    fail <- h$kitCases(b, "kit-lifecycle-metab-fail")[[1L]]
    expect_true(fail$expect_failure)
    expect_true(fail$record) # the launch record permits the guard injection (attempt 1 of the lipid failure member)
    expect_identical(fail$fault, "guard_after_qc_finalize")
    expect_identical(fail$route, "metabolomics")
    expect_identical(fail$selector, "kit:lifecycle:metab:fail")
    # The step map is part of the binding; a route without its map refuses.
    other <- b; other$harness$step_map_lipid_sha256 <- strrep("n", 64)
    expect_false(identical(h$kitCases(other, "kit-lifecycle-lipid")[[1L]]$binding, case$binding))
    none <- b; none$harness$step_map_metab_sha256 <- NULL
    expect_error(h$kitCases(none, "kit-lifecycle-metab"), "metab step map")
    unbound <- b; unbound$lifecycle_omics <- NULL
    expect_error(h$kitCases(unbound, "kit-lifecycle-lipid"), "omics lifecycle")
    # The DIA lifecycle member is unchanged in shape.
    expect_identical(h$kitCases(b, "kit-lifecycle-complete")[[1L]]$selector, "kit:lifecycle:complete")
    # Decisions on the in-memory member: the route's steps, no populations rule, the ceiling of the omics stage; the
    # failure member's head rule is not applicable without a store and is recorded, never a pass by default elsewhere.
    receipt_for <- function(cases, edit = NULL) {
        receipt <- list(cases = list())
        for (case in cases) {
            operation <- list(status = "complete", swap_bytes = 0, peak_bytes = 4e9, elapsed_seconds = 900,
                max_sample_interval = 0.06, cleanup = list(slice_stopped = TRUE), resident_end_bytes = 1e9,
                populations_complete = FALSE, head_unchanged = NA,
                phase_status = as.list(stats::setNames(rep("complete", length(ids)), ids)),
                phase_peak_anon_bytes = as.list(stats::setNames(rep(2e9, length(ids)), ids)))
            receipt$cases[[case$id]] <- list(passed = TRUE, operation = operation, verifier_completion = list(accepted = TRUE))
        }
        if (is.function(edit)) receipt <- edit(receipt)
        receipt
    }
    decision <- h$kitLifecycleDecision(receipt_for(complete), complete, b)
    expect_identical(decision$status, "lifecycle_complete_observed")
    expect_true(decision$in_memory)
    expect_identical(decision$route, "lipidomics")
    slow <- receipt_for(complete, edit = function(r) { r$cases$member$operation$elapsed_seconds <- 1801; r })
    expect_identical(h$kitLifecycleDecision(slow, complete, b)$failing_rule, "elapsed_within_ceiling")
    missing_step <- receipt_for(complete, edit = function(r) { r$cases$member$operation$phase_status$report <- NULL; r })
    expect_identical(h$kitLifecycleDecision(missing_step, complete, b)$failing_rule, "phase_timing")
    fail_cases <- list(fail)
    failed <- receipt_for(fail_cases, edit = function(r) {
        r$cases$member$operation$status <- "journey_failed"
        r$cases$member$operation$phase_status$normalisation <- "failed"
        r$cases$member$operation$guard <- list(reason = "injected_guard_stop")
        r
    })
    failure <- h$kitLifecycleDecision(failed, fail_cases, b)
    expect_identical(failure$status, "lifecycle_failure_observed")
    expect_match(failure$head_rule, "not applicable")
    expect_identical(h$kitLifecycleDecision(receipt_for(fail_cases), fail_cases, b)$failing_rule, "failure_recorded")
    wrong_step <- failed; wrong_step$cases$member$operation$phase_status$normalisation <- "complete"
    expect_identical(h$kitLifecycleDecision(wrong_step, fail_cases, b)$failing_rule, "phase_timing")
    # The DIA failure member still needs the head rule.
    dia_fail <- h$kitCases(b, "kit-lifecycle-fail")
    dia_failed <- kitSyntheticReceipt(dia_fail, h, edit = function(r) {
        r$cases$member$operation$status <- "journey_failed"
        r$cases$member$operation$phase_status$sample <- "failed"
        r$cases$member$operation$guard <- list(reason = "injected_guard_stop")
        r
    })
    expect_identical(h$kitLifecycleDecision(dia_failed, dia_fail, b)$failing_rule, "head_unchanged")
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
    # The KIT-002 record field `declaration` is storage format, not science (KIT-011 correction 1).
    with_declaration <- list(data = list(schema = "dia.component_checkpoint.v2", class = "ProteinQuantitativeData",
        declaration = list(class = "ProteinQuantitativeData", version = 1L), slots = list(sample_id = "Run")))
    without <- with_declaration; without$data$declaration <- NULL
    expect_identical(h$kitCanonical(with_declaration), h$kitCanonical(without))
    expect_false(identical(h$kitCanonical(a), h$kitCanonical(c(a, list(extra = 1)))))
    fingerprint <- function(sha) list(class = "ProteinQuantitativeData", rows = 10L, summary_sha256 = "s",
        components = list(protein_quant_table = list(kind = "table", sha256 = sha)), support = list(), products = list())
    first <- list(step = list(fingerprint = fingerprint("x")), other = list(fingerprint = NULL))
    expect_length(h$kitReproduceMismatches(first, list(step = list(fingerprint = fingerprint("x")))), 0L)
    expect_identical(h$kitReproduceMismatches(first, list(step = list(fingerprint = fingerprint("y")))), "step/components/protein_quant_table")
    expect_identical(h$kitReproduceMismatches(first, list()), "step/missing")
    # Operation identity products differ across builds by construction; a cross-build check reports them.
    with_identity <- function(sha, suffix) {
        products <- list(normalization_science = list(sha256 = "n"))
        products[[paste0("dia-normalization-ruv:", suffix)]] <- list(sha256 = sha)
        list(step = list(fingerprint = list(class = "ProteinQuantitativeData", rows = 10L, summary_sha256 = "s",
            components = list(protein_quant_table = list(kind = "table", sha256 = "x")), support = list(), products = products)))
    }
    # The missing oracle product mismatches; the other build's identity product rides along as an additional product (KIT-017).
    expect_identical(as.character(h$kitReproduceMismatches(with_identity("a", strrep("1", 12)), with_identity("b", strrep("2", 12)))), "step/products/names")
    expect_length(h$kitReproduceMismatches(with_identity("a", strrep("1", 12)), with_identity("b", strrep("2", 12)),
        ignore_products = h$kitOperationIdentityProducts), 0L)
    science_changed <- with_identity("b", strrep("2", 12))
    science_changed$step$fingerprint$products$normalization_science$sha256 <- "m"
    expect_identical(h$kitReproduceMismatches(with_identity("a", strrep("1", 12)), science_changed,
        ignore_products = h$kitOperationIdentityProducts), "step/products/normalization_science")
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
    expect_identical(vapply(cases, `[[`, "", "corpus"), c("full", "full"))
    # KIT-017: a pair bound to the screen corpus is its lifecycle run (both members on the screen corpus, not counted).
    screen <- b
    screen$tickets[["011"]]$corpus <- "screen"
    expect_identical(vapply(h$kitCases(screen, "kit-journey-011"), `[[`, "", "corpus"), c("screen", "screen"))
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
    worse_step$cases[["candidate-1"]]$operation$phase_peak_anon_bytes$qvalue <- 6.4e9 * 1.03
    expect_identical(h$kitJourneyTicketDecision(worse_step, cases, b)$failing_rule, "step_not_worse")
    # Import and design are reported, not judged (user decision 2026-10-11).
    noisy_import <- receipt
    noisy_import$cases[["candidate-1"]]$operation$phase_peak_anon_bytes$import <- 2.7e9 * 1.07
    reported <- h$kitJourneyTicketDecision(noisy_import, cases, b)
    expect_identical(reported$status, "accepted")
    expect_identical(reported$steps_guarded, steps[3:8])
    expect_equal(reported$step_ratios$import, 1.07, tolerance = 1e-6)
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

test_that("a ticket journey pair bound with the resident-end rule is decided on resident_end_bytes with the peak not worse (KIT-012)", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$builds$kit_012 <- list(library = "/lib/t", build = strrep("t", 64), exact_sum = "sum-t")
    b$tickets[["012"]] <- list(kind = "journey", control = "kit_006", candidate = "kit_012", record = FALSE, rule = "resident_end")
    cases <- h$kitCases(b, "kit-journey-012")
    expect_identical(vapply(cases, `[[`, "", "library"), c("/lib/k", "/lib/t"))
    expect_identical(cases[[2]]$selector, "kit:journey-012:candidate:1")
    steps <- h$kitStepIds()
    receipt <- kitSyntheticReceipt(cases, h, edit = function(r) {
        control <- r$cases[["control-1"]]$operation
        control$sampled_anon_bytes <- 5.5e9; control$resident_end_bytes <- 3.0e9
        r$cases[["control-1"]]$operation <- control
        candidate <- r$cases[["candidate-1"]]$operation
        candidate$sampled_anon_bytes <- 5.6e9; candidate$resident_end_bytes <- 2.85e9
        r$cases[["candidate-1"]]$operation <- candidate
        r
    })
    decision <- h$kitJourneyTicketDecision(receipt, cases, b)
    expect_identical(decision$status, "accepted")
    expect_identical(decision$rule, "resident_end")
    expect_equal(decision$resident_end_ratio, 0.95)
    expect_equal(decision$whole_slice_anon_ratio, 5.6 / 5.5, tolerance = 1e-6)
    expect_identical(decision$steps_guarded, steps[3:18])
    not_reduced <- receipt
    not_reduced$cases[["candidate-1"]]$operation$resident_end_bytes <- 2.86e9
    expect_identical(h$kitJourneyTicketDecision(not_reduced, cases, b)$failing_rule, "required_reduction")
    worse_peak <- receipt
    worse_peak$cases[["candidate-1"]]$operation$sampled_anon_bytes <- 5.5e9 * 1.03
    expect_identical(h$kitJourneyTicketDecision(worse_peak, cases, b)$failing_rule, "whole_slice_not_worse")
    worse_step <- receipt
    worse_step$cases[["candidate-1"]]$operation$phase_peak_anon_bytes$da <- 3e9 * 1.03
    expect_identical(h$kitJourneyTicketDecision(worse_step, cases, b)$failing_rule, "step_not_worse")
    # A control that ended at a guard has no resident-end reading: the pair is invalid, never rejected.
    guarded <- receipt
    guarded$cases[["control-1"]]$operation$status <- "memory_guard"
    guarded$cases[["control-1"]]$operation$outer_guard <- "memory_guard"
    guarded$cases[["control-1"]]$operation$resident_end_bytes <- NA_real_
    guarded$cases[["control-1"]]$operation$populations_complete <- FALSE
    invalid <- h$kitJourneyTicketDecision(guarded, cases, b)
    expect_identical(invalid$status, "invalid")
    expect_identical(invalid$failing_rule, "resident_end_unobserved")
    # The KIT-011 rule is the default when a ticket binds no rule.
    b$tickets[["012"]]$rule <- NULL
    expect_identical(h$kitJourneyTicketDecision(receipt, cases, b)$rule, "whole_slice_anon")
    expect_identical(h$kitJourneyTicketDecision(receipt, cases, b)$failing_rule, "required_reduction")
})

test_that("an operation pair may bind each arm's registry variant, its parameters and a journey oracle (KIT-009)", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$builds$kit_009 <- list(library = "/lib/i", build = strrep("i", 64), exact_sum = "sum-i")
    b$populations$journey_replicate <- list(step = "replicate", state = "replicate_filtered", store = "/store/journey",
        head_sha256 = strrep("7", 64), head_entity = "dia-state-head7", entity_id = strrep("8", 64), rows = 5515154L)
    b$oracles$journey_imputation <- list(class = "PeptideQuantitativeData", rows = 11826409L,
        components = list(peptide_data = list(sha256 = "p"), design_matrix = list(sha256 = "d")),
        support = list(`peptideMissingValueImputation/support_table` = list(sha256 = "s")), products = list(qc_progress = list(sha256 = "q")),
        summary_sha256 = "sum")
    b$tickets[["009"]] <- list(control = "kit_009", candidate = "kit_009",
        pairs = list(list(operation = "peptide_imputation", population = "journey_replicate", oracle = "journey_imputation",
            parameters = list(proportion_missing_values = 0.5, exclusion_column = NULL),
            variants = list(control = "grouped", candidate = "partitioned"))))
    cases <- h$kitCases(b, "kit-009")
    expect_identical(vapply(cases, `[[`, "", "id"), c("control-1", "candidate-1"))
    expect_identical(vapply(cases, `[[`, "", "variant"), c("grouped", "partitioned"))
    expect_identical(vapply(cases, `[[`, "", "library"), c("/lib/i", "/lib/i"))
    expect_identical(cases[[1]]$oracle, "journey_imputation")
    expect_identical(cases[[1]]$input$state, "replicate_filtered")
    expect_identical(cases[[2]]$parameters$proportion_missing_values, 0.5)
    expect_false(identical(cases[[1]]$binding, cases[[2]]$binding))
    # The KIT-006 pairs carry no variant: the coordinator keeps their production variants.
    expect_null(h$kitCases(b, "kit-006")[[1]]$variant)
    actual <- b$oracles$journey_imputation
    expect_identical(h$kitOracleStateMismatches(actual, b$oracles$journey_imputation), character())
    changed <- actual
    changed$components$peptide_data$sha256 <- "x"
    changed$summary_sha256 <- "other"
    changed$products <- list()
    expect_identical(h$kitOracleStateMismatches(changed, b$oracles$journey_imputation), c("components/peptide_data", "summary"))
    missing <- actual
    missing$support <- list()
    expect_identical(h$kitOracleStateMismatches(missing, b$oracles$journey_imputation), "support/names")
})

test_that("a journey population's input component is read from its snapshot head and differs from the state descriptor (KIT-009 attempt 1)", {
    h <- kitHarness()
    api <- function(name) get(name, envir = asNamespace("MultiScholaR"))
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "kit009-input", "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$saveState("replicate_filtered", object, list(), "fixture")
    state_id <- getDIAStateReference(manager)$entity_id
    closeDIAStateAdapter(manager)
    store_dir <- file.path(root, ".workflow-entities")
    component <- h$kitPopulationInputComponent(store_dir, "replicate_filtered", api)
    expect_true(is.character(component) && nzchar(component))
    expect_false(identical(component, state_id))
    expect_error(h$kitPopulationInputComponent(store_dir, "imputed", api))
    # A protein state's assay is its wide table, named by the class declaration (KIT-010).
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    root2 <- withr::local_tempdir()
    manager2 <- newDIAStateAdapter(newWorkflowEntityStore(root2, "kit010-input", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager2$.context$reference_states <- TRUE
    manager2$saveState("duplicates_removed", protein, list(), "fixture")
    expected <- readDIAStateMetadata(manager2)$data$components$protein_quant_table$entity_id
    closeDIAStateAdapter(manager2)
    expect_identical(h$kitPopulationInputComponent(file.path(root2, ".workflow-entities"), "duplicates_removed", api), expected)
})

test_that("the adapters vehicle names a member's input component through the class declaration, not a fixed table (KIT-010 audit PH-01)", {
    root <- if (exists("workflow_test_root")) workflow_test_root else here::here()
    vehicle <- readLines(file.path(root, "tests/testthat/test-workflow-execution-adapters.R"))
    # The KIT member path: the lines around the input-state selection (the frozen PERF-106 member above it keeps
    # its own peptide check).
    at <- grep("request\\$input_state", vehicle)[[1L]]
    timed <- vehicle[max(1L, at - 6L):min(length(vehicle), at + 8L)]
    expect_false(any(grepl("components\\$peptide_data\\$entity_id", timed)))
    expect_true(any(grepl("workflowClassDeclaration", timed)))
    # The worker entry point dispatches every compute role the registry validator admits.
    worker <- readLines(file.path(root, "R/utils_workflow_execution.R"))
    switch_line <- worker[grep("peptide = , protein = computeWorkflowPeptideAction", worker)]
    expect_length(switch_line, 1L)
})

test_that("a ticket's gui-mode member runs the screen journey with no record on its build and is decided on completion, dispatch and reproduction (KIT-008)", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$root <- withr::local_tempdir()
    cases <- h$kitCases(b, "kit-gui-008")
    expect_length(cases, 1L)
    case <- cases[[1L]]
    expect_identical(case$kind, "journey")
    expect_identical(case$role, "candidate")
    expect_identical(case$library, "/lib/k")
    expect_identical(case$corpus, "screen")
    expect_false(case$record)
    expect_identical(case$selector, "kit:gui-008:candidate:1")
    expect_identical(h$kitMemberSeconds("kit-gui-008", b), 600L)
    # The verifier must find the stage's one case by the id the coordinator gave it (KIT-008 correction 1).
    expect_identical(h$kitSelectorStage("kit:gui-008:candidate:1"), list(stage = "kit-gui-008", id = case$id))
    expect_error(h$kitCases(b, "kit-gui-009"), "no bound gui member")
    # Member directory with the driver's populations, two completed worker requests and the reference member.
    member <- file.path(b$root, "kit-gui-008", "member")
    reference <- file.path(b$root, "kit-lifecycle-complete", "member")
    dir.create(member, recursive = TRUE); dir.create(reference, recursive = TRUE)
    populations <- list(rollup = list(fingerprint = list(components = list(peptide_data = list(sha256 = "a")),
        support = list(), products = list(`dia-rollup-progress:1` = list(sha256 = "p")), summary_sha256 = "s", rows = 10L)))
    saveRDS(list(status = "complete", populations = populations), file.path(reference, "journey-result.rds"))
    saveRDS(list(status = "complete", populations = populations), file.path(member, "journey-result.rds"))
    for (operation in c("precursor_rollup", "minimum_peptides_per_sample")) {
        wfx <- file.path(member, "project", ".workflow-entities", ".workflow-export", paste0("wfx-", operation))
        dir.create(wfx, recursive = TRUE)
        saveRDS(list(schema = 1L, operation = operation), file.path(wfx, "request.rds"))
        writeLines('{"status":"complete"}', file.path(wfx, "journal.json"))
    }
    receipt <- kitSyntheticReceipt(cases, h)
    decision <- h$kitGuiMemberDecision(receipt, cases, b, member_dir = member)
    expect_identical(decision$status, "gui_member_observed")
    expect_identical(decision$dispatched$precursor_rollup$complete, 1L)
    expect_identical(decision$dispatched$minimum_peptides_per_sample$complete, 1L)
    expect_identical(decision$reproduction$fingerprint_mismatches, 0L)
    expect_false(decision$promotion_authorized)
    # Without the sample-filter dispatch the member is rejected; a changed state is a reference mismatch.
    unlink(file.path(member, "project", ".workflow-entities", ".workflow-export", "wfx-minimum_peptides_per_sample"), recursive = TRUE)
    rejected <- h$kitGuiMemberDecision(receipt, cases, b, member_dir = member)
    expect_identical(rejected$status, "rejected")
    expect_identical(rejected$failing_rule, "worker_dispatch_missing")
    wfx <- file.path(member, "project", ".workflow-entities", ".workflow-export", "wfx-minimum_peptides_per_sample")
    dir.create(wfx, recursive = TRUE)
    saveRDS(list(schema = 1L, operation = "minimum_peptides_per_sample"), file.path(wfx, "request.rds"))
    writeLines('{"status":"complete"}', file.path(wfx, "journal.json"))
    changed <- populations
    changed$rollup$fingerprint$components$peptide_data$sha256 <- "b"
    saveRDS(list(status = "complete", populations = changed), file.path(member, "journey-result.rds"))
    mismatch <- h$kitGuiMemberDecision(receipt, cases, b, member_dir = member)
    expect_identical(mismatch$status, "rejected")
    expect_identical(mismatch$failing_rule, "reference_mismatch")
    expect_identical(mismatch$reproduction$mismatches, "rollup/components/peptide_data")
    # An incomplete journey or a rejected verifier is invalid, never rejected.
    partial <- kitSyntheticReceipt(cases, h, function(r) { r$cases$member$operation$populations_complete <- FALSE; r })
    expect_identical(h$kitGuiMemberDecision(partial, cases, b, member_dir = member)$status, "invalid")
    unverified <- kitSyntheticReceipt(cases, h, function(r) { r$cases$member$verifier_completion$accepted <- FALSE; r })
    expect_identical(unname(h$kitGuiMemberDecision(unverified, cases, b, member_dir = member)$failing_rule), "verifier_accepted")
})

test_that("the canonical form drops the recorded worker specification, the args-embedded support copies and the step description (KIT-008 correction 2)", {
    h <- kitHarness()
    native <- list(filterMinNumPeptidesPerSample = list(peptides_per_sample_cutoff = 3L,
        core_utilisation = list(schema = 1L, backend = "multidplyr", workers = 8L), support_table = list(entity_id = "a",
            payload = list(sha256 = "s")), removal_ledger = list(entity_id = "b", payload = list(sha256 = "r")),
        filter_summary = list(kept = 2L)))
    declaration <- list(filterMinNumPeptidesPerSample = list(peptides_per_sample_cutoff = 3L, core_utilisation = NA,
        filter_summary = list(kept = 2L)))
    expect_identical(h$kitCanonical(native), h$kitCanonical(declaration))
    changed <- declaration
    changed$filterMinNumPeptidesPerSample$filter_summary$kept <- 1L
    expect_false(identical(h$kitCanonical(native), h$kitCanonical(changed)))
})

test_that("a ticket may bind single members: one operation member per step on one build, exactness and completion only (KIT-014)", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$builds$kit_014 <- list(library = "/lib/q", build = strrep("q", 64), exact_sum = "sum-q")
    population <- function(step, state, n) list(step = step, state = state, store = "/store/journey", head_sha256 = strrep(n, 64),
        head_entity = paste0("dia-state-head", n), entity_id = strrep(n, 64), rows = 100L)
    b$populations$journey_design <- population("design", "raw_data_s4", "1")
    b$populations$journey_qvalue <- population("qvalue", "qvalue_filtered", "2")
    b$tickets[["014"]] <- list(kind = "single", build = "kit_014", members = list(
        list(operation = "peptide_qvalue_filter", population = "journey_design", oracle = "journey_qvalue", variant = "streamed",
            parameters = list(qvalue_threshold = 0.01, proteotypic_only = TRUE)),
        list(operation = "peptide_protein_evidence_filter", population = "journey_qvalue", oracle = "journey_protein_peptide",
            variant = "partitioned", parameters = list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 1))))
    cases <- h$kitCases(b, "kit-014")
    expect_identical(vapply(cases, `[[`, "", "id"), c("single-1", "single-2"))
    expect_identical(vapply(cases, `[[`, "", "role"), c("single", "single"))
    expect_false(any(vapply(cases, `[[`, logical(1), "control")))
    expect_identical(vapply(cases, `[[`, "", "kind"), c("operation", "operation"))
    expect_identical(vapply(cases, `[[`, "", "variant"), c("streamed", "partitioned"))
    expect_identical(vapply(cases, `[[`, "", "oracle"), c("journey_qvalue", "journey_protein_peptide"))
    expect_identical(vapply(cases, `[[`, "", "library"), c("/lib/q", "/lib/q"))
    expect_identical(cases[[1]]$input$state, "raw_data_s4")
    expect_identical(cases[[2]]$parameters$num_peptides_per_protein_thresh, 1)
    expect_identical(cases[[2]]$selector, "kit:014:single:2")
    expect_false(identical(cases[[1]]$binding, cases[[2]]$binding))
    expect_identical(h$kitSelectorStage("kit:014:single:2"), list(stage = "kit-014", id = "single-2"))
    expect_error(h$kitSelectorStage("kit:014:single:6"))
    expect_identical(h$kitMemberSeconds("kit-014", b), 600L)
    receipt <- kitSyntheticOperationReceipt(cases, h, edit = function(r) {
        for (id in names(r$cases)) r$cases[[id]]$operation$sampled_anon_bytes <- 3e9
        r$cases[["single-2"]]$operation$science_execution <- list(admission = list(batches = 12L))
        r
    })
    decision <- h$kitSingleDecision(receipt, cases, b)
    expect_identical(decision$status, "accepted")
    expect_identical(decision$rule, "exactness_and_completion")
    expect_null(decision$failing_rule)
    expect_null(decision$pairs)
    expect_true(all(vapply(decision$members, `[[`, logical(1), "exact")))
    expect_identical(decision$members[["single-1"]]$peak_bytes, 5e9)
    expect_identical(decision$members[["single-1"]]$peak_anon_bytes, 3e9)
    expect_identical(decision$members[["single-2"]]$batches, 12L)
    expect_identical(decision$members[["single-2"]]$oracle, "journey_protein_peptide")
    expect_false(decision$promotion_authorized)
    # A member whose verifier named an oracle mismatch rejects; a verifier refusal for a harness reason, or any
    # other rule failing, is invalid (KIT-014 audit PH-01).
    inexact <- receipt
    inexact$cases[["single-2"]]$verifier_completion <- list(accepted = FALSE, reason = "Error: oracle mismatch: components/peptide_data, summary")
    expect_identical(h$kitSingleDecision(inexact, cases, b)$status, "rejected")
    expect_identical(h$kitSingleDecision(inexact, cases, b)$failing_rule, "oracle_mismatch")
    expect_false(h$kitSingleDecision(inexact, cases, b)$members[["single-2"]]$exact)
    refused <- receipt
    refused$cases[["single-2"]]$verifier_completion <- list(accepted = FALSE, reason = "Error: identical(operation$binding, case$binding) is not TRUE")
    expect_identical(h$kitSingleDecision(refused, cases, b)$status, "invalid")
    expect_identical(h$kitSingleDecision(refused, cases, b)$failing_rule, "verifier_accepted")
    late <- receipt
    late$cases[["single-1"]]$operation$elapsed_seconds <- 601
    expect_identical(h$kitSingleDecision(late, cases, b)$status, "invalid")
    expect_identical(h$kitSingleDecision(late, cases, b)$failing_rule, "elapsed_within_ceiling")
    # A member's bound input file must still carry its bound hash (KIT-015 audit PE-03).
    fasta <- withr::local_tempfile(fileext = ".RDS")
    saveRDS(data.frame(database_id = "P1"), fasta)
    hashed <- b
    hashed$tickets[["014"]]$members[[2]]$parameters$fasta_table <- fasta
    hashed$tickets[["014"]]$members[[2]]$inputs <- list(fasta_table = list(path = fasta, sha256 = digest::digest(file = fasta, algo = "sha256")))
    expect_identical(h$kitCases(hashed, "kit-014")[[2]]$parameters$fasta_table, fasta)
    saveRDS(data.frame(database_id = "P2"), fasta)
    expect_error(h$kitCases(hashed, "kit-014"), "missing or changed")
    # A member without a bound variant or oracle refuses at expansion; more than five members too.
    bare <- b
    bare$tickets[["014"]]$members[[1]]$oracle <- NULL
    expect_error(h$kitCases(bare, "kit-014"), "variant or oracle")
    many <- b
    many$tickets[["014"]]$members <- rep(b$tickets[["014"]]$members, 3)
    expect_error(h$kitCases(many, "kit-014"))
})

test_that("a ticket may bind a single journey member that stops after a step: exactness and completion only (KIT-018)", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$builds$kit_018 <- list(library = "/lib/r", build = strrep("r", 64), exact_sum = "sum-r")
    b$tickets[["018"]] <- list(kind = "single_journey", build = "kit_018", corpus = "full", stop_after = "design")
    b <- kitSyntheticOracles(b, h, "DIA", "journey")
    cases <- h$kitCases(b, "kit-journey-018")
    expect_length(cases, 1L)
    case <- cases[[1L]]
    expect_identical(case$id, "candidate-1")
    expect_identical(h$kitSelectorStage(case$selector)$id, case$id)
    expect_identical(case$kind, "journey")
    expect_identical(case$role, "candidate")
    expect_false(case$record)
    expect_identical(case$stop_after, "design")
    expect_identical(case$corpus, "full")
    expect_identical(case$library, "/lib/r")
    expect_identical(case$selector, "kit:journey-018:candidate:1")
    expect_identical(h$kitSelectorStage(case$selector), list(stage = "kit-journey-018", id = "candidate-1"))
    # The stop step is part of the binding; an unknown step or a missing build refuses at expansion.
    later <- b
    later$tickets[["018"]]$stop_after <- "qvalue"
    expect_false(identical(h$kitCases(later, "kit-journey-018")[[1L]]$binding, case$binding))
    bad <- b
    bad$tickets[["018"]]$stop_after <- "nowhere"
    expect_error(h$kitCases(bad, "kit-journey-018"), "stop_after")
    bare <- b
    bare$tickets[["018"]]$build <- NULL
    expect_error(h$kitCases(bare, "kit-journey-018"), "build or stop_after")
    expect_identical(h$kitMemberSeconds("kit-journey-018", b), as.integer(b$journey$ceiling_seconds))
    receipt <- kitSyntheticReceipt(cases, h, edit = function(r) {
        op <- r$cases[["candidate-1"]]$operation
        op$status <- "complete"
        op$populations_complete <- TRUE
        op$phase_status <- list(import = "complete", design = "complete")
        op$phase_peak_anon_bytes <- list(import = 2e9, design = 2.5e9)
        op$phase_seconds <- list(import = 100, design = 50)
        op$sampled_anon_bytes <- 2.5e9
        op$resident_end_bytes <- 1.5e9
        r$cases[["candidate-1"]]$operation <- op
        r$cases[["candidate-1"]]$verifier_completion <- list(accepted = TRUE)
        r
    })
    decision <- h$kitSingleJourneyDecision(receipt, cases, b)
    expect_identical(decision$status, "accepted")
    expect_identical(decision$rule, "exactness_and_completion")
    expect_null(decision$failing_rule)
    expect_identical(decision$steps, c("import", "design"))
    expect_identical(decision$stop_after, "design")
    expect_true(decision$member$exact)
    expect_identical(decision$member$peak_anon_bytes, 2.5e9)
    expect_identical(names(decision$member$phase_status), c("import", "design"))
    expect_false(decision$promotion_authorized)
    # A step after the stop is not expected; a missing step before it is a phase failure; an oracle mismatch rejects;
    # a harness refusal by the verifier is invalid.
    short <- receipt
    short$cases[["candidate-1"]]$operation$phase_status <- list(import = "complete")
    expect_identical(h$kitSingleJourneyDecision(short, cases, b)$status, "invalid")
    expect_identical(h$kitSingleJourneyDecision(short, cases, b)$failing_rule, "phase_timing")
    inexact <- receipt
    inexact$cases[["candidate-1"]]$verifier_completion <- list(accepted = FALSE, reason = "Error: oracle mismatch: design/components/peptide_data")
    expect_identical(h$kitSingleJourneyDecision(inexact, cases, b)$status, "rejected")
    expect_identical(h$kitSingleJourneyDecision(inexact, cases, b)$failing_rule, "oracle_mismatch")
    expect_false(h$kitSingleJourneyDecision(inexact, cases, b)$member$exact)
    refused <- receipt
    refused$cases[["candidate-1"]]$verifier_completion <- list(accepted = FALSE, reason = "Error: identical(operation$binding, case$binding) is not TRUE")
    expect_identical(h$kitSingleJourneyDecision(refused, cases, b)$status, "invalid")
    expect_identical(h$kitSingleJourneyDecision(refused, cases, b)$failing_rule, "verifier_accepted")
    incomplete <- receipt
    incomplete$cases[["candidate-1"]]$operation$populations_complete <- FALSE
    expect_identical(h$kitSingleJourneyDecision(incomplete, cases, b)$failing_rule, "populations_incomplete")
})

test_that("a product the oracle does not carry is reported as additional, a missing oracle product mismatches (KIT-017)", {
    h <- kitHarness()
    entry <- function(name) list(sha256 = strrep(substr(name, 1, 1), 64))
    fp <- function(products) list(fingerprint = list(class = "ProteinQuantitativeData", rows = 3L,
        components = list(protein_quant_table = entry("a")), support = list(),
        products = stats::setNames(lapply(products, entry), products)))
    oracle <- list(da = fp(c("da_fit:0ac83c189c9d", "qc_progress")))
    superset <- list(da = fp(c("da_fit:0ac83c189c9d", "qc_progress", "rendered_asset:96b0df0ef478")))
    result <- h$kitReproduceMismatches(oracle, superset)
    expect_length(result, 0L)
    expect_identical(attr(result, "additional_products"), "da/rendered_asset:96b0df0ef478")
    subset <- list(da = fp("qc_progress"))
    expect_identical(as.character(h$kitReproduceMismatches(oracle, subset)), "da/products/names")
    expect_length(attr(h$kitReproduceMismatches(oracle, subset), "additional_products"), 0L)
})

test_that("KIT-027: fasta_sha256 is not dropped by name and the verifier flags a state whose FASTA differs from its input", {
    h <- kitHarness()
    nested <- list(a = list(fasta_sha256 = "x"))
    expect_identical(h$kitCanonical(nested)$a$fasta_sha256, "x")
    fingerprint <- list(class = "ProteinQuantitativeData", rows = 2L, components = list(), support = list(), products = list(),
        summary_sha256 = strrep("e", 64))
    reference <- list(accession = list(fingerprint = fingerprint))
    matching <- list(accession = list(fingerprint = c(fingerprint, list(fasta_sha256 = strrep("f", 64), fasta_matches_input = TRUE))))
    differing <- list(accession = list(fingerprint = c(fingerprint, list(fasta_sha256 = strrep("0", 64), fasta_matches_input = FALSE))))
    expect_length(h$kitReproduceMismatches(reference, matching), 0L)
    expect_identical(as.vector(h$kitReproduceMismatches(reference, differing)), "accession/fasta_sha256")
})

test_that("KIT-027 audit PH-05: one FASTA comparison serves the journey and single-member verifiers", {
    h <- kitHarness()
    snapshot <- list(file_dependencies = list(fasta = list(sha256 = strrep("f", 64))))
    expect_null(h$kitFastaMatchesInput(list(class = "x"), snapshot))
    expect_true(h$kitFastaMatchesInput(list(fasta_sha256 = strrep("f", 64)), snapshot))
    expect_false(h$kitFastaMatchesInput(list(fasta_sha256 = strrep("0", 64)), snapshot))
    expect_true(h$kitFastaMatchesInput(list(fasta_sha256 = NA_character_), list()))
    source_text <- paste(readLines(file.path(if (exists("workflow_test_root")) workflow_test_root else here::here(),
        "tools/profiling/workflow_kit.R")), collapse = "\n")
    # Both verifiers call it: the store fingerprint (journey members) and the single-member oracle check.
    expect_identical(lengths(regmatches(source_text, gregexpr("kitFastaMatchesInput\\(", source_text))), 2L)
})

test_that("KIT-029: the analysis record bundle leaves the project file fingerprint unchanged", {
    h <- kitHarness()
    project <- withr::local_tempdir()
    dir.create(file.path(project, "results", "Study_report"), recursive = TRUE)
    writeLines("a\tb", file.path(project, "results", "table.tsv"))
    without <- h$kitFingerprintFiles(project)
    for (name in c("analysis_record.json", "analysis_record.docx", "analysis_record.xlsx", "analysis_record.html"))
        writeLines(name, file.path(project, "results", "Study_report", name))
    with <- h$kitFingerprintFiles(project)
    expect_identical(with[c("path", "bytes", "sha256", "method")], without[c("path", "bytes", "sha256", "method")])
})

test_that("KIT-054: the real-data keys bind a whole-journey single member with oracles on the new lane and without oracles on the legacy lane", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$builds$kit_054 <- list(library = "/lib/r", build = strrep("r", 64), exact_sum = "sum-r")
    b$corpora$kc_lfq <- b$corpora$full
    b$corpora$kc_lfq$report_file <- "proteinGroups.txt"; b$corpora$kc_lfq$format <- "maxquant"; b$corpora$kc_lfq$route <- "protein"
    b$corpora$kc_lfq$synthetic <- FALSE; b$corpora$kc_lfq$sensitive <- TRUE
    b$corpora$kc_lfq_dotted <- b$corpora$kc_lfq
    b$harness$step_map_protein_sha256 <- strrep("m", 64)
    b$tickets[["054"]] <- list(kind = "single_journey", build = "kit_054", corpus = "kc_lfq", route = "protein", whole_journey = TRUE,
        oracles = "kc_lfq", ceiling_seconds = 14400L)
    # KIT-092 (protocol v1.7): the legacy lane declares itself in memory.
    b$tickets[["054l"]] <- list(kind = "single_journey", build = "core_018", corpus = "kc_lfq_dotted", route = "protein", whole_journey = TRUE,
        ceiling_seconds = 14400L, in_memory = TRUE)
    b$tickets[["054m"]] <- list(kind = "single_journey", build = "core_018", corpus = "full", route = "protein", whole_journey = TRUE,
        ceiling_seconds = 14400L, in_memory = TRUE)
    b <- kitSyntheticOracles(b, h, "protein", "kc_lfq")
    new <- h$kitCases(b, "kit-journey-054")[[1L]]
    expect_identical(new$selector, "kit:journey-054:candidate:1")
    expect_identical(new$oracles, "kc_lfq")
    expect_identical(new$route, "protein")
    expect_false(new$record)
    legacy <- h$kitCases(b, "kit-journey-054l")[[1L]]
    expect_identical(legacy$selector, "kit:journey-054l:candidate:1")
    expect_null(legacy$oracles)
    expect_null(legacy$stop_after)
    expect_identical(legacy$corpus, "kc_lfq_dotted")
    expect_false(identical(legacy$binding, new$binding))
    expect_identical(h$kitCases(b, "kit-journey-054m")[[1L]]$corpus, "full")
    for (stage in c("kit-journey-054", "kit-journey-054l", "kit-journey-054m")) expect_identical(h$kitMemberSeconds(stage, b), 14400L)
    expect_identical(h$kitSelectorStage("kit:journey-054m:candidate:1"), list(stage = "kit-journey-054m", id = "candidate-1"))
    # The legacy lane's decision: complete with measures, no populations (no store), exactness not judged.
    protein <- h$kitStepIds("protein")
    legacy_cases <- h$kitCases(b, "kit-journey-054l")
    receipt <- kitSyntheticReceipt(legacy_cases, h, edit = function(r) {
        op <- r$cases[["candidate-1"]]$operation
        op$status <- "complete"; op$populations_complete <- FALSE
        op$phase_status <- stats::setNames(as.list(rep("complete", length(protein))), protein)
        op$phase_peak_anon_bytes <- stats::setNames(as.list(rep(2e9, length(protein))), protein)
        op$phase_seconds <- stats::setNames(as.list(rep(10, length(protein))), protein)
        op$sampled_anon_bytes <- 2e9; op$resident_end_bytes <- 1.5e9
        r$cases[["candidate-1"]]$operation <- op
        r$cases[["candidate-1"]]$verifier_completion <- list(accepted = TRUE)
        r
    })
    decision <- h$kitSingleJourneyDecision(receipt, legacy_cases, b)
    expect_identical(decision$status, "accepted")
    expect_true(decision$member$in_memory)
    expect_true(is.na(decision$member$exact))
    # The new lane still requires its populations.
    new_cases <- h$kitCases(b, "kit-journey-054")
    receipt2 <- kitSyntheticReceipt(new_cases, h, edit = function(r) {
        op <- r$cases[["candidate-1"]]$operation
        op$status <- "complete"; op$populations_complete <- FALSE
        op$phase_status <- stats::setNames(as.list(rep("complete", length(protein))), protein)
        op$phase_peak_anon_bytes <- stats::setNames(as.list(rep(2e9, length(protein))), protein)
        op$phase_seconds <- stats::setNames(as.list(rep(10, length(protein))), protein)
        op$sampled_anon_bytes <- 2e9; op$resident_end_bytes <- 1.5e9
        r$cases[["candidate-1"]]$operation <- op
        r$cases[["candidate-1"]]$verifier_completion <- list(accepted = TRUE)
        r
    })
    expect_identical(h$kitSingleJourneyDecision(receipt2, new_cases, b)$failing_rule, "populations_incomplete")
})

test_that("KIT-053: the DIA_LIMPA route's step list, its single journey member with the route's step map, oracles and normalisation settings", {
    h <- kitHarness()
    limpa <- h$kitStepIds("DIA_LIMPA")
    expect_identical(limpa, c("import", "design", "qvalue", "rollup", "peptide_log2", "peptide_accession", "dpc_quant",
        "protein_normalisation_limpa", "peptide_count", "correlation", "ruv_covariates", "session_export", "da", "report"))
    steps <- h$kitJourneySteps("DIA_LIMPA")
    expect_identical(steps[[7]]$state, "protein_dpc_quant")
    expect_identical(steps[[11]]$state, "ruv_covariates")
    expect_identical(steps[[11]]$action, "apply_ruv_covariates")
    expect_identical(steps[[13]]$state, "ruv_covariates")
    # KIT-056: the correction mode's states on the route; the DIA route ignores the argument.
    corrected <- h$kitJourneySteps("DIA_LIMPA", "correction")
    expect_identical(corrected[[11]]$state, "ruv_corrected")
    expect_identical(corrected[[13]]$state, "ruv_corrected")
    expect_identical(corrected[[11]]$action, "apply_ruv_covariates")
    expect_identical(vapply(corrected, `[[`, "", "id"), limpa)
    expect_identical(h$kitJourneySteps("DIA", "correction"), h$kitJourneySteps("DIA"))
    expect_identical(h$kitStepIds("DIA"), h$kitStepIds())
    b <- kitSyntheticBindings()
    b$builds$kit_053 <- list(library = "/lib/l", build = strrep("l", 64), exact_sum = "sum-l")
    b$harness$step_map_limpa_sha256 <- strrep("q", 64)
    b$tickets[["053"]] <- list(kind = "single_journey", build = "kit_053", corpus = "full", route = "DIA_LIMPA", whole_journey = TRUE,
        oracles = "limpa_full", ceiling_seconds = 14400L, normalisation = list(norm_method = "none", ruv_mode = "covariates", correlation_threshold = 0.5))
    b <- kitSyntheticOracles(kitSyntheticOracles(b, h, "DIA_LIMPA", "limpa_full"), h, "DIA_LIMPA", "limpa_corr_screen")
    cases <- h$kitCases(b, "kit-journey-053")
    expect_length(cases, 1L)
    case <- cases[[1L]]
    expect_identical(case$route, "DIA_LIMPA")
    expect_identical(case$oracles, "limpa_full")
    expect_identical(case$normalisation$ruv_mode, "covariates")
    expect_identical(case$selector, "kit:journey-053:candidate:1")
    expect_identical(h$kitSelectorStage("kit:journey-053:candidate:1"), list(stage = "kit-journey-053", id = "candidate-1"))
    expect_identical(h$kitMemberSeconds("kit-journey-053", b), 14400L)
    # The route's step map and its normalisation settings are part of the binding; without the map the case refuses.
    other <- b; other$harness$step_map_limpa_sha256 <- strrep("r", 64)
    expect_false(identical(h$kitCases(other, "kit-journey-053")[[1L]]$binding, case$binding))
    other <- b; other$tickets[["053"]]$normalisation$ruv_mode <- "skip"
    expect_false(identical(h$kitCases(other, "kit-journey-053")[[1L]]$binding, case$binding))
    none <- b; none$harness$step_map_limpa_sha256 <- NULL
    expect_error(h$kitCases(none, "kit-journey-053"), "limpa step map")
    # KIT-056: the correction-mode lifecycle member binds the mode and its own oracle prefix.
    b$builds$kit_056 <- list(library = "/lib/c", build = strrep("c", 64), exact_sum = "sum-c")
    b$tickets[["056"]] <- list(kind = "single_journey", build = "kit_056", corpus = "screen", route = "DIA_LIMPA", whole_journey = TRUE,
        oracles = "limpa_corr_screen", ceiling_seconds = 1800L, normalisation = list(norm_method = "none", ruv_mode = "correction", correlation_threshold = 0.5))
    corr_case <- h$kitCases(b, "kit-journey-056")[[1L]]
    expect_identical(corr_case$normalisation$ruv_mode, "correction")
    expect_identical(corr_case$oracles, "limpa_corr_screen")
    expect_identical(h$kitMemberSeconds("kit-journey-056", b), 1800L)
    # The oracles are looked up under the route prefix for the route's steps; the decision takes the whole step list.
    b$oracles$limpa_full_dpc_quant <- list(class = "ProteinQuantitativeData", components = list(protein_quant_table = "x"))
    oracles <- h$kitOracleJourneyPopulations(b, "DIA_LIMPA", "limpa_full")
    expect_identical(names(oracles), limpa)
    expect_identical(oracles$dpc_quant$fingerprint$class, "ProteinQuantitativeData")
    receipt <- kitSyntheticReceipt(cases, h, edit = function(r) {
        op <- r$cases[["candidate-1"]]$operation
        op$status <- "complete"; op$populations_complete <- TRUE
        op$phase_status <- stats::setNames(as.list(rep("complete", length(limpa))), limpa)
        op$phase_peak_anon_bytes <- stats::setNames(as.list(rep(2e9, length(limpa))), limpa)
        op$phase_seconds <- stats::setNames(as.list(rep(10, length(limpa))), limpa)
        op$sampled_anon_bytes <- 2e9; op$resident_end_bytes <- 1.5e9
        r$cases[["candidate-1"]]$operation <- op
        r$cases[["candidate-1"]]$verifier_completion <- list(accepted = TRUE)
        r
    })
    decision <- h$kitSingleJourneyDecision(receipt, cases, b)
    expect_identical(decision$status, "accepted")
    expect_identical(decision$steps, limpa)
    expect_true(decision$member$exact)
    expect_false(decision$member$in_memory)
})

test_that("KIT-057: the limpa route's RUV covariates product is an operation identity product across builds", {
    h <- kitHarness()
    pattern <- h$kitOperationIdentityProducts
    expect_true(all(grepl(pattern, c("dia-ruv-covariates:006bd48ece37", "ruv_covariates_science", "dia-pearson-selection:493ebf534ed2", "dia-normalization-ruv:0123456789ab"))))
    expect_false(any(grepl(pattern, c("qvalue_progress", "qc_progress", "dia_normalization_pca:479b0c17354b", "ruv_covariates_science_x"))))
    fp <- function(sha) list(components = list(), support = list(), products = list(ruv_covariates_science = list(kind = "scientific_product", sha256 = sha),
        qc_progress = list(kind = "value", sha256 = "same")), summary_sha256 = "s", rows = 1L)
    a <- list(ruv = list(fingerprint = fp("a"))); b <- list(ruv = list(fingerprint = fp("b")))
    expect_length(h$kitReproduceMismatches(a, b), 1L)
    expect_length(h$kitReproduceMismatches(a, b, ignore_products = pattern), 0L)
})

test_that("KIT-057: a journey pair binds its route, oracles and normalisation settings for both members", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$builds$kit_053 <- list(library = "/lib/l", build = strrep("l", 64), exact_sum = "sum-l")
    b$builds$kit_057 <- list(library = "/lib/m", build = strrep("m", 64), exact_sum = "sum-m")
    b$harness$step_map_limpa_sha256 <- strrep("q", 64)
    b$tickets[["057"]] <- list(kind = "journey", control = "kit_053", candidate = "kit_057", record = FALSE, corpus = "full",
        route = "DIA_LIMPA", oracles = "limpa_full", normalisation = list(norm_method = "cyclicloess", ruv_mode = "covariates", correlation_threshold = 0.5))
    cases <- h$kitCases(b, "kit-journey-057")
    expect_length(cases, 2L)
    expect_identical(vapply(cases, `[[`, character(1), "selector"), c("kit:journey-057:control:1", "kit:journey-057:candidate:1"))
    for (case in cases) {
        expect_identical(case$route, "DIA_LIMPA")
        expect_identical(case$oracles, "limpa_full")
        expect_identical(case$normalisation$norm_method, "cyclicloess")
        expect_false(case$record)
    }
    expect_identical(cases[[1L]]$library, "/lib/l")
    expect_identical(cases[[2L]]$library, "/lib/m")
    # The binding carries the route's settings; a pair without them (the DIA pairs) is unchanged in shape.
    other <- b; other$tickets[["057"]]$normalisation$norm_method <- "none"
    expect_false(identical(h$kitCases(other, "kit-journey-057")[[2L]]$binding, cases[[2L]]$binding))
    none <- b; none$harness$step_map_limpa_sha256 <- NULL
    expect_error(h$kitCases(none, "kit-journey-057"), "limpa step map")
    # Oracles by role: the candidate's fingerprints recorded on its own build (a different shipped engine).
    byrole <- b; byrole$tickets[["057"]]$oracles <- list(control = "limpa_full", candidate = "limpa57_full")
    pair <- h$kitCases(byrole, "kit-journey-057")
    expect_identical(pair[[1L]]$oracles, "limpa_full")
    expect_identical(pair[[2L]]$oracles, "limpa57_full")
    expect_false(identical(pair[[2L]]$binding, cases[[2L]]$binding))
    dia <- b; dia$tickets[["057"]] <- list(kind = "journey", control = "kit_053", candidate = "kit_057", record = FALSE)
    plain <- h$kitCases(dia, "kit-journey-057")
    expect_identical(plain[[1L]]$route, "DIA")
    expect_null(plain[[1L]]$oracles)
})

# KIT-081: the stage-3b instruments (protocol v1.6 section 8).
test_that("KIT-081: a suffixed operation key binds its ceiling and selectors, and unsuffixed keys keep theirs", {
    h <- new.env(parent = globalenv()); sys.source(testthat::test_path("../../tools/profiling/workflow_kit.R"), h)
    b <- list(tickets = list(`081a` = list(ceiling_seconds = 1800L), `081b` = list(ceiling_seconds = 2400L), `081c` = list()))
    expect_identical(h$kitMemberSeconds("kit-081a", b), 1800L)
    expect_error(h$kitMemberSeconds("kit-081b", b), "1,800")
    expect_error(h$kitMemberSeconds("kit-081c", b), "1,800")
    expect_identical(h$kitMemberSeconds("kit-081", b), 600L)
    expect_identical(h$kitSelectorStage("kit:081a:control:2"), list(stage = "kit-081a", id = "control-2"))
    expect_identical(h$kitSelectorStage("kit:081:candidate:1"), list(stage = "kit-081", id = "candidate-1"))
    expect_error(h$kitSelectorStage("kit:081ab:control:1"))
})

test_that("KIT-081: a list assay is identified by the digest of its sorted element ids, whatever their order", {
    h <- new.env(parent = globalenv()); sys.source(testthat::test_path("../../tools/profiling/workflow_kit.R"), h)
    one <- list(entity_id = "abc")
    expect_identical(h$kitAssayComponentId(one), "abc")
    a <- list(pos = list(entity_id = "b"), neg = list(entity_id = "a"))
    b <- list(neg = list(entity_id = "a"), pos = list(entity_id = "b"))
    expect_identical(h$kitAssayComponentId(a), h$kitAssayComponentId(b))
    expect_identical(h$kitAssayComponentId(a), digest::digest(c("a", "b"), algo = "sha256"))
})

test_that("KIT-081: the phase measures carry the controller's peak when the trace has it, and only then", {
    h <- new.env(parent = globalenv()); sys.source(testthat::test_path("../../tools/profiling/workflow_kit.R"), h)
    phases <- data.frame(uptime = c(1, 2, 3, 4), step = c("a", "a", "b", "b"), event = c("start", "complete", "start", "complete"))
    trace <- data.frame(uptime = c(1, 1.5, 2, 3, 3.5, 4), current = 1:6, anon = 11:16)
    expect_null(h$kitPhaseMeasures(phases, trace)$phase_controller_peak_anon_bytes)
    trace$controller_anon <- c(5, 7, 6, NA, 2, 3)
    expect_identical(h$kitPhaseMeasures(phases, trace)$phase_controller_peak_anon_bytes, list(a = 7, b = 3))
})

#' A receipt and cases from a recorded pair's members, as `kitJourneyTicketDecision` reads them.
kit081Replay <- function(members) {
    h <- new.env(parent = globalenv()); sys.source(testthat::test_path("../../tools/profiling/workflow_kit.R"), h)
    cases <- list(); receipt <- list(cases = list())
    for (role in c("control", "candidate")) {
        member <- members[[paste0(role, "-1")]] %||% members[[role]]
        id <- paste0(role, "-1")
        steps <- names(member$phase_peak_anon_bytes)
        cases[[length(cases) + 1L]] <- list(id = id, role = role, selector = "kit:journey-081:candidate:1")
        receipt$cases[[id]] <- list(verifier_completion = list(accepted = TRUE), operation = list(status = member$status %||% "complete",
            max_sample_interval = 0.05, peak_bytes = member$peak_bytes, swap_bytes = 0, cleanup = list(slice_stopped = TRUE), sampled_anon_bytes = member$peak_anon_bytes,
            populations_complete = TRUE, elapsed_seconds = member$elapsed_seconds, resident_end_bytes = member$resident_end_bytes,
            phase_peak_anon_bytes = member$phase_peak_anon_bytes, phase_seconds = member$phase_seconds,
            phase_status = member$phase_status %||% stats::setNames(rep(list("complete"), length(steps)), steps)))
    }
    list(h = h, receipt = receipt, cases = cases)
}

test_that("KIT-081: the recorded pair decisions of KIT-017, KIT-036 and KIT-057 reproduce unchanged", {
    for (file in c("kit-017/kit-journey-017-stage.json", "kit-036/kit-journey-036-stage.json", "kit-057/kit-journey-057-stage.json")) {
        path <- testthat::test_path("../../docs/qa/workflow-execution/kit", file)
        skip_if_not(file.exists(path))
        recorded <- jsonlite::read_json(path)
        replay <- kit081Replay(recorded$members)
        bindings <- list(journey = list(ceiling_seconds = 14400L), tickets = list(`081` = list(rule = recorded$decision$rule, ceiling_seconds = 14400L)))
        decision <- replay$h$kitJourneyTicketDecision(replay$receipt, replay$cases, bindings)
        expect_identical(decision$status, recorded$decision$status, info = file)
        expect_identical(as.character(decision$failing_rule %||% character()), as.character(unlist(recorded$decision$failing_rule)), info = file)
        expect_equal(decision$whole_slice_anon_ratio, recorded$decision$whole_slice_anon_ratio, info = file)
        expect_equal(decision$step_ratios[names(recorded$decision$step_ratios)], recorded$decision$step_ratios, info = file)
    }
})

test_that("KIT-081: a worker pair decides its two claims separately, the worker steps reported on the whole slice and judged on the application", {
    steps <- c("import", "design", "intensity", "qc_finalize", "normalisation", "report")
    member <- function(whole, step_anon, controller) list(status = "complete", outcome = "complete", peak_bytes = whole, peak_anon_bytes = whole, elapsed_seconds = 10,
        resident_end_bytes = 1, phase_seconds = stats::setNames(as.list(rep(1, 6)), steps), phase_peak_anon_bytes = stats::setNames(as.list(step_anon), steps))
    decide <- function(candidate_whole = 1.86, candidate_steps = c(1, 1, 1.24, 1, 1.7, 1), candidate_controller = c(0.9, 0.9), control_controller = c(1, 2.5)) {
        replay <- kit081Replay(list(`control-1` = member(2.5, c(1, 1, 1, 1, 2.5, 1)), `candidate-1` = member(candidate_whole, candidate_steps)))
        receipt <- replay$receipt
        receipt$cases$`control-1`$operation$phase_controller_peak_anon_bytes <- list(intensity = control_controller[[1]], normalisation = control_controller[[2]])
        receipt$cases$`candidate-1`$operation$phase_controller_peak_anon_bytes <- list(intensity = candidate_controller[[1]], normalisation = candidate_controller[[2]])
        replay$h$kitJourneyTicketDecision(receipt, replay$cases, list(journey = list(ceiling_seconds = 14400L),
            tickets = list(`081` = list(ceiling_seconds = 14400L, worker_steps = c("intensity", "normalisation")))))
    }
    good <- decide()
    expect_identical(good$status, "accepted")
    expect_identical(good$rule_decisions$whole_slice_anon$status, "accepted")
    expect_identical(good$rule_decisions$app_process$status, "accepted")
    # The intensity step's 1.24 on the whole slice is reported, not guarded.
    expect_equal(good$rule_decisions$whole_slice_anon$step_ratios_reported$intensity, 1.24)
    expect_false("intensity" %in% good$rule_decisions$whole_slice_anon$steps_guarded)
    # A non-worker step is still guarded.
    worse <- decide(candidate_steps = c(1, 1, 1.24, 1.5, 1.7, 1))
    expect_identical(worse$rule_decisions$whole_slice_anon$failing_rule, "step_not_worse")
    expect_identical(worse$rule_decisions$app_process$status, "accepted")
    expect_identical(worse$status, "rejected")
    # The total's reduction is required.
    expect_identical(decide(candidate_whole = 2.45)$rule_decisions$whole_slice_anon$failing_rule, "required_reduction")
    # The application's: every worker step at most 0.95.
    app <- decide(candidate_controller = c(0.99, 0.9))
    expect_identical(app$rule_decisions$app_process$failing_rule, "app_process_reduction")
    expect_identical(app$rule_decisions$whole_slice_anon$status, "accepted")
    expect_identical(app$failing_rule, "app_process_reduction")
    # A worker step with no controller reading cannot be decided.
    unobserved <- decide(control_controller = c(NA_real_, 2.5))
    expect_identical(unobserved$status, "invalid")
    expect_identical(unobserved$failing_rule, "worker_steps_unobserved")
})

test_that("KIT-081: an operation member's case carries its population's route, which marks a list assay's recorded state", {
    h <- new.env(parent = globalenv()); sys.source(testthat::test_path("../../tools/profiling/workflow_kit.R"), h)
    b <- list(protocol_sha256 = "p", builds = list(kit_081 = list(library = "l", build = "b", exact_sum = "e")),
        populations = list(lipid81_design = list(store = "s", head_sha256 = "h", rows = 1L, entity_id = "state", state = "lipid_raw_data_s4", route = "lipidomics"),
            plain = list(store = "s", head_sha256 = "h", rows = 1L, entity_id = "table")),
        tickets = list(`081a` = list(control = "kit_081", candidate = "kit_081", ceiling_seconds = 1800L, pairs = list(
            list(operation = "lipid_intensity_filter", population = "lipid81_design", oracle = "o", parameters = list(a = 1L), variants = list(control = "whole", candidate = "partitioned")),
            list(operation = "precursor_rollup", population = "plain", oracle = "o", parameters = list(a = 1L), variants = list(control = "whole", candidate = "whole"))))))
    cases <- h$kitCases(b, "kit-081a")
    expect_identical(vapply(cases, `[[`, "", "selector"), c("kit:081a:control:1", "kit:081a:candidate:1", "kit:081a:control:2", "kit:081a:candidate:2"))
    expect_identical(cases[[1L]]$input$route, "lipidomics")
    expect_identical(cases[[2L]]$variant, "partitioned")
    expect_null(cases[[3L]]$input$route)
})

kitTrxCorpusTool <- function() {
    h <- new.env(parent = globalenv()); sys.source(testthat::test_path("../../tools/profiling/kit_omics_corpus.R"), h)
    h
}

kitTrxCountCells <- function(directory) {
    utils::read.delim(file.path(directory, "counts.tsv"), colClasses = "character", na.strings = character(), check.names = FALSE, quote = "")
}

test_that("KIT-040: the trx_screen corpus regenerates byte for byte with the notebook's groups, technical pairs, shuffled columns, edge genes and one-value fixtures", {
    h <- kitTrxCorpusTool()
    recipe <- h$kitOmicsRecipe("trx_screen", 3L)
    expect_identical(recipe$seed, 20261106L)
    expect_identical(recipe$features_per_assay, 6000L)
    a <- withr::local_tempdir(); b <- withr::local_tempdir()
    ma <- h$writeKitTrxCorpus(a, recipe); mb <- h$writeKitTrxCorpus(b, recipe)
    files <- c("counts.tsv", "design_matrix.tab", "contrast_strings.tab", "design.tsv", "truth.json")
    expect_setequal(names(ma$files), files)
    expect_identical(ma$samples, 12L)
    for (f in files) expect_identical(unname(tools::md5sum(file.path(a, f))), unname(tools::md5sum(file.path(b, f))), info = f)
    expect_identical(ma$files, mb$files)
    # The notebook's levels and file layout: Run, group, replicates; the contrast file's header line and one entry.
    dm <- utils::read.delim(file.path(a, "design_matrix.tab"), check.names = FALSE)
    expect_identical(names(dm), c("Run", "group", "replicates"))
    expect_setequal(unique(dm$group), c("RPMI", "Sera"))
    expect_identical(readLines(file.path(a, "contrast_strings.tab")), c("contrasts", "Sera_vs_RPMI=groupSera-groupRPMI"))
    design <- utils::read.delim(file.path(a, "design.tsv"))
    expect_identical(names(design), c("sample", "group", "batch", "tech_rep_group", "replicates"))
    expect_identical(design$sample, dm$Run)
    expect_true(all(table(design$tech_rep_group) == 2L))
    expect_identical(length(unique(design$tech_rep_group)), 6L)
    expect_true(all(tapply(design$group, design$tech_rep_group, function(g) length(unique(g))) == 1L))
    expect_setequal(unique(design$batch), c("B1", "B2", "B3"))
    # The count columns are a seeded permutation of the design order, recorded in truth.json.
    truth <- jsonlite::read_json(file.path(a, "truth.json"), simplifyVector = TRUE)
    counts <- kitTrxCountCells(a)
    runs <- names(counts)[-(1:2)]
    expect_identical(names(counts)[1:2], c("GeneID", "GeneName"))
    expect_identical(runs, truth$column_order)
    expect_identical(runs, design$sample[truth$column_permutation])
    expect_false(identical(runs, design$sample))
    # DE is checked by Run ID: the per-group means taken by the design's Run IDs recover the declared fold changes.
    x <- matrix(as.numeric(unlist(counts[, runs])), ncol = length(runs), dimnames = list(counts$GeneID, runs))
    cpm <- t(t(x) / colSums(x)) * 1e6
    sera <- design$sample[design$group == "Sera"]; rpmi <- design$sample[design$group == "RPMI"]
    de <- truth$differential
    observed <- log2(rowMeans(cpm[de$GeneID, sera]) + 1) - log2(rowMeans(cpm[de$GeneID, rpmi]) + 1)
    expressed <- rowMeans(cpm[de$GeneID, ]) > 50
    expect_gt(mean(sign(observed[expressed]) == sign(de$log2fc[expressed])), 0.9)
    positional <- log2(rowMeans(cpm[de$GeneID, runs[seq_along(sera)]]) + 1) - log2(rowMeans(cpm[de$GeneID, runs[-seq_along(sera)]]) + 1)
    expect_lt(abs(stats::cor(positional[expressed], de$log2fc[expressed])), abs(stats::cor(observed[expressed], de$log2fc[expressed])))
    expect_identical(nrow(de), 600L)
    expect_true(all(abs(de$log2fc) >= 1 & abs(de$log2fc) <= 2))
    # Edge content: 50 all-zero genes, one constant gene of 100, 1% each of NA, empty and duplicated GeneName.
    # The 50 declared all-zero genes; low-mean genes may also draw all zeros by chance, which truth.json counts.
    expect_identical(length(truth$zero_genes), 50L)
    expect_true(all(truth$zero_genes %in% rownames(x)[rowSums(x) == 0]))
    expect_identical(sum(rowSums(x) == 0), truth$observed_all_zero_genes)
    constant <- rownames(x)[apply(x, 1, function(r) all(r == 100))]
    expect_identical(constant, truth$constant_genes)
    expect_identical(sum(counts$GeneName == "NA"), 60L)
    expect_identical(sum(counts$GeneName == ""), 60L)
    named <- counts$GeneName[!counts$GeneName %in% c("", "NA")]
    expect_identical(sum(duplicated(named)), 60L)
    expect_true(all(x == round(x)) && all(x >= 0))
    # The fixtures: each differs from its parent in one value.
    fixtures <- h$writeKitTrxFixtures(a, withr::local_tempdir())
    expect_setequal(names(fixtures), c("counts_fractional", "counts_negative", "counts_na", "counts_whole_double"))
    parent <- as.matrix(counts)
    expected <- c(counts_fractional = "20.37", counts_negative = "-1", counts_na = "NA", counts_whole_double = "20.0")
    for (name in names(fixtures)) {
        changed <- as.matrix(kitTrxCountCells(fixtures[[name]]))
        expect_identical(dim(changed), dim(parent))
        expect_identical(sum(changed != parent), 1L, info = name)
        expect_identical(changed[changed != parent], unname(expected[[name]]), info = name)
        for (f in c("design_matrix.tab", "contrast_strings.tab", "design.tsv")) {
            expect_identical(unname(tools::md5sum(file.path(fixtures[[name]], f))), unname(tools::md5sum(file.path(a, f))))
        }
    }
})

test_that("KIT-040: the notebook reference driver's preflight refuses invalid counts, a missing input and a missing package by class", {
    h <- kitTrxCorpusTool()
    corpus <- withr::local_tempdir()
    h$writeKitTrxCorpus(corpus, h$kitOmicsRecipe("trx_screen", 3L))
    fixtures <- h$writeKitTrxFixtures(corpus, withr::local_tempdir())
    d <- new.env(parent = globalenv()); sys.source(testthat::test_path("../../tools/profiling/trx_notebook_reference.R"), d)
    # The package list is the notebook's own p_load list (and pacman), so the notebook never installs at render time.
    packages <- d$trxNotebookPackages(testthat::test_path("../../dev/distributable/plans/transcriptomics_workflow.rmd"))
    expect_true(all(c("pacman", "edgeR", "limma", "RUVSeq", "ruv", "Glimma", "here", "vroom") %in% packages))
    for (name in c("counts_fractional", "counts_negative", "counts_na")) {
        expect_error(d$trxPreflight(fixtures[[name]], packages = "stats"), class = "kit_trx_invalid_counts", info = name)
    }
    expect_no_error(d$trxPreflight(fixtures[["counts_whole_double"]], packages = "stats"))
    expect_no_error(d$trxPreflight(corpus, packages = "stats"))
    for (missing in c("design_matrix.tab", "contrast_strings.tab")) {
        partial <- withr::local_tempdir()
        file.copy(list.files(corpus, full.names = TRUE), partial)
        unlink(file.path(partial, missing))
        expect_error(d$trxPreflight(partial, packages = "stats"), class = "kit_trx_missing_input", info = missing)
    }
    expect_error(d$trxPreflight(corpus, packages = c("stats", "kitTrxNoSuchPackage")), class = "kit_trx_missing_package")
})

# KIT-092 (protocol v1.7, user decision 2026-11-06): content fingerprints for tables inside products; KIT098-F02.
test_that("KIT-092: a table nested in a product is hashed by its content, so the same content in different bytes matches", {
    skip_if_not_installed("fastduckdb")
    skip_if(inherits(try(workflowFastWriterInfo("fast_r"), silent = TRUE), "try-error"), "the direct table writer is not qualified here")
    h <- kitHarness()
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, basename(root), "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    identity <- function(tag) newSavedWorkflowScientificIdentity("input", "design", tag, "all", "kit092")
    table <- data.frame(id = letters[1:4], quantity = as.double(1:4))
    whole <- saveWorkflowTable(store, table, "whole", identity("whole"), "id")
    # The same rows through the other writer: the same content in different bytes, as the partitioned intensity filter gives.
    other <- saveWorkflowTable(store, table, "other", identity("other"), "id", write_fn = function(payload, metadata) {
        writer <- openWorkflowTableWriter(payload, metadata, NULL, backend = "fast_r")
        appendWorkflowTableWriter(writer, table, seq_len(nrow(table)))
        finishWorkflowTableWriter(writer)
    })
    changed_table <- table
    changed_table$quantity[[3L]] <- 30
    changed <- saveWorkflowTable(store, changed_table, "changed", identity("changed"), "id")
    expect_false(identical(whole$payload$sha256, other$payload$sha256))
    product <- function(reference) list(product_type = "lipid.intensity_filter_result", counts = c(206L, 189L), filtered = reference)
    canonical <- function(reference) h$kitCanonical(product(reference), store, api, memo = new.env(parent = emptyenv()))
    expect_identical(canonical(whole), canonical(other))
    expect_false(identical(canonical(whole), canonical(changed)))
    expect_identical(canonical(whole)$filtered$content, h$kitFingerprintTable(store, whole, api)$sha256)
    # Without a store (the args summary, an .rds publication file) a nested table keeps its payload hash, as before.
    expect_identical(h$kitCanonical(product(whole))$filtered$sha256, whole$payload$sha256)
})

test_that("KIT-092: the verifier rebuilds the products the driver attached to each step from the recomputed store", {
    h <- kitHarness()
    product <- function(type, sha) list(kind = "scientific_product", product_type = type, sha256 = sha)
    products <- list(p1 = product("lipid.intensity_filter_result", strrep("b", 64)), p2 = product(NULL, strrep("a", 64)),
        p3 = product("omics-figure", strrep("c", 64)))
    products$p2$kind <- "rendered_asset"
    populations <- list(import = list(state = NULL, product_ids = character()),
        intensity = list(state = "lipid_intensity_filtered", product_ids = c("p1", "p2")),
        duplicates = list(state = "lipid_duplicates_resolved", product_ids = c("p1", "p2", "p3", "gone")))
    attached <- h$kitStepProducts(populations, products)
    # First appearance, ordered by hash, named <type>:<sha12> as the driver names them; an id the head no longer holds is skipped.
    expect_identical(names(attached$intensity), c(paste0("rendered_asset:", strrep("a", 12)), paste0("lipid.intensity_filter_result:", strrep("b", 12))))
    expect_identical(names(attached$duplicates), paste0("omics-figure:", strrep("c", 12)))
    expect_length(attached$import, 0L)
    # The in-process recorder leaves rendered pages out.
    expect_identical(names(h$kitStepProducts(populations, products, rendered = FALSE)$intensity), paste0("lipid.intensity_filter_result:", strrep("b", 12)))
})

test_that("KIT-092 (KIT098-F02): a single journey with an unbound step oracle, or an undeclared in-memory member, is refused", {
    h <- kitHarness()
    b <- kitSyntheticBindings()
    b$builds$kit_092 <- list(library = "/lib/q", build = strrep("q", 64), exact_sum = "sum-q")
    b$corpora$lipid_screen <- b$corpora$screen
    b$harness <- list(step_map_lipid_sha256 = strrep("l", 64))
    b$tickets[["092a"]] <- list(kind = "single_journey", build = "kit_092", corpus = "lipid_screen", route = "lipidomics",
        whole_journey = TRUE, oracles = "lipid92", ceiling_seconds = 1800L)
    b <- kitSyntheticOracles(b, h, "lipidomics", "lipid92")
    case <- h$kitCases(b, "kit-journey-092a")[[1L]]
    expect_false(case$in_memory)
    # An unbound step oracle is refused at binding, never skipped at verification; an oracle without components is unbound.
    gap <- b; gap$oracles$lipid92_duplicates <- NULL
    expect_error(h$kitCases(gap, "kit-journey-092a"), "duplicates")
    empty <- b; empty$oracles$lipid92_intensity$components <- list()
    expect_error(h$kitCases(empty, "kit-journey-092a"), "intensity")
    # A member stopped after a step needs the oracles of the steps it runs only.
    early <- gap; early$tickets[["092a"]]$whole_journey <- NULL; early$tickets[["092a"]]$stop_after <- "intensity"
    expect_identical(h$kitCases(early, "kit-journey-092a")[[1L]]$stop_after, "intensity")
    # A member without oracles is refused unless it declares in_memory (DIA on the full corpus keeps the KIT-005 oracles).
    bare <- b; bare$tickets[["092a"]]$oracles <- NULL
    expect_error(h$kitCases(bare, "kit-journey-092a"), "in_memory")
    screen <- b; screen$tickets[["092d"]] <- list(kind = "single_journey", build = "kit_092", corpus = "screen", whole_journey = TRUE, ceiling_seconds = 1800L)
    expect_error(h$kitCases(screen, "kit-journey-092d"), "in_memory")
    declared <- bare; declared$tickets[["092a"]]$in_memory <- TRUE
    memory_case <- h$kitCases(declared, "kit-journey-092a")[[1L]]
    expect_true(memory_case$in_memory)
    # Decisions: in memory only where declared; an undeclared member needs its populations and its exactness is judged.
    ids <- h$kitStepIds("lipidomics")
    receipt_for <- function(case, populations_complete) list(cases = stats::setNames(list(list(passed = TRUE,
        verifier_completion = list(accepted = TRUE), operation = list(status = "complete", swap_bytes = 0, peak_bytes = 4e9,
            elapsed_seconds = 900, max_sample_interval = 0.06, cleanup = list(slice_stopped = TRUE), resident_end_bytes = 1e9,
            populations_complete = populations_complete, head_unchanged = NA,
            phase_status = as.list(stats::setNames(rep("complete", length(ids)), ids)),
            phase_peak_anon_bytes = as.list(stats::setNames(rep(2e9, length(ids)), ids))))), case$id))
    judged <- h$kitSingleJourneyDecision(receipt_for(case, FALSE), list(case), b)
    expect_identical(judged$failing_rule, "populations_incomplete")
    expect_false(judged$member$in_memory)
    exact <- h$kitSingleJourneyDecision(receipt_for(case, TRUE), list(case), b)
    expect_identical(exact$status, "accepted")
    expect_true(exact$member$exact)
    memory <- h$kitSingleJourneyDecision(receipt_for(memory_case, FALSE), list(memory_case), declared)
    expect_identical(memory$status, "accepted")
    expect_true(is.na(memory$member$exact))
    # An omics lifecycle member is in memory only when the binding declares it: undeclared, the populations and head rules apply.
    b$lifecycle_omics <- list(ceiling_seconds = 1800L, build = "kit_092",
        fault = list(name = "guard_after_qc_finalize", after_step = "qc_finalize", failed_step = "normalisation"),
        normalisation = list(ruv_mode = "automatic"))
    complete <- h$kitCases(b, "kit-lifecycle-lipid")
    expect_false(complete[[1L]]$in_memory)
    expect_identical(h$kitLifecycleDecision(receipt_for(complete[[1L]], FALSE), complete, b)$failing_rule, "populations_incomplete")
    fail <- h$kitCases(b, "kit-lifecycle-lipid-fail")
    failed <- receipt_for(fail[[1L]], TRUE)
    failed$cases$member$operation$status <- "journey_failed"
    failed$cases$member$operation$phase_status$normalisation <- "failed"
    failed$cases$member$operation$guard <- list(reason = "injected_guard_stop")
    expect_identical(h$kitLifecycleDecision(failed, fail, b)$failing_rule, "head_unchanged")
    failed$cases$member$operation$head_unchanged <- TRUE
    expect_identical(h$kitLifecycleDecision(failed, fail, b)$status, "lifecycle_failure_observed")
})

# KIT-098: the transcriptomics route (protocol v1.8, kit-098/plan.md). Each test failed before its change.
kitSyntheticTrxBindings <- function(h) {
    b <- kitSyntheticBindings()
    file <- function(name) list(path = paste0("/corpus/", name), bytes = 1, sha256 = strrep("a", 64))
    b$builds$kit_098 <- list(library = "/lib/t", build = strrep("t", 64), exact_sum = "sum-t")
    b$corpora$trx_screen <- list(profile = "trx_screen", samples = 24L, route = "transcriptomics", format = "trx_counts",
        vendor = "trx_counts", report_file = "counts.tsv", assay_files = list(counts = "counts.tsv"), design_file = "design.tsv",
        files = list(`counts.tsv` = file("counts.tsv"), `design_matrix.tab` = file("design_matrix.tab"),
            `contrast_strings.tab` = file("contrast_strings.tab"), `truth.json` = file("truth.json"), `design.tsv` = file("design.tsv")))
    b$harness <- list(step_map_trx_sha256 = strrep("x", 64))
    norm <- list(ruv_mode = "automatic", auto_percentage_min = 5, auto_percentage_max = 30, auto_percentage_step = 5, k_penalty_weight = 0.05)
    settings <- list(da = list(q_value = 0.05, treat_lfc = 1L), ruviii_m_source = "group")
    workers <- c("trx_expression_filter", "trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_da_contrast", "trx_ruviii_export")
    b$lifecycle_trx <- list(ceiling_seconds = 1800L, build = "kit_098", corpus = "trx_screen", oracles = "trx98a",
        normalisation = norm, driver_settings = settings,
        fault = list(name = "guard_after_qc_finalize", after_step = "qc_finalize", failed_step = "normalisation"))
    b$tickets[["098a"]] <- list(kind = "single_journey", build = "kit_098", corpus = "trx_screen", route = "transcriptomics",
        whole_journey = TRUE, oracles = "trx98a", ceiling_seconds = 1800L, normalisation = norm, driver_settings = settings,
        worker_steps = as.list(workers))
    b <- kitSyntheticOracles(b, h, "transcriptomics", "trx98a")
    b$oracles$trx98a_ruviii_export$products <- list(`trx.ruviii:0123456789ab` = list(kind = "value", sha256 = "r"))
    b
}

test_that("KIT-098: the transcriptomics step list follows KIT-097's GUI, and an unknown route is refused", {
    h <- kitHarness()
    expect_identical(h$kitStepIds("transcriptomics"), c("import", "design", "expression_filter", "qc_finalize", "normalisation",
        "da", "ruviii_export", "report", "session_export"))
    steps <- h$kitJourneySteps("transcriptomics")
    expect_identical(lapply(steps, `[[`, "state"), list(NULL, "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete",
        "trx_ruvs", "trx_ruvs", "trx_ruvs", NULL, NULL))
    # A route the harness does not know is refused, never driven or bound as DIA.
    b <- kitSyntheticTrxBindings(h)
    b$tickets[["098z"]] <- list(kind = "single_journey", build = "kit_098", corpus = "trx_screen", route = "proteomics",
        whole_journey = TRUE, in_memory = TRUE, ceiling_seconds = 1800L)
    expect_error(h$kitCases(b, "kit-journey-098z"), "route")
    unknown <- list(id = "member", route = "proteomics", in_memory = FALSE, expect_failure = FALSE)
    expect_error(h$kitLifecycleDecision(list(cases = list(member = list(operation = list()))), list(unknown), b), "route")
})

test_that("KIT-098: a trx case binds its step map, its oracles and its own settings; the lifecycle pair takes 1,800 s", {
    h <- kitHarness()
    b <- kitSyntheticTrxBindings(h)
    case <- h$kitCases(b, "kit-journey-098a")[[1L]]
    expect_identical(case$route, "transcriptomics")
    expect_identical(case$oracles, "trx98a")
    expect_false(case$in_memory)
    expect_identical(case$normalisation$auto_percentage_max, 30)
    expect_identical(case$driver_settings$ruviii_m_source, "group")
    other <- b; other$harness$step_map_trx_sha256 <- strrep("y", 64)
    expect_false(identical(h$kitCases(other, "kit-journey-098a")[[1L]]$binding, case$binding))
    none <- b; none$harness$step_map_trx_sha256 <- NULL
    expect_error(h$kitCases(none, "kit-journey-098a"), "trx step map")
    # The lifecycle pair: its own spec, oracles bound, never in memory, the failure member with the launch record.
    complete <- h$kitCases(b, "kit-lifecycle-trx")
    expect_identical(complete[[1L]]$selector, "kit:lifecycle:trx:complete")
    expect_identical(complete[[1L]]$oracles, "trx98a")
    expect_false(complete[[1L]]$in_memory)
    expect_false(complete[[1L]]$record)
    fail <- h$kitCases(b, "kit-lifecycle-trx-fail")
    expect_true(fail[[1L]]$record)
    expect_identical(fail[[1L]]$fault_spec$failed_step, "normalisation")
    expect_identical(h$kitSelectorStage("kit:lifecycle:trx:complete"), list(stage = "kit-lifecycle-trx", id = "member"))
    expect_identical(h$kitSelectorStage("kit:lifecycle:trx:fail"), list(stage = "kit-lifecycle-trx-fail", id = "member"))
    expect_identical(h$kitMemberSeconds("kit-lifecycle-trx", b), 1800L)
    expect_identical(h$kitMemberSeconds("kit-lifecycle-trx-fail", b), 1800L)
    long <- b; long$lifecycle_trx$ceiling_seconds <- 3600L
    expect_error(h$kitMemberSeconds("kit-lifecycle-trx", long), "1,800")
    memory <- b; memory$lifecycle_trx$in_memory <- TRUE
    expect_error(h$kitCases(memory, "kit-lifecycle-trx"), "in memory")
    ids <- h$kitStepIds("transcriptomics")
    receipt_for <- function(case, edit = NULL) {
        operation <- list(status = "complete", swap_bytes = 0, peak_bytes = 4e9, elapsed_seconds = 1700, max_sample_interval = 0.06,
            cleanup = list(slice_stopped = TRUE), resident_end_bytes = 1e9, populations_complete = TRUE, head_unchanged = NA,
            phase_status = as.list(stats::setNames(rep("complete", length(ids)), ids)),
            phase_peak_anon_bytes = as.list(stats::setNames(rep(2e9, length(ids)), ids)))
        r <- list(cases = stats::setNames(list(list(passed = TRUE, operation = operation, verifier_completion = list(accepted = TRUE))), case$id))
        if (is.function(edit)) r <- edit(r)
        r
    }
    # 1,700 s is inside the trx lifecycle ceiling (the DIA lifecycle's 300 s would make it invalid).
    decision <- h$kitLifecycleDecision(receipt_for(complete[[1L]]), complete, b)
    expect_identical(decision$status, "lifecycle_complete_observed")
    expect_false(decision$in_memory)
    failed <- receipt_for(fail[[1L]], function(r) {
        r$cases$member$operation$status <- "journey_failed"
        r$cases$member$operation$phase_status$normalisation <- "failed"
        r$cases$member$operation$guard <- list(reason = "injected_guard_stop")
        r
    })
    expect_identical(h$kitLifecycleDecision(failed, fail, b)$failing_rule, "head_unchanged")
    failed$cases$member$operation$head_unchanged <- TRUE
    expect_identical(h$kitLifecycleDecision(failed, fail, b)$status, "lifecycle_failure_observed")
})

test_that("KIT-098: a trx case without an oracle for a step, or without the RUVIII product, is refused (oracle_missing)", {
    h <- kitHarness()
    b <- kitSyntheticTrxBindings(h)
    gap <- b; gap$oracles$trx98a_normalisation <- NULL
    expect_error(h$kitCases(gap, "kit-journey-098a"), "oracle_missing.*normalisation")
    expect_error(h$kitCases(gap, "kit-lifecycle-trx"), "oracle_missing.*normalisation")
    product <- b; product$oracles$trx98a_ruviii_export$products <- list()
    expect_error(h$kitCases(product, "kit-journey-098a"), "oracle_missing.*RUVIII")
    expect_error(h$kitCases(product, "kit-lifecycle-trx-fail"), "oracle_missing.*RUVIII")
    bare <- b; bare$tickets[["098a"]]$oracles <- NULL; bare$tickets[["098a"]]$in_memory <- TRUE
    expect_error(h$kitCases(bare, "kit-journey-098a"), "oracles")
    unbound <- b; unbound$lifecycle_trx$oracles <- NULL
    expect_error(h$kitCases(unbound, "kit-lifecycle-trx"), "oracles")
})

test_that("KIT-098: a trx product outside the identity set is judged; only the trx *_science products are reported", {
    h <- kitHarness()
    expect_true(grepl(h$kitOperationIdentityProducts, "trx_ruvs_science"))
    expect_true(grepl(h$kitOperationIdentityProducts, "trx_expression_filter_science"))
    expect_false(grepl(h$kitOperationIdentityProducts, "trx.ruvs_result:0123456789ab"))
    expect_false(grepl(h$kitOperationIdentityProducts, "trx.ruviii:0123456789ab"))
    fingerprint <- function(science, result) list(fingerprint = list(class = "TranscriptomicsAssayData", rows = 4565L,
        components = list(count_data = list(sha256 = "c")), support = list(), summary_sha256 = "s",
        products = list(trx_ruvs_science = list(sha256 = science), `trx.ruvs_result:0123456789ab` = list(sha256 = result))))
    reference <- list(normalisation = fingerprint("a", "w"))
    expect_length(h$kitReproduceMismatches(reference, list(normalisation = fingerprint("b", "w")), h$kitOperationIdentityProducts), 0L)
    expect_identical(as.vector(h$kitReproduceMismatches(reference, list(normalisation = fingerprint("a", "v")), h$kitOperationIdentityProducts)),
        "normalisation/products/trx.ruvs_result:0123456789ab")
    # The state's history is judged by name against the oracle's recorded history.
    b <- kitSyntheticTrxBindings(h)
    b$oracles$trx98a_history <- list("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs")
    case <- h$kitCases(b, "kit-journey-098a")[[1L]]
    expect_length(h$kitHistoryMismatches(b, case, c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs")), 0L)
    expect_identical(h$kitHistoryMismatches(b, case, c("initial", "trx_raw_data_s4", "trx_qc_complete", "trx_ruvs")), "history")
    nohistory <- b; nohistory$oracles$trx98a_history <- NULL
    expect_identical(h$kitHistoryMismatches(nohistory, case, "initial"), "history")
})

test_that("KIT-098: a trx journey needs a completed worker journal for each listed worker step (worker_dispatch_missing)", {
    h <- kitHarness()
    b <- kitSyntheticTrxBindings(h)
    case <- h$kitCases(b, "kit-journey-098a")[[1L]]
    ids <- h$kitStepIds("transcriptomics")
    member <- withr::local_tempdir()
    journal <- function(operation, status) {
        dir <- file.path(member, "project", ".dia-native-x", ".workflow-entities", ".workflow-export", paste0("wfx-", operation))
        dir.create(dir, recursive = TRUE)
        saveRDS(list(operation = operation), file.path(dir, "request.rds"))
        jsonlite::write_json(list(status = status), file.path(dir, "journal.json"), auto_unbox = TRUE)
    }
    for (operation in c("trx_expression_filter", "trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_da_contrast")) journal(operation, "complete")
    journal("trx_ruviii_export", "worker_failed")
    receipt <- list(cases = list(`candidate-1` = list(passed = TRUE, verifier_completion = list(accepted = TRUE),
        operation = list(status = "complete", swap_bytes = 0, peak_bytes = 4e9, elapsed_seconds = 900, max_sample_interval = 0.06,
            cleanup = list(slice_stopped = TRUE), resident_end_bytes = 1e9, populations_complete = TRUE,
            phase_status = as.list(stats::setNames(rep("complete", length(ids)), ids)),
            phase_peak_anon_bytes = as.list(stats::setNames(rep(2e9, length(ids)), ids))))))
    missing <- h$kitSingleJourneyDecision(receipt, list(case), b, member_dir = member)
    expect_identical(missing$status, "invalid")
    expect_identical(missing$failing_rule, "worker_dispatch_missing")
    expect_identical(missing$dispatched$trx_ruviii_export, list(requests = 1L, complete = 0L))
    journal("trx_ruviii_export_again", "complete")
    file.rename(file.path(member, "project/.dia-native-x/.workflow-entities/.workflow-export/wfx-trx_ruviii_export_again"),
        file.path(member, "project/.dia-native-x/.workflow-entities/.workflow-export/wfx-2"))
    saveRDS(list(operation = "trx_ruviii_export"), file.path(member, "project/.dia-native-x/.workflow-entities/.workflow-export/wfx-2/request.rds"))
    accepted <- h$kitSingleJourneyDecision(receipt, list(case), b, member_dir = member)
    expect_identical(accepted$status, "accepted")
    expect_identical(accepted$worker_steps, c("trx_expression_filter", "trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_da_contrast", "trx_ruviii_export"))
    # Without the member's directory the listed steps cannot be shown, so the member is not accepted.
    expect_identical(h$kitSingleJourneyDecision(receipt, list(case), b)$failing_rule, "worker_dispatch_missing")
})

test_that("KIT-098: the trx driver runs its own module ids; the lipid-or-metab steps refuse any other route", {
    d <- new.env(parent = globalenv())
    sys.source(testthat::test_path("../../tools/profiling/workflow_kit.R"), d)
    sys.source(testthat::test_path("../../tools/profiling/workflow_kit_driver.R"), d)
    expect_error(d$runKitOmicsImportActions(NULL, NULL, tempdir(), list(assays = list(counts = "x")), "transcriptomics"), "lipidomics or metabolomics")
    expect_error(d$runKitOmicsStep(NULL, NULL, tempdir(), list(), list(route = "transcriptomics"), list(id = "intensity")), "lipidomics or metabolomics")
    seen <- new.env(parent = emptyenv()); seen$ids <- character(); seen$values <- list()
    h <- list(checkpointSetInput = function(driver, id, value) { seen$ids <- c(seen$ids, id); seen$values[[id]] <- value },
        checkpointClickControl = function(driver, id) seen$ids <- c(seen$ids, id),
        checkpointWaitForJs = function(...) TRUE, checkpointWaitBudget = function(...) 1000L, checkpointWaitForIdle = function(...) TRUE,
        checkpointFileSelection = function(path, roots) list(files = path),
        checkpointClickAndWaitText = function(driver, id, selector, text) { seen$ids <- c(seen$ids, id); seen$values[[paste0(id, "#wait")]] <- text })
    driver <- list(get_js = function(script, timeout = NULL) {
        id <- regmatches(script, regexpr("getElementById\\(\"[^\"]+\"\\)\\.click", script))
        if (length(id)) seen$ids <- c(seen$ids, sub("^getElementById\\(\"([^\"]+)\"\\)\\.click$", "\\1", id))
        NULL
    }, set_inputs = function(...) invisible(NULL))
    journey <- list(route = "transcriptomics", normalisation = list(ruv_mode = "automatic", auto_percentage_min = 5, auto_percentage_max = 30,
        auto_percentage_step = 5, k_penalty_weight = 0.05), da = list(q_value = 0.05, treat_lfc = 1L), ruviii_m_source = "technical_replicate",
        experiment_label = "kit_omics_lifecycle")
    corpus <- list(assays = list(counts = "/corpus/counts.tsv"), vendor = "trx_counts", design = "/corpus/design.tsv")
    for (id in c("import", "expression_filter", "qc_finalize", "normalisation", "da", "ruviii_export", "report", "session_export"))
        d$runKitTrxStep(h, driver, tempdir(), corpus, journey, list(id = id))
    expect_true("tile-transcriptomics" %in% seen$ids)
    expect_true(all(startsWith(setdiff(seen$ids, c("tile-transcriptomics", "btn-start-analysis", "btn-confirm-setup")), "transcriptomics_workflow-")))
    for (id in c("import-process_import", "qc-filter-apply_filter", "qc-s4_finalize-finalize_qc", "norm-run_normalization", "de-run_da_analysis",
        "summary-run_ruviii_export", "summary-generate_report", "summary-export_session_state"))
        expect_true(paste0("transcriptomics_workflow-", id) %in% seen$ids, info = id)
    expect_identical(seen$values[["transcriptomics_workflow-import-vendor_format"]], "trx_counts")
    expect_identical(seen$values[["transcriptomics_workflow-norm-auto_percentage_max"]], 30)
    expect_identical(seen$values[["transcriptomics_workflow-norm-k_penalty_weight"]], 0.05)
    expect_identical(seen$values[["transcriptomics_workflow-de-treat_lfc_cutoff"]], 1L)
    expect_identical(seen$values[["transcriptomics_workflow-summary-ruviii_m_source"]], "technical_replicate")
    expect_identical(seen$values[["transcriptomics_workflow-summary-run_ruviii_export#wait"]], "M = technical_replicate")
    expect_false(any(grepl("lipid|metab", seen$ids)))
    expect_error(d$runKitTrxStep(h, driver, tempdir(), corpus, journey, list(id = "intensity")), "unknown transcriptomics step")
})

test_that("KIT-098 correction: the technical-replicate summary is matched by its header and its sample set, in any run order", {
    d <- new.env(parent = globalenv())
    sys.source(testthat::test_path("../../tools/profiling/workflow_kit.R"), d)
    sys.source(testthat::test_path("../../tools/profiling/workflow_kit_driver.R"), d)
    # The builder lists a pair's runs in the design matrix's order: the counts file's order on the transcriptomics route
    # (kit-lifecycle-trx attempt 1 waited for design.tsv's order and timed out).
    text <- paste("Group: RPMI, Biological Replicate: 3\n  Samples: RPMI_03_r2, RPMI_03_r1\n  Technical Replicates: 2, 1",
        "Group: Sera, Biological Replicate: 1\n  Samples: Sera_01_r1, Sera_01_r2\n  Technical Replicates: 1, 2", sep = "\n\n")
    expect_true(d$kitTechRepSummaryMatches(text, "RPMI", 3, c("RPMI_03_r1", "RPMI_03_r2")))
    expect_true(d$kitTechRepSummaryMatches(text, "Sera", 1, c("Sera_01_r1", "Sera_01_r2")))
    expect_false(d$kitTechRepSummaryMatches(text, "RPMI", 3, c("RPMI_03_r1", "RPMI_04_r2")))
    expect_false(d$kitTechRepSummaryMatches(text, "RPMI", 1, c("RPMI_03_r1", "RPMI_03_r2")))
    expect_false(d$kitTechRepSummaryMatches("No technical replicates assigned yet.", "RPMI", 3, c("RPMI_03_r1", "RPMI_03_r2")))
    # The lipid and metab routes' summaries, in design.tsv's order, still match.
    expect_true(d$kitTechRepSummaryMatches("Group: A, Biological Replicate: 1\n  Samples: s1, s2\n  Technical Replicates: 1, 2", "A", 1, c("s1", "s2")))
})
