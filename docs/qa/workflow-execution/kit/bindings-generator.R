# Tracked generator for docs/qa/workflow-execution/kit/bindings.json (KIT-001; modelled on core/bindings-generator.R).
# Run from the repository root after any harness, protocol, step-map, corpus or build change; commit both files.
# Usage: Rscript docs/qa/workflow-execution/kit/bindings-generator.R [kit_001=<install run dir> ...] [kit_journey=<kit-journey stage dir>]
# The screen corpus (8 samples of the PERF-107 calibration recipe) is generated once under the proof
# root if absent; its hashes are then bound. The 300-sample corpora are PERF-107's generated cohorts.
args <- commandArgs(TRUE)
stopifnot(all(grepl("^kit_([0-9]{3}|journey)=.+", args)))
values <- stats::setNames(normalizePath(sub("^[^=]+=", "", args), mustWork = TRUE), sub("=.*$", "", args))
# kit_journey=<kit-journey stage directory>: ingest the control journey's per-step populations and oracles (KIT-005).
`%||%` <- function(x, y) if (is.null(x)) y else x
journey_stage <- if ("kit_journey" %in% names(values)) values[["kit_journey"]] else NULL
installs <- values[setdiff(names(values), "kit_journey")]
repo <- normalizePath(".")
stopifnot(file.exists(file.path(repo, "DESCRIPTION")))
proof <- "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5"
hash <- function(path) digest::digest(file = path, algo = "sha256")
core <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/core/bindings.json"))
path <- file.path(repo, "docs/qa/workflow-execution/kit/bindings.json")
b <- list(schema = 1L, owner = "KIT-001; extended by later KIT tickets", created = format(Sys.Date()),
    repository = core$repository)
# Builds: core_018 (the control until KIT-006 accepts) carried verbatim from the CORE bindings.
b$builds <- list(core_018 = core$builds$core_018)
b$outcomes <- list(list(ticket = "CORE-018", outcome = "superseded", build = "core_018",
    note = "post-closure CORE build (source 906c50dc, install browser-smoke.gQ40C7); superseded as the KIT control by kit_001 (user decision 2026-10-09: launch hook, RUV replicate structure, quantile ceiling removal)"))
b$control <- "core_018"
probe <- function(library) {
    out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(paste0(
        '.libPaths(c("', library, '", .libPaths())); ns <- asNamespace("MultiScholaR"); ',
        'cat(get("workflowExecutionBuildIdentity", ns)(), get("diaRollupNativeInfo", ns)()$manifest$binary_sha256, ',
        'normalizePath(getNamespaceInfo(ns, "path")), sep = "\\n")'))), stdout = TRUE)
    out <- tail(out, 3L)
    stopifnot(identical(out[[3L]], normalizePath(file.path(library, "MultiScholaR"))))
    list(build = out[[1L]], exact_sum = out[[2L]])
}
control_probe <- probe(b$builds$core_018$library)
stopifnot(identical(control_probe$build, b$builds$core_018$build), identical(control_probe$exact_sum, b$builds$core_018$exact_sum))
for (key in names(installs)) {
    install <- installs[[key]]
    library <- file.path(install, "library")
    new <- probe(library)
    b$builds[[key]] <- list(ticket = sub("^kit_", "KIT-", key), library = library, build = new$build, exact_sum = new$exact_sum,
        install_member = basename(install), source_commit = trimws(readLines(file.path(install, "source-commit.txt"))[[1L]]),
        runtime_patch_bytes = file.size(file.path(install, "runtime.patch")))
    if (key == "kit_001") {
        b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-001", outcome = "control", build = "kit_001",
            note = "core_018 source plus the KIT journey launch record, RUV-III technical-replicate structure and quantile-ceiling removal (user decisions 2026-10-09); the KIT control until KIT-006 accepts")
        b$control <- "kit_001"
    }
}
# CORE populations and oracles, copied verbatim (the P14 page oracle and the sample-filter oracle).
b$populations <- core$populations[c("screen_p14", "full_p14", "full_primary")]
b$oracles <- core$oracles
# Corpora: PERF-107 synthetic DIA-NN cohorts. Full = calibration (300 samples); heldout stays reserved.
synthetic <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/perf-107/synthetic-corpus.json"))
bind_corpus <- function(directory, manifest) {
    files <- c("report.tsv", "design.tsv", "proteins.fasta")
    list(profile = manifest$profile, synthetic = TRUE, samples = manifest$samples, rows = manifest$rows,
        recipe = manifest$recipe, directory = directory, expected_retained_samples = manifest$expected_retained_samples,
        files = stats::setNames(lapply(files, function(name) {
            file <- file.path(directory, name)
            stopifnot(identical(hash(file), manifest$files[[name]]$sha256))
            list(path = file, bytes = file.size(file), sha256 = manifest$files[[name]]$sha256)
        }), files))
}
full <- synthetic$profiles$calibration
b$corpora <- list(full = bind_corpus(dirname(full$files[["report.tsv"]]$path), full))
# Screen corpus v2 (2026-10-09): eight samples of the calibration recipe with dense counts, so the
# group-count intensity filter (two runs in two groups) keeps about 56% of features and the 3000
# cutoff drops exactly one sample. v1 (counts 500...60000) left 365 features and the sample filter
# emptied the table (attempt 3 of kit-lifecycle-complete); it stays on disk as the record.
screen_counts <- c(3000L, 20000L, 30000L, 40000L, 45000L, 50000L, 55000L, 60000L)
screen_dir <- file.path(proof, "kit/corpus/screen-v2")
if (!dir.exists(screen_dir)) {
    code <- new.env(parent = globalenv())
    sys.source(file.path(repo, "tests/testthat/helper-workflow-checkpoint-fixtures.R"), code)
    sys.source(file.path(repo, "tools/profiling/workflow_perf107_synthetic.R"), code)
    recipe <- code$perf107SyntheticRecipe("calibration")
    recipe$samples <- 8L
    recipe$counts <- screen_counts
    recipe$profile <- "calibration"
    dir.create(dirname(screen_dir), recursive = TRUE, showWarnings = FALSE)
    code$writePerf107SyntheticCohort(screen_dir, recipe)
}
screen <- jsonlite::read_json(file.path(screen_dir, "manifest.json"))
stopifnot(identical(screen$samples, 8L), identical(screen$recipe$cutoff, 3000L),
    identical(as.integer(unlist(screen$recipe$counts)), screen_counts))
b$corpora$screen <- bind_corpus(screen_dir, screen)
b$corpora$screen$note <- "screen corpus v2: eight samples of the calibration recipe (same universe, proteins and cutoff) with counts 3000, 20000, 30000, 40000, 45000, 50000, 55000, 60000; about 56% of features pass the group-count intensity filter and the 3000 cutoff drops one sample"
# Journey: frozen settings (audit A13/A14). GUI defaults except the sample-filter cutoff, which is the
# corpus's designed cutoff; the ruv_pearson profile of runCheckpointNormalizationActions; DA as runCheckpointDaActions.
b$journey <- list(
    experiment_label = "dia_baseline_smoke",
    state_manager_mode = "dia_native_sessions",
    member = list(programme = "kit-journey-v1", deadline = 5400),
    controller_deadline_seconds = 5400,
    claim_limits = c(
        "The journey runs the production native DIA session path (`initializeDIANativeImport`: on-disk reference states from import onward, in-process reference kernels) switched on by the KIT journey launch record, so every step is snapshotted after its own step.",
        "A production-default journey (no launch record, in-memory WorkflowState) has no per-step references; KIT-005 runs one such member as the memory baseline and cannot fingerprint it per step."),
    diagnostic_settings = FALSE,
    settings = list(sample = list("sample_filter-min_peptides_per_sample" = 3000)),
    normalisation = list(norm_method = "cyclicloess", ruv_mode = "automatic",
        ruv = list(auto_percentage_min = 1, auto_percentage_max = 20, separation_metric = "max_difference",
            k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_grouping_variable = "group"),
        correlation_threshold = 0.5, profile = "ruv_pearson"),
    da = list(formula = "~ 0 + group", q_value = 0.05, treat_lfc = 0),
    steps = c("import", "design", "qvalue", "protein_peptide", "rollup", "intensity", "sample", "replicate", "imputation",
        "iq", "accession", "protein_intensity", "duplicates", "protein_replicate", "normalisation_ruv", "correlation",
        "session_export", "da"),
    # KIT-005: the KIT-001 screen journey took 269 s for 320,872 report rows (8 samples); the full corpus has
    # 5,846,578 rows (18.2 x), so a linear scaling gives 4,900 s; the ceiling is the protocol's 90 min maximum.
    ceiling_seconds = 5400L,
    ceiling_note = "KIT-005: KIT-001 screen journey 269 s at 320,872 rows scaled linearly to 5,846,578 rows = 4,900 s; capped at the protocol maximum of 90 min",
    # The baseline stage: two control journey members on kit_001 with the journey record (the first produces
    # the per-step populations and oracles, the second reproduces them) and the production-default member.
    pair = list(control = b$control, candidate = NULL, reps = 2L, production_default = TRUE))
# Operation-member lifecycle pair (KIT-005): the precursor rollup on screen_p14 through the installed action
# dispatch; the failure member fails inside the worker after the first ordered result stage (CORE-004's fault).
b$operation_lifecycle <- list(ceiling_seconds = 300L,
    fault = list(name = "after_aggregate", operation = "precursor_rollup", population = "screen_p14",
        failed_phase = "ordered_result_stage", error_pattern = "CORE-004 forced failure after the aggregate stage",
        mechanism = "the worker's phase hook raises after the ordered result stage completes; the action completes as worker_failed and the head stays on the input"))
# Operation pairs per ticket (control and candidate builds, operation, population): added by KIT-006 onward.
b$tickets <- list()
if (!is.null(journey_stage)) {
    result <- readRDS(file.path(journey_stage, "control-1", "journey-result.rds"))
    stopifnot(identical(result$status, "complete"), isTRUE(result$populations_complete))
    receipt <- readRDS(file.path(journey_stage, "stage.rds"))
    stopifnot(identical(receipt$status, "complete"), identical(receipt$decision$status, "baseline_recorded"))
    b$journey$baseline <- list(stage = journey_stage, member = "control-1", final_head_sha256 = result$final_head_sha256,
        verified_by = "control-2 (reproduced)", decision = receipt$decision$status)
    for (step in names(result$populations)) {
        population <- result$populations[[step]]
        if (is.null(population$state) || is.null(population$fingerprint)) next
        key <- paste0("journey_", step)
        b$populations[[key]] <- list(step = step, state = population$state, store = population$store,
            head_sha256 = population$head_sha256, entity_id = population$entity_id, rows = population$rows,
            snapshot_time = population$snapshot_time, parameters = NULL)
        b$oracles[[key]] <- population$fingerprint
    }
}
b$lifecycle <- list(ceiling_seconds = 300L,
    fault = list(name = "guard_after_intensity", after_step = "intensity", failed_step = "sample",
        mechanism = "application-guard-inject.rds written after the intensity filter (the first post-rollup step) completes; the application monitor stops the slice (injected_guard_stop); the sample-filter step fails on disconnection; the head must equal the intensity-filter snapshot"))
files <- c("tools/profiling/workflow_kit.R", "tools/profiling/workflow_kit_driver.R",
    "tools/profiling/workflow_kit_coordinator.R", "tools/profiling/workflow_kit_member.R", "tools/profiling/run_workflow_kit.sh",
    "tools/profiling/workflow_checkpoint_dia_driver.R", "tests/testthat/test-workflow-kit-harness.R",
    # The operation-member vehicle (KIT-005): the installed adapters test and its inputs, CORE's files unchanged.
    "tests/testthat/test-workflow-execution-adapters.R", "tests/testthat/helper-workflow-checkpoint-fixtures.R",
    "tests/testthat/helper-workflow-bulk-rollup.R", "tests/testdata/workflow-execution/qualification.json",
    "tests/testdata/workflow-checkpoints/scenarios.json", "tools/profiling/workflow_perf106.R",
    "tools/profiling/workflow_core.R")
core_pinned <- c("tools/profiling/run_workflow_checkpoint_proof.R", "tools/profiling/run_workflow_checkpoint_proof.sh",
    "tools/profiling/workflow_checkpoint_dia_driver.R", "tools/profiling/workflow_perf106_repair.R")
for (file in core_pinned) stopifnot(identical(hash(file.path(repo, file)), core$harness$files[[file]]))
b$harness <- list(root = file.path(proof, "kit/runtime"),
    coordinator = "tools/profiling/workflow_kit_coordinator.R", driver = "tools/profiling/workflow_kit_driver.R",
    member = "tools/profiling/workflow_kit_member.R",
    runner = "tools/profiling/run_workflow_kit.sh", protocol = "docs/qa/workflow-execution/kit/protocol.md",
    protocol_sha256 = hash(file.path(repo, "docs/qa/workflow-execution/kit/protocol.md")),
    step_map = "docs/qa/workflow-execution/kit/kit-001/step-map.json",
    step_map_sha256 = hash(file.path(repo, "docs/qa/workflow-execution/kit/kit-001/step-map.json")),
    core_pinned_unchanged = TRUE,
    files = as.list(stats::setNames(vapply(file.path(repo, files), hash, character(1), USE.NAMES = FALSE), files)))
jsonlite::write_json(b, path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
cat("bindings written:", path, "\nprotocol", b$harness$protocol_sha256, "\nstep map", b$harness$step_map_sha256, "\n")
