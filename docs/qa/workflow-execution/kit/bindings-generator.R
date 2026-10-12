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
            note = "core_018 source plus the KIT journey launch record, RUV-III technical-replicate structure and quantile-ceiling removal (user decisions 2026-10-09) plus the KIT-005 ceiling constant (launch-record deadline up to 14400 s, protocol v1.2) the once-per-sweep RUV canonical-correlation input, the sweep's k validity rule and the RUV-III guards (user decisions 2026-10-10; branch kit-control-ceiling 8e77d92d); the earlier installs browser-smoke.y8rDOi (9075b361), browser-smoke.sMOGTK (0f96397b) and browser-smoke.3nBUjx (42a7a4bb) are superseded and stay recorded; the KIT control until KIT-006 accepts")
        b$control <- "kit_001"
        # Protocol section 2 names the control's 2026-10-09 composition; the installed control also carries
        # the three 2026-10-10 user decisions (KIT-005 audit PH-05). The protocol text itself is not edited
        # (a protocol edit is a new version that re-runs controls).
        b$control_note <- "kit_001 = core_018 source (906c50dc) + KIT journey launch record, RUV-III technical-replicate structure, quantile-ceiling removal (2026-10-09) + launch-record deadline up to 14400 s (protocol v1.2), RUV sweep canonical-correlation input prepared once per sweep, sweep k validity rule and RUV-III guards (2026-10-10); source 8e77d92d on branch kit-control-ceiling; install browser-smoke.X9aj9P"
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
    member = list(programme = "kit-journey-v1", deadline = 14400),
    controller_deadline_seconds = 14400,
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
    # KIT-005 (protocol v1.2, user decision 2026-10-10): the first 300-sample attempt showed the data steps scaling
    # about 20 x from 8 to 300 samples (q-value 14 s -> 296 s, evidence 9 s -> 190 s), about 5,000 s for the whole
    # journey at the KIT-001 host speed and more on a shared host, so the 90 min cap was infeasible; the ceiling
    # is 4 h, the protocol v1.2 maximum.
    ceiling_seconds = 14400L,
    ceiling_note = "KIT-005 / protocol v1.2: data steps scale about 20 x from 8 to 300 samples (attempt 1: q-value 296 s, evidence 190 s); 4 h, the protocol maximum",
    # The baseline stage: two control journey members on kit_001 with the journey record (the first produces
    # the per-step populations and oracles, the second reproduces them) and the production-default member.
    pair = list(control = b$control, candidate = NULL, reps = 2L, production_default = TRUE))
# Operation-member lifecycle pair (KIT-005): the precursor rollup on screen_p14 through the installed action
# dispatch; the failure member fails inside the worker after the first ordered result stage (CORE-004's fault).
b$operation_lifecycle <- list(ceiling_seconds = 300L,
    fault = list(name = "after_aggregate", operation = "precursor_rollup", population = "screen_p14",
        failed_phase = "ordered_result_stage", error_pattern = "CORE-004 forced failure after the aggregate stage",
        mechanism = "the worker's phase hook raises after the ordered result stage completes; the action completes as worker_failed and the head stays on the input"))
# KIT-007: the registry changes the operation launch path, so its lifecycle pair runs on the kit_007 build.
if (!is.null(b$builds$kit_007)) b$operation_lifecycle$build <- "kit_007"
# Operation pairs per ticket (control and candidate builds; one or two (operation, population) pairs):
# KIT-006 binds the P14 rollup on full_p14 (P14 page oracle) and the sample filter on full_primary
# (the native filter oracle), both at the kit_001 control against the kit_006 candidate once installed.
b$tickets <- list()
if (!is.null(b$builds$kit_006)) {
    # KIT-006 accepted 2026-10-11 (stage kit-006: P14 0.991 x, sample filter 0.996 x the kit_001 control, both exact):
    # the KIT-006 build is the control from here (protocol section 2); the ticket binding keeps the control it ran against.
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-006", outcome = "control", build = "kit_006",
        note = "complete-sample partitioner (partitionBySamples) under the P14 rollup and the physical sample filter; accepted against kit_001 (stage kit-006, 2026-10-11): exact on both oracles, slice peak 0.991 x (P14) and 0.996 x (sample filter) the control; the control for KIT-007 onward")
    b$control <- "kit_006"
}
# The journey pair's control is the current control (KIT-006 audit PX-02: the pair was built above, before the
# promotion); KIT-005's own stage is recorded under journey$baseline with its control named there.
b$journey$pair$control <- b$control
if (!is.null(b$builds$kit_006)) b$tickets[["006"]] <- list(control = "kit_001", candidate = "kit_006",
    pairs = list(list(operation = "precursor_rollup", population = "full_p14"),
        list(operation = "minimum_peptides_per_sample", population = "full_primary")))
if (!is.null(journey_stage)) {
    harness <- new.env(parent = globalenv())
    sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
    result <- readRDS(file.path(journey_stage, "control-1", "journey-result.rds"))
    stopifnot(identical(result$status, "complete"), isTRUE(result$populations_complete))
    receipt <- readRDS(file.path(journey_stage, "stage.rds"))
    stopifnot(identical(receipt$status, "complete"), identical(receipt$decision$status, "baseline_recorded"))
    b$journey$baseline <- list(stage = journey_stage, member = "control-1", control = "kit_001",
        final_head_sha256 = result$final_head_sha256, verified_by = "control-2 (reproduced)", decision = receipt$decision$status)
    # The oracle fingerprints are recomputed from control-1's immutable store under the harness's current canonical
    # form (KIT-008 correction 2 dropped `core_utilisation`, the args-embedded support copies and the step description;
    # KIT-011 correction 1 dropped `declaration`), overlaid on the in-run record as the verifier overlays a candidate's:
    # class, rows, components, support and the args summary from the store, the state record's products recomputed,
    # the per-step assets the driver attached kept from the in-run fingerprint. Otherwise a candidate whose steps run
    # as declarations mismatches the oracles on bookkeeping, as the KIT-008 member did (KIT-012 pre-step).
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    stores <- unique(vapply(Filter(function(p) !is.null(p$state), result$populations), `[[`, character(1), "store"))
    stopifnot(length(stores) == 1L, dir.exists(stores))
    recomputed <- harness$kitFingerprintStore(stores, api)
    stopifnot(identical(recomputed$head_sha256, result$final_head_sha256))
    oracle_changes <- list()
    for (step in names(result$populations)) {
        population <- result$populations[[step]]
        if (is.null(population$state) || is.null(population$fingerprint)) next
        fresh <- recomputed$states[[population$state]]
        stopifnot(!is.null(fresh), identical(population$fingerprint$class, fresh$class), identical(population$fingerprint$rows, fresh$rows),
            setequal(names(population$fingerprint$components), names(fresh$components)),
            setequal(names(population$fingerprint$support), names(fresh$support)))
        changed <- character()
        for (field in c("components", "support", "summary_sha256"))
            if (!identical(population$fingerprint[[field]], fresh[[field]])) changed <- c(changed, field)
        for (name in names(fresh$products))
            if (!identical(population$fingerprint$products[[name]], fresh$products[[name]])) changed <- c(changed, paste0("products/", name))
        merged <- population$fingerprint
        for (field in c("components", "support", "summary_sha256")) merged[[field]] <- fresh[[field]]
        for (name in names(fresh$products)) merged$products[[name]] <- fresh$products[[name]]
        result$populations[[step]]$fingerprint <- merged
        if (length(changed)) oracle_changes[[step]] <- changed
    }
    b$journey$baseline$oracle_canonical_form <- list(
        recomputed_from = stores, date = format(Sys.Date()), harness = "tools/profiling/workflow_kit.R (kitFingerprintStore)",
        note = "KIT-012 pre-step: state fingerprints recomputed from control-1's store under the current canonical form (KIT-011 correction 1, KIT-008 correction 2) and overlaid on the in-run record; the in-run per-step assets are kept",
        changed = oracle_changes)
    for (step in names(result$populations)) {
        population <- result$populations[[step]]
        if (is.null(population$state) || is.null(population$fingerprint)) next
        key <- paste0("journey_", step)
        # The head after a step is the descriptor of a checkpoint entity in the immutable store; it is the
        # population's identity once the live head has moved on (KIT-005 audit PH-01). Records from the
        # attempt 6 driver predate the field, so it is looked up here.
        head_entity <- population$head_entity %||% harness$kitFindHeadEntity(population$store, population$head_sha256)
        if (is.null(head_entity)) stop("no checkpoint descriptor in ", population$store, " hashes to the head of step ", step)
        b$populations[[key]] <- list(step = step, state = population$state, store = population$store,
            head_sha256 = population$head_sha256, head_entity = head_entity, entity_id = population$entity_id, rows = population$rows,
            snapshot_time = population$snapshot_time, parameters = NULL,
            note = if (identical(step, "import")) "initial state: no table component, no products; nothing to verify" else
                if (identical(step, "session_export")) "the session export writes no state: same entity as journey_correlation under a later head" else NULL)
        b$oracles[[key]] <- population$fingerprint
    }
}
# KIT-005: the screen lifecycle ceiling is 600 s (was 300 s) for the harness proof of the per-member driver,
# the member-ceiling wait budget and the notification log: on the shared host of 2026-10-10 every step ran
# about twice as slowly as during KIT-001 (269 s), so 300 s stopped the check in the normalisation step.
b$lifecycle <- list(ceiling_seconds = 600L,
    fault = list(name = "guard_after_intensity", after_step = "intensity", failed_step = "sample",
        mechanism = "application-guard-inject.rds written after the intensity filter (the first post-rollup step) completes; the application monitor stops the slice (injected_guard_stop); the sample-filter step fails on disconnection; the head must equal the intensity-filter snapshot"))
# KIT-011: the journey pair as users run the GUI (no launch record): the kit_006 control (in-memory default,
# expected to end at the memory guard as KIT-005's production default did) against the kit_011 candidate
# (on-disk default), verified against the KIT-005 per-step oracles; whole-slice reduction with a per-step guard.
if (!is.null(b$builds$kit_011)) {
    b$tickets[["011"]] <- list(kind = "journey", control = "kit_006", candidate = "kit_011", record = FALSE)
    # The GUI's default state manager changed: the screen lifecycle member runs on kit_011 (audit PX-02).
    b$lifecycle$build <- "kit_011"
    # KIT-011 accepted 2026-10-11 (stage kit-journey-011: the kit_011 production default completed the 300-sample
    # journey at 5.51 GiB anonymous against the kit_006 default's 15.40 GiB at the memory guard, every state equal to
    # the KIT-005 oracles): the KIT-011 build is the control from here (protocol section 2).
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-011", outcome = "control", build = "kit_011",
        note = "on-disk native DIA state manager as the production default; accepted against kit_006 run as users run it (stage kit-journey-011, 2026-10-11): whole-slice anonymous peak 0.358 x the control's guard reading, every data step 0.23 to 0.48 x; the control for KIT-013 onward")
    b$control <- "kit_011"
    b$journey$pair$control <- b$control
}

# KIT-008: one `gui`-mode managed member on the KIT-008 build (production worker dispatch for P14 and the sample
# filter): the screen journey launched with no record under a managed application, verified against the KIT-011
# screen member (`kit-lifecycle-complete`), with the two worker dispatches required; 15 min ceiling (user amendment
# 2026-10-12, audit PE-03: the in-process screen reference took 486 s and each worker dispatch adds admission and acknowledgement).
if (!is.null(b$builds$kit_008)) {
    b$tickets[["008"]] <- list(kind = "gui", build = "kit_008", corpus = "screen", reference = "kit-lifecycle-complete",
        operations = c("precursor_rollup", "minimum_peptides_per_sample"), ceiling_seconds = 900L)
}

# KIT-009: the imputation step as an operation on the feature-group partitioner, control and candidate on the
# KIT-005 imputation population (input: the state after the replicate filter, `journey_replicate`; oracle: the
# state after imputation, `journey_imputation`), with the journey's recorded parameters (proportion 0.5, no
# exclusion column). Both arms run the one kit_009 build: `grouped` is the per-group loop as it stood before
# KIT-009, `partitioned` the partitioner (user decision pending at generation; see kit-009/outcome.md).
if (!is.null(b$builds$kit_009)) {
    b$tickets[["009"]] <- list(control = "kit_009", candidate = "kit_009",
        pairs = list(list(operation = "peptide_imputation", population = "journey_replicate", oracle = "journey_imputation",
            parameters = list(proportion_missing_values = 0.5, exclusion_column = NULL),
            variants = list(control = "grouped", candidate = "partitioned"))))
}

# KIT-012: the journey pair as users run the GUI (no launch record): the kit_008 control (the on-disk default with
# the production worker dispatch) against the kit_012 candidate (holders to references), decided on
# `resident_end_bytes` (candidate <= 0.95 x control) with the whole-slice anonymous peak not worse than 1.02 x and
# the per-step guard; verified against the KIT-005 per-step oracles under the current canonical form (above).
if (!is.null(b$builds$kit_012)) {
    b$tickets[["012"]] <- list(kind = "journey", control = "kit_008", candidate = "kit_012", record = FALSE, rule = "resident_end")
}

files <- c("tools/profiling/workflow_kit.R", "tools/profiling/workflow_kit_driver.R",
    "tools/profiling/workflow_kit_coordinator.R", "tools/profiling/workflow_kit_member.R", "tools/profiling/run_workflow_kit.sh",
    "tools/profiling/workflow_checkpoint_dia_driver.R", "tests/testthat/test-workflow-kit-harness.R",
    # The operation-member vehicle (KIT-005): the installed adapters test and its inputs, CORE's files unchanged.
    "tests/testthat/test-workflow-execution-adapters.R", "tests/testthat/helper-workflow-checkpoint-fixtures.R",
    "tests/testthat/helper-workflow-bulk-rollup.R", "tests/testdata/workflow-execution/qualification.json",
    "tests/testdata/workflow-checkpoints/scenarios.json", "tools/profiling/workflow_perf106.R",
    "tools/profiling/workflow_core.R")
# The checkpoint driver left the CORE pin at KIT-013 (the entry fixture's input became the generic `entry_hook`;
# its PERF-107 lifecycle scenarios are retired): it is hashed with the kit harness files above instead.
core_pinned <- c("tools/profiling/run_workflow_checkpoint_proof.R", "tools/profiling/run_workflow_checkpoint_proof.sh",
    "tools/profiling/workflow_perf106_repair.R")
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
