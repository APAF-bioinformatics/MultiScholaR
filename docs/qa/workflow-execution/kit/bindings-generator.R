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
# KIT-021: a corpus binds a file list (its report file first), the import format and the route its journey runs.
bind_corpus <- function(directory, manifest, files = c("report.tsv", "design.tsv", "proteins.fasta"), format = "diann", route = "DIA") {
    list(profile = manifest$profile, synthetic = TRUE, samples = manifest$samples, rows = manifest$rows,
        recipe = manifest$recipe, directory = directory, expected_retained_samples = manifest$expected_retained_samples,
        report_file = files[[1L]], format = format, route = route, notes = manifest$notes,
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
# KIT-021: the protein-level corpora in the real tool output formats (tools/profiling/kit_protein_corpus.R): a MaxQuant
# proteinGroups.txt and a Proteome Discoverer protein export, 300 samples x 8,000 proteins each (seeded), and the two
# e2e fixtures copied as the screen lifecycle corpora of the protein route. Generated once under the proof root.
protein_corpora <- list(
    lfq_300 = list(dir = "lfq-300", kind = "lfq", report = "proteinGroups.txt", format = "maxquant"),
    tmt_300 = list(dir = "tmt-300", kind = "tmt", report = "pd_proteins.tsv", format = "pd_tmt"),
    lfq_fixture = list(dir = "lfq-fixture", kind = "lfq_fixture", report = "proteinGroups.txt", format = "maxquant"),
    tmt_fixture = list(dir = "tmt-fixture", kind = "tmt_fixture", report = "pd_tmt_peptides.tsv", format = "pd_tmt"))
for (key in names(protein_corpora)) {
    spec <- protein_corpora[[key]]
    directory <- file.path(proof, "kit/corpus", spec$dir)
    if (!file.exists(file.path(directory, "manifest.json"))) {
        status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(repo, "tools/profiling/kit_protein_corpus.R")), spec$kind, shQuote(directory)))
        stopifnot(identical(status, 0L))
    }
    manifest <- jsonlite::read_json(file.path(directory, "manifest.json"))
    b$corpora[[key]] <- bind_corpus(directory, manifest, files = c(spec$report, "design.tsv", "proteins.fasta"), format = spec$format, route = "protein")
    b$corpora[[key]]$generator <- "tools/profiling/kit_protein_corpus.R"
    # A five-protein fixture cannot carry the automatic RUV sweep (no valid negative-control percentage): its
    # lifecycle members and oracle run the normalisation with RUV skipped; the 300-sample corpora keep the journey's settings.
    if (grepl("_fixture$", key)) b$corpora[[key]]$normalisation <- list(ruv_mode = "skip")
    b$corpora[[key]]$fasta_file <- "proteins.fasta"; b$corpora[[key]]$design_file <- "design.tsv"
}
# KIT-054: the user's real MaxQuant export (sensitive: bound by path and sha256 from its own directory, never copied into
# the repository; the manifest carries hashes and counts only, written by kit-054/bind-real-corpus.R) and, for the
# legacy lane on kit_036 (whose reader knows only the dotted column names), the header-renamed copy under the proof root
# (kit-054/dotted-copy.R). Both are skipped when absent so the bindings regenerate on any machine.
bind_real <- function(key, directory) {
    manifest_path <- file.path(directory, "manifest.json")
    if (!file.exists(manifest_path)) return(invisible(NULL))
    manifest <- jsonlite::read_json(manifest_path)
    files <- c(manifest$report_file, manifest$design_file, manifest$fasta_file)
    b$corpora[[key]] <<- bind_corpus(directory, manifest, files = files, format = manifest$format, route = "protein")
    b$corpora[[key]]$synthetic <<- FALSE; b$corpora[[key]]$sensitive <<- TRUE; b$corpora[[key]]$kind <<- manifest$kind
    b$corpora[[key]]$fasta_file <<- manifest$fasta_file; b$corpora[[key]]$design_file <<- manifest$design_file
    b$corpora[[key]]$design_groups <<- manifest$design_groups; b$corpora[[key]]$notes <<- manifest$note
}
bind_real("kc_lfq", Sys.getenv("KIT_054_LFQ_DIR", "/home/doktersmol/Projects/APAF/Projects/KC/MaxQuent - 12 Oct 2026"))
bind_real("kc_lfq_dotted", file.path(proof, "kit/corpus/kc-lfq-dotted"))
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
    # KIT-009 accepted 2026-10-12 (stage kit-009: both arms exact against journey_imputation, the partitioned arm
    # 0.962 x the grouped arm's slice peak and 0.952 x its anonymous peak, 625 batches each): the KIT-009 build is
    # the control from here (protocol section 2).
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-009", outcome = "control", build = "kit_009",
        note = "feature-group partitioner under the peptide-group family and the peptide_imputation operation; accepted against the pre-KIT-009 loop on the same build (stage kit-009, 2026-10-12): exact, slice peak 0.962 x, anonymous 0.952 x; the control for KIT-010 onward")
    b$control <- "kit_009"
    b$journey$pair$control <- b$control
}

# KIT-010: the protein replicate filter as an operation on the row partitioner, control and candidate on the
# KIT-005 protein population (input: the state after duplicate removal, `journey_duplicates`; oracle: the state
# after the replicate filter, `journey_protein_replicate`), with the journey's recorded parameters (grouping
# `group`, 4 cores recorded). Both arms run the one kit_010 build: `whole` is the GUI's in-process path, `rows`
# the partitioner (the KIT-009 precedent: the control build carries no such operation).
if (!is.null(b$builds$kit_010)) {
    b$tickets[["010"]] <- list(control = "kit_010", candidate = "kit_010",
        pairs = list(list(operation = "protein_replicate_filter", population = "journey_duplicates", oracle = "journey_protein_replicate",
            parameters = list(grouping_variable = "group", parallel_cores = 4L),
            variants = list(control = "whole", candidate = "rows"))))
    # KIT-010 accepted 2026-10-13 (stage kit-010, third launch: both arms exact against journey_protein_replicate, the
    # rows arm 0.994 x the whole arm's slice peak and 0.998 x its anonymous peak): the KIT-010 build is the control
    # from here (protocol section 2).
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-010", outcome = "control", build = "kit_010",
        note = "row partitioner and wide-assay sample projection under the protein replicate filter (protein_replicate_filter operation, protein worker roles); accepted against the GUI's in-process path on the same build (stage kit-010, 2026-10-13): exact, slice peak 0.994 x, anonymous 0.998 x; the control for KIT-014 onward")
    b$control <- "kit_010"
    b$journey$pair$control <- b$control
}

# KIT-014: the four peptide QC steps not measured by KIT-009 as single members on the one kit_014 build (protocol
# section 2: exactness and completion only, peak recorded, no ratio, no control member). Each runs on the KIT-005
# population the journey fed it with the parameters its oracle state records (the state's configuration, integers as the
# journey's inputs arrived; the driver's settings list names other values the GUI did not take), verified against the state the journey
# produced from it: the Q-value filter (raw design state -> `journey_qvalue`), the protein identification evidence
# filter (`journey_qvalue` -> `journey_protein_peptide`), the intensity filter (`journey_rollup` -> `journey_intensity`)
# and the replicate filter (`journey_sample` -> `journey_replicate`); the imputation step was measured by KIT-009.
if (!is.null(b$builds$kit_014)) {
    b$tickets[["014"]] <- list(kind = "single", build = "kit_014", members = list(
        list(operation = "peptide_qvalue_filter", population = "journey_design", oracle = "journey_qvalue", variant = "streamed",
            parameters = list(qvalue_threshold = 0.01, global_qvalue_threshold = 0.01, global_pg_qvalue_threshold = 0.01,
                proteotypic_only = TRUE, confidence_provenance_mode = "diann_main_report")),
        list(operation = "peptide_protein_evidence_filter", population = "journey_qvalue", oracle = "journey_protein_peptide",
            variant = "partitioned", parameters = list(num_peptides_per_protein_thresh = 1L, num_peptidoforms_per_protein_thresh = 2L)),
        list(operation = "peptide_intensity_filter", population = "journey_rollup", oracle = "journey_intensity", variant = "partitioned",
            parameters = list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L,
                peptides_intensity_cutoff_percentile = 1L, grouping_variable = "group")),
        list(operation = "peptide_replicate_filter", population = "journey_sample", oracle = "journey_replicate", variant = "partitioned",
            parameters = list(replicate_group_column = "replicates"))))
    # KIT-014 accepted 2026-10-13 (stage kit-014: all four single members exact against their journey oracles; peaks
    # recorded, no ratio): the KIT-014 build is the control from here (protocol section 2, latest accepted build).
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-014", outcome = "control", build = "kit_014",
        note = "peptide QC declarations (peptide_qvalue_filter, peptide_protein_evidence_filter, peptide_intensity_filter, peptide_replicate_filter) and the managed-step dispatch in the five peptide QC modules; accepted as four single members (stage kit-014, 2026-10-13): exact, anonymous peaks 2.47 / 2.37 / 2.24 / 2.23 GiB; the control for KIT-015 onward")
    b$control <- "kit_014"
    b$journey$pair$control <- b$control
}

# KIT-015: the four protein QC steps not measured by KIT-010 as single members on the one kit_015 build (protocol
# section 2: exactness and completion only). Each runs on the KIT-005 population the journey fed it with the parameters
# its oracle state records, verified against the journey's next state: the IQ rollup (`journey_imputation` ->
# `journey_iq`), the accession cleanup (`journey_iq` -> `journey_accession`; the journey's FASTA table from its project,
# read-only), the protein intensity filter (`journey_accession` -> `journey_protein_intensity`) and the duplicate removal
# (`journey_protein_intensity` -> `journey_duplicates`); the replicate filter was measured by KIT-010.
if (!is.null(b$builds$kit_015)) {
    if (is.null(journey_stage)) stop("kit_015 needs the kit_journey stage for the FASTA table")
    fasta <- file.path(journey_stage, "control-1/project/dia_baseline_smoke/scripts/proteomics/aa_seq_tbl_final.RDS")
    if (!file.exists(fasta)) stop("the control journey's FASTA table is missing: ", fasta)
    b$tickets[["015"]] <- list(kind = "single", build = "kit_015", members = list(
        list(operation = "protein_iq_rollup", population = "journey_imputation", oracle = "journey_iq", variant = "whole",
            parameters = list()),
        list(operation = "protein_accession_cleanup", population = "journey_iq", oracle = "journey_accession", variant = "whole",
            parameters = list(delimiter = ";", aggregation_method = "mean", fasta_table = fasta),
            # The bound file's hash: the case expansion refuses a changed FASTA table before any member runs (audit PE-03).
            inputs = list(fasta_table = list(path = fasta, sha256 = hash(fasta)))),
        list(operation = "protein_intensity_filter", population = "journey_accession", oracle = "journey_protein_intensity",
            variant = "whole", parameters = list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L,
                proteins_intensity_cutoff_percentile = 1L)),
        list(operation = "protein_duplicate_removal", population = "journey_protein_intensity", oracle = "journey_duplicates",
            variant = "whole", parameters = list(aggregation_method = "mean"))))
    # KIT-015 accepted 2026-10-14 (stage kit-015: all four single members exact against their journey oracles; peaks
    # recorded, no ratio): the KIT-015 build is the control from here (protocol section 2, latest accepted build).
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-015", outcome = "control", build = "kit_015",
        note = "protein QC declarations (protein_iq_rollup, protein_accession_cleanup, protein_intensity_filter, protein_duplicate_removal) and the managed-step dispatch in the five protein modules; accepted as four single members (stage kit-015, 2026-10-14): exact, anonymous peaks 2.57 / 1.19 / 1.18 / 1.34 GiB; the control for KIT-016 onward")
    b$control <- "kit_015"
    b$journey$pair$control <- b$control
}

# KIT-018: one single journey member on the kit_018 build (protocol section 2: exactness and completion only, no
# control, no ratio): the journey driven from import through the design step on the 300-sample corpus as users run
# it (no launch record), the import and design states verified against the KIT-005 `journey_import` and
# `journey_design` oracles, the peak recorded.
# `KIT_018_CORPUS=screen` binds the screen corpus for the lifecycle run of the stop-after path (pre-execution audit
# PE-02; not counted); the counted member binds the full corpus (the default).
if (!is.null(b$builds$kit_018)) {
    corpus_018 <- Sys.getenv("KIT_018_CORPUS", "full")
    stopifnot(corpus_018 %in% c("full", "screen"))
    b$tickets[["018"]] <- list(kind = "single_journey", build = "kit_018", corpus = corpus_018, stop_after = "design")
    # KIT-018 accepted 2026-10-14 (stage kit-journey-018: the single journey member complete through design, exact against
    # journey_design; 225 s, 3.16 GiB anonymous peak, no ratio): the KIT-018 build is the control from here (protocol section 2).
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-018", outcome = "control", build = "kit_018",
        note = "design default verified, reports from the accepted normalisation product, RUV and enrichment tables as paged display tables, legacy adoption verified; accepted as one single journey member (stage kit-journey-018, 2026-10-14): import and design exact, 225 s, 3.16 GiB anonymous peak; the control for KIT-019 onward")
    b$control <- "kit_018"
    b$journey$pair$control <- b$control
}

# KIT-028: the replay origin (user decision 2026-10-18): one single journey member on the kit_027 build (the first
# build that records step provenance with values) on the screen corpus as users run it (no launch record), import
# through DA, ceiling 900 s, one attempt. It measures nothing: its store is the record's origin, verified for
# self-consistency as a screen member is; the replay is in process (KIT-028 items 1-2).
if (!is.null(b$builds$kit_027)) {
    b$tickets[["028"]] <- list(kind = "single_journey", build = "kit_027", corpus = "screen", stop_after = "da",
        ceiling_seconds = 900L, purpose = "replay_origin")
}

# KIT-021: the LFQ and TMT single journey members on the kit_021 build (protocol v1.5: the whole journey of the protein
# route, no stop step; exactness against in-process oracles recorded on the same build and corpus; completion and the
# peaks recorded, no ratio, no control). Ticket key `021` is the LFQ member on the MaxQuant corpus and `021t` the TMT
# member on the Proteome Discoverer corpus (user instruction 2026-10-20: both corpora in the real tool output formats),
# each with its own ceiling (4 h, the protocol maximum). `KIT_021_CORPUS=lfq_fixture` / `KIT_021T_CORPUS=tmt_fixture`
# bind the fixture corpora for the screen lifecycle runs (not counted; ceiling 1,800 s).
if (!is.null(b$builds$kit_021)) {
    corpus_021 <- Sys.getenv("KIT_021_CORPUS", "lfq_300"); stopifnot(corpus_021 %in% c("lfq_300", "lfq_fixture"))
    corpus_021t <- Sys.getenv("KIT_021T_CORPUS", "tmt_300"); stopifnot(corpus_021t %in% c("tmt_300", "tmt_fixture"))
    b$tickets[["021"]] <- list(kind = "single_journey", build = "kit_021", corpus = corpus_021, route = "protein", whole_journey = TRUE,
        oracles = corpus_021, ceiling_seconds = if (identical(corpus_021, "lfq_300")) 14400L else 1800L)
    b$tickets[["021t"]] <- list(kind = "single_journey", build = "kit_021", corpus = corpus_021t, route = "protein", whole_journey = TRUE,
        oracles = corpus_021t, ceiling_seconds = if (identical(corpus_021t, "tmt_300")) 14400L else 1800L)
}
# KIT-054: the real-data members (3 single journeys): the new lane on the reader-repaired build with in-process oracles
# on the real corpus; the legacy in-memory lane on kit_036 on the header-renamed copy and on the synthetic tmt_300
# (no oracles: the harness's in-memory branch, completion and measures only).
if (!is.null(b$builds$kit_054) && !is.null(b$corpora$kc_lfq))
    b$tickets[["054"]] <- list(kind = "single_journey", build = "kit_054", corpus = "kc_lfq", route = "protein", whole_journey = TRUE,
        oracles = "kc_lfq", ceiling_seconds = 14400L)
if (!is.null(b$corpora$kc_lfq_dotted))
    b$tickets[["054l"]] <- list(kind = "single_journey", build = "kit_036", corpus = "kc_lfq_dotted", route = "protein", whole_journey = TRUE,
        ceiling_seconds = 14400L)
b$tickets[["054m"]] <- list(kind = "single_journey", build = "kit_036", corpus = "tmt_300", route = "protein", whole_journey = TRUE,
    ceiling_seconds = 14400L)
# KIT-053: the DIA_LIMPA single journey member (protocol v1.5: the whole journey of the route; exactness against the
# route's in-process oracles recorded on the same build and corpus by kit-053/oracle-record.R; completion and the peaks
# recorded, no ratio, no control). `KIT_053_CORPUS=screen` binds the screen corpus for the lifecycle run (not counted,
# ceiling 1,800 s); the counted member runs on the 300-sample DIA-NN corpus (`full`, ceiling 14,400 s). The route's
# normalisation tab settings are bound with the ticket: RUV covariates on both corpora; no between-sample normalisation
# on the screen corpus (the route's shipped default, the script's choice), cyclic loess on the full corpus. On the full
# synthetic corpus every sample of the `contrast` group carries a global offset (about 0.19 log2), so with `none` the
# group ANOVA is significant for every protein and the script's optimiser finds no control set at any percentage
# (kit-053/ruv-sweep-full-none.json, a valid unfavourable result kept); the route's optional between-sample
# normalisation step removes the offset and the sweep then finds its control sets (the DIA route's journey on the
# same corpus also normalises before its sweep). Kernels, engine and the counted claim are unchanged by the binding.
if (!is.null(b$builds$kit_053)) {
    corpus_053 <- Sys.getenv("KIT_053_CORPUS", "full"); stopifnot(corpus_053 %in% c("full", "screen"))
    b$tickets[["053"]] <- list(kind = "single_journey", build = "kit_053", corpus = corpus_053, route = "DIA_LIMPA", whole_journey = TRUE,
        oracles = paste0("limpa_", corpus_053), ceiling_seconds = if (identical(corpus_053, "full")) 14400L else 1800L,
        normalisation = list(norm_method = if (identical(corpus_053, "full")) "cyclicloess" else "none", ruv_mode = "covariates",
            correlation_threshold = 0.5))
}
# KIT-053 oracles: the route's per-step fingerprints on the screen and full DIA-NN corpora (the same corpora as the DIA
# route), recorded by kit-053/oracle-record.R into oracles-limpa_<corpus>.rds and bound under limpa_<corpus>_<step>.
for (key in c("screen", "full")) {
    if (is.null(b$corpora[[key]])) next
    oracle_file <- file.path(repo, "docs/qa/workflow-execution/kit/kit-053", paste0("oracles-limpa_", key, ".rds"))
    if (!file.exists(oracle_file)) next
    recorded <- readRDS(oracle_file)
    stopifnot(identical(recorded$corpus, key), identical(recorded$route, "DIA_LIMPA"), identical(recorded$corpus_sha256, b$corpora[[key]]$files[[1L]]$sha256))
    for (step in names(recorded$populations)) {
        population <- recorded$populations[[step]]
        b$populations[[paste0("limpa_", key, "_", step)]] <- list(step = step, state = population$state, store = recorded$store,
            head_sha256 = recorded$head_sha256, entity_id = population$entity_id, rows = population$rows, parameters = NULL,
            oracle_source = "in-process (protocol v1.5)", build = recorded$build, recorded_at = recorded$recorded_at, route = "DIA_LIMPA")
        b$oracles[[paste0("limpa_", key, "_", step)]] <- population$fingerprint
    }
}
# KIT-021 oracles (protocol v1.5 item 1): per-step fingerprints of an in-process run of the bound build on each protein
# corpus, recorded read-only by docs/qa/workflow-execution/kit/kit-021/oracle-record.R into oracles-<corpus>.rds; bound
# under <corpus>_<step> as the KIT-005 oracles are bound under journey_<step>. KIT-054's real corpus records into kit-054/.
# KIT-054 data handling rule: a sensitive corpus's fingerprints must not carry its tables' column names (the sample
# names); each column list becomes its length and sha256. The verifier judges components by their sha256 only.
redact_columns <- function(value) {
    if (is.list(value)) {
        if (is.character(value$columns)) value$columns <- list(n = length(value$columns), sha256 = digest::digest(value$columns, algo = "sha256"))
        return(lapply(value, redact_columns))
    }
    value
}
oracle_dirs <- c(stats::setNames(rep("kit-021", length(protein_corpora)), names(protein_corpora)), kc_lfq = "kit-054")
for (key in names(oracle_dirs)) {
    if (is.null(b$corpora[[key]])) next
    oracle_file <- file.path(repo, "docs/qa/workflow-execution/kit", oracle_dirs[[key]], paste0("oracles-", key, ".rds"))
    if (!file.exists(oracle_file)) next
    recorded <- readRDS(oracle_file)
    stopifnot(identical(recorded$corpus, key), identical(recorded$corpus_sha256, b$corpora[[key]]$files[[1L]]$sha256))
    for (step in names(recorded$populations)) {
        population <- recorded$populations[[step]]
        b$populations[[paste0(key, "_", step)]] <- list(step = step, state = population$state, store = recorded$store,
            head_sha256 = recorded$head_sha256, entity_id = population$entity_id, rows = population$rows, parameters = NULL,
            oracle_source = "in-process (protocol v1.5)", build = recorded$build, recorded_at = recorded$recorded_at)
        b$oracles[[paste0(key, "_", step)]] <- if (isTRUE(b$corpora[[key]]$sensitive)) redact_columns(population$fingerprint) else population$fingerprint
    }
}

# KIT-019: the DA publication through the managed worker (no counted member; one install recorded). The build is
# bound as the final stage-1 candidate for KIT-017's closure pair; the control stays the latest accepted build.
if (!is.null(b$builds$kit_019)) {
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-019", outcome = "candidate", build = "kit_019",
        note = "DA publication as the da_publication worker operation (the limma kernel stays in process); no counted member (probe recorded: Glimma absent from the application-side process); the final stage-1 build, KIT-017's candidate")
}

# KIT-022: kit debt before closure (no counted member; one install recorded). The build supersedes kit_019 as the
# final stage-1 candidate for KIT-017's closure pair; the control stays the latest accepted build.
if (!is.null(b$builds$kit_022)) {
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-022", outcome = "candidate", build = "kit_022",
        note = "DA arguments committed before the contrasts (manifest source = head), the active owner from the state scope, the Q-value admission estimate, the kernel contracts named; no counted member; the final stage-1 build, KIT-017's candidate (supersedes kit_019)")
}

# KIT-024: kit consistency before closure (no counted member; one install recorded). The build supersedes kit_022
# as the final stage-1 candidate for KIT-017's closure pair; the control stays the latest accepted build.
if (!is.null(b$builds$kit_024)) {
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-024", outcome = "candidate", build = "kit_024",
        note = "every declaration registered, server-side paged renders resolved from the accepted state, baseline test quarantine, KIT-023's global-free publication phases; no counted member; the final stage-1 build, KIT-017's candidate (supersedes kit_022)")
}

# KIT-017: the DIA stage closure pair (protocol section 4, ticket journey pair; ticket revision 1.1): the control is the
# latest accepted build (kit_018) and the candidate the final stage-1 build (kit_024), each one journey member as users
# run the GUI (no launch record) on the 300-sample corpus in one session, the candidate's states verified against the
# KIT-005 oracles. Rule `resident_end` (the one memory mechanism between the two builds is KIT-019's DA publication
# through the worker, whose claim is end-of-journey resident memory): candidate resident_end <= 0.95 x the control's,
# the whole-slice anonymous peak not worse than 1.02 x, the per-step guard from the first saved state. The stage's
# whole-slice reduction is reported against the production-default baselines (KIT-005 default-1 and KIT-011's
# control, both guard readings), not judged by this pair.
# `KIT_017_CORPUS=screen` binds the screen corpus for the pair's lifecycle run before the counted launch (both members
# self-verified on the screen corpus; not counted); the counted pair binds the full corpus (the default).
# The first screen run (2026-10-15, attempt 1) failed both members: the candidate's RUV grid page request found no
# state owner (Shiny serves the widget's page requests outside any reactive domain; the browser's error dialog then
# blocked the driver) and the control's application monitor died reading the worker registration file in the window
# after the controller unlinked it. Both are repaired in the package, so the candidate is the repaired install
# `kit_017` (the kit_024 source plus the two repairs) when it is bound; the control build keeps the monitor race.
if (!is.null(b$builds$kit_024)) {
    corpus_017 <- Sys.getenv("KIT_017_CORPUS", "full")
    stopifnot(corpus_017 %in% c("full", "screen"))
    candidate_017 <- if (!is.null(b$builds$kit_017)) "kit_017" else "kit_024"
    b$tickets[["017"]] <- list(kind = "journey", control = b$control, candidate = candidate_017, record = FALSE,
        rule = "resident_end", corpus = corpus_017)
    if (!is.null(b$builds$kit_017)) b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-017", outcome = "candidate",
        build = "kit_017", note = "kit_024 plus the screen-run repairs (page readers take the owner from the render and the served page is unnamed as DT serves it; the application monitor tolerates the worker registration vanishing); the installs core/browser-smoke.UpIfdA (9293655c) and core/browser-smoke.sOtQzu (c9635938) are superseded by the install from 7efb6ccd, which adds the managed step's completion context for the peptide role only (counted pair attempt 1, CP-02); the closure pair's candidate")
}

# KIT-036: self-contained pages without pandoc (the KIT-017 audit's root cause: pandoc's pass over the Glimma page is
# the DA step's peak, a 2.84 GiB child process on the 300-sample page; the R writer peaks at 1.13 GiB with no child).
# Journey pair as KIT-017 bound it: control kit_018, candidate kit_036 (kit_017's source plus KIT-037 and KIT-036),
# rule `resident_end` with the per-step guard; two journey members from the user (2026-10-16). `KIT_036_CORPUS=screen`
# binds the screen corpus for the pair's lifecycle run (not counted; KIT-037's live check).
if (!is.null(b$builds$kit_036)) {
    corpus_036 <- Sys.getenv("KIT_036_CORPUS", "full")
    stopifnot(corpus_036 %in% c("full", "screen"))
    b$tickets[["036"]] <- list(kind = "journey", control = b$control, candidate = "kit_036", record = FALSE,
        rule = "resident_end", corpus = corpus_036)
    # KIT-036 accepted 2026-10-16 (stage kit-journey-036: both members complete and verified; resident end 0.921,
    # whole slice 1.005, DA step 0.675 against kit_018): the KIT-036 build is the control from here (protocol section 2),
    # and stage 1 ends on an accepted build.
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-036", outcome = "control", build = "kit_036",
        note = "self-contained widget pages written in R (no pandoc pass) and the DA managed step acknowledged on its status text (KIT-037); accepted against kit_018 (stage kit-journey-036, 2026-10-16): resident end 0.921, whole slice 1.005, DA step 0.675")
    b$control <- "kit_036"
    b$journey$pair$control <- b$control
}

# KIT-035: technical replicates declared once (science; no counted member; one install recorded). The DA model and
# RUV-III take their replicate block from the object's declared `technical_replicate_id`; on the synthetic corpora the
# block equals the former `group_replicates` string, so the journey_da oracle is unchanged. The control stays kit_036.
if (!is.null(b$builds$kit_035)) {
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-035", outcome = "recorded", build = "kit_035",
        note = "the DA block and RUV-III's replicate structure from the one declaration (technical_replicate_id); a spanned block is standard lmFit with the reason recorded; DA tables identical on the screen project, the 300-sample oracle unchanged; no counted member, control unchanged (kit_036)")
}

# KIT-025: step provenance recorded at commit (no counted member; one install recorded with the source commit in the
# installed DESCRIPTION). The head gains step_provenance and software_environments; state records unchanged. Control
# stays kit_036.
if (!is.null(b$builds$kit_025)) {
    b$outcomes[[length(b$outcomes) + 1L]] <- list(ticket = "KIT-025", outcome = "recorded", build = "kit_025",
        note = "the head snapshot records each state's operation, variant, parameter names, build, source commit and defaults hash (legacy for states committed before); software environment per build; IQ fixed parameters declared; no counted member, control unchanged (kit_036)")
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
    # KIT-021 (protocol v1.5 item 5): the protein route's own step map.
    step_map_protein = "docs/qa/workflow-execution/kit/kit-001/step-map-protein.json",
    step_map_protein_sha256 = hash(file.path(repo, "docs/qa/workflow-execution/kit/kit-001/step-map-protein.json")),
    # KIT-053: the DIA_LIMPA route's step map (the route's rows and the shared steps), part of the route's case binding.
    step_map_limpa = "docs/qa/workflow-execution/kit/kit-001/step-map-limpa.json",
    step_map_limpa_sha256 = hash(file.path(repo, "docs/qa/workflow-execution/kit/kit-001/step-map-limpa.json")),
    core_pinned_unchanged = TRUE,
    files = as.list(stats::setNames(vapply(file.path(repo, files), hash, character(1), USE.NAMES = FALSE), files)))
jsonlite::write_json(b, path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
cat("bindings written:", path, "\nprotocol", b$harness$protocol_sha256, "\nstep map", b$harness$step_map_sha256, "\n")
