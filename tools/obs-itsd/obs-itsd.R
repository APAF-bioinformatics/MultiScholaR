# Observation for finding K242-P01 (tickets/kit/audits/AUDIT-KIT-242-plan-20261125.md): does the internal-standard
# (ITSD) step change any value on the package's own lipidomics and metabolomics fixtures and simulators, which hold
# `ISTD` in the annotation column while the method matches the slot regex `^IS_` against that column?
# No package code. Branch obs-itsd-20261125, never merged. Decides only whether P01's reading holds per generator.
#
# Usage (run-tests scripts lane, from the package root; fetch_refs kit-190-ruv-control-selection and
# kit-188-ruv-default-calibration):
#   Rscript tools/obs-itsd/obs-itsd.R <out dir>                          the job
#   Rscript tools/obs-itsd/obs-itsd.R --group <group> <tree> <out rds>   one group of cases (child process)
#
# Each case builds the object as its journey does, runs the journey's steps before ITSD as the journey runs them
# (warnings and messages suppressed there, as the journeys do), then the ITSD step with the same function, arguments
# and route, with every warning and message caught (withCallingHandlers) and kept. Per assay it records: the largest
# absolute difference of a value before (state <prefix>_qc_complete) against after (<prefix>_itsd_norm), on the
# features both hold; the rows matching the object's slot regex in each non-sample column before and after; the
# method's own record (`@args$ITSDNormalization`: applied, method_type, pattern columns, features per assay); the
# warnings; and the analysis record's `itsd_applied` for the ITSD state.
#
# Groups:
#   staging  this tree, pkgload::load_all + the test helpers (as docs/qa/.../kit-172/t1b-t11.R loads them): the
#            declarations fixtures (lipidDeclarationsJourney / metabDeclarationsJourney, in process), KIT-172's F172
#            (T1a, T11; seed 172001) and T11's concentration fixture through kit172Journey's steps, and the POSITIVE
#            control: F172 lipid with the annotation column (LipidClass) of the five standard rows set from "ISTD" to
#            the row's own id ("IS_1" .. "IS_5"); nothing else changed.
#   worker   this tree: the declarations fixtures through omicsWorkerJourney's steps (helper-workflow-omics-worker.R),
#            the worker's entry point and acceptance in this process with the test's two mocks; the helper's
#            `omicsWorkerRun` is replaced by a copy whose only change is that warnings are not suppressed there (they
#            reach this script's handler, which records and muffles them).
#   k190     KIT-190's counted stage-1 tree (e344c691, effect.R SHA-256 7d76467a...), loaded as effect.R's loadTree
#            does; e2Data from that effect.R (by git show); journeys as kit172Journey (E2 cell N-r r 1, both classes;
#            aligned cell N-w r 1, lipid, as the bounded-k measurement ran it).
#   k188     this tree, installed, loaded as sim-k188.R loads it; the generator from sim-k188.R at the head of
#            kit-188-ruv-default-calibration (by git show); journeyL's steps (defaults from routeDefaults): L12 N1
#            data set 401 at masks none and vd20, L64 N1 401 none, K172c lipid and metab.

args <- commandArgs(trailingOnly = TRUE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))

K190_COMMIT <- "e344c69142269110b7c9918f8e01940d78a9120b"
K190_EFFECT <- "docs/qa/workflow-execution/kit/kit-190/effect.R"
K190_EFFECT_SHA <- "7d76467a37c9fcf890a557efa0ba8b26498052301f18f15f9231f87177cc5b5a"
K188_REFS <- c("origin/kit-188-ruv-default-calibration", "kit-188-ruv-default-calibration")
K188_SIM <- "docs/qa/workflow-execution/kit/kit-188/sim-k188.R"

oneLine <- function(x) gsub("[\t\r\n]+", " ", paste(x, collapse = " "))
quietly <- function(expr) suppressWarnings(suppressMessages(expr))

# ---- measurement (shared by every group) ----------------------------------------------------------------------------

#' A state's assays: the sample-value matrix, the ids and the regex count per non-sample column.
obsState <- function(owner, binding, state) {
  object <- owner$getState(state)
  if (is.null(object)) stop("state ", state, " absent")
  id <- methods::slot(object, binding$id_slot)
  runs <- as.character(object@design_matrix[[object@sample_id]])
  regex <- object@internal_standard_regex
  assays <- lapply(methods::slot(object, binding$data_slot), function(table) {
    table <- as.data.frame(table, check.names = FALSE, stringsAsFactors = FALSE)
    samples <- intersect(runs, names(table))
    m <- as.matrix(table[, samples, drop = FALSE])
    storage.mode(m) <- "double"
    others <- setdiff(names(table), samples)
    counts <- vapply(others, function(column) sum(grepl(regex, as.character(table[[column]]))), integer(1L))
    list(m = m, ids = as.character(table[[id]]), counts = counts, rows = nrow(table))
  })
  list(object = object, id = id, annotation = object@annotation_id_column, regex = regex, assays = assays)
}

#' Before against after on the features both hold (matched by id) and the samples both hold.
obsCompare <- function(before, after) {
  common <- intersect(after$ids, before$ids)
  common <- common[!is.na(common)]
  columns <- intersect(colnames(after$m), colnames(before$m))
  b <- before$m[match(common, before$ids), columns, drop = FALSE]
  a <- after$m[match(common, after$ids), columns, drop = FALSE]
  d <- a - b
  finite <- is.finite(d)
  positive <- is.finite(a) & is.finite(b) & a > 0 & b > 0
  list(common_rows = length(common), cells_compared = sum(finite), cells_changed = sum(finite & d != 0),
    max_abs_diff = if (any(finite)) max(abs(d[finite])) else NA_real_,
    max_abs_log2_ratio = if (any(positive)) max(abs(log2(a[positive] / b[positive]))) else NA_real_,
    missing_pattern_changed = sum(is.na(a) != is.na(b)))
}

#' Evaluate `expr` with every warning and message kept (and muffled), an error returned as the value.
obsCapture <- function(expr) {
  warnings <- character()
  messages <- character()
  value <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(value = value, warnings = warnings, messages = messages)
}

#' The analysis record's `itsd_applied` for the ITSD state.
obsRecordFlag <- function(owner, state) tryCatch({
  record <- quietly(buildWorkflowAnalysisRecord(owner))
  states <- vapply(record$steps, function(s) as.character(s$state %||% NA_character_), character(1L))
  at <- match(state, states)
  if (is.na(at)) "no step for the state" else {
    flag <- record$steps[[at]]$outcomes$itsd_applied
    if (is.null(flag)) "absent" else paste(format(flag), collapse = ",")
  }
}, error = function(e) paste("error:", oneLine(conditionMessage(e))))

#' Run the ITSD call on a prepared owner and measure it. Returns rows (one per assay), warnings and messages.
obsItsd <- function(case, generator, route, tree, owner, binding, itsd_call, note = "") {
  prefix <- binding$prefix
  before <- obsState(owner, binding, paste0(prefix, "_qc_complete"))
  captured <- obsCapture(itsd_call())
  step_error <- if (inherits(captured$value, "error")) oneLine(conditionMessage(captured$value)) else ""
  after <- tryCatch(obsState(owner, binding, paste0(prefix, "_itsd_norm")), error = function(e) NULL)
  itsd_args <- if (is.null(after)) NULL else after$object@args$ITSDNormalization
  flag <- obsRecordFlag(owner, paste0(prefix, "_itsd_norm"))
  field <- function(x) if (is.null(x)) "absent" else oneLine(format(x))
  counts <- function(state, assay) {
    if (is.null(state) || is.null(state$assays[[assay]])) return(c(id = NA_integer_, annotation = NA_integer_, all = NA_character_))
    k <- state$assays[[assay]]$counts
    c(id = unname(k[state$id]), annotation = unname(k[state$annotation]), all = paste(names(k), k, sep = "=", collapse = ";"))
  }
  rows <- lapply(names(before$assays), function(assay) {
    cmp <- if (is.null(after) || is.null(after$assays[[assay]])) NULL else obsCompare(before$assays[[assay]], after$assays[[assay]])
    cb <- counts(before, assay)
    ca <- counts(after, assay)
    features <- itsd_args$itsd_features_per_assay[[assay]]
    mine <- captured$warnings[grepl(sprintf("'%s'", assay), captured$warnings, fixed = TRUE)]
    data.frame(case = case, generator = generator, route = route, tree = tree, class = class(before$object)[[1L]],
      assay = assay, slot_regex = before$regex, id_column = before$id, annotation_column = before$annotation,
      rows_before = before$assays[[assay]]$rows, rows_after = if (is.null(after)) NA_integer_ else after$assays[[assay]]$rows %||% NA_integer_,
      regex_rows_id_before = cb[["id"]], regex_rows_id_after = ca[["id"]],
      regex_rows_annotation_before = cb[["annotation"]], regex_rows_annotation_after = ca[["annotation"]],
      regex_rows_by_column_before = cb[["all"]], regex_rows_by_column_after = ca[["all"]],
      common_rows = cmp$common_rows %||% NA_integer_, cells_compared = cmp$cells_compared %||% NA_integer_,
      cells_changed = cmp$cells_changed %||% NA_integer_, max_abs_diff = cmp$max_abs_diff %||% NA_real_,
      max_abs_log2_ratio = cmp$max_abs_log2_ratio %||% NA_real_,
      missing_pattern_changed = cmp$missing_pattern_changed %||% NA_integer_,
      value_changed = if (is.null(cmp)) NA else isTRUE(cmp$cells_changed > 0L),
      method_applied = field(itsd_args$applied), method_type = field(itsd_args$method_type),
      method_pattern_columns = field(itsd_args$itsd_pattern_columns), method_removed_itsd = field(itsd_args$removed_itsd),
      method_features_per_assay = if (is.null(features)) "absent" else paste(features, collapse = ";"),
      method_assays_with_features = if (is.null(itsd_args$itsd_features_per_assay)) "absent" else
        paste(names(itsd_args$itsd_features_per_assay), collapse = ";"),
      warnings_naming_assay = length(mine), warnings_total = length(captured$warnings),
      record_itsd_applied = flag, step_error = step_error, note = note, stringsAsFactors = FALSE)
  })
  keep <- grepl("ITSD|Processing assay|annotation column|identify", captured$messages, ignore.case = TRUE)
  list(rows = do.call(rbind, rows),
    warnings = if (length(captured$warnings)) data.frame(case = case, index = seq_along(captured$warnings),
      text = vapply(captured$warnings, oneLine, "", USE.NAMES = FALSE), stringsAsFactors = FALSE),
    messages = if (any(keep)) data.frame(case = case, text = vapply(captured$messages[keep], oneLine, "", USE.NAMES = FALSE), stringsAsFactors = FALSE))
}

obsFailed <- function(case, generator, route, tree, e) {
  list(rows = data.frame(case = case, generator = generator, route = route, tree = tree, class = NA_character_, assay = NA_character_,
    step_error = paste("case failed before measurement:", oneLine(conditionMessage(e))), stringsAsFactors = FALSE))
}

obsRun <- function(case, generator, route, tree, body) {
  cat(sprintf("case %s ...\n", case))
  t0 <- proc.time()[[3L]]
  out <- tryCatch(body(), error = function(e) obsFailed(case, generator, route, tree, e))
  cat(sprintf("case %s done in %.1f s\n", case, proc.time()[[3L]] - t0))
  out
}

#' Evaluate the top-level assignments of `names` from `file` into `env`, in file order; stop if one is absent.
obsDefinitions <- function(file, names, env = globalenv()) {
  found <- character()
  for (e in parse(file, keep.source = FALSE)) {
    if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) && as.character(e[[2L]]) %in% names) {
      eval(e, env)
      found <- c(found, as.character(e[[2L]]))
    }
  }
  absent <- setdiff(names, found)
  if (length(absent)) stop("definitions absent from ", file, ": ", paste(absent, collapse = ", "))
  invisible(found)
}

# ---- the journeys' steps before ITSD --------------------------------------------------------------------------------

#' kit172Journey (helper-kit172.R) up to ITSD: intensity filter (1, 0.5), QC finalise, ITSD (median, NULL).
kit172Case <- function(case, generator, tree, values, class, root, positive = FALSE, note = "") {
  object <- kit172Object(values, class)
  binding <- kit172Binding(class)
  if (positive) {
    # The positive control: the standard rows' annotation value set from "ISTD" to the row's own id; nothing else.
    assays <- methods::slot(object, binding$data_slot)
    id <- methods::slot(object, binding$id_slot)
    annotation <- object@annotation_id_column
    for (a in names(assays)) {
      hit <- startsWith(as.character(assays[[a]][[id]]), "IS_")
      assays[[a]][[annotation]][hit] <- as.character(assays[[a]][[id]][hit])
    }
    methods::slot(object, binding$data_slot) <- assays
    methods::validObject(object)
  }
  fixture <- quietly(kit172Owner(root, object, class))
  owner <- fixture$owner
  on.exit(try(closeDIAStateAdapter(owner), silent = TRUE), add = TRUE)
  config <- fixture$workflow$config_list
  prefix <- binding$prefix
  quietly({
    runWorkflowOmicsReferenceIntensityFilter(binding, owner, 1, 0.5, config)
    commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
  })
  obsItsd(case, generator, "in process: kit172Journey steps; runWorkflowOmicsReferenceItsdNormalisation(binding, owner, 'median', NULL, config)",
    tree, owner, binding, function() runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config), note)
}

#' lipidDeclarationsJourney / metabDeclarationsJourney up to ITSD, in process.
declarationsCase <- function(case, tree, class, root) {
  if (identical(class, "lipid")) {
    fixture <- quietly(lipidDeclarationsOwner(root))
    binding <- lipidomicsDeclarationBinding
  } else {
    fixture <- quietly(metabDeclarationsOwner(root))
    binding <- metabolomicsDeclarationBinding
  }
  owner <- fixture$owner
  on.exit(try(closeDIAStateAdapter(owner), silent = TRUE), add = TRUE)
  config <- fixture$workflow$config_list
  if (identical(class, "lipid")) {
    quietly({
      runLipidReferenceIntensityFilter(owner, 1, 0.5, config, hydrate = FALSE)
      runLipidReferenceDuplicateResolution(owner, config, hydrate = FALSE)
      commitLipidReferenceStoreStep(owner, "lipid_qc_finalize", "lipid_qc_complete", config, "QC processing complete - ready for normalization")
    })
    call <- function() runLipidReferenceItsdNormalisation(owner, "median", NULL, config, hydrate = FALSE)
    generator <- "helper-workflow-lipid-declarations.R lipidDeclarationsObject()"
    route <- "in process: lipidDeclarationsJourney steps; runLipidReferenceItsdNormalisation(owner, 'median', NULL, config, hydrate = FALSE)"
  } else {
    quietly({
      runMetabReferenceIntensityFilter(owner, 1, 0.5, config, hydrate = FALSE)
      runMetabReferenceDuplicateResolution(owner, config, hydrate = FALSE)
      commitWorkflowOmicsStoreStep(owner, "metab_qc_finalize", "metab_qc_complete", config, "QC processing complete - ready for normalization")
    })
    call <- function() runMetabReferenceItsdNormalisation(owner, "median", NULL, config, hydrate = FALSE)
    generator <- "helper-workflow-metab-declarations.R metabDeclarationsObject()"
    route <- "in process: metabDeclarationsJourney steps; runMetabReferenceItsdNormalisation(owner, 'median', NULL, config, hydrate = FALSE)"
  }
  obsItsd(case, generator, route, tree, owner, binding, call)
}

#' omicsWorkerJourney (helper-workflow-omics-worker.R) up to ITSD: the worker's entry point and acceptance in process.
workerCase <- function(case, tree, class, root) {
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
    finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) }, .package = "MultiScholaR")
  if (identical(class, "lipid")) {
    fixture <- quietly(lipidDeclarationsOwner(file.path(root, "lipid")))
    binding <- lipidomicsDeclarationBinding
  } else {
    fixture <- quietly(metabDeclarationsOwner(file.path(root, "metab")))
    binding <- metabolomicsDeclarationBinding
  }
  owner <- fixture$owner
  on.exit(try(closeDIAStateAdapter(owner), silent = TRUE), add = TRUE)
  config <- fixture$workflow$config_list
  controller <- omicsWorkerController(fixture)
  commit <- function(step) omicsWorkerAction(fixture, controller, step$operation, step$parameters,
    list(state_name = step$state_name, config = config, description = step$description))
  quietly({
    commit(omicsIntensityFilterStep(binding, 1, 0.5))
    commit(omicsDuplicateResolutionStep(binding))
    commitWorkflowOmicsStoreStep(owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"), config,
      "QC processing complete - ready for normalization")
  })
  obsItsd(case, sprintf("helper-workflow-%s-declarations.R %sDeclarationsObject() via omicsWorkerFixtures", class, class),
    "worker entry point + acceptance in process (omicsWorkerJourney steps): omicsWorkerAction(omicsItsdNormalisationStep(binding, 'median', NULL))",
    tree, owner, binding, function() commit(omicsItsdNormalisationStep(binding, "median", NULL)))
}

#' sim-k188.R journeyL up to ITSD (filter at the route's defaults, duplicate resolution, QC finalise, ITSD d$itsd NULL).
k188Case <- function(case, tree, sim, mask, root) {
  omic <- sim$p$omic
  d <- tryCatch(routeDefaults(omic), error = function(e) e)
  defaults_note <- "defaults: routeDefaults()"
  if (inherits(d, "error")) {
    defaults_note <- paste("routeDefaults failed (", oneLine(conditionMessage(d)), "); EXPECTED_DEFAULTS used")
    d <- EXPECTED_DEFAULTS[[omic]]
  }
  old <- options(MultiScholaR.ruv_design_panel = FALSE)
  on.exit(options(old), add = TRUE)
  values <- maskedValues(sim, mask)
  binding <- bindingOf(omic)
  prefix <- binding$prefix
  object <- omicObject(sim, values)
  fixture <- quietly(if (identical(omic, "lipid")) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object))
  fixture$workflow$contrasts_tbl <- contrastsOf(sim)
  owner <- fixture$owner
  on.exit(try(closeDIAStateAdapter(owner), silent = TRUE), add = TRUE)
  config <- fixture$workflow$config_list
  quietly({
    runWorkflowOmicsReferenceIntensityFilter(binding, owner, d$percentile, d$proportion, config)
    runWorkflowOmicsReferenceDuplicateResolution(binding, owner, config)
    commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
  })
  obsItsd(case, sprintf("sim-k188.R %s (cell %s, scenario %s, data set %s, mask %s)",
      if (identical(sim$p$idx, 0L)) "k172cSim" else "generate", sim$cell, sim$scenario, sim$r, mask),
    sprintf("in process: journeyL steps; runWorkflowOmicsReferenceItsdNormalisation(binding, owner, '%s', NULL, config)", d$itsd),
    tree, owner, binding, function() runWorkflowOmicsReferenceItsdNormalisation(binding, owner, d$itsd, NULL, config),
    sprintf("%s: percentile %s, proportion %s, itsd %s", defaults_note, d$percentile, d$proportion, d$itsd))
}

# ---- groups (child process) -----------------------------------------------------------------------------------------

#' A copy of helper-workflow-omics-worker.R's omicsWorkerRun whose only change is that `suppressWarnings` is removed.
obsWorkerRun <- function(fixture, operation, parameters = list(), directory = NULL, token = list(request_id = "wfx-test"), config = list(mode = "adapters")) {
  owner <- fixture$owner
  store <- owner$.context$store
  if (is.null(directory)) directory <- tempfile("request-", tmpdir = fixture$source_dir)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
    input_references = list(state = getWorkflowStateReference(owner)), parameters = parameters,
    token = token, owner = store$identity, project_root = dirname(store$directory),
    directory = directory, config = config, worker_cap = 12 * 1024^3)
  path <- file.path(directory, "request.rds")
  writeWorkflowExecutionMessage(path, request)
  withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"), suppressMessages(utils::capture.output(try(
    workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))), silent = TRUE))))
  result <- readRDS(file.path(directory, "worker-result.rds"))
  result$directory <- directory
  result
}

runGroup <- function(group, tree) {
  work <- file.path(tempdir(), paste0("obs-", group))
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  root <- function(case) file.path(work, case)
  here <- system2("git", c("-C", shQuote(tree), "rev-parse", "HEAD"), stdout = TRUE)
  results <- list()
  add <- function(x) results[[length(results) + 1L]] <<- x
  if (group %in% c("staging", "worker")) {
    suppressMessages(pkgload::load_all(tree, export_all = TRUE, helpers = FALSE, attach_testthat = FALSE, quiet = TRUE))
    helpers <- new.env(parent = globalenv())
    suppressMessages(testthat::source_test_helpers(file.path(tree, "tests", "testthat"), env = helpers))
    attach(helpers, name = "obs_helpers", warn.conflicts = FALSE)
    if (identical(group, "worker")) {
      # VAL-54 instrument repair 1 (run 45579225993: sourcing every helper under the namespace hit a locked
      # environment): the three worker functions alone re-homed under the namespace, as in a test, so that the
      # test's two mocks reach them; their bodies unchanged except omicsWorkerRun's suppressWarnings (see above).
      worker_env <- new.env(parent = asNamespace("MultiScholaR"))
      for (name in c("omicsWorkerAction", "omicsWorkerController")) {
        fn <- get(name, envir = helpers)
        environment(fn) <- worker_env
        assign(name, fn, envir = worker_env)
      }
      environment(obsWorkerRun) <- worker_env
      assign("omicsWorkerRun", obsWorkerRun, envir = worker_env)
      attach(worker_env, name = "obs_worker", warn.conflicts = FALSE)
    }
  }
  if (identical(group, "staging")) {
    f172 <- kit172Values(172001L)$values
    concentration <- kit172ConcentrationValues()
    for (class in c("lipid", "metab")) {
      add(obsRun(paste0("decl_", class), "declarations", "in process", here, function() declarationsCase(paste0("decl_", class), here, class, root(paste0("decl_", class)))))
    }
    for (class in c("lipid", "metab")) {
      case <- paste0("k172_f172_", class)
      add(obsRun(case, "F172", "kit172Journey", here, function() kit172Case(case,
        "helper-kit172.R kit172Values(172001L)$values (F172; KIT-172 T1a at s = 1, T11, KIT-190 E1 t11_f172)", here, f172, class, root(case))))
      case <- paste0("k172_conc_", class)
      add(obsRun(case, "T11 concentration", "kit172Journey", here, function() kit172Case(case,
        "helper-kit172.R kit172ConcentrationValues() (KIT-172 T11; KIT-190 E1 t11_concentration)", here, concentration, class, root(case))))
    }
    case <- "positive_k172_f172_lipid"
    add(obsRun(case, "F172 positive control", "kit172Journey", here, function() kit172Case(case,
      "POSITIVE CONTROL: F172 (seed 172001), lipid, LipidClass of the 5 IS_ rows set from 'ISTD' to the row's own id", here, f172,
      "lipid", root(case), positive = TRUE, note = "positive control")))
  }
  if (identical(group, "worker")) {
    for (class in c("lipid", "metab")) {
      case <- paste0("worker_", class)
      add(obsRun(case, "worker", "worker", here, function() workerCase(case, here, class, root(case))))
    }
  }
  if (identical(group, "k190")) {
    effect <- Sys.getenv("OBS_K190_EFFECT")
    # effect.R's loadTree (the tree compiled first, as effect.R does before its sides)
    pkgbuild::compile_dll(tree, quiet = TRUE)
    suppressPackageStartupMessages(pkgload::load_all(tree, export_all = TRUE, helpers = FALSE, quiet = TRUE))
    logger::log_threshold(logger::FATAL)
    helpers <- new.env(parent = globalenv())
    suppressMessages(testthat::source_test_helpers(file.path(tree, "tests", "testthat"), env = helpers))
    attach(helpers, name = "k190_helpers", warn.conflicts = FALSE)
    obsDefinitions(effect, c("e2PickCells", "e2Values", "E2_CELLS", "ALIGNED_CELLS", "%||%", "e2Data"))
    for (spec in list(list("N-r", "lipid"), list("N-r", "metab"), list("N-w", "lipid"))) {
      case <- sprintf("k190_%s_r1_%s", gsub("-", "", spec[[1L]]), spec[[2L]])
      data <- e2Data(spec[[1L]], 1L)
      add(obsRun(case, "KIT-190", "kit172Journey", here, function() kit172Case(case,
        sprintf("KIT-190 effect.R (SHA-256 %s) e2Data('%s', 1) (seed %s)", substr(K190_EFFECT_SHA, 1L, 8L), spec[[1L]], data$seed),
        here, data$values, spec[[2L]], root(case))))
    }
  }
  if (identical(group, "k188")) {
    sim_file <- Sys.getenv("OBS_K188_SIM")
    # sim-k188.R's loading, at its lines 41-49
    suppressPackageStartupMessages(library(MultiScholaR))
    attach(asNamespace("MultiScholaR"), name = "MultiScholaR:namespace", warn.conflicts = FALSE)
    for (helper in c("helper-workflow-omics-owner.R", "helper-workflow-lipid-declarations.R", "helper-workflow-metab-declarations.R",
        "helper-kit172.R")) {
      sys.source(file.path(tree, "tests", "testthat", helper), envir = globalenv())
    }
    logger::log_threshold(logger::FATAL)
    options(MultiScholaR.ruv_design_panel = FALSE)
    obsDefinitions(sim_file, c("GROUPS", "SCENARIOS", "MED_VAR_L", "MED_VAR_P", "cellDef", "CELLS", "standardise", "residStd",
      "generate", "solveLod", "maskOf", "maskedValues", "contrastsOf", "omicObject", "bindingOf", "builderValue", "routeDefaults",
      "EXPECTED_DEFAULTS", "k172cSim"))
    # sim-k188.R's next line after EXPECTED_DEFAULTS (not a plain assignment, so obsDefinitions leaves it)
    defaults <- get("EXPECTED_DEFAULTS", envir = globalenv())
    defaults$metab <- defaults$lipid
    assign("EXPECTED_DEFAULTS", defaults, envir = globalenv())
    specs <- list(list("L12", "N1", 401L, "none"), list("L12", "N1", 401L, "vd20"), list("L64", "N1", 401L, "none"))
    for (spec in specs) {
      case <- sprintf("k188_%s_%s_%d_%s", spec[[1L]], spec[[2L]], spec[[3L]], spec[[4L]])
      add(obsRun(case, "KIT-188", "journeyL", here, function() k188Case(case, here, generate(spec[[1L]], spec[[2L]], spec[[3L]]),
        spec[[4L]], root(case))))
    }
    for (class in c("lipid", "metab")) {
      case <- paste0("k188_K172c_", class)
      add(obsRun(case, "KIT-188", "journeyL", here, function() k188Case(case, here, k172cSim(class), "none", root(case))))
    }
  }
  results
}

if (length(args) >= 1L && identical(args[[1L]], "--group")) {
  stopifnot(length(args) == 4L)
  result <- tryCatch(runGroup(args[[2L]], normalizePath(args[[3L]])), error = function(e) list(error = oneLine(conditionMessage(e))))
  saveRDS(result, args[[4L]])
  quit(save = "no", status = 0L)
}

# ---- the job (parent) -----------------------------------------------------------------------------------------------
stopifnot(length(args) %in% 1:2)
out <- args[[1L]]
# optional: a comma list of the groups to run (default all; the VAL-54 rerun runs only the failed group)
wanted <- if (length(args) == 2L) strsplit(args[[2L]], ",", fixed = TRUE)[[1L]] else c("staging", "worker", "k188", "k190")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
git <- function(...) suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
head_commit <- git("rev-parse", "HEAD")[[1L]]
sha256 <- function(f) digest::digest(file = f, algo = "sha256")
provenance <- c(sprintf("head\t%s", head_commit), sprintf("R\t%s", R.version.string), sprintf("started\t%s", format(Sys.time(), tz = "UTC")))

# KIT-190: the counted stage-1 tree and its effect.R, by git show
k190_tree <- file.path(tempdir(), "k190-tree")
k190_add <- git("worktree", "add", "--detach", shQuote(k190_tree), K190_COMMIT)
k190_effect <- file.path(tempdir(), "effect-k190.R")
system2("git", c("show", paste0(K190_COMMIT, ":", K190_EFFECT)), stdout = k190_effect)
k190_sha <- sha256(k190_effect)
provenance <- c(provenance, sprintf("k190_commit\t%s", K190_COMMIT), sprintf("k190_worktree\t%s", oneLine(k190_add)),
  sprintf("k190_effect_sha256\t%s", k190_sha), sprintf("k190_effect_sha256_expected\t%s", K190_EFFECT_SHA))

# KIT-188: sim-k188.R at the head of its branch, by git show
k188_ref <- NA_character_
for (ref in K188_REFS) {
  hit <- git("rev-parse", "--verify", "--quiet", paste0(ref, "^{commit}"))
  if (length(hit) && grepl("^[0-9a-f]{40}$", hit[[1L]])) { k188_ref <- hit[[1L]]; break }
}
k188_sim <- file.path(tempdir(), "sim-k188.R")
if (!is.na(k188_ref)) system2("git", c("show", paste0(k188_ref, ":", K188_SIM)), stdout = k188_sim)
provenance <- c(provenance, sprintf("k188_commit\t%s", k188_ref),
  sprintf("k188_sim_sha256\t%s", if (file.exists(k188_sim)) sha256(k188_sim) else NA))

groups <- list(list(name = "staging", tree = getwd(), minutes = 10),
  list(name = "worker", tree = getwd(), minutes = 8),
  list(name = "k188", tree = getwd(), minutes = 10),
  list(name = "k190", tree = k190_tree, minutes = 18))
rscript <- file.path(R.home("bin"), "Rscript")
status <- list()
collected <- list()
for (g in Filter(function(g) g$name %in% wanted, groups)) {
  rds <- file.path(tempdir(), paste0("child-", g$name, ".rds"))
  log <- file.path(out, paste0("child-", g$name, ".txt"))
  t0 <- proc.time()[[3L]]
  ok_inputs <- switch(g$name, k190 = dir.exists(g$tree) && identical(k190_sha, K190_EFFECT_SHA), k188 = file.exists(k188_sim), TRUE)
  if (!ok_inputs) {
    status[[g$name]] <- data.frame(group = g$name, exit = NA_integer_, seconds = 0, error = "inputs absent or SHA-256 mismatch")
    next
  }
  exit <- withr::with_envvar(c(OBS_K190_EFFECT = k190_effect, OBS_K188_SIM = k188_sim),
    system2("timeout", c("--kill-after=30", paste0(g$minutes, "m"), shQuote(rscript), shQuote(script), "--group", g$name,
      shQuote(g$tree), shQuote(rds)), stdout = log, stderr = log))
  result <- if (file.exists(rds)) readRDS(rds) else list(error = "no result file (timeout or crash)")
  status[[g$name]] <- data.frame(group = g$name, exit = exit, seconds = round(proc.time()[[3L]] - t0, 1),
    error = oneLine(result$error %||% ""), stringsAsFactors = FALSE)
  if (is.null(result$error)) collected <- c(collected, result)
  cat(sprintf("group %s: exit %s, %.0f s, %s\n", g$name, exit, proc.time()[[3L]] - t0, result$error %||% "ok"))
}

bind <- function(part) {
  frames <- Filter(Negate(is.null), lapply(collected, `[[`, part))
  if (!length(frames)) return(NULL)
  columns <- unique(unlist(lapply(frames, names)))
  do.call(rbind, lapply(frames, function(f) { for (c in setdiff(columns, names(f))) f[[c]] <- NA; f[columns] }))
}
write_tsv <- function(x, name) if (!is.null(x)) utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA")
rows <- bind("rows")
write_tsv(rows, "itsd-observation.tsv")
write_tsv(bind("warnings"), "itsd-warnings.tsv")
write_tsv(bind("messages"), "itsd-messages.tsv")
write_tsv(do.call(rbind, status), "groups.tsv")
writeLines(c(provenance, sprintf("finished\t%s", format(Sys.time(), tz = "UTC"))), file.path(out, "provenance.txt"))
if (!is.null(rows)) print(rows[, intersect(c("case", "assay", "value_changed", "max_abs_diff", "regex_rows_id_before", "regex_rows_id_after",
  "regex_rows_annotation_before", "method_applied", "method_features_per_assay", "warnings_naming_assay", "record_itsd_applied",
  "step_error"), names(rows))], right = FALSE)
git("worktree", "remove", "--force", shQuote(k190_tree))
