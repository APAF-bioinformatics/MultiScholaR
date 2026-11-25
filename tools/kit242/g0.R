# KIT-242 gate G0 (ticket r1.2, section 6): observe, at the base commit, each case in which a lipidomics or
# metabolomics correction does not happen but is recorded and reported as done. No package code: branch kit-242-g0
# carries only this script and is never merged. Nothing here judges a case: every case records what the step returns,
# the conditions it raises (warning, message, error, logger line), what the tab says, what the record says and what the
# differential analysis (DA) then runs on; the comparison with the ticket's predictions is made by reading the tables.
#
# Usage (run-tests scripts lane, from the package root; fetch_refs kit-188-ruv-default-calibration):
#   Rscript tools/kit242/g0.R <out dir>                        the job: each group in its own child process
#   Rscript tools/kit242/g0.R --group <group> <tree> <rds>     one group (child)
#
# Groups: core (C1 declaration and worker arms, C2, C3, C3m, C4, C6, C7, C9, C11, C12 a/b and a-decl, C14),
#   tab (the managed tab of C1, C3, C3m, C12 a-decl), replay (C8), lanes (C5, C13), k188 (C12 c).
# C10 is not run here: it reuses run 45579225993 (obs-itsd, package commit 4a82ddff = 3eeaf100 plus its script).
# Seeds: case k draws with set.seed(242000 + k, "Mersenne-Twister", "Inversion", "Rejection") before any package call.

args <- commandArgs(trailingOnly = TRUE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
K188_REFS <- c("origin/kit-188-ruv-default-calibration", "kit-188-ruv-default-calibration")
K188_SIM <- "docs/qa/workflow-execution/kit/kit-188/sim-k188.R"

# ---- the case definitions, evaluated in an environment whose parent is the package namespace (as a test file's) -----
caseDefinitions <- quote({

oneLine <- function(x) gsub("[\t\r\n]+", " ", paste(as.character(x), collapse = " "))
SEED <- function(k) set.seed(242000L + k, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
OUT <- new.env()
OUT$obs <- list(); OUT$cond <- list(); OUT$log <- list(); OUT$notes <- list()
LOGBUF <- new.env()
LOGBUF$lines <- character()

txt <- function(x) {
  if (is.null(x)) return("NULL")
  if (inherits(x, "condition")) return(oneLine(paste(class(x)[[1L]], conditionMessage(x), sep = ": ")))
  if (is.list(x)) {
    return(tryCatch(oneLine(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "string", digits = NA, force = TRUE)),
      error = function(e) oneLine(utils::capture.output(utils::str(x, max.level = 3L)))))
  }
  if (!length(x)) return("<empty>")
  oneLine(paste(format(x), collapse = ";"))
}
obs <- function(case, class, arm, item, value) {
  OUT$obs[[length(OUT$obs) + 1L]] <- data.frame(case = case, class = class, arm = arm, item = item, value = txt(value),
    stringsAsFactors = FALSE)
}
cond <- function(case, class, arm, step, kind, cls, text) {
  OUT$cond[[length(OUT$cond) + 1L]] <- data.frame(case = case, class = class, arm = arm, step = step, kind = kind, cls = cls,
    text = oneLine(text), stringsAsFactors = FALSE)
}
LOGGER_KEEP <- paste("^(WARN|ERROR)|Skipping|RUV correction applied|Dimensions before cleaning|No assays were|Returning original",
  "Using [0-9]+ control|Invalid K|control features|Cleaning design matrix", sep = "|")

#' Evaluate `expr`, keeping every warning, the relevant messages, the logger lines and an error (returned as the value).
capture <- function(expr, case, class, arm, step) {
  LOGBUF$lines <- character()
  warns <- list()
  msgs <- character()
  value <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warns[[length(warns) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  for (w in warns) cond(case, class, arm, step, "warning", paste(class(w), collapse = ","), conditionMessage(w))
  for (m in msgs[grepl("ITSD|Identified|Removed|Kept|Skipping|unnormalized|normalization complete", msgs)]) {
    cond(case, class, arm, step, "message", "message", m)
  }
  lines <- LOGBUF$lines
  for (l in lines[grepl(LOGGER_KEEP, lines)]) cond(case, class, arm, step, "logger", sub(" .*$", "", l), l)
  if (inherits(value, "error")) cond(case, class, arm, step, "error", paste(class(value), collapse = ","), conditionMessage(value))
  list(value = value, failed = inherits(value, "error"), warnings = vapply(warns, conditionMessage, ""), logger = lines)
}

installLogger <- function() {
  appender <- function(lines) LOGBUF$lines <- c(LOGBUF$lines, lines)
  for (ns in c("global", "MultiScholaR")) {
    logger::log_appender(appender, namespace = ns)
    logger::log_layout(logger::layout_simple, namespace = ns)
    logger::log_threshold(logger::INFO, namespace = ns)
  }
}

# ---- synthetic objects (ticket section 4) ---------------------------------------------------------------------------
RUNS <- kit172Design()$Run
S <- "B_6"

#' Two assays of `n` features and `n_is` standards `IS_<assay>_<i>`; 2 groups x 6 samples (KIT-172's design); levels
#' about 20 log2; one loading per sample N(0, 0.5^2) on log2, carried by the features and the standards; noise N(0, 0.25^2).
g0Draw <- function(k, n = c(A1 = 40L, A2 = 40L), n_is = c(A1 = 5L, A2 = 5L)) {
  SEED(k)
  loading <- stats::rnorm(length(RUNS), 0, 0.5)
  out <- list()
  for (a in names(n)) {
    level <- stats::rnorm(n[[a]], 20, 1.5)
    noise <- matrix(stats::rnorm(n[[a]] * length(RUNS), 0, 0.25), n[[a]])
    y <- outer(level, rep(1, length(RUNS))) + outer(rep(1, n[[a]]), loading) + noise
    ids <- sprintf("%s_F%03d", a, seq_len(n[[a]]))
    if (n_is[[a]] > 0L) {
      s <- matrix(16 + stats::rnorm(n_is[[a]] * length(RUNS), 0, 0.1), n_is[[a]]) + outer(rep(1, n_is[[a]]), loading)
      y <- rbind(y, s)
      ids <- c(ids, sprintf("IS_%s_%d", a, seq_len(n_is[[a]])))
    }
    m <- 2^y
    dimnames(m) <- list(ids, RUNS)
    out[[a]] <- m
  }
  out
}

#' The annotation column holds the standards' own `IS_` names (as helper-module-ci-lipid-norm.R), unless given.
g0Table <- function(m, class, annotation = NULL) {
  ids <- rownames(m)
  ann <- if (is.null(annotation)) ifelse(startsWith(ids, "IS_"), ids, "PC") else annotation
  table <- if (identical(class, "lipid")) data.frame(LipidName = ids, LipidClass = ann, stringsAsFactors = FALSE) else
    data.frame(Metabolite = ids, Ontology = ann, stringsAsFactors = FALSE)
  cbind(table, stats::setNames(as.data.frame(unname(m)), colnames(m)))
}

g0Object <- function(values, class, annotations = list()) {
  assays <- lapply(stats::setNames(names(values), names(values)), function(a) g0Table(values[[a]], class, annotations[[a]]))
  design <- kit172Design()
  args <- list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
  if (identical(class, "lipid")) {
    methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = design, lipid_id_column = "LipidName",
      annotation_id_column = "LipidClass", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  } else {
    methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = design, metabolite_id_column = "Metabolite",
      annotation_id_column = "Ontology", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  }
}

istdAnnotation <- function(m) ifelse(startsWith(rownames(m), "IS_"), "ISTD", "PC")

# ---- reading states -------------------------------------------------------------------------------------------------
g0State <- function(owner, state) tryCatch(owner$getState(state), error = function(e) NULL)

#' An assay's sample-value matrix, rows by feature id (sorted), columns the design's runs present in the table.
g0Matrix <- function(object, binding, assay) {
  if (is.null(object)) return(NULL)
  table <- methods::slot(object, binding$data_slot)[[assay]]
  if (is.null(table)) return(NULL)
  table <- as.data.frame(table, check.names = FALSE, stringsAsFactors = FALSE)
  runs <- intersect(RUNS, names(table))
  m <- as.matrix(table[, runs, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- as.character(table[[methods::slot(object, binding$id_slot)]])
  m[order(rownames(m)), , drop = FALSE]
}

#' Equal (all.equal, tolerance 0) on the rows and columns both hold; with the dimensions.
g0Equal <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NA)
  rows <- intersect(rownames(a), rownames(b))
  cols <- intersect(colnames(a), colnames(b))
  isTRUE(all.equal(unname(a[rows, cols, drop = FALSE]), unname(b[rows, cols, drop = FALSE]), tolerance = 0))
}

#' Per assay: dimensions, equality with the input (non-`IS_` rows), finite values of sample S, IS_ rows.
g0Compare <- function(case, class, arm, before, after, binding, label) {
  if (is.null(after)) {
    obs(case, class, arm, paste(label, "state"), "absent")
    return(invisible())
  }
  for (assay in names(methods::slot(before, binding$data_slot))) {
    mb <- g0Matrix(before, binding, assay)
    ma <- g0Matrix(after, binding, assay)
    if (is.null(ma)) {
      obs(case, class, arm, sprintf("%s %s", label, assay), "assay absent from the returned object")
      next
    }
    feat_b <- mb[!startsWith(rownames(mb), "IS_"), , drop = FALSE]
    feat_a <- ma[!startsWith(rownames(ma), "IS_"), , drop = FALSE]
    obs(case, class, arm, sprintf("%s %s equal to input (non-IS rows; all.equal tol 0)", label, assay), g0Equal(feat_b, feat_a))
    obs(case, class, arm, sprintf("%s %s dims in -> out", label, assay), sprintf("%dx%d -> %dx%d", nrow(mb), ncol(mb), nrow(ma), ncol(ma)))
    obs(case, class, arm, sprintf("%s %s IS_ rows in -> out", label, assay),
      sprintf("%d -> %d", sum(startsWith(rownames(mb), "IS_")), sum(startsWith(rownames(ma), "IS_"))))
    obs(case, class, arm, sprintf("%s %s sample %s finite values in -> out", label, assay, S),
      sprintf("%s -> %s", if (S %in% colnames(mb)) sum(is.finite(mb[, S])) else "column absent",
        if (S %in% colnames(ma)) sum(is.finite(ma[, S])) else "column absent"))
  }
}

g0Design <- function(object) if (is.null(object)) "absent" else as.character(object@design_matrix[[object@sample_id]])

#' A state's design samples and, per assay, its columns and sample S's finite values.
g0StateSummary <- function(case, class, arm, owner, state, binding, label) {
  object <- g0State(owner, state)
  if (is.null(object)) {
    obs(case, class, arm, paste(label, "state"), "absent")
    return(invisible(NULL))
  }
  design <- g0Design(object)
  obs(case, class, arm, paste(label, "design samples"), sprintf("%d; %s in design: %s", length(design), S, S %in% design))
  for (assay in names(methods::slot(object, binding$data_slot))) {
    m <- g0Matrix(object, binding, assay)
    obs(case, class, arm, sprintf("%s %s", label, assay), sprintf("%d features x %d samples; %s column: %s", nrow(m), ncol(m), S,
      if (S %in% colnames(m)) sprintf("present, %d finite", sum(is.finite(m[, S]))) else "absent"))
  }
  obs(case, class, arm, paste(label, "assays"), names(methods::slot(object, binding$data_slot)))
  invisible(object)
}

g0Record <- function(owner) tryCatch(suppressWarnings(suppressMessages(buildWorkflowAnalysisRecord(owner))), error = function(e) e)

g0StepObs <- function(case, class, arm, record, state) {
  if (inherits(record, "error")) {
    obs(case, class, arm, paste("record", state), paste("record failed:", conditionMessage(record)))
    return(invisible())
  }
  step <- NULL
  for (s in record$steps) if (identical(s$state, state)) step <- s
  if (is.null(step)) {
    obs(case, class, arm, paste("record", state), "no step for the state")
    return(invisible())
  }
  params <- vapply(step$parameters, function(p) paste0(p$name, "=", tryCatch(oneLine(workflowRecordValueText(p$value)),
    error = function(e) txt(p$value))), "")
  obs(case, class, arm, paste("record", state, "operation"), step$operation)
  obs(case, class, arm, paste("record", state, "parameters"), paste(params, collapse = "; "))
  obs(case, class, arm, paste("record", state, "outcomes"), step$outcomes)
}

#' Find every element named `name` in a nested list.
findNamed <- function(x, name) {
  found <- list()
  walk <- function(y) {
    if (!is.list(y)) return(invisible())
    if (!is.null(names(y)) && name %in% names(y)) found[[length(found) + 1L]] <<- y[[name]]
    for (z in y) walk(z)
  }
  walk(x)
  found
}

g0Owner <- function(root, object, class) suppressWarnings(suppressMessages(kit172Owner(root, object, class)))

g0Pre <- function(binding, owner, config, offset = 0) {
  commitWorkflowOmicsStoreStep(owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"), config, "QC complete")
  runWorkflowOmicsReferenceLog2(binding, owner, offset, config)
}

DA_SETTINGS <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)

g0Da <- function(case, class, arm, binding, owner) {
  da <- capture(runWorkflowOmicsReferenceDaAnalysis(binding, owner, getWorkflowStateReference(owner), DA_SETTINGS, kit172Contrasts()),
    case, class, arm, "da")
  if (da$failed) {
    obs(case, class, arm, "DA", paste("failed:", conditionMessage(da$value)))
    return(invisible(da))
  }
  manifest <- da$value$manifest
  obs(case, class, arm, "DA source state", owner$current_state)
  for (contrast in names(manifest$results)) {
    entry <- manifest$results[[contrast]]
    obs(case, class, arm, paste("DA", contrast, "assays with a table"), names(entry$assays))
    obs(case, class, arm, paste("DA", contrast, "rows per assay"), vapply(entry$assays, function(e) as.integer(e$rows), 1L))
    obs(case, class, arm, paste("DA", contrast, "errors"), if (length(entry$errors)) entry$errors else "none")
    obs(case, class, arm, paste("DA", contrast, "qvalue_warnings"), if (length(entry$qvalue_warnings)) entry$qvalue_warnings else "none")
  }
  invisible(da)
}

closeQuietly <- function(owner) try(closeDIAStateAdapter(owner), silent = TRUE)

# ---- C1: an empty sample column, cyclic loess (declaration and worker arms) ----------------------------------------
objectC1 <- function(class) {
  v <- g0Draw(1L)
  v$A2[, S] <- NA_real_
  g0Object(v, class)
}

caseC1 <- function(root) {
  for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    object <- objectC1(class)
    local({
      fixture <- g0Owner(file.path(root, "c1", class, "decl"), object, class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      capture(g0Pre(binding, owner, config), "C1", class, "declaration", "qc_finalize+log2")
      before <- g0State(owner, paste0(prefix, "_log2"))
      run <- capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config), "C1", class,
        "declaration", "between_sample")
      after <- g0State(owner, paste0(prefix, "_normalized"))
      g0Compare("C1", class, "declaration", before, after, binding, "between-sample")
      obs("C1", class, "declaration", "kernel product applied", if (run$failed) "step failed" else run$value$value$applied)
      obs("C1", class, "declaration", "S4 args normalisation_method", if (is.null(after)) "absent" else after@args$normalisation_method)
      g0StepObs("C1", class, "declaration", g0Record(owner), paste0(prefix, "_normalized"))
      g0Da("C1", class, "declaration", binding, owner)
    })
    local({
      testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
      fixture <- g0Owner(file.path(root, "c1", class, "worker"), object, class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      capture(g0Pre(binding, owner, config), "C1", class, "worker", "qc_finalize+log2")
      handoff <- capture(omicsWorkerRun(fixture, paste0(prefix, "_between_sample_normalisation"), list(normalisation_method = "cyclicloess")),
        "C1", class, "worker", "omicsWorkerRun")
      h <- handoff$value
      if (handoff$failed) {
        obs("C1", class, "worker", "handoff", paste("omicsWorkerRun failed:", conditionMessage(h)))
      } else {
        strings <- unlist(rapply(h, function(x) x, classes = "character", how = "unlist"), use.names = FALSE)
        obs("C1", class, "worker", "handoff status", h$status)
        obs("C1", class, "worker", "handoff error", h$error %||% "none")
        obs("C1", class, "worker", "handoff result fields", names(h$result))
        obs("C1", class, "worker", "handoff 'applied' values found", findNamed(h, "applied"))
        obs("C1", class, "worker", "handoff holds the text 'Skipping normalization' or 'only NA'",
          any(grepl("Skipping normalization|only NA", strings)))
        obs("C1", class, "worker", "handoff holds any 'warning' field", length(findNamed(h, "warnings")) + length(findNamed(h, "warning")))
      }
    })
  }
}

# ---- C2: the limma method fails (quantile with one value; a planted cyclic loess error) ------------------------------
caseC2 <- function(root) {
  for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    v <- g0Draw(2L)
    one <- v$A2
    keep <- one[1L, S]
    one[, S] <- NA_real_
    one[1L, S] <- keep
    v1 <- v
    v1$A2 <- one
    for (setup in c("quantile_one_value", "planted_cyclicloess")) {
      object <- g0Object(if (identical(setup, "quantile_one_value")) v1 else v, class)
      method <- if (identical(setup, "quantile_one_value")) "quantile" else "cyclicloess"
      arm <- paste("declaration", setup)
      local({
        fixture <- g0Owner(file.path(root, "c2", class, setup), object, class)
        owner <- fixture$owner
        on.exit(closeQuietly(owner), add = TRUE)
        config <- fixture$workflow$config_list
        capture(g0Pre(binding, owner, config), "C2", class, arm, "qc_finalize+log2")
        before <- g0State(owner, paste0(prefix, "_log2"))
        if (identical(setup, "planted_cyclicloess")) {
          testthat::local_mocked_bindings(normalizeCyclicLoess = function(...) stop("planted"), .package = "limma")
        }
        run <- capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, method, config), "C2", class, arm, "between_sample")
        after <- g0State(owner, paste0(prefix, "_normalized"))
        g0Compare("C2", class, arm, before, after, binding, "between-sample")
        obs("C2", class, arm, "kernel product applied", if (run$failed) "step failed" else run$value$value$applied)
        if (identical(setup, "planted_cyclicloess")) obs("C2", class, arm, "instrument: a warning carries 'planted'", any(grepl("planted", run$warnings)))
        g0StepObs("C2", class, arm, g0Record(owner), paste0(prefix, "_normalized"))
      })
    }
  }
}

# ---- C3, C3m: ITSD matches nothing in an assay; manual ids for A1 only -------------------------------------------------
g0ItsdObs <- function(case, class, arm, owner, binding, run) {
  prefix <- binding$prefix
  before <- g0State(owner, paste0(prefix, "_qc_complete"))
  after <- g0State(owner, paste0(prefix, "_itsd_norm"))
  g0Compare(case, class, arm, before, after, binding, "ITSD")
  itsd <- if (is.null(after)) NULL else after@args$ITSDNormalization
  obs(case, class, arm, "args ITSDNormalization applied", itsd$applied)
  obs(case, class, arm, "args ITSDNormalization method_type", itsd$method_type)
  obs(case, class, arm, "args itsd_features_per_assay (assays)", names(itsd$itsd_features_per_assay) %||% "none")
  for (a in names(itsd$itsd_features_per_assay)) obs(case, class, arm, paste("args itsd_features_per_assay", a), itsd$itsd_features_per_assay[[a]])
  obs(case, class, arm, "args itsd_counts_per_assay", itsd$itsd_counts_per_assay %||% "absent (no such field at base)")
  g0StepObs(case, class, arm, g0Record(owner), paste0(prefix, "_itsd_norm"))
}

caseC3 <- function(root) {
  for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    v <- g0Draw(3L)
    setups <- list(A2_unmatched = list(A2 = istdAnnotation(v$A2)),
      both_unmatched = list(A1 = istdAnnotation(v$A1), A2 = istdAnnotation(v$A2)))
    for (setup in names(setups)) {
      object <- g0Object(v, class, setups[[setup]])
      arm <- paste("declaration", setup)
      local({
        fixture <- g0Owner(file.path(root, "c3", class, setup), object, class)
        owner <- fixture$owner
        on.exit(closeQuietly(owner), add = TRUE)
        config <- fixture$workflow$config_list
        commitWorkflowOmicsStoreStep(owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"), config, "QC complete")
        run <- capture(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config), "C3", class, arm, "itsd")
        g0ItsdObs("C3", class, arm, owner, binding, run)
      })
    }
    # C3m: standards matchable in both assays, manual ids for A1 only
    object <- g0Object(v, class)
    local({
      fixture <- g0Owner(file.path(root, "c3m", class), object, class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      commitWorkflowOmicsStoreStep(owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"), config, "QC complete")
      ids <- list(A1 = sprintf("IS_A1_%d", 1:5))
      run <- capture(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", ids, config), "C3m", class, "declaration", "itsd")
      g0ItsdObs("C3m", class, "declaration", owner, binding, run)
    })
  }
}

# ---- C4: a sample without standard signal in A2, then every later step ----------------------------------------------
caseC4 <- function(root) {
  for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    v <- g0Draw(4L)
    v$A2[startsWith(rownames(v$A2), "IS_"), S] <- NA_real_
    object <- g0Object(v, class)
    arm <- "declaration chain"
    local({
      fixture <- g0Owner(file.path(root, "c4", class), object, class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      steps <- list(
        qc_finalize = function() commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete"),
        itsd = function() runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config),
        log2 = function() runWorkflowOmicsReferenceLog2(binding, owner, 0, config),
        between_sample = function() runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config),
        ruv_legacy_automatic = function() runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config),
        correlation_0.5 = function() runWorkflowOmicsReferenceCorrelation(binding, owner, "group", 0.5, config))
      states <- c(qc_finalize = "_qc_complete", itsd = "_itsd_norm", log2 = "_log2", between_sample = "_normalized",
        ruv_legacy_automatic = "_ruv_corrected", correlation_0.5 = "_correlation_filtered")
      for (step in names(steps)) {
        run <- capture(steps[[step]](), "C4", class, arm, step)
        if (run$failed) {
          obs("C4", class, arm, paste(step, "result"), paste("failed:", conditionMessage(run$value)))
          break
        }
        obs("C4", class, arm, paste(step, "result"), "completed")
        g0StateSummary("C4", class, arm, owner, paste0(prefix, states[[step]]), binding, step)
        if (identical(step, "itsd")) g0ItsdObs("C4", class, arm, owner, binding, run)
        if (identical(step, "between_sample")) {
          g0Compare("C4", class, arm, g0State(owner, paste0(prefix, "_log2")), g0State(owner, paste0(prefix, "_normalized")), binding, "between-sample")
          g0StepObs("C4", class, arm, g0Record(owner), paste0(prefix, "_normalized"))
        }
        if (identical(step, "ruv_legacy_automatic")) {
          obs("C4", class, arm, "RUV best_k", run$value$value$best_k)
          obs("C4", class, arm, "RUV control_count", run$value$value$control_count)
          dims <- grep("Dimensions before cleaning", run$logger, value = TRUE)
          obs("C4", class, arm, "RUV logged dimensions before/after cleaning (:713-726)", if (length(dims)) dims else "none logged")
          g0Compare("C4", class, arm, g0State(owner, paste0(prefix, "_normalized")), g0State(owner, paste0(prefix, "_ruv_corrected")), binding, "RUV")
          g0StepObs("C4", class, arm, g0Record(owner), paste0(prefix, "_ruv_corrected"))
        }
      }
      if (identical(owner$current_state, paste0(prefix, "_correlation_filtered"))) g0Da("C4", class, arm, binding, owner)
    })
  }
}

# ---- C6, C11: KIT-145's guard on -------------------------------------------------------------------------------------
guardObs <- function(case, class, arm, run) {
  e <- run$value
  if (!run$failed) {
    obs(case, class, arm, "guard", "no error")
    return(invisible())
  }
  obs(case, class, arm, "error classes", class(e))
  obs(case, class, arm, "error fields", list(assay = e$assay, sample = e$sample, fewest = e$fewest, minimum = e$minimum))
  obs(case, class, arm, "warnings raised before the error", if (length(run$warnings)) run$warnings else "none")
}

caseC6C11 <- function(root) {
  for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    small <- g0Object(g0Draw(6L, n = c(A1 = 40L, A2 = 3L), n_is = c(A1 = 0L, A2 = 0L)), class)
    local({
      withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
      fixture <- g0Owner(file.path(root, "c6", class, "decl"), small, class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      capture(g0Pre(binding, owner, config), "C6", class, "declaration", "qc_finalize+log2")
      guardObs("C6", class, "declaration", capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config),
        "C6", class, "declaration", "between_sample"))
    })
    local({
      withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
      testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
      fixture <- g0Owner(file.path(root, "c6", class, "worker"), small, class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      capture(g0Pre(binding, owner, config), "C6", class, "worker", "qc_finalize+log2")
      handoff <- capture(omicsWorkerRun(fixture, paste0(prefix, "_between_sample_normalisation"), list(normalisation_method = "cyclicloess")),
        "C6", class, "worker", "omicsWorkerRun")
      h <- handoff$value
      obs("C6", class, "worker", "handoff status", if (handoff$failed) paste("omicsWorkerRun failed:", conditionMessage(h)) else h$status)
      if (!handoff$failed) obs("C6", class, "worker", "handoff error", h$error %||% "none")
    })
    local({
      withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
      fixture <- g0Owner(file.path(root, "c11", class), objectC1(class), class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      capture(g0Pre(binding, owner, config), "C11", class, "declaration", "qc_finalize+log2")
      guardObs("C11", class, "declaration", capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config),
        "C11", class, "declaration", "between_sample"))
    })
  }
}

# ---- C7: the q-value fallback and the record ------------------------------------------------------------------------
caseC7 <- function(root) {
  for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    v <- g0Draw(7L, n = c(A1 = 40L, A2 = 10L), n_is = c(A1 = 0L, A2 = 0L))
    v$A2[1:8, startsWith(RUNS, "B_")] <- NA_real_
    object <- g0Object(v, class)
    local({
      fixture <- g0Owner(file.path(root, "c7", class), object, class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      capture(g0Pre(binding, owner, config), "C7", class, "declaration", "qc_finalize+log2")
      g0Da("C7", class, "declaration", binding, owner)
      record <- g0Record(owner)
      if (inherits(record, "error")) {
        obs("C7", class, "declaration", "record", paste("failed:", conditionMessage(record)))
      } else {
        for (entry in record$analyses$differential_abundance$results) {
          obs("C7", class, "declaration", sprintf("record DA entry %s / %s fields", entry$contrast, entry$assay %||% "-"), names(entry))
          obs("C7", class, "declaration", sprintf("record DA entry %s / %s fdr_method", entry$contrast, entry$assay %||% "-"),
            entry$fdr_method %||% "absent")
        }
      }
    })
  }
}

# ---- C9: proteomics and transcriptomics parity ----------------------------------------------------------------------
caseC9 <- function(root) {
  SEED(9L)
  loading <- stats::rnorm(length(RUNS), 0, 0.5)
  y <- outer(stats::rnorm(60L, 20, 1.5), rep(1, length(RUNS))) + outer(rep(1, 60L), loading) +
    matrix(stats::rnorm(60L * length(RUNS), 0, 0.25), 60L)
  counts <- matrix(stats::rpois(100L * 6L, rep(exp(stats::runif(100L, 2, 8)), 6L)), 100L, dimnames = list(sprintf("G%03d", 1:100), sprintf("S%d", 1:6)))
  dimnames(y) <- list(sprintf("P%03d", 1:60), RUNS)
  protObject <- function(m) methods::new("ProteinQuantitativeData",
    protein_quant_table = data.frame(Protein.Ids = rownames(m), as.data.frame(m, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE),
    protein_id_column = "Protein.Ids", design_matrix = kit172Design(), sample_id = "Run", group_id = "group",
    technical_replicate_id = "tech_rep_group", args = list())
  empty <- y
  empty[, S] <- NA_real_
  one <- y
  one[, S] <- NA_real_
  one[1L, S] <- y[1L, S]
  setups <- list(list("empty column", empty, "cyclicloess"), list("empty column", empty, "quantile"), list("empty column", empty, "scale"),
    list("one value", one, "quantile"))
  for (s in setups) {
    arm <- paste("proteomics", s[[1L]], s[[3L]])
    input <- s[[2L]]
    run <- capture(normaliseBetweenSamples(protObject(input), normalisation_method = s[[3L]]), "C9", "prot", arm, "normaliseBetweenSamples")
    if (run$failed) {
      obs("C9", "prot", arm, "result", paste("stopped:", conditionMessage(run$value)))
      next
    }
    out <- run$value@protein_quant_table
    m <- as.matrix(out[, RUNS])
    rownames(m) <- out$Protein.Ids
    m <- m[rownames(input), , drop = FALSE]
    others <- setdiff(RUNS, S)
    obs("C9", "prot", arm, "result", "completed")
    obs("C9", "prot", arm, "S4 args normalisation_method", run$value@args$normalisation_method %||% "absent")
    obs("C9", "prot", arm, "finite values in the whole table (of 720)", sum(is.finite(m)))
    obs("C9", "prot", arm, "every value missing", all(!is.finite(m)))
    obs("C9", "prot", arm, "finite values in the other samples (of 660)", sum(is.finite(m[, others])))
    obs("C9", "prot", arm, "finite values in sample S", sum(is.finite(m[, S])))
    obs("C9", "prot", arm, "other samples changed (max abs diff)", if (any(is.finite(m[, others]))) max(abs(m[, others] - input[, others]), na.rm = TRUE) else NA)
  }
  group <- factor(rep(c("A", "B"), each = 3L))
  control <- capture(trxTmmLogCpm(counts, group), "C9", "trx", "TMM control", "trxTmmLogCpm")
  obs("C9", "trx", "TMM control", "result", if (control$failed) paste("stopped:", conditionMessage(control$value)) else "completed")
  local({
    testthat::local_mocked_bindings(calcNormFactors = function(...) stop("planted"), .package = "edgeR")
    planted <- capture(trxTmmLogCpm(counts, group), "C9", "trx", "TMM planted error", "trxTmmLogCpm")
    obs("C9", "trx", "TMM planted error", "result", if (planted$failed) paste("stopped:", conditionMessage(planted$value)) else "completed (instrument did not plant)")
  })
}

# ---- C12 (a), (b), (a-decl): legacy RUV-III-C cannot correct an assay -----------------------------------------------
caseC12 <- function(root) {
  for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    v <- g0Draw(12L)
    object <- g0Object(v, class)
    local({
      fixture <- g0Owner(file.path(root, "c12", class, "s4"), object, class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      capture(g0Pre(binding, owner, config), "C12", class, "S4", "qc_finalize+log2")
      capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config), "C12", class, "S4", "between_sample")
      input <- g0State(owner, paste0(prefix, "_normalized"))
      ids <- sprintf("A1_F%03d", 1:8)
      ctrls <- list(a_A1_only = list(A1 = ids, A2 = c("absent_1", "absent_2")), b_neither = list(A1 = c("absent_1", "absent_2"), A2 = c("absent_3", "absent_4")))
      for (setup in names(ctrls)) {
        arm <- paste("S4", setup)
        run <- capture(ruvIII_C_Varying(input, ruv_grouping_variable = "group", ruv_number_k = list(A1 = 1L, A2 = 1L), ctrl = ctrls[[setup]]),
          "C12", class, arm, "ruvIII_C_Varying")
        if (run$failed) {
          obs("C12", class, arm, "result", paste("stopped:", conditionMessage(run$value)))
          next
        }
        obs("C12", class, arm, "assays returned", names(methods::slot(run$value, binding$data_slot)))
        g0Compare("C12", class, arm, input, run$value, binding, "RUV")
        obs("C12", class, arm, "warning() conditions raised", length(run$warnings))
        obs("C12", class, arm, "logger WARN lines", grep("^WARN", run$logger, value = TRUE))
      }
    })
    # (a-decl): every A2 feature has one missing value (no complete feature in A2), legacy Automatic through the declaration
    va <- v
    for (j in seq_len(nrow(va$A2))) va$A2[j, RUNS[((j - 1L) %% length(RUNS)) + 1L]] <- NA_real_
    objectA <- g0Object(va, class)
    local({
      fixture <- g0Owner(file.path(root, "c12", class, "decl"), objectA, class)
      owner <- fixture$owner
      on.exit(closeQuietly(owner), add = TRUE)
      config <- fixture$workflow$config_list
      arm <- "declaration a-decl (A2 without a complete feature)"
      capture(g0Pre(binding, owner, config), "C12", class, arm, "qc_finalize+log2")
      capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config), "C12", class, arm, "between_sample")
      run <- capture(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config), "C12", class, arm, "ruv")
      if (run$failed) {
        obs("C12", class, arm, "result", paste("stopped:", conditionMessage(run$value)))
      } else {
        obs("C12", class, arm, "RUV best_k", run$value$value$best_k)
        obs("C12", class, arm, "RUV control_count", run$value$value$control_count)
        obs("C12", class, arm, "optimisation success per assay", lapply(run$value$optimisation, function(o) o$success %||% NA))
        g0Compare("C12", class, arm, g0State(owner, paste0(prefix, "_normalized")), g0State(owner, paste0(prefix, "_ruv_corrected")), binding, "RUV")
        g0StepObs("C12", class, arm, g0Record(owner), paste0(prefix, "_ruv_corrected"))
        g0Da("C12", class, arm, binding, owner)
      }
    })
  }
}

# ---- C14: MS-DIAL's default pattern on a synthetic metabolomics object ---------------------------------------------
MSDIAL_PATTERN <- "^IS_|_d[0-9]+$|ISTD|IS-|w/o MS2:"

objectC14 <- function() {
  v <- g0Draw(14L, n_is = c(A1 = 0L, A2 = 0L))
  assays <- lapply(stats::setNames(names(v), names(v)), function(a) {
    m <- v[[a]]
    ann <- sprintf("Metabolite %s %d", a, seq_len(nrow(m)))
    if (identical(a, "A1")) {
      ann[1:5] <- sprintf("w/o MS2: Feature %d", 1:5)
      ann[6:7] <- c("PC 15:0_18:1(d7)", "LPC 18:1(d7)")
    }
    table <- data.frame(`Alignment ID` = rownames(m), `Metabolite name` = ann, check.names = FALSE, stringsAsFactors = FALSE)
    cbind(table, stats::setNames(as.data.frame(unname(m)), colnames(m)))
  })
  methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = kit172Design(), metabolite_id_column = "Alignment ID",
    annotation_id_column = "Metabolite name", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
    database_identifier_type = "Unknown", internal_standard_regex = MSDIAL_PATTERN,
    args = list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0)))
}

caseC14 <- function(root) {
  class <- "metab"
  binding <- kit172Binding(class)
  prefix <- binding$prefix
  object <- objectC14()
  wo <- sprintf("A1_F%03d", 1:5)
  d7 <- sprintf("A1_F%03d", 6:7)
  local({
    fixture <- g0Owner(file.path(root, "c14", "itsd"), object, class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    arm <- "ITSD regex mode, median"
    commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
    before <- g0State(owner, paste0(prefix, "_qc_complete"))
    obs("C14", class, arm, "object regex / annotation column", c(before@internal_standard_regex, before@annotation_id_column))
    run <- capture(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config), "C14", class, arm, "itsd")
    after <- g0State(owner, paste0(prefix, "_itsd_norm"))
    if (run$failed || is.null(after)) {
      obs("C14", class, arm, "result", paste("failed:", if (run$failed) conditionMessage(run$value) else "no state"))
    } else {
      ma <- g0Matrix(after, binding, "A1")
      mb <- g0Matrix(before, binding, "A1")
      obs("C14", class, arm, "'w/o MS2:' rows present after ITSD", intersect(wo, rownames(ma)))
      obs("C14", class, arm, "'(d7)' rows present after ITSD", intersect(d7, rownames(ma)))
      obs("C14", class, arm, "A1 rows in -> out", sprintf("%d -> %d", nrow(mb), nrow(ma)))
      keep <- setdiff(rownames(ma), c(wo, d7))
      obs("C14", class, arm, "A1 other features changed by ITSD", !g0Equal(mb[keep, , drop = FALSE], ma[keep, , drop = FALSE]))
      itsd <- after@args$ITSDNormalization
      obs("C14", class, arm, "args itsd_features_per_assay A1", itsd$itsd_features_per_assay$A1 %||% "absent")
      obs("C14", class, arm, "args itsd_features_per_assay (assays)", names(itsd$itsd_features_per_assay) %||% "none")
      obs("C14", class, arm, "args itsd_counts_per_assay", itsd$itsd_counts_per_assay %||% "absent (no such field at base)")
      g0StepObs("C14", class, arm, g0Record(owner), paste0(prefix, "_itsd_norm"))
    }
  })
  local({
    fixture <- g0Owner(file.path(root, "c14", "integration"), object, class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    arm <- "integration input, ITSD unticked"
    capture(g0Pre(binding, owner, config), "C14", class, arm, "qc_finalize+log2")
    capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config), "C14", class, arm, "between_sample")
    src <- capture(withWorkflowStateScope(owner, function() omicsIntegrationInputSource(binding, owner)), "C14", class, arm, "omicsIntegrationInputSource")
    if (src$failed || is.null(src$value)) {
      obs("C14", class, arm, "integration source", if (src$failed) paste("failed:", conditionMessage(src$value)) else "NULL")
    } else {
      for (view in names(src$value$features)) {
        f <- src$value$features[[view]]
        obs("C14", class, arm, paste("features table columns", view), names(f))
        marked <- as.character(f$feature_id[f$is_internal_standard %in% TRUE])
        obs("C14", class, arm, paste("rows marked is_internal_standard", view), if (length(marked)) marked else "none")
      }
      obs("C14", class, arm, "manifest standards rule", findNamed(src$value$manifest, "internal_standards"))
    }
  })
}

# ---- the managed tab (C1, C3, C3m, C12 a-decl) ------------------------------------------------------------------------
tabRun <- function(case, class, object, root, inputs, selections = NULL, after = NULL) {
  omic <- if (identical(class, "lipid")) "lipidomics" else "metabolomics"
  binding <- kit172Binding(class)
  prefix <- binding$prefix
  fixture <- g0Owner(root, object, class)
  owner <- fixture$owner
  on.exit(closeQuietly(owner), add = TRUE)
  probe <- new.env()
  notes <- new.env()
  notes$updates <- list(); notes$notes <- list(); notes$modals <- integer(); notes$choices <- list()
  res <- new.env()
  arm <- "managed tab"
  run <- capture(local({
    guiEntryInstallBrowserDouble(notes)
    omicsManagedChainMocks(fixture, probe)
    utils::capture.output(shiny::testServer(omicsTabServer(fixture, omic, probe, TRUE), {
      session$userData$probe_owner <- owner
      base <- omicsTabInputs()
      base[names(inputs)] <- inputs
      do.call(session$setInputs, base)
      commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), list(),
        "QC processing complete - ready for normalization")
      if (!is.null(selections)) {
        session$userData$probe_norm$itsd_selections <- selections$rows
        session$userData$probe_norm$itsd_candidate_ids <- selections$candidates
      }
      session$setInputs(run_normalization = 1)
      for (i in seq_len(10L)) {
        session$elapse(200)
        session$flushReact()
      }
      res$normalization_complete <- shiny::isolate(session$userData$probe_norm$normalization_complete)
      res$ruv_complete <- shiny::isolate(session$userData$probe_norm$ruv_complete)
    }))
    NULL
  }), case, class, arm, "run_normalization")
  if (run$failed) obs(case, class, arm, "tab run", paste("failed:", conditionMessage(run$value)))
  for (i in seq_along(probe$log)) {
    OUT$log[[length(OUT$log) + 1L]] <- data.frame(case = case, class = class, arm = arm, index = i, line = oneLine(probe$log[[i]]),
      stringsAsFactors = FALSE)
  }
  for (n in notes$notes) {
    OUT$notes[[length(OUT$notes) + 1L]] <- data.frame(case = case, class = class, arm = arm, type = n$type, message = oneLine(n$message),
      stringsAsFactors = FALSE)
  }
  obs(case, class, arm, "notifications shown", length(notes$notes))
  obs(case, class, arm, "normalization_complete", res$normalization_complete %||% "unset")
  obs(case, class, arm, "ruv_complete", res$ruv_complete %||% "unset")
  obs(case, class, arm, "worker operations", probe$operations %||% "none")
  obs(case, class, arm, "history", unlist(owner$getHistory()))
  if (is.function(after)) after(owner, binding, arm)
}

caseTab <- function(root) {
  for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    tabRun("C1", class, objectC1(class), file.path(root, "c1", class),
      list(apply_itsd = FALSE, norm_method = "cyclicloess", ruv_mode = "skip", log_offset = 0),
      after = function(owner, binding, arm) {
        g0Compare("C1", class, arm, g0State(owner, paste0(prefix, "_log2")), g0State(owner, paste0(prefix, "_normalized")), binding, "between-sample")
        g0StepObs("C1", class, arm, g0Record(owner), paste0(prefix, "_normalized"))
      })
    v <- g0Draw(3L)
    tabRun("C3", class, g0Object(v, class, list(A1 = istdAnnotation(v$A1), A2 = istdAnnotation(v$A2))), file.path(root, "c3", class),
      list(apply_itsd = TRUE, itsd_aggregation = "median", norm_method = "cyclicloess", ruv_mode = "skip", log_offset = 0),
      after = function(owner, binding, arm) g0ItsdObs("C3", class, arm, owner, binding, NULL))
    tabRun("C3m", class, g0Object(v, class), file.path(root, "c3m", class),
      list(apply_itsd = TRUE, itsd_aggregation = "median", norm_method = "cyclicloess", ruv_mode = "skip", log_offset = 0),
      selections = list(rows = list(A1 = 1:5), candidates = list(A1 = sprintf("IS_A1_%d", 1:5))),
      after = function(owner, binding, arm) g0ItsdObs("C3m", class, arm, owner, binding, NULL))
    va <- g0Draw(12L)
    for (j in seq_len(nrow(va$A2))) va$A2[j, RUNS[((j - 1L) %% length(RUNS)) + 1L]] <- NA_real_
    tabRun("C12", class, g0Object(va, class), file.path(root, "c12", class),
      list(apply_itsd = FALSE, norm_method = "cyclicloess", ruv_mode = "automatic", log_offset = 0),
      after = function(owner, binding, arm) {
        g0Compare("C12", class, arm, g0State(owner, paste0(prefix, "_normalized")), g0State(owner, paste0(prefix, "_ruv_corrected")), binding, "RUV")
        g0StepObs("C12", class, arm, g0Record(owner), paste0(prefix, "_ruv_corrected"))
      })
  }
}

# ---- C8: C1's record saved and replayed -----------------------------------------------------------------------------
caseC8 <- function(root, tree) {
  harness <- new.env(parent = globalenv())
  sys.source(file.path(tree, "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  for (class in c("lipid", "metab")) {
    binding <- kit172Binding(class)
    local({
      fixture <- g0Owner(file.path(root, "c8", class), objectC1(class), class)
      owner <- fixture$owner
      config <- fixture$workflow$config_list
      arm <- "record and replay"
      capture(g0Pre(binding, owner, config), "C8", class, arm, "qc_finalize+log2")
      origin <- capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config), "C8", class, arm, "between_sample (origin)")
      obs("C8", class, arm, "origin: empty-column warnings", sum(grepl("only NA values", origin$warnings)))
      capture(runWorkflowOmicsReferenceDaAnalysis(binding, owner, getWorkflowStateReference(owner), DA_SETTINGS, kit172Contrasts()), "C8", class, arm, "da (origin)")
      json <- file.path(root, "c8", paste0(class, "-record.json"))
      capture(writeWorkflowAnalysisRecord(owner, json), "C8", class, arm, "write record")
      history <- unlist(owner$getHistory())
      closeQuietly(owner)
      source_dir <- fixture$source_dir
      record <- jsonlite::read_json(json, simplifyVector = FALSE)
      obs("C8", class, arm, "record valid", validateWorkflowAnalysisRecord(record)$valid)
      copy <- file.path(root, "c8", paste0(class, "-copy"))
      dir.create(copy, recursive = TRUE)
      stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
      store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
      copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
      operations <- vapply(record$steps, function(s) s$operation, "")
      states <- vapply(record$steps, function(s) s$state, "")
      pre <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
      obs("C8", class, arm, "replay preconditions", pre)
      files <- c(list.files(file.path(source_dir, "import_sources"), full.names = TRUE), file.path(source_dir, "design_matrix.tab"))
      sums <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "")
      inputs <- character()
      for (input in record$inputs) {
        hit <- files[sums == (input$sha256 %||% "")]
        if (length(hit)) inputs[[input$role]] <- hit[[1L]]
      }
      obs("C8", class, arm, "replay inputs (roles)", names(inputs))
      prepared <- capture(harness$kitPrepareReplayProject(copy_root, pre, api), "C8", class, arm, "prepare replay")
      if (prepared$failed) {
        obs("C8", class, arm, "replay", paste("prepare failed:", conditionMessage(prepared$value)))
        return(invisible())
      }
      scratch <- file.path(root, "c8", paste0(class, "-scratch"))
      dir.create(scratch, recursive = TRUE)
      replay <- capture(harness$kitReplayJourney(record, prepared$value, store_root, inputs, fixture$paths, scratch, api),
        "C8", class, arm, "replay")
      if (replay$failed) {
        obs("C8", class, arm, "replay", paste("failed:", conditionMessage(replay$value)))
      } else {
        obs("C8", class, arm, "replay: empty-column warnings", sum(grepl("only NA values", replay$warnings)))
        obs("C8", class, arm, "replay states equal history", identical(replay$value$states, history))
        obs("C8", class, arm, "replay mismatches", if (length(replay$value$mismatches)) as.character(replay$value$mismatches) else "none")
        obs("C8", class, arm, "replayed operations", vapply(replay$value$per_step, function(x) x$operation, ""))
      }
    })
  }
}

# ---- C5, C13: the GUI-entry lanes as committed -----------------------------------------------------------------------
laneCapture <- function(manager, short, omic) {
  binding <- kit172Binding(short)
  prefix <- binding$prefix
  out <- list()
  withWorkflowStateScope(manager, function() {
    history <- unlist(manager$.context$snapshot$state_history)
    out$history <<- history
    record <- g0Record(manager)
    for (state in paste0(prefix, c("_itsd_norm", "_normalized"))) {
      step <- NULL
      if (!inherits(record, "error")) for (s in record$steps) if (identical(s$state, state)) step <- s
      out[[paste("record", state, "outcomes")]] <<- if (is.null(step)) "no step" else step$outcomes
      out[[paste("record", state, "parameters")]] <<- if (is.null(step)) "no step" else
        vapply(step$parameters, function(p) paste0(p$name, "=", tryCatch(oneLine(workflowRecordValueText(p$value)), error = function(e) "?")), "")
    }
    qc <- g0State(manager, paste0(prefix, "_qc_complete"))
    itsd <- g0State(manager, paste0(prefix, "_itsd_norm"))
    if (!is.null(itsd)) {
      out$`args ITSDNormalization` <<- itsd@args$ITSDNormalization[intersect(c("applied", "method_type", "itsd_pattern_columns"),
        names(itsd@args$ITSDNormalization))]
      out$`args itsd_features_per_assay (assays)` <<- names(itsd@args$ITSDNormalization$itsd_features_per_assay) %||% "none"
    }
    if (!is.null(qc)) {
      regex <- qc@internal_standard_regex
      annotation <- qc@annotation_id_column
      id_col <- methods::slot(qc, binding$id_slot)
      out$`object regex` <<- regex
      out$`object annotation column` <<- annotation
      build <- omicGuiConfig(omic)$norm$itsd_table
      for (assay in names(methods::slot(qc, binding$data_slot))) {
        table <- as.data.frame(methods::slot(qc, binding$data_slot)[[assay]], check.names = FALSE)
        matched <- if (!is.na(annotation) && annotation %in% names(table)) sum(stringr::str_detect(as.character(table[[annotation]]), regex) %in% TRUE) else NA
        out[[paste("C13", assay, "rows matched by the object's regex in the annotation column")]] <<- matched
        candidates <- tryCatch({
          selection <- build(table, id_col, annotation)
          sprintf("%d of %d (%s)", sum(selection$is_candidate %in% TRUE), nrow(selection),
            paste(utils::head(as.character(selection$feature_id[selection$is_candidate %in% TRUE]), 8L), collapse = ", "))
        }, error = function(e) paste("builder failed:", conditionMessage(e)))
        out[[paste("C13", assay, "candidates flagged by the ITSD Selection builder")]] <<- candidates
      }
    }
  })
  base <- tryCatch(get("project_dirs", envir = globalenv())[[omic]]$base_dir, error = function(e) NULL)
  if (!is.null(base)) {
    files <- list.files(base, pattern = "\\.(txt|md)$", recursive = TRUE, full.names = TRUE)
    for (f in files) {
      lines <- readLines(f, warn = FALSE)
      hit <- grep("ITSD|Internal Standard|Normali[sz]ation", lines, ignore.case = TRUE)
      if (!length(hit)) next
      keep <- sort(unique(unlist(lapply(hit, function(i) i:min(length(lines), i + 3L)))))
      out[[paste("file", sub(paste0("^", base, "/?"), "", f))]] <- paste(lines[keep], collapse = " | ")
    }
  }
  out
}

caseLane <- function(lane_id, short, omic, record_dir) {
  res <- new.env()
  original <- closeDIAStateAdapter
  local({
    testthat::local_mocked_bindings(closeDIAStateAdapter = function(manager) {
      history <- tryCatch(unlist(manager$.context$snapshot$state_history), error = function(e) character())
      if (is.null(res$captured) && any(endsWith(history, "_normalized"))) {
        res$captured <- tryCatch(laneCapture(manager, short, omic), error = function(e) list(error = conditionMessage(e)))
      }
      original(manager)
    }, .package = "MultiScholaR")
    withr::local_envvar(MULTISCHOLAR_GUI_ENTRY_RECORD_DIR = record_dir)
    run <- capture(guiEntryRunLane(lane_id), "C5", short, lane_id, "guiEntryRunLane")
    res$run <- run
  })
  run <- res$run
  if (run$failed) obs("C5", short, lane_id, "lane", paste("failed:", conditionMessage(run$value)))
  for (name in names(res$captured)) {
    case <- if (startsWith(name, "C13 ") || name %in% c("object regex", "object annotation column")) "C13" else "C5"
    obs(case, short, lane_id, sub("^C13 ", "", name), res$captured[[name]])
  }
  if (!run$failed) {
    for (step in names(run$value$steps)) {
      outputs <- run$value$steps[[step]]$outputs
      for (name in names(outputs)) {
        value <- outputs[[name]]
        if (is.character(value) && any(grepl("ITSD|ormalization complete|ormalisation complete", value))) {
          obs("C5", short, lane_id, paste("output at step", step, name), substr(paste(value, collapse = " "), 1L, 1500L))
        }
      }
      for (note in run$value$steps[[step]]$notifications) {
        OUT$notes[[length(OUT$notes) + 1L]] <- data.frame(case = "C5", class = short, arm = paste(lane_id, step), type = note$type,
          message = oneLine(note$message), stringsAsFactors = FALSE)
      }
    }
    obs("C5", short, lane_id, "lane seconds", run$value$seconds)
  }
}

# ---- the group runner -------------------------------------------------------------------------------------------------
runCase <- function(name, fn) {
  cat(sprintf("case %s ...\n", name))
  t0 <- proc.time()[[3L]]
  outcome <- tryCatch({ fn(); "ok" }, error = function(e) paste("case failed:", oneLine(conditionMessage(e))))
  if (!identical(outcome, "ok")) obs(name, "-", "-", "case", outcome)
  cat(sprintf("case %s: %s in %.1f s\n", name, outcome, proc.time()[[3L]] - t0))
}

runGroup <- function(group, tree, out_dir) {
  installLogger()
  root <- file.path(tempdir(), paste0("g0-", group))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (identical(group, "core")) {
    runCase("C1", function() caseC1(root))
    runCase("C2", function() caseC2(root))
    runCase("C3+C3m", function() caseC3(root))
    runCase("C4", function() caseC4(root))
    runCase("C6+C11", function() caseC6C11(root))
    runCase("C7", function() caseC7(root))
    runCase("C9", function() caseC9(root))
    runCase("C12", function() caseC12(root))
    runCase("C14", function() caseC14(root))
  }
  if (identical(group, "tab")) runCase("tab", function() caseTab(root))
  if (identical(group, "replay")) runCase("C8", function() caseC8(root, tree))
  if (identical(group, "lanes")) {
    record_dir <- file.path(out_dir, "lanes")
    dir.create(record_dir, recursive = TRUE, showWarnings = FALSE)
    runCase("C5 lipid_canonical", function() caseLane("lipid_canonical", "lipid", "lipidomics", record_dir))
    runCase("C5 metab_lc", function() caseLane("metab_lc", "metab", "metabolomics", record_dir))
  }
  bind <- function(x) if (length(x)) do.call(rbind, x)
  list(obs = bind(OUT$obs), cond = bind(OUT$cond), log = bind(OUT$log), notes = bind(OUT$notes))
}
})

# ---- KIT-188's journey (C12 c), loaded as sim-k188.R loads the package ----------------------------------------------
k188Group <- function(tree, sim_file) {
  suppressPackageStartupMessages(library(MultiScholaR))
  attach(asNamespace("MultiScholaR"), name = "MultiScholaR:namespace", warn.conflicts = FALSE)
  for (helper in c("helper-workflow-omics-owner.R", "helper-workflow-lipid-declarations.R", "helper-workflow-metab-declarations.R",
      "helper-kit172.R")) {
    sys.source(file.path(tree, "tests", "testthat", helper), envir = globalenv())
  }
  invisible(RUVIIIC::omp_set_num_threads(1L))
  options(MultiScholaR.ruv_design_panel = FALSE)
  constants <- c("GROUPS", "SCENARIOS", "MED_VAR_L", "MED_VAR_P", "CELLS", "MASKS", "REFUSALS", "BLIND", "warning_count",
    "EXPECTED_DEFAULTS", "ARMS_OF", "ROW_COLUMNS")
  for (e in parse(sim_file, keep.source = FALSE)) {
    if (!is.call(e) || !identical(e[[1L]], as.name("<-"))) next
    target <- e[[2L]]
    rhs <- e[[3L]]
    is_function <- is.call(rhs) && identical(rhs[[1L]], as.name("function"))
    named <- is.name(target) && as.character(target) %in% constants
    member <- is.call(target) && identical(target[[1L]], as.name("$")) && as.character(target[[2L]]) %in% c("warning_count", "EXPECTED_DEFAULTS")
    if ((is.name(target) && is_function) || named || member) eval(e, globalenv())
  }
  lines <- new.env()
  lines$all <- character()
  appender <- function(x) lines$all <- c(lines$all, x)
  for (ns in c("global", "MultiScholaR")) {
    logger::log_appender(appender, namespace = ns)
    logger::log_layout(logger::layout_simple, namespace = ns)
    logger::log_threshold(logger::INFO, namespace = ns)
  }
  d <- tryCatch(routeDefaults("lipid"), error = function(e) EXPECTED_DEFAULTS$lipid)
  sim <- generate("L12", "N1", 401L)
  values <- maskedValues(sim, "vd20")
  j <- journeyL(sim, values, "LA", REQUESTS_L(d)$LA, d, tempfile("k188-g0-"))
  row <- function(item, value) data.frame(case = "C12", class = "lipid", arm = "c KIT-188 L12 N1 401 vd20 LA", item = item,
    value = gsub("[\t\r\n]+", " ", paste(format(value), collapse = ";")), stringsAsFactors = FALSE)
  obs <- list(row("journey status (KIT-188 rule)", j$status), row("journey class", j$class %||% NA), row("journey step", j$step %||% NA),
    row("missing share", mean(is.na(values[sim$feat_ids, ]))),
    row("RUV best_k", unlist(j$ruv_value$best_k)), row("RUV best_k assays", names(j$ruv_value$best_k)),
    row("RUV control_count", unlist(j$ruv_value$control_count)),
    row("optimisation success per assay", vapply(j$optimisation, function(o) as.character(o$success %||% NA), "")),
    row("optimisation error per assay", vapply(j$optimisation, function(o) gsub("[\r\n\t]+", " ", as.character(o$error %||% "none"))[[1L]], "")),
    row("assays into RUV", names(j$m_in)), row("assays out of RUV", names(j$m_out) %||% "none"),
    row("assays equal in -> out", vapply(names(j$m_in), function(a) identical(j$m_in[[a]], j$m_out[[a]]), logical(1L))),
    row("logger lines total", length(lines$all)))
  keep <- grep("^(WARN|ERROR)|Skipping|No assays were|Returning original", lines$all, value = TRUE)
  cond <- if (length(keep)) data.frame(case = "C12", class = "lipid", arm = "c KIT-188 L12 N1 401 vd20 LA", step = "journeyL",
    kind = "logger", cls = sub(" .*$", "", keep), text = gsub("[\t\r\n]+", " ", keep), stringsAsFactors = FALSE)
  list(obs = do.call(rbind, obs), cond = cond)
}

# ---- child -----------------------------------------------------------------------------------------------------------
if (length(args) >= 1L && identical(args[[1L]], "--group")) {
  stopifnot(length(args) == 5L)
  group <- args[[2L]]
  tree <- normalizePath(args[[3L]])
  result <- tryCatch({
    if (identical(group, "k188")) {
      k188Group(tree, args[[5L]])
    } else {
      suppressMessages(pkgload::load_all(tree, export_all = TRUE, helpers = TRUE, quiet = TRUE))
      env <- new.env(parent = asNamespace("MultiScholaR"))
      eval(caseDefinitions, env)
      env$runGroup(group, tree, dirname(args[[4L]]))
    }
  }, error = function(e) list(error = gsub("[\t\r\n]+", " ", conditionMessage(e))))
  saveRDS(result, args[[4L]])
  quit(save = "no", status = 0L)
}

# ---- the job (parent) ------------------------------------------------------------------------------------------------
stopifnot(length(args) %in% 1:2)
out <- args[[1L]]
wanted <- if (length(args) == 2L) strsplit(args[[2L]], ",", fixed = TRUE)[[1L]] else c("core", "tab", "replay", "lanes", "k188")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
git <- function(...) suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
provenance <- c(sprintf("head\t%s", git("rev-parse", "HEAD")[[1L]]), sprintf("R\t%s", R.version.string),
  sprintf("started\t%s", format(Sys.time(), tz = "UTC")))
k188_ref <- NA_character_
for (ref in K188_REFS) {
  hit <- git("rev-parse", "--verify", "--quiet", paste0(ref, "^{commit}"))
  if (length(hit) && grepl("^[0-9a-f]{40}$", hit[[1L]])) { k188_ref <- hit[[1L]]; break }
}
k188_sim <- file.path(tempdir(), "sim-k188.R")
if (!is.na(k188_ref)) system2("git", c("show", paste0(k188_ref, ":", K188_SIM)), stdout = k188_sim)
provenance <- c(provenance, sprintf("k188_commit\t%s", k188_ref),
  sprintf("k188_sim_sha256\t%s", if (file.exists(k188_sim)) digest::digest(file = k188_sim, algo = "sha256") else NA))
minutes <- c(core = 30, tab = 15, replay = 12, lanes = 30, k188 = 12)
rscript <- file.path(R.home("bin"), "Rscript")
status <- list()
collected <- list()
for (g in intersect(names(minutes), wanted)) {
  rds <- file.path(out, paste0("child-", g, ".rds"))
  log <- file.path(out, paste0("child-", g, ".txt"))
  t0 <- proc.time()[[3L]]
  if (identical(g, "k188") && !file.exists(k188_sim)) {
    status[[g]] <- data.frame(group = g, exit = NA_integer_, seconds = 0, error = "sim-k188.R absent", stringsAsFactors = FALSE)
    next
  }
  exit <- system2("timeout", c("--kill-after=30", paste0(minutes[[g]], "m"), shQuote(rscript), shQuote(script), "--group", g,
    shQuote(getwd()), shQuote(rds), shQuote(k188_sim)), stdout = log, stderr = log)
  result <- if (file.exists(rds)) readRDS(rds) else list(error = "no result file (timeout or crash)")
  if (file.exists(rds)) unlink(rds)
  status[[g]] <- data.frame(group = g, exit = exit, seconds = round(proc.time()[[3L]] - t0, 1),
    error = if (is.null(result$error)) "" else result$error, stringsAsFactors = FALSE)
  if (is.null(result$error)) collected[[g]] <- result
  cat(sprintf("group %s: exit %s, %.0f s, %s\n", g, exit, proc.time()[[3L]] - t0, if (is.null(result$error)) "ok" else result$error))
}
bindPart <- function(part) {
  frames <- Filter(Negate(is.null), lapply(collected, `[[`, part))
  if (!length(frames)) return(NULL)
  do.call(rbind, frames)
}
writeTsv <- function(x, name) if (!is.null(x)) utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeTsv(bindPart("obs"), "observations.tsv")
writeTsv(bindPart("cond"), "conditions.tsv")
writeTsv(bindPart("log"), "tab-log.tsv")
writeTsv(bindPart("notes"), "notifications.tsv")
writeTsv(do.call(rbind, status), "groups.tsv")
writeLines(c(provenance, sprintf("finished\t%s", format(Sys.time(), tz = "UTC"))), file.path(out, "provenance.txt"))
cat(readLines(file.path(out, "groups.tsv")), sep = "\n")
