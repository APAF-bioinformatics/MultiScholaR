# KIT-260 gate G0 (ticket r1.2, section 6.1): observe, at the base commit, what the internal-standard pattern proposed
# for MS-DIAL files matches and what each consumer then does. No package code: branch kit-260-g0 carries only this
# script and its records and is never merged. The predictions (PREDICTIONS below, from the ticket's 6.1 table) are
# written to predictions.tsv before any case runs; nothing here judges a case except D0's per-pattern count check
# (observed against predicted, written beside each other). C14 is not repeated (KIT-242 run 45584784285).
#
# Usage (run-tests scripts lane, from the package root):
#   Rscript tools/kit260/g0.R <out dir>                        the job: each group in its own child process
#   Rscript tools/kit260/g0.R --group <group> <tree> <rds>     one group (child)
#
# Groups: d0 (D0, D0b), core (D1, D2, D2n, D3, D4, D5), d6 (D6), d7 (D7).
# Seeds: case k draws with set.seed(260000 + k, "Mersenne-Twister", "Inversion", "Rejection") before any package call:
#   D1, D4, D5 260001 (one object); D2, D2n 260002; D3 260003; D6 260006; D7 260007. D0 and D0b draw nothing.

args <- commandArgs(trailingOnly = TRUE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))

# ---- predictions (ticket r1.2, 6.1), written before any value is read -------------------------------------------------
BASE_MSDIAL <- "^IS_|_d[0-9]+$|ISTD|IS-|w/o MS2:"
OPTIONS <- c(
  base = BASE_MSDIAL,
  a = "^IS_|_d[0-9]+$|ISTD|IS-",
  b = "^IS_|_d[0-9]+$|ISTD|^IS-",
  c = "^IS_|_d[0-9]+$|ISTD|^IS-|\\(d[0-9]+\\)",
  c_hyphen = "^IS_|[_-]d[0-9]+$|ISTD|^IS-|\\(d[0-9]+\\)",
  e_on_c = "^((w/o|no) MS2: ?|low score: )?IS[_-]|_d[0-9]+$|ISTD|\\(d[0-9]+\\)")
# Section 4's named rows (20; r1.0's 14 plus r1.1's 6), with their role.
NAMED <- data.frame(stringsAsFactors = FALSE,
  annotation = c("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine",
    "w/o MS2: Alanine", "w/o MS2: Serine", "w/o MS2: Valine", "w/o MS2: Citric acid", "w/o MS2: Lactic acid",
    "CIS-ACONITIC ACID", "HIS-LEU",
    "cis-Aconitic acid", "PC 15:0_18:1(d7)", "Unknown",
    "w/o MS2:Glucose", "no MS2: Glucose", "low score: Glucose", "w/o MS2:IS_Caffeine", "Caffeine-d9", NA),
  role = c(rep("intended", 4), rep("wo_space", 5), rep("capitalised", 2), rep("control", 3),
    "wo_nospace", "ms5_no", "ms5_low", "wo_prefixed_standard", "hyphen_label", "missing_annotation"),
  r10 = c(rep(TRUE, 14), rep(FALSE, 6)))
pick <- function(...) { x <- c(...); sort(x) }
PRED_SETS <- list(
  base = pick("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine", NAMED$annotation[5:9], "w/o MS2:Glucose",
    "w/o MS2:IS_Caffeine", "CIS-ACONITIC ACID", "HIS-LEU"),
  a = pick("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine", "CIS-ACONITIC ACID", "HIS-LEU"),
  b = pick("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine"),
  c = pick("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine", "PC 15:0_18:1(d7)"),
  c_hyphen = pick("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine", "PC 15:0_18:1(d7)", "Caffeine-d9"),
  e_on_c = pick("IS_Caffeine", "Caffeine_d9", "Caffeine-ISTD", "IS-Leucine", "PC 15:0_18:1(d7)", "w/o MS2:IS_Caffeine"),
  # The other formats "match only the intended labels they name" (ticket D0); custom has no pattern (NA).
  progenesis = pick("IS_Caffeine", "Caffeine-ISTD"), xcms = pick("IS_Caffeine", "Caffeine-ISTD"),
  compound_discoverer = pick("IS_Caffeine", "Caffeine-ISTD"), lipidsearch = pick("IS_Caffeine", "Caffeine-ISTD"),
  custom = character())
# Per alternative of base's MS-DIAL default (the ticket's 2.4 table): counts on the 20 names.
PRED_BASE_ALT <- c(`^IS_` = 1L, `_d[0-9]+$` = 1L, ISTD = 1L, `IS-` = 3L, `w/o MS2:` = 7L)
LIPIDSEARCH <- "^IS_|ISTD|Internal|_IS"

PREDICTIONS <- data.frame(stringsAsFactors = FALSE, case = character(), item = character(), predicted = character())
pred <- function(case, item, predicted) PREDICTIONS[nrow(PREDICTIONS) + 1L, ] <<- c(case, item, paste(predicted, collapse = "; "))
for (p in names(PRED_SETS)) pred("D0", sprintf("pattern %s: names matched (TRUE) of the 20", p), sprintf("%d: %s", length(PRED_SETS[[p]]), paste(PRED_SETS[[p]], collapse = " | ")))
pred("D0", "every pattern: the missing annotation", "NA (not TRUE)")
for (a in names(PRED_BASE_ALT)) pred("D0", sprintf("base MS-DIAL alternative '%s': count", a), PRED_BASE_ALT[[a]])
pred("D0", "lipid and metab MS-DIAL, progenesis, xcms, compound_discoverer, custom defaults: identical strings per format", "TRUE (both readers hold the same strings)")
pred("D0b", "metab, every layout (msdial, custom, malformed, progenesis, xcms, compound_discoverer): updates$is_pattern", BASE_MSDIAL)
pred("D0b", "lipid, every layout except LipidSearch (msdial, custom, malformed, progenesis, xcms, compound_discoverer): updates$is_pattern", BASE_MSDIAL)
pred("D0b", "lipid, LipidSearch layout: detected_format / updates$is_pattern", sprintf("lipidsearch / %s", LIPIDSEARCH))
pred("D1", "itsd_counts_per_assay (declaration and worker)", "13 per assay")
pred("D1", "rows absent after ITSD", "the 13 matched rows and the missing-annotation row")
pred("D1", "rows present after ITSD", "PC 15:0_18:1(d7), Caffeine-d9, no MS2: Glucose, low score: Glucose, cis-Aconitic acid, Unknown, the 40 ordinary")
pred("D1", "console line 'Identified %d ITSD features via regex'", "prints NA (INFERRED); 'Using 13 ITSD features' after it")
pred("D1", "study-parameters file, ITSD Features Per Assay", "13 identifiers per assay (integers), first 5 shown, '... +8 more'")
pred("D1", "record", "itsd_applied TRUE; no names in the record")
pred("D2", "lipid twin, annotation 'Name'", "as D1")
pred("D2n", "lipid, name column 'Metabolite name'", "annotation NA (the preview proposes none); ITSD skips both assays; the integration input marks nothing")
pred("D3", "LipidSearch-built lipid object", "normalisation matches nothing (class column); the integration input marks 'IS_PC 15:0' (id column) in each view")
pred("D4", "D1's object, ITSD unticked: is_internal_standard", "TRUE on the 13 rows per view, FALSE on the missing-annotation row; the model's drop (integrationInternalStandards) 13 per view")
pred("D5", "QC tab analysis, box empty", "the 13 rows per assay reported as internal standards (total over the two assays 26)")
pred("D6", "complete version, A against B, method none (RUV off)", "every logFC differs by one common constant = minus the group contrast of log2(a_A / a_B), to 1e-9")
pred("D6", "complete version, A against B, cyclicloess and quantile (RUV off)", "identical to 1e-9")
pred("D6", "complete version, A against B, scale (RUV off)", "nonzero differences")
pred("D6", "missing-value version, every RUV-on row, B against C", "not predicted (recorded)")
pred("D7", "where the pattern lives; what replay reads", "column_mapping.json and the object slot hold the typed pattern; the design object takes the stored pattern; replay reads the stored state (INFERRED), not today's default")

# ---- the case definitions, evaluated in an environment whose parent is the package namespace (as a test file's) -----
caseDefinitions <- quote({

oneLine <- function(x) gsub("[\t\r\n]+", " ", paste(as.character(x), collapse = " "))
SEED <- function(k) set.seed(260000L + k, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
OUT <- new.env()
OUT$obs <- list(); OUT$cond <- list(); OUT$tables <- list()
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
  oneLine(paste(ifelse(is.na(x), "NA", format(x)), collapse = ";"))
}
obs <- function(case, class, arm, item, value) {
  OUT$obs[[length(OUT$obs) + 1L]] <- data.frame(case = case, class = class, arm = arm, item = item, value = txt(value),
    stringsAsFactors = FALSE)
}
cond <- function(case, class, arm, step, kind, cls, text) {
  OUT$cond[[length(OUT$cond) + 1L]] <- data.frame(case = case, class = class, arm = arm, step = step, kind = kind, cls = cls,
    text = oneLine(text), stringsAsFactors = FALSE)
}
addTable <- function(name, frame) OUT$tables[[name]] <- rbind(OUT$tables[[name]], frame)
LOGGER_KEEP <- "^(WARN|ERROR)|Skipping|Found IS|Analyzing internal|RUV|Invalid K|control features|No assays were|Returning original"

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
  for (m in msgs[grepl("ITSD|Identified|Removed|Kept|Skipping|unnormalized|normalization complete|annotation column", msgs)]) {
    cond(case, class, arm, step, "message", "message", m)
  }
  lines <- LOGBUF$lines
  for (l in lines[grepl(LOGGER_KEEP, lines)]) cond(case, class, arm, step, "logger", sub(" .*$", "", l), l)
  if (inherits(value, "error")) cond(case, class, arm, step, "error", paste(class(value), collapse = ","), conditionMessage(value))
  list(value = value, failed = inherits(value, "error"), warnings = vapply(warns, conditionMessage, ""), messages = msgs, logger = lines)
}

installLogger <- function() {
  appender <- function(lines) LOGBUF$lines <- c(LOGBUF$lines, lines)
  for (ns in c("global", "MultiScholaR")) {
    logger::log_appender(appender, namespace = ns)
    logger::log_layout(logger::layout_simple, namespace = ns)
    logger::log_threshold(logger::INFO, namespace = ns)
  }
}

closeQuietly <- function(owner) try(closeDIAStateAdapter(owner), silent = TRUE)
g0Owner <- function(root, object, class) suppressWarnings(suppressMessages(kit172Owner(root, object, class)))
g0State <- function(owner, state) tryCatch(owner$getState(state), error = function(e) NULL)
g0Record <- function(owner) tryCatch(suppressWarnings(suppressMessages(buildWorkflowAnalysisRecord(owner))), error = function(e) e)
RUNS <- kit172Design()$Run
DA_SETTINGS <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)

# ---- synthetic objects (ticket section 4) ---------------------------------------------------------------------------
#' Log2 values of two assays, `n` rows each: row levels N(20, 1.5^2), one loading per sample N(0, 0.5^2) shared by the
#' assays and every row, noise N(0, 0.25^2). Draws only; the seed is set by the caller's case.
drawLogs <- function(n) {
  loading <- stats::rnorm(length(RUNS), 0, 0.5)
  out <- list()
  for (a in c("A1", "A2")) {
    level <- stats::rnorm(n, 20, 1.5)
    noise <- matrix(stats::rnorm(n * length(RUNS), 0, 0.25), n)
    out[[a]] <- outer(level, rep(1, length(RUNS))) + outer(rep(1, n), loading) + noise
  }
  out
}

#' Annotations per assay: 40 ordinary ("Metabolite 001" ...) then the named rows (all 20, or r1.0's 14).
annotations <- function(full = TRUE) c(sprintf("Metabolite %03d", 1:40), NAMED$annotation[if (full) TRUE else NAMED$r10])

#' An assay table: the id column (integers 1..n), the name column, the samples (2^y).
assayTable <- function(y, ann, id_col, name_col, ids = seq_len(nrow(y))) {
  table <- data.frame(ids, ann, check.names = FALSE, stringsAsFactors = FALSE)
  names(table) <- c(id_col, name_col)
  cbind(table, stats::setNames(as.data.frame(2^unname(y)), RUNS))
}

ARGS <- list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
metabObject <- function(tables, pattern, id_col = "Alignment ID", ann_col = "Metabolite name") {
  methods::new("MetaboliteAssayData", metabolite_data = tables, design_matrix = kit172Design(), metabolite_id_column = id_col,
    annotation_id_column = ann_col, sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
    database_identifier_type = "Unknown", internal_standard_regex = pattern, args = ARGS)
}
lipidObject <- function(tables, pattern, id_col, ann_col) {
  methods::new("LipidomicsAssayData", lipid_data = tables, design_matrix = kit172Design(), lipid_id_column = id_col,
    annotation_id_column = ann_col, sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
    database_identifier_type = "Unknown", internal_standard_regex = pattern, args = ARGS)
}

#' D1's metab object (full list, 60 rows per assay); the pattern from the reader's default at this commit, not typed.
objectD1 <- function() {
  SEED(1L)
  y <- drawLogs(60L)
  pattern <- getMetabolomicsColumnDefaults("msdial")$is_pattern
  metabObject(lapply(y, assayTable, ann = annotations(), id_col = "Alignment ID", name_col = "Metabolite name"), pattern)
}
#' D2 / D2n's lipid object: the same layout, name column `Name` (D2) or `Metabolite name` with no annotation (D2n).
objectD2 <- function(name_col, annotation_column) {
  SEED(2L)
  y <- drawLogs(60L)
  pattern <- getLipidomicsColumnDefaults("msdial")$is_pattern
  lipidObject(lapply(y, assayTable, ann = annotations(), id_col = "Alignment ID", name_col = name_col), pattern,
    "Alignment ID", annotation_column)
}

# ---- reading states -------------------------------------------------------------------------------------------------
idsOf <- function(object, binding, assay) as.character(methods::slot(object, binding$data_slot)[[assay]][[methods::slot(object, binding$id_slot)]])
annOf <- function(object, binding, assay, column) as.character(methods::slot(object, binding$data_slot)[[assay]][[column]])

#' Per assay: which named rows (by their annotation) are present after a step, by id.
namedPresence <- function(case, class, arm, before, after, binding, name_col) {
  for (assay in names(methods::slot(before, binding$data_slot))) {
    ids_b <- idsOf(before, binding, assay)
    ann_b <- annOf(before, binding, assay, name_col)
    ids_a <- if (is.null(after)) character() else idsOf(after, binding, assay)
    present <- ids_b %in% ids_a
    obs(case, class, arm, sprintf("%s rows in -> out", assay), sprintf("%d -> %d", length(ids_b), length(ids_a)))
    named <- !grepl("^Metabolite [0-9]{3}$", ann_b) | is.na(ann_b)
    obs(case, class, arm, sprintf("%s named rows absent after", assay), ifelse(is.na(ann_b[named & !present]), "<NA annotation>", ann_b[named & !present]))
    obs(case, class, arm, sprintf("%s named rows present after", assay), ifelse(is.na(ann_b[named & present]), "<NA annotation>", ann_b[named & present]))
    obs(case, class, arm, sprintf("%s ordinary rows present after (of 40)", assay), sum(present & !named))
  }
}

itsdArgs <- function(case, class, arm, object) {
  itsd <- if (is.null(object)) NULL else object@args$ITSDNormalization
  obs(case, class, arm, "args ITSDNormalization applied", itsd$applied %||% "absent")
  obs(case, class, arm, "args itsd_counts_per_assay", itsd$itsd_counts_per_assay %||% "absent")
  obs(case, class, arm, "args itsd_features_per_assay (assays)", names(itsd$itsd_features_per_assay) %||% "none")
  for (a in names(itsd$itsd_features_per_assay)) obs(case, class, arm, paste("args itsd_features_per_assay", a), itsd$itsd_features_per_assay[[a]])
}

stepObs <- function(case, class, arm, record, state) {
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
  text <- tryCatch(as.character(jsonlite::toJSON(record, auto_unbox = TRUE, force = TRUE, null = "null", na = "string")), error = function(e) "")
  obs(case, class, arm, "record JSON contains 'w/o MS2' / 'CIS-ACONITIC'", c(grepl("w/o MS2", text, fixed = TRUE), grepl("CIS-ACONITIC", text, fixed = TRUE)))
}

#' The study-parameters file's ITSD block for an object (createWorkflowArgsFromConfig, no workflow data).
studyParameters <- function(case, class, arm, object, root) {
  dir <- file.path(root, paste0("sp-", case, "-", class))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  run <- capture(utils::capture.output(createWorkflowArgsFromConfig("kit260_g0", source_dir_path = dir, final_s4_object = object)),
    case, class, arm, "study parameters")
  file <- file.path(dir, "study_parameters.txt")
  if (!file.exists(file)) {
    obs(case, class, arm, "study_parameters.txt", if (run$failed) paste("not written:", conditionMessage(run$value)) else "not written")
    return(invisible())
  }
  lines <- readLines(file, warn = FALSE)
  start <- grep("^ITSD Normali[sz]ation:", lines)
  if (!length(start)) {
    obs(case, class, arm, "study_parameters.txt ITSD block", "absent")
    return(invisible())
  }
  end <- start[[1L]] + which(!nzchar(trimws(lines[-seq_len(start[[1L]])])))[1L] - 1L
  if (is.na(end)) end <- length(lines)
  obs(case, class, arm, "study_parameters.txt ITSD block", paste(lines[start[[1L]]:end], collapse = " | "))
}

qcFinalize <- function(binding, owner, config) {
  commitWorkflowOmicsStoreStep(owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"), config, "QC complete")
}

#' The integration input of the current owner (after log2 and cyclic loess, ITSD unticked): rows marked per view and
#' the count the model's rule drops (integrationInternalStandards, the summary's `internal_standards_dropped`).
integrationObs <- function(case, class, arm, binding, owner, config, ann_col) {
  capture(runWorkflowOmicsReferenceLog2(binding, owner, 0, config), case, class, arm, "log2")
  capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config), case, class, arm, "between_sample")
  src <- capture(withWorkflowStateScope(owner, function() omicsIntegrationInputSource(binding, owner)), case, class, arm, "omicsIntegrationInputSource")
  if (src$failed || is.null(src$value)) {
    obs(case, class, arm, "integration source", if (src$failed) paste("failed:", conditionMessage(src$value)) else "NULL")
    return(invisible())
  }
  for (view in names(src$value$features)) {
    f <- src$value$features[[view]]
    marked <- f$is_internal_standard %in% TRUE
    label <- if (!is.na(ann_col) && ann_col %in% names(f)) ifelse(is.na(f[[ann_col]]), "<NA annotation>", f[[ann_col]]) else f$feature_id
    obs(case, class, arm, paste(view, "rows marked is_internal_standard"), sum(marked))
    obs(case, class, arm, paste(view, "marked rows (annotation, else id)"), if (any(marked)) label[marked] else "none")
    if (!is.na(ann_col) && ann_col %in% names(f)) {
      obs(case, class, arm, paste(view, "missing-annotation row marked"), if (any(is.na(f[[ann_col]]))) f$is_internal_standard[is.na(f[[ann_col]])] else "no such row")
    }
    obs(case, class, arm, paste(view, "integrationInternalStandards() count (the model's drop)"), sum(integrationInternalStandards(f)))
    as_text <- f
    as_text$is_internal_standard <- as.character(as_text$is_internal_standard)
    obs(case, class, arm, paste(view, "integrationInternalStandards() count, column read back as text"),
      sum(integrationInternalStandards(as_text)))
  }
  obs(case, class, arm, "manifest standards rule", src$value$manifest$internal_standards %||% "absent")
}

# ---- D0: every default and every OD-1 option on the 20 names ---------------------------------------------------------
splitTop <- function(p) {
  chars <- strsplit(p, "")[[1L]]
  depth <- 0L; cur <- ""; out <- character(); esc <- FALSE; inclass <- FALSE
  for (ch in chars) {
    if (esc) { cur <- paste0(cur, ch); esc <- FALSE; next }
    if (ch == "\\") { cur <- paste0(cur, ch); esc <- TRUE; next }
    if (!inclass && ch == "[") inclass <- TRUE else if (inclass && ch == "]") inclass <- FALSE
    if (!inclass && ch == "(") depth <- depth + 1L
    if (!inclass && ch == ")") depth <- depth - 1L
    if (ch == "|" && depth == 0L && !inclass) { out <- c(out, cur); cur <- ""; next }
    cur <- paste0(cur, ch)
  }
  c(out, cur)
}

caseD0 <- function() {
  names20 <- NAMED$annotation
  lipid <- c("msdial", "progenesis", "xcms", "compound_discoverer", "lipidsearch", "custom")
  metab <- c("msdial", "progenesis", "xcms", "compound_discoverer", "custom")
  patterns <- list()
  for (f in lipid) patterns[[paste0("lipid_", f)]] <- getLipidomicsColumnDefaults(f)$is_pattern
  for (f in metab) patterns[[paste0("metab_", f)]] <- getMetabolomicsColumnDefaults(f)$is_pattern
  for (f in intersect(lipid, metab)) obs("D0", "-", "defaults", sprintf("lipid and metab '%s' strings identical", f),
    identical(patterns[[paste0("lipid_", f)]], patterns[[paste0("metab_", f)]]))
  obs("D0", "-", "defaults", "lipid msdial string at this commit", patterns$lipid_msdial)
  obs("D0", "-", "defaults", "metab msdial string at this commit", patterns$metab_msdial)
  obs("D0", "-", "defaults", "base MS-DIAL string equals the ticket's", identical(patterns$metab_msdial, OPTIONS[["base"]]))
  for (o in names(OPTIONS)) patterns[[paste0("option_", o)]] <- OPTIONS[[o]]
  matches <- list(); counts <- list()
  for (id in names(patterns)) {
    p <- patterns[[id]]
    key <- sub("^(lipid|metab|option)_", "", id)
    key <- if (identical(key, "msdial")) "base" else key
    alternatives <- if (is.na(p)) "whole" else c("whole", splitTop(p))
    for (alt in unique(alternatives)) {
      regex <- if (identical(alt, "whole")) p else alt
      hit <- tryCatch(stringr::str_detect(names20, regex), error = function(e) rep(NA, length(names20)))
      matches[[length(matches) + 1L]] <- data.frame(pattern_id = id, pattern = ifelse(is.na(p), "NA", p), alternative = alt,
        name = ifelse(is.na(names20), "<NA>", names20), role = NAMED$role, matched = ifelse(is.na(hit), "NA", as.character(hit)),
        stringsAsFactors = FALSE)
      observed <- sort(names20[hit %in% TRUE])
      predicted <- if (identical(alt, "whole")) PRED_SETS[[key]] else if (key == "base" && alt %in% names(PRED_BASE_ALT)) PRED_BASE_ALT[[alt]] else NULL
      counts[[length(counts) + 1L]] <- data.frame(pattern_id = id, alternative = alt, observed_n = length(observed),
        observed = paste(observed, collapse = " | "), na_on_missing = is.na(hit[is.na(names20)]),
        predicted_n = if (is.null(predicted)) NA_integer_ else if (is.character(predicted)) length(predicted) else as.integer(predicted),
        agrees = if (is.null(predicted)) NA else if (is.character(predicted)) identical(observed, predicted) else identical(length(observed), as.integer(predicted)),
        stringsAsFactors = FALSE)
    }
  }
  addTable("d0-matches.tsv", do.call(rbind, matches))
  addTable("d0-counts.tsv", do.call(rbind, counts))
}

# ---- D0b: what the import proposes per layout and omic -----------------------------------------------------------------
layoutTables <- function() {
  s <- data.frame(WT_1 = c(100, 110, 120), WT_2 = c(101, 111, 121), KO_1 = c(200, 210, 220), KO_2 = c(201, 211, 221))
  list(
    progenesis = cbind(data.frame(Compound = c("C1", "C2", "C3"), `Neutral mass (Da)` = c(100.1, 200.2, 300.3), `m/z` = c(101.1, 201.2, 301.3),
      Charge = 1, `Retention time (min)` = c(1.2, 2.3, 3.4), `Max Abundance` = c(1e4, 2e4, 3e4), check.names = FALSE), s),
    xcms = cbind(data.frame(featureid = c("FT0001", "FT0002", "FT0003"), mz = c(101.1, 201.2, 301.3), mzmin = c(101, 201, 301), mzmax = c(101.2, 201.3, 301.4),
      rt = c(72, 138, 204), rtmin = c(70, 136, 202), rtmax = c(74, 140, 206), npeaks = 4, check.names = FALSE), s),
    compound_discoverer = cbind(data.frame(`Compound Name` = c("Alanine", "Lactate", "Citrate"), `Molecular Formula` = c("C3H7NO2", "C3H6O3", "C6H8O7"),
      `Molecular Weight` = c(89.05, 90.03, 192.03), `RT [min]` = c(1.2, 2.3, 3.4), MZ = c(90.05, 91.04, 193.03), `Delta Mass [ppm]` = c(0.5, -0.3, 0.1),
      check.names = FALSE), s))
}
writeTsv <- function(data, path) { utils::write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE); path }

caseD0b <- function() {
  rows <- list()
  for (class in c("lipid", "metab")) {
    dir <- file.path(tempdir(), paste0("d0b-", class))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    files <- if (identical(class, "lipid")) c(
      lipidsearch = module_ci_lipid_write_table(module_ci_lipid_lipidsearch_assay(), file.path(dir, "ls.txt")),
      msdial = module_ci_lipid_write_table(module_ci_lipid_msdial_assay(), file.path(dir, "msdial_positive.tsv")),
      custom = module_ci_lipid_write_table(module_ci_lipid_custom_assay(), file.path(dir, "custom.tsv")),
      malformed = module_ci_lipid_write_malformed(file.path(dir, "malformed.tsv"))) else c(
      msdial = module_ci_metab_write_table(module_ci_metab_msdial_assay(), file.path(dir, "msdial_positive.tsv")),
      custom = module_ci_metab_write_table(module_ci_metab_custom_assay(), file.path(dir, "custom.tsv")),
      malformed = module_ci_lipid_write_malformed(file.path(dir, "malformed.tsv")))
    layouts <- layoutTables()
    names_on_disk <- c(progenesis = "progenesis_layout.tsv", xcms = "xcms_layout.tsv", compound_discoverer = "cd_layout.tsv")
    for (l in names(layouts)) files[[l]] <- writeTsv(layouts[[l]], file.path(dir, names_on_disk[[l]]))
    import <- omicGuiConfig(if (identical(class, "lipid")) "lipidomics" else "metabolomics")$import
    for (layout in names(files)) {
      run <- capture(previewOmicsImport(files[[layout]], import), "D0b", class, layout, "previewOmicsImport")
      v <- run$value
      rows[[length(rows) + 1L]] <- data.frame(class = class, layout = layout, file = basename(files[[layout]]),
        result = if (run$failed) paste("failed:", oneLine(conditionMessage(v))) else "ok",
        detected_format = if (run$failed) NA else as.character(v$detected_format %||% NA),
        confidence = if (run$failed) NA else as.character(round(v$format_confidence %||% NA, 3)),
        is_pattern = if (run$failed) NA else if (is.null(v$updates$is_pattern)) "NULL" else v$updates$is_pattern,
        id_selected = if (run$failed) NA else txt(v$updates$id$selected), annotation_selected = if (run$failed) NA else txt(v$updates$annotation$selected),
        stringsAsFactors = FALSE)
    }
  }
  addTable("d0b-previews.tsv", do.call(rbind, rows))
}

# ---- D1, D2: regex-mode ITSD through the declaration and the worker -------------------------------------------------
itsdCase <- function(case, class, object, name_col, root) {
  binding <- kit172Binding(class)
  prefix <- binding$prefix
  obs(case, class, "object", "regex / annotation column / id column", c(object@internal_standard_regex, object@annotation_id_column,
    methods::slot(object, binding$id_slot)))
  local({
    fixture <- g0Owner(file.path(root, tolower(case), class, "decl"), object, class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    arm <- "declaration (runSavedOperation)"
    qcFinalize(binding, owner, config)
    before <- g0State(owner, paste0(prefix, "_qc_complete"))
    run <- capture(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config), case, class, arm, "itsd")
    after <- g0State(owner, paste0(prefix, "_itsd_norm"))
    obs(case, class, arm, "step", if (run$failed) paste("failed:", conditionMessage(run$value)) else if (is.null(after)) "no state" else "completed")
    obs(case, class, arm, "console 'Identified ... via regex' lines", grep("Identified .* ITSD features via regex", run$messages, value = TRUE))
    obs(case, class, arm, "console 'Using ... ITSD features' lines", grep("Using [0-9NA]+ ITSD features", run$messages, value = TRUE))
    if (!is.null(after)) {
      namedPresence(case, class, arm, before, after, binding, name_col)
      itsdArgs(case, class, arm, after)
      studyParameters(case, class, arm, after, root)
    }
    stepObs(case, class, arm, g0Record(owner), paste0(prefix, "_itsd_norm"))
  })
  local({
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
    fixture <- g0Owner(file.path(root, tolower(case), class, "worker"), object, class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    arm <- "worker (omicsWorkerRun)"
    qcFinalize(binding, owner, config)
    before <- g0State(owner, paste0(prefix, "_qc_complete"))
    handoff <- capture(omicsWorkerRun(fixture, paste0(prefix, "_itsd_normalisation"), list(itsd_aggregation = "median", itsd_feature_ids = NULL)),
      case, class, arm, "omicsWorkerRun")
    h <- handoff$value
    if (handoff$failed) {
      obs(case, class, arm, "handoff", paste("omicsWorkerRun failed:", conditionMessage(h)))
      return(invisible())
    }
    obs(case, class, arm, "handoff status / error", c(h$status %||% "NULL", h$error %||% "none"))
    product <- tryCatch(rebuildDIACheckpointObject(owner$.context, h$result$checkpoint, TRUE, list()), error = function(e) e)
    if (inherits(product, "error")) {
      obs(case, class, arm, "product", paste("could not rebuild:", conditionMessage(product)))
    } else {
      namedPresence(case, class, arm, before, product, binding, name_col)
      itsdArgs(case, class, arm, product)
    }
  })
}

caseD1 <- function(root) itsdCase("D1", "metab", objectD1(), "Metabolite name", root)
caseD2 <- function(root) itsdCase("D2", "lipid", objectD2("Name", "Name"), "Name", root)

# ---- D2n: lipid with MS-DIAL's `Metabolite name` -------------------------------------------------------------------
caseD2n <- function(root) {
  class <- "lipid"
  binding <- kit172Binding(class)
  prefix <- binding$prefix
  # What the import proposes for a lipid table with MS-DIAL's header (Alignment ID, Metabolite name; no top rows).
  dir <- file.path(root, "d2n-preview")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  header <- data.frame(`Alignment ID` = 0:2, `Average Rt(min)` = c(1.2, 2.3, 3.4), `Average Mz` = c(760.6, 874.8, 703.6),
    `Metabolite name` = c("PC 34:1", "w/o MS2:TG 52:3", "SM 34:1"), `Adduct type` = "[M+H]+", Ontology = c("PC", "TG", "SM"),
    WT_1 = c(100, 110, 120), WT_2 = c(101, 111, 121), KO_1 = c(200, 210, 220), KO_2 = c(201, 211, 221), check.names = FALSE)
  file <- writeTsv(header, file.path(dir, "lipid_alignment_export.tsv"))
  pv <- capture(previewOmicsImport(file, omicGuiConfig("lipidomics")$import), "D2n", class, "preview", "previewOmicsImport")
  if (pv$failed) obs("D2n", class, "preview", "preview", paste("failed:", conditionMessage(pv$value))) else {
    obs("D2n", class, "preview", "detected_format / is_pattern", c(pv$value$detected_format, pv$value$updates$is_pattern %||% "NULL"))
    obs("D2n", class, "preview", "id selected / annotation selected", c(txt(pv$value$updates$id$selected), txt(pv$value$updates$annotation$selected)))
  }
  object <- objectD2("Metabolite name", NA_character_)
  obs("D2n", class, "object", "regex / annotation column", c(object@internal_standard_regex, object@annotation_id_column))
  local({
    fixture <- g0Owner(file.path(root, "d2n", "itsd"), object, class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    arm <- "declaration ITSD"
    qcFinalize(binding, owner, config)
    before <- g0State(owner, paste0(prefix, "_qc_complete"))
    run <- capture(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config), "D2n", class, arm, "itsd")
    after <- g0State(owner, paste0(prefix, "_itsd_norm"))
    obs("D2n", class, arm, "step", if (run$failed) paste("failed:", conditionMessage(run$value)) else if (is.null(after)) "no state" else "completed")
    if (!is.null(after)) {
      namedPresence("D2n", class, arm, before, after, binding, "Metabolite name")
      itsdArgs("D2n", class, arm, after)
    }
  })
  local({
    fixture <- g0Owner(file.path(root, "d2n", "integration"), object, class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    qcFinalize(binding, owner, config)
    integrationObs("D2n", class, "integration input, ITSD unticked", binding, owner, config, NA_character_)
  })
}

# ---- D3: a LipidSearch-built lipid object ----------------------------------------------------------------------------
caseD3 <- function(root) {
  class <- "lipid"
  binding <- kit172Binding(class)
  prefix <- binding$prefix
  SEED(3L)
  y <- drawLogs(42L)
  ids <- c(sprintf("PC %d:%d F%02d", 30 + (1:40 %% 10), 1:40 %% 4, 1:40), "PC(15:0_18:1(d7))", "IS_PC 15:0")
  tables <- lapply(y, function(m) {
    table <- data.frame(LipidName = ids, LipidClass = "PC", stringsAsFactors = FALSE)
    cbind(table, stats::setNames(as.data.frame(2^unname(m)), RUNS))
  })
  object <- lipidObject(tables, getLipidomicsColumnDefaults("lipidsearch")$is_pattern, "LipidName", "LipidClass")
  obs("D3", class, "object", "regex / annotation column", c(object@internal_standard_regex, object@annotation_id_column))
  local({
    fixture <- g0Owner(file.path(root, "d3", "itsd"), object, class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    arm <- "declaration ITSD"
    qcFinalize(binding, owner, config)
    before <- g0State(owner, paste0(prefix, "_qc_complete"))
    run <- capture(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config), "D3", class, arm, "itsd")
    after <- g0State(owner, paste0(prefix, "_itsd_norm"))
    obs("D3", class, arm, "step", if (run$failed) paste("failed:", conditionMessage(run$value)) else if (is.null(after)) "no state" else "completed")
    for (assay in c("A1", "A2")) {
      out_ids <- if (is.null(after)) character() else idsOf(after, binding, assay)
      obs("D3", class, arm, paste(assay, "rows in -> out"), sprintf("%d -> %d", length(idsOf(before, binding, assay)), length(out_ids)))
      obs("D3", class, arm, paste(assay, "'IS_PC 15:0' / 'PC(15:0_18:1(d7))' present after"), c("IS_PC 15:0", "PC(15:0_18:1(d7))") %in% out_ids)
    }
    itsdArgs("D3", class, arm, after)
  })
  local({
    fixture <- g0Owner(file.path(root, "d3", "integration"), object, class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    qcFinalize(binding, owner, config)
    integrationObs("D3", class, "integration input, ITSD unticked", binding, owner, config, "LipidClass")
  })
}

# ---- D4: D1's object with ITSD unticked, the integration input --------------------------------------------------------
caseD4 <- function(root) {
  class <- "metab"
  binding <- kit172Binding(class)
  local({
    fixture <- g0Owner(file.path(root, "d4"), objectD1(), class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    qcFinalize(binding, owner, config)
    integrationObs("D4", class, "integration input, ITSD unticked", binding, owner, config, "Metabolite name")
  })
}

# ---- D5: the QC tab's analysis with the box empty ----------------------------------------------------------------------
caseD5 <- function(root) {
  class <- "metab"
  binding <- kit172Binding(class)
  local({
    fixture <- g0Owner(file.path(root, "d5"), objectD1(), class)
    owner <- fixture$owner
    on.exit(closeQuietly(owner), add = TRUE)
    config <- fixture$workflow$config_list
    qcFinalize(binding, owner, config)
    run <- capture(analyzeWorkflowOmicsItsd(binding, owner, "", omicGuiConfig("metabolomics")$qc$itsd), "D5", class, "QC analysis", "analyzeWorkflowOmicsItsd")
    if (run$failed) {
      obs("D5", class, "QC analysis", "result", paste("failed:", conditionMessage(run$value)))
      return(invisible())
    }
    v <- run$value
    obs("D5", class, "QC analysis", "pattern used", v$pattern)
    obs("D5", class, "QC analysis", "total", v$total)
    for (a in names(v$metrics_by_assay)) {
      obs("D5", class, "QC analysis", paste(a, "standards reported"), nrow(v$metrics_by_assay[[a]]))
      obs("D5", class, "QC analysis", paste(a, "reported ids"), sort(as.character(v$metrics_by_assay[[a]]$is_id)))
    }
    obs("D5", class, "QC analysis", "logger 'Found IS in column' lines", grep("Found IS in column", run$logger, value = TRUE))
  })
}

# ---- D6: the factor and the deletion, recorded ------------------------------------------------------------------------
D6_ROUTES <- data.frame(method = c("cyclicloess", "quantile", "scale", "none", "cyclicloess"), ruv = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE)

#' Remove `share` of the cells of `rows`, lowest first within each row: total round(share * cells); each row gets
#' floor(total / rows), and the remainder goes one each to rows drawn at random (the seeded stream).
removeLowest <- function(m, rows, share) {
  total <- round(share * length(rows) * ncol(m))
  per <- rep(total %/% length(rows), length(rows))
  extra <- total - sum(per)
  if (extra > 0) { chosen <- sample.int(length(rows), extra); per[chosen] <- per[chosen] + 1L }
  for (i in seq_along(rows)) if (per[[i]] > 0) m[rows[[i]], order(m[rows[[i]], ])[seq_len(per[[i]])]] <- NA_real_
  m
}

caseD6 <- function(root) {
  class <- "metab"
  binding <- kit172Binding(class)
  prefix <- binding$prefix
  # Draws first (seed 260006): values, then the missing-value cells.
  SEED(6L)
  ann <- annotations(full = FALSE)
  y <- drawLogs(length(ann))
  ordinary_matched <- which(ann %in% NAMED$annotation[NAMED$role %in% c("wo_space", "capitalised")])
  intended <- which(ann %in% NAMED$annotation[NAMED$role == "intended"])
  group_b <- which(kit172Design()$group == "B")
  for (a in names(y)) y[[a]][ordinary_matched, group_b] <- y[[a]][ordinary_matched, group_b] + 1
  missing <- lapply(y, function(m) {
    other <- setdiff(seq_len(nrow(m)), ordinary_matched)
    removeLowest(removeLowest(m, ordinary_matched, 0.20), other, 0.05)
  })
  pattern <- getMetabolomicsColumnDefaults("msdial")$is_pattern
  obs("D6", class, "design", "rows per assay / ordinary matched rows / intended rows", c(length(ann), length(ordinary_matched), length(intended)))
  obs("D6", class, "design", "missing cells per assay (ordinary matched rows; other rows)",
    vapply(missing, function(m) sprintf("%d; %d", sum(is.na(m[ordinary_matched, ])), sum(is.na(m[-ordinary_matched, ]))), ""))
  versions <- list(complete = y, missing = missing)
  arms <- list()
  for (version in names(versions)) {
    full <- lapply(versions[[version]], assayTable, ann = ann, id_col = "Alignment ID", name_col = "Metabolite name")
    minus7 <- lapply(full, function(t) t[-ordinary_matched, , drop = FALSE])
    ids4 <- lapply(full, function(t) as.character(t[["Alignment ID"]][intended]))
    setups <- list(A = list(tables = full, ids = NULL), B = list(tables = minus7, ids = ids4), C = list(tables = full, ids = ids4))
    for (arm in names(setups)) {
      object <- metabObject(setups[[arm]]$tables, pattern)
      # The realised a_i: the median over observed raw values of the rows the arm uses as standards.
      a_i <- lapply(setups[[arm]]$tables, function(t) {
        used <- if (is.null(setups[[arm]]$ids)) which(stringr::str_detect(t[["Metabolite name"]], pattern) %in% TRUE) else
          which(t[["Metabolite name"]] %in% NAMED$annotation[NAMED$role == "intended"])
        apply(as.matrix(t[used, RUNS]), 2L, stats::median, na.rm = TRUE)
      })
      for (r in seq_len(nrow(D6_ROUTES))) {
        method <- D6_ROUTES$method[[r]]
        ruv <- D6_ROUTES$ruv[[r]]
        label <- sprintf("%s %s %s RUV %s", version, arm, method, if (ruv) "on" else "off")
        res <- local({
          fixture <- g0Owner(file.path(root, "d6", version, arm, paste0(method, if (ruv) "-ruv" else "")), object, class)
          owner <- fixture$owner
          on.exit(closeQuietly(owner), add = TRUE)
          config <- fixture$workflow$config_list
          qcFinalize(binding, owner, config)
          out <- list(status = "ok", k = list(), controls = list())
          step <- function(expr, name) {
            run <- capture(expr, "D6", class, label, name)
            if (run$failed && identical(out$status, "ok")) out$status <<- paste(name, "failed:", oneLine(conditionMessage(run$value)))
            run
          }
          itsd <- step(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", setups[[arm]]$ids, config), "itsd")
          after <- g0State(owner, paste0(prefix, "_itsd_norm"))
          out$itsd_counts <- if (is.null(after)) NULL else after@args$ITSDNormalization$itsd_counts_per_assay
          out$kept <- if (is.null(after)) NULL else vapply(c("A1", "A2"), function(a) length(idsOf(after, binding, a)), 1L)
          if (identical(out$status, "ok")) step(runWorkflowOmicsReferenceLog2(binding, owner, 0, config), "log2")
          if (identical(out$status, "ok")) step(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, method, config), "between_sample")
          if (ruv && identical(out$status, "ok")) {
            rv <- step(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config), "ruv")
            if (!rv$failed) {
              out$k <- rv$value$value$best_k
              out$controls <- rv$value$value$control_count
            }
          }
          if (identical(out$status, "ok")) {
            da <- step(runWorkflowOmicsReferenceDaAnalysis(binding, owner, getWorkflowStateReference(owner), DA_SETTINGS, kit172Contrasts()), "da")
            if (!da$failed) out$da <- as.data.frame(da$value$results[[binding$long_table]])
          }
          out
        })
        res$a_i <- a_i
        arms[[label]] <- res
        kOf <- function(a) { v <- res$k[[a]]; if (is.null(v) || !length(v)) NA else as.character(v[[1L]]) }
        cOf <- function(a) { v <- res$controls[[a]]; if (is.null(v) || !length(v)) NA else as.character(v[[1L]]) }
        addTable("d6-arms.tsv", data.frame(version = version, arm = arm, method = method, ruv = ruv, assay = c("A1", "A2"),
          status = res$status, itsd_count = vapply(c("A1", "A2"), function(a) as.character(res$itsd_counts[[a]] %||% NA), ""),
          rows_kept_after_itsd = vapply(c("A1", "A2"), function(a) as.character(res$kept[[a]] %||% NA), ""),
          ruv_k = if (ruv) vapply(c("A1", "A2"), kOf, "") else NA, ruv_controls = if (ruv) vapply(c("A1", "A2"), cOf, "") else NA,
          da_rows = vapply(c("A1", "A2"), function(a) if (is.null(res$da)) NA_integer_ else sum(res$da$assay == a), 1L),
          stringsAsFactors = FALSE))
      }
    }
  }
  # Pairs: A against B (the factor alone), B against C (the deletion alone); per route and assay.
  qcol <- function(d) intersect(c("fdr_qvalue", "q.mod", "qvalue"), names(d))[1L]
  idcol <- function(d) intersect(c("metabolite_id", "lipid_id", "feature"), names(d))[1L]
  rows <- list()
  for (version in names(versions)) for (pair in list(c("A", "B"), c("B", "C"))) for (r in seq_len(nrow(D6_ROUTES))) {
    method <- D6_ROUTES$method[[r]]
    ruv <- D6_ROUTES$ruv[[r]]
    lab <- function(arm) sprintf("%s %s %s RUV %s", version, arm, method, if (ruv) "on" else "off")
    first <- arms[[lab(pair[[1L]])]]
    second <- arms[[lab(pair[[2L]])]]
    for (assay in c("A1", "A2")) {
      row <- data.frame(version = version, pair = paste(pair, collapse = "-"), method = method, ruv = ruv, assay = assay,
        n_common = NA_integer_, max_abs_dlogfc = NA_real_, median_abs_dlogfc = NA_real_, constant = NA_real_, constant_max_dev = NA_real_,
        predicted_constant = NA_real_, up_first = NA_integer_, down_first = NA_integer_, up_second = NA_integer_, down_second = NA_integer_,
        ruv_k_first = NA_character_, ruv_controls_first = NA_character_, ruv_k_second = NA_character_, ruv_controls_second = NA_character_,
        note = "", stringsAsFactors = FALSE)
      if (ruv) {
        pickv <- function(x, a) { v <- x[[a]]; if (is.null(v) || !length(v)) NA_character_ else as.character(v[[1L]]) }
        row$ruv_k_first <- pickv(first$k, assay); row$ruv_controls_first <- pickv(first$controls, assay)
        row$ruv_k_second <- pickv(second$k, assay); row$ruv_controls_second <- pickv(second$controls, assay)
      }
      if (is.null(first$da) || is.null(second$da)) {
        row$note <- sprintf("no DA table (%s; %s)", first$status, second$status)
      } else {
        d1 <- first$da[first$da$assay == assay & first$da$comparison == "groupB-groupA", , drop = FALSE]
        d2 <- second$da[second$da$assay == assay & second$da$comparison == "groupB-groupA", , drop = FALSE]
        if (!nrow(d1)) d1 <- first$da[first$da$assay == assay, , drop = FALSE]
        if (!nrow(d2)) d2 <- second$da[second$da$assay == assay, , drop = FALSE]
        k1 <- as.character(d1[[idcol(d1)]]); k2 <- as.character(d2[[idcol(d2)]])
        common <- intersect(k1, k2)
        d <- d1$logFC[match(common, k1)] - d2$logFC[match(common, k2)]
        row$n_common <- length(common)
        row$max_abs_dlogfc <- max(abs(d), na.rm = TRUE)
        row$median_abs_dlogfc <- stats::median(abs(d), na.rm = TRUE)
        row$constant <- stats::median(d, na.rm = TRUE)
        row$constant_max_dev <- max(abs(d - row$constant), na.rm = TRUE)
        q1 <- d1[[qcol(d1)]]; q2 <- d2[[qcol(d2)]]
        row$up_first <- sum(q1 < 0.05 & d1$logFC > 0, na.rm = TRUE); row$down_first <- sum(q1 < 0.05 & d1$logFC < 0, na.rm = TRUE)
        row$up_second <- sum(q2 < 0.05 & d2$logFC > 0, na.rm = TRUE); row$down_second <- sum(q2 < 0.05 & d2$logFC < 0, na.rm = TRUE)
        g <- kit172Design()$group
        contrast <- function(v) mean(log2(v[g == "B"])) - mean(log2(v[g == "A"]))
        row$predicted_constant <- -(contrast(first$a_i[[assay]]) - contrast(second$a_i[[assay]]))
        row$note <- sprintf("DA columns: id %s, q %s", idcol(d1), qcol(d1))
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  addTable("d6-shift.tsv", do.call(rbind, rows))
}

# ---- D7: where the pattern lives, and what replay reads ---------------------------------------------------------------
caseD7 <- function(root, tree) {
  class <- "metab"
  binding <- kit172Binding(class)
  prefix <- binding$prefix
  default <- getMetabolomicsColumnDefaults("msdial")$is_pattern
  original_metab <- getMetabolomicsColumnDefaults
  original_lipid <- getLipidomicsColumnDefaults
  CHANGED <- "^KIT260_CHANGED_DEFAULT$"
  changedMetab <- function(format) { d <- original_metab(format); d$is_pattern <- CHANGED; d }
  changedLipid <- function(format) { d <- original_lipid(format); d$is_pattern <- CHANGED; d }
  # (a) the module-CI import with its MS-DIAL fixture, the box holding the proposed default.
  local({
    env <- environment()
    workflow <- new.env(parent = emptyenv())
    workflow$data_tbl <- NULL; workflow$data_format <- NULL; workflow$data_type <- NULL; workflow$column_mapping <- NULL
    workflow$processing_log <- list(setup_import = NULL)
    workflow$tab_status <- list(setup_import = "incomplete", design_matrix = "disabled")
    workflow <- omicsTestOwner(workflow = workflow, omic_type = "metabolomics", bind = FALSE, env = env)
    paths <- module_ci_metab_import_paths()
    assay <- module_ci_metab_msdial_assay(include_name = TRUE, duplicate_feature = FALSE, zero_missing = FALSE)
    file <- module_ci_metab_write_table(assay, file.path(paths$raw_dir, "msdial_positive.tsv"))
    preview <- previewOmicsImport(file, omicGuiConfig("metabolomics")$import)
    obs("D7", class, "import", "the box's proposed pattern (preview)", preview$updates$is_pattern %||% "NULL")
    run <- capture(module_ci_metab_import_run(assay = preview$import_result$data, paths = paths, workflow = workflow, assay_name = "LCMS_Pos",
      assay_file = file, vendor_format = "msdial", detected_format = preview$detected_format, metabolite_col = "Peak ID", annotation_col = "Name",
      sample_cols = preview$import_result$sample_columns, is_pattern = preview$updates$is_pattern), "D7", class, "import", "import")
    obs("D7", class, "import", "status", if (run$failed) paste("failed:", conditionMessage(run$value)) else run$value$result$status)
    if (run$failed) return(invisible())
    obs("D7", class, "import", "workflow$column_mapping$is_pattern", workflow$column_mapping$is_pattern %||% "absent")
    mapping_file <- module_ci_metab_import_artifacts(paths, "LCMS_Pos")[["column_mapping"]]
    mapping <- if (file.exists(mapping_file)) jsonlite::read_json(mapping_file, simplifyVector = TRUE) else NULL
    obs("D7", class, "import", "column_mapping.json is_pattern", mapping$is_pattern %||% "absent")
    files <- list.files(paths$base_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    holding <- files[vapply(files, function(f) tryCatch(any(grepl("w/o MS2:", readLines(f, warn = FALSE), fixed = TRUE)), error = function(e) FALSE), TRUE)]
    obs("D7", class, "import", "project files holding the pattern text", sub(paste0("^", paths$base_dir, "/?"), "", holding))
    # The design state's object from the stored mapping, with today's default replaced (instrument shown).
    testthat::local_mocked_bindings(getMetabolomicsColumnDefaults = changedMetab, getLipidomicsColumnDefaults = changedLipid, .package = "MultiScholaR")
    obs("D7", class, "design object", "instrument: the mocked default", getMetabolomicsColumnDefaults("msdial")$is_pattern)
    data <- getWorkflowOmicsBusTables(workflow, "data_tbl")
    samples <- preview$import_result$sample_columns
    design <- data.frame(Run = samples, group = sub("_[0-9]+$", "", samples), tech_rep_group = samples, stringsAsFactors = FALSE)
    built <- capture(metabDesignObject(list(data_cln = data, design_matrix = design, config_list = list()), mapping %||% workflow$column_mapping),
      "D7", class, "design object", "metabDesignObject")
    obs("D7", class, "design object", "internal_standard_regex of the design object",
      if (built$failed) paste("failed:", conditionMessage(built$value)) else built$value@internal_standard_regex)
  })
  # (b) a project with the stored base default: ITSD, log2, cyclic loess; its record; replayed with today's default changed.
  harness <- new.env(parent = globalenv())
  sys.source(file.path(tree, "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  SEED(7L)
  y <- drawLogs(60L)
  object <- metabObject(lapply(y, assayTable, ann = annotations(), id_col = "Alignment ID", name_col = "Metabolite name"), default)
  local({
    fixture <- g0Owner(file.path(root, "d7"), object, class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    arm <- "record and replay"
    qcFinalize(binding, owner, config)
    capture(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config), "D7", class, arm, "itsd (origin)")
    origin <- g0State(owner, paste0(prefix, "_itsd_norm"))
    obs("D7", class, arm, "origin itsd_counts_per_assay", origin@args$ITSDNormalization$itsd_counts_per_assay %||% "absent")
    capture(runWorkflowOmicsReferenceLog2(binding, owner, 0, config), "D7", class, arm, "log2 (origin)")
    capture(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config), "D7", class, arm, "between (origin)")
    json <- file.path(root, "d7", "record.json")
    capture(writeWorkflowAnalysisRecord(owner, json), "D7", class, arm, "write record")
    history <- unlist(owner$getHistory())
    closeQuietly(owner)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    obs("D7", class, arm, "record valid", validateWorkflowAnalysisRecord(record)$valid)
    obs("D7", class, arm, "record file contains the pattern text", any(grepl("w/o MS2:", readLines(json, warn = FALSE), fixed = TRUE)))
    for (s in record$steps) {
      hits <- Filter(function(p) grepl("pattern|regex", p$name, ignore.case = TRUE) ||
        grepl("w/o MS2", tryCatch(oneLine(workflowRecordValueText(p$value)), error = function(e) txt(p$value)), fixed = TRUE), s$parameters)
      obs("D7", class, arm, sprintf("record step %s (%s): parameters naming a pattern", s$state, s$operation),
        if (length(hits)) vapply(hits, function(p) paste0(p$name, "=", txt(p$value)), "") else "none")
    }
    source_dir <- fixture$source_dir
    copy <- file.path(root, "d7", "copy")
    dir.create(copy, recursive = TRUE)
    stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
    store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
    copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
    operations <- vapply(record$steps, function(s) s$operation, "")
    states <- vapply(record$steps, function(s) s$state, "")
    pre <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
    obs("D7", class, arm, "replay preconditions (taken from the stored project)", pre)
    files <- c(list.files(file.path(source_dir, "import_sources"), full.names = TRUE), file.path(source_dir, "design_matrix.tab"))
    sums <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "")
    inputs <- character()
    for (input in record$inputs) {
      hit <- files[sums == (input$sha256 %||% "")]
      if (length(hit)) inputs[[input$role]] <- hit[[1L]]
    }
    testthat::local_mocked_bindings(getMetabolomicsColumnDefaults = changedMetab, getLipidomicsColumnDefaults = changedLipid, .package = "MultiScholaR")
    obs("D7", class, arm, "instrument: the mocked default during replay", getMetabolomicsColumnDefaults("msdial")$is_pattern)
    prepared <- capture(harness$kitPrepareReplayProject(copy_root, pre, api), "D7", class, arm, "prepare replay")
    if (prepared$failed) {
      obs("D7", class, arm, "replay", paste("prepare failed:", conditionMessage(prepared$value)))
      return(invisible())
    }
    scratch <- file.path(root, "d7", "scratch")
    dir.create(scratch, recursive = TRUE)
    replay <- capture(harness$kitReplayJourney(record, prepared$value, store_root, inputs, fixture$paths, scratch, api), "D7", class, arm, "replay")
    if (replay$failed) {
      obs("D7", class, arm, "replay", paste("failed:", conditionMessage(replay$value)))
    } else {
      obs("D7", class, arm, "replay states equal history", identical(replay$value$states, history))
      obs("D7", class, arm, "replay mismatches", if (length(replay$value$mismatches)) as.character(replay$value$mismatches) else "none")
      obs("D7", class, arm, "replayed operations", vapply(replay$value$per_step, function(x) x$operation, ""))
      obs("D7", class, arm, "replay ITSD 'Identified/Using' console lines", grep("Using [0-9NA]+ ITSD features", replay$messages, value = TRUE))
    }
  })
}

# ---- the group runner -------------------------------------------------------------------------------------------------
runCase <- function(name, fn) {
  cat(sprintf("case %s ...\n", name))
  t0 <- proc.time()[[3L]]
  outcome <- tryCatch({ fn(); "ok" }, error = function(e) paste("case failed:", oneLine(conditionMessage(e))))
  if (!identical(outcome, "ok")) obs(name, "-", "-", "case", outcome)
  cat(sprintf("case %s: %s in %.1f s\n", name, outcome, proc.time()[[3L]] - t0))
}

runGroup <- function(group, tree) {
  installLogger()
  root <- file.path(tempdir(), paste0("g0-", group))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (identical(group, "d0")) {
    runCase("D0", caseD0)
    runCase("D0b", caseD0b)
  }
  if (identical(group, "core")) {
    runCase("D1", function() caseD1(root))
    runCase("D2", function() caseD2(root))
    runCase("D2n", function() caseD2n(root))
    runCase("D3", function() caseD3(root))
    runCase("D4", function() caseD4(root))
    runCase("D5", function() caseD5(root))
  }
  if (identical(group, "d6")) runCase("D6", function() caseD6(root))
  if (identical(group, "d7")) runCase("D7", function() caseD7(root, tree))
  bind <- function(x) if (length(x)) do.call(rbind, x)
  list(obs = bind(OUT$obs), cond = bind(OUT$cond), tables = OUT$tables)
}
})

# ---- child -----------------------------------------------------------------------------------------------------------
if (length(args) >= 1L && identical(args[[1L]], "--group")) {
  stopifnot(length(args) == 4L)
  group <- args[[2L]]
  tree <- normalizePath(args[[3L]])
  result <- tryCatch({
    suppressMessages(pkgload::load_all(tree, export_all = TRUE, helpers = TRUE, quiet = TRUE))
    env <- new.env(parent = asNamespace("MultiScholaR"))
    for (n in c("OPTIONS", "NAMED", "PRED_SETS", "PRED_BASE_ALT", "BASE_MSDIAL", "LIPIDSEARCH")) assign(n, get(n), envir = env)
    eval(caseDefinitions, env)
    env$runGroup(group, tree)
  }, error = function(e) list(error = gsub("[\t\r\n]+", " ", conditionMessage(e))))
  saveRDS(result, args[[4L]])
  quit(save = "no", status = 0L)
}

# ---- the job (parent) ------------------------------------------------------------------------------------------------
stopifnot(length(args) %in% 1:2)
out <- args[[1L]]
wanted <- if (length(args) == 2L) strsplit(args[[2L]], ",", fixed = TRUE)[[1L]] else c("d0", "core", "d6", "d7")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
writeTsv <- function(x, name) if (!is.null(x)) utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeTsv(PREDICTIONS, "predictions.tsv")
git <- function(...) suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
provenance <- c(sprintf("head\t%s", git("rev-parse", "HEAD")[[1L]]), sprintf("R\t%s", R.version.string),
  sprintf("stringr\t%s", as.character(utils::packageVersion("stringr"))), sprintf("stringi_icu\t%s", tryCatch(stringi::stri_info()$ICU.version, error = function(e) NA)),
  sprintf("limma\t%s", as.character(utils::packageVersion("limma"))), sprintf("started\t%s", format(Sys.time(), tz = "UTC")))
minutes <- c(d0 = 4, core = 12, d6 = 18, d7 = 8)
rscript <- file.path(R.home("bin"), "Rscript")
status <- list()
collected <- list()
for (g in intersect(names(minutes), wanted)) {
  rds <- file.path(out, paste0("child-", g, ".rds"))
  log <- file.path(out, paste0("child-", g, ".txt"))
  t0 <- proc.time()[[3L]]
  exit <- system2("timeout", c("--kill-after=30", paste0(minutes[[g]], "m"), shQuote(rscript), shQuote(script), "--group", g,
    shQuote(getwd()), shQuote(rds)), stdout = log, stderr = log)
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
writeTsv(bindPart("obs"), "observations.tsv")
writeTsv(bindPart("cond"), "conditions.tsv")
for (g in names(collected)) for (name in names(collected[[g]]$tables)) writeTsv(collected[[g]]$tables[[name]], name)
writeTsv(do.call(rbind, status), "groups.tsv")
writeLines(c(provenance, sprintf("finished\t%s", format(Sys.time(), tz = "UTC"))), file.path(out, "provenance.txt"))
cat(readLines(file.path(out, "groups.tsv")), sep = "\n")
