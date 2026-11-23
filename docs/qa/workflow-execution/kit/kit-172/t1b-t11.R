# KIT-172 T1b and T11 (ticket section 6.4): measured, no pass rule. Run on the commit under test (plain log2).
#
# Usage (run-tests scripts lane, from the package root; `sva` among the request's extra packages):
#   Rscript docs/qa/workflow-execution/kit/kit-172/t1b-t11.R <out dir> [<projection limit, minutes>]
#
# T11: T7's concentration fixture and F172 (seed 172001) through the default route (ITSD on, cyclic loess, legacy RUV
# automatic, ~ 0 + group) at offset 1 and offset 0, lipid and metab: the largest change of a log2 value, of a logFC
# and of t, the calls at q < 0.05 that change, the chosen k and control share. -> t11.tsv
# T1b: F172 with one planted unwanted factor and planted effects, with and without missing cells, seeds 172101 to
# 172120 (lipid; seed 172101 also metab), rows `scale` (ITSD on, RUV skipped), legacy RUV `automatic` and `manual`
# (ITSD on, cyclic loess, the panel's defaults), and the new designs `groups_all_features` and `first_pass_nulls`
# (gate on, k rule `num_sv_be`, seed 1); factors 2^-10, 2^10 (with the intensity filter) and 0.001, 1000 (after it,
# on the features kept at s = 1), each against s = 1. -> t1b-seeds.tsv (one row per seed, missing variant, row and
# factor) and t1b-summary.tsv (changed seeds per row and the rule of section 6.4 applied). The first seed is timed and
# the whole projected before the rest run; a projection above the limit (default 69 minutes of wall time for the
# parallel part, the T1b row's ceiling of 4.6 core-hours on 4 cores) stops with status `stop`.
# Each T1b seed runs in its own Rscript process (`--seed <seed> <class> <rds>`), four at a time, so that no session
# owner (duckdb) is forked.
args <- commandArgs(trailingOnly = TRUE)
seed_mode <- identical(args[[1L]], "--seed")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
if (seed_mode) {
  out <- dirname(args[[4L]])
  status <- function(step, outcome) invisible(NULL)
} else {
  dir.create(args[[1L]], recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(args[[1L]])
  limit <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 69
  status_file <- file.path(out, "status.tsv")
  cat("step\toutcome\n", file = status_file)
  status <- function(step, outcome) cat(sprintf("%s\t%s\n", step, gsub("[\t\n]", " ", outcome)), file = status_file, append = TRUE)
}

suppressMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, attach_testthat = FALSE, quiet = TRUE))
env <- new.env(parent = globalenv())
suppressMessages(testthat::source_test_helpers("tests/testthat", env = env))
options(MultiScholaR.ruv_design_panel = TRUE)
if (!seed_mode) status("start", sprintf("commit %s; cores %s; sva %s", paste(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  collapse = ""), parallel::detectCores(), requireNamespace("sva", quietly = TRUE)))

local(envir = env, {
  factors <- c(`2^-10` = 2^-10, `2^10` = 2^10, `0.001` = 0.001, `1000` = 1000)

  calls <- function(da) ifelse(!is.na(da$fdr_qvalue) & da$fdr_qvalue < 0.05, sign(da$logFC), 0)
  ruvSummary <- function(journey) {
    value <- journey$ruv$value
    if (is.null(value)) return(list(k = NA_character_, share = NA_character_, controls = NA_character_))
    one <- function(x) if (is.null(x)) NA_character_ else paste(names(x), unlist(x), sep = "=", collapse = ";")
    list(k = one(value$best_k), share = one(value$percentage), controls = one(value$control_count))
  }
  # One journey's numbers: log2 and normalised matrices, the DA table, the RUV choice.
  run <- function(root, values, class, row, offset = 0, filter = TRUE) {
    spec <- switch(row,
      scale = list(method = "scale", ruv = NULL),
      legacy_automatic = list(method = "cyclicloess", ruv = kit172LegacyRuv("automatic")),
      legacy_manual = list(method = "cyclicloess", ruv = kit172LegacyRuv("manual")),
      groups_all_features = , first_pass_nulls = list(method = "cyclicloess", ruv = list(ruv_design = row,
        ruv_grouping_variable = "group", ruv_k = NULL, ruv_k_rule = "num_sv_be", ruv_k_seed = 1L, ruv_control_ids = NULL)))
    journey <- kit172Journey(root, kit172Object(values, class), class, itsd = TRUE, method = spec$method, filter = filter,
      offset = offset, ruv = spec$ruv)
    on.exit(closeDIAStateAdapter(journey$owner), add = TRUE)
    prefix <- journey$prefix
    list(log2 = kit172StateMatrix(journey, paste0(prefix, "_log2")),
      kept = rownames(kit172StateMatrix(journey, paste0(prefix, if (filter) "_intensity_filtered" else "_qc_complete"))),
      da = kit172DaTable(journey), ruv = ruvSummary(journey))
  }
  compare <- function(a, b) {
    key <- function(d) paste(d$assay, d$comparison, d$feature)
    at <- match(key(b), key(a))
    list(max_dlogfc = max(abs(a$logFC[at] - b$logFC), na.rm = TRUE), max_dt = max(abs(a$t[at] - b$t), na.rm = TRUE),
      call_changes = sum(calls(a)[at] != calls(b), na.rm = TRUE))
  }

  rows <- c("scale", "legacy_automatic", "legacy_manual", "groups_all_features", "first_pass_nulls")
  merge_rows <- function(frames) {
    columns <- unique(unlist(lapply(frames, names)))
    do.call(rbind, lapply(frames, function(f) { for (c in setdiff(columns, names(f))) f[[c]] <- NA; f[columns] }))
  }
  oneSeed <- function(seed, class) {
    set.seed(seed)
    out_rows <- list()
    for (missing in c(TRUE, FALSE)) {
      values <- kit172Values(seed, planted = TRUE, missing = missing)$values
      for (row in rows) {
        root <- tempfile(sprintf("t1b-%d-", seed))
        reference <- tryCatch(kit172Quiet(run(file.path(root, "s1"), values, class, row)), error = function(e) e)
        if (inherits(reference, "error")) {
          out_rows[[length(out_rows) + 1L]] <- data.frame(seed = seed, class = class, missing = missing, row = row, factor = "1",
            error = conditionMessage(reference), stringsAsFactors = FALSE)
          next
        }
        for (name in names(factors)) {
          s <- factors[[name]]
          power <- name %in% c("2^-10", "2^10")
          input <- if (power) values * s else values[reference$kept, , drop = FALSE] * s
          got <- tryCatch(kit172Quiet(run(file.path(root, name), input, class, row, filter = power)), error = function(e) e)
          if (inherits(got, "error")) {
            out_rows[[length(out_rows) + 1L]] <- data.frame(seed = seed, class = class, missing = missing, row = row, factor = name,
              error = conditionMessage(got), stringsAsFactors = FALSE)
            next
          }
          d <- compare(got$da, reference$da)
          out_rows[[length(out_rows) + 1L]] <- data.frame(seed = seed, class = class, missing = missing, row = row, factor = name,
            error = "", max_dlogfc = d$max_dlogfc, max_dt = d$max_dt, call_changes = d$call_changes,
            k_1 = reference$ruv$k, k_s = got$ruv$k, share_1 = reference$ruv$share, share_s = got$ruv$share,
            controls_1 = reference$ruv$controls, controls_s = got$ruv$controls,
            changed = d$call_changes > 0 || !identical(reference$ruv$k, got$ruv$k) || !identical(reference$ruv$share, got$ruv$share),
            stringsAsFactors = FALSE)
        }
        unlink(root, recursive = TRUE)
      }
    }
    merge_rows(out_rows)
  }
  if (seed_mode) {
    saveRDS(oneSeed(as.integer(args[[2L]]), args[[3L]]), args[[4L]])
    return(invisible(NULL))
  }
  seedProcess <- function(seed, class) {
    rds <- file.path(out, sprintf("seed-%d-%s.rds", seed, class))
    code <- system2("Rscript", c(shQuote(script), "--seed", seed, class, shQuote(rds)), stdout = FALSE, stderr = FALSE)
    if (code == 0L && file.exists(rds)) {
      x <- readRDS(rds)
      unlink(rds)
      x
    } else data.frame(seed = seed, class = class, error = sprintf("seed process exit %d", code), stringsAsFactors = FALSE)
  }

  # --- T11 -------------------------------------------------------------------------------------------------------
  t11 <- list()
  for (fixture in c("t7_concentration", "f172")) for (class in c("lipid", "metab")) {
    values <- if (identical(fixture, "f172")) kit172Values(172001L)$values else kit172ConcentrationValues()
    root <- tempfile("t11-")
    r1 <- tryCatch(kit172Quiet(run(file.path(root, "o1"), values, class, "legacy_automatic", offset = 1)), error = function(e) e)
    r0 <- tryCatch(kit172Quiet(run(file.path(root, "o0"), values, class, "legacy_automatic", offset = 0)), error = function(e) e)
    if (inherits(r1, "error") || inherits(r0, "error")) {
      status(paste("t11", fixture, class), paste("mechanism:", conditionMessage(if (inherits(r1, "error")) r1 else r0)))
      next
    }
    common <- intersect(rownames(r0$log2), rownames(r1$log2))
    d <- compare(r0$da, r1$da)
    t11[[length(t11) + 1L]] <- data.frame(fixture = fixture, class = class,
      max_dlog2 = max(abs(r0$log2[common, ] - r1$log2[common, ]), na.rm = TRUE), max_dlogfc = d$max_dlogfc, max_dt = d$max_dt,
      call_changes = d$call_changes, calls_offset1 = sum(calls(r1$da) != 0), calls_offset0 = sum(calls(r0$da) != 0),
      k_offset1 = r1$ruv$k, k_offset0 = r0$ruv$k, share_offset1 = r1$ruv$share, share_offset0 = r0$ruv$share,
      controls_offset1 = r1$ruv$controls, controls_offset0 = r0$ruv$controls, stringsAsFactors = FALSE)
  }
  utils::write.table(do.call(rbind, t11), file.path(out, "t11.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  status("t11", sprintf("pass: %d rows", length(t11)))

  # --- T1b -------------------------------------------------------------------------------------------------------
  started <- Sys.time()
  first <- seedProcess(172101L, "lipid")
  first_minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  projection <- first_minutes * 20 / 4
  status("t1b_first_seed", sprintf("%.1f minutes; projected parallel part %.1f minutes (limit %.0f)", first_minutes, projection, limit))
  if (projection > limit) {
    utils::write.table(first, file.path(out, "t1b-seeds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    status("t1b", sprintf("stop: projection %.1f minutes above the limit %.0f", projection, limit))
    return(invisible(NULL))
  }
  tasks <- c(lapply(172102:172120, function(seed) list(seed = seed, class = "lipid")), list(list(seed = 172101L, class = "metab")))
  rest <- parallel::mclapply(tasks, function(task) seedProcess(task$seed, task$class), mc.cores = 4L, mc.preschedule = FALSE)
  seeds <- merge_rows(c(list(first), rest))
  utils::write.table(seeds, file.path(out, "t1b-seeds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  status("t1b_seeds", sprintf("%d rows, %d with an error, %.1f minutes in all", nrow(seeds), sum(nzchar(seeds$error %||% "")),
    as.numeric(difftime(Sys.time(), started, units = "mins"))))

  # The rule of section 6.4 (D6, D8), lipid seeds 172101 to 172120.
  lipid <- seeds[seeds$class == "lipid" & !nzchar(seeds$error), , drop = FALSE]
  summary <- do.call(rbind, lapply(rows, function(row) {
    at <- lipid[lipid$row == row, , drop = FALSE]
    call_seeds <- unique(at$seed[at$call_changes > 0])
    changed_seeds <- unique(at$seed[at$changed])
    rule <- if (row == "scale") length(call_seeds) > 1L else if (startsWith(row, "legacy_")) length(changed_seeds) > 1L else NA
    data.frame(row = row, seeds_run = length(unique(at$seed)), seeds_with_call_change = length(call_seeds),
      seeds_changed = length(changed_seeds), max_dt = if (nrow(at)) max(at$max_dt, na.rm = TRUE) else NA_real_,
      max_dlogfc = if (nrow(at)) max(at$max_dlogfc, na.rm = TRUE) else NA_real_,
      rule = if (is.na(rule)) "measured (no rule: new design)" else if (rule) "own ticket (more than 1 seed of 20)" else
        "stated limitation (at most 1 seed of 20)", stringsAsFactors = FALSE)
  }))
  utils::write.table(summary, file.path(out, "t1b-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  # Seed 172101 on metab must give the lipid tables.
  pick <- function(class) {
    x <- seeds[seeds$seed == 172101L & seeds$class == class, setdiff(names(seeds), "class"), drop = FALSE]
    x[order(x$missing, x$row, x$factor), , drop = FALSE]
  }
  same <- isTRUE(all.equal(pick("lipid"), pick("metab"), check.attributes = FALSE))
  status("t1b_metab_seed", if (same) "pass: identical to lipid" else "differs from lipid (see t1b-seeds.tsv)")
  status("t1b", "pass")
})
