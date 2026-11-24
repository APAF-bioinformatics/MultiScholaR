# Fixtures of the replay tests (KIT-159 G5, 7.3; audit I07, re-check E5) and of the entry-script
# tests (re-check E2).

# A ticket on one omic route's generated corpus (the end-to-end tickets' form, one unit): the
# workflow arm alone, a data-set truth and a per-unit truth (E4).
replayTicket <- function(dir, route, pkg_root) {
    code <- file.path(dir, "ticket.R")
    writeLines(c(
        sprintf("ROUTE <- %s", deparse(route)),
        sprintf("PKG_ROOT <- %s", deparse(pkg_root)),
        "KEYS <- NULL",
        "prepareInputs <- function(data_dir, work_dir) {",
        "    corpus <- msrval::valWriteCorpus(msrval::valCorpusFor(ROUTE),",
        "                                     file.path(work_dir, 'corpus'), PKG_ROOT)",
        "    cell <- msrval::valCorpusCell(ROUTE, corpus)",
        "    feats <- do.call(rbind, lapply(seq_along(cell$assays), function(i) {",
        "        x <- utils::read.delim(cell$files[[i]], check.names = FALSE,",
        "                               stringsAsFactors = FALSE)",
        "        id <- intersect(c('Alignment ID', 'LipidName'), names(x))[[1L]]",
        "        data.frame(assay = cell$assays[[i]], feature_id = unique(as.character(x[[id]])),",
        "                   stringsAsFactors = FALSE)",
        "    }))",
        "    KEYS <<- unique(feats$feature_id)",
        "    one <- list(files = cell$files, design = cell$design,",
        "                contrasts = vapply(cell$contrasts, paste, '', collapse = '-'),",
        "                contrast_pairs = cell$contrasts, assays = cell$assays,",
        "                import_tool = 'msdial', forced_inputs = cell$forced_inputs)",
        "    list(units = list(u01 = one), features = feats)",
        "}",
        "truthOf <- function() data.frame(key = KEYS, class = 'all_features', expected = 0,",
        "                                 truth_status = 'scored')",
        "buildTruth <- function(data_dir) list(primary = truthOf())",
        "unitTruth <- function(unit, inputs) truthOf()",
        "featureKey <- function(feature_id, arm) {",
        "    data.frame(key = as.character(feature_id), status = 'ok', reason = NA_character_)",
        "}",
        "nullArm <- function(unit, inputs) stop('not run in this test')"
    ), code)
    code
}

replayDeclaration <- function(route, replay = TRUE) {
    list(
        ticket = "KIT-REPLAY", stage = "lifecycle",
        datasets = list(list(id = "d1", route = route, contrasts = list("unknown"), seed = 159600L,
                             lifecycle_units = list("u01"),
                             not_under_test = list(arms = list(), units = list()),
                             truth = list("primary", list(id = "unit", per_unit = "unitTruth")))),
        arms = list(list(id = "wf", kind = "workflow", builder = "valArmShipped", route = route),
                    list(id = "null", kind = "custom", `function` = "nullArm")),
        hooks = list(post_cell = NULL, summary = FALSE),
        metrics = list(list(id = "mae", `function` = "valMedianAbsError", args = list(n_min = 1L),
                            class = "accuracy", views = list("intersection"), truth = "primary",
                            call_rule = "own", higher_is_better = FALSE, contrast = "unknown",
                            margin = 0.5, conf = 0.95, interval = list(method = "paired_t"),
                            comparisons = list(list("wf", "null")))),
        replay = if (replay) list(list(arms = list("wf"), units = list("u01"))),
        refusal_tolerance = list(), blind = list(seed = 159601L),
        ceilings = list(cell_minutes = 30)
    )
}

# The installed msrval library (val-tests.R installs it into $RUNNER_TEMP/val-lib) offered to a
# child Rscript through a fresh RUNNER_TEMP, so each child has its own work folders.
childRunnerTemp <- function(envir = parent.frame()) {
    lib <- dirname(find.package("msrval"))
    tmp <- withr::local_tempdir(.local_envir = envir)
    if (!identical(basename(lib), "val-lib")) return(NA_character_)
    file.symlink(lib, file.path(tmp, "val-lib"))
    tmp
}

# Run an entry script from a folder; returns list(status, lines).
runEntryScript <- function(script, args, wd, runner_temp) {
    res <- processx::run(file.path(R.home("bin"), "Rscript"), c(script, args),
                         wd = wd, error_on_status = FALSE, timeout = 1800,
                         env = c("current", RUNNER_TEMP = runner_temp, VAL_TEMPLATE_FILE = ""))
    list(status = res$status, lines = strsplit(paste(res$stdout, res$stderr), "\n")[[1L]])
}

# The template lines an entry script printed.
statusLine <- function(lines) grep("^\\[msrval\\] ", trimws(lines), value = TRUE)
