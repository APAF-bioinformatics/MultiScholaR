# The phases that run cells (KIT-159 section 7.8 steps 5 to 9): the cell table, the ticket's
# prepare.R, the cell runner by arm kind (workflow, kernel, baseline, published, custom), the
# kernel identity gate inside the runner (7.5, A14), per-unit metrics, the run-level intervals and
# classification, blind scoring in the lifecycle, and the run folder a phase writes.

toyTicket <- function(dir) {
    code <- file.path(dir, "ticket.R")
    writeLines(c(
        "prepareInputs <- function(data_dir, work_dir) {",
        "    units <- c('u01', 'u02', 'u03', 'u04', 'u05')",
        "    out <- lapply(seq_along(units), function(i) {",
        "        set.seed(100 + i)",
        "        m <- matrix(rnorm(80, 20, 0.2), 20, 4,",
        "                    dimnames = list(sprintf('k%02d', 1:20), c('A1', 'A2', 'B1', 'B2')))",
        "        m[1:5, 3:4] <- m[1:5, 3:4] + 1",
        "        f <- file.path(work_dir, paste0(units[[i]], '.tsv'))",
        "        utils::write.table(data.frame(id = rownames(m), m), f, sep = '\\t',",
        "                           quote = FALSE, row.names = FALSE)",
        "        list(files = list(quant = f),",
        "             design = data.frame(sample = colnames(m), group = c('A', 'A', 'B', 'B')),",
        "             contrasts = c('B-A' = 'B-A'))",
        "    })",
        "    list(units = stats::setNames(out, units),",
        "         features = data.frame(assay = 'all', feature_id = sprintf('k%02d', 1:20)))",
        "}",
        "buildTruth <- function(data_dir) {",
        "    list(primary = data.frame(key = sprintf('k%02d', 1:20),",
        "                              class = rep(c('changed', 'null'), c(5, 15)),",
        "                              expected = rep(c(1, 0), c(5, 15)), truth_status = 'scored'))",
        "}",
        "featureKey <- function(feature_id, arm) {",
        "    data.frame(key = feature_id, status = 'ok', reason = NA_character_)",
        "}",
        "customArm <- function(unit, inputs) {",
        "    q <- utils::read.delim(inputs$files$quant)",
        "    data.frame(feature_id = q$id, assay = 'all', contrast = 'B-A',",
        "               log2fc = rowMeans(q[, 4:5]) - rowMeans(q[, 2:3]), p = 0.5, s_own = 0.5,",
        "               call_statistic = 'none', called_own = FALSE, mean_expression = 20)",
        "}",
        "kernelSame <- function(unit, inputs) customArm(unit, inputs)",
        "kernelOther <- function(unit, inputs) {",
        "    d <- customArm(unit, inputs); d$log2fc <- d$log2fc + 1e-6; d",
        "}",
        "refusingArm <- function(unit, inputs) {",
        "    if (unit == 'u03') stop('too few replicates')",
        "    customArm(unit, inputs)",
        "}"
    ), code)
    code
}

toyDeclaration <- function(kernel = "kernelSame") {
    list(
        ticket = "KIT-TOY", stage = "counted",
        datasets = list(list(id = "d1", route = "none", contrasts = list("B-A"), seed = 159001L,
                             lifecycle_units = list("u01", "u02"),
                             not_under_test = list(arms = list(), units = list()),
                             truth = list("primary"))),
        arms = list(
            list(id = "baseline", kind = "baseline", normalise = "none", input_role = "quant",
                 column = "id"),
            list(id = "custom", kind = "custom", `function` = "customArm"),
            list(id = "kernel", kind = "kernel", `function` = kernel, identity_of = "custom",
                 identity_seed = 159102L)
        ),
        hooks = list(post_cell = NULL, summary = FALSE),
        metrics = list(list(
            id = "mae", `function` = "valMedianAbsError", args = list(n_min = 1L),
            class = "accuracy", views = list("intersection"), truth = "primary",
            call_rule = "own", higher_is_better = FALSE, contrast = "B-A", margin = 0.05,
            conf = 0.95, interval = list(method = "paired_t", form = NULL),
            comparisons = list(list("custom", "baseline")), primary_family = NULL)),
        refusal_tolerance = list(mae = 0L),
        blind = list(seed = 159101L),
        ceilings = list(cell_minutes = 5)
    )
}

toyUnits <- function() {
    data.frame(unit = sprintf("u%02d", 1:5), unit_role = "independent", form = "unit",
               stringsAsFactors = FALSE)
}

test_that("T-PH-1: the cell table crosses units and arms, one run per (run, arm), fixed seeds", {
    decl <- toyDeclaration()
    units <- toyUnits()
    cells <- valPhaseCells(decl, "d1", units, phase = "counted")
    expect_identical(nrow(cells), 15L)
    expect_true(all(c("cell", "unit", "run", "arm", "kind", "route", "seed") %in% names(cells)))
    expect_false(anyDuplicated(cells$run) > 0L)
    expect_identical(unique(cells$seed[cells$unit == "u02"]), 159002L)
    units$run <- c("u01", "u01", "u03", "u04", "u05")
    shared <- valPhaseCells(decl, "d1", units, phase = "counted", arms = "custom")
    expect_identical(shared$run[[1L]], shared$run[[2L]])
    life <- valPhaseCells(decl, "d1", toyUnits(), phase = "lifecycle")
    expect_setequal(unique(life$unit), c("u01", "u02"))
    expect_error(valPhaseCells(decl, "d1", toyUnits(), phase = "counted", arms = "nope"),
                 class = "val_spec_invalid")
})

test_that("T-PH-2: prepare.R runs in an isolated environment; universe and prepared hashes", {
    dir <- withr::local_tempdir()
    env <- valTicketEnv(toyTicket(dir))
    expect_false(exists("prepareInputs", envir = globalenv(), inherits = FALSE))
    work <- file.path(dir, "work")
    prep <- valPrepare(env, data_dir = dir, work_dir = work)
    expect_setequal(names(prep$units), sprintf("u%02d", 1:5))
    expect_true(file.exists(file.path(work, "universe.tsv")))
    expect_identical(nrow(prep$universe), 20L)
    expect_identical(nrow(prep$prepared), 5L)
    expect_true(all(grepl("^[0-9a-f]{64}$", prep$prepared$sha256)))
    bad <- new.env()
    bad$prepareInputs <- env$prepareInputs
    bad$featureKey <- env$featureKey
    bad$buildTruth <- function(data_dir) {
        list(primary = data.frame(key = c("k01", "k01"), class = "null", expected = 0,
                                  truth_status = "scored"))
    }
    expect_error(valPrepare(bad, data_dir = dir, work_dir = file.path(dir, "w2")),
                 class = "val_truth_invalid")
})

test_that("T-PH-3: the cell runner by kind writes da.tsv; a refusal is val_arm_refused", {
    dir <- withr::local_tempdir()
    env <- valTicketEnv(toyTicket(dir))
    prep <- valPrepare(env, data_dir = dir, work_dir = file.path(dir, "work"))
    decl <- toyDeclaration()
    ctx <- list(decl = decl, dataset = "d1", units = prep$units, env = env,
                work_dir = file.path(dir, "cells"), pkg_root = ".")
    cells <- valPhaseCells(decl, "d1", toyUnits(), phase = "counted")
    for (arm in c("baseline", "custom", "kernel")) {
        cell <- as.list(cells[cells$unit == "u01" & cells$arm == arm, ])
        res <- valRunCell(cell, ctx)
        da <- utils::read.delim(res$da)
        expect_true(all(c("feature_id", "assay", "contrast", "log2fc", "p", "s_own") %in%
                            names(da)))
        expect_identical(nrow(da), 20L)
    }
    decl$arms[[2L]][["function"]] <- "refusingArm"
    ctx$decl <- decl
    cell <- as.list(cells[cells$unit == "u03" & cells$arm == "custom", ])
    err <- tryCatch(valRunCell(cell, ctx), val_arm_refused = function(e) e)
    expect_s3_class(err, "val_arm_refused")
    expect_identical(err$refusal_class, "simpleError")
})

test_that("T-PH-4: kernel arms pass the identity gate before their cells, or are all refused", {
    dir <- withr::local_tempdir()
    env <- valTicketEnv(toyTicket(dir))
    prep <- valPrepare(env, data_dir = dir, work_dir = file.path(dir, "work"))
    units <- toyUnits()
    ok_ctx <- list(decl = toyDeclaration(), dataset = "d1", units = prep$units, env = env,
                   work_dir = file.path(dir, "ok"), pkg_root = ".")
    cells <- valPhaseCells(ok_ctx$decl, "d1", units, phase = "counted")
    gate <- valKernelGate(cells, ok_ctx)
    expect_true(gate$kernel$passed)
    expect_identical(length(gate$kernel$units), 3L)
    expect_identical(gate$refused, character())
    bad_ctx <- ok_ctx
    bad_ctx$decl <- toyDeclaration("kernelOther")
    bad_ctx$work_dir <- file.path(dir, "bad")
    gate <- valKernelGate(cells, bad_ctx)
    expect_false(gate$kernel$passed)
    expect_setequal(gate$refused, cells$run[cells$arm == "kernel"])
})

test_that("T-PH-5: per-unit metric rows carry the fixed columns; blinded in the lifecycle", {
    dir <- withr::local_tempdir()
    env <- valTicketEnv(toyTicket(dir))
    prep <- valPrepare(env, data_dir = dir, work_dir = file.path(dir, "work"))
    decl <- toyDeclaration()
    ctx <- list(decl = decl, dataset = "d1", units = prep$units, env = env,
                work_dir = file.path(dir, "cells"), pkg_root = ".")
    cells <- valPhaseCells(decl, "d1", toyUnits(), phase = "counted")
    u1 <- cells[cells$unit == "u01" & cells$arm %in% c("custom", "baseline"), ]
    da <- lapply(seq_len(nrow(u1)), function(i) valRunCell(as.list(u1[i, ]), ctx)$da)
    names(da) <- u1$arm
    rows <- valUnitMetrics(decl, "d1", "u01", da, prep$truth, prep$universe, env,
                           blind_seed = NULL)
    expect_identical(names(rows), msrval:::EVIDENCE_COLUMNS[["metrics.tsv"]])
    expect_setequal(rows$arm, c("custom", "baseline"))
    expect_true(all(is.finite(rows$value)))
    expect_false(any(rows$blinded))
    blind <- valUnitMetrics(decl, "d1", "u01", da, prep$truth, prep$universe, env,
                            blind_seed = 159101L)
    expect_true(all(blind$blinded))
})

test_that("T-PH-6: run-level intervals over units, refusals, and blinded outcomes", {
    decl <- toyDeclaration()
    m <- function(arm, values) {
        data.frame(unit = sprintf("u%02d", seq_along(values)), metric = "mae", arm = arm,
                   view = "intersection", contrast = "B-A", truth_id = "primary",
                   value = values, stringsAsFactors = FALSE)
    }
    metrics <- rbind(m("custom", c(0.10, 0.11, 0.09, 0.10, 0.12)),
                     m("baseline", c(0.10, 0.10, 0.10, 0.11, 0.11)))
    regimes <- rbind(regimeRows("paired_t", 5L), regimeRows("paired_t", 4L))
    cmp <- valRunLevel(decl, metrics, toyUnits(), regimes = regimes, phase = "counted")
    expect_identical(names(cmp), msrval:::EVIDENCE_COLUMNS[["comparisons.tsv"]])
    expect_identical(cmp$outcome[[1L]], "match")
    blind <- valRunLevel(decl, metrics, toyUnits(), regimes = regimes, phase = "lifecycle")
    expect_identical(blind$outcome[[1L]], "blinded")
    expect_true(blind$blinded[[1L]])
    metrics$value[metrics$arm == "custom" & metrics$unit == "u03"] <- NA
    refused <- valRunLevel(decl, metrics, toyUnits(), regimes = regimes, phase = "counted")
    expect_identical(refused$n_refused[[1L]], 1L)
})

test_that("T-PH-7: a phase writes its run folder (cells, metrics, comparisons, provenance)", {
    dir <- withr::local_tempdir()
    code <- toyTicket(dir)
    decl <- toyDeclaration()
    out <- valRunPhase(phase = "lifecycle", decl = decl, dataset = "d1", units = toyUnits(),
                       label = "toy", code_files = code, data_dir = dir,
                       work_dir = file.path(dir, "work"), out_dir = file.path(dir, "val-out"),
                       chunks = 2L, regimes = regimeRows("paired_t", 2L))
    run_dir <- file.path(dir, "val-out", "runs", "toy")
    expect_identical(out, run_dir)
    for (f in c("cells.tsv", "metrics.tsv", "comparisons.tsv", "provenance.json",
                "outputs.sha256")) {
        expect_true(file.exists(file.path(run_dir, f)), info = f)
    }
    cells <- utils::read.delim(file.path(run_dir, "cells.tsv"))
    expect_identical(names(cells), msrval:::EVIDENCE_COLUMNS[["cells.tsv"]])
    expect_setequal(unique(cells$unit), c("u01", "u02"))
    metrics <- utils::read.delim(file.path(run_dir, "metrics.tsv"))
    expect_true(all(metrics$blinded))
    expect_true(valCheckProvenance(file.path(run_dir, "provenance.json")))
    prov <- jsonlite::read_json(file.path(run_dir, "provenance.json"))
    expect_true(prov$blinded)
    expect_identical(prov$phase, "lifecycle")
})

test_that("T-PH-8: within-run bootstrap rows per unit, widening, clopper_pearson and a within-run verdict", {
    dir <- withr::local_tempdir()
    env <- valTicketEnv(toyTicket(dir))
    prep <- valPrepare(env, data_dir = dir, work_dir = file.path(dir, "work"))
    decl <- toyDeclaration()
    decl$metrics[[1L]]$within_run <- list(B = 2000L, seed = 159103L, strata = "class")
    ctx <- list(decl = decl, dataset = "d1", units = prep$units, env = env,
                work_dir = file.path(dir, "cells"), pkg_root = ".")
    cells <- valPhaseCells(decl, "d1", toyUnits(), phase = "counted")
    u1 <- cells[cells$unit == "u01" & cells$arm %in% c("custom", "baseline"), ]
    da <- lapply(seq_len(nrow(u1)), function(i) valRunCell(as.list(u1[i, ]), ctx)$da)
    names(da) <- u1$arm
    wr <- valWithinRunRows(decl, "u01", da, prep$truth, prep$universe, env, blind_seed = NULL)
    expect_identical(names(wr), msrval:::EVIDENCE_COLUMNS[["comparisons.tsv"]])
    expect_identical(wr$scope, "within_run")
    expect_identical(wr$method, "bootstrap_within_run")
    expect_true(is.finite(wr$ci_low) && is.finite(wr$ci_high))
    expect_match(wr$flags, "unit=u01", fixed = TRUE)
    m <- function(arm, values) {
        data.frame(unit = sprintf("u%02d", seq_along(values)), metric = "mae", arm = arm,
                   view = "intersection", contrast = "B-A", truth_id = "primary",
                   value = values, stringsAsFactors = FALSE)
    }
    metrics <- rbind(m("custom", c(0.10, 0.11, 0.09, 0.10, 0.12)),
                     m("baseline", c(0.10, 0.10, 0.10, 0.11, 0.11)))
    wide <- wr
    wide$ci_low <- -0.5
    wide$ci_high <- 0.5
    cmp <- valRunLevel(decl, metrics, toyUnits(), regimes = regimeRows("paired_t", 5L),
                       phase = "counted", within = wide)
    run <- cmp[cmp$scope == "run_level", ]
    # O8 (VAL-109): over independent units the unit's within-run interval is reported beside
    # the verdict, not used to widen it.
    expect_false(isTRUE(all.equal(run$ci_low[[1L]], -0.5)))
    expect_match(run$flags[[1L]], "within_run_reported=u01:[-0.5,0.5]", fixed = TRUE)
    cp <- toyDeclaration()
    cp$metrics[[1L]]$interval <- list(method = "clopper_pearson")
    cp$metrics[[1L]]$comparisons <- list(list("custom"))
    cp$metrics[[1L]]$fixed <- 0.5
    shares <- m("custom", c(1, 1, 0, 1, 1))
    res <- valRunLevel(cp, shares, toyUnits(), regimes = regimeRows("clopper_pearson", 5L,
                                                                   status = "conservative"),
                       phase = "counted")
    want <- valInterval(4L, method = "clopper_pearson", n = 5L, conf = 0.95)
    expect_equal(c(res$ci_low[[1L]], res$ci_high[[1L]]), c(want$ci_low, want$ci_high))
    wd <- toyDeclaration()
    wd$metrics[[1L]]$interval <- list(method = "bootstrap_within_run")
    wd$metrics[[1L]]$within_run <- list(B = 2000L, seed = 159103L)
    wd$metrics[[1L]]$feature_interval_decides <- TRUE
    verdict <- valRunLevel(wd, metrics, toyUnits(),
                           regimes = regimeRows("bootstrap_within_run", 20L), phase = "counted",
                           within = wr)
    expect_identical(verdict$scope[[1L]], "within_run")
    expect_false(is.na(verdict$outcome[[1L]]))
})

test_that("T-PH-9: test mode cannot be forced (r1.12, 7.3)", {
    dir <- withr::local_tempdir()
    env <- valTicketEnv(toyTicket(dir))
    prep <- valPrepare(env, data_dir = dir, work_dir = file.path(dir, "work"))
    decl <- toyDeclaration()
    decl$arms[[4L]] <- list(id = "shipped", kind = "workflow", builder = "valArmShipped",
                            route = "prot_lfq",
                            forced = list(list(setting = "multischolar.test_mode", value = FALSE,
                                               reason = "x")))
    ctx <- list(decl = decl, dataset = "d1", units = prep$units, env = env,
                work_dir = file.path(dir, "cells"), pkg_root = ".")
    cell <- as.list(valPhaseCells(decl, "d1", toyUnits(), "counted", arms = "shipped")[1L, ])
    expect_error(valRunCell(cell, ctx), class = "val_spec_invalid")
})
