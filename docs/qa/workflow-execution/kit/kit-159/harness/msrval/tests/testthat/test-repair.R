# The pre-release repairs (audit register AUDIT-KIT-159-impl-prerelease-20261124), on toy tickets
# shaped so that each audited fault would show: I01 to I09, I11, I14 to I18, I20. The end-to-end
# cell-to-verdict tests on generated inputs are in test-e2e.R; the route values in test-tab-repair.R.

prepRepair <- function(dir, ...) {
    env <- valTicketEnv(writeRepairTicket(dir, ...))
    list(env = env, prep = valPrepare(env, data_dir = dir, work_dir = file.path(dir, "work")))
}

test_that("R-I01: featureKey(feature_id, arm) is called per arm; unmapped is excluded for all arms and counted", {
    dir <- withr::local_tempdir()
    p <- prepRepair(dir, unmapped_in = "good")
    da <- list(good = p$env$goodArm("u01", NULL), null = p$env$nullArm("u01", NULL))
    rows <- valUnitMetrics(repairDeclaration(), "d1", "u01", da, p$prep$truth, p$prep$universe,
                           p$env)
    expect_true(all(is.finite(rows$value)))
    kc <- attr(rows, "key_counts")
    expect_identical(kc$n[kc$arm == "good" & kc$status == "unmapped"], 1L)
    expect_identical(kc$reason[kc$arm == "good" & kc$status == "unmapped"], "no_accession")
    tabs <- msrval:::valArmTables(da, p$env)
    s <- valScoreTable(tabs, p$prep$truth$primary, p$prep$universe, unit = "u01")
    v <- valViews(s, c("good", "null"), "intersection", "accuracy")
    expect_false("20" %in% v$key)
    expect_identical(attr(v, "excluded")$n, 1L)
    bare <- new.env()
    bare$featureKey <- function(feature_id, arm) feature_id
    expect_error(msrval:::valFeatureKeys(bare, c("a", "b"), "x"), class = "val_spec_invalid")
    harness <- Sys.getenv("VAL_HARNESS_DIR")
    skip_if(!nzchar(harness), "the harness folder is not set")
    tpl <- new.env()
    sys.source(file.path(harness, "templates", "prepare.R"), envir = tpl)
    fk <- msrval:::valFeatureKeys(tpl, c("P1", "P2"), "shipped")
    expect_identical(fk$key, c("P1", "P2"))
    expect_identical(fk$status, c("ok", "ok"))
})

test_that("R-I02: the universe is the input role's features with their assays; BH per assay; outside rows stop", {
    dir <- withr::local_tempdir()
    p <- prepRepair(dir, assays = c("LCMS_Pos", "GCMS"))
    expect_identical(nrow(p$prep$universe), 40L)
    expect_setequal(unique(p$prep$universe$assay), c("LCMS_Pos", "GCMS"))
    expect_match(p$prep$universe_sha256, "^[0-9a-f]{64}$")
    da <- list(good = p$env$goodArm("u01", NULL), null = p$env$nullArm("u01", NULL))
    rows <- valUnitMetrics(repairDeclaration(), "d1", "u01", da, p$prep$truth, p$prep$universe,
                           p$env)
    expect_true(all(is.finite(rows$value)))
    tabs <- msrval:::valArmTables(da, p$env)
    s <- valScoreTable(tabs, p$prep$truth$primary, p$prep$universe, unit = "u01")
    g <- s[s$arm == "good" & s$assay == "GCMS", ]
    expect_identical(nrow(g), 20L)
    expect_equal(g$s_bh, stats::p.adjust(g$p, "BH"))
    expect_identical(sum(is.finite(s$s_bh[s$arm == "good"])), 40L)
    da_all <- list(good = p$env$allAssayArm("u01", NULL), null = p$env$nullArm("u01", NULL))
    expect_error(valUnitMetrics(repairDeclaration(), "d1", "u01", da_all, p$prep$truth,
                                p$prep$universe, p$env), class = "val_spec_invalid")
    decl <- repairDeclaration()
    decl$datasets[[1L]]$expected_not_in_universe <- list(good = 20L)
    expect_silent(valUnitMetrics(decl, "d1", "u01", da_all, p$prep$truth, p$prep$universe,
                                 p$env))
    no_features <- new.env()
    for (f in c("buildTruth", "featureKey")) no_features[[f]] <- p$env[[f]]
    no_features$prepareInputs <- function(data_dir, work_dir) list(units = list(u01 = list()))
    expect_error(valPrepare(no_features, dir, file.path(dir, "w2")), class = "val_spec_invalid")
})

test_that("R-I03: ticket coverage regimes run, are pinned by hash and reach the classifier", {
    dir <- withr::local_tempdir()
    writeTicketRegimeGenerator(file.path(dir, "gen.R"))
    reg <- list(id = "pt3", method = "paired_t", class = "accuracy", g = 3L, levels = list(0.95),
                generator = "gen.R", simulations = 2000L, seed = 159201L)
    decl <- repairDeclaration()
    decl$coverage_regimes <- list(reg)
    tc <- valTicketCoverage(decl, dir)
    expect_identical(names(tc$coverage), msrval:::EVIDENCE_COLUMNS[["ticket-coverage.tsv"]])
    expect_setequal(tc$coverage$bound, c("two_sided", "lower", "upper"))
    expect_true(all(tc$coverage$pass))
    expect_equal(tc$coverage$coverage[tc$coverage$bound == "two_sided"], 0.95, tolerance = 0.02)
    expect_identical(names(tc$cells), msrval:::EVIDENCE_COLUMNS[["ticket-coverage-cells.tsv"]])
    strict <- decl
    strict$coverage_regimes[[1L]]$pass <- list(bound = "lower", max_miss = 0.001)
    expect_false(any(valTicketCoverage(strict, dir)$coverage$pass))
    few <- decl
    few$coverage_regimes[[1L]]$simulations <- 1999L
    expect_error(valTicketCoverage(few, dir), class = "val_spec_invalid")
    ev <- file.path(dir, "ev")
    dir.create(ev)
    path <- file.path(ev, "ticket-coverage.tsv")
    utils::write.table(tc$coverage, path, sep = "\t", quote = FALSE, row.names = FALSE)
    expect_error(valClassifierRegimes(decl, ev), class = "val_spec_invalid")
    decl$ticket_coverage <- list(path = "ticket-coverage.tsv", sha256 = sha256File(path))
    regimes <- valClassifierRegimes(decl, ev)
    expect_true(any(regimes$source == "ticket" & regimes$n == 3L & regimes$status == "pass"))
    metrics <- rbind(repairMetrics("mae", "good", c(0.10, 0.11, 0.09)),
                     repairMetrics("mae", "null", c(0.10, 0.10, 0.10)))
    cmp <- valRunLevel(decl, metrics, repairUnits(3L), regimes = regimes)
    expect_identical(cmp$regime[[1L]], "pt3")
    expect_false(is.na(cmp$outcome[[1L]]))
    failing <- tc$coverage
    failing$pass <- FALSE
    utils::write.table(failing, path, sep = "\t", quote = FALSE, row.names = FALSE)
    decl$ticket_coverage$sha256 <- sha256File(path)
    cmp2 <- valRunLevel(decl, metrics, repairUnits(3L), regimes = valClassifierRegimes(decl, ev))
    expect_true(is.na(cmp2$outcome[[1L]]))
    expect_match(cmp2$regime[[1L]], "^refused")
    decl$ticket_coverage$sha256 <- strrep("0", 64)
    expect_error(valClassifierRegimes(decl, ev), class = "val_declaration_changed")
})

test_that("R-I04: refusals are recorded per 4.7 and O7 (b); nothing aborts the scoring", {
    decl <- repairDeclaration()
    units21 <- repairUnits(21L)
    set.seed(1)
    v <- 0.1 + stats::rnorm(21, 0, 0.01)
    metrics <- rbind(repairMetrics("mae", "good", v), repairMetrics("mae", "null", v - 0.001))
    cells <- data.frame(unit = rep(units21$unit, 2), arm = rep(c("good", "null"), each = 21),
                        exit = "completed", stringsAsFactors = FALSE)
    cells$exit[cells$unit == "u21" & cells$arm == "good"] <- "refused"
    # (1) one refused unit of 21 with the harness's own table: no regime at 20, no abort
    cmp <- valRunLevel(decl, metrics, units21, cells = cells)
    expect_identical(cmp$n_refused[[1L]], 1L)
    expect_true(is.na(cmp$outcome[[1L]]))
    expect_match(cmp$flags[[1L]], "not_evaluable", fixed = TRUE)
    expect_match(cmp$regime[[1L]], "^not_evaluable:.*O7")
    # (2) two refused of fewer than 30: not_evaluable in flags
    cells2 <- cells
    cells2$exit[cells2$unit == "u20" & cells2$arm == "null"] <- "did_not_complete"
    cmp2 <- valRunLevel(decl, metrics, units21, cells = cells2)
    expect_match(cmp2$flags[[1L]], "not_evaluable", fixed = TRUE)
    expect_identical(cmp2$n_refused[[1L]], 2L)
    # with a registered regime at n - k the verdict stands, the substitute arm required
    cmp3 <- valRunLevel(decl, metrics, units21, cells = cells,
                        regimes = regimeRows("paired_t", 20L, class = "accuracy"))
    expect_false(is.na(cmp3$outcome[[1L]]))
    expect_match(cmp3$flags[[1L]], "refusal_substitute_required", fixed = TRUE)
    # (3) a one-arm comparison takes the rule
    one <- repairDeclaration()
    one$metrics[[1L]]$comparisons <- list(list("good"))
    one$metrics[[1L]]$fixed <- 0.1
    m1 <- repairMetrics("mae", "good", v, view = "intersection(good)")
    cmp4 <- valRunLevel(one, m1, units21, cells = cells[cells$arm == "good", ],
                        regimes = regimeRows("paired_t", 20L, class = "accuracy"))
    expect_identical(cmp4$n_refused[[1L]], 1L)
    expect_match(cmp4$flags[[1L]], "refusal_substitute_required", fixed = TRUE)
    # (4) draw units keep their Monte Carlo scope through the fill
    cmp5 <- valRunLevel(decl, metrics, repairUnits(21L, "draw"), cells = cells,
                        regimes = regimeRows("paired_t", 20L, class = "accuracy"))
    expect_identical(cmp5$scope[[1L]], "monte_carlo")
    # (5) a refused deletion run leaves the jackknife not evaluable, not undefined_deletion
    jk <- repairDeclaration()
    jk$metrics[[1L]]$interval <- list(method = "jackknife", form = "replicate")
    ju <- data.frame(unit = sprintf("u%02d", 1:5),
                     unit_role = c("full", rep("deletion", 4)), stringsAsFactors = FALSE)
    jm <- rbind(repairMetrics("mae", "good", c(0.1, 0.11, 0.12, 0.1, NA)),
                repairMetrics("mae", "null", c(0.1, 0.1, 0.1, 0.1, 0.1)))
    jc <- data.frame(unit = rep(ju$unit, 2), arm = rep(c("good", "null"), each = 5),
                     exit = "completed", stringsAsFactors = FALSE)
    jc$exit[jc$unit == "u05" & jc$arm == "good"] <- "refused"
    cmp6 <- valRunLevel(jk, jm, ju, cells = jc)
    expect_match(cmp6$flags[[1L]], "not_evaluable", fixed = TRUE)
    expect_false(grepl("undefined_deletion", cmp6$flags[[1L]], fixed = TRUE))
    # (6) a completed unit with a data NA is not a refusal
    md <- metrics
    md$value[md$unit == "u03" & md$arm == "good"] <- NA
    cmp7 <- valRunLevel(decl, md, units21, cells = data.frame(
        unit = rep(units21$unit, 2), arm = rep(c("good", "null"), each = 21), exit = "completed"),
        regimes = regimeRows("paired_t", 20L, class = "accuracy"))
    expect_identical(cmp7$n_refused[[1L]], 0L)
    expect_match(cmp7$flags[[1L]], "data_na=1", fixed = TRUE)
})

test_that("R-I05, I07, I08: a counted phase writes the run folder of 5.7, results.md with coverage, and passes both checks", {
    dir <- withr::local_tempdir()
    code <- writeRepairTicket(dir)
    decl <- repairDeclaration(conf = 0.975, family = TRUE)
    regimes <- rbind(regimeRows("paired_t", 5L, level = 0.975, class = "accuracy"))
    regimes$coverage <- 1
    out <- valRunPhase(phase = "counted", decl = decl, dataset = "d1", units = repairUnits(),
                       label = "toy", code_files = code, data_dir = dir,
                       work_dir = file.path(dir, "work"), out_dir = file.path(dir, "val-out"),
                       chunks = 2L, regimes = regimes, manifest = repairManifest())
    for (f in c("cells.tsv", "calls.tsv", "metrics.tsv", "comparisons.tsv", "coverage.tsv",
                "tables.sha256", "not-committed.sha256", "scored.tsv.gz", "results.md",
                "attribution.tsv", "allowed-paths.txt", "provenance.json", "outputs.sha256")) {
        expect_true(file.exists(file.path(out, f)), info = f)
    }
    calls <- utils::read.delim(file.path(out, "calls.tsv"))
    expect_identical(names(calls), msrval:::EVIDENCE_COLUMNS[["calls.tsv"]])
    expect_identical(nrow(calls), 10L)
    cov <- utils::read.delim(file.path(out, "coverage.tsv"))
    expect_identical(names(cov), msrval:::EVIDENCE_COLUMNS[["coverage.tsv"]])
    expect_setequal(unique(cov$class), c("changed", "null"))
    txt <- readLines(file.path(out, "results.md"))
    expect_true(any(grepl("measured coverage 1 ", txt, fixed = TRUE)))
    expect_true(any(grepl("fam: m = 2, nominal family level 0.95", txt, fixed = TRUE)))
    expect_true(any(grepl("entered at 3/simulations", txt, fixed = TRUE)))
    cmp <- utils::read.delim(file.path(out, "comparisons.tsv"), stringsAsFactors = FALSE)
    expect_true(all(is.na(cmp$outcome_family[!is.na(cmp$member)])))
    expect_true(all(grepl("family=fam", cmp$flags[!is.na(cmp$member)], fixed = TRUE)))
    prov <- jsonlite::read_json(file.path(out, "provenance.json"))
    expect_identical(length(prov$cells), 10L)
    expect_true(all(vapply(prov$cells, function(x) is.numeric(x$seed), logical(1))))
    expect_match(prov$inputs$universe$sha256, "^[0-9a-f]{64}$")
    expect_match(prov$scoring$metrics_sha256, "^[0-9a-f]{64}$")
    expect_true(valCheckProvenance(file.path(out, "provenance.json")))
    status <- checkEvidenceShell(out)
    if (!is.na(status)) expect_identical(status, 0L)
})

test_that("R-I08: the guard in production: a universe key as a label, a lifecycle call count, a stray file", {
    dir <- withr::local_tempdir()
    code <- writeRepairTicket(dir, extra = c(
        "LABEL_KEY <- TRUE",
        "callsMetric <- function(scored, ...) list(value = sum(scored$called_own %in% TRUE))"))
    body <- readLines(code)
    body[grepl("^    list\\(units = stats::setNames", body)] <- paste(
        "    list(result_tables = list('ticket/keys.tsv' = data.frame(label = '7', n = 1)),",
        "         units = stats::setNames(out, UNITS),")
    writeLines(body, code)
    decl <- repairDeclaration()
    decl$outputs <- list(list(path = "ticket/keys.tsv", key = "none", max_rows = 10L,
                              lineage = list("corpus"),
                              columns = list(list(name = "label", type = "label"),
                                             list(name = "n", type = "count")),
                              labels = list(label = list("7"))))
    expect_error(valRunPhase(phase = "counted", decl = decl, dataset = "d1",
                             units = repairUnits(), label = "k", code_files = code,
                             data_dir = dir, work_dir = file.path(dir, "w1"),
                             out_dir = file.path(dir, "o1"), chunks = 2L,
                             regimes = regimeRows("paired_t", 5L),
                             manifest = repairManifest("not_permitted")),
                 class = "val_publish_forbidden")
    life <- repairDeclaration(family = TRUE, conf = 0.975)
    life$metrics[[2L]]$id <- "calls"
    life$metrics[[2L]]$`function` <- "callsMetric"
    life$primary_family$fam$members$b$all_of <- list("calls:good-null")
    dir.create(file.path(dir, "t2"))
    code2 <- writeRepairTicket(file.path(dir, "t2"), extra = "callsMetric <- function(scored, ...) {
        list(value = sum(scored$called_own %in% TRUE)) }")
    expect_error(valRunPhase(phase = "lifecycle", decl = life, dataset = "d1",
                             units = repairUnits(), label = "l", code_files = code2,
                             data_dir = dir, work_dir = file.path(dir, "w2"),
                             out_dir = file.path(dir, "o2"), chunks = 2L,
                             manifest = repairManifest()),
                 class = "val_publish_forbidden")
    stray <- file.path(dir, "o3")
    dir.create(stray)
    writeLines("x", file.path(stray, "loose.txt"))
    expect_error(msrval:::valCheckOutsideRuns(stray), class = "val_publish_forbidden")
})

test_that("R-I06: phase score reads the requests' within-run rows; jackknife widens to the full unit's same view", {
    run <- valNewInterval(-0.01, 0.02, "jackknife", n_units = 4, conf = 0.975, form = "replicate")
    within <- data.frame(arm = "legacy", comparator = "baseline", metric = "mae",
                         scope = "within_run", view = c("intersection(legacy+baseline)", "union",
                                                        "intersection(legacy+baseline)"),
                         ci_low = c(-0.03, -0.5, -0.2), ci_high = c(0.01, 0.5, 0.2),
                         flags = c("unit=f1;n=10;status=ok", "unit=f1;n=10;status=ok",
                                   "unit=d1;n=10;status=ok"))
    w <- valWidenFromRequests(run, within, "legacy", "baseline", "mae",
                              view = "intersection(legacy+baseline)", unit = "f1")
    expect_identical(c(w$ci_low, w$ci_high), c(-0.03, 0.02))
    run8 <- valNewInterval(-0.01, 0.02, "jackknife", n_units = 8, conf = 0.975, form = "replicate")
    w8 <- valWidenFromRequests(run8, within, "legacy", "baseline", "mae", widen = FALSE,
                               view = "union")
    expect_true(any(grepl("^within_run_reported", w8$flags)))
    dir <- withr::local_tempdir()
    code <- writeRepairTicket(dir)
    decl <- repairDeclaration()
    decl$metrics[[1L]]$within_run <- list(B = 2000L, seed = 159103L)
    ev <- file.path(dir, "ev")
    for (lab in c("r1", "r2")) {
        ids <- if (lab == "r1") c("u01", "u02") else c("u03", "u04", "u05")
        o <- valRunPhase(phase = "counted", decl = decl, dataset = "d1", units = repairUnits(),
                         label = lab, code_files = code, data_dir = dir,
                         work_dir = file.path(dir, paste0("w", lab)),
                         out_dir = file.path(dir, paste0("o", lab)), chunks = 2L,
                         unit_ids = ids, manifest = repairManifest())
        dir.create(file.path(ev, "runs"), recursive = TRUE, showWarnings = FALSE)
        file.copy(o, file.path(ev, "runs"), recursive = TRUE)
    }
    s <- valScorePhase(decl, "d1", repairUnits(), ev, c("r1", "r2"), "s1",
                       out_dir = file.path(dir, "os"), manifest = repairManifest(),
                       regimes = regimeRows("paired_t", 5L, class = "accuracy"))
    cmp <- utils::read.delim(file.path(s, "comparisons.tsv"), stringsAsFactors = FALSE)
    expect_false(is.na(cmp$outcome[[1L]]))
    expect_identical(lengths(regmatches(cmp$flags[[1L]],
                                        gregexpr("within_run_reported=", cmp$flags[[1L]]))), 5L)
    expect_true(file.exists(file.path(s, "results.md")))
    prov <- jsonlite::read_json(file.path(s, "provenance.json"))
    expect_identical(length(prov$scoring$combined_runs), 2L)
})

test_that("R-I09, I18: declared states follow the arm's RUV run; a diagnostic step that changes state is not ok", {
    o <- valStepOutcome("diagnostic", before = "a", after = "b", history_changed = TRUE)
    expect_identical(o$outcome, "did_not_complete")
    expect_identical(valStepOutcome("diagnostic", "a", "a", FALSE)$outcome, "ok")
    norm <- function(route, ruv) {
        s <- valActionScript(route, ruv = ruv)$steps
        s[[which(vapply(s, `[[`, "", "id") == "norm")]]$declared$head
    }
    expect_identical(norm("prot_lfq", "skip"), "normalized")
    expect_identical(norm("prot_lfq", "applied"), "ruv_corrected")
    expect_setequal(norm("prot_lfq", "either"), c("normalized", "ruv_corrected"))
    expect_identical(norm("lipid_canonical", "skip"), "lipid_normalized")
    limpa <- valActionScript("prot_dia_limpa", ruv = "skip")$steps
    cov <- limpa[[which(vapply(limpa, `[[`, "", "id") == "ruv_covariates")]]
    expect_identical(cov$declared$head, "correlation_filtered")
    expect_identical(valArmRuvRun(valArmRuvOff()), "skip")
    expect_identical(valArmRuvRun(valArmRuvDesign("none")), "skip")
    expect_identical(valArmRuvRun(valArmRuvLegacy()), "applied")
    expect_identical(valArmRuvRun(valArmShipped()), "either")
})

test_that("R-I11: the declaration's statistical schema is refused before any cell", {
    ok <- repairDeclaration(family = TRUE, conf = 0.975)
    expect_true(valCheckDeclarationSchema(ok, repairUnits()))
    bad <- function(f) {
        d <- ok
        d <- f(d)
        expect_error(valCheckDeclarationSchema(d, repairUnits()), class = "val_spec_invalid")
    }
    bad(function(d) { d$metrics[[1L]]$margin <- NULL; d })
    bad(function(d) { d$metrics[[1L]]$contrast <- NULL; d })
    bad(function(d) { d$metrics[[1L]]$conf <- 0.95; d })
    bad(function(d) { d$primary_family$fam$m <- 5L; d })
    bad(function(d) { d$metrics[[1L]]$feature_interval_decides <- TRUE; d })
    bad(function(d) { d$metrics[[1L]]$interval$accept_flagged <- list(label = "x"); d })
    bad(function(d) { d$metrics[[1L]]$widen_to_within_run <- FALSE; d })
    bad(function(d) { d$metrics[[1L]]$interval$method <- "wilson"; d })
    bad(function(d) { d$metrics[[1L]]$interval$method <- "jackknife"; d })
    bad(function(d) { d$metrics[[1L]]$within_run <- list(B = 100L, seed = 1L); d })
    bad(function(d) { d$metrics[[1L]]$within_run <- list(B = 2000L); d })
    bad(function(d) { d$datasets[[1L]]$seed <- NULL; d })
    bad(function(d) { d$blind <- NULL; d })
    bad(function(d) { d$metrics[[1L]]$views <- list("union"); d })
    bad(function(d) {
        d$primary_family$fam$members$a$guards <- list(list(component = "mae:good-null"))
        d$primary_family$fam$members$a$worse_when <- list(list("mae:good-null"))
        d
    })
    bad(function(d) { d$metrics[[1L]]$`function` <- "valPartialAuc"; d })
})

test_that("R-I16: a member without a verdict is not_evaluable; worse_when never names a guard", {
    oc <- function(o, nw = FALSE) list(outcome = o, no_worse = nw)
    mem <- list(all_of = c("a", "b"))
    res <- valClassifyMember(mem, list(a = oc("match", TRUE), b = oc(NA_character_)))
    expect_identical(res$outcome, "not_evaluable")
    expect_identical(res$missing, "b")
    g <- list(all_of = c("a", "b"), guards = list(list(component = "b")),
              worse_when = list(list("b")))
    expect_error(valClassifyMember(g, list(a = oc("match", TRUE), b = oc("worse"))),
                 class = "val_spec_invalid")
    ns <- list(all_of = c("a", "c"), controls = list(c = "ctl"), worse_when = list(list("c")))
    res3 <- valClassifyMember(ns, list(a = oc("match", TRUE), c = oc("worse"),
                                       ctl = oc("match")))
    expect_identical(res3$outcome, "match")
})

test_that("R-I15: G1(a) has no limpa_da_last_digits class and no admission path (VAL-109 (5))", {
    expect_false("limpa_da_last_digits" %in% msrval:::G1A_CLASSES)
    ref <- data.frame(lane = "prot_dia_limpa", step = "da", class = "limpa_da_last_digits",
                      fields = "log2FC", keys = "1")
    live <- data.frame(class = "limpa_da_last_digits", digest = "8a4c2dd5ab")
    expect_identical(valReadChain("prot_dia_limpa", list(), ref, live_diffs = live)$status, "stop")
    # KIT-199's plan audit, ruling C1 (a): the DIA-with-limpa link is read as a declared change
    # (test-repair3.R, R-K199); the harness keeps no separate admission of the class
    expect_false("admit" %in% names(formals(valReadChain)))
})

test_that("R-I17: the chunk backstop ends a cell blocked outside R as did_not_complete (timeout)", {
    skip_if(!nzchar(Sys.which("sleep")), "no sleep command")
    dir <- withr::local_tempdir()
    runner <- file.path(dir, "runner.R")
    writeLines("runCell <- function(cell) { system2('sleep', '60'); list() }", runner)
    cells <- data.frame(cell = "u01:a", unit = "u01", run = "u01@a", arm = "a", seed = 1L)
    t0 <- Sys.time()
    res <- valRunCells(cells, runner, k = 1L, cell_minutes = 5, work_dir = dir,
                       chunk_minutes = 0.1)
    expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 45)
    expect_identical(res$exit, "did_not_complete")
    expect_identical(res$stop_reason, "timeout")
    expect_identical(msrval:::valChunkMinutes(n_runs = 9L, chunks = 4L, cell_minutes = 30), 100)
    # The fault (I17) is valRunPhase calling valRunCells without the backstop; R0 showed the
    # expectations above passing on the tests-only commit (valRunCells already had it), so the
    # call valRunPhase makes is captured.
    got <- NULL
    testthat::local_mocked_bindings(valRunCells = function(...) {
        got <<- list(...)
        stop("captured")
    }, .package = "msrval")
    d2 <- withr::local_tempdir()
    code <- writeRepairTicket(d2)
    try(valRunPhase(phase = "counted", decl = repairDeclaration(), dataset = "d1",
                    units = repairUnits(), label = "t", code_files = code, data_dir = d2,
                    work_dir = file.path(d2, "w"), out_dir = file.path(d2, "o"), chunks = 2L,
                    manifest = repairManifest()), silent = TRUE)
    expect_false(is.null(got[["chunk_minutes"]]))
    n_runs <- length(unique(got[[1L]]$run))
    expect_identical(got[["chunk_minutes"]], msrval:::valChunkMinutes(n_runs, 2L, 5))
})

test_that("R-I20: anchor refusal at step anchor; within-run seed required; confusion-table metrics", {
    arm <- list(id = "pub", kind = "published", fn = function(unit, inputs) data.frame(x = 1),
                anchor = function(da) FALSE)
    err <- tryCatch(valRunCustomArm(arm, "u01", list()), val_arm_refused = function(e) e)
    expect_s3_class(err, "val_arm_refused")
    expect_identical(err$step, "anchor")
    dir <- withr::local_tempdir()
    p <- prepRepair(dir)
    da <- list(good = p$env$goodArm("u01", NULL), null = p$env$nullArm("u01", NULL))
    noseed <- repairDeclaration()
    noseed$metrics[[1L]]$within_run <- list(B = 2000L)
    expect_error(valWithinRunRows(noseed, "u01", da, p$prep$truth, p$prep$universe, p$env),
                 class = "val_spec_invalid")
    conf <- repairDeclaration()
    conf$metrics[[1L]]$`function` <- "valFdp"
    conf$metrics[[1L]]$metric_input <- "confusion"
    conf$metrics[[1L]]$class <- "threshold_count"
    conf$metrics[[1L]]$args <- list()
    rows <- valUnitMetrics(conf, "d1", "u01", da, p$prep$truth, p$prep$universe, p$env)
    expect_true(all(is.finite(rows$value)))
})
