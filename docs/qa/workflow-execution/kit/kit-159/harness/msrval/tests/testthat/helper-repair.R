# Fixtures of the pre-release repair tests (audit register AUDIT-KIT-159-impl-prerelease-20261124,
# I01 to I20): a toy ticket under the template's own contracts (featureKey(feature_id, arm)
# returning key, status and reason; prepareInputs() returning the input role's features), and
# the end-to-end tickets on the generated corpora.

# A toy ticket: UNITS units, features "1".."20" in each of ASSAYS (one compound in every assay),
# truth 5 changed and 15 null; feature "20" unmapped in arm UNMAPPED_IN; arms of kind custom.
writeRepairTicket <- function(dir, n_units = 5L, assays = "all", unmapped_in = "none",
                              extra = character()) {
    code <- file.path(dir, "ticket.R")
    writeLines(c(
        sprintf("UNITS <- sprintf('u%%02d', seq_len(%d))", as.integer(n_units)),
        sprintf("ASSAYS <- c(%s)", paste(sprintf("'%s'", assays), collapse = ", ")),
        sprintf("UNMAPPED_IN <- '%s'", unmapped_in),
        "KEYS <- as.character(1:20)",
        "prepareInputs <- function(data_dir, work_dir) {",
        "    out <- lapply(UNITS, function(u) list(files = list(), contrasts = c('B-A' = 'B-A')))",
        "    list(units = stats::setNames(out, UNITS),",
        "         features = expand.grid(feature_id = KEYS, assay = ASSAYS,",
        "                                stringsAsFactors = FALSE))",
        "}",
        "buildTruth <- function(data_dir) {",
        "    list(primary = data.frame(key = KEYS, class = rep(c('changed', 'null'), c(5, 15)),",
        "                              expected = rep(c(1, 0), c(5, 15)), truth_status = 'scored'))",
        "}",
        "featureKey <- function(feature_id, arm) {",
        "    bad <- !is.na(arm) & arm == UNMAPPED_IN & feature_id == '20'",
        "    data.frame(key = feature_id, status = ifelse(bad, 'unmapped', 'ok'),",
        "               reason = ifelse(bad, 'no_accession', NA_character_))",
        "}",
        "armRows <- function(unit, shift, assays = ASSAYS) {",
        "    i <- match(unit, UNITS)",
        "    d <- expand.grid(feature_id = KEYS, assay = assays, stringsAsFactors = FALSE)",
        "    truth <- rep(c(1, 0), c(5, 15))[match(d$feature_id, KEYS)]",
        "    j <- seq_len(nrow(d))",
        "    data.frame(feature_id = d$feature_id, assay = d$assay, contrast = 'B-A',",
        "               log2fc = truth + shift + 0.01 * i * (j %% 3),",
        "               p = ifelse(truth > 0, 0.001, 0.2 + 0.03 * (j %% 20)),",
        "               s_own = ifelse(truth > 0, 0.01, 0.6), call_statistic = 'none',",
        "               called_own = truth > 0, mean_expression = 20)",
        "}",
        "goodArm <- function(unit, inputs) armRows(unit, 0)",
        "nullArm <- function(unit, inputs) armRows(unit, 0.05)",
        "allAssayArm <- function(unit, inputs) armRows(unit, 0, assays = 'all')",
        "refusingArm <- function(unit, inputs) {",
        "    if (unit == 'u04') stop('too few replicates')",
        "    armRows(unit, 0)",
        "}",
        extra), code)
    code
}

repairDeclaration <- function(arms = c("good", "null"), conf = 0.95, family = FALSE) {
    fn <- c(good = "goodArm", null = "nullArm", refusing = "refusingArm", all = "allAssayArm")
    metric <- function(id, cmp) {
        list(id = id, `function` = "valMedianAbsError", args = list(n_min = 1L),
             class = "accuracy", views = list("intersection"), truth = "primary",
             call_rule = "own", higher_is_better = FALSE, contrast = "B-A", margin = 0.5,
             conf = conf, interval = list(method = "paired_t"), comparisons = list(as.list(cmp)),
             primary_family = if (family) "fam")
    }
    decl <- list(
        ticket = "KIT-TOY", stage = "counted",
        datasets = list(list(id = "d1", route = "none", contrasts = list("B-A"), seed = 159001L,
                             lifecycle_units = list("u01", "u02"),
                             not_under_test = list(arms = list(), units = list()),
                             truth = list("primary"))),
        arms = lapply(arms, function(a) list(id = a, kind = "custom", `function` = fn[[a]])),
        hooks = list(post_cell = NULL, summary = FALSE),
        metrics = list(metric("mae", arms[1:2])),
        refusal_tolerance = list(), blind = list(seed = 159101L),
        ceilings = list(cell_minutes = 5)
    )
    if (family) {
        decl$metrics[[2L]] <- metric("mae2", arms[1:2])
        decl$metrics[[2L]]$higher_is_better <- TRUE
        decl$primary_family <- list(fam = list(m = 2L, members = list(
            a = list(all_of = list(paste0("mae:", arms[1], "-", arms[2]))),
            b = list(all_of = list(paste0("mae2:", arms[1], "-", arms[2]))))))
    }
    decl
}

repairUnits <- function(n = 5L, role = "independent") {
    data.frame(unit = sprintf("u%02d", seq_len(n)), unit_role = role, stringsAsFactors = FALSE)
}

repairManifest <- function(redistribution = "permitted") {
    data.frame(role = "corpus", kind = "data", licence = "generated",
               redistribution = redistribution, citation = "synthetic", sha256 = strrep("0", 64),
               commit = FALSE, file = "none", stringsAsFactors = FALSE)
}

# Per-unit metric rows for valRunLevel.
repairMetrics <- function(metric, arm, values, view = "intersection(good+null)") {
    data.frame(unit = sprintf("u%02d", seq_along(values)), metric = metric, arm = arm,
               view = view, contrast = "B-A", truth_id = "primary", value = values,
               stringsAsFactors = FALSE)
}

# A coverage-regime generator: the paired t interval over 3 normal units (2,000 simulations).
writeTicketRegimeGenerator <- function(path) {
    writeLines(c(
        "simulate <- function(nsim, seed, levels) {",
        "    set.seed(seed)",
        "    x <- matrix(stats::rnorm(3 * nsim), 3)",
        "    m <- colMeans(x)",
        "    se <- apply(x, 2, stats::sd) / sqrt(3)",
        "    do.call(rbind, lapply(levels, function(lv) {",
        "        q <- stats::qt(1 - (1 - lv) / 2, 2)",
        "        data.frame(level = lv, covered = m - q * se <= 0 & 0 <= m + q * se,",
        "                   miss_low = m - q * se > 0, miss_high = m + q * se < 0)",
        "    }))",
        "}"), path)
    path
}

# The shell pre-push check on a folder: its exit status (NA where bash or the harness is absent).
checkEvidenceShell <- function(dir) {
    harness <- Sys.getenv("VAL_HARNESS_DIR")
    if (!nzchar(harness) || !nzchar(Sys.which("bash"))) return(NA_integer_)
    status <- system2("bash", c(file.path(harness, "check-evidence.sh"), dir), stdout = FALSE,
                      stderr = FALSE)
    as.integer(status)
}
