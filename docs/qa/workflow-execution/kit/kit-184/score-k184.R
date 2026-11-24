# KIT-183: score the simulation once (ticket 6.3 to 6.5; declaration.md). Base R only.
#
# Usage: Rscript score-k183.R <jobs dir> <out dir>
#   <jobs dir>  the evidence jobs' artifacts (k183-s1 ... k183-s8, k183-det), each holding w*/<table>.tsv
#   <out dir>   written: the tables below, results-summary.md (the outcome, the deciding numbers, the placeholders)

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
jobs_dir <- arguments[[1L]]
out <- arguments[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# KIT-183's criteria, as data (ticket 6.4): a later ticket (KIT-184) scores its own cells by passing another list.
CRITERIA <- list(criterion = 0.075, cells_total = 800L, views = c("transcriptomics", "proteomics", "metabolomics_LCMS"),
    rhos = c("0.1", "0.3", "0.6"),
    e1 = list(test = "camera", cells = c(12L, 20L), kinds = c("F_main", "F_ind"), classes = c("C0", "C1", "C2")),
    e2 = list(tests = c("gsea", "camera_pr"), cell = 12L, kinds = c("F_main", "F_ind"), classes = c("C0", "C2")),
    e3 = list(tests = c("camera", "gsea"), cell = 12L, kind = "F_main", accuracy = 0.99, minimum_called = 50L))

scoreSimulation <- function(jobs_dir, out, criteria = CRITERIA) {
    CRITERION <- criteria$criterion
    VIEWS <- criteria$views
    RHOS <- criteria$rhos
    readTable <- function(name, jobs) {
        files <- unlist(lapply(jobs, function(job) Sys.glob(file.path(jobs_dir, job, "w*", paste0(name, ".tsv")))))
        if (!length(files)) return(NULL)
        tables <- lapply(files, function(f) {
            header <- names(utils::read.delim(f, nrows = 1L, check.names = FALSE))
            classes <- stats::setNames(ifelse(header %in% c("rho_key", "message", "warnings"), "character", NA_character_), header)
            t <- utils::read.delim(f, stringsAsFactors = FALSE, colClasses = classes, na.strings = "NA")
            if ("rho_key" %in% names(t)) t$rho_key[is.na(t$rho_key)] <- ""
            t
        })
        out <- do.call(rbind, tables)
        rownames(out) <- NULL
        out
    }
    write <- function(table, name) utils::write.table(table, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE,
        row.names = FALSE)

    all_jobs <- basename(Sys.glob(file.path(jobs_dir, "k183-*")))
    main <- setdiff(all_jobs, "k183-det")
    status <- readTable("status", main)
    errors <- readTable("errors", main)
    power <- readTable("power", main)
    familywise <- readTable("familywise", main)
    diagnostic <- readTable("diagnostic", main)
    formed <- readTable("formed", main)
    counts <- readTable("counts", main)
    losses <- readTable("losses", main)
    factors <- readTable("factors", main)
    orientation <- readTable("orientation", main)
    errors$rho_key[is.na(errors$rho_key)] <- ""

    # --- Completion per cell (6.5: more than 5 % not completed: the cell is not decided) ---
    completion <- do.call(rbind, lapply(split(status, status$n), function(d) data.frame(n = d$n[[1L]], replicates = nrow(d),
        completed = sum(d$status == "completed"), did_not_complete = sum(d$status == "did_not_complete"),
        not_started = sum(d$status == "not_started"), median_total_seconds = stats::median(d$total_seconds, na.rm = TRUE),
        median_train_seconds = stats::median(d$train_seconds, na.rm = TRUE))))
    completion$decided <- (completion$completed / criteria$cells_total) >= 0.95
    write(completion, "completion")

    # --- Error rates (6.3): s_r per replicate, T over the replicates with a set, the t interval ---
    errors$s <- errors$sig / errors$sets
    rate <- do.call(rbind, lapply(split(errors, list(errors$n, errors$view, errors$kind, errors$test, errors$class,
        errors$rho_key), drop = TRUE), function(d) {
        R <- nrow(d)
        T <- mean(d$s)
        se <- if (R > 1L) stats::sd(d$s) / sqrt(R) else NA_real_
        q <- if (R > 1L) stats::qt(0.975, R - 1L) else NA_real_
        data.frame(n = d$n[[1L]], view = d$view[[1L]], kind = d$kind[[1L]], test = d$test[[1L]], class = d$class[[1L]],
            rho = d$rho_key[[1L]], replicates = R, sets = sum(d$sets), T = T, mc_se = se, lower = T - q * se, upper = T + q * se,
            stringsAsFactors = FALSE)
    }))
    rate <- rate[order(rate$n, rate$test, rate$kind, rate$view, rate$class, rate$rho), ]
    write(rate, "error-rates")

    pick <- function(n, test, kinds, classes) {
        d <- rate[rate$n == n & rate$test == test & rate$kind %in% kinds & rate$class %in% classes, , drop = FALSE]
        d
    }
    deciding <- function(n, test, classes) {
        # Every deciding check must be present: per kind (F_main, F_ind), class (C1 per rho) and data set.
        wanted <- expand.grid(kind = criteria$e1$kinds, view = VIEWS, class = classes, stringsAsFactors = FALSE)
        wanted <- do.call(rbind, lapply(seq_len(nrow(wanted)), function(i) {
            w <- wanted[i, ]
            if (w$class == "C1") data.frame(w, rho = RHOS, row.names = NULL) else data.frame(w, rho = "", row.names = NULL)
        }))
        d <- merge(wanted, rate[rate$n == n & rate$test == test, ], by = c("kind", "view", "class", "rho"), all.x = TRUE)
        d$pass <- !is.na(d$upper) & d$upper <= CRITERION
        d
    }

    # --- E1, E2 (6.4) ---
    e1 <- lapply(criteria$e1$cells, function(n) deciding(n, criteria$e1$test, criteria$e1$classes))
    names(e1) <- criteria$e1$cells
    e2_gsea <- deciding(criteria$e2$cell, criteria$e2$tests[[1L]], criteria$e2$classes)
    e2_cpr <- deciding(criteria$e2$cell, criteria$e2$tests[[2L]], criteria$e2$classes)
    write(do.call(rbind, lapply(names(e1), function(n) cbind(cell = n, e1[[n]]))), "e1-checks")
    write(rbind(cbind(test = "gsea", e2_gsea), cbind(test = "camera_pr", e2_cpr)), "e2-checks")

    # Operating characteristic of a cell, from the measured standard errors (6.3).
    oc <- function(d, true) {
        q <- stats::qt(0.975, d$replicates - 1L)
        prod(stats::pnorm((CRITERION - q * d$mc_se - true) / d$mc_se), na.rm = FALSE)
    }

    # --- E3: direction (scorer's own orientation from the job's independent sign and Z1's correlation) ---
    power$key <- paste(power$cell, power$r)
    orientation$key <- paste(orientation$cell, orientation$r)
    power$expected_up_if_planted_up <- sign(orientation$sign_main * orientation$raw_cor_main_z1)[match(power$key,
        orientation$key)] > 0
    power$expected <- ifelse((power$class == "Pup") == power$expected_up_if_planted_up, "Up", "Down")
    power$correct <- ifelse(power$expected == "Up", power$called_up, power$called_down)
    main_power <- power[power$kind == criteria$e3$kind, , drop = FALSE]
    e3 <- do.call(rbind, lapply(criteria$e3$tests, function(test) {
        d <- main_power[main_power$n == criteria$e3$cell & main_power$test == test, , drop = FALSE]
        called <- sum(d$called)
        data.frame(test = test, called = called, correct = sum(d$correct), accuracy = if (called) sum(d$correct) / called else NA,
            decided = called >= criteria$e3$minimum_called, pass = called < criteria$e3$minimum_called || sum(d$correct) / called >= criteria$e3$accuracy, stringsAsFactors = FALSE)
    }))
    write(e3, "e3-direction")
    power_table <- do.call(rbind, lapply(split(main_power, list(main_power$n, main_power$test, main_power$view,
        main_power$class, main_power$nominal), drop = TRUE), function(d) data.frame(n = d$n[[1L]], test = d$test[[1L]],
        view = d$view[[1L]], class = d$class[[1L]], size = d$nominal[[1L]], sets = sum(d$sets), power = sum(d$correct) / sum(d$sets),
        stringsAsFactors = FALSE)))
    write(power_table, "power")

    # --- Reported: family-wise, diagnostic, factors, losses, pools, completion ---
    fw <- do.call(rbind, lapply(split(familywise, list(familywise$n, familywise$view, familywise$kind, familywise$test),
        drop = TRUE), function(d) {
        x <- sum(d$any_bh)
        ci <- stats::binom.test(x, nrow(d))$conf.int
        data.frame(n = d$n[[1L]], view = d$view[[1L]], kind = d$kind[[1L]], test = d$test[[1L]], replicates = nrow(d),
            familywise = x / nrow(d), lower = ci[[1L]], upper = ci[[2L]], stringsAsFactors = FALSE)
    }))
    write(fw, "familywise")
    diagnostic$rho_key[is.na(diagnostic$rho_key)] <- ""
    write(stats::aggregate(ratio ~ n + view + kind + class + rho_key, data = diagnostic, FUN = mean), "diagnostic")
    strength <- do.call(rbind, lapply(split(factors, list(factors$n, factors$kind, factors$view), drop = TRUE), function(d)
        data.frame(n = d$n[[1L]], kind = d$kind[[1L]], view = d$view[[1L]], median_r2 = stats::median(d$r2_percent),
            min_r2 = min(d$r2_percent), max_r2 = max(d$r2_percent), median_largest = stats::median(d$largest_r2_percent),
            min_largest = min(d$largest_r2_percent), max_largest = max(d$largest_r2_percent),
            median_cor_z1 = stats::median(d$cor_z1), median_cor_z2 = stats::median(d$cor_z2), stringsAsFactors = FALSE)))
    write(strength, "factor-strength")
    write(stats::aggregate(cbind(zero_lost, loaded_lost, missing_share) ~ n + view, data = losses, FUN = mean), "censoring-losses")
    pools <- stats::aggregate(cbind(wanted, formed) ~ n + view + kind + class + size, data = formed, FUN = sum)
    pools$share_formed <- pools$formed / pools$wanted
    write(pools, "sets-formed")
    write(stats::aggregate(cbind(complete, incomplete, collapsed, universe) ~ n + view, data = counts, FUN = mean), "counts")
    write(as.data.frame(table(orientation$warnings)), "mofa-warnings")

    # --- Determinism (6.3): replicates 1 to 10 of n = 12 in the separate job against the same replicates here ---
    det <- NULL
    if ("k183-det" %in% all_jobs) {
        compare <- function(name) {
            a <- readTable(name, "k183-det")
            b <- readTable(name, main)
            if (is.null(a) || is.null(b)) return(NA)
            b <- b[b$n == 12L & b$r %in% a$r, , drop = FALSE]
            key <- function(d) d[do.call(order, unname(as.list(d))), , drop = FALSE]
            a <- key(a)
            b <- key(b[names(a)])
            rownames(a) <- rownames(b) <- NULL
            identical(a, b)
        }
        det <- data.frame(table = c("errors", "power", "familywise", "diagnostic", "formed", "counts", "losses", "factors",
            "orientation"), stringsAsFactors = FALSE)
        det$identical <- vapply(det$table, compare, logical(1))
        write(det, "determinism")
    }

    # --- Outcome (6.5) ---
    cell_ok <- function(n) isTRUE(completion$decided[completion$n == n])
    e1_hold <- vapply(names(e1), function(n) all(e1[[n]]$pass), logical(1))
    e2_gsea_hold <- all(e2_gsea$pass)
    e2_cpr_hold <- all(e2_cpr$pass)
    e3_hold <- all(e3$pass)
    decided <- all(vapply(c(12L, 20L), cell_ok, logical(1)))
    outcome <- if (!decided) "not decided (a deciding cell lost more than 5 % of its replicates): to the owner" else
        if (!e3_hold) "O5" else if (!e2_gsea_hold) "O4" else
        if (e1_hold[["12"]] && e1_hold[["20"]]) "O1" else
        if (!e1_hold[["12"]] && e1_hold[["20"]] && e2_cpr_hold) "O2a" else
        if (!e1_hold[["20"]] && e2_cpr_hold) "O2b" else "O3"

    # --- The numbers written into the texts (5.5) ---
    largest <- function(d) d[which.max(d$T), , drop = FALSE]
    fmt <- function(x) sprintf("%.1f %%", 100 * x)
    interval <- function(d) sprintf("%.1f to %.1f %%", 100 * d$lower, 100 * d$upper)
    r_cam <- largest(pick(12L, "camera", c("F_main", "F_ind"), c("C0", "C1", "C2")))
    s_cam <- largest(pick(12L, "camera", c("F_weak", "F_none"), c("C4", "C4b", "C1out")))
    gsea_03 <- largest(rate[rate$n == 12L & rate$test == "gsea" & rate$kind == "F_main" & rate$class == "C1" & rate$rho == "0.3", ])
    cpr_03 <- largest(rate[rate$n == 12L & rate$test == "camera_pr" & rate$kind == "F_main" & rate$class == "C1" &
        rate$rho == "0.3", ])
    r_cam_6 <- largest(pick(6L, "camera", c("F_main", "F_ind"), c("C0", "C1", "C2")))
    placeholders <- data.frame(name = c("R_cam", "I_cam", "S_cam", "R_gsea_0.3", "I_gsea_0.3", "R_cpr_0.3", "R_cam_6"),
        value = c(fmt(r_cam$T), interval(r_cam), fmt(s_cam$T), fmt(gsea_03$T), interval(gsea_03), fmt(cpr_03$T), fmt(r_cam_6$T)),
        source = c(rep(paste(r_cam$class, r_cam$rho, r_cam$view, r_cam$kind), 2L), paste(s_cam$class, s_cam$rho, s_cam$view,
            s_cam$kind), rep(paste(gsea_03$view, gsea_03$kind), 2L), paste(cpr_03$view, cpr_03$kind),
            paste(r_cam_6$class, r_cam_6$rho, r_cam_6$view, r_cam_6$kind)), stringsAsFactors = FALSE)
    write(placeholders, "placeholders")

    lines <- c("# KIT-183 simulation: outcome", "",
        sprintf("Outcome: **%s**", outcome), "",
        sprintf("E1 (camera, 30 checks per cell, upper limit <= %.3f): n = 12 %s (largest upper %.4f), n = 20 %s (largest upper %.4f)",
            CRITERION, if (e1_hold[["12"]]) "holds" else "fails", max(e1[["12"]]$upper, na.rm = TRUE),
            if (e1_hold[["20"]]) "holds" else "fails", max(e1[["20"]]$upper, na.rm = TRUE)),
        sprintf("E2 at n = 12: gsea %s (largest upper %.4f); camera_pr %s (largest upper %.4f)",
            if (e2_gsea_hold) "holds" else "fails", max(e2_gsea$upper, na.rm = TRUE), if (e2_cpr_hold) "holds" else "fails",
            max(e2_cpr$upper, na.rm = TRUE)),
        sprintf("E3 at n = 12: camera %d of %d called in the expected direction; gsea %d of %d", e3$correct[[1L]], e3$called[[1L]],
            e3$correct[[2L]], e3$called[[2L]]),
        sprintf("Operating characteristic (measured SEs), P(cell passes E1) if every check's true rate were 0.050 / 0.055: n = 12 %.4f / %.4f; n = 20 %.4f / %.4f",
            oc(e1[["12"]], 0.05), oc(e1[["12"]], 0.055), oc(e1[["20"]], 0.05), oc(e1[["20"]], 0.055)),
        "", "Placeholders:", sprintf("- %s = %s (%s)", placeholders$name, placeholders$value, placeholders$source), "",
        "Completion:", sprintf("- n = %d: %d of %d completed, %d did not complete, %d not started", completion$n,
            completion$completed, completion$replicates, completion$did_not_complete, completion$not_started),
        if (!is.null(det)) c("", sprintf("Determinism: %s", if (all(det$identical %in% TRUE)) "identical" else
            paste("differs in", paste(det$table[!det$identical %in% TRUE], collapse = ", ")))))
    writeLines(lines, file.path(out, "results-summary.md"))
    cat(lines, sep = "\n")
}

# ======================================================================================================================
# KIT-184 block (ticket 6.3 to 6.5; declaration docs/qa/workflow-execution/kit/kit-184/declaration.md). Everything above
# is KIT-183's score-k183.R at its declared commit, unchanged but for its last line (the call), which is made below
# unless K184_SCORE=ruv. With K184_SCORE=ruv the jobs of KIT-184's cells are scored: Rscript score-k184.R <jobs dir>
# <out dir>, <jobs dir> holding the evidence jobs' artifacts (k184-*), each with w*/<table>.tsv.
# Environment, for the blind lifecycle extension L1b only: K184_BLIND (the summary is not printed),
# K184_DETERMINISM_SELFTEST (the determinism comparison repeated with one value altered), K184_L1_REFERENCE (the path of
# L1's error-rates.tsv: cell A's rates with every size summed are compared with it, TRUE or FALSE).
# ======================================================================================================================

CRITERIA184 <- list(criterion = 0.075, declared = c(A = 800L, B = 800L, D = 400L), deciding_cells = c("A", "B"),
    views = c("transcriptomics", "proteomics", "metabolomics_LCMS"), rhos = c("0.1", "0.3", "0.6"),
    test_of_record = "camera", direction_accuracy = 0.99, minimum_called = 50L, canonical_stop = 0.9, pool_share = 0.9,
    # Redefinitions after the lifecycle run L1 (declaration, "Changes forced by L1"), made in this scorer only: the
    # simulation draws and tests every set as the scripts committed before L1 do (declaration audit K184-D03).
    # Sizes left out of the reported rate of a reported class (C4), per data set and model:
    left_out = data.frame(view = c(rep("transcriptomics", 4L), rep("metabolomics_LCMS", 2L)), class = "C4",
        model = c("corrected", "default", "default", "default", "corrected", "default"),
        size = c(100L, 10L, 30L, 100L, 20L, 20L), stringsAsFactors = FALSE),
    # Decision O1 (the owner's, 2026-11-24, methods ledger VAL-150): every check, number and outcome that reads
    # proteomics C5 reads it in two compositions, sizes 10 and 30 only, and every size formed in the replicate (as the
    # ticket declared and as the scorer committed before L1 computes). Each entry names the sizes left out; a check holds
    # only when it holds in every composition listed, and a printed number takes the larger T. Size 100 alone is in
    # error-rates-by-size.tsv and decides nothing. (The options not taken: sizes 10 and 30 only,
    # `list(sizes_10_30 = 100L)`; as declared only, `list(declared = integer())`.)
    c5_proteomics = list(sizes_10_30 = 100L, declared = integer()),
    # The feature-level rate on the main factor is not defined where the model keeps every feature of the data set.
    feature_rate_undefined = c(view = "metabolomics_LCMS", kind = "F_main"),
    # The worker tables the determinism job must reproduce (`status` without its seconds).
    determinism_tables = c("errors", "power", "formed", "counts", "feature", "identity", "floor", "exact", "caught",
        "factors", "orientation", "status"))

scoreRuv <- function(jobs_dir, out, criteria = CRITERIA184) {
    limit <- criteria$criterion
    # The counted jobs' tables. The determinism job's files (`k184-det`) enter no rate; they are read again below.
    readTable <- function(name) {
        files <- Sys.glob(file.path(jobs_dir, "k184-*", "w*", paste0(name, ".tsv")))
        files <- files[!grepl("k184-det", files, fixed = TRUE)]
        if (!length(files)) return(NULL)
        tables <- lapply(files, function(f) {
            header <- names(utils::read.delim(f, nrows = 1L, check.names = FALSE))
            classes <- stats::setNames(ifelse(header %in% c("rho_key", "message", "warnings", "arm", "cell"), "character",
                NA_character_), header)
            t <- utils::read.delim(f, stringsAsFactors = FALSE, colClasses = classes, na.strings = "NA")
            if ("rho_key" %in% names(t)) t$rho_key[is.na(t$rho_key)] <- ""
            t
        })
        table <- do.call(rbind, tables)
        rownames(table) <- NULL
        table
    }
    write <- function(table, name) if (!is.null(table)) utils::write.table(table, file.path(out, paste0(name, ".tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE)
    bind <- function(parts) do.call(rbind, parts)
    interval <- function(x) {
        R <- length(x)
        se <- if (R > 1L) stats::sd(x) / sqrt(R) else NA_real_
        q <- if (R > 1L) stats::qt(0.975, R - 1L) else NA_real_
        data.frame(replicates = R, T = mean(x), mc_se = se, lower = mean(x) - q * se, upper = mean(x) + q * se)
    }
    status <- readTable("status")
    errors <- readTable("errors")
    power <- readTable("power")
    formed <- readTable("formed")
    counts <- readTable("counts")
    feature <- readTable("feature")
    identity <- readTable("identity")
    floor <- readTable("floor")
    exact <- readTable("exact")
    caught <- readTable("caught")
    orientation <- readTable("orientation")
    versions <- readTable("versions")
    status <- status[status$cell %in% names(criteria$declared), , drop = FALSE]

    # --- Completion per cell (6.3: more than 5 % not completed: the cell is not decided) ---
    completion <- bind(lapply(split(status, status$cell), function(d) data.frame(cell = d$cell[[1L]], replicates = nrow(d),
        completed = sum(d$status == "completed"), did_not_complete = sum(d$status == "did_not_complete"),
        not_started = sum(d$status == "not_started"), median_total_seconds = stats::median(d$total_seconds, na.rm = TRUE),
        largest_total_seconds = max(d$total_seconds, na.rm = TRUE),
        median_train_seconds = stats::median(d$train_seconds, na.rm = TRUE))))
    completion$declared <- criteria$declared[completion$cell]
    completion$decided <- completion$completed / completion$declared >= 0.95
    write(completion, "completion")
    write(status[status$status != "completed", , drop = FALSE], "not-completed")

    # --- Error rates (6.3): the share of a replicate's class sets with p <= 0.05, averaged over replicates ---
    # The simulation writes one row per replicate, data set, factor, arm, test, class, rho_block and set size (`nominal`).
    # Every declared estimand is over a replicate's sets of the class, all sizes together, so the sizes are summed first.
    # A size in `criteria$left_out` is left out of the reported rate of its class, and proteomics C5 is read in every
    # composition of `criteria$c5_proteomics`; `composition` is "single" for every other row.
    if (!"nominal" %in% names(errors)) stop("the errors table has no set size (`nominal`): written by an older simulation script")
    key <- c("cell", "view", "kind", "arm", "test", "class", "rho_key")
    perReplicate <- function(d) {
        if (!nrow(d)) return(NULL)
        group <- paste(d$cell, d$r, d$view, d$kind, d$arm, d$test, d$class, d$rho_key, sep = "\r")
        sums <- rowsum(cbind(sets = d$sets, sig = d$sig), group, reorder = FALSE)
        per <- d[!duplicated(group), c("cell", "r", key[-1L]), drop = FALSE]
        per$sets <- unname(sums[, "sets"])
        per$sig <- unname(sums[, "sig"])
        per$s <- per$sig / per$sets
        rownames(per) <- NULL
        per
    }
    rateOf <- function(per) {
        if (is.null(per)) return(NULL)
        bind(lapply(split(per, per[key], drop = TRUE), function(d) cbind(d[1L, key], sets = sum(d$sets), interval(d$s))))
    }
    left <- rep(FALSE, nrow(errors))
    for (i in seq_len(nrow(criteria$left_out))) {
        rule <- criteria$left_out[i, ]
        left <- left | (errors$view == rule$view & errors$class == rule$class &
            ifelse(errors$arm %in% c("D0", "D1"), "default", "corrected") == rule$model & errors$nominal == rule$size)
    }
    is_c5p <- errors$view == "proteomics" & errors$class == "C5"
    compositions <- c(list(single = NULL), criteria$c5_proteomics)
    scored <- lapply(names(compositions), function(name) {
        rows <- if (name == "single") !left & !is_c5p else is_c5p & !(errors$nominal %in% compositions[[name]])
        per <- perReplicate(errors[rows, , drop = FALSE])
        if (is.null(per)) return(NULL)
        per$composition <- name
        rate <- rateOf(per)
        rate$composition <- name
        list(per = per, rate = rate)
    })
    scored <- Filter(Negate(is.null), scored)
    per_replicate <- bind(lapply(scored, `[[`, "per"))
    rate <- bind(lapply(scored, `[[`, "rate"))
    rate <- rate[do.call(order, rate[c(key, "composition")]), ]
    write(rate, "error-rates")
    # Per set size, over the replicates in which the size formed (their number is `replicates`): reported, decides nothing.
    # A size that does not form in every replicate is conditional on its pool.
    by_size <- bind(lapply(split(errors, errors[c(key, "nominal")], drop = TRUE), function(d) cbind(d[1L, c(key, "nominal")],
        sets = sum(d$sets), interval(d$sig / d$sets))))
    write(by_size[do.call(order, by_size[c(key, "nominal")]), ], "error-rates-by-size")
    # The checks of a criterion: the test of record under one arm; C0, C1 (each rho) and C2 on the main factor and on the
    # independent covariate, and C5 on the independent covariate, in every data set (thirty-three per cell).
    wanted <- bind(lapply(criteria$views, function(view) rbind(
        expand.grid(view = view, kind = c("F_main", "F_ind"), class = c("C0", "C2"), rho_key = "", stringsAsFactors = FALSE),
        expand.grid(view = view, kind = c("F_main", "F_ind"), class = "C1", rho_key = criteria$rhos, stringsAsFactors = FALSE),
        data.frame(view = view, kind = "F_ind", class = "C5", rho_key = "", stringsAsFactors = FALSE))))
    # Every check with the composition it is read in: "single", and for proteomics C5 one row per composition of
    # `criteria$c5_proteomics`. A row the data cannot fill stays in the table with missing values and does not pass.
    withCompositions <- function(table) {
        c5p <- table$view == "proteomics" & table$class == "C5"
        single <- table[!c5p, , drop = FALSE]
        single$composition <- rep("single", nrow(single))
        rbind(single, bind(lapply(names(criteria$c5_proteomics), function(name) {
            d <- table[c5p, , drop = FALSE]
            d$composition <- rep(name, nrow(d))
            d
        })))
    }
    checks <- function(arm, table = wanted) bind(lapply(criteria$deciding_cells, function(cell) {
        d <- merge(cbind(cell = cell, withCompositions(table)), rate[rate$arm == arm & rate$test == criteria$test_of_record, ],
            by = c("cell", "view", "kind", "class", "rho_key", "composition"), all.x = TRUE)
        d$pass <- !is.na(d$upper) & d$upper <= limit
        d
    }))
    # A check holds when it holds in every composition it is read in (one for every class but proteomics C5).
    held <- function(d) stats::aggregate(d["pass"], d[c("cell", "view", "kind", "class", "rho_key")], all)
    ke1 <- checks("3")
    ke5 <- checks("2")
    ke4 <- checks("3", data.frame(view = criteria$views, kind = "F_main", class = "C5", rho_key = "", stringsAsFactors = FALSE))
    c5_arm2 <- checks("2", data.frame(view = criteria$views, kind = "F_main", class = "C5", rho_key = "",
        stringsAsFactors = FALSE))
    write(ke1, "ke1-checks")
    write(ke5, "ke5-checks")
    write(ke4, "ke4-checks")
    # Operating characteristic of the thirty-three checks of a cell from the measured standard errors (6.3). Every row
    # counts: a check read in two compositions enters twice, which can only understate the probability.
    oc <- function(d, true) prod(stats::pnorm((limit - stats::qt(0.975, d$replicates - 1L) * d$mc_se - true) / d$mc_se))
    operating <- bind(lapply(split(ke1, ke1$cell), function(d) data.frame(cell = d$cell[[1L]], checks = nrow(d),
        pass_at_0.050 = oc(d, 0.050), pass_at_0.055 = oc(d, 0.055), pass_at_0.060 = oc(d, 0.060))))
    write(operating, "operating-characteristic")

    # --- K-E3: direction accuracy of arm 3 on the main factor (the scorer's own expected direction) ---
    direction <- NULL
    power_table <- NULL
    if (!is.null(power)) {
        power$key <- paste(power$cell, power$r)
        orientation$key <- paste(orientation$cell, orientation$r)
        up <- sign(orientation$sign_main * orientation$raw_cor_main_z1)[match(power$key, orientation$key)] > 0
        power$expected <- ifelse((power$class == "Pup") == up, "Up", "Down")
        power$correct <- ifelse(power$expected == "Up", power$called_up, power$called_down)
        own <- power[power$kind == "F_main" & power$arm %in% c("1", "2", "3", "3o"), , drop = FALSE]
        direction <- bind(lapply(split(own, own[c("cell", "arm", "test")], drop = TRUE), function(d) {
            called <- sum(d$called)
            data.frame(cell = d$cell[[1L]], arm = d$arm[[1L]], test = d$test[[1L]], called = called, correct = sum(d$correct),
                accuracy = if (called) sum(d$correct) / called else NA_real_, decided = called >= criteria$minimum_called,
                pass = called < criteria$minimum_called || sum(d$correct) / called >= criteria$direction_accuracy)
        }))
        write(direction, "ke3-direction")
        # Power: planted sets called at BH <= 0.05 in the expected direction (the arms of the model trained on corrected
        # values, on its main factor); for D0 and D1 the sets called in either direction, on the default model's factors.
        power$found <- ifelse(power$arm %in% c("D0", "D1"), power$called, power$correct)
        pkey <- c("cell", "kind", "arm", "test", "view", "class", "nominal")
        power_table <- bind(lapply(split(power, power[pkey], drop = TRUE), function(d) cbind(d[1L, pkey],
            sets = sum(d$sets), power = sum(d$found) / sum(d$sets))))
        power_table <- power_table[do.call(order, power_table[pkey]), ]
        write(power_table, "power")
    }

    # --- Paired differences on the same replicates and sets (6.3): arm 1 minus arm 3, arm 2 minus arm 3 ---
    own <- per_replicate[per_replicate$arm %in% c("1", "2", "3", "3o") & per_replicate$test == criteria$test_of_record, ,
        drop = FALSE]
    wide <- stats::reshape(own[c("cell", "r", "view", "kind", "class", "rho_key", "composition", "arm", "s")],
        idvar = c("cell", "r", "view", "kind", "class", "rho_key", "composition"), timevar = "arm", direction = "wide")
    paired <- bind(lapply(c("1", "2", "3o"), function(arm) {
        column <- paste0("s.", arm)
        if (!column %in% names(wide) || !"s.3" %in% names(wide)) return(NULL)
        wide$difference <- wide[[column]] - wide[["s.3"]]
        d <- wide[!is.na(wide$difference), , drop = FALSE]
        pkey <- c("cell", "view", "kind", "class", "rho_key", "composition")
        bind(lapply(split(d, d[pkey], drop = TRUE), function(x) cbind(x[1L, pkey], comparison = paste0("arm ", arm,
            " minus arm 3"), interval(x$difference))))
    }))
    write(paired, "paired-differences")

    # --- Feature-level rates (6.3) ---
    feature_rates <- NULL
    if (!is.null(feature)) {
        feature$s <- feature$significant / feature$features
        fkey <- c("cell", "view", "kind", "arm")
        feature_rates <- bind(lapply(split(feature, feature[fkey], drop = TRUE), function(d) cbind(d[1L, fkey],
            features = sum(d$features), interval(d$s[is.finite(d$s)]))))
        # Not defined where the model keeps every feature (no feature is left that it did not select): stated, not NaN.
        undefined <- feature_rates$view == criteria$feature_rate_undefined[["view"]] &
            feature_rates$kind == criteria$feature_rate_undefined[["kind"]]
        feature_rates[undefined, c("replicates", "T", "mc_se", "lower", "upper")] <- NA
        feature_rates$defined <- !undefined
        write(feature_rates[do.call(order, feature_rates[fkey]), ], "feature-rates")
    }

    # --- Exactness inside the shipped path; what the removed factors caught; pools; the floor ---
    exact_table <- bind(lapply(split(exact, exact[c("cell", "view")], drop = TRUE), function(d) data.frame(cell = d$cell[[1L]],
        view = d$view[[1L]], replicates = nrow(d), exact = sum(d$status == "exact"), share_exact = mean(d$status == "exact"),
        largest_ratio = max(d$ratio), smallest_share = min(d$share), largest_relative_residual = max(d$largest_relative_residual),
        median_features = stats::median(d$features), median_rows = stats::median(d$rows))))
    write(exact_table, "exactness")
    identity_table <- bind(lapply(split(identity, identity[c("cell", "view", "kind")], drop = TRUE), function(d) data.frame(
        cell = d$cell[[1L]], view = d$view[[1L]], kind = d$kind[[1L]], replicates = nrow(d),
        largest_correlation_difference = max(d$max_correlation_difference),
        largest_pvalue_difference = max(d$max_pvalue_difference),
        replicates_above_1e_8 = sum(d$max_correlation_difference > 1e-8))))
    write(identity_table, "identity")
    main_counts <- counts[counts$model == "corrected" & counts$kind == "F_main", , drop = FALSE]
    caught_table <- bind(lapply(split(caught, caught[c("cell", "view")], drop = TRUE), function(d) {
        r2 <- main_counts$factor_r2[main_counts$cell == d$cell[[1L]] & main_counts$view == d$view[[1L]]]
        data.frame(cell = d$cell[[1L]], view = d$view[[1L]], replicates = nrow(d),
            median_smallest_canonical_correlation = stats::median(d$smallest_canonical_correlation),
            smallest = min(d$smallest_canonical_correlation), largest = max(d$smallest_canonical_correlation),
            below_stop = stats::median(d$smallest_canonical_correlation) < criteria$canonical_stop,
            median_factor_r2 = stats::median(r2), smallest_factor_r2 = min(r2), largest_factor_r2 = max(r2))
    }))
    write(caught_table, "caught")
    fkey <- c("cell", "model", "view", "kind", "class", "size")
    pools <- bind(lapply(split(formed, formed[fkey], drop = TRUE), function(d) cbind(d[1L, fkey], wanted = sum(d$wanted),
        formed = sum(d$formed))))
    pools$share <- pools$formed / pools$wanted
    pools$below <- pools$share < criteria$pool_share
    write(pools[do.call(order, pools[fkey]), ], "pools")
    floor_table <- if (!is.null(floor)) bind(lapply(split(floor, floor[c("cell", "view", "kind")], drop = TRUE), function(d)
        data.frame(cell = d$cell[[1L]], view = d$view[[1L]], kind = d$kind[[1L]], sets = sum(d$sets),
            largest_camera_minus_floor = max(d$largest_camera_minus_floor))))
    write(floor_table, "floor")
    write(unique(versions[c("package", "version")]), "versions")

    # --- Determinism (6.3; declaration audit K184-D04): the determinism job against the same replicates of the counted
    # jobs. Read as text, so that what is compared is what the workers wrote; ordered by every column; `identical()`. ---
    readText <- function(name, det, keep = NULL) {
        files <- Sys.glob(file.path(jobs_dir, "k184-*", "w*", paste0(name, ".tsv")))
        files <- files[grepl("k184-det", files, fixed = TRUE) == det]
        if (!length(files)) return(NULL)
        table <- do.call(rbind, lapply(files, function(f) {
            t <- utils::read.delim(f, colClasses = "character", check.names = FALSE, quote = "", na.strings = character())
            if (!is.null(keep)) t <- t[paste(t$cell, t$r) %in% keep, , drop = FALSE]
            t
        }))
        rownames(table) <- NULL
        table
    }
    sorted <- function(d) {
        d <- d[do.call(order, unname(as.list(d))), , drop = FALSE]
        rownames(d) <- NULL
        d
    }
    determinism <- NULL
    det_status <- readText("status", det = TRUE)
    if (!is.null(det_status)) {
        again <- paste(det_status$cell, det_status$r)
        pairOf <- function(name) {
            a <- readText(name, det = TRUE)
            b <- readText(name, det = FALSE, keep = again)
            if (name == "status") {
                seconds <- c("train_seconds", "total_seconds")
                if (!is.null(a)) a <- a[setdiff(names(a), seconds)]
                if (!is.null(b)) b <- b[setdiff(names(b), seconds)]
            }
            list(a = a, b = b)
        }
        compare <- function(label, p) data.frame(table = label, rows_determinism = if (is.null(p$a)) 0L else nrow(p$a),
            rows_counted = if (is.null(p$b)) 0L else nrow(p$b), identical = !is.null(p$a) && !is.null(p$b) &&
                identical(names(p$a), names(p$b)) && identical(sorted(p$a), sorted(p$b)), stringsAsFactors = FALSE)
        determinism <- bind(lapply(criteria$determinism_tables, function(name) compare(name, pairOf(name))))
        write(determinism, "determinism")
        # The comparison shown to bite (asked of the lifecycle extension L1b): the same comparison with one value of one
        # table altered must read FALSE.
        if (nzchar(Sys.getenv("K184_DETERMINISM_SELFTEST", unset = ""))) {
            p <- pairOf("errors")
            if (!is.null(p$a) && nrow(p$a)) p$a$sig[[1L]] <- paste0(p$a$sig[[1L]], "0")
            write(compare("errors, one value altered", p), "determinism-selftest")
        }
    }

    # --- The lifecycle extension's check that writing the set size moved no p-value (declaration audit K184-D03, D05):
    # the rates of cell A with every size summed and nothing left out, against the table L1's scorer wrote. ---
    reference_path <- Sys.getenv("K184_L1_REFERENCE", unset = "")
    if (nzchar(reference_path)) {
        reference <- utils::read.delim(reference_path, stringsAsFactors = FALSE,
            colClasses = c(cell = "character", arm = "character", rho_key = "character"))
        reference$rho_key[is.na(reference$rho_key)] <- ""
        reference <- reference[reference$cell == "A", , drop = FALSE]
        drawn <- rateOf(perReplicate(errors[errors$cell == "A", , drop = FALSE]))
        # Counts only, so that a FALSE can be read without opening a rate (declaration audit K184-R02): whether the
        # number of sets and of replicates agrees in every row, and in how many rows T differs. Equal sets with a few
        # differing T would point at a p-value on the other side of 0.05 on another processor, not at the summing.
        equal <- FALSE
        matched <- 0L
        sets_equal <- FALSE
        replicates_equal <- FALSE
        t_differing <- NA_integer_
        if (!is.null(drawn)) {
            both <- merge(drawn, reference, by = key, suffixes = c("_here", "_l1"))
            close <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE))
            matched <- nrow(both)
            all_matched <- nrow(drawn) == nrow(reference) && matched == nrow(reference)
            sets_equal <- all_matched && all(both$sets_here == both$sets_l1)
            replicates_equal <- all_matched && all(both$replicates_here == both$replicates_l1)
            t_differing <- sum(xor(is.na(both$T_here), is.na(both$T_l1)) |
                (!is.na(both$T_here) & !is.na(both$T_l1) & abs(both$T_here - both$T_l1) > 1e-12))
            equal <- sets_equal && replicates_equal && t_differing == 0L &&
                close(both$T_here, both$T_l1) && close(both$mc_se_here, both$mc_se_l1) &&
                close(both$lower_here, both$lower_l1) && close(both$upper_here, both$upper_l1)
        }
        write(data.frame(check = "cell A, every size summed and nothing left out, against L1's error rates",
            rows_here = if (is.null(drawn)) 0L else nrow(drawn), rows_l1 = nrow(reference), rows_matched = matched,
            sets_equal = sets_equal, replicates_equal = replicates_equal, rows_with_T_differing = t_differing, equal = equal,
            stringsAsFactors = FALSE), "l1-reference-check")
    }

    # --- The numbers of the texts (6.4) ---
    worst <- function(d) d[which.max(d$T), , drop = FALSE]
    r3 <- worst(ke1[!is.na(ke1$T), , drop = FALSE])
    r5 <- worst(ke4[!is.na(ke4$T), , drop = FALSE])
    r6 <- worst(c5_arm2[!is.na(c5_arm2$T), , drop = FALSE])
    r5m <- ke4[which(ke4$cell == r6$cell & ke4$view == r6$view & ke4$composition == r6$composition), , drop = FALSE]
    fmt <- function(d) if (nrow(d)) sprintf("%.1f %%", 100 * d$T) else NA_character_
    ival <- function(d) if (nrow(d)) sprintf("%.1f to %.1f %%", 100 * d$lower, 100 * d$upper) else NA_character_
    where <- function(d) if (nrow(d)) sprintf("cell %s, %s, %s, %s%s%s", d$cell, d$view, d$kind, d$class, ifelse(nzchar(d$rho_key),
        paste0(" rho ", d$rho_key), ""), ifelse(d$composition == "single", "", paste0(", composition ", d$composition))) else
        NA_character_
    powerRange <- function(test) {
        if (is.null(power_table)) return(NA_character_)
        d <- power_table[power_table$arm == "3" & power_table$test == test & power_table$kind == "F_main" &
            power_table$nominal == 10L & power_table$view %in% c("transcriptomics", "proteomics") &
            power_table$cell %in% criteria$deciding_cells, , drop = FALSE]
        pooled <- vapply(split(d, d[c("cell", "view")], drop = TRUE), function(x) sum(x$power * x$sets) / sum(x$sets), numeric(1))
        if (!length(pooled)) NA_character_ else sprintf("%.0f to %.0f %%", 100 * min(pooled), 100 * max(pooled))
    }
    numbers <- data.frame(placeholder = c("R3", "I3", "R5", "I5", "R6", "I6", "R5m", "PF", "PC"),
        value = c(fmt(r3), ival(r3), fmt(r5), ival(r5), fmt(r6), ival(r6), fmt(r5m), powerRange("camera_floor"),
            powerRange("camera")),
        from = c(where(r3), where(r3), where(r5), where(r5), where(r6), where(r6), where(r5m),
            "arm 3, camera_floor, sets of 10, main factor, transcriptomics and proteomics, cells A and B",
            "arm 3, camera, sets of 10, main factor, transcriptomics and proteomics, cells A and B"), stringsAsFactors = FALSE)
    write(numbers, "numbers")

    # --- The outcome (6.5) ---
    cells_decided <- all(completion$decided[completion$cell %in% criteria$deciding_cells])
    ke3 <- if (is.null(direction)) NULL else direction[direction$arm == "3" & direction$test == criteria$test_of_record &
        direction$cell %in% criteria$deciding_cells, , drop = FALSE]
    countLine <- function(d) paste(vapply(split(d, d$cell), function(x) {
        h <- held(x)
        sprintf("cell %s %d of %d hold, largest upper limit %.4f", x$cell[[1L]], sum(h$pass), nrow(h), max(x$upper, na.rm = TRUE))
    }, character(1)), collapse = "; ")
    lines <- c("# KIT-184 simulation: the scorer's summary", "",
        sprintf("- Cells decided (at least 95 %% of the declared replicates completed): %s", paste(sprintf("%s %s (%d of %d)",
            completion$cell, completion$decided, completion$completed, completion$declared), collapse = "; ")),
        sprintf("- K-E1 (arm 3, %s; 33 checks per cell, upper limit at most %.3f): %s", criteria$test_of_record, limit,
            countLine(ke1)),
        sprintf("- K-E3 (direction accuracy of arm 3 on the main factor, at least %.2f when at least %d planted sets were called): %s",
            criteria$direction_accuracy, criteria$minimum_called, if (is.null(ke3)) "no planted set scored" else
                paste(sprintf("cell %s %d called, %d correct, pass %s", ke3$cell, ke3$called, ke3$correct, ke3$pass),
                    collapse = "; ")),
        sprintf("- K-E4 (declared; arm 3, C5 on the main factor): %s", paste(sprintf("cell %s %s (%s) T %.4f upper %.4f pass %s",
            ke4$cell, ke4$view, ke4$composition, ke4$T, ke4$upper, ke4$pass), collapse = "; ")),
        sprintf("- K-E5 (arm 2, the form without the adjustment; the same 33 checks): %s", countLine(ke5)),
        sprintf("- Written factors exact: %s", paste(sprintf("cell %s %s %d of %d", exact_table$cell, exact_table$view,
            exact_table$exact, exact_table$replicates), collapse = "; ")),
        sprintf("- Identity check (largest difference of arm 3's correlation from the same model on the corrected values): %s",
            paste(sprintf("cell %s %s %s %.3g", identity_table$cell, identity_table$view, identity_table$kind,
                identity_table$largest_correlation_difference), collapse = "; ")),
        sprintf("- What the removed factors caught (median smallest canonical correlation with the planted factors; stop below %.1f): %s",
            criteria$canonical_stop, paste(sprintf("cell %s %s %.3f", caught_table$cell, caught_table$view,
                caught_table$median_smallest_canonical_correlation), collapse = "; ")),
        sprintf("- Pools: %d of %d class-and-size rows formed fewer than %.0f %% of their declared sets", sum(pools$below),
            nrow(pools), 100 * criteria$pool_share),
        sprintf("- Seconds per replicate (median; largest): %s", paste(sprintf("cell %s %.0f; %.0f", completion$cell,
            completion$median_total_seconds, completion$largest_total_seconds), collapse = "; ")),
        sprintf("- Determinism (the determinism job against the same replicates of the counted jobs, %d tables): %s",
            length(criteria$determinism_tables), if (is.null(determinism)) "no determinism job among the artifacts" else
                sprintf("%d of %d identical%s", sum(determinism$identical), nrow(determinism), if (all(determinism$identical)) "" else
                    paste0("; different: ", paste(determinism$table[!determinism$identical], collapse = ", ")))),
        sprintf("- Outcome rows of 6.5: R1 holds (K-E1 in A and B, K-E3, cells decided) %s; R7 (K-E4 fails in a data set) %s; R8 (K-E5 fails in a cell) %s",
            cells_decided && all(ke1$pass) && (is.null(ke3) || all(ke3$pass)), any(!ke4$pass), any(!ke5$pass)),
        sprintf("- Outcome row R6, read from completion.tsv, exactness.tsv and identity.tsv: a cell not decided %s; the written factors not exact in more than 1 %% of a data set's replicates %s; the identity check above 1e-8 in any replicate %s",
            !all(completion$decided), any(exact_table$share_exact < 0.99), any(identity_table$replicates_above_1e_8 > 0)))
    writeLines(lines, file.path(out, "results-summary.md"))
    # A blind lifecycle run (K184_BLIND set) prints nothing of the result; its job decides which files leave it.
    if (nzchar(Sys.getenv("K184_BLIND", unset = ""))) cat("scored blind: the summary is not printed\n") else cat(lines, sep = "\n")
}

if (identical(Sys.getenv("K184_SCORE", unset = ""), "ruv")) scoreRuv(jobs_dir, out) else scoreSimulation(jobs_dir, out)
