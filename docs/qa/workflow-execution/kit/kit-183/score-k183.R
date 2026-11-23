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

scoreSimulation(jobs_dir, out)
