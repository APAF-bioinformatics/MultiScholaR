# KIT-129 final method (decisions M29 and M30): the package's new-design kernels against the exactness round's arms,
# on the simulators' own data sets. A recorded check, not a test: the simulators live on run-only branches.
#
# Usage (from the work tree's root; R_LIBS with sva, limpa and fastlimpadpc; one BLAS/OpenMP thread):
#   Rscript docs/qa/workflow-execution/kit/kit-129/final-method/prove-arms.R <stage> <out dir>
# <stage> is the commit being proven:
#   m29   arm 25 (each comparison's standard error exact): classic proteomics (S1), DIA with limpa (S2),
#         lipidomics and metabolomics (S4)
#   m30a  arms 31 and 32's count (200 permutations; the biological units removed on the technical-duplicate cell):
#         the package's estimated count (S1: `estimateRuvFactorCount`; S2: the declared kernel's own count)
#   m30b  arm 26 (classic proteomics: treat() given the trend and robust settings) and arm 32 (all together) in every
#         omic with arm 32
# Writes <out dir>/<stage>-simulator.rds and <stage>-package.rds (named lists: per data set, arm and omic, the table
# keyed by feature and contrast, and k), compares them with packaging/verify/compare-results.R (relative 1e-6,
# absolute floor 1e-9, discrete values exact) into <stage>-compare.tsv, and writes <stage>-summary.tsv (one row per
# data set, arm and omic: the largest relative difference of the doubles and the verdict).
#
# The simulators are read with `git show` at the commits that launched the round: s1_classic.R and s2_limpa.R at
# sim-k124-exact 2ba5cf87, s4_lipid_metab.R at sim-k129-exact 2d03b103 (the kit-128 tree, where the extract is declared
# as the technical replicate). Only their definitions are evaluated (assignments before their driver; nothing that
# reads the command line or the environment). The data sets are each simulator's `simulate`/`simulateData` +
# `quantify`/`generateData` at the round's seed rule (100000 + 10000 x cell index + run), run 1 of cells C0 and C6.
# The arms' tables are the simulators' own test functions (`da`, `da_variant`; `routeTest`; `testWithW`) with E, W and
# the arm's changes formed as each simulator's `one`/`runOne` forms them (the lines are cited below).
#
# The package side is the code under proof: classic proteomics `computeProtRuvModelFactors` +
# `attachProtRuvModelFactors` + `differentialAbundanceAnalysis` (k as given: the classic route does not call the
# count; the count is `estimateRuvFactorCount`, called here); DIA with limpa the declared kernel
# (`diaRuvCovariatesOperation("declared")`, all proteins as controls, its own count) + `differentialAbundanceAnalysis`;
# lipidomics and metabolomics `runOmicsRuvFactors` + `runLipidsDA` / `runMetabolitesDA` with the simulator's W in the
# assay records (S4 computes W on the features complete in study and pooled-QC injections, which the package does not
# hold; the proof is of the DA kernel, the part M29 changes).

stage <- commandArgs(TRUE)[[1L]]
out_dir <- commandArgs(TRUE)[[2L]]
stopifnot(stage %in% c("m29", "m30a", "m30b"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
logger::log_threshold(logger::FATAL)

silently <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

# The simulator's definitions, from the branch commit, in an environment of their own.
simulator <- function(file, commit, driver) {
    text <- system2("git", c("show", paste0(commit, ":docs/qa/workflow-execution/kit/kit-124/sim/", file)), stdout = TRUE)
    stop_at <- grep(driver, text)[[1L]]
    exprs <- parse(text = text[seq_len(stop_at - 1L)], keep.source = FALSE)
    env <- new.env(parent = globalenv())
    for (e in exprs) {
        if (!is.call(e) || !identical(e[[1L]], as.name("<-"))) next
        is_function <- is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))
        if (!is_function && grepl("args|commandArgs|Sys.getenv|out_dir|load_all|arms_wanted|cells_wanted|first_run", paste(deparse(e), collapse = " "))) next
        eval(e, env)
    }
    env
}

# A table keyed by feature (and contrast), sorted, with the compared columns.
keyed <- function(feature, contrast, estimate, p, q) {
    t <- data.frame(key = paste(feature, contrast, sep = "|"), estimate = as.double(estimate), p = as.double(p),
        q = as.double(q), stringsAsFactors = FALSE)
    t <- t[order(t$key), , drop = FALSE]
    rownames(t) <- NULL
    t
}

sim_side <- list()
pkg_side <- list()
notes <- list()
add <- function(name, simulator_value, package_value, note = "") {
    sim_side[[name]] <<- simulator_value
    pkg_side[[name]] <<- package_value
    notes[[name]] <<- note
}
# Tables restricted to their common keys; the keys only one side holds are counted in the note.
add_tables <- function(name, simulator_table, package_table) {
    common <- intersect(simulator_table$key, package_table$key)
    note <- sprintf("keys: %d common, %d simulator only, %d package only", length(common),
        length(setdiff(simulator_table$key, common)), length(setdiff(package_table$key, common)))
    pick <- function(t) {
        t <- t[t$key %in% common, , drop = FALSE]
        rownames(t) <- NULL
        t
    }
    add(name, pick(simulator_table), pick(package_table), note)
}

cells <- c("C0", "C6")

# ---- S1, classic proteomics ------------------------------------------------------------------------------------------
s1 <- simulator("s1_classic.R", "2ba5cf87", "^# ---- the study")
s1_contrast <- "treated_vs_control=grouptreated-groupcontrol"
for (cell in cells) {
    p <- s1$cell_setting(cell)
    seed <- 100000L + 10000L * s1$cell_index(cell) + 1L
    sim <- s1$simulate(p, seed)
    x <- sim$x
    design <- sim$design
    xc <- x[sim$complete, , drop = FALSE]
    # E as `one` forms it (s1_classic.R lines 468-488): stream 0, the group model or the units (`tech_rep_group`)
    count <- function(B, units) {
        model <- if (units) stats::model.matrix(~ tech_rep_group, design) else stats::model.matrix(~ group, design)
        as.integer(silently(sva::num.sv(xc, model, method = "be", B = B, seed = seed * 100L)))
    }
    # W as `in_model(all_controls, k, "no_tech")` forms it (lines 538-553): every complete protein a control, the groups
    # as replicate sets
    factors <- function(k) silently(computeRuvCovariates(s1$make_object(xc, design, "no_tech"), k = k,
        controls = rep(TRUE, nrow(xc)), ruv_grouping_variable = "group"))$W
    arm_table <- function(k, exact, applied) {
        set.seed(seed * 100L + 7L)
        tab <- silently(s1$da_variant(x, design, factors(k), exact = exact, applied = applied))
        keyed(tab$id, "contrast", tab$estimate, tab$p, tab$q)
    }
    # the package: the declared state on the normalised object (the true technical-replicate declaration)
    normalised <- s1$make_object(x, design, "tech_rep_group")
    state <- function(k) {
        factors <- computeProtRuvModelFactors(normalised, ruv_k = k, ruv_k_rule = "user", ruv_grouping_variable = "group")
        attachProtRuvModelFactors(normalised, normalised, factors)
    }
    entry <- function(k) {
        result <- silently(differentialAbundanceAnalysis(state(k), contrasts_tbl = data.frame(contrasts = s1_contrast,
            stringsAsFactors = FALSE), formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
            eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Group", qvalue_column = "fdr_qvalue",
            raw_pvalue_column = "raw_pvalue"))
        tab <- result$contrasts_results$results[[1L]]
        keyed(rownames(tab), "contrast", tab$logFC, tab$raw_pvalue, tab$fdr_qvalue)
    }
    units <- if (p$duplicates) design$tech_rep_group
    package_count <- function() estimateRuvFactorCount(x, design$group, seed = seed * 100L, units = units)$k
    E20 <- count(20L, FALSE)
    if (stage == "m29") {
        stopifnot(E20 > 0L)
        add(sprintf("S1/%s/arm25/k", cell), E20, E20, "k given to the package: the arm's E (groups, B = 20)")
        add_tables(sprintf("S1/%s/arm25/classic", cell), arm_table(E20, TRUE, FALSE), entry(E20))
    }
    if (stage == "m30a") {
        add(sprintf("S1/%s/arm31/k", cell), count(200L, FALSE), estimateRuvFactorCount(x, design$group, seed = seed * 100L)$k,
            "arm 31's E (groups, B = 200) against estimateRuvFactorCount's default")
        add(sprintf("S1/%s/arm32/k", cell), count(200L, p$duplicates), package_count(),
            "arm 32's E (B = 200; the units removed on C6) against estimateRuvFactorCount's default (units on C6)")
    }
    if (stage == "m30b") {
        stopifnot(E20 > 0L)
        ruv_model <- resolveProtRuvModelInput(state(E20), "~ 0 + group")
        design_rows <- state(E20)@design_matrix
        rownames(design_rows) <- design_rows$Run
        approximate <- silently(runProtRuvModelContrasts(ruv_model$matrix, s1_contrast, design_rows, ruv_model$formula_used,
            treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE, group_id = "group",
            technical_replicate_id = "tech_rep_group", exact = FALSE))$results[[1L]]
        add_tables(sprintf("S1/%s/arm26/classic", cell), arm_table(E20, FALSE, TRUE),
            keyed(rownames(approximate), "contrast", approximate$logFC, approximate$raw_pvalue, approximate$fdr_qvalue))
        E32 <- count(200L, p$duplicates)
        k32 <- package_count()
        add(sprintf("S1/%s/arm32/k", cell), E32, k32, "arm 32's E against the package's count")
        if (k32 > 0L) add_tables(sprintf("S1/%s/arm32/classic", cell), arm_table(E32, TRUE, TRUE), entry(k32))
    }
}

# ---- S2, DIA with limpa ------------------------------------------------------------------------------------------------
if (stage %in% c("m29", "m30a", "m30b")) {
    s2 <- simulator("s2_limpa.R", "2ba5cf87", "^if \\(!is.null\\(cells_wanted\\)\\)")
    table_s2 <- s2$cellTable()
    for (cell_name in cells) {
        cell <- table_s2[[cell_name]]
        seed <- 100000L + 10000L * cell$index + 1L
        sim <- s2$simulateData(cell, seed)
        quantified <- s2$quantify(sim)
        elist <- quantified@args$limpa_dpc_quant_results$quantified_elist
        design <- sim$design
        rownames(design) <- design$Run
        # E as `runOne` forms it (s2_limpa.R lines 479-481 and 544-554): seed 100 x seed + 99
        count <- function(B, units) {
            model <- if (units) stats::model.matrix(~ tech_rep_group, design) else stats::model.matrix(~ group, design)
            as.integer(silently(sva::num.sv(elist$E, model, method = "be", B = B, seed = 100L * seed + 99L)))
        }
        # `inModel(k, all_controls, exact)` (lines 531-540)
        arm_table <- function(k) {
            set.seed(100L * seed + 7L)
            W <- silently(computeRuvCovariates(quantified, k, rep(TRUE, nrow(elist$E)), "group"))$W
            tab <- silently(s2$routeTest(elist, design, W, exact = TRUE))
            keyed(tab$id, "contrast", tab$logFC, tab$p, tab$q)
        }
        # the package: the declared kernel (all proteins as controls), then the DA entry point on its final state
        kernel <- function(...) {
            saved <- list()
            store <- list(checkpoint = function(value, label, from = NULL) {
                saved[[label]] <<- value
                label
            })
            parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features",
                ruv_k_seed = 100L * seed + 99L, ...)
            value <- silently(diaRuvCovariatesOperation("declared")$kernel(quantified, parameters, store))
            list(value = value, state = saved[[value$final]])
        }
        entry <- function(state) {
            result <- silently(differentialAbundanceAnalysis(state, contrasts_tbl = data.frame(contrasts = "B_vs_A=groupB-groupA",
                stringsAsFactors = FALSE), formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
                eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue",
                raw_pvalue_column = "raw_pvalue"))
            tab <- result$contrasts_results$results[[1L]]
            keyed(rownames(tab), "contrast", tab$logFC, tab$raw_pvalue, tab$fdr_qvalue)
        }
        units <- cell$units > 0L
        if (stage == "m29") {
            E20 <- count(20L, FALSE)
            stopifnot(E20 > 0L)
            # arm 25's k is E with the groups removed; on C6 the package's count removes the units (M28), so the arm's k
            # is given as the user's k there
            package <- if (units) kernel(ruv_k_rule = "user", ruv_k = E20) else kernel()
            add(sprintf("S2/%s/arm25/k", cell_name), E20, package$value$best_k,
                if (units) "k given to the package: the arm's E" else "the declared kernel's own count (B = 20 at this commit)")
            add_tables(sprintf("S2/%s/arm25/dia", cell_name), arm_table(E20), entry(package$state))
        }
        if (stage == "m30a") {
            package <- kernel()
            add(sprintf("S2/%s/arm%s/k", cell_name, if (units) "32" else "31"), count(200L, units), package$value$ruv_k_estimate,
                "arm 31's E (C0) or arm 32's E (C6, the units removed), B = 200, against the declared kernel's own count")
        }
        if (stage == "m30b") {
            E32 <- count(200L, units)
            package <- kernel()
            add(sprintf("S2/%s/arm32/k", cell_name), E32, package$value$best_k, "arm 32's E against the declared kernel's own count")
            if (E32 > 0L) add_tables(sprintf("S2/%s/arm32/dia", cell_name), arm_table(E32), entry(package$state))
        }
    }
}

# ---- S4, lipidomics and metabolomics ---------------------------------------------------------------------------------
if (stage %in% c("m29", "m30b")) {
    s4 <- simulator("s4_lipid_metab.R", "2d03b103", "^# ---- driver")
    table_s4 <- s4$cellTable()
    for (cell_name in cells) {
        cell <- table_s4[[cell_name]]
        seed <- 100000L + 10000L * cell$index + 1L
        data <- s4$generateData(cell, seed)
        # E and its variants as `runOne` forms them (s4_lipid_metab.R lines 836-852)
        estimate <- function(B, units) lapply(data$assays, function(a) {
            dat <- a$y[s4$completeFeatures(a), , drop = FALSE]
            if (nrow(dat) < 20L) return(NA_integer_)
            mod <- stats::model.matrix(~ factor(if (units) data$design$extract else data$design$group))
            as.integer(silently(sva::num.sv(dat, mod, method = "be", B = B, seed = seed * 100L + 99L)))
        })
        view <- s4$makeView(data, estimate(20L, FALSE))
        k_of <- if (stage == "m29") estimate(20L, FALSE) else estimate(200L, cell$duplicates)
        arm <- if (stage == "m29") "arm25" else "arm32"
        # `modelArm(view, allComplete, k)` (lines 572-583): W on the complete features, the groups as replicate sets
        W <- lapply(stats::setNames(s4$ASSAYS, s4$ASSAYS), function(assay) {
            a <- view$assays[[assay]]
            y <- s4$completeY(a)
            s4$ruvW(y, s4$groupM(view$design), rep(TRUE, ncol(y)), k_of[[assay]], keep = view$design$Sample_ID,
                residual_df = s4$residualDf(view))
        })
        stopifnot(!any(vapply(W, s4$isRefusal, logical(1L))))
        mats <- lapply(view$assays, `[[`, "y")
        for (omic in c("lipid", "metab")) {
            object <- if (omic == "lipid") s4$lipidObject(mats, view$design) else s4$metabObject(mats, view$design)
            run <- silently(runOmicsRuvFactors(object, "groups_all_features", lapply(k_of, as.integer), "user"))
            for (assay in s4$ASSAYS) {
                if (is.null(W[[assay]])) {
                    run$factors[[assay]]["W"] <- list(NULL)
                } else {
                    stopifnot(identical(dimnames(run$factors[[assay]]$W), dimnames(W[[assay]])))
                    notes[[sprintf("S4/%s/%s/%s/W_package_vs_simulator", cell_name, arm, assay)]] <- sprintf("max abs %.3g",
                        max(abs(run$factors[[assay]]$W - W[[assay]])))
                    run$factors[[assay]]$W <- W[[assay]]
                }
            }
            da <- silently(if (omic == "lipid") runLipidsDA(object, s4$contrastsTable(view$contrasts), "~ 0 + group", 0.05, 0,
                TRUE, TRUE, ruv_factors = run$factors) else runMetabolitesDA(object, s4$contrastsTable(view$contrasts),
                "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors))
            long <- if (omic == "lipid") da$da_lipids_long else da$da_metabolites_long
            id <- if (omic == "lipid") "lipid_id" else "metabolite_id"
            for (assay in s4$ASSAYS) {
                set.seed(seed * 100L + 7L)
                ref <- s4$testWithW(view, assay, W[[assay]], exact = TRUE)
                got <- long[long$assay == assay, ]
                add(sprintf("S4/%s/%s/%s/%s/k", cell_name, arm, omic, assay), as.integer(k_of[[assay]]), as.integer(run$factors[[assay]]$k))
                add_tables(sprintf("S4/%s/%s/%s/%s", cell_name, arm, omic, assay), keyed(ref$feature, ref$contrast, ref$est, ref$p, ref$q),
                    keyed(got[[id]], got$comparison, got$logFC, got$P.Value, got$fdr_qvalue))
            }
        }
    }
}

# ---- compare under the fixed rule ------------------------------------------------------------------------------------
a <- file.path(out_dir, paste0(stage, "-simulator.rds"))
b <- file.path(out_dir, paste0(stage, "-package.rds"))
compared <- file.path(out_dir, paste0(stage, "-compare.tsv"))
saveRDS(sim_side, a)
saveRDS(pkg_side, b)
system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "packaging/verify/compare-results.R", a, b, compared))
rows <- utils::read.delim(compared, stringsAsFactors = FALSE)
summary <- do.call(rbind, lapply(names(sim_side), function(name) {
    mine <- rows[rows$path == name | startsWith(rows$path, paste0(name, "$")) | startsWith(rows$path, paste0(name, "@")), , drop = FALSE]
    data.frame(entry = name, leaves = nrow(mine), max_rel = if (any(is.finite(mine$max_rel))) max(mine$max_rel, na.rm = TRUE) else NA_real_,
        verdict = if (all(mine$verdict == "agree")) "pass" else "fail", note = notes[[name]], stringsAsFactors = FALSE)
}))
w_notes <- notes[grepl("W_package_vs_simulator$", names(notes))]
if (length(w_notes)) summary <- rbind(summary, data.frame(entry = names(w_notes), leaves = NA_integer_, max_rel = NA_real_,
    verdict = "info", note = unlist(w_notes), stringsAsFactors = FALSE))
utils::write.table(summary, file.path(out_dir, paste0(stage, "-summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("stage:", stage), paste("tree:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste("uncommitted changes under R/:", length(system2("git", c("status", "--porcelain", "R"), stdout = TRUE)) > 0L),
    R.version.string, paste("BLAS:", extSoftVersion()[["BLAS"]]),
    paste("versions:", paste(vapply(c("limma", "sva", "limpa", "fastlimpadpc", "ruv", "qvalue"), \(x)
        paste(x, as.character(utils::packageVersion(x))), ""), collapse = "; ")),
    paste("finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))), file.path(out_dir, paste0(stage, "-meta.txt")))
print(summary, right = FALSE)
