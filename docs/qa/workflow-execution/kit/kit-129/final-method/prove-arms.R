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
#   count KIT-129 item 3: each omic's new-design kernel computes k itself by the estimated count, and that k equals the
#         simulator's arm 32 E (arm 31 in transcriptomics; B = 200, the biological units removed on the
#         technical-duplicate cell): classic proteomics `computeProtRuvModelFactors(ruv_k_rule = "num_sv_be")` (S1),
#         DIA with limpa the declared kernel (S2), transcriptomics `trxChooseRuvDesign` with the count (S3, run 1 of C0
#         and C6), lipidomics and metabolomics `runOmicsRuvFactors(ruv_k_rule = "num_sv_be")` (S4, on the view of the
#         features complete in study and pooled-QC injections, the features S4's E is computed on)
#   trx   KIT-129 integration item 4 (decision M29 in transcriptomics): the package's test of a comparison with W in the
#         model (`trxDaApply`: `edgeR::glmQLFTest(fit, contrast = )`) against the same test with the comparison fitted
#         as the first coefficient of the reparameterised model (`limma::contrastAsCoef`, arm 25's construction in S1,
#         S2, S4), on S3's data (run 1 of C0 and C6) with arm 31's W (RUVs, all genes, E with B = 200). S3 ran no arm
#         25 (the exactness round's declaration: nothing is approximate in this route), so the reference is the
#         reparameterised fit, not an arm's table
#   runpath KIT-129 integration item 5: the run path (the compute halves the tab and the managed worker call, in a
#         reference store) against arm 32: classic proteomics `computeDIAReferenceNormalization` with the new design
#         (the estimated count) then `computeDIAReferenceDaContrast` (S1); DIA with limpa `computeDIAReferenceRuvCovariates`
#         with all proteins as controls (its own count) then `computeDIAReferenceDaContrast` on the limpa route (S2)
#   techrep KIT-129 integration (second pass) item 1: the technical-replicate design at the kernels (the biological
#         units as RUV-III's replicate sets: `ruv_grouping_variable` = the unit column), run 1 of C6 and C8, S1 and S4:
#         the kernel's own count against arm 29's model with M30's 200 permutations (= arm 33's E); its W against arm
#         16's W (S1); its DA against arm 33 (arm 16 with arm 25's standard error, arm 26's test and that count; at k = 1
#         where the count is 0, so that the design is exercised); its corrected matrix, tested as arm 15 tests, against
#         arm 15 (arm 15's own k, 1 where E is 0). S4 on the view of the features complete in study and pooled-QC
#         injections on both sides (the simulator's `complete` subset), so that the package computes its own W
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
stopifnot(stage %in% c("m29", "m30a", "m30b", "count", "trx", "runpath", "techrep"))
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

# The run path's store (item 5): an in-process state owner holding `object` as the step's source; the compute halves of
# the tab and the managed worker run on it (`computeDIAReference*`), exactly as in the application.
runpath_store <- function(object) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(tempfile("runpath-"), "kit-129-runpath", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_replicate_filtered", object, list(), "source")
    manager
}
# One contrast's DA through the compute half on `state`, its contrasts table keyed.
runpath_da <- function(manager, state, contrast, route, id_column) {
    spec <- readDIAStateMetadata(manager, state)$data
    step <- list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0, weights = NA, contrast = contrast,
        annotations = spec$components$protein_id_table$entity_id, name = sub("=.*$", "", contrast), comparison = sub("=.*$", "", contrast))
    run <- computeDIAReferenceDaContrast(manager, state, step, route)$run
    tab <- as.data.frame(materializeDIATable(manager$.context, run$value$tables$contrasts_results_table))
    list(table = keyed(tab[[id_column]], "contrast", tab$logFC, tab$raw_pvalue, tab$fdr_qvalue), model = run$value$model)
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

cells <- if (stage == "techrep") c("C6", "C8") else c("C0", "C6")

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
    entry <- function(k) entry_state(state(k))
    entry_state <- function(declared) {
        result <- silently(differentialAbundanceAnalysis(declared, contrasts_tbl = data.frame(contrasts = s1_contrast,
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
    if (stage == "runpath") {
        E32 <- count(200L, p$duplicates)
        manager <- runpath_store(normalised)
        parameters <- list(norm_method = "none", ruv_mode = "automatic", auto_percentage_min = 5, auto_percentage_max = 30,
            separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_percentage = 50,
            ruv_k = NULL, grouping = "group", ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = seed * 100L)
        package <- silently(withDIAStateScope(manager, function() {
            run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), parameters)
            c(list(k = run$value$best_k), runpath_da(manager, run$value$final, s1_contrast, "iq", "Protein.Group"))
        }))
        closeDIAStateAdapter(manager)
        add(sprintf("S1/%s/arm32/runpath_k", cell), E32, as.integer(package$k), "arm 32's E against the run path's own count")
        if (E32 > 0L) add_tables(sprintf("S1/%s/arm32/runpath", cell), arm_table(E32, TRUE, TRUE), package$table)
    }
    if (stage == "count") {
        own <- silently(computeProtRuvModelFactors(normalised, ruv_k_rule = "num_sv_be", ruv_grouping_variable = "group",
            ruv_k_seed = seed * 100L))
        add(sprintf("S1/%s/arm32/kernel_k", cell), count(200L, p$duplicates), own$k,
            "arm 32's E against the classic kernel's own count (computeProtRuvModelFactors, num_sv_be)")
    }
    if (stage == "techrep") {
        # the kernel handed the units as replicate sets (the confirmation round's arm 38 call)
        units_factors <- function(...) silently(computeProtRuvModelFactors(normalised, ..., ruv_grouping_variable = "tech_rep_group"))
        own <- units_factors(ruv_k_rule = "num_sv_be", ruv_k_seed = seed * 100L)
        E29 <- count(200L, TRUE)
        add(sprintf("S1/%s/arm29_B200/kernel_k", cell), E29, own$k, sprintf(paste("arm 29's E (the units removed) with M30's",
            "200 permutations (= arm 33's E) against the kernel's own count; arm 29 itself (B = 20): %d; the kernel's limit",
            "k_max = %d, replicate structure %s"), count(20L, TRUE), own$k_max, own$replicate_structure))
        # arm 16's W (s1_classic.R line 626, `in_model(all_controls, E, "tech_rep_group")`, lines 538-553): RUV-III on the
        # complete proteins, every one a control, the declared units as the replicate sets
        tech_W <- function(k) silently(computeRuvCovariates(s1$make_object(xc, design, "tech_rep_group"), k = k,
            controls = rep(TRUE, nrow(xc)), ruv_grouping_variable = "group"))$W
        kk <- if (E29 > 0L) E29 else 1L
        package <- units_factors(ruv_k = kk, ruv_k_rule = "user")
        add(sprintf("S1/%s/arm16/W_k%d", cell, kk), unname(tech_W(kk)), unname(package$W), "arm 16's W against the kernel's W")
        # arm 33 (lines 147, 633: arm 16 with arm 25's standard error, arm 26's test and arm 29's E at 200 permutations)
        set.seed(seed * 100L + 16L)
        tab <- silently(s1$da_variant(x, design, tech_W(kk), exact = TRUE, applied = TRUE))
        add_tables(sprintf("S1/%s/%s/classic", cell, if (E29 > 0L) "arm33" else "arm33_at_k1"), keyed(tab$id, "contrast", tab$estimate,
            tab$p, tab$q), entry_state(attachProtRuvModelFactors(normalised, normalised, package)))
        # arm 15 (line 625, `corrected(all_controls, E, "tech_rep_group")`, lines 554-559) against the corrected matrix the
        # new design forms (`diaNormalizationOperation("declared")`: ruvIII_C_Varying on the object without the
        # declaration, with the kernel's replicate sets, k and controls), each tested as arm 15 tests
        k15 <- if (E20 > 0L) E20 else 1L
        sim15 <- s1$as_matrix(silently(ruvIII_C_Varying(s1$make_object(x, design, "tech_rep_group"), ruv_grouping_variable = "group",
            ruv_number_k = k15, ctrl = stats::setNames(rep(TRUE, nrow(x)), rownames(x)))))
        pkg15 <- s1$as_matrix(silently(ruvIII_C_Varying(protRuvGroupReplicateObject(normalised), ruv_grouping_variable = "tech_rep_group",
            ruv_number_k = k15, ctrl = units_factors(ruv_k = k15, ruv_k_rule = "user")$controls)))
        label15 <- if (E20 > 0L) "arm15" else "arm15_at_k1"
        add(sprintf("S1/%s/%s/corrected", cell, label15), sim15, pkg15[rownames(sim15), colnames(sim15), drop = FALSE],
            "arm 15's corrected matrix against the new design's")
        t15 <- silently(s1$da(sim15, design))
        p15 <- silently(s1$da(pkg15, design))
        add_tables(sprintf("S1/%s/%s/classic", cell, label15), keyed(t15$id, "contrast", t15$estimate, t15$p, t15$q),
            keyed(p15$id, "contrast", p15$estimate, p15$p, p15$q))
    }
}

# ---- S2, DIA with limpa ------------------------------------------------------------------------------------------------
if (stage %in% c("m29", "m30a", "m30b", "count", "runpath")) {
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
        if (stage == "count") {
            package <- kernel()
            add(sprintf("S2/%s/arm%s/kernel_k", cell_name, if (units) "32" else "31"), count(200L, units), package$value$best_k,
                "arm 31's E (C0) or arm 32's E (C6, the units removed) against the declared kernel's own k")
        }
        if (stage == "runpath") {
            E32 <- count(200L, units)
            manager <- runpath_store(quantified)
            manager$setWorkflowType("DIA_LIMPA")
            parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_seed = 100L * seed + 99L)
            package <- silently(withDIAStateScope(manager, function() {
                run <- computeDIAReferenceRuvCovariates(manager, getDIAStateReference(manager), parameters)
                c(list(k = run$value$best_k), runpath_da(manager, run$value$final, "B_vs_A=groupB-groupA", "limpa", "Protein.Ids"))
            }))
            closeDIAStateAdapter(manager)
            add(sprintf("S2/%s/arm32/runpath_k", cell_name), E32, as.integer(package$k), "arm 32's E against the run path's own count")
            if (E32 > 0L) add_tables(sprintf("S2/%s/arm32/runpath", cell_name), arm_table(E32), package$table)
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

# ---- S4, lipidomics and metabolomics (second-pass stages: techrep) ---------------------------------------------------
if (stage %in% c("techrep")) {
    s4 <- simulator("s4_lipid_metab.R", "2d03b103", "^# ---- driver")
    table_s4 <- s4$cellTable()
    for (cell_name in cells) {
        cell <- table_s4[[cell_name]]
        seed <- 100000L + 10000L * cell$index + 1L
        data <- s4$generateData(cell, seed)
        # E and its variants as `runOne` forms them (s4_lipid_metab.R lines 836-852), on the features complete in study and
        # pooled-QC injections
        estimate <- function(B, units) lapply(data$assays, function(a) {
            dat <- a$y[s4$completeFeatures(a), , drop = FALSE]
            if (nrow(dat) < 20L) return(NA_integer_)
            mod <- stats::model.matrix(~ factor(if (units) data$design$extract else data$design$group))
            as.integer(silently(sva::num.sv(dat, mod, method = "be", B = B, seed = seed * 100L + 99L)))
        })
        E20 <- estimate(20L, FALSE)
        variants <- list(E_200 = estimate(200L, FALSE))
        if (cell$duplicates) variants$E_units_200 <- estimate(200L, TRUE)
        # both sides on the view of those features (the simulator's `complete` subset): the package holds the study
        # injections only, and on this view computes the simulator's W itself
        view <- s4$makeView(s4$restrictData(data, lapply(data$assays, s4$completeFeatures)), E20, variants)
        mats <- lapply(view$assays, `[[`, "y")
        objects <- list(lipid = s4$lipidObject(mats, view$design), metab = s4$metabObject(mats, view$design))
        contrasts <- s4$contrastsTable(view$contrasts)
        at_least_one <- function(k) lapply(k, function(v) if (is.na(v) || v == 0L) 1L else as.integer(v))
        # an arm's tables (`modelArm`/`subtractArm`) against the package's DA with the kernel's factors
        compare_arm <- function(label, omic, reference, factors) {
            da <- silently(if (omic == "lipid") runLipidsDA(objects[[omic]], contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE,
                ruv_factors = factors) else runMetabolitesDA(objects[[omic]], contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE,
                ruv_factors = factors))
            long <- if (omic == "lipid") da$da_lipids_long else da$da_metabolites_long
            id <- if (omic == "lipid") "lipid_id" else "metabolite_id"
            for (assay in s4$ASSAYS) {
                ref <- reference$tables[[assay]]
                name <- sprintf("S4/%s/%s/%s/%s", cell_name, label, omic, assay)
                if (s4$isRefusal(ref)) {
                    notes[[paste0(name, "/W_package_vs_simulator")]] <<- paste("the simulator refused:", ref$reason)
                    next
                }
                got <- long[long$assay == assay, ]
                add_tables(name, keyed(ref$feature, ref$contrast, ref$est, ref$p, ref$q),
                    keyed(got[[id]], got$comparison, got$logFC, got$P.Value, got$fdr_qvalue))
            }
        }
        if (stage == "techrep") {
            for (omic in c("lipid", "metab")) {
                units_run <- function(...) silently(runOmicsRuvFactors(objects[[omic]], "groups_all_features", ...,
                    ruv_grouping_variable = "extract"))
                own <- units_run(ruv_k_rule = "num_sv_be", ruv_k_seed = seed * 100L + 99L)
                arm29 <- estimate(20L, TRUE)
                for (assay in s4$ASSAYS) {
                    add(sprintf("S4/%s/arm29_B200/%s/%s/kernel_k", cell_name, omic, assay), variants$E_units_200[[assay]],
                        as.integer(own$factors[[assay]]$k), sprintf(paste("arm 29's E (the extracts removed) with M30's 200",
                            "permutations (= arm 33's E) against runOmicsRuvFactors' own count; arm 29 itself (B = 20): %d; status %s"),
                            arm29[[assay]], own$factors[[assay]]$status))
                }
                # arm 33 (lines 129, 706: arm 16, the extracts as replicate sets, with arm 25's standard error and that count)
                kk <- at_least_one(variants$E_units_200)
                package <- units_run(kk, "user")
                set.seed(seed * 100L + 16L)
                compare_arm(if (all(unlist(variants$E_units_200) > 0L)) "arm33" else "arm33_at_k1_where_0", omic,
                    s4$modelArm(view, s4$allComplete, function(a, assay) kk[[assay]], m_of = s4$extractM, exact = TRUE),
                    package$factors)
                # arm 15 (line 677: the extracts as replicate sets, RUV-III-C subtracted, then the plain test) against the
                # kernel's corrected matrix tested the same way (its feature levels restored: a constant per feature)
                k15 <- at_least_one(E20)
                set.seed(seed * 100L + 15L)
                reference <- s4$subtractArm(view, function(a, assay) k15[[assay]], s4$extractM(view, NULL))
                corrected <- units_run(k15, "user")$corrected
                assays <- if (omic == "lipid") corrected@lipid_data else corrected@metabolite_data
                for (assay in s4$ASSAYS) {
                    name <- sprintf("S4/%s/%s/%s/%s", cell_name, if (all(unlist(E20) > 0L)) "arm15" else "arm15_at_k1_where_0", omic, assay)
                    ref <- reference$tables[[assay]]
                    if (s4$isRefusal(ref)) {
                        notes[[paste0(name, "/W_package_vs_simulator")]] <- paste("the simulator refused:", ref$reason)
                        next
                    }
                    tested <- silently(s4$daKernel(s4$assayMatrix(assays[[assay]], view$design$Sample_ID), view$design, "~ 0 + group",
                        view$contrasts))
                    add_tables(name, keyed(ref$feature, ref$contrast, ref$est, ref$p, ref$q),
                        keyed(tested$feature, tested$contrast, tested$est, tested$p, tested$q))
                }
            }
        }
    }
}

# ---- S3, transcriptomics (stage count) -------------------------------------------------------------------------------
if (stage == "count") {
    s3 <- simulator("s3_trx.R", "2ba5cf87", "^selected <- if")
    for (cell_name in cells) {
        cl <- s3$CELLS[[cell_name]]
        seed <- 100000L + 10000L * cl$idx + 1L
        # the shared preparation of every arm as `oneRun` forms it (s3_trx.R lines 320-330); arm 31's E (line 450)
        prep <- silently({
            d <- s3$generate(cl, seed)
            y <- trxExpressionFilter(d$counts, d$design$group)$filtered_counts
            n <- trxTmmLogCpm(y, d$design$group)
            list(d = d, y = y, n = n)
        })
        E200 <- as.integer(silently(sva::num.sv(prep$n$log_cpm_tmm, stats::model.matrix(~ group, data = prep$d$design), method = "be",
            B = 200, seed = seed)))
        choice <- silently(trxChooseRuvDesign(list(counts = prep$y, design = prep$d$design),
            trxRuvDesignParameters("groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = seed), log_cpm = prep$n$log_cpm_tmm))
        add(sprintf("S3/%s/arm31/kernel_k", cell_name), E200, choice$ruv_design$k,
            "arm 31's E (B = 200; C6's libraries are not summed on this simulator tree) against trxChooseRuvDesign's own count")
    }
}

# ---- S3, transcriptomics (stage trx: item 4) -------------------------------------------------------------------------
if (stage == "trx") {
    s3 <- simulator("s3_trx.R", "2ba5cf87", "^selected <- if")
    # edgeR's own account: the QL F-test of a contrast reformulates the design so that the contrast is a coefficient
    # (`glmLRT`, which `glmQLFTest` calls for the deviance, uses `contrastAsCoef`), then refits without it per gene.
    notes[["S3/edgeR_contrast_reformulation"]] <- sprintf("glmLRT calls contrastAsCoef: %s; glmQLFTest calls glmLRT: %s",
        any(grepl("contrastAsCoef", deparse(edgeR::glmLRT))), any(grepl("glmLRT", deparse(edgeR::glmQLFTest))))
    for (cell_name in cells) {
        cl <- s3$CELLS[[cell_name]]
        seed <- 100000L + 10000L * cl$idx + 1L
        prep <- silently({
            d <- s3$generate(cl, seed)
            y <- trxExpressionFilter(d$counts, d$design$group)$filtered_counts
            n <- trxTmmLogCpm(y, d$design$group)
            list(d = d, y = y, n = n)
        })
        E200 <- as.integer(silently(sva::num.sv(prep$n$log_cpm_tmm, stats::model.matrix(~ group, data = prep$d$design), method = "be",
            B = 200, seed = seed)))
        stopifnot(E200 > 0L)
        # arm 31's W (s3_trx.R line 449: arm 7 with E from 200 permutations): RUVs on upper-quartile counts, all genes
        s <- silently(trxRuvFactors(prep$y, prep$n$log_cpm_tmm, prep$d$design, rownames(prep$y), E200, "RUVs"))
        formula <- trxDaFormula(s$pdata)
        # the package's route: the contrast tested by glmQLFTest(contrast = )
        route <- silently(trxDaApply(trxDaFit(s$ruvs_counts, s$pdata, formula), s3$ENTRY))$results[[1L]]$top_table
        # the reparameterised fit: the comparison as the first coefficient, tested as a coefficient
        model <- trxDaModelMatrix(s$pdata, formula)
        contrast <- trxContrastMatrix(trxParseContrastEntry(s3$ENTRY), model)[, 1L]
        q <- limma::contrastAsCoef(model, contrast, first = TRUE)$design
        exact <- silently({
            dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = s$ruvs_counts, group = s$pdata$group), method = "TMM")
            dge <- edgeR::estimateDisp(dge, q, robust = TRUE)
            edgeR::topTags(edgeR::glmQLFTest(edgeR::glmQLFit(dge, q, robust = TRUE), coef = 1L), n = Inf, sort.by = "PValue")$table
        })
        add(sprintf("S3/%s/arm31/k", cell_name), E200, as.integer(ncol(s$W)), "arm 31's E, the W the route was given")
        add_tables(sprintf("S3/%s/arm31/contrast_vs_reparameterised", cell_name),
            keyed(rownames(exact), "contrast", exact$logFC, exact$PValue, exact$FDR),
            keyed(rownames(route), "contrast", route$logFC, route$PValue, route$FDR))
    }
}

# ---- S4, lipidomics and metabolomics (stage count) -------------------------------------------------------------------
if (stage == "count") {
    s4 <- simulator("s4_lipid_metab.R", "2d03b103", "^# ---- driver")
    table_s4 <- s4$cellTable()
    for (cell_name in cells) {
        cell <- table_s4[[cell_name]]
        seed <- 100000L + 10000L * cell$index + 1L
        data <- s4$generateData(cell, seed)
        # arm 32's E as `runOne` forms it (s4_lipid_metab.R lines 843-848): the features complete in study and QC
        E32 <- lapply(data$assays, function(a) {
            dat <- a$y[s4$completeFeatures(a), , drop = FALSE]
            mod <- stats::model.matrix(~ factor(if (cell$duplicates) data$design$extract else data$design$group))
            as.integer(silently(sva::num.sv(dat, mod, method = "be", B = 200, seed = seed * 100L + 99L)))
        })
        # the package holds the study injections only: it is given the view of those same features
        view <- s4$restrictData(data, lapply(data$assays, s4$completeFeatures))
        mats <- lapply(view$assays, `[[`, "y")
        for (omic in c("lipid", "metab")) {
            object <- if (omic == "lipid") s4$lipidObject(mats, view$design) else s4$metabObject(mats, view$design)
            run <- silently(runOmicsRuvFactors(object, "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = seed * 100L + 99L))
            for (assay in s4$ASSAYS) {
                add(sprintf("S4/%s/arm32/%s/%s/kernel_k", cell_name, omic, assay), E32[[assay]], as.integer(run$factors[[assay]]$k),
                    "arm 32's E against runOmicsRuvFactors' own count (num_sv_be)")
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
w_notes <- notes[grepl("W_package_vs_simulator$|edgeR_contrast_reformulation$", names(notes))]
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
