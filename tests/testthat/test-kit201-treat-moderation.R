# KIT-201 (methods ledger VAL-127, VAL-131): TREAT takes the moderated test's trend and robust settings, so the
# reported TREAT statistics come from the moderation the record states. limma's `treat()` re-estimates the prior
# itself and defaults to trend = FALSE, robust = FALSE (limma 3.68.4 R/treat.R:3, :25-38); with those defaults the
# trended robust prior of `eBayes` is computed and then discarded. A DA made before the change replays as it was made
# (replay rule marker `treat_moderation`, D3).
# Route (a), ticket section 9: T-A1 to T-A7. The legacy-design digests of test-prot-ruv-factors-in-model.R and the
# Legacy figure of test-omic-norm-ruv-design.R are rewritten in their own files (ticket 6.2 (a)(ii)).
library(testthat)

kit201Contrasts <- c("groupB-groupA", "groupC-groupA")

# Ticket section 4: 200 features x 12 samples (3 groups of 4), levels N(22, 2^2), residual variance
# 0.1 x exp(-0.541 z) with z the standardised level (a 4-fold range between the 10th and 90th percentile), 10
# hypervariable features (variance x 8) and 20 other features with an effect of +-1.5 in group B.
kit201Matrix <- function(seed = 201001L) {
    withr::local_seed(seed, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    n <- 200L
    groups <- rep(c("A", "B", "C"), each = 4L)
    samples <- sprintf("S%02d", seq_along(groups))
    level <- stats::rnorm(n, 22, 2)
    z <- (level - mean(level)) / stats::sd(level)
    variance <- 0.1 * exp(-0.541 * z)
    hyper <- sample.int(n, 10L)
    variance[hyper] <- variance[hyper] * 8
    effect <- numeric(n)
    effect[sample(setdiff(seq_len(n), hyper), 20L)] <- rep(c(1.5, -1.5), 10L)
    x <- matrix(level, n, length(groups)) + outer(effect, as.numeric(groups == "B")) +
        matrix(stats::rnorm(n * length(groups)), n, length(groups)) * sqrt(variance)
    dimnames(x) <- list(sprintf("F%03d", seq_len(n)), samples)
    design <- data.frame(Run = samples, group = groups, replicates = seq_along(groups), stringsAsFactors = FALSE)
    rownames(design) <- samples
    list(x = x, design = design)
}

kit201Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

# The contrasted fit as the kernels form it (`~ 0 + group`, limma::contrasts.fit), without moderation.
kit201ContrastFit <- function(fixture, contrasts = kit201Contrasts) {
    design_m <- stats::model.matrix(~ 0 + group, fixture$design)
    fit <- limma::lmFit(fixture$x[, rownames(design_m)], design_m)
    limma::contrasts.fit(fit, limma::makeContrasts(contrasts = contrasts, levels = colnames(design_m)))
}

kit201Prot <- function(fixture, lfc, trend = TRUE, robust = TRUE, contrasts = kit201Contrasts) {
    kit201Quiet(runTestsContrasts(fixture$x, contrast_strings = contrasts, design_matrix = fixture$design,
        formula_string = "~ 0 + group", treat_lfc_cutoff = lfc, eBayes_trend = trend, eBayes_robust = robust))
}

kit201Kernels <- list(lipid = list(fit = "fitLipidDaAssayModel", apply = "applyLipidDaContrasts"),
    metab = list(fit = "fitMetabDaAssayModel", apply = "applyMetabDaContrasts"))

kit201Omics <- function(kernel, fixture, lfc, trend = TRUE, robust = TRUE) {
    fit <- kit201Quiet(get(kernel$fit)(fixture$x, fixture$design, "~ 0 + group", sample_id_col = "Run"))
    kit201Quiet(get(kernel$apply)(fit, kit201Contrasts, lfc, trend, robust))
}

# A table's p-values in the fit's feature order.
kit201P <- function(table, fit, column) table[[column]][match(rownames(fit$coefficients), rownames(table))]

test_that("T-A1: the lipid and metabolomics TREAT step takes the trended robust prior", {
    fixture <- kit201Matrix()
    cfit <- kit201ContrastFit(fixture)
    expected <- limma::treat(limma::eBayes(cfit, trend = TRUE, robust = TRUE), lfc = 0.5, trend = TRUE, robust = TRUE)
    for (name in names(kit201Kernels)) {
        result <- kit201Omics(kit201Kernels[[name]], fixture, 0.5)
        # Fails on base at 1: treat's defaults give one prior variance for all features.
        expect_gt(length(unique(result$fit.eb$s2.prior)), 1L, label = paste(name, "distinct prior variances"))
        expect_identical(result$fit.eb$s2.prior, expected$s2.prior, label = paste(name, "s2.prior"))
        expect_length(result$fit.eb$df.prior, nrow(fixture$x))
    }
})

test_that("T-A2: proteomics' legacy test at a cutoff of 0 gives the moderated test's p-values", {
    fixture <- kit201Matrix()
    result <- kit201Prot(fixture, 0)
    # Fails on base at 1 (asserted first).
    expect_gt(length(unique(result$fit.eb$s2.prior)), 1L)
    expected <- limma::eBayes(kit201ContrastFit(fixture), trend = TRUE, robust = TRUE)
    for (contrast in kit201Contrasts) {
        expect_equal(kit201P(result$results[[contrast]], expected, "raw_pvalue"), unname(expected$p.value[, contrast]),
            tolerance = 1e-10, label = contrast)
    }
})

test_that("T-A3: proteomics at a cutoff of 1, and the lipid and metabolomics S4 methods (TREAT at 0), take the trended robust prior", {
    fixture <- kit201Matrix()
    cfit <- kit201ContrastFit(fixture)
    result <- kit201Prot(fixture, 1)
    expected <- limma::treat(limma::eBayes(cfit, trend = TRUE, robust = TRUE), lfc = 1, trend = TRUE, robust = TRUE)
    expect_gt(length(unique(result$fit.eb$s2.prior)), 1L, label = "proteomics at 1: distinct prior variances")
    expect_identical(result$fit.eb$s2.prior, expected$s2.prior)
    expect_length(result$fit.eb$df.prior, nrow(fixture$x))
    # The S4 methods test through runTestsContrasts at their default cutoff of 0 (TREAT runs at any non-NA cutoff).
    table <- data.frame(Name = rownames(fixture$x), fixture$x, check.names = FALSE, stringsAsFactors = FALSE)
    design <- fixture$design[, c("Run", "group")]
    rownames(design) <- NULL
    lipid <- createLipidomicsAssayData(lipid_data = list(Assay1 = table), design_matrix = design, lipid_id_column = "Name",
        sample_id = "Run", group_id = "group")
    metab <- methods::new("MetaboliteAssayData", metabolite_data = list(Assay1 = table), metabolite_id_column = "Name",
        annotation_id_column = "Name", database_identifier_type = "Name", internal_standard_regex = "", design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = NA_character_, args = list())
    at0 <- limma::treat(limma::eBayes(kit201ContrastFit(fixture, "groupB-groupA"), trend = TRUE, robust = TRUE), lfc = 0,
        trend = TRUE, robust = TRUE)
    for (object in list(lipid = lipid, metab = metab)) {
        s4 <- kit201Quiet(differentialAbundanceAnalysis(list(Assay1 = object),
            contrasts_tbl = data.frame(contrasts = "groupB-groupA", stringsAsFactors = FALSE), formula_string = "~ 0 + group",
            treat_lfc_cutoff = 0))
        fit <- s4$Assay1@fit.eb
        expect_gt(length(unique(fit$s2.prior)), 1L, label = paste(class(object), "distinct prior variances"))
        expect_identical(fit$s2.prior, at0$s2.prior, label = paste(class(object), "s2.prior"))
    }
})

test_that("T-A4: guards: trend and robust off, lipid/metab at 0, the declared designs and limpa are unchanged", {
    fixture <- kit201Matrix()
    cfit <- kit201ContrastFit(fixture)
    earlier <- list(MultiScholaR.treat_moderation_rule = "constant")
    # (a) both settings off: today's TREAT statistics exactly, under either rule.
    constant <- limma::treat(limma::eBayes(cfit, trend = FALSE, robust = FALSE), lfc = 0.5)
    for (name in names(kit201Kernels)) {
        off <- kit201Omics(kit201Kernels[[name]], fixture, 0.5, FALSE, FALSE)
        for (contrast in kit201Contrasts) {
            expect_identical(kit201P(off$results[[contrast]], cfit, "P.Value"),
                kit201P(limma::topTreat(constant, coef = contrast, n = Inf), cfit, "P.Value"), label = paste(name, contrast))
        }
        expect_identical(withr::with_options(earlier, kit201Omics(kit201Kernels[[name]], fixture, 0.5, FALSE, FALSE)), off,
            label = paste(name, "off, either rule"))
    }
    off <- kit201Prot(fixture, 0, FALSE, FALSE)
    constant0 <- limma::treat(limma::eBayes(cfit, trend = FALSE, robust = FALSE), lfc = 0)
    for (contrast in kit201Contrasts) {
        expect_identical(kit201P(off$results[[contrast]], cfit, "raw_pvalue"),
            kit201P(limma::topTreat(constant0, coef = contrast, n = Inf), cfit, "P.Value"), label = paste("proteomics", contrast))
    }
    expect_identical(withr::with_options(earlier, kit201Prot(fixture, 0, FALSE, FALSE)), off)
    # (b) lipid and metabolomics at 0 run no TREAT: the moderated test, under either rule.
    moderated <- limma::eBayes(cfit, trend = TRUE, robust = TRUE)
    for (name in names(kit201Kernels)) {
        at0 <- kit201Omics(kit201Kernels[[name]], fixture, 0)
        for (contrast in kit201Contrasts) {
            expect_identical(kit201P(at0$results[[contrast]], cfit, "P.Value"), unname(moderated$p.value[, contrast]),
                label = paste(name, contrast))
        }
        expect_identical(withr::with_options(earlier, kit201Omics(kit201Kernels[[name]], fixture, 0)), at0)
    }
    # (c) the declared proteomics designs already pass the settings (M30): unchanged, and the rule is not read there.
    declared <- function() kit201Quiet(runProtRuvModelContrasts(fixture$x, kit201Contrasts, fixture$design, "~ 0 + group",
        treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE, exact = TRUE))
    designs <- declared()
    expect_gt(length(unique(designs$fit.eb$s2.prior)), 1L)
    expect_identical(withr::with_options(earlier, declared()), designs)
    # (d) DIA with limpa has no TREAT step: its contrast half is unchanged under either rule.
    design_m <- stats::model.matrix(~ 0 + group, fixture$design)
    limpa_fit <- limma::lmFit(fixture$x[, rownames(design_m)], design_m)
    limpa <- function() kit201Quiet(convertDpcDAToStandardFormat(limpa_fit, "BvsA=groupB-groupA", design_m))
    expect_identical(withr::with_options(earlier, limpa()), limpa())
})

# ---- T-A5: replay (D3) ----------------------------------------------------------------------------------------------

# A long table in a canonical order: columns by name, rows by every character column, no row names (VAL-79: a replayed
# table is compared by column name with its rows aligned, as KIT-177 T5 does).
kit201Canon <- function(table) {
    table <- as.data.frame(table)
    table <- table[, sort(names(table)), drop = FALSE]
    keys <- names(table)[vapply(table, is.character, logical(1))]
    if (length(keys)) table <- table[do.call(order, unname(as.list(table[keys]))), , drop = FALSE]
    rownames(table) <- NULL
    table
}

# The DA's long tables of the head's accepted manifest, by contrast (and assay for lipidomics), and its model record.
kit201DaTables <- function(manager, long_table) withDIAStateScope(manager, function() {
    context <- manager$.context
    manifest <- resolveWorkflowStateValue(context, context$snapshot$da_results)
    out <- list()
    for (contrast in names(manifest$results)) {
        entry <- manifest$results[[contrast]]
        if (!is.null(entry$assays)) {
            for (assay in names(entry$assays))
                out[[paste(contrast, assay)]] <- kit201Canon(materializeWorkflowTable(context, entry$assays[[assay]]$tables[[long_table]]))
        } else out[[contrast]] <- kit201Canon(readDIAContrastOutput(entry, long_table)[[long_table]])
    }
    list(model = manifest$model, tables = out)
})

kit201PValues <- function(tables) lapply(tables, function(table) table[[intersect(c("raw_pvalue", "P.Value"), names(table))[[1L]]]])

# A session continued in a later build: every saved identity carries the build (`context$software`, the installed
# binary's digest in the application; R1-04), so the re-run computes rather than reuse the earlier build's products.
kit201NextBuild <- function(manager, da) {
    manager$.context$software <- paste(manager$.context$software, "kit201-next-build")
    kit201Quiet(da())
}

# A proteomics LFQ owner: the design state from the fixture, then the DA (legacy route, RUV skipped by construction) at
# a cutoff of 0, made under `rule` (NULL: as shipped; "constant": as before the change).
kit201ProtOwner <- function(root, rule = NULL) {
    fixture <- kit201Matrix()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, basename(root), "LFQ"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("LFQ")
    wide <- data.frame(Protein.Ids = rownames(fixture$x), fixture$x, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(wide) <- NULL # the protein table carries its ids as a column only (R1 instrument defect, run 45328062374)
    design <- fixture$design
    rownames(design) <- NULL
    design$tech_rep_group <- design$Run
    object <- ProteinQuantitativeData(protein_quant_table = wide, protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = wide$Protein.Ids, stringsAsFactors = FALSE), design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group", args = list())
    setDIAStepProvenance(manager, "design_builder", fixed = list(pivot = "protein x run", transform = "log2(x); 0 missing"))
    manager$saveState("protein_s4_initial", object, list(), "design")
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    contrasts <- data.frame(contrasts = "B.vs.A", full_format = "B.vs.A=groupB-groupA", friendly_names = "B.vs.A",
        stringsAsFactors = FALSE)
    da <- function() withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        workflow$contrasts_tbl <- contrasts
        runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
            da_data, workflow, list(), contrasts = contrasts, publish = FALSE)
    })
    kit201Quiet(if (is.null(rule)) da() else withr::with_options(list(MultiScholaR.treat_moderation_rule = rule), da()))
    list(manager = manager, rerun = function() kit201NextBuild(manager, da), long_table = "da_proteins_long",
        preconditions = c("initial", "protein_s4_initial"), inputs = character())
}

# A lipidomics owner (the KIT-031 fixture's two assays without their duplicate ids, on the log2 scale): the design state,
# then the DA at a cutoff of 0.5, made under `rule`.
kit201LipidOwner <- function(root, rule = NULL) {
    object <- lipidDeclarationsObject()
    for (assay in names(object@lipid_data)) {
        table <- object@lipid_data[[assay]]
        table <- table[!duplicated(table$LipidName), , drop = FALSE]
        table[object@design_matrix$Run] <- log2(table[object@design_matrix$Run])
        object@lipid_data[[assay]] <- table
    }
    fixture <- lipidDeclarationsOwner(root, object = object)
    owner <- fixture$owner
    da <- function() runLipidReferenceDaAnalysis(owner, getWorkflowStateReference(owner),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0.5), fixture$workflow$contrasts_tbl)
    kit201Quiet(if (is.null(rule)) da() else withr::with_options(list(MultiScholaR.treat_moderation_rule = rule), da()))
    store_root <- dirname(dirname(list.files(fixture$source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE,
        full.names = TRUE)[[1L]]))
    list(manager = owner, rerun = function() kit201NextBuild(owner, da), long_table = "da_lipids_long", paths = fixture$paths,
        preconditions = c("initial", "lipid_raw_data_s4"), source_dir = fixture$source_dir, store_root = store_root,
        inputs = c(lipid_assay_1 = file.path(fixture$source_dir, "import_sources", "LCMS_Neg_export.txt"),
            lipid_assay_2 = file.path(fixture$source_dir, "import_sources", "LCMS_Pos_export.txt"),
            design_file = file.path(fixture$source_dir, "design_matrix.tab")))
}

# Write the owner's record, replay it into a copy holding only the preconditions (as the harness does) and return the
# record, the original DA tables and the replayed ones.
kit201Replay <- function(session, root) {
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(session$manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    origin <- kit201DaTables(session$manager, session$long_table)
    source_root <- session$source_dir %||% root
    store_root <- session$store_root %||% root
    copy <- withr::local_tempdir(.local_envir = parent.frame())
    stopifnot(system2("cp", c("-a", file.path(source_root, "."), copy)) == 0L)
    copy_root <- if (identical(store_root, root)) copy else file.path(copy, substring(store_root, nchar(source_root) + 2L))
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- h$kitPrepareReplayProject(copy_root, session$preconditions, api)
    replayed <- NULL
    kit201Quiet(replayWorkflowAnalysisRecord(record, prepared, session$inputs, session$paths %||% list(),
        withr::local_tempdir(.local_envir = parent.frame()), on_step = function(step, replay_manager) {
            if (identical(step$operation, record$analyses$differential_abundance$operation))
                replayed <<- kit201DaTables(replay_manager, session$long_table)
        }))
    list(record = record, origin = origin, replayed = replayed)
}

test_that("T-A5: a DA recorded before the change replays with the earlier TREAT moderation, one made after with the current; a re-run computes the current", {
    skip_if_not(workflowOwnerCapable())
    # Fails on base at: the marker is not defined (NULL).
    expect_false(is.null(workflowReplayRuleMarkers$treat_moderation))
    for (omic in c("proteomics", "lipidomics")) {
        make <- if (identical(omic, "proteomics")) kit201ProtOwner else kit201LipidOwner
        # (a) made before the change: no marker in the DA model; the replay reproduces it exactly.
        root_a <- withr::local_tempdir()
        a <- make(root_a, rule = "constant")
        withr::defer(closeDIAStateAdapter(a$manager))
        replay_a <- kit201Replay(a, root_a)
        expect_null(replay_a$record$analyses$differential_abundance$model$treat_moderation, label = paste(omic, "(a) record"))
        expect_null(replay_a$replayed$model$treat_moderation, label = paste(omic, "(a) replayed model"))
        expect_identical(replay_a$replayed$tables, replay_a$origin$tables, label = paste(omic, "(a) replayed tables"))
        # (b) made as shipped: the marker, and the replay reproduces it exactly.
        root_b <- withr::local_tempdir()
        b <- make(root_b)
        withr::defer(closeDIAStateAdapter(b$manager))
        replay_b <- kit201Replay(b, root_b)
        expect_identical(unlist(replay_b$record$analyses$differential_abundance$model$treat_moderation), "as_main_test",
            label = paste(omic, "(b) record"))
        expect_identical(replay_b$replayed$model$treat_moderation, "as_main_test", label = paste(omic, "(b) replayed model"))
        expect_identical(replay_b$replayed$tables, replay_b$origin$tables, label = paste(omic, "(b) replayed tables"))
        # The two rules give different statistics (else the replay could not tell them apart).
        expect_false(isTRUE(all.equal(kit201PValues(replay_a$origin$tables), kit201PValues(replay_b$origin$tables))),
            label = paste(omic, "(a) and (b) differ"))
        # (c) the session made before the change, re-run now: the current test, recorded with the marker.
        a$rerun()
        rerun <- kit201DaTables(a$manager, a$long_table)
        expect_identical(rerun$model$treat_moderation, "as_main_test", label = paste(omic, "(c) model"))
        expect_equal(kit201PValues(rerun$tables), kit201PValues(replay_b$origin$tables), tolerance = 1e-12,
            label = paste(omic, "(c) p-values"))
    }
})

test_that("T-A6: the legacy result identity names the TREAT moderation only under the current rule", {
    base_string <- "full-population-qvalue-BH-v1:raw_pvalue:fdr_qvalue:fdr_value_bh_adjustment"
    fixture <- kit201Matrix()
    captured <- NULL
    fake <- list(.context = list(da_reference_execution = TRUE, da_base_fit_reference = NULL, software = "kit201", store = NULL))
    local_mocked_bindings(getDIAActiveManager = function(...) fake, usesDIAReferenceStates = function(...) TRUE,
        fitDIAExecutionBaseModel = function(...) list(schema = "prot_da_base_fit_v1", fit = NULL),
        newWorkflowDiaResultIdentity = function(...) {
            captured <<- list(...)[[6L]]
            rlang::abort("identity captured", class = "kit201_identity_captured")
        }, .package = "MultiScholaR")
    adjustment <- function() {
        captured <<- NULL
        expect_error(kit201Prot(fixture, 0), class = "kit201_identity_captured")
        captured
    }
    expect_identical(withr::with_options(list(MultiScholaR.treat_moderation_rule = "constant"), adjustment()), base_string)
    # Fails on base here: the strings are equal.
    expect_identical(adjustment(), paste0(base_string, ":treat-moderated-v2"))
})

test_that("T-A7: the in-app Legacy row is restated from KIT-124 S1 arm 28, with the reworded sentence (D2)", {
    legacy <- stats::setNames(omicGuiRuvTableClassic()$options, vapply(omicGuiRuvTableClassic()$options, `[[`, "", "id"))$legacy
    expect_identical(legacy$fdp, 0.389)
    expect_match(legacy$note, "10.1%", fixed = TRUE)
    expect_match(legacy$note, "For comparison with the approach of earlier versions; analyses saved with earlier versions replay unchanged.",
        fixed = TRUE)
})

# ---- T-196: TREAT with features without values (KIT-196 audit K196-I04) ---------------------------------------------

# KIT-196's fixture (`kit196Object` of test-kit196-da-no-value-features.R at 48bfe085, copied): the shared contrasts
# fixture with, in assay `pos`, the scatter (NA at F1_2 S1, F1_3 S5 and S9, F1_4 S12) and optionally the empty feature
# F1_empty (NA in every sample).
kit201Kit196Object <- function(omic, empty = TRUE) {
    object <- omicsDaContrastsFixture(omic)
    slot_name <- if (identical(omic, "lipid")) "lipid_data" else "metabolite_data"
    id <- if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
    tables <- methods::slot(object, slot_name)
    pos <- tables$pos
    samples <- object@design_matrix$Run
    pos[pos[[id]] == "F1_2", "S1"] <- NA_real_
    pos[pos[[id]] == "F1_3", c("S5", "S9")] <- NA_real_
    pos[pos[[id]] == "F1_4", "S12"] <- NA_real_
    if (empty) {
        row <- pos[1L, , drop = FALSE]
        row[[id]] <- "F1_empty"
        row[[setdiff(names(pos), c(id, samples))[[1L]]]] <- "Name empty"
        row[, samples] <- NA_real_
        pos <- rbind(pos, row)
        rownames(pos) <- NULL
    }
    tables$pos <- pos
    methods::slot(object, slot_name) <- tables
    object
}

kit201Kit196Long <- function(omic, object) {
    result <- kit201Quiet(omicsDaKernel(omic)(theObject = object, contrasts_tbl = omicsDaContrastsTable(),
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0.5, eBayes_trend = TRUE,
        eBayes_robust = TRUE))
    expect_null(result$failed_assays)
    long <- as.data.frame(result[[if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"]])
    id <- if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
    long[order(long$assay, long$comparison, long[[id]]), , drop = FALSE]
}

test_that("T-196 (guard): TREAT above 0 with a feature without values completes, leaves it untested and tests the rest as without it", {
    for (omic in c("lipid", "metab")) {
        id <- if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
        long <- kit201Kit196Long(omic, kit201Kit196Object(omic))
        reference <- kit201Kit196Long(omic, kit201Kit196Object(omic, empty = FALSE))
        empty <- long[long[[id]] %in% "F1_empty", , drop = FALSE]
        expect_identical(sort(as.character(empty$comparison)), sort(omicsDaContrastsTable()$contrasts), label = paste(omic, "F1_empty once per contrast"))
        statistics <- intersect(c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "fdr_qvalue", "fdr_value_bh", "raw_pvalue"), names(long))
        for (column in statistics) expect_true(all(is.na(empty[[column]])), label = paste(omic, column, "NA for F1_empty"))
        tested <- long[!long[[id]] %in% "F1_empty", , drop = FALSE]
        rownames(tested) <- rownames(reference) <- NULL
        expect_identical(tested, reference, label = paste(omic, "tested rows"))
    }
})
