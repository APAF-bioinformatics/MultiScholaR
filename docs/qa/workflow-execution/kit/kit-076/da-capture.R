# KIT-076: run the omics DA on one install over the frozen scenarios and keep everything the old-vs-new comparison
# reads: the public kernel's list, the wrapper's list, the per-contrast fits of a v2 manifest, and the exported files.
# Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-076/da-capture.R <inputs rds> <out dir>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); inputs <- readRDS(args[[1L]]); out <- args[[2L]]
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
unlink(out, recursive = TRUE); dir.create(out, recursive = TRUE)
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
contrast_table <- function(contrasts, friendly) {
    data.frame(contrasts = contrasts, friendly_names = friendly, full_format = paste0(friendly, "=", contrasts), stringsAsFactors = FALSE)
}
with_batch <- function(object, design_in) {
    design <- object@design_matrix
    design$batch <- design_in$batch[match(design[[object@sample_id]], design_in$sample)]
    stopifnot(!anyNA(design$batch))
    object@design_matrix <- design
    object
}
with_group_batch <- function(object, design_in) {
    object <- with_batch(object, design_in)
    design <- object@design_matrix
    design[[object@group_id]] <- paste(design[[object@group_id]], design$batch, sep = "_")
    object@design_matrix <- design
    object
}
scenarios <- function(omic) {
    screen <- inputs[[omic]]$screen
    groups <- paste0("group", c("KO", "WT"))
    pair <- function(a, b) paste0("group", a, "-group", b)
    list(
        screen = list(object = screen$object, formula = "~ 0 + group", lfc = 0, rule = "identical",
            contrasts = contrast_table(pair("KO", "WT"), "KO_vs_WT")),
        multi = list(object = with_group_batch(screen$object, screen$design_in), formula = "~ 0 + group", lfc = 0, rule = "identical",
            contrasts = contrast_table(c(pair("KO_B1", "WT_B1"), pair("KO_B2", "WT_B2"), pair("KO_B1", "KO_B2")), c("B1_KO_vs_WT", "B2_KO_vs_WT", "KO_B1_vs_B2"))),
        multi_treat = list(object = with_group_batch(screen$object, screen$design_in), formula = "~ 0 + group", lfc = 0.5, rule = "identical",
            contrasts = contrast_table(c(pair("KO_B1", "WT_B1"), pair("KO_B2", "WT_B2"), pair("KO_B1", "KO_B2")), c("B1_KO_vs_WT", "B2_KO_vs_WT", "KO_B1_vs_B2"))),
        nonorth = list(object = with_batch(screen$object, screen$design_in), formula = "~ 0 + group + batch", lfc = 0, rule = "tolerance",
            contrasts = contrast_table(c(pair("KO", "WT"), "batchB2-batchB3"), c("KO_vs_WT", "B2_vs_B3"))),
        e2e = list(object = inputs[[omic]]$lane, formula = "~ 0 + group", lfc = 0, rule = "identical",
            contrasts = contrast_table(pair("KO", "WT"), "KO_vs_WT")))
}
session_owner <- function(object, omic, root) {
    source_dir <- file.path(root, "source"); dir.create(source_dir, recursive = TRUE)
    session <- new.env(parent = emptyenv()); session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    workflow <- new.env(parent = emptyenv()); workflow$config_list <- list()
    omic_type <- if (identical(omic, "lipid")) "lipidomics" else "metabolomics"
    workflow$state_manager <- api("createDIASessionOwner")(workflow, list(source_dir = source_dir), session, omic_type = omic_type)
    shiny::withReactiveDomain(session, api("bindWorkflowOwnerTypeAtDesign")(workflow, toupper(omic_type)))
    workflow$state_manager$saveState("da_input", object, list(), "KIT-076 comparison input")
    workflow$state_manager
}
for (omic in c("lipid", "metab")) {
    kernel <- if (identical(omic, "lipid")) api("runLipidsDA") else api("runMetabolitesDA")
    wrapper <- api(if (identical(omic, "lipid")) "runLipidReferenceDaAnalysis" else "runMetabReferenceDaAnalysis")
    writer <- api(if (identical(omic, "lipid")) "outputLipidDaResultsAllContrasts" else "outputMetabDaResultsAllContrasts")
    all <- scenarios(omic)
    for (name in names(all)) {
        scenario <- all[[name]]
        started <- Sys.time()
        root <- file.path(out, omic, name); dir.create(root, recursive = TRUE)
        direct <- suppressWarnings(suppressMessages(kernel(theObject = scenario$object, contrasts_tbl = scenario$contrasts,
            formula_string = scenario$formula, da_q_val_thresh = 0.05, treat_lfc_cutoff = scenario$lfc, eBayes_trend = TRUE, eBayes_robust = TRUE)))
        manager <- session_owner(scenario$object, omic, root)
        source <- api("getWorkflowStateReference")(manager)
        run <- suppressWarnings(suppressMessages(wrapper(manager, source,
            list(formula_string = scenario$formula, da_q_val_thresh = 0.05, treat_lfc_cutoff = scenario$lfc), scenario$contrasts)))
        fits <- list()
        if (grepl("\\.v2$", run$manifest$schema)) for (contrast in names(run$manifest$results)) {
            product <- run$manifest$results[[contrast]]$products$da_contrast
            value <- api("withWorkflowStateScope")(manager, function() api("hydrateSavedWorkflowEntity")(manager$.context$store, product)$value)
            fits[[contrast]] <- lapply(value$assays, `[[`, "fit.eb")
        }
        results <- run$results
        results$theObject <- scenario$object
        files <- file.path(root, c("da", "pub")); for (d in files) dir.create(d)
        written <- suppressWarnings(suppressMessages(writer(da_results_list = results, da_output_dir = files[[1L]], publication_graphs_dir = files[[2L]],
            da_q_val_thresh = 0.05, lfc_threshold = 1)))
        results$theObject <- NULL
        direct$theObject <- NULL
        saveRDS(list(omic = omic, scenario = name, formula = scenario$formula, rule = scenario$rule, contrasts = scenario$contrasts, schema = run$manifest$schema,
            direct = direct, results = results, fits = fits, written = written), file.path(root, "capture.rds"), version = 3)
        api("closeDIAStateAdapter")(manager)
        cat(sprintf("%s/%s: schema %s, %d long rows, written %s, %.0f s\n", omic, name, run$manifest$schema, nrow(results[[run$manifest$long_table %||% if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"]]),
            written, as.numeric(difftime(Sys.time(), started, units = "secs"))))
    }
}
