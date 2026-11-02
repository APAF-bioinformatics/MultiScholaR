# KIT-076 Simulation Setup: the fixed inputs of the old-vs-new DA comparison, made once on kit_075 so both builds read
# the same bytes. The screen head states come from the recorders' stores (their head hashes are checked against the
# kit_075 oracles), the e2e lanes from their committed fixtures (log2 of the sample columns, first row per id).
# Usage: KIT076_INSTALLED=1 R_LIBS=<kit_075 library> Rscript docs/qa/workflow-execution/kit/kit-076/da-inputs.R <out rds>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); out <- args[[1L]]
`%||%` <- function(a, b) if (is.null(a)) b else a
api <- function(name) get(name, asNamespace("MultiScholaR"))
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
bindings <- jsonlite::read_json("docs/qa/workflow-execution/kit/bindings.json")
root <- bindings$harness$root
head_state <- function(work, oracle, corpus_key) {
    store_dir <- h$kitFindStore(file.path(root, "oracles", work, "project"))
    head <- h$kitReadHead(store_dir, api)
    stopifnot(identical(head$head_sha256, readRDS(oracle)$head_sha256))
    identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
    manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(dirname(store_dir), identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(api("closeDIAStateAdapter")(manager))
    object <- api("withDIAStateScope")(manager, function() manager$getState())
    corpus <- bindings$corpora[[corpus_key]]
    list(state = head$current_state, head_sha256 = head$head_sha256, object = object,
        design_in = utils::read.delim(corpus$files[[corpus$design_file %||% "design.tsv"]]$path, stringsAsFactors = FALSE, check.names = FALSE))
}
fixture_assay <- function(path, id, annotation, samples) {
    table <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    table <- table[!duplicated(table[[id]]), c(id, annotation, samples)]
    table[samples] <- lapply(table[samples], function(values) log2(as.numeric(values) + 1))
    table
}
fixture_design <- function(path) {
    design <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    names(design)[names(design) == "sample"] <- "Run"
    design
}
lipid_design <- fixture_design("tests/testdata/e2e/lipid_canonical/design.tsv")
lipid_lane <- methods::new("LipidomicsAssayData", lipid_data = list(
        LCMS_Pos = fixture_assay("tests/testdata/e2e/lipid_canonical/lipidsearch_lcms_pos.txt", "LipidName", "LipidClass", lipid_design$Run),
        LCMS_Neg = fixture_assay("tests/testdata/e2e/lipid_canonical/lipidsearch_lcms_neg.txt", "LipidName", "LipidClass", lipid_design$Run)),
    design_matrix = lipid_design, sample_id = "Run", group_id = "group", lipid_id_column = "LipidName", annotation_id_column = "LipidClass")
metab_design <- fixture_design("tests/testdata/e2e/metab_combined/design.tsv")
metab_id <- "Feature.Name"; metab_annotation <- "Platform"
metab_lane <- methods::new("MetaboliteAssayData", metabolite_data = list(
        LCMS_Pos = fixture_assay("tests/testdata/e2e/metab_combined/combined_lcms_features.tsv", metab_id, metab_annotation, metab_design$Run),
        GCMS = fixture_assay("tests/testdata/e2e/metab_combined/combined_gcms_features.tsv", metab_id, metab_annotation, metab_design$Run)),
    design_matrix = metab_design, sample_id = "Run", group_id = "group", metabolite_id_column = metab_id, annotation_id_column = metab_annotation)
inputs <- list(
    lipid = list(screen = head_state("lipid31_lipid_screen", "docs/qa/workflow-execution/kit/kit-075/oracles-lipid75_lipid_screen.rds", "lipid_screen"), lane = lipid_lane),
    metab = list(screen = head_state("metab32_metab_screen", "docs/qa/workflow-execution/kit/kit-075/oracles-metab75_metab_screen.rds", "metab_screen"), lane = metab_lane))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(inputs, out, version = 3)
cat(sprintf("inputs %s sha256 %s; lipid head %s (%s); metab head %s (%s)\n", out, digest::digest(file = out, algo = "sha256"),
    inputs$lipid$screen$state, substr(inputs$lipid$screen$head_sha256, 1, 12), inputs$metab$screen$state, substr(inputs$metab$screen$head_sha256, 1, 12)))
