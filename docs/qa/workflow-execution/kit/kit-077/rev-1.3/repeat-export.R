# KIT-077 rev 1.3 acceptance on a finished journey member, in a new process with no session attached: the reopened owner
# knows its omic, and the DA export from its saved results, with the app's settings, builds no plot and adds no asset. Every
# figure it writes is byte-identical to a stored asset. Counts only; works on a copy of the member's project.
# Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-077/rev-1.3/repeat-export.R <member dir> <lipid|metab>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); member <- args[[1L]]; omic <- args[[2L]]
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
work <- tempfile("kit077-repeat-"); dir.create(work); on.exit(unlink(work, recursive = TRUE), add = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project"), work)) == 0L)
store_dir <- h$kitFindStore(file.path(work, "project"))
identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(dirname(store_dir), identity$project_id, identity$workflow_id), reopen = TRUE)
omic_type <- if (identical(omic, "lipid")) "lipidomics" else "metabolomics"
cat(sprintf("%s: reopened with no session; omic_type %s; figure owner in scope: %s\n", omic, manager$.context$omic_type,
    api("withWorkflowStateScope")(manager, function() !is.null(api("workflowOmicsFigureOwner")(omic_type)))))
stopifnot(identical(manager$.context$omic_type, omic_type))
store <- manager$.context$store
assets <- function() {
    found <- api("findSavedWorkflowEntities")(store, list(kind = "rendered_asset"), max_matches = 2000L)$matches
    vapply(found, function(r) digest::digest(api("hydrateSavedWorkflowEntity")(store, r)$bytes, serialize = FALSE, algo = "sha256"), "")
}
before <- assets()
binding <- api(if (identical(omic, "lipid")) "lipidomicsDeclarationBinding" else "metabolomicsDeclarationBinding")
names_of <- if (identical(omic, "lipid")) c("generateLipidDAVolcanoStatic", "generateLipidDAHeatmap") else c("generateMetabDAVolcanoStatic", "generateMetabDAHeatmap")
builds <- 0L
for (name in names_of) local({
    original <- get(name, ns); force(original)
    unlockBinding(name, ns)
    assign(name, function(...) { builds <<- builds + 1L; original(...) }, envir = ns)
})
pub <- file.path(work, "repeat-pub"); da <- file.path(work, "repeat-da")
export <- api(if (identical(omic, "lipid")) "outputLipidDaResultsAllContrasts" else "outputMetabDaResultsAllContrasts")
reader <- api(if (identical(omic, "lipid")) "readLipidWorkflowState" else "readMetabWorkflowState")
ok <- suppressWarnings(suppressMessages(api("withWorkflowStateScope")(manager, function() {
    context <- manager$.context
    manifest <- api("resolveWorkflowStateValue")(context, context$snapshot$da_results)
    results <- api("readWorkflowOmicsDaResults")(binding, manager, manifest)
    results$theObject <- reader(manager, manifest$source)
    export(da_results_list = results, da_output_dir = da, publication_graphs_dir = pub, da_q_val_thresh = manifest$parameters$qvalue,
        lfc_threshold = manifest$parameters$lfc, heatmap_top_n = 50, heatmap_clustering = "both", heatmap_color_scheme = "RdBu")
})))
files <- list.files(pub, recursive = TRUE, full.names = TRUE, pattern = "\\.(png|pdf)$")
hashes <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "", USE.NAMES = FALSE)
after <- assets()
originals <- list.files(file.path(work, "project"), recursive = TRUE, pattern = "\\.(png|pdf)$")
absent <- setdiff(basename(files), basename(originals))
cat(sprintf("%s: export %s; plots built %d; figures written %d; byte-identical to a stored asset %d; assets %d before, %d after; names not in the journey's project %d\n",
    omic, ok, builds, length(files), sum(hashes %in% before), length(before), length(after), length(absent)))
api("closeDIAStateAdapter")(manager)
pass <- isTRUE(ok) && builds == 0L && length(files) > 0L && all(hashes %in% before) && length(before) == length(after) && !length(absent)
cat(if (pass) "RESULT a repeat export builds no plot and writes stored bytes\n" else "RESULT differs\n")
quit(status = if (pass) 0L else 1L)
