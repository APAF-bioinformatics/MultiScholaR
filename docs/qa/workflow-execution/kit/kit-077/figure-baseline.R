# KIT-077 Simulation Setup: the figure baseline of a screen journey member. The figure files and sizes under the project,
# the QC_figures copy list, and a hash of each QC plot's data per state and assay, from a copy of the member's store.
# Hashes and paths only. Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-077/figure-baseline.R <member dir> <lipid|metab> <out rds>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); member <- args[[1L]]; omic <- args[[2L]]; out <- args[[3L]]
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
work <- tempfile("kit077-baseline-"); dir.create(work); on.exit(unlink(work, recursive = TRUE), add = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project"), work)) == 0L)
project <- file.path(work, "project")
store_dir <- h$kitFindStore(project)
identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(dirname(store_dir), identity$project_id, identity$workflow_id), reopen = TRUE)
files <- list.files(project, pattern = "\\.(png|pdf|html|svg)$", recursive = TRUE)
files <- files[!grepl("^\\.|/\\.", files)]
sizes <- file.size(file.path(project, files))
states <- if (identical(omic, "lipid")) c("lipid_qc_complete", "lipid_normalized", "lipid_ruv_corrected") else c("metab_qc_complete", "metab_normalized", "metab_ruv_corrected")
data_hash <- function(plots) vapply(plots, function(p) if (inherits(p, "ggplot")) digest::digest(as.data.frame(p$data), algo = "sha256") else NA_character_, "")
plot_data <- list()
for (state in intersect(states, unlist(manager$getHistory()))) {
    object <- api("withWorkflowStateScope")(manager, function() manager$getState(state))
    group <- object@group_id
    plot_data[[state]] <- suppressWarnings(suppressMessages(list(
        pca = data_hash(plotPca(object, grouping_variable = group, label_column = "", shape_variable = NULL, title = "", font_size = 8)),
        rle = data_hash(plotRle(object, grouping_variable = group, yaxis_limit = c(-6, 6))),
        density = data_hash(plotDensity(object, grouping_variable = group, title = "", font_size = 8)),
        correlation = data_hash(plotPearson(object, correlation_group = group)))))
}
api("closeDIAStateAdapter")(manager)
saveRDS(list(member = member, omic = omic, files = data.frame(path = files, bytes = sizes, stringsAsFactors = FALSE), plot_data = plot_data), out, version = 3)
cat(sprintf("%s: %d figure files; %d empty; plot data for states %s; plot classes hashed %d of %d\n", omic, length(files), sum(sizes == 0),
    paste(names(plot_data), collapse = ", "), sum(!is.na(unlist(plot_data))), length(unlist(plot_data))))
