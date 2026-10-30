# KIT-071: read a lifecycle member's owner store and export journal and print the closure facts:
# the head's current state and history, the provenance kinds (a DA Open that reopened the owner adds
# none; an adoption would add a `legacy` entry), the export journal's state, and the exported session's
# state name against the head's. Usage: Rscript verify-member.R <member project dir> <lipid|metab>
args <- commandArgs(trailingOnly = TRUE)
project <- args[[1]]; omic <- args[[2]]
dir <- list.files(project, pattern = if (omic == "lipid") "^lipidomics$" else "^metabolomics$", recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
dir <- dir[grepl("/scripts/", dir) & dir.exists(dir)][[1]]
store <- list.files(dir, pattern = "^\\.dia-native-", all.files = TRUE, full.names = TRUE, include.dirs = TRUE)
suppressPackageStartupMessages(library(MultiScholaR))
ns <- asNamespace("MultiScholaR")
entities <- file.path(store, ".workflow-entities")
identity <- ns$readSavedStoreIdentity(file.path(entities, "store.json"))
opened <- ns$newWorkflowEntityStore(store, identity$project_id, identity$workflow_id, max_cache_bytes = 0)
snapshot <- ns$hydrateSavedWorkflowCheckpoint(opened, readRDS(file.path(entities, "dia-state-head.rds")))
ns$closeWorkflowEntityStore(opened)
cat("store:", basename(store), "\n")
cat("head current_state:", snapshot$current_state, "\n")
cat("head history:", paste(unlist(snapshot$state_history), collapse = " > "), "\n")
prov <- snapshot$provenance %||% snapshot$step_provenance
if (!is.null(prov)) cat("provenance entries:", paste(vapply(prov, function(p) paste0(p$kind %||% p$operation %||% p$step %||% "?", if (isTRUE(p$adopted)) "(adopted)" else ""), ""), collapse = ", "), "\n")
cat("snapshot fields:", paste(names(snapshot), collapse = ", "), "\n")
head <- snapshot
journal <- jsonlite::fromJSON(list.files(dir, "^export_complete_.*\\.json$", full.names = TRUE)[[1]])
cat("export status:", journal$status, "session:", journal$session_name, "latest:", journal$latest_name, "\n")
latest <- readRDS(file.path(dir, journal$latest_name))
cat("exported r6_current_state_name:", latest$r6_current_state_name, " equals head:", identical(latest$r6_current_state_name, head$current_state), "\n")
cat("exported sha256 equals journal:", identical(unname(as.character(digest::digest(file.path(dir, journal$latest_name), algo = "sha256", file = TRUE))), journal$latest_sha256), "\n")
rds <- list.files(project, "\\.(rds|RDS)$", recursive = TRUE, full.names = TRUE)
rds <- rds[!grepl("\\.workflow-(entities|export)|dia-native", rds)]
cat("RDS products:\n"); print(data.frame(file = sub(paste0(project, "/"), "", rds), bytes = file.size(rds)), row.names = FALSE)
