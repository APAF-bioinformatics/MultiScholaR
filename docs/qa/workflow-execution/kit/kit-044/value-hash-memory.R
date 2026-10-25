# KIT-044: the memory cost of the in-memory value hash (the fallback without bash/dd/sha256sum) against the streamed
# GNU path, on a saved peptide table's columns: each column hashed both ways in one fresh process per mode, VmHWM
# and time recorded, the hex compared. Usage: Rscript docs/qa/workflow-execution/kit/kit-044/value-hash-memory.R <project> <state> <mode: gnu|memory> <out json>
suppressMessages(devtools::load_all(quiet = TRUE))
ns <- asNamespace("MultiScholaR"); api <- function(n) get(n, ns); h <- new.env(); sys.source("tools/profiling/workflow_kit.R", h)
a <- commandArgs(TRUE); project <- a[[1]]; state <- a[[2]]; mode <- a[[3]]; out <- a[[4]]
if (identical(mode, "memory")) Sys.setenv(MULTISCHOLAR_DISABLE_TOOLS = "bash,dd,sha256sum")
rss <- function() as.numeric(strsplit(grep("VmHWM", readLines("/proc/self/status"), value = TRUE), "[[:space:]]+")[[1]][2]) / 2^20
store <- h$kitFindStore(project); root <- dirname(store); id <- api("readSavedStoreIdentity")(file.path(root, ".workflow-entities", "store.json"))
st <- api("newWorkflowEntityStore")(root, id$project_id, id$workflow_id); m <- api("newDIAStateAdapter")(st, reopen = TRUE)
result <- api("withDIAStateScope")(m, function() {
    ref <- api("getDIAStateReference")(m, state); before <- api("readDIAStateMetadata")(m, ref)
    scope <- api("newWorkflowTableReadScope")(st, before$data$components$peptide_data); on.exit(api("releaseWorkflowTableReadScope")(scope))
    rows <- scope$metadata$rows; columns <- scope$metadata$columns
    base <- rss(); hashes <- list(); seconds <- list(); peaks <- list()
    for (column in columns) {
        pages <- list()
        api("streamWorkflowTableRows")(scope, column, function(page, offset, row_names) { pages[[length(pages) + 1L]] <<- page[[column]]; invisible(NULL) })
        value <- do.call(c, pages); pages <- NULL
        attributes(value) <- NULL
        t0 <- Sys.time(); hashes[[column]] <- api("hashWorkflowSerializedValue")(list(value = value, names = NULL), st)
        seconds[[column]] <- as.numeric(difftime(Sys.time(), t0, units = "secs")); peaks[[column]] <- rss(); value <- NULL; gc()
    }
    list(mode = mode, contracts = api("workflowCapabilityContracts")(), rows = rows, columns = columns, base_hwm_gib = base, hashes = hashes, seconds = seconds, hwm_after_gib = peaks, hwm_gib = rss()) })
jsonlite::write_json(result, out, auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(mode, "| rows", result$rows, "| columns", length(result$columns), "| HWM", round(result$hwm_gib, 2), "GiB | base", round(result$base_hwm_gib, 2), "| seconds", round(sum(unlist(result$seconds))), "\n")
