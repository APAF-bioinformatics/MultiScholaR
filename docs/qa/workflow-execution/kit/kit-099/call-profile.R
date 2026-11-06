# KIT-099 call profile (declaration section 5, item 6): in process on kit_098, not a member. Two parts.
# (a) The filter and TMM/logCPM declarations, `whole` and `partitioned`, each on a fresh reflinked copy of its large pair
#     population store (trx99lf_design, trx99lq_qc_finalize; recorded to the input step, so the product is computed). Two
#     phases each: the pre-pass, from the start to the first per-row edgeR call (`filterByExpr` for the filter, `cpm` for
#     logCPM), and the rest. Each phase has its time and its peak RSS (`VmHWM` of this process, reset at each phase start
#     through /proc/self/clear_refs). The partitioned pre-pass is the TMM-effective library size (declaration section 7).
# (b) `Rprof` with memory profiling over a whole in-process run of the KIT-099 recorder at large (ticket 099l's settings on
#     trx_60k_200, label trx99lp, its record written to the runtime only): total time and allocation per candidate call.
# Usage, from a clean worktree at the bound commit: KIT099_INSTALLED=1 R_LIBS=<kit_098 library> Rscript \
#   docs/qa/workflow-execution/kit/kit-099/call-profile.R <scratch dir> <output txt>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE)
scratch <- args[[1L]]; output <- args[[2L]]
repo <- normalizePath(".")
ns <- asNamespace("MultiScholaR")
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
lines <- character()
say <- function(...) { line <- sprintf(...); lines <<- c(lines, line); cat(line, "\n") }
hwm <- function() as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE))) * 1024
reset_hwm <- function() { gc(); writeLines("5", "/proc/self/clear_refs") }
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
say("# KIT-099 call profile, %s, library %s, build %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), find.package("MultiScholaR"),
    utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit`)

# (a) the two declarations by phase.
mark <- new.env()
marker <- function() if (is.null(mark$t1)) { mark$t1 <- Sys.time(); mark$hwm1 <- hwm(); writeLines("5", "/proc/self/clear_refs") }
phase_case <- function(operation, population_key, variant, per_row) {
    population <- bindings$populations[[population_key]]
    stopifnot(!is.null(population), is.null(population$head_entity))
    dest <- file.path(scratch, paste(operation, variant, sep = "-"))
    unlink(dest, recursive = TRUE); dir.create(dest)
    stopifnot(identical(system2("cp", c("-a", "--reflink=auto", shQuote(population$store), shQuote(file.path(dest, ".workflow-entities")))), 0L))
    identity <- readSavedStoreIdentity(file.path(dest, ".workflow-entities", "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(dest, identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    stopifnot(identical(manager$.context$snapshot$current_state, population$state))
    declaration <- workflowOperationEntry(operation)$declaration
    rm(list = ls(mark), envir = mark)
    suppressMessages(trace(per_row, tracer = quote(marker()), where = asNamespace("edgeR"), print = FALSE))
    on.exit(suppressMessages(untrace(per_row, where = asNamespace("edgeR"))), add = TRUE)
    reset_hwm(); t0 <- Sys.time()
    run <- withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager)
        runSavedOperation(manager, declaration, source, list(variant = variant, worker_cap = 12 * 1024^3))
    })
    t2 <- Sys.time(); hwm2 <- hwm()
    stopifnot(identical(run$status, "computed"), !is.null(mark$t1))
    say("%-22s %-11s | pre-pass %6.1f s peak RSS %5.2f GiB | rest %6.1f s peak RSS %5.2f GiB | status %s", operation, variant,
        as.numeric(difftime(mark$t1, t0, units = "secs")), mark$hwm1 / 1024^3, as.numeric(difftime(t2, mark$t1, units = "secs")), hwm2 / 1024^3, run$status)
}
say("## (a) pre-pass and batched phases, in process, each on a fresh copy of its pair population store")
for (variant in c("whole", "partitioned")) phase_case("trx_expression_filter", "trx99lf_design", variant, "filterByExpr")
for (variant in c("whole", "partitioned")) phase_case("trx_tmm_logcpm", "trx99lq_qc_finalize", variant, "cpm")

# (b) Rprof over a whole in-process run of the recorder at large.
say("## (b) Rprof (interval 0.05 s, memory profiling) over the KIT-099 recorder at large (label trx99lp, ticket 099l, trx_60k_200)")
profile <- file.path(scratch, "rprof.out")
record <- file.path(scratch, "trx99lp.rds")
Sys.setenv(KIT099_INSTALLED = "1")
recorder_args <- c("trx99lp", "099l", record, "trx_60k_200")
started <- Sys.time()
Rprof(profile, interval = 0.05, memory.profiling = TRUE)
local({
    commandArgs <- function(trailingOnly = FALSE) recorder_args
    sys.source(file.path(repo, "docs/qa/workflow-execution/kit/kit-099/oracle-record.R"), envir = environment())
})
Rprof(NULL)
say("recorder run %.0f s", as.numeric(difftime(Sys.time(), started, units = "secs")))
summary <- summaryRprof(profile, memory = "both")$by.total
candidates <- c("svd", "La.svd", "cancor", "ruv_cancorplot", "ruv::ruv_cancorplot", "trxGridMember", "trxPrelimOrder", "estimateDisp", "glmQLFit",
    "glmQLFTest", "RUVs", "calcNormFactors", "filterByExpr", "cpm", "trxRuvGrid", "trxRuvs", "trxDaFit", "trxDaApply", "RUVIII",
    "runTrxReferenceRuvs", "runTrxReferenceDaAnalysis", "runTrxReferenceExpressionFilter", "runTrxReferenceRuviiiExport")
rows <- summary[gsub('"', "", rownames(summary)) %in% candidates, , drop = FALSE]
rows <- rows[order(-rows$total.time), , drop = FALSE]
say("%-34s %10s %8s %12s", "call", "total s", "total %", "alloc MB")
for (i in seq_len(nrow(rows))) say("%-34s %10.1f %8.1f %12.0f", gsub('"', "", rownames(rows)[[i]]), rows$total.time[[i]], rows$total.pct[[i]], rows$mem.total[[i]])
writeLines(lines, output)
