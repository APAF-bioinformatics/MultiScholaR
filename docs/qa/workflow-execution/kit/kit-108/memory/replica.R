# KIT-108 W2, instrument I3: the in-process replica of the report step's two calls (D7).
#
# On a copy of a completed member's project tree, in one R process, it runs what the summary tab runs in the
# application's controller after the save:
#   copyToResultsSummary(omic_type, label, contrasts_tbl, design_matrix, force = TRUE)
#   RenderReport(omic_type, label, rmd_filename)
# with project_dirs built by setupDirectories(base_dir, omic, label = NULL, reuse_existing = TRUE) (the application's
# "reuse existing" path, app_server.R:363-369), assigned into .GlobalEnv as the handler does when absent, and the
# contrasts and design read from source_dir as the handler's fallback does (mod_lipid_summary_workflow_helpers.R:234-247).
# time_dir/timestamp are pointed at the member's own filtering_qc/<timestamp> directory (setupDirectories mints a new
# one; the copy reads time_dir for the filtering grids).
#
# Around each call: gc(reset = TRUE) before; after, /proc/self/status raw, then gc() (used and max used since the reset),
# /proc/self/status and mallinfo2() again; an Rprof memory profile over the call. At the end: lobstr::obj_size of every
# .GlobalEnv object and of the calls' return values, then one malloc_trim(0) probe with /proc/self/status after it.
# Every row is appended and flushed as it is taken; on any failure on.exit writes a final `failed` row.
#
# Usage (R_LIBS = the build's library; TMPDIR on disk):
#   Rscript replica.R --project <member>/project --out <dir> [--omic lipidomics] [--label <experiment label>]
#       [--template lipidomics_report.rmd] [--template-dir <dir>] [--remove-template] [--keep-work] [--load-all <tree>]
#   --project       the member's project directory (one experiment directory inside), or the member directory
#   --label         default: the driver request's journey$experiment_label next to the project, else "replica"
#   --template-dir  where the template is looked up when scripts/<omic>/ has none (default: the package's reports/<dir>)
#   --remove-template  the failure probe: delete scripts/<omic>/<template> from the copy before the render
#   --keep-work     keep the copied tree (default: removed at the end, after the rows are written)
#   --load-all      load the package from a source tree with pkgload instead of the installed build (the W4/W5 trees)
# Exit 0 when both calls complete; 1 otherwise (rows up to the failure plus a final `failed` row are written).

args <- commandArgs(TRUE)
opt <- function(name, default = NULL) {
    i <- match(paste0("--", name), args)
    if (is.na(i)) return(default)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) return(TRUE)
    args[[i + 1L]]
}
project <- normalizePath(opt("project"), mustWork = TRUE)
if (dir.exists(file.path(project, "project"))) project <- file.path(project, "project")
out <- opt("out")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
rowsPath <- file.path(out, "replica-rows.tsv")

state <- new.env()
state$phase <- "init"; state$done <- FALSE; state$error <- NA_character_; state$gc <- NULL; state$t0 <- Sys.time()
allocator <- NULL

procStatus <- function() {
    lines <- readLines("/proc/self/status")
    pick <- function(key) {
        l <- grep(paste0("^", key, ":"), lines, value = TRUE)
        if (!length(l)) return(NA_real_)
        as.numeric(strsplit(trimws(sub(paste0("^", key, ":"), "", l)), "[[:space:]]+")[[1L]][[1L]])
    }
    c(vm_rss_kb = pick("VmRSS"), rss_anon_kb = pick("RssAnon"), rss_file_kb = pick("RssFile"), vm_hwm_kb = pick("VmHWM"),
        vm_data_kb = pick("VmData"))
}
uptime <- function() as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]])
allocatorValues <- function(trim = FALSE) {
    blank <- c(arena_bytes = NA, in_use_arena_bytes = NA, free_arena_bytes = NA, mmap_bytes = NA,
        main_top_free_bytes = NA, free_chunk_count = NA, trim_return = NA)
    if (is.null(allocator)) return(blank)
    v <- tryCatch(allocator(trim), error = function(e) blank)
    v[names(blank)]
}
gcValues <- function(g) {
    if (is.null(g)) return(c(gc_used_ncells_mb = NA, gc_used_vcells_mb = NA, gc_max_ncells_mb = NA, gc_max_vcells_mb = NA))
    # gc()'s matrix: rows Ncells, Vcells; column 2 = used (Mb), the "max used" (Mb) column is the 6th.
    c(gc_used_ncells_mb = g[1, 2], gc_used_vcells_mb = g[2, 2], gc_max_ncells_mb = g[1, 6], gc_max_vcells_mb = g[2, 6])
}
row <- function(phase, status = "ok", g = NULL, trim = FALSE, seconds = NA_real_, note = "") {
    a <- allocatorValues(trim)
    s <- procStatus()
    values <- c(list(epoch = sprintf("%.3f", as.numeric(Sys.time())), uptime = sprintf("%.2f", uptime()), phase = phase,
        status = status, seconds = round(seconds, 3)), as.list(s), as.list(round(gcValues(g), 1)), as.list(a),
        list(note = gsub("[\t\n]", " ", note)))
    first <- !file.exists(rowsPath)
    utils::write.table(as.data.frame(values, stringsAsFactors = FALSE), rowsPath, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = first, append = !first)
    invisible(NULL)
}

profileSummary <- function(file, name) {
    if (!file.exists(file) || file.info(file)$size == 0) return(invisible(NULL))
    p <- tryCatch(utils::summaryRprof(file, memory = "both"), error = function(e) NULL)
    if (is.null(p) || !nrow(p$by.total)) return(invisible(NULL))
    # The replica's own frames carry every sample; drop them so the ranking starts at the package's code.
    wrappers <- c("main", "measured", "force", "tryCatch", "tryCatchList", "tryCatchOne", "doTryCatch",
        "withCallingHandlers", "ns$copyToResultsSummary", "ns$RenderReport")
    tot <- p$by.total
    tot <- tot[!gsub('"', "", rownames(tot)) %in% wrappers, , drop = FALSE]
    tot <- tot[order(-tot$mem.total), , drop = FALSE]
    tot <- cbind(fn = rownames(tot), tot)
    utils::write.table(utils::head(tot, 40L), file.path(out, paste0("rprof-", name, "-by-total.tsv")), sep = "\t",
        quote = FALSE, row.names = FALSE)
    self <- p$by.self
    self <- self[order(-self$mem.total), , drop = FALSE]
    self <- cbind(fn = rownames(self), self)
    utils::write.table(utils::head(self, 40L), file.path(out, paste0("rprof-", name, "-by-self.tsv")), sep = "\t",
        quote = FALSE, row.names = FALSE)
    # file:line when the package has srcrefs (load_all trees); an installed build gives function names only.
    lines <- tryCatch(utils::summaryRprof(file, memory = "both", lines = "show")$by.self, error = function(e) NULL)
    if (!is.null(lines) && nrow(lines)) {
        lines <- lines[order(-lines$mem.total), , drop = FALSE]
        utils::write.table(cbind(where = rownames(lines), utils::head(lines, 40L)),
            file.path(out, paste0("rprof-", name, "-by-line.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    invisible(NULL)
}

measured <- function(name, expr) {
    state$phase <- name
    g0 <- gc(reset = TRUE)
    row(paste0("before_", name), g = g0)
    prof <- file.path(out, paste0("rprof-", name, ".out"))
    utils::Rprof(prof, interval = 0.05, memory.profiling = TRUE, line.profiling = TRUE)
    started <- Sys.time()
    value <- tryCatch(force(expr), finally = {
        utils::Rprof(NULL)
        state$elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    })
    row(paste0("after_", name, "_raw"), seconds = state$elapsed, note = "no explicit gc yet")
    g1 <- gc()
    row(paste0("after_", name), g = g1, seconds = state$elapsed, note = "after gc(); max used = peak since the reset")
    profileSummary(prof, name)
    value
}

main <- function() {
    # The failure path: whatever stops the run, the rows taken so far stay and a final `failed` row is appended.
    on.exit({
        if (!isTRUE(state$done)) {
            try(utils::Rprof(NULL), silent = TRUE)
            try(row("end", status = "failed", note = paste0("in ", state$phase, ": ", state$error)), silent = TRUE)
            state$done <- TRUE
        }
    }, add = TRUE)
    row("start", note = paste("R", getRversion(), "pid", Sys.getpid(), "libs", paste(.libPaths(), collapse = ":")))
    # mallinfo2 and malloc_trim through the probe of tools/profiling/workflow_checkpoint_gui_diagnostic.R:41-60,
    # compiled in this process (scratch cache under TMPDIR).
    state$phase <- "probe_compile"
    native <- new.env(parent = baseenv())
    ok <- tryCatch({
        Rcpp::cppFunction(code = '
Rcpp::NumericVector allocatorSnapshot(bool trim) {
    int released = trim ? malloc_trim(0) : -1;
    struct mallinfo2 info = mallinfo2();
    return Rcpp::NumericVector::create(
        Rcpp::Named("arena_bytes") = (double) info.arena,
        Rcpp::Named("in_use_arena_bytes") = (double) info.uordblks,
        Rcpp::Named("free_arena_bytes") = (double) info.fordblks,
        Rcpp::Named("mmap_bytes") = (double) info.hblkhd,
        Rcpp::Named("main_top_free_bytes") = (double) info.keepcost,
        Rcpp::Named("free_chunk_count") = (double) info.ordblks,
        Rcpp::Named("trim_return") = released);
}', includes = "#include <malloc.h>", env = native, cacheDir = file.path(tempdir(), "kit108-probe"))
        TRUE
    }, error = function(e) conditionMessage(e))
    if (isTRUE(ok)) allocator <<- native$allocatorSnapshot
    row("probe_compiled", status = if (isTRUE(ok)) "ok" else "warn", note = if (isTRUE(ok)) "mallinfo2 available" else ok)

    state$phase <- "load"
    # --load-all <tree>: a source tree under pkgload (the W4/W5 trees, D8); srcrefs give the profile file:line.
    tree <- opt("load-all", NULL)
    if (!is.null(tree)) suppressPackageStartupMessages(pkgload::load_all(tree, quiet = TRUE, export_all = FALSE))
    else suppressPackageStartupMessages(loadNamespace("MultiScholaR"))
    ns <- asNamespace("MultiScholaR")
    row("package_loaded", note = paste("MultiScholaR", utils::packageVersion("MultiScholaR"), dirname(getNamespaceInfo(ns, "path"))))

    state$phase <- "copy_tree"
    experiments <- list.dirs(project, recursive = FALSE)
    if (length(experiments) != 1L) stop("expected one experiment directory under ", project, ", found ", length(experiments))
    work <- file.path(out, "work")
    unlink(work, recursive = TRUE)
    dir.create(work)
    t <- Sys.time()
    status <- system2("cp", c("-a", shQuote(experiments), shQuote(work)))
    if (!identical(status, 0L)) stop("cp -a failed with status ", status)
    base <- file.path(work, basename(experiments))
    omic <- opt("omic", basename(list.dirs(file.path(base, "results"), recursive = FALSE))[[1L]])
    reportDir <- omic
    template <- opt("template", paste0(omic, "_report.rmd"))
    label <- opt("label", NULL)
    if (is.null(label)) {
        request <- file.path(dirname(project), "driver-request.rds")
        label <- if (file.exists(request)) readRDS(request)$journey$experiment_label else "replica"
    }
    row("tree_copied", seconds = as.numeric(difftime(Sys.time(), t, units = "secs")),
        note = sprintf("omic=%s label=%s template=%s base=%s", omic, label, template, base))

    state$phase <- "project_dirs"
    dirs <- ns$setupDirectories(base_dir = base, omic_types = omic, label = NULL, force = FALSE, reuse_existing = TRUE)
    stamps <- list.dirs(dirs[[omic]]$qc_dir, recursive = FALSE, full.names = FALSE)
    stamps <- setdiff(stamps[grepl("^[0-9]{8}_[0-9]{6}$", stamps)], dirs[[omic]]$timestamp)
    if (length(stamps) == 1L) {
        unlink(dirs[[omic]]$time_dir, recursive = TRUE)
        dirs[[omic]]$timestamp <- stamps
        dirs[[omic]]$time_dir <- file.path(dirs[[omic]]$qc_dir, stamps)
    }
    if (!exists("project_dirs", envir = .GlobalEnv)) assign("project_dirs", dirs, envir = .GlobalEnv)
    source <- dirs[[omic]]$source_dir
    contrasts <- if (file.exists(file.path(source, "contrasts_tbl.tab"))) readr::read_tsv(file.path(source, "contrasts_tbl.tab"), show_col_types = FALSE)
    design <- if (file.exists(file.path(source, "design_matrix.tab"))) readr::read_tsv(file.path(source, "design_matrix.tab"), show_col_types = FALSE)
    row("project_dirs_set", note = sprintf("setupDirectories(reuse_existing = TRUE); time_dir=%s; contrasts_tbl.tab %s; design_matrix.tab %s",
        basename(dirs[[omic]]$time_dir), if (is.null(contrasts)) "absent" else "read", if (is.null(design)) "absent" else "read"))

    copyValue <- measured("copy", ns$copyToResultsSummary(omic_type = omic, experiment_label = label,
        contrasts_tbl = contrasts, design_matrix = design, force = TRUE))

    # The report handler's template step (mod_lipid_summary_workflow_helpers.R:398-447): the project's
    # scripts/<omic>/ copy, else the package's reports/<dir>/ copy (--template-dir overrides the lookup).
    state$phase <- "template"
    target <- file.path(base, "scripts", omic, template)
    if (isTRUE(opt("remove-template", FALSE)) && file.exists(target)) file.remove(target)
    if (!file.exists(target)) {
        lookup <- opt("template-dir", system.file("reports", reportDir, package = "MultiScholaR"))
        candidate <- file.path(lookup, template)
        if (nzchar(lookup) && file.exists(candidate)) {
            file.copy(candidate, target)
            row("template_ready", note = paste("copied from", candidate))
        } else {
            row("template_missing", status = "warn", note = paste("not in scripts/ and not at", candidate))
        }
    } else row("template_ready", note = paste("already at", target))

    renderValue <- measured("render", ns$RenderReport(omic_type = omic, experiment_label = label, rmd_filename = template))

    state$phase <- "object_sizes"
    sizes <- data.frame(where = character(), object = character(), bytes = numeric(), stringsAsFactors = FALSE)
    for (name in ls(.GlobalEnv, all.names = TRUE)) {
        where <- if (name %in% replicaNames) "replica-own" else ".GlobalEnv"
        sizes[nrow(sizes) + 1L, ] <- list(where, name, as.numeric(lobstr::obj_size(get(name, envir = .GlobalEnv))))
    }
    sizes[nrow(sizes) + 1L, ] <- list("return", "copyToResultsSummary", as.numeric(lobstr::obj_size(copyValue)))
    sizes[nrow(sizes) + 1L, ] <- list("return", "RenderReport", as.numeric(lobstr::obj_size(renderValue)))
    sizes[nrow(sizes) + 1L, ] <- list("replica", "contrasts_tbl+design_matrix", as.numeric(lobstr::obj_size(contrasts, design)))
    utils::write.table(sizes, file.path(out, "object-sizes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    app <- sizes[sizes$where != "replica-own", ]
    row("object_sizes", note = sprintf("%d objects outside the replica's own; largest %s (%.1f MB)", nrow(app),
        app$object[which.max(app$bytes)], max(app$bytes) / 1e6))

    state$phase <- "final"
    row("final_gc", g = gc(), note = "full gc after both calls")
    row("after_malloc_trim", trim = TRUE, note = "malloc_trim(0) called; trim_return 1 = memory was released")
    if (!isTRUE(opt("keep-work", FALSE))) unlink(work, recursive = TRUE)
    state$done <- TRUE
    row("end", status = "complete", seconds = as.numeric(difftime(Sys.time(), state$t0, units = "secs")),
        note = paste("report:", if (is.character(renderValue)) renderValue else "<non-character>"))
}

replicaNames <- c(ls(.GlobalEnv, all.names = TRUE), "replicaNames", "result")
result <- tryCatch({
    withCallingHandlers(main(), error = function(e) state$error <- conditionMessage(e))
    0L
}, error = function(e) {
    message("REPLICA FAILED in ", state$phase, ": ", conditionMessage(e))
    1L
})
quit(save = "no", status = result)
