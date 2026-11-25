# KIT-259 diagnostic G0 (triage of KIT-247 G0-4, run 45604291093): why the label-free normalisation's Pearson
# correlation figure fails with "Problem while converting geom to grob" in the browser lane. No package code.
# Binding: v0.9-staging 63143aa3 plus this folder (the drawing code is unchanged from 16112917, the lane's binding).
#
# Part 1 (plain R process, shiny::testServer, as KIT-247's g0-lanes.R and helper-gui-entry.R drive it): lanes prot_lfq
#   (MaxQuant label-free, 400 x 10) and prot_dia (DIA-NN, 8 runs, the control) up to their "norm" step. Spies on
#   saveProtNormQcPlotArtifact (status and full condition chain of every figure save; the plot object saved as RDS)
#   and on buildProtNormQcPlotProduct (the S4 object each correlation figure is drawn from, saved as RDS).
# Part 2 (the same plain process, outside any Shiny session, warn = 0): each captured S4 object, through the
#   application's own producers called directly (buildProtNormQcPlotProduct for pca, rle, density, correlation; and
#   plotPearson's native plot), rendered by ggsave (the device it picks for .png) and by grDevices::png(type = "cairo");
#   then the correlation figures again under options(warn = 2). Data facts of each correlation figure.
# Part 3 (launch mode): the application started by shinytest2 exactly as the browser lane starts it (run_app(test_mode
#   = TRUE) in shinytest2's background process, Chrome), with a probe that renders the captured correlation figures
#   with ggsave in that process: once before the app object is built (app_top) and once inside an observer under
#   shiny::withProgress (app_observer), as the tab does. The process's options and devices are recorded.
#
# Usage (package root, on a runner, NOT_CRAN=true, shinytest2 installed): Rscript <this file> <out dir>

k259Args <- commandArgs(trailingOnly = TRUE)
k259Out <- if (length(k259Args) >= 1L) k259Args[[1L]] else "kit259-g0"
dir.create(k259Out, recursive = TRUE, showWarnings = FALSE)
k259Out <- normalizePath(k259Out)
k259Rds <- file.path(tempdir(), "k259-rds")
dir.create(k259Rds, showWarnings = FALSE)

# ---------------------------------------------------------------------------------------------------------------------
# PREDICTIONS (written before any value is read; regular expressions applied to the observed value). By reading at
# 63143aa3: the histogram's fixed y limits c(0, 4) (rebuildWorkflowDiaQcPlot, utils_workflow_dia_products.R:84;
# plotPearson, func_prot_s4_qc_correlation_methods.R:68) with 0.001-wide bins put any bin holding more than 4 pairs
# outside the scale; its bar's ymax becomes NA and ggplot2 removes it with a warning. A Python estimate of the
# label-free fixture's raw log2 within-group correlations gave 20 pairs with a fullest bin of 5; the DIA-NN fixture has
# 12 pairs. The frozen gui-entry oracle (tests/testdata/gui-entry/prot_lfq.json) lists pre_norm_correlation.png and
# post_norm_correlation.png under testServer, so testServer draws them. The best-guess mechanism for the browser-only
# error is that the warning is promoted to an error in the application's process (warn >= 2).
# ---------------------------------------------------------------------------------------------------------------------
k259P <- function(key, predicted, pattern) data.frame(key = key, predicted = predicted, pattern = pattern,
    stringsAsFactors = FALSE)
k259Predictions <- rbind(
    k259P("p1.prot_lfq.saved_ok", "testServer: the 8 pre/post figures of the label-free lane are saved without error", "^8/8$"),
    k259P("p1.prot_dia.saved_ok", "testServer: the 8 pre/post figures of the DIA-NN lane are saved without error", "^8/8$"),
    k259P("p1.prot_lfq.norm_head", "the label-free normalisation completes under testServer", "^normalized$"),
    k259P("p1.prot_dia.norm_head", "the DIA-NN normalisation completes under testServer", "^normalized$"),
    k259P("p1.prot_lfq.norm_notes", "with its completion notice", "Normalization completed"),
    k259P("p1.warn_option", "the testServer lane runs at warn = 1 (helper-gui-entry.R:277)", "^1$"),
    k259P("facts.prot_lfq.pre_norm.n_pairs", "20 within-group pairs (5 + 5 samples)", "^20$"),
    k259P("facts.prot_lfq.post_norm.n_pairs", "20 within-group pairs", "^20$"),
    k259P("facts.prot_dia.pre_norm.n_pairs", "12 within-group pairs (4 + 4 runs)", "^12$"),
    k259P("facts.prot_lfq.pre_norm.max_bin_over_4", "a 0.001 bin holds more than 4 pairs", "^TRUE$"),
    k259P("facts.prot_lfq.post_norm.max_bin_over_4", "also after normalisation", "^TRUE$"),
    k259P("facts.prot_dia.pre_norm.max_bin_over_4", "no DIA-NN bin holds more than 4", "^FALSE$"),
    k259P("facts.prot_dia.post_norm.max_bin_over_4", "no DIA-NN bin holds more than 4", "^FALSE$"),
    k259P("facts.prot_lfq.pre_norm.censored_bars", "ggplot2 sets at least one bar's ymax to NA", "^[1-9]"),
    k259P("facts.prot_dia.pre_norm.censored_bars", "none for DIA-NN", "^0$"),
    k259P("p2.prot_lfq.ok", "plain process, warn 0: every label-free figure renders on both devices", "^20/20$"),
    k259P("p2.prot_dia.ok", "plain process, warn 0: every DIA-NN figure renders on both devices", "^20/20$"),
    k259P("p2.prot_lfq.pre_norm.correlation.ggsave.warnings", "with the removed-rows warning", "outside the scale range"),
    k259P("p2w2.prot_lfq.pre_norm.correlation.status", "warn 2: the label-free figure fails as in the browser",
        "Problem while converting geom to grob.*converted from warning.*outside the scale range"),
    k259P("p2w2.prot_lfq.post_norm.correlation.status", "warn 2: and after normalisation",
        "Problem while converting geom to grob.*converted from warning.*outside the scale range"),
    k259P("p2w2.prot_dia.pre_norm.correlation.status", "warn 2: the DIA-NN figure renders", "^ok$"),
    k259P("p3.app_observer.env.warn", "the application's process under shinytest2 runs at warn >= 2", "^[2-9]"),
    k259P("p3.app_observer.prot_lfq-pre_norm_correlation.png.status", "in the app's observer the label-free figure fails",
        "^error: .*Problem while converting geom to grob"),
    k259P("p3.app_observer.prot_lfq-post_norm_correlation.png.status", "and the post-normalisation one",
        "^error: .*Problem while converting geom to grob"),
    k259P("p3.app_observer.prot_dia-pre_norm_correlation.png.status", "the DIA-NN figure renders in the app", "^ok$"),
    k259P("p3.app_top.prot_lfq-pre_norm_correlation.png.status", "before the app object is built it fails too",
        "^error: .*Problem while converting geom to grob"))

# ---------------------------------------------------------------------------------------------------------------------
k259Log <- function(...) cat(sprintf("[KIT-259 G0 %s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)), file = stderr())
k259Rows <- new.env(parent = emptyenv())
k259Rows$rows <- list()
source(file.path("docs", "qa", "workflow-execution", "kit", "kit-259", "k259-probe.R"))
k259Put <- function(key, value) {
    value <- k259Text(value)
    k259Rows$rows[[length(k259Rows$rows) + 1L]] <- data.frame(key = key, value = value, stringsAsFactors = FALSE)
    k259Log(key, " = ", substr(value, 1L, 400L))
}
k259Write <- function() {
    extra <- list.files(k259Out, pattern = "^p3-.*[.]tsv$", full.names = TRUE)
    observed <- do.call(rbind, c(k259Rows$rows, lapply(extra, function(f) utils::read.delim(f, quote = "",
        stringsAsFactors = FALSE, colClasses = "character"))))
    if (is.null(observed)) observed <- data.frame(key = character(), value = character())
    utils::write.table(observed, file.path(k259Out, "g0-observed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    compare <- k259Predictions
    compare$observed <- vapply(compare$key, function(key) {
        hit <- observed$value[observed$key == key]
        if (length(hit)) hit[[length(hit)]] else "<not observed>"
    }, "")
    compare$held <- ifelse(compare$observed == "<not observed>", NA,
        mapply(function(p, o) grepl(p, o, perl = TRUE), compare$pattern, compare$observed))
    utils::write.table(compare, file.path(k259Out, "g0-compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    k259Log("written: ", sum(compare$held %in% TRUE), " held, ", sum(compare$held %in% FALSE), " not held, ",
        sum(is.na(compare$held)), " not observed")
}

k259Log("loading the package and the helpers")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
k259Env <- k259Environment()
for (name in names(k259Env)) k259Put(paste0("env.", name), k259Env[[name]])

# ---------------------------------------------------------------------------------------------------------------------
# Part 1: the lanes under testServer up to "norm", with the spies.
# ---------------------------------------------------------------------------------------------------------------------
k259Part1 <- function(lane_id) {
    ns <- asNamespace("MultiScholaR")
    o_save <- get("saveProtNormQcPlotArtifact", ns)
    o_build <- get("buildProtNormQcPlotProduct", ns)
    spy <- new.env(parent = emptyenv())
    spy$s4 <- NULL
    spy$saves <- character()
    testthat::local_mocked_bindings(
        buildProtNormQcPlotProduct = function(s4Object, aesthetics, plotType, ...) {
            if (identical(plotType, "correlation")) spy$s4 <- s4Object
            o_build(s4Object, aesthetics, plotType, ...)
        },
        saveProtNormQcPlotArtifact = function(qcDir, filename, plotObject, width, height, dpi = 150,
                                              savePlotFn = ggplot2::ggsave) {
            id <- paste0(lane_id, "-", filename)
            saveRDS(plotObject, file.path(k259Rds, paste0(id, ".rds")))
            if (grepl("_correlation[.]png$", filename) && !is.null(spy$s4)) {
                saveRDS(spy$s4, file.path(k259Rds, paste0(lane_id, "-", sub("_correlation[.]png$", "", filename), "-s4.rds")))
            }
            k259Put(paste0("p1.", id, ".warn"), getOption("warn"))
            value <- NULL
            got <- k259Capture(value <- o_save(qcDir, filename, plotObject, width, height, dpi, savePlotFn))
            k259Put(paste0("p1.", id, ".status"), if (identical(got$status, "ok")) "ok" else paste("error:", got$chain))
            k259Put(paste0("p1.", id, ".warnings"), got$warnings)
            if (!identical(got$status, "ok")) {
                k259Put(paste0("p1.", id, ".backtrace"), got$backtrace)
                k259Put(paste0("p1.", id, ".stored_traces"), got$stored_traces)
            }
            spy$saves <- c(spy$saves, paste0(filename, "=", got$status))
            if (!identical(got$status, "ok")) stop(got$error)
            value
        },
        .package = "MultiScholaR")
    # The lane's steps up to and including "norm" (the helper's own list, cut; nothing else changes).
    lane_env <- environment(guiEntryRunLane)
    o_steps <- get("guiEntrySteps", envir = lane_env)
    cut_steps <- function(lane) {
        steps <- o_steps(lane)
        ids <- vapply(steps, `[[`, "", "id")
        if (!"norm" %in% ids) stop("the lane has no norm step: ", paste(ids, collapse = ", "))
        steps[seq_len(match("norm", ids))]
    }
    locked <- bindingIsLocked("guiEntrySteps", lane_env)
    if (locked) unlockBinding("guiEntrySteps", lane_env)
    assign("guiEntrySteps", cut_steps, envir = lane_env)
    on.exit({
        assign("guiEntrySteps", o_steps, envir = lane_env)
        if (locked) lockBinding("guiEntrySteps", lane_env)
    }, add = TRUE)
    record <- guiEntryRunLane(lane_id)
    norm <- record$steps$norm
    k259Put(paste0("p1.", lane_id, ".steps"), paste(names(record$steps), collapse = ","))
    k259Put(paste0("p1.", lane_id, ".norm_head"), norm$head)
    k259Put(paste0("p1.", lane_id, ".norm_notes"), vapply(norm$notifications, function(n) paste0("[", n$type, "] ",
        n$message), ""))
    k259Put(paste0("p1.", lane_id, ".pre_qc_notes"), vapply(record$steps$norm_refusal$notifications, function(n)
        paste0("[", n$type, "] ", n$message), ""))
    pngs <- grep("_norm_(pca|rle|density|correlation)[.]png=", spy$saves, value = TRUE)
    k259Put(paste0("p1.", lane_id, ".saves"), spy$saves)
    k259Put(paste0("p1.", lane_id, ".saved_ok"), sprintf("%d/%d", sum(grepl("=ok$", pngs)), length(pngs)))
    invisible(record)
}

k259Put("p1.warn_option_outside", getOption("warn"))
for (lane_id in c("prot_lfq", "prot_dia")) {
    k259Log("part 1: lane ", lane_id)
    got <- k259Capture(k259Part1(lane_id))
    if (!identical(got$status, "ok")) k259Put(paste0("p1.", lane_id, ".lane_error"), got$chain)
    k259Write()
}
# The helper sets warn = 1 for the lane (withr::local_options); the value seen inside a figure save:
k259Put("p1.warn_option", {
    seen <- vapply(k259Rows$rows, function(r) if (grepl("^p1[.].*[.]warn$", r$key)) r$value else NA_character_, "")
    paste(unique(stats::na.omit(seen)), collapse = ",")
})

# ---------------------------------------------------------------------------------------------------------------------
# Part 2: the producers called directly in this plain process (no Shiny session), both devices, then warn = 2.
# ---------------------------------------------------------------------------------------------------------------------
k259Facts <- function(prefix, plot) {
    r <- plot$data$pearson_correlation
    k259Put(paste0(prefix, ".n_pairs"), length(r))
    k259Put(paste0(prefix, ".values"), paste(format(sort(r), digits = 6), collapse = ","))
    k259Put(paste0(prefix, ".n_nonfinite"), sum(!is.finite(r)))
    k259Put(paste0(prefix, ".n_above_1"), sum(r > 1, na.rm = TRUE))
    breaks <- seq(min(round(r - 0.5, 2), na.rm = TRUE), 1, 0.001)
    inside <- r[is.finite(r) & r >= min(breaks) & r <= max(breaks)]
    counts <- graphics::hist(inside, breaks = breaks, plot = FALSE)$counts
    k259Put(paste0(prefix, ".max_bin"), max(counts))
    k259Put(paste0(prefix, ".max_bin_over_4"), max(counts) > 4)
    built <- k259Capture(ld <- ggplot2::layer_data(plot, 1L))
    k259Put(paste0(prefix, ".layer_data"), built$status)
    if (identical(built$status, "ok")) {
        k259Put(paste0(prefix, ".censored_bars"), sum(is.na(ld$ymax) | is.na(ld$y)))
        k259Put(paste0(prefix, ".layer_count_max"), max(ld$count, na.rm = TRUE))
    }
    k259Put(paste0(prefix, ".build_warnings"), built$warnings)
}

k259Render <- function(prefix, plot) {
    target <- tempfile(fileext = ".png")
    a <- k259Capture(ggplot2::ggsave(target, plot, width = 10, height = 8, dpi = 150))
    k259Put(paste0(prefix, ".ggsave.status"), if (identical(a$status, "ok")) "ok" else paste("error:", a$chain))
    k259Put(paste0(prefix, ".ggsave.warnings"), a$warnings)
    if (!identical(a$status, "ok")) k259Put(paste0(prefix, ".ggsave.backtrace"), a$backtrace)
    b <- k259Capture({
        grDevices::png(target, width = 1500, height = 1200, res = 150, type = "cairo")
        print(plot)
        grDevices::dev.off()
    })
    while (grDevices::dev.cur() > 1L) grDevices::dev.off()
    k259Put(paste0(prefix, ".cairo.status"), if (identical(b$status, "ok")) "ok" else paste("error:", b$chain))
    if (!identical(b$status, "ok")) k259Put(paste0(prefix, ".cairo.backtrace"), b$backtrace)
    c(a$status, b$status)
}

aesthetics <- list(color_var = "group", shape_var = "group")
for (lane_id in c("prot_lfq", "prot_dia")) {
    statuses <- character()
    for (stage in c("pre_norm", "post_norm")) {
        path <- file.path(k259Rds, paste0(lane_id, "-", stage, "-s4.rds"))
        if (!file.exists(path)) {
            k259Put(paste0("p2.", lane_id, ".", stage, ".s4"), "<not captured>")
            next
        }
        s4 <- readRDS(path)
        k259Put(paste0("p2.", lane_id, ".", stage, ".s4"), sprintf("%s, %d x %d", class(s4)[[1L]],
            nrow(s4@protein_quant_table), ncol(s4@protein_quant_table) - 1L))
        for (type in c("pca", "rle", "density", "correlation")) {
            prefix <- paste0("p2.", lane_id, ".", stage, ".", type)
            made <- k259Capture(plot <- buildProtNormQcPlotProduct(s4, aesthetics, type)$plot)
            k259Put(paste0(prefix, ".build"), if (identical(made$status, "ok")) "ok" else paste("error:", made$chain))
            if (!identical(made$status, "ok")) { statuses <- c(statuses, "error", "error"); next }
            if (identical(type, "correlation")) k259Facts(paste0("facts.", lane_id, ".", stage), plot)
            statuses <- c(statuses, k259Render(prefix, plot))
        }
        # plotPearson's own (native) plot, the legacy producer.
        prefix <- paste0("p2.", lane_id, ".", stage, ".pearson_native")
        made <- k259Capture(native <- plotPearson(s4, correlation_group = "group", exclude_pool_samples = TRUE))
        k259Put(paste0(prefix, ".build"), if (identical(made$status, "ok")) "ok" else paste("error:", made$chain))
        if (identical(made$status, "ok")) statuses <- c(statuses, k259Render(prefix, native)) else
            statuses <- c(statuses, "error", "error")
        # The application's own correlation figure under warn = 2.
        made <- k259Capture(plot <- buildProtNormQcPlotProduct(s4, aesthetics, "correlation")$plot)
        if (identical(made$status, "ok")) {
            w2 <- withr::with_options(list(warn = 2), k259Capture(ggplot2::ggsave(tempfile(fileext = ".png"), plot,
                width = 10, height = 8, dpi = 150)))
            k259Put(paste0("p2w2.", lane_id, ".", stage, ".correlation.status"), if (identical(w2$status, "ok")) "ok" else
                paste("error:", w2$chain))
            if (!identical(w2$status, "ok")) k259Put(paste0("p2w2.", lane_id, ".", stage, ".correlation.backtrace"),
                w2$backtrace)
        }
    }
    k259Put(paste0("p2.", lane_id, ".ok"), sprintf("%d/%d", sum(statuses == "ok"), length(statuses)))
    k259Write()
}

# ---------------------------------------------------------------------------------------------------------------------
# Part 3: the application started by shinytest2 as the browser lane starts it, with the probe in its process.
# ---------------------------------------------------------------------------------------------------------------------
k259Plots <- file.path(k259Rds, c("prot_lfq-pre_norm_correlation.png.rds", "prot_lfq-post_norm_correlation.png.rds",
    "prot_dia-pre_norm_correlation.png.rds", "prot_dia-post_norm_correlation.png.rds"))
k259Plots <- k259Plots[file.exists(k259Plots)]
k259Put("p3.plots", basename(k259Plots))
k259AppFn <- local({
    pkg <- normalizePath(".")
    out <- k259Out
    plots <- k259Plots
    probe <- normalizePath(file.path("docs", "qa", "workflow-execution", "kit", "kit-259", "k259-probe.R"))
    function() {
        setwd(pkg)
        options(multischolar.test_mode = TRUE)
        pkgload::load_all(export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
        source(probe, local = TRUE)
        k259ProbeRun("app_top", plots, out)
        app <- run_app(test_mode = TRUE)
        inner <- app$serverFuncSource
        app$serverFuncSource <- function() {
            server <- inner()
            function(input, output, session) {
                server(input, output, session)
                done <- FALSE
                shiny::observe({
                    if (done) return(invisible(NULL))
                    done <<- TRUE
                    shiny::withProgress(message = "KIT-259 probe", value = 0.5, k259ProbeRun("app_observer", plots, out))
                    writeLines("done", file.path(out, "p3-done.txt"))
                })
            }
        }
        app
    }
})
if (length(k259Plots)) {
    k259Log("part 3: starting the application through shinytest2")
    k259Put("p3.chrome", tryCatch(chromote::find_chrome(), error = function(e) paste("<error>", conditionMessage(e))))
    driver <- NULL
    started <- k259Capture(driver <- shinytest2::AppDriver$new(app_dir = k259AppFn, name = "kit259", seed = 20260428L,
        timeout = 60000L, load_timeout = 240000L, width = 1440L, height = 1100L,
        options = list(multischolar.test_mode = TRUE, shiny.testmode = TRUE)))
    k259Put("p3.driver", if (identical(started$status, "ok")) "ok" else paste("error:", started$chain))
    if (!is.null(driver)) {
        deadline <- Sys.time() + 240
        while (!file.exists(file.path(k259Out, "p3-done.txt")) && Sys.time() < deadline) Sys.sleep(2)
        k259Put("p3.observer_done", file.exists(file.path(k259Out, "p3-done.txt")))
        logs <- tryCatch(as.data.frame(driver$get_logs()), error = function(e) data.frame(error = conditionMessage(e)))
        logs[] <- lapply(logs, function(x) gsub("[\t\r\n]+", " ", as.character(x)))
        utils::write.table(logs, file.path(k259Out, "p3-applog.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
        tryCatch(driver$stop(), error = function(e) NULL)
    }
}
k259Write()
k259Log("done")
