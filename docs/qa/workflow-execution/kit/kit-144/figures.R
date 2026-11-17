# KIT-144: supplementary figures and tables for the RUV validation study (KIT-124).
#
# Reads only committed files: the scored tables under docs/qa/workflow-execution/kit/kit-124/scored/ (outputs of
# kit-124/aggregate.R); for the figures that need every simulated data set (S6, S7), kit-144/data/ (written by
# kit-144/extract-runs.R from the simulators' run-level outputs); and kit-144/definitions-source.psv, the one table
# of design, scenario and simulator definitions. Writes vector PDF, 600 dpi PNG and CSV to kit-144/out/. No result
# is typed in here.
#
# Usage, from the repository root (one process, seconds):
#   R_LIBS=<simulation libraries> nice -n 19 Rscript docs/qa/workflow-execution/kit/kit-144/figures.R
#
# To add a design (arm): one row in `designs` below. To point at a new scoring: `omics$core_dir`; folders with
# extra-hard scenarios are found by themselves (a folder named "*partial*" is used only when no other has them).

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

# ---- Configuration ---------------------------------------------------------------------------------------------

readSpec <- function(text) {
    utils::read.table(
        text = text, sep = "|", header = TRUE, strip.white = TRUE, quote = "", comment.char = "",
        colClasses = "character", na.strings = character()
    )
}

# One row per design, in the order of every figure. `colour`: Okabe-Ito hue for the designs of the option set
# (lighter tint of the same hue for a variant), grey for the rest. `highlight`: "panel" gives the design
# its own column of panels in the error-against-power figure; "ring" also marks it in every panel (the ceiling).
designs <- readSpec("
arm | design                                                | highlight | tier                    | colour
1   | No RUV                                                | panel     | Reference               | #7F7F7F
2   | True factor known (ceiling)                           | ring      | Reference               | #000000
7   | All features as controls, factors in model            | panel     | Recommended             | #009E73
16  | Technical replicates, factors in model                |           | Valid with conditions   | #56B4E9
15  | Technical replicates, corrected data tested           |           | Valid with conditions   | #A6D8F3
8   | Known null controls, factors in model                 | panel     | Valid with conditions   | #0072B2
5   | First-pass nulls (top 25% dropped), factors in model  | panel     | Valid with conditions   | #E69F00
3   | Legacy design                                         | panel     | Legacy                  | #D55E00
4   | Legacy controls and k, factors in model               |           | Other RUV designs       | #E9A27A
6   | First-pass nulls (q > 0.5), factors in model          |           | Other RUV designs       | #F2CD80
9   | Control list with 20% true effects, factors in model  |           | Other RUV designs       | #9A9A9A
10  | Known null controls, corrected data tested            |           | Other RUV designs       | #7FB3D9
11  | All features as controls, k = 1, factors in model     |           | Other RUV designs       | #7FCEB6
R   | RUVr (residuals)                                      |           | Other RUV designs       | #9A9A9A
Qs  | Pooled QC as replicates, corrected data tested        |           | Other RUV designs       | #9A9A9A
Qm  | Pooled QC as replicates, factors in model             |           | Other RUV designs       | #9A9A9A
12  | RUV-4, known null controls                            |           | Other methods           | #CC79A7
13  | RUVinv, known null controls                           |           | Other methods           | #CC79A7
14  | Surrogate variable analysis                           |           | Other methods           | #CC79A7
14E | Surrogate variable analysis, k as for RUV             |           | Other methods           | #CC79A7
I   | Internal-standard normalisation only                  |           | Other methods           | #CC79A7
L   | QC-based LOESS drift correction                       |           | Other methods           | #CC79A7
17  | Known batch in the model                              |           | Other methods           | #CC79A7
")

omics <- readSpec("
simulator | omic                         | core_dir
S1        | Proteomics, classic          | s1-core
S2        | Proteomics, DIA (limpa)      | s2-core
S3        | Transcriptomics              | s3-core
S4        | Lipidomics and metabolomics  | s4-core
")

# Scenario names as in kit-124/design-review-20261118.md and the KIT-124 cells table. `simulator` "*" = every
# simulator; a row with a simulator overrides it. `mark` picks the scenarios the forest plot singles out.
scenarios <- readSpec("
cell | simulator | scenario                                   | mark
C0   | *         | Centre                                     | centre
C1   | *         | 12 per group                               |
C2   | *         | 10% of features changed                    |
C3   | *         | No unwanted factor                         |
C4   | *         | No true effects (global null)              |
C5   | *         | Factor tied to groups (ρ 0.3)         | tied
C6   | *         | Technical duplicates                       |
C7   | *         | Three groups                               |
S1   | *         | 3 per group                                |
S2   | *         | 40% of features changed                    |
S3   | *         | Factor strongly tied to groups (ρ 0.6) |
S4   | *         | Three factors                              |
S5   | *         | No missing values                          |
S5   | S3        | Low counts                                 |
S6   | *         | 50% missing                                |
S7   | *         | Heavy-tailed noise                         |
S8   | *         | Known unbalanced batch                     |
S9   | *         | 80% of effects in one direction            |
S10  | *         | Loadings tied to abundance                 |
S11  | *         | 3 per group, factor tied to groups         |
S12  | *         | 3 per group, three factors                 |
S13  | *         | No unwanted factor, 3 per group            |
S14  | *         | No unwanted factor, 50% missing            |
S14  | S3        | No unwanted factor, low counts             |
")

# "One choice changed at a time": reference design against the design that differs in the named choice.
# `stand_in`: "simulator:arm" drawn in place of `changed` where the simulator did not run it separately
# (transcriptomics: KIT-124 notes "arm 4 not run (arm 4 = arm 3)").
choice_pairs <- readSpec("
choice                                  | reference | changed | stand_in
Test form (same controls, same k)       | 8         | 10      |
Control set (factors in model in both)  | 7         | 4       | S3:3
")

# The two rules for the number of unwanted factors (kit-124/design-review-20261118.md, section 5), each drawn from
# one design that ran it. The estimated count is one value per data set, shared by every design that uses it.
k_rules <- readSpec("
rule             | arm
Estimated count  | 7
Separation rule  | 3
")

# KIT-124 W3 decision rule: the nominal level and the limits of rules 1 and 2 (upper bounds at most 0.075); rule 3
# counts the data sets with a false discovery proportion above 0.10.
NOMINAL <- 0.05
NULL_RATE_LIMIT <- 0.075
FDP_LIMIT <- 0.075
FDP_TAIL <- 0.10
Z <- stats::qnorm(0.975)

# Figure S6: the quantiles drawn (three nested intervals and the median). Figure S7: k above K_TOP is one row.
QUANTILES <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
K_TOP <- 5

DOUBLE_COLUMN <- 7
HEAT_ROW_IN <- 0.15
HEAT_MAX_IN <- 9.2
BASE_PT <- 7
PNG_DPI <- 600

scored_dir <- here::here("docs", "qa", "workflow-execution", "kit", "kit-124", "scored")
kit_dir <- here::here("docs", "qa", "workflow-execution", "kit", "kit-144")
runs_path <- file.path(kit_dir, "data", "k124-main-runs.tsv.xz")
definitions_path <- file.path(kit_dir, "definitions-source.psv")
out_dir <- file.path(kit_dir, "out")
font_family <- if ("Liberation Sans" %in% systemfonts::system_fonts()$family) "Liberation Sans" else "sans"

# ---- Reading ---------------------------------------------------------------------------------------------------

cell_columns <- c(
    "simulator", "cell", "arm", "subset", "runs", "ok", "refused_designed", "refused_method", "crashed", "scored",
    "mean_k", "null05", "null05_se", "null01", "null01_se", "mean_fdp", "mean_fdp_se", "fdp10_share", "fdp10_lower",
    "fdp10_upper", "anycall_share", "anycall_lower", "anycall_upper", "power", "power_se", "bias", "bias_se", "rmse",
    "rmse_se", "rule1", "rule2", "rule3_limit", "rule3", "rule4", "rule5", "rule6", "rule7", "verdict", "void",
    "cell_role"
)
arm_columns <- c("simulator", "arm", "subset", "primary")

readScored <- function(folder, file, required) {
    path <- file.path(scored_dir, folder, file)
    if (!file.exists(path)) {
        cli::cli_abort("scored table not found: {.path {path}}", class = "kit144_missing_file")
    }
    x <- utils::read.delim(path, quote = "", comment.char = "", stringsAsFactors = FALSE)
    missing_cols <- setdiff(required, names(x))
    if (length(missing_cols)) {
        cli::cli_abort(
            "{.path {path}} lacks column{?s} {.val {missing_cols}}",
            class = "kit144_missing_column"
        )
    }
    x$arm <- as.character(x$arm)
    x[required]
}

# The scored rows of one folder with the `primary` flag of its design table (k124-arms.tsv).
readFolder <- function(folder) {
    cells <- readScored(folder, "k124-cells.tsv", cell_columns)
    arms <- readScored(folder, "k124-arms.tsv", arm_columns)
    x <- merge(cells, arms, by = c("simulator", "arm", "subset"), all.x = TRUE)
    if (anyNA(x$primary)) {
        cli::cli_abort("{.path {folder}}: rows with no design row in k124-arms.tsv", class = "kit144_inconsistent")
    }
    x$folder <- folder
    x
}

stressFolders <- function() {
    folders <- basename(list.dirs(scored_dir, recursive = FALSE))
    has_stress <- vapply(folders, \(f) {
        path <- file.path(scored_dir, f, "k124-cells.tsv")
        file.exists(path) && any(readScored(f, "k124-cells.tsv", "cell_role")$cell_role == "stress")
    }, logical(1))
    folders <- folders[has_stress]
    final <- folders[!grepl("partial", folders)]
    if (length(final)) final else folders
}

# Labels, order and intervals; stops on a design or scenario that the tables above do not name.
labelRows <- function(x, state) {
    unknown <- setdiff(unique(x$arm), designs$arm)
    if (length(unknown)) {
        cli::cli_abort(
            "design{?s} {.val {unknown}} in the scored tables ha{?s/ve} no row in {.code designs}; add one row each",
            class = "kit144_unknown_design"
        )
    }
    if (anyDuplicated(x[c("simulator", "cell", "arm", "subset")])) {
        cli::cli_abort("the same scored row comes from two folders", class = "kit144_inconsistent")
    }
    x$omic <- factor(omics$omic[match(x$simulator, omics$simulator)], levels = omics$omic)
    x$design <- factor(designs$design[match(x$arm, designs$arm)], levels = designs$design)
    own <- match(paste(x$cell, x$simulator), paste(scenarios$cell, scenarios$simulator))
    any_sim <- match(paste(x$cell, "*"), paste(scenarios$cell, scenarios$simulator))
    row <- ifelse(is.na(own), any_sim, own)
    if (anyNA(row) || anyNA(x$omic)) {
        cli::cli_abort("a scored row names a scenario or simulator not listed here", class = "kit144_inconsistent")
    }
    x$scenario <- scenarios$scenario[row]
    x$mark <- scenarios$mark[any_sim]
    x$cell_index <- match(x$cell, unique(scenarios$cell))
    x$state <- state
    x$null05_lower <- pmax(x$null05 - Z * x$null05_se, 0)
    x$null05_upper <- x$null05 + Z * x$null05_se
    x$mean_fdp_lower <- pmax(x$mean_fdp - Z * x$mean_fdp_se, 0)
    x$mean_fdp_upper <- x$mean_fdp + Z * x$mean_fdp_se
    x$power_lower <- pmax(x$power - Z * x$power_se, 0)
    x$power_upper <- pmin(x$power + Z * x$power_se, 1)
    x[order(x$omic, x$cell_index, x$design, x$subset), ]
}

core <- do.call(rbind, lapply(omics$core_dir, readFolder))
core <- labelRows(core[core$cell_role != "stress", ], "final")

stress_folders <- stressFolders()
stress_partial <- any(grepl("partial", stress_folders))
stress <- do.call(rbind, lapply(stress_folders, readFolder))
stress <- if (is.null(stress)) core[0, ] else {
    labelRows(stress[stress$cell_role == "stress", ], if (stress_partial) "partial" else "final")
}

# The scored data sets of the committed per-run table (kit-144/extract-runs.R): one row per omic, main scenario,
# data set and design, on the design's primary subset. A crashed run has no measures and is left out, as in the
# scoring.
run_columns <- c("simulator", "cell", "run", "arm", "subset", "class", "k", "fdp")

readRuns <- function() {
    if (!file.exists(runs_path)) {
        cli::cli_abort("per-run table not found: {.path {runs_path}}", class = "kit144_missing_file")
    }
    x <- utils::read.delim(runs_path, quote = "", comment.char = "", stringsAsFactors = FALSE,
        colClasses = c(arm = "character"))
    missing_cols <- setdiff(run_columns, names(x))
    if (length(missing_cols)) {
        cli::cli_abort(
            "{.path {runs_path}} lacks column{?s} {.val {missing_cols}}",
            class = "kit144_missing_column"
        )
    }
    x[x$class != "crashed" & !is.na(x$fdp), run_columns]
}

# One row per omic, main scenario and design: the quantiles of the false discovery proportion over the data sets
# (q01 .. q99), its mean and the share above FDP_TAIL, with the labels and the "not judged" mark of the scored
# table. Stops when the per-run table does not give the scored table's count, mean, share and mean k.
summariseRuns <- function(runs) {
    keys <- c("simulator", "cell", "arm", "subset")
    x <- do.call(rbind, lapply(split(runs, runs[keys], drop = TRUE), \(r) {
        quantiles <- stats::quantile(r$fdp, QUANTILES, names = FALSE)
        data.frame(
            r[1L, keys], data_sets = nrow(r), t(stats::setNames(quantiles, sprintf("q%02.0f", 100 * QUANTILES))),
            run_mean_fdp = mean(r$fdp), tail_share = mean(r$fdp > FDP_TAIL), run_mean_k = mean(r$k, na.rm = TRUE)
        )
    }))
    scored <- c("omic", "design", "scenario", "cell_index", "void", "scored", "mean_fdp", "fdp10_share", "mean_k")
    merged <- merge(x, core[c(keys, scored)], by = keys)
    same <- function(a, b) {
        (!is.finite(a) & !is.finite(b)) | (is.finite(a) & is.finite(b) & abs(a - b) <= 1e-5 * pmax(abs(b), 1e-6))
    }
    agree <- merged$data_sets == merged$scored & same(merged$run_mean_fdp, merged$mean_fdp) &
        same(merged$tail_share, merged$fdp10_share) & same(merged$run_mean_k, merged$mean_k)
    if (nrow(merged) != nrow(x) || !all(agree)) {
        cli::cli_abort(
            "the per-run table and the scored tables disagree; run {.file extract-runs.R} on the scored inputs",
            class = "kit144_inconsistent"
        )
    }
    merged[order(merged$omic, merged$cell_index, merged$design), ]
}

runs <- readRuns()
run_summary <- summariseRuns(runs)

# The definitions (kit-144/definitions-source.psv): one row per table (design, scenario, simulator), key (arm,
# cell or parameter), simulator and field. "*" stands for every key or every simulator.
if (!file.exists(definitions_path)) {
    cli::cli_abort("definitions not found: {.path {definitions_path}}", class = "kit144_missing_file")
}
definitions <- readSpec(paste(readLines(definitions_path, encoding = "UTF-8"), collapse = "\n"))

# The value of `field` for each key and simulator: the row of that key and simulator, else of the key in every
# simulator, else of every key in that simulator, else of every key in every simulator; NA when there is none.
definitionValue <- function(table, key, simulator, field) {
    rows <- definitions[definitions$table == table & definitions$field == field, ]
    have <- paste(rows$key, rows$simulator)
    at <- cbind(
        match(paste(key, simulator), have), match(paste(key, "*"), have), match(paste("*", simulator), have),
        match(rep("* *", length(key)), have)
    )
    rows$value[at[cbind(seq_along(key), max.col(!is.na(at), ties.method = "first"))]]
}

# ---- Shared plotting pieces ------------------------------------------------------------------------------------

pt <- function(points) points / ggplot2::.pt

themeJpr <- function() {
    theme_bw(base_size = BASE_PT, base_family = font_family) +
        theme(
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "grey30", linewidth = 0.3),
            axis.ticks = element_line(colour = "grey30", linewidth = 0.25),
            axis.text = element_text(colour = "black", size = BASE_PT),
            axis.title = element_text(size = BASE_PT),
            strip.background = element_blank(),
            strip.text = element_text(size = BASE_PT + 1, face = "bold", hjust = 0),
            legend.title = element_text(size = BASE_PT),
            legend.text = element_text(size = BASE_PT),
            legend.key.size = grid::unit(9, "pt"),
            legend.margin = margin(0, 0, 0, 0),
            plot.caption = element_text(size = BASE_PT, hjust = 0),
            plot.margin = margin(3, 4, 3, 3)
        )
}

# Row positions of the designs present in `x`, top to bottom in table order, with a gap between tiers.
designRows <- function(x) {
    rows <- designs[designs$arm %in% x$arm, ]
    tier_index <- match(rows$tier, unique(rows$tier))
    rows$y <- -(seq_len(nrow(rows)) + 0.5 * (tier_index - 1))
    rows$design <- factor(rows$design, levels = designs$design)
    rows
}

noLeadingZero <- function(v, digits) {
    ifelse(is.na(v), "", sub("^0\\.", ".", formatC(v, format = "f", digits = digits)))
}

tickLabel <- function(v) sub("^0\\.", ".", trimws(formatC(v, format = "fg", digits = 3)))

sqrtAxis <- function(axis, limit, breaks, name) {
    scale <- if (axis == "x") scale_x_continuous else scale_y_continuous
    scale(name, transform = "sqrt", breaks = breaks, labels = tickLabel, limits = c(0, limit),
        expand = expansion(mult = c(0, 0.03)))
}

savePair <- function(plot, name, width, height) {
    ggsave(file.path(out_dir, paste0(name, ".pdf")), plot, width = width, height = height,
        device = grDevices::cairo_pdf)
    ggsave(file.path(out_dir, paste0(name, ".png")), plot, width = width, height = height,
        device = ragg::agg_png, dpi = PNG_DPI)
    cli::cli_inform("written: {name}.pdf, {name}.png ({width} x {round(height, 2)} in)")
}

# ---- Heatmaps (S2, S8) -----------------------------------------------------------------------------------------

fdpFill <- function() {
    # Diverging on log2(value / nominal): white at the nominal level, clipped at a quarter and at 16 times it.
    fold <- c(-2, 0, 1, 2.5, 4)
    scale_fill_gradientn(
        "Mean false discovery proportion",
        colours = c("#2166AC", "#F7F7F7", "#FDDBC7", "#D6604D", "#67001F"),
        values = scales::rescale(fold), limits = range(fold), oob = scales::oob_squish_any,
        breaks = -2:4, labels = \(b) tickLabel(NOMINAL * 2^b)
    )
}

powerFill <- function() {
    scale_fill_gradientn("Power", colours = viridisLite::viridis(64, direction = -1), limits = c(0, 1),
        breaks = c(0, 0.5, 1), labels = tickLabel)
}

# One omic's heatmap: the designs that have a value in this omic in rows, scenarios in columns, the value written
# in each cell. A dagger over a column marks a scenario the scoring did not use to judge designs (`void`).
# `n_columns` reserves the same number of columns in every panel of a file, so cells have one size.
heatPanel <- function(x, measure, n_columns, show_scenarios) {
    rows <- designRows(x)
    cells <- stats::aggregate(void ~ cell_index + cell + scenario, data = x, FUN = any)
    cells <- cells[order(cells$cell_index), ]
    cells$x <- seq_len(nrow(cells))
    x$x <- cells$x[match(x$cell, cells$cell)]
    x$y <- rows$y[match(x$arm, rows$arm)]
    if (measure == "fdp") {
        x$fill <- log2(x$value / NOMINAL)
        x$label <- noLeadingZero(x$value, 3)
        x$dark <- x$fill > 2.2 | x$fill < -1.5
    } else {
        x$fill <- x$value
        x$label <- noLeadingZero(x$value, 2)
        x$dark <- x$value > 0.45
    }
    ggplot(x, aes(x, y)) +
        geom_tile(aes(fill = fill), colour = "white", linewidth = 0.3) +
        geom_text(aes(label = label, colour = I(ifelse(dark, "white", "black"))), size = pt(BASE_PT),
            family = font_family) +
        geom_point(aes(x = 0.3, y = y, colour = I(colour)), data = rows, shape = 15, size = 1.6) +
        geom_text(aes(x = x, y = max(rows$y) + 0.9), data = cells[cells$void, ], label = "†",
            size = pt(BASE_PT), family = font_family) +
        (if (measure == "fdp") fdpFill() else powerFill()) +
        scale_x_continuous(NULL, breaks = cells$x, labels = cells$scenario, limits = c(0.15, n_columns + 0.5),
            expand = expansion(0)) +
        scale_y_continuous(NULL, breaks = rows$y, labels = as.character(rows$design),
            limits = heatRowLimits(rows), expand = expansion(0)) +
        labs(title = as.character(x$omic[1])) +
        themeJpr() +
        theme(
            panel.border = element_blank(), axis.ticks = element_blank(),
            axis.text.x = if (show_scenarios) element_text(angle = 40, hjust = 1, vjust = 1) else element_blank(),
            plot.title = element_text(size = BASE_PT + 1, face = "bold", margin = margin(0, 0, 0, 0)),
            legend.position = "bottom", legend.title.position = "top",
            legend.key.width = grid::unit(34, "pt"), legend.key.height = grid::unit(6, "pt")
        )
}

# Room for the designs' rows, half a row below and a line for the daggers above.
heatRowLimits <- function(rows) c(min(rows$y) - 0.5, max(rows$y) + 1.2)

primaryRows <- function(x) x[x$primary, ]

# Heatmap files: the omics of `x` stacked, each with the designs it has, panel heights in proportion to their rows
# and the same cell size everywhere. Omics go to one file while it stays under HEAT_MAX_IN, then to the next
# ("-part<i>"). Scenario names go under the last panel of a file when all omics share their scenarios, else under
# each panel. `banner`: a line set above the figure (the "partial" mark).
heatFiles <- function(x, measure, stem, banner = NULL) {
    x$value <- if (measure == "fdp") x$mean_fdp else x$power
    x <- x[!is.na(x$value), ]
    parts <- lapply(levels(droplevels(x$omic)), \(o) x[x$omic == o, ])
    scenario_sets <- lapply(parts, \(part) unique(part$scenario[order(part$cell_index)]))
    shared <- length(unique(scenario_sets)) == 1
    n_columns <- max(lengths(scenario_sets))
    rows <- vapply(parts, \(part) diff(heatRowLimits(designRows(part))), numeric(1))
    # Height of a block of scenario names set at 40 degrees, from the longest name.
    names_in <- vapply(scenario_sets, \(set) 0.15 + max(nchar(set)) * 0.031, numeric(1))
    panel_in <- rows * HEAT_ROW_IN + 0.3 + if (shared) 0 else names_in
    fixed_in <- 0.6 + (if (shared) names_in[1] else 0) + 0.25 * !is.null(banner)

    file <- integer(length(parts))
    used <- Inf
    for (i in seq_along(parts)) {
        new_file <- used + panel_in[i] > HEAT_MAX_IN - fixed_in
        used <- panel_in[i] + if (new_file) 0 else used
        file[i] <- max(file) + new_file
    }
    for (f in unique(file)) {
        members <- which(file == f)
        panels <- lapply(members, \(i) {
            heatPanel(parts[[i]], measure, n_columns, show_scenarios = !shared || i == max(members))
        })
        banner_theme <- theme(
            legend.position = "bottom",
            plot.title = element_text(size = BASE_PT + 2, face = "bold", colour = "#D55E00", hjust = 0)
        )
        plot <- wrap_plots(panels, ncol = 1, heights = rows[members], guides = "collect") +
            plot_annotation(title = banner, theme = themeJpr() + banner_theme)
        name <- if (max(file) > 1) sprintf("%s-part%d", stem, f) else stem
        savePair(plot, name, DOUBLE_COLUMN, sum(panel_in[members]) + fixed_in)
    }
}

figCalibrationHeatmap <- function() {
    heatFiles(primaryRows(core), "fdp", "figS2-calibration-heatmap")
}

figStress <- function() {
    if (!nrow(stress)) {
        cli::cli_warn("no extra-hard scenarios in any scored folder: figure S8 not drawn")
        return(invisible())
    }
    x <- primaryRows(stress)
    banner <- if (stress_partial) {
        sprintf(
            "PARTIAL: extra-hard scenarios scored for %d of %d omics (input: %s)",
            length(unique(x$omic)), nrow(omics), paste(stress_folders, collapse = ", ")
        )
    }
    heatFiles(x, "fdp", "figS8a-extra-hard-fdp", banner)
    heatFiles(x, "power", "figS8b-extra-hard-power", banner)
}

# ---- Forest plot of the type I error (S3) ----------------------------------------------------------------------

figTypeOneForest <- function() {
    x <- primaryRows(core)
    rows <- designRows(x)
    x$colour <- designs$colour[match(x$arm, designs$arm)]
    n_cells <- length(unique(x$cell))
    x$y <- rows$y[match(x$arm, rows$arm)] - (match(x$cell, unique(x$cell)) - (n_cells + 1) / 2) * 0.8 / n_cells
    kinds <- c("Centre scenario", "Unwanted factor tied to groups", "Other main scenario",
        "Not judged (ceiling miscalibrated)")
    x$kind <- factor(
        ifelse(x$void, kinds[4], ifelse(x$mark == "centre", kinds[1], ifelse(x$mark == "tied", kinds[2], kinds[3]))),
        levels = kinds
    )
    bands <- rows[seq(1, nrow(rows), by = 2), ]
    plot <- ggplot(x, aes(null05, y)) +
        geom_rect(aes(xmin = 0, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5), data = bands, inherit.aes = FALSE,
            fill = "grey94") +
        geom_vline(xintercept = NOMINAL, linewidth = 0.3) +
        geom_vline(xintercept = NULL_RATE_LIMIT, linewidth = 0.3, linetype = "22") +
        geom_linerange(aes(xmin = null05_lower, xmax = null05_upper, colour = I(colour)), linewidth = 0.35) +
        geom_point(aes(shape = kind, colour = I(colour), size = kind), stroke = 0.35, fill = "white") +
        scale_shape_manual(NULL, values = c(16, 17, 21, 4), drop = FALSE) +
        scale_size_manual(NULL, values = c(1.3, 1.3, 0.9, 0.9), drop = FALSE) +
        sqrtAxis("x", max(x$null05_upper) * 1.12, c(0.01, 0.05, 0.1, 0.2, 0.3),
            "Type I error rate: share of true null features with p < 0.05 (square-root axis)") +
        scale_y_continuous(NULL, breaks = rows$y, labels = as.character(rows$design),
            limits = c(min(rows$y) - 0.5, max(rows$y) + 0.5), expand = expansion(0)) +
        facet_grid(. ~ omic, labeller = label_wrap_gen(18)) +
        themeJpr() +
        theme(axis.ticks.y = element_blank(), legend.position = "bottom", legend.justification = "right",
            panel.spacing.x = grid::unit(7, "pt"), strip.text = element_text(vjust = 0)) +
        guides(shape = guide_legend(nrow = 2, override.aes = list(colour = "black")))
    savePair(plot, "figS3-type-one-error-forest", DOUBLE_COLUMN, 1 + 0.27 * nrow(rows))
}

# ---- Error against power (S4) ----------------------------------------------------------------------------------

figErrorPower <- function() {
    x <- primaryRows(core)
    x <- x[!x$void & !is.na(x$power), ]
    centre <- x[x$mark == "centre", ]
    # Each design's judged main scenario with the highest mean false discovery proportion, set at the power of
    # the centre scenario.
    worst <- x[order(-x$mean_fdp), ]
    worst <- worst[!duplicated(worst[c("omic", "arm")]), ]
    at <- match(paste(worst$omic, worst$arm), paste(centre$omic, centre$arm))
    worst[c("power", "power_lower", "power_upper")] <- centre[at, c("power", "power_lower", "power_upper")]
    worst <- worst[!is.na(at), ]

    # One column of panels per highlighted design: that design in its colour, the ceiling as a ring, the rest grey.
    shown <- designs[nzchar(designs$highlight), ]
    roles <- c(
        "Design named above the panel, centre scenario",
        "The same design in the main scenario where its false discovery proportion is highest (joined)",
        "True factor known (ceiling), centre scenario",
        "All other designs and methods, centre scenario"
    )
    ring <- centre[centre$arm %in% shown$arm[shown$highlight == "ring"], ]
    columns <- do.call(rbind, lapply(seq_len(nrow(shown)), \(i) {
        own <- centre$arm == shown$arm[i]
        column <- rbind(
            transform(centre[own, ], role = roles[1], colour = shown$colour[i]),
            transform(worst[worst$arm == shown$arm[i], ], role = roles[2], colour = shown$colour[i]),
            transform(ring, role = roles[3], colour = "black"),
            transform(centre[!own, ], role = roles[4], colour = "grey70")
        )
        column$panel <- factor(shown$design[i], levels = shown$design)
        column
    }))
    columns$role <- factor(columns$role, levels = roles)
    # Drawn in this order: the other designs, the ring, the highest value, the design itself on top.
    columns <- columns[order(match(columns$role, roles[c(4, 3, 2, 1)])), ]
    joins <- merge(
        columns[columns$role == roles[1], c("omic", "panel", "colour", "power", "mean_fdp")],
        columns[columns$role == roles[2], c("omic", "panel", "mean_fdp")],
        by = c("omic", "panel"), suffixes = c("", "_worst")
    )
    lines <- data.frame(
        y = c(NOMINAL, FDP_LIMIT),
        line = c(sprintf("Nominal level, %s", NOMINAL), sprintf("Limit of the decision rule, %s", FDP_LIMIT))
    )
    lines$line <- factor(lines$line, levels = lines$line)

    plot <- ggplot(columns, aes(power, mean_fdp)) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = FDP_LIMIT, fill = "#009E73", alpha = 0.1) +
        geom_hline(aes(yintercept = y, linetype = line), data = lines, linewidth = 0.3) +
        geom_segment(aes(xend = power, yend = mean_fdp_worst, colour = I(colour)), data = joins,
            linewidth = 0.35) +
        geom_linerange(aes(xmin = power_lower, xmax = power_upper, colour = I(colour)),
            data = columns[columns$role %in% roles[c(1, 4)], ], linewidth = 0.3) +
        geom_linerange(aes(ymin = mean_fdp_lower, ymax = mean_fdp_upper, colour = I(colour)),
            data = columns[columns$role != roles[3], ], linewidth = 0.3) +
        geom_point(aes(shape = role, size = role, colour = I(colour)), stroke = 0.45, fill = "white") +
        scale_shape_manual(NULL, values = stats::setNames(c(16, 21, 1, 16), roles)) +
        scale_size_manual(NULL, values = stats::setNames(c(1.5, 1.2, 3.2, 0.7), roles)) +
        scale_linetype_manual(NULL, values = stats::setNames(c("solid", "22"), levels(lines$line))) +
        sqrtAxis("y", max(columns$mean_fdp_upper) * 1.06, c(0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6),
            "Mean false discovery proportion (square-root axis)") +
        scale_x_continuous("Power in the centre scenario", limits = c(0, 1), breaks = c(0, 0.5, 1),
            labels = tickLabel, expand = expansion(mult = 0.06)) +
        facet_grid(omic ~ panel, labeller = labeller(panel = label_wrap_gen(17), omic = label_wrap_gen(16))) +
        guides(
            shape = guide_legend(order = 1, ncol = 1,
                override.aes = list(colour = c("grey25", "grey25", "black", "grey70"))),
            size = guide_legend(order = 1, ncol = 1),
            linetype = guide_legend(order = 2, ncol = 1, theme = theme(legend.key.width = grid::unit(18, "pt")))
        ) +
        themeJpr() +
        theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "top",
            legend.spacing.x = grid::unit(14, "pt"), legend.key.height = grid::unit(10, "pt"),
            legend.key.spacing.y = grid::unit(0, "pt"),
            strip.text.x = element_text(size = BASE_PT, vjust = 0), strip.text.y = element_text(size = BASE_PT),
            panel.spacing = grid::unit(5, "pt"))
    savePair(plot, "figS4-error-against-power", DOUBLE_COLUMN, 6.6)
}

# ---- One choice changed at a time (S5) -------------------------------------------------------------------------

choicePanel <- function(pair, show_scenarios) {
    x <- primaryRows(core)
    stand_in <- strsplit(pair$stand_in, ":", fixed = TRUE)[[1]]
    changed_arm <- ifelse(length(stand_in) == 2 & x$simulator == stand_in[1], stand_in[2], pair$changed)
    x <- x[x$arm == pair$reference | x$arm == changed_arm, ]
    roles <- designs$design[match(c(pair$reference, pair$changed), designs$arm)]
    x$role <- factor(ifelse(x$arm == pair$reference, roles[1], roles[2]), levels = roles)
    cells <- unique(core[order(core$cell_index), c("cell", "scenario")])
    x$y <- -match(x$cell, cells$cell) + ifelse(x$arm == pair$reference, 0.16, -0.16)
    colours <- stats::setNames(designs$colour[match(c(pair$reference, pair$changed), designs$arm)], roles)
    # A line joins the two designs of a scenario where the omic has both.
    link_columns <- c("omic", "cell", "void", "mean_fdp", "y")
    links <- merge(x[x$arm == pair$reference, link_columns], x[x$arm != pair$reference, link_columns],
        by = c("omic", "cell", "void"), suffixes = c("", "_changed"))
    ggplot(x, aes(mean_fdp, y, colour = role, alpha = I(ifelse(void, 0.35, 1)))) +
        geom_vline(xintercept = NOMINAL, linewidth = 0.3) +
        geom_vline(xintercept = FDP_LIMIT, linewidth = 0.3, linetype = "22") +
        geom_segment(aes(mean_fdp, y, xend = mean_fdp_changed, yend = y_changed, alpha = I(ifelse(void, 0.35, 1))),
            data = links, inherit.aes = FALSE, colour = "grey60", linewidth = 0.3) +
        geom_linerange(aes(xmin = mean_fdp_lower, xmax = mean_fdp_upper), linewidth = 0.4) +
        geom_point(aes(shape = role), size = 1.3, stroke = 0.5, fill = "white") +
        scale_colour_manual(pair$choice, values = colours) +
        scale_shape_manual(pair$choice, values = stats::setNames(c(16, 21), roles)) +
        sqrtAxis("x", max(primaryRows(core)$mean_fdp_upper[primaryRows(core)$arm %in% unlist(choice_pairs[2:3])]) *
            1.04, c(0.01, 0.05, 0.2, 0.4, 0.6),
            if (show_scenarios) "Mean false discovery proportion (square-root axis)") +
        scale_y_continuous(NULL, breaks = -seq_len(nrow(cells)), labels = cells$scenario,
            limits = c(-nrow(cells) - 0.5, -0.5), expand = expansion(0)) +
        facet_grid(. ~ omic, drop = FALSE, labeller = label_wrap_gen(18)) +
        themeJpr() +
        theme(
            legend.position = "top", legend.justification = "left", legend.title = element_text(face = "bold"),
            legend.box.margin = margin(0, 0, -6, 0), panel.spacing.x = grid::unit(7, "pt"),
            panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
            axis.ticks.y = element_blank(),
            strip.text = if (show_scenarios) element_blank() else element_text(vjust = 0)
        )
}

figOneChoice <- function() {
    panels <- lapply(seq_len(nrow(choice_pairs)), \(i) choicePanel(choice_pairs[i, ], show_scenarios = i > 1))
    savePair(wrap_plots(panels, ncol = 1), "figS5-one-choice-changed", DOUBLE_COLUMN, 4.6)
}

# ---- Per-scenario results table --------------------------------------------------------------------------------

writeResults <- function() {
    x <- rbind(core, stress)
    x$scenario_class <- ifelse(x$cell_role == "stress", "extra-hard", "main")
    x$tier <- designs$tier[match(x$arm, designs$arm)]
    x$not_judged <- x$void
    columns <- c(
        "omic", "simulator", "scenario_class", "cell", "scenario", "design", "tier", "arm", "subset", "primary",
        "state", "not_judged", "runs", "scored", "refused_designed", "refused_method", "crashed", "mean_k",
        "null05", "null05_se", "null05_lower", "null05_upper", "null01", "null01_se", "mean_fdp", "mean_fdp_se",
        "mean_fdp_lower", "mean_fdp_upper", "fdp10_share", "fdp10_lower", "fdp10_upper", "anycall_share",
        "anycall_lower", "anycall_upper", "power", "power_se", "power_lower", "power_upper", "bias", "bias_se",
        "rmse", "rmse_se", "rule1", "rule2", "rule3_limit", "rule3", "rule4", "rule5", "rule6", "rule7", "verdict",
        "folder"
    )
    path <- file.path(out_dir, "tableS-per-scenario-results.csv")
    utils::write.csv(x[columns], path, row.names = FALSE, na = "")
    cli::cli_inform("written: {basename(path)} ({nrow(x)} rows)")
}

# ---- Distribution of the false discovery proportion over the data sets (S6) ------------------------------------

# One design's row of panels (one per omic). For each main scenario: three nested percentile intervals of the
# false discovery proportion over the data sets, the median and the mean; in the column at the right of each
# panel, the share of data sets above FDP_TAIL. `strips`: omic names above the panels; `axis`: the axis below.
distributionPanel <- function(x, design, cells, strips, axis) {
    x$y <- -match(x$cell, cells$cell)
    x$alpha <- ifelse(x$void, 0.35, 1)
    intervals <- data.frame(
        lower = c("q25", "q05", "q01"), upper = c("q75", "q95", "q99"),
        interval = c("25th to 75th percentile", "5th to 95th", "1st to 99th")
    )
    bands <- do.call(rbind, lapply(seq_len(nrow(intervals)), \(i) {
        data.frame(x[c("omic", "y", "alpha")], lower = x[[intervals$lower[i]]], upper = x[[intervals$upper[i]]],
            interval = intervals$interval[i])
    }))
    bands$interval <- factor(bands$interval, levels = intervals$interval)
    marks <- rbind(
        data.frame(x[c("omic", "y", "alpha")], value = x$q50, mark = "Median"),
        data.frame(x[c("omic", "y", "alpha")], value = x$run_mean_fdp, mark = "Mean")
    )
    marks$mark <- factor(marks$mark, levels = c("Median", "Mean"))
    lines <- data.frame(
        x = c(NOMINAL, FDP_TAIL),
        line = c(
            sprintf("Nominal level, %s", NOMINAL),
            sprintf("%.2f: the share of data sets above it is the number at the right", FDP_TAIL)
        )
    )
    lines$line <- factor(lines$line, levels = lines$line)
    # The column of shares lies beyond 1 on the axis: text ends at `share_x`, the axis at `axis_end`. Its tick label
    # is set one line below the others, which it would touch.
    share_x <- 1.5
    axis_end <- 1.58
    ggplot(x, aes(y = y)) +
        geom_vline(aes(xintercept = x, linetype = line), data = lines, linewidth = 0.3) +
        geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey30") +
        geom_linerange(aes(xmin = lower, xmax = upper, linewidth = interval, alpha = I(alpha)), data = bands,
            colour = design$colour, key_glyph = "path") +
        geom_point(aes(value, shape = mark, size = mark, alpha = I(alpha)), data = marks, colour = "black",
            fill = "white", stroke = 0.4) +
        geom_text(aes(share_x, label = noLeadingZero(tail_share, 2), alpha = I(alpha)), hjust = 1,
            size = pt(BASE_PT), family = font_family) +
        scale_linewidth_manual(NULL, values = stats::setNames(c(2.3, 1, 0.3), intervals$interval)) +
        scale_shape_manual(NULL, values = c(Median = 124, Mean = 23)) +
        scale_size_manual(NULL, values = c(Median = 2.1, Mean = 1.2)) +
        scale_linetype_manual(NULL, values = stats::setNames(c("solid", "13"), levels(lines$line))) +
        scale_x_continuous(
            if (axis) "False discovery proportion of a simulated data set (square-root axis)",
            transform = "sqrt", breaks = c(0, 0.05, 0.2, 0.5, 1, 1.27),
            labels = c(tickLabel(c(0, 0.05, 0.2, 0.5, 1)), paste0("\n>", noLeadingZero(FDP_TAIL, 2))),
            limits = c(0, axis_end), expand = expansion(mult = c(0.03, 0))
        ) +
        scale_y_continuous(NULL, breaks = -seq_len(nrow(cells)), labels = cells$scenario,
            limits = c(-nrow(cells) - 0.5, -0.5), expand = expansion(0)) +
        facet_grid(. ~ omic, drop = FALSE, labeller = label_wrap_gen(22)) +
        labs(title = design$design) +
        guides(
            linewidth = guide_legend(order = 1, ncol = 1, override.aes = list(colour = "grey25"),
                theme = theme(legend.key.width = grid::unit(16, "pt"))),
            shape = guide_legend(order = 2, ncol = 1), size = guide_legend(order = 2, ncol = 1),
            linetype = guide_legend(order = 3, ncol = 1, theme = theme(legend.key.width = grid::unit(16, "pt")))
        ) +
        themeJpr() +
        theme(
            plot.title = element_text(size = BASE_PT + 1, face = "bold", margin = margin(0, 0, 1, 0)),
            panel.spacing.x = grid::unit(7, "pt"), axis.ticks.y = element_blank(),
            panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
            strip.text = if (strips) element_text(size = BASE_PT, face = "plain", vjust = 0) else element_blank(),
            axis.text.x = if (axis) element_text() else element_blank(),
            axis.ticks.x = if (axis) element_line() else element_blank(),
            legend.key.height = grid::unit(10, "pt"), legend.key.spacing.y = grid::unit(0, "pt")
        )
}

figFdpDistribution <- function() {
    shown <- designs[nzchar(designs$highlight), ]
    x <- run_summary[run_summary$arm %in% shown$arm, ]
    cells <- unique(core[order(core$cell_index), c("cell", "scenario")])
    panels <- lapply(seq_len(nrow(shown)), \(i) {
        distributionPanel(x[x$arm == shown$arm[i], ], shown[i, ], cells, strips = i == 1, axis = i == nrow(shown))
    })
    legend_theme <- theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "top",
        legend.spacing.x = grid::unit(14, "pt"))
    plot <- wrap_plots(panels, ncol = 1, guides = "collect") +
        plot_annotation(theme = themeJpr() + legend_theme)
    savePair(plot, "figS6-fdp-distribution", DOUBLE_COLUMN, 1.2 + nrow(shown) * (0.22 + 0.115 * nrow(cells)))
}

# ---- Number of unwanted factors: estimated against true (S7) ---------------------------------------------------

shareLabel <- function(v) ifelse(v == 0, "", ifelse(v < 0.005, "<.01", noLeadingZero(v, 2)))

# The share of the data sets of each main scenario at each number of factors k, for one omic and each rule of
# `k_rules`, with the true number from the definitions. k above K_TOP is counted in the row of K_TOP.
factorCounts <- function(simulator) {
    x <- merge(runs[runs$simulator == simulator & !is.na(runs$k), ], k_rules, by = "arm")
    cells <- core[core$simulator == simulator & core$primary & core$arm %in% k_rules$arm, ]
    cells <- stats::aggregate(void ~ cell_index + cell + scenario + omic, data = cells, FUN = any)
    cells <- cells[order(cells$cell_index), ]
    cells$x <- seq_len(nrow(cells))
    cells$true_k <- as.numeric(definitionValue("scenario", cells$cell, simulator, "unwanted_factors"))
    if (anyNA(cells$true_k)) {
        cli::cli_abort("no number of unwanted factors defined for a scenario of {simulator}",
            class = "kit144_inconsistent")
    }
    levels_k <- sort(unique(c(pmin(x$k, K_TOP), cells$true_k)))
    counts <- table(
        rule = factor(x$rule, levels = k_rules$rule), cell = factor(x$cell, levels = cells$cell),
        k = factor(pmin(x$k, K_TOP), levels = levels_k)
    )
    shares <- as.data.frame(prop.table(counts, c(1, 2)), responseName = "share", stringsAsFactors = FALSE)
    shares <- shares[is.finite(shares$share), ]
    shares$rule <- factor(shares$rule, levels = k_rules$rule)
    shares$x <- cells$x[match(shares$cell, cells$cell)]
    shares$y <- match(shares$k, levels_k)
    labels_k <- tickLabel(levels_k)
    labels_k[levels_k == K_TOP & any(x$k > K_TOP)] <- paste0("≥", K_TOP)
    truth <- merge(cells, data.frame(rule = factor(k_rules$rule, levels = k_rules$rule)))
    truth$y <- match(truth$true_k, levels_k)
    list(shares = shares, cells = cells, truth = truth, labels_k = labels_k)
}

# One omic's panels (one per rule): scenarios in columns, k in rows, the share written in each cell that has data
# sets, the true number outlined. A dagger over a column marks a scenario not judged (as in the other figures).
countPanel <- function(counts, strips, show_scenarios) {
    cells <- counts$cells
    n_k <- length(counts$labels_k)
    truth_name <- "True number of unwanted factors"
    rule_names <- stats::setNames(
        sprintf("%s (as run in: %s)", k_rules$rule, designs$design[match(k_rules$arm, designs$arm)]), k_rules$rule
    )
    ggplot(counts$shares, aes(x, y)) +
        geom_tile(aes(fill = share), colour = "white", linewidth = 0.3) +
        geom_text(aes(label = shareLabel(share), colour = I(ifelse(share > 0.55, "white", "black"))),
            size = pt(BASE_PT), family = font_family) +
        geom_tile(aes(linetype = truth_name), data = counts$truth, fill = NA, colour = "black", linewidth = 0.45,
            width = 0.94, height = 0.9) +
        geom_text(aes(x = x, y = n_k + 0.9), data = cells[cells$void, ], label = "†", size = pt(BASE_PT),
            family = font_family) +
        scale_fill_gradientn("Share of the simulated data sets",
            colours = c("#F3F3F3", "#C6DBEF", "#4292C6", "#08306B"), limits = c(0, 1), breaks = c(0, 0.5, 1),
            labels = tickLabel) +
        scale_linetype_manual(NULL, values = stats::setNames("solid", truth_name)) +
        scale_x_continuous(NULL, breaks = seq_len(nrow(cells)), labels = cells$scenario,
            limits = c(0.5, nrow(cells) + 0.5), expand = expansion(0)) +
        scale_y_continuous("Factors removed, k", breaks = seq_len(n_k), labels = counts$labels_k,
            limits = c(0.5, n_k + 1.2), expand = expansion(0)) +
        facet_grid(. ~ rule, drop = FALSE, labeller = labeller(rule = as_labeller(rule_names))) +
        labs(title = as.character(cells$omic[1])) +
        guides(fill = guide_colourbar(order = 1), linetype = guide_legend(order = 2,
            theme = theme(legend.key.width = grid::unit(14, "pt"), legend.key.height = grid::unit(10, "pt")))) +
        themeJpr() +
        theme(
            panel.border = element_blank(), axis.ticks = element_blank(), panel.spacing.x = grid::unit(10, "pt"),
            axis.text.x = if (show_scenarios) element_text(angle = 40, hjust = 1, vjust = 1) else element_blank(),
            plot.title = element_text(size = BASE_PT + 1, face = "bold", margin = margin(0, 0, 0, 0)),
            strip.text = if (strips) element_text(size = BASE_PT, face = "plain", vjust = 0) else element_blank(),
            legend.position = "bottom", legend.title.position = "top",
            legend.key.width = grid::unit(34, "pt"), legend.key.height = grid::unit(6, "pt")
        )
}

figFactorCount <- function() {
    counts <- lapply(omics$simulator, factorCounts)
    rows <- vapply(counts, \(part) length(part$labels_k) + 1.2, numeric(1))
    panels <- lapply(seq_along(counts), \(i) {
        countPanel(counts[[i]], strips = i == 1, show_scenarios = i == length(counts))
    })
    names_in <- 0.15 + max(nchar(core$scenario)) * 0.031
    plot <- wrap_plots(panels, ncol = 1, heights = rows, guides = "collect") +
        plot_annotation(theme = themeJpr() + theme(legend.position = "bottom", legend.box.just = "bottom"))
    savePair(plot, "figS7-factors-estimated-against-true", DOUBLE_COLUMN,
        sum(rows * HEAT_ROW_IN + 0.3) + names_in + 0.85)
}

# ---- Definition tables -----------------------------------------------------------------------------------------

# The versions a simulator recorded, from the sidecar of its scoring.
softwareLine <- function(folder) {
    lines <- readLines(file.path(scored_dir, folder, "k124-aggregate.meta.txt"))
    from <- grep("^simulator S[1-4] meta", lines)
    if (length(from) != 1L) {
        cli::cli_abort("{.path {folder}}: no simulator block in the scoring's sidecar", class = "kit144_inconsistent")
    }
    block <- lines[seq(from + 1L, length(lines))]
    paste(trimws(grep("^  (R version|R:|packages:|DPC engine:)", block, value = TRUE)), collapse = "; ")
}

# Three tables from the definitions: designs and scenarios, one row per omic, with one column per field; the
# simulators' parameters, one row per omic and parameter. A design or scenario without a definition stops the run.
writeDefinitions <- function() {
    scored <- rbind(core, stress)
    write <- function(x, name) {
        utils::write.csv(x, file.path(out_dir, name), row.names = FALSE, na = "")
        cli::cli_inform("written: {name} ({nrow(x)} rows)")
    }
    complete <- function(x, table, key, fields) {
        for (field in fields) {
            x[[field]] <- definitionValue(table, x[[key]], x$simulator, field)
        }
        if (anyNA(x[fields])) {
            cli::cli_abort("a {table} of the scored tables has no definition in {.path {definitions_path}}",
                class = "kit144_inconsistent")
        }
        x
    }

    design_rows <- unique(scored[c("omic", "simulator", "arm", "design")])
    design_rows$tier <- designs$tier[match(design_rows$arm, designs$arm)]
    design_rows <- complete(design_rows, "design", "arm",
        c("controls", "replicate_sets", "number_of_factors", "test_on"))
    design_rows$scenarios_scored <- vapply(seq_len(nrow(design_rows)), \(i) {
        own <- scored$simulator == design_rows$simulator[i] & scored$arm == design_rows$arm[i]
        paste(unique(scored$cell[own]), collapse = " ")
    }, character(1))
    write(design_rows[order(design_rows$omic, design_rows$design), ], "tableS-design-definitions.csv")

    scenario_rows <- merge(omics[c("simulator", "omic")], data.frame(cell = unique(scenarios$cell)))
    scenario_rows$cell_index <- match(scenario_rows$cell, unique(scenarios$cell))
    own <- match(paste(scenario_rows$cell, scenario_rows$simulator), paste(scenarios$cell, scenarios$simulator))
    any_sim <- match(paste(scenario_rows$cell, "*"), paste(scenarios$cell, scenarios$simulator))
    scenario_rows$scenario <- scenarios$scenario[ifelse(is.na(own), any_sim, own)]
    scenario_rows$scenario_class <- ifelse(startsWith(scenario_rows$cell, "C"), "main", "extra-hard")
    scenario_rows <- complete(scenario_rows, "scenario", "cell", c(
        "groups", "units_per_group", "technical_replicates", "effects_share", "effects_up_share",
        "unwanted_factors", "factor_group_correlation", "missing_share", "noise"
    ))
    for (field in c("note", "not_run")) {
        scenario_rows[[field]] <- definitionValue("scenario", scenario_rows$cell, scenario_rows$simulator, field)
    }
    scenario_rows <- scenario_rows[order(match(scenario_rows$simulator, omics$simulator), scenario_rows$cell_index), ]
    write(scenario_rows[names(scenario_rows) != "cell_index"], "tableS-scenario-definitions.csv")

    parameters <- definitions[definitions$table == "simulator", ]
    parameter_rows <- do.call(rbind, lapply(omics$simulator, \(sim) {
        mine <- parameters[parameters$simulator %in% c(sim, "*"), ]
        data.frame(simulator = sim, parameter = mine$key, value = mine$value, source = mine$source)
    }))
    software <- data.frame(
        simulator = omics$simulator, parameter = "software", value = vapply(omics$core_dir, softwareLine, ""),
        source = sprintf("kit-124/scored/%s/k124-aggregate.meta.txt", omics$core_dir)
    )
    parameter_rows <- rbind(parameter_rows, software)
    parameter_rows$omic <- omics$omic[match(parameter_rows$simulator, omics$simulator)]
    parameter_rows <- parameter_rows[order(match(parameter_rows$simulator, omics$simulator)), ]
    write(parameter_rows[c("omic", "simulator", "parameter", "value", "source")],
        "tableS-simulator-parameters.csv")
}

# ---- Run -------------------------------------------------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Outputs of an earlier run under other file names (a figure split differently) do not stay behind.
unlink(list.files(out_dir, pattern = "^(fig|table)S.*[.](pdf|png|csv)$", full.names = TRUE))
writeResults()
figCalibrationHeatmap()
figTypeOneForest()
figErrorPower()
figOneChoice()
figStress()
figFdpDistribution()
figFactorCount()
writeDefinitions()
