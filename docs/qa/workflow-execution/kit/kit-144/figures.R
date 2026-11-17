# KIT-144: supplementary figures and the per-scenario results table for the RUV validation study (KIT-124).
#
# Reads only the committed scored tables under docs/qa/workflow-execution/kit/kit-124/scored/ (outputs of
# kit-124/aggregate.R) and writes vector PDF, 600 dpi PNG and one CSV to kit-144/out/. No result is typed in here.
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
# (lighter tint of the same hue for a variant), grey for the rest. `short`: the direct label in the error-against-
# power figure; a design with a short label is highlighted there.
designs <- readSpec("
arm | design                                                | short                 | tier                    | colour
1   | No RUV                                                | No RUV                | Reference               | #7F7F7F
2   | True factor known (ceiling)                           | Ceiling               | Reference               | #000000
7   | All features as controls, factors in model            | All features          | Recommended             | #009E73
16  | Technical replicates, factors in model                |                       | Valid with conditions   | #56B4E9
15  | Technical replicates, corrected data tested           |                       | Valid with conditions   | #A6D8F3
8   | Known null controls, factors in model                 | Known nulls           | Valid with conditions   | #0072B2
5   | First-pass nulls (top 25% dropped), factors in model  | First-pass nulls      | Valid with conditions   | #E69F00
3   | Legacy design                                         | Legacy                | Legacy                  | #D55E00
4   | Legacy controls and k, factors in model               |                       | Other RUV designs       | #E9A27A
6   | First-pass nulls (q > 0.5), factors in model          |                       | Other RUV designs       | #F2CD80
9   | Control list with 20% true effects, factors in model  |                       | Other RUV designs       | #9A9A9A
10  | Known null controls, corrected data tested            |                       | Other RUV designs       | #7FB3D9
11  | All features as controls, k = 1, factors in model     |                       | Other RUV designs       | #7FCEB6
R   | RUVr (residuals)                                      |                       | Other RUV designs       | #9A9A9A
Qs  | Pooled QC as replicates, corrected data tested        |                       | Other RUV designs       | #9A9A9A
Qm  | Pooled QC as replicates, factors in model             |                       | Other RUV designs       | #9A9A9A
12  | RUV-4, known null controls                            |                       | Other methods           | #CC79A7
13  | RUVinv, known null controls                           |                       | Other methods           | #CC79A7
14  | Surrogate variable analysis                           |                       | Other methods           | #CC79A7
14E | Surrogate variable analysis, k as for RUV             |                       | Other methods           | #CC79A7
I   | Internal-standard normalisation only                  |                       | Other methods           | #CC79A7
L   | QC-based LOESS drift correction                       |                       | Other methods           | #CC79A7
17  | Known batch in the model                              |                       | Other methods           | #CC79A7
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

# KIT-124 W3 decision rule: the nominal level and the limits of rules 1 and 2 (upper bounds at most 0.075).
NOMINAL <- 0.05
NULL_RATE_LIMIT <- 0.075
FDP_LIMIT <- 0.075
Z <- stats::qnorm(0.975)

SINGLE_COLUMN <- 3.33
DOUBLE_COLUMN <- 7
BASE_PT <- 7
PNG_DPI <- 600

scored_dir <- here::here("docs", "qa", "workflow-execution", "kit", "kit-124", "scored")
out_dir <- here::here("docs", "qa", "workflow-execution", "kit", "kit-144", "out")
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

tickLabel <- function(v) sub("^0\\.", ".", formatC(v, format = "fg", digits = 3))

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
        breaks = c(0, 0.5, 1))
}

# One omic's heatmap: designs in rows, scenarios in columns, the value written in each cell. A dagger over a
# column marks a scenario the scoring did not use to judge designs (`void`).
heatPanel <- function(x, rows, measure, show_designs, show_scenarios) {
    cells <- stats::aggregate(void ~ cell_index + cell + scenario, data = x, FUN = any)
    cells <- cells[order(cells$cell_index), ]
    cells$x <- seq_len(nrow(cells))
    x$x <- cells$x[match(x$cell, cells$cell)]
    x$y <- rows$y[match(x$arm, rows$arm)]
    if (measure == "fdp") {
        x$fill <- log2(x$mean_fdp / NOMINAL)
        x$label <- noLeadingZero(x$mean_fdp, 3)
        x$dark <- x$fill > 2.2 | x$fill < -1.5
    } else {
        x$fill <- x$power
        x$label <- noLeadingZero(x$power, 2)
        x$dark <- !is.na(x$power) & x$power > 0.45
    }
    top <- max(rows$y) + 1
    ggplot(x, aes(x, y)) +
        geom_tile(aes(fill = fill), colour = "white", linewidth = 0.3) +
        geom_text(aes(label = label, colour = I(ifelse(dark, "white", "black"))), size = pt(BASE_PT - 1),
            family = font_family) +
        geom_point(aes(x = 0.2, y = y, colour = I(colour)), data = rows, shape = 15, size = 1.5) +
        geom_text(aes(x = x, y = top - 0.15), data = cells[cells$void, ], label = "†", size = pt(BASE_PT),
            family = font_family) +
        (if (measure == "fdp") fdpFill() else powerFill()) +
        scale_x_continuous(NULL, breaks = cells$x, labels = cells$scenario, expand = expansion(add = c(0.25, 0.05)),
            limits = c(0.05, nrow(cells) + 0.5)) +
        scale_y_continuous(NULL, breaks = rows$y, labels = as.character(rows$design),
            limits = c(min(rows$y) - 0.5, top + 0.2), expand = expansion(0)) +
        labs(title = as.character(x$omic[1])) +
        themeJpr() +
        theme(
            panel.border = element_blank(), axis.ticks = element_blank(),
            axis.text.x = if (show_scenarios) element_text(angle = 40, hjust = 1, vjust = 1) else element_blank(),
            axis.text.y = if (show_designs) element_text() else element_blank(),
            plot.title = element_text(size = BASE_PT + 1, face = "bold", margin = margin(0, 0, 1, 0)),
            legend.position = "bottom", legend.title.position = "top",
            legend.key.width = grid::unit(34, "pt"), legend.key.height = grid::unit(6, "pt")
        )
}

primaryRows <- function(x) x[x$primary, ]

# Heatmaps of every omic in `x`, two per row, on shared design rows; scenario names under each panel when the
# omics do not share their scenarios, else under the last row only.
heatGrid <- function(x, measure, own_scenarios, watermark = NULL) {
    rows <- designRows(x)
    present <- levels(droplevels(x$omic))
    panels <- lapply(seq_along(present), \(i) {
        heatPanel(x[x$omic == present[i], ], rows, measure, show_designs = i %% 2 == 1,
            show_scenarios = own_scenarios || i > length(present) - 2) + watermark
    })
    n_rows <- ceiling(length(present) / 2)
    n_cells <- vapply(present, \(o) length(unique(x$cell[x$omic == o])), integer(1))
    widths <- c(max(n_cells[c(TRUE, FALSE)]), max(n_cells[c(FALSE, TRUE)], 1L))
    list(
        plot = wrap_plots(panels, ncol = 2, guides = "collect", widths = widths) &
            theme(legend.position = "bottom"),
        height = 0.55 + n_rows * (0.3 + 0.15 * nrow(rows)) + (if (own_scenarios) n_rows else 1) * 1.25
    )
}

figCalibrationHeatmap <- function() {
    grid <- heatGrid(primaryRows(core), "fdp", own_scenarios = FALSE)
    savePair(grid$plot, "figS2-calibration-heatmap", DOUBLE_COLUMN, grid$height)
}

figStress <- function() {
    if (!nrow(stress)) {
        cli::cli_warn("no extra-hard scenarios in any scored folder: figure S8 not drawn")
        return(invisible())
    }
    x <- primaryRows(stress)
    watermark <- caption <- NULL
    if (stress_partial) {
        watermark <- annotation_custom(grid::textGrob("PARTIAL", rot = 30, gp = grid::gpar(
            col = grDevices::adjustcolor("black", 0.13), fontsize = 40, fontface = "bold", fontfamily = font_family
        )))
        caption <- sprintf(
            "PARTIAL: extra-hard scenarios scored for %d of %d omics (input: %s).",
            length(unique(x$omic)), nrow(omics), paste(stress_folders, collapse = ", ")
        )
    }
    names <- c(fdp = "figS8a-extra-hard-fdp", power = "figS8b-extra-hard-power")
    for (measure in names(names)) {
        grid <- heatGrid(x, measure, own_scenarios = TRUE, watermark = watermark)
        plot <- grid$plot +
            plot_annotation(caption = caption, theme = themeJpr() + theme(legend.position = "bottom"))
        savePair(plot, names[[measure]], DOUBLE_COLUMN, grid$height + 0.15)
    }
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
        sqrtAxis("x", max(x$null05_upper) * 1.08, c(0.01, 0.05, 0.1, 0.2, 0.3),
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
    centre$short <- designs$short[match(centre$arm, designs$arm)]
    centre$key <- nzchar(centre$short)
    centre$colour <- ifelse(centre$key, designs$colour[match(centre$arm, designs$arm)], "grey70")
    worst <- stats::aggregate(mean_fdp ~ omic + arm, data = x, FUN = max)
    centre$worst_fdp <- worst$mean_fdp[match(paste(centre$omic, centre$arm), paste(worst$omic, worst$arm))]
    centre <- centre[order(centre$key), ]
    key <- centre[centre$key, ]
    # Labels keep clear of every drawn point: unlabelled designs and the worst-scenario markers repel them too.
    # The line from each centre point to its worst-scenario marker counts as three more points.
    along <- lapply(c(1 / 3, 2 / 3, 1), \(f) {
        transform(key, mean_fdp = (sqrt(mean_fdp) + f * (sqrt(worst_fdp) - sqrt(mean_fdp)))^2, short = "")
    })
    obstacles <- do.call(rbind, c(list(key, transform(centre[!centre$key, ], short = "")), along))
    lines <- data.frame(
        omic = factor(omics$omic[1], levels = omics$omic), y = c(NOMINAL, FDP_LIMIT),
        label = c("nominal", "limit")
    )
    plot <- ggplot(centre, aes(power, mean_fdp)) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = FDP_LIMIT, fill = "#009E73", alpha = 0.1) +
        geom_hline(yintercept = NOMINAL, linewidth = 0.3) +
        geom_hline(yintercept = FDP_LIMIT, linewidth = 0.3, linetype = "22") +
        geom_text(aes(x = 0.01, y = y, label = label), data = lines, hjust = 0, vjust = -0.35, size = pt(BASE_PT),
            family = font_family) +
        geom_segment(aes(xend = power, yend = worst_fdp, colour = I(colour)), data = key, linewidth = 0.3) +
        geom_point(aes(y = worst_fdp, colour = I(colour)), data = key, shape = 21, fill = "white", size = 1.3,
            stroke = 0.4) +
        geom_linerange(aes(xmin = power_lower, xmax = power_upper, colour = I(colour)), linewidth = 0.3) +
        geom_linerange(aes(ymin = mean_fdp_lower, ymax = mean_fdp_upper, colour = I(colour)), linewidth = 0.3) +
        geom_point(aes(colour = I(colour), size = I(ifelse(key, 1.5, 0.9)))) +
        ggrepel::geom_text_repel(aes(label = short, colour = I(colour)), data = obstacles, size = pt(BASE_PT),
            family = font_family, segment.size = 0.2, min.segment.length = 0, box.padding = 0.6,
            point.padding = 0.35, force = 8, force_pull = 0.1, max.overlaps = Inf, max.iter = 20000, seed = 7) +
        sqrtAxis("y", max(centre$worst_fdp, centre$mean_fdp_upper) * 1.06, c(0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6),
            "Mean false discovery proportion (square-root axis)") +
        scale_x_continuous("Power in the centre scenario", limits = c(0, 1), breaks = seq(0, 1, 0.25),
            labels = c("0", "0.25", "0.5", "0.75", "1"), expand = expansion(mult = 0.03)) +
        facet_wrap(~omic, ncol = 2) +
        themeJpr()
    savePair(plot, "figS4-error-against-power", DOUBLE_COLUMN, 5.6)
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
    ggplot(x, aes(mean_fdp, y, colour = role, alpha = I(ifelse(void, 0.35, 1)))) +
        geom_vline(xintercept = NOMINAL, linewidth = 0.3) +
        geom_vline(xintercept = FDP_LIMIT, linewidth = 0.3, linetype = "22") +
        geom_line(aes(group = cell), colour = "grey60", linewidth = 0.3) +
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

# ---- Run -------------------------------------------------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeResults()
figCalibrationHeatmap()
figTypeOneForest()
figErrorPower()
figOneChoice()
figStress()
