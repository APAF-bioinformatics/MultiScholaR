# The baseline arm (KIT-159 5.5; section 10.1: T-BASE, T-BASE-filter, T-BASE-batch).

baseDesign <- function(batch = NULL) {
    d <- data.frame(sample = c("A1", "A2", "B1", "B2"), group = c("A", "A", "B", "B"))
    if (!is.null(batch)) d$batch <- batch
    d
}

baseQuant <- function() {
    log2_values <- rbind(f1 = c(10, 10, 11, 11), f2 = c(12, 12, 12, 12), f3 = c(8, 9, 8, 9))
    colnames(log2_values) <- c("A1", "A2", "B1", "B2")
    2^log2_values
}

test_that("T-BASE: median and no normalisation; t equals a direct limma fit; zeros become NA", {
    q <- baseQuant()
    med <- suppressWarnings(valBaselineArm(q, baseDesign(), c(B_vs_A = "B-A"),
                                           normalise = "median"))
    da <- med$da[order(med$da$key), ]
    expect_equal(da$estimate, c(0, -1, -1), tolerance = 1e-12)
    expect_equal(unname(med$sample_medians), c(10, 10, 11, 11), tolerance = 1e-12)
    non <- suppressWarnings(valBaselineArm(q, baseDesign(), c(B_vs_A = "B-A"), normalise = "none"))
    dn <- non$da[order(non$da$key), ]
    expect_equal(dn$estimate, c(1, 0, 0), tolerance = 1e-12)
    group <- factor(c("A", "A", "B", "B"))
    x <- stats::model.matrix(~ 0 + group)
    colnames(x) <- levels(group)
    fit <- limma::lmFit(med$normalised, x)
    fit <- limma::contrasts.fit(fit, limma::makeContrasts(B - A, levels = x))
    fit <- suppressWarnings(limma::eBayes(fit))
    expect_equal(unname(da$t), unname(fit$t[da$key, 1]), tolerance = 1e-12)
    q0 <- q
    q0["f3", "A1"] <- 0
    z <- suppressWarnings(valBaselineArm(q0, baseDesign(), c(B_vs_A = "B-A"), normalise = "none"))
    expect_true(is.na(z$normalised["f3", "A1"]))
})

test_that("T-BASE-filter: a key enters a contrast by that contrast's two groups only", {
    samples <- c("B1", "B2", "C1", "C2", "E1", "E2")
    design <- data.frame(sample = samples, group = rep(c("B", "C", "E"), each = 2))
    withr::local_seed(159011)
    q <- matrix(2^(10 + stats::rnorm(5 * 6)), nrow = 5,
                dimnames = list(paste0("k", 1:5), samples))
    q["k1", c("B1", "B2")] <- NA
    contrasts <- c(E_vs_C = "E-C", E_vs_B = "E-B", C_vs_B = "C-B")
    out <- suppressWarnings(valBaselineArm(q, design, contrasts, normalise = "none"))
    k1 <- out$da[out$da$key == "k1", ]
    expect_identical(k1$contrast, "E_vs_C")
    f <- out$filter
    expect_identical(f$n_keys[f$contrast == "E_vs_C"], 5L)
    expect_identical(f$n_keys_r11[f$contrast == "E_vs_C"], 4L)
    expect_identical(f$n_keys[f$contrast == "E_vs_B"], 4L)
})

test_that("T-BASE-batch: a balanced batch term leaves the fold changes unchanged", {
    q <- baseQuant()
    b <- suppressWarnings(valBaselineArm(q, baseDesign(batch = c(1, 2, 1, 2)), c(B_vs_A = "B-A"),
                                         normalise = "none", batch = "batch"))
    expect_true(all(c("group", "batch") %in% names(b$design_columns)))
    db <- b$da[order(b$da$key), ]
    expect_equal(db$estimate, c(1, 0, 0), tolerance = 1e-12)
})

test_that("the baseline refuses a contrast group below min_per_group", {
    d <- baseDesign()
    d$group <- c("A", "A", "A", "B")
    expect_error(valBaselineArm(baseQuant(), d, c(B_vs_A = "B-A"), normalise = "none"),
                 class = "val_baseline_refused")
})

test_that("the counts baseline runs voom-limma after filterByExpr and TMM", {
    withr::local_seed(159012)
    counts <- matrix(stats::rpois(200 * 6, lambda = 50), nrow = 200,
                     dimnames = list(paste0("g", 1:200), paste0("s", 1:6)))
    design <- data.frame(sample = colnames(counts), group = rep(c("A", "B"), each = 3))
    out <- valBaselineArm(counts, design, c(B_vs_A = "B-A"), data_type = "counts")
    expect_true(nrow(out$da) > 0)
    expect_true(all(c("key", "contrast", "estimate", "p", "s_own", "called_own", "t") %in%
                        names(out$da)))
})
