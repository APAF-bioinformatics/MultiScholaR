# KIT-172: the claim of the ticket (section 3), tested through the declared operations at the application's default
# offset: T1a (unit equivariance of the log step and of the steps after it on the configurations listed, tolerances
# fixed before any run) and T2 (a two-fold change between 0.1 and 0.2 has logFC 1). Fixture F172, seed 172001.
# Excess in t or p alone, with (1) to (3) holding, is a numerical finding (the exceedance rule): it is printed, with
# its size, under the marker KIT172-T1A-TP and does not fail the test.
library(testthat)

kit172Compare <- function(at_s, at_1, shift, tolerance, what) {
  common <- intersect(rownames(at_1), rownames(at_s))
  a <- at_s[common, , drop = FALSE]
  b <- at_1[common, , drop = FALSE]
  expect_identical(is.na(a), is.na(b), info = paste(what, "missing cells"))
  expect_lte(max(abs(a - b - shift), na.rm = TRUE), tolerance, label = what)
  invisible(length(common))
}

kit172Factors <- c(`2^-10` = 2^-10, `1` = 1, `2^10` = 2^10, `0.001` = 0.001, `1000` = 1000)

for (class in c("lipid", "metab")) {
  test_that(paste("KIT-172 T1a:", class, "plain log2 is unit-equivariant through ITSD, between-sample normalisation and DA"), {
    skip_if_not(workflowOwnerCapable())
    values <- kit172Values(172001L)$values
    prefix <- kit172Binding(class)$prefix
    root <- withr::local_tempdir()
    owners <- list()
    withr::defer(for (owner in owners) closeDIAStateAdapter(owner))
    for (itsd in c(TRUE, FALSE)) for (method in c("none", "cyclicloess", "quantile")) {
      config <- paste(class, if (itsd) "itsd" else "noitsd", method, sep = "-")
      reference <- kit172Journey(file.path(root, paste0(config, "-1")), kit172Object(values, class), class, itsd = itsd,
        method = method)
      owners[[length(owners) + 1L]] <- reference$owner
      kept <- rownames(kit172StateMatrix(reference, paste0(prefix, "_intensity_filtered")))
      expect_gte(sum(!startsWith(kept, "IS_")), 500L, label = paste(config, "features after the intensity filter"))
      log_1 <- kit172StateMatrix(reference, paste0(prefix, "_log2"))
      norm_1 <- kit172StateMatrix(reference, paste0(prefix, "_normalized"))
      da_1 <- kit172DaTable(reference)
      for (name in setdiff(names(kit172Factors), "1")) {
        s <- kit172Factors[[name]]
        # The intensity filter's percentile is exact under a power of two only; the decimal factors start after it,
        # on the features it kept at s = 1.
        power_of_two <- name %in% c("2^-10", "2^10")
        input <- if (power_of_two) values * s else values[kept, , drop = FALSE] * s
        journey <- kit172Journey(file.path(root, paste0(config, "-", name)), kit172Object(input, class), class, itsd = itsd,
          method = method, filter = power_of_two)
        owners[[length(owners) + 1L]] <- journey$owner
        info <- paste(config, "s =", name)
        kit172Compare(kit172StateMatrix(journey, paste0(prefix, "_log2")), log_1, log2(s), 1e-12, paste(info, "(1) log2"))
        kit172Compare(kit172StateMatrix(journey, paste0(prefix, "_normalized")), norm_1, log2(s), 1e-9, paste(info, "(2) normalised"))
        da_s <- kit172DaTable(journey)
        key <- function(d) paste(d$assay, d$comparison, d$feature)
        at <- match(key(da_1), key(da_s))
        expect_false(anyNA(at), info = paste(info, "DA rows"))
        expect_lte(max(abs(da_s$logFC[at] - da_1$logFC), na.rm = TRUE), 1e-9, label = paste(info, "(3) logFC"))
        t_excess <- max(abs(da_s$t[at] - da_1$t) / abs(da_1$t), na.rm = TRUE)
        p_1 <- da_1$P.Value
        p_s <- da_s$P.Value[at]
        both_zero <- p_1 == 0 & p_s == 0
        p_excess <- max((abs(log(p_s) - log(p_1)) - (1e-6 * abs(log(p_1)) + 1e-12))[!both_zero], na.rm = TRUE)
        if (t_excess > 1e-6 || p_excess > 0) {
          cat(sprintf("KIT172-T1A-TP\t%s\tmax relative t change %.3g (tolerance 1e-6)\tlargest log-p excess over tolerance %.3g\n",
            info, t_excess, p_excess))
        }
        succeed()
      }
    }
  })
}

test_that("KIT-172 T2: a feature at 0.1 against 0.2 has logFC 1 (ITSD off, no between-sample normalisation, no RUV)", {
  skip_if_not(workflowOwnerCapable())
  values <- kit172Values(172001L)$values
  pattern <- c(0.9, 1.1, 1, 1, 0.95, 1.05)
  values["F001", ] <- c(0.1 * pattern, 0.2 * pattern)
  root <- withr::local_tempdir()
  owners <- list()
  withr::defer(for (owner in owners) closeDIAStateAdapter(owner))
  for (class in c("lipid", "metab")) {
    journey <- kit172Journey(file.path(root, class), kit172Object(values, class), class, itsd = FALSE, method = "none",
      filter = FALSE)
    owners[[length(owners) + 1L]] <- journey$owner
    da <- kit172DaTable(journey)
    expect_lte(abs(da$logFC[da$feature == "F001"] - 1), 1e-12, label = paste(class, "logFC of F001"))
  }
})
