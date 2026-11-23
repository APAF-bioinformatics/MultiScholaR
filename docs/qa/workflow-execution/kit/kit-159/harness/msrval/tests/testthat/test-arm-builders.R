# The workflow arm builders (KIT-159 7.5) and their inputs per route (valArmInputs). The values
# are pinned from the product's normalisation panels (R/utils_omic_gui_config.R, R/mod_omic_norm_panels.R)
# and are re-read by G3 in the lifecycle.

test_that("shipped: no option and no override on any route", {
    for (route in driverRoutes()) {
        x <- valArmInputs(valArmShipped(), route)
        expect_length(x$options, 0L)
        expect_length(x$overrides, 0L)
    }
    expect_identical(valArmShipped()$kind, "workflow")
    expect_identical(valArmShipped()$id, "shipped")
})

test_that("legacy: the gate off and the shipped mode and parameters written explicitly", {
    dia <- valArmInputs(valArmRuvLegacy(), "prot_dia")
    expect_identical(dia$options, list(MultiScholaR.ruv_design_panel = FALSE))
    expect_identical(dia$overrides, list(`normalization-ruv_mode` = "automatic"))
    limpa <- valArmInputs(valArmRuvLegacy(), "prot_dia_limpa")
    expect_identical(limpa$overrides, list(`normalization-ruv_mode` = "covariates"))
    lipid <- valArmInputs(valArmRuvLegacy(), "lipid_canonical")
    expect_identical(lipid$overrides, list(`norm-ruv_mode` = "automatic", `norm-ruv_percentage` = 5,
                                           `norm-ruv_k` = 2))
    trx <- valArmInputs(valArmRuvLegacy(), "trx_screen")
    expect_identical(trx$overrides, list(`norm-ruv_mode` = "automatic"))
})

test_that("calibrated: the gate on; the design the server selects, recorded as expected", {
    x <- valArmInputs(valArmRuvCalibrated(), "metab_lc")
    expect_identical(x$options, list(MultiScholaR.ruv_design_panel = TRUE))
    expect_length(x$overrides, 0L)
    expect_identical(x$expect, list(`norm-ruv_design` = "groups_all_features"))
    tech <- valArmInputs(valArmRuvCalibrated(), "trx_screen", technical_replicates = TRUE)
    expect_identical(tech$expect, list(`norm-ruv_design` = "technical_replicates"))
    prot <- valArmInputs(valArmRuvCalibrated(), "prot_tmt")
    expect_identical(prot$expect, list(`normalization-ruv_design` = "groups_all_features"))
})

test_that("off: skip with the gate off; transcriptomics design none with the gate on", {
    dia <- valArmInputs(valArmRuvOff(), "prot_lfq")
    expect_identical(dia$options, list(MultiScholaR.ruv_design_panel = FALSE))
    expect_identical(dia$overrides, list(`normalization-ruv_mode` = "skip"))
    trx <- valArmInputs(valArmRuvOff(), "trx_screen")
    expect_identical(trx$options, list(MultiScholaR.ruv_design_panel = TRUE))
    expect_identical(trx$overrides, list(`norm-ruv_design` = "none"))
})

test_that("design and supplied controls: the gate on, the design and k rule; the list one per line", {
    x <- valArmInputs(valArmRuvDesign("first_pass_nulls", k_rule = "sweep_bounded"), "lipid_canonical")
    expect_identical(x$options, list(MultiScholaR.ruv_design_panel = TRUE))
    expect_identical(x$overrides, list(`norm-ruv_design` = "first_pass_nulls",
                                       `norm-ruv_k_rule` = "sweep_bounded"))
    expect_error(valArmRuvDesign("made_up"), class = "val_spec_invalid")
    s <- valArmRuvSupplied("controls")
    expect_identical(s$controls_role, "controls")
    y <- valArmInputs(s, "prot_dia", roles = list(controls = c("P1", "P2", "P3")))
    expect_identical(y$overrides, list(`normalization-ruv_design` = "supplied_controls",
                                       `normalization-ruv_control_list` = "P1\nP2\nP3"))
    expect_identical(y$mask, "normalization-ruv_control_list")
    expect_error(valArmInputs(s, "prot_dia", roles = list()), class = "val_spec_invalid")
})

test_that("the matrix stage of a cell follows its route and arm (7.3; A12, revision 1.8)", {
    expect_identical(valCellRuv("prot_dia", valArmRuvOff()), "off")
    expect_identical(valCellRuv("prot_dia", valArmShipped()), "correction")
    expect_identical(valCellRuv("prot_dia_limpa", valArmShipped()), "factors_in_model")
    expect_identical(valCellRuv("trx_screen", valArmShipped()), "factors_in_model")
    expect_identical(valCellRuv("trx_screen", valArmRuvOff()), "off")
    expect_identical(valCellRuv("lipid_canonical", valArmRuvLegacy()), "correction")
    expect_identical(valCellRuv("metab_lc", valArmRuvDesign("groups_all_features")),
                     "factors_in_model")
    expect_identical(valCellRuv("metab_lc", valArmRuvDesign("legacy")), "correction")
    expect_identical(valCellRuv("lipid_canonical", valArmRuvSupplied("c")), "factors_in_model")
})
