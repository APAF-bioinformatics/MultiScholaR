# Arms run by a ticket function (KIT-159 7.5: kinds `custom` and `published`; valRunCustomArm),
# with their labels, the published anchor check (anchor_on_fail, KIT-163 N6) and refusals (4.7).

armDa <- function() {
    data.frame(feature_id = c("k1", "k2"), assay = "all", contrast = "B_vs_A", log2fc = c(1, 0),
               p = c(0.01, 0.9), s_own = c(0.02, 0.9), call_statistic = "BH",
               called_own = c(TRUE, FALSE), mean_expression = c(5, 6))
}

test_that("a custom arm is labelled 'not the harness baseline'", {
    arm <- list(id = "welch", kind = "custom", fn = function(unit, inputs) armDa())
    res <- valRunCustomArm(arm, unit = "u01", inputs = list())
    expect_false(res$refused)
    expect_identical(res$label, "not the harness baseline")
    expect_identical(res$da, armDa())
})

test_that("a published arm: anchor reproduced, stopped, or kept as context", {
    pass <- list(id = "cox", kind = "published", fn = function(unit, inputs) armDa(),
                 anchor = function(da) TRUE)
    res <- valRunCustomArm(pass, unit = "u01", inputs = list())
    expect_identical(res$label, "published processing re-run")
    expect_identical(res$anchor, "reproduced")
    expect_true(res$primary)
    fail <- utils::modifyList(pass, list(anchor = function(da) FALSE))
    expect_error(valRunCustomArm(fail, unit = "u01", inputs = list()), class = "val_arm_refused")
    context <- utils::modifyList(fail, list(anchor_on_fail = "context"))
    res <- valRunCustomArm(context, unit = "u01", inputs = list())
    expect_identical(res$label, "anchor not reproduced: context only")
    expect_false(res$primary)
    expect_error(valRunCustomArm(list(id = "p", kind = "published", fn = function(unit, inputs) armDa()),
                                 unit = "u01", inputs = list()),
                 class = "val_spec_invalid")
})

test_that("a ticket function that refuses is recorded as refused with its class and digest", {
    arm <- list(id = "c", kind = "custom", fn = function(unit, inputs) {
        cli::cli_abort("too few replicates: 2", class = c("ticket_too_few", "rlang_error"))
    })
    res <- valRunCustomArm(arm, unit = "u07", inputs = list())
    expect_true(res$refused)
    expect_identical(res$class, "ticket_too_few")
    expect_match(res$message_sha256, "^[0-9a-f]{64}$")
    expect_null(res$da)
    expect_error(valRunCustomArm(list(id = "x", kind = "workflow"), unit = "u", inputs = list()),
                 class = "val_spec_invalid")
})
