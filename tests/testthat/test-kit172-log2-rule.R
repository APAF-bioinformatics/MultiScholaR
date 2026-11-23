# KIT-172: the log2 rule per offset (ticket section 6.1), the replay of a recorded offset, the defaults and the bound on
# intensity-scale data: tests T3, T4, T5 and T6 of section 6.4. Fixtures from helper-kit172.R.
library(testthat)

kit172SmallValues <- function() {
  values <- kit172Values(172001L)$values[c(1:40, 516:520), ]
  values[is.na(values)] <- 5
  values
}

kit172AsTibbles <- function(object, class) {
  slot <- kit172Binding(class)$data_slot
  methods::slot(object, slot) <- lapply(methods::slot(object, slot), tibble::as_tibble)
  object
}

kit172Assay <- function(object, class, assay = "LCMS_Pos") methods::slot(object, kit172Binding(class)$data_slot)[[assay]]

for (class in c("lipid", "metab")) {
  test_that(paste("KIT-172 T3:", class, "direct callers at offset 0: non-positive and non-finite cells missing, counted, one warning"), {
    values <- kit172SmallValues()
    bad <- values
    bad["F001", "A_1"] <- 0
    bad["F001", "A_2"] <- -3
    bad["F002", "A_1"] <- NaN
    bad["F002", "A_2"] <- Inf
    bad["F003", "A_1"] <- -Inf
    bad["F003", "A_2"] <- NA_real_
    object <- kit172AsTibbles(kit172Object(list(LCMS_Pos = bad, LCMS_Neg = values), class), class)
    caught <- list()
    result <- withCallingHandlers(suppressMessages(logTransformAssays(object, offset = 0)),
      omics_log2_nonpositive = function(w) {
        caught[[length(caught) + 1L]] <<- w
        invokeRestart("muffleWarning")
      })
    expect_length(caught, 1L)
    expect_identical(caught[[1L]]$counts$LCMS_Pos, c(nonpositive = 3L, nonfinite = 2L))
    expect_identical(caught[[1L]]$counts$LCMS_Neg, c(nonpositive = 0L, nonfinite = 0L))
    expect_identical(result@args$log_transform_nonpositive$LCMS_Pos, c(nonpositive = 3L, nonfinite = 2L))
    expect_identical(result@args$log_transform_nonpositive$LCMS_Neg, c(nonpositive = 0L, nonfinite = 0L))
    expect_identical(result@args$log_transform_offset, 0)
    assay <- kit172Assay(result, class)
    id <- methods::slot(result, kit172Binding(class)$id_slot)
    at <- function(feature, run) assay[[run]][match(feature, assay[[id]])]
    for (cell in list(c("F001", "A_1"), c("F001", "A_2"), c("F002", "A_1"), c("F002", "A_2"), c("F003", "A_1"), c("F003", "A_2"))) {
      expect_true(is.na(at(cell[[1L]], cell[[2L]])) && !is.nan(at(cell[[1L]], cell[[2L]])), info = paste(cell, collapse = " "))
    }
    runs <- kit172Design()$Run
    logged <- as.matrix(as.data.frame(assay)[, runs])
    expect_false(any(is.nan(logged) | is.infinite(logged)))
    # Every other cell is plain log2 of its value.
    plain <- suppressWarnings(log2(bad[assay[[id]], runs]))
    plain[!is.finite(plain)] <- NA_real_
    expect_identical(unname(logged), unname(plain))
  })

  test_that(paste("KIT-172 T3:", class, "an invalid offset is refused and a failed transform stops"), {
    object <- kit172AsTibbles(kit172Object(kit172SmallValues(), class), class)
    for (offset in list(-1, NA, c(1, 2), "1")) {
      expect_error(suppressMessages(logTransformAssays(object, offset = offset)), class = "omics_log2_offset_invalid",
        info = deparse(offset))
    }
    broken <- object
    slot <- kit172Binding(class)$data_slot
    assay <- methods::slot(broken, slot)$LCMS_Pos
    # A column of functions: the transform cannot coerce it to numbers (R1, run 45121678056: a column of one-element
    # lists was coerced without an error).
    assay$A_1 <- lapply(assay$A_1, function(value) function() value)
    methods::slot(broken, slot)$LCMS_Pos <- assay
    expect_error(suppressWarnings(suppressMessages(logTransformAssays(broken, offset = 0))), class = "omics_log2_failed")
  })

  test_that(paste("KIT-172 T4:", class, "an offset above 0 returns the object of earlier releases (method)"), {
    values <- kit172SmallValues()
    values["F001", "A_1"] <- 0
    values["F001", "A_2"] <- -3
    values["F002", "A_1"] <- NA_real_
    object <- kit172AsTibbles(kit172Object(values, class), class)
    runs <- kit172Design()$Run
    for (offset in c(1, 2)) {
      result <- suppressMessages(logTransformAssays(object, offset = offset))
      expected <- object
      slot <- kit172Binding(class)$data_slot
      assay <- methods::slot(expected, slot)$LCMS_Pos
      for (run in runs) assay[[run]] <- log2(ifelse(!is.na(assay[[run]]) & assay[[run]] < 0, 0, assay[[run]]) + offset)
      methods::slot(expected, slot)$LCMS_Pos <- assay
      expected@args$log_transformed <- TRUE
      expected@args$log_transform_offset <- offset
      expect_identical(result, expected, info = paste("offset", offset))
      expect_null(result@args$log_transform_nonpositive)
    }
  })

  test_that(paste("KIT-172 T4:", class, "the declared operation replays an offset above 0 and is the same in both variants at 0"), {
    skip_if_not(workflowOwnerCapable())
    values <- kit172SmallValues()
    values["F001", "A_1"] <- 0
    values["F001", "A_2"] <- -3
    values["F002", "A_1"] <- NA_real_
    object <- kit172Object(values, class)
    binding <- kit172Binding(class)
    runs <- kit172Design()$Run
    root <- withr::local_tempdir()
    owners <- list()
    withr::defer(for (owner in owners) closeDIAStateAdapter(owner))
    run <- function(offset, variant) {
      fixture <- kit172Owner(file.path(root, paste(class, offset, variant, sep = "-")), object, class)
      owners[[length(owners) + 1L]] <<- fixture$owner
      step <- omicsLog2Step(binding, offset)
      step$parameters$variant <- variant
      result <- kit172Quiet(runWorkflowOmicsReferenceStep(fixture$owner, step, fixture$workflow$config_list))
      list(result = result, state = fixture$owner$getState(paste0(binding$prefix, "_log2")))
    }
    for (offset in c(1, 2)) for (variant in c("whole", "partitioned")) {
      got <- run(offset, variant)
      info <- paste("offset", offset, variant)
      expect_null(got$result$value$nonpositive, info = info)
      expect_null(got$state@args$log_transform_nonpositive, info = info)
      expect_identical(got$state@args$log_transform_offset, offset, info = info)
      expect_true(got$state@args$log_transformed, info = info)
      assay <- methods::slot(got$state, binding$data_slot)$LCMS_Pos
      id <- methods::slot(got$state, binding$id_slot)
      input <- values[assay[[id]], , drop = FALSE]
      for (r in runs) {
        expect_identical(assay[[r]], log2(ifelse(!is.na(input[, r]) & input[, r] < 0, 0, input[, r]) + offset) |> unname(),
          info = paste(info, r))
      }
    }
    whole <- run(0, "whole")
    partitioned <- run(0, "partitioned")
    expect_identical(omicsDeclarationsScrub(whole$state), omicsDeclarationsScrub(partitioned$state))
    expect_null(whole$state@args$log_transform_nonpositive)
    expect_null(partitioned$state@args$log_transform_nonpositive)
    expect_identical(whole$result$value$nonpositive, partitioned$result$value$nonpositive)
    expect_identical(whole$result$value$negative_log, partitioned$result$value$negative_log)
    expect_identical(partitioned$result$value$nonpositive$LCMS_Pos, c(nonpositive = 2L, nonfinite = 0L))
    expect_identical(partitioned$result$value$negative_log$LCMS_Pos, sum(values > 0 & values < 1, na.rm = TRUE))
  })

  test_that(paste("KIT-172 T5:", class, "the defaults are plain log2 and the tab offers no offset"), {
    method <- methods::getMethod("logTransformAssays", if (identical(class, "lipid")) "LipidomicsAssayData" else "MetaboliteAssayData")
    expect_identical(formals(method)$offset, 0)
    expect_identical(formals(methods::getGeneric("logTransformAssays"))$offset, 0)
    config <- if (identical(class, "lipid")) omicGuiLipidNorm() else omicGuiMetabNorm()
    expect_identical(config$panels$log2$offset, 0)
    ui <- if (identical(class, "lipid")) mod_lipid_norm_ui("norm") else mod_metab_norm_ui("norm")
    html <- as.character(htmltools::renderTags(ui)$html)
    expect_false(grepl("norm-log_offset", html, fixed = TRUE))
    expect_true(grepl("norm-log2_recorded_offset", html, fixed = TRUE))
    settings_of <- function(input, ...) omicNormRunSettings(input, list(color_var = "group", shape_var = "group"),
      list(assay_names = "LCMS_Pos"), list(qc = tempdir()), list(qc_dir_role = "qc"), ...)
    expect_identical(settings_of(list(norm_method = "cyclicloess", ruv_mode = "skip"))$log_offset, 0)
    expect_identical(settings_of(list(norm_method = "cyclicloess", ruv_mode = "skip"), logOffset = config$panels$log2$offset)$log_offset, 0)
    expect_identical(settings_of(list(norm_method = "cyclicloess", ruv_mode = "skip", log_offset = 1))$log_offset, 1)
  })

  test_that(paste("KIT-172 T5:", class, "a default run's record holds offset 0"), {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    journey <- kit172Journey(file.path(root, "j"), kit172Object(kit172SmallValues(), class), class, itsd = FALSE,
      method = "none", filter = FALSE, da = FALSE)
    withr::defer(closeDIAStateAdapter(journey$owner))
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(journey$owner, json)
    overlay <- workflowRecordReportParameters(json, list(log_transform_offset = 99), list())
    expect_equal(overlay$log_transform_offset, 0)
    provenance <- journey$owner$.context$snapshot$step_provenance[[paste0(journey$prefix, "_log2")]]
    expect_identical(provenance$values$offset, 0)
  })

  test_that(paste("KIT-172 T6:", class, "on intensity-scale data plain log2 moves each value by log2(1 + 1/x), at most 1.4427e-4"), {
    set.seed(172006L)
    values <- kit172SmallValues()
    values[] <- 1e4 * (1 + stats::runif(length(values), 0, 99))
    object <- kit172AsTibbles(kit172Object(values, class), class)
    new <- kit172Assay(suppressMessages(logTransformAssays(object)), class)
    old <- kit172Assay(suppressMessages(logTransformAssays(object, offset = 1)), class)
    runs <- kit172Design()$Run
    difference <- as.matrix(as.data.frame(old)[, runs]) - as.matrix(as.data.frame(new)[, runs])
    expect_lte(max(abs(difference)), 1.4427e-4)
    expect_lte(max(abs(difference - log2(1 + 1 / values[, runs]))), 1e-12)
  })
}
