library(testthat)

source(test_path("..", "..", "R", "mod_lipid_design_builder_display_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design_builder_action_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design_builder_state_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design_builder_ui.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design_builder_server.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design_builder.R"), local = environment())

makeReactiveValMock <- function(initialValue) {
    currentValue <- initialValue

    function(newValue) {
        if (missing(newValue)) {
            return(currentValue)
        }

        currentValue <<- newValue
        invisible(newValue)
    }
}

test_that("formatLipidDesignTechRepSummary keeps grouped summary text stable", {
    design_matrix <- data.frame(
        Run = c("Sample_1", "Sample_2", "Sample_3"),
        group = c("Control", "Control", "Treatment"),
        replicates = c(1L, 1L, 2L),
        tech_reps = c(1L, 2L, NA_integer_),
        stringsAsFactors = FALSE
    )

    expect_identical(
        formatLipidDesignTechRepSummary(design_matrix),
        paste(
            "Group: Control, Biological Replicate: 1\n  Samples: Sample_1, Sample_2\n  Technical Replicates: 1, 2"
        )
    )
})

test_that("formatLipidDesignRemovedSamplesDisplay keeps empty and sorted output stable", {
    expect_identical(
        formatLipidDesignRemovedSamplesDisplay(character(0)),
        "No samples have been removed."
    )

    expect_identical(
        formatLipidDesignRemovedSamplesDisplay(c("Sample_10", "Sample_2")),
        "Removed 2 sample(s):\nSample_2\nSample_10"
    )
})

test_that("formatLipidDesignContrastFactorsInfo keeps the formula guidance stable", {
    expect_identical(
        formatLipidDesignContrastFactorsInfo("~ 0 + group"),
        "Note: Contrasts will use 'group' prefix based on formula: ~ 0 + group"
    )

    expect_identical(
        formatLipidDesignContrastFactorsInfo("~ group + batch"),
        "Note: Contrasts will use group names as-is"
    )
})

test_that("formatLipidDesignAvailableFactorsDisplay keeps empty and joined output stable", {
    expect_identical(
        formatLipidDesignAvailableFactorsDisplay(character(0)),
        "No factors defined yet (use the 'Factors' tab)."
    )

    expect_identical(
        formatLipidDesignAvailableFactorsDisplay(c("Treatment", "Timepoint")),
        "Treatment, Timepoint"
    )
})

test_that("formatLipidDesignDefinedContrastLines keeps contrast rendering lines stable", {
    contrast_data <- data.frame(
        contrast_name = c("Treatment.vs.Control", "Late.vs.Early"),
        numerator = c("Treatment", "Late"),
        denominator = c("Control", "Early"),
        stringsAsFactors = FALSE
    )

    expect_identical(
        formatLipidDesignDefinedContrastLines(contrast_data, "~ 0 + group"),
        c(
            "Treatment_vs_Control=groupTreatment-groupControl",
            "Late_vs_Early=groupLate-groupEarly"
        )
    )

    expect_identical(
        formatLipidDesignDefinedContrastLines(contrast_data, "~ group + batch"),
        c(
            "Treatment_vs_Control=Treatment-Control",
            "Late_vs_Early=Late-Early"
        )
    )
})

test_that("formatLipidDesignRangePreview keeps success and error text stable", {
    expect_identical(
        formatLipidDesignRangePreview(
            samplesToTransform = c("Sample_A", "Sample_B"),
            rangeStart = 1L,
            rangeEnd = 2L,
            extractExperimentFn = function(sampleName, mode, start, end) {
                paste(sampleName, mode, start, end, sep = "|")
            }
        ),
        "\"Sample_A\" -> \"Sample_A|range|1|2\""
    )

    expect_identical(
        formatLipidDesignRangePreview(
            samplesToTransform = c("Sample_A", "Sample_B"),
            rangeStart = 1L,
            rangeEnd = 2L,
            extractExperimentFn = function(...) {
                stop("bad range")
            }
        ),
        "Error: bad range"
    )
})

test_that("buildLipidDesignReplicateInputLabel keeps the replicate input label stable", {
    expect_identical(
        buildLipidDesignReplicateInputLabel(c("S1", "S2", "S3")),
        "Starting replicate number for 3 selected samples:"
    )
})

test_that("appendLipidDesignFactorName keeps trimming and duplicate checks stable", {
    expect_identical(
        appendLipidDesignFactorName(c("Treatment"), "  Timepoint  "),
        c("Treatment", "Timepoint")
    )

    expect_identical(
        appendLipidDesignFactorName(c("Treatment"), " Treatment "),
        c("Treatment")
    )

    expect_identical(
        appendLipidDesignFactorName(c("Treatment"), "   "),
        c("Treatment")
    )
})

test_that("buildLipidDesignMetadataReplicateNumbers keeps replicate numbering stable", {
    expect_identical(
        buildLipidDesignMetadataReplicateNumbers(c("S1", "S2", "S3"), 4L),
        c(4L, 5L, 6L)
    )

    expect_identical(
        buildLipidDesignMetadataReplicateNumbers(c("S1", "S2"), NULL),
        NA_integer_
    )
})

test_that("buildLipidDesignMetadataGroupName keeps factor-based group labels stable", {
    expect_identical(
        buildLipidDesignMetadataGroupName("Treatment", "", "Day1"),
        "Treatment_Day1"
    )

    expect_identical(
        buildLipidDesignMetadataGroupName("", "", NULL),
        NA_character_
    )
})

test_that("listLipidDesignAssignedGroups keeps non-empty group collection stable", {
    design_matrix <- data.frame(
        group = c(NA_character_, "Control", "", "Treatment"),
        stringsAsFactors = FALSE
    )

    expect_identical(
        listLipidDesignAssignedGroups(design_matrix),
        c("Control", "Treatment")
    )
})

test_that("applyLipidDesignMetadataAssignment keeps metadata bookkeeping stable", {
    design_matrix <- data.frame(
        Run = c("Sample_1", "Sample_2", "Sample_3"),
        factor1 = c(NA_character_, NA_character_, NA_character_),
        factor2 = c(NA_character_, NA_character_, NA_character_),
        factor3 = c(NA_character_, NA_character_, NA_character_),
        group = c(NA_character_, NA_character_, NA_character_),
        replicates = c(NA_integer_, NA_integer_, NA_integer_),
        stringsAsFactors = FALSE
    )

    result <- applyLipidDesignMetadataAssignment(
        designMatrix = design_matrix,
        selectedRuns = c("Sample_1", "Sample_3"),
        factor1Selection = "Control",
        factor2Selection = "",
        factor3Selection = "Day1",
        replicateStart = 4L
    )

    expect_identical(result$factor1, c("Control", NA_character_, "Control"))
    expect_identical(result$factor2, c("", NA_character_, ""))
    expect_identical(result$factor3, c("Day1", NA_character_, "Day1"))
    expect_identical(result$replicates, c(4L, NA_integer_, 5L))
    expect_identical(result$group, c("Control_Day1", NA_character_, "Control_Day1"))
})

test_that("registerLipidDesignSummaryOutputShells keeps the summary render handoff stable", {
    output <- new.env(parent = emptyenv())
    render_calls <- character()

    result <- registerLipidDesignSummaryOutputShells(
        output = output,
        designMatrix = function() {
            data.frame(
                Run = c("Sample_1", "Sample_2"),
                group = c("Control", "Control"),
                replicates = c(1L, 1L),
                tech_reps = c(1L, 2L),
                stringsAsFactors = FALSE
            )
        },
        removedSamples = function() c("Sample_10", "Sample_2"),
        formulaString = function() "~ 0 + group",
        renderTextFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(
                list(kind = "mock_render", index = length(render_calls)),
                class = "mock_render"
            )
        },
        reqFn = function(...) invisible(NULL)
    )

    expect_identical(result, output)
    expect_identical(
        render_calls,
        c(
            "Group: Control, Biological Replicate: 1\n  Samples: Sample_1, Sample_2\n  Technical Replicates: 1, 2",
            "Removed 2 sample(s):\nSample_2\nSample_10",
            "Note: Contrasts will use 'group' prefix based on formula: ~ 0 + group"
        )
    )
    expect_s3_class(output$tech_rep_summary, "mock_render")
    expect_s3_class(output$removed_samples_display, "mock_render")
    expect_s3_class(output$contrast_factors_info, "mock_render")
})

test_that("registerLipidDesignAdjacentOutputShells keeps adjacent render handoff stable", {
    output <- new.env(parent = emptyenv())
    ui_calls <- list()
    text_calls <- character()

    result <- registerLipidDesignAdjacentOutputShells(
        output = output,
        session = list(ns = function(id) paste0("ns-", id)),
        factors = function() c("Treatment", "Timepoint"),
        contrasts = function() {
            data.frame(
                contrast_name = "Treatment.vs.Control",
                numerator = "Treatment",
                denominator = "Control",
                stringsAsFactors = FALSE
            )
        },
        formulaString = function() "~ 0 + group",
        samplesToTransform = function() c("Sample_A", "Sample_B"),
        rangeStart = function() 1L,
        rangeEnd = function() 2L,
        selectedRuns = function() c("Sample_A", "Sample_B", "Sample_C"),
        renderUiFn = function(expr) {
            ui_calls[[length(ui_calls) + 1]] <<- force(expr)
            structure(
                list(kind = "mock_ui", index = length(ui_calls)),
                class = "mock_ui"
            )
        },
        renderTextFn = function(expr) {
            text_calls <<- c(text_calls, force(expr))
            structure(
                list(kind = "mock_text", index = length(text_calls)),
                class = "mock_text"
            )
        },
        reqFn = function(...) invisible(NULL),
        paragraphFn = function(value) list(kind = "p", value = value),
        tagListFn = function(value) list(kind = "tagList", value = value),
        codeTagFn = function(value) list(kind = "code", value = value),
        extractExperimentFn = function(sampleName, mode, start, end) {
            paste(sampleName, mode, start, end, sep = "|")
        },
        numericInputFn = function(inputId, label, value, min) {
            list(
                kind = "numericInput",
                inputId = inputId,
                label = label,
                value = value,
                min = min
            )
        }
    )

    expect_identical(result, output)
    expect_identical(
        ui_calls,
        list(
            list(kind = "p", value = "Treatment, Timepoint"),
            list(
                kind = "tagList",
                value = list(
                    list(
                        kind = "p",
                        value = list(
                            kind = "code",
                            value = "Treatment_vs_Control=groupTreatment-groupControl"
                        )
                    )
                )
            ),
            list(
                kind = "numericInput",
                inputId = "ns-replicate_start",
                label = "Starting replicate number for 3 selected samples:",
                value = 1,
                min = 1
            )
        )
    )
    expect_identical(text_calls, "\"Sample_A\" -> \"Sample_A|range|1|2\"")
    expect_s3_class(output$available_factors_display, "mock_ui")
    expect_s3_class(output$defined_contrasts_display, "mock_ui")
    expect_s3_class(output$range_preview, "mock_text")
    expect_s3_class(output$replicate_inputs, "mock_ui")
})

test_that("applyLipidDesignSampleRenameMap keeps shared rename bookkeeping stable", {
    design_matrix <- data.frame(
        Run = c("Sample_1", "Sample_2", "Sample_3"),
        stringsAsFactors = FALSE
    )
    assay_list <- list(
        LCMS_Pos = data.frame(
            feature = c("L1", "L2"),
            Sample_1 = c(10, 20),
            Sample_3 = c(30, 40),
            check.names = FALSE
        ),
        LCMS_Neg = data.frame(
            feature = c("L3", "L4"),
            Sample_2 = c(50, 60),
            Sample_3 = c(70, 80),
            check.names = FALSE
        )
    )
    sample_names <- c("Sample_1", "Sample_2", "Sample_3")

    result <- applyLipidDesignSampleRenameMap(
        designMatrix = design_matrix,
        assayList = assay_list,
        sampleNames = sample_names,
        renameMap = c(Sample_1 = "Renamed_1", Sample_3 = "Range_3")
    )

    expect_identical(
        result$designMatrix$Run,
        c("Renamed_1", "Sample_2", "Range_3")
    )
    expect_identical(
        names(result$assayList$LCMS_Pos),
        c("feature", "Renamed_1", "Range_3")
    )
    expect_identical(
        names(result$assayList$LCMS_Neg),
        c("feature", "Sample_2", "Range_3")
    )
    expect_identical(
        result$sampleNames,
        c("Renamed_1", "Sample_2", "Range_3")
    )
})

test_that("buildLipidDesignBulkRenameMap keeps transform dispatch stable", {
    extract_mock <- function(sampleName, mode, start = NULL, end = NULL) {
        if (identical(mode, "range")) {
            return(paste(sampleName, mode, start, end, sep = "|"))
        }

        paste(sampleName, mode, sep = "|")
    }

    expect_identical(
        buildLipidDesignBulkRenameMap(
            samplesToTransform = c("Sample_A", "Sample_B"),
            transformMode = "range",
            rangeStart = 1L,
            rangeEnd = 2L,
            extractExperimentFn = extract_mock,
            reqFn = function(...) invisible(NULL)
        ),
        c(Sample_A = "Sample_A|range|1|2", Sample_B = "Sample_B|range|1|2")
    )

    expect_identical(
        buildLipidDesignBulkRenameMap(
            samplesToTransform = "Sample_A",
            transformMode = "before_underscore",
            extractExperimentFn = extract_mock,
            reqFn = function(...) invisible(NULL)
        ),
        c(Sample_A = "Sample_A|start")
    )

    expect_identical(
        buildLipidDesignBulkRenameMap(
            samplesToTransform = "Sample_A",
            transformMode = "after_underscore",
            extractExperimentFn = extract_mock,
            reqFn = function(...) invisible(NULL)
        ),
        c(Sample_A = "Sample_A|end")
    )
})

test_that("registerLipidDesignSampleRenameShells keeps rename observer handoff stable", {
    observed_labels <- character()
    updated_inputs <- list()

    design_matrix <- makeReactiveValMock(data.frame(
        Run = c("Sample_1", "Sample_2"),
        stringsAsFactors = FALSE
    ))
    assay_list <- makeReactiveValMock(list(
        LCMS_Pos = data.frame(
            feature = "L1",
            Sample_1 = 10,
            Sample_2 = 20,
            check.names = FALSE
        )
    ))
    sample_names <- makeReactiveValMock(c("Sample_1", "Sample_2"))

    registerLipidDesignSampleRenameShells(
        input = list(
            rename_sample = 1L,
            sample_to_rename = "Sample_1",
            new_sample_name = "Renamed_1",
            bulk_rename = 0L,
            samples_to_transform = c("Sample_2"),
            transform_mode = "before_underscore",
            range_start = NULL,
            range_end = NULL
        ),
        session = "mock-session",
        designMatrix = design_matrix,
        dataClnReactive = assay_list,
        sampleNamesReactive = sample_names,
        observeEventFn = function(eventExpr, handlerExpr) {
            event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
            observed_labels <<- c(observed_labels, event_label)

            if (identical(event_label, "input$rename_sample")) {
                eval(substitute(handlerExpr), envir = parent.frame())
            }

            invisible(NULL)
        },
        reqFn = function(...) invisible(NULL),
        updateTextInputFn = function(session, inputId, value) {
            updated_inputs[[length(updated_inputs) + 1]] <<- list(
                session = session,
                inputId = inputId,
                value = value
            )
            invisible(NULL)
        }
    )

    expect_identical(observed_labels, c("input$rename_sample", "input$bulk_rename"))
    expect_identical(design_matrix()$Run, c("Renamed_1", "Sample_2"))
    expect_identical(names(assay_list()$LCMS_Pos), c("feature", "Renamed_1", "Sample_2"))
    expect_identical(sample_names(), c("Renamed_1", "Sample_2"))
    expect_identical(
        updated_inputs,
        list(list(session = "mock-session", inputId = "new_sample_name", value = ""))
    )
})

test_that("registerLipidDesignSampleRenameShells keeps bulk rename observer handoff stable", {
    observed_labels <- character()
    build_bulk_rename_map_calls <- list()

    design_matrix <- makeReactiveValMock(data.frame(
        Run = c("Sample_1", "Sample_2", "Sample_3"),
        stringsAsFactors = FALSE
    ))
    assay_list <- makeReactiveValMock(list(
        LCMS_Pos = data.frame(
            feature = "L1",
            Sample_1 = 10,
            Sample_2 = 20,
            Sample_3 = 30,
            check.names = FALSE
        )
    ))
    sample_names <- makeReactiveValMock(c("Sample_1", "Sample_2", "Sample_3"))

    registerLipidDesignSampleRenameShells(
        input = list(
            rename_sample = 0L,
            sample_to_rename = "Sample_3",
            new_sample_name = "Ignored",
            bulk_rename = 1L,
            samples_to_transform = c("Sample_1", "Sample_2"),
            transform_mode = "range",
            range_start = 1L,
            range_end = 2L
        ),
        session = "mock-session",
        designMatrix = design_matrix,
        dataClnReactive = assay_list,
        sampleNamesReactive = sample_names,
        observeEventFn = function(eventExpr, handlerExpr) {
            event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
            observed_labels <<- c(observed_labels, event_label)

            if (identical(event_label, "input$bulk_rename")) {
                eval(substitute(handlerExpr), envir = parent.frame())
            }

            invisible(NULL)
        },
        reqFn = function(...) invisible(NULL),
        updateTextInputFn = function(...) invisible(NULL),
        buildBulkRenameMapFn = function(samplesToTransform, transformMode, rangeStart, rangeEnd) {
            build_bulk_rename_map_calls[[length(build_bulk_rename_map_calls) + 1]] <<- list(
                samplesToTransform = samplesToTransform,
                transformMode = transformMode,
                rangeStart = rangeStart,
                rangeEnd = rangeEnd
            )

            c(Sample_1 = "Range_1", Sample_2 = "Range_2")
        }
    )

    expect_identical(observed_labels, c("input$rename_sample", "input$bulk_rename"))
    expect_identical(
        build_bulk_rename_map_calls,
        list(list(
            samplesToTransform = c("Sample_1", "Sample_2"),
            transformMode = "range",
            rangeStart = 1L,
            rangeEnd = 2L
        ))
    )
    expect_identical(design_matrix()$Run, c("Range_1", "Range_2", "Sample_3"))
    expect_identical(names(assay_list()$LCMS_Pos), c("feature", "Range_1", "Range_2", "Sample_3"))
    expect_identical(sample_names(), c("Range_1", "Range_2", "Sample_3"))
})

test_that("registerLipidDesignMetadataAssignmentShells keeps add-factor observer handoff stable", {
    observed_labels <- character()
    updated_inputs <- list()

    design_matrix <- makeReactiveValMock(data.frame(
        Run = c("Sample_1", "Sample_2"),
        stringsAsFactors = FALSE
    ))
    factors_rv <- makeReactiveValMock(c("Treatment"))
    groups_rv <- makeReactiveValMock(character(0))

    registerLipidDesignMetadataAssignmentShells(
        input = list(
            add_factor = 1L,
            new_factor = "  Timepoint  ",
            assign_metadata = 0L,
            selected_runs = c("Sample_1"),
            factor1_select = "Ignored",
            factor2_select = "",
            factor3_select = "",
            replicate_start = 1L
        ),
        session = "mock-session",
        designMatrix = design_matrix,
        factors = factors_rv,
        groups = groups_rv,
        observeEventFn = function(eventExpr, handlerExpr) {
            event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
            observed_labels <<- c(observed_labels, event_label)

            if (identical(event_label, "input$add_factor")) {
                eval(substitute(handlerExpr), envir = parent.frame())
            }

            invisible(NULL)
        },
        reqFn = function(...) invisible(NULL),
        updateTextInputFn = function(session, inputId, value) {
            updated_inputs[[length(updated_inputs) + 1]] <<- list(
                session = session,
                inputId = inputId,
                value = value
            )
            invisible(NULL)
        }
    )

    expect_identical(observed_labels, c("input$add_factor", "input$assign_metadata"))
    expect_identical(factors_rv(), c("Treatment", "Timepoint"))
    expect_identical(groups_rv(), character(0))
    expect_identical(
        updated_inputs,
        list(list(session = "mock-session", inputId = "new_factor", value = ""))
    )
})

test_that("registerLipidDesignMetadataAssignmentShells keeps assign-metadata observer handoff stable", {
    observed_labels <- character()

    design_matrix <- makeReactiveValMock(data.frame(
        Run = c("Sample_1", "Sample_2", "Sample_3"),
        factor1 = c(NA_character_, NA_character_, NA_character_),
        factor2 = c(NA_character_, NA_character_, NA_character_),
        factor3 = c(NA_character_, NA_character_, NA_character_),
        group = c(NA_character_, NA_character_, "Existing_Group"),
        replicates = c(NA_integer_, NA_integer_, 9L),
        stringsAsFactors = FALSE
    ))
    factors_rv <- makeReactiveValMock(c("Control", "Batch1"))
    groups_rv <- makeReactiveValMock("Existing_Group")

    registerLipidDesignMetadataAssignmentShells(
        input = list(
            add_factor = 0L,
            new_factor = "Ignored",
            assign_metadata = 1L,
            selected_runs = c("Sample_1", "Sample_2"),
            factor1_select = "Control",
            factor2_select = "Batch1",
            factor3_select = "",
            replicate_start = 5L
        ),
        session = "mock-session",
        designMatrix = design_matrix,
        factors = factors_rv,
        groups = groups_rv,
        observeEventFn = function(eventExpr, handlerExpr) {
            event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
            observed_labels <<- c(observed_labels, event_label)

            if (identical(event_label, "input$assign_metadata")) {
                eval(substitute(handlerExpr), envir = parent.frame())
            }

            invisible(NULL)
        },
        reqFn = function(...) invisible(NULL),
        updateTextInputFn = function(...) invisible(NULL)
    )

    expect_identical(observed_labels, c("input$add_factor", "input$assign_metadata"))
    expect_identical(design_matrix()$factor1, c("Control", "Control", NA_character_))
    expect_identical(design_matrix()$factor2, c("Batch1", "Batch1", NA_character_))
    expect_identical(design_matrix()$factor3, c(NA_character_, NA_character_, NA_character_))
    expect_identical(design_matrix()$replicates, c(5L, 6L, 9L))
    expect_identical(design_matrix()$group, c("Control_Batch1", "Control_Batch1", "Existing_Group"))
    expect_identical(groups_rv(), c("Control_Batch1", "Existing_Group"))
})

test_that("applyLipidDesignTechRepAssignment keeps consolidation and validation stable", {
    design_matrix <- data.frame(
        Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
        group = c("Control", "Control", "Control", "Treatment"),
        replicates = c(1L, 2L, 3L, 1L),
        tech_reps = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_),
        stringsAsFactors = FALSE
    )

    result <- applyLipidDesignTechRepAssignment(
        designMatrix = design_matrix,
        techRepSamples = c("Sample_1", "Sample_2"),
        assignmentMode = "lowest"
    )

    expect_true(result$ok)
    expect_identical(result$type, "message")
    expect_identical(result$message, "Assigned 2 samples as technical replicates")
    expect_equal(result$designMatrix$replicates, c(1, 1, 2, 1))
    expect_identical(result$designMatrix$tech_reps, c(1L, 2L, NA_integer_, NA_integer_))

    warning_result <- applyLipidDesignTechRepAssignment(
        designMatrix = design_matrix,
        techRepSamples = "Sample_1",
        assignmentMode = "lowest"
    )

    expect_false(warning_result$ok)
    expect_identical(warning_result$type, "warning")
    expect_identical(warning_result$message, "Please select at least two samples.")
})

test_that("registerLipidDesignTechRepAssignmentShells keeps warning observer handoff stable", {
    observed_labels <- character()
    notifications <- list()
    design_matrix <- makeReactiveValMock(data.frame(
        Run = c("Sample_1", "Sample_2"),
        group = c("Control", "Treatment"),
        replicates = c(1L, 1L),
        tech_reps = c(NA_integer_, NA_integer_),
        stringsAsFactors = FALSE
    ))

    registerLipidDesignTechRepAssignmentShells(
        input = list(
            assign_tech_reps = 1L,
            tech_rep_samples = "Sample_1",
            tech_rep_assignment_mode = "lowest",
            manual_replicate_number = NULL
        ),
        designMatrix = design_matrix,
        observeEventFn = function(eventExpr, handlerExpr) {
            event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
            observed_labels <<- c(observed_labels, event_label)
            eval(substitute(handlerExpr), envir = parent.frame())
            invisible(NULL)
        },
        reqFn = function(...) invisible(NULL),
        showNotificationFn = function(message, type) {
            notifications[[length(notifications) + 1]] <<- list(message = message, type = type)
            invisible(NULL)
        }
    )

    expect_identical(observed_labels, "input$assign_tech_reps")
    expect_identical(
        notifications,
        list(list(message = "Please select at least two samples.", type = "warning"))
    )
    expect_identical(design_matrix()$replicates, c(1L, 1L))
    expect_identical(design_matrix()$tech_reps, c(NA_integer_, NA_integer_))
})

test_that("registerLipidDesignTechRepAssignmentShells keeps success observer handoff stable", {
    observed_labels <- character()
    notifications <- list()
    apply_calls <- list()
    design_matrix <- makeReactiveValMock(data.frame(
        Run = c("Sample_1", "Sample_2", "Sample_3"),
        group = c("Control", "Control", "Control"),
        replicates = c(1L, 2L, 3L),
        tech_reps = c(NA_integer_, NA_integer_, NA_integer_),
        stringsAsFactors = FALSE
    ))

    registerLipidDesignTechRepAssignmentShells(
        input = list(
            assign_tech_reps = 1L,
            tech_rep_samples = c("Sample_1", "Sample_2"),
            tech_rep_assignment_mode = "manual",
            manual_replicate_number = 5L
        ),
        designMatrix = design_matrix,
        observeEventFn = function(eventExpr, handlerExpr) {
            event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
            observed_labels <<- c(observed_labels, event_label)
            eval(substitute(handlerExpr), envir = parent.frame())
            invisible(NULL)
        },
        reqFn = function(...) invisible(NULL),
        showNotificationFn = function(message, type) {
            notifications[[length(notifications) + 1]] <<- list(message = message, type = type)
            invisible(NULL)
        },
        applyTechRepAssignmentFn = function(designMatrix, techRepSamples, assignmentMode, manualReplicateNumber) {
            apply_calls[[length(apply_calls) + 1]] <<- list(
                designMatrix = designMatrix,
                techRepSamples = techRepSamples,
                assignmentMode = assignmentMode,
                manualReplicateNumber = manualReplicateNumber
            )

            list(
                ok = TRUE,
                type = "message",
                message = "Assigned 2 samples as technical replicates",
                designMatrix = transform(
                    designMatrix,
                    replicates = c(5L, 5L, 2L),
                    tech_reps = c(1L, 2L, NA_integer_)
                )
            )
        }
    )

    expect_identical(observed_labels, "input$assign_tech_reps")
    expect_identical(length(apply_calls), 1L)
    expect_identical(apply_calls[[1]]$techRepSamples, c("Sample_1", "Sample_2"))
    expect_identical(apply_calls[[1]]$assignmentMode, "manual")
    expect_identical(apply_calls[[1]]$manualReplicateNumber, 5L)
    expect_identical(design_matrix()$replicates, c(5L, 5L, 2L))
    expect_identical(design_matrix()$tech_reps, c(1L, 2L, NA_integer_))
    expect_identical(
        notifications,
        list(list(message = "Assigned 2 samples as technical replicates", type = "message"))
    )
})

test_that("appendLipidDesignContrast keeps duplicate suppression and row append stable", {
    existing_contrasts <- data.frame(
        contrast_name = "Treatment.vs.Control",
        numerator = "Treatment",
        denominator = "Control",
        stringsAsFactors = FALSE
    )

    expect_identical(
        appendLipidDesignContrast(
            currentContrasts = existing_contrasts,
            group1 = "Treatment",
            group2 = "Control"
        ),
        existing_contrasts
    )

    expect_identical(
        appendLipidDesignContrast(
            currentContrasts = existing_contrasts,
            group1 = "Treatment",
            group2 = "Treatment"
        ),
        existing_contrasts
    )

    expect_identical(
        appendLipidDesignContrast(
            currentContrasts = existing_contrasts,
            group1 = "Responder",
            group2 = "NonResponder"
        ),
        rbind(
            existing_contrasts,
            data.frame(
                contrast_name = "Responder.vs.NonResponder",
                numerator = "Responder",
                denominator = "NonResponder",
                stringsAsFactors = FALSE
            )
        )
    )
})

test_that("registerLipidDesignContrastManagementShells keeps contrast observer handoff stable", {
    observed_labels <- character()
    append_calls <- list()
    contrasts_rv <- makeReactiveValMock(data.frame(
        contrast_name = "Treatment.vs.Control",
        numerator = "Treatment",
        denominator = "Control",
        stringsAsFactors = FALSE
    ))

    registerLipidDesignContrastManagementShells(
        input = list(
            add_contrast = 1L,
            contrast_group1 = "Responder",
            contrast_group2 = "NonResponder"
        ),
        contrasts = contrasts_rv,
        observeEventFn = function(eventExpr, handlerExpr) {
            event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
            observed_labels <<- c(observed_labels, event_label)
            eval(substitute(handlerExpr), envir = parent.frame())
            invisible(NULL)
        },
        reqFn = function(...) invisible(NULL),
        appendContrastFn = function(currentContrasts, group1, group2) {
            append_calls[[length(append_calls) + 1]] <<- list(
                currentContrasts = currentContrasts,
                group1 = group1,
                group2 = group2
            )

            rbind(
                currentContrasts,
                data.frame(
                    contrast_name = "Responder.vs.NonResponder",
                    numerator = group1,
                    denominator = group2,
                    stringsAsFactors = FALSE
                )
            )
        }
    )

    expect_identical(observed_labels, "input$add_contrast")
    expect_identical(length(append_calls), 1L)
    expect_identical(append_calls[[1]]$group1, "Responder")
    expect_identical(append_calls[[1]]$group2, "NonResponder")
    expect_identical(
        contrasts_rv(),
        data.frame(
            contrast_name = c("Treatment.vs.Control", "Responder.vs.NonResponder"),
            numerator = c("Treatment", "Responder"),
            denominator = c("Control", "NonResponder"),
            stringsAsFactors = FALSE
        )
    )
})

test_that("appendLipidDesignRemovedSamples keeps accumulation order and uniqueness stable", {
    expect_identical(
        appendLipidDesignRemovedSamples(
            currentRemovedSamples = c("Sample_1", "Sample_2"),
            selectedSamples = c("Sample_2", "Sample_3")
        ),
        c("Sample_1", "Sample_2", "Sample_3")
    )
})

test_that("registerLipidDesignSampleRemovalShells keeps remove-samples observer handoff stable", {
    observed_labels <- character()
    append_calls <- list()
    update_calls <- list()
    notifications <- list()
    removed_samples_rv <- makeReactiveValMock(c("Sample_1"))

    registerLipidDesignSampleRemovalShells(
        input = list(
            remove_samples = 1L,
            samples_to_remove = c("Sample_2", "Sample_3")
        ),
        removedSamples = removed_samples_rv,
        session = "session-token",
        observeEventFn = function(eventExpr, handlerExpr) {
            event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
            observed_labels <<- c(observed_labels, event_label)
            eval(substitute(handlerExpr), envir = parent.frame())
            invisible(NULL)
        },
        reqFn = function(...) invisible(NULL),
        updateSelectizeInputFn = function(session, inputId, selected) {
            update_calls[[length(update_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                selected = selected
            )
            invisible(NULL)
        },
        showNotificationFn = function(message, type) {
            notifications[[length(notifications) + 1]] <<- list(message = message, type = type)
            invisible(NULL)
        },
        appendRemovedSamplesFn = function(currentRemovedSamples, selectedSamples) {
            append_calls[[length(append_calls) + 1]] <<- list(
                currentRemovedSamples = currentRemovedSamples,
                selectedSamples = selectedSamples
            )

            c(currentRemovedSamples, selectedSamples)
        }
    )

    expect_identical(observed_labels, "input$remove_samples")
    expect_identical(length(append_calls), 1L)
    expect_identical(append_calls[[1]]$currentRemovedSamples, "Sample_1")
    expect_identical(append_calls[[1]]$selectedSamples, c("Sample_2", "Sample_3"))
    expect_identical(
        removed_samples_rv(),
        c("Sample_1", "Sample_2", "Sample_3")
    )
    expect_identical(
        update_calls,
        list(list(session = "session-token", inputId = "samples_to_remove", selected = ""))
    )
    expect_identical(
        notifications,
        list(list(message = "Removed 2 sample(s) from analysis.", type = "message"))
    )
})

test_that("buildLipidDesignResetConfirmationBody keeps reset messaging stable", {
    expect_identical(
        buildLipidDesignResetConfirmationBody(
            resetScope = "formula",
            htmlFn = function(value) list(kind = "html", value = value)
        ),
        list(
            kind = "html",
            value = "<p>This will revert <strong>formula</strong> to their initial state.</p>"
        )
    )
})

test_that("buildLipidDesignResetConfirmationModal keeps modal assembly stable", {
    expect_identical(
        buildLipidDesignResetConfirmationModal(
            resetScope = "all",
            nsFn = function(id) paste0("ns-", id),
            modalDialogFn = function(..., title, footer, easyClose) {
                list(
                    kind = "modalDialog",
                    body = list(...),
                    title = title,
                    footer = footer,
                    easyClose = easyClose
                )
            },
            tagListFn = function(...) list(kind = "tagList", value = list(...)),
            modalButtonFn = function(label) list(kind = "modalButton", label = label),
            actionButtonFn = function(inputId, label, class) {
                list(kind = "actionButton", inputId = inputId, label = label, class = class)
            },
            buildResetConfirmationBodyFn = function(resetScope) {
                list(kind = "body", resetScope = resetScope)
            }
        ),
        list(
            kind = "modalDialog",
            body = list(list(kind = "body", resetScope = "all")),
            title = "Confirm Reset",
            footer = list(
                kind = "tagList",
                value = list(
                    list(kind = "modalButton", label = "Cancel"),
                    list(
                        kind = "actionButton",
                        inputId = "ns-confirm_reset",
                        label = "Reset",
                        class = "btn-danger"
                    )
                )
            ),
            easyClose = TRUE
        )
    )
})

test_that("registerLipidDesignResetRequestShells keeps reset-request observer handoff stable", {
    observed_labels <- character()
    modal_calls <- list()
    shown_modals <- list()

    result <- registerLipidDesignResetRequestShells(
        input = list(
            reset_changes = 1L,
            reset_scope = "removed_samples"
        ),
        session = list(ns = function(id) paste0("ns-", id)),
        observeEventFn = function(eventExpr, handlerExpr) {
            event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
            observed_labels <<- c(observed_labels, event_label)
            eval(substitute(handlerExpr), envir = parent.frame())
            invisible(NULL)
        },
        showModalFn = function(modal) {
            shown_modals[[length(shown_modals) + 1]] <<- modal
            invisible(NULL)
        },
        buildResetConfirmationModalFn = function(resetScope, nsFn) {
            modal_calls[[length(modal_calls) + 1]] <<- list(
                resetScope = resetScope,
                confirmId = nsFn("confirm_reset")
            )

            list(kind = "modal", resetScope = resetScope)
        }
    )

    expect_null(result)
    expect_identical(observed_labels, "input$reset_changes")
    expect_identical(
        modal_calls,
        list(list(resetScope = "removed_samples", confirmId = "ns-confirm_reset"))
    )
    expect_identical(
        shown_modals,
        list(list(kind = "modal", resetScope = "removed_samples"))
    )
})

test_that("applyLipidDesignResetState keeps scope-based reset bookkeeping stable", {
    initial_design_matrix <- data.frame(
        Run = "Imported_1",
        group = "ImportedGroup",
        factor1 = "ImportedFactor1",
        factor2 = "ImportedFactor2",
        factor3 = "ImportedFactor3",
        batch = "Batch_A",
        tech_reps = 7L,
        stringsAsFactors = FALSE
    )
    design_matrix_rv <- makeReactiveValMock(data.frame(
        Run = "Current_1",
        group = "CurrentGroup",
        factor1 = "CurrentFactor1",
        factor2 = "CurrentFactor2",
        factor3 = "CurrentFactor3",
        batch = "Batch_Current",
        tech_reps = 3L,
        stringsAsFactors = FALSE
    ))
    data_cln_rv <- makeReactiveValMock(list(current = data.frame(value = 1)))
    sample_names_rv <- makeReactiveValMock("Current_1")
    removed_samples_rv <- makeReactiveValMock(c("Current_2", "Current_3"))
    factors_rv <- makeReactiveValMock("CurrentFactor")
    groups_rv <- makeReactiveValMock("CurrentGroup")
    contrasts_rv <- makeReactiveValMock(data.frame(
        contrast_name = "Current",
        numerator = "A",
        denominator = "B",
        stringsAsFactors = FALSE
    ))
    formula_updates <- character()

    result <- applyLipidDesignResetState(
        scope = "all",
        initialState = list(
            design_matrix = initial_design_matrix,
            data_cln = list(imported = data.frame(value = 2)),
            sample_names = "Imported_1",
            factors = c("ImportedFactor1", "ImportedFactor2"),
            groups = c("ImportedGroup", "ReferenceGroup"),
            contrasts = data.frame(
                contrast_name = "Imported",
                numerator = "Treatment",
                denominator = "Control",
                stringsAsFactors = FALSE
            ),
            formula = "~ group + batch"
        ),
        designMatrix = design_matrix_rv,
        dataClnReactive = data_cln_rv,
        sampleNamesReactive = sample_names_rv,
        removedSamples = removed_samples_rv,
        factors = factors_rv,
        groups = groups_rv,
        contrasts = contrasts_rv,
        updateFormulaFn = function(value) {
            formula_updates <<- c(formula_updates, value)
            invisible(value)
        }
    )

    expect_null(result)
    expect_identical(
        design_matrix_rv(),
        transform(
            initial_design_matrix,
            group = NA_character_,
            factor1 = NA_character_,
            factor2 = NA_character_,
            factor3 = NA_character_,
            tech_reps = NA_integer_
        )
    )
    expect_identical(data_cln_rv(), list(imported = data.frame(value = 2)))
    expect_identical(sample_names_rv(), "Imported_1")
    expect_identical(removed_samples_rv(), character(0))
    expect_identical(factors_rv(), c("ImportedFactor1", "ImportedFactor2"))
    expect_identical(groups_rv(), c("ImportedGroup", "ReferenceGroup"))
    expect_identical(
        contrasts_rv(),
        data.frame(
            contrast_name = "Imported",
            numerator = "Treatment",
            denominator = "Control",
            stringsAsFactors = FALSE
        )
    )
    expect_identical(formula_updates, "~ group + batch")
})

test_that("runLipidDesignResetConfirmationShell keeps apply and notification handoff stable", {
    apply_calls <- list()
    formula_updates <- list()
    remove_modal_calls <- 0L
    notifications <- list()

    result <- runLipidDesignResetConfirmationShell(
        scope = "formula",
        initialState = list(formula = "~ group"),
        designMatrix = "design-matrix-rv",
        dataClnReactive = "data-cln-rv",
        sampleNamesReactive = "sample-names-rv",
        removedSamples = "removed-samples-rv",
        factors = "factors-rv",
        groups = "groups-rv",
        contrasts = "contrasts-rv",
        session = "session-token",
        applyResetStateFn = function(
            scope,
            initialState,
            designMatrix,
            dataClnReactive,
            sampleNamesReactive,
            removedSamples,
            factors,
            groups,
            contrasts,
            updateFormulaFn
        ) {
            apply_calls[[length(apply_calls) + 1]] <<- list(
                scope = scope,
                initialState = initialState,
                designMatrix = designMatrix,
                dataClnReactive = dataClnReactive,
                sampleNamesReactive = sampleNamesReactive,
                removedSamples = removedSamples,
                factors = factors,
                groups = groups,
                contrasts = contrasts
            )
            updateFormulaFn(initialState$formula)
            invisible(NULL)
        },
        removeModalFn = function() {
            remove_modal_calls <<- remove_modal_calls + 1L
            invisible(NULL)
        },
        showNotificationFn = function(message, type) {
            notifications[[length(notifications) + 1]] <<- list(message = message, type = type)
            invisible(NULL)
        },
        updateTextInputFn = function(session, inputId, value) {
            formula_updates[[length(formula_updates) + 1]] <<- list(
                session = session,
                inputId = inputId,
                value = value
            )
            invisible(NULL)
        }
    )

    expect_null(result)
    expect_identical(
        apply_calls,
        list(list(
            scope = "formula",
            initialState = list(formula = "~ group"),
            designMatrix = "design-matrix-rv",
            dataClnReactive = "data-cln-rv",
            sampleNamesReactive = "sample-names-rv",
            removedSamples = "removed-samples-rv",
            factors = "factors-rv",
            groups = "groups-rv",
            contrasts = "contrasts-rv"
        ))
    )
    expect_identical(
        formula_updates,
        list(list(session = "session-token", inputId = "formula_string", value = "~ group"))
    )
    expect_identical(remove_modal_calls, 1L)
    expect_identical(
        notifications,
        list(list(message = "Reset of formula completed.", type = "message"))
    )
})

test_that("registerLipidDesignResetConfirmationShells keeps confirm-reset observer handoff stable", {
    observed_labels <- character()
    shell_calls <- list()

    result <- registerLipidDesignResetConfirmationShells(
        input = list(confirm_reset = 1L, reset_scope = "contrasts"),
        initialState = function() list(marker = "initial-state"),
        designMatrix = "design-matrix-rv",
        dataClnReactive = "data-cln-rv",
        sampleNamesReactive = "sample-names-rv",
        removedSamples = "removed-samples-rv",
        factors = "factors-rv",
        groups = "groups-rv",
        contrasts = "contrasts-rv",
        session = "session-token",
        observeEventFn = function(eventExpr, handlerExpr) {
            observed_labels <<- c(
                observed_labels,
                paste(deparse(substitute(eventExpr)), collapse = "")
            )
            eval(substitute(handlerExpr), envir = parent.frame())
            invisible(NULL)
        },
        runResetConfirmationShellFn = function(
            scope,
            initialState,
            designMatrix,
            dataClnReactive,
            sampleNamesReactive,
            removedSamples,
            factors,
            groups,
            contrasts,
            session
        ) {
            shell_calls[[length(shell_calls) + 1]] <<- list(
                scope = scope,
                initialState = initialState,
                designMatrix = designMatrix,
                dataClnReactive = dataClnReactive,
                sampleNamesReactive = sampleNamesReactive,
                removedSamples = removedSamples,
                factors = factors,
                groups = groups,
                contrasts = contrasts,
                session = session
            )
            invisible(NULL)
        }
    )

    expect_null(result)
    expect_identical(observed_labels, "input$confirm_reset")
    expect_identical(
        shell_calls,
        list(list(
            scope = "contrasts",
            initialState = list(marker = "initial-state"),
            designMatrix = "design-matrix-rv",
            dataClnReactive = "data-cln-rv",
            sampleNamesReactive = "sample-names-rv",
            removedSamples = "removed-samples-rv",
            factors = "factors-rv",
            groups = "groups-rv",
            contrasts = "contrasts-rv",
            session = "session-token"
        ))
    )
})

test_that("runLipidDesignSampleSelectionInputShell keeps sample input refresh stable", {
    update_calls <- list()

    result <- runLipidDesignSampleSelectionInputShell(
        input = list(
            sample_to_rename = "Sample_10",
            selected_runs = c("Sample_2", "Sample_10"),
            samples_to_transform = c("Sample_3", "Missing"),
            tech_rep_samples = c("Sample_1", "Sample_10"),
            samples_to_remove = c("Sample_10", "Sample_2")
        ),
        session = "session-token",
        designMatrix = function() {
            data.frame(
                Run = c("Sample_10", "Sample_2", "Sample_1", "Sample_3"),
                stringsAsFactors = FALSE
            )
        },
        removedSamples = function() "Sample_10",
        reqFn = function(...) invisible(NULL),
        mixedsortFn = function(values) {
            expect_identical(values, c("Sample_10", "Sample_2", "Sample_1", "Sample_3"))
            c("Sample_1", "Sample_2", "Sample_3", "Sample_10")
        },
        isolateFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
        },
        updateSelectizeInputFn = function(session, inputId, choices, selected) {
            update_calls[[length(update_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                choices = choices,
                selected = selected
            )
        }
    )

    expect_null(result)
    expect_identical(
        update_calls,
        list(
            list(
                session = "session-token",
                inputId = "sample_to_rename",
                choices = c("Sample_1", "Sample_2", "Sample_3"),
                selected = ""
            ),
            list(
                session = "session-token",
                inputId = "selected_runs",
                choices = c("Sample_1", "Sample_2", "Sample_3"),
                selected = "Sample_2"
            ),
            list(
                session = "session-token",
                inputId = "samples_to_transform",
                choices = c("Sample_1", "Sample_2", "Sample_3"),
                selected = "Sample_3"
            ),
            list(
                session = "session-token",
                inputId = "tech_rep_samples",
                choices = c("Sample_1", "Sample_2", "Sample_3"),
                selected = "Sample_1"
            ),
            list(
                session = "session-token",
                inputId = "samples_to_remove",
                choices = c("Sample_1", "Sample_2", "Sample_3"),
                selected = "Sample_2"
            )
        )
    )
})

test_that("registerLipidDesignSampleSelectionInputShells keeps observe handoff stable", {
    observe_calls <- 0L
    shell_calls <- list()

    result <- registerLipidDesignSampleSelectionInputShells(
        input = "input-token",
        session = "session-token",
        designMatrix = "design-matrix-rv",
        removedSamples = "removed-samples-rv",
        observeFn = function(expr) {
            observe_calls <<- observe_calls + 1L
            eval(substitute(expr), envir = parent.frame())
            invisible(NULL)
        },
        runSampleSelectionInputShellFn = function(input, session, designMatrix, removedSamples) {
            shell_calls[[length(shell_calls) + 1]] <<- list(
                input = input,
                session = session,
                designMatrix = designMatrix,
                removedSamples = removedSamples
            )
            invisible(NULL)
        }
    )

    expect_null(result)
    expect_identical(observe_calls, 1L)
    expect_identical(
        shell_calls,
        list(list(
            input = "input-token",
            session = "session-token",
            designMatrix = "design-matrix-rv",
            removedSamples = "removed-samples-rv"
        ))
    )
})

test_that("runLipidDesignFactorGroupDropdownShell keeps dropdown update handoff stable", {
    update_calls <- list()

    result <- runLipidDesignFactorGroupDropdownShell(
        input = list(
            factor1_select = "Treatment",
            factor2_select = "Day1",
            factor3_select = "",
            contrast_group1 = "Case",
            contrast_group2 = "Control"
        ),
        session = "session-token",
        factors = function() c("Treatment", "Day1"),
        groups = function() c("Case", "Control"),
        isolateFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
        },
        updateSelectInputFn = function(session, inputId, choices, selected) {
            update_calls[[length(update_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                choices = choices,
                selected = selected
            )
            invisible(NULL)
        }
    )

    expect_null(result)
    expect_identical(
        update_calls,
        list(
            list(
                session = "session-token",
                inputId = "factor1_select",
                choices = c("", "Treatment", "Day1"),
                selected = "Treatment"
            ),
            list(
                session = "session-token",
                inputId = "factor2_select",
                choices = c("", "Treatment", "Day1"),
                selected = "Day1"
            ),
            list(
                session = "session-token",
                inputId = "factor3_select",
                choices = c("", "Treatment", "Day1"),
                selected = ""
            ),
            list(
                session = "session-token",
                inputId = "contrast_group1",
                choices = c("", "Case", "Control"),
                selected = "Case"
            ),
            list(
                session = "session-token",
                inputId = "contrast_group2",
                choices = c("", "Case", "Control"),
                selected = "Control"
            )
        )
    )
})

test_that("registerLipidDesignFactorGroupDropdownShells keeps observe handoff stable", {
    observe_calls <- 0L
    shell_calls <- list()

    result <- registerLipidDesignFactorGroupDropdownShells(
        input = "input-token",
        session = "session-token",
        factors = "factors-rv",
        groups = "groups-rv",
        observeFn = function(expr) {
            observe_calls <<- observe_calls + 1L
            eval(substitute(expr), envir = parent.frame())
            invisible(NULL)
        },
        runFactorGroupDropdownShellFn = function(input, session, factors, groups) {
            shell_calls[[length(shell_calls) + 1]] <<- list(
                input = input,
                session = session,
                factors = factors,
                groups = groups
            )
            invisible(NULL)
        }
    )

    expect_null(result)
    expect_identical(observe_calls, 1L)
    expect_identical(
        shell_calls,
        list(list(
            input = "input-token",
            session = "session-token",
            factors = "factors-rv",
            groups = "groups-rv"
        ))
    )
})

test_that("getLipidDesignSampleColumns keeps sample-column detection stable", {
    assay_list <- list(
        LCMS_Pos = data.frame(
            lipid_id = c("L1", "L2"),
            annotation = c("A1", "A2"),
            Sample_2 = c(10, 20),
            Average.Rt = c(1, 2),
            Sample_10 = c(30, 40),
            check.names = FALSE
        )
    )

    detected <- getLipidDesignSampleColumns(
        assayList = assay_list,
        colMap = list(lipid_id_col = "lipid_id", annotation_col = "annotation"),
        messageFn = function(...) invisible(NULL)
    )
    fallback <- getLipidDesignSampleColumns(
        assayList = assay_list,
        colMap = NULL,
        messageFn = function(...) invisible(NULL)
    )

    expect_identical(detected, c("Sample_2", "Sample_10"))
    expect_identical(fallback, c("Sample_2", "Average.Rt", "Sample_10"))
})

test_that("buildLipidDesignInitialState keeps imported-state reuse stable", {
    assay_list <- list(
        LCMS_Pos = data.frame(
            Sample_2 = c(10, 20),
            Sample_10 = c(30, 40),
            check.names = FALSE
        )
    )
    existing_design <- data.frame(
        Run = c("Sample_2", "Sample_10"),
        group = c("Control", "Treatment"),
        factor1 = c("Control", "Treatment"),
        factor2 = c("Day1", "Day2"),
        stringsAsFactors = FALSE
    )
    existing_contrasts <- data.frame(
        contrasts = "groupTreatment-groupControl",
        stringsAsFactors = FALSE
    )

    result <- buildLipidDesignInitialState(
        assayList = assay_list,
        configList = list(deAnalysisParameters = list(formula_string = "~ group + batch")),
        colMap = list(sample_columns = c("Sample_2", "Sample_10")),
        existingDesignMatrix = existing_design,
        existingContrasts = existing_contrasts,
        getSampleColumnsFn = function(assayList, colMap) colMap$sample_columns,
        messageFn = function(...) invisible(NULL)
    )

    expect_identical(result$design_matrix, existing_design)
    expect_identical(result$data_cln, assay_list)
    expect_identical(result$sample_names, c("Sample_2", "Sample_10"))
    expect_identical(result$groups, c("Control", "Treatment"))
    expect_identical(result$factors, c("Control", "Treatment", "Day1", "Day2"))
    expect_identical(result$formula, "~ group + batch")
    expect_identical(result$contrasts$contrast_name, "Treatment.vs.Control")
    expect_identical(result$contrasts$numerator, "Treatment")
    expect_identical(result$contrasts$denominator, "Control")
})

test_that("buildLipidDesignInitialState keeps fresh-state defaults stable", {
    assay_list <- list(
        LCMS_Pos = data.frame(
            Sample_10 = c(10, 20),
            Sample_2 = c(30, 40),
            check.names = FALSE
        )
    )

    result <- buildLipidDesignInitialState(
        assayList = assay_list,
        configList = list(deAnalysisParameters = list(formula_string = "~ 0 + group")),
        colMap = list(sample_columns = c("Sample_10", "Sample_2")),
        getSampleColumnsFn = function(assayList, colMap) colMap$sample_columns,
        mixedsortFn = function(samples) c("Sample_2", "Sample_10"),
        messageFn = function(...) invisible(NULL)
    )

    expect_identical(result$design_matrix$Run, c("Sample_2", "Sample_10"))
    expect_true(all(is.na(result$design_matrix$group)))
    expect_identical(result$data_cln, assay_list)
    expect_identical(result$sample_names, c("Sample_10", "Sample_2"))
    expect_identical(result$groups, character(0))
    expect_identical(result$factors, character(0))
    expect_identical(result$formula, "~ 0 + group")
    expect_identical(
        result$contrasts,
        data.frame(
            contrast_name = character(),
            numerator = character(),
            denominator = character(),
            stringsAsFactors = FALSE
        )
    )
})

test_that("runLipidDesignInitialStateShell keeps reset state handoff stable", {
    design_matrix_rv <- makeReactiveValMock(NULL)
    data_cln_rv <- makeReactiveValMock(NULL)
    sample_names_rv <- makeReactiveValMock(NULL)
    groups_rv <- makeReactiveValMock(NULL)
    factors_rv <- makeReactiveValMock(NULL)
    contrasts_rv <- makeReactiveValMock(NULL)
    removed_samples_rv <- makeReactiveValMock("stale")
    selectize_calls <- list()
    text_calls <- list()
    state <- list(
        design_matrix = data.frame(
            Run = c("Sample_2", "Sample_10"),
            stringsAsFactors = FALSE
        ),
        data_cln = list(LCMS_Pos = data.frame(Sample_2 = 1, check.names = FALSE)),
        sample_names = c("Sample_2", "Sample_10"),
        groups = c("Control", "Treatment"),
        factors = c("Treatment", "Day1"),
        contrasts = data.frame(
            contrast_name = "Treatment.vs.Control",
            numerator = "Treatment",
            denominator = "Control",
            stringsAsFactors = FALSE
        ),
        formula = "~ 0 + group"
    )

    result <- runLipidDesignInitialStateShell(
        state = state,
        session = "session-token",
        designMatrix = design_matrix_rv,
        dataClnReactive = data_cln_rv,
        sampleNamesReactive = sample_names_rv,
        groups = groups_rv,
        factors = factors_rv,
        contrasts = contrasts_rv,
        removedSamples = removed_samples_rv,
        updateSelectizeInputFn = function(session, inputId, choices, selected) {
            selectize_calls[[length(selectize_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                choices = choices,
                selected = selected
            )
            invisible(NULL)
        },
        updateTextInputFn = function(session, inputId, value) {
            text_calls[[length(text_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                value = value
            )
            invisible(NULL)
        }
    )

    expect_identical(result, state)
    expect_identical(design_matrix_rv(), state$design_matrix)
    expect_identical(data_cln_rv(), state$data_cln)
    expect_identical(sample_names_rv(), state$sample_names)
    expect_identical(groups_rv(), state$groups)
    expect_identical(factors_rv(), state$factors)
    expect_identical(contrasts_rv(), state$contrasts)
    expect_identical(removed_samples_rv(), character(0))
    expect_identical(
        vapply(selectize_calls, `[[`, character(1), "inputId"),
        c(
            "sample_to_rename",
            "selected_runs",
            "samples_to_transform",
            "tech_rep_samples",
            "samples_to_remove"
        )
    )
    expect_true(all(vapply(selectize_calls, function(call) identical(call$choices, state$design_matrix$Run), logical(1))))
    expect_true(all(vapply(selectize_calls, function(call) identical(call$selected, ""), logical(1))))
    expect_identical(
        text_calls,
        list(list(
            session = "session-token",
            inputId = "formula_string",
            value = "~ 0 + group"
        ))
    )
})

test_that("registerLipidDesignInitialStateShells keeps bind-event handoff stable", {
    observe_calls <- 0L
    bind_calls <- list()
    req_calls <- list()
    shell_calls <- list()

    result <- registerLipidDesignInitialStateShells(
        dataTbl = function() "data-token",
        input = list(reset_changes = 7L),
        initialState = function() "state-token",
        session = "session-token",
        designMatrix = "design-matrix-rv",
        dataClnReactive = "data-cln-rv",
        sampleNamesReactive = "sample-names-rv",
        groups = "groups-rv",
        factors = "factors-rv",
        contrasts = "contrasts-rv",
        removedSamples = "removed-samples-rv",
        observeFn = function(expr) {
            observe_calls <<- observe_calls + 1L
            eval(substitute(expr), envir = parent.frame())
            "observer-token"
        },
        bindEventFn = function(observer, ...) {
            bind_calls[[length(bind_calls) + 1]] <<- list(
                observer = observer,
                labels = vapply(
                    as.list(substitute(list(...)))[-1],
                    function(expr) paste(deparse(expr), collapse = ""),
                    character(1)
                )
            )
            observer
        },
        reqFn = function(...) {
            req_calls[[length(req_calls) + 1]] <<- list(...)
            invisible(NULL)
        },
        runInitialStateShellFn = function(
            state,
            session,
            designMatrix,
            dataClnReactive,
            sampleNamesReactive,
            groups,
            factors,
            contrasts,
            removedSamples
        ) {
            shell_calls[[length(shell_calls) + 1]] <<- list(
                state = state,
                session = session,
                designMatrix = designMatrix,
                dataClnReactive = dataClnReactive,
                sampleNamesReactive = sampleNamesReactive,
                groups = groups,
                factors = factors,
                contrasts = contrasts,
                removedSamples = removedSamples
            )
            invisible(NULL)
        }
    )

    expect_null(result)
    expect_identical(observe_calls, 1L)
    expect_identical(req_calls, list(list("state-token")))
    expect_identical(
        shell_calls,
        list(list(
            state = "state-token",
            session = "session-token",
            designMatrix = "design-matrix-rv",
            dataClnReactive = "data-cln-rv",
            sampleNamesReactive = "sample-names-rv",
            groups = "groups-rv",
            factors = "factors-rv",
            contrasts = "contrasts-rv",
            removedSamples = "removed-samples-rv"
        ))
    )
    expect_identical(
        bind_calls,
        list(list(
            observer = "observer-token",
            labels = c("dataTbl()", "input$reset_changes")
        ))
    )
})

test_that("buildLipidDesignActiveDataTable keeps removed-run filtering stable", {
    design_matrix <- data.frame(
        Run = c("Sample_1", "Sample_2", "Sample_3"),
        group = c("Control", "Treatment", "Control"),
        stringsAsFactors = FALSE
    )

    result <- buildLipidDesignActiveDataTable(
        designMatrix = design_matrix,
        currentRemovedSamples = c("Sample_2", "Missing")
    )

    expect_identical(result$Run, c("Sample_1", "Sample_3"))
    expect_identical(result$group, c("Control", "Control"))
})

test_that("registerLipidDesignDataTableShells keeps render and refresh handoff stable", {
    output <- new.env(parent = emptyenv())
    render_calls <- list()
    observe_calls <- 0L
    replace_calls <- list()

    result <- registerLipidDesignDataTableShells(
        output = output,
        proxyDataTable = "proxy-token",
        designMatrix = function() {
            data.frame(
                Run = c("Sample_1", "Sample_2", "Sample_3"),
                group = c("Control", "Treatment", "Control"),
                stringsAsFactors = FALSE
            )
        },
        removedSamples = function() "Sample_2",
        renderDtFn = function(expr, selection, options) {
            render_calls[[length(render_calls) + 1]] <<- list(
                data = force(expr),
                selection = selection,
                options = options
            )
            structure(list(kind = "mock_dt"), class = "mock_dt")
        },
        observeFn = function(expr) {
            observe_calls <<- observe_calls + 1L
            eval(substitute(expr), envir = parent.frame())
            invisible(NULL)
        },
        reqFn = function(...) invisible(NULL),
        replaceDataFn = function(proxy, data, resetPaging) {
            replace_calls[[length(replace_calls) + 1]] <<- list(
                proxy = proxy,
                data = data,
                resetPaging = resetPaging
            )
            invisible(NULL)
        }
    )

    expect_null(result)
    expect_identical(length(render_calls), 1L)
    expect_identical(render_calls[[1]]$data$Run, c("Sample_1", "Sample_3"))
    expect_identical(render_calls[[1]]$selection, "none")
    expect_identical(
        render_calls[[1]]$options,
        list(pageLength = 10, scrollX = TRUE, server = FALSE)
    )
    expect_identical(observe_calls, 1L)
    expect_identical(
        replace_calls,
        list(list(
            proxy = "proxy-token",
            data = data.frame(
                Run = c("Sample_1", "Sample_3"),
                group = c("Control", "Control"),
                stringsAsFactors = FALSE
            ),
            resetPaging = FALSE
        ))
    )
    expect_s3_class(output$data_table, "mock_dt")
})

test_that("buildLipidDesignSaveResultsDesignMatrix keeps assignment filtering stable", {
    design_matrix <- data.frame(
        Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
        group = c("Control", "", "Treatment", NA_character_),
        replicates = c(1L, 2L, 3L, 4L),
        stringsAsFactors = FALSE
    )

    result <- buildLipidDesignSaveResultsDesignMatrix(
        designMatrix = design_matrix,
        currentRemovedSamples = "Sample_3"
    )

    expect_identical(result$Run, "Sample_1")
    expect_identical(result$group, "Control")
    expect_identical(result$tech_rep_group, "Control_1")
})

test_that("buildLipidDesignSaveResultsDataList keeps assay column selection stable", {
    assay_list <- list(
        LCMS_Pos = data.frame(
            lipid_id = c("L1", "L2"),
            annotation = c("A1", "A2"),
            Sample_1 = c(10, 20),
            Sample_3 = c(30, 40),
            Other = c("x", "y"),
            check.names = FALSE
        )
    )

    result <- buildLipidDesignSaveResultsDataList(
        dataCln = assay_list,
        assignedSamples = c("Sample_3", "Sample_2"),
        colMap = list(lipid_id_col = "lipid_id", annotation_col = "annotation")
    )

    expect_identical(
        names(result$LCMS_Pos),
        c("lipid_id", "annotation", "Sample_3")
    )
})

test_that("buildLipidDesignSaveResultsPayload keeps payload assembly stable", {
    design_matrix <- data.frame(
        Run = c("Sample_1", "Sample_2"),
        group = c("Control", "Treatment"),
        replicates = c(1L, 2L),
        stringsAsFactors = FALSE
    )
    assay_list <- list(
        LCMS_Pos = data.frame(
            lipid_id = c("L1", "L2"),
            annotation = c("A1", "A2"),
            Sample_1 = c(10, 20),
            Sample_2 = c(30, 40),
            check.names = FALSE
        )
    )
    contrast_data <- data.frame(
        contrast_name = "Treatment.vs.Control",
        numerator = "Treatment",
        denominator = "Control",
        stringsAsFactors = FALSE
    )
    config_list <- list(deAnalysisParameters = list(formula_string = "~ group"))

    result <- buildLipidDesignSaveResultsPayload(
        designMatrix = design_matrix,
        currentRemovedSamples = "Sample_2",
        dataCln = assay_list,
        colMap = list(lipid_id_col = "lipid_id", annotation_col = "annotation"),
        contrastData = contrast_data,
        configList = config_list,
        formulaString = "~ 0 + group"
    )

    expect_identical(result$design_matrix$Run, "Sample_1")
    expect_identical(result$data_cln$LCMS_Pos$Sample_1, c(10, 20))
    expect_identical(result$contrasts_tbl$full_format, "Treatment_vs_Control=groupTreatment-groupControl")
    expect_identical(
        result$config_list$deAnalysisParameters$formula_string,
        "~ 0 + group"
    )
})

test_that("runLipidDesignSaveResultsShell keeps payload and notification handoff stable", {
    result_calls <- list()
    notifications <- list()

    result <- runLipidDesignSaveResultsShell(
        designMatrix = "design-matrix",
        currentRemovedSamples = "removed-samples",
        dataCln = "data-list",
        colMap = "column-map",
        contrastData = "contrast-data",
        configList = "config-list",
        formulaString = "~ 0 + group",
        resultSetter = function(value) {
            result_calls[[length(result_calls) + 1]] <<- value
            invisible(NULL)
        },
        buildSaveResultsPayloadFn = function(
            designMatrix,
            currentRemovedSamples,
            dataCln,
            colMap,
            contrastData,
            configList,
            formulaString
        ) {
            expect_identical(designMatrix, "design-matrix")
            expect_identical(currentRemovedSamples, "removed-samples")
            expect_identical(dataCln, "data-list")
            expect_identical(colMap, "column-map")
            expect_identical(contrastData, "contrast-data")
            expect_identical(configList, "config-list")
            expect_identical(formulaString, "~ 0 + group")

            list(marker = "final-result")
        },
        showNotificationFn = function(message, type, duration = NULL) {
            notifications[[length(notifications) + 1]] <<- list(
                message = message,
                type = type,
                duration = duration
            )
            invisible(NULL)
        }
    )

    expect_identical(result, list(marker = "final-result"))
    expect_identical(result_calls, list(list(marker = "final-result")))
    expect_identical(
        notifications,
        list(list(
            message = "Design saved successfully!",
            type = "message",
            duration = 5
        ))
    )
})

test_that("registerLipidDesignSaveResultsShells keeps save observer handoff stable", {
    observed_labels <- character()
    shell_calls <- list()

    result <- registerLipidDesignSaveResultsShells(
        input = list(save_results = 1L, formula_string = "~ group + batch"),
        designMatrix = function() "design-matrix",
        removedSamples = function() "removed-samples",
        dataClnReactive = function() "data-list",
        columnMapping = function() "column-map",
        contrasts = function() "contrast-data",
        configList = function() "config-list",
        resultRv = "result-rv",
        observeEventFn = function(eventExpr, handlerExpr) {
            observed_labels <<- c(
                observed_labels,
                paste(deparse(substitute(eventExpr)), collapse = "")
            )
            eval(substitute(handlerExpr), envir = parent.frame())
            invisible(NULL)
        },
        runSaveResultsShellFn = function(
            designMatrix,
            currentRemovedSamples,
            dataCln,
            colMap,
            contrastData,
            configList,
            formulaString,
            resultSetter
        ) {
            shell_calls[[length(shell_calls) + 1]] <<- list(
                designMatrix = designMatrix,
                currentRemovedSamples = currentRemovedSamples,
                dataCln = dataCln,
                colMap = colMap,
                contrastData = contrastData,
                configList = configList,
                formulaString = formulaString,
                resultSetter = resultSetter
            )
            invisible(NULL)
        }
    )

    expect_null(result)
    expect_identical(observed_labels, "input$save_results")
    expect_identical(
        shell_calls,
        list(list(
            designMatrix = "design-matrix",
            currentRemovedSamples = "removed-samples",
            dataCln = "data-list",
            colMap = "column-map",
            contrastData = "contrast-data",
            configList = "config-list",
            formulaString = "~ group + batch",
            resultSetter = "result-rv"
        ))
    )
})
