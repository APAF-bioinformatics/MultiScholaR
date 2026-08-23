.formatSupportInventory <- function() {
    jsonlite::read_json(
        testthat::test_path("..", "testdata", "omics-capabilities.json"),
        simplifyVector = FALSE
    )
}

.formatSupportReaderResult <- function(idColumn) {
    data <- data.frame(
        feature_id = "F1",
        sample_a = 1,
        check.names = FALSE
    )
    names(data)[[1L]] <- idColumn
    list(
        data = data,
        detected_columns = list(
            metabolite_id = idColumn,
            lipid_id = idColumn,
            annotation = ""
        ),
        sample_columns = "sample_a",
        is_pattern = NA_character_
    )
}

test_that("runtime format support matches the versioned inventory", {
    inventory <- .formatSupportInventory()
    runtime <- workflowFormatSupportCatalogue()
    inventoryRows <- lapply(inventory$formats, function(entry) {
        list(
            format_id = entry$format_id,
            omic_type = entry$omic_type,
            ui_value = entry$ui_value,
            support_status = entry$support_status
        )
    })
    runtimeRows <- lapply(runtime, function(entry) {
        entry[c("format_id", "omic_type", "ui_value", "support_status")]
    })

    expect_identical(runtimeRows, inventoryRows)
    expect_identical(anyDuplicated(vapply(
        runtime,
        `[[`,
        character(1),
        "format_id"
    )), 0L)
    expect_match(inventory$inventory_version, "^2026-08-22[.]")
})

test_that("format support seams preserve existing positional arguments", {
    expect_identical(
        head(names(formals(importProtImportDataByFormat)), 9L),
        c(
            "format", "searchResultsPath", "input", "importDiann",
            "importSpectronaut", "importFragpipe", "importMaxquant",
            "importPdTmt", "logError"
        )
    )
    expect_identical(
        tail(names(formals(registerProtImportObservers)), 2L),
        c("messageFn", "resolveFormatSupportFn")
    )
    expect_identical(
        head(names(formals(buildMetabImportWorkflowPayload)), 16L),
        c(
            "assay1Name", "assay1Data", "assay1File", "assay2File",
            "assay2Name", "vendorFormat", "detectedFormat", "metaboliteCol",
            "annotationCol", "sampleCols", "sanitizeNames", "isPattern",
            "assay2Importer", "cleanNamesFn", "mapAssaysFn", "timestampFn"
        )
    )
    expect_identical(
        deparse(formals(buildMetabImportWorkflowPayload)$assay2Importer),
        "importMetabMSDIALData"
    )
    expect_identical(
        head(names(formals(prepareMetabImportAssaySelectionState)), 5L),
        c(
            "assay1File", "detectFormatFn", "defaultImporter", "importers",
            "readHeadersFn"
        )
    )
    expect_identical(
        deparse(formals(prepareMetabImportAssaySelectionState)$defaultImporter),
        "importMetabMSDIALData"
    )
    expect_identical(
        head(names(formals(loadLipidImportAssayPreview)), 6L),
        c(
            "assay1File", "readHeaders", "detectFormat", "importMsdial",
            "importLipidSearch", "logInfo"
        )
    )
    expect_identical(
        head(names(formals(assembleLipidImportAssayList)), 10L),
        c(
            "assay1Name", "assay1Data", "assay2File", "assay2Name",
            "dataFormat", "lipidIdCol", "annotationCol", "importSecondAssay",
            "resolveSecondAssayReader", "callSecondAssayReader"
        )
    )
})

test_that("all advertised formats resolve or reject by audited support status", {
    for (entry in workflowFormatSupportCatalogue()) {
        detected <- if (identical(entry$ui_value, "custom")) {
            "msdial"
        } else {
            entry$ui_value
        }
        resolve <- function() {
            resolveWorkflowFormatSupport(
                omicType = entry$omic_type,
                requestedFormat = entry$ui_value,
                detectedFormat = detected,
                detectionConfidence = 1
            )
        }
        if (entry$support_status %in% .WORKFLOW_FORMAT_ALLOWED_STATUSES) {
            decision <- resolve()
            expect_identical(decision$format, entry$ui_value, info = entry$format_id)
            expect_identical(
                decision$support_status,
                entry$support_status,
                info = entry$format_id
            )
        } else {
            expectedClass <- workflowFormatSupportStatusClass(entry$support_status)
            expect_error(resolve(), class = expectedClass, info = entry$format_id)
        }
    }
})

test_that("format support distinguishes unknown, mismatch, and ambiguous input", {
    expect_error(
        resolveWorkflowFormatSupport("proteomics", "unknown_vendor", "diann", 1),
        class = "multischolar_format_unknown"
    )
    expect_error(
        resolveWorkflowFormatSupport("proteomics", "diann", "maxquant", 1),
        class = "multischolar_format_mismatch"
    )
    expect_error(
        resolveWorkflowFormatSupport("proteomics", "auto", "diann", 0.1),
        class = "multischolar_format_ambiguous"
    )
    expect_error(
        resolveWorkflowFormatSupport("lipidomics", "auto", "unknown", 0),
        class = "multischolar_format_unknown"
    )

    custom <- resolveWorkflowFormatSupport(
        "metabolomics",
        "custom",
        "progenesis",
        1
    )
    expect_identical(custom$format, "custom")
    expect_identical(custom$detected_format, "progenesis")
})

test_that("format support exposes unsupported and duplicate-catalogue classes", {
    unsupported <- newWorkflowFormatSupport(
        "proteomics.retired",
        "proteomics",
        "retired",
        "Retired format",
        "unsupported"
    )
    expect_error(
        resolveWorkflowFormatSupport(
            "proteomics",
            "retired",
            "retired",
            1,
            catalogue = c(workflowFormatSupportCatalogue(), list(unsupported))
        ),
        class = "multischolar_format_unsupported"
    )

    diann <- findWorkflowFormatSupport("proteomics", "diann")
    expect_error(
        resolveWorkflowFormatSupport(
            "proteomics",
            "diann",
            "diann",
            1,
            catalogue = c(workflowFormatSupportCatalogue(), list(diann))
        ),
        class = "multischolar_format_ambiguous"
    )
})

test_that("Spectronaut rejection occurs before its reader is invoked", {
    readerCalled <- FALSE
    expect_error(
        importProtImportDataByFormat(
            format = "spectronaut",
            searchResultsPath = "spectronaut.tsv",
            input = list(spectronaut_quantity = "pg"),
            importSpectronaut = function(...) {
                readerCalled <<- TRUE
            },
            logError = function(...) NULL
        ),
        class = "multischolar_format_unverified"
    )
    expect_false(readerCalled)
})

test_that("metabolomics preview rejects fallback and preserves explicit custom identity", {
    readerCalls <- character()
    detector <- function(headers, filename) {
        list(format = "progenesis", confidence = 0.9)
    }
    readers <- list(
        msdial = function(path) readerCalls <<- c(readerCalls, "msdial"),
        custom = function(path) {
            readerCalls <<- c(readerCalls, "custom")
            .formatSupportReaderResult("Feature")
        }
    )

    expect_error(
        prepareMetabImportAssaySelectionState(
            assay1File = "progenesis.tsv",
            detectFormatFn = detector,
            importers = readers,
            readHeadersFn = function(path) c("Compound", "Sample A"),
            vendorFormat = "progenesis"
        ),
        class = "multischolar_format_detection_only"
    )
    expect_length(readerCalls, 0L)

    expect_error(
        prepareMetabImportAssaySelectionState(
            assay1File = "mismatch.tsv",
            detectFormatFn = function(headers, filename) {
                list(format = "xcms", confidence = 0.9)
            },
            importers = readers,
            readHeadersFn = function(path) c("Feature", "sample_a"),
            vendorFormat = "msdial"
        ),
        class = "multischolar_format_mismatch"
    )
    expect_length(readerCalls, 0L)

    custom <- prepareMetabImportAssaySelectionState(
        assay1File = "custom.tsv",
        detectFormatFn = detector,
        importers = readers,
        readHeadersFn = function(path) c("Feature", "sample_a"),
        vendorFormat = "custom"
    )
    expect_identical(readerCalls, "custom")
    expect_identical(custom$formatInfo$format, "custom")
    expect_identical(custom$formatInfo$observed_format, "progenesis")
})

test_that("lipidomics preview rejects fallback and preserves explicit custom identity", {
    readerCalls <- character()
    detector <- function(headers, filename) {
        list(format = "xcms", confidence = 0.9)
    }
    importMsdial <- function(path) {
        readerCalls <<- c(readerCalls, "generic")
        .formatSupportReaderResult("Lipid")
    }

    expect_error(
        loadLipidImportAssayPreview(
            assay1File = "xcms.tsv",
            readHeaders = function(path) c("featureid", "sample_a"),
            detectFormat = detector,
            importMsdial = importMsdial,
            importLipidSearch = function(...) readerCalls <<- c(readerCalls, "lipidsearch"),
            logInfo = function(...) NULL,
            vendorFormat = "xcms"
        ),
        class = "multischolar_format_detection_only"
    )
    expect_length(readerCalls, 0L)

    expect_error(
        loadLipidImportAssayPreview(
            assay1File = "mismatch.tsv",
            readHeaders = function(path) c("Lipid", "sample_a"),
            detectFormat = function(headers, filename) {
                list(format = "xcms", confidence = 0.9)
            },
            importMsdial = importMsdial,
            importLipidSearch = function(...) {
                readerCalls <<- c(readerCalls, "lipidsearch")
            },
            logInfo = function(...) NULL,
            vendorFormat = "msdial"
        ),
        class = "multischolar_format_mismatch"
    )
    expect_length(readerCalls, 0L)

    custom <- loadLipidImportAssayPreview(
        assay1File = "custom.tsv",
        readHeaders = function(path) c("Lipid", "sample_a"),
        detectFormat = detector,
        importMsdial = importMsdial,
        importLipidSearch = function(...) readerCalls <<- c(readerCalls, "lipidsearch"),
        logInfo = function(...) NULL,
        vendorFormat = "custom"
    )
    expect_identical(readerCalls, "generic")
    expect_identical(custom$detectedFormat, "custom")
})

test_that("processing gates reject before state, files, or notifications change", {
    metabolomicsWorkflow <- new.env(parent = emptyenv())
    metabolomicsWorkflow$sentinel <- "unchanged"
    metabCalls <- character()
    expect_error(
        runMetabImportProcessing(
            assay1Data = data.frame(Feature = "F1", sample_a = 1),
            assay1Name = "LCMS_Pos",
            assay2File = NULL,
            assay2Name = "",
            vendorFormat = "progenesis",
            detectedFormat = "progenesis",
            sanitizeNames = FALSE,
            isPattern = "",
            getMetaboliteIdColFn = function() metabCalls <<- c(metabCalls, "id"),
            getAnnotationColFn = function() metabCalls <<- c(metabCalls, "annotation"),
            getSampleColumnsFn = function() metabCalls <<- c(metabCalls, "samples"),
            workflowData = metabolomicsWorkflow,
            showNotificationFn = function(...) metabCalls <<- c(metabCalls, "notify")
        ),
        class = "multischolar_format_detection_only"
    )
    expect_identical(metabCalls, character())
    expect_identical(as.list(metabolomicsWorkflow), list(sentinel = "unchanged"))

    lipidomicsWorkflow <- new.env(parent = emptyenv())
    lipidomicsWorkflow$sentinel <- "unchanged"
    lipidCalls <- character()
    expect_error(
        runLipidImportProcessing(
            workflowData = lipidomicsWorkflow,
            assay1Name = "LCMS_Pos",
            assay1Data = data.frame(Lipid = "L1", sample_a = 1),
            vendorFormat = "xcms",
            detectedFormat = "xcms",
            lipidIdCol = "Lipid",
            annotationCol = "",
            sampleColumns = "sample_a",
            isPattern = "",
            sanitizeNames = FALSE,
            notify = function(...) lipidCalls <<- c(lipidCalls, "notify")
        ),
        class = "multischolar_format_detection_only"
    )
    expect_identical(lipidCalls, character())
    expect_identical(as.list(lipidomicsWorkflow), list(sentinel = "unchanged"))
})
