# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.WORKFLOW_FORMAT_ALLOWED_STATUSES <- c(
    "scientifically_supported",
    "reader_characterized"
)

newWorkflowFormatSupport <- function(
    formatId,
    omicType,
    uiValue,
    uiLabel,
    supportStatus
) {
    list(
        format_id = formatId,
        omic_type = omicType,
        ui_value = uiValue,
        ui_label = uiLabel,
        support_status = supportStatus
    )
}

.WORKFLOW_FORMAT_SUPPORT_CATALOGUE <- list(
    newWorkflowFormatSupport(
        "proteomics.diann",
        "proteomics",
        "diann",
        "DIA-NN",
        "scientifically_supported"
    ),
    newWorkflowFormatSupport(
        "proteomics.spectronaut",
        "proteomics",
        "spectronaut",
        "Spectronaut DIA (unverified)",
        "advertised_unverified"
    ),
    newWorkflowFormatSupport(
        "proteomics.fragpipe",
        "proteomics",
        "fragpipe",
        "FragPipe LFQ",
        "scientifically_supported"
    ),
    newWorkflowFormatSupport(
        "proteomics.maxquant",
        "proteomics",
        "maxquant",
        "MaxQuant LFQ",
        "scientifically_supported"
    ),
    newWorkflowFormatSupport(
        "proteomics.pd_tmt",
        "proteomics",
        "pd_tmt",
        "Proteome Discoverer TMT",
        "scientifically_supported"
    ),
    newWorkflowFormatSupport(
        "metabolomics.msdial",
        "metabolomics",
        "msdial",
        "MS-DIAL",
        "reader_characterized"
    ),
    newWorkflowFormatSupport(
        "metabolomics.progenesis",
        "metabolomics",
        "progenesis",
        "Progenesis QI (unavailable)",
        "detection_only"
    ),
    newWorkflowFormatSupport(
        "metabolomics.xcms",
        "metabolomics",
        "xcms",
        "XCMS (unavailable)",
        "detection_only"
    ),
    newWorkflowFormatSupport(
        "metabolomics.compound_discoverer",
        "metabolomics",
        "compound_discoverer",
        "Compound Discoverer (unavailable)",
        "detection_only"
    ),
    newWorkflowFormatSupport(
        "metabolomics.custom",
        "metabolomics",
        "custom",
        "Other/Custom",
        "scientifically_supported"
    ),
    newWorkflowFormatSupport(
        "lipidomics.msdial",
        "lipidomics",
        "msdial",
        "MS-DIAL",
        "reader_characterized"
    ),
    newWorkflowFormatSupport(
        "lipidomics.progenesis",
        "lipidomics",
        "progenesis",
        "Progenesis QI (unavailable)",
        "detection_only"
    ),
    newWorkflowFormatSupport(
        "lipidomics.xcms",
        "lipidomics",
        "xcms",
        "XCMS (unavailable)",
        "detection_only"
    ),
    newWorkflowFormatSupport(
        "lipidomics.compound_discoverer",
        "lipidomics",
        "compound_discoverer",
        "Compound Discoverer (unavailable)",
        "detection_only"
    ),
    newWorkflowFormatSupport(
        "lipidomics.lipidsearch",
        "lipidomics",
        "lipidsearch",
        "LipidSearch",
        "scientifically_supported"
    ),
    newWorkflowFormatSupport(
        "lipidomics.custom",
        "lipidomics",
        "custom",
        "Other/Custom",
        "reader_characterized"
    )
)

workflowFormatSupportCatalogue <- function() {
    .WORKFLOW_FORMAT_SUPPORT_CATALOGUE
}

workflowFormatSupportAbort <- function(
    message,
    class,
    omicType,
    requestedFormat,
    detectedFormat = NULL,
    supportStatus = NULL
) {
    rlang::abort(
        message,
        class = c(
            class,
            "multischolar_format_support_error",
            "multischolar_workflow_error"
        ),
        omic_type = omicType,
        requested_format = requestedFormat,
        detected_format = detectedFormat,
        support_status = supportStatus
    )
}

workflowFormatSupportScalar <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

findWorkflowFormatSupport <- function(
    omicType,
    inputFormat,
    catalogue = workflowFormatSupportCatalogue()
) {
    matches <- vapply(catalogue, function(entry) {
        identical(entry$omic_type, omicType) &&
            identical(entry$ui_value, inputFormat)
    }, logical(1))
    if (sum(matches) > 1L) {
        workflowFormatSupportAbort(
            "format support catalogue contains an ambiguous format key",
            "multischolar_format_ambiguous",
            omicType,
            inputFormat
        )
    }
    if (!any(matches)) NULL else catalogue[[which(matches)]]
}

workflowFormatSupportStatusClass <- function(status) {
    switch(
        status,
        advertised_unverified = "multischolar_format_unverified",
        detection_only = "multischolar_format_detection_only",
        unsupported = "multischolar_format_unsupported",
        "multischolar_format_unknown"
    )
}

workflowFormatDetectionThreshold <- function(omicType) {
    switch(omicType, proteomics = 0.3, metabolomics = 0.2, lipidomics = 0.2, 1)
}

resolveWorkflowFormatSupport <- function(
    omicType,
    requestedFormat,
    detectedFormat = NULL,
    detectionConfidence = NULL,
    catalogue = workflowFormatSupportCatalogue()
) {
    if (!workflowFormatSupportScalar(omicType) ||
        !workflowFormatSupportScalar(requestedFormat)) {
        workflowFormatSupportAbort(
            "omic and requested format must be non-empty scalar strings",
            "multischolar_format_unknown",
            as.character(omicType %||% "unknown"),
            as.character(requestedFormat %||% "unknown"),
            detectedFormat
        )
    }

    if (!identical(requestedFormat, "auto")) {
        requestedEntry <- findWorkflowFormatSupport(
            omicType,
            requestedFormat,
            catalogue
        )
        if (is.null(requestedEntry)) {
            workflowFormatSupportAbort(
                sprintf("Unknown %s import format: %s", omicType, requestedFormat),
                "multischolar_format_unknown",
                omicType,
                requestedFormat,
                detectedFormat
            )
        }
        if (!requestedEntry$support_status %in% .WORKFLOW_FORMAT_ALLOWED_STATUSES) {
            workflowFormatSupportAbort(
                sprintf(
                    "%s import format '%s' is %s and cannot be imported",
                    omicType,
                    requestedFormat,
                    gsub("_", " ", requestedEntry$support_status)
                ),
                workflowFormatSupportStatusClass(requestedEntry$support_status),
                omicType,
                requestedFormat,
                detectedFormat,
                requestedEntry$support_status
            )
        }
        if (identical(requestedFormat, "custom")) {
            return(c(
                requestedEntry,
                list(
                    format = "custom",
                    detected_format = detectedFormat,
                    detection_confidence = detectionConfidence
                )
            ))
        }
    }

    if (!workflowFormatSupportScalar(detectedFormat) ||
        identical(detectedFormat, "unknown")) {
        workflowFormatSupportAbort(
            sprintf("Could not identify a supported %s import format", omicType),
            "multischolar_format_unknown",
            omicType,
            requestedFormat,
            detectedFormat
        )
    }
    if (!is.null(detectionConfidence)) {
        confidence <- suppressWarnings(as.numeric(detectionConfidence))
        validConfidence <- length(confidence) == 1L && is.finite(confidence) &&
            confidence >= workflowFormatDetectionThreshold(omicType) &&
            confidence <= 1
        if (!validConfidence) {
            workflowFormatSupportAbort(
                sprintf("%s format detection is ambiguous or low confidence", omicType),
                "multischolar_format_ambiguous",
                omicType,
                requestedFormat,
                detectedFormat
            )
        }
    }

    if (!identical(requestedFormat, "auto") &&
        !identical(requestedFormat, detectedFormat)) {
        workflowFormatSupportAbort(
            sprintf(
                "Selected %s format '%s' does not match detected format '%s'",
                omicType,
                requestedFormat,
                detectedFormat
            ),
            "multischolar_format_mismatch",
            omicType,
            requestedFormat,
            detectedFormat
        )
    }

    activeFormat <- if (identical(requestedFormat, "auto")) {
        detectedFormat
    } else {
        requestedFormat
    }
    entry <- findWorkflowFormatSupport(omicType, activeFormat, catalogue)
    if (is.null(entry)) {
        workflowFormatSupportAbort(
            sprintf("Unknown %s import format: %s", omicType, activeFormat),
            "multischolar_format_unknown",
            omicType,
            requestedFormat,
            detectedFormat
        )
    }
    if (!entry$support_status %in% .WORKFLOW_FORMAT_ALLOWED_STATUSES) {
        workflowFormatSupportAbort(
            sprintf(
                "%s import format '%s' is %s and cannot be imported",
                omicType,
                activeFormat,
                gsub("_", " ", entry$support_status)
            ),
            workflowFormatSupportStatusClass(entry$support_status),
            omicType,
            requestedFormat,
            detectedFormat,
            entry$support_status
        )
    }

    c(
        entry,
        list(
            format = activeFormat,
            detected_format = detectedFormat,
            detection_confidence = detectionConfidence
        )
    )
}

workflowFormatUiChoices <- function(omicType, includeAuto = FALSE) {
    entries <- Filter(
        function(entry) identical(entry$omic_type, omicType),
        workflowFormatSupportCatalogue()
    )
    choices <- stats::setNames(
        vapply(entries, `[[`, character(1), "ui_value"),
        vapply(entries, `[[`, character(1), "ui_label")
    )
    if (isTRUE(includeAuto)) {
        choices <- c("Auto-detect" = "auto", choices)
    }
    choices
}
