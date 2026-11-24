# KIT-103 W6: the proteomics builder's two edits that differ from the shared builder's (C4), for the
# `dia_long` adapter (mod_omic_design_builder.R). Its other edit observers are the shared
# builder's. Both keep today's proteomics behaviour (this file's observers at de2d891c).

#' Assign factor levels, group and replicate numbers to the selected runs, as the proteomics
#' builder does: the replicate numbers are integers (`seq.int`), where the shared C4 edit's are
#' doubles (`seq`); otherwise `applyOmicsDesignMetadataAssignment()`.
#' @return The design matrix.
#' @noRd
applyProtDesignMetadataAssignment <- function(designMatrix, selectedRuns, factor1Selection, factor2Selection = NULL,
                                              factor3Selection = NULL, replicateStart = NULL) {
    replicateNumbers <- if (!is.null(replicateStart) && !is.na(replicateStart)) {
        seq.int(replicateStart, length.out = length(selectedRuns))
    } else {
        NA_integer_
    }
    selected <- which(designMatrix$Run %in% selectedRuns)
    designMatrix$factor1[selected] <- factor1Selection
    designMatrix$factor2[selected] <- factor2Selection
    designMatrix$factor3[selected] <- if (!is.null(factor3Selection) && factor3Selection != "") {
        factor3Selection
    } else {
        NA_character_
    }
    designMatrix$replicates[selected] <- replicateNumbers
    designMatrix$group[selected] <- buildOmicsDesignGroupName(factor1Selection, factor2Selection, factor3Selection)
    # KIT-222 (OD-1): a new assignment ends the runs' technical-replicate marking.
    if ("tech_reps" %in% names(designMatrix)) designMatrix$tech_reps[selected] <- NA
    designMatrix
}

#' The shared technical-replicate edit (`applyOmicsDesignTechRepAssignment()`) with the proteomics
#' builder's notification texts (KIT-222: the shared `number_held` and `partial_sample` refusals keep the shared text).
#' @return A list: `ok`, `type`, `message` and `designMatrix` (and `code` on a refusal).
#' @noRd
applyProtDesignTechRepAssignment <- function(designMatrix, techRepSamples, assignmentMode,
                                             manualReplicateNumber = NULL) {
    integerNumbers <- is.integer(designMatrix$replicates)
    assignment <- applyOmicsDesignTechRepAssignment(designMatrix, techRepSamples, assignmentMode,
        manualReplicateNumber)
    if (length(techRepSamples) < 2) {
        assignment$message <- "Please select at least two samples to assign as technical replicates."
    } else if (identical(assignment$code, "groups")) {
        assignment$message <- "All selected samples must belong to the same group. Please assign metadata first."
    } else if (isTRUE(assignment$ok)) {
        # KIT-222: proteomics' replicate numbers stay integers (a typed number is a double).
        if (integerNumbers) assignment$designMatrix$replicates <- as.integer(assignment$designMatrix$replicates)
        # The biological replicate number the selected runs now share.
        baseReplicateNumber <- assignment$designMatrix$replicates[assignment$designMatrix$Run == techRepSamples[1]]
        assignment$message <- paste("Assigned", length(techRepSamples),
            "samples as technical replicates with biological replicate number", baseReplicateNumber,
            "and adjusted subsequent replicate numbers")
    }
    # KIT-222: the shared refusals (`number_held`, `partial_sample`) keep their own text.
    assignment
}
