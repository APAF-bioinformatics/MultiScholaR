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
    designMatrix
}

#' The shared technical-replicate edit (`applyOmicsDesignTechRepAssignment()`) with the proteomics
#' builder's notification texts.
#' @return A list: `ok`, `type`, `message` and `designMatrix`.
#' @noRd
applyProtDesignTechRepAssignment <- function(designMatrix, techRepSamples, assignmentMode,
                                             manualReplicateNumber = NULL) {
    assignment <- applyOmicsDesignTechRepAssignment(designMatrix, techRepSamples, assignmentMode,
        manualReplicateNumber)
    if (length(techRepSamples) < 2) {
        assignment$message <- "Please select at least two samples to assign as technical replicates."
    } else if (!isTRUE(assignment$ok)) {
        assignment$message <- "All selected samples must belong to the same group. Please assign metadata first."
    } else {
        # The biological replicate number the selected runs now share.
        baseReplicateNumber <- assignment$designMatrix$replicates[assignment$designMatrix$Run == techRepSamples[1]]
        assignment$message <- paste("Assigned", length(techRepSamples),
            "samples as technical replicates with biological replicate number", baseReplicateNumber,
            "and adjusted subsequent replicate numbers")
    }
    assignment
}
