proteomicsPublicationNegativeBase <- function(format) {
    if (identical(format, "diann")) {
        return(data.frame(
            Protein.Group = c("P1", "P2"),
            Stripped.Sequence = c("PEPTIDEA", "PEPTIDEB"),
            Run = c("S1", "S2"),
            Precursor.Quantity = c(100, 200),
            Precursor.Normalised = c(90, 180),
            check.names = FALSE,
            stringsAsFactors = FALSE
        ))
    }
    if (identical(format, "maxquant")) {
        return(data.frame(
            Protein.IDs = c("P1", "P2"),
            Potential.contaminant = c("", ""),
            Reverse = c("", ""),
            LFQ.intensity.S1 = c(100, 200),
            LFQ.intensity.S2 = c(110, 210),
            check.names = FALSE,
            stringsAsFactors = FALSE
        ))
    }
    if (identical(format, "fragpipe")) {
        return(data.frame(
            `Protein ID` = c("P1", "P2"),
            `S1 MaxLFQ Intensity` = c(100, 200),
            `S2 MaxLFQ Intensity` = c(110, 210),
            check.names = FALSE,
            stringsAsFactors = FALSE
        ))
    }
    if (identical(format, "pd_tmt")) {
        return(data.frame(
            Annotated.Sequence = c("PEPTIDEA", "PEPTIDEB"),
            Master.Protein.Accessions = c("P1", "P2"),
            `Abundance: F1: 126, S1` = c(100, 200),
            `Abundance: F1: 127N, S2` = c(110, 210),
            check.names = FALSE,
            stringsAsFactors = FALSE
        ))
    }
    proteomicsPublicationAbort("negative payload format is unsupported")
}

proteomicsPublicationWriteNegative <- function(data, path) {
    proteomicsPublicationWriteTable(data, path, append = FALSE)
    invisible(path)
}

proteomicsPublicationMutateDiann <- function(data, mutation, path) {
    if (identical(mutation, "missing_protein_group")) {
        data$Protein.Group <- NULL
    } else if (identical(mutation, "missing_quantity")) {
        data$Precursor.Quantity <- NULL
        data$Precursor.Normalised <- NULL
    } else if (identical(mutation, "nonrectangular")) {
        header <- paste(names(data), collapse = "\t")
        first <- paste(unlist(data[1L, ], use.names = FALSE), collapse = "\t")
        second <- paste(
            unlist(data[2L, -ncol(data)], use.names = FALSE),
            collapse = "\t"
        )
        writeLines(c(header, first, second), path, useBytes = TRUE)
        return(invisible(path))
    } else if (identical(mutation, "zero_and_nonfinite_quantity")) {
        data$Precursor.Quantity <- c(0, Inf)
        data$Precursor.Normalised <- c(-Inf, 0)
    }
    proteomicsPublicationWriteNegative(data, path)
}

proteomicsPublicationMutateMaxQuant <- function(data, mutation, path) {
    if (identical(mutation, "missing_protein_id")) {
        data$Protein.IDs <- NULL
    } else if (identical(mutation, "missing_intensity")) {
        data <- data[, !startsWith(names(data), "LFQ.intensity."), drop = FALSE]
    } else if (identical(mutation, "native_space_lfq_columns")) {
        names(data) <- sub("^LFQ\\.intensity\\.", "LFQ intensity ", names(data))
    } else if (identical(mutation, "duplicate_and_contaminant_rows")) {
        duplicate <- data[1L, , drop = FALSE]
        contaminant <- data[2L, , drop = FALSE]
        contaminant$Potential.contaminant <- "+"
        data <- rbind(data, duplicate, contaminant)
    }
    proteomicsPublicationWriteNegative(data, path)
}

proteomicsPublicationMutateFragPipe <- function(data, mutation, path) {
    if (identical(mutation, "missing_protein_id")) {
        data[["Protein ID"]] <- NULL
    } else if (identical(mutation, "missing_intensity")) {
        data <- data[, !grepl("Intensity$", names(data)), drop = FALSE]
    } else if (identical(mutation, "regular_intensity_without_maxlfq")) {
        names(data) <- sub(" MaxLFQ Intensity$", " Intensity", names(data))
    } else if (identical(mutation, "duplicate_protein_rows")) {
        data <- rbind(data, data[1L, , drop = FALSE])
    }
    proteomicsPublicationWriteNegative(data, path)
}

proteomicsPublicationMutatePdTmt <- function(data, mutation, path) {
    if (identical(mutation, "missing_accession")) {
        data$Master.Protein.Accessions <- NULL
    } else if (identical(mutation, "missing_abundance")) {
        data <- data[, !startsWith(names(data), "Abundance:"), drop = FALSE]
    } else if (identical(mutation, "duplicate_normalized_channel")) {
        names(data)[[4L]] <- "Abundance: F2: 126, S1"
    } else if (identical(mutation, "native_grouped_abundance_columns")) {
        names(data)[[2L]] <- "Master Protein Accessions"
        names(data)[3:4] <- c(
            "Abundances (Grouped): F1, 126",
            "Abundances (Grouped): F1, 127N"
        )
    }
    proteomicsPublicationWriteNegative(data, path)
}

proteomicsPublicationPrepareNegative <- function(case, path) {
    proteomicsPublicationValidateNegativeCase(case)
    if (file.exists(path) || dir.exists(path)) {
        proteomicsPublicationAbort("negative output already exists")
    }
    data <- proteomicsPublicationNegativeBase(case$format)
    switch(
        case$format,
        diann = proteomicsPublicationMutateDiann(data, case$mutation, path),
        maxquant = proteomicsPublicationMutateMaxQuant(data, case$mutation, path),
        fragpipe = proteomicsPublicationMutateFragPipe(data, case$mutation, path),
        pd_tmt = proteomicsPublicationMutatePdTmt(data, case$mutation, path)
    )
    list(
        path = path,
        sha256 = publicationFileDigest(path),
        size_bytes = as.numeric(file.info(path)$size)
    )
}

proteomicsPublicationEvaluateNegative <- function(case, path) {
    importer <- proteomicsPublicationImporter(case$format)
    observed <- tryCatch(
        suppressMessages(suppressWarnings(importer(path))),
        error = \(error) error
    )
    rejected <- inherits(observed, "error")
    expected_rejection <- identical(case$expected_outcome, "classed_rejection")
    if (!identical(rejected, expected_rejection)) {
        proteomicsPublicationAbort("negative outcome differs")
    }
    list(
        case_id = case$case_id,
        expected_outcome = case$expected_outcome,
        observed_outcome = if (rejected) {
            "classed_rejection"
        } else {
            "accepted_edge_characterization"
        },
        condition_class = if (rejected) class(observed)[[1L]] else NULL,
        imported_rows = if (rejected) NULL else as.numeric(nrow(observed$data)),
        performance_authority = FALSE,
        promotion_authority = FALSE
    )
}
