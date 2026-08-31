proteomicsPublicationTmtChannels <- function() {
    c(
        "126", "127N", "127C", "128N", "128C", "129N", "129C",
        "130N", "130C", "131N", "131C", "132N", "132C", "133N",
        "133C", "134N", "134C", "135N"
    )
}

proteomicsPublicationFormatSchema <- function(format, sample_count) {
    if (identical(format, "diann")) {
        required <- c(
            "Protein.Group", "Protein.Ids", "Genes", "Stripped.Sequence",
            "Modified.Sequence", "Precursor.Id", "Precursor.Charge", "Run",
            "Q.Value", "Global.Q.Value", "PG.Q.Value", "Global.PG.Q.Value",
            "Proteotypic", "Precursor.Quantity", "Precursor.Normalised"
        )
        orientation <- "long"
        rule <- "Precursor.Normalised"
    } else if (identical(format, "maxquant")) {
        required <- c(
            "Protein.IDs", "Gene.names", "Protein.names",
            "Potential.contaminant", "Reverse"
        )
        orientation <- "wide"
        rule <- "LFQ.intensity.<sample_id>"
    } else if (identical(format, "fragpipe")) {
        required <- c("Protein ID", "Gene")
        orientation <- "wide"
        rule <- "<sample_id> MaxLFQ Intensity"
    } else if (identical(format, "pd_tmt")) {
        required <- c("Annotated.Sequence", "Master.Protein.Accessions")
        orientation <- "wide_peptide_input"
        rule <- "Abundance: F<index>: <channel>, <sample_id>"
    } else {
        proteomicsPublicationAbort("serializer format is unsupported")
    }
    list(
        schema_id = paste(
            "multischolar.proteomics",
            format,
            "publication",
            sep = "."
        ),
        schema_version = "1.0.0",
        orientation = orientation,
        delimiter = "tab",
        required_columns = as.list(required),
        quantity_column_rule = rule,
        line_ending = "LF",
        encoding = "UTF-8",
        input_column_count = as.integer(length(required) + if (
            identical(format, "diann")
        ) 0L else sample_count)
    )
}

proteomicsPublicationBlockRanges <- function(entity_count, chunk_rows) {
    starts <- seq.int(1L, entity_count, by = chunk_rows)
    lapply(starts, \(start) {
        seq.int(start, min(entity_count, start + chunk_rows - 1L))
    })
}

proteomicsPublicationFormatNumeric <- function(value) {
    output <- rep.int("", length(value))
    finite <- is.finite(value)
    output[finite] <- formatC(
        value[finite],
        format = "f",
        digits = 6L,
        decimal.mark = "."
    )
    output[is.infinite(value) & value > 0] <- "Inf"
    output[is.infinite(value) & value < 0] <- "-Inf"
    output
}

proteomicsPublicationWriteTable <- function(data, path, append) {
    numeric <- vapply(data, is.double, logical(1))
    data[numeric] <- lapply(data[numeric], proteomicsPublicationFormatNumeric)
    utils::write.table(
        data,
        file = path,
        append = append,
        quote = FALSE,
        sep = "\t",
        eol = "\n",
        na = "",
        row.names = FALSE,
        col.names = !append,
        fileEncoding = "UTF-8"
    )
    invisible(path)
}

proteomicsPublicationDiannBlock <- function(plan, block) {
    sample_count <- nrow(plan$design)
    entity_count <- nrow(block$entities)
    protein_index <- rep(
        block$entities$protein_index,
        each = sample_count
    )
    peptide_id <- rep(block$entities$peptide_id, each = sample_count)
    precursor_id <- rep(block$entities$entity_id, each = sample_count)
    charge <- rep(block$entities$charge, each = sample_count)
    values <- as.vector(t(block$values))
    data.frame(
        Protein.Group = plan$proteins$protein_id[protein_index],
        Protein.Ids = plan$proteins$protein_id[protein_index],
        Genes = plan$proteins$gene_id[protein_index],
        Stripped.Sequence = peptide_id,
        Modified.Sequence = paste0("_", peptide_id, "_"),
        Precursor.Id = precursor_id,
        Precursor.Charge = charge,
        Run = rep(plan$design$sample_id, times = entity_count),
        Q.Value = 0.001,
        Global.Q.Value = 0.001,
        PG.Q.Value = 0.001,
        Global.PG.Q.Value = 0.001,
        Proteotypic = 1L,
        Precursor.Quantity = values,
        Precursor.Normalised = values,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

proteomicsPublicationMaxQuantBlock <- function(plan, block) {
    protein_index <- block$entities$protein_index
    annotations <- data.frame(
        Protein.IDs = plan$proteins$protein_id[protein_index],
        Gene.names = plan$proteins$gene_id[protein_index],
        Protein.names = paste0("Synthetic protein ", protein_index),
        Potential.contaminant = "",
        Reverse = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(block$values, check.names = FALSE)
    names(quantities) <- paste0("LFQ.intensity.", plan$design$sample_id)
    cbind(annotations, quantities)
}

proteomicsPublicationFragPipeBlock <- function(plan, block) {
    protein_index <- block$entities$protein_index
    annotations <- data.frame(
        `Protein ID` = plan$proteins$protein_id[protein_index],
        Gene = plan$proteins$gene_id[protein_index],
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(block$values, check.names = FALSE)
    names(quantities) <- paste(plan$design$sample_id, "MaxLFQ Intensity")
    cbind(annotations, quantities)
}

proteomicsPublicationPdTmtBlock <- function(plan, block) {
    protein_index <- block$entities$protein_index
    annotations <- data.frame(
        Annotated.Sequence = block$entities$peptide_id,
        Master.Protein.Accessions = plan$proteins$protein_id[protein_index],
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(block$values, check.names = FALSE)
    channels <- proteomicsPublicationTmtChannels()[seq_len(nrow(plan$design))]
    names(quantities) <- sprintf(
        "Abundance: F%d: %s, %s",
        seq_along(channels),
        channels,
        plan$design$sample_id
    )
    cbind(annotations, quantities)
}

proteomicsPublicationSerializeBlock <- function(plan, block) {
    switch(
        plan$contract$capability$input_format,
        diann = proteomicsPublicationDiannBlock(plan, block),
        maxquant = proteomicsPublicationMaxQuantBlock(plan, block),
        fragpipe = proteomicsPublicationFragPipeBlock(plan, block),
        pd_tmt = proteomicsPublicationPdTmtBlock(plan, block),
        proteomicsPublicationAbort("serializer dispatch is unsupported")
    )
}

proteomicsPublicationSerialize <- function(plan, path, observer = NULL) {
    if (file.exists(path) || dir.exists(path)) {
        proteomicsPublicationAbort("payload output already exists")
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    entity_count <- nrow(proteomicsPublicationEntityMap(plan))
    ranges <- proteomicsPublicationBlockRanges(
        entity_count,
        as.integer(plan$contract$generator$chunk_rows)
    )
    for (index in seq_along(ranges)) {
        block <- proteomicsPublicationGenerateBlock(plan, ranges[[index]])
        payload <- proteomicsPublicationSerializeBlock(plan, block)
        proteomicsPublicationWriteTable(payload, path, append = index > 1L)
        if (is.function(observer)) observer(block, index)
        rm(block, payload)
        if (index %% 10L == 0L) gc(verbose = FALSE)
    }
    list(
        path = path,
        sha256 = publicationFileDigest(path),
        size_bytes = as.numeric(file.info(path)$size),
        block_count = as.integer(length(ranges)),
        entity_count = as.integer(entity_count)
    )
}
