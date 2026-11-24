# KIT-167 sections 3 (items 3 and 4), 5.6 and 8.3: the reference mapping in the application, behind the option
# `MultiScholaR.reference_mapping` (off by default; nothing here is rendered or registered when it is off).
#   - The import step: an optional "Map to a reference organism" panel under the proteomics source panel, and an
#     optional annotation step at the transcriptomics count import (GFF3/GTF plus a strain protein FASTA to map here or
#     a saved mapping table of the same strain). `runReferenceMappingImport()` and `runReferenceGenesImport()` do the
#     work and store the tables on the head; they return the plain sentence of a refusal instead of raising it.
#   - The enrichment tab: the 2 x 2 counts per contrast (changed and unchanged by mapped and unmapped, with both
#     shares), the unmapped changed features as a download, and the taxon pre-filled from the reference.
# Long work runs in the session, under a progress bar, as the import step's own parsing does (a mapping of a bacterial
# strain takes under a minute at the shipped settings, gate F2); no ExtendedTask is added for it.

#' The import step's reference mapping controls (gate on): proteomics maps a strain FASTA; transcriptomics reads the
#' annotation and takes its mapping from a FASTA mapped here or from a saved mapping table
#' @param ns The import module's namespace function.
#' @param config The omic's GUI configuration.
#' @noRd
referenceMappingImportUi <- function(ns, config) {
    organisms <- readReferenceOrganisms()
    choices <- stats::setNames(organisms$uniprot_upid, organisms$name)
    transcriptomics <- identical(config$shell$omic_type, "transcriptomics")
    shiny::wellPanel(
        shiny::h4(if (transcriptomics) "Gene annotation and reference organism (optional)" else
            "Map to a reference organism (optional)"),
        shiny::helpText(paste("For an individually sequenced strain: its proteins are mapped to their orthologues in a",
            "well-annotated reference of the same species, so enrichment can use the reference's annotation.")),
        shiny::helpText(referenceMappingTexts$applicability),
        if (transcriptomics) shiny::fileInput(ns("refmap_annotation"), "Genome annotation of the strain (GFF3 or GTF):",
            accept = c(".gff", ".gff3", ".gtf", ".gz")),
        shiny::fileInput(ns("refmap_strain_fasta"), "Strain protein FASTA:", accept = c(".fasta", ".fa", ".faa", ".gz")),
        if (transcriptomics) shiny::fileInput(ns("refmap_mapping_table"),
            "Or a saved mapping table of the same strain (TSV):", accept = c(".tsv", ".txt")),
        shiny::selectInput(ns("refmap_reference"), "Reference organism:", choices = choices),
        shiny::numericInput(ns("refmap_coverage"), "Coverage of both sequences (%):", value = 80, min = 1, max = 100),
        shiny::numericInput(ns("refmap_evalue"), "E-value:", value = 1e-5, min = 0),
        shiny::numericInput(ns("refmap_min_identity"), "Minimum identity (%), empty for none:", value = NA,
            min = 0, max = 100),
        shiny::textInput(ns("refmap_strain_taxon"), "Strain taxon (only for UniProt-style files of several organisms):",
            value = ""),
        shiny::actionButton(ns("refmap_run"), if (transcriptomics) "Match genes" else "Map to the reference",
            icon = shiny::icon("project-diagram")),
        shiny::verbatimTextOutput(ns("refmap_summary"), placeholder = TRUE))
}

#' The import step's reference mapping server (gate on)
#' @noRd
referenceMappingImportServer <- function(input, output, session, workflow_data, experiment_paths, config) {
    transcriptomics <- identical(config$shell$omic_type, "transcriptomics")
    summary <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$refmap_run, {
        manager <- workflow_data$state_manager
        keep <- function(file) referenceMappingKeepInput(file, experiment_paths)
        fasta <- keep(input$refmap_strain_fasta)
        settings <- list(coverage = input$refmap_coverage, evalue = input$refmap_evalue,
            min_identity = if (is.null(input$refmap_min_identity) || is.na(input$refmap_min_identity)) NULL else
                input$refmap_min_identity,
            strain_taxon = if (nzchar(trimws(input$refmap_strain_taxon %||% ""))) trimws(input$refmap_strain_taxon))
        done <- shiny::withProgress(message = "Mapping to the reference organism...", value = 0.2, {
            mapped <- if (!is.null(fasta)) runReferenceMappingImport(manager, fasta, input$refmap_reference, settings) else
                list(ok = TRUE, result = NULL)
            if (!isTRUE(mapped$ok) || !transcriptomics) mapped else {
                genes <- runReferenceGenesImport(manager, keep(input$refmap_annotation),
                    mapping_table = keep(input$refmap_mapping_table))
                if (isTRUE(genes$ok)) genes$mapping_result <- mapped$result
                genes
            }
        })
        if (!isTRUE(done$ok)) {
            shiny::showNotification(done$message, type = "error", duration = NULL)
            return(invisible(NULL))
        }
        summary(c(if (!is.null(done$mapping_result %||% (if (!transcriptomics) done$result)))
            referenceMappingSummaryText(done$mapping_result %||% done$result),
            if (transcriptomics) referenceGenesSummaryText(done$result)))
    })
    output$refmap_summary <- shiny::renderText(paste(summary() %||% "", collapse = "\n"))
    invisible(NULL)
}

#' Copy an uploaded input once into the project's `reference_mapping/inputs/` under its own name (section 4)
#' @param file A `fileInput` value, or NULL.
#' @param experiment_paths The project paths (`source_dir`).
#' @return The kept file's path, or NULL.
#' @noRd
referenceMappingKeepInput <- function(file, experiment_paths) {
    if (is.null(file) || !nrow(file)) return(NULL)
    if (is.null(experiment_paths$source_dir)) return(file$datapath[[1L]])
    folder <- file.path(experiment_paths$source_dir, "reference_mapping", "inputs")
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    target <- file.path(folder, basename(file$name[[1L]]))
    file.copy(file$datapath[[1L]], target, overwrite = TRUE)
    target
}

#' Map a strain at the import step and store the mapping on the head; a refusal is returned, not raised
#' @param manager The session owner.
#' @param strain_fasta The kept strain FASTA.
#' @param reference A shipped reference's UPID, or a supplied-reference list.
#' @param settings `coverage`, `evalue`, `min_identity`, `strain_taxon` (missing entries take the defaults).
#' @return `list(ok, result)` or `list(ok = FALSE, class, message)`.
#' @noRd
runReferenceMappingImport <- function(manager, strain_fasta, reference, settings = list()) {
    tryCatch({
        result <- mapStrainToReference(strain_fasta, reference, coverage = settings$coverage %||% 80,
            evalue = settings$evalue %||% 1e-5, strain_taxon = settings$strain_taxon,
            min_identity = settings$min_identity)
        storeReferenceMapping(manager, result, inputs = c(reference_mapping_strain_fasta = strain_fasta))
        list(ok = TRUE, result = result)
    }, reference_mapping_refusal = function(e) list(ok = FALSE, class = class(e)[[1L]], message = conditionMessage(e)))
}

#' Read a saved mapping table (TSV with `strain_accession`, `reference_accession_used`, `map_status`)
#' @noRd
readReferenceMappingTable <- function(path) {
    table <- utils::read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, colClasses = "character",
        na.strings = c("NA", ""))
    needed <- c("strain_accession", "reference_accession_used", "map_status")
    if (!all(needed %in% names(table))) {
        referenceMappingRefuse("reference_mapping_reference_invalid", sprintf("the mapping table lacks the column(s) %s",
            paste(setdiff(needed, names(table)), collapse = ", ")))
    }
    table
}

#' The transcriptomics annotation step: genes matched to proteins and stored; a refusal is returned, not raised
#' @param manager The session owner.
#' @param annotation The kept GFF3 or GTF file.
#' @param mapping_table NULL, or a saved mapping table's path (else the project's stored mapping, if any).
#' @return `list(ok, result)` or `list(ok = FALSE, class, message)`.
#' @noRd
runReferenceGenesImport <- function(manager, annotation, mapping_table = NULL) {
    tryCatch({
        if (is.null(annotation)) referenceMappingRefuse("reference_annotation_missing")
        mapping <- if (!is.null(mapping_table)) readReferenceMappingTable(mapping_table) else
            readReferenceMapping(manager)$mapping
        result <- withCallingHandlers(matchStrainGenesToProteins(annotation, mapping = mapping),
            reference_annotation_not_mapped_share = function(w) {
                if (!is.null(shiny::getDefaultReactiveDomain())) {
                    shiny::showNotification(conditionMessage(w), type = "warning", duration = NULL)
                    invokeRestart("muffleWarning")
                }
            })
        inputs <- c(reference_mapping_annotation = annotation,
            if (!is.null(mapping_table)) c(reference_mapping_table = mapping_table))
        storeReferenceGenes(manager, result, inputs = inputs)
        list(ok = TRUE, result = result)
    }, reference_mapping_refusal = function(e) list(ok = FALSE, class = class(e)[[1L]], message = conditionMessage(e)))
}

#' The import step's summary of a mapping: reference, counts and shares per status, identity summary, caution and
#' applicability (sections 5.3, 5.6)
#' @noRd
referenceMappingSummaryText <- function(result) {
    total <- max(1L, result$counts$strain_proteins)
    statuses <- vapply(referenceMappingStatuses, function(s) {
        n <- result$counts$by_status[[s]] %||% 0L
        sprintf("%s: %d (%.1f%%)", s, n, 100 * n / total)
    }, character(1))
    identity <- result$identity_summary
    c(sprintf("Reference: %s (UniProt release %s)", result$provenance$reference_name %||% "",
            result$provenance$reference_release %||% ""),
        sprintf("Strain proteins: %d", result$counts$strain_proteins), unname(statuses),
        sprintf("Median identity of the one-to-one pairs: %.1f%% (5th percentile %.1f%%; below 40%%: %.1f%%)",
            identity$median, identity$p05, 100 * identity$share_below_40),
        if (isTRUE(identity$caution)) sprintf(referenceMappingTexts$identity_caution, identity$median),
        referenceMappingTexts$applicability)
}

#' The annotation step's summary: genes per status, genes without protein, proteins not in the mapping
#' @noRd
referenceGenesSummaryText <- function(result) {
    counts <- result$counts
    present <- Filter(function(n) n > 0L, counts$by_status)
    c(sprintf("Genes: %d (%s, key %s)", counts$genes, result$provenance$format, result$provenance$gene_key),
        sprintf("%s: %d", names(present), unlist(present)),
        sprintf("Genes without a protein: %d; protein accessions on several genes: %d", counts$genes_without_protein,
            counts$accessions_on_several_genes),
        if (!is.na(counts$protein_not_mapped_share)) sprintf("Annotation proteins not in the mapping: %d of %d",
            counts$protein_not_mapped_proteins, counts$annotation_proteins))
}

#' The reference's UniProt taxon from the stored mapping, for the enrichment tab's pre-fill, or NULL
#' @noRd
referenceMappingTaxon <- function(manager) {
    stored <- tryCatch(readReferenceMapping(manager), error = function(e) NULL)
    if (is.null(stored)) return(NULL)
    taxon <- stored$record$value[stored$record$field == "uniprot_taxon"]
    if (length(taxon) == 1L && nzchar(taxon)) taxon else NULL
}

#' One contrast's 2 x 2 table (section 5.6): changed and unchanged features by mapped and unmapped (ambiguous groups
#' counted as unmapped), with each column's shares
#' @param counts One contrast's counts (`attachReferenceMapping()`).
#' @noRd
referenceMappingContrastTable <- function(counts) {
    changed <- c(counts$changed_mapped, counts$changed_unmapped)
    unchanged <- c(counts$unchanged_mapped, counts$unchanged_unmapped)
    table <- data.frame(features = c("mapped", "unmapped"), changed = as.integer(changed),
        unchanged = as.integer(unchanged), share_of_changed = changed / max(1, sum(changed)),
        share_of_unchanged = unchanged / max(1, sum(unchanged)), stringsAsFactors = FALSE)
    # Audit K167-I09 (KIT-205): rows without a p-value in the comparison are in neither list nor the background.
    untested <- counts$untested %||% 0L
    if (untested > 0L) {
        table <- rbind(table, data.frame(features = "untested (no p-value in this comparison)", changed = NA_integer_,
            unchanged = as.integer(untested), share_of_changed = NA_real_, share_of_unchanged = NA_real_,
            stringsAsFactors = FALSE))
    }
    table
}

#' The unmapped changed features of one contrast, for download
#' @noRd
referenceMappingUnmappedTable <- function(counts) {
    data.frame(protein_group = as.character(counts$unmapped_changed_groups %||% character()), stringsAsFactors = FALSE)
}

#' The enrichment tab's reference mapping outputs (gate on): the 2 x 2 tables, the download and the taxon pre-fill
#' @noRd
registerReferenceMappingEnrichOutputs <- function(input, output, session, workflowData) {
    output$reference_mapping_counts <- shiny::renderTable({
        counts <- workflowData$enrichment_reference_mapping
        shiny::req(length(counts))
        do.call(rbind, lapply(names(counts), function(contrast)
            cbind(contrast = contrast, referenceMappingContrastTable(counts[[contrast]]))))
    }, digits = 3)
    output$download_unmapped_changed <- shiny::downloadHandler(
        filename = function() "unmapped_changed_features.tsv",
        content = function(file) {
            counts <- workflowData$enrichment_reference_mapping %||% list()
            rows <- lapply(names(counts), function(contrast)
                cbind(contrast = rep(contrast, length(counts[[contrast]]$unmapped_changed_groups)),
                    referenceMappingUnmappedTable(counts[[contrast]])))
            table <- if (length(rows)) do.call(rbind, rows) else data.frame(contrast = character(), protein_group = character())
            utils::write.table(table, file, sep = "\t", quote = FALSE, row.names = FALSE)
        })
    shiny::observeEvent(input$id_source, {
        if (!identical(input$id_source, "reference")) return(invisible(NULL))
        taxon <- referenceMappingTaxon(workflowData$state_manager)
        if (!is.null(taxon)) shiny::updateTextInput(session, "organism_taxid", value = taxon)
    })
    invisible(NULL)
}
