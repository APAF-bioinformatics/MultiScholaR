# KIT-167 sections 7, 8.3 and 8.5: the reference mapping in the workflow.
#   - The gate: `referenceMappingEnabled()` (option `MultiScholaR.reference_mapping`, default FALSE).
#   - The mapping's identity (section 7), its record sheet and Methods paragraph.
#   - Store and read on the head: the mapping table and its record sheet as owned tables
#     (`saveWorkflowOwnedValue`, labels `reference_mapping` and `reference_mapping_record`); the input files by role
#     (`recordWorkflowInputs`). No declared operation (8.5).
#   - Replay (the stored table, no engine) and compute-again (the pinned engine on the stored inputs, compared by the
#     canonical text's SHA-256).
#   - Proteomics enrichment with the reference identifiers (8.3): the protein-group join (`attachReferenceMapping`),
#     the reference's annotation table (`referenceAnnotationTable`), the 2 x 2 unmapped counts.

#' Whether the reference mapping tool is switched on (option `MultiScholaR.reference_mapping`, default FALSE)
#' @noRd
referenceMappingEnabled <- function() isTRUE(getOption("MultiScholaR.reference_mapping", FALSE))

#' The elements of a mapping's identity (section 7), from its provenance
#' @noRd
referenceMappingIdentityFields <- c("strain_selection_sha256", "reference_fasta_sha256", "reference_xref_sha256",
    "annotation_sha256", "engine_version", "engine_sha256", "engine_arguments", "coverage", "evalue", "min_identity",
    "header_classes", "strain_taxon", "ranking_key", "resolver_version")

#' The identity of a mapping: the section 7 elements of its provenance, and their SHA-256
#' @param provenance The `provenance` list of `mapStrainToReference()` (with `annotation_sha256` when a gene table
#'   was made from it).
#' @return `list(identity, sha256)`.
#' @noRd
referenceMappingIdentity <- function(provenance) {
    identity <- lapply(stats::setNames(referenceMappingIdentityFields, referenceMappingIdentityFields),
        function(field) provenance[[field]] %||% "none")
    list(identity = identity, sha256 = digest::digest(identity, algo = "sha256", serializeVersion = 2L))
}

#' The record sheet of a mapping (section 7, A14): one row per field, every value as text
#' @param result `mapStrainToReference()`'s list.
#' @return Data frame (field, value).
#' @noRd
referenceMappingRecordSheet <- function(result) {
    p <- result$provenance
    counts <- result$counts
    identity <- result$identity_summary
    text <- function(x) if (is.null(x) || !length(x)) "none" else paste(as.character(x), collapse = ";")
    statuses <- vapply(referenceMappingStatuses, function(s) text(counts$by_status[[s]]), character(1))
    shares <- vapply(referenceMappingStatuses, function(s)
        formatC(counts$by_status[[s]] / max(1L, counts$strain_proteins), digits = 17, format = "g"), character(1))
    arguments <- vapply(p$engine_arguments, paste, character(1), collapse = " ")
    rows <- c(
        reference = p$reference_name, uniprot_upid = p$reference_upid, uniprot_release = p$reference_release,
        reference_supplied = text(p$reference_supplied), uniprot_taxon = text(p$reference_uniprot_taxon),
        string_species = text(p$reference_string_species), string_status = p$reference_string_status,
        kegg_code = p$reference_kegg_code, strain_fasta_sha256 = p$strain_fasta_sha256,
        strain_selection_sha256 = p$strain_selection_sha256, reference_fasta_sha256 = p$reference_fasta_sha256,
        reference_xref_sha256 = p$reference_xref_sha256, engine = p$engine, engine_version = p$engine_version,
        engine_sha256 = p$engine_sha256, engine_platform = p$engine_platform,
        stats::setNames(arguments, paste0("engine_arguments_", names(arguments))),
        coverage = text(p$coverage), evalue = text(p$evalue), min_identity = text(p$min_identity),
        header_classes = text(p$header_classes), strain_taxon = text(p$strain_taxon), ranking_key = p$ranking_key,
        resolver_version = p$resolver_version, strain_proteins = text(counts$strain_proteins),
        identical_sequence_groups = text(counts$identical_sequence_groups),
        stats::setNames(statuses, paste0("count_", referenceMappingStatuses)),
        stats::setNames(shares, paste0("share_", referenceMappingStatuses)),
        one_direction_pairs = text(counts$one_direction_pairs),
        used_with_string = text(counts$used_with_string), used_with_kegg = text(counts$used_with_kegg),
        identity_pairs = text(identity$pairs), identity_median = text(identity$median),
        identity_p05 = text(identity$p05), identity_share_below_40 = text(identity$share_below_40),
        identity_caution = text(identity$caution), applicability = referenceMappingTexts$applicability,
        mapping_sha256 = p$mapping_sha256, identity_sha256 = referenceMappingIdentity(p)$sha256)
    data.frame(field = names(rows), value = unname(rows), stringsAsFactors = FALSE)
}

#' The Methods paragraph of a mapping (section 7), built from its record sheet
#' @param sheet `referenceMappingRecordSheet()`'s table.
#' @return One string.
#' @noRd
referenceMappingMethodsText <- function(sheet) {
    v <- stats::setNames(sheet$value, sheet$field)
    one_to_one <- as.numeric(v[["count_one_to_one"]])
    total <- as.numeric(v[["strain_proteins"]])
    reference <- if (identical(v[["reference_supplied"]], "TRUE")) {
        sprintf("%s (UniProt proteome supplied by the user)", v[["reference"]])
    } else {
        sprintf("%s (UniProt reference proteome %s, UniProtKB release %s, shipped with MultiScholaR)",
            v[["reference"]], v[["uniprot_upid"]], v[["uniprot_release"]])
    }
    floor <- if (identical(v[["min_identity"]], "none")) "no identity floor" else
        sprintf("a minimum identity of %s %%", v[["min_identity"]])
    paste0("Strain proteins were mapped to their orthologues in ", reference, " by reciprocal best hits of DIAMOND ",
        v[["engine_version"]], " (Buchfink et al. 2021) searches in both directions (ultra-sensitive mode, repeat ",
        "masking off, full Smith-Waterman extension, composition-based statistics mode 1, e-value <= ", v[["evalue"]],
        ", at least ", v[["coverage"]], " % of both sequences covered by one alignment, ", floor, "). A strain protein ",
        "took the reference accession only when each protein was the other's best candidate (one-to-one); ties, ",
        "one-way hits, many-to-one hits and proteins without a candidate were left unmapped. ",
        sprintf("%d of %d strain proteins (%.1f %%) were mapped one-to-one; ", as.integer(one_to_one), as.integer(total),
            100 * one_to_one / max(1, total)),
        sprintf("the median identity of the pairs was %s %%. ", formatC(as.numeric(v[["identity_median"]]),
            digits = 4, format = "g")),
        v[["applicability"]], " The reference files are kept in the project and should be deposited with the data.")
}

#' Store a mapping on the head (section 7): the mapping table and its record sheet as owned tables, the input files
#' by role. Raw DIAMOND outputs are not stored.
#' @param manager The session owner.
#' @param result `mapStrainToReference()`'s list.
#' @param inputs Named character vector, role to file (`reference_mapping_strain_fasta`, ...).
#' @return `list(mapping, record)`: the saved references.
#' @noRd
storeReferenceMapping <- function(manager, result, inputs = character()) {
    withWorkflowStateScope(manager, function() {
        mapping <- saveWorkflowOwnedValue(result$mapping, "reference_mapping", manager)
        record <- saveWorkflowOwnedValue(referenceMappingRecordSheet(result), "reference_mapping_record", manager)
        if (length(inputs)) recordWorkflowInputs(manager, inputs)
        list(mapping = mapping, record = record)
    })
}

#' Read the mapping stored on the head, or NULL when none is
#' @param manager The session owner.
#' @return `list(mapping, record)`, the tables, or NULL.
#' @noRd
readReferenceMapping <- function(manager) {
    if (!usesWorkflowReferenceStates(manager)) return(NULL)
    references <- manager$.context$snapshot$gui_references
    if (is.null(references$reference_mapping) || is.null(references$reference_mapping_record)) return(NULL)
    list(mapping = resolveWorkflowOwnedValue(references$reference_mapping, manager),
        record = resolveWorkflowOwnedValue(references$reference_mapping_record, manager))
}

#' A stored mapping is reused only when every element of its identity matches (section 7)
#' @param stored `readReferenceMapping()`'s list, or NULL.
#' @param identity_sha256 The identity SHA-256 of the mapping about to be made.
#' @return TRUE when the stored mapping has that identity.
#' @noRd
referenceMappingReusable <- function(stored, identity_sha256) {
    if (is.null(stored)) return(FALSE)
    identical(stored$record$value[stored$record$field == "identity_sha256"], identity_sha256)
}

#' Replay: the stored mapping table, checked against its recorded canonical SHA-256; no engine is run
#' @param manager The session owner.
#' @return The mapping table.
#' @noRd
replayReferenceMapping <- function(manager) {
    stored <- readReferenceMapping(manager)
    if (is.null(stored)) cli::cli_abort("the project holds no reference mapping", class = "reference_mapping_replay")
    recorded <- stored$record$value[stored$record$field == "mapping_sha256"]
    if (!identical(canonicalReferenceMappingSha256(stored$mapping), recorded)) {
        cli::cli_abort("the stored reference mapping does not match its recorded SHA-256",
            class = "reference_mapping_replay")
    }
    stored$mapping
}

#' Compute again: the pinned engine on the stored inputs with the recorded settings; compared with the stored table by
#' the canonical text's SHA-256. Without the engine it is refused in words (`reference_mapping_engine_missing`).
#' @param manager The session owner.
#' @param strain_fasta,reference The stored inputs (the reference as `mapStrainToReference()` takes it).
#' @return `list(identical, stored_sha256, computed_sha256, result)`.
#' @noRd
computeReferenceMappingAgain <- function(manager, strain_fasta, reference) {
    stored <- readReferenceMapping(manager)
    if (is.null(stored)) cli::cli_abort("the project holds no reference mapping", class = "reference_mapping_replay")
    v <- stats::setNames(stored$record$value, stored$record$field)
    number <- function(x) if (identical(x, "none")) NULL else as.numeric(x)
    result <- mapStrainToReference(strain_fasta, reference, coverage = as.numeric(v[["coverage"]]),
        evalue = as.numeric(v[["evalue"]]), header_classes = strsplit(v[["header_classes"]], ";", fixed = TRUE)[[1]],
        strain_taxon = number(v[["strain_taxon"]]), min_identity = number(v[["min_identity"]]))
    computed <- result$provenance$mapping_sha256
    list(identical = identical(computed, v[["mapping_sha256"]]), stored_sha256 = v[["mapping_sha256"]],
        computed_sha256 = computed, result = result)
}

#' The reference accession of each protein group by the protein-group rule (section 8.3 items 1 and 8)
#'
#' A group's members are split on `;` and a `ref|...|` wrapper removed; the members found among the mapping's strain
#' accessions are kept (strain entries by entry class, never by prefix: other members never match). If every kept
#' member with a `reference_accession_used` agrees, that accession; if they disagree, NA (`ambiguous_group`); if none
#' has one, NA (`unmapped`).
#' @param groups Character vector of protein-group identifiers.
#' @param mapping The mapping table (strain_accession, reference_accession_used).
#' @return Data frame: group, reference_accession_used, group_class (`mapped`, `ambiguous_group`, `unmapped`).
#' @noRd
referenceAccessionForGroups <- function(groups, mapping) {
    members <- strsplit(as.character(groups), ";", fixed = TRUE)
    lookup <- stats::setNames(mapping$reference_accession_used, mapping$strain_accession)
    resolved <- vapply(members, function(m) {
        m <- sub("^ref\\|([^|]+)\\|?.*$", "\\1", trimws(m))
        used <- unique(lookup[m[m %in% names(lookup)]])
        used <- used[!is.na(used)]
        if (length(used) == 1L) used else if (length(used) > 1L) "\r ambiguous" else NA_character_
    }, character(1))
    class <- ifelse(is.na(resolved), "unmapped", ifelse(resolved == "\r ambiguous", "ambiguous_group", "mapped"))
    resolved[class != "mapped"] <- NA_character_
    data.frame(group = as.character(groups), reference_accession_used = resolved, group_class = class,
        stringsAsFactors = FALSE)
}

#' The protein-group join of the enrichment (section 8.3 items 1 and 4): each contrast's table gains
#' `reference_accession_used`; rows without one are removed, so the query lists and the background are "quantified
#' and mapped". The tables on disk are not changed.
#' @param da_results The `da_results_for_enrichment` object the run read.
#' @param mapping The mapping table.
#' @param id_column The tables' protein-group column.
#' @param q_cutoff,up_cutoff,down_cutoff The run's cutoffs, for the 2 x 2 counts (changed: q below the cutoff and
#'   log2 fold change at or beyond a cutoff).
#' @return `list(da_results, counts)`: counts per contrast (groups by class; the 2 x 2 changed/unchanged by
#'   mapped/unmapped) and the unmapped changed groups.
#' @noRd
attachReferenceMapping <- function(da_results, mapping, id_column, q_cutoff = 0.05, up_cutoff = 0, down_cutoff = 0) {
    counts <- list()
    for (contrast in names(da_results@da_data)) {
        table <- da_results@da_data[[contrast]]
        if (!is.data.frame(table) || !id_column %in% names(table)) next
        joined <- referenceAccessionForGroups(table[[id_column]], mapping)
        changed <- referenceMappingChanged(table, q_cutoff, up_cutoff, down_cutoff)
        mapped <- joined$group_class == "mapped"
        counts[[contrast]] <- list(
            groups = nrow(table), mapped = sum(mapped), ambiguous_group = sum(joined$group_class == "ambiguous_group"),
            unmapped = sum(joined$group_class == "unmapped"),
            changed_mapped = sum(changed & mapped), changed_unmapped = sum(changed & !mapped),
            unchanged_mapped = sum(!changed & mapped), unchanged_unmapped = sum(!changed & !mapped),
            unmapped_changed_groups = sort(as.character(table[[id_column]][changed & !mapped]), method = "radix"))
        table$reference_accession_used <- joined$reference_accession_used
        da_results@da_data[[contrast]] <- table[mapped, , drop = FALSE]
    }
    list(da_results = da_results, counts = counts)
}

#' Which rows of a DA table are changed (q below the cutoff and log2 fold change at or beyond the up or down cutoff)
#' @noRd
referenceMappingChanged <- function(table, q_cutoff, up_cutoff, down_cutoff) {
    q <- table[["fdr_qvalue"]] %||% rep(NA_real_, nrow(table))
    fc <- table[["log2FC"]] %||% rep(NA_real_, nrow(table))
    !is.na(q) & !is.na(fc) & q < q_cutoff & (fc >= up_cutoff | fc <= -down_cutoff)
}

#' The shipped GO capture of a curated reference, converted as `getUniprotAnnotations` converts a download (8.3 item 3)
#'
#' Read as `directUniprotDownload` reads a reply; `go_ids_not_in_go_db` counted on the raw ids (KIT-173); GO ids to
#' terms (`uniprotGoIdToTerm`) and the standard columns (`standardizeUniprotColumns`); the attribute set last.
#' @param path The shipped `<UPID>_go.tsv`.
#' @return The converted annotation table, with attribute `go_ids_not_in_go_db`.
#' @noRd
readReferenceGoTable <- function(path) {
    table <- utils::read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE)
    names(table) <- gsub(" ", ".", names(table))
    table$From <- table$Entry
    lost <- goIdsNotInGoDb(unlist(strsplit(as.character(table$Gene.Ontology.IDs), "; ", fixed = TRUE)))
    processed <- uniprotGoIdToTerm(table, uniprot_id_column = Entry, go_id_column = Gene.Ontology.IDs, sep = "; ")
    processed <- standardizeUniprotColumns(processed)
    attr(processed, "go_ids_not_in_go_db") <- lost
    processed
}

#' The annotation table of the reference identifiers (8.3 item 3): the shipped GO capture for a curated reference; for
#' a user-supplied or current-release reference one `getUniprotAnnotations()` fetch on the mapping's
#' `reference_accession_used`, cached under the mapping's SHA-256, its record rows labelled with the reference.
#' @param mapping The stored mapping table.
#' @param record The stored record sheet.
#' @param project_dir The project folder (the cache's parent), for a fetched table.
#' @param getUniprotAnnotationsFn The fetch (tests pass a stub).
#' @return `list(table, go_table_sha256, source)`.
#' @noRd
referenceAnnotationTable <- function(mapping, record, project_dir = NULL,
                                     getUniprotAnnotationsFn = getUniprotAnnotations) {
    v <- stats::setNames(record$value, record$field)
    if (!identical(v[["reference_supplied"]], "TRUE")) {
        path <- referenceOrganismFiles(v[["uniprot_upid"]])$go
        return(list(table = readReferenceGoTable(path), go_table_sha256 = digest::digest(file = path, algo = "sha256"),
            source = "shipped"))
    }
    used <- sort(unique(mapping$reference_accession_used[!is.na(mapping$reference_accession_used)]), method = "radix")
    cache_dir <- file.path(project_dir %||% tempdir(), "reference_mapping", "annotation", v[["mapping_sha256"]])
    label <- sprintf("reference annotation: %s (%s)", v[["reference"]],
        if (identical(v[["uniprot_upid"]], "supplied")) "supplied" else v[["uniprot_upid"]])
    table <- getUniprotAnnotationsFn(input_tbl = data.frame(Protein.Ids = used, stringsAsFactors = FALSE),
        cache_dir = cache_dir, taxon_id = as.numeric(v[["uniprot_taxon"]]), label = label)
    sha <- if (is.data.frame(table) && nrow(table)) canonicalReferenceMappingSha256(as.data.frame(table), "Entry") else
        NA_character_
    list(table = table, go_table_sha256 = sha, source = "fetched")
}
