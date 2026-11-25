# KIT-167 section 5.5: a strain's genes (locus tags) matched to its proteins from the genome annotation (GFF3 or GTF),
# without sequence search, and joined with a protein mapping by the gene-level rule.
#   - `matchStrainGenesToProteins()` (exported).
#   - `readAnnotationRows()`, `readGff3Cds()`, `readGtfCds()`: a base-R parser; malformed lines refused by number.
#   - `resolveGeneStatuses()`: `no_protein`, `protein_not_mapped`, `ambiguous_gene`, else the proteins' status.
#   - `joinGenesToCountTable()`: the exact join to a count table's gene identifier, with both sides' absences counted.

#' The gene statuses of 5.5 besides the protein statuses of 5.3
#' @noRd
referenceAnnotationStatuses <- c("no_protein", "protein_not_mapped", "ambiguous_gene")

#' The annotation format from the first feature line's attribute syntax (`key=value` GFF3, `key "value";` GTF)
#' @param path The annotation file.
#' @return `"gff3"` or `"gtf"`.
#' @noRd
detectAnnotationFormat <- function(path) {
    connection <- file(path, open = "r")
    on.exit(close(connection), add = TRUE)
    repeat {
        lines <- readLines(connection, n = 200L, warn = FALSE)
        if (!length(lines)) break
        lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
        if (length(lines)) {
            attributes <- vapply(strsplit(lines[[1L]], "\t", fixed = TRUE), function(f) f[9L] %||% "", character(1))
            if (grepl("^[[:space:]]*[A-Za-z_][A-Za-z0-9_]*[[:space:]]+\"", attributes)) return("gtf")
            return("gff3")
        }
    }
    if (grepl("\\.gtf(\\.gz)?$", path, ignore.case = TRUE)) "gtf" else "gff3"
}

#' Read the feature lines of a GFF3 or GTF file; a line without 9 tab-separated fields is refused with its number
#' @param path The annotation file (plain or gzip).
#' @return Data frame: seqid, type, start, end, attributes, line.
#' @noRd
readAnnotationRows <- function(path) {
    connection <- file(path, open = "r")
    on.exit(close(connection), add = TRUE)
    lines <- sub("\r$", "", readLines(connection, warn = FALSE))
    fasta <- which(lines == "##FASTA")
    if (length(fasta)) lines <- lines[seq_len(fasta[[1L]] - 1L)]
    number <- seq_along(lines)
    keep <- nzchar(lines) & !startsWith(lines, "#")
    fields <- strsplit(lines[keep], "\t", fixed = TRUE)
    widths <- lengths(fields)
    if (any(widths != 9L)) referenceMappingRefuse("reference_annotation_malformed", number[keep][which(widths != 9L)[1L]],
        fields = list(line = number[keep][which(widths != 9L)[1L]]))
    column <- function(i) vapply(fields, `[[`, character(1), i)
    data.frame(seqid = column(1L), type = column(3L), start = as.numeric(column(4L)), end = as.numeric(column(5L)),
        attributes = column(9L), line = number[keep], stringsAsFactors = FALSE)
}

#' One attribute of each row (GFF3 values percent-decoded)
#' @noRd
annotationAttribute <- function(attributes, key, format) {
    pattern <- if (identical(format, "gtf")) sprintf("(^|;)[[:space:]]*%s[[:space:]]+\"([^\"]*)\"", key) else
        sprintf("(^|;)%s=([^;]*)", key)
    found <- regmatches(attributes, regexec(pattern, attributes))
    value <- vapply(found, function(m) if (length(m) == 3L) m[[3L]] else NA_character_, character(1))
    if (!identical(format, "gtf")) {
        encoded <- !is.na(value) & grepl("%", value, fixed = TRUE)
        value[encoded] <- vapply(value[encoded], utils::URLdecode, character(1), USE.NAMES = FALSE)
    }
    value[!is.na(value) & !nzchar(value)] <- NA_character_
    value
}

#' Total bases in the union of intervals
#' @noRd
annotationIntervalUnion <- function(start, end) {
    order_at <- order(start, end)
    start <- start[order_at]
    end <- end[order_at]
    total <- 0
    current_start <- start[[1L]]
    current_end <- end[[1L]]
    for (i in seq_along(start)[-1L]) {
        if (start[[i]] > current_end + 1) {
            total <- total + current_end - current_start + 1
            current_start <- start[[i]]
            current_end <- end[[i]]
        } else {
            current_end <- max(current_end, end[[i]])
        }
    }
    total + current_end - current_start + 1
}

#' The CDS, gene and exon rows of a GFF3 file, keyed by locus tag
#' @noRd
readGff3Cds <- function(path) {
    rows <- readAnnotationRows(path)
    rows$key <- annotationAttribute(rows$attributes, "locus_tag", "gff3")
    rows$protein_id <- annotationAttribute(rows$attributes, "protein_id", "gff3")
    rows$gene <- annotationAttribute(rows$attributes, "gene", "gff3")
    rows$biotype <- annotationAttribute(rows$attributes, "gene_biotype", "gff3")
    rows$pseudo <- annotationAttribute(rows$attributes, "pseudo", "gff3") %in% "true"
    rows$locus_tag <- rows$key
    list(rows = rows, format = "gff3", gene_key = "locus_tag")
}

#' The CDS, gene and exon rows of a GTF file, keyed by `locus_tag` when every CDS row has one, else `gene_id`
#' @noRd
readGtfCds <- function(path) {
    rows <- readAnnotationRows(path)
    locus <- annotationAttribute(rows$attributes, "locus_tag", "gtf")
    gene_id <- annotationAttribute(rows$attributes, "gene_id", "gtf")
    cds <- rows$type == "CDS"
    gene_key <- if (any(cds) && all(!is.na(locus[cds]))) "locus_tag" else "gene_id"
    rows$key <- if (identical(gene_key, "locus_tag")) locus else gene_id
    rows$locus_tag <- locus
    rows$protein_id <- annotationAttribute(rows$attributes, "protein_id", "gtf")
    rows$gene <- annotationAttribute(rows$attributes, "gene", "gtf")
    rows$biotype <- annotationAttribute(rows$attributes, "gene_biotype", "gtf")
    rows$pseudo <- annotationAttribute(rows$attributes, "pseudo", "gtf") %in% "true"
    list(rows = rows, format = "gtf", gene_key = gene_key)
}

#' The gene status of each gene from its proteins' statuses (section 5.5)
#'
#' `no_protein`: no CDS with a protein (pseudogenes, RNA genes); `protein_not_mapped`: a protein the annotation names
#' is absent from the mapped FASTA; `ambiguous_gene`: several mapped proteins whose `reference_accession_used` differ
#' (one with an accession and one without included); otherwise the proteins' status (when they agree on the accession
#' but differ in status, the status of the first protein in C-locale order).
#' @param proteins List of each gene's protein accessions.
#' @param mapping NULL or the mapping table (strain_accession, reference_accession_used, map_status).
#' @return Data frame: reference_accession_used, map_status.
#' @noRd
resolveGeneStatuses <- function(proteins, mapping = NULL) {
    out <- data.frame(reference_accession_used = rep(NA_character_, length(proteins)),
        map_status = rep(NA_character_, length(proteins)), stringsAsFactors = FALSE)
    for (i in seq_along(proteins)) {
        accessions <- proteins[[i]]
        if (!length(accessions)) {
            out$map_status[i] <- "no_protein"
            next
        }
        if (is.null(mapping)) next
        at <- match(accessions, mapping$strain_accession)
        if (anyNA(at)) {
            out$map_status[i] <- "protein_not_mapped"
            next
        }
        used <- mapping$reference_accession_used[at]
        if (length(unique(used)) > 1L) {
            out$map_status[i] <- "ambiguous_gene"
            next
        }
        out$reference_accession_used[i] <- used[[1L]]
        out$map_status[i] <- mapping$map_status[at][[1L]]
    }
    out
}

#' Match a strain's genes to its proteins from the genome annotation
#'
#' Reads a GFF3 or GTF file with base R (no sequence search): one row per gene (locus tag) with its protein
#' accession(s), gene name, biotype and gene length. With a protein mapping (`mapStrainToReference()`), each gene takes
#' `reference_accession_used` and a status by the gene-level rule: `no_protein`, `protein_not_mapped`,
#' `ambiguous_gene`, or the status its proteins share.
#'
#' Gene length is the bases in the union of the gene's exon rows, else of its CDS rows, else its gene row's span
#' (`length_source`); a GTF's CDS excludes the stop codon and a GFF3's includes it.
#'
#' @param annotation Path to the GFF3 or GTF file (plain or gzip).
#' @param format `"auto"` (from the attribute syntax), `"gff3"` or `"gtf"`.
#' @param mapping NULL, a mapping table, or the list `mapStrainToReference()` returns.
#' @param max_not_mapped_share Share of the annotation's proteins absent from the mapping above which a warning is
#'   given (default 0.05).
#' @return `list(genes, counts, provenance)`: `genes` (gene_id, locus_tag, protein_accession, gene_name, biotype,
#'   gene_length, length_source, reference_accession_used, map_status), counts (genes per status summing to the
#'   genes, CDS rows, genes without protein, accessions on several genes, proteins not mapped and their share), and
#'   the file's identity, format and gene key.
#' @export
matchStrainGenesToProteins <- function(annotation, format = c("auto", "gff3", "gtf"), mapping = NULL,
                                       max_not_mapped_share = 0.05) {
    format <- match.arg(format)
    checkmate::assert_file_exists(annotation)
    checkmate::assert_number(max_not_mapped_share, lower = 0, upper = 1)
    if (is.list(mapping) && !is.data.frame(mapping)) mapping <- mapping$mapping
    if (identical(format, "auto")) format <- detectAnnotationFormat(annotation)
    parsed <- if (identical(format, "gtf")) readGtfCds(annotation) else readGff3Cds(annotation)
    rows <- parsed$rows
    rows <- rows[!is.na(rows$key), , drop = FALSE]
    cds <- rows[rows$type == "CDS", , drop = FALSE]
    if (!any(!is.na(cds$protein_id))) referenceMappingRefuse("reference_annotation_no_protein_id")
    gene_rows <- rows[rows$type %in% c("gene", "pseudogene"), , drop = FALSE]
    exons <- rows[rows$type == "exon", , drop = FALSE]
    keys <- sort(unique(c(gene_rows$key, cds$key, exons$key)), method = "radix")
    pseudo_genes <- unique(c(gene_rows$key[gene_rows$pseudo | gene_rows$type == "pseudogene"],
        gene_rows$key[gene_rows$biotype %in% "pseudogene"]))
    proteins <- lapply(keys, function(k) {
        if (k %in% pseudo_genes) return(character())
        ids <- cds$protein_id[cds$key == k & !cds$pseudo]
        sort(unique(ids[!is.na(ids)]), method = "radix")
    })
    firstValue <- function(table, k, field) {
        v <- table[[field]][table$key == k]
        v <- v[!is.na(v)]
        if (length(v)) v[[1L]] else NA_character_
    }
    lengths_and_source <- lapply(keys, function(k) {
        for (part in list(list(exons, "exon"), list(cds, "cds"))) {
            at <- part[[1L]]$key == k
            if (any(at)) return(list(annotationIntervalUnion(part[[1L]]$start[at], part[[1L]]$end[at]),
                paste0(part[[2L]], "_", parsed$format)))
        }
        at <- which(gene_rows$key == k)[1L]
        list(gene_rows$end[at] - gene_rows$start[at] + 1, paste0("gene_", parsed$format))
    })
    statuses <- resolveGeneStatuses(proteins, mapping)
    genes <- data.frame(gene_id = keys,
        locus_tag = vapply(keys, function(k) firstValue(rows, k, "locus_tag"), character(1), USE.NAMES = FALSE),
        protein_accession = vapply(proteins, function(p) if (length(p)) paste(p, collapse = ";") else NA_character_,
            character(1)),
        gene_name = vapply(keys, function(k) {
            name <- firstValue(gene_rows, k, "gene")
            if (is.na(name)) firstValue(cds, k, "gene") else name
        }, character(1), USE.NAMES = FALSE),
        biotype = vapply(keys, function(k) firstValue(gene_rows, k, "biotype"), character(1), USE.NAMES = FALSE),
        gene_length = vapply(lengths_and_source, function(x) as.numeric(x[[1L]]), numeric(1)),
        length_source = vapply(lengths_and_source, `[[`, character(1), 2L),
        reference_accession_used = statuses$reference_accession_used, map_status = statuses$map_status,
        stringsAsFactors = FALSE)
    annotation_proteins <- sort(unique(unlist(proteins)), method = "radix")
    not_mapped <- if (is.null(mapping)) NA_integer_ else sum(!annotation_proteins %in% mapping$strain_accession)
    share <- if (is.null(mapping) || !length(annotation_proteins)) NA_real_ else not_mapped / length(annotation_proteins)
    status_names <- c(referenceAnnotationStatuses, referenceMappingStatuses)
    counted <- ifelse(is.na(genes$map_status), "no_mapping", genes$map_status)
    on_genes <- table(unlist(lapply(proteins, unique)))
    counts <- list(genes = nrow(genes),
        by_status = lapply(stats::setNames(c(status_names, "no_mapping"), c(status_names, "no_mapping")),
            function(s) sum(counted == s)),
        cds_rows = nrow(cds), genes_without_protein = sum(lengths(proteins) == 0L),
        annotation_proteins = length(annotation_proteins), accessions_on_several_genes = sum(on_genes > 1L),
        protein_not_mapped = sum(genes$map_status %in% "protein_not_mapped"), protein_not_mapped_proteins = not_mapped,
        protein_not_mapped_share = share, max_not_mapped_share = max_not_mapped_share)
    if (!is.na(share) && share > max_not_mapped_share) {
        rlang::warn(sprintf(referenceMappingTexts$reference_annotation_not_mapped_share, not_mapped,
            length(annotation_proteins), 100 * share), class = "reference_annotation_not_mapped_share")
    }
    provenance <- list(annotation_sha256 = digest::digest(file = annotation, algo = "sha256"),
        annotation_bytes = file.size(annotation), format = parsed$format, gene_key = parsed$gene_key)
    list(genes = genes, counts = counts, provenance = provenance)
}

#' Join a gene table to a count table's gene identifiers, exactly (no prefix stripping), counting absences
#' @param genes `matchStrainGenesToProteins()`'s gene table.
#' @param gene_ids The count table's gene identifiers.
#' @return `list(in_both, count_only, annotation_only, match)`: `match` gives each count-table gene's row in `genes`.
#' @noRd
joinGenesToCountTable <- function(genes, gene_ids) {
    ids <- unique(as.character(gene_ids))
    list(in_both = sum(ids %in% genes$gene_id), count_only = sum(!ids %in% genes$gene_id),
        annotation_only = sum(!genes$gene_id %in% ids), match = match(as.character(gene_ids), genes$gene_id))
}
