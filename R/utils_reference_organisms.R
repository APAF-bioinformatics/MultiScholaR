# KIT-167 section 5.7: the four curated references shipped with the package (UniProt release 2026_04, CC BY 4.0, see
# inst/extdata/reference_mapping/ATTRIBUTION.txt), their service identifiers, and the parsing of their multi-valued
# cross-reference fields. `downloadReferenceProteome()` is the explicit "current UniProt release" choice.

#' The shipped reference list (`inst/extdata/reference_organisms.tsv`)
#' @return Data frame, one row per curated reference, in the order of section 5.7 (the reference index).
#' @noRd
readReferenceOrganisms <- function() {
    path <- system.file("extdata", "reference_organisms.tsv", package = "MultiScholaR")
    utils::read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings = "NA",
        colClasses = c(name = "character", uniprot_upid = "character", uniprot_taxon = "integer",
            uniprot_release = "character", entries = "integer", fasta_sha256 = "character", xref_sha256 = "character",
            go_sha256 = "character", string_species = "character", string_status = "character",
            string_recognised_share = "numeric", kegg_code = "character", kegg_taxon = "integer",
            locus_pattern = "character"))
}

#' The shipped files of one curated reference
#' @param upid The UniProt proteome id.
#' @return `list(fasta, xref, go)`: installed paths.
#' @noRd
referenceOrganismFiles <- function(upid) {
    folder <- system.file("extdata", "reference_mapping", "references", package = "MultiScholaR")
    files <- list(fasta = file.path(folder, paste0(upid, ".fasta")), xref = file.path(folder, paste0(upid, "_xref.tsv")),
        go = file.path(folder, paste0(upid, "_go.tsv")))
    if (!nzchar(folder) || !all(file.exists(unlist(files)))) {
        referenceMappingRefuse("reference_mapping_reference_invalid",
            sprintf("the shipped files of %s are missing from the installation", upid))
    }
    files
}

#' Read a cross-reference table (UniProt TSV: Entry, Gene Names (primary), Gene Names (ordered locus), STRING, KEGG)
#' @param path The table.
#' @return Data frame, every column as text, empty fields as "".
#' @noRd
readReferenceCrossReferenceTable <- function(path) {
    table <- utils::read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, colClasses = "character",
        na.strings = character(), check.names = FALSE)
    needed <- c("Entry", "Gene Names (primary)", "Gene Names (ordered locus)", "STRING", "KEGG")
    if (!all(needed %in% names(table))) {
        referenceMappingRefuse("reference_mapping_reference_invalid",
            sprintf("the cross-reference table lacks the column(s) %s", paste(setdiff(needed, names(table)),
                collapse = ", ")))
    }
    table
}

#' Split a `;`-separated UniProt cross-reference field and keep the entries with a prefix, in C-locale order
#' @noRd
referenceXrefEntries <- function(field, prefix) {
    lapply(strsplit(field, ";", fixed = TRUE), function(values) {
        values <- trimws(values)
        values <- values[nzchar(values) & startsWith(values, prefix)]
        sort(unique(values), method = "radix")
    })
}

#' Parse the multi-valued cross-reference fields of a reference (section 5.7)
#'
#' `reference_kegg_id`: the entries of the KEGG field whose prefix is the reference's KEGG code;
#' `reference_string_id`: the entry of the STRING field whose prefix is the STRING species; `reference_locus`: the
#' locus of `reference_kegg_id`, else the ordered-locus tokens matching the reference's locus pattern. Several values are
#' `;`-joined in C-locale order and counted.
#' @param xref The cross-reference table (`readReferenceCrossReferenceTable`).
#' @param kegg_code The reference's KEGG organism code (`eco`, `sao`, `kpn`, `spy`).
#' @param string_species The STRING species id, or NA when STRING is unavailable for the reference.
#' @param locus_pattern Regular expression of the reference's ordered-locus names, anchored.
#' @return Data frame: accession, reference_gene, reference_locus, reference_string_id, reference_kegg_id,
#'   reference_locus_values, reference_kegg_values.
#' @noRd
parseReferenceCrossReferences <- function(xref, kegg_code, string_species, locus_pattern) {
    joined <- function(values) vapply(values, function(v) if (length(v)) paste(v, collapse = ";") else NA_character_,
        character(1))
    kegg <- referenceXrefEntries(xref[["KEGG"]], paste0(kegg_code, ":"))
    string <- if (is.na(string_species)) rep(list(character()), nrow(xref)) else
        lapply(referenceXrefEntries(xref[["STRING"]], paste0(string_species, ".")), utils::head, 1L)
    oln <- strsplit(xref[["Gene Names (ordered locus)"]], "[ ;,/]+")
    locus <- lapply(seq_len(nrow(xref)), function(i) {
        if (length(kegg[[i]])) return(sort(unique(sub("^[^:]+:", "", kegg[[i]])), method = "radix"))
        tokens <- oln[[i]][grepl(locus_pattern, oln[[i]])]
        sort(unique(tokens), method = "radix")
    })
    gene <- xref[["Gene Names (primary)"]]
    data.frame(accession = xref[["Entry"]], reference_gene = ifelse(nzchar(gene), gene, NA_character_),
        reference_locus = unname(joined(locus)), reference_string_id = unname(joined(string)),
        reference_kegg_id = unname(joined(kegg)), reference_locus_values = lengths(locus),
        reference_kegg_values = lengths(kegg), stringsAsFactors = FALSE, row.names = NULL)
}

#' The reference a mapping is made against: a curated name or UPID, or a user-supplied list
#' @param reference A curated reference's name or UniProt proteome id, or
#'   `list(fasta, xref, name, uniprot_taxon, string_species, kegg_code)` (with optional `locus_pattern`, `release`).
#' @return `list(name, upid, release, fasta, xref, go, uniprot_taxon, string_species, string_status, kegg_code,
#'   locus_pattern, supplied)`.
#' @noRd
resolveMappingReference <- function(reference) {
    if (is.character(reference) && length(reference) == 1L) {
        organisms <- readReferenceOrganisms()
        row <- organisms[organisms$uniprot_upid == reference | organisms$name == reference, , drop = FALSE]
        if (nrow(row) != 1L) {
            referenceMappingRefuse("reference_mapping_reference_invalid",
                sprintf("'%s' is not one of the shipped references (%s)", reference,
                    paste(organisms$name, collapse = "; ")))
        }
        files <- referenceOrganismFiles(row$uniprot_upid)
        return(list(name = row$name, upid = row$uniprot_upid, release = row$uniprot_release, fasta = files$fasta,
            xref = files$xref, go = files$go, uniprot_taxon = row$uniprot_taxon, string_species = row$string_species,
            string_status = row$string_status, kegg_code = row$kegg_code, locus_pattern = row$locus_pattern,
            supplied = FALSE))
    }
    if (!is.list(reference)) {
        referenceMappingRefuse("reference_mapping_reference_invalid", "it is neither a shipped reference nor a list")
    }
    needed <- c("fasta", "xref", "name", "uniprot_taxon", "string_species", "kegg_code")
    missing <- setdiff(needed, names(reference))
    if (length(missing)) {
        referenceMappingRefuse("reference_mapping_reference_invalid",
            sprintf("the supplied reference lacks %s", paste(missing, collapse = ", ")))
    }
    problems <- c(
        if (!is.character(reference$fasta) || !file.exists(reference$fasta)) "its FASTA file does not exist",
        if (!is.character(reference$xref) || !file.exists(reference$xref)) "its cross-reference table does not exist",
        if (!grepl("^[0-9]+$", as.character(reference$uniprot_taxon))) "its UniProt taxon is not a number",
        if (!is.na(reference$string_species) && !grepl("^[0-9]+$", as.character(reference$string_species)))
            "its STRING species is not a number",
        if (!grepl("^[a-z]{3,4}$", as.character(reference$kegg_code))) "its KEGG code is not 3 or 4 lower-case letters")
    if (length(problems)) {
        referenceMappingRefuse("reference_mapping_reference_invalid", paste(problems, collapse = "; "))
    }
    list(name = as.character(reference$name), upid = reference$upid %||% "supplied",
        release = reference$release %||% "supplied", fasta = reference$fasta, xref = reference$xref,
        go = reference$go %||% NA_character_, uniprot_taxon = as.integer(reference$uniprot_taxon),
        string_species = as.character(reference$string_species),
        string_status = if (is.na(reference$string_species)) "unavailable" else "supplied",
        kegg_code = reference$kegg_code, locus_pattern = reference$locus_pattern %||% "^$", supplied = TRUE)
}
