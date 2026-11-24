# KIT-167: map an individually sequenced strain's proteins to their orthologues in one well-annotated reference of the
# same species, by reciprocal best hits of DIAMOND searches in both directions (sections 5.1 to 5.4).
#   - `mapStrainToReference()` (exported): strain FASTA and one reference -> one row per strain protein with
#     `map_status` and `reference_accession_used` (an accession for `one_to_one` only), counts, the identity summary
#     and the provenance that is the mapping's identity (section 7).
#   - `readStrainProteinFasta()`, `classifyStrainHeaders()`, `prepareStrainEntities()`: the strain entries by header
#     class and the input refusals of section 8.4.
#   - `collapseIdenticalSequences()`: one search entity per distinct sequence (5.4).
#   - `resolveReciprocalBestHits()`: the status rules of 5.2 and 5.3 on two hit tables, in R.
#   - `assembleReferenceMapping()`, `summariseReferenceMapping()`, `canonicalReferenceMappingText()`.

#' The six statuses of a strain protein, in the order the rules are evaluated (section 5.3)
#' @noRd
referenceMappingStatuses <- c("no_hit", "below_threshold", "tie", "one_to_one", "many_to_one", "one_way")

#' The version of the resolver's rules, part of a mapping's identity (section 7)
#' @noRd
referenceMappingResolverVersion <- "kit167-rbh-1"

#' The ranking key (section 5.2): the raw `score` if gate F2 shows it is an integer, else the bit score as printed.
#' @noRd
referenceMappingRankingKey <- "score"

#' The strain header classes (section 8.4)
#' @noRd
referenceMappingHeaderClasses <- c("refseq", "ncbi_cds", "plain", "uniprot")

#' The residues a protein sequence may hold: the twenty amino acids and X U O B Z J (section 8.4)
#' @noRd
referenceMappingAlphabet <- "^[ACDEFGHIKLMNPQRSTVWYXUOBZJ]*$"

#' Read a protein FASTA file with base R (plain or gzip)
#' @param path The file.
#' @return Data frame: header (without `>`), sequence (as written, lines joined).
#' @noRd
readStrainProteinFasta <- function(path) {
    if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
        referenceMappingRefuse("reference_mapping_too_few", 0L)
    }
    connection <- file(path, open = "r")
    on.exit(close(connection), add = TRUE)
    lines <- sub("\r$", "", readLines(connection, warn = FALSE))
    is_header <- startsWith(lines, ">")
    entry <- cumsum(is_header)
    keep <- !is_header & entry > 0L & nzchar(lines)
    sequences <- character(sum(is_header))
    if (any(keep)) {
        joined <- tapply(gsub("[[:space:]]", "", lines[keep]), entry[keep], paste, collapse = "")
        sequences[as.integer(names(joined))] <- unname(joined)
    }
    data.frame(header = substring(lines[is_header], 2L), sequence = sequences, stringsAsFactors = FALSE)
}

#' Class each FASTA header and select the strain entries (section 8.4)
#'
#' `>ref|WP_...|` (accession the WP_), NCBI translated CDS `>lcl|... [protein_id=...]` (accession the protein_id),
#' plain `>ACCESSION ...`, UniProt `>sp|`/`>tr|` (accession the second field). MaxQuant contaminant entries are excluded
#' by class. UniProt-style entries are strain entries only when `strain_taxon` is given and equals their `OX`, or when
#' the file has no entry of another strain class and exactly one `OX`; otherwise they are excluded and counted per `OX`.
#' A file of UniProt-style entries only, with several `OX` and no `strain_taxon`, is refused.
#' @param entries `readStrainProteinFasta()`'s table.
#' @param header_classes The classes that may hold strain entries.
#' @param strain_taxon The strain's NCBI taxon, or NULL.
#' @return `list(entries, counts)`: the entries with `class`, `accession`, `locus_tag`, `ox`, `selected`; counts per
#'   class (`by_class`) and the UniProt-style entries excluded per `OX` (`excluded_by_ox`).
#' @noRd
classifyStrainHeaders <- function(entries, header_classes = referenceMappingHeaderClasses, strain_taxon = NULL) {
    header <- entries$header
    first <- sub("[[:space:]].*$", "", header)
    tag <- "(SWISS-PROT|TREMBL|REFSEQ|ENSEMBL|H-INV):"
    contaminant <- grepl("^(CON__|CONTAMINANT_|cRAP-)", first) | grepl(paste0("^", tag), first) |
        grepl(paste0("^[^[:space:]]+[[:space:]]+", tag), header) | grepl("^Streptavidin", header)
    refseq <- grepl("^ref\\|[^|]+\\|", header)
    cds <- grepl("^lcl\\|", header) & grepl("\\[protein_id=[^]]+\\]", header)
    uniprot <- grepl("^(sp|tr)\\|[^|]+\\|", header)
    plain <- nzchar(first) & !grepl("[|:]", first)
    class <- ifelse(contaminant, "contaminant", ifelse(refseq, "refseq", ifelse(cds, "ncbi_cds",
        ifelse(uniprot, "uniprot", ifelse(plain, "plain", "other")))))
    accession <- rep(NA_character_, length(header))
    accession[class == "refseq"] <- sub("^ref\\|([^|]+)\\|.*$", "\\1", header[class == "refseq"])
    accession[class == "ncbi_cds"] <- sub("^.*\\[protein_id=([^]]+)\\].*$", "\\1", header[class == "ncbi_cds"])
    accession[class == "uniprot"] <- sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", header[class == "uniprot"])
    accession[class == "plain"] <- first[class == "plain"]
    locus_tag <- ifelse(grepl("\\[locus_tag=[^]]+\\]", header), sub("^.*\\[locus_tag=([^]]+)\\].*$", "\\1", header),
        NA_character_)
    ox <- ifelse(class == "uniprot" & grepl("[[:space:]]OX=[0-9]+", header),
        sub("^.*[[:space:]]OX=([0-9]+).*$", "\\1", header), NA_character_)
    strain_class <- class %in% intersect(header_classes, c("refseq", "ncbi_cds", "plain"))
    candidate <- class == "uniprot" & "uniprot" %in% header_classes
    taxa <- sort(unique(ox[candidate & !is.na(ox)]), method = "radix")
    chosen <- if (!is.null(strain_taxon)) {
        candidate & !is.na(ox) & ox == as.character(strain_taxon)
    } else if (!any(strain_class) && length(taxa) == 1L) {
        candidate & !is.na(ox) & ox == taxa
    } else if (!any(strain_class) && length(taxa) > 1L) {
        per_ox <- table(ox[candidate & !is.na(ox)])
        referenceMappingRefuse("reference_mapping_taxon_ambiguous",
            paste(sprintf("OX=%s: %d entries", names(per_ox), as.integer(per_ox)), collapse = "; "),
            fields = list(entries_per_ox = as.list(c(per_ox))))
    } else {
        rep(FALSE, length(header))
    }
    excluded <- candidate & !chosen
    excluded_ox <- ifelse(is.na(ox[excluded]), "none", ox[excluded])
    entries$class <- class
    entries$accession <- accession
    entries$locus_tag <- locus_tag
    entries$ox <- ox
    entries$selected <- strain_class | chosen
    classes <- c("refseq", "ncbi_cds", "plain", "uniprot", "contaminant", "other")
    counts <- list(entries = nrow(entries),
        by_class = lapply(stats::setNames(classes, classes), function(name) sum(class == name)),
        excluded_by_ox = as.list(lengths(split(excluded_ox, excluded_ox))), selected = sum(entries$selected))
    list(entries = entries, counts = counts)
}

#' One search entity per distinct sequence (section 5.4)
#'
#' Key: SHA-256 of the upper-case sequence without a terminal `*`. The entity's id is its members' first accession in
#' C-locale order; every member inherits the entity's status.
#' @param entries Data frame with `accession` and `sequence`.
#' @return `list(entities, members, identical_sequence_groups)`: entities (entity_id, sequence, key, n_members),
#'   members (strain_accession, entity_id).
#' @noRd
collapseIdenticalSequences <- function(entries) {
    sequence <- sub("\\*+$", "", toupper(entries$sequence))
    key <- vapply(sequence, digest::digest, character(1), algo = "sha256", serialize = FALSE, USE.NAMES = FALSE)
    members <- unique(data.frame(strain_accession = entries$accession, key = key, stringsAsFactors = FALSE))
    by_key <- order(members$key, members$strain_accession, method = "radix")
    first <- by_key[!duplicated(members$key[by_key])]
    entity_of_key <- stats::setNames(members$strain_accession[first], members$key[first])
    members$entity_id <- unname(entity_of_key[members$key])
    sizes <- table(members$entity_id)
    entities <- data.frame(entity_id = unname(entity_of_key), sequence = sequence[match(names(entity_of_key), key)],
        key = names(entity_of_key), stringsAsFactors = FALSE)
    entities$n_members <- as.integer(sizes[entities$entity_id])
    entities <- entities[order(entities$entity_id, method = "radix"), , drop = FALSE]
    members <- members[order(members$strain_accession, method = "radix"), c("strain_accession", "entity_id")]
    rownames(entities) <- NULL
    rownames(members) <- NULL
    list(entities = entities, members = members, identical_sequence_groups = sum(entities$n_members > 1L))
}

#' The strain entries of a FASTA file, validated and collapsed (section 8.4)
#' @param path The strain protein FASTA.
#' @param header_classes,strain_taxon As `classifyStrainHeaders()`.
#' @param min_entities The fewest distinct strain proteins accepted (50).
#' @return `list(entities, members, identical_sequence_groups, counts, selection_sha256)`.
#' @noRd
prepareStrainEntities <- function(path, header_classes = referenceMappingHeaderClasses, strain_taxon = NULL,
                                  min_entities = 50L) {
    classified <- classifyStrainHeaders(readStrainProteinFasta(path), header_classes, strain_taxon)
    selected <- classified$entries[classified$entries$selected, c("accession", "sequence"), drop = FALSE]
    selected$sequence <- sub("\\*$", "", toupper(selected$sequence))
    stopped <- grepl("*", selected$sequence, fixed = TRUE)
    foreign <- !grepl(referenceMappingAlphabet, selected$sequence)
    if (any(stopped | foreign)) {
        bad <- selected$accession[stopped | foreign]
        referenceMappingRefuse("reference_mapping_not_protein",
            sprintf("%d entries with an internal stop or a character outside the protein alphabet, for example %s",
                length(bad), paste(utils::head(bad, 5L), collapse = ", ")), fields = list(accessions = bad))
    }
    residues <- sum(nchar(selected$sequence))
    nucleotide <- sum(nchar(gsub("[^ACGTN]", "", selected$sequence)))
    if (residues > 0 && nucleotide / residues > 0.9) {
        referenceMappingRefuse("reference_mapping_not_protein",
            sprintf("%.0f%% of the residues are A, C, G, T or N", 100 * nucleotide / residues))
    }
    distinct <- unique(selected)
    duplicated_ids <- sort(unique(distinct$accession[duplicated(distinct$accession)]), method = "radix")
    if (length(duplicated_ids)) {
        referenceMappingRefuse("reference_mapping_duplicate_id", paste(utils::head(duplicated_ids, 20L), collapse = ", "),
            fields = list(accessions = duplicated_ids))
    }
    collapsed <- collapseIdenticalSequences(distinct)
    if (nrow(collapsed$entities) < min_entities) {
        referenceMappingRefuse("reference_mapping_too_few", nrow(collapsed$entities))
    }
    ordered <- distinct[order(distinct$accession, method = "radix"), , drop = FALSE]
    selection_text <- paste0(paste(ordered$accession, ordered$sequence, sep = "\t"), "\n", collapse = "")
    c(collapsed, list(counts = classified$counts,
        selection_sha256 = digest::digest(selection_text, algo = "sha256", serialize = FALSE)))
}

#' Best sets of one search direction (section 5.2)
#' @param hits The direction's hit table.
#' @param coverage,evalue,min_identity,ranking_key The settings.
#' @return `list(best, top, max_all, max_candidate)`: per query, the best candidate set (C-locale order; empty when
#'   the top hit is filtered), the top hit row, and the maximum keys.
#' @noRd
referenceMappingBestSets <- function(hits, coverage, evalue, min_identity, ranking_key) {
    hits <- hits[!is.na(hits$evalue) & hits$evalue <= evalue, , drop = FALSE]
    hits$key <- as.numeric(hits[[ranking_key]])
    hits$candidate <- hits$qcovhsp >= coverage & hits$scovhsp >= coverage &
        (if (is.null(min_identity)) TRUE else hits$pident >= min_identity)
    if (!nrow(hits)) return(list(best = list(), top = hits, max_all = numeric(), max_candidate = numeric()))
    max_all <- tapply(hits$key, hits$qseqid, max)
    candidates <- hits[hits$candidate, , drop = FALSE]
    max_candidate <- if (nrow(candidates)) tapply(candidates$key, candidates$qseqid, max) else numeric()
    best_rows <- candidates[candidates$key == max_candidate[candidates$qseqid], , drop = FALSE]
    best <- lapply(split(best_rows$sseqid, best_rows$qseqid), sort, method = "radix")
    filtered <- names(max_candidate)[max_all[names(max_candidate)] > max_candidate]
    best[filtered] <- NULL
    top_order <- order(hits$qseqid, -hits$key, !hits$candidate, hits$sseqid, method = "radix")
    top <- hits[top_order, , drop = FALSE]
    top <- top[!duplicated(top$qseqid), , drop = FALSE]
    list(best = best, top = top, max_all = max_all, max_candidate = max_candidate)
}

#' Resolve reciprocal best hits (sections 5.2 and 5.3)
#'
#' A candidate is a hit with e <= `evalue`, both coverages >= `coverage` and, when `min_identity` is set, identity >=
#' it, in both directions alike. The top hit is ranked by the key before coverage; when it is not a candidate and its
#' key is strictly greater than the best candidate's, the best set is empty (`top_hit_filtered`). Ties are exact
#' equality of the key. Statuses: `no_hit`, `below_threshold`, `tie`, `one_to_one`, `many_to_one`, `one_way`.
#' @param forward Hits of the strain entities against the reference (`diamondOutputFields`).
#' @param reverse Hits of the reference proteins against the strain entities.
#' @param entity_ids Every strain search entity.
#' @param coverage Percent of both sequences an HSP must cover (80).
#' @param evalue E-value (1e-5).
#' @param min_identity NULL (no floor) or a percent identity a candidate needs in both directions.
#' @param ranking_key `"score"` or `"bitscore"`.
#' @return Data frame, one row per entity in the order given, with `map_status`, `status_detail`,
#'   `reference_accession`, `reference_accession_used`, the reported hit's statistics and the top-hit columns;
#'   attribute `one_direction_pairs` (pairs whose e-value passes in one direction only).
#' @noRd
resolveReciprocalBestHits <- function(forward, reverse, entity_ids, coverage = 80, evalue = 1e-5, min_identity = NULL,
                                      ranking_key = referenceMappingRankingKey) {
    ranking_key <- match.arg(ranking_key, c("score", "bitscore"))
    fwd <- referenceMappingBestSets(forward, coverage, evalue, min_identity, ranking_key)
    rev <- referenceMappingBestSets(reverse, coverage, evalue, min_identity, ranking_key)
    n <- length(entity_ids)
    status <- character(n)
    detail <- rep(NA_character_, n)
    reference <- rep(NA_character_, n)
    partner <- rep(NA_character_, n)
    pending <- logical(n)
    for (i in seq_len(n)) {
        s <- entity_ids[[i]]
        best <- fwd$best[[s]]
        if (is.na(fwd$max_all[s])) {
            status[i] <- "no_hit"
        } else if (!length(best)) {
            status[i] <- "below_threshold"
            detail[i] <- if (is.na(fwd$max_candidate[s])) "no_candidate" else "top_hit_filtered"
        } else if (length(best) > 1L) {
            status[i] <- "tie"
            detail[i] <- "forward_tie"
            reference[i] <- paste(best, collapse = ";")
        } else {
            reference[i] <- best
            back <- rev$best[[best]]
            if (identical(back, s)) {
                status[i] <- "one_to_one"
                partner[i] <- best
            } else if (length(back) > 1L && s %in% back) {
                status[i] <- "tie"
                detail[i] <- "reverse_tie"
            } else {
                pending[i] <- TRUE
            }
        }
    }
    partnered <- partner[!is.na(partner)]
    status[pending] <- ifelse(reference[pending] %in% partnered, "many_to_one", "one_way")
    resolved <- data.frame(entity_id = entity_ids, reference_accession = reference, reference_accession_used = partner,
        map_status = status, status_detail = detail, stringsAsFactors = FALSE)
    hits <- forward[!is.na(forward$evalue) & forward$evalue <= evalue, , drop = FALSE]
    reported <- ifelse(status %in% c("one_to_one", "many_to_one", "one_way") |
        (status == "tie" & detail == "reverse_tie"), reference, NA_character_)
    top_index <- match(entity_ids, fwd$top$qseqid)
    reported[status == "below_threshold"] <- fwd$top$sseqid[top_index[status == "below_threshold"]]
    first_tied <- sub(";.*$", "", reference)
    reported[status == "tie" & detail == "forward_tie"] <- first_tied[status == "tie" & detail == "forward_tie"]
    row <- match(paste(entity_ids, reported, sep = "\r"), paste(hits$qseqid, hits$sseqid, sep = "\r"))
    for (field in c("pident", "length", "qcovhsp", "scovhsp", "evalue", "bitscore", "score", "qlen", "slen")) {
        resolved[[field]] <- hits[[field]][row]
    }
    resolved$top_hit_accession <- fwd$top$sseqid[top_index]
    resolved$top_hit_score <- fwd$top$key[top_index]
    resolved$top_hit_is_candidate <- fwd$top$candidate[top_index]
    pairs <- function(h, flip) {
        h <- h[!is.na(h$evalue) & h$evalue <= evalue, , drop = FALSE]
        unique(if (flip) paste(h$sseqid, h$qseqid, sep = "\r") else paste(h$qseqid, h$sseqid, sep = "\r"))
    }
    forward_pairs <- pairs(forward, FALSE)
    reverse_pairs <- pairs(reverse, TRUE)
    attr(resolved, "one_direction_pairs") <- length(setdiff(forward_pairs, reverse_pairs)) +
        length(setdiff(reverse_pairs, forward_pairs))
    resolved
}

#' One row per strain protein: the members of each entity with the entity's status and the reference's cross-references
#' @param resolved `resolveReciprocalBestHits()`'s table.
#' @param members Data frame (strain_accession, entity_id).
#' @param xref NULL, or `parseReferenceCrossReferences()`'s table.
#' @return The mapping table of section 8.2, ordered by strain accession (C locale).
#' @noRd
assembleReferenceMapping <- function(resolved, members, xref = NULL) {
    at <- match(members$entity_id, resolved$entity_id)
    mapping <- data.frame(strain_accession = members$strain_accession, entity_id = members$entity_id,
        reference_accession = resolved$reference_accession[at],
        reference_accession_used = resolved$reference_accession_used[at], stringsAsFactors = FALSE)
    used <- if (is.null(xref)) rep(NA_integer_, nrow(mapping)) else match(mapping$reference_accession_used, xref$accession)
    for (field in c("reference_gene", "reference_locus", "reference_string_id", "reference_kegg_id")) {
        mapping[[field]] <- if (is.null(xref)) rep(NA_character_, nrow(mapping)) else xref[[field]][used]
    }
    for (field in c("pident", "length", "qcovhsp", "scovhsp", "evalue", "bitscore", "score", "qlen", "slen",
        "top_hit_accession", "top_hit_score", "top_hit_is_candidate", "map_status", "status_detail")) {
        mapping[[field]] <- resolved[[field]][at]
    }
    mapping <- mapping[order(mapping$strain_accession, method = "radix"), , drop = FALSE]
    rownames(mapping) <- NULL
    mapping
}

#' Counts by status and the identity summary of a mapping (sections 5.3 and 5.6)
#'
#' Identity over `one_to_one` pairs (one per entity): median, 5th percentile (`quantile` type 7) and share below 40 %;
#' `caution` when the median is below 95 % (VAL-38).
#' @param mapping The mapping table.
#' @return `list(counts, identity_summary)`.
#' @noRd
summariseReferenceMapping <- function(mapping) {
    by_status <- lapply(stats::setNames(referenceMappingStatuses, referenceMappingStatuses),
        function(s) sum(mapping$map_status == s))
    group_sizes <- table(mapping$entity_id)
    pairs <- mapping[mapping$map_status == "one_to_one" & !duplicated(mapping$entity_id), , drop = FALSE]
    identity <- as.numeric(pairs$pident)
    median_identity <- if (length(identity)) stats::median(identity) else NA_real_
    list(counts = list(strain_proteins = nrow(mapping), entities = length(group_sizes),
            identical_sequence_groups = sum(group_sizes > 1L), by_status = by_status,
            used_with_string = sum(!is.na(mapping$reference_string_id) & !is.na(mapping$reference_accession_used)),
            used_with_kegg = sum(!is.na(mapping$reference_kegg_id) & !is.na(mapping$reference_accession_used))),
        identity_summary = list(pairs = length(identity), median = median_identity,
            p05 = if (length(identity)) unname(stats::quantile(identity, 0.05, type = 7)) else NA_real_,
            share_below_40 = if (length(identity)) mean(identity < 40) else NA_real_,
            caution = isTRUE(median_identity < 95)))
}

#' The canonical text of a table (section 7): rows in the given order, statistics printed with 17 significant digits,
#' UTF-8 and `\n`. Its SHA-256 is the table's identity.
#' @param table A data frame (the mapping or the gene table).
#' @param order_by The column the rows are ordered by (C locale, radix), or NULL to keep the order.
#' @return One string.
#' @noRd
canonicalReferenceMappingText <- function(table, order_by = "strain_accession") {
    if (!is.null(order_by)) table <- table[order(table[[order_by]], method = "radix"), , drop = FALSE]
    columns <- lapply(table, function(column) {
        text <- if (is.double(column)) formatC(column, digits = 17, format = "g") else as.character(column)
        text[is.na(column)] <- "NA"
        enc2utf8(text)
    })
    body <- if (nrow(table)) do.call(paste, c(columns, sep = "\t")) else character()
    paste0(paste(c(paste(names(table), collapse = "\t"), body), collapse = "\n"), "\n")
}

#' SHA-256 of a table's canonical text
#' @noRd
canonicalReferenceMappingSha256 <- function(table, order_by = "strain_accession") {
    digest::digest(canonicalReferenceMappingText(table, order_by), algo = "sha256", serialize = FALSE)
}

#' Read a reference proteome FASTA (UniProt headers) for the search
#' @noRd
readReferenceProteome <- function(path) {
    entries <- readStrainProteinFasta(path)
    uniprot <- grepl("^(sp|tr)\\|[^|]+\\|", entries$header)
    if (!nrow(entries) || !all(uniprot)) {
        referenceMappingRefuse("reference_mapping_reference_invalid",
            sprintf("%d of %d FASTA headers are not UniProt headers (>sp| or >tr|)", sum(!uniprot), nrow(entries)))
    }
    entries$accession <- sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", entries$header)
    if (anyDuplicated(entries$accession)) {
        referenceMappingRefuse("reference_mapping_reference_invalid", "its FASTA repeats an accession")
    }
    entries$sequence <- sub("\\*+$", "", toupper(entries$sequence))
    entries
}

#' Map a strain's proteins to their orthologues in one reference organism
#'
#' DIAMOND v2.2.8 searches the strain against the reference and the reference against the strain at the shipped
#' settings (`--ultra-sensitive --masking none --ext full --comp-based-stats 1 -k0 --max-hsps 1`), and reciprocal
#' best hits are resolved in R. A strain protein is paired with a reference protein (`one_to_one`) only when each is
#' the other's best candidate; only those pairs carry `reference_accession_used`. Validated for a strain mapped to a
#' reference of the same species; for any other reference the mapping is reported but is not covered by the
#' validation.
#'
#' @param strain_fasta Path to the strain's protein FASTA (RefSeq `>ref|WP_...|`, NCBI translated CDS, plain or
#'   UniProt headers; MaxQuant contaminants and other organisms' UniProt entries are left out).
#' @param reference A shipped reference (its name or UniProt proteome id, see `inst/extdata/reference_organisms.tsv`),
#'   or a list `list(fasta, xref, name, uniprot_taxon, string_species, kegg_code)` for a UniProt proteome the user
#'   supplies.
#' @param coverage Percent of both the strain and the reference protein one alignment must cover (default 80).
#' @param evalue E-value threshold (default 1e-5).
#' @param header_classes FASTA header classes that may hold strain entries.
#' @param strain_taxon The strain's NCBI taxon, needed only when UniProt-style entries of several organisms are given.
#' @param min_identity NULL (default, no identity floor) or a percent identity both directions must reach.
#' @param threads DIAMOND threads, or NULL for its default (results do not depend on it).
#' @param scratch A folder for the search files, or NULL for a temporary one; removed after the run.
#' @return `list(mapping, counts, identity_summary, provenance)`: `mapping` has one row per strain protein (columns
#'   of the ticket's section 8.2); `counts` per status (summing to the strain proteins); `identity_summary` over the
#'   one-to-one pairs; `provenance` the mapping's identity (inputs' SHA-256, engine, every argument, settings).
#' @export
mapStrainToReference <- function(strain_fasta, reference, coverage = 80, evalue = 1e-5,
                                 header_classes = c("refseq", "ncbi_cds", "plain", "uniprot"), strain_taxon = NULL,
                                 min_identity = NULL, threads = NULL, scratch = NULL) {
    checkmate::assert_number(coverage, lower = 0, upper = 100)
    checkmate::assert_number(evalue, lower = 0)
    checkmate::assert_subset(header_classes, referenceMappingHeaderClasses, empty.ok = FALSE)
    checkmate::assert_number(min_identity, lower = 0, upper = 100, null.ok = TRUE)
    checkmate::assert_count(threads, positive = TRUE, null.ok = TRUE)
    started <- Sys.time()
    ref <- resolveMappingReference(reference)
    strain <- prepareStrainEntities(strain_fasta, header_classes, strain_taxon)
    proteome <- readReferenceProteome(ref$fasta)
    xref <- parseReferenceCrossReferences(readReferenceCrossReferenceTable(ref$xref), ref$kegg_code,
        ref$string_species, ref$locus_pattern)
    engine <- checkDiamond()
    folder <- scratch %||% tempfile("reference-mapping-")
    dir.create(file.path(folder, "tmp"), recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(folder, recursive = TRUE), add = TRUE)
    strain_ids <- sprintf("s%07d", seq_len(nrow(strain$entities)))
    reference_ids <- sprintf("r%07d", seq_len(nrow(proteome)))
    writeLines(paste0(">", strain_ids, "\n", strain$entities$sequence), file.path(folder, "strain.fasta"))
    writeLines(paste0(">", reference_ids, "\n", proteome$sequence), file.path(folder, "reference.fasta"))
    left <- function() max(1, 600 - as.numeric(difftime(Sys.time(), started, units = "secs")))
    diamondCommands <- function(threads) {
        thread_args <- if (!is.null(threads)) c("--threads", as.character(threads))
        list(makedb_reference = c("makedb", "--in", "reference.fasta", "--db", "reference", thread_args),
            makedb_strain = c("makedb", "--in", "strain.fasta", "--db", "strain", thread_args),
            forward = diamondSearchArguments("strain.fasta", "reference", "forward.tsv", evalue, threads),
            reverse = diamondSearchArguments("reference.fasta", "strain", "reverse.tsv", evalue, threads))
    }
    for (command in diamondCommands(threads)) runDiamondCommand(engine$path, command, scratch = folder, timeout = left())
    forward <- readDiamondHits(file.path(folder, "forward.tsv"))
    reverse <- readDiamondHits(file.path(folder, "reverse.tsv"))
    to_strain <- stats::setNames(strain$entities$entity_id, strain_ids)
    to_reference <- stats::setNames(proteome$accession, reference_ids)
    forward$qseqid <- unname(to_strain[forward$qseqid])
    forward$sseqid <- unname(to_reference[forward$sseqid])
    reverse$qseqid <- unname(to_reference[reverse$qseqid])
    reverse$sseqid <- unname(to_strain[reverse$sseqid])
    resolved <- resolveReciprocalBestHits(forward, reverse, strain$entities$entity_id, coverage = coverage,
        evalue = evalue, min_identity = min_identity, ranking_key = referenceMappingRankingKey)
    mapping <- assembleReferenceMapping(resolved, strain$members, xref)
    summary <- summariseReferenceMapping(mapping)
    summary$counts$one_direction_pairs <- attr(resolved, "one_direction_pairs")
    summary$counts$strain_file <- strain$counts
    summary$counts$reference_proteins <- nrow(proteome)
    summary$counts$used_with_several_kegg <- sum(xref$reference_kegg_values[match(mapping$reference_accession_used,
        xref$accession)] > 1L, na.rm = TRUE)
    provenance <- list(
        strain_fasta_sha256 = digest::digest(file = strain_fasta, algo = "sha256"),
        strain_fasta_bytes = file.size(strain_fasta), strain_selection_sha256 = strain$selection_sha256,
        reference_name = ref$name, reference_upid = ref$upid, reference_release = ref$release,
        reference_supplied = ref$supplied, reference_uniprot_taxon = ref$uniprot_taxon,
        reference_string_species = ref$string_species, reference_string_status = ref$string_status,
        reference_kegg_code = ref$kegg_code,
        reference_fasta_sha256 = digest::digest(file = ref$fasta, algo = "sha256"),
        reference_xref_sha256 = digest::digest(file = ref$xref, algo = "sha256"),
        engine = "DIAMOND", engine_version = engine$version, engine_sha256 = engine$sha256,
        engine_platform = engine$platform,
        engine_arguments = diamondCommands(NULL),
        coverage = coverage, evalue = evalue, min_identity = min_identity, header_classes = header_classes,
        strain_taxon = strain_taxon, ranking_key = referenceMappingRankingKey,
        resolver_version = referenceMappingResolverVersion,
        mapping_sha256 = canonicalReferenceMappingSha256(mapping))
    list(mapping = mapping, counts = summary$counts, identity_summary = summary$identity_summary,
        provenance = provenance)
}
