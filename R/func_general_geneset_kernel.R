# KIT-173: Gene Ontology parent terms for gene sets built from direct annotations.
#   - `propagateGoAnnotations()`: one aspect's direct (term, gene) pairs propagated to every ancestor term over a named
#     relation set read from GO.db's parent maps. The shared function KIT-168 and KIT-181 reuse.
#   - `buildUniprotGoGeneSets()`: the proteomics clusterProfiler route's gene sets from the UniProt annotation table,
#     direct (`"none"`, today's tables exactly) or propagated with the aspect roots removed. KIT-179 and KIT-168's
#     proteomics adapter reuse it.
#   - `goIdsNotInGoDb()`: the raw GO ids the installed GO.db does not hold, counted by class where raw ids first meet
#     GO.db (`getUniprotAnnotations`).
# The names below are the ticket's contract with KIT-167, KIT-168, KIT-170, KIT-179 and KIT-181 (KIT-173 section 14),
# so the constants keep their camelCase names.

#' GO.db relation labels followed for each propagating `go_propagation` value (KIT-173, VAL-28)
#'
#' The labels are those GO.db 3.23.1's parent maps (`GOBPPARENTS`, `GOMFPARENTS`, `GOCCPARENTS`) return (KIT-173 F1).
#' `is_a_part_of_regulates` is the set g:Profiler and GO.db's ancestor tables use; `is_a_part_of` is the GO
#' Consortium's set that is safe for grouping annotations.
#' @noRd
goPropagationRelations <- list(
    is_a_part_of = c("isa", "part of"),
    is_a_part_of_regulates = c("isa", "part of", "regulates", "positively regulates", "negatively regulates")
)

#' The shipped `go_propagation` value (KIT-173, VAL-28): the five relations.
#' @noRd
goPropagationShipped <- "is_a_part_of_regulates"

#' The three aspect roots (biological process, molecular function, cellular component).
#' @noRd
goAspectRoots <- c(BP = "GO:0008150", MF = "GO:0003674", CC = "GO:0005575")

# GO.db reads, made once per R session (the parent maps per aspect, the term aspects, the release).
goDbCache <- new.env(parent = emptyenv())

#' One aspect's parent edges from GO.db, read once per session
#' @param ontology `"BP"`, `"MF"` or `"CC"`.
#' @return Data frame `child`, `parent`, `relation` (the parent map's relation labels).
#' @noRd
goParentTable <- function(ontology) {
    key <- paste0("parents_", ontology)
    if (is.null(goDbCache[[key]])) {
        map <- switch(ontology, BP = GO.db::GOBPPARENTS, MF = GO.db::GOMFPARENTS, CC = GO.db::GOCCPARENTS)
        edges <- AnnotationDbi::toTable(map)
        goDbCache[[key]] <- data.frame(child = edges[[1L]], parent = edges[[2L]], relation = edges[[3L]],
            stringsAsFactors = FALSE)
    }
    goDbCache[[key]]
}

#' The aspect of every term GO.db holds, read once per session
#' @return Named character vector: GO id to `"BP"`, `"MF"` or `"CC"`.
#' @noRd
goDbTermAspects <- function() {
    if (is.null(goDbCache$aspects)) {
        aspects <- AnnotationDbi::Ontology(GO.db::GOTERM)
        goDbCache$aspects <- aspects[aspects %in% c("BP", "MF", "CC")]
    }
    goDbCache$aspects
}

#' The installed GO.db's version and the date of the ontology it was built from
#' @return `list(version, date)`.
#' @noRd
goDbRelease <- function() {
    if (is.null(goDbCache$release)) {
        info <- GO.db::GO_dbInfo()
        goDbCache$release <- list(version = as.character(utils::packageVersion("GO.db")),
            date = info$value[match("GOSOURCEDATE", info$name)])
    }
    goDbCache$release
}

#' Every (term, ancestor) pair reachable over the relations, each term its own first entry
#' @param terms GO ids of one aspect.
#' @param parents `goParentTable()` of that aspect.
#' @param relations Relation labels to follow.
#' @return Data frame `term`, `node`: each term with itself and each ancestor; the root `all` is never a node.
#' @noRd
goAncestorPairs <- function(terms, parents, relations) {
    edges <- parents[parents$relation %in% relations & parents$parent != "all", c("child", "parent")]
    seen <- data.frame(term = terms, node = terms, stringsAsFactors = FALSE)
    seen_keys <- paste(seen$term, seen$node)
    frontier <- seen
    while (nrow(frontier) > 0L) {
        step <- merge(frontier, edges, by.x = "node", by.y = "child")
        found <- data.frame(term = step$term, node = step$parent, stringsAsFactors = FALSE)
        keys <- paste(found$term, found$node)
        fresh <- !duplicated(keys) & !keys %in% seen_keys
        frontier <- found[fresh, , drop = FALSE]
        seen <- rbind(seen, frontier)
        seen_keys <- c(seen_keys, keys[fresh])
    }
    seen
}

#' Propagate one aspect's direct GO annotations to their parent terms (KIT-173)
#'
#' Applies the true-path rule: a gene annotated to a term is annotated to every ancestor of that term reachable over
#' the named relations, within the term's aspect, in the installed GO.db. The root `all` is never a term; the aspect
#' roots are kept (a caller removes them by its own rule). Input terms GO.db does not hold, or that belong to another
#' aspect, are dropped and counted.
#' @param term2gene_direct Data frame with character columns `TERM` (GO ids) and `GENE` (any identifier): the direct
#'   annotations of one aspect.
#' @param ontology One of `"BP"`, `"MF"`, `"CC"`.
#' @param relations Non-empty character vector of GO.db relation labels (`goPropagationRelations[[value]]`); no default.
#' @return Data frame `TERM`, `GENE`, distinct, sorted by `TERM` then `GENE`, with attributes `relations`,
#'   `go_db_version`, `go_source_date`, `not_in_go_db` (count of distinct ids), `not_in_go_db_examples` (up to 20 ids),
#'   `not_in_aspect` (count of distinct ids), `pairs_direct` (distinct direct pairs used) and `pairs_propagated`.
#' @noRd
propagateGoAnnotations <- function(term2gene_direct, ontology, relations) {
    requireOptionalPackage("GO.db", "GO parent terms")
    if (!is.character(ontology) || length(ontology) != 1L || !ontology %in% names(goAspectRoots)) {
        cli::cli_abort("{.arg ontology} must be one of BP, MF or CC, not {.val {ontology}}.",
            class = "go_propagation_invalid")
    }
    labels <- unique(unlist(goPropagationRelations, use.names = FALSE))
    if (!is.character(relations) || !length(relations) || anyNA(relations) || !all(relations %in% labels)) {
        cli::cli_abort("{.arg relations} must be GO.db relation labels ({.val {labels}}), not {.val {relations}}.",
            class = "go_propagation_invalid")
    }
    if (!is.data.frame(term2gene_direct) || !all(c("TERM", "GENE") %in% names(term2gene_direct))) {
        cli::cli_abort("{.arg term2gene_direct} must be a data frame with columns TERM and GENE.",
            class = "go_propagation_invalid")
    }
    direct <- unique(data.frame(TERM = as.character(term2gene_direct$TERM), GENE = as.character(term2gene_direct$GENE),
        stringsAsFactors = FALSE))
    direct <- direct[!is.na(direct$TERM) & !is.na(direct$GENE), , drop = FALSE]
    aspects <- goDbTermAspects()
    term_aspect <- unname(aspects[direct$TERM])
    not_in_db <- sort(unique(direct$TERM[is.na(term_aspect)]))
    not_in_aspect <- unique(direct$TERM[!is.na(term_aspect) & term_aspect != ontology])
    kept <- direct[!is.na(term_aspect) & term_aspect == ontology, , drop = FALSE]
    pairs <- goAncestorPairs(unique(kept$TERM), goParentTable(ontology), relations)
    joined <- merge(kept, pairs, by.x = "TERM", by.y = "term")
    propagated <- unique(data.frame(TERM = as.character(joined$node), GENE = as.character(joined$GENE),
        stringsAsFactors = FALSE))
    propagated <- propagated[order(propagated$TERM, propagated$GENE, method = "radix"), , drop = FALSE]
    rownames(propagated) <- NULL
    release <- goDbRelease()
    attr(propagated, "relations") <- relations
    attr(propagated, "go_db_version") <- release$version
    attr(propagated, "go_source_date") <- release$date
    attr(propagated, "not_in_go_db") <- length(not_in_db)
    attr(propagated, "not_in_go_db_examples") <- utils::head(not_in_db, 20L)
    attr(propagated, "not_in_aspect") <- length(not_in_aspect)
    attr(propagated, "pairs_direct") <- nrow(kept)
    attr(propagated, "pairs_propagated") <- nrow(propagated)
    propagated
}

#' Check a `go_propagation` value
#' @param go_propagation `"none"` or a name of `goPropagationRelations`.
#' @return `go_propagation`, invisibly; aborts with class `go_propagation_invalid` otherwise.
#' @noRd
checkGoPropagation <- function(go_propagation) {
    allowed <- c("none", names(goPropagationRelations))
    if (!is.character(go_propagation) || length(go_propagation) != 1L || !go_propagation %in% allowed) {
        cli::cli_abort("{.arg go_propagation} must be one of {.val {allowed}}, not {.val {go_propagation}}.",
            class = "go_propagation_invalid")
    }
    invisible(go_propagation)
}

#' One aspect's direct (Entry, TERM) pairs from the UniProt annotation table's `go_id_*` column
#' @noRd
uniprotGoAspectPairs <- function(go_annotations, column) {
    go_annotations |>
        dplyr::filter(!is.na(.data[[column]])) |>
        tidyr::separate_rows(dplyr::all_of(column), sep = "; ") |>
        dplyr::select(dplyr::all_of(c("Entry", column))) |>
        dplyr::rename(TERM = dplyr::all_of(column))
}

#' The proteomics clusterProfiler route's GO gene sets from the UniProt annotation table (KIT-173)
#'
#' With `"none"`, the tables `processEnrichments` built before KIT-173, exactly: the direct annotations of the
#' `go_id_go_*` columns, the three aspects pooled. Otherwise each aspect's direct annotations are propagated to their
#' parent terms over the value's relations (`propagateGoAnnotations()`), the three aspect roots are removed (in one
#' family pooled over the aspects a root's test reflects how well proteins are studied, not their biology), and
#' `term2gene` is sorted by term then protein, so the result does not depend on the annotation table's row order.
#' @param go_annotations The UniProt annotation table (`Entry`, `go_id_go_biological_process`,
#'   `go_id_go_molecular_function`, `go_id_go_cellular_compartment`; `Reviewed` when present).
#' @param go_propagation `"none"`, `"is_a_part_of"` or `"is_a_part_of_regulates"`.
#' @return `list(term2gene, term2name, all_terms)`: `term2gene` (`TERM`, `Entry`) and `term2name` (`TERM`, `NAME`)
#'   as `clusterProfiler::enricher` takes them, and `all_terms` (`Entry`, `TERM`, `ONTOLOGY`). On a propagating value
#'   the list carries the attributes `go_propagation`, `relations`, `go_db_version`, `go_source_date` (once),
#'   `not_in_go_db`, `not_in_go_db_examples`, `not_in_aspect`, `pairs_direct` (summed over the aspects),
#'   `pairs_propagated` (the pairs of `term2gene`, roots removed), `root_only_proteins` (the proteins whose direct GO
#'   ids are all aspect roots, which leave the background) and, when `Reviewed` is present, `unreviewed_pairs` (the
#'   pairs of `term2gene` whose protein is unreviewed).
#' @noRd
buildUniprotGoGeneSets <- function(go_annotations, go_propagation) {
    checkGoPropagation(go_propagation)
    if (identical(go_propagation, "none")) return(buildUniprotGoGeneSetsDirect(go_annotations))
    columns <- c(BP = "go_id_go_biological_process", MF = "go_id_go_molecular_function",
        CC = "go_id_go_cellular_compartment")
    relations <- goPropagationRelations[[go_propagation]]
    direct <- lapply(columns, function(column) uniprotGoAspectPairs(go_annotations, column))
    propagated <- lapply(names(columns), function(aspect) propagateGoAnnotations(
        data.frame(TERM = direct[[aspect]]$TERM, GENE = direct[[aspect]]$Entry, stringsAsFactors = FALSE),
        aspect, relations))
    names(propagated) <- names(columns)
    all_terms <- do.call(rbind, lapply(names(columns), function(aspect) {
        pairs <- propagated[[aspect]]
        pairs <- pairs[pairs$TERM != goAspectRoots[[aspect]], , drop = FALSE]
        data.frame(Entry = pairs$GENE, TERM = pairs$TERM, ONTOLOGY = rep(aspect, nrow(pairs)),
            stringsAsFactors = FALSE)
    }))
    term2gene <- unique(all_terms[, c("TERM", "Entry")])
    term2gene <- term2gene[order(term2gene$TERM, term2gene$Entry, method = "radix"), , drop = FALSE]
    rownames(term2gene) <- NULL
    terms <- unique(term2gene$TERM)
    term2name <- data.frame(TERM = terms, NAME = if (length(terms)) unname(AnnotationDbi::Term(terms)) else character(),
        stringsAsFactors = FALSE)
    direct_all <- do.call(rbind, direct)
    root_only <- tapply(direct_all$TERM %in% goAspectRoots, direct_all$Entry, all)
    sets <- list(term2gene = term2gene, term2name = term2name, all_terms = all_terms)
    release <- goDbRelease()
    count <- function(name) sum(vapply(propagated, function(p) attr(p, name), numeric(1)))
    attr(sets, "go_propagation") <- go_propagation
    attr(sets, "relations") <- relations
    attr(sets, "go_db_version") <- release$version
    attr(sets, "go_source_date") <- release$date
    attr(sets, "not_in_go_db") <- count("not_in_go_db")
    attr(sets, "not_in_go_db_examples") <- utils::head(sort(unique(unlist(lapply(propagated, attr,
        "not_in_go_db_examples")))), 20L)
    attr(sets, "not_in_aspect") <- count("not_in_aspect")
    attr(sets, "pairs_direct") <- count("pairs_direct")
    attr(sets, "pairs_propagated") <- nrow(term2gene)
    attr(sets, "root_only_proteins") <- sort(names(root_only)[root_only])
    reviewed_column <- intersect(c("Reviewed", "reviewed"), names(go_annotations))
    if (length(reviewed_column)) {
        status <- stats::setNames(as.character(go_annotations[[reviewed_column[[1L]]]]), go_annotations$Entry)
        attr(sets, "unreviewed_pairs") <- sum(status[term2gene$Entry] %in% "unreviewed")
    }
    sets
}

#' The direct gene sets `processEnrichments` built before KIT-173, unchanged
#' @noRd
buildUniprotGoGeneSetsDirect <- function(go_annotations) {
    bp_terms <- go_annotations |>
      dplyr::filter(!is.na(.data$go_id_go_biological_process)) |>
      tidyr::separate_rows(.data$go_id_go_biological_process, sep = "; ") |>
      dplyr::select(.data$Entry, .data$go_id_go_biological_process) |>
      dplyr::rename(TERM = .data$go_id_go_biological_process)

    mf_terms <- go_annotations |>
      dplyr::filter(!is.na(.data$go_id_go_molecular_function)) |>
      tidyr::separate_rows(.data$go_id_go_molecular_function, sep = "; ") |>
      dplyr::select(.data$Entry, .data$go_id_go_molecular_function) |>
      dplyr::rename(TERM = .data$go_id_go_molecular_function)

    cc_terms <- go_annotations |>
      dplyr::filter(!is.na(.data$go_id_go_cellular_compartment)) |>
      tidyr::separate_rows(.data$go_id_go_cellular_compartment, sep = "; ") |>
      dplyr::select(.data$Entry, .data$go_id_go_cellular_compartment) |>
      dplyr::rename(TERM = .data$go_id_go_cellular_compartment)

    all_terms <- rbind(
      cbind(bp_terms, ONTOLOGY = "BP"),
      cbind(mf_terms, ONTOLOGY = "MF"),
      cbind(cc_terms, ONTOLOGY = "CC")
    )

    term2gene <- all_terms |>
      dplyr::select(.data$TERM, .data$Entry) |>
      dplyr::distinct()

    term2name <- data.frame(
      TERM = unique(all_terms$TERM),
      NAME = purrr::map_chr(unique(all_terms$TERM),
                            ~tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[.x]]),
                                      error = function(e) .x))
    )

    list(term2gene = term2gene, term2name = term2name, all_terms = all_terms)
}

#' Raw GO ids the installed GO.db does not hold, by class (KIT-173)
#'
#' Called where raw GO ids first meet GO.db (`getUniprotAnnotations`, before `uniprotGoIdToTerm` drops them). An id is
#' obsolete when GO.db lists it as obsolete, secondary when it is a secondary id of a live term, unknown otherwise
#' (newer than the GO.db release, or not a GO id).
#' @param raw_ids Character vector of GO ids (repeats allowed; NA and empty strings ignored).
#' @return `list(obsolete, secondary, unknown, total, examples, go_db_version, go_source_date)`: counts of distinct
#'   ids, up to 20 example ids, and the GO.db release in force.
#' @noRd
goIdsNotInGoDb <- function(raw_ids) {
    requireOptionalPackage("GO.db", "GO identifier counts")
    ids <- unique(as.character(raw_ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    missing <- sort(setdiff(ids, names(goDbTermAspects())))
    con <- GO.db::GO_dbconn()
    obsolete <- DBI::dbGetQuery(con, "SELECT go_id FROM go_obsolete")$go_id
    secondary <- DBI::dbGetQuery(con, "SELECT secondary FROM go_synonym WHERE secondary IS NOT NULL")$secondary
    is_obsolete <- missing %in% obsolete
    is_secondary <- !is_obsolete & missing %in% secondary
    release <- goDbRelease()
    list(obsolete = sum(is_obsolete), secondary = sum(is_secondary), unknown = sum(!is_obsolete & !is_secondary),
        total = length(missing), examples = utils::head(missing, 20L), go_db_version = release$version,
        go_source_date = release$date)
}
