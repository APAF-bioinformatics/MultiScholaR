# ----------------------------------------------------------------------------
# goIdToTerm
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Convert a list of Gene Ontology IDs to their respective human readable name (e.g. GO term).
#' @param go_string A string consisting of a list of Gene Ontology ID, separated by a delimiter
#' @param goterms Output from running \code{goterms <- Term(GOTERM)} from the GO.db library.
#' @param gotypes Output from running \code{gotypes <- Ontology(GOTERM)} from the GO.db library.
#' @return A table with three columns. go_biological_process, go_celluar_compartment, and go_molecular_function. Each column is a list of gene ontology terms, separated by '; '.
#' @export
#' @examples
#' go_string <- "GO:0016021; GO:0030659; GO:0031410; GO:0035915; GO:0042742; GO:0045087; GO:0045335; GO:0050829; GO:0050830"
#' go_id_to_term(go_string)
goIdToTerm <- function(go_string, sep = "; ", goterms, gotypes) {

  if (!is.na(go_string)) {
    go_string_tbl <- data.frame(go_id = go_string) |>
      separate_rows(go_id, sep = sep) |>
      mutate(go_term = purrr::map_chr(go_id, function(x) { if (x %in% names(goterms)) { return(goterms[[x]]) }; return(NA) })) |>
      mutate(go_type = purrr::map_chr(go_id, function(x) { if (x %in% names(gotypes)) { return(gotypes[[x]]) }; return(NA) })) |>
      group_by(go_type) |>
      summarise(go_term = paste(go_term, collapse = "; ")) |>
      ungroup() |>
      mutate(go_type = case_when(go_type == "BP" ~ "go_biological_process",
                                 go_type == "CC" ~ "go_cellular_compartment",
                                 go_type == "MF" ~ "go_molecular_function")) |>
      pivot_wider(names_from = go_type,
                  values_from = go_term)

    return(go_string_tbl)
  }
  return(NA)
}

# ----------------------------------------------------------------------------
# uniprotGoIdToTerm
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Convert UniProt Accession to GO Term
#' @param uniprot_dat  a table with uniprot accessions and a column with GO-ID
#' @param uniprot_id_column The name of the column with the uniprot accession, as a tidyverse header format, not a string
#' @param go_id_column The name of the column with the GO-ID, as a tidyverse header format, not a string
#' @param goterms Output from running \code{goterms <- Term(GOTERM)} from the GO.db library.
#' @param gotypes Output from running \code{gotypes <- Ontology(GOTERM)} from the GO.db library.
#' @return A table with three columns. go_biological_process, go_celluar_compartment, and go_molecular_function. Each column is a list of gene ontology terms, separated by '; '.
#' @export
uniprotGoIdToTerm <- function(uniprot_dat, uniprot_id_column = UNIPROTKB
                              , go_id_column = `GO-IDs`,  sep = "; "
                              , goterms = AnnotationDbi::Term(GO.db::GOTERM)
                              , gotypes = AnnotationDbi::Ontology(GO.db::GOTERM)) {

  uniprot_acc_to_go_id <- uniprot_dat |>
    dplyr::distinct({{uniprot_id_column}}, {{go_id_column}}) |>
    separate_rows({{go_id_column}}, sep = sep) |>
    dplyr::distinct({{uniprot_id_column}}, {{go_id_column}}) |>
    dplyr::filter(!is.na({{go_id_column}}))

  go_term_temp <- uniprot_acc_to_go_id |>
    dplyr::distinct({{go_id_column}}) |>
    mutate(go_term = purrr::map_chr({{go_id_column}}, function(x) { if (x %in% names(goterms)) { return(goterms[[x]]) }; return(NA) })) |>
    mutate(go_type = purrr::map_chr({{go_id_column}}, function(x) { if (x %in% names(gotypes)) { return(gotypes[[x]]) }; return(NA) })) |>
    mutate(go_type = case_when(go_type == "BP" ~ "go_biological_process",
                               go_type == "CC" ~ "go_cellular_compartment",
                               go_type == "MF" ~ "go_molecular_function"))

  uniprot_acc_to_go_term <- uniprot_acc_to_go_id |>
    left_join(go_term_temp, by = join_by({{go_id_column}} == {{go_id_column}})) |>
    dplyr::filter(!is.na(go_term)) |>
    group_by({{uniprot_id_column}}, go_type) |>
    summarise(go_id = paste({{go_id_column}}, collapse = "; ")
              , go_term = paste(go_term, collapse = "; ")
              , .groups = 'drop' ) |>
    pivot_wider(id_cols = {{uniprot_id_column}},
                names_from = go_type,
                values_from = c(go_term, go_id))


  output_uniprot_dat <- uniprot_dat |>
    left_join(uniprot_acc_to_go_term, by = join_by({{uniprot_id_column}} == {{uniprot_id_column}})) |>
    dplyr::select( -{{go_id_column}})

  return(output_uniprot_dat)

}

