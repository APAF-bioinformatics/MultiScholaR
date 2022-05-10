
# Function used for parsing a list of minimum or maximum gene set size from command line
parseNumList <-  function ( input_text ) {
  if( str_detect( input_text, "[.,;:]")) {
    str_split( input_text, "[.,;:]")[[1]] %>% purrr::map_int( as.integer)
  } else {
    return( as.integer( input_text ))
  }
}

convertIdToAnnotation <- function( id, id_to_annotation_dictionary, go_aspect ) {

  # if( !is.na( go_aspect)) {
  #   return( ifelse( !is.null(id_to_annotation_dictionary[[id]] ),
  #                   Term( id_to_annotation_dictionary[[id]] ),
  #                   NA_character_))
  # } else {
    return( ifelse( !is.null(id_to_annotation_dictionary[[id]] ),
                    id_to_annotation_dictionary[[id]] ,
                    NA_character_))
  # }

}


one_go_enrichment <- function(go_annot, background_list, go_aspect, query_list, id_to_annotation_dictionary,
                              annotation_id, protein_id, aspect_column, p_val_thresh, min_gene_set_size,  max_gene_set_size  ) {

  join_condition <- rlang::set_names( c( colnames(background_list)[1]),
                                      c( quo_name(enquo( protein_id)) ) )

  #print("reached here 1")
  if ( !is.na( go_aspect)) {
    go_annot_filt <- go_annot %>%
      dplyr::filter( {{aspect_column}} == go_aspect)
  } else {
    go_annot_filt <- go_annot
  }

  filtered_go_terms <- go_annot_filt %>%
    inner_join( background_list, by =join_condition )  %>%
    group_by( {{annotation_id}} ) %>%
    summarise(counts =n()) %>%
    ungroup() %>%
    arrange(desc(counts)) %>%
    dplyr::filter( counts <= max_gene_set_size & counts >= min_gene_set_size ) %>%
    dplyr::select(-counts)

  term_to_gene_tbl_filt <- go_annot_filt %>%
    inner_join( background_list, by =join_condition )  %>%
    dplyr::inner_join( filtered_go_terms, by = quo_name(enquo( annotation_id))  )  %>%
    dplyr::rename( gene = quo_name(enquo(protein_id )) ,
                   term = quo_name(enquo( annotation_id)) ) %>%
    dplyr::select(term, gene) %>%
    dplyr::distinct( term, gene )

  # print(nrow(term_to_gene_tbl_filt))

  ## Avoid singleton GO terms
  terms_to_avoid <- term_to_gene_tbl_filt %>%
    distinct() %>%
    dplyr::inner_join(data.frame(uniprot_acc = query_list), by=c("gene" = "uniprot_acc") )  %>%
    distinct() %>%
    group_by( term) %>%
    summarise(counts =n()) %>%
    ungroup() %>%
    dplyr::filter( counts < 2)

  term_to_gene_tbl_filt_no_singleton <- term_to_gene_tbl_filt %>%
    dplyr::anti_join( terms_to_avoid, by="term")

  enrichment_result <- enricher(
    intersect( query_list ,
               term_to_gene_tbl_filt_no_singleton %>% distinct(gene) %>% pull(gene)),
    pvalueCutoff = p_val_thresh,
    pAdjustMethod = "BH",
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    qvalueCutoff = p_val_thresh,
    TERM2GENE =term_to_gene_tbl_filt_no_singleton
  )


  output_table <-  as.data.frame( enrichment_result ) %>%
    dplyr::mutate( term = purrr::map_chr( ID,
                                          function(id) {
                                            convertIdToAnnotation(id,
                                                                  id_to_annotation_dictionary,
                                                                  go_aspect) } )) %>%
    dplyr::relocate( term, .before="Description") %>%
    dplyr::mutate( min_gene_set_size = min_gene_set_size,
                   max_gene_set_size = max_gene_set_size )


  output_table_with_go_aspect <- NA
  if ( !is.na( go_aspect)) {
    output_table_with_go_aspect <- output_table %>%
      dplyr::mutate( {{aspect_column}} := go_aspect)
  } else {
    output_table_with_go_aspect <- output_table
  }

  return(output_table_with_go_aspect)

}




runOneGoEnrichmentInOutFunction <- function(input_table,
                                            comparison_column,
                                            protein_id_column,
                                            ## These are the parameters that are usually presented as different combinations via the cross functions
                                            go_aspect,
                                            input_comparison,
                                            min_gene_set_size,
                                            max_gene_set_size) {

  # print("start enrichment")
  query_list <- input_table %>%
    dplyr::filter( {{comparison_column}} == input_comparison) %>%
    pull( {{protein_id_column}})

  enrichment_result <- one_go_enrichment_partial(
    go_aspect = go_aspect,
    query_list=query_list,
    min_gene_set_size=min_gene_set_size,
    max_gene_set_size=max_gene_set_size ) %>%
    dplyr::mutate(  {{comparison_column}} := input_comparison )

  return(enrichment_result )

}


convertProteinAccToGeneSymbol <- function( gene_id_list, dictionary ) {

  purrr::map_chr( gene_id_list,
                  ~{ ifelse( . %in% names(dictionary ),
                             dictionary[[.]],
                             NA_character_)   } )  %>%
    paste( collapse="/")
}
