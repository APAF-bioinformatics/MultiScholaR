
# Function used for parsing a list of minimum or maximum gene set size from command line
#'@export
parseNumList <-  function ( input_text ) {
  if( str_detect( input_text, "[.,;:]")) {
    str_split( input_text, "[.,;:]")[[1]] %>% purrr::map_int( as.integer)
  } else {
    return( as.integer( input_text ))
  }
}

#'@export
convertIdToAnnotation <- function( id, id_to_annotation_dictionary) {

    return( ifelse( !is.null(id_to_annotation_dictionary[[id]] ),
                    id_to_annotation_dictionary[[id]] ,
                    NA_character_))

}


#'@param go_annot: Go annotation table.
#'@export
oneGoEnrichment <- function(go_annot, background_list, go_aspect, query_list, id_to_annotation_dictionary,
                              annotation_id, protein_id, aspect_column, p_val_thresh, min_gene_set_size,  max_gene_set_size  ) {

  join_condition <- rlang::set_names( c( colnames(background_list)[1]),
                                      c( quo_name(enquo( protein_id)) ) )

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

  no_singleton_terms_query_gene_list <- intersect( query_list ,
                                                   term_to_gene_tbl_filt_no_singleton %>%
                                                     dplyr::distinct(gene) %>%
                                                     dplyr::pull(gene))

  # print(quo_name(enquo(aspect_column)))
  # print(go_aspect)
  # print(nrow( go_annot))
  # print(nrow( go_annot_filt))
  # print(nrow( filtered_go_terms))
  # print(nrow(term_to_gene_tbl_filt))
  # print( length( no_singleton_terms_query_gene_list) )
  # print( nrow( term_to_gene_tbl_filt_no_singleton))
  # print(p_val_thresh )
  # print( min_gene_set_size)
  # print( max_gene_set_size)

  enrichment_result <- enricher(
    no_singleton_terms_query_gene_list,
    pvalueCutoff = p_val_thresh,
    pAdjustMethod = "BH",
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    qvalueCutoff = p_val_thresh,
    TERM2GENE =term_to_gene_tbl_filt_no_singleton
  )

  if(!is.null(enrichment_result) ) {

    output_table <-  as.data.frame( enrichment_result ) %>%
      dplyr::mutate( term = purrr::map_chr( ID,
                                            function(id) {
                                              convertIdToAnnotation(id,
                                                                    id_to_annotation_dictionary) } )) %>%
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

  } else {

    return( NULL)
  }

}



#'@export
runOneGoEnrichmentInOutFunction <- function(comparison_column,
                                            protein_id_column,
                                            go_annot = go_annot,
                                            background_list,
                                            id_to_annotation_dictionary,
                                            annotation_id,
                                            protein_id,
                                            aspect_column,
                                            p_val_thresh,
                                            ## These are the parameters that are usually presented as different combinations via the cross functions
                                            names_of_genes_list,
                                            input_table,
                                            go_aspect,
                                            input_comparison,
                                            min_gene_set_size,
                                            max_gene_set_size) {

  oneGoEnrichmentPartial <- NA
  if(!is.null(aspect_column )) {
    oneGoEnrichmentPartial <- purrr::partial( oneGoEnrichment,
                                              go_annot = go_annot,
                                              background_list = background_list,
                                              id_to_annotation_dictionary=id_to_annotation_dictionary,
                                              annotation_id= {{annotation_id}},
                                              protein_id={{protein_id}},
                                              aspect_column=!!rlang::sym(aspect_column),
                                              p_val_thresh=p_val_thresh)
  } else {
    oneGoEnrichmentPartial <- purrr::partial( oneGoEnrichment,
                                              go_annot = go_annot,
                                              background_list = background_list,
                                              id_to_annotation_dictionary=id_to_annotation_dictionary,
                                              annotation_id= {{annotation_id}},
                                              protein_id={{protein_id}},
                                              aspect_column=aspect_column,
                                              p_val_thresh=p_val_thresh)
  }

  query_list <- input_table %>%
    dplyr::filter( {{comparison_column}} == input_comparison) %>%
    pull( {{protein_id_column}})

  # print( paste( "size of query list = ", length( query_list)) )

  enrichment_temp <- oneGoEnrichmentPartial(
    go_aspect = go_aspect,
    query_list=query_list,
    min_gene_set_size=min_gene_set_size,
    max_gene_set_size=max_gene_set_size )

  if( !is.null( enrichment_temp)) {
    enrichment_result <- enrichment_temp %>%
      dplyr::mutate(  {{comparison_column}} := input_comparison ) %>%
      dplyr::mutate( names_of_genes_list = names_of_genes_list)

    return(enrichment_result )

  } else {
    return (NULL)
  }

}

#'@export
convertProteinAccToGeneSymbol <- function( gene_id_list, dictionary ) {

  purrr::map_chr( gene_id_list,
                  ~{ ifelse( . %in% names(dictionary ),
                             dictionary[[.]],
                             NA_character_)   } )  %>%
    paste( collapse="/")
}


#'@export
buildAnnotationIdToAnnotationNameDictionary <- function(input_table, annotation_column, annotation_id_column) {

  id_to_annotation_dictionary <- NA

  dictionary_pair <- input_table %>%
    distinct({{annotation_column}},
             {{annotation_id_column}})

  id_to_annotation_dictionary <- dictionary_pair %>%
    pull({{annotation_column}} )

  names(id_to_annotation_dictionary ) <-  dictionary_pair %>%
    pull( {{annotation_id_column}})

  id_to_annotation_dictionary

}


#'@export
buildOneProteinToAnnotationList <- function( input_table, annotation_id, protein_id ) {
  temp_table <- input_table %>%
    group_by( {{annotation_id}}) %>%
    nest( ) %>%
    ungroup()  %>%
    mutate( gene_set = purrr::map( data, ~{ (.) %>% pull( {{protein_id}} )} ))

  gene_set_list <- temp_table %>% pull(gene_set)

  names(gene_set_list ) <-temp_table %>% pull({{annotation_id}})

  gene_set_list
}

#'@export
listifyTableByColumn  <- function(input_table, column_name) {

  nested_table <- input_table %>%
    dplyr::filter(!is.na( {{column_name}})) %>%
    group_by({{column_name}}) %>%
    nest( ) %>%
    ungroup()

  list_of_tables <- nested_table %>%
    pull(data)

  names( list_of_tables) <- nested_table %>%
    distinct({{column_name}})  %>%
    pull( {{column_name}})

  list_of_tables
}
