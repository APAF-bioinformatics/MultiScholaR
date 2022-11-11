
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
                                      c( as_name(enquo( protein_id)) ) )

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
    dplyr::inner_join( filtered_go_terms, by = as_name(enquo( annotation_id))  )  %>%
    dplyr::rename( gene = as_name(enquo(protein_id )) ,
                   term = as_name(enquo( annotation_id)) ) %>%
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

  # print(as_name(enquo(aspect_column)))
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


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
cmriCamera <- function( abundance_mat, design_mat, contrast_name, index_name, index,
                        contrast, min_set_size = 4, max_set_size=300) {

  print(paste("contrast_name =", contrast_name))
  print(paste("index_name =", index_name))

  print(contrast)

  this_contrast <- contrast
  msigdb_ids <- index

  #convert gene sets into a list of gene indices
  camera_indices <- ids2indices(msigdb_ids,
                                rownames(abundance_mat))

  ## At least two genes in the gene set
  set_sizes <- purrr::map(camera_indices, length)
  camera_indices_filt <- camera_indices[ set_sizes >= min_set_size &
                                           set_sizes <= max_set_size]



  camera_result <- NA
  if (length(camera_indices_filt) > 0) {
    camera_result <- camera(y = abundance_mat,
                            design = design_mat,
                            index = camera_indices_filt,
                            contrast = this_contrast) %>%
      mutate( rank = row_number())
  }


  info_list <- list(camera = camera_result,
                    y = abundance_mat,
                    design = design_mat,

                    index_name = index_name,
                    index = camera_indices_filt,

                    contrast_name = contrast_name,
                    contrast = this_contrast,
                    min_set_size = min_set_size,
                    max_set_size= max_set_size)

  return(info_list)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
runGsea <- function(index_name, contrast_name, list_of_de_proteins, list_of_gene_sets, min_set_size = 4, max_set_size = 300) {

  gene_list <- list_of_de_proteins[[contrast_name]]

  msigdb_gene_set <- geneIds(list_of_gene_sets[[index_name]])

  query_gene_list <- data.frame(gene = names(gene_list))

  term_to_gene_tab <- tibble(term = names(msigdb_gene_set), gene = msigdb_gene_set) %>%
    unnest(gene) %>%
    dplyr::inner_join(query_gene_list, by = c("gene"))

  terms_to_keep <- term_to_gene_tab %>%
    group_by(term) %>%
    summarise(counts = n()) %>%
    ungroup() %>%
    dplyr::filter( counts >= min_set_size &
                     counts <= max_set_size) %>%
    dplyr::select(-counts)

  term_to_gene_tab_filt <- term_to_gene_tab %>%
    inner_join(terms_to_keep, by = "term") %>%
    mutate(gene = as.character(gene))

  ## Check that there is overlap
  # intersect( names( gene_list_final) ,  unique( term_to_gene_tab_filt$gene )) %>% length


  gsea_results <- GSEA(geneList = gene_list, TERM2GENE = as.data.frame(term_to_gene_tab_filt))

  return(gsea_results)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
runEnricher <- function(index_name, contrast_name, list_of_de_proteins, list_of_gene_sets, min_set_size = 4, max_set_size = 300) {

  gene_list <- list_of_de_proteins[[contrast_name]]

  msigdb_gene_set <- geneIds(list_of_gene_sets[[index_name]])

  query_gene_list <- data.frame(gene = gene_list)

  term_to_gene_tab <- tibble(term = names(msigdb_gene_set), gene = msigdb_gene_set) %>%
    unnest(gene) %>%
    dplyr::inner_join(query_gene_list, by = c("gene"))

  terms_to_keep <- term_to_gene_tab %>%
    group_by(term) %>%
    summarise(counts = n()) %>%
    ungroup() %>%
    dplyr::filter( counts >= min_set_size &
                     counts <= max_set_size  ) %>%
    dplyr::select(-counts)

  term_to_gene_tab_filt <- term_to_gene_tab %>%
    inner_join(terms_to_keep, by = "term") %>%
    mutate(gene = as.character(gene))

  ## Check that there is overlap
  # intersect( names( gene_list_final) ,  unique( term_to_gene_tab_filt$gene )) %>% length

  print(intersect(gene_list, unique(term_to_gene_tab_filt$gene)) %>% length)

  gsea_results <- enricher(gene = gene_list, TERM2GENE = as.data.frame(term_to_gene_tab_filt))

  return(gsea_results)

}

#'@export
getUniprotAccToGeneSymbolDictionary <- function( input_table,
              protein_id_lookup_column,
          gene_symbol_column,
          protein_id) {

  # Clean up protein ID to gene sybmol table
  uniprot_to_gene_symbol <- input_table  %>%
    dplyr::select( {{protein_id_lookup_column}},
                   {{gene_symbol_column}}) %>%
    dplyr::rename( {{protein_id}} := as_name(enquo(protein_id_lookup_column)) ) %>%
    dplyr::rename( gene_symbol = as_name(enquo( gene_symbol_column)) ) %>%
    dplyr::mutate( gene_symbol = str_split(  gene_symbol , " " ) %>%
                     purrr::map_chr( 1)) %>%
    dplyr::distinct( {{protein_id}}, gene_symbol)

  ## Convert to lookup dictionary
  uniprot_to_gene_symbol_dict <- uniprot_to_gene_symbol %>%
    pull( gene_symbol)

  names( uniprot_to_gene_symbol_dict )  <- uniprot_to_gene_symbol %>%
    pull( {{protein_id}} )

  uniprot_to_gene_symbol_dict

}

#######################################################
### Query revigo
#'@export
queryRevigo <- function( input_list,
                         cutoff=0.5,
                         speciesTaxon = 10090,
                         temp_file=NA) {

  userData <-  paste(input_list,  collapse= "\n")

  httr::POST(
    url = "http://revigo.irb.hr/Revigo.aspx",
    body = list(
      cutoff = as.character(cutoff),
      valueType = "pvalue",
      speciesTaxon = as.character(speciesTaxon),
      measure = "SIMREL",
      goList = userData
    ),
    # application/x-www-form-urlencoded
    encode = "form"
  )  -> res

  # print(res)

  dat <- httr::content(res, encoding = "UTF-8")


  dat <- stri_replace_all_fixed(dat, "\r", "")

  if(is.na( temp_file) |
     is.null(temp_file)) {
    temp_file <- tempfile(pattern = "temp_revigo",
                          tmpdir = tempdir(),
                          fileext = "html")
  }

  cat(dat, file=temp_file , fill = FALSE)

  html_doc <- rvest::read_html(dat, as.data.frame=T, stringsAsFactors = FALSE)

  revigo_tbl <- html_doc  %>%
    html_nodes("table") %>%
    purrr::map( ~html_table(.)) %>%
    discard( ~{ nrow(.) ==0 }) %>%
    bind_rows()

  if( file.exists( temp_file) ) {
    file.remove(temp_file)
  }

  revigo_tbl
}



#'@export
clusterPathways <- function ( input_table, added_columns, remove_duplicted_entries = TRUE ) {

    duplicated_entries <- input_table %>%
      mutate(set_type = case_when( str_detect( gene_set, "positive") ~"positive",
                                   str_detect( gene_set, "negative") ~ "negative",
                                   TRUE ~ "neutral")) %>%
      group_by( across(c( any_of(added_columns), comparison, set_type, annotation_id) ) ) %>%
      dplyr::summarise( temp_qvalue = min(qvalue )) %>%
      ungroup() %>%
      dplyr::group_by( across(c( any_of(added_columns), comparison, annotation_id) ) ) %>%
      dplyr::summarise(counts = n(),
                best_p_adj_value = min(temp_qvalue)) %>%
      ungroup() %>%
      dplyr::filter( counts > 1)

  if( remove_duplicted_entries == TRUE |
      remove_duplicted_entries == "delete" ) {
    input_table <- input_table  %>%
      anti_join( duplicated_entries, by =c("comparison" = "comparison",
                                           "annotation_id" = "annotation_id",
                                           added_columns))
  } else if( remove_duplicted_entries == "merge" ) {

    duplicates_tbl <- input_table %>%
      inner_join( duplicated_entries, by =c("comparison" = "comparison",
                                            "annotation_id" = "annotation_id",
                                            added_columns)) %>%
      dplyr::filter( qvalue == best_p_adj_value ) %>%
      mutate( gene_set = "shared" )

    input_table <- input_table  %>%
      anti_join( duplicated_entries, by =c("comparison" = "comparison",
                                           "annotation_id" = "annotation_id",
                                            added_columns)) %>%
      bind_rows( duplicates_tbl )

  }

    scores_for_clustering <- input_table %>%
      mutate( neg_log_p_value = -log10( p.adjust) ) %>%
      mutate(score = case_when( str_detect( gene_set, "positive") ~neg_log_p_value,
                                str_detect( gene_set, "negative") ~ -1* neg_log_p_value,
                                TRUE ~ neg_log_p_value))  %>%
      pivot_wider( id_cols = c(annotation_id),
                   names_from = c(any_of(added_columns), comparison, gene_set, go_type) ,
                   values_from = score,
                   values_fill = 0 )    %>%
      column_to_rownames("annotation_id") %>%
      as.matrix()

  pathways_clustered <- hclust(dist(scores_for_clustering))

  pathways_sorting <- cutree(pathways_clustered, k=1:ncol(scores_for_clustering)) %>%
    as.data.frame %>%
    rownames_to_column("Term") %>%
    arrange( across( matches("\\d+"))) %>%
    mutate( ordering = row_number()) %>%
  arrange(ordering)

   annot_heat_map_ordering <-  input_table %>%
     mutate( neg_log_p_value = -log10( p.adjust) )  %>%
     dplyr::select(  c(any_of(added_columns), comparison, annotation_id, term,  neg_log_p_value,  gene_set, go_type )) %>%
     mutate( annotation_id = as.character(annotation_id)) %>%
     left_join(pathways_sorting, by=c("annotation_id" = "Term")) %>%
     arrange(ordering)


   annot_heat_map_ordered <- annot_heat_map_ordering %>%
     mutate(term = factor( term,  levels = unique(annot_heat_map_ordering$term)))

  annot_heat_map_ordered
}

########################

#'@export
getEnrichmentHeatmap <- function( input_table, x_axis, input_go_type, input_plot_title,
                                  facet_by_column = NA, xaxis_levels=NA,
                                  scales="fixed") {

  get_shape <- list( negative_list = 25,
                     positive_list=24,
                     positive_only = 24,
                     negative_only = 25,
                     positive_sum_sig_phosphosites = 24,
                     negative_sum_sig_phosphosites = 25,
                     shared = 1,
                     positive_plus_overlap = 24,
                     negative_plus_overlap = 25,
                     all_significant = 1,
                     overlap_only = 1)

  my_get_shape <-function(x) {
    if(x %in% names(get_shape)) {
      return( get_shape[[x]])
    } else  {
      return( 16)
    }
  }

  get_colour <- list( negative_list = "blue", positive_list = "red",
                      positive_only = "red",
                      negative_only = "blue",
                      positive_sum_sig_phosphosites = "red",
                      negative_sum_sig_phosphosites = "blue",
                      shared = "black",
                      positive_plus_overlap = "red",
                      negative_plus_overlap = "blue",
                      all_significant = "black",
                      overlap_only = "black")

  my_get_colour <-function(x) {
    if(x %in% names(get_colour)) {
      return( get_colour[[x]])
    } else  {
      return( "black")
    }
  }

  table_filtering <- NA
  if(!is.na( input_go_type)) {
    table_filtering <- input_table %>%
      dplyr::filter(  go_type == input_go_type)
  } else {
    table_filtering <- input_table
  }

  table_shape_colour <- table_filtering %>%
    mutate( use_shape = purrr::map_dbl( gene_set, my_get_shape)) %>%
    mutate( use_colour = purrr::map_chr( gene_set, my_get_colour )) %>%
    mutate(term = factor( term,  levels = unique(input_table$term)))

  if( length(xaxis_levels) > 1  ) {

    # If we are manually ordering the x axis labels from left to right,
    # We need to make sure the factor levels in the input covers all the things we need to label.
    all_x_axis_labels <- table_shape_colour %>%
      distinct( {{x_axis}} ) %>%
      pull({{x_axis}})


    if( length(setdiff( all_x_axis_labels, xaxis_levels)) ==0) {
      table_shape_colour <- table_shape_colour %>%
        mutate( {{x_axis}} := factor( {{x_axis}}, levels=xaxis_levels))
    } else {
      print(setdiff( all_x_axis_labels, xaxis_levels))
      print( "Cannot locate x_axis ordering.")
      stop()
    }

  } else {
    table_shape_colour <- table_shape_colour %>%
      mutate( {{x_axis}} := purrr::map_chr( {{x_axis}}, as.character))
  }

  output_heat_map <- table_shape_colour %>%
    ggplot( aes(  {{x_axis}}, term,
                  fill = use_colour,
                  col = use_colour,
                  shape=use_shape,
                  size = neg_log_p_value)) +
    geom_point() +
    scale_size_continuous( name = "-log10(p-value)"  ) +
    scale_shape_identity() +
    scale_color_identity() +
    scale_fill_identity()  +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y  = element_text(face = "bold")) +
    theme(strip.text.y = element_text(angle = 0))  +
    scale_x_discrete(labels = function(input){ str_wrap(input, width = 15) }) +
    labs(title=input_plot_title)

  if( length(which( c(24, 25)  %in% c(table_shape_colour %>% pull( use_shape)) ) > 0 ) ) {
    output_heat_map <- output_heat_map +
      guides(size = guide_legend(override.aes = list(shape=17)),
             shape = guide_legend(override.aes = list( size = 5   )))
  }

  if( !is.na( quo_get_expr(enquo(facet_by_column) ) ) ) {
    if( as_name(enquo(facet_by_column)) %in% colnames(table_filtering )) {
      print("Using faceting")
      print( as_name( enquo( facet_by_column)) )

      output_heat_map <- output_heat_map  +
        facet_wrap( vars({{facet_by_column}}), scales=scales  )
    }
  }


  output_heat_map

}


#' @export
readEnrichmentResultFiles <- function( table_of_files, file_names_column=file_name, go_type="KEGG") {

  list_of_files <- table_of_files %>%
    pull( {{file_names_column}})

  added_columns <- setdiff(colnames(table_of_files),
          as_name(enquo(file_names_column)))

  ## Gets error if input table have zero rows, so need filtering to remove table with zero rows
  list_of_tables <- purrr::map( list_of_files, vroom::vroom)

  num_lines <- purrr::map_int( list_of_tables, nrow)

  list_of_tables_with_rows <- purrr::keep( list_of_tables, ~{nrow(.) > 0})

  names(list_of_tables_with_rows) <- list_of_files[num_lines > 0]

  cleaned_tbl <-  list_of_tables_with_rows %>%
    bind_rows( .id=as_name(enquo(file_names_column)))

  enriched_results_tbl <- cleaned_tbl %>%
    dplyr::rename(annotation_id = "ID", gene_set = "names_of_genes_list",
           min_set_size = "min_gene_set_size",
           max_set_size = "max_gene_set_size") %>%
    left_join( table_of_files, by = as_name(enquo(file_names_column)) ) %>%
    relocate( any_of(added_columns ) ,
              .before=as_name(enquo(file_names_column))) %>%
    dplyr::select(-{{file_names_column}})

  if ( ! "go_type" %in% colnames(enriched_results_tbl) ) {

    enriched_results_tbl <- enriched_results_tbl %>%
      dplyr::mutate( go_type = go_type)
  }

  return( enriched_results_tbl )

}

# enriched_results_tbl <- readEnrichmentResultFiles( table_of_files, go_type="KEGG")

#'@export
filterResultsWithRevigo <- function(enriched_results_tbl,  added_columns, is_run_revigo=TRUE, revigo_cutoff=0.7) {

  enrich_revigo <- NA

  if ( is_run_revigo == TRUE) {

    annotation_list <-  enriched_results_tbl %>%
      group_by( across( c(any_of(added_columns), comparison, gene_set, go_type) )) %>%
      nest() %>%
      ungroup() %>%
      mutate( annot_id_list = purrr::map( data, ~{ pull(., annotation_id)} ))

    annotation_list_revigo <-  annotation_list %>%
      mutate( revigo_results = purrr::map( annot_id_list,
                                           function(x){queryRevigo(x,
                                                                   cutoff=revigo_cutoff,
                                                                   speciesTaxon = species_taxon,
                                                                   temp_file=NA )}))

    # annotation_list_revigo %>%
    #   unnest(revigo_results) %>%
    #   colnames %>% print
    #
    # annotation_list_revigo %>%
    #   unnest(revigo_results) %>% head %>% print

    revigo_tbl <- annotation_list_revigo %>%
      unnest(revigo_results)  %>%
      dplyr::select(-data, - annot_id_list) %>%
      dplyr::rename(annotation_id = "Term ID")

    join_condition <- rlang::set_names(c("annotation_id", "comparison", "go_type", "gene_set", added_columns),
                                       c("annotation_id", "comparison", "go_type", "gene_set", added_columns))



    enrich_revigo <- enriched_results_tbl %>%
      dplyr::mutate( annotation_id = as.character( annotation_id)) %>%
      left_join( revigo_tbl %>%
                   dplyr::select(-Name),
                 by = join_condition) %>%
      dplyr::filter( Eliminated == "False" |
                       is.na(Eliminated))

  } else {

    enrich_revigo <- enriched_results_tbl
  }

  return( enrich_revigo)
}

#'@export
saveFilteredFunctionalEnrichmentTable <- function( enriched_results_tbl,
                                                   set_size_min,
                                                   set_size_max,
                                                   results_dir,
                                                   file_name  ) {

  max_excel_cell_length <- 32760

  vroom::vroom_write( enriched_results_tbl %>%
                        dplyr::filter( min_set_size == set_size_min,
                                       max_set_size == set_size_max),
                      file.path(results_dir,
                                paste0( file_name, ".tab" )))

  writexl::write_xlsx( enriched_results_tbl %>%
                         dplyr::filter( min_set_size == set_size_min,
                                        max_set_size == set_size_max) %>%
                         mutate_at( c( "gene_symbol"), ~substr(., 1, max_excel_cell_length)) ,
                       path=file.path(results_dir,
                                      paste0( file_name, ".xlsx" ) ))

  vroom::vroom_write( enriched_results_tbl ,
                      file.path(results_dir,
                                paste0( file_name, "_unfiltered.tab" )))

  writexl::write_xlsx( enriched_results_tbl %>%
                         mutate_at( c( "gene_symbol"), ~substr(., 1, max_excel_cell_length)) ,
                       path=file.path(results_dir,
                                      paste0( file_name, "_unfiltered.xlsx" ) ))

}


#'@export
evaluateBestMinMaxGeneSetSize <- function(enrichment_results_tble, added_columns) {

  plotting_data <- enrichment_results_tble %>%
    group_by(  across( c(any_of(added_columns), comparison, min_set_size, max_set_size, gene_set, go_type) ) ) %>%
    summarise( counts =n()) %>%
    ungroup() %>%
    mutate( set_size = paste(min_set_size, max_set_size, sep="-" ) ) %>%
    dplyr::mutate( gene_set_mod = ifelse(!is.na(go_type),
                                         paste(  gene_set, go_type, sep="-"),
                                         gene_set) )

  plotting_data %>%
    unite(  custom_comparison , comparison, any_of( added_columns ) ) %>%
  ggplot( aes( set_size, counts, group=custom_comparison)) +
    geom_line(aes(col=custom_comparison)) +
    theme (axis.text.x = element_text (angle = 90, vjust = 1))  +
    facet_wrap( . ~ gene_set_mod    , scales="free_y")

}


#'@export
drawListOfFunctionalEnrichmentHeatmaps <- function(enriched_results_tbl,
                                                   added_columns,
                                                   set_size_min,
                                                   set_size_max,
                                                   x_axis = Analysis_Type,
                                                   analysis_column = Analysis_Type,
                                                   facet_by_column = NA,
                                                   remove_duplicted_entries = TRUE,
                                                   xaxis_levels=NA,
                                                   scales="fixed") {

  input_table <- enriched_results_tbl %>%
    dplyr::filter( min_set_size == set_size_min,
                   max_set_size == set_size_max) %>%
    group_by(  across( c(any_of(added_columns), comparison, gene_set, go_type) )) %>%
    arrange( comparison, pvalue) %>%
    mutate(  ranking = row_number() ) %>%
    ungroup()

  if ( nrow(input_table) == 0 ) {
    stop("drawListOfFunctionalEnrichmentHeatmaps: No more rows for clustering analysis after gene set size filtering.")
  }

  list_of_columns_to_exclude <- c(as_name(enquo(x_axis)))
  if(! is.na( quo_get_expr(enquo(facet_by_column) ) ) ) {
    list_of_columns_to_exclude <- c(as_name(enquo(x_axis)), as_name(enquo(facet_by_column)))
  }

  annot_heat_map_ordered <- clusterPathways( input_table,
                                             added_columns,
                                             remove_duplicted_entries = remove_duplicted_entries) %>%
    unite(  {{analysis_column}} , comparison, any_of( c(setdiff(added_columns, list_of_columns_to_exclude))) )

  combinations <- annot_heat_map_ordered %>%
    distinct(  go_type)

  list_of_heatmaps <- purrr::pmap( combinations, function( go_type){
    print( paste(  go_type) )
    getEnrichmentHeatmap( input_table=annot_heat_map_ordered,
                          x_axis={{x_axis}},
                          input_go_type=go_type,
                          input_plot_title=go_type,
                          facet_by_column = {{facet_by_column}},
                          xaxis_levels = xaxis_levels,
                          scales=scales) } )

  names( list_of_heatmaps) <- annot_heat_map_ordered %>%
    distinct(  go_type) %>%
    mutate( output_name = go_type ) %>%
    pull(output_name)

  return(list_of_heatmaps)

}



#'@export
saveListOfFunctionalEnrichmentHeatmaps <- function(list_of_heatmaps,
                                                   results_dir,
                                                   file_name,
                                                   plot_width = 10,
                                                   plot_height = 10 ) {

  purrr::walk2(names( list_of_heatmaps),
               list_of_heatmaps,
               function(output_name, plot){
                 ggsave( filename=file.path( results_dir,
                                             paste0(file_name, "_", output_name   ,".pdf")),
                         plot=plot,
                         width = plot_width,
                         height=plot_height)   }  )

  purrr::walk2(names( list_of_heatmaps),
               list_of_heatmaps,
               function(output_name, plot){
                 ggsave( filename=file.path( results_dir,
                                             paste0(file_name, "_", output_name, ".png")),
                         plot=plot,
                         width = plot_width,
                         height=plot_height)   }  )


}

###------------------------------------------------------------------------------------------------------------------------

#' @export
enrichedPathwayBarPlot <- function( input_table, input_go_type = NA, remove_duplicted_entries = "merge") {

  duplicated_entries <- input_table %>%
    mutate(set_type = case_when( str_detect( gene_set, "positive") ~"positive",
                                 str_detect( gene_set, "negative") ~ "negative",
                                 TRUE ~ "neutral")) %>%
    group_by( across(c( any_of(added_columns), comparison, set_type, annotation_id) ) ) %>%
    dplyr::summarise( temp_qvalue = min(qvalue )) %>%
    ungroup() %>%
    dplyr::group_by( across(c( any_of(added_columns), comparison, annotation_id) ) ) %>%
    dplyr::summarise(counts = n(),
                     best_p_adj_value = min(temp_qvalue)) %>%
    ungroup() %>%
    dplyr::filter( counts > 1)

  if( remove_duplicted_entries == TRUE |
      remove_duplicted_entries == "delete" ) {
    input_table <- input_table  %>%
      anti_join( duplicated_entries, by =c("comparison" = "comparison",
                                           "annotation_id" = "annotation_id"))
  } else if( remove_duplicted_entries == "merge" ) {

    duplicates_tbl <- input_table %>%
      inner_join( duplicated_entries, by =c("comparison" = "comparison",
                                            "annotation_id" = "annotation_id")) %>%
      dplyr::filter( qvalue == best_p_adj_value ) %>%
      mutate( gene_set = "shared" )

    input_table <- input_table  %>%
      anti_join( duplicated_entries, by =c("comparison" = "comparison",
                                           "annotation_id" = "annotation_id")) %>%
      bind_rows( duplicates_tbl )

  }


  if (!is.na(input_go_type )) {
    if(! (input_go_type %in% (input_table %>% distinct(go_type) %>% pull(go_type)) ) ) {
      stop( paste0( "input_go_type = ", input_go_type, ", is not in the input_table" ) )
    }

    filt_input_table <-   input_table %>%
      dplyr::filter( go_type == input_go_type )
  } else {
    filt_input_table <- input_table
  }


  bar_plot_data <- filt_input_table  %>%
    mutate( neg_log10_qvalue =  -log10(qvalue))  %>%
    arrange( neg_log10_qvalue) %>%
    mutate( gene_set = case_when( gene_set == "positive_list" ~ "Up-regulated"
                                  , gene_set == "negative_list" ~ "Down-regulated"
                                  , gene_set == "shared" ~ "Both")) %>%
    group_by(gene_set, term ) %>%
    summarise( neg_log10_qvalue = max(neg_log10_qvalue)) %>%
    ungroup()

  term_ordering <- bar_plot_data %>%
    dplyr::select( gene_set, neg_log10_qvalue, term) %>%
    distinct(gene_set, neg_log10_qvalue, term) %>%
    arrange( (gene_set), (neg_log10_qvalue), term ) %>%
    pull(term )

  x_label <- "Pathway"
  if( !is.na(input_go_type) ) {
    x_label <- "GO Term"
  }

  bar_plot_data %>%
    mutate ( term = factor(term, levels= term_ordering)) %>%
    mutate ( gene_set = factor( gene_set, levels =c( "Up-regulated"
                                                     , "Down-regulated"
                                                     , "Both"))) %>%
    arrange( (gene_set), (neg_log10_qvalue) ) %>%
    ggplot( aes( term, neg_log10_qvalue, fill=gene_set  )) +
    geom_col() +
    scale_fill_manual(values = c( "Up-regulated" = "#F8766D", "Down-regulated" = "#00BFC4", "Both" = "grey")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))  +
    #facet_grid(  gene_set ~ ., scales = "free_y", space = "free_y") +
    coord_flip() +
    #theme(strip.text.y.right = element_text(angle = 0)) +
    labs(fill='Query Proteins') +
    ylab( expression("Significance, -log"[10]*"(q-value)") ) +
    xlab( x_label)

}



#'@description given input table, draw a bar plot representing the GO enrichment results.
#'The height of each bar represents the negative log (base 10) q-values of the query proteins.
#'@export
enrichedGoTermBarPlot <- function( input_table, output_dir,
                                   analysis_type = "GO", file_suffix, width=10, height = 7) {

  partial_go_term_bar_plot <- partial( enrichedPathwayBarPlot,
                                       input_table = filtered_enrich_revigo)

  list_of_go_type <- filtered_enrich_revigo %>%
    distinct( go_type) %>%
    arrange(go_type) %>%
    pull(go_type)

  list_of_barplots <- purrr::map( list_of_go_type,
                                  ~partial_go_term_bar_plot(input_go_type = .))

  names( list_of_barplots) <- list_of_go_type

  suffix_list <- file_suffix

  plotGOBarPlotWithSuffix <- function(suffix, input_plot, plot_name, width, height) {
    ggsave(plot = input_plot,
           filename=file.path( output_dir,
                               paste0(analysis_type, "_",
                                      plot_name,
                                      "_bar_plot.",
                                      suffix )),
           width= width,
           height = height) }

  plotBarPlot <- function(plot_name, input_plot, width, height ) {
    purrr::walk( suffix_list,
                 ~plotGOBarPlotWithSuffix(., input_plot, plot_name, width, height) ) }

  purrr::walk2(names( list_of_barplots), list_of_barplots,
               ~plotBarPlot(.x, .y, width=width, height=height) )


}

#'@description Create a word frequency distribution table for Word Cloud generation.
#'Based on article by Céline Van den Rul, How to Generate Word Clouds in R, Simple Steps on How and When to Use Them,
#' https://towardsdatascience.com/create-a-word-cloud-with-r-bde3e7422e8a (accessed 7th November 2022)
#'@export
#'@param text_list, a vector of text (e.g. a list of GO terms name)
createWordCloudDataFrame <- function( text_list) {

  docs <- Corpus(VectorSource(text_list))

  docs <- docs %>%
    tm_map(removeNumbers) %>%
    tm_map(removePunctuation) %>%
    tm_map(stripWhitespace)
  docs <- tm_map(docs, content_transformer(tolower))
  docs <- tm_map(docs, removeWords, stopwords("english"))

  dtm <- TermDocumentMatrix(docs)
  matrix <- as.matrix(dtm)
  words <- sort(rowSums(matrix),decreasing=TRUE)
  df <- data.frame(word = names(words),freq=words)
}

########################
