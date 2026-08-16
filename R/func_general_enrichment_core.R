# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# func_general_enrichment_core.R
# ============================================================================
# Purpose: Core general enrichment analysis and annotation helpers
# 
# This file owns shared GO enrichment, GSEA, annotation lookup, and table
# conversion helpers used across the omics workflows.
#
# Dependencies:
# - clusterProfiler, gProfiler2, GO.db, ReactomePA
# - func_general_plotting_*_helpers.R (for visualization)
# - func_general_helpers.R (for utility functions)
# ============================================================================

#'@export
parseNumList <-  function ( input_text ) {
  if( str_detect( input_text, "[.,;:]")) {
    str_split( input_text, "[.,;:]")[[1]] %>% purrr::map_int( as.integer)
  } else {
    return( as.integer( input_text ))
  }
}


# ----------------------------------------------------------------------------
# convertIdToAnnotation
# ----------------------------------------------------------------------------
#'@export
convertIdToAnnotation <- function( id, id_to_annotation_dictionary) {

    return( ifelse( !is.null(id_to_annotation_dictionary[[id]] ),
                    id_to_annotation_dictionary[[id]] ,
                    NA_character_))

}


# ----------------------------------------------------------------------------
# oneGoEnrichment
# ----------------------------------------------------------------------------
#'@title Run one GO enrichment
#'@export
#' @param go_annot,background_list,go_aspect,query_list,id_to_annotation_dictionary,annotation_id,protein_id,aspect_column,p_val_thresh,min_gene_set_size,max_gene_set_size,get_cluster_profiler_object Runtime inputs used by this function; see the usage section for accepted values.
oneGoEnrichment <- function( go_annot
                             , background_list
                             , go_aspect
                             , query_list
                             , id_to_annotation_dictionary
                             ,  annotation_id
                             , protein_id
                             , aspect_column
                             , p_val_thresh
                             , min_gene_set_size
                             , max_gene_set_size
                             , get_cluster_profiler_object = FALSE) {

  join_condition <- rlang::set_names( c( colnames(background_list)[1]),
                                      c( as_name(enquo( protein_id)) ) )

  if ( !is.na( go_aspect)) {
    go_annot_filt <- go_annot %>%
      dplyr::filter( {{aspect_column}} == go_aspect) |>
      mutate( {{protein_id}} := purrr::map_chr( {{protein_id}}, as.character))
  } else {
    go_annot_filt <- go_annot |>
      mutate( {{protein_id}} := purrr::map_chr( {{protein_id}}, as.character))
  }

  filtered_go_terms <- go_annot_filt %>%
    inner_join( background_list, by =join_condition )  %>%
    group_by( {{annotation_id}} ) %>%
    summarise(counts =n()) %>%
    ungroup() %>%
    arrange(desc(counts)) %>%
    dplyr::filter( counts <= max_gene_set_size & counts >= min_gene_set_size ) %>%
    dplyr::select(-counts)


  # print(head( filtered_go_terms))

  term_to_gene_tbl_filt <- go_annot_filt %>%
    inner_join( background_list, by =join_condition )  %>%
    dplyr::inner_join( filtered_go_terms
                       , by = as_name(enquo( annotation_id)) )  %>%
    dplyr::rename( gene = as_name(enquo(protein_id )) ,
                   term = as_name(enquo( annotation_id)) ) %>%
    dplyr::select(term, gene) %>%
    dplyr::distinct( term, gene )

  # print(nrow(term_to_gene_tbl_filt))

  ## Avoid singleton GO terms
  terms_to_avoid <- term_to_gene_tbl_filt %>%
    distinct() %>%
    dplyr::inner_join( data.frame(uniprot_acc = query_list)
                       , by=c("gene" = "uniprot_acc") )  %>%
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

  enrichment_result <- clusterProfiler::enricher(
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
                                              convertIdToAnnotation( id
                                                                     , id_to_annotation_dictionary) } )) %>%
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

    if( get_cluster_profiler_object == TRUE) {
      return( list( output_table = output_table_with_go_aspect
                    , cluster_profiler_object = enrichment_result))
    } else {
      return(output_table_with_go_aspect)
    }

  } else {

    return( NA)
  }

}


# ----------------------------------------------------------------------------
# runOneGoEnrichmentInOutFunction
# ----------------------------------------------------------------------------
#'@export
runOneGoEnrichmentInOutFunction <- function(significant_proteins,
                                            background_proteins,
                                            go_annotations,
                                            uniprot_data,
                                            p_val_thresh = 0.05,
                                            min_gene_set_size = 4,
                                            max_gene_set_size = 200,
                                            min_sig_gene_set_size = 2) {

  # Debug information
  print("Starting enrichment analysis")
  print(paste("Number of significant proteins:", length(significant_proteins)))
  print(paste("Number of background proteins:", length(background_proteins)))
  print("UniProt data columns:")
  print(colnames(uniprot_data))


  significant_df <- data.frame(uniprot_acc = significant_proteins)
  background_df <- data.frame(uniprot_acc = background_proteins)


  filtered_go_annotations <- go_annotations  |>
    inner_join( background_df
                , by =join_by( uniprot_acc == uniprot_acc  ) )  |>
    group_by( go_id ) |>
    dplyr::summarise(counts =n()) |>
    ungroup() |>
    dplyr::filter( counts <= max_gene_set_size & counts >= min_gene_set_size ) |>
    dplyr::select(-counts)

  ## There should be at least that many significnat proteins in a go term before it is accepted for analysis
  filtered_go_annotations_by_sig_proteins <- go_annotations |>
    inner_join( significant_df
                , by =join_by( uniprot_acc == uniprot_acc  ) )  |>
    group_by( go_id ) |>
    dplyr::summarise(counts =n()) |>
    ungroup() |>
    dplyr::filter( counts >= min_sig_gene_set_size ) |>
    dplyr::select(-counts)


  # Create term2gene and term2name for enricher
  term2gene <- go_annotations |>
    inner_join( filtered_go_annotations
                , by= join_by( go_id) ) |>
    inner_join( filtered_go_annotations_by_sig_proteins
                , by= join_by( go_id) ) |>
    dplyr::select(go_id, uniprot_acc, go_type) |>
    dplyr::distinct()

  term2name <- go_annotations |>
    inner_join( filtered_go_annotations
                , by= join_by( go_id) ) |>
    inner_join( filtered_go_annotations_by_sig_proteins
                , by= join_by( go_id) ) |>
    dplyr::select(go_id, go_term, go_type) |>
    dplyr::distinct()


  # Run enricher
  enricher_result <- purrr::map( c("Biological Process", "Cellular Compartment", "Molecular Function")
                                 , \(x) {
                                   clusterProfiler::enricher(
                                     gene = significant_proteins,
                                     universe = background_proteins,
                                     TERM2GENE = term2gene |>
                                       dplyr::filter( go_type == x) |>
                                       dplyr::select(go_id, uniprot_acc) ,
                                     TERM2NAME = term2name |>
                                       dplyr::filter( go_type == x) |>
                                       dplyr::select(go_id, go_term),
                                     pvalueCutoff = p_val_thresh,
                                     minGSSize = min_gene_set_size,
                                     maxGSSize = max_gene_set_size
                                   ) }) |>
    purrr::discard(is.null) |>  # Remove NULL results
    purrr::map(\(x) {
      if (is(x, "enrichResult") && nrow(x@result) > 0) {
        as.data.frame(x@result)
      } else {
        NULL
      }
    }) |>
    purrr::discard(is.null) |>  # Remove empty results
    purrr::reduce(bind_rows, .init = NULL)  # Safely combine results


  enricher_result_filt <- enricher_result
  # # Run Revigo
  # if( run_revigo == TRUE) {
  #
  #   revigo_input_list <- enricher_result |>
  #     distinct(ID) |>
  #     dplyr::pull( ID)
  #
  #   print("Running Revigo")
  #
  #   print( length(revigo_input_list))
  #   revigo_output <- queryRevigo( revigo_input_list,
  #                                 cutoff=revigo_cutoff,
  #                                 speciesTaxon = revigo_taxon,
  #                                 temp_file=NA)  |>
  #     dplyr::rename( go_id = "Term ID")  |>
  #     dplyr::filter( Eliminated == "False" |
  #                      is.na(Eliminated)) |>
  #     dplyr::distinct( go_id )
  #
  #   enricher_result_filt <- enricher_result |>
  #     inner_join( revigo_output, by=join_by( go_id))
  #
  # }

  # If no enrichment found, return NULL
  if (is.null(enricher_result_filt)) {
    print("No enrichment results found")
    return(NULL)
  }

  print("Processing enrichment results")

  # Get the column name for gene names in uniprot_data
  gene_name_col <- if ("Gene.Names" %in% colnames(uniprot_data)) {
    "Gene.Names"
  } else if ("GENENAME" %in% colnames(uniprot_data)) {
    "GENENAME"
  } else if ("gene_names" %in% colnames(uniprot_data)) {
    "gene_names"
  } else {
    stop("Could not find gene names column in uniprot_data. Available columns: ",
         paste(colnames(uniprot_data), collapse = ", "))
  }

  # Get the column name for UniProt accessions in uniprot_data
  uniprot_acc_col <- if ("Entry" %in% colnames(uniprot_data)) {
    "Entry"
  } else if ("UNIPROTKB" %in% colnames(uniprot_data)) {
    "UNIPROTKB"
  } else if ("uniprot_acc" %in% colnames(uniprot_data)) {
    "uniprot_acc"
  } else {
    stop("Could not find UniProt accession column in uniprot_data. Available columns: ",
         paste(colnames(uniprot_data), collapse = ", "))
  }

  message("List enricher_reult header")
  message( print( paste( colnames( enricher_result)) ))

  # Convert to data frame and add gene symbols
  enrichment_results <- enricher_result_filt |>
    dplyr::left_join(term2name, by = c("ID" = "go_id")) |>
    # Split gene list and get corresponding gene symbols
    dplyr::mutate(
      uniprot_list = purrr::map(
        geneID,
        \(x) str_split(x, "/")[[1]]
      )
    ) |>
    mutate( gene_names = purrr::map_chr(uniprot_list, \(acc_list) {
      genes <- uniprot_data |>
        dplyr::filter( Entry %in% acc_list) |>
        dplyr::mutate( gene_names_first = purrr::map_chr( gene_names, \(x) str_split(x, ";")[[1]][1])) |>
        dplyr::pull(gene_names_first)

      if (length(genes) == 0) {
        return(NA_character_)
      }

      unique_genes <- unique(unlist(stringr::str_split(genes, "\\s+")))
      paste(unique_genes[!is.na(unique_genes)], collapse = ";")
    })
    )

  print(paste("Found", nrow(enrichment_results), "enriched terms"))
  return(enrichment_results)
}


# ----------------------------------------------------------------------------
# convertProteinAccToGeneSymbol
# ----------------------------------------------------------------------------
#'@export
convertProteinAccToGeneSymbol <- function( gene_id_list, dictionary ) {

  purrr::map_chr( gene_id_list,
                  ~{ ifelse( . %in% names(dictionary ),
                             dictionary[[.]],
                             NA_character_)   } )  %>%
    paste( collapse="/")
}


# ----------------------------------------------------------------------------
# buildAnnotationIdToAnnotationNameDictionary
# ----------------------------------------------------------------------------
#'@export
buildAnnotationIdToAnnotationNameDictionary <- function(input_table, annotation_column, annotation_id_column) {

  id_to_annotation_dictionary <- NA

  dictionary_pair <- input_table %>%
    dplyr::filter( !is.na({{annotation_column}}) & !is.na({{annotation_id_column}})) %>%
    distinct({{annotation_column}},
             {{annotation_id_column}})

  id_to_annotation_dictionary <- dictionary_pair %>%
    dplyr::pull({{annotation_column}} )

  names(id_to_annotation_dictionary ) <-  dictionary_pair %>%
    dplyr::pull( {{annotation_id_column}})

  id_to_annotation_dictionary

}


# ----------------------------------------------------------------------------
# buildOneProteinToAnnotationList
# ----------------------------------------------------------------------------
#'@export
buildOneProteinToAnnotationList <- function( input_table, annotation_id, protein_id ) {

  temp_table <- input_table %>%
    dplyr::filter( !is.na( {{annotation_id}} ) & !is.na( {{protein_id}} )) %>%
    group_by( {{annotation_id}}) %>%
    nest( ) %>%
    ungroup()  %>%
    mutate( gene_set = purrr::map( data, ~{ (.) %>% dplyr::pull( {{protein_id}} )} ))

  gene_set_list <- temp_table %>% dplyr::pull(gene_set)

  names(gene_set_list ) <-temp_table %>% dplyr::pull({{annotation_id}})

  gene_set_list
}


# ----------------------------------------------------------------------------
# listifyTableByColumn
# ----------------------------------------------------------------------------
#'@export
listifyTableByColumn  <- function(input_table, column_name) {

  nested_table <- input_table %>%
    dplyr::filter(!is.na( {{column_name}})) %>%
    group_by({{column_name}}) %>%
    nest( ) %>%
    ungroup()

  list_of_tables <- nested_table %>%
    dplyr::pull(data)

  names( list_of_tables) <- nested_table %>%
    distinct({{column_name}})  %>%
    dplyr::pull( {{column_name}})

  list_of_tables
}


# ----------------------------------------------------------------------------
# runGsea
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
runGsea <- function(index_name, contrast_name, list_of_de_proteins, list_of_gene_sets, min_set_size = 4, max_set_size = 300) {

  gene_list <- list_of_de_proteins[[contrast_name]]

  msigdb_gene_set <- GSEABase::geneIds(list_of_gene_sets[[index_name]])

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


  gsea_results <- clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = as.data.frame(term_to_gene_tab_filt)
  )

  return(gsea_results)

}


# ----------------------------------------------------------------------------
# runEnricher
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
runEnricher <- function(index_name, contrast_name, list_of_de_proteins, list_of_gene_sets, min_set_size = 4, max_set_size = 300) {

  gene_list <- list_of_de_proteins[[contrast_name]]

  msigdb_gene_set <- GSEABase::geneIds(list_of_gene_sets[[index_name]])

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

  gsea_results <- clusterProfiler::enricher(
    gene = gene_list,
    TERM2GENE = as.data.frame(term_to_gene_tab_filt)
  )

  return(gsea_results)

}


# ----------------------------------------------------------------------------
# getUniprotAccToGeneSymbolDictionary
# ----------------------------------------------------------------------------
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
    dplyr::pull( gene_symbol)

  names( uniprot_to_gene_symbol_dict )  <- uniprot_to_gene_symbol %>%
    dplyr::pull( {{protein_id}} )

  uniprot_to_gene_symbol_dict

}


# ----------------------------------------------------------------------------
