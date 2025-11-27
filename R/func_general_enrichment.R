# ============================================================================
# func_general_enrichment.R
# ============================================================================
# Purpose: General enrichment analysis functions
# 
# This file contains general-purpose enrichment analysis functions used
# across all omics types, including GO enrichment, KEGG enrichment,
# Reactome enrichment, GSEA, and enrichment visualization. Functions in
# this file are used by enrichment modules across proteomics, metabolomics,
# and multiomics workflows.
#
# Functions to extract here:
# - GO enrichment functions (oneGoEnrichment, runOneGoEnrichmentInOutFunction)
# - KEGG enrichment functions (runKeggEnrichment)
# - Reactome enrichment functions (runReactomeEnrichment)
# - GSEA functions (runGsea)
# - Enricher functions (runEnricher)
# - Pathway enrichment functions (enrichProteinsPathways, etc.)
# - Revigo functions (queryRevigo, filterResultsWithRevigo)
# - Enrichment visualization functions (cnetplotEdited, etc.)
# - GO term conversion functions (goIdToTerm, uniprotGoIdToTerm, etc.)
# - Gene set functions (extract_geneSets, evaluateBestMinMaxGeneSetSize)
# - Enrichment result processing functions (cleanDuplicatesEnrichment, etc.)
# - Additional enrichment helper functions
#
# Dependencies:
# - clusterProfiler, gProfiler2, GO.db, ReactomePA
# - func_general_plotting.R (for visualization)
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# === GO Enrichment Functions ===

# Function 1: oneGoEnrichment()
# Current location: R/enrichment_functions.R
# Description: Runs one GO enrichment analysis
# oneGoEnrichment <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 2: runOneGoEnrichmentInOutFunction()
# Current location: R/enrichment_functions.R
# Description: Runs GO enrichment with input/output handling
# runOneGoEnrichmentInOutFunction <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === KEGG Enrichment Functions ===

# Function 3: runKeggEnrichment()
# Current location: R/multiomics_enrichment_functions.R
# Description: Runs KEGG pathway enrichment
# runKeggEnrichment <- function(...) {
#   # Extract from R/multiomics_enrichment_functions.R
# }

# === Reactome Enrichment Functions ===

# Function 4: runReactomeEnrichment()
# Current location: R/multiomics_enrichment_functions.R
# Description: Runs Reactome pathway enrichment
# runReactomeEnrichment <- function(...) {
#   # Extract from R/multiomics_enrichment_functions.R
# }

# === GSEA Functions ===

# Function 5: runGsea()
# Current location: R/enrichment_functions.R
# Description: Runs Gene Set Enrichment Analysis
# runGsea <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Enricher Functions ===

# Function 6: runEnricher()
# Current location: R/enrichment_functions.R
# Description: Runs enricher analysis
# runEnricher <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Pathway Enrichment Functions ===

# Function 7: enrichProteinsPathways()
# Current location: R/enrichment_functions.R
# Description: Enriches proteins pathways
# enrichProteinsPathways <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 8: enrichProteinsPathwaysHelper()
# Current location: R/enrichment_functions.R
# Description: Helper for enriching protein pathways
# enrichProteinsPathwaysHelper <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Revigo Functions ===

# Function 9: queryRevigo()
# Current location: R/enrichment_functions.R
# Description: Queries Revigo for GO term reduction
# queryRevigo <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 10: filterResultsWithRevigo()
# Current location: R/enrichment_functions.R
# Description: Filters enrichment results with Revigo
# filterResultsWithRevigo <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 11: filterResultsWithRevigoScholar()
# Current location: R/enrichment_functions.R
# Description: Filters results with Revigo (Scholar version)
# filterResultsWithRevigoScholar <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Enrichment Visualization Functions ===

# Function 12: clusterPathways()
# Current location: R/enrichment_functions.R
# Description: Clusters pathways
# clusterPathways <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 13: cnetplotEdited()
# Current location: R/enrichment_functions.R
# Description: Creates edited cnetplot
# cnetplotEdited <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 14: add_node_label()
# Current location: R/enrichment_functions.R
# Description: Adds node labels to plot
# add_node_label <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 15: node_add_alpha()
# Current location: R/enrichment_functions.R
# Description: Adds alpha to nodes
# node_add_alpha <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === GO Term Conversion Functions ===

# Function 16: goIdToTerm()
# Current location: R/enrichment_functions.R
# Description: Converts GO ID to term
# goIdToTerm <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 17: uniprotGoIdToTerm()
# Current location: R/enrichment_functions.R
# Description: Converts UniProt GO ID to term
# uniprotGoIdToTerm <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 18: uniprotGoIdToTermSimple()
# Current location: R/enrichment_functions.R
# Description: Simple UniProt GO ID to term conversion
# uniprotGoIdToTermSimple <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Gene Set Functions ===

# Function 19: extract_geneSets()
# Current location: R/enrichment_functions.R
# Description: Extracts gene sets
# extract_geneSets <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 20: evaluateBestMinMaxGeneSetSize()
# Current location: R/enrichment_functions.R
# Description: Evaluates best min/max gene set size
# evaluateBestMinMaxGeneSetSize <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Enrichment Result Processing Functions ===

# Function 21: cleanDuplicatesEnrichment()
# Current location: R/enrichment_functions.R
# Description: Cleans duplicate enrichment results
# cleanDuplicatesEnrichment <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 22: saveFilteredFunctionalEnrichmentTable()
# Current location: R/enrichment_functions.R
# Description: Saves filtered functional enrichment table
# saveFilteredFunctionalEnrichmentTable <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 23: readEnrichmentResultFiles()
# Current location: R/enrichment_functions.R
# Description: Reads enrichment result files
# readEnrichmentResultFiles <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Additional Enrichment Helper Functions ===

# Function 24: createDEResultsForEnrichment()
# Current location: R/functional_enrichment.R
# Description: Creates DE results for enrichment
# createDEResultsForEnrichment <- function(...) {
#   # Extract from R/functional_enrichment.R
# }

# Function 25: createEnrichmentResults()
# Current location: R/functional_enrichment.R
# Description: Creates enrichment results structure
# createEnrichmentResults <- function(...) {
#   # Extract from R/functional_enrichment.R
# }

# Function 26: perform_enrichment()
# Current location: R/functional_enrichment.R
# Description: Performs enrichment analysis
# perform_enrichment <- function(...) {
#   # Extract from R/functional_enrichment.R
# }

# Function 27: generate_enrichment_plots()
# Current location: R/functional_enrichment.R
# Description: Generates enrichment plots
# generate_enrichment_plots <- function(...) {
#   # Extract from R/functional_enrichment.R
# }

# Function 28: summarize_enrichment()
# Current location: R/functional_enrichment.R
# Description: Summarizes enrichment results
# summarize_enrichment <- function(...) {
#   # Extract from R/functional_enrichment.R
# }

# Function 29: processEnrichments()
# Current location: R/functional_enrichment.R
# Description: Processes enrichment results
# processEnrichments <- function(...) {
#   # Extract from R/functional_enrichment.R
# }

# Function 30: getEnrichmentResult()
# Current location: R/functional_enrichment.R
# Description: Gets enrichment result
# getEnrichmentResult <- function(...) {
#   # Extract from R/functional_enrichment.R
# }

# Function 31: getEnrichmentSummary()
# Current location: R/functional_enrichment.R
# Description: Gets enrichment summary
# getEnrichmentSummary <- function(...) {
#   # Extract from R/functional_enrichment.R
# }

# Function 32: getEnrichmentPlotly()
# Current location: R/functional_enrichment.R
# Description: Gets interactive enrichment plot
# getEnrichmentPlotly <- function(...) {
#   # Extract from R/functional_enrichment.R
# }

# Function 33: getEnrichmentHeatmap()
# Current location: R/enrichment_functions.R
# Description: Gets enrichment heatmap
# getEnrichmentHeatmap <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 34: drawListOfFunctionalEnrichmentHeatmaps()
# Current location: R/enrichment_functions.R
# Description: Draws list of functional enrichment heatmaps
# drawListOfFunctionalEnrichmentHeatmaps <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 35: drawListOfFunctionalEnrichmentHeatmapsScholar()
# Current location: R/enrichment_functions.R
# Description: Draws list of enrichment heatmaps (Scholar version)
# drawListOfFunctionalEnrichmentHeatmapsScholar <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 36: saveListOfFunctionalEnrichmentHeatmaps()
# Current location: R/enrichment_functions.R
# Description: Saves list of functional enrichment heatmaps
# saveListOfFunctionalEnrichmentHeatmaps <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 37: enrichedGoTermBarPlot()
# Current location: R/enrichment_functions.R
# Description: Creates enriched GO term bar plot
# enrichedGoTermBarPlot <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 38: enrichedPathwayBarPlot()
# Current location: R/enrichment_functions.R
# Description: Creates enriched pathway bar plot
# enrichedPathwayBarPlot <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 39: plotEnrichmentBarplot()
# Current location: R/enrichment_functions.R
# Description: Plots enrichment barplot
# plotEnrichmentBarplot <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 40: fc_readable()
# Current location: R/enrichment_functions.R
# Description: Makes fold change readable
# fc_readable <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 41: update_n()
# Current location: R/enrichment_functions.R
# Description: Updates n value for plots
# update_n <- function(...) {
#   # Extract from R/enrichment_functions.R
# }


# ----------------------------------------------------------------------------
# parseNumList
# ----------------------------------------------------------------------------
# Function used for parsing a list of minimum or maximum gene set size from command line
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
#'@param go_annot: Go annotation table.
#'@export
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


# ----------------------------------------------------------------------------
# runEnricher
# ----------------------------------------------------------------------------
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
# queryRevigo
# ----------------------------------------------------------------------------
#######################################################
### Query revigo
#'@export
queryRevigo <- function( input_list,
                         cutoff=0.5,
                         speciesTaxon = 10090,
                         temp_file=NA) {

  userData <-  paste(input_list,  collapse= "\n")

  print("userData")
  print(userData)

  httr::POST(
    url = "http://revigo.irb.hr/Revigo", # .aspx
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

  print("This is dat")
  print(dat)


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

  print("QueryRevigoExit")

  revigo_tbl
}


# ----------------------------------------------------------------------------
# clusterPathways
# ----------------------------------------------------------------------------
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
      anti_join( duplicated_entries
                 , by =c("comparison" = "comparison",
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


  if ( nrow( scores_for_clustering ) >= 2 ) {
    pathways_clustered <- hclust(dist(scores_for_clustering))

    pathways_sorting <- cutree(pathways_clustered, k=1:nrow(scores_for_clustering)) %>%
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
  } else {

    input_table %>%
      mutate( neg_log_p_value = -log10( p.adjust) )  %>%
      dplyr::select(  c(any_of(added_columns), comparison, annotation_id, term,  neg_log_p_value,  gene_set, go_type )) %>%
      mutate( annotation_id = as.character(annotation_id))
  }


}


# ----------------------------------------------------------------------------
# getEnrichmentHeatmap
# ----------------------------------------------------------------------------
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

  get_colour <- list( negative = "blue",
                      positive = "red",
                      negative_list = "blue", positive_list = "red",
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
      dplyr::pull({{x_axis}})


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

  if( length(which( c(24, 25)  %in% c(table_shape_colour %>% dplyr::pull( use_shape)) ) > 0 ) ) {
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


# ----------------------------------------------------------------------------
# readEnrichmentResultFiles
# ----------------------------------------------------------------------------
#' @export
readEnrichmentResultFiles <- function( table_of_files, file_names_column=file_name, go_type="KEGG") {

  list_of_files <- table_of_files %>%
    dplyr::pull( {{file_names_column}})

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


# ----------------------------------------------------------------------------
# filterResultsWithRevigo
# ----------------------------------------------------------------------------
#'@export
filterResultsWithRevigo <- function( enriched_results_tbl
                                     , added_columns
                                     , is_run_revigo=TRUE
                                     , revigo_cutoff=0.7
                                     , species_taxon = 9606 # Human
) {

  enrich_revigo <- NA

  if ( is_run_revigo == TRUE) {

    annotation_list <-  enriched_results_tbl %>%
      group_by( across( c(any_of(added_columns), comparison, gene_set, go_type) )) %>%
      nest() %>%
      ungroup() %>%
      mutate( annot_id_list = purrr::map( data, ~{ dplyr::pull(., annotation_id)} ))

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
      dplyr::select(-data, - annot_id_list)

    if(nrow(revigo_tbl) > 0 )  {
      revigo_tbl <- revigo_tbl %>%
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
      warning("filterResultsWithRevigo: Revigo summarizatio did not return any useful GO terms, return original input table.")
      enrich_revigo <- enriched_results_tbl
    }


  } else {

    enrich_revigo <- enriched_results_tbl
  }

  return( enrich_revigo)
}


# ----------------------------------------------------------------------------
# filterResultsWithRevigoScholar
# ----------------------------------------------------------------------------
#'@export
filterResultsWithRevigoScholar <- function( enriched_results_tbl
                                            , added_columns
                                            , is_run_revigo=TRUE
                                            , revigo_cutoff=0.7
                                            , species_taxon = 9606 # Human
) {

  enrich_revigo <- NA

  if ( is_run_revigo == TRUE) {

    annotation_list <-  enriched_results_tbl %>%
      group_by( across( c(any_of(added_columns),  go_type) )) %>%
      nest() %>%
      ungroup() %>%
      mutate( annot_id_list = purrr::map( data, ~{ dplyr::pull(., ID)} ))

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
      dplyr::select(-data, - annot_id_list)

    if(nrow(revigo_tbl) > 0 )  {
      revigo_tbl <- revigo_tbl %>%
        dplyr::rename(go_id = "Term ID")

      join_condition <- rlang::set_names(c("go_id", "go_type", added_columns),
                                         c("go_id", "go_type", added_columns))

      enrich_revigo <- enriched_results_tbl |>
        dplyr::rename( go_id = ID) |>
        left_join( revigo_tbl |>
                     dplyr::select(-Name),
                   by = join_condition) %>%
        dplyr::filter( Eliminated == "False" |
                         is.na(Eliminated))

    } else {
      warning("filterResultsWithRevigo: Revigo summarizatio did not return any useful GO terms, return original input table.")
      enrich_revigo <- enriched_results_tbl
    }


  } else {

    enrich_revigo <- enriched_results_tbl
  }

  return( enrich_revigo)
}


# ----------------------------------------------------------------------------
# saveFilteredFunctionalEnrichmentTable
# ----------------------------------------------------------------------------
#'@export
saveFilteredFunctionalEnrichmentTable <- function( enriched_results_tbl,
                                                   set_size_min,
                                                   set_size_max,
                                                   results_dir,
                                                   file_name,
                                                   list_of_columns_to_trim = c("gene_symbol")) {

  max_excel_cell_length <- 32760

  vroom::vroom_write( enriched_results_tbl %>%
                        dplyr::filter( min_set_size == set_size_min,
                                       max_set_size == set_size_max),
                      file.path(results_dir,
                                paste0( file_name, ".tab" )))

  writexl::write_xlsx( enriched_results_tbl %>%
                         dplyr::filter( min_set_size == set_size_min,
                                        max_set_size == set_size_max) %>%
                         mutate( across( one_of(list_of_columns_to_trim ), \(x)substr(x, 1, max_excel_cell_length)) ),
                       path=file.path(results_dir,
                                      paste0( file_name, ".xlsx" ) ))

  vroom::vroom_write( enriched_results_tbl ,
                      file.path(results_dir,
                                paste0( file_name, "_unfiltered.tab" )))

  writexl::write_xlsx( enriched_results_tbl %>%
                         mutate( across( one_of(list_of_columns_to_trim ), \(x)substr(x, 1, max_excel_cell_length)) ),
                       path=file.path(results_dir,
                                      paste0( file_name, "_unfiltered.xlsx" ) ))

}


# ----------------------------------------------------------------------------
# evaluateBestMinMaxGeneSetSize
# ----------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------
# drawListOfFunctionalEnrichmentHeatmaps
# ----------------------------------------------------------------------------
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

  added_columns <- unique( added_columns)

  input_table <- enriched_results_tbl |>
    distinct() %>%
    dplyr::filter( min_set_size == set_size_min,
                   max_set_size == set_size_max) %>%
    group_by(  across(  c(any_of(added_columns), comparison, gene_set, go_type)  )) %>%
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
    dplyr::pull(output_name)

  return(list_of_heatmaps)

}


# ----------------------------------------------------------------------------
# drawListOfFunctionalEnrichmentHeatmapsScholar
# ----------------------------------------------------------------------------
#'@export
drawListOfFunctionalEnrichmentHeatmapsScholar <- function(enriched_results_tbl,
                                                          added_columns,
                                                          x_axis = Analysis_Type,
                                                          analysis_column = Analysis_Type,
                                                          facet_by_column = NA,
                                                          remove_duplicted_entries = TRUE,
                                                          xaxis_levels=NA,
                                                          scales="fixed") {

  added_columns <- unique( added_columns)

  input_table <- enriched_results_tbl |>
    distinct() %>%
    group_by(  across(  c(any_of(added_columns), comparison, gene_set, go_type)  )) %>%
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
    dplyr::pull(output_name)

  return(list_of_heatmaps)

}


# ----------------------------------------------------------------------------
# saveListOfFunctionalEnrichmentHeatmaps
# ----------------------------------------------------------------------------
#'@export
saveListOfFunctionalEnrichmentHeatmaps <- function(list_of_heatmaps,
                                                   results_dir,
                                                   file_name,
                                                   plot_width = 10,
                                                   plot_height = 10 ) {

  if ( length(list_of_heatmaps) == length( plot_width)
       & length(list_of_heatmaps) == length( plot_height) ) {

  } else if( length( plot_width) == 1
             & length( plot_height) == 1) {
    plot_width <- rep(plot_width,  length(list_of_heatmaps))
    plot_height <- rep(plot_height,  length(list_of_heatmaps))

  } else {
    stop("Length of plot_width and plot_height should be one or same as the lenght of list of heatmaps.")
  }


  purrr::pwalk( list( output_name = names( list_of_heatmaps),
                      plot = list_of_heatmaps ,
                      plot_width = plot_width,
                      plot_height = plot_height),
                function(output_name, plot, plot_width, plot_height){
                  savePlot( plot=plot,
                            base_dir=results_dir,
                            plot_name = paste0( file_name, "_", output_name),
                            width = plot_width,
                            height=plot_height,
                            limitsize=FALSE) }
  )


}


# ----------------------------------------------------------------------------
# enrichedPathwayBarPlot
# ----------------------------------------------------------------------------
#' @export
enrichedPathwayBarPlot <- function( input_table, input_go_type = NA, remove_duplicted_entries = "merge", added_columns = "comparison") {

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
    if(! (input_go_type %in% (input_table %>% distinct(go_type) %>% dplyr::pull(go_type)) ) ) {
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
    dplyr::pull(term )

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


# ----------------------------------------------------------------------------
# enrichedGoTermBarPlot
# ----------------------------------------------------------------------------
#'@title Draw list of functional enrichment heatmaps
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
    dplyr::pull(go_type)

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


# ----------------------------------------------------------------------------
# createWordCloudDataFrame
# ----------------------------------------------------------------------------
#'@title Create a word frequency distribution table for Word Cloud generation.
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


# ----------------------------------------------------------------------------
# cleanDuplicatesEnrichment
# ----------------------------------------------------------------------------
#'@export
cleanDuplicatesEnrichment <- function( input_table
                                       , pathway_column = term
                                       , fdr_column = qvalue
                                       , gene_set_column = gene_set) {

  duplicated_entries <- input_table |>
    group_by( across(c( {{gene_set_column}}, {{pathway_column}}) ) ) |>
    dplyr::summarise( temp_qvalue = min({{fdr_column}} )) |>
    ungroup() |>
    dplyr::group_by( across(c( {{pathway_column}}) ) ) |>
    dplyr::summarise(counts = n(),
                     best_p_adj_value = min(temp_qvalue)) |>
    ungroup() |>
    dplyr::filter( counts > 1)

  duplicates_tbl <- input_table |>
    inner_join( duplicated_entries, by =join_by( {{pathway_column}} )) |>
    dplyr::filter( qvalue == best_p_adj_value ) |>
    mutate( gene_set = "Both" )

  input_table_cln <- input_table  %>%
    anti_join( duplicated_entries, by = join_by( {{pathway_column}} ) ) |>
    bind_rows( duplicates_tbl )

  positive_label <- input_table_cln |>
    distinct({{gene_set_column}}) |>
    dplyr::filter( str_detect( {{gene_set_column}}, "[P|p]ositive")) |>
    dplyr::pull({{gene_set_column}})
  negative_label <- input_table_cln |>
    distinct({{gene_set_column}}) |>
    dplyr::filter( str_detect( {{gene_set_column}}, "[N|n]egative")) |>
    dplyr::pull({{gene_set_column}})
  both_label <- input_table_cln |>
    distinct({{gene_set_column}}) |>
    dplyr::filter( str_detect( {{gene_set_column}}, "[B|b]oth")) |>
    dplyr::pull({{gene_set_column}})

  proteomics_go_helper <- input_table_cln |>
    mutate( neg_log_10_fdr = -log10({{fdr_column}}))    |>
    mutate( {{gene_set_column}} := factor( {{gene_set_column}}, levels=c( positive_label
                                                                          , negative_label
                                                                          , both_label ))) |>
    arrange( {{gene_set_column}}, desc({{gene_set_column}}), desc(neg_log_10_fdr))


   proteomics_go_helper
}


# ----------------------------------------------------------------------------
# plotEnrichmentBarplot
# ----------------------------------------------------------------------------
#'@export
plotEnrichmentBarplot <- function( input_table
                                   , pathway_column = term
                                   , fdr_column = qvalue
                                   , gene_set_column = gene_set
                                   , xlab_string = expression(-log[10](FDR))
                                   , ylab_string = "Enriched Terms"
                                   , legend_title = "Gene Set"
                                   , legend_colours = c( "red", "blue", "black")) {

  proteomics_go_helper <- cleanDuplicatesEnrichment( input_table
                                                                , pathway_column = {{pathway_column}}
                                                                , fdr_column = {{fdr_column}}
                                                                , gene_set_column = {{gene_set_column}} )

  proteomics_go_helper2 <- proteomics_go_helper |>
    mutate( {{pathway_column}} := factor( {{pathway_column}}, levels = rev( proteomics_go_helper |> distinct({{pathway_column}}) |> dplyr::pull( {{pathway_column}}) )) )

  # proteomics_go_helper_2 <- proteomics_go_helper |>
  #   mutate( {{pathway_column}} := factor( {{pathway_column}}, levels = rev( unique( proteomics_go_helper |> dplyr::pull( {{pathway_column}}) ))  ) )

  output_barplot <- proteomics_go_helper2 |>
    ggplot( aes( neg_log_10_fdr,   {{pathway_column}},  fill={{gene_set_column}})) +
    geom_bar( stat="identity" ) +
    scale_fill_manual(legend_title, values = legend_colours) +
    xlab ( xlab_string) +
    ylab ( ylab_string) +
    theme_bw() +
    scale_x_continuous(limits = c(0,-log10( min ( proteomics_go_helper|> dplyr::pull({{fdr_column}}) ) )*1.05), expand = c(0, 0))
  #+
    #facet_grid( rows = vars({{gene_set_column}})  , scales="free_y", space = "free_y")

  output_barplot
}


# ----------------------------------------------------------------------------
# list2df
# ----------------------------------------------------------------------------
#' @export
list2df <- function(inputList) {
  # ldf <- lapply(1:length(inputList), function(i) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })

  do.call('rbind', ldf)
}


# ----------------------------------------------------------------------------
# list2graph
# ----------------------------------------------------------------------------
#' @export
list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}


# ----------------------------------------------------------------------------
# get_param_change_message
# ----------------------------------------------------------------------------
#' @export
get_param_change_message <- function(parameter, params_df) {
  paste0("Use '", params_df[parameter, "listname"],
         " = list(", params_df[parameter, "present"],
         " = your_value)' instead of '", params_df[parameter, "original"],
         "'.\n The ", params_df[parameter, "original"],
         " parameter will be removed in the next version.")
}


# ----------------------------------------------------------------------------
# node_add_alpha
# ----------------------------------------------------------------------------
#' @export
node_add_alpha <- function(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight) {
  alpha_node <- rep(1, nrow(p$data))
  if (!is.null(hilight_category)) {
    alpha_node <- rep(alpha_nohilight, nrow(p$data))
    hilight_node <- c(hilight_category, hilight_gene)
    alpha_node[match(hilight_node, p$data$name)] <- alpha_hilight
  }
  p$data$alpha <- alpha_node
  return(p)
}


# ----------------------------------------------------------------------------
# get_enrichplot_color
# ----------------------------------------------------------------------------
#' @export
get_enrichplot_color <- function(n = 2) {
  colors <- getOption("enrichplot.colours")
  if (!is.null(colors)) return(colors)

  if (n != 2 && n != 3) stop("'n' should be 2 or 3")

  colors = c("#e06663", "#327eba")
  if (n == 2) return(colors)

  if (n == 3) return(c(colors[1], "white", colors[2]))
}


# ----------------------------------------------------------------------------
# set_enrichplot_color
# ----------------------------------------------------------------------------
#' @export
set_enrichplot_color <- function(colors = get_enrichplot_color(2),
                                 type = "color", name = NULL, .fun = NULL, ...) {

  type <- match.arg(type, c("color", "colour", "fill"))

  n <- length(colors)
  if (n < 2) {
    stop("'colors' should be of length >= 2")
  } else if (n == 2) {
    params <- list(low = colors[1], high = colors[2])
    fn_suffix <- "continuous"
  } else if (n == 3) {
    params <- list(low = colors[1], mid = colors[2], high = colors[3])
    fn_suffix <- "gradient2"
  } else {
    params <- list(colors = colors)
    fn_suffix <- "gradientn"
  }

  if (!is.null(.fun)) {
    if (n == 3) {
      # should determine parameter for user selected functions: 'gradient2' or 'gradientn'
      fn_type <- which_scale_fun(.fun)
      if (fn_type == "gradientn") {
        params <- list(colors = colors)
      } else {
        params <- list(low = colors[1], mid = colors[2], high = colors[3])
      }
    }
  } else {
    fn <- sprintf("scale_%s_%s", type, fn_suffix)
    .fun <- getFromNamespace(fn, "ggplot2")
  }

  params$guide <- guide_colorbar(reverse=TRUE, order=1)
  params$name <- name # no legend name setting by default as 'name = NULL'

  params <- modifyList(params, list(...))

  do.call(.fun, params)
}


# ----------------------------------------------------------------------------
# add_node_label
# ----------------------------------------------------------------------------
#' @export
add_node_label <- function(p, data, label_size_node, cex_label_node, shadowtext) {
  # If use 'aes_(alpha =~I(alpha))' will put an error for AsIs object.
  # But I(alpha) is necessory, so use 'alpha = I(data$alpha)'.
  segment.size <- get_ggrepel_segsize()
  if (is.null(data)) {
    data <- p$data
  }
  if (shadowtext) {
    p <- p + geom_node_text(aes_(label=~name), data = data,
                            alpha = I(data$alpha),
                            size = label_size_node * cex_label_node, bg.color = "white",
                            repel=TRUE, segment.size = segment.size)
  } else {
    p <- p + geom_node_text(aes_(label=~name), data = data,
                            alpha = I(data$alpha),
                            size = label_size_node * cex_label_node, repel=TRUE,
                            segment.size = segment.size)
  }
  return(p)
}


# ----------------------------------------------------------------------------
# get_ggrepel_segsize
# ----------------------------------------------------------------------------
#' @export
get_ggrepel_segsize <- function(default = 0.2) {
  getOption("ggrepel.segment.size", default = default)
}


# ----------------------------------------------------------------------------
# cnetplotEdited
# ----------------------------------------------------------------------------
#' @export
cnetplotEdited <- function(
    geneSets,
    showCategory = 5,
    foldChange = NULL,
    layout = "kk",
    colorEdge = FALSE,
    circular = FALSE,
    node_label = "all",
    cex_category = 1,
    cex_gene = 1,
    cex_label_category = 1,
    cex_label_gene = 1,
    color_category = "#E5C494",
    color_gene = "#B3B3B3",
    shadowtext = "all",
    color.params=list(
      foldChange = NULL,
      edge = FALSE,
      category = "#E5C494",
      gene = "#B3B3B3"
    ),
    cex.params=list(
      category_node = 1,
      gene_node = 1,
      category_label = 1,
      gene_label = 1
    ),
    hilight.params=list(
      category = NULL,
      alpha_hilight = 1,
      alpha_no_hilight = 0.3
    ),
    ...) {

  label_size_category <- 5
  label_size_gene <- 5
  node_label <- match.arg(node_label, c("category", "gene", "all", "none"))

  params_df <- as.data.frame(rbind(
    c("foldChange", "color.params", "foldChange"),
    c("colorEdge", "color.params", "edge"),
    c("color_category", "color.params", "category"),
    c("color_gene", "color.params", "gene"),

    c("cex_category", "cex.params", "category_node"),
    c("cex_gene", "cex.params", "gene_node"),
    c("cex_label_category", "cex.params", "category_label"),
    c("cex_label_gene", "cex.params", "gene_label"))
  )
  colnames(params_df) <- c("original", "listname", "present")
  rownames(params_df) <- params_df$original

  default.color.params <- list(
    foldChange = NULL,
    edge = FALSE,
    category = "#E5C494",
    gene = "#B3B3B3"
  )
  default.cex.params <- list(
    category_node = 1,
    gene_node = 1,
    category_label = 1,
    gene_label = 1
  )
  default.hilight.params <- list(
    category = NULL,
    alpha_hilight = 1,
    alpha_no_hilight = 0.3
  )

  # use modifyList to change the values of parameter
  color.params <- modifyList(default.color.params, color.params)
  cex.params <- modifyList(default.cex.params, cex.params)
  hilight.params <- modifyList(default.hilight.params, hilight.params)
  params_list <- list( showCategory = showCategory,
                       foldChange = foldChange,
                       layout = layout,
                       colorEdge = colorEdge,
                       circular = circular,
                       node_label = node_label,
                       cex_category = cex_category,
                       cex_gene = cex_gene,
                       cex_label_category = cex_label_category,
                       cex_label_gene = cex_label_gene,
                       color_category = color_category,
                       color_gene = color_gene,
                       shadowtext = shadowtext,
                       color.params = color.params,
                       cex.params = cex.params,
                       hilight.params = hilight.params
  )

  # get all parameters value
  args <- as.list(match.call())
  removed_params <- intersect(params_df$original, names(args))
  if (length(removed_params) > 0) {
    for (i in removed_params) {
      params_list[[params_df[i, 2]]][[params_df[i, 3]]] <- get(i)
      warn <- get_param_change_message(i, params_df)
      warning(warn)
    }
  }

  color.params <- params_list[["color.params"]]
  cex.params <- params_list[["cex.params"]]
  hilight.params <- params_list[["hilight.params"]]

  foldChange <- color.params[["foldChange"]]
  colorEdge <- color.params[["edge"]]
  color_category <- color.params[["category"]]
  color_gene <- color.params[["gene"]]

  cex_category <- cex.params[["category_node"]]
  cex_gene <- cex.params[["gene_node"]]
  cex_label_category <- cex.params[["category_label"]]
  cex_label_gene <- cex.params[["gene_label"]]

  hilight_category <- hilight.params[["category"]]
  alpha_hilight <- hilight.params[["alpha_hilight"]]
  alpha_nohilight <- hilight.params[["alpha_no_hilight"]]

  if (circular) {
    layout <- "linear"
    geom_edge <- geom_edge_arc
  } else {
    geom_edge <- geom_edge_link
  }
  if (is.logical(shadowtext)) {
    shadowtext <- ifelse(shadowtext, "all", "none")
  }
  shadowtext_category <- shadowtext_gene <- FALSE
  if (shadowtext == "all") shadowtext_category <- shadowtext_gene <- TRUE
  if (shadowtext == "category") shadowtext_category <- TRUE
  if (shadowtext == "gene") shadowtext_gene <- TRUE

  g <- list2graph(geneSets)

  # if (!inherits(x,  "list")) {
  #     foldChange <- fc_readable(x, foldChange)
  # }

  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2
  n <- length(geneSets)
  V(g)$size[1:n] <- size
  node_scales <- c(rep(cex_category, n), rep(cex_gene, (length(V(g)) - n)))

  # add edge alpha
  hilight_category <- intersect(hilight_category, names(geneSets))

  if (!is.null(hilight_category) && length(hilight_category) > 0) {
    edges <- attr(E(g), "vnames")
    E(g)$alpha <- rep(alpha_nohilight, length(E(g)))
    hilight_edge <- grep(paste(hilight_category, collapse = "|"), edges)
    hilight_gene <- edges[hilight_edge]
    hilight_gene <- gsub(".*\\|", "", hilight_gene)
    E(g)$alpha[hilight_edge] <- min(0.8, alpha_hilight)
  } else {
    E(g)$alpha <- rep(0.8, length(E(g)))
  }

  show_legend <- c(FALSE, TRUE)
  names(show_legend) <- c("alpha", "color")
  if (colorEdge) {
    E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    edge_layer <- geom_edge(aes_(color = ~category, alpha = ~I(alpha)),
                            show.legend = show_legend)
  } else {
    edge_layer <- geom_edge(aes_(alpha = ~I(alpha)), colour='darkgrey',
                            show.legend = FALSE)
  }

  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
    V(g)$color <- NA
    # V(g)$color[1:n] <- color_category
    V(g)$color[(n+1):length(V(g))] <- fc
    show_legend <- c(TRUE, FALSE)
    names(show_legend) <- c("color", "size")
    p <- ggraph(g, layout=layout, circular = circular)
    p$data[-(1:n), "size"] <- 3 * cex_gene

    p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)

    alpha_category <- c(rep(1, n), rep(0, nrow(p$data)-n))
    alpha_gene <- c(rep(0, n), rep(1, nrow(p$data)-n))

    if (!is.null(hilight_category) && length(hilight_category) > 0) {
      alpha_category <- c(rep(alpha_nohilight, n), rep(0, nrow(p$data)-n))
      alpha_gene <- c(rep(0, n), rep(alpha_nohilight, nrow(p$data)-n))
      alpha_gene[match(hilight_gene, p$data$name)] <- alpha_hilight
      alpha_gene[match(hilight_category, p$data$name)] <- alpha_hilight
    }

    p <- p + edge_layer +
      geom_node_point(aes_(size=~size), color=I(color_category),
                      data = NULL, show.legend = show_legend,
                      alpha = I(alpha_category)) +
      ggnewscale::new_scale_color() +
      geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size),
                      data = NULL, alpha = I(alpha_gene)) +
      scale_size(range=c(3, 8) * cex_category) +
      # scale_colour_gradient2(name = "fold change") +
      set_enrichplot_color(colors = rev(get_enrichplot_color(3)), name = "fold change")


  } else {
    V(g)$color <- color_gene
    V(g)$color[1:n] <- color_category
    p <- ggraph(g, layout=layout, circular=circular)
    p$data[-(1:n), "size"] <- 3 * cex_gene
    p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)
    p <- p + edge_layer +
      geom_node_point(aes_(color=~I(color), size=~size, alpha=~I(alpha)))+
      scale_size(range=c(3, 8) * cex_category)

  }

  p <- p + theme_void()

  if (node_label == "category") {
    p$data[-c(1:n), "name"] <- NA
    p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
                        cex_label_node = cex_label_category, shadowtext = shadowtext_category)
  } else if (node_label == "gene") {
    p$data[1:n, "name"] <- NA
    p <- add_node_label(p = p, data = NULL, label_size_node = label_size_gene,
                        cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
  } else if (node_label == "all") {
    p <- add_node_label(p = p, data = NULL,
                        label_size_node = c(rep(label_size_category, n), rep(label_size_gene, nrow(p$data)-n)),
                        cex_label_node = c(rep(cex_label_category, n), rep(cex_label_gene, nrow(p$data)-n)),
                        shadowtext = shadowtext_gene)
  }
  if (!is.null(foldChange)) {
    p <- p + guides(size  = guide_legend(order = 1),
                    color = guide_colorbar(order = 2))
  }
  return(p + guides(alpha = "none"))
}


# ----------------------------------------------------------------------------
# fc_readable
# ----------------------------------------------------------------------------
fc_readable <- function(x, foldChange = NULL) {
  if (is.null(foldChange))
    return(NULL)

  if (x@readable && x@keytype != "SYMBOL") {
    gid <- names(foldChange)
    if (is(x, 'gseaResult')) {
      ii <- gid %in% names(x@geneList)
    } else {
      ii <- gid %in% x@gene
    }
    gid[ii] <- x@gene2Symbol[gid[ii]]
    names(foldChange) <- gid
  }
  return(foldChange)
}


# ----------------------------------------------------------------------------
# update_n
# ----------------------------------------------------------------------------
#' @export
update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    if (inherits(x, 'list')) {
      showCategory <- showCategory[showCategory %in% names(x)]
    } else {
      showCategory <- intersect(showCategory, x$Description)
    }
    return(showCategory)
  }

  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (inherits(x, 'list')) {
    nn <- length(x)
  } else {
    nn <- nrow(x)
  }
  if (nn < n) {
    n <- nn
  }

  return(n)
}


# ----------------------------------------------------------------------------
# extract_geneSets
# ----------------------------------------------------------------------------
#' @export
extract_geneSets <- function(x, n) {
  n <- update_n(x, n)

  if (inherits(x, 'list')) {
    geneSets <- x
  } else {
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description
  }
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}


# ----------------------------------------------------------------------------
# enrichProteinsPathwaysHelper
# ----------------------------------------------------------------------------
#' Helper function to perform protein pathway enrichment analysis for a single contrast
#'
#' @description
#' This function performs pathway enrichment analysis on protein data using GO terms and gene symbols
#' downloaded from UniProt. The data is cached for future use to improve performance.
#'
#' @param de_analysis_results Output from deAnalysisWrapperFunction containing differential expression results
#' @param organism_taxid NCBI taxonomy ID for the organism (e.g., "9606" for human)
#' @param min_gene_set_size Minimum number of genes in a gene set (default: 4)
#' @param max_gene_set_size Maximum number of genes in a gene set (default: 200)
#' @param p_val_thresh P-value threshold for enrichment significance (default: 0.05)
#' @param protein_p_val_thresh P-value threshold for protein significance (default: 0.05)
#' @param cache_dir Directory to store cached UniProt data (default: "cache")
#' @param output_dir Directory for output files (default: "proteins_pathways_enricher")
#' @param use_cached Whether to use cached data if available (default: TRUE)
#' @param protein_id_delimiter Delimiter used in protein IDs (default: ":")
#' @param protein_id_column Name of the protein ID column (default: "Protein.Ids")
#'
#' @return A data frame containing enrichment results
#'
#' @import UniProt.ws
#' @import clusterProfiler
#' @import GO.db
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import plotly
#' @importFrom purrr map map_chr walk
#' @importFrom stringr str_split str_replace_all
#'
#' @export
enrichProteinsPathwaysHelper <- function(de_analysis_results,
                                  organism_taxid,
                                  min_gene_set_size = 4,
                                  max_gene_set_size = 200,
                                  p_val_thresh = 0.05,
                                  protein_p_val_thresh = 0.05,
                                  cache_dir = "cache",
                                  cache_file = "uniprot_annotations.RDS",
                                  output_dir = "proteins_pathways_enricher",
                                  use_cached = TRUE,
                                  protein_id_delimiter = ":",
                                  protein_id_column = "Protein.Ids"
                                  ) {

  # Create directories if they don't exist
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Initialize UniProt.ws handle
  up <- UniProt.ws(taxId = organism_taxid)

  # Cache file paths
  go_cache_file <- file.path(cache_dir, cache_file)

  # Extract protein IDs from all rows and clean them
  protein_data <- de_analysis_results$de_proteins_wide |>
    dplyr::select(!!sym(protein_id_column)) |>
    distinct() |>
    dplyr::mutate(!!sym(protein_id_column) := purrr::map_chr(!!sym(protein_id_column), \(x) {
      str_split(x, protein_id_delimiter)[[1]][1]
    })) |>
    dplyr::mutate(uniprot_acc = .cleanProteinIds(!!sym(protein_id_column)))

  # Get cached or fresh data
  uniprot_data <-  download_uniprot_data(protein_ids = protein_data
                                , cache_file = go_cache_file
                                , uniprot_handle = up
                                , protein_id_delimiter = protein_id_delimiter)
  # Process protein data for each comparison
  enrichment_results <- list()

    # Get significant proteins for this comparison and clean their IDs
    sig_data <- de_analysis_results$significant_rows |>
      dplyr::filter(analysis_type == "RUV applied") |>
      dplyr::mutate(!!sym(protein_id_column) := purrr::map_chr(!!sym(protein_id_column), \(x) {
        str_split(x, protein_id_delimiter)[[1]][1]
      })) |>
      dplyr::mutate(uniprot_acc = .cleanProteinIds(!!sym(protein_id_column)))

    # Get positive and negative proteins
    pos_prots <- sig_data |>
      dplyr::filter(fdr_qvalue < protein_p_val_thresh,
                   log2FC > 0) |>
      dplyr::pull(uniprot_acc)

    neg_prots <- sig_data |>
      dplyr::filter(fdr_qvalue < protein_p_val_thresh,
                   log2FC < 0) |>
      dplyr::pull(uniprot_acc)

    # Background proteins - use all proteins from the dataset
    background_proteins <- protein_data |>
      dplyr::pull(uniprot_acc)


    # Create GO annotation data frame from UniProt data
    uniprot_data_updated <- uniprot_data
    if( "Entry" %in% colnames(uniprot_data) & !("uniprot_acc" %in%  colnames(uniprot_data)) ) {
      uniprot_data_updated <- uniprot_data|>
        dplyr::rename( uniprot_acc = "Entry")
    }

   goterms <- AnnotationDbi::Term(GO.db::GOTERM)
   gotypes <- AnnotationDbi::Ontology(GO.db::GOTERM)

    go_annot_temp <- uniprot_data_updated |>
      pivot_longer( cols=matches("go_id"), names_to = "go_type", values_to = "go_id") |>
      dplyr::select(-go_type) |>
      separate_rows(go_id, sep = "; ")

    go_id_to_term <- go_annot_temp |>
      distinct(go_id)    |>
      dplyr::mutate(
        go_term = purrr::map_chr(go_id, \(x) {
          if (x %in% names(goterms)) {
            return(goterms[[x]])
          }
          return(NA_character_)
        }),
        go_type = purrr::map_chr(go_id, \(x) {
          if (x %in% names(gotypes)) {
            type <- gotypes[[x]]
            return(dplyr::case_when(
              type == "BP" ~ "Biological Process",
              type == "CC" ~ "Cellular Compartment",
              type == "MF" ~ "Molecular Function",
              TRUE ~ NA_character_
            ))
          }
          return(NA_character_)
        }) )

    go_annot <- go_annot_temp |>
      left_join( go_id_to_term, by = "go_id") |>
      dplyr::distinct(uniprot_acc, go_id, go_term, go_type)  |>
      dplyr::filter(!is.na(go_id))

  # Create empty result tibble
  empty_result <- tibble::tibble(
    ID = character(0),
    Description = character(0),
    GeneRatio = character(0),
    BgRatio = character(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    qvalue = numeric(0),
    geneID = character(0),
    Count = integer(0),
    term = character(0),
    go_type = character(0)
  )

  # Perform enrichment analysis
  pos_enrich <- runOneGoEnrichmentInOutFunction(
    significant_proteins = pos_prots,
    background_proteins = background_proteins,
    go_annotations = go_annot,
    uniprot_data = uniprot_data,
    p_val_thresh = p_val_thresh,
    min_gene_set_size = min_gene_set_size,
    max_gene_set_size = max_gene_set_size
  )

  neg_enrich <- runOneGoEnrichmentInOutFunction(
    significant_proteins = neg_prots,
    background_proteins = background_proteins,
    go_annotations = go_annot,
    uniprot_data = uniprot_data,
    p_val_thresh = p_val_thresh,
    min_gene_set_size = min_gene_set_size,
    max_gene_set_size = max_gene_set_size
  )

  # Handle NULL results
  pos_enrich <- if(is.null(pos_enrich)) empty_result else pos_enrich
  neg_enrich <- if(is.null(neg_enrich)) empty_result else neg_enrich

  # Add directionality
  enrichment_results <- bind_rows(
    pos_enrich |> mutate(directionality = "positive"),
    neg_enrich |> mutate(directionality = "negative")
  )

  return(enrichment_results)
}


# ----------------------------------------------------------------------------
# enrichProteinsPathways
# ----------------------------------------------------------------------------
#' Perform protein pathway enrichment analysis across multiple contrasts
#'
#' @description
#' This function performs pathway enrichment analysis across multiple contrasts in a proteomics dataset.
#' It processes each contrast separately and combines the results into a single data frame.
#'
#' @param de_analysis_results_list List of differential expression results for each contrast
#' @param taxon_id NCBI taxonomy ID for the organism (e.g., "9606" for human)
#' @param pathway_dir Directory for storing pathway analysis results
#' @param protein_id_delimiter Delimiter used in protein IDs (default: ":")
#' @param protein_p_val_thresh P-value threshold for protein significance (default: 0.05)
#' @param min_gene_set_size Minimum number of genes in a gene set (default: 4)
#' @param max_gene_set_size Maximum number of genes in a gene set (default: 200)
#' @param p_val_thresh P-value threshold for enrichment significance (default: 0.05)
#' @param cache_dir Directory to store cached UniProt data (default: "cache")
#' @param cache_file Name of the cache file for UniProt annotations (default: "uniprot_annotations.RDS")
#' @param use_cached Whether to use cached data if available (default: TRUE)
#'
#' @return A data frame containing combined enrichment results across all contrasts
#'
#' @import dplyr
#' @importFrom purrr map set_names
#'
#' @export
enrichProteinsPathways <- function(de_analysis_results_list,
                                 taxon_id,
                                 protein_id_delimiter = ":",
                                 protein_p_val_thresh = 0.05,
                                 min_gene_set_size = 4,
                                 max_gene_set_size = 200,
                                 p_val_thresh = 0.05,
                                 cache_dir = "cache",
                                 cache_file = "uniprot_annotations.RDS",
                                 use_cached = TRUE   ) {

  # Create a list to store all enrichment results
  all_enrichment_results_by_group <- names(de_analysis_results_list) |>
    purrr::set_names() |>  # Keep the contrast names
    purrr::map(\(contrast_name) {
      message(paste("Processing enrichment for contrast:", contrast_name))

      # Get the DE results for this contrast
      de_results <- de_analysis_results_list[[contrast_name]]

      # Run enrichment analysis
      enrichment_result <- enrichProteinsPathwaysHelper(
        de_analysis_results = de_results,
        organism_taxid = as.character(taxon_id),
        protein_p_val_thresh = protein_p_val_thresh,
        min_gene_set_size = min_gene_set_size,
        max_gene_set_size = max_gene_set_size,
        p_val_thresh = p_val_thresh,
        cache_dir = cache_dir,
        cache_file = cache_file,
        use_cached = use_cached,
        protein_id_delimiter = protein_id_delimiter )

      return(enrichment_result)
    })

  # Combine results from all contrasts
  go_results_table_by_group <- bind_rows(all_enrichment_results_by_group, .id = "comparison")

  return(go_results_table_by_group)
}


# ----------------------------------------------------------------------------
# download_uniprot_data
# ----------------------------------------------------------------------------
#'@export
download_uniprot_data <- function(protein_ids, cache_file, uniprot_handle, protein_id_delimiter = ":") {
  # Ensure cache directory exists
  cache_dir <- dirname(cache_file)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  # Check if cache exists
  if (file.exists(cache_file)) {
    tryCatch({
      cached_data <- readRDS(cache_file)

      # Find which proteins are not in cache
      missing_proteins <- protein_ids |>
        dplyr::filter(!uniprot_acc %in% cached_data$Entry)

      if (nrow(missing_proteins) == 0) {
        return(cached_data)
      }

      # Download only missing proteins
      message("Downloading data for ", nrow(missing_proteins), " new protein IDs")
      new_data <- batchQueryEvidenceGeneId(
        data.frame(uniprot_acc = unique(missing_proteins$uniprot_acc)  ),
        gene_id_column = "uniprot_acc",
        uniprot_handle = uniprot_handle,
        uniprot_columns = c(
          "protein_existence",
          "annotation_score",
          "reviewed",
          "gene_names",
          "protein_name",
          "length",
          "xref_ensembl",
          "go_id",  # This is the correct column name from UniProt.ws
          "keyword"
        )
      )



      # Process new data
      if (!is.null(new_data) && nrow(new_data) > 0) {

        new_data_processed <- new_data |>
          uniprotGoIdToTerm(
            uniprot_id_column = Entry,
            go_id_column = Gene.Ontology.IDs,
            sep = "; "
          ) |>
          dplyr::rename(
            Protein_existence = "Protein.existence",
            Protein_names = "Protein.names"
          ) |>
          mutate( Protein_existence = purrr::map_dbl(Protein_existence, as.double ))


        # Combine with cached data
        combined_data <- dplyr::bind_rows(cached_data, new_data_processed)

        # Try to update cache
        tryCatch({
          saveRDS(combined_data, cache_file)
        }, error = function(e) {
          warning("Could not update cache file: ", e$message)
        })

        return(combined_data)
      }
      return(cached_data)
    }, error = function(e) {
      warning("Error reading cache: ", e$message, ". Downloading all data fresh.")
      # Fall through to download all data
    })
  }

  # If no cache exists or cache read failed, download all data
  message("Downloading data for all protein IDs")

  all_data <- batchQueryEvidenceGeneId(
    data.frame(uniprot_acc = unique(protein_ids$uniprot_acc) ),
    gene_id_column = "uniprot_acc",
    uniprot_handle = uniprot_handle,
    uniprot_columns = c(
      "protein_existence",
      "annotation_score",
      "reviewed",
      "gene_names",
      "protein_name",
      "length",
      "xref_ensembl",
      "go_id",  # This is the correct column name from UniProt.ws
      "keyword"
    )
  )

  if (!is.null(all_data) && nrow(all_data) > 0) {
    processed_data <- all_data |>
      uniprotGoIdToTerm(
        uniprot_id_column = Entry,
        go_id_column = Gene.Ontology.IDs,
        sep = "; "
      ) |>
      dplyr::rename(
        Protein_existence = "Protein.existence",
        Protein_names = "Protein.names"
      )

    # Try to save to cache
    tryCatch({
      saveRDS(processed_data, cache_file)
    }, error = function(e) {
      warning("Could not save to cache file: ", e$message)
    })

    return(processed_data)
  }

  return(NULL)
}


# ----------------------------------------------------------------------------
# uniprotGoIdToTermSimple
# ----------------------------------------------------------------------------
#' Convert UniProt GO IDs to terms without grouping or pivoting
#' @param uniprot_dat  a table with uniprot accessions and a column with GO-ID
#' @param uniprot_id_column The name of the column with the uniprot accession, as a tidyverse header format, not a string
#' @param go_id_column The name of the column with the GO-ID, as a tidyverse header format, not a string
#' @param goterms Output from running \code{goterms <- Term(GOTERM)} from the GO.db library.
#' @param gotypes Output from running \code{gotypes <- Ontology(GOTERM)} from the GO.db library.
#' @return A table with columns for uniprot_id, go_id, go_term, and go_type
#' @export
uniprotGoIdToTermSimple <- function(uniprot_dat
                                   , uniprot_id_column = UNIPROTKB
                                   , go_id_column = `GO-IDs`
                                   , gene_name_column = Gene.Names
                                   , sep = "; "
                                   , goterms = AnnotationDbi::Term(GO.db::GOTERM)
                                   , gotypes = AnnotationDbi::Ontology(GO.db::GOTERM)) {

  print("Run uniprotGoIdToTermSimple")

  uniprot_acc_to_go_term <- uniprot_dat |>
    dplyr::distinct({{uniprot_id_column}}, {{gene_name_column}}, {{go_id_column}}) |>
    tidyr::separate_rows({{go_id_column}}, sep = sep) |>
    dplyr::filter(!is.na({{go_id_column}})) |>
    dplyr::mutate(
      go_term = purrr::map_chr({{go_id_column}}, \(x) {
        if (x %in% names(goterms)) {
          return(goterms[[x]])
        }
        return(NA_character_)
      }),
      go_type = purrr::map_chr({{go_id_column}}, \(x) {
        if (x %in% names(gotypes)) {
          type <- gotypes[[x]]
          return(dplyr::case_when(
            type == "BP" ~ "Biological Process",
            type == "CC" ~ "Cellular Compartment",
            type == "MF" ~ "Molecular Function",
            TRUE ~ NA_character_
          ))
        }
        return(NA_character_)
      })
    ) |>
    dplyr::filter(!is.na(go_term)) |>
    dplyr::rename(
      uniprot_acc = {{uniprot_id_column}},
      go_id = {{go_id_column}}
    )

  return(uniprot_acc_to_go_term)
}


# ----------------------------------------------------------------------------
# createDEResultsForEnrichment
# ----------------------------------------------------------------------------
#' Create DE Results For Enrichment
#'
#' @param contrasts_tbl A tibble containing contrast information
#' @param design_matrix A data frame containing the design matrix
#' @param de_output_dir Directory containing DE results files
#' @return An S4 object of class de_results_for_enrichment
#' @export
createDEResultsForEnrichment <- function(contrasts_tbl, design_matrix, de_output_dir) {
  # Helper function to format contrast filename
  format_contrast_filename <- function(contrast_string) {
    contrast_name <- stringr::str_split(contrast_string, "=")[[1]][1] |>
      stringr::str_replace_all("\\.", "_")

    paste0("de_proteins_", contrast_name, "_long_annot.tsv")
  }

  # Create new S4 object
  de_results <- new("de_results_for_enrichment")

  # Convert contrasts_tbl to tibble if it isn't already
  contrasts_tbl <- tibble::as_tibble(contrasts_tbl)

  # Fill slots
  de_results@contrasts <- contrasts_tbl
  de_results@design_matrix <- design_matrix
  de_results@de_data <- contrasts_tbl$contrasts |>
    purrr::set_names() |>
    purrr::map(function(contrast) {
      filename <- format_contrast_filename(contrast)
      filepath <- file.path(de_output_dir, filename)

      if (!file.exists(filepath)) {
        warning("File not found: ", filepath)
        return(NULL)
      }

      readr::read_tsv(filepath, show_col_types = FALSE)
    })

  return(de_results)
}


# ----------------------------------------------------------------------------
# createEnrichmentResults
# ----------------------------------------------------------------------------
# Constructor function
#' @export
createEnrichmentResults <- function(contrasts_tbl) {
  new("EnrichmentResults",
      contrasts = contrasts_tbl,
      enrichment_data = list(),
      enrichment_plots = list(),
      enrichment_plotly = list(),
      enrichment_summaries = list())
}


# ----------------------------------------------------------------------------
# perform_enrichment
# ----------------------------------------------------------------------------
#' @export
perform_enrichment <- function(data_subset,
                               species,
                               threshold,
                               sources,
                               domain_scope,
                               custom_bg,
                               exclude_iea = FALSE,
                               max_retries = 5,
                               wait_time = 5,
                               protein_id_column,
                               correction_method = "gSCS") {
  
  message("--- Entering perform_enrichment ---")
  message(sprintf("   perform_enrichment Arg: species = %s", species))
  message(sprintf("   perform_enrichment Arg: threshold = %s", threshold))
  message(sprintf("   perform_enrichment Arg: sources = %s", paste(sources, collapse = ", ")))
  message(sprintf("   perform_enrichment Arg: domain_scope = %s", domain_scope))
  message(sprintf("   perform_enrichment Arg: exclude_iea = %s", exclude_iea))
  message(sprintf("   perform_enrichment Arg: protein_id_column = %s", protein_id_column))
  message(sprintf("   perform_enrichment Arg: max_retries = %s", max_retries))
  
  message("      Data State (data_subset): Checking input data...")
  message(sprintf("      Data State (data_subset): Dims = %d rows, %d cols", nrow(data_subset), ncol(data_subset)))
  message("      Data State (data_subset) Structure:")
  utils::str(data_subset)
  message("      Data State (data_subset) Head:")
  print(head(data_subset))
  
  if (nrow(data_subset) == 0) {
    message("   perform_enrichment Step: Data subset is empty, returning NULL")
    message("--- Exiting perform_enrichment. Returning: NULL ---")
    return(NULL)
  }

  message("   perform_enrichment Step: Extracting protein IDs from data...")
  # Clean data before enrichment
  protein_ids <- data_subset[[protein_id_column]]
  
  message(sprintf("      Data State (protein_ids): Length = %d", length(protein_ids)))
  message(sprintf("      Data State (protein_ids): Class = %s", class(protein_ids)))
  message(sprintf("      Data State (protein_ids): First 5 values = %s", paste(head(protein_ids, 5), collapse = ", ")))
  
  # ✅ DEBUG66: Comprehensive NA analysis
  na_count_initial <- sum(is.na(protein_ids))
  empty_count_initial <- sum(protein_ids == "", na.rm = TRUE)
  valid_count_initial <- length(protein_ids) - na_count_initial - empty_count_initial
  
  message(sprintf("      ╔══════════════════════════════════════════════════════════════════╗"))
  message(sprintf("      ║ DEBUG66: PROTEIN ID AVAILABILITY ANALYSIS                       ║"))
  message(sprintf("      ╠══════════════════════════════════════════════════════════════════╣"))
  message(sprintf("      ║ Total proteins in data_subset:     %5d (100.0%%)                 ║", length(protein_ids)))
  message(sprintf("      ║ NA values in %s column:   %5d (%5.1f%%)                 ║", 
                 protein_id_column, na_count_initial, (na_count_initial/length(protein_ids))*100))
  message(sprintf("      ║ Empty string values:                %5d (%5.1f%%)                 ║", 
                 empty_count_initial, (empty_count_initial/length(protein_ids))*100))
  message(sprintf("      ║ Valid (non-NA, non-empty) values:  %5d (%5.1f%%)                 ║", 
                 valid_count_initial, (valid_count_initial/length(protein_ids))*100))
  message(sprintf("      ╚══════════════════════════════════════════════════════════════════╝"))

  if (any(is.na(protein_ids))) {
    na_count <- sum(is.na(protein_ids))
    message(sprintf("   ⚠️  WARNING: Filtering out %d proteins with NA %s values!", na_count, protein_id_column))
    message(sprintf("   ⚠️  %d proteins will be EXCLUDED from enrichment analysis!", na_count))
    
    if (na_count > (length(protein_ids) * 0.5)) {
      message("   ╔═══════════════════════════════════════════════════════════════════════════╗")
      message("   ║ ⚠️  CRITICAL WARNING: > 50% of proteins have NA gene names!              ║")
      message("   ║ Consider using 'Protein.Ids' instead of 'gene_name' for enrichment!      ║")
      message("   ║ OR ensure gene names are properly annotated in your data.                ║")
      message("   ╚═══════════════════════════════════════════════════════════════════════════╝")
    }
    
    warning(paste("NA values found in", protein_id_column, "column"))
    data_subset <- data_subset |> dplyr::filter(!is.na(.data[[protein_id_column]]))
    protein_ids <- data_subset[[protein_id_column]]
    message(sprintf("      Data State (protein_ids after NA filter): Length = %d", length(protein_ids)))
  }

  protein_ids <- unique(protein_ids[!is.na(protein_ids) & protein_ids != ""])
  
  message(sprintf("      DEBUG66: Final protein IDs to submit to gprofiler2: %d unique IDs", length(protein_ids)))
  if (length(protein_ids) > 0 && length(protein_ids) <= 20) {
    message(sprintf("      DEBUG66: All protein IDs being submitted: %s", paste(protein_ids, collapse = ", ")))
  } else if (length(protein_ids) > 20) {
    message(sprintf("      DEBUG66: First 20 protein IDs being submitted: %s", paste(head(protein_ids, 20), collapse = ", ")))
  }

  message("   perform_enrichment Step: Checking custom background...")
  message(sprintf("      Data State (custom_bg): Length = %d", length(custom_bg)))
  message(sprintf("      Data State (custom_bg): Class = %s", class(custom_bg)))
  message(sprintf("      Data State (custom_bg): First 5 values = %s", paste(head(custom_bg, 5), collapse = ", ")))

  if (any(is.na(custom_bg))) {
    na_bg_count <- sum(is.na(custom_bg))
    message(sprintf("   perform_enrichment Step: Found %d NA values in custom background, filtering...", na_bg_count))
    warning("NA values found in custom background IDs")
    custom_bg <- custom_bg[!is.na(custom_bg)]
    message(sprintf("      Data State (custom_bg after NA filter): Length = %d", length(custom_bg)))
  }

  custom_bg <- unique(custom_bg)

  result <- NULL
  attempt <- 1

  message(sprintf("   perform_enrichment Step: Starting gprofiler2 gost() retry loop (max %d attempts)...", max_retries))

  while (is.null(result) && attempt <= max_retries) {
    message(sprintf("   perform_enrichment Attempt %d: Calling gprofiler2::gost()...", attempt))
    
    tryCatch({
      message("   perform_enrichment Step: About to call gprofiler2::gost() with parameters:")
      message(sprintf("      gost query: %d protein IDs", length(protein_ids)))
      message(sprintf("      gost organism: %s", species))
      message(sprintf("      gost sources: %s", paste(sources, collapse = ", ")))
      message(sprintf("      gost user_threshold: %s", threshold))
      message(sprintf("      gost domain_scope: %s", domain_scope))
      message(sprintf("      gost custom_bg: %d background IDs", length(custom_bg)))
      message(sprintf("      gost exclude_iea: %s", exclude_iea))
      
      result <- gprofiler2::gost(
        query = protein_ids,
        organism = species,
        ordered_query = FALSE,
        sources = sources,
        user_threshold = threshold,
        correction_method = correction_method,
        exclude_iea = exclude_iea,
        evcodes = TRUE,
        domain_scope = domain_scope,
        custom_bg = custom_bg,
        significant = TRUE,
        highlight = TRUE
      )
      
      message("   perform_enrichment Step: gprofiler2::gost() completed successfully")
      message("      Data State (gost result): Checking gost result structure...")
      
      if (is.null(result)) {
        message("      Data State (gost result): Result is NULL")
      } else {
        message("      Data State (gost result) Structure:")
        utils::str(result)
        
        if ("result" %in% names(result)) {
          if (is.null(result$result)) {
            message("      Data State (gost result$result): Is NULL")
          } else {
            message(sprintf("      Data State (gost result$result): Dims = %d rows, %d cols", nrow(result$result), ncol(result$result)))
            message("      Data State (gost result$result) Column names:")
            print(names(result$result))
            if (nrow(result$result) > 0) {
              message("      Data State (gost result$result) Head:")
              print(head(result$result))
            } else {
              message("      Data State (gost result$result): NO ROWS - Empty results!")
            }
          }
        } else {
          message("      Data State (gost result): No 'result' component found")
        }
      }
      
    }, error = function(e) {
      message(sprintf("   perform_enrichment Attempt %d: ERROR in gprofiler2::gost(): %s", attempt, e$message))
      message(sprintf("   perform_enrichment Attempt %d: Will wait %d seconds before retry...", attempt, wait_time))
      Sys.sleep(wait_time)
    })
    
    attempt <- attempt + 1
  }

  if (is.null(result)) {
    message("   perform_enrichment Step: All retry attempts failed, returning NULL")
    message("--- Exiting perform_enrichment. Returning: NULL ---")
    return(NULL)
  }

  message("   perform_enrichment Step: Successfully obtained gprofiler2 result")
  message("--- Exiting perform_enrichment. Returning: gost result object ---")
  return(result)
}


# ----------------------------------------------------------------------------
# generate_enrichment_plots
# ----------------------------------------------------------------------------
# Plot generation function
#' @export
generate_enrichment_plots <- function(enrichment_result, contrast, direction, pathway_dir) {
  # Defensive check for empty or NULL results
  if (is.null(enrichment_result) || is.null(enrichment_result$result) || nrow(enrichment_result$result) == 0) {
    return(list(static = NULL, interactive = NULL))
  }

  # Extract the significance threshold from the result object metadata
  significance_threshold <- enrichment_result$meta$query_metadata$user_threshold

  # Prepare data for plotting
  plot_data <- enrichment_result$result %>%
    dplyr::mutate(
      neg_log10_p = -log10(p_value),
      # Ensure 'source' is a factor with a consistent level order for plotting
      source = factor(source, levels = c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC"))
    ) %>%
    # Drop any levels that are not actually present in the data to avoid empty spaces on the plot
    dplyr::mutate(source = forcats::fct_drop(source))

  # Identify the top significant terms to add labels for
  top_terms <- plot_data %>%
    dplyr::group_by(source) %>%
    dplyr::arrange(p_value) %>%
    dplyr::slice_head(n = 3) %>%
    dplyr::ungroup()

  # Create the static ggplot object
  static_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = source, y = neg_log10_p, text = paste0(
    "Term: ", term_name, "\n",
    "ID: ", term_id, "\n",
    "P-value: ", signif(p_value, 3), "\n",
    "Genes: ", intersection_size
  ))) +
    # Add a line for the significance threshold
    ggplot2::geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "red") +
    ggplot2::geom_jitter(ggplot2::aes(color = source, size = term_size), alpha = 0.7, width = 0.2) +
    ggrepel::geom_text_repel(
      data = top_terms,
      ggplot2::aes(label = term_name),
      size = 3, max.overlaps = 15, box.padding = 0.5, point.padding = 0.3, force = 5
    ) +
    ggplot2::scale_color_manual(
      values = c(`GO:BP` = "#ff9900", `GO:CC` = "#109618", `GO:MF` = "#dc3912", 
                 KEGG = "#dd4477", REAC = "#3366cc"),
      name = "Source", drop = FALSE 
    ) +
    ggplot2::scale_size_continuous(name = "Term Size", range = c(3, 10)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = paste0("Enrichment for ", contrast, " (", direction, "-regulated)"),
      x = "Annotation Source", y = "-log10(p-value)"
    ) +
    ggplot2::guides(color = "none")

  # Convert to an interactive plotly object
  interactive_plot <- tryCatch({
      plotly::ggplotly(static_plot, tooltip = "text")
  }, error = function(e) {
      warning(paste("Failed to convert ggplot to plotly for", contrast, direction, ":", e$message))
      return(NULL)
  })

  # Save the results data table
  tryCatch({
    result_table <- enrichment_result$result
    result_table$parents <- sapply(result_table$parents, paste, collapse = ", ")
    readr::write_tsv(result_table,
      file = file.path(pathway_dir, paste0(contrast, "_", direction, "_enrichment_results.tsv"))
    )
  }, error = function(e) {
    warning(paste("Failed to write enrichment results table for", contrast, direction, ":", e$message))
  })
  
  return(list(
    static = static_plot,
    interactive = interactive_plot
  ))
}


# ----------------------------------------------------------------------------
# summarize_enrichment
# ----------------------------------------------------------------------------
# Summary function
#' @export
summarize_enrichment <- function(enrichment_result) {
  if (is.null(enrichment_result) || length(enrichment_result$result) == 0) {
    return(data.frame(
      total = 0,
      GO_BP = 0,
      GO_CC = 0,
      GO_MF = 0,
      KEGG = 0,
      REAC = 0
    ))
  }

  result_df <- enrichment_result$result

  data.frame(
    total = nrow(result_df),
    GO_BP = sum(result_df$source == "GO:BP"),
    GO_CC = sum(result_df$source == "GO:CC"),
    GO_MF = sum(result_df$source == "GO:MF"),
    KEGG = sum(result_df$source == "KEGG"),
    REAC = sum(result_df$source == "REAC")
  )
}


# ----------------------------------------------------------------------------
# processEnrichments
# ----------------------------------------------------------------------------
#' Process Enrichments
#'
#' @param de_results S4 object containing differential expression results
#' @param taxon_id NCBI taxonomy ID for the organism
#' @param up_cutoff Log2 fold change cutoff for up-regulated proteins (default: 0)
#' @param down_cutoff Log2 fold change cutoff for down-regulated proteins (default: 0)
#' @param q_cutoff FDR q-value threshold for enrichment significance (default: 0.05)
#' @param pathway_dir Directory for saving pathway results
#' @param go_annotations UniProt GO annotations (required for unsupported organisms)
#' @param exclude_iea Whether to exclude IEA (Inferred Electronic Annotation) terms (default: FALSE)
#' @param protein_id_column Name of the protein ID column (default: "Protein.IDs")
#' @param contrast_names Vector of contrast names for output labeling
#' @param correction_method Method for FDR correction (default: "gSCS")
#'
#' @return S4 EnrichmentResults object containing enrichment data, plots, and summaries
#'
#' @export
processEnrichments <- function(de_results,
                               taxon_id,
                               up_cutoff = 0,
                               down_cutoff = 0,
                               q_cutoff = 0.05,
                               pathway_dir,
                               go_annotations = NULL,
                               exclude_iea = FALSE,
                               protein_id_column = "Protein.IDs",
                               contrast_names = NULL,
                               correction_method = "gSCS") {

  message("--- RUNNING processEnrichments VERSION [Timestamp: ", Sys.time(), "] ---")

  # Validate exclude_iea parameter
  if (is.null(exclude_iea)) {
    stop("exclude_iea must be explicitly set to TRUE or FALSE")
  }

  if (!is.logical(exclude_iea)) {
    stop("exclude_iea must be a logical value (TRUE or FALSE)")
  }

  # Common model organisms lookup
  supported_organisms <- tibble::tribble(
    ~taxid,     ~id,            ~name,
    "9606",     "hsapiens",     "Homo sapiens",
    "10090",    "mmusculus",    "Mus musculus",
    "10116",    "rnorvegicus",  "Rattus norvegicus",
    "7227",     "dmelanogaster", "Drosophila melanogaster",
    "6239",     "celegans",     "Caenorhabditis elegans",
    "4932",     "scerevisiae",  "Saccharomyces cerevisiae",
    "3702",     "athaliana",    "Arabidopsis thaliana",
    "7955",     "drerio",       "Danio rerio",
    "9031",     "ggallus",      "Gallus gallus",
    "9823",     "sscrofa",      "Sus scrofa",
    "9913",     "btaurus",      "Bos taurus",
    "9544",     "mmulatta",     "Macaca mulatta",
    "9598",     "ptroglodytes", "Pan troglodytes"
  )

  is_supported <- as.character(taxon_id) %in% supported_organisms$taxid

  if(is_supported) {
    message(sprintf("Taxon ID %s found in supported organisms. Proceeding with gprofiler2 analysis...", taxon_id))

    # Convert taxon_id to species
    species <- supported_organisms |>
      dplyr::filter(.data$taxid == as.character(taxon_id)) |>
      dplyr::pull(.data$id)

    enrichment_results <- createEnrichmentResults(de_results@contrasts)

    # Process each contrast
    results <- de_results@de_data |>
      purrr::map(function(de_data) {
        tryCatch({
        if(is.null(de_data)) {
          warning("No DE data found for contrast")
          return(NULL)
        }

        # Split data into up/down regulated
        message(sprintf("Total proteins before filtering: %d", nrow(de_data)))
        
        # ✅ DIAGNOSTIC: Check protein_id_column and data structure
        message(sprintf("DEBUG: protein_id_column = '%s'", protein_id_column))
        message(sprintf("DEBUG: Available columns in de_data: %s", paste(names(de_data), collapse = ", ")))
        
        if (!protein_id_column %in% names(de_data)) {
          stop(sprintf("ERROR: Column '%s' not found in DE data. Available columns: %s", 
                      protein_id_column, paste(names(de_data), collapse = ", ")))
        }
        
        # Check if the column has valid data
        protein_col_data <- de_data[[protein_id_column]]
        message(sprintf("DEBUG: First 3 values in protein_id_column: %s", 
                       paste(head(protein_col_data, 3), collapse = ", ")))
        message(sprintf("DEBUG: Column class: %s", class(protein_col_data)))

                 message("   processEnrichments Step: About to apply protein ID splitting and FDR filtering...")
         message(sprintf("      Data State (de_data Before Modify): Dims=%d rows, %d cols. First 5 IDs: %s", nrow(de_data), ncol(de_data), paste(head(de_data[[protein_id_column]], 5), collapse=", ")))
         
         # ✅ DEBUG66: Check NA values BEFORE filtering
         na_count_before <- sum(is.na(de_data[[protein_id_column]]))
         message(sprintf("      DEBUG66: NA values in %s column BEFORE filtering: %d out of %d (%.1f%%)", 
                        protein_id_column, na_count_before, nrow(de_data), (na_count_before/nrow(de_data))*100))
         
         # ✅ DEBUG66: Check how many pass q-value filter (before NA filter)
         passes_q_filter <- de_data |> dplyr::filter(.data$fdr_qvalue < q_cutoff)
         message(sprintf("      DEBUG66: Proteins passing q-value filter (< %.3f): %d", q_cutoff, nrow(passes_q_filter)))
         if (nrow(passes_q_filter) > 0) {
           na_in_sig <- sum(is.na(passes_q_filter[[protein_id_column]]))
           message(sprintf("      DEBUG66: NA values in %s among q-significant proteins: %d out of %d (%.1f%%)",
                          protein_id_column, na_in_sig, nrow(passes_q_filter), (na_in_sig/nrow(passes_q_filter))*100))
           message(sprintf("      DEBUG66: First 10 %s values in q-significant proteins: %s",
                          protein_id_column, paste(head(passes_q_filter[[protein_id_column]], 10), collapse = ", ")))
         }

         subset_sig <- de_data |>
            dplyr::mutate(
              !!rlang::sym(protein_id_column) := stringr::str_remove(.data[[protein_id_column]], ":.*")
            ) |>
           dplyr::filter(.data$fdr_qvalue < q_cutoff)
         
         # ✅ DEBUG66: Check if NAs are being filtered somewhere
         na_count_after_q <- sum(is.na(subset_sig[[protein_id_column]]))
         message(sprintf("      DEBUG66: NA values in %s AFTER q-filter (before any NA removal): %d out of %d",
                        protein_id_column, na_count_after_q, nrow(subset_sig)))

         message("   processEnrichments Step: Protein ID splitting and filtering complete.")
         message(sprintf("      Data State (subset_sig After Modify): Dims=%d rows, %d cols.", nrow(subset_sig), ncol(subset_sig)))
         if(nrow(subset_sig) > 0) {
            message(sprintf("      Data State (subset_sig After Modify): First 5 IDs: %s", paste(head(subset_sig[[protein_id_column]], 5), collapse=", ")))
            # Also check for the literal "Protein.Ids" to be sure
            if(any(subset_sig[[protein_id_column]] == "Protein.Ids", na.rm = TRUE)) {
              message("      WARNING: Literal 'Protein.Ids' found in the protein ID column after modification!")
            }
         }

        message(sprintf("Proteins passing q-value cutoff (%.3f): %d", q_cutoff, nrow(subset_sig)))
        
        # ✅ DEBUG 66: Extensive logging for subset_sig
        message("      Data State (subset_sig) Structure:")
        utils::str(subset_sig)
        if (nrow(subset_sig) > 0) {
          message("      Data State (subset_sig) Head:")
          print(head(subset_sig))
          message("      Data State (subset_sig fdr_qvalue range):")
          print(range(subset_sig$fdr_qvalue, na.rm = TRUE))
          message("      Data State (subset_sig log2FC range):")
          print(range(subset_sig$log2FC, na.rm = TRUE))
        }

        # No longer need subset_for_enrichment, as subset_sig is already filtered
        message(sprintf("Proteins available for enrichment: %d", nrow(subset_sig)))
        
        # ✅ DEBUG 66: Check subset_for_enrichment
        if (nrow(subset_sig) == 0) {
          message("      Data State (subset_sig): NO PROTEINS PASS ENRICHMENT CUTOFF!")
        }

        up_matrix <- subset_sig |>
          dplyr::filter(.data$log2FC > up_cutoff)
        message(sprintf("Up-regulated proteins (log2FC > %g): %d", up_cutoff, nrow(up_matrix)))
        
        # ✅ DEBUG66: Check NA values in up_matrix
        if (nrow(up_matrix) > 0) {
          na_in_up <- sum(is.na(up_matrix[[protein_id_column]]))
          message(sprintf("      DEBUG66: NA values in %s for UP-regulated proteins: %d out of %d (%.1f%%)",
                         protein_id_column, na_in_up, nrow(up_matrix), (na_in_up/nrow(up_matrix))*100))
          message(sprintf("      DEBUG66: First 10 UP gene names: %s",
                         paste(head(up_matrix[[protein_id_column]], 10), collapse = ", ")))
        }
        
        # ✅ DEBUG 66: Check up_matrix details
        if (nrow(up_matrix) > 0) {
          message("      Data State (up_matrix) Structure:")
          utils::str(up_matrix)
          message("      Data State (up_matrix) Head:")
          print(head(up_matrix))
          message("      Data State (up_matrix log2FC values):")
          print(summary(up_matrix$log2FC))
        } else {
          message("      Data State (up_matrix): NO UP-REGULATED PROTEINS FOUND!")
          message(sprintf("      Debug: up_cutoff = %g, max log2FC in subset_sig = %g", 
                         up_cutoff, if(nrow(subset_sig) > 0) max(subset_sig$log2FC, na.rm = TRUE) else NA))
        }

        down_matrix <- subset_sig |>
          dplyr::filter(.data$log2FC < -down_cutoff)
        message(sprintf("Down-regulated proteins (log2FC < -%g): %d", down_cutoff, nrow(down_matrix)))
        
        # ✅ DEBUG66: Check NA values in down_matrix
        if (nrow(down_matrix) > 0) {
          na_in_down <- sum(is.na(down_matrix[[protein_id_column]]))
          message(sprintf("      DEBUG66: NA values in %s for DOWN-regulated proteins: %d out of %d (%.1f%%)",
                         protein_id_column, na_in_down, nrow(down_matrix), (na_in_down/nrow(down_matrix))*100))
          message(sprintf("      DEBUG66: First 10 DOWN gene names: %s",
                         paste(head(down_matrix[[protein_id_column]], 10), collapse = ", ")))
        }
        
        # ✅ DEBUG 66: Check down_matrix details
        if (nrow(down_matrix) > 0) {
          message("      Data State (down_matrix) Structure:")
          utils::str(down_matrix)
          message("      Data State (down_matrix) Head:")
          print(head(down_matrix))
          message("      Data State (down_matrix log2FC values):")
          print(summary(down_matrix$log2FC))
        } else {
          message("      Data State (down_matrix): NO DOWN-REGULATED PROTEINS FOUND!")
          message(sprintf("      Debug: down_cutoff = %g, min log2FC in subset_sig = %g", 
                         down_cutoff, if(nrow(subset_sig) > 0) min(subset_sig$log2FC, na.rm = TRUE) else NA))
        }

        # Get background IDs from the full de_data for this contrast
        custom_bg <- de_data[[protein_id_column]] |>
          unique()

        message(sprintf("Using %d unique proteins as background for enrichment analysis", length(custom_bg)))
        
        # ✅ DEBUG 66: Check background details
        message("      Data State (custom_bg): First 5 background proteins:")
        message(paste(head(custom_bg, 5), collapse = ", "))

        # Process up and down regulated genes
        list(
          up = tryCatch({
            if(nrow(up_matrix) > 0) {
              # ✅ DEBUG 66: Log call to perform_enrichment for up-regulated
              message("   processEnrichments Step: CALLING perform_enrichment for UP-REGULATED proteins...")
              message(sprintf("      About to analyze %d up-regulated proteins", nrow(up_matrix)))
              
              protein_col <- protein_id_column
              message(sprintf("      Using protein_id_column: %s", protein_col))
              message(sprintf("      Using species: %s", species))
              message(sprintf("      Using threshold: %s", q_cutoff))
              message(sprintf("      Using exclude_iea: %s", exclude_iea))
              
              up_result <- perform_enrichment(
                data_subset = up_matrix,
                species = species,
                threshold = q_cutoff,
                sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                domain_scope = "custom",
                custom_bg = custom_bg,
                exclude_iea = exclude_iea,
                protein_id_column = protein_col,
                correction_method = correction_method
              )
              
              # ✅ DEBUG 66: Log result from perform_enrichment
              message("   processEnrichments Step: UP-REGULATED enrichment completed")
              message("      Data State (up_result): Checking perform_enrichment result...")
              
              if (is.null(up_result)) {
                message("      Data State (up_result): Result is NULL")
              } else {
                message("      Data State (up_result) Structure:")
                utils::str(up_result)
                
                if ("result" %in% names(up_result) && !is.null(up_result$result)) {
                  message(sprintf("      Data State (up_result$result): Found %d enriched terms", nrow(up_result$result)))
                  if (nrow(up_result$result) > 0) {
                    message("      Data State (up_result$result) Head:")
                    print(head(up_result$result))
                  }
                } else {
                  message("      Data State (up_result): No 'result' component or result is NULL")
                }
              }
              
              up_result
            } else {
              message("   processEnrichments Step: SKIPPING UP-REGULATED enrichment (no proteins)")
              NULL
            }
          }, error = function(e) {
            message(sprintf("   processEnrichments Step: ERROR in up-regulated enrichment: %s", e$message))
            warning(sprintf("Error processing up-regulated genes: %s", e$message))
            message("Debug info for up-regulated genes:")
            message("Number of proteins: ", nrow(up_matrix))
            NULL
          }),

          down = tryCatch({
            if(nrow(down_matrix) > 0) {
              # ✅ DEBUG 66: Log call to perform_enrichment for down-regulated
              message("   processEnrichments Step: CALLING perform_enrichment for DOWN-REGULATED proteins...")
              message(sprintf("      About to analyze %d down-regulated proteins", nrow(down_matrix)))
              
              protein_col <- protein_id_column
              message(sprintf("      Using protein_id_column: %s", protein_col))
              message(sprintf("      Using species: %s", species))
              message(sprintf("      Using threshold: %s", q_cutoff))
              message(sprintf("      Using exclude_iea: %s", exclude_iea))
              
              down_result <- perform_enrichment(
                data_subset = down_matrix,
                species = species,
                threshold = q_cutoff,
                sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                domain_scope = "custom",
                custom_bg = custom_bg,
                exclude_iea = exclude_iea,
                protein_id_column = protein_col,
                correction_method = correction_method
              )
              
              # ✅ DEBUG 66: Log result from perform_enrichment
              message("   processEnrichments Step: DOWN-REGULATED enrichment completed")
              message("      Data State (down_result): Checking perform_enrichment result...")
              
              if (is.null(down_result)) {
                message("      Data State (down_result): Result is NULL")
              } else {
                message("      Data State (down_result) Structure:")
                utils::str(down_result)
                
                if ("result" %in% names(down_result) && !is.null(down_result$result)) {
                  message(sprintf("      Data State (down_result$result): Found %d enriched terms", nrow(down_result$result)))
                  if (nrow(down_result$result) > 0) {
                    message("      Data State (down_result$result) Head:")
                    print(head(down_result$result))
                  }
                } else {
                  message("      Data State (down_result): No 'result' component or result is NULL")
                }
              }
              
              down_result
            } else {
              message("   processEnrichments Step: SKIPPING DOWN-REGULATED enrichment (no proteins)")
              NULL
            }
          }, error = function(e) {
            message(sprintf("   processEnrichments Step: ERROR in down-regulated enrichment: %s", e$message))
            warning(sprintf("Error processing down-regulated genes: %s", e$message))
            message("Debug info for down-regulated genes:")
            message("Number of proteins: ", nrow(down_matrix))
            NULL
          })
        )
        }, error = function(e) {
          message(sprintf("*** ERROR in contrast processing: %s", e$message))
          message(sprintf("*** ERROR details: %s", class(e)))
          message(sprintf("*** Call stack: %s", paste(deparse(sys.calls()), collapse = " -> ")))
          return(NULL)
        })
      })

    # Store enrichment results
    enrichment_results@enrichment_data <- results

    # Generate and store both static and interactive plots
    plot_results <- purrr::map(names(results), function(contrast) {
      list(
        up = if(!is.null(results[[contrast]]$up)) {
          generate_enrichment_plots(results[[contrast]]$up, contrast, "up", pathway_dir)
        } else NULL,

        down = if(!is.null(results[[contrast]]$down)) {
          generate_enrichment_plots(results[[contrast]]$down, contrast, "down", pathway_dir)
        } else NULL
      )
    }) |>
      purrr::set_names(names(results))

    # Store static plots
    enrichment_results@enrichment_plots <- purrr::map(plot_results, function(x) {
      list(
        up = if(!is.null(x$up)) x$up$static else NULL,
        down = if(!is.null(x$down)) x$down$static else NULL
      )
    })

    # Store interactive plotly objects
    enrichment_results@enrichment_plotly <- purrr::map(plot_results, function(x) {
      list(
        up = if(!is.null(x$up)) x$up$interactive else NULL,
        down = if(!is.null(x$down)) x$down$interactive else NULL
      )
    })

    # Generate and store summaries
    enrichment_results@enrichment_summaries <- purrr::map(names(results), function(contrast) {
      list(
        up = if(!is.null(results[[contrast]]$up)) summarize_enrichment(results[[contrast]]$up) else NULL,
        down = if(!is.null(results[[contrast]]$down)) summarize_enrichment(results[[contrast]]$down) else NULL
      )
    }) |>
      purrr::set_names(names(results))

    # Explicitly ensure the names of plot_results match contrast_names_to_use
    # This guards against issues if names were lost or mismatched during creation
    if (!identical(names(plot_results), contrast_names)) {
       message("DEBUG: Mismatch detected or names missing in plot_results. Reassigning names.")
       # Check if the number of elements still matches before trying to assign names
       if(length(plot_results) == length(contrast_names)) {
         names(plot_results) <- contrast_names
       } else {
         stop("Critical error: Number of plot results does not match the number of contrast names.")
       }
    }

    # Save plots using the desired contrast names for files
    purrr::iwalk(contrast_names, function(contrast, i) {
      plots <- plot_results[[contrast]]

      # Double-check if plots object is NULL, which might happen if naming failed
      if(is.null(plots)) {
          warning(paste("Could not find plot data for contrast:", contrast, "- Skipping save for this contrast."))
          return(NULL) # Skip to the next iteration
      }

      # Simple sanitization (should be redundant now)
      sanitized_contrast <- gsub("[^A-Za-z0-9_.-]", "_", contrast)

      message(paste("Loop", i, "- Saving plots for contrast:", contrast, " (Sanitized:", sanitized_contrast, ")")) # Debug

      purrr::walk(c("up", "down"), function(direction) {
        # Check if the plot component exists and is not NULL
        if (!is.null(plots[[direction]]) && !is.null(plots[[direction]]$interactive)) {

          file_name <- paste0(sanitized_contrast, "_", direction, "_enrichment_plot.html")
          file_path <- file.path(pathway_dir, file_name)
          
          # Add PNG file path
          png_file_name <- paste0(sanitized_contrast, "_", direction, "_enrichment_plot.png")
          png_file_path <- file.path(pathway_dir, png_file_name)

          message(paste("  Attempting to save:", file_path)) # Debug

          # Ensure the directory exists before saving
          target_dir_for_file <- dirname(file_path)
          if (!dir.exists(target_dir_for_file)) {
             message(paste("  Creating directory:", target_dir_for_file)) # Debug
             dir.create(target_dir_for_file, recursive = TRUE)
          }

          tryCatch({
            # Save the interactive Plotly HTML file
            htmlwidgets::saveWidget(
              plots[[direction]]$interactive,
              file_path,
              selfcontained = TRUE
            )
            message(paste("  Successfully saved:", file_path)) # Debug success
            
            # Save the static ggplot as PNG
            ggplot2::ggsave(
              filename = png_file_path,
              plot = plots[[direction]]$static,
              width = 10, 
              height = 8,
              dpi = 300,
              bg = "white"
            )
            message(paste("  Successfully saved:", png_file_path)) # Debug success
          }, error = function(e) {
            # Print more detailed error context
            warning(paste("  ERROR saving plots for contrast:", contrast,
                          "direction:", direction,
                          "path:", file_path,
                          "- Error message:", e$message))
          })
        } else {
           message(paste("  Skipping save for direction:", direction, "- Plot component is NULL or not interactive."))
        }
      })
    })

    return(enrichment_results)

  } else {
    if(is.null(go_annotations)) {
      stop("For unsupported organisms, GO annotations must be provided")
    }

    message(sprintf("Using custom GO annotations for taxon ID %s", taxon_id))

    # Prepare GO term mappings
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

    # Combine all terms
    all_terms <- rbind(
      cbind(bp_terms, ONTOLOGY = "BP"),
      cbind(mf_terms, ONTOLOGY = "MF"),
      cbind(cc_terms, ONTOLOGY = "CC")
    )

    # Create term mappings with explicit dplyr namespace
    term2gene <- all_terms |>
      dplyr::select(.data$TERM, .data$Entry) |>
      dplyr::distinct()

    term2name <- data.frame(
      TERM = unique(all_terms$TERM),
      NAME = purrr::map_chr(unique(all_terms$TERM),
                            ~tryCatch(GO.db::Term(GO.db::GOTERM[[.x]]),
                                      error = function(e) .x))
    )

    # Get the internal long names (only used initially if short names aren't provided or needed for mapping)
    internal_contrast_names <- names(de_results@de_data)

    # Determine which names to use (Prefer explicitly passed short names)
    if (is.null(contrast_names)) {
      warning("Explicit contrast_names not provided, using internal names from de_results@de_data which might be long or contain invalid characters.")
      contrast_names_to_use <- internal_contrast_names
    } else {
      if(length(contrast_names) != length(internal_contrast_names)) {
        stop("Length of provided 'contrast_names' does not match the number of contrasts in 'de_results@de_data'.")
      }
      contrast_names_to_use <- contrast_names # Use the short names
    }
    message("DEBUG: Using the following contrast names for processing and output:")
    print(contrast_names_to_use)

    # Initialize lists
    results_list_long_names <- list() # Temporary list to hold results with original names
    enrichment_results <- createEnrichmentResults(de_results@contrasts)

    # --- Step 1: Process enrichment using internal loop mapped to short names ---
    message("--- Starting Enrichment Processing Loop ---")
    results_list_long_names <- purrr::map2(contrast_names_to_use, internal_contrast_names, function(short_name, long_name) {
        message(paste("Processing contrast:", short_name, "(original:", long_name, ")"))

        de_data <- de_results@de_data[[long_name]] # Access input data using long name

        if(is.null(de_data)) {
          warning(paste("No DE data found for internal contrast:", long_name))
          # Return NULL for both up and down to keep list structure aligned
          return(list(up = NULL, down = NULL))
        }

        # Split data into up/down regulated
        message(sprintf("Total proteins before filtering: %d", nrow(de_data)))

        subset_sig <- de_data |>
          dplyr::filter(.data$fdr_qvalue < q_cutoff)

        message(sprintf("Proteins passing FDR cutoff (%g): %d", q_cutoff, nrow(subset_sig)))

        up_genes <- subset_sig |>
          dplyr::filter(.data$log2FC > up_cutoff) |>
          dplyr::pull({{protein_id_column}})

        down_genes <- subset_sig |>
          dplyr::filter(.data$log2FC < -down_cutoff) |>
          dplyr::pull({{protein_id_column}})

        # Background genes
        background_IDs <- unique(de_data |> dplyr::pull({{protein_id_column}}))

        # Perform enrichment for up-regulated genes
        up_enrich <- tryCatch({
          if(length(up_genes) > 0) {
            clusterProfiler::enricher(
              gene = up_genes,
              universe = background_IDs,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = q_cutoff,
              pAdjustMethod = "BH"
            )
          } else NULL
        }, error = function(e) {
          warning(sprintf("Error processing up-regulated genes: %s", e$message))
          NULL
        })

        # Perform enrichment for down-regulated genes
        down_enrich <- tryCatch({
          if(length(down_genes) > 0) {
            clusterProfiler::enricher(
              gene = down_genes,
              universe = background_IDs,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = q_cutoff,
              pAdjustMethod = "BH"
            )
          } else NULL
        }, error = function(e) {
          warning(sprintf("Error processing down-regulated genes: %s", e$message))
          NULL
        })

        # Return results for this contrast
        list(up = up_enrich, down = down_enrich)
    }) |> purrr::set_names(internal_contrast_names)
    message("--- Finished Enrichment Processing Loop ---")

    # --- Step 2: Assign results to the S4 object with SHORT names ---
    if(length(results_list_long_names) == length(contrast_names_to_use)) {
        names(results_list_long_names) <- contrast_names_to_use # Rename the temporary list
        enrichment_results@enrichment_data <- results_list_long_names # Assign renamed list
        message("DEBUG: Assigned enrichment data to S4 object with short names.")
    } else {
        stop("Mismatch between number of processed results and expected contrast names.")
    }
    # Now enrichment_results@enrichment_data uses SHORT names

    # Create GO term mappings once (moved outside the plotting function)
    go_term_map <- dplyr::bind_rows(
      go_annotations |>
        tidyr::separate_rows(.data$go_id_go_biological_process, .data$go_term_go_biological_process, sep = "; ") |>
        dplyr::select(.data$go_id_go_biological_process, .data$go_term_go_biological_process) |>
        dplyr::rename(ID = .data$go_id_go_biological_process, term = .data$go_term_go_biological_process),

      go_annotations |>
        tidyr::separate_rows(.data$go_id_go_molecular_function, .data$go_term_go_molecular_function, sep = "; ") |>
        dplyr::select(.data$go_id_go_molecular_function, .data$go_term_go_molecular_function) |>
        dplyr::rename(ID = .data$go_id_go_molecular_function, term = .data$go_term_go_molecular_function),

      go_annotations |>
        tidyr::separate_rows(.data$go_id_go_cellular_compartment, .data$go_term_go_cellular_compartment, sep = "; ") |>
        dplyr::select(.data$go_id_go_cellular_compartment, .data$go_term_go_cellular_compartment) |>
        dplyr::rename(ID = .data$go_id_go_cellular_compartment, term = .data$go_term_go_cellular_compartment)
    ) |>
      dplyr::distinct()

    # Create category mapping once
    go_category_map <- all_terms |>
      dplyr::distinct(.data$TERM, .data$ONTOLOGY) |>
      dplyr::mutate(
        source = dplyr::case_when(
          .data$ONTOLOGY == "BP" ~ "GO:BP",
          .data$ONTOLOGY == "CC" ~ "GO:CC",
          .data$ONTOLOGY == "MF" ~ "GO:MF"
        )
      )

    # --- Step 3: Generate Plots using SHORT names ---
    message("--- Starting Plot Generation Loop ---")
    plot_results_list <- purrr::map(contrast_names_to_use, function(contrast) {
        message(paste("Generating plots for contrast:", contrast))
        # Access enrichment data using SHORT name
        current_enrich_data <- enrichment_results@enrichment_data[[contrast]]

        # Process up and down directions
        direction_results <- purrr::map(c("up", "down"), function(direction) {
            result_data <- current_enrich_data[[direction]] # Access up/down list

            if(!is.null(result_data) && nrow(result_data@result) > 0) {
                message(sprintf("Processing %s-regulated genes for contrast %s", direction, contrast))

                # Prepare data for plotting and tables
                plot_data <- result_data@result |>
                  dplyr::left_join(go_category_map, by = c("ID" = "TERM")) |>
                  dplyr::left_join(go_term_map, by = "ID") |>
                  dplyr::mutate(
                    source = dplyr::coalesce(.data$source, "Other"),
                    source = factor(.data$source, levels = c("GO:BP", "GO:CC", "GO:MF", "Other")),
                    neg_log10_q = -log10(.data$qvalue),  # Using qvalue directly from clusterProfiler output
                    gene_count = .data$Count,
                    significant = .data$qvalue < q_cutoff  # Add significance flag based on q_cutoff
                  ) |>
                  dplyr::mutate(
                    term = dplyr::coalesce(.data$term, .data$Description)
                  )

                # Save results table
                readr::write_tsv(
                  plot_data,
                  file.path(pathway_dir,
                            paste0(contrast, "_", direction, "_enrichment_results.tsv"))
                )

                # Generate static plot with q-value threshold line
                # First identify top significant terms to label
                top_terms <- plot_data %>%
                  dplyr::filter(.data$qvalue < q_cutoff) %>%
                  dplyr::arrange(.data$qvalue) %>%
                  dplyr::slice_head(n = 5) # Label top 5 most significant terms
                
                static <- ggplot2::ggplot(plot_data,
                                          ggplot2::aes(x = .data$source,
                                                       y = .data$neg_log10_q,
                                                       text = paste0(
                                                         "Term: ", .data$term, "\n",
                                                         "ID: ", .data$ID, "\n",
                                                         "Genes: ", .data$Count, "\n",
                                                         "Gene Ratio: ", .data$GeneRatio, "\n",
                                                         "Background Ratio: ", .data$BgRatio, "\n",
                                                         "Q-value: ", signif(.data$qvalue, 3)
                                                       ))) +
                  ggplot2::geom_hline(yintercept = -log10(q_cutoff),
                                      linetype = "dashed",
                                      color = "darkgrey") +
                  ggplot2::geom_jitter(ggplot2::aes(size = .data$gene_count,
                                                    color = -log10(.data$qvalue)),
                                       alpha = 0.7,
                                       width = 0.2) +
                  # Add labels for significant terms
                  ggrepel::geom_text_repel(
                    data = top_terms,
                    ggplot2::aes(label = .data$term),
                    size = 3,
                    max.overlaps = 15,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    force = 5
                  ) +
                  ggplot2::scale_color_gradient(low = "#FED976",
                                                high = "#800026",
                                                name = "-log10(q-value)") +
                  ggplot2::scale_size_continuous(name = "Gene Count",
                                                 range = c(3, 12)) +
                  ggplot2::theme_minimal() +
                  ggplot2::theme(
                    axis.text.x = ggplot2::element_text(size = 10, angle = 0),
                    axis.text.y = ggplot2::element_text(size = 8),
                    plot.title = ggplot2::element_text(size = 12, face = "bold"),
                    legend.title = ggplot2::element_text(size = 10),
                    legend.text = ggplot2::element_text(size = 8)
                  ) +
                  ggplot2::labs(
                    title = paste0(contrast, " ", tools::toTitleCase(direction), "-regulated"),
                    x = "GO Category",
                    y = "-log10(q-value)"
                  )

                list(
                  static = static,
                  interactive = plotly::ggplotly(static, tooltip = "text")
                )
            } else {
                NULL # Return NULL if no results
            }
        }) |> purrr::set_names(c("up", "down"))
        
        # Return plot components for this contrast
        direction_results
    }) |> purrr::set_names(contrast_names_to_use)
    message("--- Finished Plot Generation Loop ---")

    # --- Step 4: Store plots in S4 object using SHORT names ---
    enrichment_results@enrichment_plots <- purrr::map(plot_results_list, function(x) {
      list(up = x$up$static, down = x$down$static) # Handle NULLs gracefully if needed
    })
    enrichment_results@enrichment_plotly <- purrr::map(plot_results_list, function(x) {
      list(up = x$up$interactive, down = x$down$interactive) # Handle NULLs
    })
    message("DEBUG: Assigned plot data to S4 object slots.")
    message("DEBUG: Names in enrichment_results@enrichment_plotly:")
    print(names(enrichment_results@enrichment_plotly))

    # --- Step 5: Save interactive plots using SHORT names ---
    message("--- Starting Final Save Loop ---")
    purrr::walk(contrast_names_to_use, function(contrast) {
        # Access plot data using SHORT name (now stored with short names)
        plots <- enrichment_results@enrichment_plotly[[contrast]]

        if(is.null(plots)) {
            # This warning should NOT trigger if steps above worked
            warning(paste("Could not find plot data for contrast:", contrast, "- Skipping save. (This indicates an earlier issue)"))
            return(NULL)
        }

        # Sanitize (optional if short names are clean)
        sanitized_contrast <- gsub("[^A-Za-z0-9_.-]", "_", contrast)
        message(paste("Saving plots for contrast:", contrast))

        # Process both directions
        purrr::walk(c("up", "down"), function(direction) {
            # Check if interactive plot exists
            if (!is.null(plots[[direction]])) { # Access up/down list
                interactive_plot <- plots[[direction]] # The actual plotly object
                static_plot <- enrichment_results@enrichment_plots[[contrast]][[direction]] # The static ggplot

                file_name <- paste0(sanitized_contrast, "_", direction, "_enrichment_plot.html")
                file_path <- file.path(pathway_dir, file_name)
                
                # Add PNG file path
                png_file_name <- paste0(sanitized_contrast, "_", direction, "_enrichment_plot.png")
                png_file_path <- file.path(pathway_dir, png_file_name)
                
                message(paste("  Attempting to save:", file_path))

                # Ensure directory exists
                target_dir_for_file <- dirname(file_path)
                if (!dir.exists(target_dir_for_file)) {
                   dir.create(target_dir_for_file, recursive = TRUE)
                }

                # Save the widget
                tryCatch({
                  # Save HTML interactive plot
                  htmlwidgets::saveWidget(
                    interactive_plot, # Pass the plot object directly
                    file_path,
                    selfcontained = TRUE
                  )
                  message(paste("  Successfully saved:", file_path))
                  
                  # Save static ggplot as PNG
                  if (!is.null(static_plot)) {
                    ggplot2::ggsave(
                      filename = png_file_path,
                      plot = static_plot,
                      width = 10, 
                      height = 8,
                      dpi = 300,
                      bg = "white"
                    )
                    message(paste("  Successfully saved:", png_file_path))
                  }
                }, error = function(e) {
                  warning(paste("  ERROR saving plots:", file_path, "-", e$message))
                })
            } else {
                message(paste("  Skipping save for direction:", direction, "- Plot component is NULL."))
            }
        })
    })
    message("--- Finished Final Save Loop ---")

    return(enrichment_results)
  }
}


# ----------------------------------------------------------------------------
# getEnrichmentResult
# ----------------------------------------------------------------------------
# Helper function to access results
#' @export
getEnrichmentResult <- function(enrichment_results, contrast, direction) {
  if(!contrast %in% names(enrichment_results@enrichment_data)) {
    stop("Contrast not found")
  }
  if(!direction %in% c("up", "down")) {
    stop("Direction must be 'up' or 'down'")
  }
  enrichment_results@enrichment_data[[contrast]][[direction]]
}


# ----------------------------------------------------------------------------
# getEnrichmentPlotly
# ----------------------------------------------------------------------------
# Helper function to access plotly objects
#' @export
getEnrichmentPlotly <- function(enrichment_results, contrast, direction) {
  if(!contrast %in% names(enrichment_results@enrichment_plotly)) {
    stop("Contrast not found")
  }
  if(!direction %in% c("up", "down")) {
    stop("Direction must be 'up' or 'down'")
  }
  enrichment_results@enrichment_plotly[[contrast]][[direction]]
}


# ----------------------------------------------------------------------------
# getEnrichmentSummary
# ----------------------------------------------------------------------------
# Helper function to get summary
#' @export
getEnrichmentSummary <- function(enrichment_results) {
  summaries <- purrr::map_df(names(enrichment_results@enrichment_summaries), function(contrast) {
    summary <- enrichment_results@enrichment_summaries[[contrast]]

    data.frame(
      contrast = contrast,
      up_total = if(!is.null(summary$up)) summary$up$total else 0,
      down_total = if(!is.null(summary$down)) summary$down$total else 0,
      up_GO_BP = if(!is.null(summary$up)) summary$up$GO_BP else 0,
      down_GO_BP = if(!is.null(summary$down)) summary$down$GO_BP else 0,
      up_GO_CC = if(!is.null(summary$up)) summary$up$GO_CC else 0,
      down_GO_CC = if(!is.null(summary$down)) summary$down$GO_CC else 0,
      up_GO_MF = if(!is.null(summary$up)) summary$up$GO_MF else 0,
      down_GO_MF = if(!is.null(summary$down)) summary$down$GO_MF else 0,
      up_KEGG = if(!is.null(summary$up)) summary$up$KEGG else 0,
      down_KEGG = if(!is.null(summary$down)) summary$down$KEGG else 0,
      up_REAC = if(!is.null(summary$up)) summary$up$REAC else 0,
      down_REAC = if(!is.null(summary$down)) summary$down$REAC else 0
    )
  })

  return(summaries)
}

