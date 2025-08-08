#' @title Parse a Delimited String of Numbers
#' @description This utility function parses a character string containing numbers
#' separated by common delimiters (`,`, `.`, `;`, `:`) into an integer vector.
#' If no delimiters are found, it assumes the string is a single number.
#'
#' @param input_text A character string containing one or more numbers.
#'
#' @return An integer vector of the parsed numbers.
#'
#' @importFrom stringr str_detect str_split
#' @importFrom purrr map_int
#' @export
#'
#' @examples
#' parseNumList("1,2,3")
#' parseNumList("10;20;30")
#' parseNumList("5")
parseNumList <-  function ( input_text ) {
  if( str_detect( input_text, "[.,;:]")) {
    str_split( input_text, "[.,;:]")[[1]] %>% purrr::map_int( as.integer)
  } else {
    return( as.integer( input_text ))
  }
}

#' @title Convert an ID to an Annotation using a Dictionary
#' @description A simple key-value lookup function that retrieves an annotation for a
#' given ID from a named list or vector (acting as a dictionary).
#'
#' @param id The identifier (key) to look up.
#' @param id_to_annotation_dictionary A named list or vector where names are the
#'   IDs and values are the annotations.
#'
#' @return The corresponding annotation value if the ID is found, otherwise `NA_character_`.
#' @export
#'
#' @examples
#' my_dict <- c("ID1" = "Annotation A", "ID2" = "Annotation B")
#' convertIdToAnnotation("ID1", my_dict)
#' convertIdToAnnotation("ID3", my_dict)
convertIdToAnnotation <- function( id, id_to_annotation_dictionary) {

    return( ifelse( !is.null(id_to_annotation_dictionary[[id]] ),
                    id_to_annotation_dictionary[[id]] ,
                    NA_character_))

}


#' @title Perform Gene Ontology (GO) Enrichment for a Single Aspect
#' @description This function runs a GO enrichment analysis for a specific GO aspect
#' (e.g., Biological Process) using `clusterProfiler::enricher`. It filters the
#' GO terms by size and ensures that only terms with at least two significant
#' genes in the query set are considered.
#'
#' @details
#' The workflow is as follows:
#' 1.  Filters the main GO annotation table for the specified `go_aspect`.
#' 2.  Filters the GO terms to include only those that fall within the specified
#'     `min_gene_set_size` and `max_gene_set_size`, based on the background list.
#' 3.  Removes "singleton" GO terms, which are terms that are associated with only
#'     one gene in the user's `query_list`. This avoids spurious, single-gene enrichments.
#' 4.  Performs the enrichment test using `clusterProfiler::enricher`.
#' 5.  Formats the results into a data frame and adds back the term descriptions.
#'
#' @param go_annot A data frame of GO annotations, containing columns for protein IDs,
#'   GO IDs, GO terms, and GO aspects.
#' @param background_list A data frame or vector of all background protein IDs used in the experiment.
#' @param go_aspect A character string specifying the GO aspect to test (e.g., "BP", "CC", "MF").
#'   If `NA`, no aspect filtering is performed.
#' @param query_list A character vector of the protein IDs to be tested for enrichment (the "query set").
#' @param id_to_annotation_dictionary A named list/vector for mapping GO IDs to GO term names.
#' @param annotation_id The unquoted column name for GO identifiers in `go_annot`.
#' @param protein_id The unquoted column name for protein identifiers in `go_annot`.
#' @param aspect_column The unquoted column name for the GO aspect in `go_annot`.
#' @param p_val_thresh The p-value cutoff for significance.
#' @param min_gene_set_size The minimum number of genes a GO term must have to be included in the analysis.
#' @param max_gene_set_size The maximum number of genes a GO term can have to be included.
#' @param get_cluster_profiler_object Logical. If `TRUE`, returns a list containing both
#'   the results data frame and the raw `enrichResult` object from clusterProfiler.
#'   Defaults to `FALSE`.
#'
#' @return If `get_cluster_profiler_object` is `FALSE` (default), a data frame of
#'   the enrichment results. If `TRUE`, a list containing the results data frame
#'   and the `enrichResult` object. Returns `NA` if no enrichment results are found.
#'
#' @importFrom rlang set_names as_name enquo
#' @importFrom clusterProfiler enricher
#' @export
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



#' @title Run GO Enrichment Analysis
#' @description A wrapper function to perform Gene Ontology (GO) enrichment analysis
#' for a set of significant proteins against a background set. It runs the analysis
#' for all three main GO aspects (Biological Process, Cellular Compartment, Molecular Function).
#'
#' @details
#' The function filters GO terms based on size and the number of significant proteins
#' they contain. It then uses `clusterProfiler::enricher` to perform a hypergeometric
#' test for each GO aspect. The results are combined and annotated with gene symbols.
#'
#' @param significant_proteins A character vector of significant protein identifiers (e.g., UniProt accessions).
#' @param background_proteins A character vector of all protein identifiers used as the background/universe.
#' @param go_annotations A data frame containing GO annotations, mapping protein IDs to GO IDs, terms, and types.
#' @param uniprot_data A data frame containing UniProt data, used to map protein accessions to gene names.
#' @param p_val_thresh The p-value cutoff for enrichment significance. Defaults to 0.05.
#' @param min_gene_set_size The minimum number of genes a GO term must have in the background to be tested. Defaults to 4.
#' @param max_gene_set_size The maximum number of genes a GO term can have. Defaults to 200.
#' @param min_sig_gene_set_size The minimum number of significant genes a GO term must have to be included in the results. Defaults to 2.
#'
#' @return A data frame of enrichment results for all three GO domains, or `NULL` if no enrichment is found.
#'
#' @importFrom clusterProfiler enricher
#' @importFrom dplyr inner_join group_by summarise ungroup filter select distinct mutate left_join pull case_when
#' @importFrom purrr map map_chr discard reduce
#' @importFrom stringr str_split
#' @export
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

#' @title Convert Protein Accessions to a String of Gene Symbols
#' @description Looks up a list of protein accessions in a dictionary and returns a
#' single string of the corresponding gene symbols, concatenated with a forward slash.
#'
#' @param gene_id_list A character vector of protein accessions to convert.
#' @param dictionary A named vector or list where names are protein accessions and
#'   values are gene symbols.
#'
#' @return A single character string with gene symbols separated by "/".
#' @export
convertProteinAccToGeneSymbol <- function( gene_id_list, dictionary ) {

  purrr::map_chr( gene_id_list,
                  ~{ ifelse( . %in% names(dictionary ),
                             dictionary[[.]],
                             NA_character_)   } )  %>%
    paste( collapse="/")
}


#' @title Build a Dictionary Mapping Annotation IDs to Names
#' @description Creates a named vector to serve as a lookup dictionary, mapping
#' unique annotation identifiers (e.g., GO IDs) to their descriptive names (e.g., GO terms).
#'
#' @param input_table A data frame containing at least two columns: one for annotation IDs
#'   and one for annotation names.
#' @param annotation_column The unquoted column name for the annotation descriptions/names.
#' @param annotation_id_column The unquoted column name for the annotation identifiers.
#'
#' @return A named character vector where names are the annotation IDs and values are the annotation names.
#' @export
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


#' @title Create a Gene Set List from an Annotation Table
#' @description Transforms a data frame of protein-to-annotation mappings into a
#' named list, where each element represents a gene set. This format is required
#' by many enrichment analysis tools.
#'
#' @details The function groups the input table by the annotation ID. For each
#' annotation ID, it collects all associated protein IDs into a vector.
#'
#' @param input_table A data frame containing protein-to-annotation mappings.
#' @param annotation_id The unquoted column name for the annotation identifier (e.g., pathway ID).
#' @param protein_id The unquoted column name for the protein identifier.
#'
#' @return A named list where names are the annotation IDs and each element is a
#'   character vector of associated protein IDs.
#' @export
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

#' @title Convert a Table to a List of Tables by a Column
#' @description Splits a data frame into a list of data frames based on the unique
#' values in a specified column.
#'
#' @param input_table The data frame to be split.
#' @param column_name The unquoted name of the column to group by. The unique values
#'   in this column will become the names of the list elements.
#'
#' @return A named list of data frames.
#' @export
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



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Run Gene Set Enrichment Analysis (GSEA)
#' @description This function performs Gene Set Enrichment Analysis (GSEA) on a
#' ranked list of genes against a specified collection of gene sets (e.g., from MSigDB).
#' It filters gene sets by size before running the analysis.
#'
#' @param index_name The name of the gene set collection to use from `list_of_gene_sets`.
#' @param contrast_name The name of the contrast (comparison) to use from `list_of_de_proteins`,
#'   which contains the ranked gene lists.
#' @param list_of_de_proteins A named list where each element is a ranked list of genes
#'   (numeric vector with gene names as names) for a specific contrast.
#' @param list_of_gene_sets A named list where each element is a `GSEABase` gene set
#'   collection object.
#' @param min_set_size The minimum number of genes a gene set must have to be included.
#' @param max_set_size The maximum number of genes a gene set can have to be included.
#'
#' @return A `gseaResult` object from `clusterProfiler::GSEA`.
#'
#' @importFrom GSEABase geneIds
#' @importFrom tibble tibble
#' @importFrom tidyr unnest
#' @importFrom dplyr inner_join group_by summarise ungroup filter select mutate
#' @importFrom clusterProfiler GSEA
#' @export
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
#' @title Run Over-Representation Analysis (ORA) with clusterProfiler
#' @description This function performs an Over-Representation Analysis (ORA) using
#' `clusterProfiler::enricher` on a list of significant genes against a specified
#' collection of gene sets. It filters gene sets by size.
#'
#' @param index_name The name of the gene set collection to use from `list_of_gene_sets`.
#' @param contrast_name The name of the contrast to use from `list_of_de_proteins`,
#'   which contains the list of significant genes.
#' @param list_of_de_proteins A named list where each element is a character vector
#'   of significant gene identifiers for a specific contrast.
#' @param list_of_gene_sets A named list where each element is a `GSEABase` gene set
#'   collection object.
#' @param min_set_size The minimum number of genes a gene set must have to be included.
#' @param max_set_size The maximum number of genes a gene set can have to be included.
#'
#' @return An `enrichResult` object from `clusterProfiler::enricher`.
#'
#' @importFrom GSEABase geneIds
#' @importFrom tibble tibble
#' @importFrom tidyr unnest
#' @importFrom dplyr inner_join group_by summarise ungroup filter select mutate
#' @importFrom clusterProfiler enricher
#' @export
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

#' @title Create a UniProt Accession to Gene Symbol Dictionary
#' @description This function processes a table of UniProt annotations to create a
#' named vector that serves as a lookup dictionary from UniProt accessions to gene symbols.
#' It handles cases where multiple gene symbols are present by taking the first one.
#'
#' @param input_table A data frame containing UniProt data.
#' @param protein_id_lookup_column The unquoted column name for UniProt accessions in `input_table`.
#' @param gene_symbol_column The unquoted column name for gene symbols in `input_table`.
#' @param protein_id The desired unquoted name for the protein ID column in the output dictionary (will be the names of the vector).
#'
#' @return A named character vector where names are UniProt accessions and values are gene symbols.
#'
#' @importFrom dplyr select rename mutate distinct pull
#' @importFrom stringr str_split
#' @importFrom purrr map_chr
#' @importFrom rlang enquo as_name
#' @export
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

#######################################################
#' @title Query the Revigo Web Service
#' @description Submits a list of Gene Ontology (GO) terms to the Revigo web service
#' to summarize and reduce redundancy. It parses the resulting HTML table into a data frame.
#'
#' @param input_list A character vector of GO term IDs.
#' @param cutoff A numeric value between 0.4 and 0.9 representing the similarity
#'   cutoff for Revigo. Allowed values: 0.9, 0.7, 0.5, 0.4. Defaults to 0.5.
#' @param speciesTaxon The NCBI taxon ID for the species (e.g., 9606 for Human, 10090 for Mouse).
#'   Defaults to 10090 (Mouse).
#' @param temp_file An optional file path to save the intermediate HTML response from Revigo.
#'   If `NA` or `NULL`, a temporary file is created and deleted automatically.
#'
#' @return A data frame (tibble) containing the summarized GO terms from Revigo,
#'   including columns like `Term ID`, `Name`, `Frequency`, `Plot X`, `Plot Y`,
#'   `Plot Size`, `Uniqueness`, `Dispensability`, and `Eliminated`.
#'
#' @importFrom httr POST content
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom rvest read_html html_nodes html_table
#' @importFrom purrr map discard
#' @importFrom dplyr bind_rows
#' @export
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



#' @title Cluster Enriched Pathways
#' @description This function takes a table of enrichment results, resolves duplicate
#' entries for the same pathway, calculates a signed score based on p-values and
#' regulation direction, and performs hierarchical clustering on the pathways.
#'
#' @details
#' The scoring logic is as follows:
#' - `score = -log10(p.adjust)` for up-regulated or neutral gene sets.
#' - `score = -1 * -log10(p.adjust)` for down-regulated gene sets.
#'
#' Duplicate entries (same pathway enriched in both up- and down-regulated sets for the same comparison)
#' can be handled by either removing them (`remove_duplicted_entries = "delete"`) or merging them
#' into a "shared" category (`remove_duplicted_entries = "merge"`).
#'
#' The function returns an ordered data frame based on the clustering, which can be
#' used for creating ordered heatmaps.
#'
#' @param input_table A data frame of enrichment results. Must contain columns for
#'   `comparison`, `annotation_id`, `p.adjust`, `gene_set`, and `term`.
#' @param added_columns A character vector of additional column names to include in grouping
#'   when identifying and clustering pathways.
#' @param remove_duplicted_entries A string indicating how to handle duplicate entries.
#'   Can be `"delete"` (or `TRUE`), `"merge"`. Defaults to `TRUE`.
#'
#' @return A data frame of enrichment results, ordered according to the hierarchical
#'   clustering of pathways. An `ordering` column is added. Returns a slightly
#'   different, un-clustered table if there are fewer than 2 rows to cluster.
#'
#' @importFrom dplyr mutate case_when group_by summarise ungroup filter anti_join inner_join bind_rows select left_join arrange row_number
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom stats hclust dist cutree
#' @export
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

########################

#' @title Generate an Enrichment Heatmap
#' @description Creates a ggplot2 heatmap of functional enrichment results. The heatmap
#' displays enriched terms on the y-axis and experimental comparisons or other
#' categories on the x-axis. The size of the points corresponds to the statistical
#' significance, and the shape/color can represent the direction of regulation.
#'
#' @param input_table A data frame of enrichment results, typically the output of
#'   `clusterPathways` to ensure correct ordering. It must contain columns for the
#'   term, the x-axis variable, `neg_log_p_value`, and `gene_set`.
#' @param x_axis The unquoted column name to be used for the x-axis of the heatmap.
#' @param input_go_type An optional character string to filter the results to a single
#'   GO type (e.g., "Biological Process"). If `NA`, all types are used.
#' @param input_plot_title The title for the plot.
#' @param facet_by_column Optional unquoted column name to facet the plot by.
#' @param xaxis_levels An optional character vector to specify the order of items
#'   on the x-axis.
#' @param scales The `scales` argument passed to `facet_wrap` (e.g., "free", "fixed").
#'   Defaults to "fixed".
#'
#' @return A ggplot object representing the enrichment heatmap.
#'
#' @import ggplot2
#' @importFrom dplyr filter mutate pull
#' @importFrom purrr map_dbl map_chr
#' @importFrom stringr str_wrap
#' @importFrom rlang enquo quo_get_expr as_name
#' @export
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


#' @title Read and Combine Enrichment Result Files
#' @description Reads multiple enrichment result files (in a tabular format like .tab or .csv)
#' specified in a manifest data frame, combines them into a single data frame, and
#' adds metadata columns from the manifest.
#'
#' @param table_of_files A data frame where one column contains the file paths to the
#'   enrichment results, and other columns contain metadata to be joined to the results.
#' @param file_names_column The unquoted column name in `table_of_files` that holds the
#'   file paths. Defaults to `file_name`.
#' @param go_type A character string to assign as the `go_type` if this column is not
#'   present in the result files. Defaults to "KEGG".
#'
#' @return A single data frame containing the combined and annotated enrichment results.
#'
#' @importFrom vroom vroom
#' @importFrom purrr map map_int keep
#' @importFrom dplyr pull bind_rows rename left_join relocate select mutate
#' @importFrom rlang enquo as_name
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

# enriched_results_tbl <- readEnrichmentResultFiles( table_of_files, go_type="KEGG")

#' @title Filter Enrichment Results Using Revigo
#' @description This function takes a table of enrichment results, groups them by
#' specified columns, sends the GO terms for each group to Revigo for summarization,
#' and then filters the original table to keep only the non-redundant terms
#' returned by Revigo.
#'
#' @param enriched_results_tbl A data frame of enrichment results.
#' @param added_columns A character vector of column names to group by before sending
#'   data to Revigo. `comparison`, `gene_set`, and `go_type` are also used for grouping.
#' @param is_run_revigo A logical flag. If `TRUE` (default), the Revigo filtering is
#'   performed. If `FALSE`, the original table is returned unmodified.
#' @param revigo_cutoff The similarity cutoff for Revigo (e.g., 0.7).
#' @param species_taxon The NCBI species taxon ID for Revigo (e.g., 9606 for Human).
#'
#' @return A data frame of enrichment results filtered to include only the representative
#'   terms identified by Revigo. If Revigo returns no terms for a group, a warning is
#'   issued, and the original data for that group may be affected depending on the join.
#'
#' @importFrom dplyr group_by nest ungroup mutate select unnest rename left_join filter
#' @importFrom purrr map
#' @importFrom rlang set_names
#' @export
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



#' @title Filter GSEA/ORA Results Using Revigo (clusterProfiler `scholar` specific)
#' @description A variant of `filterResultsWithRevigo` tailored for results from
#' `clusterProfiler::GSEA` or `clusterProfiler::enricher` (often used with the `scholar` package).
#' It groups results, sends term IDs to Revigo, and filters based on the summary.
#'
#' @param enriched_results_tbl A data frame of enrichment results, typically from `clusterProfiler`.
#'   Must contain an `ID` column with term identifiers.
#' @param added_columns A character vector of column names to group by before sending data to Revigo.
#'   `go_type` is also used for grouping.
#' @param is_run_revigo A logical flag. If `TRUE` (default), Revigo filtering is performed.
#' @param revigo_cutoff The similarity cutoff for Revigo.
#' @param species_taxon The NCBI species taxon ID for Revigo.
#'
#' @return A data frame of enrichment results filtered by Revigo.
#'
#' @importFrom dplyr group_by nest ungroup mutate select unnest rename left_join filter
#' @importFrom purrr map
#' @importFrom rlang set_names
#' @export
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


#' @title Save Filtered Functional Enrichment Table
#' @description Saves a functional enrichment results table to both .tab (tab-separated)
#' and .xlsx (Excel) formats. It saves both a filtered version (based on gene set size)
#' and the complete, unfiltered table.
#'
#' @details
#' The function truncates cell content for specified columns in the Excel file to avoid
#' exceeding Excel's cell character limit (32,760 characters).
#'
#' Four files will be created:
#' - `[file_name].tab`: Filtered data.
#' - `[file_name].xlsx`: Filtered data, with long columns truncated.
#' - `[file_name]_unfiltered.tab`: Complete data.
#' - `[file_name]_unfiltered.xlsx`: Complete data, with long columns truncated.
#'
#' @param enriched_results_tbl The data frame of enrichment results to save.
#' @param set_size_min The minimum gene set size to filter by for the primary output files.
#' @param set_size_max The maximum gene set size to filter by for the primary output files.
#' @param results_dir The directory where the files will be saved.
#' @param file_name The base name for the output files (without extension).
#' @param list_of_columns_to_trim A character vector of column names whose content should
#'   be truncated for the Excel export. Defaults to `c("gene_symbol")`.
#'
#' @return This function is called for its side effect of writing files and does not return a value.
#'
#' @importFrom vroom vroom_write
#' @importFrom writexl write_xlsx
#' @importFrom dplyr filter mutate across one_of
#' @export
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


#' @title Evaluate Optimal Gene Set Size Parameters
#' @description Generates a line plot to help evaluate the effect of different
#' minimum and maximum gene set size parameters on the number of significantly
#' enriched terms found.
#'
#' @details
#' The plot shows the number of significant pathways (y-axis) found for each
#' combination of `min_set_size` and `max_set_size` (x-axis). Lines are grouped by
#' comparison and faceted by gene set type (e.g., "positive_list-BP"). This visualization
#' helps in choosing size parameters that yield a reasonable number of results.
#'
#' @param enrichment_results_tble A data frame of enrichment results, containing data
#'   from runs with multiple `min_set_size` and `max_set_size` values.
#' @param added_columns A character vector of additional columns to use for grouping comparisons.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise ungroup mutate if_else
#' @importFrom tidyr unite
#' @export
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


#' @title Draw a List of Functional Enrichment Heatmaps
#' @description A wrapper function that processes an enrichment results table, clusters
#' the pathways, and generates a list of heatmaps, one for each GO type (or other category).
#'
#' @details
#' This function first filters the results by the specified `set_size_min` and
#' `set_size_max`. It then calls `clusterPathways` to handle duplicates and order the
#' pathways. Finally, it iterates through each `go_type` and calls `getEnrichmentHeatmap`
#' to generate a plot for it.
#'
#' @param enriched_results_tbl The data frame of enrichment results.
#' @param added_columns A character vector of column names to be used in clustering and analysis.
#' @param set_size_min The minimum gene set size to filter the results by.
#' @param set_size_max The maximum gene set size to filter the results by.
#' @param x_axis The unquoted column name to use for the x-axis of the heatmaps.
#' @param analysis_column The unquoted column name to unite comparison columns into.
#' @param facet_by_column Optional unquoted column name to facet the heatmaps by.
#' @param remove_duplicted_entries How to handle duplicate pathway entries (passed to `clusterPathways`).
#' @param xaxis_levels Optional character vector to specify the order of x-axis items.
#' @param scales The `scales` argument for faceting (e.g., "free", "fixed").
#'
#' @return A named list of ggplot objects, where each element is a heatmap for a specific GO type.
#'
#' @importFrom dplyr filter group_by arrange mutate ungroup distinct
#' @importFrom purrr pmap
#' @importFrom tidyr unite
#' @export
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





#' @title Draw a List of Functional Enrichment Heatmaps (Scholar Version)
#' @description A wrapper function that processes an enrichment results table from
#' `clusterProfiler` (often used with `scholar`), clusters the pathways, and generates
#' a list of heatmaps, one for each GO type or category. This version does not
#' filter by gene set size, assuming it has been done upstream.
#'
#' @details
#' This function calls `clusterPathways` to handle duplicates and order the
#' pathways. It then iterates through each `go_type` and calls `getEnrichmentHeatmap`
#' to generate a plot.
#'
#' @param enriched_results_tbl The data frame of enrichment results.
#' @param added_columns A character vector of column names to be used in clustering and analysis.
#' @param x_axis The unquoted column name to use for the x-axis of the heatmaps.
#' @param analysis_column The unquoted column name to unite comparison columns into.
#' @param facet_by_column Optional unquoted column name to facet the heatmaps by.
#' @param remove_duplicted_entries How to handle duplicate pathway entries (passed to `clusterPathways`).
#' @param xaxis_levels Optional character vector to specify the order of x-axis items.
#' @param scales The `scales` argument for faceting (e.g., "free", "fixed").
#'
#' @return A named list of ggplot objects, where each element is a heatmap for a specific GO type.
#'
#' @importFrom dplyr group_by arrange mutate ungroup distinct
#' @importFrom purrr pmap
#' @importFrom tidyr unite
#' @export
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



#' @title Save a List of ggplot Heatmaps to Files
#' @description Iterates through a list of ggplot objects and saves each one to a
#' file, typically in PNG or PDF format. It allows specifying different dimensions for each plot.
#'
#' @param list_of_heatmaps A named list of ggplot objects. The names of the list
#'   elements are used in the output filenames.
#' @param results_dir The directory where the plots will be saved.
#' @param file_name The base name for the output plot files. The final filename will be
#'   `[file_name]_[plot_name]`.
#' @param plot_width The width of the plot(s). Can be a single numeric value (applied
#'   to all plots) or a numeric vector with the same length as `list_of_heatmaps`.
#' @param plot_height The height of the plot(s). Can be a single value or a vector.
#'
#' @return This function is called for its side effect of writing files and does not return a value.
#'
#' @importFrom purrr pwalk
#' @export
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

###------------------------------------------------------------------------------------------------------------------------

#' @title Create a Bar Plot for Enriched Pathways
#' @description Generates a horizontal bar plot of enriched pathways or GO terms,
#' with bar length representing significance (`-log10(qvalue)`). It handles
#' up-regulated, down-regulated, and shared terms, displaying them with different colors.
#'
#' @details
#' This function first resolves duplicate entries (e.g., a pathway found in both
#' up- and down-regulated lists) based on the `remove_duplicted_entries` parameter.
#' It then orders the terms by significance and regulation status for a clean visual presentation.
#'
#' @param input_table A data frame of enrichment results.
#' @param input_go_type An optional character string to filter the results to a single GO type.
#' @param remove_duplicted_entries How to handle duplicates (`"merge"` or `"delete"`).
#' @param added_columns Additional columns to consider when identifying duplicates.
#'
#' @return A ggplot object representing the bar plot.
#'
#' @import ggplot2
#' @importFrom dplyr mutate case_when group_by summarise ungroup filter anti_join inner_join bind_rows distinct pull arrange
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



#' @title Generate and Save GO Term Enrichment Bar Plots
#' @description A wrapper function that creates and saves bar plots for GO enrichment
#' results. It generates a separate plot for each GO type (e.g., BP, CC, MF) found
#' in the input data.
#'
#' @param input_table The data frame of filtered and Revigo-summarized enrichment results.
#' @param output_dir The directory to save the plots in.
#' @param analysis_type A string used as a prefix in the output filename (e.g., "GO").
#' @param file_suffix A character vector of file extensions (e.g., `c("png", "pdf")`)
#'   for saving the plots.
#' @param width The width of the output plots.
#' @param height The height of the output plots.
#'
#' @return This function is called for its side effect of writing files and does not return a value.
#'
#' @importFrom dplyr distinct arrange pull
#' @importFrom purrr map walk2
#' @importFrom ggplot2 ggsave
#' @importFrom rlang partial
#' @export
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

#' @title Create a Word Cloud Data Frame
#' @description Processes a vector of text (e.g., GO term names) to generate a
#' data frame of word frequencies, suitable for creating a word cloud.
#'
#' @details This function performs standard text mining pre-processing steps:
#' - Removes numbers and punctuation.
#' - Converts text to lower case.
#' - Removes common English stopwords.
#' The process is based on an article by Céline Van den Rul.
#'
#' @param text_list A character vector of text documents (e.g., GO term descriptions).
#'
#' @return A data frame with two columns: `word` and `freq` (frequency of the word).
#'
#' @importFrom tm Corpus VectorSource tm_map removeNumbers removePunctuation stripWhitespace content_transformer removeWords stopwords TermDocumentMatrix
#' @export
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

#' @title Clean Duplicate Enrichment Results
#' @description This function resolves duplicate entries in an enrichment result table
#' where the same pathway/term is enriched in both up- and down-regulated gene sets.
#' It merges these duplicates into a single "Both" category, keeping the most significant p-value.
#'
#' @param input_table The data frame of enrichment results.
#' @param pathway_column The unquoted column name for the pathway/term description.
#' @param fdr_column The unquoted column name for the false discovery rate or q-value.
#' @param gene_set_column The unquoted column name for the gene set (e.g., indicating up/down regulation).
#'
#' @return A data frame with duplicates resolved and an added `neg_log_10_fdr` column,
#'   ordered by significance.
#'
#' @importFrom dplyr group_by summarise ungroup filter inner_join mutate anti_join bind_rows distinct pull arrange desc
#' @importFrom stringr str_detect
#' @export
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

#' @title Plot Enrichment Bar Plot
#' @description Creates a horizontal bar plot from enrichment results. The bars represent
#' enriched terms, their length corresponds to significance (`-log10(FDR)`), and
#' they are colored by gene set (e.g., up-regulated, down-regulated).
#'
#' @details This function first uses `cleanDuplicatesEnrichment` to resolve duplicate
#' entries. It then creates a ggplot bar plot with pathways ordered by significance.
#'
#' @param input_table The data frame of enrichment results.
#' @param pathway_column The unquoted column name for the pathway/term description.
#' @param fdr_column The unquoted column name for the false discovery rate or q-value.
#' @param gene_set_column The unquoted column name for the gene set.
#' @param xlab_string The label for the x-axis.
#' @param ylab_string The label for the y-axis.
#' @param legend_title The title for the fill legend.
#' @param legend_colours A character vector of colors for the fill aesthetic.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom dplyr pull distinct
#' @export
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


########################

## Functions taken from clusterProfiler package to enable drawing of cnet (e.g. Pathways and proteins association network).
# Reference: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
#            https://github.com/YuLab-SMU/clusterProfiler
## Functions taken on 1st November 2023



#' @title Convert a Named List to a Two-Column Data Frame
#' @description This internal helper function converts a named list (like a list of
#' gene sets) into a two-column data frame. The first column holds the names of the
#' list elements (categories), and the second column holds the unlisted values (genes).
#'
#' @param inputList A named list.
#'
#' @return A data frame with two columns: `categoryID` and `Gene`.
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


#' @title Convert a List to an igraph Graph Object
#' @description This function takes a named list (representing categories and items,
#' such as pathways and genes) and converts it into an undirected `igraph` graph object.
#'
#' @param inputList A named list, which will be converted to a two-column data frame
#'   by `list2df`.
#'
#' @return An `igraph` graph object.
#'
#' @importFrom igraph graph.data.frame
#' @export
list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}


#' @title Generate a Deprecation Message for a Parameter
#' @description Creates a standardized warning message for a deprecated parameter,
#' guiding the user on the new, correct syntax.
#'
#' @param parameter The name of the deprecated parameter (string).
#' @param params_df A data frame mapping old parameter names to their new list-based syntax.
#'   It must contain columns: `listname`, `present`, `original`.
#'
#' @return A character string containing the formatted warning message.
#' @export
get_param_change_message <- function(parameter, params_df) {
  paste0("Use '", params_df[parameter, "listname"],
         " = list(", params_df[parameter, "present"],
         " = your_value)' instead of '", params_df[parameter, "original"],
         "'.\n The ", params_df[parameter, "original"],
         " parameter will be removed in the next version.")
}

#' @title Add Alpha/Transparency to Nodes in a ggraph Plot
#' @description Modifies a `ggraph` plot object to set the alpha (transparency)
#' of nodes. It can highlight specified categories and genes by giving them a
#' different alpha value than other nodes.
#'
#' @param p A `ggraph` plot object.
#' @param hilight_category A character vector of category node names to highlight.
#' @param hilight_gene A character vector of gene node names to highlight.
#' @param alpha_nohilight The alpha value for non-highlighted nodes.
#' @param alpha_hilight The alpha value for highlighted nodes.
#'
#' @return The modified `ggraph` plot object with an `alpha` column added to its data.
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


#' @title Get Default enrichplot Colors
#' @description Retrieves the default color palette used in `enrichplot` visualizations.
#' It returns a two-color vector for gradients or a three-color vector for divergent scales.
#' Users can override the default by setting the `enrichplot.colours` option.
#'
#' @param n The number of colors to return (either 2 or 3).
#'
#' @return A character vector of hex color codes.
#' @export
get_enrichplot_color <- function(n = 2) {
  colors <- getOption("enrichplot.colours")
  if (!is.null(colors)) return(colors)

  if (n != 2 && n != 3) stop("'n' should be 2 or 3")

  colors = c("#e06663", "#327eba")
  if (n == 2) return(colors)

  if (n == 3) return(c(colors[1], "white", colors[2]))
}

#' @title Set a Color Scale for enrichplot
#' @description A flexible helper function to create a ggplot2 color scale (e.g.,
#' `scale_color_continuous`, `scale_fill_gradientn`) for use in `enrichplot`
#' visualizations. It can create continuous, 2-point gradient, or n-point gradient scales.
#'
#' @param colors A character vector of color hex codes. The number of colors determines
#'   the type of scale created (2 for continuous, 3 for diverging, >3 for n-color gradient).
#' @param type The type of scale to create: "color", "colour", or "fill".
#' @param name The title for the color scale legend.
#' @param .fun An optional custom scaling function to use instead of the ggplot2 defaults.
#' @param ... Additional arguments passed to the ggplot2 scaling function.
#'
#' @return A ggplot2 scale object.
#'
#' @importFrom ggplot2 guide_colorbar
#' @importFrom utils getFromNamespace
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

#' @title Add Node Labels to a ggraph Plot
#' @description Adds text labels to nodes in a `ggraph` plot, using `ggrepel` to
#' prevent overlapping text. It can apply a "shadow" effect by drawing a white
#' background for the text.
#'
#' @param p A `ggraph` plot object.
#' @param data The data frame containing node information (usually `p$data`).
#' @param label_size_node The base font size for the labels.
#' @param cex_label_node A scaling factor for the label size.
#' @param shadowtext A logical value indicating whether to apply the shadow text effect.
#'
#' @return The modified `ggraph` plot object with a `geom_node_text` layer.
#'
#' @importFrom ggraph geom_node_text
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

#' @title Get ggrepel Segment Size Option
#' @description A helper function to retrieve the global option for `ggrepel.segment.size`,
#' which controls the thickness of the lines connecting labels to their points.
#'
#' @param default The default value to return if the option is not set. Defaults to 0.2.
#'
#' @return The value of the `ggrepel.segment.size` option, or the default.
#' @export
get_ggrepel_segsize <- function(default = 0.2) {
  getOption("ggrepel.segment.size", default = default)
}

#' @title Create a Category-Gene Network Plot (Edited)
#' @description An edited version of `clusterProfiler::cnetplot` for visualizing the
#' relationships between gene sets (categories) and genes. It creates a network
#' where nodes are either categories or genes, and an edge connects a gene to a
#' category it belongs to.
#'
#' @details This function is highly customizable, allowing control over layout, colors,
#' node sizes, and labels. It can also map fold change values to gene node colors.
#' It handles parameter deprecation gracefully by mapping old parameters to a new,
#' list-based system (`color.params`, `cex.params`, `hilight.params`).
#'
#' @param geneSets A named list of character vectors, where each element is a gene set
#'   and the name is the category.
#' @param showCategory The number of categories to display, or a character vector of
#'   category names to display.
#' @param foldChange A named numeric vector of fold changes for genes.
#' @param layout The layout algorithm to use for the network (e.g., "kk", "fr").
#' @param colorEdge A logical value. If `TRUE`, edges are colored by category.
#' @param circular A logical value. If `TRUE`, a circular layout is used.
#' @param node_label Which nodes to label: "all", "gene", "category", or "none".
#' @param cex_category,cex_gene,cex_label_category,cex_label_gene Deprecated parameters for node and label sizes.
#' @param color_category,color_gene Deprecated parameters for node colors.
#' @param shadowtext Which labels to apply a shadow effect to: "all", "gene", "category", or "none".
#' @param color.params A list controlling color aesthetics:
#'   - `foldChange`: Named numeric vector for gene colors.
#'   - `edge`: Logical, whether to color edges by category.
#'   - `category`: Color for category nodes.
#'   - `gene`: Default color for gene nodes.
#' @param cex.params A list controlling size aesthetics:
#'   - `category_node`, `gene_node`: Scaling factors for node sizes.
#'   - `category_label`, `gene_label`: Scaling factors for label sizes.
#' @param hilight.params A list controlling highlighting:
#'   - `category`: Character vector of categories to highlight.
#'   - `alpha_hilight`, `alpha_no_hilight`: Alpha values for highlighted and non-highlighted nodes/edges.
#' @param ... Additional parameters.
#'
#' @return A `ggraph` plot object.
#'
#' @import ggraph
#' @importFrom ggplot2 theme_void guides guide_legend
#' @importFrom ggraph geom_edge_arc geom_edge_link geom_node_point scale_edge_width
#' @importFrom ggnewscale new_scale_color
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

##################################S


#' @title Make Fold Change Vector Readable
#' @description A helper function to translate gene identifiers in a fold change vector
#' to gene symbols if the input enrichment object `x` is in a "readable" format
#' (i.e., contains a mapping from gene IDs to symbols).
#'
#' @param x An enrichment result object (e.g., `enrichResult` or `gseaResult`) which
#'   may contain a `gene2Symbol` mapping.
#' @param foldChange A named numeric vector of fold changes, where names are gene IDs.
#'
#' @return A named numeric vector of fold changes with names translated to gene symbols,
#'   if possible. Otherwise, the original vector is returned.
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


#' @title Update `showCategory` Parameter
#' @description A helper function to process the `showCategory` parameter. If it's a
#' number, it ensures it's not larger than the total number of categories. If it's a
#' character vector of category names, it filters them to ensure they exist in the
#' results.
#'
#' @param x An enrichment result object or a list of gene sets.
#' @param showCategory A numeric value indicating the number of top categories to show,
#'   or a character vector of category names to show.
#'
#' @return An integer or a character vector, representing the validated number or names
#'   of categories to display.
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

#' @title Extract Gene Sets from an Enrichment Object
#' @description Extracts the gene sets (e.g., pathways and their associated genes) from
#' an enrichment result object or a list. It can select a specific number of top
#' categories or a specific set of categories by name.
#'
#' @param x An enrichment result object (`enrichResult`, `gseaResult`) or a named list of gene sets.
#' @param n The number of top gene sets to extract, or a character vector of gene set names to extract.
#'
#' @return A named list where each element is a character vector of genes, and the
#'   names are the gene set descriptions.
#'
#' @importFrom DOSE geneInCategory
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

########

#' @title Protein Pathway Enrichment Analysis Helper
#' @description This function performs pathway enrichment analysis for a single
#' differential expression contrast. It identifies significant up- and down-regulated
#' proteins and runs GO enrichment analysis against a background set.
#'
#' @details
#' This is a helper for `enrichProteinsPathways`. Its main steps are:
#' 1.  Clean and prepare protein IDs.
#' 2.  Download or load cached UniProt annotation data.
#' 3.  Identify significant proteins based on specified thresholds.
#' 4.  Run GO enrichment for positive and negative sets separately using `runOneGoEnrichmentInOutFunction`.
#' 5.  Combine and return the results.
#'
#' @param de_analysis_results A list object with `de_proteins_wide` and `significant_rows`.
#' @param organism_taxid NCBI taxonomy ID for the organism (e.g., "9606").
#' @param min_gene_set_size Minimum number of genes in a gene set.
#' @param max_gene_set_size Maximum number of genes in a gene set.
#' @param p_val_thresh P-value threshold for enrichment significance.
#' @param protein_p_val_thresh FDR/q-value threshold for significant proteins.
#' @param cache_dir Directory to store cached UniProt data.
#' @param cache_file Name of the cache file.
#' @param use_cached Logical, whether to use cached data.
#' @param protein_id_delimiter Delimiter in the protein ID column.
#' @param protein_id_column Name of the protein ID column.
#'
#' @return A data frame with enrichment results for the contrast, including a `directionality` column.
#'
#' @importFrom UniProt.ws UniProt.ws
#' @importFrom GO.db GOTERM
#' @importFrom AnnotationDbi Term Ontology
#' @importFrom dplyr select distinct mutate filter pull left_join case_when bind_rows sym
#' @importFrom tidyr pivot_longer separate_rows
#' @importFrom purrr map_chr
#' @importFrom stringr str_split
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

#' @title Perform Protein Pathway Enrichment Across Multiple Contrasts
#' @description Main wrapper to run protein GO pathway enrichment. It iterates over a list
#' of differential expression results and performs enrichment analysis for each contrast.
#'
#' @param de_analysis_results_list A named list of DE results objects.
#' @param taxon_id NCBI taxonomy ID for the organism.
#' @param protein_id_delimiter Delimiter in the protein ID column.
#' @param protein_p_val_thresh FDR/q-value threshold for significant proteins.
#' @param min_gene_set_size Minimum gene set size.
#' @param max_gene_set_size Maximum gene set size.
#' @param p_val_thresh P-value threshold for enrichment.
#' @param cache_dir Directory for cached data.
#' @param cache_file Filename for the cache file.
#' @param use_cached Logical, whether to use cached data.
#'
#' @return A single data frame of combined enrichment results from all contrasts.
#'
#' @importFrom dplyr bind_rows
#' @importFrom purrr map set_names
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

#' @title Download and Cache UniProt Annotation Data
#' @description Downloads protein annotation from UniProt for a list of protein IDs,
#' with caching to avoid re-downloads. If a cache exists, it only fetches data for new IDs.
#'
#' @param protein_ids A data frame with a `uniprot_acc` column.
#' @param cache_file Full path to the RDS cache file.
#' @param uniprot_handle An active `UniProt.ws` handle.
#' @param protein_id_delimiter Not used directly, for consistency.
#'
#' @return A data frame of UniProt annotations, also saved to the cache.
#'
#' @importFrom dplyr filter bind_rows rename mutate
#' @importFrom purrr map_dbl
#' @export
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

#' @title Convert UniProt GO IDs to Terms and Types
#' @description Parses a column of semicolon-separated GO IDs, expands the table to have
#' one GO ID per row, and annotates it with the corresponding term and type (BP, CC, MF).
#'
#' @param uniprot_dat A data frame from a UniProt query.
#' @param uniprot_id_column Unquoted column name for UniProt accessions.
#' @param go_id_column Unquoted column name for GO IDs.
#' @param gene_name_column Unquoted column name for gene names.
#' @param sep Separator for GO IDs.
#' @param goterms A GO term dictionary from `AnnotationDbi::Term(GO.db::GOTERM)`.
#' @param gotypes A GO type dictionary from `AnnotationDbi::Ontology(GO.db::GOTERM)`.
#'
#' @return A long-format data frame with columns `uniprot_acc`, `go_id`, `go_term`, `go_type`.
#'
#' @importFrom tidyr separate_rows
#' @importFrom dplyr distinct filter mutate rename
#' @importFrom purrr map_chr
#' @importFrom GO.db GOTERM
#' @importFrom AnnotationDbi Term Ontology
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

#' @title Clean Protein Identifiers
#' @description Internal helper to clean protein IDs, typically by removing
#' isoform suffixes (e.g., "-1") from UniProt accessions.
#'
#' @param ids A character vector of protein IDs.
#'
#' @return A character vector of cleaned protein IDs.
#'
#' @importFrom stringr str_split str_replace_all
#' @importFrom purrr map_chr
#' @export
.cleanProteinIds <- function(ids) {
  ids |>
    stringr::str_split("-") |>
    purrr::map_chr(1) |>
    stringr::str_replace_all("-\\d+$", "")
}
