

# protein_quant_matrix <- ruv_normalized_for_de_analysis_obj@protein_data
#
# log_fold_change_table <- de_analysis_results_list$de_proteins_long
#
# uniprot_tbl <- uniprot_dat_cln
#
# annotation_column <- "go_term"
# min_gene_set_size <- 4
# max_gene_set_size <- 300
# protein_id <- "uniprot_acc"
# annotation_id <- "go_term"
# aspect_column <- "go_type"
# annotation_type <- "go_enrichment"
# protein_id_lookup_column <- "Entry"
# gene_symbol_column <- "Gene.Names"
# protein_id_or_gene_symbol <- "gene_symbol"
# fdr_column_name <- "q.mod"
# log_fc_column_name <- "log2FC"
# protein_p_val_thresh <- 0.05
# id_column_for_enrichment_test <- "uniprot_acc_first"



# runFunctionalEnrichment (protein_quant_matrix = ruv_normalized_for_de_analysis_obj@protein_data
#                          , de_analysis_results_list$de_proteins_long
#                          , uniprot_dat_cln)


runFunctionalEnrichment <- function(protein_quant_matrix
                                    , log_fold_change_table
                                    , uniprot_tbl
                                    ,output_dir="../../results/proteomics/de_proteins_list"
                                    ,annotation_file="../../data/proteomics/UniProt/go_terms_table_python_all.tab"
                                    ,dictionary_file="../../data/proteomics/UniProt/go_terms_table_python_all.tab"
                                    ,annotation_column="go_term"
                                    ,min_gene_set_size=4
                                    ,max_gene_set_size=300
                                    ,protein_id="uniprot_acc"
                                    ,annotation_id="go_term"
                                    ,aspect_column="go_type"
                                    ,annotation_type="go_enrichment"
                                    ,protein_id_lookup_column="Entry"
                                    ,gene_symbol_column="Gene.Names"
                                    ,protein_id_or_gene_symbol="gene_symbol"

                                    ) {

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print ("Read proteomics abundance data.")

  proteins_tbl_orig <- log_fold_change_table |>
    mutate( uniprot_acc_first = purrr::map_chr( !!sym( protein_id), ~str_split(., ":") |> map_chr(1))) |>
    mutate( uniprot_acc_first = str_replace_all( uniprot_acc_first, "-\\d+$", ""))   # Strip away isoform information




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print ("Compile positive proteins list.")

uniprot_tbl <- NA
proteins_tbl_updated <- proteins_tbl_orig
if(isArgumentDefined(args, "uniprot_to_gene_symbol_file")) {
  # uniprot_to_gene_symbol_file <- "/home/ubuntu/Workings/2021/ALPK1_BMP_06/Data/UniProt/data.tab"

  # protein_id_lookup_column <- "Entry"
  # gene_symbol_column <- "Gene names"

  # Clean up protein ID to gene symbol table

  proteins_tbl_updated <- proteins_tbl_orig |>
    left_join( uniprot_tbl
               , by = join_by(  !!sym(protein_id) == !!sym(protein_id_lookup_column)) ) |>
    dplyr::rename( UNIPROT_GENENAME = !!sym(gene_symbol_column) ) |>
    dplyr::rename( `PROTEIN-NAMES` = "Protein_names" ) |>
    dplyr::mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x) { ifelse( !is.na(x), str_split_1(x, "\\s+"), NA_character_) }) ) |>
    mutate( gene_name_first = purrr::map_chr( UNIPROT_GENENAME, ~str_split(., ":") |> map_chr(1)))  |>
    mutate( protein_name_first = purrr::map_chr( `PROTEIN-NAMES`, ~str_split(., ":") |> map_chr(1))) |>
    mutate( gene_name_first = ifelse( is.na(gene_name_first), uniprot_acc_first, gene_name_first) )

}


  if( "UNIPROT_GENENAME" %in% colnames( proteins_tbl_updated))  {

    positive_proteins <- proteins_tbl_updated |>
      dplyr::filter( !!rlang::sym(fdr_column_name) < protein_p_val_thresh &   !!rlang::sym(log_fc_column_name) > 0 ) |>
      group_by(comparison, uniprot_acc_first, gene_name_first, protein_name_first) |>
      summarise( max_norm_logFC = max(!!rlang::sym(log_fc_column_name))) |>
      ungroup() |>
      arrange( comparison, desc(max_norm_logFC  ) )


  } else {

    positive_proteins <- proteins_tbl_updated |>
      dplyr::filter( !!rlang::sym(fdr_column_name) < protein_p_val_thresh &   !!rlang::sym(log_fc_column_name) > 0 ) |>
      group_by(comparison, uniprot_acc_first ) |>
      summarise( max_norm_logFC = max(!!rlang::sym(log_fc_column_name))) |>
      ungroup() |>
      arrange( comparison, desc(max_norm_logFC  ) )
  }



  vroom::vroom_write( positive_proteins,
                      file.path( output_dir,
                                 "all_proteins_with_positive_logFC.tab" ),
                      col_names=FALSE)



  list_of_comparisons <- positive_proteins |> distinct( comparison) |> pull( comparison)



  purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(output_dir,. )) )



  purrr::walk( list_of_comparisons, function( input_comparison){

    positive_proteins |>
      dplyr::filter( comparison == input_comparison) |>
      dplyr::select( !!sym( id_column_for_enrichment_test)) |>
      vroom::vroom_write( file.path( output_dir,
                                     input_comparison,
                                     "all_proteins_with_positive_logFC.tab" ),
                          col_names=FALSE)

  } )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# norm_phos_logFC
print ("Compile negative proteins list.")

if( "UNIPROT_GENENAME" %in% colnames( proteins_tbl_updated))  {

  negative_proteins <- proteins_tbl_updated |>
    dplyr::filter( !!rlang::sym(fdr_column_name) < protein_p_val_thresh & !!rlang::sym(log_fc_column_name)  < 0 ) |>
    group_by(comparison, uniprot_acc_first, gene_name_first, protein_name_first) |>
    summarise( min_norm_logFC = min(!!rlang::sym(log_fc_column_name))) |>
    ungroup() |>
    arrange( comparison, min_norm_logFC  )

} else {

  negative_proteins <- proteins_tbl_updated |>
    dplyr::filter( !!rlang::sym(fdr_column_name) < protein_p_val_thresh & !!rlang::sym(log_fc_column_name)  < 0 ) |>
    group_by(comparison, uniprot_acc_first ) |>
    summarise( min_norm_logFC = min(!!rlang::sym(log_fc_column_name))) |>
    ungroup() |>
    arrange( comparison, min_norm_logFC  )
}

vroom::vroom_write( negative_proteins,
                    file.path( output_dir,
                               "all_proteins_with_negative_logFC.tab" ),
                    col_names=FALSE )

list_of_comparisons <- negative_proteins |> distinct( comparison) |> pull( comparison)

purrr::walk( list_of_comparisons, ~createDirIfNotExists( file.path(output_dir,. )) )

purrr::walk( list_of_comparisons, function( input_comparison){

  negative_proteins |>
    dplyr::filter( comparison == input_comparison) |>
    dplyr::select( !!sym( id_column_for_enrichment_test) ) |>
    vroom::vroom_write( file.path( output_dir,
                                   input_comparison,
                                   "all_proteins_with_negative_logFC.tab" ),
                        col_names=FALSE)

} )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print ("Compile background proteins list.")


background_proteins <- proteins_tbl_updated |>
  dplyr::distinct(!!sym( id_column_for_enrichment_test))


vroom::vroom_write( background_proteins,
                    file.path(output_dir, "background_proteins.tab" ),
                    col_names=FALSE )


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print ("Compile annotation ID to annotation term name dictionary.")

## Tidy up the annotation ID to annotation term name dictionary

dictionary <- vroom::vroom( dictionary_file )

id_to_annotation_dictionary <- buildAnnotationIdToAnnotationNameDictionary( input_table=dictionary,
                                                                            annotation_column = !!rlang::sym(annotation_column),
                                                                            annotation_id_column = !!rlang::sym(annotation_id))

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print ("Preparing the enrichment test.")

## preparing the enrichment test
go_annot_orig <- vroom::vroom(  annotation_file   )

go_annot <- go_annot_orig
if(isArgumentDefined(args, "uniprot_to_gene_symbol_file")) {

  # Clean up protein ID to gene symbol table

  go_annot <- go_annot_orig |>
    left_join( uniprot_tbl |>
                 dplyr::select( !!sym(protein_id_lookup_column), !!sym(gene_symbol_column))
               , by = join_by(  !!sym(protein_id) == !!sym(protein_id_lookup_column)) ) |>
    dplyr::rename( UNIPROT_GENENAME = !!sym(gene_symbol_column) ) |>
    dplyr::mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x) { ifelse( !is.na(x), str_split_1(x, "\\s+"), NA_character_) }) ) |>
    mutate( gene_name_first = purrr::map_chr( UNIPROT_GENENAME, ~str_split(., ":") |> map_chr(1)))

}


background_list <- background_proteins

#print( min_gene_set_size)
min_gene_set_size_list <- parseNumList(min_gene_set_size)
max_gene_set_size_list <- parseNumList(max_gene_set_size)

list_of_comparisons <- negative_proteins |>
  bind_rows(positive_proteins) |>
  distinct(comparison) |>
  pull(comparison)

# Tidy up GO aspect list, marked as null if not using GO terms
go_aspect_list <- NA
if(!is.null(aspect_column)) {
  go_aspect_list <- go_annot |>
    dplyr::filter( !is.na(!!rlang::sym( aspect_column) )) |>
    distinct( !!rlang::sym( aspect_column) ) |>
    pull( !!rlang::sym( aspect_column) ) # c("C", "F", "P")
} else {
  go_aspect_list <- NA
}

# print( paste("is.na(go_aspect_list) =", is.na(go_aspect_list)) )

list_of_genes_list <- list( negative_list=negative_proteins,
                            positive_list=positive_proteins)

input_params <- expand_grid(
  names_of_genes_list = names( list_of_genes_list),
  go_aspect=go_aspect_list,
  input_comparison = list_of_comparisons,
  min_size = min_gene_set_size_list,
  max_size = max_gene_set_size_list)

input_params_updated <- input_params |>
  mutate( input_table = purrr::map(names_of_genes_list,  \(x) list_of_genes_list[[x]] ))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print ("Run Enrichment.")

runOneGoEnrichmentInOutFunctionPartial <- purrr::partial ( runOneGoEnrichmentInOutFunction,
                                                           comparison_column = comparison,
                                                           protein_id_column = !!sym( id_column_for_enrichment_test),
                                                           go_annot = go_annot,
                                                           background_list = background_list,
                                                           id_to_annotation_dictionary=id_to_annotation_dictionary,
                                                           annotation_id=!!rlang::sym(annotation_id),
                                                           protein_id=!!rlang::sym(id_column_for_enrichment_test), # protein_id
                                                           aspect_column=aspect_column,
                                                           p_val_thresh=p_val_thresh
                                                           , get_cluster_profiler_object = TRUE )

i <- 1
runOneGoEnrichmentInOutFunctionPartial(
  names_of_genes_list = input_params_updated$names_of_genes_list[[i]],
  input_table=input_params_updated$input_table[[i]],
  go_aspect=input_params_updated$go_aspect[[i]],
  input_comparison=input_params_updated$input_comparison[[i]],
  min_gene_set_size=input_params_updated$min_size[[i]],
  max_gene_set_size=input_params_updated$max_size[[i]])

enrichment_combined <- input_params_updated  |>
  mutate( enrichment_temp = purrr::pmap( list( names_of_genes_list
                                               , input_table
                                               , go_aspect
                                               , input_comparison
                                               , min_size
                                               , max_size) ,
                                         \(names_of_genes_list
                                           ,input_table
                                           , go_aspect
                                           , input_comparison
                                           , min_size
                                           , max_size) {runOneGoEnrichmentInOutFunctionPartial(
                                             names_of_genes_list = names_of_genes_list,
                                             input_table=input_table,
                                             go_aspect=go_aspect,
                                             input_comparison=input_comparison,
                                             min_gene_set_size=min_size,
                                             max_gene_set_size=max_size)} ) )    |>
  dplyr::mutate( is_logical  = map_chr(enrichment_temp, \(x) { typeof(x) } )  ) |>
  dplyr::filter( is_logical == "list" )  |>
  mutate( enrichment_objects = purrr::map(enrichment_temp, \(x){ x$enrichment_object} )  ) |>
  mutate( enrichment_results = purrr::map(enrichment_temp, \(x){ x$enrichment_result} )  )  |>
  dplyr::select(-enrichment_temp)


enrichment_objects <-  enrichment_combined |>
  dplyr::select(names_of_genes_list, go_aspect, input_comparison, min_size, max_size, enrichment_objects)

enrichment_result <- enrichment_combined |>
  dplyr::select(-names_of_genes_list) |>
  unnest(cols = c("enrichment_results"))  |>
  dplyr::select(-enrichment_objects)

if(is.null(enrichment_result) |
   nrow(enrichment_result) == 0 ) {
  warnings("No enriched terms were identified.")
} else {


  enrichment_result_add_gene_symbol <- NA
  ## Convert Uniprot accession to gene names
  if(isArgumentDefined(args, "uniprot_to_gene_symbol_file") &&
     protein_id_or_gene_symbol== "protein_id" ) {
    # uniprot_to_gene_symbol_file <- "/home/ubuntu/Workings/2021/ALPK1_BMP_06/Data/UniProt/data.tab"

    # protein_id_lookup_column <- "Entry"
    # gene_symbol_column <- "Gene names"

    # Clean up protein ID to gene symbol table
    uniprot_tbl <- vroom::vroom( file.path( uniprot_to_gene_symbol_file))

    uniprot_to_gene_symbol_dict <- getUniprotAccToGeneSymbolDictionary( uniprot_tbl,
                                                                        !!rlang::sym(protein_id_lookup_column),
                                                                        !!rlang::sym(gene_symbol_column),
                                                                        !!rlang::sym(protein_id) )

    convertProteinAccToGeneSymbolPartial <- purrr::partial( convertProteinAccToGeneSymbol,
                                                            dictionary = uniprot_to_gene_symbol_dict)

    print( head(enrichment_result))

    enrichment_result_add_gene_symbol <- enrichment_result |>
      mutate( gene_id_list = str_split( geneID, "/") ) |>
      mutate( gene_symbol = purrr::map_chr( gene_id_list,
                                            convertProteinAccToGeneSymbolPartial)) |>
      dplyr::select(-gene_id_list)

  } else {
    enrichment_result_add_gene_symbol <- enrichment_result
  }

  ## generate the output files
  purrr::walk(list_of_comparisons, function(input_comparison) {
    output_file <- paste0( annotation_type, "_table_", input_comparison, ".tab" )

    vroom::vroom_write(enrichment_result_add_gene_symbol |>
                         dplyr::filter( comparison == input_comparison),
                       file=file.path( output_dir,
                                       input_comparison,
                                       output_file ) ) })

  saveRDS( enrichment_objects, file.path(output_dir, "enrichment_objects.rds"))

}

## Draw Cnet
enrichment_objects |>
  mutate( file_name = paste0( paste( input_comparison, go_aspect, names_of_genes_list, sep="-"), ".pdf")) |>
  dplyr::select( input_comparison, file_name, enrichment_objects) |>
  purrr::pmap( \(input_comparison, file_name, enrichment_objects) {
    if(nrow(as.data.frame(enrichment_objects))>0) {
      pdf( file.path( output_dir, input_comparison, file_name))
      print( cnetplot(enrichment_objects) )
      dev.off()  }
  })


enrichment_objects |>
  mutate( file_name = paste0( paste( input_comparison, go_aspect, names_of_genes_list, sep="-"), ".png")) |>
  dplyr::select( input_comparison, file_name, enrichment_objects) |>
  purrr::pmap( \(input_comparison, file_name, enrichment_objects) {
    if(nrow(as.data.frame(enrichment_objects))>0) {
      png( file.path( output_dir, input_comparison, file_name))
      print( cnetplot(enrichment_objects) )
      dev.off()  }
  })


}
