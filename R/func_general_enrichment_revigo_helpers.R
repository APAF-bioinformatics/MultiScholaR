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

