# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
addColumnsToEvidenceTbl <- function(evidence_tbl, phospho_site_prob_col = phospho_sty_probabilities) {
  evidence_tbl_cleaned <- evidence_tbl %>%
      mutate( evidence_id = (row_number() - 1))  %>%
      # dplyr:::select(one_of(c("evidence_id", evidence_col_to_use %>% pull(Columns)))) %>%
      mutate( cleaned_peptide = str_replace_all({{phospho_site_prob_col}}, "[\\(\\)0-9\\.]", ""))

  return( evidence_tbl_cleaned)

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
getMaxProb <- function(phosphopeptide, num_sites=1) {

  pass_thresh <- str_match_all( phosphopeptide,
                                "\\((\\d+\\.*\\d*)\\)") %>%
    .[[1]]  %>%
    .[,2] %>%
    as.numeric

  # Try to preserve the order in which the probability is listed in the peptide,
  # while using sort to find the top 'num_sites'

  if( length(pass_thresh) == 0 ) {
    return( c())
  }

  # Sort from maximum to minimum and then take the first few numbers according to number of sites required
  # Please ensure decreasing is set to TRUE.
  top_site_index <-  sort.int(pass_thresh,
                              index.return=TRUE,
                              decreasing=TRUE)$ix[seq_len(num_sites)]


  pass_thresh[sort( top_site_index)]  %>%
    keep( ~{!is.na(.)}  )

}

#' @export
getMaxProbFutureMap <- function(phosphopeptide, num_sites=1 ) {
  furrr::future_map2( phosphopeptide, num_sites,
                      ~{getMaxProb(.x, .y)}  )
}





## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
getBestPosition <- function(phosphopeptide, num_sites=1 ) {

 if(str_detect(phosphopeptide, "p" ) ) {
   stop("Input phosphopetide string should not have little 'p' as characters.")
 }


pass_thresh <- str_match_all( phosphopeptide,
                "\\((\\d+\\.*\\d*)\\)") %>%
                          .[[1]]  %>%
                          .[,2] %>%
                          as.numeric

 prob_list <- getMaxProb(phosphopeptide, num_sites)

 ## I might need to fix this line as if we have two poisition sharing the same maximum score,
 ## we currently only use the first one as best position
 selected_pos <- which( pass_thresh %in% prob_list)

 little_p_position <- str_replace_all( phosphopeptide,
                                       "\\(\\d+\\.*\\d*\\)", "p" ) %>%
                      str_locate_all("p") %>%
                      .[[1]] %>%
                      .[,1]

  to_adj_pos <- seq_along( little_p_position)

  clean_pos <- little_p_position - to_adj_pos

  return( clean_pos[selected_pos] )

}

#' @export
getBestPositionFutureMap <- function(phosphopeptide, num_sites=1  ) {
  furrr::future_map2( phosphopeptide, num_sites,
                     ~{getBestPosition(.x, .y)}  )
}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Get PTM position string for the modified peptide
#' @description get_pos_string
#' @param peptide_start_position: Start position of peptide relative to the entire protein sequence
#' @param site_relative_position: Position of the modification site relative to the N-terminus of the peptide sequence
#' @return A string. If the position are limited to one unique peptide in the protein, all the phosphosite on that peptide. If the positions are found in multiple repeated peptides in the protein, the phosphosites in each peptide will be contained in round bracket (e.g. (144,148),(170,174),(183,187) ).
#' @export
getPosString <-  function(peptide_start_position, site_relative_position) {

  a <- peptide_start_position
  b <- site_relative_position

  pos_group <- cross2( a, b) %>%
    map_dbl( ~{sum(unlist(.))-1} )

  pos_mat <-  matrix( pos_group,
                      ncol= length(b),
                      nrow= length(a),
                      byrow=FALSE)

  # print( pos_mat)

  if ( length( a) > 1) {
    pos_string <- list()

    for( i in seq_len(nrow(pos_mat)) ) {
      pos_string <- c( pos_string, paste0( "(", paste(  pos_mat[i,], collapse=";"), ")"  ) )
    }
    return( paste( pos_string, collapse="|") )

  } else {

    pos_string <-  paste(  pos_mat[1,], collapse=";")
    return( pos_string )
  }

}





## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Get X-mer string given the primary sequence and the position of the PTM
#' @description Given the bioiString object, the uniprot accession, and the position of the phosphorylation site, return the 15mer sequence with the phosphorylation site at the middle.
#' @export
getXMerString <- function(seq, uniprot_acc, position, padding_length=7 ) {

    start <- position - padding_length
    end <- position + padding_length
    seq_end <- str_length( seq )

    # print( seq_end)

    end_padding <- 0
    start_padding <- 0

    if ( seq_end < end ) {
      end_padding <-   position + padding_length - seq_end
      end <- seq_end
    }

    if ( start < 1) {
      start_padding <- abs( start - 1)
      start <- 1
    }

    start_padding_string <- paste0(rep("_", start_padding), collapse="")
    end_padding_string <- paste0(rep("_", end_padding), collapse="" )

    my_X_mer_partial <- str_sub(seq, start, end )[[1]]

    my_X_mer <- paste0( start_padding_string ,
                         my_X_mer_partial,
                         end_padding_string   )

    return( my_X_mer)

}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
getXMersList <-  function(seq, uniprot_acc,
                          peptide_start_position, site_relative_position, padding_length=7 ) {

  a <- peptide_start_position
  b <- site_relative_position

  pos_group <- cross2( a, b) %>%
    map_dbl( ~{sum(unlist(.))-1} )

  pos_mat <-  matrix( pos_group,
                      ncol= length(b),
                      nrow= length(a),
                      byrow=FALSE)

  # print( uniprot_acc)
  # print( pos_mat)
  # print(as.vector(pos_mat[1,] ) )
  # print( paste( "peptide_start_position = ", paste( peptide_start_position, collapse=";")))
  # print( paste( "site_relative_position = ", paste( site_relative_position, collapse=";")))
  my_Xmers_list <- purrr::map_chr( as.vector(pos_mat[1,] ),
                                ~{getXMerString(seq, uniprot_acc, ., padding_length=padding_length)}) %>%
    paste( collapse=";")

  return( my_Xmers_list )

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title formatPhosphositePosition
#' @param gene_symbol The gene symbol as a string
#' @param position The positions of the phosphorylation sites as a string, multiple sites are separated by the delimiter (e.g. delim parameter). Multiple positions are in the same corresponding order as the residue parameter.
#' @param residue The amino acid residue where the phosphorylation is at, multiple residues are separate by the delimiter (e.g. delim paramter). Multiple residues are in the same corresponding order as the position parameter.
#' @return A string with the gene symbol, followed by a semi-colon, and the amino acid residue and the position, multiple residue and positions separated by commas (e.g. YFG1;S199,T203 )
#' @export

formatPhosphositePosition <- function( gene_symbol, position, residue, delim=";") {
  position_cln <- str_split( position, "\\|" )[[1]][1] |>
    str_replace_all( "\\(|\\)", "") |>
    str_split(delim) |>
    unlist()

  residue_cln <- residue |>
    str_split( ";") |>
    unlist()

  list_of_positions <- purrr::map2_chr(residue_cln,  position_cln, \(x,y){paste0(x, y)})

  formatted_positions <- paste( gene_symbol, paste0(list_of_positions, collapse=","), sep=";" )

  return(formatted_positions)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
removePeptidesWithoutAbundances <- function(evidence_tbl_cleaned, col_pattern) {

  sites_to_accept <- evidence_tbl_cleaned %>%
    mutate( across( matches(col_pattern, perl=TRUE), ~.==0 )) %>%
    dplyr::filter( !if_all( matches(col_pattern, perl=TRUE), ~. ==TRUE )) %>%
    dplyr::select( evidence_id)

  ## Removing entries where all the "Reporter intensity corrected" rows are zero
  evidence_tbl_filtered <- evidence_tbl_cleaned %>%
    inner_join(sites_to_accept, by=c("evidence_id" = "evidence_id") )

  return( evidence_tbl_filtered )
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
filterPeptideAndExtractProbabilities <- function(evidence_tbl_cleaned, accession_gene_name_tbl,
                                                 col_pattern="corrected",
                                                 accession_col = leading_proteins,
                                                 phospho_site_prob_col = phospho_sty_probabilities,
                                                 num_phospho_site_col = phospho_sty) {


  sites_probability_tbl <- evidence_tbl_cleaned %>%
    ## Must be a phosphopeptide (at least one site)
    dplyr::filter( {{num_phospho_site_col}} >=1) %>%
    ## Remove REV_ and CON_
    dplyr::filter( !str_detect( {{accession_col}}, "^REV__|^CON__" )  ) %>%
    dplyr::mutate( best_phos_prob = getMaxProbFutureMap({{phospho_site_prob_col}},
                                                        {{num_phospho_site_col}})) %>%
    dplyr::filter( map_lgl(best_phos_prob, ~{length(.) > 0} )) %>%
    dplyr::mutate( best_phos_pos = getBestPositionFutureMap({{phospho_site_prob_col}},
                                                            {{num_phospho_site_col}})) %>%
    ## Avoid cases where there are multiple positions having the same top scores
    dplyr::filter( map2_lgl(best_phos_prob, best_phos_pos, ~{length(.x) == length(.y)} )) %>%
    left_join( accession_gene_name_tbl, by="evidence_id") %>%
    dplyr::select( one_of( c( "best_phos_prob", "best_phos_pos",


                              colnames(accession_gene_name_tbl),
                              colnames( evidence_tbl_cleaned))),
                   {{phospho_site_prob_col}},
                   {{num_phospho_site_col}}
                   )


  number_of_rows_without_uniprot_acc <- sites_probability_tbl %>%
    dplyr::filter( is.na(uniprot_acc) ) %>%
    nrow()


  if( length(number_of_rows_without_uniprot_acc) > 0) {
    warnings( paste("There are", number_of_rows_without_uniprot_acc, "proteins do not have sequence information in FASTA file. Removing them from analysis"))
  }

  sites_probability_filt <- sites_probability_tbl %>%
    dplyr::filter( !is.na(uniprot_acc) )


  return(sites_probability_filt )

}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
addPeptideStartAndEnd <- function(sites_probability_tbl, aa_seq_tbl ) {
  peptide_start_and_end <- sites_probability_tbl %>%
  left_join( aa_seq_tbl %>% dplyr::select(uniprot_acc, seq), by=c("uniprot_acc" = "uniprot_acc")) %>%
  mutate( peptide_location =  str_locate_all(seq, cleaned_peptide)) %>%
  mutate( pep_start  = map_chr ( peptide_location, ~paste(.[,"start"], collapse="|" )     )   ) %>%
  mutate( pep_end  = map_chr ( peptide_location, ~paste(.[,"end"], collapse="|" )     )   )

  return( peptide_start_and_end)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
addPhosphositesPositionsString <- function(peptide_start_and_end ) {

  phosphosite_pos_string_tbl <- peptide_start_and_end %>%
    mutate( best_phos_pos_string = map_chr(best_phos_pos, ~paste(., collapse=";") )) %>%
    mutate( temp_check_pos =  map2(peptide_location, best_phos_pos, ~{cross2( .x[,"start"] , .y) } )   ) %>%
    mutate( check_pos =  purrr::map(temp_check_pos, ~{ map_dbl(., function(x){sum(unlist(x)) -1} )}   ) ) %>%
    mutate( protein_site_positions = map2_chr(peptide_location, best_phos_pos, ~{getPosString(.x[, "start"] , .y) } )  )  %>%
    mutate( best_phos_prob_string = map_chr(best_phos_prob, ~paste(., collapse=";") ))

  return( phosphosite_pos_string_tbl)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
addXMerStrings <- function (phosphosite_pos_string_tbl, padding_length=7) {


  my_get_X_mers_list <-  function(uniprot_acc, peptide_location, best_phos_pos, seq) {
    getXMersList(seq, uniprot_acc, peptide_location, best_phos_pos, padding_length=padding_length)
  }

  get_15_mer_tbl <-   phosphosite_pos_string_tbl %>%
    mutate( phos_15mer_seq = furrr::future_pmap_chr(  list( uniprot_acc = uniprot_acc,
                                                            peptide_location = peptide_location,
                                                            best_phos_pos = best_phos_pos,
                                                            seq = seq),
                                                     my_get_X_mers_list )  ) %>%
    distinct()

  return( get_15_mer_tbl)

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
temp <- function ( myinput = `lotr`) {
  # print(as_string( {{myinput}}))

  my_tab <- data.frame(lort=rep(10,10))

  my_tab %>%
    dplyr::select( {{myinput}})

   my_string <- as_name( enquo(myinput))

   print(typeof(my_string))
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
filterByScoreAndGetSimilarPeptides <- function(get_15_mer_tbl, site_prob_threshold, secondary_site_prob_threshold = 0.5, num_phospho_site_col = phospho_sty ) {

  ## Find peptide in which at least one phosphosite has one position >= site probability threshold
  all_peptide_and_sites_pass_filter <- get_15_mer_tbl %>%
    dplyr::filter(  map2_lgl ( best_phos_prob,  {{num_phospho_site_col}},
                              function(x, y) {  return( length(which( x >= site_prob_threshold)) >=  y) }  ) )  %>%
    dplyr::select( uniprot_acc, protein_site_positions) %>%
    distinct() %>%
    dplyr::mutate( protein_site_positions = str_split(protein_site_positions, "[;\\|]")  )   %>%
    unnest( protein_site_positions) %>%
    dplyr::mutate( protein_site_positions = as.integer( str_replace_all( protein_site_positions, "\\(|\\)", "")  )) %>%
    distinct() %>%
    arrange( uniprot_acc, protein_site_positions) %>%
    group_by(uniprot_acc) %>%
    nest(high_quality_sites = protein_site_positions ) %>%
    ungroup()  %>%
    mutate( high_quality_sites = purrr::map( high_quality_sites, ~{ .$protein_site_positions }))


  ## Find peptide where all major sites probability is > secondary_site_prob_threshold
  peptide_and_pos_pass_filt <- get_15_mer_tbl %>%
    dplyr::filter(   map2_lgl ( best_phos_prob,   {{num_phospho_site_col}},
                              function(x, y) {  return(  length(which( x >  secondary_site_prob_threshold)) >= y  ) }  ) )  %>%
    dplyr::select( uniprot_acc, protein_site_positions) %>%
    distinct() %>%
    dplyr::mutate( protein_site_positions = str_split(protein_site_positions, "[;\\|]")  ) %>%
    dplyr::mutate( protein_site_positions = purrr::map( protein_site_positions, ~{ str_replace_all( ., "\\(|\\)", "") %>% purrr::map_int(as.integer) }   )   ) %>%
    distinct()


  all_filtered_peptide <- peptide_and_pos_pass_filt %>%
    dplyr::inner_join( all_peptide_and_sites_pass_filter, by = "uniprot_acc") %>%
    dplyr::filter( map2_lgl(protein_site_positions,
                            high_quality_sites,
                            function(to_check, hq_sites) { length(which(to_check %in%  hq_sites) ) == length(to_check)   }))   %>%
    dplyr::mutate( protein_site_positions =  map_chr( protein_site_positions,
                                                      ~paste(., collapse=";")))

  ## Check that all sites in peptide has been found at least one in set "all_peptide_and_sites_pass_filter"
  # Keep peptide if there is another peptide that has >0.75 at all of these sites and secondary site with highest probability in the same position.
  get_15_mer_tbl_filt <- get_15_mer_tbl %>%
    dplyr::mutate( temp_protein_site_positions = purrr::map_chr( protein_site_positions,
                                                            ~{ str_replace_all(., "\\(|\\)", "") %>%
                                                               str_replace_all( "\\|", ";" ) }  )) %>%
    dplyr::inner_join( all_filtered_peptide, by =c( "uniprot_acc" = "uniprot_acc",
                                                         "temp_protein_site_positions" = "protein_site_positions")) %>%
    dplyr::select(-temp_protein_site_positions)

  return( get_15_mer_tbl_filt)

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
allPhosphositesPivotLonger <- function(get_15_mer_tbl,
                                       additional_cols = c("experiment"),
                                       col_pattern = "Reporter intensity corrected",
                                       pattern_suffix = "_\\d+",
                                       extract_patt_suffix="_(\\d+)",
                                       phospho_site_prob_col = phospho_sty_probabilities,
                                       num_phospho_site_col = phospho_sty
                                           ) {

  usual_columns <- c( "evidence_id", "uniprot_acc", "gene_name", "sequence", # "gene_names",
                      "protein_site_positions", "phos_15mer_seq", as_name(enquo( phospho_site_prob_col)), as_name(enquo( num_phospho_site_col)) )

  cols_to_use <- usual_columns

  if ( !is.na( additional_cols) & additional_cols != "") {
    cols_to_use <- c( usual_columns, additional_cols)
  }

  all_sites_long <- get_15_mer_tbl %>%
    dplyr::select( {{cols_to_use}},
                   matches( paste0(tolower(col_pattern), pattern_suffix ), perl = TRUE) ) %>%
    pivot_longer( cols = matches(c( paste0(tolower(col_pattern), pattern_suffix )), perl = TRUE) ,
                  names_to = "replicate",
                  values_to = "value")

  # print(head(all_sites_long))
  # print( col_pattern)

  if ( extract_patt_suffix != "") {
    all_sites_long <- all_sites_long %>%
      dplyr::mutate( replicate = str_replace(replicate,
                                             paste0(tolower(col_pattern), extract_patt_suffix ),
                                             "\\1") %>%
                       map_int(as.integer) )
  }


  return(all_sites_long)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
groupParalogPeptides <- function(all_sites_long,
                                 additional_cols = c("experiment"),
                                 phospho_site_prob_col = phospho_sty_probabilities,
                                 num_phospho_site_col = phospho_sty) {

  grouping_variables <- c( "evidence_id", "replicate",
                           "value", "sequence",
                           as_name(enquo( phospho_site_prob_col)),
                           as_name(enquo( num_phospho_site_col)))

  if ( !is.na( additional_cols) & additional_cols != "" ) {
    grouping_variables <- c( grouping_variables, additional_cols)
  }

  final_select_var <- c( "evidence_id", "uniprot_acc", "gene_names",
                         "protein_site_positions", "phos_15mer_seq",
                         "replicate", "value",
                           as_name(enquo( phospho_site_prob_col)),
                           as_name(enquo( num_phospho_site_col)) )

  if ( !is.na( additional_cols)  & additional_cols != "" ) {
    final_select_var <- c( final_select_var,
                           additional_cols )
  }

  ## Group homolog gene names, uniprot_acc, site posiitons
  paralog_sites_long <- all_sites_long %>%
    group_by( across({{ grouping_variables }}) ) %>%
    summarise( gene_names = paste(gene_name, collapse = ":"),
               uniprot_acc = paste(uniprot_acc, collapse=":"),
               protein_site_positions = paste( protein_site_positions, collapse=":"),
               phos_15mer_seq = paste( phos_15mer_seq, collapse=":"),
    ) %>%
    ungroup() %>%
    dplyr::select( one_of( final_select_var ) )

  return(paralog_sites_long)

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
allPhosphositesPivotWider <- function(all_phos_sites_long_tbl,
                                      additional_cols = c("experiment"),
                                      phospho_site_prob_col = phospho_sty_probabilities,
                                      num_phospho_site_col = phospho_sty ) {
 cols_to_use <- "replicate"

  temp_tbl <- all_phos_sites_long_tbl

  if ( !is.na( additional_cols) & additional_cols != "" ) {
    cols_to_use <-c(  "replicate", additional_cols)

    temp_tbl <- all_phos_sites_long_tbl %>%
      mutate_at( additional_cols, toupper )
  }

  all_phos_sites_wide_tbl <-  temp_tbl %>%
    pivot_wider( id_cols = c( evidence_id, uniprot_acc, gene_names,
                              protein_site_positions, phos_15mer_seq,
                           as_name(enquo( phospho_site_prob_col)),
                           as_name(enquo( num_phospho_site_col))),
                 names_from = all_of(cols_to_use) ,
                 values_from = value) %>%
    arrange( uniprot_acc, gene_names, protein_site_positions,
             phos_15mer_seq, evidence_id)

  return( all_phos_sites_wide_tbl)
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
uniquePhosphositesSummariseLongList <- function(all_phos_sites_long_tbl,
                                                additional_cols = c("experiment") ) {

  ## Summarise the input table with a summarisation function
  group_summary <- function( input_tbl, additional_cols, method=mean ) {

    usual_columns <- c( "uniprot_acc",  "gene_names", "protein_site_positions",  "phos_15mer_seq", "replicate" )

    cols_to_use <- usual_columns

    if ( !is.na( additional_cols) & additional_cols != "" ) {
      cols_to_use <- c( usual_columns, additional_cols)
    }

    temp_tbl <- input_tbl %>%
      group_by( across({{ cols_to_use }}) ) %>%
      # uniprot_acc, gene_names, protein_site_positions, phos_15mer_seq, experiment, replicate
      summarise( value =  method( value) ,
                 maxquant_row_ids= paste0(evidence_id, collapse=";") ) %>%
      ungroup  %>%
      mutate( replicate = toupper(replicate))


    output_tbl <- temp_tbl
    if ( !is.na( additional_cols) & additional_cols != "" ) {

      output_tbl <- temp_tbl %>%
      mutate_at(  additional_cols, toupper)

    }

    return( output_tbl)
  }

  summary_funcs <- list( mean=mean, median=median, sum=sum)

  summarised_long_tbl_list <- purrr::map( summary_funcs, ~group_summary(all_phos_sites_long_tbl, additional_cols, .))

  return( summarised_long_tbl_list)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
uniquePhosphositesSummariseWideList <- function(summarised_long_tbl_list,
                                                additional_cols=c("experiment")) {


  cols_to_use <- c("replicate")

  summarised_wide_tbl_list_edited <- NA
  if ( !is.na( additional_cols) & additional_cols != "" ) {
    cols_to_use <- c( "replicate", additional_cols)

    experiment_col <- additional_cols[[1]]

    summarised_wide_tbl_list_edited <- purrr::map( summarised_long_tbl_list,
                                                   function(input_table){ output_table <- input_table %>%
      ## When there is additional cols use the first additional cols and add it to the maxquant_row_ids
      mutate( maxquant_row_ids = paste0( paste(!!rlang::sym(experiment_col ), sep="_") , "(", maxquant_row_ids, ")") ) %>%
      pivot_wider( id_cols = c("uniprot_acc", "gene_names", "protein_site_positions", "phos_15mer_seq", "maxquant_row_ids"),
                   names_from = all_of(cols_to_use),
                   values_from=c("value") )%>%
      unite( "sites_id", uniprot_acc, gene_names, protein_site_positions, phos_15mer_seq, sep="!" )

    return(output_table)}  )
  } else {

    summarised_wide_tbl_list_edited <- purrr::map( summarised_long_tbl_list, function(input_table){  output_table <- input_table %>%
      ## When there is additional cols use the first additional cols and add it to the maxquant_row_ids
      pivot_wider( id_cols = c("uniprot_acc", "gene_names", "protein_site_positions", "phos_15mer_seq", "maxquant_row_ids"),
                   names_from = all_of(cols_to_use),
                   values_from=c("value") )%>%
      unite( "sites_id", uniprot_acc, gene_names, protein_site_positions, phos_15mer_seq, sep="!" )

    return(output_table)}  )
    }



    ## Summarize MaxQuant evidence IDs from different multiplex experiment
    clean_maxquant_ids <- function(input_tab ) {
      maxquant_ids_tbl <- input_tab  %>%
        group_by( sites_id) %>%
        summarise( maxquant_row_ids = paste(maxquant_row_ids, collapse=",")  ) %>%
        ungroup()


      values_tbl <- input_tab %>%
        dplyr::select(-maxquant_row_ids) %>%
        group_by( sites_id) %>%
        summarise_all( ~sum(., na.rm=TRUE)) %>%
        ungroup()

      output_tab <- values_tbl %>%
        left_join( maxquant_ids_tbl, by="sites_id") %>%
        relocate( maxquant_row_ids, .after="sites_id")

      colnames( output_tab) <- colnames( output_tab ) %>%
        toupper( ) %>%
        str_replace_all( "SITES_ID", "sites_id")  %>%
        str_replace_all( "MAXQUANT_ROW_IDS", "maxquant_row_ids")

      return( output_tab)

    }

    summarised_wide_tbl_cln_list <- purrr::map( summarised_wide_tbl_list_edited, clean_maxquant_ids)

    return( summarised_wide_tbl_cln_list)


  }
# The values for sum is way too large, so I think it is going to be median or mean



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
processMultisiteEvidence <- function(fasta_file,
                                     evidence_tbl,
                                     accession_col = leading_proteins,
                                     group_id,
                                     additional_cols = c(experiment),
                                     col_pattern="corrected",
                                     extract_pattern = "Reporter intensity corrected",
                                     col_to = "",
                                     site_prob_threshold = 0.75,
                                     columns_to_use = NA) {
  ## Read fasta file

  print( "Step 1: Reading the fasta file.")
  aa_seq_tbl <- parseFastaFile(fasta_file)

  ## Add the row id column and create a column containing the cleaned  peptide
  print("Step 2: Get row ID and get cleaned peptide sequence.")
  evidence_tbl_cleaned <- addColumnsToEvidenceTbl(evidence_tbl, evidence_col_to_use )

  ## Get best accession per entry, work out peptides mapped to multiple genes
  print("Step 3: Use decision tree to get best accession per phosphosite evidence entry")
  accession_gene_name_tbl <- chooseBestAccession(evidence_tbl_cleaned,
                                                 aa_seq_tbl,
                                                 {{accession_col}},
                                                 {{group_id}})

  ## Remove peptides without abundance values at all
  print("Step 4: Remove peptides without abundance values at all")
  evidence_tbl_filt <- removePeptidesWithoutAbundances(evidence_tbl_cleaned, col_pattern)

  ## For all the multi-phosphosites peptide extract their intensity, filter peptide with no intensity across all samples, extract site probabilities
  print("Step 5: Filter peptides with no intensity across all samples, extract intensity data, extract sites")
  sites_probability_tbl <- filterPeptideAndExtractProbabilities (evidence_tbl_filt,
                                                                 accession_gene_name_tbl,
                                                                 col_pattern,
                                                                 accession_col = {{accession_col}} )

  ## Get the peptide start and end position for each peptide
  print("Step 6: Add peptide start and end position")
  peptide_start_and_end <- addPeptideStartAndEnd(sites_probability_tbl , aa_seq_tbl )


  ## Get the phosphosites position string
  print("Step 7: Add string listing the positions of phosphosites")
  phosphosite_pos_string_tbl <- addPhosphositesPositionsString(peptide_start_and_end )

  ## Get the string listing all the 15-mer sequences, each sequence has the phosphorylation site in the middle
  print("Step 8: Add string listing all 15-mer sequences, each sequence has phosphosite in the center")
  get_15_mer_tbl <- addXMerStrings(phosphosite_pos_string_tbl, padding_length=7)

  ## Get peptides with at least one phosphosite over threshold. Find all peptides with same sites as another peptides that contained at least one phosphosies >= threshold.
  print("Step 9: Get high conf. peptides (e.g. phosphosites >= threshold). Get peptide W/ same sites as high conf. peptide.")
  get_15_mer_tbl_filt <- filterByScoreAndGetSimilarPeptides(get_15_mer_tbl, site_prob_threshold)

  ## Pivot the phosphosites to a longer table
  print("Step 10: Pivot phosphosites/phosphopeptide table to long format")
  all_phos_sites_long_tbl <- allPhosphositesPivotLonger(get_15_mer_tbl_filt,
                                                        additional_cols ,
                                                        col_pattern )

  ## Group peptides from paralog proteins
  print("Step 11: Group peptides from paralog proteins ")
  paralog_sites_long <- groupParalogPeptides (all_phos_sites_long_tbl, additional_cols )

  ## Pivot the phosphosites data to a wide format
  print("Step 12: Pivot phosphosites/phosphopeptide table to wide format")
  all_phos_sites_wide_tbl <- allPhosphositesPivotWider(paralog_sites_long,
                                                       additional_cols )



  ## Summarise the abundance values for each unique phosphosites (mean, median, sum), return table in long format
  print("Step 13: Summarise abundance values for each unique phosphosites, long format")
  summarised_long_tbl_list <- uniquePhosphositesSummariseLongList(paralog_sites_long,
                                                                  additional_cols )

  ## Summarise the abundance values for each unique phosphosites (mean, median, sum), return table in wide format
  print("Step 14: Summarise abundance values for each unique phosphosites, wide format")
  summarised_wide_tbl_list <- uniquePhosphositesSummariseWideList(summarised_long_tbl_list,
                                                                  additional_cols)

  ## The values for sum is way too large, so I think it is going to be median or mean


  return( list( summarised_wide_list = summarised_wide_tbl_list,
                summarised_long_list = summarised_long_tbl_list,
                all_phos_sites_wide  = all_phos_sites_wide_tbl,
                all_phos_sites_long  = all_phos_sites_long_tbl ))

}

###--------------------------------------------------------------------------------------------------------------------------------

#' @export
#' @description Given an input table with sites_id and uniprot_acc columns, work out the ranking of the uniprot_acc within the sites_id
#' @param input_table, an input table with sites_id and uniprot_acc columns
#' @param uniprot_acc, name of Uniprot accession column tidyverse style column name input
#' @param sites_id, name of sites_id column, in tidyverse style column name
getUniprotAccRankFromSitesId <- function( input_table, uniprot_acc, sites_id) {
  input_table %>%
    dplyr::mutate( uniprot_acc_split = str_split({{uniprot_acc}}, ":" ) %>% purrr::map_chr(1) ) %>%
    dplyr::mutate( uniprot_list = str_split( {{sites_id}}, "!") %>% purrr::map_chr(1) %>% str_split( ":")) %>%
    dplyr::mutate( gene_list_position = purrr::map2_int( uniprot_acc_split, uniprot_list, ~{ which(.x == .y)[1]}))  %>%
    relocate(uniprot_acc_split, .after=lazyeval::as_name(enquo(sites_id)) ) %>%
    relocate(uniprot_list, .after="uniprot_acc_split" ) %>%
    relocate(gene_list_position, .after="uniprot_list" )
}

###--------------------------------------------------------------------------------------------------------------------------------
