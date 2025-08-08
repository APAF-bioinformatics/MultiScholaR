#' Add Helper Columns to Evidence Table
#'
#' @description
#' Adds an `evidence_id` and a `cleaned_peptide` column to the evidence table.
#'
#' @param evidence_tbl A data frame or tibble, typically the `evidence.txt` output from MaxQuant.
#' @param phospho_site_prob_col An unquoted column name specifying the column that
#'   contains phosphopeptide sequences with probability scores (e.g., `phospho_sty_probabilities`).
#'   This column is used to generate the `cleaned_peptide` by removing probability scores.
#'
#' @return The input tibble with two new columns: `evidence_id` (a unique row identifier)
#'   and `cleaned_peptide` (the peptide sequence with probability scores removed).
#'
#' @importFrom dplyr mutate row_number
#' @importFrom stringr str_replace_all
#' @export
addColumnsToEvidenceTbl <- function(evidence_tbl, phospho_site_prob_col = phospho_sty_probabilities) {
  evidence_tbl_cleaned <- evidence_tbl %>%
      mutate( evidence_id = (row_number() - 1))  %>%
      # dplyr:::select(one_of(c("evidence_id", evidence_col_to_use %>% pull(Columns)))) %>%
      mutate( cleaned_peptide = str_replace_all({{phospho_site_prob_col}}, "[\\(\\)0-9\\.]", ""))

  return( evidence_tbl_cleaned)

}



#' Extract Highest Phosphorylation Site Probabilities
#'
#' @description
#' From a phosphopeptide string containing modification probabilities in parentheses,
#' this function extracts and returns the highest probability values.
#'
#' @param phosphopeptide A character string for a single phosphopeptide, where
#'   probabilities are enclosed in parentheses, e.g., `_VSDSG(0.99)YSSG(0.01)SLSGR_`.
#' @param num_sites The number of top probability scores to return. For example, if
#'   a peptide has three phosphosites but `num_sites` is 2, only the top two
#'   probabilities will be returned. Default is 1.
#'
#' @return A numeric vector containing the top `num_sites` probability scores,
#'   sorted in the order they appeared in the original peptide string. Returns an
#'   empty vector if no probabilities are found.
#'
#' @importFrom stringr str_match_all
#' @importFrom purrr keep
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

#' Map `getMaxProb` over Lists using `furrr`
#'
#' @description
#' A wrapper that applies the `getMaxProb` function over parallel lists of
#' phosphopeptides and their corresponding number of sites to consider.
#'
#' @param phosphopeptide A list or vector of phosphopeptide character strings.
#' @param num_sites A list or vector of integers specifying the number of top
#'   sites to extract for each corresponding phosphopeptide.
#'
#' @return A list of numeric vectors, where each vector contains the top
#'   probability scores for the corresponding input phosphopeptide.
#'
#' @importFrom furrr future_map2
#' @export
getMaxProbFutureMap <- function(phosphopeptide, num_sites=1 ) {
  furrr::future_map2( phosphopeptide, num_sites,
                      ~{getMaxProb(.x, .y)}  )
}





#' Get Best Phosphorylation Site Position
#'
#' @description
#' From a phosphopeptide string, this function identifies the positions of the
#' modification sites corresponding to the highest probability scores.
#'
#' @details
#' The function first identifies all modification sites (marked by parentheses)
#' and their probabilities. It then determines the highest probabilities using
#' `getMaxProb` and finds the relative positions of these top sites within the
#' peptide sequence. The positions are adjusted to account for the characters
#' removed (the probability scores themselves).
#'
#' @param phosphopeptide A character string for a single phosphopeptide, e.g.,
#'   `_VSDSG(0.99)YSSG(0.01)SLSGR_`.
#' @param num_sites The number of top sites (and their positions) to return.
#'   Default is 1.
#'
#' @return A numeric vector of the relative positions of the best-scoring
#'   phosphorylation sites within the peptide.
#'
#' @importFrom stringr str_detect str_match_all str_replace_all str_locate_all
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

#' Map `getBestPosition` over Lists using `furrr`
#'
#' @description
#' A wrapper that applies the `getBestPosition` function over parallel lists of
#' phosphopeptides and their corresponding number of sites.
#'
#' @param phosphopeptide A list or vector of phosphopeptide character strings.
#' @param num_sites A list or vector of integers specifying the number of top
#'   sites for which to find positions.
#'
#' @return A list of numeric vectors, where each vector contains the relative
#'   positions of the best-scoring sites for the corresponding phosphopeptide.
#'
#' @importFrom furrr future_map2
#' @export
getBestPositionFutureMap <- function(phosphopeptide, num_sites=1  ) {
  furrr::future_map2( phosphopeptide, num_sites,
                     ~{getBestPosition(.x, .y)}  )
}




#' Calculate Absolute PTM Positions in Protein Sequence
#'
#' @description
#' Calculates the absolute positions of modification sites within the full protein
#' sequence, given the start position(s) of the peptide and the relative positions
#' of the sites within the peptide.
#'
#' @param peptide_start_position A numeric vector of start positions of the peptide
#'   within the protein. Can contain multiple values if the peptide sequence is repeated.
#' @param site_relative_position A numeric vector of the positions of modification
#'   sites relative to the start of the peptide (1-indexed).
#'
#' @return A character string of absolute site positions.
#'   - If the peptide occurs once, positions are semi-colon-separated (e.g., "144;148").
#'   - If the peptide is repeated, positions for each peptide instance are grouped in
#'     parentheses and pipe-separated (e.g., "(144;148)|(170;174)").
#'
#' @importFrom purrr cross2 map_dbl
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





#' Extract X-mer Sequence Motif
#'
#' @description
#' Extracts a sequence motif (X-mer) of a specified length centered around a
#' modification site. If the site is too close to the N or C terminus, the
#' sequence is padded with underscores.
#'
#' @param seq The full protein sequence as a character string.
#' @param uniprot_acc The UniProt accession of the protein (currently unused but
#'   kept for signature consistency).
#' @param position The absolute position of the modification site in the protein.
#' @param padding_length The number of amino acids to include on each side of the
#'   central modification site. The total length of the X-mer will be `2 * padding_length + 1`.
#'   Default is 7 (for a 15-mer).
#'
#' @return A character string representing the X-mer sequence, padded with `_`
#'   if necessary.
#'
#' @importFrom stringr str_length str_sub
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




#' Get List of X-mer Motifs for a Peptide
#'
#' @description
#' A wrapper function that calculates absolute PTM positions and then extracts
#' the corresponding X-mer sequence motifs for all modification sites within a peptide.
#'
#' @param seq The full protein sequence.
#' @param uniprot_acc The UniProt accession of the protein.
#' @param peptide_start_position A numeric vector of peptide start positions.
#' @param site_relative_position A numeric vector of PTM positions relative to the peptide start.
#' @param padding_length The padding length to use for the X-mer extraction (see `getXMerString`).
#'   Default is 7.
#'
#' @return A single character string with all X-mer sequences for the peptide,
#'   separated by semicolons.
#'
#' @importFrom purrr cross2 map_dbl map_chr
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



#' Format Phosphosite Position String
#'
#' @description
#' Creates a standardized string representing phosphosites, combining the gene
#' symbol with the residue and position of each site.
#'
#' @param gene_symbol The gene symbol as a string.
#' @param position A character string of site positions, potentially with complex
#'   formatting (e.g., "(1986)|(1998)"). The function processes the first
#'   pipe-separated element.
#' @param residue A character string of modified amino acid residues (e.g., "S;Y"),
#'   separated by semicolons.
#' @param delim The delimiter used within the `position` string. Default is ";".
#'
#' @return A formatted string, e.g., "YFG1;S199,T203".
#'
#' @importFrom stringr str_split str_replace_all
#' @importFrom purrr map2_chr
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

# # formatPhosphositePosition( "MFG1", "(1986)|(1998)|(2010)	", "S")
# # formatPhosphositePosition( "MFG1", "1082;1087", "S;Y")

#' Remove Peptides with No Abundance Across All Samples
#'
#' @description
#' Filters an evidence table to remove rows (peptides) where all abundance
#' columns (matching a specified pattern) are zero.
#'
#' @param evidence_tbl_cleaned The input evidence table, which must have an `evidence_id` column.
#' @param col_pattern A regex pattern to identify the abundance columns to check.
#'
#' @return A filtered tibble containing only rows with at least one non-zero
#'   abundance value in the specified columns.
#'
#' @importFrom dplyr mutate across matches filter if_all select inner_join
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



#' Filter and Process Phosphopeptides from Evidence Table
#'
#' @description
#' This function performs initial filtering and processing of phosphopeptide
#' evidence. It filters for actual phosphopeptides, removes contaminants,
#' extracts localization probabilities and positions, and joins with gene/protein
#' annotation information.
#'
#' @param evidence_tbl_cleaned The input evidence table.
#' @param accession_gene_name_tbl A table mapping evidence IDs to UniProt accessions
#'   and gene names.
#' @param col_pattern (Currently unused) A regex pattern for columns.
#' @param accession_col The unquoted column name for protein accessions (e.g., `leading_proteins`).
#' @param phospho_site_prob_col The unquoted column name for phosphosite probabilities.
#' @param num_phospho_site_col The unquoted column name for the number of phosphosites.
#'
#' @return A tibble with processed and filtered phosphopeptide data, including
#'   extracted probabilities (`best_phos_prob`) and positions (`best_phos_pos`).
#'
#' @importFrom dplyr filter mutate select one_of left_join
#' @importFrom stringr str_detect
#' @importFrom purrr map_lgl map2_lgl
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




#' Add Peptide Start and End Positions
#'
#' @description
#' Maps peptide sequences to their corresponding full protein sequences to find
#' their start and end locations.
#'
#' @param sites_probability_tbl A tibble containing processed phosphopeptide data,
#'   including `uniprot_acc` and `cleaned_peptide` columns.
#' @param aa_seq_tbl A tibble containing full protein sequences, with `uniprot_acc`
#'   and `seq` columns.
#'
#' @return The input tibble with new columns: `peptide_location` (a list of
#'   stringr locate matrices), `pep_start` (character string of start positions),
#'   and `pep_end` (character string of end positions).
#'
#' @importFrom dplyr left_join select mutate
#' @importFrom stringr str_locate_all
#' @importFrom purrr map_chr
#' @export
addPeptideStartAndEnd <- function(sites_probability_tbl, aa_seq_tbl ) {
  peptide_start_and_end <- sites_probability_tbl %>%
  left_join( aa_seq_tbl %>% dplyr::select(uniprot_acc, seq), by=c("uniprot_acc" = "uniprot_acc")) %>%
  mutate( peptide_location =  str_locate_all(seq, cleaned_peptide)) %>%
  mutate( pep_start  = map_chr ( peptide_location, ~paste(.[,"start"], collapse="|" )     )   ) %>%
  mutate( pep_end  = map_chr ( peptide_location, ~paste(.[,"end"], collapse="|" )     )   )

  return( peptide_start_and_end)

}


#' Add Formatted Phosphosite Position Strings
#'
#' @description
#' Calculates and adds several formatted string columns related to phosphosite
#' positions and probabilities to the data table.
#'
#' @param peptide_start_and_end A tibble from `addPeptideStartAndEnd`.
#'
#' @return The input tibble with new columns: `best_phos_pos_string`,
#'   `protein_site_positions` (absolute positions), and `best_phos_prob_string`.
#'
#' @importFrom dplyr mutate
#' @importFrom purrr map_chr map2 map_dbl cross2
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


#' Add X-mer Sequence Motif Strings
#'
#' @description
#' For each phosphosite, this function extracts the X-mer sequence motif
#' (e.g., a 15-mer) centered on the site.
#'
#' @param phosphosite_pos_string_tbl A tibble from `addPhosphositesPositionsString`.
#' @param padding_length The padding length for the X-mer (see `getXMerString`).
#'   Default is 7 for a 15-mer.
#'
#' @return The input tibble with a new column `phos_15mer_seq` containing a
#'   semicolon-separated string of all X-mer motifs for the peptide.
#'
#' @importFrom dplyr mutate distinct
#' @importFrom furrr future_pmap_chr
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



#' Filter Phosphosites by Probability Score and Similarity
#'
#' @description
#' Performs a two-step filtering process to retain high-confidence phosphosites.
#' First, it identifies all sites on proteins that have at least one high-probability
#' site. Second, it filters the full peptide list to keep only those peptides where
#' all identified sites meet a secondary probability threshold and belong to the
#' set of high-confidence sites.
#'
#' @param get_15_mer_tbl The processed phosphosite table.
#' @param site_prob_threshold The primary probability threshold (e.g., 0.75) to
#'   define a "high-quality" site.
#' @param secondary_site_prob_threshold A more lenient threshold (e.g., 0.5) that
#'   all sites in a kept peptide must pass.
#' @param num_phospho_site_col An unquoted column name for the number of phosphosites.
#'
#' @return A filtered tibble containing only high-confidence phosphopeptide evidence.
#'
#' @importFrom dplyr filter select distinct mutate group_by nest ungroup inner_join arrange
#' @importFrom tidyr unnest
#' @importFrom purrr map map2_lgl map_chr map_int
#' @importFrom stringr str_split str_replace_all
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



#' Pivot Phosphosite Table to Long Format
#'
#' @description
#' Converts a wide-format phosphosite table (with samples as columns) into a
#' long format, making it suitable for analysis and plotting with tools like ggplot2.
#'
#' @param get_15_mer_tbl The wide-format input tibble.
#' @param additional_cols A character vector of additional columns to preserve
#'   during pivoting. Default is `c("experiment")`.
#' @param col_pattern A regex pattern to identify the sample/abundance columns.
#' @param pattern_suffix A regex pattern for the suffix of sample columns.
#' @param extract_patt_suffix A regex pattern with a capture group to extract
#'   the replicate number from the column name.
#' @param phospho_site_prob_col Unquoted column name for phosphosite probabilities.
#' @param num_phospho_site_col Unquoted column name for the number of phosphosites.
#'
#' @return A long-format tibble with columns `replicate` and `value`.
#'
#' @importFrom dplyr select matches mutate
#' @importFrom tidyr pivot_longer
#' @importFrom rlang as_name enquo
#' @importFrom stringr str_replace
#' @importFrom purrr map_int
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


#' Group Peptides from Paralogous Proteins
#'
#' @description
#' Collapses rows that represent peptides from paralogous proteins. It groups
#' by evidence ID and other identifying columns, then pastes together the
#' differing gene names, UniProt accessions, and site positions.
#'
#' @param all_sites_long A long-format phosphosite table.
#' @param additional_cols Additional columns to include in the grouping.
#' @param phospho_site_prob_col Unquoted column name for phosphosite probabilities.
#' @param num_phospho_site_col Unquoted column name for the number of phosphosites.
#'
#' @return A tibble where rows from paralogous proteins have been collapsed.
#'
#' @importFrom dplyr group_by summarise ungroup select one_of across
#' @importFrom rlang as_name enquo
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



#' Pivot Phosphosite Table to Wide Format
#'
#' @description
#' Converts a long-format phosphosite table back to a wide format, where each
#' sample/replicate becomes a separate column.
#'
#' @param all_phos_sites_long_tbl The long-format input tibble.
#' @param additional_cols Additional columns to use in creating the new column names.
#' @param phospho_site_prob_col Unquoted column name for phosphosite probabilities.
#' @param num_phospho_site_col Unquoted column name for the number of phosphosites.
#'
#' @return A wide-format tibble.
#'
#' @importFrom dplyr mutate_at pivot_wider arrange all_of
#' @importFrom rlang as_name enquo
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



#' Summarize Unique Phosphosites (Long Format)
#'
#' @description
#' Groups the data by unique phosphosites and summarizes the abundance values
#' (e.g., by mean, median, sum). Returns a list of summarized tibbles in long format.
#'
#' @param all_phos_sites_long_tbl A long-format phosphosite table.
#' @param additional_cols Additional columns to include in the grouping.
#'
#' @return A named list of tibbles (`mean`, `median`, `sum`), where each tibble
#'   contains the summarized abundance values in long format.
#'
#' @importFrom dplyr group_by summarise ungroup mutate mutate_at across
#' @importFrom purrr map
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


#' Summarize Unique Phosphosites (Wide Format)
#'
#' @description
#' Takes the list of long-format summarized tables and pivots them to wide format.
#' It also collapses evidence IDs from different experiments into a single identifier.
#'
#' @param summarised_long_tbl_list A list of long-format tibbles from
#'   `uniquePhosphositesSummariseLongList`.
#' @param additional_cols Additional columns used for pivoting.
#'
#' @return A named list of wide-format tibbles (`mean`, `median`, `sum`) with
#'   unique phosphosites as rows and samples as columns.
#'
#' @importFrom dplyr group_by summarise ungroup select left_join relocate mutate mutate_at summarise_all
#' @importFrom tidyr pivot_wider unite
#' @importFrom purrr map
#' @importFrom stringr str_replace_all
#' @importFrom rlang sym
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



#' Process Multi-site Phosphoproteomics Evidence
#'
#' @description
#' A high-level wrapper function that orchestrates the entire phosphoproteomics
#' data processing workflow, from reading raw files to generating summarized,
#' wide-format tables of unique phosphosites.
#'
#' @param fasta_file Path to the FASTA file containing protein sequences.
#' @param evidence_tbl The `evidence.txt` data table from MaxQuant.
#' @param accession_col Unquoted column name for protein accessions.
#' @param group_id Unquoted column name for the grouping variable.
#' @param additional_cols Unquoted column name(s) for additional grouping variables.
#' @param col_pattern Regex pattern for abundance columns.
#' @param extract_pattern (Currently unused) Regex pattern for extraction.
#' @param col_to (Currently unused) Column name for extracted values.
#' @param site_prob_threshold Primary probability threshold for phosphosites.
#' @param columns_to_use (Currently unused) Specific columns to use.
#'
#' @return A list containing four main data structures:
#'   - `summarised_wide_list`: A list of wide-format tables (mean, median, sum).
#'   - `summarised_long_list`: A list of long-format tables (mean, median, sum).
#'   - `all_phos_sites_wide`: A wide table of all filtered phosphosites before summarization.
#'   - `all_phos_sites_long`: A long table of all filtered phosphosites before summarization.
#'
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

#' Get UniProt Accession Rank from Site ID
#'
#' @description Given an input table with `sites_id` and `uniprot_acc` columns,
#'   this function determines the rank (position) of a specific UniProt accession
#'   within the list of accessions contained in the `sites_id`.
#'
#' @details The `sites_id` is expected to be a composite identifier where the
#'   first element (before the "!" separator) is a colon-separated list of
#'   UniProt accessions. This function finds which position the primary
#'   `uniprot_acc` for the row occupies in that list.
#'
#' @param input_table An input data frame with `sites_id` and `uniprot_acc` columns.
#' @param uniprot_acc The unquoted column name for the primary UniProt accession of the row.
#' @param sites_id The unquoted column name for the composite site identifier.
#'
#' @return The input table with new columns: `uniprot_acc_split` (the primary accession),
#'   `uniprot_list` (the list of accessions from `sites_id`), and `gene_list_position`
#'   (the 1-based index of the primary accession in the list).
#'
#' @importFrom dplyr mutate relocate
#' @importFrom stringr str_split
#' @importFrom purrr map_chr map2_int
#' @importFrom lazyeval as_name
#' @importFrom rlang enquo
#' @export
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
