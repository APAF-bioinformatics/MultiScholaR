# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children's Medical Research Institute, finding cures for childhood genetic diseases

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
cleanIsoformNumber <- function(string ) {
  # "Q8K4R4-2"
  str_replace( string, "-\\d+$", "")

}

# clean_isoform_number("Q8K4R4-2")



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
getFastaFields <- function(string, pattern) {

  field_found <- str_detect( {{string}}, paste0(pattern, "="))

  extract_data <- NA_character_
  if( field_found ) {

    extract_data <- str_replace_all( {{string}},
                                     paste0("(.*)",
                                            pattern,
                                            "=(.*?)(\\s..=.*|$)"), "\\2")
  }

  case_when(  field_found ~ extract_data,
              TRUE ~ NA_character_  )
}

#  get_fasta_fields( "6PGL_RAT 6-phosphogluconolactonase OS=Rattus norvegicus GN=Pgls PE=1 SV=1", "GN")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Parse FASTA object from seqinr
#' @description parse_fasta_object: Parse FASTA headers
#' @param aa_seq AAStringSet object, output from running seqinr
#' @return A table containing the protein evidence, isoform number, uniprot accession without isoform number in the uniprot_acc column, gene name
#' @export
parseFastaObject <- function(aa_seq ) {

  accession_tab <-  data.frame( header=names(aa_seq)) %>%
    separate( header, into=c("db", "uniprot_acc", "description"), sep="\\|") %>%
    mutate( uniprot_id = str_replace( description, "(.*?)\\s(.*)", "\\1" ) ) %>%
    mutate( OS = purrr::map_chr(description, ~getFastaFields(., "OS")))  %>%
    mutate( OX = purrr::map_int(description, ~as.integer(getFastaFields(., "OX")))) %>%
    mutate( GN = purrr::map_chr(description, ~getFastaFields(., "GN"))) %>%
    mutate( GN = ifelse( is.na(GN), "", GN)) %>%
    mutate( PE = purrr::map_int(description, ~as.integer(getFastaFields(., "PE")))) %>%
    mutate( SV = purrr::map_int(description, ~as.integer(getFastaFields(., "SV")))) %>%
    dplyr::select(-description) %>%
    dplyr::rename( species = "OS",
                   tax_id = "OX",
                   gene_name = "GN",
                   protein_evidence = "PE",
                   sequence_version = "SV")

  acc_detail_tab <- accession_tab %>%
    mutate( is_isoform = case_when( str_detect( uniprot_acc, "-\\d+") ~ "Isoform",
                                    TRUE ~ "Canonical") ) %>%
    mutate (isoform_num = case_when ( is_isoform == "Isoform" ~ str_replace_all( uniprot_acc,
                                                                                 "(.*)(-)(\\d{1,})",
                                                                                 "\\3") %>%
                                        as.numeric,
                                      is_isoform == "Canonical" ~ 0,
                                      TRUE ~ NA_real_ ) ) %>%
    mutate( cleaned_acc = cleanIsoformNumber(uniprot_acc)) %>%
    mutate( protein_evidence  = factor(protein_evidence, levels =1:5 )) %>%
    mutate( status = factor( db, levels =c( "sp", "tr"), labels=c("reviewed", "unreviewed"))) %>%
    mutate( is_isoform = factor(is_isoform, levels =c("Canonical", "Isoform")))

  return(acc_detail_tab )

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



#'Parse the headers of a Uniprot FASTA file and extract the headers and sequences into a data frame
#'Use seqinr object instead as it seems to be a lot faster to run substring
#' @param path to input faster file with header format described in https://www.uniprot.org/help/fasta-headers
#' @return A table containing the following columns:
#' db  sp for Swiss-Prot, tr for TrEMBL
#' uniprot_acc Uniprot Accession
#' uniprot_id  Uniprot ID
#' species     Species
#' tax_id      Taxonomy ID
#' gene_name   Gene symbol
#' protein_evidence 1 to 5, the lower the value, the more evidence that supports the existence of this protein
#' sequence_version Sequence version
#' is_isoform  Is it a protein isoform (not the canonical form)
#' isoform_num     Isoform number.
#' cleaned_acc Cleaned accession without isoform number.
#' status  Reviewed or unreviewed.
#' seq     Amino acid sequence.
#' seq_length      Sequence length (integer).
#' @export
parseFastaFile <- function(fasta_file) {

  aa_seqinr <-  read.fasta( file = fasta_file,
                            seqtype="AA",
                            whole.header	=TRUE,
                            as.string=TRUE)

  acc_detail_tab <- parseFastaObject(aa_seqinr)

  names(aa_seqinr) <- str_match( names(aa_seqinr), "(sp|tr)\\|(.+?)\\|(.*)\\s+" )[,3]

  aa_seq_tbl <- acc_detail_tab %>%
    mutate(seq = map_chr( aa_seqinr, 1)) %>%
    mutate(seq_length = purrr::map_int(seq, str_length) )
}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
chooseBestPhosphositeAccession <- function(input_tbl, acc_detail_tab, accessions_column, group_id) {

  resolve_acc_helper <- input_tbl %>%
    dplyr::select( {{group_id}}, {{accessions_column}}, cleaned_peptide) %>%
    mutate( uniprot_acc = str_split( {{accessions_column}}, ";") ) %>%
    unnest( uniprot_acc )   %>%
    mutate( cleaned_acc = cleanIsoformNumber(uniprot_acc)) %>%
    left_join( acc_detail_tab,
               by=c("uniprot_acc" = "uniprot_acc",
                    "cleaned_acc" = "cleaned_acc") ) %>%
    ## Just a sanity check that the peptide is actually in the sequence
    dplyr::filter( str_detect( seq, cleaned_peptide  )) %>%
    dplyr::select({{group_id}}, one_of(c( "uniprot_acc", "gene_name", "cleaned_acc",
                                          "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length"  ))) %>%

    distinct %>%
    arrange( {{group_id}}, desc(protein_evidence), status, is_isoform, desc(seq_length), isoform_num )

  # print( colnames(head(resolve_acc_helper)) )


  score_isoforms <- resolve_acc_helper %>%
    mutate( gene_name = ifelse( is.na(gene_name) | gene_name == "", "NA", gene_name)) %>%
    group_by( {{group_id}},  gene_name ) %>%
    arrange( {{group_id}},  desc(protein_evidence),
             status, is_isoform, desc(seq_length), isoform_num, cleaned_acc )  %>%
    mutate(ranking = row_number()) %>%
    ungroup


  # print( colnames(head(score_isoforms)) )

  ## For each gene name find the uniprot_acc with the lowest ranking
  group_gene_names_and_uniprot_accs <- score_isoforms %>%
    distinct( {{group_id}}, gene_name, ranking ) %>%
    dplyr::filter( ranking == 1) %>%
    left_join( score_isoforms %>%
                 dplyr::select( {{group_id}}, ranking, gene_name, uniprot_acc),
               by = join_by( {{group_id}} == {{group_id}}
                             , ranking == ranking
                             , gene_name == gene_name ) )   %>%
    dplyr::select(-ranking)

  # %>%
  #   group_by({{group_id}}) %>%
  #   summarise( num_gene_names = n(),
  #              gene_names = paste( gene_name, collapse=":"),
  #              uniprot_acc = paste( uniprot_acc, collapse=":")) %>%
  #   ungroup() %>%
  #   mutate( is_unique = case_when( num_gene_names == 1 ~ "Unique",
  #                                  TRUE ~ "Multimapped"))


  return( group_gene_names_and_uniprot_accs )

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@description From a list of UniProt accessions, choose the best accession to use based on the UniProt score for quality of annotation for the protein entries
#'@param input_tbl Contain the following columns, 'group_id' which is the Id for each protein group, 'accessions_column' which is the column with the accession of the protein
#'@param acc_detail_tabl The out table from running the function 'parseFastaFile'
#'@param accessions_column The name of the column with the list of protein accessions, separated by ';' semi-colon. No need to quote the name as we are using tidyverse programming quosure.
#'@param group_id The name of the column with the group ID for each protein group. No need to quote the name as we are using tidyverse programming quosure.
#' @returns A table with the following columns:
#'  maxquant_row_id: Row ID
#'  num_gene_names: Number of gene names associated with this row ID
#'  gene_names: The gene names
#'  uniprot_acc: List of uniprot accessions, but with the list ordered by the best one to less useful one to use
#'  is_unique: Is the protein group assined to a unique UniProt accession or multiple UniProt accessions
#'@export
chooseBestProteinAccessionHelper <- function(input_tbl
                                             , acc_detail_tab
                                             , accessions_column
                                             , row_id_column = "uniprot_acc"
                                             , group_id
                                             , delim= ";") {


  resolve_acc_temp <- input_tbl |>
    dplyr::select( { { group_id } }, { { accessions_column } }) |>
    mutate(row_id_column_with_isoform = str_split({ { accessions_column } }, delim)) |>
    unnest( row_id_column_with_isoform ) |>
    mutate( !!sym(row_id_column) := cleanIsoformNumber( row_id_column_with_isoform))

  print(head( resolve_acc_temp))

  print(head( acc_detail_tab))

  resolve_acc_helper <- resolve_acc_temp |>
    left_join( acc_detail_tab ,
               by = join_by( !!sym(row_id_column) == !!sym(row_id_column) ),
               copy = TRUE,
               keep = NULL)  |>
    dplyr::select( { { group_id } }, one_of(c(row_id_column, "gene_name", "cleaned_acc",
                                              "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length"))) |>
    distinct() |>
    arrange( { { group_id } }, desc(protein_evidence), status, is_isoform, desc(seq_length), isoform_num)


  score_isoforms <- resolve_acc_helper |>
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", "NA", gene_name)) |>
    group_by({ { group_id } }, gene_name) |>
    arrange( { { group_id } }, desc(protein_evidence),
             status, is_isoform, desc(seq_length), isoform_num, cleaned_acc) |>
    mutate( ranking = row_number()) |>
    ungroup()


  ## For each gene name find the uniprot_acc with the lowest rankinG
  group_gene_names_and_uniprot_accs <- score_isoforms |>
    distinct( { { group_id } }, gene_name, ranking) |>
    dplyr::filter(ranking == 1) |>
    left_join(score_isoforms |>
                dplyr::select({ { group_id } }, ranking, gene_name, !!sym(row_id_column)),
              by = join_by( {{ group_id }} == {{ group_id }}
                            , ranking == ranking
                            , gene_name == gene_name)) |>

    dplyr::select(-ranking) |>
    group_by({ { group_id } }) |>
    summarise(num_gene_names = n(),
              gene_names = paste(gene_name, collapse = ":"),
              !!sym(row_id_column) := paste(!!sym(row_id_column), collapse = ":")) |>
    ungroup() |>
    mutate(is_unique = case_when(num_gene_names == 1 ~ "Unique",
                                 TRUE ~ "Multimapped"))


  return(group_gene_names_and_uniprot_accs)

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@description From a list of UniProt accessions, rank the accession to use based on the UniProt score for quality of annotation for the protein entries
#'@param input_tbl Contain the following columns, 'group_id' which is the Id for each protein group, 'accessions_column' which is the column with the accession of the protein
#'@param acc_detail_tabl The out table from running the function 'parseFastaFile'
#'@param accessions_column The name of the column with the list of protein accessions, separated by ';' semi-colon. No need to quote the name as we are using tidyverse programming quosure.
#'@param group_id The name of the column with the group ID for each protein group. No need to quote the name as we are using tidyverse programming quosure.
#' @returns A table with the following columns:
#'  maxquant_row_id: Row ID
#'  num_gene_names: Number of gene names associated with this row ID
#'  gene_names: The gene names
#'  uniprot_acc: List of uniprot accessions, but with the list ordered by the best one to less useful one to use
#'  is_unique: Is the protein group assined to a unique UniProt accession or multiple UniProt accessions
#'@export
rankProteinAccessionHelper <- function(input_tbl
                                       , acc_detail_tab
                                       , accessions_column
                                       , row_id_column = "uniprot_acc"
                                       , group_id
                                       , delim= ";") {

  resolve_acc_helper <- input_tbl |>
    dplyr::select( { { group_id } }, { { accessions_column } }) |>
    mutate( !!sym(row_id_column) := str_split({ { accessions_column } }, delim)) |>
    unnest( !!sym(row_id_column)) |>
    mutate( cleaned_acc = cleanIsoformNumber(row_id_column))   |>
    left_join( acc_detail_tab ,
               by = join_by( cleaned_acc == !!sym(row_id_column) ),
               copy = TRUE,
               keep = NULL)  |>
    dplyr::select( { { group_id } }, one_of(c(row_id_column, "gene_name", "cleaned_acc",
                                              "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length"))) |>
    distinct() |>
    arrange( { { group_id } }, desc(protein_evidence), status, is_isoform, desc(seq_length), isoform_num)


  score_isoforms <- resolve_acc_helper |>
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", "NA", gene_name)) |>
    group_by({ { group_id } }, gene_name) |>
    arrange( { { group_id } }, desc(protein_evidence),
             status, is_isoform, desc(seq_length), isoform_num, cleaned_acc) |>
    mutate( ranking = row_number()) |>
    ungroup()


  ## For each gene name find the uniprot_acc with the lowest rankinG
  group_gene_names_and_uniprot_accs <- score_isoforms |>
    distinct( { { group_id } }, gene_name, ranking) |>
    left_join(score_isoforms |>
                dplyr::select({ { group_id } }, ranking, gene_name, !!sym(row_id_column)),
              by = join_by( {{ group_id }} == {{ group_id }}
                            , ranking == ranking
                            , gene_name == gene_name)) |>

    dplyr::select(-ranking) |>
    group_by({ { group_id } }) |>
    summarise(num_gene_names = n(),
              gene_names = paste(gene_name, collapse = ":"),
              !!sym(row_id_column) := paste(!!sym(row_id_column), collapse = ":")) |>
    ungroup() |>
    mutate(is_unique = case_when(num_gene_names == 1 ~ "Unique",
                                 TRUE ~ "Multimapped"))


  return(group_gene_names_and_uniprot_accs)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
processFastaFile <- function(fasta_file_path, uniprot_search_results = NULL, uniparc_search_results = NULL, fasta_meta_file, organism_name) {
  # Properly suppress all vroom messages
  withr::local_options(list(
    vroom.show_col_types = FALSE,
    vroom.show_progress = FALSE
  ))

  startsWith <- function(x, prefix) {
    substr(x, 1, nchar(prefix)) == prefix
  }

  parseFastaFileStandard <- function(fasta_file) {
    message("Reading FASTA file with seqinr...")
    utils::flush.console()

    aa_seqinr <- seqinr::read.fasta(file = fasta_file, seqtype = "AA",
                                    whole.header = TRUE, as.string = TRUE)
    headers <- names(aa_seqinr)
    total_entries <- length(headers)

    message(sprintf("\nProcessing %d FASTA entries...", total_entries))
    utils::flush.console()

    # Create a text progress bar
    pb <- utils::txtProgressBar(min = 0, max = total_entries, style = 3, width = 50)

    parsed_headers <- vector("list", length(headers))

    for(i in seq_along(headers)) {
      header <- headers[i]
      parsed_headers[[i]] <- {
        parts <- strsplit(substr(header, 2, nchar(header)), " ", fixed = TRUE)[[1]]
        id_parts <- strsplit(parts[1], "|", fixed = TRUE)[[1]]

        # Extract protein evidence level
        protein_evidence <- stringr::str_extract(header, "PE=[0-9]") |>
          stringr::str_extract("[0-9]") |>
          as.integer()

        # Determine status based on entry type
        status <- if(startsWith(header, ">sp|")) "reviewed" else "unreviewed"

        # Extract gene name (GN=)
        gene_name <- stringr::str_extract(header, "GN=\\S+") |>
          stringr::str_remove("GN=")

        # For entries without isoforms, set defaults
        is_isoform <- FALSE
        isoform_num <- 0L
        cleaned_acc <- id_parts[2]

        list(
          accession = id_parts[2],
          database_id = id_parts[2],
          cleaned_acc = cleaned_acc,
          gene_name = gene_name,
          protein = paste(parts[-1], collapse = " "),
          attributes = paste(parts[-1], collapse = " "),
          protein_evidence = protein_evidence,
          status = status,
          is_isoform = is_isoform,
          isoform_num = isoform_num
        )
      }

      # Update progress bar every 100 entries
      if(i %% 100 == 0 || i == total_entries) {
        utils::setTxtProgressBar(pb, i)
      }
    }

    close(pb)

    message("\nBinding rows and creating final table...")
    utils::flush.console()

    acc_detail_tab <- dplyr::bind_rows(parsed_headers)
    aa_seq_tbl <- acc_detail_tab |>
      dplyr::mutate(
        seq = purrr::map_chr(aa_seqinr, 1),
        seq_length = stringr::str_length(seq),
        description = headers
      )

    return(aa_seq_tbl)
  }

  parseFastaHeader <- function(header) {
    parts <- strsplit(substr(header, 2, nchar(header)), " ", fixed = TRUE)[[1]]
    id_parts <- strsplit(parts[1], "|", fixed = TRUE)[[1]]
    accession <- id_parts[2]
    attributes <- paste(parts[-1], collapse = " ")
    locus_tag <- stringr::str_extract(attributes, "(?<=\\[locus_tag=)[^\\]]+")
    protein <- stringr::str_extract(attributes, "(?<=\\[protein=)[^\\]]+")
    ncbi_refseq <- stringr::str_extract(attributes, "(?<=\\[protein_id=)WP_[^\\]]+")
    list(
      accession = accession,
      protein_id = locus_tag,
      protein = protein,
      ncbi_refseq = ncbi_refseq,
      attributes = attributes
    )
  }

  parseFastaFileNonStandard <- function(fasta_file) {
    message("Reading FASTA file with seqinr...")
    utils::flush.console()

    aa_seqinr <- seqinr::read.fasta(file = fasta_file, seqtype = "AA",
                                    whole.header = TRUE, as.string = TRUE)
    headers <- names(aa_seqinr)
    total_entries <- length(headers)

    message(sprintf("\nProcessing %d non-standard FASTA entries...", total_entries))
    utils::flush.console()

    # Create a text progress bar
    pb <- utils::txtProgressBar(min = 0, max = total_entries, style = 3, width = 50)

    parsed_headers <- vector("list", length(headers))

    for(i in seq_along(headers)) {
      header <- headers[i]
      parsed_headers[[i]] <- parseFastaHeader(header)

      # Update progress bar every 100 entries
      if(i %% 100 == 0 || i == total_entries) {
        utils::setTxtProgressBar(pb, i)
      }
    }

    close(pb)

    message("\nBinding rows and creating final table...")
    utils::flush.console()

    acc_detail_tab <- dplyr::bind_rows(parsed_headers)
    aa_seq_tbl <- acc_detail_tab |>
      dplyr::mutate(
        seq = purrr::map_chr(aa_seqinr, 1),
        seq_length = stringr::str_length(seq),
        description = headers
      )

    return(aa_seq_tbl)
  }

  matchAndUpdateDataFrames <- function(aa_seq_tbl, uniprot_search_results, uniparc_search_results, organism_name) {
    message("Matching and updating dataframes...")
    flush.console()

    uniprot_filtered <- uniprot_search_results |>
      dplyr::filter(Organism == organism_name) |>
      dplyr::select("ncbi_refseq", "uniprot_id")

    uniparc_prepared <- uniparc_search_results |>
      dplyr::select(ncbi_refseq, uniparc_id = uniprot_id)

    aa_seq_tbl_updated <- aa_seq_tbl |>
      dplyr::left_join(uniprot_filtered, by = "ncbi_refseq") |>
      dplyr::left_join(uniparc_prepared, by = "ncbi_refseq") |>
      dplyr::mutate(database_id = dplyr::coalesce(uniprot_id, uniparc_id)) |>
      dplyr::select(-uniprot_id, -uniparc_id)

    return(aa_seq_tbl_updated)
  }

  message("Reading FASTA file...")
  flush.console()

  suppressMessages({
    fasta_file_raw <- vroom::vroom(fasta_file_path, delim = "\n", col_names = FALSE, progress = FALSE)
  })
  first_line <- fasta_file_raw$X1[1]

  if (startsWith(first_line, ">sp|") || startsWith(first_line, ">tr|")) {
    message("Processing standard UniProt FASTA format...")
    flush.console()
    aa_seq_tbl <- parseFastaFileStandard(fasta_file_path)
    message("Saving results...")
    flush.console()
    saveRDS(aa_seq_tbl, fasta_meta_file)
    return(aa_seq_tbl)
  } else {
    message("Processing non-standard FASTA format...")
    flush.console()
    aa_seq_tbl <- parseFastaFileNonStandard(fasta_file_path)

    if (!is.null(uniprot_search_results) && !is.null(uniparc_search_results)) {
      aa_seq_tbl_final <- matchAndUpdateDataFrames(aa_seq_tbl, uniprot_search_results, uniparc_search_results, organism_name)
    } else {
      aa_seq_tbl_final <- aa_seq_tbl |>
        dplyr::mutate(database_id = NA_character_)
    }

    message("Writing results...")
    flush.console()

    vroom::vroom_write(aa_seq_tbl_final,
                       file = "aa_seq_tbl.tsv",
                       delim = "\t",
                       na = "",
                       quote = "none",
                       progress = FALSE)

    saveRDS(aa_seq_tbl_final, fasta_meta_file)
    return(aa_seq_tbl_final)
  }
}

################################################################################################################################################################################################

processFastaFile_deprecated <- function(fasta_file_path, uniprot_search_results, uniparc_search_results, fasta_meta_file) {
  startsWith <- function(x, prefix) {
    substr(x, 1, nchar(prefix)) == prefix
  }
  fasta_file_raw <- vroom::vroom(fasta_file_path, delim = "\n", col_names = FALSE)
  first_line <- fasta_file_raw$X1[1]

  if (startsWith(first_line, ">sp|") || startsWith(first_line, ">tr|")) {
    aa_seq_tbl <- parseFastaFile(fasta_file_path)
    saveRDS(aa_seq_tbl, fasta_meta_file)
    return(aa_seq_tbl)
  } else {  # Custom parsing for non-standard headers
    parseFastaHeader <- function(header) {
      parts <- strsplit(substr(header, 2, nchar(header)), " ", fixed = TRUE)[[1]]
      id_parts <- strsplit(parts[1], "|", fixed = TRUE)[[1]]
      accession <- id_parts[2]
      attributes <- paste(parts[-1], collapse = " ")
      locus_tag <- str_extract(attributes, "(?<=\\[locus_tag=)[^\\]]+")
      protein <- str_extract(attributes, "(?<=\\[protein=)[^\\]]+")
      ncbi_refseq <- str_extract(attributes, "(?<=\\[protein_id=)WP_[^\\]]+")
      list(
        accession = accession,
        protein_id = locus_tag,
        protein = protein,
        ncbi_refseq = ncbi_refseq,
        attributes = attributes
      )
    }

    parseFastaFile <- function(fasta_file) {
      aa_seqinr <- read.fasta(file = fasta_file, seqtype = "AA",
                              whole.header = TRUE, as.string = TRUE)
      headers <- names(aa_seqinr)
      parsed_headers <- lapply(headers, parseFastaHeader)
      acc_detail_tab <- bind_rows(parsed_headers)
      aa_seq_tbl <- acc_detail_tab |>
        mutate(seq = map_chr(aa_seqinr, 1),
               seq_length = map_int(seq, str_length),
               description = headers)

      return(aa_seq_tbl)
    }

    aa_seq_tbl <- parseFastaFile(fasta_file_path)

    matchAndUpdateDataFrames <- function(aa_seq_tbl, uniprot_search_results, uniparc_search_results) {
      uniprot_filtered <- uniprot_search_results |>
        dplyr::filter(Organism == "Klebsiella variicola") |>
        dplyr::select("ncbi_refseq", "uniprot_id")

      uniparc_prepared <- uniparc_search_results |>
        dplyr::select(ncbi_refseq, uniparc_id = uniprot_id)

      aa_seq_tbl_updated <- aa_seq_tbl |>
        dplyr::left_join(uniprot_filtered, by = "ncbi_refseq") |>
        dplyr::left_join(uniparc_prepared, by = "ncbi_refseq") |>
        dplyr::mutate(database_id = dplyr::coalesce(uniprot_id, uniparc_id)) |>
        dplyr::select(-uniprot_id, -uniparc_id)

      return(aa_seq_tbl_updated)
    }

    aa_seq_tbl_final <- matchAndUpdateDataFrames(aa_seq_tbl, uniprot_search_results, uniparc_search_results)

    vroom::vroom_write(aa_seq_tbl_final,
                       file = "aa_seq_tbl.tsv",
                       delim = "\t",
                       na = "",
                       quote = "none")

    saveRDS(aa_seq_tbl_final, fasta_meta_file)
    return(aa_seq_tbl_final)
  }
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export

updateProteinIDs <- function(protein_data, aa_seq_tbl_final) {
  # Check if ncbi_refseq column exists in aa_seq_tbl_final
  if (!"ncbi_refseq" %in% colnames(aa_seq_tbl_final)) {
    message("No ncbi_refseq column found in aa_seq_tbl_final. Returning original data unchanged.")
    return(protein_data)
  }

  # Generic NCBI protein ID patterns - escaped special characters
  ncbi_patterns <- c(
    "WP_\\d+\\.?\\d*",                     # WP_123456789.1
    "[A-Z]{2}_\\d+\\.?\\d*",              # NP_123456.1, XP_123456.1
    "[A-Z]{3}\\d+\\.?\\d*",               # ABC12345.1
    "\\w+\\.\\d+_prot_\\w+_\\d+"          # Assembly specific patterns like NZ_LR130543.1_prot_ABC_123
  )

  pattern <- paste0("(", paste(ncbi_patterns, collapse = "|"), ")")

  protein_data <- protein_data |>
    dplyr::mutate(matching_id = stringr::str_extract(Protein.Ids, pattern))

  lookup_table <- aa_seq_tbl_final |>
    dplyr::mutate(matching_id = stringr::str_extract(accession, pattern)) |>
    dplyr::select(matching_id, database_id, ncbi_refseq)

  updated_protein_data <- protein_data |>
    dplyr::left_join(lookup_table, by = "matching_id") |>
    dplyr::mutate(Protein.Ids_new = coalesce(database_id, ncbi_refseq, Protein.Ids)) |>
    dplyr::select(-matching_id, -database_id, -ncbi_refseq)

  changes <- sum(updated_protein_data$Protein.Ids_new != updated_protein_data$Protein.Ids)
  cat("Number of Protein.Ids that would be updated:", changes, "\n")

  if (changes > 0) {
    cat("\nSample of changes:\n")
    changed <- which(updated_protein_data$Protein.Ids_new != updated_protein_data$Protein.Ids)
    sample_changes <- head(changed, 5)
    for (i in sample_changes) {
      cat("Old:", updated_protein_data$Protein.Ids[i], "-> New:", updated_protein_data$Protein.Ids_new[i], "\n")
    }
  }

  # Replace old Protein.Ids with new ones
  updated_protein_data$Protein.Ids <- updated_protein_data$Protein.Ids_new
  updated_protein_data$Protein.Ids_new <- NULL

  return(updated_protein_data)
}



#' Clean MaxQuant Protein Data
#'
#' This function processes and cleans protein data from MaxQuant output,
#' filtering based on peptide counts and removing contaminants.
#'
#' @param fasta_file Path to input FASTA file
#' @param raw_counts_file Path to MaxQuant proteinGroups.txt file
#' @param output_counts_file Name of cleaned counts table output file
#' @param accession_record_file Name of cleaned accession to protein group mapping file
#' @param column_pattern Pattern to match intensity columns (e.g., "Reporter intensity corrected")
#' @param group_pattern Pattern to identify experimental groups (default: "")
#' @param razor_unique_peptides_group_thresh Threshold for razor + unique peptides (default: 0)
#' @param unique_peptides_group_thresh Threshold for unique peptides (default: 1)
#' @param fasta_meta_file Name of FASTA metadata RDS file (default: "aa_seq_tbl.RDS")
#' @param output_dir Directory for results (default: "results/proteomics/clean_proteins")
#' @param tmp_dir Directory for temporary files (default: "cache")
#' @param log_file Name of log file (default: "output.log")
#' @param debug Enable debug output (default: FALSE)
#' @param silent Only print critical information (default: FALSE)
#' @param no_backup Deactivate backup of previous run (default: FALSE)
#' @return List containing cleaned data and statistics
#' @import tidyverse vroom magrittr knitr rlang optparse seqinr ProteomeRiver janitor tictoc configr logging
#' @export
cleanMaxQuantProteins <- function(
    fasta_file,
    raw_counts_file,
    output_counts_file = "counts_table_cleaned.tab",
    accession_record_file = "cleaned_accession_to_protein_group.tab",
    column_pattern = "Reporter intensity corrected",
    group_pattern = "",
    razor_unique_peptides_group_thresh = 0,
    unique_peptides_group_thresh = 1,
    fasta_meta_file = "aa_seq_tbl.RDS",
    output_dir = "results/proteomics/clean_proteins",
    tmp_dir = "cache",
    log_file = "output.log",
    debug = FALSE,
    silent = FALSE,
    no_backup = FALSE
) {
  tic()

  # Initialize argument list with direct parameters
  args <- list(
    fasta_file = fasta_file,
    raw_counts_file = raw_counts_file,
    output_counts_file = output_counts_file,
    accession_record_file = accession_record_file,
    column_pattern = column_pattern,
    group_pattern = group_pattern,
    razor_unique_peptides_group_thresh = razor_unique_peptides_group_thresh,
    unique_peptides_group_thresh = unique_peptides_group_thresh,
    fasta_meta_file = fasta_meta_file,
    output_dir = output_dir,
    tmp_dir = tmp_dir,
    log_file = log_file,
    debug = debug,
    silent = silent,
    no_backup = no_backup
  )

  # Create directories
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }
  if (!dir.exists(args$tmp_dir)) {
    dir.create(args$tmp_dir, recursive = TRUE)
  }

  # Configure logging
  logReset()
  addHandler(writeToConsole)
  addHandler(writeToFile, file = file.path(args$output_dir, args$log_file))

  level <- ifelse(args$debug, loglevels["DEBUG"], loglevels["INFO"])
  setLevel(level = ifelse(args$silent, loglevels["ERROR"], level))

  # Log start of processing
  loginfo("Starting protein data cleaning")

  # Validate required files
  required_files <- c(args$fasta_file, args$raw_counts_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
  }

  # Set default values for pattern suffixes
  args$pattern_suffix <- "_\\d+"
  args$extract_patt_suffix <- "_(\\d+)"
  args$remove_more_peptides <- FALSE

  # Read counts file
  loginfo("Reading the counts file")
  dat_tbl <- vroom::vroom(args$raw_counts_file)

  # Clean counts table header
  loginfo("Cleaning counts table header")
  dat_cln <- janitor::clean_names(dat_tbl)
  colnames(dat_cln) <- str_replace(colnames(dat_cln), "_i_ds", "_ids")

  # Prepare regular expressions
  pattern_suffix <- args$pattern_suffix
  if (args$group_pattern != "") {
    pattern_suffix <- paste(args$pattern_suffix, tolower(args$group_pattern), sep = "_")
  }

  extract_patt_suffix <- args$extract_patt_suffix
  if (args$group_pattern != "") {
    extract_patt_suffix <- paste0(args$extract_patt_suffix, "_(", tolower(args$group_pattern), ")")
  }

  column_pattern <- tolower(paste0(make_clean_names(args$column_pattern), pattern_suffix))
  extract_replicate_group <- tolower(paste0(make_clean_names(args$column_pattern), extract_patt_suffix))

  # Prepare peptide count columns
  razor_unique_peptides_group_col <- "razor_unique_peptides"
  unique_peptides_group_col <- "unique_peptides"

  if (args$group_pattern != "") {
    razor_unique_peptides_group_col <- paste0("razor_unique_peptides_", tolower(args$group_pattern))
    unique_peptides_group_col <- paste0("unique_peptides_", tolower(args$group_pattern))
  }

  # Process FASTA file
  fasta_meta_file <- file.path(args$tmp_dir, args$fasta_meta_file)
  loginfo("Processing FASTA file")

  if (file.exists(fasta_meta_file)) {
    aa_seq_tbl <- readRDS(fasta_meta_file)
  } else {
    aa_seq_tbl <- parseFastaFile(args$fasta_file)

    saveRDS(aa_seq_tbl, fasta_meta_file)
  }

  # Process and filter data
  evidence_tbl <- dat_cln %>%
    mutate(maxquant_row_id = id)

  print(aa_seq_tbl)

  # Filter and clean data
  loginfo("Identify best UniProt accession per entry, extract sample number and simplify column header")

  filtered_data <- processAndFilterData(
    evidence_tbl,
    args,
    razor_unique_peptides_group_col,
    unique_peptides_group_col,
    column_pattern,
    aa_seq_tbl,
    extract_replicate_group
  )

  # Save results
  saveResults(filtered_data, args)

  # Log completion and session info
  te <- toc(quiet = TRUE)
  loginfo("%f sec elapsed", te$toc - te$tic)
  writeLines(capture.output(sessionInfo()), file.path(args$output_dir, "sessionInfo.txt"))

  return(filtered_data)
}

#' Helper function to process and filter data
#' @noRd
processAndFilterData <- function(
    evidence_tbl,
    args,
    razor_unique_peptides_group_col,
    unique_peptides_group_col,
    column_pattern,
    aa_seq_tbl,
    extract_replicate_group
) {
  # Initialize tracking of protein numbers
  num_proteins_remaining <- numeric(3)
  names(num_proteins_remaining) <- c(
    "Number of proteins in raw unfiltered file",
    "Number of proteins after removing reverse decoy and contaminant proteins",
    paste0(
      "Number of proteins after removing proteins with no. of razor + unique peptides < ",
      args$razor_unique_peptides_group_thresh,
      " and no. of unique peptides < ",
      args$unique_peptides_group_thresh
    )
  )

  # Filter and process data
  select_columns <- evidence_tbl %>%
    dplyr::select(
      maxquant_row_id,
      protein_ids,
      !!rlang::sym(razor_unique_peptides_group_col),
      !!rlang::sym(unique_peptides_group_col),
      reverse,
      potential_contaminant,
      matches(column_pattern)
    )

  num_proteins_remaining[1] <- nrow(select_columns)

  remove_reverse_and_contaminant <- select_columns  %>%
    dplyr::filter( is.na(reverse) &
                     is.na(potential_contaminant)) %>%
    dplyr::filter( !str_detect(protein_ids, "^CON__") &
                     !str_detect(protein_ids, "^REV__") )

  remove_reverse_and_contaminant_more_hits <- remove_reverse_and_contaminant

  # Remove reverse decoy peptides and contaminant peptides even if it is not the first ranked Protein IDs (e.g. it is lower down in the list of protein IDs)
  if( args$remove_more_peptides == TRUE) {
    remove_reverse_and_contaminant_more_hits <- remove_reverse_and_contaminant  %>%
      dplyr::filter( is.na(reverse) &
                       is.na(potential_contaminant)) %>%
      dplyr::filter( !str_detect(protein_ids, "CON__") &
                       !str_detect(protein_ids, "REV__") )
  }

  # Record the number of proteins after removing reverse decoy and contaminant proteins
  # The numbers will be saved into the file 'number_of_proteins_remaining_after_each_filtering_step.tab'
  num_proteins_remaining[2] <- nrow(remove_reverse_and_contaminant_more_hits)

  helper_unnest_unique_and_razor_peptides <- remove_reverse_and_contaminant_more_hits %>%
    dplyr::mutate(protein_ids = str_split(protein_ids, ";")) %>%
    dplyr::mutate(!!rlang::sym(razor_unique_peptides_group_col) := str_split(!!rlang::sym(razor_unique_peptides_group_col), ";")) %>%
    dplyr::mutate(!!rlang::sym(unique_peptides_group_col) := str_split(!!rlang::sym(unique_peptides_group_col), ";")) %>%
    unnest(cols = c(protein_ids,
                    !!rlang::sym(razor_unique_peptides_group_col),
                    !!rlang::sym(unique_peptides_group_col)))


  evidence_tbl_cleaned <- helper_unnest_unique_and_razor_peptides %>%
    dplyr::filter(!!rlang::sym(razor_unique_peptides_group_col) >= args$razor_unique_peptides_group_thresh &
                    !!rlang::sym(unique_peptides_group_col) >= args$unique_peptides_group_thresh)


  ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  print(head(evidence_tbl_cleaned))


  accession_gene_name_tbl <- chooseBestProteinAccessionHelper(input_tbl = evidence_tbl_cleaned,
                                                        acc_detail_tab = aa_seq_tbl,
                                                        accessions_column = protein_ids,
                                                        row_id_column = "uniprot_acc",
                                                        group_id = maxquant_row_id)


  accession_gene_name_tbl_record <- accession_gene_name_tbl %>%
    left_join(evidence_tbl %>% dplyr::select(maxquant_row_id, protein_ids), by = c("maxquant_row_id"))


  evidence_tbl_filt <- evidence_tbl_cleaned |>
    inner_join(accession_gene_name_tbl |>
                 dplyr::select(maxquant_row_id, uniprot_acc), by = "maxquant_row_id") |>
    dplyr::select(uniprot_acc, matches(column_pattern), -contains(c("razor", "unique"))) |>
    distinct()

  # Record the number of proteins after removing proteins with low no. of razor + unique peptides and low no. of unique peptides
  num_proteins_remaining[3] <- nrow( evidence_tbl_filt)

  # Record the number of proteins remaining after each filtering step into the file 'number_of_proteins_remaining_after_each_filtering_step.tab'
  num_proteins_remaining_tbl <- data.frame( step=names( num_proteins_remaining), num_proteins_remaining=num_proteins_remaining)

  #TODO: This part need improvement. There is potential for bugs.
  extraction_pattern <- "\\1"
  if (args$group_pattern != "") {
    extraction_pattern <- "\\1_\\2"
  }

  colnames(evidence_tbl_filt) <- str_replace_all(colnames(evidence_tbl_filt), tolower(extract_replicate_group), extraction_pattern) %>%
    toupper( ) %>%
    str_replace_all( "UNIPROT_ACC", "uniprot_acc")


  return(list(
    evidence_tbl_filt = evidence_tbl_filt,
    num_proteins_remaining = num_proteins_remaining,
    accession_gene_name_tbl_record = accession_gene_name_tbl_record
  ))
}

#' Helper function to save results
#' @noRd
saveResults <- function(filtered_data, args) {
  # Save cleaned counts
  vroom::vroom_write(
    filtered_data$evidence_tbl_filt,
    file.path(args$output_dir, args$output_counts_file)
  )

  # Save accession records
  vroom::vroom_write(
    filtered_data$accession_gene_name_tbl_record,
    file.path(args$output_dir, args$accession_record_file)
  )

  # Save protein numbers
  vroom::vroom_write(
    data.frame(
      step = names(filtered_data$num_proteins_remaining),
      num_proteins_remaining = filtered_data$num_proteins_remaining
    ),
    file.path(args$output_dir, "number_of_proteins_remaining_after_each_filtering_step.tab")
  )

  # Save sample names
  sample_names <- colnames(filtered_data$evidence_tbl_filt)[-1]
  vroom::vroom_write(
    data.frame(sample_names = t(t(sample_names))),
    file.path(args$output_dir, "sample_names.tab")
  )
}
