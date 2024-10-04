# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases

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
#' @return A table containing the protein evidence, isoform number, uniprot accession without isoform number, gene name
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
    arrange( {{group_id}}, protein_evidence, status, is_isoform, desc(seq_length), isoform_num )

  # print( colnames(head(resolve_acc_helper)) )


  score_isoforms <- resolve_acc_helper %>%
    mutate( gene_name = ifelse( is.na(gene_name) | gene_name == "", "NA", gene_name)) %>%
    group_by( {{group_id}},  gene_name ) %>%
    arrange( {{group_id}},  protein_evidence,
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
    arrange( { { group_id } }, protein_evidence, status, is_isoform, desc(seq_length), isoform_num)


  score_isoforms <- resolve_acc_helper |>
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", "NA", gene_name)) |>
    group_by({ { group_id } }, gene_name) |>
    arrange( { { group_id } }, protein_evidence,
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
