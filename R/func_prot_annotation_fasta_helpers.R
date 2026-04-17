# ----------------------------------------------------------------------------
# getFastaFields
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# parseFastaObject
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# parseFastaFile
# ----------------------------------------------------------------------------
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

