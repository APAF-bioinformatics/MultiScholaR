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
parseFastaFile <- function(fasta_file_path, uniprot_search_results = NULL, uniparc_search_results = NULL, fasta_meta_file, organism_name) {
  startsWith <- function(x, prefix) {
    substr(x, 1, nchar(prefix)) == prefix
  }
  
  parseFastaFile <- function(fasta_file) {
  aa_seqinr <- seqinr::read.fasta(file = fasta_file, seqtype = "AA", 
                          whole.header = TRUE, as.string = TRUE)
  headers <- names(aa_seqinr)
  
  parsed_headers <- lapply(headers, function(header) {
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
  })
  
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
    aa_seqinr <- seqinr::read.fasta(file = fasta_file, seqtype = "AA", 
                            whole.header = TRUE, as.string = TRUE)
    headers <- names(aa_seqinr)
    parsed_headers <- lapply(headers, parseFastaHeader)
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

  fasta_file_raw <- vroom::vroom(fasta_file_path, delim = "\n", col_names = FALSE)
  first_line <- fasta_file_raw$X1[1]

  if (startsWith(first_line, ">sp|") || startsWith(first_line, ">tr|")) {
    aa_seq_tbl <- parseFastaFile(fasta_file_path)
    saveRDS(aa_seq_tbl, fasta_meta_file)
    return(aa_seq_tbl)
  } else {
    aa_seq_tbl <- parseFastaFileNonStandard(fasta_file_path)
    
    if (!is.null(uniprot_search_results) && !is.null(uniparc_search_results)) {
      aa_seq_tbl_final <- matchAndUpdateDataFrames(aa_seq_tbl, uniprot_search_results, uniparc_search_results, organism_name)
    } else {
      aa_seq_tbl_final <- aa_seq_tbl |>
        dplyr::mutate(database_id = NA_character_)
    }

    vroom::vroom_write(aa_seq_tbl_final,
                      file = "aa_seq_tbl.tsv",
                      delim = "\t",
                      na = "",
                      quote = "none")

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