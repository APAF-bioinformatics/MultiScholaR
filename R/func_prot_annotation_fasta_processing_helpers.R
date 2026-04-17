# ----------------------------------------------------------------------------
# processFastaFile
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
# ----------------------------------------------------------------------------
# processFastaFile
# ----------------------------------------------------------------------------
#' @title Process FASTA file and map to UniProt/UniParc
#' @description Parses a FASTA file and maps accessions to UniProt/UniParc search results
#' to create a standardized sequence mapping table.
#' @param fasta_file_path Path to the FASTA file.
#' @param uniprot_search_results UniProt search results data frame (optional).
#' @param uniparc_search_results UniParc search results data frame (optional).
#' @param fasta_meta_file Path to save/load parsed fasta metadata.
#' @param organism_name Name of the organism for metadata.
#' @return Mapping table between accessions and standardized IDs.
#' @export
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
          isoform_num = isoform_num,
          annotation_score = NA_real_
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
      attributes = attributes,
      annotation_score = NA_real_
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
    
    # Create metadata for standard UniProt format
    fasta_metadata <- list(
      fasta_format = "standard_uniprot",
      available_columns = colnames(aa_seq_tbl),
      has_protein_evidence = "protein_evidence" %in% colnames(aa_seq_tbl),
      has_gene_names = "gene_name" %in% colnames(aa_seq_tbl),
      has_isoform_info = "is_isoform" %in% colnames(aa_seq_tbl),
      has_status_info = "status" %in% colnames(aa_seq_tbl),
      num_sequences = nrow(aa_seq_tbl),
      processing_timestamp = Sys.time()
    )
    
    message("Saving results...")
    flush.console()
    saveRDS(aa_seq_tbl, fasta_meta_file)
    
    # Return list with data and metadata
    return(list(
      aa_seq_tbl_final = aa_seq_tbl,
      fasta_metadata = fasta_metadata
    ))
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
    
    # Create metadata for non-standard format
    fasta_metadata <- list(
      fasta_format = "non_standard",
      available_columns = colnames(aa_seq_tbl_final),
      has_protein_evidence = "protein_evidence" %in% colnames(aa_seq_tbl_final),
      has_gene_names = "gene_name" %in% colnames(aa_seq_tbl_final),
      has_isoform_info = "is_isoform" %in% colnames(aa_seq_tbl_final),
      has_status_info = "status" %in% colnames(aa_seq_tbl_final),
      num_sequences = nrow(aa_seq_tbl_final),
      processing_timestamp = Sys.time()
    )

    message("Writing results...")
    flush.console()

    vroom::vroom_write(aa_seq_tbl_final,
                       file = "aa_seq_tbl.tsv",
                       delim = "\t",
                       na = "",
                       quote = "none",
                       progress = FALSE)

    saveRDS(aa_seq_tbl_final, fasta_meta_file)
    
    # Return list with data and metadata
    return(list(
      aa_seq_tbl_final = aa_seq_tbl_final,
      fasta_metadata = fasta_metadata
    ))
  }
}

