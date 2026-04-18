# ----------------------------------------------------------------------------
# chooseBestPhosphositeAccession
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
chooseBestPhosphositeAccession <- function(input_tbl, acc_detail_tab, accessions_column, group_id) {

  # Join with FASTA data
  resolve_acc_joined <- input_tbl %>%
    dplyr::select( {{group_id}}, {{accessions_column}}, cleaned_peptide) %>%
    mutate( uniprot_acc = str_split( {{accessions_column}}, ";") ) %>%
    unnest( uniprot_acc )   %>%
    mutate( cleaned_acc = cleanIsoformNumber(uniprot_acc)) %>%
    left_join( acc_detail_tab,
               by=c("uniprot_acc" = "uniprot_acc",
                    "cleaned_acc" = "cleaned_acc") )
  
  # Detect which columns are available and add missing ones with defaults
  available_cols <- colnames(resolve_acc_joined)
  
  if (!"gene_name" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(gene_name = NA_character_)
  }
  
  if (!"protein_evidence" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(protein_evidence = 3)
  }
  
  if (!"status" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(status = "unknown")
  }
  
  if (!"is_isoform" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(is_isoform = "Canonical")
  }
  
  if (!"isoform_num" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(isoform_num = 0)
  }
  
  if (!"seq_length" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(seq_length = NA_integer_)
  }
  
  if (!"seq" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(seq = NA_character_)
  }
  
  # Now select, filter, and arrange using the complete set of columns
  resolve_acc_helper <- resolve_acc_joined %>%
    ## Just a sanity check that the peptide is actually in the sequence (skip if seq not available)
    {if ("seq" %in% colnames(.) && all(!is.na(.$seq))) filter(., str_detect( seq, cleaned_peptide  )) else .} %>%
    dplyr::select({{group_id}}, one_of(c( "uniprot_acc", "gene_name", "cleaned_acc",
                                          "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length"  ))) %>%
    distinct %>%
    mutate(annotation_score = if_else(is.na(annotation_score), 0, annotation_score)) %>%
    arrange( {{group_id}}, desc(annotation_score), protein_evidence, status, is_isoform, desc(seq_length), isoform_num )

  # print( colnames(head(resolve_acc_helper)) )


  score_isoforms <- resolve_acc_helper %>%
    mutate( gene_name = ifelse( is.na(gene_name) | gene_name == "", "NA", gene_name)) %>%
    mutate(annotation_score = if_else(is.na(annotation_score), 0, annotation_score)) %>%
    group_by( {{group_id}},  gene_name ) %>%
    arrange( {{group_id}},  desc(annotation_score), protein_evidence,
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

