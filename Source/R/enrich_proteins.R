## Packages installation and loading

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# load pacman package manager
if(!require(pacman)){
    install.packages("pacman")
    library(pacman)
}


## Load packages
p_load(tidyverse)
p_load(plotly)
p_load(vroom)
p_load(ggplot2)
p_load(ggpubr)
p_load(writexl)
p_load(magrittr)
p_load(knitr)
p_load(rlang)


p_load(limma)
p_load(qvalue)

p_load(UniProt.ws)
p_load(msigdb)
p_load(ExperimentHub)
p_load(GSEABase)
p_load(ProteomeRiver)
p_load(configr)
p_load(logging)
p_load(optparse)
p_load(tictoc)

## Parameters


## Directories management
base_dir <- here::here()

# output_dir <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/enrich_proteins"
# group_pattern <- "\\d+"
# de_proteins_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/de_proteins/de_proteins_long.tsv"
# contrasts_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Source/N155S/contrast_strings.tab"
# design_matrix_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Source/N155S/design_matrix.tab"
# counts_table_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/de_proteins/normalized_counts_after_ruv.tsv"
# formula_string <- "~ 0 + group "
# my_sample_id <- "Sample_ID"
# my_group_id <- "group"
# my_row_id <- "uniprot_acc"
# taxonomy_id <- 10090


command_line_options <- commandArgs(trailingOnly = TRUE)

  parser <- OptionParser(add_help_option =TRUE)


  parser <- add_option(parser, c( "--tax-id"), type="integer", default="", dest = "taxonomy_id",
                       help="The NCBI taxonomy ID of the organism being investigated (e.g. M. musculus=10090, H. sapien=9606).",
                       metavar="integer")

  parser <- add_option(parser, c( "--group-pattern"), type="character", default="", dest = "group_pattern",
                       help="Regular expression pattern to identify columns with abundance values belonging to the experiment. [default %default]",
                       metavar="string")

  parser <- add_option(parser, c("-c", "--counts"), type="character", default="", dest = "counts_table_file",
                       help="Input file with the protein abundance values",
                       metavar="string")

  parser <- add_option(parser, c("-c", "--de-proteins"), type="character", default="", dest = "de_proteins_file",
                       help="Input file with the list of diffierentiall expressed protein log fold-change and q-values for every contrasts.",
                       metavar="string")

  parser <- add_option(parser, c("--contrasts"), type="character", default="", dest = "contrasts_file",
                       help="Input file with a table listing all comparisons to be made in string, one comparison per line (e.g. groupB.vs.group_A = groupB - groupA).",
                       metavar="string")

  parser <- add_option(parser, c("--formula"), type="character", default="", dest = "formula_string",
                       help="A string representing the formula for input into the model.frame function. (e.g. ~ 0 + group).",
                       metavar="string")

  parser <- add_option(parser, c("-d", "--design-matrix"), type="character", default="", dest = "design_matrix_file",
                       help="Input file with the design matrix",
                       metavar="string")

  parser <- add_option(parser, c("-o", "--output-dir"), type="character", default="", dest = "output_dir",
                       help="Directory path for all results files.",
                       metavar="string")

  parser <- add_option(parser, c("--sample-id"), type="character", default="Sample_ID", dest = "my_sample_id",
                       help="A string describing the sample ID. This must be a column that exists in the design matrix.",
                       metavar="string")

  parser <- add_option(parser, c("-g", "--group-id"), type="character", default="group", dest = "my_group_id",
                       help="A string describing the experimental group ID. This must be a column that exists in the design matrix.",
                       metavar="string")

  parser <- add_option(parser, c("-r", "--row-id"), type="character", default="uniprot_acc", dest = "my_row_id",
                       help="A string describing the row id.",
                       metavar="string")


  print(commandArgs(trailingOnly = TRUE))

  args <- parse_args(parser)


  #parse and merge the configuration file options.
  if (args$config != "") {
    args <- config.list.merge(eval.config(file = args$config, config="enrich_camera_proteins"), args)
  }

  args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="enrich_camera_proteins" )

  createOutputDir(args$output_dir, args$no_backup)


  ## Logger configuration
  logReset()
  addHandler(writeToConsole , formatter = cmriFormatter)
  addHandler(writeToFile, file = file.path(args$output_dir, args$log_file), formatter = cmriFormatter)

  level <- ifelse(args$debug, loglevels["DEBUG"], loglevels["INFO"])
  setLevel(level = ifelse(args$silent, loglevels["ERROR"], level))

  cmriWelcome("ProteomeRiver", c("Ignatius Pang", "Pablo Galaviz"))
  loginfo("Reading configuration file %s", args$config)
  loginfo("Argument: Value")
  loginfo("----------------------------------------------------")
  for (v in names(args))
  {
    loginfo("%s : %s", v, args[v])
  }
  loginfo("----------------------------------------------------")

  print(args)

  testRequiredArguments(args, c(
    args$group_pattern,
    args$output_dir,
    args$my_sample_id,
    args$my_group_id,
    args$my_row_id,
    args$taxonomy_id,
  ))

  testRequiredFiles(c(
    args$design_matrix_file,
    args$counts_table_file,
    args$de_proteins_file,
    args$contrasts_file
  ))

  group_pattern <- args$group_pattern
  design_matrix_file <- args$design_matrix_file
  counts_table_file <- args$counts_table_file
  de_proteins_file <- args$de_proteins_file
  contrasts_file <- args$contrasts_file
  output_dir <- args$output_dir
  my_sample_id <- args$my_sample_id
  my_group_id <-  args$my_group_id
  my_row_id <- args$my_row_id
  taxonomy_id <- args$taxonomy_id

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Create directories
createDirectoryIfNotExists( output_dir)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Read the Contrasts file
contrasts_tbl <- NA
if( contrasts_file != "") {
  loginfo("Read file with lists of experimental contrasts to test %s", args$contrasts_file)
  captured_output<-capture.output(
    contrasts_tbl <- vroom::vroom( contrasts_file, delim="\t"),
    type = "message"
  )
  logdebug(captured_output)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## List of read abundance tables after normalization with RUVIII
loginfo("Read abundance tables after normalization with RUVIII.")
captured_output <- capture.output(
norm_abundance <- vroom::vroom( counts_table_file )  %>%
  mutate( best_uniprot_acc = purrr::map_chr( uniprot_acc, ~str_split(. , ":"   )[[1]][1] )),
type = "message"
)
logdebug(captured_output)

norm_abundance_mat <- norm_abundance %>%
  dplyr::select(best_uniprot_acc, matches( group_pattern ) ) %>%
  column_to_rownames("best_uniprot_acc")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Functions to create the replicate matrix
loginfo("Create the replicate matrix.")

ruvIII_replicates_matrix <- getRuvIIIReplicateMatrix( design_mat_cln,
                                                         !!rlang::sym(my_sample_id),
                                                         !!rlang::sym(my_group_id))

ruvIII_replicates_matrix


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Set up list of contrasts

loginfo("Create the list of contrasts.")

ff <- as.formula(  formula_string)
mod_frame <- model.frame( ff, design_mat_cln)
design_m <- model.matrix( ff, mod_frame)

contr.matrix <- makeContrasts( contrasts = contrasts_tbl %>% pull(contrasts),
                                 levels = colnames(design_m))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Obatain the Mol. Sig DB data.")

eh <- ExperimentHub()
# query(eh , 'msigdb')

org_abbrev <- case_when(taxonomy_id == 10090  ~ "mm",
          taxonomy_id == 9606 ~ "hs",
          TRUE ~ NA_character_)

if( is.na( org_abbrev ) ) {
  stop("Organism not yet supported by this code. Only human and mouse are supported.")
}

msigdb.EZID <- getMsigdb(org_abbrev, 'EZID')
msigdb.EZID <- appendKEGG(msigdb.EZID)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## List of all proteins

loginfo("Get the list of all DE proteins.")
captured_output <- capture.output(
proteins_cln <- vroom::vroom(de_proteins_file) %>%
  mutate( best_uniprot_acc = purrr::map_chr( uniprot_acc, ~str_split(. , ":"   )[[1]][1] )) %>%
  dplyr::select( -uniprot_acc ) %>%
  dplyr::rename( uniprot_acc = "best_uniprot_acc")
,type = "message"
)
logdebug(captured_output)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Convert Uniprot accession to Entrez gene ID

loginfo("Convert Uniprot accession to Entrez gene ID.")

if( file.exists(file.path( output_dir, "uniprot_acc_to_entrez_id.RDS") )) {

    uniprot_acc_to_entrez_id <- readRDS( file.path( output_dir, "uniprot_acc_to_entrez_id.RDS")) %>%
      mutate( ENTREZ_GENE= purrr::map_chr( ENTREZ_GENE,   ~str_replace_all(., "(^\\s+|\\s+$)", "")))

} else {
    up <- UniProt.ws(taxId=taxonomy_id )

    #keytypes(up)

    uniprot_acc_to_entrez_id <- batchQueryEvidence( proteins_cln %>% distinct(uniprot_acc), uniprot_acc, up,
                                                      uniprot_columns = c("ENTREZ_GENE"))  %>%
      mutate( ENTREZ_GENE= purrr::map_chr( ENTREZ_GENE,   ~str_replace_all(., "(^\\s+|\\s+$)", "")))

    saveRDS( uniprot_acc_to_entrez_id, file.path( output_dir, "uniprot_acc_to_entrez_id.RDS") )
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Convert Uniprot accession to Entrez ID.")

de_prot_for_camera_output <- mergeWithEntrezId( proteins_cln,
                                                   uniprot_acc_to_entrez_id  )

de_prot_for_camera_tab <-  de_prot_for_camera_output$de_proteins %>%
  group_by( comparison) %>%
  nest() %>%
  ungroup()

de_prot_for_camera_output$de_proteins %>%
  group_by(comparison, ENTREZ_GENE) %>%
  summarise(counts=n()) %>%
  ungroup() %>%
  dplyr::filter( counts > 1)

de_prot_for_camera_mapped <- de_prot_for_camera_tab %>%
     mutate(  data_edited = purrr::map( data, ~{ as.matrix(column_to_rownames(.,  "ENTREZ_GENE" ))   }) )

de_prot_for_camera_list <- de_prot_for_camera_tab$data
names( de_prot_for_camera_list) <- de_prot_for_camera_tab$comparison


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


norm_abundance_helper <- norm_abundance_mat %>%
  rownames_to_column("uniprot_acc") %>%
  as.data.frame %>%
  mutate ( temp =1) %>%
  inner_join( data.frame( comparison=unique(names(de_prot_for_camera_list)), temp=1 ), by="temp"        )  %>%
  dplyr::select( -temp) %>%
  anti_join(de_prot_for_camera_output$without_entrez_id, by=c("uniprot_acc" = "uniprot_acc",
                                                              "comparison" = "comparison") ) %>%
  anti_join(de_prot_for_camera_output$excluded_duplicates, by=c("uniprot_acc" = "uniprot_acc",
                                                                "comparison" = "comparison") )  %>%
  left_join( uniprot_acc_to_entrez_id %>%
               dplyr::filter(!is.na(ENTREZ_GENE)),
             by=c("uniprot_acc" = "UNIPROTKB") ) %>%
  dplyr::filter( !is.na(ENTREZ_GENE))  %>%
  dplyr::relocate(ENTREZ_GENE, "uniprot_acc" ) %>%
  dplyr::select(-uniprot_acc) %>%
  mutate( ENTREZ_GENE  =trimws( ENTREZ_GENE))


norm_abundance_mat_list <- norm_abundance_helper %>%
   group_by(comparison) %>%
   nest() %>%
  ungroup %>%
  mutate( data =  purrr::map(data,  function(x) {  as.matrix(column_to_rownames(x, "ENTREZ_GENE")) }))

## norm_abundance_mat_list$data[[1]]



## Prepare gene ID to gene set dictionary

list_of_collections <- purrr::map(msigdb.EZID, ~bcCategory(collectionType(.)) )

list_of_sub_collections <- purrr::map(msigdb.EZID, ~bcSubCategory(collectionType(.)) )


collections_tab <- data.frame( collection =unlist(list_of_collections),
                               sub_collection = unlist( list_of_sub_collections ) ) %>%
  distinct() %>%
  dplyr::filter( collection !="archived")

gene_set_idx_list_file <- file.path(output_dir,  "gene_set_idx_list_file.RDS")

if(file.exists(gene_set_idx_list_file)) {
  gene_set_idx_list <- readRDS(gene_set_idx_list_file)

} else {


  gene_set_idx_list <- purrr::pmap(collections_tab, function( collection, sub_collection) {

    if( is.na(sub_collection )) {
      subsetCollection(msigdb.EZID,
                       collection=collection )

    } else {
      subsetCollection(msigdb.EZID,
                       collection=collection,
                       subcollection=sub_collection)

    }


  } )


  names( gene_set_idx_list) <- collections_tab %>%
    mutate( all_collection = paste( collection, sub_collection, sep=", ")) %>%
    pull(all_collection)

  saveRDS(gene_set_idx_list, gene_set_idx_list_file)

}




## Example usage

abundance_mat_trimmed <- norm_abundance_mat_list[1,"data"][[1]][[1]]

msigdb_ids <- geneIds(gene_set_idx_list[["h, NA"]])

#convert gene sets into a list of gene indices
camera_indices <- ids2indices(msigdb_ids,
                              rownames(abundance_mat_trimmed))


camera( y=abundance_mat_trimmed,
        index=camera_indices,
        design =ruvIII_replicates_matrix,
        contrast=c(1, -1, 0, 0 ))


## List of contrasts


# colnames(ruvIII_replicates_matrix)
#
# "RPE_B" "RPE_A" "RPE_X" "RPE_Y"
#
# contrasts_tbl %>%
#   separate( contrasts, sep="=", into=c("contrast", "calculation")) %>%
#   pull(contrast)

contrasts_tbl <- vroom::vroom(contrasts_file, delim = "\t")
contrast_strings <- contrasts_tbl[, 1][[1]]
  ## Make contrasts
  contr.matrix <- makeContrasts(contrasts = contrast_strings,
                                levels = colnames(design_m))

lists_of_contrasts <-   map(colnames(contr.matrix), ~contr.matrix[,.] )
names(lists_of_contrasts) <-  colnames(contr.matrix) %>% str_split( "=") %>% purrr::map_chr(1)

# lists_of_contrasts <- list( RPE_Y.vs.RPE_X=c(0, 0, -1, 1 ) ,
#                             RPE_B.vs.RPE_A=c(1, -1, 0, 0 ) ,
#                             RPE_B.vs.RPE_X=c(1, 0, -1, 0 ) ,
#                             RPE_Y.vs.RPE_A=c(0, -1, 0, 1 ) ,
#                             RPE_A.vs.RPE_X=c(0, 1, -1, 0 ) ,
#                             RPE_B.vs.RPE_Y=c(1, 0, 0, -1 )   )







contrast_name <- "SS.vs.WT"
index_name <- "c5, GO:BP"
abundance_mat = norm_abundance_mat_list
replicates_mat = ruvIII_replicates_matrix
lists_of_contrasts = lists_of_contrasts
list_of_gene_sets = gene_set_idx_list
min_set_size=4

  print(paste("contrast_name =", contrast_name))
  print(paste("index_name =", index_name))


  groupA <- str_replace(contrast_name, "(.*)\\.vs\\.(.*)", "\\1")
  groupB <- str_replace(contrast_name, "(.*)\\.vs\\.(.*)", "\\2")

  design_choose_column <- replicates_mat[, c(groupA, groupB)]

  design_trimmed <- design_choose_column[rowSums(design_choose_column) > 0,]

  # print(design_trimmed)
  #
  # print(head( abundance_mat[[contrast_name]][ , rownames(design_trimmed)]))

  abundance_mat_trimmed <- abundance_mat %>%
    dplyr::filter( comparison == contrast_name) %>%
    dplyr::pull( data) %>%
    .[[1]] %>%
    .[,rownames(design_trimmed)]

  contrast_mat_trimmed <- lists_of_contrasts[[contrast_name]][colnames(replicates_mat) %in% colnames(design_trimmed)]

  index <- list_of_gene_sets[[index_name]]

  msigdb_ids <- geneIds(index)

  #convert gene sets into a list of gene indices
  camera_indices <- ids2indices(msigdb_ids,
                                rownames(abundance_mat_trimmed))

  ## At least two genes in the gene set
  camera_indices_filt <- camera_indices[purrr::map(camera_indices, length) >= min_set_size]


  camera_result <- NA
  if (length(camera_indices_filt) > 0) {
    camera_result <- camera(y = abundance_mat_trimmed, design = design_trimmed, index = camera_indices_filt, contrast = contrast_mat_trimmed)
  }

  info_list <- list(camera = camera_result,
                    y = abundance_mat_trimmed,
                    design = design_trimmed,

                    index_name = index_name,
                    index = camera_indices_filt,

                    contrast_name = contrast_name,
                    contrast = contrast_mat_trimmed)



## Run the camera test


    createDirectoryIfNotExists(output_dir)

    camera_results_file <- file.path(output_dir,  "camera_results.RDS")

if( file.exists(camera_results_file) ) {
    camera_results <- readRDS( camera_results_file )
} else {

    index_tab <- data.frame( index_name = names(gene_set_idx_list) ) %>%
        mutate( temp = 1 )


    contrast_tab <-    data.frame( contrast_name = names(lists_of_contrasts)) %>%
        mutate( temp = 1 )

    combination_tab <- contrast_tab %>%
        full_join( index_tab, by="temp") %>%
        dplyr::select(-temp)

    # Pre-fill the data, other data to be filled with the pmap function
    my_partial_camera <- partial( cmriCamera,
                                  abundance_mat = norm_abundance_mat_list ,
                                  replicates_mat = ruvIII_replicates_matrix,
                                  lists_of_contrasts = lists_of_contrasts,
                                  list_of_gene_sets = gene_set_idx_list,
                                  min_set_size=4 )

    camera_results <- purrr::pmap (combination_tab ,
                                   my_partial_camera)


    names(norm_abundance_mat_list)



    saveRDS( camera_results, camera_results_file)
}




## Covert the camera results that are stored in list structures into a table


camera_results_filt <- camera_results [purrr::map_lgl( camera_results, ~{!is.na(.[["camera"]][[1]][1])})]


camera_results_cln <-  purrr::map( camera_results_filt,
  ~ { rownames_to_column(.[["camera"]], "pathway")   }  )

camera_results_tbl <- tibble( temp = camera_results_cln)  %>%
    bind_cols(  data.frame( comparison =  purrr::map_chr( camera_results_filt,
   ~ {  .[["contrast_name"]]  } ) )  ) %>%
    bind_cols(  data.frame( gene_set =  purrr::map_chr( camera_results_filt,
   ~ {  .[["index_name"]]  } ) )  ) %>%
    unnest(temp)

camera_results_filt <- camera_results_tbl %>%
    dplyr::filter( FDR <0.05)


camera_results_tbl %>%
    arrange(FDR)


vroom::vroom_write( camera_results_filt,
                    file.path( output_dir,  "filtered_camera_results.tab"))


writexl::write_xlsx( camera_results_filt,
                     file.path(output_dir,  "filtered_camera_results.xlsx"))



