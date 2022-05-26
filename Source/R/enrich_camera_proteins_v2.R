#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


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
p_load(furrr)

## Parameters
tic()


## Directories management
base_dir <- here::here()

# output_dir <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/enrich_proteins"
# args$group_pattern <- "\\d+"
# de_proteins_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/de_proteins/de_proteins_long.tsv"
# contrasts_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Source/N155S/contrast_strings.tab"
# design_matrix_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Source/N155S/design_matrix.tab"
# counts_table_file <- "/home/ignatius/PostDoc/2022/PYROXD1_BMP_13/Results/N155S/de_proteins/normalized_counts_after_ruv.tsv"
# formula_string <- "~ 0 + group "
# sample_id <- "Sample_ID"
# group_id <- "group"


command_line_options <- commandArgs(trailingOnly = TRUE)

  parser <- OptionParser(add_help_option =TRUE)

  parser <- add_option(parser, c("-c", "--config"), type = "character", default =  "/home/ubuntu/Workings/2022/EStim_Brain_Organoids_BMP_17/Source/config_prot.ini", # "/home/ignatius/PostDoc/2021/ALPK1_BMP_06/Source/P90/config_prot_NR.ini",
                       help = "Configuration file.  [default %default]",
                       metavar = "string")

  parser <- add_option(parser, c("-l", "--log_file"), type = "character", default = "output.log",
                       help = "Name of the logging file.  [default %default]",
                       metavar = "string")

  parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                       help = "Deactivate backup of previous run.  [default %default]")

  parser <- add_option(parser, c("-t","--tmp_dir"), type = "character", default = "cache", dest = "tmp_dir",
                       help = "Directory path for temporary files.",
                       metavar = "string")

  parser <- add_option(parser, c( "--group-pattern"), type="character", dest = "group_pattern",
                       help="Regular expression pattern to identify columns with abundance values belonging to the experiment. [default %default]",
                       metavar="string")

  parser <- add_option(parser, c( "--counts"), type="character", dest = "counts_table_file",
                       help="Input file with the protein abundance values",
                       metavar="string")

  parser <- add_option(parser, c( "--de-proteins"), type="character", dest = "de_proteins_file",
                       help="Input file with the list of diffierentiall expressed protein log fold-change and q-values for every contrasts.",
                       metavar="string")

  parser <- add_option(parser, c("--contrasts"), type="character", dest = "contrasts_file",
                       help="Input file with a table listing all comparisons to be made in string, one comparison per line (e.g. groupB.vs.group_A = groupB - groupA).",
                       metavar="string")

  parser <- add_option(parser, c("--formula"), type="character", dest = "formula_string",
                       help="A string representing the formula for input into the model.frame function. (e.g. ~ 0 + group).",
                       metavar="string")

  parser <- add_option(parser, c( "--design-matrix"), type="character", dest = "design_matrix_file",
                       help="Input file with the design matrix",
                       metavar="string")

  parser <- add_option(parser, c("-o", "--output-dir"), type="character", dest = "output_dir",
                       help="Directory path for all results files.",
                       metavar="string")

  parser <- add_option(parser, c("--sample-id"), type="character",dest = "sample_id",
                       help="A string describing the sample ID. This must be a column that exists in the design matrix.",
                       metavar="string")

  parser <- add_option(parser, c( "--group-id"), type="character", dest = "group_id",
                       help="A string describing the experimental group ID. This must be a column that exists in the design matrix.",
                       metavar="string")

  parser <- add_option(parser, c("--dictionary_file"), type="character", dest = "dictionary_file",
                       help="File that converts annotation ID to annotation name.",
                       metavar="string")

  parser <- add_option(parser, c("--annotation_file"), type="character", dest = "annotation_file",
                       help="File with protein accession and functional annotation data.",
                       metavar="string")

  parser <- add_option(parser, c("--protein_id"), type="character", dest = "protein_id",
                       help="The protein id in the annotation file.",
                       metavar="string")

  parser <- add_option(parser, c("--annotation_id"), type="character", dest = "annotation_id",
                       help="The annotatio id in the annotation file.",
                       metavar="string")

  parser <- add_option(parser, c("--annotation_column"), type="character", dest = "annotation_column",
                       help="Column with the long name and biological details of the functional annotation.",
                       metavar="string")

  parser <- add_option(parser, c("--annotation_type"), type="character", dest = "annotation_type",
                       help="Name of the functional annotation data, added to the output file name.",
                       metavar="string")

  parser <- add_option(parser, c("--aspect_column"), type="character", dest = "aspect_column",
                       help="The aspect of the GO term.",
                       metavar="string")

  parser <- add_option(parser, c("--max_gene_set_size"), type="character", dest = "max_gene_set_size",
                       help="The maximum number of genes associaed with each gene set.",
                       metavar="string")

  parser <- add_option(parser, c("--min_gene_set_size"), type="character", dest = "min_gene_set_size",
                       help="The minimum number of genes associaed with each gene set.",
                       metavar="string")

  parser <- add_option(parser, "--p_val_thresh", type = "double",
                       help = "Adjusted p-value threshold below which a GO term is significantly enriched",
                       metavar = "double")

  parser <- add_option(parser, "--protein_q_val_thresh", type = "double",
                       help = "q-value threshold below which a GO term is significantly enriched",
                       metavar = "double")

  parser <- add_option(parser, "--log_fc_column_name", type="character",
                       help="Column name in the input file that contains the log fold-change values of the proteins.",
                       metavar="string")

  parser <- add_option(parser, "--fdr_column_name", type="character",
                       help="Column name in the input file that contains the false discovery rate values of the proteins.",
                       metavar="string")

  parser <- add_option(parser, "--uniprot_to_gene_symbol_file", type = "character",
                       help = "A file that contains the dictionary to convert uniprot accession to gene symbol. Uses the column specified in 'protein_id_lookup_column' flag for protein ID. Uses the column specified in 'gene_symbol_column' for the gene symbol column.",
                       metavar = "string")

  parser <- add_option(parser, "--protein_id_lookup_column", type = "character",
                       help = "The  name of the column that contained the protein ID to convert into gene symobl.",
                       metavar = "string")

  parser <- add_option(parser, "--gene_symbol_column", type = "character",
                       help = "The  name of the column that contained the gene symobl.",
                       metavar = "string")

  parser <- add_option(parser, c("--num_cores"), type="integer", dest = "num_cores",
                       help="The number of CPU cores to used for computation.",
                       metavar="string")


  print(commandArgs(trailingOnly = TRUE))

  args <- parse_args(parser)

  #parse and merge the configuration file options.
  if (args$config != "") {
    args <- config.list.merge(eval.config(file = args$config, config="enrich_camera_proteins"), args)
  }

  args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="enrich_camera_proteins" )

  args <- setArgsDefault(args, "max_gene_set_size", as_func=as.character, default_val="200" )
  args <- setArgsDefault(args, "min_gene_set_size", as_func=as.character, default_val="4" )
  args <- setArgsDefault(args, "protein_id", as_func=as.character, default_val="uniprot_acc" )
  args <- setArgsDefault(args, "annotation_id", as_func=as.character, default_val="go_id" )
  args <- setArgsDefault(args, "aspect_column", as_func=as.character, default_val=NULL )
  if( isArgumentDefined(args, "aspect_column" )) {
    if( args$aspect_column == "NULL" ) {
      args$aspect_column <- NULL
    }
  }
  args <- setArgsDefault(args, "annotation_column", as_func=as.character, default_val="term" )
  args <- setArgsDefault(args, "annotation_type", as_func=as.character, default_val="annotation" )
  args <- setArgsDefault(args, "p_val_thresh", as_func=as.double , default_val=0.05 )
  args <- setArgsDefault(args, "protein_q_val_thresh", as_func=as.double, default_val=0.05 )
  args <- setArgsDefault(args, "log_fc_column_name", as_func=as.character, default_val="log2FC" )
  args <- setArgsDefault(args, "fdr_column_name", as_func=as.character, default_val="q.mod" )

  if( isArgumentDefined( args, "annotation_file" )) {

    testRequiredArguments(args, c(
      "annotation_file",
      "annotation_type",
      "dictionary_file",
      "annotation_column"
    ))

    testRequiredFiles(c(
      args$annotation_file,
      args$dictionary_file))
  }

  if(isArgumentDefined(args, "uniprot_to_gene_symbol_file")) {
    testRequiredFiles(c(
      args$uniprot_to_gene_symbol_file))

    args <- setArgsDefault( args,
                            "protein_id_lookup_column",
                            as_func=as.character,
                            default_val="Entry" )
    args <- setArgsDefault( args,
                            "gene_symbol_column",
                            as_func=as.character,
                            default_val="Gene names" )
  }

  testRequiredArguments(args, c(
    "formula_string",
    "group_pattern",
    "output_dir",
    "sample_id",
    "group_id"
  ))

  testRequiredFiles(c(
    args$design_matrix_file,
    args$counts_table_file,
    args$de_proteins_file,
    args$contrasts_file,
    args$config,
    args$tmp_dir
  ))

  args<-parseString(args,
                    c(  "formula_string",
                        "group_pattern",
                        "output_dir",
                        "sample_id",
                        "group_id",
                        "tmp_dir",
                        "dictionary_file",
                        "annotation_file",
                        "protein_id",
                        "annotation_id",
                        "annotation_column",
                        "annotation_type",
                        "aspect_column",
                        "max_gene_set_size",
                        "min_gene_set_size",
                        "log_fc_column_name",
                        "fdr_column_name" ))

  args <- parseType( args, c("num_cores"),
                     as.integer)

  args <- parseType( args, c("p_val_thresh",
                             "protein_q_val_thresh"),
                     as.double)



  #############################################
  createOutputDir(args$output_dir, args$no_backup)
  createDirectoryIfNotExists(args$tmp_dir)

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

####################

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Create directories
createDirectoryIfNotExists( args$output_dir)


  ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  loginfo("Read design matrix file.")

  captured_output <- capture.output(
    design_mat_cln <- vroom::vroom(args$design_matrix_file) %>%
      as.data.frame() %>%
      dplyr::mutate(!!rlang::sym(args$sample_id) := as.character(!!rlang::sym(args$sample_id)))
    ,type = "message"
  )
  logdebug(captured_output)

  rownames(design_mat_cln) <- design_mat_cln %>% pull(as.name(args$sample_id))

  ## Check that the sample ID and group ID name is found in the design matrix
  if (length(which(c(args$sample_id, args$group_id) %in% colnames(design_mat_cln))) != 2) {
    logerror("Sample ID and group ID are not matching to the column names used in the design matrix.")
    q()
  }

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Read the Contrasts file
contrasts_tbl <- NA
if( args$contrasts_file != "") {
  loginfo("Read file with lists of experimental contrasts to test %s", args$contrasts_file)
  captured_output<-capture.output(
    contrasts_tbl <- vroom::vroom( args$contrasts_file, delim="\t"),
    type = "message"
  )
  logdebug(captured_output)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## List of read abundance tables after normalization with RUVIII
loginfo("Read abundance tables after normalization with RUVIII.")
captured_output <- capture.output(
norm_abundance <- vroom::vroom( args$counts_table_file )  %>%
  mutate( best_uniprot_acc = purrr::map_chr( uniprot_acc, ~str_split(. , ":"   )[[1]][1] )),
type = "message"
)
logdebug(captured_output)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Prepare Normalized Abundance Matrix.")

norm_abundance_mat <- norm_abundance %>%
  dplyr::select(best_uniprot_acc, matches( args$group_pattern ) ) %>%
  column_to_rownames("best_uniprot_acc")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Functions to create the replicate matrix
loginfo("Create the replicate matrix.")

ruvIII_replicates_matrix <- getRuvIIIReplicateMatrix( design_mat_cln,
                                                         !!rlang::sym(args$sample_id),
                                                         !!rlang::sym(args$group_id))

ruvIII_replicates_matrix


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Set up list of contrasts

loginfo("Create the list of contrasts.")

ff <- as.formula(  args$formula_string)
mod_frame <- model.frame( ff, design_mat_cln)
design_m <- model.matrix( ff, mod_frame)

contr.matrix <- makeContrasts( contrasts = contrasts_tbl %>% pull(contrasts),
                                 levels = colnames(design_m))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## List of all proteins

loginfo("Get the list of all DE proteins.")
captured_output <- capture.output(
proteins_cln <- vroom::vroom(args$de_proteins_file) %>%
  mutate( best_uniprot_acc = purrr::map_chr( uniprot_acc, ~str_split(. , ":"   )[[1]][1] )) %>%
  dplyr::select( -uniprot_acc ) %>%
  dplyr::rename( uniprot_acc = "best_uniprot_acc")
,type = "message"
)
logdebug(captured_output)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Prepare gene ID to gene set dictionary.")

## Tidy up the annotation ID to annotation term name dictionary

dictionary <- vroom::vroom( args$dictionary_file )

id_to_annotation_dictionary <- buildAnnotationIdToAnnotationNameDictionary( input_table=dictionary,
                                                                            annotation_column = !!rlang::sym(args$annotation_column),
                                                                            annotation_id_column = !!rlang::sym(args$annotation_id))


annotation_types_list <- NA
if( !is.null( args$aspect_column )){

  if( args$aspect_column  %in% colnames( dictionary) ) {

    annotation_types_list <- listifyTableByColumn( dictionary,
                                                   !!rlang::sym(args$aspect_column))

  } else {
    annotation_types_list <- list( dictionary )
    names( annotation_types_list) <- args$annotation_type
  }

} else {
  annotation_types_list <- list( dictionary )
  names( annotation_types_list) <- args$annotation_type
}

annotation_gene_set_list <- purrr::map( annotation_types_list,
                                        ~{ buildOneProteinToAnnotationList(.,
                                                                           !!rlang::sym(args$annotation_id),
                                                                           !!rlang::sym(args$protein_id))})

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Gene set size maximum and minimum
min_gene_set_size_list <- parseNumList(args$min_gene_set_size)
max_gene_set_size_list <- parseNumList(args$max_gene_set_size)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## List of contrasts
loginfo("Prepare list of contrasts.")

contrasts_tbl <- vroom::vroom(args$contrasts_file, delim = "\t")
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

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Run the camera test.")

createDirectoryIfNotExists(args$output_dir)

camera_results_file <- file.path(args$output_dir,  "camera_results.RDS")

if( file.exists(camera_results_file) ) {
  loginfo("Read file with camera test results %s", camera_results_file)
  captured_output<-capture.output(
    camera_results <- readRDS( camera_results_file ),
    type = "message"
  )
  logdebug(captured_output)
} else {

    combination_tab <- cross_df( list( index_name = names( annotation_gene_set_list ),
                                    contrast_name = names( lists_of_contrasts),
                                    min_set_size = min_gene_set_size_list,
                                    max_set_size = max_gene_set_size_list )) %>%
      mutate( contrast= purrr::map( contrast_name, ~lists_of_contrasts[[.]])  ) %>%
      mutate( index= purrr::map(index_name,  ~annotation_gene_set_list[[.]]  ) )

    # combination_tab[6,"contrast"][[1]][[1]]
    # combination_tab[6,"index"][[1]][[1]]

    # Pre-fill the data, other data to be filled with the pmap function
    my_partial_camera <- partial( cmriCamera,
                                  abundance_mat = norm_abundance_mat ,
                                  design_mat = design_m)

    # camera( norm_abundance_mat, index = annotation_gene_set_list[[1]],
    #         design = design_m, contrast = lists_of_contrasts[[1]])

    plan(multisession, workers = args$num_cores)

    camera_results <- furrr::future_pmap (combination_tab,
                                   function(index_name,
                                            contrast_name,
                                            index,
                                            contrast,
                                            min_set_size,
                                            max_set_size) {
                                     my_partial_camera( index_name=index_name,
                                                        contrast_name=contrast_name,
                                                        index=index,
                                                        contrast=contrast,
                                                        min_set_size=min_set_size,
                                                        max_set_size=max_set_size ) })

    loginfo("Save file with camera test results %s", camera_results_file)
    captured_output<-capture.output(
      saveRDS( camera_results, camera_results_file),
      type = "message" )
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Covert the camera results that are stored in list structures into a table.")

camera_results_temp <- camera_results [purrr::map_lgl( camera_results, ~{!is.na(.[["camera"]][[1]][1])})]

camera_results_cln <-  purrr::map( camera_results_temp,
  ~ { rownames_to_column(.[["camera"]], args$annotation_id)   }  )

camera_results_tbl <- tibble( temp = camera_results_cln)  %>%
  bind_cols(  data.frame( comparison = purrr::map_chr( camera_results_temp,
                                                       ~ {  .[["contrast_name"]]  } ) )  ) %>%
  bind_cols(  data.frame( gene_set = purrr::map_chr( camera_results_temp,
                                                     ~ {  .[["index_name"]]  } ) )  ) %>%
  bind_cols(  data.frame( min_set_size = purrr::map_chr( camera_results_temp,
                                                         ~ {  .[["min_set_size"]]  } ) )  ) %>%
  bind_cols(  data.frame( max_set_size = purrr::map_chr( camera_results_temp,
                                                         ~ {  .[["max_set_size"]]  } ) )  ) %>%
  unnest(temp) %>%
  mutate( term = purrr::map_chr( !!rlang::sym(args$annotation_id),  ~id_to_annotation_dictionary[[.]]))

if( !is.null( args$aspect_column)) {
  camera_results_tbl <- camera_results_tbl %>%
    dplyr::rename( !!rlang::sym(args$aspect_column) := "gene_set")
}

rm( camera_results)
gc()

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Add UniProt accession of proteins that had significant log fold-change and also associated with annotation.")

join_condition <- NA
columns_included  <- NA
if( !is.null( args$aspect_column) ) {

  join_condition <- rlang::set_names( c(args$annotation_id, args$aspect_column),
                                      c(args$annotation_id, args$aspect_column ) )

  columns_included <- c(args$annotation_id, args$aspect_column, "uniprot_acc")

} else {
  join_condition <- rlang::set_names( c(args$annotation_id),
                                      c(args$annotation_id ) )

  columns_included <- c(args$annotation_id, "uniprot_acc")

}

join_condition_two <- rlang::set_names( c("comparison", "uniprot_acc"),
                                        c("comparison", "uniprot_acc") )


camera_results_with_uniprot_acc <- camera_results_tbl %>%
  mutate(  !!rlang::sym( args$annotation_id) := as.character(  !!rlang::sym( args$annotation_id))) %>%
  left_join( dictionary %>%
               dplyr::rename( uniprot_acc =  args$protein_id) %>%
               mutate(  !!rlang::sym( args$annotation_id) := as.character(  !!rlang::sym( args$annotation_id))) %>%
               dplyr::select( one_of(columns_included)),
             by=join_condition ) %>%
  left_join( proteins_cln  %>%
             dplyr::filter( !!rlang::sym( args$fdr_column_name) < args$protein_q_val_thresh) %>%
             dplyr::select( !!rlang::sym( args$fdr_column_name),
                               !!rlang::sym( args$log_fc_column_name),
                               uniprot_acc,
                               comparison),
              by = join_condition_two ) %>%
  dplyr::filter ( is.na(!!rlang::sym( args$log_fc_column_name)) |
                  ( Direction == "Down" & !!rlang::sym( args$log_fc_column_name) < 0) |
                  ( Direction == "Up" & !!rlang::sym( args$log_fc_column_name) > 0) ) %>%
  dplyr::select( -one_of(c( args$fdr_column_name,
                           args$log_fc_column_name)))


rm( camera_results_tbl)
gc()

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Add gene symbol for the proteins that had significant log fold-change and also associated with annotation.")

camera_results_with_gene_symbol <- NA
## Convert Uniprot accession to gene names
if(isArgumentDefined(args, "uniprot_to_gene_symbol_file")) {
  # args$uniprot_to_gene_symbol_file <- "/home/ubuntu/Workings/2021/ALPK1_BMP_06/Data/UniProt/data.tab"

  # args$protein_id_lookup_column <- "Entry"
  # args$gene_symbol_column <- "Gene names"

  # Clean up protein ID to gene sybmol table

  uniprot_tab_delimited_tbl <- vroom::vroom( file.path( args$uniprot_to_gene_symbol_file))

  uniprot_to_gene_names <- uniprot_tab_delimited_tbl %>%
    dplyr::mutate( gene_symbol = str_split(!!rlang::sym(args$gene_symbol_column), " |;|/")) %>%
    unnest(gene_symbol) %>%
    distinct( !!rlang::sym(args$protein_id_lookup_column), gene_symbol) %>%
    dplyr::filter( gene_symbol != "" & !is.na(gene_symbol)) %>%
    dplyr::rename( uniprot_acc = args$protein_id_lookup_column)

  uniprot_acc_to_gene_symbol_join_condition <- rlang::set_names(c("uniprot_acc"),
                                                                c("uniprot_acc"))

  camera_results_with_gene_symbol <- camera_results_with_uniprot_acc %>%
    left_join( uniprot_to_gene_names, by=uniprot_acc_to_gene_symbol_join_condition)

  # uniprot_to_gene_symbol_dict <- getUniprotAccToGeneSymbolDictionary( uniprot_tab_delimited_tbl,
  #                                      !!rlang::sym(args$protein_id_lookup_column),
  #                                      !!rlang::sym(args$gene_symbol_column),
  #                                      !!rlang::sym(args$protein_id) )
  #
  # camera_results_with_gene_symbol <- camera_results_with_uniprot_acc %>%
  #   mutate( gene_symbol = furrr::future_map_chr(  !!rlang::sym(args$protein_id),
  #                                       function(x){ ifelse(!is.na(x) & (x %in% names( uniprot_to_gene_symbol_dict)),
  #                                                           uniprot_to_gene_symbol_dict[[x]],
  #                                                           NA_character_) }))

} else {
  camera_results_with_gene_symbol <- camera_results_with_uniprot_acc
}



rm( camera_results_with_uniprot_acc)
gc()

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Convert gene symbols and UniProt accession into comma separated spring.")

camera_results_unfilt <- camera_results_with_gene_symbol %>%
  group_by_at( vars( !contains("uniprot_acc") & !contains("gene_symbol")) ) %>%
  nest() %>%
  ungroup() %>%
  mutate( accession_list = furrr::future_map_chr( data, function(x) { x %>%
      arrange (uniprot_acc) %>%
      pull( uniprot_acc ) %>%
      sort %>%
      paste( collapse=", ")    } )) %>%
  mutate( gene_symbol = furrr::future_map_chr( data, function(x) { x %>%
      arrange(gene_symbol) %>%
      pull( gene_symbol) %>%
      sort %>%
      paste( collapse=", ")    } ) ) %>%
  dplyr::select(-data)

camera_results_filt <- camera_results_unfilt %>%
  dplyr::filter( FDR < args$p_val_thresh)

# camera_results_tbl %>%
#     arrange(FDR)
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


loginfo("Save camera results table in tab-separated table %s", file.path( args$output_dir,  "filtered_camera_results.tab"))
captured_output<-capture.output(
  vroom::vroom_write( camera_results_filt,
                      file.path( args$output_dir,  "filtered_camera_results.tab")),
  type = "message" )

loginfo("Save camera results table in Excel format %s", file.path( args$output_dir,  "filtered_camera_results.xlsx"))
captured_output<-capture.output(
writexl::write_xlsx( camera_results_filt,
                     file.path( args$output_dir,  "filtered_camera_results.xlsx")),
type = "message" )





loginfo("Save camera results table in tab-separated table %s", file.path( args$output_dir,  "unfiltered_camera_results.tab"))
captured_output<-capture.output(
  vroom::vroom_write( camera_results_unfilt,
                      file.path( args$output_dir,  "unfiltered_camera_results.tab")),
  type = "message" )

loginfo("Save camera results table in Excel format %s", file.path( args$output_dir,  "unfiltered_camera_results.xlsx"))
captured_output<-capture.output(
  writexl::write_xlsx( camera_results_unfilt,
                       file.path( args$output_dir,  "unfiltered_camera_results.xlsx")),
  type = "message" )



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))


