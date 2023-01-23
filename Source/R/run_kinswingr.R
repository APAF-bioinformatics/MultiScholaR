#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases

## ---------------------------------------------------------------------------------------------------------------------------------------------------

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
   BiocManager::install(version = "3.12")
}

# load pacman package manager
if(!require(pacman)){
    install.packages("pacman")
    library(pacman)
}

p_load(tidyverse)
p_load(vroom)
p_load(BiocParallel)
p_load(magrittr)
p_load(KinSwingR)
p_load(RColorBrewer)
p_load(ggrepel)

p_load(optparse)
p_load(configr)
p_load(logging)
p_load(svglite)
p_load(tictoc)
p_load(ProteomeRiver)
p_load(UniProt.ws)
p_load(furrr)

tic()

## ---------------------------------------------------------------------------------------------------------------------------------------------------

command_line_options <- commandArgs(trailingOnly = TRUE)
#Note: options with default values are ignored in the configuration file parsing.

parser <- OptionParser(add_help_option = TRUE)

#Note: options with default values are ignored in the configuration file parsing.
parser <- add_option(parser, c("-d", "--debug"), action = "store_true", default = FALSE,
                     help = "Print debugging output")

parser <- add_option(parser, c("-s", "--silent"), action = "store_true", default = FALSE,
                     help = "Only print critical information to the console.")

parser <- add_option(parser, c("-c", "--config"), type = "character", default = "config_phos.ini",
                     help = "Configuration file.",
                     metavar = "string")

parser <- add_option(parser, c("-t","--tmp_dir"), type = "character",  dest = "tmp_dir",
                     help = "Directory path for temporary files.",
                     metavar = "string")

parser <- add_option(parser, c( "--num_cores"), type = "integer",
                     help = "The number of cores used for the computation.",
                     metavar = "integer")

parser <- add_option(parser, c( "--random_seed"), type = "integer",
                     help = "The number of cores used for the computation.",
                     metavar = "integer")

parser <- add_option(parser, c( "--motif_score_iteration"), type = "integer",
                     help = "The number iterations for scoring each substrate against kinase motifs.",
                     metavar = "integer")

parser <- add_option(parser, c( "--swing_iteration"), type = "integer",
                     help = "The number iterations for swing score calcualtions.",
                     metavar = "integer")

parser <- add_option(parser, c( "--min_num_sites_per_kinase"), type = "integer",
                     help = "The minimum number of known substrates for each kinase",
                     metavar = "integer")

parser <- add_option(parser, c( "--p_value_cutoff"), type = "double",
                     help = "The p-value cutoff for identifying a kinase as significant by KinSwingR",
                     metavar = "integer")

parser <- add_option(parser, c( "--ggrepel_p_value_cutoff"), type = "double",
                     help = "The p-value cutoff for showing the kinase name in the volcano plot",
                     metavar = "integer")

parser <- add_option(parser, c( "--kinase_specificity"), type = "character",
                     help = "Specificity of the kinase. One of Ser/Thr kinase= ST, Tyr kinase = Y, or Ser/Thr/Tyr kinase = STY",
                     metavar = "integer")

parser <- add_option(parser, c( "--phosphosite_db_dir"), type = "character",
                     help = "The directory in which PhosphositePlus data is stored",
                     metavar = "character")

parser <- add_option(parser, c( "--uniprot_kinase_file"), type = "character",
                     help = "The path to the UniProt tab-separated data file.",
                     metavar = "character")

parser <- add_option(parser, "--norm_phos_logfc_file", type="character",
                     help="Results table in which the phosphorylation log fold-change is normalized by protein log fold-change.",
                     metavar="string")

parser <- add_option(parser, "--uniprot_other_kinase_file", type="character",
                     help="File listing the UniProt accession of atypical and other kinases and their UniProt Keywords.",
                     metavar="string")

parser <- add_option(parser, "--log_fc_column_name", type="character",
                     help="Column name in the input file that contains the log fold-change values of the phosphosites.",
                     metavar="string")

parser <- add_option(parser,  "--taxonomy_id", type="integer", dest = "taxonomy_id",
                     help="The NCBI taxonomy ID of the organism being investigated (e.g. M. musculus=10090, H. sapien=9606).",
                     metavar="integer")

parser <- add_option(parser, "--fdr_column_name", type="character",
                     help="Column name in the input file that contains the false discovery rate values of the phosphosites.",
                     metavar="string")

parser <- add_option(parser, c("-o", "--output_dir"), type = "character", dest = "output_dir",
                     help = "Directory path for all results files.",
                     metavar = "string")

parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                     help = "Deactivate backup of previous run.")

parser <- add_option(parser, c("-l", "--log_file"), type = "character", default = "output.log",
                     help = "Name of the logging file.",
                     metavar = "string")

parser <- add_option(parser, "--plots_format", type = "character",
                     help = "A comma separated strings to indicate the fortmat output for the plots [pdf,png,svg].",
                     metavar = "string")

parser <- add_option(parser, "--is_multisite", type = "logical",
                     help = "If NA, include both single-site and multi-site phosphorylation. If TRUE, include only multi-site phosphorylation. If FALSE, include only single-site phosphorylation.",
                     metavar = "string")

parser <- add_option(parser, "--reuse_old", type = "logical",
                     help = "If TRUE, reuse result from files saved from previous run.")

#parse command line arguments first.
args <- parse_args(parser)


# knitr::purl(input="/home/ignatius/PostDoc/2021/proteomeriver/Source/Rmd/run_kinswingr.Rmd",
#             output="/home/ignatius/PostDoc/2021/proteomeriver/Source/R/run_kinswingr.R" )

## ----------------------------------------------------------------------------------------------------------------------------------------
# args$debug <- TRUE
# if (args$debug == TRUE) {
#
#   args$tmp_dir <- "/home/ignatius/PostDoc/2021/Podocyte_2021/Results/Part_1/cache"
#   args <- config.list.merge(eval.config(file = args$config, config = "phos_kinase_substrate"), args)
#
# }

## ----------------------------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config = "phos_kinase_substrate"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="phos_kinswingr" )


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

testRequiredArguments(args, c(
  "kinase_specificity",
  "phosphosite_db_dir",
  "uniprot_kinase_file",
  "norm_phos_logfc_file",
  "uniprot_other_kinase_file",
  "output_dir",
  "p_value_cutoff",
  "ggrepel_p_value_cutoff",
  "plots_format",
  "taxonomy_id"
))


ks_file <-  file.path( args$phosphosite_db_dir, "Kinase_Substrate_Dataset")

testRequiredFiles(c( args$norm_phos_logfc_file,
                     args$uniprot_kinase_file
))

## Set default values

if(isArgumentDefined(args,"plots_format"))
{
  args <- parseList(args,
                    c("plots_format"))
}else {
  logwarn("plots_format is undefined, default output set to pdf.")
  args$plots_format <- list("pdf")
}

args <- setArgsDefault(args, "random_seed", as_func=as.integer, default_val=123456 )
args <- setArgsDefault(args, "motif_score_iteration", as_func=as.integer, default_val=1000 )
args <- setArgsDefault(args, "swing_iteration", as_func=as.integer, default_val=1000 )
args <- setArgsDefault(args, "min_num_sites_per_kinase", as_func=as.integer, default_val=10 )
args <- setArgsDefault(args, "num_cores", as_func=as.integer, default_val=1 )
args <- setArgsDefault(args, "log_fc_column_name", as_func=as.character, default_val="norm_phos_logFC" )
args <- setArgsDefault(args, "fdr_column_name", as_func=as.character, default_val="combined_q_mod" )
args <- setArgsDefault(args, "p_value_cutoff", as_func=as.double, default_val=0.05 )
args <- setArgsDefault(args, "is_multisite", as_func=as.logical, default_val=NA )
args <- setArgsDefault(args, "tmp_dir", as_func=as.character, default_val="cache" )
args <- setArgsDefault(args, "reuse_old", as_func=as.logical, default_val=FALSE )


createDirectoryIfNotExists(args$tmp_dir)


## ----------------------------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Read phosphosite log fold-change normalized by protein log fold-change file")
captured_output<-capture.output(
  de_phos <- vroom::vroom( args$norm_phos_logfc_file )
  ,type = "message"
)
logdebug(captured_output)

## ---------------------------------------------------------------------------------------------------------------------------------------------------

# uniprot_kinase_file<-file.path(args$tmp_dir,args$uniprot_kinase_file)
# if(!file.exists(uniprot_kinase_file))
# {
#   logwarn("Download UniProt kinases list file.")
#   status <- download.file(url="https://www.uniprot.org/docs/pkinfam.txt", destfile=uniprot_kinase_file)
#   loginfo(status)
# }
loginfo("Reading UniProt kinases list file.")
captured_output<-capture.output(
  uniprot_kinase_tbl <- vroom::vroom(  args$uniprot_kinase_file   )
    , type = "message"
)
logdebug(captured_output)

## ---------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Download information from UniProt.")
uniprot_other_kinase_file<-file.path(args$tmp_dir,args$uniprot_other_kinase_file)
if( ! file.exists( uniprot_other_kinase_file )) {
  human_taxonomy_id <- 9606
  up <- UniProt.ws(taxId=human_taxonomy_id)
  list_of_sp_columns <- c("EXISTENCE"
                          , "SCORE"
                          , "REVIEWED"
                          , "GENENAME"
                          , "PROTEIN-NAMES"
                          , "LENGTH"
                          , "ENSEMBL"
                          , "GO-ID"
                          , "KEYWORDS"

                          ,"protein_existence"
                          ,"annotation_score"#?
                          ,"reviewed"
                          ,"gene_names"
                          ,"protein_name"
                          ,"length"
                          ,"xref_ensembl"
                          , "go_id"
                          , "keyword")

  up_cls<-unlist(columns(up))
  list_intersect<-intersect(list_of_sp_columns, up_cls)

  if(length(setdiff( list_of_sp_columns,list_intersect)) > 0) {
    logwarn("UniProt fields not found: %s",paste(setdiff( list_of_sp_columns,list_intersect),sep=", "))
  }

  my_keytype <- "UniProtKB"
  if( "UNIPROTKB" %in% keytypes(up) ) {
    my_keytype <- "UNIPROTKB"
  }


  uniprot_acc_tbl <- uniprot_kinase_tbl |>
    dplyr::filter( Family == "Other" |
                     str_detect(Family, "Atypical")) |>
    dplyr::distinct(uniprot_acc_human) |>
    dplyr::rename(uniprot_acc = uniprot_acc_human)

  atypical_and_other <- batchQueryEvidence(uniprot_acc_tbl, uniprot_acc_column="uniprot_acc", uniprot_handle=up,
                                           uniprot_columns = list_of_sp_columns, uniprot_keytype=my_keytype)

  saveRDS( atypical_and_other, file.path( uniprot_other_kinase_file))

}

loginfo("Reading UniProt data from %s.",uniprot_other_kinase_file)
atypical_and_other <- readRDS( uniprot_other_kinase_file )

## ----------------------------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Read PhosphoSitePlus (PSP) kinase-substrate table.")
captured_output<-capture.output(
  ks_tbl <- vroom::vroom( ks_file, skip=3 ) |>
  mutate( SUB_MOD_RSD_CLN = str_replace_all(SUB_MOD_RSD, "([A-Z])(\\d+)", "\\1 \\2")) |>
  separate( SUB_MOD_RSD_CLN, into=c("residue", "position"))
  ,type = "message"
)
logdebug(captured_output)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
# As an example of control over multi-core processing
# load BiocParallel library
# finally set/register the number of cores to use

loginfo("Set number of cores for parallel computation.")
register(SnowParam(workers = args$num_cores))

## ---------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Check log fold-change and FDR column names exists in input file.")


if ( ! args$log_fc_column_name %in% colnames(de_phos ) ) {
  logerror("Column '%s' is not found in the input table.",args$log_fc_column_name )
}

if ( ! args$fdr_column_name %in% colnames(de_phos ) ) {
  logerror("Column '%s' is not found in the input table.",args$fdr_column_name)
}

if(  ! args$log_fc_column_name %in% colnames(de_phos )  |
     ! args$fdr_column_name %in% colnames(de_phos )  ) {
  stop()
}

## ---------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Filter the correct subset of kinases for the analysis.")

uniprot_kinases <- NA

if( args$kinase_specificity == "ST") {
   uniprot_kinases <- uniprot_kinase_tbl |>
     dplyr::filter( str_detect( Family, "Ser/Thr" )  | str_detect( Family, "Atypical") | str_detect( Family, "Other" ) ) |>
     left_join( atypical_and_other |>
                  dplyr::select(UNIPROTKB, KEYWORDS),
                by=c("uniprot_acc_human" = "UNIPROTKB")) |>
      dplyr::filter( str_detect( Family, "Ser/Thr" ) |
                       (   str_detect(KEYWORDS, "Serine/threonine-protein kinase") &
                          (!str_detect(KEYWORDS, "Tyrosine-protein kinase") ) ) )

} else if ( args$kinase_specificity == "Y") {
  uniprot_kinases <- uniprot_kinase_tbl |>
    dplyr::filter( str_detect( Family, "Tyr" )  | str_detect( Family, "Atypical") | str_detect( Family, "Other" ) ) |>
    left_join( atypical_and_other |>
                  dplyr::select(UNIPROTKB, KEYWORDS),
               by=c("uniprot_acc_human" = "UNIPROTKB")) |>
    dplyr::filter( str_detect( Family, "Tyr" ) |
                     ( (!str_detect(KEYWORDS, "Serine/threonine-protein kinase")) &
                         str_detect(KEYWORDS, "Tyrosine-protein kinase") ) )


} else if ( args$kinase_specificity == "STY") {
  uniprot_kinases <- uniprot_kinase_tbl |>
    dplyr::filter(  str_detect( Family, "Atypical") | str_detect( Family, "Other" ) ) |>
    left_join( atypical_and_other |>
                  dplyr::select(UNIPROTKB, KEYWORDS),
               by=c("uniprot_acc_human" = "UNIPROTKB")) |>
    dplyr::filter( str_detect(KEYWORDS, "Serine/threonine-protein kinase") &
                         str_detect(KEYWORDS, "Tyrosine-protein kinase")  )
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------


phosphositeplus <- ks_tbl |>
  dplyr::select( GENE, KINASE, `SITE_+/-7_AA`, SUB_MOD_RSD) |>
  distinct() |>
  dplyr::rename( kinase = "KINASE",
                 substrate = "SITE_+/-7_AA") |>
  dplyr::mutate( substrate = toupper(substrate)) |>
  dplyr::mutate( kinase = str_replace( kinase, " ", "_")) |>
  arrange( GENE, kinase, substrate)  |>
  mutate( GENE = toupper(GENE),
          kinase = toupper(kinase),
          residue = purrr::map_chr( SUB_MOD_RSD , \(x){str_sub(., 1,1 )}   ) ) |>
  dplyr::select( - SUB_MOD_RSD )

kinases_to_include <-  phosphositeplus |>
  left_join( uniprot_kinases |>
               dplyr::select(gene_name) |>
               dplyr::mutate( is_uniprot_a= 1),
             by=c( "GENE" = "gene_name") ) |>
  left_join( uniprot_kinases |>
               dplyr::select(gene_name) |>
               dplyr::mutate( is_uniprot_b= 1),
             by=c( "kinase" = "gene_name") ) |>
  dplyr::filter( is_uniprot_a == 1 | is_uniprot_b == 1) |>
  dplyr::select( -is_uniprot_a, - is_uniprot_b ) |>
  dplyr::filter(  str_detect(  args$kinase_specificity, residue )  ) |>
  group_by(kinase )  |>
  summarise( counts =n()) |>
  ungroup() |>
  dplyr::filter( counts >= args$min_num_sites_per_kinase)

loginfo("Filter the correct subset of substrates for the analysis.")
phosphositeplus_filt <- phosphositeplus |>
  mutate( kinase = toupper(kinase) ) |>
  inner_join( kinases_to_include, by="kinase") |>
  dplyr::filter(  str_detect(  args$kinase_specificity, residue )  ) |>
  dplyr::select(-counts, - GENE, - residue) |>
  distinct() |>
  as.matrix()




## ---------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Build the position-specific scoring matrices (PSSM).")

pwms <- buildPWM(as.matrix(phosphositeplus_filt))

# pwms$kinase

## ---------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Clean up phoshopsite log fold-change table.")

annotated_data_pre_residue_filter <- de_phos |>
  mutate( peptide = sequence) |>
  mutate( positions_peptide = position   ) |>
  mutate(uniprot_acc = str_split( uniprot_acc, ":") |> purrr::map_chr(1) ) |>
  mutate(gene_name = str_split( gene_name, ":") |> purrr::map_chr(1) )  |>
  mutate(position = str_split( position, ":") |> purrr::map_chr(1) ) |>
  mutate( peptide =  str_split(peptide, ":"  )|> purrr::map_chr(1) ) |>
  mutate( residue =  str_split(residue, ":"  )|> purrr::map_chr(1) ) |>
  mutate( position = str_split( position, "\\|") |> purrr::map_chr(1) |> str_replace_all( "\\(|\\)", "") ) |>
  mutate ( is_multisite = case_when ( str_detect( position, ";") ~ TRUE,
                                   TRUE ~ FALSE)) |>
  separate_rows(peptide, position, residue, sep=";") |>
  dplyr::mutate( peptide_copy = peptide) |>
  dplyr::filter(  !str_detect( peptide, "X"))   |>
  dplyr::select(sites_id, comparison, uniprot_acc, gene_name, peptide, peptide_copy, position, residue,
                one_of( c( as.character(args$log_fc_column_name), as.character(args$fdr_column_name))), is_multisite) |>
  dplyr::filter(  str_detect(  args$kinase_specificity, residue )  )  |>
  dplyr::filter(!is.na(peptide)) |>
  dplyr::filter( str_sub( peptide, 8, 8) == residue) |>
  dplyr::rename( fc = args$log_fc_column_name,
                 pval = args$fdr_column_name)

# annotated_data_pre_residue_filter |> colnames()

----------------------------------------------------------------------------------------
loginfo("Filtering single-site or multisite phosphorylation")

if( !is.na(args$is_multisite)   ) {
  annotated_data_multisite_filtered <- annotated_data_pre_residue_filter |>
    dplyr::filter( is_multisite == args$is_multisite)

  if( args$is_multisite) {
    loginfo("Keeping only multi-site phosphorylation.")
  } else {
    loginfo("Keeping only single-site phosphorylation.")
  }

} else {
  loginfo("Keeping both single-site and multi-site phosphorylation.")
  annotated_data_multisite_filtered <- annotated_data_pre_residue_filter
}

----------------------------------------------------------------------------------------

loginfo("For sites with multisites information, get best p-value (and log fc if there are ties) for each comparison, UniProt accession, and position combination")

best_p_value <- annotated_data_multisite_filtered |>
  group_by( uniprot_acc, position, comparison) |>
  dplyr::summarise( best_p_val = min(  pval )) |>
  ungroup()

best_log_fc_pval_site_join <- c( "uniprot_acc",
                                 "position",
                                 "best_p_val",
                                 "comparison")
names( best_log_fc_pval_site_join) <- c( "uniprot_acc",
                                         "position",
                                         "pval",
                                         "comparison")

best_logfc_pval_site <- annotated_data_pre_residue_filter |>
  inner_join(best_p_value, by = best_log_fc_pval_site_join) |>
  group_by( uniprot_acc, position, comparison, pval) |>
  dplyr::summarise( best_abs_log_fc =  max(abs( fc  ))) |>
  ungroup()

best_site <- annotated_data_pre_residue_filter |>
  mutate( abs_log_fc = abs( fc)) |>
  inner_join(best_logfc_pval_site, by = c("uniprot_acc" = "uniprot_acc",
                                          "position" = "position",
                                          "pval" = "pval",
                                          "abs_log_fc" = "best_abs_log_fc",
                                          "comparison" = "comparison") ) |>
  dplyr::select(-abs_log_fc)


annotated_data <- best_site   |>
    unite( annotation,  uniprot_acc, gene_name, position, peptide_copy , sep="|"  )  |>
    distinct()

## ---------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Perform analysis of data from each contrast")

grouped_annotated_data <- annotated_data  |>
  dplyr::select( -is_multisite, -sites_id, -residue ) |>
  distinct() |>
  group_by( comparison) |>
  nest() |>
  ungroup

## ---------------------------------------------------------------------------------------------------------------------------------------------------
# grouped_annotated_data$data[[1]]


## ----------------------------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Score each phosphorylation site for similarity to kinase-specific motif.")

rds_dir <- args$output_dir

if (!args$no_backup) {
  rds_dir <- paste(args$output_dir, "_prev", sep = "")
}

if( args$reuse_old==TRUE & file.exists(file.path(rds_dir,  paste0("kinswingr_scores_list_", args$kinase_specificity, ".RDS")))) {

  scores_list <- readRDS( file.path(rds_dir, paste0("kinswingr_scores_list_", args$kinase_specificity, ".RDS")))

  captured_output<-capture.output(
    copy_attempt <- file.copy( from=  file.path(rds_dir, paste0("kinswingr_scores_list_", args$kinase_specificity, ".RDS")),
                               to =  file.path(args$output_dir, paste0("kinswingr_scores_list_", args$kinase_specificity, ".RDS")) ),
    type = "message"
  )
  logdebug(captured_output)


} else {
  # set seed for reproducible results
  set.seed(args$random_seed )

  loginfo(paste( "Number of sites per contrast: ",
                 paste( names(grouped_annotated_data$data), collapse=", "),
                 paste( purrr::map_int(grouped_annotated_data$data, \(x){ nrow(as.data.frame(x) )}), collapse=", ") ) )

  ## Use the number of data to determine the number of iterations
  num_of_iterations <- purrr::map_dbl(grouped_annotated_data$data,
                                      \(x){ num_of_rows <- nrow(as.data.frame(x) )
                                            if (num_of_rows <  args$motif_score_iteration) {
                                              motif_score_iteration_to_use <- num_of_rows - ( num_of_rows %% 10 )
                                            } else {
                                              args$motif_score_iteration
                                            }

                                      } )

  scores_list <-  grouped_annotated_data |>
    dplyr::mutate( pwms_scores = purrr::map2( data, num_of_iterations, \(.x, .y){scoreSequences(input_data = as.data.frame(.x),
                                                                   pwm_in = pwms,
                                                                   n = .y )} ))


    saveRDS( scores_list, file.path(args$output_dir,  paste0("kinswingr_scores_list_", args$kinase_specificity, ".RDS")))

}


# scores_list <- readRDS(  file.path(args$output_dir,  paste0("kinswingr_scores_list_", args$kinase_specificity, ".RDS")))

## ----------------------------------------------------------------------------------------------------------------------------------------
purrr::walk2( scores_list$data,
              scores_list$comparison,
              \(.x, .y){ vroom::vroom_write( .x  |>
                                     mutate( comparison = .y),
                                   file.path( args$output_dir,
                                              paste0( "input_data_", args$kinase_specificity, "_",  .y , ".tsv" )) ) } )

## ---------------------------------------------------------------------------------------------------------------------------------------------------

purrr::walk2( scores_list$pwms_scores,
              scores_list$comparison,
              \(.x, .y){vroom::vroom_write( .x$peptide_scores  |>
                                     mutate( comparison = .y),
                                   file.path( args$output_dir,
                                              paste0( "peptide_scores_", args$kinase_specificity, "_",  .y , ".tsv" )) ) } )

## ---------------------------------------------------------------------------------------------------------------------------------------------------

purrr::walk2( scores_list$pwms_scores,
              scores_list$comparison,
              \(.x, .y){vroom::vroom_write( .x$peptide_p  |>
                                     mutate( comparison = .y),
                                   file.path( args$output_dir,
                                              paste0( "peptide_p_", args$kinase_specificity, "_",  .y , ".tsv" )) ) } )

## ---------------------------------------------------------------------------------------------------------------------------------------------------

purrr::walk2( scores_list$pwms_scores,
              scores_list$comparison,
              \(.x, .y){ vroom::vroom_write( .x$background  |>
                                     mutate( comparison = .y),
                                   file.path( args$output_dir,
                                              paste0( "peptide_background_", args$kinase_specificity, "_",  .y , ".tsv" )) ) } )



## ----------------------------------------------------------------------------------------------------------------------------------------
loginfo("Perform randomization analysis with Swing.")

# set seed for reproducible results
set.seed(args$random_seed)


if( args$reuse_old==TRUE &  file.exists(file.path(rds_dir,  paste0("kinswingr_swing_out_list_", args$kinase_specificity, ".RDS")))) {

  swing_out_list <- readRDS( file.path(rds_dir,
                                       paste0("kinswingr_swing_out_list_", args$kinase_specificity, ".RDS")))

  captured_output<-capture.output(
    copy_attempt <- file.copy( from=  file.path(rds_dir, paste0("kinswingr_swing_out_list_", args$kinase_specificity, ".RDS")),
                               to =  file.path(args$output_dir, paste0("kinswingr_swing_out_list_", args$kinase_specificity, ".RDS")) ),
    type = "message"
  )
  logdebug(captured_output)

} else {


  swing_out_list  <-  scores_list  |>
    dplyr::mutate( swing_result = purrr::map2( data, pwms_scores, \(.x, .y) swing(input_data = as.data.frame( .x),
                                                                        pwm_in = pwms,
                                                                        pwm_scores = .y,
                                                                        permutations = args$swing_iteration,
                                                                        return_network= TRUE)))


  # This will produce two tables, one is a network for use with e.g. Cytoscape and the other is the scores. To access the scores:
  head(swing_out_list$swing_result[[1]]$scores)


  saveRDS( swing_out_list, file.path( args$output_dir,
                                      paste0("kinswingr_swing_out_list_", args$kinase_specificity, ".RDS")))
}


# swing_out_list<- readRDS(  file.path( args$output_dir,
#                                     paste0("kinswingr_swing_out_list_", args$kinase_specificity, ".RDS")))


## ---------------------------------------------------------------------------------------------------------------------------------------------------

purrr::walk2( swing_out_list$swing_result,
              swing_out_list$comparison,
              \(.x, .y){vroom::vroom_write( .x$scores |>
                                     mutate( comparison = .y),
                                   file.path( args$output_dir,
                                              paste0( "KinSwingR_", args$kinase_specificity, "_",  .y , ".tsv" )) ) })

## ---------------------------------------------------------------------------------------------------------------------------------------------------

purrr::walk2( swing_out_list$swing_result,
              swing_out_list$comparison,
              \(.x, .y){vroom::vroom_write( .x$network  |>
                                     mutate( comparison = .y),
                                   file.path( args$output_dir,
                                              paste0( "network_", args$kinase_specificity, "_",  .y , ".tsv" )) ) } )

## ----------------------------------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------------------------------
loginfo("Plot kinase-level volcano plot.")

kinase_volcano_plot <- function(
    x_label_name,
    kinase_type = "ST",
    list_position,
    brewer_pal_set = "Set1",
    p_value_cutoff = 0.05,
    ggrepel_p_value_cutoff = 0.05,
    n_cutoff = 10 ) {

  swing_scores_tbl <- swing_out_list$swing_result[[list_position]]$scores

  kinase_type_string <- case_when( kinase_type == "ST" ~ "predicted Ser/Thr PK activity",
                                   kinase_type == "Y" ~ "predicted Tyr PK activity",
                                   kinase_type == "STY" ~ "predicted Ser/Thr/Tyr PK activity")

  my_colours <- setdiff( brewer.pal(length(swing_out_list$comparison)+1, brewer_pal_set), "#FFFF33")

  significance_tbl <- swing_scores_tbl |>
    mutate( p_value = case_when( p_greater < p_less ~ p_greater,
                                 p_less <= p_greater  ~ p_less )) |>
    mutate( colour = case_when( p_value < p_value_cutoff ~ "Significant",
                                TRUE ~ "Not significant")) |>
    mutate( colour = factor( colour, levels=c(  "Not significant", "Significant") )) |>
    mutate( colour_vector = case_when( colour == "Not significant" ~ "grey",
                                       colour == "Significant" ~ my_colours[list_position] ) )

  scale_colour_vector <- significance_tbl |>
    arrange( colour) |>
    distinct( colour_vector) |>
    pull( colour_vector)

  volcano_plot_dat <- significance_tbl|>
    mutate(kinase_label = case_when( p_value < ggrepel_p_value_cutoff  &
                                       n > n_cutoff ~ kinase,
                                     TRUE ~ ""))

 if (nrow(volcano_plot_dat) > 0 ) {

    volcano_plot <- volcano_plot_dat |>
      ggplot( aes( x=swing, y = -log10(p_value + .Machine$double.eps ), size=n , col=colour, label=kinase_label, alpha=0.5))  +
      theme_bw () +
      geom_point( ) +
      scale_color_manual( values =scale_colour_vector , name="Is significant" ) +
      scale_size_continuous(name = "Num. substrates") +
      ylab( expression("Significance, -log"[10]*"(P)") ) +
      xlab( paste0( x_label_name[list_position], ", ", kinase_type_string)) +
      geom_hline(aes(yintercept = -log10(p_value_cutoff)), linetype=3)  +
      geom_text_repel( show.legend = FALSE, size = 5 ) +
      theme(text = element_text(size = 20))

    return(volcano_plot)

  } else {

    return(NA)
  }

}

partial_kinase_volcano_plot <-  partial( kinase_volcano_plot,
           x_label_name = swing_out_list$comparison,
           kinase_type = args$kinase_specificity,
           p_value_cutoff = args$p_value_cutoff,
           n_cutoff = args$min_num_sites_per_kinase,
           ggrepel_p_value_cutoff = args$ggrepel_p_value_cutoff)

kinase_volcano_plot <- purrr::map( seq_along(swing_out_list$comparison),
                                   \(x){partial_kinase_volcano_plot(list_position = x)})

for( format_ext in args$plots_format) {
  purrr::walk2( kinase_volcano_plot, swing_out_list$comparison,
                \(.x, .y) { if(length(.x) > 1 ) { ggsave( filename=file.path(args$output_dir,
                                                                    paste0(  "kinase_volcano_plot_",
                                                                             args$kinase_specificity, "_",  .y, ".", format_ext )),
                         plot=.x)}} )

}


## ----------------------------------------------------------------------------------------------------------------------------------------

plotAvgLogfcKnownSites <- function( list_position,
                                    y_label_name,
                                    brewer_pal_set= "Set1",
                                    kinase_type,
                                    p_value_cutoff = 0.05,
                                    ggrepel_p_value_cutoff = 0.05) {

  kinase_type_string <- case_when( kinase_type == "ST" ~ "predicted Ser/Thr PK activity",
                            kinase_type == "Y" ~ "predicted Tyr PK activity",
                            kinase_type == "STY" ~ "predicted Ser/Thr/Tyr PK activity")

  my_colours <- brewer.pal(length(swing_out_list$comparison), brewer_pal_set)

  my_comparison <- swing_out_list$comparison[[list_position]]

  fc_pval_tab <- scores_list$data[[list_position]]

  up_or_down_kinases <- swing_out_list$swing_result[[list_position]]$scores |>
    dplyr::filter(  p_greater < p_value_cutoff | p_less < p_value_cutoff )

  motif_score_table <- scores_list$pwms_scores[[list_position]]$peptide_scores |>
    pivot_longer( cols= !(contains("annotation") | contains("peptide")),
                  names_to="kinase",
                  values_to = "motif.score")

   peptide_p_long_tbl <- scores_list$pwms_scores[[list_position]]$peptide_p |>
    pivot_longer( cols= !(contains("annotation") | contains("peptide")),
                  names_to="kinase",
                  values_to = "motif.p.value")

  selected_scores_list <- peptide_p_long_tbl |>
     left_join( motif_score_table, by=c("annotation" = "annotation",
                                       "peptide" = "peptide",
                                       "kinase" = "kinase")) |>
    inner_join( up_or_down_kinases, by=c( "kinase" = "kinase")) |>
    inner_join( fc_pval_tab, by=c("annotation" = "annotation",
                                  "peptide" = "peptide")) |>
    left_join ( annotated_data |>
                  dplyr::select( annotation, sites_id),
                by=c("annotation" = "annotation")) |>
    left_join( de_phos |>
                 dplyr::select(sites_id, KINASE), by = c("sites_id" = "sites_id")) |>
    distinct() |>
    separate_rows( KINASE , sep= "//")  |>
    dplyr::mutate( KINASE = toupper(KINASE) ) |>
    dplyr::filter( kinase == KINASE)

  avg_logfc_vs_swing_score <- selected_scores_list |>
    group_by( kinase, swing, p_greater, p_less ) |>
    summarise( avg_fc = mean(fc),
               median_fc = median(fc),
               total_fc = sum(fc),
               counts = n()) |>
    ungroup() |>
    dplyr::filter( counts >= 5 )

  avg_logfc_vs_swing_tbl <- avg_logfc_vs_swing_score |>
    dplyr::mutate( p_value = case_when( swing > 0     ~ p_greater,
                                        swing < 0     ~ p_less,
                                        TRUE ~ 1))  |>
    mutate( colour = case_when( ( p_value < p_value_cutoff |
                                    p_value < p_value_cutoff ) ~ "Significant",
                                TRUE ~ "Not significant" )) |>
    dplyr::mutate( my_kinase = case_when( counts >= 10 |
                                            ( p_greater < ggrepel_p_value_cutoff |
                                                p_less < ggrepel_p_value_cutoff ) ~ kinase,
                                          TRUE ~ "") )

   if (nrow(avg_logfc_vs_swing_tbl) > 0 ) {

     avg_logfc_vs_swing_score_plot <- avg_logfc_vs_swing_tbl |>
       ggplot(aes( y = avg_fc, x = swing, label = my_kinase, size = counts,  color=colour, alpha=0.5 )) +
       geom_point() +
       geom_text_repel(show.legend = FALSE, size = 5) +
       scale_color_manual( values = c("grey", my_colours[list_position]), name="Is significant" ) +
       ylab( substitute( paste(a , ", Mean substrate log"[2],"(intensity) difference"), list(a=y_label_name[list_position]) ) ) +
       xlab( paste0( y_label_name[list_position], ", ", kinase_type_string)) +
       theme(text = element_text(size = 20))

     return( avg_logfc_vs_swing_score_plot)

   }

  return(NA)
}

# scores_list <- readRDS( "/home/ignatius/PostDoc/2021/Tautomycetin_2021/Results/phos_kinswingr_ST_exp1/kinswingr_scores_list_ST.RDS" )
# swing_out_list <- readRDS( "/home/ignatius/PostDoc/2021/Tautomycetin_2021/Results/phos_kinswingr_ST_exp1/kinswingr_swing_out_list_ST.RDS" )

## Do this only if the species is Human or Mouse
if( args$taxonomy_id %in% c(9606, 10090) ) {
  loginfo("Plot average log fold-change of known phoshosites.")


  partial_plotAvgLogfcKnownSites <-  partial( plotAvgLogfcKnownSites,
                                              y_label_name = swing_out_list$comparison,
                                              brewer_pal_set= "Set1",
                                              kinase_type = args$kinase_specificity,
                                              p_value_cutoff = args$p_value_cutoff,
                                              ggrepel_p_value_cutoff = args$ggrepel_p_value_cutoff)


  kinase_avg_logfc <- purrr::map( seq_along(swing_out_list$comparison),
                                  \(x) partial_plotAvgLogfcKnownSites(x) )

  for( format_ext in args$plots_format) {
    purrr::walk2( kinase_avg_logfc,
                  swing_out_list$comparison,
                  \(.x, .y){ if ( length(.x) > 1) {
                    ggsave( filename=file.path(args$output_dir,
                                               paste0( "kinase_avg_logfc_of_known_sites_",
                                                       args$kinase_specificity, "_", .y,
                                                       ".", format_ext )),
                            plot=.x) }} )
  }
}

## ---------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Compile KinSwingR results together.")

list_of_kinswinger_columns <- c(  "substrate_uniprot_acc" ,
                                 "subsrate_gene_symbol" ,
                                 "phosphosite_position" ,
                                 "sequence_context",
                                 "kinase_uniprot_acc",
                                 "kinase_gene_symbol",
                                 "kinase_family",
                                 "phosphosite_log2FC",
                                 "phosphosite_fdr_value",
                                 "motif.score",
                                 "motif.p.value",
                                 "swing_kinswingr",
                                 "p_value_kinswingr",
                                 "pos_kinswingr",
                                 "neg_kinswingr",
                                 "all_kinswingr",
                                 "pk_kinswingr",
                                 "nk_kinswingr",
                                 "swing_raw_kinswingr",
                                 "n_kinswingr",
                                 "p_greater_kinswingr",
                                 "p_less_kinswingr",
                                 "is_kinase_phosphorylated",
                                 "known_upstream_kinase",
                                 "prediction_match_known_kinase",
                                 "reactome_term",
                                 "ON_FUNCTION",
                                 "ON_PROCESS",
                                 "ON_PROT_INTERACT",
                                 "ON_OTHER_INTERACT",
                                 "REG_SITES_NOTES",
                                 "kinase_uniprot_id_human",
                                 "kinase_uniprot_acc_human",
                                 "kinase_uniprot_id_mouse",
                                 "kinase_uniprot_acc_mouse",
                                 "sites_id",
                                 "substrate_name")

compileKinswingerResults <- function( list_position ) {

  my_comparison <- swing_out_list$comparison[[list_position]]

  fc_pval_tab <- scores_list$data[[list_position]]

  up_or_down_kinases <- swing_out_list$swing_result[[list_position]]$scores |>
    dplyr::filter( p_greater < 0.2 | p_less < 0.2 )

  motif_score_table <- scores_list$pwms_scores[[list_position]]$peptide_scores |>
    pivot_longer( cols= !(contains("annotation") | contains("peptide")),
                  names_to="kinase",
                  values_to = "motif.score")

  is_phosphoylated_tbl <- annotated_data |>
    dplyr::filter( pval < args$p_value_cutoff ) |>
    dplyr::mutate( uniprot_acc = str_split(annotation, "\\|") |> purrr::map_chr(1) ) |>
    distinct( uniprot_acc) |>
    mutate(is_kinase_phosphorylated = 1)

  selected_columns <- intersect( colnames(de_phos),
                                 c( "sites_id", "PROTEIN-NAMES", "reactome_term",
                                    "KINASE",
                                    "ON_FUNCTION",
                                    "ON_PROCESS",
                                    "ON_PROT_INTERACT",
                                    "ON_OTHER_INTERACT",
                                    "REG_SITES_NOTES"))

  step_1 <- scores_list$pwms_scores[[list_position]]$peptide_p |>
    pivot_longer( cols= !(contains("annotation") | contains("peptide")),
                  names_to="kinase",
                  values_to = "motif.p.value") |>
    left_join( motif_score_table, by=c("annotation" = "annotation",
                                       "peptide" = "peptide",
                                       "kinase" = "kinase")) |>
    inner_join( up_or_down_kinases, by=c( "kinase" = "kinase")) |>
    dplyr::filter( motif.p.value < 0.2 )


  step_2 <- step_1 |>
    left_join( ks_tbl |>
                 mutate( KINASE = toupper(KINASE)) |>
                 dplyr::distinct( KINASE, KIN_ACC_ID),
               by = c("kinase" = "KINASE")) |>
    left_join( ks_tbl |>
                 mutate( GENE = toupper(GENE)) |>
                 dplyr::distinct( GENE, KIN_ACC_ID),
               by = c("kinase" = "GENE")) |>
    mutate( kinase_uniprot_acc = ifelse( is.na( KIN_ACC_ID.x),
                                         KIN_ACC_ID.y,
                                         KIN_ACC_ID.x)) |>
    dplyr::select(-KIN_ACC_ID.x, -KIN_ACC_ID.y) |>
    inner_join( fc_pval_tab |>
                  dplyr::filter( pval < args$p_value_cutoff ),
                by=c("annotation" = "annotation",
                     "peptide" = "peptide")) |>
    left_join ( annotated_data |>
                  dplyr::filter( comparison == swing_out_list$comparison[[list_position]] ) |>
                  dplyr::select( annotation, sites_id, peptide),
                by=c("annotation" = "annotation",
                     "peptide" = "peptide"))

  rm(step_1)
  gc()

  step_3 <- step_2 |>
    left_join( phosphositeplus |>
                 distinct( GENE, kinase), by=c("kinase" = "kinase") )

  rm(step_2)
  gc()

  plan(multisession, workers = args$num_cores)


  step_4 <- step_3 |>
    dplyr::mutate( substrate_gene_name  = furrr::future_map(annotation,
                                                            \(x){ str_split(x, "\\|") |> purrr::map_chr(2)}) )

  rm(step_3)
  gc()


  step_5 <-  step_4 |>
    left_join( de_phos |>
                 dplyr::filter( comparison == swing_out_list$comparison[[list_position]]  ) |>
                 dplyr::select( one_of( selected_columns ) ),
               by=c("sites_id" = "sites_id"))


  rm(step_4)
  gc()


  selected_scores_list_help <- step_5 |>
    left_join( uniprot_kinases |>
                 dplyr::select(-KEYWORDS),
               by= c("GENE" = "gene_name")) |>
    dplyr::rename( kinase_gene_name = "GENE") |>
    distinct()

  rm(step_5)
  gc()

  if( "KINASE" %in% selected_columns) {
    selected_scores_list_help <- selected_scores_list_help |>
      dplyr::mutate( kinase_copy = KINASE ) |>
      dplyr::rename( known_upstream_kinase = "kinase_copy") |>
      dplyr::rename( one_known_kinase = "KINASE") |>
      separate_rows( one_known_kinase , sep= "//")  |>
      dplyr::mutate( one_known_kinase = toupper(one_known_kinase) )  |>
      dplyr::mutate( prediction_match_known_kinase = case_when( kinase == one_known_kinase ~ TRUE,
                                                                TRUE ~ FALSE) )  |>
      dplyr::select(-one_known_kinase)
  }

  ## Use mouse or human uniprot accession if it makes sense to do so.
  selected_scores_list <- selected_scores_list_help
  if( length( intersect(   is_phosphoylated_tbl$uniprot_acc, uniprot_kinases$uniprot_acc_human  ) ) >0 )  {
    selected_scores_list <- selected_scores_list_help  |>
      left_join( is_phosphoylated_tbl, by=c("uniprot_acc_human" = "uniprot_acc"))  |>
      distinct()
  } else if ( length( intersect(   is_phosphoylated_tbl$uniprot_acc, uniprot_kinases$uniprot_acc_mouse  ) > 0 ))  {
    selected_scores_list <- selected_scores_list_help  |>
      left_join( is_phosphoylated_tbl, by=c("uniprot_acc_mouse" = "uniprot_acc"))  |>
      distinct()
  } else {

    list_of_kinswinger_columns <- setdiff(list_of_kinswinger_columns,
                                          c( "is_kinase_phosphorylated"
                                             ,"kinase_family"
                                             ,"kinase_uniprot_id_human"
                                             ,"kinase_uniprot_acc_human"
                                             ,"kinase_uniprot_id_mouse"
                                             ,"kinase_uniprot_acc_mouse"
                                             ,"known_upstream_kinase"
                                             ,"prediction_match_known_kinase"
                                             ,"reactome_term"
                                             ,"ON_FUNCTION"
                                             ,"ON_PROCESS"
                                             ,"ON_PROT_INTERACT"
                                             ,"ON_OTHER_INTERACT",
                                             "REG_SITES_NOTES"  ) )
  }

  selected_scores_list_cln_step_1 <- selected_scores_list |>
    separate( annotation, into=c( "substrate_uniprot_acc", "subsrate_gene_symbol", "phosphosite_position", "sequence_context") ) |>
    dplyr::mutate( p_value_kinswingr = case_when( swing > 0     ~ p_greater,
                                        swing < 0     ~ p_less,
                                        TRUE ~ 1)) |>
    dplyr::rename(kinase_gene_symbol = "kinase",
                  substrate_name = "PROTEIN-NAMES",
                  phosphosite_log2FC = "fc",
                  phosphosite_fdr_value = "pval",
                  swing_kinswingr = "swing",
                  pos_kinswingr = "pos",
                  neg_kinswingr = "neg",
                  all_kinswingr = "all",
                  pk_kinswingr = "pk",
                  nk_kinswingr = "nk",
                  swing_raw_kinswingr = "swing_raw",
                  n_kinswingr = "n",
                  p_greater_kinswingr = "p_greater",
                  p_less_kinswingr = "p_less")

  selected_scores_list_cln_step_2 <- selected_scores_list_cln_step_1
  if( "Family" %in% colnames( selected_scores_list_cln_step_1)) {

    selected_scores_list_cln_step_2 <- selected_scores_list_cln_step_1 |>
      dplyr::rename( kinase_family = "Family",
                     kinase_uniprot_id_human = "uniprot_id_human",
                     kinase_uniprot_acc_human = "uniprot_acc_human",
                     kinase_uniprot_id_mouse = "uniprot_id_mouse",
                     kinase_uniprot_acc_mouse = "uniprot_acc_mouse") |>
      dplyr::filter ( (args$taxonomy_id == 9606 &  kinase_uniprot_acc_human == kinase_uniprot_acc) | ## human only
                      (args$taxonomy_id == 10090 &  kinase_uniprot_acc_mouse == kinase_uniprot_acc) |  ## mouse only
                      !(args$taxonomy_id %in% c(9606, 10090)))

  }

  print( setdiff(list_of_kinswinger_columns, colnames( selected_scores_list_cln_step_2)) )

  selected_scores_list_cln_final <- selected_scores_list_cln_step_2 |>
    dplyr::select(all_of( list_of_kinswinger_columns) )

  return( selected_scores_list_cln_final)
}


gc()

selected_scores_list <- purrr::map( seq_along(swing_out_list$comparison),
                                    \(x) compileKinswingerResults(x))

names( selected_scores_list) <- swing_out_list$comparison

purrr::walk2( selected_scores_list,
              swing_out_list$comparison,
              \(.x, .y) { vroom::vroom_write( .x  |>
                                     mutate( comparison = .y) |>
                                     relocate( comparison, .before="substrate_uniprot_acc"),
                                   file.path( args$output_dir,
                                              paste0( "selected_kinase_substrate_",
                                                      args$kinase_specificity, "_",
                                                      .y ,
                                                      ".tsv" )) )} )

## ----------------------------------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#  RColorBrewer::display.brewer.all()

## ---------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))

# y_label_name <- c("hello")
# list_position <- 1
# iris |>
# ggplot( aes(x=Petal.Width, y=Sepal.Length)) +
#   geom_point() +
#   labs (title = substitute( paste(a , ", Mean substrate log"[2],"(intensity) difference"), list(a=y_label_name[list_position]) ) )
