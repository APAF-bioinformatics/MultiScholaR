# Package-level imports that are shared across multiple workflow modules.

#' @importFrom graphics text title
#' @importFrom grDevices as.raster cairo_pdf colorRampPalette
#' @importFrom methods callNextMethod canCoerce
#' @importFrom igraph E<- V<-
#' @importFrom limma duplicateCorrelation plotSA
#' @importFrom purrr cross2 discard keep map2_chr map2_lgl
#' @importFrom RUVIIIC RUVIII_C_Varying
#' @importFrom stats contrasts density median pt reorder rnorm runif step
#' @importFrom stringr str_locate_all str_match_all str_sub str_wrap
#' @importFrom utils capture.output combn data download.file flush.console getFromNamespace install.packages modifyList object.size read.delim sessionInfo tail unzip
NULL

# Symbols evaluated inside dplyr, tidyr, ggplot2, and Shiny data masks.
utils::globalVariables(c(
  ".", ".impute_design_order", ".impute_excluded", ".impute_group", ".impute_group_size", ".impute_peptide",
  ".impute_protein", ".impute_quantity", ".impute_run", ".impute_source_order", ".new_itsd_id", ".original_row_id",
  ".peptide_feature_key", "abs_factor1", "Abundance", "acc_order_id", "accession", "adj.P.Val",
  "adjusted_min_reps", "All", "analysis_type", "Analysis_Type", "annot_id_list", "annotation",
  "annotation_id", "annotation_score", "assay", "AveExpr", "average_value", "Averagine",
  "avg_detected", "avg_intensity", "avg_peptides", "avg_proteins", "Avg.Log2.Protein.Imputed", "batch",
  "Batch", "best_k", "best_p_adj_value", "best_phos_pos", "best_phos_prob", "best_protein",
  "best_uniprot_acc", "category", "cc", "ChEBI", "chebi_id", "clean_acc",
  "cleaned_acc", "cleaned_id", "cleaned_peptide", "coefficient_of_variation", "Cohort_Sample", "col_name",
  "collaborator_patient_id", "colour", "compare_column", "comparision_short", "comparison", "composite_score",
  "contrast", "contrast_name", "Control", "core_enrichment", "Count", "counts",
  "CScore", "custom_comparison", "cv", "cv_percent", "cv_status", "database_id",
  "database_identifier", "db", "description", "Description", "detail", "direction",
  "Direction", "directionality", "display_name", "eligible_for_imputation", "Eliminated", "enrichmentScore",
  "Entry", "Evidence", "evidence_col_to_use", "evidence_id", "experiment", "experiment_label",
  "experiment_paths", "Factor1", "falseDiscoveryRate", "FDR", "fdr_qvalue", "fdr_safe",
  "fdr_value_bh_adjustment", "feature", "feature_id", "feature_key", "featureset", "file_name",
  "File.Name", "First.Protein.Description", "formatted", "Fragment.Correlations", "Fragment.Quant.Raw", "friendly_name",
  "From", "frozen_protein_table", "gene", "gene_list_position", "gene_name", "gene_name_significant",
  "gene_names", "gene_names_first", "gene_set", "gene_symbol", "Gene.Names", "Gene.Ontology.IDs",
  "geneID", "general_sample_info", "Genes", "Genes.MaxLFQ", "Genes.MaxLFQ.Unique", "Genes.Normalised",
  "Genes.Quantity", "genesMapped", "GeneSymbol", "genotype_group", "GG.Q.Value", "Global.PG.Q.Value",
  "Global.Q.Value", "GN", "go_id", "go_term", "go_type", "GO-IDs",
  "group", "group_has_technical_replicates", "has_value", "header", "high_quality_sites", "hsa_pathway_id",
  "id", "ID", "iIM", "IM", "index", "intensity",
  "Intensity", "intersection_size", "iRT", "is_below", "is_candidate", "is_HEK",
  "is_id", "IS_ID", "is_isoform", "is_replicate_temp", "is_selected", "is_unique",
  "isoform_num", "join_uniprot_acc", "K", "KEGG", "kegg_id", "label",
  "leading_proteins", "left", "left_group", "Levels", "Lib.PG.Q.Value", "Lib.Q.Value",
  "lipid_id", "lipid_name", "log_intensity", "log_values", "log2_intensity", "log2_value",
  "Log2.Protein.Imputed", "log2FC", "Log2Intensity", "log2norm", "logFC", "logFC_safe",
  "lqm", "map_id", "mappedIDs", "Mass.Evidence", "matching_id", "max_missing",
  "max_set_size", "maxquant_row_id", "maxquant_row_ids", "mean_intensity", "mean_na_percentage", "mean_observed_raw_quantity",
  "metabolite", "metabolite_id", "metabolite_identification", "metabolite_name", "min_set_size", "Missing_Count",
  "missing_distinct_runs", "missing_fraction", "missing_percent", "missingness", "Modified.Sequence", "ms_filename",
  "ms_filename.x", "ms_filename.y", "Ms1.Area", "Ms1.Normalised", "Ms1.Profile.Corr", "MS2.Scan",
  "my_value", "n_detected", "n_lipids", "n_metabolites", "n_peptides", "n_proteins",
  "n_reps", "na_percentage", "name", "Name", "ncbi_refseq", "neg_log_10_fdr",
  "neg_log_p_value", "neg_log10_p", "neg_log10_q", "neg_log10_qvalue", "negLog10FDR", "NES",
  "new_chebi_id", "No. of Missing Values", "No. of Values", "norm_name", "Normalisation.Factor", "num_below_per_group",
  "num_below_threshold", "num_gene_names", "num_groups", "num_groups_failed", "num_missing", "num_missing_per_group",
  "num_missing_values", "num_missing.left", "num_missing.right", "num_observed_above", "num_peptides_after_impute", "num_per_group",
  "num_present", "num_present.left", "num_present.right", "num_proteins_with_values", "num_tech_rep", "num_values",
  "numeric_id", "observed_distinct_runs", "omic_type", "ordering", "organism", "Organism",
  "organism_name", "output_name", "p_value", "P.Value", "parameter_name", "pathway",
  "pathway_id", "pathway_name", "pathway_num", "PC1", "PC2", "pearson_correlation",
  "PEP", "Peptide", "peptide_count", "peptide_location", "Peptide.Imputed", "Peptide.Normalised",
  "Peptide.RawQuantity", "peptides_for_protein_count", "peptides_status", "peptidoform_count", "peptidoforms_for_protein_count", "perc_below_per_group",
  "perc_missing", "perc_missing_per_group", "perc_present", "percent_groups_failed", "percent_missing", "percentage",
  "percentage_bin", "percentage_requested", "PG.MaxLFQ", "PG.Normalised", "PG.Q.Value", "PG.Quantity",
  "phos_15mer_seq", "phospho_sty", "phospho_sty_probabilities", "potential_contaminant", "Precursor.Charge", "Precursor.Id",
  "Precursor.Normalised", "Precursor.Quantity", "Predicted.iIM", "Predicted.IM", "Predicted.iRT", "Predicted.RT",
  "proportion_below_threshold", "Protein", "protein_count", "protein_evidence", "Protein_existence", "protein_id",
  "protein_id_column", "protein_ids", "Protein_Ids", "protein_site_positions", "Protein.Group", "Protein.Ids",
  "Protein.Imputed", "Protein.Names", "Protein.Normalised", "Protein.Q.Value", "proteins", "Proteotypic",
  "pvalue", "q.mod", "Q.Value", "Quantiles", "Quantity.Quality", "rank_suffix",
  "ranking", "raw_pvalue", "reactome_id", "reactome_term", "replicate_number", "replicates",
  "results_dir", "reverse", "revigo_results", "right", "right_group", "rle.x.factor",
  "row_id", "row_id_column_with_isoform", "RT", "RT.Start", "RT.Stop", "Run",
  "Sample", "sample_id", "Sample_ID", "sample_type", "samples", "Samples ID",
  "score", "sd_intensity", "search_field", "separation_score", "seq_length", "set_size",
  "set_type", "setSize", "significant", "significant_label", "sites_id", "sites_id_short",
  "source_dir", "species_id", "Spectrum.Similarity", "Stage", "status", "Stripped.Sequence",
  "sum_intensity", "tables", "taxid", "taxon_id", "tech_rep_numbers", "tech_reps",
  "techRepNumbers", "temp", "temp_base", "temp_check_pos", "temp_column", "temp_column_base",
  "temp_protein_site_positions", "temp_qvalue", "temp_sample_id", "temporary_values_choose_accession", "term", "term_id",
  "term_name", "term_size", "termDescription", "termDescriptionAbbrev", "termID", "total_distinct_runs",
  "total_lipids", "total_metabolites", "total_num_tech_rep", "total_peptides", "Translated.Q.Value", "type",
  "uniparc_id", "uniprot_acc", "uniprot_acc_cleaned", "uniprot_acc_copy", "uniprot_acc_first", "uniprot_acc_split",
  "UNIPROT_GENENAME", "uniprot_id", "uniprot_list", "uniprot_tbl", "UNIPROTKB", "updated_tables",
  "use_colour", "use_shape", "V1", "V2", "value", "values",
  "view"
))
