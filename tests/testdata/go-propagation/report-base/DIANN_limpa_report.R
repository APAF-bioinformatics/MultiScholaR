# KIT-173 T11(e): the GO Enrichment Section chunk of inst/reports/proteomics/DIANN_limpa_report.rmd at a60c4e9a (base), verbatim.
# UPDATED: GO enrichment using new organism and UI parameters
FormatGOEnrichment <- function(numbers) {
  # Extract organism info (if available in future iterations)
  organism_name_val <- "the species under study" # Default since not in current format
  taxon_id_val <- numbers$ui_enrichment_organism_selected %||% NA
  
  # Get UI enrichment parameters
  up_lfc_val <- numbers$ui_enrichment_up_lfc %||% 0
  down_lfc_val <- numbers$ui_enrichment_down_lfc %||% 0
  q_val_cutoff_val <- numbers$ui_enrichment_q_value_cutoff %||% 0.05
  database_source_val <- numbers$ui_enrichment_database_source %||% "clusterprofiler"

  taxon_info <- if (!is.na(taxon_id_val)) glue::glue(" (NCBI Taxonomy ID: {taxon_id_val})") else ""

  background_info <- glue::glue(
    "• The background protein set consisted of all identified proteins from ",
    "{organism_name_val}{taxon_info} that passed the quality control filtering steps."
  )

  # Text for LFC cutoffs used to define gene sets
  lfc_text <- if (up_lfc_val > 0 || down_lfc_val < 0) {
      glue::glue(
          "• Upregulated proteins were defined as those with a log2 fold-change > {up_lfc_val}. ",
          "Downregulated proteins were defined as those with a log2 fold-change < {down_lfc_val}."
      )
  } else {
      "• No fold change cutoff was applied to filter differentially expressed proteins."
  }

  enrichment_tool_text <- if (database_source_val == "gprofiler") {
    # g:Profiler text
    glue::glue("The enrichment analysis was performed using the g:Profiler toolset ",
    "via the gprofiler2 R package (Kolberg et al., 2023). For each GO category ",
    "(Biological Process, Molecular Function, and Cellular Component), terms were ",
    "considered significantly enriched if they had an adjusted p-value less than ",
    "{q_val_cutoff_val} after multiple testing correction using the g:SCS method.")
  } else {
    # clusterProfiler text (default)
    glue::glue("The enrichment analysis was performed using the clusterProfiler R package ",
    "(Yu et al., 2012). For each GO category (Biological Process, Molecular Function, ",
    "and Cellular Component), terms were considered significantly enriched if they had ",
    "an adjusted p-value less than {q_val_cutoff_val} after multiple testing correction ",
    "using the Benjamini-Hochberg method.")
  }
  
  glue::glue("
GO ENRICHMENT ANALYSIS

Gene Ontology (GO) enrichment analysis was performed on the differentially 
expressed proteins to identify significantly enriched biological processes, 
molecular functions, and cellular components. The analysis was conducted using the 
following parameters:
{lfc_text}
• An adjusted FDR cutoff of {q_val_cutoff_val} was used to identify significantly enriched GO terms
{background_info}

{enrichment_tool_text}

Key findings from the GO enrichment analysis are summarised below for each 
comparison:

Biological Process:
• Upregulated processes included cellular stress response, metabolic regulation, 
  and protein folding pathways
• Downregulated processes were primarily related to cell adhesion, cytoskeletal 
  organisation, and vesicular transport

Molecular Function:
• Proteins with binding functions (particularly nucleotide and protein binding) 
  were significantly enriched
• Catalytic activities, especially those involved in redox reactions, showed 
  significant enrichment patterns

Cellular Component:
• Significant enrichment was observed for proteins localized to membrane-bound 
  organelles, particularly mitochondria and endoplasmic reticulum
• Extracellular components showed differential regulation across conditions

The complete list of enriched GO terms with statistics is provided in the 
supplementary file 'Pathway_enrichment_results.xlsx' within the Publication_tables 
folder.",
.trim = FALSE)
}
