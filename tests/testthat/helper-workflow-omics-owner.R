# KIT-070/071: the owner fixtures shared by test-workflow-omics-owner.R and test-workflow-omics-sessions.R: a fake
# Shiny session, a four-sample design and one small object per omic class.

omicsOwnerSession <- function() {
  session <- new.env(parent = emptyenv())
  session$userData <- new.env(parent = emptyenv())
  session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
  session
}

omicsOwnerDesign <- function() {
  data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
}

omicsOwnerLipidObject <- function(n = 6L) {
  assay <- function(prefix) {
    data.frame(lipid_id = paste0(prefix, "_L", seq_len(n)), lipid_name = paste("Lipid", seq_len(n)),
      S1 = seq_len(n) + 0.5, S2 = seq_len(n) + 1.5, S3 = seq_len(n) + 2.5, S4 = seq_len(n) + 3.5,
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  methods::new("LipidomicsAssayData", lipid_data = list(pos = assay("pos"), neg = assay("neg")),
    design_matrix = omicsOwnerDesign(), sample_id = "Run", group_id = "group",
    lipid_id_column = "lipid_id", annotation_id_column = "lipid_name")
}

omicsOwnerMetabObject <- function(n = 5L) {
  assay <- function(prefix) {
    data.frame(metabolite_id = paste0(prefix, "_M", seq_len(n)), metabolite_name = paste("Met", seq_len(n)),
      S1 = seq_len(n) * 2, S2 = seq_len(n) * 3, S3 = seq_len(n) * 4, S4 = seq_len(n) * 5,
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  methods::new("MetaboliteAssayData", metabolite_data = list(LCMS_Pos = assay("pos"), GCMS = assay("gc")),
    design_matrix = omicsOwnerDesign(), sample_id = "Run", group_id = "group",
    metabolite_id_column = "metabolite_id", annotation_id_column = "metabolite_name")
}

