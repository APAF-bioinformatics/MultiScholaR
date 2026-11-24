# KIT-167: every user-facing sentence of the reference mapping tool in one place (the owner's wording list, section 12).
# The refusal texts are `sprintf` templates; the class of each refusal is its name.

#' The reference mapping tool's user-facing texts (KIT-167, wording list)
#' @noRd
referenceMappingTexts <- list(
    reference_mapping_engine_missing = paste(
        "The sequence search program DIAMOND is not available in this installation, so a new mapping cannot be made.",
        "Saved mappings and the demo tables still work."),
    reference_mapping_engine_mismatch = paste(
        "The sequence search program at %s is not the verified DIAMOND %s for this computer (its checksum or version",
        "differs), so a new mapping cannot be made. Saved mappings and the demo tables still work."),
    reference_mapping_engine_runtime_missing = paste(
        "The sequence search program DIAMOND could not start because the Microsoft Visual C++ runtime (Visual C++ v14,",
        "64-bit) is missing or too old on this computer. Install it from Microsoft's download page",
        "(https://learn.microsoft.com/cpp/windows/latest-supported-vc-redist) and try again. Saved mappings and the demo",
        "tables still work."),
    reference_mapping_engine_failed = "The sequence search program DIAMOND stopped with an error (exit status %s): %s",
    reference_mapping_timeout = paste(
        "The mapping did not finish within %d minutes and was stopped; nothing was kept. A strain and a bacterial",
        "reference normally take a few minutes."),
    reference_mapping_taxon_ambiguous = paste(
        "The protein file holds UniProt entries of several organisms (%s) and no strain taxon was given, so the strain's",
        "own proteins cannot be told apart. Give the strain's taxon."),
    reference_mapping_duplicate_id = "The protein file gives different sequences under the same identifier: %s.",
    reference_mapping_not_protein = paste(
        "The selected entries do not look like protein sequences (%s). The tool needs a protein FASTA file, not",
        "nucleotides."),
    reference_mapping_too_few = paste(
        "Only %d distinct strain proteins were found in the file (at least 50 are needed). Check that the file is the",
        "strain's protein FASTA."),
    reference_mapping_reference_invalid = "The reference cannot be used: %s",
    reference_mapping_download_failed = "The current UniProt release of %s could not be downloaded: %s",
    reference_mapping_download_unavailable = paste(
        "Downloading the current UniProt release needs the httr2 package, which is not installed. The shipped",
        "references work without it."),
    reference_annotation_missing = "Choose the strain's genome annotation file (GFF3 or GTF) to match its genes.",
    reference_annotation_malformed = "The annotation file is not GFF3 or GTF: line %d does not have 9 tab-separated fields.",
    reference_annotation_no_protein_id = paste(
        "The annotation file names no protein for any coding sequence (no protein_id attribute), so genes cannot be",
        "matched to proteins. Use the GFF3 file of the same genome, which carries protein_id."),
    reference_annotation_not_mapped_share = paste(
        "%d of %d proteins named in the annotation (%.1f%%) are not in the mapped protein file; the annotation and the",
        "protein file may come from different versions of the genome."),
    applicability = paste(
        "Validated for a strain mapped to a reference of the same species. For any other reference the mapping is",
        "reported but is not covered by the validation."),
    identity_caution = paste(
        "The median identity of the one-to-one pairs is %.1f%%, below 95%%. Strains of one species are usually above",
        "it; check that the reference is of the strain's species."),
    current_release_label = "Current UniProt release (not the validated reference)",
    organism_filter_not_needed = paste(
        "With orthologue identifiers only the mapped strain's proteins are used; the organism filter is not needed.")
)

#' A refusal of the reference mapping tool: a classed error with the plain sentence of `referenceMappingTexts`.
#' @param class The refusal class (a name of `referenceMappingTexts`).
#' @param ... Values for the sentence's template.
#' @param fields Named list of fields carried on the condition.
#' @return Does not return.
#' @noRd
referenceMappingRefuse <- function(class, ..., fields = list()) {
    text <- referenceMappingTexts[[class]]
    if (is.null(text)) cli::cli_abort("unknown reference mapping refusal {class}", class = "reference_mapping_internal")
    message <- if (...length()) sprintf(text, ...) else text
    do.call(rlang::abort, c(list(message = message, class = c(class, "reference_mapping_refusal")), fields))
}
