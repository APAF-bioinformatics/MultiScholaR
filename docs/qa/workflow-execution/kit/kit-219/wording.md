# KIT-219 user-facing texts: DRAFTS for the owner's wording pass (VAL-41; OD-4)

Where used: `R/utils_workflow_summary.R` (`workflowEnrichmentTexts`, the Copy notice and the Report refusal) and the two
report templates. The ticket's drafts (7) with the deviations below.

- T1 (Copy, warning, proteomics): "Publication files copied. <n> enrichment file(s) in the pathway_enrichment folder were
  not copied, because the current enrichment did not make them (they come from an earlier enrichment or an earlier
  analysis). They are still on disk: <first three names> [and <m> more]." Deviation from the ticket's draft: the trailing
  ellipsis is "and <m> more" (the form of KIT-212's refusal text).
- T1b (Copy, warning, lipidomics, metabolomics, transcriptomics): "Publication files copied. <n> enrichment file(s) in the
  pathway_enrichment folder were not copied: this application has no enrichment analysis for <omic>, so they were not made
  by this analysis. They are still on disk: <first three names> [and <m> more]."
- With optional artifacts missing the same sentence follows "Publication files copied with optional missing artifacts."
  in place of the first words; after a failed combined workbook it follows that notice.
- T2 (Report refused; the notice reads "Report generation failed: <T2>"): "The report was not generated: the publication
  folder holds enrichment results that the current enrichment did not make. Press 'Copy to Publication Directory' again,
  then generate the report." The notice then says "Report generation failed" and "The report was not generated": one of
  the two should go at the wording pass.
- Lipidomics and metabolomics templates: "No enrichment analysis was run for this project." under the heading "Pathway
  Enrichment Analysis" (the proteomics sentence). T3 (the proteomics sentence for an enrichment that is no longer current) is
  not built: it needs a value passed to the report (F7); the one sentence above is printed.
- NEWS (`### Fixes`): the ticket's draft (7) with the OD-5 (i) sentence and the template-copy limit sentence.

## Added after the implementation audit (K219-I08, I12; for the owner, not built)

- T3 (when the enrichment stops being current: the DA re-run, a Revert or Reset, no product, another machine) the proteomics report
  prints "No enrichment analysis was run for this project.", which is untrue: one was run, on an earlier analysis. Before this change
  the same report showed the earlier enrichment's figures as if current. The templates are NOT changed here. The register's three
  options: (a) keep the sentence until the wording pass and state the case in NEWS meanwhile; (b) one sentence for both cases now,
  "No enrichment result of the current differential analysis is included." (needs no new mechanism; edits four templates and one
  hand-written assertion, `tests/testthat/test-prot-summary-report-templates.R:43-56`); (c) T3 as drafted, "The enrichment analysis of
  this project was run on an earlier differential analysis and is not included; run it again to include it." (needs a value passed to
  the report, F7). Register's recommendation: (b), or (c) at the wording pass, before the release commit.
- T2 has no variant for omics without an enrichment analysis (T1b's form): for lipidomics, metabolomics and transcriptomics it says
  "the current enrichment"; a T2b would say "this application has no enrichment analysis for <omic>".
- T2's doubled opening ("Report generation failed: The report was not generated: ...").
