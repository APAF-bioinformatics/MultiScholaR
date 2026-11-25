# KIT-214: user-facing text, drafts for the owner (VAL-41)

The owner sets the final wording. Nothing here is a decision.

| Where | Text (draft) | Status |
|---|---|---|
| Refusal, `omicsDuplicateResolutionSampleColumns` in `R/utils_workflow_omics_declarations.R` (shipped as written) | "Duplicate features cannot be resolved in assay '<assay>': the design's sample column(s) <up to three> are missing from it or do not hold numbers. Check the design matrix and the imported table." | open: the owner's wording; with no sample left it reads "(none left)" (AUDIT-KIT-214-impl-20261125 I11) |
| NEWS entry, sentence added after "...could decide which row was kept." (in `NEWS.md`) | "Where the duplicated rows were measured in different numbers of samples, the row whose measured values average highest is now kept, even if it was measured in fewer samples." | open: draft (I01); depends on the owner's N1 answer |
| NEWS entry, the refusal sentence (not changed in `NEWS.md`) | now: "A design sample that is missing from a table with a duplicated identifier, or is not numeric, is now refused with a message." Proposed: "When an identifier is listed more than once, a design sample that is missing from that table or does not hold numbers is now refused with a message." | open (I03) |
| NEWS entry, "Not affected" list (changed in `NEWS.md`) | "the metabolomics workbooks" replaced by "scripts that call the S4 method `resolveDuplicateFeatures()` (it already averaged the design's samples)" (VAL-209: no workbook ships) | applied as drafted (I03); final wording the owner's |
