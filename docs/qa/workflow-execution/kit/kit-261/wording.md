# KIT-261: texts for the owner's wording pass (VAL-41)

All drafts. Nothing users are told about the science changes.

| Draft | Where | Code |
|---|---|---|
| "Please wait: a step is still running..." | the Run button's label on the lipidomics, metabolomics and transcriptomics Normalization tab, managed application, while any step of the omic runs | `omicNormRunButtonTexts()`, `R/mod_omic_norm.R` |
| NEWS entry "Lipidomics, metabolomics and transcriptomics normalisation, managed application: 'Run Normalization Pipeline' is greyed out ..." | `NEWS.md`, `### Fixes` | |
| (unchanged, KIT-194's) "A step is still running. Press "Run Normalization Pipeline" again when it has finished." | the refusal, kept as the server-side guard; already in the owner's list | `omicNormGuardTexts()$busy` |
