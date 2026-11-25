# KIT-235 wording for the owner's pass (VAL-41)

Drafts only; the owner sets the final wording.

## N1. NEWS entry (NEWS.md, unreleased heading, "Fixes")

The entry is a DRAFT. Its fourth and fifth sentences are the implementation audit's text (AUDIT-KIT-235-impl-20261125,
K235-I01 (b)), which replaced "The Methods now state the number of assays, the internal-standard normalisation, … as
recorded." so that NEWS does not read as if the report were now correct. On both observed GUI lanes the
internal-standard sentence the fix makes visible is false (no internal standard matched; KIT-242), and the report prints
the fixed and unsupported sentences KIT-236 corrects.

Draft as committed:

> Reports (lipidomics, metabolomics): the report now reads the study-parameters file the application writes. Since
> March 2026 the application wrote this file in a form these two reports did not read, so a report generated in the
> application showed 'Error running Main function: …' in place of its Sample Details and no Methods, Results, Comments
> or differential-abundance text. The Methods now print again, with the number of assays, the log transformation, the
> between-sample normalisation and, for each assay, the number of unwanted components removed and of control features,
> as the study-parameters file records them. Some of the text these two reports now print again is not computed from
> the data, or can be false for the analysis that was run (for example, internal-standard normalisation is reported as
> applied when no internal standard matched); that text is corrected in this release (KIT-236). A report generated
> after a save that failed (the short study-parameters file the application writes on an error) still shows 'Error
> running Main function: …'; that file records no assays. A project keeps the copy of the report template made at its
> first report; projects that generated a lipidomics or metabolomics report before this version keep the earlier
> template unless that copy (`scripts/<omic>/<omic>_report.rmd`) is removed (KIT-235).

## Report text this ticket makes users read again

Unchanged template text, listed in ticket section 11 (V1 to V7). What it becomes is KIT-236's owner question, plus
KIT-219, KIT-227's wording pass and the wording pass's SUM-7.
