# KIT-273 user-facing text: drafts for the owner's wording pass (VAL-41)

Substance and the two file names are decided (VAL-332); every sentence below is a draft.

| Where | Text |
|---|---|
| File, correction mode | `Publication_tables/ruv_corrected_results.tsv` |
| File, covariates mode | `Study_report/design_matrix_with_ruv_factors.tab` |
| Display name, correction | RUV-III Corrected Results TSV |
| Display name, covariates | Design Matrix with RUV-III Factors |
| Copy console line, correction | COPY: DIA with limpa, RUV-III correction: the corrected matrix the analysis used is delivered as Publication_tables/ruv_corrected_results.tsv |
| Copy console line, covariates | COPY: DIA with limpa, RUV-III factors as covariates: the design with the factors the analysis used is delivered as Study_report/design_matrix_with_ruv_factors.tab |
| Failed-copy record (Copy's summary) | display name above, the error text, destination folder; type `ruv_delivery` |
| NEWS | the bullet beginning "Proteomics, DIA with limpa: the results package now holds what the differential analysis used" in NEWS.md (Fixes); it adds one sentence the ticket's draft did not carry: after a RUV run and "Skip RUV" the next analysis still uses the factors or the corrected matrix (KIT-276 has not landed), so the file is then still delivered |
| Failed-copy record, generic display name (AUDIT-KIT-273-impl I09) | RUV-III delivery table (used when the failure comes before the mode is decided; destination is then the folder only) |
| Draft for the owner's pass (merge train 34): the two file names | `Publication_tables/ruv_corrected_results.tsv` (correction) and `Study_report/design_matrix_with_ruv_factors.tab` (covariates): decided in VAL-332; the wording pass confirms them |
| Draft for the owner's pass (merge train 34): the two console lines | the two `COPY:` lines of the rows above, and the failed-copy summary entry (display name, error, destination folder) |
| Draft for the owner's pass (merge train 34): NEWS "Skip RUV" sentence | The sentence "pressing 'Skip RUV' after a RUV run does not yet undo that run ... so the file is then still delivered" is true today (VAL-334 (2)) and becomes untrue when KIT-276 lands (Skip then removes the RUV result, VAL-346 Q1); it is then removed from NEWS (AUDIT-KIT-273-impl I07; KIT-276's A7 is to carry it) |
