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
