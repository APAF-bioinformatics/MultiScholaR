# KIT-271 wording drafts (VAL-41: the owner sets the final wording)

The NEWS bullet below is the text as the branch had it. The merge commit of merge train 26 carries a version with four precision points of
AUDIT-KIT-271-impl-20261126 I06, marked here as DRAFTS for the wording pass: (1) the Skip mode computes no RUV factors; (2) in the correction mode the
delivered `normalised_results.tsv` is the uncorrected table; (3) a normalisation press also removes the old file; (4) with the limpa default
(method none) "normalised" means "before RUV", not that a between-sample normalisation ran. The owner sets the final wording.

## NEWS (Fixes, under "## MultiScholaR (unreleased; becomes 0.9.1, the first alpha for testing on real data)")

Proteomics, DIA with limpa: the table written by the correlation filter, and the results package's copy of it, no longer
carry the name of an RUV-normalised table. The table is the protein table before any RUV step (RUV factors are computed
after the filter), but it was stored as `ruv_normalised_results_cln_with_replicates.tsv` and delivered as
`Publication_tables/RUV_normalised_results.tsv`; it is now stored as `normalised_results_cln_with_replicates.tsv` and
delivered as `normalised_results.tsv`, as on the other proteomics routes when RUV did not correct the table. The values are
unchanged. In the RUV correction mode the matrix the differential analysis used (RUV-III's corrected matrix) is not
delivered in the results package, and in the factors-as-covariates mode neither are the factors' design columns; both stay
as before. Projects saved by an earlier build need no re-run: Copy delivers the table they hold as `normalised_results.tsv`
and the session export works; the old file name stays in their protein QC folder until the correlation filter is run
again. This route is new in 0.9 (tag v0.9.0 and later). (KIT-271)

## Console (the correlation filter on DIA with limpa)

Before: `*** CORRELATION: RUV was applied - using 'ruv_normalised' file names ***`
After (draft): `*** CORRELATION: DIA with limpa - the correlation filter's table precedes any RUV step; using 'normalised' file names ***`

## Console (Copy, DIA with limpa)

Draft: `COPY: Using normalized file names (DIA with limpa: the table precedes any RUV step)`

Labels on Copy's list: "Normalized Results TSV" and "Normalized Results RDS" (the labels the other routes already use).
