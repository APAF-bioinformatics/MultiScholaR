# The limpa script mapped onto the kit: the `DIA_LIMPA` workflow (2026-10-19)

Source: `APAF-bioinformatics/e_35036_MorpheusTex_GarryNunn_20260329_b` (private), `main`, `scripts/proteomics_isoseq/DIA_workflow_limpa_isoseq_tissue.Rmd` (2,175 lines); the twins `DIA_workflow_limpa_isoseq_ev.Rmd`, `DIA_workflow_limpa_isoseq_juice.Rmd` and `scripts/proteomics_metagenomics/DIA_workflow_limpa_metagenomics_tissue.Rmd` differ only in project constants (experiment label, sample regex, contrasts, FASTA, taxon, the import reader and an annotation block). The script ran against MultiScholaR branch `workshop` at `d9af2dde`, where `deAnalysisWrapperFunction` is a deprecated wrapper over `protein_deAnalysisWrapperFunction`; on `performance-v2` the equivalent limpa detection lives in `daAnalysisWrapperFunction` (`R/func_prot_da_model_wrapper.R`), which the GUI does not call. The one customised function is `findBestNegCtrlPercentageCustom` (script lines 1086–1172); no other `Custom` function exists in the four scripts.

Tickets: KIT-049 (the DIA route's optional imputation and the entry-point labels), KIT-046 (the type), KIT-047 (peptide stage), KIT-041 (protein stage), KIT-053 (DA, record, report, journey). Kernels are never rewritten; the two additions named below are the script's own code (the capped optimiser) and the existing dpcDE branch reached on the new route only.

## Chunk by chunk

| # | Script chunk | `DIA_LIMPA` step (state, step label) | Kernel | Unit | Ticket |
|---|---|---|---|---|---|
| 4 | Import DIA-NN report (parquet or tsv; `^Cont_` contaminants dropped) | `initial`, owned table | import module as DIA | store | KIT-046 |
| 5–8 | FASTA processing, UniProt lookup, design matrix | `raw_data_s4` (`design_builder`) | as DIA | store | KIT-046 |
| 9–10 | `PeptideQuantitativeData`; raw QC (`1_Raw Data`) | as DIA | — | — | — |
| 11 | q-value and proteotypic filter (`2_qval Filtered`) | `peptide_qvalue_filter` (shared) | `srlQvalueProteotypicPeptideClean` | feature group | shared |
| 12 | precursor → peptide rollup (`3_peptidoform count`) | `precursor_rollup` (shared, P14) | `rollUpPrecursorToPeptide` | complete sample | shared |
| 13 | peptide matrix and log2; RLE/PCA/density/Pearson QC (`4_QC_raw`) | `peptide_log2` (`4_log2_matrix`) | `calcPeptideMatrix`, `log2TransformPeptideMatrix` | row | KIT-047 |
| 16 | best accession per peptide (`7_annotation_cleanup`) | `peptide_accession` (`7_annotation_cleanup`) | `chooseBestProteinAccession` (peptide method; `delim = ";"`, `aggregation_method = "mean"`, FASTA table) | feature group | KIT-047 |
| 17 | Limpa DPC imputation and protein quantification; limpa QC plots; protein RLE/PCA/density/Pearson | `dpc_quant` (`9_protein_dpc_quant`) | `proteinMissingValueImputationLimpa` (peptide method: DPC fit, `dpcQuant`; `engine = fastlimpadpc` default, `limpa` selectable); `generateLimpaQCPlots` as a product | whole matrix | KIT-047 |
| 17b | protein cyclic loess (the script skips it: `protein_normalised_obj <- protein_imputed_obj`) | `protein_normalisation` (`10_protein_normalised`), default `none` | `normaliseBetweenSamples` | whole matrix | KIT-041 |
| 18 | ≥ 1 unique peptide and ≥ 2 peptidoforms per protein (`8_peptideperprotein_filter`) | `peptide_count` (`11_peptide_count_filtered`) | `filterMinNumPeptidesPerProtein` (protein method) | row | KIT-041 |
| 19a | Pearson sample filter, threshold 0.5 (`9_samplecorr_filter`) | `correlation` (`12_samplecorr_filtered`; shared `pearson-pairs`, `pearson-selection`) | `pearsonCorForSamplePairs`, `filterSamplesByProteinCorrelationThreshold` | whole matrix | shared |
| 19b | protein-level RUV-III as covariates | `ruv_covariates` (`13_ruv_covariates`) | `findBestNegCtrlPercentageCapped` (the script's optimiser), `ruv::RUVIII`, W computation; the adjusted matrix as a product for QC plots | whole matrix | KIT-041 |
| 20–21 | audit-trail files; composite QC figure | `session_export`; products | as DIA | store | shared |
| 22 | DE with `formula + W1 + … + Wk`; dpcDE when limpa data are present | `da` (`DA-contrast-output`, `model = dpcDE`) | `limpa::dpcDE`, `convertDpcDAToStandardFormat` | whole matrix | KIT-053 |
| 23, 23a, 22b | DE outputs, Glimma volcano, summary bar plot | `da_publication` (shared) | existing | — | shared |
| 24–26 | STRING, g:Profiler / clusterProfiler enrichment | `enrichment` (shared) | existing | — | shared |
| 27 | study parameters | the analysis record | — | — | KIT-053 |
| 30 | report (`DIANN_report.rmd` in the script) | `DIANN_limpa_report.rmd` from the record (deviation: the package's dedicated template) | — | — | KIT-053 |

## The RUV-as-covariates block (chunk 19b), as the script does it

1. `Y_protein = t(protein matrix)`; `M = getRuvIIIReplicateMatrixHelper(design, sample_id_column = "Run", grouping_variable = "group")`.
2. `findBestNegCtrlPercentageCustom(obj, percentage_range = 1:30, num_components_to_impute = 5, ruv_grouping_variable = "group", ruv_fdr_method = "BH", ruv_qval_cutoff = 0.05, separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE)`: per percentage, `getNegCtrlProtAnova` → skip if fewer than 5 controls → `ruvCancor` → cap `K ≤ min(max_acceptable_k, n_controls − 1)` (adaptive `max_acceptable_k` from `calculateAdaptiveMaxK(sample_size)`) → `separation_score = calculateSeparationScore` → `best_k` = the K at the maximum of (all − control) canonical correlation → `composite_score = calculateCompositeScore`; the best percentage maximises the composite score. Differences from the package's `findBestNegCtrlPercentage`: the K cap, the max-difference k rule instead of `findBestKElbow`, the 1–30 range (package 1–20), `BH` (package `qvalue`); the package's complete-controls skip rule is not in the script.
3. `ruv_res = ruv::RUVIII(Y, M, ctl, k, return.info = TRUE)`; `alpha = fullalpha[1:k, ]`; `ac = alpha[, ctl]`; `W = Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))`; `W1..Wk` appended to `design_matrix`.
4. The quantified matrix is not changed; `t(ruv_res$newY)` is used for the QC plots only.
5. The DE formula becomes `config formula + " + W1 + … + Wk"`.

## What the script does not run (and the route omits)

Peptide intensity filter, minimum peptides per sample, peptide replicate filter, peptide imputation (replicate mean), IQ MaxLFQ, protein-level accession cleanup, protein intensity, duplicate and replicate filters, RUV-III-C correction of the matrix (`ruvIII_C_Varying`).

## Engine

`fastlimpadpc` (`APAF-bioinformatics/fast_limpa_dpc`, `main` at 2026-09-28; exports `dpc`, `dpcCN`, `dpcImpute`, `dpcQuant`, `dpcQuantByRow`, `peptides2Proteins`; imports `limpa::estimateDPCIntercept`; Cargo, rustc ≥ 1.65, xz; README: max |ΔE| 1.65e-8, max |ΔSE| 8.5e-5 against limpa on 50,000 × 16; Rayon-parallel per protein group). It does not provide `dpcDE`, which stays limpa's. Determinism across runs and thread counts is checked in KIT-047, not assumed; the engine and its version are recorded on every record.
