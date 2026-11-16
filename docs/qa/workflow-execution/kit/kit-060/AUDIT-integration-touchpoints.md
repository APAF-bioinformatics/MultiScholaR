# KIT-060: integration touchpoints, the inputs contract, the GUI shape, the test corpus and the MOFA2 source check

Audit only. Nothing was run, installed or changed. Read from the tree at `015d6577` (branch `performance-v2`) on 2026-11-17.
The companion file `outcome.md` holds the decisions for the user, the re-estimate and the acceptance status.

Short names used below:

- "the workbook" is `inst/workbooks/integration/integration_workflow_starter.rmd` (3,374 lines, 30 code chunks). Line numbers `wb:N` point into it.
- "the project" is the folder the GUI sets up with `setupDirectories()` (`R/func_general_filemgmt.R:66-77`). In the GUI the label is `NULL`, so each omic's folders are named by the omic alone.
- "the record" is `scripts/<omic>/analysis_record.json`, written at the Summary tab (`R/utils_workflow_summary.R:112-118`).

Row ids: I (inputs), A (alignment), F (feature naming), M (model), DS (downstream), E (enrichment), C (corrected matrices). D1 to D15 are the decisions in `outcome.md`.

Closing tickets: KIT-061 (per-omic input), KIT-062 (owner, model, GUI import and model tabs, corpus), KIT-063 (downstream, enrichment, record, report, workbook rewrite). "User" means a decision listed in `outcome.md`. "Workbook" means the answer waits for the user's verified methods workbook.

## 0. What changed since the ticket's inventory of 2026-10-26

| Inventory statement | State today |
|---|---|
| Transcriptomics is on its own branch | On the state owner, with its own binding (`R/utils_workflow_trx_declarations.R:700-737`). |
| Lipidomics and metabolomics run in memory | On the state owner (`createWorkflowSessionOwner`, `R/utils_workflow_dia_sessions.R:875`; bindings in `utils_workflow_lipid_declarations.R:17`, `utils_workflow_metab_declarations.R`). |
| Nothing in the package mentions integration | Every omic's Summary tab already writes an "integration file" into `<project>/integration/` (`R/utils_workflow_summary.R:119-139`). The four files have four different shapes (section 2). The home page already has an Integration tile, disabled, with no module behind it (`R/app_ui.R:471-479`). |
| The enrichment helpers are untouched by the series | The STRING submission already records an external resource on the active owner (`R/func_multiomics_enrich_stringdb_api.R:169`); KEGG is guarded as an optional package (`R/func_multiomics_enrich_kegg_helpers.R:20`). |
| No synthetic corpus on shared samples | Still true. `tests/testdata/e2e/cross_omic_packs.json` names three "MOFA" packs, but they only list lanes whose samples are unrelated; its tests check the JSON's shape and never fit a model (`tests/testthat/test-e2e-cross-omic.R`). |

## 1. The project folder as it is today

```
<project>/
  data/<omic>/                      data_dir (user inputs); data/UniProt/ for proteomics
  scripts/<omic>/                   source_dir: design_matrix.tab, contrast_strings.tab, data_cln*.tab,
                                    config.ini, study_parameters.txt, analysis_record.json, session exports,
                                    the owner's store (.dia-native-*/.workflow-entities/)
  results/<omic>/<subfolders>       protein_qc, gene_qc, lipid_qc, metabolite_qc, da_*, publication_graphs, ...
  results_summary/<omic>/           QC_figures, Publication_figures, Publication_tables, Study_report
  integration/                      shared; each omic's Summary tab writes one file here
  results/integration/              multiomic_inputs, mofa/{inputs,plots,tables}, mixomics/{...},
                                    integration_enrichment/{inputs,plots,tables}   (created only when "integration" is set up)
  data/integration/, scripts/integration/   the same; the two starter workbooks are copied to scripts/integration/
```

Sources: `R/func_general_filemgmt_setupdirectories_helpers.R:43-135` (per-omic subfolders), `:153-166` (`integration_dir`), `:452-472` (roots); `R/app_server.R:326-345` (workbook copy).

## 2. What each omic exposes today

All four omics build their object in the GUI with the same identity columns: sample id `Run`, group `group`, technical replicates `tech_rep_group` (which is `<group>_<replicates>`; `R/func_omics_design_shared.R:683-688`). Sources: proteomics `R/utils_workflow_dia_state.R:1419-1421`, lipid `R/utils_workflow_lipid_declarations.R:230-231`, metab `R/utils_workflow_metab_declarations.R:215`, trx `R/func_trx_design.R:113-116`.

The design file is `scripts/<omic>/design_matrix.tab` in every omic (`R/func_general_applets.R:378-382`, `R/utils_workflow_omics_declarations.R:1103`). It is the design as built. Samples removed later by QC and the RUV `W` columns are only in the states; the record's `design` block carries the last assay state's design rows with `sample_id`, `group_id` and `technical_replicate_id` (`R/utils_workflow_analysis_record.R:304-311`).

### 2.1 Per omic and route

| Omic, route | Normalised matrix (RUV off) | Corrected matrix (RUV on) | Files in the project today | Feature id and annotation | Integration file written at Summary |
|---|---|---|---|---|---|
| Proteomics DIA (IQ rollup), LFQ, TMT ("classic") | State `normalized` (RUV skipped; `R/utils_workflow_dia_products.R:602`), then `correlation_filtered`. | State `ruv_corrected`, then `correlation_filtered`. Kernel: section 2.2, row C1. | `results/proteomics/protein_qc/normalized_protein_matrix_pre_ruv.tsv` (always; `utils_workflow_dia_products.R:641-642`). Final, after the correlation filter: `ruv_normalised_results_cln_with_replicates.tsv` when RUV ran, `normalised_results_cln_with_replicates.tsv` when it was skipped (`R/utils_workflow_dia_norm.R:645-656`, written at `utils_workflow_dia_products.R:1109-1112`). Wide TSV: the protein id column, then one column per run. | Id: the import's protein column (`protein_id_column`). `protein_id_table` holds only the distinct ids (`utils_workflow_dia_state.R:1417`). UniProt annotation: `scripts/proteomics/uniprot_dat_cln.RDS`; the query is recorded as an external resource (`R/func_prot_annotation_uniprot_helpers.R:207`). | `integration/proteomics_<label>_final.msdia.tar` on a native session (the whole session archive), else `proteomics_<label>_final_s4.RDS` (`R/utils_workflow_dia_summary.R:34-55`). `<label>` is the free text typed at the Summary tab. |
| Proteomics DIA_LIMPA | State `correlation_filtered` (the DPC-quantified, normalised matrix). | Two modes (`utils_workflow_dia_products.R:942-956`). Covariates (default): state `ruv_covariates`; the matrix is unchanged, `W1..Wk` join the design; the adjusted matrix is a product only. Correction: state `ruv_corrected`; the adjusted matrix becomes the state's table. Kernel: row C2. | The correlation file as above. Covariates step: `protein_qc/ruv_adjusted_protein_matrix.tsv`, `ruv_covariates_design.tsv`, `ruv_covariates_optimization.tsv` (`utils_workflow_dia_products.R:981-986`). | As classic. | As classic. The Summary takes `ruv_covariates` first (`utils_workflow_dia_summary.R:84-90`), whose table is the uncorrected matrix. |
| Transcriptomics | TMM logCPM (`edgeR::cpm(log = TRUE)`, `R/func_trx_normalisation.R:11-15`). A saved product of step `trx_tmm_logcpm` only: no state holds it and no file is written. RUV cannot be skipped (modes are automatic or manual; `utils_workflow_trx_declarations.R:234`). | Two different matrices: rows C3a and C3b. State `trx_ruvs` keeps the counts and gains `W_1..W_k` in its design (`utils_workflow_trx_declarations.R:268-281`). | `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv`: C3a, genes by runs with a `GeneID` column (`:557-558`). | `GeneID`, `GeneName`, both columns of the count table (`R/func_trx_s4_class.R:14-19`). | `integration/transcriptomics_normalized_data.tsv`: C3b, runs by genes with a `Run` column (`:643-648`). Refused until the user has run the separate "RUVIII Integration Export" step (`R/utils_omic_gui_config.R:2072-2074`). No label in the name. |
| Lipidomics | State `lipid_normalized` (log2, between-sample normalised; `R/utils_workflow_omics_declarations.R:973-975`), then `lipid_correlation_filtered`. | State `lipid_ruv_corrected` (`:981-983`), then `lipid_correlation_filtered`. Kernel: row C4. One table per assay. | No matrix file. `scripts/lipidomics/lipid_filtered_session_data_latest.rds` and a timestamped copy (a session list; `R/utils_workflow_saved_exports.R:171-190`), `lipid_ruv_optimization_results.RDS`. | Id: the mapped id column; annotation: the mapped annotation column; both are columns of each assay table (`utils_workflow_lipid_declarations.R:228-230`). | `integration/lipidomics_<label>_final_s4.RDS`: the whole S4 object of the first state found among `lipid_correlation_filtered`, `lipid_norm_complete`, `lipid_ruv_corrected`, `lipid_normalized`, ... (`utils_omic_gui_config.R:1618-1623`; written by the worker, `utils_workflow_omics_declarations.R:1390-1400`). |
| Metabolomics | State `metab_normalized`, then `metab_correlation_filtered`. | State `metab_ruv_corrected`, then `metab_correlation_filtered`. Kernel: row C4. One table per assay. | As lipid, prefix `metab_`. | As lipid; defaults `database_identifier`, `metabolite_identification` (`R/func_metab_s4_objects.R:83-84`). | `integration/metabolomics_<label>_final_s4.RDS`, as lipid (`utils_omic_gui_config.R:1835-1840`). |

### 2.2 "Corrected" names five matrices, not four

The ticket's revision 1.3 lists four. Transcriptomics has two, so the register lists five.

| Id | Where | What it is | Replicate structure | Levels | Where it lands |
|---|---|---|---|---|---|
| C1 | Classic proteomics. S4 method `ruvIII_C_Varying` for `ProteinQuantitativeData`, `R/func_prot_s4_norm_methods.R:204-287`; kernel `RUVIIIC::RUVIII_C_Varying` at `:255`; run from `R/utils_workflow_dia_products.R:365`. | RUV-III-C on the normalised matrix. | Technical replicates when the design has them, else the grouping variable (`:245-254`); recorded as `ruv_replicate_structure`. | Each protein's level restored after correction (`restoreRuvFeatureLevels`, `:272`, KIT-123). | State `ruv_corrected`; file `ruv_normalised_results_cln_with_replicates.tsv`. |
| C2 | DIA_LIMPA. `computeRuvCovariates`, `R/func_prot_norm_optimization_helpers.R:1174-1205`: `ruv::RUVIII` at `:1192`, `W` recomputed from centred controls, `adjusted <- t(Y - W %*% alpha)` at `:1202`. Operation `protein_ruv_covariates`, `utils_workflow_dia_products.R:786`. | RUV-III `newY` with the same fit, k and controls as the covariates mode. | As C1; recorded as `replicate_structure`, `block_reason`. | Kept (centred `W`). | Covariates mode: a product (`ruv-adjusted-science`, `:819`) and `ruv_adjusted_protein_matrix.tsv`; the state's matrix is uncorrected. Correction mode: state `ruv_corrected`. |
| C3a | Transcriptomics display matrix. `limma::removeBatchEffect(log_cpm_tmm, covariates = W, design = ~ group)`, `R/func_trx_ruv.R:152`; rebuilt for display at `utils_workflow_trx_declarations.R:505`. | TMM logCPM with the RUVs factors regressed out, the group term protected. "Visualisation only" in the code's own words (`func_trx_ruv.R:136-137`). | RUVs `scIdx` from the groups (`RUVSeq::makeGroups(design$group)`, `:145`). | Kept. | Saved table `corrected` of product `trx.ruvs_result`; `gene_qc/ruv_normalised_log_values.tsv`. |
| C3b | Transcriptomics RUVIII export. `ruv::RUVIII(Y = t(log_cpm_tmm), M, ctl, k)`, `R/func_trx_ruv.R:171`; operation `trx_ruviii_export`, `utils_workflow_trx_declarations.R:394-404`. | RUV-III `newY` on the TMM logCPM with the RUVs controls and k. Made for integration only; the DA does not use it. | The group by default, or technical replicates (`ruviii_m_source`; `func_trx_ruv.R:180-191`). | Not restored (plain `ruv::RUVIII`). | Product `trx.ruviii`, recorded under `snapshot$exports`; `integration/transcriptomics_normalized_data.tsv`. |
| C4 | Lipidomics and metabolomics. S4 method `ruvIII_C_Varying` for `LipidomicsAssayData` / `MetaboliteAssayData`, `R/func_lipid_s4_ruv_methods.R:397-400`, `R/func_metab_s4_ruv_methods.R:397-400`; kernel `RUVIII_C_Varying` at `:806` in both. | RUV-III-C per assay on the log2 normalised table. | The RUV grouping variable (a declared parameter). | Not restored today (no `restoreRuvFeatureLevels` call in either file; KIT-127 is pending). | State `<prefix>_ruv_corrected`. No file. |

Consequences for the input decision of revision 1.2:

- The workbook reads the matrix C3b for transcriptomics (`wb:338`, a file with a `Run` column) and a file named like C1 for proteomics (`wb:337`). Revision 1.3 names C3a for transcriptomics. The two differ in method and in what they keep. This is a decision (outcome D3).
- Transcriptomics has no "RUV off" route, and its normalised matrix is never written to a file.
- On DIA_LIMPA in covariates mode the corrected matrix is not the state the Summary exports.
- KIT-126 and KIT-127 (pending) move the classic proteomics and the lipid/metab DA to "factors in the model". The corrected matrix must then still be produced as a product for integration, as C2 is today. KIT-061 depends on how those two tickets land.
- A per-feature Z-scaling (what the workbook does) removes a constant level shift per feature, so the missing level restoration in C3b and C4 does not change a Z-scaled input. It does change any view that is not centred per feature.

## 3. The touchpoint register

Columns: where; what it reads or calls; owning omic; what the owner exposes today; the gap; closing ticket.

### 3.1 Inputs per omics

| # | Workbook chunk / function | Reads | Owner | Exposed today | Gap | Closes |
|---|---|---|---|---|---|---|
| I1 | Chunk 3, `wb:254-263` | `setupDirectories(omic_types = c(..., "integration"))`; every later path is `project_dirs[["integration_<label>"]]$data_dir`, i.e. `data/integration/`. | Integration | The folders exist (section 1). In the GUI the key is `integration` (no label). | The workbook reads every input from `data/integration/`. No omic writes there. The omics write to `<project>/integration/`. `results/integration/multiomic_inputs/` exists and is unused. Three candidate places, one contract needed. | KIT-061; user (D1) |
| I2 | Chunk 4, `wb:337` | `ruv_normalized_results_cln_with_replicates.tsv` (proteins by samples; column `Protein.Ids`). | Proteomics | `protein_qc/ruv_normalised_results_cln_with_replicates.tsv` (British spelling), or `normalised_...` when RUV was skipped. | The name the workbook reads is not the name the package writes. The id column is fixed as `Protein.Ids` in the workbook (`wb:603`); the package uses the import's mapped column. The file has no hash in the record. On DIA_LIMPA covariates mode this file is the uncorrected matrix. | KIT-061 |
| I3 | Chunk 4, `wb:338`; chunk 6, `wb:594-600` | `transcriptomics_normalized_data.tsv` (samples by genes, column `Run`), transposed by hand. | Transcriptomics | Exactly this file, at `<project>/integration/` (matrix C3b). | Orientation differs from every other view (samples in rows). Only written after a separate export step. It is C3b, not the C3a that revision 1.3 names. No RUV-off form. | KIT-061; user (D3) |
| I4 | Chunk 5, `wb:452-481` | `readRDS("workshop_data_final_metabolomics_s4_object.RDS")`; `@metabolite_data[[1]]`, `[[2]]`; columns `metabolite`, `database_identifier`, `metabolite_identification`, `chemical_formula`, `smiles`, `inchi`. | Metabolomics | `integration/metabolomics_<label>_final_s4.RDS`, the whole S4 object of the final state. | The integration would have to know the S4 class and which columns are samples. The workbook names the feature column `metabolite`; the package's id column is the mapped one (default `database_identifier`). The file does not say whether its state is corrected or only normalised; only the record does. Assays are taken by position (1 = LC, 2 = GC). | KIT-061 |
| I5 | Chunk 5, `wb:455-456`, `wb:470-471` | Drops rows whose `database_identifier` or `metabolite_identification` starts with `ITSD`. | Metabolomics | The object carries `internal_standard_regex`; the ITSD normalisation uses it. Internal-standard rows stay in the final tables. | Whether internal standards are removed from the integration input, and by which rule, is not declared. The workbook's rule is a fixed prefix. | KIT-061; workbook |
| I6 | (no workbook chunk) | Lipidomics is not read by the workbook at all. | Lipidomics | `integration/lipidomics_<label>_final_s4.RDS`. | No precedent in the workbook for lipid views (one per assay, which names). | KIT-061; workbook |
| I7 | Chunk 10, `wb:1258-1268` | `design_matrix.tab` (columns `Run`, `group`), from the proteomics design only; the group factor `Y` with levels `RPMI`, `Sera` fixed in code. | Proteomics (as the reference design) | `scripts/<omic>/design_matrix.tab` per omic; the final design in the record's `design` block. | One omic's design stands in for all. Group levels are fixed in code. The design file is the one before QC. | KIT-061 (per-omic design in the input), KIT-062 (integration design) |
| I8 | (implicit) | Sample identity column and original ids. | Each omic | `Run` in every omic; the original file names are the `Run` values unless the import sanitised them. | No column says which biological sample a run is, across omics (section 3.2). | KIT-061; user (D5) |
| I9 | (implicit) | The states' entity ids and the record. | Each omic | `analysis_record.json`: steps with `state`, `operation`, parameters, `result_entity_id`, the RUV outcomes (k, controls, replicate structure, mode), `analyses` (DA, the trx RUVIII export), `inputs` by role with sha256, `external_resources` (`R/utils_workflow_analysis_record.R:372-396`). | The record names entity ids, not the files in the project. Nothing ties an integration file's bytes to a state. The record is written only when the user saves at Summary. No field names "the integration input". | KIT-061 |
| I10 | (implicit) | File naming. | Each omic | Proteomics, lipid, metab: `<omic>_<label>_final...`; trx: a fixed name. | The label is free text, so a scan can find several files per omic and cannot tell which is current. | KIT-061; user (D2) |

### 3.2 Sample alignment

| # | Where | What it does | Exposed today | Gap | Closes |
|---|---|---|---|---|---|
| A1 | Chunk 5, `wb:457-458`, `wb:472-473`; chunk 6, `wb:600` | Removes the sample the proteomics QC dropped from the other views, by a sample id typed in the code: `51591` for metabolomics, `51519` for transcriptomics. | Each omic's final design lists its surviving runs. | The two ids differ. One of them may be a typing error; this audit did not open the demo data to check. The rule "a sample dropped by one omic is dropped from all" is not declared, and MOFA2 itself accepts missing values (section 7). | KIT-062; user (D4) |
| A2 | Chunk 5, `wb:483-492`; chunk 6, `wb:608-625` | Sorts each view's sample ids, then renames the columns by position to `R1..R6`, `S1..S5`. Alignment across views is "same position after sorting". | `Run`, `group`, `replicates`, `tech_rep_group` in every design. | Position after sorting is not a key. No shared biological-sample key exists across omics. `tech_rep_group` (`<group>_<replicates>`) could serve only if the user numbered replicates the same way in every omic. | KIT-062; user (D5) |
| A3 | Chunk 6, `wb:618-622`; chunk 10, `wb:1260-1262` | Builds `sample_names_lookup` (original id per view to the new id) and joins it to the proteomics design. | Nothing. | This is the integration design: one table mapping each view's run to the biological sample and the groups. It has no owner, file or hash today. | KIT-062 |
| A4 | (not in the workbook) | Technical replicates. | Trx, lipid and metab designs can hold technical replicates (`tech_rep_group`). | The demo has none, so the workbook has no rule. A model needs one column per biological sample per view: collapse, choose one, or refuse. A methods question. | KIT-062; workbook |

### 3.3 Feature naming

| # | Where | Reads | Owner | Exposed today | Gap | Closes |
|---|---|---|---|---|---|---|
| F1 | Chunk 7, `wb:741-756` | `Galaxy4-[Diamond on data 1_ Blast Tabular].tabular` (a DIAMOND BLAST table, RefSeq protein to UniProt). | External, organism-specific | Nothing. | A hand-made lookup for one organism. Not a general input. | KIT-063; user (D11) |
| F2 | Chunk 7, `wb:759-760`, `wb:786-787` | `idmapping_2025_05_08.tsv`, `idmapping_2024_10_07.tsv` (UniProt id mapping downloads). | External | Proteomics has its own UniProt annotation table with the query recorded. | Two dated downloads with no hash. The proteomics annotation the package already holds is not used. | KIT-061 (publish the annotation with the input), KIT-063 |
| F3 | Chunk 7, `wb:783` | `kv_uniparc.csv` (UniParc lookup). | External | Nothing. | As F1. | KIT-063; user (D11) |
| F4 | Chunk 8, `wb:947-1009` | `k_variicola.gff3`, parsed by hand into locus tag to gene name. | External | Transcriptomics carries `GeneName` beside `GeneID` in its count table. | The package's own annotation column would replace this for display names. The GFF3 is still needed where the workbook maps transcripts to protein ids for STRING (E2). | KIT-061, KIT-063 |
| F5 | Chunks 7 to 9, `wb:778-779`, `wb:810-811`, `wb:1008-1009`, `wb:1118-1120` | Writes and re-reads `prot_refseq_id_to_gene_name.tsv`, `uniprot_id_to_gene_name.tsv`, `transcript_id_to_gene_name.tsv` in `data/integration/`. | Integration | Nothing. | Derived tables written into the input folder, unhashed, re-read later by the enrichment (E1, E2). | KIT-063 |
| F6 | Chunk 9, `wb:1123-1140` | Replaces row names by display names; a duplicate gets a `_<n>` suffix (`wb:771-775`). | Integration | Nothing. | Display names become the model's feature ids. The stable id is lost inside the model, and the enrichment has to map names back (`wb:2443-2445`). A contract should keep the stable id as the model's feature id and carry the display name beside it. | KIT-062 |
| F7 | `plotMofaWeights`, `R/func_multiomics_mofa.R:112-113` | Strips `_transcriptome` and `_proteome` suffixes from feature names. | Package | n/a | MOFA2 adds a view suffix when feature names repeat across views. The function knows two view names only. | KIT-063 |

### 3.4 The model

| # | Where | Call | Exposed today | Gap | Closes |
|---|---|---|---|---|---|
| M1 | Chunk 9, `wb:1147-1150` | `t(scale(t(x)))` per view: each feature centred and scaled to unit variance across samples. Views: `transcriptome`, `proteome`, `metabolome_lc`, `metabolome_gc`. | Nothing. | The scaling is by hand, outside MOFA2, and unrecorded. A feature with zero variance becomes `NaN`. MOFA2's own documentation speaks of scaling views, not features (section 7). | KIT-062; workbook |
| M2 | Chunk 11, `wb:1329` | `create_mofa(mofa_matrix_input)`: a list of feature-by-sample matrices, one group. | Nothing. | The view names and their order are the list's; no declaration. | KIT-062 |
| M3 | Chunk 12, `wb:1538-1551` | `get_default_data_options`, with `center_groups <- FALSE` and `scale_views` left at its default (the line setting it is commented out); `get_default_model_options` with `num_factors <- 2`; `get_default_training_options` unchanged. | Nothing. | Three option sets, two changed values, the seed and convergence mode left at MOFA2's defaults and never written down. The number of factors is fixed for the demo. | KIT-062; workbook |
| M4 | Chunk 13, `wb:1629-1634` | `prepare_mofa(object, data_options, model_options, training_options)`. | Nothing. | n/a beyond M3. | KIT-062 |
| M5 | Chunk 14, `wb:1715-1738` | The training is commented out. The live line is `model <- readRDS("model.RDS")` from `data/integration/`. The commented code calls `run_mofa(..., use_basilisk = FALSE)` inside `try()`, then `run_mofa(..., use_basilisk = TRUE)`, writes `model.hdf5` and `model.RDS`. | `MOFA2` and `reticulate` are in `Suggests` (`DESCRIPTION:417`, `:434`); `basilisk` is in the dependency list of `loadDependencies` (`R/func_general_filemgmt.R:248`). | The shipped workbook does not train a model, so there is no run to be exact against. The Python runtime (`mofapy2`) is reached two ways with a deliberate failed first call. No probe, no stated absence behaviour. With basilisk, the first run downloads an environment. | KIT-062; user (D6, D12) |
| M6 | Chunk 14 | The trained model file. | The folder `results/integration/mofa/inputs/` exists (`mofa_inputs_dir`). | The model is an RDS of a MOFA object plus an HDF5 file, neither hashed nor recorded. By revision 1.3 a model trained in Python is compared numerically, not by hash. | KIT-062 |
| M7 | Chunk 15, `wb:1804-1811` | `samples_metadata(model) <- data.frame(sample, condition = Y)`. | Nothing. | The group labels are attached after training, by position. | KIT-062 |
| M8 | Chunk 29, `wb:3276` | `set.seed(42)` before `plot_dimred` (UMAP, t-SNE). | Nothing. | The only seed the workbook sets; it is for the display, not the model. | KIT-063 |

### 3.5 Downstream

| # | Where | Call | Output today | Gap | Closes |
|---|---|---|---|---|---|
| DS1 | Chunk 16, `wb:1897-1903`; chunk 17, `wb:2015-2072` | `model@cache$variance_explained` (`r2_total`, `r2_per_factor`), `plot_variance_explained`, two hand-made bar plots; `savePlot` to `mofa_plots_dir`. | PNG and PDF in `results/integration/mofa/plots/`. | The view order is fixed in code. The plot named "total variance" is drawn from the per-factor table (`wb:2061`), not from `variance_total_tbl`. Outputs are not recorded. | KIT-063 |
| DS2 | Chunk 18, `wb:2102-2103`; chunk 19, `wb:2192-2231` | `plot_variance_explained(plot_total = TRUE)`, `plot_factor`, `plot_factors` coloured by `condition`, colours fixed for `RPMI` and `Sera`. | PNG and PDF. | Group names and colours fixed in code. This is the plot the claim limit of revision 1.3 governs: factors against the design are descriptive. | KIT-063 |
| DS3 | Chunk 20, `wb:2294-2318` | `plot_weights(view, factor = 1, nfeatures = 10, scale = TRUE)` per view. | PNG and PDF. | Factor 1 only. | KIT-063 |
| DS4 | Chunk 25, `wb:2802-2807`; `plotMofaWeights`, `R/func_multiomics_mofa.R:97-136` | `MOFA2::get_weights(model)[[view]]`, top 20 by absolute weight. | Saved by the function itself through `getProjectPaths(omic_type, experiment_label)$mofa_plots_dir` (`:130`). | The function reads `omic_type` and `experiment_label` as free variables from the caller's environment and `project_dirs` from the global environment. Its bars and axis label use `Factor1` whatever `factor_level` is (`:118`, `:124`). It cannot run inside a module as it is. Kernels are not rewritten, so KIT-063 wraps or repairs it by decision. | KIT-063 |
| DS5 | Chunk 26, `wb:2877-2916`; chunk 28, `wb:3144-3208` | `plot_data_heatmap`, `plot_data_scatter` per view, factor 1. | Scatter plots saved; heatmaps shown only. | Views fixed in code. | KIT-063 |
| DS6 | Chunk 27, `wb:2924-3076` | A helper that searches a view's row names for a pattern (`dimethylglycine|betaine`) and draws box plots by group, the group taken from the first letter of the renamed sample id. | Shown only. | Demo-specific. | KIT-063; workbook |
| DS7 | Chunk 29, `wb:3280-3290` | `plot_dimred` (UMAP, t-SNE). | Shown only. | With eleven samples. | KIT-063; workbook |
| DS8 | Chunk 30, `wb:3355-3369` | `get_factors`, `get_weights`, `get_data` as long tables, written to `mofa_tables_dir` as `factors.tab`, `weights.tab`, `data.tab`. | Three TSV files. | Unrecorded. These three tables are the natural saved products of the model step. | KIT-062 (products), KIT-063 (export) |
| DS9 | (not in the workbook) | A report. | Each omic has a report template and the record renders to DOCX or HTML (`R/utils_workflow_analysis_record.R:765`). | No integration report template, no integration record. | KIT-063 |

### 3.6 Enrichment

| # | Where | Call | External resource | Gap | Closes |
|---|---|---|---|---|---|
| E1 | Chunk 21, `wb:2434-2475`; `runOneStringDbRankEnrichmentMofa`, `R/func_multiomics_enrich_stringdb_rank_helpers.R:163-205` | `get_expectations(model, "W")`, factor 1 of view `proteome`, joined back to RefSeq ids through F5, sent to STRING as a ranked list. | STRING API v12 (`https://version-12-0.string-db.org/api`, `R/func_multiomics_enrich_stringdb_api.R:114`), with an API key and a STRING species id for one organism. | An API key is written in the workbook in six places (`wb:2454`, `:2496`, `:2513`, `:2525`, `:2536`, `:2650`); it should be removed from the shipped file and never recorded. Network needed. Factor 1 only. The result is cached as an RDS and re-read. | KIT-063; user (D11) |
| E2 | Chunk 22, `wb:2633-2675` | The same for view `transcriptome`, mapped through the GFF3 lookup to protein ids. | STRING. | Transcripts without a protein id are dropped silently (`wb:2640`). | KIT-063 |
| E3 | Unnamed chunk, `wb:2489-2540`; `runStringDbEnrichmentFromDAResults`, `..._rank_helpers.R:273`; `runStringDbEnrichmentFromDAResultsMultiple`, `:441`; `runStringDbEnrichmentAllContrasts`, `R/func_multiomics_enrich_stringdb_reporting_helpers.R:340`; `plotStringDbEnrichmentResults`, `:571` | STRING enrichment from a proteomics DA result object, not from the model. | STRING. | The workbook passes `de_results_for_enrichment =`; the function's argument is `da_results_for_enrichment` (`:273`), so the chunk fails as written. It reads an S4 class `da_results_for_enrichment`, not the DA manifest the owner holds. Not an integration step. | KIT-063; workbook |
| E4 | Chunk 23, `wb:2719-2748`; `runMetabolomicsPathwayEnrichment`, `R/func_multiomics_enrich_metabolomics_helpers.R:162`; `runMetabolomicsEnrichmentAnalysis`, `:16` | Weights of views `metabolome_lc` and `metabolome_gc`, factor 1; ChEBI ids from the metabolomics S4 object's assay tables (`@metabolite_data[[i]]`, columns `metabolite`, `database_identifier`); `metaboliteIDmapping.txt`. | A mapping table file (from AnnotationHub, by the function's documentation); KEGG through `KEGGREST` (`R/func_multiomics_enrich_kegg_helpers.R:260`, `:583`); Reactome through `clusterProfiler::GSEA` with fgsea (`R/func_multiomics_enrich_reactome_helpers.R:302`). | Reads the metabolomics S4 object directly: the one place downstream that reaches into an omic's object. View names and assay positions are fixed in the function (`:50-53`, `:207`, `:223`). KEGG species (`kpn`) and Reactome organism fixed in the call. The mapping file is unhashed. GSEA is seeded nowhere in the workbook. Paths come from `getProjectPaths()` and the global `project_dirs`. | KIT-061 (publish the feature annotation with its database ids), KIT-063 |
| E5 | `runKeggEnrichment`, `R/func_multiomics_enrich_kegg_helpers.R:13`; `runReactomeEnrichment`, `R/func_multiomics_enrich_reactome_helpers.R:52`; `resolveReactomeOrganism`, `:12` | Called by E4. | KEGG REST, Reactome. | Neither records an external resource (STRING does). | KIT-063 |
| E6 | `submitStringDBEnrichment`, `R/func_multiomics_enrich_stringdb_api.R:31`; `retrieveStringDBEnrichmentResults`, `:447`; `downloadStringDBGraph`, `:259`; `downloadStringDBResultsFile`, `:327`; `runOneStringDbRankEnrichment`, `..._rank_helpers.R:53` | Called by E1 to E3. | STRING. | The submission records the resource on `getDIAActiveManager()` (`:169`): on whichever proteomics owner is active, not on an integration owner. Polling with fixed waits. | KIT-063 |
| E7 | `printStringDbFunctionalEnrichmentBarGraph`, `R/func_multiomics_enrich_stringdb_reporting_helpers.R:3`; chunk 24, `wb:2755-2772` | Bar graphs per view and combined. | n/a | The combined plot labels the comparison by view name in code. | KIT-063 |
| E8 | `getStringDbSpecies`, `..._reporting_helpers.R:80`; `searchStringDbSpecies`, `:235` | Species lookup. | STRING. | Not called by the workbook; the GUI would need it to replace the fixed species id. | KIT-063 |

Coverage check (the ticket's simulation setup): every `readRDS`, `vroom`, `create_mofa`, `prepare_mofa`, `run_mofa`, `get_*` and plotting call of the workbook has a row (I2 to I4, I7, F1 to F5, M1 to M8, DS1 to DS8, E1 to E4), and each of the 19 top-level functions in `R/func_multiomics_*.R` is named in DS4 or E1 to E8 (`plotMofaWeights` and 18 enrichment helpers; the ticket's inventory counted 19 helpers). The second file, `inst/workbooks/integration/MultiOmics_workflow_starter.rmd` (72 lines), only installs the package and calls a function named `multisetupAndShowDirectories`; it touches no omics result.

## 4. The per-omics inputs contract

Written against the folder layout of section 1 and the hook that already exists: each omic's Summary tab writes its integration file through `config$summary$integration` (`R/utils_workflow_summary.R:152-172`), on a managed session through the binding's `integration_export(manager, parameters, path)` (`R/utils_workflow_omics_declarations.R:1390-1400`). KIT-061 would change what is written there, not add a second mechanism.

### 4.1 The product

One folder per omic, with a fixed name and no free-text label. Its place is decision D1; the names below assume `<project>/integration/<omic>/`.

| File | Content |
|---|---|
| `integration_input.json` | The manifest (fields in 4.2). The integration import reads this file first and nothing else by guesswork. |
| `matrix_<view>.tsv` | One per view. Features in rows, samples in columns. First column `feature_id` (the omic's stable id), then one column per run, named by the omic's `Run` values. Values on the log2 scale as the omic's final state holds them. Missing values stay missing. No annotation columns. |
| `samples.tsv` | The final design of the state the matrix came from: `Run`, `group`, `replicates`, `tech_rep_group`, any factor and batch columns, the RUV `W` columns when the state has them, and the biological-sample key (decision D5). Only the runs that are in the matrix. |
| `features_<view>.tsv` | One row per feature of the view: `feature_id`, `display_name`, and the omic's annotation columns (for proteomics the UniProt accession and gene name; for transcriptomics `GeneName`; for lipid and metab the annotation column and any database identifiers such as ChEBI or KEGG). Internal standards flagged in a column, not removed (I5). |

Plain tables are proposed because the integration must not know an omic's S4 class or store (the ticket's point 4), and because the four files written today are a tar archive, two S4 objects and one TSV in a different orientation. Whether today's files stay beside the new product is decision D2.

### 4.2 The manifest's fields

| Field | Meaning |
|---|---|
| `schema` | `multischolar.integration_input.v1`. |
| `omic_type`, `workflow` | The omic and its route (`DIA`, `DIA_LIMPA`, `LFQ`, `TMT`, `TRANSCRIPTOMICS`, `LIPIDOMICS`, `METABOLOMICS`), from the record. |
| `views` | One entry per view: `name`, `matrix_file`, `features_file`, `rows`, `columns`, `sha256` of each file. Lipid and metab: one view per assay, named by the assay. |
| `matrix_role` | `corrected` or `normalised`. Named by role, not as "the DA input" (revision 1.2). |
| `correction` | When `matrix_role` is `corrected`: the method id (one of the rows C1 to C4 of section 2.2), `k`, the control count, the replicate structure (`technical_replicates` or the grouping variable, and its name), the grouping variable's levels. This is "the RUV design of every view" that revision 1.2 asks the integration record to name. |
| `scale` | `log2` (or `log2_cpm` for transcriptomics); whether per-feature levels were kept by the correction. |
| `source_state` | The state's name and `result_entity_id`, as the record lists them. For a product that is not a state (C2 in covariates mode, C3a, C3b): the product's entity id and type. |
| `store_identity` | `project_id`, `workflow_id` from the record. |
| `sample_id_column`, `group_column`, `technical_replicate_column`, `biological_sample_column` | The column names in `samples.tsv`. |
| `feature_id_column`, `display_name_column` | The column names in the features file. |
| `samples_file`, with `sha256` | The design. |
| `analysis_record` | The relative path of the omic's `analysis_record.json` and its sha256 at the time of writing. |
| `software` | Package version and build, from the record. |
| `written_at` | Time of writing. |

### 4.3 One corrected matrix per omic

| Omic, route | RUV on: the matrix of the input | RUV off: the matrix of the input | Open |
|---|---|---|---|
| Proteomics classic | C1 after the correlation filter (state `correlation_filtered` whose parent is `ruv_corrected`). | State `correlation_filtered` whose parent is `normalized`. | After KIT-126, C1 must remain available as a product. |
| Proteomics DIA_LIMPA | C2 restricted to the samples of the final state. In correction mode this is the state's table; in covariates mode it is the product `ruv-adjusted-science`. | State `correlation_filtered`. | The Summary's final state on this route is `ruv_covariates`, whose table is uncorrected; the writer must take the product. |
| Transcriptomics | C3a or C3b: decision D3. | TMM logCPM. Not reachable today (no skip mode); needed only if RUV becomes optional for transcriptomics. | D3; omics parity for "RUV off". |
| Lipidomics, metabolomics | C4 after the correlation filter, one view per assay. | `<prefix>_normalized` after the correlation filter. | KIT-127 (levels kept) changes C4's values. |

### 4.4 Record fields each omic adds

- In the omic's own record: an `analyses$integration_input` entry (the same mechanism that lists the trx RUVIII export, `R/utils_workflow_analysis_record.R:351`), with the manifest's path and sha256, `matrix_role` and `source_state`.
- The native session export carries the folder as members with their hashes (the proteomics export already lists its files, `R/utils_workflow_dia_sessions.R:1608-1645`).

### 4.5 What the contract does not settle

The handling of technical replicates (A4), feature filtering before the model, and whether internal standards are dropped are methods questions. The contract carries what is needed to apply any answer (the replicate column, the flag column) and applies none.

## 5. The GUI shape (design notes only)

### 5.1 How an omic is hosted today

- The home page lets the user pick omics; `values$selected_omics` drives everything (`R/app_server.R:126`).
- For each selected omic the app calls `mod_<omic>_ui(id)` and `mod_<omic>_server(id, project_dirs, omic_type, experiment_label, volumes)` by name (`R/app_server.R:437`, `:576`), after `setupDirectories()` has made its folders and copied `inst/workbooks/<omic>/` to `scripts/<omic>/` (`:326-355`, `:526-565`).
- The four `mod_<omic>_*` functions are thin wrappers over one shell, `mod_omic_workflow_ui/server(id, ..., config = omicGuiConfig(<omic>))` (`R/mod_omic_workflow.R:32`, `:71`; wrapper example `R/mod_transcriptomics.R`).
- `omicGuiConfig()` accepts exactly four omics (`R/utils_omic_gui_config.R:100-110`) and its validator requires the groups `shell`, `binding`, `import`, `design`, `qc`, `da`, `summary` (`norm` and `enrichment` may be `NULL`; `:120-140`).
- The shell builds a pill tabset from `shell$tabs`, a stepper, the per-session bus, the state owner (`binding$owner`), the runner, and the tab chain that unlocks the next tab (`R/mod_omic_workflow.R:104-166`). Tab modules come from a fixed registry (`omicGuiComponents`, `:210-232`) and are called through one of eight call shapes (`:262-277`).
- An Integration tile is already on the home page, disabled, enabled by JavaScript once two omics are selected, with the comment that no module backs it (`R/app_ui.R:471-479`, `:553`, `:570-582`). The side menu already has an icon for it (`R/app_server.R:409`). `setupDirectories()` already accepts `integration`.

### 5.2 What an integration workflow would look like in that frame

Tabs, in order: Import (scan), Alignment, Model, Results, Enrichment, Summary & Export. The science of Model, Results and Enrichment waits for the workbook; the tabs are named here only to size the frame.

- **Import.** No file upload. On entry it scans the project for each omic's `integration_input.json` (section 4), reads the manifests, verifies the hashes, and lists what it found: omic, route, views with their sizes, `matrix_role`, the correction's replicate structure, when it was written, and whether the omic's record has changed since. The user ticks the views to include. Two or more are needed to continue (decision D7). An omic that is processed but has not written its input is listed with the instruction to save at its Summary tab. This is the form the user asked for on 2026-11-17: data already processed through MultiScholaR, from the project's folders.
- **Alignment.** Shows the integration design: one row per biological sample, one column per view with that view's run, the groups. Built from the `samples.tsv` files by the biological-sample key. Shows which samples are missing in which view and applies the declared rule (D4). The user confirms. This replaces chunks 5, 6 and 10 of the workbook.
- **Model.** The options the workbook sets, a run button, the Python capability's state. A long step, so it runs on the managed worker like the other omics' heavy steps.
- **Results, Enrichment, Summary.** As KIT-063 defines them.

### 5.3 What can be reused

| Piece | Where | Reuse |
|---|---|---|
| Dispatch by name, folder setup, workbook copy | `R/app_server.R:326-365`, `:437`, `:526-600` | As is: defining `mod_integration_ui` and `mod_integration_server` is enough for the app to host it. The watchdog passes the whole `project_dirs` (all omics' paths), which the scan needs. |
| The home tile and menu entry | `R/app_ui.R:471-479`, `R/app_server.R:409` | Enable the tile (remove the placeholder state). The existing JavaScript rule ties it to "two omics selected in this session"; with the new form the rule should be "two processed omics in the project", which the home page cannot know before a project is chosen. A small change, listed as new. |
| Tab panel, stepper, tab chain, initial status | `omicWorkflowTabPanel`, `render_workflow_stepper`, `installOmicWorkflowChain`, `omicWorkflowInitialBus`, `omicGuiTab` | As is; they read only the `shell` group. |
| The shell functions `mod_omic_workflow_ui/server` | `R/mod_omic_workflow.R:32-102` | Not as is: both call `validateOmicGuiConfig`, which demands `design`, `qc`, `da` and an `import$source_panel` of `feature_table` or `proteomics`. Either the validator gains a second, smaller schema for a workflow that is not a single-omic pipeline, or `mod_integration_*` composes the shell-level helpers itself with a config validated by `validateOmicGuiShell` only. The second is the smaller change and leaves the four omics' contract untouched. |
| The state owner and the worker | `createWorkflowSessionOwner(workflow_data, experiment_paths, session, omic_type)`, `installWorkflowOmicsActions` | The owner needs only `experiment_paths$source_dir`, which exists for integration. Its documentation lists proteomics, lipidomics and metabolomics as the `omic_type`; transcriptomics already passes its own. An `integration` owner is new but follows the same call. |
| Operations with saved products | `newWorkflowOperation` (`R/utils_workflow_operations.R:30`), as proteomics declares `protein_ruv_covariates` | The model step, the alignment and the scaling fit this form. `newOmicsOperation` and the omics binding assume one S4 assay class with a data slot; integration has several views, so it either gets a small S4 class of its own or declares its operations directly. KIT-062's design decides. |
| Publishing files with the action | `publishWorkflowActionOutputs` (`R/utils_workflow_saved_exports.R:424`) | For plots and tables written to `results/integration/mofa/`. |
| Capability probes | `requireOptionalPackage` (`R/utils_optional_packages.R:6`), the `workflowToolPath` pattern (`R/utils_workflow_saved_entities.R:678`) | `MOFA2` and `reticulate` by the first; the Python module `mofapy2` needs a new probe in the second's style (memoised, can be disabled for tests). |
| External resources in the record | `recordWorkflowExternalResource` (`R/utils_workflow_dia_state.R:2929`) | For STRING, KEGG, Reactome and the mapping files. STRING's existing call must point at the integration owner. |
| The record and its rendering | `buildWorkflowAnalysisRecord` and the rest of `R/utils_workflow_analysis_record.R` | The step, parameter, product and software parts apply. The `inputs` list by role with sha256 is where the consumed integration inputs go. The `design` block assumes one design; the integration design fits. `analyses` assumes DA manifests and exports; the model's outputs would be a new entry. |
| The Summary tab | `mod_omic_summary_ui/server` (`R/mod_omic_summary.R:37`, `:162`) | Its config group requires `final_state`, `copy`, `integration`, `report` (`R/utils_omic_gui_config.R:336-350`). The save, record and report flow is reusable; the "integration export" part has no meaning here. Reuse depends on the validator question above. |

### 5.4 What is new

1. `mod_integration_ui/server` and its shell config (tabs, bus, status keys).
2. The Import tab: the project scan, manifest validation, the view list. No existing tab reads another omic's folders.
3. The Alignment tab and the integration design table.
4. The Model tab and its declarations (scaling, `create_mofa`, `prepare_mofa`, `run_mofa`), the model file as a saved entity, the three long tables as products.
5. The Python capability probe and its absence text.
6. Results and Enrichment tabs over functions that today read global variables (`plotMofaWeights`, the metabolomics enrichment helpers).
7. An integration report template (`inst/reports` has one per omic; none for integration).
8. Test ids, a module-CI lane and an e2e lane for the new tabs, on the corpus of section 6.

## 6. A synthetic multi-omics corpus on shared samples with known shared factors

None exists. This is its specification; KIT-062 builds it before its first run.

### 6.1 The generators that exist

| Generator | File | Writes | Design | What varies |
|---|---|---|---|---|
| Lipidomics, metabolomics | `tools/profiling/kit_omics_corpus.R`: `kitOmicsRecipe` (`:30`), `kitOmicsDesign` (`:53`), `kitOmicsAssayMatrix` (`:65`), `writeKitLipidCorpus` (`:103`), `writeKitMetabCorpus` (`:161`) | LipidSearch and MS-DIAL exports, `design.tsv`, `manifest.json` with sha256 | Groups `WT`, `KO`; biological samples each measured twice; three batches; `tech_rep_group` is the biological sample | Group effect on 15% of features (KO), a batch effect with a per-feature loading, biological noise per biological sample, sample scale, technical noise; all drawn inside the function from the assay's seed |
| Transcriptomics | same file: `kitTrxRecipe` (`:198`), `kitTrxCounts` (`:216`), `writeKitTrxCorpus` (`:249`) | `counts.tsv`, `design_matrix.tab`, `truth.json`, manifest | Groups `RPMI`, `Sera`; the same replicate and batch layout | Gamma-Poisson counts: batch effect on 30% of genes, a Sera effect on 10%, a gamma multiplier per gene and biological sample |
| Proteomics LFQ, TMT | `tools/profiling/kit_protein_corpus.R`: `kitProteinRecipe` (`:14`), `kitProteinMatrix` (`:33`), `writeKitLfqCorpus` (`:75`), `writeKitTmtCorpus` (`:118`) | MaxQuant `proteinGroups.txt`, Proteome Discoverer export, FASTA, manifest | Groups `WT`, `KO`, alternating; batches of ten consecutive samples; no technical replicates | Group effect on 10% of proteins, batch, sample scale, noise, missing values |
| Proteomics DIA (and the e2e lanes) | `tests/testdata/e2e/generate_fixtures.R`: `protFixtureMatrix` (`:138`), `generateProtDiannFixture` (`:197`) | A DIA-NN `report.tsv` expanded from `kitProteinMatrix` | Samples `WT_1..`, `KO_1..` | As the protein recipe |

Why they cannot be put side by side as they are: each draws its own biological noise from its own seed, so nothing but the group is shared; the group labels differ (`WT`/`KO` against `RPMI`/`Sera`, and the effect is tied to the literal `KO` or `Sera` in the code, `kit_omics_corpus.R:77`, `:237`); the sample names and the replicate layout differ; no generator takes an effect from outside.

### 6.2 The corpus

**One design, shared.** `kitOmicsDesign()` with one recipe gives the biological samples (`tech_rep_group`), groups and batches. Proposed size: two groups, ten biological samples per group (twenty biological samples, above the "more than 15" that MOFA2's FAQ asks for; section 7). Technical replicates as each generator has them: two for lipid, metab and trx, none for proteomics. This difference is kept on purpose: it exercises the technical-replicate rule (A4).

**Known factors, per biological sample.**

| Factor | Scores | Active in | Purpose |
|---|---|---|---|
| Z1 | The group (0 or 1) | Every view, on each generator's declared differential features | The factor every view shares and the design explains. |
| Z2 | Continuous, standard normal, drawn once per biological sample, independent of the group | Every view, on a declared share of features per view (proposed 10%), with declared loadings | A shared biological factor the design does not explain: what an integration should find beyond the groups. |
| Z3 | Continuous, drawn once per biological sample | Two views only (proposed proteomics and transcriptomics) | A partly shared factor: the variance-explained table has a known pattern of views. |
| Batch | Each omic's own batch assignment, different per omic | One view each | Technical, not shared. With RUV on it should be reduced; with RUV off it should appear as a single-view factor. |

**How the generators take the factors.** Each needs one additive, default-preserving argument: a matrix of effects on the log2 scale (features by samples), added before the values are rounded (lipid, metab, protein) or before the Poisson draw (trx). With the argument absent the existing corpora regenerate byte for byte, the rule KIT-086 and KIT-099 followed for their additions. The effect group's label also becomes a recipe field with today's value as default, so one pair of group names serves every omic. These are changes to test tooling under `tools/` and `tests/testdata/`, made by KIT-062, not by this audit.

**A composer** (new, proposed `tools/profiling/kit_integration_corpus.R`) that: builds the design; draws Z2 and Z3 and the loadings from one declared seed; calls the four generators with the same biological samples and the effect matrices; gives each omic its own run names (as real exports have) and writes the mapping; removes one biological sample from one omic's files (the QC-drop case, A1); writes `truth.json` (factor scores per biological sample, loadings per feature and view, the active views per factor, the differential and batch feature sets, the dropped sample, every recipe) and a manifest with a sha256 per file.

**Two tiers.**

- Tier A, direct: the composer also writes each omic's integration input in the contract's form (section 4) straight from the generated log2 matrices, with `matrix_role = "normalised"` and no correction. Small and fast; the integration owner's unit tests (scan, alignment, scaling, model call, record) run on it without running any omic's journey.
- Tier B, through the journeys: the raw files are run through each omic's journey to the Summary export in one project, with the kit driver's existing routes (`kitJourneySteps` accepts `DIA`, `protein`, `DIA_LIMPA`, `lipidomics`, `metabolomics`, `transcriptomics`; `tools/profiling/workflow_kit.R:12-15`), once with RUV on and once with RUV off. This is the corpus for the import scan, for the one e2e lane, and for comparing factors on corrected against uncorrected inputs. Four small journeys per variant; sizes as the e2e fixtures (300 to 400 features per view).

**What the truth supports.**

- Recovery of Z2 and Z3 is measured by the absolute correlation between a recovered factor and the known scores, and by the canonical correlations between the recovered and the known factor spaces. Factors have no fixed order or sign, so nothing is compared by factor index. The threshold is set from a pilot on Tier A and written into KIT-062 before its first counted run; this audit sets none.
- The pattern of views per factor in the variance-explained table against the declared pattern.
- The claim limit of revision 1.3 holds on the corpus too: on a matrix corrected with the groups as replicate structure, Z1's separation of the groups is not a test result. The corpus can show the effect the caveat describes (compare Z1 on RUV-on and RUV-off inputs), and that comparison is reported as a description.
- The corpus is synthetic and small. It tests that the plumbing recovers what was put in. It supports no statement about performance on real data.

## 7. MOFA2 source check

Two pages were fetched on 2026-11-17 and read as text.

- FAQ: <https://biofam.github.io/MOFA2/faq.html>
- Bioconductor vignette "MOFA2: training a model in R" (MOFA2 1.22.1): <https://www.bioconductor.org/packages/release/bioc/vignettes/MOFA2/inst/doc/getting_started_R.html>

### 7.1 Removing technical variation before fitting

The FAQ says it, in these words:

> (Q) Should I remove undesired sources of variability (i.e. batch effects) before fitting the model?
> Yes. If you have clear technical factors, we strongly encourage to regress it out a priori using for example a linear model like limma. This is very important for the feature selection step. If you do not remove technical variability then MOFA will "focus" on capturing this huge variability driven by the technical factors, and smaller sources of (biological) variability could be missed.

So the sentence "remove technical factors first" has a source and can stay, cited to the FAQ. What the source does and does not cover:

- It speaks of "clear technical factors" regressed out with a linear model, limma as the example. It does not mention RUV, estimated factors, or negative controls.
- It says nothing about a correction whose replicate structure is the experimental groups. The caveat and the claim limit of revisions 1.2 and 1.3 are the project's own and are not contradicted or supported by the page.
- The vignette does not mention batch effects or technical variation at all.

### 7.2 Normalisation and scaling per view

The FAQ:

> (Q) How do I normalise the data?
> Proper normalisation is critical for the model to work. First, one needs to remove library size effects. For count-based data such as RNA-seq or ATAC-seq we recommend size factor normalisation + variance stabilisation. If this is not done correctly, the model will learn a very strong Factor 1 that will capture differences in the total expression per sample, and more subtle sources of variation will be downweighted.

> Importantly, the use of non-gaussian likelihoods require statistical approximations and are not as accurate as the gaussian likelihood. If your data can be safely transformed to match the gaussian likelihood assumptions, this is ALWAYS recommended. For example RNA-seq data is expected to be normalised and modelled with a gaussian distribution, do not input the counts directly.

The vignette, section 3.1 and the data options:

> Proper normalisation of the data is critical. [...] Non-gaussian likelihoods give non-optimal results, we recommend the user to apply data transformations to obtain continuous measurements. For example, for count-based data such as RNA-seq or ATAC-seq we recommend size factor normalisation + variance stabilisation (i.e. a log transformation).

> scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance. Default is FALSE

> scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance. Default is FALSE

Neither page says to scale each feature to unit variance. The workbook's `t(scale(t(x)))` (M1) is therefore the workbook's own choice, not MOFA2's documented practice, and it leaves `scale_views` at `FALSE`. It is recorded here as what the workbook does; whether it stays is for the user's verified workbook.

### 7.3 Other statements that bear on the design

| Topic | What the FAQ says | Bearing |
|---|---|---|
| Feature filtering | "It is strongly recommended that you filter highly variable features (HVGs) per assay. When doing multi-group inference, you have to regress out the group effect before selecting HVGs." | The workbook filters nothing. |
| View sizes | "Bigger data modalities will tend to be overrepresented in the factors. It is good practice to filter uninformative features (based for example on a minimum variance threshold) in order to have the different views within the same order of magnitudes." | A transcriptome of thousands of genes beside a few hundred metabolites. |
| Sample count | "Factor Analysis models are only useful with large sample sizes, at least more than 15." | The demo has eleven samples (decision D7). |
| Missing values | "Yes. It simply ignores them from the likelihood, there is no hidden imputation step." | Dropping a sample from every view because one omic lost it (A1) is not required by the model (decision D4). |
| Known covariates in the model | "We recommend that you learn the factors in a completely unsupervised manner and then relate them to the biological covariates a posteriori". | Agrees with the claim limit's direction: factors are related to the design afterwards, descriptively. |
| Multi-group mode | "The aim of the multi-group framework is not to capture differential changes in mean levels between the groups [...] If your aim is to find a factor that 'separates' the groups, you DO NOT want to use the multi-group framework." | The workbook uses one group, so this mode is not in play. If the GUI ever offers groups, this warning goes with it. |
| The Python runtime | Vignette: "in this step the MOFA2 R package connets with the mofapy2 Python package using reticulate. This is the source of most problems when running MOFA." And the training log: "Connecting to the mofapy2 package using basilisk. Set 'use_basilisk' to FALSE if you prefer to manually set the python binary using 'reticulate'." | The capability question (decision D6). |
| Seed and run-to-run differences | The pages as rendered do not address them. (The FAQ's HTML source has a commented-out answer on rotation invariance and on initialisation; it is not published text and is not relied on here.) | The seed is recorded from the training options; exactness is numerical (revision 1.3). |

## 8. Other things seen, recorded and left untouched

- An API key for STRING is written in the shipped workbook (lines in E1).
- `wb:2061` draws the "total variance" plot from the per-factor table.
- `wb:2490`: the argument name does not match the function (E3).
- `plotMofaWeights` ignores `factor_level` in its bars and label (DS4).
- The two hard-coded dropped-sample ids differ (A1).
- `inst/workbooks/integration/MultiOmics_workflow_starter.rmd` calls `multisetupAndShowDirectories`, which is not defined under `R/` nor listed in `NAMESPACE` (grep of both).
- `results/integration/mixomics/` is created by `setupDirectories()`; nothing in the workbook or the package uses it.
