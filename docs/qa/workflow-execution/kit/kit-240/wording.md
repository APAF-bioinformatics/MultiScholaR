# KIT-240: user-facing texts (drafts awaiting the owner)

Every text below is a **draft awaiting the owner's wording** (VAL-41 allows implementing proposed wording; ticket
KIT-240 r1.4, OD-1 to OD-3; implementation audit AUDIT-KIT-240-impl-20261125). The values the marks describe are decided (VAL-221); only the words are open.
Ledger IDs in the marks refer to `docs/methods-decision-ledger.md` in the repository, which the installed application
does not carry (AUDIT-KIT-240-plan-20261125, Re-check 1, R1-07). The predicate `integrationParameterProvisional()`
marks a default as waiting when its mark contains "awaiting the owner": a reworded waiting mark must keep that token, or
the predicate changes with it (stop rule 7).

## OD-1: the Model tab when the software is missing (`R/utils_integration_gui.R`)

Each is followed by the probe's reason under "Technical detail". In the application the commands are plain text (no
code formatting).

| Key | Draft |
|---|---|
| `readiness_missing_mofa2` | The model cannot be trained on this computer: the R package MOFA2 is not installed. The MultiScholaR installers include it. From R, install it with BiocManager::install("MOFA2"), then press Check again. |
| `readiness_missing_reticulate` | The model cannot be trained on this computer: the R package reticulate, which MOFA2 uses to run Python, is not installed. From R, install it with install.packages("reticulate"), then press Check again. |
| `readiness_missing_mofapy2` | The model cannot be trained on this computer: no Python with the module mofapy2 was found. The MultiScholaR installers include one. From R: install into a Python the mofapy2 version your MOFA2 asks for (0.7.3 for MOFA2 1.22.0), for example pip install mofapy2==0.7.3; close MultiScholaR and R; set the environment variable RETICULATE_PYTHON to that Python; then start R and MultiScholaR again. |

## OD-2: README, NEWS and the testing guide

| Where | Draft |
|---|---|
| `README.md`, key features | Multi-Omics Integration: combine proteomics, metabolomics, lipidomics and transcriptomics processed in one project (MOFA2). |
| `README.md`, diagram | node "Multi-Omics Integration (MOFA2)", its edges from the three omics' QC and normalisation nodes (the integration reads the normalised matrices saved at each Summary tab) |
| `README.md`, table | Integration: MOFA2. Shared factors across two or more omics processed in the same project (the Integration tabs); for each data set, pathway enrichment of the features that follow a factor: STRING for protein and gene data sets (your own STRING key), KEGG and Reactome for metabolite and lipid data sets |
| `NEWS.md` | the entry under the unreleased heading, "New features" (KIT-240) |
| `docs/testing/data-testing-guide.md` section 2 | When it is on, the window title ends with **(test features on)** and each omic's Normalisation tab shows the RUV design choices. The Integration tile is on for everyone and needs no setting. |

## OD-3: the marks, the notes and the stability sentences

### Marks (`settled_by`; the record's "Settled by" column and `integrationSettings()`)

| Parameter | Draft |
|---|---|
| `center_groups` | VAL-221 (decided 2026-11-25: TRUE kept; with one sample group, each feature is centred over all samples) |
| `num_factors` | VAL-221 (decided 2026-11-25: 5 factors) |
| `seed` (model) | VAL-221 (decided 2026-11-25: 20261119, the seed of the reported model) |
| `training_seeds` | INT-14 (several seeds, 2026-11-21); VAL-221 (decided 2026-11-25: 5 seeds, the declared seed and the next four) |
| `technical_replicates` | INT-2 (decided 2026-11-19: averaged per biological sample) |
| `internal_standards` | INT-3 (decided 2026-11-19: dropped) |
| `missing_samples` | INT-010 (KIT-060 D4, decided 2026-11-17: kept as missing in that data set) |
| the 17 Results settings | VAL-221 (decided 2026-11-25: the default kept) |
| `string_fdr`, `string_rank_direction` | KIT-141 (decided 2026-11-20 as recorded best practice; methods ledger, the STRING ranked enrichment client) |
| `correlation_method`, `min_samples` | KIT-141 (chosen 2026-11-21 under INT-18; awaiting the owner's confirmation; methods ledger, the ranked list) |
| `kegg_organism`, `reactome_organism`, `min_set_size`, `max_set_size`, `seed` (enrichment) | KIT-141 (chosen 2026-11-21 under INT-19; awaiting the owner's confirmation; methods ledger, KEGG and Reactome rules) |
| `string_species` | KIT-141 (human, 9606, as the starting value; no decision recorded; awaiting the owner's decision) |

The other marks (`input_matrix`, `feature_selection`, `features_per_view`, `min_observed_fraction`, `scaling`,
`scale_views`, `convergence_mode`, `maxiter`, `precision`, and the enrichment `test`, `permutations`,
`adjust_ruv_factors`) are unchanged.

### Notes and labels

| Key | Where shown | Draft |
|---|---|---|
| `align_provisional` (r1.4, I03) | Sample Matching, under the rules | The averaging of runs, the removal of internal standards and the handling of missing samples are decided defaults (methods ledger INT-2, INT-3 and INT-010); the analysis record names the decision behind each. |
| `fact_sparse` (r1.4, I02; %s data set, %d features, the fewest samples, the samples kept) | Sample Matching, the rules list | %s: %d features measured in fewer than %d of the %d biological samples are left out. |
| `col_sparse` (r1.4, I02) | Sample Matching, the data-set table's column | Features left out (measured in too few samples) |
| `settings_provisional` | above the Model settings | Each default below is a decision recorded in the methods ledger, named in the analysis record. A value you change is recorded as set by you. |
| `mark_provisional` | in brackets after the label of each Enrichment setting whose default is not yet confirmed (the eight above: correlation, smallest number of samples, STRING species, KEGG organism, Reactome organism, smallest and largest pathway, seed); none of the Model or Results settings (r1.3, coordinator's ruling on R1-04) | default not yet confirmed |
| `results_provisional` | Results tab, above its settings and above Save | These defaults are decisions recorded in the methods ledger (VAL-221); a value you change is recorded as set by you. |
| `legend_provisional_default` (record; r1.4, I05) | the report's status legend (always printed) | a default not yet confirmed |
| `status_provisional_default` (record; r1.4, I05) | the record's status of such a setting, and the legend | default not yet confirmed |

`methods_provisional` keeps its text; it is shown only for a default that still waits, which no Model or Results
setting does.

### The stability sentences (VAL-221 (1); Re-check 1 R1-01)

| Key | Where shown | Draft |
|---|---|---|
| `stability_elbo_rank` | Results tab, after the stability introduction, when every seed's final ELBO is known | By the final value of MOFA2's training objective (the ELBO, a lower bound on the model evidence; higher is better), the model shown ranks %d of %d: %s below the highest. The final values of all seeds span %s. |
| `stability_limit` | Results tab, after the stability introduction, whenever two or more seeds were trained | The table describes only the factors of the model shown: a factor that other seeds find but this model lacks is not listed. |
| `methods_seeds_elbo` (record; r1.4, I09, I10) | the methods paragraph, after the seeds sentences, when the rank is known | By its final evidence lower bound (ELBO, MOFA2's training objective; higher is better), the reported model ranked %d of the %d models, %s below the highest (range over the models %s, in the objective's units). |
| `methods_seeds_limit` (record; r1.4, I09) | the methods paragraph, after the seeds sentence, whenever two or more seeds were trained | The stability figures describe only the reported model's factors: a factor that other seeds find but the reported model lacks is not listed. |

The two %s are ELBO differences written with two decimals (`integrationElboText`); a positive difference below 0.005
reads "less than 0.01" (r1.4, I01), so the rank line never prints "0.00 below the highest" at rank 2 or lower. No
interpretive clause is added unless the owner adds one with a stated basis.

### Open for the owner: how the rank is presented (AUDIT-KIT-240-impl I01)

The seeds' final ELBOs can differ by less than training's own stopping slack (KIT-151c n20weak cap1000 slow: five seeds
span 0.0146, the declared seed ranks 5 of 5; KIT-062: one solution stopped at 186 and 191 rounds ended 0.676 apart), so
the rank can order one solution. Options, none implemented:

- (a) recommended by the audit: one sentence after `stability_elbo_rank`, and the same clause in `methods_seeds_elbo`,
  with its basis in the ledger row: "Training stops by a rule, not at an exact optimum, so a small difference in the
  final value can come from where training stopped rather than from a better solution."
- (b) a per-analysis figure instead: the reported model's rise over its last ELBO evaluations, read from
  `value$elbo$trace`, for example "In its last evaluations the reported model's value still rose by %s." (a new field
  and a test).
- (c) drop the ordinal and keep the distance and the span: "the model shown is %s below the highest final value; the
  final values of all seeds span %s." (changes VAL-221 (1)'s "where ... ranks can be shown").
- or keep it as is.

## Two facts for the owner (Re-check 1, R1-07)

- "The owner" means nothing to a user, and `legend_provisional_default` is printed in every report: something like
  "a default not yet confirmed" reads plainly.
- The ledger IDs the marks cite are in the repository's `docs/methods-decision-ledger.md`, not in the installed
  application; the paper or README should say where to read them.
