# KIT-128 W0: technical replicates handled the same in every omic (design and authoring evidence)

Read-only. Written 2026-11-17 on `performance-v2` at `934e9de9`. No code changed, no test or simulation run.
R was used only to print installed function sources and help pages. Nothing here is measured unless a file is named.

**Bounded item.** Outcome: the facts and the smallest design KIT-128 needs before it starts. Owner: `tickets/kit/KIT-128.md`
(revision 1.1), parked rows 23, 24, 25 (`parked-20261113.md`). Files touched: this document only. Exit: the seven
sections below.

**Read first.** `AGENTS.md`; `conventions/R.md`, `R-shiny.md`, `R-golem.md`, `R-bioinformatics.md`; the ticket; the three rows.

**Sources read for section 3** (the installed copies and the current Bioconductor release copies carry the same revision dates):

- limma User's Guide, "Last revised 19 May 2026" (limma 3.68.4 installed):
  <https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf>
- edgeR User's Guide, "Last revised 25 April 2026" (edgeR 4.10.1 installed):
  <https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf>
- limpa vignette (limpa 1.4.0 installed; 1.4.2 is the release):
  <https://bioconductor.org/packages/release/bioc/vignettes/limpa/inst/doc/limpa.html>
- Help pages of the installed packages: `edgeR::sumTechReps`, `edgeR::voomLmFit`, `limma::duplicateCorrelation`.

---

## 0. The one finding that changes the ticket's size

Every synthetic corpus of the three omics declares technical replicates: `tools/profiling/kit_omics_corpus.R:12-16, 22, 36, 202`
("twelve biological samples each measured twice"; `technical_replicates = 2L` for lipid, metab and trx), and the design files
it writes carry the biological replicate number in `replicates` (:152, :189, :269), which the design step turns into a
duplicated `tech_rep_group`. So the change is not confined to a new fixture: **every stored lipid, metab and trx reference
taken on those corpora from the DA step on (trx: from the design state on) will move.** Section 4.5 lists them. The
ticket's "no re-record here (KIT-129)" therefore needs a rule for the time between the two tickets (decision D6).

---

## 1. Proteomics today

### 1.1 The column and how it is declared

- The design builder keeps `replicates` (the biological replicate number within the group) and `tech_reps`. Assigning
  technical replicates gives the chosen runs one shared `replicates` number (`applyOmicsDesignTechRepAssignment`,
  `R/func_omics_design_shared.R:629-652`).
- At save, the design gets `tech_rep_group = paste(group, replicates, sep = "_")`:
  proteomics `R/utils_workflow_dia_state.R:1238` (builder) and `:1341-1342` (import); the shared builder
  `R/func_omics_design_shared.R:684-688`; the shared import `R/utils_workflow_omics_design.R:286-287`.
- The object declares it in the slot `technical_replicate_id = "tech_rep_group"` (`R/utils_workflow_dia_state.R:1421`, `:1515`).
  The class default is `"replicates"` (`R/func_pept_s4_core.R:44`).
- A design without technical replicates has a unique `tech_rep_group` per run; with them, runs of one biological unit share it.

### 1.2 One rule, two readers

`resolveProtTechnicalReplicateBlock(design_matrix, group_id, technical_replicate_id)`,
`R/func_prot_norm_optimization_helpers.R:399-415`:

- not declared or column absent: `no_declaration`; any NA: `missing_values`;
- block = the column's value, prefixed by the group when a value recurs across groups (:409-410);
- no duplicated block: `no_duplicates`;
- every group is one block: `spanned_by_design` (not technical; limma could not estimate a correlation);
- else `technical = TRUE`, reason `technical_replicates`.

The function reads only a data frame and two column names. Nothing in it is specific to proteins.

### 1.3 (a) The DA

- Classic route (DIA with IQ, LFQ, TMT): `differentialAbundanceAnalysisHelper` calls `runTestsContrasts(...,
  group_id = theObject@group_id, technical_replicate_id = theObject@technical_replicate_id)`,
  `R/func_prot_da_model_methods.R:408-417` (the two names at :411). `runTestsContrasts`
  (`R/func_prot_da_model_stats.R:419-431`) builds the model matrix from the formula and hands the names on at :480-482
  to `fitDIAExecutionBaseModel` (`R/utils_workflow_dia_da.R:74-106`), which calls `fitProtDaBaseModel`
  (`R/func_prot_da_model_stats.R:227-266`).
- The kernel (:234-252): the block from the rule; when `technical`,
  `limma::duplicateCorrelation(data_subset, design = design_m, block = block)` (:246-247), then
  `limma::lmFit(data_subset, design = design_m, block = block, correlation = consensus.correlation)` (:251-252).
  Otherwise plain `lmFit` with the reason in a message (:253-260). The artifact records `block`, `block_reason`,
  `uses_duplicate_correlation`, `consensus_correlation` (:261-266). Then `contrasts.fit`, `eBayes`, `treat` as before
  (`applyProtDaContrastModel`, :283-312).
- DIA_LIMPA route: the same rule, the block passed to `limpa::dpcDE(elist, design_for_dpcde, plot = FALSE, block = block)`
  (`R/func_prot_da_model_methods.R:388-392`); `dpcDE` forwards `...` to `voomaLmFitWithImputation`, which runs
  `duplicateCorrelation` itself. The record: `block_reason`, `uses_duplicate_correlation` (:404-406).
- Test: `tests/testthat/test-tech-reps-limma.R`.
- Measured (KIT-124 pilot, 20 runs, cell C6: 4 units per group, 2 injections each, 75% of the residual variance
  biological; arm 2, the true factor in the model): null rate at 0.05 = 0.049 classic (`kit-124/s1-pilot/s1-runs.tsv`),
  0.053 DIA_LIMPA (`kit-124/s2-pilot/s2-runs.tsv`). One configuration, a pilot.

### 1.4 (b) RUV-III's replicate structure

`resolveProtRuvReplicateBlocks(design_matrix, sample_id, group_id, technical_replicate_id, ruv_grouping_variable)`,
`R/func_prot_norm_optimization_helpers.R:417-433`: when the rule says `technical`, the block becomes the column
`.ruv_replicate_block` and the structure is `technical_replicates`; else the RUV grouping variable, structure `group`.
Readers: `ruvIII_C_Varying` (`R/func_prot_s4_norm_methods.R:244-253`; M from `getRuvIIIReplicateMatrixHelper(...,
.ruv_replicate_block)`, the structure recorded in `ruv_replicate_structure`) and `computeRuvCovariates`
(`R/func_prot_norm_optimization_helpers.R:1187-1191`, returned at :1204). Test: `tests/testthat/test-prot-ruv-replicate-structure.R`.

---

## 2. The other three omics today

All three declare the same column: `technical_replicate_id = "tech_rep_group"` on the design object
(lipid `R/utils_workflow_lipid_declarations.R:231`; metab `R/utils_workflow_metab_declarations.R:216`;
trx `R/func_trx_design.R:113-117`). The column is in the design of every state. It is dropped at the DA kernel's door.

### 2.1 Lipidomics (kernel: `limma::lmFit`, `contrasts.fit`, `eBayes(trend, robust)`, optional `treat`)

| Step | Where | What happens to the column |
|---|---|---|
| All-contrast orchestrator | `runLipidsDA`, `R/func_lipid_da_model.R:372`; slots read at :396-401 | reads `sample_id`, `group_id`; never `technical_replicate_id` |
| Kernel wrapper | `runTestsContrastsLipidDA`, :220-229, call at :470-479 | no block argument |
| Fit half | `fitLipidDaAssayModel(data, design_matrix, formula_string, sample_id_col)`, :6-36 | `limma::lmFit(data_subset, design = design_m)` at :35 |
| Per-contrast declaration (the GUI's path) | `omicsDaContrastOperation`, `R/utils_workflow_omics_declarations.R:445-503`; `omicsDaAssayFit`, :418-438; binding `da_fit`, `R/utils_workflow_lipid_declarations.R:90` | `binding$da_fit(data, design, formula_string, sample_id)` at :424: the design is passed whole, the two names are not |
| Fit identity | `newWorkflowOmicsDaFitIdentity`, `R/utils_workflow_omics_declarations.R:405-411` | hashes values, formula, design; no block |
| RUV-III | `R/func_lipid_s4_ruv_methods.R:413` reads the slot "for context"; M at :690-696 from `ruv_grouping_variable` only | read, not used |

### 2.2 Metabolomics (the same kernel; the twin files)

| Step | Where |
|---|---|
| Orchestrator | `runMetabolitesDA`, `R/func_metab_da_orchestration.R:155`; slots at :179-184; kernel call at :244-253 |
| Kernel wrapper and fit half | `runTestsContrastsMetabDA`, `R/func_metab_da_core_engine.R:221`; `fitMetabDaAssayModel`, :6-36, `limma::lmFit(data_subset, design = design_m)` at :35 |
| Per-contrast declaration | the shared one above; binding `da_fit`, `R/utils_workflow_metab_declarations.R:78` |
| RUV-III | `R/func_metab_s4_ruv_methods.R:413` (read), M at :690-696 from the grouping variable only |

Not the GUI's path, recorded only: the legacy S4 method `differentialAbundanceAnalysisHelper` for `MetaboliteAssayData`
already passes the declaration to the proteomics kernel (`R/func_metab_s4_da_methods.R:162-171`, KIT-035); the lipid one
(`R/func_lipid_s4_da_analysis_methods.R:159-167`) passes nothing, so the kernel's default column `replicates` applies.

### 2.3 Transcriptomics (kernel: edgeR `calcNormFactors(TMM)`, `estimateDisp(robust)`, `glmQLFit(robust)`, `glmQLFTest`)

| Step | Where | What happens to the column |
|---|---|---|
| Design object | `trxDesignObject`, `R/func_trx_design.R:113-117`; class slot `R/func_trx_s4_class.R:15-20` | declared |
| Every kernel's input | `trxAssayInputs`, `R/func_trx_s4_class.R:32-48`; `validateTrxDesign`, `R/func_trx_design.R:15-41` | reads `Run` and `group` only; one column per library |
| RUV grid | `trxRuvGrid`, `R/func_trx_ruv.R:16-36` | `~ group`; no replicates |
| RUVs | `trxRuvs`, `R/func_trx_ruv.R:140-155`: `RUVSeq::makeGroups(design$group)` at :145 | replicate sets = the group, always |
| DA design | `trxDaDesign`, `R/utils_workflow_trx_declarations.R:290-293`; formula fixed to `~ 0 + group + W_1..W_k` (:298-301, refused otherwise :346-350) | not used |
| DA fit | `trxDaFit`, :323-331 (`estimateDisp` :329, `glmQLFit` :330); the kernel form `trxDaContrasts`, `R/func_trx_da.R:60-83` (:74-75) | each library is a sample |
| RUVIII export (integration only) | `trxRuviii` and `trxRuviiiM`, `R/func_trx_ruv.R:165-191`; operation `R/utils_workflow_trx_declarations.R:394-403`; choice offered by `trxHeadHasTechnicalReplicates`, `R/utils_workflow_summary.R:458-472` | the only reader: M from technical replicates when the user picks it |

`glmQLFit` takes no block or correlation argument (ticket revision 1.1; consistent with the edgeR guide, which has no
random-effect or `duplicateCorrelation` section).

---

## 3. The statistically correct treatment, from the guides

### 3.1 Log-intensity data with limma (proteomics, lipidomics, metabolomics)

limma User's Guide, section 11.3 "A Correlation Approach to Technical Replication":

> "The first and second and third and fourth arrays are technical replicates. It would not be correct to treat this
> experiment as comprising four replicate arrays because the technical replicate pairs are not independent, in fact they
> are likely to be positively correlated."

> "One way to analyze these data is the following:
> `biolrep <- c(1, 1, 2, 2)`
> `corfit <- duplicateCorrelation(MA, ndups = 1, block = biolrep)`
> `fit <- lmFit(MA, block = biolrep, cor = corfit$consensus)`
> `fit <- eBayes(fit)`"

> "The vector biolrep indicates the two blocks corresponding to biological replicates. The value corfit$consensus
> estimates the average correlation within the blocks and should be positive. This analysis is analogous to mixed model
> analysis of variance [...] except that information has been borrowed between genes."

> "Beware however that treating technical replicates as biological gives p-values that are smaller than they should be."

Section 9.7 "Multi-level Experiments" gives the same two calls with a design matrix
(`duplicateCorrelation(eset, design, block = targets$Subject)`, then `lmFit(eset, design, block = targets$Subject,
correlation = corfit$consensus)`). Section 2.1 names the reference for "duplicate spots or technical replication":
Smyth, Michaud and Scott (2005), Bioinformatics 21, 2067-2075.

`?duplicateCorrelation`: "Typically the blocks are biological replicates and repeated observations on the same block may
be correlated." and "For this function to return statistically useful results, there must be at least two more arrays
than the number of coefficients to be estimated".

limpa vignette, "Differential expression": "The `dpcDE()` function accepts any argument that `voomaLmFitWithImputation`
or `voomaLmFit` do. [...] Or `dpc(y.protein, design, block=subject)` could be used to model the correlation between
repeated observations on the same subject." (The vignette writes `dpc(`; the argument is `dpcDE`'s.)

**So:** block = the biological unit, `duplicateCorrelation` then `lmFit(block, correlation)`. This is what proteomics
does, and it is the documented treatment for lipid and metab log-intensities fitted with limma. Limits to state, not
to repair: one consensus correlation for all features; the residual degrees of freedom are not reduced.

### 3.2 Count data (transcriptomics)

edgeR User's Guide, section 4.9.1 (case study 4.9, "Time course RNA-seq experiments of Drosophila melanogaster"):

> "Each biological sample was sequenced several times and we view these as technical replicates."

> "The data has 127 columns of counts but only 30 biologically independent samples [...] We therefore sum the technical
> replicates to get total genewise read counts for each biological sample:
> `PooledCounts <- sumTechReps(Counts, ID=Samples$stage)`"

The summing comes first: the `DGEList` (4.9.2), filtering and normalisation (4.9.4) and the dispersion (4.9.7) all follow it.

limma User's Guide, section 18.2.5 "Assemble DGEList and sum counts for technical replicates":

> "Then we collapse the data so that there is a single column for each GEO sample by summing the counts over the
> technical replicates: `y <- sumTechReps(y.all, ID=SRA$GEO)`"

`?sumTechReps` (edgeR): "Condense the columns of a matrix or DGEList object so that counts are summed over technical
replicate samples."

The only other mention in the edgeR guide is section 2.13 (no replicates), which gives a typical BCV of "0.01 for
technical replicates"; it is not a method for this case.

The documented alternative that keeps one column per library is not an edgeR QL fit: limma guide section 15.7,
"we now recommend voomLmFit as the preferred interface for RNA-seq analyses [...] voomLmFit also automates the
estimation of sample weights and intrablock correlations", and `?voomLmFit` (edgeR): "`block`: vector or factor
specifying a blocking variable on the samples. [...] Samples within each block are assumed to be correlated." and
"If `block` is specified, then the intra-block correlation is estimated using `duplicateCorrelation`".

**So:** for the route's edgeR QL kernel the guides' treatment is to sum the technical replicates' counts per biological
sample before filtering, normalisation and the fit. A block for counts exists only in `voomLmFit` (limma-voom), which
would replace the kernel for those designs. Nothing else is proposed here.

---

## 4. The smallest change per omic

Common to all three: reuse `resolveProtTechnicalReplicateBlock` as it is (same rule, same reasons, no rename, no move;
`AGENTS.md`: reuse existing helpers). A design the rule calls `no_duplicates`, `spanned_by_design`, `no_declaration` or
`missing_values` keeps today's computation exactly.

### 4.1 Lipidomics and metabolomics: the block in the fit half

1. `fitLipidDaAssayModel` and `fitMetabDaAssayModel` gain `group_id = NULL, technical_replicate_id = NULL`. After the
   existing subsetting to `common_samples` (:29-31), resolve the block on the design rows of those samples; when
   `technical`, `duplicateCorrelation(data_subset, design = design_m, block = block)` and `lmFit(..., block = block,
   correlation = consensus)`; else the current line :35. The returned `MArrayLM` carries `block_reason`,
   `uses_duplicate_correlation`, `consensus_correlation` as extra components (limma keeps unknown components through
   `contrasts.fit`/`eBayes`; to be asserted by the test). With the defaults (`NULL`) the function is unchanged, so direct
   callers and their tests do not move.
2. Callers pass the two slots: `runTestsContrastsLipidDA` (:220-241) and `runLipidsDA` (:470-479);
   `runTestsContrastsMetabDA` (`func_metab_da_core_engine.R:221-242`) and `runMetabolitesDA` (:244-253).
3. The per-contrast path: `omicsDaAssayFit` (`utils_workflow_omics_declarations.R:418-438`) and its call at :473 pass
   `object@group_id` and `object@technical_replicate_id`; the lipid and metab bindings' `da_fit` already forward `...`;
   the trx binding's `da_fit` (`utils_workflow_trx_declarations.R:726`) names its arguments and ignores the new ones.
4. The fit identity (`newWorkflowOmicsDaFitIdentity`, :405-411) must include the block (or the rule's version), or a
   store holding a fit saved before the change would return it as a hit for the same values, formula and design.
5. The model record: `accepted$model` (:390-391 and the per-contrast manifest) gains `block_reason` and
   `uses_duplicate_correlation`, as the proteomics record has.

Files: `R/func_lipid_da_model.R`, `R/func_metab_da_core_engine.R`, `R/func_metab_da_orchestration.R`,
`R/utils_workflow_omics_declarations.R`. No new file in `R/`, no new dependency (limma is imported), no NAMESPACE change
(`duplicateCorrelation` is already imported, `R/func_general_package_bindings.R:7`; called as `limma::`).

### 4.2 Transcriptomics: sum at the design object (recommended; decision D1)

1. In `trxDesignObject` (`R/func_trx_design.R:113-117`), the one function both the builder and the import use to make the
   design state's object: resolve the block on the builder's design; when `technical`, the count table's run columns
   are replaced by `edgeR::sumTechReps(counts, ID = block)` and the design by one row per block (`Run` = the block
   value, `group`, `replicates`, `tech_rep_group`; first occurrence order, as `sumTechReps` orders its columns). The
   runs that were summed are recorded in `args` (unit to runs) so the record and the report can say it.
2. Nothing downstream changes: the filter, TMM, the grid, RUVs, the DA and the figures already work on whatever columns
   the object holds. After summing no replicate set has two runs, so `trxHeadHasTechnicalReplicates` is FALSE and the
   RUVIII export offers the group only (decision D2).
3. To check at W1, not assumed here: readers of the builder's library-level design on the bus after the design step
   (`workflowData$design_matrix`), and the design tab's summary text, which should say how many runs became how many
   samples.

Files: `R/func_trx_design.R` (and its test). edgeR is already required by every trx step (`requireOptionalPackage`).

Why not inside the DA kernel only: RUVs returns one W per library and the DA formula needs W per fitted column;
carrying W from libraries to summed units has no documented rule. Both guides sum before filtering and normalisation.

### 4.3 Tests (failing first)

- **Fixture** (one per omic, seeded, generated in the test, no stored file): two groups, 4 biological units per group,
  2 injections or libraries per unit, no true effect, residual variance mostly biological (KIT-124's C6 recipe:
  `kit-124/sim/s1_classic.R:50-52`, `s4_lipid_metab.R:69`, `s3_trx.R:48-49`), `tech_rep_group` naming the unit.
- **Limit**: null rate at 0.05 over all features and seeds at most 0.075 (ticket revision 1.1). KIT-124's rule 1 is
  written on the one-sided 97.5% upper bound; the test should apply it the same way and state its size (for example
  20 seeds by 500 features). Before the change the same test fails: 0.137 measured for lipid and metab
  (`kit-124/s4-pilot/s4-runs.tsv`, C6, arm 2, full), 0.179 for trx (`kit-124/s3-pilot/s3-runs.tsv`, C6, arm 2).
- **Lipid, metab**: through `runLipidsDA` / `runMetabolitesDA` and through the per-contrast operation (both paths must
  agree, as `test-workflow-omics-da-contrasts.R:19` already requires); plus: a design with unique `tech_rep_group` gives
  a fit identical to today's; a spanned design records `spanned_by_design` and fits as today; the record carries the reason.
- **Trx**: the object built from a duplicated design equals `edgeR::sumTechReps` of the counts with the collapsed design;
  a design without duplicates gives the object it gives today; then the null rate through the route's kernels in order
  (filter, TMM, grid, RUVs, DA).
- **Runtime**: limma and edgeR on 8 to 16 columns; seconds per seed. No long run, no worker journey.
- Affected regressions (targeted files, then the dependents grep, `R CMD check`, install, per the test rule):
  `test-lipid-02-da-core-helpers.R`, `test-metab-02g-*`, `test-metab-da-core-engine-direct-shared.R`,
  `test-metab-da-orchestration-direct-shared.R`, `test-metab-02h-*`, `test-workflow-omics-da-contrasts.R`,
  `test-workflow-lipid-declarations.R`, `test-workflow-metab-declarations.R`, `test-module-ci-lipid-da.R`,
  `test-module-ci-metab-da.R`, `test-trx-kernels.R`, `test-workflow-trx-declarations.R`, `test-module-ci-trx-design.R`,
  `test-tech-reps-limma.R`, `test-prot-ruv-replicate-structure.R`.
- To check at W1: unit fixtures that reuse replicate labels inside a group (`tests/testthat/helper-module-ci-lipid-da.R:10`,
  `helper-module-ci-metab-da.R:10`: `replicates = rep(c("R1", "R2"), ...)`) become technical replicates only if their
  object's declared column is built from them.

### 4.4 Interaction with RUV (for sequencing with KIT-126 and KIT-127)

- **Classic proteomics (KIT-126)**: no interaction. Both readers already take the rule.
- **Lipid, metab RUV-III replicate structure**: today the grouping variable only (`func_lipid_s4_ruv_methods.R:690-696`,
  metab the same lines). Parity is one call: `resolveProtRuvReplicateBlocks(design_matrix_filtered, sample_id, group_id,
  technical_replicate_id, ruv_grouping_variable_final)` and `.ruv_replicate_block` in place of the grouping column, with the
  structure recorded. It is result-changing for every state after RUV on a design with technical replicates. KIT-127
  rewrites the same method (levels kept, factors in the model). **Recommended: the DA block here, the RUV replicate
  structure in KIT-127**, so the method is edited once and its outputs move once (decision D3).
- **Order**: KIT-128 before KIT-127. Both edit the lipid and metab fit halves (KIT-128 adds the block, KIT-127 adds W
  to the model); the block and W covariates combine the way the DIA_LIMPA path already does (design with W, block
  passed; `func_prot_da_model_methods.R:382-392`). Not in parallel.
- **KIT-124's C6 cells for trx, lipid and metab** were all measured with the unblocked test
  (`kit-124/sim/s4_lipid_metab.R:128-129`: "the 16 injections of C6 are tested as 16 independent samples in every arm,
  the ceiling included"). Arms 15 and 16 there (technical replicates as the RUV replicate sets) read 0.16 and 0.137 for
  that reason, so they say nothing yet about the replicate structure. Any gate conclusion drawn from C6 in S3 or S4
  needs those cells re-measured after KIT-128.
- **Trx RUVs**: after summing there are no technical-replicate sets, so RUVs' replicate sets stay the group over
  biological units (`func_trx_ruv.R:145`). The rule "technical replicates when declared" is met by the summing; it has
  no second application in RUVs (decision D2).

### 4.5 Stored references that would move

Only on designs the rule calls `technical_replicates`. All three kit corpora are such designs (section 0).

| Omic | Moves | Does not move |
|---|---|---|
| Proteomics | nothing | everything |
| Lipid, metab | the DA step and everything after it on `lipid_screen` / `metab_screen` journeys: t, P, q, BH, B, the significance call and counts, the saved `da_fit` entities (new identity), the DA manifests' model record, DA figures, exports, the report text naming the model. Expected and to be confirmed at the re-record: logFC and AveExpr unchanged on a balanced complete design. | import, QC, normalisation, RUV, correlation filter (until KIT-127 changes the RUV replicate structure); kernels called directly without the two names |
| Trx | everything from the design state on, on `trx_screen` corpora: the design state (rows, sample names), filter counts, TMM/logCPM, grid, RUVs, DA, figures, exports. The KIT-040 notebook references (`docs/qa/workflow-execution/kit/kit-040/reference-<corpus>.rds`) no longer describe the route on these corpora: the unchanged notebook treats each library as a sample. `test-workflow-trx-declarations.R` (against those references) and `test-module-ci-trx-design.R:44-57` (the committed design's partition) would fail. | kernels on matrices (`test-trx-kernels.R`); any design without technical replicates |

The frozen declaration fixtures (`tests/testthat/fixtures/kit073-omics-declarations-frozen.rds`,
`kit094-characterisation.rds`) were not opened; whether the model record's two new fields touch them is a W1 check.

---

## 5. Where the kernel error is swallowed (row 24 b)

Three sites, one behaviour: a failure of the assay's fit is printed and logged, and the assay is dropped.

| Site | File:line | What it does |
|---|---|---|
| `runLipidsDA` | `R/func_lipid_da_model.R:468-515`; the handler at :510-514 | `cat("[ERROR] ...")`, `log_error`, `per_assay_results[[assay_name]] <<- NULL`; the loop goes on; `bind_rows` of nothing returns no rows (:519) |
| `runMetabolitesDA` (the twin) | `R/func_metab_da_orchestration.R:242-292`; the handler at :287-291 | the same |
| `omicsDaAssayFit` (the per-contrast path the GUI runs; lipid, metab and trx) | `R/utils_workflow_omics_declarations.R:424-429`, and `if (is.null(fit)) next` at :474 | the same for the fit half; its header (:413-415) records it as "the orchestrator's silent drop, kept for the fit half only" (KIT-076) |

The contrast half on the per-contrast path already does the right thing: the error is recorded under
`accepted$errors[[assay]]` (:480-486) and the DA tab warns (`test-workflow-omics-da-contrasts.R:68`, :146).

Smallest repair (recommended; decision D4):

- `omicsDaAssayFit` returns the message instead of NULL (`list(error = conditionMessage(e))`), and the caller writes it
  to `accepted$errors[[assay]]` before `next`: the existing channel, the existing warning.
- `runLipidsDA` and `runMetabolitesDA` keep the loop, collect `failed_assays` (assay and message) in the returned list,
  raise one `warning()` naming them, and stop with the messages when every assay failed (today: an empty table).
  `omicsDaOperation` (:373-392) copies `failed_assays` into its manifest.
- A classed condition and `cli` wording per `conventions/R.md` (errors: "NEVER silently swallow"); the surrounding
  legacy `cat`/`message` lines are left as they are (existing-code scope rule).
- Test, failing first: an assay with a feature that has no value in one group (row 24's case) returns a named failure
  on all three sites; an all-failed run stops.

---

## 6. Row 25: `dpcDE` on complete data, "No residual degrees of freedom in linear model fits"

From the code only (limpa 1.4.0 installed; the same lines are in the 1.4.2 source tarball, `R/dpcDE.R:7-8`,
`R/voomaLmFitWithImputation.R:136-139` and :234-237). Not run.

- Our call: `limpa::dpcDE(elist, design_for_dpcde, plot = FALSE, block = block)`
  (`R/func_prot_da_model_methods.R:391-392`), on the EList `dpcQuant` returned, subset to the model's samples. These are
  the documented arguments. The message itself is raised later by `limma::eBayes` in our converter
  (`R/func_prot_limpa_da_helpers.R:52`), on a fit that is already empty.
- `dpcDE` passes `predictor = log(y$other$standard.error + 1e-6)` and `imputed = !y$other$n.observations`. With no
  missing value, `imputed` is dropped (`max(NImputed) <= 1L`).
- In `voomaLmFitWithImputation`: `sxc <- rowMeans(predictor)`; `vartrend <- lm.fit(cbind(1, sx, sxc), sy)`;
  `beta <- coef(vartrend)`; `mu <- beta[1] + beta[2]*mu + beta[3]*predictor`. When the standard errors are one constant,
  `sxc` is collinear with the intercept; `lm.fit` pivots it out and returns `NA` for `beta[3]`; `mu` becomes all `NA`;
  the weights `1/f(mu)^4` are all `NA`; `lmFit` treats an observation with a non-finite weight as missing
  (`limma::lm.series`: `M[!is.finite(weights)] <- NA`); every fit has zero residual degrees of freedom; `eBayes` stops
  with exactly this message. There is no guard on `beta` in limpa.
- "2 of 5 runs, by floating-point noise" fits this: `lm.fit` drops the column only when the collinearity is inside its
  tolerance (1e-7); when the standard errors differ in the last digits the coefficient is finite (and arbitrary) and the
  run completes.
- Why the standard errors were constant in the pilot: every protein has 4 precursors with one noise sd
  (`kit-124/sim/s2_limpa.R:52, 57`), no value is missing (cell S5), and `dpcQuant`'s per-protein sigma is squeezed
  towards a common value (`dpcQuantHyperparam`: `squeezeVar`). The posterior standard error on complete data depends
  only on sigma, the precursor count and the priors.

**Conclusion: upstream, in limpa (present in the current release 1.4.2), not in our call.** It needs a degenerate
predictor: no missing values, the same precursor count for every protein, homogeneous variances. Real DIA data has
varying precursor counts, so the standard errors vary even when nothing is missing. Per the ticket it is recorded,
not repaired here. A guard in our call (a clear message when the standard errors have no spread) would be new behaviour
on our side; it is listed as decision D5, not proposed as part of this ticket.

---

## 7. Decisions for the user, and the estimate

### Decisions

| # | Decision | Recommendation |
|---|---|---|
| D1 | **Trx method.** (a) Sum the technical replicates' counts per biological sample at the design step, then the route as it is (edgeR guide 4.9.1, limma guide 18.2.5). (b) Keep libraries and replace the test by `voomLmFit(block =)` for designs with technical replicates (limma guide 15.7): a second test engine, other result columns. | (a) |
| D2 | **Trx consequences of (a).** On such designs the route no longer equals the unchanged notebook on library counts (the reference becomes the notebook rendered on summed counts, at the re-record); sample names become the unit names; RUVs' replicate sets are the group over biological units; the RUVIII export loses its "technical replicates" M. | accept |
| D3 | **Lipid and metab RUV replicate structure** (technical replicates as RUV-III's replicate sets): in KIT-128, or in KIT-127 which rewrites the same method. | KIT-127 |
| D4 | **Kernel error behaviour.** A failed assay is named in the result and warned about, the other assays are kept, and the run stops only when every assay failed; or the run stops at the first failed assay. | named and kept; stop when all fail |
| D5 | **Row 25.** Record only; and whether you want it reported to the limpa authors (your account, your call), and a guard message on our side in a later ticket. | record; report upstream; guard later if it ever shows on real data |
| D6 | **Between KIT-128 and KIT-129.** The lipid, metab and trx references on the kit corpora fail from the moment KIT-128 lands until the re-record (section 0). Hold KIT-128 on its branch until KIT-129's re-record, or land it with those reference tests listed as known-moved. | hold on branch `kit-128`, merge with the re-record |

Not a decision, to fix when the ticket is authored: the ticket's estimate (10 h) predates section 0; and KIT-124's
C6 cells in S3 and S4 need re-measuring after this ticket (section 4.4).

### Estimate

| Step | Content | Hours |
|---|---|---|
| W1 | Authoring: files, declarations, budgets; the three W1 checks named above (bus design readers in trx, frozen fixtures, unit fixtures with repeated labels) | 1.5 |
| W2 | Lipid and metab: failing-first fixture and null-rate test; the block in both fit halves, the callers, the per-contrast path, the fit identity, the model record | 4 |
| W3 | Kernel error surfaced at the three sites, with its test | 1.5 |
| W4 | Trx: failing-first tests; summing at the design object, the record of summed runs, the design summary text | 4 |
| W5 | Row 25 recorded in the parked file with section 6's reading | 0.5 |
| W6 | Checks: targeted files, dependents grep, `R CMD check`, install; outcome record | 2.5 |
| | **Total** | **14** |
| optional | Lipid and metab RUV replicate structure, if D3 puts it here | +2.5 |
| not here | Re-record of the moved references (KIT-129) | 0 |

### Convention notes

- `resolveProtTechnicalReplicateBlock` keeps its name and file though three omics will call it (no rename, `AGENTS.md`).
- New and changed code follows the files' existing idiom (base pipe, `limma::` qualified calls, roxygen `@noRd`); the
  surrounding legacy diagnostics are not reformatted.
- The trx class still does not inherit `SummarizedExperiment` (recorded divergence, `R/func_trx_s4_class.R:1-4`); unchanged.
