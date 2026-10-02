# DIA saved-entity compatibility contract

Version 1.0 — DIA-001, 2026-09-11. Status: source/ownership contract recorded;
diagnostic failures retained. This is G0 preparation, **not G0 acceptance**.

Authority: [restart specification](../reversion-document.md),
[execution plan](../../../tickets/openai/01-DIA-SAVED-ENTITIES-PLAN.md),
[ticket](../../../tickets/openai/DIA-001.md).
The machine-readable companion is [coverage.json](coverage.json): source hashes,
line-indexed ownership inventory, route/scenario/consumer ownership, exact commands,
returned diagnostic logs, failed attempts and outstanding qualification work.

## 1. Frozen source and evidence boundary

B0-original is janitor commit c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017.
Work is on performance-v2, created at that same commit. R/, DESCRIPTION, NAMESPACE,
man/, inst/ and tests/ remain unchanged against that commit. Ticket tooling,
the restart specification, plan, audit/repair records and DIA tickets are
uncommitted additions; they are not an integrated storage candidate. No package
installation, large fixture, GUI benchmark or runtime optimization was made here.

The source inventory binds the Git index entries for the package surface, SHA-256
of 194 relevant R files, test members, fixture files and the four conventions.
The line indexes cover state access across all omics and DIA payload owners,
file I/O and completion consumers. These are reproducible lexical indexes,
supplemented by the call traces below, not proof that every dynamically captured
allocation has been measured. Profiling and actual process ownership remain
DIA-003/004 requirements.

Diagnostics used R 4.6.1, pkgload/devtools source loading, and the existing library
at /home/doktersmol/R/x86_64-pc-linux-gnu-library/4.6 plus /usr/lib/R/library.
Relevant installed versions: iq 2.0.1, limpa 1.4.0, limma 3.68.4, qvalue 2.44.0,
mixOmics 6.36.0, testthat 3.3.2. The receipt contains versions, library paths,
BLAS/LAPACK and sessionInfo. MultiScholaR is not installed in those libraries.
Thus these runs cannot qualify installed startup, dispatch, completion time or
combined peak memory. DIA-003 must build and bind the control library.

PQC-001–008, the historical board, plan, package source and tests are preserved.
Their earlier pass counts are historical assertions, not this run's results.
The archived performance branch and failed evidence are not acceptance evidence
for this branch. The restart specification and original audit/repair records
remain immutable.

## 2. Corrections to inherited assumptions

These corrections freeze the code actually present; none authorizes an estimator
change.

1. **No-TREAT return:** runTestsContrasts initializes t.fit to NA, then assigns
   eb.fit on its no-TREAT branch. Its public fit.eb returns that moderated fit,
   or the TREAT fit on the other branch; it never means the uncontrasted base fit.
   The previous plan/DIA-012 statement that it stays NA was wrong. The active
   descriptions are corrected, with predecessor hashes/text in coverage.json.
   A real limma execution below verifies the no-TREAT return exactly.
2. **GUI DA is not the wrapper:** the registered ProteinQuantitativeData
   differentialAbundanceAnalysis method calls differentialAbundanceAnalysisHelper
   in R/func_prot_da_model_methods.R. That helper extracts the displayed protein
   matrix directly and calls runTestsContrasts. It does not call
   daAnalysisWrapperFunction or dpcDE and has no Limpa dispatch condition.
   The separate exported wrapper in R/func_prot_da_model_wrapper.R has the
   uncertainty-bearing Limpa branch. The loaded method-body check confirms this
   registration/call distinction; it is not a real Limpa GUI execution.
3. **Density population is call-specific:** current GUI normalization PCA and
   buildProtNormDensityPlot both call plotPca with the same arguments; the latter
   calls plotPcaBox on the returned plot. At this source, both use the helper's
   default CV population. Other direct S4 plotPcaBox/getPcaMatrix callers have
   their own populations. Preserve each effective call; do not impose a historical
   population split where this build does not have one, or merge distinct calls
   merely because both are labeled PCA/density.
4. **Fixture prose is stale:** both tracked E2E report.tsv files currently contain
   60 rows, five Protein.Group values and six runs. The Limpa README describes
   72 rows/six proteins. It is not a valid dimensions or scientific oracle.
5. **Browser sequence differs from UI order:** the current peptide tab order is
   q-value → protein evidence → rollup → intensity → sample → replicate →
   imputation. The historical browser helper runs q-value → rollup → intensity →
   protein evidence → sample → replicate → imputation, and still selects the old
   “Protein Peptides” tab label. Freeze the actual chosen order and correct the
   harness to address existing controls; never silently reorder the scientific
   pipeline to make a benchmark run.

Sources: [statistics](../../../R/func_prot_da_model_stats.R),
[registered methods](../../../R/func_prot_da_model_methods.R),
[wrapper](../../../R/func_prot_da_model_wrapper.R),
[QC generation](../../../R/mod_prot_norm_qc_generation_helpers.R),
[QC support](../../../R/mod_prot_norm_qc_support_helpers.R),
[peptide coordinator](../../../R/mod_prot_qc_peptide.R),
[browser helper](../../../tests/testthat/helper-e2e-browser.R).

## 3. Producer → consumer → completion ledger

The active GUI is synchronous at these boundaries. A returned S4, completion
flag, progress notification or existing file alone is not durable publication.

| Stage / named state | Actual producer and scientific input | Consumers / retained owners / current completion |
| --- | --- | --- |
| Import | importDIANNData and import processing helpers read the whole TSV/Parquet; retain column mapping and sanitized run-name mapping | workflow_data data_tbl and data_cln, import module reactives, design builder; setup_import complete enables design |
| Design / raw_data_s4 | buildProtDesignStateCheckpoint constructs PeptideQuantitativeDataDiann; tech_rep_group = paste(group, replicates, sep = "_") | state manager, design matrix/config, QC modules; annotation lookup/save completes or warns, then design_matrix complete and qc_trigger(TRUE) |
| qvalue_filtered | runPeptideQvalueApplyStep → srlQvalueProteotypicPeptideClean; fixed identification confidence/exclusion and evidence freeze | full peptide S4 in states; audit lineage/support, module result/grid and global FilteringProgress |
| protein_peptide_filtered | runProteinPeptideApplyStep → filterMinNumPeptidesPerProtein using frozen identification evidence | peptide S4, audit record and filtering grids; state name/observer return |
| precursor_rollup | runPeptideRollupApplyStep → rollUpPrecursorToPeptide | peptide S4 with changed long table and quantity-column selectors; raw and normalized linear sums |
| intensity_filtered | runPeptideIntensityApplyStep → peptideIntensityFiltering | filtered S4, threshold/support summaries, grids |
| sample_filtered | runPeptideSampleApplyStep → filterMinNumPeptidesPerSample | filtered data/design, retained support/exclusion evidence and plots |
| replicate_filtered | runPeptideReplicateApplyStep → removePeptidesWithOnlyOneReplicate | filtered S4/design and technical-replicate support |
| imputed | runPeptideImputationStep → peptideMissingValueImputation | Peptide.Imputed and is_imputed in long data; args support table; no GUI peptide-matrix construction appears in this path |
| protein_s4_created (IQ) | runProteinIqRollupApplyStep writes IQ input, runs actual iq::process_long_format, reads output, constructs ProteinQuantitativeData | state manager plus returned protein object/table, IQ files, plotting/audit owners |
| protein_s4_created (Limpa) | runProteinLimpaRollupApplyStep → proteinMissingValueImputationLimpa(PeptideQuantitativeData); real package if available, test-mode fallback otherwise | protein display table, complete uncertainty EList and peptide inputs in args; Limpa-specific qualification held |
| protein_accession_cleaned | runProteinAccessionCleanupStep → chooseBestProteinAccession | identifiers/annotations; Protein.Group route preserves quantitative grouping instead of regrouping accessions; legacy accession path differs |
| protein_intensity_filtered | protein intensity module → removeRowsWithMissingValuesPercent | full filtered protein state plus missingness/intensity support and plots |
| duplicates_removed | runProteinDuplicateRemovalStep performs grouped aggregation using the selected function | protein table, resulting S4 and grids; preliminary duplicate counting still names Protein.Ids literally |
| protein_replicate_filtered | replicate module → removeProteinsWithOnlyOneReplicate | checkpoint and TSV; coordinator directly observes states$protein_replicate_filtered and enables normalization |
| Pre-QC / between-sample normalization | current getState() → pre-QC; normaliseBetweenSamples → post-QC and normalized TSV | norm_data normalized_protein_obj, current_s4/run_context closures, path lists, temporary plots; normalization_complete flag is local, not output acknowledgement |
| normalized (RUV skip) | applyProtNormSkippedRuvState | same normalized S4 also in ruv_normalized_obj; skip metadata file and state save can warn and continue |
| ruv_corrected (manual/automatic) | resolveProtNormRuvParameters → ruvIII_C_Varying → removeRowsWithMissingValuesPercent | normalized and corrected inputs, controls/cancor/optimization tables, cleanup object; checkpoint save errors can be swallowed |
| correlation_filtered (Apply/Skip) | correlation helper computes/filter samples or passes RUV object unchanged; final QC, TSV/RDS precede completeProtNormCorrelationWorkflow | norm_data correlation_filtered_obj and workflow_data ruv_normalised_for_da_analysis_obj; saveState then tab normalization complete/DA pending; write failure can still lead to success |
| Export filtered session | collectProtNormExportSessionData → saveProtNormExportArtifacts | all states plus two active S4 fields and metadata; timestamped and latest each serialized separately; explicit export notification after summary; optional metadata writer swallows errors |
| Load filtered session | da_server_load_session_handler → readRDS → validateProtDaFilteredSession → public states assignment | restored full history, da_data current_s4_object, design/config/contrasts/annotation globals and workflow aliases; tab statuses/trigger updated |
| DA | run-analysis observer loops contrasts; registered S4 methods → runTestsContrasts; current state overwritten with UI-parameter metadata | per-contrast result lists each contain theObject and plots; combined long table, first-result aggregate, workflow_data full result list; analysis_complete/tab complete are set before required output writer |
| DA table/volcano/heatmap | result lists → current view subset, plots/ComplexHeatmap and clustering | table currently slices a fully resident table; volcano and heatmap callbacks execute substantial work; current_heatmap_plot/cluster vectors retained |
| Existing summary/file consumers | final getState, save args, integration S4, publication copy, template/report and summary-session exporters | existing files and formats remain dependencies; enrichment/new report features are held, not added by this ticket |

The broad state-access index also names LFQ/TMT, metabolomics and lipidomics
callers. Their default resident behavior must remain unchanged. R/func_prot_da_model_wrapper.R
also handles MetaboliteAssayData; shared file/publication and WorkflowState changes
cannot be enabled merely because a store exists.

## 4. Route-specific scientific input contract

### 4.1 Common DIA import and peptide lineage

The importer chooses Precursor.Normalised in its mapping when requested/present,
otherwise Precursor.Quantity. The DIA S4 constructor nevertheless explicitly sets
raw_quantity_column = Precursor.Quantity and norm_quantity_column =
Precursor.Normalised. Therefore the import UI mapping is not, by itself, the
identity of every downstream estimator input. Bind the selected slots and columns
at each stage. Run sanitization uses janitor::make_clean_names with uniqueness;
retain the original-to-clean map and design alignment.

Protein.Group is the active biological key when present, with Protein.Ids as the
legacy fallback. Stripped.Sequence is the peptide key; Modified.Sequence and
precursor/charge identities support identification counts. Preserve reversible
feature-key maps, multi-accession annotations and ambiguity. Never parse an
opaque matrix row label as a protein identifier.

Identification uses Q.Value, Global.Q.Value and Global.PG.Q.Value <= 0.01,
proteotypic selection and explicit decoy/contaminant handling; strict main-report
and explicit upstream_prefiltered provenance modes differ. Frozen experiment-wide
distinct peptide/peptidoform evidence survives later quantitative filtering.
Missing confidence cannot be silently invented to qualify a fixture.

Precursor rollup groups by active protein/stripped peptide/run and sums linear raw
and normalized quantities across modified/charge states. All-missing groups stay
NA; partial-missing groups sum observations. Logged inputs are rejected. It sets
raw_quantity_column = Peptide.RawQuantity and norm_quantity_column =
Peptide.Normalised. Intensity filtering uses its effective quantity selector and
biological-group support; sample and replicate filtering use distinct identities
and declared technical groups, not duplicated long-table row counts.

Ordinary imputation reads the active raw quantity column, not an inferred numeric
column. Eligible missing measurements use the mean of observed raw quantities
within the same declared technical group; this is not stochastic downshift imputation.
The default maximum missing fraction is 0.50 within distinct design runs
of a technical group (eligibility <=). Explicit exclusion groups are respected.
Without a declared technical-replicate column, quantities pass through and
is_imputed is FALSE, with an explicit skip reason. Preserve the observed values,
imputation mask, group support and args. This GUI path is not peptide-level Limpa
normalization/RUV from a workbook.

### 4.2 IQ quantification and standard GUI DA

IQ consumes Peptide.Imputed from the long table with active protein/peptide/run
selectors. Before writing peptide_values_imputed.tsv the adapter maps design
runs, in design order, to S_001 etc.; changes Q.Value to 0.0009 and PG.Q.Value to
0.009; and changes missing input intensity to zero. Preserve this existing
compatibility behavior and its audit record; it is not authorization to alter
the earlier identification population.

The exact explicit IQ call is:

- sample_id = active sample column; primary_id = active protein column;
  secondary_id = active peptide sequence column;
- intensity_col = Peptide.Imputed;
  filter_double_less = c(Q.Value = "0.01", PG.Q.Value = "0.01");
- normalization = "none"; input/output paths from project peptide/protein QC.

Bind inherited iq 2.0.1 defaults too: na_string = "0",
log2_intensity_cutoff = 0, pdf_out = "qc-plots.pdf", show_boxplot = TRUE,
pdf_width = 12, pdf_height = 8, intensity_col_sep = ";". The library preprocess
log2-transforms intensity and applies its cutoff; MaxLFQ output is log2-scale.
The QC PDF is an implicit working-directory output, not free or optional merely
because the adapter does not name it. DIA-003 must isolate/capture that file in
the run boundary without rewriting the producer. No double log transform.

IQ result rows/order come from the actual output file; aliased sample names are
restored. Surviving design samples are the intersection in original design order.
Missing required input columns, unknown design runs or missing IQ output are
errors (existing output wait bound 30 seconds). Do not promise real IQ parity
from a mock processLongFormatFn.

The GUI normalizes the protein display matrix, then applies the selected RUV and
correlation branch. Registered standard DA converts all non-ID columns with
column_to_rownames/as.matrix. The separate wrapper calls getDataMatrix(), which
uses the same protein extraction, but that is not the GUI call edge. Both reach
runTestsContrasts(weights = NA). Custom pre-existing data-frame row names are
rejected by column_to_rownames in the diagnostic library; preserve that rejection
instead of silently resetting imported scientific data.

runTestsContrasts constructs model.frame(formula, design), then model.matrix,
and selects assay columns by model row names in that order. Freeze na.action,
factor levels/contrasts, omitted samples, model columns and rank. GUI preflight
uses na.pass to validate the formula, which is not proof of the effective
na.omit fitting population. A missing/invalid design is not a cache hit.

With a replicates column, the blocking vector is paste(group, replicates, "_")
for model-selected samples; duplicated blocks invoke duplicateCorrelation and
then correlated lmFit. Otherwise use ordinary lmFit. Preserve the existing
literal group/replicates behavior, even when the S4 technical-replicate selector
has another name. Non-NA weights are not a new supported storage optimization:
the current conditional is not a general observation-weight matrix contract.

The exact order is lmFit → contrasts.fit → eBayes(trend, robust) → optional
treat(lfc). Direct runTestsContrasts defaults trend/robust FALSE and cutoff NA;
GUI/helper defaults are trend/robust TRUE and cutoff 0 (thus TREAT runs even at
zero). Each GUI contrast is currently analyzed separately. Reusing a base fit
must preserve each downstream contrast's tested population and adjustment order.

topTreat(n = Inf) / topTable(number = Inf) provide complete contrast tables.
Q-values use only finite non-NA P.Value entries, with BH fallback on the existing
edge/failure conditions; the separate BH column uses that same valid mask.
Invalid entries remain NA. Preserve sorted result identifiers, contrast names,
warnings and the return list results/fit.eb/qvalue_warnings. Table filtering and
pagination occur after statistics, never before FDR estimation.

### 4.3 Limpa: separate effective representations and a held GUI route

The real peptide-to-protein method requires an existing peptide_matrix and its
feature-key map. It log2(x + 1)-transforms when is_logged_data is FALSE, replaces
infinite/NaN with NA, and removes all-NA rows. DPC parameters come from the explicit
argument, stored DPC object/parameters, or limpa::dpc; on estimation failure the
existing fallback uses dpc.slope (default 0.8) in dpcQuant. Preserve verbose/chunk
(default 1000), EList E/genes ordering and protein/peptide identity mapping.

dpcQuant returns log2-scale E plus genes and other$standard.error /
other$n.observations. The method saves the **whole** quantified_elist, duplicate
support arrays, peptide matrix for DPC plots, peptide counts and DPC provenance
in args$limpa_dpc_quant_results. If input is not logged, the separate displayed
protein table is transformed back with 2^E - 1 and clipped below zero. The stored
EList is not replaced by that table.

| Subsequent operation | Display table/design | Stored quantified_elist / support |
| --- | --- | --- |
| chooseBestProteinAccession | Protein.Group preserves quantitative grouping; legacy accession path may rename/aggregate | no general EList synchronization in this method |
| ordinary protein intensity / duplicate / replicate filtering | filters or changes table/ID table; cleanDesignMatrix may reorder/drop design rows | no corresponding general EList rewrite |
| filterMinNumPeptidesPerProtein, Protein method | filters via Limpa peptide summary | explicitly subsets EList rows via genes$protein.id; separately stored support/summary fields still need their own alignment check |
| normaliseBetweenSamples | changes display values, nonfinite → NA, cleans design | leaves EList unchanged |
| ruvIII_C_Varying and missingness cleanup | correct/filter display population and clean design | leaves EList unchanged |
| correlation Apply | drops display sample columns, cleans design | leaves EList sample population unchanged |
| correlation Skip | display passed through | EList passed through |
| separate daAnalysisWrapperFunction, with Limpa installed | design matrix for dpcDE built from S4 design | dpcDE consumes quantified_elist directly; no last-moment alignment or display substitution |
| actual registered GUI DA helper | uses displayed protein matrix and standard limma | does not select dpcDE at this build |

No GUI calcPeptideMatrix call was found in the normal peptide sequence. The real
rollup's missing-matrix rejection is therefore a baseline risk to test, not
permission to insert a new transform. Missing Limpa at rollup is an error outside
test mode; the test-mode fallback is not real Limpa evidence. The separate wrapper
falls back to standard limma if Limpa is unavailable, even with DPC metadata.

DIA-017 must resolve actual GUI routing, matrix construction, scale, sample/feature
alignment and uncertainty semantics before a real Limpa G0/G3 pass. Any needed
scientific correction requires a separately authorized B0 versus B0-original
decision. This cannot be “fixed” incidentally by hydration, nor can a displayed
matrix key authorize an uncertainty-weighted fit. Same E with different uncertainty,
or same dimensions with different support order, is a different input.

### 4.4 Normalization, plots, RUV and view populations

Supported GUI normalization choices are cyclicloess, quantile, scale and none.
The Protein method calls normalizeCyclicLoess, normalizeQuantiles,
normalizeMedianAbsValues or no transform, respectively, after nonfinite → NA;
it restores a data frame and calls cleanDesignMatrix. Global populations remain
global. The method parameter is normalisation_method, not the config/UI spelling
method; record effective resolution as well as UI input.

Automatic RUV defaults: percentages 1:20 by 1, max_difference separation,
k penalty 0.5, adaptive penalty TRUE, grouping group. Other UI separation choices
mean_difference/auc/weighted_difference remain covered. Manual uses selected
percentage (UI 5), selected k or fallback 3, getNegCtrlProtAnova with q cutoff 0.05
and BH, then ruvCancor with num_components_to_impute = 2. Skip has NA k/controls
and explicit skip provenance, not a fitted zero-factor correction.
ruvIII_C_Varying selects samples in design order, constructs the existing replicate
matrix, corrects all selected features using named control membership, orients
the returned matrix, removes all-nonfinite rows/columns and cleans design.
Subsequent missingness/intensity cleanup remains part of the population.

PCA helper presence cutoff is floor(0.05 * sample count), followed by CV quantile
(default 0.90), then transposed input to mixOmics::pca. Preserve complete-variable
prcomp fallback and zero-score degenerate fallback. The Protein plotPca method
currently does not forward its cv_percentile formal to the helper. Density's
GUI path runs plotPca again, then plotPcaBox on its scores; direct
plotPcaBox(ProteinQuantitativeData) is a separate operation. Do not use parameter
names alone as a cache identity when the current producer ignores them.

RLE uses the full effective matrix and existing NA behavior; GUI y limits -6:6.
Pearson QC excludes pool samples. Correlation Apply uses
pearsonCorForSamplePairs(tech_rep_remove_regex = "pool",
correlation_group = selected grouping), then
filterSamplesByProteinCorrelationThreshold(min_pearson_correlation_threshold).
Default threshold 0.75; Skip retains all input samples. Apply retains samples belonging to at least one passing
pair, plus sample columns absent from the pair table; it is not an all-pairs rule.
Reapply starts from
norm_data$ruv_normalized_obj, not the last correlation-filtered result.
Style changes rerun available pre/post/RUV QC. Pre-QC resolves current getState,
not necessarily the original pre-normalization checkpoint: bind actual source
identity and treat that distinction as a baseline correctness risk.

Tables select contrast/significance first, take the leading table_max_rows, show
25-row pages and format numeric columns to four decimals. Summary counts use
the unpaged contrast population. DT CSV/Excel are view exports, not refitted DA.
The existing download_da_results ZIP handler writes placeholder text, not a
real archive; do not claim it as a complete scientific export.

Heatmap significance filtering precedes top-N absolute fold-change selection;
expression rows retain source-table order via membership, and samples use the
table/design intersection. Row/column/both/none scaling, distance/clustering,
tree cutting, ties and zero replacement after scaling remain unchanged.
No-significant-input returns NULL/placeholder. Volcano projections use the
accepted full-population result columns and existing thresholds/labels.

### 4.5 Equality rules

Storage/rehydration must be identical for values, storage types, classes, column
and row order, dimnames, factor/ordered-factor levels, attributes, annotations,
NA/NaN/+Inf/-Inf, zeros, nested args and checkpoint metadata. Codec comparisons
must not normalize away real nonfinite differences.

Existing expect_identical assertions are exact. Existing testthat edition-3
expect_equal defaults to sqrt(.Machine$double.eps), approximately
1.4901161193847656e-8 on this host; preserve the existing comparison semantics,
not a new loose percent tolerance. Explicit support-summary tolerances include
1e-8 for numeric support and 1e-5 for rounded missingness percentages. Plot
dimension tolerance 10 in reporting tests is **not** a scientific tolerance.
These assertions and mocked methods do not establish a universal IQ/Limpa/PCA/RUV
oracle. DIA-002 must bind method-specific expected arrays, seeds, justified PCA
sign/subspace handling and existing numerical comparisons before timing.

Only recorded volatile fields may be excluded in scientific replay.
.peptideQcSubstantiveAudit removes timestamp_utc, retaining record identity,
parent lineage, thresholds, support, input/output identity and removal ledger.
Do not drop args or provenance wholesale. Equality exceptions must be explicit
per field and route, with whole-object exact storage round-trips still required.

## 5. State ownership and compatible cold-history strategy

WorkflowState currently exposes an ordinary list. Entries contain timestamp,
data, config, description, previous_state, s4_class and audit_metadata (initial
has a smaller field set). saveState assigns the object; copy-on-modify does not
prove that every state is an independent full allocation. Shared allocations
and data.table/environment aliases require measurement and negative controls.

getState(NULL) reads current_state; missing names currently return NULL.
revertToState truncates state_history, sets current_state and returns data, but
does **not** delete later states or immutable audit_records. Saving an existing
name replaces its entry without duplicating that active history name. A branch
after revert retains inactive records in states. Reset chooses the latest
matching pre-normalization name from active history; it clears selected norm
objects/flags and disables DA, not the whole states list or every downstream alias.

Direct public consumers requiring real records are:

- mod_prot_qc.R's final-state observer;
- collectProtNormExportSessionData's full list read;
- DA run-analysis metadata lookup and saveState replacement;
- DA load-session whole-list assignment, with separate history/current assignment;
- utils_test_mode's state-name enumeration;
- WorkflowState methods themselves, including getStateAudit.

There is no reset() method on WorkflowState. getStateConfig is referenced by the
general study-parameter writer but is not a defined WorkflowState method; do not
mistake that guarded historical call for a supported accessor.

The implementation strategy to prove in DIA-006/008/018 is:

1. Keep the default public WorkflowState resident and its accepted read/write
   behavior unchanged. Cold operation is explicit DIA opt-in only, after gates.
2. Keep descriptors private. getState(name) hydrates a real S4 and current
   metadata is bounded. Public states access must return a normal, complete list
   of historical records with real S4 data, not a proxy disguised as that list.
3. A cold manager may implement a compatibility materialization boundary for
   states reads/writes. A full read or R replacement expression can hydrate the
   entire list. That is an explicitly measured compatibility scope, not a
   bounded lookup; avoid it in ordinary internal GUI metadata checks.
4. Whole-list assignment and nested replacement must work for valid legacy
   records; commit replacement only after validation/persistence, preserve active
   state on failure, and do not secretly fit to reconstruct history. Exercise
   list extraction, lapply, names/length, external assignment and full-list
   saveRDS/readRDS using the same examples as the resident control.
5. Legacy filtered-session export must produce the complete independent list
   envelope, including inactive states and both current S4 aliases, with no
   project/store dependency. Direct serialization of a cold R6 environment does
   not automatically invoke active-binding getters; it is not an acceptable
   substitute. Preserve default whole-manager serialization; prove any opt-in
   serialization surface explicitly before adoption.
6. Release export materializations and owned connections after use. Caller-held
   public snapshots remain caller-owned and cannot be reclaimed by the store.
   The 64 MiB optional cache cap is not a cap on required active/export payloads;
   those must fit the unchanged total service limit.

This permits normal cold operation without claiming every compatibility access
is memory bounded. It does not yet demonstrate that a cold facade passes all
R6 replacement/clone/serialization cases or the performance gates. If maintaining
those public semantics requires permanent full residency, reject that integration;
do not weaken the API. If full-history export alone exceeds 8 GiB, report the
measured conflict. Independent product experiments remain possible under the plan.

### Large-object owners that must be released or accounted for

| Owner | Retention / release requirement for the later proof |
| --- | --- |
| WorkflowState states, inactive branch records, exported public snapshots | Persist before eviction; retain checkpoint reachability, measure unique allocations |
| workflow_data data_tbl/data_cln and design builder reactive closures | Imported raw long table can outlive peptide/protein rollup; deleting a state alone will not release it |
| Module preview/apply/revert reactives and local result/grid objects | A result plot or closure can retain a full table/S4; inspect the actual closure environment |
| norm_data normalized_protein_obj / ruv_normalized_obj / correlation_filtered_obj | Three explicit scientific aliases plus step-local before/after objects; releasing only QC_composite_figure is insufficient |
| norm_data final/filtering plots, QC path lists, optimization result | Path lists are small; plots/cancor/optimization objects may embed data/environments |
| DA current_s4_object, per-contrast theObject, da_results_list and combined tables | Metadata-only @args updates may create distinct wrappers; list object.size sums overcount shared allocations |
| Current heatmap, row/column clusters, volcano plot/widget payloads | Persistent UI outputs/closures count until invalidated or session end |
| Globals filtering_progress, config_list, contrasts_tbl, uniprot_dat_cln, aa_seq_tbl_final and project context | Existing compatibility aliases are real owners; do not remove blindly or assume package-scoped state |
| Serialization/composite buffers, PNG raster loads, IQ and vroom inputs | Account for overlap within one action and release before later heavy work where semantics permit |
| GUI diagnostic/test-mode captures and full RDS oracles | Must be off or in separate processes for timed runs; they cannot retain baseline histories invisibly |

## 6. Required files and downstream lookup contract

Paths below are relative to named experiment_paths members, not guessed project
roots. Freeze both those resolved directories and all conditional obligations.
An optional unavailable panel needs an explicit reason; unexpected absence is
not an accepted optimization. Every newly written file, required copy and
compatibility serialization belongs in the later completion boundary.

| Directory / file(s) | Producer → consumer / condition |
| --- | --- |
| source_dir/aa_seq_tbl_final.RDS, fasta_metadata.RDS; FASTA cache aa_seq_tbl.RDS | Import with FASTA → accession cleanup, resumed session, report metadata |
| source_dir/uniprot_dat_cln.RDS and annotation cache | Design annotation step when available → DA reload, table/volcano labels and existing report consumers; frozen offline input must be identical between timed members |
| source_dir/design_matrix.tab, contrasts_tbl.tab (accepted legacy aliases contrast_strings.tab/contrasts.tab on read), study_parameters.txt | Design/summary → report and publication copy |
| peptide_qc_dir/peptide_values_imputed.tsv | IQ adapter's aliased/compatibility-adjusted long input → IQ |
| protein_qc_dir/iq_output_file.txt; working-directory qc-plots.pdf | Actual IQ writer/default QC → protein construction / required run inventory |
| protein_qc_dir/remove_proteins_with_only_one_rep.tsv | Protein replicate filtering output |
| time_dir/{step}_{plot}.{png,pdf} and {step}_combined_plots.{png,pdf} | updateProteinFiltering; plot names proteins_total, proteins_per_run, peptides_total, peptides_per_protein, peptides_per_run; availability depends on peptide/protein statistics and valid project path |
| protein_qc_dir/normalized_protein_matrix_pre_ruv.tsv | Between-sample normalization, before RUV |
| protein_qc_dir/{pre_norm,post_norm}_{pca,density,rle,correlation}.png | Every generated normalization QC stage; 150 dpi, PCA/density 8×6 in, RLE 10×6, correlation 10×8 → panel images and composites |
| protein_qc_dir/ruv_corrected_{pca,density,rle,correlation,cancor}.png | RUV performed, cancor when available → RUV panels/composite; cancor 8×6, 150 dpi |
| protein_qc_dir/proteomics_composite_QC_figure.{rds,pdf,png} | generateProtNormCompositeQcFigure → savePlot always writes plot RDS, then default PDF/PNG; dimensions determined by actual two/three-column layout |
| source_dir/ruv_optimization_results.RDS | Automatic/manual/skip metadata → reload, summary/export |
| protein_qc_dir/normalised_results_cln_with_replicates.{tsv,RDS} | Final correlation Apply/Skip with RUV skipped; RDS is the full final S4 |
| protein_qc_dir/ruv_normalised_results_cln_with_replicates.{tsv,RDS} | Same with RUV applied; name selection is based on stored skip metadata |
| source_dir/filtered_session_data_YYYYmmdd_HHMMSS.rds, filtered_session_data_latest.rds, filtered_session_summary.txt | Explicit filtered-session export → DA loader; complete full-history envelope, two serializations at baseline |
| source_dir/accession_cleanup_results.RDS, qc_params.RDS, mixed_species_analysis.RDS | Existing available metadata → export/report; FASTA/RUV metadata only written by export if absent, other named metadata rewritten |
| source_dir/limpa_parameters.RDS, limpa_dpc_quant_results.RDS | Limpa metadata present → downstream/report; full EList/support payload must be preserved |
| da_output_dir/da_proteins_{safe_contrast}_long_annot.{tsv,xlsx} | GUI outputDaResultsAllContrasts; safe name replaces non A–Z/a–z/0–9/_/- with underscore → DA publication/report consumers |
| publication_graphs_dir/Volcano_Plots/{safe_contrast}.{png,pdf}, all_volcano_plots_combined.pdf | Per-contrast GUI automatic writer; complete result population, no new fit |
| publication_graphs_dir/NumSigDaMolecules/da_proteins_num_sig_da_molecules.{tab,xlsx,png,pdf} | Registered writer's conditional summary products; file inventory must use the loaded method, not another duplicate definition |
| publication_graphs_dir/Heatmap/{safe_prefix}_heatmap.{png,pdf}, {safe_prefix}_heatmap_methods.md, {safe_prefix}_clusters.csv | Explicit Save Heatmap; prefix sanitizes prot_{contrast} to alphanumeric/underscore, CSV conditional on row clusters; images 10×8 in, PNG 300 dpi |
| results_summary_dir/Study_report/{design_matrix.tab,contrasts_tbl.tab,study_parameters.txt} | Existing publication copies/write-outs |
| results_summary_dir/Publication_tables/{normalised_results.tsv,normalised_results.RDS} or {RUV_normalised_results.tsv,ruv_normalised_results.RDS} | Copy from final QC files; existing RUV-skip/legacy fallback resolution |
| results_summary_dir/Publication_figures/{Volcano_Plots,Interactive_Volcano_Plots,NumSigDaMolecules,Heatmap}/… and QC_figures/… | Existing publication copy contract; only produced/selected capabilities, no new enrichment |
| results_summary_dir/Publication_tables/DA_results_proteomics.xlsx | Existing DA table aggregation consumer |
| source_dir/DIANN_report.rmd or DIANN_limpa_report.rmd, rendered report/download; source_dir/session_state_{Date}.RDS | Existing browser/report-file compatibility surfaces; summary session state is metadata, NOT the full-history filtered session |
| integration_dir/proteomics_{experiment_label}_final_s4.RDS | Existing Save Workflow Arguments final S4 export |

Filtering step labels are not numeric checkpoint identifiers: current labels
include 2_qval_filtered, 3_protein_peptide_filtered, 4_precursor_rollup,
5_intensity_filtered, 6_sample_filtered, 7_replicate_filtered, 8_imputed,
9_protein_s4_created, 10_protein_accession_cleaned, 11_protein_intensity_filtered,
12_duplicates_removed, 13_protein_replicate_filtered, 11_RUV_filtered and
12_correlation_filtered. Preserve names even where numbers overlap.

The standalone legacy outputDaAnalysisResults writer additionally has PCA/RLE,
num_of_values, plotSA_after_ruvIII, fit.eb.RDS, lfc_qval_long, wide/long annotated
tables, p-value plots and plot RDS products. Do not require that entire separate
inventory for the GUI merely because the function exists. Conversely, don't
claim GUI automatic outputs passed by testing that different writer.

Publication risks to qualify: the copy helper expects unprefixed
composite_QC_figure.{png,pdf}, while current normalization writes a prefixed name;
the main report resolves both normalized/RUV final RDS names, whereas the Limpa
report's correlation reader names only ruv_normalised_results.RDS. Bare existence
of a stub report or placeholder ZIP does not prove those consumers work.
DIA-003/004 must enumerate the actual output set after each original action;
DIA-014/020/021/022 retain complete contents and later tests close these obligations.
No new report/enrichment implementation is authorized here.

## 7. Failures, publication and recovery boundaries

The existing tests explicitly characterize several success-after-write-failure
paths. Preserving scientific/public data does not mean approving that completion
behavior for the candidate.

- Normalized TSV, skip/selected RUV metadata, RUV checkpoint save and QC figures
  can fail with only logged warnings.
- Both correlation workflows catch required save errors and return the S4;
  their observers then save correlation_filtered and enable DA.
- The filtered-session writer can partially publish timestamped/latest/summary,
  with no multi-file acknowledgement; metadata saves may be incomplete.
- DA sets analysis_complete and differential_expression complete before writing
  outputDaResultsAllContrasts; its write exception is caught.
- Heatmap image/method/cluster saves catch individual errors before the observer
  reports success. Summary can create an error-bearing study-parameters fallback.

Freeze this as the original behavior, not a target success contract. DIA-022 owns
required-output status and journaled publication across these included consumers.
A successful candidate action requires committed science plus all required outputs
acknowledged. Output failure preserves committed science and prior complete
outputs, keeps dependent actions disabled, and retry does output work only.
Failure before checkpoint commit leaves the previous state usable. Cancellation,
worker loss, stale project incarnation, partial publication and loss before
acknowledgement require fresh-process recovery tests; none was run here.

## 8. Diagnostic evidence and remaining gates

| Receipt | Actual outcome | What it establishes |
| --- | --- | --- |
| D001-primary | 709 passed; 0 failures/warnings/skips | Existing peptide module/identification assertions; no matching workflow-state test file exists |
| D001-expanded | 4 errors, 1 warning, 1 skip; command exit 0 | Wider source-level contracts were executed; exit 0 did not mean test success |
| D001-supplement | 636 passed, 3 failures, 5 errors, 2 warnings, 1 skip | Additional rollup, peptide normalization, Limpa characterization, PCA and summary contracts; not a passed suite |
| D001-public-v1 | Failed: attempted unexported constructor; sandbox timezone warning | Diagnostic authoring failure retained; no package modification |
| D001-public-v2 | Failed: matrix-derived custom data-frame row names rejected | Example corrected to use valid ordinary row names; original failure retained |
| D001-public-v3 | Passed: real S4/list read/write, revert/branch, full legacy RDS and real no-TREAT return | Bounded source-mode compatibility smoke only; sandbox timedatectl warning retained |
| D001-dispatch / defaults | Passed method-registration and library-default inspection | Actual GUI/helper vs wrapper distinction, IQ inherited options, testthat tolerance |

Expanded errors: invalid FakeProteinIqInputState lacks protein_id_column;
FakeSharedProteinRollupProteinState lacks that slot; direct protein-normalization
fixture omits required protein_id_table; Limpa converter fails qvalue smoothing.
The IQ timeout test skips based on its historical seam condition. One wrapper
test warns of a missing raw_pvalue column for its second contrast.

Additional IQ baseline risk: some later protein GUI checks still name Protein.Ids
literally (deduplication pre-count, RUV cleanup/final counts/export), whereas the
new active key can be Protein.Group. This can cause failure or incorrect counts;
no real IQ GUI run has qualified those paths. DIA-004 must resolve the measured
control outcome rather than treat the earlier small suite as a whole-route pass.

Supplemental failures all arise in the existing Limpa characterization suite:
five “Too few rows of data” errors and three failed assertions about captured/
returned protein IDs and annotations. Its expected doubles did not establish
real Limpa correctness. Deprecated ruvCancorFast warnings remain; the cp03 snapshot
is only a Git LFS pointer, so that assertion skipped. These requirements stay held,
not passed. Do not fetch or fabricate a binary implicitly.

The full returned logs and exact commands are in coverage.json. One early
expanded-output chunk was tool-truncated; that limitation is explicitly recorded
and the final failure/warning/skip summary retained. No failed run was discarded
or replaced by a later successful record. A first oversized read-only inventory
response also exceeded the tool output limit; it was replaced by compact line
indexes, not by omitting failing evidence.

Required future work:

- DIA-002: realistic valid IQ recipes first, all C01–C26 scenarios, actual expected
  quantities/order/masks, representative scale and numeric protocol freeze.
  Existing tiny fixtures are not representative workloads.
- DIA-003/004: installed original build, original browser actions, matched frozen
  topology/threads/native dependencies/filesystem, full output receipts, profiling
  separately from serial timed runs, valid IQ control.
- DIA-008/009/015/016: exact codec/public behavior, cold-history/product tradeoffs,
  complete exports and actual IQ science before GUI integration.
- DIA-017: separately qualify real Limpa inputs, routing/scale/uncertainty and
  failed characterization cases after the accepted IQ pilot.
- DIA-023/024: all included consumer/file contents, external restore, faults and
  ten successive cycles, plus both fixed performance gates.

G0, G1, G2_history, G2, G3 and G4 are **not accepted** by this document.
No elapsed saving, CPU improvement, combined-peak reduction, scaling claim or
cold-store compatibility pass is demonstrated. Existing 8 GiB / zero swap /
900-second limits and 64 MiB optional cache ceilings are unchanged.

## 9. Rollback and convention compliance

This ticket changes documentation/metadata only. Rollback means retaining these
evidence files and their hashes, reverting only this ticket's documentation or
marking its contract superseded, and restoring the two corrected ticket-plan
sentences from the recorded predecessor if required. Do not reset performance-v2,
delete diagnostics, restore archived fixtures wholesale, or rewrite PQC history.

The original resident path remains present and the bounded real-object
serialization/state smoke passes. There is no store-backed session to migrate
or hydrate before disabling anything. The same full-history envelope must remain
recoverable if a later store implementation is rolled back.

All four conventions were read. No unrelated legacy style, globals, synchronous
callback, S4 inheritance, dependency locking or namespace-export convention was
“fixed” here. Preserving established architecture/public classes takes precedence
over new-project patterns. The examples use the actual exported S4 class through
methods::new, not a new public constructor or fake class.

## 10. Executable compatibility example

Run from the repository root with Rscript --vanilla. This is a small deterministic
source-mode diagnostic, not a timed or installed-package acceptance member.
It serializes the real legacy envelope through RDS raw connections (no generated
fixture file), retains inactive branch records, restores a new manager, and
compares the no-TREAT fit to the actual limma composition. Disk round-trip,
fresh-process restore, all state classes and cold-store replacement/serialization
remain required by later tickets.

```r
pkgload::load_all(export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
local({
    set.seed(1001)
    quantities <- matrix(rnorm(40 * 8, 22, 1), nrow = 40,
        dimnames = list(paste0("P", seq_len(40)), paste0("s", seq_len(8))))
    design <- data.frame(Run = colnames(quantities),
        group = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
        replicates = rep(seq_len(4), 2), row.names = colnames(quantities))
    protein_ids <- data.frame(Protein.Ids = rownames(quantities))
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(protein_ids, quantities, check.names = FALSE, row.names = NULL),
        protein_id_column = "Protein.Ids", protein_id_table = protein_ids,
        design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", args = list(globalParameters = list(workflow_type = "DIA")))
    stopifnot(methods::validObject(object), identical(MultiScholaR:::getDataMatrix(object), quantities))
    manager <- WorkflowState$new()
    manager$setWorkflowType("DIA")
    manager$saveState("raw_data_s4", object, list(route = "iq"), "fixture")
    manager$saveState("normalized", object, list(method = "none"), "unchanged quantities")
    stopifnot(identical(manager$states[["raw_data_s4"]]$data, object),
        identical(manager$states$normalized$data, manager$getState("normalized")))
    snapshot <- manager$states
    snapshot$normalized$config$external <- TRUE
    manager$states <- snapshot
    manager$states$normalized$description <- "external nested replacement"
    stopifnot(isTRUE(manager$states$normalized$config$external),
        identical(manager$states$normalized$description, "external nested replacement"))
    manager$revertToState("raw_data_s4")
    stopifnot(identical(manager$getHistory(), c("initial", "raw_data_s4")),
        "normalized" %in% names(manager$states), identical(manager$getState("normalized"), object))
    manager$saveState("branch", object, list(branch = TRUE), "new branch")
    stopifnot(identical(manager$getHistory(), c("initial", "raw_data_s4", "branch")),
        "normalized" %in% names(manager$states), is.null(manager$getState("missing")))
    manager$saveState("branch", object, list(reapply = TRUE), "replace existing name")
    stopifnot(sum(manager$getHistory() == "branch") == 1L)
    workflow <- list(state_manager = manager, design_matrix = design,
        config_list = list(globalParameters = list(workflow_type = "DIA")),
        contrasts_tbl = data.frame(contrasts = "groupB-groupA",
            full_format = "B_vs_A=groupB-groupA", friendly_names = "B_vs_A"))
    norm <- list(correlation_filtered_obj = object, best_k = NA,
        correlation_threshold = NULL)
    envelope <- MultiScholaR:::collectProtNormExportSessionData(
        workflow, norm, list(norm_method = "none", ruv_mode = "skip"),
        timeFn = function() as.POSIXct("2026-09-11", tz = "UTC"),
        messageFn = function(...) invisible(NULL))
    connection <- rawConnection(raw(), "wb")
    saveRDS(envelope, connection)
    serialized <- rawConnectionValue(connection)
    close(connection)
    connection <- rawConnection(serialized, "rb")
    restored <- readRDS(connection)
    close(connection)
    stopifnot(identical(restored, envelope),
        identical(restored$r6_complete_states$normalized$data, object))
    new_manager <- WorkflowState$new()
    new_manager$states <- restored$r6_complete_states
    new_manager$state_history <- restored$r6_state_history
    new_manager$current_state <- restored$r6_current_state_name
    stopifnot(identical(new_manager$getState(), object),
        identical(new_manager$revertToState("raw_data_s4"), object))
    contrast <- "B_vs_A=groupB-groupA"
    actual <- runTestsContrasts(quantities, contrast, design, "~ 0 + group",
        weights = NA, treat_lfc_cutoff = NA, eBayes_trend = FALSE, eBayes_robust = FALSE)
    model <- model.matrix(~ 0 + group, model.frame(~ 0 + group, design))
    contrasts <- limma::makeContrasts(contrasts = contrast, levels = colnames(model))
    expected <- limma::eBayes(limma::contrasts.fit(
        limma::lmFit(quantities[, rownames(model)], design = model), contrasts = contrasts),
        trend = FALSE, robust = FALSE)
    stopifnot(identical(names(actual), c("results", "fit.eb", "qvalue_warnings")),
        inherits(actual$fit.eb, "MArrayLM"), identical(actual$fit.eb, expected),
        identical(names(actual$results), contrast))
    message("DIA-001 PUBLIC CONTRACT PASS: actual S4/state read/write/revert/branch/full legacy RDS and no-TREAT return.")
    packages <- c("MultiScholaR", "R6", "devtools", "pkgload", "testthat", "iq",
        "limpa", "limma", "qvalue", "mixOmics", "dplyr", "data.table", "vroom",
        "shiny", "ggplot2", "ComplexHeatmap", "jsonlite")
    installed <- installed.packages()
    metadata <- list(r_version = R.version.string, libraries = .libPaths(),
        package_load_mode = "pkgload source (NOT installed candidate)",
        versions = lapply(setNames(packages, packages), function(name) {
            indices <- which(installed[, "Package"] == name)
            if (!length(indices)) return(list(status = "missing"))
            as.list(installed[indices[[1]], c("Version", "LibPath", "Built")])
        }),
        rng = RNGkind(), options = options()[c("na.action", "contrasts")],
        session_info = capture.output(sessionInfo()), legacy_bytes = length(serialized))
    cat(jsonlite::toJSON(metadata, auto_unbox = TRUE, pretty = TRUE), "\n")
})
```

## Shared table I/O repair (DIA-025, 2026-09-25)

Small correctness, lifecycle and full-journey science/output checks pass. The
small comparison fails performance acceptance: only ~1.7% lower I/O counters and
11.64% higher combined peak. These helpers remain experimental; see the complete
receipts and remaining anonymous-memory problem in gui-capacity-results.md.

The shared helpers require a `WorkflowEntityStore`, references and explicit column
roles; they do not require a DIA manager or scientific class. Public S4 methods,
legacy payloads and native archive layout remain unchanged. Restoration v3 keeps
the v2 types/attributes but derives row count from the Parquet footer; descriptor
shape and the final writer count must agree. Existing v1/v2 files remain readable.

- `beginWorkflowTableAction(store)` / `endWorkflowTableAction(store, action)` share
  action-local validation receipts and at most one DuckDB connection. They retain
  no assay. Receipts expire at action exit; external recovery still validates fully.
- `newWorkflowTableReadScope(..., metadata_only = TRUE)` and
  `readSavedWorkflowTablePrototype(scope)` inspect descriptors/footer metadata
  without payload hashing or database construction. Actual reads validate once
  per action and check the file identity throughout. Ordinary validation never
  performs durable recovery flushing.
- `streamSavedWorkflowTable()` reads explicit projections sequentially through
  Arrow, or uses the existing fastduckdb range bridge for an explicitly selected
  fast scope. SQL remains explicit for aggregation, ranking and joins.
  `streamWorkflowTableRows()` assembles wide output rows from <=64-column requests.
  Existing 10,000-row / 64-column / 64-MiB per-read ceilings remain in force.
- `openWorkflowTableWriter(path, metadata)`, `appendWorkflowTableWriter(writer,
  table, at)` and `finishWorkflowTableWriter(writer)` encode bounded rows to one
  Parquet file. Always register `abortWorkflowTableWriter(writer)` for failure.
  Use these inside `saveWorkflowTable(..., metadata = metadata, write_fn = ...)`
  so the established descriptor, hash, sync and commit sequence remains authoritative.
  With v3, unknown initial row count is `metadata["rows"] <- list(NULL)`; finishing
  supplies the physical footer count before a descriptor can be committed.
- `openWorkflowTableGroups(scope, columns, group)` exposes one projected cursor
  of complete independence groups. Its `next_group()` returns an exact group plus
  `workflow_source_order`; `close()` is mandatory on failure. No assay copy/index
  is created. Scientific callbacks still receive their complete declared design.
- `writeWorkflowTableChanges()` merges a narrow source-row/change map and optional
  generated rows with one sequential source pass. Callers explicitly name the
  source-row and generated-row ordering columns; unchanged columns are copied
  only into the final output. Maps cannot duplicate or invent source-row IDs.

Normal store information now contains fixed-size I/O counters for payload hash
bytes/counts, recovery flushes, encoded outputs, staged streams, cursor use and
sampled maximum scratch size. These are operation counters, not physical device
I/O or complete heap accounting. The existing guarded runner supplies independent
combined peak, anonymous/file memory, physical I/O and elapsed observations.

Fastduckdb's pinned API still has no incremental writer/custom metadata parameter;
Arrow handles that role. No engine replacement, persistent integrity bypass,
global cache eviction or other-omics migration is implied. New large runs remain
held until the small evidence is reviewed.

The current audited repair also counts fetched group rows and the restored group
buffer bytes once per delivered buffer. These remain materialization counters,
not allocated/resident memory. Diagnostic phase records additionally expose scalar
Arrow pool allocation/high-water bytes and DuckDB query dispatch counts; they do
not log SQL or assay values.

Writer row groups are independent of read batches: at most 122,880 rows and a
32-MiB pending Arrow-data target, plus the incoming bounded batch. A single incoming
batch may itself exceed the pending target (up to the existing 64-MiB production
limit), so this is not a total-process memory guarantee. Flush drops references to its pending
Arrow tables; their reclamation follows Arrow/R ownership. Local Arrow scans disable Parquet prebuffering; cancellation
destroys the private reader without invoking its draining Close operation.

`workflowTableChangesAreIdentity()` permits component reuse only for a complete,
unique source-row map, no changed/generated columns, and identical logical
restoration metadata. Format version and hidden physical column names do not
force an unchanged assay rewrite. `snapshotWorkflowFile()` uses independent GNU
copy-on-write snapshots where supported, with normal copying as the filesystem
fallback. Its existing consumers still hash and durably publish every required
artifact. Neither helper weakens external integrity validation.

The destination-check repair and complete small journey pass independent science
and required-output verification. Observed read/write counters fall against the
prior repair, but the small elapsed/peak outcome does not meet the performance
objective; see DIA-025. Live Arrow allocation is zero at the composite boundaries
while RSS remains elevated. The full-journey SQL counter missed late namespace
loading and must be treated as unavailable, not zero. No 96/300 capacity or
production-default claim follows from these results.


The continued DIA-025 diagnosis identifies an owner outside the named GUI holders:
ggplot2's global last-plot store can retain the raster-backed normalization composite
after its local producer returns. The retained-image replay releases 542.413 MiB
live R when that root alone is cleared. Freed glibc pages also remain resident after
GC; diagnostic-only release demonstrates this in both the real GUI and plot replay.
These are diagnoses, not production fixes. Preserve prior last-plot semantics in
any future scoped repair, and do not treat an empty assay holder or zero Arrow pool
allocation as complete process-memory accounting. See the current evidence record.


The scoped composite-export repair now restores the caller's previous ggplot2
last-plot reference on exit, including failures. This preserves prior-plot semantics
and prevents the export's raster composite from surviving through ggplot2's global
store. The retained production-image replay measures zero last-plot retention and
542.561 MiB lower post-GC live R, with identical PNG/RDS output and PDF changes
limited to creation-time metadata. No production allocator trim or GC was added;
allocator residency and full-workflow capacity remain separate acceptance decisions.


Allocator audit b8UwdY locates the repaired composite's remaining reclaimable
residency primarily in free interior pages of the main glibc heap. The secondary
arena is only 1.078 MiB, while diagnostic page release removes 317.285 MiB anonymous
RSS. Fastduckdb metadata-cache clearing and DuckDB lease closure already occur;
qs2 cache-file deletion is not an allocator-release operation. These are diagnostic
findings, not a production trimming policy or full-workflow capacity proof. Arena
virtual/free counters, R live bytes and resident mapping bytes overlap and must not
be added into a fabricated total.


Stock QC products now reconstruct their display through `rebuildWorkflowDiaQcPlot`
after extracting native calculated data and bounded metadata. The reconstruction
frame receives no original S4, scientific matrix or fit. The factory restores the
caller's prior last plot on exit. Installed checks prove collectability of the
scientific input for PCA/RLE/density/correlation while keeping the returned product
usable, with native layer-data and exact tick-label parity. Explicitly injected
custom producers preserve their callback semantics and are not covered by this
stock-producer lifetime guarantee. Releasing these live roots does not itself
return all free allocator pages to the operating system.
