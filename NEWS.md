# MultiScholaR NEWS

## MultiScholaR (unreleased; becomes 0.9.1, the first alpha for testing on real data)

### New features

 * Reference mapping for individually sequenced strains (KIT-167; the application's panels and enrichment choice are behind the option
   `MultiScholaR.reference_mapping`, off by default until the tool's validation passes; the two exported functions can
   be called but are not yet validated): `mapStrainToReference()` maps a
   strain's proteins to their orthologues in one well-annotated reference of the same species by reciprocal best hits
   of DIAMOND v2.2.8 searches in both directions (repeat masking off, full extension, 80 % coverage of both sequences,
   e-value 1e-5, no identity floor by default; every status reported and counted; a strain protein takes a reference
   accession only when one-to-one). Four UniProt reference proteomes (release 2026_04, CC BY 4.0) ship with the
   package: E. coli K-12 MG1655, S. aureus NCTC 8325, K. pneumoniae MGH 78578 and S. pyogenes SF370 (STRING
   unavailable for SF370). `matchStrainGenesToProteins()` matches the strain's genes (locus tags) to its proteins from
   the GFF3 or GTF annotation without sequence search. Proteomics enrichment can use the orthologues' accessions, with
   the reference's own GO table, unmapped groups left out of the query lists and the background, and the counts of
   changed and unchanged features by mapped and unmapped; the integration features gain the mapping columns.

 * Integration enrichment (KEGG and Reactome, for metabolite and lipid data sets): besides the correlation-adjusted
   test (the default, unchanged) and the pre-ranked test, two more tests can be chosen: the correlation-adjusted test
   with negative correlations counted as zero, and a permutation test of the factor's values (PCGSE). The tab and every
   record show how each test behaved in a check on simulated data of the same kind (false calls and power; no lipid
   data were simulated), so the choice can be made on those figures. With a permutation result, the tab and the
   record's Methods state the smallest attainable p-value as an exact fraction (for example 1/1,001) and the smallest
   number of terms that can be called together at a 5 % false discovery rate, or that none can be called when that
   p-value is above 0.05 (VAL-115). Results of the default test are unchanged. The
   pre-ranked test now draws its random numbers with R's default generator whatever the session's setting; results in
   a session with default settings are unchanged (KIT-191).

 * Integration inputs record the unwanted factors a correction for unwanted variation (RUV) removed: where an omic
   applied RUV to a data set, saving its Summary writes one more file, `removed_<data set>.parquet`, beside the
   normalised and the corrected matrix (a basis of the difference between the two over the features measured in every
   run, checked against the correction's number of factors; contract v2, an addition that earlier versions read as
   before). An input saved earlier must be saved again at the omic's Summary step before the correlation-adjusted
   enrichment test can run with the adjustment for removed RUV factors (KIT-184).

 * Integration enrichment, an optional adjustment for removed RUV factors: the setting "Adjust for removed RUV factors"
   (`adjust_ruv_factors`), in the Enrichment tab and in `integrationEnrichFactor(overrides = )`, off unless chosen.
   Switched on, the correlation-adjusted tests (the default test and the form with negative correlations counted as
   zero) put the unwanted factors the RUV correction removed from the data set into their model (intercept, factor,
   removed factors), so a feature counts for what it shares with the factor beyond those factors; for a model trained
   on RUV-corrected values and features measured in every run this is identical to fitting that model to the corrected
   values. In a check on simulated data (12 biological samples, 2 removed factors) the adjusted test called sets with
   no link to what was tested, of the kinds the check decided on, in at most 5.4 % of cases (95 % interval of that
   largest figure 4.5 to 6.4 %; sets made only of the most variable features: up to 9.4 %, against 6.8 % without the
   adjustment; sets in which half the features truly followed the main factor, tested against an unrelated covariate:
   up to 8.6 %, against up to 7.0 %) and found more true terms than the test without it. In that check it did not
   protect sets of features that share unwanted variation when that variation follows the groups: such sets were called
   on the model's main factor in 8.1 to 29.6 % of cases with the adjustment, against 2.0 to 11.2 % without it and 1.6
   to 2.8 % had the true unwanted factors been known. No
   result changes unless the setting is switched on: every model, one trained on RUV-corrected values included, is
   tested as before on the values before RUV without the factors. Every result of a correlation-adjusted test now
   carries the line "Adjusted for removed RUV factors: yes (n factor(s))." or "... no." in the tab, in
   `integrationEnrichFactor()`'s summaries and in records. The adjusted test is refused when fewer than max(4, the
   smallest number of samples per feature) biological samples are left after one per removed factor, when the factor is
   a combination of the removed factors, and for a data set without usable removed factors. The pre-ranked test, the
   permutation test and STRING do not use the removed factors. A recorded result replays in the form it was saved in; a
   result recorded before this version replays without the adjustment (VAL-53, VAL-58, VAL-62, VAL-160, VAL-178;
   KIT-184).

### Fixes

 * Every omic: when the table changed after the differential analysis (a Reset and a new run, a Revert, the
   correlation filter applied again, or Quality Control reverted) and the analysis was not run again, the Summary's
   'Save Workflow Arguments', 'Generate Report' and 'Copy to Publication Directory' now refuse with a message; a
   refused Save also removes this omic's earlier integration input. 'Generate Report' and 'Copy' also refuse when the
   results folder holds differential analysis files that no analysis of this session made (from an earlier session,
   an opened session file, or a contrast the last analysis did not produce), and name them. The differential-analysis
   tab of lipidomics, metabolomics and transcriptomics now clears its results when the table changes, as proteomics
   already did, and after a session file is loaded and the table then moves, the next analysis runs on the current
   table. Before, the saved parameters, the report and the publication folder could combine one table's parameters
   with another table's differential analysis. They still can if the analysis is run again and Copy or the Report is
   pressed without pressing 'Save Workflow Arguments' first; the Summary does not check this. 'Generate Report' after
   the analysis was run again since the last 'Copy' now asks for 'Copy' first. A saved analysis record whose
   differential analysis was run on another table than its final one no longer replays. To check an earlier save: in
   analysis_record.json, the differential_abundance entry's source_state should be the last step's state, or, for
   proteomics with enrichment, the last state before enrichment_completed; if it is empty or names another state, run
   the differential analysis again and save again (KIT-212).
 * Proteomics enrichment: a list of up- or down-regulated proteins that could not be tested (on g:Profiler none of its
   proteins has a gene name; on the Gene Ontology route none has a Gene Ontology annotation that the analysis tests) is
   now shown in the tab, the download and the report as 'Not tested', with the reason, and is recorded in the analysis
   record; a direction in which no protein passed the cutoffs now says 'no list'. Before, both looked like a list in
   which no term was enriched. A list that was tested and in which no tested term holds one of its proteins is shown,
   as before, as tested with no term. On the Gene Ontology route a comparison that has a list but no tested protein
   with an annotation, or no Gene Ontology term of 10 to 500 tested proteins, is now refused with a message, as is a
   run in which the test of a list stops with an error (on both routes); records of such runs made before this change
   no longer replay. On g:Profiler, a list with no significant term cannot be told apart from a list in which the
   service recognised none of the gene names; the report says so (KIT-215).

 * The help for the RUV grouping variable in proteomics advised selecting 'batch' for TMT. That advice was wrong: with
   'batch' the correction is built from the differences between your groups and the differences between batches stay.
   To check an earlier analysis, look for 'RUV grouping variable' in its study_parameters.txt (or the normalisation step
   of its analysis record); if it is not 'group', run the normalisation again with 'group' and, if you know the batches,
   add batch to the model formula. The normalisation tab now warns when the chosen column puts samples of different
   groups together (KIT-174).
 * Proteomics: 'Reset to Pre-Normalization' now returns to the table the normalisation started from. Before, it
   returned to the newest of a fixed list of step names. When the normalisation had started before the Protein
   Replicate Filter, that was, on label-free (MaxQuant, FragPipe, Spectronaut) and TMT data, the design table (so the
   next normalisation ran without the quality-control filters), and on DIA with iq a peptide table (undoing the protein
   roll-up and protein Quality Control; the next normalisation then failed). On DIA with limpa every Reset returned the
   precursor roll-up (undoing Log2 Matrix, Peptide Accession Cleanup, DPC Quantification, Between-Sample Normalisation
   and the Peptide Count Filter). Reset pressed before any normalisation also undid Quality Control steps; it now leaves
   Quality Control as it is. When the starting table can no longer be found, Reset now changes nothing and says so
   (KIT-225).

 * Proteomics: the enrichment tab's 'Download All Results' now writes its archive itself. Before, it needed a separate
   `zip` program, which the installed Windows application does not have (and which some Linux computers lack), so on
   those computers the download produced no archive. This change does not alter the files in the archive (KIT-241).

 * Proteomics: 'Export Filtered Data Session', and the session archive that 'Save Workflow Arguments' refreshes in the
   Summary tab, now refuse with a message when the current table is not the normalisation tab's last result, for
   example after a quality-control step was reverted or applied again following normalisation. Before, the session
   file named that table as current while carrying the earlier normalisation's results. 'Export S4 compatibility RDS',
   which runs only in a guarded diagnostic session, refuses in the same case. The session export now takes its results
   file and its statement of whether RUV was applied from the normalisation that was run, not from the current RUV
   setting; a missing file is reported plainly and no longer blocks the next steps until the application is restarted.
   A Summary save refused because the current table is not the normalisation tab's last result also removes the
   proteomics input that the Integration tab reads from an earlier save, as in the other omics. Exports made at the
   end of the normalisation tab, and projects reopened from them, are unchanged, except that the session file now
   states the RUV mode of the run: on DIA with limpa the RUV covariates step's mode, and 'none' where the chosen rule
   applied no RUV, rather than the RUV setting shown. In a session file exported before this version, the 'RUV Mode'
   line of filtered_session_summary.txt is the setting shown at export. The RUV the run applied is in the project's
   analysis_record.json (the normalisation step's ruv_applied) and in the normalisation tab's RUV summary. No result,
   report or Methods text used that line (KIT-206).

 * Lipidomics and metabolomics: the differential analysis no longer stops with 'prior.weights contain NA values' (or
   'covariate contains NA or infinite values') when an assay holds a feature with no value in any sample. Examples are a
   feature below the detection limit everywhere, which the intensity filter keeps, or a feature measured only in samples
   the correlation filter removed. Such a feature is not tested: it appears in each comparison's table without a fold
   change or p-value and is counted in the summary's 'Not tested' line, and it no longer affects any other feature's
   result. Before, the assay's analysis failed, and the whole analysis stopped when every assay failed. A saved analysis
   that failed this way now completes when it is replayed. With the default settings every result produced before is
   unchanged, except in an assay of only one or two features, one of them without values, where the B statistic changes;
   an assay that failed this way while the other assays completed now has results too, also when the saved analysis is
   replayed. The older S4 method `differentialAbundanceAnalysis()` on lipidomics and metabolomics objects runs the proteomics test and is repaired with it (see the proteomics entry below, KIT-200). With the trend option off (available only when calling the functions directly), on data that contain such
   features the result tables change: those features move to the end of each table without an average intensity. With
   the robust option also on, the other features' results change, to those of the same data without them. With RUV
   factors in the model these features were already left out of the tables and listed separately; that is unchanged. An
   assay in which no feature has any value is now refused with 'no feature has a value in any sample, so nothing can be
   tested'. `runLipidsDA()` and `runMetabolitesDA()` return `features_without_values` (per assay, the ids) when there
   are such features (KIT-196).
 * Proteomics: the differential analysis no longer stops with 'prior.weights contain NA values' (or 'covariate contains NA
   or infinite values') when a protein has no value in any of the analysed samples, for example a protein measured only
   in samples the correlation filter removed. As for lipidomics and metabolomics, such a protein is not tested: it
   appears in each comparison's table without a fold change or p-value and does not change any other protein's
   differential-abundance result. Like every row without a p-value, it is left out of that comparison's enrichment
   background (KIT-205). This applies to every
   proteomics analysis route except DIA with limpa, which was not affected, and to the older
   `differentialAbundanceAnalysis()` on lipidomics and metabolomics objects. Results of analyses that completed before
   with the default settings are unchanged, except for data of only one or two proteins, one of them without values,
   where the B statistic changes. When the trend option is off (possible only in scripted calls) on data that contain
   such proteins, their rows move to the end of each table without an average intensity and are no longer in `fit.eb`; with the robust option also on and unequal residual degrees of freedom, the other proteins' results can change slightly, to those of the same data without them. Data in which no protein has any value is
   refused with 'no feature has a value in any sample, so nothing can be tested'. `runTestsContrasts()` returns
   `features_without_values` (the ids) when there are such proteins (KIT-200).
 * Lipidomics and metabolomics: after Quality Control was reverted past the normalisation, or after 'Reset to
   Pre-Normalization', 'Apply Correlation Filter', 'Skip Filtering & Proceed to DE', 'Export Normalized Data' and the
   differential analysis could run on data that had not been normalised (no internal-standard step, no log2, no
   normalisation), and reported the analysis complete. They now stop with a message asking to run the normalisation
   again. Lipidomics, metabolomics and transcriptomics: pressing 'Run Normalization Pipeline' a second time without
   Reset normalised the already normalised table (log2 applied twice), and in transcriptomics a second run with fewer
   factors stopped with an error and kept the first run's factors (seen with 2 then 1), and the expression filter
   applied after normalisation kept factors estimated on the unfiltered genes; it now runs again from the table the
   first run started from, and it is refused on a table loaded from a session file and while a step is still running.
   Quality-control steps can no longer be applied to a normalised table; the message names the step to Revert to. A
   table that carries the name of a normalisation step is treated as normalised even when it was made by a script, or
   loaded from a session file, rather than by the tab: quality control and the normalisation refuse it, and the message
   asks to Revert to the last quality-control step (or, when none is left, to import the data again). Scripts that call
   the package's step functions directly are not guarded (in proteomics the exported quality-control step functions
   refuse a normalised table, the internal computations do not); a script decides which table it passes. The
   differential analysis no longer runs on a session file named after a design, quality-control, reset,
   internal-standard or log2 state (such a file was exported after Quality Control was reverted past the normalisation).
   Results change only for sessions that took these routes; analysis records are not changed, and replaying one repeats
   the steps it recorded. To check an earlier lipidomics or metabolomics result, open its analysis record
   (analysis_record.xlsx, exported with the results). In the sheet 'analyses', the column source_state names the table
   the differential analysis tested. If that name ends in _intensity_filtered, _duplicates_resolved, _qc_complete,
   _reset, _itsd_norm or _log2, the analysis tested values that were not normalised, or not fully. If it ends in
   _correlation_filtered or _norm_complete, look at the sheet 'steps': unless a log2 step (state ending _log2) and a
   normalisation step (state ending _normalized) stand between the quality-control steps and that state, the analysis
   tested values that were not normalised. If the sheet 'steps' starts from a session file loaded with 'Load Filtered
   Session' (one step, operation 'legacy'), the record cannot show how that file was made: this version refuses files
   named after a design, quality-control, reset, internal-standard or log2 state, and a file named after the correlation
   step is affected only if Quality Control was reverted before it was exported. A second normalisation run without
   Reset, and quality control applied after normalisation, do not show in the sheet 'steps'. They show in
   analysis_record.json: the first internal-standard or log2 step's input (input_entity_ids), or a quality-control
   step's input, is the result (result_entity_id) of a normalisation, RUV or correlation step. The data were then
   normalised twice, or filtered after normalisation. Run such analyses again with this version (KIT-194).

 * Lipidomics, metabolomics and transcriptomics: the Summary tab's 'Save Workflow Arguments', 'Generate Report' and
   'Export Session State' now refuse with a message when the current table is not the result of the last
   normalisation, for example after Quality Control was reverted or the normalisation was reset after it ran. Before,
   the saved final object (<omic>_<label>_final_s4.RDS), the study parameters and the analysis record could describe
   the quality-control table (not normalised), an empty object, or an earlier normalisation, while the differential
   analysis shown came from the normalised data. A refused save writes nothing and removes the omic's integration
   input folder (integration/<omic>/), as a save without a final table already did, so the Integration tab does not
   read an earlier input after a refused save. This includes an input saved earlier in the same session: after a
   Revert and a refused save, run the normalisation again and save before integrating. 'Generate Report' and 'Export
   Session State' refuse without removing anything. After a Reset, press Revert once in a Quality Control sub-tab (or
   run the normalisation again) before saving; after one Revert to the RUV-corrected table that table is saved; a save
   right after a run, before the correlation step, is allowed; the Summary does not check that the differential
   analysis shown was run on the table it saves. Transcriptomics no longer writes an earlier normalisation's RUVIII
   matrix as the integration file: after a second normalisation without a new RUVIII export the save says so. In every
   omic, proteomics included, a save whose integration input is refused for another reason now names the reason. Saves
   at the end of the normalisation, after its correlation step or Skip, and after the differential analysis are
   unchanged. To check an earlier save: open analysis_record.json in the folder of study_parameters.txt (the same save
   wrote both). If the last step it lists is not a normalisation result (a state ending in _normalized,
   _ruv_corrected, _correlation_filtered or _norm_complete, or loaded_for_de; for transcriptomics trx_ruvs), for
   example a quality-control, Reset, internal-standard or log2 state, that save's final object, study parameters and
   (transcriptomics) integration file do not describe the normalised analysis: save again with this version. A
   transcriptomics save that showed 'Study parameters saved with warnings' after a second normalisation may have
   written the earlier RUVIII export as transcriptomics_normalized_data.tsv: run the RUVIII export and save again
   (KIT-207).

 * Transcriptomics: the differential analysis now runs after a normalisation without removal of unwanted variation (the
   RUV option 'No removal of unwanted variation', or an RUV design whose rule finds no unwanted factor, k = 0). The
   model is then ~ 0 + group, tested with the edgeR quasi-likelihood F-test. Before, the Differential Analysis tab
   refused these runs with 'Run the normalisation first', although the normalisation had run. Runs that applied RUV are
   unchanged. Both cases are offered only when the RUV design panel is switched on. For every user, the message shown
   when the Differential Analysis tab is pressed before the normalisation now reads 'Run the normalisation first: the
   differential analysis runs on the normalised data.', and the tab's model help says that the RUV factors enter the
   model only when the normalisation added them. The report then states that no correction for unwanted variation was
   applied.

 * Transcriptomics: on runs of the RUV design panel (the option 'No removal of unwanted variation', an RUV design whose
   rule finds no unwanted factor, or any RUV design that adds factors), the Differential Analysis tab no longer reports
   'Could not save results to disk: non-numeric argument to mathematical function' after writing the result tables and
   figures. The RUV grid summary (`ruv_optimization_summary.txt`) is written only for runs that used the RUV grid
   search; the other runs have no grid to summarise, and their RUV choices are in the analysis record and the report.
   Files written by an earlier run of the same project that used the RUV grid search (this summary, the
   canonical-correlation panels and the grid workbook) are not removed, and describe that earlier run. Runs that used
   the RUV grid are unchanged. The RUV design panel is off by default.
 * Proteomics normalisation run again after changes in Quality Control now uses the current quality-control result.
   Before, once normalisation had run, every later run reused the protein table the first run had read until Reset was
   pressed, and in some cases even the first run used a table from before the latest quality-control changes, so
   filters changed or completed afterwards were silently ignored (DIA, label-free and TMT; DIA with limpa was not
   affected). Results of such a run change: normalisation, RUV, the correlation filter and the differential analysis are
   now computed from the filtered table the Quality Control tab shows, the pre-normalisation figures and the
   pre-normalisation column of the composite QC figure are redrawn from that same table, and a notice says that the
   table differs from the previous run's. A run with no quality-control change in between gives the same results as
   before. The protein quality-control steps can no longer be applied to a normalised table; press Revert in the Quality
   Control tab to change quality control after normalisation. If a project saved earlier already holds quality-control
   steps applied to a normalised table, normalisation now stops with a message naming the step to go back to, instead
   of normalising an earlier table or normalising the data twice. "Apply Correlation Filter" and "Skip Filtering &
   Proceed to DE" now stop with a message when quality control has changed since normalisation, instead of filtering
   the earlier normalisation's result. Saved projects and analysis records are not changed; in a project saved before this release in which
   normalisation was run again after a quality-control step was reverted or applied again, that run used the earlier
   table, and its analysis record names the quality-control result as the step's input; to recompute such a project,
   press Revert in the Quality Control tab back to the step you changed, apply the steps again and press "Run
   Normalization & RUV"; replaying a record repeats the steps it recorded (KIT-182).

 * Proteomics: when protein annotation from UniProt is empty or incomplete the app now says why, keeps the notice
   until it is closed, writes the fact into the analysis record, and no longer keeps an empty or incomplete answer
   for the next attempt. An incomplete annotation stays marked as incomplete when the project is saved and reopened.
   Projects whose annotation was empty under an earlier version: save the design again. Enrichment now refuses to run
   when the project has no protein annotation (previously: the UniProt GO route stopped with a 'column not found'
   error, and the g:Profiler route was sent no protein, since data imported in the app carry no gene names of their
   own). An earlier analysis record whose enrichment ran without any protein annotation is no longer replayed: the
   replay is refused, because that run had nothing to test. Each request to UniProt is now tried up to three times,
   and asks for each protein by its accession only (`accession:` terms; an isoform is asked for by its canonical
   accession, as the design step already did), so an answer holds exactly the requested entries; a protein UniProt
   has merged into another entry (a secondary accession) still gets no annotation. DA tables written before the
   annotation was retrieved keep their blank gene names until the DA is run again. (KIT-175)
 * Proteomics: when the protein table built in a script already carries gene names and the protein annotation is
   also present, the differential analysis tables lost every gene name (the two name columns collided); they now take
   the annotation's name and fall back to the table's own. In such tables a column that both carry is now named
   `<name>.input` (the table's) and `<name>` (the annotation's), where it was `<name>.x` and `<name>.y`. On the
   g:Profiler enrichment route, a run in which a comparison has proteins in its up or down list but none of its
   proteins has a gene name is now refused with a message instead of being saved as a run with no terms. A direction
   (up or down) whose proteins all lack gene names is no longer sent empty: it is reported as not tested, by a warning
   when the run finishes and in the analysis record (`enrichment_lists_not_tested`); the enrichment tab, the download
   and the report show such a list as 'Not tested' (KIT-215). The analysis record also states,
   for each comparison, how many of its proteins have no gene name and were therefore left out of the lists and the
   background (`enrichment_unnamed_proteins`). An analysis record of a run that is now refused, made before this
   change, no longer replays; it is refused with the same message. The first change does not affect data imported in
   the application, whose names come from the annotation; the others change only g:Profiler runs in which some
   protein has no gene name (KIT-208).

### Scientific changes

 * Lipidomics and metabolomics: the intensity filter now counts only the samples of the design. Before, it took every
   numeric column of the table as a sample, both when it set the intensity threshold at the chosen percentile and when
   it counted in how many samples a feature was below it. An MS-DIAL table imported with its default identifier column
   ('Peak ID' or 'Alignment ID', whichever the export has), whose identifiers are whole numbers, had that column counted
   as a sample, so the threshold could be set by identifier numbers and whether a feature was removed could depend on
   its row number. For example, on an MS-DIAL metabolite table of 8 features and 4 samples at the 10th percentile and
   0.5, the threshold was 4.9 instead of 16,410 and no feature was removed instead of 2. A custom-format table with a
   numeric identifier or annotation column was affected in the same way, and so were the metabolomics workbooks, which
   keep retention time and m/z columns, and designs imported from a folder whose tables hold MS-DIAL's numeric columns
   (retention time, m/z, scores). Tables with no such column (for example LipidSearch imports, and custom imports whose
   identifier and annotation columns are text) give the same results as before. Analyses saved before this change
   replay as they were made. A design sample that is missing from a table or is not numeric is now refused with a
   message. The helpers `lipidIntensityFilteringHelper()` and `metaboliteIntensityFilteringHelper()` take a new optional
   argument `sample_columns`; without it they count every numeric column except the identifier column and warn which
   columns they counted. (KIT-202, VAL-118)

 * Proteomics enrichment: the background (universe) of each comparison is now the proteins tested in that comparison,
   those with a p-value in its differential-abundance table, on both the g:Profiler and the clusterProfiler routes.
   Before, it also held proteins that could not be tested (no value in a compared group, or none in the analysed
   samples), which can make a term look more or less enriched than it is (in a worked example, a term of 50 proteins
   none of which was untested, 300 of 3,000 proteins untested and lists of 150 drawn with no real change, a test at the
   5 % level called the term on 5.6 % of lists). Comparisons in which every protein was tested give the same results
   as before. Analyses recorded before this change replay with the earlier background, and their reports describe it
   (KIT-205).

 * Gene Ontology enrichment for organisms outside g:Profiler's set now calls a term significant by its
   Benjamini-Hochberg adjusted p-value over every term of 10 to 500 tested proteins, including terms with none of the
   listed proteins. The earlier flag, Storey's q-value over the terms with at least one listed protein, is still
   offered; in simulations at a cutoff of 0.05 it called terms on up to 42.9 % of random protein lists without parent
   terms and up to 22.9 % with them. Analyses recorded before this change replay with the q-value, unchanged (KIT-203).

 * Proteomics, MaxQuant, FragPipe, Spectronaut and Proteome Discoverer TMT: protein quantities are now
   log2-transformed without first adding 1. Adding 1 changed each value v by log2(1 + 1/v): 0.00014 at 10,000, 0.0014
   at 1,000, 0.14 at 10; a fold change of 10 against 80 was reported as 2.88 instead of 3. Measured on simulated data
   (2,000 proteins, two groups of five, 10 % of proteins changed by a factor of 2; 20 simulated data sets per setting):
   where typical values are near 10, fold changes were shrunk to 0.81 of their true size with the offset and 0.94
   without it; near 40, 0.96 and 1.00; from about 150 on the two agree. Low-abundance calls stay unreliable either
   way: in the lowest third of abundance, with values near 10, 14.5 % of the true changes were found with the offset
   and 8.4 % without it with RUV Automatic (5.2 % and 3.1 % with RUV skipped), and 65 % against 71 % of the calls
   there were false with RUV Automatic (50 % against 42 % with RUV skipped; the whole list's false share was 0.19
   against 0.17, and 0.05 against 0.04). The Benjamini-Hochberg procedure controls the false discovery rate of the
   whole list, not of an abundance stratum. The test used by default does not model the larger variance of log values
   at low abundance. On the package's five test data sets (values of 500 and above) no call changed (fold changes
   moved by at most 0.0001). Some Proteome Discoverer exports report signal-to-noise values, which are small; there
   the change is largest. A value of 0 is read as not measured (missing) and counted; for Spectronaut tables and for
   designs re-imported from files saved before 7 November 2026 this changes which values are tested. A negative or
   infinite value now stops the design step with a message naming the run. Saved projects keep their results; saving
   the design again applies the new transform, and a project whose design was saved by an earlier version says so
   once when it is opened. Sessions from versions before the session owner (filtered sessions adopted for DA) are not
   told. Scale normalisation (limma's two-colour M-value scaling) multiplies each sample's log2 values by a factor, so
   its result depends on the unit of the table on any data; on log2 values below 0 (values below 1) it also treats -2
   like +2, and a warning now says so when the data hold such values. Cyclic Loess and Quantile do not depend on the
   unit.

 * Proteomics (DIA, label-free, TMT): the protein missing-value filter ('Protein Intensity Filter', and the clean-up
   after RUV) now decides by counting missing values. Before, it compared each group's share of missing values with a
   percentage rounded to three decimals. With three samples per group and 'Min Replicates per Group' 2, one missing
   value of three (33.3333 %) exceeded the rounded 33.333 %, so at the tab's defaults only proteins measured in every
   sample were kept. The same happened for some other group sizes (for example 6, 7 or 11 samples with 1 required; 11
   to 14 with 2) and for three groups with 'Min Groups Required' 2. Group sizes whose allowed share is exact or rounds
   up (for example 2, 4, 5, 6, 7, 8, 9 or 10 samples with 2 required) give the same results as before. In simulated
   data (KIT-185's generator with no unwanted variation and no RUV, 2 groups of 3, 3,000 proteins, R = 100 draws per
   condition) the filter now keeps, at 5 %, 10 % and 25 % value-dependent missingness, 2,838, 2,668 and 2,086 proteins
   instead of 2,442, 2,060 and 1,261; the differential analysis finds 53, 81 and 110 more proteins (the share of truly
   changed proteins found rises from 0.65 to 0.73, 0.55 to 0.68 and 0.33 to 0.51), with the false discovery proportion
   unchanged at about 0.05 for the whole list (among the newly kept proteins about 6 to 8 % of the proteins called
   were false; a follow-up simulation on new draws, pooled over 2, 4 and 5 groups of 3 and 3 groups of 4, traced this
   higher share to lower power among them (truly changed proteins found 0.61 against 0.80 among the proteins kept under
   both rules), while unchanged proteins among them reached p below 0.05 about as often (5.1 % against 5.0 %, at most
   3.4 % more often), except the few that hold a value at or below the intensity threshold, which the filter counts as
   missing but the test uses (2.5 times as often; about 2 % of them, 15 % with four groups of three and 5 % of values
   missing); it also appears, at about half the size, when the same share of values is missing at random (KIT-197));
   these are simulated figures, and the size on real data is not yet measured. Designs with unequal group
   sizes keep the rule of the largest group, as before, so a smaller group can pass with fewer values than 'Min
   Replicates per Group' (for example one of three when the largest group has six). Analysis records made before this
   version replay with the earlier rule; each step of a record replays with the rule it was made with. In an analysis
   record, the protein intensity filter and normalisation steps made with this version list the fixed parameter
   'missing_value_comparison = counts'; a step without it was made with the earlier rule. A session saved before this
   version and continued after it applies the new rule from its next missing-value step. The peptide intensity
   filter's percentage-only call (deprecated) reads its percentages the same way. The proteomics reports now say that
   a protein failed a group when its share of missing values exceeded the allowed share, and was kept when the share
   of groups it failed did not exceed the allowed share; they said 'or more' and 'at most' of the rounded percentages.
   (KIT-187, VAL-71)
 * Integration enrichment of a factor against KEGG and Reactome pathways now uses a correlation-adjusted test by
   default (limma's camera on each feature's regression on the factor, allowing for the correlation measured between
   a pathway's features, on the features measured in every biological sample; of several features on one compound
   identifier, the one with the largest mean value is kept). The pre-ranked test used before (GSEA on each feature's
   correlation with the factor) assumes features are independent, which gives p-values that are too small for
   pathways whose features vary together; it stays selectable, with that assumption stated. Existing scripts that call
   `integrationEnrichFactor()` without naming a test now get the new test and different p-values (the size of the
   change is given by the check on simulated data below; on the synthetic test corpus, 40 pathways of randomly chosen
   measured features were called by neither test). In a check on simulated data shaped like a 12-sample study, on the
   main factor, unrelated sets whose features vary together at a correlation of 0.3 were called at the 5 % level by
   the earlier test in 13.5 % (sets of 5 to 10 metabolites) to 43.0 % (sets of 10 to 100 genes) of cases and by the
   new test in at most 0.6 %. Over every kind of unrelated set checked, the new test called at most 4.4 % on the main
   factor and at most 6.6 % against values unrelated to the data, except sets made only of the most variable features
   (up to 10.1 %); on weaker factors both tests call unrelated sets of the model's features too often (the new test up
   to 12.8 %, the earlier one up to 60.1 %). To get the earlier test, pass `overrides = list(test = "gsea")`; its
   results are unchanged. In the Enrichment tab the new test is pre-selected. Results made before this change keep their test when replayed. Every
   KEGG and Reactome summary now carries a caution (the factor was learned from the same data, so the results are
   exploratory) and the factor's variance explained in the data set. STRING is unchanged and ignores the setting.
   Further test options follow in KIT-191 (VAL-85 (2)): the correlation-adjusted test with negative correlations
   floored, a permutation test and the pre-ranked test, each shown with its measured false-positive rate and power from
   this check.

 * Proteomics enrichment for organisms outside the 13 served by g:Profiler: Gene Ontology sets now include each
   term's more specific terms (parent-term propagation over is a, part of, regulates, positively and negatively
   regulates, as g:Profiler does; is a and part of only is available to scripts as an option). Terms that had fewer
   than 10 measured proteins can now be tested, terms with more than 500 drop out, and parent terms of significant
   terms tend to be significant too. The three top-level terms are not tested. `processEnrichments()` propagates by
   default (`go_propagation = "is_a_part_of_regulates"`; `"none"` gives the previous sets). Analysis records made
   before this version replay without propagation; records made with this version need this version or later. GO
   identifiers that the installed GO release does not contain are counted and reported. Proteins annotated only to
   the three top-level terms leave the background on such a run; their number is recorded. Measured on one proteome,
   the reviewed E. coli K-12 entries of UniProt 2026_04 (4,024 annotated proteins, the background for every setting),
   with synthetic lists, 20 per scenario, the same lists for none, is a and part of, and the five relations: 314,
   1,242 and 1,249 terms in the testable range (10 to 500 annotated proteins); median significant terms (q-value
   below 0.05) for lists made of half of a specific, a mid-level or a broad process term's proteins plus as many
   random proteins 1.5, 20 and 19; 10, 53.5 and 46.5; 12.5, 129.5 and 129. With parent terms the planted term was
   significant in 20 of 20 lists in every scenario; without them the specific and the broad term are not tested at
   all (too few directly annotated proteins) and the mid-level term was significant in 20 of 20. Random lists: median
   0 significant terms; under each setting 1 of 20 lists had a term significant by q-value (by Benjamini-Hochberg 0,
   1 and 1 of 20); without parent terms that list had 120 of its 156 tested terms significant by q-value and none by
   Benjamini-Hochberg, with parent terms 4 by either. These counts describe this proteome and these lists; they are
   not estimates of power or of the false-discovery rate (1 of 20 is compatible with any rate from 0.1 % to 25 %).
   (KIT-173, VAL-65)
 * Proteomics reports: for the 13 organisms served by g:Profiler, the enrichment paragraph now describes g:Profiler.
   It described clusterProfiler with Benjamini-Hochberg correction, because the report did not recognise the recorded
   method. The correction applied is named from the analysis record. On both routes the fold-change bullet now states
   the cutoffs applied (it gave the down cutoff with the wrong sign, or said no cutoff was applied when only the down
   cutoff was set). No result changes. (KIT-173)
 * Proteomics enrichment for organisms outside the 13 served by g:Profiler: the `geneID` cells of the result tables
   now list their proteins alphabetically. Before, the order inside a cell was set by the R process (enrichit's hash
   set over `Rcpp::unique()`, ordered by memory address), so a replay could list the same proteins in another order. No
   number changes: counts, set sizes, ratios, p-values, q-values and the sets themselves are unchanged. (KIT-173,
   VAL-61)


### Scientific changes

 * Lipidomics and metabolomics normalisation: the log2 step is now plain log2(x) by default; before, it was
   log2(x + 1), with the offset described as protecting zeros. Since KIT-118 the import already reads every zero as
   not detected (missing, not imputed) and refuses negative values, so the offset protected nothing; its only effect
   was on positive values, and it depended on the unit of the table: a value v moved by log2(1 + 1/v) (1.44e-4 at
   v = 10,000, 0.0014 at v = 1,000, 1 at v = 1), and a true two-fold change between 0.1 and 0.2 was reported as
   0.126 instead of 1. After this change the log step is unit-equivariant (multiplying a table by a constant shifts
   every log2 value by the same amount). Intensity tables (values of 10,000 and above) change by at most 1.44e-4 log2;
   concentration tables with values near or below 1 change materially, and their earlier fold changes were compressed
   toward 0. Through the default normalisation with RUV-III Automatic, results on such tables change far more than
   the log step alone: on two simulated concentration tables with no true difference between the groups (one data
   set each), 165 and 0 features were called significant with the offset and 16 and 37 with plain log2, all of them
   false. The RUV-III Automatic and Manual modes and the Scale option depend on the unit of the table (a change of
   unit alone changed the significant calls in 20 of 20 simulated data sets) and are followed up separately. The
   "Offset for zeros" field is removed from the normalisation tab. A project whose record used an offset
   replays with that offset unchanged, and the tab says so. Called directly, `logTransformAssays()` now defaults to
   offset 0: a value of 0 or below, or a value that is not finite, becomes missing and is counted in one warning
   (class `omics_log2_nonpositive`); an offset above 0 runs the earlier transformation unchanged; an offset below 0 is
   refused (class `omics_log2_offset_invalid`); an error of the transformation now stops (class `omics_log2_failed`)
   instead of returning untransformed values labelled as log2.
 * Statistical position, stated: zeros in these tables are often values below a detection limit. They are read as
   not detected and not imputed, so a feature is tested only where it was measured (limma on the available values),
   and a fold change can be underestimated when values fall below the detection limit in one group; a feature with
   no value in one of the compared groups is not tested, and the differential analysis summary now counts these
   ("Not tested (no value in a compared group)"). No variance-stabilising transform (generalised log) is offered; the
   test models variance against mean level (limma-trend with robust hyperparameters). Unchanged: a table that is
   already on the log scale is logged again, and tables with negative values are refused at import.
 * The `scale` between-sample option is limma's two-colour M-value scaling (`limma::normalizeMedianAbsValues`): it
   divides each sample's log2 values so that their median absolute value is the same in every sample, so its result
   depends on the unit of the table, and on log2 values below 0 the absolute value treats -2 like +2. It stays
   selectable; when it runs on log2 values below 0, the normalisation log and a notification now say so. Cyclic
   Loess and Quantile do not depend on the unit.
 * Lipidomics and metabolomics filtering figures on log2 values (the workbooks draw them after the log step): a value
   is counted as detected when it is not missing, and the step's missingness counts missing cells only (before, a log2
   value of 0 or below, that is a measured value of 1 or below, was counted as not detected). On intensities a zero
   still counts as not detected.
 * The interactive volcano draws a sample's missing value at 0 (it cannot draw missing values); it now says so when
   the assay holds missing values.


### Scientific changes

 * Label-free (MaxQuant, FragPipe, Spectronaut) and TMT (Proteome Discoverer) analyses now go through the
   protein-level quality control that DIA analyses use. After the design is saved, the step guide moves to Quality
   Control, which has four steps: Accession Cleanup, Protein Intensity Filter, Duplicate Removal and Protein
   Replicate Filter. The guide moves on to Normalization once the Protein Replicate Filter has been applied.
   Normalization can still be opened earlier, and then shows a warning that the quality control steps are not
   complete. Before this change, the app marked quality control as done as soon as the design was saved, so an
   analysis that followed the step guide normalised every protein of the imported table. With RUV at its default, the
   only proteins removed before testing were those failing the missing-value rule in every group; with RUV skipped,
   none. Results made by following the step guide can therefore change. A protein is now tested only if at least two
   groups each have at least two values above the lowest 1 % of all intensities (for groups of equal size). Protein
   identifiers follow the accession the clean-up chooses. Normalisation, RUV and the test are computed on the
   filtered table. With RUV at its default, the clean-up after RUV now applies the same rule to the corrected values,
   so it can remove a few more proteins. On the package's test data (one simulated data set per format: 400 proteins,
   two groups of five, 2 % to 8 % of values missing at random; 40 proteins with a true difference), with RUV at its
   default: proteins tested 398 to 392 (TMT), 398 to 394 (MaxQuant), 397 to 393 (FragPipe); significant at q < 0.05,
   20 to 18 (TMT; 17 to 16 of the 40 true differences, 3 to 2 of the 360 others), 16 to 18 (MaxQuant; 15 to 15, 1 to
   3), 25 to 25 (FragPipe; 22 to 22, 3 to 3). With RUV skipped: proteins tested 400 to 397, 400 to 397 and 400 to
   396; significant 19 to 18 (18 to 17 true, 1 to 1 other), 17 to 16 (16 to 15, 1 to 1) and 25 to 22 (24 to 21, 1 to
   1). Each figure comes from that one data set: a change of one to three proteins is not a measured change in power
   or in the false-discovery rate, and a protein the filters remove counts as not significant. Real data usually have
   more missing values, mostly in low-abundance proteins, so the filters can remove a much larger share. The Quality
   Control tab shows the number of proteins left after each step. Projects and analysis records made before this
   change are not altered; replaying a record repeats the steps it recorded.


### Distribution

 * MultiScholaR now opens normally when it is started from a shell that an editor extension or a coding agent
   opened; before, it closed at once without a window. (Linux and Windows.)
 * Help, About names the commit of the analysis code that was installed, beside the commit the installer was built
   from. Reports and analysis records made with the installed application name the same commit; before, they named
   none.


### Scientific changes

 * MaxQuant import: contaminant protein groups are now removed from proteinGroups.txt files written by MaxQuant 1.4
   and earlier, which name the column `Contaminant` instead of `Potential contaminant` (the column is `Potential
   contaminant` from MaxQuant 1.5.2.8 at the latest); before, these rows were kept although 'Filter contaminants' was
   ticked (44 protein groups in the Cox et al. 2014 benchmark, PXD000279, 22 of them sharing a human accession such
   as keratin or actin with a contaminant entry and removed whole, as MaxQuant instructs). Reports generated before
   this change from such files stated that contaminant rows had been removed; they had not. Until a later change
   corrects it, the label-free report still states this removal when 'Filter contaminants' is unticked or the file
   has no contaminant column. The rows removed per flag are now logged. Existing projects are not changed; import the
   file again to apply the correction. On spike-in standards whose proteins overlap MaxQuant's contaminant list, such
   as UPS1 and UPS2 (gelsolin, serum albumin, transferrins), this rule can remove true analytes; untick 'Filter
   contaminants' to keep them, which also keeps the reverse (decoy) rows.

## MultiScholaR 0.9.0 (staging)

First staging build with every omic in the distributable (Linux and Windows) and the multi-omics integration module (MOFA2) behind its gate; versions 0.9.0.x follow per larger feature, 1.0.0 at acceptance of the JPR paper.

## Development version (before 0.9.0)

### Scientific changes

 * Lipidomics and metabolomics differential abundance: technical replicates declared in the design
   (several runs of one biological unit, the design step's `tech_rep_group`) are now fitted with the
   biological unit as a block (`limma::duplicateCorrelation`, then `lmFit(block, correlation)`; limma
   User's Guide 11.3), the rule proteomics already applies. Before, each injection was tested as an
   independent sample, which made P-values too small (null rate 0.16 at a nominal 0.05 on duplicated
   injections; at most 0.075 now). The DA record names the model (`block_reason`,
   `uses_duplicate_correlation`). Designs without technical replicates give the same results as before.
   The RUV-III replicate structure of these omics is unchanged here.
 * Transcriptomics: technical replicates declared in the design are summed per biological unit at the
   design step (`edgeR::sumTechReps`; edgeR User's Guide 4.9.1), before the expression filter, TMM, RUVs and
   the edgeR QL test, which takes no block. Before, each library was tested as an independent sample (null
   rate 0.29 at a nominal 0.05 on duplicated libraries through the route; at most 0.075 now). On such designs
   the samples are now the units (named by `tech_rep_group`), the summed runs are recorded with the design
   state, results no longer equal the notebook on library counts, and the RUVIII export no longer offers a
   technical-replicate M. Designs without technical replicates are unchanged.
 * A design whose every group is one biological unit measured more than once is refused in every omic: no
   biological replication. Proteomics, lipidomics and metabolomics stop the differential analysis with an error
   (class `omics_da_no_biological_replication`) naming the groups; before, the runs were tested as independent
   samples. Import, design saving and normalisation are unchanged.
 * Transcriptomics design check: a group needs two biological units, as in every omic. A group of one unit
   sequenced twice (two libraries, one `tech_rep_group`) is now refused at save; before, it passed as two samples.
 * Lipidomics and metabolomics differential abundance: an error on one assay is no longer dropped silently.
   `runLipidsDA` and `runMetabolitesDA` name the failed assay in `failed_assays` with a warning and keep the
   other assays, and stop when every assay failed; the DA tab records and warns about a failed assay fit.
 * DIA-NN import: a precursor quantity of exactly 0 (a precursor DIA-NN did not quantify) is now read as a
   missing value in both quantity columns (`Precursor.Quantity` and `Precursor.Normalised`),
   with a notice of the count per column, as the LFQ and TMT readers already did. Before, a 0 counted as an
   observed value in the peptide rollup and lowered the technical-replicate imputation's mean. A negative
   or infinite quantity now refuses the import. Results change only for reports that contain zeros.
 * Lipidomics and metabolomics import: an intensity of exactly 0 in a sample column (a feature the vendor
   software did not detect) is now read as a missing value, with a warning of the count per assay, as
   proteomics does. Before, a 0 was kept as a measured intensity. A negative intensity now refuses the
   import. Transcriptomics counts keep their zeros. Results change only for tables that contain zeros.
 * Enrichment through g:Profiler (`processEnrichments()`): the custom background ids are now stripped of a
   colon suffix as the query ids always were (`P12345:iso1` is read as `P12345` in both). Before, the
   background kept the suffix, so such a query id was missing from its own background. Results change only
   for ids that contain a colon; standard UniProt accessions and gene names are unchanged.
 * `runReactomeEnrichment()`: the organism name is now matched in full, ignoring case, against the species of
   the Reactome mapping file, and a name that matches no species (or several) stops with an error that
   lists the supported names. Before, the name was matched as a substring and the first hit taken
   ("saccharomyces" picked one of two yeasts), and an unmatched name fell back to *Homo sapiens* without
   notice. Results change only where the name given was not a species' full name.
 * Transcriptomics RUV (KIT-129): three designs and a mode without RUV are available to the RUV step beside
   today's grid (the kernels only; no screen offers them yet and no default changes). The factors are
   estimated on upper-quartile normalised counts (`EDASeq::betweenLaneNormalization(which = "upper")`, the
   RUVSeq vignette's usage) by `RUVSeq::RUVs` with every gene as a control and the groups as replicate sets,
   or by `RUVSeq::RUVg` with a supplied list of control genes or with the genes outside the 25% most
   significant of the first test of the groups. The factors enter the edgeR model and the test reads the
   filtered counts, as today; the display matrix for figures and export is still `removeBatchEffect` with
   the factors. The number of factors k is given to the step with the rule that produced it; a k of 0, or the
   mode without RUV, fits `~ 0 + group` and leaves the display matrix uncorrected. A k the data cannot carry
   is refused, never reduced (RUVSeq itself returns fewer factors without notice). A supplied list is taken
   once per id and in the matrix's gene order, as are the first test's controls; an id that is NA or empty is
   refused, every id not found among the filtered genes is warned of, a list with fewer than two genes found
   is refused, and the sha256 of the sorted list is recorded (`ruv_control_list_sha256`). Today's grid, when it is computed, refuses a request that also carries a
   parameter of these designs (`trx_ruv_grid_declared_parameters`) instead of dropping it. Today's design (the grid's
   controls and k, RUVs on the filtered counts) is the default and every result is as before.
 * Proteomics RUV-III correction (`ruvIII_C_Varying()` for protein data, and the adjusted matrix of the
   DIA_LIMPA route): each protein keeps its mean abundance. Before, the corrected values of a protein were
   shifted by a constant that depended on its loadings (a median of 27 log2 on one synthetic DIA dataset).
   Differences between samples within a protein are the same as before in the classic correction. On the
   DIA_LIMPA route the factors are now computed from controls centred on their means, so the factors, the
   adjusted matrix and the results that use them change. Abundance-dependent steps after RUV can change:
   the low-intensity filter that follows the classic correction may now remove a different protein.
   Lipidomics, metabolomics and the peptide method are not changed.
 * RUV, the number of unwanted factors k (KIT-129): an estimated count is available to the RUV step of every
   omic (the kernel only; no screen offers it yet and no default changes). It is `sva::num.sv()` (method "be",
   20 permutations, a fixed and recorded seed) on the features observed in every sample: the groups are
   removed, and a component of the residuals is counted while at most 10% of the permutations reach its share
   of the residual variance. `sva` is a suggested package; without it the estimate is refused with the reason
   and k has to be given, never replaced by another rule. The estimate is also refused when fewer than 20
   features are observed in every sample or the samples leave fewer than two residual components; a count
   of 0 means no RUV for the run, with a notice. A k above what the design can carry (samples minus replicate
   sets, a residual degree of freedom beside the groups, fewer factors than controls) is refused, not
   reduced. The estimate has its own diagnostic plot (observed share of residual variance per component
   against the permutation null, the count marked). Today's rule on the canonical-correlation plot is
   unchanged and stays the rule of today's designs; it is not defined when all features are the controls.
   The estimate draws under a fixed and recorded random-number kind, so it does not depend on the session's
   `RNGkind()`; a sample without a group and a single group are refusals with their reason; a given k is
   used whatever error the estimate raises, the reason recorded. The count removes the groups only: on a
   design with declared technical replicates its record says so (`technical_replicates_ignored`).
 * DIA_LIMPA RUV covariates (KIT-129): the control rule and the rule for k are now parameters of the step
   (the kernel only; no screen offers them yet and no default changes). Besides today's sweep (the limpa
   script's rule: the least significant proteins as controls, k from the canonical-correlation plot), the
   step can take all proteins as controls, a supplied list of control proteins, or the proteins outside the
   25% most significant of a first test of the groups (`limpa::dpcDE` without RUV, the F test of the
   differences between the groups); k is the estimated count, today's estimate where a separate control set
   exists, or a given number. These designs depart from the limpa script; they run in the covariates mode
   only, the quantified matrix is unchanged and the factors enter the `dpcDE` model as today. Declared
   technical replicates stay RUV-III's replicate sets. When the estimated count is 0, or the data cannot
   give a k, the run has no RUV and says why. With the defaults every result is as before.
   Such a run is recorded as one: mode `none` and `ruv_applied = FALSE` in the product, the state and the
   commit, with the reason; no percentage, controls or control count are given as used, and the integration
   input labels its matrix `normalised`. A declared run's commit names its control rule, its control count,
   k and the rule that gave it; the sweep's own record, when the sweep gave k, is kept apart
   (`ruv_result$k_sweep`). Supplied control ids that are not in the protein table raise a warning with the
   matched and listed counts. Today's declaration refuses a request that carries a rule parameter.
 * DIA_LIMPA differential abundance (KIT-129): a model that leaves no residual degree of freedom (the
   formula's terms and the RUV factors against the samples) is refused with a message that names k, instead
   of failing inside the fit with "NA/NaN/Inf in 'y'".
 * Proteomics automatic RUV parameters (`findBestNegCtrlPercentage()`): the number of factors k is now chosen
   among the values up to the sample-size limit (`max_acceptable_k`). Before, the limit was applied only to
   the penalty, and a larger k could be returned. Results change only where the canonical-correlation curve
   pointed above the limit.
 * Classic proteomics (DIA with IQ, LFQ, TMT) can account for RUV in the test itself (KIT-126): the
   differential-abundance model reads the normalised, uncorrected matrix with the RUV-III factors as
   covariates, the groups as replicate sets and every protein as a negative control. The factors are
   computed on the proteins observed in every sample (validated to 25% missing values); the RUV-III-C
   corrected matrix is still produced for figures, export and integration. Every protein of the normalised
   matrix is tested, also the few that RUV-III-C cannot correct: those are named with the result and have no
   corrected intensity. The test is refused, with the reason, when the samples left cannot carry the model. The number of factors k is the
   automatic sweep's estimate, reduced to what the samples and complete proteins can carry, or the value
   given; its source is recorded with the result. The design is available to the normalisation step
   (`ruv_design = "groups_all_features"`) and is not selected by default or from the application yet: every
   result is unchanged until it is chosen. Today's design (the corrected matrix tested) stays as it is.
 * Lipidomics and metabolomics RUV: a second design is available to `runLipidsDA()` and `runMetabolitesDA()`
   through the new argument `ruv_factors` (the kernel only; no screen offers it yet and no default changes).
   The RUV-III factors are computed per assay from the features observed in every sample, with the groups as
   replicate sets and either all those features or a supplied list (internal standards) as controls, and the
   differential test reads the uncorrected values with the factors in the limma model instead of testing the
   corrected matrix. Supplied controls are not tested; a feature with too few values to estimate the
   factors' coefficients is not tested and is named in the result; a feature with no value in one whole group
   is tested between the groups where it is observed and listed as not observed in that group. The number of
   factors has no default: it is given, or (with supplied controls only) taken from the automatic sweep and
   limited by the sample-size limit, and is recorded with its source. An assay whose factors cannot be
   computed is analysed without RUV, with a warning; a model left without a residual degree of freedom, or the
   corrected object given together with the factors, stops. The corrected matrix of this design
   (for figures, export and integration) keeps each feature's mean abundance. Without `ruv_factors` every
   result is as before: `ruvIII_C_Varying()` for lipid and metabolite data is not changed and still shifts
   each feature by a constant.

### Other changes

 * The fast DIA_LIMPA engine (`fastlimpadpc`) is pinned at fast_limpa_dpc 0ad8c527 (KIT-157), which compiles its
   Rust core for the x86-64-v3 baseline (AVX2, FMA) instead of the building machine's CPU
   (`-C target-cpu=native`). A Linux bundle built on an AVX-512 machine no longer stops with "illegal instruction"
   on CPUs without AVX-512. Windows builds had no CPU flag and are unchanged.
 * The managed worker no longer needs systemd (KIT-131). A managed application (the internal launcher a desktop
   shell starts R through) runs each heavy step as a child process of the application where systemd is absent (Windows, or Linux without `systemctl`), or when the
   option `multischolar.workflow_launcher = "process"` is set. The worker keeps the per-step deadline, stops with the
   application and with its whole process tree, and returns its result the same way. Without systemd there is no
   memory limit on the worker and no memory measurement. Linux with systemd is unchanged. `processx (>= 3.9.0)` and `ps`
   are now declared imports (both were already installed through `callr`).
 * The packaged desktop application (KIT-134) has its own entry, `runPackagedApplication()` (internal), which the
   bundle's shell starts. It takes the launcher the shell's start-up probe chose and tries, in order: the systemd
   worker; the plain launcher, after one trial child process; the application inside the entry's own R process. It
   tells the shell on one `MSR-LAUNCHER` line which launcher runs and why, stops R when the shell is gone, and never
   installs or downloads anything. Under systemd, a unit started from inside the bundle now also receives the bundle's
   `RHOME`, `LD_LIBRARY_PATH`, library, Python, pandoc, `PATH` and locale; outside a bundle the unit's arguments are
   unchanged. `later` is now a declared import (it was already installed through `shiny`).
 * `fastduckdb (>= 0.2.0)` is now a declared import (it was in Suggests). The lipidomics and metabolomics intensity
   filter, the DIA precursor rollup and the direct table writer need it, so an installation without it could not run
   those steps. The pinned source is unchanged (DESCRIPTION's Remotes). No function or result changes.
 * The on-disk state owner compares paths with forward slashes on every platform, so its containment checks no longer
   refuse valid project paths on Windows. Nothing changes on Linux or macOS.
 * Windows: a session archive written without GNU tar (every DIA session export, and the archive refresh of the Summary
   save) no longer leaves a full copy of its files in the temporary folder. The staging folder was removed while it was
   still R's working directory, which Windows refuses (KIT-156, D-1). The archive itself is unchanged.

## Version 0.5.0 (Release Candidate)

This release candidate captures the large janitor-branch stabilization campaign:
substantial de-monolithing of oversized GUI and helper modules, parity auditing
against the pinned `main` baseline, coverage-backed verification across the refactor
surface, and a correctness-focused audit of automatic RUV parameter selection.

Final promotion from release candidate to stable release is intended to follow the
remaining end-to-end GUI suite validation across the default proteomics,
metabolomics, and lipidomics workflows.

### Architecture, De-Monolithing, and Helper Extraction

The codebase underwent a broad stabilization pass focused on reducing god modules,
surfacing hidden logic, and making wrapper behavior auditable.

-- Wrapper and module families stabilized
 * Proteomics wrappers and workflows were stabilized across import, design, builder,
   normalization, QC, enrichment, and summary flows.
 * Metabolomics wrappers and workflows were stabilized across import, QC,
   normalization, differential abundance, design, and summary flows.
 * Lipidomics wrappers and workflows were stabilized across import, QC, design,
   normalization, differential abundance, and summary flows.
 * Shared infrastructure was stabilized in the general file-management and plotting
   families, alongside S4 object surfaces in proteomics, peptide, metabolomics, and
   lipidomics code paths.

-- Functions-within-functions reduced
 * Refactor work lifted implicit nested logic into named helper surfaces where
   possible, making behavior easier to test, replay, and compare against the
   baseline branch.
 * Legacy wrapper shells and helper families were reorganized into explicit bundle
   seams so behavior could be reviewed at the contract, helper, and lineage level
   rather than only as monolithic files.

### Refactor Fidelity, Parity, and Coverage Campaign

This branch was not treated as a cosmetic refactor. It was audited against a pinned
`main` baseline and then driven through a dedicated fidelity workflow.

-- Closeout and parity gates
 * The integrated refactor closeout reported `proven parity: true` against pinned
   baseline `main@326c049`.
 * Surface, manifest, behavior, contract, and high-severity gates all passed in the
   latest closeout run.
 * Open exceptions closed to `0`, with `0` high-severity open exceptions at closeout.

-- Bundle- and contract-level verification
 * The fidelity bundle map covered `233` bundles spanning wrapper entrypoints, helper
   surfaces, S4 method surfaces, and lineage families.
 * The latest targeted coverage probe executed `1161` tests across `124` test files
   and recorded `55983 / 70076` covered lines (`79.9%`) in the audited target set.
 * Coverage evidence recorded `233` passing bundle gates, `227` bundles explicitly
   justified by tests, `0` low-coverage bundles, and `0` regressed bundles.

-- Per-omic coverage characterization added
 * Coverage characterization was added for metabolomics normalization, import, and
   differential abundance lineages.
 * Coverage characterization was added for proteomics import and design-builder
   lineages.
 * Coverage characterization was added for lipidomics normalization, differential
   abundance, and design-builder lineages.

### Automatic RUV Parameter Optimization — Integrity Audit and Corrections

This release candidate completes a correctness overhaul of automatic RUV-III-C
parameter optimization across the omics pipelines. **Automatic-mode results may
differ from previous versions for the same dataset.** Manual mode semantics are
intended to remain unchanged.

-- Audit scope and outcome
 * Automatic RUV behavior was audited across proteomics, peptide, metabolomics, and
   lipidomics implementations.
 * The audit verified multiple correctness defects in shared K-selection logic,
   objective alignment, percentage-search behavior, and peptide optimization flow.
 * The branch includes the follow-through fix sequence plus a regression matrix and
   release-gate pass dedicated to the RUV work.

-- Automatic K selection (`findBestKElbow`)
 * K selection now uses a first-plateau rule instead of argmax-style winner
   selection. The optimizer chooses the smallest K whose separation score lies
   within 5% of the best observed value, favoring the earliest stable correction
   level over late, marginal gains.
 * Shared K-selection, scoring, and curve-alignment logic were hardened so the
   optimizer behaves deterministically and the selected K is consistent with the
   actual comparison curve being evaluated.
 * `findBestK()` is deprecated in favor of `findBestKElbow()`. The legacy entrypoint
   delegates and emits a deprecation warning.

-- Proteomics automatic mode
 * Proteomics automatic optimization now uses deterministic tie-breaking and richer
   optimization traces, improving reproducibility and post-run diagnostics.
 * Shared scoring updates reduce ambiguity when multiple candidate settings produce
   near-identical separation behavior.

-- Metabolomics and lipidomics automatic mode
 * Automatic mode now genuinely searches the user-requested percentage range. Each
   candidate percentage is evaluated as a whole-object optimization step instead of
   silently collapsing to `percentage_max`.
 * This change aligns actual behavior with the documented automatic-mode contract and
   removes the previous "single-point search disguised as optimization" behavior.

-- Peptide automatic mode
 * Peptide automatic optimization now evaluates candidates with the full
   `ruvCancor()` objective rather than the `ruvCancorFast()` surrogate during winner
   selection.
 * This reconciles the optimizer with the objective used to interpret the final
   result and reduces the chance of selecting controls on one scoring regime and
   applying them under another.
 * `ruvCancorFast()` is deprecated for optimizer selection use.

-- Deprecations and user-facing behavior
 * The `weighted_difference` separation metric is now deprecated because it pushes
   high-K selection in the opposite direction from the composite score penalty.
 * UI labels and documentation now identify deprecated RUV choices explicitly.
 * Workbook and release documentation were updated to describe the corrected
   automatic-mode behavior and its implications for result comparability.

### Documentation and Release-Engineering Cleanup

-- Documentation refresh
 * Manual documentation was reorganized into `docs/manual/`, replacing the previous
   root-level `HowTo/` placement.
 * RUV release notes, workbook text, and deprecation language were refreshed to
   match current optimizer behavior.

-- Repository hygiene
 * Local workflow artifacts, refactor audit outputs, development planning material,
   `renv/`, `renv.lock`, and local workbook scratch content were moved out of Git
   tracking and into ignore rules where appropriate.
 * This keeps the repository focused on source, docs, and releasable assets rather
   than machine-local refactor state.

## Version 0.4.1.2

-- All Modules & Infrastructure
 * Project-wide relicensing to **LGPL v3.0** for open-source compliance
 * MultiScholaR Launcher scripts now correctly skip branch selection when started in `--local` mode
 * Standardized differential abundance parameters (using `da_` prefix globally)
 * Fixed numerical parsing errors where q-value thresholds were incorrectly evaluated during Volcano plot generation
 * **New Feature**: Added a "Sanitize Sample Names" option in the **Setup & Import** step across all pipelines. This uses `janitor::make_clean_names()` to proactively clean sample IDs (e.g., "Sample-A!" -> "sample_a") at the point of entry for better compatibility.
 * Standardized composite QC figure generation across all modules to use the `savePlot()` framework, ensuring simultaneous `.pdf`, `.png`, and `.rds` exports.

-- Proteomics
 * improved interactive volcano plot
 * Parquet format input
 * Peptide filtering - Group-aware peptide intensity filtering instead of % of samples higher than intensity theshold filtering
 * General bug fixes and refactoring
 * Included functions for peptide-base missing values imputation using the `limpa` R library
 * **Robust IQ Rollup**: Implemented internal sample aliasing and dynamic design matrix filtering to prevent name mangling and handle samples dropped during quality filtering.

-- Metabolomics
 * improved interactive volcano plot
 * Fixed a design matrix creation bug related to sample assignment indexing
 * Improved composite QC figure generation using `savePlot()`.

-- Lipidomics
 * improved interactive volcano plot
 * Improved composite QC figure generation using `savePlot()`.

## Version 0.4.1.1

### Heatmap Visualization Improvements

#### Standardized Heatmap Saving

- **Manual Save Requirement**: Heatmaps are now saved manually via a new "Save Heatmap" button in the Heatmap tab, rather than automatically during DE analysis. This allows users to fine-tune clustering and aesthetics before exporting.
- **Warning Banner**: Added an informative banner to the results tab in all omics modules (Proteomics, Metabolomics, Lipidomics) to remind users that heatmaps must be saved manually.
- **Reproducibility Artifacts**: The save function now generates:
  - High-resolution `.png` and `.pdf` images.
  - A `_heatmap_methods.md` file documenting the exact parameters used (clustering method, distance metric, etc.).
  - A `_clusters.csv` file containing the cluster membership for every molecule in the heatmap.

#### Enhanced Clustering Options

- **Tree Cutting**: Added explicit "Tree Cutting" controls to the heatmap UI, allowing users to define clusters by:
  - **K Clusters**: Arbitrary number of clusters.
  - **Height Cutoff**: Dendrogram height.
  - **Dynamic**: Dynamic tree cutting (if available).
- **Cluster Summary**: Added a "Cluster Information" panel that displays the size and members of each cluster interactively.

#### Summary Report Integration

- **Automated Copying**: Updated `copyToResultsSummary` path management to ensure that manually saved heatmaps in `publication_graphs/Heatmap` are automatically included in the final "Publication Figures" export directory.

## Version 0.4.1 (Previous Stable Release)

### Major Feature: Lipidomics & Metabolomics GUI Integration

This release marks a significant milestone, introducing full GUI support for both Lipidomics and Metabolomics workflows, bringing them to parity with the existing Proteomics pipeline.

#### Comprehensive Workflow Chain

Both Lipidomics and Metabolomics now follow the standard MultiScholaR analysis chain:
**Import -> Quality Control -> Normalization -> Differential Expression -> Enrichment -> Reporting**

#### Data Structures & S4 Classes

- **New S4 Classes**: Introduced `LipidomicsAssayData` and `MetabolomicsAssayData` classes to robustly handle assay-specific data, metadata, and processing states.
- **Class Definitions**: These classes enforce strict structural validity, ensuring data integrity throughout the analysis pipeline.

#### Data Import & Vendor Support

- **Native Vendor Support**:
  - **Metabolomics**: MS-Dial output formats.
  - **Lipidomics**: MS-Dial and LipidSearch output formats.
- **Universal Import Engine**: A flexible new import system allows users to map custom annotation fields (e.g. retention time, adduct, formula) to internal standardized columns, enabling support for virtually any vendor output (Beta).
- **Import Validation**: The import module performs rigorous integrity checks, ensuring:
  - Valid unique identifiers (LipidName/MetaboliteName).
  - Consistency between intensity data and sample lists.
  - Automatic cleaning of invalid or empty rows.

#### Quality Control Chain

The dedicated **QC & Filtering** module provides a comprehensive, interactive cleaning pipeline:

1.  **Intensity Filtering**: Interactive filtering based on feature coverage (missing values) within experimental groups or across the entire dataset.
2.  **Duplicate Resolution**: Tools to identify and resolve duplicate feature entries (e.g., isomers or multiple retention times for the same feature) by retaining the most abundant or highest quality peak.
3.  **Internal Standard (ITSD) QC**: Evaluation of internal standard performance, allowing users to select appropriate standards for normalization and visualize their stability across samples.
4.  **Finalization**: Creation of the cleaned, rigorous S4 object ready for normalization.

#### Multi-Omics Integration

- **Automated Export**: The `mod_lipid_summary` and `mod_metab_summary` modules now automatically export the final processed S4 objects to a dedicated `integration/` directory.
- **Integration-Ready**: These exported objects are pre-formatted for seamless loading into downstream multi-omics integration tools (e.g., MOFA2).

### Improvements & Bug Fixes

- **Report Generation**: Fixed a bug in `lipidomics_report.rmd` where the wrong parameter extraction function was called, ensuring all study parameters are now correctly reflected in the generated reports.
- **Documentation**: Comprehensive refactor of `README.md` to reflect the expanded multi-omics capabilities.

---

## Version 0.3.6.3 (Current Release)

### Enrichment Analysis Module Fixes

#### "No Contrasts Available" Bug Fix

- **Fixed strict state checking**: Enrichment module was only accepting `"correlation_filtered"` state, rejecting other valid states like `"normalized"`, `"ruv_corrected"`, `"protein_replicate_filtered"`
- **Added backup observer**: New observer watches `workflow_data$de_analysis_results_list` directly, ensuring contrasts populate as soon as DE analysis completes regardless of tab state
- **Improved state condition**: Now uses `current_state %in% valid_states_for_enrichment` instead of exact string match

#### S4 Object NULL Crash Fix

- **Multiple S4 source fallbacks**: When `state_manager$getState()` returns NULL, now tries:
  1. DE results `theObject` (`de_results_list[[1]]$theObject`)
  2. Combined DE results structure (`de_results_list$theObject`)
- **Safe slot access**: All `@` slot accesses now wrapped with NULL guards and `tryCatch` with sensible fallbacks
- **Prevents**: `"no applicable method for '@' applied to an object of class 'NULL'"` errors

#### Multi-Species Organism Filtering Fix

- **Fixed organism name to taxon ID mapping**: Now handles `Organism` column containing names like `"Mus musculus (Mouse)"` instead of numeric taxon IDs
- **Built-in organism lookup table**: Maps common organism names to NCBI taxon IDs (human, mouse, rat, zebrafish, yeast, etc.)
- **Improved column detection**: Checks multiple column name variations for both organism and accession columns
- **Debug logging**: Logs available columns and mapping results for easier troubleshooting
- **Single-species fallback**: If no organism column found, assumes all proteins are from target organism

### Mixed Species FASTA Workflow Integration

#### Session Export/Restore

- **Export**: `mod_prot_norm.R` now includes `mixed_species_analysis` in exported session data and saves as separate `mixed_species_analysis.RDS` file
- **Restore**: `mod_prot_de.R` restores `mixed_species_analysis` from loaded session, maintaining full workflow continuity
- **Auto-enable filtering**: Enrichment module automatically enables organism filter checkbox when mixed species analysis was enabled at import

#### Workflow Parameter Reporting

- **New report sections**: `createWorkflowArgsFromConfig()` in `func_general_filemgmt.R` now includes:
  - "Mixed Species FASTA Analysis" section with organism distribution
  - "Enrichment Analysis - Organism Filtering" section with protein counts before/after filtering

### DE Analysis Output Improvements

#### Enhanced Long Format Output

- **Clearer column naming**: Renamed `left_group` and `right_group` columns to `numerator` and `denominator` respectively, explicitly indicating the contrast direction.
- **Sample-specific columns**: For single-contrast results, sample columns are now named with the specific Sample ID and Group name (e.g., `log2norm.Sample1.GroupA`) instead of generic indices, improving readability.

### Report Generation Improvements

#### Dynamic DE Summary in Reports

- **Automated Text Generation**: Added logic to read the `num_sig_de_molecules.tab` file and generate a natural language summary paragraph in the "Differential Expression Analysis" section of reports.
- **Summary Details**: Paragraph now dynamically lists the number of upregulated and downregulated proteins for each comparison (e.g., "For the comparison X vs Y, there were 10 upregulated and 5 downregulated proteins.").
- **File Management**: Updated `copyToResultsSummary()` to ensure the `de_proteins_num_sig_de_molecules.tab` file is correctly copied to the `Publication_tables` directory for report access.
- **Affected Reports**: Implemented in `LFQ_report.rmd`, `TMT_report.rmd`, and `DIANN_report.rmd`.
- **Cleanup**: Removed static placeholder text from GO Enrichment sections to prevent misleading interpretations in generated reports.

---

## Version 0.3.6.2

### DE Analysis Output Fixes

#### Windows Path Normalization

- **Fixed volcano plots not saving**: Resolved critical bug where paths with double slashes (`C://`) caused `dir.create()` to silently fail on Windows
- **Added path sanitization**: `outputDeResultsAllContrasts()` now normalizes paths using `gsub("//+", "/")` and `normalizePath()` before directory creation
- **Improved error handling**: Directory creation now verifies success and provides informative error messages

#### NumSigDeMolecules Figure Output

- **Fixed missing NumSigDE barplot**: The `outputDeResultsAllContrasts()` function was missing the code to save the "Number of Significant DE Molecules" figure
- **Complete implementation**: Added logic to aggregate `num_sig_de_molecules_first_go` tables from all contrasts
- **Output files**: Now writes combined table as `.tab` and `.xlsx`, plus faceted barplot as PNG and PDF to `publication_graphs/NumSigDeMolecules/`

#### Mixed FASTA Gene Name Lookup

- **Fixed empty gene names**: Enhanced `getGlimmaVolcanoProteomicsWidget()` to handle mixed-species FASTA databases
- **Base accession matching**: Now strips isoform suffixes (e.g., `P12345-2` → `P12345`) before annotation lookup, improving match rates
- **Improved fallback logic**: Fixed to handle empty strings `""` in addition to `NA` values
- **Smart fallback**: Uses protein accession as display name when no gene name annotation is found

---

## Version 0.3.6.1

### Critical Memory Optimization

#### Root Cause Analysis: Environment Capture

- **Identified hidden memory bloat**: R's internal memory reporting (`pryr::mem_used()`, `gc()`) was showing stable memory while actual process memory (Task Manager) spiked 20GB+
- **Root cause discovered**: ggplot objects and tidyverse pipe operations capture their **entire parent environment**, creating massive hidden memory consumption that R doesn't report
- **S4 slot assignment triggers full object copies**: Modifying any S4 slot can duplicate the entire object in memory

#### New Memory Tracking Utilities (`R/func_general_helpers.R`)

- **`getProcessMemoryMB()`**: Returns ACTUAL process memory from OS (what Task Manager shows)
- **`getRHeapMemoryMB()`**: Returns R's internal heap memory
- **`checkMemoryBoth(step, context)`**: Logs BOTH values + calculates "hidden" memory from environment capture
- **`reportMemoryDelta(baseline, step, context)`**: Shows memory change since baseline with automatic warnings for >500MB or >1GB increases
- **Cross-platform support**: Works on Windows (via `tasklist`) and Unix/Mac (via `ps`)

#### Pearson Correlation Optimization

- **Rewrote `pearsonCorForSamplePairs()`**: Replaced memory-intensive `pivot_longer` + `left_join` approach with efficient matrix algebra using `cor()`
- **Performance improvement**: Correlation calculation now completes in <1 second instead of 20-30 seconds
- **Memory reduction**: From ~23GB intermediate objects to <100MB

#### Plot Memory Management

- **Disk-first approach**: QC plots now save to disk immediately after generation, then clear from memory
- **New `generateCompositeFromFiles()`**: Stitches saved PNG files using `cowplot::draw_image()` instead of holding all plot objects in memory
- **Removed QC_composite_figure S4 storage**: Eliminated massive S4 object that was duplicating all plots

#### Base R Optimization

- **Rewrote `cleanDesignMatrix()`**: Replaced tidyverse pipe chains with base R operations to prevent environment capture
- **Explicit memory cleanup**: Added `rm()` and `gc()` calls after large object operations

#### Instrumented Functions

The following functions now have comprehensive memory tracking showing both R-heap and actual process memory:

- `plotPcaHelper` (R/func_general_plotting.R)
- `removeRowsWithMissingValuesPercent` (R/func_prot_s4_objects.R)
- `filterSamplesByProteinCorrelationThreshold` (R/func_prot_s4_objects.R)
- `cleanDesignMatrix` (R/func_prot_s4_objects.R)
- `pearsonCorForSamplePairs` (R/func_prot_s4_objects.R)

#### Result

- **Eliminated `std::bad_alloc` crashes**: Normalization workflow now completes successfully without memory exhaustion
- **Stable memory profile**: Memory usage stays within reasonable bounds throughout the workflow
- **Session export works**: `saveRDS()` no longer triggers 50GB+ serialization attempts

---

## Version 0.3.6

### Data Import Enhancements

#### Mixed Species FASTA Support (All Workflows)

- Added **mixed species FASTA analysis** for non-specific search databases
- **Available for all proteomics workflows**: DIA-NN, TMT-PD, LFQ-Fragpipe
- **"Multiple organism" checkbox**: Users can enable mixed species analysis by checking the "Mixed species FASTA (analyze organism distribution)" checkbox during import
- **Smart distribution detection**: When the checkbox is enabled, automatically analyzes organism distribution from FASTA headers (OS/OX fields) and protein matches
- **Interactive organism selection**: Modal dialog displays organism distribution table with protein counts and percentages for each detected organism
- **Smart organism detection**: Automatically identifies the most abundant organism based on protein matches
- **Optional data filtering**: Users can choose to filter data to keep only proteins from the selected primary organism
- **UI improvements**: Organism information inputs (Step 3) are automatically disabled when mixed species mode is enabled via the checkbox
- **Processing feedback**: Enhanced processing modal with spinner and status updates during organism distribution analysis
- **Debug output**: Comprehensive terminal logging for troubleshooting import issues

### Design Matrix Builder Improvements

#### New "Remove Samples" Tab

- Added ability to **exclude samples from analysis** directly in the Design Matrix Builder
- Select multiple samples for removal with a single action
- Removed samples are excluded from the design matrix display and all downstream analysis
- **Reversible**: Restore removed samples via Settings tab → Reset Scope → "Removed Samples Only"
- Follows the same reset pattern as other builder operations for consistency

#### Improved Save Design UX

- Added **immediate visual feedback** when clicking "Save Design" button
- Spinning icon with "Saving design matrix and preparing data..." message appears instantly
- Eliminates confusion from perceived UI freeze during file writing and S4 object creation
- Modal transitions smoothly to UniProt annotation progress bar

#### Batch Auto-Assignment

- **Automatic batch detection**: Detects "Batch" column from TMT workflow imports and pre-populates design matrix
- **Cross-batch comparisons enabled**: Batch is kept separate from experimental factors, allowing comparisons across batches (e.g., Treatment vs Control across different batches)
- **Persistent during resets**: Batch column is preserved during design matrix resets since it originates from data import
- Eliminates manual batch assignment for TMT workflows

#### Technical Replicate Group Creation

- **Automatic tech_rep_group identifier**: Creates unique technical replicate group identifiers combining group and replicate numbers
- **Proper imputation grouping**: Enables correct technical replicate grouping for missing value imputation
- **Smart logic**: When technical replicates exist, samples share the same tech_rep_group; when no tech reps exist, each sample has a unique tech_rep_group (prevents inappropriate cross-sample imputation)
- Applied automatically during both design matrix creation and import operations

### RUV Normalization Improvements

#### Pool/QC Sample Handling

- **Smart Pool/QC exclusion**: Pool/QC samples are now **excluded from ANOVA** (negative control selection) and **Pearson correlation** (QC plots) but **INCLUDED in RUV-III correction** (noise estimation)
- **Rationale**: Pool/QC samples should not influence negative control selection or within-group correlation calculations but can help estimate unwanted variation patterns
- **Automatic detection**: Case-insensitive detection of groups containing "Pool" or "QC" in group names
- **Comprehensive logging**: Detailed messages indicating which samples are excluded and why
- **Functions affected**: `getNegCtrlProtAnovaHelper()`, `pearsonCorForSamplePairs()`, `plotPearson()`, `findBestNegCtrlPercentage()`, `ruvCancor()`, `ruvIII_C_Varying()`
- **Validation**: Ensures sufficient non-Pool samples remain for analysis

#### Improved Logging and User Feedback

- **Detailed exclusion logging**: Logs show exactly which samples are excluded and the reason (NA groups, Pool/QC samples, etc.)
- **Sample name reporting**: Lists excluded sample names (up to 5) with count if more samples are excluded
- **Progress information**: Provides counts of samples before and after filtering operations
- **Workflow context**: Log messages include function names and workflow context for easier debugging
- **User-friendly notifications**: Clear notifications in Shiny UI when samples are excluded or operations complete

### Bug Fixes

#### UniProt Annotation

- **Fixed API timeout issue**: Resolved a critical bug where malformed protein IDs (containing full FASTA headers) were being sent to the UniProt API, causing 400 Bad Request errors and timeouts during the annotation step.
- **Robust ID extraction**: Implemented improved parsing logic to correctly extract clean accession IDs from various FASTA header formats (UniProt, generic pipe, space-separated, and colon-prefixed Ensembl IDs) before API queries.

---

## Version 0.3.5

### Major Release: End-to-End GUI Implementation

MultiScholaR v0.3.5 represents a **transformative release** that transitions the package from a collection of R Markdown workflows to a **production-grade, end-to-end graphical user interface** for multi-omics analysis. This release makes advanced proteomics analysis accessible to users without programming expertise while maintaining the transparency and educational value that defines MultiScholaR.

### New Features

#### End-to-End Multiomics GUI

- Complete Shiny application built using the `{golem}` framework
- Production-grade modular architecture following "app-as-a-package" philosophy
- Centralized state management using R6/S4 hybrid architecture
- Workflow-based navigation through all analysis stages
- Interactive data import for multiple proteomics platforms
- Real-time quality control visualizations
- Interactive experimental design builder
- Automated normalization with visual feedback
- Integrated differential expression analysis
- Enrichment analysis with multiple annotation sources

#### Three Complete Proteomics Workflows

- **DIA-NN Workflow**: Full GUI implementation with peptide-to-protein rollup using IQ tool, comprehensive QC pipeline, automated normalization with RUV-IIIc optimization
- **TMT-PD Workflow**: Complete TMT (Tandem Mass Tag) analysis pipeline with protein-level analysis from start, TMT-specific QC procedures
- **LFQ-Fragpipe Workflow**: Label-free quantification workflow with FragPipe-specific data handling and optimization

#### Automated RUV-IIIc Parameter Optimization

- Novel algorithm (`findBestNegCtrlPercentage` and `findBestNegCtrlPercentagePeptides`) eliminates manual parameter tuning
- Systematic parameter testing: Evaluates 1-30% negative control percentages
- Optimal factor selection using canonical correlation analysis
- Composite scoring system balancing separation quality with complexity
- Fully automated and dataset-specific adaptation

#### Publication-Ready Reports

- Automated report generation (HTML and Word formats)
- Workflow-specific templates: `DIANN_report.rmd`, `TMT_report.rmd`, `LFQ_report.rmd`
- Comprehensive analysis documentation with all parameters and results
- Includes QC plots, statistical results, enrichment analyses, and methodology

#### Metabolomics Workflow (R Markdown)

- Fully functional metabolomics analysis workflow
- LC-MS and GC-MS support
- Comprehensive QC pipeline
- Multiple normalization methods including RUVIIIc
- Metabolite pathway enrichment
- Two workflow templates: Starter and experienced versions
- Note: GUI integration planned for v0.4

#### Enhanced R Markdown Workflows

- DIA-NN with limpa imputation: Advanced missing value handling
- Workflow variants: Starter (educational) and experienced (streamlined) versions
- Improved documentation and educational content

### Architecture Improvements

- Modular Shiny framework with self-contained modules for each workflow stage
- State management system using WorkflowState R6 class
- Enhanced project structure with `inst/app/`, `inst/workbooks/`, `inst/reports/`, and `dev/` directories
- Function-based organization with reusable core functions
- Easy extension capability for new omics types

### User Experience Enhancements

- MultiScholaR Launcher: Platform-specific launchers (Windows/macOS) with automatic setup
- Comprehensive README update with GUI instructions
- Two usage options: GUI (recommended) and R Markdown (advanced)
- Clear workflow diagrams using Mermaid
- Detailed architecture documentation

### Technical Improvements

- Improved code organization with modular architecture
- Enhanced quality control with interactive visualizations
- Advanced normalization with automated optimization
- Better error handling in workflow state management
- Enhanced data validation during import
- Improved visualization quality and performance
- Enhanced logging and debugging capabilities

### Breaking Changes

- Banner image location moved from `images/MultiScholaR.png` to `inst/shiny/www/MultiScholaR.png`
- Some workflow files renamed for consistency (e.g., `DIA_workflow_beginner.rmd` → `DIA_workflow_starter.rmd`)
- New package structure with `inst/app/` and `inst/reports/` directories

### Dependencies Added

- `{golem}`: Shiny application framework
- `DT`: Interactive data tables
- `shinyjqui`: Enhanced UI components
- Additional Shiny ecosystem packages

### Documentation

- Comprehensive README rewrite with GUI focus
- Architecture documentation with detailed explanation of modular design
- Workflow diagrams showing visual representation of analysis pipeline
- Usage instructions for both GUI and R Markdown users

---

## Version 0.2

### Initial Multiomics Functionality

- Core multiomics workflow implementation
- Cross-omics data linking
- Factor analysis methods (MOFA) with visualization
- Pathway/StringDB integration
- Multiomic pathway visualisations
- Basic R Markdown workflows for proteomics
