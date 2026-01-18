# MultiScholaR NEWS

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

## Version 0.4.1 (Current Release)

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
