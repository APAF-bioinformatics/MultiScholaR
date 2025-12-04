# MultiScholaR NEWS

## Version 0.3.6 (Development)

### Data Import Enhancements

#### Mixed Species FASTA Support
- Added **mixed species FASTA analysis** for non-specific search databases
- When enabled, automatically analyzes organism distribution from FASTA headers (OS/OX fields)
- **Interactive organism selection**: Modal dialog displays organism distribution table with protein counts and percentages
- **Smart organism detection**: Automatically identifies the most abundant organism based on protein matches
- **Optional data filtering**: Users can choose to filter data to keep only proteins from the selected primary organism
- **UI improvements**: Organism information inputs (Step 3) are automatically disabled when mixed species mode is enabled
- **Processing feedback**: Enhanced processing modal with spinner and status updates during data import
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

---

## Version 0.3.5 (Current Release)

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

