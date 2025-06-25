# MultiScholaR Shiny Application

## Overview

This is a modular Shiny application that transforms the MultiScholaR R Markdown workflows into an interactive web interface. The app supports multiple omics types (proteomics, metabolomics, transcriptomics, lipidomics) and allows for integrated multi-omics analysis. 

**NEW ARCHITECTURE (2025)**: The application has been restructured to support **format-specific workflows** through component-based design, enabling different analysis pipelines for DIA, DDA/LFQ, and specialized workflows like glycoproteomics within the same omics type.

## Architectural Design Philosophy

### Component-Based Architecture
The application follows a **clean separation of concerns** with three distinct layers:

1. **UI Components** (`ui/`): Pure UI functions with no server logic
2. **Server Components** (`server/`): Pure server logic with no UI definitions  
3. **Module Coordinators** (`modules/`): Lightweight wrappers that combine UI+server and provide public APIs

### Format-Specific Workflow Support
Each omics type can support multiple analysis formats:
- **Proteomics**: DIA (peptide-centric), DDA/LFQ (protein-only), Glycoproteomics (custom pipeline)
- **Metabolomics**: LC-MS, GC-MS, NMR (different preprocessing needs)
- **Future**: Easy extensibility for new formats within existing omics types

## Directory Structure

```
R/shiny/
├── app.R                    # Main entry point
│
├── ui/                      # Pure UI Components (No server logic)
│   ├── common/             # Shared UI components across omics
│   ├── proteomics/         # Proteomics-specific UI
│   │   ├── workflow_ui.R   # Main proteomics workflow UI
│   │   └── qc/             # QC component UIs
│   │       ├── peptideQCApplet_ui.R      # ✅ Peptide filtering UI
│   │       ├── proteinQCApplet_ui.R      # Protein filtering UI
│   │       ├── diaQC_ui.R               # DIA-specific QC workflow
│   │       ├── ddaQC_ui.R               # DDA/LFQ-specific QC workflow
│   │       └── glycoQC_ui.R             # Glycoproteomics QC workflow
│   ├── metabolomics/       # Metabolomics-specific UI
│   ├── transcriptomics/    # Transcriptomics-specific UI
│   ├── lipidomics/         # Lipidomics-specific UI
│   └── integration/        # Multi-omics integration UI
│
├── server/                  # Pure Server Components (No UI definitions)
│   ├── common/             # Shared server logic across omics
│   ├── proteomics/         # Proteomics-specific server logic
│   │   ├── workflow_server.R # Main proteomics workflow server
│   │   └── qc/             # QC component servers
│   │       ├── peptideQCApplet_server.R  # ✅ Peptide filtering server
│   │       ├── proteinQCApplet_server.R  # Protein filtering server
│   │       ├── diaQC_server.R           # DIA-specific QC server
│   │       ├── ddaQC_server.R           # DDA/LFQ-specific QC server
│   │       └── glycoQC_server.R         # Glycoproteomics QC server
│   ├── metabolomics/       # Metabolomics-specific server logic
│   ├── transcriptomics/    # Transcriptomics-specific server logic
│   ├── lipidomics/         # Lipidomics-specific server logic
│   └── integration/        # Multi-omics integration server logic
│
├── modules/                 # Module Coordinators & Public APIs
│   ├── common/             # Shared module coordinators
│   │   ├── workflowState.R # Workflow state management
│   │   └── designMatrix*.R # Design matrix components
│   └── proteomics/         # Proteomics module coordinators
│       ├── qualityControlApplet.R    # 🔄 QC coordinator (lightweight wrapper)
│       ├── normalizationApplet.R     # Normalization coordinator
│       ├── setupImportApplet.R       # Import coordinator
│       └── differentialExpressionApplet.R # DE coordinator
│
└── www/                     # Static assets (CSS, JS, images)
```

## Component Architecture Details

### UI Components (`ui/`)
Pure UI functions that return Shiny UI elements:
```r
# Example: R/shiny/ui/proteomics/qc/peptideQCApplet_ui.R
peptideQCAppletUI <- function(id) {
  ns <- NS(id)
  # UI definition only - no server logic
  tabsetPanel(...)
}
```

### Server Components (`server/`)
Pure server functions that handle reactive logic:
```r
# Example: R/shiny/server/proteomics/qc/peptideQCApplet_server.R
peptideQCAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  moduleServer(id, function(input, output, session) {
    # Server logic only - no UI definitions
    observeEvent(input$apply_filter, {...})
  })
}
```

### Module Coordinators (`modules/`)
Lightweight wrappers that combine components and provide stable public APIs:
```r
# Example: R/shiny/modules/proteomics/qualityControlApplet.R
source("R/shiny/ui/proteomics/qc/peptideQCApplet_ui.R")
source("R/shiny/server/proteomics/qc/peptideQCApplet_server.R")

qualityControlAppletUI <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    tabPanel("Raw Data Overview", rawDataQCUI(ns("raw"))),
    tabPanel("Peptide Filtering", peptideQCAppletUI(ns("peptide"))),
    tabPanel("Protein Filtering", proteinQCAppletUI(ns("protein")))
  )
}

qualityControlAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  moduleServer(id, function(input, output, session) {
    peptideQCAppletServer("peptide", workflow_data, experiment_paths, omic_type, experiment_label)
    proteinQCAppletServer("protein", workflow_data, experiment_paths, omic_type, experiment_label)
  })
}
```

## Format-Specific Workflow Implementation

### Current Status: Proteomics Format Support

**✅ Completed Components:**
- `peptideQCApplet_ui.R` - All 7 peptide filtering tabs extracted
- `peptideQCApplet_server.R` - Server skeleton created  

**🔄 In Progress:**
- Server logic extraction from existing monolithic `qualityControlApplet.R`
- `proteinQCApplet_ui.R` and `proteinQCApplet_server.R` extraction

**📋 Planned Format Routers:**
1. **DIA Workflow Router**: Raw Data → Peptide QC → Protein QC → Normalization → DE
2. **DDA/LFQ Workflow Router**: Raw Data → Protein QC → Normalization → DE (skip peptide level)
3. **Glycoproteomics Router**: Raw Data → Glycopeptide QC → Site Assignment → Glycoprotein QC → Normalization → DE

### Format Detection & Routing
```r
# Future implementation in workflow_ui.R
output$format_specific_workflow <- renderUI({
  req(workflow_data$detected_format)
  
  switch(workflow_data$detected_format,
    "diann" = diaQCUI(ns("dia_workflow")),
    "spectronaut" = diaQCUI(ns("dia_workflow")), 
    "fragpipe" = ddaQCUI(ns("dda_workflow")),
    "maxquant" = ddaQCUI(ns("dda_workflow")),
    "glyco" = glycoQCUI(ns("glyco_workflow"))
  )
})
```

## Running the App

To run the application:

```r
# Load MultiScholaR
library(MultiScholaR)

# Run the app
shiny::runApp("R/shiny/app.R")
```

## Features

### Landing Page
- Select one or more omics types to analyze
- Integration analysis automatically enabled when 2+ omics are selected
- Clean, intuitive interface for omics selection

### Workflow Tabs
Each omics workflow is divided into sequential tabs:

1. **Setup & Import**: Upload data files, specify organism information, **format detection**
2. **Design Matrix**: Build experimental design using interactive interface (format-independent)
3. **Quality Control**: **Format-specific QC pipelines** with real-time visualization
4. **Normalization**: **Format-aware normalization** methods with before/after comparison
5. **Differential Expression**: Statistical analysis with volcano plots (shared)
6. **Enrichment Analysis**: Pathway and GO term enrichment (shared)
7. **Report**: Generate comprehensive analysis reports (shared)

### Key Features
- **Format-Specific Workflows**: Automatic routing based on detected data format
- **Progress Tracking**: Visual indicators show completed steps
- **Interactive Parameters**: Adjust analysis parameters in real-time
- **Data Persistence**: Results persist between tabs within a session
- **Export Options**: Download processed data, reports, and session info
- **Error Handling**: Graceful error messages with troubleshooting hints
- **Component Reusability**: Mix and match UI/server components for custom workflows

## Development Guide

### Adding a New Format to Existing Omics Type

1. **Create format-specific components**:
```bash
# Example: Adding NMR support to metabolomics
touch R/shiny/ui/metabolomics/qc/nmrQC_ui.R
touch R/shiny/server/metabolomics/qc/nmrQC_server.R
```

2. **Create format router**:
```r
# R/shiny/modules/metabolomics/nmrQCRouter.R
nmrQCRouterUI <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    tabPanel("Raw Data Overview", rawDataQCUI(ns("raw"))),
    tabPanel("NMR-Specific QC", nmrQCUI(ns("nmr_qc"))),
    tabPanel("Metabolite QC", metaboliteQCUI(ns("metabolite_qc")))
  )
}
```

3. **Update main workflow router** to include new format detection

### Adding a New Omics Type

1. Create directory structure:
```bash
mkdir -p R/shiny/ui/new_omics/qc
mkdir -p R/shiny/server/new_omics/qc  
mkdir -p R/shiny/modules/new_omics
```

2. Follow the component pattern established for proteomics

### Component Extraction Pattern

**When extracting monolithic modules:**

1. **UI Extraction**: Move all UI elements to `ui/[omics]/[component]_ui.R`
2. **Server Extraction**: Move all server logic to `server/[omics]/[component]_server.R`  
3. **Create Coordinator**: Update `modules/[omics]/[component].R` to import and coordinate
4. **Maintain API**: Ensure public function signatures remain unchanged
5. **Test Isolation**: Verify components work independently
6. **Update References**: Update any imports/sources in other files

### Module Communication

Modules communicate through:
- `workflow_data`: Shared reactive values across workflow steps
- `experiment_paths`: Standardized directory structure
- Return values from server modules  
- R6 state managers for complex state persistence

### Git Subtree Management

The modular structure supports git subtree workflows:

```bash
# Add a subtree for proteomics components
git subtree add --prefix=R/shiny/ui/proteomics https://github.com/org/proteomics-ui.git main

# Pull updates to specific components
git subtree pull --prefix=R/shiny/server/proteomics https://github.com/org/proteomics-server.git main

# Push component changes
git subtree push --prefix=R/shiny/ui/proteomics https://github.com/org/proteomics-ui.git main
```

## Integration with Existing Functions

The app reuses existing MultiScholaR functions by:
1. Sourcing workflow functions from the appropriate directories
2. Wrapping functions in reactive contexts within server components
3. Handling progress updates and error messages in coordinators
4. Managing intermediate results in reactive values and R6 state managers

## Migration Status & Roadmap

### ✅ Phase 1: Foundation (Completed)
- [x] Component-based architecture established
- [x] Directory structure created (`ui/`, `server/`, `modules/`)
- [x] Peptide QC UI component extracted (7 filtering tabs)
- [x] Peptide QC server skeleton created

### 🔄 Phase 2: Component Extraction (In Progress) 
- [ ] Extract peptide QC server logic from monolithic module
- [ ] Extract protein QC UI and server components
- [ ] Create lightweight coordinator wrappers
- [ ] Test backward compatibility

### 📋 Phase 3: Format-Specific Routing (Planned)
- [ ] Implement DIA workflow router (current full pipeline)
- [ ] Implement DDA/LFQ workflow router (protein-only pipeline)
- [ ] Implement glycoproteomics workflow router (specialized pipeline)
- [ ] Add format detection to Setup & Import step
- [ ] Update main workflow to route based on detected format

### 🚀 Phase 4: Advanced Features (Future)
- [ ] Dynamic workflow composition based on user preferences
- [ ] Custom pipeline builder interface
- [ ] Cross-format comparison capabilities
- [ ] Advanced caching and performance optimization

## Troubleshooting

### Common Issues

1. **"Module not found" errors**
   - Ensure UI/server component files exist in correct directories
   - Check that coordinator modules properly source component files
   - Verify function names match between components and coordinators

2. **Data not persisting between tabs**
   - Verify reactive values are properly initialized in coordinators
   - Check that server components receive workflow_data correctly
   - Ensure state managers are passed to appropriate components

3. **Format-specific workflow not loading**
   - Verify format detection logic in Setup & Import
   - Check that format routers are properly implemented
   - Ensure workflow_data$detected_format is set correctly

4. **Component communication failures**
   - Verify module IDs are consistent between UI and server
   - Check that coordinators pass parameters correctly to components
   - Ensure reactive contexts are properly established

## Future Enhancements

### Near-term (2025)
- [x] Component-based architecture implementation
- [ ] Complete format-specific workflow support
- [ ] Enhanced state management with R6 patterns
- [ ] Improved error handling and user feedback

### Medium-term  
- [ ] Save/load analysis sessions with component state
- [ ] Batch processing capabilities for multiple datasets
- [ ] Advanced visualization options with interactive plots
- [ ] Plugin architecture for custom components

### Long-term
- [ ] Cloud deployment support with scalable architecture
- [ ] API endpoints for programmatic access to components  
- [ ] Real-time collaboration features
- [ ] Machine learning integration for parameter optimization