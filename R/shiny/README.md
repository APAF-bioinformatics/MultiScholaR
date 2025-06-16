# MultiScholaR Shiny Application

## Overview

This is a modular Shiny application that transforms the MultiScholaR R Markdown workflows into an interactive web interface. The app supports multiple omics types (proteomics, metabolomics, transcriptomics, lipidomics) and allows for integrated multi-omics analysis.

## Directory Structure

```
R/shiny/
├── app.R                    # Main entry point
├── ui/                      # UI modules
│   ├── common/             # Shared UI components
│   ├── proteomics/         # Proteomics-specific UI
│   ├── metabolomics/       # Metabolomics-specific UI
│   ├── transcriptomics/    # Transcriptomics-specific UI
│   ├── lipidomics/         # Lipidomics-specific UI
│   └── integration/        # Multi-omics integration UI
├── server/                  # Server modules
│   ├── common/             # Shared server logic
│   ├── proteomics/         # Proteomics-specific server
│   ├── metabolomics/       # Metabolomics-specific server
│   ├── transcriptomics/    # Transcriptomics-specific server
│   ├── lipidomics/         # Lipidomics-specific server
│   └── integration/        # Multi-omics integration server
├── modules/                 # Reusable Shiny modules
│   ├── common/             # Shared modules
│   └── [omics]/            # Omics-specific modules
└── www/                     # Static assets (CSS, JS, images)
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

1. **Setup & Import**: Upload data files, specify organism information
2. **Design Matrix**: Build experimental design using interactive interface
3. **Quality Control**: Apply QC filters with real-time visualization
4. **Normalization**: Multiple normalization methods with before/after comparison
5. **Differential Expression**: Statistical analysis with volcano plots
6. **Enrichment Analysis**: Pathway and GO term enrichment
7. **Report**: Generate comprehensive analysis reports

### Key Features
- **Progress Tracking**: Visual indicators show completed steps
- **Interactive Parameters**: Adjust analysis parameters in real-time
- **Data Persistence**: Results persist between tabs within a session
- **Export Options**: Download processed data, reports, and session info
- **Error Handling**: Graceful error messages with troubleshooting hints

## Development Guide

### Adding a New Omics Type

1. Create directory structure:
```bash
mkdir -p R/shiny/ui/new_omics
mkdir -p R/shiny/server/new_omics
mkdir -p R/shiny/modules/new_omics
mkdir -p R/new_omics/workflow
mkdir -p R/new_omics/utils
```

2. Create UI module: `R/shiny/ui/new_omics/workflow_ui.R`
```r
newOmicsWorkflowUi <- function(id) {
  ns <- NS(id)
  # Define UI for workflow tabs
}
```

3. Create server module: `R/shiny/server/new_omics/workflow_server.R`
```r
newOmicsWorkflowServer <- function(id, workflow_state) {
  moduleServer(id, function(input, output, session) {
    # Define server logic
  })
}
```

4. Update `app.R` to include the new omics type in the selection interface

### Module Communication

Modules communicate through:
- `workflow_state`: Shared reactive values across omics types
- Return values from server modules
- Global reactive values in the main app

### Git Subtree Management

The modular structure supports git subtree workflows:

```bash
# Add a subtree for proteomics
git subtree add --prefix=R/proteomics https://github.com/org/proteomics-workflow.git main

# Pull updates
git subtree pull --prefix=R/proteomics https://github.com/org/proteomics-workflow.git main

# Push changes
git subtree push --prefix=R/proteomics https://github.com/org/proteomics-workflow.git main
```

## Integration with Existing Functions

The app reuses existing MultiScholaR functions by:
1. Sourcing workflow functions from the appropriate directories
2. Wrapping functions in reactive contexts
3. Handling progress updates and error messages
4. Managing intermediate results in reactive values

## Troubleshooting

### Common Issues

1. **"Module not found" errors**
   - Ensure all required `.R` files are in the correct directories
   - Check that `sourceDir()` is finding your modules

2. **Data not persisting between tabs**
   - Verify reactive values are properly initialized
   - Check that values are being assigned correctly

3. **Design matrix applet not launching**
   - Ensure `RunApplet()` function is available
   - Check that required objects exist in the environment

## Future Enhancements

- [ ] Save/load analysis sessions
- [ ] Batch processing capabilities
- [ ] Advanced visualization options
- [ ] Cloud deployment support
- [ ] API endpoints for programmatic access
- [ ] Real-time collaboration features