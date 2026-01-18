# MultiScholaR <img src="https://img.shields.io/badge/Version-0.4.1-orange?style=for-the-badge" alt="Version 0.4.1">

![Banner](inst/shiny/www/MultiScholaR.png "MultiScholaR Banner")

> **⚠️ Disclaimer:** `MultiScholaR` is under active development. The **end-to-end Multiomics GUI** currently supports Proteomics, Metabolomics, and Lipidomics workflows.

## Overview

**MultiScholaR** bridges the gap between complex multi-omics data and accessible, reproducible analysis. Unlike "black box" tools, it provides transparent, well-documented templates alongside a production-grade GUI, empowering researchers to understand the _how_ and _why_ behind their results.

**Key Features:**

- **Modular Architecture:** Flexible object-oriented design for easy tool integration.
- **Reproducibility:** Standardized workflows with automated, publication-ready reporting.
- **Advanced Quality Control:** Stringent measures including FDR thresholds, intensity filtering, and automated batch effect correction (RUV-IIIc).
- **Multi-Omics Integration:** Seamlessly integrate Proteomics, Metabolomics, and Lipidomics data.

## Analysis Workflow

```mermaid
graph TD
    %% Proteomics
    subgraph Proteomics ["Proteomics Workflows"]
        DIA["DIA-NN"] --> Pep["Peptide Rollup (IQ)"]
        Pep --> ProtQC
        TMT["TMT / PD"] --> ProtQC
        LFQ["LFQ / FragPipe"] --> ProtQC

        ProtQC["QC & Normalization<br/>(RUV-IIIc / cyclic-loess)"] --> ProtDE["Differential Expression"]
        ProtDE --> ProtEnrich["Enrichment Analysis"]
        ProtEnrich --> ProtReport["Proteomics Report"]
    end

    %% Metabolomics
    subgraph Metabolomics ["Metabolomics Workflow"]
        MetInput["LC-MS / GC-MS"] --> MetQC["QC & Normalization<br/>(RUV-IIIc / cyclic-loess)"]
        MetQC --> MetDE["Differential Expression"]
        MetDE --> MetEnrich["Enrichment Analysis"]
        MetEnrich --> MetReport["Metabolomics Report"]
    end

    %% Lipidomics
    subgraph Lipidomics ["Lipidomics Workflow"]
        LipInput["LipidSearch / MS-Dial"] --> LipQC["QC & Normalization<br/>(RUV-IIIc / cyclic-loess)"]
        LipQC --> LipDE["Differential Expression"]
        LipDE --> LipEnrich["Enrichment Analysis"]
        LipEnrich --> LipReport["Lipidomics Report"]
    end

    %% Multi-Omics Integration (Union)
    ProtDE -.-> Integration["Multi-Omics Integration<br/>(MOFA2 / MixOmics)"]
    MetDE -.-> Integration
    LipDE -.-> Integration

    %% Styling
    style Integration fill:#f3e5f5,stroke:#4a148c,stroke-width:2px
    style Pep fill:#fff3e0,stroke:#f57c00

    style ProtReport fill:#e8f5e9,stroke:#1b5e20
    style MetReport fill:#e8f5e9,stroke:#1b5e20
    style LipReport fill:#e8f5e9,stroke:#1b5e20
```

## Getting Started

### 1. The GUI Experience (Recommended)

For a **code-free, interactive experience**, use the standalone launcher.

[![Download Launcher](https://img.shields.io/badge/Download-MultiScholaR_Launcher-blue?style=for-the-badge&logo=github)](https://github.com/APAF-bioinformatics/MultiScholaR-launcher)

1.  Download the appropriate launcher for your OS (Windows/Mac).
2.  Double-click to install/update and launch the application.
3.  Follow the guided wizard for Import, QC, and Analysis.

### 2. R Markdown Workflows (Advanced)

For **programmatic, reproducible analyses**, typically for advanced users or batch processing.

1.  **Install Prerequisites**: R (4.4.3+), RStudio, and Rtools (Windows).
2.  **Setup Project**: Download and run the setup script to generate a structured project with all necessary templates.

    <a href="https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/project_setup.R" download="project_setup.R">
        <img src="https://img.shields.io/badge/Download-Setup_Script-276DC3?style=for-the-badge&logo=r" alt="Download Setup Script">
    </a>

3.  **Run Workflows**: Open the generated `.Rmd` templates (Proteomics, Metabolomics, etc.) in `inst/workbooks/` and execute chunks sequentially.

> **Need Help?** Visit the [Multiomics Masterclass](https://zenodo.org/records/15573343) for a comprehensive guide and tutorial data.

## Supported Workflows

| Omics Type       | Methods Supported            | Key Features                                              |
| :--------------- | :--------------------------- | :-------------------------------------------------------- |
| **Proteomics**   | DIA-NN, TMT-PD, LFQ-Fragpipe | Peptide rollup (IQ), Auto-RUV-IIIc, Mixed Species support |
| **Metabolomics** | LC-MS, GC-MS (MS-Dial)       | Intensity QC, Internal Standard Normalization             |
| **Lipidomics**   | MS-Dial, LipidSearch         | Lipid-class specific analysis, QC & filtering             |
| **Integration**  | MOFA2, MixOmics              | Multi-view factor analysis, Cross-omics enrichment        |

## Contributors

- Ignatius Pang (ignatius.pang@mq.edu.au)
- Will Klare (william.klare@mq.edu.au)

## Attributions

Derived from [ProteomeRiver](https://bitbucket.org/cmri-bioinformatics/proteomeriver/src/main/) (LGPLv3). Significant modifications by APAF-bioinformatics.
AI assistance used responsibly in development.
