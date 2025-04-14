# MultiScholaR <img src="https://img.shields.io/badge/Version-0.1-orange?style=for-the-badge" alt="Version 0.1">

>**⚠️ Disclaimer:** `MultiScholaR` is currently under active development and is **not yet ready for general use**. The core package structure and initial proteomics workflow are being established. Functionality described in the roadmap is planned for future releases.

## Overview

`MultiScholaR` aims to provide a standardized, reproducible, and user-friendly R environment for the analysis of various omics data types (proteomics, metabolomics, transcriptomics) and their integration. It leverages S4 objects and established Bioconductor principles for robust data handling and analysis.

## Quick Start & Setup

*(Currently unavailable - See Roadmap)*

Details on installation and setup will be provided as the project matures towards stable releases.

## Contributors
* Ignatius Pang (ignatius.pang@mq.edu.au)
* Will Klare (william.klare@mq.edu.au)

## Development Roadmap

```mermaid
graph TD
    A[v0.1: Foundation & Proteomics Core] --> B[v0.2: Metabolomics Module];
    B --> C[v0.3: Transcriptomics Module];
    C --> D[v0.4: Basic Multiomics Integration];
    D --> E[v0.5: Advanced Multiomics Analysis & Visualization];

    subgraph "v0.1"
        A
        A1(✓ Refactor Core Package)
        A2(✓ Solidify Proteomics Workflow)
        A3(✓ Establish S4 Structure)
    end
    subgraph "v0.2"
        B
        B1(Metabolomics Data Input)
        B2(Metabolomics QC & Normalization)
        B3(Basic Metabolomics Analysis)
    end
    subgraph "v0.3"
        C
        C1(Transcriptomics Data Input)
        C2(Transcriptomics QC & Normalization)
        C3(Basic Transcriptomics Analysis - DEGs)
    end
     subgraph "v0.4"
        D
        D1(MultiAssayExperiment Integration)
        D2(Cross-Omics Data Linking)
        D3(Initial Integrated Visualization)
    end
     subgraph "v0.5"
        E
        E1(Factor Analysis Methods - MOFA+)
        E2(Advanced Integration Visualization)
        E3(Pathway/Network Integration)
    end

    style A fill:#ADD8E6,color:#000000
    style B fill:#FFB6C1,color:#000000
    style C fill:#90EE90,color:#000000
    style D fill:#FFDAB9,color:#000000
    style E fill:#E6E6FA,color:#000000
```