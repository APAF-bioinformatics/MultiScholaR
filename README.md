# MultiScholaR <img src="https://img.shields.io/badge/Version-0.1-orange?style=for-the-badge" alt="Version 0.1">

>**⚠️ Disclaimer:** `MultiScholaR` is currently under active development and is **not yet ready for general use**. The core package structure and initial proteomics workflow are being established. Functionality described in the roadmap is planned for future releases.

## Overview

Modern multi-omics datasets require substantial programming expertise and practical knowledge to analyse. Analysis is skill-gated, making it challenging for many researchers to apply modern statistical best practices. MultiScholaR is a package in R that addresses this challenge by providing a novel pipeline aiming to enable users to perform comprehensive multi-omics analyses, including single-omic analyses (e.g. transcriptomics, proteomics, phosphoproteomics and metabolomics datasets) and integrative multi-omics analysis. Through well-documented workflow templates, researchers can systematically apply best-practice to all steps of integrative multi-omics analysis.

MultiScholaR implements stringent quality control measures for multi-omics analysis by incorporating criteria such as false discovery rate thresholds, filtering criterias, and missing value limitations across samples. It integrates several sophisticated analytical tools: the IQ tool for peptide-to-protein quantitative data summarization, RUVIII-C for removing unwanted variation, and edgeR and/or limma for sample normalization and linear modelling. Pathway analysis can be performed either using user-supplied annotations via clusterProfiler or through automated analysis with gProfiler2. Multivariate and integrative multi-omics analyses are implemented using MOFA and MixOmics.

Structured on modular, object-oriented components, MultiScholaR's architecture facilitates easy integration of new tools as they emerge. The inclusion of comprehensive, documented workbooks that guide users through each analytical step, facilitates reproducibility and enabling public sharing of analyses.

By streamlining complex multi-omics analyses, MultiScholaR makes advanced analytical techniques accessible to researchers across all levels of programming expertise.

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
        D3(Factor Analysis Methods - MOFA+Visualisation)
    end
     subgraph "v0.5"
        E
        E1(Advanced Integration Visualization)
        E2(Pathway/Network Integration)
        E3(StringDB Integration)
    end

    style A fill:#ADD8E6,color:#000000
    style B fill:#FFB6C1,color:#000000
    style C fill:#90EE90,color:#000000
    style D fill:#FFDAB9,color:#000000
    style E fill:#E6E6FA,color:#000000
```