# ProteomeScholaR <img src="https://img.shields.io/badge/Version-0.8-green?style=for-the-badge" alt="Version 0.8">

## Quick Start

### 1. Required Software

<a href="https://cran.r-project.org/bin/windows/base/" target="_blank">
    <img src="https://img.shields.io/badge/Download-R_(Windows)-276DC3?style=for-the-badge&logo=r" alt="Download R for Windows">
</a>
<a href="https://cran.r-project.org/bin/macosx/" target="_blank">
    <img src="https://img.shields.io/badge/Download-R_(macOS)-276DC3?style=for-the-badge&logo=r" alt="Download R for macOS">
</a>
<a href="https://posit.co/download/rstudio-desktop/" target="_blank">
    <img src="https://img.shields.io/badge/Download-RStudio_Desktop-75AADB?style=for-the-badge&logo=rstudio" alt="Download RStudio">
</a>

### 2. Setup Script

<a href="https://raw.githubusercontent.com/APAF-bioinformatics/ProteomeScholaR/dev-jr/project_setup.R" download="project_setup.R">
    <img src="https://img.shields.io/badge/Download-Setup_Script-blue?style=for-the-badge&logo=r" alt="Download Setup Script">
</a>

### 3. Tutorial Data (Optional)

<a href="https://drive.google.com/file/d/1qeS2X1uA_Y7HFGMVru0_tAbEQmVjsdlD/view?usp=drive_link" target="_blank">
    <img src="https://img.shields.io/badge/Download-Tutorial_Data-orange?style=for-the-badge&logo=google-drive" alt="Download Tutorial Data">
</a>

This tutorial dataset contains example data from *Klebsiella variicola*, including:
- Example DIA-NN search results
- Example organism FASTA file
- NCBI annotation protein data searched against UniProt and UniParc databases

The data is derived from the publicly available dataset published in [Mu, Klare, Baines, Pang et al., (2023) Nature Communications](https://www.nature.com/articles/s41467-023-37200-w), which performed integrative omics analysis on sepsis-causing bacteria.

## Setup Instructions

1. Install RStudio Desktop and R if you haven't already (use button above)
2. Download the setup script using the button above (right click + save as)
3. Open the downloaded file in RStudio
4. Change the project name at the top of the script (optional: specify custom directory)
5. Run the entire script (Ctrl+A then Ctrl+Enter)
6. A new RStudio project will open automatically with all required files and structure

## What Gets Set Up

- Complete directory structure for proteomics analysis
- Latest version of DIA workflow
- Configuration files
- R project file

## Using the Workflow

1. Run the DIA_workflow.rmd file (this should automatically have opened if you did the above correctly!)
2. Please copy your organism .fasta to the data/UniProt subdirectory
3. Please copy your searched data to the data/proteomics subdirectory
4. Proceed chunk by chunk
5. Use the embedded Shiny app to define your:
   - Experimental design
   - Contrasts
   - Linear model
6. Find all results in the summary_results folder
7. Run DIA_report_1.0.rmd to generate a shareable report (current not finished)

## Contributors 
* Ignatius Pang (ignatius.pang@mq.edu.au) 
* Will Klare (william.klare@mq.edu.au) 

## Version Information

v0.8 provides a feature-complete workflow for processing searched DIA-NN proteomics data.

## Roadmap


```mermaid
graph TD
A[Current Version 0.8] --> B[Version 0.9]
B --> C[Version 1.0]
B --> D[Interactive Volcano Plots]
B --> E[Enhanced Dynamic Report Generation]
C --> F[Unsupervised Analysis & Interactive Plots]
C --> G[Enhanced Enrichment Analysis]
G --> H[Greater Control Over Enrichment Platforms]
G --> I[Improved Data Visualisation for Enrichments]
style A fill:#90EE90,color:#000000
style B fill:#FFE4B5,color:#000000
style C fill:#FFB6C1,color:#000000
```

## Need Help?

If you encounter any issues:
1. Check the [Issues](https://github.com/APAF-bioinformatics/ProteomeScholaR/issues) page
2. Contact the contributors
3. Submit a new issue

## We very much encourage you to report any and all bugs using the issues page - we are commited to making this tool as robust as possible and we would love your help to do this :)

Enjoy! 🧬🔬