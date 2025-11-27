# Function Files Index

This document provides an overview of all `func_*.R` files created for golem-ification of the MultiScholaR package.

## Overview

The `func_*.R` files serve as structured containers for functions organized by:
- **Omics type**: proteomics, metabolomics, lipidomics, multiomics
- **Analysis stage**: import, QC, normalization, DE, enrichment
- **General utilities**: plotting, file management, helpers, annotation, enrichment, design, S4 objects

These files are designed to enable frictionless manual extraction of functions from monolithic files.

## File Structure

### Proteomics Functions

1. **`func_prot_import.R`** - Proteomics data import functions
   - DIA-NN, TMT-PD, LFQ-Fragpipe, MaxQuant, Spectronaut import
   - Format detection and conversion functions

2. **`func_prot_qc_peptide.R`** - Peptide-level QC filtering
   - Intensity filtering, missing value filtering, replicate filtering
   - Q-value and proteotypic peptide filtering

3. **`func_prot_qc.R`** - Protein-level QC filtering
   - Intensity filtering, missing value filtering, replicate filtering
   - Protein cleanup operations

4. **`func_prot_rollup.R`** - Peptide-to-protein rollup
   - Precursor-to-peptide rollup
   - Counting functions for peptides and proteins

5. **`func_prot_norm.R`** - Protein normalization
   - RUV-III-C normalization, negative control selection
   - Normalization parameter optimization

6. **`func_pept_norm.R`** - Peptide normalization
   - RUV-III-C normalization for peptides
   - Log2 transformation

7. **`func_peptide_qc_imputation.R`** - Missing value imputation
   - Limpa-based imputation for peptides and proteins
   - Imputation validation functions

8. **`func_prot_de.R`** - Protein differential expression
   - Limma-based DE analysis
   - Result formatting and visualization

### Metabolomics Functions

9. **`func_metab_import.R`** - Metabolomics data import
   - MetaboliteAssayData S4 object creation
   - Data extraction functions

10. **`func_metab_qc.R`** - Metabolomics QC filtering
    - Intensity filtering, missing value analysis
    - CV calculations, internal standard metrics

11. **`func_metab_norm.R`** - Metabolomics normalization
    - Log transformation, between-sample normalization
    - Negative control selection

12. **`func_metab_de.R`** - Metabolomics differential abundance
    - Limma-based DE analysis for metabolites
    - Result formatting

### Lipidomics Functions (Placeholders)

13. **`func_lipid_import.R`** - Placeholder for future lipidomics import
14. **`func_lipid_qc.R`** - Placeholder for future lipidomics QC
15. **`func_lipid_norm.R`** - Placeholder for future lipidomics normalization
16. **`func_lipid_de.R`** - Placeholder for future lipidomics DE

### Multiomics Functions

17. **`func_multiomics_mofa.R`** - MOFA integration
    - MOFA model building and analysis
    - Visualization functions

18. **`func_multiomics_enrich.R`** - Multiomics enrichment
    - StringDB enrichment
    - Metabolomics pathway enrichment
    - Rank-based enrichment methods

### S4 Object Files

19. **`func_general_s4_generics.R`** - S4 generic function definitions
    - All setGeneric() definitions (replaces allGenerics.R)
    - Must be loaded early before setMethod() calls

20. **`func_prot_s4_objects.R`** - Proteomics S4 classes and methods
    - ProteinQuantitativeData, PeptideQuantitativeData classes
    - Proteomics-specific S4 methods

21. **`func_metab_s4_objects.R`** - Metabolomics S4 classes and methods
    - MetaboliteAssayData class
    - Metabolomics-specific S4 methods

22. **`func_lipid_s4_objects.R`** - Lipidomics S4 classes (placeholder)

23. **`func_general_s4_objects.R`** - Shared S4 classes
    - FilteringProgress classes
    - Shared S4 utilities

### General Functions

24. **`func_general_plotting.R`** - General plotting utilities
    - Volcano plots, PCA plots, RLE plots
    - Correlation plots, enrichment plots
    - Heatmaps, interactive plots
    - Plot saving and theme functions

25. **`func_general_filemgmt.R`** - File management
    - Directory setup and management
    - Configuration file handling
    - Project path management

26. **`func_general_helpers.R`** - Helper utilities
    - Parameter checking, parsing functions
    - Data manipulation, correlation calculations
    - Result extraction, filtering progress tracking

27. **`func_general_annotation.R`** - Annotation functions
    - UniProt annotation retrieval
    - Protein ID cleaning and selection
    - FASTA processing, phosphoproteomics annotation

28. **`func_general_enrichment.R`** - General enrichment analysis
    - GO, KEGG, Reactome enrichment
    - GSEA, enricher functions
    - Enrichment visualization and processing

29. **`func_general_design.R`** - Design matrix functions
    - Design matrix cleaning and manipulation
    - Category mapping and column creation

## Extraction Guide

### Step-by-Step Extraction Process

1. **Identify Function Location**
   - Each `func_*.R` file contains TODO comments listing functions to extract
   - Each function entry includes:
     - Function name
     - Current file location
     - Line numbers (approximate)
     - Brief description

2. **Extract Function**
   - Open the source file listed in the TODO comment
   - Locate the function using the provided line numbers
   - Copy the entire function (including roxygen documentation)
   - Paste into the appropriate `func_*.R` file

3. **Update Dependencies**
   - Check function dependencies
   - Ensure required packages are imported
   - Update function calls if needed

4. **Test Function**
   - After extraction, test the function in isolation
   - Verify it works with extracted dependencies
   - Update documentation if needed

5. **Update Source File**
   - Remove extracted function from source file
   - Add `source()` call or ensure function is loaded from `func_*.R`
   - Update any direct calls to use the new location

### Best Practices

- **Extract one function at a time** to avoid breaking dependencies
- **Test after each extraction** to catch issues early
- **Maintain function signatures** - don't change parameter names or order
- **Preserve roxygen documentation** - keep all `@param`, `@return`, `@export` tags
- **Check for S4 methods** - ensure `setMethod()` calls are properly extracted
- **Update NAMESPACE** - ensure exported functions remain exported

### Dependency Management

Functions may depend on:
- Other functions in the same `func_*.R` file
- Functions in other `func_*.R` files
- Functions from external packages
- S4 class definitions

When extracting:
1. Identify all dependencies
2. Extract dependencies first (or ensure they're available)
3. Update import statements in the `func_*.R` file
4. Test that all dependencies resolve correctly

## File Naming Conventions

- **`func_[omic]_[stage].R`** - Omics-specific functions
  - Example: `func_prot_qc_peptide.R`
- **`func_general_[category].R`** - General utility functions
  - Example: `func_general_plotting.R`

## Classification Reference

See `FUNC_CLASSIFICATION.md` for detailed function-to-file mappings.

## Next Steps

1. Begin extracting functions from monolithic files
2. Start with low-dependency functions (e.g., plotting utilities)
3. Progress to higher-level functions (e.g., DE analysis)
4. Test each extraction thoroughly
5. Update module files to use extracted functions

## Notes

- All `func_*.R` files are currently placeholders with TODO comments
- Functions should be extracted manually as described above
- The structure is designed to minimize friction during extraction
- Each file includes comprehensive function lists with source locations

