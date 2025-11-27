# Reorganization Summary

This document summarizes the changes made to address organizational concerns.

## Changes Made

### 1. Omic-Specific S4 Objects

**Created separate files for omic-specific S4 classes:**

- **`func_prot_s4_objects.R`** - Proteomics S4 classes and methods
  - `ProteinQuantitativeData` class
  - `PeptideQuantitativeData` class
  - All proteomics-specific S4 methods

- **`func_metab_s4_objects.R`** - Metabolomics S4 classes and methods
  - `MetaboliteAssayData` class
  - `MetabolomicsDifferentialAbundanceResults` class
  - All metabolomics-specific S4 methods

- **`func_lipid_s4_objects.R`** - Placeholder for future lipidomics S4 classes

- **`func_general_s4_objects.R`** - Updated to contain only shared S4 classes
  - `FilteringProgress` class (shared)
  - `FilteringProgressMetabolomics` class
  - Shared S4 utilities

### 2. S4 Generics File

**Created `func_general_s4_generics.R`:**
- Contains all `setGeneric()` definitions from `R/allGenerics.R`
- Must be loaded early before any `setMethod()` calls
- Replaces the need for a separate `allGenerics.R` file (or can coexist)

### 3. File Naming Simplification

**Renamed files for clarity:**

- `func_prot_qc_protein.R` → **`func_prot_qc.R`**
  - Rationale: "protein" suffix was redundant since it's already in the proteomics context
  - Peptide QC remains in `func_prot_qc_peptide.R` to distinguish from protein QC

- `func_prot_qc_rollup.R` → **`func_prot_rollup.R`**
  - Rationale: Rollup is a distinct operation, not strictly QC, so removing "qc" from name

### 4. Peptide Functions Clarification

**Clarified in `func_prot_qc_peptide.R`:**
- Added note explaining that peptides are part of proteomics workflow
- The "prot" prefix is correct because peptides are within the proteomics context
- Peptide-specific functions are correctly placed in proteomics files

## File Structure Summary

### Proteomics Functions (8 files)
1. `func_prot_import.R` - Data import
2. `func_prot_qc_peptide.R` - Peptide QC (peptides are part of proteomics)
3. `func_prot_qc.R` - Protein QC (renamed from func_prot_qc_protein.R)
4. `func_prot_rollup.R` - Peptide-to-protein rollup (renamed from func_prot_qc_rollup.R)
5. `func_prot_norm.R` - Protein normalization
6. `func_pept_norm.R` - Peptide normalization
7. `func_peptide_qc_imputation.R` - Missing value imputation
8. `func_prot_de.R` - Differential expression

### S4 Object Files (5 files)
1. `func_general_s4_generics.R` - All S4 generic definitions (NEW)
2. `func_prot_s4_objects.R` - Proteomics S4 classes and methods (NEW)
3. `func_metab_s4_objects.R` - Metabolomics S4 classes and methods (NEW)
4. `func_lipid_s4_objects.R` - Lipidomics S4 classes (placeholder, NEW)
5. `func_general_s4_objects.R` - Shared S4 classes (updated)

### Other Files
- Metabolomics: 4 files (unchanged)
- Lipidomics: 4 placeholder files (unchanged)
- Multiomics: 2 files (unchanged)
- General: 7 files (unchanged)

## Rationale

1. **Omic-specific S4 objects**: Separates concerns - each omics type has its own S4 classes and methods, making it easier to maintain and understand.

2. **S4 generics file**: Centralizes all generic definitions, ensuring they're loaded before methods. This is critical for S4 method dispatch.

3. **Simplified naming**: 
   - `func_prot_qc.R` is clearer than `func_prot_qc_protein.R` (redundant "protein")
   - `func_prot_rollup.R` better reflects that rollup is a distinct operation, not just QC

4. **Peptide functions**: Peptides are correctly placed in proteomics files since they're part of the proteomics workflow. The distinction between peptide and protein QC is maintained through separate files.

## Next Steps

1. Extract functions according to the updated file structure
2. Ensure `func_general_s4_generics.R` is loaded early in package initialization
3. Update any references to renamed files
4. Test that S4 method dispatch works correctly after reorganization

