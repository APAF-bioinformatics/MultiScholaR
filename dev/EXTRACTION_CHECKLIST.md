# Function Extraction Checklist

**Purpose**: Guide for extracting functions from monolithic R files to organized `func_*.R` files.

**Tool**: `dev/extract_functions.R` - Run with `source("dev/extract_functions.R")` in R.

---

## Quick Start for New Agent Session

```r
# 1. Load the extraction tool
source("dev/extract_functions.R")

# 2. Preview what will be extracted (always do this first!)
extract_all(dry_run = TRUE)

# 3. Execute extraction
extract_all(dry_run = FALSE)

# 4. Verify package still builds
devtools::document()

# 5. Check logs
# - dev/extraction_log.md (full log)
# - dev/extraction_log_errors.md (functions not found)
```

---

## Session Summary: 2025-11-27 (Session 4 - Phase 3)

**Phase 3 Accomplishments:**
- **Verified all standalone functions** extracted to func_*.R or mod_*.R files
- **Extracted 25+ additional functions** missed in earlier phases:
  - `calcHtSize`, `updateProteinFiltering`, `generateMetaboliteFilteringPlots`
  - `plotStringDbEnrichmentResults`, `runKeggEnrichment`, `runReactomeEnrichment`
  - `writeInteractiveVolcanoPlotProteomics*`, `getDataMatrix`, `plotPca`
  - `PeptideQuantitativeDataDiann`, `compareTwoPeptideDataObjects`, etc.
- **Fixed extraction tool** comment detection bug (was falsely matching commented code)
- **Moved RunApplet** to `launch_app.R` (consolidated app launchers)
- **DELETED 11 redundant source files** (all functions extracted, no S4 content):
  - `annotation.R`, `enrichment_functions.R`, `file_management.R`
  - `multiomics_functions_MOFA.R`, `de_proteins_functions.R`
  - `string_enrichment_functions_refactored.R`, `de_analysis_function_wrapper.R`
  - `multiomics_enrichment_functions.R`, `shiny_applets.R`
  - `get_best_accession_helper.R`, `phosphoproteomics_helper.R`
- **Package builds successfully** with `devtools::document()`
- **Total functions extracted: 340+**

**Previous Session (2025-11-27 - Session 3 - Phase 2):**
- Created 3 new func_*.R files: `func_phospho_annotation.R`, `func_prot_limpa.R`, `func_pept_limpa.R`
- Extracted **44 additional functions** from 4 overlooked source files

**Previous Session (2025-11-27 - Session 2):**
- Deleted 7 redundant/empty placeholder files
- Extracted 63 setGeneric() definitions from `allGenerics.R` to `func_general_s4_generics.R`

**Original Session (2025-11-27 - Session 1):**
- Fixed 49 original extraction errors by correcting source file mappings
- Extracted 195 total functions (36 in Phase 1 + 159 in Phases 2-5)
- Audited and configured 5 major source files

---

## Current File Structure Summary

**Populated func_*.R files:** 22 files, 340+ functions
**S4 Generics:** 63 setGeneric() definitions in `func_general_s4_generics.R`
**Deleted source files:** 11 (all standalone functions extracted)
**Remaining S4 source files:** 13 (contain S4 classes/methods that must stay)

### Files Deleted in Phase 3

These files contained only standalone functions (now in func_*.R):
- `annotation.R`, `enrichment_functions.R`, `file_management.R`
- `multiomics_functions_MOFA.R`, `de_proteins_functions.R`
- `string_enrichment_functions_refactored.R`, `de_analysis_function_wrapper.R`
- `multiomics_enrichment_functions.R`, `shiny_applets.R`
- `get_best_accession_helper.R`, `phosphoproteomics_helper.R`

### Files Kept (S4 Content)

These files contain S4 classes or methods and must be preserved:
- `peptideVsSamplesS4Objects.R` - PeptideQuantitativeData class + 23 methods
- `proteinVsSamplesS4Objects.R` - ProteinQuantitativeData methods
- `metaboliteVsSamplesS4Objects.R` - MetaboliteAssayData class + 16 methods
- `functional_enrichment.R` - EnrichmentResults classes + methods
- `limpa_functions.R` - S4 methods for limpa imputation
- `allGenerics.R` - setGeneric() definitions (needed for load order)
- `helper_functions.R`, `qc_and_rollup.R`, `QC_visualisation.R` - FilteringProgress classes
- `metabolite_*.R` files - S4 methods for metabolomics
- `protein_de_analysis_wrapper.R` - S4 methods for DE analysis

---

## Source Files Status (After Phase 3)

| File | Status | Notes |
|------|--------|-------|
| `helper_functions.R` | KEEP (S4) | Contains DirectoryManager setClass |
| `de_proteins_functions.R` | **DELETED** | All functions extracted |
| `multiomics_enrichment_functions.R` | **DELETED** | All functions extracted |
| `metaboliteVsSamplesS4Objects.R` | KEEP (S4) | MetaboliteAssayData class + 16 methods |
| `qc_and_rollup.R` | KEEP (S4) | Contains FilteringProgress setClass |
| `peptideVsSamplesS4Objects.R` | KEEP (S4) | PeptideQuantitativeData class + 23 methods |
| `enrichment_functions.R` | **DELETED** | All functions extracted |
| `QC_visualisation.R` | KEEP (S4) | Contains FilteringProgressMetabolomics setClass |
| `limpa_functions.R` | KEEP (S4) | S4 methods for limpa imputation |
| `functional_enrichment.R` | KEEP (S4) | EnrichmentResults classes + methods |
| `de_analysis_function_wrapper.R` | **DELETED** | All functions extracted |
| `protein_de_analysis_wrapper.R` | KEEP (S4) | Contains S4 methods |
| `get_best_accession_helper.R` | **DELETED** | All functions extracted |
| `annotation.R` | **DELETED** | All functions extracted |
| `file_management.R` | **DELETED** | All functions extracted |
| `phosphoproteomics_helper.R` | **DELETED** | All functions extracted |
| `metabolite_qc.R` | KEEP (S4) | S4 methods for metabolite QC |
| `string_enrichment_functions_refactored.R` | **DELETED** | All functions extracted |
| `allGenerics.R` | KEEP (S4) | setGeneric() definitions (load order) |
| `metabolite_normalization.R` | KEEP (S4) | S4 methods for normalization |
| `github_managment.R` | NOT APPLICABLE | Git utilities (keep as-is) |
| `metabolites_de_analysis_wrapper.R` | KEEP (S4) | S4 class + methods |
| `metabolite_de_analysis_wrapper.R` | KEEP (S4) | S4 class + methods |
| `proteinVsSamplesS4Objects.R` | KEEP (S4) | S4 method for protein imputation |
| `multiomics_functions_MOFA.R` | **DELETED** | All functions extracted |
| `shiny_applets.R` | **DELETED** | RunApplet moved to launch_app.R |

---

## Target func_*.R Files

### Proteomics Files

| Target File | Status | Source Files | Notes |
|------------|--------|--------------|-------|
| `func_prot_import.R` | [x] COMPLETE | mod_prot_import.R, file_management.R | 9 functions |
| `func_prot_qc_peptide.R` | [x] COMPLETE | peptideVsSamplesS4Objects.R, qc_and_rollup.R, helper_functions.R | All helpers extracted |
| `func_prot_qc.R` | [x] COMPLETE | proteinVsSamplesS4Objects.R, qc_and_rollup.R, de_proteins_functions.R | All functions extracted |
| `func_prot_rollup.R` | [x] COMPLETE | peptideVsSamplesS4Objects.R, qc_and_rollup.R | 10+ functions |
| `func_prot_norm.R` | [x] COMPLETE | proteinVsSamplesS4Objects.R, de_proteins_functions.R, qc_and_rollup.R | All functions extracted |
| `func_pept_norm.R` | [x] COMPLETE | peptideVsSamplesS4Objects.R, qc_and_rollup.R | log2Transformation extracted |
| `func_peptide_qc_imputation.R` | [x] COMPLETE | peptideVsSamplesS4Objects.R, limpa_functions.R, qc_and_rollup.R, helper_functions.R | All functions |
| `func_prot_de.R` | [x] COMPLETE | protein_de_analysis_wrapper.R, de_proteins_functions.R | 20+ functions |
| `func_prot_annotation.R` | [x] COMPLETE | annotation.R, de_proteins_functions.R, get_best_accession_helper.R | 29+ functions |
| `func_prot_limpa.R` | [x] COMPLETE | limpa_functions.R | 2 functions (QC plots, DE conversion) |
| `func_pept_limpa.R` | [ ] PLACEHOLDER | - | For peptide limpa helpers (future) |
| `func_phospho_annotation.R` | [x] COMPLETE | phosphoproteomics_helper.R | 22 phospho annotation functions |
| `func_prot_s4_objects.R` | [ ] PLACEHOLDER | peptideVsSamplesS4Objects.R, proteinVsSamplesS4Objects.R | S4 classes (future) |

### Metabolomics Files

| Target File | Status | Source Files | Notes |
|------------|--------|--------------|-------|
| `func_metab_import.R` | [x] COMPLETE | metaboliteVsSamplesS4Objects.R, QC_visualisation.R | 2 functions |
| `func_metab_qc.R` | [x] COMPLETE | metabolite_qc.R, QC_visualisation.R | 15 functions |
| `func_metab_norm.R` | [ ] PLACEHOLDER | metabolite_normalization.R | Target for future extraction |
| `func_metab_de.R` | [x] COMPLETE | metaboliteVsSamplesS4Objects.R | 3 functions |
| `func_metab_s4_objects.R` | [ ] PLACEHOLDER | metaboliteVsSamplesS4Objects.R | S4 classes (future) |

### Multiomics Files

| Target File | Status | Source Files | Notes |
|------------|--------|--------------|-------|
| `func_multiomics_mofa.R` | [x] COMPLETE | multiomics_functions_MOFA.R | 1 function |
| `func_multiomics_enrich.R` | [x] COMPLETE | multiomics_enrichment_functions.R | 14 functions |

### General Files

| Target File | Status | Source Files | Notes |
|------------|--------|--------------|-------|
| `func_general_plotting.R` | [x] COMPLETE | qc_and_rollup.R, de_proteins_functions.R, helper_functions.R | 30+ plotting functions |
| `func_general_filemgmt.R` | [x] COMPLETE | file_management.R, qc_and_rollup.R, helper_functions.R | 25+ file/config functions |
| `func_general_helpers.R` | [x] COMPLETE | qc_and_rollup.R, helper_functions.R | 15+ utility functions |
| `func_general_enrichment.R` | [x] COMPLETE | enrichment_functions.R, functional_enrichment.R | 52 functions |
| `func_general_design.R` | [x] COMPLETE | proteinVsSamplesS4Objects.R, peptideVsSamplesS4Objects.R, qc_and_rollup.R | 5 functions |
| `func_general_s4_generics.R` | [x] COMPLETE | allGenerics.R | 63 setGeneric() definitions |
| `func_general_s4_objects.R` | [ ] PLACEHOLDER | qc_and_rollup.R | For shared S4 classes (future) |

### Deleted Files (Redundant Placeholders)

The following files were deleted on 2025-11-27 as redundant or premature:
- `func_prot_qc_protein.R` - Functions in func_prot_qc.R
- `func_prot_qc_rollup.R` - Functions in func_prot_rollup.R
- `func_lipid_import.R` - Create when lipidomics implemented
- `func_lipid_qc.R` - Create when lipidomics implemented
- `func_lipid_norm.R` - Create when lipidomics implemented
- `func_lipid_de.R` - Create when lipidomics implemented
- `func_lipid_s4_objects.R` - Create when lipidomics implemented

---

## How to Add New Extractions

### Step 1: Find Functions to Extract

Use grep to find all function definitions in a source file:

```bash
# Find all function definitions
grep -n "^[a-zA-Z_][a-zA-Z0-9_]* *<- *function" R/SOURCE_FILE.R

# Find S4 method definitions
grep -n "setMethod" R/SOURCE_FILE.R

# Find S4 generic definitions  
grep -n "setGeneric" R/SOURCE_FILE.R

# Find S4 class definitions
grep -n "setClass" R/SOURCE_FILE.R
```

### Step 2: Add Configuration to extract_functions.R

Open `dev/extract_functions.R` and add entries to the `EXTRACTIONS` list:

```r
, list(
    target = "R/func_TARGET.R"
    , source = "R/SOURCE_FILE.R"
    , functions = c(
      "functionName1"
      , "functionName2"
      , "functionName3"
    )
  )
```

### Step 3: Test with Dry Run

```r
source("dev/extract_functions.R")
extract_all(dry_run = TRUE)
```

### Step 4: Execute and Verify

```r
extract_all(dry_run = FALSE)
devtools::document()
```

---

## Remaining Work (Future Sessions)

### S4 Class Extraction (Optional)
1. Extract S4 class definitions to `func_prot_s4_objects.R` and `func_metab_s4_objects.R`
2. Note: S4 methods are tightly coupled to classes - extraction may not be beneficial
3. S4 classes in `functional_enrichment.R` (`de_results_for_enrichment`, `EnrichmentResults`) remain there

### Cleanup (Optional)
1. Remove extracted functions from source files (functions are now duplicated)
2. Consider removing `allGenerics.R` since content is now in `func_general_s4_generics.R`
3. Update any internal references
4. Final `devtools::document()` and testing

---

## Reference Documents

- `FUNC_CLASSIFICATION.md` - Original function mappings (many corrections applied)
- `FUNC_FILES_INDEX.md` - Overview of func_*.R file structure
- `REORGANIZATION_SUMMARY.md` - Reorganization rationale
- `dev/extraction_log.md` - Full extraction log
- `dev/extraction_log_errors.md` - Functions that need investigation

---

## Notes for Handover

1. **The extraction tool works well** - all errors resolved, 312+ functions extracted
2. **312+ functions successfully extracted** across 22 populated target files
3. **63 S4 generics extracted** to func_general_s4_generics.R
4. **7 redundant placeholder files deleted** - cleaner R/ directory
5. **Package builds correctly** - devtools::document() passes
6. **Source file mappings corrected** - many functions were in different files than documented
7. **Always run dry_run first** - prevents duplicate extractions
8. **Phase 2 complete** - All 4 overlooked files now have functions extracted

---

*Last Updated: 2025-11-27*
*Session 3 (Phase 2): 44 functions extracted from 4 overlooked files*
