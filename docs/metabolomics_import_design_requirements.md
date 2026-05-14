# Metabolomics Import Design Requirements

For the MultiScholaR Metabolomics module's "Import Existing Design" functionality, the app expects you to select a **folder** containing specific files. 

Here is what is required and how to construct them:

## 1. Required Files
The import script explicitly checks for the following mandatory files in your selected folder:
*   **`design_matrix.tab`**: Must be present.
*   **`assay_manifest.txt`**: Must be present. 
*   **Assay Data Files (`data_cln_*.tab`)**: The script reads the names listed in `assay_manifest.txt` and explicitly checks that a corresponding `data_cln_[name].tab` file exists for each entry.

## 2. Optional Files
*   **`contrast_strings.tab`**: If present, it loads the predefined differential comparisons.
*   **`column_mapping.json`**: Explicitly maps columns for features and sample names. If missing, the app tries to auto-detect them.
*   **`config.ini`**: App configuration. Falls back to a default if missing.

## 3. The Role of `assay_manifest.txt` vs `manifest.json`
The Metabolomics module allows you to analyze multiple separate assay tables simultaneously (e.g., Positive and Negative modes). 
Instead of hardcoding them into the general `manifest.json`, the script reads `assay_manifest.txt` to see exactly which assays exist. 
For example, if `assay_manifest.txt` contains the lines `LCMS_Pos` and `LCMS_Neg`, the script dynamically constructs the filenames `data_cln_LCMS_Pos.tab` and `data_cln_LCMS_Neg.tab` and strictly requires them to be in the folder.

## 4. How to Construct These Files Manually

### `assay_manifest.txt`
A simple plain text file with one assay identifier per line (NO file extensions or prefixes).
```text
LCMS_Neg
LCMS_Pos
```

### `design_matrix.tab`
A tab-separated file that defines your experiment. It *must* include matching columns used internally: `Run`, `group`, and `replicates`. 
```text
Run	group	replicates
Sample1_Neg	Control	1
Sample2_Neg	Control	2
Sample3_Neg	Treated	1
```

### `data_cln_LCMS_Neg.tab` (and other assays)
A tab-separated data matrix. The sample column headers must **exactly match** the values in the `Run` column of your `design_matrix.tab`. It also needs a unique feature ID column (e.g., "Peak ID") and an annotation column (e.g., "Name").
```text
Peak ID	Name	Sample1_Neg	Sample2_Neg	Sample3_Neg
101	Glucose	14500	15200	9800
102	Fructose	8000	8100	4000
```

### `contrast_strings.tab`
A plain text file containing the exact statistical comparisons you want Limma to run, formatted as `group[Level1]-group[Level2]`.
```text
groupTreated-groupControl
```

If you format a folder with these generated text files, you can point the app to it via the "Import Existing Design" button to bypass the UI setup entirely!

---

<!-- APAF Bioinformatics | metabolomics_import_design_requirements.md | Approved | 2026-03-19 -->
