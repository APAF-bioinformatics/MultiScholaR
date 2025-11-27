# Functions Not Found - Investigation Required

**Generated:** 2025-11-27 16:10:21

**Total functions not found:** 13

These functions were listed in FUNC_CLASSIFICATION.md but could not be
located in their expected source files. Possible reasons:

1. Function name differs from classification (typo, case sensitivity)
2. Function is in a different file than documented
3. Function was renamed or removed
4. Function is defined inline within another function
5. Function uses a different definition pattern (e.g., = instead of <-)

---

## Functions by Target File

### func_prot_qc.R

| Function | Expected Source | Action Needed |
|----------|-----------------|---------------|
| `proteinIntensityFiltering` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |
| `removeProteinsWithOnlyOneReplicate` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |
| `removeRowsWithMissingValuesPercent` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |

### func_prot_norm.R

| Function | Expected Source | Action Needed |
|----------|-----------------|---------------|
| `normaliseBetweenSamples` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |
| `normaliseUntransformedData` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |
| `ruvIII_C_Varying` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |
| `ruvCancor` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |
| `ruvCancorFast` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |
| `getNegCtrlProtAnova` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |

### func_metab_qc.R

| Function | Expected Source | Action Needed |
|----------|-----------------|---------------|
| `resolveDuplicateFeatures` | metabolite_qc.R | [ ] Locate/Reclassify |

### func_metab_de.R

| Function | Expected Source | Action Needed |
|----------|-----------------|---------------|
| `differentialAbundanceAnalysis` | metaboliteVsSamplesS4Objects.R | [ ] Locate/Reclassify |
| `differentialAbundanceAnalysisHelper` | metaboliteVsSamplesS4Objects.R | [ ] Locate/Reclassify |

### func_general_design.R

| Function | Expected Source | Action Needed |
|----------|-----------------|---------------|
| `cleanDesignMatrix` | proteinVsSamplesS4Objects.R | [ ] Locate/Reclassify |

---

## Functions by Source File (for grep investigation)

Use these commands to search for functions in the codebase:

```bash
# Search for a function definition
grep -n -E '(proteinIntensityFiltering|removeProteinsWithOnlyOneReplicate|removeRowsWithMissingValuesPercent|normaliseBetweenSamples|normaliseUntransformedData|ruvIII_C_Varying|ruvCancor|ruvCancorFast|getNegCtrlProtAnova|cleanDesignMatrix)' R/*.R
grep -n -E '(resolveDuplicateFeatures)' R/*.R
grep -n -E '(differentialAbundanceAnalysis|differentialAbundanceAnalysisHelper)' R/*.R
```

## Quick Reference: All Missing Functions

```
proteinIntensityFiltering -> func_prot_qc.R (expected in proteinVsSamplesS4Objects.R)
removeProteinsWithOnlyOneReplicate -> func_prot_qc.R (expected in proteinVsSamplesS4Objects.R)
removeRowsWithMissingValuesPercent -> func_prot_qc.R (expected in proteinVsSamplesS4Objects.R)
normaliseBetweenSamples -> func_prot_norm.R (expected in proteinVsSamplesS4Objects.R)
normaliseUntransformedData -> func_prot_norm.R (expected in proteinVsSamplesS4Objects.R)
ruvIII_C_Varying -> func_prot_norm.R (expected in proteinVsSamplesS4Objects.R)
ruvCancor -> func_prot_norm.R (expected in proteinVsSamplesS4Objects.R)
ruvCancorFast -> func_prot_norm.R (expected in proteinVsSamplesS4Objects.R)
getNegCtrlProtAnova -> func_prot_norm.R (expected in proteinVsSamplesS4Objects.R)
resolveDuplicateFeatures -> func_metab_qc.R (expected in metabolite_qc.R)
differentialAbundanceAnalysis -> func_metab_de.R (expected in metaboliteVsSamplesS4Objects.R)
differentialAbundanceAnalysisHelper -> func_metab_de.R (expected in metaboliteVsSamplesS4Objects.R)
cleanDesignMatrix -> func_general_design.R (expected in proteinVsSamplesS4Objects.R)
```

