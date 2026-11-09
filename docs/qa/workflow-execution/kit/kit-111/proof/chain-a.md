# KIT-111 W6, chain A: Simulation Setup, run 1, run 3 (2026-11-10)

Tree: the merged worktree at `bcd3396b` (R code as `341dab60`), `load_all`, R 4.6.1, strict skips, `R_LIBS` = `kit_108`.

## Simulation Setup (`kit-105/lifecycle/lifecycle.R`, traced)

| Lane | Identical to oracle | Calls traced |
|---|---|---|
| `prot_dia` (370 s) | TRUE | runner, compute, commit, publish: 1 each for normalisation and correlation; `computeWorkflowDiaNormalizationAction` 0 (in-process session: the runner runs in place); `installWorkflowOmicsActions` 1 |
| `prot_dia_limpa` (260 s) | TRUE | pass-through 1; correlation runner/compute 1/1; RUV covariates runner, compute, commit, publish 1 each; `computeWorkflowDiaRuvCovariatesAction` 0 |

Traces: `sim-prot_dia.txt`, `sim-prot_dia_limpa.txt`.

## Run 1: the 49 targeted files

49 files, 15,753 passed, 0 failed, 0 errors (`run1.tsv`). The DIA memory set within it: `test-workflow-native-session.R` 733 passed / 2 skipped; `test-workflow-dia-products.R` 189/0; `test-module-ci-prot-dia-limpa.R` 33/2; `test-module-ci-prot-norm.R` 60/4; `test-module-ci-prot-da.R` 218/1.

## Run 3: the entry-point suite with records

trx 29, lipid 35, metab 105, proteomics 193 passed, 0 failed, 0 skipped. Every lane's record against its oracle (`jsondiff.py`): 0 step differences in all ten lanes.
