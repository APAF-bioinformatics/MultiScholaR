# KIT-112 W5, chain A: Simulation Setup, run 1, run 3 (2026-11-10)

**Setup:** the merged worktree at `a53f9060`, `load_all`, R 4.6.1, strict skips, `NOT_CRAN=true`, `R_LIBS` = `kit_111`.

- **Simulation Setup:** `prot_dia` (241 s) and `prot_dia_limpa` (162 s) are both identical to their oracles. In each, `runDIAReferenceDaAnalysis`, `diaDaRequest`, `computeDIAReferenceDaContrast` and `acceptDIADaManifest` are called once (one contrast). `computeWorkflowDiaContrastAction` is called 0 times, as expected for an in-process session. Traces: `sim-prot_dia*.txt`.
- **Run 1:** 54 files, 16,994 passed, 0 failed, 0 errors ([../targeted.txt](../targeted.txt)).
- **Run 4, inside run 1:** the DIA session and memory set gives 733/2, 191/1, 38/1, 174/4 and 80/11, each equal ([../dia-memory-set.txt](../dia-memory-set.txt)).
- **Run 3:** trx 29, lipid 35, metab 105, proteomics 193 passed, with 0 failed and 0 skipped. All ten lanes show 0 step differences against their oracles. No re-freeze.
