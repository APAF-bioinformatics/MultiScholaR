# KIT-212 outcome (ticket r1.4)

Branch `kit-212` from `v0.9-staging` 32766889 (no KIT-216 code in the base: E6's and G4's "base without KIT-216" columns
apply). No R was run on the laptop; every run is listed in `runs.tsv`.

## Step 0 (before R1)

### The survey (KIT-207's form, R2's scope)

Observed by reading at 32766889 (no R run):

- **W1 Summary captures.** The nine cases of section 6 (r1.4) move: `lipid_main`, `lipid_variants`, `lipid_github_fallback`,
  `lipid_global_another`, `metab_main`, `metab_variants`, `metab_github_fallback` (module-CI owners, DA tables on disk, no DA),
  `prot_owner`, `prot_owner_variants` (a session owner, no DA, DA tables on disk). The six proteomics cases on the in-memory
  manager are not checked (S1). `test-omics-summary-core.R` runs lipid, metab and trx cases only.
- **Save, Copy or Report after a table change after a DA, in merged tests.** One move, declared as E14: KIT-207 E11
  (`test-kit207-summary-head-guard.R:604-633`, Copy after a Revert to `lipid_qc_complete`; the DA's source left the history,
  so Copy refuses with `summary_da_stale` before it checks the publication artifacts). Every other KIT-207 row that saves
  after a table change is refused first by KIT-207's own guard (unchanged); its G9 (Reset, then one Revert) keeps the DA's
  source in the history with the same reference, so it stays "current".
- **E15 is not made: KIT-206 E11b does not move.** Its DA is `k206Da` (`test-kit206-prot-norm-export-guard.R:83-86`), which
  only re-commits the DA arguments (`commitDIADaArguments`) and accepts no manifest, so `snapshot$da_results` is NULL, the
  status is "none", and KIT-212's Save guard passes; KIT-206's refusal is reached as before. R2 checks this (the file must be
  `same`). Re-check 3's (m2) assumed a DA manifest.
- **Other files with a DA and a Summary call** (`test-kit172-downstream.R`, `test-workflow-analysis-record.R`,
  `test-workflow-omics-application-reads.R`, `test-workflow-figure-paths.R`, `test-workflow-memory-mechanism-member.R`,
  `test-workflow-checkpoint-harness.R`, `test-workflow-kit-harness.R`, `test-module-ui-wrapper-direct-shared.R`, the GUI-entry
  lanes): the Summary runs at the DA's head, or on a manager that is not a session owner, or without DA files. Not a move by
  reading; R2 judges them.
- **Replay** (the auditor's re-check 3, confirmed): every replaying test record has its DA source among its steps.
- **DA tab observer:** the DA module captures' buses are plain environments, so the observer cannot fire there.
- **F172 threshold:** not needed. E4, E9 to E9c, E10, G11 and G12 do not depend on a sample being removed: a Revert or a
  second run takes the DA's source out of the active history ("stale" by history), and E9c uses the DA on `lipid_log2`
  followed by the normalisation (another table at the head). G11 and G12 print whether the table is the same.

### G0: which DA files the files check judges, per omic

The candidates are the files Copy collects whose names follow the omic's DA publication naming; the made list is computed
from the manifest by the writers' own (now callable) naming. Anything else in the folders (e.g. a fixture's `artifact.txt`)
is "other" and never refused.

| Omic | Folder (Copy) | Candidate naming | Made by (callable) | Kind |
|---|---|---|---|---|
| lipid | `da_output_dir` | `(da\|de)_.+_long_annot.xlsx` | `omicsDaFileNames("lipids", assay, contrast)$table` | per contrast, assay |
| lipid | `Heatmaps` | `<mode>_<contrast>_heatmap.(png\|pdf)`; `all_heatmaps_combined.pdf` | `omicsDaFileNames(...)$plot`; per run | per contrast, assay; per run |
| lipid | `NumSigDeMolecules` | `lipids_num_sig_de_molecules.(tab\|xlsx)`, `lipids_num_sig_barplot.png` | per run | per run |
| metab | `da_output_dir`, `Heatmaps` | as lipid (`metabolites`) | `omicsDaFileNames("metabolites", ...)` | as lipid |
| trx | `da_output_dir` | `de_results_<comparison>.csv` | `trxDaPublicationNames(comparison)$table` | per contrast |
| trx | `Volcano_Plots`, `NumSigDeMolecules` | `volcano_plot_<c>.(png\|pdf)`, `summary_barplot_<c>.(png\|pdf)` | `trxDaPublicationNames` | per contrast |
| prot | `da_output_dir`, `Volcano_Plots`, `NumSigDaMolecules` | tables, `<safe>.(png\|pdf)`, `all_volcano_plots_combined.pdf`, `da_proteins_num_sig_da_molecules.*` | `diaDaPublicationFiles` (the publication's own list, `assertDIADaPublicationFiles`), plus the S4 kernel's `.png`/`.pdf` count plots (`R/func_prot_da_results_methods.R:377-383`) | per contrast; per run |
| prot | `Interactive_Volcano_Plots` | `<friendly>.html`, `<friendly>.unavailable.txt` | `diaInteractiveVolcanoFileName`, `diaInteractiveVolcanoUnavailablePath` | per contrast |
| prot | `Heatmap` (Save Heatmap) | `<prefix>_(heatmap.(png\|pdf)\|heatmap_methods.md\|clusters.csv)` | `omicDaHeatmapSavePrefix("prot_", contrast)` for the manifest's contrasts (any of its three names) | user-saved, by contrast |

Lipid and metab Copy do not collect `Volcano_Plots`; their user-saved heatmaps go to `Heatmap`, which only proteomics' Copy
collects. The lipid and metab names are pinned by G1 and E12 (the real publication's files, none unmade), trx's and
proteomics' by G1.

## Runs

See `runs.tsv`.

## Not checked

(Filled in at the audit point.)
