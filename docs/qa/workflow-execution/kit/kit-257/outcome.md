# KIT-257 outcome: each normalisation first removes its own earlier files

Ticket `tickets/kit/KIT-257.md` r1.2 (step 0 applied Re-check 1's R1-01 to R1-05). Implements VAL-250. Branch `kit-257`
from `v0.9-staging` 63143aa3. KIT-243 had not merged at the branch point, so the removal core is this ticket's
`removeWorkflowOwnFiles` (`unlink(files, expand = FALSE)`). No R was run on the laptop; every run below is a GitHub
runner (`run-tests`).

## Commits

| Commit | What |
|---|---|
| 1e1b0cee | tests only: `tests/testthat/test-kit257-norm-own-files.R` (nothing under `R/`) |
| 413d739a | the fix (`R/`), NEWS, ledger row `VAL-<merge>` |
| 5ad937e4 | records: R1, wording drafts |
| (next) | records: R2, outcome, runs |

## What changed (at 413d739a)

- `removeWorkflowOwnFiles(files, message)` (`R/utils_workflow_saved_exports.R:674`): the present, non-directory files
  removed with `unlink(expand = FALSE)`; any still present raises `workflow_figure_unwritten`; one logger line names
  the removed files.
- Transcriptomics, lipidomics, metabolomics: `omicNormApplyRunStart` gains `remove` (`R/mod_omic_norm.R:999`), called
  after `omicNormRunStart`'s refusals and before the previous start is selected, inside the same `tryCatch` (a failure
  is returned as a refusal); one tab-log line when anything was removed. The lists: `omicNormOwnFiles`
  (`R/mod_omic_norm.R:1025`): trx CCA panels (png, pdf) and grid workbook in `publication_graphs_dir`; lipid/metab
  `<assay>_ruv_corrected_{pca,density,rle,correlation}.png` in the omic's QC folder, assay names from the tab, else
  the head's record (R1-03). Callers: `runOmicNormPipeline` (`R/mod_omic_norm.R:456`), `runTrxNormalisation`
  (`R/utils_workflow_trx_norm.R:142`).
- Proteomics: `removeProtNormOwnFiles(paths, stage)` (`R/utils_workflow_dia_norm.R:670`): at a normalisation the five
  RUV figures and both results-table names (TSV and RDS); at a correlation the four table files. Positions:
  in process `runDIAReferenceNormalization` directly after the request (`R/utils_workflow_dia_products.R:748`) and
  `runDIAReferenceCorrelation` after `diaCorrelationSource` (`:1812`); managed normalisation in the actions function
  resolved when the actions start (`R/mod_prot_norm_workflow_helpers.R:638`, R1-01); managed correlation in its first
  action function (`:732`); DIA with limpa in process in the tab's entry before the pass-through (`:510`), managed in
  `dispatchProtNormPassThrough`'s action function (`:270`). The pass-through itself is not changed: its body is pinned
  by digest (`tests/testthat/test-prot-norm-core.R:45`).
- Texts (`omicNormGuardTexts`): `norm_unremoved`, `correlation_unremoved` (drafts, `wording.md`).

## Runs

See `runs.tsv`.

### R1 (45618988736, tests-only 1e1b0cee, 1 job, 5.1 test minutes)

Observed (`r1-tests.tsv`, `r1-failures.txt`, `r1-evidence.txt`): 51 passed, 36 failed, 0 errors, 0 skipped. Every
failure is an assertion. Failing tests: E1 to E9, G3, G5 (proteomics, DIA with limpa, lipid, metab), each at its base
value:

| Row | Base value observed |
|---|---|
| E1 | after the declared run the three CCA files present |
| E2 | at the first dispatch the three planted files present |
| E3 | RUV table present beside the current one; Copy: "Using RUV-normalized file names", `RUV_normalised_results.tsv` and `ruv_normalised_results.RDS` delivered |
| E4 | the five planted RUV figures present |
| E5 | lipid and metab: the four figures present; Copy put all four in `QC_figures/` |
| E6 | lipid and metab: present at the first dispatch |
| E7 | the DA ran at `normalized` (reachable; not refused); the first run's RUV table present; Copy delivered it as `RUV_normalised_results.tsv` |
| E8 | after the second pass-through all four table files present |
| E9 | the correlation's table present when the queued normalisation started |
| G3 | the run proceeded, committed and failed at the publication ("the publication's files could not be copied into the project"); `removeWorkflowOwnFiles` absent |
| G5 | every route proceeded: no refusal; the head moved where the route commits |

G1 (trx, proteomics, lipid) and G2 (trx and proteomics user files; lipid users' files, another omic's folder, the
`Study [pilot]` / `Study p` sibling) passed.

Inferred: the reading of every omic was right (stop rule 1 not met); E7's sequence is reachable in the application
(the proteomics DA accepts `normalized`), so the proteomics defect is not hypothetical.

### R2

(to be filled from the run)

## Acceptance criteria

(to be filled from R2)

## Not checked

- The managed worker paths on a real worker (E2, E6, E9 and the managed G5 cases test the dispatch in the session with
  a recording runner).
- Windows file locks (the case the refusal exists for); the directory stands in for a file that cannot be removed.
- Metabolomics with RUV applied (G1): F172 runs RUV for lipidomics only (KIT-194, KIT-207 precedent); metabolomics
  shares the code path and its removal is shown by E5, E6 and G5.
- G4: KIT-243's file is not on staging at the branch point (not applicable).
- Projects made by earlier builds on real disks.

## For the change report

### Main part (legacy; confirmed with `git show main:<path>`, main 436cdbe8)

**Proteomics: Copy delivered an earlier analysis's normalised table: wrong data in the results package.** The
correlation step writes its results table under a name chosen by the RUV setting
(`main:R/mod_prot_norm.R:2078-2079`, `:2273-2274`: `ruv_normalised_results_cln_with_replicates.tsv` with RUV,
`normalised_results_cln_with_replicates.tsv` without), and nothing removed the other name. Copy prefers the RUV name
whenever that file exists (`main:R/func_general_filemgmt.R:2241-2243`; at the base commit
`R/func_general_filemgmt_results_summary_helpers.R:381-387`). After a normalisation with RUV and the correlation, then
a normalisation without RUV (or with another k, when the differential analysis is then run without a new correlation
step), Copy put the earlier analysis's RUV-normalised abundance table into the results package as
`Publication_tables/RUV_normalised_results.tsv` (and, outside native sessions, its `.RDS`): the package's only
normalised table, while the differential-analysis tables, the report and the study parameters describe the later
analysis. Nothing marked it as earlier. Who was affected: every proteomics route (DIA, LFQ, TMT; DIA with limpa when
the correlation is applied again after the RUV covariates), in every earlier release. The report's text was not
affected (no proteomics report uses the table's values). How to recognise an affected package: a
`RUV_normalised_results.tsv` beside a report or study parameters saying RUV-III was not applied. Remedy: run the
normalisation (and the correlation step) again, then Copy. The earlier RUV-corrected figures (`ruv_corrected_*.png`,
`ruv_corrected_cancor.png`) also stayed and the native session archive packed them.

**Metabolomics: Copy collected earlier RUV-corrected QC figures.** Copy copies every
`_(pre_norm|post_norm|ruv_corrected)_(pca|density|rle|correlation).png` of the QC folder
(`main:R/func_general_filemgmt.R:2285`; base `R/func_general_filemgmt_results_summary_helpers.R:429`). The RUV
figures are drawn only when RUV runs, so after a run with RUV skipped `QC_figures/` carried the earlier run's
RUV-corrected figures. The composite and the report did not use them.

### Appendix (new to 0.9)

- **Transcriptomics** (no transcriptomics on `main`): the legacy RUV grid writes the CCA panels and the grid workbook;
  a declared design does not (base `R/utils_workflow_trx_declarations.R:829-847`). Copy collected them by name
  (`R/func_general_filemgmt_results_summary_helpers.R:935`, `:939`) and the report showed the panel whenever the file
  existed (`inst/reports/transcriptomics/transcriptomics_report.rmd:131-133`), beside Methods describing the declared
  design.
- **Lipidomics**: the same Copy pattern as metabolomics (base `:527-531`); no such Copy on `main` (git grep).

### Both parts

- What changed in results: no number changes. What Copy and the report deliver now matches the analysis run: after
  the fix the earlier files are absent (Copy names a missing file where nothing current replaces it).
- Seams: base 63143aa3 had none; fix 413d739a at the lines listed under "What changed".
- Projects processed before the fix keep their earlier files until the stage is run again. Revert after a re-run
  leaves the removed files absent until the step is run again; the same after a proteomics normalisation or correlation
  that fails after it has started.
- A file that cannot be removed now refuses the run with a message naming it.
- Evidence: R1 45618988736 (base values above); R2 (below).
- Not in this ticket: F1, the transcriptomics RUVIII export is never cleared, so after a later normalisation without
  RUV the report can still state that RUVIII was applied with an earlier k and show its figure. Routed to its own
  ticket (parked-20261125 row 38 (2)); until it has an owner this report states it.
