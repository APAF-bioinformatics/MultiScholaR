# KIT-097: stopped and superseded (2026-11-06)

**Stopped by the user's decision before closure.** The normalisation, DA and summary tabs were built by adding transcriptomics `profile` hooks to the lipid modules: about 67 branches across the lipid DA, normalisation and summary files, a 455-line change to `mod_lipid_norm_ui_helpers.R`, and about 1,070 lines of trx module code. The user judged that "lazy and incomplete". The lipid modules become an accidental shared base, metabolomics keeps its own twin copy, and lipid assumptions leak into transcriptomics. The e2e defect below is one such leak.

## State at the stop

**Committed:**
- `efb8b605`: the plan.
- `55a53d44`: the tabs.
- `8e97b55e`: the coordinator-review fix for a function-name clash. `trxDaContrasts` in `mod_trx_da.R` masked KIT-093's kernel and broke DA, and `R CMD check` showed a fourth WARNING.

**Install:** `kit_097` = QqMykR, from 8e97b55e.

**Evidence in this folder** (checks run on the stopped code):
- `R CMD check` at the 3-WARNING baseline after the rename;
- the trx kernel, design, binding, module-CI and application-reads tests pass;
- the lipid, metab and proteomics module CI pass with their test files unchanged.

**The e2e failed** on its counted attempt (`e2e-failed-attempt1.txt`): `Timed out waiting for transcriptomics/differential_analysis to become complete; last status was <NULL>`, after import, design, QC and normalisation completed, in 1,264 s. The likely cause is that the lipid DA helpers' tab-status writes don't respect the transcriptomics key. The agent had started a reproduction test when it was stopped; that unfinished change is kept as `stopped-wip-test-module-ci-trx-app.patch`, and the test file is restored.

## What supersedes it

A ticket that does for the GUI what KIT-073 did for the declarations: shared omics GUI modules driven by each omic's binding, with lipid, metab and transcriptomics as configurations. KIT-096's and KIT-097's transcriptomics tabs get reworked onto it. The committed hook code stays until that ticket replaces it. It serves as a working behavioural reference only, not the final structure.
