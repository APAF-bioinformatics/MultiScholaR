# KIT-106 W3: the five W1 values that D3 and D4 change, re-frozen once

Lead decision, 2026-11-09: re-freeze these five W1 expectations once, with attribution. The stop rule is not applied.
All five are in `tests/testthat/test-omic-da-module.R` (`omicDaW1Frozen`), and each carries a comment that cites this
decision. No other value changed.

## The values

| Case, step, field | Before (W1, c3a35150) | After (W3 tree) | Cause |
|---|---|---|---|
| `lipid_managed` select `notes` | `[default NULL lipidomics_action_working] Applying lipidomics step ...`, `remove lipidomics_action_working` | `character(0)` | D3 |
| `lipid_managed` select `files` | the run step's files plus the page's store entries (`omics-volcano-html` 1; `dia-state` 14 -> 15) | `"<unchanged>"` | D3 |
| `lipid_managed` select `operations` | `"lipid_volcano_page"` | `character(0)` | D3 |
| `lipid_no_friendly` run `calls` | the table widget only | the export's 4 `generateLipidDAVolcanoStatic` and 4 `generateLipidDAHeatmap` calls (KO and RES, LCMS_Pos and GCMS), then the table widget | D4 |
| `lipid_no_friendly` run `files` | 4 TSVs and their XLSX (tables without `friendly_name`), NumSigDE; store: `dia-state` 13, no figure entities | the 4 TSVs with `friendly_name` (new digests), their XLSX, NumSigDE, every per-assay volcano and heatmap PNG/PDF and both combined PDFs; store: `omics-figure` 11, `omics-figure-asset` 19, `dia-state` 43 | D4 |

- **D3, the volcano on demand.** The box is unticked by default (user unification, 2026-11-06), so the managed page
  is not requested at the select step. Under `testServer`, outputs render when the session flushes, even if the test
  never reads them. At c3a35150 the box was ticked by default, so selecting an assay (GCMS) rendered
  `volcano_glimma`, and the managed factory asked the worker for the page. Now the render stops at
  `shiny::req(isTRUE(input$volcano_interactive))`.
- **D4, the in-process export through the publication runner.** The in-process export now goes through
  `runWorkflowOmicsPublicationInProcess(..., "lipid_da_publication")`. That runner reads the head's accepted DA manifest,
  not the tab's list. The case mocks `friendly_name` out of the reference DA's returned list only, so the tab's list
  lacks the column while the manifest keeps it. At c3a35150 the exporter ran on the tab's list, and its volcano and
  heatmap kernels failed on the missing column, so it wrote the tables and NumSigDE only. Now the export reads the
  manifest's rows, with `friendly_name`, and draws every figure.

## Proof that each cause is the only one

Each proof is a scratch run, not committed, through `run1.R` with `load_all`, `MULTISCHOLAR_STRICT_SKIPS=1` and
`TMPDIR` on disk.

- **D3:** `test-attr-d3.R` runs the W3 tree with one change, the box's default set back to ticked in
  `omicDaFeatureLayout`. Result: **65/0**. With the old default, the frame reproduces W1's frozen `lipid_managed`
  capture exactly, at every step and in every field.
- **D4:** `test-attr-d4.R` runs the W3 tree unchanged, with `readWorkflowOmicsDaResults` mocked to drop
  `friendly_name`, so the publication sees the same rows as the tab. In `lipid_no_friendly`'s run step, every kernel
  call and every project file (names and table digests) then equals W1's frozen value. Two store-entity counts still
  differ:
  - `omics-figure` 1 and `omics-figure-asset` 1 are added;
  - `dia-state` goes from 13 to 15.

### The store-count difference that remains

- **The entity:** one `omics-figure` scientific product and its `omics-figure-asset` rendered PNG. Both are
  `publication_graphs/NumSigDeMolecules/lipids_num_sig_barplot.png`, written through `writeWorkflowOmicsFigure`
  (`func_lipid_da_export.R:412`). It is the one figure the exporter can still draw on rows without `friendly_name`.
- **Why the new path creates it:** `writeWorkflowOmicsFigure` keeps a figure in the store only when it finds an owner
  (`workflowOmicsFigureOwner`). It looks first on the active state scope, then on the session's registered owners.
  `runWorkflowOmicsPublicationInProcess` runs the export inside `withWorkflowStateScope(manager, ...)`, so the owner is
  always found. In this case the owner was made by `omicsTestOwner` in a separate owner session and is not registered
  on the module's session. The start commit's exporter ran outside any scope, found no owner, wrote the PNG only, and
  kept nothing in the store.
- **The `dia-state` +2:** the two `rootWorkflowProduct` records that root the product and its asset under the head.
- **In the re-frozen value:** the same mechanism, applied to every figure the manifest's rows allow, gives
  `omics-figure` 11, `omics-figure-asset` 19 and `dia-state` 13 + 30 = 43. These are the same figure counts that
  `lipid_main` and `lipid_managed` record.
- **`lipid_main`, the realistic case, shows no store-count difference.** Its run step's frozen capture already holds
  `omics-figure` 11 and `omics-figure-asset` 19: its owner is adopted by the tab's session load on the module's session,
  so the start commit's exporter found it there too. `lipid_main` passes identically on the W3 tree (every step and
  every field, the store counts included).

## Freezing

The new values come from a run of the W3 tree (43517efa's R/, which `scratch attr2/` copies). They were dumped as R
code by a scratch test that writes the five fields with `omicDaW1Code`. Two dumps in separate processes were
byte-identical (`cmp`; 6,053 bytes), and the values were pasted from the first.

## Result

`test-omic-da-module.R` (`run1.R`, `load_all`, strict skips, `TMPDIR` on disk): `COUNTS test-omic-da-module.R pass 1218
fail 0 error 0 skip 0`.

## Draft quarantine rows (W3; for `quarantined-tests.md`)

| Id | Test | What it covered | Owner | Replacement coverage |
|---|---|---|---|---|
| Q-130 | `test-omics-da-core.R`: "(W2, D4; comparison side) lipid: the in-process publication equals the start commit's write" | It compared the in-process publication with `outputLipidDaResultsAllContrasts` on `lipidDaResultsWithSource` (mod_lipid_da_session_helpers.R, deleted by W3) | KIT-106 W3 | None expected. The core-side block stays live, as does `test-omic-da-module.R` (b) lipid_main's run step |
| Q-133 | `test-lipid-02b-da-module-contracts.R`: the whole file (line-1 skip; its `source()` lines name deleted files), 77 blocks | It sourced the five lipid DA helper files and mod_lipid_da.R's registrations, and tested their seams | KIT-106 W3 | None expected. The shared frame mod_omic_da.R and panel library mod_omic_da_panels.R over the W2 core, covered by `test-omic-da-module.R` (a) to (c) and its W3 blocks, MCI-022, `test-omics-da-core.R` and the lipid_canonical DA steps |
| Q-134 | `test-lipid-02ab-da-module-characterization.R`: "mod_lipid_da_server preserves run-analysis public success behavior" | It counted calls to a mocked `outputLipidDaResultsAllContrasts` from the in-process run, on a mocked `runLipidReferenceDaAnalysis` that accepts no manifest | KIT-106 W3 | None expected. The export is the binding's `lipid_da_publication` through `runWorkflowOmicsPublicationInProcess` (D4), covered by `test-omic-da-module.R` (b) lipid_main's run step and `test-omics-da-core.R` |
| Q-135 | `test-lipid-02ab-da-module-characterization.R`: "mod_lipid_da_server preserves render and heatmap-save public behavior" | It read `volcano_glimma` without ticking the box and expected a mocked `generateLipidDAVolcanoPlotGlimma` on the tab's results | KIT-106 W3 | None expected. The volcano is on demand, as the binding's `lipid_volcano_page` (D3), covered by the W3 D3 blocks and (b) lipid_main (heatmap, cluster summary, Save Heatmap's params and prefix, the table) |
