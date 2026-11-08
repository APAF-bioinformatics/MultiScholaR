# KIT-106 W4: the five metab W1 values that D3 and D4 change, re-frozen once

Lead decision, 2026-11-09 (pre-authorised for W4, as for W3): `metab_managed`'s select step and `metab_no_friendly`'s run
step may be re-frozen once, only for D3 or D4, with attribution. The stop rule is not applied. All five values are in
`tests/testthat/test-omic-da-module.R` (`omicDaW1Frozen`), and each carries a comment that cites this decision. No other
value changed. They are the metab twins of W3's five lipid values ([../w3/w1-refreeze.md](../w3/w1-refreeze.md)).

## The values

| Case, step, field | Before (W1, frozen) | After (W4 tree) | Cause |
|---|---|---|---|
| `metab_managed` select `notes` | `[default NULL metabolomics_action_working] Applying metabolomics step ...`, `remove metabolomics_action_working` | `character(0)` | D3 |
| `metab_managed` select `files` | the run step's files plus the page's store entries (`omics-volcano-html` 1; `dia-state` 14 -> 15) | `"<unchanged>"` | D3 |
| `metab_managed` select `operations` | `"metab_volcano_page"` | `character(0)` | D3 |
| `metab_no_friendly` run `calls` | the two table widgets only | the export's 4 `generateMetabDAVolcanoStatic` and 4 `generateMetabDAHeatmap` calls (KO and RES, LCMS_Pos and GCMS), then the two table widgets | D4 |
| `metab_no_friendly` run `files` | 4 TSVs and their XLSX (tables without `friendly_name`), NumSigDE; store: `dia-state` 13, no figure entities | the 4 TSVs with `friendly_name` (new digests, equal to `metab_managed`'s), their XLSX, NumSigDE, every per-assay volcano and heatmap PNG/PDF and both combined PDFs; store: `omics-figure` 11, `omics-figure-asset` 19, `dia-state` 43 | D4 |

Before the re-freeze, `test-omic-da-module.R` on the W4 tree: `pass 1229 fail 5`; the five failures were exactly these
five fields (`kit-106-w4/after1/test-omic-da-module.R.txt` in the agent's scratch).

- **D3, the volcano on demand.** The box is unticked by default (user unification, 2026-11-06), so the managed page is
  not requested at the select step. Under `testServer`, outputs render when the session flushes. With the old module the
  box was ticked by default, so selecting GCMS rendered `volcano_glimma` and the managed factory asked the worker for
  the page. On the frame the render stops at `shiny::req(isTRUE(input$volcano_interactive))`.
- **D4, the in-process export through the publication runner.** The frame's in-process export is
  `runWorkflowOmicsPublicationInProcess(..., "metab_da_publication")`, which reads the head's accepted DA manifest, not
  the tab's list. The case mocks `friendly_name` out of the reference DA's returned list only, so the tab's list lacks the
  column while the manifest keeps it. The old module's `writeMetabDaResultsArtifacts` ran `outputMetabDaResultsAllContrasts`
  on the tab's list, whose volcano and heatmap kernels failed on the missing column, so it wrote the tables and NumSigDE
  only. Now the export reads the manifest's rows, with `friendly_name`, and draws every figure.

## Proof that each cause is the only one

Scratch runs, not committed, `load_all`, `MULTISCHOLAR_STRICT_SKIPS=1`, `TMPDIR` on disk, in testthat's reproducible
context (collation C).

- **D3:** a copy of the W4 tree with one change, the box's default set back to ticked in `omicDaFeatureLayout`
  (`R/mod_omic_da_panels.R:146`, `value = TRUE`). The per-field comparison of `metab_managed` against W1's frozen capture
  (`case_diff.R`): **0 differing fields**, at every step. With the old default the frame reproduces `metab_managed`
  exactly, store counts included.
- **D4:** the W4 tree unchanged, with `readWorkflowOmicsDaResults` mocked to drop `friendly_name` from
  `da_metabolites_long`, so the publication sees the same rows as the tab (`attr_d4.R`). In `metab_no_friendly`, every
  step and field equals W1's frozen value except one: the run step's `files`. Every kernel call and every project file
  (names and table digests) equals W1's. What remains is two store-entity counts:
  - `omics-figure` 1 and `omics-figure-asset` 1 are added;
  - `dia-state` goes from 13 to 15.

### The store-count difference that remains

The mechanism W3 recorded for lipid, unchanged for metab:

- **The entity:** one `omics-figure` product and its `omics-figure-asset` PNG,
  `publication_graphs/NumSigDeMolecules/metabolites_num_sig_barplot.png`, written through `writeWorkflowOmicsFigure` by
  the exporter. It is the one figure the exporter can draw on rows without `friendly_name`.
- **Why the new path keeps it:** `writeWorkflowOmicsFigure` keeps a figure in the store only when it finds an owner
  (`workflowOmicsFigureOwner`): the active state scope first, then the session's registered owners.
  `runWorkflowOmicsPublicationInProcess` runs the export inside `withWorkflowStateScope(manager, ...)`, so the owner is
  always found. In this case the owner was made by `omicsTestOwner` in a separate owner session and is not registered on
  the module's session, so the old exporter, outside any scope, found none and wrote the PNG only.
- **The `dia-state` +2:** the two `rootWorkflowProduct` records that root the product and its asset under the head.
- **In the re-frozen value:** the same mechanism over every figure the manifest's rows allow gives `omics-figure` 11,
  `omics-figure-asset` 19 and `dia-state` 13 + 30 = 43, the figure counts `metab_managed`'s run step records.
- **`metab_main`, the realistic case, shows no difference.** It passes identically on the W4 tree at every step and
  field, store counts included: its owner is adopted by the tab's session load on the module's session.

## Freezing

The five values were dumped as R code (`omicDaW1Code(value, 8L)`) from the W4 tree by a scratch script (`dump.R`) in a
test context, in two separate processes. The two dumps were byte-identical (`cmp`; 6,229 bytes), and the values were
pasted from the first, each with its comment.

## Result

`test-omic-da-module.R` on the W4 tree after the re-freeze (`run1.R`, `load_all`, strict skips, `TMPDIR` on disk):
`COUNTS test-omic-da-module.R pass 1234 fail 0 error 0 skip 0` (W1's 1218 plus the W4 blocks' 16).
