# KIT-120: proteomics' normalisation tab on the shared frame (2026-11-12)

**Outcome:** done except one open finding (memory on the DIA full journey, under Exit runs). The tab runs on the shared frame for all five proteomics routes. No numeric result, file product, input id or harness file changed. All counts below are `load_all`, one process per file, unless a row says otherwise.

## What the join is

Proteomics' normalisation tab is `mod_omic_norm_ui`/`mod_omic_norm_server` under `omicGuiConfig("proteomics")$norm`. The shared frame builds the whole display: the option column, the QC tabset, the image grid and the correlation tab. Proteomics' actions (run, RUV, correlation, covariates, export, reset, restore, the S4 download) are its existing handlers, registered on the frame's state by `protNormFrameServer`.

**The handlers were not unified with the frame's own pipeline.** Their field names and file products are what the worker and the journeys are exact on; merging them is a second, larger change ([sites.md](sites.md), section 3). So the join retires about 760 lines of display code, and about 2,700 lines of handlers stay.

## Work items

| W | Item | Result | Commits | Evidence |
|---|---|---|---|---|
| W0 | Sites and map | [sites.md](sites.md). Two of its findings were wrong and are corrected: the "RUV grouping default" (a misreading: both frames select `group`), and the member list (the DIA `full` member covers the DIA route, but by a before-and-after pair: its verifier decision is KIT-065's known mismatch set on every build). | `33ab4c55`, `c58e22d7` | reading |
| W1 | Contract | The `norm` contract gains: an omic's own server (by name, as the DA plugins are); the Post-Normalisation QC tab; an S4 download section; a section's own texts; the stage grid's label and file-name options and QC directory; a slot under the correlation summary; a manual k that may be unset. With `server`, the frame's own run, correlation and export fields are not required. | `6b1011b3` | `test-omic-gui-config.R` 1892/0 (new blocks failed first) |
| W3 | Frame entries | The library entries and the server hook for those fields. Lipid's, metab's and trx's frozen UI and captures are unchanged. | `fa60b6a5` | `test-omic-norm-module.R` 592/0 (new blocks failed first); lipid, metab and trx entry lanes 35/0, 105/0, 29/0 |
| W2 | Proteomics' group | `omicGuiProteomicsNorm()` and `protNormFrameServer()`, not yet mounted. The frame's proteomics UI has every input of the tab's own UI with the same starting value (frozen in a test). | `7a0b3b60`, `c668802d` | `test-omic-gui-config.R` 1906/0; ten proteomics normalisation files and the entry lane at their counts. **A missed file:** `test-prot-norm-core.R` held KIT-111's assertion that the group is NULL and was not in W2's run; corrected in `c668802d`. |
| W4 | The switch | `mod_prot_norm_ui` and `mod_prot_norm_server` are thin wrappers over the frame, with unchanged exported signatures. **Stop rule not triggered:** the five re-frozen oracles differ from the stored ones only by the twelve plot ids (1,392 removed, 1,392 added), the added `ruv_qc_ui` container (116) and the oracle's date and commit. No state, table, file hash, notification or other output differs. Each plot slot shows a plot exactly when it did (444 shown, 948 empty). | `dc2a05b2`; oracles `c18cb435` | `test-gui-entry-proteomics.R` 193/0 on the new oracles; 31 dependent files at their counts; Q-184 |
| W5 | Removals | The three server shells, the image payload, renderer and registrar, and the state constructor (284 lines; with W4's UI body, about 760). | `9a16bfa3` | dependents equal the post-switch results except the quarantined blocks; Q-185, Q-186 |
| W6 | Browser lanes | No selector needed changing: no browser test named an old plot id. | none | the e2e lanes below |
| W7 | Records | This file, [sites.md](sites.md), [exit/suite.md](exit/suite.md), the quarantine register, P-32. | | |

Order: W3 was done before W2, because proteomics' group names its server function, which plugs into the hook W3 adds.

## What a user sees

- The plot outputs and the RUV QC container have the shared ids (tests and custom CSS only).
- The option column has the library's order (Plot Aesthetics above Normalization Method), heading levels and icons. Proteomics' help texts, labels, inputs and defaults are its own, unchanged.
- The correlation tab has the shared wording. Three help sentences differ slightly from proteomics' former ones (for example "will be flagged" for "will be removed from the analysis").
- No Log tab and no assay label row (user, 2026-11-12).

## Quarantined tests

Nine blocks, each in `quarantined-tests.md` with what covers it now: Q-184 (the server seam, W4), Q-185 (three shell blocks, W5), Q-186 (five image and state blocks, W5).

## Left as it is, recorded

- `updateProtNormDesignDrivenChoices` and `getProtNormPlotAesthetics` are no longer called by the application. Their tests sit in blocks that also cover live helpers. `protNormFrameServer` passes the registrar a no-op for the first, since the frame registers the design's choices. Removing both is a small follow-up.
- Eight handler functions were already unreferenced in `R/` before this ticket (`applyProtNormRuvCorrectionStep`, `finalizeProtNormRuvCleanupStep`, `prepareProtNormNormalizationRun`, `runProtNormBetweenSamplesStep`, `runProtNormPostNormalizationQcStep`, `applyProtNormSkippedRuvState`, `resolveProtNormStep6QcObject`, `runProtNormStep6RuvQc`). Not touched.
- The step map's line references into `mod_prot_norm_*` files (`kit-001/step-map.json`) may have shifted. That is KIT-117's kind of correction.

## Exit runs

- **Full suite, both modes, at `9a16bfa3`** ([exit/suite.md](exit/suite.md)): `load_all` 36,769 passed, 0 failed; check mode 31,710 passed, 149 failed (the known check-mode set). No regression against KIT-109's exit.
- **`R CMD check --no-tests`:** 3 WARNINGs, 3 NOTEs (baseline). **Install:** kit_923 (`browser-smoke.KaJpsV`).
- **e2e lanes** ([exit/e2e.txt](exit/e2e.txt)): seven of seven pass on kit_923, one attempt each, with KIT-109's counts: lipid 88/0, metab LC-GC 190/0, trx screen 16/0, proteomics DIA 42/0, DIA_LIMPA 86/0, enrichment 98/0, TMT/LFQ 153/0. The four proteomics lanes drive the switched tab in a browser.
- **Module CI `--omic all`** ([exit/module-ci-all.log](exit/module-ci-all.log)): passed, 110 scenarios.
- **Journeys** ([exit/journeys/](exit/journeys/)), one attempt each:

| Member | Journey | Build | Verifier | Evidence |
|---|---|---|---|---|
| 923m | DIA_LIMPA, covariates | kit_923 | accepted | exact against `limpa65_screen` |
| 923l | protein, LFQ | kit_923 | accepted | exact against `lfq101_fixture` |
| 923t | protein, TMT | kit_923 | accepted | exact against `tmt101_fixture` |
| 922r | DIA_LIMPA, correction | kit_922 (before the ticket) | rejected, `oracle_mismatch` | the oracle is kit_056's (`limpa_corr_screen`), never re-recorded |
| 923r | DIA_LIMPA, correction | kit_923 | rejected, `oracle_mismatch` | **store against store with 922r: 11 states, 75 comparisons, 0 different** (`compare-922r-923r.txt`) |
| 922q, 923q | DIA (IQ), full | kit_922, kit_923 | both rejected, `oracle_mismatch` (KIT-065's known set, on every build) | **store against store: 16 states, 105 comparisons, 0 different** (`compare-922q-923q.txt`) |

- **The correction journey's mismatch predates this ticket:** it is the same on the build before it, and the two builds' stores are identical. **Finding for KIT-118:** the correction-mode oracle is stale; that journey has had no exact reference since kit_056. Which step first differs was not established.
- **Memory, DIA_LIMPA (922r against 923r, the same journey before and after):** the normalisation steps' slice peaks are `protein_normalisation_limpa` 1.53 and 1.53 GiB, `correlation` 1.56 and 1.57, `ruv_covariates` 1.75 and 1.76; whole member 1.94 and 1.94. No change.
- **Memory, DIA (IQ) full (922q against 923q): OPEN FINDING.** 923q's application process sits 0.3 to 0.54 GiB above 922q's from the `design` step to the end (`exit/journeys/memory-922q-923q.txt`). The declared rule was "no worse beyond the recorded member-to-member spread"; the recorded spread is up to 0.37 GiB (KIT-112, 111q against 112q, on steps no ticket changed), so this is outside it and is not accepted on one pair. What the traces show:
  - **It starts before the normalisation tab is used.** `import` is equal (peak 2.316 and 2.300 GiB). During `design` 922q's application falls to 1.323 GiB and 923q's only to 2.048: the release after the import did not happen in 923q. From there the offset is constant to the last step.
  - **The normalisation step's own rise is the same:** 922q 1.036 to 2.084 GiB (+1.05), 923q 1.574 to 2.547 (+0.97).
  - **The screen pair shows no offset at all:** 922r against 923r, the application is within 0.016 GiB at every step (`exit/journeys/memory-922r-923r.txt`); 923m's `rollup` minimum is 0.729 GiB, in the released group of P-34.
  - **Reading:** this has the form of the run-level release floor of P-34 and KIT-112 (an import-time floor that is sometimes not released, on code no ticket changed), not of a cost in the tab. One pair cannot separate the two on the large corpus. **Decision needed from the user:** one more 923q run (about an hour). If its application releases after the import as 922q's did, the offset is run-level and the finding closes; if not, it is this ticket's and is diagnosed before closure.

## Acceptance

- [x] W0 recorded; the three decisions and the three later ones answered.
- [x] `omicGuiConfig("proteomics")$norm` is valid; the tab is the shared module for the five routes.
- [x] Every proteomics entry-lane data capture identical; the UI ids re-frozen once, attributed.
- [x] Lipid, metab and trx normalisation captures unchanged.
- [x] Proteomics journeys: three exact against their oracles; the correction and DIA full journeys identical store against store to the build before the ticket. Harness files unchanged.
- [ ] **Memory no worse: open** on the DIA full pair (above). Level on the DIA_LIMPA pair.
- [x] e2e lanes, module CI, `R CMD check` at baseline, install, full suite in both modes.
- [x] `mod_prot_norm_ui` and `mod_prot_norm_server` exported and working as wrappers.
- [ ] The series closure (28 of 28 cells) is stated when the memory finding is settled.
