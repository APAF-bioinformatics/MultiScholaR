# KIT-106 outcome: one shared DA tab for every omic; the volcano on demand (2026-11-09)

**Done.** The four omics now run their DA tab on one shared frame (`mod_omic_da.R`) and panel library (`mod_omic_da_panels.R`), configured by each omic's `da` group in `omicGuiConfig`.
- **Lipid, metab and trx** run fully on the frame. Proteomics shares the UI sections and renders; its init, load, run and archive handlers are plugins with their names and signatures kept. Proteomics' kernel placement is unchanged (KIT-111 to KIT-113).
- **The interactive volcano** is off by default and built on demand for every omic, as the user decided on 2026-11-06 (D3).
- **Metab's tree-cut controls** were already fixed by KIT-101 F-03; their parity guard held throughout.
- **Back-references:** 6 → **2**. Only the exported design-builder wrappers (#4, #5) remain.
- **Code:** across `R/`, 4,000 lines fewer (2,065 added, 6,065 deleted, against `b75561de`; [lines.tsv](lines.tsv)).

  | Area | Before | After |
  |---|---|---|
  | Lipid DA | 2,569 | 58 |
  | Metab DA | 2,483 | 56 |
  | Proteomics DA | 1,868 | 783 |
  | Trx DA | 171 | 0 |
  | Shared frame | 0 | 1,111 |
- **Builds:**
  - `kit_106` = `core/browser-smoke.m4DXz7` (`9152e613`), for runs 5 to 7;
  - the corrected build `J57eHO` (`e1a5b960`, bound as `kit_906`), for the correction pass's reruns.

## Work items

| Item | Commit(s) | What |
|---|---|---|
| W1 | `22bd937a` | The `da` group and validator. `test-omic-da-module.R`: 1,183 captures over four omics, frozen on the start commit. |
| W2 | `31062b13`, `c3a35150` | The core moves (back-references #29 to #32), the in-process volcano page factory and the DA session loader. The D4 export equivalence held. The lead wired the config. |
| W3 | `d7a9b96e`, `43517efa`, `96307945` | The shared frame and panel library, with lipid and trx on it. D3 for both. F106-05 fixed for lipid. Five W1 values moved with D3/D4 and were re-frozen with attribution (lead decision). |
| W4 | `3de87970`, `315959a0` | Metab on the frame, with its preflight through a config adapter. D3. F106-05 for metab. The twin W1 values re-frozen. |
| W5 | `f7aeeaf0`, `8ebd07d6` | Proteomics' UI sections and renders on the frame; its handlers as plugins; D3. Two `volcano_glimma` values re-frozen. The user's stop rule did not fire. |
| W6 | `0034c4aa`, `d4ca601f`, `9152e613` | The driver's volcano tick (D10, P-23; its test failed first). The wired function fields made required. Quarantine rows. Merge checks. |
| Lane re-freezes | `8420bb06`, `de484cee` | lipid_canonical and the three metab lanes, once each, for D3: exactly 6 `volcano_glimma` values per lane ([lane-changes.md](lane-changes.md)) |
| Correction | `b9820d9e`, `e30ca0b2` | F106-06: the volcano iframe is sandboxed without `allow-modals`. Test-only: a check-mode source-reference guard. |

## Proof

| Check | Result |
|---|---|
| Merge checks ([merge-checks.txt](merge-checks.txt)) | `Collate` 353/353. NAMESPACE byte-identical; `man/`, the app files and the shell unchanged. `tools/` changes only in D10's two lines. Back-references 2. The dependents grep found no live caller of the 120 deleted names or the 14 files. |
| Run 1, targeted ([targeted.txt](targeted.txt)) | 49 files: 9,441 passed, 0 failed, 0 errors, under strict skips |
| Run 2, Simulation Setup ([entry-suite.md](entry-suite.md)) | `metab_combined` and `prot_dia` alone: identical, with the expected traced calls |
| Run 3, entry suite ([entry-suite.md](entry-suite.md)) | All 10 lanes identical. Four oracles were re-frozen for D3, once each. |
| Run 4, DIA session and memory set ([dia-memory-set.txt](dia-memory-set.txt)) | 733/2, 191/1, 38/1, 174/4, 80/11, each equal |
| Run 5, `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs. Sections and code NOTE identical to KIT-110, on `9152e613` and again on `e1a5b960` after the correction. |
| Run 7, journeys ([member-results.txt](member-results.txt), [da-peaks.txt](da-peaks.txt)) | 106a to 106t: all accepted and exact. Dispatches identical to KIT-110: 106a and 106b each build `<prefix>_volcano_page` once (D10's tick). DA-phase peaks within ±2.4%. After the correction: 106x (as 106a) and 106y (as 106b) accepted and exact, dispatches identical. |
| Run 8, e2e ([e2e-final.txt](e2e-final.txt)) | On the corrected build: lipid 88, metab 190, **trx 17** (14 plus D3's 3), DIA 42, DIA_LIMPA 79, enrichment 98, TMT/LFQ 153. 0 failures, 0 skips. |
| Run 9, skip comparison ([skip-compare.md](skip-compare.md)) | Every change is a recorded quarantine, a new file, MCI-016's recorded retarget, or (check mode) one more source-tree-read error in the harness test file |

## The correction pass (used once)

**F106-06.** Run 8's trx lane hung after the volcano tick and failed after 1,378 s (pass 1, error 1).
- **Cause** ([correction/trx-e2e.md](correction/trx-e2e.md)): not the page's size (958 KB). Trx's Glimma table asks DataTables for a `GeneID` column it does not have, so the page calls `alert()` twice. In a headless browser that modal blocks the page. Users would have seen two warning dialogs.
- **Fix:** the DA frame sandboxes the volcano iframe without `allow-modals`, with a test that fails first. Page bytes, dispatches and results are unchanged, and no capture changed.
- **Reruns, once:** install and check; the two journeys that tick the volcano (106x, 106y); all e2e lanes (trx 17/0); the skip comparison.
- **The trigger itself** (F106-07, the `GeneID` display column pinned to the KIT-040 notebook) is parked as P-27.

**Instrumentation failures, recorded and rerun once:**
- Run 8's first attempt: the on-disk TMPDIR path was too long for Chrome's socket, so every lane skipped.
- The lead's first metab record run used a pinned tree whose checkout had failed.
- The correction's first bind was refused because a build key must be `kit_` plus three digits.

## Findings, decisions and parked items

- **F106-01:** the proteomics lanes never select a contrast, so proteomics' DA renders are covered by W1 and MCI-009 only.
- **F106-02, F106-03** (proteomics display): queued to KIT-111, with P-25, for one attributed proteomics re-freeze.
- **F106-04:** parked as P-26.
- **F106-05:** fixed for lipid and metab.
- **F106-06:** fixed.
- **F106-07:** parked as P-27.
- **Lead decisions:**
  - re-freeze the W1 values that D3 and D4 move, with proofs that nothing else moved (W3, W4, W5);
  - the driver tick (P-23, for review);
  - required config fields.
- **Divergences:**
  - console-only log lines;
  - the frame's registration order;
  - an export writer that returns FALSE now warns;
  - the frame keeps no locals.

  No capture changed beyond the attributed values.
- **Parked from this ticket:** P-23 (lead decision, for review), P-24, P-25, P-26, P-27.
- **Hand-offs:**
  - **KIT-107:** enrichment takes proteomics' DA hand-off, unchanged.
  - **KIT-109:** remove `profile` from the lipid wrappers, and add impact-map rules for `mod_omic_da*.R`.
  - **KIT-111:** F106-02, F106-03 and P-25 (in its ticket).

## Rollback

`git revert` of the KIT-106 commits, in reverse order. Restore `tests/testdata/gui-entry/{lipid_canonical,metab_lc,metab_gc,metab_combined}.json` from `8420bb06^` and `de484cee^`. Remove the Q-130 onwards skip lines. Revert `0034c4aa` (the driver), then regenerate the bindings.
