# KIT-122: proteomics normalisation figures on the worker, and the tab's leftovers (2026-11-16)

Scope: the executive call of 2026-11-15 (audit E1, F2, F3, F5 to F8, F10; [../gui-journey-audit-20261115.md](../gui-journey-audit-20261115.md)). Implemented in a separate worktree, brought onto `performance-v2` as eight commits, plus one correction.

**Result: every item landed and the journeys are exact, but E1's memory aim is not met.** The figures moved to the worker were not the cost; the composite figure is, and it is still built in the app (parked row 17).

| Item | Change | Commit |
|---|---|---|
| F2 | The RUV summary reads mode and grouping from the run's saved product each time the result changes (`resolveDIARuvRunDescription`). | `3a9f59f9` |
| F10 | The reset clears the state the summaries render, not the outputs. | `2431bfec` |
| F5 | The RUV tab's single-result content is a library entry (`omicNormRuvSingleResult`, `norm$ruv_single_result`). | `a1e34950` |
| F6 | `QC_composite_figure` and `qc_plots` leave the tab's state fields. | `46beeb49` |
| F7 | The module's start line names the tab's omic. | `963b0846` |
| F8 | The inert design-choices observer and its no-op are gone. | `fba00242` |
| Dead code | Seven unreferenced functions removed. Three of the listed handlers are kept: the profiling diagnostic names them. | `4c9a5690` |
| E1 / F3 | A new worker operation, `normalization-qc-figures`, saves the four stage figures (products and images) and the correlation tab's PCA. The tab's generators then write the saved bytes. With no runner they compute as before. | `5bb57f2e` |
| Correction | The tab-entry and aesthetics dispatch declines while another module's action holds the application, and the tab draws in process. | `01fe465e` |

## Journeys

- **kit_937** (`587f92b5`): 937m, a, b, c accepted. 937q (DIA full) ends in its known `oracle_mismatch`; against 935q, store against store: 16 states, 105 comparisons, 0 different (`compare-935q-937q.tsv`).
- **937l and 937t (LFQ, TMT) failed.** The driver opened the normalisation tab while the last QC step's action was still active. The new figures dispatch was refused ("another heavy action is active") and the journey read the notice as the normalisation's failure. A user who opens the tab that fast would have seen the same notice.
- **The one correction pass, kit_938** (`522ad32c`): 938l and 938t accepted.

## Memory: the measure the ticket names

App memory in `normalisation_ruv`, one run each (`app-memory-by-step.tsv`):

| Journey | Before | After |
|---|---|---|
| DIA full | 0.88 to 1.79 GiB (934q) | 0.91 to 1.84 GiB (937q) |
| LFQ screen | 0.57 to 1.15 GiB (936l) | 0.57 to 1.02 GiB (938l) |

- DIA full is unchanged. The whole rise falls in the step's last 20 s, after the workers finish, where the app builds the composite figure and the filtering-progress grid.
- The composite alone, built from 937q's 13 images in a fresh R process, takes the process from 0.37 to 1.29 GiB. It also writes a 29.6 MB `.rds` of the plot.
- The LFQ difference is one run against one run and is not a finding.

## Byte identity (tests, `test-prot-norm-core.R`)

- Managed equals unmanaged: every stage image's md5 after every step on the DIA fixture (six steps) and the DIA_LIMPA fixture (three steps), and the correlation PCA's data and title.
- Managed, the app saves no product and no image; the worker saves 16 products and 24 images.
- KIT-111 W1's frozen captures are unchanged. Re-frozen with attribution: the managed drives' operation lists, working-notice counts and action positions. The 116 `ruv_qc_ui` lines of the five proteomics entry oracles are deleted (F5: the output no longer exists).

## Tests (passed / failed / skipped on the final tree)

`test-prot-norm-core.R` 1,872/0/0; `test-omic-norm-module.R` 639/0/0; `test-prot-05b-norm-module-contracts.R` 439/0/22; `test-module-ci-prot-norm.R` 60/0/4. Quarantine entries Q-187 (eight blocks on the deleted functions) and Q-188 (20 assertions removed from live blocks).

## Noticed, left untouched

- More worker round trips per click: normalisation +1 or +2, correlation +1, tab entry +1, an aesthetics change up to +3.
- Controller-side identity lookups do not see entities a worker saved until the store is reopened. It is fixed at figures delivery only (`indexWorkflowHeadProducts`).
- By reading only, not reproduced: a managed retry after a failed publication may be refused (`assertDIAExportReady` with no `allow_action`).
- DIA_LIMPA gets no pre-normalisation QC on tab entry. With RUV skipped the RUV summary still reads "Manual RUV Parameters" with NA values.

## Final tree

Full suite in both modes, `R CMD check` and seven e2e lanes at `ed99a0a8`: no regression ([final/suite.md](final/suite.md)). `test-prot-norm-core.R` now runs for about 1,600 s and reached the suite's 1,800 s per-file ceiling once in check mode (it passes alone). The two drives this ticket added to it should move to a file of their own.

## Acceptance

Met: F2, F5 to F8, F10, the dead code, figures byte-identical, journeys exact, no numeric result changed. Not met: E1's memory aim on the DIA full journey.

## Revision 1.1 (2026-11-16): the composite figure is built by the worker

E1's memory aim was missed because the app still built the composite QC figure (0.37 to 1.29 GiB alone). Revision 1.1 moves it into the figures operation (`figures = "composite"`: `saveWorkflowDiaCompositeFigure`, `R/utils_workflow_dia_products.R`); the app reads the saved image. Commits `134af313`, `23148110`, `4e5320ae`, `296367e8`, the test patch after them, `561eb2b1` (kit_939 = `rYf7jA`).

**App memory in `normalisation_ruv`, one run each:**

| Journey | Before revision 1.1 | After |
|---|---|---|
| DIA full | 0.91 to 1.84 GiB (937q) | 0.90 to 0.94 GiB (939q) |
| LFQ screen | 0.57 to 1.02 GiB (938l) | 0.58 to 0.58 GiB (939l) |

The exit of revision 1.1 (a rise under 0.3 GiB on the DIA full journey) is met.

**Journeys (kit_939):** 939l (LFQ), 939t (TMT) and 939m (DIA_LIMPA) accepted, exact. 939q (DIA full) ends `oracle_mismatch` as that route always does and is judged store against store with 935q: 104 of 105 fields identical ([compare-935q-939q.tsv](compare-935q-939q.tsv)); the one difference is the head's rendered page count, which is the composite's three assets (PNG, PDF, RDS) now saved by the worker. No numeric result differs.

- 939m's first attempt was `invalid` (`measurement_valid`: a sampling gap while simulator pilots loaded the machine to 18). It is kept (`runtime-keep/kit-journey-939m-invalid-under-load`); the rerun on a quiet machine is the one counted.
- 939q's decision file was written after the member script's 4,000 s wait ended; its records were written from the finished stage afterwards.

**Store growth:** each managed normalisation now roots three more rendered assets per composite identity (keyed on the product entity and the sha256 of the RUV mode and stage images).

**Tests:** the worker-figures drives are in `test-prot-norm-worker-figures.R` (63 passed); `test-prot-norm-core.R` 1,854 passed alone and is off the suite's time ceiling in check mode.

**Final tree:** full suite in both modes, `R CMD check` and seven e2e lanes at `ea236689`: no regression ([final2/suite.md](final2/suite.md)).

**Acceptance, revised:** every criterion met, E1's memory aim included. Parked row 17 is closed.
