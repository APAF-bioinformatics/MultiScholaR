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

## Acceptance

Met: F2, F5 to F8, F10, the dead code, figures byte-identical, journeys exact, no numeric result changed. Not met: E1's memory aim on the DIA full journey. The full suite, `R CMD check` and the e2e lanes for the final tree are in [final/suite.md](final/suite.md).
