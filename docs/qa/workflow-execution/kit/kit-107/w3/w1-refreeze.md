# KIT-107 W3: the one re-freeze of W1 (b)'s captures

`tests/testthat/test-omic-enrich-module.R`, the block between `OMIC_ENRICH_FROZEN_START` and `OMIC_ENRICH_FROZEN_END`, re-frozen
once on the W3 tree (base `eacd4a5b`), with an attribution comment in the block. W1 (a) (the UI) did not move.

## How

- Two dumps (`MULTISCHOLAR_KIT107_W1_DUMP`) in separate processes on the W3 tree, same runner and library as W1
  (`kit/scripts-078/run1.R`, `load_all`, strict skips, `core/browser-smoke.J57eHO/library`, R 4.6.1):
  both `sha256 5660b57486249a0b99973f3a7b57acd68d6dc337154411df02307c7af536e271` (byte-identical).
- The dump's case entries replaced the frozen ones. The textual diff of the block (old 2,812 lines, new 2,756) has 18 hunks,
  and every changed value is one listed below. Every other case, step and field is byte-identical to W1's freeze
  (`prefill` and every `start`, `da`, `open`, `select`, filter and download step: no hunk).
- Before: the re-frozen file on the start tree fails exactly in these cases (`test-first.txt`, "re-frozen file on the start
  tree"). After: 975/0 on the W3 tree (W1's 947 plus W3's 28).

## What moved, and why

| Fix | Case / step | Field | Start commit | W3 |
|---|---|---|---|---|
| F107-01 (W3: the tab keeps the manifest it holds) | `gprofiler/run`, `clusterprofiler/run`, `kernel_gprofiler/run`, `kernel_clusterprofiler/run`, `organism_filter no_mapping/run`, `restore/run` | `returned$da_results_data` | `<list> 8a3eda2205cf9527` (`manifest$results`, written over the holder) | `<list> 00c45cf144cf6614` (the manifest, as held before the Run) |
| F107-01 | `organism_filter mapped/run`, `refusals pre_da/run_taxon_abc`, `refusals q_cutoff/run_q_cutoff` (failed runs) | `returned` | the 13 fields with `da_results_data` = `8a3eda…` | `<unchanged>` (the failed run no longer touches the holder) |
| F107-03 (the lead's decision) | `refusals table_removed/run_table_removed` | every field the run wrote | the run completed on a NULL DA table: completion notice, progress to 1.0, `tab_status` complete, state `enrichment_completed` (61dae2e69108addd), a product, the external resource, the kernel call (`da_rows=<NULL>`), the bus values, the display tables, the pathway TSVs, the 18 returned fields | refused before any write: `[error 10 -] Error in enrichment analysis: The DA results table is missing for contrast groupB-groupA.`, progress stops at 0.3, no kernel call, and `tab_status`, history, states, product, resources, bus, display and `returned` are `<unchanged>`; no pathway TSV and no new store entry (dia-state 15, dia-table 8) |
| F107-05 (W3: the live kernel's results made portable before the commit) | `kernel_gprofiler/run` | `history`, `states`, store entries | `<unchanged>` (the commit was refused: "payload contains active or lazy bindings"); dia-state 17 | `initial, correlation_filtered, enrichment_completed`; `enrichment_completed` = `dbde2059c38a85f0`; dia-state 22 |
| F107-05 | `kernel_clusterprofiler/run` | `history`, `states`, store entries | `<unchanged>` (refused: lazy bindings, then a live readr parser pointer); dia-state 17 | `initial, correlation_filtered, enrichment_completed`; `enrichment_completed` = `9421e4861d18927e`; dia-state 21 |

Not moved by any fix: the kernel's arguments and the fake services' query and background, the product and its parameters,
the external resources, the bus values, the display tables, the pathway files and the outputs of every completed run; the
`gprofiler`/`clusterprofiler` (double) state digests (`61dae2e69108addd`, `a62a41b234214817`) and every other state digest.
The F107-05 helper is a no-op on the double's results.
