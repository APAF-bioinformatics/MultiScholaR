# KIT-111 W5: the one attributed re-freeze (D7)

**Lead decision, 2026-11-10.** The values that F106-02 and F106-03 move are re-frozen once, with attribution, after W3's merge. Nothing else changed.

## How

- **The tree:** `df9c601f`, which holds W1, W2, W2b, F111-03/04, W3, and the two fixes (`5f7b6771`). The settings: R 4.6.1, `load_all`, strict skips, `R_LIBS` = `kit_108`.
- **The entry lanes:** `test-gui-entry-proteomics.R` was recorded twice in separate processes (`MULTISCHOLAR_GUI_ENTRY_RECORD_DIR`).
  - The two records differ only in each lane's wall-clock `seconds`, which is not a captured value.
  - Each run gave 131 passed and 62 failed.
  - The records were compared with each oracle (`jsondiff.py`). Apart from the header fields records and oracles carry differently (`frozen`, `load_mode`, `r_version`, `source_sha`, `omic`, `seconds`), the only differences are 62 `cluster_summary` outputs.
  - All 62 equal the oracle with the trailing `NULL` removed: 13 in `prot_dia`, 13 in `prot_dia_limpa`, 12 in `prot_lfq`, 12 in `prot_lfq_fragpipe`, 12 in `prot_tmt`.
  - The oracles were then written with `MULTISCHOLAR_GUI_ENTRY_FREEZE=1`.
- **The W1 module captures:** `test-omic-da-module.R` (`MULTISCHOLAR_KIT106_W1_DUMP`) and `test-omic-enrich-module.R` (`MULTISCHOLAR_KIT107_W1_DUMP`).
  - Each was dumped twice in separate processes, and the dumps are byte-identical (`cmp`).
  - Each dump was evaluated and compared with the file's frozen values path by path. That gives 24 differences, every one attributed:

  | File | F106-03 (`cluster_summary` loses `NULL`) | F106-02 (the run's notes gain `remove da_working`) | Other |
  |---|---|---|---|
  | `test-omic-da-module.R` | 11 | 3 (`prot_load`, `prot_main`, `prot_no_contrasts` run steps) | 0 |
  | `test-omic-enrich-module.R` | 0 | 10 (every proteomics DA run the enrichment cases drive) | 0 |

  - The counts equal the files' failure counts before the edit (14 and 10).
  - The frozen R code was edited at exactly those values: `"remove da_working"` appended after the completion notice, and `heatmap.NULL` → `heatmap.`. Each file cites this decision once, above its frozen list.

## Not changed

- **The lipid, metab and trx lanes:** their DA steps do not use these handlers.
- **Every proteomics state, product, file and digest.**
- **Notifications shown:** the lanes do not record removals (`helper-gui-entry.R:242`), so F106-02 moves no lane value. This corrects KIT-106's note, as the ticket's D7 says.
