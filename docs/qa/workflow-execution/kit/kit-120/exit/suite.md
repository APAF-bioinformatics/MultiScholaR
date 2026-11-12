# KIT-120 exit: full suite at `9a16bfa3` (2026-11-12)

Same protocol as KIT-109's exit: `tests/testthat` (404 files), one process per file, six in parallel, 1800 s per file, both load modes. Build `browser-smoke.KaJpsV` (kit_923) at `9a16bfa3`, all of W1 to W5 in, installed from kit_922's library (`kit/scripts-k120/verify.sh`). The baseline is KIT-109's exit ([../../kit-109/exit/suite.md](../../kit-109/exit/suite.md), `63ce6d6c`): nothing but KIT-109's own later test correction (`a15d64c0`) and records lies between that tree and this ticket's first commit. Check mode's aggregate again stops at its header step after the per-file table is written.

`R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs (the series baseline).

| Mode | Files | Passed | Failed | Error blocks | No result |
|---|---|---|---|---|---|
| `load_all`, KIT-109 exit | 404 | 36,582 | 0 | 2 | 0 |
| `load_all`, this run | 404 | 36,769 | 0 | 2 | 0 |
| check, KIT-109 exit | 404 | 31,465 | 154 | 237 | |
| check, this run | 404 | 31,710 | 149 | 237 | |

## Every difference (`kit-105/skipcmp.py`, per file)

- **`test-prot-05b-norm-module-contracts.R`, both modes:** 509 passes to 476, skips 9 to 14. The five blocks on the deleted image payload, renderer, registrar and state constructor are quarantined (Q-186, W5).
- **`test-prot-05c-norm-module-seams.R`, `load_all`:** 53 passes to 0, skips 0 to 4. Its four blocks all tested the tab's own server shell: Q-184 (W4) and Q-185 (W5). In check mode the file was already skipped (it sources `R/` files).
- **`test-omics-module-globals.R`:** `load_all` 19 passes to 14, check mode 6 failures to 1. This is KIT-109's correction `a15d64c0` (six per-family expectations became one), made after KIT-109's exit run and not part of this ticket.
- **More passes elsewhere, no file with fewer blocks run:** the new tests of W1 to W4 (`test-omic-gui-config.R`, `test-omic-norm-module.R`) and `test-prot-norm-core.R`, 1,690 to 1,851. The last was counted per block at `c58e22d7` and at `9a16bfa3`: +1 in the adaptor block (a corrected assertion), +80 in each of the two KIT-111 blocks that compare every field of the tab's state, which now has the frame's fields too. The other 22 blocks are unchanged.
- No new failure, no new file error, no new file-level skip. The two error blocks are the diagnostic ones, as before.

**Decision:** no regression. The only tests no longer run are the nine quarantined blocks (Q-184, Q-185, Q-186), each with its register entry.

Tables: `skips-load_all.tsv`, `skips-check.tsv`, `header-load_all.tsv`; `check.log`.
