# GUI audit remediation: full suite at the merged tree (2026-11-11)

The whole of `tests/testthat` (404 files), one process per file, six in parallel, 1800 s per file (`tools/testing/test-baseline.R suite`), in both load modes, on `cc8764a4` (R1 to R6 and R9 merged). Build `nJSJyS` (`kit/scripts-audit-20261111/verify.sh`, from kit_098's library). `load_all` ran on a pinned worktree; check mode ran on the installed package against the `R CMD check` directory's copy of `tests/`, as KIT-105 did. The check-mode aggregate stopped at its header step (the check directory is not a git repository), as at KIT-105; the per-file table was written before it.

`R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs (as every build of the series).

| Mode | Files | Passed | Failed | Error blocks | File-level skips |
|---|---|---|---|---|---|
| `load_all`, overlay (KIT-105, then 110, 106, 107, 108, 111, 112, 113) | 401 | 37,666 | 0 | 2 | |
| `load_all`, this run | 404 | 35,793 | 1 | 2 | 32 |
| check, overlay | 401 | 27,699 | 149 | 411 | |
| check, this run | 404 | 29,106 | 149 | 237 | 37 |

The overlay is the series' latest recorded row per file (`kit-105/skipcmp.py` over the overlaid tables); KIT-105's run was the last whole-suite run before this one.

## Every difference, accounted for

- **One failure, pre-existing, a stale guard.** `test-omics-module-globals.R:25` expects more than 100 lipid and metab source files; the series deleted the duplicate copies, so the count was 90 at `ec7db20e` already (104 at KIT-105). Its global-environment assertions pass. Already listed in KIT-109's scoping ("the `> 100` count fails").
- **Three per-file timeouts under the parallel load,** rerun alone without the ceiling: `test-gui-entry-proteomics.R` 193/0 in both modes (1558 s, 1446 s); `test-prot-norm-core.R` in check mode 1603/0, 24 of 24 blocks (1549 s). The entry lanes also assert a 900 s ceiling per lane (see the remediation record).
- **Check mode, error blocks 411 to 237:** R6 (`tests/testthat/setup-here-root.R`). The checkpoint-fixture blocks now run under the installed layout. `test-prot-norm-core.R` went from 67 passes with 14 error blocks to 1603 passes with none.
- **Fewer tests in three files, all from the P-43 revert (R3):** `test-omic-enrich-module.R` 2296 to 1060 and `test-omics-enrich-core.R` 1341 to 241 (KIT-113's managed blocks deleted with the code they tested); `test-workflow-operation-registry.R` 1044 to 1040 (the file returned to its state at `132e5445`).
- **Three new files:** `test-omics-design-import-refusal.R` 39/0, `test-omics-norm-itsd-selection.R` 22/0, `test-prot-import-quantity-text.R` 313/0 (R4).
- **New file-level skips since KIT-105 (11 under `load_all`, 10 under check):** each a recorded quarantine of KIT-106 to KIT-108 (Q-133, Q-143, Q-144, Q-151 to Q-154, Q-163, Q-164, Q-170 to Q-172, Q-180, Q-181). None without a Q-entry.
- **Check-mode file errors, pre-existing and unchanged:** `test-metab-01-qc-metrics-characterization.R`, `test-metab-01e-qc-variability-helper-characterization.R` (a source-tree file read), `test-tech-reps-limma.R` (`load_all` inside the test). The 149 check-mode failures are the same count as the overlay.

**Decision:** no regression from the remediation in either mode. Not covered by this run: the commits after `cc8764a4` (A112-06's `383dd37e`, the RUV follow-up), which ran their targeted files; journeys; e2e browser lanes; any real managed worker.

Tables: `skips-load_all.tsv`, `skips-check.tsv` (this directory); `check.log`.
