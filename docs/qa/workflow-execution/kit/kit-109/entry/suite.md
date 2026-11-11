# KIT-109 entry: full suite at `35769a54` (2026-11-12)

The whole of `tests/testthat` (404 files), one process per file, six in parallel, 1800 s per file (`tools/testing/test-baseline.R suite`), in both load modes, on `35769a54` (KIT-115 and KIT-116 merged, before any KIT-109 change). Build `browser-smoke.yhcfUf`, from kit_921's library (`kit/scripts-k109/verify.sh`, a copy of the remediation's `verify.sh` with that library and these paths). As at KIT-105 and the remediation, check mode's aggregate stops at its header step (the check directory is not a git repository) after the per-file table is written.

`R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs (the series baseline).

| Mode | Files | Passed | Failed | Error blocks | No result |
|---|---|---|---|---|---|
| `load_all`, remediation (`cc8764a4`) | 404 | 35,793 | 1 | 2 | 1 |
| `load_all`, this run | 404 | 34,673 | 4 | 2 | 2 |
| check, remediation | 404 | 29,106 | 149 | 237 | |
| check, this run | 404 | 29,587 | 152 | 237 | |

## Every difference (`kit-105/skipcmp.py`, per file, against `gui-audit-20261111/entry/`)

- **`test-omic-qc-module.R`, 3 failures in both modes:** two frozen DIA_LIMPA QC UI digests changed by KIT-116's P-02 label (`08e9b25e`), a file KIT-116 did not run. Shown by running the file at `08e9b25e~1` (158/0) and `08e9b25e` (155/3). Re-frozen with that attribution in `347e5dd5`; 157/0 after KIT-109 W1.
- **`test-prot-norm-core.R`, `load_all`, no result:** the 1800 s per-file ceiling under the six-way load, as the remediation recorded for it in check mode. Its last run alone on the same code (KIT-116 F116-01, `5b3a93c9`) is 1690/0. Check mode ran it (one hidden block more than the remediation's row, the A111-05 guard block).
- **Unchanged from the remediation:** the stale `> 100` count in `test-omics-module-globals.R` (fixed in W6, `4c56d2f0`); `test-gui-entry-proteomics.R` over the ceiling (no result in both runs; run alone at the exit); the two diagnostic error blocks; the 149 check-mode failures and three check-mode file errors.
- No new file-level skip, no file with fewer blocks run apart from the timeout above, no new file error.

**Decision:** the entry baseline holds. The one regression in it (the QC digests) belongs to KIT-116 and is fixed. The removals start from this table.

Tables: `skips-load_all.tsv`, `skips-check.tsv`, `header-load_all.tsv` (this directory); `check.log`.
