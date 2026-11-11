# KIT-109 exit: full suite at `63ce6d6c` (2026-11-12)

Same protocol as the entry ([../entry/suite.md](../entry/suite.md)): `tests/testthat` (404 files), one process per file, six in parallel, 1800 s per file, both load modes, build `browser-smoke.iZxvdt` (kit_922) at `63ce6d6c`, all of W1 to W7 in. Check mode's aggregate again stops at its header step after the per-file table is written.

`R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs (the series baseline).

| Mode | Files | Passed | Failed | Error blocks | No result |
|---|---|---|---|---|---|
| `load_all`, entry (`35769a54`) | 404 | 34,673 | 4 | 2 | 2 |
| `load_all`, exit | 404 | 36,582 | 0 | 2 | 0 |
| check, entry | 404 | 29,587 | 152 | 237 | |
| check, exit | 404 | 31,465 | 154 | 237 | |

The pass counts rise mainly because the two files that hit the ceiling at the entry (`test-prot-norm-core.R`, `test-gui-entry-proteomics.R`) finished this time.

## Every difference (`kit-105/skipcmp.py`, per file, against `../entry/`)

- **`test-omic-qc-module.R`:** 3 failures to 0 in both modes (the re-freeze `347e5dd5`).
- **`test-omics-module-globals.R`:** `load_all` 1 failure to 0 (W6). **Check mode 1 to 6: a W6 regression.** The installed tests have no `R/` tree, so each of W6's six per-family expectations failed where the old count failed once. Corrected in `a15d64c0` (the families are one expectation): 14/0 under `load_all`, and one failure in check mode, as at the entry. Check mode was not rerun for it, so the check row above still counts 6.
- **`test-omic-gui-config.R`:** 2 fewer passes in both modes (W1: the two expectations on `trxOmicProfile()` itself, which no longer exists. W3's two `expect_null` checks replaced two value checks one for one).
- **`test-omics-summary-core.R`:** 4 fewer passes (`load_all`) and 2 fewer (check) (W2: the two deleted functions' expectations).
- **`test-prot-norm-core.R`, check mode:** one skip, `skip_other`. At the entry this file hit the ceiling and gave no result. The skip is the file's check-mode guard (line 166, "source references are kept under load_all only"), not a new skip.
- No new file error, no new file skip, no file with fewer blocks run.

**Decision:** no regression stands. The one regression (check-mode globals count) is corrected by `a15d64c0`.

Tables: `skips-load_all.tsv`, `skips-check.tsv`, `header-load_all.tsv`; `check.log`.
