# Merge train 7 of 2026-11-25: KIT-252, KIT-250

This train follows the coordinator's instruction after train 6 landed (staging c4505c7f, ledger rows through VAL-261).
The procedure is the same as train 6 (`merge-train-20261125f.md`).

## Train

Base: `c4505c7f0fdda79692f1baaafa636b2c710d7a19` (`origin/v0.9-staging`, fetched before building).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-252 | `origin/kit-252` bd87aa17 | 13c295cb | `VAL-<train7-1>` | AUDIT-KIT-252-impl-20261125 (section 6; I07 status form) |
| 2 | KIT-250 | `origin/kit-250` aee51510 | fbf42a9c | `VAL-<train7-2>` | AUDIT-KIT-250-impl-20261125 (closure check) |

The run binds head `fbf42a9c` (full hash in the request). This file and the records commit after the run touch only
`docs/qa/`. The two ledger rows are numbered by the coordinator at landing.

### Conflicts met

| File | Merges | Resolution |
|---|---|---|
| `docs/methods-decision-ledger.md` | both | Staging's rows through VAL-261 kept. Each member's row follows as `VAL-<train7-k>` with a clean status cell. Each `outcome.md` placeholder replaced to match. |
| `NEWS.md` | both | Each member's entry against the one before it, at the top of "### Fixes" under the one unreleased heading. Every entry kept, earlier members first. |

There was no conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`. `DESCRIPTION`, `NAMESPACE`
and `man/` are unchanged.

### Path lengths

Checked by script before sending: 112 paths are new since the base. The longest is 99 bytes with the `MultiScholaR/`
prefix (`kit-250/h1-45611275368/g0-h1a-reimport_dialog_open.json`). None exceeds 100.

## Declaration (written before the run is sent)

- **Purpose:** each member's check on the merged tree, done once on the train.
- **Scope:** 22 files, base and head, with R CMD check on both sides and `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` on both sides.
  The run uses 4 shards per side, with no extra packages.
  - KIT-252 (15 files):
    - `test-kit252-managed-failure-text.R`;
    - the plan's 13: `test-prot-norm-core.R`, `test-prot-norm-input-current.R`, `test-omic-da-module.R`,
      `test-workflow-omics-application-reads.R`, `test-workflow-execution-process.R`, `test-workflow-execution.R`,
      `test-perf-107.R`, `test-perf-003.R`, `test-perf107-audit.R`, `test-prot-07b-da-handlers-compat.R`,
      `test-workflow-core-harness.R`, `test-workflow-kit-harness.R`, `test-workflow-execution-adapters.R`;
    - `test-kit251-replicate-filter-one-click.R`, because KIT-251 is now on staging (register section 6).
  - KIT-250 (7 files, closure check): `test-kit250-organism-dialog.R`, `test-prot-01c-import-module-contracts.R`,
    `test-prot-01d-import-module-characterization.R`, `test-module-ci-prot-import.R`, `test-module-ci-trx-app.R`,
    `test-workflow-gui-memory-lifetimes.R`, `test-workflow-server-shared.R`.
- **Ceiling:** 1,200 job-minutes, the sum of the members' merge-check ceilings (KIT-252's R2 600 and KIT-250's R2 600,
  each ticket section 10).
- **Base values:** read from the base as it now is, that is from the latest train run whose tree equals c4505c7f in
  that file's code:
  - train 6, 45608483144: `test-prot-norm-core.R` 1548/306 (already with KIT-224's two digests);
    `test-perf107-audit.R` 191/0/1 skipped (KIT-251's edit); `test-kit251-replicate-filter-one-click.R` 9/0;
    `test-prot-norm-input-current.R` 321/0; `test-workflow-omics-application-reads.R` 230/0;
    `test-workflow-kit-harness.R` 734/1.
  - train 4, 45586240983: `test-omic-da-module.R` 1316/53; `test-prot-07b-da-handlers-compat.R` 35/0/2;
    `test-module-ci-prot-import.R` 107/0; `test-module-ci-trx-app.R` 50/0; `test-workflow-gui-memory-lifetimes.R`
    81/0/11.
  - Files no train run has measured since are predicted `same`, base against head, with the members' R2 values as
    the expected counts.

### Predictions, KIT-252 (register section 6 and the implementer's report)

- `test-kit252-managed-failure-text.R`: only on head, 72 passed, 0 failed.
- Every other file `same`:
  - `test-prot-norm-core.R` (1548/306) and `test-omic-da-module.R` (1316/53): exactly the base's failures, test by test.
    The A111-03 and KIT-112 W3 blocks pass on both sides.
  - `test-workflow-execution.R`: 114/0, 10 skipped.
  - `test-workflow-omics-application-reads.R`: 230/0.
  - `test-workflow-kit-harness.R`: 734/1.
  - `test-perf107-audit.R`: 191/0, 1 skipped, KIT-251's numbers.
  - `test-kit251-replicate-filter-one-click.R`: 9/0.
  - `test-workflow-execution-process.R` 191/0, `test-perf-107.R` 21/0, `test-perf-003.R` 137/0/1,
    `test-prot-07b-da-handlers-compat.R` 35/0/2, `test-workflow-core-harness.R` 131/0,
    `test-workflow-execution-adapters.R` 0/0/6 skipped and `test-prot-norm-input-current.R` 321/0, as in R2.

### Predictions, KIT-250 (closure check; any other verdict for these seven is a stop)

| File | Head (passed/failed) | Verdict |
|---|---|---|
| `test-kit250-organism-dialog.R` | 119/0 | only on head |
| `test-prot-01c-import-module-contracts.R` | 437/0 | `same` |
| `test-prot-01d-import-module-characterization.R` | 79/0 | `same` |
| `test-module-ci-prot-import.R` | 107/0 | `same` |
| `test-module-ci-trx-app.R` | 50/0 | `same` |
| `test-workflow-gui-memory-lifetimes.R` | 81/0, 11 skipped | `same` |
| `test-workflow-server-shared.R` | 226/0 | `same` |

### Predictions, R CMD check

- Findings identical on both sides, with nothing only on head.
- The portable-file-names NOTE lists the same number of paths on both sides. No `kit-250` or `kit-252` path is in it,
  since every new path is at most 99 bytes.

### Decision and attempts

- If every prediction holds, the train is GREEN.
- A member whose prediction fails is named with its files and counts, and is not repaired here.
- One run is sent. One re-send is allowed, only for an instrument defect where no result was judged.

## Result

(filled after the run)
