# Merge train 22 of 2026-11-26: KIT-265

**This train is NOT stacked.** It is branched from `origin/v0.9-staging` at
`176f9d3ffa100f0480fd5100a32c002c5a5cdbe5` (fetched before building; no code merge since 176f9d3f). Trains 19, 20 and 21
(`merge-train-20261126l`, `m`, `n`) are other agents' and are still out. If they land first, this train lands by a merge
whose only conflicts are `NEWS.md` and the ledger (the audit's `git merge-tree`: no code conflict against any of them).
Because `kit-265` already contains staging 15d86bd4, staging is the right base.

## Train

Base: `176f9d3ffa100f0480fd5100a32c002c5a5cdbe5`.

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-265 | `origin/kit-265` b94e62bb | 75810a2f | `VAL-330` (replaces `VAL-<merge>`) | AUDIT-KIT-265-impl-20261126 with its closure check (verdict ready) |

Head (bound by the run): `75810a2f8e22df7875775153aa0407190ed0da7d`. This file is in the declaration commit after it,
under `docs/` only.

### Conflicts met and their resolutions

- `docs/methods-decision-ledger.md` only: staging's VAL-329 kept, then KIT-265's row as `VAL-330` with the closure
  check's C1 edits (the notice clause; correction pass 45711795705 in the Status cell), `<date>` = 2026-11-26,
  `<t>` = 22. The train's run id in the Status cell is filled at landing.
- `NEWS.md` merged automatically (no conflict).
- C2: `wording.md` "(48 words." became "(58 words; longer than the shipped text: for the wording pass to shorten.".
  `kit-265/outcome.md` line 92 now also names `VAL-330`.
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`.

### Path lengths

Every path new since the base is at most 87 bytes with the `MultiScholaR/` prefix
(`docs/qa/workflow-execution/kit/kit-265/r3-45711795705/compare-failures.txt`); under 100.

## Declaration (written before the run is sent)

**Purpose.** KIT-265's merge check on the merged tree.

**Scope recount.** The ticket's grep gives 9 test files on the base 176f9d3f and 10 on the head (the new
`test-kit265-organism-routes.R`): no change from the ticket's 9 plus the new file. Scope = **10 files** + R CMD check
on both sides.

**Settings.** 1 shard per side; no env switches; no extra packages; reference BLAS; check on both sides; no scripts.
Ceiling 120 job-minutes (= 8.0 core-hours); expected about 55 (R2 56.58, correction run 53.30).

### Predictions (from the closure check, C-6)

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit265-organism-routes.R | absent | 129/0/0/0 | only-on-head |
| test-kit250-organism-dialog.R | 119/0 | 119/0 | same |
| test-prot-01c-import-module-contracts.R | 437/0 | 438/0 | count-changed |
| test-e2e-test-mode-integration.R | 61/0 | 61/0 | same |
| test-module-ci-prot-import.R | 107/0 | 107/0 | same |
| test-prot-01d-import-module-characterization.R | 79/0 | 79/0 | same |
| test-prot-import-core.R | 38/0 | 38/0 | same |
| test-workflow-gui-memory-lifetimes.R | 81/0, 11 skipped | 81/0, 11 skipped | same |
| test-workflow-server-shared.R | 226/0 | 226/0 | same |
| test-workflow-ui-shared.R | 178/0 | 178/0 | same |
| R CMD check | | | identical findings, none only on head |

### Decision and attempts

- Every prediction holds: TRAIN GREEN.
- Any other verdict is a result: recorded as it is, not re-sent.
- One run. An HTTP 500/502 at dispatch with no run created: started once by hand, recorded.

## What the coordinator resolves at landing

- The ledger: `VAL-330` takes its final number if trains 19 to 21 land first; its `<id>` fills; VAL-291 and VAL-304
  Status cells and the other C-7 records (ticket closure lines, plan row, candidates) as the closure check lists them.
- `NEWS.md` and the ledger are the only expected conflicts against trains 19 to 21.
