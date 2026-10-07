# PERF-107 deferred qualification handoff

PERF-107 is deliberately incomplete and deferred at the user's request. Architecture
performance improvements come next; revisit the harness against their final source
and installed builds. No automatic default activation or further qualification runs
are authorized by this handoff.

## Included implementation

The package changes add peptide-stage checkpoint Save/Open using the existing
archive and publication machinery, transactional restore, supported-route checks,
Design-cache cleanup, and restored numeric-input handling. The proof runner includes
bounded entry, checkpoint and synthetic-fixture cases, including a fresh-only replay
that consumes an explicitly bound retained archive. Scientific kernels and checkpoint
format remain unchanged.

## Retained verification

- Installed build `5uSiUc`, check `49dWk9`: 234 assertions passed.
- Fresh-application case `uU1Tuk`: exact history and cutoff100 restored; two continued
  Apply operations matched native outputs; continued Save and rejected-invalid-Open
  preservation passed. Complete verification141.0045s; application/browser/worker
  peak4,310,474,752 bytes (4.014GiB), zero swap. Cleanup left no managed jobs running.
- Earlier same-application checkpoint case `OXizbB` passed on build `T22M8I`.
  Its failed fresh-application case remains historical evidence, not a pass.
- Subsequent component checks against5uSiUc passed73 assertions, including rejection
  of stale designs and skipped imputation. This is component sensitivity evidence,
  not a complete installed ownership inventory.
- Synthetic calibration/heldout corpora each contain300 samples; generation and
  bounded readback were verified. This does not establish real-study applicability
  or a large-input performance result.

Raw receipts, original failures and exact input/build bindings remain in the local
ignored `docs/qa/workflow-execution/perf-107/` evidence directory and external proof
store. They are not distributed by this commit. Existing runner selectors still
require those inputs and source-bound authorities; a cloned checkout alone cannot
replay the retained archive. The current local ticket is `tickets/openai/PERF-107.md`,
with status `blocked` meaning user-deferred, not completed.

## Shelved diagnostics

[Recoverable diagnostic draft](checkpoint-ownership-diagnostics.patch) is inactive
and is not part of the executable runner or test suite. Its direct-materialization
probe passed51 checks. The following holder-inventory run had55 passing assertions
and1 fixture error: `roots$field <- NULL` deletes the required named root. Remaining
holder checks did not run. The draft has not been installed-GUI qualified.

The working runner and consumer test were restored byte-for-byte to the earlier
73-check source bindings before this commit. Review the patch against the future
architecture; do not apply it blindly. Preserve named NULL entries when repairing
the fixture, then qualify non-hydrating observation, direct assay reads, cleanup
and detector sensitivity. A zero `full_s4_reads` counter alone does not prove the
absence of full-assay materialization.

## Resume gates

1. Reconcile future source/build identities, input bindings, spent allowances and a
   finite new verification sequence. Reuse evidence only where bindings remain valid.
2. Finish installed H ownership/history coverage, F adoption/publication failure
   coverage and G interruption/guard coverage, including outstanding legacy/recovery
   obligations. Resolve or explicitly scope the unavailable cp04 LFS snapshot.
3. Only then run the planned matched large-input comparisons, followed by the
   corpus-bound crossover, heldout and routing decisions. No old timing exception
   changes PERF-107's frozen thresholds.

This commit preserves working repairs and a reviewable deferred handoff. It does
not claim whole-workflow performance qualification or ticket completion.
