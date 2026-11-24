# KIT-159 outcome, builder B (replay, the entry script, G1(a)'s job pieces, E4)

Branch `kit-159-replay` from `kit-159-validation-harness` at b8683b3d, worktree
`/home/doktersmol/Projects/APAF/MultiScholaR-kit-159-b`. Builder A works on the ticket branch at
the same time (coverage and results code, the main records, one diagnostic cell); this file and
`runs-b.tsv` are builder B's records and are merged at the end. Budget: builder B's 6 of the
REPAIR line's 25 core-hours (14.81 used before builder B; builder A shares the rest).

## Items

| Item | Register | What was built | Tests |
|---|---|---|---|
| 1. Replay of declared cells (G5) | I07, E5 | `R/replay.R`: the cell writes its analysis record inside its session (`writeWorkflowAnalysisRecord`); after the session, a copy of its store is pruned to the record's import and design states by the package's `kitPrepareReplayProject` (extracted from `tools/profiling/workflow_kit.R` with its closure, never edited), the record replayed there with `replayWorkflowAnalysisRecord` on fresh experiment paths, under the cell's thread variables and seed; the two stores compared by the package's own fingerprints (every state, as KIT-028's `kitReplayJourney`, and the head's DA result). Statuses `identical`, `differs` (the differing states and `da_results` named), `unsupported` (class); `not_declared` and `not_run` for cells not replayed. `replay = not_built` removed; `valCheckProvenance` refuses a cell without a status of this set and a `differs` without its differences | T-REP-1 to T-REP-4, T-PROV-2 |
| 2. `val-run.R` driven itself | E2 | The cell and score phases read the committed declaration through `valCheckDeclaration` first (the working tree is no longer read for it), and any `val_error` ends the job with one template line (class) and exit 3 | T-RUN-E2 |
| 3. G1(a)'s job pieces | I15 remaining, VAL-120 (2), VAL-133 (1) | `R/g1a-entry.R`: `valJobEnvironment` (CPU model, image, the ten D7 keys, lockfile SHA-256; derived as KIT-199's `writeEnvironment`), `valLockfileUnchanged`, `valG1aDaInputs` (live table, link, mode), `valG1aGatesMd`, `valParityL1Tables` (L1's DA tables for G13 read by path and checked by SHA-256); the entry script `val-g1a.R` (environment, part (1), part (2), `gates.md`) | T-G1A-ENV-1, -ENV-2, -GATES, -ENTRY, T-PAR-5; the one-lane smoke run |
| 4. E4 | E4 | `helpers_extracted` filled by `valRunPhase` whenever its caller (`val-run.R`) passes none: the GUI-entry helpers where a phase runs workflow cells, the replay kit where it replays. Per-unit truths: builder A's `truth.per_unit` (unit -> truth id -> hash) kept at the merge; builder B's parallel version removed | in T-REP-4, T-PROV-2 |

## Decisions and divergences for the owner (they change what a G5 status means)

- **G5 compares stores, not the DA table's text.** 7.3 says "the replayed DA table compared by
  masked digest, then by value". `replayWorkflowAnalysisRecord` re-runs the DA without publishing
  its tables (proteomics `publish = FALSE`; the omics DA writes no files), so there is no replayed
  DA table to read. G5 here compares the package's content fingerprints of every state and of the
  head's DA result (`snapshot$da_results`), the comparison KIT-028's replay tests use. `identical`
  is therefore at least as strict as the DA-table digest (every state and the DA result); a
  `differs` names what differs but carries no largest value difference. Writing the replayed DA
  tables would need a product change or a re-implementation of the publication; neither is in
  this ticket's scope.
- **L1's DA tables for G13 and the feature-blind lifecycle.** G13's chain needs L1's driver DA
  table per route (keys, q and fold-change columns) committed; 5.2 says no feature-level table
  leaves a lifecycle job. The reader and its hash check are built (`valParityL1Tables`); what L1
  writes, and whether its generated-input DA tables are exempt from the feature-blind rule, is the
  owner's to rule before L1a.

## Runs

See `runs-b.tsv`.

- RB0 (45335018540, e2485409, tests only): observed, every new test fails at its fault (table in
  `runs-b.tsv`). Also observed: `val-run.R` lifecycle, counted and score ran to completion on the
  toy ticket in a temporary repository on the unchanged script; only its refusal path failed.
- Smoke (k159b-g1a-smoke-27059497, val-g1a.R part (1) on lipid_canonical) and RB1
  (k159b-rb1-1080f82a, the fix with --lint): SENT, results not yet read.

## Resume here (handed over)

1. Read the smoke run and RB1 (fetch.sh); record both in `runs-b.tsv`; one correction per item.
2. DONE early (coordinator's note): merged builder A's head 67594694 (conflicts in `phases.R` and `provenance.R` resolved to A's `truth.per_unit`; tests adapted). Still to do: re-merge A's later head if any, (builder A's head) into this branch, run the changed test
   files once (replay, records, valrun, g1a-entry, phases, repair plus A's), push, report the merge.

## Not checked

- Part (2) of G1(a) on a real lane (needs G1(r)'s `g1a-reference.tsv`; not built here).
