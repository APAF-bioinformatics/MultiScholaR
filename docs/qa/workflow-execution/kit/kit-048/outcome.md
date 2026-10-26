# KIT-048: the in-memory proteomics path retired

Status 2026-10-26: **in progress**. Base build `kit_044` = `core/browser-smoke.qD1Mrl`.

## Mechanism (one): deletion, plus the behaviour change that makes it possible

- **The owner from session start.** `mod_proteomics_server` creates the state owner on the project's source directory as soon as the experiment paths are known (`createDIASessionOwner`: a store with the pending identity `PROTEOMICS-pending` / `PENDING` and the placeholder `initial` state, attached to the session as the import-created owner was). The import binds the format and the workflow type to it (`bindDIASessionOwnerType`: the store marker rewritten with the identity the import-time creation used to give, the store's and the head's identity following, the manager's type set); a second import of another route is refused as before. `WorkflowState$new` no longer appears in the proteomics modules; the class stays (the adapter keeps one inside; the legacy ingest and the q-value template construct their own).
- **Branches deleted.** Every `usesDIAReferenceStates` site in the proteomics modules is gone (65 at the start of the ticket, counted on 2026-10-26; the ticket's 55 of 2026-10-19 had grown with KIT-053 and KIT-056): the apply steps keep their reference call alone; the managed-step guards lose their in-memory half; the design, DA, enrichment, normalisation, summary and rollup branches keep the owner path. The direct state calls become the owner's named entry points in `utils_workflow_dia_state.R`: `diaStateHistory` (the head's history, no hydration), `saveDIAOwnerState` (the adapter's own checkpoint of an in-process object, asserted on the owner: the protein-level design step, the enrichment results, the RUV helper, the peptide audit) and `readDIAStateObject` (the adapter's materialising read, for the few sites that need the whole object). `selectDIAStateReference` refuses a non-owner instead of reverting an in-memory manager. The fallback protein S4 creation tab (`mod_prot_qc_protein_s4.R`) is removed with its UI entry and Collate line; the legacy RDS session write (`saveProtNormExportArtifacts` without a publication request) is retired.
- **`usesDIAReferenceStates`** stays for the harness and the stage-3 code; its roxygen says no proteomics module calls it.
- **Carried from KIT-042 (F-05)**: the Norm tab re-derives `normalization_input_reference`, `normalized_protein_obj`, `ruv_normalized_obj` and `correlation_filtered_obj` from the restored head when it observes `native_session_restored` (`mod_prot_norm_observer_helpers.R`), so the references into a replaced store are overwritten there; confirmed, nothing further released.

## Acceptance greps (2026-10-26, after the deletion)

| Grep | Count |
|---|---|
| `usesDIAReferenceStates(` in `R/mod_prot_*.R`, `R/mod_proteomics.R`, `R/func_prot_qc_peptide_audit_helpers.R` | 0 |
| `WorkflowState$new(` in `R/mod_prot*.R`, `R/mod_proteomics.R` | 0 |
| `$saveState(`, `$getState(`, `$getHistory(`, `$revertToState(` in those files | 0 |
| Lines of `R/mod_prot_*.R` + `R/mod_proteomics.R` | 27,167 → 24,423 (2,744 fewer, the removed module included) |
| `DESCRIPTION` Collate | `mod_prot_qc_protein_s4.R` removed |

**Platform limit stated:** DIA's precursor rollup on the owner is the P14 native engine, x86_64 Linux only (`diaRollupNativeInfo` aborts elsewhere and no caller selects `backend = "legacy"`); the DIA and `DIA_LIMPA` routes are unsupported off x86_64 Linux until KIT-051/052 port the engine or connect the legacy backend. The LFQ and TMT routes carry no such limit.

## Tests (item 4) and a recorded deviation

- The three entry points `diaStateHistory`, `saveDIAOwnerState` and `readDIAStateObject` **forward** to the manager's commit surface rather than asserting a `DIAStateAdapter`: in the application the manager is always the owner (from session start), but the KIT-042 legacy ingest stages a legacy RDS session on a `WorkflowState` of its own before adopting it into an owner (the ticket allows that construction), and the module contract tests pass fakes exposing the same three methods. Asserting the adapter broke the legacy ingest (`buildProtDesignStateCheckpoint` on the staging manager); forwarding keeps one named entry point per operation and the grep criteria, at the cost that the entry points do not reject a non-owner. Recorded as the chosen bounded form; the ticket's item 4 (an owner fixture for the twenty mock-based module tests) is deferred: the tests that only drove the in-memory apply/revert path are quarantined (Q-22, 155 tests in 36 files; Q-21 for the removed fallback S4 tab's tests), the module tests whose logic is independent of the manager's class run unchanged on their fakes.
- The DA session validation accepts a saved reference or the `ProteinQuantitativeData` S4 a legacy RDS session carries (the KIT-042 ingest), as before; my first cut refused the S4 and the legacy adoption tests caught it.
- The server tests stub the session-owner creation point (`createDIASessionOwner`) in place of `WorkflowState$new`.

## What the GUI journeys found (2026-10-26)

The in-process checks (the module tests, the recorders, the replay) passed before the GUI journeys ran; the four screen journeys through the driver then failed at the first QC step on every route, for two defects the in-process paths cannot see:

1. **A store's identity is part of every descriptor it issues.** The first form of the owner-from-session-start rewrote the placeholder store's identity at import (the marker, the handle, the head). The placeholder's `initial` state kept the pending owner in its descriptor; the first managed worker reopened the store and refused it (`descriptor belongs to another store`). The binding now creates the bound store at import exactly as the import used to, attaches it in the placeholder's place, and removes the placeholder (`bindDIASessionOwnerType`); the test asserts the bound owner, the removed placeholder and a reopen.
2. **The Q-value compute used the module's in-memory step as its template.** `computeDIAReferenceQvalue` (KIT-014) runs the existing R method on a zero-row prototype through a private in-memory `WorkflowState` to obtain parameter resolution, validation and the output schema; it did that by calling the module's apply step, whose in-memory body this ticket deleted, so the worker recursed into the reference runner with no store (`invalid saved-entity handle`). The kernel sequence is kept as `runPeptideQvalueTemplateStep` beside the compute (a move from the module, not a rewrite); the module's step stays the reference call. No other compute leaned on a module step.

The LFQ and TMT GUI journeys completed all ten states on the fixed build; their verification against KIT-021's fixture oracles mismatched every state's summary and one identity product only: the shipped configuration gained keys since KIT-021 recorded those oracles (`peptideMissingValueImputation.skip`, `limpaDpcQuant`, `limpaProteinNormalisation`, `ruvCovariates`), carried in every state's `args` and `config`; the components are identical. The fixture oracles are re-recorded on the final build for the verdict.

## Evidence

- Final code commit ff755910; install `kit_048` = `core/browser-smoke.2GA8g6` (over `kit_044`'s library); bindings at 6092b7eb, the re-recorded oracles bound at 10acea56.
- Pre-ticket full suite on a clean worktree of d6d4a385 ($P/kit/scratch-048-pre): 31,092 pass, 0 fail, 2 known errors (the two diagnostic files recorded since KIT-017), 140 skips, 346 files.
- Targeted tests on the final commit (clean worktree, the installed library): `test-workflow-checkpoint-dia-state.R` 148 pass / 4 skips (Q-22), `test-workflow-prot-declarations.R` 303 / 2, `test-workflow-partitions.R` 260 / 2; `test-dia-baseline-compatibility.R` after the two quarantines 26 / 2.
- `R CMD check --no-tests --no-manual` at ff755910 (git-archive export, `_R_CHECK_FORCE_SUGGESTS_=false`, the installed library on `R_LIBS`): 3 WARNINGs, 3 NOTEs, the baseline.
- The four screen journeys through the GUI on `kit_048` (recorded, not counted; `systemd` user slices, the driver, the verifier):

| Member | Stage archive | States | Verification | Elapsed | Peak (anon) |
|---|---|---|---|---|---|
| DIA lifecycle (`kit-lifecycle-complete`, screen corpus) | `kit-lifecycle-complete.kit048-dia-135034` | 18 complete | 0 mismatches against the KIT-005 oracles | 326 s | 2.42 GiB (2.06 anon) |
| DIA_LIMPA (`kit-journey-048m`, screen) | `kit-journey-048m.kit048r-limpa-141036` | 13 complete | accepted, 0 mismatches (`limpa48_screen`) | 234 s | 2.22 GiB |
| LFQ (`kit-journey-048l`, MaxQuant fixture) | `kit-journey-048l.kit048r-lfq-141346` | 10 complete | accepted, 0 mismatches (`lfq48_fixture`) | 125 s | 1.83 GiB |
| TMT (`kit-journey-048t`, PD fixture) | `kit-journey-048t.kit048r-tmt-141656` | 10 complete | accepted, 0 mismatches (`tmt48_fixture`) | 123 s | 1.81 GiB |

  Summaries: [lifecycle-complete.json](lifecycle-complete.json), [kit-journey-048m-stage.json](kit-journey-048m-stage.json), [kit-journey-048l-stage.json](kit-journey-048l-stage.json), [kit-journey-048t-stage.json](kit-journey-048t-stage.json). The first runs of the three screen members on this install (archives `kit-journey-048m.kit048-limpa-135513`, `kit-journey-048l.kit048-lfq`, `kit-journey-048t` refused) completed every state and mismatched the older records on the config summaries only: KIT-057's screen record of d4f0fbc8 predates KIT-056's `ruv_mode` (three summaries), KIT-021's fixture records of 88df9314 predate the limpa, RUV sweep and imputation defaults of KIT-046/047/041 (every summary, plus the `normalization_science` product whose identity carries the build). Every component was identical (compared field by field against the records made on this build). The oracles were re-recorded on `kit_048` by the KIT-053 and KIT-021 recorders into this directory and bound under their own prefixes; the three members were rerun against them.
- Post-ticket full suite, first pass on a clean worktree of ff755910 (the installed library): 29,067 pass, 3 fail, 13 errors, 295 skips, 346 files. Beyond the two known diagnostic errors, four files failed, all green before the ticket: `test-dia-baseline-compatibility.R` (2), `test-perf107-checkpoints.R` (1), `test-workflow-gui-memory-repairs.R` (1) and `test-workflow-gui-memory-lifetimes.R` (3 failures, 7 errors in 8 tests). Each failing test drove a module step on an in-memory manager or on an owner that was never activated (the adapter's legacy mode: `usesDIAReferenceStates` false, legacy `saveState` checkpoints), i.e. the retired branch (`state requires explicit legacy conversion`, `invalid saved-entity "input"`, `table read scope requires a Parquet descriptor`, the IQ step's export request without a source directory); none exercised the owner path. The twelve are quarantined under Q-22 with the file-level reasons in the registry; the four files rerun on the quarantined tree: 26/2, 38/1, 0/6, 80/11 (pass/skip). The second full pass on the quarantined tree is recorded below.
- Post-ticket full suite, second pass on the quarantined tree: (filled on completion).
- The counted member (`kit:journey-048:candidate:1`, the 300-sample DIA-NN corpus, ceiling 14,400 s, the KIT-005 oracles): declared below, launched only on the user's go.

Unrelated discovery, left untouched: the FASTA processing helper writes `aa_seq_tbl.tsv` into the process's working directory (`R/func_prot_annotation_fasta_processing_helpers.R`, the `vroom_write` beside the `aa_seq_tbl.RDS` metadata file); an in-process recorder run from the repository root left the file in the tree and the harness refused the next launch as dirty. Removed by hand; not a KIT-048 change.
