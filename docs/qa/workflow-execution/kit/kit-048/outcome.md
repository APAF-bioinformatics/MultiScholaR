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

- The three entry points `diaStateHistory`, `saveDIAOwnerState` and `readDIAStateObject` **forward** to the manager's commit surface rather than asserting a `DIAStateAdapter`: in the application the manager is always the owner (from session start), but the KIT-042 legacy ingest stages a legacy RDS session on a `WorkflowState` of its own before adopting it into an owner (the ticket allows that construction), and the module contract tests pass fakes exposing the same three methods. Asserting the adapter broke the legacy ingest (`buildProtDesignStateCheckpoint` on the staging manager); forwarding keeps one named entry point per operation and the grep criteria, at the cost that the entry points do not reject a non-owner. Recorded as the chosen bounded form; the ticket's item 4 (an owner fixture for the twenty mock-based module tests) is deferred: the tests that only drove the in-memory apply/revert path are quarantined (Q-22, 53 tests in twelve files; Q-21 for the removed fallback S4 tab's tests), the module tests whose logic is independent of the manager's class run unchanged on their fakes.
- The DA session validation accepts a saved reference or the `ProteinQuantitativeData` S4 a legacy RDS session carries (the KIT-042 ingest), as before; my first cut refused the S4 and the legacy adoption tests caught it.
- The server tests stub the session-owner creation point (`createDIASessionOwner`) in place of `WorkflowState$new`.

## Evidence

(filled as the checks run)
