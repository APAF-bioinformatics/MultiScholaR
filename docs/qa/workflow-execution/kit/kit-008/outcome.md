# KIT-008: production worker dispatch for P14 and the sample filter

Status 2026-10-12: **complete.** Stage `kit-gui-008` decided `gui_member_observed` after two verifier-only corrections (member measurement untouched). Placement only: no memory decision, no promotion; the control stays `kit_011`. Authority: `tickets/kit/00-KIT-PLAN.md` (run ledger KIT-008: one `gui`-mode managed member on the screen corpus, 15 min after the user's amendment of 2026-10-12); protocol v1.3. Control: `kit_011` (`core/browser-smoke.AFaqC1`). Predecessor cases: PERF-107 (Apply, next consumer, Revert routing), PERF-113; audit findings B18, B19, B20, B21.

## What was built

1. **Dispatch rule (B18).** `workflowSessionApplication(session)` returns the managed application a session carries (set by `workflowManagedServer` under `launchWorkflowApplication`: Linux, systemd) or NULL. A module runs a registered operation on the managed worker when the application is present, and the same declaration in process otherwise (plain `run_app()`, macOS, Windows). No qualification flag takes part: the registry entry now carries the production variant an unqualified launch runs (`default_variant`: `p14` for the rollup, `physical` for the sample filter), and `workflowQualificationVariant` returns it when no qualification is bound.
2. **One reference module pattern.** `installWorkflowManagedStep(session, input, workflow_data, operation, publication, ui, deliver)` in `R/utils_workflow_execution.R` is the managed-step state machine both modules use: it reserves the application's one heavy action, starts the operation on the worker, polls it, delivers the accepted product to the module's outputs, requires the browser's acknowledgement of the delivered content before closing the action, retries a pending publication without dispatching science again, and releases the action on Revert, project replacement and session end. `readManagedPeptideResult` reads the accepted product for both operations.
3. **P14 in production (B21, hand-off item 8).** The rollup module dispatches `precursor_rollup` through the step under a managed application, so the GUI rollup gets the admission cap and recovery. In process, `runDIAReferencePrecursorRollup` now runs the declaration's backend (P14, the exact ordered-sum aggregate the worker runs) instead of the legacy row-indexed path. Backend switch checked before the member: [in-process-rollup-check.R](in-process-rollup-check.R), [in-process-rollup-check.json](in-process-rollup-check.json): on `screen_p14` (960,000 rows) legacy and P14 produce the same 48 output pages (480,000 rows), equal to the bound P14 screen oracle; legacy 16.5 s, P14 6.4 s in process (rerun at `9fe8028d`).
4. **Sample filter in production.** The sample module dispatches `minimum_peptides_per_sample` through the step without the `perf107-phase-a` gate; in process it runs `runDIAReferenceSampleFilter` (the declaration's physical mode, KIT-006's adaptive complete-sample batches, then the shared finishing step), no longer the native per-group loop. `workflowPeptideSampleQualification`, `runWorkflowNativeSampleControl`, the frozen-member ownership check and the mid-journey staging branch are gone (B19: after KIT-011 the owner is on disk from import).
5. **Qualification surface retired.** `workflowGuiLaunchRecord` accepts the KIT journey record only; `launchWorkflowApplication`'s gui-member qualification validation, the PERF-107 branch of `workflowQualificationVariant`, `workflowGuiQualifiedOperation`, the entry fixture (`runWorkflowGUIEntryFixture`, `workflowEntryHook`) and `finishWorkflowGUIAction`'s PERF-107 browser-membership requirement are removed. The application keeps its generic `entry_hook` slot (KIT-013) with nothing installed. `runPerf107GuiLifecycle` in the CORE-pinned driver is untouched, retired harness history (its scenarios were retired at KIT-013).
6. **Harness.** Stage `kit-gui-0NN`: one journey member on the ticket's bound build and corpus, launched with no record (the build's production default under a managed application), ceiling 600 s. Decision `kitGuiMemberDecision`: every step complete within the ceiling, populations complete and verified, every state, support entry and summary reproducing the reference screen member (`kit-lifecycle-complete`, the KIT-011 screen member on the control build; the per-build identity products reported, not judged), and a completed worker request for both operations in the member's project store. `bindings$tickets$008 = {kind gui, build kit_008, corpus screen}`.
7. **Tests.** `test-perf107-entry.R` retired (fixture tests; its hook-dispatch test moved to `test-workflow-execution-launch-record.R`); `test-perf107-audit.R` re-pointed (pending-publication retry on the dispatch rule; the native control, lifecycle classification and scientific-verification tests of the retired scenarios removed; a rollup dispatch test added: worker under a managed application, in-process declaration otherwise); `test-perf107-consumers.R` gains the in-process sample filter exactness check (physical declaration equals the native method, state committed); `test-perf107-checkpoints.R` and `test-perf107-synthetic.R` unchanged (checkpoints and the synthetic corpus are live); `test-perf-107.R` re-pointed to the dispatch rule and the shared acknowledgement machinery; registry and launch-record tests updated; harness tests cover the `kit-gui` stage and decision.

Grep evidence (acceptance 3): no `perf107-phase-a`, `workflowPeptideSampleQualification`, `qualification$entry` or `entry_fixture` reference remains in `R/`; `qualification` survives only for the adapters vehicle (`perf106-repair-v2` timed members) and in documentation.

## Results

| Check | Evidence | Result |
|---|---|---|
| Legacy vs P14 page identity on `screen_p14` (before the member) | `in-process-rollup-check.json` | identical 48 pages, 480,000 rows; both equal the P14 screen oracle |
| In-process sample filter through the declaration equals the native method | `test-perf107-consumers.R` | passes |
| Rollup: worker under a managed application, in-process declaration otherwise; pending publication retried without a second dispatch | `test-perf107-audit.R` | passes |
| Harness stage and decision | `test-workflow-kit-harness.R` (224) | passes |

Test files run after the change (working tree, proof library on `R_LIBS`): registry 35, launch record 25, perf-107 21, perf107-audit 63, perf107-consumers 29, perf107-checkpoints 49, perf107-synthetic 12, kit harness 224, operations 66, execution 64 (10 skipped installed-only), adapters (6 skipped installed-only), peptide sample shared 48, orchestrator shared 15, peptide module contracts 772, saved QC rollup 59, perf-106 193, memory mechanism member (13 skipped installed-only); all green.

Pre-execution audit 2026-10-12: [KIT-008-PREEXEC-AUDIT-20261012.md](KIT-008-PREEXEC-AUDIT-20261012.md), remediation [KIT-008-PREEXEC-REMEDIATION-20261012.md](KIT-008-PREEXEC-REMEDIATION-20261012.md) (the managed completion now carries the bus's QC parameters and prior progress so the worker's progress product equals the in-process one; the verifier takes the self-consistency branch for a screen member without a record; managed non-DIA sessions run in process; ceiling 900 s by user amendment).

## Member (counted, pending the user's go)

Stage `kit-gui-008` (launched 2026-10-12 11:56 after the pre-execution audit repairs, user "go"): one journey on the screen corpus, build `kit_008` = `core/browser-smoke.otE3Ck` (source `9fe8028d`), no launch record, managed application, ceiling 900 s.

| Measure | Value |
|---|---|
| Decision | `gui_member_observed` (after corrections 1 and 2) |
| Elapsed | 303 s |
| Whole-slice peak | 2.462 GiB |
| Anonymous peak | 2.119 GiB |
| Resident end | 1.745 GiB |
| Worker dispatches | `precursor_rollup` 1 request, complete; `minimum_peptides_per_sample` 1 request, complete |
| Journey | 18 steps complete, 18 populations, in-run fingerprints equal to the immutable store, 527 publication files |
| Reproduction of the KIT-011 screen member (`kit-lifecycle-complete`) | 0 mismatches over every state's tables, support entries, summary and products (the operation identity products reported, not judged) |

Against the KIT-011 screen member (in process, 486 s, 2.53 GiB slice, 2.10 GiB anonymous) the member ran faster with the same states; time is reported, memory recorded, nothing judged beyond completion, dispatch and exactness.

Corrections (verifier and decision only; the member was never rerun; the launch decision and the first recomputation are kept in `stage.rds$decisions_before_correction`; script [member-verifier-correction.R](member-verifier-correction.R)):

1. The verifier never ran at launch: `kitSelectorStage` mapped `kit:gui-008:candidate:1` to case id `candidate-1` while the stage's one case is `member`, so `verifyKitMember` aborted before writing a record (`member/verification-launch.log`), and the decision read the reference member under a bindings field the real bindings do not carry. Both fixed; the verifier ran from the corrected harness against the frozen bindings; the recomputed decision was `rejected` with 13 summary mismatches from the sample-filter state onward.
2. Those mismatches were bookkeeping, not science (user decision "1"): every table, support entry and product was equal; the sample-filter state's args recorded `core_utilisation = NA` (declaration run) where the native loop recorded its eight-worker multidplyr spec and copied the support references into the args, and the two word the step description differently; every later state's summary carried that args section forward. The canonical form now drops `core_utilisation`, the args-embedded `support_table`/`removal_ledger` copies and the step description (KIT-011 correction 1 did the same for the `declaration` field), and the gui decision recomputes both members' state fingerprints from their immutable stores under the current form. Recomputed: 0 mismatches.

## Claim limits

- **Hand-off to KIT-012.** The KIT-005 oracle fingerprints bound in `bindings.json` were computed under the earlier canonical form (native-loop args, descriptions). Before KIT-012's journey pair, the generator must recompute them from control-1's immutable store under the corrected form (`kitFingerprintStore`), or the pair's candidate (which now runs the two steps as declarations) would mismatch on bookkeeping as this member did. Recorded in the ledger.

- The in-process sample filter changed mechanism (native per-group loop to the physical declaration); exactness rests on KIT-006's accepted stage (full_primary and screen) and the new unit check, not on a new counted run.
- The managed path's browser acknowledgement keeps the PERF-107 message name (`multischolar-perf107-ready`) and JavaScript globals as the wire contract; only the server-side qualification logic was retired.
- Memory refusals mention free host memory through the existing admission messages (CORE-010); no GUI message text was changed in this ticket.
