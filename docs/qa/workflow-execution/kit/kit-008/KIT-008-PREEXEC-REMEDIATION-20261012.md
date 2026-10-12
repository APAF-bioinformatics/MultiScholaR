# KIT-008 pre-execution remediation (2026-10-12)

Closes findings of [KIT-008-PREEXEC-AUDIT-20261012.md](KIT-008-PREEXEC-AUDIT-20261012.md) (user: "ill take all your recommendations - fix findings"). The audit register is unchanged. Code at `9fe8028d`; install `kit_008` = `core/browser-smoke.otE3Ck` (replaces `Hag4dl`, which was never run).

| ID | Action | Closure evidence | State |
|---|---|---|---|
| PE-01 | `kitSelectorStage` accepts `kit:gui-0NN:candidate:1`; the verifier's KIT-005 oracle branch applies only to full-corpus members without a record; a screen member without a record takes the self-consistency branch and is reproduced against its reference screen member by `kitGuiMemberDecision`. | `test-workflow-kit-harness.R` 225 pass (selector test added); stage `kit-gui-008` validates on `otE3Ck`. | closed |
| PE-02 | The GUI step hands the action a bounded completion context (`workflowCompletionContext`: the bus's `peptide_qc` parameters and prior progress slots, ≤ 4 MiB, no live owners); `runWorkflowAction` keeps it on the job; `completeWorkflowAction` finishes the worker result on `newWorkflowCompletionWorkflow(manager, context)`, so the `qc_progress` product equals the in-process one; `readManagedPeptideResult` merges the step's `peptide_qc` entry into the bus instead of replacing the bus. | `test-perf107-audit.R` 78 pass: context built from the bus, validated, seeded; the reader keeps the earlier entries; the rollup dispatch carries the context. | closed |
| PE-03 | User amendment: the member's ceiling is 900 s (ledger row and KIT-008 ticket rev 1.2; bindings `tickets$008$ceiling_seconds`). | bindings, ledger | closed |
| PE-04 | Both modules dispatch to the worker only when the session carries a managed application **and** the owner is reference-backed; otherwise the declaration runs in process. | `test-perf107-audit.R`: a managed launch with an in-memory owner runs in process. | closed |
| PE-08 | `in-process-rollup-check.R` rerun at `9fe8028d`: identical 48 pages, equal to the oracle (legacy 16.5 s, P14 6.4 s). | `in-process-rollup-check.json` | closed |
| PE-05 to PE-07 | No change. | — | accepted |
