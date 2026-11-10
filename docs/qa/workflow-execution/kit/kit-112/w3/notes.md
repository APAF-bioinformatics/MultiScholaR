# KIT-112 W3: the run plugin dispatches through the runner (D3), Phase C (2026-11-10)

Bounded item: W3 of `tickets/kit/KIT-112.md`. Files: `R/mod_prot_da_handlers_run_analysis.R` (the run observer and the
runner lookup before it), `tests/testthat/test-omic-da-module.R` (7 blocks appended), this folder. Exit: the W3 blocks
fail first on the start commit and pass after; the unmanaged plugin identical (prot_main, prot_load and every other
frozen capture unchanged); the targeted set and the five proteomics lanes at their counts.

Start commit 6f41a5dd (performance-v2, W2 merged). Commit: `661186bc`. R 4.6.1, kit_111's library, load_all, strict
skips, `TMPDIR=/home/doktersmol/.cache/k112w3`, one process per file through `tools/testing/test-baseline.R file`, at most
4 at once.

## What changed

- **Runner lookup (before the observer).** On a Shiny, proxy or mock session the plugin takes
  `session$userData$workflow_omics_actions[["proteomics"]]` (KIT-111's runner, installed by the normalisation shell, which
  starts first in tab order; one per session), else installs one with KIT-111's call and guard
  (`installWorkflowOmicsActions(session, workflow_data, "proteomics", experiment_label %||% "dia_baseline_smoke",
  experiment_paths)`). Unmanaged, `installWorkflowOmicsActions` returns NULL, so there is no runner.
- **The click, managed with a runner** (`workflowSessionApplication` non-NULL and a runner), after today's no-data
  refusal, `da_working` and settings:
  - F112-02 guard 1: if the `da_publication` step or the runner's step is active, "A DA publication action is already
    running" / "A proteomics step action is already running" (the step's own text, type warning), `da_working` removed,
    nothing else (no request half, no commit).
  - F112-02 guard 2: if the publication holds a pending retry, `da_working` removed and today's
    `managed_step$apply(list(annotations = ...))` (the retry), with no request half and no dispatch.
  - Otherwise `proteomicsDiaBinding$da$request(...)` (diaDaRequest) in the controller; errors take today's tryCatch
    path and text. Then `runner$run` with one `protNormManagedAction("DA-contrast-output", c(step, route =
    diaActiveRoute(manager)), request$source, NULL, deliver)` per step. `deliver` keeps `result$value` by contrast name
    and adopts `result$display$base_fit` into `da_base_fit_reference` only when not NULL (R1). `done`: inside
    `withDIAStateScope(manager, ...)` the accept half (`proteomicsDiaBinding$da$accept`), `da_working` removed, then the
    publication call (today's `managed_step$apply`, now through a local `publication()` shared with the retry and the
    legacy branch). `failed`: today's error report ("Error in DA analysis: <message>", error, 10 s, plus today's console
    line) and `da_working` removed.
- **Unmanaged:** today's `runDIAReferenceDaAnalysis(..., publish = TRUE)` call, notices and errors, unchanged in effect.
  The error handler is the same body moved into a local `report_error` (shared with `failed`).
- **Managed with no runner** (not reachable when the application is on the session at server start, as the shell and
  KIT-111 assume): today's managed branch is kept (kernel in the controller, then the publication).
- Lines: 77 -> 127 (+50, at the ticket's upper estimate).

## Tests (test-omic-da-module.R, appended; prefix `omicDaK112`)

A fixture as `omicDaW1ProtSetup` (MCI-009's object on an on-disk DIA owner, the same bus and tab) made before the
session so that `omicsManagedChainMocks` holds the owner; each case's fixture and mocks live for its drive only. Driven by
`omicDaW1Drive` (its notification doubles and captures). Around the chain: the dispatch also runs the real gates
`assertDIAExportReady(owner$manager)`, `assertWorkflowOperationDispatchable` and `assertWorkflowActionInputSelected`
before delegating (the F111-03 lesson, so far as the chain allows); the request halves and runner installs are counted;
the poll can hold an action `running` or end the publication `publication_pending`; the kernel can be forced to fail.
`da_publication` is acknowledged as the browser does (`workflow_ack` with the chain's reservation id and the step's
plot identity).

1. One contrast: operations `DA-contrast-output`, `da_publication`; notes `da_working`, the runner's working notice and
   removal, `remove da_working`, the publication's working notice; after the ack "DA publication applied
   successfully". Values, updates, tab_status, history, states, tab (volcano identity masked), the manifest's rows and
   the published files equal prot_main's frozen in-process run.
2. Two contrasts (B_vs_A, A_vs_B): `DA-contrast-output` x 2 then `da_publication`; each request carries the in-process
   parameter names plus `route = "iq"` on the head; the bus's, the tab's and the head's manifests are identical; the
   same shared fields equal the same fixture run unmanaged.
3. Worker failure then a click: "Error applying proteomics step: worker failed: ..." (runner, 15 s), then "Error in DA
   analysis: worker failed: ..." (10 s), `remove da_working`; no manifest, tab unchanged; the next click completes.
4. Click while the runner's action runs: refused before the request half (1 request half in total), states and history
   unchanged; the held run then completes.
5. Click (with a changed q-value) while the publication awaits its acknowledgement: refused, no commit; the
   publication then completes.
6. Click while the publication holds a pending retry: only `da_publication` is reserved again; no request half; it
   completes.
7. Mounted alone (managed): one runner installed. Full session: the pre-installed runner (on the `normalization`
   scope, as KIT-111's shell) is the plugin's, no second install, and its `run` is the one used. Unmanaged: no runner,
   the run completes in process.

## Notices on the managed path (as D3 records)

Today's managed set plus the runner's "Applying proteomics step ..." working notice per contrast; a worker failure also
shows the runner's "Error applying proteomics step: ...". The BH-fallback warning and the completion notice are still
not shown on the managed path (F112-01, queued to KIT-109; not added here). F112-05 (the `da_publication` completion's
pending provenance) is present and not addressed.

## Draft findings (for the lead)

- **Runner retry state not guarded (draft).** F112-02's guards cover the runner's step being active and the publication
  step's pending retry, as the ticket names them. If a `DA-contrast-output` completion itself ended
  `publication_pending` on the runner's step (its `omics` publication), the next click would run the request half (a
  possible arguments commit) and then `runner$run` would resume the earlier click's actions with the earlier callbacks
  (`runner$retrying()`, utils_workflow_omics_actions.R:63-70). Not reproduced; no lane or journey reaches it. A third
  guard (`runner$retrying()` -> `runner$run(list())` without the request half) would mirror guard 2. Not done (outside the
  ticket's two guards).
- The plugin reuses `protNormManagedAction` (mod_prot_norm_workflow_helpers.R) by name, as D3 says ("as KIT-111's
  protNormManagedAction builds it"); its name is normalisation-specific.
