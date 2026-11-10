# KIT-113 W3: the frame dispatches through the runner (D3), Phase C (2026-11-10)

Bounded item: W3 of `tickets/kit/KIT-113.md`. Files: `R/mod_omic_enrich.R` (`runOmicEnrichAnalysis` only),
`tests/testthat/test-omic-enrich-module.R` (15 blocks appended, 8 from one loop; nothing existing edited), this folder.
Exit: the W3 blocks fail first on the start commit and pass after; the unmanaged frame identical (every W1 (b), F107 and
KIT-113 W1 frozen capture unchanged); the targeted set and the five proteomics lanes at their start counts.

Start commit 1b5d7f44 (performance-v2, W2 merged). Code and tests: `6eacd081`. R 4.6.1, kit_112's library, load_all,
strict skips, `TMPDIR=kit/tmp-113`, one process per file through `tools/testing/test-baseline.R file`, at most 4 at once.

## What changed (`runOmicEnrichAnalysis`, 692 -> 759 lines in the file, +67, insertions only)

- The `req()`s and the working notice are unchanged. Then the runner lookup:
  `workflowOmicsManagedActions(workflowActionSpec(config$operation()$name)$omic_type)` ("proteomics"). The frame installs
  no runner.
- **No runner** (unmanaged, or the frame mounted alone): the start commit's code follows, byte-identical (the diff
  against 1b5d7f44 deletes no line).
- **A runner:**
  - F112-02's guard first: `runner$step$active()` or `runner$step$pending()` (isolated) refuses the click with the
    step's own warning (`paste("A", omic, "step action is already running")`, type warning, default duration), removes
    `enrichment_working` and returns: nothing resolved, written or dispatched.
  - The progress is the runner's (`runner$progress("Running enrichment analysis...")`), passed to both halves.
  - The request half (`requestWorkflowEnrichment`, the same parameters list the in-process call builds) inside
    `tryCatch`; an error takes today's error notice (text, type error, duration 10, today's console lines), then the
    tab holds the source reference and `enrichment_working` is removed. R7's check is the request half's (W2).
  - `runner$run` with one action: `operation = "enrichment"`, `parameters = request$parameters`,
    `input = manifest$source`, `deliver` keeps the delivered value, `context = list(state_name = "enrichment_completed",
    config = list(), description = "enrichment")` (non-NULL, as `protNormManagedAction` builds a state-producing one).
  - `done` (isolated: the runner calls it from its step's poll observer): the accept half
    (`acceptWorkflowEnrichment(..., request$context, delivered, progress)`), `holdOmicEnrichRun`, then the progress
    closed, `tab$current_s4_object`, the completion notice and `enrichment_working` removed. An accept-half error takes
    today's error notice.
  - `failed` (isolated): today's error notice, the progress closed, `tab$current_s4_object`, `enrichment_working`
    removed.

## Tests (test-omic-enrich-module.R, appended; prefix `omicEnrichK113`)

KIT-107 W1 (b)'s fixture and driver unchanged (`omicEnrichW1CaseSpec`, `omicEnrichW1Drive`). After W1 (b)'s mount the
session's owner exists, and the wrapper server then installs, scoped to the case: `omicsManagedChainMocks` on that owner;
`runWorkflowAction` replaced by the real dispatch's gates (`assertDIAExportReady`, `assertWorkflowOperationDispatchable`,
`assertWorkflowActionInputSelected`, `validateWorkflowActionInput`) and a named-input dispatch through
`workflowWorkerMain`, the controller's acceptance and omics completion (KIT-111 W2's `protNormW2Dispatch`, reproduced);
a poll that can hold the action `running` or end it `publication_pending`; a counted request half; a compute half that
can be forced to fail; then the proteomics runner on the `normalization` scope (as the shell installs it), its progress
display recorded in place of shiny's (as `omicsTabServer` replaces it). Each step records the runner's progress, the
dispatches' outcomes, the reserved operations and the request halves (`k113`).

Design decision: the DA step runs with the session's application withheld (in process), so W1 (b)'s DA captures stay
frozen; every later step is managed. The managed DA is KIT-112 W3's, proven there.

1. KIT-107's eight cases managed: every field of every step equals W1 (b)'s frozen value, and each session's KIT-113 W1
   D6 fields equal theirs, apart from D6's recorded exclusions, each asserted exactly as W2 asserts them (store
   transfer values; the request half's 0.5 before a worker refusal; `uniprot_dat_cln` before the missing-table refusal;
   `enrichment_organism_filter` unwritten after a worker failure and after the taxon refusal, F113-04), plus D3's managed
   notices (the runner's working notice and removal before the completion notice; for a worker refusal also the
   runner's "Error applying proteomics step: <the same message>" before today's error) and the runner's progress
   (today's sequence and "close"; shiny's progress empty). Dispatches: input = head on a first Run, request < 65,536 B,
   value <= 8 MiB, the completion's provenance pending before the accept half; acceptances per case as W2.
2. A second Run: dispatched on the manifest's source while the head is `enrichment_completed`; the first Run's values.
3. A worker failure: today's error and the runner's, `enrichment_working` removed, no state, product, resource, provenance
   or pathway file, the annotation table on the bus (D6); then a Run gives W1 (b)'s run step and the D6 fields.
4. A missing DA table: today's refusal text (and the runner's), no state, product or pathway file, `uniprot_dat_cln`
   written (D6's divergence).
5. R7: the tab's manifest made to differ from the head's: today's error notice with the request half's message, nothing
   dispatched, committed or saved to the store (`uniprot_dat_cln` written first, see draft F113-06).
6. A click while the runner's step is active (a held Run): the step's warning, `enrichment_working` removed, no request
   half, nothing moved; the held Run then completes with W1 (b)'s run values. A click while the step holds a pending
   publication: the same refusal.
7. The frame mounted alone on a managed session (no runner): W1 (b)'s run step exactly, nothing dispatched.

The unmanaged frame is covered by the 14 existing blocks, unchanged and passing (997/997, as at W2).

## Divergences recorded (conventions, existing-code scope rule)

- The Run stays a managed step, not a `shiny::ExtendedTask` (the frame's header names KIT-113 for this; the managed step is
  the series' mechanism).
- The managed branch keeps today's `cat()` console lines and `e$message` notices, as the in-process branch does.
- The managed branch is inline in `runOmicEnrichAnalysis` (the work item allows that one function only), so the function
  is longer than the conventions' guidance; the in-process code is untouched below it.

## Draft findings for the lead (numbers are the lead's)

- **F113-06 (R7's refusal point, information or a D6 row).** The request half checks the held manifest against the head
  after the annotation resolve and the mapping (W2's order), so a refused Run has already written `uniprot_dat_cln` to
  the bus (when it had none) and shown the no-mapping warning (when the filter is on without a mapping), as D6 records
  for the missing-table refusal. Nothing is saved to the store or dispatched. In process there is no such refusal at
  all (the stale manifest is run and committed). Proposed: one more of D6's recorded exclusions, the same kind.
- **F113-07 (the pending-publication guard, information for P-42 or the parked list).** The guard refuses an enrichment
  click while the runner holds any pending publication, including another proteomics tab's, and only that tab's own
  button retries it; the notice is the step's text ("already running"), not a "press that step's button again" message.
  The enrichment declaration stages nothing, so an enrichment action itself should not leave one. D3 asks for exactly
  this; the wording is today's step text. No change made.
- **Information, the worker-failure text.** In these tests the chain raises the worker's own message, so the managed
  notices carry the in-process text exactly. On a real managed session the controller sees a terminal poll status, so
  the text would be the step's ("managed action ended: worker_failed"); the notice's format, type and duration are
  today's. This is D3's "notices and progress" exclusion.
- **Information, found and fixed inside W3 before the commit.** The first draft of `done` passed the accept half as a
  lazy argument to `finish`, which closed the progress first; the runner-progress expectation caught it, and
  `force(completed)` fixed it (test-first.txt).
