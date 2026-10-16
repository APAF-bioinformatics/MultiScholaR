# KIT-037: DA managed step acknowledges on what it shows

Status 2026-10-16: **complete.** Implemented at `9a886fbc`; observed live on KIT-036's screen pair (the candidate's `da_publication` record `status complete`, `browser_acknowledged TRUE`, `measurement_valid TRUE`); full suite at the commit green.

## The defect (KIT-017 post-hoc audit PH-04)

After the worker publishes and the controller completes, the managed step waits up to 30 s for the browser to acknowledge the accepted content: the status text and the plot named in the step's `ui`. The DA step named the static volcano (`volcano_plot_static`), which renders only when its sub-tab is shown, so in a journey and for any user who does not open that tab within 30 s the acknowledgement never arrived; `fail()` cancelled the action, logged an error and showed a 15-second error notification after a successful publication, while the status text read "[OK] Analysis Complete". Both counted KIT-017 attempts and the screen attempt 4 recorded `browser_ack_timeout`, `measurement_valid FALSE`.

## The repair (one mechanism)

- A managed step may show no plot at completion: `ui$plot = NULL`. The ready message (now `workflowManagedReadyPayload`) then carries no plot id; the browser's acknowledgement script treats the absent plot as ready and acknowledges on the text alone after the same 200 ms stable interval; the server-side handler is unchanged.
- The DA module binds `plot = NULL`; its `plot_identity` hook still sets the volcano's identity as the plot's alt text, an observation, never the acknowledgement.
- The peptide and protein QC steps are unchanged: their acknowledged content is the plot grid they show.

## Evidence

| Check | Result |
|---|---|
| JavaScript contract (`tests/testthat/fixtures/perf107-ack-contract.cjs`, node): the existing lifecycle plus the plotless case (text must match, the stable interval still applies, `plot_ready` sent) | `ack lifecycle contract passed` |
| `test-workflow-operation-registry.R` (new: the payload without a plot; the DA handler binds no plot) | 134 pass |
| `test-perf107-audit.R`, `test-perf-107.R` (the acknowledgement state machine, the head script) | 193, 21 pass |
| `test-module-ci-prot-da.R`, `test-prot-07b-da-handlers-compat.R`, `test-workflow-saved-da.R` | 245, 55, 109 pass |
| Full suite at `9a886fbc` from a clean worktree ([logs](check-logs/)) | 343 files, 29,645 pass, 0 fail, 140 skip; the one error is the PERF-101 snapshot test, the known runner artefact (passes alone); `R CMD check` 3 warnings, 3 notes, nothing new |
| Live: KIT-036 screen pair (2026-10-16 10:02–10:25, candidate `kit_036` which carries this ticket), the `da_publication` completion record | `status complete`, `browser_acknowledged TRUE`, `measurement_valid TRUE`; no error notification; the control (`kit_018`, before this ticket) still times out |

## Claim limits

The acknowledgement now attests the status text the step shows; the volcano's rendering is not attested by the step (it was never reliably attested before either, since the step timed out). The 30 s window and the 200 ms stable interval are unchanged.
