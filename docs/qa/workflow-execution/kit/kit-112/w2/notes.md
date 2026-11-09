# KIT-112 W2: the split, the handler, the registry, the adaptor (D3, D4), and F112-03 (2026-11-10)

Bounded item: W2 of `tickets/kit/KIT-112.md` (Phase B) plus F112-03 as its own commit. Files: `R/utils_workflow_dia_da.R`,
`R/utils_workflow_operations.R` (the registry loops only), `R/utils_workflow_prot_declarations.R` (the adaptor's `da`
entry), `tests/testthat/test-workflow-saved-da.R` (appended blocks), `tests/testthat/test-workflow-operation-registry.R`
(D11's two retargets). Exit: W1 (a) identical, W1 (b) at core level equal, D11 failing first then passing, the
targeted set and the DIA set at their counts.

Start commit 525cbd78 (performance-v2). Commits: `8085d7d0` (W2), `9e2ee21d` (F112-03). R 4.6.1, kit_111's library,
load_all, strict skips, `TMPDIR=/home/doktersmol/.cache/k112w2`, one process per file through
`tools/testing/test-baseline.R file`, at most 4 at once. "Before" runs used a `git archive` copy of the start commit (or
of `8085d7d0`) with only the test file replaced.

## What changed

- `runDIAReferenceDaAnalysis` (name and signature kept) is the request half, the flags as today, the compute half per
  contrast with today's per-computed `da_dependencies` publish, the accept half, then today's publication and display.
  - `diaDaRequest`: today's :345-389 (gate, source, defaults, design-only validation, manifest skeleton, arguments
    commit) and the steps. The arguments commit is moved verbatim into `commitDIADaArguments` to keep both under ~50
    lines (a sixth private name; the lifecycle trace list of the Simulation Setup is unaffected).
  - `computeDIAReferenceDaContrast`: sets `da_reference_execution`/`da_source_reference`, restores the previous values on
    exit, runs the declaration; returns `list(run, base_fit)`, `base_fit` = the context's base-fit reference when the run
    was `computed` (the worker's owner is fresh per action), else NULL.
  - `acceptDIADaManifest`: today's :407-423.
  - `computeWorkflowDiaContrastAction`: route from the request (refused if absent or not iq/limpa), stripped; compute
    half on the named input; a hit rooted with the source and the value's native references; returns
    `list(checkpoint = source record, reference, native, display = list(base_fit))`.
- Registry: `DA-contrast-output` moved from the in-process loop to KIT-111's omics loop (handler, omics completion and
  publication, `named_input = TRUE`, `section = NULL` by the loop's rule). `enrichment` stays in process.
- Adaptor: `proteomicsDiaBinding$da` as D4.
- Request-half reads: exactly today's (state records via `readDIAStateMetadata`, the design via
  `readDIAExecutionComponents(fields = "design_matrix")`, the head's record for the arguments commit, the commit's
  reference validation). No read added. The protein matrix is first read in `runSavedOperation` inside the compute half.

## F112-03

Cause: `writeDIAInteractiveVolcanoOutputs` handed `generateProtDAVolcanoPlotGlimma`'s result to `saveSelfContainedWidget`;
when a page is saved for the contrast's result and threshold (KIT-016) that result is the saved page's iframe tag, and
`htmlwidgets::saveWidget` stops ("invalid type/length (symbol/0) in vector allocation").

Fix (utils_workflow_dia_da.R only):
1. `writeDIAInteractiveVolcanoOutputs` looks the page up first (the lookup `generateProtDAVolcanoPlotGlimma` does first
   anyway); a saved page is written as the file (`writeBin(charToRaw(page))`, the saved bytes), not rebuilt,
   re-serialised or saved again. Without a saved page the code path is today's.
2. `saveDIAInteractiveVolcanoAsset` keeps an asset already saved under the identity. Needed for the managed path:
   with fix 1 only, the worker reuses the page, and the controller's `publishWorkflowDaResult` (utils_workflow_execution.R,
   unchanged) saves every staged page as an asset, so a second asset with the same identity was saved and every later
   lookup of the page became "multiple" (ambiguous). Probe: `probe-f11203-worker.R`, output
   `probe-f11203-worker-on-writer-fix-only.txt` ("asset lookup for the reused page: multiple"). In process the guard is
   always a miss (the writer no longer calls it for a reused page).

Not chosen: `saveSelfContainedWidget` (shared by the enrichment and omics writers) would have had to unwrap a shiny tag
to recover the page; the reuse belongs where the page is resolved.

## Test first

See `test-first.txt`.

## W1 (b) at core level (the equivalence)

Blocks appended to `test-workflow-saved-da.R` (`savedDaW2*`): the managed click = `diaDaRequest` in the controller,
one `DA-contrast-output` request per contrast through `workflowWorkerMain` on the named input (with `route`), the
acceptance and completion, `display$base_fit` adopted only when not NULL, `acceptDIADaManifest` in the controller, then
`da_publication` through the worker and the controller's publication, `readManagedDaResult` and `completeDIADaDisplay`.

- Every D6 case passes: one contrast, two (base-fit hit), repeat, reopened rerun (R1: base-fit fields NULL, as frozen),
  non-head (R7), changed arguments, q only, refusal, worker failure then retry, DIA_LIMPA.
- The DIA block also runs a twin owner through today's managed path (runner with `publish = FALSE`, then the same
  publication step): worker path = today's managed path in every captured field except the `saveDIAStateValue` count.
- Worker path = W1 (a) frozen in-process values in every field except, by design:
  - the `saveDIAStateValue` count (head snapshots; D6);
  - `pending`: on both managed paths the `da_publication` completion leaves `da_publication` pending
    (`completeWorkflowAction` sets it; KIT-019, unchanged), where in process the accept half leaves none; after a
    worker failure the controller has none (the worker's owner held it; D6);
  - `provenance` (committed step provenance) in `dia_qvalue_only` and `dia_refused` only: the DA-arguments commit that
    follows a managed publication records that pending `da_publication` entry as its step provenance, on today's
    managed path too (draft finding below);
  - non-head: the publication is refused on both managed paths ("DA manifest does not belong to the head state"), so
    `outputs`, `files`, `publication`, `display`, `error`, the run value (`manifest`, `model_unset`) and the
    `publishDIAActionOutputs` count differ from in process; the head manifest, states, products, fits, entities,
    dependencies and record equal W1 (a)'s.
  - `products`/`dependencies`: no difference observed (the capture already drops `product_dependencies`).
- R5: 13 contrast values, each < 8 MiB, no live owner, names `checkpoint, reference, native, display` (+ the worker's
  `operation, roots, owned`); base fit delivered exactly on the contrasts that ran the kernel.
- Route stripping: a unit block shows the compute half receives exactly the step without `route`, for iq and limpa,
  and a request without a route is refused (`workflow_operation_parameters`).
- Adaptor: `proteomicsDiaBinding$da` names the registered omics operation with the handler, `named_input`, the
  `da_publication` step and existing functions; the runner's formals are unchanged.

## Draft findings (for the lead)

- **F112-05 (draft).** On the managed path (today's and the worker's alike) the `da_publication` completion leaves
  `da_publication` as the controller's pending step provenance, and the next state commit (a DA-arguments commit on a
  later DA click) records it as that state's step provenance; in process that commit records the provenance it would
  without it. Pre-existing since KIT-019/KIT-025, not changed by KIT-112 (the twin comparison shows both managed paths
  equal). Recorded, not fixed (utils_workflow_execution.R is outside W2).
- The KIT-019 comment above the `da_publication` registration ("the DA kernel (limma) stays in process",
  utils_workflow_operations.R:683-685) is now stale; outside W2's lines, left.
