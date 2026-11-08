# KIT-106 correction pass: the trx e2e lane (run 8), 2026-11-09

## The failure

Run 8 (`kit_106`, 9152e613): `test-e2e-transcriptomics-screen.R` pass 1, error 1, 1,378 s. The error was "app_wait_for_idle:
An error occurred while waiting for Shiny to be stable". It was raised by the `wait_for_idle` inside
`e2e_set_input_and_idle` straight after D3's tick, before the `wait_for_js` for the iframe. The browser log stops after
`shiny:value ...-de-volcano_glimma` and one `test_state_digest`.

## Diagnosis (the trx lane alone, under `systemd-inhibit --what=sleep`, at fcbb06b4 (= 9152e613's code), a scratch copy of the lane)

The scratch copy ([probe/diag-lane.R.txt](probe/diag-lane.R.txt)) runs the lane up to DA. It then ticks the box twice:

- the first tick with the value captured and its rendering prevented (`shiny:value` default prevented);
- the second tick with the value rendered.

- **What `volcano_glimma` renders:** an `iframe` with a `srcdoc` (the Glimma `glimmaXY` page). There is no error, no
  "Could not generate Glimma plot" and no "Preparing" placeholder ([probe/phase1-k.json](probe/phase1-k.json)). The value
  arrives once, 0.5 s after the tick.
- **In process:** the lane's app is `run_app(test_mode = TRUE)`, which has no managed application
  (`workflowSessionApplication()` is NULL). The page is built in the app process by
  `workflowOmicsInProcessVolcanoPageFactory`. A single value with no placeholder agrees with that.
- **Size:** the page is 958,462 bytes. The iframe's HTML is 1,130,297 bytes and the JSON value 1,134,191 characters. The
  trx_screen lane has 300 genes. Size is not the cause: with rendering prevented, the browser stays responsive and `get_js`
  answers at once.
- **Rendered:** after the second tick, no browser console event is logged after the value. Every `get_js` probe returns no
  value for 280 s ([probe/phase2-polls.txt](probe/phase2-polls.txt)). `wait_for_idle` fails
  ([probe/phase2-idle.txt](probe/phase2-idle.txt)). No `shiny:disconnected` is seen, and the app process stays up. The
  browser's page is blocked; the websocket is not closed.
- **The cause:** the page alone in a fresh headless browser:
  - With `alert` recorded instead of shown, the page calls `alert()` twice: "DataTables warning: table id=DataTables_Table_0 -
    Requested unknown parameter 'GeneID' for row 0, column 1" ([probe/page-probe2-record.json](probe/page-probe2-record.json)).
  - As shipped, `Runtime.evaluate` times out on every call after the load
    ([probe/page-probe2-plain.json](probe/page-probe2-plain.json)). The modal dialog blocks the page.
  - The page's `cols` are `gene, GeneID, GeneName, DisplaySymbol, logFC, negLog10PValue`. `trxGlimmaInputs` passes
    `display_columns = c("GeneID", "GeneName", "DisplaySymbol")`, but its `anno` moved `GeneID` into the row names, which
    Glimma shows as `gene`. So the table has no `GeneID` column, and DataTables (errMode `alert`) reports it.
- **The same page in a sandboxed iframe (no `allow-modals`):** it renders (one vega SVG, 10 table rows) and the page answers
  every call ([probe/wrapper-sandbox.json](probe/wrapper-sandbox.json)). Unsandboxed, it blocks as above
  ([probe/wrapper-plain.json](probe/wrapper-plain.json)).

KIT-097's e2e hang (`kit-097/outcome.md`, defect 2: "`de-volcano_glimma` then arrived, and every later wait for Shiny to be
stable failed") has the same signature on the same 300-gene lane. That widget sat in the page itself. KIT-097 attributed the
hang to the widget's size; the alert is the likelier cause (not re-run).

## The fix (b9820d9e; fa239458 before the rebase onto ec6f40da)

`R/mod_omic_da.R`: `omicDaVolcanoPageFactory` passes the managed or in-process factory's result through
`omicDaVolcanoPageFrame`. That function adds `sandbox="allow-scripts allow-same-origin allow-downloads"` to an iframe and
returns anything else unchanged. The page's scripts, origin, storage and downloads (Glimma's Save Plot / Save Data) work as
before. Only modal dialogs, popups, forms and top navigation are refused, so an `alert()` is ignored.

The change covers lipid, metab and trx, managed and in process. Proteomics' DA-time page does not pass through this
factory and is unchanged. The page's bytes, the asset, the dispatches and every DA result, threshold, default and recorded
state are unchanged.

**Failing first:** a new block in `test-omic-da-module.R`, "(KIT-106 correction, F106-06)", covers trx in process and
lipid managed. Without the R change it fails 2 and errors 1 ([test-first.txt](test-first.txt)). With the change it passes
(in the module file's run below).

## Reruns on the fixed build (scratch library over m4DXz7, the app from the worktree source)

See [reruns.txt](reruns.txt).

## Findings (proposed, for the lead's findings.md)

- **F106-06 (resolution).** The cause is a DataTables `alert()` in trx's page, a modal that blocks the browser page. Size,
  the websocket and the test's waits are not the cause. Fixed by the sandboxed iframe (b9820d9e) for every omic's
  on-demand page.
- **F106-07 (new, not fixed here).** The trigger itself stays:
  - trx's interactive volcano names a `GeneID` table column that Glimma's table does not have.
  - `trxGlimmaInputs` (`func_trx_publication.R:229`) keeps the notebook's `display_columns`. That value is pinned by
    `test-trx-publication.R:152` against the KIT-040 notebook.
  - In the sandboxed frame the table's `GeneID` column is blank, and the ids show under `gene`. Before the fix, a user saw
    two DataTables warning dialogs.
  - Changing the columns is a display change against a notebook-parity pin, outside this pass's files. The lead decides.
