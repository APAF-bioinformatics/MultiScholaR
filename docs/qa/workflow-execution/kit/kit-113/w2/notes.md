# KIT-113 W2 notes (the W2 agent, 2026-11-10)

Bounded item: D3's core. Split `runWorkflowEnrichment` and `runDIAReferenceEnrichment` into helpers called in today's
order (names and signatures kept), add the managed request, compute and accept halves and the worker handler, carry
the annotation match results by reference (F113-02), and move `enrichment` from the in-process registry loop to the
omics loop. Start commit `dcdc4a3f`. The frame (`R/mod_omic_enrich.R`) is W3's and is unchanged.

## Commits

- `2984636c` code and the W2 test blocks;
- `c2db9053` F113-02 (its test first, then the fix);
- the evidence commit (this file, `test-first.txt`, `targeted.txt`).

## What changed

- `R/utils_workflow_enrichment.R`: `runWorkflowEnrichment` is now `workflowEnrichmentContext`,
  `readWorkflowEnrichmentDa`, the annotation resolve (in place), `workflowEnrichmentFilter` (the match, then
  `workflowEnrichmentMapping` with its warning, `applyWorkflowEnrichmentFilter`, `persistWorkflowEnrichmentFilter`),
  progress 0.5 and the method (in place), `workflowEnrichmentRun` and `publishWorkflowEnrichment`. Each is a contiguous
  run of the start commit's statements, in the same order. The managed halves: `requestWorkflowEnrichment`,
  `computeWorkflowEnrichment`, `acceptWorkflowEnrichment`, with `workflowEnrichmentKernelName` /
  `workflowEnrichmentKernel` (by name, looked up when called) and the two transfer helpers
  (`saveWorkflowEnrichmentValue` / `readWorkflowEnrichmentValue`, codec `rds`, NULL stays NULL).
- `R/utils_workflow_dia_da.R`: `runDIAReferenceEnrichment` is `computeDIAReferenceEnrichment` then
  `recordDIAEnrichmentResources`; `computeWorkflowDiaEnrichmentAction` beside it (`commit = NULL`, `display` holds
  references and the filter's outcome).
- `R/utils_workflow_operations.R`: the in-process loop is gone; its comment moved to the omics loop, which gains
  `enrichment`. `diaEnrichmentOperation` is unchanged.

## Decisions taken inside D3

- `workflowEnrichmentRun` takes `kernel = NULL` and `run = config$run` beyond D3's listed signature: `kernel %||%
  config$kernel$select()` keeps the selection after the argument validation (today's order), and the worker passes
  `computeDIAReferenceEnrichment` as the run (D3: "with computeDIAReferenceEnrichment as the run").
- The records' query summary takes the contrast count from the declared `results` and the taxon from the declared
  `taxon_id` on the managed path; in process `runDIAReferenceEnrichment` passes `arguments$taxon_id` as today. They are
  the same values (the declared parameters are `arguments[...]` plus `results`).
- The request checks the held manifest `identical()` to the head snapshot's `da_results` and refuses otherwise
  (`workflow_enrichment_input`), as D3 and R7 say. The request half runs outside the owner's scope, as W1's emulation
  did; the accept half runs inside `withDIAStateScope`.
- `requestWorkflowEnrichment` returns the context with the method and the run's inputs, so the accept half's signature
  is D3's (`workflow, config, context, delivered, progress`).

## Test design (core level)

The new blocks replace only the frame's core call (`runWorkflowEnrichment`, mocked in the namespace) by the managed
halves: request half, one `enrichment` action through `workflowWorkerMain` on the request's named input (a local copy of
KIT-111 W2's `protNormW2Dispatch` on the in-process test worker `omicsWorkerController`), the controller's acceptance
and omics completion, `resumeWorkflowAction` / `readManagedOmicsResult`, then the accept half. KIT-107's eight W1 (b)
cases then run their frozen steps unchanged, so every captured field of every step is compared with the frozen
in-process value, plus the KIT-113 W1 block's fields after each session. The comparison allows only D6's recorded
exclusions, each asserted exactly:
- the store's `dia-state` entries: +2 per completed Run (the annotation table and the results; the match results are
  NULL on this fixture and travel as NULL, asserted per dispatch), +2 for `organism_filter/mapped` (the annotation table
  and the mapping, then the worker fails on today's `clean_acc` error), +1 for `q_cutoff` and `table_removed` (the
  annotation table), 0 for `run_taxon_abc` (the request half refuses before saving);
- progress gains the request half's 0.5 before a worker refusal (`mapped/run`, `table_removed/run_table_removed`);
- `uniprot_dat_cln` is written by the request half before the worker refuses a missing DA table (`table_removed`);
- `enrichment_organism_filter` is not written when the Run fails after the point where in process writes it
  (`q_cutoff/run_q_cutoff`, a worker failure, D6; `pre_da/run_taxon_abc`, a request-half failure: draft F113-04).
Everything else is identical, including every notice (the worker's error text is the in-process message), the states,
history, products, displays, tab outputs, returned values, files and the W1 block's record, product fingerprint,
roots, provenance and analysis record step. `prefill` has no Run; it runs unchanged and dispatches nothing.

Also: a second managed Run (input = the manifest's source while the head is `enrichment_completed`; the same values; the
record's previous state is `enrichment_completed`, as in process); a forced worker failure (nothing committed,
recorded or published, no provenance pending, today's error notice) then a Run that completes with W1 (b)'s run step
and the W1 block's values; the kernel's name both ways, looked up when called (a recorder in the double's binding maps,
a double captured before it is refused), and refusal of any other function or name; the sizes (R3) per dispatch
(request 16,824 to 22,480 B by object.size, value about 65 KB, both within their limits).

## Divergences recorded (conventions, existing-code scope rule)

- The moved statements keep their `*ProtEnrich*` callees, `cat()` logging and `@` slot reads.
- `workflowEnrichmentRun` has 11 parameters and `publishWorkflowEnrichment` 9 (conventions: 5 to 7): D3 names the
  first signature; both keep today's statements without a new parameter object.
- `utils_workflow_enrichment.R` grew by 235 lines against the ticket's "+70 to +140": most of it is the roxygen of 15
  new internal functions (@noRd); no logic is duplicated.

## Draft findings for the lead (numbers are the lead's)

- **F113-03 (F113-02's mechanism, corrected).** A value over 8 MiB is refused earlier than F113-02 says: the worker's
  own handoff (`writeWorkflowExecutionMessage`, `utils_workflow_execution.R:1300-1302`, the same 8 MiB bound) fails to
  write `worker-result.rds`, so the controller sees the worker fail (no result), not an `invalid_result` acceptance.
  The fix is the same and is made (`c2db9053`). Information only.
- **F113-04 (D6's list, one more case of a recorded exclusion).** D6 records that the bus's `enrichment_organism_filter`
  is not written after a worker failure. It is also not written when the request half fails after the point where in
  process writes it: the taxon refusal (`refusals/pre_da/run_taxon_abc`, "organism_taxid must be a single positive
  integer NCBI taxonomy ID.", raised by the method) comes after the filter metadata in process and before any
  metadata on the managed path (the metadata is the accept half's). Same nature as D6's row; the notice is today's.
  Proposed: add it to D6's recorded exclusions.
- **F113-05 (lines).** `utils_workflow_enrichment.R` 1,134 -> 1,369 (+235) against the ticket's "+70 to +140"; the
  excess is the roxygen of the new internal functions. `utils_workflow_dia_da.R` 1,116 -> 1,154 (+38, inside
  "+25 to +50"); `utils_workflow_operations.R` 766 -> 760 (-6, against "-5 to +2").
- **Registration order (information).** `enrichment` is now registered after `DA-contrast-output`, so
  `workflowOperationNames()` lists it later; no test or consumer depends on the order (registry, partitions,
  binding-steps and the omics registry tests are identical).
