# KIT-228 step 0 (2026-11-26)

Base: `origin/v0.9-staging` 0bdf72131 (the ticket's r1.1 was written at 3deba98d). Nothing run; every statement below is by reading.

## Bounded work item

Outcome: a Summary Save that does not write a new integration input leaves no integration input for that omic
(`integration/<omic>/`): the earlier one is removed once, at the start of the Save, and the user is told (VAL-201 (1)).
Owning ticket KIT-228 r1.1 with Re-check 1 (R01 to R05). Files: `R/utils_workflow_summary.R` (the Save path only; Copy and
Report belong to KIT-219), `R/utils_workflow_integration_input.R` (the one text at the failed publish), `NEWS.md`, the
ledger row, the new test file `tests/testthat/test-kit228-integration-input-refused-save.R`, this folder. Checks: the
ticket's rows E1 to E6, G0 to G4 plus R01's text row, then R2. Exit: the rows pass on head, fail on base for their stated
reason, R2 shows no regression. KIT-277's part (the organism Save) is not implemented.

## Owner's answers (ledger VAL-201, read in `docs/methods-decision-ledger.md`)

(1) the earlier input is removed at the start of every Save, in one place for all four omics, and the user is told; a removal
that itself fails is named and the Save continues (fail closed). (2) Limit accepted: an omic's input is the one its last
successful Save wrote. (3) Merged before the release commit. The ticket's recommendations and the row agree.

## The identity script ("the I05 edits (now applied)")

DISPATCH-STATE (line 1336) says the merge check as written was blocked by the retired identity script until the I05 edits
were applied. Meaning: `docs/qa/workflow-execution/kit/kit-240/default-identity.R` was retired at KIT-267's merge (VAL-340); the
ticket's closing note (KIT-228.md:174) now reads "re-runs `tests/testthat/test-kit240-integration-default.R` beside its own
files". For my runs: no identity script is run; `test-kit240-integration-default.R` is in the comparison run. (KIT-228 changes
no integration-GUI or record file, so no textual neighbour is expected; KIT-240's R2 file list is included for the check.)

## What moved since the ticket (r1.1 at 3deba98d; head 0bdf72131)

| Ticket statement | At 0bdf72131 | Effect |
|---|---|---|
| Save at `R/utils_workflow_summary.R:124-199`; `tryCatch({` `:163`; head guard `:164` | `saveWorkflowSummaryParameters` `:254-335`; `tryCatch({` `:287`, `assertWorkflowSummaryHeadFinal` `:288`, then `assertWorkflowSummaryDaCurrent` (KIT-212; its class `workflow_summary_da_stale` is now in the branch with KIT-207's) | place unchanged: first statement inside `tryCatch({` |
| per-class removals `:135`, `:150` call `removeIntegrationInput(<folder>)` | `:267`, `:283` call `removeOmicIntegrationInput(dir, omicType)` (KIT-247, VAL-239 (iv); `R/utils_workflow_omics_design.R:443`): the omic's folder AND the top-level `<omic>_*` file the Save writes beside it (`lipidomics_<label>_final_s4.RDS`, `transcriptomics_normalized_data.tsv`); test KIT-247 E23 (`test-kit247-design-resave.R:404`) asserts a refused Save removes both | **changes what is built, see below** |
| KIT-216 I10 raises KIT-207's class for proteomics | `assertProtSummaryHeadFinal` (`:109-134`) is called from `assertWorkflowSummaryHeadFinal` | none |
| KIT-216 SD-4 inside `protSummaryIntegrationWrite` | `R/utils_workflow_dia_summary.R:51-75`; the writer removes the folder itself with no source (`:74`) | none |
| worker failure loses its class (`R/utils_workflow_execution.R:964`) | `fail(action, simpleError(workflowManagedFailureText(result, minutes)), result$status)` at `:1011` (and `:998`, `:1024`, `:1028`, `:1044`) | same mechanism, line moved |
| gate `:216-219` | `startWorkflowSummaryIntegrationExport` `:353-360`, `done(FALSE)` | line moved |
| trx refusal `trx_declarations.R:947` | `:990` and `:1070` (same text) | line moved |
| lipid/metab writers `omics_declarations.R:1854-1859`, `:1886` | `omicsSummaryIntegrationWrite` `:2591`, `writeOmicsIntegrationInput` `:2621-2625` (removes with no source) | line moved |
| retry path `utils_workflow_omics_actions.R:73-81` | `:73-82` (`runner$retrying()`) | none |
| `writeIntegrationInput` `:407-484`, text `:470`, `:478` | `R/utils_workflow_integration_input.R` `:594` ("it is kept"), `:599`, `:602` ("the previous input is kept") | R01's text is at `:602` |
| only caller `R/mod_omic_summary.R:247` | `:274` (`omicSummarySave`) | none |
| the integration workbook (T5) | no workbook ships (KIT-234, VAL-209; `inst/` has none) | T5 dropped, as the ticket's addendum says |
| KIT-243 (part D) | new assert `assertWorkflowDaPublicationComplete`, used by Copy/Report, not Save | none |

Nothing else of the ticket's section 1 and 2 table changed in kind: rows 1 to 12 still describe the paths.

### The one departure: which helper removes at the start

The ticket says `removeIntegrationInput(<folder>)` and that the per-class calls are deleted. Since KIT-247 the per-class calls
remove folder plus top-level file (E23 asserts it). Deleting them and removing only the folder at the start would stop
removing the top-level file on a refused Save and fail E23. So the start step calls `removeOmicIntegrationInput(integrationDir,
omicType)`: one removal, the omic's whole earlier integration input as KIT-247 defines it. Consequence: on every Save the
earlier top-level `<omic>_*` file is removed at the start and the Save rewrites it (in process: `saveRDSFn`; worker: the
action) when it succeeds; a leftover with another label disappears on any Save. That is fail-closed and inside VAL-201 (1)
("the earlier input is removed"); the integration reads only the folder (`R/utils_workflow_integration_input.R:17-19`), so no
computed value changes. This is recorded as a departure in outcome.md and put to the coordinator; if the coordinator rules
folder-only, the per-class calls must then stay, which contradicts "one removal". "Earlier" for the notices (T1 to T4) is the
folder's existence, as in the ticket.

## Re-check 1 repairs R01 to R05: what each changes in what I build or test

- **R01 (major).** `R/utils_workflow_integration_input.R:602` says "the previous input is kept" on every failed publish.
  Under the start removal `previous` is NULL, so it would be untrue. Build: the text at `:602` is conditional: without a
  previous input "the integration input could not be published at %s"; with one restored, as today. `:594` unchanged (reached
  only when a folder exists at publish, i.e. after a failed start removal). Test: row R01 in the new file (a failed publish
  with no previous input: no "kept"; failing first on base). `test-workflow-integration-input.R:98` has a previous input and
  stays as is. File added to the change list.
- **R02.** T3's text replaced by the register's; ticket section 3 item 8 (retry) and section 9 line are text only; no test
  row: the retry path is not exercised (not checked, listed). Build: T3 as worded in R02.
- **R03.** T4 replaced ("may still read it"). Build and G3 expectation use the register's text.
- **R04.** (a) T1 is emitted in both branches of `summaryFailed`; (b) T2 is keyed on the gate's branch (the proteomics writer also
  returns FALSE with no source: that is T3); build: the gate's `done(FALSE)` carries an attribute `gate`, read by the `done`
  wrapper; (c) citations only; the `:491` citation is moot at this head (`test-workflow-integration-input.R` line numbers
  moved) and is not used.
- **R05.** The workbook is not shipped, so T5 and its edit are gone; the workbook test files are not in R2 for the T5 reason
  (`test-integration-workbook*.R` match R2's grep only if they reference the listed names). Files changed: as in the bounded
  item above; no `NAMESPACE`, `man/` or version change.

## Predictions made at step 0 (before any run)

- Stored references predicted to move: none, provided G4 holds (ticket section 6), plus one thing the ticket did not see:
  tests that assert the top-level file of a Save that failed. KIT-247 E23 expects it removed on a refused Save: unchanged.
- `test-omic-summary-module.R` and the GUI-entry lanes: a Save with an earlier folder that ends without a new input is the
  only thing that changes; G0 is the inventory (see predictions.md).
