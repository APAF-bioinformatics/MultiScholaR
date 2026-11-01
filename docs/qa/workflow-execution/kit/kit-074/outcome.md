# KIT-074 outcome: omic-neutral names for the kit functions the omics code calls

Status 2026-11-01: **completed**. Stage 3b, ticket 3 of 10.

- **Commits** on `performance-v2`: code 4f4a6b02, bindings and records ec39d8bc, and this evidence commit.
- **Base:** 47ee3113 (KIT-073 closed).
- **Install:** `kit_074` = `core/browser-smoke.E9Srqc` (4f4a6b02), over the kit_073 library.
- **Ticket:** `tickets/kit/KIT-074.md`, revision 1.1.

## What changed

1. **Inventory.** [inventory.sh](inventory.sh) lists every DIA-named function in scope, with its neutral name and its callers by omic. It was run before the change ([inventory-before.tsv](inventory-before.tsv), [residual-before.txt](residual-before.txt)) and after it ([inventory.tsv](inventory.tsv), [residual-after.txt](residual-after.txt)). It has 19 rows:
   - **Rule (a), called by the omics code (6):** `saveWorkflowStateValue`, `recordWorkflowInputs`, `clearWorkflowStepProvenance`, `collectWorkflowNativeReferences`, `runWorkflowOwnedSessionExport`, `createWorkflowSessionOwner`.
   - **Rule (b), named by KIT-075 and KIT-078 (8):** `readWorkflowStateGuiMetadata`, `getWorkflowOwnedPayload`, `setWorkflowOwnedPayload`, `getWorkflowOwnedInput`, `resolveWorkflowOwnedValue`, `getWorkflowReactiveReference`, `assertWorkflowExportReady`, `publishWorkflowActionOutputs`.
   - **Already neutral since KIT-045 (5):** `saveWorkflowOwnedValue`, `publishWorkflowStateHead`, `withWorkflowStateScope`, `selectWorkflowStateReference`, `closeWorkflowStateOwner`.
2. **Neutral definitions.** Each of the 14 new names was placed at the function's existing definition, in `utils_workflow_dia_state.R`, `utils_workflow_dia_sessions.R`, `utils_workflow_dia_da.R` and `utils_workflow_saved_exports.R`. The bodies are unchanged and still call the DIA names internally. Each DIA name follows its definition as `<dia> <- <neutral>`, with the KIT-045 roxygen line.
3. **Omics call sites.** These now use the neutral names:
   - the core (`recordWorkflowInputs`, `collectWorkflowNativeReferences`, `publishWorkflowStateHead`, `clearWorkflowStepProvenance`, `saveWorkflowStateValue`);
   - the lipid and metab norm export sites (`withWorkflowStateScope`, `runWorkflowOwnedSessionExport`);
   - the two omic servers (`createWorkflowSessionOwner`).

   The record's omics branches had no DIA call. After the change, the residual list of DIA-named calls in the omics files is empty.
4. **Unchanged:**
   - proteomics call sites, the harness and the tests (except the two server tests below);
   - the store directory, entity kinds, envelope classes, the `DIAStateAdapter` class string (still tested by `inherits` in the core), the `dia.*` schemas and the `dia_capacity_manager` slot.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Aliases | `tests/testthat/test-workflow-neutral-names.R` | All 19 rows: the neutral name is a function and the DIA alias is `identical()` to it. 39 pass. |
| Affected tests | [affected-tests.txt](affected-tests.txt) | 178 files, one process per file: **16,799 pass, 0 fail, 0 error**. |
| R CMD check | [check-item1.log](check-item1.log) | Baseline: 3 WARNINGs, 3 NOTEs. The code note is line-for-line identical to KIT-073's. |
| Exactness set (1)(3)(6) | [exactness-gates.txt](exactness-gates.txt) | The `lipid74` and `metab74` records are valid with no legacy step, and **identical** to `lipid31` and `metab32`. The KIT-028 replay is **exact** (16 states). The LFQ, TMT and DIA_LIMPA records are **identical** to the kit_065 oracles. |
| Exactness set (2) | [074l](kit074-kit-journey-074l-stage.json), [074m](kit074-kit-journey-074m-stage.json) | Both accepted and exact, with 0 oracle mismatches. `074l`: 241 s, 2.53 GiB peak. `074m`: 240 s, 2.61 GiB peak. No sampler gaps. |
| Exactness set (4) | [074l](member-record-replay-074l.json), [074m](member-record-replay-074m.json) | Both replays are exact (10 states), and a tampered input is refused. |
| Exactness set (5) | [lipid](kit074-kit-lifecycle-lipid-stage.json), [metab](kit074-kit-lifecycle-metab-stage.json) | Both observed. All four members were launched one at a time. Lipid: 242 s, 2.54 GiB peak. Metab: 242 s, 2.61 GiB peak. |

**Memory (reported, not claimed).** The lipid lifecycle peak is back at 2.54 GiB, against kit_032's 2.53. This supports reading kit_073's 2.68 GiB as run-to-run variance.

## Deviations and discoveries

- **Files outside the ticket list.**
  - `R/utils_workflow_dia_da.R` holds `readDIAStateGuiMetadata`.
  - `tests/testthat/test-workflow-server-shared.R`: the lipid and metab server tests mock the owner constructor by name, so they follow the call site to `createWorkflowSessionOwner`. The two proteomics tests keep the DIA mock.
  - `bindings-generator.R` needed no edit: its stage-3b loop already covers 074.
- **Mocks.** A test that mocks a DIA name no longer intercepts an omics call made through the neutral name. KIT-075 and later tickets mock the neutral names for omics paths.
- **Tracing.** `tools/profiling/workflow_checkpoint_gui_diagnostic.R` traces `saveDIAStateValue` and `runDIAOwnedSessionExport` by name, so it no longer sees omics calls. It is a proteomics diagnostic, so I recorded this and left it untouched.
- **Pre-existing warning, recorded and left untouched.** Both omics recorders print the `copyToResultsSummary` warning "N files/objects/directories failed to copy correctly" (lipid 3, metab 2). The kit_073 build prints the same lipid warning. It belongs to KIT-077's results-summary area.

## Backward compatibility and rollback

- Every DIA name resolves to the same function.
- No saved format, directory or class string changed. Records and stores of earlier builds replay unchanged (KIT-028; the stage-3 oracles).
- Rollback: `git revert` the ticket's commits.

## Handoff

KIT-075 writes against the neutral names: the owned-payload helpers, `resolveWorkflowOwnedValue`, `getWorkflowReactiveReference` and `readWorkflowStateGuiMetadata`. A name a later ticket finds missing is added there under item 1(b) and recorded.
