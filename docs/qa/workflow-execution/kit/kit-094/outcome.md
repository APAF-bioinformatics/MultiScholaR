# KIT-094 outcome: code part (2026-11-06)

**Scope of this run:** the code, the characterisation snapshot, the tests, the dependents grep and `R CMD check`, in an agent worktree branched from `performance-v2` at `86a32799`. **Not run:** the install, the `kit_094` build binding (`bindings-generator.R` and `bindings.json` untouched) and the two screen journeys `094a`/`094b`. They need KIT-092 (protocol v1.7, re-fingerprinted `lipid86s`/`metab86s`), which is still pending. They will be run later from the main checkout, under the ticket's decision table.

## What changed

The shared core now reads every per-step list from the binding. Lipidomics and metabolomics declare their current values.

| Binding field | Read by | Lipid/metab value |
|---|---|---|
| `steps` | `omicsWorkflowOperations` | `omicsRouteSteps`: the 13 declarations, same order |
| `worker_steps` | registry loop | the nine names, written literally per binding (`omicsWorkerSteps` removed) |
| `default_variants` | registry loop | `c(<prefix>_intensity_filter = "partitioned")` |
| `publications` | registry loop | `omicsRoutePublications`: the same nine, same handlers, same order |
| `filtering_steps` | compute handler | `<prefix>_intensity_filter`, `<prefix>_duplicate_resolution` |
| `native_steps` | compute handler | `<prefix>_da_contrast` |
| `dependency_parameters` | compute handler | `pairs` |
| `replay` | `replayWorkflowOmicsStep` | `omicsRouteReplay`: one entry per former `switch` arm, same wrappers and arguments |
| `record_outcomes` | `workflowRecordStepOutcomes` | `omicsRouteRecordOutcomes`: the former RUV and ITSD code, moved unchanged |

These are derived from the bindings, with no new field:
- The replay's input roles and regex come from `design_input_prefix`.
- The reopen allowlist, the `setWorkflowType` allowlist and its message, and the omics entries of `workflowRoutesByType` come from `omic_type`, through `omicsWorkflowTypes()`.
- `workflowRoutesByType` is now a function, so a binding injected by a mock reaches it. Its three call sites were updated.
- The proteomics entries stay literal.

Not changed, as the ticket decided:
- `workflowRecordReportSources`;
- `kitOperationIdentityProducts` (`tools/` is untouched);
- the name-keyed section, default and label tables.

`omicsDeclarationBindings()` still returns the two bindings. There is no registration API. No exported function changed.

## Evidence

- **Snapshot:** taken on `86a32799` before any `R/` edit by `characterise.R`. Two independent runs gave the same sha256 (`c69a6008…57bdfa8`). This is `tests/testthat/fixtures/kit094-characterisation.rds`.
- **Test first:** see `test-first.txt`.
  - Run 1: every stub test fails on the base commit, at registration.
  - Run 2 (`test-first-sites.R`, diagnostic): with only `steps` honoured, each generalised site fails on its own. The sites are worker steps, default variant, publications, filtering, native, dependency, replay, record hook, input role, route, reopen and session type.
- **After the change:** stub test plus characterisation 48/0/0/0, so the snapshot is `identical`. The frozen dump ran 68/0/0 with **0 skips**, and its test file and fixture have an empty `git diff` against the base. Lipid and metab are 259/0/0/0 each, with literal worker steps and a tamper case per input-role prefix. Dependents and targeted files are in `dependents.txt`: 12 files, 3522 passed, 0 failed.
- **Perturbation (shown once, reverted):** I dropped `lipid_log2` from the lipid `worker_steps`. The binding-steps file then gave 47 passed, 1 failed. Only the characterisation failed: `actual$lipidomics$lipid_log2$compute` was `"in_process"` against the expected `"omics"`. The file was restored with `git checkout`, and the tree was clean afterwards.
- **`R CMD check`** (`--no-tests --no-manual --no-build-vignettes`, `check.log`): 3 WARNINGs, 3 NOTEs. The same sections as the KIT-090 baseline, and the code-problems NOTE body is identical.

## Divergences and notes

- `workflowRecordStepOutcomes` now calls a binding's hook only when `isTRUE(startsWith(operation, "<prefix>_"))`. Before, a `NULL` operation raised an error in the `%in%` test, and that can no longer happen. The record builder always supplies a string, falling back to `"legacy"`, so there is no observable change.
- Under `load_all`, the owner's software identity hashes the body of `newWorkflowStateOwner`, which this ticket edits. The characterisation therefore removes entity ids, which carry that identity, as the ticket's contract allows.
- `docs/qa/workflow-execution/kit/kit-046/route-check.R:135` reads `workflowRoutesByType` as a value. It is a historical script outside the allowlist, recorded here and left unedited.

## Remaining for the ticket

1. KIT-092 lands.
2. Then, from the main checkout: install, bind `kit_094` in `bindings-generator.R`, and regenerate `bindings.json`.
3. Run `094a` and `094b` (one attempt each, 1,800 s) against the re-fingerprinted `lipid86s`/`metab86s`, and record the decision table's branch.

## The screen journeys (2026-11-06, after KIT-092)

`kit_094` = te0soO, installed from a worktree pinned at e6919066: the merged KIT-094 code, KIT-093's trx code and KIT-092's protocol v1.7, with no uncommitted work from other agents. One attempt each, ceiling 1,800 s.

| Member | Corpus | Oracle (v1.7) | Decision | Exact | Files | Mismatches | Time | Whole slice | Application |
|---|---|---|---|---|---|---|---|---|---|
| `094a` | lipid_screen_missing | `lipid86s` | accepted | yes | 208 | 0 | 195 s | 1.95 GiB | 0.94 GiB |
| `094b` | metab_screen_missing | `metab86s` | accepted | yes | 199 | 0 | 185 s | 1.95 GiB | 0.94 GiB |

Both are exact under the decision table's first row, so no attribution member was needed. The step-list generalisation changes nothing in a whole lipid or metab journey.

(The verifier's oracle label still reads "protocol v1.5", a fixed string in the harness. The oracles are the v1.7 re-fingerprints bound by KIT-092.)
