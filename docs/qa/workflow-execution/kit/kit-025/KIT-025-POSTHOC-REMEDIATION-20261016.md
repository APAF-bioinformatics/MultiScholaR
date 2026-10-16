# KIT-025 post-hoc audit remediation (2026-10-16)

Audit: [KIT-025-POSTHOC-AUDIT-20261016.md](KIT-025-POSTHOC-AUDIT-20261016.md) (unchanged). User: "fix all findings".

| ID | Repair | Closure evidence |
|---|---|---|
| PH-01 | `clearDIAStepProvenance(manager)` after the DA wrapper's contrast loop (the one product-only user of `runSavedOperation`; the enrichment runs `runProjected` and sets nothing). | `test-workflow-saved-da.R`: after the DA wrapper nothing is pending and the head state's entry stays what produced it; `test-workflow-dia-products.R`: a commit after a cleared slot is `legacy`. |
| PH-02 | `workflowSectionOperations` maps every registered worker step's kernel section to its operation; the three in-process finishing steps (`finishDIAReferencePeptideStep`, `finishDIAReferenceQvalueStep`, `finishDIAReferenceProteinStep`) name the operation from the section unless the worker completion already set one. | `test-workflow-operation-registry.R`: the map covers exactly the worker operations; `test-workflow-dia-products.R`: a section names the operation, a pending operation wins; the check script's in-process replicate filter on the fixture. |
| PH-03 | The check runs the protein replicate filter in process on the fixture from `duplicates_removed` and verifies every recorded parameter name is a name in the state's `slots$args[[section]]` or `config`. | `provenance-check-screen.json`: entry `protein_replicate_filter` (`rows`), `grouping_variable` and `parallel_cores` both recorded in the state (`in_process_names_recorded` TRUE, TRUE). |
| PH-04 | `workflowSoftwareIdentity()` memoised per process. | Code; the products test's identity assertions unchanged. |
| PH-05 | The scan follows `do.call` targets; the remaining limits stated in the outcome. | `test-workflow-operation-registry.R`. |
| PH-06 | Ticket revision 1.4: allowlist corrected; the superseded install recorded. | Ticket text. |
| PH-07 | Plan rule: the four pinned files never in an allowlist; installs through `run_workflow_kit.sh install`. | `00-KIT-PLAN.md`. |

All seven findings closed. Full suite at `c25de1ee` (clean worktree): 343 files, 29,724 pass, 0 fail, 140 skip, the one known PERF-101 runner artefact; `R CMD check` at baseline (3 warnings, 3 notes); logs in `check-logs/`. The remediated build is installed as `kit_025` (installed source commit `c25de1ee` verified) and the bindings regenerated at `d500bc72`.
