# Stage 3b remediation, 2026-11-01

Closure record for [STAGE3B-AUDIT-20261101](STAGE3B-AUDIT-20261101.md); the audit itself is unchanged. User instruction: "fix all please" (2026-11-01). Mode `audit_and_repair` on the ticket texts only: no package, harness or protocol file was changed. `protocol.md` stays v1.5 until KIT-081 freezes v1.6 on the user's go.

## Decisions taken under "fix all" (the audit's recommended option in each case)

| Finding | Decision | Where it lands |
|---|---|---|
| S3B-19 (B): the worker's objective had no observable | KIT-078 proves the controller never materialises an assay with a hydration probe (implementation conformance) and makes no memory claim. KIT-081 samples the controller service's cgroup. KIT-080 binds two pre-declared rules, `whole_slice_anon` for total memory and `app_process` (≤ 0.95, the protocol's existing threshold) for the application process, each deciding a named claim. | KIT-078 item 7, KIT-081 item 2, KIT-080 item 3, protocol v1.6 proposal item 3 |
| S3B-23 (B): ITSD partition not exact with the unchanged kernel | ITSD removed from KIT-079 and stays `whole`. Partitioning it would need an additive kernel parameter, which is a separate user decision and not taken. | KIT-079 item 3; KIT-080 counts 3 partitioned steps |
| S3B-26 (B): operation pairs impossible as scoped | The tooling moves to a new ticket, KIT-081, between KIT-079 and KIT-066 (lists-assay input components, suffixed operation keys of two pairs each, bound ceilings). KIT-080's default switch is keyed on `reduction` AND `accepted`. | KIT-081 items 3–4, KIT-080 item 4 |
| S3B-27 (B): control and corpus outside v1.5 | A protocol v1.6 proposal (`protocol-v1.6-proposal.md`) names `kit_032` as the stage-3b control and allows real-corpus pairs with per-build in-process oracles. KIT-081 freezes it only on the user's go, before any counted member, and KIT-080 cannot start without it. | proposal items 1–2, KIT-081 item 1, KIT-080 blocked state |

## Closure per finding

| Finding | Repair | Evidence (revision 1.1 unless stated) |
|---|---|---|
| S3B-01 (H) | The stage-3b exactness set is defined once in the plan: oracles of the build, members against them, the cross-build comparison against the stage-3 oracles, the records, the lifecycle members, and the proteomics part. Every ticket references it, adds its `0NNl`/`0NNm` keys, and lists `bindings-generator.R`. | `00-KIT-PLAN.md` "Stage-3b exactness set"; KIT-072–081 files |
| S3B-02 (H) | Forward references kept. The order (KIT-073 creates the core) is in each ticket's dependencies, and KIT-075 notes the current lines "(the KIT-073 core after that ticket)". | KIT-075 item 2; index order |
| S3B-03 (L) | Sites and tests re-cited (`restoreLipidDaCurrentS4FromSession`; the exports already one path; the kernel-only tests excluded; the exact `02aa`/`04` files; the R6 stub named). | KIT-072 items 1, 3 |
| S3B-04 (M) | The omitted module sites, the DA Open fallbacks, the design import tolerance and the two shared session helpers added. A grep script before/after is required. | KIT-072 items 1–2, Simulation Setup |
| S3B-05 (M) | The omitted tests and lanes added. The pass-through tests are recorded as contract changes. | KIT-072 item 3, `files_to_modify` |
| S3B-06 (M) | The binding carries data and per-omic kernel adaptors (duplicate resolution). The full list of differences, the `…Fn` call-site names, both replay pairs, the labels and the corrected range are recorded. | KIT-073 items 1–4 |
| S3B-07 (M) | The scope rule is kept, and the examples are regenerated from the tree. The already-neutral names are listed and the forward clause (b) added. | KIT-074 item 1 |
| S3B-08 (H) | Per-wrapper `hydrate` and the module `$object` readers listed; the DA `theObject` lines named. `revertToState` is unchanged and its GUI callers move to `selectWorkflowStateReference`. `getDIAOwnedInput` is added, with a cross-omic misroute test. | KIT-075 items 1, 2, 6 |
| S3B-09 (H) | The residency walk names the real holders; `qcData` is removed. Dead slots are recorded and the receipt probe re-cited with its change. | KIT-075 item 3, Simulation Setup |
| S3B-10 (M) | The omitted reads added to items 3–5. | KIT-075 |
| S3B-11 (M) | The memory observable, expectation and outcomes declared before the run, including a stop-and-report if the peak is above `kit_032`. | KIT-075 Experiment budget |
| S3B-12 (H) | The equality rule frozen before outcomes: `identical()` for `~ 0 + group`, and ≤ 1e-12 absolute with identical calls, rows and order otherwise, recorded as a deviation. A non-orthogonal fixture is added. | KIT-076 Scientific contract |
| S3B-13 (M) | Failure isolation, ordering, `significant_counts` aggregation and the product-hash change declared. The metab `tryCatch` is cited. | KIT-076 items 3, 6 |
| S3B-14 (M) | Consumers listed. No F/`F.p.value` consumer exists, so neither is produced per contrast. No stage-5 reader exists. | KIT-076 item 5 |
| S3B-15 (L) | The kernel halves re-cited. The omics fit identity (`newWorkflowOmicsDaFitIdentity`) is new code; the proteomics identity is not reused. | KIT-076 Authority, item 2 |
| S3B-16 (H) | A new omics figure product (`R/utils_workflow_omics_figures.R`): plot data for the four rebuildable types, rendered bytes for every figure (grids and PDFs bytes only). The proteomics functions are untouched and their identities byte-identical. A PNG determinism check is added. | KIT-077 item 1, Simulation Setup |
| S3B-17 (M) | Defects re-stated precisely. The metab repair is "look for the written name"; the lipid repair is "delete the fallback". The fixtures that hide the defect are corrected. | KIT-077 item 3 |
| S3B-18 (M) | The PDFs, the combined PDFs, the composite's three files and `updateMetaboliteFiltering` added. The Glimma interactive volcanos are out of scope, recorded. | KIT-077 items 2, 4 |
| S3B-19 (B) | See the decisions above. | KIT-078, KIT-081, KIT-080, proposal |
| S3B-20 (H) | The `gui`-mode managed proteomics member `078p` after each chain item and at closure. | KIT-078 item 9; plan exactness set item 6 |
| S3B-21 (H) | Each contract decided:<br>• the section rule changed (`section = NULL` with provenance-resolvable parameters), recorded as a pin change;<br>• deadlines kept, so no pin change;<br>• `input_references` unchanged, with ids and contrasts as parameters and the pairs as an entity-id parameter (the Pearson kernel reads the entity);<br>• product-only results return the source checkpoint;<br>• publication registers through `publishDIAActionOutputs`, and the export-ready state list comes from the omic's entries;<br>• the peptide-shaped completion functions are left alone, with omics getting its own. | KIT-078 items 1, 3, 4 |
| S3B-22 (L) | `:80-82` unchanged; `:793` unchanged; the driver re-cited (`:536-603`); the managed launch described as already true, with dispatch as the gap. | KIT-078 items 1, 3, 6 |
| S3B-23 (B) | See the decisions above. | KIT-079 |
| S3B-24 (M) | Merge axes per step (row-binds in source or kernel order; ITSD out). The design is carried; counts are recomputed after the merge. | KIT-079 item 2 |
| S3B-25 (L) | The codec tests re-cited (`test-workflow-saved-projection.R`, `test-perf-103*.R`). | KIT-079 Regression Checks |
| S3B-26 (B) | See the decisions above. | KIT-081, KIT-080 |
| S3B-27 (B) | See the decisions above. | proposal, KIT-081, KIT-080 |
| S3B-28 (M) | Stated as KIT-080's blocked state (corpora, v1.6, the user's go). | KIT-080 "Blocked state" |
| S3B-29 (L) | The consequence recorded in the plan's decision row and in KIT-066. | plan 3.2 decision row; KIT-066 note |

## Re-audit

- Ticket schema validation after the repairs: KIT-066 and KIT-072–081 pass both validators. The selector returns KIT-072. The inventory finds 79 tickets, no cycle, no duplicate and no missing dependency, with every new ticket status-consistent and carrying all required sections. Its 2 remaining "missing links" are the known false positives (`[[omic]]`, `[[field]]`, R subscripts).
- Not re-run: the semantic reference passes. The repairs cite lines taken from those passes, each already verified at its cited line, and add no new reference that was not checked there except the new files the tickets create.

## States

- `implementation_verified`: false. Nothing is implemented.
- `execution_complete`: true for the repairs.
- `decision_supported`: true for the ticket texts. The four decisions above were taken under the user's instruction, and the protocol v1.6 text still needs the user's go at KIT-081.
- `promotion_authorized`: false.
- `claim_scope`: the stage-3b ticket texts, plan 3.2.
