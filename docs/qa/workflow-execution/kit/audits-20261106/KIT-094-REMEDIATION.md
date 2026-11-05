# KIT-094 remediation (audit_and_repair, 2026-11-06)

**Authority:** the user, 2026-11-06: "fix all recommendations", relayed by the coordinator with five decisions (decision 3: parameterise rather than copy; decision 5: KIT-092 is protocol v1.7).
**Audit:** [KIT-094-AUDIT.md](KIT-094-AUDIT.md), left unchanged.
**Ticket:** `tickets/kit/KIT-094.md` revision 1.0 → 1.1, date-updated 2026-11-06. Only this file and the ticket were edited.

## Index changes for the coordinator to sync

- **Dependencies:** `["KIT-093"]` becomes `["KIT-093", "KIT-092"]` (F03). Until the index is synced, the selector reports `frontmatter dependencies do not match index dependencies`; this is the only validation error. The OpenAI validation is valid, the status is consistent, and there are 4 acceptance items.
- **Title:** unchanged.
- **Description:** revised in the frontmatter. Sync it if the index mirrors descriptions.
- **Time estimate:** 5h → 6h, for the added core sites.
- **Plan text (F11, not editable here):**
  - the approved plan's registry-loop citation `R/utils_workflow_operations.R:682-718` should read `:682-713`;
  - its run table's "Tests only" for KIT-094 should read "two screen journeys (094a, 094b)";
  - its KIT-094 list still says `kitOperationIdentityProducts` and `workflowRecordReportSources` move to binding fields. Revision 1.1 decides otherwise (F07 and F09: a harness constant extended in KIT-098, and an ordered table with KIT-095 rows).

## Closure register

| Finding | Severity | Disposition | Where in the ticket | Remaining |
|---|---|---|---|---|
| KIT094-F01: four hard-coded core sites missing | should-fix | **Closed.** The default variant, the publication list, the compute handler's three checks and the replay switch are now binding fields: `default_variants`, `publications`, `filtering_steps`, `native_steps`, `dependency_parameters`, `replay`. | Outcome, field table | Implementation |
| KIT094-F02: allowlist omissions | should-fix | **Closed.** Added `R/utils_workflow_lipid_declarations.R`, `R/utils_workflow_metab_declarations.R`, `bindings-generator.R` and `bindings.json`. `tests/testthat/` is replaced by the named test files and the new fixture. `tools/profiling/workflow_kit.R` is removed (F07). | Frontmatter | — |
| KIT094-F03: journey control unnamed; validity after KIT-092 | should-fix | **Closed.** Dependency on KIT-092 made explicit. The members are named: `094a` (lipid_screen_missing, oracle `lipid86s`) and `094b` (metab_screen_missing, oracle `metab86s`), on build `kit_094`, 1,800 s each, under protocol v1.7 with the oracles as KIT-092 re-fingerprints them. | Frontmatter, experiment budget | The index sync |
| KIT094-F04: no attribution or terminal outcomes | should-fix | **Closed.** A decision table is fixed before any run: exact; invalid (one rerun after the cause is shown); mismatch (one attribution member `kit_094c` on the pre-KIT-094 commit, then stop either way). | Experiment budget | — |
| KIT094-F05: moved lists not falsifiable; silent input-role skip | should-fix | **Closed.** A characterisation snapshot is captured on the base commit before any edit (`fixtures/kit094-characterisation.rds`), covering the report sources, allowlists and routes, the input-role vector, the fixture records' steps and outcomes, and the touched registry fields. It must be `identical` after the change, and a perturbation must fail it. Tamper cases are added per role prefix. | Simulation setup, acceptance | — |
| KIT094-F06: frozen dump can skip or be edited | should-fix | **Closed.** The test must run with 0 skips, shown by reporter output. `git diff` must be empty for `test-workflow-omics-declarations.R` and the KIT-073 fixture, so both are kept out of the allowlist. The lipid and metab tests assert the nine worker steps literally. | Acceptance | — |
| KIT094-F07: identity regex in a binding would let the build set the verifier's ignore-list | should-fix | **Closed by routing** (agreed by the coordinator). `kitOperationIdentityProducts` stays a harness constant, and KIT-098 extends it. KIT-094 does not touch `tools/profiling/`. | Outcome, dispositions | KIT-098 to carry it (coordinator) |
| KIT094-F08: stub test mechanism and coverage | should-fix | **Closed.** The stub is injected with `local_mocked_bindings(omicsDeclarationBindings = ...)`, with no registration API. The test asserts every generalised site, including the allowlists, route, input role, record outcomes and replay. `omicsDeclarationBindings()` still returns exactly two bindings. | Simulation setup, backward compatibility | — |
| KIT094-F09: binding field versus additive entry, by list | minor | **Closed.** Each list's disposition is decided in the ticket. **Binding fields:** steps, worker steps, default variants, publications, the compute handler's checks, replay and record outcomes. **Derived from existing fields:** input roles (from `design_input_prefix`), type allowlists and routes. **Unchanged tables with additive KIT-095 rows:** report sources, sections, defaults and labels. **Harness constant:** the identity regex. | Outcome, dispositions | The plan text (see the index changes above) |
| KIT094-F10: wrong scientific contract | minor | **Closed.** Replaced with KIT-094's own no-change contract. | Scientific contract | — |
| KIT094-F11: plan citation and run-table inconsistency | minor | **Routed to the coordinator** (the plan is out of this pass's scope). The ticket now cites the correct lines (`:683`, `:697`, `:703-712`). | Field table | Plan edit (coordinator) |
| KIT095-F04 (part, routed in): compute handler checks and replay switch | — | **Closed** in this ticket as the `filtering_steps`, `native_steps`, `dependency_parameters` and `replay` fields (`R/utils_workflow_omics_declarations.R:1138,1144,1148,1321-1352`). | Field table | — |
| KIT095-F07 (routed in): per-binding publication list | — | **Closed.** The `publications` field: lipid and metab keep the same nine with the same handlers, and transcriptomics adds CCA diagnostics in KIT-095 without changing them. | Field table | — |
| KIT095-F15 / KIT099-F13 (routed in): the `*_intensity_filter` default-variant rule | — | **Closed.** The `default_variants` field replaces the name rule at `R/utils_workflow_operations.R:697`. Lipid and metab keep `partitioned` on the intensity filter. KIT-099 sets the transcriptomics defaults through this field. | Field table | — |

**Counts:** 15 rows (11 audit findings and 4 routed-in). Of the 11 audit findings:
- 9 closed in the ticket;
- 1 closed by routing to KIT-098 (F07);
- 1 routed to the coordinator for plan text (F11).

All 4 routed-in findings are closed. None is unclosable. The one validation error clears when the index takes the new dependency.
