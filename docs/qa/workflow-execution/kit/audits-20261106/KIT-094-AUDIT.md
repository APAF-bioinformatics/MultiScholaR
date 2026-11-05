# KIT-094 audit (audit_only, 2026-11-06)

**Target:** `tickets/kit/KIT-094.md` rev 1.0, status `pending`, `needs_planning: true`.
**Mode:** `audit_only`. Nothing was repaired, edited, run or committed apart from this file. No tests, installs or journeys were run. Nothing under `kit/runtime` and no KC data was read.
**Baseline:** branch `performance-v2` at `bbe1fae9`. The tree is dirty only with four untracked `kit-090/*.txt` files, which belong to KIT-090.
**Authorities read:**
- `AGENTS.md` and `conventions/R.md` (scope rule, "Minimal", functions, testing);
- `tickets/kit/00-KIT-PLAN.md` v3.8 (rows 231–240, decision 2026-11-06);
- the approved plan `~/.claude/plans/write-a-series-of-atomic-dragon.md`;
- `docs/qa/workflow-execution/kit/protocol.md` v1.6;
- neighbours KIT-040, KIT-089 (outcome), KIT-090, KIT-092, KIT-093, KIT-095 to KIT-099.

**Validator:** `next-openai-ticket.sh KIT-094` passes: valid, 0 errors, 0 warnings, status consistent. The dependency `KIT-093` matches the index. Planning is explicit.

## Finding register

| ID | Severity | Dimension | Blocked state |
|---|---|---|---|
| KIT094-F01 | should-fix | 1 Completeness | implementation |
| KIT094-F02 | should-fix | 1 Completeness | implementation, execution |
| KIT094-F03 | should-fix | 2 Scientific rigor | execution |
| KIT094-F04 | should-fix | 2 Scientific rigor | decision |
| KIT094-F05 | should-fix | 2 Scientific rigor | implementation |
| KIT094-F06 | should-fix | 2 Scientific rigor | decision |
| KIT094-F07 | should-fix | 3 Best practice / 2 independence | implementation |
| KIT094-F08 | should-fix | 3 Test-first | implementation |
| KIT094-F09 | minor | 3 Smallest sufficient | implementation (planning) |
| KIT094-F10 | minor | 3 Ticket hygiene | none |
| KIT094-F11 | minor | 4 Reference accuracy | none |

### KIT094-F01: hard-coded step-name sites in the shared core are not listed (should-fix)

**Evidence.** These sites name lipid/metab step suffixes for every binding, and none is in the ticket's list:
- `R/utils_workflow_operations.R:697`: the registry loop sets `default_variant` to `partitioned` only for `paste0(prefix, "_intensity_filter")`.
- `R/utils_workflow_operations.R:703-712`: one fixed list of nine publication operations, the same for every binding. KIT-095 needs a publication the list doesn't have ("CCA diagnostics").
- `R/utils_workflow_omics_declarations.R:1138,1144,1148`: the worker's compute handler hard-codes three things:
  - the `pairs` dependency;
  - native references for `_da_contrast` only;
  - filtering metrics for `_intensity_filter` and `_duplicate_resolution` only.

  Transcriptomics needs filtering metrics for `trx_expression_filter` (KIT-096 QC tab).
- `R/utils_workflow_omics_declarations.R:1321-1352`: `replayWorkflowOmicsStep` is a `switch` over the lipid/metab suffixes. Any other step aborts `workflow_replay_unsupported`. Both the record's replay (`R/utils_workflow_analysis_record.R:580-583`) and KIT-095's "replay exact" reach it.

**Impact.** "Each binding declares its own steps" is only partly delivered. KIT-095 would have to reopen the core (its allowlist does include the file, but the generalisation is KIT-094's outcome) or special-case `trx_` names. KIT-094's stub test would not find this gap unless it covers these paths.

**Repair.** Add the four sites to the bounded outcome as binding-declared data:
- the default variant per step;
- the publication list;
- the filtering-display steps;
- replay dispatch, as a per-step replay entry or a binding function.

Alternatively, assign them to KIT-095 by name in both tickets.

**Closure evidence.** A revised ticket lists the sites. The stub-binding test (F08) exercises each site and fails before the change.

### KIT094-F02: `files_to_modify` omits files the outcome needs (should-fix)

**Evidence.**
- The binding objects the new fields go on are defined in `R/utils_workflow_lipid_declarations.R:17-92` (`lipidomicsDeclarationBinding`) and the metab equivalent. Neither file is in `files_to_modify`. Today the bindings have no `steps` field.
- The two journeys need a build binding (`b$builds$kit_094`) and member entries in `docs/qa/workflow-execution/kit/bindings-generator.R`, regenerating `bindings.json`. KIT-089 listed this file for the same kind of members (`tickets/kit/KIT-089.md:15`). KIT-094 does not.
- `tests/testthat/` is listed as a whole directory, which is broad under AGENTS.md's allowlist rule.

**Impact.** The implementation can't add `binding$steps` within the allowlist, and the journeys can't be declared or bound.

**Repair.**
- Add both binding files, `bindings-generator.R` and `bindings.json`.
- Replace `tests/testthat/` with the named files:
  - `test-workflow-omics-declarations.R`;
  - `test-workflow-lipid-declarations.R`;
  - `test-workflow-metab-declarations.R`;
  - `test-workflow-operation-registry.R`;
  - `test-workflow-analysis-record.R`;
  - the new stub test.

**Closure evidence.** The revised frontmatter validates, and the selector shows the new list.

### KIT094-F03: the journey control and its validity on the current build are not specified (should-fix)

**Evidence.**
- "Both screen journeys exact against their existing oracles" names no oracle, corpus, build label or protocol version. Two lipid/metab screen oracle families exist:
  - `lipid82`/`metab82` (complete screen corpora, `kit_082`);
  - `lipid86s`/`metab86s` (`*_screen_missing`, `kit_086`) (`bindings-generator.R:397-398, 432-445`).
- Since KIT-089 the intensity filter defaults to `partitioned`. Against `lipid86s`/`metab86s`, the verifier marked both KIT-089 members `invalid` on `products/names`, although their content was identical (`kit-089/outcome.md`).
- KIT-092 is the fix. It re-fingerprints "at least `lipid86s`, `metab86s`, `metab86`, `lipid86`", not `lipid82`/`metab82`.
- KIT-094 reaches KIT-092 only through KIT-093 → KIT-040 → KIT-092.

**Impact.** Run as written, the members may be `invalid` for a reason unrelated to KIT-094, or be judged against an oracle that was never re-fingerprinted. "Exact" isn't decidable from the ticket text.

**Repair.** In the experiment budget, name:
- the oracles: `lipid86s`/`metab86s` as re-fingerprinted by KIT-092;
- the corpora: `lipid_screen_missing`/`metab_screen_missing`;
- the members, for example `094a`/`094b`, on build `kit_094`;
- the protocol version in force after KIT-092.

Add `KIT-092` to `dependencies` explicitly.

**Closure evidence.** The revised ticket, and the bindings-generator entries resolving to existing oracle keys (`b$oracles$lipid86s_import` present).

### KIT094-F04: no attribution rule or terminal outcomes for a mismatch (should-fix)

**Evidence.**
- The last lipid/metab journeys ran on `kit_089`.
- After that, KIT-090 replaces NIPALS. NIPALS is used by the lipid/metab RUV and plot paths (`R/func_lipid_s4_ruv_methods.R`, `R/func_metab_norm.R`, `R/func_metab_s4_ruv_methods.R`, `R/func_*_plotting_methods.R`). KIT-090 is proved by component tests only ("No journeys", `KIT-090.md` budget).
- KIT-093 then changes DESCRIPTION and the class registry.
- The KIT-094 members are therefore the first end-to-end lipid/metab check of everything since `kit_089`. The missing-value corpora exercise NIPALS.
- The budget allows "one attempt each … one bounded correction pass" but defines no outcomes: exact, mismatch, or invalid instrumentation.

**Impact.** A mismatch can't be attributed to KIT-094 rather than KIT-090 or KIT-093. Without a decision rule, the correction pass could be spent re-running instead of diagnosing.

**Repair.** Declare the terminal outcomes before running:
- exact and accepted: pass;
- `invalid` (instrumentation): stop, and record it as distinct from a mismatch;
- mismatch: the correction pass may be spent on one member on the pre-KIT-094 commit (same corpus and oracle) to attribute it; otherwise stop and report.

**Closure evidence.** The decision table is in the ticket before the run, and the evidence folder records which branch was taken.

### KIT094-F05: "no behaviour change" is not falsifiable for the lists that move (should-fix)

**Evidence.**
- What the frozen dump covers: registry entries and declarations (kernel excluded) and the two in-process wrapper journeys (`helper-workflow-omics-declarations.R:20-53`).
- What it does not cover: any of the five lists the ticket moves:
  - `workflowRecordStepOutcomes` (`R/utils_workflow_analysis_record.R:118-160`);
  - the replay input-role filter (`:463-466`);
  - `workflowRecordReportSources` (`:735-779`, ordered: "the first the record carries is taken", `:738-739`);
  - the type allowlists (`R/utils_workflow_dia_state.R:1330,1410`; `workflowRoutesByType`, `R/utils_workflow_dia_sessions.R:817`);
  - `kitOperationIdentityProducts` (`tools/profiling/workflow_kit.R:729`).
- The replay input filter fails silently. A role it doesn't recognise is skipped (`next`), not refused (`:464-466`). If a refactor dropped a lipid/metab role, sha256 checking would stop for that input without failing a test, unless a tamper test covers that role. `test-workflow-lipid-declarations.R:216` covers only the tampered file it uses.

**Impact.** The acceptance can pass while record, report or replay behaviour has changed.

**Repair.** Before editing, capture a characterisation snapshot on the base commit:
- the realised `workflowRecordReportSources`;
- the allowlist and route sets;
- a fixed role vector through the input filter;
- the identity regex applied to a fixed product-name vector;
- the lipid/metab fixture records' steps and outcomes, with software identity removed.

After the change, assert each one is `identical`. Add a tamper case per omic role prefix.

**Closure evidence.** The snapshot script and its output are in `kit-094/`, the test fails on a deliberately perturbed list, and it passes after the change.

### KIT094-F06: the frozen-dump acceptance can pass by skip or by editing its normalisation (should-fix)

**Evidence.**
- The test skips when `workflowOwnerCapable()` is false or the fixture is missing (`test-workflow-omics-declarations.R:7,9`).
- Its lines 13-40 already strip the later tickets' changes (KIT-076/078/079/089). The `tests/testthat/` allowlist would permit more stripping.
- `test-workflow-lipid-declarations.R:25` and `test-workflow-metab-declarations.R:25` use the global `omicsWorkerSteps`, which the ticket makes per-binding. They must change, and they could be changed to agree with whatever the code does.

**Impact.** The acceptance evidence could be weakened inside the ticket that it is meant to judge.

**Repair.** The acceptance requires:
- reporter output showing the frozen-dump test ran with 0 skips;
- an empty `git diff` for `tests/testthat/fixtures/kit073-omics-declarations-frozen.rds` and for lines 13-40 of that test;
- the lipid/metab tests asserting the nine worker steps as literal names, not through the new field.

**Closure evidence.** The reporter log and the diff, in `kit-094/`.

### KIT094-F07: moving `kitOperationIdentityProducts` into a binding field lets the build under test set the verifier's ignore-list (should-fix)

**Evidence.**
- `tools/profiling/workflow_kit.R` is sourced standalone ("Pure functions only", lines 1-6; `sys.source` in `test-workflow-kit-harness.R:7`). It has no package references (`MultiScholaR::`, `asNamespace`: none).
- The regex decides which products the cross-build verifier does not judge (`:1003, 1531, 1535`).
- Deriving it from the candidate build's binding means the build being verified chooses what is exempt from verification.
- KIT-098 already owns `workflow_kit.R` for the transcriptomics route.

**Impact.** The verifier is no longer independent of the build it judges, and the harness becomes coupled to package internals.

**Repair.** Keep the regex as a harness constant and remove it from KIT-094's list. KIT-098 extends it for `trx_*_science` under protocol v1.7.

**Closure evidence.** The revised KIT-094 and KIT-098 texts, and `workflow_kit.R:729` unchanged by KIT-094.

### KIT094-F08: the test-first stub test is under-specified (should-fix)

**Evidence.**
- Simulation Setup says only "a test of a stub third binding with its own steps fails before the change".
- `omicsDeclarationBindings()` is a constant function (`R/utils_workflow_omics_declarations.R:580-582`), and `test-workflow-omics-declarations.R:60` asserts exactly the two bindings. Registering a stub in the package would break that test, or need a registration API that no consumer requires.
- A test of `omicsWorkflowOperations(stub)` alone would not falsify F01's sites or the five moved lists.

**Impact.** The test-first evidence could prove only the easiest site.

**Repair.**
- Inject the stub with `testthat::local_mocked_bindings(omicsDeclarationBindings = ...)`.
- Assert that the stub's own steps and worker steps register.
- Assert that its default variant, publications, record outcome hook, input role, workflow type and route are honoured, and that replay dispatches its steps.
- Add no public registration function.

**Closure evidence.** `kit-094/test-first.txt` shows the failures before and the pass after.

### KIT094-F09: the per-list decision between a binding field and an additive entry is not justified (minor)

**Evidence.**
- `workflowRecordReportSources` is one ordered table of report variables to record paths, shared with proteomics (`:735-779`). Distributing it over bindings needs reassembly that preserves order.
- KIT-095 also lists "record outcomes and report sources" in its registry outcome, so it isn't clear which ticket owns the mechanism and which owns the content.
- On the other side, KIT-095's allowlist omits `utils_workflow_dia_state.R` and `utils_workflow_dia_sessions.R`, so the type allowlists do have to be derived from bindings in KIT-094. This is consistent.

**Impact.** There is a risk of machinery beyond what's needed (`conventions/R.md` "Minimal"; AGENTS.md's smallest-sufficient rule), and the KIT-094/095 split is unclear.

**Repair.** In the planning step, record for each list whether it becomes binding-derived or takes an additive `trx` entry in KIT-095, with the reason. KIT-095 then states that it adds content only.

**Closure evidence.** The KIT-094 file-by-file plan names each list's disposition.

### KIT094-F10: the "Scientific contract" section is the series' transcriptomics contract (minor)

**Evidence.** `KIT-094.md:34-36` repeats the notebook, edgeR and count contract. None of it applies to this ticket. KIT-094's own contract appears only in the acceptance and backward-compatibility sections.

**Impact.** The claim scope is unclear to the implementer and the reviewer.

**Repair.** Replace the section with KIT-094's own contract: lipid/metab registry, declarations, records (software identity aside), replay and journey content unchanged; no oracle, fixture or protocol change.

**Closure evidence.** The ticket revision.

### KIT094-F11: small reference and plan inconsistencies (minor)

**Evidence.**
- The approved plan cites the registry loop as `R/utils_workflow_operations.R:682-718`. The loop is 682-713, `invisible(...)` is at 714, and the function closes at 715.
- The approved plan's run table says KIT-094 is "Tests only". Its KIT-094 body, the ticket and `00-KIT-PLAN.md:233` say "two screen journeys".

**Impact.** Cosmetic. The journeys are the authoritative intent.

**Repair.** Correct the plan text when the ticket is next revised.

**Closure evidence.** The plan diff.

## Checked and excluded

- `workflowOmicsSessionClass` (`R/utils_workflow_dia_sessions.R:1318`) is legacy-session adoption only. Transcriptomics has no legacy sessions, and an unknown omic is refused by class (`:1340`).
- `setupDirectories` already accepts `transcriptomics` (`R/func_general_filemgmt_setupdirectories_helpers.R:25,78`).
- The harness route lists belong to KIT-098 or KIT-040, which list those files:
  - `workflow_kit.R:15,22,111,178,251,482,622`;
  - `workflow_kit_coordinator.R:18`;
  - `workflow_kit_driver.R:196`;
  - `kit_omics_corpus.R:22`.
- `workflowSectionOperationsByType`, `workflowOperationDefaults` and `workflowParameterLabels` are additive tables keyed by name, filled in KIT-095.
- The record's DA binding lookup (`R/utils_workflow_analysis_record.R:604`) already works through `omicsDeclarationBindings()` by the `_da_contrast` suffix.

## Reference verification (dimension 4)

Every reference below was checked against the code at `bbe1fae9` and is accurate:
- `omicsWorkerSteps` (`R/utils_workflow_omics_declarations.R:13`);
- `omicsWorkflowOperations` (`:567`, 13 declarations, matching "keep their 13");
- `omicsDeclarationBindings` (`:580`);
- `registerBuiltinWorkflowOperations` (`R/utils_workflow_operations.R:520`, binding loop `:682`);
- `workflowSectionOperationsByType` (`:55`), `workflowOperationDefaults` (`:265`), `workflowParameterLabels` (`:383`);
- `workflowRecordStepOutcomes` (`R/utils_workflow_analysis_record.R:118`, lists at `:143,158`), input roles (`:465-466`), `workflowRecordReportSources` (`:735`);
- `kitOperationIdentityProducts` (`tools/profiling/workflow_kit.R:729`);
- the allowlists (`R/utils_workflow_dia_state.R:1330,1410`);
- `workflowRoutesByType` (`R/utils_workflow_dia_sessions.R:817`);
- the frozen fixture `tests/testthat/fixtures/kit073-omics-declarations-frozen.rds` (exists);
- the `R CMD check` baseline of 3 WARNINGs and 3 NOTEs (`kit-089/outcome.md`).

The ticket itself cites no line numbers. The only line inaccuracy is the plan's (F11).

## Verdict per dimension

1. **Completeness: insufficient as written.** All five listed sites are real. However, four core sites that a third binding hits are missing (F01), and the allowlist omits the binding files and the bindings generator (F02).
2. **Scientific rigor: not yet falsifiable.** The frozen dump is real, but it covers only registry, declarations and in-process journeys, and it can skip (F05, F06). The journey control has no named oracle or build. It is valid only after KIT-092 (F03), and it has no attribution rule against the KIT-090 and KIT-093 changes it would be the first to exercise end to end (F04).
3. **Best practices: acceptable intent, with two risks.** The generalisation serves a named current consumer (KIT-095), under the user's decision (1). The risks are the verifier's ignore-list coming from the build under test (F07), and moving lists where an additive entry would do (F09). The test-first stub needs a mechanism and coverage (F08).
4. **Reference accuracy: accurate.** Every cited function and file exists as described. The only discrepancies are in the plan text (F11).

**Overall:** no blocking finding. Eight should-fix findings gate implementation (F01, F02, F05, F07, F08), execution (F02, F03) and the acceptance decision (F04, F06). They can be resolved in the ticket's mandatory planning step without new authority, except adding KIT-092 as an explicit dependency (F03), which is an index edit.
