# KIT-098 audit (audit_only), 2026-11-06

**Target:** `tickets/kit/KIT-098.md` revision 1.0, status `pending`, `needs_planning: true`, depends on KIT-097.
**Mode:** `audit_only`. Nothing was repaired, run, installed or committed. This file is the only write.
**Baseline:** branch `performance-v2` at `bbe1fae9`. Untracked files present beforehand: four `kit-090/*.txt`, not touched.

Short sha256 of the files read: `KIT-098.md` 9fb6ecbe682379f4, `workflow_kit.R` c0a0ff3e7549f5c7, `workflow_kit_driver.R` d5fb1a2e5becc352, `workflow_kit_coordinator.R` c58e1b921af67797, `run_workflow_kit.sh` fe209d6a852137a9, `bindings-generator.R` 18f5e9bcb75c7671, `protocol.md` 1531ad5b75e7d0b5 (equal to `bindings.json` `harness.protocol_sha256`), notebook `transcriptomics_workflow.rmd` 44d5f643946212c1.

**Read for context:** AGENTS.md; the four `conventions/R*.md` files (listed only, since no R code is under review); KIT-040, KIT-092 to KIT-097, KIT-099 and KIT-088; plan 3.8; the approved plan `write-a-series-of-atomic-dragon.md`; `protocol.md` v1.6; the kit README; the KIT-089 outcome; the KIT-031 and KIT-032 recorders; `kit-053/stage-summary.R`; and the harness test file.

**Mechanical checks:** `next-openai-ticket.sh KIT-098` validates: every required section is present, status is consistent and there are no warnings. The `audit_inventory.py` check finds no cycles, duplicates or missing dependencies across the 97 tickets. Its one link finding is in the plan (`omic`) and has nothing to do with KIT-098. In the index, KIT-098 is `pending` and depends on `KIT-097`, the same as the ticket.

## Finding register

| ID | Severity | Dimension | Blocked state |
|---|---|---|---|
| KIT098-F01 | blocking | 2 rigor, 1 completeness | decision, execution |
| KIT098-F02 | blocking | 2 rigor, 1 completeness | decision |
| KIT098-F03 | should-fix | 1 completeness | implementation |
| KIT098-F04 | should-fix | 2 rigor | decision |
| KIT098-F05 | should-fix | 2 rigor | decision (claim scope) |
| KIT098-F06 | should-fix | 3 practice, 1 completeness | execution |
| KIT098-F07 | should-fix | 2 rigor, 3 practice | decision |
| KIT098-F08 | should-fix | 3 practice | execution |
| KIT098-F09 | should-fix | 1 completeness | implementation, execution |
| KIT098-F10 | minor | 3 practice | decision (provenance) |
| KIT098-F11 | minor | 1 completeness, 3 practice | execution |
| KIT098-F12 | minor | 4 references | none |
| KIT098-F13 | minor | 1 completeness, 3 practice | implementation |

### KIT098-F01 (blocking): the comparison with the notebook reference is undefined

**Evidence.** The acceptance criterion asks for "final states identical to the notebook reference" (`KIT-098.md:51`). The ticket names no comparator, no mapping and no settings binding:

- **Products.** KIT-040 freezes these as R objects in `kit-040/reference-<corpus>.rds` (`KIT-040.md:28-32`): filter mask, TMM factors, logCPM, grid scores, best k and percentage, W, the edgeR table, the RUVs-corrected logCPM, RUVIII `newY`, PCA and RLE inputs. The journey leaves parquet-backed states and products in a store.
- **Mapping.** The ticket does not say which product is compared with which state, component or product, or in what projection. `identical()` is strict about attributes: matrix versus data.frame, dimnames, factor levels.
- **Step list.** The plan's step list (`write-a-series-of-atomic-dragon.md:139`: import, design, filter, qc_finalize, normalisation, da, session_export, report) has no step that leaves the RUVIII export (`trx_ruviii_export`, KIT-095). Nothing in the journey would hold `newY` to compare.
- **Settings.** The notebook fixes group levels `c("RPMI","Sera")` (rmd l.124-140), takes its contrast from `contrast_strings.tab` (l.540) and fixes the RUV percentage range and the `max_k` cap (l.195-258). It uses `M` from `group` (l.704). The ticket does not bind the driver and recorder to these. The KIT-032 recorder model hard-codes "the second group minus the first" and `~ 0 + group` (`kit-032/oracle-record.R` design block).
- **The two journeys.** They are not specified (`KIT-098.md:42`). If one of them uses `M` from technical replicates (plan decision 2), the notebook has no matching product (KIT-093 records it "as a distinct product"), so "identical to the notebook" cannot be met for `newY`.

**Impact.** The only check in this ticket that is independent of the package cannot be decided as written. It could be satisfied by a comparator chosen after the outcomes are seen, which would be a post-hoc loosening.

**Smallest repair.** In the KIT-098 plan, before any run:
1. Name the comparator script (`kit-098/notebook-compare.R`).
2. Give a table mapping each frozen reference product to the journey step, state and component or product it reads from the **journey member's store**, with the projection applied (for example `as.matrix`, dimnames, level order).
3. Declare any tolerance and its reason, or state `identical` throughout.
4. Bind the design levels, contrast, RUV grid range, k cap, automatic selection and `M = group` for both journeys.
5. Name the two journeys. Either both use the notebook settings (a reproducibility repeat), or the second is stated to be compared with the in-process oracle only, with no notebook claim for `newY`.
6. Give `newY` a journey step or product, or remove it from the claim.

**Closure evidence.** A KIT-098 `plan.md` containing this mapping, committed before the first member's bindings commit. The comparator's output records every reference product as compared, or as excluded with a reason.

### KIT098-F02 (blocking): exactness and the lifecycle head rule can pass without checking anything

**Evidence.**
- **Missing oracles are skipped.** `kitOracleJourneyPopulations` (`workflow_kit.R:806-813`) returns `fingerprint = NULL` for a step with no bound oracle. `kitReproduceMismatches` then skips that step (`workflow_kit.R:737`, `if (is.null(expected)) next`). The verifier (`workflow_kit.R:1527-1535`) asserts only that the mismatch vector is empty. A step the recorder did not record therefore passes.
- **Members without oracles are treated as in-memory.** `kitSingleJourneyDecision` sets `in_memory <- is.null(case$oracles) && !identical(kitRouteOf(case), "DIA")` (`:464`). That makes `exact = NA`, skips the populations rule, and still returns `accepted`. `kitLifecycleDecision` uses the same test (`:619`) and waives `populations_complete` (`:630`) and `head_unchanged` (`:641`).
- **Lifecycle cases bind no oracles.** The omics lifecycle cases (`:251-268`) bind none. A trx lifecycle pair modelled on them would be classified in-memory, even though the transcriptomics module holds an on-disk owner from the start (KIT-096). The lipid lifecycle has had a store since KIT-070 (`kit-070/outcome.md`), so this misclassification already exists there.

**Impact.**
- "Both screen journeys exact" (`KIT-098.md:51`) can be accepted with steps never judged, or with exactness `NA`.
- "Lifecycle pair observed (completion and failure)" (`:50`) can be accepted without the head-unchanged rule. That rule is what the failure member exists to show.

**Smallest repair.**
1. Bind `oracles` on every KIT-098 case, the lifecycle pair included.
2. Make the trx path assert oracle coverage: every step whose `kitJourneySteps("transcriptomics")` state is non-NULL has a bound oracle, otherwise the result is `invalid` with `oracle_missing`.
3. For this route, decide in-memory status from whether a store exists, not from `is.null(oracles)`.
4. Add a harness test that fails before the change: a trx case with one oracle step missing, and a trx lifecycle case with a store.

**Closure evidence.** The harness test fails before the change and passes after. The stage JSONs for the lifecycle pair and the journeys show `in_memory = FALSE`, `exact = TRUE` and `head_unchanged = TRUE` on the failure member, and a per-step list of the oracles judged.

### KIT098-F03 (should-fix): the route-specific harness sites are not enumerated, and three would silently fall back

**Evidence.** The ticket lists, at a coarse level, steps, the driver, coordinator fields, runner and stage regexes, and tests (`KIT-098.md:26-29`). Found by grep:

**`workflow_kit.R`:**
- route allowlist `:15`;
- `kitMemberSeconds` `:111`;
- `kitCases` stage regex `:178`;
- step-map guards `:194-199`, which have no trx guard;
- step-map `switch` `:205-207`, **which defaults to the DIA step map for an unknown route**, against protocol v1.5 item 5;
- omics lifecycle block `:251-268`;
- `kitSelectorStage` `:482-486`;
- `kitLifecycleDecision` stage mapping `:622`, **which maps any non-lipid/metab route to `kit-lifecycle-complete`, a 600 s ceiling in the bindings, not 1,800 s**;
- `kitOperationIdentityProducts` `:729`;
- `kitGuiMemberDecision` `:957`, which uses DIA `kitStepIds()` (not used here).

**Driver:**
- dispatch `:196`;
- `runKitOmicsImportActions` and `runKitOmicsStep` map the route with `short <- if (identical(route, "lipidomics")) "lipid" else "metab"` (`:506`, `:538`), **so any other route is driven as metabolomics**.

**Coordinator:**
- stage regex `:18`;
- corpus list `:298-307`;
- normalisation merge `:313-316`;
- fault `:318`.

**Runner:** usage comment `:6-8` and stage regex `:85`.

**Bindings generator:**
- corpora (pattern at `:117-157`);
- `lifecycle_omics` (`:159-185`). It is a single spec shared by lipid and metab, with ITSD, cyclic loess and one build, which does not fit trx;
- harness `step_map_*` (`:996-1011`), which has no `step_map_trx` entry for KIT-095's `kit-001/step-map-trx.json`;
- oracle ingestion (pattern at `:393-412`);
- the ticket keys.

**Unchanged:** `workflow_kit_member.R` and `kit-053/stage-summary.R` are route-agnostic.

**Impact.**
- A trx case could carry the DIA step-map hash in its binding digest.
- A trx lifecycle member would be judged against the wrong ceiling and become `invalid` on `elapsed_within_ceiling`.
- A wrong import or tab path would fail at run time and use up attempts.

**Smallest repair.**
1. List every site above in the KIT-098 plan.
2. Replace the three silent defaults with explicit refusals for unknown routes.
3. Give trx its own lifecycle spec (or a per-route `lifecycle_omics` entry), its fault (`after_step` and `failed_step`, preferably a worker step) and its step-map binding.
4. Cover each site with a harness test that fails before the route exists.

**Closure evidence.**
- Harness tests: unknown route refused, trx step map bound into the digest, trx lifecycle ceiling 1,800 s, trx driver prefixes.
- A grep in `dependents.txt` that shows no remaining `lipid|metab` binary.

### KIT098-F04 (should-fix): intermediate transcriptomics states would go unjudged

**Evidence.** The verifier compares the state that is current after each journey step (`workflow_kit.R:1513-1530`; the driver attaches one state per step, `workflow_kit_driver.R:322-337`). The plan's step list folds TMM and logCPM, the RUV grid and RUVs into one `normalisation` step (`write-a-series-of-atomic-dragon.md:139`). KIT-095 declares them as separate worker steps (`KIT-095.md`, step table). Their states are fingerprinted in the store but never compared with the oracle.

**Impact.** Grid scores, best k and percentage, and W, the parts of the notebook that are most sensitive scientifically, are judged only through the final state. If the DA and correction outputs agree, a compensating difference upstream is unlikely but not excluded.

**Smallest repair.** Either give each worker step its own journey step (the driver marks one step per tab action), or have the trx verifier compare every history state against the oracle's history. F01's comparator also reads these intermediates from the journey store.

**Closure evidence.** A receipt that lists every trx history state as judged.

### KIT098-F05 (should-fix): "control" is the wrong term, and the claim limits are not stated

**Evidence.**
- The ticket calls the same-build in-process oracle "the control" (`KIT-098.md:32`). In the protocol a control is a build run in the same session for a memory ratio (§2, §3). The in-process oracle is an exactness oracle, and §7.1 already allows it for a route without KIT-005 oracles: it "proves that the owner path reproduces the in-process path of the same kernels on the same build; it is not a cross-build claim" (`protocol.md:53`).
- The recorder model drives the same functions with hand-coded defaults (`kit-032/oracle-record.R`).

**What each comparison proves:**
- **Journey against the same-build oracle.** The GUI, owner, worker and store path is equivalent to the in-process path for the same inputs and settings. It proves nothing about scientific correctness. A defect shared by both (a kernel, a declaration, or a setting both copy) passes.
- **Journey store against the notebook reference.** The package reproduces independently written code (the notebook calls edgeR, RUVSeq and ruv directly and loads no MultiScholaR) on one synthetic corpus, for the frozen products. This is the only check that is not circular. It is transitive with KIT-093 and KIT-095, but it is only new evidence if it reads the journey member's store. Reading the oracle's store would repeat KIT-095.
- **Lifecycle pair.** The instrumentation observes completion and failure.

**Impact.** Wording that could be read as a memory or correctness control overstates the claim. v1.7 may also re-authorise what §7.1 already allows.

**Smallest repair.** Word v1.7 as:
- the oracle under §7.1;
- the notebook reference as a second, implementation-independent exactness check, compared from the journey store;
- no memory claim (no control build exists);
- no cross-build or real-data claim;
- the trx stages and selectors.

**Closure evidence.** The v1.7 text, with the user's decision quoted, and the outcome's claim-scope line.

### KIT098-F06 (should-fix): the protocol version and its hash conflict with KIT-092 and with existing bindings

**Evidence.**
- KIT-092 is a user-decided change to "the frozen protocol (v1.6)", with "the protocol version and the change ... recorded in the plan" (`KIT-092.md`). It runs before KIT-098 (via KIT-040) and may itself become v1.7.
- Every case binding digest includes `protocol = bindings$protocol_sha256` (`workflow_kit.R:204`, `:237`). The verifier asserts `identical(operation$binding, case$binding)` (`:1432`).

**Impact.**
- The version number can collide.
- Editing `protocol.md` changes the binding of every existing route. Older stages can then no longer be re-verified under the new bindings, so the evidence for "Other routes unchanged" (`KIT-098.md:59`) has to be defined differently.

**Smallest repair.**
1. Take the version number after KIT-092 closes.
2. Amend before any KIT-098 member is bound.
3. State that earlier receipts stand as recorded.
4. Evidence backward compatibility with a `bindings.json` diff limited to the trx keys and the protocol hash, plus the harness tests for the other routes.

**Closure evidence.** The protocol header shows the version sequence after KIT-092. The bindings diff is in `kit-098/`.

### KIT098-F07 (should-fix): the listed worker steps are checked by nothing

**Evidence.**
- v1.7 "lists the worker steps" (`KIT-098.md:33`). Only `kitWorkerPairDecision` reads `worker_steps` (`workflow_kit.R:896-929`), and that is for pairs. `kitSingleJourneyDecision` (`:450-478`) ignores them.
- Since KIT-089 the in-process path runs the registry default. A journey that silently ran a step in process would therefore still be exact against the oracle.
- `kitGuiMemberDecision` already shows how to check dispatch from worker journals (`:967-985`).

**Impact.** The plan's architectural claim ("compute in the worker") is not verified by this ticket's evidence.

**Smallest repair.** Either verify a completed worker journal for each listed step in the trx journeys' receipts, reusing the `kitGuiMemberDecision` journal scan, or drop the list from v1.7 and state that worker placement is not judged here.

**Closure evidence.** The receipt's `dispatched` counts per worker step, or the v1.7 text that excludes the check.

### KIT098-F08 (should-fix): the experiment budget leaves out runs, member kinds and terminal outcomes

**Evidence.** The budget is "The lifecycle pair, then two screen journeys (24 by 6,000), one attempt each, 1,800 s; one bounded correction pass" (`KIT-098.md:42`). It is silent on:
- the in-process recorder runs (how many, on which settings, ceiling);
- the notebook comparison;
- whether the members are counted (plan run ledger row: "lifecycle, two journeys", `00-KIT-PLAN.md:237`);
- the member kinds and selectors (for example `single_journey` with `whole_journey = TRUE`, `kit:journey-098[a-z]:candidate:1`, `kit:lifecycle:trx:(complete|fail)`);
- the terminal outcomes. A valid `rejected` is recorded and never retried. The correction pass is only for a demonstrated instrumentation or harness failure (protocol §1).
- what happens if the user declines v1.7.

**Impact.** AGENTS.md requires purpose, bindings, ceiling, attempts and decision to be declared before execution. Without them, a mismatch could be "corrected" by rerunning or by changing the recorder.

**Smallest repair.** Put these items in the KIT-098 plan and the ledger row before the bindings commit.

**Closure evidence.** The plan and ledger row are committed before the first launch. The outcome reports attempts against the budget.

### KIT098-F09 (should-fix): binding the corpus and reference inputs is left implicit

**Evidence.**
- KIT-040 adds the `trx` kind to `kit_omics_corpus.R` but does not modify `bindings-generator.R` (`KIT-040.md` `files_to_modify`). So binding the trx corpora falls to KIT-098 and is not mentioned.
- `bind_corpus` for omics expects `assay_files` plus `design.tsv` (`bindings-generator.R:148-157`). The trx files are `counts`, `design_matrix.tab`, `contrast_strings.tab` and `truth.json`.
- The coordinator builds `corpus` from `report_file`, `fasta`, `design`, `assays` and `vendor` (`workflow_kit_coordinator.R:298-307`).
- The coordinator merges the shared `bindings$journey$normalisation` (the proteomics cyclic loess and RUV 1-20 % settings) into every route (`:312-316`).

**Impact.** The journey could run on a corpus whose bytes differ from the one the notebook reference used, or with normalisation settings that do not belong to trx, and the mismatch would surface only after a run.

**Smallest repair.**
1. Bind the trx screen corpus with sha256 values equal to the inputs recorded in `kit-040/reference-<corpus>.rds`, and assert that equality in the generator, as the oracle loop does with `corpus_sha256` (`:405`).
2. Bind trx's own normalisation and DA settings on the ticket keys, not the shared journey block.

**Closure evidence.** The generator assertion, and the corpus and settings entries in `bindings.json`.

### KIT098-F10 (minor): the recorder's provenance is incomplete

**Evidence.**
- The model recorder records the build, corpus sha and settings. It does not record the `workflow_kit.R` hash it fingerprinted with, the BLAS or LAPACK library, or the variant run per step (`kit-032/oracle-record.R`, the `recorded <- list(...)` block).
- The oracle names its products `<type>:<sha12>` from the canonical form (same file). A canonical-form change after recording (KIT-092 is one) silently changes those names.
- The bindings and verifier label every in-process oracle "protocol v1.5" (`bindings-generator.R:410`, `workflow_kit.R:1539`).
- Host check: R 4.6.1 with reference BLAS `libblas.so.3.12.0`. `identical` depends on this staying fixed.

**Impact.** If the oracle and the verifier disagree after a harness edit, the record cannot show why.

**Smallest repair.** In the trx recorder:
- record the `workflow_kit.R` sha256, the BLAS and LAPACK paths and the per-step variants;
- read the settings from the same binding the driver uses rather than hard-coding them;
- label the protocol in force.

**Closure evidence.** These fields appear in `kit-098/oracles-*.rds` and in its bindings entry.

### KIT098-F11 (minor): the README's pre-launch chain check is not referenced

**Evidence.**
- The README's "Before a counted member exercises a new operation role or class" check (items 1-6: registry, compute helper, worker dispatch, vehicle checks, publication, verifier branch) is not cited by the ticket.
- Transcriptomics introduces a new class (`TranscriptomicsAssayData`) and new compute roles.

**Impact.** KIT-010's three launches are the precedent for the cost of skipping it. The 1,800 s single attempts are at risk.

**Smallest repair.** Run the chain in process on the screen corpus before the lifecycle pair, reusing KIT-095's tests where they already cover an item, and record which items they cover.

**Closure evidence.** A chain-check record in `kit-098/`.

### KIT098-F12 (minor): small reference inaccuracies

**Evidence.**
- "793 lines" holds by line count. The file has no trailing newline, so `wc -l` reports 792. sha256 44d5f643….
- `bindings.json`, which is generated and committed, is not in `files_to_modify`. That is consistent with the series convention (KIT-089 and KIT-092 omit it too).
- The plan's step list puts `da` before `session_export`. On lipid and metab the DA tab loads the export (`workflow_kit.R:36-38`). The order has to be confirmed against KIT-097's GUI.
- **Out of scope, recorded only:** the plan table row for KIT-088 (`00-KIT-PLAN.md:240`) lists dependencies "KIT-040, KIT-087". The ticket and the index say KIT-087, 089, 090, 092 and 099.

**Impact.** Cosmetic, except the step order, which F01 and F03 depend on.

**Smallest repair.** Settle these in the KIT-098 plan.

**Closure evidence.** The plan states the step order with GUI evidence.

### KIT098-F13 (minor): harness ownership overlaps with KIT-094

**Evidence.**
- KIT-094's description moves "the harness" lists to binding fields. Its files include `tools/profiling/workflow_kit.R`, and the plan names `kitOperationIdentityProducts` (`workflow_kit.R:729`) among them.
- KIT-098 adds route alternations to the same file.
- Whether trx `*_science` products are judged or ignored is undecided. KIT-031 ignored the lipid and metab ones because of an `@args` timestamp in their bytes.

**Impact.** Two tickets editing the same mechanism. Adding a trx alternation to the ignore regex without a stated cause would hide products from judgement.

**Smallest repair.**
1. Consume KIT-094's binding fields.
2. Judge trx products by default.
3. Add a trx pattern to the identity-product set only with the cause shown, in the KIT-031 way.

**Closure evidence.** The plan's decision, and a harness test that a trx product is judged.

## Verdicts

1. **Completeness: not sufficient as written.** The coarse outcome list misses several required sites:
   - the lifecycle fault and spec;
   - the step-map binding;
   - the corpus and settings binding;
   - the selector and ceiling mappings;
   - the identity-product decision;
   - three silent fallbacks (step map to DIA, unknown omic to metab, lifecycle ceiling to 600 s).

   The recorder and protocol amendment are named but their contents are not. The notebook comparison has no comparator or mapping (F01, F03, F09, F11, F13). The planning gate is the right place to fix this, and the plan must cover these items before any run.
2. **Scientific rigor: sound in principle, not decidable as written.** The same-build oracle is an equivalence check for the GUI and worker path and is circular for science. The notebook comparison is the only independent check, and it is undefined (F01). Exactness can pass without judging anything through existing harness defaults (F02). Intermediate RUV states would go unjudged (F04). "Control" overstates the claim (F05). Worker placement is asserted but not checked (F07). KIT-092 removes the partitioned-variant byte issue only if its canonical form also covers what the trx states and products hold. KIT-098 should pin the harness hash it records with (F10). Partitioned against whole remains KIT-099's.
3. **Implementation practice: mostly aligned, with gaps.**
   - Aligned: test-first, the lifecycle pair before the journeys, one attempt each, one correction pass.
   - Missing: the recorder and comparator runs, the counted or uncounted status, selectors and terminal outcomes (F08).
   - The protocol version and hash interplay with KIT-092 and with the other routes needs sequencing (F06).
   - Reusing KIT-094's binding fields would keep the harness change smallest (F13).
4. **Reference accuracy: accurate.** Everything cited by the ticket and the brief exists as described:
   - `kitJourneySteps` `:12`, `kitOracleJourneyPopulations` `:806`, `kitCanonical` `:1061`, `kitOperationIdentityProducts` `:729`;
   - the oracle branch of `verifyKitMember` around `:1495-1543`;
   - `runKitOmicsStep` `:536` and `runKitOmicsImportActions` `:505`, which hard-code lipid and metab;
   - protocol §7.1;
   - the R CMD check baseline of 3 WARNINGs and 3 NOTEs (`kit-089/check.log`);
   - the KIT-089 byte-hash finding.

   Only the minor nits in F12 remain.
