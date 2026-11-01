# KIT-073 outcome: one shared omics declaration core

Status 2026-11-01: **completed**. Stage 3b, ticket 2 of 10.

- **Commits** on `performance-v2`: code 52a2c179, bindings and records 54a5d53a, and this evidence commit.
- **Base:** 67a1a3d9 (KIT-072 closed).
- **Install:** `kit_073` = `core/browser-smoke.R4tSbp` (52a2c179), over the kit_072 library.
- **Ticket:** `tickets/kit/KIT-073.md`, revision 1.1.

## What changed

1. **The core: `R/utils_workflow_omics_declarations.R` (new, 454 lines).** One constructor per step builds the twelve declarations of an omic from its binding. The same file holds:
   - `omicsWorkflowOperations(binding)` and `omicsDeclarationBindings()`;
   - the wrappers `runWorkflowOmicsReference*` (intensity, duplicates, ITSD, log2, between-sample, RUV, correlation, DA);
   - `recordWorkflowOmicsDesignProvenance` and `replayWorkflowOmicsStep`;
   - the omic-neutral functions moved out of the lipid file: `runWorkflowOmicsReferenceOperation`, `markWorkflowOmicsStoreStep` and `commitWorkflowOmicsStoreStep`.
2. **The bindings.** `lipidomicsDeclarationBinding` and `metabolomicsDeclarationBinding` hold the names and strings as data: class, slots, prefix, feature noun, long table, design input prefix, abort scopes, intensity parameter and threshold names, and step labels. Where the two omics call their kernels differently, each binding supplies its own adaptor:
   - `filter_intensity`: lipid's local `…Fn` alias, metab's bare name;
   - `resolve_duplicates`: lipid's per-assay loop with counts, metab's `resolveMetabDuplicateAssayData` with stats;
   - `optimise_ruv`, `extract_best_k`, `extract_controls`;
   - `filter_correlation`;
   - `run_da`.

   The omic files fell from 471 + 389 to 164 + 147 lines. They keep every name the modules, the record, the harness and the tests use: the `run…Reference*` wrappers, the `…Operation()` constructors, `…WorkflowOperations()`, `…StepLabels`, `…FeatureCounts`, the manifest predicates, the design recorders, the replay steps and the lipid aliases of the wrapper core.
3. **The registry and the record.** The registry builds both omics in one loop over the bindings. The record replays steps through one omics branch keyed by `omic_type`, and re-runs DA through one branch keyed by the DA operation. DESCRIPTION Collate puts the core before the omic files.

## Results

| Gate | Evidence | Result |
|---|---|---|
| **Nothing observable changed** | [freeze-dump.R](freeze-dump.R), `tests/testthat/test-workflow-omics-declarations.R`, `tests/testthat/fixtures/kit073-omics-declarations-frozen.rds` | The dump was frozen on the base commit. It holds the 24 registry entries and every declaration field except the kernel, plus both omics' wrapper journeys on the declaration fixtures: every state, the kernel values, the correlation tables and the DA manifest and results, with times scrubbed. Two base runs were identical, so the dump is deterministic. The current tree's dump is **identical**: 40 expectations, 0 failures. |
| Affected tests | [affected-tests.txt](affected-tests.txt) | 148 files, one process per file: **12,878 pass, 0 fail, 0 error**. |
| R CMD check | [check-item1.log](check-item1.log) | Baseline: 3 WARNINGs, 3 NOTEs, with the code note identical to KIT-072's. |
| Exactness set (1)(3)(6) | [exactness-gates.txt](exactness-gates.txt) | The `lipid73` and `metab73` records are valid with no legacy step, and **identical** to `lipid31` and `metab32`. The KIT-028 replay is **exact** (16 states, 265 s). The LFQ, TMT and DIA_LIMPA records are **identical** to the kit_065 oracles. |
| Exactness set (2) | [073l](kit073-kit-journey-073l-stage.json), [073m](kit073-kit-journey-073m-stage.json) | Both accepted and exact. `073l`: 241 s, 2.52 GiB peak. `073m`: 239 s, 2.62 GiB peak. No sampler gaps. |
| Exactness set (4) | [073l](member-record-replay-073l.json), [073m](member-record-replay-073m.json) | Both records are valid with no legacy step. Both replays are exact (10 states), and a tampered input is refused. |
| Exactness set (5) | [lipid](kit073-kit-lifecycle-lipid-stage.json), [metab](kit073-kit-lifecycle-metab-stage.json) | Both observed, launched one at a time. Lipid: 240 s, 2.68 GiB peak. Metab: 240 s, 2.59 GiB peak. |

**Memory (reported, not claimed).** The lipid lifecycle peak is 0.15 GiB above kit_032's 2.53 GiB. The lipid journey on the same build and corpus peaked at 2.52 GiB, and the frozen dump shows that no step's computation changed. Earlier lipid lifecycle peaks were 2.46 GiB (kit_070), 2.53 GiB (kit_032) and 2.53 GiB (kit_072), so the reading is run-to-run variance of the sampled peak. It is recorded so later tickets can compare against it.

## Backward compatibility and rollback

- Every name resolves; the frozen-dump test checks the identity of the wrapper aliases.
- Records and stores of kit_031, kit_032 and kit_072 replay unchanged.
- No saved format changed. Proteomics is identical.
- Rollback: `git revert` the ticket's commits.

## Handoff

KIT-074 (neutral names) has one declaration file to move. KIT-040 (transcriptomics) starts from the core. The exactness set's references stay `lipid31`/`metab32`, which equal `lipid72`/`metab72` and `lipid73`/`metab73`.
