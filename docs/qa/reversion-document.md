# DIA saved-entity reuse and memory proof: provenance and implementation basis

- **Document version:** 1.0, 2026-09-11.
- **Status:** planning basis; no new implementation or performance acceptance.
- **Initial scope:** reusable package utilities and standalone DIA fixtures for
  history, scientific products, DA, plots, table lookup and exports; then DIA/IQ
  GUI integration and real Limpa verification.
- **Restored baseline:** local `janitor` at
  `c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017`.
- **First deliverable:** tested save, hydrate, lookup-by-factors, projected-read
  and cache-release functions in a flat, golem-compliant package R file.
- **Primary question:** where can saved entities reduce resident memory,
  reconstruction, repeated science and rendering across the DIA workflow?
- **Authority:** the user's checkpoint discussion and fixture-first request.
  This document reorganizes that discussion; it does not authorize branch
  deletion, erase earlier evidence, or certify a storage technology.

**Navigation:** [decision and sequence](#1-read-this-first),
[provenance](#2-provenance-observations-and-unresolved-hypotheses),
[API contracts](#5-reusable-package-foundation-before-gui-work),
[consumer comparisons](#56-consumer-lookup-and-reuse-map),
[fixture coverage](#63-required-scenario-coverage),
[measurement and gates](#8-measurement-protocol-and-decision-gates),
[implementation checklist](#10-implementation-handoff-checklist).

## 1. Read this first

Build and measure the reusable primitives **outside Shiny before changing GUI
consumers**. The standalone fixtures must call the same package functions that
the GUI will eventually use. A benchmark-only storage path is not acceptable.

Keep the current S4 working state in memory. Store older checkpoints on physical
disk and hydrate only the requested state on revert. Keep small metadata and
bounded derived products available without reconstructing an assay.

The user's follow-up explicitly includes finding saved scientific entities for
DA and plot generation. Checkpoint eviction is one use case of the reusable
foundation, not the entire objective. Compare complete S4/RDS restoration with
loading an exact saved fit, selected result table, scientific plot product or
ready-rendered asset. Compare each with the best available resident-data path
as well; saved lookup is not automatically preferable to memory access.

A saved fit or scientific product may itself be an RDS file. The question is
whether loading that smaller, sufficient entity avoids reconstructing a complete
workflow S4 or recomputing its results. Entity granularity and storage format are
separate choices and must be measured separately.

Start with whole-object RDS as the simplest compatible storage comparator.
Parquet/Arrow/DuckDB are candidates for workloads that benefit from selective
table reads; they are not a prerequisite for saving or finding a checkpoint.
Compare alternative formats only after the basic lifecycle works.

The intended order is:

1. Preserve provenance and establish a clean baseline from local janitor.
2. Specify and implement the small reusable saved-entity API and checkpoint use.
3. Prove exactness, identity, ownership and failure behavior with small fixtures.
4. Demonstrate the declared memory and performance benefit with large,
   GUI-compatible synthetic DIA fixtures outside the GUI.
5. Replay actual DIA science and consumers outside the GUI, including DA, plots,
   table browsing and exports; test each reuse mechanism separately.
6. Integrate only accepted operations into one complete DIA/IQ GUI journey.
7. Verify real Limpa before claiming DIA coverage.
8. Extend to other omics only after their own compatibility and benefit checks.

### 1.1 DIA-only: out-of-the-box baseline, then targeted refactors

The first execution establishes how the selected janitor DIA workflow behaves
**out of the box**, with its existing scientific methods, ordinary execution
path, checkpoints and required outputs. Record its elapsed time, CPU, retained
memory, combined peak and scientific/output results before performance refactors.
Begin with DIA/IQ; real Limpa is the second DIA path required for a DIA-wide claim.

Next, prove the shared save/hydrate/lookup utilities and their DIA consumers on
large compatible fixtures outside Shiny. Use the unmodified scientific functions
as the reference. Only after that evidence passes should targeted GUI refactors
replace a named consumer's access or reuse operation. Each refactor must identify
the measured cost it addresses and retain a directly comparable baseline.

"Ecosystem" here means the DIA workflow's connected operations: import/QC,
checkpoint history, normalization/RUV, DA, plotting, lookup and required exports.
It does not authorize concurrent metabolomics, lipidomics, LFQ or TMT migrations.
Design the utilities for reuse, but implement and accept their consumers in DIA
first. Improvements in other omics require their own later verification.

"All possible scenarios" means an explicit, versioned coverage matrix of supported
cases, boundaries and failures. Finite tests cannot prove every possible dataset,
host or interleaving. Unsupported cases and unexecuted cases must be visible.
Correctness must hold across supported cases; performance claims must state
exactly which workloads demonstrated the claimed benefit.

## 2. Provenance, observations and unresolved hypotheses

The original user-created file contained two pasted discussion responses, with
terminal-rendered tables and truncated paragraphs. This version restores the
structure and missing context, retains their substantive observations, and adds
the user's requirement to build reusable utilities and large fixtures first.
The subsequent request broadens the benefit proof to saved-entity lookup for
DA, plot generation and the wider workflow ecosystem, still using DIA first.

SHA-256 of that original 376-line file before this edit:
`6c6e4fe88cb6566bf9e26dc621ab0cb8e0ebfcf5fb73729fa4712496f5ba185b`.

### 2.1 What was observed

| Observation | Implication | Evidence and limit |
| --- | --- | --- |
| Memory-mode reversion shortens active history but leaves later records in `states`. | Inactive objects may remain resident. Inspect all owners before expecting eviction to free memory. | `WorkflowState$revertToState()` in [the current implementation](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/R/utils_workflow_state.R); also observed in local janitor. |
| S4 objects can share unchanged allocations. | Summing separate `object.size()` values does not establish physical duplication. | [Original memory audit](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/docs/qa/dia-memory-audit.md). Measure unique retained allocations and process memory. |
| All seven policy-09 GUI routes missed their fixed elapsed and combined-memory gates. | The full Parquet/DuckDB migration has not demonstrated acceptance. | [Policy-09 results](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/docs/qa/omics-gui-repair-121-results.md). Controls are corrected memory builds, not clean janitor. |
| Later pilot memory decreased while elapsed regressed. | A memory reduction alone is not an accepted end-to-end performance improvement. | [Policy-10 journal](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/docs/qa/omics-gui-integration-117-121.md). Aggregate observations do not identify the precise cause. |
| FragPipe profiling found substantial composite/raster/RDS work. | Plot and export reuse are separate, evidence-backed optimization candidates. | [Profile record](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/docs/qa/omics-gui-repair-116-121.md#116-profile-and-comparator-freeze): 16.95 of 23.69 sampled normalization CPU seconds attributed to composite extras. Inclusive costs overlap; do not sum them. |
| Protein and metabolomics export serialize the same session twice for timestamped and latest files. | Compare one serialization plus a required file copy against two serializations. | [Protein export](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/R/mod_prot_norm_session_helpers.R), [metabolomics export](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/R/mod_metab_norm_session_helpers.R); also observed in janitor. |
| Full protein legacy export materializes historical states together. | A cold-history design may reduce ordinary residency but still peak during export. | [History materialization](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/R/mod_prot_norm_session_artifact_helpers.R). Retain this boundary in acceptance measurements. |
| This host's `/tmp` was RAM-backed `tmpfs`. | Writing there does not demonstrate physical-disk memory relief. | Filesystem observation on 2026-09-11; recheck on every measurement host. |
| `checkMemoryBoth()` is already disabled; the protein composite builder still explicitly calls `gc()`. | Do not claim savings from removing already-disabled diagnostics. | [Runtime helpers](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/R/func_general_runtime_helpers.R), [composite builder](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/source/R/mod_prot_norm_qc_support_helpers.R). |

### 2.2 Failure lessons to carry forward

The earlier branch also exposed correctness defects: product identities omitted
component selections/patches, a batched product path lost unavailable-panel
handling, and caught-error publication rollback was mistaken for worker-crash
recovery. Its export cache did not cover every component decoder.

Carry those failure cases forward as tests. Do not transplant the corresponding
implementations unchanged. A successful GUI assertion or a prior build's
scientific comparison does not certify a newly changed codec, key or consumer.

The prospective checkpoint proof remains **unproven**. It is not a claim that
RDS, Parquet, DuckDB or a background worker will be faster on this workload.

### 2.3 Preserve the restart boundary

Before any branch transition, preserve this document, the dirty worktree,
untracked implementation files, selected scientific fixtures, benchmark inputs,
receipts and relevant external evidence. Verify the archive can be read.

A Git branch or tag alone does not preserve uncommitted files or evidence under
`/tmp`. Local janitor was 32 commits ahead of `origin/janitor`;
resolve the intended local commit explicitly. Prefer a separate clean worktree.

The restart plan must review scientific fixes made after janitor individually.
A pure-janitor baseline and a baseline with essential correctness repairs must
have different source identities and labels.

Archive-relative evidence locations must be recorded before leaving this branch.
Several links in this document point to files added on performance and will not
exist automatically in a clean janitor checkout.

### 2.4 Janitor restoration completed on 2026-09-11

The working checkout now points to local janitor at
`c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017`. Tracked runtime, package and
test files match that commit. This provenance document is carried forward as an
untracked planning addition; no new scientific or performance implementation
has started.

Recovery archive:
`/home/doktersmol/Projects/APAF/MultiScholaR-performance-archive-20260911-5ywjc71f`.

- [Recovery receipt](../../../MultiScholaR-performance-archive-20260911-5ywjc71f/receipt.json) records the source identity,
  verified file count and every original-to-archive path mapping.
- `source/` preserves 4,057 verified tracked/untracked source and evidence files.
- `source.bundle` preserves complete Git history for performance, janitor and
  the recovery stash. Bundle verification passed.
- Recovery stash: `8b45761f6d70a4a50d19d47ea517015d69c6453a`.
- `generated-in-repository/` holds approximately 133 GiB of former workload,
  comparator, browser-run and build outputs.
- `temporary-runs/` holds 145 former `/tmp/multischolar-*` directories,
  approximately 8.4 GiB.
- 101 already-empty experiment directories were removed during Git stash cleanup.
  Populated artifacts were relocated, not permanently deleted.
- Manuscript files, local R-library files, repository conventions and editor
  state were preserved in place. Local Git exclusions preserve their prior
  ignored status without changing janitor's tracked ignore file.

Historical reports retain their original embedded paths and hashes. Resolve moved
paths through the recovery receipt; do not rewrite archived evidence. Local
evidence links in this working document now point to the archived source.

## 3. Scope and fixed compatibility contracts

### In scope first

- Saving complete DIA checkpoint records and hydrating one selected state.
- Finding saved records by typed metadata factors without reading payloads.
- Saving/retrieving exact scientific products, fitted models, result tables and
  rendered assets for the existing DIA consumers.
- Comparing metadata-only lookup, whole-object hydration, projected table reads
  and direct saved-product/asset use, including the cost of producing the entity.
- Bounded resident caching, safe release, reopening and failure recovery.
- Real `PeptideQuantitativeData` and `ProteinQuantitativeData` fixtures.
- Large synthetic data, small edge cases, and real IQ/Limpa scientific adapters.
- Exact history, reversion, session-export and automatic-output compatibility.
- Independent measurement of memory, elapsed time, CPU and I/O.

### Held until separately justified

- Broad migration of the remaining omics workflows.
- General delta stores, arbitrary query engines, plugin systems or new packages.
- Replacing statistical methods, changing populations, reducing precision or
  dropping required output files.
- Mandatory persistence of every plot product or every small parameter leaf.
- HDF5/delayed-array conversion unless a required active matrix cannot fit.
- Rebuilding enrichment or reports beyond existing required file consumers.

### Non-negotiable behavior

1. Preserve one public R package and the existing S4 classes, public methods and
   dispatch. Apply all four [repository conventions](../../conventions/README.md).
   Their existing-code compatibility rule takes precedence over an unrelated
   S4 architecture redesign.
2. A normal scientific `getState()` call returns a real compatible S4 object.
   Do not substitute paths into S4 slots that require matrices or data frames.
3. Direct use of public `WorkflowState$states`, serialized histories and
   compatibility helpers must be inventoried. An internal descriptor is not
   equivalent to the public historical record.
4. Preserve exact IDs, order, row names, attributes, factor levels, annotations,
   missing/nonfinite behavior, configuration and substantive provenance.
5. Keep global normalization, PCA, density's distinct PCA population, RUV,
   contrasts, empirical-Bayes and FDR populations unchanged. Preserve
   `runTestsContrasts(weights = NA)`.
6. Keep checkpoint names/order and revert/reset behavior. Cold storage must not
   require rerunning science to recover a checkpoint.
7. Preserve complete self-contained legacy session RDS files and required
   timestamped/latest files, TSVs and figures. Changing the public export format
   is outside this proof.
8. Keep the existing 8 GiB / zero-swap / 900-second service limits and registry
   limits. Optional resident caches remain at or below their existing 64 MiB
   ceilings. Active scientific inputs and serialization working memory must
   also fit within the total process/service budget.
9. Scientific failure, checkpoint-write failure, output-publication failure and
   cancellation are distinct outcomes. Output retry creates no new scientific
   checkpoint and performs no model fitting.

## 4. Working-state and checkpoint ownership

| Category | Intended owner and residency | Release rule |
| --- | --- | --- |
| Current S4 | One scientific execution owner; resident while active | Replace only after the next state is accepted |
| Current operation's inputs/output | Same owner; scoped to the action | Release obsolete inputs before later serialization/rendering where safe |
| Older checkpoint payloads | Immutable files in project-owned disk storage | Evict resident copies only after successful persistence |
| History index and saved factors | Small, bounded metadata | Keep while the corresponding checkpoint is retained |
| Recently hydrated history | Optional byte-bounded cache | Evict by declared policy; persistent checkpoint remains |
| Scientific plot products | Small scores/statistics keyed by exact dependencies | Reuse within budget; persist only when measured reuse justifies it |
| Rendered assets | Files keyed separately by product, style and dimensions | Invalidate only dependent assets |
| Legacy export objects | Explicit export scope | Release after publication; include this scope in peak measurements |

### Save and forward execution

1. Validate the object and checkpoint request before changing active state.
2. Write an immutable checkpoint record to a request-owned staging location.
3. Close the file and complete the declared integrity checks.
4. Publish its record/reference; only then make it eligible for resident eviction.
5. Advance the current state through the coordinator's existing acceptance rule.
6. Release older resident payloads and obsolete action inputs.

If persistence fails, preserve the last valid checkpoint and recoverable resident
object. Never report a safely persisted checkpoint when only a temporary file
exists.

### Revert and reopening

1. Resolve an exact checkpoint reference from metadata.
2. Hydrate and validate the requested object.
3. Replace the active state only after successful hydration.
4. Update active history according to existing semantics.
5. Evict obsolete resident objects within the declared byte budget.

A failed hydrate leaves the previously active state usable. Reopening must work
without the previous R6 instance or its resident cache. Unreferenced staging files
must not become the latest valid state.

### Ownership pitfalls

- Clearing the R6 reference is insufficient if module state, a closure, plot,
  promise, global binding or diagnostic collector still retains the object.
- Metadata-only states may share allocations. Test shared and independently
  transformed histories separately.
- By-reference mutation must not alter an older checkpoint or a shared current
  input. Persist/copy at explicit ownership boundaries.
- Background writers cannot accumulate an unbounded queue of full S4 objects.
  Any queued payload and transport copy counts against memory and elapsed cost.
- Do not add a second scientific owner just to serialize a snapshot. The first
  utility proof can run in one R process. A later reusable-worker integration
  must measure startup and transport and return bounded GUI results.

## 5. Reusable package foundation before GUI work

### 5.1 Initial file surface

The names below are a proposed implementation contract, not existing functions.

| Planned file | Responsibility |
| --- | --- |
| `R/utils_workflow_saved_entities.R` | Shared store handle, save/hydrate, typed metadata lookup, projected-read capabilities, descriptor validation and bounded cache lifecycle; checkpoint wrappers initially live here |
| `tests/testthat/helper-workflow-checkpoint-fixtures.R` | Deterministic GUI-compatible fixture recipes and small builders |
| `tests/testthat/test-workflow-checkpoints.R` | Exactness, lookup, ownership, round-trip and failure tests |
| `tests/testthat/test-workflow-checkpoint-science.R` | Original-method versus hydrated-method DIA comparisons |
| `tests/testthat/test-workflow-saved-consumers.R` | Saved-product/fit/table/asset reuse versus original DIA consumers, including hit/miss/invalidation cases |
| `tests/testdata/workflow-checkpoints/` | Small scenario specifications, seeds, expected metadata and fixture receipts |
| `tools/profiling/run_workflow_checkpoint_proof.R` | Explicit large standalone runner using the installed package and existing measurement tools |

Keep `R/` flat. Start with one cohesive utility file; split by a concrete
responsibility only if needed. This is package business logic, not a new Shiny
module, and does not need a UI/server pair.

Functions use camelCase; arguments/list fields use snake_case. Document internal
functions with roxygen `@noRd`. Qualify dependencies, maintain required package
declarations, and regenerate NAMESPACE/collation from their source annotations.
Preserve native registration. Update the dependency lock only when dependencies
actually change. No runtime installation, global store singleton, `setwd()`,
or package-wide reformatting.

### 5.2 Minimal API and responsibilities

| Proposed entry point | Contract |
| --- | --- |
| `newWorkflowEntityStore(project_root, identity, policy)` | Create/open one explicitly owned store; policy supplies limits, supported entity kinds, codec and integrity rules |
| `saveWorkflowEntity(store, value, specification)` | Save one checkpoint, scientific product, fit, result table or owned rendered asset using its declared contract; return an accepted immutable reference |
| `describeSavedWorkflowEntity(store, entity_ref)` | Return metadata, saved factors and declared shape/class without payload hydration |
| `findSavedWorkflowEntities(store, factors, limit)` | Return bounded matching descriptors using metadata only; no payload reads or implicit computation |
| `hydrateSavedWorkflowEntity(store, entity_ref)` | Validate exact ownership/integrity and reconstruct one compatible object; account for its resident lifetime |
| `readSavedWorkflowTable(store, entity_ref, selection)` | Select declared rows/columns before collection where the codec supports it; reject unsupported projection or oversized reads explicitly |
| `resolveSavedWorkflowAsset(store, entity_ref)` | Return a validated owned file reference for display/export; no plot refitting or construction |
| `releaseWorkflowEntityCache(store, entity_id = NULL)` | Drop owned resident cache entries; never delete persistent history or required products |
| `workflowEntityCacheInfo(store)` | Report entries, accounted bytes, hits/misses and owned resources without hydrating |
| `closeWorkflowEntityStore(store)` | Idempotently release owned resources; accepted entities remain reopenable |

Thin `saveWorkflowCheckpoint()` and `hydrateWorkflowCheckpoint()` wrappers can
use this shared contract with `entity_kind = "checkpoint"`. Do not build a
separate store/index/recovery implementation for each consumer. The initial
entity kinds are package-owned and explicit; arbitrary classes, callbacks or
user-defined query code do not become a new extension framework.

The implementation may use one small R6 owner for mutable resources. Keep the
scientific transformations in their existing functions/S4 methods. No package
store method should fit PCA, normalize an assay or select a new RUV population.

Saving, finding and hydrating are deliberately explicit operations. Do not begin
with a convenience function that silently recomputes science on a lookup miss.

### 5.3 Checkpoint and reference contracts

A checkpoint record must preserve the original S4 plus its historical metadata.
The bounded reference should include:

- Schema/codec version and opaque checkpoint ID.
- Project and workflow identity, plus a project-incarnation identity that changes
  on replacement rather than merely hashing an unchanged directory name.
- Logical state name, parent checkpoint, and committed/checkpoint revision.
- Project-relative payload location, byte size and integrity binding.
- S4 class and small row/sample/assay summaries.
- Typed lookup factors and relevant scientific dependency identities.
- Writer/package versions and preserved provenance references.

The store must reject cross-project references, unsupported codecs, malformed
metadata, ambiguous ownership, path escapes and corrupted payloads. Full external
restore/import validation remains required. Any reuse of an internally verified
payload must have an explicit ownership/trust scope; a matching size alone is
not an external-integrity check.

The same ownership, version and integrity envelope applies to all saved entities.
Entity-specific metadata distinguishes a checkpoint from a fit, table, scientific
plot product or rendered asset; matching method labels alone cannot equate them.

The first codec should use one object record and a small index rather than
externalizing every metadata leaf. No thousands-of-files design before a workload
demonstrates why it is needed.

### 5.4 Lookup by factors: exact metadata semantics

Here "factors" means named retrieval criteria, not necessarily R `factor`
columns. Examples: workflow, state role, assay, source checkpoint, normalization
method, RUV settings, design identity, contrast or product kind.

Freeze which factors determine compatibility for each record kind before storing
or querying it:

- Unknown fields and invalid types are errors.
- Missing criteria mean "unspecified"; an explicitly supplied `NA` has a
  defined missing-value meaning and is not a wildcard.
- Preserve types, sample order and factor levels where they affect science.
  Do not collapse them into an untyped concatenated string.
- An exact checkpoint ID is authoritative for historical restore.
- A multi-match query returns bounded matches; the caller must explicitly resolve
  ambiguity. A cache miss is a normal no-match result, not a failed restore.
- Do not silently pick a filename containing `latest` or reuse a product
  merely because dimensions and method labels match.
- Scientific product keys include effective data, selections/value patches,
  relevant design, parameters and implementation version.
- Render keys additionally include style, format, dimensions and renderer
  version. A style change must not invalidate an unchanged scientific product.
- Start with conservative checkpoint identities. Cross-generation semantic reuse
  is an additional tested feature, not a prerequisite for the first proof.

A metadata lookup is different from a row-level query on a saved assay. The
`readSavedWorkflowTable()` contract declares which codecs can select before
collection. Its implementation is required for any selective-query treatment
being evaluated. Hydrating a whole RDS and then slicing it must never be described
as a bounded selective disk read. Keep an explicit whole-object comparator for
codecs without projection support.

### 5.5 Publication and recovery expectations

Use explicit immutable IDs and the smallest publication protocol that preserves
the existing checkpoint contract. A checkpoint is readable only when its payload
and accepted record agree.

Test interrupted writes, incomplete metadata, restart, retry and failure after
publication but before acknowledgement. Recovery must work from persistent
evidence, not a live closure containing per-file progress. Retrying the same
request must not duplicate a checkpoint. Conflicting reuse of an ID is an error.

Resident eviction, persistent deletion and project cleanup are separate
operations. This initial API needs cache release, not a general history-deletion
facility. Required outputs, parents and pending-recovery inputs cannot be
removed to satisfy a memory target.

### 5.6 Consumer lookup and reuse map

Design these consumers and standalone comparisons up front. They use the same
utility foundation and do not require Shiny to execute. Scientific entity
production remains in the existing scientific methods; the store persists and
retrieves their outputs.

| Consumer | Smallest sufficient saved entity | Required comparison |
| --- | --- | --- |
| Revert/reset | Exact checkpoint and its history metadata | Resident S4 versus one complete saved-state hydrate |
| QC overview/table preview | Saved summary or bounded rows/columns | Recompute/restore whole state versus summary lookup/projected read |
| PCA/RLE/correlation panels | Scores, explained variance, RLE statistics or sample-pair statistics | Full input hydrate plus preparation versus saved scientific product plus identical rendering |
| Density panel | Its own correctly defined PCA/summary product | Preserve its different scientific population; do not substitute the PCA panel's cached scores |
| Composite QC/export | Matching individual rendered assets or bounded products | Rebuild panels/composite versus asset reuse, including any required copies |
| DA analysis | Exact model inputs, reusable fit or complete accepted results, depending on request | Resident-input fit, full-S4 hydrate/fit, projected-input fit and compatible saved-fit/result reuse |
| DA table/volcano | Selected contrast results and required annotation columns | Whole DA object restoration versus bounded results lookup; fit count must remain zero on a compatible hit |
| DA heatmap | The declared feature/sample matrix and clustering inputs/results | Whole S4 restoration versus exact selected matrix and compatible clustering reuse |
| Existing report/session file consumers | Existing required files and exact metadata | Rebuild/reserialize versus validated existing output reuse, without omitting files |

Reuse at the highest compatible level. An unchanged request may reuse a rendered
asset; changed aesthetics may reuse a scientific product; a changed model may
reuse only unchanged input data. A miss returns explicit metadata so the consumer
can call the existing producer, then save its result. A stale or corrupt record
must not be quietly treated as a trustworthy hit.

DA compatibility factors include effective quantity identity, feature/sample
order, design/formula, blocking/replicate structure, weights, contrast definitions,
normalization/RUV state, fitting parameters and relevant software versions.
Define separate dependencies for the base fit, contrast statistics, multiple-test
adjustment and rendered view. A changed contrast is not automatically eligible
for every previously saved statistical product.

Filtering a displayed results table may reuse already computed statistics.
Changing the tested feature/sample population requires the appropriate scientific
recomputation. Never fit or adjust p-values only on the displayed page and present
them as the original full-population result. Keep ordinary whole-object callers
working through the existing public API.

The first implementation may cover entity kinds sequentially, but all consumers
being claimed at the GUI gate must already have standalone scientific and
performance evidence. No unmeasured generic auto-selector decides the access
route: freeze explicit consumer policies from observed results.

## 6. Synthetic fixtures: large, realistic and independent of the GUI

### 6.1 Reuse the actual scientific contracts

Build valid package S4 objects with the constructors and field names present in
the selected janitor baseline. Do not wrap random matrices in fake S4 classes.

The first fixture family must cover:

- DIA-NN-like long peptide/precursor data and associated peptide matrices.
- Protein quantification tables, protein-ID/annotation tables and sample design.
- Biological groups and technical replicates with enough observations for the
  actual supported methods.
- Q-value fields, peptide/protein mappings, quantity scales and meaningful
  missingness required by the existing DIA fixtures.
- Small inline args, substantial provenance/support members, and real Limpa
  support objects where that route produces them.
- Real class transitions from peptide to protein objects.
- Scientific plot products, actual fitted models, contrast-result tables,
  annotation projections and rendered assets generated from these same inputs.
  Test both small derived products and large result/matrix entities.

Review existing small fixtures for compatibility, including
[protein normalization fixtures](../../tests/testthat/helper-module-ci-prot-norm.R)
and [peptide QC fixtures](../../tests/testthat/helper-module-ci-prot-peptide-qc.R).
The performance branch's component fixtures are comparison/salvage material, not
dependencies to drag into the janitor baseline.

Every generated object must pass class validity and the relevant reader/design
alignment checks. Unsupported edge inputs belong in explicit rejection tests.

### 6.2 Scale by allocation and shape, not row count alone

| Tier | Purpose | Required design |
| --- | --- | --- |
| Small exactness | Fast default tests | Empty/singleton boundaries, rich types, small histories and failures |
| Representative large | Demonstrate useful memory relief | One current S4 fits comfortably, but several independently changed historical states retain substantial memory |
| Shape contrast | Expose layout-dependent costs | Tall peptide tables, wide protein matrices, annotation-heavy records and varying sample counts |
| History stress | Test bounded residency | Fixed current-state scale while history grows; exercise at least ten successive operation/revert/style/export cycles |
| Capacity boundary | Establish limits | Largest declared supported case, oversized single object, and history that makes the baseline exceed its limit |

Concrete dimensions, seeds, expected decoded bytes and history depths must be
frozen in scenario files before timed comparisons. Select large cases from
representative DIA shapes and a resource preflight, not from whichever sizes make
a codec look fast. Exact sizes remain to be set on the chosen baseline/host.

Do not allocate the full stress history merely to generate fixtures. Create inputs
incrementally, preserve seed/recipe receipts, and separate fixture generation
from timed execution. Include input loading in any cold-start boundary that
claims to measure it.

Generate both shared-history controls and genuinely independently transformed
quantities. `rep(list(object), n)` does not prove the cost of n independent
states. Avoid highly repetitive synthetic values that compress unrealistically;
vary annotation repetition and numeric entropy deliberately.

Do not leave the oracle, fixture builder or all expected states in a timed
process. Verify outputs in a separate untimed process with exact fixture bindings.
Default tests generate only small data. Large runs are explicitly invoked on a
validated disk-backed output root; generated datasets and run directories stay
out of package source and Git.

### 6.3 Required scenario coverage

| ID | Scenario | Required evidence |
| --- | --- | --- |
| C01 | Save, close/reopen, hydrate each supported S4 class | Exact object and metadata round-trip; usable S4 dispatch |
| C02 | Find by method/design/assay/source factors | Correct match/no-match/ambiguity; zero payload reads |
| C03 | Small versus large objects; zero optional cache | Correct uncached behavior; no silent cap increase |
| C04 | Unchanged/shared metadata versus changed quantities | Memory accounting distinguishes shared allocations; changed science never hits an old product |
| C05 | Row/column subsets, reordering, duplicate selectors and patches | Defined selector semantics, exact IDs/order and safe cache invalidation |
| C06 | Factors/ordered factors, names, row names, annotations and nested args | Preservation under the codec; explicit unsupported-value errors |
| C07 | NA, NaN, infinities, zero, all-missing and constant features | Preserve values and existing unavailable-panel/rejection behavior |
| C08 | Forward history, immediate revert, deep revert, reset and reapply | Same state/history semantics; only requested state hydrated |
| C09 | Revert then create a new branch of history | Old required history remains recoverable; no stale active refs |
| C10 | Cache hit, eviction, reload and changed style | No refit on valid hits; correct invalidation; byte ceiling respected |
| C11 | Original owners released versus deliberately retained aliases | Memory savings occur only when all relevant references are released |
| C12 | Many checkpoints and metadata matches | Bounded metadata lookup and resident cache; no full-history hydrate |
| C13 | Failed save, corruption, truncation, missing payload or invalid codec | Classed failure; previous active state/history remains usable |
| C14 | Process exit during save or before acknowledgement | Persistent recovery and idempotent retry in a fresh process |
| C15 | Cancellation and project replacement | Stale results cannot activate/publish into the replacement project |
| C16 | Two stores/projects and conflicting writers | No cross-project reuse; reject or serialize conflicts explicitly |
| C17 | Permission failure and disk-full injection | No falsely accepted checkpoint, lost previous output or orphaned authority |
| C18 | Complete legacy export and timestamped/latest outputs | Full historical contents and consumer compatibility; export peak included |
| C19 | Ten normalize/revert/style/export cycles in one process | Resources/cache bounds checked each cycle; retain every RSS-guard failure |
| C20 | Restart with no resident cache | Cold restore cost and compatibility measured independently of cache hits |
| C21 | Same DA request, fit/result reuse and changed contrast | No fitting on a fully compatible result hit; correct invalidation/recomputation at each statistical layer |
| C22 | Changed design, weights, sample order or tested feature population | Reject incompatible fits/results even when dimensions or filenames match |
| C23 | Plot style-only change, scientific change and asset hit | Separate zero-fit, zero-source-hydration and zero-render expectations according to the actual hit level |
| C24 | DA result pages, volcano projections and heatmap selections | Correct ordered values and full-population statistics; actual bounded reads for capable codecs |
| C25 | One-off access versus frequent reuse; hit/miss mixtures | Include initial production/persistence, lookup, cold read, invalidation and eviction costs |
| C26 | Combined consumer sequence | Saving for one consumer cannot shift more cost or peak memory into another consumer without reporting it |

Run all small supported correctness cases. For expensive large cases, freeze a
coverage plan spanning size, shape, history depth, mutation, cache state and
access pattern, with explicit high-risk combinations. Do not invent an exhaustive
cross-product or omit slow cases after seeing results.

## 7. Standalone scientific replay before Shiny integration

The utility harness and the scientific harness answer different questions.

- The utility harness proves storage, identity, metadata lookup, lifecycle and
  memory behavior using valid objects.
- The scientific harness runs existing package methods on original versus
  saved/hydrated/projected inputs and compares their outputs. It also compares
  saved fits/results and plot products against fresh production-method outputs.
- The GUI harness later proves that actual consumers release references and use
  the accepted APIs correctly.

Begin with DIA/IQ, then real Limpa. Cover all existing normalization choices,
automatic/manual/skipped RUV, correlation Apply/Skip/reapply, class transitions,
checkpoints and revert/reset before claiming full DIA scientific support. Freeze
the route matrix from the actual baseline methods rather than guessing defaults.

Use real IQ/Limpa for scientific acceptance. Mocks are appropriate for injecting
storage failure, not for replacing the estimator whose compatibility is claimed.
Seeds, tolerances, feature/sample populations and substantive provenance remain
fixed. Storage-only round-trips should be exact; floating scientific comparisons
use the existing documented tolerances.

Standalone consumer scenarios must include first DA execution, repeated analysis
with unchanged parameters, contrast changes, table paging, static/interactive
plot selection, heatmap selection, plot style changes and all required automatic
exports. Classify each comparison as access-only, scientific reuse, render reuse
or a combined sequence. Do not call skipping a required computation an I/O gain.

A shared transformed input may be reused across adjacent operations where its
lifetime and bytes are accounted for. Do not preserve unrelated assays merely
to pass a read-count assertion.

## 8. Measurement protocol and decision gates

### 8.1 Comparators and boundaries

Use the same fixture and operation sequence for:

- **B0-original:** selected janitor DIA workflow out of the box, measured before
  performance changes and preserved even if correctness checks fail.
- **B0:** qualified janitor DIA control. This is B0-original when it passes the
  scientific contract; otherwise it includes only separately recorded essential
  correctness repairs, with a new source identity and a new baseline measurement.
  Never silently replace B0-original or use a scientifically invalid control to
  support a speedup.
- **B1:** current-in-memory / cold-history candidate using RDS, isolating history
  residency from downstream reuse.
- **B2:** compatible saved scientific entity or rendered-asset reuse for one
  named consumer against B0/B1's best equivalent path.
- **B3, only if justified:** one alternative codec or selective-query adapter
  against the accepted mechanism.

For DA/plot/table consumers, explicitly compare applicable access paths:

1. Operate on an already-resident S4 or result object.
2. Load a complete RDS/S4 and run the existing consumer.
3. Load a smaller saved fit, result table or scientific plot product.
4. Read only the required data from a projection-capable saved table.
5. Reuse the exact already-rendered asset when the display/export request matches.

Resident access is the control where the application already has the object.
Do not force it through an artificial RDS round-trip to make lookup appear faster.
Likewise, rendering a saved image and fitting/rendering a new plot answer different
questions unless the scientific/display requests are equivalent and production
cost is included in a separate full-sequence comparison.

For each pair, keep active working state, mandatory outputs and process topology
the same. If topology changes, measure that as a separate treatment, including
worker startup, input transport and shutdown. Do not give the memory baseline
extra retained oracle objects or spare full histories.

Before measurements, record build hashes, fixture recipes/input hashes,
dependencies, R/native-library versions, storage device/filesystem, settings,
thread counts, resource limits, cache state and exact boundaries. Use fresh
installed-package execution. Reuse existing profiling/measurement tools only
after verifying they run on the selected baseline without importing the whole
performance architecture.

Freeze replication count, member order, aggregate/statistical summary and retry
rules before outcomes. Run timed members serially with profiling disabled and no
competing owned builds/tests. Separate diagnostic profiling from acceptance;
include unattributed time and overlapping peaks explicitly.

### 8.2 Measure the costs separately

| Metric | Required interpretation |
| --- | --- |
| Checkpoint save time and CPU | Includes encoding, compression, required verification and acceptance |
| Metadata lookup latency | No payload hydration; hit, miss and ambiguous match cases |
| Cold hydrate/revert time | Includes read, decode, validation and activation; not comparable to a resident hit without labeling |
| Resident hit time | Reports the benefit and cost of the declared cache |
| Saved consumer elapsed/CPU | Compare accessing input, retrieving a fit/product, rendering and completing the requested DA/plot/table action |
| Reuse break-even | Initial production/write cost versus savings per compatible reuse; include miss/invalidation frequency |
| Forward-sequence elapsed/CPU | Includes every required checkpoint write and automatic output |
| Retained memory by step | Unique reachable allocations where measurable; do not sum duplicate bindings |
| Peak process/service memory | Includes workers, buffers and relevant filesystem memory |
| Export elapsed and peak | Complete legacy history and all required copies remain inside the boundary |
| Disk bytes/files and I/O | Show whether lower RAM was exchanged for excessive writes/storage |
| Lifecycle counters | Saves, actual decodes, reads, cache hits, fits, open resources and checkpoint count |

Respect 8 GiB, zero swap, 900 seconds and existing registry/cache limits. Confirm
the output root is physical-disk-backed; on this host `/tmp` is unsuitable for
claiming disk offload. Record OS filesystem cache effects. "Cold" must specify
whether the R cache, process or OS page cache is cold; do not flush shared host
caches or alter system policy as a benchmark shortcut.

### 8.3 Stage decisions

| Gate | Required outcome | What it permits |
| --- | --- | --- |
| G0: baseline/fixture freeze | Recorded out-of-the-box DIA execution, qualified control, realistic objects, resource-safe cases and declared boundaries | Utility implementation and focused tests |
| G1: utility correctness | All supported exactness, lookup, lifecycle and failure contracts pass | Large standalone performance runs |
| G2: standalone benefit | Demonstrated memory benefit and acceptable elapsed/CPU for checkpoint sequences and each claimed DA/plot/table reuse treatment; every tradeoff reported | Candidate becomes eligible for scientific replay; no GUI performance claim |
| G3: scientific replay | Existing methods agree on original/hydrated/projected inputs and reused fits/results/products, including real IQ/Limpa for claimed scope | Narrow GUI integration |
| G4: DIA GUI acceptance | Complete journey, all outputs/restore, science, repeated-cycle stability and fixed performance gates | DIA-specific acceptance |
| G5: next omic | Its own fixtures, consumers and performance gates pass | Extension to that omic only |

The existing final journey target remains **at least 10% lower elapsed completion
time and 10% lower combined peak memory than the qualified control**. Any immediate
accepted comparator also requires no regression. This document does not relax
those thresholds or substitute retained heap for combined peak.

Before G2, freeze its numeric retained-memory target, allowable save/cold-revert
latencies, CPU criterion and target workload classes. Until those are specified,
G2 cannot be called accepted. A storage-only write need not be faster than an
in-memory assignment; the claim is evaluated over the declared complete sequence.
Memory relief with unacceptable elapsed cost is a tradeoff, not a dual win.

For repeated consumers, report initial entity production and each subsequent use.
Compare total measured cost over a predeclared realistic access sequence. Do not
assume unlimited cache hits to amortize an expensive save. Report the access count
at which reuse first pays back, if it does, and the resident memory/disk cost of
achieving that point. A one-off request may properly keep the original path.

Overall ecosystem gains mean lower total cost for the declared DIA producer and
consumer sequence, not an average of selectively chosen fast lookup calls. Show
per-operation regressions alongside the total, preserve required outputs, and
verify that composing individually accepted mechanisms does not introduce new
memory overlap or invalidation bugs. Later omics need separate measurements.

If the baseline exceeds its resource limit, report capacity improvement only;
do not invent a finite baseline runtime or compute a speedup from an aborted run.
If both fit, include all added storage costs. Small, low-history and frequent
cold-revert cases are necessary negative controls, not cases to hide.

Keep implementation conformance, scientific agreement, runtime completion and
performance acceptance separate. Utility gains do not predict the GUI's gains.
Synthetic gains do not establish representative biological-workload performance.

### 8.4 Failure and correction policy

Retain failed/unfavourable receipts. Diagnose before changing code. Each bounded
correction must name the observed defect or measured cost, its implementation
scope, affected tests, and the cases invalidated by that change.

Do not repeat unchanged cases to hunt for a pass, change thresholds after
outcomes, prewarm outside the boundary, remove required files, or repeatedly
change the architecture until a favourable number appears. If a checkpoint
format cannot meet the declared tradeoff, record a no-go for that use case.

## 9. Ecosystem treatments and evaluation priority

The saved-product, DA and selective-read cases below are part of the upfront
standalone design requested by the user. Implement them as separately measured
treatments on the shared foundation. Algorithm changes, new public formats and
other-omics GUI migrations remain outside the initial DIA proof.

| Candidate | Intended benefit | Compatibility boundary |
| --- | --- | --- |
| Small scientific products and rendered-asset reuse | Avoid repeated fits, plot builds and composite reconstruction | Core standalone consumer proof; exact dependencies and density's separate population |
| Serialize session once, copy completed bytes to latest | Reduce duplicate traversal/compression | Preserve both complete public files and measure the required copy/verification |
| Required-column projections and fewer matrix conversions | Reduce action working memory and repeated allocation | Same inputs/populations/annotations and substantive provenance |
| Carefully scoped data.table operations | Reduce measured join/aggregation costs | By-reference mutation cannot corrupt historical objects |
| DA fit/results reuse | Feed volcano, heatmap, tables and exports from accepted results | Core standalone consumer proof; preserve full contrast/FDR populations and invalidation |
| Selective Arrow/DuckDB queries | Avoid collecting large unused table portions | Upfront access-path comparison where applicable; preserve ordering/missingness and resource bounds |
| Lazy package loading and controlled parallelism | Reduce startup residency and redundant process/thread work | Keep required dependencies and S4/native registration; count startup |
| HDF5Array/DelayedArray | Address a required active matrix that exceeds memory | Requires compatible methods; no silent materialization or S4 slot substitution |

Evaluate history first to establish the storage/ownership lifecycle, then the
scientific-product, DA-result and rendered-asset reuse cases. Add a projected-read
codec where a consumer can benefit; measure serialization/copy and working-data
optimizations separately. Finally run the combined DIA sequence. A promising idea
or a component gain on the old branch is salvage evidence, not acceptance for
the restarted package.

## 10. Implementation handoff checklist

This document is the starting specification. Administrative preservation and
branch restoration are complete; scientific implementation and proof remain open.

The bounded execution breakdown is
[DIA-001–024](../../tickets/openai/01-DIA-SAVED-ENTITIES-PLAN.md).
All tickets begin pending; their creation is not implementation, scientific
agreement or performance acceptance. The first ticket freezes the actual DIA
consumer and public compatibility contract before any runtime change.

- [x] Archive dirty worktree, this document, selected fixtures and prior evidence;
  verify preservation before branch changes.
- [x] Resolve local janitor and restore the clean tracked development baseline.
- [ ] Measure and preserve out-of-the-box DIA/IQ behavior before performance
  refactors; qualify any essential correctness-repaired control separately.
- [ ] Inventory direct `states` consumers and every large-object owner.
- [ ] Freeze metadata factors, reference schema, codec and integrity contract.
- [ ] Map the checkpoint/product/fit/table/asset entities and their exact DA,
  plotting, QC, export and existing report consumers.
- [ ] Define concrete large fixture dimensions, seeds, memory limits and coverage.
- [ ] Freeze G2 numeric targets and repetition/order rules before timed runs.
- [ ] Implement the shared utility file and small exactness/failure fixtures.
- [ ] Validate a fresh installed package, exports, S4 dispatch and native routines.
- [ ] Run large standalone utility comparisons and record all costs.
- [ ] Compare resident access, full RDS/S4 restore, saved-fit/product lookup,
  projected reads and exact asset reuse for applicable consumers.
- [ ] Record hit/miss/invalidation behavior, reuse break-even and combined
  producer-consumer costs before selecting access policies.
- [ ] Run standalone scientific comparisons against unchanged methods.
- [ ] Integrate accepted primitives into DIA/IQ and measure the complete journey.
- [ ] Keep each GUI refactor tied to a named measured cost and its accepted
  standalone consumer; compare each treatment and the combined DIA result.
- [ ] Verify real Limpa and repeated-cycle DIA behavior before claiming DIA scope.
- [ ] Decide acceptance from the declared gates; hold other omics until their
  own evidence is available.

Every run receipt should include scenario ID, source/build/fixture bindings,
process/storage/cache settings, boundaries, repetitions, measurements, actual
output inventory, scientific comparison result, failures and claim scope.

Use a small decision record with: question, observation, chosen action, rejected
alternatives, supporting receipt, unresolved limitations and affected consumers.
Keep prior records immutable; append corrections with their reason.

## 11. Reference material

These sources informed the discussion. They are not benchmarks of MultiScholaR.
Recheck version-sensitive details against the implementation environment.

- [R single-object serialization](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html):
  RDS save/restore and compression choices.
- [R serialization internals](https://stat.ethz.ch/CRAN/doc/manuals/r-devel/R-ints.html):
  object/reference serialization behavior.
- [qs2](https://github.com/qsbase/qs2): optional R-serialization-compatible internal
  codec candidate; its restricted qdata format is not a general S4 replacement.
- [Arrow datasets](https://arrow.apache.org/docs/r/articles/dataset.html):
  filtering/projection before collection.
- [DuckDB memory management](https://duckdb.org/2024/07/09/memory-management):
  streaming, buffering and the cost of spilling.
- [DuckDB workload tuning](https://duckdb.org/docs/current/guides/performance/how_to_tune_workloads):
  workload-dependent thread and storage tradeoffs.
- [Shiny caching](https://shiny.posit.co/r/articles/improve/caching/):
  reuse and the correctness/cost of cache keys.
- [Shiny non-blocking operations](https://shiny.posit.co/r/articles/improve/nonblocking/):
  responsiveness, execution ownership and ExtendedTask.
- [data.table reference semantics](https://rdatatable.gitlab.io/data.table/articles/datatable-reference-semantics.html):
  mutation and explicit copying boundaries.
- [HDF5Array](https://bioconductor.org/packages/release/bioc/html/HDF5Array.html):
  disk-backed arrays and delayed/block-processed methods.
