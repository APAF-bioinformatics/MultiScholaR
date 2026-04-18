# Refactor Fidelity Audit Implementation Plan

## Objective

Build a deterministic, resumable audit system that can prove or clearly classify
parity between a refactored branch and a baseline branch for MultiScholaR.

The audit must produce:

- an authoritative machine-readable store suitable for automation
- append-only event history for inter-session resume
- clear human-readable reports for review and curation

The system must distinguish:

- exact-source parity
- normalized-text parity
- AST/code-shape parity
- behavior parity
- curated exceptions

It must not collapse those into one generic "match" result.

## Validated Findings From This Planning Pass

### Surface/API Prototype

Prototype comparison: `main` (`326c049`) vs current `janitor` `HEAD`.

Observed:

- `R/` files: `92` on `main`, `200` on `HEAD`
- unique top-level symbols: `615` on `main`, `1066` on `HEAD`
- unique `setMethod` targets: `65` on both sides
- unique `setClass` targets: `13` on both sides
- exported symbols in `NAMESPACE`: `482` on both sides

Implication:

- `NAMESPACE` parity is necessary but not sufficient
- internal helper surface changed massively even while exports remained stable
- audit must track internal symbol inventory, S4 inventory, and file/load-order
  metadata in addition to public API

### Manifest-Aware Block Fidelity Prototype

Prototypes covered five selector kinds:

- `symbol`
- `setMethod`
- `setClass`
- `expr_index`
- `anchor_range`

Prototype sample:

- `14` manifest entries across `5` manifests
- `12/14` were `rawExact`
- `14/14` were `normalizedTextExact`
- `14/14` were `exprExact`

Critical finding:

- `anchor_range` selectors were not resolvable against the target via the same
  selector logic used for the source
- those blocks required `normalized_text_search` in the target file
- the observed drift was whitespace-only, not semantic

Implication:

- source resolution and target resolution must be treated asymmetrically
- fidelity must be recorded at multiple levels:
  - `raw_exact`
  - `normalized_text_exact`
  - `ast_exact`

### Differential Behavior Prototype

Prototype differential replay against `main` vs `HEAD` passed for:

- `count_num_peptides`
- `countProteinsPerRun` on peptide-style input
- `countProteinsPerRun` on protein-quant-table-style input
- `findBestK`

Observed:

- pure/helper function replay is cheap and deterministic
- ggplot-backed helper replay is also feasible with lightweight fixtures

### S4 Isolated Replay Prototype

Prototype differential replay in isolated R processes confirmed parity for:

- `GridPlotData`
- `InitialiseGrid`

Critical finding:

- isolated processes are required for reliable S4 generic/method comparison
- R startup noise contaminated stdout ahead of JSON output

Implication:

- the final runner needs strict JSON-channel hygiene
- S4 behavior checks should run in subprocesses, not a shared interpreter

### Parsing/Infrastructure Findings

Observed during prototypes:

- naive regex parsing for `DESCRIPTION` `Collate` is not reliable enough
- CLI tools `extract_blocks.R` and `verify_refactor.R` are useful, but they are
  not reusable as libraries because they always execute `main()`

Implication:

- audit implementation needs a shared importable parser library
- `DESCRIPTION` and other package metadata should be parsed with real parsers,
  not regex scraping

## Implementation Checkpoint: 2026-04-17

`FA-P0` is now implemented.

Implemented artifacts:

- `tools/refactor/fidelity_inventory.R`
- `tools/refactor/fidelity_manifest.R`
- `tools/refactor/fidelity_inventory_snapshot.R`
- `tools/refactor/fidelity_audit.py`
- `tools/refactor/tests/fidelity-audit.test.cjs`

Existing tools now refit onto the shared libraries:

- `tools/refactor/extract_blocks.R`
- `tools/refactor/verify_refactor.R`

What this checkpoint now guarantees:

- parser-backed `DESCRIPTION` `Collate` extraction via `read.dcf`
- parser-backed `NAMESPACE` export inventory
- parser-backed top-level inventory for:
  - symbol assignments
  - `setGeneric`
  - `setMethod`
  - `setClass`
- durable audit bootstrap at `.refactor-fidelity-audit/` with:
  - SQLite at `.refactor-fidelity-audit/audit.db`
  - append-only JSONL at `.refactor-fidelity-audit/events.jsonl`
  - materialized summary and exception report placeholders
- machine-readable inventory capture via:
  - `python3 tools/refactor/fidelity_audit.py inventory`

Validation completed at this checkpoint:

- `node --test tools/refactor/tests/fidelity-audit.test.cjs`
- `python3 tools/refactor/fidelity_audit.py inventory --repo-root . --side target --target-ref HEAD`
- `Rscript --vanilla tools/refactor/fidelity_inventory_snapshot.R --repo-root .`

Regression coverage added in this checkpoint:

- inventory snapshot captures exports, `Collate`, symbols, and S4 entities
- audit bootstrap persists SQLite plus JSONL state and writes reports
- shared selector library is exercised through the real:
  - `verify_refactor.R`
  - `extract_blocks.R`
  entrypoints across:
  - `symbol`
  - `expr_index`
  - `setMethod` with signature
  - `anchor_range`

Next implementation target:

- `FA-P1`: inventory and package contract audit implemented

`FA-P1` is now implemented.

Implemented additions in this checkpoint:

- `tools/refactor/fidelity_audit.py surface`
- persisted file surface storage in `inventory_files`
- persisted export surface storage in `inventory_exports`
- persisted diff evidence in `inventory_diffs.evidence_json`

What this checkpoint now guarantees:

- baseline-vs-target surface comparison from either:
  - git refs
  - direct snapshot paths
- drift classification across:
  - missing definitions
  - extra definitions
  - duplicate definitions
  - moved definitions
  - definition hash drift
  - missing exports
  - extra exports
  - missing files
  - extra files
  - collate drift
- materialized human-readable summary for the latest surface run
- first-class diff rows in SQLite for automation and later curation

Validation completed at this checkpoint:

- `node --test tools/refactor/tests/fidelity-audit.test.cjs`
- `python3 tools/refactor/fidelity_audit.py surface --repo-root . --baseline-ref main --target-ref WORKTREE`

Live surface smoke result on this repo:

- baseline: `main@326c04904e504339a7ff15af00240a51cf337343`
- target: `WORKTREE`
- total diffs: `928`
- open diffs: `723`
- diff buckets:
  - `collate_drift`: `90`
  - `definition_drift`: `29`
  - `duplicate_definition`: `57`
  - `extra_definition`: `439`
  - `extra_file`: `108`
  - `moved_definition`: `205`

Immediate implication:

- the integrated branch is structurally far from the pinned pre-refactor
  baseline even where public exports still match
- `FA-P2` needs manifest-aware classification next so intentional moves can be
  separated from real drift

`FA-P2` is now implemented.

Implemented additions in this checkpoint:

- `tools/refactor/fidelity_manifest_compare.R`
- manifest catalog persistence in SQLite
- manifest comparison persistence in SQLite
- manifest exception candidate emission into `exceptions`
- `python3 tools/refactor/fidelity_audit.py manifest`

What this checkpoint now guarantees:

- manifest entries are cataloged into `manifest_entries`
- per-entry comparison evidence is recorded in `manifest_comparisons`
- comparison fidelity is tracked separately as:
  - `raw_exact`
  - `normalized_text_exact`
  - `ast_exact`
- target resolution is asymmetric and can use:
  - `selector`
  - `raw_text_search`
  - `normalized_text_search`
  - `ast_search`
  - file-level raw/normalized/AST fallback for dedicated targets
- non-exact or unresolved entries emit first-class exception candidates

Validation completed at this checkpoint:

- `node --test tools/refactor/tests/fidelity-audit.test.cjs`
- `python3 tools/refactor/fidelity_audit.py manifest --repo-root . --baseline-ref main --target-ref WORKTREE --manifest-path tools/refactor/manifest-wave1.yml --manifest-path tools/refactor/manifest-pept-s4-wave4.yml`

Live manifest smoke result on selected manifests:

- baseline: `main@326c04904e504339a7ff15af00240a51cf337343`
- target: `WORKTREE`
- manifests: `2`
- entries: `47`
- exceptions emitted: `31`
- status buckets:
  - `raw_exact`: `16`
  - `ast_only`: `1`
  - `content_missing`: `29`
  - `manual_merge_expected`: `1`

Important limitation discovered by the live smoke:

- many early-wave manifest entries now report `content_missing` because their
  original target files were later split again into second-generation helper
  files
- this is not a parser failure; it is a lineage-tracking gap in the manifest
  layer as currently designed
- that gap must be handled through later exception curation, or by a future
  transitive manifest-chain enhancement if the exception volume is too high

Immediate implication:

- `FA-P2` now separates exact matches from unresolved lineage, but it does not
  yet collapse multi-wave extraction ancestry into one proof chain
- `FA-P3` can proceed, but the later closeout will need either curated
  exceptions or chained manifest ancestry for these older waves

`FA-P3` is now implemented.

Implemented additions in this checkpoint:

- `tools/refactor/fidelity_behavior_replay.R`
- `tools/refactor/behavior-cases-smoke.json`
- `python3 tools/refactor/fidelity_audit.py behavior`
- behavior case catalog persistence in SQLite
- behavior result persistence in SQLite
- behavior exception candidate emission from replay mismatches

What this checkpoint now guarantees:

- behavior case catalogs are loaded from JSON into `behavior_cases`
- each baseline/target replay runs in an isolated `Rscript` process
- replay comparison tracks separately:
  - raw artifact equality
  - normalized result equality
  - warning equality
  - message equality
  - error equality
- behavior outcomes are classified distinctly as:
  - `exact_match`
  - `normalized_match`
  - `error_match`
  - `mismatch`
  - `runner_failure`
- latest behavior summary and exception reports are materialized beside the
  structural and manifest reports

Validation completed at this checkpoint:

- `python3 -m py_compile tools/refactor/fidelity_audit.py`
- `node --test tools/refactor/tests/fidelity-audit.test.cjs`
- `python3 tools/refactor/fidelity_audit.py behavior --repo-root . --baseline-ref main --target-ref WORKTREE --case-path tools/refactor/behavior-cases-smoke.json`

Regression coverage added in this checkpoint:

- exact scalar replay
- normalized-only replay for reordered tabular output
- warning mismatch detection
- message mismatch detection
- shared error equivalence
- isolated S4 replay via deterministic synthetic class/generic/method fixtures

Live behavior smoke result on this repo:

- baseline: `main@326c04904e504339a7ff15af00240a51cf337343`
- target: `WORKTREE`
- case files: `1`
- cases: `4`
- exceptions emitted: `0`
- status buckets:
  - `exact_match`: `4`

Important implementation nuance:

- isolated S4 replay is validated in the deterministic synthetic suite
- the current live smoke catalog is intentionally helper-focused, because
  sourcing whole legacy S4 owner files pulls unrelated package-load behavior
  into the replay surface and is better handled by curated case definitions in
  later phases

Immediate implication:

- `FA-P3` now gives the audit a decision-grade behavior layer for replayable
  helper surfaces
- `FA-P4` should focus on module contract integration and broader testthat
  mapping next

`FA-P4` is now implemented.

Implemented additions in this checkpoint:

- `tools/refactor/fidelity_contract_runner.R`
- `python3 tools/refactor/fidelity_audit.py contracts`
- testing matrix persistence in SQLite
- contract test file and case catalog persistence in SQLite
- optional contract execution result persistence in SQLite

What this checkpoint now guarantees:

- existing `tests/testthat` surfaces are cataloged into:
  - family
  - module family
  - target certainty tier
  - scenario tags
- the implementation plan testing matrix is materialized into
  `testing_matrix_entries`
- contract execution can run selected files in isolated `Rscript` processes and
  persist:
  - per-file status
  - per-file test count
  - per-file failure and skip counts
  - machine-readable result payloads
- latest contract summaries are materialized beside the structural, manifest,
  and behavior reports

Validation completed at this checkpoint:

- `python3 -m py_compile tools/refactor/fidelity_audit.py`
- `node --test tools/refactor/tests/fidelity-audit.test.cjs`
- `python3 tools/refactor/fidelity_audit.py contracts --repo-root . --target-path .`
- `python3 tools/refactor/fidelity_audit.py contracts --repo-root <synthetic-fixture> --target-path <synthetic-fixture> --execute`

Regression coverage added in this checkpoint:

- testthat family classification:
  - `module_contract`
  - `characterization`
  - `compat`
  - `golden_master`
  - `general_unit`
- scenario extraction for:
  - initialization
  - invalid input
  - export/report side effects
- execution result persistence for:
  - passing files
  - failing files
  - skipped tests inside otherwise passing files

Important implementation nuance:

- the original contract runner bug was a scope error in R:
  `reporter_results <<-` skipped the current function frame and left the local
  result object `NULL`
- the runner now uses local assignment and correctly records real `testthat`
  results instead of zero-test false passes
- full live execution on this repo remains optional, because many package-level
  tests require dependency-complete `load_all()` semantics that are outside the
  audit engine's deterministic default path

Live contract catalog result on this repo:

- target: `WORKTREE`
- status: `completed`
- matrix surfaces: `7`
- test files cataloged: `25`
- test cases cataloged: `747`
- family buckets:
  - `module_contract`: `6`
  - `characterization`: `2`
  - `compat`: `1`
  - `golden_master`: `1`
  - `general_unit`: `15`
- surface buckets:
  - `demonolithed_wrappers_and_modules`: `6`
  - `cross_module_and_workflow_segments`: `4`
  - `unspecified_general`: `15`

Synthetic execution validation result:

- status: `completed_with_failures`
- test files executed: `5`
- test cases cataloged: `7`
- execution buckets:
  - `passed`: `4`
  - `failed`: `1`
- skip handling is preserved inside the persisted file result payloads

Immediate implication:

- `FA-P4` now gives the audit a decision-grade contract and scenario mapping
  layer over the existing `testthat` estate
- `FA-P5` should focus on exception curation, disposition state, and closeout
  review workflows next

`FA-P5` is now implemented.

Implemented additions in this checkpoint:

- first-class exception keys and lineage metadata in the SQLite store
- aggregated exception ledger materialization in:
  - `.refactor-fidelity-audit/reports/latest-exceptions.json`
  - `.refactor-fidelity-audit/reports/latest-exceptions.md`
- `python3 tools/refactor/fidelity_audit.py exceptions`
- `python3 tools/refactor/fidelity_audit.py curate`
- exception emission from:
  - surface inventory diffs
  - manifest fidelity
  - behavior replay
  - contract execution failures

What this checkpoint now guarantees:

- exceptions are no longer mode-local blobs; they are durable, queryable ledger
  rows with:
  - stable exception keys
  - baseline and target lineage
  - curation status
  - disposition
  - owner
  - review note
  - updated timestamp
- auto-curation can safely close low-risk exception classes for:
  - whitespace-only drift
  - comment-only drift
  - target-resolution asymmetry
- manual curation can update exceptions by:
  - `exception_id`
  - `exception_key`
- latest exception reports can be materialized either as:
  - full ledger views
  - open-only closeout views
  - run-scoped review views

Validation completed at this checkpoint:

- `python3 -m py_compile tools/refactor/fidelity_audit.py`
- `node --test tools/refactor/tests/fidelity-audit.test.cjs`
- `python3 tools/refactor/fidelity_audit.py manifest --repo-root . --baseline-ref main --target-ref WORKTREE --manifest-path tools/refactor/manifest-wave1.yml --manifest-path tools/refactor/manifest-pept-s4-wave4.yml --run-id fa-p5-manifest-live`
- `python3 tools/refactor/fidelity_audit.py exceptions --repo-root . --run-id fa-p5-manifest-live --auto-curate`

Regression coverage added in this checkpoint:

- surface diffs emit curatable exceptions for open structural drift
- manifest exceptions support stable-key auto-curation
- contract failures emit curatable exceptions
- exception report generation captures:
  - total count
  - open count
  - high-severity open count
  - per-layer, per-type, per-status buckets
- manual curate updates by `exception_key`

Important implementation nuance:

- legacy audit stores needed an in-place schema migration fix because the new
  `exception_key` index cannot be created before the column exists
- the migration is now ordered correctly through `ensure_column(...)` followed
  by explicit index creation

Live curation smoke result on this repo:

- run id: `fa-p5-manifest-live`
- scope: selected manifests
  - `tools/refactor/manifest-wave1.yml`
  - `tools/refactor/manifest-pept-s4-wave4.yml`
- total exceptions: `31`
- open exceptions: `30`
- high-severity open exceptions: `29`
- auto-curated exceptions: `1`
- type buckets:
  - `missing_target_block`: `29`
  - `manual_merge_expected`: `1`
  - `comment_only_drift`: `1`

Immediate implication:

- the curation layer now makes the unresolved lineage gap explicit instead of
  burying it in prose
- `FA-P6` should focus on the integrated closeout run and readiness criteria
  for the post-parity coverage campaign

`FA-P6` is now implemented.

Implemented additions in this checkpoint:

- `python3 tools/refactor/fidelity_audit.py closeout`
- dedicated closeout report materialization in:
  - `.refactor-fidelity-audit/reports/latest-closeout.json`
  - `.refactor-fidelity-audit/reports/latest-closeout.md`
- integrated unresolved-exception materialization for closeout-scoped runs
- deterministic subcommand orchestration across:
  - surface
  - manifest
  - behavior
  - contracts

What this checkpoint now guarantees:

- one command can now run the integrated branch closeout audit against:
  - an explicit pinned baseline
  - an explicit mainline comparator
  - the target branch/worktree
- the closeout path auto-curates low-risk exception classes only for its own
  component runs
- the final closeout summary now records:
  - component run ids
  - per-layer readiness gates
  - blockers
  - open and high-severity exception counts
  - coverage-campaign readiness

Validation completed at this checkpoint:

- `python3 -m py_compile tools/refactor/fidelity_audit.py`
- `node --test tools/refactor/tests/fidelity-audit.test.cjs`
- synthetic integrated closeout fixture through the regression suite
- `python3 tools/refactor/fidelity_audit.py closeout --repo-root . --baseline-ref main --mainline-ref main --target-ref WORKTREE --manifest-path tools/refactor/manifest-wave1.yml --manifest-path tools/refactor/manifest-pept-s4-wave4.yml --case-path tools/refactor/behavior-cases-smoke.json --run-id fa-p6-closeout-live`

Regression coverage added in this checkpoint:

- integrated closeout orchestration across all audit layers
- auto-curated manifest drift closing a synthetic ready-to-ship fixture
- readiness report materialization in dedicated closeout artifacts
- final unresolved-exception report emission for the integrated run scope

Important implementation nuance:

- closeout intentionally emits a separate dedicated report instead of relying on
  whichever subcommand last touched `latest-summary.*`
- the real-repo smoke is expected to remain blocked until the structural and
  manifest lineage gaps are actually resolved; the closeout command is the
  decision artifact, not a paper-over layer

Live closeout smoke result on this repo:

- run id: `fa-p6-closeout-live`
- baseline: `main@326c04904e504339a7ff15af00240a51cf337343`
- mainline: `main@326c04904e504339a7ff15af00240a51cf337343`
- target: `WORKTREE`
- status: `blocked`
- component highlights:
  - baseline surface exceptions: `723`
  - manifest exceptions: `31`
  - behavior exceptions: `0`
  - contract exceptions: `0` with catalog-only execution
- readiness blockers:
  - open structural surface drift
  - open manifest lineage gaps
  - contract execution not requested
  - unresolved high-severity exceptions

Immediate implication:

- the fidelity system now has its final decision-grade closeout command
- the current branch is not yet ready to start the `>= 80%` coverage campaign
- the next practical work is no longer audit infrastructure; it is resolving the
  remaining metabolomics/general parity gaps that the closeout report now makes
  explicit

## Required Audit Layers

### Layer 0: Baseline Resolution

The audit must pin comparison anchors explicitly.

Required baseline modes:

- `branch_tip`: compare against a named branch, typically `main`
- `merge_base`: compare against merge-base of target branch and baseline branch
- `pinned_commit`: compare against an explicit commit SHA recorded in audit
  metadata

Rules:

- lane work should compare against a pinned baseline or merge-base
- post-merge closeout should compare the integrated refactor branch against both:
  - the pinned pre-refactor baseline
  - the active shared mainline branch

### Layer 1: Surface and Package Contract Inventory

This layer answers: "is the declared package surface still complete?"

Required surfaces:

- `NAMESPACE` exports
- `DESCRIPTION` `Collate`
- `R/` file inventory
- top-level function symbol inventory
- top-level `setMethod` inventory
- top-level `setClass` inventory
- top-level `setGeneric` inventory

Keying rules:

- symbol: symbol name
- setClass: class name
- setGeneric: generic name
- setMethod: generic name plus normalized signature payload

Outputs:

- additions
- removals
- moved definitions
- duplicate definitions
- signature drift
- export drift
- collate/load-order drift

### Layer 2: Manifest-Aware Block Fidelity

This layer answers: "for manifest-managed extractions, did the moved content
stay faithful?"

Input sources:

- manifest files under `tools/refactor/manifest-*.yml`
- baseline branch content
- current working tree content

Supported selector kinds:

- `symbol`
- `setMethod`
- `setClass`
- `setGeneric`
- `expr_index`
- `anchor_range`

Resolution rules:

- source side:
  - always resolve via manifest selector
- target side:
  - first try selector-based resolution where valid
  - then try exact raw text search
  - then try normalized-text search
  - record which target-resolution strategy succeeded

Required fidelity measures:

- raw text hash
- normalized-text hash
- AST/code fingerprint hash

Result classes:

- `raw_exact`
- `normalized_only`
- `ast_only`
- `content_missing`
- `selector_ambiguous`
- `manual_merge_expected`
- `semantic_drift_suspected`

### Layer 3: Differential Behavior Replay

This layer answers: "for code not proven by exact-source parity, does it still
behave the same?"

Behavior families:

- pure helpers
- dataframe/tibble helpers
- ggplot-backed helpers
- S4 methods
- module helper functions
- Shiny module contracts
- end-to-end workflow segments

Execution modes:

- same-process replay for pure helpers
- isolated subprocess replay for S4 methods and anything that mutates global
  method/class state
- `testServer` or characterization harness for module/server behavior

Comparison strategy:

- compare normalized return values
- compare warnings/messages where they are contractually relevant
- compare error class and message text when failure is expected
- compare side-effect artifacts for file-writing helpers

### Layer 4: Exception Curation

This layer answers: "if parity is not proven automatically, what exactly is the
reason and who needs to act?"

Exceptions must not live only in Markdown notes.

Required exception types:

- whitespace-only drift
- comment/roxygen-only drift
- target resolution asymmetry
- expected manual merge
- intentional API change
- intentional collate change
- unresolved semantic drift
- behavior mismatch
- fixture gap
- unsupported audit shape

Each exception record must include:

- stable key
- severity
- audit layer
- baseline target and current target
- exact reason
- evidence pointers
- curation status
- disposition

### Layer 5: Human and Machine Reporting

The audit must emit:

- machine-readable snapshots for automation
- human-readable summary for review
- diff-focused drilldowns for curated items

Report tiers:

- run summary
- branch summary
- manifest/wave summary
- unresolved exception report
- proven parity report

## Certainty Model

The system should explicitly track certainty tier per entity or block.

Proposed tiers:

- `T0_unseen`: not audited yet
- `T1_surface_present`: inventory/accounting only
- `T2_raw_exact`: raw content exact
- `T3_normalized_or_ast_exact`: formatting drift only, AST exact
- `T4_behavior_replayed`: differential behavior matched
- `T5_contract_replayed`: characterization/contract harness matched
- `T6_end_to_end_verified`: workflow-level parity verified

Release-closeout target:

- no unresolved high-severity exceptions
- every manifest-managed extraction reaches at least `T2` or an approved
  exception disposition
- every curated non-exact item reaches at least `T4` or an approved exception
  disposition
- critical module and workflow surfaces reach `T5` or `T6`

## Storage and Inter-Session Memory Design

### Authoritative Store

Use SQLite as the authoritative local audit database.

Proposed path:

- `.refactor-fidelity-audit/audit.db`

Rationale:

- local and durable
- queryable from Python/R/shell
- easy to checkpoint across sessions
- suitable for CI artifact publishing

### Append-Only Event Log

Use JSONL for append-only operational memory.

Proposed path:

- `.refactor-fidelity-audit/events.jsonl`

Each event should record:

- timestamp
- run id
- phase
- status transition
- touched artifact ids
- summary note

### Materialized Reports

Use deterministic materialization outputs:

- `.refactor-fidelity-audit/reports/latest-summary.json`
- `.refactor-fidelity-audit/reports/latest-summary.md`
- `.refactor-fidelity-audit/reports/latest-exceptions.json`
- `.refactor-fidelity-audit/reports/latest-exceptions.md`

### Progress Manifest

Keep a separate small machine-readable progress file for compaction/resume.

Proposed path:

- `tools/refactor/fidelity-audit-plan.json`

This file should track:

- implementation phase status
- last validated prototypes
- next implementation checkpoint
- blockers and open questions

## Proposed SQLite Schema

### `audit_runs`

One row per audit execution.

Fields:

- `run_id`
- `started_at`
- `completed_at`
- `repo_root`
- `baseline_ref`
- `target_ref`
- `mode`
- `status`
- `summary_json`

### `inventory_entities`

One row per discovered entity per side.

Fields:

- `entity_id`
- `run_id`
- `side` (`baseline` or `target`)
- `entity_kind`
- `entity_name`
- `signature_key`
- `file_path`
- `line_start`
- `line_end`
- `exported`
- `collate_index`
- `hash_raw`
- `hash_ast`

### `inventory_diffs`

Surface-level reconciliation records.

Fields:

- `diff_id`
- `run_id`
- `entity_kind`
- `entity_key`
- `diff_class`
- `baseline_entity_id`
- `target_entity_id`
- `severity`
- `status`

### `manifest_entries`

Normalized manifest entry catalog.

Fields:

- `manifest_entry_id`
- `run_id`
- `manifest_path`
- `entry_id`
- `selector_kind`
- `selector_payload_json`
- `source_path`
- `target_path`
- `group_name`
- `action`

### `manifest_comparisons`

One row per manifest entry comparison.

Fields:

- `comparison_id`
- `run_id`
- `manifest_entry_id`
- `baseline_source_resolver`
- `target_resolver`
- `raw_exact`
- `normalized_text_exact`
- `ast_exact`
- `hash_raw_baseline`
- `hash_raw_target`
- `hash_ast_baseline`
- `hash_ast_target`
- `status`

### `behavior_cases`

Catalog of replayable behavior fixtures.

Fields:

- `behavior_case_id`
- `family`
- `entity_key`
- `fixture_kind`
- `fixture_path`
- `normalizer_key`
- `risk_level`

### `behavior_results`

One row per behavior replay.

Fields:

- `behavior_result_id`
- `run_id`
- `behavior_case_id`
- `baseline_status`
- `target_status`
- `normalized_equal`
- `warning_equal`
- `error_equal`
- `artifact_equal`
- `result_json`

### `exceptions`

Curated deviations and unresolved issues.

Fields:

- `exception_id`
- `run_id`
- `audit_layer`
- `entity_key`
- `severity`
- `exception_type`
- `reason`
- `evidence_json`
- `curation_status`
- `disposition`
- `owner`
- `resolved_at`

## Implementation Architecture

### Shared Parser Library

Create an importable parser module instead of re-implementing selector logic in
multiple scripts.

Recommended modules:

- `tools/refactor/fidelity_inventory.R`
- `tools/refactor/fidelity_manifest.R`
- `tools/refactor/fidelity_compare.R`
- `tools/refactor/fidelity_behavior.py`
- `tools/refactor/fidelity_report.py`

Refactor existing selector logic out of:

- `extract_blocks.R`
- `verify_refactor.R`

so the audit can import the same logic instead of shelling out blindly.

### CLI Entry Points

Recommended commands:

- `python3 tools/refactor/fidelity_audit.py inventory`
- `python3 tools/refactor/fidelity_audit.py manifest`
- `python3 tools/refactor/fidelity_audit.py behavior`
- `python3 tools/refactor/fidelity_audit.py report`
- `python3 tools/refactor/fidelity_audit.py run`

### Skill Layout

Recommended skill name:

- `refactor-fidelity-audit`

Recommended shape:

- `SKILL.md`
- `references/schema.md`
- `references/curation.md`
- `scripts/run_audit.py`

The skill should orchestrate deterministic repo-local tools rather than
re-implement audit logic in prompt text.

### Skill Basis and Division of Responsibility

Yes, these audit layers form the basis of the skill.

The skill should be intentionally thin.

The skill owns:

- selecting the right audit mode
- choosing baseline strategy
- invoking repo-local audit engines
- deciding whether a finding is structural, behavioral, or curated
- pointing the user to the right human and machine reports

The skill should not own:

- parsing branch state ad hoc in prompt text
- re-implementing manifest resolution in prompt text
- re-implementing behavior replay logic in prompt text
- storing long-lived audit memory inside chat history only

The repo-local engines own:

- parser-backed inventory extraction
- manifest-aware block comparison
- behavior replay and normalization
- exception persistence and curation state
- report materialization

This separation is required so that:

- the skill remains stable across sessions and compaction
- the evidence remains machine-readable
- automation and CI can use the same underlying engines
- human review can trust that the same logic is used every time

## Rollout Plan

### Phase 0: Foundation

Deliver:

- parser-backed inventory engine
- parser-backed `DESCRIPTION` handling
- importable manifest resolver
- SQLite bootstrap
- JSONL event log bootstrap

Acceptance:

- inventory runs against `main` vs `HEAD`
- `Collate` parsed deterministically
- manifest selector resolution works for all current selector kinds on source

### Phase 1: Surface Audit

Deliver:

- inventory diff tables and reports
- export/collate drift detection
- duplicate and ambiguity detection

Acceptance:

- current `main` vs `HEAD` diff report emitted in DB, JSON, and Markdown

### Phase 2: Manifest Fidelity

Deliver:

- manifest entry normalization
- asymmetric target resolver
- raw/normalized/AST comparison
- exception generation for non-exact cases

Acceptance:

- sampled manifests prove the prototype result pattern observed in this planning
  pass
- all active selector kinds are handled

### Phase 3: Differential Behavior

Deliver:

- behavior case catalog
- replay harness for pure helpers, ggplot helpers, and isolated S4 methods
- result normalization and evidence capture

Acceptance:

- current sample behavior cases can be replayed automatically
- JSON output is clean even under repo startup noise

### Phase 4: Contract and Module Replay

Deliver:

- integration with existing `testthat` characterization, contract, compat, and
  golden-master tests
- module family selection manifest

Acceptance:

- representative module families can be mapped into fidelity tiers `T5` or `T6`

### Phase 5: Exception Curation Workflow

Deliver:

- curated exception DB flow
- human review report
- disposition rules

Acceptance:

- whitespace-only and target-resolution exceptions can be auto-curated
- semantic/behavioral mismatches are clearly escalated

### Phase 6: Final Closeout Audit

Deliver:

- integrated branch audit against pinned baseline and mainline
- final exception report
- readiness summary for coverage push

Acceptance:

- all high-severity fidelity exceptions resolved or explicitly accepted
- audit ready to gate the start of the 80 percent coverage campaign

## Existing Test Families To Reuse

The audit should not invent a parallel testing taxonomy. Reuse current families:

- characterization tests
- module contract tests
- compat tests
- golden master tests

Observed examples:

- `test-prot-01c-import-module-contracts.R`
- `test-prot-07b-da-handlers-compat.R`
- `test-prot-07c-da-results-characterization.R`
- `test-protein-da-golden-master.R`

These already describe the seam between structural parity and behavioral
parity. The audit should catalogue and reuse them.

## Testing Matrix and Coverage Policy

Coverage is required, but it is one evidence channel inside the fidelity
system, not the whole system.

The final system should record test evidence in the audit database and reports
alongside structural and behavior-parity evidence.

### Package-Level Targets

- package-wide post-closeout line coverage target: `>= 80%`
- touched extracted helper families target: `>= 90%`
- touched preserved helper families target: `>= 85%`
- wrapper/module files: no coverage-only acceptance; require scenario and
  contract completeness in addition to measured coverage

### Surface 1: Extracted Pure Helpers

Examples:

- count helpers
- summarizers
- table transformers
- small orchestration-free utility functions

Primary evidence:

- differential replay against baseline
- direct unit tests on the current branch

Fixture strategy:

- generated minimal fixtures
- representative edge-case fixtures
- NA/empty/invalid input fixtures

Acceptance:

- deterministic normalized output equality against baseline
- touched family coverage `>= 90%`
- no unresolved high-severity behavior mismatches

### Surface 2: Preserved but Touched Helpers

Examples:

- preserved functions whose code changed during seam extraction
- helpers that now delegate to extracted owners

Primary evidence:

- differential replay against baseline
- characterization tests on current branch
- unit tests for new branching paths

Fixture strategy:

- existing characterization fixtures
- regression fixtures for preserved bugfixes
- backward-compatibility fixtures where old aliases or old data shapes still
  matter

Acceptance:

- normalized output parity or approved exception
- touched family coverage `>= 85%`
- explicit exception record for any intentional drift

### Surface 3: S4 Classes, Generics, and Methods

Examples:

- `setClass`
- `setMethod`
- grid and QC S4 helpers

Primary evidence:

- inventory parity
- manifest-aware block fidelity
- isolated-process behavior replay

Fixture strategy:

- subprocess-per-variant replay
- slot shape assertions
- dispatch and return-class assertions
- representative object fixtures

Acceptance:

- class and dispatch surfaces preserved
- behavior parity proven in isolated replay where applicable
- touched S4 family coverage `>= 85%`
- JSON output channel remains clean under subprocess execution

### Surface 4: Extracted Wrapper Helpers

Examples:

- helpers split out of large module/server files
- action/render/state helper families

Primary evidence:

- helper-level unit tests
- differential replay where baseline extraction target is available
- side-effect assertions for notifications, file writes, and workflow-state
  mutations

Fixture strategy:

- minimal reactive stubs
- mock state managers
- file-system tempdirs
- notification/log capture

Acceptance:

- touched helper family coverage `>= 90%`
- explicit side-effect assertions for all public helper responsibilities
- no hidden dependency on original monolith-global state unless recorded as
  compat behavior

### Surface 5: De-Monolithed Wrappers and Modules

Examples:

- `mod_*` wrappers
- large server/ui modules
- workflow coordinators

Primary evidence:

- `testServer` contract suites
- compatibility tests
- scenario-based workflow tests

Required scenario classes:

- initialization and default state
- happy path
- missing prerequisite path
- invalid user input path
- state restore/session reload path
- reset/cleanup path
- export/report side-effect path
- error handling and notification path

Coverage policy:

- wrapper/module coverage is measured, but not used alone as acceptance
- scenario completeness is mandatory
- if wrapper coverage is low, that is acceptable only when helper coverage and
  module contract coverage prove the public behavior comprehensively

Acceptance:

- module contract suite green
- representative scenarios green
- no unresolved high-severity behavioral drift
- wrapper reaches certainty tier `T5` or `T6`

### Surface 6: Cross-Module and Workflow Segments

Examples:

- import to design handoff
- normalization to DA handoff
- report/export workflows

Primary evidence:

- characterization tests
- compat tests
- golden-master tests
- workflow smoke tests

Fixture strategy:

- pinned workflow fixtures
- captured checkpoints and RDS snapshots where appropriate
- representative omics-specific flows

Acceptance:

- critical workflows reach `T6_end_to_end_verified`
- golden-master or characterization evidence exists for all critical user
  journeys touched by the refactor

### Surface 7: Exception Curation Flow

Primary evidence:

- automated exception emission from structural and behavior layers
- human review disposition

Required exception buckets:

- whitespace-only drift
- comment-only drift
- target-resolution asymmetry
- intentional behavior change
- unsupported automation shape
- missing fixture
- unresolved mismatch

Acceptance:

- every non-exact or non-matching result is captured as a first-class exception
  record
- no important deviation lives only in prose or in chat history

### How Coverage Fits Into the Audit DB

Coverage should be recorded as evidence, not as a standalone dashboard only.

Recommended additional schema:

- `coverage_runs`
- `coverage_entities`

Suggested fields:

- package coverage percentage
- file coverage percentage
- touched file flag
- helper family label
- wrapper family label
- source test family

This allows the final report to answer:

- is the package at `>= 80%`
- are touched helper families at their stricter targets
- are low-coverage wrappers compensated by strong contract and scenario
  evidence

### Recommended Gating Order

Use this order for closeout:

1. structural inventory and manifest fidelity
2. differential helper replay
3. S4 isolated replay
4. module contract and scenario replay
5. package-wide coverage pass
6. final exception curation and closeout report

Do not start the broad `>= 80%` coverage campaign before the structural and
behavioral fidelity layers are in place, or the test effort risks locking in
accidental drift.

## Progress Checkpoints For Future Sessions

Use these stable checkpoint ids when updating progress:

- `FA-P0`: parser library extracted and reusable
- `FA-P1`: inventory and package contract audit implemented
- `FA-P2`: manifest fidelity engine implemented
- `FA-P3`: behavior replay engine implemented
- `FA-P4`: module/contract harness integration implemented
- `FA-P5`: exception curation flow implemented
- `FA-P6`: final integrated audit run completed

This plan file and the adjacent JSON manifest should be updated together as
implementation progresses.
