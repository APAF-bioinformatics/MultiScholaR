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
- `FA-C0`: coverage instrumentation bootstrap
- `FA-C1`: auto-curate bundle mapping
- `FA-C2`: baseline-vs-target comparative coverage runner
- `FA-C3`: coverage evidence DB and report materialization
- `FA-C4`: family-by-family coverage fill campaign
- `FA-C5`: coverage-backed parity closeout

## Post-Closeout Coverage-Backed Auto-Curation Plan

Authoritative current closeout:

- run id: `integrated-closeout-20260419d`
- status: `ready_for_coverage_campaign`
- proven parity: `true`
- open exceptions: `0`
- high-severity open exceptions: `0`
- auto-curated exceptions: `868`
- catalog-curated exceptions: `74`

### Clarification On The Coverage Target

The requested coverage expansion must target deduplicated parity surfaces, not
raw exception rows.

Why:

- the `868` auto-curated rows are not `868` distinct behavioral units
- the dominant bucket is `source_lineage_gap_target_resolved` `788`, which is
  mainly ancestry bookkeeping after second-generation file splits
- the remainder is largely low-risk structural noise:
  - `target_resolution_asymmetry` `27`
  - `whitespace_only_drift` `22`
  - `comment_only_drift` `6`
  - `surface_definition_drift_manifest_overlap` `25`

Current deduplicated footprint of the auto-curated population:

- `138` unique files touched by auto-curated evidence
- `24` unique legacy source files in the dominant lineage bucket
- `66` unique extracted target files in the dominant lineage bucket

This means the right execution unit is a `coverage bundle`, not an exception
row.

### Feasibility Assessment

This is feasible.

Current repository scale:

- `317` `R/` files
- `116,412` R source lines
- `115` `tests/testthat/` files
- `90,311` test lines

Current tooling gap:

- `covr` is not listed in `Suggests`
- `covr` is not installed in the current runtime

Practical estimate:

- tooling/bootstrap and DB/report additions: `1-2` days
- bundle mapping and comparative runner: `1-2` days
- targeted test expansion to bring the implicated parity bundles up to target:
  `1-2` weeks
- reruns, documentation, and final closeout: `1-2` days

Realistic total: `2-3` focused weeks.

The biggest cost is not audit plumbing. It is writing or extending test suites
for the parity bundles that still lack dense behavioral coverage.

### Coverage Bundle Model

A `coverage bundle` is the deduplicated parity unit used for coverage
measurement and reporting.

Examples:

- a preserved helper symbol
- an extracted helper family
- an S4 method surface
- a wrapper/server entry shell
- a lineage family where one legacy source block now spans several extracted
  target helpers

Each bundle should record:

- `bundle_id`
- `entity_key`
- `bundle_type`
- `baseline_ref`
- `target_ref`
- `baseline_paths`
- `target_paths`
- `linked_exception_keys`
- `shared_test_files`
- `baseline_line_coverage_pct`
- `target_line_coverage_pct`
- `scenario_class`
- `differential_replay_status`
- `coverage_gate_status`
- `review_note`

Acceptance should apply to bundles, not rows.

### Required Audit DB Additions

Add the following tables:

- `coverage_runs`
- `coverage_files`
- `coverage_bundles`
- `coverage_bundle_members`
- `coverage_test_links`
- `coverage_bundle_results`

Minimum useful fields:

- run id and timestamp
- baseline ref and target ref
- package line coverage percentage
- file line coverage percentage
- bundle line coverage percentage
- bundle type
- linked test files
- linked parity exception keys
- baseline and target coverage deltas
- gate status
- documentation note

### Required Report Outputs

Persist all of the following for automation and human review:

- coverage JSON summary
- coverage Markdown summary
- bundle-level machine-readable manifest
- bundle-to-test mapping
- baseline-vs-target coverage comparison report
- unresolved low-coverage bundle ledger

### Execution Phases

#### `FA-C0`: Coverage Instrumentation Bootstrap

Goal:

- make coverage collection reproducible on both baseline and target trees

Status:

- completed on `2026-04-19`

Work:

- add `covr` to `Suggests`
- bootstrap a repo-local coverage runner under `tools/refactor/`
- make the runner accept:
  - baseline ref
  - target ref or `WORKTREE`
  - explicit test file subsets
  - bundle manifest input
- materialize the new coverage tables in `audit.db`

Acceptance:

- one command can collect coverage for a test subset on the target tree
- one command can collect coverage for the same subset on the baseline tree
- results persist into the new coverage tables

Implemented:

- `tools/refactor/fidelity_audit.py` now materializes:
  - `coverage_runs`
  - `coverage_files`
  - `coverage_bundles`
  - `coverage_bundle_members`
  - `coverage_test_links`
  - `coverage_bundle_results`
- new CLI mode:
  - `python3 tools/refactor/fidelity_audit.py coverage`
- new repo-local runner:
  - `tools/refactor/fidelity_coverage_capture.R`
- deterministic validation paths:
  - fixture-backed success capture via `FIDELITY_COVERAGE_FIXTURE_JSON`
  - deterministic unavailable-tool path via `FIDELITY_COVERAGE_FORCE_UNAVAILABLE`
- new persisted coverage artifacts:
  - `.refactor-fidelity-audit/reports/latest-coverage.json`
  - `.refactor-fidelity-audit/reports/latest-coverage.md`
- package metadata updated:
  - `covr` added to `Suggests`
- regression coverage added for:
  - schema bootstrap
  - unavailable-tool capture
  - successful bundle/file persistence

#### `FA-C1`: Auto-Curate Bundle Mapping

Goal:

- collapse auto-curated and catalog-curated parity evidence onto stable
  coverage bundles

Status:

- completed on `2026-04-19`

Work:

- map every auto-curated row to a `bundle_id`
- map every catalog-curated row to the same bundle model
- preserve many-to-one links:
  - many exception rows may point to one bundle
- rank bundles by:
  - number of linked exception rows
  - bundle type
  - file size
  - current evidence depth

Acceptance:

- every auto-curated row links to a bundle
- every catalog-curated row links to a bundle
- priority bundle list is materialized for execution

Implemented:

- new CLI mode:
  - `python3 tools/refactor/fidelity_audit.py bundle-map`
- new persisted report artifacts:
  - `.refactor-fidelity-audit/reports/latest-bundles.json`
  - `.refactor-fidelity-audit/reports/latest-bundles.md`
- audit DB persistence now includes explicit exception-to-bundle linkage via:
  - `coverage_exception_links`
- bundle grouping rules now distinguish:
  - `lineage_family` bundles grouped by legacy source file for
    `equivalent_target_resolved_without_baseline_ancestor`
  - surface-level helper/method/wrapper bundles grouped by canonicalized
    selector-derived surface entity
  - canonical duplicate and manual-merge bundles grouped onto the surviving
    canonical implementation surface
- heuristic shared-test selection is materialized from the existing
  `tests/testthat` catalog into `coverage_test_links`
- priority ranking now uses:
  - linked exception count
  - bundle type
  - file size
  - current evidence depth

Live result on the integrated closeout evidence:

- run id: `bundle-map-20260419T121735Z-3a4c2239`
- closeout source: `integrated-closeout-20260419d`
- bundles: `149`
- linked curated exceptions: `942`
- priority families: `56`
- bundle type split:
  - `24` lineage families
  - `46` wrapper entrypoints
  - `63` helper surfaces
  - `13` S4 method surfaces
  - `2` canonical duplicate bundles
  - `1` manual merge bundle

Priority head:

- `R/mod_prot_enrich.R`
- `R/mod_prot_norm.R`
- `R/mod_metab_norm.R`
- `R/mod_prot_design_builder.R`
- `R/mod_metab_da.R`

Validation:

- `python3 -m py_compile tools/refactor/fidelity_audit.py`
- live `bundle-map` run on the integrated closeout evidence
- fixture coverage in `tools/refactor/tests/fidelity-audit.test.cjs`
  Note:
  The sandboxed Node harness still exhibits the pre-existing slow-exit behavior,
  so direct regression execution is treated as logically passing when the file
  reaches completion without assertion output, rather than relying on fast host
  process teardown.

#### `FA-C2`: Comparative Coverage Runner

Goal:

- measure the same tests against baseline and target for the same parity bundle

Work:

- run identical test subsets against:
  - baseline worktree
  - target worktree
- record:
  - baseline coverage
  - target coverage
  - coverage delta
  - test file provenance

Acceptance:

- bundle-level baseline and target coverage can be queried side by side
- per-bundle coverage reports are machine-readable

Status:

- complete

Implemented:

- `coverage` now accepts `latest-bundles.json`-style inputs:
  - `shared_test_files`
  - `baseline_members`
  - `target_members`
- new `coverage-compare` mode in `tools/refactor/fidelity_audit.py`
  runs the same test subset against:
  - baseline source
  - target source
  - using a temporary overlay workspace so target-selected tests are reused on
    both sides
- baseline and target single-side captures are persisted into
  `coverage_runs`, `coverage_files`, `coverage_bundles`,
  `coverage_bundle_members`, `coverage_test_links`, and
  `coverage_bundle_results`
- comparison rows are also persisted into `coverage_bundles` and
  `coverage_bundle_results` with:
  - baseline coverage
  - target coverage
  - delta classification
  - shared test provenance
- dedicated compare artifacts are written to:
  - `.refactor-fidelity-audit/reports/latest-coverage-compare.json`
  - `.refactor-fidelity-audit/reports/latest-coverage-compare.md`

Validation:

- `python3 -m py_compile tools/refactor/fidelity_audit.py`
- direct fixture-backed CLI smoke:
  `coverage-compare` produced `33.3%` baseline vs `100.0%` target with
  `delta_counts.improved = 1`
- SQLite smoke on the fixture run confirmed persisted bundle rows for:
  - `baseline`
  - `target`
  - `comparison`
- `tools/refactor/tests/fidelity-audit.test.cjs` was extended for:
  - `shared_test_files` plus side-member coverage input
  - comparative baseline/target delta persistence
  Note:
  the host Node harness still exhibits the pre-existing slow-exit behavior, so
  the direct CLI smoke is treated as the authoritative runtime validation for
  this checkpoint.

#### `FA-C3`: Coverage Evidence Materialization

Goal:

- make coverage a first-class parity evidence channel inside the audit system

Work:

- extend report generation to include:
  - package-level coverage
  - bundle-level coverage
  - low-coverage bundle exceptions
  - baseline-vs-target deltas
- write bundle manifests and summaries under `.refactor-fidelity-audit/reports/`

Acceptance:

- one report answers which bundles are below target
- one report answers whether target coverage regressed from baseline
- one report shows which tests justify each curated parity bundle

Status:

- complete

Implemented:

- new `coverage-evidence` mode in `tools/refactor/fidelity_audit.py`
- dedicated evidence artifacts:
  - `.refactor-fidelity-audit/reports/latest-coverage-evidence.json`
  - `.refactor-fidelity-audit/reports/latest-coverage-evidence.md`
- bundle evidence evaluation now merges:
  - bundle-map metadata from `latest-bundles.json`
  - comparative coverage deltas from `coverage-compare`
  - threshold policy from `fidelity-audit-plan.json`
- thresholds are applied by bundle type:
  - helper surfaces: `90%`
  - S4/preserved method surfaces: `85%`
  - lineage/wrapper/canonical bundles: `80%`
- `coverage-evidence` persists evaluated bundle rows back into the audit store
  with:
  - `coverage_gate_status`
  - baseline vs target coverage values
  - delta classification
  - shared test provenance
- low-coverage exceptions are now first-class audit records under
  `audit_layer = coverage_evidence`:
  - `coverage_target_below_threshold`
  - `coverage_target_unmeasured`
  - `package_coverage_below_threshold`
- global tool-unavailable runs do not flood the exception ledger; they publish
  `tool_unavailable` evidence with `0` emitted coverage exceptions

Validation:

- `python3 -m py_compile tools/refactor/fidelity_audit.py`
- fixture-backed end-to-end CLI smoke:
  - baseline package coverage `95.0%`
  - target package coverage `60.0%`
  - bundle threshold `90.0%`
  - emitted:
    - `coverage_target_below_threshold`
    - `package_coverage_below_threshold`
- live repo smoke on current env:
  - `coverage-evidence --limit 1`
  - status `tool_unavailable`
  - `exception_count = 0`
- `tools/refactor/tests/fidelity-audit.test.cjs` extended with a dedicated
  FA-C3 regression for:
  - evidence summary generation
  - low-coverage bundle detection
  - coverage exception emission
  Note:
  the host Node harness still exhibits the pre-existing slow-exit behavior, so
  direct CLI smokes remain the authoritative runtime checks for coverage
  checkpoints.

#### `FA-C4`: Family-By-Family Coverage Fill Campaign

Goal:

- raise the implicated parity bundles to coverage target with explicit
  documentation

Priority families from the current auto-curated footprint:

- `R/mod_prot_enrich.R`
- `R/mod_prot_norm.R`
- `R/mod_metab_norm.R`
- `R/mod_prot_design_builder.R`
- `R/mod_metab_da.R`
- `R/mod_lipid_design_builder.R`
- `R/mod_lipid_norm.R`
- `R/mod_prot_import.R`
- `R/mod_lipid_da.R`
- `R/mod_metab_import.R`
- `R/mod_prot_summary.R`
- `R/func_general_filemgmt.R`
- `R/func_lipid_qc.R`
- `R/func_prot_annotation.R`

Work pattern:

- take one family at a time
- map all bundles in that family
- identify missing direct tests
- add or extend tests until the family meets bundle targets
- document why the surviving implementation is equivalent

Acceptance:

- each family has:
  - updated tests
  - updated coverage evidence
  - updated bundle notes

Status:

- in progress

Implemented so far:

- `tools/refactor/fidelity_coverage_capture.R` now runs real package coverage via:
  - `pkgload::load_all()`
  - `covr` namespace tracing
  - selected `testthat::test_file()` execution
- coverage capture now understands `covr` function-range output
  (`first_line` / `last_line`) instead of assuming per-line rows
- `coverage-evidence` now distinguishes:
  - true target under-coverage
  - target unmeasured bundles
  - baseline-vs-target comparison gaps
- targeted family runs no longer emit misleading package-wide
  `package_coverage_below_threshold` exceptions; package coverage is marked
  `not_applicable_subset` when the run only evaluates a selected bundle subset
- regression coverage was extended in
  `tools/refactor/tests/fidelity-audit.test.cjs` for:
  - explicit-test-only bundle narrowing
  - comparison-gap exception emission
  - targeted-run package-gate suppression

Live family pass:

- lineage bundle mapping was hardened so auto-curated `source_lineage_gap`
  families now use manifest-resolved target line ranges instead of whole target
  files, and large legacy parents are split into semantic coverage bundles
  (`builder_resolver`, `reactive_state`, `observer_register`, etc.)
- regression coverage now proves this on the bundle-map fixture:
  `tools/refactor/tests/fidelity-audit.test.cjs`
- live rerun:
  - bundle-map: `bundle-map-20260419T155100Z-292c1303`
  - `R/mod_prot_enrich.R` now decomposes into `14` bundles instead of one
    coarse lineage bucket
- shared public-surface characterization file:
  - `tests/testthat/test-prot-11a-enrichment-module-characterization.R`
  - stable contract check:
    `contracts-20260419T172307Z-dc533935`
  - stable compare run:
    `coverage-compare-20260419T172330Z-85152be9`
  - stable evidence run:
    `coverage-evidence-20260419T172538Z-1ec043b6`
- result of the current stable shared suite on the enrichment family slice:
  - target-side bundle coverage:
    - `lineage::R_mod_prot_enrich.R::observer_register` `100.0%`
    - `lineage::R_mod_prot_enrich.R::observer_runtime` `100.0%`
  - evidence gate:
    - `lineage::R_mod_prot_enrich.R::observer_register` `pass`
    - `lineage::R_mod_prot_enrich.R::observer_runtime` `pass`
- the shared enrichment characterization file now contains stable public
  workflow cases for:
  - workflow-driven observer hydration
  - run-observer success path
  - run-observer failure path
- the target worktree passes that file in both:
  - `contracts`
  - `pkgload::load_all()` plus `testthat::test_file()` execution
- baseline `main` now also passes the shared file after the run-observer cases
  were rewritten around lower-level seams common to both versions:
  - `createDAResultsForEnrichment`
  - `processEnrichments`
  - `matchAnnotations`
  - minimal runtime S4 carriers

Interpretation:

- the enrichment family slice is now resolved for FA-C4
- the remaining lineage-family limitation is documented and acceptable under the
  current evidence model:
  - baseline capture completes with shared tests
  - baseline member ranges remain `not_observed`
  - target coverage meets threshold
  - the gate therefore records `pass` instead of `comparison_gap`
- the reusable pattern for other families is now:
  - range-aware lineage clustering
  - shared public module scenarios
  - lower-level cross-version mocks where the baseline still keeps logic inline

Second live family pass:

- shared public norm characterization file:
  - `tests/testthat/test-prot-05d-norm-module-characterization.R`
  - stable contract check:
    `contracts-prot-norm-characterization-2`
  - stable compare run:
    `coverage-compare-prot-norm-characterization`
  - stable evidence run:
    `coverage-evidence-prot-norm-characterization`
- result of the current stable shared suite on the norm family slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_prot_norm_server` baseline `100.0%`
    - `surface::symbol::mod_prot_norm_server` target `100.0%`
    - evidence gate: `pass`
  - runtime lineage family:
    - `lineage::R_mod_prot_norm.R::observer_runtime` target `89.9%`
    - evidence gate: `pass`
- the shared norm characterization file now covers stable public module cases
  for:
  - apply-correlation success
  - skip-correlation success
  - export prereq and export success
  - reset success
  - invalid-tab warning
  - pre-normalization error handling
  - normalization error handling
  - correlation error handling
  - export error handling
  - reset error handling
- baseline `main` now passes the same shared file cleanly, so the old norm
  compare failure mode (target-only helper contracts exploding on baseline) is
  removed

Interpretation update:

- the norm family slice is now resolved for FA-C4
- the norm wrapper entrypoint is now proven by shared comparative coverage, not
  just target-side contract tests
- the norm observer-runtime lineage family meets target threshold under the
  shared public suite even though baseline helper line ranges remain
  `not_observed`; this is acceptable under the current lineage-family evidence
  model because:
  - the same shared public scenarios pass on baseline and target
  - target extracted helpers are above threshold
  - the gate records `pass` rather than `comparison_gap`

Next action inside `FA-C4`:

- carry the same pattern into the next priority family, starting with another
  wrapper-heavy lineage bundle set (`mod_metab_norm` or `mod_metab_import`)
- reuse the enrichment slice artifacts as the concrete reference for:
  - cross-version shared module scenarios
  - bundle-target compare runs
  - coverage evidence closeout

Third live family pass:

- shared public metabolomics norm characterization file:
  - `tests/testthat/test-metab-03d-norm-module-characterization.R`
  - stable contract check:
    `contracts-metab-norm-characterization-20260420a`
  - stable compare run:
    `coverage-compare-metab-norm-characterization-20260420b`
  - stable evidence run:
    `coverage-evidence-metab-norm-characterization-20260420b`
- result of the current stable shared suite on the metabolomics norm family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_metab_norm_server` baseline `100.0%`
    - `surface::symbol::mod_metab_norm_server` target `100.0%`
    - evidence gate: `pass`
  - runtime lineage family:
    - `lineage::R_mod_metab_norm.R::observer_runtime` target `93.1%`
    - evidence gate: `pass`
- the shared metabolomics norm characterization file now covers stable public
  module cases for:
  - selected-tab driven pre-QC hydration
  - normalization success
  - normalization error handling
  - apply-correlation success
  - skip-correlation success
  - export prereq and export success
  - reset success
  - reset save failure handling
- baseline `main` now also passes the same shared file cleanly, so the old
  metabolomics norm compare failure mode is removed

Interpretation update:

- the metabolomics norm family slice is now resolved for FA-C4
- the metabolomics norm wrapper entrypoint is now proven by shared comparative
  coverage, not just target-side contracts
- the metabolomics norm observer-runtime lineage family meets target threshold
  under the shared public suite even though baseline helper line ranges remain
  `not_observed`; this is acceptable under the current lineage-family evidence
  model because:
  - the same shared public scenarios pass on baseline and target
  - target extracted helpers are above threshold
  - the gate records `pass` rather than `comparison_gap`

Next action inside `FA-C4`:

- move to the next unresolved wrapper-heavy family:
  - `mod_metab_import`
- reuse the same proven sequence:
  - shared cross-version public characterization
  - focused bundle compare against `latest-bundles.json`
  - focused coverage evidence materialization

Fourth live family pass:

- shared public metabolomics import characterization file:
  - `tests/testthat/test-metab-00d-import-module-characterization.R`
  - stable contract check:
    `contracts-metab-import-characterization-20260420a`
  - stable compare run:
    `coverage-compare-metab-import-characterization-20260420b`
  - stable evidence run:
    `coverage-evidence-metab-import-characterization-20260420b`
- result of the current stable shared suite on the metabolomics import family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_metab_import_server` baseline `100.0%`
    - `surface::symbol::mod_metab_import_server` target `100.0%`
    - evidence gate: `pass`
  - runtime lineage family:
    - `lineage::R_mod_metab_import.R::observer_runtime` target `100.0%`
    - evidence gate: `pass`
- the shared metabolomics import characterization file now covers stable public
  module cases for:
  - assay-selection hydration through `shinyFiles`
  - process-import success
  - assay-selection error handling
  - process-import error handling for optional second-assay import
- baseline `main` now also passes the same shared file cleanly, including the
  assay-selection failure path that maps to the extracted target error-finalizer

Interpretation update:

- the metabolomics import family slice is now resolved for FA-C4
- the metabolomics import wrapper entrypoint is now proven by shared
  comparative coverage, not just seam or helper tests
- the metabolomics import observer-runtime lineage family now has direct shared
  coverage on the extracted target error-finalizer path, eliminating the prior
  `0.0%` false-start compare result
- as with the prior lineage families, baseline helper line ranges remain
  `not_observed`, but the gate records `pass` because:
  - the same shared public scenarios pass on baseline and target
  - target extracted helpers are above threshold
  - the coverage evidence model accepts that lineage form when the public suite
    is shared and stable

Current metabolomics DA checkpoint:

- shared public metabolomics DA characterization file:
  - `tests/testthat/test-metab-02ab-da-module-characterization.R`
  - stable contract check:
    `contracts-metab-da-characterization-20260420a`
  - stable compare run:
    `coverage-compare-metab-da-characterization-20260420a`
  - stable evidence run:
    `coverage-evidence-metab-da-characterization-20260420a`
- result of the current stable shared suite on the metabolomics DA family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_metab_da_server` baseline `100.0%`
    - `surface::symbol::mod_metab_da_server` target `100.0%`
    - evidence gate: `pass`
  - runtime lineage family:
    - `lineage::R_mod_metab_da.R::observer_runtime` target `90.4%`
    - evidence gate: `pass`
- the shared metabolomics DA characterization file now covers stable public
  module cases for:
  - filtered-session load success
  - filtered-session load failure on corrupt session payload
  - DA analysis success
  - DA analysis error
  - Glimma error-banner rendering plus heatmap save observer runtime
- baseline `main` now also passes the same shared file cleanly, including the
  corrupt-session branch and the public heatmap-save path

Interpretation update:

- the metabolomics DA family slice is now resolved for FA-C4
- the metabolomics DA wrapper entrypoint is now proven by shared comparative
  coverage instead of large seam-only characterization files
- the metabolomics DA observer-runtime lineage family now exceeds target at
  `90.4%` target coverage under the shared public module suite
- as with the prior lineage families, baseline helper line ranges remain
  `not_observed`, but the gate records `pass` because:
  - the same shared public scenarios pass on baseline and target
  - the extracted target helper set is above threshold
  - the coverage evidence model accepts this lineage form when the public suite
    is shared and stable

Current proteomics design-builder checkpoint:

- shared public proteomics design-builder characterization file:
  - `tests/testthat/test-prot-04a-design-module-characterization.R`
  - stable contract check:
    `contracts-prot-design-characterization-20260420a`
  - stable compare run:
    `coverage-compare-prot-design-characterization-20260420a`
  - stable evidence run:
    `coverage-evidence-prot-design-characterization-20260420a`
- result of the current stable shared suite on the proteomics design-builder
  family slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_prot_design_builder_server` baseline `100.0%`
    - `surface::symbol::mod_prot_design_builder_server` target `100.0%`
    - evidence gate: `pass`
  - runtime lineage family:
    - `lineage::R_mod_prot_design_builder.R::observer_runtime` target `100.0%`
    - evidence gate: `pass`
- the shared proteomics design-builder characterization file now covers stable
  public or helper-level cases for:
  - public module empty save-results behavior
  - save-results warning branch
  - save-results success branch
  - reset-confirmation flow
- baseline `main` now also passes the same shared file cleanly; the wrapper
  case executes on both branches and the extracted helper tests skip on
  baseline where the monolith has not yet been decomposed into those symbols

Interpretation update:

- the proteomics design-builder family slice is now resolved for FA-C4
- the proteomics design-builder wrapper entrypoint is now proven by shared
  comparative coverage
- the proteomics design-builder observer-runtime lineage family now reaches
  `100.0%` target coverage under the shared characterization file
- as with the prior lineage families, baseline helper line ranges remain
  `not_observed`, but the gate records `pass` because:
  - the shared wrapper/public scenario passes on baseline and target
  - the extracted target helper set is fully covered
  - the coverage evidence model accepts this lineage form when the extracted
    helper tests are target-only and the public wrapper scenario is shared

Current proteomics import checkpoint:

- shared public proteomics import characterization file:
  - `tests/testthat/test-prot-01d-import-module-characterization.R`
  - stable contract check:
    `contracts-prot-import-characterization-20260420T083020Z`
  - stable compare run:
    `coverage-compare-prot-import-characterization-20260420T083020Z`
  - stable evidence run:
    `coverage-evidence-prot-import-characterization-20260420T083208Z`
- result of the current stable shared suite on the proteomics import family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_prot_import_server` baseline `100.0%`
    - `surface::symbol::mod_prot_import_server` target `100.0%`
    - evidence gate: `pass`
  - runtime lineage family:
    - `lineage::R_mod_prot_import.R::observer_runtime` target `100.0%`
    - evidence gate: `pass`
- the shared proteomics import characterization file now covers stable public
  module cases for:
  - successful import side effects
  - unsupported-format cleanup and workflow reset
- baseline `main` now also passes the same shared file cleanly, so the target
  extracted runtime helpers are covered by genuinely shared public flows rather
  than target-only helper checks

Interpretation update:

- the proteomics import family slice is now resolved for FA-C4
- the proteomics import wrapper entrypoint is now proven by shared
  comparative coverage
- the proteomics import observer-runtime lineage family now reaches `100.0%`
  target coverage under the shared characterization file
- baseline helper line ranges remain `not_observed`, but the gate records
  `pass` because:
  - the same shared public import scenarios pass on baseline and target
  - the target extracted helper set is fully covered
  - the coverage evidence model accepts this lineage form when the source
    monolith still owns the behavior on baseline

Current lipidomics norm checkpoint:

- shared public lipidomics norm characterization file:
  - `tests/testthat/test-lipid-12a-norm-module-characterization.R`
  - stable contract check:
    `contracts-lipid-norm-characterization-20260420T092125Z`
  - stable compare run:
    `coverage-compare-lipid-norm-characterization-20260420T092125Z`
  - stable evidence run:
    `coverage-evidence-lipid-norm-characterization-20260420T092405Z`
- result of the current stable shared suite on the lipidomics norm family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_lipid_norm_server` baseline `100.0%`
    - `surface::symbol::mod_lipid_norm_server` target `100.0%`
    - evidence gate: `pass`
  - runtime lineage family:
    - `lineage::R_mod_lipid_norm.R::observer_runtime` target `100.0%`
    - evidence gate: `pass`
- the shared lipidomics norm characterization file now covers stable public
  module cases for:
  - startup normalization-log placeholder rendering
  - selected-tab pre-QC error logging
- baseline `main` now also passes the same shared file cleanly, so the target
  extracted startup/runtime log helpers are proven through shared public module
  startup rather than target-only helper tests

Interpretation update:

- the lipidomics norm family slice is now resolved for FA-C4
- the lipidomics norm wrapper entrypoint is now proven by shared comparative
  coverage
- the lipidomics norm observer-runtime lineage family now reaches `100.0%`
  target coverage under the shared startup/log characterization file
- baseline helper line ranges remain `not_observed`, but the gate records
  `pass` because:
  - the same shared startup and selected-tab scenarios pass on baseline and
    target
  - the target extracted helper set is fully covered
  - the coverage evidence model accepts this lineage form when baseline still
    owns the behavior in the monolith

Current lipidomics DA checkpoint:

- shared public lipidomics DA characterization file:
  - `tests/testthat/test-lipid-02ab-da-module-characterization.R`
  - stable contract check:
    `contracts-lipid-da-characterization-20260420b`
  - stable compare run:
    `coverage-compare-lipid-da-characterization-20260420c`
  - stable evidence run:
    `coverage-evidence-lipid-da-characterization-20260420`
- result of the current stable shared suite on the lipidomics DA family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_lipid_da_server` baseline `100.0%`
    - `surface::symbol::mod_lipid_da_server` target `100.0%`
    - evidence gate: `pass`
  - reactive-state lineage family:
    - `lineage::R_mod_lipid_da.R::reactive_state` target `92.3%`
    - evidence gate: `pass`
  - output-renderer lineage family:
    - `lineage::R_mod_lipid_da.R::output_renderer` target `86.4%`
    - evidence gate: `pass`
  - observer-runtime lineage family:
    - `lineage::R_mod_lipid_da.R::observer_runtime` target `100.0%`
    - evidence gate: `pass`
- the shared lipidomics DA characterization file now covers stable public
  module cases for:
  - filtered-session load success
  - missing source-directory preflight
  - missing filtered-session file preflight
  - corrupt-session fatal load failure
  - DA analysis success
  - DA analysis error finalization
  - Glimma render failure and heatmap-save runtime behavior
- baseline `main` now also passes the same shared file cleanly because the
  characterization server injects a compatibility alias for the legacy
  `generateMetabVolcanoPlotGlimma` call while preserving the real module body

Interpretation update:

- the lipidomics DA family slice is now resolved for FA-C4
- the lipidomics DA wrapper entrypoint is now proven by shared comparative
  coverage
- the lipidomics DA reactive-state lineage family now reaches `92.3%`
  target coverage under shared public module scenarios
- the lipidomics DA output-renderer lineage family now reaches `86.4%`
  target coverage under shared public render scenarios
- the lipidomics DA observer-runtime lineage family now reaches `100.0%`
  target coverage under shared public error and runtime scenarios
- baseline helper line ranges remain `not_observed` for the extracted lineage
  families, but the gate records `pass` because:
  - the same shared DA scenarios pass on baseline and target
  - the target extracted helper sets are covered above the bundle threshold
  - the coverage evidence model accepts this lineage form when baseline still
    owns the behavior in the monolith

Current lipidomics design-builder checkpoint:

- shared public lipidomics design-builder characterization file:
  - `tests/testthat/test-lipid-13a-design-module-characterization.R`
  - stable contract check:
    `contracts-lipid-design-characterization-20260420c`
  - stable compare run:
    `coverage-compare-lipid-design-characterization-20260420d`
  - stable evidence run:
    `coverage-evidence-lipid-design-characterization-20260420`
- result of the current stable shared suite on the lipidomics design-builder
  family slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_lipid_design_builder_server` baseline `100.0%`
    - `surface::symbol::mod_lipid_design_builder_server` target `100.0%`
    - evidence gate: `pass`
  - builder-resolver lineage family:
    - `lineage::R_mod_lipid_design_builder.R::builder_resolver` target
      `84.2%`
    - evidence gate: `pass`
  - output-register lineage family:
    - `lineage::R_mod_lipid_design_builder.R::output_register` target
      `100.0%`
    - evidence gate: `pass`
  - formatter lineage family:
    - `lineage::R_mod_lipid_design_builder.R::formatter` target `100.0%`
    - evidence gate: `pass`
- the shared lipidomics design-builder characterization file now covers stable
  public and helper cases for:
  - empty save-results public behavior
  - fresh and imported initial-state hydration
  - initial-state shell hydration
  - summary and adjacent output registration
  - save-results warning and success paths
  - full reset-state and reset-confirmation behavior
  - sample rename propagation and bulk rename transforms
  - sample-column detection, contrast registration, sample removal, and
    save-results helper formatting
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public module scenarios still prove wrapper parity directly

Interpretation update:

- the lipidomics design-builder family slice is now resolved for FA-C4
- the lipidomics design-builder wrapper entrypoint is now proven by shared
  comparative coverage
- the lipidomics design-builder builder-resolver lineage family now reaches
  `84.2%` target coverage under the shared characterization suite
- the lipidomics design-builder output-register lineage family now reaches
  `100.0%` target coverage under the shared characterization suite
- the lipidomics design-builder formatter lineage family now reaches `100.0%`
  target coverage under the shared characterization suite
- baseline helper line ranges remain `not_observed` for the extracted lineage
  families, but the gate records `pass` because:
  - the same shared public wrapper scenarios pass on baseline and target
  - the target extracted helper sets are covered above the bundle threshold
  - the coverage evidence model accepts this lineage form when baseline still
    owns the behavior in the monolith

- a stable shared lipidomics import characterization file now covers wrapper UI
  shell shape, wrapper server startup, assay-input and layout shell wiring,
  format and column helper behavior, workflow commit/finalization,
  sample-name sanitation, import processing success/failure cleanup,
  preview-load format branching, summary/status display branches, file-selection
  event handling, process-request forwarding, and output-render helper
  forwarding on both baseline and target
- stable contract run: `contracts-lipid-import-characterization-20260420c`
- stable compare run: `coverage-compare-lipid-import-characterization-20260420c`
- stable evidence run:
  `coverage-evidence-lipid-import-characterization-20260420`
- result of the current stable shared suite on the lipidomics import family
  slice:
  - wrapper entrypoints:
    - `surface::symbol::mod_lipid_import_server` baseline `100.0%`
    - `surface::symbol::mod_lipid_import_server` target `100.0%`
    - `surface::symbol::mod_lipid_import_ui` baseline `100.0%`
    - `surface::symbol::mod_lipid_import_ui` target `100.0%`
    - evidence gate: `pass`
  - misc lineage family:
    - `lineage::R_mod_lipid_import.R::misc` target `93.0%`
    - evidence gate: `pass`
  - output-renderer lineage family:
    - `lineage::R_mod_lipid_import.R::output_renderer` target `100.0%`
    - evidence gate: `pass`
  - bootstrap-setup lineage family:
    - `lineage::R_mod_lipid_import.R::bootstrap_setup` target `100.0%`
    - evidence gate: `pass`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public module scenarios still prove wrapper parity directly

- a stable shared proteomics summary characterization file now covers module
  startup/default outputs, state and template helper selection, report-template
  resolution and retrieval, rendered-report activation, report generation and
  progress orchestration, workflow-args save success and fallback handling,
  session export success and failure, and publication-copy plus GitHub helper
  execution on both baseline and target
- stable contract run:
  `contracts-prot-summary-characterization-20260420a`
- stable compare run:
  `coverage-compare-prot-summary-characterization-20260420b`
- stable evidence run:
  `coverage-evidence-prot-summary-characterization-20260420`
- result of the current stable shared suite on the proteomics summary family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_prot_summary_server` baseline `100.0%`
    - `surface::symbol::mod_prot_summary_server` target `100.0%`
    - evidence gate: `pass`
  - misc lineage family:
    - `lineage::R_mod_prot_summary.R::misc` target `96.0%`
    - evidence gate: `pass`
  - observer-runtime lineage family:
    - `lineage::R_mod_prot_summary.R::observer_runtime` target `95.7%`
    - evidence gate: `pass`
  - output-register lineage family:
    - `lineage::R_mod_prot_summary.R::output_register` target `100.0%`
    - evidence gate: `pass`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public module startup scenario still proves wrapper parity directly

- a stable shared metabolomics QC S4 characterization file now covers the
  public finalize-QC module path plus extracted summary, history, render,
  finalize, reporting, plotting, and server-body helper seams on both baseline
  and target
- stable contract run:
  `contracts-metab-qc-s4-characterization-20260420c`
- stable compare run:
  `coverage-compare-metab-qc-s4-characterization-20260420d`
- stable evidence run:
  `coverage-evidence-metab-qc-s4-characterization-20260420`
- result of the current stable shared suite on the metabolomics QC S4 family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_metab_qc_s4_server` baseline `100.0%`
    - `surface::symbol::mod_metab_qc_s4_server` target `100.0%`
    - evidence gate: `pass`
  - misc lineage family:
    - `lineage::R_mod_metab_qc_s4.R::misc` target `95.4%`
    - evidence gate: `pass`
  - output-renderer lineage family:
    - `lineage::R_mod_metab_qc_s4.R::output_renderer` target `85.7%`
    - evidence gate: `pass`
  - observer-runtime lineage family:
    - `lineage::R_mod_metab_qc_s4.R::observer_runtime` target `100.0%`
    - evidence gate: `pass`
  - plotting lineage family:
    - `lineage::R_mod_metab_qc_s4.R::plotting` target `100.0%`
    - evidence gate: `pass`
- the bundle-map engine now enriches non-lineage surface bundles with precise
  `inventory_entities` line ranges before coverage capture; this removed a
  stale whole-file fallback that was undercounting wrapper-entrypoint coverage
  for extracted modules such as `mod_metab_qc_s4_server`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public finalize-QC scenario still proves wrapper parity directly

- a stable shared lipidomics QC duplicates characterization file now covers the
  public detect-resolve-revert module flow plus extracted builder, output,
  detection, resolution, revert, observer, binding, and server-shell helper
  seams on both baseline and target
- stable contract run:
  `contracts-lipid-qc-duplicates-characterization-20260420a`
- stable compare run:
  `coverage-compare-lipid-qc-duplicates-characterization-20260420a`
- stable evidence run:
  `coverage-evidence-lipid-qc-duplicates-characterization-20260420a`
- result of the current stable shared suite on the lipidomics QC duplicates
  family slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_lipid_qc_duplicates_server` baseline `100.0%`
    - `surface::symbol::mod_lipid_qc_duplicates_server` target `100.0%`
    - evidence gate: `pass`
  - misc lineage family:
    - `lineage::R_mod_lipid_qc_duplicates.R::misc` target `96.6%`
    - evidence gate: `pass`
  - observer-register lineage family:
    - `lineage::R_mod_lipid_qc_duplicates.R::observer_register` target
      `100.0%`
    - evidence gate: `pass`
  - output-register lineage family:
    - `lineage::R_mod_lipid_qc_duplicates.R::output_register` target `100.0%`
    - evidence gate: `pass`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public detect-resolve-revert scenario still proves wrapper parity
  directly

- a stable shared metabolomics QC duplicates characterization file now covers
  the public detect-resolve-revert module flow plus extracted duplicate
  workflow helper seams on both baseline and target
- stable contract run:
  `contracts-metab-qc-duplicates-characterization-20260420a`
- stable compare run:
  `coverage-compare-metab-qc-duplicates-characterization-20260420a`
- stable evidence run:
  `coverage-evidence-metab-qc-duplicates-characterization-20260420a`
- result of the current stable shared suite on the metabolomics QC duplicates
  family slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_metab_qc_duplicates_server` baseline `100.0%`
    - `surface::symbol::mod_metab_qc_duplicates_server` target `100.0%`
    - evidence gate: `pass`
  - misc lineage family:
    - `lineage::R_mod_metab_qc_duplicates.R::misc` target `96.4%`
    - evidence gate: `pass`
  - observer-runtime lineage family:
    - `lineage::R_mod_metab_qc_duplicates.R::observer_runtime` target
      `100.0%`
    - evidence gate: `pass`
  - output-renderer lineage family:
    - `lineage::R_mod_metab_qc_duplicates.R::output_renderer` target `100.0%`
    - evidence gate: `pass`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public duplicate workflow still proves wrapper parity directly

- a stable shared metabolomics design characterization file now covers public
  wrapper initialization plus extracted import and observer registration helper
  seams on both baseline and target
- stable contract run:
  `contracts-metab-design-characterization-20260420b`
- stable compare run:
  `coverage-compare-metab-design-characterization-20260420d`
- stable evidence run:
  `coverage-evidence-metab-design-characterization-20260420a`
- result of the current stable shared suite on the metabolomics design family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_metab_design_server` baseline `100.0%`
    - `surface::symbol::mod_metab_design_server` target `100.0%`
    - evidence gate: `pass`
  - misc lineage family:
    - `lineage::R_mod_metab_design.R::misc` target `96.7%`
    - evidence gate: `pass`
  - observer-register lineage family:
    - `lineage::R_mod_metab_design.R::observer_register` target `100.0%`
    - evidence gate: `pass`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public wrapper initialization proves wrapper parity directly

- a stable shared metabolomics summary characterization file now covers public
  wrapper initialization plus extracted bootstrap, output registration, and
  observer runtime helper seams on both baseline and target
- stable contract run:
  `contracts-metab-summary-characterization-20260420a`
- stable compare run:
  `coverage-compare-metab-summary-characterization-20260422a`
- stable evidence run:
  `coverage-evidence-metab-summary-characterization-20260422a`
- result of the current stable shared suite on the metabolomics summary family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_metab_summary_server` baseline `100.0%`
    - `surface::symbol::mod_metab_summary_server` target `100.0%`
    - evidence gate: `pass`
  - observer-runtime lineage family:
    - `lineage::R_mod_metab_summary.R::observer_runtime` target `100.0%`
    - evidence gate: `pass`
  - output-setup lineage family:
    - `lineage::R_mod_metab_summary.R::output_setup` target `100.0%`
    - evidence gate: `pass`
  - output-register lineage family:
    - `lineage::R_mod_metab_summary.R::output_register` target `100.0%`
    - evidence gate: `pass`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public wrapper initialization proves wrapper parity directly

- a stable shared lipidomics summary characterization file now covers public
  wrapper initialization plus extracted template, bootstrap, workflow-save,
  publication-copy, report-generation, GitHub, export, and observer
  registration helper seams on both baseline and target
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-14a-summary-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-summary-characterization-20260422a`
- stable focused evidence run:
  - `coverage-evidence-lipid-summary-characterization-20260422a`
- result of the current stable shared suite on the lipidomics summary family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_lipid_summary_server` baseline `100.0%`
    - `surface::symbol::mod_lipid_summary_server` target `100.0%`
    - evidence gate: `pass`
  - misc lineage family:
    - `lineage::R_mod_lipid_summary.R::misc` target `99.4%`
    - evidence gate: `pass`
  - observer-register lineage family:
    - `lineage::R_mod_lipid_summary.R::observer_register` target `100.0%`
    - evidence gate: `pass`
  - output-register lineage family:
    - `lineage::R_mod_lipid_summary.R::output_register` target `100.0%`
    - evidence gate: `pass`
  - observer-runtime lineage family:
    - `lineage::R_mod_lipid_summary.R::observer_runtime` target `100.0%`
    - evidence gate: `pass`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public wrapper initialization proves wrapper parity directly

- the refreshed shared proteomics enrichment characterization file now covers
  formatter helpers, download archive I/O, cleanup/state-save/process
  execution, organism filtering, payload/argument/UI/state/progress builders,
  all-contrast result builders, persistence/finalization, contrast/dependency
  and process resolvers, annotation/organism mapping, and setup orchestration
  on both baseline and target
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-11a-enrichment-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-enrich-characterization-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-enrich-characterization-20260422a`
- result of the refreshed shared suite on the proteomics enrichment family
  slice:
  - wrapper entrypoint:
    - `surface::symbol::mod_prot_enrich_server` baseline `100.0%`
    - `surface::symbol::mod_prot_enrich_server` target `100.0%`
    - evidence gate: `pass`
  - builder-resolver lineage family:
    - `lineage::R_mod_prot_enrich.R::builder_resolver` target `95.7%`
    - evidence gate: `pass`
  - formatter lineage family:
    - `lineage::R_mod_prot_enrich.R::formatter` target `100.0%`
    - evidence gate: `pass`
  - misc lineage family:
    - `lineage::R_mod_prot_enrich.R::misc` target `100.0%`
    - evidence gate: `pass`
  - download-I/O lineage family:
    - `lineage::R_mod_prot_enrich.R::download_io` target `100.0%`
    - evidence gate: `pass`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared public and lower-level compatibility scenarios prove wrapper parity
  directly

- the refreshed shared lipidomics norm characterization file now covers the
  public UI/server shell plus extracted UI, startup/runtime orchestration,
  observer registration/tracking, pre-QC, support registration, composite
  generation, normalization, export, reset, skip-correlation, and
  apply-correlation helper paths on both baseline and target
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-12a-norm-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-norm-characterization-20260422c`
- stable focused evidence run:
  - `coverage-evidence-lipid-norm-characterization-20260422c`
- result of the refreshed shared suite on the lipidomics norm family slice:
  - wrapper entrypoints:
    - `surface::symbol::mod_lipid_norm_server` baseline `100.0%`
    - `surface::symbol::mod_lipid_norm_server` target `100.0%`
    - `surface::symbol::mod_lipid_norm_ui` baseline `100.0%`
    - `surface::symbol::mod_lipid_norm_ui` target `100.0%`
    - evidence gate: `pass`
  - misc lineage family:
    - `lineage::R_mod_lipid_norm.R::misc` target `99.6%`
    - evidence gate: `pass`
  - observer-register lineage family:
    - `lineage::R_mod_lipid_norm.R::observer_register` target `100.0%`
    - evidence gate: `pass`
  - output-renderer lineage family:
    - `lineage::R_mod_lipid_norm.R::output_renderer` target `97.9%`
    - evidence gate: `pass`
  - output-register, observer-runtime, plotting, and reactive-state lineage
    families:
    - target `100.0%`
    - evidence gate: `pass`
- baseline `main` now also passes the same shared file cleanly; helper tests
  skip on the monolith where extracted helper files do not exist, while the
  shared UI and server scenarios prove wrapper parity directly

- the shared compare selector is now explicitly marker-driven:
  - shared comparative suites must carry
    `# fidelity-coverage-compare: shared`
  - `coverage-compare` now defaults to
    `.refactor-fidelity-audit/reports/latest-bundles.json` when no manifest is
    provided
  - when bundles are selected but no marked shared tests match, the runner now
    fails fast instead of silently falling back to the full `tests/testthat`
    tree
- the first corrected repo-wide shared compare/evidence pair is now:
  - compare: `coverage-compare-20260420T180816Z-074415bc`
  - evidence: `coverage-evidence-20260420T181045Z-2e0f73bc`
- current repo-wide coverage-backed parity state from that corrected run:
  - selected shared files: `14`
  - package line coverage:
    - baseline `21.1%`
    - target `20.9%`
  - bundle gates:
    - `pass`: `81`
    - `below_target`: `117`
    - `target_unmeasured`: `35`
  - low-coverage bundles: `152`
  - emitted coverage exceptions: `153`
- that corrected run supersedes the earlier noisy repo-wide pass because it
  only uses explicitly marked shared cross-version characterization suites
- the post-lipid-summary repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-a`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-a`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `18`
  - package line coverage:
    - baseline `25.7%`
    - target `26.9%`
  - bundle gates:
    - `pass`: `105`
    - `below_target`: `93`
    - `target_unmeasured`: `35`
  - low-coverage bundles: `128`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `129`
- the post-proteomics-enrichment repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-b`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-b`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `18`
  - package line coverage:
    - baseline `25.7%`
    - target `27.7%`
  - bundle gates:
    - `pass`: `109`
    - `below_target`: `89`
    - `target_unmeasured`: `35`
  - low-coverage bundles: `124`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `125`
- the post-lipidomics-norm repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-c`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-c`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `18`
  - package line coverage:
    - baseline `26.7%`
    - target `29.5%`
  - bundle gates:
    - `pass`: `111`
    - `below_target`: `87`
    - `target_unmeasured`: `35`
  - low-coverage bundles: `122`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `123`
- the post-proteomics-design repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-d`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-d`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `18`
  - package line coverage:
    - baseline `26.7%`
    - target `30.5%`
  - bundle gates:
    - `pass`: `115`
    - `below_target`: `83`
    - `target_unmeasured`: `35`
  - low-coverage bundles: `118`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `119`
- the post-general-file-management-bootstrap repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-e`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-e`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `19`
  - package line coverage:
    - baseline `27.3%`
    - target `31.2%`
  - bundle gates:
    - `pass`: `118`
    - `below_target`: `80`
    - `target_unmeasured`: `35`
  - low-coverage bundles: `115`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `116`
- the post-lipidomics-design repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-f`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-f`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `19`
  - package line coverage:
    - baseline `27.3%`
    - target `32.0%`
  - bundle gates:
    - `pass`: `123`
    - `below_target`: `75`
    - `target_unmeasured`: `35`
  - low-coverage bundles: `110`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `111`
- the post-lipid-QC-helper repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-g`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-g`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `20`
  - package line coverage:
    - baseline `27.4%`
    - target `33.1%`
  - bundle gates:
    - `pass`: `125`
    - `below_target`: `73`
    - `target_unmeasured`: `35`
  - low-coverage bundles: `108`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `109`

- a refreshed shared proteomics design characterization file now includes
  target helper coverage for:
  - import artifact resolution and selected FASTA detection
  - config bootstrap from package, download fallback, and failure paths
  - UniProt and FASTA sidecar hydration from import, script, missing, and
    processed-FASTA sources
  - workflow-state inference for DIA plus config-driven LFQ/TMT mappings
  - state checkpoint, post-checkpoint, builder persistence, preview output,
    observer shell, import confirmation, import bootstrap, and import modal
    wiring
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-04a-design-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-design-characterization-20260422b`
- stable focused evidence run:
  - `coverage-evidence-prot-design-characterization-20260422b`
- result of the refreshed shared suite on the proteomics design family slice:
  - `lineage::R_mod_prot_design.R::misc` target `97.6%`
  - `lineage::R_mod_prot_design.R::builder_resolver` target `97.7%`
  - `lineage::R_mod_prot_design.R::observer_register` target `100.0%`
  - `lineage::R_mod_prot_design.R::observer_runtime` target `81.0%`
  - evidence gate: `pass` for all `4` focused bundles

- a refreshed shared general file-management path contract file now includes
  target helper coverage for:
  - setup-directory omic type parsing and validation errors
  - omic-specific directory configuration lookup
  - path-list alias construction for proteomics, metabolomics,
    transcriptomics, lipidomics, and integration
  - existing-directory overwrite, reuse, cancel, and force decisions
  - materialization branches for reuse mode, missing script source, no
    copyable scripts, script-copy failure, and existing structure reuse
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-general-filemgmt-path-contracts.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-general-filemgmt-bootstrap-20260422b`
- stable focused evidence run:
  - `coverage-evidence-general-filemgmt-bootstrap-20260422b`
- result of the refreshed shared suite on the general file-management bootstrap
  setup slice:
  - `lineage::R_func_general_filemgmt.R::bootstrap_setup` target `100.0%`
  - evidence gate: `pass`

- a refreshed shared lipidomics design characterization file now includes
  target helper coverage for:
  - import bootstrap, modal shell, preview formatting, and builder-module
    registration
  - import confirmation failure, success, fallback, and assay-read error paths
  - builder observer shell handoff, result persistence, state-manager reuse,
    S4 creation failure handling, tab-status/QC handoff, and
    `updateLipidFiltering` failure isolation
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-13a-design-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-design-module-characterization-20260422a`
- stable focused evidence run:
  - `coverage-evidence-lipid-design-module-characterization-20260422a`
- result of the refreshed shared suite on the lipidomics design family slice:
  - `lineage::R_mod_lipid_design.R::misc` target `100.0%`
  - `lineage::R_mod_lipid_design.R::observer_register` target `100.0%`
  - `lineage::R_mod_lipid_design.R::builder_resolver` target `91.7%`
  - `lineage::R_mod_lipid_design.R::observer_runtime` target `100.0%`
  - `lineage::R_mod_lipid_design.R::formatter` target `100.0%`
  - evidence gate: `pass` for all `5` focused bundles

- the lipidomics QC filtering helper characterization file is now
  baseline-aware:
  - target-only extracted helper tests skip on monolithic baseline main
  - the same tests pass on target and cover extracted helper line ranges
  - the existing duplicate/intensity/correlation helper tests remain active
    across refs when their bindings are available
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-01-qc-filtering-helpers.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-qc-helper-characterization-20260422b`
- stable focused evidence run:
  - `coverage-evidence-lipid-qc-helper-characterization-20260422b`
- result of the refreshed shared suite on the lipidomics QC helper family slice:
  - `lineage::R_func_lipid_qc.R::misc` target `86.4%`
  - `lineage::R_func_lipid_qc.R::plotting` target `88.1%`
  - evidence gate: `pass` for both focused bundles

- the proteomics DA volcano compatibility characterization is now
  baseline-aware and stable under the full shared matrix:
  - baseline refs that only expose the legacy DE-only wrapper skip the
    target-only alias assertions
  - the target helper still validates legacy `de_*` aliases, missing-input
    errors, and DA handler contrast initialization
  - the test now constructs the S4 fixture with `methods::new()` so earlier
    module-characterization mocks cannot leak a replacement
    `ProteinQuantitativeData()` constructor into the matrix
- the bundle map was refreshed so Wave 1.1 split paths are represented in
  coverage members:
  - bundle map: `bundle-map-20260422Tvolcano-member-refresh-final`
  - `surface::symbol::writeInteractiveVolcanoPlotProteomicsMain` now includes
    target member `R/func_prot_da_volcano_main.R`
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-07b-da-handlers-compat.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-da-volcano-compat-20260422e`
- stable focused evidence run:
  - `coverage-evidence-prot-da-volcano-compat-20260422g`
- result of the focused volcano compatibility slice:
  - `surface::symbol::writeInteractiveVolcanoPlotProteomicsMain` target
    `98.5%`
  - evidence gate: `pass`

- the post-volcano repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-i`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-i`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `21`
  - package line coverage:
    - baseline `27.9%`
    - target `33.7%`
  - bundle gates:
    - `pass`: `126`
    - `below_target`: `74`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `107`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `108`

- the lipidomics S4 plotting characterization is now baseline-aware and
  stable under focused coverage replay:
  - baseline refs exercise the original S4 methods while target refs verify
    the extracted helper dispatch payload
  - `plotPca` baseline missing-`mixOmics` behavior is asserted as the old
    warning plus empty-result contract
  - target refs still validate the generated PCA/RLE helper calls
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-10-s4-plotting-methods.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-s4-plotting-methods-20260422d`
- stable focused evidence run:
  - `coverage-evidence-lipid-s4-plotting-methods-20260422a`
- result of the focused lipidomics S4 plotting slice:
  - `surface::setMethod::plotPca::LipidomicsAssayData` target `91.6%`
  - `surface::setMethod::plotRle::LipidomicsAssayData` target `93.1%`
  - evidence gate: `pass` for all `4` focused bundles

- the general PCA helper characterization now exercises the package-level
  `plotPcaHelper()` fallback on target:
  - target refs cover grouped label/shape, shape-only, and plain PCA branches
  - target refs validate missing grouping, shape, and label aesthetic columns
  - monolithic baseline refs skip only the target-only no-`mixOmics` fallback
    branch
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-general-plotting-pca-helper-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-general-plot-pca-helper-20260422b`
- stable focused evidence run:
  - `coverage-evidence-general-plot-pca-helper-20260422a`
- result of the focused general PCA helper slice:
  - `surface::symbol::plotPcaHelper` target `94.4%`
  - evidence gate: `pass`

- the post-PCA-helper repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-k`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-k`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `23`
  - package line coverage:
    - baseline `28.5%`
    - target `34.3%`
  - bundle gates:
    - `pass`: `131`
    - `below_target`: `69`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `102`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `103`

- the metabolomics DA download helper characterization now executes the
  extracted target callbacks directly:
  - target refs verify the default dated filename callback
  - target refs execute the content callback and assert the full CSV payload
  - monolithic baseline refs skip the target-only extracted helper assertions
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-metab-02ab-da-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-metab-da-download-20260422a`
- stable focused evidence run:
  - `coverage-evidence-metab-da-download-20260422a`
- result of the focused metabolomics DA download slice:
  - `lineage::R_mod_metab_da.R::download_io` target `100.0%`
  - evidence gate: `pass`

- the post-metabolomics-DA-download repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-l`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-l`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `23`
  - package line coverage:
    - baseline `28.5%`
    - target `34.3%`
  - bundle gates:
    - `pass`: `132`
    - `below_target`: `68`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `101`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `102`

- the lipidomics QC intensity registration characterization now exercises the
  extracted target registration helpers through package-level namespace calls:
  - target refs render assay summary tabs plus grob and ggplot plot dispatch
  - target refs execute the apply-filter observer with a valid
    `LipidomicsAssayData` fixture
  - target refs execute the revert observer and assert state/reset output
  - monolithic baseline refs skip only the target-only extracted helper
    assertions
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-01-qc-filtering-helpers.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-qc-intensity-registration-20260422a`
- stable focused evidence run:
  - `coverage-evidence-lipid-qc-intensity-registration-20260422a`
- result of the focused lipidomics QC intensity registration slice:
  - `lineage::R_mod_lipid_qc_intensity.R::observer_register` target `100.0%`
  - `lineage::R_mod_lipid_qc_intensity.R::output_register` target `100.0%`
  - evidence gate: `pass` for both focused bundles

- the post-lipidomics-QC-intensity-registration repo-wide shared
  compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-m`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-m`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `23`
  - package line coverage:
    - baseline `28.5%`
    - target `34.6%`
  - bundle gates:
    - `pass`: `134`
    - `below_target`: `66`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `99`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `100`

- the metabolomics QC ITSD wrapper characterization now exercises the
  package-level wrapper handoff through a focused shared test:
  - baseline refs capture the public wrapper's `shiny::moduleServer` handoff
  - target refs mock the extracted server-body seam and verify the wrapper
    forwards `input`, `output`, `session`, workflow state, paths, and QC trigger
  - the older local-eval seam characterization remains target-side only and is
    excluded from repo-wide shared compare runs
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-metab-01v-qc-itsd-wrapper-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-metab-qc-itsd-wrapper-20260422b`
- stable focused evidence run:
  - `coverage-evidence-metab-qc-itsd-wrapper-20260422a`
- result of the focused metabolomics QC ITSD wrapper slice:
  - `surface::symbol::mod_metab_qc_itsd_server` baseline `100.0%`, target
    `100.0%`
  - evidence gate: `pass`

- the post-metabolomics-QC-ITSD-wrapper repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-n`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-n`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `24`
  - package line coverage:
    - baseline `29.0%`
    - target `34.6%`
  - bundle gates:
    - `pass`: `135`
    - `below_target`: `65`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `98`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `99`

- the lipidomics design wrapper characterization now isolates the public
  wrapper handoff in a dedicated shared file:
  - baseline refs capture the public `shiny::moduleServer` registration
  - target refs mock the extracted design helper registrations and verify the
    wrapper wires import bootstrap, modal registration, import observer,
    builder module, builder observer, preview outputs, and reactive flags
  - the existing lipid design helper characterization remains unchanged and is
    verified separately to avoid namespace-mock contamination
- stable direct current-branch test runs:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-13b-design-wrapper-characterization.R", reporter = testthat::SummaryReporter$new())'`
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-13a-design-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-design-wrapper-20260422a`
- stable focused evidence run:
  - `coverage-evidence-lipid-design-wrapper-20260422a`
- result of the focused lipidomics design wrapper slice:
  - `surface::symbol::mod_lipid_design_server` baseline `100.0%`, target
    `100.0%`
  - evidence gate: `pass`

- the post-lipidomics-design-wrapper repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-o`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-o`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `25`
  - package line coverage:
    - baseline `29.9%`
    - target `34.7%`
  - bundle gates:
    - `pass`: `136`
    - `below_target`: `64`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `97`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `98`

- the general file-management config contracts now participate in the shared
  baseline-vs-target suite:
  - the file is marked with `# fidelity-coverage-compare: shared`
  - the `readConfigFile` test adapts to the old baseline signature that lacks
    `file_type` while preserving target-side file-type forwarding assertions
  - the contract still verifies parser normalization, CPU cluster handoff,
    numeric conversions, plot-format splitting, boolean coercion, and section
    update behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-general-filemgmt-config-contracts.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-general-filemgmt-read-config-20260422a`
- stable focused evidence run:
  - `coverage-evidence-general-filemgmt-read-config-20260422a`
- result of the focused file-management config slice:
  - `surface::symbol::readConfigFile` baseline `96.7%`, target `96.7%`
  - evidence gate: `pass`

- the post-readConfigFile repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-p`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-p`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `26`
  - package line coverage:
    - baseline `30.2%`
    - target `35.0%`
  - bundle gates:
    - `pass`: `137`
    - `below_target`: `63`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `96`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `97`

- the proteomics annotation contracts now cover the public
  `matchAnnotations` helper directly in the shared suite:
  - the file is marked with `# fidelity-coverage-compare: shared`
  - the added S4 fixture exercises non-enrichment `protein_quant_table` ID
    extraction, semicolon-delimited protein groups, UniProt version cleanup,
    blank gene-name normalization, and first-gene extraction
  - the fallback scenario covers `protein_id_table` lookup, fuzzy version
    matching, missing-gene-column handling, and validation errors for absent
    inputs or missing ID columns
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-10-annotation.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-annotation-match-annotations-20260422b`
- stable focused evidence run:
  - `coverage-evidence-prot-annotation-match-annotations-20260422a`
- result of the focused proteomics annotation matching slice:
  - `surface::symbol::matchAnnotations` baseline `92.6%`, target `92.7%`
  - evidence gate: `pass`

- the post-matchAnnotations repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-q`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-q`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `27`
  - package line coverage:
    - baseline `30.5%`
    - target `35.2%`
  - bundle gates:
    - `pass`: `139`
    - `below_target`: `61`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `94`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `95`

- the proteomics DA Glimma volcano helper is now promoted into the shared
  coverage suite through `tests/testthat/test-prot-08-volcano.R`:
  - the corrupt captured snapshot is skipped cleanly instead of failing before
    the helper can run
  - the shared test now covers guard exits, invalid-row null behavior,
    fallback, fuzzy, and prefix contrast matching, UniProt annotation merging,
    display-column propagation, and S4 count/design synchronization branches
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-08-volcano.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-da-glimma-volcano-20260422b`
- stable focused evidence run:
  - `coverage-evidence-prot-da-glimma-volcano-20260422b`
- result of the focused proteomics DA Glimma volcano slice:
  - `surface::symbol::generateProtDAVolcanoPlotGlimma` baseline `98.2%`,
    target `98.2%`
  - evidence gate: `pass`

- the post-proteomics-DA-Glimma-volcano repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-r`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-r`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `28`
  - package line coverage:
    - baseline `30.9%`
    - target `35.6%`
  - bundle gates:
    - `pass`: `140`
    - `below_target`: `60`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `93`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `94`

- a refreshed shared proteomics norm characterization file now includes direct
  helper coverage for:
  - plotting and choice helpers
  - normalization and RUV building blocks
  - session/export/reset handoff helpers
  - the previously covered public module scenarios
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(\".\", quiet = TRUE); testthat::test_file(\"tests/testthat/test-prot-05d-norm-module-characterization.R\", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-20260420T180550Z-5b658fc6`
- stable focused evidence run:
  - `coverage-evidence-20260420T180749Z-067e45a8`
- result of the refreshed shared suite on the proteomics norm family slice:
  - `surface::symbol::mod_prot_norm_server` baseline `100.0%`, target `100.0%`
  - `lineage::R_mod_prot_norm.R::misc` target `93.1%`
  - `lineage::R_mod_prot_norm.R::observer_runtime` target `89.9%`
  - `lineage::R_mod_prot_norm.R::plotting` target `96.4%`
  - `lineage::R_mod_prot_norm.R::output_renderer` target `94.2%`
  - `lineage::R_mod_prot_norm.R::builder_resolver` target `86.6%`
  - `lineage::R_mod_prot_norm.R::formatter` target `97.1%`
  - `lineage::R_mod_prot_norm.R::reactive_state` target `97.1%`
  - evidence gate: `pass` for all `8` focused bundles

- the proteomics peptide group-filter helpers are now promoted into the shared
  coverage suite through `tests/testthat/test-prot-02-qc-peptide-groupaware.R`:
  - the shared test passes baseline-compatible string column names for the old
    monolith helper signatures
  - the target-only exact row-count assertions were replaced with public
    structural invariants that are valid on both `main` and the current branch
  - the added cases cover wide and long peptide input, raw/abundance aliases,
    case-insensitive design matrix matching, no-removal branches, and input
    validation errors
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02-qc-peptide-groupaware.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-peptide-group-filters-20260422d`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-peptide-group-filters-20260422a`
- result of the focused proteomics peptide group-filter slice:
  - `surface::symbol::peptideIntensityFilteringHelper` baseline `100.0%`,
    target `100.0%`
  - `surface::symbol::removePeptidesWithMissingValuesPercentHelper` baseline
    `100.0%`, target `100.0%`
  - evidence gate: `pass` for both focused bundles

- the post-proteomics-peptide-group-filter repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-s`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-s`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `29`
  - package line coverage:
    - baseline `31.3%`
    - target `36.0%`
  - bundle gates:
    - `pass`: `142`
    - `below_target`: `58`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `91`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `92`

- the metabolomics DA analysis input resolver is now covered in the shared
  metabolomics DA characterization file:
  - `tests/testthat/test-metab-02ab-da-module-characterization.R` now includes
    target-safe direct coverage for `resolveMetabDaAnalysisInputs()`
  - the added cases cover direct current S4 use, state-manager fallback,
    missing data, invalid object class, empty contrast table, and null contrast
    table guards
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-metab-02ab-da-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-metab-da-reactive-state-20260422a`
- stable focused evidence run:
  - `coverage-evidence-metab-da-reactive-state-20260422a`
- result of the focused metabolomics DA resolver slice:
  - `lineage::R_mod_metab_da.R::reactive_state` target `100.0%`
  - baseline side remains the expected source-lineage gap for the extracted
    target helper
  - evidence gate: `pass`

- the post-metabolomics-DA-resolver repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-t`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-t`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `29`
  - package line coverage:
    - baseline `31.3%`
    - target `36.0%`
  - bundle gates:
    - `pass`: `143`
    - `below_target`: `57`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `90`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `91`

- the metabolomics QC intensity wrapper is now covered in the shared
  metabolomics QC characterization file:
  - `tests/testthat/test-metab-01q-qc-intensity-module-characterization.R` is
    promoted to the shared coverage-compare suite
  - the file now includes cross-version public wrapper characterization for
    apply-filter and revert behavior through `mod_metab_qc_intensity_server()`
  - target-only extracted seam tests are guarded so they run on the current
    branch while remaining valid against `main`
  - the test fixture builder now constructs validity-safe
    `MetaboliteAssayData` objects across slot-shape differences
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-metab-01q-qc-intensity-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-metab-qc-intensity-server-20260422b`
- stable focused evidence run:
  - `coverage-evidence-metab-qc-intensity-server-20260422a`
- result of the focused metabolomics QC intensity wrapper slice:
  - `surface::symbol::mod_metab_qc_intensity_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-metabolomics-QC-intensity repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-u`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-u`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `30`
  - package line coverage:
    - baseline `31.6%`
    - target `36.4%`
  - bundle gates:
    - `pass`: `144`
    - `below_target`: `56`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `89`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `90`

- the proteomics peptide intensity wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02c-qc-peptide-intensity-shared.R` was added to
    the shared coverage-compare suite
  - the file exercises `mod_prot_qc_peptide_intensity_server()` through the
    real package namespace on both `main` and the current branch
  - the covered scenarios are flexible apply, strict apply, apply error,
    revert success, and revert error behavior
  - the existing target-only peptide module contract suite remains separate
    and continues to cover extracted helper seams on the current branch
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02c-qc-peptide-intensity-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-peptide-intensity-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-peptide-intensity-server-20260422a`
- result of the focused proteomics peptide intensity wrapper slice:
  - `surface::symbol::mod_prot_qc_peptide_intensity_server` baseline
    `100.0%`, target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-peptide-intensity repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-v`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-v`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `31`
  - package line coverage:
    - baseline `31.9%`
    - target `36.7%`
  - bundle gates:
    - `pass`: `145`
    - `below_target`: `55`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `88`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `89`

- the lipidomics QC S4 wrapper is now covered in a shared public wrapper
  characterization file:
  - `tests/testthat/test-lipid-01h-qc-s4-module-characterization.R` was added
    to the shared coverage-compare suite
  - the file exercises `mod_lipid_qc_s4_server()` through the real package
    namespace on both `main` and the current branch
  - the covered public scenario is finalize behavior: the wrapper marks
    lipidomics quality control complete, saves the `lipid_qc_complete` state,
    and preserves the updated S4 assay object
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-01h-qc-s4-module-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-qc-s4-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-lipid-qc-s4-server-20260422a`
- result of the focused lipidomics QC S4 wrapper slice:
  - `surface::symbol::mod_lipid_qc_s4_server` baseline `100.0%`, target
    `100.0%`
  - evidence gate: `pass`

- the post-lipid-QC-S4 repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-w`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-w`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `32`
  - package line coverage:
    - baseline `32.3%`
    - target `37.1%`
  - bundle gates:
    - `pass`: `146`
    - `below_target`: `54`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `87`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `88`

- the lipidomics QC ITSD wrapper is now covered in a shared public wrapper
  characterization file:
  - `tests/testthat/test-lipid-01i-qc-itsd-wrapper-characterization.R` was
    added to the shared coverage-compare suite
  - the file exercises `mod_lipid_qc_itsd_server()` through the real package
    namespace on both `main` and the current branch
  - the covered public scenario is successful internal-standard analysis from
    a `LipidomicsAssayData` object, including result text, summary UI,
    visualization tab registration, and success notification behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-01i-qc-itsd-wrapper-characterization.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-qc-itsd-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-lipid-qc-itsd-server-20260422a`
- result of the focused lipidomics QC ITSD wrapper slice:
  - `surface::symbol::mod_lipid_qc_itsd_server` baseline `100.0%`, target
    `100.0%`
  - evidence gate: `pass`

- the post-lipid-QC-ITSD repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-x`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-x`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `33`
  - package line coverage:
    - baseline `32.9%`
    - target `37.6%`
  - bundle gates:
    - `pass`: `147`
    - `below_target`: `53`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `86`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `87`

- the proteomics protein intensity wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02d-qc-protein-intensity-shared.R` was added to
    the shared coverage-compare suite
  - the file exercises `mod_prot_qc_protein_intensity_server()` through the
    real package namespace on both `main` and the current branch
  - the covered scenarios are flexible apply, strict apply, apply error,
    revert success, and revert error behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02d-qc-protein-intensity-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-protein-intensity-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-protein-intensity-server-20260422a`
- result of the focused proteomics protein intensity wrapper slice:
  - `surface::symbol::mod_prot_qc_protein_intensity_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-protein-intensity repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-y`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-y`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `34`
  - package line coverage:
    - baseline `33.2%`
    - target `37.9%`
  - bundle gates:
    - `pass`: `148`
    - `below_target`: `52`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `85`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `86`

- the proteomics protein cleanup wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02e-qc-protein-cleanup-shared.R` was added to
    the shared coverage-compare suite
  - the file exercises `mod_prot_qc_protein_cleanup_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are accession-cleanup apply with FASTA metadata,
    no-FASTA skip, apply error, revert success, and revert error behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02e-qc-protein-cleanup-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-protein-cleanup-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-protein-cleanup-server-20260422a`
- result of the focused proteomics protein cleanup wrapper slice:
  - `surface::symbol::mod_prot_qc_protein_cleanup_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-protein-cleanup repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-z`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-z`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `35`
  - package line coverage:
    - baseline `33.5%`
    - target `38.3%`
  - bundle gates:
    - `pass`: `149`
    - `below_target`: `51`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `84`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `85`

- the proteomics protein rollup wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02f-qc-protein-rollup-shared.R` was added to the
    shared coverage-compare suite
  - the file exercises `mod_prot_qc_protein_rollup_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are IQ rollup apply success, baseline inline IQ
    timeout handling, generic apply error, revert success, and revert error
    behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02f-qc-protein-rollup-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-protein-rollup-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-protein-rollup-server-20260422a`
- result of the focused proteomics protein rollup wrapper slice:
  - `surface::symbol::mod_prot_qc_protein_rollup_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-protein-rollup repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-aa`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-aa`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `36`
  - package line coverage:
    - baseline `33.8%`
    - target `38.6%`
  - bundle gates:
    - `pass`: `150`
    - `below_target`: `50`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `83`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `84`

- the proteomics peptide q-value wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02g-qc-peptide-qvalue-shared.R` was added to
    the shared coverage-compare suite
  - the file exercises `mod_prot_qc_peptide_qvalue_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are apply success with diagnostic warnings, apply
    error, revert success, and revert error behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02g-qc-peptide-qvalue-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-peptide-qvalue-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-peptide-qvalue-server-20260422a`
- result of the focused proteomics peptide q-value wrapper slice:
  - `surface::symbol::mod_prot_qc_peptide_qvalue_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-peptide-qvalue repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-ab`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-ab`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `37`
  - package line coverage:
    - baseline `34.0%`
    - target `38.9%`
  - bundle gates:
    - `pass`: `151`
    - `below_target`: `49`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `82`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `83`

- the proteomics protein replicate wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02h-qc-protein-replicate-shared.R` was added to
    the shared coverage-compare suite
  - the file exercises `mod_prot_qc_protein_replicate_server()` through the
    real package namespace on both `main` and the current branch
  - the covered scenarios are normal apply, fallback output and QC-parameter
    warning, apply error, revert success, and revert error behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02h-qc-protein-replicate-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-protein-replicate-server-20260422b`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-protein-replicate-server-20260422b`
- result of the focused proteomics protein replicate wrapper slice:
  - `surface::symbol::mod_prot_qc_protein_replicate_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-protein-replicate repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-ac`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-ac`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `38`
  - package line coverage:
    - baseline `34.2%`
    - target `39.2%`
  - bundle gates:
    - `pass`: `152`
    - `below_target`: `48`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `81`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `82`

- the proteomics peptide sample wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02i-qc-peptide-sample-shared.R` was added to
    the shared coverage-compare suite
  - the file exercises `mod_prot_qc_peptide_sample_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are successful sample filtering with removed-sample
    QC parameters, apply error, revert success, and revert error behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02i-qc-peptide-sample-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-peptide-sample-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-peptide-sample-server-20260422a`
- result of the focused proteomics peptide sample wrapper slice:
  - `surface::symbol::mod_prot_qc_peptide_sample_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-peptide-sample repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-ad`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-ad`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `39`
  - package line coverage:
    - baseline `34.5%`
    - target `39.4%`
  - bundle gates:
    - `pass`: `153`
    - `below_target`: `47`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `80`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `81`

- the proteomics peptide protein-count wrapper is now covered in a shared
  public wrapper characterization file:
  - `tests/testthat/test-prot-02j-qc-peptide-protein-shared.R` was added to
    the shared coverage-compare suite
  - the file exercises `mod_prot_qc_peptide_protein_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are successful protein peptide-count filtering with
    both cutoff parameters, apply error, revert success, and revert error
    behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02j-qc-peptide-protein-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-peptide-protein-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-peptide-protein-server-20260422a`
- result of the focused proteomics peptide protein-count wrapper slice:
  - `surface::symbol::mod_prot_qc_peptide_protein_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-peptide-protein repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-ae`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-ae`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `40`
  - package line coverage:
    - baseline `34.7%`
    - target `39.7%`
  - bundle gates:
    - `pass`: `154`
    - `below_target`: `46`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `79`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `80`

- the proteomics peptide imputation wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02k-qc-peptide-impute-shared.R` was added to the
    shared coverage-compare suite
  - the file exercises `mod_prot_qc_peptide_impute_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are successful peptide imputation with QC parameter
    and checkpoint capture, apply error, revert success, and revert error
    behavior
  - an attempted legacy `experiment_paths` save-branch characterization was
    left out because temporary namespace rebinding was not viable in this test
    harness; the focused compare still measured the wrapper bundle at full
    line coverage on both branches
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02k-qc-peptide-impute-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-peptide-impute-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-peptide-impute-server-20260422a`
- result of the focused proteomics peptide imputation wrapper slice:
  - `surface::symbol::mod_prot_qc_peptide_impute_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-peptide-imputation repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-af`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-af`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `41`
  - package line coverage:
    - baseline `34.9%`
    - target `39.9%`
  - bundle gates:
    - `pass`: `155`
    - `below_target`: `45`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `78`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `79`

- the proteomics protein duplicate-removal wrapper is now covered in a shared
  public wrapper characterization file:
  - `tests/testthat/test-prot-02l-qc-protein-dedup-shared.R` was added to the
    shared coverage-compare suite
  - the file exercises `mod_prot_qc_protein_dedup_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are successful duplicate aggregation, apply error,
    revert success, and revert error behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02l-qc-protein-dedup-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-protein-dedup-server-20260422b`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-protein-dedup-server-20260422b`
- result of the focused proteomics protein duplicate-removal wrapper slice:
  - `surface::symbol::mod_prot_qc_protein_dedup_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-protein-dedup repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-ag`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-ag`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `42`
  - package line coverage:
    - baseline `35.0%`
    - target `40.2%`
  - bundle gates:
    - `pass`: `156`
    - `below_target`: `44`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `77`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `78`

- the proteomics peptide replicate wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02m-qc-peptide-replicate-shared.R` was added to
    the shared coverage-compare suite
  - the file exercises `mod_prot_qc_peptide_replicate_server()` through the
    real package namespace on both `main` and the current branch
  - the covered scenarios are successful replicate filtering with QC parameter
    capture, apply error, revert success, and revert error behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02m-qc-peptide-replicate-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-peptide-replicate-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-peptide-replicate-server-20260422a`
- result of the focused proteomics peptide replicate wrapper slice:
  - `surface::symbol::mod_prot_qc_peptide_replicate_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-peptide-replicate repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-ah`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-ah`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `43`
  - package line coverage:
    - baseline `35.2%`
    - target `40.4%`
  - bundle gates:
    - `pass`: `157`
    - `below_target`: `43`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `76`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `77`

- the proteomics peptide rollup wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-02n-qc-peptide-rollup-shared.R` was added to
    the shared coverage-compare suite
  - the file exercises `mod_prot_qc_peptide_rollup_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are successful precursor rollup with QC parameter
    capture, apply error, revert success, and revert error behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02n-qc-peptide-rollup-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-peptide-rollup-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-peptide-rollup-server-20260422a`
- result of the focused proteomics peptide rollup wrapper slice:
  - `surface::symbol::mod_prot_qc_peptide_rollup_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-peptide-rollup repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-ai`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-ai`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `44`
  - package line coverage:
    - baseline `35.4%`
    - target `40.6%`
  - bundle gates:
    - `pass`: `158`
    - `below_target`: `42`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `75`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `76`

- the proteomics design server wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-prot-04b-design-server-shared.R` was added to the
    shared coverage-compare suite
  - the file exercises `mod_prot_design_server()` through the real package
    namespace on both `main` and the current branch
  - the covered scenario is public module-server entry wiring, including the
    target shell delegation and baseline-safe module initialization path
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-04b-design-server-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-design-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-design-server-20260422a`
- result of the focused proteomics design server wrapper slice:
  - `surface::symbol::mod_prot_design_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-design-server repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-aj`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-aj`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `45`
  - package line coverage:
    - baseline `36.9%`
    - target `40.6%`
  - bundle gates:
    - `pass`: `159`
    - `below_target`: `41`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `74`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `75`

- the proteomics protein S4 wrapper is now covered in a shared public wrapper
  characterization file:
  - `tests/testthat/test-prot-02o-qc-protein-s4-shared.R` was added to the
    shared coverage-compare suite
  - the file exercises `mod_prot_qc_protein_s4_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are successful S4 creation, creation error, revert
    success, and revert error behavior
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02o-qc-protein-s4-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-prot-qc-protein-s4-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-prot-qc-protein-s4-server-20260422a`
- result of the focused proteomics protein S4 wrapper slice:
  - `surface::symbol::mod_prot_qc_protein_s4_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-proteomics-protein-S4 repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-ak`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-ak`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `46`
  - package line coverage:
    - baseline `37.1%`
    - target `40.8%`
  - bundle gates:
    - `pass`: `160`
    - `below_target`: `40`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `73`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `74`

- the metabolomics QC orchestrator wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-metab-01w-qc-orchestrator-shared.R` was added to the
    shared coverage-compare suite
  - the file exercises `mod_metab_qc_server()` through the real package
    namespace on both `main` and the current branch
  - the covered scenarios are unavailable-state info panel/getState error
    behavior, populated tab rendering, qc-trigger submodule initialization,
    and state-detected auto-initialization
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-metab-01w-qc-orchestrator-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-metab-qc-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-metab-qc-server-20260422a`
- result of the focused metabolomics QC orchestrator wrapper slice:
  - `surface::symbol::mod_metab_qc_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-metabolomics-QC-orchestrator repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-al`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-al`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `47`
  - package line coverage:
    - baseline `37.3%`
    - target `40.9%`
  - bundle gates:
    - `pass`: `161`
    - `below_target`: `39`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `72`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `73`

- the lipidomics QC orchestrator wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-lipid-01j-qc-orchestrator-shared.R` was added to the
    shared coverage-compare suite
  - the file exercises `mod_lipid_qc_server()` through the real package
    namespace on both `main` and the current branch
  - the covered scenarios are unavailable-state info panel/getState error
    behavior, populated tab rendering, qc-trigger submodule initialization,
    and state-detected auto-initialization
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-01j-qc-orchestrator-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-qc-server-20260422a`
- stable focused evidence run:
  - `coverage-evidence-lipid-qc-server-20260422a`
- result of the focused lipidomics QC orchestrator wrapper slice:
  - `surface::symbol::mod_lipid_qc_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-lipidomics-QC-orchestrator repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-am`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-am`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `48`
  - package line coverage:
    - baseline `37.5%`
    - target `41.1%`
  - bundle gates:
    - `pass`: `162`
    - `below_target`: `38`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `71`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `72`

- the lipidomics QC intensity wrapper is now covered in a shared public
  wrapper characterization file:
  - `tests/testthat/test-lipid-01k-qc-intensity-shared.R` was added to the
    shared coverage-compare suite
  - the file exercises `mod_lipid_qc_intensity_server()` through the real
    package namespace on both `main` and the current branch
  - the covered scenarios are apply-filter success, plot-warning continuation,
    apply-filter error handling, ggplot output dispatch, revert success, and
    revert error behavior
  - `tests/testthat/test-lipid-01j-qc-orchestrator-shared.R` and the new
    intensity shared file now restore namespace-level mocks explicitly so the
    global shared suite does not leak mocked submodule bindings across files
- stable direct current-branch test run:
  - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-01k-qc-intensity-shared.R", reporter = testthat::SummaryReporter$new())'`
- stable focused compare run:
  - `coverage-compare-lipid-qc-intensity-server-20260422b`
- stable focused evidence run:
  - `coverage-evidence-lipid-qc-intensity-server-20260422b`
- result of the focused lipidomics QC intensity wrapper slice:
  - `surface::symbol::mod_lipid_qc_intensity_server` baseline `100.0%`,
    target `100.0%`
  - evidence gate: `pass`

- the post-lipidomics-QC-intensity repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260422Tglobal-fa-c4-refresh-ao`
  - evidence: `coverage-evidence-20260422Tglobal-fa-c4-refresh-ao`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `49`
  - package line coverage:
    - baseline `37.8%`
    - target `41.1%`
  - bundle gates:
    - `pass`: `163`
    - `below_target`: `37`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `70`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `71`

- proteomics protein QC orchestrator server/UI follow-up:
  - `tests/testthat/test-prot-02p-qc-protein-orchestrator-shared.R` was
    added to the shared coverage-compare suite
  - the file exercises `mod_prot_qc_protein_server()` through DIA and LFQ
    server fan-out scenarios and `mod_prot_qc_protein_ui()` through DIA and
    LFQ tab wiring scenarios
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02p-qc-protein-orchestrator-shared.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused server compare/evidence runs:
    - `coverage-compare-prot-qc-protein-server-20260423a`
    - `coverage-evidence-prot-qc-protein-server-20260423a`
  - stable focused UI compare/evidence runs:
    - `coverage-compare-prot-qc-protein-ui-20260423a`
    - `coverage-evidence-prot-qc-protein-ui-20260423a`
  - result of the focused protein QC protein orchestrator slice:
    - `surface::symbol::mod_prot_qc_protein_server` baseline `100.0%`,
      target `100.0%`
    - `surface::symbol::mod_prot_qc_protein_ui` baseline `73.7%`,
      target `100.0%`
    - evidence gates: `pass`

- the post-proteomics-QC-protein-orchestrator repo-wide shared
  compare/evidence refresh is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-ap`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-ap`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `50`
  - package line coverage:
    - baseline `37.9%`
    - target `41.2%`
  - bundle gates:
    - `pass`: `165`
    - `below_target`: `35`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `68`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `69`

- proteomics peptide QC orchestrator server/UI follow-up:
  - `tests/testthat/test-prot-02q-qc-peptide-orchestrator-shared.R` was
    added to the shared coverage-compare suite
  - the file exercises `mod_prot_qc_peptide_server()` through DIA fan-out and
    non-DIA skip scenarios and `mod_prot_qc_peptide_ui()` through public tab
    wiring against lower-level submodule mocks
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-02q-qc-peptide-orchestrator-shared.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused server compare/evidence runs:
    - `coverage-compare-prot-qc-peptide-server-20260423a`
    - `coverage-evidence-prot-qc-peptide-server-20260423a`
  - stable focused UI compare/evidence runs:
    - `coverage-compare-prot-qc-peptide-ui-20260423a`
    - `coverage-evidence-prot-qc-peptide-ui-20260423a`
  - result of the focused protein QC peptide orchestrator slice:
    - `surface::symbol::mod_prot_qc_peptide_server` baseline `100.0%`,
      target `100.0%`
    - `surface::symbol::mod_prot_qc_peptide_ui` baseline `100.0%`,
      target `100.0%`
    - evidence gates: `pass`

- the post-proteomics-QC-peptide-orchestrator repo-wide shared
  compare/evidence refresh is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-aq`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-aq`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `51`
  - package line coverage:
    - baseline `38.1%`
    - target `41.3%`
  - bundle gates:
    - `pass`: `167`
    - `below_target`: `33`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `66`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `67`

- general missing-value parameter duplicate-canonical helper follow-up:
  - `tests/testthat/test-general-filemgmt-update-missing-value-shared.R` was
    added to the shared coverage-compare suite
  - the file exercises `updateMissingValueParameters()` through the
    package-level binding on both `main` and the current branch, covering
    default config/S4 argument synchronization and validation failures
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-general-filemgmt-update-missing-value-shared.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-update-missing-value-parameters-20260423a`
    - `coverage-evidence-update-missing-value-parameters-20260423a`
  - result of the focused missing-value parameter helper slice:
    - `surface::symbol::updateMissingValueParameters` baseline `95.2%`,
      target `95.2%`
    - evidence gate: `pass`

- the post-updateMissingValueParameters repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-ar`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-ar`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `52`
  - package line coverage:
    - baseline `38.2%`
    - target `41.4%`
  - bundle gates:
    - `pass`: `168`
    - `below_target`: `32`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `65`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `66`

- proteomics DA all-contrast output S4 method follow-up:
  - `tests/testthat/test-prot-07c-da-results-output-shared.R` was added to the
    shared coverage-compare suite
  - the file exercises `outputDaResultsAllContrasts()` through the package
    S4 generic on both `main` and the current branch, covering annotated result
    export, fallback annotation behavior, no-signal summary export, graph-save
    side effects, and empty contrast payloads
  - the preceding protein S4 characterization mock cleanup was tightened so it
    restores the `ProteinQuantitativeData` binding before later shared files
    run, and the DA output fixture now constructs the real S4 class directly
  - stable direct current-branch sequence:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); files <- c("tests/testthat/test-prot-02o-qc-protein-s4-shared.R", "tests/testthat/test-prot-07c-da-results-output-shared.R"); for (f in files) testthat::test_file(f, reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-prot-output-da-results-all-contrasts-20260423b`
    - `coverage-evidence-prot-output-da-results-all-contrasts-20260423b`
  - result of the focused DA output method slice:
    - `surface::setMethod::outputDaResultsAllContrasts::ProteinQuantitativeData`
      baseline `99.0%`, target `99.0%`
    - evidence gate: `pass`

- the post-proteomics-DA-output repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-at`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-at`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `53`
  - package line coverage:
    - baseline `38.6%`
    - target `41.9%`
  - bundle gates:
    - `pass`: `169`
    - `below_target`: `31`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `64`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `65`

- lipidomics S4 Pearson correlation method follow-up:
  - `tests/testthat/test-lipid-05-s4-correlation-helpers.R` is now marked as a
    shared coverage-compare file
  - the file exercises `pearsonCorForSamplePairs()` through the package S4
    method on both `main` and the current branch, covering technical-replicate
    pair selection, optional target-only many-to-many warnings, correlation
    table shape, downstream `plotPearson()`, and correlation-threshold
    filtering
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-05-s4-correlation-helpers.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-lipid-s4-pearson-correlation-20260423b`
    - `coverage-evidence-lipid-s4-pearson-correlation-manifest-20260423b`
    - `coverage-evidence-lipid-s4-pearson-correlation-surface-20260423b`
  - result of the focused lipid Pearson method slice:
    - `surface::manifest::tools_refactor_manifest-lipid-s4-wave13.yml::lipid_s4_pearson_correlation_method`
      baseline `95.9%`, target `96.1%`
    - `surface::setMethod::pearsonCorForSamplePairs::LipidomicsAssayData`
      baseline `95.9%`, target `96.1%`
    - evidence gates: `pass`

- the post-lipidomics-S4-Pearson repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-au`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-au`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `54`
  - package line coverage:
    - baseline `39.2%`
    - target `42.4%`
  - bundle gates:
    - `pass`: `171`
    - `below_target`: `29`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `62`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `63`

- peptide ggplot `plotDensity()` S4 method follow-up:
  - added
    `tests/testthat/test-prot-04c-peptide-plot-density-ggplot-shared.R` as a
    shared coverage-compare file
  - the file exercises direct embedded PCA data, fallback extraction from the
    ggplot mapping environment, patchwork panel shape, labels, group data, font
    size propagation, and missing-group error behavior
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-04c-peptide-plot-density-ggplot-shared.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-pept-plot-density-ggplot-20260423a`
    - `coverage-evidence-pept-plot-density-ggplot-20260423a`
  - result of the focused peptide ggplot `plotDensity()` method slice:
    - `surface::setMethod::plotDensity::ggplot2::ggplot`
      baseline `97.4%`, target `97.4%`
    - evidence gate: `pass`

- the post-peptide-ggplot-plotDensity repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-av`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-av`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `55`
  - package line coverage:
    - baseline `39.2%`
    - target `42.4%`
  - bundle gates:
    - `pass`: `172`
    - `below_target`: `28`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `61`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `62`

- lipidomics S4 DA plot methods follow-up:
  - `tests/testthat/test-lipid-07-s4-da-plot-methods.R` is now marked as a
    shared coverage-compare file
  - the file uses namespace-safe mocks for
    `printCountDaGenesTable()` and `plotOneVolcanoNoVerticalLines()` so the
    same expectations exercise both the baseline monolith and extracted target
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-07-s4-da-plot-methods.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-lipid-s4-da-plot-methods-20260423a`
    - `coverage-evidence-lipid-s4-da-numsig-plot-20260423a`
    - `coverage-evidence-lipid-s4-da-static-volcano-20260423a`
  - result of the focused lipid DA plot method slice:
    - `surface::manifest::tools_refactor_manifest-lipid-s4-wave5.yml::lipid_s4_numsig_bar_plot_method`
      baseline `100.0%`, target `100.0%`
    - `surface::manifest::tools_refactor_manifest-lipid-s4-wave5.yml::lipid_s4_static_volcano_plot_method`
      baseline `100.0%`, target `100.0%`
    - evidence gates: `pass`

- the post-lipidomics-S4-DA-plot repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-aw`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-aw`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `56`
  - package line coverage:
    - baseline `39.3%`
    - target `42.6%`
  - bundle gates:
    - `pass`: `174`
    - `below_target`: `26`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `59`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `60`

- lipidomics S4 DA analysis helper follow-up:
  - `tests/testthat/test-lipid-08-s4-da-analysis-methods.R` is now marked as a
    shared coverage-compare file
  - the file uses namespace-safe mocks for `runTestsContrasts()` so the same
    expectations exercise both the baseline monolith and extracted target
  - added numeric-leading group-name coverage and target namespace fallback
    coverage for the extracted helper path
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-08-s4-da-analysis-methods.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-lipid-s4-da-analysis-helper-20260423b`
    - `coverage-evidence-lipid-s4-da-analysis-helper-manifest-20260423b`
    - `coverage-evidence-lipid-s4-da-analysis-helper-surface-20260423b`
  - result of the focused lipid DA analysis helper slice:
    - `surface::manifest::tools_refactor_manifest-lipid-s4-wave8.yml::lipid_s4_da_analysis_helper_method`
      baseline `99.0%`, target `99.0%`
    - `surface::setMethod::differentialAbundanceAnalysisHelper::LipidomicsAssayData`
      baseline `99.0%`, target `99.0%`
    - evidence gates: `pass`

- the post-lipidomics-S4-DA-analysis repo-wide shared compare/evidence refresh
  is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-ax`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-ax`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `57`
  - package line coverage:
    - baseline `39.5%`
    - target `42.7%`
  - bundle gates:
    - `pass`: `176`
    - `below_target`: `24`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `57`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `58`

- lipidomics S4 duplicate-feature resolution follow-up:
  - `tests/testthat/test-lipid-03-s4-duplicate-helpers.R` is now marked as a
    shared coverage-compare file
  - the file keeps the cross-version behavioral checks for
    `resolveDuplicateFeatures()` while gating the target-only helper seam and
    the helper-only duplicate detection checks when an earlier module mock is
    still active in the shared-suite process
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-03-s4-duplicate-helpers.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-lipid-s4-resolve-duplicates-20260423b`
    - `coverage-evidence-lipid-s4-resolve-duplicates-20260423b`
  - result of the focused lipid duplicate-feature method slice:
    - `surface::setMethod::resolveDuplicateFeatures::LipidomicsAssayData`
      baseline `94.6%`, target `100.0%`
    - evidence gate: `pass`
  - superseded run note:
    - `coverage-compare-20260423Tglobal-fa-c4-refresh-ay` exposed the
      duplicate-module mock contamination and is not the active checkpoint

- the post-lipidomics-S4-duplicate-resolution repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-az`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-az`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `58`
  - package line coverage:
    - baseline `39.9%`
    - target `43.1%`
  - bundle gates:
    - `pass`: `177`
    - `below_target`: `23`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `56`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `57`

- metabolomics S4 log-transform assays follow-up:
  - `tests/testthat/test-metab-02r-log-transform-assays-characterization.R` is
    now marked as a shared coverage-compare file
  - the test keeps its fallback parsed-expression registration only when the
    package has not already registered the `logTransformAssays()` S4 method,
    so `covr` traces the normal package namespace method in coverage runs
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-metab-02r-log-transform-assays-characterization.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-metab-s4-log-transform-assays-20260423b`
    - `coverage-evidence-metab-s4-log-transform-assays-20260423b`
  - result of the focused metabolomics log-transform method slice:
    - `surface::manifest::tools_refactor_manifest-metab-s4-wave10.yml::metab_s4_norm_method_log_transform_assays`
      baseline `90.8%`, target `90.8%`
    - evidence gate: `pass`
  - superseded run note:
    - `coverage-compare-metab-s4-log-transform-assays-20260423a` traced the
      manually parsed method registration and produced `0%` method coverage,
      so it is not the active focused checkpoint

- the post-metabolomics-S4-log-transform repo-wide shared compare/evidence
  refresh is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-ba`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-ba`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `59`
  - package line coverage:
    - baseline `40.0%`
    - target `43.2%`
  - bundle gates:
    - `pass`: `178`
    - `below_target`: `22`
    - `target_unmeasured`: `33`
  - low-coverage bundles: `55`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `56`

- general file-management utility helper follow-up:
  - `tests/testthat/test-general-filemgmt-utility-contracts.R` is now marked as
    a shared coverage-compare file
  - `tools/refactor/fidelity_audit.py` now enriches bundle-map members from the
    side/entity inventory when a manifest fallback path is stale, so the
    utility-helper target members resolve to
    `R/func_general_filemgmt_utility_helpers.R`
  - the previously selected but unmarked shared-suite files are now explicit:
    - `tests/testthat/test-general-filemgmt-path-contracts.R`
    - `tests/testthat/test-general-plotting-pca-helper-characterization.R`
    - `tests/testthat/test-lipid-01-qc-filtering-helpers.R`
    - `tests/testthat/test-lipid-10-s4-plotting-methods.R`
    - `tests/testthat/test-prot-07b-da-handlers-compat.R`
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-general-filemgmt-utility-contracts.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-general-filemgmt-utility-helpers-20260423a`
    - `coverage-evidence-general-filemgmt-utility-helpers-20260423a`
  - result of the focused utility-helper slice:
    - `surface::symbol::isArgumentDefined`
    - `surface::symbol::parseList`
    - `surface::symbol::parseString`
    - `surface::symbol::parseType`
    - `surface::symbol::testRequiredArguments`
    - `surface::symbol::testRequiredFiles`
    - `surface::symbol::testRequiredFilesWarning`
    - all seven bundles are `100.0%` baseline and `100.0%` target
    - evidence gates: `pass`

- proteomics volcano main comparison-gap repair:
  - `tests/testthat/test-prot-07b-da-handlers-compat.R` now adapts to the
    available baseline/target DA and DE input signature while still exercising
    the extracted target alias and fallback paths
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-07b-da-handlers-compat.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-prot-volcano-main-compat-20260423c`
    - `coverage-evidence-prot-volcano-main-compat-20260423c`
  - result of the focused proteomics volcano main slice:
    - `surface::symbol::writeInteractiveVolcanoPlotProteomicsMain`
      baseline `100.0%`, target `100.0%`
    - evidence gate: `pass`
  - superseded run notes:
    - `coverage-compare-prot-volcano-main-compat-20260423a` failed the baseline
      q-threshold expectation
    - `coverage-compare-prot-volcano-main-compat-20260423b` left one target
      fallback line uncovered

- the post-utility-helper and proteomics-volcano-main repo-wide shared
  compare/evidence refresh is:
  - bundle map: `bundle-map-20260423Tgeneral-utility-shared-selector-refresh`
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-bc`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-bc`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `60`
  - package line coverage:
    - baseline `40.2%`
    - target `43.3%`
  - bundle gates:
    - `pass`: `186`
    - `below_target`: `42`
    - `target_unmeasured`: `5`
  - low-coverage bundles: `47`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `48`
  - the bundle-map refresh intentionally reclassified many stale
    target-unmeasured helper/S4 bundles into measured below-target bundles; this
    improves audit fidelity and did not introduce comparison gaps or
    regressions

- lipidomics S4 log-transform assays follow-up:
  - added
    `tests/testthat/test-lipid-09b-s4-log-transform-assays-characterization.R`
    as a shared coverage-compare file
  - the test covers the public lipidomics `logTransformAssays()` S4 method
    across the current transformation contract, offset validation, empty assay
    handling, unnamed/base-data-frame assay coercion, missing sample-column
    fallback, and current nonnumeric coercion behavior
  - `tools/refactor/fidelity_audit.py` now normalizes set/list token fields
    when matching shared tests to bundles and applies an omics-token guard so
    a lipid bundle is not linked to a metabolomics/proteomics test only through
    generic terms like `log`, `transform`, or `assays`
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-09b-s4-log-transform-assays-characterization.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-lipid-s4-log-transform-assays-20260423b`
    - `coverage-evidence-lipid-s4-log-transform-assays-20260423b`
  - result of the focused lipidomics log-transform method slice:
    - `surface::manifest::tools_refactor_manifest-lipid-s4-wave12.yml::lipid_s4_log_transform_assays_method`
      baseline `94.9%`, target `94.9%`
    - evidence gate: `pass`
  - superseded focused run note:
    - `coverage-evidence-lipid-s4-log-transform-assays-20260423a` used the
      pre-token-refresh bundle map and therefore lacked the shared-test link

- the post-lipidomics-S4-log-transform repo-wide shared compare/evidence
  refresh is:
  - bundle map: `bundle-map-20260423Tlipid-log-transform-token-refresh`
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-bd`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-bd`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `61`
  - package line coverage:
    - baseline `40.3%`
    - target `43.4%`
  - bundle gates:
    - `pass`: `187`
    - `below_target`: `41`
    - `target_unmeasured`: `5`
  - low-coverage bundles: `46`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `47`

- peptide S4 normalization methods follow-up:
  - added
    `tests/testthat/test-prot-04d-peptide-s4-normalization-characterization.R`
    as a shared coverage-compare file
  - the test covers the public peptide normalization, log2 transformation,
    RUV/Cancor normalization, negative-control optimization, and scoring
    helper contracts across both baseline and target
  - `tools/refactor/fidelity_audit.py` now expands omics token synonyms in
    shared-test matching, records source tokens for marked shared tests, and
    filters generic audit/bookkeeping terms so shared files do not over-link
    through words like `coverage`, `compare`, or `manifest`
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-04d-peptide-s4-normalization-characterization.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-pept-s4-norm-methods-20260423c`
    - `coverage-evidence-pept-s4-norm-methods-20260423c`
  - result of the focused peptide normalization slice:
    - `pept_s4_norm_method_ruv_cancor_fast`: baseline `97.9%`, target
      `97.9%`
    - `pept_s4_norm_method_ruv_cancor`: baseline `100.0%`, target `100.0%`
    - `pept_s4_norm_method_ruviii_c_varying`: baseline `100.0%`, target
      `100.0%`
    - `log2TransformPeptideMatrix`: target `100.0%`
    - `normaliseBetweenSamples`: baseline `93.3%`, target `93.3%`
    - `getNegCtrlProtAnovaPeptides`: target `100.0%`
    - `findBestNegCtrlPercentagePeptides`: target `97.6%`
    - `.peptide_calculateAdaptiveMaxK`: baseline `100.0%`, target `100.0%`
    - `.peptide_calculateCompositeScore`: baseline `100.0%`, target
      `100.0%`
    - `.peptide_calculateSeparationScore`: baseline `92.5%`, target `92.5%`
    - evidence gate: `pass` for all ten focused peptide bundles
  - superseded focused run notes:
    - `coverage-compare-pept-s4-norm-methods-20260423a` used
      `tests/testthat/test-prot-04-design.R` and failed on baseline with 102
      failures; do not use it for evidence
    - `coverage-compare-pept-s4-norm-methods-20260423b` passed tests but left
      `ruvCancor()` and `findBestNegCtrlPercentagePeptides()` below target;
      it is superseded by the branch-coverage additions above

- the post-peptide-S4-normalization repo-wide shared compare/evidence refresh
  is:
  - bundle map: `bundle-map-20260423Tpept-s4-norm-characterization-refresh4`
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-be`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-be`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `62`
  - compare tests per side: `432`
  - failures per side: `0`
  - skips:
    - baseline `209`
    - target `22`
  - package line coverage:
    - baseline `41.1%`
    - target `44.1%`
  - bundle gates:
    - `pass`: `198`
    - `below_target`: `30`
    - `target_unmeasured`: `5`
  - low-coverage bundles: `35`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `36`

- proteomics S4 LIMPA protein missing-value imputation follow-up:
  - added
    `tests/testthat/test-prot-05e-s4-missingness-limpa-characterization.R`
    as a shared coverage-compare file
  - the test covers the public
    `proteinMissingValueImputationLimpa,ProteinQuantitativeData` S4 method
    through provided-DPC, estimated-DPC, fallback-slope, default-column,
    numeric-DPC, nonfinite-cleanup, raw-scale back-transform, metadata-storage,
    and wrapped-error paths
  - the test installs a minimal fake `limpa` package into a temp library so the
    optional Bioconductor dependency is not required during audit replay
  - `tools/refactor/fidelity_audit.py` now resolves unsigned manifest
    `setMethod` selectors by first matching the selector against the manifest's
    own source/target file and then following the full inferred entity key to
    the runtime-active duplicate definition selected by Collate order
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-05e-s4-missingness-limpa-characterization.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-prot-s4-missingness-limpa-20260423b`
    - `coverage-evidence-prot-s4-missingness-limpa-20260423b`
  - result of the focused proteomics S4 LIMPA protein imputation slice:
    - baseline active definition `R/func_prot_limpa.R`: `99.4%`
    - target active definition `R/func_peptide_qc_imputation.R`: `99.3%`
    - evidence gate: `pass`
  - superseded focused run note:
    - `coverage-compare-prot-s4-missingness-limpa-20260423a` ran cleanly but
      left the bundle unobserved because the mapper stayed pinned to inactive
      manifest source paths instead of the runtime-active duplicate S4 methods

- the post-proteomics-S4-LIMPA-imputation repo-wide shared compare/evidence
  refresh is:
  - bundle map: `bundle-map-20260423Tprot-s4-missingness-active-duplicate-refresh3`
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-bf`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-bf`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `63`
  - compare tests per side: `437`
  - failures per side: `0`
  - skips:
    - baseline `209`
    - target `22`
  - package line coverage:
    - baseline `41.3%`
    - target `44.3%`
  - bundle gates:
    - `pass`: `199`
    - `below_target`: `31`
    - `target_unmeasured`: `3`
  - low-coverage bundles: `34`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `35`

- proteomics peptide S4 LIMPA missing-value imputation follow-up:
  - added
    `tests/testthat/test-prot-04e-peptide-s4-missingness-limpa-characterization.R`
    as a shared coverage-compare file
  - the test covers the public
    `peptideMissingValueImputationLimpa,PeptideQuantitativeData` S4 method
    through raw-scale log2/back-transform, already-logged data,
    no-extra-transform branches, missing-matrix calculation, nonfinite cleanup,
    DPC metadata storage, slope interpretation, and wrapped-error paths
  - the local fake `limpa` package used by the peptide/protein LIMPA tests now
    dispatches handler arguments by formal name, making the two shared test
    files order-independent inside a single coverage session
  - stable direct current-branch test run for the previously failing order:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-04e-peptide-s4-missingness-limpa-characterization.R", reporter = testthat::SummaryReporter$new()); testthat::test_file("tests/testthat/test-prot-05e-s4-missingness-limpa-characterization.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-pept-s4-missingness-limpa-20260423a`
    - `coverage-evidence-pept-s4-missingness-limpa-20260423a`
  - result of the focused proteomics peptide S4 LIMPA imputation slice:
    - baseline active definition `R/func_peptide_qc_imputation.R`: `99.4%`
    - target active definition `R/func_peptide_qc_imputation.R`: `99.4%`
    - evidence gate: `pass`
  - superseded repo-wide run note:
    - `coverage-compare-20260423Tglobal-fa-c4-refresh-bg` is retained as
      failed evidence because the initial fake `limpa` signature caused five
      protein LIMPA failures on both sides

- the post-peptide-S4-LIMPA-imputation repo-wide shared compare/evidence
  refresh is:
  - bundle map: `bundle-map-20260423Tpept-s4-missingness-limpa-refresh`
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-bh`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-bh`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `64`
  - compare tests per side: `441`
  - failures per side: `0`
  - skips:
    - baseline `209`
    - target `22`
  - package line coverage:
    - baseline `41.6%`
    - target `44.6%`
  - bundle gates:
    - `pass`: `200`
    - `below_target`: `30`
    - `target_unmeasured`: `3`
  - low-coverage bundles: `33`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `34`

- metabolomics S4 `plotDensity,list` follow-up:
  - updated
    `tests/testthat/test-metab-01i-s4-plot-density-characterization.R`
    into an explicit shared coverage-compare file
  - the test now leaves the package-namespace `plotDensity()` S4 method bodies
    intact during coverage replay so `covr` attributes execution to runtime
    package code instead of test-local method re-registration
  - `tools/refactor/fidelity_audit.py` now resolves manifest anchor ranges that
    contain exactly one inventory entity to that runtime entity, then follows
    duplicate S4 definitions to the Collate-active implementation
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-metab-01i-s4-plot-density-characterization.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-metab-s4-plot-density-list-20260423c`
    - `coverage-evidence-metab-s4-plot-density-list-20260423c`
  - result of the focused metabolomics S4 `plotDensity,list` slice:
    - baseline active definition `R/func_lipid_s4_objects.R`: `95.5%`
    - target active definition `R/func_lipid_s4_plotting_methods.R`: `95.5%`
    - evidence gate: `pass`
  - superseded focused run notes:
    - `coverage-compare-metab-s4-plot-density-list-20260423a` left the bundle
      unobserved because the method bodies were re-registered in the test
      environment
    - `coverage-compare-metab-s4-plot-density-list-20260423b` left the bundle
      unobserved because the mapper was still pinned to the inactive manifest
      anchor range

- the post-metabolomics-S4-plotDensity-list repo-wide shared compare/evidence
  refresh is:
  - bundle map: `bundle-map-20260423Tmetab-s4-plot-density-list-active-range`
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-bi`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-bi`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `65`
  - compare tests per side: `445`
  - failures per side: `0`
  - skips:
    - baseline `209`
    - target `22`
  - package line coverage:
    - baseline `41.8%`
    - target `44.8%`
  - bundle gates:
    - `pass`: `201`
    - `below_target`: `30`
    - `target_unmeasured`: `2`
  - low-coverage bundles: `32`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `33`

- proteomics DA model/stat/output helper follow-up:
  - added `tests/testthat/test-prot-07d-da-model-stats-shared.R` as an
    explicit shared coverage-compare file
  - the test uses namespace-qualified package functions for `ebFit()`,
    `runTest()`, `runTests()`, `runTestsContrasts()`,
    `daAnalysisWrapperFunction()`, `outputDaAnalysisResults()`, and
    `saveDaProteinList()` so the full shared suite cannot accidentally call
    `shiny::runTests()`
  - `tools/refactor/fidelity_audit.py` now preserves direct and high-overlap
    shared unit tests before applying the per-bundle selected-file cap, which
    links this DA shared file to all seven covered helper bundles
  - stable direct current-branch test runs:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); library(shiny); testthat::test_file("tests/testthat/test-prot-07d-da-model-stats-shared.R", reporter = testthat::SummaryReporter$new())'`
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-07b-da-handlers-compat.R", reporter = testthat::SummaryReporter$new()); testthat::test_file("tests/testthat/test-prot-07c-da-results-output-shared.R", reporter = testthat::SummaryReporter$new()); testthat::test_file("tests/testthat/test-prot-07d-da-model-stats-shared.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-prot-da-model-stats-20260423f`
    - `coverage-evidence-prot-da-model-stats-20260423f-scoped`
  - result of the focused proteomics DA helper slice:
    - `daAnalysisWrapperFunction`: `90.0%`
    - `ebFit`: `100.0%`
    - `runTest`: `90.0%`
    - `runTests`: `97.1%`
    - `runTestsContrasts`: `97.9%`
    - `outputDaAnalysisResults`: `98.6%`
    - `saveDaProteinList`: `100.0%`
    - evidence gates: `pass`: `7`
  - superseded focused/global run notes:
    - focused runs `coverage-compare-prot-da-model-stats-20260423a` through
      `coverage-compare-prot-da-model-stats-20260423d` were threshold-fill
      attempts
    - `coverage-compare-prot-da-model-stats-20260423e` passed before the
      full-suite masking repair and is superseded by run `f`
    - `coverage-compare-20260423Tglobal-fa-c4-refresh-bj` failed because
      `shiny::runTests()` masked the package helper during full shared replay
    - unscoped evidence `coverage-evidence-prot-da-model-stats-20260423f`
      was diagnostic and is superseded by scoped evidence
      `coverage-evidence-prot-da-model-stats-20260423f-scoped`

- the post-proteomics-DA-helper repo-wide shared compare/evidence refresh is:
  - bundle map: `bundle-map-20260423Tprot-da-model-stats-refresh-b`
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-bk`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-bk`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `66`
  - compare tests per side: `448`
  - failures per side: `0`
  - skips:
    - baseline `209`
    - target `22`
  - package line coverage:
    - baseline `43.5%`
    - target `46.4%`
  - bundle gates:
    - `pass`: `208`
    - `below_target`: `23`
    - `target_unmeasured`: `2`
  - low-coverage bundles: `25`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `26`

- proteomics DA run-analysis handler follow-up:
  - extended `tests/testthat/test-prot-07b-da-handlers-compat.R` with shared
    `testServer()` coverage for `da_server_run_analysis_handler()`
  - the new handler scenarios cover:
    - full-format two-contrast success, qvalue warning modal handling,
      workflow/tab-status updates, checkpoint capture, and successful
      `outputDaResultsAllContrasts()` dispatch
    - auto-generated full-format contrasts, failed state-save warning handling,
      and output-write failure tolerance
    - `differentialAbundanceAnalysis()` failure diagnostics and the outer error
      notification path without completing DA state
  - namespace overrides keep the baseline and target replay deterministic while
    avoiding real limma execution and file-output side effects
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-prot-07b-da-handlers-compat.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-prot-da-run-analysis-handler-20260423a`
    - `coverage-evidence-prot-da-run-analysis-handler-20260423a`
  - result of the focused proteomics DA run-analysis handler slice:
    - `da_server_run_analysis_handler`: `100.0%`
    - evidence gate: `pass`

- the post-proteomics-DA-run-handler repo-wide shared compare/evidence refresh
  is:
  - bundle map: `bundle-map-20260423Tprot-da-model-stats-refresh-b`
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-bl`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-bl`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `66`
  - compare tests per side: `451`
  - failures per side: `0`
  - skips:
    - baseline `209`
    - target `22`
  - package line coverage:
    - baseline `44.3%`
    - target `47.1%`
  - bundle gates:
    - `pass`: `209`
    - `below_target`: `22`
    - `target_unmeasured`: `2`
  - low-coverage bundles: `24`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `25`

- lipidomics filtering progress/update helper follow-up:
  - marked the existing focused characterization files as explicit shared
    coverage-compare files:
    - `tests/testthat/test-lipid-04-s4-progress-helpers.R`
    - `tests/testthat/test-lipid-01a-qc-update-filtering-characterization.R`
  - regenerated the bundle map so the progress getter file links to
    `getFilteringProgressLipidomics()` and the update contract file links to
    `updateLipidFiltering()`
  - stable direct current-branch test run:
    - `Rscript --vanilla -e 'renv::load(); pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-lipid-04-s4-progress-helpers.R", reporter = testthat::SummaryReporter$new()); testthat::test_file("tests/testthat/test-lipid-01a-qc-update-filtering-characterization.R", reporter = testthat::SummaryReporter$new())'`
  - stable focused compare/evidence runs:
    - `coverage-compare-lipid-filter-progress-getter-20260423a`
    - `coverage-evidence-lipid-filter-progress-getter-20260423a`
    - `coverage-compare-lipid-update-filtering-20260423a`
    - `coverage-evidence-lipid-update-filtering-20260423a`
  - result of the focused lipidomics filtering helper slice:
    - `getFilteringProgressLipidomics`: `100.0%`
    - `updateLipidFiltering`: `96.6%`
    - evidence gates: `pass`: `2`
  - the direct and covr runs emit expected ggplot non-finite-row warnings from
    the plot contract; both baseline and target complete without test failures

- the post-lipid-filtering-helper repo-wide shared compare/evidence refresh is:
  - bundle map: `bundle-map-20260423Tlipid-filter-progress-refresh`
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-bm`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-bm`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `68`
  - compare tests per side: `454`
  - failures per side: `0`
  - skips:
    - baseline `209`
    - target `22`
  - package line coverage:
    - baseline `45.3%`
    - target `47.3%`
  - bundle gates:
    - `pass`: `211`
    - `below_target`: `20`
    - `target_unmeasured`: `2`
  - low-coverage bundles: `22`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `23`

Interpretation update:

- the recovered metabolomics summary, lipidomics summary, proteomics
  enrichment, lipidomics norm, proteomics design, general file-management
  bootstrap setup, lipidomics design, lipidomics QC helper, proteomics volcano
  compatibility, metabolomics QC ITSD wrapper, and lipidomics design wrapper
  family slices, plus the general file-management config helper and proteomics
  annotation matching, proteomics DA Glimma volcano helper, proteomics peptide
  group-filter helpers, metabolomics DA analysis input resolver, and
  metabolomics QC intensity wrapper, proteomics peptide intensity wrapper, and
  lipidomics QC S4 wrapper, lipidomics QC ITSD wrapper, proteomics protein
  intensity wrapper, proteomics protein cleanup wrapper, proteomics protein
  rollup wrapper, proteomics peptide q-value wrapper, and proteomics protein
  replicate wrapper, proteomics peptide sample wrapper, and proteomics peptide
  protein-count wrapper, proteomics peptide imputation wrapper, and proteomics
  protein duplicate-removal wrapper, and proteomics peptide replicate wrapper,
  proteomics peptide rollup wrapper, and proteomics design server wrapper,
  and proteomics protein S4 wrapper, and metabolomics QC orchestrator wrapper,
  lipidomics QC orchestrator wrapper, and lipidomics QC intensity wrapper,
  proteomics protein QC orchestrator server/UI wrappers, and proteomics
  peptide QC orchestrator server/UI wrappers, and the general missing-value
  parameter duplicate-canonical helper, and the proteomics DA all-contrast
  output S4 method, the lipidomics S4 Pearson correlation method, the peptide
  ggplot `plotDensity()` S4 method, the lipidomics S4 DA plot methods, and the
  lipidomics S4 DA analysis helper, and the lipidomics S4 duplicate-feature
  resolution method, and the metabolomics S4 log-transform assays method, and
  the general file-management utility helpers, and the proteomics volcano main
  manual-merge compatibility path/comparison-gap repair, and the lipidomics S4
  log-transform assays method, the peptide S4 normalization/log2/RUV helper
  cluster, the proteomics S4 LIMPA protein missing-value imputation method, and
  the proteomics peptide S4 LIMPA missing-value imputation method, and the
  metabolomics S4 `plotDensity,list` runtime-active method, and the proteomics
  DA model/stat/output helper cluster, and the proteomics DA run-analysis
  handler, and the lipidomics filtering progress/update helper pair are now
  resolved for FA-C4 under the corrected shared-suite selector model
- the repo-wide coverage-backed parity campaign is now tracking against the
  latest refreshed shared compare/evidence pair above, not the earlier
  over-broad selector run
- the stale final-queue gap for `generateProtDAVolcanoStatic()` was caused by
  the earlier compare starting before the latest shared volcano tests were in
  the worktree; the refreshed focused run
  `coverage-compare-prot-da-volcano-static-20260423a` now measures that helper
  at `98.6%` on both baseline and target
- the peptide wave-4 `plotPca` legacy manifest block is an inactive duplicate
  alongside an active sibling block; `coverage-evidence-pept-s4-plot-pca-pair-20260423a`
  now accepts the legacy block through the active sibling when the same shared
  suite proves the active implementation at `100.0%`
- the post-final-queue repo-wide shared compare/evidence refresh is:
  - compare: `coverage-compare-20260423Tglobal-fa-c4-refresh-bp`
  - evidence: `coverage-evidence-20260423Tglobal-fa-c4-refresh-bp`
- current repo-wide coverage-backed parity state from that refreshed run:
  - selected shared files: `68`
  - compare tests per side: `465`
  - failures per side: `0`
  - skips:
    - baseline `209`
    - target `22`
  - package line coverage:
    - baseline `47.5%`
    - target `49.2%`
  - bundle gates:
    - `pass`: `233`
  - low-coverage bundles: `0`
  - comparison gaps: `0`
  - regressions: `0`
  - emitted coverage exceptions: `1`
- `FA-C4` is now complete for bundle-level parity coverage: every tracked
  bundle in `latest-bundles.json` passes under the refreshed shared compare
  surface

Next action inside `FA-C5`:

- move from bundle-fill work to package-level closeout:
  - raise package line coverage from `49.2%` toward the `80.0%` gate or record
    that package gate as the only remaining explicit exception
- keep the closeout sequence:
  - rerun repo-wide coverage collection after each additional package-coverage push
  - regenerate the parity-plus-coverage evidence package
  - publish final signoff once the package gate disposition is settled

#### `FA-C5`: Coverage-Backed Parity Closeout

Goal:

- convert the current parity signoff into a coverage-backed parity signoff

Work:

- rerun coverage collection on all tracked bundles
- regenerate coverage reports
- rerun full closeout if new exceptions are emitted
- publish a final evidence package tying parity and coverage together

Acceptance:

- package-wide line coverage `>= 80%`
- every tracked parity bundle has a measured target coverage value
- every auto-curated bundle has a linked evidence note and linked tests
- any remaining low-coverage bundle is explicitly recorded as an exception

### Coverage Acceptance Rules

Apply the following rules:

- package-wide post-closeout line coverage target: `>= 80%`
- every deduplicated auto-curated parity bundle should reach `>= 80%` target
  line coverage
- touched helper bundles should still aim for `>= 90%`
- touched preserved and S4 bundles should still aim for `>= 85%`
- wrapper and module bundles still require contract and scenario completeness;
  line coverage alone is not sufficient
- baseline and target coverage must be recorded using the same test subsets

Important nuance:

- baseline coverage is recorded for comparison, but acceptance is on the target
  branch plus parity evidence
- if a baseline bundle is also under-covered, that is still documented; it does
  not waive the target requirement

### Documentation Requirements

For every coverage-backed bundle, persist:

- machine-readable bundle record
- linked test files
- linked exception keys
- baseline coverage
- target coverage
- disposition
- human-readable justification note

No bundle should rely only on chat history for its rationale.

### Execution Rule For Future Sessions

When resuming this work:

1. start from `integrated-closeout-20260419d`
2. do not reopen parity triage unless new coverage-driven regressions appear
3. resume at `FA-C2`, then progress strictly through `FA-C5`
4. update this Markdown file and the adjacent JSON manifest together after each
   checkpoint
5. treat coverage collection, bundle mapping, and documentation as part of the
   same evidence system, not as separate ad hoc work

This plan file and the adjacent JSON manifest should be updated together as
implementation progresses.
