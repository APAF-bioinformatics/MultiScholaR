# Refactor Tooling

Operational note:

- Codex acts as the outer supervisor for stabilization work in this repo.
- The inner loop runs one bounded checkpoint at a time.
- After each checkpoint, Codex decides go/no-go for the next loop based on
  verification and only stops for real strategic blockers.
- Production commits are made at module boundaries, not after every small seam.
- In default mode, this supervision happens in-chat.
- In AFK mode, if the user explicitly authorizes unattended progress, one
  bounded local supervisor may perform the same supervision cycle on Codex's
  behalf.
- Keep the runtime bounded to one supervisor, one inner loop, and one verifier
  at a time.

Current post-design state:

- Proteomics design/builder is complete and committed in `77fb6c7`
  (`Stabilize proteomics design and builder wrappers`).
- Proteomics completion is now `4/8` buckets (`50.0%`).
- Repo-wide completion is now `4/13` buckets (`30.8%`).
- There is no active supervisor or loop by default after the design boundary.
- If parallelism is used from this point forward, it must be done only in
  separate git worktrees with isolated `.god-module-stabilization/` state.
- Same-worktree parallel loops are still not safe because they compete on
  `DESCRIPTION`, handover/backlog docs, commits, and runtime state.
- The first recommended separate-worktree lipidomics pilot target is
  `R/func_lipid_qc.R`; it is lower-touch than `R/mod_lipid_norm.R` because it
  is direct-extraction-ready and has no Shiny wrapper shell.

This directory contains manifest-driven tooling for splitting large `R/` files
into smaller files without hand-rewriting function bodies.

The design goal is simple:

1. Humans decide destination files and ordering.
2. Tooling moves exact source ranges.
3. Verification checks that selectors, targets, and parse state still make sense.

## Files

- `extract_blocks.R`: manifest-driven source extractor
- `verify_refactor.R`: manifest and parse verifier
- `check_wave_apply.R`: post-apply checker for a wave already materialized into `R/`
- `audit_refactor_coupling.R`: filename-coupling audit for `R/*.R` references outside the manifest flow
- `audit_file_sizes.R`: file-size budget audit for the whole repo or a single wave manifest
- `stabilization-status.py`: machine-readable target/backlog progress estimator for the stabilization workflow
- `stabilization-reviewer.py`: local checkpoint reviewer that reruns reported test commands and verifies status
- `stabilization-codex-runner.py`: default noninteractive `codex exec` runner for loop iterations
- `stabilization-loop.py`: resumable outer loop controller for one-checkpoint-at-a-time stabilization
- `stabilization-context-poller.py`: passive poller that records loop status plus an explicit supervision reminder for later Codex turns
- `PLAYBOOK.md`: operational protocol for later waves
- `GOD_MODULE_STABILIZATION_BACKLOG.md`: prioritized queue for stabilization-first refactors
- `manifest-wave1.yml`: initial Wave 1 proteomics DA scaffold

## Manifest Shape

Top-level keys:

- `version`: schema version
- `wave`: human label
- `source_root`: usually `R`
- `collate_mode`: informational for now
- `defaults`: default extractor behavior
- `collate_targets`: desired target file order
- `entries`: ordered extraction plan

Entry keys:

- `id`: stable identifier
- `action`: one of `extract`, `manual_merge`, `skip`
- `source`: source file path
- `selector`: how to locate the block
- `target`: destination file path
- `group`: logical grouping label
- `note`: optional human note

## Selector Kinds

### `symbol`

Matches a top-level assignment to a symbol.

```yaml
selector:
  kind: symbol
  value: mod_prot_da_server
```

Supported assignment forms:

- `name <- function(...)`
- `name = function(...)`

### `setMethod`

Matches a top-level `setMethod()` call by method name.

```yaml
selector:
  kind: setMethod
  value: differentialAbundanceAnalysis
```

The extractor matches the second argument to `setMethod(...)` after normalizing
either symbol or string forms.

### `setGeneric`

Matches a top-level `setGeneric()` call by generic name.

```yaml
selector:
  kind: setGeneric
  value: plotPca
```

### `setClass`

Matches a top-level `setClass()` call by class name.

```yaml
selector:
  kind: setClass
  value: MetaboliteAssayData
```

### `expr_index`

Matches a top-level parse expression by 1-based parse order.

Use this sparingly for `NULL` doc blocks or other non-symbol expressions.

```yaml
selector:
  kind: expr_index
  value: 1
```

### `anchor_range`

Matches a raw line range between two regex anchors.

```yaml
selector:
  kind: anchor_range
  start: "^# Handler 4: Volcano Plot Rendering$"
  end: "^# Handler 5: Heatmap Rendering$"
  end_inclusive: false
```

This is the fallback when the source block is not its own top-level expression.

## Actions

### `extract`

Normal automatic extraction. The block is copied into the target file.

### `manual_merge`

Do not auto-copy into the target. The selector still must resolve, but the
extractor reports it as unresolved for human review.

Use this for real symbol collisions or intentionally divergent duplicates.

Important: `manual_merge` does not remove the old block from the source file.
If the merge is resolved in staging, the source cleanup still needs an explicit
follow-up step when the wave is applied to live `R/`.

### `skip`

Validate the selector exists, but do not extract it.

## Leading Comment Behavior

For parse-based selectors, the extractor includes contiguous comment lines
immediately above the matched top-level expression. This preserves nearby
roxygen and section comments without swallowing the file header.

It does not attempt to synthesize new roxygen blocks.

## Usage

Dry run:

```bash
Rscript tools/refactor/extract_blocks.R --manifest tools/refactor/manifest-wave1.yml --dry-run
```

Write target files only:

```bash
Rscript tools/refactor/extract_blocks.R --manifest tools/refactor/manifest-wave1.yml --write-targets
```

Write target files into a staging directory instead of the repo root:

```bash
Rscript tools/refactor/extract_blocks.R \
  --manifest tools/refactor/manifest-wave1.yml \
  --write-targets \
  --output-root tools/refactor/staging/wave1
```

Write target files and rewrite sources:

```bash
Rscript tools/refactor/extract_blocks.R --manifest tools/refactor/manifest-wave1.yml --write-targets --rewrite-sources
```

Emit collate list:

```bash
Rscript tools/refactor/extract_blocks.R --manifest tools/refactor/manifest-wave1.yml --emit-collate tools/refactor/collate-wave1.txt
```

Verify:

```bash
Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-wave1.yml
```

Post-apply check:

```bash
Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-wave1.yml
```

Filename-coupling audit:

```bash
Rscript tools/refactor/audit_refactor_coupling.R
Rscript tools/refactor/audit_refactor_coupling.R --output tools/refactor/AUDIT-filename-coupling.md
```

File-size audit:

```bash
Rscript tools/refactor/audit_file_sizes.R
Rscript tools/refactor/audit_file_sizes.R --output tools/refactor/AUDIT-file-sizes.md
Rscript tools/refactor/audit_file_sizes.R --manifest tools/refactor/manifest-wave1.yml --output tools/refactor/AUDIT-wave1-file-sizes.md
```

Stabilization status:

```bash
python3 tools/refactor/stabilization-status.py --target R/mod_prot_norm.R --backlog tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md --json
```

Run one stabilization loop iteration:

```bash
python3 tools/refactor/stabilization-loop.py run --iteration-limit 1
```

Run the passive reminder poller:

```bash
python3 tools/refactor/stabilization-context-poller.py run
```

Inspect current loop state:

```bash
python3 tools/refactor/stabilization-loop.py status
```

Recover a stale active run before relaunching:

```bash
python3 tools/refactor/stabilization-loop.py recover
```

Current loop caveat:

- use the stabilization loop for explicit, bounded targets and status polling
- explicit `--target` / `--item-id` override is supported; use it when backlog
  state may be stale
- do not rely on backlog text alone for auto-queueing until readiness checks
  and backlog/live drift detection are added
- `python3 tools/refactor/stabilization-loop.py status` now returns a safe JSON
  error if the loop state file does not exist yet
- `status` now separates workflow state from runtime liveness:
  a target can still be `in_progress` while runtime is `idle` between bounded
  loop invocations
- active loop runs now record PID/heartbeat metadata and executor/reviewer
  stdout/stderr artifacts under `.god-module-stabilization/`
- `status` and `recover` both reconcile stale `activeRun` metadata into a
  blocked terminal transition instead of leaving ambiguous live-state markers
- `stabilization-reviewer.py` now replays sanitized direct commands instead of
  `bash -lc`, strips common checkpoint prose suffixes, rejects shell-control
  operators, and deduplicates duplicate replay commands
- the preferred replay contract is now typed:
  `verification.replayCommands[].argv`
  with human commentary kept separately in `summary`, `notes`, or
  `verification.display`
- `stabilization-supervisor.py` now treats both `blocked` and `stopped` as hold
  states and exits instead of thrashing relaunches
- `python3 tools/refactor/stabilization-loop.py run --clear-stop-file ...`
  exists for deliberate resumes after a real stop marker
- `python3 tools/refactor/stabilization-loop.py recover`
  can now rerun a stored blocked reviewer payload and append an explicit
  recovery event when the original block was a false positive
- real `codex exec` loop runs still require network access to the Codex backend;
  DNS/network failures are classified as blocked infrastructure runs
- the loop now supports repo-local target overrides in
  `tools/refactor/stabilization-target-overrides.json`
  so a backlog bucket can keep its public wrapper target while the executor is
  aimed at a different late-stage structural work file
- the replay contract now distinguishes replay-safe verification from
  one-shot operational commands:
  - `replayMode: "verify"` means the reviewer should rerun the command
  - `replayMode: "operate"` means the reviewer should record but not rerun it
  - when omitted, the reviewer infers replay safety from `phase` and command
    content, so `apply_wave.py` and similar mutation steps no longer cause
    false blocked states on successful wave applies
- the loop prompt now emits a compact checkpoint brief with the latest
  checkpoint, latest summary, explicit primary work file, focused gate
  commands, and a trimmed handover excerpt
- current practical caveat: legacy `testsRun` fallback still exists for older
  stored payloads, but new loop iterations should emit structured
  `verification.replayCommands` instead of free-form replay strings

Current live state at compact time:

- active bucket: `5-proteomics-design-and-builder`
- loop status: `blocked`
- runtime status: `idle`
- iteration: `189 / 200`
- latest clean checkpoint: `builder_state_setup_checkpoint`
- current hold reason: `blocked_executor_timeout`
- the latest timeout exposed late-stage prompt drift rather than a code/test
  regression: the bucket target stayed on `R/mod_prot_design.R` while the real
  remaining work had moved into `R/mod_prot_design_builder.R`
- the new override layer now pins:
  - handover: `tools/refactor/HANDOVER-prot-design-seams.md`
  - primary work file: `R/mod_prot_design_builder.R`
  - focused gate: `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
- next safe step: reopen the blocked design bucket under the hardened compact
  prompt and prove one clean bounded builder-file checkpoint before re-enabling
  unattended supervisor mode

Practical replay-policy update:

- the design bucket later hit a false blocked reviewer replay on
  `builder_pre_server_wave2_applied_checkpoint` because the old reviewer tried
  to rerun `apply_wave.py` after the manifest had already been applied
- that is now fixed:
  - the reviewer skips non-replay-safe operational commands
  - `stabilization-loop.py recover` can reopen the old blocked state cleanly
  - the design bucket recovered back to `in_progress` without manual state
    edits

## Safety Notes

- Default mode is dry-run.
- Source rewriting is opt-in.
- `manual_merge` entries are never auto-applied.
- `manual_merge` entries also are not auto-removed from their original source files.
- The scripts preserve exact source text for matched blocks.
- They do not regenerate code from semantics.
