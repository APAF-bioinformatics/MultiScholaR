# God Module Stabilization Loop Spec

## Operating Mode

Default operating mode for this repo is now:

1. Codex is the outer supervisor.
2. The inner stabilization loop runs one bounded checkpoint at a time.
3. After each inner-loop checkpoint, Codex decides whether to relaunch the next
   checkpoint without asking the user again.
4. Codex only stops when there is a real blocker that needs strategic
   elevation, redesign, or explicit user direction.
5. Production commits happen at module boundaries, not after every checkpoint.
6. In default mode, Codex supervises from chat.
7. In AFK mode, if the user explicitly authorizes unattended progress, Codex
   may launch one bounded repo-local supervisor to perform the same control
   cycle until the repo is complete or genuinely blocked.
8. Keep the runtime bounded to one supervisor, one inner loop, and one
   verifier at a time.

This means the user retains visibility through status artifacts and the chat
history, but does not need to manually approve every ordinary next checkpoint.

## Goal

Provide a resumable outer-loop controller for the stabilization workflow so the
user does not need to manually watch every checkpoint.

The loop should:

1. pick the next active or pending stabilization target from the backlog
2. run exactly one bounded executor checkpoint at a time
3. persist machine-readable progress after every checkpoint
4. run an independent reviewer pass at every checkpoint
5. emit an estimated `%` progress line that can be polled by another process or
   by a later Codex turn

This spec is intentionally repo-local so the implementation survives context
compaction and can evolve independently of the installed skill body.

## Scope

This loop is for `god-module-stabilization`, not for generic refactoring.

It is designed around the existing MultiScholaR workflow:

- `tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md`
- `tools/refactor/HANDOVER-*.md`
- the installed `god-module-stabilization` skill
- the repo’s existing gate and wave tooling

## Non-Goals

- exact semantic completion percentages
- replacing the existing refactor scripts
- autonomous broad repo edits without checkpoint boundaries
- claiming a target is complete based only on line count reduction

## Core Model

The loop follows the same broad architecture family as the `bug-hunter` outer
loop:

1. durable loop state
2. append-only checkpoint history
3. one bounded iteration per executor turn
4. resumable execution
5. stop-file support
6. explicit reviewer step
7. machine-readable status command

Unlike `bug-hunter`, the work unit is not a domain scan. The work unit here is
one stabilization checkpoint on one target.

## Files

Default runtime files under the repo root:

- `.god-module-stabilization/loop-state.json`
- `.god-module-stabilization/loop.jsonl`
- `.god-module-stabilization/loop.stop`

These are runtime artifacts, not source-of-truth documentation.

`loop.stop` is a deliberate finish-current-iteration pause marker. The loop
should notice it, let the active iteration complete normally, persist the
checkpoint/reviewer result, and then park in `status=stopped`. Resume should
use the supported clear-stop path instead of editing state by hand.

The bounded supervisor should also support an explicit clear-stop resume path
so unattended lanes can be restarted without a separate manual loop run.

## Loop State

`loop-state.json` is the current snapshot. It should be rewritten after each
checkpoint.

Minimum fields:

```json
{
  "schemaVersion": 2,
  "createdAt": "2026-04-12T00:00:00Z",
  "updatedAt": "2026-04-12T00:00:00Z",
  "repoRoot": "/repo",
  "backlogPath": "/repo/tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md",
  "runtimeDir": "/repo/.god-module-stabilization",
  "status": "in_progress",
  "iteration": 3,
  "maxIterations": 20,
  "currentItemId": "3-proteomics-normalization",
  "items": [
    {
      "id": "3-proteomics-normalization",
      "bucketNumber": 3,
      "title": "Proteomics Normalization",
      "targetPath": "/repo/R/mod_prot_norm.R",
      "handoverPath": "/repo/tools/refactor/HANDOVER-norm-server-seams.md",
      "supplementalReviewCommands": [
        ["Rscript", "tools/test_with_renv.R", "tests/testthat/test-prot-07-da-analysis.R"]
      ],
      "status": "in_progress",
      "startedAt": "2026-04-12T00:00:00Z",
      "completedAt": null,
      "lastCheckpoint": "late-observer-seam",
      "lastSummary": "Wrapped correlation/export/reset observers.",
      "lastProgress": {},
      "lastReview": {}
    }
  ],
  "history": [],
  "activeRun": {
    "runId": "3-proteomics-normalization-iter-004",
    "loopPid": 12345,
    "itemId": "3-proteomics-normalization",
    "targetPath": "/repo/R/mod_prot_norm.R",
    "startedAt": "2026-04-12T00:00:00Z",
    "heartbeatAt": "2026-04-12T00:00:10Z",
    "phase": "executor",
    "executor": {
      "payloadPath": "/repo/.god-module-stabilization/3-proteomics-normalization-iter-004-executor-payload.json",
      "stdoutPath": "/repo/.god-module-stabilization/3-proteomics-normalization-iter-004-executor.stdout.log",
      "stderrPath": "/repo/.god-module-stabilization/3-proteomics-normalization-iter-004-executor.stderr.log",
      "pid": 12346,
      "startedAt": "2026-04-12T00:00:01Z",
      "completedAt": null,
      "returncode": null
    },
    "reviewer": {
      "payloadPath": "/repo/.god-module-stabilization/3-proteomics-normalization-iter-004-reviewer-payload.json",
      "stdoutPath": "/repo/.god-module-stabilization/3-proteomics-normalization-iter-004-reviewer.stdout.log",
      "stderrPath": "/repo/.god-module-stabilization/3-proteomics-normalization-iter-004-reviewer.stderr.log",
      "pid": null,
      "startedAt": null,
      "completedAt": null,
      "returncode": null
    }
  }
}
```

Item status values:

- `pending`
- `in_progress`
- `done`
- `blocked`
- `failed`

Loop status values:

- `in_progress`
- `completed`
- `paused_max_iterations`
- `stopped`
- `blocked`

## Append-Only History

`loop.jsonl` is append-only and records every checkpoint outcome.

Each line should contain:

- timestamp
- iteration
- target id
- executor result
- reviewer result
- fresh progress snapshot

This is for auditability and for status reconstruction if the snapshot becomes
suspect.

## Progress Model

The loop must print an estimated `%`, but the estimate must be honest about its
basis.

The primary progress estimate is structural, not semantic.

### Active Target Estimate

For a target file family such as `R/mod_prot_norm.R`:

- `legacy_top_level_functions`
  - top-level function count remaining in the legacy target file
- `refactored_top_level_functions`
  - top-level function count in sibling extracted/helper files matching the
    target family prefix

The target estimate is:

```text
target_refactored_pct_estimate =
  refactored_top_level_functions /
  (refactored_top_level_functions + legacy_top_level_functions)
```

This is only an estimate. It intentionally measures “how much of the function
surface has moved out of the legacy file family” rather than “true semantic
completion”.

### Backlog Estimate

Backlog percentages are bucket-based:

- total backlog buckets
- completed backlog buckets
- active backlog buckets
- pending backlog buckets

Repo estimate:

```text
repo_completed_pct_estimate = completed_buckets / total_buckets
```

Proteomics estimate:

- same calculation, filtered to buckets whose titles contain `Proteomics`

### Shape Metrics

Status output must also include current shape metrics:

- legacy file lines
- legacy file max function lines
- legacy file observers count
- legacy file renderers count
- helper file count
- helper function count

## Reviewer Fidelity Contract

The reviewer should prefer typed replay commands:

- `verification.replayCommands[].argv`

The loop may also supply bucket-level supplemental reviewer commands through
item state or target overrides:

- `supplementalReviewCommands`

These should be merged into reviewer replay so fidelity checks do not depend
only on executor-reported commands. This is the preferred way to enforce extra
golden-master, exported-helper, or bucket-specific regression checks.

These are often more informative than the single `%`.

## Status Emission

Every successful loop iteration should emit one single-line progress marker:

```text
GODMOD_PROGRESS target=R/mod_prot_norm.R target_pct=78 repo_pct=18 proteomics_pct=31 refactored_functions=42 legacy_functions=21 review=approved checkpoint=late-observer-seam
```

This line is for terminal monitoring only.

The source of truth remains the JSON state and status command.

## Status Command

The status command should return machine-readable JSON that includes:

- current loop summary
- current item
- current runtime/liveness snapshot
- any stale-run reconciliation event performed during the status call
- latest progress snapshot
- bucket counts
- progress marker line

If no loop state exists, `stabilization-loop.py status` should return a safe
JSON error rather than a traceback. Direct target snapshots remain available via
`stabilization-status.py`.

If the state contains an `activeRun` whose loop PID and child PID are both dead,
or whose heartbeat has expired, `status` should reconcile that stale run into a
terminal blocked transition and report that reconciliation explicitly.

## Executor Contract

One executor iteration equals one bounded stabilization checkpoint.

The executor prompt must instruct the agent to:

- use `$god-module-stabilization`
- read backlog and handover first
- do exactly one bounded seam, staged wave, or equivalent checkpoint
- rerun the focused gate
- update handover/backlog
- stop at the checkpoint
- return structured JSON only

Executor response schema:

```json
{
  "status": "completed",
  "target": "/repo/R/mod_prot_norm.R",
  "checkpoint": "late-observer-seam",
  "targetStatusAfterIteration": "in_progress",
  "summary": "Wrapped the remaining late observers.",
  "reasonCode": "completed_checkpoint",
  "verification": {
    "replayCommands": [
      {
        "argv": [
          "Rscript",
          "tools/test_with_renv.R",
          "tests/testthat/test-prot-05b-norm-module-contracts.R"
        ],
        "label": "focused normalization gate",
        "phase": "post"
      }
    ],
    "display": [
      "focused normalization gate passed"
    ]
  },
  "filesChanged": [
    "R/mod_prot_norm.R",
    "tools/refactor/HANDOVER-norm-server-seams.md"
  ],
  "notes": []
}
```

`targetStatusAfterIteration` values:

- `in_progress`
- `done`
- `blocked`

This field decides whether the loop advances to the next backlog item or stays
on the same target.

## Reviewer Contract

Each executor checkpoint must be followed by a reviewer pass.

The reviewer must be independent from the executor. In the MVP, that means:

- rerun each command reported by `verification.replayCommands`
- compute a fresh status snapshot for the target
- check that claimed files exist
- reject a claimed `done` state if the status snapshot still clearly shows
  `needs-seam-introduction`
- assign an explicit reason code when blocking or requesting changes

Reviewer response schema:

```json
{
  "status": "approved",
  "summary": "Checkpoint verified and gate rerun stayed green.",
  "reasonCode": "approved",
  "testsRun": [
    "Rscript tools/test_with_renv.R tests/testthat/test-prot-05b-norm-module-contracts.R"
  ],
  "verification": {
    "replayedCommands": [
      {
        "argv": [
          "Rscript",
          "tools/test_with_renv.R",
          "tests/testthat/test-prot-05b-norm-module-contracts.R"
        ],
        "display": "Rscript tools/test_with_renv.R tests/testthat/test-prot-05b-norm-module-contracts.R",
        "label": "focused normalization gate",
        "phase": "post",
        "source": "structured"
      }
    ],
    "legacyFallbackUsed": false
  },
  "issues": [],
  "notes": []
}
```

Reviewer status values:

- `approved`
- `needs_changes`
- `blocked`

If review is not approved, the loop must stop and mark the item blocked.
If a blocked review later proves to be a false positive, `recover` may rerun
the stored reviewer payload and append an explicit recovery event instead of
requiring manual state edits.

## Queue Construction

Initial queue is derived from the backlog:

- completed buckets are imported as `done`
- the bucket containing `next active target` starts as `in_progress`
- remaining non-completed buckets start as `pending`

The loop uses the bucket’s first listed file as the primary target for the MVP.

This is acceptable for the current backlog because the active proteomics
targets already list the controlling wrapper file first.

## Stop Conditions

The loop stops when:

- the stop file exists
- max iterations are reached
- the reviewer returns `needs_changes` or `blocked`
- the executor returns `blocked` or `failed`
- all queued items are `done`

## Polling Model

The loop writes state and prints progress markers.

Another process, agent, or later Codex turn should poll via:

```bash
python3 tools/refactor/stabilization-loop.py status
python3 tools/refactor/stabilization-status.py --target R/mod_prot_norm.R --json
```

This is intentionally pull-based. The chat session should not rely on
unsolicited background messages.

## MVP Implementation

The repo-local MVP consists of:

1. `tools/refactor/stabilization-status.py`
   - live target and backlog progress estimator
2. `tools/refactor/stabilization-reviewer.py`
   - local reviewer that reruns reported test commands
3. `tools/refactor/stabilization-codex-runner.py`
   - default noninteractive `codex exec` executor backend
4. `tools/refactor/stabilization-loop.py`
   - loop controller with state, JSONL history, stop file, executor hook,
     reviewer step, and status command

The first implementation may use a runner hook rather than a fully native Codex
SDK backend. That is acceptable as long as the state model, status model, and
review flow are real.

## Recommended Next Step After MVP

After the MVP is stable:

1. add a native `codex exec` backend to the loop controller
2. add a second-agent reviewer backend in addition to the local reviewer
3. surface the progress marker in a small watcher command or tmux-friendly view
4. optionally add line-mass-weighted backlog estimates next to bucket counts

## Operational Gaps

The repo-local MVP is real, but it is not yet safe to treat as a fully
autonomous queue runner for live repository work.

The following work still needs to be implemented explicitly before the loop
should drive real backlog slices without manual supervision:

1. readiness gate before queueing work
   - do not trust backlog `done`/`in_progress` markers by themselves
   - recalculate readiness from live code shape, handover state, and git state
   - refuse to auto-queue a target if the live wrapper still classifies as
     `high-risk-wrapper` / `needs-seam-introduction`
2. backlog/live drift detection
   - flag when handover or backlog says a slice is complete but the current live
     classifier disagrees
   - import is the concrete example discovered during this session:
     `mod_prot_import.R` is a breadcrumb stub, but
     `mod_prot_import_server.R` still classifies as a high-risk wrapper
3. one real end-to-end live loop run
   - the fixture harness is green, but the controller still needs a verified run
     against an actual repository target using the real `codex exec` runner
4. reviewer result classification improvements
   - distinguish environment failures from code regressions
   - for example, missing optional packages such as `mixOmics` should be called
     out separately from behavior failures
5. operator-facing watch path
   - add a tiny `watch` or polling helper so the loop status can be monitored
     without manually retyping the status command

Implemented on April 12, 2026:

- explicit `--target` and `--item-id` override support in
  [stabilization-loop.py](/home/doktersmol/Documents/MultiScholaR/tools/refactor/stabilization-loop.py:1)
- safe `status` output when no loop state exists yet, so polling can happen
  before the first run without a traceback
- network/DNS executor failures are now classified as blocked infrastructure
  runs rather than target failures
- reviewer replay hardening in
  [stabilization-reviewer.py](/home/doktersmol/Documents/MultiScholaR/tools/refactor/stabilization-reviewer.py:1)
  now strips checkpoint prose suffixes, parses commands with `shlex`, reruns
  direct argv instead of `bash -lc`, rejects shell-control operators, and
  deduplicates repeated replay commands
- supervisor hold-state hardening in
  [stabilization-supervisor.py](/home/doktersmol/Documents/MultiScholaR/tools/refactor/stabilization-supervisor.py:1)
  now treats both `blocked` and `stopped` as deliberate hold states instead of
  relaunch triggers
- loop prompt/path hardening in
  [stabilization-loop.py](/home/doktersmol/Documents/MultiScholaR/tools/refactor/stabilization-loop.py:1)
  now resolves runtime paths relative to `--repo-root`, supports
  `--clear-stop-file`, and explicitly forbids inner-loop edits to control
  artifacts like `loop-state.json`, `loop.jsonl`, and `loop.stop`

Current live-run note:

- the loop and bounded supervisor have both now been exercised against the live
  proteomics design bucket
- the earlier false reviewer block on `iter-052` was cleared via
  `stabilization-loop.py recover`, which reran the stored reviewer payload and
  restored the bucket to `in_progress`
- the executor schema was then tightened to an OpenAI-compatible nested-object
  form for structured replay commands
- `iter-055` completed successfully against the live design bucket using the new
  structured replay contract
- overnight April 12-13, 2026 the bounded unattended supervisor then remained
  stable across many sequential live design iterations, auto-extended the loop
  cap from `125` to `150`, and continued launching bounded checkpoints without
  blocked/stopped thrash
- current live state at last check is still the design bucket:
  `in_progress`, runtime `running`, iteration `126 / 150`, run id
  `5-proteomics-design-and-builder-iter-127`
- the latest successful checkpoint is now
  `ui_default_fallback_preview_guidance_review_checkpoint`
- the focused design gate stayed green and the reviewer approved the checkpoint
  using `verification.replayCommands[].argv` rather than legacy `testsRun`
- the main remaining drag is no longer harness stability; it is conservative
  late-stage closeout in the sibling
  [R/mod_prot_design_builder.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_design_builder.R:1),
  which is now a large review-stage file rather than a reactive/render
  god-wrapper

Liveness and recovery status as of April 12, 2026:

Implemented:

1. durable executor liveness tracking
   - `loop-state.json` now records loop PID, child PID, and heartbeat metadata
     in `activeRun`
2. stale-run reconciliation in `status`
   - `python3 tools/refactor/stabilization-loop.py status` now detects a stale
     `activeRun`, writes a blocked recovery transition, and returns a
     `reconciled` event in the response
3. durable executor log capture
   - executor and reviewer stdout/stderr are written under
     `.god-module-stabilization/`
4. explicit recovery command
   - `python3 tools/refactor/stabilization-loop.py recover` now performs the
     same stale-run reconciliation path directly

Remaining practical caveat:

- real `codex exec` runs still depend on backend network access; infrastructure
  outages correctly land as blocked runs with durable logs, but still block
  autonomous progress
- legacy `testsRun` fallback still exists for old stored payloads and recovery
  cases, but new live iterations should emit structured
  `verification.replayCommands` with `argv`, `label`, and `phase`

At this stage, the loop is best treated as:

- a verified checkpoint/status framework
- a reviewer-backed executor harness
- a bounded supervisor harness that can hold correctly on `blocked` and
  `stopped`
- not yet a blind backlog autopilot until the remaining replay-suffix
  normalization gap is closed

Practical update as of April 13, 2026:

- the replay-suffix gap has now been closed and live structured replay is
  working
- the next late-stage failure was a genuine `blocked_executor_timeout` on the
  design bucket, but the failure mode was control-plane drift, not test or
  code breakage
- the bucket target remained `R/mod_prot_design.R` while the real remaining
  structural tail had moved into `R/mod_prot_design_builder.R`
- the timeout run spent its turn rediscovering backlog/handover context instead
  of taking one bounded seam, which is exactly the failure mode the new
  hardening now targets
- the loop now has a repo-local target override layer in
  `tools/refactor/stabilization-target-overrides.json`
  so each bucket can carry:
  - `handoverPath`
  - `workTargetPath`
  - `focusedGateCommands`
  - `promptHints`
- the replay contract now distinguishes replay-safe verification from one-shot
  operational commands:
  - `verification.replayCommands[].replayMode = "verify"` means the reviewer
    should rerun it
  - `verification.replayCommands[].replayMode = "operate"` means the reviewer
    should record it but never rerun it
  - when `replayMode` is omitted, the reviewer infers it from `phase` and
    command content, treating `apply`, `stage`, `commit`, and similar
    operations as non-replay-safe
- the loop prompt now emits a compact checkpoint brief with:
  - latest checkpoint
  - latest reason code
  - latest summary
  - bucket progress
  - optional work-file progress
  - focused gate replay commands
  - a trimmed handover excerpt
- this keeps bucket identity stable while aiming the executor at the real
  structural file and reduces the risk of giant late-stage rediscovery turns
- this also prevents false reviewer blocks after a successful apply checkpoint,
  because the reviewer no longer tries to rerun manifest apply/stage commands
- live wave applies now need a stronger safety contract for shared helper
  targets:
  - apply must preserve any pre-existing top-level target symbols that are not
    part of the current wave
  - apply should be transactional and restore touched files automatically if
    post-apply verification fails
  - the repo-local `tools/refactor/apply_wave.py` now enforces that contract by
    backing up touched files, calling
    `extract_blocks.R --preserve-existing-targets`, validating preserved
    target symbols, and rolling back on failure
- staged-wave collate artifacts now have a matching path contract too:
  - `extract_blocks.R --emit-collate` treats paths with directory components as
    repo-relative outputs instead of nesting them under `--output-root`
  - this prevents ghost staging paths that previously caused false reviewer
    `needs_changes_missing_changed_file` holds
- reviewer changed-file validation is now stricter about real repo files and
  looser about generated staging bookkeeping:
  - missing generated `tools/refactor/staging/**/collate-*.txt` artifacts are
    recorded as notes
  - they no longer fail an otherwise green staged checkpoint
- blocked-state recovery now covers one more false-positive class:
  - `stabilization-loop.py recover` can reopen reviewer `needs_changes` holds
    when a rerun of the stored reviewer payload approves the checkpoint
- module-boundary commits still have not happened for design because the
  sibling builder file remains open

## Current Recommendation

Use the loop now for:

- status polling
- one-iteration bounded runs on explicitly chosen targets
- reviewer-backed checkpoint verification

Do not use the loop yet for:

- automatic selection of the next real target from backlog text alone
- auto-closing slices whose docs say `done` without rechecking the live wrapper
- unattended continuation across a reviewer replay false positive without first
  checking whether the blocked state came from malformed replay commands

Use the loop override layer for:

- late-stage buckets where the public wrapper is mostly frozen but the
  remaining structural work has moved into a sibling helper/builder file
- buckets whose backlog text does not carry a reliable handover link
- buckets that need explicit focused gate replay commands in the executor brief

Practical update as of April 14, 2026:

- the first overnight dual-lane run exposed two separate control-plane issues,
  neither of which was a target-code regression:
  - the main proteomics annotation lane failed one executor launch during
    session startup with `Failed to create session: Read-only file system`
  - the lipid worktree lane completed `func_lipid_qc` successfully, then
    drifted into a proteomics bucket because generic backlog fallback had no
    lane-family filter
- the runner/session issue is now hardened in
  `tools/refactor/stabilization-codex-runner.py` by:
  - building a sanitized child environment for `codex exec`
  - removing inherited `CODEX_THREAD_ID`
  - removing inherited `CODEX_SANDBOX_NETWORK_DISABLED`
  - launching child sessions from a temp cwd rather than the repo root
  - retrying once on known session-init failure markers
- the lane-drift issue is now hardened in
  `tools/refactor/stabilization-loop.py` and
  `tools/refactor/stabilization-supervisor.py` by:
  - adding persistent `scopePrefixes` to loop state
  - filtering `next_open_item()` through scope membership
  - forwarding `--scope-prefix` from the supervisor into loop launches
- a follow-on usability gap is now hardened too:
  - the lipid lane could still park truthfully as `completed` after exhausting
    just two manually seeded items because there was no first-class way to add
    more scoped targets into an existing lane state
  - `stabilization-loop.py queue` now seeds additional manual targets into an
    existing state file while preserving lane scope and resumability
- as a result, the current live topology is now intentionally scoped:
  - main worktree: proteomics/peptide families only
  - lipid worktree: lipid families only
- this is now the required rule for any future parallel lanes:
  - every lane must have its own worktree
  - every lane must have its own `.god-module-stabilization/` state
  - every lane must declare explicit scope prefixes before unattended fallback
    selection is trusted
  - if a lane starts from manual targets rather than backlog-native items, it
    must either be seeded upfront or extended later with
    `stabilization-loop.py queue` so completion state reflects real workflow
    exhaustion rather than a tiny initial seed set
