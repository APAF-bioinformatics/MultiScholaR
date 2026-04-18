# Outer Supervision Protocol

This file is the explicit reminder for future Codex turns.

## Authority Boundary

1. Default mode is Codex-in-chat as the outer supervisor.
2. If the user explicitly authorizes unattended progress, Codex may launch one
   repo-local bounded supervisor that watches the inner loop, decides ordinary
   go/no-go transitions, and commits at module boundaries.
3. Do not run more than one repo-local supervisor at a time.
4. Keep the process budget bounded:
   - one supervisor
   - one inner loop
   - one verifier when active
5. Do not spin up additional autonomous drivers outside that bounded set.

## Required Cycle

For each module or active target:

1. poll the inner loop state
2. inspect the checkpoint result and changed files
3. verify fidelity with the focused gate and any obvious sanity checks
4. decide go/no-go for the next bounded inner loop
5. launch the next bounded inner loop if safe
6. commit only when a real module boundary is reached

When AFK/autonomous mode is enabled, the bounded local supervisor performs this
cycle on Codex's behalf and records its decisions in repo-local state/logs.

## Stop Rule

Keep going until:

- the module boundary is complete and committed, or
- there is a genuine blocker that requires strategic elevation, redesign, or
  explicit user direction

Do not stop for ordinary gate repair, routine seam work, manifest review/apply,
or verification reruns.

## Pause And Resume Rule

- Creating `.god-module-stabilization/loop.stop` means: let the current
  bounded iteration finish normally, persist the checkpoint/reviewer result,
  then park the lane in `status=stopped` / `runtimeStatus=idle`.
- A stopped lane is an intentional hold state, not a failed lane.
- Resume a stopped lane with the supported clear-stop flow before or alongside
  restarting the supervisor:
  - `python3 tools/refactor/stabilization-loop.py run --clear-stop-file ...`
  - `python3 tools/refactor/stabilization-supervisor.py run --clear-stop-file`

## Commit Rule

At module boundaries:

- commit production files only
- prefer family-scoped `R/` files and `DESCRIPTION` when needed
- leave docs/tests/tooling uncommitted unless explicitly requested

## Current Focus

- Proteomics lane is now genuinely complete for the current omics scope.
- Lipidomics lane is now genuinely complete for the current omics scope.
- Metabolomics is the only active omics lane still open.
- `origin/janitor` now contains the merged lipidomics + proteomics production
  history via:
  - `aab75ec` `Merge local proteomics janitor history into shared janitor branch`
- Key recent proteomics closeout commits included in that push:
  - `946f546` `Stabilize mod_prot_summary`
  - `e924175` `Finalize residual proteomics S4 helper extractions`
- Key lipidomics closeout commit already included in that push chain:
  - `2ddbde5` `Finalize residual lipidomics helper extractions after lane completion`
- Current omics completion estimate:
  - proteomics: complete for the scoped omics lane
  - lipidomics: complete for the scoped omics lane
  - metabolomics: still in progress
- Current active lane at this documentation update:
  - repo root: `/home/doktersmol/Documents/MultiScholaR-metab-safe`
  - target: `R/mod_metab_import.R`
  - lane state: clean blocked hold
  - current status snapshot: `blocked`, `runtimeStatus=idle`
  - latest completed checkpoint:
    `metab-import-module-annotation-status-custom-render-registration-seam`
  - blocker class: false reviewer replay failure from malformed structured
    `verification.replayCommands[].argv`, not a bad metabolomics refactor

## Parallelism Policy

- Same-worktree parallel loops remain unsafe.
- Parallel stabilization is allowed only in separate git worktrees with
  isolated `.god-module-stabilization/` state per lane.
- Shared merge points must remain serialized:
  - `DESCRIPTION`
  - backlog and handover docs under `tools/refactor/`
  - final production commits onto the main branch
- `dev/test_lipid_*` scripts remain optional manual-dev helpers only and are
  not a blocker for mirrored lipid/metabolite waves.
- Each parallel lane must carry explicit scope prefixes so it cannot drift into
  another omics family when it falls back to generic backlog item selection.
- Each lane should also run the hardened Codex runner that isolates child
  session state from the parent chat environment before launching `codex exec`.
- Each lane should also be able to inject bucket-level supplemental reviewer
  commands so fidelity checks do not depend only on executor-reported replay
  commands.

## Next Candidate Lanes

- Immediate next operational step:
  - recover and relaunch metabolomics from the current blocked reviewer hold
- Recovery note for the next session:
  - the replay-contract fix is already implemented and verified in the main
    repo tooling
  - regression suite passed `30/30`
  - what still needs doing is to sync the patched
    `stabilization-loop.py` / `stabilization-reviewer.py` into the metabolomics
    worktree, run `recover`, and restart the metabolomics supervisor
- After metabolomics closes:
  - only the general shared-code lane should remain
  - likely targets:
    - `R/func_general_filemgmt.R`
    - `R/func_general_plotting.R`

## Live Lane Setup

- Main proteomics lane:
  - repo root: `/home/doktersmol/Documents/MultiScholaR`
  - lane state: `completed`, `runtimeStatus=idle`
  - scope prefixes:
    - `R/func_prot_`
    - `R/mod_prot_`
    - `R/func_pept_`
- Parallel lipidomics lane:
  - repo root: `/home/doktersmol/Documents/MultiScholaR-lipid-lane`
  - branch: `stabilization-lipid-lane`
  - lane state: `completed`, `runtimeStatus=idle`
  - scope prefixes:
    - `R/func_lipid_`
    - `R/mod_lipid_`
- Parallel metabolomics lane:
  - repo root: `/home/doktersmol/Documents/MultiScholaR-metab-safe`
  - lane state: `blocked`, `runtimeStatus=idle`
  - active target at hold: `R/mod_metab_import.R`
  - current item id: `manual-mod-metab-import`
  - scope prefixes:
    - `R/func_metab_`
    - `R/mod_metab_`
- Worktree note:
  - sibling worktrees do not inherit repo-local untracked `tools/`, so the
    current local `tools/` directory must be copied into each worktree when the
    harness changes
- Recent hardening that made the current multi-lane topology viable:
  - `stabilization-codex-runner.py` now launches child `codex exec` from an
    isolated temp cwd with a sanitized environment, avoiding parent-session
    contamination and the observed read-only session-init failure
  - `stabilization-loop.py` and `stabilization-supervisor.py` now persist and
    honor lane `scopePrefixes`, preventing a lipid lane from drifting into a
    proteomics bucket after completing a manual target
  - `stabilization-loop.py` now also supports explicit queue seeding so a
    scoped lane can continue through the rest of a workflow instead of
    reporting `completed` after exhausting one or two manual seed targets
  - `stabilization-loop.py` and `stabilization-reviewer.py` now support
    `supplementalReviewCommands`, so loop-state overrides can force extra
    bucket-level fidelity checks even if the executor omits them
  - `stabilization-loop.py` and `stabilization-reviewer.py` now defensively
    trim malformed structured replay-command metadata tails when bad executor
    output leaks `label` / `phase` / `replayMode` tokens into `argv`
  - this last hardening is what the metabolomics lane needs in its worktree
    before the current blocked reviewer hold is recovered

## AFK Mode

When the user is away and explicitly authorizes unattended progress:

1. stop any standalone passive poller so only one supervisor remains
2. launch one bounded local supervisor
3. let it poll the loop state, relaunch ordinary next checkpoints, and commit
   only production files at module boundaries
4. leave docs, tests, and tooling uncommitted unless explicitly requested
5. stop only for a genuine strategic blocker that the supervisor cannot repair
