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

- Proteomics design is complete and committed in `77fb6c7`
  (`Stabilize proteomics design and builder wrappers`).
- Proteomics closed buckets are now:
  - DA follow-up
  - import
  - normalization
  - QC/rollup
  - design/builder
- Current completion estimate:
  - proteomics: `4/8` done = `50.0%`
  - repo-wide: `4/13` done = `30.8%`
- Active main-worktree lane at this documentation update:
  - target: `R/mod_prot_enrich.R`
  - lane state: deliberately paused after current iteration
  - current status snapshot: `stopped`, `runtimeStatus=idle`
  - latest completed checkpoint: `prot_enrich_stage_observer_shell_wave11`
- Active parallel lipid worktree lane at this documentation update:
  - target: `R/mod_lipid_qc.R`
  - lane state: deliberately paused after current iteration
  - current status snapshot: `stopped`, `runtimeStatus=idle`
  - latest completed checkpoint: `bounded_wrapper_init_seam`

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

- Recommended serial next bucket in the current worktree:
  - `6. Proteomics Annotation`
- Recommended first parallel lipidomics candidate in a separate worktree:
  - `R/func_lipid_import.R`
- Why `R/func_lipid_import.R` now:
  - the first two parallel proof targets are already complete in the lipid
    worktree:
    - `1e55c0f` `Stabilize func_lipid_qc`
    - `d5add9d` `Stabilize func_lipid_da`
  - the lipid lane now continues through a queued family workflow rather than
    parking after one manual target
  - `func_lipid_import.R` is `direct-extraction-ready` and lower-risk than the
    remaining lipid wrapper modules
- Why not start with `R/mod_lipid_norm.R`:
  - it is still `high-risk-wrapper, needs-seam-introduction`
  - current direct metric is `2107` lines, `7` observers, `8` renderers
  - better kept for a later lane once the parallel pattern is proven

## Live Lane Setup

- Current worktree serial proteomics lane:
  - repo root: `/home/doktersmol/Documents/MultiScholaR`
  - active target: `R/mod_prot_enrich.R`
  - loop state: `.god-module-stabilization/loop-state.json`
  - supervisor state: `.god-module-stabilization/supervisor-state.json`
  - scope prefixes:
    - `R/func_prot_`
    - `R/mod_prot_`
    - `R/func_pept_`
- Parallel lipidomics pilot lane:
  - repo root: `/home/doktersmol/Documents/MultiScholaR-lipid-lane`
  - branch: `stabilization-lipid-lane`
  - active target: `R/mod_lipid_qc.R`
  - loop state: `.god-module-stabilization/loop-state.json`
  - supervisor state: `.god-module-stabilization/supervisor-state.json`
  - scope prefixes:
    - `R/func_lipid_`
    - `R/mod_lipid_`
- Worktree note:
  - the lipid worktree does not inherit `tools/` from git because the refactor
    harness is repo-local and untracked, so the current local `tools/`
    directory was copied into the worktree deliberately
- Recent hardening that made the current dual-lane topology viable:
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

## AFK Mode

When the user is away and explicitly authorizes unattended progress:

1. stop any standalone passive poller so only one supervisor remains
2. launch one bounded local supervisor
3. let it poll the loop state, relaunch ordinary next checkpoints, and commit
   only production files at module boundaries
4. leave docs, tests, and tooling uncommitted unless explicitly requested
5. stop only for a genuine strategic blocker that the supervisor cannot repair
