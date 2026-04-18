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
- There is no active supervisor or inner loop by default after the design
  commit boundary; the next lane should be chosen deliberately.

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

## Next Candidate Lanes

- Recommended serial next bucket in the current worktree:
  - `6. Proteomics Annotation`
- Recommended first parallel lipidomics candidate in a separate worktree:
  - `R/func_lipid_qc.R`
- Why `R/func_lipid_qc.R` first:
  - `direct-extraction-ready`
  - no `moduleServer()`, observers, or renderers
  - lower orchestration risk than `R/mod_lipid_norm.R`
- Why not start with `R/mod_lipid_norm.R`:
  - it is still `high-risk-wrapper, needs-seam-introduction`
  - current direct metric is `2107` lines, `7` observers, `8` renderers
  - better kept for a later lane once the parallel pattern is proven

## AFK Mode

When the user is away and explicitly authorizes unattended progress:

1. stop any standalone passive poller so only one supervisor remains
2. launch one bounded local supervisor
3. let it poll the loop state, relaunch ordinary next checkpoints, and commit
   only production files at module boundaries
4. leave docs, tests, and tooling uncommitted unless explicitly requested
5. stop only for a genuine strategic blocker that the supervisor cannot repair
