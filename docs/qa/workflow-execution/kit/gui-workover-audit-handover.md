# GUI workover series: handover for a post-hoc audit (2026-11-11)

You are auditing the GUI workover series (KIT-101 to KIT-113) on branch `performance-v2` of MultiScholaR. The lead ran it unattended under the user's instruction of 2026-11-07:

> "continue with this series agentically until the GUI is done. any fixes along the way - do it as long as the scientific methodology doesnt change. anything that needs my input park for end of the chain."

The user wants an independent check of the work, the evidence and the lead's reports before deciding the parked questions.

**Mode:** `post_hoc`, audit only. Use the `audit-evidence-contracts` skill if you have it. Make no repairs, start no runs, and change no ticket status. Write findings to a new file, `docs/qa/workflow-execution/kit/gui-workover-audit-<date>.md`; do not edit the evidence you audit.

## 1. Where things are

| What | Where |
|---|---|
| Repository rules | `AGENTS.md`; `conventions/R.md`, `R-shiny.md`, `R-golem.md`, `R-bioinformatics.md` |
| Tickets | `tickets/kit/KIT-1NN.md` (git-ignored, read them from the main checkout); index `tickets/kit/tickets-index.json` |
| Series plan and design | `tickets/kit/00-KIT-PLAN.md`; `docs/qa/workflow-execution/kit/kit-101/design.md`; scoping `docs/qa/workflow-execution/kit/gui-workover-scoping/` |
| Per-ticket evidence | `docs/qa/workflow-execution/kit/kit-1NN/`: `outcome.md`, `findings.md`, `memory-outcome.md` (where measured), `skip-compare.md`, `entry-suite.md`, `merge-checks.txt`, `lines.tsv`, `proof/`, `w<N>/`, `kit1NN-kit-journey-*-stage.json` |
| Decisions parked for the user | `docs/qa/workflow-execution/kit/gui-workover-parked-decisions.md` (P-01 to P-43; nothing parked there has been acted on, except where a row says "done") |
| Quarantined tests | `docs/qa/workflow-execution/kit/quarantined-tests.md` |
| Run protocol | `docs/qa/workflow-execution/kit/protocol.md` (v1.9); `bindings-generator.R`; `bindings.json` |
| Proof machinery | `tools/profiling/workflow_kit*.R`; `tools/testing/test-baseline.R`; `kit-105/skipcmp.py`; `kit-105/lifecycle/lifecycle.R` |
| Builds and run outputs (outside the repo) | `/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/` (`core/browser-smoke.*`, `kit/runtime/kit-journey-*`, `kit/scripts-1NN/`) |

**Ticket order run:** KIT-101 → 102 → 103 → 104 → 105 → 110 → 106 → 107 → 108 → 111 → 112 → 113. KIT-109 (closure: removals, quarantine manifest, full suite) is a planning-gated stub and has not started. KIT-114 is methodology, parked for the user.

**Builds:**

| Build | Install | Source |
|---|---|---|
| `kit_107` | `66sPmE` | |
| `kit_108` | `X6BzNV` | |
| `kit_111` | `MlQroA` | `c568dfa5` |
| `kit_112` | `OftnRb` | `5ecc9fd1` |
| `kit_113` | `cmBMvS` | `95600989` |

## 2. State at handover

- **KIT-101 to KIT-112:** marked `completed`, with their evidence committed and pushed.
- **KIT-113 is not closed.**
  - Runs 1 to 8 and the memory measurement are done and committed.
  - The ticket is still `in_progress`, with its boxes unticked.
  - `kit-113/outcome.md` still says `RUN9` for the skip comparison.
- **Run 9 was interrupted and restarted.** The check-mode skip run was cut off by a session restart that also cleared the lead's scratch folder.
  - The lead rebuilt the overlaid skip baseline from the committed tables, in ticket order: KIT-105, then 110, 106, 107, 108, 111, 112.
  - The rebuild was checked by reproducing KIT-112's recorded comparison exactly in both modes.
  - Then only the check-mode run was restarted. **Check that rebuild.**
- **The test rule** (user, 2026-10-18): each ticket runs its targeted files plus a dependents grep, `R CMD check` and an install. The full suite runs only at checkpoints, so none ran inside this series; KIT-109 runs it. **Find when the last full suite ran.**

## 3. Standing constraints (do not break them while auditing)

- Never read raw real-data (KC) datasets or their outputs; aggregates and hashes only.
- R is pinned at 4.6.1. Never install into or upgrade the default library.
- Never run `roxygen2::document()`.
- Never use `pgrep -f` / `pkill -f`.
- Start no `multischolar*` systemd units.
- Commits carry no Co-Authored-By or Claude attribution.

## 4. What to assess

Each item below states the claim and what to check. They are the lead's own list of weak points, ordered by how much a wrong answer would matter. Add what you find beyond them.

### A. Scientific record and correctness

1. **F112-05: provenance on the managed DA path** (`kit-112/findings.md`).
   - **Claim:** since KIT-019/025, a managed `da_publication` completion leaves `da_publication` pending as the controller's step provenance, and the next state commit records it as its own.
   - **Check:** whether any recorded analysis record (journeys, oracles) already carries wrong provenance. Assess whether deferring it to KIT-109 was sound.
2. **F112-01: missing notices on the managed DA path.** The BH-fallback warning and the completion notice exist only on the in-process branch. The managed path is the default on Linux. Assess whether a missing BH-fallback warning is a scientific-reporting defect rather than "a notice text".
3. **P-35 / F111-02: repeated RUV covariates.** A second DIA_LIMPA "Apply RUV covariates" click applies RUV to its own output. It was parked, not fixed, as a methodology-adjacent change. Assess that classification.
4. **F112-04: validation of all-NA design terms.** Validation accepts a formula term on an all-NA design column. Recorded, not fixed.
5. **The "identical" claims rest on frozen captures plus recorded exclusions.**
   - KIT-111 to KIT-113 compare worker against in process with listed exclusions: D6 lists, F113-04, F113-06, entity counts, notices and progress.
   - **Check:** that each exclusion is asserted exactly in the tests (not masked loosely), and that no exclusion hides a state, product, provenance or file difference.
   - Tests: `test-workflow-saved-da.R` (KIT-112 W1/W2 blocks), `test-omics-enrich-core.R` and `test-omic-enrich-module.R` (KIT-113 blocks), `test-prot-norm-core.R` (KIT-111).
6. **Re-frozen oracles.** Few were re-frozen; they are attributed in `kit-104/w7/p16-diff.txt`, `kit-111/lane-changes.md` and the `kit-10N/lane-changes.md` files. Check that every re-frozen value is attributed to a named fix, and that no re-freeze absorbed an unexplained change.

### B. Memory and performance claims

7. **KIT-112's "better" claim** (`kit-112/memory-outcome.md`). `da` 2.208 → 1.575 GiB, one run per member.
   - 112q's controller sat 0.15 to 0.37 GiB lower from `import` onwards, in steps KIT-112 does not change. The lead therefore states the claim on the step's own rise (+0.500 → +0.015 GiB above `session_export`).
   - **Check:** that reasoning against the traces (`kit/runtime/kit-journey-11{1,2}q/candidate-1/memory-trace.tsv`, 8 columns, controller anon is the last).
   - **Also unexplained:** early steps of 112q ran 1.5 to 1.7× slower, and the whole journey 1.065×.
8. **KIT-113's placement despite a projected regression** (`kit-113/w0/enrichment-heap.md`, `kit-113/findings.md` W0 decision).
   - The W0 probe projects the managed controller **+0.66 GiB worse** at enrichment. The lead followed D8's declared outcome for "unreachable": move made, P-43 raised.
   - **Assess the lead's judgement:** D8 anticipated "not better"; was a projected "worse" inside that outcome or a gate?
   - **Assess the probe's validity:**
     - same-process emulation of the worker;
     - the warm-pass choice, declared by the probe agent, not the ticket;
     - an order effect the agent reported;
     - the hydration size measured with `lobstr::obj_size`.
9. **KIT-111's correction pass** (`kit-111/memory-outcome.md`). 111m exceeded a "no worse" rule by +0.13 GiB. The rerun 111n passed and was used to hold the rule. Check this was a legitimate use of the one correction pass (a diagnosed one-off release floor), not a retry until pass. Both members are reported.
10. **KIT-108 D8 missed** (P-30). Check the diagnosis and the measured figures.
11. **Host conditions.**
    - During 111q, two read-only agents were working; this is recorded.
    - Check `kit-1NN/power-state.txt` for every measured member, and that no measured window overlapped heavy work.

### C. Test evidence and coverage

12. **Check-mode error blocks are a growing class.** New test blocks that use the DIA checkpoint fixture `readWorkflowCheckpointRecipe` error under the `R CMD check` layout, because `tests/testdata/workflow-checkpoints/scenarios.json` is read from the repository root and is not shipped. That is KIT-111's 11 and KIT-112's 8 (see each `skip-compare.md`). Those blocks therefore never run under check. Assess whether "the known source-tree read" is an acceptable classification, or whether it hides lost coverage.
13. **Skip comparisons.** Every ticket compares its targeted files against an overlaid baseline. Check the overlay order and the filtering (only targeted files are reported), and that no file lost blocks or gained skips unseen.
14. **Coverage gaps the series recorded but did not close:**
    - no run ever reaches enrichment on a managed session (P-40);
    - no lane, journey or e2e renders proteomics' DA tab (F106-01);
    - the managed enrichment path rests only on module tests with an in-process test worker.
15. **Agent-written code.** W1 to W3 of KIT-111, KIT-112 and KIT-113 were implemented by sub-agents in worktrees. The lead verified by test counts, merge checks and spot reads of diffs, not by a line-by-line review.
    - Review the `R/` diffs: KIT-111 `6c04767f..c568dfa5`, KIT-112 `c568dfa5..5ecc9fd1`, KIT-113 `132e5445..95600989`.
    - Check them for duplicated logic, convention breaches beyond the recorded divergences, and behaviour on the unmanaged path.
16. **Lines over estimate:**
    - KIT-111: `mod_prot_norm_workflow_helpers.R` +153;
    - KIT-112: `mod_prot_da_handlers_run_analysis.R` +51;
    - KIT-113: `utils_workflow_enrichment.R` +235 and `mod_omic_enrich.R` +67.

    Assess against the conventions' size rules and the series' aim of removing duplication.

### D. Series goal and closure

17. **Has the series met its goal?** The goal was one shared GUI layer of light reactive modules for every omic (plan: `/home/doktersmol/.claude/plans/write-a-series-of-atomic-dragon.md`, `kit-101/design.md`).
    - Proteomics' normalisation tab is still its own (P-32).
    - The profile hooks and the duplicate lipid/metab copies are removed only in KIT-109.
    - Assess how far the series is from its stated outcome.
18. **Items queued to KIT-109:**
    - F108-01: failures in closures are lost;
    - F108-04: git info depends on the checkout;
    - F111-05: QC image paths point into the staging copy;
    - F112-01, F112-04, F112-05, F112-06;
    - the stale KIT-019 comment at `utils_workflow_operations.R`;
    - the W2/W3 convention divergences;
    - the diagnostic exit boundary (`workflow_checkpoint_gui_diagnostic.R`).

    Check none was lost and none should have been fixed inside its own ticket.
19. **Worktrees and housekeeping.** There are about 59 agent worktrees under `.claude/worktrees/` and about 139 entries in `git worktree list`. Report which are stale. Remove nothing.

## 5. Report back

For each item:
- a finding id;
- severity;
- the observed evidence (file:line or commit);
- the impact on a claim;
- the smallest sufficient repair;
- whether it needs the user's decision.

**Keep separate:**
- implementation verified;
- evidence valid;
- claim supported;
- what the user must decide.

An inconclusive result is acceptable. Say so rather than inferring.
