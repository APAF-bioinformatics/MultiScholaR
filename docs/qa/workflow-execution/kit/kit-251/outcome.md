# KIT-251 outcome: one click on Protein Replicate Filter is handled once

Ticket `tickets/kit/KIT-251.md` r1.2 (plan audit `AUDIT-KIT-251-plan-20261125`, Re-check 1: may start now).
Branch `kit-251` from `v0.9-staging` 16112917329a63ab062655f04db5bf1dd7df9729. No R was run on the laptop; every
run below is a GitHub-hosted run-tests run.

## Commits

| Commit | What |
|---|---|
| c0e1ee89 | G0 script `kit-251/g0-base-values.R` (staging + this file only; nothing under `R/`) |
| 1a60a4c7 | G0 records (`g0-45593374307/`) |
| 6155d2a0 | tests only: new `tests/testthat/test-kit251-replicate-filter-one-click.R` (E1, E2, E3, G1); the declared edit of `tests/testthat/test-perf107-audit.R` (the comment at `:545-546` and `in_process_calls = 2L` at `:574` deleted) |
| 26bc37fb | the fix: `R/mod_prot_qc_protein_replicate.R` `:310-321` at 16112917 deleted (the root-session observer of the Apply input and the blank line before it) |
| (records commit) | records, NEWS, ledger row `VAL-<merge>` |

## Runs (`runs.tsv`)

### G0, run 45593374307, commit c0e1ee89 (base values, observation only)

Instrument: KIT-233's journey driver (`journey.R` at cce7b04c, sourced), network namespace with loopback only, test
mode, the 8-run `prot_dia` fixture, steps S1 to S5.4 by the driver, then ONE real mouse click (chromote
`Input.dispatchMouseEvent`) on Protein Replicate Filter Apply; the lane's permissive retry for that tab removed.

OBSERVED (`g0-45593374307/g0-base-values.tsv`):

| | (a) managed: `entry.R`, plain launcher (`MSR-LAUNCHER {"launcher":"process"}`) | (b) in process: `run_app(test_mode = TRUE)` |
|---|---|---|
| driver mouse clicks; element at the point is the button | 1; yes | 1; yes |
| click events on the button (trusted) | 1 (1) | 1 (1) |
| `shiny:inputchanged` for the button's id | 1 (value 1) | 1 (value 1) |
| worker requests for `protein_replicate_filter` (`request.rds` under the project's `.workflow-export/`) | **1** (12 requests in all, one per operation; `requests-A.tsv`) | 0 |
| notifications from the click on | "Applying protein replicate filter ..." (default); **"A protein replicate filter action is already running" (warning)**; "protein replicate filter applied successfully" (message) | "Applying protein replicate filter..."; "Protein replicate filter applied successfully. Pipeline complete!" |
| step outcome; head (digest) | complete; `protein_replicate_filtered` | complete; `protein_replicate_filtered` |

Stop rule 1 (two worker requests from one click) does **not** fire. The R-log line counts ("observed via module
input" / "via root input") read 0 on both sessions: they are **not values**. The application's logger writes to the
in-app log panel, not to the console (`R/utils_shiny_logging.R:62-70`), so the instrument could not see them. The
managed session's warning is the direct evidence that the second observer ran; in process the second call is silent
(`already_complete`). By VAL-54 a run that produced a result is not re-run; the coordinator accepted G0 as reported.

### R1, run 45596850061, commit 6155d2a0 (tests only; failing first)

OBSERVED (`r1-45596850061/`): `test-kit251-replicate-filter-one-click.R` 6 passed, 3 failed; `test-perf107-audit.R`
190 passed, 1 failed, 1 skipped. Each failure at its predicted value:

| Row | Expected (head) | Actual on the tests-only commit |
|---|---|---|
| E1 managed: `apply` calls | 1 | 2 |
| E2 managed, normalised head: refusal notifications | 1 | 2 |
| E3 in process: in-process apply path calls | 1 | 2 |
| declared edit, `test-perf107-audit.R` replicate module in-process count | 1 | 2 |
| G1 in process: step once, head `protein_replicate_filtered`, success notice once | passes | passed |

### R2, run 45599010969, fix 26bc37fb against tests-only 6155d2a0

OBSERVED (`r2-45599010969/compare.tsv`, `check-compare.txt`): 10 files `same` (identical counts on both sides:
`test-workflow-memory-mechanism-member.R`, `test-workflow-checkpoint-dia-state.R`, `test-workflow-gui-memory-repairs.R`,
`test-prot-norm-input-current.R`, `test-kit200-prot-da-no-value-proteins.R`, `test-e2e-browser-harness.R`,
`test-prot-02c-qc-protein-module-contracts.R`, `test-prot-02h-qc-protein-replicate-shared.R`,
`test-prot-02p-qc-protein-orchestrator-shared.R`, `test-module-ci-prot-protein-qc.R`); 2 files `fixed`, as predicted:
`test-perf107-audit.R` (190/1 failed on base, 191/0 on head) and `test-kit251-replicate-filter-one-click.R` (6/3 on
base, 9/0 on head). R CMD check on both sides: 0 ERRORs, 3 WARNINGs, 3 NOTEs, findings identical (none only on head,
none only on base). Base is the tests-only commit, so the new file shows `fixed`, not `only-on-head` (predicted).
Recounted grep at the branch point: 11 files, as the ticket.

## Acceptance criteria

| # | Criterion | Result | Evidence |
|---|---|---|---|
| A1 | G0 base values (a) managed and (b) in process | met | G0 above; `g0-45593374307/g0-base-values.tsv` |
| A2 | E1, E2, E3 and the declared edit fail at G0's values on the tests-only commit; G1 passes; all pass on the fix | met | R1 above; R2 |
| A3 | R2 base against head: every file `same` except the new file and the declared edit | met | R2 `compare.tsv` |
| A4 | NEWS, ledger, this outcome with "For the change report", the browser step handed to KIT-233 J | met | NEWS (Fixes); ledger `VAL-<merge>`; `tickets/kit/KIT-233.md` addendum "2026-11-25, KIT-251 hand-off" |

## Core-hours

G0 20.2 + R1 11.8 + R2 139.2 = 171.2 job-minutes x 4 / 60 = **11.41 core-hours**, against 23.33 expected and 102.0
ceiling. No correction pass, no instrument rerun.

## Inferred (not executed)

- In process on base, the second call returned `already_complete` (`R/mod_prot_qc_protein_replicate.R:143-146`):
  inferred from one success notice and one step completion; pinned by R1's E3 and `test-perf107-audit.R`.
- The fix cannot lose a click: the removed listener read the same reactive input as the remaining one (ticket
  section 2); G0 shows one browser input change per click, which the module observer receives.

## Not checked

- R-log observer lines in the running application (not on the console; see G0).
- A click on head in a browser (handed to KIT-233's journey J; stop rule 2).
- Whether a first-call failure before `R/utils_workflow_execution.R:945` let the second call start a second request
  on base (not a normal path).
- The TMT and label-free e2e lanes (same button, helper and module; stated, not run).
- No full suite (the merge train's).

## For the change report (appendix: code new to 0.9, not in the legacy pipeline)

**What was wrong, who was affected.** In the 0.9 development line the proteomics Protein Replicate Filter tab listened
to its Apply button twice: through the module's own input and through the root session's copy of the same input (an
observer added with the DIA browser test lane, 132e1359, 2026-05-05; `main` has one observer,
`main:R/mod_prot_qc_protein_replicate.R:85`). Every proteomics user of the installed application on the DIA (iq),
label-free and TMT routes saw, at every Apply, the warning "A protein replicate filter action is already running"
beside the working notice, although nothing else was running; on a normalised table KIT-182's refusal was shown
twice. In process (no managed worker) the second call was skipped silently.

**Code.** Base 16112917: `R/mod_prot_qc_protein_replicate.R:303-309` (module observer, kept) and `:311-321` (root
observer, removed); the guard that produced the warning `R/utils_workflow_execution.R:896-898`. Fix 26bc37fb: the
module observer alone (`R/mod_prot_qc_protein_replicate.R:303-309`).

**What changed in results.** Nothing: the step ran once per click before and after (G0: one worker request from one
real click under the plain launcher; in process one step run). One false warning disappears; a refusal is shown once.

**Evidence.** G0 45593374307 (one real mouse click in headless Chrome against the application started through
`packaging/electron/entry.R`, network closed, the 8-run fixture; and the same in process); R1 45596850061 (E1 to E3 and
`test-perf107-audit.R` at 2 on the base, G1 passing); R2 45599010969 (base against head). Limits: synthetic 8-run fixture;
the browser check on head is KIT-233 J's step.
