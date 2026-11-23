# KIT-176 bundle runs (declared 2026-11-23, before any run)

Ticket `tickets/kit/KIT-176.md` revision 1.4, section 11. Two run branches; the workflow's concurrency group is per ref,
so they do not cancel each other.

## Bindings

| | |
|---|---|
| Staging commit (the head is cut from it) | `b052b82441ef2f11c6067420dbeff8e95f8c8505` (`v0.9-staging`) |
| R1 branch | `run-k176-base` = b052b824 + 568c19c2 (the test-only change, also the first commit of `kit-176-launcher`) + the T3 instrumentation commit (this branch's workflow copy only) + this declaration |
| R2 branch | `run-k176-head` = `kit-176-launcher` head + its own declaration commit |
| Workflow | `KIT-134 distributable`, `workflow_dispatch`, stage `full`, platform `both` (ubuntu-22.04, ubuntu-24.04, windows-2025) |
| Read token | `DEPS_READ_TOKEN` only in the library install step, as the existing workflow; no cache; no upload; nothing built leaves a job |
| Declared package commit (R1) | `MSR_PACKAGE_COMMIT` = b052b824 (full sha above), set in this branch's workflow `env:`; T2 (lifecycle pass 1) and T4 expect it |
| T3 pairs | `t3-pairs.tsv` beside this file: (85770930, 22a18b63) gives 22a18b63 true; (22a18b63, 12560b8a) gives 22a18b63 false; (9a53d791, 9a53d791) gives 9a53d791 true. Re-derived before the run with read-only git (merge bases 22a18b63, 12560b8a, 9a53d791, one each; package paths equal, different, equal). 85770930 is reachable only from `origin/run-k158-install`, which exists at this declaration (GitHub API: branch head 85770930d1f74b4b61e936ac0d28d5e18d14aff4); if it is gone by the run, the presence check reports an instrumentation failure |

## Expected (R1, base application plus the test-only change)

- Lifecycle, every job: pass 1 fails only on T2's checks (`VERSION.json names the package commit ...`, `... package_on_release_line is a JSON boolean`, `the About lines ...`); pass 2 passes; pass 3 (`ELECTRON_RUN_AS_NODE=1`) fails with `FAIL the shell started` and `FAIL the window loaded the application and a live session rendered the first tab`, the variable check failing too, then `TypeError: Invalid URL`, exit 1. ubuntu-22.04: exit 0 of the shell; ubuntu-24.04: non-zero (Node rejects `--no-sandbox`). windows-2025: the inference to be tested; if pass 3 prints `ok  the shell started` and `ok  the window loaded ...`, finding 1 does not exist on Windows (recorded; NEWS names Linux only; the fuse stays).
- Fuse read of the installed binary: `RunAsNode` enabled on all three.
- T4: fails on all three (field absent, `NA`).
- T3: the three pairs as declared, the shallow copy stops at the shallow check, the two dirty trees (` M R/app_server.R`, `?? inst/kit176-t3-untracked.txt`) stop at the dirty-tree check, on all three jobs; the sha256 of `packaging/package-commit.sh` printed (compared with R2's rule step, B13).
- Every other step as in run 44998302733; the GUI-entry lanes step fails on its recorded difference families as there; the families are R2's comparison base (A5).

## Expected (R2, head)

All three lifecycle passes pass on all three jobs; the variable check passes in pass 3; `RunAsNode` disabled; `package_commit` = b052b824, `package_on_release_line` true, `commit` = the head branch's commit (different); T4 passes; the rule step's sha256 of `packaging/package-commit.sh` equals R1's; lane difference families and e2e tables identical to R1's, job by job (`trx_screen` repeats as the union).

## Ceiling and attempts

Per run 3 jobs x 300 min x 4 cores = 60 core-hours; expected 13.2 (re-derived from run 44998302733). R1 once, R2 once, one
correction pass in total (repeating R2, or R1 if R1 was the invalid run). An instrumentation failure is recorded as such.
A runner stall (more than 30 minutes in a set-up step) is cancelled and its failed jobs rerun; that is not a correction pass.

Results below this line after the runs.
