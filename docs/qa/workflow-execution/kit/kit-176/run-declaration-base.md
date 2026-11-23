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

Added after the implementation audit (finding KIT-176-I03): T3's step is the workflow copy in 368f9fab on `run-k176-base` only; it is not merged because its first pair needs 85770930, reachable only from `run-k158-install`. Both run branches and `run-k158-install` are kept so that the step, the declarations and the pairs stay reachable.

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

## Results (read from the job logs 2026-11-23; no artifact, no cache, the token never printed)

R1 [45112185913](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/45112185913) on 001e8e4d and
R2 [45112365197](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/45112365197) on e3717079 (declared in
`run-declaration-head.md`). Per job: `runs.tsv`; families: `lane-families.tsv`; stamps: `stamp.json`; fuses: `fuse-reads.txt`.

- **Failing first (A2), as declared on all three systems.** Pass 3 on the base: no `start` event, `FAIL the shell started`,
  `FAIL the window loaded ...` with exit 0 (ubuntu-22.04, windows-2025) and exit 9 (ubuntu-24.04, Node rejecting
  `--no-sandbox`), then `TypeError: Invalid URL`. **Windows is affected** on the base, so the NEWS bullet keeps
  "(Linux and Windows.)". Pass 1 failed only on T2's three checks; pass 2 passed. `RunAsNode is Enabled` on each binary.
  T4 failed on each (installed field `NA`). T3: every call as declared on all three jobs (pairs; shallow stop; dirty
  stops naming ` M R/app_server.R` and `?? inst/kit176-t3-untracked.txt`).
- **Head (A3, A4, A7).** The rule step on each job: `package commit b052b824...; on the release line (v0.9-staging):
  true; build commit e3717079...` (the release line had moved to 271cd3c4 by then, a ledger-only descendant of b052b824;
  the merge base stayed b052b824). All three lifecycle passes passed on all three jobs; pass 3 recorded
  `runAsNodeVariable "1"`; `RunAsNode is Disabled` on each binary; T2's About lines `Analysis code: commit b052b824...` /
  `Test build from: commit e3717079...` matched; T4 passed on each.
- **B13.** sha256 of `packaging/package-commit.sh` 30c217ce... in R1's T3 step and in R2's rule step on all six jobs: equal.
- **A5.** Every other step's outcome identical, job by job; the GUI-entry difference blocks identical (290, 290 and 65
  families with equal counts, no family in one run only, so none involves the software block and D4 stands); the e2e
  tables identical (ubuntu-22.04 7 lanes, windows-2025 4 lanes, all at the laptop's counts).
- **Checkout** with `fetch-depth: 0`: 6 to 14 s against 3 to 9 s in run 44998302733.
- **Budget.** R1 197 job-minutes (13.1 core-hours), R2 209 (13.9); 27.1 core-hours of the 180 ceiling (expected 26.4;
  the R2 Windows job took 120 min against R1's 86, inside its timeout). No correction pass used.
- **Unrelated discovery, left untouched.** `lockfile_sha256` in `VERSION.json` is 483c964b... on Linux and b911b381...
  on Windows in both runs (the base too): on Windows the lock file is hashed as checked out (line endings, presumably).
