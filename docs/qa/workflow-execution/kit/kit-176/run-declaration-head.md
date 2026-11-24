# KIT-176 R2 on run-k176-head (declared 2026-11-23, before the run)

Ticket `tickets/kit/KIT-176.md` revision 1.4, section 11; R1 is declared on `run-k176-base` (`run-declaration.md`).

| | |
|---|---|
| Branch | `run-k176-head` = `kit-176-launcher` 93f96917 (568c19c2 test-only change, 5700c8a1 the fix, 93f96917 NEWS; cut from `v0.9-staging` b052b824) + this declaration (a file under `docs/` only) |
| Workflow | `KIT-134 distributable`, `workflow_dispatch`, stage `full`, platform `both`; `DEPS_READ_TOKEN` as the existing workflow; no cache, no upload |
| Declared package commit | `b052b82441ef2f11c6067420dbeff8e95f8c8505` (the staging commit the head was cut from; no commit since touches the package paths), `package_on_release_line` true, `commit` = this branch's head (different) |
| Expected | the rule step: no stop; the sha256 of `packaging/package-commit.sh` equal to R1's T3 line; all three lifecycle passes pass on ubuntu-22.04, ubuntu-24.04 and windows-2025, pass 3's variable check passes; `RunAsNode` disabled on each installed binary; T2 and T4 pass; every other step, the GUI-entry lanes' difference families and the e2e tables as in R1, job by job (A5) |
| Ceiling | 60 core-hours (3 jobs x 300 min x 4 cores); expected 13.2. One correction pass in total for the ticket |

Results are recorded in `runs.tsv` and `lane-families.tsv` on `kit-176-launcher`.
