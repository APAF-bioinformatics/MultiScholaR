# KIT-017 pre-execution audit remediation (2026-10-15)

Audit: [KIT-017-PREEXEC-AUDIT-20261015.md](KIT-017-PREEXEC-AUDIT-20261015.md) (unchanged).

| ID | Repair | Closure evidence |
|---|---|---|
| PE-01 | Bindings regenerated at the closure tree (driver, `workflow_kit.R` and harness test hashes current); KIT-034's skipped regeneration recorded here. | `bindings.json` harness file hashes equal the working tree's; `validate kit-journey-017` passed. |
| PE-02 | Ticket revision 1.1 (control `kit_018`, candidate `kit_024`, rule `resident_end`, stage whole-slice reduction reported against the production-default baselines). | Ticket text; `bindings$tickets$017$rule = "resident_end"`; the user's launch decision. |
| PE-03 | The pair's lifecycle run on the screen corpus before the counted launch (`KIT_017_CORPUS=screen`); the harness pair branch takes the bound corpus key. | `test-workflow-kit-harness.R` 328 pass, 0 fail; screen stage result recorded in `outcome.md` once run. |
| PE-04 | Elapsed projection from the screen run recorded with the launch declaration. | `outcome.md`. |
| PE-05, PE-06, PE-10 | None (informational; recorded in the audit and the checklist). | — |
| PE-07 | Screen stage archived as `kit-journey-017-screen-lifecycle` before the counted launch. | Runtime directory listing in `outcome.md`. |
| PE-08 | Recorded in revision 1.1 `files_to_modify`. | Ticket frontmatter. |
| PE-09 | Full suite and `R CMD check` from a clean worktree at the binding commit, in the background. | `full-suite.txt`, `r-cmd-check.log` in this directory, compared with KIT-034's suite record and the KIT-002 check record. |
