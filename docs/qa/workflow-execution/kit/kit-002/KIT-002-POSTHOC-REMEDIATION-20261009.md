# KIT-002 post-hoc audit: remediation record

Register: [KIT-002-POSTHOC-AUDIT-20261009.md](KIT-002-POSTHOC-AUDIT-20261009.md) (not rewritten). Remediation commit `80af6cdb`.

| Finding | State | Closure evidence |
|---|---|---|
| PB-01 | closed | `require_complete_design_runs` documented in the roxygen block of `peptideIntensityFilteringHelper`; `man/peptideIntensityFilteringHelper.Rd` regenerated (only that Rd kept; roxygen also wanted to touch `ProteinQuantitativeData-class.Rd` and delete `runTestsContrasts.Rd`, pre-existing drift left untouched). `R CMD check` from a clean worktree of `80af6cdb` under `~/.cache`, proof library on `R_LIBS`, suggests not forced: 3 warnings, 3 notes ([r-cmd-check-remediated.log](r-cmd-check-remediated.log)), the CORE-018 set (`cli` `::` use, Makevars flags, GNU extensions). |
| PB-02 | closed | The remediation check was built and run under `~/.cache`; the README now says so. |
| PB-03 | closed | Check recipe recorded in `docs/qa/workflow-execution/kit/README.md` ("Checks that need the proof environment"). |
| PB-04 | closed | `checkpointDIAExecutionResult` resolves the declaration first, so an undeclared class is refused before any component write; test "an undeclared S4 result is refused with a typed class before any component is written" (declaration file now 50 assertions; checkpoint state file 238 pass). |
| PB-05 | informational | KIT-003 decides which callers may see the new classes when it adds them. |

Pre-existing drift noted, not touched: roxygen comments and committed Rd differ for `ProteinQuantitativeData-class.Rd` and `runTestsContrasts.Rd`; a future documentation pass should regenerate them deliberately.
