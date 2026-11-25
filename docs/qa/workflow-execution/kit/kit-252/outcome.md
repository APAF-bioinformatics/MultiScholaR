# KIT-252 outcome: a managed step that ends without its own error, said in plain words

Ticket KIT-252 r1.2 (plan audit AUDIT-KIT-252-plan-20261125 and its Re-check 1). Branch `kit-252` from
`v0.9-staging` 16112917. No R was run on the laptop; every run below is a GitHub-hosted runner (R 4.6.1, reference
BLAS). `data/benchmark/` was not opened.

| Commit | What |
|---|---|
| ce9d359d | tests only: `tests/testthat/test-kit252-managed-failure-text.R` and the five declared edits in four files |
| 0f34035a | the fix: `workflowManagedFailureText(result, minutes)` and the explicit log line (`R/utils_workflow_execution.R`); NEWS; ledger row `VAL-262` |
| (this commit) | records |

## What was run (observed)

| Run | Binds | Purpose, ceiling | Observed |
|---|---|---|---|
| 45593919125 (R1) | ce9d359d, head only | failing-first; 1 job, ceiling 45 job-minutes (the request file stated 75, a precaution before the files' times were known; the run stayed within the ticket's 45) | 40.9 job-minutes. `test-kit252-managed-failure-text.R` 12 passed, 60 failed: every failure in E1a, E1b, E2a, E2b, E3, E4, E5, E5b, E6, E7 and "E1a in a chain", each at the base value "managed action ended: <status>" (E7 also: no log line "managed action <operation> ended: memory_guard"); G1, G2, G3 passed. `test-prot-norm-core.R` 307 failed (the branch point's 306 plus the A111-03 block's `expect_identical`, now showing base "managed action ended: memory_guard" against WP-1). `test-omic-da-module.R` 54 (53 plus the KIT-112 W3 expectation, base "managed action ended: publication_pending"). `test-workflow-execution.R` 2 (KIT-115 at the old `:834`, `:840`). `test-workflow-omics-application-reads.R` 1 (KIT-082, the old `:567`, base "managed action ended: failed"). `r1-tests.tsv`, `r1-failures-kit252.txt` |
| 45600826132 (R2) | head 0f34035a against base ce9d359d | regression scope (13 files by `grep -rlE "managed action ended\|installWorkflowManagedStep\|worker_failed\|pollWorkflowAction" tests/testthat`, helpers excluded, recounted at 16112917) plus the new file; 3 shards per side, check; ceiling 600 job-minutes | 224.2 job-minutes. `fixed`: the new file (72 passed, 0 failed), `test-prot-norm-core.R` (307 to 306), `test-omic-da-module.R` (54 to 53), `test-workflow-execution.R` (2 to 0), `test-workflow-omics-application-reads.R` (1 to 0). `same`: the other nine. R CMD check: 0 errors, 3 warnings, 3 notes on both sides, findings identical (`r2-check-compare.txt`). `r2-compare.tsv`, `r2-compare-failures.txt`, `r2-tests.tsv`, `r2-tests-base.tsv` |

**The base counts at the branch point** were read from merge train 3's run 45573711299 (head 8bce32aa, whose `R/` and
`tests/` are byte-identical to 16112917): `test-prot-norm-core.R` 306 failures, `test-omic-da-module.R` 53; neither
the A111-03 block nor the KIT-112 W3 block among them (stop rule 4 did not apply). Comparing that run's
`failures.txt` with R1's and R2's head, test by test with temporary paths masked: R1 holds exactly the branch point's
failures plus the one declared expectation in each file; R2's head holds exactly the branch point's failures (306 and
53, the same messages). `test-workflow-kit-harness.R`'s one failure is the same on R2's two sides and at the branch
point.

## Acceptance criteria

| # | Criterion | Met | Evidence |
|---|---|---|---|
| A1 | Base values from KIT-233 L3w and the existing tests; no new run | met | KIT-233 L3w ("Error applying precursor rollup: managed action ended: worker_failed"); R1 shows the base value for every status in the new tests |
| A2 | Tests-only commit: E rows and the five declared edits fail at the base values by assertion, the two files that already fail judged per test (base count plus one); G1 to G3 pass; all pass on the fix, the two files at the base count with the same messages | met | R1 (45593919125) and R2 (45600826132) as above |
| A3 | R2 base against head: every file `same` except the new file and the declared edits | met | `r2-compare.tsv`: the five `fixed` files are exactly the new file and the four declared-edit files; nine `same`; check identical |
| A4 | NEWS, ledger, `outcome.md` with "For the change report" | met | `NEWS.md` (### Fixes, under the one unreleased heading), `docs/methods-decision-ledger.md` row `VAL-262` after VAL-229, this file, `wording.md` |

Stop rules: (2) no E or G row showed the head moved on a single action's failure (G2: `protein_replicate_filtered`
after E7's click; G3: the lipidomics head stays at `lipid_normalized`, `lipid_ruv_corrected` absent); (3) no G row
failed, no unpredicted R2 verdict; (4) did not apply.

Core-hours: R1 2.73 + R2 14.95 = **17.68** (expected 21.33, ceiling 86.0). No correction pass, no instrument rerun.

## What is inferred (not executed)

- The test harnesses mock the poll's receipt; the statuses' real causes (2.1 of the ticket) are by reading. Under
  systemd a worker stopped at `RuntimeMaxSec` is taken to report `Result` "timeout" (systemd documentation).
- The sentences are true for their cause by reading the code cited in `wording.md`; they are DRAFTS for the owner's
  wording pass (VAL-41, row QC-8), each listed there.

## Repair after the implementation audit

AUDIT-KIT-252-impl-20261125 (I01 to I09) was applied in one commit of texts and records, without a run of its own
(the merge train's run executes the new file and the four declared-edit files on the merged head): WP-5 "Nothing was
run." became "It did not run." (I02); WP-8's "the application's log file has the details" became "the application's
log has the details" (I03), in `R/utils_workflow_execution.R`, the new test file's constants and the declared edit in
`test-workflow-omics-application-reads.R`; NEWS, the ledger row and `wording.md` no longer say a user is told of low
memory on the computer (I01); `wording.md` lists every frame of the sentences and the owner's open points (I04 to I06,
I08). Convention divergence (I10): 25 added lines in `R/utils_workflow_execution.R` exceed R.md's 100 characters
(longest 126) and the test file uses `function(...)` lambdas, following the existing file's idiom (411 lines over 100
at base).

## For the change report

- **What was wrong, who was affected.** Every user of the installed application (every omic, user and test mode):
  when a step run on the managed worker ended without an error of its own (the application's memory limit for one
  step, a step that could not start, the per-step time limit of 30 minutes, a process killed or crashed, a result
  committed whose files were not written, an internal check), the notice showed an internal code, for example
  "Error applying precursor rollup: managed action ended: worker_failed" (KIT-233 run L3w), with no statement of what
  happened, whether anything was saved, or what to do. For `publication_pending` a generic sentence would have said
  nothing was saved when the result had been committed.
- **Legacy or new.** New to 0.9: `main` has neither `R/utils_workflow_execution.R` nor `R/utils_workflow_omics_actions.R`
  (`git show main:<path>`: not in `main`), so the managed worker and this notice do not exist in the legacy pipeline.
  The old literal was in the `v0.9.0` staging-build tag (12560b8a, `R/utils_workflow_execution.R:963`), which is not on
  `main` (AUDIT-KIT-252-impl I12).
- **Code.** Base 16112917: `R/utils_workflow_execution.R:960-966` (`paste("managed action ended:", result$status)`).
  Fix 0f34035a: `workflowManagedFailureText` at `R/utils_workflow_execution.R:824`; the log line
  "managed action <operation> ended: <status>" at `:1008`; the call at `:1011`; the action records its operation at
  `:945`.
- **What changed in results.** None: texts only. The prefix "Error applying <step>:", the notice's type and duration,
  the worker's own error when it wrote one (A111-03), the statuses, commits, journals and the `publication_pending`
  retry are unchanged (R2: every other test in the 14 files unchanged; G1 to G3).
- **Evidence.** R1 45593919125, R2 45600826132 (tables in this folder); branch-point counts from 45573711299.
  Limits: synthetic fixtures (the proteomics "small" recipe, the lipidomics declaration fixture), the poll mocked;
  no installed application or real worker was run.

## Not checked

- Which statuses a real user meets, and the notice on the installed application (no journey run; KIT-233 L3w is the
  only observation of the base value).
- The Windows plain launcher's statuses (KIT-237's runs may show them); systemd's `Result` at `RuntimeMaxSec`.
- The full suite (rides the merge train).
- The owner's final words (QC-8).
- Under systemd the application's own memory guards stop the whole application (AUDIT-KIT-252-impl I01); the
  low-memory sentence is mapped but not shown.

## Correction (coordinator, 2026-11-26)

A4's evidence cell says the ledger row `VAL-262` sits "after VAL-229". It sits after VAL-261: it
was numbered when merge train 7 landed (42846c0e, 2026-11-26).
