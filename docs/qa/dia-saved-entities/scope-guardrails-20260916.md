# Scope guardrail findings — 2026-09-16

Read-only review findings recorded before the user-authorized documentation repair.
Scope: AGENTS.md, the DIA execution plan/index and DIA-003. No runtime work or
acceptance promotion is authorized by this review.

| ID | Severity | Evidence and impact | Smallest repair / closure evidence |
| --- | --- | --- | --- |
| SG-01 | High | AGENTS.md protects scientific compatibility but lacks an executable scope/stop checklist. Local improvements can expand without a delivery endpoint. | Add bounded-work, necessity, verification and stop rules; inspect the resulting instructions. |
| SG-02 | High | DIA-003 revision 1.6 combines harness work with several approved repair amendments and a broad file list. This makes completed repairs appear continuously open. | Name the current comparison, its active files and exit condition; explicitly prohibit reopening unrelated repairs. |
| SG-03 | High | The plan has correction rules, but execution accumulated diagnostic modes and instrumentation repairs before useful comparisons. | Require small lifecycle preflights and fixed execution budgets; prohibit resetting budgets by relabeling work. Preserve all existing gates and failed evidence. |
| SG-04 | Low | Index plan_version is 1.1 while the plan is 1.2. | Synchronize the index with the new guardrail revision, without changing ticket statuses or dependencies. |

These are workflow-control findings, not scientific or performance findings.
The user authorized project instructions to prevent scope creep, speculative
optimization and overengineering, then explicitly requested continuation.
Original experiment protocols, numeric thresholds and historical results remain
unchanged. Full-series empirical acceptance is outside this documentation repair.
