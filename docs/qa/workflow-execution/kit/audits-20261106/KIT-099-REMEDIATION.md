# KIT-099 remediation (2026-11-06)

Mode: `audit_and_repair`, on the user's "fix all recommendations" (2026-11-06), relayed by the coordinator. The audit is [KIT-099-AUDIT.md](KIT-099-AUDIT.md), unchanged. The only files edited are `tickets/kit/KIT-099.md` (revision 1.0 → 1.1) and this file.

## For the coordinator to sync (index and plan)

- **Title changed:**
  - **was:** "Transcriptomics at size: large journeys, partitioned defaults decided, and the four KIT-088 sizes fixed"
  - **now:** "Transcriptomics at size: the four KIT-088 sizes fixed, large journeys and the variant pairs recorded".
  - **Index title:** suggested "Transcriptomics at size: the four KIT-088 sizes, large journeys, the variant pairs recorded".
- **Description changed:** measurement only; the default change is KIT-100's.
- **Dependencies changed:** `["KIT-098"]` → `["KIT-098", "KIT-090"]` (F18). The validator now reports `frontmatter dependencies do not match index dependencies` until the index entry is synced. That is the only validation error; the OpenAI schema check is valid.
- **Plan rows:**
  - KIT-099's row: dependencies KIT-098 and KIT-090; runs "under its run gate (declaration, ledger row, user's go)".
  - A ledger row for KIT-099, when its declaration is written.
  - F14 and F15 are yours, as instructed.
- **KIT-100:** the new stub you create. KIT-099 hands it the pairs.
- **KIT-091:** the scope extension you are making (transcriptomics RUV grid: SVD and `cancor` loops). KIT-099 routes any transcriptomics RUV or grid ceiling to it.
- **KIT-098:** v1.8 must extend §8.2 and §8.7–8.9 to transcriptomics operation members. If it doesn't, KIT-099's run gate stops the work (F02).

## Findings

| ID | Severity | Repair | Where | State |
|---|---|---|---|---|
| F01 | blocking | A run gate: no member runs, the stepped trial included, without the full declaration, a ledger row and the user's go. The trial is declared with its points, run kind, ceilings, attempts, host record and the rule from trial to sizes. | KIT-099 "Run gate", "Declaration" 2 | closed in ticket; the ledger row waits on authoring (coordinator, plan) |
| F02 | blocking | The ticket binds protocol v1.8 (KIT-098) and lists the §8 rules it needs for transcriptomics operation members; if v1.8 lacks them, stop. | KIT-099 "Run gate" | routed to KIT-098 (v1.8); closed in KIT-099 |
| F03 | should-fix | Size rule: raw genes held at a declared realistic count (about 20k–60k), a declared low-count fraction, samples stepped; per size, raw and retained genes, samples, raw and retained cells, compared with the two-assay omics. | KIT-099 "Declaration" 3 | closed |
| F04 | should-fix | Terminal outcomes: a limit stop is a valid not-finishing result (step and peak reached recorded, never retried, no ratio, even though `kitIneligible` labels it `invalid`). Correction only for instrumentation. Sizes are never reduced. Rule for ultra. | "Decision rules and terminal outcomes"; "Declaration" 4 | closed |
| F05 | should-fix | Split: KIT-099 measures only. The default change moved to KIT-100. `R/utils_workflow_operations.R` removed from the file list. "Defaults change only on the user's decision" replaced by "no default changes (KIT-100's)". | Frontmatter, outcome, backward compatibility | closed; KIT-100 stub by the coordinator |
| F06 | should-fix | Pairs run at large size; `kitOperationDecision` recorded verbatim on §3 (≤ 1.02 kernel and whole-slice anonymous; ≤ 0.95 kernel, reported); KIT-080's rule shown with its real-data divergence; the logCPM branch explicit; one attempt, control first, cache caveat. Acceptance now records decisions and hands them to KIT-100. | "Decision rules"; "Declaration" 5.4; acceptance | closed |
| F07 | should-fix | Each filter pair reports the pre-pass phase separately. The pre-pass must be the TMM-effective library size (notebook :166-170; `filterByExpr.DGEList` uses `lib.size * norm.factors`). If KIT-095 declares raw library sizes, stop. | "Decision rules" | closed in KIT-099; the pre-pass definition stays with KIT-093/095 (the coordinator may carry it to their audits) |
| F08 | should-fix | The report question is now "which step and which call set peak and time per size", with the candidates listed. A transcriptomics RUV or grid ceiling goes to KIT-091. | Outcome "Report"; hand-offs | closed; KIT-091 scope by the coordinator |
| F09 | should-fix | Journey small test (completion and forced failure) and a screen operation stage (`081a` pattern) before any large member. | "Declaration" 5.1–5.2; simulation setup | closed |
| F10 | should-fix | The KIT-088 hand-off table is specified. Large here = KIT-088's large. No `kit_086` control (single-build curve). No missing values. Reuse only when build, protocol and corpus hash all match. | Outcome "Sizes"; hand-offs | closed |
| F11 | should-fix | Claim limits: synthetic, recorded, uncounted; no real-data, safe-limit or RUV-validity claim; TRX-007 not discharged; the KIT ledger governs, not the TRX allowance; stated in the outcome's headline. | Scientific contract "Claim limits" | closed in KIT-099; a one-line note in the TRX plan that its run allowance is superseded is KIT-040's (archiving) or the coordinator's |
| F12 | should-fix | One large in-process oracle per bound build (14,400 s, one attempt); pair populations from its stored states. | "Declaration" 5.3 | closed |
| F13 | minor | File list: the operations registry removed, the harness test added. The default-variant field is KIT-094's (binding field) and its change is KIT-100's. | Frontmatter | closed |
| F14 | minor | Index title: yours to sync (suggested above). | — | routed to the coordinator |
| F15 | minor | Plan dependency rows for KIT-088 and KIT-050: yours to sync. | — | routed to the coordinator |
| F16 | minor | Reuse named: `single_084`, `pair_084`/`stage_084` generalised (not copied), `kitOperationDecision`; the operation runner script and its hash recorded in the bindings. | "Evidence and handoff" | closed |
| F17 | minor | BLAS provider, threads and host state recorded per member. Machine time is declared in the declaration; the effort estimate (6 h) is unchanged. | "Declaration" 1, 2, 6 | closed |
| F18 | minor | KIT-090 added as a dependency, so the bound build carries the clean-room NIPALS. | Frontmatter | closed; index sync by the coordinator |

Counts: 18 findings. 13 closed in the ticket (F01, F03–F06, F08–F13, F16, F17). 3 closed in the ticket with the owning ticket's part routed (F02 → KIT-098 v1.8, F07 → KIT-093/095, F18 → index sync). 2 routed entirely (F14, F15 → coordinator). None is unclosable. F05's KIT-100 and F08's KIT-091 extension depend on the coordinator's edits.

Validation: `next-openai-ticket.sh KIT-099` returns OpenAI schema valid, status consistent. Its only error is the expected dependency mismatch with the index until it is synced.
