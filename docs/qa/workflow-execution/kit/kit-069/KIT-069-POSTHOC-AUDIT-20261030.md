# KIT-069 post-hoc audit (2026-10-30, audit_only)

Scope: the completed KIT-069 (`tickets/kit/KIT-069.md` rev 1.0, completed at 00bc9550; code b1ff6015; install `kit_069` = `core/browser-smoke.4nx7p3`) against the ticket, [kit-064/AUDIT-ALIGNMENT-20261028.md](../kit-064/AUDIT-ALIGNMENT-20261028.md) and [REVIEW-20261028-stage3-recut.md](../REVIEW-20261028-stage3-recut.md). Audit only: no artefact of the ticket was changed; this register is additive. Ticket status is workflow metadata and counts as no evidence here.

## Mechanical and provenance (verified)

- `audit_inventory.py` on `tickets/kit` with `--ticket KIT-069`: 0 findings, exit 0 (schema, status, dependencies, graph).
- Bindings `builds.kit_069`: `source_commit` b1ff6015 = the install's `source-commit.txt`; library `core/browser-smoke.4nx7p3/library`; `lifecycle_omics.build` still `kit_030` (no member ran on kit_069; see A-01). The kit harness test on the regenerated bindings: 483 pass.
- `git show b1ff6015 --stat`: 16 R files, the test file, two evidence scripts. None of the 16 is named as a kernel in `kit-001/step-map-lipid.json` or `step-map-metab.json` (the kernels are the import readers, constructors, intensity filters, duplicate resolvers, the four normalisation kernels, correlation, RUV, DA and the report); the two S4 plotting-method files hold plot-list filters, the class KIT-064 item B1 treated as non-kernel for proteomics.
- R CMD check on the git-archive export of b1ff6015: 3 WARNINGs, 3 NOTEs = baseline.

## Implementation and tests (verified)

- BUG-10: `shiny::isolate(modules_initialized())` in both orchestrators. The real `workflow_data` is a `shiny::reactiveValues` (`mod_lipidomics.R:139`, `mod_metabolomics.R:139`), so the observer's remaining dependency is the `state_manager` field; `WorkflowState$getState()` reads no reactive. The `once = TRUE` trigger path keeps its own guard. Outcomes unchanged in every ordering.
- BUG-16: `seq_len` at both sites; identical to the lipid and proteomics copies.
- BUG-17: the floor is computed on `contrast_data$fdr_qvalue` before the `mutate`, the same column the old expression read; with any finite positive FDR the stand-in equals the old one (`is.finite` also drops `Inf`, which `> 0` kept, immaterial for a q-value); with none it is `.Machine$double.xmin`; NA FDR still maps to `1e-10`. The metab twin is the same six lines.
- BUG-23: 19 sites. `vapply` coerces an integer result into `numeric(1)` (checked: `vapply(1:2, function(i) 1L, numeric(1))` gives doubles) and `vapply(list(), nrow, integer(1))` is `integer(0)`, so every converted site returns what `sapply` simplified to on non-empty input and a typed empty vector on empty input.
- The test file reproduces every class on the 59f26056 worktree (7 failures, 3 errors) and passes on b1ff6015 (18); no existing test changed.

## Findings

| ID | Severity | Observed | Affected claim / impact | Smallest sufficient repair | Blocked state |
|---|---|---|---|---|---|
| A-01 | medium | The series plan's run ledger (`00-KIT-PLAN.md`, the "KIT-030, KIT-069, KIT-070, KIT-071" row) expects "the two KIT-068 lifecycle members on each build, recorded"; the ticket's Experiment budget names the install only; no member ran on kit_069 and the outcome's claim limits say so. | The ledger's expectation is unmet; the outcome claims nothing about the GUI journey on kit_069, so no recorded claim is false. KIT-070 runs its members on kit_070 whatever is decided. | User decision: run the two lifecycle members on kit_069 (about 4 minutes each) and file their stage summaries under kit-069/, or amend the ledger row to the tickets' wording. | decision (the ledger row), not implementation |
| A-02 | low | The module-CI comparison ran with `pkgload::load_all` of the two source trees against the kit_030 dependency library, not "on the installed build" as the ticket's budget words it. Same source commit as the install (b1ff6015); the install's own smoke member passed; KIT-030 used the same procedure. | The before/after exactness claim rests on the source trees; the installed library was not the test subject. | Record the reading in the outcome, or rerun the 12 files with `library(MultiScholaR)` from `4nx7p3/library` and file the counts. | none (claim scope) |
| A-03 | low | Installed over the kit_030 library (`U9nsh0`); the ticket says "over `kit_068`'s library", which does not exist (KIT-068 was tools-only and ran on kit_065). | None; kit_030 is the predecessor build. | Record in the outcome. | none |
| A-04 | low | `files_to_modify` names 4 R files and one test file; b1ff6015 touches 12 R files and a new test file. Item 4's text authorises the BUG-23 sweep over `R/mod_lipid*`, `R/mod_metab*`, `R/func_lipid*`, `R/func_metab*`, so the edits are in scope; the outcome records only the metab BUG-17 twin as a file-list divergence. | Traceability of the file list. | One sentence in the outcome naming the frontmatter list as the divergence source for the other files. | none |
| A-05 | low | Commit message b1ff6015 says "22 sapply sites"; the count is 19. Recorded in the outcome and in the closure commit's message. | None; the commit is immutable and pushed. | None. | none |
| A-06 | low (test design) | The BUG-10 test asserts total state reads = 2 (the tab renderer's one plus the observer's one) and sub-module initialisation = 1; the invalidation claim relies on the renderer's read being constant, stated in the test's comment (3 reads before, 2 after, measured). | The test is a proxy for "the observer ran once"; a renderer change would move the count without an observer defect. | Count the observer's reads alone (a state manager that records the calling frame, or a flag set inside the observer's branch). Not required for the ticket's acceptance. | none |
| A-07 | info | With no plot files the composite generator (all three omic copies) ends in its error handler (`split` of a zero-length grouping against `paste0(character(0), ")")`, which is `")"`). Recorded in the outcome; unreachable from the modules. | None. | None here; a shared-helper ticket if ever wanted. | none |
| A-08 | info | The 124-file lipid/metab sweep in one process showed 31 failures in 7 files that pass in their own processes (cross-file state leakage); the outcome records both. | The per-process result is the valid one; the single-process harness is weaker than the per-file runs the test rule names. | None. | none |

## States

- `implementation_verified`: yes (the four repairs match their findings; kernels untouched; tests reproduce and pass; check at baseline).
- `execution_complete`: yes for the ticket's budget (one install, the module-CI comparison); **no** for the series ledger's lifecycle members on kit_069 (A-01).
- `decision_supported`: the behaviour-preservation claim is supported by the module-CI content comparison (230 paths; every data product identical after masking timestamps, paths, ids and the provenance block) and the per-process suites; the all-zero-FDR product change is supported by the unit test only (no fixture reaches it), as the outcome states.
- `promotion_authorized`: not applicable (no counted member; nothing promoted).
- `claim_scope`: the module-CI lanes and the unit tests on the source trees; no GUI member on kit_069; the kernel `sapply` sites remain on the audit's list for the oracle tickets.
