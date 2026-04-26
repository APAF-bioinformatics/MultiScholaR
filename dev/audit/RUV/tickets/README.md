# RUV Fix-Plan Ticket Backlog

These tickets are the repo-local translation of `dev/audit/RUV/fix-plan.md` into actionable work items. There is no `goyoke` or `/ticket` scaffold in this repository today, so the backlog uses a structured Markdown format intended to be easy to import or restate into that apparatus.

## Authority and reconciliation

- Primary source of truth: `dev/audit/RUV/fix-plan.md`
- Verification sources: `dev/audit/RUV/math-audit.md`, `dev/audit/RUV/router-evaluation.md`, `dev/audit/RUV/staff-architect-handoff.md`
- Reconciled against live repo state on `2026-04-26`
- Important correction: the peptide `ctrl`-override claim in `staff-architect-review.md` is retracted and is not a valid implementation target

## Queue

| ID | Priority | Depends on | Summary |
| --- | --- | --- | --- |
| [RUV-001](./RUV-001-shared-k-selection-and-scoring-hardening.md) | P0 | none | Harden shared curve parsing, K selection, and composite scoring contracts |
| [RUV-002](./RUV-002-proteomics-automatic-optimizer-hardening.md) | P0 | RUV-001 | Bring proteomics automatic optimization onto the hardened contract and trace schema |
| [RUV-003](./RUV-003-metabolomics-automatic-search-rewrite.md) | P0 | RUV-001 | Replace metabolomics automatic one-shot behavior with real percentage search |
| [RUV-004](./RUV-004-lipidomics-automatic-search-rewrite.md) | P0 | RUV-001 | Replace lipidomics automatic one-shot behavior with real percentage search |
| [RUV-005](./RUV-005-peptide-optimizer-reconciliation.md) | P0 | RUV-001 | Guard peptide grouping-column handling and move automatic optimization onto the full objective |
| [RUV-006](./RUV-006-docs-ui-deprecations-and-release-notes.md) | P1 | RUV-001, RUV-002, RUV-003, RUV-004, RUV-005 | Align UIs, workbooks, roxygen, NEWS, and versioning with the corrected behavior |
| [RUV-007](./RUV-007-regression-matrix-and-release-gate.md) | P0 | RUV-001, RUV-002, RUV-003, RUV-004, RUV-005, RUV-006 | Update characterization tests, add regression coverage, and enforce the release gate |

## Global current-state notes

- `R/func_prot_norm_optimization_helpers.R` still uses argmax `findBestK()` logic and still allows vector `best_k` on ties.
- `R/func_metab_norm.R` and `R/func_lipid_norm_ruv_helpers.R` still create `percentage_range` in automatic mode and then ignore it, always returning `best_percentage = percentage_max`.
- `R/func_pept_s4_norm_methods.R` still runs peptide optimization through `ruvCancorFast()`, still drops `group_id` unconditionally in the ANOVA-prep path, and still mirrors the old separation/composite helper math.
- `weighted_difference` is still exposed as a normal metric in proteomics, metabolomics, and lipidomics UIs.
- Several tests currently encode the broken behavior and must be updated as part of the work, not treated as ground truth.

## Intended usage

Each ticket file is already shaped for direct execution:

- the live-repo defect is restated plainly
- the file scope is explicit
- acceptance criteria are concrete
- the affected tests are named up front
- non-goals are called out where the audit history was noisy

Each ticket's "Tests to update or add" section describes the minimum test changes needed to keep the test suite passing after that ticket's code changes. RUV-007 then performs the comprehensive test matrix audit, fills coverage gaps, and gates the release.

If these are fed into an external `/ticket` workflow, keep the dependencies intact. `RUV-001` is the shared foundation and should land before any omics-specific automatic-mode rewrite.
