# Ticket Index

Central index of all tickets across the repository. See `dev/ticket-schema.md` for the frontmatter contract.

## LIP — Lipidomics Scaffolding

Created: 2026-01-15. Cloned metabolomics infrastructure for lipidomics support.

| ID | Title | Status | Priority | Depends on |
|----|-------|--------|----------|------------|
| [LIP-001](TICKET-001_S4_Objects.md) | Lipidomics S4 objects and helpers | completed | P0 | — |
| [LIP-002](TICKET-002_Core_Modules.md) | Core data acquisition modules | completed | P0 | LIP-001 |
| [LIP-003](TICKET-003_QC_Norm.md) | Quality control and normalization modules | completed | P0 | LIP-001, LIP-002 |
| [LIP-004](TICKET-004_DE.md) | Differential expression | completed | P1 | LIP-003 |
| [LIP-005](TICKET-005_Orchestrator.md) | Orchestrator and summary | completed | P0 | LIP-001, LIP-002, LIP-003, LIP-004 |

## RUV — Automatic Parameter Selection Fix

Created: 2026-04-25. Derived from `dev/audit/RUV/fix-plan.md`. See `dev/audit/RUV/tickets/README.md` for series-specific notes.

| ID | Title | Status | Priority | Depends on |
|----|-------|--------|----------|------------|
| [RUV-001](audit/RUV/tickets/RUV-001-shared-k-selection-and-scoring-hardening.md) | Shared K-selection and scoring hardening | pending | P0 | — |
| [RUV-002](audit/RUV/tickets/RUV-002-proteomics-automatic-optimizer-hardening.md) | Proteomics automatic optimizer hardening | pending | P0 | RUV-001 |
| [RUV-003](audit/RUV/tickets/RUV-003-metabolomics-automatic-search-rewrite.md) | Metabolomics automatic search rewrite | pending | P0 | RUV-001 |
| [RUV-004](audit/RUV/tickets/RUV-004-lipidomics-automatic-search-rewrite.md) | Lipidomics automatic search rewrite | pending | P0 | RUV-001 |
| [RUV-005](audit/RUV/tickets/RUV-005-peptide-optimizer-reconciliation.md) | Peptide optimizer reconciliation | pending | P0 | RUV-001 |
| [RUV-006](audit/RUV/tickets/RUV-006-docs-ui-deprecations-and-release-notes.md) | Docs, UI deprecations, and release notes | pending | P1 | RUV-001–005 |
| [RUV-007](audit/RUV/tickets/RUV-007-regression-matrix-and-release-gate.md) | Regression matrix and release gate | pending | P0 | RUV-001–006 |
