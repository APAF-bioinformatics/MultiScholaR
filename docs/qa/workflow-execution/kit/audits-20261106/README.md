# Transcriptomics re-cut audit, 2026-11-06

Eight independent auditors, one per ticket (KIT-040, KIT-093 to 099), each using the `audit-evidence-contracts` skill in `audit_only` mode. They covered completeness, scientific rigour, implementation practice and codebase reference accuracy. Nothing was repaired; each report is in this folder.

| Ticket | Blocking | Should-fix | Minor | Verdict |
|---|---|---|---|---|
| KIT-040 corpus and notebook reference | 1 | 12 | 4 | not ready |
| KIT-093 class and kernels | 1 | 11 | 4 | not ready |
| KIT-094 per-binding step list | 0 | 8 | 3 | revise in planning |
| KIT-095 binding and declarations | 2 | 10 | 4 | not ready |
| KIT-096 GUI 1 | 1 | 7 | 7 | blocked on a decision |
| KIT-097 GUI 2 | 2 | 9 | 7 | not ready |
| KIT-098 kit route | 2 | 7 | 4 | not ready |
| KIT-099 at size | 2 | 10 | 6 | gated stub; not executable |
| **Total** | **11** | **74** | **39** | |

## Blocking findings, checked against the code by the specifier

| ID | Finding | Checked |
|---|---|---|
| KIT040-F01, KIT093-F01 | RUVSeq and EDASeq (and EDASeq's dependencies) are installed in no library, so the unchanged notebook can't run | `requireNamespace` FALSE in the default and kit libraries |
| KIT095-F01, KIT097-F01 | The DA products accept only limma: `newSavedWorkflowDaFit` needs an `MArrayLM`, and `model = "limma"` is hard-coded (`utils_workflow_omics_declarations.R:368,1090`). Reuse would save an edgeR fit as limma or abort | confirmed |
| KIT095-F02 | No ticket owns the transcriptomics renderers (CCA panels and workbook, PCA/RLE grid, volcano, bar plot, Glimma, DE tables, RUVIII grid) | from the ticket texts |
| KIT096-F01 | There is no shared design builder (lipid 6 files, metab 2, prot 7), so "reuse, not copy" and "other omics unchanged" conflict | confirmed |
| KIT097-F02 | RUVIII placement: the notebook runs it after DA for export; as a state it breaks the DA head rule, and as a product the record can't carry `M` | from the code and the notebook |
| KIT098-F01 | "Identical to the notebook reference" has no comparator, product-to-state map, or bound run settings | from the ticket text |
| KIT098-F02 | Exactness can pass unchecked: a step with no bound oracle is skipped (`workflow_kit.R:737`), and a non-DIA member with no oracles is "in memory", with exactness not judged (`:464`, `:619`) | confirmed |
| KIT099-F01 | The stepped trial isn't declared before it runs, and has no ledger row and no user go | from the ticket text |
| KIT099-F02 | Protocol v1.6 §8 covers lipid and metab only; nothing extends it to transcriptomics operation members | from the protocol |

## Errors in the specifier's plan that the audit found

- `ruv` is already in Imports (`DESCRIPTION:522`); only edgeR, RUVSeq and EDASeq are new.
- `filterByExpr` takes library size × TMM factors from the unfiltered counts, not plain library sizes (edgeR 4.10.1 source); the zero-variance filter and the second TMM pass were left out.
- The corpus's WT/KO groups become NA in the notebook, which hard-codes RPMI/Sera.
- "Prepared once, as KIT-085" points at NIPALS imputation, which doesn't apply to counts.
- Plan rows for KIT-088 and KIT-050 still show KIT-040, and the full-suite checkpoint line is stale.
