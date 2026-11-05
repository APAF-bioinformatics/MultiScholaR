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

## Repairs (2026-11-06)

The user decided five points and said "fix all recommendations":
1. RUVSeq and EDASeq go into the kit library only (KIT-040).
2. The shared DA saving accepts an edgeR QL fit and records the real method (KIT-095).
3. The lipid modules are parameterised for a new omic, not copied (KIT-096, KIT-097).
4. RUVIII is a product-only export after DA, with `M` recorded (KIT-095, KIT-097).
5. KIT-092 is protocol v1.7 and also takes the KIT098-F02 exactness gaps; transcriptomics is v1.8 (KIT-098).

Each auditor repaired its own ticket, and the remediation record sits beside each audit. Every finding is closed in the ticket text or routed to a named owner; none is left open. The evidence for each closure arrives when its ticket runs.

**Specifier's share:**
- KIT-092 revision 1.1 (the exactness gaps, v1.7).
- KIT-091 revision 1.1 (the transcriptomics RUV grid at size).
- KIT-100, new: the default change split out of KIT-099.
- Plan 3.9: the stale rows, the TRX paragraph and the decision row.
- The index synced from each ticket's frontmatter.

**Consistency pass:**
- KIT-040 freezes the products KIT-093 names and writes `kit-040/product-map.tsv` mapping kernel names to notebook objects. It adds the derived filter values, the k cap, the annotation join and the QL table.
- KIT-095 drops the `skip` RUV mode, since the notebook has none and KIT-097 offers none.
- The class-declaration entry is KIT-093's.
- The Glimma outcome is recorded explicitly.
- KIT-088's stub waits for KIT-099.

**Validation:** all twelve touched tickets pass `next-openai-ticket.sh`, and `audit_inventory.py` finds no defects.

**Open with the user:** KIT-093 refuses designs with more than two groups (`trx_design_groups_unsupported`), because the notebook's `coef = 2` control ranking defines a two-group test only. Allowing more changes the science.
