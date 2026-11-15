# KIT-124 audit of 2026-11-16: remediation record

**Authority for the repairs:** user, 2026-11-16, on applying the repairs and on `sva` as a comparator: "those decisions should follow best practices and best practices only". The audit register (`audit-20261116.md`) is unchanged. The audited revision 1.0 is kept as `KIT-124-rev1.0-audited.md` (sha256 `e3eb05d41b98c9fa…`). All repairs are edits to tickets; no package code, test or evidence file was changed.

| Finding | Repair, in KIT-124 revision 1.1 unless said otherwise | Closed by |
|---|---|---|
| F-01 | Fixed limits in place of "two standard errors": null rate at most 0.075, excess false discovery proportion at most 0.02, power loss at most 0.02, each judged by a one-sided 97.5% Monte Carlo bound. Repetitions set from a 50-run pilot to a stated precision, 200 to 1,000. Pass, fail and inconclusive defined. | "W3. Decision rule"; "Repetitions and seeds" |
| F-02 | The runs of 2026-11-16 are labelled exploratory in the ticket; confirming runs start at seed 100001. | "Exploratory evidence"; "Repetitions and seeds" |
| F-03 | Comparator arms added: the true factor in the model, `ruv::RUV4`, `ruv::RUVinv`, limma with surrogate variables, the known batch in the model; `svaseq` and the RUVg vignette design for counts. `sva` in a scratch library. | "Methods (arms)" 2, 14 to 18 |
| F-04 | Candidates that do not use the result: all features as controls, a supplied list (50 known nulls), technical replicates as the replicate structure; k fixed at 1 or 2 or estimated. The shared objective of the sweeps is named. Pooled-QC replicates are inventoried in W0 and enter S4's specification. | Arms 8 to 13; "The five implementations today"; W0 item 3 |
| F-05 | Primary cells added: the factor correlated with group at 0.3 and 0.6, three factors, heavy-tailed noise, a known unbalanced batch; drift for S4. | "Primary cells" |
| F-06 | An arm whose k may be 0; a no-unwanted-factor cell; a separate "no harm" rule with a fixed limit. | Arm 13; rule 4 |
| F-07 | The ticket ends at the user's gate. W4 and W5 are removed; implementation is authored after the gate from W0's consumer table. Estimate 40 hours. | "Bounded outcome"; "Gate" |
| F-08 | Measures added: null rate at 0.05 and 0.01, share of runs above 10% false discoveries, bias and error of the log2 differences, failure rate, each with its Monte Carlo standard error. | "Performance measures" |
| F-09 | The split is of biological units with technical replicates kept together, at least 3 per false group; measures and limit stated; the permutation fallback marked weaker. | "Null split"; rule 6 |
| F-10 | 16 primary cells per simulator (a centre, one setting moved at a time); 60 core-hours per simulator; the order in which arms are dropped. | "Primary cells"; "Compute" |
| F-11 | A consumer table is a W0 deliverable and the basis of the later implementation ticket. | W0 item 4 |
| F-12 | The first column is "Reference usage" with a source per cell and "example" where it is one; the references read through summaries and the four standards are re-read at source in W0. | "The five implementations today"; W0 item 1; "Standards" |
| F-13 | W0 inventories a known batch in the model, replicate and pooled-QC designs and the peptide method's use; a batch cell and a technical-replicate cell are primary; one end-to-end null split through the worker operations. | W0 item 3; "Primary cells"; "Null split" |
| F-14 | S4's model and comparators are fixed in a revision after W0, before S4 runs. | "Data-generating mechanisms" |
| F-15 | The series' five sections added to KIT-124 and KIT-123 (revision 1.2); the index entries' description, priority and estimate corrected; KIT-123 retitled. Inventory script: 0 findings. | `tickets/kit/KIT-123.md`, `tickets-index.json` |
| F-16 | A sidecar per result file: commit, library, versions, command, seeds. | "Repetitions and seeds" |
| F-17 | Retitled: "... up to the user's decision". | Title |
| F-18 | Closed in the audit. The checked fact is now in the table. | "The five implementations today" |

**Added beyond the findings, from the audit's stated limit (same author and auditor):** an independent read of the ticket and simulators before W3 runs ("Independent read").

**Not established by these repairs:** that any RUV design is calibrated or is best practice. The limits 0.075 and 0.02 and the choice of primary cells are stated judgements, fixed before any run; the four standards are cited from memory until W0 re-reads them.

**States:** implementation verified: no. Execution complete: no. Decision supported: no. Promotion authorised: no. Claim scope: the ticket's design is repaired against the recorded findings.
